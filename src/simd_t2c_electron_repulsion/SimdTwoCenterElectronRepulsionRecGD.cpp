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



#include "SimdTwoCenterElectronRepulsionRecGD.hpp"

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
compute_gd_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 2))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_gd_electron_repulsion: Basis functions must be of angular momenta four and two"));
    }

    if (harmonics.size() < 6)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_gd_electron_repulsion: Harmonics must reach angular momentum six"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_gd_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 8);

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
    // call, which fills the orders 0 to 6 of every row. The terms read the
    // orders 4 to 6 alone, and the orders below them are formed on the
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

            const auto ff_0 = fbase * bexp * bexp / (fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * bexp * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * bexp * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(5, i * nprim_b + j);

            const auto *bv_1 = boys.data(6, i * nprim_b + j);

            const auto *bv_2 = boys.data(7, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, pe_2, bv_0, bv_1, bv_2 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
                pe_2[k] += ff_2 * bv_2[k];
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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_10_7 = 10.0 / 7.0;
    const auto f_12_7 = 12.0 / 7.0;
    const auto f_17_14 = 17.0 / 14.0;
    const auto f_17_77 = 17.0 / 77.0;
    const auto f_18_7 = 18.0 / 7.0;
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_1_2 = 1.0 / 2.0;
    const auto f_1_21 = 1.0 / 21.0;
    const auto f_20_77 = 20.0 / 77.0;
    const auto f_2_1 = 2.0;
    const auto f_2_7 = 2.0 / 7.0;
    const auto f_3_1 = 3.0;
    const auto f_3_4 = 3.0 / 4.0;
    const auto f_3_7 = 3.0 / 7.0;
    const auto f_4_11 = 4.0 / 11.0;
    const auto f_4_21 = 4.0 / 21.0;
    const auto f_4_33 = 4.0 / 33.0;
    const auto f_4_7 = 4.0 / 7.0;
    const auto f_5_11 = 5.0 / 11.0;
    const auto f_8_77 = 8.0 / 77.0;
    const auto f_9_2 = 9.0 / 2.0;
    const auto fs_10_77_3 = std::sqrt(300.0 / 5929.0);
    const auto fs_1_11_14 = std::sqrt(14.0 / 121.0);
    const auto fs_1_11_5 = std::sqrt(5.0 / 121.0);
    const auto fs_1_11_6 = std::sqrt(6.0 / 121.0);
    const auto fs_1_11_7 = std::sqrt(7.0 / 121.0);
    const auto fs_1_14_30 = std::sqrt(15.0 / 98.0);
    const auto fs_1_21_15 = std::sqrt(5.0 / 147.0);
    const auto fs_1_21_30 = std::sqrt(10.0 / 147.0);
    const auto fs_1_21_35 = std::sqrt(5.0 / 63.0);
    const auto fs_1_22_110 = std::sqrt(5.0 / 22.0);
    const auto fs_1_22_2 = std::sqrt(1.0 / 242.0);
    const auto fs_1_2_6 = std::sqrt(3.0 / 2.0);
    const auto fs_1_33_105 = std::sqrt(35.0 / 363.0);
    const auto fs_1_33_15 = std::sqrt(5.0 / 363.0);
    const auto fs_1_33_165 = std::sqrt(5.0 / 33.0);
    const auto fs_1_33_210 = std::sqrt(70.0 / 363.0);
    const auto fs_1_33_35 = std::sqrt(35.0 / 1089.0);
    const auto fs_1_33_70 = std::sqrt(70.0 / 1089.0);
    const auto fs_1_42_10 = std::sqrt(5.0 / 882.0);
    const auto fs_1_42_70 = std::sqrt(5.0 / 126.0);
    const auto fs_1_66_10 = std::sqrt(5.0 / 2178.0);
    const auto fs_1_66_2 = std::sqrt(1.0 / 2178.0);
    const auto fs_1_66_330 = std::sqrt(5.0 / 66.0);
    const auto fs_1_66_70 = std::sqrt(35.0 / 2178.0);
    const auto fs_1_77_30 = std::sqrt(30.0 / 5929.0);
    const auto fs_1_7_21 = std::sqrt(3.0 / 7.0);
    const auto fs_2_11_3 = std::sqrt(12.0 / 121.0);
    const auto fs_2_21_5 = std::sqrt(20.0 / 441.0);
    const auto fs_2_33_30 = std::sqrt(40.0 / 363.0);
    const auto fs_2_33_42 = std::sqrt(56.0 / 363.0);
    const auto fs_2_77_21 = std::sqrt(12.0 / 847.0);
    const auto fs_3_14_10 = std::sqrt(45.0 / 98.0);
    const auto fs_3_14_21 = std::sqrt(27.0 / 28.0);
    const auto fs_3_14_70 = std::sqrt(45.0 / 14.0);
    const auto fs_3_2_5 = std::sqrt(45.0 / 4.0);
    const auto fs_3_4_15 = std::sqrt(135.0 / 16.0);
    const auto fs_3_4_30 = std::sqrt(135.0 / 8.0);
    const auto fs_3_4_35 = std::sqrt(315.0 / 16.0);
    const auto fs_3_77_21 = std::sqrt(27.0 / 847.0);
    const auto fs_3_7_15 = std::sqrt(135.0 / 49.0);
    const auto fs_3_7_30 = std::sqrt(270.0 / 49.0);
    const auto fs_3_7_35 = std::sqrt(45.0 / 7.0);
    const auto fs_3_8_10 = std::sqrt(45.0 / 32.0);
    const auto fs_3_8_70 = std::sqrt(315.0 / 32.0);
    const auto fs_4_33_7 = std::sqrt(112.0 / 1089.0);
    const auto fs_5_154_42 = std::sqrt(75.0 / 1694.0);
    const auto fs_5_28_42 = std::sqrt(75.0 / 56.0);
    const auto fs_5_33_7 = std::sqrt(175.0 / 1089.0);
    const auto fs_5_7_3 = std::sqrt(75.0 / 49.0);
    const auto fs_6_77_15 = std::sqrt(540.0 / 5929.0);
    const auto fs_6_7_5 = std::sqrt(180.0 / 49.0);
    const auto fs_9_154_6 = std::sqrt(243.0 / 11858.0);
    const auto fs_9_28_6 = std::sqrt(243.0 / 392.0);

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_p2, ph4_m4, ph4_m3, ph4_p2, ph4_p3, ph6_m5, ph6_m4, ph6_m3, ph6_p2, ph6_p3, ph6_p5, ph6_p6, ab_2, pc_0, pc_1, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];

        pc_0[k] = e_0 * fs_3_4_35 * h2_p2 + e_1 * fs_1_7_21 * h4_p2 - e_1 * fs_3_7_35 * r_2 * h2_p2 + e_2 * fs_1_66_2 * h6_p2 - e_2 * fs_1_22_110 * h6_p6 - e_2 * fs_2_77_21 * r_2 * h4_p2 + e_2 * fs_1_21_35 * r_4 * h2_p2;

        pc_1[k] = - e_1 * fs_1_2_6 * h4_p3 - e_2 * fs_1_22_2 * h6_p3 - e_2 * fs_1_66_330 * h6_p5 + e_2 * fs_1_11_6 * r_2 * h4_p3;

        pc_2[k] = e_1 * f_2_1 * h4_m4 + e_2 * fs_1_11_5 * h6_m4 - e_2 * f_4_11 * r_2 * h4_m4;

        pc_3[k] = - e_1 * fs_1_2_6 * h4_m3 + e_2 * fs_1_66_330 * h6_m5 - e_2 * fs_1_22_2 * h6_m3 + e_2 * fs_1_11_6 * r_2 * h4_m3;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_p1, ph2_p2, ph4_m2, ph4_p1, ph4_p2, ph4_p4, ph6_m6, ph6_m2, ph6_p1, ph6_p2, ph6_p4, ph6_p5, ab_2, pc_4, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];

        pc_4[k] = e_0 * fs_3_4_35 * h2_m2 + e_1 * fs_1_7_21 * h4_m2 - e_1 * fs_3_7_35 * r_2 * h2_m2 + e_2 * fs_1_22_110 * h6_m6 + e_2 * fs_1_66_2 * h6_m2 - e_2 * fs_2_77_21 * r_2 * h4_m2 + e_2 * fs_1_21_35 * r_4 * h2_m2;

        pc_5[k] = e_0 * fs_3_8_70 * h2_p1 + e_1 * fs_3_14_21 * h4_p1 - e_1 * fs_3_14_70 * r_2 * h2_p1 + e_2 * fs_1_66_10 * h6_p1 - e_2 * fs_1_33_165 * h6_p5 - e_2 * fs_3_77_21 * r_2 * h4_p1 + e_2 * fs_1_42_70 * r_4 * h2_p1;

        pc_6[k] = e_0 * fs_3_8_70 * h2_p2 - e_1 * fs_5_28_42 * h4_p2 + e_1 * fs_1_2_6 * h4_p4 - e_1 * fs_3_14_70 * r_2 * h2_p2 - e_2 * f_4_33 * h6_p2 - e_2 * fs_2_33_30 * h6_p4 + e_2 * fs_5_154_42 * r_2 * h4_p2 - e_2 * fs_1_11_6 * r_2 * h4_p4 + e_2 * fs_1_42_70 * r_4 * h2_p2;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ab_2, pc_7, pc_8, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];

        pc_7[k] = e_1 * f_1_2 * h4_m3 + e_2 * fs_2_11_3 * h6_m3 - e_2 * f_1_11 * r_2 * h4_m3;

        pc_8[k] = e_0 * fs_3_8_70 * h2_m2 - e_1 * fs_1_2_6 * h4_m4 - e_1 * fs_5_28_42 * h4_m2 - e_1 * fs_3_14_70 * r_2 * h2_m2 + e_2 * fs_2_33_30 * h6_m4 - e_2 * f_4_33 * h6_m2 + e_2 * fs_1_11_6 * r_2 * h4_m4 + e_2 * fs_5_154_42 * r_2 * h4_m2 + e_2 * fs_1_42_70 * r_4 * h2_m2;

        pc_9[k] = e_0 * fs_3_8_70 * h2_m1 + e_1 * fs_3_14_21 * h4_m1 - e_1 * fs_3_14_70 * r_2 * h2_m1 + e_2 * fs_1_33_165 * h6_m5 + e_2 * fs_1_66_10 * h6_m1 - e_2 * fs_3_77_21 * r_2 * h4_m1 + e_2 * fs_1_42_70 * r_4 * h2_m1;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_0, ph2_p1, ph4_m2, ph4_0, ph4_p1, ph4_p3, ph4_p4, ph6_m2, ph6_0, ph6_p1, ph6_p3, ph6_p4, ab_2, pc_10, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];

        pc_10[k] = e_0 * fs_3_4_15 * h2_0 + e_1 * fs_3_7_15 * h4_0 - e_1 * fs_1_7_21 * h4_p4 - e_1 * fs_3_7_15 * r_2 * h2_0 + e_2 * fs_1_33_15 * h6_0 - e_2 * fs_1_33_105 * h6_p4 - e_2 * fs_6_77_15 * r_2 * h4_0 + e_2 * fs_2_77_21 * r_2 * h4_p4 + e_2 * fs_1_21_15 * r_4 * h2_0;

        pc_11[k] = e_0 * fs_3_2_5 * h2_p1 - e_1 * fs_9_28_6 * h4_p1 + e_1 * fs_5_28_42 * h4_p3 - e_1 * fs_6_7_5 * r_2 * h2_p1 - e_2 * fs_1_33_35 * h6_p1 - e_2 * fs_1_11_14 * h6_p3 + e_2 * fs_9_154_6 * r_2 * h4_p1 - e_2 * fs_5_154_42 * r_2 * h4_p3 + e_2 * fs_2_21_5 * r_4 * h2_p1;

        pc_12[k] = e_0 * fs_3_4_15 * h2_m2 - e_1 * f_4_7 * h4_m2 - e_1 * fs_3_7_15 * r_2 * h2_m2 + e_2 * fs_2_33_42 * h6_m2 + e_2 * f_8_77 * r_2 * h4_m2 + e_2 * fs_1_21_15 * r_4 * h2_m2;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m1, ph2_p1, ph4_m4, ph4_m3, ph4_m1, ph4_p1, ph4_p3, ph6_m4, ph6_m3, ph6_m1, ph6_p1, ph6_p3, ab_2, pc_13, pc_14, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];

        pc_13[k] = e_0 * fs_3_2_5 * h2_m1 - e_1 * fs_5_28_42 * h4_m3 - e_1 * fs_9_28_6 * h4_m1 - e_1 * fs_6_7_5 * r_2 * h2_m1 + e_2 * fs_1_11_14 * h6_m3 - e_2 * fs_1_33_35 * h6_m1 + e_2 * fs_5_154_42 * r_2 * h4_m3 + e_2 * fs_9_154_6 * r_2 * h4_m1 + e_2 * fs_2_21_5 * r_4 * h2_m1;

        pc_14[k] = e_1 * fs_1_7_21 * h4_m4 + e_2 * fs_1_33_105 * h6_m4 - e_2 * fs_2_77_21 * r_2 * h4_m4;

        pc_15[k] = - e_0 * fs_3_8_10 * h2_p1 - e_1 * fs_5_7_3 * h4_p1 - e_1 * fs_3_14_21 * h4_p3 + e_1 * fs_3_14_10 * r_2 * h2_p1 - e_2 * fs_1_66_70 * h6_p1 - e_2 * fs_1_11_7 * h6_p3 + e_2 * fs_10_77_3 * r_2 * h4_p1 + e_2 * fs_3_77_21 * r_2 * h4_p3 - e_2 * fs_1_42_10 * r_4 * h2_p1;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph2_0, ph2_p2, ph4_m2, ph4_m1, ph4_0, ph4_p2, ph6_m2, ph6_m1, ph6_0, ph6_p2, ab_2, pc_16, pc_17, pc_18 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];

        pc_16[k] = e_0 * fs_3_4_30 * h2_0 + e_0 * fs_3_8_10 * h2_p2 - e_1 * fs_1_14_30 * h4_0 + e_1 * fs_9_28_6 * h4_p2 - e_1 * fs_3_7_30 * r_2 * h2_0 - e_1 * fs_3_14_10 * r_2 * h2_p2 - e_2 * fs_2_33_30 * h6_0 - e_2 * fs_4_33_7 * h6_p2 + e_2 * fs_1_77_30 * r_2 * h4_0 - e_2 * fs_9_154_6 * r_2 * h4_p2 + e_2 * fs_1_21_30 * r_4 * h2_0 + e_2 * fs_1_42_10 * r_4 * h2_p2;

        pc_17[k] = e_0 * fs_3_4_30 * h2_m1 - e_1 * f_17_14 * h4_m1 - e_1 * fs_3_7_30 * r_2 * h2_m1 + e_2 * fs_1_33_210 * h6_m1 + e_2 * f_17_77 * r_2 * h4_m1 + e_2 * fs_1_21_30 * r_4 * h2_m1;

        pc_18[k] = - e_0 * fs_3_8_10 * h2_m2 - e_1 * fs_9_28_6 * h4_m2 + e_1 * fs_3_14_10 * r_2 * h2_m2 + e_2 * fs_4_33_7 * h6_m2 + e_2 * fs_9_154_6 * r_2 * h4_m2 - e_2 * fs_1_42_10 * r_4 * h2_m2;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph2_0, ph4_m3, ph4_m2, ph4_m1, ph4_0, ph6_m3, ph6_m2, ph6_m1, ph6_0, ab_2, pc_19, pc_20, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];

        pc_19[k] = e_0 * fs_3_8_10 * h2_m1 + e_1 * fs_3_14_21 * h4_m3 + e_1 * fs_5_7_3 * h4_m1 - e_1 * fs_3_14_10 * r_2 * h2_m1 + e_2 * fs_1_11_7 * h6_m3 + e_2 * fs_1_66_70 * h6_m1 - e_2 * fs_3_77_21 * r_2 * h4_m3 - e_2 * fs_10_77_3 * r_2 * h4_m1 + e_2 * fs_1_42_10 * r_4 * h2_m1;

        pc_20[k] = e_0 * f_3_4 * h2_m2 + e_1 * fs_3_7_15 * h4_m2 - e_1 * f_3_7 * r_2 * h2_m2 + e_2 * fs_1_33_70 * h6_m2 - e_2 * fs_6_77_15 * r_2 * h4_m2 + e_2 * f_1_21 * r_4 * h2_m2;

        pc_21[k] = - e_0 * f_3_1 * h2_m1 - e_1 * fs_1_14_30 * h4_m1 + e_1 * f_12_7 * r_2 * h2_m1 + e_2 * fs_5_33_7 * h6_m1 + e_2 * fs_1_77_30 * r_2 * h4_m1 - e_2 * f_4_21 * r_4 * h2_m1;

        pc_22[k] = e_0 * f_9_2 * h2_0 - e_1 * f_10_7 * h4_0 - e_1 * f_18_7 * r_2 * h2_0 + e_2 * f_5_11 * h6_0 + e_2 * f_20_77 * r_2 * h4_0 + e_2 * f_2_7 * r_4 * h2_0;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m1, ph2_p1, ph2_p2, ph4_m3, ph4_m1, ph4_p1, ph4_p2, ph6_m3, ph6_m1, ph6_p1, ph6_p2, ab_2, pc_23, pc_24, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];

        pc_23[k] = - e_0 * f_3_1 * h2_p1 - e_1 * fs_1_14_30 * h4_p1 + e_1 * f_12_7 * r_2 * h2_p1 + e_2 * fs_5_33_7 * h6_p1 + e_2 * fs_1_77_30 * r_2 * h4_p1 - e_2 * f_4_21 * r_4 * h2_p1;

        pc_24[k] = e_0 * f_3_4 * h2_p2 + e_1 * fs_3_7_15 * h4_p2 - e_1 * f_3_7 * r_2 * h2_p2 + e_2 * fs_1_33_70 * h6_p2 - e_2 * fs_6_77_15 * r_2 * h4_p2 + e_2 * f_1_21 * r_4 * h2_p2;

        pc_25[k] = - e_0 * fs_3_8_10 * h2_m1 + e_1 * fs_3_14_21 * h4_m3 - e_1 * fs_5_7_3 * h4_m1 + e_1 * fs_3_14_10 * r_2 * h2_m1 + e_2 * fs_1_11_7 * h6_m3 - e_2 * fs_1_66_70 * h6_m1 - e_2 * fs_3_77_21 * r_2 * h4_m3 + e_2 * fs_10_77_3 * r_2 * h4_m1 - e_2 * fs_1_42_10 * r_4 * h2_m1;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_0, ph2_p1, ph2_p2, ph4_m2, ph4_0, ph4_p1, ph4_p2, ph6_m2, ph6_0, ph6_p1, ph6_p2, ab_2, pc_26, pc_27, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];

        pc_26[k] = - e_0 * fs_3_8_10 * h2_m2 - e_1 * fs_9_28_6 * h4_m2 + e_1 * fs_3_14_10 * r_2 * h2_m2 + e_2 * fs_4_33_7 * h6_m2 + e_2 * fs_9_154_6 * r_2 * h4_m2 - e_2 * fs_1_42_10 * r_4 * h2_m2;

        pc_27[k] = e_0 * fs_3_4_30 * h2_p1 - e_1 * f_17_14 * h4_p1 - e_1 * fs_3_7_30 * r_2 * h2_p1 + e_2 * fs_1_33_210 * h6_p1 + e_2 * f_17_77 * r_2 * h4_p1 + e_2 * fs_1_21_30 * r_4 * h2_p1;

        pc_28[k] = e_0 * fs_3_4_30 * h2_0 - e_0 * fs_3_8_10 * h2_p2 - e_1 * fs_1_14_30 * h4_0 - e_1 * fs_9_28_6 * h4_p2 - e_1 * fs_3_7_30 * r_2 * h2_0 + e_1 * fs_3_14_10 * r_2 * h2_p2 - e_2 * fs_2_33_30 * h6_0 + e_2 * fs_4_33_7 * h6_p2 + e_2 * fs_1_77_30 * r_2 * h4_0 + e_2 * fs_9_154_6 * r_2 * h4_p2 + e_2 * fs_1_21_30 * r_4 * h2_0 - e_2 * fs_1_42_10 * r_4 * h2_p2;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m1, ph2_p1, ph4_m4, ph4_m3, ph4_m1, ph4_p1, ph4_p3, ph6_m4, ph6_m3, ph6_m1, ph6_p1, ph6_p3, ab_2, pc_29, pc_30, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];

        pc_29[k] = - e_0 * fs_3_8_10 * h2_p1 - e_1 * fs_5_7_3 * h4_p1 + e_1 * fs_3_14_21 * h4_p3 + e_1 * fs_3_14_10 * r_2 * h2_p1 - e_2 * fs_1_66_70 * h6_p1 + e_2 * fs_1_11_7 * h6_p3 + e_2 * fs_10_77_3 * r_2 * h4_p1 - e_2 * fs_3_77_21 * r_2 * h4_p3 - e_2 * fs_1_42_10 * r_4 * h2_p1;

        pc_30[k] = e_1 * fs_1_7_21 * h4_m4 + e_2 * fs_1_33_105 * h6_m4 - e_2 * fs_2_77_21 * r_2 * h4_m4;

        pc_31[k] = - e_0 * fs_3_2_5 * h2_m1 - e_1 * fs_5_28_42 * h4_m3 + e_1 * fs_9_28_6 * h4_m1 + e_1 * fs_6_7_5 * r_2 * h2_m1 + e_2 * fs_1_11_14 * h6_m3 + e_2 * fs_1_33_35 * h6_m1 + e_2 * fs_5_154_42 * r_2 * h4_m3 - e_2 * fs_9_154_6 * r_2 * h4_m1 - e_2 * fs_2_21_5 * r_4 * h2_m1;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph6_p4, ab_2, pc_32, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];

        pc_32[k] = e_0 * fs_3_4_15 * h2_p2 - e_1 * f_4_7 * h4_p2 - e_1 * fs_3_7_15 * r_2 * h2_p2 + e_2 * fs_2_33_42 * h6_p2 + e_2 * f_8_77 * r_2 * h4_p2 + e_2 * fs_1_21_15 * r_4 * h2_p2;

        pc_33[k] = e_0 * fs_3_2_5 * h2_p1 - e_1 * fs_9_28_6 * h4_p1 - e_1 * fs_5_28_42 * h4_p3 - e_1 * fs_6_7_5 * r_2 * h2_p1 - e_2 * fs_1_33_35 * h6_p1 + e_2 * fs_1_11_14 * h6_p3 + e_2 * fs_9_154_6 * r_2 * h4_p1 + e_2 * fs_5_154_42 * r_2 * h4_p3 + e_2 * fs_2_21_5 * r_4 * h2_p1;

        pc_34[k] = e_0 * fs_3_4_15 * h2_0 + e_1 * fs_3_7_15 * h4_0 + e_1 * fs_1_7_21 * h4_p4 - e_1 * fs_3_7_15 * r_2 * h2_0 + e_2 * fs_1_33_15 * h6_0 + e_2 * fs_1_33_105 * h6_p4 - e_2 * fs_6_77_15 * r_2 * h4_0 - e_2 * fs_2_77_21 * r_2 * h4_p4 + e_2 * fs_1_21_15 * r_4 * h2_0;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph4_m4, ph4_m2, ph4_m1, ph4_p3, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph6_p3, ab_2, pc_35, pc_36, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p3 = ph6_p3[k];

        pc_35[k] = - e_0 * fs_3_8_70 * h2_m1 - e_1 * fs_3_14_21 * h4_m1 + e_1 * fs_3_14_70 * r_2 * h2_m1 + e_2 * fs_1_33_165 * h6_m5 - e_2 * fs_1_66_10 * h6_m1 + e_2 * fs_3_77_21 * r_2 * h4_m1 - e_2 * fs_1_42_70 * r_4 * h2_m1;

        pc_36[k] = - e_0 * fs_3_8_70 * h2_m2 - e_1 * fs_1_2_6 * h4_m4 + e_1 * fs_5_28_42 * h4_m2 + e_1 * fs_3_14_70 * r_2 * h2_m2 + e_2 * fs_2_33_30 * h6_m4 + e_2 * f_4_33 * h6_m2 + e_2 * fs_1_11_6 * r_2 * h4_m4 - e_2 * fs_5_154_42 * r_2 * h4_m2 - e_2 * fs_1_42_70 * r_4 * h2_m2;

        pc_37[k] = e_1 * f_1_2 * h4_p3 + e_2 * fs_2_11_3 * h6_p3 - e_2 * f_1_11 * r_2 * h4_p3;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_p1, ph2_p2, ph4_m2, ph4_p1, ph4_p2, ph4_p4, ph6_m6, ph6_m2, ph6_p1, ph6_p2, ph6_p4, ph6_p5, ab_2, pc_38, pc_39, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];

        pc_38[k] = e_0 * fs_3_8_70 * h2_p2 - e_1 * fs_5_28_42 * h4_p2 - e_1 * fs_1_2_6 * h4_p4 - e_1 * fs_3_14_70 * r_2 * h2_p2 - e_2 * f_4_33 * h6_p2 + e_2 * fs_2_33_30 * h6_p4 + e_2 * fs_5_154_42 * r_2 * h4_p2 + e_2 * fs_1_11_6 * r_2 * h4_p4 + e_2 * fs_1_42_70 * r_4 * h2_p2;

        pc_39[k] = e_0 * fs_3_8_70 * h2_p1 + e_1 * fs_3_14_21 * h4_p1 - e_1 * fs_3_14_70 * r_2 * h2_p1 + e_2 * fs_1_66_10 * h6_p1 + e_2 * fs_1_33_165 * h6_p5 - e_2 * fs_3_77_21 * r_2 * h4_p1 + e_2 * fs_1_42_70 * r_4 * h2_p1;

        pc_40[k] = - e_0 * fs_3_4_35 * h2_m2 - e_1 * fs_1_7_21 * h4_m2 + e_1 * fs_3_7_35 * r_2 * h2_m2 + e_2 * fs_1_22_110 * h6_m6 - e_2 * fs_1_66_2 * h6_m2 + e_2 * fs_2_77_21 * r_2 * h4_m2 - e_2 * fs_1_21_35 * r_4 * h2_m2;
    }

    // NOTE: the angular components are formed in 14 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_p2, ph4_m3, ph4_p2, ph4_p3, ph4_p4, ph6_m5, ph6_m3, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ph6_p6, ab_2, pc_41, pc_42, pc_43, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];

        pc_41[k] = e_1 * fs_1_2_6 * h4_m3 + e_2 * fs_1_66_330 * h6_m5 + e_2 * fs_1_22_2 * h6_m3 - e_2 * fs_1_11_6 * r_2 * h4_m3;

        pc_42[k] = e_1 * f_2_1 * h4_p4 + e_2 * fs_1_11_5 * h6_p4 - e_2 * f_4_11 * r_2 * h4_p4;

        pc_43[k] = - e_1 * fs_1_2_6 * h4_p3 - e_2 * fs_1_22_2 * h6_p3 + e_2 * fs_1_66_330 * h6_p5 + e_2 * fs_1_11_6 * r_2 * h4_p3;

        pc_44[k] = e_0 * fs_3_4_35 * h2_p2 + e_1 * fs_1_7_21 * h4_p2 - e_1 * fs_3_7_35 * r_2 * h2_p2 + e_2 * fs_1_66_2 * h6_p2 + e_2 * fs_1_22_110 * h6_p6 - e_2 * fs_2_77_21 * r_2 * h4_p2 + e_2 * fs_1_21_35 * r_4 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[45] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44};

    for (size_t n = 0; n < 45; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
