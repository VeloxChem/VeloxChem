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



#include "SimdTwoCenterElectronRepulsionRecPI.hpp"

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
compute_pi_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 1) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_pi_electron_repulsion: Basis functions must be of angular momenta one and six"));
    }

    if (harmonics.size() < 7)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_pi_electron_repulsion: Harmonics must reach angular momentum 7"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_pi_electron_repulsion: Number of values exceeds number of atom pairs"));
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
    // orders 6 to 7 alone, and the orders below them are formed on the
    // way to them by the recursion the Boys function is evaluated with.

    simdfunc::compute_boys_function(boys);

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time. Every exponent factor is
    // used with one order of the Boys function alone, so the order follows from
    // the factor and one accumulator per factor suffices.

    auto buffer = CSimdMatrix(2, nvalues);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);

    std::fill(pe_0, pe_0 + nvalues, 0.0);
    std::fill(pe_1, pe_1 + nvalues, 0.0);

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

            const auto ff_0 = fbase * aexp * aexp * aexp * aexp * aexp / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * aexp * aexp * aexp * aexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(7, i * nprim_b + j);

            const auto *bv_1 = boys.data(8, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, bv_0, bv_1 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
            }
        }
    }

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_3_1 = 3.0;
    const auto f_6_13 = 6.0 / 13.0;
    const auto f_7_13 = 7.0 / 13.0;
    const auto fs_1_13_11 = std::sqrt(11.0 / 169.0);
    const auto fs_1_13_13 = std::sqrt(1.0 / 13.0);
    const auto fs_1_13_14 = std::sqrt(14.0 / 169.0);
    const auto fs_1_13_15 = std::sqrt(15.0 / 169.0);
    const auto fs_1_13_21 = std::sqrt(21.0 / 169.0);
    const auto fs_1_13_3 = std::sqrt(3.0 / 169.0);
    const auto fs_1_13_33 = std::sqrt(33.0 / 169.0);
    const auto fs_1_13_35 = std::sqrt(35.0 / 169.0);
    const auto fs_1_13_39 = std::sqrt(3.0 / 13.0);
    const auto fs_1_13_5 = std::sqrt(5.0 / 169.0);
    const auto fs_1_1_5 = std::sqrt(5.0);
    const auto fs_1_26_110 = std::sqrt(55.0 / 338.0);
    const auto fs_1_26_182 = std::sqrt(7.0 / 26.0);
    const auto fs_1_26_2 = std::sqrt(1.0 / 338.0);
    const auto fs_1_26_30 = std::sqrt(15.0 / 338.0);
    const auto fs_1_26_6 = std::sqrt(3.0 / 338.0);
    const auto fs_1_2_11 = std::sqrt(11.0 / 4.0);
    const auto fs_1_2_14 = std::sqrt(7.0 / 2.0);
    const auto fs_1_2_15 = std::sqrt(15.0 / 4.0);
    const auto fs_1_2_21 = std::sqrt(21.0 / 4.0);
    const auto fs_1_2_3 = std::sqrt(3.0 / 4.0);
    const auto fs_1_2_33 = std::sqrt(33.0 / 4.0);
    const auto fs_1_2_35 = std::sqrt(35.0 / 4.0);
    const auto fs_1_2_5 = std::sqrt(5.0 / 4.0);
    const auto fs_1_4_110 = std::sqrt(55.0 / 8.0);
    const auto fs_1_4_2 = std::sqrt(1.0 / 8.0);
    const auto fs_1_4_6 = std::sqrt(3.0 / 8.0);
    const auto fs_2_13_10 = std::sqrt(40.0 / 169.0);
    const auto fs_2_13_5 = std::sqrt(20.0 / 169.0);
    const auto fs_2_13_6 = std::sqrt(24.0 / 169.0);
    const auto fs_2_13_7 = std::sqrt(28.0 / 169.0);
    const auto fs_2_1_2 = std::sqrt(8.0);
    const auto fs_3_13_2 = std::sqrt(18.0 / 169.0);
    const auto fs_3_13_3 = std::sqrt(27.0 / 169.0);
    const auto fs_3_13_5 = std::sqrt(45.0 / 169.0);
    const auto fs_3_26_10 = std::sqrt(45.0 / 338.0);
    const auto fs_3_2_2 = std::sqrt(9.0 / 2.0);
    const auto fs_3_2_3 = std::sqrt(27.0 / 4.0);
    const auto fs_3_4_10 = std::sqrt(45.0 / 8.0);
    const auto fs_4_13_2 = std::sqrt(32.0 / 169.0);
    const auto fs_4_13_3 = std::sqrt(48.0 / 169.0);

    // NOTE: the angular components are formed in 8 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ph7_p1, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ph7_p7, ab_2, pc_0, pc_1, pc_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_0[k] = e_0 * fs_1_2_33 * h5_p5 + e_1 * fs_1_26_2 * h7_p5 + e_1 * fs_1_26_182 * h7_p7 - e_1 * fs_1_13_33 * r_2 * h5_p5;

        pc_1[k] = e_0 * fs_1_4_110 * h5_p4 + e_1 * fs_1_26_6 * h7_p4 + e_1 * fs_1_13_39 * h7_p6 - e_1 * fs_1_26_110 * r_2 * h5_p4;

        pc_2[k] = e_0 * fs_3_4_10 * h5_p3 + e_0 * fs_1_4_2 * h5_p5 + e_1 * fs_1_13_3 * h7_p3 + e_1 * fs_1_13_33 * h7_p5 - e_1 * fs_3_26_10 * r_2 * h5_p3 - e_1 * fs_1_26_2 * r_2 * h5_p5;

        pc_3[k] = e_0 * fs_3_2_2 * h5_p2 + e_0 * fs_1_4_6 * h5_p4 + e_1 * fs_1_13_5 * h7_p2 + e_1 * fs_1_26_110 * h7_p4 - e_1 * fs_3_13_2 * r_2 * h5_p2 - e_1 * fs_1_26_6 * r_2 * h5_p4;

        pc_4[k] = e_0 * fs_1_2_14 * h5_p1 + e_0 * fs_1_2_3 * h5_p3 + e_1 * fs_1_26_30 * h7_p1 + e_1 * fs_3_26_10 * h7_p3 - e_1 * fs_1_13_14 * r_2 * h5_p1 - e_1 * fs_1_13_3 * r_2 * h5_p3;
    }

    // NOTE: the angular components are formed in 8 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph5_0, ph5_p2, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ph7_0, ph7_p2, ab_2, pc_5, pc_6, pc_7, pc_8, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];

        pc_5[k] = e_0 * fs_1_2_21 * h5_0 + e_0 * fs_1_2_5 * h5_p2 + e_1 * fs_1_13_21 * h7_0 + e_1 * fs_3_13_2 * h7_p2 - e_1 * fs_1_13_21 * r_2 * h5_0 - e_1 * fs_1_13_5 * r_2 * h5_p2;

        pc_6[k] = - e_0 * fs_1_2_15 * h5_m1 - e_1 * fs_2_13_7 * h7_m1 + e_1 * fs_1_13_15 * r_2 * h5_m1;

        pc_7[k] = - e_0 * fs_1_2_5 * h5_m2 - e_1 * fs_3_13_2 * h7_m2 + e_1 * fs_1_13_5 * r_2 * h5_m2;

        pc_8[k] = - e_0 * fs_1_2_3 * h5_m3 - e_0 * fs_1_2_14 * h5_m1 - e_1 * fs_3_26_10 * h7_m3 - e_1 * fs_1_26_30 * h7_m1 + e_1 * fs_1_13_3 * r_2 * h5_m3 + e_1 * fs_1_13_14 * r_2 * h5_m1;

        pc_9[k] = - e_0 * fs_1_4_6 * h5_m4 - e_0 * fs_3_2_2 * h5_m2 - e_1 * fs_1_26_110 * h7_m4 - e_1 * fs_1_13_5 * h7_m2 + e_1 * fs_1_26_6 * r_2 * h5_m4 + e_1 * fs_3_13_2 * r_2 * h5_m2;
    }

    // NOTE: the angular components are formed in 8 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph5_m5, ph5_m4, ph5_m3, ph7_m7, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ab_2, pc_10, pc_11, pc_12, pc_13, pc_14, pc_15, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];

        pc_10[k] = - e_0 * fs_1_4_2 * h5_m5 - e_0 * fs_3_4_10 * h5_m3 - e_1 * fs_1_13_33 * h7_m5 - e_1 * fs_1_13_3 * h7_m3 + e_1 * fs_1_26_2 * r_2 * h5_m5 + e_1 * fs_3_26_10 * r_2 * h5_m3;

        pc_11[k] = - e_0 * fs_1_4_110 * h5_m4 - e_1 * fs_1_13_39 * h7_m6 - e_1 * fs_1_26_6 * h7_m4 + e_1 * fs_1_26_110 * r_2 * h5_m4;

        pc_12[k] = - e_0 * fs_1_2_33 * h5_m5 - e_1 * fs_1_26_182 * h7_m7 - e_1 * fs_1_26_2 * h7_m5 + e_1 * fs_1_13_33 * r_2 * h5_m5;

        pc_13[k] = - e_1 * fs_1_13_13 * h7_m6;

        pc_14[k] = e_0 * fs_1_2_11 * h5_m5 - e_1 * fs_2_13_6 * h7_m5 - e_1 * fs_1_13_11 * r_2 * h5_m5;

        pc_15[k] = e_0 * fs_1_1_5 * h5_m4 - e_1 * fs_1_13_33 * h7_m4 - e_1 * fs_2_13_5 * r_2 * h5_m4;

        pc_16[k] = e_0 * fs_3_2_3 * h5_m3 - e_1 * fs_2_13_10 * h7_m3 - e_1 * fs_3_13_3 * r_2 * h5_m3;
    }

    // NOTE: the angular components are formed in 8 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph5_m2, ph5_m1, ph5_0, ph5_p1, ph5_p2, ph7_m2, ph7_m1, ph7_0, ph7_p1, ph7_p2, ab_2, pc_17, pc_18, pc_19, pc_20, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_17[k] = e_0 * fs_2_1_2 * h5_m2 - e_1 * fs_3_13_5 * h7_m2 - e_1 * fs_4_13_2 * r_2 * h5_m2;

        pc_18[k] = e_0 * fs_1_2_35 * h5_m1 - e_1 * fs_4_13_3 * h7_m1 - e_1 * fs_1_13_35 * r_2 * h5_m1;

        pc_19[k] = e_0 * f_3_1 * h5_0 - e_1 * f_7_13 * h7_0 - e_1 * f_6_13 * r_2 * h5_0;

        pc_20[k] = e_0 * fs_1_2_35 * h5_p1 - e_1 * fs_4_13_3 * h7_p1 - e_1 * fs_1_13_35 * r_2 * h5_p1;

        pc_21[k] = e_0 * fs_2_1_2 * h5_p2 - e_1 * fs_3_13_5 * h7_p2 - e_1 * fs_4_13_2 * r_2 * h5_p2;
    }

    // NOTE: the angular components are formed in 8 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph5_m5, ph5_p3, ph5_p4, ph5_p5, ph7_m7, ph7_m5, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ab_2, pc_22, pc_23, pc_24, pc_25, pc_26 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m5 = ph5_m5[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_22[k] = e_0 * fs_3_2_3 * h5_p3 - e_1 * fs_2_13_10 * h7_p3 - e_1 * fs_3_13_3 * r_2 * h5_p3;

        pc_23[k] = e_0 * fs_1_1_5 * h5_p4 - e_1 * fs_1_13_33 * h7_p4 - e_1 * fs_2_13_5 * r_2 * h5_p4;

        pc_24[k] = e_0 * fs_1_2_11 * h5_p5 - e_1 * fs_2_13_6 * h7_p5 - e_1 * fs_1_13_11 * r_2 * h5_p5;

        pc_25[k] = - e_1 * fs_1_13_13 * h7_p6;

        pc_26[k] = e_0 * fs_1_2_33 * h5_m5 - e_1 * fs_1_26_182 * h7_m7 + e_1 * fs_1_26_2 * h7_m5 - e_1 * fs_1_13_33 * r_2 * h5_m5;
    }

    // NOTE: the angular components are formed in 8 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ab_2, pc_27, pc_28, pc_29, pc_30, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_27[k] = e_0 * fs_1_4_110 * h5_m4 - e_1 * fs_1_13_39 * h7_m6 + e_1 * fs_1_26_6 * h7_m4 - e_1 * fs_1_26_110 * r_2 * h5_m4;

        pc_28[k] = - e_0 * fs_1_4_2 * h5_m5 + e_0 * fs_3_4_10 * h5_m3 - e_1 * fs_1_13_33 * h7_m5 + e_1 * fs_1_13_3 * h7_m3 + e_1 * fs_1_26_2 * r_2 * h5_m5 - e_1 * fs_3_26_10 * r_2 * h5_m3;

        pc_29[k] = - e_0 * fs_1_4_6 * h5_m4 + e_0 * fs_3_2_2 * h5_m2 - e_1 * fs_1_26_110 * h7_m4 + e_1 * fs_1_13_5 * h7_m2 + e_1 * fs_1_26_6 * r_2 * h5_m4 - e_1 * fs_3_13_2 * r_2 * h5_m2;

        pc_30[k] = - e_0 * fs_1_2_3 * h5_m3 + e_0 * fs_1_2_14 * h5_m1 - e_1 * fs_3_26_10 * h7_m3 + e_1 * fs_1_26_30 * h7_m1 + e_1 * fs_1_13_3 * r_2 * h5_m3 - e_1 * fs_1_13_14 * r_2 * h5_m1;

        pc_31[k] = - e_0 * fs_1_2_5 * h5_m2 - e_1 * fs_3_13_2 * h7_m2 + e_1 * fs_1_13_5 * r_2 * h5_m2;
    }

    // NOTE: the angular components are formed in 8 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ph7_0, ph7_p1, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ab_2, pc_32, pc_33, pc_34, pc_35, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_32[k] = - e_0 * fs_1_2_15 * h5_p1 - e_1 * fs_2_13_7 * h7_p1 + e_1 * fs_1_13_15 * r_2 * h5_p1;

        pc_33[k] = e_0 * fs_1_2_21 * h5_0 - e_0 * fs_1_2_5 * h5_p2 + e_1 * fs_1_13_21 * h7_0 - e_1 * fs_3_13_2 * h7_p2 - e_1 * fs_1_13_21 * r_2 * h5_0 + e_1 * fs_1_13_5 * r_2 * h5_p2;

        pc_34[k] = e_0 * fs_1_2_14 * h5_p1 - e_0 * fs_1_2_3 * h5_p3 + e_1 * fs_1_26_30 * h7_p1 - e_1 * fs_3_26_10 * h7_p3 - e_1 * fs_1_13_14 * r_2 * h5_p1 + e_1 * fs_1_13_3 * r_2 * h5_p3;

        pc_35[k] = e_0 * fs_3_2_2 * h5_p2 - e_0 * fs_1_4_6 * h5_p4 + e_1 * fs_1_13_5 * h7_p2 - e_1 * fs_1_26_110 * h7_p4 - e_1 * fs_3_13_2 * r_2 * h5_p2 + e_1 * fs_1_26_6 * r_2 * h5_p4;

        pc_36[k] = e_0 * fs_3_4_10 * h5_p3 - e_0 * fs_1_4_2 * h5_p5 + e_1 * fs_1_13_3 * h7_p3 - e_1 * fs_1_13_33 * h7_p5 - e_1 * fs_3_26_10 * r_2 * h5_p3 + e_1 * fs_1_26_2 * r_2 * h5_p5;
    }

    // NOTE: the angular components are formed in 8 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph5_p4, ph5_p5, ph7_p4, ph7_p5, ph7_p6, ph7_p7, ab_2, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_37[k] = e_0 * fs_1_4_110 * h5_p4 + e_1 * fs_1_26_6 * h7_p4 - e_1 * fs_1_13_39 * h7_p6 - e_1 * fs_1_26_110 * r_2 * h5_p4;

        pc_38[k] = e_0 * fs_1_2_33 * h5_p5 + e_1 * fs_1_26_2 * h7_p5 - e_1 * fs_1_26_182 * h7_p7 - e_1 * fs_1_13_33 * r_2 * h5_p5;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[39] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};

    for (size_t n = 0; n < 39; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
