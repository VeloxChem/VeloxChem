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



#include "SimdTwoCenterElectronRepulsionRecDP.hpp"

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
compute_dp_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 1))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_dp_electron_repulsion: Basis functions must be of angular momenta two and one"));
    }

    if (harmonics.size() < 3)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_dp_electron_repulsion: Harmonics must reach angular momentum three"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_dp_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 5);

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
    // call, which fills the orders 0 to 3 of every row. The terms read the
    // orders 2 to 3 alone, and the orders below them are formed on the
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

            const auto ff_0 = fbase * bexp / (fexp * fexp);

            const auto ff_1 = fbase * bexp * fmu / (fexp * fexp);

            const auto *bv_0 = boys.data(3, i * nprim_b + j);

            const auto *bv_1 = boys.data(4, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, bv_0, bv_1 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_1_2 = 1.0 / 2.0;
    const auto f_1_5 = 1.0 / 5.0;
    const auto f_2_5 = 2.0 / 5.0;
    const auto f_3_5 = 3.0 / 5.0;
    const auto fs_1_10_2 = std::sqrt(1.0 / 50.0);
    const auto fs_1_10_30 = std::sqrt(3.0 / 10.0);
    const auto fs_1_2_3 = std::sqrt(3.0 / 4.0);
    const auto fs_1_5_3 = std::sqrt(3.0 / 25.0);
    const auto fs_1_5_5 = std::sqrt(1.0 / 5.0);
    const auto fs_1_5_6 = std::sqrt(6.0 / 25.0);
    const auto fs_2_5_2 = std::sqrt(8.0 / 25.0);

    // NOTE: the angular components are formed in 3 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph1_m1, ph1_0, ph1_p1, ph3_m3, ph3_m2, ph3_m1, ph3_0, ph3_p1, ph3_p2, ph3_p3, ab_2, pc_0, pc_1, pc_2, pc_3, pc_4, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];

        pc_0[k] = - e_0 * fs_1_2_3 * h1_p1 - e_1 * fs_1_10_2 * h3_p1 - e_1 * fs_1_10_30 * h3_p3 + e_1 * fs_1_5_3 * r_2 * h1_p1;

        pc_1[k] = e_1 * fs_1_5_5 * h3_m2;

        pc_2[k] = - e_0 * fs_1_2_3 * h1_m1 + e_1 * fs_1_10_30 * h3_m3 - e_1 * fs_1_10_2 * h3_m1 + e_1 * fs_1_5_3 * r_2 * h1_m1;

        pc_3[k] = - e_0 * fs_1_2_3 * h1_0 - e_1 * fs_1_5_3 * h3_0 - e_1 * fs_1_5_5 * h3_p2 + e_1 * fs_1_5_3 * r_2 * h1_0;

        pc_4[k] = - e_0 * fs_1_2_3 * h1_m1 + e_1 * fs_2_5_2 * h3_m1 + e_1 * fs_1_5_3 * r_2 * h1_m1;

        pc_5[k] = e_1 * fs_1_5_5 * h3_m2;

        pc_6[k] = e_0 * f_1_2 * h1_m1 + e_1 * fs_1_5_6 * h3_m1 - e_1 * f_1_5 * r_2 * h1_m1;
    }

    // NOTE: the angular components are formed in 3 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph1_m1, ph1_0, ph1_p1, ph3_m3, ph3_m2, ph3_m1, ph3_0, ph3_p1, ph3_p2, ab_2, pc_7, pc_8, pc_9, pc_10, pc_11, pc_12, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];

        pc_7[k] = - e_0 * h1_0 + e_1 * f_3_5 * h3_0 + e_1 * f_2_5 * r_2 * h1_0;

        pc_8[k] = e_0 * f_1_2 * h1_p1 + e_1 * fs_1_5_6 * h3_p1 - e_1 * f_1_5 * r_2 * h1_p1;

        pc_9[k] = e_1 * fs_1_5_5 * h3_m2;

        pc_10[k] = - e_0 * fs_1_2_3 * h1_p1 + e_1 * fs_2_5_2 * h3_p1 + e_1 * fs_1_5_3 * r_2 * h1_p1;

        pc_11[k] = - e_0 * fs_1_2_3 * h1_0 - e_1 * fs_1_5_3 * h3_0 + e_1 * fs_1_5_5 * h3_p2 + e_1 * fs_1_5_3 * r_2 * h1_0;

        pc_12[k] = e_0 * fs_1_2_3 * h1_m1 + e_1 * fs_1_10_30 * h3_m3 + e_1 * fs_1_10_2 * h3_m1 - e_1 * fs_1_5_3 * r_2 * h1_m1;

        pc_13[k] = e_1 * fs_1_5_5 * h3_p2;
    }

    // NOTE: the angular components are formed in 3 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, ph1_p1, ph3_p1, ph3_p3, ab_2, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];

        pc_14[k] = - e_0 * fs_1_2_3 * h1_p1 - e_1 * fs_1_10_2 * h3_p1 + e_1 * fs_1_10_30 * h3_p3 + e_1 * fs_1_5_3 * r_2 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

    for (size_t n = 0; n < 15; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
