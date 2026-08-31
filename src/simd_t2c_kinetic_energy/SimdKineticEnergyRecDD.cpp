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



#include "SimdKineticEnergyRecDD.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_dd_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 2))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_dd_kinetic_energy: Basis functions must be of angular momenta two and two"));
    }

    if (harmonics.size() < 4)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_dd_kinetic_energy: Harmonics must reach angular momentum four"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_dd_kinetic_energy: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_kinetic_energy_primitive_bound,
        threshold / static_cast<double>(nprims));

    // NOTE: the buffer spans the atom pairs reached by the pair of primitives
    // reaching furthest, which is searched for rather than assumed. The
    // primitives are sorted by descending exponent, but the bound of a pair of
    // primitives carries their prefactor as well as their decay, so a tighter
    // pair with a larger prefactor reaches further than a more diffuse pair with
    // a smaller one, and the last pair is not always the furthest reaching.

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 25 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time.

    auto buffer = CSimdMatrix(4, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);

    const auto *ab_2 = coordinates.data(6);

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the prefactor of each exponent factor over the pairs of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto bexp = b_exps[j];

            const auto fexp = aexp + bexp;

            const auto fmu = aexp * bexp / fexp;

            const auto fovl = fpi / fexp;

            const auto fbase = anorm * b_norms[j] * fovl * std::sqrt(fovl);

            const auto ff_0 = fbase * fmu / (fexp * fexp);

            const auto ff_1 = fbase * fmu * fmu / (fexp * fexp);

            const auto ff_2 = fbase * fmu * fmu * fmu / (fexp * fexp);

            const auto ff_3 = fbase * fmu * fmu * fmu * fmu / (fexp * fexp);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fterm = std::exp(-fmu * ab_2[k]);

                pe_0[k] += ff_0 * fterm;
                pe_1[k] += ff_1 * fterm;
                pe_2[k] += ff_2 * fterm;
                pe_3[k] += ff_3 * fterm;
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

    auto *pc_0 = values + 0 * nvalues;
    auto *pc_1 = values + 1 * nvalues;
    auto *pc_2 = values + 2 * nvalues;
    auto *pc_3 = values + 3 * nvalues;
    auto *pc_4 = values + 4 * nvalues;
    auto *pc_5 = values + 6 * nvalues;
    auto *pc_6 = values + 7 * nvalues;
    auto *pc_7 = values + 8 * nvalues;
    auto *pc_8 = values + 9 * nvalues;
    auto *pc_9 = values + 12 * nvalues;
    auto *pc_10 = values + 13 * nvalues;
    auto *pc_11 = values + 14 * nvalues;
    auto *pc_12 = values + 18 * nvalues;
    auto *pc_13 = values + 19 * nvalues;
    auto *pc_14 = values + 24 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_132_35 = 132.0 / 35.0;
    const auto f_18_7 = 18.0 / 7.0;
    const auto f_198_35 = 198.0 / 35.0;
    const auto f_21_2 = 21.0 / 2.0;
    const auto f_21_4 = 21.0 / 4.0;
    const auto f_21_5 = 21.0 / 5.0;
    const auto f_24_35 = 24.0 / 35.0;
    const auto f_2_5 = 2.0 / 5.0;
    const auto f_2_7 = 2.0 / 7.0;
    const auto f_33_35 = 33.0 / 35.0;
    const auto f_36_35 = 36.0 / 35.0;
    const auto f_36_7 = 36.0 / 7.0;
    const auto f_4_7 = 4.0 / 7.0;
    const auto f_6_35 = 6.0 / 35.0;
    const auto f_9_1 = 9.0;
    const auto f_9_2 = 9.0 / 2.0;
    const auto fs_12_35_5 = std::sqrt(144.0 / 245.0);
    const auto fs_18_7_3 = std::sqrt(972.0 / 49.0);
    const auto fs_2_7_3 = std::sqrt(12.0 / 49.0);
    const auto fs_33_35_15 = std::sqrt(3267.0 / 245.0);
    const auto fs_33_35_30 = std::sqrt(6534.0 / 245.0);
    const auto fs_33_35_35 = std::sqrt(1089.0 / 35.0);
    const auto fs_33_70_10 = std::sqrt(1089.0 / 490.0);
    const auto fs_33_70_70 = std::sqrt(1089.0 / 70.0);
    const auto fs_3_35_10 = std::sqrt(18.0 / 245.0);
    const auto fs_3_35_70 = std::sqrt(18.0 / 35.0);
    const auto fs_66_35_5 = std::sqrt(4356.0 / 245.0);
    const auto fs_6_35_15 = std::sqrt(108.0 / 245.0);
    const auto fs_6_35_30 = std::sqrt(216.0 / 245.0);
    const auto fs_6_35_35 = std::sqrt(36.0 / 35.0);
    const auto fs_9_2_3 = std::sqrt(243.0 / 4.0);

    // NOTE: the angular components are formed in 2 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph2_0, ph2_p1, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph4_0, ph4_p1, ph4_p3, ph4_p4, ab_2, pc_0, pc_1, pc_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
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
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];

        pc_0[k] = e_0 * f_21_4 + e_1 * f_9_1 * h2_0 - e_1 * f_21_2 * r_2 + e_2 * f_33_35 * h4_0 - e_2 * fs_33_35_35 * h4_p4 - e_2 * f_36_7 * r_2 * h2_0 + e_2 * f_21_5 * r_4 - e_3 * f_6_35 * r_2 * h4_0 + e_3 * fs_6_35_35 * r_2 * h4_p4 + e_3 * f_4_7 * r_4 * h2_0 - e_3 * f_2_5 * r_6;

        pc_1[k] = - e_1 * fs_9_2_3 * h2_p1 - e_2 * fs_33_70_10 * h4_p1 - e_2 * fs_33_70_70 * h4_p3 + e_2 * fs_18_7_3 * r_2 * h2_p1 + e_3 * fs_3_35_10 * r_2 * h4_p1 + e_3 * fs_3_35_70 * r_2 * h4_p3 - e_3 * fs_2_7_3 * r_4 * h2_p1;

        pc_2[k] = e_1 * f_9_1 * h2_m2 + e_2 * fs_33_35_15 * h4_m2 - e_2 * f_36_7 * r_2 * h2_m2 - e_3 * fs_6_35_15 * r_2 * h4_m2 + e_3 * f_4_7 * r_4 * h2_m2;

        pc_3[k] = - e_1 * fs_9_2_3 * h2_m1 + e_2 * fs_33_70_70 * h4_m3 - e_2 * fs_33_70_10 * h4_m1 + e_2 * fs_18_7_3 * r_2 * h2_m1 - e_3 * fs_3_35_70 * r_2 * h4_m3 + e_3 * fs_3_35_10 * r_2 * h4_m1 - e_3 * fs_2_7_3 * r_4 * h2_m1;

        pc_4[k] = e_2 * fs_33_35_35 * h4_m4 - e_3 * fs_6_35_35 * r_2 * h4_m4;
    }

    // NOTE: the angular components are formed in 2 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph2_0, ph2_p1, ph2_p2, ph4_m3, ph4_m2, ph4_m1, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ab_2, pc_5, pc_6, pc_7, pc_8, pc_9, pc_10, pc_11, pc_12, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
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
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];

        pc_5[k] = e_0 * f_21_4 - e_1 * f_9_2 * h2_0 + e_1 * fs_9_2_3 * h2_p2 - e_1 * f_21_2 * r_2 - e_2 * f_132_35 * h4_0 - e_2 * fs_66_35_5 * h4_p2 + e_2 * f_18_7 * r_2 * h2_0 - e_2 * fs_18_7_3 * r_2 * h2_p2 + e_2 * f_21_5 * r_4 + e_3 * f_24_35 * r_2 * h4_0 + e_3 * fs_12_35_5 * r_2 * h4_p2 - e_3 * f_2_7 * r_4 * h2_0 + e_3 * fs_2_7_3 * r_4 * h2_p2 - e_3 * f_2_5 * r_6;

        pc_6[k] = - e_1 * f_9_2 * h2_m1 + e_2 * fs_33_35_30 * h4_m1 + e_2 * f_18_7 * r_2 * h2_m1 - e_3 * fs_6_35_30 * r_2 * h4_m1 - e_3 * f_2_7 * r_4 * h2_m1;

        pc_7[k] = - e_1 * fs_9_2_3 * h2_m2 + e_2 * fs_66_35_5 * h4_m2 + e_2 * fs_18_7_3 * r_2 * h2_m2 - e_3 * fs_12_35_5 * r_2 * h4_m2 - e_3 * fs_2_7_3 * r_4 * h2_m2;

        pc_8[k] = e_1 * fs_9_2_3 * h2_m1 + e_2 * fs_33_70_70 * h4_m3 + e_2 * fs_33_70_10 * h4_m1 - e_2 * fs_18_7_3 * r_2 * h2_m1 - e_3 * fs_3_35_70 * r_2 * h4_m3 - e_3 * fs_3_35_10 * r_2 * h4_m1 + e_3 * fs_2_7_3 * r_4 * h2_m1;

        pc_9[k] = e_0 * f_21_4 - e_1 * f_9_1 * h2_0 - e_1 * f_21_2 * r_2 + e_2 * f_198_35 * h4_0 + e_2 * f_36_7 * r_2 * h2_0 + e_2 * f_21_5 * r_4 - e_3 * f_36_35 * r_2 * h4_0 - e_3 * f_4_7 * r_4 * h2_0 - e_3 * f_2_5 * r_6;

        pc_10[k] = - e_1 * f_9_2 * h2_p1 + e_2 * fs_33_35_30 * h4_p1 + e_2 * f_18_7 * r_2 * h2_p1 - e_3 * fs_6_35_30 * r_2 * h4_p1 - e_3 * f_2_7 * r_4 * h2_p1;

        pc_11[k] = e_1 * f_9_1 * h2_p2 + e_2 * fs_33_35_15 * h4_p2 - e_2 * f_36_7 * r_2 * h2_p2 - e_3 * fs_6_35_15 * r_2 * h4_p2 + e_3 * f_4_7 * r_4 * h2_p2;

        pc_12[k] = e_0 * f_21_4 - e_1 * f_9_2 * h2_0 - e_1 * fs_9_2_3 * h2_p2 - e_1 * f_21_2 * r_2 - e_2 * f_132_35 * h4_0 + e_2 * fs_66_35_5 * h4_p2 + e_2 * f_18_7 * r_2 * h2_0 + e_2 * fs_18_7_3 * r_2 * h2_p2 + e_2 * f_21_5 * r_4 + e_3 * f_24_35 * r_2 * h4_0 - e_3 * fs_12_35_5 * r_2 * h4_p2 - e_3 * f_2_7 * r_4 * h2_0 - e_3 * fs_2_7_3 * r_4 * h2_p2 - e_3 * f_2_5 * r_6;

        pc_13[k] = - e_1 * fs_9_2_3 * h2_p1 - e_2 * fs_33_70_10 * h4_p1 + e_2 * fs_33_70_70 * h4_p3 + e_2 * fs_18_7_3 * r_2 * h2_p1 + e_3 * fs_3_35_10 * r_2 * h4_p1 - e_3 * fs_3_35_70 * r_2 * h4_p3 - e_3 * fs_2_7_3 * r_4 * h2_p1;

        pc_14[k] = e_0 * f_21_4 + e_1 * f_9_1 * h2_0 - e_1 * f_21_2 * r_2 + e_2 * f_33_35 * h4_0 + e_2 * fs_33_35_35 * h4_p4 - e_2 * f_36_7 * r_2 * h2_0 + e_2 * f_21_5 * r_4 - e_3 * f_6_35 * r_2 * h4_0 - e_3 * fs_6_35_35 * r_2 * h4_p4 + e_3 * f_4_7 * r_4 * h2_0 - e_3 * f_2_5 * r_6;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[25] = {0, 1, 2, 3, 4, 1, 6, 7, 8, 9, 2, 7, 12, 13, 14, 3, 8, 13, 18, 19, 4, 9, 14, 19, 24};

    for (size_t n = 0; n < 25; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
