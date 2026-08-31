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



#include "SimdKineticEnergyRecDF.hpp"

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
compute_df_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_df_kinetic_energy: Basis functions must be of angular momenta two and three"));
    }

    if (harmonics.size() < 5)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_df_kinetic_energy: Harmonics must reach angular momentum five"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_df_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 35 * nvalues, 0.0);

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

            const auto ff_0 = fbase * aexp * fmu / (fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * fmu * fmu / (fexp * fexp * fexp);

            const auto ff_2 = fbase * aexp * fmu * fmu * fmu / (fexp * fexp * fexp);

            const auto ff_3 = fbase * aexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp);

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_130_21 = 130.0 / 21.0;
    const auto f_18_35 = 18.0 / 35.0;
    const auto f_20_21 = 20.0 / 21.0;
    const auto f_22_3 = 22.0 / 3.0;
    const auto f_22_5 = 22.0 / 5.0;
    const auto f_243_10 = 243.0 / 10.0;
    const auto f_243_35 = 243.0 / 35.0;
    const auto f_2_3 = 2.0 / 3.0;
    const auto f_2_5 = 2.0 / 5.0;
    const auto f_33_2 = 33.0 / 2.0;
    const auto f_66_5 = 66.0 / 5.0;
    const auto f_81_4 = 81.0 / 4.0;
    const auto f_88_15 = 88.0 / 15.0;
    const auto f_8_15 = 8.0 / 15.0;
    const auto f_99_10 = 99.0 / 10.0;
    const auto fs_10_21_2 = std::sqrt(200.0 / 441.0);
    const auto fs_11_15_30 = std::sqrt(242.0 / 15.0);
    const auto fs_11_3_2 = std::sqrt(242.0 / 9.0);
    const auto fs_12_35_2 = std::sqrt(288.0 / 1225.0);
    const auto fs_13_21_105 = std::sqrt(845.0 / 21.0);
    const auto fs_13_21_35 = std::sqrt(845.0 / 63.0);
    const auto fs_13_21_42 = std::sqrt(338.0 / 21.0);
    const auto fs_13_21_5 = std::sqrt(845.0 / 441.0);
    const auto fs_13_42_14 = std::sqrt(169.0 / 126.0);
    const auto fs_13_42_2 = std::sqrt(169.0 / 882.0);
    const auto fs_13_42_210 = std::sqrt(845.0 / 42.0);
    const auto fs_13_42_30 = std::sqrt(845.0 / 294.0);
    const auto fs_13_7_10 = std::sqrt(1690.0 / 49.0);
    const auto fs_13_7_7 = std::sqrt(169.0 / 7.0);
    const auto fs_162_35_2 = std::sqrt(52488.0 / 1225.0);
    const auto fs_1_15_30 = std::sqrt(2.0 / 15.0);
    const auto fs_1_21_14 = std::sqrt(2.0 / 63.0);
    const auto fs_1_21_2 = std::sqrt(2.0 / 441.0);
    const auto fs_1_21_210 = std::sqrt(10.0 / 21.0);
    const auto fs_1_21_30 = std::sqrt(10.0 / 147.0);
    const auto fs_1_3_2 = std::sqrt(2.0 / 9.0);
    const auto fs_22_15_2 = std::sqrt(968.0 / 225.0);
    const auto fs_22_15_5 = std::sqrt(484.0 / 45.0);
    const auto fs_26_21_14 = std::sqrt(1352.0 / 63.0);
    const auto fs_26_21_3 = std::sqrt(676.0 / 147.0);
    const auto fs_26_21_7 = std::sqrt(676.0 / 63.0);
    const auto fs_27_2_2 = std::sqrt(729.0 / 2.0);
    const auto fs_27_4_3 = std::sqrt(2187.0 / 16.0);
    const auto fs_27_4_5 = std::sqrt(3645.0 / 16.0);
    const auto fs_27_4_6 = std::sqrt(2187.0 / 8.0);
    const auto fs_27_8_2 = std::sqrt(729.0 / 32.0);
    const auto fs_27_8_30 = std::sqrt(10935.0 / 32.0);
    const auto fs_2_15_2 = std::sqrt(8.0 / 225.0);
    const auto fs_2_15_5 = std::sqrt(4.0 / 45.0);
    const auto fs_2_21_105 = std::sqrt(20.0 / 21.0);
    const auto fs_2_21_35 = std::sqrt(20.0 / 63.0);
    const auto fs_2_21_42 = std::sqrt(8.0 / 21.0);
    const auto fs_2_21_5 = std::sqrt(20.0 / 441.0);
    const auto fs_2_7_10 = std::sqrt(40.0 / 49.0);
    const auto fs_2_7_7 = std::sqrt(4.0 / 7.0);
    const auto fs_33_10_2 = std::sqrt(1089.0 / 50.0);
    const auto fs_33_10_5 = std::sqrt(1089.0 / 20.0);
    const auto fs_33_20_30 = std::sqrt(3267.0 / 40.0);
    const auto fs_33_4_2 = std::sqrt(1089.0 / 8.0);
    const auto fs_33_5_3 = std::sqrt(3267.0 / 25.0);
    const auto fs_33_5_5 = std::sqrt(1089.0 / 5.0);
    const auto fs_3_35_2 = std::sqrt(18.0 / 1225.0);
    const auto fs_3_35_30 = std::sqrt(54.0 / 245.0);
    const auto fs_44_15_3 = std::sqrt(1936.0 / 75.0);
    const auto fs_44_15_5 = std::sqrt(1936.0 / 45.0);
    const auto fs_4_15_3 = std::sqrt(16.0 / 75.0);
    const auto fs_4_15_5 = std::sqrt(16.0 / 45.0);
    const auto fs_4_21_14 = std::sqrt(32.0 / 63.0);
    const auto fs_4_21_3 = std::sqrt(16.0 / 147.0);
    const auto fs_4_21_7 = std::sqrt(16.0 / 63.0);
    const auto fs_52_21_5 = std::sqrt(13520.0 / 441.0);
    const auto fs_65_21_2 = std::sqrt(8450.0 / 441.0);
    const auto fs_6_35_3 = std::sqrt(108.0 / 1225.0);
    const auto fs_6_35_5 = std::sqrt(36.0 / 245.0);
    const auto fs_6_35_6 = std::sqrt(216.0 / 1225.0);
    const auto fs_81_10_3 = std::sqrt(19683.0 / 100.0);
    const auto fs_81_10_5 = std::sqrt(6561.0 / 20.0);
    const auto fs_81_10_6 = std::sqrt(19683.0 / 50.0);
    const auto fs_81_20_2 = std::sqrt(6561.0 / 200.0);
    const auto fs_81_20_30 = std::sqrt(19683.0 / 40.0);
    const auto fs_81_35_3 = std::sqrt(19683.0 / 1225.0);
    const auto fs_81_35_5 = std::sqrt(6561.0 / 245.0);
    const auto fs_81_35_6 = std::sqrt(39366.0 / 1225.0);
    const auto fs_81_5_2 = std::sqrt(13122.0 / 25.0);
    const auto fs_81_70_2 = std::sqrt(6561.0 / 2450.0);
    const auto fs_81_70_30 = std::sqrt(19683.0 / 490.0);
    const auto fs_8_21_5 = std::sqrt(320.0 / 441.0);

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_0, ph1_p1, ph3_m2, ph3_0, ph3_p1, ph3_p3, ph5_m2, ph5_0, ph5_p1, ph5_p3, ph5_p4, ph5_p5, ab_2, pc_0, pc_1, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];

        pc_0[k] = e_0 * fs_27_8_30 * h1_p1 + e_1 * fs_33_10_5 * h3_p1 - e_1 * fs_81_20_30 * r_2 * h1_p1 + e_2 * fs_13_42_2 * h5_p1 - e_2 * fs_13_21_105 * h5_p5 - e_2 * fs_22_15_5 * r_2 * h3_p1 + e_2 * fs_81_70_30 * r_4 * h1_p1 - e_3 * fs_1_21_2 * r_2 * h5_p1 + e_3 * fs_2_21_105 * r_2 * h5_p5 + e_3 * fs_2_15_5 * r_4 * h3_p1 - e_3 * fs_3_35_30 * r_6 * h1_p1;

        pc_1[k] = e_0 * fs_27_4_5 * h1_0 + e_1 * fs_33_5_5 * h3_0 - e_1 * fs_81_10_5 * r_2 * h1_0 + e_2 * fs_13_21_5 * h5_0 - e_2 * fs_13_7_7 * h5_p4 - e_2 * fs_44_15_5 * r_2 * h3_0 + e_2 * fs_81_35_5 * r_4 * h1_0 - e_3 * fs_2_21_5 * r_2 * h5_0 + e_3 * fs_2_7_7 * r_2 * h5_p4 + e_3 * fs_4_15_5 * r_4 * h3_0 - e_3 * fs_6_35_5 * r_6 * h1_0;

        pc_2[k] = - e_0 * fs_27_8_2 * h1_p1 - e_1 * fs_33_5_3 * h3_p1 - e_1 * fs_33_10_5 * h3_p3 + e_1 * fs_81_20_2 * r_2 * h1_p1 - e_2 * fs_13_42_30 * h5_p1 - e_2 * fs_13_21_35 * h5_p3 + e_2 * fs_44_15_3 * r_2 * h3_p1 + e_2 * fs_22_15_5 * r_2 * h3_p3 - e_2 * fs_81_70_2 * r_4 * h1_p1 + e_3 * fs_1_21_30 * r_2 * h5_p1 + e_3 * fs_2_21_35 * r_2 * h5_p3 - e_3 * fs_4_15_3 * r_4 * h3_p1 - e_3 * fs_2_15_5 * r_4 * h3_p3 + e_3 * fs_3_35_2 * r_6 * h1_p1;

        pc_3[k] = e_1 * fs_33_5_5 * h3_m2 + e_2 * fs_13_21_35 * h5_m2 - e_2 * fs_44_15_5 * r_2 * h3_m2 - e_3 * fs_2_21_35 * r_2 * h5_m2 + e_3 * fs_4_15_5 * r_4 * h3_m2;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph3_m3, ph3_m1, ph3_p2, ph5_m5, ph5_m4, ph5_m3, ph5_m1, ph5_p2, ph5_p4, ab_2, pc_4, pc_5, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];

        pc_4[k] = - e_0 * fs_27_8_2 * h1_m1 + e_1 * fs_33_10_5 * h3_m3 - e_1 * fs_33_5_3 * h3_m1 + e_1 * fs_81_20_2 * r_2 * h1_m1 + e_2 * fs_13_21_35 * h5_m3 - e_2 * fs_13_42_30 * h5_m1 - e_2 * fs_22_15_5 * r_2 * h3_m3 + e_2 * fs_44_15_3 * r_2 * h3_m1 - e_2 * fs_81_70_2 * r_4 * h1_m1 - e_3 * fs_2_21_35 * r_2 * h5_m3 + e_3 * fs_1_21_30 * r_2 * h5_m1 + e_3 * fs_2_15_5 * r_4 * h3_m3 - e_3 * fs_4_15_3 * r_4 * h3_m1 + e_3 * fs_3_35_2 * r_6 * h1_m1;

        pc_5[k] = e_2 * fs_13_7_7 * h5_m4 - e_3 * fs_2_7_7 * r_2 * h5_m4;

        pc_6[k] = - e_0 * fs_27_8_30 * h1_m1 - e_1 * fs_33_10_5 * h3_m1 + e_1 * fs_81_20_30 * r_2 * h1_m1 + e_2 * fs_13_21_105 * h5_m5 - e_2 * fs_13_42_2 * h5_m1 + e_2 * fs_22_15_5 * r_2 * h3_m1 - e_2 * fs_81_70_30 * r_4 * h1_m1 - e_3 * fs_2_21_105 * r_2 * h5_m5 + e_3 * fs_1_21_2 * r_2 * h5_m1 - e_3 * fs_2_15_5 * r_4 * h3_m1 + e_3 * fs_3_35_30 * r_6 * h1_m1;

        pc_7[k] = - e_1 * fs_33_4_2 * h3_p2 - e_2 * fs_13_42_14 * h5_p2 - e_2 * fs_13_21_42 * h5_p4 + e_2 * fs_11_3_2 * r_2 * h3_p2 + e_3 * fs_1_21_14 * r_2 * h5_p2 + e_3 * fs_2_21_42 * r_2 * h5_p4 - e_3 * fs_1_3_2 * r_4 * h3_p2;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph1_0, ph1_p1, ph3_m1, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_m1, ph5_0, ph5_p1, ph5_p2, ph5_p3, ab_2, pc_8, pc_9, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];

        pc_8[k] = e_0 * fs_27_4_5 * h1_p1 - e_1 * fs_33_20_30 * h3_p1 + e_1 * fs_33_4_2 * h3_p3 - e_1 * fs_81_10_5 * r_2 * h1_p1 - e_2 * fs_26_21_3 * h5_p1 - e_2 * fs_26_21_14 * h5_p3 + e_2 * fs_11_15_30 * r_2 * h3_p1 - e_2 * fs_11_3_2 * r_2 * h3_p3 + e_2 * fs_81_35_5 * r_4 * h1_p1 + e_3 * fs_4_21_3 * r_2 * h5_p1 + e_3 * fs_4_21_14 * r_2 * h5_p3 - e_3 * fs_1_15_30 * r_4 * h3_p1 + e_3 * fs_1_3_2 * r_4 * h3_p3 - e_3 * fs_6_35_5 * r_6 * h1_p1;

        pc_9[k] = e_0 * fs_27_2_2 * h1_0 - e_1 * fs_33_10_2 * h3_0 + e_1 * fs_33_20_30 * h3_p2 - e_1 * fs_81_5_2 * r_2 * h1_0 - e_2 * fs_65_21_2 * h5_0 - e_2 * fs_13_42_210 * h5_p2 + e_2 * fs_22_15_2 * r_2 * h3_0 - e_2 * fs_11_15_30 * r_2 * h3_p2 + e_2 * fs_162_35_2 * r_4 * h1_0 + e_3 * fs_10_21_2 * r_2 * h5_0 + e_3 * fs_1_21_210 * r_2 * h5_p2 - e_3 * fs_2_15_2 * r_4 * h3_0 + e_3 * fs_1_15_30 * r_4 * h3_p2 - e_3 * fs_12_35_2 * r_6 * h1_0;

        pc_10[k] = - e_0 * fs_27_4_3 * h1_m1 - e_1 * fs_33_10_2 * h3_m1 + e_1 * fs_81_10_3 * r_2 * h1_m1 + e_2 * fs_52_21_5 * h5_m1 + e_2 * fs_22_15_2 * r_2 * h3_m1 - e_2 * fs_81_35_3 * r_4 * h1_m1 - e_3 * fs_8_21_5 * r_2 * h5_m1 - e_3 * fs_2_15_2 * r_4 * h3_m1 + e_3 * fs_6_35_3 * r_6 * h1_m1;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph1_0, ph3_m3, ph3_m2, ph3_m1, ph3_0, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph5_0, ab_2, pc_11, pc_12, pc_13, pc_14, pc_15, pc_16, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];

        pc_11[k] = - e_1 * fs_33_20_30 * h3_m2 + e_2 * fs_13_42_210 * h5_m2 + e_2 * fs_11_15_30 * r_2 * h3_m2 - e_3 * fs_1_21_210 * r_2 * h5_m2 - e_3 * fs_1_15_30 * r_4 * h3_m2;

        pc_12[k] = - e_0 * fs_27_4_5 * h1_m1 - e_1 * fs_33_4_2 * h3_m3 + e_1 * fs_33_20_30 * h3_m1 + e_1 * fs_81_10_5 * r_2 * h1_m1 + e_2 * fs_26_21_14 * h5_m3 + e_2 * fs_26_21_3 * h5_m1 + e_2 * fs_11_3_2 * r_2 * h3_m3 - e_2 * fs_11_15_30 * r_2 * h3_m1 - e_2 * fs_81_35_5 * r_4 * h1_m1 - e_3 * fs_4_21_14 * r_2 * h5_m3 - e_3 * fs_4_21_3 * r_2 * h5_m1 - e_3 * fs_1_3_2 * r_4 * h3_m3 + e_3 * fs_1_15_30 * r_4 * h3_m1 + e_3 * fs_6_35_5 * r_6 * h1_m1;

        pc_13[k] = e_1 * fs_33_4_2 * h3_m2 + e_2 * fs_13_21_42 * h5_m4 + e_2 * fs_13_42_14 * h5_m2 - e_2 * fs_11_3_2 * r_2 * h3_m2 - e_3 * fs_2_21_42 * r_2 * h5_m4 - e_3 * fs_1_21_14 * r_2 * h5_m2 + e_3 * fs_1_3_2 * r_4 * h3_m2;

        pc_14[k] = e_1 * f_33_2 * h3_m3 + e_2 * fs_26_21_7 * h5_m3 - e_2 * f_22_3 * r_2 * h3_m3 - e_3 * fs_4_21_7 * r_2 * h5_m3 + e_3 * f_2_3 * r_4 * h3_m3;

        pc_15[k] = e_2 * fs_13_7_7 * h5_m2 - e_3 * fs_2_7_7 * r_2 * h5_m2;

        pc_16[k] = e_0 * fs_27_4_6 * h1_m1 - e_1 * f_99_10 * h3_m1 - e_1 * fs_81_10_6 * r_2 * h1_m1 + e_2 * fs_13_7_10 * h5_m1 + e_2 * f_22_5 * r_2 * h3_m1 + e_2 * fs_81_35_6 * r_4 * h1_m1 - e_3 * fs_2_7_10 * r_2 * h5_m1 - e_3 * f_2_5 * r_4 * h3_m1 - e_3 * fs_6_35_6 * r_6 * h1_m1;

        pc_17[k] = e_0 * f_81_4 * h1_0 - e_1 * f_66_5 * h3_0 - e_1 * f_243_10 * r_2 * h1_0 + e_2 * f_130_21 * h5_0 + e_2 * f_88_15 * r_2 * h3_0 + e_2 * f_243_35 * r_4 * h1_0 - e_3 * f_20_21 * r_2 * h5_0 - e_3 * f_8_15 * r_4 * h3_0 - e_3 * f_18_35 * r_6 * h1_0;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_p1, ph3_m2, ph3_p1, ph3_p3, ph5_m4, ph5_m2, ph5_p1, ph5_p2, ph5_p3, ab_2, pc_18, pc_19, pc_20, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];

        pc_18[k] = e_0 * fs_27_4_6 * h1_p1 - e_1 * f_99_10 * h3_p1 - e_1 * fs_81_10_6 * r_2 * h1_p1 + e_2 * fs_13_7_10 * h5_p1 + e_2 * f_22_5 * r_2 * h3_p1 + e_2 * fs_81_35_6 * r_4 * h1_p1 - e_3 * fs_2_7_10 * r_2 * h5_p1 - e_3 * f_2_5 * r_4 * h3_p1 - e_3 * fs_6_35_6 * r_6 * h1_p1;

        pc_19[k] = e_2 * fs_13_7_7 * h5_p2 - e_3 * fs_2_7_7 * r_2 * h5_p2;

        pc_20[k] = e_1 * f_33_2 * h3_p3 + e_2 * fs_26_21_7 * h5_p3 - e_2 * f_22_3 * r_2 * h3_p3 - e_3 * fs_4_21_7 * r_2 * h5_p3 + e_3 * f_2_3 * r_4 * h3_p3;

        pc_21[k] = - e_1 * fs_33_4_2 * h3_m2 + e_2 * fs_13_21_42 * h5_m4 - e_2 * fs_13_42_14 * h5_m2 + e_2 * fs_11_3_2 * r_2 * h3_m2 - e_3 * fs_2_21_42 * r_2 * h5_m4 + e_3 * fs_1_21_14 * r_2 * h5_m2 - e_3 * fs_1_3_2 * r_4 * h3_m2;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph1_p1, ph3_m3, ph3_m2, ph3_m1, ph3_p1, ph5_m3, ph5_m2, ph5_m1, ph5_p1, ab_2, pc_22, pc_23, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];

        pc_22[k] = e_0 * fs_27_4_5 * h1_m1 - e_1 * fs_33_4_2 * h3_m3 - e_1 * fs_33_20_30 * h3_m1 - e_1 * fs_81_10_5 * r_2 * h1_m1 + e_2 * fs_26_21_14 * h5_m3 - e_2 * fs_26_21_3 * h5_m1 + e_2 * fs_11_3_2 * r_2 * h3_m3 + e_2 * fs_11_15_30 * r_2 * h3_m1 + e_2 * fs_81_35_5 * r_4 * h1_m1 - e_3 * fs_4_21_14 * r_2 * h5_m3 + e_3 * fs_4_21_3 * r_2 * h5_m1 - e_3 * fs_1_3_2 * r_4 * h3_m3 - e_3 * fs_1_15_30 * r_4 * h3_m1 - e_3 * fs_6_35_5 * r_6 * h1_m1;

        pc_23[k] = - e_1 * fs_33_20_30 * h3_m2 + e_2 * fs_13_42_210 * h5_m2 + e_2 * fs_11_15_30 * r_2 * h3_m2 - e_3 * fs_1_21_210 * r_2 * h5_m2 - e_3 * fs_1_15_30 * r_4 * h3_m2;

        pc_24[k] = - e_0 * fs_27_4_3 * h1_p1 - e_1 * fs_33_10_2 * h3_p1 + e_1 * fs_81_10_3 * r_2 * h1_p1 + e_2 * fs_52_21_5 * h5_p1 + e_2 * fs_22_15_2 * r_2 * h3_p1 - e_2 * fs_81_35_3 * r_4 * h1_p1 - e_3 * fs_8_21_5 * r_2 * h5_p1 - e_3 * fs_2_15_2 * r_4 * h3_p1 + e_3 * fs_6_35_3 * r_6 * h1_p1;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_0, ph1_p1, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ab_2, pc_25, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];

        pc_25[k] = e_0 * fs_27_2_2 * h1_0 - e_1 * fs_33_10_2 * h3_0 - e_1 * fs_33_20_30 * h3_p2 - e_1 * fs_81_5_2 * r_2 * h1_0 - e_2 * fs_65_21_2 * h5_0 + e_2 * fs_13_42_210 * h5_p2 + e_2 * fs_22_15_2 * r_2 * h3_0 + e_2 * fs_11_15_30 * r_2 * h3_p2 + e_2 * fs_162_35_2 * r_4 * h1_0 + e_3 * fs_10_21_2 * r_2 * h5_0 - e_3 * fs_1_21_210 * r_2 * h5_p2 - e_3 * fs_2_15_2 * r_4 * h3_0 - e_3 * fs_1_15_30 * r_4 * h3_p2 - e_3 * fs_12_35_2 * r_6 * h1_0;

        pc_26[k] = e_0 * fs_27_4_5 * h1_p1 - e_1 * fs_33_20_30 * h3_p1 - e_1 * fs_33_4_2 * h3_p3 - e_1 * fs_81_10_5 * r_2 * h1_p1 - e_2 * fs_26_21_3 * h5_p1 + e_2 * fs_26_21_14 * h5_p3 + e_2 * fs_11_15_30 * r_2 * h3_p1 + e_2 * fs_11_3_2 * r_2 * h3_p3 + e_2 * fs_81_35_5 * r_4 * h1_p1 + e_3 * fs_4_21_3 * r_2 * h5_p1 - e_3 * fs_4_21_14 * r_2 * h5_p3 - e_3 * fs_1_15_30 * r_4 * h3_p1 - e_3 * fs_1_3_2 * r_4 * h3_p3 - e_3 * fs_6_35_5 * r_6 * h1_p1;

        pc_27[k] = - e_1 * fs_33_4_2 * h3_p2 - e_2 * fs_13_42_14 * h5_p2 + e_2 * fs_13_21_42 * h5_p4 + e_2 * fs_11_3_2 * r_2 * h3_p2 + e_3 * fs_1_21_14 * r_2 * h5_p2 - e_3 * fs_2_21_42 * r_2 * h5_p4 - e_3 * fs_1_3_2 * r_4 * h3_p2;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph3_m3, ph3_m1, ph3_p2, ph5_m5, ph5_m4, ph5_m3, ph5_m1, ph5_p2, ab_2, pc_28, pc_29, pc_30, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];

        pc_28[k] = e_0 * fs_27_8_30 * h1_m1 + e_1 * fs_33_10_5 * h3_m1 - e_1 * fs_81_20_30 * r_2 * h1_m1 + e_2 * fs_13_21_105 * h5_m5 + e_2 * fs_13_42_2 * h5_m1 - e_2 * fs_22_15_5 * r_2 * h3_m1 + e_2 * fs_81_70_30 * r_4 * h1_m1 - e_3 * fs_2_21_105 * r_2 * h5_m5 - e_3 * fs_1_21_2 * r_2 * h5_m1 + e_3 * fs_2_15_5 * r_4 * h3_m1 - e_3 * fs_3_35_30 * r_6 * h1_m1;

        pc_29[k] = e_2 * fs_13_7_7 * h5_m4 - e_3 * fs_2_7_7 * r_2 * h5_m4;

        pc_30[k] = e_0 * fs_27_8_2 * h1_m1 + e_1 * fs_33_10_5 * h3_m3 + e_1 * fs_33_5_3 * h3_m1 - e_1 * fs_81_20_2 * r_2 * h1_m1 + e_2 * fs_13_21_35 * h5_m3 + e_2 * fs_13_42_30 * h5_m1 - e_2 * fs_22_15_5 * r_2 * h3_m3 - e_2 * fs_44_15_3 * r_2 * h3_m1 + e_2 * fs_81_70_2 * r_4 * h1_m1 - e_3 * fs_2_21_35 * r_2 * h5_m3 - e_3 * fs_1_21_30 * r_2 * h5_m1 + e_3 * fs_2_15_5 * r_4 * h3_m3 + e_3 * fs_4_15_3 * r_4 * h3_m1 - e_3 * fs_3_35_2 * r_6 * h1_m1;

        pc_31[k] = e_1 * fs_33_5_5 * h3_p2 + e_2 * fs_13_21_35 * h5_p2 - e_2 * fs_44_15_5 * r_2 * h3_p2 - e_3 * fs_2_21_35 * r_2 * h5_p2 + e_3 * fs_4_15_5 * r_4 * h3_p2;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_0, ph1_p1, ph3_0, ph3_p1, ph3_p3, ph5_0, ph5_p1, ph5_p3, ph5_p4, ph5_p5, ab_2, pc_32, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];

        pc_32[k] = - e_0 * fs_27_8_2 * h1_p1 - e_1 * fs_33_5_3 * h3_p1 + e_1 * fs_33_10_5 * h3_p3 + e_1 * fs_81_20_2 * r_2 * h1_p1 - e_2 * fs_13_42_30 * h5_p1 + e_2 * fs_13_21_35 * h5_p3 + e_2 * fs_44_15_3 * r_2 * h3_p1 - e_2 * fs_22_15_5 * r_2 * h3_p3 - e_2 * fs_81_70_2 * r_4 * h1_p1 + e_3 * fs_1_21_30 * r_2 * h5_p1 - e_3 * fs_2_21_35 * r_2 * h5_p3 - e_3 * fs_4_15_3 * r_4 * h3_p1 + e_3 * fs_2_15_5 * r_4 * h3_p3 + e_3 * fs_3_35_2 * r_6 * h1_p1;

        pc_33[k] = e_0 * fs_27_4_5 * h1_0 + e_1 * fs_33_5_5 * h3_0 - e_1 * fs_81_10_5 * r_2 * h1_0 + e_2 * fs_13_21_5 * h5_0 + e_2 * fs_13_7_7 * h5_p4 - e_2 * fs_44_15_5 * r_2 * h3_0 + e_2 * fs_81_35_5 * r_4 * h1_0 - e_3 * fs_2_21_5 * r_2 * h5_0 - e_3 * fs_2_7_7 * r_2 * h5_p4 + e_3 * fs_4_15_5 * r_4 * h3_0 - e_3 * fs_6_35_5 * r_6 * h1_0;

        pc_34[k] = e_0 * fs_27_8_30 * h1_p1 + e_1 * fs_33_10_5 * h3_p1 - e_1 * fs_81_20_30 * r_2 * h1_p1 + e_2 * fs_13_42_2 * h5_p1 + e_2 * fs_13_21_105 * h5_p5 - e_2 * fs_22_15_5 * r_2 * h3_p1 + e_2 * fs_81_70_30 * r_4 * h1_p1 - e_3 * fs_1_21_2 * r_2 * h5_p1 - e_3 * fs_2_21_105 * r_2 * h5_p5 + e_3 * fs_2_15_5 * r_4 * h3_p1 - e_3 * fs_3_35_30 * r_6 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[35] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};

    for (size_t n = 0; n < 35; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
