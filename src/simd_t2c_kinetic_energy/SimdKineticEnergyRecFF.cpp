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



#include "SimdKineticEnergyRecFF.hpp"

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
compute_ff_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 3) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ff_kinetic_energy: Basis functions must be of angular momenta three and three"));
    }

    if (harmonics.size() < 6)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ff_kinetic_energy: Harmonics must reach angular momentum six"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ff_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 49 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time.

    auto buffer = CSimdMatrix(5, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);
    std::fill(pe_4, pe_4 + nmax, 0.0);

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

            const auto ff_0 = fbase * fmu / (fexp * fexp * fexp);

            const auto ff_1 = fbase * fmu * fmu / (fexp * fexp * fexp);

            const auto ff_2 = fbase * fmu * fmu * fmu / (fexp * fexp * fexp);

            const auto ff_3 = fbase * fmu * fmu * fmu * fmu / (fexp * fexp * fexp);

            const auto ff_4 = fbase * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fterm = std::exp(-fmu * ab_2[k]);

                pe_0[k] += ff_0 * fterm;
                pe_1[k] += ff_1 * fterm;
                pe_2[k] += ff_2 * fterm;
                pe_3[k] += ff_3 * fterm;
                pe_4[k] += ff_4 * fterm;
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
    auto *pc_7 = values + 8 * nvalues;
    auto *pc_8 = values + 9 * nvalues;
    auto *pc_9 = values + 10 * nvalues;
    auto *pc_10 = values + 11 * nvalues;
    auto *pc_11 = values + 12 * nvalues;
    auto *pc_12 = values + 13 * nvalues;
    auto *pc_13 = values + 16 * nvalues;
    auto *pc_14 = values + 17 * nvalues;
    auto *pc_15 = values + 18 * nvalues;
    auto *pc_16 = values + 19 * nvalues;
    auto *pc_17 = values + 20 * nvalues;
    auto *pc_18 = values + 24 * nvalues;
    auto *pc_19 = values + 25 * nvalues;
    auto *pc_20 = values + 26 * nvalues;
    auto *pc_21 = values + 27 * nvalues;
    auto *pc_22 = values + 32 * nvalues;
    auto *pc_23 = values + 33 * nvalues;
    auto *pc_24 = values + 34 * nvalues;
    auto *pc_25 = values + 40 * nvalues;
    auto *pc_26 = values + 41 * nvalues;
    auto *pc_27 = values + 48 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_10_21 = 10.0 / 21.0;
    const auto f_10_231 = 10.0 / 231.0;
    const auto f_117_14 = 117.0 / 14.0;
    const auto f_117_7 = 117.0 / 7.0;
    const auto f_135_8 = 135.0 / 8.0;
    const auto f_150_77 = 150.0 / 77.0;
    const auto f_165_4 = 165.0 / 4.0;
    const auto f_18_77 = 18.0 / 77.0;
    const auto f_198_7 = 198.0 / 7.0;
    const auto f_200_231 = 200.0 / 231.0;
    const auto f_20_77 = 20.0 / 77.0;
    const auto f_234_77 = 234.0 / 77.0;
    const auto f_24_77 = 24.0 / 77.0;
    const auto f_25_77 = 25.0 / 77.0;
    const auto f_27_1 = 27.0;
    const auto f_297_14 = 297.0 / 14.0;
    const auto f_2_7 = 2.0 / 7.0;
    const auto f_312_77 = 312.0 / 77.0;
    const auto f_33_1 = 33.0;
    const auto f_33_7 = 33.0 / 7.0;
    const auto f_36_7 = 36.0 / 7.0;
    const auto f_36_77 = 36.0 / 77.0;
    const auto f_375_77 = 375.0 / 77.0;
    const auto f_39_14 = 39.0 / 14.0;
    const auto f_39_2 = 39.0 / 2.0;
    const auto f_44_7 = 44.0 / 7.0;
    const auto f_45_1 = 45.0;
    const auto f_468_77 = 468.0 / 77.0;
    const auto f_495_14 = 495.0 / 14.0;
    const auto f_500_77 = 500.0 / 77.0;
    const auto f_50_77 = 50.0 / 77.0;
    const auto f_55_7 = 55.0 / 7.0;
    const auto f_6_11 = 6.0 / 11.0;
    const auto f_6_77 = 6.0 / 77.0;
    const auto f_78_11 = 78.0 / 11.0;
    const auto f_78_7 = 78.0 / 7.0;
    const auto f_78_77 = 78.0 / 77.0;
    const auto f_8_21 = 8.0 / 21.0;
    const auto f_99_4 = 99.0 / 4.0;
    const auto fs_100_77_14 = std::sqrt(20000.0 / 847.0);
    const auto fs_10_231_105 = std::sqrt(500.0 / 2541.0);
    const auto fs_10_231_14 = std::sqrt(200.0 / 7623.0);
    const auto fs_10_231_210 = std::sqrt(1000.0 / 2541.0);
    const auto fs_10_231_231 = std::sqrt(100.0 / 231.0);
    const auto fs_10_231_462 = std::sqrt(200.0 / 231.0);
    const auto fs_10_77_21 = std::sqrt(300.0 / 847.0);
    const auto fs_117_14_3 = std::sqrt(41067.0 / 196.0);
    const auto fs_117_14_7 = std::sqrt(13689.0 / 28.0);
    const auto fs_11_14_30 = std::sqrt(1815.0 / 98.0);
    const auto fs_11_7_2 = std::sqrt(242.0 / 49.0);
    const auto fs_11_7_5 = std::sqrt(605.0 / 49.0);
    const auto fs_125_77_14 = std::sqrt(31250.0 / 847.0);
    const auto fs_12_77_5 = std::sqrt(720.0 / 5929.0);
    const auto fs_150_77_7 = std::sqrt(22500.0 / 847.0);
    const auto fs_156_77_5 = std::sqrt(121680.0 / 5929.0);
    const auto fs_165_8_2 = std::sqrt(27225.0 / 32.0);
    const auto fs_18_77_3 = std::sqrt(972.0 / 5929.0);
    const auto fs_18_77_7 = std::sqrt(324.0 / 847.0);
    const auto fs_1_21_30 = std::sqrt(10.0 / 147.0);
    const auto fs_20_231_21 = std::sqrt(400.0 / 2541.0);
    const auto fs_20_77_7 = std::sqrt(400.0 / 847.0);
    const auto fs_22_7_3 = std::sqrt(1452.0 / 49.0);
    const auto fs_22_7_5 = std::sqrt(2420.0 / 49.0);
    const auto fs_234_77_3 = std::sqrt(164268.0 / 5929.0);
    const auto fs_234_77_7 = std::sqrt(54756.0 / 847.0);
    const auto fs_25_154_14 = std::sqrt(625.0 / 1694.0);
    const auto fs_25_154_210 = std::sqrt(9375.0 / 1694.0);
    const auto fs_25_77_105 = std::sqrt(9375.0 / 847.0);
    const auto fs_25_77_14 = std::sqrt(1250.0 / 847.0);
    const auto fs_25_77_210 = std::sqrt(18750.0 / 847.0);
    const auto fs_25_77_231 = std::sqrt(1875.0 / 77.0);
    const auto fs_25_77_462 = std::sqrt(3750.0 / 77.0);
    const auto fs_2_21_2 = std::sqrt(8.0 / 441.0);
    const auto fs_2_21_5 = std::sqrt(20.0 / 441.0);
    const auto fs_33_2_3 = std::sqrt(3267.0 / 4.0);
    const auto fs_33_2_5 = std::sqrt(5445.0 / 4.0);
    const auto fs_33_4_2 = std::sqrt(1089.0 / 8.0);
    const auto fs_33_4_5 = std::sqrt(5445.0 / 16.0);
    const auto fs_33_8_30 = std::sqrt(16335.0 / 32.0);
    const auto fs_39_14_15 = std::sqrt(22815.0 / 196.0);
    const auto fs_39_14_21 = std::sqrt(4563.0 / 28.0);
    const auto fs_39_14_3 = std::sqrt(4563.0 / 196.0);
    const auto fs_39_14_35 = std::sqrt(7605.0 / 28.0);
    const auto fs_39_14_7 = std::sqrt(1521.0 / 28.0);
    const auto fs_39_7_5 = std::sqrt(7605.0 / 49.0);
    const auto fs_40_231_14 = std::sqrt(3200.0 / 7623.0);
    const auto fs_495_28_2 = std::sqrt(245025.0 / 392.0);
    const auto fs_4_21_3 = std::sqrt(16.0 / 147.0);
    const auto fs_4_21_5 = std::sqrt(80.0 / 441.0);
    const auto fs_50_231_14 = std::sqrt(5000.0 / 7623.0);
    const auto fs_50_77_21 = std::sqrt(7500.0 / 847.0);
    const auto fs_55_14_2 = std::sqrt(3025.0 / 98.0);
    const auto fs_5_21_2 = std::sqrt(50.0 / 441.0);
    const auto fs_5_231_14 = std::sqrt(50.0 / 7623.0);
    const auto fs_5_231_210 = std::sqrt(250.0 / 2541.0);
    const auto fs_6_77_15 = std::sqrt(540.0 / 5929.0);
    const auto fs_6_77_21 = std::sqrt(108.0 / 847.0);
    const auto fs_6_77_3 = std::sqrt(108.0 / 5929.0);
    const auto fs_6_77_35 = std::sqrt(180.0 / 847.0);
    const auto fs_6_77_7 = std::sqrt(36.0 / 847.0);
    const auto fs_75_77_21 = std::sqrt(16875.0 / 847.0);
    const auto fs_78_77_15 = std::sqrt(91260.0 / 5929.0);
    const auto fs_78_77_21 = std::sqrt(18252.0 / 847.0);
    const auto fs_78_77_3 = std::sqrt(18252.0 / 5929.0);
    const auto fs_78_77_35 = std::sqrt(30420.0 / 847.0);
    const auto fs_78_77_7 = std::sqrt(6084.0 / 847.0);
    const auto fs_99_14_2 = std::sqrt(9801.0 / 98.0);
    const auto fs_99_14_5 = std::sqrt(49005.0 / 196.0);
    const auto fs_99_28_30 = std::sqrt(147015.0 / 392.0);
    const auto fs_99_7_3 = std::sqrt(29403.0 / 49.0);
    const auto fs_99_7_5 = std::sqrt(49005.0 / 49.0);

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ab_2, pc_0, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];

        pc_0[k] = e_0 * f_135_8 + e_1 * f_165_4 * h2_0 - e_1 * f_45_1 * r_2 + e_2 * f_117_14 * h4_0 - e_2 * f_495_14 * r_2 * h2_0 + e_2 * f_27_1 * r_4 + e_3 * f_25_77 * h6_0 + e_3 * fs_25_77_462 * h6_p6 - e_3 * f_234_77 * r_2 * h4_0 + e_3 * f_55_7 * r_4 * h2_0 - e_3 * f_36_7 * r_6 - e_4 * f_10_231 * r_2 * h6_0 - e_4 * fs_10_231_462 * r_2 * h6_p6 + e_4 * f_18_77 * r_4 * h4_0 - e_4 * f_10_21 * r_6 * h2_0 + e_4 * f_2_7 * r_8;

        pc_1[k] = - e_1 * fs_165_8_2 * h2_p1 - e_2 * fs_39_14_15 * h4_p1 + e_2 * fs_495_28_2 * r_2 * h2_p1 - e_3 * fs_25_154_14 * h6_p1 + e_3 * fs_25_77_231 * h6_p5 + e_3 * fs_78_77_15 * r_2 * h4_p1 - e_3 * fs_55_14_2 * r_4 * h2_p1 + e_4 * fs_5_231_14 * r_2 * h6_p1 - e_4 * fs_10_231_231 * r_2 * h6_p5 - e_4 * fs_6_77_15 * r_4 * h4_p1 + e_4 * fs_5_21_2 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_p2, ph4_m4, ph4_m3, ph4_m2, ph4_p2, ph4_p4, ph6_m4, ph6_m3, ph6_m2, ph6_p2, ph6_p4, ab_2, pc_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];

        pc_2[k] = e_1 * fs_33_4_5 * h2_p2 + e_2 * fs_117_14_3 * h4_p2 + e_2 * fs_39_14_21 * h4_p4 - e_2 * fs_99_14_5 * r_2 * h2_p2 + e_3 * fs_25_77_14 * h6_p2 + e_3 * fs_25_77_105 * h6_p4 - e_3 * fs_234_77_3 * r_2 * h4_p2 - e_3 * fs_78_77_21 * r_2 * h4_p4 + e_3 * fs_11_7_5 * r_4 * h2_p2 - e_4 * fs_10_231_14 * r_2 * h6_p2 - e_4 * fs_10_231_105 * r_2 * h6_p4 + e_4 * fs_18_77_3 * r_4 * h4_p2 + e_4 * fs_6_77_21 * r_4 * h4_p4 - e_4 * fs_2_21_5 * r_6 * h2_p2;

        pc_3[k] = - e_2 * fs_117_14_7 * h4_m3 - e_3 * fs_50_77_21 * h6_m3 + e_3 * fs_234_77_7 * r_2 * h4_m3 + e_4 * fs_20_231_21 * r_2 * h6_m3 - e_4 * fs_18_77_7 * r_4 * h4_m3;

        pc_4[k] = e_1 * fs_33_4_5 * h2_m2 - e_2 * fs_39_14_21 * h4_m4 + e_2 * fs_117_14_3 * h4_m2 - e_2 * fs_99_14_5 * r_2 * h2_m2 - e_3 * fs_25_77_105 * h6_m4 + e_3 * fs_25_77_14 * h6_m2 + e_3 * fs_78_77_21 * r_2 * h4_m4 - e_3 * fs_234_77_3 * r_2 * h4_m2 + e_3 * fs_11_7_5 * r_4 * h2_m2 + e_4 * fs_10_231_105 * r_2 * h6_m4 - e_4 * fs_10_231_14 * r_2 * h6_m2 - e_4 * fs_6_77_21 * r_4 * h4_m4 + e_4 * fs_18_77_3 * r_4 * h4_m2 - e_4 * fs_2_21_5 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph4_0, ph4_p4, ph6_m6, ph6_m5, ph6_m1, ph6_0, ph6_p4, ab_2, pc_5, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];

        pc_5[k] = - e_1 * fs_165_8_2 * h2_m1 - e_2 * fs_39_14_15 * h4_m1 + e_2 * fs_495_28_2 * r_2 * h2_m1 - e_3 * fs_25_77_231 * h6_m5 - e_3 * fs_25_154_14 * h6_m1 + e_3 * fs_78_77_15 * r_2 * h4_m1 - e_3 * fs_55_14_2 * r_4 * h2_m1 + e_4 * fs_10_231_231 * r_2 * h6_m5 + e_4 * fs_5_231_14 * r_2 * h6_m1 - e_4 * fs_6_77_15 * r_4 * h4_m1 + e_4 * fs_5_21_2 * r_6 * h2_m1;

        pc_6[k] = - e_3 * fs_25_77_462 * h6_m6 + e_4 * fs_10_231_462 * r_2 * h6_m6;

        pc_7[k] = e_0 * f_135_8 - e_1 * f_45_1 * r_2 - e_2 * f_39_2 * h4_0 - e_2 * fs_39_14_35 * h4_p4 + e_2 * f_27_1 * r_4 - e_3 * f_150_77 * h6_0 + e_3 * fs_150_77_7 * h6_p4 + e_3 * f_78_11 * r_2 * h4_0 + e_3 * fs_78_77_35 * r_2 * h4_p4 - e_3 * f_36_7 * r_6 + e_4 * f_20_77 * r_2 * h6_0 - e_4 * fs_20_77_7 * r_2 * h6_p4 - e_4 * f_6_11 * r_4 * h4_0 - e_4 * fs_6_77_35 * r_4 * h4_p4 + e_4 * f_2_7 * r_8;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph2_p1, ph4_m3, ph4_m2, ph4_m1, ph4_p1, ph4_p3, ph6_m3, ph6_m2, ph6_m1, ph6_p1, ph6_p3, ab_2, pc_8, pc_9, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];

        pc_8[k] = - e_1 * fs_33_8_30 * h2_p1 + e_2 * f_78_7 * h4_p1 - e_2 * fs_39_14_7 * h4_p3 + e_2 * fs_99_28_30 * r_2 * h2_p1 + e_3 * fs_25_154_210 * h6_p1 + e_3 * fs_75_77_21 * h6_p3 - e_3 * f_312_77 * r_2 * h4_p1 + e_3 * fs_78_77_7 * r_2 * h4_p3 - e_3 * fs_11_14_30 * r_4 * h2_p1 - e_4 * fs_5_231_210 * r_2 * h6_p1 - e_4 * fs_10_77_21 * r_2 * h6_p3 + e_4 * f_24_77 * r_4 * h4_p1 - e_4 * fs_6_77_7 * r_4 * h4_p3 + e_4 * fs_1_21_30 * r_6 * h2_p1;

        pc_9[k] = e_1 * fs_33_2_5 * h2_m2 - e_2 * fs_39_14_3 * h4_m2 - e_2 * fs_99_7_5 * r_2 * h2_m2 - e_3 * fs_100_77_14 * h6_m2 + e_3 * fs_78_77_3 * r_2 * h4_m2 + e_3 * fs_22_7_5 * r_4 * h2_m2 + e_4 * fs_40_231_14 * r_2 * h6_m2 - e_4 * fs_6_77_3 * r_4 * h4_m2 - e_4 * fs_4_21_5 * r_6 * h2_m2;

        pc_10[k] = - e_1 * fs_33_8_30 * h2_m1 + e_2 * fs_39_14_7 * h4_m3 + e_2 * f_78_7 * h4_m1 + e_2 * fs_99_28_30 * r_2 * h2_m1 - e_3 * fs_75_77_21 * h6_m3 + e_3 * fs_25_154_210 * h6_m1 - e_3 * fs_78_77_7 * r_2 * h4_m3 - e_3 * f_312_77 * r_2 * h4_m1 - e_3 * fs_11_14_30 * r_4 * h2_m1 + e_4 * fs_10_77_21 * r_2 * h6_m3 - e_4 * fs_5_231_210 * r_2 * h6_m1 + e_4 * fs_6_77_7 * r_4 * h4_m3 + e_4 * f_24_77 * r_4 * h4_m1 + e_4 * fs_1_21_30 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph2_0, ph2_p2, ph4_m4, ph4_m1, ph4_0, ph4_p2, ph6_m5, ph6_m4, ph6_m1, ph6_0, ph6_p2, ab_2, pc_11, pc_12, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];

        pc_11[k] = e_2 * fs_39_14_35 * h4_m4 - e_3 * fs_150_77_7 * h6_m4 - e_3 * fs_78_77_35 * r_2 * h4_m4 + e_4 * fs_20_77_7 * r_2 * h6_m4 + e_4 * fs_6_77_35 * r_4 * h4_m4;

        pc_12[k] = e_1 * fs_165_8_2 * h2_m1 + e_2 * fs_39_14_15 * h4_m1 - e_2 * fs_495_28_2 * r_2 * h2_m1 - e_3 * fs_25_77_231 * h6_m5 + e_3 * fs_25_154_14 * h6_m1 - e_3 * fs_78_77_15 * r_2 * h4_m1 + e_3 * fs_55_14_2 * r_4 * h2_m1 + e_4 * fs_10_231_231 * r_2 * h6_m5 - e_4 * fs_5_231_14 * r_2 * h6_m1 + e_4 * fs_6_77_15 * r_4 * h4_m1 - e_4 * fs_5_21_2 * r_6 * h2_m1;

        pc_13[k] = e_0 * f_135_8 - e_1 * f_99_4 * h2_0 + e_1 * fs_33_2_3 * h2_p2 - e_1 * f_45_1 * r_2 + e_2 * f_39_14 * h4_0 - e_2 * fs_39_7_5 * h4_p2 + e_2 * f_297_14 * r_2 * h2_0 - e_2 * fs_99_7_3 * r_2 * h2_p2 + e_2 * f_27_1 * r_4 + e_3 * f_375_77 * h6_0 + e_3 * fs_25_77_210 * h6_p2 - e_3 * f_78_77 * r_2 * h4_0 + e_3 * fs_156_77_5 * r_2 * h4_p2 - e_3 * f_33_7 * r_4 * h2_0 + e_3 * fs_22_7_3 * r_4 * h2_p2 - e_3 * f_36_7 * r_6 - e_4 * f_50_77 * r_2 * h6_0 - e_4 * fs_10_231_210 * r_2 * h6_p2 + e_4 * f_6_77 * r_4 * h4_0 - e_4 * fs_12_77_5 * r_4 * h4_p2 + e_4 * f_2_7 * r_6 * h2_0 - e_4 * fs_4_21_3 * r_6 * h2_p2 + e_4 * f_2_7 * r_8;

        pc_14[k] = - e_1 * fs_33_4_2 * h2_m1 + e_2 * fs_39_14_15 * h4_m1 + e_2 * fs_99_14_2 * r_2 * h2_m1 - e_3 * fs_125_77_14 * h6_m1 - e_3 * fs_78_77_15 * r_2 * h4_m1 - e_3 * fs_11_7_2 * r_4 * h2_m1 + e_4 * fs_50_231_14 * r_2 * h6_m1 + e_4 * fs_6_77_15 * r_4 * h4_m1 + e_4 * fs_2_21_2 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ab_2, pc_15, pc_16, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];

        pc_15[k] = - e_1 * fs_33_2_3 * h2_m2 + e_2 * fs_39_7_5 * h4_m2 + e_2 * fs_99_7_3 * r_2 * h2_m2 - e_3 * fs_25_77_210 * h6_m2 - e_3 * fs_156_77_5 * r_2 * h4_m2 - e_3 * fs_22_7_3 * r_4 * h2_m2 + e_4 * fs_10_231_210 * r_2 * h6_m2 + e_4 * fs_12_77_5 * r_4 * h4_m2 + e_4 * fs_4_21_3 * r_6 * h2_m2;

        pc_16[k] = e_1 * fs_33_8_30 * h2_m1 + e_2 * fs_39_14_7 * h4_m3 - e_2 * f_78_7 * h4_m1 - e_2 * fs_99_28_30 * r_2 * h2_m1 - e_3 * fs_75_77_21 * h6_m3 - e_3 * fs_25_154_210 * h6_m1 - e_3 * fs_78_77_7 * r_2 * h4_m3 + e_3 * f_312_77 * r_2 * h4_m1 + e_3 * fs_11_14_30 * r_4 * h2_m1 + e_4 * fs_10_77_21 * r_2 * h6_m3 + e_4 * fs_5_231_210 * r_2 * h6_m1 + e_4 * fs_6_77_7 * r_4 * h4_m3 - e_4 * f_24_77 * r_4 * h4_m1 - e_4 * fs_1_21_30 * r_6 * h2_m1;

        pc_17[k] = - e_1 * fs_33_4_5 * h2_m2 - e_2 * fs_39_14_21 * h4_m4 - e_2 * fs_117_14_3 * h4_m2 + e_2 * fs_99_14_5 * r_2 * h2_m2 - e_3 * fs_25_77_105 * h6_m4 - e_3 * fs_25_77_14 * h6_m2 + e_3 * fs_78_77_21 * r_2 * h4_m4 + e_3 * fs_234_77_3 * r_2 * h4_m2 - e_3 * fs_11_7_5 * r_4 * h2_m2 + e_4 * fs_10_231_105 * r_2 * h6_m4 + e_4 * fs_10_231_14 * r_2 * h6_m2 - e_4 * fs_6_77_21 * r_4 * h4_m4 - e_4 * fs_18_77_3 * r_4 * h4_m2 + e_4 * fs_2_21_5 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph6_0, ph6_p1, ph6_p2, ph6_p3, ab_2, pc_18, pc_19, pc_20, pc_21, pc_22, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];

        pc_18[k] = e_0 * f_135_8 - e_1 * f_33_1 * h2_0 - e_1 * f_45_1 * r_2 + e_2 * f_117_7 * h4_0 + e_2 * f_198_7 * r_2 * h2_0 + e_2 * f_27_1 * r_4 - e_3 * f_500_77 * h6_0 - e_3 * f_468_77 * r_2 * h4_0 - e_3 * f_44_7 * r_4 * h2_0 - e_3 * f_36_7 * r_6 + e_4 * f_200_231 * r_2 * h6_0 + e_4 * f_36_77 * r_4 * h4_0 + e_4 * f_8_21 * r_6 * h2_0 + e_4 * f_2_7 * r_8;

        pc_19[k] = - e_1 * fs_33_4_2 * h2_p1 + e_2 * fs_39_14_15 * h4_p1 + e_2 * fs_99_14_2 * r_2 * h2_p1 - e_3 * fs_125_77_14 * h6_p1 - e_3 * fs_78_77_15 * r_2 * h4_p1 - e_3 * fs_11_7_2 * r_4 * h2_p1 + e_4 * fs_50_231_14 * r_2 * h6_p1 + e_4 * fs_6_77_15 * r_4 * h4_p1 + e_4 * fs_2_21_2 * r_6 * h2_p1;

        pc_20[k] = e_1 * fs_33_2_5 * h2_p2 - e_2 * fs_39_14_3 * h4_p2 - e_2 * fs_99_7_5 * r_2 * h2_p2 - e_3 * fs_100_77_14 * h6_p2 + e_3 * fs_78_77_3 * r_2 * h4_p2 + e_3 * fs_22_7_5 * r_4 * h2_p2 + e_4 * fs_40_231_14 * r_2 * h6_p2 - e_4 * fs_6_77_3 * r_4 * h4_p2 - e_4 * fs_4_21_5 * r_6 * h2_p2;

        pc_21[k] = - e_2 * fs_117_14_7 * h4_p3 - e_3 * fs_50_77_21 * h6_p3 + e_3 * fs_234_77_7 * r_2 * h4_p3 + e_4 * fs_20_231_21 * r_2 * h6_p3 - e_4 * fs_18_77_7 * r_4 * h4_p3;

        pc_22[k] = e_0 * f_135_8 - e_1 * f_99_4 * h2_0 - e_1 * fs_33_2_3 * h2_p2 - e_1 * f_45_1 * r_2 + e_2 * f_39_14 * h4_0 + e_2 * fs_39_7_5 * h4_p2 + e_2 * f_297_14 * r_2 * h2_0 + e_2 * fs_99_7_3 * r_2 * h2_p2 + e_2 * f_27_1 * r_4 + e_3 * f_375_77 * h6_0 - e_3 * fs_25_77_210 * h6_p2 - e_3 * f_78_77 * r_2 * h4_0 - e_3 * fs_156_77_5 * r_2 * h4_p2 - e_3 * f_33_7 * r_4 * h2_0 - e_3 * fs_22_7_3 * r_4 * h2_p2 - e_3 * f_36_7 * r_6 - e_4 * f_50_77 * r_2 * h6_0 + e_4 * fs_10_231_210 * r_2 * h6_p2 + e_4 * f_6_77 * r_4 * h4_0 + e_4 * fs_12_77_5 * r_4 * h4_p2 + e_4 * f_2_7 * r_6 * h2_0 + e_4 * fs_4_21_3 * r_6 * h2_p2 + e_4 * f_2_7 * r_8;

        pc_23[k] = - e_1 * fs_33_8_30 * h2_p1 + e_2 * f_78_7 * h4_p1 + e_2 * fs_39_14_7 * h4_p3 + e_2 * fs_99_28_30 * r_2 * h2_p1 + e_3 * fs_25_154_210 * h6_p1 - e_3 * fs_75_77_21 * h6_p3 - e_3 * f_312_77 * r_2 * h4_p1 - e_3 * fs_78_77_7 * r_2 * h4_p3 - e_3 * fs_11_14_30 * r_4 * h2_p1 - e_4 * fs_5_231_210 * r_2 * h6_p1 + e_4 * fs_10_77_21 * r_2 * h6_p3 + e_4 * f_24_77 * r_4 * h4_p1 + e_4 * fs_6_77_7 * r_4 * h4_p3 + e_4 * fs_1_21_30 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p4, ph6_0, ph6_p1, ph6_p2, ph6_p4, ph6_p5, ab_2, pc_24, pc_25, pc_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];

        pc_24[k] = e_1 * fs_33_4_5 * h2_p2 + e_2 * fs_117_14_3 * h4_p2 - e_2 * fs_39_14_21 * h4_p4 - e_2 * fs_99_14_5 * r_2 * h2_p2 + e_3 * fs_25_77_14 * h6_p2 - e_3 * fs_25_77_105 * h6_p4 - e_3 * fs_234_77_3 * r_2 * h4_p2 + e_3 * fs_78_77_21 * r_2 * h4_p4 + e_3 * fs_11_7_5 * r_4 * h2_p2 - e_4 * fs_10_231_14 * r_2 * h6_p2 + e_4 * fs_10_231_105 * r_2 * h6_p4 + e_4 * fs_18_77_3 * r_4 * h4_p2 - e_4 * fs_6_77_21 * r_4 * h4_p4 - e_4 * fs_2_21_5 * r_6 * h2_p2;

        pc_25[k] = e_0 * f_135_8 - e_1 * f_45_1 * r_2 - e_2 * f_39_2 * h4_0 + e_2 * fs_39_14_35 * h4_p4 + e_2 * f_27_1 * r_4 - e_3 * f_150_77 * h6_0 - e_3 * fs_150_77_7 * h6_p4 + e_3 * f_78_11 * r_2 * h4_0 - e_3 * fs_78_77_35 * r_2 * h4_p4 - e_3 * f_36_7 * r_6 + e_4 * f_20_77 * r_2 * h6_0 + e_4 * fs_20_77_7 * r_2 * h6_p4 - e_4 * f_6_11 * r_4 * h4_0 + e_4 * fs_6_77_35 * r_4 * h4_p4 + e_4 * f_2_7 * r_8;

        pc_26[k] = - e_1 * fs_165_8_2 * h2_p1 - e_2 * fs_39_14_15 * h4_p1 + e_2 * fs_495_28_2 * r_2 * h2_p1 - e_3 * fs_25_154_14 * h6_p1 - e_3 * fs_25_77_231 * h6_p5 + e_3 * fs_78_77_15 * r_2 * h4_p1 - e_3 * fs_55_14_2 * r_4 * h2_p1 + e_4 * fs_5_231_14 * r_2 * h6_p1 + e_4 * fs_10_231_231 * r_2 * h6_p5 - e_4 * fs_6_77_15 * r_4 * h4_p1 + e_4 * fs_5_21_2 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 9 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph6_p6, ab_2, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];

        pc_27[k] = e_0 * f_135_8 + e_1 * f_165_4 * h2_0 - e_1 * f_45_1 * r_2 + e_2 * f_117_14 * h4_0 - e_2 * f_495_14 * r_2 * h2_0 + e_2 * f_27_1 * r_4 + e_3 * f_25_77 * h6_0 - e_3 * fs_25_77_462 * h6_p6 - e_3 * f_234_77 * r_2 * h4_0 + e_3 * f_55_7 * r_4 * h2_0 - e_3 * f_36_7 * r_6 - e_4 * f_10_231 * r_2 * h6_0 + e_4 * fs_10_231_462 * r_2 * h6_p6 + e_4 * f_18_77 * r_4 * h4_0 - e_4 * f_10_21 * r_6 * h2_0 + e_4 * f_2_7 * r_8;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[49] = {0, 1, 2, 3, 4, 5, 6, 1, 8, 9, 10, 11, 12, 13, 2, 9, 16, 17, 18, 19, 20, 3, 10, 17, 24, 25, 26, 27, 4, 11, 18, 25, 32, 33, 34, 5, 12, 19, 26, 33, 40, 41, 6, 13, 20, 27, 34, 41, 48};

    for (size_t n = 0; n < 49; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
