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



#include "SimdKineticEnergyRecHD.hpp"

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
compute_hd_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 2))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hd_kinetic_energy: Basis functions must be of angular momenta five and two"));
    }

    if (harmonics.size() < 7)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hd_kinetic_energy: Harmonics must reach angular momentum 7"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hd_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 55 * nvalues, 0.0);

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

            const auto ff_0 = fbase * bexp * bexp * bexp * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * bexp * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * bexp * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * bexp * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

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

    const auto f_100_13 = 100.0 / 13.0;
    const auto f_1071_143 = 1071.0 / 143.0;
    const auto f_10_13 = 10.0 / 13.0;
    const auto f_126_143 = 126.0 / 143.0;
    const auto f_130_11 = 130.0 / 11.0;
    const auto f_150_13 = 150.0 / 13.0;
    const auto f_15_1 = 15.0;
    const auto f_195_2 = 195.0 / 2.0;
    const auto f_20_33 = 20.0 / 33.0;
    const auto f_20_39 = 20.0 / 39.0;
    const auto f_25_1 = 25.0;
    const auto f_2_39 = 2.0 / 39.0;
    const auto f_45_2 = 45.0 / 2.0;
    const auto f_4_13 = 4.0 / 13.0;
    const auto f_5_2 = 5.0 / 2.0;
    const auto f_60_13 = 60.0 / 13.0;
    const auto f_65_1 = 65.0;
    const auto f_6_13 = 6.0 / 13.0;
    const auto f_75_2 = 75.0 / 2.0;
    const auto f_90_13 = 90.0 / 13.0;
    const auto fs_102_143_105 = std::sqrt(1092420.0 / 20449.0);
    const auto fs_102_143_5 = std::sqrt(52020.0 / 20449.0);
    const auto fs_102_143_55 = std::sqrt(52020.0 / 1859.0);
    const auto fs_10_13_15 = std::sqrt(1500.0 / 169.0);
    const auto fs_10_13_21 = std::sqrt(2100.0 / 169.0);
    const auto fs_10_13_5 = std::sqrt(500.0 / 169.0);
    const auto fs_10_33_2 = std::sqrt(200.0 / 1089.0);
    const auto fs_10_39_2 = std::sqrt(200.0 / 1521.0);
    const auto fs_10_39_3 = std::sqrt(100.0 / 507.0);
    const auto fs_117_4_10 = std::sqrt(68445.0 / 8.0);
    const auto fs_117_4_7 = std::sqrt(95823.0 / 16.0);
    const auto fs_12_143_105 = std::sqrt(15120.0 / 20449.0);
    const auto fs_12_143_5 = std::sqrt(720.0 / 20449.0);
    const auto fs_12_143_55 = std::sqrt(720.0 / 1859.0);
    const auto fs_13_11_105 = std::sqrt(17745.0 / 121.0);
    const auto fs_13_11_35 = std::sqrt(5915.0 / 121.0);
    const auto fs_13_11_42 = std::sqrt(7098.0 / 121.0);
    const auto fs_13_11_5 = std::sqrt(845.0 / 121.0);
    const auto fs_13_1_14 = std::sqrt(2366.0);
    const auto fs_13_1_3 = std::sqrt(507.0);
    const auto fs_13_1_7 = std::sqrt(1183.0);
    const auto fs_13_22_14 = std::sqrt(1183.0 / 242.0);
    const auto fs_13_22_2 = std::sqrt(169.0 / 242.0);
    const auto fs_13_22_210 = std::sqrt(17745.0 / 242.0);
    const auto fs_13_22_30 = std::sqrt(2535.0 / 242.0);
    const auto fs_13_2_105 = std::sqrt(17745.0 / 4.0);
    const auto fs_13_2_35 = std::sqrt(5915.0 / 4.0);
    const auto fs_13_2_42 = std::sqrt(3549.0 / 2.0);
    const auto fs_13_2_5 = std::sqrt(845.0 / 4.0);
    const auto fs_13_4_14 = std::sqrt(1183.0 / 8.0);
    const auto fs_13_4_2 = std::sqrt(169.0 / 8.0);
    const auto fs_13_4_210 = std::sqrt(17745.0 / 8.0);
    const auto fs_13_4_30 = std::sqrt(2535.0 / 8.0);
    const auto fs_153_143_14 = std::sqrt(327726.0 / 20449.0);
    const auto fs_153_143_30 = std::sqrt(702270.0 / 20449.0);
    const auto fs_153_143_5 = std::sqrt(117045.0 / 20449.0);
    const auto fs_153_286_110 = std::sqrt(117045.0 / 3718.0);
    const auto fs_15_13_30 = std::sqrt(6750.0 / 169.0);
    const auto fs_15_4_30 = std::sqrt(3375.0 / 8.0);
    const auto fs_18_143_14 = std::sqrt(4536.0 / 20449.0);
    const auto fs_18_143_30 = std::sqrt(9720.0 / 20449.0);
    const auto fs_18_143_5 = std::sqrt(1620.0 / 20449.0);
    const auto fs_195_4_2 = std::sqrt(38025.0 / 8.0);
    const auto fs_1_13_30 = std::sqrt(30.0 / 169.0);
    const auto fs_1_33_14 = std::sqrt(14.0 / 1089.0);
    const auto fs_1_33_2 = std::sqrt(2.0 / 1089.0);
    const auto fs_1_33_210 = std::sqrt(70.0 / 363.0);
    const auto fs_1_33_30 = std::sqrt(10.0 / 363.0);
    const auto fs_204_143_15 = std::sqrt(624240.0 / 20449.0);
    const auto fs_204_143_21 = std::sqrt(873936.0 / 20449.0);
    const auto fs_204_143_5 = std::sqrt(208080.0 / 20449.0);
    const auto fs_20_13_14 = std::sqrt(5600.0 / 169.0);
    const auto fs_20_13_35 = std::sqrt(14000.0 / 169.0);
    const auto fs_24_143_15 = std::sqrt(8640.0 / 20449.0);
    const auto fs_24_143_21 = std::sqrt(12096.0 / 20449.0);
    const auto fs_24_143_5 = std::sqrt(2880.0 / 20449.0);
    const auto fs_25_2_2 = std::sqrt(625.0 / 2.0);
    const auto fs_25_2_3 = std::sqrt(1875.0 / 4.0);
    const auto fs_26_11_14 = std::sqrt(9464.0 / 121.0);
    const auto fs_26_11_3 = std::sqrt(2028.0 / 121.0);
    const auto fs_26_11_7 = std::sqrt(4732.0 / 121.0);
    const auto fs_26_1_5 = std::sqrt(3380.0);
    const auto fs_2_11_10 = std::sqrt(40.0 / 121.0);
    const auto fs_2_11_7 = std::sqrt(28.0 / 121.0);
    const auto fs_2_33_105 = std::sqrt(140.0 / 363.0);
    const auto fs_2_33_35 = std::sqrt(140.0 / 1089.0);
    const auto fs_2_33_42 = std::sqrt(56.0 / 363.0);
    const auto fs_2_33_5 = std::sqrt(20.0 / 1089.0);
    const auto fs_2_39_15 = std::sqrt(20.0 / 507.0);
    const auto fs_2_39_21 = std::sqrt(28.0 / 507.0);
    const auto fs_2_39_5 = std::sqrt(20.0 / 1521.0);
    const auto fs_306_143_10 = std::sqrt(936360.0 / 20449.0);
    const auto fs_357_143_5 = std::sqrt(637245.0 / 20449.0);
    const auto fs_35_13_6 = std::sqrt(7350.0 / 169.0);
    const auto fs_35_4_6 = std::sqrt(3675.0 / 8.0);
    const auto fs_36_143_10 = std::sqrt(12960.0 / 20449.0);
    const auto fs_39_11_10 = std::sqrt(15210.0 / 121.0);
    const auto fs_39_11_7 = std::sqrt(10647.0 / 121.0);
    const auto fs_39_1_5 = std::sqrt(7605.0);
    const auto fs_39_2_10 = std::sqrt(7605.0 / 2.0);
    const auto fs_39_2_14 = std::sqrt(10647.0 / 2.0);
    const auto fs_39_2_3 = std::sqrt(4563.0 / 4.0);
    const auto fs_39_2_7 = std::sqrt(10647.0 / 4.0);
    const auto fs_39_4_105 = std::sqrt(159705.0 / 16.0);
    const auto fs_39_4_35 = std::sqrt(53235.0 / 16.0);
    const auto fs_39_4_42 = std::sqrt(31941.0 / 8.0);
    const auto fs_39_4_5 = std::sqrt(7605.0 / 16.0);
    const auto fs_39_8_14 = std::sqrt(10647.0 / 32.0);
    const auto fs_39_8_2 = std::sqrt(1521.0 / 32.0);
    const auto fs_39_8_210 = std::sqrt(159705.0 / 32.0);
    const auto fs_39_8_30 = std::sqrt(22815.0 / 32.0);
    const auto fs_3_143_10 = std::sqrt(90.0 / 20449.0);
    const auto fs_3_143_1430 = std::sqrt(90.0 / 143.0);
    const auto fs_3_143_2 = std::sqrt(18.0 / 20449.0);
    const auto fs_3_143_2002 = std::sqrt(126.0 / 143.0);
    const auto fs_3_143_22 = std::sqrt(18.0 / 1859.0);
    const auto fs_3_143_30 = std::sqrt(270.0 / 20449.0);
    const auto fs_42_143_5 = std::sqrt(8820.0 / 20449.0);
    const auto fs_4_33_14 = std::sqrt(224.0 / 1089.0);
    const auto fs_4_33_3 = std::sqrt(16.0 / 363.0);
    const auto fs_4_33_7 = std::sqrt(112.0 / 1089.0);
    const auto fs_4_39_14 = std::sqrt(224.0 / 1521.0);
    const auto fs_4_39_35 = std::sqrt(560.0 / 1521.0);
    const auto fs_50_13_2 = std::sqrt(5000.0 / 169.0);
    const auto fs_50_13_3 = std::sqrt(7500.0 / 169.0);
    const auto fs_51_143_105 = std::sqrt(273105.0 / 20449.0);
    const auto fs_51_143_143 = std::sqrt(2601.0 / 143.0);
    const auto fs_51_143_165 = std::sqrt(39015.0 / 1859.0);
    const auto fs_51_143_210 = std::sqrt(546210.0 / 20449.0);
    const auto fs_51_143_35 = std::sqrt(91035.0 / 20449.0);
    const auto fs_51_143_66 = std::sqrt(15606.0 / 1859.0);
    const auto fs_51_286_10 = std::sqrt(13005.0 / 40898.0);
    const auto fs_51_286_1430 = std::sqrt(13005.0 / 286.0);
    const auto fs_51_286_2 = std::sqrt(2601.0 / 40898.0);
    const auto fs_51_286_2002 = std::sqrt(18207.0 / 286.0);
    const auto fs_51_286_22 = std::sqrt(2601.0 / 3718.0);
    const auto fs_51_286_30 = std::sqrt(39015.0 / 40898.0);
    const auto fs_52_11_5 = std::sqrt(13520.0 / 121.0);
    const auto fs_5_1_14 = std::sqrt(350.0);
    const auto fs_5_1_35 = std::sqrt(875.0);
    const auto fs_5_2_15 = std::sqrt(375.0 / 4.0);
    const auto fs_5_2_21 = std::sqrt(525.0 / 4.0);
    const auto fs_5_2_5 = std::sqrt(125.0 / 4.0);
    const auto fs_65_11_2 = std::sqrt(8450.0 / 121.0);
    const auto fs_65_2_2 = std::sqrt(4225.0 / 2.0);
    const auto fs_6_143_105 = std::sqrt(3780.0 / 20449.0);
    const auto fs_6_143_143 = std::sqrt(36.0 / 143.0);
    const auto fs_6_143_165 = std::sqrt(540.0 / 1859.0);
    const auto fs_6_143_210 = std::sqrt(7560.0 / 20449.0);
    const auto fs_6_143_35 = std::sqrt(1260.0 / 20449.0);
    const auto fs_6_143_66 = std::sqrt(216.0 / 1859.0);
    const auto fs_7_39_6 = std::sqrt(98.0 / 507.0);
    const auto fs_8_33_5 = std::sqrt(320.0 / 1089.0);
    const auto fs_9_143_110 = std::sqrt(810.0 / 1859.0);

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_m5, ph5_m4, ph5_p3, ph5_p4, ph7_m6, ph7_m5, ph7_m4, ph7_p3, ph7_p4, ph7_p6, ph7_p7, ab_2, pc_0, pc_1, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_0[k] = - e_0 * fs_39_4_105 * h3_p3 - e_1 * fs_5_2_15 * h5_p3 + e_1 * fs_13_2_105 * r_2 * h3_p3 - e_2 * fs_51_286_2 * h7_p3 + e_2 * fs_51_286_2002 * h7_p7 + e_2 * fs_10_13_15 * r_2 * h5_p3 - e_2 * fs_13_11_105 * r_4 * h3_p3 + e_3 * fs_3_143_2 * r_2 * h7_p3 - e_3 * fs_3_143_2002 * r_2 * h7_p7 - e_3 * fs_2_39_15 * r_4 * h5_p3 + e_3 * fs_2_33_105 * r_6 * h3_p3;

        pc_1[k] = e_1 * fs_15_4_30 * h5_p4 + e_2 * fs_51_286_22 * h7_p4 + e_2 * fs_51_143_143 * h7_p6 - e_2 * fs_15_13_30 * r_2 * h5_p4 - e_3 * fs_3_143_22 * r_2 * h7_p4 - e_3 * fs_6_143_143 * r_2 * h7_p6 + e_3 * fs_1_13_30 * r_4 * h5_p4;

        pc_2[k] = - e_1 * f_75_2 * h5_m5 - e_2 * fs_51_143_66 * h7_m5 + e_2 * f_150_13 * r_2 * h5_m5 + e_3 * fs_6_143_66 * r_2 * h7_m5 - e_3 * f_10_13 * r_4 * h5_m5;

        pc_3[k] = e_1 * fs_15_4_30 * h5_m4 - e_2 * fs_51_143_143 * h7_m6 + e_2 * fs_51_286_22 * h7_m4 - e_2 * fs_15_13_30 * r_2 * h5_m4 + e_3 * fs_6_143_143 * r_2 * h7_m6 - e_3 * fs_3_143_22 * r_2 * h7_m4 + e_3 * fs_1_13_30 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_p2, ph3_p3, ph5_m3, ph5_p2, ph5_p3, ph5_p5, ph7_m7, ph7_m3, ph7_p2, ph7_p3, ph7_p5, ph7_p6, ab_2, pc_4, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_4[k] = - e_0 * fs_39_4_105 * h3_m3 - e_1 * fs_5_2_15 * h5_m3 + e_1 * fs_13_2_105 * r_2 * h3_m3 - e_2 * fs_51_286_2002 * h7_m7 - e_2 * fs_51_286_2 * h7_m3 + e_2 * fs_10_13_15 * r_2 * h5_m3 - e_2 * fs_13_11_105 * r_4 * h3_m3 + e_3 * fs_3_143_2002 * r_2 * h7_m7 + e_3 * fs_3_143_2 * r_2 * h7_m3 - e_3 * fs_2_39_15 * r_4 * h5_m3 + e_3 * fs_2_33_105 * r_6 * h3_m3;

        pc_5[k] = - e_0 * fs_117_4_7 * h3_p2 - e_1 * f_15_1 * h5_p2 + e_1 * fs_39_2_7 * r_2 * h3_p2 - e_2 * fs_51_286_10 * h7_p2 + e_2 * fs_51_286_1430 * h7_p6 + e_2 * f_60_13 * r_2 * h5_p2 - e_2 * fs_39_11_7 * r_4 * h3_p2 + e_3 * fs_3_143_10 * r_2 * h7_p2 - e_3 * fs_3_143_1430 * r_2 * h7_p6 - e_3 * f_4_13 * r_4 * h5_p2 + e_3 * fs_2_11_7 * r_6 * h3_p2;

        pc_6[k] = - e_0 * fs_39_4_42 * h3_p3 + e_1 * fs_35_4_6 * h5_p3 - e_1 * fs_15_4_30 * h5_p5 + e_1 * fs_13_2_42 * r_2 * h3_p3 + e_2 * fs_102_143_5 * h7_p3 + e_2 * fs_102_143_55 * h7_p5 - e_2 * fs_35_13_6 * r_2 * h5_p3 + e_2 * fs_15_13_30 * r_2 * h5_p5 - e_2 * fs_13_11_42 * r_4 * h3_p3 - e_3 * fs_12_143_5 * r_2 * h7_p3 - e_3 * fs_12_143_55 * r_2 * h7_p5 + e_3 * fs_7_39_6 * r_4 * h5_p3 - e_3 * fs_1_13_30 * r_4 * h5_p5 + e_3 * fs_2_33_42 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ab_2, pc_7, pc_8, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_7[k] = - e_1 * f_15_1 * h5_m4 - e_2 * fs_51_143_165 * h7_m4 + e_2 * f_60_13 * r_2 * h5_m4 + e_3 * fs_6_143_165 * r_2 * h7_m4 - e_3 * f_4_13 * r_4 * h5_m4;

        pc_8[k] = - e_0 * fs_39_4_42 * h3_m3 + e_1 * fs_15_4_30 * h5_m5 + e_1 * fs_35_4_6 * h5_m3 + e_1 * fs_13_2_42 * r_2 * h3_m3 - e_2 * fs_102_143_55 * h7_m5 + e_2 * fs_102_143_5 * h7_m3 - e_2 * fs_15_13_30 * r_2 * h5_m5 - e_2 * fs_35_13_6 * r_2 * h5_m3 - e_2 * fs_13_11_42 * r_4 * h3_m3 + e_3 * fs_12_143_55 * r_2 * h7_m5 - e_3 * fs_12_143_5 * r_2 * h7_m3 + e_3 * fs_1_13_30 * r_4 * h5_m5 + e_3 * fs_7_39_6 * r_4 * h5_m3 + e_3 * fs_2_33_42 * r_6 * h3_m3;

        pc_9[k] = - e_0 * fs_117_4_7 * h3_m2 - e_1 * f_15_1 * h5_m2 + e_1 * fs_39_2_7 * r_2 * h3_m2 - e_2 * fs_51_286_1430 * h7_m6 - e_2 * fs_51_286_10 * h7_m2 + e_2 * f_60_13 * r_2 * h5_m2 - e_2 * fs_39_11_7 * r_4 * h3_m2 + e_3 * fs_3_143_1430 * r_2 * h7_m6 + e_3 * fs_3_143_10 * r_2 * h7_m2 - e_3 * f_4_13 * r_4 * h5_m2 + e_3 * fs_2_11_7 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_p1, ph3_p2, ph5_m3, ph5_p1, ph5_p2, ph5_p4, ph5_p5, ph7_m3, ph7_p1, ph7_p2, ph7_p4, ph7_p5, ab_2, pc_10, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_10[k] = - e_0 * fs_39_4_35 * h3_p1 - e_1 * fs_5_1_14 * h5_p1 + e_1 * fs_5_2_15 * h5_p5 + e_1 * fs_13_2_35 * r_2 * h3_p1 - e_2 * fs_51_286_30 * h7_p1 + e_2 * fs_153_286_110 * h7_p5 + e_2 * fs_20_13_14 * r_2 * h5_p1 - e_2 * fs_10_13_15 * r_2 * h5_p5 - e_2 * fs_13_11_35 * r_4 * h3_p1 + e_3 * fs_3_143_30 * r_2 * h7_p1 - e_3 * fs_9_143_110 * r_2 * h7_p5 - e_3 * fs_4_39_14 * r_4 * h5_p1 + e_3 * fs_2_39_15 * r_4 * h5_p5 + e_3 * fs_2_33_35 * r_6 * h3_p1;

        pc_11[k] = - e_0 * fs_39_2_14 * h3_p2 + e_1 * fs_25_2_2 * h5_p2 - e_1 * fs_35_4_6 * h5_p4 + e_1 * fs_13_1_14 * r_2 * h3_p2 + e_2 * fs_153_143_5 * h7_p2 + e_2 * fs_153_286_110 * h7_p4 - e_2 * fs_50_13_2 * r_2 * h5_p2 + e_2 * fs_35_13_6 * r_2 * h5_p4 - e_2 * fs_26_11_14 * r_4 * h3_p2 - e_3 * fs_18_143_5 * r_2 * h7_p2 - e_3 * fs_9_143_110 * r_2 * h7_p4 + e_3 * fs_10_39_2 * r_4 * h5_p2 - e_3 * fs_7_39_6 * r_4 * h5_p4 + e_3 * fs_4_33_14 * r_6 * h3_p2;

        pc_12[k] = - e_0 * fs_39_2_7 * h3_m3 + e_1 * f_5_2 * h5_m3 + e_1 * fs_13_1_7 * r_2 * h3_m3 - e_2 * fs_153_143_30 * h7_m3 - e_2 * f_10_13 * r_2 * h5_m3 - e_2 * fs_26_11_7 * r_4 * h3_m3 + e_3 * fs_18_143_30 * r_2 * h7_m3 + e_3 * f_2_39 * r_4 * h5_m3 + e_3 * fs_4_33_7 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_m1, ph5_m5, ph5_m4, ph5_m2, ph5_m1, ph7_m5, ph7_m4, ph7_m2, ph7_m1, ab_2, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_13[k] = - e_0 * fs_39_2_14 * h3_m2 + e_1 * fs_35_4_6 * h5_m4 + e_1 * fs_25_2_2 * h5_m2 + e_1 * fs_13_1_14 * r_2 * h3_m2 - e_2 * fs_153_286_110 * h7_m4 + e_2 * fs_153_143_5 * h7_m2 - e_2 * fs_35_13_6 * r_2 * h5_m4 - e_2 * fs_50_13_2 * r_2 * h5_m2 - e_2 * fs_26_11_14 * r_4 * h3_m2 + e_3 * fs_9_143_110 * r_2 * h7_m4 - e_3 * fs_18_143_5 * r_2 * h7_m2 + e_3 * fs_7_39_6 * r_4 * h5_m4 + e_3 * fs_10_39_2 * r_4 * h5_m2 + e_3 * fs_4_33_14 * r_6 * h3_m2;

        pc_14[k] = - e_0 * fs_39_4_35 * h3_m1 - e_1 * fs_5_2_15 * h5_m5 - e_1 * fs_5_1_14 * h5_m1 + e_1 * fs_13_2_35 * r_2 * h3_m1 - e_2 * fs_153_286_110 * h7_m5 - e_2 * fs_51_286_30 * h7_m1 + e_2 * fs_10_13_15 * r_2 * h5_m5 + e_2 * fs_20_13_14 * r_2 * h5_m1 - e_2 * fs_13_11_35 * r_4 * h3_m1 + e_3 * fs_9_143_110 * r_2 * h7_m5 + e_3 * fs_3_143_30 * r_2 * h7_m1 - e_3 * fs_2_39_15 * r_4 * h5_m5 - e_3 * fs_4_39_14 * r_4 * h5_m1 + e_3 * fs_2_33_35 * r_6 * h3_m1;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph3_p3, ph5_0, ph5_p1, ph5_p3, ph5_p4, ph7_0, ph7_p1, ph7_p3, ph7_p4, ab_2, pc_15, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_15[k] = - e_0 * fs_39_4_35 * h3_0 - e_1 * fs_5_1_35 * h5_0 + e_1 * f_15_1 * h5_p4 + e_1 * fs_13_2_35 * r_2 * h3_0 - e_2 * fs_51_143_35 * h7_0 + e_2 * fs_51_143_165 * h7_p4 + e_2 * fs_20_13_35 * r_2 * h5_0 - e_2 * f_60_13 * r_2 * h5_p4 - e_2 * fs_13_11_35 * r_4 * h3_0 + e_3 * fs_6_143_35 * r_2 * h7_0 - e_3 * fs_6_143_165 * r_2 * h7_p4 - e_3 * fs_4_39_35 * r_4 * h5_0 + e_3 * f_4_13 * r_4 * h5_p4 + e_3 * fs_2_33_35 * r_6 * h3_0;

        pc_16[k] = - e_0 * fs_39_8_210 * h3_p1 - e_0 * fs_39_8_14 * h3_p3 + e_1 * fs_5_2_21 * h5_p1 - e_1 * fs_25_2_2 * h5_p3 + e_1 * fs_13_4_210 * r_2 * h3_p1 + e_1 * fs_13_4_14 * r_2 * h3_p3 + e_2 * fs_204_143_5 * h7_p1 + e_2 * fs_204_143_15 * h7_p3 - e_2 * fs_10_13_21 * r_2 * h5_p1 + e_2 * fs_50_13_2 * r_2 * h5_p3 - e_2 * fs_13_22_210 * r_4 * h3_p1 - e_2 * fs_13_22_14 * r_4 * h3_p3 - e_3 * fs_24_143_5 * r_2 * h7_p1 - e_3 * fs_24_143_15 * r_2 * h7_p3 + e_3 * fs_2_39_21 * r_4 * h5_p1 - e_3 * fs_10_39_2 * r_4 * h5_p3 + e_3 * fs_1_33_210 * r_6 * h3_p1 + e_3 * fs_1_33_14 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph3_m1, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ab_2, pc_17, pc_18, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_17[k] = - e_0 * fs_117_4_7 * h3_m2 + e_1 * f_15_1 * h5_m2 + e_1 * fs_39_2_7 * r_2 * h3_m2 - e_2 * fs_306_143_10 * h7_m2 - e_2 * f_60_13 * r_2 * h5_m2 - e_2 * fs_39_11_7 * r_4 * h3_m2 + e_3 * fs_36_143_10 * r_2 * h7_m2 + e_3 * f_4_13 * r_4 * h5_m2 + e_3 * fs_2_11_7 * r_6 * h3_m2;

        pc_18[k] = e_0 * fs_39_8_14 * h3_m3 - e_0 * fs_39_8_210 * h3_m1 + e_1 * fs_25_2_2 * h5_m3 + e_1 * fs_5_2_21 * h5_m1 - e_1 * fs_13_4_14 * r_2 * h3_m3 + e_1 * fs_13_4_210 * r_2 * h3_m1 - e_2 * fs_204_143_15 * h7_m3 + e_2 * fs_204_143_5 * h7_m1 - e_2 * fs_50_13_2 * r_2 * h5_m3 - e_2 * fs_10_13_21 * r_2 * h5_m1 + e_2 * fs_13_22_14 * r_4 * h3_m3 - e_2 * fs_13_22_210 * r_4 * h3_m1 + e_3 * fs_24_143_15 * r_2 * h7_m3 - e_3 * fs_24_143_5 * r_2 * h7_m1 + e_3 * fs_10_39_2 * r_4 * h5_m3 + e_3 * fs_2_39_21 * r_4 * h5_m1 - e_3 * fs_1_33_14 * r_6 * h3_m3 + e_3 * fs_1_33_210 * r_6 * h3_m1;

        pc_19[k] = - e_1 * f_15_1 * h5_m4 - e_2 * fs_51_143_165 * h7_m4 + e_2 * f_60_13 * r_2 * h5_m4 + e_3 * fs_6_143_165 * r_2 * h7_m4 - e_3 * f_4_13 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph7_0, ph7_p1, ph7_p2, ph7_p3, ab_2, pc_20, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_20[k] = e_0 * fs_39_8_30 * h3_p1 + e_0 * fs_39_8_2 * h3_p3 + e_1 * fs_25_2_3 * h5_p1 + e_1 * fs_5_1_14 * h5_p3 - e_1 * fs_13_4_30 * r_2 * h3_p1 - e_1 * fs_13_4_2 * r_2 * h3_p3 + e_2 * fs_51_143_35 * h7_p1 + e_2 * fs_51_143_105 * h7_p3 - e_2 * fs_50_13_3 * r_2 * h5_p1 - e_2 * fs_20_13_14 * r_2 * h5_p3 + e_2 * fs_13_22_30 * r_4 * h3_p1 + e_2 * fs_13_22_2 * r_4 * h3_p3 - e_3 * fs_6_143_35 * r_2 * h7_p1 - e_3 * fs_6_143_105 * r_2 * h7_p3 + e_3 * fs_10_39_3 * r_4 * h5_p1 + e_3 * fs_4_39_14 * r_4 * h5_p3 - e_3 * fs_1_33_30 * r_6 * h3_p1 - e_3 * fs_1_33_2 * r_6 * h3_p3;

        pc_21[k] = - e_0 * fs_39_1_5 * h3_0 - e_0 * fs_39_2_3 * h3_p2 + e_1 * fs_5_2_5 * h5_0 - e_1 * fs_5_2_21 * h5_p2 + e_1 * fs_26_1_5 * r_2 * h3_0 + e_1 * fs_13_1_3 * r_2 * h3_p2 + e_2 * fs_357_143_5 * h7_0 + e_2 * fs_51_143_210 * h7_p2 - e_2 * fs_10_13_5 * r_2 * h5_0 + e_2 * fs_10_13_21 * r_2 * h5_p2 - e_2 * fs_52_11_5 * r_4 * h3_0 - e_2 * fs_26_11_3 * r_4 * h3_p2 - e_3 * fs_42_143_5 * r_2 * h7_0 - e_3 * fs_6_143_210 * r_2 * h7_p2 + e_3 * fs_2_39_5 * r_4 * h5_0 - e_3 * fs_2_39_21 * r_4 * h5_p2 + e_3 * fs_8_33_5 * r_6 * h3_0 + e_3 * fs_4_33_3 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph3_m1, ph3_0, ph5_m3, ph5_m2, ph5_m1, ph5_0, ph7_m3, ph7_m2, ph7_m1, ph7_0, ab_2, pc_22, pc_23, pc_24, pc_25, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];

        pc_22[k] = - e_0 * fs_117_4_10 * h3_m1 + e_1 * f_45_2 * h5_m1 + e_1 * fs_39_2_10 * r_2 * h3_m1 - e_2 * fs_102_143_105 * h7_m1 - e_2 * f_90_13 * r_2 * h5_m1 - e_2 * fs_39_11_10 * r_4 * h3_m1 + e_3 * fs_12_143_105 * r_2 * h7_m1 + e_3 * f_6_13 * r_4 * h5_m1 + e_3 * fs_2_11_10 * r_6 * h3_m1;

        pc_23[k] = e_0 * fs_39_2_3 * h3_m2 + e_1 * fs_5_2_21 * h5_m2 - e_1 * fs_13_1_3 * r_2 * h3_m2 - e_2 * fs_51_143_210 * h7_m2 - e_2 * fs_10_13_21 * r_2 * h5_m2 + e_2 * fs_26_11_3 * r_4 * h3_m2 + e_3 * fs_6_143_210 * r_2 * h7_m2 + e_3 * fs_2_39_21 * r_4 * h5_m2 - e_3 * fs_4_33_3 * r_6 * h3_m2;

        pc_24[k] = - e_0 * fs_39_8_2 * h3_m3 - e_0 * fs_39_8_30 * h3_m1 - e_1 * fs_5_1_14 * h5_m3 - e_1 * fs_25_2_3 * h5_m1 + e_1 * fs_13_4_2 * r_2 * h3_m3 + e_1 * fs_13_4_30 * r_2 * h3_m1 - e_2 * fs_51_143_105 * h7_m3 - e_2 * fs_51_143_35 * h7_m1 + e_2 * fs_20_13_14 * r_2 * h5_m3 + e_2 * fs_50_13_3 * r_2 * h5_m1 - e_2 * fs_13_22_2 * r_4 * h3_m3 - e_2 * fs_13_22_30 * r_4 * h3_m1 + e_3 * fs_6_143_105 * r_2 * h7_m3 + e_3 * fs_6_143_35 * r_2 * h7_m1 - e_3 * fs_4_39_14 * r_4 * h5_m3 - e_3 * fs_10_39_3 * r_4 * h5_m1 + e_3 * fs_1_33_2 * r_6 * h3_m3 + e_3 * fs_1_33_30 * r_6 * h3_m1;

        pc_25[k] = - e_0 * fs_39_4_5 * h3_m2 - e_1 * fs_5_1_35 * h5_m2 + e_1 * fs_13_2_5 * r_2 * h3_m2 - e_2 * fs_153_143_14 * h7_m2 + e_2 * fs_20_13_35 * r_2 * h5_m2 - e_2 * fs_13_11_5 * r_4 * h3_m2 + e_3 * fs_18_143_14 * r_2 * h7_m2 - e_3 * fs_4_39_35 * r_4 * h5_m2 + e_3 * fs_2_33_5 * r_6 * h3_m2;

        pc_26[k] = e_0 * fs_195_4_2 * h3_m1 + e_1 * fs_5_2_5 * h5_m1 - e_1 * fs_65_2_2 * r_2 * h3_m1 - e_2 * fs_204_143_21 * h7_m1 - e_2 * fs_10_13_5 * r_2 * h5_m1 + e_2 * fs_65_11_2 * r_4 * h3_m1 + e_3 * fs_24_143_21 * r_2 * h7_m1 + e_3 * fs_2_39_5 * r_4 * h5_m1 - e_3 * fs_10_33_2 * r_6 * h3_m1;

        pc_27[k] = - e_0 * f_195_2 * h3_0 + e_1 * f_25_1 * h5_0 + e_1 * f_65_1 * r_2 * h3_0 - e_2 * f_1071_143 * h7_0 - e_2 * f_100_13 * r_2 * h5_0 - e_2 * f_130_11 * r_4 * h3_0 + e_3 * f_126_143 * r_2 * h7_0 + e_3 * f_20_39 * r_4 * h5_0 + e_3 * f_20_33 * r_6 * h3_0;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m1, ph3_p1, ph3_p2, ph5_m3, ph5_m1, ph5_p1, ph5_p2, ph7_m3, ph7_m1, ph7_p1, ph7_p2, ab_2, pc_28, pc_29, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_28[k] = e_0 * fs_195_4_2 * h3_p1 + e_1 * fs_5_2_5 * h5_p1 - e_1 * fs_65_2_2 * r_2 * h3_p1 - e_2 * fs_204_143_21 * h7_p1 - e_2 * fs_10_13_5 * r_2 * h5_p1 + e_2 * fs_65_11_2 * r_4 * h3_p1 + e_3 * fs_24_143_21 * r_2 * h7_p1 + e_3 * fs_2_39_5 * r_4 * h5_p1 - e_3 * fs_10_33_2 * r_6 * h3_p1;

        pc_29[k] = - e_0 * fs_39_4_5 * h3_p2 - e_1 * fs_5_1_35 * h5_p2 + e_1 * fs_13_2_5 * r_2 * h3_p2 - e_2 * fs_153_143_14 * h7_p2 + e_2 * fs_20_13_35 * r_2 * h5_p2 - e_2 * fs_13_11_5 * r_4 * h3_p2 + e_3 * fs_18_143_14 * r_2 * h7_p2 - e_3 * fs_4_39_35 * r_4 * h5_p2 + e_3 * fs_2_33_5 * r_6 * h3_p2;

        pc_30[k] = - e_0 * fs_39_8_2 * h3_m3 + e_0 * fs_39_8_30 * h3_m1 - e_1 * fs_5_1_14 * h5_m3 + e_1 * fs_25_2_3 * h5_m1 + e_1 * fs_13_4_2 * r_2 * h3_m3 - e_1 * fs_13_4_30 * r_2 * h3_m1 - e_2 * fs_51_143_105 * h7_m3 + e_2 * fs_51_143_35 * h7_m1 + e_2 * fs_20_13_14 * r_2 * h5_m3 - e_2 * fs_50_13_3 * r_2 * h5_m1 - e_2 * fs_13_22_2 * r_4 * h3_m3 + e_2 * fs_13_22_30 * r_4 * h3_m1 + e_3 * fs_6_143_105 * r_2 * h7_m3 - e_3 * fs_6_143_35 * r_2 * h7_m1 - e_3 * fs_4_39_14 * r_4 * h5_m3 + e_3 * fs_10_39_3 * r_4 * h5_m1 + e_3 * fs_1_33_2 * r_6 * h3_m3 - e_3 * fs_1_33_30 * r_6 * h3_m1;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_0, ph3_p1, ph3_p2, ph5_m2, ph5_0, ph5_p1, ph5_p2, ph7_m2, ph7_0, ph7_p1, ph7_p2, ab_2, pc_31, pc_32, pc_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_31[k] = e_0 * fs_39_2_3 * h3_m2 + e_1 * fs_5_2_21 * h5_m2 - e_1 * fs_13_1_3 * r_2 * h3_m2 - e_2 * fs_51_143_210 * h7_m2 - e_2 * fs_10_13_21 * r_2 * h5_m2 + e_2 * fs_26_11_3 * r_4 * h3_m2 + e_3 * fs_6_143_210 * r_2 * h7_m2 + e_3 * fs_2_39_21 * r_4 * h5_m2 - e_3 * fs_4_33_3 * r_6 * h3_m2;

        pc_32[k] = - e_0 * fs_117_4_10 * h3_p1 + e_1 * f_45_2 * h5_p1 + e_1 * fs_39_2_10 * r_2 * h3_p1 - e_2 * fs_102_143_105 * h7_p1 - e_2 * f_90_13 * r_2 * h5_p1 - e_2 * fs_39_11_10 * r_4 * h3_p1 + e_3 * fs_12_143_105 * r_2 * h7_p1 + e_3 * f_6_13 * r_4 * h5_p1 + e_3 * fs_2_11_10 * r_6 * h3_p1;

        pc_33[k] = - e_0 * fs_39_1_5 * h3_0 + e_0 * fs_39_2_3 * h3_p2 + e_1 * fs_5_2_5 * h5_0 + e_1 * fs_5_2_21 * h5_p2 + e_1 * fs_26_1_5 * r_2 * h3_0 - e_1 * fs_13_1_3 * r_2 * h3_p2 + e_2 * fs_357_143_5 * h7_0 - e_2 * fs_51_143_210 * h7_p2 - e_2 * fs_10_13_5 * r_2 * h5_0 - e_2 * fs_10_13_21 * r_2 * h5_p2 - e_2 * fs_52_11_5 * r_4 * h3_0 + e_2 * fs_26_11_3 * r_4 * h3_p2 - e_3 * fs_42_143_5 * r_2 * h7_0 + e_3 * fs_6_143_210 * r_2 * h7_p2 + e_3 * fs_2_39_5 * r_4 * h5_0 + e_3 * fs_2_39_21 * r_4 * h5_p2 + e_3 * fs_8_33_5 * r_6 * h3_0 - e_3 * fs_4_33_3 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p3, ph5_m4, ph5_p1, ph5_p3, ph7_m4, ph7_p1, ph7_p3, ab_2, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];

        pc_34[k] = e_0 * fs_39_8_30 * h3_p1 - e_0 * fs_39_8_2 * h3_p3 + e_1 * fs_25_2_3 * h5_p1 - e_1 * fs_5_1_14 * h5_p3 - e_1 * fs_13_4_30 * r_2 * h3_p1 + e_1 * fs_13_4_2 * r_2 * h3_p3 + e_2 * fs_51_143_35 * h7_p1 - e_2 * fs_51_143_105 * h7_p3 - e_2 * fs_50_13_3 * r_2 * h5_p1 + e_2 * fs_20_13_14 * r_2 * h5_p3 + e_2 * fs_13_22_30 * r_4 * h3_p1 - e_2 * fs_13_22_2 * r_4 * h3_p3 - e_3 * fs_6_143_35 * r_2 * h7_p1 + e_3 * fs_6_143_105 * r_2 * h7_p3 + e_3 * fs_10_39_3 * r_4 * h5_p1 - e_3 * fs_4_39_14 * r_4 * h5_p3 - e_3 * fs_1_33_30 * r_6 * h3_p1 + e_3 * fs_1_33_2 * r_6 * h3_p3;

        pc_35[k] = - e_1 * f_15_1 * h5_m4 - e_2 * fs_51_143_165 * h7_m4 + e_2 * f_60_13 * r_2 * h5_m4 + e_3 * fs_6_143_165 * r_2 * h7_m4 - e_3 * f_4_13 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m1, ph3_p2, ph5_m3, ph5_m1, ph5_p2, ph7_m3, ph7_m1, ph7_p2, ab_2, pc_36, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_36[k] = e_0 * fs_39_8_14 * h3_m3 + e_0 * fs_39_8_210 * h3_m1 + e_1 * fs_25_2_2 * h5_m3 - e_1 * fs_5_2_21 * h5_m1 - e_1 * fs_13_4_14 * r_2 * h3_m3 - e_1 * fs_13_4_210 * r_2 * h3_m1 - e_2 * fs_204_143_15 * h7_m3 - e_2 * fs_204_143_5 * h7_m1 - e_2 * fs_50_13_2 * r_2 * h5_m3 + e_2 * fs_10_13_21 * r_2 * h5_m1 + e_2 * fs_13_22_14 * r_4 * h3_m3 + e_2 * fs_13_22_210 * r_4 * h3_m1 + e_3 * fs_24_143_15 * r_2 * h7_m3 + e_3 * fs_24_143_5 * r_2 * h7_m1 + e_3 * fs_10_39_2 * r_4 * h5_m3 - e_3 * fs_2_39_21 * r_4 * h5_m1 - e_3 * fs_1_33_14 * r_6 * h3_m3 - e_3 * fs_1_33_210 * r_6 * h3_m1;

        pc_37[k] = - e_0 * fs_117_4_7 * h3_p2 + e_1 * f_15_1 * h5_p2 + e_1 * fs_39_2_7 * r_2 * h3_p2 - e_2 * fs_306_143_10 * h7_p2 - e_2 * f_60_13 * r_2 * h5_p2 - e_2 * fs_39_11_7 * r_4 * h3_p2 + e_3 * fs_36_143_10 * r_2 * h7_p2 + e_3 * f_4_13 * r_4 * h5_p2 + e_3 * fs_2_11_7 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph3_p3, ph5_0, ph5_p1, ph5_p3, ph5_p4, ph7_0, ph7_p1, ph7_p3, ph7_p4, ab_2, pc_38, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_38[k] = - e_0 * fs_39_8_210 * h3_p1 + e_0 * fs_39_8_14 * h3_p3 + e_1 * fs_5_2_21 * h5_p1 + e_1 * fs_25_2_2 * h5_p3 + e_1 * fs_13_4_210 * r_2 * h3_p1 - e_1 * fs_13_4_14 * r_2 * h3_p3 + e_2 * fs_204_143_5 * h7_p1 - e_2 * fs_204_143_15 * h7_p3 - e_2 * fs_10_13_21 * r_2 * h5_p1 - e_2 * fs_50_13_2 * r_2 * h5_p3 - e_2 * fs_13_22_210 * r_4 * h3_p1 + e_2 * fs_13_22_14 * r_4 * h3_p3 - e_3 * fs_24_143_5 * r_2 * h7_p1 + e_3 * fs_24_143_15 * r_2 * h7_p3 + e_3 * fs_2_39_21 * r_4 * h5_p1 + e_3 * fs_10_39_2 * r_4 * h5_p3 + e_3 * fs_1_33_210 * r_6 * h3_p1 - e_3 * fs_1_33_14 * r_6 * h3_p3;

        pc_39[k] = - e_0 * fs_39_4_35 * h3_0 - e_1 * fs_5_1_35 * h5_0 - e_1 * f_15_1 * h5_p4 + e_1 * fs_13_2_35 * r_2 * h3_0 - e_2 * fs_51_143_35 * h7_0 - e_2 * fs_51_143_165 * h7_p4 + e_2 * fs_20_13_35 * r_2 * h5_0 + e_2 * f_60_13 * r_2 * h5_p4 - e_2 * fs_13_11_35 * r_4 * h3_0 + e_3 * fs_6_143_35 * r_2 * h7_0 + e_3 * fs_6_143_165 * r_2 * h7_p4 - e_3 * fs_4_39_35 * r_4 * h5_0 - e_3 * f_4_13 * r_4 * h5_p4 + e_3 * fs_2_33_35 * r_6 * h3_0;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_m1, ph3_p3, ph5_m5, ph5_m4, ph5_m2, ph5_m1, ph5_p3, ph7_m5, ph7_m4, ph7_m2, ph7_m1, ph7_p3, ab_2, pc_40, pc_41, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p3 = ph7_p3[k];

        pc_40[k] = e_0 * fs_39_4_35 * h3_m1 - e_1 * fs_5_2_15 * h5_m5 + e_1 * fs_5_1_14 * h5_m1 - e_1 * fs_13_2_35 * r_2 * h3_m1 - e_2 * fs_153_286_110 * h7_m5 + e_2 * fs_51_286_30 * h7_m1 + e_2 * fs_10_13_15 * r_2 * h5_m5 - e_2 * fs_20_13_14 * r_2 * h5_m1 + e_2 * fs_13_11_35 * r_4 * h3_m1 + e_3 * fs_9_143_110 * r_2 * h7_m5 - e_3 * fs_3_143_30 * r_2 * h7_m1 - e_3 * fs_2_39_15 * r_4 * h5_m5 + e_3 * fs_4_39_14 * r_4 * h5_m1 - e_3 * fs_2_33_35 * r_6 * h3_m1;

        pc_41[k] = e_0 * fs_39_2_14 * h3_m2 + e_1 * fs_35_4_6 * h5_m4 - e_1 * fs_25_2_2 * h5_m2 - e_1 * fs_13_1_14 * r_2 * h3_m2 - e_2 * fs_153_286_110 * h7_m4 - e_2 * fs_153_143_5 * h7_m2 - e_2 * fs_35_13_6 * r_2 * h5_m4 + e_2 * fs_50_13_2 * r_2 * h5_m2 + e_2 * fs_26_11_14 * r_4 * h3_m2 + e_3 * fs_9_143_110 * r_2 * h7_m4 + e_3 * fs_18_143_5 * r_2 * h7_m2 + e_3 * fs_7_39_6 * r_4 * h5_m4 - e_3 * fs_10_39_2 * r_4 * h5_m2 - e_3 * fs_4_33_14 * r_6 * h3_m2;

        pc_42[k] = - e_0 * fs_39_2_7 * h3_p3 + e_1 * f_5_2 * h5_p3 + e_1 * fs_13_1_7 * r_2 * h3_p3 - e_2 * fs_153_143_30 * h7_p3 - e_2 * f_10_13 * r_2 * h5_p3 - e_2 * fs_26_11_7 * r_4 * h3_p3 + e_3 * fs_18_143_30 * r_2 * h7_p3 + e_3 * f_2_39 * r_4 * h5_p3 + e_3 * fs_4_33_7 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph5_p4, ph5_p5, ph7_p1, ph7_p2, ph7_p4, ph7_p5, ab_2, pc_43, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_43[k] = - e_0 * fs_39_2_14 * h3_p2 + e_1 * fs_25_2_2 * h5_p2 + e_1 * fs_35_4_6 * h5_p4 + e_1 * fs_13_1_14 * r_2 * h3_p2 + e_2 * fs_153_143_5 * h7_p2 - e_2 * fs_153_286_110 * h7_p4 - e_2 * fs_50_13_2 * r_2 * h5_p2 - e_2 * fs_35_13_6 * r_2 * h5_p4 - e_2 * fs_26_11_14 * r_4 * h3_p2 - e_3 * fs_18_143_5 * r_2 * h7_p2 + e_3 * fs_9_143_110 * r_2 * h7_p4 + e_3 * fs_10_39_2 * r_4 * h5_p2 + e_3 * fs_7_39_6 * r_4 * h5_p4 + e_3 * fs_4_33_14 * r_6 * h3_p2;

        pc_44[k] = - e_0 * fs_39_4_35 * h3_p1 - e_1 * fs_5_1_14 * h5_p1 - e_1 * fs_5_2_15 * h5_p5 + e_1 * fs_13_2_35 * r_2 * h3_p1 - e_2 * fs_51_286_30 * h7_p1 - e_2 * fs_153_286_110 * h7_p5 + e_2 * fs_20_13_14 * r_2 * h5_p1 + e_2 * fs_10_13_15 * r_2 * h5_p5 - e_2 * fs_13_11_35 * r_4 * h3_p1 + e_3 * fs_3_143_30 * r_2 * h7_p1 + e_3 * fs_9_143_110 * r_2 * h7_p5 - e_3 * fs_4_39_14 * r_4 * h5_p1 - e_3 * fs_2_39_15 * r_4 * h5_p5 + e_3 * fs_2_33_35 * r_6 * h3_p1;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m5, ph5_m3, ph5_m2, ph5_p4, ph7_m6, ph7_m5, ph7_m3, ph7_m2, ph7_p4, ab_2, pc_45, pc_46, pc_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p4 = ph7_p4[k];

        pc_45[k] = e_0 * fs_117_4_7 * h3_m2 + e_1 * f_15_1 * h5_m2 - e_1 * fs_39_2_7 * r_2 * h3_m2 - e_2 * fs_51_286_1430 * h7_m6 + e_2 * fs_51_286_10 * h7_m2 - e_2 * f_60_13 * r_2 * h5_m2 + e_2 * fs_39_11_7 * r_4 * h3_m2 + e_3 * fs_3_143_1430 * r_2 * h7_m6 - e_3 * fs_3_143_10 * r_2 * h7_m2 + e_3 * f_4_13 * r_4 * h5_m2 - e_3 * fs_2_11_7 * r_6 * h3_m2;

        pc_46[k] = e_0 * fs_39_4_42 * h3_m3 + e_1 * fs_15_4_30 * h5_m5 - e_1 * fs_35_4_6 * h5_m3 - e_1 * fs_13_2_42 * r_2 * h3_m3 - e_2 * fs_102_143_55 * h7_m5 - e_2 * fs_102_143_5 * h7_m3 - e_2 * fs_15_13_30 * r_2 * h5_m5 + e_2 * fs_35_13_6 * r_2 * h5_m3 + e_2 * fs_13_11_42 * r_4 * h3_m3 + e_3 * fs_12_143_55 * r_2 * h7_m5 + e_3 * fs_12_143_5 * r_2 * h7_m3 + e_3 * fs_1_13_30 * r_4 * h5_m5 - e_3 * fs_7_39_6 * r_4 * h5_m3 - e_3 * fs_2_33_42 * r_6 * h3_m3;

        pc_47[k] = - e_1 * f_15_1 * h5_p4 - e_2 * fs_51_143_165 * h7_p4 + e_2 * f_60_13 * r_2 * h5_p4 + e_3 * fs_6_143_165 * r_2 * h7_p4 - e_3 * f_4_13 * r_4 * h5_p4;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_p2, ph3_p3, ph5_m3, ph5_p2, ph5_p3, ph5_p5, ph7_m7, ph7_m3, ph7_p2, ph7_p3, ph7_p5, ph7_p6, ab_2, pc_48, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_48[k] = - e_0 * fs_39_4_42 * h3_p3 + e_1 * fs_35_4_6 * h5_p3 + e_1 * fs_15_4_30 * h5_p5 + e_1 * fs_13_2_42 * r_2 * h3_p3 + e_2 * fs_102_143_5 * h7_p3 - e_2 * fs_102_143_55 * h7_p5 - e_2 * fs_35_13_6 * r_2 * h5_p3 - e_2 * fs_15_13_30 * r_2 * h5_p5 - e_2 * fs_13_11_42 * r_4 * h3_p3 - e_3 * fs_12_143_5 * r_2 * h7_p3 + e_3 * fs_12_143_55 * r_2 * h7_p5 + e_3 * fs_7_39_6 * r_4 * h5_p3 + e_3 * fs_1_13_30 * r_4 * h5_p5 + e_3 * fs_2_33_42 * r_6 * h3_p3;

        pc_49[k] = - e_0 * fs_117_4_7 * h3_p2 - e_1 * f_15_1 * h5_p2 + e_1 * fs_39_2_7 * r_2 * h3_p2 - e_2 * fs_51_286_10 * h7_p2 - e_2 * fs_51_286_1430 * h7_p6 + e_2 * f_60_13 * r_2 * h5_p2 - e_2 * fs_39_11_7 * r_4 * h3_p2 + e_3 * fs_3_143_10 * r_2 * h7_p2 + e_3 * fs_3_143_1430 * r_2 * h7_p6 - e_3 * f_4_13 * r_4 * h5_p2 + e_3 * fs_2_11_7 * r_6 * h3_p2;

        pc_50[k] = e_0 * fs_39_4_105 * h3_m3 + e_1 * fs_5_2_15 * h5_m3 - e_1 * fs_13_2_105 * r_2 * h3_m3 - e_2 * fs_51_286_2002 * h7_m7 + e_2 * fs_51_286_2 * h7_m3 - e_2 * fs_10_13_15 * r_2 * h5_m3 + e_2 * fs_13_11_105 * r_4 * h3_m3 + e_3 * fs_3_143_2002 * r_2 * h7_m7 - e_3 * fs_3_143_2 * r_2 * h7_m3 + e_3 * fs_2_39_15 * r_4 * h5_m3 - e_3 * fs_2_33_105 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 19 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_m4, ph5_p3, ph5_p4, ph5_p5, ph7_m6, ph7_m4, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ph7_p7, ab_2, pc_51, pc_52, pc_53, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_51[k] = - e_1 * fs_15_4_30 * h5_m4 - e_2 * fs_51_143_143 * h7_m6 - e_2 * fs_51_286_22 * h7_m4 + e_2 * fs_15_13_30 * r_2 * h5_m4 + e_3 * fs_6_143_143 * r_2 * h7_m6 + e_3 * fs_3_143_22 * r_2 * h7_m4 - e_3 * fs_1_13_30 * r_4 * h5_m4;

        pc_52[k] = - e_1 * f_75_2 * h5_p5 - e_2 * fs_51_143_66 * h7_p5 + e_2 * f_150_13 * r_2 * h5_p5 + e_3 * fs_6_143_66 * r_2 * h7_p5 - e_3 * f_10_13 * r_4 * h5_p5;

        pc_53[k] = e_1 * fs_15_4_30 * h5_p4 + e_2 * fs_51_286_22 * h7_p4 - e_2 * fs_51_143_143 * h7_p6 - e_2 * fs_15_13_30 * r_2 * h5_p4 - e_3 * fs_3_143_22 * r_2 * h7_p4 + e_3 * fs_6_143_143 * r_2 * h7_p6 + e_3 * fs_1_13_30 * r_4 * h5_p4;

        pc_54[k] = - e_0 * fs_39_4_105 * h3_p3 - e_1 * fs_5_2_15 * h5_p3 + e_1 * fs_13_2_105 * r_2 * h3_p3 - e_2 * fs_51_286_2 * h7_p3 - e_2 * fs_51_286_2002 * h7_p7 + e_2 * fs_10_13_15 * r_2 * h5_p3 - e_2 * fs_13_11_105 * r_4 * h3_p3 + e_3 * fs_3_143_2 * r_2 * h7_p3 + e_3 * fs_3_143_2002 * r_2 * h7_p7 - e_3 * fs_2_39_15 * r_4 * h5_p3 + e_3 * fs_2_33_105 * r_6 * h3_p3;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[55] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54};

    for (size_t n = 0; n < 55; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
