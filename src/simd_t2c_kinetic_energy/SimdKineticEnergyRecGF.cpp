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



#include "SimdKineticEnergyRecGF.hpp"

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
compute_gf_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gf_kinetic_energy: Basis functions must be of angular momenta four and three"));
    }

    if (harmonics.size() < 7)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gf_kinetic_energy: Harmonics must reach angular momentum 7"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gf_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 63 * nvalues, 0.0);

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

            const auto ff_0 = fbase * bexp * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * bexp * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_117_2 = 117.0 / 2.0;
    const auto f_117_4 = 117.0 / 4.0;
    const auto f_120_13 = 120.0 / 13.0;
    const auto f_132_1 = 132.0;
    const auto f_13_11 = 13.0 / 11.0;
    const auto f_13_2 = 13.0 / 2.0;
    const auto f_14_33 = 14.0 / 33.0;
    const auto f_14_39 = 14.0 / 39.0;
    const auto f_150_7 = 150.0 / 7.0;
    const auto f_15_1 = 15.0;
    const auto f_165_2 = 165.0 / 2.0;
    const auto f_176_21 = 176.0 / 21.0;
    const auto f_26_1 = 26.0;
    const auto f_273_4 = 273.0 / 4.0;
    const auto f_2975_429 = 2975.0 / 429.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_2_33 = 2.0 / 33.0;
    const auto f_30_1 = 30.0;
    const auto f_350_429 = 350.0 / 429.0;
    const auto f_35_2 = 35.0 / 2.0;
    const auto f_396_7 = 396.0 / 7.0;
    const auto f_39_1 = 39.0;
    const auto f_39_11 = 39.0 / 11.0;
    const auto f_39_2 = 39.0 / 2.0;
    const auto f_39_4 = 39.0 / 4.0;
    const auto f_40_91 = 40.0 / 91.0;
    const auto f_4_11 = 4.0 / 11.0;
    const auto f_4_13 = 4.0 / 13.0;
    const auto f_52_11 = 52.0 / 11.0;
    const auto f_600_91 = 600.0 / 91.0;
    const auto f_60_13 = 60.0 / 13.0;
    const auto f_70_13 = 70.0 / 13.0;
    const auto f_78_11 = 78.0 / 11.0;
    const auto f_8_13 = 8.0 / 13.0;
    const auto f_8_21 = 8.0 / 21.0;
    const auto f_8_33 = 8.0 / 33.0;
    const auto f_91_11 = 91.0 / 11.0;
    const auto f_91_2 = 91.0 / 2.0;
    const auto fs_100_7_3 = std::sqrt(30000.0 / 49.0);
    const auto fs_100_91_7 = std::sqrt(10000.0 / 1183.0);
    const auto fs_10_13_3 = std::sqrt(300.0 / 169.0);
    const auto fs_10_143_10 = std::sqrt(1000.0 / 20449.0);
    const auto fs_10_143_55 = std::sqrt(500.0 / 1859.0);
    const auto fs_10_143_70 = std::sqrt(7000.0 / 20449.0);
    const auto fs_10_143_77 = std::sqrt(700.0 / 1859.0);
    const auto fs_10_1_3 = std::sqrt(300.0);
    const auto fs_10_273_105 = std::sqrt(500.0 / 3549.0);
    const auto fs_10_39_3 = std::sqrt(100.0 / 507.0);
    const auto fs_10_429_14 = std::sqrt(1400.0 / 184041.0);
    const auto fs_10_429_165 = std::sqrt(500.0 / 5577.0);
    const auto fs_10_429_210 = std::sqrt(7000.0 / 61347.0);
    const auto fs_10_429_231 = std::sqrt(700.0 / 5577.0);
    const auto fs_10_429_42 = std::sqrt(1400.0 / 61347.0);
    const auto fs_10_429_462 = std::sqrt(1400.0 / 5577.0);
    const auto fs_10_429_7 = std::sqrt(700.0 / 184041.0);
    const auto fs_10_429_858 = std::sqrt(200.0 / 429.0);
    const auto fs_10_7_105 = std::sqrt(1500.0 / 7.0);
    const auto fs_10_91_105 = std::sqrt(1500.0 / 1183.0);
    const auto fs_10_91_15 = std::sqrt(1500.0 / 8281.0);
    const auto fs_10_91_210 = std::sqrt(3000.0 / 1183.0);
    const auto fs_10_91_7 = std::sqrt(100.0 / 1183.0);
    const auto fs_110_91_10 = std::sqrt(121000.0 / 8281.0);
    const auto fs_117_4_3 = std::sqrt(41067.0 / 16.0);
    const auto fs_117_4_7 = std::sqrt(95823.0 / 16.0);
    const auto fs_11_273_70 = std::sqrt(1210.0 / 10647.0);
    const auto fs_125_28_10 = std::sqrt(78125.0 / 392.0);
    const auto fs_125_91_10 = std::sqrt(156250.0 / 8281.0);
    const auto fs_13_11_15 = std::sqrt(2535.0 / 121.0);
    const auto fs_13_11_21 = std::sqrt(3549.0 / 121.0);
    const auto fs_13_11_3 = std::sqrt(507.0 / 121.0);
    const auto fs_13_11_35 = std::sqrt(5915.0 / 121.0);
    const auto fs_13_11_7 = std::sqrt(1183.0 / 121.0);
    const auto fs_13_1_5 = std::sqrt(845.0);
    const auto fs_13_2_15 = std::sqrt(2535.0 / 4.0);
    const auto fs_13_2_21 = std::sqrt(3549.0 / 4.0);
    const auto fs_13_2_3 = std::sqrt(507.0 / 4.0);
    const auto fs_13_2_35 = std::sqrt(5915.0 / 4.0);
    const auto fs_13_2_7 = std::sqrt(1183.0 / 4.0);
    const auto fs_150_91_7 = std::sqrt(22500.0 / 1183.0);
    const auto fs_160_91_7 = std::sqrt(25600.0 / 1183.0);
    const auto fs_165_16_2 = std::sqrt(27225.0 / 128.0);
    const auto fs_165_16_30 = std::sqrt(408375.0 / 128.0);
    const auto fs_165_16_42 = std::sqrt(571725.0 / 128.0);
    const auto fs_165_16_6 = std::sqrt(81675.0 / 128.0);
    const auto fs_165_4_3 = std::sqrt(81675.0 / 16.0);
    const auto fs_165_8_10 = std::sqrt(136125.0 / 32.0);
    const auto fs_165_8_14 = std::sqrt(190575.0 / 32.0);
    const auto fs_165_8_15 = std::sqrt(408375.0 / 64.0);
    const auto fs_165_8_6 = std::sqrt(81675.0 / 32.0);
    const auto fs_165_8_7 = std::sqrt(190575.0 / 64.0);
    const auto fs_170_143_14 = std::sqrt(404600.0 / 20449.0);
    const auto fs_170_143_22 = std::sqrt(57800.0 / 1859.0);
    const auto fs_170_429_210 = std::sqrt(2023000.0 / 61347.0);
    const auto fs_170_429_42 = std::sqrt(404600.0 / 61347.0);
    const auto fs_170_429_6 = std::sqrt(57800.0 / 61347.0);
    const auto fs_190_91_6 = std::sqrt(216600.0 / 8281.0);
    const auto fs_198_7_3 = std::sqrt(117612.0 / 49.0);
    const auto fs_1_21_2 = std::sqrt(2.0 / 441.0);
    const auto fs_1_21_30 = std::sqrt(10.0 / 147.0);
    const auto fs_1_21_42 = std::sqrt(2.0 / 21.0);
    const auto fs_1_21_6 = std::sqrt(2.0 / 147.0);
    const auto fs_20_13_15 = std::sqrt(6000.0 / 169.0);
    const auto fs_20_13_5 = std::sqrt(2000.0 / 169.0);
    const auto fs_20_143_14 = std::sqrt(5600.0 / 20449.0);
    const auto fs_20_143_22 = std::sqrt(800.0 / 1859.0);
    const auto fs_20_273_7 = std::sqrt(400.0 / 10647.0);
    const auto fs_20_429_210 = std::sqrt(28000.0 / 61347.0);
    const auto fs_20_429_42 = std::sqrt(5600.0 / 61347.0);
    const auto fs_20_429_6 = std::sqrt(800.0 / 61347.0);
    const auto fs_20_91_21 = std::sqrt(1200.0 / 1183.0);
    const auto fs_20_91_7 = std::sqrt(400.0 / 1183.0);
    const auto fs_215_28_2 = std::sqrt(46225.0 / 392.0);
    const auto fs_215_91_2 = std::sqrt(92450.0 / 8281.0);
    const auto fs_22_21_2 = std::sqrt(968.0 / 441.0);
    const auto fs_22_21_30 = std::sqrt(4840.0 / 147.0);
    const auto fs_22_21_42 = std::sqrt(968.0 / 21.0);
    const auto fs_22_21_6 = std::sqrt(968.0 / 147.0);
    const auto fs_22_273_10 = std::sqrt(4840.0 / 74529.0);
    const auto fs_25_14_105 = std::sqrt(9375.0 / 28.0);
    const auto fs_25_273_10 = std::sqrt(6250.0 / 74529.0);
    const auto fs_25_2_3 = std::sqrt(1875.0 / 4.0);
    const auto fs_25_7_7 = std::sqrt(625.0 / 7.0);
    const auto fs_26_11_5 = std::sqrt(3380.0 / 121.0);
    const auto fs_2_11_3 = std::sqrt(12.0 / 121.0);
    const auto fs_2_11_7 = std::sqrt(28.0 / 121.0);
    const auto fs_2_21_10 = std::sqrt(40.0 / 441.0);
    const auto fs_2_21_14 = std::sqrt(8.0 / 63.0);
    const auto fs_2_21_15 = std::sqrt(20.0 / 147.0);
    const auto fs_2_21_6 = std::sqrt(8.0 / 147.0);
    const auto fs_2_21_7 = std::sqrt(4.0 / 63.0);
    const auto fs_2_273_105 = std::sqrt(20.0 / 3549.0);
    const auto fs_2_273_15 = std::sqrt(20.0 / 24843.0);
    const auto fs_2_273_210 = std::sqrt(40.0 / 3549.0);
    const auto fs_2_33_15 = std::sqrt(20.0 / 363.0);
    const auto fs_2_33_21 = std::sqrt(28.0 / 363.0);
    const auto fs_2_33_3 = std::sqrt(4.0 / 363.0);
    const auto fs_2_33_35 = std::sqrt(140.0 / 1089.0);
    const auto fs_2_33_7 = std::sqrt(28.0 / 1089.0);
    const auto fs_2_39_3 = std::sqrt(4.0 / 507.0);
    const auto fs_300_91_7 = std::sqrt(90000.0 / 1183.0);
    const auto fs_32_273_7 = std::sqrt(1024.0 / 10647.0);
    const auto fs_33_1_10 = std::sqrt(10890.0);
    const auto fs_33_1_14 = std::sqrt(15246.0);
    const auto fs_33_1_15 = std::sqrt(16335.0);
    const auto fs_33_1_6 = std::sqrt(6534.0);
    const auto fs_33_1_7 = std::sqrt(7623.0);
    const auto fs_33_2_2 = std::sqrt(1089.0 / 2.0);
    const auto fs_33_2_30 = std::sqrt(16335.0 / 2.0);
    const auto fs_33_2_42 = std::sqrt(22869.0 / 2.0);
    const auto fs_33_2_6 = std::sqrt(3267.0 / 2.0);
    const auto fs_340_429_30 = std::sqrt(1156000.0 / 61347.0);
    const auto fs_340_429_70 = std::sqrt(8092000.0 / 184041.0);
    const auto fs_38_273_6 = std::sqrt(2888.0 / 24843.0);
    const auto fs_39_11_3 = std::sqrt(4563.0 / 121.0);
    const auto fs_39_11_7 = std::sqrt(10647.0 / 121.0);
    const auto fs_39_2_3 = std::sqrt(4563.0 / 4.0);
    const auto fs_39_2_5 = std::sqrt(7605.0 / 4.0);
    const auto fs_39_2_7 = std::sqrt(10647.0 / 4.0);
    const auto fs_39_4_15 = std::sqrt(22815.0 / 16.0);
    const auto fs_39_4_21 = std::sqrt(31941.0 / 16.0);
    const auto fs_39_4_3 = std::sqrt(4563.0 / 16.0);
    const auto fs_39_4_35 = std::sqrt(53235.0 / 16.0);
    const auto fs_39_4_7 = std::sqrt(10647.0 / 16.0);
    const auto fs_400_91_3 = std::sqrt(480000.0 / 8281.0);
    const auto fs_40_13_3 = std::sqrt(4800.0 / 169.0);
    const auto fs_40_429_30 = std::sqrt(16000.0 / 61347.0);
    const auto fs_40_429_70 = std::sqrt(112000.0 / 184041.0);
    const auto fs_40_7_7 = std::sqrt(1600.0 / 7.0);
    const auto fs_40_91_105 = std::sqrt(24000.0 / 1183.0);
    const auto fs_425_429_42 = std::sqrt(2528750.0 / 61347.0);
    const auto fs_43_273_2 = std::sqrt(3698.0 / 74529.0);
    const auto fs_44_21_10 = std::sqrt(19360.0 / 441.0);
    const auto fs_44_21_14 = std::sqrt(3872.0 / 63.0);
    const auto fs_44_21_15 = std::sqrt(9680.0 / 147.0);
    const auto fs_44_21_6 = std::sqrt(3872.0 / 147.0);
    const auto fs_44_21_7 = std::sqrt(1936.0 / 63.0);
    const auto fs_45_14_21 = std::sqrt(6075.0 / 28.0);
    const auto fs_45_14_35 = std::sqrt(10125.0 / 28.0);
    const auto fs_4_21_3 = std::sqrt(16.0 / 147.0);
    const auto fs_4_273_21 = std::sqrt(16.0 / 3549.0);
    const auto fs_4_33_5 = std::sqrt(80.0 / 1089.0);
    const auto fs_4_39_15 = std::sqrt(80.0 / 507.0);
    const auto fs_4_39_5 = std::sqrt(80.0 / 1521.0);
    const auto fs_50_13_3 = std::sqrt(7500.0 / 169.0);
    const auto fs_50_429_42 = std::sqrt(35000.0 / 61347.0);
    const auto fs_50_91_105 = std::sqrt(37500.0 / 1183.0);
    const auto fs_55_14_10 = std::sqrt(15125.0 / 98.0);
    const auto fs_55_28_70 = std::sqrt(15125.0 / 56.0);
    const auto fs_55_91_70 = std::sqrt(30250.0 / 1183.0);
    const auto fs_595_429_15 = std::sqrt(1770125.0 / 61347.0);
    const auto fs_595_429_3 = std::sqrt(354025.0 / 61347.0);
    const auto fs_5_143_10 = std::sqrt(250.0 / 20449.0);
    const auto fs_5_143_110 = std::sqrt(250.0 / 1859.0);
    const auto fs_5_143_2 = std::sqrt(50.0 / 20449.0);
    const auto fs_5_143_286 = std::sqrt(50.0 / 143.0);
    const auto fs_5_14_105 = std::sqrt(375.0 / 28.0);
    const auto fs_5_14_15 = std::sqrt(375.0 / 196.0);
    const auto fs_5_14_210 = std::sqrt(375.0 / 14.0);
    const auto fs_5_1_15 = std::sqrt(375.0);
    const auto fs_5_1_5 = std::sqrt(125.0);
    const auto fs_5_2_3 = std::sqrt(75.0 / 4.0);
    const auto fs_5_429_2 = std::sqrt(50.0 / 184041.0);
    const auto fs_5_429_6006 = std::sqrt(350.0 / 429.0);
    const auto fs_5_7_21 = std::sqrt(75.0 / 7.0);
    const auto fs_5_91_30 = std::sqrt(750.0 / 8281.0);
    const auto fs_66_1_3 = std::sqrt(13068.0);
    const auto fs_6_91_21 = std::sqrt(108.0 / 1183.0);
    const auto fs_6_91_35 = std::sqrt(180.0 / 1183.0);
    const auto fs_70_429_15 = std::sqrt(24500.0 / 61347.0);
    const auto fs_70_429_3 = std::sqrt(4900.0 / 61347.0);
    const auto fs_75_14_7 = std::sqrt(5625.0 / 28.0);
    const auto fs_75_28_30 = std::sqrt(84375.0 / 392.0);
    const auto fs_75_7_7 = std::sqrt(5625.0 / 7.0);
    const auto fs_75_91_30 = std::sqrt(168750.0 / 8281.0);
    const auto fs_80_273_3 = std::sqrt(6400.0 / 24843.0);
    const auto fs_85_143_10 = std::sqrt(72250.0 / 20449.0);
    const auto fs_85_143_55 = std::sqrt(36125.0 / 1859.0);
    const auto fs_85_143_70 = std::sqrt(505750.0 / 20449.0);
    const auto fs_85_143_77 = std::sqrt(50575.0 / 1859.0);
    const auto fs_85_286_10 = std::sqrt(36125.0 / 40898.0);
    const auto fs_85_286_110 = std::sqrt(36125.0 / 3718.0);
    const auto fs_85_286_2 = std::sqrt(7225.0 / 40898.0);
    const auto fs_85_286_286 = std::sqrt(7225.0 / 286.0);
    const auto fs_85_429_14 = std::sqrt(101150.0 / 184041.0);
    const auto fs_85_429_165 = std::sqrt(36125.0 / 5577.0);
    const auto fs_85_429_210 = std::sqrt(505750.0 / 61347.0);
    const auto fs_85_429_231 = std::sqrt(50575.0 / 5577.0);
    const auto fs_85_429_42 = std::sqrt(101150.0 / 61347.0);
    const auto fs_85_429_462 = std::sqrt(101150.0 / 5577.0);
    const auto fs_85_429_7 = std::sqrt(50575.0 / 184041.0);
    const auto fs_85_429_858 = std::sqrt(14450.0 / 429.0);
    const auto fs_85_858_2 = std::sqrt(7225.0 / 368082.0);
    const auto fs_85_858_6006 = std::sqrt(50575.0 / 858.0);
    const auto fs_88_21_3 = std::sqrt(7744.0 / 147.0);
    const auto fs_8_273_105 = std::sqrt(320.0 / 3549.0);
    const auto fs_8_39_3 = std::sqrt(64.0 / 507.0);
    const auto fs_90_91_21 = std::sqrt(24300.0 / 1183.0);
    const auto fs_90_91_35 = std::sqrt(40500.0 / 1183.0);
    const auto fs_95_14_6 = std::sqrt(27075.0 / 98.0);
    const auto fs_99_14_2 = std::sqrt(9801.0 / 98.0);
    const auto fs_99_14_30 = std::sqrt(147015.0 / 98.0);
    const auto fs_99_14_42 = std::sqrt(29403.0 / 14.0);
    const auto fs_99_14_6 = std::sqrt(29403.0 / 98.0);
    const auto fs_99_7_10 = std::sqrt(98010.0 / 49.0);
    const auto fs_99_7_14 = std::sqrt(19602.0 / 7.0);
    const auto fs_99_7_15 = std::sqrt(147015.0 / 49.0);
    const auto fs_99_7_6 = std::sqrt(58806.0 / 49.0);
    const auto fs_99_7_7 = std::sqrt(9801.0 / 7.0);

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph7_p6, ph7_p7, ab_2, pc_0, pc_1 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_0[k] = - e_0 * fs_165_8_14 * h1_p1 - e_1 * fs_39_4_21 * h3_p1 + e_1 * fs_33_1_14 * r_2 * h1_p1 - e_2 * fs_5_14_210 * h5_p1 + e_2 * fs_13_2_21 * r_2 * h3_p1 - e_2 * fs_99_7_14 * r_4 * h1_p1 - e_3 * fs_85_858_2 * h7_p1 - e_3 * fs_85_858_6006 * h7_p7 + e_3 * fs_10_91_210 * r_2 * h5_p1 - e_3 * fs_13_11_21 * r_4 * h3_p1 + e_3 * fs_44_21_14 * r_6 * h1_p1 + e_4 * fs_5_429_2 * r_2 * h7_p1 + e_4 * fs_5_429_6006 * r_2 * h7_p7 - e_4 * fs_2_273_210 * r_4 * h5_p1 + e_4 * fs_2_33_21 * r_6 * h3_p1 - e_4 * fs_2_21_14 * r_8 * h1_p1;

        pc_1[k] = e_1 * fs_39_4_35 * h3_p2 + e_2 * fs_5_1_5 * h5_p2 - e_2 * fs_13_2_35 * r_2 * h3_p2 + e_3 * fs_85_286_2 * h7_p2 - e_3 * fs_85_286_286 * h7_p6 - e_3 * fs_20_13_5 * r_2 * h5_p2 + e_3 * fs_13_11_35 * r_4 * h3_p2 - e_4 * fs_5_143_2 * r_2 * h7_p2 + e_4 * fs_5_143_286 * r_2 * h7_p6 + e_4 * fs_4_39_5 * r_4 * h5_p2 - e_4 * fs_2_33_35 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_p3, ph5_m5, ph5_m4, ph5_m3, ph5_p3, ph5_p5, ph7_m5, ph7_m4, ph7_m3, ph7_p3, ph7_p5, ab_2, pc_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];

        pc_2[k] = - e_1 * fs_39_4_21 * h3_p3 - e_2 * fs_10_1_3 * h5_p3 - e_2 * fs_5_1_15 * h5_p5 + e_2 * fs_13_2_21 * r_2 * h3_p3 - e_3 * fs_85_286_10 * h7_p3 - e_3 * fs_85_286_110 * h7_p5 + e_3 * fs_40_13_3 * r_2 * h5_p3 + e_3 * fs_20_13_15 * r_2 * h5_p5 - e_3 * fs_13_11_21 * r_4 * h3_p3 + e_4 * fs_5_143_10 * r_2 * h7_p3 + e_4 * fs_5_143_110 * r_2 * h7_p5 - e_4 * fs_8_39_3 * r_4 * h5_p3 - e_4 * fs_4_39_15 * r_4 * h5_p5 + e_4 * fs_2_33_21 * r_6 * h3_p3;

        pc_3[k] = e_2 * f_30_1 * h5_m4 + e_3 * fs_85_429_165 * h7_m4 - e_3 * f_120_13 * r_2 * h5_m4 - e_4 * fs_10_429_165 * r_2 * h7_m4 + e_4 * f_8_13 * r_4 * h5_m4;

        pc_4[k] = - e_1 * fs_39_4_21 * h3_m3 + e_2 * fs_5_1_15 * h5_m5 - e_2 * fs_10_1_3 * h5_m3 + e_2 * fs_13_2_21 * r_2 * h3_m3 + e_3 * fs_85_286_110 * h7_m5 - e_3 * fs_85_286_10 * h7_m3 - e_3 * fs_20_13_15 * r_2 * h5_m5 + e_3 * fs_40_13_3 * r_2 * h5_m3 - e_3 * fs_13_11_21 * r_4 * h3_m3 - e_4 * fs_5_143_110 * r_2 * h7_m5 + e_4 * fs_5_143_10 * r_2 * h7_m3 + e_4 * fs_4_39_15 * r_4 * h5_m5 - e_4 * fs_8_39_3 * r_4 * h5_m3 + e_4 * fs_2_33_21 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m7, ph7_m6, ph7_m2, ph7_m1, ab_2, pc_5, pc_6 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_5[k] = e_1 * fs_39_4_35 * h3_m2 + e_2 * fs_5_1_5 * h5_m2 - e_2 * fs_13_2_35 * r_2 * h3_m2 + e_3 * fs_85_286_286 * h7_m6 + e_3 * fs_85_286_2 * h7_m2 - e_3 * fs_20_13_5 * r_2 * h5_m2 + e_3 * fs_13_11_35 * r_4 * h3_m2 - e_4 * fs_5_143_286 * r_2 * h7_m6 - e_4 * fs_5_143_2 * r_2 * h7_m2 + e_4 * fs_4_39_5 * r_4 * h5_m2 - e_4 * fs_2_33_35 * r_6 * h3_m2;

        pc_6[k] = - e_0 * fs_165_8_14 * h1_m1 - e_1 * fs_39_4_21 * h3_m1 + e_1 * fs_33_1_14 * r_2 * h1_m1 - e_2 * fs_5_14_210 * h5_m1 + e_2 * fs_13_2_21 * r_2 * h3_m1 - e_2 * fs_99_7_14 * r_4 * h1_m1 + e_3 * fs_85_858_6006 * h7_m7 - e_3 * fs_85_858_2 * h7_m1 + e_3 * fs_10_91_210 * r_2 * h5_m1 - e_3 * fs_13_11_21 * r_4 * h3_m1 + e_3 * fs_44_21_14 * r_6 * h1_m1 - e_4 * fs_5_429_6006 * r_2 * h7_m7 + e_4 * fs_5_429_2 * r_2 * h7_m1 - e_4 * fs_2_273_210 * r_4 * h5_m1 + e_4 * fs_2_33_21 * r_6 * h3_m1 - e_4 * fs_2_21_14 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p5, ph7_0, ph7_p1, ph7_p5, ph7_p6, ab_2, pc_7, pc_8 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_7[k] = - e_0 * fs_165_8_7 * h1_0 - e_1 * fs_117_4_7 * h3_0 + e_1 * fs_33_1_7 * r_2 * h1_0 - e_2 * fs_75_14_7 * h5_0 + e_2 * fs_39_2_7 * r_2 * h3_0 - e_2 * fs_99_7_7 * r_4 * h1_0 - e_3 * fs_85_429_7 * h7_0 - e_3 * fs_85_429_858 * h7_p6 + e_3 * fs_150_91_7 * r_2 * h5_0 - e_3 * fs_39_11_7 * r_4 * h3_0 + e_3 * fs_44_21_7 * r_6 * h1_0 + e_4 * fs_10_429_7 * r_2 * h7_0 + e_4 * fs_10_429_858 * r_2 * h7_p6 - e_4 * fs_10_91_7 * r_4 * h5_0 + e_4 * fs_2_11_7 * r_6 * h3_0 - e_4 * fs_2_21_7 * r_8 * h1_0;

        pc_8[k] = - e_0 * fs_165_16_42 * h1_p1 + e_1 * fs_39_4_7 * h3_p1 + e_1 * fs_33_2_42 * r_2 * h1_p1 + e_2 * fs_55_28_70 * h5_p1 + e_2 * fs_25_2_3 * h5_p5 - e_2 * fs_13_2_7 * r_2 * h3_p1 - e_2 * fs_99_14_42 * r_4 * h1_p1 + e_3 * fs_170_429_6 * h7_p1 - e_3 * fs_170_143_22 * h7_p5 - e_3 * fs_55_91_70 * r_2 * h5_p1 - e_3 * fs_50_13_3 * r_2 * h5_p5 + e_3 * fs_13_11_7 * r_4 * h3_p1 + e_3 * fs_22_21_42 * r_6 * h1_p1 - e_4 * fs_20_429_6 * r_2 * h7_p1 + e_4 * fs_20_143_22 * r_2 * h7_p5 + e_4 * fs_11_273_70 * r_4 * h5_p1 + e_4 * fs_10_39_3 * r_4 * h5_p5 - e_4 * fs_2_33_7 * r_6 * h3_p1 - e_4 * fs_1_21_42 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph3_p2, ph5_m4, ph5_m3, ph5_m2, ph5_p2, ph5_p4, ph7_m4, ph7_m3, ph7_m2, ph7_p2, ph7_p4, ab_2, pc_9, pc_10, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];

        pc_9[k] = e_1 * fs_39_4_7 * h3_p2 - e_2 * f_35_2 * h5_p2 + e_2 * fs_5_2_3 * h5_p4 - e_2 * fs_13_2_7 * r_2 * h3_p2 - e_3 * fs_85_143_10 * h7_p2 - e_3 * fs_85_143_55 * h7_p4 + e_3 * f_70_13 * r_2 * h5_p2 - e_3 * fs_10_13_3 * r_2 * h5_p4 + e_3 * fs_13_11_7 * r_4 * h3_p2 + e_4 * fs_10_143_10 * r_2 * h7_p2 + e_4 * fs_10_143_55 * r_2 * h7_p4 - e_4 * f_14_39 * r_4 * h5_p2 + e_4 * fs_2_39_3 * r_4 * h5_p4 - e_4 * fs_2_33_7 * r_6 * h3_p2;

        pc_10[k] = - e_1 * fs_117_4_7 * h3_m3 + e_2 * f_15_1 * h5_m3 + e_2 * fs_39_2_7 * r_2 * h3_m3 + e_3 * fs_340_429_30 * h7_m3 - e_3 * f_60_13 * r_2 * h5_m3 - e_3 * fs_39_11_7 * r_4 * h3_m3 - e_4 * fs_40_429_30 * r_2 * h7_m3 + e_4 * f_4_13 * r_4 * h5_m3 + e_4 * fs_2_11_7 * r_6 * h3_m3;

        pc_11[k] = e_1 * fs_39_4_7 * h3_m2 - e_2 * fs_5_2_3 * h5_m4 - e_2 * f_35_2 * h5_m2 - e_2 * fs_13_2_7 * r_2 * h3_m2 + e_3 * fs_85_143_55 * h7_m4 - e_3 * fs_85_143_10 * h7_m2 + e_3 * fs_10_13_3 * r_2 * h5_m4 + e_3 * f_70_13 * r_2 * h5_m2 + e_3 * fs_13_11_7 * r_4 * h3_m2 - e_4 * fs_10_143_55 * r_2 * h7_m4 + e_4 * fs_10_143_10 * r_2 * h7_m2 - e_4 * fs_2_39_3 * r_4 * h5_m4 - e_4 * f_14_39 * r_4 * h5_m2 - e_4 * fs_2_33_7 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m6, ph7_m5, ph7_m1, ab_2, pc_12, pc_13 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];

        pc_12[k] = - e_0 * fs_165_16_42 * h1_m1 + e_1 * fs_39_4_7 * h3_m1 + e_1 * fs_33_2_42 * r_2 * h1_m1 - e_2 * fs_25_2_3 * h5_m5 + e_2 * fs_55_28_70 * h5_m1 - e_2 * fs_13_2_7 * r_2 * h3_m1 - e_2 * fs_99_14_42 * r_4 * h1_m1 + e_3 * fs_170_143_22 * h7_m5 + e_3 * fs_170_429_6 * h7_m1 + e_3 * fs_50_13_3 * r_2 * h5_m5 - e_3 * fs_55_91_70 * r_2 * h5_m1 + e_3 * fs_13_11_7 * r_4 * h3_m1 + e_3 * fs_22_21_42 * r_6 * h1_m1 - e_4 * fs_20_143_22 * r_2 * h7_m5 - e_4 * fs_20_429_6 * r_2 * h7_m1 - e_4 * fs_10_39_3 * r_4 * h5_m5 + e_4 * fs_11_273_70 * r_4 * h5_m1 - e_4 * fs_2_33_7 * r_6 * h3_m1 - e_4 * fs_1_21_42 * r_8 * h1_m1;

        pc_13[k] = e_3 * fs_85_429_858 * h7_m6 - e_4 * fs_10_429_858 * r_2 * h7_m6;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p4, ph5_p5, ph7_0, ph7_p1, ph7_p4, ph7_p5, ab_2, pc_14, pc_15 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_14[k] = e_0 * fs_165_16_2 * h1_p1 + e_1 * fs_117_4_3 * h3_p1 - e_1 * fs_33_2_2 * r_2 * h1_p1 + e_2 * fs_75_28_30 * h5_p1 - e_2 * fs_75_14_7 * h5_p5 - e_2 * fs_39_2_3 * r_2 * h3_p1 + e_2 * fs_99_14_2 * r_4 * h1_p1 + e_3 * fs_85_429_14 * h7_p1 - e_3 * fs_85_429_462 * h7_p5 - e_3 * fs_75_91_30 * r_2 * h5_p1 + e_3 * fs_150_91_7 * r_2 * h5_p5 + e_3 * fs_39_11_3 * r_4 * h3_p1 - e_3 * fs_22_21_2 * r_6 * h1_p1 - e_4 * fs_10_429_14 * r_2 * h7_p1 + e_4 * fs_10_429_462 * r_2 * h7_p5 + e_4 * fs_5_91_30 * r_4 * h5_p1 - e_4 * fs_10_91_7 * r_4 * h5_p5 - e_4 * fs_2_11_3 * r_6 * h3_p1 + e_4 * fs_1_21_2 * r_8 * h1_p1;

        pc_15[k] = - e_0 * fs_165_4_3 * h1_0 - e_1 * fs_39_4_3 * h3_0 + e_1 * fs_66_1_3 * r_2 * h1_0 + e_2 * fs_100_7_3 * h5_0 + e_2 * fs_10_7_105 * h5_p4 + e_2 * fs_13_2_3 * r_2 * h3_0 - e_2 * fs_198_7_3 * r_4 * h1_0 + e_3 * fs_595_429_3 * h7_0 - e_3 * fs_85_143_77 * h7_p4 - e_3 * fs_400_91_3 * r_2 * h5_0 - e_3 * fs_40_91_105 * r_2 * h5_p4 - e_3 * fs_13_11_3 * r_4 * h3_0 + e_3 * fs_88_21_3 * r_6 * h1_0 - e_4 * fs_70_429_3 * r_2 * h7_0 + e_4 * fs_10_143_77 * r_2 * h7_p4 + e_4 * fs_80_273_3 * r_4 * h5_0 + e_4 * fs_8_273_105 * r_4 * h5_p4 + e_4 * fs_2_33_3 * r_6 * h3_0 - e_4 * fs_4_21_3 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_m2, ph3_p1, ph3_p3, ph5_m2, ph5_p1, ph5_p3, ph7_m2, ph7_p1, ph7_p3, ab_2, pc_16, pc_17 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];

        pc_16[k] = - e_0 * fs_165_16_30 * h1_p1 + e_1 * fs_39_2_5 * h3_p1 - e_1 * fs_117_4_3 * h3_p3 + e_1 * fs_33_2_30 * r_2 * h1_p1 - e_2 * fs_215_28_2 * h5_p1 + e_2 * fs_45_14_21 * h5_p3 - e_2 * fs_13_1_5 * r_2 * h3_p1 + e_2 * fs_39_2_3 * r_2 * h3_p3 - e_2 * fs_99_14_30 * r_4 * h1_p1 - e_3 * fs_85_429_210 * h7_p1 - e_3 * fs_85_143_70 * h7_p3 + e_3 * fs_215_91_2 * r_2 * h5_p1 - e_3 * fs_90_91_21 * r_2 * h5_p3 + e_3 * fs_26_11_5 * r_4 * h3_p1 - e_3 * fs_39_11_3 * r_4 * h3_p3 + e_3 * fs_22_21_30 * r_6 * h1_p1 + e_4 * fs_10_429_210 * r_2 * h7_p1 + e_4 * fs_10_143_70 * r_2 * h7_p3 - e_4 * fs_43_273_2 * r_4 * h5_p1 + e_4 * fs_6_91_21 * r_4 * h5_p3 - e_4 * fs_4_33_5 * r_6 * h3_p1 + e_4 * fs_2_11_3 * r_6 * h3_p3 - e_4 * fs_1_21_30 * r_8 * h1_p1;

        pc_17[k] = - e_1 * fs_39_4_3 * h3_m2 - e_2 * fs_5_7_21 * h5_m2 + e_2 * fs_13_2_3 * r_2 * h3_m2 + e_3 * fs_170_429_210 * h7_m2 + e_3 * fs_20_91_21 * r_2 * h5_m2 - e_3 * fs_13_11_3 * r_4 * h3_m2 - e_4 * fs_20_429_210 * r_2 * h7_m2 - e_4 * fs_4_273_21 * r_4 * h5_m2 + e_4 * fs_2_33_3 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph5_m5, ph5_m4, ph5_m3, ph5_m1, ph7_m5, ph7_m4, ph7_m3, ph7_m1, ab_2, pc_18, pc_19, pc_20 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];

        pc_18[k] = - e_0 * fs_165_16_30 * h1_m1 + e_1 * fs_117_4_3 * h3_m3 + e_1 * fs_39_2_5 * h3_m1 + e_1 * fs_33_2_30 * r_2 * h1_m1 - e_2 * fs_45_14_21 * h5_m3 - e_2 * fs_215_28_2 * h5_m1 - e_2 * fs_39_2_3 * r_2 * h3_m3 - e_2 * fs_13_1_5 * r_2 * h3_m1 - e_2 * fs_99_14_30 * r_4 * h1_m1 + e_3 * fs_85_143_70 * h7_m3 - e_3 * fs_85_429_210 * h7_m1 + e_3 * fs_90_91_21 * r_2 * h5_m3 + e_3 * fs_215_91_2 * r_2 * h5_m1 + e_3 * fs_39_11_3 * r_4 * h3_m3 + e_3 * fs_26_11_5 * r_4 * h3_m1 + e_3 * fs_22_21_30 * r_6 * h1_m1 - e_4 * fs_10_143_70 * r_2 * h7_m3 + e_4 * fs_10_429_210 * r_2 * h7_m1 - e_4 * fs_6_91_21 * r_4 * h5_m3 - e_4 * fs_43_273_2 * r_4 * h5_m1 - e_4 * fs_2_11_3 * r_6 * h3_m3 - e_4 * fs_4_33_5 * r_6 * h3_m1 - e_4 * fs_1_21_30 * r_8 * h1_m1;

        pc_19[k] = - e_2 * fs_10_7_105 * h5_m4 + e_3 * fs_85_143_77 * h7_m4 + e_3 * fs_40_91_105 * r_2 * h5_m4 - e_4 * fs_10_143_77 * r_2 * h7_m4 - e_4 * fs_8_273_105 * r_4 * h5_m4;

        pc_20[k] = - e_0 * fs_165_16_2 * h1_m1 - e_1 * fs_117_4_3 * h3_m1 + e_1 * fs_33_2_2 * r_2 * h1_m1 + e_2 * fs_75_14_7 * h5_m5 - e_2 * fs_75_28_30 * h5_m1 + e_2 * fs_39_2_3 * r_2 * h3_m1 - e_2 * fs_99_14_2 * r_4 * h1_m1 + e_3 * fs_85_429_462 * h7_m5 - e_3 * fs_85_429_14 * h7_m1 - e_3 * fs_150_91_7 * r_2 * h5_m5 + e_3 * fs_75_91_30 * r_2 * h5_m1 - e_3 * fs_39_11_3 * r_4 * h3_m1 + e_3 * fs_22_21_2 * r_6 * h1_m1 - e_4 * fs_10_429_462 * r_2 * h7_m5 + e_4 * fs_10_429_14 * r_2 * h7_m1 + e_4 * fs_10_91_7 * r_4 * h5_m5 - e_4 * fs_5_91_30 * r_4 * h5_m1 + e_4 * fs_2_11_3 * r_6 * h3_m1 - e_4 * fs_1_21_2 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p2, ph3_p3, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph7_p1, ph7_p2, ph7_p3, ph7_p4, ab_2, pc_21, pc_22 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_21[k] = - e_1 * fs_39_4_15 * h3_p2 - e_2 * fs_25_14_105 * h5_p2 - e_2 * fs_45_14_35 * h5_p4 + e_2 * fs_13_2_15 * r_2 * h3_p2 - e_3 * fs_85_429_42 * h7_p2 - e_3 * fs_85_429_231 * h7_p4 + e_3 * fs_50_91_105 * r_2 * h5_p2 + e_3 * fs_90_91_35 * r_2 * h5_p4 - e_3 * fs_13_11_15 * r_4 * h3_p2 + e_4 * fs_10_429_42 * r_2 * h7_p2 + e_4 * fs_10_429_231 * r_2 * h7_p4 - e_4 * fs_10_273_105 * r_4 * h5_p2 - e_4 * fs_6_91_35 * r_4 * h5_p4 + e_4 * fs_2_33_15 * r_6 * h3_p2;

        pc_22[k] = e_0 * fs_165_16_6 * h1_p1 + e_1 * f_39_1 * h3_p1 + e_1 * fs_39_4_15 * h3_p3 - e_1 * fs_33_2_6 * r_2 * h1_p1 - e_2 * fs_125_28_10 * h5_p1 + e_2 * fs_5_14_105 * h5_p3 - e_2 * f_26_1 * r_2 * h3_p1 - e_2 * fs_13_2_15 * r_2 * h3_p3 + e_2 * fs_99_14_6 * r_4 * h1_p1 - e_3 * fs_170_429_42 * h7_p1 - e_3 * fs_170_143_14 * h7_p3 + e_3 * fs_125_91_10 * r_2 * h5_p1 - e_3 * fs_10_91_105 * r_2 * h5_p3 + e_3 * f_52_11 * r_4 * h3_p1 + e_3 * fs_13_11_15 * r_4 * h3_p3 - e_3 * fs_22_21_6 * r_6 * h1_p1 + e_4 * fs_20_429_42 * r_2 * h7_p1 + e_4 * fs_20_143_14 * r_2 * h7_p3 - e_4 * fs_25_273_10 * r_4 * h5_p1 + e_4 * fs_2_273_105 * r_4 * h5_p3 - e_4 * f_8_33 * r_6 * h3_p1 - e_4 * fs_2_33_15 * r_6 * h3_p3 + e_4 * fs_1_21_6 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph1_0, ph3_m1, ph3_0, ph3_p2, ph5_m1, ph5_0, ph5_p2, ph7_m1, ph7_0, ph7_p2, ab_2, pc_23, pc_24 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];

        pc_23[k] = - e_0 * fs_165_8_15 * h1_0 + e_1 * fs_39_4_15 * h3_0 - e_1 * f_39_1 * h3_p2 + e_1 * fs_33_1_15 * r_2 * h1_0 - e_2 * fs_5_14_15 * h5_0 + e_2 * fs_40_7_7 * h5_p2 - e_2 * fs_13_2_15 * r_2 * h3_0 + e_2 * f_26_1 * r_2 * h3_p2 - e_2 * fs_99_7_15 * r_4 * h1_0 - e_3 * fs_595_429_15 * h7_0 - e_3 * fs_85_143_70 * h7_p2 + e_3 * fs_10_91_15 * r_2 * h5_0 - e_3 * fs_160_91_7 * r_2 * h5_p2 + e_3 * fs_13_11_15 * r_4 * h3_0 - e_3 * f_52_11 * r_4 * h3_p2 + e_3 * fs_44_21_15 * r_6 * h1_0 + e_4 * fs_70_429_15 * r_2 * h7_0 + e_4 * fs_10_143_70 * r_2 * h7_p2 - e_4 * fs_2_273_15 * r_4 * h5_0 + e_4 * fs_32_273_7 * r_4 * h5_p2 - e_4 * fs_2_33_15 * r_6 * h3_0 + e_4 * f_8_33 * r_6 * h3_p2 - e_4 * fs_2_21_15 * r_8 * h1_0;

        pc_24[k] = - e_0 * fs_165_8_10 * h1_m1 + e_1 * fs_39_4_15 * h3_m1 + e_1 * fs_33_1_10 * r_2 * h1_m1 - e_2 * fs_95_14_6 * h5_m1 - e_2 * fs_13_2_15 * r_2 * h3_m1 - e_2 * fs_99_7_10 * r_4 * h1_m1 + e_3 * fs_340_429_70 * h7_m1 + e_3 * fs_190_91_6 * r_2 * h5_m1 + e_3 * fs_13_11_15 * r_4 * h3_m1 + e_3 * fs_44_21_10 * r_6 * h1_m1 - e_4 * fs_40_429_70 * r_2 * h7_m1 - e_4 * fs_38_273_6 * r_4 * h5_m1 - e_4 * fs_2_33_15 * r_6 * h3_m1 - e_4 * fs_2_21_10 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m2, ph3_m1, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ab_2, pc_25, pc_26, pc_27, pc_28, pc_29, pc_30 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
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

        pc_25[k] = e_1 * f_39_1 * h3_m2 - e_2 * fs_40_7_7 * h5_m2 - e_2 * f_26_1 * r_2 * h3_m2 + e_3 * fs_85_143_70 * h7_m2 + e_3 * fs_160_91_7 * r_2 * h5_m2 + e_3 * f_52_11 * r_4 * h3_m2 - e_4 * fs_10_143_70 * r_2 * h7_m2 - e_4 * fs_32_273_7 * r_4 * h5_m2 - e_4 * f_8_33 * r_6 * h3_m2;

        pc_26[k] = - e_0 * fs_165_16_6 * h1_m1 - e_1 * fs_39_4_15 * h3_m3 - e_1 * f_39_1 * h3_m1 + e_1 * fs_33_2_6 * r_2 * h1_m1 - e_2 * fs_5_14_105 * h5_m3 + e_2 * fs_125_28_10 * h5_m1 + e_2 * fs_13_2_15 * r_2 * h3_m3 + e_2 * f_26_1 * r_2 * h3_m1 - e_2 * fs_99_14_6 * r_4 * h1_m1 + e_3 * fs_170_143_14 * h7_m3 + e_3 * fs_170_429_42 * h7_m1 + e_3 * fs_10_91_105 * r_2 * h5_m3 - e_3 * fs_125_91_10 * r_2 * h5_m1 - e_3 * fs_13_11_15 * r_4 * h3_m3 - e_3 * f_52_11 * r_4 * h3_m1 + e_3 * fs_22_21_6 * r_6 * h1_m1 - e_4 * fs_20_143_14 * r_2 * h7_m3 - e_4 * fs_20_429_42 * r_2 * h7_m1 - e_4 * fs_2_273_105 * r_4 * h5_m3 + e_4 * fs_25_273_10 * r_4 * h5_m1 + e_4 * fs_2_33_15 * r_6 * h3_m3 + e_4 * f_8_33 * r_6 * h3_m1 - e_4 * fs_1_21_6 * r_8 * h1_m1;

        pc_27[k] = e_1 * fs_39_4_15 * h3_m2 + e_2 * fs_45_14_35 * h5_m4 + e_2 * fs_25_14_105 * h5_m2 - e_2 * fs_13_2_15 * r_2 * h3_m2 + e_3 * fs_85_429_231 * h7_m4 + e_3 * fs_85_429_42 * h7_m2 - e_3 * fs_90_91_35 * r_2 * h5_m4 - e_3 * fs_50_91_105 * r_2 * h5_m2 + e_3 * fs_13_11_15 * r_4 * h3_m2 - e_4 * fs_10_429_231 * r_2 * h7_m4 - e_4 * fs_10_429_42 * r_2 * h7_m2 + e_4 * fs_6_91_35 * r_4 * h5_m4 + e_4 * fs_10_273_105 * r_4 * h5_m2 - e_4 * fs_2_33_15 * r_6 * h3_m2;

        pc_28[k] = e_1 * f_117_4 * h3_m3 + e_2 * fs_75_7_7 * h5_m3 - e_2 * f_39_2 * r_2 * h3_m3 + e_3 * fs_85_429_210 * h7_m3 - e_3 * fs_300_91_7 * r_2 * h5_m3 + e_3 * f_39_11 * r_4 * h3_m3 - e_4 * fs_10_429_210 * r_2 * h7_m3 + e_4 * fs_20_91_7 * r_4 * h5_m3 - e_4 * f_2_11 * r_6 * h3_m3;

        pc_29[k] = - e_1 * f_273_4 * h3_m2 + e_2 * fs_25_7_7 * h5_m2 + e_2 * f_91_2 * r_2 * h3_m2 + e_3 * fs_85_143_70 * h7_m2 - e_3 * fs_100_91_7 * r_2 * h5_m2 - e_3 * f_91_11 * r_4 * h3_m2 - e_4 * fs_10_143_70 * r_2 * h7_m2 + e_4 * fs_20_273_7 * r_4 * h5_m2 + e_4 * f_14_33 * r_6 * h3_m2;

        pc_30[k] = e_0 * fs_165_8_6 * h1_m1 + e_1 * f_39_4 * h3_m1 - e_1 * fs_33_1_6 * r_2 * h1_m1 - e_2 * fs_55_14_10 * h5_m1 - e_2 * f_13_2 * r_2 * h3_m1 + e_2 * fs_99_7_6 * r_4 * h1_m1 + e_3 * fs_425_429_42 * h7_m1 + e_3 * fs_110_91_10 * r_2 * h5_m1 + e_3 * f_13_11 * r_4 * h3_m1 - e_3 * fs_44_21_6 * r_6 * h1_m1 - e_4 * fs_50_429_42 * r_2 * h7_m1 - e_4 * fs_22_273_10 * r_4 * h5_m1 - e_4 * f_2_33 * r_6 * h3_m1 + e_4 * fs_2_21_6 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph3_p2, ph5_0, ph5_p1, ph5_p2, ph7_0, ph7_p1, ph7_p2, ab_2, pc_31, pc_32, pc_33 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_31[k] = - e_0 * f_165_2 * h1_0 + e_1 * f_117_2 * h3_0 + e_1 * f_132_1 * r_2 * h1_0 - e_2 * f_150_7 * h5_0 - e_2 * f_39_1 * r_2 * h3_0 - e_2 * f_396_7 * r_4 * h1_0 + e_3 * f_2975_429 * h7_0 + e_3 * f_600_91 * r_2 * h5_0 + e_3 * f_78_11 * r_4 * h3_0 + e_3 * f_176_21 * r_6 * h1_0 - e_4 * f_350_429 * r_2 * h7_0 - e_4 * f_40_91 * r_4 * h5_0 - e_4 * f_4_11 * r_6 * h3_0 - e_4 * f_8_21 * r_8 * h1_0;

        pc_32[k] = e_0 * fs_165_8_6 * h1_p1 + e_1 * f_39_4 * h3_p1 - e_1 * fs_33_1_6 * r_2 * h1_p1 - e_2 * fs_55_14_10 * h5_p1 - e_2 * f_13_2 * r_2 * h3_p1 + e_2 * fs_99_7_6 * r_4 * h1_p1 + e_3 * fs_425_429_42 * h7_p1 + e_3 * fs_110_91_10 * r_2 * h5_p1 + e_3 * f_13_11 * r_4 * h3_p1 - e_3 * fs_44_21_6 * r_6 * h1_p1 - e_4 * fs_50_429_42 * r_2 * h7_p1 - e_4 * fs_22_273_10 * r_4 * h5_p1 - e_4 * f_2_33 * r_6 * h3_p1 + e_4 * fs_2_21_6 * r_8 * h1_p1;

        pc_33[k] = - e_1 * f_273_4 * h3_p2 + e_2 * fs_25_7_7 * h5_p2 + e_2 * f_91_2 * r_2 * h3_p2 + e_3 * fs_85_143_70 * h7_p2 - e_3 * fs_100_91_7 * r_2 * h5_p2 - e_3 * f_91_11 * r_4 * h3_p2 - e_4 * fs_10_143_70 * r_2 * h7_p2 + e_4 * fs_20_273_7 * r_4 * h5_p2 + e_4 * f_14_33 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph3_p3, ph5_m4, ph5_m2, ph5_p3, ph7_m4, ph7_m2, ph7_p3, ab_2, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];

        pc_34[k] = e_1 * f_117_4 * h3_p3 + e_2 * fs_75_7_7 * h5_p3 - e_2 * f_39_2 * r_2 * h3_p3 + e_3 * fs_85_429_210 * h7_p3 - e_3 * fs_300_91_7 * r_2 * h5_p3 + e_3 * f_39_11 * r_4 * h3_p3 - e_4 * fs_10_429_210 * r_2 * h7_p3 + e_4 * fs_20_91_7 * r_4 * h5_p3 - e_4 * f_2_11 * r_6 * h3_p3;

        pc_35[k] = - e_1 * fs_39_4_15 * h3_m2 + e_2 * fs_45_14_35 * h5_m4 - e_2 * fs_25_14_105 * h5_m2 + e_2 * fs_13_2_15 * r_2 * h3_m2 + e_3 * fs_85_429_231 * h7_m4 - e_3 * fs_85_429_42 * h7_m2 - e_3 * fs_90_91_35 * r_2 * h5_m4 + e_3 * fs_50_91_105 * r_2 * h5_m2 - e_3 * fs_13_11_15 * r_4 * h3_m2 - e_4 * fs_10_429_231 * r_2 * h7_m4 + e_4 * fs_10_429_42 * r_2 * h7_m2 + e_4 * fs_6_91_35 * r_4 * h5_m4 - e_4 * fs_10_273_105 * r_4 * h5_m2 + e_4 * fs_2_33_15 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m2, ph3_m1, ph5_m3, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ab_2, pc_36, pc_37 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_36[k] = e_0 * fs_165_16_6 * h1_m1 - e_1 * fs_39_4_15 * h3_m3 + e_1 * f_39_1 * h3_m1 - e_1 * fs_33_2_6 * r_2 * h1_m1 - e_2 * fs_5_14_105 * h5_m3 - e_2 * fs_125_28_10 * h5_m1 + e_2 * fs_13_2_15 * r_2 * h3_m3 - e_2 * f_26_1 * r_2 * h3_m1 + e_2 * fs_99_14_6 * r_4 * h1_m1 + e_3 * fs_170_143_14 * h7_m3 - e_3 * fs_170_429_42 * h7_m1 + e_3 * fs_10_91_105 * r_2 * h5_m3 + e_3 * fs_125_91_10 * r_2 * h5_m1 - e_3 * fs_13_11_15 * r_4 * h3_m3 + e_3 * f_52_11 * r_4 * h3_m1 - e_3 * fs_22_21_6 * r_6 * h1_m1 - e_4 * fs_20_143_14 * r_2 * h7_m3 + e_4 * fs_20_429_42 * r_2 * h7_m1 - e_4 * fs_2_273_105 * r_4 * h5_m3 - e_4 * fs_25_273_10 * r_4 * h5_m1 + e_4 * fs_2_33_15 * r_6 * h3_m3 - e_4 * f_8_33 * r_6 * h3_m1 + e_4 * fs_1_21_6 * r_8 * h1_m1;

        pc_37[k] = e_1 * f_39_1 * h3_m2 - e_2 * fs_40_7_7 * h5_m2 - e_2 * f_26_1 * r_2 * h3_m2 + e_3 * fs_85_143_70 * h7_m2 + e_3 * fs_160_91_7 * r_2 * h5_m2 + e_3 * f_52_11 * r_4 * h3_m2 - e_4 * fs_10_143_70 * r_2 * h7_m2 - e_4 * fs_32_273_7 * r_4 * h5_m2 - e_4 * f_8_33 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph3_p2, ph5_0, ph5_p1, ph5_p2, ph7_0, ph7_p1, ph7_p2, ab_2, pc_38, pc_39 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_38[k] = - e_0 * fs_165_8_10 * h1_p1 + e_1 * fs_39_4_15 * h3_p1 + e_1 * fs_33_1_10 * r_2 * h1_p1 - e_2 * fs_95_14_6 * h5_p1 - e_2 * fs_13_2_15 * r_2 * h3_p1 - e_2 * fs_99_7_10 * r_4 * h1_p1 + e_3 * fs_340_429_70 * h7_p1 + e_3 * fs_190_91_6 * r_2 * h5_p1 + e_3 * fs_13_11_15 * r_4 * h3_p1 + e_3 * fs_44_21_10 * r_6 * h1_p1 - e_4 * fs_40_429_70 * r_2 * h7_p1 - e_4 * fs_38_273_6 * r_4 * h5_p1 - e_4 * fs_2_33_15 * r_6 * h3_p1 - e_4 * fs_2_21_10 * r_8 * h1_p1;

        pc_39[k] = - e_0 * fs_165_8_15 * h1_0 + e_1 * fs_39_4_15 * h3_0 + e_1 * f_39_1 * h3_p2 + e_1 * fs_33_1_15 * r_2 * h1_0 - e_2 * fs_5_14_15 * h5_0 - e_2 * fs_40_7_7 * h5_p2 - e_2 * fs_13_2_15 * r_2 * h3_0 - e_2 * f_26_1 * r_2 * h3_p2 - e_2 * fs_99_7_15 * r_4 * h1_0 - e_3 * fs_595_429_15 * h7_0 + e_3 * fs_85_143_70 * h7_p2 + e_3 * fs_10_91_15 * r_2 * h5_0 + e_3 * fs_160_91_7 * r_2 * h5_p2 + e_3 * fs_13_11_15 * r_4 * h3_0 + e_3 * f_52_11 * r_4 * h3_p2 + e_3 * fs_44_21_15 * r_6 * h1_0 + e_4 * fs_70_429_15 * r_2 * h7_0 - e_4 * fs_10_143_70 * r_2 * h7_p2 - e_4 * fs_2_273_15 * r_4 * h5_0 - e_4 * fs_32_273_7 * r_4 * h5_p2 - e_4 * fs_2_33_15 * r_6 * h3_0 - e_4 * f_8_33 * r_6 * h3_p2 - e_4 * fs_2_21_15 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p2, ph3_p3, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph7_p1, ph7_p2, ph7_p3, ph7_p4, ab_2, pc_40, pc_41 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_40[k] = e_0 * fs_165_16_6 * h1_p1 + e_1 * f_39_1 * h3_p1 - e_1 * fs_39_4_15 * h3_p3 - e_1 * fs_33_2_6 * r_2 * h1_p1 - e_2 * fs_125_28_10 * h5_p1 - e_2 * fs_5_14_105 * h5_p3 - e_2 * f_26_1 * r_2 * h3_p1 + e_2 * fs_13_2_15 * r_2 * h3_p3 + e_2 * fs_99_14_6 * r_4 * h1_p1 - e_3 * fs_170_429_42 * h7_p1 + e_3 * fs_170_143_14 * h7_p3 + e_3 * fs_125_91_10 * r_2 * h5_p1 + e_3 * fs_10_91_105 * r_2 * h5_p3 + e_3 * f_52_11 * r_4 * h3_p1 - e_3 * fs_13_11_15 * r_4 * h3_p3 - e_3 * fs_22_21_6 * r_6 * h1_p1 + e_4 * fs_20_429_42 * r_2 * h7_p1 - e_4 * fs_20_143_14 * r_2 * h7_p3 - e_4 * fs_25_273_10 * r_4 * h5_p1 - e_4 * fs_2_273_105 * r_4 * h5_p3 - e_4 * f_8_33 * r_6 * h3_p1 + e_4 * fs_2_33_15 * r_6 * h3_p3 + e_4 * fs_1_21_6 * r_8 * h1_p1;

        pc_41[k] = - e_1 * fs_39_4_15 * h3_p2 - e_2 * fs_25_14_105 * h5_p2 + e_2 * fs_45_14_35 * h5_p4 + e_2 * fs_13_2_15 * r_2 * h3_p2 - e_3 * fs_85_429_42 * h7_p2 + e_3 * fs_85_429_231 * h7_p4 + e_3 * fs_50_91_105 * r_2 * h5_p2 - e_3 * fs_90_91_35 * r_2 * h5_p4 - e_3 * fs_13_11_15 * r_4 * h3_p2 + e_4 * fs_10_429_42 * r_2 * h7_p2 - e_4 * fs_10_429_231 * r_2 * h7_p4 - e_4 * fs_10_273_105 * r_4 * h5_p2 + e_4 * fs_6_91_35 * r_4 * h5_p4 + e_4 * fs_2_33_15 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph5_m5, ph5_m4, ph5_m3, ph5_m1, ph7_m5, ph7_m4, ph7_m3, ph7_m1, ab_2, pc_42, pc_43, pc_44 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];

        pc_42[k] = e_0 * fs_165_16_2 * h1_m1 + e_1 * fs_117_4_3 * h3_m1 - e_1 * fs_33_2_2 * r_2 * h1_m1 + e_2 * fs_75_14_7 * h5_m5 + e_2 * fs_75_28_30 * h5_m1 - e_2 * fs_39_2_3 * r_2 * h3_m1 + e_2 * fs_99_14_2 * r_4 * h1_m1 + e_3 * fs_85_429_462 * h7_m5 + e_3 * fs_85_429_14 * h7_m1 - e_3 * fs_150_91_7 * r_2 * h5_m5 - e_3 * fs_75_91_30 * r_2 * h5_m1 + e_3 * fs_39_11_3 * r_4 * h3_m1 - e_3 * fs_22_21_2 * r_6 * h1_m1 - e_4 * fs_10_429_462 * r_2 * h7_m5 - e_4 * fs_10_429_14 * r_2 * h7_m1 + e_4 * fs_10_91_7 * r_4 * h5_m5 + e_4 * fs_5_91_30 * r_4 * h5_m1 - e_4 * fs_2_11_3 * r_6 * h3_m1 + e_4 * fs_1_21_2 * r_8 * h1_m1;

        pc_43[k] = - e_2 * fs_10_7_105 * h5_m4 + e_3 * fs_85_143_77 * h7_m4 + e_3 * fs_40_91_105 * r_2 * h5_m4 - e_4 * fs_10_143_77 * r_2 * h7_m4 - e_4 * fs_8_273_105 * r_4 * h5_m4;

        pc_44[k] = e_0 * fs_165_16_30 * h1_m1 + e_1 * fs_117_4_3 * h3_m3 - e_1 * fs_39_2_5 * h3_m1 - e_1 * fs_33_2_30 * r_2 * h1_m1 - e_2 * fs_45_14_21 * h5_m3 + e_2 * fs_215_28_2 * h5_m1 - e_2 * fs_39_2_3 * r_2 * h3_m3 + e_2 * fs_13_1_5 * r_2 * h3_m1 + e_2 * fs_99_14_30 * r_4 * h1_m1 + e_3 * fs_85_143_70 * h7_m3 + e_3 * fs_85_429_210 * h7_m1 + e_3 * fs_90_91_21 * r_2 * h5_m3 - e_3 * fs_215_91_2 * r_2 * h5_m1 + e_3 * fs_39_11_3 * r_4 * h3_m3 - e_3 * fs_26_11_5 * r_4 * h3_m1 - e_3 * fs_22_21_30 * r_6 * h1_m1 - e_4 * fs_10_143_70 * r_2 * h7_m3 - e_4 * fs_10_429_210 * r_2 * h7_m1 - e_4 * fs_6_91_21 * r_4 * h5_m3 + e_4 * fs_43_273_2 * r_4 * h5_m1 - e_4 * fs_2_11_3 * r_6 * h3_m3 + e_4 * fs_4_33_5 * r_6 * h3_m1 + e_4 * fs_1_21_30 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p2, ph3_p3, ph5_p1, ph5_p2, ph5_p3, ph7_p1, ph7_p2, ph7_p3, ab_2, pc_45, pc_46 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];

        pc_45[k] = - e_1 * fs_39_4_3 * h3_p2 - e_2 * fs_5_7_21 * h5_p2 + e_2 * fs_13_2_3 * r_2 * h3_p2 + e_3 * fs_170_429_210 * h7_p2 + e_3 * fs_20_91_21 * r_2 * h5_p2 - e_3 * fs_13_11_3 * r_4 * h3_p2 - e_4 * fs_20_429_210 * r_2 * h7_p2 - e_4 * fs_4_273_21 * r_4 * h5_p2 + e_4 * fs_2_33_3 * r_6 * h3_p2;

        pc_46[k] = - e_0 * fs_165_16_30 * h1_p1 + e_1 * fs_39_2_5 * h3_p1 + e_1 * fs_117_4_3 * h3_p3 + e_1 * fs_33_2_30 * r_2 * h1_p1 - e_2 * fs_215_28_2 * h5_p1 - e_2 * fs_45_14_21 * h5_p3 - e_2 * fs_13_1_5 * r_2 * h3_p1 - e_2 * fs_39_2_3 * r_2 * h3_p3 - e_2 * fs_99_14_30 * r_4 * h1_p1 - e_3 * fs_85_429_210 * h7_p1 + e_3 * fs_85_143_70 * h7_p3 + e_3 * fs_215_91_2 * r_2 * h5_p1 + e_3 * fs_90_91_21 * r_2 * h5_p3 + e_3 * fs_26_11_5 * r_4 * h3_p1 + e_3 * fs_39_11_3 * r_4 * h3_p3 + e_3 * fs_22_21_30 * r_6 * h1_p1 + e_4 * fs_10_429_210 * r_2 * h7_p1 - e_4 * fs_10_143_70 * r_2 * h7_p3 - e_4 * fs_43_273_2 * r_4 * h5_p1 - e_4 * fs_6_91_21 * r_4 * h5_p3 - e_4 * fs_4_33_5 * r_6 * h3_p1 - e_4 * fs_2_11_3 * r_6 * h3_p3 - e_4 * fs_1_21_30 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p4, ph5_p5, ph7_0, ph7_p1, ph7_p4, ph7_p5, ab_2, pc_47, pc_48 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_47[k] = - e_0 * fs_165_4_3 * h1_0 - e_1 * fs_39_4_3 * h3_0 + e_1 * fs_66_1_3 * r_2 * h1_0 + e_2 * fs_100_7_3 * h5_0 - e_2 * fs_10_7_105 * h5_p4 + e_2 * fs_13_2_3 * r_2 * h3_0 - e_2 * fs_198_7_3 * r_4 * h1_0 + e_3 * fs_595_429_3 * h7_0 + e_3 * fs_85_143_77 * h7_p4 - e_3 * fs_400_91_3 * r_2 * h5_0 + e_3 * fs_40_91_105 * r_2 * h5_p4 - e_3 * fs_13_11_3 * r_4 * h3_0 + e_3 * fs_88_21_3 * r_6 * h1_0 - e_4 * fs_70_429_3 * r_2 * h7_0 - e_4 * fs_10_143_77 * r_2 * h7_p4 + e_4 * fs_80_273_3 * r_4 * h5_0 - e_4 * fs_8_273_105 * r_4 * h5_p4 + e_4 * fs_2_33_3 * r_6 * h3_0 - e_4 * fs_4_21_3 * r_8 * h1_0;

        pc_48[k] = e_0 * fs_165_16_2 * h1_p1 + e_1 * fs_117_4_3 * h3_p1 - e_1 * fs_33_2_2 * r_2 * h1_p1 + e_2 * fs_75_28_30 * h5_p1 + e_2 * fs_75_14_7 * h5_p5 - e_2 * fs_39_2_3 * r_2 * h3_p1 + e_2 * fs_99_14_2 * r_4 * h1_p1 + e_3 * fs_85_429_14 * h7_p1 + e_3 * fs_85_429_462 * h7_p5 - e_3 * fs_75_91_30 * r_2 * h5_p1 - e_3 * fs_150_91_7 * r_2 * h5_p5 + e_3 * fs_39_11_3 * r_4 * h3_p1 - e_3 * fs_22_21_2 * r_6 * h1_p1 - e_4 * fs_10_429_14 * r_2 * h7_p1 - e_4 * fs_10_429_462 * r_2 * h7_p5 + e_4 * fs_5_91_30 * r_4 * h5_p1 + e_4 * fs_10_91_7 * r_4 * h5_p5 - e_4 * fs_2_11_3 * r_6 * h3_p1 + e_4 * fs_1_21_2 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m2, ph3_m1, ph5_m5, ph5_m4, ph5_m2, ph5_m1, ph7_m6, ph7_m5, ph7_m4, ph7_m2, ph7_m1, ab_2, pc_49, pc_50, pc_51 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_49[k] = e_3 * fs_85_429_858 * h7_m6 - e_4 * fs_10_429_858 * r_2 * h7_m6;

        pc_50[k] = e_0 * fs_165_16_42 * h1_m1 - e_1 * fs_39_4_7 * h3_m1 - e_1 * fs_33_2_42 * r_2 * h1_m1 - e_2 * fs_25_2_3 * h5_m5 - e_2 * fs_55_28_70 * h5_m1 + e_2 * fs_13_2_7 * r_2 * h3_m1 + e_2 * fs_99_14_42 * r_4 * h1_m1 + e_3 * fs_170_143_22 * h7_m5 - e_3 * fs_170_429_6 * h7_m1 + e_3 * fs_50_13_3 * r_2 * h5_m5 + e_3 * fs_55_91_70 * r_2 * h5_m1 - e_3 * fs_13_11_7 * r_4 * h3_m1 - e_3 * fs_22_21_42 * r_6 * h1_m1 - e_4 * fs_20_143_22 * r_2 * h7_m5 + e_4 * fs_20_429_6 * r_2 * h7_m1 - e_4 * fs_10_39_3 * r_4 * h5_m5 - e_4 * fs_11_273_70 * r_4 * h5_m1 + e_4 * fs_2_33_7 * r_6 * h3_m1 + e_4 * fs_1_21_42 * r_8 * h1_m1;

        pc_51[k] = - e_1 * fs_39_4_7 * h3_m2 - e_2 * fs_5_2_3 * h5_m4 + e_2 * f_35_2 * h5_m2 + e_2 * fs_13_2_7 * r_2 * h3_m2 + e_3 * fs_85_143_55 * h7_m4 + e_3 * fs_85_143_10 * h7_m2 + e_3 * fs_10_13_3 * r_2 * h5_m4 - e_3 * f_70_13 * r_2 * h5_m2 - e_3 * fs_13_11_7 * r_4 * h3_m2 - e_4 * fs_10_143_55 * r_2 * h7_m4 - e_4 * fs_10_143_10 * r_2 * h7_m2 - e_4 * fs_2_39_3 * r_4 * h5_m4 + e_4 * f_14_39 * r_4 * h5_m2 + e_4 * fs_2_33_7 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p4, ph7_p2, ph7_p3, ph7_p4, ab_2, pc_52, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_52[k] = - e_1 * fs_117_4_7 * h3_p3 + e_2 * f_15_1 * h5_p3 + e_2 * fs_39_2_7 * r_2 * h3_p3 + e_3 * fs_340_429_30 * h7_p3 - e_3 * f_60_13 * r_2 * h5_p3 - e_3 * fs_39_11_7 * r_4 * h3_p3 - e_4 * fs_40_429_30 * r_2 * h7_p3 + e_4 * f_4_13 * r_4 * h5_p3 + e_4 * fs_2_11_7 * r_6 * h3_p3;

        pc_53[k] = e_1 * fs_39_4_7 * h3_p2 - e_2 * f_35_2 * h5_p2 - e_2 * fs_5_2_3 * h5_p4 - e_2 * fs_13_2_7 * r_2 * h3_p2 - e_3 * fs_85_143_10 * h7_p2 + e_3 * fs_85_143_55 * h7_p4 + e_3 * f_70_13 * r_2 * h5_p2 + e_3 * fs_10_13_3 * r_2 * h5_p4 + e_3 * fs_13_11_7 * r_4 * h3_p2 + e_4 * fs_10_143_10 * r_2 * h7_p2 - e_4 * fs_10_143_55 * r_2 * h7_p4 - e_4 * f_14_39 * r_4 * h5_p2 - e_4 * fs_2_39_3 * r_4 * h5_p4 - e_4 * fs_2_33_7 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p5, ph7_0, ph7_p1, ph7_p5, ph7_p6, ab_2, pc_54, pc_55 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_54[k] = - e_0 * fs_165_16_42 * h1_p1 + e_1 * fs_39_4_7 * h3_p1 + e_1 * fs_33_2_42 * r_2 * h1_p1 + e_2 * fs_55_28_70 * h5_p1 - e_2 * fs_25_2_3 * h5_p5 - e_2 * fs_13_2_7 * r_2 * h3_p1 - e_2 * fs_99_14_42 * r_4 * h1_p1 + e_3 * fs_170_429_6 * h7_p1 + e_3 * fs_170_143_22 * h7_p5 - e_3 * fs_55_91_70 * r_2 * h5_p1 + e_3 * fs_50_13_3 * r_2 * h5_p5 + e_3 * fs_13_11_7 * r_4 * h3_p1 + e_3 * fs_22_21_42 * r_6 * h1_p1 - e_4 * fs_20_429_6 * r_2 * h7_p1 - e_4 * fs_20_143_22 * r_2 * h7_p5 + e_4 * fs_11_273_70 * r_4 * h5_p1 - e_4 * fs_10_39_3 * r_4 * h5_p5 - e_4 * fs_2_33_7 * r_6 * h3_p1 - e_4 * fs_1_21_42 * r_8 * h1_p1;

        pc_55[k] = - e_0 * fs_165_8_7 * h1_0 - e_1 * fs_117_4_7 * h3_0 + e_1 * fs_33_1_7 * r_2 * h1_0 - e_2 * fs_75_14_7 * h5_0 + e_2 * fs_39_2_7 * r_2 * h3_0 - e_2 * fs_99_7_7 * r_4 * h1_0 - e_3 * fs_85_429_7 * h7_0 + e_3 * fs_85_429_858 * h7_p6 + e_3 * fs_150_91_7 * r_2 * h5_0 - e_3 * fs_39_11_7 * r_4 * h3_0 + e_3 * fs_44_21_7 * r_6 * h1_0 + e_4 * fs_10_429_7 * r_2 * h7_0 - e_4 * fs_10_429_858 * r_2 * h7_p6 - e_4 * fs_10_91_7 * r_4 * h5_0 + e_4 * fs_2_11_7 * r_6 * h3_0 - e_4 * fs_2_21_7 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m7, ph7_m6, ph7_m2, ph7_m1, ab_2, pc_56, pc_57 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_56[k] = e_0 * fs_165_8_14 * h1_m1 + e_1 * fs_39_4_21 * h3_m1 - e_1 * fs_33_1_14 * r_2 * h1_m1 + e_2 * fs_5_14_210 * h5_m1 - e_2 * fs_13_2_21 * r_2 * h3_m1 + e_2 * fs_99_7_14 * r_4 * h1_m1 + e_3 * fs_85_858_6006 * h7_m7 + e_3 * fs_85_858_2 * h7_m1 - e_3 * fs_10_91_210 * r_2 * h5_m1 + e_3 * fs_13_11_21 * r_4 * h3_m1 - e_3 * fs_44_21_14 * r_6 * h1_m1 - e_4 * fs_5_429_6006 * r_2 * h7_m7 - e_4 * fs_5_429_2 * r_2 * h7_m1 + e_4 * fs_2_273_210 * r_4 * h5_m1 - e_4 * fs_2_33_21 * r_6 * h3_m1 + e_4 * fs_2_21_14 * r_8 * h1_m1;

        pc_57[k] = - e_1 * fs_39_4_35 * h3_m2 - e_2 * fs_5_1_5 * h5_m2 + e_2 * fs_13_2_35 * r_2 * h3_m2 + e_3 * fs_85_286_286 * h7_m6 - e_3 * fs_85_286_2 * h7_m2 + e_3 * fs_20_13_5 * r_2 * h5_m2 - e_3 * fs_13_11_35 * r_4 * h3_m2 - e_4 * fs_5_143_286 * r_2 * h7_m6 + e_4 * fs_5_143_2 * r_2 * h7_m2 - e_4 * fs_4_39_5 * r_4 * h5_m2 + e_4 * fs_2_33_35 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_p3, ph5_m5, ph5_m3, ph5_p3, ph5_p4, ph5_p5, ph7_m5, ph7_m3, ph7_p3, ph7_p4, ph7_p5, ab_2, pc_58, pc_59, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_58[k] = e_1 * fs_39_4_21 * h3_m3 + e_2 * fs_5_1_15 * h5_m5 + e_2 * fs_10_1_3 * h5_m3 - e_2 * fs_13_2_21 * r_2 * h3_m3 + e_3 * fs_85_286_110 * h7_m5 + e_3 * fs_85_286_10 * h7_m3 - e_3 * fs_20_13_15 * r_2 * h5_m5 - e_3 * fs_40_13_3 * r_2 * h5_m3 + e_3 * fs_13_11_21 * r_4 * h3_m3 - e_4 * fs_5_143_110 * r_2 * h7_m5 - e_4 * fs_5_143_10 * r_2 * h7_m3 + e_4 * fs_4_39_15 * r_4 * h5_m5 + e_4 * fs_8_39_3 * r_4 * h5_m3 - e_4 * fs_2_33_21 * r_6 * h3_m3;

        pc_59[k] = e_2 * f_30_1 * h5_p4 + e_3 * fs_85_429_165 * h7_p4 - e_3 * f_120_13 * r_2 * h5_p4 - e_4 * fs_10_429_165 * r_2 * h7_p4 + e_4 * f_8_13 * r_4 * h5_p4;

        pc_60[k] = - e_1 * fs_39_4_21 * h3_p3 - e_2 * fs_10_1_3 * h5_p3 + e_2 * fs_5_1_15 * h5_p5 + e_2 * fs_13_2_21 * r_2 * h3_p3 - e_3 * fs_85_286_10 * h7_p3 + e_3 * fs_85_286_110 * h7_p5 + e_3 * fs_40_13_3 * r_2 * h5_p3 - e_3 * fs_20_13_15 * r_2 * h5_p5 - e_3 * fs_13_11_21 * r_4 * h3_p3 + e_4 * fs_5_143_10 * r_2 * h7_p3 - e_4 * fs_5_143_110 * r_2 * h7_p5 - e_4 * fs_8_39_3 * r_4 * h5_p3 + e_4 * fs_4_39_15 * r_4 * h5_p5 + e_4 * fs_2_33_21 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 26 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph7_p6, ph7_p7, ab_2, pc_61, pc_62 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_61[k] = e_1 * fs_39_4_35 * h3_p2 + e_2 * fs_5_1_5 * h5_p2 - e_2 * fs_13_2_35 * r_2 * h3_p2 + e_3 * fs_85_286_2 * h7_p2 + e_3 * fs_85_286_286 * h7_p6 - e_3 * fs_20_13_5 * r_2 * h5_p2 + e_3 * fs_13_11_35 * r_4 * h3_p2 - e_4 * fs_5_143_2 * r_2 * h7_p2 - e_4 * fs_5_143_286 * r_2 * h7_p6 + e_4 * fs_4_39_5 * r_4 * h5_p2 - e_4 * fs_2_33_35 * r_6 * h3_p2;

        pc_62[k] = - e_0 * fs_165_8_14 * h1_p1 - e_1 * fs_39_4_21 * h3_p1 + e_1 * fs_33_1_14 * r_2 * h1_p1 - e_2 * fs_5_14_210 * h5_p1 + e_2 * fs_13_2_21 * r_2 * h3_p1 - e_2 * fs_99_7_14 * r_4 * h1_p1 - e_3 * fs_85_858_2 * h7_p1 + e_3 * fs_85_858_6006 * h7_p7 + e_3 * fs_10_91_210 * r_2 * h5_p1 - e_3 * fs_13_11_21 * r_4 * h3_p1 + e_3 * fs_44_21_14 * r_6 * h1_p1 + e_4 * fs_5_429_2 * r_2 * h7_p1 - e_4 * fs_5_429_6006 * r_2 * h7_p7 - e_4 * fs_2_273_210 * r_4 * h5_p1 + e_4 * fs_2_33_21 * r_6 * h3_p1 - e_4 * fs_2_21_14 * r_8 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[63] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62};

    for (size_t n = 0; n < 63; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
