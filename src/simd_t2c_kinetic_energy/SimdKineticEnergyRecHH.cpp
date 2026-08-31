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



#include "SimdKineticEnergyRecHH.hpp"

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
compute_hh_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hh_kinetic_energy: Basis functions must be of angular momenta five and five"));
    }

    if (harmonics.size() < 10)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hh_kinetic_energy: Harmonics must reach angular momentum 10"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hh_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 121 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time.

    auto buffer = CSimdMatrix(7, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);
    auto *pe_6 = buffer.data(6);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);
    std::fill(pe_4, pe_4 + nmax, 0.0);
    std::fill(pe_5, pe_5 + nmax, 0.0);
    std::fill(pe_6, pe_6 + nmax, 0.0);

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

            const auto ff_0 = fbase * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_6 = fbase * fmu * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fterm = std::exp(-fmu * ab_2[k]);

                pe_0[k] += ff_0 * fterm;
                pe_1[k] += ff_1 * fterm;
                pe_2[k] += ff_2 * fterm;
                pe_3[k] += ff_3 * fterm;
                pe_4[k] += ff_4 * fterm;
                pe_5[k] += ff_5 * fterm;
                pe_6[k] += ff_6 * fterm;
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
    const auto *ph10_m10 = harmonics[9].data(0);
    const auto *ph10_m9 = harmonics[9].data(1);
    const auto *ph10_m8 = harmonics[9].data(2);
    const auto *ph10_m7 = harmonics[9].data(3);
    const auto *ph10_m6 = harmonics[9].data(4);
    const auto *ph10_m5 = harmonics[9].data(5);
    const auto *ph10_m4 = harmonics[9].data(6);
    const auto *ph10_m3 = harmonics[9].data(7);
    const auto *ph10_m2 = harmonics[9].data(8);
    const auto *ph10_m1 = harmonics[9].data(9);
    const auto *ph10_0 = harmonics[9].data(10);
    const auto *ph10_p1 = harmonics[9].data(11);
    const auto *ph10_p2 = harmonics[9].data(12);
    const auto *ph10_p3 = harmonics[9].data(13);
    const auto *ph10_p4 = harmonics[9].data(14);
    const auto *ph10_p5 = harmonics[9].data(15);
    const auto *ph10_p6 = harmonics[9].data(16);
    const auto *ph10_p7 = harmonics[9].data(17);
    const auto *ph10_p8 = harmonics[9].data(18);
    const auto *ph10_p9 = harmonics[9].data(19);
    const auto *ph10_p10 = harmonics[9].data(20);

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
    auto *pc_11 = values + 12 * nvalues;
    auto *pc_12 = values + 13 * nvalues;
    auto *pc_13 = values + 14 * nvalues;
    auto *pc_14 = values + 15 * nvalues;
    auto *pc_15 = values + 16 * nvalues;
    auto *pc_16 = values + 17 * nvalues;
    auto *pc_17 = values + 18 * nvalues;
    auto *pc_18 = values + 19 * nvalues;
    auto *pc_19 = values + 20 * nvalues;
    auto *pc_20 = values + 21 * nvalues;
    auto *pc_21 = values + 24 * nvalues;
    auto *pc_22 = values + 25 * nvalues;
    auto *pc_23 = values + 26 * nvalues;
    auto *pc_24 = values + 27 * nvalues;
    auto *pc_25 = values + 28 * nvalues;
    auto *pc_26 = values + 29 * nvalues;
    auto *pc_27 = values + 30 * nvalues;
    auto *pc_28 = values + 31 * nvalues;
    auto *pc_29 = values + 32 * nvalues;
    auto *pc_30 = values + 36 * nvalues;
    auto *pc_31 = values + 37 * nvalues;
    auto *pc_32 = values + 38 * nvalues;
    auto *pc_33 = values + 39 * nvalues;
    auto *pc_34 = values + 40 * nvalues;
    auto *pc_35 = values + 41 * nvalues;
    auto *pc_36 = values + 42 * nvalues;
    auto *pc_37 = values + 43 * nvalues;
    auto *pc_38 = values + 48 * nvalues;
    auto *pc_39 = values + 49 * nvalues;
    auto *pc_40 = values + 50 * nvalues;
    auto *pc_41 = values + 51 * nvalues;
    auto *pc_42 = values + 52 * nvalues;
    auto *pc_43 = values + 53 * nvalues;
    auto *pc_44 = values + 54 * nvalues;
    auto *pc_45 = values + 60 * nvalues;
    auto *pc_46 = values + 61 * nvalues;
    auto *pc_47 = values + 62 * nvalues;
    auto *pc_48 = values + 63 * nvalues;
    auto *pc_49 = values + 64 * nvalues;
    auto *pc_50 = values + 65 * nvalues;
    auto *pc_51 = values + 72 * nvalues;
    auto *pc_52 = values + 73 * nvalues;
    auto *pc_53 = values + 74 * nvalues;
    auto *pc_54 = values + 75 * nvalues;
    auto *pc_55 = values + 76 * nvalues;
    auto *pc_56 = values + 84 * nvalues;
    auto *pc_57 = values + 85 * nvalues;
    auto *pc_58 = values + 86 * nvalues;
    auto *pc_59 = values + 87 * nvalues;
    auto *pc_60 = values + 96 * nvalues;
    auto *pc_61 = values + 97 * nvalues;
    auto *pc_62 = values + 98 * nvalues;
    auto *pc_63 = values + 108 * nvalues;
    auto *pc_64 = values + 109 * nvalues;
    auto *pc_65 = values + 120 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_100_429 = 100.0 / 429.0;
    const auto f_10125_8 = 10125.0 / 8.0;
    const auto f_1022_2717 = 1022.0 / 2717.0;
    const auto f_1029_143 = 1029.0 / 143.0;
    const auto f_10731_286 = 10731.0 / 286.0;
    const auto f_10_429 = 10.0 / 429.0;
    const auto f_1102_187 = 1102.0 / 187.0;
    const auto f_1125_11 = 1125.0 / 11.0;
    const auto f_1125_143 = 1125.0 / 143.0;
    const auto f_1125_2 = 1125.0 / 2.0;
    const auto f_1125_8 = 1125.0 / 8.0;
    const auto f_116_561 = 116.0 / 561.0;
    const auto f_1224_143 = 1224.0 / 143.0;
    const auto f_12285_32 = 12285.0 / 32.0;
    const auto f_12285_8 = 12285.0 / 8.0;
    const auto f_1250_11 = 1250.0 / 11.0;
    const auto f_1250_143 = 1250.0 / 143.0;
    const auto f_125_11 = 125.0 / 11.0;
    const auto f_125_143 = 125.0 / 143.0;
    const auto f_125_2 = 125.0 / 2.0;
    const auto f_1260_46189 = 1260.0 / 46189.0;
    const auto f_126_46189 = 126.0 / 46189.0;
    const auto f_1368_187 = 1368.0 / 187.0;
    const auto f_13770_143 = 13770.0 / 143.0;
    const auto f_14175_16 = 14175.0 / 16.0;
    const auto f_1425_22 = 1425.0 / 22.0;
    const auto f_14490_46189 = 14490.0 / 46189.0;
    const auto f_1449_46189 = 1449.0 / 46189.0;
    const auto f_1470_2717 = 1470.0 / 2717.0;
    const auto f_15120_46189 = 15120.0 / 46189.0;
    const auto f_1520_187 = 1520.0 / 187.0;
    const auto f_1575_16 = 1575.0 / 16.0;
    const auto f_160_561 = 160.0 / 561.0;
    const auto f_16875_8 = 16875.0 / 8.0;
    const auto f_16_187 = 16.0 / 187.0;
    const auto f_1710_11 = 1710.0 / 11.0;
    const auto f_173880_46189 = 173880.0 / 46189.0;
    const auto f_1824_187 = 1824.0 / 187.0;
    const auto f_1875_11 = 1875.0 / 11.0;
    const auto f_1875_143 = 1875.0 / 143.0;
    const auto f_1875_2 = 1875.0 / 2.0;
    const auto f_1900_11 = 1900.0 / 11.0;
    const auto f_195_2 = 195.0 / 2.0;
    const auto f_196_2717 = 196.0 / 2717.0;
    const auto f_19992_2717 = 19992.0 / 2717.0;
    const auto f_204_143 = 204.0 / 143.0;
    const auto f_20580_2717 = 20580.0 / 2717.0;
    const auto f_20_143 = 20.0 / 143.0;
    const auto f_20_187 = 20.0 / 187.0;
    const auto f_21462_2717 = 21462.0 / 2717.0;
    const auto f_2280_11 = 2280.0 / 11.0;
    const auto f_228_11 = 228.0 / 11.0;
    const auto f_2295_143 = 2295.0 / 143.0;
    const auto f_2295_4 = 2295.0 / 4.0;
    const auto f_23625_16 = 23625.0 / 16.0;
    const auto f_24_143 = 24.0 / 143.0;
    const auto f_26460_46189 = 26460.0 / 46189.0;
    const auto f_2755_22 = 2755.0 / 22.0;
    const auto f_285_11 = 285.0 / 11.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_304290_46189 = 304290.0 / 46189.0;
    const auto f_3060_11 = 3060.0 / 11.0;
    const auto f_30_143 = 30.0 / 143.0;
    const auto f_31752_46189 = 31752.0 / 46189.0;
    const auto f_3375_4 = 3375.0 / 4.0;
    const auto f_365148_46189 = 365148.0 / 46189.0;
    const auto f_36_143 = 36.0 / 143.0;
    const auto f_375_1 = 375.0;
    const auto f_4116_2717 = 4116.0 / 2717.0;
    const auto f_434_2717 = 434.0 / 2717.0;
    const auto f_4557_286 = 4557.0 / 286.0;
    const auto f_456_187 = 456.0 / 187.0;
    const auto f_4590_11 = 4590.0 / 11.0;
    const auto f_4725_8 = 4725.0 / 8.0;
    const auto f_48_187 = 48.0 / 187.0;
    const auto f_4998_143 = 4998.0 / 143.0;
    const auto f_50_143 = 50.0 / 143.0;
    const auto f_5145_143 = 5145.0 / 143.0;
    const auto f_551_11 = 551.0 / 11.0;
    const auto f_5625_4 = 5625.0 / 4.0;
    const auto f_5670_46189 = 5670.0 / 46189.0;
    const auto f_570_11 = 570.0 / 11.0;
    const auto f_570_187 = 570.0 / 187.0;
    const auto f_585_1 = 585.0;
    const auto f_625_1 = 625.0;
    const auto f_64_187 = 64.0 / 187.0;
    const auto f_65205_46189 = 65205.0 / 46189.0;
    const auto f_684_11 = 684.0 / 11.0;
    const auto f_6_143 = 6.0 / 143.0;
    const auto f_70_2717 = 70.0 / 2717.0;
    const auto f_735_286 = 735.0 / 286.0;
    const auto f_750_11 = 750.0 / 11.0;
    const auto f_750_143 = 750.0 / 143.0;
    const auto f_760_11 = 760.0 / 11.0;
    const auto f_765_11 = 765.0 / 11.0;
    const auto f_765_2 = 765.0 / 2.0;
    const auto f_765_8 = 765.0 / 8.0;
    const auto f_7875_8 = 7875.0 / 8.0;
    const auto f_78_11 = 78.0 / 11.0;
    const auto f_816_143 = 816.0 / 143.0;
    const auto f_9114_2717 = 9114.0 / 2717.0;
    const auto f_912_11 = 912.0 / 11.0;
    const auto f_9180_143 = 9180.0 / 143.0;
    const auto f_952_2717 = 952.0 / 2717.0;
    const auto f_980_2717 = 980.0 / 2717.0;
    const auto fs_100_561_3 = std::sqrt(10000.0 / 104907.0);
    const auto fs_1029_143_15 = std::sqrt(15882615.0 / 20449.0);
    const auto fs_1029_2717_210 = std::sqrt(222356610.0 / 7382089.0);
    const auto fs_1029_572_210 = std::sqrt(111178305.0 / 163592.0);
    const auto fs_102_143_10 = std::sqrt(104040.0 / 20449.0);
    const auto fs_102_143_70 = std::sqrt(728280.0 / 20449.0);
    const auto fs_10_187_30 = std::sqrt(3000.0 / 34969.0);
    const auto fs_10_429_15 = std::sqrt(500.0 / 61347.0);
    const auto fs_10_429_21 = std::sqrt(700.0 / 61347.0);
    const auto fs_10_429_5 = std::sqrt(500.0 / 184041.0);
    const auto fs_11025_32_6 = std::sqrt(364651875.0 / 512.0);
    const auto fs_1125_4_14 = std::sqrt(8859375.0 / 8.0);
    const auto fs_1125_4_35 = std::sqrt(44296875.0 / 16.0);
    const auto fs_1125_8_15 = std::sqrt(18984375.0 / 64.0);
    const auto fs_1125_8_21 = std::sqrt(26578125.0 / 64.0);
    const auto fs_1125_8_5 = std::sqrt(6328125.0 / 64.0);
    const auto fs_112_2717_55 = std::sqrt(62720.0 / 671099.0);
    const auto fs_114_187_154 = std::sqrt(181944.0 / 3179.0);
    const auto fs_114_187_30 = std::sqrt(389880.0 / 34969.0);
    const auto fs_114_187_55 = std::sqrt(64980.0 / 3179.0);
    const auto fs_114_187_70 = std::sqrt(909720.0 / 34969.0);
    const auto fs_1235_44_6 = std::sqrt(4575675.0 / 968.0);
    const auto fs_125_11_15 = std::sqrt(234375.0 / 121.0);
    const auto fs_125_11_21 = std::sqrt(328125.0 / 121.0);
    const auto fs_125_11_5 = std::sqrt(78125.0 / 121.0);
    const auto fs_125_143_15 = std::sqrt(234375.0 / 20449.0);
    const auto fs_125_143_21 = std::sqrt(328125.0 / 20449.0);
    const auto fs_125_143_5 = std::sqrt(78125.0 / 20449.0);
    const auto fs_125_1_14 = std::sqrt(218750.0);
    const auto fs_125_1_35 = std::sqrt(546875.0);
    const auto fs_125_2_15 = std::sqrt(234375.0 / 4.0);
    const auto fs_125_2_21 = std::sqrt(328125.0 / 4.0);
    const auto fs_125_2_5 = std::sqrt(78125.0 / 4.0);
    const auto fs_126_2717_7 = std::sqrt(111132.0 / 7382089.0);
    const auto fs_126_46189_143 = std::sqrt(15876.0 / 14919047.0);
    const auto fs_126_46189_3003 = std::sqrt(333396.0 / 14919047.0);
    const auto fs_126_46189_30030 = std::sqrt(3333960.0 / 14919047.0);
    const auto fs_126_46189_33 = std::sqrt(47628.0 / 193947611.0);
    const auto fs_126_46189_36465 = std::sqrt(238140.0 / 877591.0);
    const auto fs_126_46189_46189 = std::sqrt(15876.0 / 46189.0);
    const auto fs_126_46189_92378 = std::sqrt(31752.0 / 46189.0);
    const auto fs_12_143_5 = std::sqrt(720.0 / 20449.0);
    const auto fs_1323_286_7 = std::sqrt(12252303.0 / 81796.0);
    const auto fs_1425_22_11 = std::sqrt(2030625.0 / 44.0);
    const auto fs_1425_22_7 = std::sqrt(14214375.0 / 484.0);
    const auto fs_1425_44_30 = std::sqrt(30459375.0 / 968.0);
    const auto fs_1449_46189_143 = std::sqrt(2099601.0 / 14919047.0);
    const auto fs_1449_46189_3003 = std::sqrt(44091621.0 / 14919047.0);
    const auto fs_1449_46189_30030 = std::sqrt(440916210.0 / 14919047.0);
    const auto fs_1449_46189_33 = std::sqrt(6298803.0 / 193947611.0);
    const auto fs_1449_46189_36465 = std::sqrt(31494015.0 / 877591.0);
    const auto fs_1449_46189_46189 = std::sqrt(2099601.0 / 46189.0);
    const auto fs_1449_46189_92378 = std::sqrt(4199202.0 / 46189.0);
    const auto fs_1449_92378_10010 = std::sqrt(73486035.0 / 29838094.0);
    const auto fs_1449_92378_2002 = std::sqrt(14697207.0 / 29838094.0);
    const auto fs_1449_92378_22 = std::sqrt(2099601.0 / 387895222.0);
    const auto fs_1470_2717_14 = std::sqrt(30252600.0 / 7382089.0);
    const auto fs_1470_2717_143 = std::sqrt(2160900.0 / 51623.0);
    const auto fs_1470_2717_286 = std::sqrt(4321800.0 / 51623.0);
    const auto fs_1470_2717_6 = std::sqrt(12965400.0 / 7382089.0);
    const auto fs_1470_2717_66 = std::sqrt(12965400.0 / 671099.0);
    const auto fs_147_11_2 = std::sqrt(43218.0 / 121.0);
    const auto fs_147_143_429 = std::sqrt(64827.0 / 143.0);
    const auto fs_147_143_858 = std::sqrt(129654.0 / 143.0);
    const auto fs_147_22_15 = std::sqrt(324135.0 / 484.0);
    const auto fs_147_26_55 = std::sqrt(1188495.0 / 676.0);
    const auto fs_147_2717_286 = std::sqrt(43218.0 / 51623.0);
    const auto fs_147_2717_6006 = std::sqrt(907578.0 / 51623.0);
    const auto fs_147_2717_66 = std::sqrt(129654.0 / 671099.0);
    const auto fs_147_286_1430 = std::sqrt(108045.0 / 286.0);
    const auto fs_147_286_5005 = std::sqrt(756315.0 / 572.0);
    const auto fs_147_572_286 = std::sqrt(21609.0 / 1144.0);
    const auto fs_147_572_6006 = std::sqrt(453789.0 / 1144.0);
    const auto fs_147_572_66 = std::sqrt(64827.0 / 14872.0);
    const auto fs_14_209_15 = std::sqrt(2940.0 / 43681.0);
    const auto fs_14_247_55 = std::sqrt(10780.0 / 61009.0);
    const auto fs_14_2717_1430 = std::sqrt(1960.0 / 51623.0);
    const auto fs_14_2717_5005 = std::sqrt(6860.0 / 51623.0);
    const auto fs_152_187_35 = std::sqrt(808640.0 / 34969.0);
    const auto fs_1530_11_5 = std::sqrt(11704500.0 / 121.0);
    const auto fs_1575_16_15 = std::sqrt(37209375.0 / 256.0);
    const auto fs_1575_16_21 = std::sqrt(52093125.0 / 256.0);
    const auto fs_1575_16_5 = std::sqrt(12403125.0 / 256.0);
    const auto fs_1575_8_14 = std::sqrt(17364375.0 / 32.0);
    const auto fs_1575_8_35 = std::sqrt(86821875.0 / 64.0);
    const auto fs_16_187_5 = std::sqrt(1280.0 / 34969.0);
    const auto fs_16_187_7 = std::sqrt(1792.0 / 34969.0);
    const auto fs_16_561_35 = std::sqrt(8960.0 / 314721.0);
    const auto fs_171_187_66 = std::sqrt(175446.0 / 3179.0);
    const auto fs_171_22_66 = std::sqrt(87723.0 / 22.0);
    const auto fs_1764_2717_70 = std::sqrt(217818720.0 / 7382089.0);
    const auto fs_1764_2717_77 = std::sqrt(21781872.0 / 671099.0);
    const auto fs_1764_46189_165 = std::sqrt(46675440.0 / 193947611.0);
    const auto fs_1890_46189_143 = std::sqrt(3572100.0 / 14919047.0);
    const auto fs_189_46189_10010 = std::sqrt(2500470.0 / 14919047.0);
    const auto fs_189_46189_110 = std::sqrt(357210.0 / 193947611.0);
    const auto fs_190_11_2 = std::sqrt(72200.0 / 121.0);
    const auto fs_190_11_35 = std::sqrt(1263500.0 / 121.0);
    const auto fs_190_187_42 = std::sqrt(1516200.0 / 34969.0);
    const auto fs_196_2717_15 = std::sqrt(576240.0 / 7382089.0);
    const auto fs_19_11_14 = std::sqrt(5054.0 / 121.0);
    const auto fs_19_11_231 = std::sqrt(7581.0 / 11.0);
    const auto fs_19_17_42 = std::sqrt(15162.0 / 289.0);
    const auto fs_19_1_5 = std::sqrt(1805.0);
    const auto fs_19_2_42 = std::sqrt(7581.0 / 2.0);
    const auto fs_20286_46189_165 = std::sqrt(6172826940.0 / 193947611.0);
    const auto fs_204_143_15 = std::sqrt(624240.0 / 20449.0);
    const auto fs_204_143_21 = std::sqrt(873936.0 / 20449.0);
    const auto fs_204_143_30 = std::sqrt(1248480.0 / 20449.0);
    const auto fs_204_143_35 = std::sqrt(1456560.0 / 20449.0);
    const auto fs_204_143_6 = std::sqrt(249696.0 / 20449.0);
    const auto fs_20_187_11 = std::sqrt(400.0 / 3179.0);
    const auto fs_20_187_7 = std::sqrt(2800.0 / 34969.0);
    const auto fs_20_429_14 = std::sqrt(5600.0 / 184041.0);
    const auto fs_20_429_35 = std::sqrt(14000.0 / 184041.0);
    const auto fs_20_561_42 = std::sqrt(5600.0 / 104907.0);
    const auto fs_21735_46189_143 = std::sqrt(472410225.0 / 14919047.0);
    const auto fs_228_11_5 = std::sqrt(259920.0 / 121.0);
    const auto fs_228_11_7 = std::sqrt(363888.0 / 121.0);
    const auto fs_2295_143_15 = std::sqrt(79005375.0 / 20449.0);
    const auto fs_2295_143_21 = std::sqrt(110607525.0 / 20449.0);
    const auto fs_2295_143_30 = std::sqrt(158010750.0 / 20449.0);
    const auto fs_2295_143_35 = std::sqrt(184345875.0 / 20449.0);
    const auto fs_2295_143_6 = std::sqrt(31602150.0 / 20449.0);
    const auto fs_2295_286_10 = std::sqrt(26335125.0 / 40898.0);
    const auto fs_2295_286_70 = std::sqrt(184345875.0 / 40898.0);
    const auto fs_2352_2717_55 = std::sqrt(27659520.0 / 671099.0);
    const auto fs_2375_22_3 = std::sqrt(16921875.0 / 484.0);
    const auto fs_247_187_6 = std::sqrt(366054.0 / 34969.0);
    const auto fs_247_22_6 = std::sqrt(183027.0 / 242.0);
    const auto fs_250_11_14 = std::sqrt(875000.0 / 121.0);
    const auto fs_250_11_35 = std::sqrt(2187500.0 / 121.0);
    const auto fs_250_143_14 = std::sqrt(875000.0 / 20449.0);
    const auto fs_250_143_35 = std::sqrt(2187500.0 / 20449.0);
    const auto fs_252_46189_1001 = std::sqrt(444528.0 / 14919047.0);
    const auto fs_252_46189_12155 = std::sqrt(317520.0 / 877591.0);
    const auto fs_252_46189_2431 = std::sqrt(63504.0 / 877591.0);
    const auto fs_252_46189_3003 = std::sqrt(1333584.0 / 14919047.0);
    const auto fs_2646_2717_7 = std::sqrt(49009212.0 / 7382089.0);
    const auto fs_26_561_6 = std::sqrt(1352.0 / 104907.0);
    const auto fs_285_11_11 = std::sqrt(81225.0 / 11.0);
    const auto fs_285_11_7 = std::sqrt(568575.0 / 121.0);
    const auto fs_285_187_30 = std::sqrt(2436750.0 / 34969.0);
    const auto fs_285_22_154 = std::sqrt(568575.0 / 22.0);
    const auto fs_285_22_30 = std::sqrt(1218375.0 / 242.0);
    const auto fs_285_22_55 = std::sqrt(406125.0 / 44.0);
    const auto fs_285_22_70 = std::sqrt(2842875.0 / 242.0);
    const auto fs_285_44_10 = std::sqrt(406125.0 / 968.0);
    const auto fs_285_44_210 = std::sqrt(8528625.0 / 968.0);
    const auto fs_2898_46189_1001 = std::sqrt(58788828.0 / 14919047.0);
    const auto fs_2898_46189_12155 = std::sqrt(41992020.0 / 877591.0);
    const auto fs_2898_46189_2431 = std::sqrt(8398404.0 / 877591.0);
    const auto fs_2898_46189_3003 = std::sqrt(176366484.0 / 14919047.0);
    const auto fs_28_209_2 = std::sqrt(1568.0 / 43681.0);
    const auto fs_28_2717_429 = std::sqrt(2352.0 / 51623.0);
    const auto fs_28_2717_858 = std::sqrt(4704.0 / 51623.0);
    const auto fs_294_209_15 = std::sqrt(1296540.0 / 43681.0);
    const auto fs_294_247_55 = std::sqrt(4753980.0 / 61009.0);
    const auto fs_294_2717_1430 = std::sqrt(864360.0 / 51623.0);
    const auto fs_294_2717_5005 = std::sqrt(3025260.0 / 51623.0);
    const auto fs_2_187_10 = std::sqrt(40.0 / 34969.0);
    const auto fs_2_187_210 = std::sqrt(840.0 / 34969.0);
    const auto fs_2_51_42 = std::sqrt(56.0 / 867.0);
    const auto fs_3024_46189_77 = std::sqrt(64012032.0 / 193947611.0);
    const auto fs_304_11_3 = std::sqrt(277248.0 / 121.0);
    const auto fs_315_46189_2002 = std::sqrt(1389150.0 / 14919047.0);
    const auto fs_329_2717_6 = std::sqrt(649446.0 / 7382089.0);
    const auto fs_3375_16_30 = std::sqrt(170859375.0 / 128.0);
    const auto fs_34776_46189_77 = std::sqrt(8465591232.0 / 193947611.0);
    const auto fs_35_2717_154 = std::sqrt(17150.0 / 671099.0);
    const auto fs_35_2717_330 = std::sqrt(36750.0 / 671099.0);
    const auto fs_35_2717_858 = std::sqrt(7350.0 / 51623.0);
    const auto fs_35_429_6 = std::sqrt(2450.0 / 61347.0);
    const auto fs_375_22_30 = std::sqrt(2109375.0 / 242.0);
    const auto fs_375_286_30 = std::sqrt(2109375.0 / 40898.0);
    const auto fs_375_4_30 = std::sqrt(2109375.0 / 8.0);
    const auto fs_378_46189_2431 = std::sqrt(142884.0 / 877591.0);
    const auto fs_378_46189_3003 = std::sqrt(3000564.0 / 14919047.0);
    const auto fs_378_46189_330 = std::sqrt(4286520.0 / 193947611.0);
    const auto fs_378_46189_4290 = std::sqrt(4286520.0 / 14919047.0);
    const auto fs_378_46189_770 = std::sqrt(10001880.0 / 193947611.0);
    const auto fs_380_187_2 = std::sqrt(288800.0 / 34969.0);
    const auto fs_38_11_210 = std::sqrt(303240.0 / 121.0);
    const auto fs_38_11_462 = std::sqrt(60648.0 / 11.0);
    const auto fs_38_17_5 = std::sqrt(7220.0 / 289.0);
    const auto fs_38_187_14 = std::sqrt(20216.0 / 34969.0);
    const auto fs_38_187_231 = std::sqrt(30324.0 / 3179.0);
    const auto fs_3_143_10 = std::sqrt(90.0 / 20449.0);
    const auto fs_3_143_70 = std::sqrt(630.0 / 20449.0);
    const auto fs_408_143_5 = std::sqrt(832320.0 / 20449.0);
    const auto fs_40_561_2 = std::sqrt(3200.0 / 314721.0);
    const auto fs_4116_2717_15 = std::sqrt(254121840.0 / 7382089.0);
    const auto fs_42_2717_10 = std::sqrt(17640.0 / 7382089.0);
    const auto fs_42_2717_165 = std::sqrt(26460.0 / 671099.0);
    const auto fs_42_2717_715 = std::sqrt(8820.0 / 51623.0);
    const auto fs_4347_46189_2431 = std::sqrt(18896409.0 / 877591.0);
    const auto fs_4347_46189_3003 = std::sqrt(396824589.0 / 14919047.0);
    const auto fs_4347_46189_330 = std::sqrt(566892270.0 / 193947611.0);
    const auto fs_4347_46189_4290 = std::sqrt(566892270.0 / 14919047.0);
    const auto fs_4347_46189_770 = std::sqrt(1322748630.0 / 193947611.0);
    const auto fs_4347_92378_10010 = std::sqrt(661374315.0 / 29838094.0);
    const auto fs_4347_92378_110 = std::sqrt(94482045.0 / 387895222.0);
    const auto fs_441_143_70 = std::sqrt(13613670.0 / 20449.0);
    const auto fs_441_143_77 = std::sqrt(1361367.0 / 1859.0);
    const auto fs_441_286_10 = std::sqrt(972405.0 / 40898.0);
    const auto fs_441_286_165 = std::sqrt(2917215.0 / 7436.0);
    const auto fs_441_286_715 = std::sqrt(972405.0 / 572.0);
    const auto fs_456_187_5 = std::sqrt(1039680.0 / 34969.0);
    const auto fs_456_187_7 = std::sqrt(1455552.0 / 34969.0);
    const auto fs_4590_143_5 = std::sqrt(105340500.0 / 20449.0);
    const auto fs_4725_32_30 = std::sqrt(334884375.0 / 512.0);
    const auto fs_475_11_2 = std::sqrt(451250.0 / 121.0);
    const auto fs_475_11_3 = std::sqrt(676875.0 / 121.0);
    const auto fs_475_22_42 = std::sqrt(4738125.0 / 242.0);
    const auto fs_49_2717_210 = std::sqrt(504210.0 / 7382089.0);
    const auto fs_4_187_154 = std::sqrt(224.0 / 3179.0);
    const auto fs_4_187_30 = std::sqrt(480.0 / 34969.0);
    const auto fs_4_187_55 = std::sqrt(80.0 / 3179.0);
    const auto fs_4_187_70 = std::sqrt(1120.0 / 34969.0);
    const auto fs_4_51_5 = std::sqrt(80.0 / 2601.0);
    const auto fs_4_561_14 = std::sqrt(224.0 / 314721.0);
    const auto fs_4_561_231 = std::sqrt(112.0 / 9537.0);
    const auto fs_504_46189_1430 = std::sqrt(2540160.0 / 14919047.0);
    const auto fs_504_46189_2145 = std::sqrt(3810240.0 / 14919047.0);
    const auto fs_504_46189_55 = std::sqrt(1270080.0 / 193947611.0);
    const auto fs_50_429_2 = std::sqrt(5000.0 / 184041.0);
    const auto fs_50_429_3 = std::sqrt(2500.0 / 61347.0);
    const auto fs_5292_46189_33 = std::sqrt(84015792.0 / 193947611.0);
    const auto fs_5625_8_2 = std::sqrt(31640625.0 / 32.0);
    const auto fs_5625_8_3 = std::sqrt(94921875.0 / 64.0);
    const auto fs_570_11_5 = std::sqrt(1624500.0 / 121.0);
    const auto fs_570_11_7 = std::sqrt(2274300.0 / 121.0);
    const auto fs_570_187_11 = std::sqrt(324900.0 / 3179.0);
    const auto fs_570_187_7 = std::sqrt(2274300.0 / 34969.0);
    const auto fs_5796_46189_1430 = std::sqrt(335936160.0 / 14919047.0);
    const auto fs_5796_46189_2145 = std::sqrt(503904240.0 / 14919047.0);
    const auto fs_5796_46189_55 = std::sqrt(167968080.0 / 193947611.0);
    const auto fs_57_11_154 = std::sqrt(45486.0 / 11.0);
    const auto fs_57_11_30 = std::sqrt(97470.0 / 121.0);
    const auto fs_57_11_55 = std::sqrt(16245.0 / 11.0);
    const auto fs_57_11_70 = std::sqrt(227430.0 / 121.0);
    const auto fs_57_187_10 = std::sqrt(32490.0 / 34969.0);
    const auto fs_57_187_210 = std::sqrt(682290.0 / 34969.0);
    const auto fs_57_22_10 = std::sqrt(16245.0 / 242.0);
    const auto fs_57_22_210 = std::sqrt(341145.0 / 242.0);
    const auto fs_588_143_55 = std::sqrt(1728720.0 / 1859.0);
    const auto fs_588_209_2 = std::sqrt(691488.0 / 43681.0);
    const auto fs_588_2717_429 = std::sqrt(1037232.0 / 51623.0);
    const auto fs_588_2717_858 = std::sqrt(2074464.0 / 51623.0);
    const auto fs_5_143_30 = std::sqrt(750.0 / 20449.0);
    const auto fs_60858_46189_33 = std::sqrt(11111088492.0 / 193947611.0);
    const auto fs_608_187_3 = std::sqrt(1108992.0 / 34969.0);
    const auto fs_625_11_2 = std::sqrt(781250.0 / 121.0);
    const auto fs_625_11_3 = std::sqrt(1171875.0 / 121.0);
    const auto fs_625_143_2 = std::sqrt(781250.0 / 20449.0);
    const auto fs_625_143_3 = std::sqrt(1171875.0 / 20449.0);
    const auto fs_625_2_2 = std::sqrt(390625.0 / 2.0);
    const auto fs_625_2_3 = std::sqrt(1171875.0 / 4.0);
    const auto fs_63_46189_10010 = std::sqrt(277830.0 / 14919047.0);
    const auto fs_63_46189_2002 = std::sqrt(55566.0 / 14919047.0);
    const auto fs_63_46189_22 = std::sqrt(7938.0 / 193947611.0);
    const auto fs_64_561_3 = std::sqrt(4096.0 / 104907.0);
    const auto fs_6909_2717_6 = std::sqrt(286405686.0 / 7382089.0);
    const auto fs_6909_572_6 = std::sqrt(143202843.0 / 163592.0);
    const auto fs_6_143_15 = std::sqrt(540.0 / 20449.0);
    const auto fs_6_143_21 = std::sqrt(756.0 / 20449.0);
    const auto fs_6_143_30 = std::sqrt(1080.0 / 20449.0);
    const auto fs_6_143_35 = std::sqrt(1260.0 / 20449.0);
    const auto fs_6_143_6 = std::sqrt(216.0 / 20449.0);
    const auto fs_6_187_66 = std::sqrt(216.0 / 3179.0);
    const auto fs_70_2717_14 = std::sqrt(68600.0 / 7382089.0);
    const auto fs_70_2717_143 = std::sqrt(4900.0 / 51623.0);
    const auto fs_70_2717_286 = std::sqrt(9800.0 / 51623.0);
    const auto fs_70_2717_6 = std::sqrt(29400.0 / 7382089.0);
    const auto fs_70_2717_66 = std::sqrt(29400.0 / 671099.0);
    const auto fs_7245_92378_2002 = std::sqrt(367430175.0 / 29838094.0);
    const auto fs_735_2717_154 = std::sqrt(7563150.0 / 671099.0);
    const auto fs_735_2717_330 = std::sqrt(16206750.0 / 671099.0);
    const auto fs_735_2717_858 = std::sqrt(3241350.0 / 51623.0);
    const auto fs_735_286_14 = std::sqrt(3781575.0 / 40898.0);
    const auto fs_735_286_143 = std::sqrt(540225.0 / 572.0);
    const auto fs_735_286_286 = std::sqrt(540225.0 / 286.0);
    const auto fs_735_286_6 = std::sqrt(1620675.0 / 40898.0);
    const auto fs_735_286_66 = std::sqrt(1620675.0 / 3718.0);
    const auto fs_735_572_154 = std::sqrt(3781575.0 / 14872.0);
    const auto fs_735_572_330 = std::sqrt(8103375.0 / 14872.0);
    const auto fs_735_572_858 = std::sqrt(1620675.0 / 1144.0);
    const auto fs_756_46189_385 = std::sqrt(20003760.0 / 193947611.0);
    const auto fs_760_11_3 = std::sqrt(1732800.0 / 121.0);
    const auto fs_765_11_15 = std::sqrt(8778375.0 / 121.0);
    const auto fs_765_11_21 = std::sqrt(12289725.0 / 121.0);
    const auto fs_765_11_30 = std::sqrt(17556750.0 / 121.0);
    const auto fs_765_11_35 = std::sqrt(20482875.0 / 121.0);
    const auto fs_765_11_6 = std::sqrt(3511350.0 / 121.0);
    const auto fs_765_16_10 = std::sqrt(2926125.0 / 128.0);
    const auto fs_765_16_70 = std::sqrt(20482875.0 / 128.0);
    const auto fs_765_22_10 = std::sqrt(2926125.0 / 242.0);
    const auto fs_765_22_70 = std::sqrt(20482875.0 / 242.0);
    const auto fs_765_4_5 = std::sqrt(2926125.0 / 16.0);
    const auto fs_765_8_15 = std::sqrt(8778375.0 / 64.0);
    const auto fs_765_8_21 = std::sqrt(12289725.0 / 64.0);
    const auto fs_765_8_30 = std::sqrt(8778375.0 / 32.0);
    const auto fs_765_8_35 = std::sqrt(20482875.0 / 64.0);
    const auto fs_765_8_6 = std::sqrt(1755675.0 / 32.0);
    const auto fs_76_11_35 = std::sqrt(202160.0 / 121.0);
    const auto fs_76_187_210 = std::sqrt(1212960.0 / 34969.0);
    const auto fs_76_187_462 = std::sqrt(242592.0 / 3179.0);
    const auto fs_7875_16_2 = std::sqrt(62015625.0 / 128.0);
    const auto fs_7875_16_3 = std::sqrt(186046875.0 / 256.0);
    const auto fs_7875_16_6 = std::sqrt(186046875.0 / 128.0);
    const auto fs_7_2717_286 = std::sqrt(98.0 / 51623.0);
    const auto fs_7_2717_6006 = std::sqrt(2058.0 / 51623.0);
    const auto fs_7_2717_66 = std::sqrt(294.0 / 671099.0);
    const auto fs_84_2717_70 = std::sqrt(493920.0 / 7382089.0);
    const auto fs_84_2717_77 = std::sqrt(49392.0 / 671099.0);
    const auto fs_855_44_66 = std::sqrt(2193075.0 / 88.0);
    const auto fs_8694_46189_385 = std::sqrt(2645497260.0 / 193947611.0);
    const auto fs_875_22_6 = std::sqrt(2296875.0 / 242.0);
    const auto fs_875_286_6 = std::sqrt(2296875.0 / 40898.0);
    const auto fs_875_4_6 = std::sqrt(2296875.0 / 8.0);
    const auto fs_882_2717_10 = std::sqrt(7779240.0 / 7382089.0);
    const auto fs_882_2717_165 = std::sqrt(11668860.0 / 671099.0);
    const auto fs_882_2717_715 = std::sqrt(3889620.0 / 51623.0);
    const auto fs_8_561_210 = std::sqrt(4480.0 / 104907.0);
    const auto fs_8_561_462 = std::sqrt(896.0 / 9537.0);
    const auto fs_950_187_3 = std::sqrt(2707500.0 / 34969.0);
    const auto fs_95_11_210 = std::sqrt(1895250.0 / 121.0);
    const auto fs_95_11_42 = std::sqrt(379050.0 / 121.0);
    const auto fs_95_11_462 = std::sqrt(379050.0 / 11.0);
    const auto fs_95_22_14 = std::sqrt(63175.0 / 242.0);
    const auto fs_95_22_231 = std::sqrt(189525.0 / 44.0);
    const auto fs_95_2_5 = std::sqrt(45125.0 / 4.0);
    const auto fs_95_4_42 = std::sqrt(189525.0 / 8.0);

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ab_2, pc_0 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];

        pc_0[k] = e_0 * f_12285_32 + e_1 * f_23625_16 * h2_0 - e_1 * f_12285_8 * r_2 + e_2 * f_2295_4 * h4_0 - e_2 * f_16875_8 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 + e_3 * f_1425_22 * h6_0 - e_3 * f_4590_11 * r_2 * h4_0 + e_3 * f_1875_2 * r_4 * h2_0 - e_3 * f_585_1 * r_6 + e_4 * f_735_286 * h8_0 - e_4 * f_285_11 * r_2 * h6_0 + e_4 * f_13770_143 * r_4 * h4_0 - e_4 * f_1875_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 + e_5 * f_1449_46189 * h10_0 + e_5 * fs_1449_46189_92378 * h10_p10 - e_5 * f_1470_2717 * r_2 * h8_0 + e_5 * f_570_187 * r_4 * h6_0 - e_5 * f_1224_143 * r_6 * h4_0 + e_5 * f_1875_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 - e_6 * f_126_46189 * r_2 * h10_0 - e_6 * fs_126_46189_92378 * r_2 * h10_p10 + e_6 * f_70_2717 * r_4 * h8_0 - e_6 * f_20_187 * r_6 * h6_0 + e_6 * f_36_143 * r_8 * h4_0 - e_6 * f_50_143 * r_10 * h2_0 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ab_2, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_1[k] = - e_1 * fs_4725_32_30 * h2_p1 - e_2 * f_2295_4 * h4_p1 + e_2 * fs_3375_16_30 * r_2 * h2_p1 - e_3 * fs_285_44_210 * h6_p1 + e_3 * f_4590_11 * r_2 * h4_p1 - e_3 * fs_375_4_30 * r_4 * h2_p1 - e_4 * fs_441_286_10 * h8_p1 + e_4 * fs_57_22_210 * r_2 * h6_p1 - e_4 * f_13770_143 * r_4 * h4_p1 + e_4 * fs_375_22_30 * r_6 * h2_p1 - e_5 * fs_1449_92378_22 * h10_p1 + e_5 * fs_1449_46189_46189 * h10_p9 + e_5 * fs_882_2717_10 * r_2 * h8_p1 - e_5 * fs_57_187_210 * r_4 * h6_p1 + e_5 * f_1224_143 * r_6 * h4_p1 - e_5 * fs_375_286_30 * r_8 * h2_p1 + e_6 * fs_63_46189_22 * r_2 * h10_p1 - e_6 * fs_126_46189_46189 * r_2 * h10_p9 - e_6 * fs_42_2717_10 * r_4 * h8_p1 + e_6 * fs_2_187_210 * r_6 * h6_p1 - e_6 * f_36_143 * r_8 * h4_p1 + e_6 * fs_5_143_30 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_2[k] = e_1 * fs_1575_16_15 * h2_p2 + e_2 * f_2295_4 * h4_p2 - e_2 * fs_1125_8_15 * r_2 * h2_p2 + e_3 * fs_475_22_42 * h6_p2 - e_3 * f_4590_11 * r_2 * h4_p2 + e_3 * fs_125_2_15 * r_4 * h2_p2 + e_4 * fs_735_286_14 * h8_p2 + e_4 * fs_735_286_143 * h8_p8 - e_4 * fs_95_11_42 * r_2 * h6_p2 + e_4 * f_13770_143 * r_4 * h4_p2 - e_4 * fs_125_11_15 * r_6 * h2_p2 + e_5 * fs_1449_46189_33 * h10_p2 + e_5 * fs_4347_46189_2431 * h10_p8 - e_5 * fs_1470_2717_14 * r_2 * h8_p2 - e_5 * fs_1470_2717_143 * r_2 * h8_p8 + e_5 * fs_190_187_42 * r_4 * h6_p2 - e_5 * f_1224_143 * r_6 * h4_p2 + e_5 * fs_125_143_15 * r_8 * h2_p2 - e_6 * fs_126_46189_33 * r_2 * h10_p2 - e_6 * fs_378_46189_2431 * r_2 * h10_p8 + e_6 * fs_70_2717_14 * r_4 * h8_p2 + e_6 * fs_70_2717_143 * r_4 * h8_p8 - e_6 * fs_20_561_42 * r_6 * h6_p2 + e_6 * f_36_143 * r_8 * h4_p2 - e_6 * fs_10_429_15 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ab_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];

        pc_3[k] = - e_2 * fs_765_8_21 * h4_p3 - e_3 * fs_1425_22_7 * h6_p3 + e_3 * fs_765_11_21 * r_2 * h4_p3 - e_4 * fs_735_572_154 * h8_p3 + e_4 * fs_735_572_858 * h8_p7 + e_4 * fs_285_11_7 * r_2 * h6_p3 - e_4 * fs_2295_143_21 * r_4 * h4_p3 - e_5 * fs_1449_46189_143 * h10_p3 + e_5 * fs_2898_46189_2431 * h10_p7 + e_5 * fs_735_2717_154 * r_2 * h8_p3 - e_5 * fs_735_2717_858 * r_2 * h8_p7 - e_5 * fs_570_187_7 * r_4 * h6_p3 + e_5 * fs_204_143_21 * r_6 * h4_p3 + e_6 * fs_126_46189_143 * r_2 * h10_p3 - e_6 * fs_252_46189_2431 * r_2 * h10_p7 - e_6 * fs_35_2717_154 * r_4 * h8_p3 + e_6 * fs_35_2717_858 * r_4 * h8_p7 + e_6 * fs_20_187_7 * r_6 * h6_p3 - e_6 * fs_6_143_21 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p4, ph6_m5, ph6_p4, ph6_p6, ph8_m5, ph8_p4, ph8_p6, ph10_m5, ph10_p4, ph10_p6, ab_2, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];

        pc_4[k] = e_2 * fs_765_8_6 * h4_p4 + e_3 * fs_1425_44_30 * h6_p4 + e_3 * fs_285_22_55 * h6_p6 - e_3 * fs_765_11_6 * r_2 * h4_p4 + e_4 * fs_735_572_330 * h8_p4 + e_4 * fs_147_286_5005 * h8_p6 - e_4 * fs_285_22_30 * r_2 * h6_p4 - e_4 * fs_57_11_55 * r_2 * h6_p6 + e_4 * fs_2295_143_6 * r_4 * h4_p4 + e_5 * fs_1449_92378_2002 * h10_p4 + e_5 * fs_2898_46189_1001 * h10_p6 - e_5 * fs_735_2717_330 * r_2 * h8_p4 - e_5 * fs_294_2717_5005 * r_2 * h8_p6 + e_5 * fs_285_187_30 * r_4 * h6_p4 + e_5 * fs_114_187_55 * r_4 * h6_p6 - e_5 * fs_204_143_6 * r_6 * h4_p4 - e_6 * fs_63_46189_2002 * r_2 * h10_p4 - e_6 * fs_252_46189_1001 * r_2 * h10_p6 + e_6 * fs_35_2717_330 * r_4 * h8_p4 + e_6 * fs_14_2717_5005 * r_4 * h8_p6 - e_6 * fs_10_187_30 * r_6 * h6_p4 - e_6 * fs_4_187_55 * r_6 * h6_p6 + e_6 * fs_6_143_6 * r_8 * h4_p4;

        pc_5[k] = - e_3 * fs_1425_22_11 * h6_m5 - e_4 * fs_735_286_286 * h8_m5 + e_4 * fs_285_11_11 * r_2 * h6_m5 - e_5 * fs_1449_46189_3003 * h10_m5 + e_5 * fs_1470_2717_286 * r_2 * h8_m5 - e_5 * fs_570_187_11 * r_4 * h6_m5 + e_6 * fs_126_46189_3003 * r_2 * h10_m5 - e_6 * fs_70_2717_286 * r_4 * h8_m5 + e_6 * fs_20_187_11 * r_6 * h6_m5;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ab_2, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];

        pc_6[k] = e_2 * fs_765_8_6 * h4_m4 - e_3 * fs_285_22_55 * h6_m6 + e_3 * fs_1425_44_30 * h6_m4 - e_3 * fs_765_11_6 * r_2 * h4_m4 - e_4 * fs_147_286_5005 * h8_m6 + e_4 * fs_735_572_330 * h8_m4 + e_4 * fs_57_11_55 * r_2 * h6_m6 - e_4 * fs_285_22_30 * r_2 * h6_m4 + e_4 * fs_2295_143_6 * r_4 * h4_m4 - e_5 * fs_2898_46189_1001 * h10_m6 + e_5 * fs_1449_92378_2002 * h10_m4 + e_5 * fs_294_2717_5005 * r_2 * h8_m6 - e_5 * fs_735_2717_330 * r_2 * h8_m4 - e_5 * fs_114_187_55 * r_4 * h6_m6 + e_5 * fs_285_187_30 * r_4 * h6_m4 - e_5 * fs_204_143_6 * r_6 * h4_m4 + e_6 * fs_252_46189_1001 * r_2 * h10_m6 - e_6 * fs_63_46189_2002 * r_2 * h10_m4 - e_6 * fs_14_2717_5005 * r_4 * h8_m6 + e_6 * fs_35_2717_330 * r_4 * h8_m4 + e_6 * fs_4_187_55 * r_6 * h6_m6 - e_6 * fs_10_187_30 * r_6 * h6_m4 + e_6 * fs_6_143_6 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_7[k] = - e_2 * fs_765_8_21 * h4_m3 - e_3 * fs_1425_22_7 * h6_m3 + e_3 * fs_765_11_21 * r_2 * h4_m3 - e_4 * fs_735_572_858 * h8_m7 - e_4 * fs_735_572_154 * h8_m3 + e_4 * fs_285_11_7 * r_2 * h6_m3 - e_4 * fs_2295_143_21 * r_4 * h4_m3 - e_5 * fs_2898_46189_2431 * h10_m7 - e_5 * fs_1449_46189_143 * h10_m3 + e_5 * fs_735_2717_858 * r_2 * h8_m7 + e_5 * fs_735_2717_154 * r_2 * h8_m3 - e_5 * fs_570_187_7 * r_4 * h6_m3 + e_5 * fs_204_143_21 * r_6 * h4_m3 + e_6 * fs_252_46189_2431 * r_2 * h10_m7 + e_6 * fs_126_46189_143 * r_2 * h10_m3 - e_6 * fs_35_2717_858 * r_4 * h8_m7 - e_6 * fs_35_2717_154 * r_4 * h8_m3 + e_6 * fs_20_187_7 * r_6 * h6_m3 - e_6 * fs_6_143_21 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_8[k] = e_1 * fs_1575_16_15 * h2_m2 + e_2 * f_2295_4 * h4_m2 - e_2 * fs_1125_8_15 * r_2 * h2_m2 + e_3 * fs_475_22_42 * h6_m2 - e_3 * f_4590_11 * r_2 * h4_m2 + e_3 * fs_125_2_15 * r_4 * h2_m2 - e_4 * fs_735_286_143 * h8_m8 + e_4 * fs_735_286_14 * h8_m2 - e_4 * fs_95_11_42 * r_2 * h6_m2 + e_4 * f_13770_143 * r_4 * h4_m2 - e_4 * fs_125_11_15 * r_6 * h2_m2 - e_5 * fs_4347_46189_2431 * h10_m8 + e_5 * fs_1449_46189_33 * h10_m2 + e_5 * fs_1470_2717_143 * r_2 * h8_m8 - e_5 * fs_1470_2717_14 * r_2 * h8_m2 + e_5 * fs_190_187_42 * r_4 * h6_m2 - e_5 * f_1224_143 * r_6 * h4_m2 + e_5 * fs_125_143_15 * r_8 * h2_m2 + e_6 * fs_378_46189_2431 * r_2 * h10_m8 - e_6 * fs_126_46189_33 * r_2 * h10_m2 - e_6 * fs_70_2717_143 * r_4 * h8_m8 + e_6 * fs_70_2717_14 * r_4 * h8_m2 - e_6 * fs_20_561_42 * r_6 * h6_m2 + e_6 * f_36_143 * r_8 * h4_m2 - e_6 * fs_10_429_15 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m10, ph10_m9, ph10_m1, ab_2, pc_9, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];

        pc_9[k] = - e_1 * fs_4725_32_30 * h2_m1 - e_2 * f_2295_4 * h4_m1 + e_2 * fs_3375_16_30 * r_2 * h2_m1 - e_3 * fs_285_44_210 * h6_m1 + e_3 * f_4590_11 * r_2 * h4_m1 - e_3 * fs_375_4_30 * r_4 * h2_m1 - e_4 * fs_441_286_10 * h8_m1 + e_4 * fs_57_22_210 * r_2 * h6_m1 - e_4 * f_13770_143 * r_4 * h4_m1 + e_4 * fs_375_22_30 * r_6 * h2_m1 - e_5 * fs_1449_46189_46189 * h10_m9 - e_5 * fs_1449_92378_22 * h10_m1 + e_5 * fs_882_2717_10 * r_2 * h8_m1 - e_5 * fs_57_187_210 * r_4 * h6_m1 + e_5 * f_1224_143 * r_6 * h4_m1 - e_5 * fs_375_286_30 * r_8 * h2_m1 + e_6 * fs_126_46189_46189 * r_2 * h10_m9 + e_6 * fs_63_46189_22 * r_2 * h10_m1 - e_6 * fs_42_2717_10 * r_4 * h8_m1 + e_6 * fs_2_187_210 * r_6 * h6_m1 - e_6 * f_36_143 * r_8 * h4_m1 + e_6 * fs_5_143_30 * r_10 * h2_m1;

        pc_10[k] = - e_5 * fs_1449_46189_92378 * h10_m10 + e_6 * fs_126_46189_92378 * r_2 * h10_m10;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_11[k] = e_0 * f_12285_32 + e_1 * f_4725_8 * h2_0 - e_1 * f_12285_8 * r_2 - e_2 * f_2295_4 * h4_0 - e_2 * f_3375_4 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 - e_3 * f_2280_11 * h6_0 + e_3 * f_4590_11 * r_2 * h4_0 + e_3 * f_375_1 * r_4 * h2_0 - e_3 * f_585_1 * r_6 - e_4 * f_4557_286 * h8_0 - e_4 * fs_441_286_715 * h8_p8 + e_4 * f_912_11 * r_2 * h6_0 - e_4 * f_13770_143 * r_4 * h4_0 - e_4 * f_750_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 - e_5 * f_14490_46189 * h10_0 + e_5 * fs_2898_46189_12155 * h10_p8 + e_5 * f_9114_2717 * r_2 * h8_0 + e_5 * fs_882_2717_715 * r_2 * h8_p8 - e_5 * f_1824_187 * r_4 * h6_0 + e_5 * f_1224_143 * r_6 * h4_0 + e_5 * f_750_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 + e_6 * f_1260_46189 * r_2 * h10_0 - e_6 * fs_252_46189_12155 * r_2 * h10_p8 - e_6 * f_434_2717 * r_4 * h8_0 - e_6 * fs_42_2717_715 * r_4 * h8_p8 + e_6 * f_64_187 * r_6 * h6_0 - e_6 * f_36_143 * r_8 * h4_0 - e_6 * f_20_143 * r_10 * h2_0 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_12[k] = - e_1 * fs_11025_32_6 * h2_p1 + e_2 * fs_7875_16_6 * r_2 * h2_p1 + e_3 * fs_95_4_42 * h6_p1 - e_3 * fs_875_4_6 * r_4 * h2_p1 + e_4 * fs_147_11_2 * h8_p1 - e_4 * fs_147_286_1430 * h8_p7 - e_4 * fs_19_2_42 * r_2 * h6_p1 + e_4 * fs_875_22_6 * r_6 * h2_p1 + e_5 * fs_4347_92378_110 * h10_p1 + e_5 * fs_1449_46189_36465 * h10_p7 - e_5 * fs_588_209_2 * r_2 * h8_p1 + e_5 * fs_294_2717_1430 * r_2 * h8_p7 + e_5 * fs_19_17_42 * r_4 * h6_p1 - e_5 * fs_875_286_6 * r_8 * h2_p1 - e_6 * fs_189_46189_110 * r_2 * h10_p1 - e_6 * fs_126_46189_36465 * r_2 * h10_p7 + e_6 * fs_28_209_2 * r_4 * h8_p1 - e_6 * fs_14_2717_1430 * r_4 * h8_p7 - e_6 * fs_2_51_42 * r_6 * h6_p1 + e_6 * fs_35_429_6 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_13[k] = e_1 * f_4725_8 * h2_p2 + e_2 * fs_765_8_15 * h4_p2 - e_2 * f_3375_4 * r_2 * h2_p2 - e_3 * fs_285_22_70 * h6_p2 - e_3 * fs_285_22_154 * h6_p6 - e_3 * fs_765_11_15 * r_2 * h4_p2 + e_3 * f_375_1 * r_4 * h2_p2 - e_4 * fs_1029_572_210 * h8_p2 + e_4 * fs_147_572_286 * h8_p6 + e_4 * fs_57_11_70 * r_2 * h6_p2 + e_4 * fs_57_11_154 * r_2 * h6_p6 + e_4 * fs_2295_143_15 * r_4 * h4_p2 - e_4 * f_750_11 * r_6 * h2_p2 - e_5 * fs_5796_46189_55 * h10_p2 + e_5 * fs_5796_46189_1430 * h10_p6 + e_5 * fs_1029_2717_210 * r_2 * h8_p2 - e_5 * fs_147_2717_286 * r_2 * h8_p6 - e_5 * fs_114_187_70 * r_4 * h6_p2 - e_5 * fs_114_187_154 * r_4 * h6_p6 - e_5 * fs_204_143_15 * r_6 * h4_p2 + e_5 * f_750_143 * r_8 * h2_p2 + e_6 * fs_504_46189_55 * r_2 * h10_p2 - e_6 * fs_504_46189_1430 * r_2 * h10_p6 - e_6 * fs_49_2717_210 * r_4 * h8_p2 + e_6 * fs_7_2717_286 * r_4 * h8_p6 + e_6 * fs_4_187_70 * r_6 * h6_p2 + e_6 * fs_4_187_154 * r_6 * h6_p6 + e_6 * fs_6_143_15 * r_8 * h4_p2 - e_6 * f_20_143 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph4_p3, ph6_m4, ph6_p3, ph6_p5, ph8_m4, ph8_p3, ph8_p5, ph10_m4, ph10_p3, ph10_p5, ab_2, pc_14, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_14[k] = - e_2 * fs_765_8_30 * h4_p3 + e_3 * fs_285_44_10 * h6_p3 - e_3 * fs_855_44_66 * h6_p5 + e_3 * fs_765_11_30 * r_2 * h4_p3 + e_4 * fs_588_143_55 * h8_p3 + e_4 * fs_147_143_429 * h8_p5 - e_4 * fs_57_22_10 * r_2 * h6_p3 + e_4 * fs_171_22_66 * r_2 * h6_p5 - e_4 * fs_2295_143_30 * r_4 * h4_p3 + e_5 * fs_1449_92378_10010 * h10_p3 + e_5 * fs_7245_92378_2002 * h10_p5 - e_5 * fs_2352_2717_55 * r_2 * h8_p3 - e_5 * fs_588_2717_429 * r_2 * h8_p5 + e_5 * fs_57_187_10 * r_4 * h6_p3 - e_5 * fs_171_187_66 * r_4 * h6_p5 + e_5 * fs_204_143_30 * r_6 * h4_p3 - e_6 * fs_63_46189_10010 * r_2 * h10_p3 - e_6 * fs_315_46189_2002 * r_2 * h10_p5 + e_6 * fs_112_2717_55 * r_4 * h8_p3 + e_6 * fs_28_2717_429 * r_4 * h8_p5 - e_6 * fs_2_187_10 * r_6 * h6_p3 + e_6 * fs_6_187_66 * r_6 * h6_p5 - e_6 * fs_6_143_30 * r_8 * h4_p3;

        pc_15[k] = e_2 * f_2295_4 * h4_m4 + e_3 * fs_570_11_5 * h6_m4 - e_3 * f_4590_11 * r_2 * h4_m4 - e_4 * fs_147_26_55 * h8_m4 - e_4 * fs_228_11_5 * r_2 * h6_m4 + e_4 * f_13770_143 * r_4 * h4_m4 - e_5 * fs_2898_46189_3003 * h10_m4 + e_5 * fs_294_247_55 * r_2 * h8_m4 + e_5 * fs_456_187_5 * r_4 * h6_m4 - e_5 * f_1224_143 * r_6 * h4_m4 + e_6 * fs_252_46189_3003 * r_2 * h10_m4 - e_6 * fs_14_247_55 * r_4 * h8_m4 - e_6 * fs_16_187_5 * r_6 * h6_m4 + e_6 * f_36_143 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_16[k] = - e_2 * fs_765_8_30 * h4_m3 + e_3 * fs_855_44_66 * h6_m5 + e_3 * fs_285_44_10 * h6_m3 + e_3 * fs_765_11_30 * r_2 * h4_m3 - e_4 * fs_147_143_429 * h8_m5 + e_4 * fs_588_143_55 * h8_m3 - e_4 * fs_171_22_66 * r_2 * h6_m5 - e_4 * fs_57_22_10 * r_2 * h6_m3 - e_4 * fs_2295_143_30 * r_4 * h4_m3 - e_5 * fs_7245_92378_2002 * h10_m5 + e_5 * fs_1449_92378_10010 * h10_m3 + e_5 * fs_588_2717_429 * r_2 * h8_m5 - e_5 * fs_2352_2717_55 * r_2 * h8_m3 + e_5 * fs_171_187_66 * r_4 * h6_m5 + e_5 * fs_57_187_10 * r_4 * h6_m3 + e_5 * fs_204_143_30 * r_6 * h4_m3 + e_6 * fs_315_46189_2002 * r_2 * h10_m5 - e_6 * fs_63_46189_10010 * r_2 * h10_m3 - e_6 * fs_28_2717_429 * r_4 * h8_m5 + e_6 * fs_112_2717_55 * r_4 * h8_m3 - e_6 * fs_6_187_66 * r_6 * h6_m5 - e_6 * fs_2_187_10 * r_6 * h6_m3 - e_6 * fs_6_143_30 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_17[k] = e_1 * f_4725_8 * h2_m2 + e_2 * fs_765_8_15 * h4_m2 - e_2 * f_3375_4 * r_2 * h2_m2 + e_3 * fs_285_22_154 * h6_m6 - e_3 * fs_285_22_70 * h6_m2 - e_3 * fs_765_11_15 * r_2 * h4_m2 + e_3 * f_375_1 * r_4 * h2_m2 - e_4 * fs_147_572_286 * h8_m6 - e_4 * fs_1029_572_210 * h8_m2 - e_4 * fs_57_11_154 * r_2 * h6_m6 + e_4 * fs_57_11_70 * r_2 * h6_m2 + e_4 * fs_2295_143_15 * r_4 * h4_m2 - e_4 * f_750_11 * r_6 * h2_m2 - e_5 * fs_5796_46189_1430 * h10_m6 - e_5 * fs_5796_46189_55 * h10_m2 + e_5 * fs_147_2717_286 * r_2 * h8_m6 + e_5 * fs_1029_2717_210 * r_2 * h8_m2 + e_5 * fs_114_187_154 * r_4 * h6_m6 - e_5 * fs_114_187_70 * r_4 * h6_m2 - e_5 * fs_204_143_15 * r_6 * h4_m2 + e_5 * f_750_143 * r_8 * h2_m2 + e_6 * fs_504_46189_1430 * r_2 * h10_m6 + e_6 * fs_504_46189_55 * r_2 * h10_m2 - e_6 * fs_7_2717_286 * r_4 * h8_m6 - e_6 * fs_49_2717_210 * r_4 * h8_m2 - e_6 * fs_4_187_154 * r_6 * h6_m6 + e_6 * fs_4_187_70 * r_6 * h6_m2 + e_6 * fs_6_143_15 * r_8 * h4_m2 - e_6 * f_20_143 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m7, ph8_m1, ph10_m9, ph10_m8, ph10_m7, ph10_m1, ab_2, pc_18, pc_19, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_18[k] = - e_1 * fs_11025_32_6 * h2_m1 + e_2 * fs_7875_16_6 * r_2 * h2_m1 + e_3 * fs_95_4_42 * h6_m1 - e_3 * fs_875_4_6 * r_4 * h2_m1 + e_4 * fs_147_286_1430 * h8_m7 + e_4 * fs_147_11_2 * h8_m1 - e_4 * fs_19_2_42 * r_2 * h6_m1 + e_4 * fs_875_22_6 * r_6 * h2_m1 - e_5 * fs_1449_46189_36465 * h10_m7 + e_5 * fs_4347_92378_110 * h10_m1 - e_5 * fs_294_2717_1430 * r_2 * h8_m7 - e_5 * fs_588_209_2 * r_2 * h8_m1 + e_5 * fs_19_17_42 * r_4 * h6_m1 - e_5 * fs_875_286_6 * r_8 * h2_m1 + e_6 * fs_126_46189_36465 * r_2 * h10_m7 - e_6 * fs_189_46189_110 * r_2 * h10_m1 + e_6 * fs_14_2717_1430 * r_4 * h8_m7 + e_6 * fs_28_209_2 * r_4 * h8_m1 - e_6 * fs_2_51_42 * r_6 * h6_m1 + e_6 * fs_35_429_6 * r_10 * h2_m1;

        pc_19[k] = e_4 * fs_441_286_715 * h8_m8 - e_5 * fs_2898_46189_12155 * h10_m8 - e_5 * fs_882_2717_715 * r_2 * h8_m8 + e_6 * fs_252_46189_12155 * r_2 * h10_m8 + e_6 * fs_42_2717_715 * r_4 * h8_m8;

        pc_20[k] = e_1 * fs_4725_32_30 * h2_m1 + e_2 * f_2295_4 * h4_m1 - e_2 * fs_3375_16_30 * r_2 * h2_m1 + e_3 * fs_285_44_210 * h6_m1 - e_3 * f_4590_11 * r_2 * h4_m1 + e_3 * fs_375_4_30 * r_4 * h2_m1 + e_4 * fs_441_286_10 * h8_m1 - e_4 * fs_57_22_210 * r_2 * h6_m1 + e_4 * f_13770_143 * r_4 * h4_m1 - e_4 * fs_375_22_30 * r_6 * h2_m1 - e_5 * fs_1449_46189_46189 * h10_m9 + e_5 * fs_1449_92378_22 * h10_m1 - e_5 * fs_882_2717_10 * r_2 * h8_m1 + e_5 * fs_57_187_210 * r_4 * h6_m1 - e_5 * f_1224_143 * r_6 * h4_m1 + e_5 * fs_375_286_30 * r_8 * h2_m1 + e_6 * fs_126_46189_46189 * r_2 * h10_m9 - e_6 * fs_63_46189_22 * r_2 * h10_m1 + e_6 * fs_42_2717_10 * r_4 * h8_m1 - e_6 * fs_2_187_210 * r_6 * h6_m1 + e_6 * f_36_143 * r_8 * h4_m1 - e_6 * fs_5_143_30 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_21[k] = e_0 * f_12285_32 - e_1 * f_1575_16 * h2_0 - e_1 * f_12285_8 * r_2 - e_2 * f_2295_4 * h4_0 + e_2 * f_1125_8 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 + e_3 * f_2755_22 * h6_0 + e_3 * fs_95_11_462 * h6_p6 + e_3 * f_4590_11 * r_2 * h4_0 - e_3 * f_125_2 * r_4 * h2_0 - e_3 * f_585_1 * r_6 + e_4 * f_10731_286 * h8_0 - e_4 * fs_147_143_858 * h8_p6 - e_4 * f_551_11 * r_2 * h6_0 - e_4 * fs_38_11_462 * r_2 * h6_p6 - e_4 * f_13770_143 * r_4 * h4_0 + e_4 * f_125_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 + e_5 * f_65205_46189 * h10_0 + e_5 * fs_4347_46189_4290 * h10_p6 - e_5 * f_21462_2717 * r_2 * h8_0 + e_5 * fs_588_2717_858 * r_2 * h8_p6 + e_5 * f_1102_187 * r_4 * h6_0 + e_5 * fs_76_187_462 * r_4 * h6_p6 + e_5 * f_1224_143 * r_6 * h4_0 - e_5 * f_125_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 - e_6 * f_5670_46189 * r_2 * h10_0 - e_6 * fs_378_46189_4290 * r_2 * h10_p6 + e_6 * f_1022_2717 * r_4 * h8_0 - e_6 * fs_28_2717_858 * r_4 * h8_p6 - e_6 * f_116_561 * r_6 * h6_0 - e_6 * fs_8_561_462 * r_6 * h6_p6 - e_6 * f_36_143 * r_8 * h4_0 + e_6 * f_10_429 * r_10 * h2_0 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_22[k] = - e_1 * fs_7875_16_2 * h2_p1 + e_2 * fs_765_8_15 * h4_p1 + e_2 * fs_5625_8_2 * r_2 * h2_p1 + e_3 * fs_95_22_14 * h6_p1 + e_3 * fs_95_22_231 * h6_p5 - e_3 * fs_765_11_15 * r_2 * h4_p1 - e_3 * fs_625_2_2 * r_4 * h2_p1 - e_4 * fs_6909_572_6 * h8_p1 - e_4 * fs_147_572_6006 * h8_p5 - e_4 * fs_19_11_14 * r_2 * h6_p1 - e_4 * fs_19_11_231 * r_2 * h6_p5 + e_4 * fs_2295_143_15 * r_4 * h4_p1 + e_4 * fs_625_11_2 * r_6 * h2_p1 - e_5 * fs_4347_46189_330 * h10_p1 + e_5 * fs_21735_46189_143 * h10_p5 + e_5 * fs_6909_2717_6 * r_2 * h8_p1 + e_5 * fs_147_2717_6006 * r_2 * h8_p5 + e_5 * fs_38_187_14 * r_4 * h6_p1 + e_5 * fs_38_187_231 * r_4 * h6_p5 - e_5 * fs_204_143_15 * r_6 * h4_p1 - e_5 * fs_625_143_2 * r_8 * h2_p1 + e_6 * fs_378_46189_330 * r_2 * h10_p1 - e_6 * fs_1890_46189_143 * r_2 * h10_p5 - e_6 * fs_329_2717_6 * r_4 * h8_p1 - e_6 * fs_7_2717_6006 * r_4 * h8_p5 - e_6 * fs_4_561_14 * r_6 * h6_p1 - e_6 * fs_4_561_231 * r_6 * h6_p5 + e_6 * fs_6_143_15 * r_8 * h4_p1 + e_6 * fs_50_429_2 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_23[k] = e_1 * fs_1575_8_14 * h2_p2 + e_2 * fs_765_8_30 * h4_p4 - e_2 * fs_1125_4_14 * r_2 * h2_p2 - e_3 * fs_95_2_5 * h6_p2 - e_3 * fs_1235_44_6 * h6_p4 - e_3 * fs_765_11_30 * r_2 * h4_p4 + e_3 * fs_125_1_14 * r_4 * h2_p2 + e_4 * fs_147_22_15 * h8_p2 - e_4 * fs_147_572_66 * h8_p4 + e_4 * fs_19_1_5 * r_2 * h6_p2 + e_4 * fs_247_22_6 * r_2 * h6_p4 + e_4 * fs_2295_143_30 * r_4 * h4_p4 - e_4 * fs_250_11_14 * r_6 * h2_p2 + e_5 * fs_4347_46189_770 * h10_p2 + e_5 * fs_4347_92378_10010 * h10_p4 - e_5 * fs_294_209_15 * r_2 * h8_p2 + e_5 * fs_147_2717_66 * r_2 * h8_p4 - e_5 * fs_38_17_5 * r_4 * h6_p2 - e_5 * fs_247_187_6 * r_4 * h6_p4 - e_5 * fs_204_143_30 * r_6 * h4_p4 + e_5 * fs_250_143_14 * r_8 * h2_p2 - e_6 * fs_378_46189_770 * r_2 * h10_p2 - e_6 * fs_189_46189_10010 * r_2 * h10_p4 + e_6 * fs_14_209_15 * r_4 * h8_p2 - e_6 * fs_7_2717_66 * r_4 * h8_p4 + e_6 * fs_4_51_5 * r_6 * h6_p2 + e_6 * fs_26_561_6 * r_6 * h6_p4 + e_6 * fs_6_143_30 * r_8 * h4_p4 - e_6 * fs_20_429_14 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m3, ph8_m3, ph10_m3, ab_2, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m3 = ph10_m3[k];

        pc_24[k] = - e_2 * f_2295_4 * h4_m3 + e_3 * fs_2375_22_3 * h6_m3 + e_3 * f_4590_11 * r_2 * h4_m3 - e_4 * fs_735_286_66 * h8_m3 - e_4 * fs_475_11_3 * r_2 * h6_m3 - e_4 * f_13770_143 * r_4 * h4_m3 - e_5 * fs_4347_46189_3003 * h10_m3 + e_5 * fs_1470_2717_66 * r_2 * h8_m3 + e_5 * fs_950_187_3 * r_4 * h6_m3 + e_5 * f_1224_143 * r_6 * h4_m3 + e_6 * fs_378_46189_3003 * r_2 * h10_m3 - e_6 * fs_70_2717_66 * r_4 * h8_m3 - e_6 * fs_100_561_3 * r_6 * h6_m3 - e_6 * f_36_143 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m4, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_25[k] = e_1 * fs_1575_8_14 * h2_m2 - e_2 * fs_765_8_30 * h4_m4 - e_2 * fs_1125_4_14 * r_2 * h2_m2 + e_3 * fs_1235_44_6 * h6_m4 - e_3 * fs_95_2_5 * h6_m2 + e_3 * fs_765_11_30 * r_2 * h4_m4 + e_3 * fs_125_1_14 * r_4 * h2_m2 + e_4 * fs_147_572_66 * h8_m4 + e_4 * fs_147_22_15 * h8_m2 - e_4 * fs_247_22_6 * r_2 * h6_m4 + e_4 * fs_19_1_5 * r_2 * h6_m2 - e_4 * fs_2295_143_30 * r_4 * h4_m4 - e_4 * fs_250_11_14 * r_6 * h2_m2 - e_5 * fs_4347_92378_10010 * h10_m4 + e_5 * fs_4347_46189_770 * h10_m2 - e_5 * fs_147_2717_66 * r_2 * h8_m4 - e_5 * fs_294_209_15 * r_2 * h8_m2 + e_5 * fs_247_187_6 * r_4 * h6_m4 - e_5 * fs_38_17_5 * r_4 * h6_m2 + e_5 * fs_204_143_30 * r_6 * h4_m4 + e_5 * fs_250_143_14 * r_8 * h2_m2 + e_6 * fs_189_46189_10010 * r_2 * h10_m4 - e_6 * fs_378_46189_770 * r_2 * h10_m2 + e_6 * fs_7_2717_66 * r_4 * h8_m4 + e_6 * fs_14_209_15 * r_4 * h8_m2 - e_6 * fs_26_561_6 * r_6 * h6_m4 + e_6 * fs_4_51_5 * r_6 * h6_m2 - e_6 * fs_6_143_30 * r_8 * h4_m4 - e_6 * fs_20_429_14 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m6, ph6_m5, ph6_m1, ph8_m6, ph8_m5, ph8_m1, ph10_m6, ph10_m5, ph10_m1, ab_2, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_26[k] = - e_1 * fs_7875_16_2 * h2_m1 + e_2 * fs_765_8_15 * h4_m1 + e_2 * fs_5625_8_2 * r_2 * h2_m1 - e_3 * fs_95_22_231 * h6_m5 + e_3 * fs_95_22_14 * h6_m1 - e_3 * fs_765_11_15 * r_2 * h4_m1 - e_3 * fs_625_2_2 * r_4 * h2_m1 + e_4 * fs_147_572_6006 * h8_m5 - e_4 * fs_6909_572_6 * h8_m1 + e_4 * fs_19_11_231 * r_2 * h6_m5 - e_4 * fs_19_11_14 * r_2 * h6_m1 + e_4 * fs_2295_143_15 * r_4 * h4_m1 + e_4 * fs_625_11_2 * r_6 * h2_m1 - e_5 * fs_21735_46189_143 * h10_m5 - e_5 * fs_4347_46189_330 * h10_m1 - e_5 * fs_147_2717_6006 * r_2 * h8_m5 + e_5 * fs_6909_2717_6 * r_2 * h8_m1 - e_5 * fs_38_187_231 * r_4 * h6_m5 + e_5 * fs_38_187_14 * r_4 * h6_m1 - e_5 * fs_204_143_15 * r_6 * h4_m1 - e_5 * fs_625_143_2 * r_8 * h2_m1 + e_6 * fs_1890_46189_143 * r_2 * h10_m5 + e_6 * fs_378_46189_330 * r_2 * h10_m1 + e_6 * fs_7_2717_6006 * r_4 * h8_m5 - e_6 * fs_329_2717_6 * r_4 * h8_m1 + e_6 * fs_4_561_231 * r_6 * h6_m5 - e_6 * fs_4_561_14 * r_6 * h6_m1 + e_6 * fs_6_143_15 * r_8 * h4_m1 + e_6 * fs_50_429_2 * r_10 * h2_m1;

        pc_27[k] = - e_3 * fs_95_11_462 * h6_m6 + e_4 * fs_147_143_858 * h8_m6 + e_4 * fs_38_11_462 * r_2 * h6_m6 - e_5 * fs_4347_46189_4290 * h10_m6 - e_5 * fs_588_2717_858 * r_2 * h8_m6 - e_5 * fs_76_187_462 * r_4 * h6_m6 + e_6 * fs_378_46189_4290 * r_2 * h10_m6 + e_6 * fs_28_2717_858 * r_4 * h8_m6 + e_6 * fs_8_561_462 * r_6 * h6_m6;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ab_2, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_28[k] = e_1 * fs_11025_32_6 * h2_m1 - e_2 * fs_7875_16_6 * r_2 * h2_m1 - e_3 * fs_95_4_42 * h6_m1 + e_3 * fs_875_4_6 * r_4 * h2_m1 + e_4 * fs_147_286_1430 * h8_m7 - e_4 * fs_147_11_2 * h8_m1 + e_4 * fs_19_2_42 * r_2 * h6_m1 - e_4 * fs_875_22_6 * r_6 * h2_m1 - e_5 * fs_1449_46189_36465 * h10_m7 - e_5 * fs_4347_92378_110 * h10_m1 - e_5 * fs_294_2717_1430 * r_2 * h8_m7 + e_5 * fs_588_209_2 * r_2 * h8_m1 - e_5 * fs_19_17_42 * r_4 * h6_m1 + e_5 * fs_875_286_6 * r_8 * h2_m1 + e_6 * fs_126_46189_36465 * r_2 * h10_m7 + e_6 * fs_189_46189_110 * r_2 * h10_m1 + e_6 * fs_14_2717_1430 * r_4 * h8_m7 - e_6 * fs_28_209_2 * r_4 * h8_m1 + e_6 * fs_2_51_42 * r_6 * h6_m1 - e_6 * fs_35_429_6 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_29[k] = - e_1 * fs_1575_16_15 * h2_m2 - e_2 * f_2295_4 * h4_m2 + e_2 * fs_1125_8_15 * r_2 * h2_m2 - e_3 * fs_475_22_42 * h6_m2 + e_3 * f_4590_11 * r_2 * h4_m2 - e_3 * fs_125_2_15 * r_4 * h2_m2 - e_4 * fs_735_286_143 * h8_m8 - e_4 * fs_735_286_14 * h8_m2 + e_4 * fs_95_11_42 * r_2 * h6_m2 - e_4 * f_13770_143 * r_4 * h4_m2 + e_4 * fs_125_11_15 * r_6 * h2_m2 - e_5 * fs_4347_46189_2431 * h10_m8 - e_5 * fs_1449_46189_33 * h10_m2 + e_5 * fs_1470_2717_143 * r_2 * h8_m8 + e_5 * fs_1470_2717_14 * r_2 * h8_m2 - e_5 * fs_190_187_42 * r_4 * h6_m2 + e_5 * f_1224_143 * r_6 * h4_m2 - e_5 * fs_125_143_15 * r_8 * h2_m2 + e_6 * fs_378_46189_2431 * r_2 * h10_m8 + e_6 * fs_126_46189_33 * r_2 * h10_m2 - e_6 * fs_70_2717_143 * r_4 * h8_m8 - e_6 * fs_70_2717_14 * r_4 * h8_m2 + e_6 * fs_20_561_42 * r_6 * h6_m2 - e_6 * f_36_143 * r_8 * h4_m2 + e_6 * fs_10_429_15 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_30[k] = e_0 * f_12285_32 - e_1 * f_4725_8 * h2_0 - e_1 * f_12285_8 * r_2 - e_2 * f_765_8 * h4_0 - e_2 * fs_765_8_35 * h4_p4 + e_2 * f_3375_4 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 + e_3 * f_1710_11 * h6_0 + e_3 * fs_570_11_7 * h6_p4 + e_3 * f_765_11 * r_2 * h4_0 + e_3 * fs_765_11_35 * r_2 * h4_p4 - e_3 * f_375_1 * r_4 * h2_0 - e_3 * f_585_1 * r_6 - e_4 * f_4998_143 * h8_0 - e_4 * fs_441_143_77 * h8_p4 - e_4 * f_684_11 * r_2 * h6_0 - e_4 * fs_228_11_7 * r_2 * h6_p4 - e_4 * f_2295_143 * r_4 * h4_0 - e_4 * fs_2295_143_35 * r_4 * h4_p4 + e_4 * f_750_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 - e_5 * f_173880_46189 * h10_0 + e_5 * fs_5796_46189_2145 * h10_p4 + e_5 * f_19992_2717 * r_2 * h8_0 + e_5 * fs_1764_2717_77 * r_2 * h8_p4 + e_5 * f_1368_187 * r_4 * h6_0 + e_5 * fs_456_187_7 * r_4 * h6_p4 + e_5 * f_204_143 * r_6 * h4_0 + e_5 * fs_204_143_35 * r_6 * h4_p4 - e_5 * f_750_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 + e_6 * f_15120_46189 * r_2 * h10_0 - e_6 * fs_504_46189_2145 * r_2 * h10_p4 - e_6 * f_952_2717 * r_4 * h8_0 - e_6 * fs_84_2717_77 * r_4 * h8_p4 - e_6 * f_48_187 * r_6 * h6_0 - e_6 * fs_16_187_7 * r_6 * h6_p4 - e_6 * f_6_143 * r_8 * h4_0 - e_6 * fs_6_143_35 * r_8 * h4_p4 + e_6 * f_20_143 * r_10 * h2_0 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_31[k] = - e_1 * fs_1575_16_21 * h2_p1 + e_2 * fs_765_16_70 * h4_p1 - e_2 * fs_765_16_10 * h4_p3 + e_2 * fs_1125_8_21 * r_2 * h2_p1 - e_3 * fs_760_11_3 * h6_p1 + e_3 * fs_285_22_30 * h6_p3 - e_3 * fs_765_22_70 * r_2 * h4_p1 + e_3 * fs_765_22_10 * r_2 * h4_p3 - e_3 * fs_125_2_21 * r_4 * h2_p1 + e_4 * fs_1323_286_7 * h8_p1 - e_4 * fs_441_286_165 * h8_p3 + e_4 * fs_304_11_3 * r_2 * h6_p1 - e_4 * fs_57_11_30 * r_2 * h6_p3 + e_4 * fs_2295_286_70 * r_4 * h4_p1 - e_4 * fs_2295_286_10 * r_4 * h4_p3 + e_4 * fs_125_11_21 * r_6 * h2_p1 + e_5 * fs_8694_46189_385 * h10_p1 + e_5 * fs_1449_46189_30030 * h10_p3 - e_5 * fs_2646_2717_7 * r_2 * h8_p1 + e_5 * fs_882_2717_165 * r_2 * h8_p3 - e_5 * fs_608_187_3 * r_4 * h6_p1 + e_5 * fs_114_187_30 * r_4 * h6_p3 - e_5 * fs_102_143_70 * r_6 * h4_p1 + e_5 * fs_102_143_10 * r_6 * h4_p3 - e_5 * fs_125_143_21 * r_8 * h2_p1 - e_6 * fs_756_46189_385 * r_2 * h10_p1 - e_6 * fs_126_46189_30030 * r_2 * h10_p3 + e_6 * fs_126_2717_7 * r_4 * h8_p1 - e_6 * fs_42_2717_165 * r_4 * h8_p3 + e_6 * fs_64_561_3 * r_6 * h6_p1 - e_6 * fs_4_187_30 * r_6 * h6_p3 + e_6 * fs_3_143_70 * r_8 * h4_p1 - e_6 * fs_3_143_10 * r_8 * h4_p3 + e_6 * fs_10_429_21 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ab_2, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];

        pc_32[k] = e_1 * fs_1575_8_35 * h2_m2 - e_2 * fs_765_8_21 * h4_m2 - e_2 * fs_1125_4_35 * r_2 * h2_m2 + e_3 * fs_475_11_2 * h6_m2 + e_3 * fs_765_11_21 * r_2 * h4_m2 + e_3 * fs_125_1_35 * r_4 * h2_m2 + e_4 * fs_735_286_6 * h8_m2 - e_4 * fs_190_11_2 * r_2 * h6_m2 - e_4 * fs_2295_143_21 * r_4 * h4_m2 - e_4 * fs_250_11_35 * r_6 * h2_m2 - e_5 * fs_34776_46189_77 * h10_m2 - e_5 * fs_1470_2717_6 * r_2 * h8_m2 + e_5 * fs_380_187_2 * r_4 * h6_m2 + e_5 * fs_204_143_21 * r_6 * h4_m2 + e_5 * fs_250_143_35 * r_8 * h2_m2 + e_6 * fs_3024_46189_77 * r_2 * h10_m2 + e_6 * fs_70_2717_6 * r_4 * h8_m2 - e_6 * fs_40_561_2 * r_6 * h6_m2 - e_6 * fs_6_143_21 * r_8 * h4_m2 - e_6 * fs_20_429_35 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_33[k] = - e_1 * fs_1575_16_21 * h2_m1 + e_2 * fs_765_16_10 * h4_m3 + e_2 * fs_765_16_70 * h4_m1 + e_2 * fs_1125_8_21 * r_2 * h2_m1 - e_3 * fs_285_22_30 * h6_m3 - e_3 * fs_760_11_3 * h6_m1 - e_3 * fs_765_22_10 * r_2 * h4_m3 - e_3 * fs_765_22_70 * r_2 * h4_m1 - e_3 * fs_125_2_21 * r_4 * h2_m1 + e_4 * fs_441_286_165 * h8_m3 + e_4 * fs_1323_286_7 * h8_m1 + e_4 * fs_57_11_30 * r_2 * h6_m3 + e_4 * fs_304_11_3 * r_2 * h6_m1 + e_4 * fs_2295_286_10 * r_4 * h4_m3 + e_4 * fs_2295_286_70 * r_4 * h4_m1 + e_4 * fs_125_11_21 * r_6 * h2_m1 - e_5 * fs_1449_46189_30030 * h10_m3 + e_5 * fs_8694_46189_385 * h10_m1 - e_5 * fs_882_2717_165 * r_2 * h8_m3 - e_5 * fs_2646_2717_7 * r_2 * h8_m1 - e_5 * fs_114_187_30 * r_4 * h6_m3 - e_5 * fs_608_187_3 * r_4 * h6_m1 - e_5 * fs_102_143_10 * r_6 * h4_m3 - e_5 * fs_102_143_70 * r_6 * h4_m1 - e_5 * fs_125_143_21 * r_8 * h2_m1 + e_6 * fs_126_46189_30030 * r_2 * h10_m3 - e_6 * fs_756_46189_385 * r_2 * h10_m1 + e_6 * fs_42_2717_165 * r_4 * h8_m3 + e_6 * fs_126_2717_7 * r_4 * h8_m1 + e_6 * fs_4_187_30 * r_6 * h6_m3 + e_6 * fs_64_561_3 * r_6 * h6_m1 + e_6 * fs_3_143_10 * r_8 * h4_m3 + e_6 * fs_3_143_70 * r_8 * h4_m1 + e_6 * fs_10_429_21 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ab_2, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];

        pc_34[k] = e_2 * fs_765_8_35 * h4_m4 - e_3 * fs_570_11_7 * h6_m4 - e_3 * fs_765_11_35 * r_2 * h4_m4 + e_4 * fs_441_143_77 * h8_m4 + e_4 * fs_228_11_7 * r_2 * h6_m4 + e_4 * fs_2295_143_35 * r_4 * h4_m4 - e_5 * fs_5796_46189_2145 * h10_m4 - e_5 * fs_1764_2717_77 * r_2 * h8_m4 - e_5 * fs_456_187_7 * r_4 * h6_m4 - e_5 * fs_204_143_35 * r_6 * h4_m4 + e_6 * fs_504_46189_2145 * r_2 * h10_m4 + e_6 * fs_84_2717_77 * r_4 * h8_m4 + e_6 * fs_16_187_7 * r_6 * h6_m4 + e_6 * fs_6_143_35 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ab_2, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_35[k] = e_1 * fs_7875_16_2 * h2_m1 - e_2 * fs_765_8_15 * h4_m1 - e_2 * fs_5625_8_2 * r_2 * h2_m1 - e_3 * fs_95_22_231 * h6_m5 - e_3 * fs_95_22_14 * h6_m1 + e_3 * fs_765_11_15 * r_2 * h4_m1 + e_3 * fs_625_2_2 * r_4 * h2_m1 + e_4 * fs_147_572_6006 * h8_m5 + e_4 * fs_6909_572_6 * h8_m1 + e_4 * fs_19_11_231 * r_2 * h6_m5 + e_4 * fs_19_11_14 * r_2 * h6_m1 - e_4 * fs_2295_143_15 * r_4 * h4_m1 - e_4 * fs_625_11_2 * r_6 * h2_m1 - e_5 * fs_21735_46189_143 * h10_m5 + e_5 * fs_4347_46189_330 * h10_m1 - e_5 * fs_147_2717_6006 * r_2 * h8_m5 - e_5 * fs_6909_2717_6 * r_2 * h8_m1 - e_5 * fs_38_187_231 * r_4 * h6_m5 - e_5 * fs_38_187_14 * r_4 * h6_m1 + e_5 * fs_204_143_15 * r_6 * h4_m1 + e_5 * fs_625_143_2 * r_8 * h2_m1 + e_6 * fs_1890_46189_143 * r_2 * h10_m5 - e_6 * fs_378_46189_330 * r_2 * h10_m1 + e_6 * fs_7_2717_6006 * r_4 * h8_m5 + e_6 * fs_329_2717_6 * r_4 * h8_m1 + e_6 * fs_4_561_231 * r_6 * h6_m5 + e_6 * fs_4_561_14 * r_6 * h6_m1 - e_6 * fs_6_143_15 * r_8 * h4_m1 - e_6 * fs_50_429_2 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_36[k] = - e_1 * f_4725_8 * h2_m2 - e_2 * fs_765_8_15 * h4_m2 + e_2 * f_3375_4 * r_2 * h2_m2 + e_3 * fs_285_22_154 * h6_m6 + e_3 * fs_285_22_70 * h6_m2 + e_3 * fs_765_11_15 * r_2 * h4_m2 - e_3 * f_375_1 * r_4 * h2_m2 - e_4 * fs_147_572_286 * h8_m6 + e_4 * fs_1029_572_210 * h8_m2 - e_4 * fs_57_11_154 * r_2 * h6_m6 - e_4 * fs_57_11_70 * r_2 * h6_m2 - e_4 * fs_2295_143_15 * r_4 * h4_m2 + e_4 * f_750_11 * r_6 * h2_m2 - e_5 * fs_5796_46189_1430 * h10_m6 + e_5 * fs_5796_46189_55 * h10_m2 + e_5 * fs_147_2717_286 * r_2 * h8_m6 - e_5 * fs_1029_2717_210 * r_2 * h8_m2 + e_5 * fs_114_187_154 * r_4 * h6_m6 + e_5 * fs_114_187_70 * r_4 * h6_m2 + e_5 * fs_204_143_15 * r_6 * h4_m2 - e_5 * f_750_143 * r_8 * h2_m2 + e_6 * fs_504_46189_1430 * r_2 * h10_m6 - e_6 * fs_504_46189_55 * r_2 * h10_m2 - e_6 * fs_7_2717_286 * r_4 * h8_m6 + e_6 * fs_49_2717_210 * r_4 * h8_m2 - e_6 * fs_4_187_154 * r_6 * h6_m6 - e_6 * fs_4_187_70 * r_6 * h6_m2 - e_6 * fs_6_143_15 * r_8 * h4_m2 + e_6 * f_20_143 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_37[k] = e_2 * fs_765_8_21 * h4_m3 + e_3 * fs_1425_22_7 * h6_m3 - e_3 * fs_765_11_21 * r_2 * h4_m3 - e_4 * fs_735_572_858 * h8_m7 + e_4 * fs_735_572_154 * h8_m3 - e_4 * fs_285_11_7 * r_2 * h6_m3 + e_4 * fs_2295_143_21 * r_4 * h4_m3 - e_5 * fs_2898_46189_2431 * h10_m7 + e_5 * fs_1449_46189_143 * h10_m3 + e_5 * fs_735_2717_858 * r_2 * h8_m7 - e_5 * fs_735_2717_154 * r_2 * h8_m3 + e_5 * fs_570_187_7 * r_4 * h6_m3 - e_5 * fs_204_143_21 * r_6 * h4_m3 + e_6 * fs_252_46189_2431 * r_2 * h10_m7 - e_6 * fs_126_46189_143 * r_2 * h10_m3 - e_6 * fs_35_2717_858 * r_4 * h8_m7 + e_6 * fs_35_2717_154 * r_4 * h8_m3 - e_6 * fs_20_187_7 * r_6 * h6_m3 + e_6 * fs_6_143_21 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p2 = ph10_p2[k];

        pc_38[k] = e_0 * f_12285_32 - e_1 * f_14175_16 * h2_0 + e_1 * fs_7875_16_3 * h2_p2 - e_1 * f_12285_8 * r_2 + e_2 * f_765_2 * h4_0 - e_2 * fs_765_4_5 * h4_p2 + e_2 * f_10125_8 * r_2 * h2_0 - e_2 * fs_5625_8_3 * r_2 * h2_p2 + e_2 * f_12285_8 * r_4 - e_3 * f_570_11 * h6_0 + e_3 * fs_95_11_210 * h6_p2 - e_3 * f_3060_11 * r_2 * h4_0 + e_3 * fs_1530_11_5 * r_2 * h4_p2 - e_3 * f_1125_2 * r_4 * h2_0 + e_3 * fs_625_2_3 * r_4 * h2_p2 - e_3 * f_585_1 * r_6 - e_4 * f_1029_143 * h8_0 - e_4 * fs_441_143_70 * h8_p2 + e_4 * f_228_11 * r_2 * h6_0 - e_4 * fs_38_11_210 * r_2 * h6_p2 + e_4 * f_9180_143 * r_4 * h4_0 - e_4 * fs_4590_143_5 * r_4 * h4_p2 + e_4 * f_1125_11 * r_6 * h2_0 - e_4 * fs_625_11_3 * r_6 * h2_p2 + e_4 * f_195_2 * r_8 + e_5 * f_304290_46189 * h10_0 + e_5 * fs_20286_46189_165 * h10_p2 + e_5 * f_4116_2717 * r_2 * h8_0 + e_5 * fs_1764_2717_70 * r_2 * h8_p2 - e_5 * f_456_187 * r_4 * h6_0 + e_5 * fs_76_187_210 * r_4 * h6_p2 - e_5 * f_816_143 * r_6 * h4_0 + e_5 * fs_408_143_5 * r_6 * h4_p2 - e_5 * f_1125_143 * r_8 * h2_0 + e_5 * fs_625_143_3 * r_8 * h2_p2 - e_5 * f_78_11 * r_10 - e_6 * f_26460_46189 * r_2 * h10_0 - e_6 * fs_1764_46189_165 * r_2 * h10_p2 - e_6 * f_196_2717 * r_4 * h8_0 - e_6 * fs_84_2717_70 * r_4 * h8_p2 + e_6 * f_16_187 * r_6 * h6_0 - e_6 * fs_8_561_210 * r_6 * h6_p2 + e_6 * f_24_143 * r_8 * h4_0 - e_6 * fs_12_143_5 * r_8 * h4_p2 + e_6 * f_30_143 * r_10 * h2_0 - e_6 * fs_50_429_3 * r_10 * h2_p2 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m2, ph10_m1, ab_2, pc_39, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_m1 = ph10_m1[k];

        pc_39[k] = - e_1 * fs_1575_16_5 * h2_m1 + e_2 * fs_765_8_6 * h4_m1 + e_2 * fs_1125_8_5 * r_2 * h2_m1 - e_3 * fs_190_11_35 * h6_m1 - e_3 * fs_765_11_6 * r_2 * h4_m1 - e_3 * fs_125_2_5 * r_4 * h2_m1 + e_4 * fs_1029_143_15 * h8_m1 + e_4 * fs_76_11_35 * r_2 * h6_m1 + e_4 * fs_2295_143_6 * r_4 * h4_m1 + e_4 * fs_125_11_5 * r_6 * h2_m1 - e_5 * fs_60858_46189_33 * h10_m1 - e_5 * fs_4116_2717_15 * r_2 * h8_m1 - e_5 * fs_152_187_35 * r_4 * h6_m1 - e_5 * fs_204_143_6 * r_6 * h4_m1 - e_5 * fs_125_143_5 * r_8 * h2_m1 + e_6 * fs_5292_46189_33 * r_2 * h10_m1 + e_6 * fs_196_2717_15 * r_4 * h8_m1 + e_6 * fs_16_561_35 * r_6 * h6_m1 + e_6 * fs_6_143_6 * r_8 * h4_m1 + e_6 * fs_10_429_5 * r_10 * h2_m1;

        pc_40[k] = - e_1 * fs_7875_16_3 * h2_m2 + e_2 * fs_765_4_5 * h4_m2 + e_2 * fs_5625_8_3 * r_2 * h2_m2 - e_3 * fs_95_11_210 * h6_m2 - e_3 * fs_1530_11_5 * r_2 * h4_m2 - e_3 * fs_625_2_3 * r_4 * h2_m2 + e_4 * fs_441_143_70 * h8_m2 + e_4 * fs_38_11_210 * r_2 * h6_m2 + e_4 * fs_4590_143_5 * r_4 * h4_m2 + e_4 * fs_625_11_3 * r_6 * h2_m2 - e_5 * fs_20286_46189_165 * h10_m2 - e_5 * fs_1764_2717_70 * r_2 * h8_m2 - e_5 * fs_76_187_210 * r_4 * h6_m2 - e_5 * fs_408_143_5 * r_6 * h4_m2 - e_5 * fs_625_143_3 * r_8 * h2_m2 + e_6 * fs_1764_46189_165 * r_2 * h10_m2 + e_6 * fs_84_2717_70 * r_4 * h8_m2 + e_6 * fs_8_561_210 * r_6 * h6_m2 + e_6 * fs_12_143_5 * r_8 * h4_m2 + e_6 * fs_50_429_3 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_41[k] = e_1 * fs_1575_16_21 * h2_m1 + e_2 * fs_765_16_10 * h4_m3 - e_2 * fs_765_16_70 * h4_m1 - e_2 * fs_1125_8_21 * r_2 * h2_m1 - e_3 * fs_285_22_30 * h6_m3 + e_3 * fs_760_11_3 * h6_m1 - e_3 * fs_765_22_10 * r_2 * h4_m3 + e_3 * fs_765_22_70 * r_2 * h4_m1 + e_3 * fs_125_2_21 * r_4 * h2_m1 + e_4 * fs_441_286_165 * h8_m3 - e_4 * fs_1323_286_7 * h8_m1 + e_4 * fs_57_11_30 * r_2 * h6_m3 - e_4 * fs_304_11_3 * r_2 * h6_m1 + e_4 * fs_2295_286_10 * r_4 * h4_m3 - e_4 * fs_2295_286_70 * r_4 * h4_m1 - e_4 * fs_125_11_21 * r_6 * h2_m1 - e_5 * fs_1449_46189_30030 * h10_m3 - e_5 * fs_8694_46189_385 * h10_m1 - e_5 * fs_882_2717_165 * r_2 * h8_m3 + e_5 * fs_2646_2717_7 * r_2 * h8_m1 - e_5 * fs_114_187_30 * r_4 * h6_m3 + e_5 * fs_608_187_3 * r_4 * h6_m1 - e_5 * fs_102_143_10 * r_6 * h4_m3 + e_5 * fs_102_143_70 * r_6 * h4_m1 + e_5 * fs_125_143_21 * r_8 * h2_m1 + e_6 * fs_126_46189_30030 * r_2 * h10_m3 + e_6 * fs_756_46189_385 * r_2 * h10_m1 + e_6 * fs_42_2717_165 * r_4 * h8_m3 - e_6 * fs_126_2717_7 * r_4 * h8_m1 + e_6 * fs_4_187_30 * r_6 * h6_m3 - e_6 * fs_64_561_3 * r_6 * h6_m1 + e_6 * fs_3_143_10 * r_8 * h4_m3 - e_6 * fs_3_143_70 * r_8 * h4_m1 - e_6 * fs_10_429_21 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m4, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_42[k] = - e_1 * fs_1575_8_14 * h2_m2 - e_2 * fs_765_8_30 * h4_m4 + e_2 * fs_1125_4_14 * r_2 * h2_m2 + e_3 * fs_1235_44_6 * h6_m4 + e_3 * fs_95_2_5 * h6_m2 + e_3 * fs_765_11_30 * r_2 * h4_m4 - e_3 * fs_125_1_14 * r_4 * h2_m2 + e_4 * fs_147_572_66 * h8_m4 - e_4 * fs_147_22_15 * h8_m2 - e_4 * fs_247_22_6 * r_2 * h6_m4 - e_4 * fs_19_1_5 * r_2 * h6_m2 - e_4 * fs_2295_143_30 * r_4 * h4_m4 + e_4 * fs_250_11_14 * r_6 * h2_m2 - e_5 * fs_4347_92378_10010 * h10_m4 - e_5 * fs_4347_46189_770 * h10_m2 - e_5 * fs_147_2717_66 * r_2 * h8_m4 + e_5 * fs_294_209_15 * r_2 * h8_m2 + e_5 * fs_247_187_6 * r_4 * h6_m4 + e_5 * fs_38_17_5 * r_4 * h6_m2 + e_5 * fs_204_143_30 * r_6 * h4_m4 - e_5 * fs_250_143_14 * r_8 * h2_m2 + e_6 * fs_189_46189_10010 * r_2 * h10_m4 + e_6 * fs_378_46189_770 * r_2 * h10_m2 + e_6 * fs_7_2717_66 * r_4 * h8_m4 - e_6 * fs_14_209_15 * r_4 * h8_m2 - e_6 * fs_26_561_6 * r_6 * h6_m4 - e_6 * fs_4_51_5 * r_6 * h6_m2 - e_6 * fs_6_143_30 * r_8 * h4_m4 + e_6 * fs_20_429_14 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_43[k] = e_2 * fs_765_8_30 * h4_m3 + e_3 * fs_855_44_66 * h6_m5 - e_3 * fs_285_44_10 * h6_m3 - e_3 * fs_765_11_30 * r_2 * h4_m3 - e_4 * fs_147_143_429 * h8_m5 - e_4 * fs_588_143_55 * h8_m3 - e_4 * fs_171_22_66 * r_2 * h6_m5 + e_4 * fs_57_22_10 * r_2 * h6_m3 + e_4 * fs_2295_143_30 * r_4 * h4_m3 - e_5 * fs_7245_92378_2002 * h10_m5 - e_5 * fs_1449_92378_10010 * h10_m3 + e_5 * fs_588_2717_429 * r_2 * h8_m5 + e_5 * fs_2352_2717_55 * r_2 * h8_m3 + e_5 * fs_171_187_66 * r_4 * h6_m5 - e_5 * fs_57_187_10 * r_4 * h6_m3 - e_5 * fs_204_143_30 * r_6 * h4_m3 + e_6 * fs_315_46189_2002 * r_2 * h10_m5 + e_6 * fs_63_46189_10010 * r_2 * h10_m3 - e_6 * fs_28_2717_429 * r_4 * h8_m5 - e_6 * fs_112_2717_55 * r_4 * h8_m3 - e_6 * fs_6_187_66 * r_6 * h6_m5 + e_6 * fs_2_187_10 * r_6 * h6_m3 + e_6 * fs_6_143_30 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];

        pc_44[k] = - e_2 * fs_765_8_6 * h4_m4 - e_3 * fs_285_22_55 * h6_m6 - e_3 * fs_1425_44_30 * h6_m4 + e_3 * fs_765_11_6 * r_2 * h4_m4 - e_4 * fs_147_286_5005 * h8_m6 - e_4 * fs_735_572_330 * h8_m4 + e_4 * fs_57_11_55 * r_2 * h6_m6 + e_4 * fs_285_22_30 * r_2 * h6_m4 - e_4 * fs_2295_143_6 * r_4 * h4_m4 - e_5 * fs_2898_46189_1001 * h10_m6 - e_5 * fs_1449_92378_2002 * h10_m4 + e_5 * fs_294_2717_5005 * r_2 * h8_m6 + e_5 * fs_735_2717_330 * r_2 * h8_m4 - e_5 * fs_114_187_55 * r_4 * h6_m6 - e_5 * fs_285_187_30 * r_4 * h6_m4 + e_5 * fs_204_143_6 * r_6 * h4_m4 + e_6 * fs_252_46189_1001 * r_2 * h10_m6 + e_6 * fs_63_46189_2002 * r_2 * h10_m4 - e_6 * fs_14_2717_5005 * r_4 * h8_m6 - e_6 * fs_35_2717_330 * r_4 * h8_m4 + e_6 * fs_4_187_55 * r_6 * h6_m6 + e_6 * fs_10_187_30 * r_6 * h6_m4 - e_6 * fs_6_143_6 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph8_0, ph8_p1, ph10_0, ph10_p1, ab_2, pc_45, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p1 = ph10_p1[k];

        pc_45[k] = e_0 * f_12285_32 - e_1 * f_7875_8 * h2_0 - e_1 * f_12285_8 * r_2 + e_2 * f_2295_4 * h4_0 + e_2 * f_5625_4 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 - e_3 * f_1900_11 * h6_0 - e_3 * f_4590_11 * r_2 * h4_0 - e_3 * f_625_1 * r_4 * h2_0 - e_3 * f_585_1 * r_6 + e_4 * f_5145_143 * h8_0 + e_4 * f_760_11 * r_2 * h6_0 + e_4 * f_13770_143 * r_4 * h4_0 + e_4 * f_1250_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 - e_5 * f_365148_46189 * h10_0 - e_5 * f_20580_2717 * r_2 * h8_0 - e_5 * f_1520_187 * r_4 * h6_0 - e_5 * f_1224_143 * r_6 * h4_0 - e_5 * f_1250_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 + e_6 * f_31752_46189 * r_2 * h10_0 + e_6 * f_980_2717 * r_4 * h8_0 + e_6 * f_160_561 * r_6 * h6_0 + e_6 * f_36_143 * r_8 * h4_0 + e_6 * f_100_429 * r_10 * h2_0 + e_6 * f_2_11 * r_12;

        pc_46[k] = - e_1 * fs_1575_16_5 * h2_p1 + e_2 * fs_765_8_6 * h4_p1 + e_2 * fs_1125_8_5 * r_2 * h2_p1 - e_3 * fs_190_11_35 * h6_p1 - e_3 * fs_765_11_6 * r_2 * h4_p1 - e_3 * fs_125_2_5 * r_4 * h2_p1 + e_4 * fs_1029_143_15 * h8_p1 + e_4 * fs_76_11_35 * r_2 * h6_p1 + e_4 * fs_2295_143_6 * r_4 * h4_p1 + e_4 * fs_125_11_5 * r_6 * h2_p1 - e_5 * fs_60858_46189_33 * h10_p1 - e_5 * fs_4116_2717_15 * r_2 * h8_p1 - e_5 * fs_152_187_35 * r_4 * h6_p1 - e_5 * fs_204_143_6 * r_6 * h4_p1 - e_5 * fs_125_143_5 * r_8 * h2_p1 + e_6 * fs_5292_46189_33 * r_2 * h10_p1 + e_6 * fs_196_2717_15 * r_4 * h8_p1 + e_6 * fs_16_561_35 * r_6 * h6_p1 + e_6 * fs_6_143_6 * r_8 * h4_p1 + e_6 * fs_10_429_5 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ab_2, pc_47, pc_48 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p3 = ph10_p3[k];

        pc_47[k] = e_1 * fs_1575_8_35 * h2_p2 - e_2 * fs_765_8_21 * h4_p2 - e_2 * fs_1125_4_35 * r_2 * h2_p2 + e_3 * fs_475_11_2 * h6_p2 + e_3 * fs_765_11_21 * r_2 * h4_p2 + e_3 * fs_125_1_35 * r_4 * h2_p2 + e_4 * fs_735_286_6 * h8_p2 - e_4 * fs_190_11_2 * r_2 * h6_p2 - e_4 * fs_2295_143_21 * r_4 * h4_p2 - e_4 * fs_250_11_35 * r_6 * h2_p2 - e_5 * fs_34776_46189_77 * h10_p2 - e_5 * fs_1470_2717_6 * r_2 * h8_p2 + e_5 * fs_380_187_2 * r_4 * h6_p2 + e_5 * fs_204_143_21 * r_6 * h4_p2 + e_5 * fs_250_143_35 * r_8 * h2_p2 + e_6 * fs_3024_46189_77 * r_2 * h10_p2 + e_6 * fs_70_2717_6 * r_4 * h8_p2 - e_6 * fs_40_561_2 * r_6 * h6_p2 - e_6 * fs_6_143_21 * r_8 * h4_p2 - e_6 * fs_20_429_35 * r_10 * h2_p2;

        pc_48[k] = - e_2 * f_2295_4 * h4_p3 + e_3 * fs_2375_22_3 * h6_p3 + e_3 * f_4590_11 * r_2 * h4_p3 - e_4 * fs_735_286_66 * h8_p3 - e_4 * fs_475_11_3 * r_2 * h6_p3 - e_4 * f_13770_143 * r_4 * h4_p3 - e_5 * fs_4347_46189_3003 * h10_p3 + e_5 * fs_1470_2717_66 * r_2 * h8_p3 + e_5 * fs_950_187_3 * r_4 * h6_p3 + e_5 * f_1224_143 * r_6 * h4_p3 + e_6 * fs_378_46189_3003 * r_2 * h10_p3 - e_6 * fs_70_2717_66 * r_4 * h8_p3 - e_6 * fs_100_561_3 * r_6 * h6_p3 - e_6 * f_36_143 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p4, ph6_p4, ph6_p5, ph8_p4, ph8_p5, ph10_p4, ph10_p5, ab_2, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];

        pc_49[k] = e_2 * f_2295_4 * h4_p4 + e_3 * fs_570_11_5 * h6_p4 - e_3 * f_4590_11 * r_2 * h4_p4 - e_4 * fs_147_26_55 * h8_p4 - e_4 * fs_228_11_5 * r_2 * h6_p4 + e_4 * f_13770_143 * r_4 * h4_p4 - e_5 * fs_2898_46189_3003 * h10_p4 + e_5 * fs_294_247_55 * r_2 * h8_p4 + e_5 * fs_456_187_5 * r_4 * h6_p4 - e_5 * f_1224_143 * r_6 * h4_p4 + e_6 * fs_252_46189_3003 * r_2 * h10_p4 - e_6 * fs_14_247_55 * r_4 * h8_p4 - e_6 * fs_16_187_5 * r_6 * h6_p4 + e_6 * f_36_143 * r_8 * h4_p4;

        pc_50[k] = - e_3 * fs_1425_22_11 * h6_p5 - e_4 * fs_735_286_286 * h8_p5 + e_4 * fs_285_11_11 * r_2 * h6_p5 - e_5 * fs_1449_46189_3003 * h10_p5 + e_5 * fs_1470_2717_286 * r_2 * h8_p5 - e_5 * fs_570_187_11 * r_4 * h6_p5 + e_6 * fs_126_46189_3003 * r_2 * h10_p5 - e_6 * fs_70_2717_286 * r_4 * h8_p5 + e_6 * fs_20_187_11 * r_6 * h6_p5;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p2 = ph10_p2[k];

        pc_51[k] = e_0 * f_12285_32 - e_1 * f_14175_16 * h2_0 - e_1 * fs_7875_16_3 * h2_p2 - e_1 * f_12285_8 * r_2 + e_2 * f_765_2 * h4_0 + e_2 * fs_765_4_5 * h4_p2 + e_2 * f_10125_8 * r_2 * h2_0 + e_2 * fs_5625_8_3 * r_2 * h2_p2 + e_2 * f_12285_8 * r_4 - e_3 * f_570_11 * h6_0 - e_3 * fs_95_11_210 * h6_p2 - e_3 * f_3060_11 * r_2 * h4_0 - e_3 * fs_1530_11_5 * r_2 * h4_p2 - e_3 * f_1125_2 * r_4 * h2_0 - e_3 * fs_625_2_3 * r_4 * h2_p2 - e_3 * f_585_1 * r_6 - e_4 * f_1029_143 * h8_0 + e_4 * fs_441_143_70 * h8_p2 + e_4 * f_228_11 * r_2 * h6_0 + e_4 * fs_38_11_210 * r_2 * h6_p2 + e_4 * f_9180_143 * r_4 * h4_0 + e_4 * fs_4590_143_5 * r_4 * h4_p2 + e_4 * f_1125_11 * r_6 * h2_0 + e_4 * fs_625_11_3 * r_6 * h2_p2 + e_4 * f_195_2 * r_8 + e_5 * f_304290_46189 * h10_0 - e_5 * fs_20286_46189_165 * h10_p2 + e_5 * f_4116_2717 * r_2 * h8_0 - e_5 * fs_1764_2717_70 * r_2 * h8_p2 - e_5 * f_456_187 * r_4 * h6_0 - e_5 * fs_76_187_210 * r_4 * h6_p2 - e_5 * f_816_143 * r_6 * h4_0 - e_5 * fs_408_143_5 * r_6 * h4_p2 - e_5 * f_1125_143 * r_8 * h2_0 - e_5 * fs_625_143_3 * r_8 * h2_p2 - e_5 * f_78_11 * r_10 - e_6 * f_26460_46189 * r_2 * h10_0 + e_6 * fs_1764_46189_165 * r_2 * h10_p2 - e_6 * f_196_2717 * r_4 * h8_0 + e_6 * fs_84_2717_70 * r_4 * h8_p2 + e_6 * f_16_187 * r_6 * h6_0 + e_6 * fs_8_561_210 * r_6 * h6_p2 + e_6 * f_24_143 * r_8 * h4_0 + e_6 * fs_12_143_5 * r_8 * h4_p2 + e_6 * f_30_143 * r_10 * h2_0 + e_6 * fs_50_429_3 * r_10 * h2_p2 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_52[k] = - e_1 * fs_1575_16_21 * h2_p1 + e_2 * fs_765_16_70 * h4_p1 + e_2 * fs_765_16_10 * h4_p3 + e_2 * fs_1125_8_21 * r_2 * h2_p1 - e_3 * fs_760_11_3 * h6_p1 - e_3 * fs_285_22_30 * h6_p3 - e_3 * fs_765_22_70 * r_2 * h4_p1 - e_3 * fs_765_22_10 * r_2 * h4_p3 - e_3 * fs_125_2_21 * r_4 * h2_p1 + e_4 * fs_1323_286_7 * h8_p1 + e_4 * fs_441_286_165 * h8_p3 + e_4 * fs_304_11_3 * r_2 * h6_p1 + e_4 * fs_57_11_30 * r_2 * h6_p3 + e_4 * fs_2295_286_70 * r_4 * h4_p1 + e_4 * fs_2295_286_10 * r_4 * h4_p3 + e_4 * fs_125_11_21 * r_6 * h2_p1 + e_5 * fs_8694_46189_385 * h10_p1 - e_5 * fs_1449_46189_30030 * h10_p3 - e_5 * fs_2646_2717_7 * r_2 * h8_p1 - e_5 * fs_882_2717_165 * r_2 * h8_p3 - e_5 * fs_608_187_3 * r_4 * h6_p1 - e_5 * fs_114_187_30 * r_4 * h6_p3 - e_5 * fs_102_143_70 * r_6 * h4_p1 - e_5 * fs_102_143_10 * r_6 * h4_p3 - e_5 * fs_125_143_21 * r_8 * h2_p1 - e_6 * fs_756_46189_385 * r_2 * h10_p1 + e_6 * fs_126_46189_30030 * r_2 * h10_p3 + e_6 * fs_126_2717_7 * r_4 * h8_p1 + e_6 * fs_42_2717_165 * r_4 * h8_p3 + e_6 * fs_64_561_3 * r_6 * h6_p1 + e_6 * fs_4_187_30 * r_6 * h6_p3 + e_6 * fs_3_143_70 * r_8 * h4_p1 + e_6 * fs_3_143_10 * r_8 * h4_p3 + e_6 * fs_10_429_21 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_53[k] = e_1 * fs_1575_8_14 * h2_p2 - e_2 * fs_765_8_30 * h4_p4 - e_2 * fs_1125_4_14 * r_2 * h2_p2 - e_3 * fs_95_2_5 * h6_p2 + e_3 * fs_1235_44_6 * h6_p4 + e_3 * fs_765_11_30 * r_2 * h4_p4 + e_3 * fs_125_1_14 * r_4 * h2_p2 + e_4 * fs_147_22_15 * h8_p2 + e_4 * fs_147_572_66 * h8_p4 + e_4 * fs_19_1_5 * r_2 * h6_p2 - e_4 * fs_247_22_6 * r_2 * h6_p4 - e_4 * fs_2295_143_30 * r_4 * h4_p4 - e_4 * fs_250_11_14 * r_6 * h2_p2 + e_5 * fs_4347_46189_770 * h10_p2 - e_5 * fs_4347_92378_10010 * h10_p4 - e_5 * fs_294_209_15 * r_2 * h8_p2 - e_5 * fs_147_2717_66 * r_2 * h8_p4 - e_5 * fs_38_17_5 * r_4 * h6_p2 + e_5 * fs_247_187_6 * r_4 * h6_p4 + e_5 * fs_204_143_30 * r_6 * h4_p4 + e_5 * fs_250_143_14 * r_8 * h2_p2 - e_6 * fs_378_46189_770 * r_2 * h10_p2 + e_6 * fs_189_46189_10010 * r_2 * h10_p4 + e_6 * fs_14_209_15 * r_4 * h8_p2 + e_6 * fs_7_2717_66 * r_4 * h8_p4 + e_6 * fs_4_51_5 * r_6 * h6_p2 - e_6 * fs_26_561_6 * r_6 * h6_p4 - e_6 * fs_6_143_30 * r_8 * h4_p4 - e_6 * fs_20_429_14 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_54[k] = - e_2 * fs_765_8_30 * h4_p3 + e_3 * fs_285_44_10 * h6_p3 + e_3 * fs_855_44_66 * h6_p5 + e_3 * fs_765_11_30 * r_2 * h4_p3 + e_4 * fs_588_143_55 * h8_p3 - e_4 * fs_147_143_429 * h8_p5 - e_4 * fs_57_22_10 * r_2 * h6_p3 - e_4 * fs_171_22_66 * r_2 * h6_p5 - e_4 * fs_2295_143_30 * r_4 * h4_p3 + e_5 * fs_1449_92378_10010 * h10_p3 - e_5 * fs_7245_92378_2002 * h10_p5 - e_5 * fs_2352_2717_55 * r_2 * h8_p3 + e_5 * fs_588_2717_429 * r_2 * h8_p5 + e_5 * fs_57_187_10 * r_4 * h6_p3 + e_5 * fs_171_187_66 * r_4 * h6_p5 + e_5 * fs_204_143_30 * r_6 * h4_p3 - e_6 * fs_63_46189_10010 * r_2 * h10_p3 + e_6 * fs_315_46189_2002 * r_2 * h10_p5 + e_6 * fs_112_2717_55 * r_4 * h8_p3 - e_6 * fs_28_2717_429 * r_4 * h8_p5 - e_6 * fs_2_187_10 * r_6 * h6_p3 - e_6 * fs_6_187_66 * r_6 * h6_p5 - e_6 * fs_6_143_30 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p4, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p6, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];

        pc_55[k] = e_2 * fs_765_8_6 * h4_p4 + e_3 * fs_1425_44_30 * h6_p4 - e_3 * fs_285_22_55 * h6_p6 - e_3 * fs_765_11_6 * r_2 * h4_p4 + e_4 * fs_735_572_330 * h8_p4 - e_4 * fs_147_286_5005 * h8_p6 - e_4 * fs_285_22_30 * r_2 * h6_p4 + e_4 * fs_57_11_55 * r_2 * h6_p6 + e_4 * fs_2295_143_6 * r_4 * h4_p4 + e_5 * fs_1449_92378_2002 * h10_p4 - e_5 * fs_2898_46189_1001 * h10_p6 - e_5 * fs_735_2717_330 * r_2 * h8_p4 + e_5 * fs_294_2717_5005 * r_2 * h8_p6 + e_5 * fs_285_187_30 * r_4 * h6_p4 - e_5 * fs_114_187_55 * r_4 * h6_p6 - e_5 * fs_204_143_6 * r_6 * h4_p4 - e_6 * fs_63_46189_2002 * r_2 * h10_p4 + e_6 * fs_252_46189_1001 * r_2 * h10_p6 + e_6 * fs_35_2717_330 * r_4 * h8_p4 - e_6 * fs_14_2717_5005 * r_4 * h8_p6 - e_6 * fs_10_187_30 * r_6 * h6_p4 + e_6 * fs_4_187_55 * r_6 * h6_p6 + e_6 * fs_6_143_6 * r_8 * h4_p4;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_56[k] = e_0 * f_12285_32 - e_1 * f_4725_8 * h2_0 - e_1 * f_12285_8 * r_2 - e_2 * f_765_8 * h4_0 + e_2 * fs_765_8_35 * h4_p4 + e_2 * f_3375_4 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 + e_3 * f_1710_11 * h6_0 - e_3 * fs_570_11_7 * h6_p4 + e_3 * f_765_11 * r_2 * h4_0 - e_3 * fs_765_11_35 * r_2 * h4_p4 - e_3 * f_375_1 * r_4 * h2_0 - e_3 * f_585_1 * r_6 - e_4 * f_4998_143 * h8_0 + e_4 * fs_441_143_77 * h8_p4 - e_4 * f_684_11 * r_2 * h6_0 + e_4 * fs_228_11_7 * r_2 * h6_p4 - e_4 * f_2295_143 * r_4 * h4_0 + e_4 * fs_2295_143_35 * r_4 * h4_p4 + e_4 * f_750_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 - e_5 * f_173880_46189 * h10_0 - e_5 * fs_5796_46189_2145 * h10_p4 + e_5 * f_19992_2717 * r_2 * h8_0 - e_5 * fs_1764_2717_77 * r_2 * h8_p4 + e_5 * f_1368_187 * r_4 * h6_0 - e_5 * fs_456_187_7 * r_4 * h6_p4 + e_5 * f_204_143 * r_6 * h4_0 - e_5 * fs_204_143_35 * r_6 * h4_p4 - e_5 * f_750_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 + e_6 * f_15120_46189 * r_2 * h10_0 + e_6 * fs_504_46189_2145 * r_2 * h10_p4 - e_6 * f_952_2717 * r_4 * h8_0 + e_6 * fs_84_2717_77 * r_4 * h8_p4 - e_6 * f_48_187 * r_6 * h6_0 + e_6 * fs_16_187_7 * r_6 * h6_p4 - e_6 * f_6_143 * r_8 * h4_0 + e_6 * fs_6_143_35 * r_8 * h4_p4 + e_6 * f_20_143 * r_10 * h2_0 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_57[k] = - e_1 * fs_7875_16_2 * h2_p1 + e_2 * fs_765_8_15 * h4_p1 + e_2 * fs_5625_8_2 * r_2 * h2_p1 + e_3 * fs_95_22_14 * h6_p1 - e_3 * fs_95_22_231 * h6_p5 - e_3 * fs_765_11_15 * r_2 * h4_p1 - e_3 * fs_625_2_2 * r_4 * h2_p1 - e_4 * fs_6909_572_6 * h8_p1 + e_4 * fs_147_572_6006 * h8_p5 - e_4 * fs_19_11_14 * r_2 * h6_p1 + e_4 * fs_19_11_231 * r_2 * h6_p5 + e_4 * fs_2295_143_15 * r_4 * h4_p1 + e_4 * fs_625_11_2 * r_6 * h2_p1 - e_5 * fs_4347_46189_330 * h10_p1 - e_5 * fs_21735_46189_143 * h10_p5 + e_5 * fs_6909_2717_6 * r_2 * h8_p1 - e_5 * fs_147_2717_6006 * r_2 * h8_p5 + e_5 * fs_38_187_14 * r_4 * h6_p1 - e_5 * fs_38_187_231 * r_4 * h6_p5 - e_5 * fs_204_143_15 * r_6 * h4_p1 - e_5 * fs_625_143_2 * r_8 * h2_p1 + e_6 * fs_378_46189_330 * r_2 * h10_p1 + e_6 * fs_1890_46189_143 * r_2 * h10_p5 - e_6 * fs_329_2717_6 * r_4 * h8_p1 + e_6 * fs_7_2717_6006 * r_4 * h8_p5 - e_6 * fs_4_561_14 * r_6 * h6_p1 + e_6 * fs_4_561_231 * r_6 * h6_p5 + e_6 * fs_6_143_15 * r_8 * h4_p1 + e_6 * fs_50_429_2 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_58[k] = e_1 * f_4725_8 * h2_p2 + e_2 * fs_765_8_15 * h4_p2 - e_2 * f_3375_4 * r_2 * h2_p2 - e_3 * fs_285_22_70 * h6_p2 + e_3 * fs_285_22_154 * h6_p6 - e_3 * fs_765_11_15 * r_2 * h4_p2 + e_3 * f_375_1 * r_4 * h2_p2 - e_4 * fs_1029_572_210 * h8_p2 - e_4 * fs_147_572_286 * h8_p6 + e_4 * fs_57_11_70 * r_2 * h6_p2 - e_4 * fs_57_11_154 * r_2 * h6_p6 + e_4 * fs_2295_143_15 * r_4 * h4_p2 - e_4 * f_750_11 * r_6 * h2_p2 - e_5 * fs_5796_46189_55 * h10_p2 - e_5 * fs_5796_46189_1430 * h10_p6 + e_5 * fs_1029_2717_210 * r_2 * h8_p2 + e_5 * fs_147_2717_286 * r_2 * h8_p6 - e_5 * fs_114_187_70 * r_4 * h6_p2 + e_5 * fs_114_187_154 * r_4 * h6_p6 - e_5 * fs_204_143_15 * r_6 * h4_p2 + e_5 * f_750_143 * r_8 * h2_p2 + e_6 * fs_504_46189_55 * r_2 * h10_p2 + e_6 * fs_504_46189_1430 * r_2 * h10_p6 - e_6 * fs_49_2717_210 * r_4 * h8_p2 - e_6 * fs_7_2717_286 * r_4 * h8_p6 + e_6 * fs_4_187_70 * r_6 * h6_p2 - e_6 * fs_4_187_154 * r_6 * h6_p6 + e_6 * fs_6_143_15 * r_8 * h4_p2 - e_6 * f_20_143 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ab_2, pc_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];

        pc_59[k] = - e_2 * fs_765_8_21 * h4_p3 - e_3 * fs_1425_22_7 * h6_p3 + e_3 * fs_765_11_21 * r_2 * h4_p3 - e_4 * fs_735_572_154 * h8_p3 - e_4 * fs_735_572_858 * h8_p7 + e_4 * fs_285_11_7 * r_2 * h6_p3 - e_4 * fs_2295_143_21 * r_4 * h4_p3 - e_5 * fs_1449_46189_143 * h10_p3 - e_5 * fs_2898_46189_2431 * h10_p7 + e_5 * fs_735_2717_154 * r_2 * h8_p3 + e_5 * fs_735_2717_858 * r_2 * h8_p7 - e_5 * fs_570_187_7 * r_4 * h6_p3 + e_5 * fs_204_143_21 * r_6 * h4_p3 + e_6 * fs_126_46189_143 * r_2 * h10_p3 + e_6 * fs_252_46189_2431 * r_2 * h10_p7 - e_6 * fs_35_2717_154 * r_4 * h8_p3 - e_6 * fs_35_2717_858 * r_4 * h8_p7 + e_6 * fs_20_187_7 * r_6 * h6_p3 - e_6 * fs_6_143_21 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_60[k] = e_0 * f_12285_32 - e_1 * f_1575_16 * h2_0 - e_1 * f_12285_8 * r_2 - e_2 * f_2295_4 * h4_0 + e_2 * f_1125_8 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 + e_3 * f_2755_22 * h6_0 - e_3 * fs_95_11_462 * h6_p6 + e_3 * f_4590_11 * r_2 * h4_0 - e_3 * f_125_2 * r_4 * h2_0 - e_3 * f_585_1 * r_6 + e_4 * f_10731_286 * h8_0 + e_4 * fs_147_143_858 * h8_p6 - e_4 * f_551_11 * r_2 * h6_0 + e_4 * fs_38_11_462 * r_2 * h6_p6 - e_4 * f_13770_143 * r_4 * h4_0 + e_4 * f_125_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 + e_5 * f_65205_46189 * h10_0 - e_5 * fs_4347_46189_4290 * h10_p6 - e_5 * f_21462_2717 * r_2 * h8_0 - e_5 * fs_588_2717_858 * r_2 * h8_p6 + e_5 * f_1102_187 * r_4 * h6_0 - e_5 * fs_76_187_462 * r_4 * h6_p6 + e_5 * f_1224_143 * r_6 * h4_0 - e_5 * f_125_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 - e_6 * f_5670_46189 * r_2 * h10_0 + e_6 * fs_378_46189_4290 * r_2 * h10_p6 + e_6 * f_1022_2717 * r_4 * h8_0 + e_6 * fs_28_2717_858 * r_4 * h8_p6 - e_6 * f_116_561 * r_6 * h6_0 + e_6 * fs_8_561_462 * r_6 * h6_p6 - e_6 * f_36_143 * r_8 * h4_0 + e_6 * f_10_429 * r_10 * h2_0 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_61[k] = - e_1 * fs_11025_32_6 * h2_p1 + e_2 * fs_7875_16_6 * r_2 * h2_p1 + e_3 * fs_95_4_42 * h6_p1 - e_3 * fs_875_4_6 * r_4 * h2_p1 + e_4 * fs_147_11_2 * h8_p1 + e_4 * fs_147_286_1430 * h8_p7 - e_4 * fs_19_2_42 * r_2 * h6_p1 + e_4 * fs_875_22_6 * r_6 * h2_p1 + e_5 * fs_4347_92378_110 * h10_p1 - e_5 * fs_1449_46189_36465 * h10_p7 - e_5 * fs_588_209_2 * r_2 * h8_p1 - e_5 * fs_294_2717_1430 * r_2 * h8_p7 + e_5 * fs_19_17_42 * r_4 * h6_p1 - e_5 * fs_875_286_6 * r_8 * h2_p1 - e_6 * fs_189_46189_110 * r_2 * h10_p1 + e_6 * fs_126_46189_36465 * r_2 * h10_p7 + e_6 * fs_28_209_2 * r_4 * h8_p1 + e_6 * fs_14_2717_1430 * r_4 * h8_p7 - e_6 * fs_2_51_42 * r_6 * h6_p1 + e_6 * fs_35_429_6 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_62[k] = e_1 * fs_1575_16_15 * h2_p2 + e_2 * f_2295_4 * h4_p2 - e_2 * fs_1125_8_15 * r_2 * h2_p2 + e_3 * fs_475_22_42 * h6_p2 - e_3 * f_4590_11 * r_2 * h4_p2 + e_3 * fs_125_2_15 * r_4 * h2_p2 + e_4 * fs_735_286_14 * h8_p2 - e_4 * fs_735_286_143 * h8_p8 - e_4 * fs_95_11_42 * r_2 * h6_p2 + e_4 * f_13770_143 * r_4 * h4_p2 - e_4 * fs_125_11_15 * r_6 * h2_p2 + e_5 * fs_1449_46189_33 * h10_p2 - e_5 * fs_4347_46189_2431 * h10_p8 - e_5 * fs_1470_2717_14 * r_2 * h8_p2 + e_5 * fs_1470_2717_143 * r_2 * h8_p8 + e_5 * fs_190_187_42 * r_4 * h6_p2 - e_5 * f_1224_143 * r_6 * h4_p2 + e_5 * fs_125_143_15 * r_8 * h2_p2 - e_6 * fs_126_46189_33 * r_2 * h10_p2 + e_6 * fs_378_46189_2431 * r_2 * h10_p8 + e_6 * fs_70_2717_14 * r_4 * h8_p2 - e_6 * fs_70_2717_143 * r_4 * h8_p8 - e_6 * fs_20_561_42 * r_6 * h6_p2 + e_6 * f_36_143 * r_8 * h4_p2 - e_6 * fs_10_429_15 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_63[k] = e_0 * f_12285_32 + e_1 * f_4725_8 * h2_0 - e_1 * f_12285_8 * r_2 - e_2 * f_2295_4 * h4_0 - e_2 * f_3375_4 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 - e_3 * f_2280_11 * h6_0 + e_3 * f_4590_11 * r_2 * h4_0 + e_3 * f_375_1 * r_4 * h2_0 - e_3 * f_585_1 * r_6 - e_4 * f_4557_286 * h8_0 + e_4 * fs_441_286_715 * h8_p8 + e_4 * f_912_11 * r_2 * h6_0 - e_4 * f_13770_143 * r_4 * h4_0 - e_4 * f_750_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 - e_5 * f_14490_46189 * h10_0 - e_5 * fs_2898_46189_12155 * h10_p8 + e_5 * f_9114_2717 * r_2 * h8_0 - e_5 * fs_882_2717_715 * r_2 * h8_p8 - e_5 * f_1824_187 * r_4 * h6_0 + e_5 * f_1224_143 * r_6 * h4_0 + e_5 * f_750_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 + e_6 * f_1260_46189 * r_2 * h10_0 + e_6 * fs_252_46189_12155 * r_2 * h10_p8 - e_6 * f_434_2717 * r_4 * h8_0 + e_6 * fs_42_2717_715 * r_4 * h8_p8 + e_6 * f_64_187 * r_6 * h6_0 - e_6 * f_36_143 * r_8 * h4_0 - e_6 * f_20_143 * r_10 * h2_0 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ab_2, pc_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_64[k] = - e_1 * fs_4725_32_30 * h2_p1 - e_2 * f_2295_4 * h4_p1 + e_2 * fs_3375_16_30 * r_2 * h2_p1 - e_3 * fs_285_44_210 * h6_p1 + e_3 * f_4590_11 * r_2 * h4_p1 - e_3 * fs_375_4_30 * r_4 * h2_p1 - e_4 * fs_441_286_10 * h8_p1 + e_4 * fs_57_22_210 * r_2 * h6_p1 - e_4 * f_13770_143 * r_4 * h4_p1 + e_4 * fs_375_22_30 * r_6 * h2_p1 - e_5 * fs_1449_92378_22 * h10_p1 - e_5 * fs_1449_46189_46189 * h10_p9 + e_5 * fs_882_2717_10 * r_2 * h8_p1 - e_5 * fs_57_187_210 * r_4 * h6_p1 + e_5 * f_1224_143 * r_6 * h4_p1 - e_5 * fs_375_286_30 * r_8 * h2_p1 + e_6 * fs_63_46189_22 * r_2 * h10_p1 + e_6 * fs_126_46189_46189 * r_2 * h10_p9 - e_6 * fs_42_2717_10 * r_4 * h8_p1 + e_6 * fs_2_187_210 * r_6 * h6_p1 - e_6 * f_36_143 * r_8 * h4_p1 + e_6 * fs_5_143_30 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 56 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ab_2, pc_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];

        pc_65[k] = e_0 * f_12285_32 + e_1 * f_23625_16 * h2_0 - e_1 * f_12285_8 * r_2 + e_2 * f_2295_4 * h4_0 - e_2 * f_16875_8 * r_2 * h2_0 + e_2 * f_12285_8 * r_4 + e_3 * f_1425_22 * h6_0 - e_3 * f_4590_11 * r_2 * h4_0 + e_3 * f_1875_2 * r_4 * h2_0 - e_3 * f_585_1 * r_6 + e_4 * f_735_286 * h8_0 - e_4 * f_285_11 * r_2 * h6_0 + e_4 * f_13770_143 * r_4 * h4_0 - e_4 * f_1875_11 * r_6 * h2_0 + e_4 * f_195_2 * r_8 + e_5 * f_1449_46189 * h10_0 - e_5 * fs_1449_46189_92378 * h10_p10 - e_5 * f_1470_2717 * r_2 * h8_0 + e_5 * f_570_187 * r_4 * h6_0 - e_5 * f_1224_143 * r_6 * h4_0 + e_5 * f_1875_143 * r_8 * h2_0 - e_5 * f_78_11 * r_10 - e_6 * f_126_46189 * r_2 * h10_0 + e_6 * fs_126_46189_92378 * r_2 * h10_p10 + e_6 * f_70_2717 * r_4 * h8_0 - e_6 * f_20_187 * r_6 * h6_0 + e_6 * f_36_143 * r_8 * h4_0 - e_6 * f_50_143 * r_10 * h2_0 + e_6 * f_2_11 * r_12;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[121] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 2, 13, 24, 25, 26, 27, 28, 29, 30, 31, 32, 3, 14, 25, 36, 37, 38, 39, 40, 41, 42, 43, 4, 15, 26, 37, 48, 49, 50, 51, 52, 53, 54, 5, 16, 27, 38, 49, 60, 61, 62, 63, 64, 65, 6, 17, 28, 39, 50, 61, 72, 73, 74, 75, 76, 7, 18, 29, 40, 51, 62, 73, 84, 85, 86, 87, 8, 19, 30, 41, 52, 63, 74, 85, 96, 97, 98, 9, 20, 31, 42, 53, 64, 75, 86, 97, 108, 109, 10, 21, 32, 43, 54, 65, 76, 87, 98, 109, 120};

    for (size_t n = 0; n < 121; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
