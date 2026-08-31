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



#include "SimdKineticEnergyRecGG.hpp"

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
compute_gg_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gg_kinetic_energy: Basis functions must be of angular momenta four and four"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gg_kinetic_energy: Harmonics must reach angular momentum 8"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gg_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 81 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time.

    auto buffer = CSimdMatrix(6, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);
    std::fill(pe_4, pe_4 + nmax, 0.0);
    std::fill(pe_5, pe_5 + nmax, 0.0);

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

            const auto ff_0 = fbase * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fterm = std::exp(-fmu * ab_2[k]);

                pe_0[k] += ff_0 * fterm;
                pe_1[k] += ff_1 * fterm;
                pe_2[k] += ff_2 * fterm;
                pe_3[k] += ff_3 * fterm;
                pe_4[k] += ff_4 * fterm;
                pe_5[k] += ff_5 * fterm;
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

    auto *pc_0 = values + 0 * nvalues;
    auto *pc_1 = values + 1 * nvalues;
    auto *pc_2 = values + 2 * nvalues;
    auto *pc_3 = values + 3 * nvalues;
    auto *pc_4 = values + 4 * nvalues;
    auto *pc_5 = values + 5 * nvalues;
    auto *pc_6 = values + 6 * nvalues;
    auto *pc_7 = values + 7 * nvalues;
    auto *pc_8 = values + 8 * nvalues;
    auto *pc_9 = values + 10 * nvalues;
    auto *pc_10 = values + 11 * nvalues;
    auto *pc_11 = values + 12 * nvalues;
    auto *pc_12 = values + 13 * nvalues;
    auto *pc_13 = values + 14 * nvalues;
    auto *pc_14 = values + 15 * nvalues;
    auto *pc_15 = values + 16 * nvalues;
    auto *pc_16 = values + 17 * nvalues;
    auto *pc_17 = values + 20 * nvalues;
    auto *pc_18 = values + 21 * nvalues;
    auto *pc_19 = values + 22 * nvalues;
    auto *pc_20 = values + 23 * nvalues;
    auto *pc_21 = values + 24 * nvalues;
    auto *pc_22 = values + 25 * nvalues;
    auto *pc_23 = values + 26 * nvalues;
    auto *pc_24 = values + 30 * nvalues;
    auto *pc_25 = values + 31 * nvalues;
    auto *pc_26 = values + 32 * nvalues;
    auto *pc_27 = values + 33 * nvalues;
    auto *pc_28 = values + 34 * nvalues;
    auto *pc_29 = values + 35 * nvalues;
    auto *pc_30 = values + 40 * nvalues;
    auto *pc_31 = values + 41 * nvalues;
    auto *pc_32 = values + 42 * nvalues;
    auto *pc_33 = values + 43 * nvalues;
    auto *pc_34 = values + 44 * nvalues;
    auto *pc_35 = values + 50 * nvalues;
    auto *pc_36 = values + 51 * nvalues;
    auto *pc_37 = values + 52 * nvalues;
    auto *pc_38 = values + 53 * nvalues;
    auto *pc_39 = values + 60 * nvalues;
    auto *pc_40 = values + 61 * nvalues;
    auto *pc_41 = values + 62 * nvalues;
    auto *pc_42 = values + 70 * nvalues;
    auto *pc_43 = values + 71 * nvalues;
    auto *pc_44 = values + 80 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_1040_99 = 1040.0 / 99.0;
    const auto f_1064_1287 = 1064.0 / 1287.0;
    const auto f_10_99 = 10.0 / 99.0;
    const auto f_1105_21 = 1105.0 / 21.0;
    const auto f_1105_7 = 1105.0 / 7.0;
    const auto f_1105_8 = 1105.0 / 8.0;
    const auto f_112_1287 = 112.0 / 1287.0;
    const auto f_1155_16 = 1155.0 / 16.0;
    const auto f_1215_14 = 1215.0 / 14.0;
    const auto f_1215_143 = 1215.0 / 143.0;
    const auto f_1215_22 = 1215.0 / 22.0;
    const auto f_1215_28 = 1215.0 / 28.0;
    const auto f_1300_21 = 1300.0 / 21.0;
    const auto f_1300_7 = 1300.0 / 7.0;
    const auto f_133_1287 = 133.0 / 1287.0;
    const auto f_135_2 = 135.0 / 2.0;
    const auto f_136_99 = 136.0 / 99.0;
    const auto f_1445_66 = 1445.0 / 66.0;
    const auto f_1485_28 = 1485.0 / 28.0;
    const auto f_14_1287 = 14.0 / 1287.0;
    const auto f_162_1001 = 162.0 / 1001.0;
    const auto f_170_33 = 170.0 / 33.0;
    const auto f_170_693 = 170.0 / 693.0;
    const auto f_18_91 = 18.0 / 91.0;
    const auto f_1925_8 = 1925.0 / 8.0;
    const auto f_200_693 = 200.0 / 693.0;
    const auto f_2080_693 = 2080.0 / 693.0;
    const auto f_260_1 = 260.0;
    const auto f_260_3 = 260.0 / 3.0;
    const auto f_260_99 = 260.0 / 99.0;
    const auto f_2_33 = 2.0 / 33.0;
    const auto f_2_9 = 2.0 / 9.0;
    const auto f_2_99 = 2.0 / 99.0;
    const auto f_324_1001 = 324.0 / 1001.0;
    const auto f_325_2 = 325.0 / 2.0;
    const auto f_34_33 = 34.0 / 33.0;
    const auto f_34_99 = 34.0 / 99.0;
    const auto f_3645_1001 = 3645.0 / 1001.0;
    const auto f_3645_154 = 3645.0 / 154.0;
    const auto f_3645_77 = 3645.0 / 77.0;
    const auto f_36_143 = 36.0 / 143.0;
    const auto f_3724_1287 = 3724.0 / 1287.0;
    const auto f_385_2 = 385.0 / 2.0;
    const auto f_392_1287 = 392.0 / 1287.0;
    const auto f_405_11 = 405.0 / 11.0;
    const auto f_405_14 = 405.0 / 14.0;
    const auto f_405_4 = 405.0 / 4.0;
    const auto f_405_91 = 405.0 / 91.0;
    const auto f_40_99 = 40.0 / 99.0;
    const auto f_4420_693 = 4420.0 / 693.0;
    const auto f_455_2 = 455.0 / 2.0;
    const auto f_455_8 = 455.0 / 8.0;
    const auto f_4_9 = 4.0 / 9.0;
    const auto f_5200_693 = 5200.0 / 693.0;
    const auto f_520_21 = 520.0 / 21.0;
    const auto f_520_7 = 520.0 / 7.0;
    const auto f_54_143 = 54.0 / 143.0;
    const auto f_55_1 = 55.0;
    const auto f_55_9 = 55.0 / 9.0;
    const auto f_578_99 = 578.0 / 99.0;
    const auto f_65_1 = 65.0;
    const auto f_65_3 = 65.0 / 3.0;
    const auto f_680_99 = 680.0 / 99.0;
    const auto f_68_9 = 68.0 / 9.0;
    const auto f_7290_1001 = 7290.0 / 1001.0;
    const auto f_7448_1287 = 7448.0 / 1287.0;
    const auto f_784_1287 = 784.0 / 1287.0;
    const auto f_80_693 = 80.0 / 693.0;
    const auto f_810_143 = 810.0 / 143.0;
    const auto f_850_33 = 850.0 / 33.0;
    const auto f_85_22 = 85.0 / 22.0;
    const auto f_85_3 = 85.0 / 3.0;
    const auto f_85_66 = 85.0 / 66.0;
    const auto f_8_99 = 8.0 / 99.0;
    const auto f_9310_1287 = 9310.0 / 1287.0;
    const auto f_980_1287 = 980.0 / 1287.0;
    const auto fs_100_693_3 = std::sqrt(10000.0 / 160083.0);
    const auto fs_108_1001_5 = std::sqrt(58320.0 / 1002001.0);
    const auto fs_10_231_21 = std::sqrt(100.0 / 2541.0);
    const auto fs_10_693_30 = std::sqrt(1000.0 / 160083.0);
    const auto fs_10_99_6 = std::sqrt(200.0 / 3267.0);
    const auto fs_1105_132_6 = std::sqrt(1221025.0 / 2904.0);
    const auto fs_1215_1001_35 = std::sqrt(7381125.0 / 143143.0);
    const auto fs_1215_154_35 = std::sqrt(7381125.0 / 3388.0);
    const auto fs_1215_77_5 = std::sqrt(7381125.0 / 5929.0);
    const auto fs_130_21_21 = std::sqrt(16900.0 / 21.0);
    const auto fs_130_77_6 = std::sqrt(101400.0 / 5929.0);
    const auto fs_130_7_15 = std::sqrt(253500.0 / 49.0);
    const auto fs_130_7_21 = std::sqrt(50700.0 / 7.0);
    const auto fs_133_1287_2310 = std::sqrt(1238230.0 / 50193.0);
    const auto fs_133_2574_330 = std::sqrt(88445.0 / 100386.0);
    const auto fs_133_2574_6006 = std::sqrt(123823.0 / 7722.0);
    const auto fs_133_429_14 = std::sqrt(247646.0 / 184041.0);
    const auto fs_133_429_286 = std::sqrt(35378.0 / 1287.0);
    const auto fs_133_429_55 = std::sqrt(88445.0 / 16731.0);
    const auto fs_133_429_715 = std::sqrt(88445.0 / 1287.0);
    const auto fs_133_858_10 = std::sqrt(88445.0 / 368082.0);
    const auto fs_133_858_1430 = std::sqrt(88445.0 / 2574.0);
    const auto fs_133_858_2 = std::sqrt(17689.0 / 368082.0);
    const auto fs_133_858_286 = std::sqrt(17689.0 / 2574.0);
    const auto fs_135_28_35 = std::sqrt(91125.0 / 112.0);
    const auto fs_135_4_5 = std::sqrt(91125.0 / 16.0);
    const auto fs_136_33_5 = std::sqrt(92480.0 / 1089.0);
    const auto fs_136_99_15 = std::sqrt(92480.0 / 3267.0);
    const auto fs_13_99_6 = std::sqrt(338.0 / 3267.0);
    const auto fs_14_1287_2310 = std::sqrt(13720.0 / 50193.0);
    const auto fs_14_429_14 = std::sqrt(2744.0 / 184041.0);
    const auto fs_14_429_286 = std::sqrt(392.0 / 1287.0);
    const auto fs_14_429_55 = std::sqrt(980.0 / 16731.0);
    const auto fs_14_429_715 = std::sqrt(980.0 / 1287.0);
    const auto fs_170_11_5 = std::sqrt(144500.0 / 121.0);
    const auto fs_170_33_15 = std::sqrt(144500.0 / 363.0);
    const auto fs_17_33_42 = std::sqrt(4046.0 / 363.0);
    const auto fs_18_1001_35 = std::sqrt(1620.0 / 143143.0);
    const auto fs_18_143_5 = std::sqrt(1620.0 / 20449.0);
    const auto fs_195_14_6 = std::sqrt(114075.0 / 98.0);
    const auto fs_195_4_15 = std::sqrt(570375.0 / 16.0);
    const auto fs_195_7_21 = std::sqrt(114075.0 / 7.0);
    const auto fs_195_8_21 = std::sqrt(798525.0 / 64.0);
    const auto fs_1_33_42 = std::sqrt(14.0 / 363.0);
    const auto fs_20_231_15 = std::sqrt(2000.0 / 17787.0);
    const auto fs_20_693_21 = std::sqrt(400.0 / 22869.0);
    const auto fs_20_99_3 = std::sqrt(400.0 / 3267.0);
    const auto fs_221_99_6 = std::sqrt(97682.0 / 3267.0);
    const auto fs_2430_1001_5 = std::sqrt(29524500.0 / 1002001.0);
    const auto fs_25_693_42 = std::sqrt(1250.0 / 22869.0);
    const auto fs_2600_693_3 = std::sqrt(6760000.0 / 160083.0);
    const auto fs_260_231_21 = std::sqrt(67600.0 / 2541.0);
    const auto fs_260_693_30 = std::sqrt(676000.0 / 160083.0);
    const auto fs_260_99_6 = std::sqrt(135200.0 / 3267.0);
    const auto fs_266_1287_858 = std::sqrt(141512.0 / 3861.0);
    const auto fs_266_429_10 = std::sqrt(707560.0 / 184041.0);
    const auto fs_266_429_70 = std::sqrt(4952920.0 / 184041.0);
    const auto fs_266_429_77 = std::sqrt(495292.0 / 16731.0);
    const auto fs_28_1287_858 = std::sqrt(1568.0 / 3861.0);
    const auto fs_28_429_10 = std::sqrt(7840.0 / 184041.0);
    const auto fs_28_429_70 = std::sqrt(54880.0 / 184041.0);
    const auto fs_28_429_77 = std::sqrt(5488.0 / 16731.0);
    const auto fs_2_33_11 = std::sqrt(4.0 / 99.0);
    const auto fs_2_33_30 = std::sqrt(40.0 / 363.0);
    const auto fs_2_99_105 = std::sqrt(140.0 / 3267.0);
    const auto fs_2_99_210 = std::sqrt(280.0 / 3267.0);
    const auto fs_2_99_42 = std::sqrt(56.0 / 3267.0);
    const auto fs_2_99_462 = std::sqrt(56.0 / 297.0);
    const auto fs_325_14_42 = std::sqrt(316875.0 / 14.0);
    const auto fs_325_16_42 = std::sqrt(2218125.0 / 128.0);
    const auto fs_325_42_42 = std::sqrt(105625.0 / 42.0);
    const auto fs_325_4_3 = std::sqrt(316875.0 / 16.0);
    const auto fs_340_99_3 = std::sqrt(115600.0 / 3267.0);
    const auto fs_34_33_11 = std::sqrt(1156.0 / 99.0);
    const auto fs_34_33_30 = std::sqrt(11560.0 / 363.0);
    const auto fs_34_99_105 = std::sqrt(40460.0 / 3267.0);
    const auto fs_34_99_210 = std::sqrt(80920.0 / 3267.0);
    const auto fs_34_99_42 = std::sqrt(16184.0 / 3267.0);
    const auto fs_34_99_462 = std::sqrt(16184.0 / 297.0);
    const auto fs_390_7_15 = std::sqrt(2281500.0 / 49.0);
    const auto fs_405_1001_35 = std::sqrt(820125.0 / 143143.0);
    const auto fs_405_143_5 = std::sqrt(820125.0 / 20449.0);
    const auto fs_405_14_5 = std::sqrt(820125.0 / 196.0);
    const auto fs_405_154_35 = std::sqrt(820125.0 / 3388.0);
    const auto fs_405_22_5 = std::sqrt(820125.0 / 484.0);
    const auto fs_405_28_35 = std::sqrt(820125.0 / 112.0);
    const auto fs_425_33_3 = std::sqrt(180625.0 / 363.0);
    const auto fs_455_8_6 = std::sqrt(621075.0 / 32.0);
    const auto fs_4_33_11 = std::sqrt(16.0 / 99.0);
    const auto fs_4_33_7 = std::sqrt(112.0 / 1089.0);
    const auto fs_4_99_30 = std::sqrt(160.0 / 3267.0);
    const auto fs_4_99_66 = std::sqrt(32.0 / 297.0);
    const auto fs_520_231_15 = std::sqrt(1352000.0 / 17787.0);
    const auto fs_520_693_21 = std::sqrt(270400.0 / 22869.0);
    const auto fs_532_429_11 = std::sqrt(283024.0 / 16731.0);
    const auto fs_54_1001_35 = std::sqrt(14580.0 / 143143.0);
    const auto fs_56_429_11 = std::sqrt(3136.0 / 16731.0);
    const auto fs_585_14_6 = std::sqrt(1026675.0 / 98.0);
    const auto fs_585_16_6 = std::sqrt(1026675.0 / 128.0);
    const auto fs_5_77_6 = std::sqrt(150.0 / 5929.0);
    const auto fs_650_21_3 = std::sqrt(422500.0 / 147.0);
    const auto fs_650_693_42 = std::sqrt(845000.0 / 22869.0);
    const auto fs_650_7_3 = std::sqrt(1267500.0 / 49.0);
    const auto fs_65_1_6 = std::sqrt(25350.0);
    const auto fs_65_21_30 = std::sqrt(42250.0 / 147.0);
    const auto fs_65_3_6 = std::sqrt(8450.0 / 3.0);
    const auto fs_65_4_21 = std::sqrt(88725.0 / 16.0);
    const auto fs_65_7_21 = std::sqrt(12675.0 / 7.0);
    const auto fs_65_7_30 = std::sqrt(126750.0 / 49.0);
    const auto fs_65_8_30 = std::sqrt(63375.0 / 32.0);
    const auto fs_665_1287_66 = std::sqrt(884450.0 / 50193.0);
    const auto fs_665_429_14 = std::sqrt(6191150.0 / 184041.0);
    const auto fs_68_33_11 = std::sqrt(4624.0 / 99.0);
    const auto fs_68_33_7 = std::sqrt(32368.0 / 1089.0);
    const auto fs_68_99_30 = std::sqrt(46240.0 / 3267.0);
    const auto fs_68_99_66 = std::sqrt(9248.0 / 297.0);
    const auto fs_70_1287_66 = std::sqrt(9800.0 / 50193.0);
    const auto fs_70_429_14 = std::sqrt(68600.0 / 184041.0);
    const auto fs_7_1287_330 = std::sqrt(490.0 / 50193.0);
    const auto fs_7_1287_6006 = std::sqrt(686.0 / 3861.0);
    const auto fs_7_429_10 = std::sqrt(490.0 / 184041.0);
    const auto fs_7_429_1430 = std::sqrt(490.0 / 1287.0);
    const auto fs_7_429_2 = std::sqrt(98.0 / 184041.0);
    const auto fs_7_429_286 = std::sqrt(98.0 / 1287.0);
    const auto fs_85_11_11 = std::sqrt(7225.0 / 11.0);
    const auto fs_85_11_7 = std::sqrt(50575.0 / 121.0);
    const auto fs_85_22_11 = std::sqrt(7225.0 / 44.0);
    const auto fs_85_22_30 = std::sqrt(108375.0 / 242.0);
    const auto fs_85_33_30 = std::sqrt(72250.0 / 363.0);
    const auto fs_85_33_66 = std::sqrt(14450.0 / 33.0);
    const auto fs_85_44_42 = std::sqrt(151725.0 / 968.0);
    const auto fs_85_66_105 = std::sqrt(252875.0 / 1452.0);
    const auto fs_85_66_210 = std::sqrt(252875.0 / 726.0);
    const auto fs_85_66_42 = std::sqrt(50575.0 / 726.0);
    const auto fs_85_66_462 = std::sqrt(50575.0 / 66.0);
    const auto fs_8_33_5 = std::sqrt(320.0 / 1089.0);
    const auto fs_8_99_15 = std::sqrt(320.0 / 3267.0);
    const auto fs_931_429_10 = std::sqrt(8667610.0 / 184041.0);
    const auto fs_931_429_2 = std::sqrt(1733522.0 / 184041.0);
    const auto fs_98_429_10 = std::sqrt(96040.0 / 184041.0);
    const auto fs_98_429_2 = std::sqrt(19208.0 / 184041.0);

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph8_0, ph8_p1, ph8_p7, ph8_p8, ab_2, pc_0, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * f_1155_16 + e_1 * f_455_2 * h2_0 - e_1 * f_1925_8 * r_2 + e_2 * f_135_2 * h4_0 - e_2 * f_260_1 * r_2 * h2_0 + e_2 * f_385_2 * r_4 + e_3 * f_170_33 * h6_0 - e_3 * f_405_11 * r_2 * h4_0 + e_3 * f_260_3 * r_4 * h2_0 - e_3 * f_55_1 * r_6 + e_4 * f_133_1287 * h8_0 - e_4 * fs_133_429_715 * h8_p8 - e_4 * f_136_99 * r_2 * h6_0 + e_4 * f_810_143 * r_4 * h4_0 - e_4 * f_1040_99 * r_6 * h2_0 + e_4 * f_55_9 * r_8 - e_5 * f_14_1287 * r_2 * h8_0 + e_5 * fs_14_429_715 * r_2 * h8_p8 + e_5 * f_8_99 * r_4 * h6_0 - e_5 * f_36_143 * r_6 * h4_0 + e_5 * f_40_99 * r_8 * h2_0 - e_5 * f_2_9 * r_10;

        pc_1[k] = - e_1 * fs_455_8_6 * h2_p1 - e_2 * fs_135_4_5 * h4_p1 + e_2 * fs_65_1_6 * r_2 * h2_p1 - e_3 * fs_85_66_42 * h6_p1 + e_3 * fs_405_22_5 * r_2 * h4_p1 - e_3 * fs_65_3_6 * r_4 * h2_p1 - e_4 * fs_133_858_2 * h8_p1 - e_4 * fs_133_858_1430 * h8_p7 + e_4 * fs_34_99_42 * r_2 * h6_p1 - e_4 * fs_405_143_5 * r_4 * h4_p1 + e_4 * fs_260_99_6 * r_6 * h2_p1 + e_5 * fs_7_429_2 * r_2 * h8_p1 + e_5 * fs_7_429_1430 * r_2 * h8_p7 - e_5 * fs_2_99_42 * r_4 * h6_p1 + e_5 * fs_18_143_5 * r_6 * h4_p1 - e_5 * fs_10_99_6 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph6_p5, ph6_p6, ph8_p2, ph8_p3, ph8_p5, ph8_p6, ab_2, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_2[k] = e_1 * fs_65_4_21 * h2_p2 + e_2 * fs_405_28_35 * h4_p2 - e_2 * fs_130_7_21 * r_2 * h2_p2 + e_3 * fs_85_33_30 * h6_p2 - e_3 * fs_85_33_66 * h6_p6 - e_3 * fs_1215_154_35 * r_2 * h4_p2 + e_3 * fs_130_21_21 * r_4 * h2_p2 + e_4 * fs_133_858_10 * h8_p2 - e_4 * fs_133_2574_6006 * h8_p6 - e_4 * fs_68_99_30 * r_2 * h6_p2 + e_4 * fs_68_99_66 * r_2 * h6_p6 + e_4 * fs_1215_1001_35 * r_4 * h4_p2 - e_4 * fs_520_693_21 * r_6 * h2_p2 - e_5 * fs_7_429_10 * r_2 * h8_p2 + e_5 * fs_7_1287_6006 * r_2 * h8_p6 + e_5 * fs_4_99_30 * r_4 * h6_p2 - e_5 * fs_4_99_66 * r_4 * h6_p6 - e_5 * fs_54_1001_35 * r_6 * h4_p2 + e_5 * fs_20_693_21 * r_8 * h2_p2;

        pc_3[k] = - e_2 * fs_135_4_5 * h4_p3 - e_3 * fs_170_33_15 * h6_p3 - e_3 * fs_85_11_11 * h6_p5 + e_3 * fs_405_22_5 * r_2 * h4_p3 - e_4 * fs_133_2574_330 * h8_p3 - e_4 * fs_133_858_286 * h8_p5 + e_4 * fs_136_99_15 * r_2 * h6_p3 + e_4 * fs_68_33_11 * r_2 * h6_p5 - e_4 * fs_405_143_5 * r_4 * h4_p3 + e_5 * fs_7_1287_330 * r_2 * h8_p3 + e_5 * fs_7_429_286 * r_2 * h8_p5 - e_5 * fs_8_99_15 * r_4 * h6_p3 - e_5 * fs_4_33_11 * r_4 * h6_p5 + e_5 * fs_18_143_5 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_m4, ph4_m3, ph6_m5, ph6_m4, ph6_m3, ph8_m5, ph8_m4, ph8_m3, ab_2, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];

        pc_4[k] = e_2 * f_135_2 * h4_m4 + e_3 * fs_170_11_5 * h6_m4 - e_3 * f_405_11 * r_2 * h4_m4 + e_4 * fs_133_429_55 * h8_m4 - e_4 * fs_136_33_5 * r_2 * h6_m4 + e_4 * f_810_143 * r_4 * h4_m4 - e_5 * fs_14_429_55 * r_2 * h8_m4 + e_5 * fs_8_33_5 * r_4 * h6_m4 - e_5 * f_36_143 * r_6 * h4_m4;

        pc_5[k] = - e_2 * fs_135_4_5 * h4_m3 + e_3 * fs_85_11_11 * h6_m5 - e_3 * fs_170_33_15 * h6_m3 + e_3 * fs_405_22_5 * r_2 * h4_m3 + e_4 * fs_133_858_286 * h8_m5 - e_4 * fs_133_2574_330 * h8_m3 - e_4 * fs_68_33_11 * r_2 * h6_m5 + e_4 * fs_136_99_15 * r_2 * h6_m3 - e_4 * fs_405_143_5 * r_4 * h4_m3 - e_5 * fs_7_429_286 * r_2 * h8_m5 + e_5 * fs_7_1287_330 * r_2 * h8_m3 + e_5 * fs_4_33_11 * r_4 * h6_m5 - e_5 * fs_8_99_15 * r_4 * h6_m3 + e_5 * fs_18_143_5 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m8, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2, pc_6, pc_7, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_6[k] = e_1 * fs_65_4_21 * h2_m2 + e_2 * fs_405_28_35 * h4_m2 - e_2 * fs_130_7_21 * r_2 * h2_m2 + e_3 * fs_85_33_66 * h6_m6 + e_3 * fs_85_33_30 * h6_m2 - e_3 * fs_1215_154_35 * r_2 * h4_m2 + e_3 * fs_130_21_21 * r_4 * h2_m2 + e_4 * fs_133_2574_6006 * h8_m6 + e_4 * fs_133_858_10 * h8_m2 - e_4 * fs_68_99_66 * r_2 * h6_m6 - e_4 * fs_68_99_30 * r_2 * h6_m2 + e_4 * fs_1215_1001_35 * r_4 * h4_m2 - e_4 * fs_520_693_21 * r_6 * h2_m2 - e_5 * fs_7_1287_6006 * r_2 * h8_m6 - e_5 * fs_7_429_10 * r_2 * h8_m2 + e_5 * fs_4_99_66 * r_4 * h6_m6 + e_5 * fs_4_99_30 * r_4 * h6_m2 - e_5 * fs_54_1001_35 * r_6 * h4_m2 + e_5 * fs_20_693_21 * r_8 * h2_m2;

        pc_7[k] = - e_1 * fs_455_8_6 * h2_m1 - e_2 * fs_135_4_5 * h4_m1 + e_2 * fs_65_1_6 * r_2 * h2_m1 - e_3 * fs_85_66_42 * h6_m1 + e_3 * fs_405_22_5 * r_2 * h4_m1 - e_3 * fs_65_3_6 * r_4 * h2_m1 + e_4 * fs_133_858_1430 * h8_m7 - e_4 * fs_133_858_2 * h8_m1 + e_4 * fs_34_99_42 * r_2 * h6_m1 - e_4 * fs_405_143_5 * r_4 * h4_m1 + e_4 * fs_260_99_6 * r_6 * h2_m1 - e_5 * fs_7_429_1430 * r_2 * h8_m7 + e_5 * fs_7_429_2 * r_2 * h8_m1 - e_5 * fs_2_99_42 * r_4 * h6_m1 + e_5 * fs_18_143_5 * r_6 * h4_m1 - e_5 * fs_10_99_6 * r_8 * h2_m1;

        pc_8[k] = e_4 * fs_133_429_715 * h8_m8 - e_5 * fs_14_429_715 * r_2 * h8_m8;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ab_2, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];

        pc_9[k] = e_0 * f_1155_16 + e_1 * f_455_8 * h2_0 - e_1 * f_1925_8 * r_2 - e_2 * f_405_4 * h4_0 - e_2 * f_65_1 * r_2 * h2_0 + e_2 * f_385_2 * r_4 - e_3 * f_1445_66 * h6_0 + e_3 * fs_85_66_462 * h6_p6 + e_3 * f_1215_22 * r_2 * h4_0 + e_3 * f_65_3 * r_4 * h2_0 - e_3 * f_55_1 * r_6 - e_4 * f_1064_1287 * h8_0 - e_4 * fs_266_1287_858 * h8_p6 + e_4 * f_578_99 * r_2 * h6_0 - e_4 * fs_34_99_462 * r_2 * h6_p6 - e_4 * f_1215_143 * r_4 * h4_0 - e_4 * f_260_99 * r_6 * h2_0 + e_4 * f_55_9 * r_8 + e_5 * f_112_1287 * r_2 * h8_0 + e_5 * fs_28_1287_858 * r_2 * h8_p6 - e_5 * f_34_99 * r_4 * h6_0 + e_5 * fs_2_99_462 * r_4 * h6_p6 + e_5 * f_54_143 * r_6 * h4_0 + e_5 * f_10_99 * r_8 * h2_0 - e_5 * f_2_9 * r_10;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ab_2, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];

        pc_10[k] = - e_1 * fs_325_16_42 * h2_p1 + e_2 * fs_135_28_35 * h4_p1 + e_2 * fs_325_14_42 * r_2 * h2_p1 + e_3 * fs_1105_132_6 * h6_p1 + e_3 * fs_85_22_11 * h6_p5 - e_3 * fs_405_154_35 * r_2 * h4_p1 - e_3 * fs_325_42_42 * r_4 * h2_p1 + e_4 * fs_133_429_14 * h8_p1 - e_4 * fs_133_429_286 * h8_p5 - e_4 * fs_221_99_6 * r_2 * h6_p1 - e_4 * fs_34_33_11 * r_2 * h6_p5 + e_4 * fs_405_1001_35 * r_4 * h4_p1 + e_4 * fs_650_693_42 * r_6 * h2_p1 - e_5 * fs_14_429_14 * r_2 * h8_p1 + e_5 * fs_14_429_286 * r_2 * h8_p5 + e_5 * fs_13_99_6 * r_4 * h6_p1 + e_5 * fs_2_33_11 * r_4 * h6_p5 - e_5 * fs_18_1001_35 * r_6 * h4_p1 - e_5 * fs_25_693_42 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_m3, ph4_p2, ph4_p4, ph6_m3, ph6_p2, ph6_p4, ph8_m3, ph8_p2, ph8_p4, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];

        pc_11[k] = e_1 * fs_195_8_21 * h2_p2 + e_2 * fs_135_28_35 * h4_p2 + e_2 * fs_135_4_5 * h4_p4 - e_2 * fs_195_7_21 * r_2 * h2_p2 - e_3 * fs_85_22_30 * h6_p2 - e_3 * f_85_22 * h6_p4 - e_3 * fs_405_154_35 * r_2 * h4_p2 - e_3 * fs_405_22_5 * r_2 * h4_p4 + e_3 * fs_65_7_21 * r_4 * h2_p2 - e_4 * fs_266_429_10 * h8_p2 - e_4 * fs_532_429_11 * h8_p4 + e_4 * fs_34_33_30 * r_2 * h6_p2 + e_4 * f_34_33 * r_2 * h6_p4 + e_4 * fs_405_1001_35 * r_4 * h4_p2 + e_4 * fs_405_143_5 * r_4 * h4_p4 - e_4 * fs_260_231_21 * r_6 * h2_p2 + e_5 * fs_28_429_10 * r_2 * h8_p2 + e_5 * fs_56_429_11 * r_2 * h8_p4 - e_5 * fs_2_33_30 * r_4 * h6_p2 - e_5 * f_2_33 * r_4 * h6_p4 - e_5 * fs_18_1001_35 * r_6 * h4_p2 - e_5 * fs_18_143_5 * r_6 * h4_p4 + e_5 * fs_10_231_21 * r_8 * h2_p2;

        pc_12[k] = - e_2 * f_405_4 * h4_m3 + e_3 * fs_425_33_3 * h6_m3 + e_3 * f_1215_22 * r_2 * h4_m3 + e_4 * fs_665_1287_66 * h8_m3 - e_4 * fs_340_99_3 * r_2 * h6_m3 - e_4 * f_1215_143 * r_4 * h4_m3 - e_5 * fs_70_1287_66 * r_2 * h8_m3 + e_5 * fs_20_99_3 * r_4 * h6_m3 + e_5 * f_54_143 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ab_2, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];

        pc_13[k] = e_1 * fs_195_8_21 * h2_m2 - e_2 * fs_135_4_5 * h4_m4 + e_2 * fs_135_28_35 * h4_m2 - e_2 * fs_195_7_21 * r_2 * h2_m2 + e_3 * f_85_22 * h6_m4 - e_3 * fs_85_22_30 * h6_m2 + e_3 * fs_405_22_5 * r_2 * h4_m4 - e_3 * fs_405_154_35 * r_2 * h4_m2 + e_3 * fs_65_7_21 * r_4 * h2_m2 + e_4 * fs_532_429_11 * h8_m4 - e_4 * fs_266_429_10 * h8_m2 - e_4 * f_34_33 * r_2 * h6_m4 + e_4 * fs_34_33_30 * r_2 * h6_m2 - e_4 * fs_405_143_5 * r_4 * h4_m4 + e_4 * fs_405_1001_35 * r_4 * h4_m2 - e_4 * fs_260_231_21 * r_6 * h2_m2 - e_5 * fs_56_429_11 * r_2 * h8_m4 + e_5 * fs_28_429_10 * r_2 * h8_m2 + e_5 * f_2_33 * r_4 * h6_m4 - e_5 * fs_2_33_30 * r_4 * h6_m2 + e_5 * fs_18_143_5 * r_6 * h4_m4 - e_5 * fs_18_1001_35 * r_6 * h4_m2 + e_5 * fs_10_231_21 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m6, ph6_m5, ph6_m1, ph8_m7, ph8_m6, ph8_m5, ph8_m1, ab_2, pc_14, pc_15, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];

        pc_14[k] = - e_1 * fs_325_16_42 * h2_m1 + e_2 * fs_135_28_35 * h4_m1 + e_2 * fs_325_14_42 * r_2 * h2_m1 - e_3 * fs_85_22_11 * h6_m5 + e_3 * fs_1105_132_6 * h6_m1 - e_3 * fs_405_154_35 * r_2 * h4_m1 - e_3 * fs_325_42_42 * r_4 * h2_m1 + e_4 * fs_133_429_286 * h8_m5 + e_4 * fs_133_429_14 * h8_m1 + e_4 * fs_34_33_11 * r_2 * h6_m5 - e_4 * fs_221_99_6 * r_2 * h6_m1 + e_4 * fs_405_1001_35 * r_4 * h4_m1 + e_4 * fs_650_693_42 * r_6 * h2_m1 - e_5 * fs_14_429_286 * r_2 * h8_m5 - e_5 * fs_14_429_14 * r_2 * h8_m1 - e_5 * fs_2_33_11 * r_4 * h6_m5 + e_5 * fs_13_99_6 * r_4 * h6_m1 - e_5 * fs_18_1001_35 * r_6 * h4_m1 - e_5 * fs_25_693_42 * r_8 * h2_m1;

        pc_15[k] = - e_3 * fs_85_66_462 * h6_m6 + e_4 * fs_266_1287_858 * h8_m6 + e_4 * fs_34_99_462 * r_2 * h6_m6 - e_5 * fs_28_1287_858 * r_2 * h8_m6 - e_5 * fs_2_99_462 * r_4 * h6_m6;

        pc_16[k] = e_1 * fs_455_8_6 * h2_m1 + e_2 * fs_135_4_5 * h4_m1 - e_2 * fs_65_1_6 * r_2 * h2_m1 + e_3 * fs_85_66_42 * h6_m1 - e_3 * fs_405_22_5 * r_2 * h4_m1 + e_3 * fs_65_3_6 * r_4 * h2_m1 + e_4 * fs_133_858_1430 * h8_m7 + e_4 * fs_133_858_2 * h8_m1 - e_4 * fs_34_99_42 * r_2 * h6_m1 + e_4 * fs_405_143_5 * r_4 * h4_m1 - e_4 * fs_260_99_6 * r_6 * h2_m1 - e_5 * fs_7_429_1430 * r_2 * h8_m7 - e_5 * fs_7_429_2 * r_2 * h8_m1 + e_5 * fs_2_99_42 * r_4 * h6_m1 - e_5 * fs_18_143_5 * r_6 * h4_m1 + e_5 * fs_10_99_6 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ab_2, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_17[k] = e_0 * f_1155_16 - e_1 * f_65_1 * h2_0 - e_1 * f_1925_8 * r_2 - e_2 * f_1485_28 * h4_0 - e_2 * fs_405_28_35 * h4_p4 + e_2 * f_520_7 * r_2 * h2_0 + e_2 * f_385_2 * r_4 + e_3 * f_85_3 * h6_0 + e_3 * fs_85_11_7 * h6_p4 + e_3 * f_405_14 * r_2 * h4_0 + e_3 * fs_1215_154_35 * r_2 * h4_p4 - e_3 * f_520_21 * r_4 * h2_0 - e_3 * f_55_1 * r_6 + e_4 * f_3724_1287 * h8_0 - e_4 * fs_266_429_77 * h8_p4 - e_4 * f_68_9 * r_2 * h6_0 - e_4 * fs_68_33_7 * r_2 * h6_p4 - e_4 * f_405_91 * r_4 * h4_0 - e_4 * fs_1215_1001_35 * r_4 * h4_p4 + e_4 * f_2080_693 * r_6 * h2_0 + e_4 * f_55_9 * r_8 - e_5 * f_392_1287 * r_2 * h8_0 + e_5 * fs_28_429_77 * r_2 * h8_p4 + e_5 * f_4_9 * r_4 * h6_0 + e_5 * fs_4_33_7 * r_4 * h6_p4 + e_5 * f_18_91 * r_6 * h4_0 + e_5 * fs_54_1001_35 * r_6 * h4_p4 - e_5 * f_80_693 * r_8 * h2_0 - e_5 * f_2_9 * r_10;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_m2, ph8_p1, ph8_p3, ab_2, pc_18, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_18[k] = - e_1 * fs_585_16_6 * h2_p1 + e_2 * fs_405_14_5 * h4_p1 - e_2 * fs_135_28_35 * h4_p3 + e_2 * fs_585_14_6 * r_2 * h2_p1 - e_3 * fs_85_44_42 * h6_p1 + e_3 * fs_85_66_105 * h6_p3 - e_3 * fs_1215_77_5 * r_2 * h4_p1 + e_3 * fs_405_154_35 * r_2 * h4_p3 - e_3 * fs_195_14_6 * r_4 * h2_p1 - e_4 * fs_931_429_2 * h8_p1 - e_4 * fs_133_1287_2310 * h8_p3 + e_4 * fs_17_33_42 * r_2 * h6_p1 - e_4 * fs_34_99_105 * r_2 * h6_p3 + e_4 * fs_2430_1001_5 * r_4 * h4_p1 - e_4 * fs_405_1001_35 * r_4 * h4_p3 + e_4 * fs_130_77_6 * r_6 * h2_p1 + e_5 * fs_98_429_2 * r_2 * h8_p1 + e_5 * fs_14_1287_2310 * r_2 * h8_p3 - e_5 * fs_1_33_42 * r_4 * h6_p1 + e_5 * fs_2_99_105 * r_4 * h6_p3 - e_5 * fs_108_1001_5 * r_6 * h4_p1 + e_5 * fs_18_1001_35 * r_6 * h4_p3 - e_5 * fs_5_77_6 * r_8 * h2_p1;

        pc_19[k] = e_1 * fs_195_4_15 * h2_m2 - e_2 * f_1485_28 * h4_m2 - e_2 * fs_390_7_15 * r_2 * h2_m2 + e_3 * f_405_14 * r_2 * h4_m2 + e_3 * fs_130_7_15 * r_4 * h2_m2 + e_4 * fs_665_429_14 * h8_m2 - e_4 * f_405_91 * r_4 * h4_m2 - e_4 * fs_520_231_15 * r_6 * h2_m2 - e_5 * fs_70_429_14 * r_2 * h8_m2 + e_5 * f_18_91 * r_6 * h4_m2 + e_5 * fs_20_231_15 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_20, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];

        pc_20[k] = - e_1 * fs_585_16_6 * h2_m1 + e_2 * fs_135_28_35 * h4_m3 + e_2 * fs_405_14_5 * h4_m1 + e_2 * fs_585_14_6 * r_2 * h2_m1 - e_3 * fs_85_66_105 * h6_m3 - e_3 * fs_85_44_42 * h6_m1 - e_3 * fs_405_154_35 * r_2 * h4_m3 - e_3 * fs_1215_77_5 * r_2 * h4_m1 - e_3 * fs_195_14_6 * r_4 * h2_m1 + e_4 * fs_133_1287_2310 * h8_m3 - e_4 * fs_931_429_2 * h8_m1 + e_4 * fs_34_99_105 * r_2 * h6_m3 + e_4 * fs_17_33_42 * r_2 * h6_m1 + e_4 * fs_405_1001_35 * r_4 * h4_m3 + e_4 * fs_2430_1001_5 * r_4 * h4_m1 + e_4 * fs_130_77_6 * r_6 * h2_m1 - e_5 * fs_14_1287_2310 * r_2 * h8_m3 + e_5 * fs_98_429_2 * r_2 * h8_m1 - e_5 * fs_2_99_105 * r_4 * h6_m3 - e_5 * fs_1_33_42 * r_4 * h6_m1 - e_5 * fs_18_1001_35 * r_6 * h4_m3 - e_5 * fs_108_1001_5 * r_6 * h4_m1 - e_5 * fs_5_77_6 * r_8 * h2_m1;

        pc_21[k] = e_2 * fs_405_28_35 * h4_m4 - e_3 * fs_85_11_7 * h6_m4 - e_3 * fs_1215_154_35 * r_2 * h4_m4 + e_4 * fs_266_429_77 * h8_m4 + e_4 * fs_68_33_7 * r_2 * h6_m4 + e_4 * fs_1215_1001_35 * r_4 * h4_m4 - e_5 * fs_28_429_77 * r_2 * h8_m4 - e_5 * fs_4_33_7 * r_4 * h6_m4 - e_5 * fs_54_1001_35 * r_6 * h4_m4;

        pc_22[k] = e_1 * fs_325_16_42 * h2_m1 - e_2 * fs_135_28_35 * h4_m1 - e_2 * fs_325_14_42 * r_2 * h2_m1 - e_3 * fs_85_22_11 * h6_m5 - e_3 * fs_1105_132_6 * h6_m1 + e_3 * fs_405_154_35 * r_2 * h4_m1 + e_3 * fs_325_42_42 * r_4 * h2_m1 + e_4 * fs_133_429_286 * h8_m5 - e_4 * fs_133_429_14 * h8_m1 + e_4 * fs_34_33_11 * r_2 * h6_m5 + e_4 * fs_221_99_6 * r_2 * h6_m1 - e_4 * fs_405_1001_35 * r_4 * h4_m1 - e_4 * fs_650_693_42 * r_6 * h2_m1 - e_5 * fs_14_429_286 * r_2 * h8_m5 + e_5 * fs_14_429_14 * r_2 * h8_m1 - e_5 * fs_2_33_11 * r_4 * h6_m5 - e_5 * fs_13_99_6 * r_4 * h6_m1 + e_5 * fs_18_1001_35 * r_6 * h4_m1 + e_5 * fs_25_693_42 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ab_2, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];

        pc_23[k] = - e_1 * fs_65_4_21 * h2_m2 - e_2 * fs_405_28_35 * h4_m2 + e_2 * fs_130_7_21 * r_2 * h2_m2 + e_3 * fs_85_33_66 * h6_m6 - e_3 * fs_85_33_30 * h6_m2 + e_3 * fs_1215_154_35 * r_2 * h4_m2 - e_3 * fs_130_21_21 * r_4 * h2_m2 + e_4 * fs_133_2574_6006 * h8_m6 - e_4 * fs_133_858_10 * h8_m2 - e_4 * fs_68_99_66 * r_2 * h6_m6 + e_4 * fs_68_99_30 * r_2 * h6_m2 - e_4 * fs_1215_1001_35 * r_4 * h4_m2 + e_4 * fs_520_693_21 * r_6 * h2_m2 - e_5 * fs_7_1287_6006 * r_2 * h8_m6 + e_5 * fs_7_429_10 * r_2 * h8_m2 + e_5 * fs_4_99_66 * r_4 * h6_m6 - e_5 * fs_4_99_30 * r_4 * h6_m2 + e_5 * fs_54_1001_35 * r_6 * h4_m2 - e_5 * fs_20_693_21 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ab_2, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];

        pc_24[k] = e_0 * f_1155_16 - e_1 * f_1105_8 * h2_0 + e_1 * fs_325_4_3 * h2_p2 - e_1 * f_1925_8 * r_2 + e_2 * f_1215_28 * h4_0 - e_2 * fs_405_14_5 * h4_p2 + e_2 * f_1105_7 * r_2 * h2_0 - e_2 * fs_650_7_3 * r_2 * h2_p2 + e_2 * f_385_2 * r_4 + e_3 * f_85_66 * h6_0 + e_3 * fs_85_66_210 * h6_p2 - e_3 * f_3645_154 * r_2 * h4_0 + e_3 * fs_1215_77_5 * r_2 * h4_p2 - e_3 * f_1105_21 * r_4 * h2_0 + e_3 * fs_650_21_3 * r_4 * h2_p2 - e_3 * f_55_1 * r_6 - e_4 * f_7448_1287 * h8_0 - e_4 * fs_266_429_70 * h8_p2 - e_4 * f_34_99 * r_2 * h6_0 - e_4 * fs_34_99_210 * r_2 * h6_p2 + e_4 * f_3645_1001 * r_4 * h4_0 - e_4 * fs_2430_1001_5 * r_4 * h4_p2 + e_4 * f_4420_693 * r_6 * h2_0 - e_4 * fs_2600_693_3 * r_6 * h2_p2 + e_4 * f_55_9 * r_8 + e_5 * f_784_1287 * r_2 * h8_0 + e_5 * fs_28_429_70 * r_2 * h8_p2 + e_5 * f_2_99 * r_4 * h6_0 + e_5 * fs_2_99_210 * r_4 * h6_p2 - e_5 * f_162_1001 * r_6 * h4_0 + e_5 * fs_108_1001_5 * r_6 * h4_p2 - e_5 * f_170_693 * r_8 * h2_0 + e_5 * fs_100_693_3 * r_8 * h2_p2 - e_5 * f_2_9 * r_10;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m3, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_25, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_25[k] = - e_1 * fs_65_8_30 * h2_m1 + e_2 * f_1215_28 * h4_m1 + e_2 * fs_65_7_30 * r_2 * h2_m1 - e_3 * fs_85_66_210 * h6_m1 - e_3 * f_3645_154 * r_2 * h4_m1 - e_3 * fs_65_21_30 * r_4 * h2_m1 + e_4 * fs_931_429_10 * h8_m1 + e_4 * fs_34_99_210 * r_2 * h6_m1 + e_4 * f_3645_1001 * r_4 * h4_m1 + e_4 * fs_260_693_30 * r_6 * h2_m1 - e_5 * fs_98_429_10 * r_2 * h8_m1 - e_5 * fs_2_99_210 * r_4 * h6_m1 - e_5 * f_162_1001 * r_6 * h4_m1 - e_5 * fs_10_693_30 * r_8 * h2_m1;

        pc_26[k] = - e_1 * fs_325_4_3 * h2_m2 + e_2 * fs_405_14_5 * h4_m2 + e_2 * fs_650_7_3 * r_2 * h2_m2 - e_3 * fs_85_66_210 * h6_m2 - e_3 * fs_1215_77_5 * r_2 * h4_m2 - e_3 * fs_650_21_3 * r_4 * h2_m2 + e_4 * fs_266_429_70 * h8_m2 + e_4 * fs_34_99_210 * r_2 * h6_m2 + e_4 * fs_2430_1001_5 * r_4 * h4_m2 + e_4 * fs_2600_693_3 * r_6 * h2_m2 - e_5 * fs_28_429_70 * r_2 * h8_m2 - e_5 * fs_2_99_210 * r_4 * h6_m2 - e_5 * fs_108_1001_5 * r_6 * h4_m2 - e_5 * fs_100_693_3 * r_8 * h2_m2;

        pc_27[k] = e_1 * fs_585_16_6 * h2_m1 + e_2 * fs_135_28_35 * h4_m3 - e_2 * fs_405_14_5 * h4_m1 - e_2 * fs_585_14_6 * r_2 * h2_m1 - e_3 * fs_85_66_105 * h6_m3 + e_3 * fs_85_44_42 * h6_m1 - e_3 * fs_405_154_35 * r_2 * h4_m3 + e_3 * fs_1215_77_5 * r_2 * h4_m1 + e_3 * fs_195_14_6 * r_4 * h2_m1 + e_4 * fs_133_1287_2310 * h8_m3 + e_4 * fs_931_429_2 * h8_m1 + e_4 * fs_34_99_105 * r_2 * h6_m3 - e_4 * fs_17_33_42 * r_2 * h6_m1 + e_4 * fs_405_1001_35 * r_4 * h4_m3 - e_4 * fs_2430_1001_5 * r_4 * h4_m1 - e_4 * fs_130_77_6 * r_6 * h2_m1 - e_5 * fs_14_1287_2310 * r_2 * h8_m3 - e_5 * fs_98_429_2 * r_2 * h8_m1 - e_5 * fs_2_99_105 * r_4 * h6_m3 + e_5 * fs_1_33_42 * r_4 * h6_m1 - e_5 * fs_18_1001_35 * r_6 * h4_m3 + e_5 * fs_108_1001_5 * r_6 * h4_m1 + e_5 * fs_5_77_6 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_28, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_28[k] = - e_1 * fs_195_8_21 * h2_m2 - e_2 * fs_135_4_5 * h4_m4 - e_2 * fs_135_28_35 * h4_m2 + e_2 * fs_195_7_21 * r_2 * h2_m2 + e_3 * f_85_22 * h6_m4 + e_3 * fs_85_22_30 * h6_m2 + e_3 * fs_405_22_5 * r_2 * h4_m4 + e_3 * fs_405_154_35 * r_2 * h4_m2 - e_3 * fs_65_7_21 * r_4 * h2_m2 + e_4 * fs_532_429_11 * h8_m4 + e_4 * fs_266_429_10 * h8_m2 - e_4 * f_34_33 * r_2 * h6_m4 - e_4 * fs_34_33_30 * r_2 * h6_m2 - e_4 * fs_405_143_5 * r_4 * h4_m4 - e_4 * fs_405_1001_35 * r_4 * h4_m2 + e_4 * fs_260_231_21 * r_6 * h2_m2 - e_5 * fs_56_429_11 * r_2 * h8_m4 - e_5 * fs_28_429_10 * r_2 * h8_m2 + e_5 * f_2_33 * r_4 * h6_m4 + e_5 * fs_2_33_30 * r_4 * h6_m2 + e_5 * fs_18_143_5 * r_6 * h4_m4 + e_5 * fs_18_1001_35 * r_6 * h4_m2 - e_5 * fs_10_231_21 * r_8 * h2_m2;

        pc_29[k] = e_2 * fs_135_4_5 * h4_m3 + e_3 * fs_85_11_11 * h6_m5 + e_3 * fs_170_33_15 * h6_m3 - e_3 * fs_405_22_5 * r_2 * h4_m3 + e_4 * fs_133_858_286 * h8_m5 + e_4 * fs_133_2574_330 * h8_m3 - e_4 * fs_68_33_11 * r_2 * h6_m5 - e_4 * fs_136_99_15 * r_2 * h6_m3 + e_4 * fs_405_143_5 * r_4 * h4_m3 - e_5 * fs_7_429_286 * r_2 * h8_m5 - e_5 * fs_7_1287_330 * r_2 * h8_m3 + e_5 * fs_4_33_11 * r_4 * h6_m5 + e_5 * fs_8_99_15 * r_4 * h6_m3 - e_5 * fs_18_143_5 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph6_0, ph6_p1, ph8_0, ph8_p1, ph8_p2, ab_2, pc_30, pc_31, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_30[k] = e_0 * f_1155_16 - e_1 * f_325_2 * h2_0 - e_1 * f_1925_8 * r_2 + e_2 * f_1215_14 * h4_0 + e_2 * f_1300_7 * r_2 * h2_0 + e_2 * f_385_2 * r_4 - e_3 * f_850_33 * h6_0 - e_3 * f_3645_77 * r_2 * h4_0 - e_3 * f_1300_21 * r_4 * h2_0 - e_3 * f_55_1 * r_6 + e_4 * f_9310_1287 * h8_0 + e_4 * f_680_99 * r_2 * h6_0 + e_4 * f_7290_1001 * r_4 * h4_0 + e_4 * f_5200_693 * r_6 * h2_0 + e_4 * f_55_9 * r_8 - e_5 * f_980_1287 * r_2 * h8_0 - e_5 * f_40_99 * r_4 * h6_0 - e_5 * f_324_1001 * r_6 * h4_0 - e_5 * f_200_693 * r_8 * h2_0 - e_5 * f_2_9 * r_10;

        pc_31[k] = - e_1 * fs_65_8_30 * h2_p1 + e_2 * f_1215_28 * h4_p1 + e_2 * fs_65_7_30 * r_2 * h2_p1 - e_3 * fs_85_66_210 * h6_p1 - e_3 * f_3645_154 * r_2 * h4_p1 - e_3 * fs_65_21_30 * r_4 * h2_p1 + e_4 * fs_931_429_10 * h8_p1 + e_4 * fs_34_99_210 * r_2 * h6_p1 + e_4 * f_3645_1001 * r_4 * h4_p1 + e_4 * fs_260_693_30 * r_6 * h2_p1 - e_5 * fs_98_429_10 * r_2 * h8_p1 - e_5 * fs_2_99_210 * r_4 * h6_p1 - e_5 * f_162_1001 * r_6 * h4_p1 - e_5 * fs_10_693_30 * r_8 * h2_p1;

        pc_32[k] = e_1 * fs_195_4_15 * h2_p2 - e_2 * f_1485_28 * h4_p2 - e_2 * fs_390_7_15 * r_2 * h2_p2 + e_3 * f_405_14 * r_2 * h4_p2 + e_3 * fs_130_7_15 * r_4 * h2_p2 + e_4 * fs_665_429_14 * h8_p2 - e_4 * f_405_91 * r_4 * h4_p2 - e_4 * fs_520_231_15 * r_6 * h2_p2 - e_5 * fs_70_429_14 * r_2 * h8_p2 + e_5 * f_18_91 * r_6 * h4_p2 + e_5 * fs_20_231_15 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ab_2, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_33[k] = - e_2 * f_405_4 * h4_p3 + e_3 * fs_425_33_3 * h6_p3 + e_3 * f_1215_22 * r_2 * h4_p3 + e_4 * fs_665_1287_66 * h8_p3 - e_4 * fs_340_99_3 * r_2 * h6_p3 - e_4 * f_1215_143 * r_4 * h4_p3 - e_5 * fs_70_1287_66 * r_2 * h8_p3 + e_5 * fs_20_99_3 * r_4 * h6_p3 + e_5 * f_54_143 * r_6 * h4_p3;

        pc_34[k] = e_2 * f_135_2 * h4_p4 + e_3 * fs_170_11_5 * h6_p4 - e_3 * f_405_11 * r_2 * h4_p4 + e_4 * fs_133_429_55 * h8_p4 - e_4 * fs_136_33_5 * r_2 * h6_p4 + e_4 * f_810_143 * r_4 * h4_p4 - e_5 * fs_14_429_55 * r_2 * h8_p4 + e_5 * fs_8_33_5 * r_4 * h6_p4 - e_5 * f_36_143 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ab_2, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];

        pc_35[k] = e_0 * f_1155_16 - e_1 * f_1105_8 * h2_0 - e_1 * fs_325_4_3 * h2_p2 - e_1 * f_1925_8 * r_2 + e_2 * f_1215_28 * h4_0 + e_2 * fs_405_14_5 * h4_p2 + e_2 * f_1105_7 * r_2 * h2_0 + e_2 * fs_650_7_3 * r_2 * h2_p2 + e_2 * f_385_2 * r_4 + e_3 * f_85_66 * h6_0 - e_3 * fs_85_66_210 * h6_p2 - e_3 * f_3645_154 * r_2 * h4_0 - e_3 * fs_1215_77_5 * r_2 * h4_p2 - e_3 * f_1105_21 * r_4 * h2_0 - e_3 * fs_650_21_3 * r_4 * h2_p2 - e_3 * f_55_1 * r_6 - e_4 * f_7448_1287 * h8_0 + e_4 * fs_266_429_70 * h8_p2 - e_4 * f_34_99 * r_2 * h6_0 + e_4 * fs_34_99_210 * r_2 * h6_p2 + e_4 * f_3645_1001 * r_4 * h4_0 + e_4 * fs_2430_1001_5 * r_4 * h4_p2 + e_4 * f_4420_693 * r_6 * h2_0 + e_4 * fs_2600_693_3 * r_6 * h2_p2 + e_4 * f_55_9 * r_8 + e_5 * f_784_1287 * r_2 * h8_0 - e_5 * fs_28_429_70 * r_2 * h8_p2 + e_5 * f_2_99 * r_4 * h6_0 - e_5 * fs_2_99_210 * r_4 * h6_p2 - e_5 * f_162_1001 * r_6 * h4_0 - e_5 * fs_108_1001_5 * r_6 * h4_p2 - e_5 * f_170_693 * r_8 * h2_0 - e_5 * fs_100_693_3 * r_8 * h2_p2 - e_5 * f_2_9 * r_10;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_36[k] = - e_1 * fs_585_16_6 * h2_p1 + e_2 * fs_405_14_5 * h4_p1 + e_2 * fs_135_28_35 * h4_p3 + e_2 * fs_585_14_6 * r_2 * h2_p1 - e_3 * fs_85_44_42 * h6_p1 - e_3 * fs_85_66_105 * h6_p3 - e_3 * fs_1215_77_5 * r_2 * h4_p1 - e_3 * fs_405_154_35 * r_2 * h4_p3 - e_3 * fs_195_14_6 * r_4 * h2_p1 - e_4 * fs_931_429_2 * h8_p1 + e_4 * fs_133_1287_2310 * h8_p3 + e_4 * fs_17_33_42 * r_2 * h6_p1 + e_4 * fs_34_99_105 * r_2 * h6_p3 + e_4 * fs_2430_1001_5 * r_4 * h4_p1 + e_4 * fs_405_1001_35 * r_4 * h4_p3 + e_4 * fs_130_77_6 * r_6 * h2_p1 + e_5 * fs_98_429_2 * r_2 * h8_p1 - e_5 * fs_14_1287_2310 * r_2 * h8_p3 - e_5 * fs_1_33_42 * r_4 * h6_p1 - e_5 * fs_2_99_105 * r_4 * h6_p3 - e_5 * fs_108_1001_5 * r_6 * h4_p1 - e_5 * fs_18_1001_35 * r_6 * h4_p3 - e_5 * fs_5_77_6 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ab_2, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_37[k] = e_1 * fs_195_8_21 * h2_p2 + e_2 * fs_135_28_35 * h4_p2 - e_2 * fs_135_4_5 * h4_p4 - e_2 * fs_195_7_21 * r_2 * h2_p2 - e_3 * fs_85_22_30 * h6_p2 + e_3 * f_85_22 * h6_p4 - e_3 * fs_405_154_35 * r_2 * h4_p2 + e_3 * fs_405_22_5 * r_2 * h4_p4 + e_3 * fs_65_7_21 * r_4 * h2_p2 - e_4 * fs_266_429_10 * h8_p2 + e_4 * fs_532_429_11 * h8_p4 + e_4 * fs_34_33_30 * r_2 * h6_p2 - e_4 * f_34_33 * r_2 * h6_p4 + e_4 * fs_405_1001_35 * r_4 * h4_p2 - e_4 * fs_405_143_5 * r_4 * h4_p4 - e_4 * fs_260_231_21 * r_6 * h2_p2 + e_5 * fs_28_429_10 * r_2 * h8_p2 - e_5 * fs_56_429_11 * r_2 * h8_p4 - e_5 * fs_2_33_30 * r_4 * h6_p2 + e_5 * f_2_33 * r_4 * h6_p4 - e_5 * fs_18_1001_35 * r_6 * h4_p2 + e_5 * fs_18_143_5 * r_6 * h4_p4 + e_5 * fs_10_231_21 * r_8 * h2_p2;

        pc_38[k] = - e_2 * fs_135_4_5 * h4_p3 - e_3 * fs_170_33_15 * h6_p3 + e_3 * fs_85_11_11 * h6_p5 + e_3 * fs_405_22_5 * r_2 * h4_p3 - e_4 * fs_133_2574_330 * h8_p3 + e_4 * fs_133_858_286 * h8_p5 + e_4 * fs_136_99_15 * r_2 * h6_p3 - e_4 * fs_68_33_11 * r_2 * h6_p5 - e_4 * fs_405_143_5 * r_4 * h4_p3 + e_5 * fs_7_1287_330 * r_2 * h8_p3 - e_5 * fs_7_429_286 * r_2 * h8_p5 - e_5 * fs_8_99_15 * r_4 * h6_p3 + e_5 * fs_4_33_11 * r_4 * h6_p5 + e_5 * fs_18_143_5 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ab_2, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_39[k] = e_0 * f_1155_16 - e_1 * f_65_1 * h2_0 - e_1 * f_1925_8 * r_2 - e_2 * f_1485_28 * h4_0 + e_2 * fs_405_28_35 * h4_p4 + e_2 * f_520_7 * r_2 * h2_0 + e_2 * f_385_2 * r_4 + e_3 * f_85_3 * h6_0 - e_3 * fs_85_11_7 * h6_p4 + e_3 * f_405_14 * r_2 * h4_0 - e_3 * fs_1215_154_35 * r_2 * h4_p4 - e_3 * f_520_21 * r_4 * h2_0 - e_3 * f_55_1 * r_6 + e_4 * f_3724_1287 * h8_0 + e_4 * fs_266_429_77 * h8_p4 - e_4 * f_68_9 * r_2 * h6_0 + e_4 * fs_68_33_7 * r_2 * h6_p4 - e_4 * f_405_91 * r_4 * h4_0 + e_4 * fs_1215_1001_35 * r_4 * h4_p4 + e_4 * f_2080_693 * r_6 * h2_0 + e_4 * f_55_9 * r_8 - e_5 * f_392_1287 * r_2 * h8_0 - e_5 * fs_28_429_77 * r_2 * h8_p4 + e_5 * f_4_9 * r_4 * h6_0 - e_5 * fs_4_33_7 * r_4 * h6_p4 + e_5 * f_18_91 * r_6 * h4_0 - e_5 * fs_54_1001_35 * r_6 * h4_p4 - e_5 * f_80_693 * r_8 * h2_0 - e_5 * f_2_9 * r_10;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph6_p5, ph6_p6, ph8_p1, ph8_p2, ph8_p5, ph8_p6, ab_2, pc_40, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_40[k] = - e_1 * fs_325_16_42 * h2_p1 + e_2 * fs_135_28_35 * h4_p1 + e_2 * fs_325_14_42 * r_2 * h2_p1 + e_3 * fs_1105_132_6 * h6_p1 - e_3 * fs_85_22_11 * h6_p5 - e_3 * fs_405_154_35 * r_2 * h4_p1 - e_3 * fs_325_42_42 * r_4 * h2_p1 + e_4 * fs_133_429_14 * h8_p1 + e_4 * fs_133_429_286 * h8_p5 - e_4 * fs_221_99_6 * r_2 * h6_p1 + e_4 * fs_34_33_11 * r_2 * h6_p5 + e_4 * fs_405_1001_35 * r_4 * h4_p1 + e_4 * fs_650_693_42 * r_6 * h2_p1 - e_5 * fs_14_429_14 * r_2 * h8_p1 - e_5 * fs_14_429_286 * r_2 * h8_p5 + e_5 * fs_13_99_6 * r_4 * h6_p1 - e_5 * fs_2_33_11 * r_4 * h6_p5 - e_5 * fs_18_1001_35 * r_6 * h4_p1 - e_5 * fs_25_693_42 * r_8 * h2_p1;

        pc_41[k] = e_1 * fs_65_4_21 * h2_p2 + e_2 * fs_405_28_35 * h4_p2 - e_2 * fs_130_7_21 * r_2 * h2_p2 + e_3 * fs_85_33_30 * h6_p2 + e_3 * fs_85_33_66 * h6_p6 - e_3 * fs_1215_154_35 * r_2 * h4_p2 + e_3 * fs_130_21_21 * r_4 * h2_p2 + e_4 * fs_133_858_10 * h8_p2 + e_4 * fs_133_2574_6006 * h8_p6 - e_4 * fs_68_99_30 * r_2 * h6_p2 - e_4 * fs_68_99_66 * r_2 * h6_p6 + e_4 * fs_1215_1001_35 * r_4 * h4_p2 - e_4 * fs_520_693_21 * r_6 * h2_p2 - e_5 * fs_7_429_10 * r_2 * h8_p2 - e_5 * fs_7_1287_6006 * r_2 * h8_p6 + e_5 * fs_4_99_30 * r_4 * h6_p2 + e_5 * fs_4_99_66 * r_4 * h6_p6 - e_5 * fs_54_1001_35 * r_6 * h4_p2 + e_5 * fs_20_693_21 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p6, ph8_0, ph8_p1, ph8_p6, ph8_p7, ab_2, pc_42, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_42[k] = e_0 * f_1155_16 + e_1 * f_455_8 * h2_0 - e_1 * f_1925_8 * r_2 - e_2 * f_405_4 * h4_0 - e_2 * f_65_1 * r_2 * h2_0 + e_2 * f_385_2 * r_4 - e_3 * f_1445_66 * h6_0 - e_3 * fs_85_66_462 * h6_p6 + e_3 * f_1215_22 * r_2 * h4_0 + e_3 * f_65_3 * r_4 * h2_0 - e_3 * f_55_1 * r_6 - e_4 * f_1064_1287 * h8_0 + e_4 * fs_266_1287_858 * h8_p6 + e_4 * f_578_99 * r_2 * h6_0 + e_4 * fs_34_99_462 * r_2 * h6_p6 - e_4 * f_1215_143 * r_4 * h4_0 - e_4 * f_260_99 * r_6 * h2_0 + e_4 * f_55_9 * r_8 + e_5 * f_112_1287 * r_2 * h8_0 - e_5 * fs_28_1287_858 * r_2 * h8_p6 - e_5 * f_34_99 * r_4 * h6_0 - e_5 * fs_2_99_462 * r_4 * h6_p6 + e_5 * f_54_143 * r_6 * h4_0 + e_5 * f_10_99 * r_8 * h2_0 - e_5 * f_2_9 * r_10;

        pc_43[k] = - e_1 * fs_455_8_6 * h2_p1 - e_2 * fs_135_4_5 * h4_p1 + e_2 * fs_65_1_6 * r_2 * h2_p1 - e_3 * fs_85_66_42 * h6_p1 + e_3 * fs_405_22_5 * r_2 * h4_p1 - e_3 * fs_65_3_6 * r_4 * h2_p1 - e_4 * fs_133_858_2 * h8_p1 + e_4 * fs_133_858_1430 * h8_p7 + e_4 * fs_34_99_42 * r_2 * h6_p1 - e_4 * fs_405_143_5 * r_4 * h4_p1 + e_4 * fs_260_99_6 * r_6 * h2_p1 + e_5 * fs_7_429_2 * r_2 * h8_p1 - e_5 * fs_7_429_1430 * r_2 * h8_p7 - e_5 * fs_2_99_42 * r_4 * h6_p1 + e_5 * fs_18_143_5 * r_6 * h4_p1 - e_5 * fs_10_99_6 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];

        pc_44[k] = e_0 * f_1155_16 + e_1 * f_455_2 * h2_0 - e_1 * f_1925_8 * r_2 + e_2 * f_135_2 * h4_0 - e_2 * f_260_1 * r_2 * h2_0 + e_2 * f_385_2 * r_4 + e_3 * f_170_33 * h6_0 - e_3 * f_405_11 * r_2 * h4_0 + e_3 * f_260_3 * r_4 * h2_0 - e_3 * f_55_1 * r_6 + e_4 * f_133_1287 * h8_0 + e_4 * fs_133_429_715 * h8_p8 - e_4 * f_136_99 * r_2 * h6_0 + e_4 * f_810_143 * r_4 * h4_0 - e_4 * f_1040_99 * r_6 * h2_0 + e_4 * f_55_9 * r_8 - e_5 * f_14_1287 * r_2 * h8_0 - e_5 * fs_14_429_715 * r_2 * h8_p8 + e_5 * f_8_99 * r_4 * h6_0 - e_5 * f_36_143 * r_6 * h4_0 + e_5 * f_40_99 * r_8 * h2_0 - e_5 * f_2_9 * r_10;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[81] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 10, 11, 12, 13, 14, 15, 16, 17, 2, 11, 20, 21, 22, 23, 24, 25, 26, 3, 12, 21, 30, 31, 32, 33, 34, 35, 4, 13, 22, 31, 40, 41, 42, 43, 44, 5, 14, 23, 32, 41, 50, 51, 52, 53, 6, 15, 24, 33, 42, 51, 60, 61, 62, 7, 16, 25, 34, 43, 52, 61, 70, 71, 8, 17, 26, 35, 44, 53, 62, 71, 80};

    for (size_t n = 0; n < 81; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
