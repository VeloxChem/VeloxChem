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



#include "SimdKineticEnergyRecGH.hpp"

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
compute_gh_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gh_kinetic_energy: Basis functions must be of angular momenta four and five"));
    }

    if (harmonics.size() < 9)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gh_kinetic_energy: Harmonics must reach angular momentum 9"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_gh_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 99 * nvalues, 0.0);

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

            const auto ff_0 = fbase * aexp * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * aexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * aexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * aexp * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * aexp * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

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
    const auto *ph9_m9 = harmonics[8].data(0);
    const auto *ph9_m8 = harmonics[8].data(1);
    const auto *ph9_m7 = harmonics[8].data(2);
    const auto *ph9_m6 = harmonics[8].data(3);
    const auto *ph9_m5 = harmonics[8].data(4);
    const auto *ph9_m4 = harmonics[8].data(5);
    const auto *ph9_m3 = harmonics[8].data(6);
    const auto *ph9_m2 = harmonics[8].data(7);
    const auto *ph9_m1 = harmonics[8].data(8);
    const auto *ph9_0 = harmonics[8].data(9);
    const auto *ph9_p1 = harmonics[8].data(10);
    const auto *ph9_p2 = harmonics[8].data(11);
    const auto *ph9_p3 = harmonics[8].data(12);
    const auto *ph9_p4 = harmonics[8].data(13);
    const auto *ph9_p5 = harmonics[8].data(14);
    const auto *ph9_p6 = harmonics[8].data(15);
    const auto *ph9_p7 = harmonics[8].data(16);
    const auto *ph9_p8 = harmonics[8].data(17);
    const auto *ph9_p9 = harmonics[8].data(18);

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
    auto *pc_63 = values + 63 * nvalues;
    auto *pc_64 = values + 64 * nvalues;
    auto *pc_65 = values + 65 * nvalues;
    auto *pc_66 = values + 66 * nvalues;
    auto *pc_67 = values + 67 * nvalues;
    auto *pc_68 = values + 68 * nvalues;
    auto *pc_69 = values + 69 * nvalues;
    auto *pc_70 = values + 70 * nvalues;
    auto *pc_71 = values + 71 * nvalues;
    auto *pc_72 = values + 72 * nvalues;
    auto *pc_73 = values + 73 * nvalues;
    auto *pc_74 = values + 74 * nvalues;
    auto *pc_75 = values + 75 * nvalues;
    auto *pc_76 = values + 76 * nvalues;
    auto *pc_77 = values + 77 * nvalues;
    auto *pc_78 = values + 78 * nvalues;
    auto *pc_79 = values + 79 * nvalues;
    auto *pc_80 = values + 80 * nvalues;
    auto *pc_81 = values + 81 * nvalues;
    auto *pc_82 = values + 82 * nvalues;
    auto *pc_83 = values + 83 * nvalues;
    auto *pc_84 = values + 84 * nvalues;
    auto *pc_85 = values + 85 * nvalues;
    auto *pc_86 = values + 86 * nvalues;
    auto *pc_87 = values + 87 * nvalues;
    auto *pc_88 = values + 88 * nvalues;
    auto *pc_89 = values + 89 * nvalues;
    auto *pc_90 = values + 90 * nvalues;
    auto *pc_91 = values + 91 * nvalues;
    auto *pc_92 = values + 92 * nvalues;
    auto *pc_93 = values + 93 * nvalues;
    auto *pc_94 = values + 94 * nvalues;
    auto *pc_95 = values + 95 * nvalues;
    auto *pc_96 = values + 96 * nvalues;
    auto *pc_97 = values + 97 * nvalues;
    auto *pc_98 = values + 98 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_102_13 = 102.0 / 13.0;
    const auto f_10_33 = 10.0 / 33.0;
    const auto f_1200_143 = 1200.0 / 143.0;
    const auto f_1260_11 = 1260.0 / 11.0;
    const auto f_12635_429 = 12635.0 / 429.0;
    const auto f_13300_2431 = 13300.0 / 2431.0;
    const auto f_13300_429 = 13300.0 / 429.0;
    const auto f_1330_143 = 1330.0 / 143.0;
    const auto f_1330_221 = 1330.0 / 221.0;
    const auto f_1365_2 = 1365.0 / 2.0;
    const auto f_1365_4 = 1365.0 / 4.0;
    const auto f_1680_143 = 1680.0 / 143.0;
    const auto f_1764_2431 = 1764.0 / 2431.0;
    const auto f_17_13 = 17.0 / 13.0;
    const auto f_18522_2431 = 18522.0 / 2431.0;
    const auto f_210_1 = 210.0;
    const auto f_2205_8 = 2205.0 / 8.0;
    const auto f_245_1 = 245.0;
    const auto f_252_2431 = 252.0 / 2431.0;
    const auto f_255_2 = 255.0 / 2.0;
    const auto f_255_26 = 255.0 / 26.0;
    const auto f_260_3 = 260.0 / 3.0;
    const auto f_260_33 = 260.0 / 33.0;
    const auto f_2646_2431 = 2646.0 / 2431.0;
    const auto f_2660_7293 = 2660.0 / 7293.0;
    const auto f_2800_7293 = 2800.0 / 7293.0;
    const auto f_280_2431 = 280.0 / 2431.0;
    const auto f_28_143 = 28.0 / 143.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_2_39 = 2.0 / 39.0;
    const auto f_300_1 = 300.0;
    const auto f_325_3 = 325.0 / 3.0;
    const auto f_325_33 = 325.0 / 33.0;
    const auto f_3325_143 = 3325.0 / 143.0;
    const auto f_390_1 = 390.0;
    const auto f_4095_16 = 4095.0 / 16.0;
    const auto f_4095_8 = 4095.0 / 8.0;
    const auto f_40_143 = 40.0 / 143.0;
    const auto f_420_1 = 420.0;
    const auto f_42_2431 = 42.0 / 2431.0;
    const auto f_441_2431 = 441.0 / 2431.0;
    const auto f_4_13 = 4.0 / 13.0;
    const auto f_50540_7293 = 50540.0 / 7293.0;
    const auto f_510_13 = 510.0 / 13.0;
    const auto f_53200_7293 = 53200.0 / 7293.0;
    const auto f_5320_2431 = 5320.0 / 2431.0;
    const auto f_56_143 = 56.0 / 143.0;
    const auto f_585_2 = 585.0 / 2.0;
    const auto f_630_11 = 630.0 / 11.0;
    const auto f_65_1 = 65.0;
    const auto f_65_11 = 65.0 / 11.0;
    const auto f_665_26 = 665.0 / 26.0;
    const auto f_675_2 = 675.0 / 2.0;
    const auto f_6825_16 = 6825.0 / 16.0;
    const auto f_6825_8 = 6825.0 / 8.0;
    const auto f_68_13 = 68.0 / 13.0;
    const auto f_700_2431 = 700.0 / 2431.0;
    const auto f_70_221 = 70.0 / 221.0;
    const auto f_735_11 = 735.0 / 11.0;
    const auto f_765_13 = 765.0 / 13.0;
    const auto f_840_143 = 840.0 / 143.0;
    const auto f_85_1 = 85.0;
    const auto f_85_4 = 85.0 / 4.0;
    const auto f_8_33 = 8.0 / 33.0;
    const auto f_8_39 = 8.0 / 39.0;
    const auto f_900_11 = 900.0 / 11.0;
    const auto f_945_2 = 945.0 / 2.0;
    const auto f_945_4 = 945.0 / 4.0;
    const auto f_975_2 = 975.0 / 2.0;
    const auto f_980_143 = 980.0 / 143.0;
    const auto f_98_429 = 98.0 / 429.0;
    const auto fs_100_143_105 = std::sqrt(1050000.0 / 20449.0);
    const auto fs_10290_2431_3 = std::sqrt(317652300.0 / 5909761.0);
    const auto fs_1029_2431_165 = std::sqrt(15882615.0 / 537251.0);
    const auto fs_105_11_3 = std::sqrt(33075.0 / 121.0);
    const auto fs_10_13_70 = std::sqrt(7000.0 / 169.0);
    const auto fs_10_143_7 = std::sqrt(700.0 / 20449.0);
    const auto fs_10_187_7 = std::sqrt(700.0 / 34969.0);
    const auto fs_10_1_21 = std::sqrt(2100.0);
    const auto fs_10_2431_154 = std::sqrt(1400.0 / 537251.0);
    const auto fs_10_2431_70 = std::sqrt(7000.0 / 5909761.0);
    const auto fs_10_2431_858 = std::sqrt(600.0 / 41327.0);
    const auto fs_10_429_105 = std::sqrt(3500.0 / 61347.0);
    const auto fs_10_561_33 = std::sqrt(100.0 / 9537.0);
    const auto fs_1125_16_10 = std::sqrt(6328125.0 / 128.0);
    const auto fs_115_2431_14 = std::sqrt(185150.0 / 5909761.0);
    const auto fs_1176_2431_2 = std::sqrt(2765952.0 / 5909761.0);
    const auto fs_12160_7293_21 = std::sqrt(1035059200.0 / 17729283.0);
    const auto fs_12348_2431_2 = std::sqrt(304946208.0 / 5909761.0);
    const auto fs_125_2_10 = std::sqrt(78125.0 / 2.0);
    const auto fs_125_561_2 = std::sqrt(31250.0 / 314721.0);
    const auto fs_130_33_6 = std::sqrt(33800.0 / 363.0);
    const auto fs_130_3_6 = std::sqrt(33800.0 / 3.0);
    const auto fs_1330_143_11 = std::sqrt(1768900.0 / 1859.0);
    const auto fs_1330_2431_22 = std::sqrt(3537800.0 / 537251.0);
    const auto fs_1330_429_14 = std::sqrt(24764600.0 / 184041.0);
    const auto fs_1330_429_165 = std::sqrt(8844500.0 / 5577.0);
    const auto fs_1330_7293_858 = std::sqrt(3537800.0 / 123981.0);
    const auto fs_135_11_21 = std::sqrt(382725.0 / 121.0);
    const auto fs_135_11_35 = std::sqrt(637875.0 / 121.0);
    const auto fs_1365_16_10 = std::sqrt(9316125.0 / 128.0);
    const auto fs_1365_16_14 = std::sqrt(13042575.0 / 128.0);
    const auto fs_1365_16_15 = std::sqrt(27948375.0 / 256.0);
    const auto fs_1365_16_2 = std::sqrt(1863225.0 / 128.0);
    const auto fs_1365_16_21 = std::sqrt(39127725.0 / 256.0);
    const auto fs_1365_16_3 = std::sqrt(5589675.0 / 256.0);
    const auto fs_1365_16_42 = std::sqrt(39127725.0 / 128.0);
    const auto fs_1365_16_6 = std::sqrt(5589675.0 / 128.0);
    const auto fs_1365_32_2 = std::sqrt(1863225.0 / 512.0);
    const auto fs_1365_32_42 = std::sqrt(39127725.0 / 512.0);
    const auto fs_1365_32_6 = std::sqrt(5589675.0 / 512.0);
    const auto fs_1365_4_6 = std::sqrt(5589675.0 / 8.0);
    const auto fs_1365_8_10 = std::sqrt(9316125.0 / 32.0);
    const auto fs_1365_8_14 = std::sqrt(13042575.0 / 32.0);
    const auto fs_1365_8_15 = std::sqrt(27948375.0 / 64.0);
    const auto fs_1365_8_21 = std::sqrt(39127725.0 / 64.0);
    const auto fs_1365_8_3 = std::sqrt(5589675.0 / 64.0);
    const auto fs_1365_8_6 = std::sqrt(5589675.0 / 32.0);
    const auto fs_140_143_3 = std::sqrt(58800.0 / 20449.0);
    const auto fs_140_1_3 = std::sqrt(58800.0);
    const auto fs_140_2431_14 = std::sqrt(274400.0 / 5909761.0);
    const auto fs_140_7293_6 = std::sqrt(39200.0 / 17729283.0);
    const auto fs_1425_2431_30 = std::sqrt(60918750.0 / 5909761.0);
    const auto fs_1425_286_66 = std::sqrt(6091875.0 / 3718.0);
    const auto fs_1425_572_30 = std::sqrt(30459375.0 / 163592.0);
    const auto fs_147_2431_1001 = std::sqrt(151263.0 / 41327.0);
    const auto fs_147_2431_1155 = std::sqrt(2268945.0 / 537251.0);
    const auto fs_147_2431_12155 = std::sqrt(108045.0 / 2431.0);
    const auto fs_147_2431_143 = std::sqrt(21609.0 / 41327.0);
    const auto fs_147_2431_33 = std::sqrt(64827.0 / 537251.0);
    const auto fs_147_2431_4290 = std::sqrt(648270.0 / 41327.0);
    const auto fs_147_4862_2 = std::sqrt(21609.0 / 11819522.0);
    const auto fs_147_4862_22 = std::sqrt(21609.0 / 1074502.0);
    const auto fs_147_4862_30030 = std::sqrt(2268945.0 / 82654.0);
    const auto fs_147_4862_330 = std::sqrt(324135.0 / 1074502.0);
    const auto fs_147_4862_6006 = std::sqrt(453789.0 / 82654.0);
    const auto fs_147_4862_770 = std::sqrt(756315.0 / 1074502.0);
    const auto fs_14_2431_1001 = std::sqrt(1372.0 / 41327.0);
    const auto fs_14_2431_1155 = std::sqrt(20580.0 / 537251.0);
    const auto fs_14_2431_12155 = std::sqrt(980.0 / 2431.0);
    const auto fs_14_2431_143 = std::sqrt(196.0 / 41327.0);
    const auto fs_14_2431_33 = std::sqrt(588.0 / 537251.0);
    const auto fs_14_2431_4290 = std::sqrt(5880.0 / 41327.0);
    const auto fs_14_429_3 = std::sqrt(196.0 / 61347.0);
    const auto fs_150_11_7 = std::sqrt(157500.0 / 121.0);
    const auto fs_150_143_30 = std::sqrt(675000.0 / 20449.0);
    const auto fs_150_1_7 = std::sqrt(157500.0);
    const auto fs_150_2431_66 = std::sqrt(135000.0 / 537251.0);
    const auto fs_1520_2431_10 = std::sqrt(23104000.0 / 5909761.0);
    const auto fs_1520_2431_165 = std::sqrt(34656000.0 / 537251.0);
    const auto fs_15485_1716_6 = std::sqrt(239785225.0 / 490776.0);
    const auto fs_15485_7293_6 = std::sqrt(479570450.0 / 17729283.0);
    const auto fs_1575_8_3 = std::sqrt(7441875.0 / 64.0);
    const auto fs_15_11_105 = std::sqrt(23625.0 / 121.0);
    const auto fs_15_11_15 = std::sqrt(3375.0 / 121.0);
    const auto fs_15_11_210 = std::sqrt(47250.0 / 121.0);
    const auto fs_15_1_10 = std::sqrt(2250.0);
    const auto fs_15_2431_210 = std::sqrt(47250.0 / 5909761.0);
    const auto fs_15_2431_4290 = std::sqrt(6750.0 / 41327.0);
    const auto fs_15_2_70 = std::sqrt(7875.0 / 2.0);
    const auto fs_168_2431_21 = std::sqrt(592704.0 / 5909761.0);
    const auto fs_168_2431_33 = std::sqrt(84672.0 / 537251.0);
    const auto fs_168_2431_55 = std::sqrt(141120.0 / 537251.0);
    const auto fs_175_1_3 = std::sqrt(91875.0);
    const auto fs_1764_2431_21 = std::sqrt(65345616.0 / 5909761.0);
    const auto fs_1764_2431_33 = std::sqrt(9335088.0 / 537251.0);
    const auto fs_1764_2431_55 = std::sqrt(15558480.0 / 537251.0);
    const auto fs_17_13_15 = std::sqrt(4335.0 / 169.0);
    const auto fs_17_13_21 = std::sqrt(6069.0 / 169.0);
    const auto fs_17_13_30 = std::sqrt(8670.0 / 169.0);
    const auto fs_17_13_35 = std::sqrt(10115.0 / 169.0);
    const auto fs_17_13_6 = std::sqrt(1734.0 / 169.0);
    const auto fs_17_26_10 = std::sqrt(1445.0 / 338.0);
    const auto fs_17_26_70 = std::sqrt(10115.0 / 338.0);
    const auto fs_180_143_21 = std::sqrt(680400.0 / 20449.0);
    const auto fs_180_143_35 = std::sqrt(1134000.0 / 20449.0);
    const auto fs_190_187_7 = std::sqrt(252700.0 / 34969.0);
    const auto fs_190_2431_154 = std::sqrt(505400.0 / 537251.0);
    const auto fs_190_2431_70 = std::sqrt(2527000.0 / 5909761.0);
    const auto fs_190_2431_858 = std::sqrt(216600.0 / 41327.0);
    const auto fs_190_561_33 = std::sqrt(36100.0 / 9537.0);
    const auto fs_1935_16_2 = std::sqrt(3744225.0 / 128.0);
    const auto fs_195_1_6 = std::sqrt(228150.0);
    const auto fs_195_2_10 = std::sqrt(190125.0 / 2.0);
    const auto fs_195_2_14 = std::sqrt(266175.0 / 2.0);
    const auto fs_195_2_15 = std::sqrt(570375.0 / 4.0);
    const auto fs_195_2_21 = std::sqrt(798525.0 / 4.0);
    const auto fs_195_2_3 = std::sqrt(114075.0 / 4.0);
    const auto fs_195_4_2 = std::sqrt(38025.0 / 8.0);
    const auto fs_195_4_42 = std::sqrt(798525.0 / 8.0);
    const auto fs_195_4_6 = std::sqrt(114075.0 / 8.0);
    const auto fs_196_2431_15 = std::sqrt(576240.0 / 5909761.0);
    const auto fs_1_11_10 = std::sqrt(10.0 / 121.0);
    const auto fs_1_33_2 = std::sqrt(2.0 / 1089.0);
    const auto fs_1_33_42 = std::sqrt(14.0 / 363.0);
    const auto fs_1_33_6 = std::sqrt(2.0 / 363.0);
    const auto fs_1_39_10 = std::sqrt(10.0 / 1521.0);
    const auto fs_1_39_70 = std::sqrt(70.0 / 1521.0);
    const auto fs_200_143_7 = std::sqrt(280000.0 / 20449.0);
    const auto fs_200_1_3 = std::sqrt(120000.0);
    const auto fs_200_2431_11 = std::sqrt(40000.0 / 537251.0);
    const auto fs_2058_2431_15 = std::sqrt(63530460.0 / 5909761.0);
    const auto fs_20_13_10 = std::sqrt(4000.0 / 169.0);
    const auto fs_20_143_105 = std::sqrt(42000.0 / 20449.0);
    const auto fs_20_143_15 = std::sqrt(6000.0 / 20449.0);
    const auto fs_20_143_210 = std::sqrt(84000.0 / 20449.0);
    const auto fs_20_143_7 = std::sqrt(2800.0 / 20449.0);
    const auto fs_20_1_105 = std::sqrt(42000.0);
    const auto fs_20_2431_2002 = std::sqrt(5600.0 / 41327.0);
    const auto fs_20_429_7 = std::sqrt(2800.0 / 184041.0);
    const auto fs_20_663_105 = std::sqrt(14000.0 / 146523.0);
    const auto fs_20_7293_210 = std::sqrt(28000.0 / 17729283.0);
    const auto fs_210_11_15 = std::sqrt(661500.0 / 121.0);
    const auto fs_210_11_5 = std::sqrt(220500.0 / 121.0);
    const auto fs_210_2431_33 = std::sqrt(132300.0 / 537251.0);
    const auto fs_215_2_2 = std::sqrt(46225.0 / 2.0);
    const auto fs_2185_2431_14 = std::sqrt(66839150.0 / 5909761.0);
    const auto fs_2185_572_14 = std::sqrt(33419575.0 / 163592.0);
    const auto fs_2185_858_30 = std::sqrt(23871125.0 / 122694.0);
    const auto fs_2185_858_70 = std::sqrt(167097875.0 / 368082.0);
    const auto fs_21_2431_10 = std::sqrt(4410.0 / 5909761.0);
    const auto fs_21_2431_110 = std::sqrt(4410.0 / 537251.0);
    const auto fs_21_2431_1430 = std::sqrt(4410.0 / 41327.0);
    const auto fs_21_2431_4290 = std::sqrt(13230.0 / 41327.0);
    const auto fs_21_2431_770 = std::sqrt(30870.0 / 537251.0);
    const auto fs_2205_2431_33 = std::sqrt(14586075.0 / 537251.0);
    const auto fs_225_11_7 = std::sqrt(354375.0 / 121.0);
    const auto fs_225_1_3 = std::sqrt(151875.0);
    const auto fs_225_22_30 = std::sqrt(759375.0 / 242.0);
    const auto fs_225_4_7 = std::sqrt(354375.0 / 16.0);
    const auto fs_225_8_105 = std::sqrt(5315625.0 / 64.0);
    const auto fs_230_7293_30 = std::sqrt(529000.0 / 17729283.0);
    const auto fs_230_7293_70 = std::sqrt(3703000.0 / 53187849.0);
    const auto fs_2375_132_2 = std::sqrt(5640625.0 / 8712.0);
    const auto fs_2375_286_7 = std::sqrt(39484375.0 / 81796.0);
    const auto fs_2375_561_2 = std::sqrt(11281250.0 / 314721.0);
    const auto fs_240_11_7 = std::sqrt(403200.0 / 121.0);
    const auto fs_250_143_10 = std::sqrt(625000.0 / 20449.0);
    const auto fs_250_2431_7 = std::sqrt(437500.0 / 5909761.0);
    const auto fs_255_13_5 = std::sqrt(325125.0 / 169.0);
    const auto fs_255_26_15 = std::sqrt(975375.0 / 676.0);
    const auto fs_255_26_21 = std::sqrt(1365525.0 / 676.0);
    const auto fs_255_26_30 = std::sqrt(975375.0 / 338.0);
    const auto fs_255_26_35 = std::sqrt(2275875.0 / 676.0);
    const auto fs_255_26_6 = std::sqrt(195075.0 / 338.0);
    const auto fs_255_52_10 = std::sqrt(325125.0 / 1352.0);
    const auto fs_255_52_70 = std::sqrt(2275875.0 / 1352.0);
    const auto fs_2565_286_6 = std::sqrt(19737675.0 / 40898.0);
    const auto fs_25_187_6 = std::sqrt(3750.0 / 34969.0);
    const auto fs_25_1_105 = std::sqrt(65625.0);
    const auto fs_25_429_10 = std::sqrt(6250.0 / 184041.0);
    const auto fs_25_7293_6006 = std::sqrt(8750.0 / 123981.0);
    const auto fs_265_7293_66 = std::sqrt(140450.0 / 1611753.0);
    const auto fs_2660_2431_14 = std::sqrt(99058400.0 / 5909761.0);
    const auto fs_2660_7293_6 = std::sqrt(14151200.0 / 17729283.0);
    const auto fs_270_2431_6 = std::sqrt(437400.0 / 5909761.0);
    const auto fs_2755_858_21 = std::sqrt(53130175.0 / 245388.0);
    const auto fs_280_143_15 = std::sqrt(1176000.0 / 20449.0);
    const auto fs_280_143_5 = std::sqrt(392000.0 / 20449.0);
    const auto fs_280_2431_11 = std::sqrt(78400.0 / 537251.0);
    const auto fs_280_7293_14 = std::sqrt(1097600.0 / 53187849.0);
    const auto fs_280_7293_165 = std::sqrt(392000.0 / 1611753.0);
    const auto fs_2850_2431_66 = std::sqrt(48735000.0 / 537251.0);
    const auto fs_285_11_6 = std::sqrt(487350.0 / 121.0);
    const auto fs_285_2431_210 = std::sqrt(17057250.0 / 5909761.0);
    const auto fs_285_2431_4290 = std::sqrt(2436750.0 / 41327.0);
    const auto fs_285_572_210 = std::sqrt(8528625.0 / 163592.0);
    const auto fs_285_572_4290 = std::sqrt(1218375.0 / 1144.0);
    const auto fs_28_2431_10 = std::sqrt(7840.0 / 5909761.0);
    const auto fs_28_2431_1001 = std::sqrt(5488.0 / 41327.0);
    const auto fs_28_2431_165 = std::sqrt(11760.0 / 537251.0);
    const auto fs_28_2431_2145 = std::sqrt(11760.0 / 41327.0);
    const auto fs_28_2431_2431 = std::sqrt(784.0 / 2431.0);
    const auto fs_28_429_15 = std::sqrt(3920.0 / 61347.0);
    const auto fs_28_429_5 = std::sqrt(3920.0 / 184041.0);
    const auto fs_290_7293_21 = std::sqrt(588700.0 / 17729283.0);
    const auto fs_294_2431_10 = std::sqrt(864360.0 / 5909761.0);
    const auto fs_294_2431_1001 = std::sqrt(605052.0 / 41327.0);
    const auto fs_294_2431_165 = std::sqrt(1296540.0 / 537251.0);
    const auto fs_294_2431_2145 = std::sqrt(1296540.0 / 41327.0);
    const auto fs_294_2431_22 = std::sqrt(172872.0 / 537251.0);
    const auto fs_294_2431_2431 = std::sqrt(86436.0 / 2431.0);
    const auto fs_2_11_2 = std::sqrt(8.0 / 121.0);
    const auto fs_2_33_10 = std::sqrt(40.0 / 1089.0);
    const auto fs_2_33_14 = std::sqrt(56.0 / 1089.0);
    const auto fs_2_33_15 = std::sqrt(20.0 / 363.0);
    const auto fs_2_33_21 = std::sqrt(28.0 / 363.0);
    const auto fs_2_33_3 = std::sqrt(4.0 / 363.0);
    const auto fs_2_39_10 = std::sqrt(40.0 / 1521.0);
    const auto fs_2_39_15 = std::sqrt(20.0 / 507.0);
    const auto fs_2_39_21 = std::sqrt(28.0 / 507.0);
    const auto fs_2_39_30 = std::sqrt(40.0 / 507.0);
    const auto fs_2_39_35 = std::sqrt(140.0 / 1521.0);
    const auto fs_2_39_6 = std::sqrt(8.0 / 507.0);
    const auto fs_2_429_105 = std::sqrt(140.0 / 61347.0);
    const auto fs_2_429_15 = std::sqrt(20.0 / 61347.0);
    const auto fs_2_429_210 = std::sqrt(280.0 / 61347.0);
    const auto fs_300_143_7 = std::sqrt(630000.0 / 20449.0);
    const auto fs_3040_429_21 = std::sqrt(64691200.0 / 61347.0);
    const auto fs_3087_2431_22 = std::sqrt(19059138.0 / 537251.0);
    const auto fs_30_11_21 = std::sqrt(18900.0 / 121.0);
    const auto fs_315_2_3 = std::sqrt(297675.0 / 4.0);
    const auto fs_315_4_15 = std::sqrt(1488375.0 / 16.0);
    const auto fs_315_4_5 = std::sqrt(496125.0 / 16.0);
    const auto fs_315_8_3 = std::sqrt(297675.0 / 64.0);
    const auto fs_320_143_7 = std::sqrt(716800.0 / 20449.0);
    const auto fs_32_429_7 = std::sqrt(7168.0 / 184041.0);
    const auto fs_3325_858_42 = std::sqrt(77389375.0 / 122694.0);
    const auto fs_34_13_5 = std::sqrt(5780.0 / 169.0);
    const auto fs_350_7293_42 = std::sqrt(1715000.0 / 17729283.0);
    const auto fs_35_1_3 = std::sqrt(3675.0);
    const auto fs_35_2431_858 = std::sqrt(7350.0 / 41327.0);
    const auto fs_375_22_10 = std::sqrt(703125.0 / 242.0);
    const auto fs_3800_2431_11 = std::sqrt(14440000.0 / 537251.0);
    const auto fs_380_143_10 = std::sqrt(1444000.0 / 20449.0);
    const auto fs_380_143_165 = std::sqrt(2166000.0 / 1859.0);
    const auto fs_380_143_6 = std::sqrt(866400.0 / 20449.0);
    const auto fs_380_2431_2002 = std::sqrt(2021600.0 / 41327.0);
    const auto fs_380_663_105 = std::sqrt(5054000.0 / 146523.0);
    const auto fs_380_7293_210 = std::sqrt(10108000.0 / 17729283.0);
    const auto fs_3895_858_6 = std::sqrt(15171025.0 / 122694.0);
    const auto fs_38_429_6 = std::sqrt(2888.0 / 61347.0);
    const auto fs_405_8_21 = std::sqrt(3444525.0 / 64.0);
    const auto fs_405_8_35 = std::sqrt(5740875.0 / 64.0);
    const auto fs_4085_286_2 = std::sqrt(16687225.0 / 40898.0);
    const auto fs_4095_16_10 = std::sqrt(83845125.0 / 128.0);
    const auto fs_4095_16_2 = std::sqrt(16769025.0 / 128.0);
    const auto fs_4095_32_10 = std::sqrt(83845125.0 / 512.0);
    const auto fs_4095_8_2 = std::sqrt(16769025.0 / 32.0);
    const auto fs_40_143_21 = std::sqrt(33600.0 / 20449.0);
    const auto fs_410_7293_6 = std::sqrt(336200.0 / 17729283.0);
    const auto fs_420_11_3 = std::sqrt(529200.0 / 121.0);
    const auto fs_42_2431_1001 = std::sqrt(12348.0 / 41327.0);
    const auto fs_42_2431_143 = std::sqrt(1764.0 / 41327.0);
    const auto fs_42_2431_2431 = std::sqrt(1764.0 / 2431.0);
    const auto fs_42_2431_70 = std::sqrt(123480.0 / 5909761.0);
    const auto fs_42_2431_715 = std::sqrt(8820.0 / 41327.0);
    const auto fs_42_2431_770 = std::sqrt(123480.0 / 537251.0);
    const auto fs_430_143_2 = std::sqrt(369800.0 / 20449.0);
    const auto fs_430_2431_2 = std::sqrt(369800.0 / 5909761.0);
    const auto fs_4370_7293_30 = std::sqrt(190969000.0 / 17729283.0);
    const auto fs_4370_7293_70 = std::sqrt(1336783000.0 / 53187849.0);
    const auto fs_43_429_2 = std::sqrt(3698.0 / 184041.0);
    const auto fs_441_2431_1001 = std::sqrt(1361367.0 / 41327.0);
    const auto fs_441_2431_143 = std::sqrt(194481.0 / 41327.0);
    const auto fs_441_2431_2431 = std::sqrt(194481.0 / 2431.0);
    const auto fs_441_2431_70 = std::sqrt(13613670.0 / 5909761.0);
    const auto fs_441_2431_715 = std::sqrt(972405.0 / 41327.0);
    const auto fs_441_2431_770 = std::sqrt(13613670.0 / 537251.0);
    const auto fs_441_4862_10 = std::sqrt(972405.0 / 11819522.0);
    const auto fs_441_4862_110 = std::sqrt(972405.0 / 1074502.0);
    const auto fs_441_4862_1430 = std::sqrt(972405.0 / 82654.0);
    const auto fs_441_4862_4290 = std::sqrt(2917215.0 / 82654.0);
    const auto fs_441_4862_770 = std::sqrt(6806835.0 / 1074502.0);
    const auto fs_450_11_7 = std::sqrt(1417500.0 / 121.0);
    const auto fs_45_1_21 = std::sqrt(42525.0);
    const auto fs_45_1_35 = std::sqrt(70875.0);
    const auto fs_45_2431_66 = std::sqrt(12150.0 / 537251.0);
    const auto fs_45_2_105 = std::sqrt(212625.0 / 4.0);
    const auto fs_45_4_21 = std::sqrt(42525.0 / 16.0);
    const auto fs_45_8_105 = std::sqrt(212625.0 / 64.0);
    const auto fs_45_8_15 = std::sqrt(30375.0 / 64.0);
    const auto fs_45_8_210 = std::sqrt(212625.0 / 32.0);
    const auto fs_4655_858_30 = std::sqrt(108345125.0 / 122694.0);
    const auto fs_4750_2431_7 = std::sqrt(157937500.0 / 5909761.0);
    const auto fs_475_1716_6006 = std::sqrt(1579375.0 / 3432.0);
    const auto fs_475_187_6 = std::sqrt(1353750.0 / 34969.0);
    const auto fs_475_286_165 = std::sqrt(3384375.0 / 7436.0);
    const auto fs_475_44_6 = std::sqrt(676875.0 / 968.0);
    const auto fs_475_7293_6006 = std::sqrt(3158750.0 / 123981.0);
    const auto fs_490_7293_30 = std::sqrt(2401000.0 / 17729283.0);
    const auto fs_495_16_70 = std::sqrt(8575875.0 / 128.0);
    const auto fs_495_8_10 = std::sqrt(1225125.0 / 32.0);
    const auto fs_4_33_6 = std::sqrt(32.0 / 363.0);
    const auto fs_4_39_5 = std::sqrt(80.0 / 1521.0);
    const auto fs_4_429_21 = std::sqrt(112.0 / 61347.0);
    const auto fs_5035_1716_66 = std::sqrt(25351225.0 / 44616.0);
    const auto fs_5035_7293_66 = std::sqrt(50702450.0 / 1611753.0);
    const auto fs_50_1_7 = std::sqrt(17500.0);
    const auto fs_50_2431_165 = std::sqrt(37500.0 / 537251.0);
    const auto fs_5130_2431_6 = std::sqrt(157901400.0 / 5909761.0);
    const auto fs_525_11_3 = std::sqrt(826875.0 / 121.0);
    const auto fs_5320_2431_11 = std::sqrt(28302400.0 / 537251.0);
    const auto fs_5320_7293_14 = std::sqrt(396233600.0 / 53187849.0);
    const auto fs_5320_7293_165 = std::sqrt(141512000.0 / 1611753.0);
    const auto fs_5510_7293_21 = std::sqrt(212520700.0 / 17729283.0);
    const auto fs_55_1_10 = std::sqrt(30250.0);
    const auto fs_55_2_70 = std::sqrt(105875.0 / 2.0);
    const auto fs_560_143_3 = std::sqrt(940800.0 / 20449.0);
    const auto fs_56_2431_210 = std::sqrt(658560.0 / 5909761.0);
    const auto fs_56_2431_30 = std::sqrt(94080.0 / 5909761.0);
    const auto fs_56_2431_429 = std::sqrt(9408.0 / 41327.0);
    const auto fs_56_2431_715 = std::sqrt(15680.0 / 41327.0);
    const auto fs_56_429_3 = std::sqrt(3136.0 / 61347.0);
    const auto fs_585_2_2 = std::sqrt(342225.0 / 2.0);
    const auto fs_585_4_10 = std::sqrt(1711125.0 / 8.0);
    const auto fs_588_2431_210 = std::sqrt(72606240.0 / 5909761.0);
    const auto fs_588_2431_30 = std::sqrt(10372320.0 / 5909761.0);
    const auto fs_588_2431_429 = std::sqrt(1037232.0 / 41327.0);
    const auto fs_588_2431_6 = std::sqrt(2074464.0 / 5909761.0);
    const auto fs_588_2431_715 = std::sqrt(1728720.0 / 41327.0);
    const auto fs_5_143_30 = std::sqrt(750.0 / 20449.0);
    const auto fs_5_1_105 = std::sqrt(2625.0);
    const auto fs_5_1_15 = std::sqrt(375.0);
    const auto fs_5_1_210 = std::sqrt(5250.0);
    const auto fs_5_2431_30030 = std::sqrt(5250.0 / 41327.0);
    const auto fs_600_11_3 = std::sqrt(1080000.0 / 121.0);
    const auto fs_600_143_7 = std::sqrt(2520000.0 / 20449.0);
    const auto fs_60_11_105 = std::sqrt(378000.0 / 121.0);
    const auto fs_6174_2431_6 = std::sqrt(228709656.0 / 5909761.0);
    const auto fs_640_7293_21 = std::sqrt(2867200.0 / 17729283.0);
    const auto fs_645_22_2 = std::sqrt(416025.0 / 242.0);
    const auto fs_65_11_2 = std::sqrt(8450.0 / 121.0);
    const auto fs_65_1_2 = std::sqrt(8450.0);
    const auto fs_65_22_10 = std::sqrt(21125.0 / 242.0);
    const auto fs_65_2_10 = std::sqrt(21125.0 / 2.0);
    const auto fs_65_33_10 = std::sqrt(42250.0 / 1089.0);
    const auto fs_65_33_14 = std::sqrt(59150.0 / 1089.0);
    const auto fs_65_33_15 = std::sqrt(21125.0 / 363.0);
    const auto fs_65_33_21 = std::sqrt(29575.0 / 363.0);
    const auto fs_65_33_3 = std::sqrt(4225.0 / 363.0);
    const auto fs_65_3_10 = std::sqrt(42250.0 / 9.0);
    const auto fs_65_3_14 = std::sqrt(59150.0 / 9.0);
    const auto fs_65_3_15 = std::sqrt(21125.0 / 3.0);
    const auto fs_65_3_21 = std::sqrt(29575.0 / 3.0);
    const auto fs_65_3_3 = std::sqrt(4225.0 / 3.0);
    const auto fs_65_66_2 = std::sqrt(4225.0 / 2178.0);
    const auto fs_65_66_42 = std::sqrt(29575.0 / 726.0);
    const auto fs_65_66_6 = std::sqrt(4225.0 / 726.0);
    const auto fs_65_6_2 = std::sqrt(4225.0 / 18.0);
    const auto fs_65_6_42 = std::sqrt(29575.0 / 6.0);
    const auto fs_65_6_6 = std::sqrt(4225.0 / 6.0);
    const auto fs_6650_7293_42 = std::sqrt(619115000.0 / 17729283.0);
    const auto fs_665_143_14 = std::sqrt(6191150.0 / 20449.0);
    const auto fs_665_2431_858 = std::sqrt(2653350.0 / 41327.0);
    const auto fs_665_286_22 = std::sqrt(442225.0 / 3718.0);
    const auto fs_665_429_6 = std::sqrt(884450.0 / 61347.0);
    const auto fs_665_572_858 = std::sqrt(1326675.0 / 1144.0);
    const auto fs_665_858_858 = std::sqrt(442225.0 / 858.0);
    const auto fs_675_16_30 = std::sqrt(6834375.0 / 128.0);
    const auto fs_675_4_7 = std::sqrt(3189375.0 / 16.0);
    const auto fs_675_8_7 = std::sqrt(3189375.0 / 64.0);
    const auto fs_6_143_21 = std::sqrt(756.0 / 20449.0);
    const auto fs_6_143_35 = std::sqrt(1260.0 / 20449.0);
    const auto fs_700_143_3 = std::sqrt(1470000.0 / 20449.0);
    const auto fs_70_1_15 = std::sqrt(73500.0);
    const auto fs_70_1_5 = std::sqrt(24500.0);
    const auto fs_70_2431_143 = std::sqrt(4900.0 / 41327.0);
    const auto fs_70_2431_22 = std::sqrt(9800.0 / 537251.0);
    const auto fs_70_2431_462 = std::sqrt(205800.0 / 537251.0);
    const auto fs_70_429_3 = std::sqrt(4900.0 / 61347.0);
    const auto fs_70_7293_858 = std::sqrt(9800.0 / 123981.0);
    const auto fs_735_2431_143 = std::sqrt(540225.0 / 41327.0);
    const auto fs_735_2431_462 = std::sqrt(22689450.0 / 537251.0);
    const auto fs_75_11_105 = std::sqrt(590625.0 / 121.0);
    const auto fs_75_1_7 = std::sqrt(39375.0);
    const auto fs_75_2431_30 = std::sqrt(168750.0 / 5909761.0);
    const auto fs_75_2_30 = std::sqrt(84375.0 / 2.0);
    const auto fs_7790_7293_6 = std::sqrt(121368200.0 / 17729283.0);
    const auto fs_7_2431_2 = std::sqrt(98.0 / 5909761.0);
    const auto fs_7_2431_22 = std::sqrt(98.0 / 537251.0);
    const auto fs_7_2431_30030 = std::sqrt(10290.0 / 41327.0);
    const auto fs_7_2431_330 = std::sqrt(1470.0 / 537251.0);
    const auto fs_7_2431_6006 = std::sqrt(2058.0 / 41327.0);
    const auto fs_7_2431_770 = std::sqrt(3430.0 / 537251.0);
    const auto fs_800_143_3 = std::sqrt(1920000.0 / 20449.0);
    const auto fs_80_143_105 = std::sqrt(672000.0 / 20449.0);
    const auto fs_80_1_7 = std::sqrt(44800.0);
    const auto fs_80_2431_10 = std::sqrt(64000.0 / 5909761.0);
    const auto fs_80_2431_165 = std::sqrt(96000.0 / 537251.0);
    const auto fs_80_429_3 = std::sqrt(6400.0 / 61347.0);
    const auto fs_815_7293_6 = std::sqrt(1328450.0 / 17729283.0);
    const auto fs_8170_2431_2 = std::sqrt(133497800.0 / 5909761.0);
    const auto fs_855_2431_66 = std::sqrt(4386150.0 / 537251.0);
    const auto fs_855_572_66 = std::sqrt(2193075.0 / 14872.0);
    const auto fs_855_8_6 = std::sqrt(2193075.0 / 32.0);
    const auto fs_85_2_5 = std::sqrt(36125.0 / 4.0);
    const auto fs_85_4_15 = std::sqrt(108375.0 / 16.0);
    const auto fs_85_4_21 = std::sqrt(151725.0 / 16.0);
    const auto fs_85_4_30 = std::sqrt(108375.0 / 8.0);
    const auto fs_85_4_35 = std::sqrt(252875.0 / 16.0);
    const auto fs_85_4_6 = std::sqrt(21675.0 / 8.0);
    const auto fs_85_8_10 = std::sqrt(36125.0 / 32.0);
    const auto fs_85_8_70 = std::sqrt(252875.0 / 32.0);
    const auto fs_8_429_105 = std::sqrt(2240.0 / 61347.0);
    const auto fs_90_1_7 = std::sqrt(56700.0);
    const auto fs_9310_7293_30 = std::sqrt(866761000.0 / 17729283.0);
    const auto fs_950_143_11 = std::sqrt(902500.0 / 1859.0);
    const auto fs_950_2431_165 = std::sqrt(13537500.0 / 537251.0);
    const auto fs_95_143_2002 = std::sqrt(126350.0 / 143.0);
    const auto fs_95_1_6 = std::sqrt(54150.0);
    const auto fs_95_22_7 = std::sqrt(63175.0 / 484.0);
    const auto fs_95_2431_30030 = std::sqrt(1895250.0 / 41327.0);
    const auto fs_95_286_154 = std::sqrt(63175.0 / 3718.0);
    const auto fs_95_286_70 = std::sqrt(315875.0 / 40898.0);
    const auto fs_95_286_858 = std::sqrt(27075.0 / 286.0);
    const auto fs_95_39_105 = std::sqrt(315875.0 / 507.0);
    const auto fs_95_429_210 = std::sqrt(631750.0 / 61347.0);
    const auto fs_95_572_30030 = std::sqrt(947625.0 / 1144.0);
    const auto fs_95_66_33 = std::sqrt(9025.0 / 132.0);
    const auto fs_980_2431_3 = std::sqrt(2881200.0 / 5909761.0);
    const auto fs_98_2431_165 = std::sqrt(144060.0 / 537251.0);

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ab_2, pc_0 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];

        pc_0[k] = e_0 * fs_4095_32_10 * h1_p1 + e_1 * fs_315_4_15 * h3_p1 - e_1 * fs_4095_16_10 * r_2 * h1_p1 + e_2 * fs_85_4_6 * h5_p1 - e_2 * fs_70_1_15 * r_2 * h3_p1 + e_2 * fs_585_4_10 * r_4 * h1_p1 + e_3 * fs_95_286_70 * h7_p1 - e_3 * fs_255_26_6 * r_2 * h5_p1 + e_3 * fs_210_11_15 * r_4 * h3_p1 - e_3 * fs_65_2_10 * r_6 * h1_p1 + e_4 * fs_147_4862_2 * h9_p1 - e_4 * fs_441_2431_2431 * h9_p9 - e_4 * fs_190_2431_70 * r_2 * h7_p1 + e_4 * fs_17_13_6 * r_4 * h5_p1 - e_4 * fs_280_143_15 * r_6 * h3_p1 + e_4 * fs_65_22_10 * r_8 * h1_p1 - e_5 * fs_7_2431_2 * r_2 * h9_p1 + e_5 * fs_42_2431_2431 * r_2 * h9_p9 + e_5 * fs_10_2431_70 * r_4 * h7_p1 - e_5 * fs_2_39_6 * r_6 * h5_p1 + e_5 * fs_28_429_15 * r_8 * h3_p1 - e_5 * fs_1_11_10 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph9_p8, ab_2, pc_1 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];

        pc_1[k] = e_0 * f_4095_16 * h1_0 + e_1 * f_945_2 * h3_0 - e_1 * f_4095_8 * r_2 * h1_0 + e_2 * f_255_2 * h5_0 - e_2 * f_420_1 * r_2 * h3_0 + e_2 * f_585_2 * r_4 * h1_0 + e_3 * f_1330_143 * h7_0 - e_3 * f_765_13 * r_2 * h5_0 + e_3 * f_1260_11 * r_4 * h3_0 - e_3 * f_65_1 * r_6 * h1_0 + e_4 * f_441_2431 * h9_0 - e_4 * fs_147_2431_12155 * h9_p8 - e_4 * f_5320_2431 * r_2 * h7_0 + e_4 * f_102_13 * r_4 * h5_0 - e_4 * f_1680_143 * r_6 * h3_0 + e_4 * f_65_11 * r_8 * h1_0 - e_5 * f_42_2431 * r_2 * h9_0 + e_5 * fs_14_2431_12155 * r_2 * h9_p8 + e_5 * f_280_2431 * r_4 * h7_0 - e_5 * f_4_13 * r_6 * h5_0 + e_5 * f_56_143 * r_8 * h3_0 - e_5 * f_2_11 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2, pc_2 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_2[k] = - e_0 * fs_1365_32_2 * h1_p1 - e_1 * fs_315_2_3 * h3_p1 + e_1 * fs_1365_16_2 * r_2 * h1_p1 - e_2 * fs_85_4_30 * h5_p1 + e_2 * fs_140_1_3 * r_2 * h3_p1 - e_2 * fs_195_4_2 * r_4 * h1_p1 - e_3 * fs_1330_429_14 * h7_p1 - e_3 * fs_665_858_858 * h7_p7 + e_3 * fs_255_26_30 * r_2 * h5_p1 - e_3 * fs_420_11_3 * r_4 * h3_p1 + e_3 * fs_65_6_2 * r_6 * h1_p1 - e_4 * fs_441_4862_10 * h9_p1 - e_4 * fs_441_2431_715 * h9_p7 + e_4 * fs_5320_7293_14 * r_2 * h7_p1 + e_4 * fs_1330_7293_858 * r_2 * h7_p7 - e_4 * fs_17_13_30 * r_4 * h5_p1 + e_4 * fs_560_143_3 * r_6 * h3_p1 - e_4 * fs_65_66_2 * r_8 * h1_p1 + e_5 * fs_21_2431_10 * r_2 * h9_p1 + e_5 * fs_42_2431_715 * r_2 * h9_p7 - e_5 * fs_280_7293_14 * r_4 * h7_p1 - e_5 * fs_70_7293_858 * r_4 * h7_p7 + e_5 * fs_2_39_30 * r_6 * h5_p1 - e_5 * fs_56_429_3 * r_8 * h3_p1 + e_5 * fs_1_33_2 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ab_2, pc_3 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];

        pc_3[k] = e_1 * fs_315_4_5 * h3_p2 + e_2 * fs_85_4_35 * h5_p2 - e_2 * fs_70_1_5 * r_2 * h3_p2 + e_3 * fs_665_143_14 * h7_p2 - e_3 * fs_95_143_2002 * h7_p6 - e_3 * fs_255_26_35 * r_2 * h5_p2 + e_3 * fs_210_11_5 * r_4 * h3_p2 + e_4 * fs_147_4862_330 * h9_p2 - e_4 * fs_441_4862_1430 * h9_p6 - e_4 * fs_2660_2431_14 * r_2 * h7_p2 + e_4 * fs_380_2431_2002 * r_2 * h7_p6 + e_4 * fs_17_13_35 * r_4 * h5_p2 - e_4 * fs_280_143_5 * r_6 * h3_p2 - e_5 * fs_7_2431_330 * r_2 * h9_p2 + e_5 * fs_21_2431_1430 * r_2 * h9_p6 + e_5 * fs_140_2431_14 * r_4 * h7_p2 - e_5 * fs_20_2431_2002 * r_4 * h7_p6 - e_5 * fs_2_39_35 * r_6 * h5_p2 + e_5 * fs_28_429_5 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_m4, ph5_p3, ph5_p5, ph7_m4, ph7_p3, ph7_p5, ph9_m4, ph9_p3, ph9_p5, ab_2, pc_4, pc_5 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];

        pc_4[k] = - e_1 * fs_45_8_210 * h3_p3 - e_2 * fs_85_4_30 * h5_p3 - e_2 * fs_85_4_6 * h5_p5 + e_2 * fs_5_1_210 * r_2 * h3_p3 - e_3 * f_3325_143 * h7_p3 - e_3 * fs_1330_143_11 * h7_p5 + e_3 * fs_255_26_30 * r_2 * h5_p3 + e_3 * fs_255_26_6 * r_2 * h5_p5 - e_3 * fs_15_11_210 * r_4 * h3_p3 - e_4 * fs_441_4862_110 * h9_p3 - e_4 * fs_147_4862_6006 * h9_p5 + e_4 * f_13300_2431 * r_2 * h7_p3 + e_4 * fs_5320_2431_11 * r_2 * h7_p5 - e_4 * fs_17_13_30 * r_4 * h5_p3 - e_4 * fs_17_13_6 * r_4 * h5_p5 + e_4 * fs_20_143_210 * r_6 * h3_p3 + e_5 * fs_21_2431_110 * r_2 * h9_p3 + e_5 * fs_7_2431_6006 * r_2 * h9_p5 - e_5 * f_700_2431 * r_4 * h7_p3 - e_5 * fs_280_2431_11 * r_4 * h7_p5 + e_5 * fs_2_39_30 * r_6 * h5_p3 + e_5 * fs_2_39_6 * r_6 * h5_p5 - e_5 * fs_2_429_210 * r_8 * h3_p3;

        pc_5[k] = e_2 * f_255_2 * h5_m4 + e_3 * fs_1330_429_165 * h7_m4 - e_3 * f_765_13 * r_2 * h5_m4 + e_4 * fs_441_2431_143 * h9_m4 - e_4 * fs_5320_7293_165 * r_2 * h7_m4 + e_4 * f_102_13 * r_4 * h5_m4 - e_5 * fs_42_2431_143 * r_2 * h9_m4 + e_5 * fs_280_7293_165 * r_4 * h7_m4 - e_5 * f_4_13 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ab_2, pc_6 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];

        pc_6[k] = - e_1 * fs_45_8_210 * h3_m3 + e_2 * fs_85_4_6 * h5_m5 - e_2 * fs_85_4_30 * h5_m3 + e_2 * fs_5_1_210 * r_2 * h3_m3 + e_3 * fs_1330_143_11 * h7_m5 - e_3 * f_3325_143 * h7_m3 - e_3 * fs_255_26_6 * r_2 * h5_m5 + e_3 * fs_255_26_30 * r_2 * h5_m3 - e_3 * fs_15_11_210 * r_4 * h3_m3 + e_4 * fs_147_4862_6006 * h9_m5 - e_4 * fs_441_4862_110 * h9_m3 - e_4 * fs_5320_2431_11 * r_2 * h7_m5 + e_4 * f_13300_2431 * r_2 * h7_m3 + e_4 * fs_17_13_6 * r_4 * h5_m5 - e_4 * fs_17_13_30 * r_4 * h5_m3 + e_4 * fs_20_143_210 * r_6 * h3_m3 - e_5 * fs_7_2431_6006 * r_2 * h9_m5 + e_5 * fs_21_2431_110 * r_2 * h9_m3 + e_5 * fs_280_2431_11 * r_4 * h7_m5 - e_5 * f_700_2431 * r_4 * h7_m3 - e_5 * fs_2_39_6 * r_6 * h5_m5 + e_5 * fs_2_39_30 * r_6 * h5_m3 - e_5 * fs_2_429_210 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ab_2, pc_7 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];

        pc_7[k] = e_1 * fs_315_4_5 * h3_m2 + e_2 * fs_85_4_35 * h5_m2 - e_2 * fs_70_1_5 * r_2 * h3_m2 + e_3 * fs_95_143_2002 * h7_m6 + e_3 * fs_665_143_14 * h7_m2 - e_3 * fs_255_26_35 * r_2 * h5_m2 + e_3 * fs_210_11_5 * r_4 * h3_m2 + e_4 * fs_441_4862_1430 * h9_m6 + e_4 * fs_147_4862_330 * h9_m2 - e_4 * fs_380_2431_2002 * r_2 * h7_m6 - e_4 * fs_2660_2431_14 * r_2 * h7_m2 + e_4 * fs_17_13_35 * r_4 * h5_m2 - e_4 * fs_280_143_5 * r_6 * h3_m2 - e_5 * fs_21_2431_1430 * r_2 * h9_m6 - e_5 * fs_7_2431_330 * r_2 * h9_m2 + e_5 * fs_20_2431_2002 * r_4 * h7_m6 + e_5 * fs_140_2431_14 * r_4 * h7_m2 - e_5 * fs_2_39_35 * r_6 * h5_m2 + e_5 * fs_28_429_5 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m9, ph9_m8, ph9_m7, ph9_m1, ab_2, pc_8, pc_9, pc_10 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];

        pc_8[k] = - e_0 * fs_1365_32_2 * h1_m1 - e_1 * fs_315_2_3 * h3_m1 + e_1 * fs_1365_16_2 * r_2 * h1_m1 - e_2 * fs_85_4_30 * h5_m1 + e_2 * fs_140_1_3 * r_2 * h3_m1 - e_2 * fs_195_4_2 * r_4 * h1_m1 + e_3 * fs_665_858_858 * h7_m7 - e_3 * fs_1330_429_14 * h7_m1 + e_3 * fs_255_26_30 * r_2 * h5_m1 - e_3 * fs_420_11_3 * r_4 * h3_m1 + e_3 * fs_65_6_2 * r_6 * h1_m1 + e_4 * fs_441_2431_715 * h9_m7 - e_4 * fs_441_4862_10 * h9_m1 - e_4 * fs_1330_7293_858 * r_2 * h7_m7 + e_4 * fs_5320_7293_14 * r_2 * h7_m1 - e_4 * fs_17_13_30 * r_4 * h5_m1 + e_4 * fs_560_143_3 * r_6 * h3_m1 - e_4 * fs_65_66_2 * r_8 * h1_m1 - e_5 * fs_42_2431_715 * r_2 * h9_m7 + e_5 * fs_21_2431_10 * r_2 * h9_m1 + e_5 * fs_70_7293_858 * r_4 * h7_m7 - e_5 * fs_280_7293_14 * r_4 * h7_m1 + e_5 * fs_2_39_30 * r_6 * h5_m1 - e_5 * fs_56_429_3 * r_8 * h3_m1 + e_5 * fs_1_33_2 * r_10 * h1_m1;

        pc_9[k] = e_4 * fs_147_2431_12155 * h9_m8 - e_5 * fs_14_2431_12155 * r_2 * h9_m8;

        pc_10[k] = - e_0 * fs_4095_32_10 * h1_m1 - e_1 * fs_315_4_15 * h3_m1 + e_1 * fs_4095_16_10 * r_2 * h1_m1 - e_2 * fs_85_4_6 * h5_m1 + e_2 * fs_70_1_15 * r_2 * h3_m1 - e_2 * fs_585_4_10 * r_4 * h1_m1 - e_3 * fs_95_286_70 * h7_m1 + e_3 * fs_255_26_6 * r_2 * h5_m1 - e_3 * fs_210_11_15 * r_4 * h3_m1 + e_3 * fs_65_2_10 * r_6 * h1_m1 + e_4 * fs_441_2431_2431 * h9_m9 - e_4 * fs_147_4862_2 * h9_m1 + e_4 * fs_190_2431_70 * r_2 * h7_m1 - e_4 * fs_17_13_6 * r_4 * h5_m1 + e_4 * fs_280_143_15 * r_6 * h3_m1 - e_4 * fs_65_22_10 * r_8 * h1_m1 - e_5 * fs_42_2431_2431 * r_2 * h9_m9 + e_5 * fs_7_2431_2 * r_2 * h9_m1 - e_5 * fs_10_2431_70 * r_4 * h7_m1 + e_5 * fs_2_39_6 * r_6 * h5_m1 - e_5 * fs_28_429_15 * r_8 * h3_m1 + e_5 * fs_1_11_10 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ab_2, pc_11 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];

        pc_11[k] = - e_1 * fs_1575_8_3 * h3_p2 - e_2 * fs_85_4_21 * h5_p2 + e_2 * fs_175_1_3 * r_2 * h3_p2 - e_3 * fs_285_572_210 * h7_p2 + e_3 * fs_255_26_21 * r_2 * h5_p2 - e_3 * fs_525_11_3 * r_4 * h3_p2 - e_4 * fs_147_4862_22 * h9_p2 - e_4 * fs_294_2431_2431 * h9_p8 + e_4 * fs_285_2431_210 * r_2 * h7_p2 - e_4 * fs_17_13_21 * r_4 * h5_p2 + e_4 * fs_700_143_3 * r_6 * h3_p2 + e_5 * fs_7_2431_22 * r_2 * h9_p2 + e_5 * fs_28_2431_2431 * r_2 * h9_p8 - e_5 * fs_15_2431_210 * r_4 * h7_p2 + e_5 * fs_2_39_21 * r_6 * h5_p2 - e_5 * fs_70_429_3 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2, pc_12 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_12[k] = e_0 * fs_4095_16_2 * h1_p1 - e_1 * fs_315_8_3 * h3_p1 - e_1 * fs_4095_8_2 * r_2 * h1_p1 - e_2 * fs_85_4_30 * h5_p1 + e_2 * fs_35_1_3 * r_2 * h3_p1 + e_2 * fs_585_2_2 * r_4 * h1_p1 - e_3 * fs_2185_572_14 * h7_p1 + e_3 * fs_665_572_858 * h7_p7 + e_3 * fs_255_26_30 * r_2 * h5_p1 - e_3 * fs_105_11_3 * r_4 * h3_p1 - e_3 * fs_65_1_2 * r_6 * h1_p1 - e_4 * fs_294_2431_10 * h9_p1 - e_4 * fs_588_2431_715 * h9_p7 + e_4 * fs_2185_2431_14 * r_2 * h7_p1 - e_4 * fs_665_2431_858 * r_2 * h7_p7 - e_4 * fs_17_13_30 * r_4 * h5_p1 + e_4 * fs_140_143_3 * r_6 * h3_p1 + e_4 * fs_65_11_2 * r_8 * h1_p1 + e_5 * fs_28_2431_10 * r_2 * h9_p1 + e_5 * fs_56_2431_715 * r_2 * h9_p7 - e_5 * fs_115_2431_14 * r_4 * h7_p1 + e_5 * fs_35_2431_858 * r_4 * h7_p7 + e_5 * fs_2_39_30 * r_6 * h5_p1 - e_5 * fs_14_429_3 * r_8 * h3_p1 - e_5 * fs_2_11_2 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ab_2, pc_13 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];

        pc_13[k] = e_0 * f_1365_4 * h1_0 + e_1 * f_945_4 * h3_0 - e_1 * f_1365_2 * r_2 * h1_0 - e_2 * f_255_2 * h5_0 - e_2 * f_210_1 * r_2 * h3_0 + e_2 * f_390_1 * r_4 * h1_0 - e_3 * f_12635_429 * h7_0 + e_3 * fs_475_1716_6006 * h7_p6 + e_3 * f_765_13 * r_2 * h5_0 + e_3 * f_630_11 * r_4 * h3_0 - e_3 * f_260_3 * r_6 * h1_0 - e_4 * f_2646_2431 * h9_0 - e_4 * fs_441_4862_4290 * h9_p6 + e_4 * f_50540_7293 * r_2 * h7_0 - e_4 * fs_475_7293_6006 * r_2 * h7_p6 - e_4 * f_102_13 * r_4 * h5_0 - e_4 * f_840_143 * r_6 * h3_0 + e_4 * f_260_33 * r_8 * h1_0 + e_5 * f_252_2431 * r_2 * h9_0 + e_5 * fs_21_2431_4290 * r_2 * h9_p6 - e_5 * f_2660_7293 * r_4 * h7_0 + e_5 * fs_25_7293_6006 * r_4 * h7_p6 + e_5 * f_4_13 * r_6 * h5_0 + e_5 * f_28_143 * r_8 * h3_0 - e_5 * f_8_33 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_14 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_14[k] = - e_0 * fs_1365_32_6 * h1_p1 - e_1 * f_2205_8 * h3_p1 + e_1 * fs_1365_16_6 * r_2 * h1_p1 + e_2 * fs_85_8_10 * h5_p1 + e_2 * fs_85_4_21 * h5_p5 + e_2 * f_245_1 * r_2 * h3_p1 - e_2 * fs_195_4_6 * r_4 * h1_p1 + e_3 * fs_3325_858_42 * h7_p1 + e_3 * fs_95_286_154 * h7_p5 - e_3 * fs_255_52_10 * r_2 * h5_p1 - e_3 * fs_255_26_21 * r_2 * h5_p5 - e_3 * f_735_11 * r_4 * h3_p1 + e_3 * fs_65_6_6 * r_6 * h1_p1 + e_4 * fs_588_2431_30 * h9_p1 - e_4 * fs_588_2431_429 * h9_p5 - e_4 * fs_6650_7293_42 * r_2 * h7_p1 - e_4 * fs_190_2431_154 * r_2 * h7_p5 + e_4 * fs_17_26_10 * r_4 * h5_p1 + e_4 * fs_17_13_21 * r_4 * h5_p5 + e_4 * f_980_143 * r_6 * h3_p1 - e_4 * fs_65_66_6 * r_8 * h1_p1 - e_5 * fs_56_2431_30 * r_2 * h9_p1 + e_5 * fs_56_2431_429 * r_2 * h9_p5 + e_5 * fs_350_7293_42 * r_4 * h7_p1 + e_5 * fs_10_2431_154 * r_4 * h7_p5 - e_5 * fs_1_39_10 * r_6 * h5_p1 - e_5 * fs_2_39_21 * r_6 * h5_p5 - e_5 * f_98_429 * r_8 * h3_p1 + e_5 * fs_1_33_6 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph3_p2, ph5_m3, ph5_p2, ph5_p4, ph7_m3, ph7_p2, ph7_p4, ph9_m3, ph9_p2, ph9_p4, ab_2, pc_15, pc_16 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_15[k] = e_1 * fs_495_16_70 * h3_p2 + e_2 * fs_85_8_10 * h5_p2 + e_2 * fs_85_4_30 * h5_p4 - e_2 * fs_55_2_70 * r_2 * h3_p2 - e_3 * f_665_26 * h7_p2 - e_3 * fs_665_286_22 * h7_p4 - e_3 * fs_255_52_10 * r_2 * h5_p2 - e_3 * fs_255_26_30 * r_2 * h5_p4 + e_3 * fs_15_2_70 * r_4 * h3_p2 - e_4 * fs_147_2431_1155 * h9_p2 - e_4 * fs_147_2431_4290 * h9_p4 + e_4 * f_1330_221 * r_2 * h7_p2 + e_4 * fs_1330_2431_22 * r_2 * h7_p4 + e_4 * fs_17_26_10 * r_4 * h5_p2 + e_4 * fs_17_13_30 * r_4 * h5_p4 - e_4 * fs_10_13_70 * r_6 * h3_p2 + e_5 * fs_14_2431_1155 * r_2 * h9_p2 + e_5 * fs_14_2431_4290 * r_2 * h9_p4 - e_5 * f_70_221 * r_4 * h7_p2 - e_5 * fs_70_2431_22 * r_4 * h7_p4 - e_5 * fs_1_39_10 * r_6 * h5_p2 - e_5 * fs_2_39_30 * r_6 * h5_p4 + e_5 * fs_1_39_70 * r_8 * h3_p2;

        pc_16[k] = - e_1 * fs_675_8_7 * h3_m3 - e_2 * f_255_2 * h5_m3 + e_2 * fs_75_1_7 * r_2 * h3_m3 + e_3 * fs_4655_858_30 * h7_m3 + e_3 * f_765_13 * r_2 * h5_m3 - e_3 * fs_225_11_7 * r_4 * h3_m3 + e_4 * fs_1764_2431_33 * h9_m3 - e_4 * fs_9310_7293_30 * r_2 * h7_m3 - e_4 * f_102_13 * r_4 * h5_m3 + e_4 * fs_300_143_7 * r_6 * h3_m3 - e_5 * fs_168_2431_33 * r_2 * h9_m3 + e_5 * fs_490_7293_30 * r_4 * h7_m3 + e_5 * f_4_13 * r_6 * h5_m3 - e_5 * fs_10_143_7 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_17 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];

        pc_17[k] = e_1 * fs_495_16_70 * h3_m2 - e_2 * fs_85_4_30 * h5_m4 + e_2 * fs_85_8_10 * h5_m2 - e_2 * fs_55_2_70 * r_2 * h3_m2 + e_3 * fs_665_286_22 * h7_m4 - e_3 * f_665_26 * h7_m2 + e_3 * fs_255_26_30 * r_2 * h5_m4 - e_3 * fs_255_52_10 * r_2 * h5_m2 + e_3 * fs_15_2_70 * r_4 * h3_m2 + e_4 * fs_147_2431_4290 * h9_m4 - e_4 * fs_147_2431_1155 * h9_m2 - e_4 * fs_1330_2431_22 * r_2 * h7_m4 + e_4 * f_1330_221 * r_2 * h7_m2 - e_4 * fs_17_13_30 * r_4 * h5_m4 + e_4 * fs_17_26_10 * r_4 * h5_m2 - e_4 * fs_10_13_70 * r_6 * h3_m2 - e_5 * fs_14_2431_4290 * r_2 * h9_m4 + e_5 * fs_14_2431_1155 * r_2 * h9_m2 + e_5 * fs_70_2431_22 * r_4 * h7_m4 - e_5 * f_70_221 * r_4 * h7_m2 + e_5 * fs_2_39_30 * r_6 * h5_m4 - e_5 * fs_1_39_10 * r_6 * h5_m2 + e_5 * fs_1_39_70 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m6, ph7_m5, ph7_m1, ph9_m6, ph9_m5, ph9_m1, ab_2, pc_18, pc_19 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_18[k] = - e_0 * fs_1365_32_6 * h1_m1 - e_1 * f_2205_8 * h3_m1 + e_1 * fs_1365_16_6 * r_2 * h1_m1 - e_2 * fs_85_4_21 * h5_m5 + e_2 * fs_85_8_10 * h5_m1 + e_2 * f_245_1 * r_2 * h3_m1 - e_2 * fs_195_4_6 * r_4 * h1_m1 - e_3 * fs_95_286_154 * h7_m5 + e_3 * fs_3325_858_42 * h7_m1 + e_3 * fs_255_26_21 * r_2 * h5_m5 - e_3 * fs_255_52_10 * r_2 * h5_m1 - e_3 * f_735_11 * r_4 * h3_m1 + e_3 * fs_65_6_6 * r_6 * h1_m1 + e_4 * fs_588_2431_429 * h9_m5 + e_4 * fs_588_2431_30 * h9_m1 + e_4 * fs_190_2431_154 * r_2 * h7_m5 - e_4 * fs_6650_7293_42 * r_2 * h7_m1 - e_4 * fs_17_13_21 * r_4 * h5_m5 + e_4 * fs_17_26_10 * r_4 * h5_m1 + e_4 * f_980_143 * r_6 * h3_m1 - e_4 * fs_65_66_6 * r_8 * h1_m1 - e_5 * fs_56_2431_429 * r_2 * h9_m5 - e_5 * fs_56_2431_30 * r_2 * h9_m1 - e_5 * fs_10_2431_154 * r_4 * h7_m5 + e_5 * fs_350_7293_42 * r_4 * h7_m1 + e_5 * fs_2_39_21 * r_6 * h5_m5 - e_5 * fs_1_39_10 * r_6 * h5_m1 - e_5 * f_98_429 * r_8 * h3_m1 + e_5 * fs_1_33_6 * r_10 * h1_m1;

        pc_19[k] = - e_3 * fs_475_1716_6006 * h7_m6 + e_4 * fs_441_4862_4290 * h9_m6 + e_4 * fs_475_7293_6006 * r_2 * h7_m6 - e_5 * fs_21_2431_4290 * r_2 * h9_m6 - e_5 * fs_25_7293_6006 * r_4 * h7_m6;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ab_2, pc_20 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];

        pc_20[k] = - e_0 * fs_4095_16_2 * h1_m1 + e_1 * fs_315_8_3 * h3_m1 + e_1 * fs_4095_8_2 * r_2 * h1_m1 + e_2 * fs_85_4_30 * h5_m1 - e_2 * fs_35_1_3 * r_2 * h3_m1 - e_2 * fs_585_2_2 * r_4 * h1_m1 - e_3 * fs_665_572_858 * h7_m7 + e_3 * fs_2185_572_14 * h7_m1 - e_3 * fs_255_26_30 * r_2 * h5_m1 + e_3 * fs_105_11_3 * r_4 * h3_m1 + e_3 * fs_65_1_2 * r_6 * h1_m1 + e_4 * fs_588_2431_715 * h9_m7 + e_4 * fs_294_2431_10 * h9_m1 + e_4 * fs_665_2431_858 * r_2 * h7_m7 - e_4 * fs_2185_2431_14 * r_2 * h7_m1 + e_4 * fs_17_13_30 * r_4 * h5_m1 - e_4 * fs_140_143_3 * r_6 * h3_m1 - e_4 * fs_65_11_2 * r_8 * h1_m1 - e_5 * fs_56_2431_715 * r_2 * h9_m7 - e_5 * fs_28_2431_10 * r_2 * h9_m1 - e_5 * fs_35_2431_858 * r_4 * h7_m7 + e_5 * fs_115_2431_14 * r_4 * h7_m1 - e_5 * fs_2_39_30 * r_6 * h5_m1 + e_5 * fs_14_429_3 * r_8 * h3_m1 + e_5 * fs_2_11_2 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph3_p3, ph5_m2, ph5_p3, ph7_m2, ph7_p3, ph7_p7, ph9_m8, ph9_m2, ph9_p3, ph9_p7, ab_2, pc_21, pc_22 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];

        pc_21[k] = e_1 * fs_1575_8_3 * h3_m2 + e_2 * fs_85_4_21 * h5_m2 - e_2 * fs_175_1_3 * r_2 * h3_m2 + e_3 * fs_285_572_210 * h7_m2 - e_3 * fs_255_26_21 * r_2 * h5_m2 + e_3 * fs_525_11_3 * r_4 * h3_m2 + e_4 * fs_294_2431_2431 * h9_m8 + e_4 * fs_147_4862_22 * h9_m2 - e_4 * fs_285_2431_210 * r_2 * h7_m2 + e_4 * fs_17_13_21 * r_4 * h5_m2 - e_4 * fs_700_143_3 * r_6 * h3_m2 - e_5 * fs_28_2431_2431 * r_2 * h9_m8 - e_5 * fs_7_2431_22 * r_2 * h9_m2 + e_5 * fs_15_2431_210 * r_4 * h7_m2 - e_5 * fs_2_39_21 * r_6 * h5_m2 + e_5 * fs_70_429_3 * r_8 * h3_m2;

        pc_22[k] = e_1 * fs_675_8_7 * h3_p3 + e_2 * f_255_2 * h5_p3 - e_2 * fs_75_1_7 * r_2 * h3_p3 + e_3 * fs_1425_572_30 * h7_p3 - e_3 * fs_95_572_30030 * h7_p7 - e_3 * f_765_13 * r_2 * h5_p3 + e_3 * fs_225_11_7 * r_4 * h3_p3 + e_4 * fs_147_2431_33 * h9_p3 - e_4 * fs_294_2431_1001 * h9_p7 - e_4 * fs_1425_2431_30 * r_2 * h7_p3 + e_4 * fs_95_2431_30030 * r_2 * h7_p7 + e_4 * f_102_13 * r_4 * h5_p3 - e_4 * fs_300_143_7 * r_6 * h3_p3 - e_5 * fs_14_2431_33 * r_2 * h9_p3 + e_5 * fs_28_2431_1001 * r_2 * h9_p7 + e_5 * fs_75_2431_30 * r_4 * h7_p3 - e_5 * fs_5_2431_30030 * r_4 * h7_p7 - e_5 * f_4_13 * r_6 * h5_p3 + e_5 * fs_10_143_7 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ab_2, pc_23 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];

        pc_23[k] = - e_1 * fs_45_2_105 * h3_p2 + e_2 * fs_85_4_15 * h5_p2 + e_2 * fs_20_1_105 * r_2 * h3_p2 + e_3 * fs_2565_286_6 * h7_p2 + e_3 * fs_95_286_858 * h7_p6 - e_3 * fs_255_26_15 * r_2 * h5_p2 - e_3 * fs_60_11_105 * r_4 * h3_p2 + e_4 * fs_147_4862_770 * h9_p2 - e_4 * fs_147_4862_30030 * h9_p6 - e_4 * fs_5130_2431_6 * r_2 * h7_p2 - e_4 * fs_190_2431_858 * r_2 * h7_p6 + e_4 * fs_17_13_15 * r_4 * h5_p2 + e_4 * fs_80_143_105 * r_6 * h3_p2 - e_5 * fs_7_2431_770 * r_2 * h9_p2 + e_5 * fs_7_2431_30030 * r_2 * h9_p6 + e_5 * fs_270_2431_6 * r_4 * h7_p2 + e_5 * fs_10_2431_858 * r_4 * h7_p6 - e_5 * fs_2_39_15 * r_6 * h5_p2 - e_5 * fs_8_429_105 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_24 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_24[k] = e_0 * fs_1365_16_14 * h1_p1 - e_1 * fs_405_8_21 * h3_p1 - e_1 * fs_1365_8_14 * r_2 * h1_p1 - e_2 * f_255_2 * h5_p5 + e_2 * fs_45_1_21 * r_2 * h3_p1 + e_2 * fs_195_2_14 * r_4 * h1_p1 + e_3 * fs_2375_132_2 * h7_p1 + e_3 * fs_5035_1716_66 * h7_p5 + e_3 * f_765_13 * r_2 * h5_p5 - e_3 * fs_135_11_21 * r_4 * h3_p1 - e_3 * fs_65_3_14 * r_6 * h1_p1 + e_4 * fs_441_2431_70 * h9_p1 - e_4 * fs_441_2431_1001 * h9_p5 - e_4 * fs_2375_561_2 * r_2 * h7_p1 - e_4 * fs_5035_7293_66 * r_2 * h7_p5 - e_4 * f_102_13 * r_4 * h5_p5 + e_4 * fs_180_143_21 * r_6 * h3_p1 + e_4 * fs_65_33_14 * r_8 * h1_p1 - e_5 * fs_42_2431_70 * r_2 * h9_p1 + e_5 * fs_42_2431_1001 * r_2 * h9_p5 + e_5 * fs_125_561_2 * r_4 * h7_p1 + e_5 * fs_265_7293_66 * r_4 * h7_p5 + e_5 * f_4_13 * r_6 * h5_p5 - e_5 * fs_6_143_21 * r_8 * h3_p1 - e_5 * fs_2_33_14 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2, pc_25 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_25[k] = e_0 * fs_1365_16_21 * h1_0 - e_1 * fs_45_4_21 * h3_0 - e_1 * fs_1365_8_21 * r_2 * h1_0 - e_2 * fs_85_4_21 * h5_0 - e_2 * fs_85_4_15 * h5_p4 + e_2 * fs_10_1_21 * r_2 * h3_0 + e_2 * fs_195_2_21 * r_4 * h1_0 + e_3 * fs_3040_429_21 * h7_0 + e_3 * fs_950_143_11 * h7_p4 + e_3 * fs_255_26_21 * r_2 * h5_0 + e_3 * fs_255_26_15 * r_2 * h5_p4 - e_3 * fs_30_11_21 * r_4 * h3_0 - e_3 * fs_65_3_21 * r_6 * h1_0 + e_4 * fs_1764_2431_21 * h9_0 - e_4 * fs_294_2431_2145 * h9_p4 - e_4 * fs_12160_7293_21 * r_2 * h7_0 - e_4 * fs_3800_2431_11 * r_2 * h7_p4 - e_4 * fs_17_13_21 * r_4 * h5_0 - e_4 * fs_17_13_15 * r_4 * h5_p4 + e_4 * fs_40_143_21 * r_6 * h3_0 + e_4 * fs_65_33_21 * r_8 * h1_0 - e_5 * fs_168_2431_21 * r_2 * h9_0 + e_5 * fs_28_2431_2145 * r_2 * h9_p4 + e_5 * fs_640_7293_21 * r_4 * h7_0 + e_5 * fs_200_2431_11 * r_4 * h7_p4 + e_5 * fs_2_39_21 * r_6 * h5_0 + e_5 * fs_2_39_15 * r_6 * h5_p4 - e_5 * fs_4_429_21 * r_8 * h3_0 - e_5 * fs_2_33_21 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_26 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_26[k] = - e_0 * fs_1365_16_3 * h1_p1 - e_1 * fs_1935_16_2 * h3_p1 - e_1 * fs_675_16_30 * h3_p3 + e_1 * fs_1365_8_3 * r_2 * h1_p1 + e_2 * fs_85_2_5 * h5_p1 + e_2 * fs_215_2_2 * r_2 * h3_p1 + e_2 * fs_75_2_30 * r_2 * h3_p3 - e_2 * fs_195_2_3 * r_4 * h1_p1 - e_3 * fs_2755_858_21 * h7_p1 + e_3 * fs_95_22_7 * h7_p3 - e_3 * fs_255_13_5 * r_2 * h5_p1 - e_3 * fs_645_22_2 * r_4 * h3_p1 - e_3 * fs_225_22_30 * r_4 * h3_p3 + e_3 * fs_65_3_3 * r_6 * h1_p1 - e_4 * fs_2058_2431_15 * h9_p1 - e_4 * fs_441_2431_770 * h9_p3 + e_4 * fs_5510_7293_21 * r_2 * h7_p1 - e_4 * fs_190_187_7 * r_2 * h7_p3 + e_4 * fs_34_13_5 * r_4 * h5_p1 + e_4 * fs_430_143_2 * r_6 * h3_p1 + e_4 * fs_150_143_30 * r_6 * h3_p3 - e_4 * fs_65_33_3 * r_8 * h1_p1 + e_5 * fs_196_2431_15 * r_2 * h9_p1 + e_5 * fs_42_2431_770 * r_2 * h9_p3 - e_5 * fs_290_7293_21 * r_4 * h7_p1 + e_5 * fs_10_187_7 * r_4 * h7_p3 - e_5 * fs_4_39_5 * r_6 * h5_p1 - e_5 * fs_43_429_2 * r_8 * h3_p1 - e_5 * fs_5_143_30 * r_8 * h3_p3 + e_5 * fs_2_33_3 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ab_2, pc_27 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];

        pc_27[k] = e_1 * fs_225_1_3 * h3_m2 - e_2 * fs_85_4_21 * h5_m2 - e_2 * fs_200_1_3 * r_2 * h3_m2 + e_3 * fs_95_429_210 * h7_m2 + e_3 * fs_255_26_21 * r_2 * h5_m2 + e_3 * fs_600_11_3 * r_4 * h3_m2 + e_4 * fs_3087_2431_22 * h9_m2 - e_4 * fs_380_7293_210 * r_2 * h7_m2 - e_4 * fs_17_13_21 * r_4 * h5_m2 - e_4 * fs_800_143_3 * r_6 * h3_m2 - e_5 * fs_294_2431_22 * r_2 * h9_m2 + e_5 * fs_20_7293_210 * r_4 * h7_m2 + e_5 * fs_2_39_21 * r_6 * h5_m2 + e_5 * fs_80_429_3 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m4, ph5_m1, ph7_m4, ph7_m3, ph7_m1, ph9_m4, ph9_m3, ph9_m1, ab_2, pc_28, pc_29 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_28[k] = - e_0 * fs_1365_16_3 * h1_m1 + e_1 * fs_675_16_30 * h3_m3 - e_1 * fs_1935_16_2 * h3_m1 + e_1 * fs_1365_8_3 * r_2 * h1_m1 + e_2 * fs_85_2_5 * h5_m1 - e_2 * fs_75_2_30 * r_2 * h3_m3 + e_2 * fs_215_2_2 * r_2 * h3_m1 - e_2 * fs_195_2_3 * r_4 * h1_m1 - e_3 * fs_95_22_7 * h7_m3 - e_3 * fs_2755_858_21 * h7_m1 - e_3 * fs_255_13_5 * r_2 * h5_m1 + e_3 * fs_225_22_30 * r_4 * h3_m3 - e_3 * fs_645_22_2 * r_4 * h3_m1 + e_3 * fs_65_3_3 * r_6 * h1_m1 + e_4 * fs_441_2431_770 * h9_m3 - e_4 * fs_2058_2431_15 * h9_m1 + e_4 * fs_190_187_7 * r_2 * h7_m3 + e_4 * fs_5510_7293_21 * r_2 * h7_m1 + e_4 * fs_34_13_5 * r_4 * h5_m1 - e_4 * fs_150_143_30 * r_6 * h3_m3 + e_4 * fs_430_143_2 * r_6 * h3_m1 - e_4 * fs_65_33_3 * r_8 * h1_m1 - e_5 * fs_42_2431_770 * r_2 * h9_m3 + e_5 * fs_196_2431_15 * r_2 * h9_m1 - e_5 * fs_10_187_7 * r_4 * h7_m3 - e_5 * fs_290_7293_21 * r_4 * h7_m1 - e_5 * fs_4_39_5 * r_6 * h5_m1 + e_5 * fs_5_143_30 * r_8 * h3_m3 - e_5 * fs_43_429_2 * r_8 * h3_m1 + e_5 * fs_2_33_3 * r_10 * h1_m1;

        pc_29[k] = e_2 * fs_85_4_15 * h5_m4 - e_3 * fs_950_143_11 * h7_m4 - e_3 * fs_255_26_15 * r_2 * h5_m4 + e_4 * fs_294_2431_2145 * h9_m4 + e_4 * fs_3800_2431_11 * r_2 * h7_m4 + e_4 * fs_17_13_15 * r_4 * h5_m4 - e_5 * fs_28_2431_2145 * r_2 * h9_m4 - e_5 * fs_200_2431_11 * r_4 * h7_m4 - e_5 * fs_2_39_15 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ab_2, pc_30 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_30[k] = - e_0 * fs_1365_16_14 * h1_m1 + e_1 * fs_405_8_21 * h3_m1 + e_1 * fs_1365_8_14 * r_2 * h1_m1 + e_2 * f_255_2 * h5_m5 - e_2 * fs_45_1_21 * r_2 * h3_m1 - e_2 * fs_195_2_14 * r_4 * h1_m1 - e_3 * fs_5035_1716_66 * h7_m5 - e_3 * fs_2375_132_2 * h7_m1 - e_3 * f_765_13 * r_2 * h5_m5 + e_3 * fs_135_11_21 * r_4 * h3_m1 + e_3 * fs_65_3_14 * r_6 * h1_m1 + e_4 * fs_441_2431_1001 * h9_m5 - e_4 * fs_441_2431_70 * h9_m1 + e_4 * fs_5035_7293_66 * r_2 * h7_m5 + e_4 * fs_2375_561_2 * r_2 * h7_m1 + e_4 * f_102_13 * r_4 * h5_m5 - e_4 * fs_180_143_21 * r_6 * h3_m1 - e_4 * fs_65_33_14 * r_8 * h1_m1 - e_5 * fs_42_2431_1001 * r_2 * h9_m5 + e_5 * fs_42_2431_70 * r_2 * h9_m1 - e_5 * fs_265_7293_66 * r_4 * h7_m5 - e_5 * fs_125_561_2 * r_4 * h7_m1 - e_5 * f_4_13 * r_6 * h5_m5 + e_5 * fs_6_143_21 * r_8 * h3_m1 + e_5 * fs_2_33_14 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m7, ph7_m6, ph7_m3, ph7_m2, ph9_m7, ph9_m6, ph9_m3, ph9_m2, ab_2, pc_31, pc_32 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_31[k] = e_1 * fs_45_2_105 * h3_m2 - e_2 * fs_85_4_15 * h5_m2 - e_2 * fs_20_1_105 * r_2 * h3_m2 - e_3 * fs_95_286_858 * h7_m6 - e_3 * fs_2565_286_6 * h7_m2 + e_3 * fs_255_26_15 * r_2 * h5_m2 + e_3 * fs_60_11_105 * r_4 * h3_m2 + e_4 * fs_147_4862_30030 * h9_m6 - e_4 * fs_147_4862_770 * h9_m2 + e_4 * fs_190_2431_858 * r_2 * h7_m6 + e_4 * fs_5130_2431_6 * r_2 * h7_m2 - e_4 * fs_17_13_15 * r_4 * h5_m2 - e_4 * fs_80_143_105 * r_6 * h3_m2 - e_5 * fs_7_2431_30030 * r_2 * h9_m6 + e_5 * fs_7_2431_770 * r_2 * h9_m2 - e_5 * fs_10_2431_858 * r_4 * h7_m6 - e_5 * fs_270_2431_6 * r_4 * h7_m2 + e_5 * fs_2_39_15 * r_6 * h5_m2 + e_5 * fs_8_429_105 * r_8 * h3_m2;

        pc_32[k] = - e_1 * fs_675_8_7 * h3_m3 - e_2 * f_255_2 * h5_m3 + e_2 * fs_75_1_7 * r_2 * h3_m3 + e_3 * fs_95_572_30030 * h7_m7 - e_3 * fs_1425_572_30 * h7_m3 + e_3 * f_765_13 * r_2 * h5_m3 - e_3 * fs_225_11_7 * r_4 * h3_m3 + e_4 * fs_294_2431_1001 * h9_m7 - e_4 * fs_147_2431_33 * h9_m3 - e_4 * fs_95_2431_30030 * r_2 * h7_m7 + e_4 * fs_1425_2431_30 * r_2 * h7_m3 - e_4 * f_102_13 * r_4 * h5_m3 + e_4 * fs_300_143_7 * r_6 * h3_m3 - e_5 * fs_28_2431_1001 * r_2 * h9_m7 + e_5 * fs_14_2431_33 * r_2 * h9_m3 + e_5 * fs_5_2431_30030 * r_4 * h7_m7 - e_5 * fs_75_2431_30 * r_4 * h7_m3 + e_5 * f_4_13 * r_6 * h5_m3 - e_5 * fs_10_143_7 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p4, ph5_p5, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ph9_p3, ph9_p4, ph9_p5, ph9_p6, ab_2, pc_33, pc_34 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_33[k] = - e_2 * f_255_2 * h5_p4 - e_3 * fs_475_286_165 * h7_p4 - e_3 * fs_285_572_4290 * h7_p6 + e_3 * f_765_13 * r_2 * h5_p4 - e_4 * fs_147_2431_143 * h9_p4 - e_4 * fs_147_4862_6006 * h9_p6 + e_4 * fs_950_2431_165 * r_2 * h7_p4 + e_4 * fs_285_2431_4290 * r_2 * h7_p6 - e_4 * f_102_13 * r_4 * h5_p4 + e_5 * fs_14_2431_143 * r_2 * h9_p4 + e_5 * fs_7_2431_6006 * r_2 * h9_p6 - e_5 * fs_50_2431_165 * r_4 * h7_p4 - e_5 * fs_15_2431_4290 * r_4 * h7_p6 + e_5 * f_4_13 * r_6 * h5_p4;

        pc_34[k] = e_1 * fs_405_8_35 * h3_p3 + e_2 * f_255_2 * h5_p5 - e_2 * fs_45_1_35 * r_2 * h3_p3 - e_3 * fs_475_44_6 * h7_p3 - e_3 * fs_855_572_66 * h7_p5 - e_3 * f_765_13 * r_2 * h5_p5 + e_3 * fs_135_11_35 * r_4 * h3_p3 - e_4 * fs_294_2431_165 * h9_p3 - e_4 * fs_294_2431_1001 * h9_p5 + e_4 * fs_475_187_6 * r_2 * h7_p3 + e_4 * fs_855_2431_66 * r_2 * h7_p5 + e_4 * f_102_13 * r_4 * h5_p5 - e_4 * fs_180_143_35 * r_6 * h3_p3 + e_5 * fs_28_2431_165 * r_2 * h9_p3 + e_5 * fs_28_2431_1001 * r_2 * h9_p5 - e_5 * fs_25_187_6 * r_4 * h7_p3 - e_5 * fs_45_2431_66 * r_4 * h7_p5 - e_5 * f_4_13 * r_6 * h5_p5 + e_5 * fs_6_143_35 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_35 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_35[k] = - e_1 * fs_45_8_105 * h3_p2 + e_2 * fs_85_4_15 * h5_p2 + e_2 * fs_5_1_105 * r_2 * h3_p2 - e_3 * fs_15485_1716_6 * h7_p2 + e_3 * fs_95_66_33 * h7_p4 - e_3 * fs_255_26_15 * r_2 * h5_p2 - e_3 * fs_15_11_105 * r_4 * h3_p2 - e_4 * fs_441_4862_770 * h9_p2 - e_4 * fs_441_2431_715 * h9_p4 + e_4 * fs_15485_7293_6 * r_2 * h7_p2 - e_4 * fs_190_561_33 * r_2 * h7_p4 + e_4 * fs_17_13_15 * r_4 * h5_p2 + e_4 * fs_20_143_105 * r_6 * h3_p2 + e_5 * fs_21_2431_770 * r_2 * h9_p2 + e_5 * fs_42_2431_715 * r_2 * h9_p4 - e_5 * fs_815_7293_6 * r_4 * h7_p2 + e_5 * fs_10_561_33 * r_4 * h7_p4 - e_5 * fs_2_39_15 * r_6 * h5_p2 - e_5 * fs_2_429_105 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_36 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_36[k] = e_0 * fs_1365_32_42 * h1_p1 - e_1 * fs_90_1_7 * h3_p1 + e_1 * fs_225_8_105 * h3_p3 - e_1 * fs_1365_16_42 * r_2 * h1_p1 + e_2 * fs_85_8_70 * h5_p1 - e_2 * fs_85_4_15 * h5_p3 + e_2 * fs_80_1_7 * r_2 * h3_p1 - e_2 * fs_25_1_105 * r_2 * h3_p3 + e_2 * fs_195_4_42 * r_4 * h1_p1 - e_3 * fs_3895_858_6 * h7_p1 + e_3 * fs_4085_286_2 * h7_p3 - e_3 * fs_255_52_70 * r_2 * h5_p1 + e_3 * fs_255_26_15 * r_2 * h5_p3 - e_3 * fs_240_11_7 * r_4 * h3_p1 + e_3 * fs_75_11_105 * r_4 * h3_p3 - e_3 * fs_65_6_42 * r_6 * h1_p1 - e_4 * fs_588_2431_210 * h9_p1 - e_4 * fs_1764_2431_55 * h9_p3 + e_4 * fs_7790_7293_6 * r_2 * h7_p1 - e_4 * fs_8170_2431_2 * r_2 * h7_p3 + e_4 * fs_17_26_70 * r_4 * h5_p1 - e_4 * fs_17_13_15 * r_4 * h5_p3 + e_4 * fs_320_143_7 * r_6 * h3_p1 - e_4 * fs_100_143_105 * r_6 * h3_p3 + e_4 * fs_65_66_42 * r_8 * h1_p1 + e_5 * fs_56_2431_210 * r_2 * h9_p1 + e_5 * fs_168_2431_55 * r_2 * h9_p3 - e_5 * fs_410_7293_6 * r_4 * h7_p1 + e_5 * fs_430_2431_2 * r_4 * h7_p3 - e_5 * fs_1_39_70 * r_6 * h5_p1 + e_5 * fs_2_39_15 * r_6 * h5_p3 - e_5 * fs_32_429_7 * r_8 * h3_p1 + e_5 * fs_10_429_105 * r_8 * h3_p3 - e_5 * fs_1_33_42 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ab_2, pc_37 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];

        pc_37[k] = e_0 * fs_1365_8_6 * h1_0 - e_1 * fs_855_8_6 * h3_0 + e_1 * fs_1125_16_10 * h3_p2 - e_1 * fs_1365_4_6 * r_2 * h1_0 + e_2 * fs_85_4_6 * h5_0 - e_2 * fs_85_8_70 * h5_p2 + e_2 * fs_95_1_6 * r_2 * h3_0 - e_2 * fs_125_2_10 * r_2 * h3_p2 + e_2 * fs_195_1_6 * r_4 * h1_0 + e_3 * fs_665_429_6 * h7_0 + e_3 * fs_2375_286_7 * h7_p2 - e_3 * fs_255_26_6 * r_2 * h5_0 + e_3 * fs_255_52_70 * r_2 * h5_p2 - e_3 * fs_285_11_6 * r_4 * h3_0 + e_3 * fs_375_22_10 * r_4 * h3_p2 - e_3 * fs_130_3_6 * r_6 * h1_0 - e_4 * fs_6174_2431_6 * h9_0 - e_4 * fs_1029_2431_165 * h9_p2 - e_4 * fs_2660_7293_6 * r_2 * h7_0 - e_4 * fs_4750_2431_7 * r_2 * h7_p2 + e_4 * fs_17_13_6 * r_4 * h5_0 - e_4 * fs_17_26_70 * r_4 * h5_p2 + e_4 * fs_380_143_6 * r_6 * h3_0 - e_4 * fs_250_143_10 * r_6 * h3_p2 + e_4 * fs_130_33_6 * r_8 * h1_0 + e_5 * fs_588_2431_6 * r_2 * h9_0 + e_5 * fs_98_2431_165 * r_2 * h9_p2 + e_5 * fs_140_7293_6 * r_4 * h7_0 + e_5 * fs_250_2431_7 * r_4 * h7_p2 - e_5 * fs_2_39_6 * r_6 * h5_0 + e_5 * fs_1_39_70 * r_6 * h5_p2 - e_5 * fs_38_429_6 * r_8 * h3_0 + e_5 * fs_25_429_10 * r_8 * h3_p2 - e_5 * fs_4_33_6 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ab_2, pc_38, pc_39 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_38[k] = - e_0 * fs_1365_16_10 * h1_m1 - e_1 * fs_45_8_15 * h3_m1 + e_1 * fs_1365_8_10 * r_2 * h1_m1 + e_2 * fs_85_4_6 * h5_m1 + e_2 * fs_5_1_15 * r_2 * h3_m1 - e_2 * fs_195_2_10 * r_4 * h1_m1 - e_3 * fs_2185_858_70 * h7_m1 - e_3 * fs_255_26_6 * r_2 * h5_m1 - e_3 * fs_15_11_15 * r_4 * h3_m1 + e_3 * fs_65_3_10 * r_6 * h1_m1 + e_4 * fs_12348_2431_2 * h9_m1 + e_4 * fs_4370_7293_70 * r_2 * h7_m1 + e_4 * fs_17_13_6 * r_4 * h5_m1 + e_4 * fs_20_143_15 * r_6 * h3_m1 - e_4 * fs_65_33_10 * r_8 * h1_m1 - e_5 * fs_1176_2431_2 * r_2 * h9_m1 - e_5 * fs_230_7293_70 * r_4 * h7_m1 - e_5 * fs_2_39_6 * r_6 * h5_m1 - e_5 * fs_2_429_15 * r_8 * h3_m1 + e_5 * fs_2_33_10 * r_10 * h1_m1;

        pc_39[k] = - e_1 * fs_1125_16_10 * h3_m2 + e_2 * fs_85_8_70 * h5_m2 + e_2 * fs_125_2_10 * r_2 * h3_m2 - e_3 * fs_2375_286_7 * h7_m2 - e_3 * fs_255_52_70 * r_2 * h5_m2 - e_3 * fs_375_22_10 * r_4 * h3_m2 + e_4 * fs_1029_2431_165 * h9_m2 + e_4 * fs_4750_2431_7 * r_2 * h7_m2 + e_4 * fs_17_26_70 * r_4 * h5_m2 + e_4 * fs_250_143_10 * r_6 * h3_m2 - e_5 * fs_98_2431_165 * r_2 * h9_m2 - e_5 * fs_250_2431_7 * r_4 * h7_m2 - e_5 * fs_1_39_70 * r_6 * h5_m2 - e_5 * fs_25_429_10 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ab_2, pc_40 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_40[k] = - e_0 * fs_1365_32_42 * h1_m1 - e_1 * fs_225_8_105 * h3_m3 + e_1 * fs_90_1_7 * h3_m1 + e_1 * fs_1365_16_42 * r_2 * h1_m1 + e_2 * fs_85_4_15 * h5_m3 - e_2 * fs_85_8_70 * h5_m1 + e_2 * fs_25_1_105 * r_2 * h3_m3 - e_2 * fs_80_1_7 * r_2 * h3_m1 - e_2 * fs_195_4_42 * r_4 * h1_m1 - e_3 * fs_4085_286_2 * h7_m3 + e_3 * fs_3895_858_6 * h7_m1 - e_3 * fs_255_26_15 * r_2 * h5_m3 + e_3 * fs_255_52_70 * r_2 * h5_m1 - e_3 * fs_75_11_105 * r_4 * h3_m3 + e_3 * fs_240_11_7 * r_4 * h3_m1 + e_3 * fs_65_6_42 * r_6 * h1_m1 + e_4 * fs_1764_2431_55 * h9_m3 + e_4 * fs_588_2431_210 * h9_m1 + e_4 * fs_8170_2431_2 * r_2 * h7_m3 - e_4 * fs_7790_7293_6 * r_2 * h7_m1 + e_4 * fs_17_13_15 * r_4 * h5_m3 - e_4 * fs_17_26_70 * r_4 * h5_m1 + e_4 * fs_100_143_105 * r_6 * h3_m3 - e_4 * fs_320_143_7 * r_6 * h3_m1 - e_4 * fs_65_66_42 * r_8 * h1_m1 - e_5 * fs_168_2431_55 * r_2 * h9_m3 - e_5 * fs_56_2431_210 * r_2 * h9_m1 - e_5 * fs_430_2431_2 * r_4 * h7_m3 + e_5 * fs_410_7293_6 * r_4 * h7_m1 - e_5 * fs_2_39_15 * r_6 * h5_m3 + e_5 * fs_1_39_70 * r_6 * h5_m1 - e_5 * fs_10_429_105 * r_8 * h3_m3 + e_5 * fs_32_429_7 * r_8 * h3_m1 + e_5 * fs_1_33_42 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph3_m2, ph5_m5, ph5_m2, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ph9_m5, ph9_m4, ph9_m3, ph9_m2, ab_2, pc_41, pc_42 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_41[k] = e_1 * fs_45_8_105 * h3_m2 - e_2 * fs_85_4_15 * h5_m2 - e_2 * fs_5_1_105 * r_2 * h3_m2 - e_3 * fs_95_66_33 * h7_m4 + e_3 * fs_15485_1716_6 * h7_m2 + e_3 * fs_255_26_15 * r_2 * h5_m2 + e_3 * fs_15_11_105 * r_4 * h3_m2 + e_4 * fs_441_2431_715 * h9_m4 + e_4 * fs_441_4862_770 * h9_m2 + e_4 * fs_190_561_33 * r_2 * h7_m4 - e_4 * fs_15485_7293_6 * r_2 * h7_m2 - e_4 * fs_17_13_15 * r_4 * h5_m2 - e_4 * fs_20_143_105 * r_6 * h3_m2 - e_5 * fs_42_2431_715 * r_2 * h9_m4 - e_5 * fs_21_2431_770 * r_2 * h9_m2 - e_5 * fs_10_561_33 * r_4 * h7_m4 + e_5 * fs_815_7293_6 * r_4 * h7_m2 + e_5 * fs_2_39_15 * r_6 * h5_m2 + e_5 * fs_2_429_105 * r_8 * h3_m2;

        pc_42[k] = - e_1 * fs_405_8_35 * h3_m3 - e_2 * f_255_2 * h5_m5 + e_2 * fs_45_1_35 * r_2 * h3_m3 + e_3 * fs_855_572_66 * h7_m5 + e_3 * fs_475_44_6 * h7_m3 + e_3 * f_765_13 * r_2 * h5_m5 - e_3 * fs_135_11_35 * r_4 * h3_m3 + e_4 * fs_294_2431_1001 * h9_m5 + e_4 * fs_294_2431_165 * h9_m3 - e_4 * fs_855_2431_66 * r_2 * h7_m5 - e_4 * fs_475_187_6 * r_2 * h7_m3 - e_4 * f_102_13 * r_4 * h5_m5 + e_4 * fs_180_143_35 * r_6 * h3_m3 - e_5 * fs_28_2431_1001 * r_2 * h9_m5 - e_5 * fs_28_2431_165 * r_2 * h9_m3 + e_5 * fs_45_2431_66 * r_4 * h7_m5 + e_5 * fs_25_187_6 * r_4 * h7_m3 + e_5 * f_4_13 * r_6 * h5_m5 - e_5 * fs_6_143_35 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m5, ph5_m4, ph5_m3, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph9_m6, ph9_m5, ph9_m4, ph9_m3, ab_2, pc_43, pc_44, pc_45, pc_46 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_43[k] = e_2 * f_255_2 * h5_m4 + e_3 * fs_285_572_4290 * h7_m6 + e_3 * fs_475_286_165 * h7_m4 - e_3 * f_765_13 * r_2 * h5_m4 + e_4 * fs_147_4862_6006 * h9_m6 + e_4 * fs_147_2431_143 * h9_m4 - e_4 * fs_285_2431_4290 * r_2 * h7_m6 - e_4 * fs_950_2431_165 * r_2 * h7_m4 + e_4 * f_102_13 * r_4 * h5_m4 - e_5 * fs_7_2431_6006 * r_2 * h9_m6 - e_5 * fs_14_2431_143 * r_2 * h9_m4 + e_5 * fs_15_2431_4290 * r_4 * h7_m6 + e_5 * fs_50_2431_165 * r_4 * h7_m4 - e_5 * f_4_13 * r_6 * h5_m4;

        pc_44[k] = e_2 * f_255_2 * h5_m5 + e_3 * fs_1425_286_66 * h7_m5 - e_3 * f_765_13 * r_2 * h5_m5 + e_4 * fs_147_2431_1001 * h9_m5 - e_4 * fs_2850_2431_66 * r_2 * h7_m5 + e_4 * f_102_13 * r_4 * h5_m5 - e_5 * fs_14_2431_1001 * r_2 * h9_m5 + e_5 * fs_150_2431_66 * r_4 * h7_m5 - e_5 * f_4_13 * r_6 * h5_m5;

        pc_45[k] = - e_2 * f_255_2 * h5_m4 + e_3 * fs_380_143_165 * h7_m4 + e_3 * f_765_13 * r_2 * h5_m4 + e_4 * fs_735_2431_143 * h9_m4 - e_4 * fs_1520_2431_165 * r_2 * h7_m4 - e_4 * f_102_13 * r_4 * h5_m4 - e_5 * fs_70_2431_143 * r_2 * h9_m4 + e_5 * fs_80_2431_165 * r_4 * h7_m4 + e_5 * f_4_13 * r_6 * h5_m4;

        pc_46[k] = e_1 * fs_675_4_7 * h3_m3 - e_2 * f_255_2 * h5_m3 - e_2 * fs_150_1_7 * r_2 * h3_m3 + e_3 * fs_2185_858_30 * h7_m3 + e_3 * f_765_13 * r_2 * h5_m3 + e_3 * fs_450_11_7 * r_4 * h3_m3 + e_4 * fs_2205_2431_33 * h9_m3 - e_4 * fs_4370_7293_30 * r_2 * h7_m3 - e_4 * f_102_13 * r_4 * h5_m3 - e_4 * fs_600_143_7 * r_6 * h3_m3 - e_5 * fs_210_2431_33 * r_2 * h9_m3 + e_5 * fs_230_7293_30 * r_4 * h7_m3 + e_5 * f_4_13 * r_6 * h5_m3 + e_5 * fs_20_143_7 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ab_2, pc_47, pc_48 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_47[k] = e_1 * fs_225_4_7 * h3_m2 - e_2 * f_85_4 * h5_m2 - e_2 * fs_50_1_7 * r_2 * h3_m2 - e_3 * fs_380_143_10 * h7_m2 + e_3 * f_255_26 * r_2 * h5_m2 + e_3 * fs_150_11_7 * r_4 * h3_m2 + e_4 * fs_735_2431_462 * h9_m2 + e_4 * fs_1520_2431_10 * r_2 * h7_m2 - e_4 * f_17_13 * r_4 * h5_m2 - e_4 * fs_200_143_7 * r_6 * h3_m2 - e_5 * fs_70_2431_462 * r_2 * h9_m2 - e_5 * fs_80_2431_10 * r_4 * h7_m2 + e_5 * f_2_39 * r_6 * h5_m2 + e_5 * fs_20_429_7 * r_8 * h3_m2;

        pc_48[k] = e_0 * fs_1365_16_15 * h1_m1 - e_1 * fs_495_8_10 * h3_m1 - e_1 * fs_1365_8_15 * r_2 * h1_m1 + e_2 * f_85_1 * h5_m1 + e_2 * fs_55_1_10 * r_2 * h3_m1 + e_2 * fs_195_2_15 * r_4 * h1_m1 - e_3 * fs_95_39_105 * h7_m1 - e_3 * f_510_13 * r_2 * h5_m1 - e_3 * fs_15_1_10 * r_4 * h3_m1 - e_3 * fs_65_3_15 * r_6 * h1_m1 + e_4 * fs_10290_2431_3 * h9_m1 + e_4 * fs_380_663_105 * r_2 * h7_m1 + e_4 * f_68_13 * r_4 * h5_m1 + e_4 * fs_20_13_10 * r_6 * h3_m1 + e_4 * fs_65_33_15 * r_8 * h1_m1 - e_5 * fs_980_2431_3 * r_2 * h9_m1 - e_5 * fs_20_663_105 * r_4 * h7_m1 - e_5 * f_8_39 * r_6 * h5_m1 - e_5 * fs_2_39_10 * r_8 * h3_m1 - e_5 * fs_2_33_15 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph7_0, ph7_p1, ph9_0, ph9_p1, ab_2, pc_49, pc_50 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];

        pc_49[k] = e_0 * f_6825_16 * h1_0 - e_1 * f_675_2 * h3_0 - e_1 * f_6825_8 * r_2 * h1_0 + e_2 * f_255_2 * h5_0 + e_2 * f_300_1 * r_2 * h3_0 + e_2 * f_975_2 * r_4 * h1_0 - e_3 * f_13300_429 * h7_0 - e_3 * f_765_13 * r_2 * h5_0 - e_3 * f_900_11 * r_4 * h3_0 - e_3 * f_325_3 * r_6 * h1_0 + e_4 * f_18522_2431 * h9_0 + e_4 * f_53200_7293 * r_2 * h7_0 + e_4 * f_102_13 * r_4 * h5_0 + e_4 * f_1200_143 * r_6 * h3_0 + e_4 * f_325_33 * r_8 * h1_0 - e_5 * f_1764_2431 * r_2 * h9_0 - e_5 * f_2800_7293 * r_4 * h7_0 - e_5 * f_4_13 * r_6 * h5_0 - e_5 * f_40_143 * r_8 * h3_0 - e_5 * f_10_33 * r_10 * h1_0;

        pc_50[k] = e_0 * fs_1365_16_15 * h1_p1 - e_1 * fs_495_8_10 * h3_p1 - e_1 * fs_1365_8_15 * r_2 * h1_p1 + e_2 * f_85_1 * h5_p1 + e_2 * fs_55_1_10 * r_2 * h3_p1 + e_2 * fs_195_2_15 * r_4 * h1_p1 - e_3 * fs_95_39_105 * h7_p1 - e_3 * f_510_13 * r_2 * h5_p1 - e_3 * fs_15_1_10 * r_4 * h3_p1 - e_3 * fs_65_3_15 * r_6 * h1_p1 + e_4 * fs_10290_2431_3 * h9_p1 + e_4 * fs_380_663_105 * r_2 * h7_p1 + e_4 * f_68_13 * r_4 * h5_p1 + e_4 * fs_20_13_10 * r_6 * h3_p1 + e_4 * fs_65_33_15 * r_8 * h1_p1 - e_5 * fs_980_2431_3 * r_2 * h9_p1 - e_5 * fs_20_663_105 * r_4 * h7_p1 - e_5 * f_8_39 * r_6 * h5_p1 - e_5 * fs_2_39_10 * r_8 * h3_p1 - e_5 * fs_2_33_15 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p4, ph7_p2, ph7_p3, ph7_p4, ph9_p2, ph9_p3, ph9_p4, ab_2, pc_51, pc_52, pc_53 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];

        pc_51[k] = e_1 * fs_225_4_7 * h3_p2 - e_2 * f_85_4 * h5_p2 - e_2 * fs_50_1_7 * r_2 * h3_p2 - e_3 * fs_380_143_10 * h7_p2 + e_3 * f_255_26 * r_2 * h5_p2 + e_3 * fs_150_11_7 * r_4 * h3_p2 + e_4 * fs_735_2431_462 * h9_p2 + e_4 * fs_1520_2431_10 * r_2 * h7_p2 - e_4 * f_17_13 * r_4 * h5_p2 - e_4 * fs_200_143_7 * r_6 * h3_p2 - e_5 * fs_70_2431_462 * r_2 * h9_p2 - e_5 * fs_80_2431_10 * r_4 * h7_p2 + e_5 * f_2_39 * r_6 * h5_p2 + e_5 * fs_20_429_7 * r_8 * h3_p2;

        pc_52[k] = e_1 * fs_675_4_7 * h3_p3 - e_2 * f_255_2 * h5_p3 - e_2 * fs_150_1_7 * r_2 * h3_p3 + e_3 * fs_2185_858_30 * h7_p3 + e_3 * f_765_13 * r_2 * h5_p3 + e_3 * fs_450_11_7 * r_4 * h3_p3 + e_4 * fs_2205_2431_33 * h9_p3 - e_4 * fs_4370_7293_30 * r_2 * h7_p3 - e_4 * f_102_13 * r_4 * h5_p3 - e_4 * fs_600_143_7 * r_6 * h3_p3 - e_5 * fs_210_2431_33 * r_2 * h9_p3 + e_5 * fs_230_7293_30 * r_4 * h7_p3 + e_5 * f_4_13 * r_6 * h5_p3 + e_5 * fs_20_143_7 * r_8 * h3_p3;

        pc_53[k] = - e_2 * f_255_2 * h5_p4 + e_3 * fs_380_143_165 * h7_p4 + e_3 * f_765_13 * r_2 * h5_p4 + e_4 * fs_735_2431_143 * h9_p4 - e_4 * fs_1520_2431_165 * r_2 * h7_p4 - e_4 * f_102_13 * r_4 * h5_p4 - e_5 * fs_70_2431_143 * r_2 * h9_p4 + e_5 * fs_80_2431_165 * r_4 * h7_p4 + e_5 * f_4_13 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph5_p5, ph7_m6, ph7_m4, ph7_p5, ph9_m6, ph9_m4, ph9_p5, ab_2, pc_54, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p5 = ph9_p5[k];

        pc_54[k] = e_2 * f_255_2 * h5_p5 + e_3 * fs_1425_286_66 * h7_p5 - e_3 * f_765_13 * r_2 * h5_p5 + e_4 * fs_147_2431_1001 * h9_p5 - e_4 * fs_2850_2431_66 * r_2 * h7_p5 + e_4 * f_102_13 * r_4 * h5_p5 - e_5 * fs_14_2431_1001 * r_2 * h9_p5 + e_5 * fs_150_2431_66 * r_4 * h7_p5 - e_5 * f_4_13 * r_6 * h5_p5;

        pc_55[k] = - e_2 * f_255_2 * h5_m4 + e_3 * fs_285_572_4290 * h7_m6 - e_3 * fs_475_286_165 * h7_m4 + e_3 * f_765_13 * r_2 * h5_m4 + e_4 * fs_147_4862_6006 * h9_m6 - e_4 * fs_147_2431_143 * h9_m4 - e_4 * fs_285_2431_4290 * r_2 * h7_m6 + e_4 * fs_950_2431_165 * r_2 * h7_m4 - e_4 * f_102_13 * r_4 * h5_m4 - e_5 * fs_7_2431_6006 * r_2 * h9_m6 + e_5 * fs_14_2431_143 * r_2 * h9_m4 + e_5 * fs_15_2431_4290 * r_4 * h7_m6 - e_5 * fs_50_2431_165 * r_4 * h7_m4 + e_5 * f_4_13 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph3_m2, ph5_m5, ph5_m2, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ph9_m5, ph9_m4, ph9_m3, ph9_m2, ab_2, pc_56, pc_57 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_56[k] = e_1 * fs_405_8_35 * h3_m3 - e_2 * f_255_2 * h5_m5 - e_2 * fs_45_1_35 * r_2 * h3_m3 + e_3 * fs_855_572_66 * h7_m5 - e_3 * fs_475_44_6 * h7_m3 + e_3 * f_765_13 * r_2 * h5_m5 + e_3 * fs_135_11_35 * r_4 * h3_m3 + e_4 * fs_294_2431_1001 * h9_m5 - e_4 * fs_294_2431_165 * h9_m3 - e_4 * fs_855_2431_66 * r_2 * h7_m5 + e_4 * fs_475_187_6 * r_2 * h7_m3 - e_4 * f_102_13 * r_4 * h5_m5 - e_4 * fs_180_143_35 * r_6 * h3_m3 - e_5 * fs_28_2431_1001 * r_2 * h9_m5 + e_5 * fs_28_2431_165 * r_2 * h9_m3 + e_5 * fs_45_2431_66 * r_4 * h7_m5 - e_5 * fs_25_187_6 * r_4 * h7_m3 + e_5 * f_4_13 * r_6 * h5_m5 + e_5 * fs_6_143_35 * r_8 * h3_m3;

        pc_57[k] = - e_1 * fs_45_8_105 * h3_m2 + e_2 * fs_85_4_15 * h5_m2 + e_2 * fs_5_1_105 * r_2 * h3_m2 - e_3 * fs_95_66_33 * h7_m4 - e_3 * fs_15485_1716_6 * h7_m2 - e_3 * fs_255_26_15 * r_2 * h5_m2 - e_3 * fs_15_11_105 * r_4 * h3_m2 + e_4 * fs_441_2431_715 * h9_m4 - e_4 * fs_441_4862_770 * h9_m2 + e_4 * fs_190_561_33 * r_2 * h7_m4 + e_4 * fs_15485_7293_6 * r_2 * h7_m2 + e_4 * fs_17_13_15 * r_4 * h5_m2 + e_4 * fs_20_143_105 * r_6 * h3_m2 - e_5 * fs_42_2431_715 * r_2 * h9_m4 + e_5 * fs_21_2431_770 * r_2 * h9_m2 - e_5 * fs_10_561_33 * r_4 * h7_m4 - e_5 * fs_815_7293_6 * r_4 * h7_m2 - e_5 * fs_2_39_15 * r_6 * h5_m2 - e_5 * fs_2_429_105 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ab_2, pc_58 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_58[k] = e_0 * fs_1365_32_42 * h1_m1 - e_1 * fs_225_8_105 * h3_m3 - e_1 * fs_90_1_7 * h3_m1 - e_1 * fs_1365_16_42 * r_2 * h1_m1 + e_2 * fs_85_4_15 * h5_m3 + e_2 * fs_85_8_70 * h5_m1 + e_2 * fs_25_1_105 * r_2 * h3_m3 + e_2 * fs_80_1_7 * r_2 * h3_m1 + e_2 * fs_195_4_42 * r_4 * h1_m1 - e_3 * fs_4085_286_2 * h7_m3 - e_3 * fs_3895_858_6 * h7_m1 - e_3 * fs_255_26_15 * r_2 * h5_m3 - e_3 * fs_255_52_70 * r_2 * h5_m1 - e_3 * fs_75_11_105 * r_4 * h3_m3 - e_3 * fs_240_11_7 * r_4 * h3_m1 - e_3 * fs_65_6_42 * r_6 * h1_m1 + e_4 * fs_1764_2431_55 * h9_m3 - e_4 * fs_588_2431_210 * h9_m1 + e_4 * fs_8170_2431_2 * r_2 * h7_m3 + e_4 * fs_7790_7293_6 * r_2 * h7_m1 + e_4 * fs_17_13_15 * r_4 * h5_m3 + e_4 * fs_17_26_70 * r_4 * h5_m1 + e_4 * fs_100_143_105 * r_6 * h3_m3 + e_4 * fs_320_143_7 * r_6 * h3_m1 + e_4 * fs_65_66_42 * r_8 * h1_m1 - e_5 * fs_168_2431_55 * r_2 * h9_m3 + e_5 * fs_56_2431_210 * r_2 * h9_m1 - e_5 * fs_430_2431_2 * r_4 * h7_m3 - e_5 * fs_410_7293_6 * r_4 * h7_m1 - e_5 * fs_2_39_15 * r_6 * h5_m3 - e_5 * fs_1_39_70 * r_6 * h5_m1 - e_5 * fs_10_429_105 * r_8 * h3_m3 - e_5 * fs_32_429_7 * r_8 * h3_m1 - e_5 * fs_1_33_42 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_m2, ph3_p1, ph5_m2, ph5_p1, ph7_m2, ph7_p1, ph9_m2, ph9_p1, ab_2, pc_59, pc_60 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p1 = ph9_p1[k];

        pc_59[k] = - e_1 * fs_1125_16_10 * h3_m2 + e_2 * fs_85_8_70 * h5_m2 + e_2 * fs_125_2_10 * r_2 * h3_m2 - e_3 * fs_2375_286_7 * h7_m2 - e_3 * fs_255_52_70 * r_2 * h5_m2 - e_3 * fs_375_22_10 * r_4 * h3_m2 + e_4 * fs_1029_2431_165 * h9_m2 + e_4 * fs_4750_2431_7 * r_2 * h7_m2 + e_4 * fs_17_26_70 * r_4 * h5_m2 + e_4 * fs_250_143_10 * r_6 * h3_m2 - e_5 * fs_98_2431_165 * r_2 * h9_m2 - e_5 * fs_250_2431_7 * r_4 * h7_m2 - e_5 * fs_1_39_70 * r_6 * h5_m2 - e_5 * fs_25_429_10 * r_8 * h3_m2;

        pc_60[k] = - e_0 * fs_1365_16_10 * h1_p1 - e_1 * fs_45_8_15 * h3_p1 + e_1 * fs_1365_8_10 * r_2 * h1_p1 + e_2 * fs_85_4_6 * h5_p1 + e_2 * fs_5_1_15 * r_2 * h3_p1 - e_2 * fs_195_2_10 * r_4 * h1_p1 - e_3 * fs_2185_858_70 * h7_p1 - e_3 * fs_255_26_6 * r_2 * h5_p1 - e_3 * fs_15_11_15 * r_4 * h3_p1 + e_3 * fs_65_3_10 * r_6 * h1_p1 + e_4 * fs_12348_2431_2 * h9_p1 + e_4 * fs_4370_7293_70 * r_2 * h7_p1 + e_4 * fs_17_13_6 * r_4 * h5_p1 + e_4 * fs_20_143_15 * r_6 * h3_p1 - e_4 * fs_65_33_10 * r_8 * h1_p1 - e_5 * fs_1176_2431_2 * r_2 * h9_p1 - e_5 * fs_230_7293_70 * r_4 * h7_p1 - e_5 * fs_2_39_6 * r_6 * h5_p1 - e_5 * fs_2_429_15 * r_8 * h3_p1 + e_5 * fs_2_33_10 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ab_2, pc_61 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];

        pc_61[k] = e_0 * fs_1365_8_6 * h1_0 - e_1 * fs_855_8_6 * h3_0 - e_1 * fs_1125_16_10 * h3_p2 - e_1 * fs_1365_4_6 * r_2 * h1_0 + e_2 * fs_85_4_6 * h5_0 + e_2 * fs_85_8_70 * h5_p2 + e_2 * fs_95_1_6 * r_2 * h3_0 + e_2 * fs_125_2_10 * r_2 * h3_p2 + e_2 * fs_195_1_6 * r_4 * h1_0 + e_3 * fs_665_429_6 * h7_0 - e_3 * fs_2375_286_7 * h7_p2 - e_3 * fs_255_26_6 * r_2 * h5_0 - e_3 * fs_255_52_70 * r_2 * h5_p2 - e_3 * fs_285_11_6 * r_4 * h3_0 - e_3 * fs_375_22_10 * r_4 * h3_p2 - e_3 * fs_130_3_6 * r_6 * h1_0 - e_4 * fs_6174_2431_6 * h9_0 + e_4 * fs_1029_2431_165 * h9_p2 - e_4 * fs_2660_7293_6 * r_2 * h7_0 + e_4 * fs_4750_2431_7 * r_2 * h7_p2 + e_4 * fs_17_13_6 * r_4 * h5_0 + e_4 * fs_17_26_70 * r_4 * h5_p2 + e_4 * fs_380_143_6 * r_6 * h3_0 + e_4 * fs_250_143_10 * r_6 * h3_p2 + e_4 * fs_130_33_6 * r_8 * h1_0 + e_5 * fs_588_2431_6 * r_2 * h9_0 - e_5 * fs_98_2431_165 * r_2 * h9_p2 + e_5 * fs_140_7293_6 * r_4 * h7_0 - e_5 * fs_250_2431_7 * r_4 * h7_p2 - e_5 * fs_2_39_6 * r_6 * h5_0 - e_5 * fs_1_39_70 * r_6 * h5_p2 - e_5 * fs_38_429_6 * r_8 * h3_0 - e_5 * fs_25_429_10 * r_8 * h3_p2 - e_5 * fs_4_33_6 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_62 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_62[k] = e_0 * fs_1365_32_42 * h1_p1 - e_1 * fs_90_1_7 * h3_p1 - e_1 * fs_225_8_105 * h3_p3 - e_1 * fs_1365_16_42 * r_2 * h1_p1 + e_2 * fs_85_8_70 * h5_p1 + e_2 * fs_85_4_15 * h5_p3 + e_2 * fs_80_1_7 * r_2 * h3_p1 + e_2 * fs_25_1_105 * r_2 * h3_p3 + e_2 * fs_195_4_42 * r_4 * h1_p1 - e_3 * fs_3895_858_6 * h7_p1 - e_3 * fs_4085_286_2 * h7_p3 - e_3 * fs_255_52_70 * r_2 * h5_p1 - e_3 * fs_255_26_15 * r_2 * h5_p3 - e_3 * fs_240_11_7 * r_4 * h3_p1 - e_3 * fs_75_11_105 * r_4 * h3_p3 - e_3 * fs_65_6_42 * r_6 * h1_p1 - e_4 * fs_588_2431_210 * h9_p1 + e_4 * fs_1764_2431_55 * h9_p3 + e_4 * fs_7790_7293_6 * r_2 * h7_p1 + e_4 * fs_8170_2431_2 * r_2 * h7_p3 + e_4 * fs_17_26_70 * r_4 * h5_p1 + e_4 * fs_17_13_15 * r_4 * h5_p3 + e_4 * fs_320_143_7 * r_6 * h3_p1 + e_4 * fs_100_143_105 * r_6 * h3_p3 + e_4 * fs_65_66_42 * r_8 * h1_p1 + e_5 * fs_56_2431_210 * r_2 * h9_p1 - e_5 * fs_168_2431_55 * r_2 * h9_p3 - e_5 * fs_410_7293_6 * r_4 * h7_p1 - e_5 * fs_430_2431_2 * r_4 * h7_p3 - e_5 * fs_1_39_70 * r_6 * h5_p1 - e_5 * fs_2_39_15 * r_6 * h5_p3 - e_5 * fs_32_429_7 * r_8 * h3_p1 - e_5 * fs_10_429_105 * r_8 * h3_p3 - e_5 * fs_1_33_42 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph3_p3, ph5_p2, ph5_p5, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ph9_p2, ph9_p3, ph9_p4, ph9_p5, ab_2, pc_63, pc_64 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];

        pc_63[k] = - e_1 * fs_45_8_105 * h3_p2 + e_2 * fs_85_4_15 * h5_p2 + e_2 * fs_5_1_105 * r_2 * h3_p2 - e_3 * fs_15485_1716_6 * h7_p2 - e_3 * fs_95_66_33 * h7_p4 - e_3 * fs_255_26_15 * r_2 * h5_p2 - e_3 * fs_15_11_105 * r_4 * h3_p2 - e_4 * fs_441_4862_770 * h9_p2 + e_4 * fs_441_2431_715 * h9_p4 + e_4 * fs_15485_7293_6 * r_2 * h7_p2 + e_4 * fs_190_561_33 * r_2 * h7_p4 + e_4 * fs_17_13_15 * r_4 * h5_p2 + e_4 * fs_20_143_105 * r_6 * h3_p2 + e_5 * fs_21_2431_770 * r_2 * h9_p2 - e_5 * fs_42_2431_715 * r_2 * h9_p4 - e_5 * fs_815_7293_6 * r_4 * h7_p2 - e_5 * fs_10_561_33 * r_4 * h7_p4 - e_5 * fs_2_39_15 * r_6 * h5_p2 - e_5 * fs_2_429_105 * r_8 * h3_p2;

        pc_64[k] = e_1 * fs_405_8_35 * h3_p3 - e_2 * f_255_2 * h5_p5 - e_2 * fs_45_1_35 * r_2 * h3_p3 - e_3 * fs_475_44_6 * h7_p3 + e_3 * fs_855_572_66 * h7_p5 + e_3 * f_765_13 * r_2 * h5_p5 + e_3 * fs_135_11_35 * r_4 * h3_p3 - e_4 * fs_294_2431_165 * h9_p3 + e_4 * fs_294_2431_1001 * h9_p5 + e_4 * fs_475_187_6 * r_2 * h7_p3 - e_4 * fs_855_2431_66 * r_2 * h7_p5 - e_4 * f_102_13 * r_4 * h5_p5 - e_4 * fs_180_143_35 * r_6 * h3_p3 + e_5 * fs_28_2431_165 * r_2 * h9_p3 - e_5 * fs_28_2431_1001 * r_2 * h9_p5 - e_5 * fs_25_187_6 * r_4 * h7_p3 + e_5 * fs_45_2431_66 * r_4 * h7_p5 + e_5 * f_4_13 * r_6 * h5_p5 + e_5 * fs_6_143_35 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m3, ph5_p4, ph7_m7, ph7_m3, ph7_p4, ph7_p6, ph9_m7, ph9_m3, ph9_p4, ph9_p6, ab_2, pc_65, pc_66 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];

        pc_65[k] = - e_2 * f_255_2 * h5_p4 - e_3 * fs_475_286_165 * h7_p4 + e_3 * fs_285_572_4290 * h7_p6 + e_3 * f_765_13 * r_2 * h5_p4 - e_4 * fs_147_2431_143 * h9_p4 + e_4 * fs_147_4862_6006 * h9_p6 + e_4 * fs_950_2431_165 * r_2 * h7_p4 - e_4 * fs_285_2431_4290 * r_2 * h7_p6 - e_4 * f_102_13 * r_4 * h5_p4 + e_5 * fs_14_2431_143 * r_2 * h9_p4 - e_5 * fs_7_2431_6006 * r_2 * h9_p6 - e_5 * fs_50_2431_165 * r_4 * h7_p4 + e_5 * fs_15_2431_4290 * r_4 * h7_p6 + e_5 * f_4_13 * r_6 * h5_p4;

        pc_66[k] = e_1 * fs_675_8_7 * h3_m3 + e_2 * f_255_2 * h5_m3 - e_2 * fs_75_1_7 * r_2 * h3_m3 + e_3 * fs_95_572_30030 * h7_m7 + e_3 * fs_1425_572_30 * h7_m3 - e_3 * f_765_13 * r_2 * h5_m3 + e_3 * fs_225_11_7 * r_4 * h3_m3 + e_4 * fs_294_2431_1001 * h9_m7 + e_4 * fs_147_2431_33 * h9_m3 - e_4 * fs_95_2431_30030 * r_2 * h7_m7 - e_4 * fs_1425_2431_30 * r_2 * h7_m3 + e_4 * f_102_13 * r_4 * h5_m3 - e_4 * fs_300_143_7 * r_6 * h3_m3 - e_5 * fs_28_2431_1001 * r_2 * h9_m7 - e_5 * fs_14_2431_33 * r_2 * h9_m3 + e_5 * fs_5_2431_30030 * r_4 * h7_m7 + e_5 * fs_75_2431_30 * r_4 * h7_m3 - e_5 * f_4_13 * r_6 * h5_m3 + e_5 * fs_10_143_7 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ab_2, pc_67 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];

        pc_67[k] = - e_1 * fs_45_2_105 * h3_m2 + e_2 * fs_85_4_15 * h5_m2 + e_2 * fs_20_1_105 * r_2 * h3_m2 - e_3 * fs_95_286_858 * h7_m6 + e_3 * fs_2565_286_6 * h7_m2 - e_3 * fs_255_26_15 * r_2 * h5_m2 - e_3 * fs_60_11_105 * r_4 * h3_m2 + e_4 * fs_147_4862_30030 * h9_m6 + e_4 * fs_147_4862_770 * h9_m2 + e_4 * fs_190_2431_858 * r_2 * h7_m6 - e_4 * fs_5130_2431_6 * r_2 * h7_m2 + e_4 * fs_17_13_15 * r_4 * h5_m2 + e_4 * fs_80_143_105 * r_6 * h3_m2 - e_5 * fs_7_2431_30030 * r_2 * h9_m6 - e_5 * fs_7_2431_770 * r_2 * h9_m2 - e_5 * fs_10_2431_858 * r_4 * h7_m6 + e_5 * fs_270_2431_6 * r_4 * h7_m2 - e_5 * fs_2_39_15 * r_6 * h5_m2 - e_5 * fs_8_429_105 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m4, ph7_m5, ph7_m4, ph7_m1, ph9_m5, ph9_m4, ph9_m1, ab_2, pc_68, pc_69 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m1 = ph9_m1[k];

        pc_68[k] = e_0 * fs_1365_16_14 * h1_m1 - e_1 * fs_405_8_21 * h3_m1 - e_1 * fs_1365_8_14 * r_2 * h1_m1 + e_2 * f_255_2 * h5_m5 + e_2 * fs_45_1_21 * r_2 * h3_m1 + e_2 * fs_195_2_14 * r_4 * h1_m1 - e_3 * fs_5035_1716_66 * h7_m5 + e_3 * fs_2375_132_2 * h7_m1 - e_3 * f_765_13 * r_2 * h5_m5 - e_3 * fs_135_11_21 * r_4 * h3_m1 - e_3 * fs_65_3_14 * r_6 * h1_m1 + e_4 * fs_441_2431_1001 * h9_m5 + e_4 * fs_441_2431_70 * h9_m1 + e_4 * fs_5035_7293_66 * r_2 * h7_m5 - e_4 * fs_2375_561_2 * r_2 * h7_m1 + e_4 * f_102_13 * r_4 * h5_m5 + e_4 * fs_180_143_21 * r_6 * h3_m1 + e_4 * fs_65_33_14 * r_8 * h1_m1 - e_5 * fs_42_2431_1001 * r_2 * h9_m5 - e_5 * fs_42_2431_70 * r_2 * h9_m1 - e_5 * fs_265_7293_66 * r_4 * h7_m5 + e_5 * fs_125_561_2 * r_4 * h7_m1 - e_5 * f_4_13 * r_6 * h5_m5 - e_5 * fs_6_143_21 * r_8 * h3_m1 - e_5 * fs_2_33_14 * r_10 * h1_m1;

        pc_69[k] = e_2 * fs_85_4_15 * h5_m4 - e_3 * fs_950_143_11 * h7_m4 - e_3 * fs_255_26_15 * r_2 * h5_m4 + e_4 * fs_294_2431_2145 * h9_m4 + e_4 * fs_3800_2431_11 * r_2 * h7_m4 + e_4 * fs_17_13_15 * r_4 * h5_m4 - e_5 * fs_28_2431_2145 * r_2 * h9_m4 - e_5 * fs_200_2431_11 * r_4 * h7_m4 - e_5 * fs_2_39_15 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ab_2, pc_70 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_70[k] = e_0 * fs_1365_16_3 * h1_m1 + e_1 * fs_675_16_30 * h3_m3 + e_1 * fs_1935_16_2 * h3_m1 - e_1 * fs_1365_8_3 * r_2 * h1_m1 - e_2 * fs_85_2_5 * h5_m1 - e_2 * fs_75_2_30 * r_2 * h3_m3 - e_2 * fs_215_2_2 * r_2 * h3_m1 + e_2 * fs_195_2_3 * r_4 * h1_m1 - e_3 * fs_95_22_7 * h7_m3 + e_3 * fs_2755_858_21 * h7_m1 + e_3 * fs_255_13_5 * r_2 * h5_m1 + e_3 * fs_225_22_30 * r_4 * h3_m3 + e_3 * fs_645_22_2 * r_4 * h3_m1 - e_3 * fs_65_3_3 * r_6 * h1_m1 + e_4 * fs_441_2431_770 * h9_m3 + e_4 * fs_2058_2431_15 * h9_m1 + e_4 * fs_190_187_7 * r_2 * h7_m3 - e_4 * fs_5510_7293_21 * r_2 * h7_m1 - e_4 * fs_34_13_5 * r_4 * h5_m1 - e_4 * fs_150_143_30 * r_6 * h3_m3 - e_4 * fs_430_143_2 * r_6 * h3_m1 + e_4 * fs_65_33_3 * r_8 * h1_m1 - e_5 * fs_42_2431_770 * r_2 * h9_m3 - e_5 * fs_196_2431_15 * r_2 * h9_m1 - e_5 * fs_10_187_7 * r_4 * h7_m3 + e_5 * fs_290_7293_21 * r_4 * h7_m1 + e_5 * fs_4_39_5 * r_6 * h5_m1 + e_5 * fs_5_143_30 * r_8 * h3_m3 + e_5 * fs_43_429_2 * r_8 * h3_m1 - e_5 * fs_2_33_3 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ab_2, pc_71 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];

        pc_71[k] = e_1 * fs_225_1_3 * h3_p2 - e_2 * fs_85_4_21 * h5_p2 - e_2 * fs_200_1_3 * r_2 * h3_p2 + e_3 * fs_95_429_210 * h7_p2 + e_3 * fs_255_26_21 * r_2 * h5_p2 + e_3 * fs_600_11_3 * r_4 * h3_p2 + e_4 * fs_3087_2431_22 * h9_p2 - e_4 * fs_380_7293_210 * r_2 * h7_p2 - e_4 * fs_17_13_21 * r_4 * h5_p2 - e_4 * fs_800_143_3 * r_6 * h3_p2 - e_5 * fs_294_2431_22 * r_2 * h9_p2 + e_5 * fs_20_7293_210 * r_4 * h7_p2 + e_5 * fs_2_39_21 * r_6 * h5_p2 + e_5 * fs_80_429_3 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_72 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_72[k] = - e_0 * fs_1365_16_3 * h1_p1 - e_1 * fs_1935_16_2 * h3_p1 + e_1 * fs_675_16_30 * h3_p3 + e_1 * fs_1365_8_3 * r_2 * h1_p1 + e_2 * fs_85_2_5 * h5_p1 + e_2 * fs_215_2_2 * r_2 * h3_p1 - e_2 * fs_75_2_30 * r_2 * h3_p3 - e_2 * fs_195_2_3 * r_4 * h1_p1 - e_3 * fs_2755_858_21 * h7_p1 - e_3 * fs_95_22_7 * h7_p3 - e_3 * fs_255_13_5 * r_2 * h5_p1 - e_3 * fs_645_22_2 * r_4 * h3_p1 + e_3 * fs_225_22_30 * r_4 * h3_p3 + e_3 * fs_65_3_3 * r_6 * h1_p1 - e_4 * fs_2058_2431_15 * h9_p1 + e_4 * fs_441_2431_770 * h9_p3 + e_4 * fs_5510_7293_21 * r_2 * h7_p1 + e_4 * fs_190_187_7 * r_2 * h7_p3 + e_4 * fs_34_13_5 * r_4 * h5_p1 + e_4 * fs_430_143_2 * r_6 * h3_p1 - e_4 * fs_150_143_30 * r_6 * h3_p3 - e_4 * fs_65_33_3 * r_8 * h1_p1 + e_5 * fs_196_2431_15 * r_2 * h9_p1 - e_5 * fs_42_2431_770 * r_2 * h9_p3 - e_5 * fs_290_7293_21 * r_4 * h7_p1 - e_5 * fs_10_187_7 * r_4 * h7_p3 - e_5 * fs_4_39_5 * r_6 * h5_p1 - e_5 * fs_43_429_2 * r_8 * h3_p1 + e_5 * fs_5_143_30 * r_8 * h3_p3 + e_5 * fs_2_33_3 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2, pc_73 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_73[k] = e_0 * fs_1365_16_21 * h1_0 - e_1 * fs_45_4_21 * h3_0 - e_1 * fs_1365_8_21 * r_2 * h1_0 - e_2 * fs_85_4_21 * h5_0 + e_2 * fs_85_4_15 * h5_p4 + e_2 * fs_10_1_21 * r_2 * h3_0 + e_2 * fs_195_2_21 * r_4 * h1_0 + e_3 * fs_3040_429_21 * h7_0 - e_3 * fs_950_143_11 * h7_p4 + e_3 * fs_255_26_21 * r_2 * h5_0 - e_3 * fs_255_26_15 * r_2 * h5_p4 - e_3 * fs_30_11_21 * r_4 * h3_0 - e_3 * fs_65_3_21 * r_6 * h1_0 + e_4 * fs_1764_2431_21 * h9_0 + e_4 * fs_294_2431_2145 * h9_p4 - e_4 * fs_12160_7293_21 * r_2 * h7_0 + e_4 * fs_3800_2431_11 * r_2 * h7_p4 - e_4 * fs_17_13_21 * r_4 * h5_0 + e_4 * fs_17_13_15 * r_4 * h5_p4 + e_4 * fs_40_143_21 * r_6 * h3_0 + e_4 * fs_65_33_21 * r_8 * h1_0 - e_5 * fs_168_2431_21 * r_2 * h9_0 - e_5 * fs_28_2431_2145 * r_2 * h9_p4 + e_5 * fs_640_7293_21 * r_4 * h7_0 - e_5 * fs_200_2431_11 * r_4 * h7_p4 + e_5 * fs_2_39_21 * r_6 * h5_0 - e_5 * fs_2_39_15 * r_6 * h5_p4 - e_5 * fs_4_429_21 * r_8 * h3_0 - e_5 * fs_2_33_21 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_74 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_74[k] = e_0 * fs_1365_16_14 * h1_p1 - e_1 * fs_405_8_21 * h3_p1 - e_1 * fs_1365_8_14 * r_2 * h1_p1 + e_2 * f_255_2 * h5_p5 + e_2 * fs_45_1_21 * r_2 * h3_p1 + e_2 * fs_195_2_14 * r_4 * h1_p1 + e_3 * fs_2375_132_2 * h7_p1 - e_3 * fs_5035_1716_66 * h7_p5 - e_3 * f_765_13 * r_2 * h5_p5 - e_3 * fs_135_11_21 * r_4 * h3_p1 - e_3 * fs_65_3_14 * r_6 * h1_p1 + e_4 * fs_441_2431_70 * h9_p1 + e_4 * fs_441_2431_1001 * h9_p5 - e_4 * fs_2375_561_2 * r_2 * h7_p1 + e_4 * fs_5035_7293_66 * r_2 * h7_p5 + e_4 * f_102_13 * r_4 * h5_p5 + e_4 * fs_180_143_21 * r_6 * h3_p1 + e_4 * fs_65_33_14 * r_8 * h1_p1 - e_5 * fs_42_2431_70 * r_2 * h9_p1 - e_5 * fs_42_2431_1001 * r_2 * h9_p5 + e_5 * fs_125_561_2 * r_4 * h7_p1 - e_5 * fs_265_7293_66 * r_4 * h7_p5 - e_5 * f_4_13 * r_6 * h5_p5 - e_5 * fs_6_143_21 * r_8 * h3_p1 - e_5 * fs_2_33_14 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph7_p6, ph7_p7, ph9_p2, ph9_p3, ph9_p6, ph9_p7, ab_2, pc_75, pc_76 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_75[k] = - e_1 * fs_45_2_105 * h3_p2 + e_2 * fs_85_4_15 * h5_p2 + e_2 * fs_20_1_105 * r_2 * h3_p2 + e_3 * fs_2565_286_6 * h7_p2 - e_3 * fs_95_286_858 * h7_p6 - e_3 * fs_255_26_15 * r_2 * h5_p2 - e_3 * fs_60_11_105 * r_4 * h3_p2 + e_4 * fs_147_4862_770 * h9_p2 + e_4 * fs_147_4862_30030 * h9_p6 - e_4 * fs_5130_2431_6 * r_2 * h7_p2 + e_4 * fs_190_2431_858 * r_2 * h7_p6 + e_4 * fs_17_13_15 * r_4 * h5_p2 + e_4 * fs_80_143_105 * r_6 * h3_p2 - e_5 * fs_7_2431_770 * r_2 * h9_p2 - e_5 * fs_7_2431_30030 * r_2 * h9_p6 + e_5 * fs_270_2431_6 * r_4 * h7_p2 - e_5 * fs_10_2431_858 * r_4 * h7_p6 - e_5 * fs_2_39_15 * r_6 * h5_p2 - e_5 * fs_8_429_105 * r_8 * h3_p2;

        pc_76[k] = e_1 * fs_675_8_7 * h3_p3 + e_2 * f_255_2 * h5_p3 - e_2 * fs_75_1_7 * r_2 * h3_p3 + e_3 * fs_1425_572_30 * h7_p3 + e_3 * fs_95_572_30030 * h7_p7 - e_3 * f_765_13 * r_2 * h5_p3 + e_3 * fs_225_11_7 * r_4 * h3_p3 + e_4 * fs_147_2431_33 * h9_p3 + e_4 * fs_294_2431_1001 * h9_p7 - e_4 * fs_1425_2431_30 * r_2 * h7_p3 - e_4 * fs_95_2431_30030 * r_2 * h7_p7 + e_4 * f_102_13 * r_4 * h5_p3 - e_4 * fs_300_143_7 * r_6 * h3_p3 - e_5 * fs_14_2431_33 * r_2 * h9_p3 - e_5 * fs_28_2431_1001 * r_2 * h9_p7 + e_5 * fs_75_2431_30 * r_4 * h7_p3 + e_5 * fs_5_2431_30030 * r_4 * h7_p7 - e_5 * f_4_13 * r_6 * h5_p3 + e_5 * fs_10_143_7 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ab_2, pc_77 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];

        pc_77[k] = - e_1 * fs_1575_8_3 * h3_m2 - e_2 * fs_85_4_21 * h5_m2 + e_2 * fs_175_1_3 * r_2 * h3_m2 - e_3 * fs_285_572_210 * h7_m2 + e_3 * fs_255_26_21 * r_2 * h5_m2 - e_3 * fs_525_11_3 * r_4 * h3_m2 + e_4 * fs_294_2431_2431 * h9_m8 - e_4 * fs_147_4862_22 * h9_m2 + e_4 * fs_285_2431_210 * r_2 * h7_m2 - e_4 * fs_17_13_21 * r_4 * h5_m2 + e_4 * fs_700_143_3 * r_6 * h3_m2 - e_5 * fs_28_2431_2431 * r_2 * h9_m8 + e_5 * fs_7_2431_22 * r_2 * h9_m2 - e_5 * fs_15_2431_210 * r_4 * h7_m2 + e_5 * fs_2_39_21 * r_6 * h5_m2 - e_5 * fs_70_429_3 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m6, ph7_m1, ph9_m7, ph9_m6, ph9_m1, ab_2, pc_78, pc_79 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m1 = ph9_m1[k];

        pc_78[k] = e_0 * fs_4095_16_2 * h1_m1 - e_1 * fs_315_8_3 * h3_m1 - e_1 * fs_4095_8_2 * r_2 * h1_m1 - e_2 * fs_85_4_30 * h5_m1 + e_2 * fs_35_1_3 * r_2 * h3_m1 + e_2 * fs_585_2_2 * r_4 * h1_m1 - e_3 * fs_665_572_858 * h7_m7 - e_3 * fs_2185_572_14 * h7_m1 + e_3 * fs_255_26_30 * r_2 * h5_m1 - e_3 * fs_105_11_3 * r_4 * h3_m1 - e_3 * fs_65_1_2 * r_6 * h1_m1 + e_4 * fs_588_2431_715 * h9_m7 - e_4 * fs_294_2431_10 * h9_m1 + e_4 * fs_665_2431_858 * r_2 * h7_m7 + e_4 * fs_2185_2431_14 * r_2 * h7_m1 - e_4 * fs_17_13_30 * r_4 * h5_m1 + e_4 * fs_140_143_3 * r_6 * h3_m1 + e_4 * fs_65_11_2 * r_8 * h1_m1 - e_5 * fs_56_2431_715 * r_2 * h9_m7 + e_5 * fs_28_2431_10 * r_2 * h9_m1 - e_5 * fs_35_2431_858 * r_4 * h7_m7 - e_5 * fs_115_2431_14 * r_4 * h7_m1 + e_5 * fs_2_39_30 * r_6 * h5_m1 - e_5 * fs_14_429_3 * r_8 * h3_m1 - e_5 * fs_2_11_2 * r_10 * h1_m1;

        pc_79[k] = - e_3 * fs_475_1716_6006 * h7_m6 + e_4 * fs_441_4862_4290 * h9_m6 + e_4 * fs_475_7293_6006 * r_2 * h7_m6 - e_5 * fs_21_2431_4290 * r_2 * h9_m6 - e_5 * fs_25_7293_6006 * r_4 * h7_m6;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ab_2, pc_80 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_80[k] = e_0 * fs_1365_32_6 * h1_m1 + e_1 * f_2205_8 * h3_m1 - e_1 * fs_1365_16_6 * r_2 * h1_m1 - e_2 * fs_85_4_21 * h5_m5 - e_2 * fs_85_8_10 * h5_m1 - e_2 * f_245_1 * r_2 * h3_m1 + e_2 * fs_195_4_6 * r_4 * h1_m1 - e_3 * fs_95_286_154 * h7_m5 - e_3 * fs_3325_858_42 * h7_m1 + e_3 * fs_255_26_21 * r_2 * h5_m5 + e_3 * fs_255_52_10 * r_2 * h5_m1 + e_3 * f_735_11 * r_4 * h3_m1 - e_3 * fs_65_6_6 * r_6 * h1_m1 + e_4 * fs_588_2431_429 * h9_m5 - e_4 * fs_588_2431_30 * h9_m1 + e_4 * fs_190_2431_154 * r_2 * h7_m5 + e_4 * fs_6650_7293_42 * r_2 * h7_m1 - e_4 * fs_17_13_21 * r_4 * h5_m5 - e_4 * fs_17_26_10 * r_4 * h5_m1 - e_4 * f_980_143 * r_6 * h3_m1 + e_4 * fs_65_66_6 * r_8 * h1_m1 - e_5 * fs_56_2431_429 * r_2 * h9_m5 + e_5 * fs_56_2431_30 * r_2 * h9_m1 - e_5 * fs_10_2431_154 * r_4 * h7_m5 - e_5 * fs_350_7293_42 * r_4 * h7_m1 + e_5 * fs_2_39_21 * r_6 * h5_m5 + e_5 * fs_1_39_10 * r_6 * h5_m1 + e_5 * f_98_429 * r_8 * h3_m1 - e_5 * fs_1_33_6 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph3_p3, ph5_m4, ph5_m2, ph5_p3, ph7_m4, ph7_m2, ph7_p3, ph9_m4, ph9_m2, ph9_p3, ab_2, pc_81, pc_82 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p3 = ph9_p3[k];

        pc_81[k] = - e_1 * fs_495_16_70 * h3_m2 - e_2 * fs_85_4_30 * h5_m4 - e_2 * fs_85_8_10 * h5_m2 + e_2 * fs_55_2_70 * r_2 * h3_m2 + e_3 * fs_665_286_22 * h7_m4 + e_3 * f_665_26 * h7_m2 + e_3 * fs_255_26_30 * r_2 * h5_m4 + e_3 * fs_255_52_10 * r_2 * h5_m2 - e_3 * fs_15_2_70 * r_4 * h3_m2 + e_4 * fs_147_2431_4290 * h9_m4 + e_4 * fs_147_2431_1155 * h9_m2 - e_4 * fs_1330_2431_22 * r_2 * h7_m4 - e_4 * f_1330_221 * r_2 * h7_m2 - e_4 * fs_17_13_30 * r_4 * h5_m4 - e_4 * fs_17_26_10 * r_4 * h5_m2 + e_4 * fs_10_13_70 * r_6 * h3_m2 - e_5 * fs_14_2431_4290 * r_2 * h9_m4 - e_5 * fs_14_2431_1155 * r_2 * h9_m2 + e_5 * fs_70_2431_22 * r_4 * h7_m4 + e_5 * f_70_221 * r_4 * h7_m2 + e_5 * fs_2_39_30 * r_6 * h5_m4 + e_5 * fs_1_39_10 * r_6 * h5_m2 - e_5 * fs_1_39_70 * r_8 * h3_m2;

        pc_82[k] = - e_1 * fs_675_8_7 * h3_p3 - e_2 * f_255_2 * h5_p3 + e_2 * fs_75_1_7 * r_2 * h3_p3 + e_3 * fs_4655_858_30 * h7_p3 + e_3 * f_765_13 * r_2 * h5_p3 - e_3 * fs_225_11_7 * r_4 * h3_p3 + e_4 * fs_1764_2431_33 * h9_p3 - e_4 * fs_9310_7293_30 * r_2 * h7_p3 - e_4 * f_102_13 * r_4 * h5_p3 + e_4 * fs_300_143_7 * r_6 * h3_p3 - e_5 * fs_168_2431_33 * r_2 * h9_p3 + e_5 * fs_490_7293_30 * r_4 * h7_p3 + e_5 * f_4_13 * r_6 * h5_p3 - e_5 * fs_10_143_7 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_83 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_83[k] = e_1 * fs_495_16_70 * h3_p2 + e_2 * fs_85_8_10 * h5_p2 - e_2 * fs_85_4_30 * h5_p4 - e_2 * fs_55_2_70 * r_2 * h3_p2 - e_3 * f_665_26 * h7_p2 + e_3 * fs_665_286_22 * h7_p4 - e_3 * fs_255_52_10 * r_2 * h5_p2 + e_3 * fs_255_26_30 * r_2 * h5_p4 + e_3 * fs_15_2_70 * r_4 * h3_p2 - e_4 * fs_147_2431_1155 * h9_p2 + e_4 * fs_147_2431_4290 * h9_p4 + e_4 * f_1330_221 * r_2 * h7_p2 - e_4 * fs_1330_2431_22 * r_2 * h7_p4 + e_4 * fs_17_26_10 * r_4 * h5_p2 - e_4 * fs_17_13_30 * r_4 * h5_p4 - e_4 * fs_10_13_70 * r_6 * h3_p2 + e_5 * fs_14_2431_1155 * r_2 * h9_p2 - e_5 * fs_14_2431_4290 * r_2 * h9_p4 - e_5 * f_70_221 * r_4 * h7_p2 + e_5 * fs_70_2431_22 * r_4 * h7_p4 - e_5 * fs_1_39_10 * r_6 * h5_p2 + e_5 * fs_2_39_30 * r_6 * h5_p4 + e_5 * fs_1_39_70 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_84 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_84[k] = - e_0 * fs_1365_32_6 * h1_p1 - e_1 * f_2205_8 * h3_p1 + e_1 * fs_1365_16_6 * r_2 * h1_p1 + e_2 * fs_85_8_10 * h5_p1 - e_2 * fs_85_4_21 * h5_p5 + e_2 * f_245_1 * r_2 * h3_p1 - e_2 * fs_195_4_6 * r_4 * h1_p1 + e_3 * fs_3325_858_42 * h7_p1 - e_3 * fs_95_286_154 * h7_p5 - e_3 * fs_255_52_10 * r_2 * h5_p1 + e_3 * fs_255_26_21 * r_2 * h5_p5 - e_3 * f_735_11 * r_4 * h3_p1 + e_3 * fs_65_6_6 * r_6 * h1_p1 + e_4 * fs_588_2431_30 * h9_p1 + e_4 * fs_588_2431_429 * h9_p5 - e_4 * fs_6650_7293_42 * r_2 * h7_p1 + e_4 * fs_190_2431_154 * r_2 * h7_p5 + e_4 * fs_17_26_10 * r_4 * h5_p1 - e_4 * fs_17_13_21 * r_4 * h5_p5 + e_4 * f_980_143 * r_6 * h3_p1 - e_4 * fs_65_66_6 * r_8 * h1_p1 - e_5 * fs_56_2431_30 * r_2 * h9_p1 - e_5 * fs_56_2431_429 * r_2 * h9_p5 + e_5 * fs_350_7293_42 * r_4 * h7_p1 - e_5 * fs_10_2431_154 * r_4 * h7_p5 - e_5 * fs_1_39_10 * r_6 * h5_p1 + e_5 * fs_2_39_21 * r_6 * h5_p5 - e_5 * f_98_429 * r_8 * h3_p1 + e_5 * fs_1_33_6 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ab_2, pc_85 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];

        pc_85[k] = e_0 * f_1365_4 * h1_0 + e_1 * f_945_4 * h3_0 - e_1 * f_1365_2 * r_2 * h1_0 - e_2 * f_255_2 * h5_0 - e_2 * f_210_1 * r_2 * h3_0 + e_2 * f_390_1 * r_4 * h1_0 - e_3 * f_12635_429 * h7_0 - e_3 * fs_475_1716_6006 * h7_p6 + e_3 * f_765_13 * r_2 * h5_0 + e_3 * f_630_11 * r_4 * h3_0 - e_3 * f_260_3 * r_6 * h1_0 - e_4 * f_2646_2431 * h9_0 + e_4 * fs_441_4862_4290 * h9_p6 + e_4 * f_50540_7293 * r_2 * h7_0 + e_4 * fs_475_7293_6006 * r_2 * h7_p6 - e_4 * f_102_13 * r_4 * h5_0 - e_4 * f_840_143 * r_6 * h3_0 + e_4 * f_260_33 * r_8 * h1_0 + e_5 * f_252_2431 * r_2 * h9_0 - e_5 * fs_21_2431_4290 * r_2 * h9_p6 - e_5 * f_2660_7293 * r_4 * h7_0 - e_5 * fs_25_7293_6006 * r_4 * h7_p6 + e_5 * f_4_13 * r_6 * h5_0 + e_5 * f_28_143 * r_8 * h3_0 - e_5 * f_8_33 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2, pc_86 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_86[k] = e_0 * fs_4095_16_2 * h1_p1 - e_1 * fs_315_8_3 * h3_p1 - e_1 * fs_4095_8_2 * r_2 * h1_p1 - e_2 * fs_85_4_30 * h5_p1 + e_2 * fs_35_1_3 * r_2 * h3_p1 + e_2 * fs_585_2_2 * r_4 * h1_p1 - e_3 * fs_2185_572_14 * h7_p1 - e_3 * fs_665_572_858 * h7_p7 + e_3 * fs_255_26_30 * r_2 * h5_p1 - e_3 * fs_105_11_3 * r_4 * h3_p1 - e_3 * fs_65_1_2 * r_6 * h1_p1 - e_4 * fs_294_2431_10 * h9_p1 + e_4 * fs_588_2431_715 * h9_p7 + e_4 * fs_2185_2431_14 * r_2 * h7_p1 + e_4 * fs_665_2431_858 * r_2 * h7_p7 - e_4 * fs_17_13_30 * r_4 * h5_p1 + e_4 * fs_140_143_3 * r_6 * h3_p1 + e_4 * fs_65_11_2 * r_8 * h1_p1 + e_5 * fs_28_2431_10 * r_2 * h9_p1 - e_5 * fs_56_2431_715 * r_2 * h9_p7 - e_5 * fs_115_2431_14 * r_4 * h7_p1 - e_5 * fs_35_2431_858 * r_4 * h7_p7 + e_5 * fs_2_39_30 * r_6 * h5_p1 - e_5 * fs_14_429_3 * r_8 * h3_p1 - e_5 * fs_2_11_2 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph3_p2, ph5_m1, ph5_p2, ph7_m1, ph7_p2, ph9_m9, ph9_m1, ph9_p2, ph9_p8, ab_2, pc_87, pc_88 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];

        pc_87[k] = - e_1 * fs_1575_8_3 * h3_p2 - e_2 * fs_85_4_21 * h5_p2 + e_2 * fs_175_1_3 * r_2 * h3_p2 - e_3 * fs_285_572_210 * h7_p2 + e_3 * fs_255_26_21 * r_2 * h5_p2 - e_3 * fs_525_11_3 * r_4 * h3_p2 - e_4 * fs_147_4862_22 * h9_p2 + e_4 * fs_294_2431_2431 * h9_p8 + e_4 * fs_285_2431_210 * r_2 * h7_p2 - e_4 * fs_17_13_21 * r_4 * h5_p2 + e_4 * fs_700_143_3 * r_6 * h3_p2 + e_5 * fs_7_2431_22 * r_2 * h9_p2 - e_5 * fs_28_2431_2431 * r_2 * h9_p8 - e_5 * fs_15_2431_210 * r_4 * h7_p2 + e_5 * fs_2_39_21 * r_6 * h5_p2 - e_5 * fs_70_429_3 * r_8 * h3_p2;

        pc_88[k] = e_0 * fs_4095_32_10 * h1_m1 + e_1 * fs_315_4_15 * h3_m1 - e_1 * fs_4095_16_10 * r_2 * h1_m1 + e_2 * fs_85_4_6 * h5_m1 - e_2 * fs_70_1_15 * r_2 * h3_m1 + e_2 * fs_585_4_10 * r_4 * h1_m1 + e_3 * fs_95_286_70 * h7_m1 - e_3 * fs_255_26_6 * r_2 * h5_m1 + e_3 * fs_210_11_15 * r_4 * h3_m1 - e_3 * fs_65_2_10 * r_6 * h1_m1 + e_4 * fs_441_2431_2431 * h9_m9 + e_4 * fs_147_4862_2 * h9_m1 - e_4 * fs_190_2431_70 * r_2 * h7_m1 + e_4 * fs_17_13_6 * r_4 * h5_m1 - e_4 * fs_280_143_15 * r_6 * h3_m1 + e_4 * fs_65_22_10 * r_8 * h1_m1 - e_5 * fs_42_2431_2431 * r_2 * h9_m9 - e_5 * fs_7_2431_2 * r_2 * h9_m1 + e_5 * fs_10_2431_70 * r_4 * h7_m1 - e_5 * fs_2_39_6 * r_6 * h5_m1 + e_5 * fs_28_429_15 * r_8 * h3_m1 - e_5 * fs_1_11_10 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m8, ph9_m7, ph9_m1, ab_2, pc_89, pc_90 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];

        pc_89[k] = e_4 * fs_147_2431_12155 * h9_m8 - e_5 * fs_14_2431_12155 * r_2 * h9_m8;

        pc_90[k] = e_0 * fs_1365_32_2 * h1_m1 + e_1 * fs_315_2_3 * h3_m1 - e_1 * fs_1365_16_2 * r_2 * h1_m1 + e_2 * fs_85_4_30 * h5_m1 - e_2 * fs_140_1_3 * r_2 * h3_m1 + e_2 * fs_195_4_2 * r_4 * h1_m1 + e_3 * fs_665_858_858 * h7_m7 + e_3 * fs_1330_429_14 * h7_m1 - e_3 * fs_255_26_30 * r_2 * h5_m1 + e_3 * fs_420_11_3 * r_4 * h3_m1 - e_3 * fs_65_6_2 * r_6 * h1_m1 + e_4 * fs_441_2431_715 * h9_m7 + e_4 * fs_441_4862_10 * h9_m1 - e_4 * fs_1330_7293_858 * r_2 * h7_m7 - e_4 * fs_5320_7293_14 * r_2 * h7_m1 + e_4 * fs_17_13_30 * r_4 * h5_m1 - e_4 * fs_560_143_3 * r_6 * h3_m1 + e_4 * fs_65_66_2 * r_8 * h1_m1 - e_5 * fs_42_2431_715 * r_2 * h9_m7 - e_5 * fs_21_2431_10 * r_2 * h9_m1 + e_5 * fs_70_7293_858 * r_4 * h7_m7 + e_5 * fs_280_7293_14 * r_4 * h7_m1 - e_5 * fs_2_39_30 * r_6 * h5_m1 + e_5 * fs_56_429_3 * r_8 * h3_m1 - e_5 * fs_1_33_2 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ab_2, pc_91 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];

        pc_91[k] = - e_1 * fs_315_4_5 * h3_m2 - e_2 * fs_85_4_35 * h5_m2 + e_2 * fs_70_1_5 * r_2 * h3_m2 + e_3 * fs_95_143_2002 * h7_m6 - e_3 * fs_665_143_14 * h7_m2 + e_3 * fs_255_26_35 * r_2 * h5_m2 - e_3 * fs_210_11_5 * r_4 * h3_m2 + e_4 * fs_441_4862_1430 * h9_m6 - e_4 * fs_147_4862_330 * h9_m2 - e_4 * fs_380_2431_2002 * r_2 * h7_m6 + e_4 * fs_2660_2431_14 * r_2 * h7_m2 - e_4 * fs_17_13_35 * r_4 * h5_m2 + e_4 * fs_280_143_5 * r_6 * h3_m2 - e_5 * fs_21_2431_1430 * r_2 * h9_m6 + e_5 * fs_7_2431_330 * r_2 * h9_m2 + e_5 * fs_20_2431_2002 * r_4 * h7_m6 - e_5 * fs_140_2431_14 * r_4 * h7_m2 + e_5 * fs_2_39_35 * r_6 * h5_m2 - e_5 * fs_28_429_5 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m5, ph5_m3, ph5_p4, ph7_m5, ph7_m3, ph7_p4, ph9_m5, ph9_m3, ph9_p4, ab_2, pc_92, pc_93 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p4 = ph9_p4[k];

        pc_92[k] = e_1 * fs_45_8_210 * h3_m3 + e_2 * fs_85_4_6 * h5_m5 + e_2 * fs_85_4_30 * h5_m3 - e_2 * fs_5_1_210 * r_2 * h3_m3 + e_3 * fs_1330_143_11 * h7_m5 + e_3 * f_3325_143 * h7_m3 - e_3 * fs_255_26_6 * r_2 * h5_m5 - e_3 * fs_255_26_30 * r_2 * h5_m3 + e_3 * fs_15_11_210 * r_4 * h3_m3 + e_4 * fs_147_4862_6006 * h9_m5 + e_4 * fs_441_4862_110 * h9_m3 - e_4 * fs_5320_2431_11 * r_2 * h7_m5 - e_4 * f_13300_2431 * r_2 * h7_m3 + e_4 * fs_17_13_6 * r_4 * h5_m5 + e_4 * fs_17_13_30 * r_4 * h5_m3 - e_4 * fs_20_143_210 * r_6 * h3_m3 - e_5 * fs_7_2431_6006 * r_2 * h9_m5 - e_5 * fs_21_2431_110 * r_2 * h9_m3 + e_5 * fs_280_2431_11 * r_4 * h7_m5 + e_5 * f_700_2431 * r_4 * h7_m3 - e_5 * fs_2_39_6 * r_6 * h5_m5 - e_5 * fs_2_39_30 * r_6 * h5_m3 + e_5 * fs_2_429_210 * r_8 * h3_m3;

        pc_93[k] = e_2 * f_255_2 * h5_p4 + e_3 * fs_1330_429_165 * h7_p4 - e_3 * f_765_13 * r_2 * h5_p4 + e_4 * fs_441_2431_143 * h9_p4 - e_4 * fs_5320_7293_165 * r_2 * h7_p4 + e_4 * f_102_13 * r_4 * h5_p4 - e_5 * fs_42_2431_143 * r_2 * h9_p4 + e_5 * fs_280_7293_165 * r_4 * h7_p4 - e_5 * f_4_13 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ab_2, pc_94 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];

        pc_94[k] = - e_1 * fs_45_8_210 * h3_p3 - e_2 * fs_85_4_30 * h5_p3 + e_2 * fs_85_4_6 * h5_p5 + e_2 * fs_5_1_210 * r_2 * h3_p3 - e_3 * f_3325_143 * h7_p3 + e_3 * fs_1330_143_11 * h7_p5 + e_3 * fs_255_26_30 * r_2 * h5_p3 - e_3 * fs_255_26_6 * r_2 * h5_p5 - e_3 * fs_15_11_210 * r_4 * h3_p3 - e_4 * fs_441_4862_110 * h9_p3 + e_4 * fs_147_4862_6006 * h9_p5 + e_4 * f_13300_2431 * r_2 * h7_p3 - e_4 * fs_5320_2431_11 * r_2 * h7_p5 - e_4 * fs_17_13_30 * r_4 * h5_p3 + e_4 * fs_17_13_6 * r_4 * h5_p5 + e_4 * fs_20_143_210 * r_6 * h3_p3 + e_5 * fs_21_2431_110 * r_2 * h9_p3 - e_5 * fs_7_2431_6006 * r_2 * h9_p5 - e_5 * f_700_2431 * r_4 * h7_p3 + e_5 * fs_280_2431_11 * r_4 * h7_p5 + e_5 * fs_2_39_30 * r_6 * h5_p3 - e_5 * fs_2_39_6 * r_6 * h5_p5 - e_5 * fs_2_429_210 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ab_2, pc_95 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];

        pc_95[k] = e_1 * fs_315_4_5 * h3_p2 + e_2 * fs_85_4_35 * h5_p2 - e_2 * fs_70_1_5 * r_2 * h3_p2 + e_3 * fs_665_143_14 * h7_p2 + e_3 * fs_95_143_2002 * h7_p6 - e_3 * fs_255_26_35 * r_2 * h5_p2 + e_3 * fs_210_11_5 * r_4 * h3_p2 + e_4 * fs_147_4862_330 * h9_p2 + e_4 * fs_441_4862_1430 * h9_p6 - e_4 * fs_2660_2431_14 * r_2 * h7_p2 - e_4 * fs_380_2431_2002 * r_2 * h7_p6 + e_4 * fs_17_13_35 * r_4 * h5_p2 - e_4 * fs_280_143_5 * r_6 * h3_p2 - e_5 * fs_7_2431_330 * r_2 * h9_p2 - e_5 * fs_21_2431_1430 * r_2 * h9_p6 + e_5 * fs_140_2431_14 * r_4 * h7_p2 + e_5 * fs_20_2431_2002 * r_4 * h7_p6 - e_5 * fs_2_39_35 * r_6 * h5_p2 + e_5 * fs_28_429_5 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2, pc_96 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_96[k] = - e_0 * fs_1365_32_2 * h1_p1 - e_1 * fs_315_2_3 * h3_p1 + e_1 * fs_1365_16_2 * r_2 * h1_p1 - e_2 * fs_85_4_30 * h5_p1 + e_2 * fs_140_1_3 * r_2 * h3_p1 - e_2 * fs_195_4_2 * r_4 * h1_p1 - e_3 * fs_1330_429_14 * h7_p1 + e_3 * fs_665_858_858 * h7_p7 + e_3 * fs_255_26_30 * r_2 * h5_p1 - e_3 * fs_420_11_3 * r_4 * h3_p1 + e_3 * fs_65_6_2 * r_6 * h1_p1 - e_4 * fs_441_4862_10 * h9_p1 + e_4 * fs_441_2431_715 * h9_p7 + e_4 * fs_5320_7293_14 * r_2 * h7_p1 - e_4 * fs_1330_7293_858 * r_2 * h7_p7 - e_4 * fs_17_13_30 * r_4 * h5_p1 + e_4 * fs_560_143_3 * r_6 * h3_p1 - e_4 * fs_65_66_2 * r_8 * h1_p1 + e_5 * fs_21_2431_10 * r_2 * h9_p1 - e_5 * fs_42_2431_715 * r_2 * h9_p7 - e_5 * fs_280_7293_14 * r_4 * h7_p1 + e_5 * fs_70_7293_858 * r_4 * h7_p7 + e_5 * fs_2_39_30 * r_6 * h5_p1 - e_5 * fs_56_429_3 * r_8 * h3_p1 + e_5 * fs_1_33_2 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph9_p8, ab_2, pc_97 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];

        pc_97[k] = e_0 * f_4095_16 * h1_0 + e_1 * f_945_2 * h3_0 - e_1 * f_4095_8 * r_2 * h1_0 + e_2 * f_255_2 * h5_0 - e_2 * f_420_1 * r_2 * h3_0 + e_2 * f_585_2 * r_4 * h1_0 + e_3 * f_1330_143 * h7_0 - e_3 * f_765_13 * r_2 * h5_0 + e_3 * f_1260_11 * r_4 * h3_0 - e_3 * f_65_1 * r_6 * h1_0 + e_4 * f_441_2431 * h9_0 + e_4 * fs_147_2431_12155 * h9_p8 - e_4 * f_5320_2431 * r_2 * h7_0 + e_4 * f_102_13 * r_4 * h5_0 - e_4 * f_1680_143 * r_6 * h3_0 + e_4 * f_65_11 * r_8 * h1_0 - e_5 * f_42_2431 * r_2 * h9_0 - e_5 * fs_14_2431_12155 * r_2 * h9_p8 + e_5 * f_280_2431 * r_4 * h7_0 - e_5 * f_4_13 * r_6 * h5_0 + e_5 * f_56_143 * r_8 * h3_0 - e_5 * f_2_11 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 69 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ab_2, pc_98 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];

        pc_98[k] = e_0 * fs_4095_32_10 * h1_p1 + e_1 * fs_315_4_15 * h3_p1 - e_1 * fs_4095_16_10 * r_2 * h1_p1 + e_2 * fs_85_4_6 * h5_p1 - e_2 * fs_70_1_15 * r_2 * h3_p1 + e_2 * fs_585_4_10 * r_4 * h1_p1 + e_3 * fs_95_286_70 * h7_p1 - e_3 * fs_255_26_6 * r_2 * h5_p1 + e_3 * fs_210_11_15 * r_4 * h3_p1 - e_3 * fs_65_2_10 * r_6 * h1_p1 + e_4 * fs_147_4862_2 * h9_p1 + e_4 * fs_441_2431_2431 * h9_p9 - e_4 * fs_190_2431_70 * r_2 * h7_p1 + e_4 * fs_17_13_6 * r_4 * h5_p1 - e_4 * fs_280_143_15 * r_6 * h3_p1 + e_4 * fs_65_22_10 * r_8 * h1_p1 - e_5 * fs_7_2431_2 * r_2 * h9_p1 - e_5 * fs_42_2431_2431 * r_2 * h9_p9 + e_5 * fs_10_2431_70 * r_4 * h7_p1 - e_5 * fs_2_39_6 * r_6 * h5_p1 + e_5 * fs_28_429_15 * r_8 * h3_p1 - e_5 * fs_1_11_10 * r_10 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[99] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98};

    for (size_t n = 0; n < 99; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
