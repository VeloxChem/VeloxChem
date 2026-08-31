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



#include "SimdKineticEnergyRecIF.hpp"

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
compute_if_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_if_kinetic_energy: Basis functions must be of angular momenta six and three"));
    }

    if (harmonics.size() < 9)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_if_kinetic_energy: Harmonics must reach angular momentum 9"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_if_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 91 * nvalues, 0.0);

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

            const auto ff_0 = fbase * bexp * bexp * bexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * bexp * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * bexp * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * bexp * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * bexp * bexp * bexp * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_1008_2431 = 1008.0 / 2431.0;
    const auto f_100_143 = 100.0 / 143.0;
    const auto f_102_13 = 102.0 / 13.0;
    const auto f_10_429 = 10.0 / 429.0;
    const auto f_1125_11 = 1125.0 / 11.0;
    const auto f_1125_2 = 1125.0 / 2.0;
    const auto f_114_17 = 114.0 / 17.0;
    const auto f_119_13 = 119.0 / 13.0;
    const auto f_126_221 = 126.0 / 221.0;
    const auto f_12_221 = 12.0 / 221.0;
    const auto f_14_39 = 14.0 / 39.0;
    const auto f_1500_11 = 1500.0 / 11.0;
    const auto f_1500_143 = 1500.0 / 143.0;
    const auto f_150_1 = 150.0;
    const auto f_168_221 = 168.0 / 221.0;
    const auto f_1764_221 = 1764.0 / 221.0;
    const auto f_1785_26 = 1785.0 / 26.0;
    const auto f_19152_2431 = 19152.0 / 2431.0;
    const auto f_2000_143 = 2000.0 / 143.0;
    const auto f_200_429 = 200.0 / 429.0;
    const auto f_20_143 = 20.0 / 143.0;
    const auto f_225_8 = 225.0 / 8.0;
    const auto f_255_2 = 255.0 / 2.0;
    const auto f_255_4 = 255.0 / 4.0;
    const auto f_25_1 = 25.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_3375_8 = 3375.0 / 8.0;
    const auto f_375_1 = 375.0;
    const auto f_450_11 = 450.0 / 11.0;
    const auto f_4788_143 = 4788.0 / 143.0;
    const auto f_4_13 = 4.0 / 13.0;
    const auto f_500_1 = 500.0;
    const auto f_50_143 = 50.0 / 143.0;
    const auto f_51_13 = 51.0 / 13.0;
    const auto f_57_2 = 57.0 / 2.0;
    const auto f_595_4 = 595.0 / 4.0;
    const auto f_600_143 = 600.0 / 143.0;
    const auto f_675_4 = 675.0 / 4.0;
    const auto f_6_17 = 6.0 / 17.0;
    const auto f_75_11 = 75.0 / 11.0;
    const auto f_765_13 = 765.0 / 13.0;
    const auto f_765_26 = 765.0 / 26.0;
    const auto fs_100_143_105 = std::sqrt(1050000.0 / 20449.0);
    const auto fs_100_143_14 = std::sqrt(140000.0 / 20449.0);
    const auto fs_100_143_210 = std::sqrt(2100000.0 / 20449.0);
    const auto fs_100_143_231 = std::sqrt(210000.0 / 1859.0);
    const auto fs_100_143_462 = std::sqrt(420000.0 / 1859.0);
    const auto fs_100_1_14 = std::sqrt(140000.0);
    const auto fs_1026_2431_165 = std::sqrt(15790140.0 / 537251.0);
    const auto fs_105_221_143 = std::sqrt(121275.0 / 3757.0);
    const auto fs_105_221_231 = std::sqrt(2546775.0 / 48841.0);
    const auto fs_105_442_154 = std::sqrt(848925.0 / 97682.0);
    const auto fs_1083_2431_110 = std::sqrt(11728890.0 / 537251.0);
    const auto fs_1083_572_110 = std::sqrt(5864445.0 / 14872.0);
    const auto fs_10887_2431_2 = std::sqrt(237053538.0 / 5909761.0);
    const auto fs_10887_572_2 = std::sqrt(118526769.0 / 163592.0);
    const auto fs_10_143_21 = std::sqrt(2100.0 / 20449.0);
    const auto fs_10_221_143 = std::sqrt(1100.0 / 3757.0);
    const auto fs_10_221_231 = std::sqrt(23100.0 / 48841.0);
    const auto fs_10_429_105 = std::sqrt(3500.0 / 61347.0);
    const auto fs_10_429_14 = std::sqrt(1400.0 / 184041.0);
    const auto fs_10_429_210 = std::sqrt(7000.0 / 61347.0);
    const auto fs_10_429_231 = std::sqrt(700.0 / 5577.0);
    const auto fs_10_429_462 = std::sqrt(1400.0 / 5577.0);
    const auto fs_111_2431_30 = std::sqrt(369630.0 / 5909761.0);
    const auto fs_1125_8_14 = std::sqrt(8859375.0 / 32.0);
    const auto fs_1140_143_10 = std::sqrt(12996000.0 / 20449.0);
    const auto fs_114_143_1430 = std::sqrt(129960.0 / 143.0);
    const auto fs_114_187_42 = std::sqrt(545832.0 / 34969.0);
    const auto fs_114_221_273 = std::sqrt(272916.0 / 3757.0);
    const auto fs_114_221_30 = std::sqrt(389880.0 / 48841.0);
    const auto fs_1197_143_5 = std::sqrt(7164045.0 / 20449.0);
    const auto fs_119_26_2 = std::sqrt(14161.0 / 338.0);
    const auto fs_125_1_14 = std::sqrt(218750.0);
    const auto fs_126_221_66 = std::sqrt(1047816.0 / 48841.0);
    const auto fs_12_221_66 = std::sqrt(9504.0 / 48841.0);
    const auto fs_12_2431_14 = std::sqrt(2016.0 / 5909761.0);
    const auto fs_1311_286_5 = std::sqrt(8593605.0 / 81796.0);
    const auto fs_135_2431_10 = std::sqrt(182250.0 / 5909761.0);
    const auto fs_1368_2431_70 = std::sqrt(130999680.0 / 5909761.0);
    const auto fs_138_2431_5 = std::sqrt(95220.0 / 5909761.0);
    const auto fs_1425_143_14 = std::sqrt(28428750.0 / 20449.0);
    const auto fs_1425_2431_66 = std::sqrt(12183750.0 / 537251.0);
    const auto fs_1425_572_66 = std::sqrt(6091875.0 / 14872.0);
    const auto fs_147_221_66 = std::sqrt(1426194.0 / 48841.0);
    const auto fs_14_221_66 = std::sqrt(12936.0 / 48841.0);
    const auto fs_150_11_21 = std::sqrt(472500.0 / 121.0);
    const auto fs_150_1_7 = std::sqrt(157500.0);
    const auto fs_15_221_66 = std::sqrt(14850.0 / 48841.0);
    const auto fs_15_2431_330 = std::sqrt(6750.0 / 537251.0);
    const auto fs_1653_2431_22 = std::sqrt(5464818.0 / 537251.0);
    const auto fs_1653_572_22 = std::sqrt(2732409.0 / 14872.0);
    const auto fs_168_221_105 = std::sqrt(2963520.0 / 48841.0);
    const auto fs_16_221_105 = std::sqrt(26880.0 / 48841.0);
    const auto fs_1710_2431_110 = std::sqrt(29241000.0 / 537251.0);
    const auto fs_1710_2431_210 = std::sqrt(614061000.0 / 5909761.0);
    const auto fs_171_13_13 = std::sqrt(29241.0 / 13.0);
    const auto fs_171_13_3 = std::sqrt(87723.0 / 169.0);
    const auto fs_171_286_10 = std::sqrt(146205.0 / 40898.0);
    const auto fs_171_286_55 = std::sqrt(146205.0 / 7436.0);
    const auto fs_1785_52_2 = std::sqrt(3186225.0 / 1352.0);
    const auto fs_17_13_14 = std::sqrt(4046.0 / 169.0);
    const auto fs_17_13_15 = std::sqrt(4335.0 / 169.0);
    const auto fs_17_13_21 = std::sqrt(6069.0 / 169.0);
    const auto fs_17_13_3 = std::sqrt(867.0 / 169.0);
    const auto fs_17_13_33 = std::sqrt(9537.0 / 169.0);
    const auto fs_17_13_35 = std::sqrt(10115.0 / 169.0);
    const auto fs_17_13_42 = std::sqrt(12138.0 / 169.0);
    const auto fs_17_13_6 = std::sqrt(1734.0 / 169.0);
    const auto fs_17_26_30 = std::sqrt(4335.0 / 338.0);
    const auto fs_17_26_6 = std::sqrt(867.0 / 338.0);
    const auto fs_17_26_66 = std::sqrt(9537.0 / 338.0);
    const auto fs_180_2431_21 = std::sqrt(680400.0 / 5909761.0);
    const auto fs_1824_2431_33 = std::sqrt(9980928.0 / 537251.0);
    const auto fs_189_221_14 = std::sqrt(500094.0 / 48841.0);
    const auto fs_18_221_14 = std::sqrt(4536.0 / 48841.0);
    const auto fs_18_2431_10 = std::sqrt(3240.0 / 5909761.0);
    const auto fs_18_2431_55 = std::sqrt(1620.0 / 537251.0);
    const auto fs_192_2431_15 = std::sqrt(552960.0 / 5909761.0);
    const auto fs_1995_143_3 = std::sqrt(11940075.0 / 20449.0);
    const auto fs_1_13_10 = std::sqrt(10.0 / 169.0);
    const auto fs_1_13_22 = std::sqrt(22.0 / 169.0);
    const auto fs_1_221_10010 = std::sqrt(770.0 / 3757.0);
    const auto fs_1_221_14 = std::sqrt(14.0 / 48841.0);
    const auto fs_1_221_182 = std::sqrt(14.0 / 3757.0);
    const auto fs_1_221_2 = std::sqrt(2.0 / 48841.0);
    const auto fs_1_221_26 = std::sqrt(2.0 / 3757.0);
    const auto fs_1_221_462 = std::sqrt(462.0 / 48841.0);
    const auto fs_1_221_6006 = std::sqrt(462.0 / 3757.0);
    const auto fs_1_39_30 = std::sqrt(10.0 / 507.0);
    const auto fs_1_39_6 = std::sqrt(2.0 / 507.0);
    const auto fs_1_39_66 = std::sqrt(22.0 / 507.0);
    const auto fs_200_143_21 = std::sqrt(840000.0 / 20449.0);
    const auto fs_20_143_7 = std::sqrt(2800.0 / 20449.0);
    const auto fs_20_429_21 = std::sqrt(2800.0 / 61347.0);
    const auto fs_2109_2431_30 = std::sqrt(133436430.0 / 5909761.0);
    const auto fs_2109_572_30 = std::sqrt(66718215.0 / 163592.0);
    const auto fs_21_221_105 = std::sqrt(46305.0 / 48841.0);
    const auto fs_21_221_14 = std::sqrt(6174.0 / 48841.0);
    const auto fs_21_221_195 = std::sqrt(6615.0 / 3757.0);
    const auto fs_21_221_231 = std::sqrt(101871.0 / 48841.0);
    const auto fs_21_221_3003 = std::sqrt(101871.0 / 3757.0);
    const auto fs_21_221_3094 = std::sqrt(6174.0 / 221.0);
    const auto fs_21_221_455 = std::sqrt(15435.0 / 3757.0);
    const auto fs_21_221_858 = std::sqrt(29106.0 / 3757.0);
    const auto fs_21_221_910 = std::sqrt(30870.0 / 3757.0);
    const auto fs_21_221_9282 = std::sqrt(18522.0 / 221.0);
    const auto fs_21_442_10010 = std::sqrt(169785.0 / 7514.0);
    const auto fs_21_442_14 = std::sqrt(3087.0 / 97682.0);
    const auto fs_21_442_182 = std::sqrt(3087.0 / 7514.0);
    const auto fs_21_442_2 = std::sqrt(441.0 / 97682.0);
    const auto fs_21_442_26 = std::sqrt(441.0 / 7514.0);
    const auto fs_21_442_462 = std::sqrt(101871.0 / 97682.0);
    const auto fs_21_442_6006 = std::sqrt(101871.0 / 7514.0);
    const auto fs_225_11_21 = std::sqrt(1063125.0 / 121.0);
    const auto fs_225_16_14 = std::sqrt(354375.0 / 128.0);
    const auto fs_225_16_210 = std::sqrt(5315625.0 / 128.0);
    const auto fs_225_2_14 = std::sqrt(354375.0 / 2.0);
    const auto fs_225_4_21 = std::sqrt(1063125.0 / 16.0);
    const auto fs_225_8_105 = std::sqrt(5315625.0 / 64.0);
    const auto fs_225_8_14 = std::sqrt(354375.0 / 32.0);
    const auto fs_225_8_210 = std::sqrt(5315625.0 / 32.0);
    const auto fs_225_8_231 = std::sqrt(11694375.0 / 64.0);
    const auto fs_225_8_462 = std::sqrt(11694375.0 / 32.0);
    const auto fs_228_2431_14 = std::sqrt(727776.0 / 5909761.0);
    const auto fs_240_2431_10 = std::sqrt(576000.0 / 5909761.0);
    const auto fs_2451_143_3 = std::sqrt(18022203.0 / 20449.0);
    const auto fs_24_221_7 = std::sqrt(4032.0 / 48841.0);
    const auto fs_24_2431_1430 = std::sqrt(5760.0 / 41327.0);
    const auto fs_252_221_7 = std::sqrt(444528.0 / 48841.0);
    const auto fs_252_2431_5 = std::sqrt(317520.0 / 5909761.0);
    const auto fs_255_13_2 = std::sqrt(130050.0 / 169.0);
    const auto fs_255_13_21 = std::sqrt(1365525.0 / 169.0);
    const auto fs_255_13_7 = std::sqrt(455175.0 / 169.0);
    const auto fs_255_26_14 = std::sqrt(455175.0 / 338.0);
    const auto fs_255_26_15 = std::sqrt(975375.0 / 676.0);
    const auto fs_255_26_21 = std::sqrt(1365525.0 / 676.0);
    const auto fs_255_26_3 = std::sqrt(195075.0 / 676.0);
    const auto fs_255_26_33 = std::sqrt(2145825.0 / 676.0);
    const auto fs_255_26_35 = std::sqrt(2275875.0 / 676.0);
    const auto fs_255_26_42 = std::sqrt(1365525.0 / 338.0);
    const auto fs_255_26_6 = std::sqrt(195075.0 / 338.0);
    const auto fs_255_4_11 = std::sqrt(715275.0 / 16.0);
    const auto fs_255_4_3 = std::sqrt(195075.0 / 16.0);
    const auto fs_255_4_5 = std::sqrt(325125.0 / 16.0);
    const auto fs_255_4_7 = std::sqrt(455175.0 / 16.0);
    const auto fs_255_52_30 = std::sqrt(975375.0 / 1352.0);
    const auto fs_255_52_6 = std::sqrt(195075.0 / 1352.0);
    const auto fs_255_52_66 = std::sqrt(2145825.0 / 1352.0);
    const auto fs_255_8_10 = std::sqrt(325125.0 / 32.0);
    const auto fs_255_8_22 = std::sqrt(715275.0 / 32.0);
    const auto fs_2565_2431_10 = std::sqrt(65792250.0 / 5909761.0);
    const auto fs_2565_286_11 = std::sqrt(6579225.0 / 7436.0);
    const auto fs_2565_572_10 = std::sqrt(32896125.0 / 163592.0);
    const auto fs_25_1_105 = std::sqrt(65625.0);
    const auto fs_25_1_14 = std::sqrt(8750.0);
    const auto fs_25_1_210 = std::sqrt(131250.0);
    const auto fs_25_1_231 = std::sqrt(144375.0);
    const auto fs_25_1_462 = std::sqrt(288750.0);
    const auto fs_25_2_14 = std::sqrt(4375.0 / 2.0);
    const auto fs_25_2_210 = std::sqrt(65625.0 / 2.0);
    const auto fs_2622_2431_5 = std::sqrt(34374420.0 / 5909761.0);
    const auto fs_270_2431_11 = std::sqrt(72900.0 / 537251.0);
    const auto fs_2793_286_5 = std::sqrt(39004245.0 / 81796.0);
    const auto fs_285_2431_330 = std::sqrt(2436750.0 / 537251.0);
    const auto fs_285_286_15 = std::sqrt(1218375.0 / 81796.0);
    const auto fs_285_286_66 = std::sqrt(243675.0 / 3718.0);
    const auto fs_285_572_330 = std::sqrt(1218375.0 / 14872.0);
    const auto fs_28_221_30 = std::sqrt(23520.0 / 48841.0);
    const auto fs_294_221_30 = std::sqrt(2593080.0 / 48841.0);
    const auto fs_294_2431_5 = std::sqrt(432180.0 / 5909761.0);
    const auto fs_2_13_11 = std::sqrt(44.0 / 169.0);
    const auto fs_2_13_3 = std::sqrt(12.0 / 169.0);
    const auto fs_2_13_5 = std::sqrt(20.0 / 169.0);
    const auto fs_2_13_7 = std::sqrt(28.0 / 169.0);
    const auto fs_2_221_105 = std::sqrt(420.0 / 48841.0);
    const auto fs_2_221_14 = std::sqrt(56.0 / 48841.0);
    const auto fs_2_221_195 = std::sqrt(60.0 / 3757.0);
    const auto fs_2_221_231 = std::sqrt(924.0 / 48841.0);
    const auto fs_2_221_3003 = std::sqrt(924.0 / 3757.0);
    const auto fs_2_221_3094 = std::sqrt(56.0 / 221.0);
    const auto fs_2_221_455 = std::sqrt(140.0 / 3757.0);
    const auto fs_2_221_858 = std::sqrt(264.0 / 3757.0);
    const auto fs_2_221_910 = std::sqrt(280.0 / 3757.0);
    const auto fs_2_221_9282 = std::sqrt(168.0 / 221.0);
    const auto fs_2_39_14 = std::sqrt(56.0 / 1521.0);
    const auto fs_2_39_15 = std::sqrt(20.0 / 507.0);
    const auto fs_2_39_21 = std::sqrt(28.0 / 507.0);
    const auto fs_2_39_3 = std::sqrt(4.0 / 507.0);
    const auto fs_2_39_33 = std::sqrt(44.0 / 507.0);
    const auto fs_2_39_35 = std::sqrt(140.0 / 1521.0);
    const auto fs_2_39_42 = std::sqrt(56.0 / 507.0);
    const auto fs_2_39_6 = std::sqrt(8.0 / 507.0);
    const auto fs_300_11_14 = std::sqrt(1260000.0 / 121.0);
    const auto fs_300_143_21 = std::sqrt(1890000.0 / 20449.0);
    const auto fs_300_2431_14 = std::sqrt(1260000.0 / 5909761.0);
    const auto fs_30_221_7 = std::sqrt(6300.0 / 48841.0);
    const auto fs_30_2431_15 = std::sqrt(13500.0 / 5909761.0);
    const auto fs_30_2431_66 = std::sqrt(5400.0 / 537251.0);
    const auto fs_315_221_7 = std::sqrt(694575.0 / 48841.0);
    const auto fs_315_442_66 = std::sqrt(3274425.0 / 97682.0);
    const auto fs_3363_286_3 = std::sqrt(33929307.0 / 81796.0);
    const auto fs_3420_2431_21 = std::sqrt(245624400.0 / 5909761.0);
    const auto fs_342_143_70 = std::sqrt(8187480.0 / 20449.0);
    const auto fs_342_2431_10 = std::sqrt(1169640.0 / 5909761.0);
    const auto fs_342_2431_55 = std::sqrt(584820.0 / 537251.0);
    const auto fs_34_13_2 = std::sqrt(2312.0 / 169.0);
    const auto fs_34_13_21 = std::sqrt(24276.0 / 169.0);
    const auto fs_34_13_7 = std::sqrt(8092.0 / 169.0);
    const auto fs_354_2431_3 = std::sqrt(375948.0 / 5909761.0);
    const auto fs_3591_286_5 = std::sqrt(64476405.0 / 81796.0);
    const auto fs_3648_2431_15 = std::sqrt(199618560.0 / 5909761.0);
    const auto fs_36_221_13 = std::sqrt(1296.0 / 3757.0);
    const auto fs_36_221_14 = std::sqrt(18144.0 / 48841.0);
    const auto fs_36_221_3 = std::sqrt(3888.0 / 48841.0);
    const auto fs_375_11_14 = std::sqrt(1968750.0 / 121.0);
    const auto fs_378_221_14 = std::sqrt(2000376.0 / 48841.0);
    const auto fs_378_2431_5 = std::sqrt(714420.0 / 5909761.0);
    const auto fs_399_26_6 = std::sqrt(477603.0 / 338.0);
    const auto fs_3_221_1430 = std::sqrt(990.0 / 3757.0);
    const auto fs_3_221_2002 = std::sqrt(1386.0 / 3757.0);
    const auto fs_3_221_26 = std::sqrt(18.0 / 3757.0);
    const auto fs_3_221_910 = std::sqrt(630.0 / 3757.0);
    const auto fs_3_2431_30030 = std::sqrt(1890.0 / 41327.0);
    const auto fs_400_143_14 = std::sqrt(2240000.0 / 20449.0);
    const auto fs_40_221_11 = std::sqrt(17600.0 / 48841.0);
    const auto fs_40_429_14 = std::sqrt(22400.0 / 184041.0);
    const auto fs_420_221_11 = std::sqrt(1940400.0 / 48841.0);
    const auto fs_420_2431_3 = std::sqrt(529200.0 / 5909761.0);
    const auto fs_42_221_1001 = std::sqrt(135828.0 / 3757.0);
    const auto fs_42_221_1547 = std::sqrt(12348.0 / 221.0);
    const auto fs_42_221_21 = std::sqrt(37044.0 / 48841.0);
    const auto fs_42_221_231 = std::sqrt(407484.0 / 48841.0);
    const auto fs_42_221_6 = std::sqrt(10584.0 / 48841.0);
    const auto fs_42_221_70 = std::sqrt(123480.0 / 48841.0);
    const auto fs_42_221_715 = std::sqrt(97020.0 / 3757.0);
    const auto fs_450_11_7 = std::sqrt(1417500.0 / 121.0);
    const auto fs_4560_2431_10 = std::sqrt(207936000.0 / 5909761.0);
    const auto fs_456_143_33 = std::sqrt(623808.0 / 1859.0);
    const auto fs_456_2431_1430 = std::sqrt(2079360.0 / 41327.0);
    const auto fs_45_187_2 = std::sqrt(4050.0 / 34969.0);
    const auto fs_45_2431_286 = std::sqrt(4050.0 / 41327.0);
    const auto fs_4788_2431_5 = std::sqrt(114624720.0 / 5909761.0);
    const auto fs_4_221_1001 = std::sqrt(1232.0 / 3757.0);
    const auto fs_4_221_1547 = std::sqrt(112.0 / 221.0);
    const auto fs_4_221_21 = std::sqrt(336.0 / 48841.0);
    const auto fs_4_221_231 = std::sqrt(3696.0 / 48841.0);
    const auto fs_4_221_70 = std::sqrt(1120.0 / 48841.0);
    const auto fs_4_221_715 = std::sqrt(880.0 / 3757.0);
    const auto fs_4_39_2 = std::sqrt(32.0 / 1521.0);
    const auto fs_4_39_21 = std::sqrt(112.0 / 507.0);
    const auto fs_4_39_7 = std::sqrt(112.0 / 1521.0);
    const auto fs_500_143_14 = std::sqrt(3500000.0 / 20449.0);
    const auto fs_50_143_14 = std::sqrt(35000.0 / 20449.0);
    const auto fs_50_143_210 = std::sqrt(525000.0 / 20449.0);
    const auto fs_50_1_21 = std::sqrt(52500.0);
    const auto fs_50_429_14 = std::sqrt(35000.0 / 184041.0);
    const auto fs_510_13_2 = std::sqrt(520200.0 / 169.0);
    const auto fs_5130_2431_11 = std::sqrt(26316900.0 / 537251.0);
    const auto fs_513_286_165 = std::sqrt(3947535.0 / 7436.0);
    const auto fs_516_2431_3 = std::sqrt(798768.0 / 5909761.0);
    const auto fs_51_13_11 = std::sqrt(28611.0 / 169.0);
    const auto fs_51_13_3 = std::sqrt(7803.0 / 169.0);
    const auto fs_51_13_5 = std::sqrt(13005.0 / 169.0);
    const auto fs_51_13_7 = std::sqrt(18207.0 / 169.0);
    const auto fs_51_26_10 = std::sqrt(13005.0 / 338.0);
    const auto fs_51_26_22 = std::sqrt(28611.0 / 338.0);
    const auto fs_54_2431_165 = std::sqrt(43740.0 / 537251.0);
    const auto fs_5586_2431_5 = std::sqrt(156016980.0 / 5909761.0);
    const auto fs_5700_2431_14 = std::sqrt(454860000.0 / 5909761.0);
    const auto fs_570_2431_15 = std::sqrt(4873500.0 / 5909761.0);
    const auto fs_570_2431_66 = std::sqrt(1949400.0 / 537251.0);
    const auto fs_573_2431_2 = std::sqrt(656658.0 / 5909761.0);
    const auto fs_57_143_14 = std::sqrt(45486.0 / 20449.0);
    const auto fs_57_221_26 = std::sqrt(6498.0 / 3757.0);
    const auto fs_57_221_910 = std::sqrt(227430.0 / 3757.0);
    const auto fs_57_22_42 = std::sqrt(68229.0 / 242.0);
    const auto fs_57_2431_110 = std::sqrt(32490.0 / 537251.0);
    const auto fs_57_2431_30030 = std::sqrt(682290.0 / 41327.0);
    const auto fs_57_26_273 = std::sqrt(68229.0 / 52.0);
    const auto fs_57_26_30 = std::sqrt(48735.0 / 338.0);
    const auto fs_57_52_26 = std::sqrt(3249.0 / 104.0);
    const auto fs_57_52_910 = std::sqrt(113715.0 / 104.0);
    const auto fs_57_572_30030 = std::sqrt(341145.0 / 1144.0);
    const auto fs_595_8_2 = std::sqrt(354025.0 / 32.0);
    const auto fs_5_221_154 = std::sqrt(3850.0 / 48841.0);
    const auto fs_5_429_14 = std::sqrt(350.0 / 184041.0);
    const auto fs_5_429_210 = std::sqrt(1750.0 / 61347.0);
    const auto fs_600_143_7 = std::sqrt(2520000.0 / 20449.0);
    const auto fs_63_221_385 = std::sqrt(1528065.0 / 48841.0);
    const auto fs_63_221_55 = std::sqrt(218295.0 / 48841.0);
    const auto fs_63_442_1430 = std::sqrt(218295.0 / 7514.0);
    const auto fs_63_442_2002 = std::sqrt(305613.0 / 7514.0);
    const auto fs_63_442_910 = std::sqrt(138915.0 / 7514.0);
    const auto fs_6726_2431_3 = std::sqrt(135717228.0 / 5909761.0);
    const auto fs_675_4_7 = std::sqrt(3189375.0 / 16.0);
    const auto fs_675_8_21 = std::sqrt(9568125.0 / 64.0);
    const auto fs_684_221_13 = std::sqrt(467856.0 / 3757.0);
    const auto fs_684_221_3 = std::sqrt(1403568.0 / 48841.0);
    const auto fs_68_13_2 = std::sqrt(9248.0 / 169.0);
    const auto fs_6_187_42 = std::sqrt(1512.0 / 34969.0);
    const auto fs_6_221_273 = std::sqrt(756.0 / 3757.0);
    const auto fs_6_221_30 = std::sqrt(1080.0 / 48841.0);
    const auto fs_6_221_385 = std::sqrt(13860.0 / 48841.0);
    const auto fs_6_221_55 = std::sqrt(1980.0 / 48841.0);
    const auto fs_7182_2431_5 = std::sqrt(257905620.0 / 5909761.0);
    const auto fs_72_2431_70 = std::sqrt(362880.0 / 5909761.0);
    const auto fs_75_11_105 = std::sqrt(590625.0 / 121.0);
    const auto fs_75_11_14 = std::sqrt(78750.0 / 121.0);
    const auto fs_75_11_210 = std::sqrt(1181250.0 / 121.0);
    const auto fs_75_11_231 = std::sqrt(118125.0 / 11.0);
    const auto fs_75_11_462 = std::sqrt(236250.0 / 11.0);
    const auto fs_75_1_21 = std::sqrt(118125.0);
    const auto fs_75_22_14 = std::sqrt(39375.0 / 242.0);
    const auto fs_75_22_210 = std::sqrt(590625.0 / 242.0);
    const auto fs_75_2431_66 = std::sqrt(33750.0 / 537251.0);
    const auto fs_765_26_11 = std::sqrt(6437475.0 / 676.0);
    const auto fs_765_26_3 = std::sqrt(1755675.0 / 676.0);
    const auto fs_765_26_5 = std::sqrt(2926125.0 / 676.0);
    const auto fs_765_26_7 = std::sqrt(4096575.0 / 676.0);
    const auto fs_765_52_10 = std::sqrt(2926125.0 / 1352.0);
    const auto fs_765_52_22 = std::sqrt(6437475.0 / 1352.0);
    const auto fs_7980_2431_3 = std::sqrt(191041200.0 / 5909761.0);
    const auto fs_798_221_6 = std::sqrt(3820824.0 / 48841.0);
    const auto fs_7_39_2 = std::sqrt(98.0 / 1521.0);
    const auto fs_84_221_273 = std::sqrt(148176.0 / 3757.0);
    const auto fs_84_221_91 = std::sqrt(49392.0 / 3757.0);
    const auto fs_855_143_21 = std::sqrt(15351525.0 / 20449.0);
    const auto fs_855_187_2 = std::sqrt(1462050.0 / 34969.0);
    const auto fs_855_2431_286 = std::sqrt(1462050.0 / 41327.0);
    const auto fs_855_286_110 = std::sqrt(3655125.0 / 3718.0);
    const auto fs_855_286_210 = std::sqrt(76757625.0 / 40898.0);
    const auto fs_855_44_2 = std::sqrt(731025.0 / 968.0);
    const auto fs_855_572_286 = std::sqrt(731025.0 / 1144.0);
    const auto fs_85_1_2 = std::sqrt(14450.0);
    const auto fs_85_2_2 = std::sqrt(7225.0 / 2.0);
    const auto fs_85_2_21 = std::sqrt(151725.0 / 4.0);
    const auto fs_85_2_7 = std::sqrt(50575.0 / 4.0);
    const auto fs_85_4_14 = std::sqrt(50575.0 / 8.0);
    const auto fs_85_4_15 = std::sqrt(108375.0 / 16.0);
    const auto fs_85_4_21 = std::sqrt(151725.0 / 16.0);
    const auto fs_85_4_3 = std::sqrt(21675.0 / 16.0);
    const auto fs_85_4_33 = std::sqrt(238425.0 / 16.0);
    const auto fs_85_4_35 = std::sqrt(252875.0 / 16.0);
    const auto fs_85_4_42 = std::sqrt(151725.0 / 8.0);
    const auto fs_85_4_6 = std::sqrt(21675.0 / 8.0);
    const auto fs_85_8_30 = std::sqrt(108375.0 / 32.0);
    const auto fs_85_8_6 = std::sqrt(21675.0 / 32.0);
    const auto fs_85_8_66 = std::sqrt(238425.0 / 32.0);
    const auto fs_87_2431_22 = std::sqrt(15138.0 / 537251.0);
    const auto fs_8_221_273 = std::sqrt(1344.0 / 3757.0);
    const auto fs_8_221_91 = std::sqrt(448.0 / 3757.0);
    const auto fs_8_39_2 = std::sqrt(128.0 / 1521.0);
    const auto fs_90_2431_110 = std::sqrt(81000.0 / 537251.0);
    const auto fs_90_2431_210 = std::sqrt(1701000.0 / 5909761.0);
    const auto fs_912_143_15 = std::sqrt(12476160.0 / 20449.0);
    const auto fs_96_2431_33 = std::sqrt(27648.0 / 537251.0);
    const auto fs_9804_2431_3 = std::sqrt(288355248.0 / 5909761.0);

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph9_p3, ph9_p4, ph9_p8, ph9_p9, ab_2, pc_0, pc_1 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h9_p9 = ph9_p9[k];

        pc_0[k] = - e_0 * fs_225_8_462 * h3_p3 - e_1 * fs_85_8_66 * h5_p3 + e_1 * fs_25_1_462 * r_2 * h3_p3 - e_2 * fs_171_286_55 * h7_p3 + e_2 * fs_255_52_66 * r_2 * h5_p3 - e_2 * fs_75_11_462 * r_4 * h3_p3 - e_3 * fs_21_442_2 * h9_p3 - e_3 * fs_21_221_9282 * h9_p9 + e_3 * fs_342_2431_55 * r_2 * h7_p3 - e_3 * fs_17_26_66 * r_4 * h5_p3 + e_3 * fs_100_143_462 * r_6 * h3_p3 + e_4 * fs_1_221_2 * r_2 * h9_p3 + e_4 * fs_2_221_9282 * r_2 * h9_p9 - e_4 * fs_18_2431_55 * r_4 * h7_p3 + e_4 * fs_1_39_66 * r_6 * h5_p3 - e_4 * fs_10_429_462 * r_8 * h3_p3;

        pc_1[k] = e_1 * fs_255_8_22 * h5_p4 + e_2 * fs_57_26_30 * h7_p4 - e_2 * fs_765_52_22 * r_2 * h5_p4 + e_3 * fs_21_442_26 * h9_p4 - e_3 * fs_21_221_3094 * h9_p8 - e_3 * fs_114_221_30 * r_2 * h7_p4 + e_3 * fs_51_26_22 * r_4 * h5_p4 - e_4 * fs_1_221_26 * r_2 * h9_p4 + e_4 * fs_2_221_3094 * r_2 * h9_p8 + e_4 * fs_6_221_30 * r_4 * h7_p4 - e_4 * fs_1_13_22 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m5, ph5_p5, ph7_m7, ph7_m6, ph7_m5, ph7_p5, ph7_p7, ph9_m7, ph9_m6, ph9_m5, ph9_p5, ph9_p7, ab_2, pc_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m5 = ph5_m5[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p7 = ph9_p7[k];

        pc_2[k] = - e_1 * fs_255_8_22 * h5_p5 - e_2 * fs_171_13_3 * h7_p5 - e_2 * fs_57_26_273 * h7_p7 + e_2 * fs_765_52_22 * r_2 * h5_p5 - e_3 * fs_21_442_182 * h9_p5 - e_3 * fs_21_221_910 * h9_p7 + e_3 * fs_684_221_3 * r_2 * h7_p5 + e_3 * fs_114_221_273 * r_2 * h7_p7 - e_3 * fs_51_26_22 * r_4 * h5_p5 + e_4 * fs_1_221_182 * r_2 * h9_p5 + e_4 * fs_2_221_910 * r_2 * h9_p7 - e_4 * fs_36_221_3 * r_4 * h7_p5 - e_4 * fs_6_221_273 * r_4 * h7_p7 + e_4 * fs_1_13_22 * r_6 * h5_p5;

        pc_3[k] = e_2 * fs_171_13_13 * h7_m6 + e_3 * fs_21_221_455 * h9_m6 - e_3 * fs_684_221_13 * r_2 * h7_m6 - e_4 * fs_2_221_455 * r_2 * h9_m6 + e_4 * fs_36_221_13 * r_4 * h7_m6;

        pc_4[k] = - e_1 * fs_255_8_22 * h5_m5 + e_2 * fs_57_26_273 * h7_m7 - e_2 * fs_171_13_3 * h7_m5 + e_2 * fs_765_52_22 * r_2 * h5_m5 + e_3 * fs_21_221_910 * h9_m7 - e_3 * fs_21_442_182 * h9_m5 - e_3 * fs_114_221_273 * r_2 * h7_m7 + e_3 * fs_684_221_3 * r_2 * h7_m5 - e_3 * fs_51_26_22 * r_4 * h5_m5 - e_4 * fs_2_221_910 * r_2 * h9_m7 + e_4 * fs_1_221_182 * r_2 * h9_m5 + e_4 * fs_6_221_273 * r_4 * h7_m7 - e_4 * fs_36_221_3 * r_4 * h7_m5 + e_4 * fs_1_13_22 * r_6 * h5_m5;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m4, ph5_m3, ph7_m4, ph7_m3, ph9_m9, ph9_m8, ph9_m4, ph9_m3, ab_2, pc_5, pc_6 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_5[k] = e_1 * fs_255_8_22 * h5_m4 + e_2 * fs_57_26_30 * h7_m4 - e_2 * fs_765_52_22 * r_2 * h5_m4 + e_3 * fs_21_221_3094 * h9_m8 + e_3 * fs_21_442_26 * h9_m4 - e_3 * fs_114_221_30 * r_2 * h7_m4 + e_3 * fs_51_26_22 * r_4 * h5_m4 - e_4 * fs_2_221_3094 * r_2 * h9_m8 - e_4 * fs_1_221_26 * r_2 * h9_m4 + e_4 * fs_6_221_30 * r_4 * h7_m4 - e_4 * fs_1_13_22 * r_6 * h5_m4;

        pc_6[k] = - e_0 * fs_225_8_462 * h3_m3 - e_1 * fs_85_8_66 * h5_m3 + e_1 * fs_25_1_462 * r_2 * h3_m3 - e_2 * fs_171_286_55 * h7_m3 + e_2 * fs_255_52_66 * r_2 * h5_m3 - e_2 * fs_75_11_462 * r_4 * h3_m3 + e_3 * fs_21_221_9282 * h9_m9 - e_3 * fs_21_442_2 * h9_m3 + e_3 * fs_342_2431_55 * r_2 * h7_m3 - e_3 * fs_17_26_66 * r_4 * h5_m3 + e_3 * fs_100_143_462 * r_6 * h3_m3 - e_4 * fs_2_221_9282 * r_2 * h9_m9 + e_4 * fs_1_221_2 * r_2 * h9_m3 - e_4 * fs_18_2431_55 * r_4 * h7_m3 + e_4 * fs_1_39_66 * r_6 * h5_m3 - e_4 * fs_10_429_462 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph7_p7, ph9_p2, ph9_p3, ph9_p7, ph9_p8, ab_2, pc_7, pc_8 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h9_p8 = ph9_p8[k];

        pc_7[k] = - e_0 * fs_225_8_231 * h3_p2 - e_1 * fs_85_4_33 * h5_p2 + e_1 * fs_25_1_231 * r_2 * h3_p2 - e_2 * fs_285_572_330 * h7_p2 + e_2 * fs_255_26_33 * r_2 * h5_p2 - e_2 * fs_75_11_231 * r_4 * h3_p2 - e_3 * fs_21_442_14 * h9_p2 - e_3 * fs_42_221_1547 * h9_p8 + e_3 * fs_285_2431_330 * r_2 * h7_p2 - e_3 * fs_17_13_33 * r_4 * h5_p2 + e_3 * fs_100_143_231 * r_6 * h3_p2 + e_4 * fs_1_221_14 * r_2 * h9_p2 + e_4 * fs_4_221_1547 * r_2 * h9_p8 - e_4 * fs_15_2431_330 * r_4 * h7_p2 + e_4 * fs_2_39_33 * r_6 * h5_p2 - e_4 * fs_10_429_231 * r_8 * h3_p2;

        pc_8[k] = - e_0 * fs_225_8_231 * h3_p3 + e_1 * fs_85_4_33 * h5_p3 + e_1 * fs_25_1_231 * r_2 * h3_p3 + e_2 * fs_1083_572_110 * h7_p3 + e_2 * fs_57_52_910 * h7_p7 - e_2 * fs_255_26_33 * r_2 * h5_p3 - e_2 * fs_75_11_231 * r_4 * h3_p3 + e_3 * f_126_221 * h9_p3 - e_3 * fs_84_221_273 * h9_p7 - e_3 * fs_1083_2431_110 * r_2 * h7_p3 - e_3 * fs_57_221_910 * r_2 * h7_p7 + e_3 * fs_17_13_33 * r_4 * h5_p3 + e_3 * fs_100_143_231 * r_6 * h3_p3 - e_4 * f_12_221 * r_2 * h9_p3 + e_4 * fs_8_221_273 * r_2 * h9_p7 + e_4 * fs_57_2431_110 * r_4 * h7_p3 + e_4 * fs_3_221_910 * r_4 * h7_p7 - e_4 * fs_2_39_33 * r_6 * h5_p3 - e_4 * fs_10_429_231 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m5, ph7_m6, ph7_m5, ph7_m4, ph7_p4, ph7_p6, ph9_m6, ph9_m5, ph9_m4, ph9_p4, ph9_p6, ab_2, pc_9, pc_10, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m5 = ph5_m5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];

        pc_9[k] = - e_2 * f_57_2 * h7_p4 - e_2 * fs_57_52_26 * h7_p6 - e_3 * fs_21_221_195 * h9_p4 - e_3 * fs_63_442_910 * h9_p6 + e_3 * f_114_17 * r_2 * h7_p4 + e_3 * fs_57_221_26 * r_2 * h7_p6 + e_4 * fs_2_221_195 * r_2 * h9_p4 + e_4 * fs_3_221_910 * r_2 * h9_p6 - e_4 * f_6_17 * r_4 * h7_p4 - e_4 * fs_3_221_26 * r_4 * h7_p6;

        pc_10[k] = - e_1 * fs_255_4_11 * h5_m5 + e_2 * fs_399_26_6 * h7_m5 + e_2 * fs_765_26_11 * r_2 * h5_m5 + e_3 * fs_84_221_91 * h9_m5 - e_3 * fs_798_221_6 * r_2 * h7_m5 - e_3 * fs_51_13_11 * r_4 * h5_m5 - e_4 * fs_8_221_91 * r_2 * h9_m5 + e_4 * fs_42_221_6 * r_4 * h7_m5 + e_4 * fs_2_13_11 * r_6 * h5_m5;

        pc_11[k] = e_2 * fs_57_52_26 * h7_m6 - e_2 * f_57_2 * h7_m4 + e_3 * fs_63_442_910 * h9_m6 - e_3 * fs_21_221_195 * h9_m4 - e_3 * fs_57_221_26 * r_2 * h7_m6 + e_3 * f_114_17 * r_2 * h7_m4 - e_4 * fs_3_221_910 * r_2 * h9_m6 + e_4 * fs_2_221_195 * r_2 * h9_m4 + e_4 * fs_3_221_26 * r_4 * h7_m6 - e_4 * f_6_17 * r_4 * h7_m4;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m7, ph7_m3, ph7_m2, ph9_m8, ph9_m7, ph9_m3, ph9_m2, ab_2, pc_12, pc_13 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_12[k] = - e_0 * fs_225_8_231 * h3_m3 + e_1 * fs_85_4_33 * h5_m3 + e_1 * fs_25_1_231 * r_2 * h3_m3 - e_2 * fs_57_52_910 * h7_m7 + e_2 * fs_1083_572_110 * h7_m3 - e_2 * fs_255_26_33 * r_2 * h5_m3 - e_2 * fs_75_11_231 * r_4 * h3_m3 + e_3 * fs_84_221_273 * h9_m7 + e_3 * f_126_221 * h9_m3 + e_3 * fs_57_221_910 * r_2 * h7_m7 - e_3 * fs_1083_2431_110 * r_2 * h7_m3 + e_3 * fs_17_13_33 * r_4 * h5_m3 + e_3 * fs_100_143_231 * r_6 * h3_m3 - e_4 * fs_8_221_273 * r_2 * h9_m7 - e_4 * f_12_221 * r_2 * h9_m3 - e_4 * fs_3_221_910 * r_4 * h7_m7 + e_4 * fs_57_2431_110 * r_4 * h7_m3 - e_4 * fs_2_39_33 * r_6 * h5_m3 - e_4 * fs_10_429_231 * r_8 * h3_m3;

        pc_13[k] = - e_0 * fs_225_8_231 * h3_m2 - e_1 * fs_85_4_33 * h5_m2 + e_1 * fs_25_1_231 * r_2 * h3_m2 - e_2 * fs_285_572_330 * h7_m2 + e_2 * fs_255_26_33 * r_2 * h5_m2 - e_2 * fs_75_11_231 * r_4 * h3_m2 + e_3 * fs_42_221_1547 * h9_m8 - e_3 * fs_21_442_14 * h9_m2 + e_3 * fs_285_2431_330 * r_2 * h7_m2 - e_3 * fs_17_13_33 * r_4 * h5_m2 + e_3 * fs_100_143_231 * r_6 * h3_m2 - e_4 * fs_4_221_1547 * r_2 * h9_m8 + e_4 * fs_1_221_14 * r_2 * h9_m2 - e_4 * fs_15_2431_330 * r_4 * h7_m2 + e_4 * fs_2_39_33 * r_6 * h5_m2 - e_4 * fs_10_429_231 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph7_p6, ph7_p7, ph9_p1, ph9_p2, ph9_p6, ph9_p7, ab_2, pc_14, pc_15 : simd::cache_line_size())
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

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_14[k] = - e_0 * fs_225_8_105 * h3_p1 - e_1 * fs_85_4_42 * h5_p1 + e_1 * fs_25_1_105 * r_2 * h3_p1 - e_2 * fs_2565_572_10 * h7_p1 - e_2 * fs_57_572_30030 * h7_p7 + e_2 * fs_255_26_42 * r_2 * h5_p1 - e_2 * fs_75_11_105 * r_4 * h3_p1 - e_3 * fs_21_221_14 * h9_p1 - e_3 * fs_42_221_1001 * h9_p7 + e_3 * fs_2565_2431_10 * r_2 * h7_p1 + e_3 * fs_57_2431_30030 * r_2 * h7_p7 - e_3 * fs_17_13_42 * r_4 * h5_p1 + e_3 * fs_100_143_105 * r_6 * h3_p1 + e_4 * fs_2_221_14 * r_2 * h9_p1 + e_4 * fs_4_221_1001 * r_2 * h9_p7 - e_4 * fs_135_2431_10 * r_4 * h7_p1 - e_4 * fs_3_2431_30030 * r_4 * h7_p7 + e_4 * fs_2_39_42 * r_6 * h5_p1 - e_4 * fs_10_429_105 * r_8 * h3_p1;

        pc_15[k] = - e_0 * fs_675_4_7 * h3_p2 + e_1 * f_255_4 * h5_p2 + e_1 * fs_150_1_7 * r_2 * h3_p2 + e_2 * fs_1140_143_10 * h7_p2 + e_2 * fs_114_143_1430 * h7_p6 - e_2 * f_765_26 * r_2 * h5_p2 - e_2 * fs_450_11_7 * r_4 * h3_p2 + e_3 * fs_21_442_462 * h9_p2 - e_3 * fs_63_442_2002 * h9_p6 - e_3 * fs_4560_2431_10 * r_2 * h7_p2 - e_3 * fs_456_2431_1430 * r_2 * h7_p6 + e_3 * f_51_13 * r_4 * h5_p2 + e_3 * fs_600_143_7 * r_6 * h3_p2 - e_4 * fs_1_221_462 * r_2 * h9_p2 + e_4 * fs_3_221_2002 * r_2 * h9_p6 + e_4 * fs_240_2431_10 * r_4 * h7_p2 + e_4 * fs_24_2431_1430 * r_4 * h7_p6 - e_4 * f_2_13 * r_6 * h5_p2 - e_4 * fs_20_143_7 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_m4, ph5_p3, ph5_p5, ph7_m4, ph7_p3, ph7_p5, ph9_m4, ph9_p3, ph9_p5, ab_2, pc_16, pc_17 : simd::cache_line_size())
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

        pc_16[k] = - e_0 * fs_225_8_105 * h3_p3 + e_1 * fs_85_4_15 * h5_p3 - e_1 * fs_255_4_3 * h5_p5 + e_1 * fs_25_1_105 * r_2 * h3_p3 - e_2 * fs_10887_572_2 * h7_p3 + e_2 * fs_1653_572_22 * h7_p5 - e_2 * fs_255_26_15 * r_2 * h5_p3 + e_2 * fs_765_26_3 * r_2 * h5_p5 - e_2 * fs_75_11_105 * r_4 * h3_p3 - e_3 * fs_63_221_55 * h9_p3 - e_3 * fs_21_221_3003 * h9_p5 + e_3 * fs_10887_2431_2 * r_2 * h7_p3 - e_3 * fs_1653_2431_22 * r_2 * h7_p5 + e_3 * fs_17_13_15 * r_4 * h5_p3 - e_3 * fs_51_13_3 * r_4 * h5_p5 + e_3 * fs_100_143_105 * r_6 * h3_p3 + e_4 * fs_6_221_55 * r_2 * h9_p3 + e_4 * fs_2_221_3003 * r_2 * h9_p5 - e_4 * fs_573_2431_2 * r_4 * h7_p3 + e_4 * fs_87_2431_22 * r_4 * h7_p5 - e_4 * fs_2_39_15 * r_6 * h5_p3 + e_4 * fs_2_13_3 * r_6 * h5_p5 - e_4 * fs_10_429_105 * r_8 * h3_p3;

        pc_17[k] = - e_1 * fs_255_4_5 * h5_m4 + e_2 * fs_456_143_33 * h7_m4 + e_2 * fs_765_26_5 * r_2 * h5_m4 + e_3 * fs_42_221_715 * h9_m4 - e_3 * fs_1824_2431_33 * r_2 * h7_m4 - e_3 * fs_51_13_5 * r_4 * h5_m4 - e_4 * fs_4_221_715 * r_2 * h9_m4 + e_4 * fs_96_2431_33 * r_4 * h7_m4 + e_4 * fs_2_13_5 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ab_2, pc_18 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];

        pc_18[k] = - e_0 * fs_225_8_105 * h3_m3 + e_1 * fs_255_4_3 * h5_m5 + e_1 * fs_85_4_15 * h5_m3 + e_1 * fs_25_1_105 * r_2 * h3_m3 - e_2 * fs_1653_572_22 * h7_m5 - e_2 * fs_10887_572_2 * h7_m3 - e_2 * fs_765_26_3 * r_2 * h5_m5 - e_2 * fs_255_26_15 * r_2 * h5_m3 - e_2 * fs_75_11_105 * r_4 * h3_m3 + e_3 * fs_21_221_3003 * h9_m5 - e_3 * fs_63_221_55 * h9_m3 + e_3 * fs_1653_2431_22 * r_2 * h7_m5 + e_3 * fs_10887_2431_2 * r_2 * h7_m3 + e_3 * fs_51_13_3 * r_4 * h5_m5 + e_3 * fs_17_13_15 * r_4 * h5_m3 + e_3 * fs_100_143_105 * r_6 * h3_m3 - e_4 * fs_2_221_3003 * r_2 * h9_m5 + e_4 * fs_6_221_55 * r_2 * h9_m3 - e_4 * fs_87_2431_22 * r_4 * h7_m5 - e_4 * fs_573_2431_2 * r_4 * h7_m3 - e_4 * fs_2_13_3 * r_6 * h5_m5 - e_4 * fs_2_39_15 * r_6 * h5_m3 - e_4 * fs_10_429_105 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m7, ph7_m6, ph7_m2, ph7_m1, ph9_m7, ph9_m6, ph9_m2, ph9_m1, ab_2, pc_19, pc_20 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_19[k] = - e_0 * fs_675_4_7 * h3_m2 + e_1 * f_255_4 * h5_m2 + e_1 * fs_150_1_7 * r_2 * h3_m2 - e_2 * fs_114_143_1430 * h7_m6 + e_2 * fs_1140_143_10 * h7_m2 - e_2 * f_765_26 * r_2 * h5_m2 - e_2 * fs_450_11_7 * r_4 * h3_m2 + e_3 * fs_63_442_2002 * h9_m6 + e_3 * fs_21_442_462 * h9_m2 + e_3 * fs_456_2431_1430 * r_2 * h7_m6 - e_3 * fs_4560_2431_10 * r_2 * h7_m2 + e_3 * f_51_13 * r_4 * h5_m2 + e_3 * fs_600_143_7 * r_6 * h3_m2 - e_4 * fs_3_221_2002 * r_2 * h9_m6 - e_4 * fs_1_221_462 * r_2 * h9_m2 - e_4 * fs_24_2431_1430 * r_4 * h7_m6 + e_4 * fs_240_2431_10 * r_4 * h7_m2 - e_4 * f_2_13 * r_6 * h5_m2 - e_4 * fs_20_143_7 * r_8 * h3_m2;

        pc_20[k] = - e_0 * fs_225_8_105 * h3_m1 - e_1 * fs_85_4_42 * h5_m1 + e_1 * fs_25_1_105 * r_2 * h3_m1 + e_2 * fs_57_572_30030 * h7_m7 - e_2 * fs_2565_572_10 * h7_m1 + e_2 * fs_255_26_42 * r_2 * h5_m1 - e_2 * fs_75_11_105 * r_4 * h3_m1 + e_3 * fs_42_221_1001 * h9_m7 - e_3 * fs_21_221_14 * h9_m1 - e_3 * fs_57_2431_30030 * r_2 * h7_m7 + e_3 * fs_2565_2431_10 * r_2 * h7_m1 - e_3 * fs_17_13_42 * r_4 * h5_m1 + e_3 * fs_100_143_105 * r_6 * h3_m1 - e_4 * fs_4_221_1001 * r_2 * h9_m7 + e_4 * fs_2_221_14 * r_2 * h9_m1 + e_4 * fs_3_2431_30030 * r_4 * h7_m7 - e_4 * fs_135_2431_10 * r_4 * h7_m1 + e_4 * fs_2_39_42 * r_6 * h5_m1 - e_4 * fs_10_429_105 * r_8 * h3_m1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_0, ph3_p1, ph5_0, ph5_p5, ph7_0, ph7_p1, ph7_p5, ph7_p6, ph9_0, ph9_p1, ph9_p5, ph9_p6, ab_2, pc_21, pc_22 : simd::cache_line_size())
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

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_21[k] = - e_0 * fs_225_4_21 * h3_0 - e_1 * fs_85_2_21 * h5_0 + e_1 * fs_50_1_21 * r_2 * h3_0 - e_2 * fs_855_143_21 * h7_0 - e_2 * fs_855_572_286 * h7_p6 + e_2 * fs_255_13_21 * r_2 * h5_0 - e_2 * fs_150_11_21 * r_4 * h3_0 - e_3 * fs_42_221_21 * h9_0 - e_3 * fs_21_442_10010 * h9_p6 + e_3 * fs_3420_2431_21 * r_2 * h7_0 + e_3 * fs_855_2431_286 * r_2 * h7_p6 - e_3 * fs_34_13_21 * r_4 * h5_0 + e_3 * fs_200_143_21 * r_6 * h3_0 + e_4 * fs_4_221_21 * r_2 * h9_0 + e_4 * fs_1_221_10010 * r_2 * h9_p6 - e_4 * fs_180_2431_21 * r_4 * h7_0 - e_4 * fs_45_2431_286 * r_4 * h7_p6 + e_4 * fs_4_39_21 * r_6 * h5_0 - e_4 * fs_20_429_21 * r_8 * h3_0;

        pc_22[k] = - e_0 * fs_675_8_21 * h3_p1 + e_1 * f_255_4 * h5_p5 + e_1 * fs_75_1_21 * r_2 * h3_p1 + e_2 * fs_855_44_2 * h7_p1 + e_2 * fs_1425_572_66 * h7_p5 - e_2 * f_765_26 * r_2 * h5_p5 - e_2 * fs_225_11_21 * r_4 * h3_p1 + e_3 * fs_42_221_70 * h9_p1 - e_3 * fs_42_221_1001 * h9_p5 - e_3 * fs_855_187_2 * r_2 * h7_p1 - e_3 * fs_1425_2431_66 * r_2 * h7_p5 + e_3 * f_51_13 * r_4 * h5_p5 + e_3 * fs_300_143_21 * r_6 * h3_p1 - e_4 * fs_4_221_70 * r_2 * h9_p1 + e_4 * fs_4_221_1001 * r_2 * h9_p5 + e_4 * fs_45_187_2 * r_4 * h7_p1 + e_4 * fs_75_2431_66 * r_4 * h7_p5 - e_4 * f_2_13 * r_6 * h5_p5 - e_4 * fs_10_143_21 * r_8 * h3_p1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_p2, ph5_m3, ph5_p2, ph5_p4, ph7_m3, ph7_p2, ph7_p4, ph9_m3, ph9_p2, ph9_p4, ab_2, pc_23, pc_24 : simd::cache_line_size())
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

        pc_23[k] = - e_0 * fs_675_8_21 * h3_p2 + e_1 * fs_255_4_3 * h5_p2 - e_1 * f_255_2 * h5_p4 + e_1 * fs_75_1_21 * r_2 * h3_p2 - e_2 * fs_2109_572_30 * h7_p2 + e_2 * fs_513_286_165 * h7_p4 - e_2 * fs_765_26_3 * r_2 * h5_p2 + e_2 * f_765_13 * r_2 * h5_p4 - e_2 * fs_225_11_21 * r_4 * h3_p2 - e_3 * fs_105_442_154 * h9_p2 - e_3 * fs_105_221_143 * h9_p4 + e_3 * fs_2109_2431_30 * r_2 * h7_p2 - e_3 * fs_1026_2431_165 * r_2 * h7_p4 + e_3 * fs_51_13_3 * r_4 * h5_p2 - e_3 * f_102_13 * r_4 * h5_p4 + e_3 * fs_300_143_21 * r_6 * h3_p2 + e_4 * fs_5_221_154 * r_2 * h9_p2 + e_4 * fs_10_221_143 * r_2 * h9_p4 - e_4 * fs_111_2431_30 * r_4 * h7_p2 + e_4 * fs_54_2431_165 * r_4 * h7_p4 - e_4 * fs_2_13_3 * r_6 * h5_p2 + e_4 * f_4_13 * r_6 * h5_p4 - e_4 * fs_10_143_21 * r_8 * h3_p2;

        pc_24[k] = - e_0 * fs_225_4_21 * h3_m3 - e_1 * fs_85_4_3 * h5_m3 + e_1 * fs_50_1_21 * r_2 * h3_m3 - e_2 * fs_171_286_10 * h7_m3 + e_2 * fs_255_26_3 * r_2 * h5_m3 - e_2 * fs_150_11_21 * r_4 * h3_m3 + e_3 * fs_420_221_11 * h9_m3 + e_3 * fs_342_2431_10 * r_2 * h7_m3 - e_3 * fs_17_13_3 * r_4 * h5_m3 + e_3 * fs_200_143_21 * r_6 * h3_m3 - e_4 * fs_40_221_11 * r_2 * h9_m3 - e_4 * fs_18_2431_10 * r_4 * h7_m3 + e_4 * fs_2_39_3 * r_6 * h5_m3 - e_4 * fs_20_429_21 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_25 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];

        pc_25[k] = - e_0 * fs_675_8_21 * h3_m2 + e_1 * f_255_2 * h5_m4 + e_1 * fs_255_4_3 * h5_m2 + e_1 * fs_75_1_21 * r_2 * h3_m2 - e_2 * fs_513_286_165 * h7_m4 - e_2 * fs_2109_572_30 * h7_m2 - e_2 * f_765_13 * r_2 * h5_m4 - e_2 * fs_765_26_3 * r_2 * h5_m2 - e_2 * fs_225_11_21 * r_4 * h3_m2 + e_3 * fs_105_221_143 * h9_m4 - e_3 * fs_105_442_154 * h9_m2 + e_3 * fs_1026_2431_165 * r_2 * h7_m4 + e_3 * fs_2109_2431_30 * r_2 * h7_m2 + e_3 * f_102_13 * r_4 * h5_m4 + e_3 * fs_51_13_3 * r_4 * h5_m2 + e_3 * fs_300_143_21 * r_6 * h3_m2 - e_4 * fs_10_221_143 * r_2 * h9_m4 + e_4 * fs_5_221_154 * r_2 * h9_m2 - e_4 * fs_54_2431_165 * r_4 * h7_m4 - e_4 * fs_111_2431_30 * r_4 * h7_m2 - e_4 * f_4_13 * r_6 * h5_m4 - e_4 * fs_2_13_3 * r_6 * h5_m2 - e_4 * fs_10_143_21 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m1, ph5_m5, ph7_m6, ph7_m5, ph7_m1, ph9_m6, ph9_m5, ph9_m1, ab_2, pc_26, pc_27 : simd::cache_line_size())
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

        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_26[k] = - e_0 * fs_675_8_21 * h3_m1 - e_1 * f_255_4 * h5_m5 + e_1 * fs_75_1_21 * r_2 * h3_m1 - e_2 * fs_1425_572_66 * h7_m5 + e_2 * fs_855_44_2 * h7_m1 + e_2 * f_765_26 * r_2 * h5_m5 - e_2 * fs_225_11_21 * r_4 * h3_m1 + e_3 * fs_42_221_1001 * h9_m5 + e_3 * fs_42_221_70 * h9_m1 + e_3 * fs_1425_2431_66 * r_2 * h7_m5 - e_3 * fs_855_187_2 * r_2 * h7_m1 - e_3 * f_51_13 * r_4 * h5_m5 + e_3 * fs_300_143_21 * r_6 * h3_m1 - e_4 * fs_4_221_1001 * r_2 * h9_m5 - e_4 * fs_4_221_70 * r_2 * h9_m1 - e_4 * fs_75_2431_66 * r_4 * h7_m5 + e_4 * fs_45_187_2 * r_4 * h7_m1 + e_4 * f_2_13 * r_6 * h5_m5 - e_4 * fs_10_143_21 * r_8 * h3_m1;

        pc_27[k] = e_2 * fs_855_572_286 * h7_m6 + e_3 * fs_21_442_10010 * h9_m6 - e_3 * fs_855_2431_286 * r_2 * h7_m6 - e_4 * fs_1_221_10010 * r_2 * h9_m6 + e_4 * fs_45_2431_286 * r_4 * h7_m6;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_28 : simd::cache_line_size())
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

        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_28[k] = e_0 * fs_225_8_14 * h3_p1 + e_1 * fs_85_4_35 * h5_p1 - e_1 * fs_85_8_6 * h5_p5 - e_1 * fs_25_1_14 * r_2 * h3_p1 + e_2 * fs_1995_143_3 * h7_p1 - e_2 * fs_2565_286_11 * h7_p5 - e_2 * fs_255_26_35 * r_2 * h5_p1 + e_2 * fs_255_52_6 * r_2 * h5_p5 + e_2 * fs_75_11_14 * r_4 * h3_p1 + e_3 * fs_21_221_105 * h9_p1 - e_3 * fs_21_442_6006 * h9_p5 - e_3 * fs_7980_2431_3 * r_2 * h7_p1 + e_3 * fs_5130_2431_11 * r_2 * h7_p5 + e_3 * fs_17_13_35 * r_4 * h5_p1 - e_3 * fs_17_26_6 * r_4 * h5_p5 - e_3 * fs_100_143_14 * r_6 * h3_p1 - e_4 * fs_2_221_105 * r_2 * h9_p1 + e_4 * fs_1_221_6006 * r_2 * h9_p5 + e_4 * fs_420_2431_3 * r_4 * h7_p1 - e_4 * fs_270_2431_11 * r_4 * h7_p5 - e_4 * fs_2_39_35 * r_6 * h5_p1 + e_4 * fs_1_39_6 * r_6 * h5_p5 + e_4 * fs_10_429_14 * r_8 * h3_p1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2, pc_29 : simd::cache_line_size())
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

        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_29[k] = - e_0 * fs_225_2_14 * h3_0 - e_1 * fs_85_4_14 * h5_0 + e_1 * fs_255_8_10 * h5_p4 + e_1 * fs_100_1_14 * r_2 * h3_0 + e_2 * fs_1425_143_14 * h7_0 + e_2 * fs_285_286_66 * h7_p4 + e_2 * fs_255_26_14 * r_2 * h5_0 - e_2 * fs_765_52_10 * r_2 * h5_p4 - e_2 * fs_300_11_14 * r_4 * h3_0 + e_3 * fs_189_221_14 * h9_0 - e_3 * fs_63_442_1430 * h9_p4 - e_3 * fs_5700_2431_14 * r_2 * h7_0 - e_3 * fs_570_2431_66 * r_2 * h7_p4 - e_3 * fs_17_13_14 * r_4 * h5_0 + e_3 * fs_51_26_10 * r_4 * h5_p4 + e_3 * fs_400_143_14 * r_6 * h3_0 - e_4 * fs_18_221_14 * r_2 * h9_0 + e_4 * fs_3_221_1430 * r_2 * h9_p4 + e_4 * fs_300_2431_14 * r_4 * h7_0 + e_4 * fs_30_2431_66 * r_4 * h7_p4 + e_4 * fs_2_39_14 * r_6 * h5_0 - e_4 * fs_1_13_10 * r_6 * h5_p4 - e_4 * fs_40_429_14 * r_8 * h3_0;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph3_p1, ph3_p3, ph5_m2, ph5_p1, ph5_p3, ph7_m2, ph7_p1, ph7_p3, ph9_m2, ph9_p1, ph9_p3, ab_2, pc_30, pc_31 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_30[k] = - e_0 * fs_225_8_210 * h3_p1 - e_0 * fs_225_8_14 * h3_p3 + e_1 * fs_85_4_21 * h5_p1 - e_1 * fs_595_8_2 * h5_p3 + e_1 * fs_25_1_210 * r_2 * h3_p1 + e_1 * fs_25_1_14 * r_2 * h3_p3 - e_2 * fs_1311_286_5 * h7_p1 + e_2 * fs_912_143_15 * h7_p3 - e_2 * fs_255_26_21 * r_2 * h5_p1 + e_2 * fs_1785_52_2 * r_2 * h5_p3 - e_2 * fs_75_11_210 * r_4 * h3_p1 - e_2 * fs_75_11_14 * r_4 * h3_p3 - e_3 * fs_315_221_7 * h9_p1 - e_3 * fs_315_442_66 * h9_p3 + e_3 * fs_2622_2431_5 * r_2 * h7_p1 - e_3 * fs_3648_2431_15 * r_2 * h7_p3 + e_3 * fs_17_13_21 * r_4 * h5_p1 - e_3 * fs_119_26_2 * r_4 * h5_p3 + e_3 * fs_100_143_210 * r_6 * h3_p1 + e_3 * fs_100_143_14 * r_6 * h3_p3 + e_4 * fs_30_221_7 * r_2 * h9_p1 + e_4 * fs_15_221_66 * r_2 * h9_p3 - e_4 * fs_138_2431_5 * r_4 * h7_p1 + e_4 * fs_192_2431_15 * r_4 * h7_p3 - e_4 * fs_2_39_21 * r_6 * h5_p1 + e_4 * fs_7_39_2 * r_6 * h5_p3 - e_4 * fs_10_429_210 * r_8 * h3_p1 - e_4 * fs_10_429_14 * r_8 * h3_p3;

        pc_31[k] = - e_0 * fs_225_2_14 * h3_m2 + e_1 * fs_85_2_2 * h5_m2 + e_1 * fs_100_1_14 * r_2 * h3_m2 - e_2 * fs_1197_143_5 * h7_m2 - e_2 * fs_255_13_2 * r_2 * h5_m2 - e_2 * fs_300_11_14 * r_4 * h3_m2 + e_3 * fs_105_221_231 * h9_m2 + e_3 * fs_4788_2431_5 * r_2 * h7_m2 + e_3 * fs_34_13_2 * r_4 * h5_m2 + e_3 * fs_400_143_14 * r_6 * h3_m2 - e_4 * fs_10_221_231 * r_2 * h9_m2 - e_4 * fs_252_2431_5 * r_4 * h7_m2 - e_4 * fs_4_39_2 * r_6 * h5_m2 - e_4 * fs_40_429_14 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m1, ph5_m4, ph5_m3, ph5_m1, ph7_m4, ph7_m3, ph7_m1, ph9_m4, ph9_m3, ph9_m1, ab_2, pc_32, pc_33 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_32[k] = e_0 * fs_225_8_14 * h3_m3 - e_0 * fs_225_8_210 * h3_m1 + e_1 * fs_595_8_2 * h5_m3 + e_1 * fs_85_4_21 * h5_m1 - e_1 * fs_25_1_14 * r_2 * h3_m3 + e_1 * fs_25_1_210 * r_2 * h3_m1 - e_2 * fs_912_143_15 * h7_m3 - e_2 * fs_1311_286_5 * h7_m1 - e_2 * fs_1785_52_2 * r_2 * h5_m3 - e_2 * fs_255_26_21 * r_2 * h5_m1 + e_2 * fs_75_11_14 * r_4 * h3_m3 - e_2 * fs_75_11_210 * r_4 * h3_m1 + e_3 * fs_315_442_66 * h9_m3 - e_3 * fs_315_221_7 * h9_m1 + e_3 * fs_3648_2431_15 * r_2 * h7_m3 + e_3 * fs_2622_2431_5 * r_2 * h7_m1 + e_3 * fs_119_26_2 * r_4 * h5_m3 + e_3 * fs_17_13_21 * r_4 * h5_m1 - e_3 * fs_100_143_14 * r_6 * h3_m3 + e_3 * fs_100_143_210 * r_6 * h3_m1 - e_4 * fs_15_221_66 * r_2 * h9_m3 + e_4 * fs_30_221_7 * r_2 * h9_m1 - e_4 * fs_192_2431_15 * r_4 * h7_m3 - e_4 * fs_138_2431_5 * r_4 * h7_m1 - e_4 * fs_7_39_2 * r_6 * h5_m3 - e_4 * fs_2_39_21 * r_6 * h5_m1 + e_4 * fs_10_429_14 * r_8 * h3_m3 - e_4 * fs_10_429_210 * r_8 * h3_m1;

        pc_33[k] = - e_1 * fs_255_8_10 * h5_m4 - e_2 * fs_285_286_66 * h7_m4 + e_2 * fs_765_52_10 * r_2 * h5_m4 + e_3 * fs_63_442_1430 * h9_m4 + e_3 * fs_570_2431_66 * r_2 * h7_m4 - e_3 * fs_51_26_10 * r_4 * h5_m4 - e_4 * fs_3_221_1430 * r_2 * h9_m4 - e_4 * fs_30_2431_66 * r_4 * h7_m4 + e_4 * fs_1_13_10 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ab_2, pc_34 : simd::cache_line_size())
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

        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_34[k] = - e_0 * fs_225_8_14 * h3_m1 + e_1 * fs_85_8_6 * h5_m5 - e_1 * fs_85_4_35 * h5_m1 + e_1 * fs_25_1_14 * r_2 * h3_m1 + e_2 * fs_2565_286_11 * h7_m5 - e_2 * fs_1995_143_3 * h7_m1 - e_2 * fs_255_52_6 * r_2 * h5_m5 + e_2 * fs_255_26_35 * r_2 * h5_m1 - e_2 * fs_75_11_14 * r_4 * h3_m1 + e_3 * fs_21_442_6006 * h9_m5 - e_3 * fs_21_221_105 * h9_m1 - e_3 * fs_5130_2431_11 * r_2 * h7_m5 + e_3 * fs_7980_2431_3 * r_2 * h7_m1 + e_3 * fs_17_26_6 * r_4 * h5_m5 - e_3 * fs_17_13_35 * r_4 * h5_m1 + e_3 * fs_100_143_14 * r_6 * h3_m1 - e_4 * fs_1_221_6006 * r_2 * h9_m5 + e_4 * fs_2_221_105 * r_2 * h9_m1 + e_4 * fs_270_2431_11 * r_4 * h7_m5 - e_4 * fs_420_2431_3 * r_4 * h7_m1 - e_4 * fs_1_39_6 * r_6 * h5_m5 + e_4 * fs_2_39_35 * r_6 * h5_m1 - e_4 * fs_10_429_14 * r_8 * h3_m1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_35 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_35[k] = - e_0 * fs_225_16_14 * h3_p2 - e_1 * fs_595_8_2 * h5_p2 - e_1 * fs_85_4_6 * h5_p4 + e_1 * fs_25_2_14 * r_2 * h3_p2 - e_2 * fs_3591_286_5 * h7_p2 - e_2 * fs_855_286_110 * h7_p4 + e_2 * fs_1785_52_2 * r_2 * h5_p2 + e_2 * fs_255_26_6 * r_2 * h5_p4 - e_2 * fs_75_22_14 * r_4 * h3_p2 - e_3 * fs_21_221_231 * h9_p2 - e_3 * fs_21_221_858 * h9_p4 + e_3 * fs_7182_2431_5 * r_2 * h7_p2 + e_3 * fs_1710_2431_110 * r_2 * h7_p4 - e_3 * fs_119_26_2 * r_4 * h5_p2 - e_3 * fs_17_13_6 * r_4 * h5_p4 + e_3 * fs_50_143_14 * r_6 * h3_p2 + e_4 * fs_2_221_231 * r_2 * h9_p2 + e_4 * fs_2_221_858 * r_2 * h9_p4 - e_4 * fs_378_2431_5 * r_4 * h7_p2 - e_4 * fs_90_2431_110 * r_4 * h7_p4 + e_4 * fs_7_39_2 * r_6 * h5_p2 + e_4 * fs_2_39_6 * r_6 * h5_p4 - e_4 * fs_5_429_14 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_36 : simd::cache_line_size())
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

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_36[k] = e_0 * fs_225_16_210 * h3_p1 + e_0 * fs_225_16_14 * h3_p3 + e_1 * fs_85_4_21 * h5_p1 + e_1 * fs_85_1_2 * h5_p3 - e_1 * fs_25_2_210 * r_2 * h3_p1 - e_1 * fs_25_2_14 * r_2 * h3_p3 - e_2 * fs_2793_286_5 * h7_p1 - e_2 * fs_285_286_15 * h7_p3 - e_2 * fs_255_26_21 * r_2 * h5_p1 - e_2 * fs_510_13_2 * r_2 * h5_p3 + e_2 * fs_75_22_210 * r_4 * h3_p1 + e_2 * fs_75_22_14 * r_4 * h3_p3 - e_3 * fs_252_221_7 * h9_p1 - e_3 * fs_126_221_66 * h9_p3 + e_3 * fs_5586_2431_5 * r_2 * h7_p1 + e_3 * fs_570_2431_15 * r_2 * h7_p3 + e_3 * fs_17_13_21 * r_4 * h5_p1 + e_3 * fs_68_13_2 * r_4 * h5_p3 - e_3 * fs_50_143_210 * r_6 * h3_p1 - e_3 * fs_50_143_14 * r_6 * h3_p3 + e_4 * fs_24_221_7 * r_2 * h9_p1 + e_4 * fs_12_221_66 * r_2 * h9_p3 - e_4 * fs_294_2431_5 * r_4 * h7_p1 - e_4 * fs_30_2431_15 * r_4 * h7_p3 - e_4 * fs_2_39_21 * r_6 * h5_p1 - e_4 * fs_8_39_2 * r_6 * h5_p3 + e_4 * fs_5_429_210 * r_8 * h3_p1 + e_4 * fs_5_429_14 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m1, ph3_0, ph3_p2, ph5_m1, ph5_0, ph5_p2, ph7_m1, ph7_0, ph7_p2, ph9_m1, ph9_0, ph9_p2, ab_2, pc_37, pc_38 : simd::cache_line_size())
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

        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];

        pc_37[k] = - e_0 * fs_1125_8_14 * h3_0 - e_0 * fs_225_16_210 * h3_p2 + e_1 * fs_85_4_14 * h5_0 - e_1 * fs_85_8_30 * h5_p2 + e_1 * fs_125_1_14 * r_2 * h3_0 + e_1 * fs_25_2_210 * r_2 * h3_p2 + e_2 * fs_57_143_14 * h7_0 + e_2 * fs_3363_286_3 * h7_p2 - e_2 * fs_255_26_14 * r_2 * h5_0 + e_2 * fs_255_52_30 * r_2 * h5_p2 - e_2 * fs_375_11_14 * r_4 * h3_0 - e_2 * fs_75_22_210 * r_4 * h3_p2 - e_3 * fs_378_221_14 * h9_0 - e_3 * fs_63_221_385 * h9_p2 - e_3 * fs_228_2431_14 * r_2 * h7_0 - e_3 * fs_6726_2431_3 * r_2 * h7_p2 + e_3 * fs_17_13_14 * r_4 * h5_0 - e_3 * fs_17_26_30 * r_4 * h5_p2 + e_3 * fs_500_143_14 * r_6 * h3_0 + e_3 * fs_50_143_210 * r_6 * h3_p2 + e_4 * fs_36_221_14 * r_2 * h9_0 + e_4 * fs_6_221_385 * r_2 * h9_p2 + e_4 * fs_12_2431_14 * r_4 * h7_0 + e_4 * fs_354_2431_3 * r_4 * h7_p2 - e_4 * fs_2_39_14 * r_6 * h5_0 + e_4 * fs_1_39_30 * r_6 * h5_p2 - e_4 * fs_50_429_14 * r_8 * h3_0 - e_4 * fs_5_429_210 * r_8 * h3_p2;

        pc_38[k] = - e_0 * fs_1125_8_14 * h3_m1 + e_1 * fs_85_4_35 * h5_m1 + e_1 * fs_125_1_14 * r_2 * h3_m1 - e_2 * fs_2451_143_3 * h7_m1 - e_2 * fs_255_26_35 * r_2 * h5_m1 - e_2 * fs_375_11_14 * r_4 * h3_m1 + e_3 * fs_168_221_105 * h9_m1 + e_3 * fs_9804_2431_3 * r_2 * h7_m1 + e_3 * fs_17_13_35 * r_4 * h5_m1 + e_3 * fs_500_143_14 * r_6 * h3_m1 - e_4 * fs_16_221_105 * r_2 * h9_m1 - e_4 * fs_516_2431_3 * r_4 * h7_m1 - e_4 * fs_2_39_35 * r_6 * h5_m1 - e_4 * fs_50_429_14 * r_8 * h3_m1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph3_m1, ph5_m3, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ph9_m3, ph9_m2, ph9_m1, ab_2, pc_39, pc_40 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_39[k] = e_0 * fs_225_16_210 * h3_m2 + e_1 * fs_85_8_30 * h5_m2 - e_1 * fs_25_2_210 * r_2 * h3_m2 - e_2 * fs_3363_286_3 * h7_m2 - e_2 * fs_255_52_30 * r_2 * h5_m2 + e_2 * fs_75_22_210 * r_4 * h3_m2 + e_3 * fs_63_221_385 * h9_m2 + e_3 * fs_6726_2431_3 * r_2 * h7_m2 + e_3 * fs_17_26_30 * r_4 * h5_m2 - e_3 * fs_50_143_210 * r_6 * h3_m2 - e_4 * fs_6_221_385 * r_2 * h9_m2 - e_4 * fs_354_2431_3 * r_4 * h7_m2 - e_4 * fs_1_39_30 * r_6 * h5_m2 + e_4 * fs_5_429_210 * r_8 * h3_m2;

        pc_40[k] = - e_0 * fs_225_16_14 * h3_m3 - e_0 * fs_225_16_210 * h3_m1 - e_1 * fs_85_1_2 * h5_m3 - e_1 * fs_85_4_21 * h5_m1 + e_1 * fs_25_2_14 * r_2 * h3_m3 + e_1 * fs_25_2_210 * r_2 * h3_m1 + e_2 * fs_285_286_15 * h7_m3 + e_2 * fs_2793_286_5 * h7_m1 + e_2 * fs_510_13_2 * r_2 * h5_m3 + e_2 * fs_255_26_21 * r_2 * h5_m1 - e_2 * fs_75_22_14 * r_4 * h3_m3 - e_2 * fs_75_22_210 * r_4 * h3_m1 + e_3 * fs_126_221_66 * h9_m3 + e_3 * fs_252_221_7 * h9_m1 - e_3 * fs_570_2431_15 * r_2 * h7_m3 - e_3 * fs_5586_2431_5 * r_2 * h7_m1 - e_3 * fs_68_13_2 * r_4 * h5_m3 - e_3 * fs_17_13_21 * r_4 * h5_m1 + e_3 * fs_50_143_14 * r_6 * h3_m3 + e_3 * fs_50_143_210 * r_6 * h3_m1 - e_4 * fs_12_221_66 * r_2 * h9_m3 - e_4 * fs_24_221_7 * r_2 * h9_m1 + e_4 * fs_30_2431_15 * r_4 * h7_m3 + e_4 * fs_294_2431_5 * r_4 * h7_m1 + e_4 * fs_8_39_2 * r_6 * h5_m3 + e_4 * fs_2_39_21 * r_6 * h5_m1 - e_4 * fs_5_429_14 * r_8 * h3_m3 - e_4 * fs_5_429_210 * r_8 * h3_m1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph5_m4, ph5_m3, ph5_m2, ph7_m4, ph7_m3, ph7_m2, ph9_m4, ph9_m3, ph9_m2, ab_2, pc_41, pc_42, pc_43 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_41[k] = e_0 * fs_225_16_14 * h3_m2 + e_1 * fs_85_4_6 * h5_m4 + e_1 * fs_595_8_2 * h5_m2 - e_1 * fs_25_2_14 * r_2 * h3_m2 + e_2 * fs_855_286_110 * h7_m4 + e_2 * fs_3591_286_5 * h7_m2 - e_2 * fs_255_26_6 * r_2 * h5_m4 - e_2 * fs_1785_52_2 * r_2 * h5_m2 + e_2 * fs_75_22_14 * r_4 * h3_m2 + e_3 * fs_21_221_858 * h9_m4 + e_3 * fs_21_221_231 * h9_m2 - e_3 * fs_1710_2431_110 * r_2 * h7_m4 - e_3 * fs_7182_2431_5 * r_2 * h7_m2 + e_3 * fs_17_13_6 * r_4 * h5_m4 + e_3 * fs_119_26_2 * r_4 * h5_m2 - e_3 * fs_50_143_14 * r_6 * h3_m2 - e_4 * fs_2_221_858 * r_2 * h9_m4 - e_4 * fs_2_221_231 * r_2 * h9_m2 + e_4 * fs_90_2431_110 * r_4 * h7_m4 + e_4 * fs_378_2431_5 * r_4 * h7_m2 - e_4 * fs_2_39_6 * r_6 * h5_m4 - e_4 * fs_7_39_2 * r_6 * h5_m2 + e_4 * fs_5_429_14 * r_8 * h3_m2;

        pc_42[k] = e_0 * f_225_8 * h3_m3 + e_1 * fs_85_2_7 * h5_m3 - e_1 * f_25_1 * r_2 * h3_m3 + e_2 * fs_855_286_210 * h7_m3 - e_2 * fs_255_13_7 * r_2 * h5_m3 + e_2 * f_75_11 * r_4 * h3_m3 + e_3 * fs_42_221_231 * h9_m3 - e_3 * fs_1710_2431_210 * r_2 * h7_m3 + e_3 * fs_34_13_7 * r_4 * h5_m3 - e_3 * f_100_143 * r_6 * h3_m3 - e_4 * fs_4_221_231 * r_2 * h9_m3 + e_4 * fs_90_2431_210 * r_4 * h7_m3 - e_4 * fs_4_39_7 * r_6 * h5_m3 + e_4 * f_10_429 * r_8 * h3_m3;

        pc_43[k] = - e_0 * f_675_4 * h3_m2 - e_1 * fs_255_4_7 * h5_m2 + e_1 * f_150_1 * r_2 * h3_m2 + e_2 * fs_342_143_70 * h7_m2 + e_2 * fs_765_26_7 * r_2 * h5_m2 - e_2 * f_450_11 * r_4 * h3_m2 + e_3 * fs_147_221_66 * h9_m2 - e_3 * fs_1368_2431_70 * r_2 * h7_m2 - e_3 * fs_51_13_7 * r_4 * h5_m2 + e_3 * f_600_143 * r_6 * h3_m2 - e_4 * fs_14_221_66 * r_2 * h9_m2 + e_4 * fs_72_2431_70 * r_4 * h7_m2 + e_4 * fs_2_13_7 * r_6 * h5_m2 - e_4 * f_20_143 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m1, ph3_0, ph3_p1, ph5_0, ph7_m1, ph7_0, ph7_p1, ph9_m1, ph9_0, ph9_p1, ab_2, pc_44, pc_45, pc_46 : simd::cache_line_size())
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

        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];

        pc_44[k] = e_0 * f_3375_8 * h3_m1 - e_1 * f_375_1 * r_2 * h3_m1 - e_2 * fs_57_22_42 * h7_m1 + e_2 * f_1125_11 * r_4 * h3_m1 + e_3 * fs_294_221_30 * h9_m1 + e_3 * fs_114_187_42 * r_2 * h7_m1 - e_3 * f_1500_143 * r_6 * h3_m1 - e_4 * fs_28_221_30 * r_2 * h9_m1 - e_4 * fs_6_187_42 * r_4 * h7_m1 + e_4 * f_50_143 * r_8 * h3_m1;

        pc_45[k] = - e_0 * f_1125_2 * h3_0 + e_1 * f_595_4 * h5_0 + e_1 * f_500_1 * r_2 * h3_0 - e_2 * f_4788_143 * h7_0 - e_2 * f_1785_26 * r_2 * h5_0 - e_2 * f_1500_11 * r_4 * h3_0 + e_3 * f_1764_221 * h9_0 + e_3 * f_19152_2431 * r_2 * h7_0 + e_3 * f_119_13 * r_4 * h5_0 + e_3 * f_2000_143 * r_6 * h3_0 - e_4 * f_168_221 * r_2 * h9_0 - e_4 * f_1008_2431 * r_4 * h7_0 - e_4 * f_14_39 * r_6 * h5_0 - e_4 * f_200_429 * r_8 * h3_0;

        pc_46[k] = e_0 * f_3375_8 * h3_p1 - e_1 * f_375_1 * r_2 * h3_p1 - e_2 * fs_57_22_42 * h7_p1 + e_2 * f_1125_11 * r_4 * h3_p1 + e_3 * fs_294_221_30 * h9_p1 + e_3 * fs_114_187_42 * r_2 * h7_p1 - e_3 * f_1500_143 * r_6 * h3_p1 - e_4 * fs_28_221_30 * r_2 * h9_p1 - e_4 * fs_6_187_42 * r_4 * h7_p1 + e_4 * f_50_143 * r_8 * h3_p1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph9_p2, ph9_p3, ab_2, pc_47, pc_48 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];

        pc_47[k] = - e_0 * f_675_4 * h3_p2 - e_1 * fs_255_4_7 * h5_p2 + e_1 * f_150_1 * r_2 * h3_p2 + e_2 * fs_342_143_70 * h7_p2 + e_2 * fs_765_26_7 * r_2 * h5_p2 - e_2 * f_450_11 * r_4 * h3_p2 + e_3 * fs_147_221_66 * h9_p2 - e_3 * fs_1368_2431_70 * r_2 * h7_p2 - e_3 * fs_51_13_7 * r_4 * h5_p2 + e_3 * f_600_143 * r_6 * h3_p2 - e_4 * fs_14_221_66 * r_2 * h9_p2 + e_4 * fs_72_2431_70 * r_4 * h7_p2 + e_4 * fs_2_13_7 * r_6 * h5_p2 - e_4 * f_20_143 * r_8 * h3_p2;

        pc_48[k] = e_0 * f_225_8 * h3_p3 + e_1 * fs_85_2_7 * h5_p3 - e_1 * f_25_1 * r_2 * h3_p3 + e_2 * fs_855_286_210 * h7_p3 - e_2 * fs_255_13_7 * r_2 * h5_p3 + e_2 * f_75_11 * r_4 * h3_p3 + e_3 * fs_42_221_231 * h9_p3 - e_3 * fs_1710_2431_210 * r_2 * h7_p3 + e_3 * fs_34_13_7 * r_4 * h5_p3 - e_3 * f_100_143 * r_6 * h3_p3 - e_4 * fs_4_221_231 * r_2 * h9_p3 + e_4 * fs_90_2431_210 * r_4 * h7_p3 - e_4 * fs_4_39_7 * r_6 * h5_p3 + e_4 * f_10_429 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_49 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];

        pc_49[k] = - e_0 * fs_225_16_14 * h3_m2 + e_1 * fs_85_4_6 * h5_m4 - e_1 * fs_595_8_2 * h5_m2 + e_1 * fs_25_2_14 * r_2 * h3_m2 + e_2 * fs_855_286_110 * h7_m4 - e_2 * fs_3591_286_5 * h7_m2 - e_2 * fs_255_26_6 * r_2 * h5_m4 + e_2 * fs_1785_52_2 * r_2 * h5_m2 - e_2 * fs_75_22_14 * r_4 * h3_m2 + e_3 * fs_21_221_858 * h9_m4 - e_3 * fs_21_221_231 * h9_m2 - e_3 * fs_1710_2431_110 * r_2 * h7_m4 + e_3 * fs_7182_2431_5 * r_2 * h7_m2 + e_3 * fs_17_13_6 * r_4 * h5_m4 - e_3 * fs_119_26_2 * r_4 * h5_m2 + e_3 * fs_50_143_14 * r_6 * h3_m2 - e_4 * fs_2_221_858 * r_2 * h9_m4 + e_4 * fs_2_221_231 * r_2 * h9_m2 + e_4 * fs_90_2431_110 * r_4 * h7_m4 - e_4 * fs_378_2431_5 * r_4 * h7_m2 - e_4 * fs_2_39_6 * r_6 * h5_m4 + e_4 * fs_7_39_2 * r_6 * h5_m2 - e_4 * fs_5_429_14 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph3_m1, ph5_m3, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ph9_m3, ph9_m2, ph9_m1, ab_2, pc_50, pc_51 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_50[k] = - e_0 * fs_225_16_14 * h3_m3 + e_0 * fs_225_16_210 * h3_m1 - e_1 * fs_85_1_2 * h5_m3 + e_1 * fs_85_4_21 * h5_m1 + e_1 * fs_25_2_14 * r_2 * h3_m3 - e_1 * fs_25_2_210 * r_2 * h3_m1 + e_2 * fs_285_286_15 * h7_m3 - e_2 * fs_2793_286_5 * h7_m1 + e_2 * fs_510_13_2 * r_2 * h5_m3 - e_2 * fs_255_26_21 * r_2 * h5_m1 - e_2 * fs_75_22_14 * r_4 * h3_m3 + e_2 * fs_75_22_210 * r_4 * h3_m1 + e_3 * fs_126_221_66 * h9_m3 - e_3 * fs_252_221_7 * h9_m1 - e_3 * fs_570_2431_15 * r_2 * h7_m3 + e_3 * fs_5586_2431_5 * r_2 * h7_m1 - e_3 * fs_68_13_2 * r_4 * h5_m3 + e_3 * fs_17_13_21 * r_4 * h5_m1 + e_3 * fs_50_143_14 * r_6 * h3_m3 - e_3 * fs_50_143_210 * r_6 * h3_m1 - e_4 * fs_12_221_66 * r_2 * h9_m3 + e_4 * fs_24_221_7 * r_2 * h9_m1 + e_4 * fs_30_2431_15 * r_4 * h7_m3 - e_4 * fs_294_2431_5 * r_4 * h7_m1 + e_4 * fs_8_39_2 * r_6 * h5_m3 - e_4 * fs_2_39_21 * r_6 * h5_m1 - e_4 * fs_5_429_14 * r_8 * h3_m3 + e_4 * fs_5_429_210 * r_8 * h3_m1;

        pc_51[k] = e_0 * fs_225_16_210 * h3_m2 + e_1 * fs_85_8_30 * h5_m2 - e_1 * fs_25_2_210 * r_2 * h3_m2 - e_2 * fs_3363_286_3 * h7_m2 - e_2 * fs_255_52_30 * r_2 * h5_m2 + e_2 * fs_75_22_210 * r_4 * h3_m2 + e_3 * fs_63_221_385 * h9_m2 + e_3 * fs_6726_2431_3 * r_2 * h7_m2 + e_3 * fs_17_26_30 * r_4 * h5_m2 - e_3 * fs_50_143_210 * r_6 * h3_m2 - e_4 * fs_6_221_385 * r_2 * h9_m2 - e_4 * fs_354_2431_3 * r_4 * h7_m2 - e_4 * fs_1_39_30 * r_6 * h5_m2 + e_4 * fs_5_429_210 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_0, ph3_p1, ph3_p2, ph5_0, ph5_p1, ph5_p2, ph7_0, ph7_p1, ph7_p2, ph9_0, ph9_p1, ph9_p2, ab_2, pc_52, pc_53 : simd::cache_line_size())
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

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];

        pc_52[k] = - e_0 * fs_1125_8_14 * h3_p1 + e_1 * fs_85_4_35 * h5_p1 + e_1 * fs_125_1_14 * r_2 * h3_p1 - e_2 * fs_2451_143_3 * h7_p1 - e_2 * fs_255_26_35 * r_2 * h5_p1 - e_2 * fs_375_11_14 * r_4 * h3_p1 + e_3 * fs_168_221_105 * h9_p1 + e_3 * fs_9804_2431_3 * r_2 * h7_p1 + e_3 * fs_17_13_35 * r_4 * h5_p1 + e_3 * fs_500_143_14 * r_6 * h3_p1 - e_4 * fs_16_221_105 * r_2 * h9_p1 - e_4 * fs_516_2431_3 * r_4 * h7_p1 - e_4 * fs_2_39_35 * r_6 * h5_p1 - e_4 * fs_50_429_14 * r_8 * h3_p1;

        pc_53[k] = - e_0 * fs_1125_8_14 * h3_0 + e_0 * fs_225_16_210 * h3_p2 + e_1 * fs_85_4_14 * h5_0 + e_1 * fs_85_8_30 * h5_p2 + e_1 * fs_125_1_14 * r_2 * h3_0 - e_1 * fs_25_2_210 * r_2 * h3_p2 + e_2 * fs_57_143_14 * h7_0 - e_2 * fs_3363_286_3 * h7_p2 - e_2 * fs_255_26_14 * r_2 * h5_0 - e_2 * fs_255_52_30 * r_2 * h5_p2 - e_2 * fs_375_11_14 * r_4 * h3_0 + e_2 * fs_75_22_210 * r_4 * h3_p2 - e_3 * fs_378_221_14 * h9_0 + e_3 * fs_63_221_385 * h9_p2 - e_3 * fs_228_2431_14 * r_2 * h7_0 + e_3 * fs_6726_2431_3 * r_2 * h7_p2 + e_3 * fs_17_13_14 * r_4 * h5_0 + e_3 * fs_17_26_30 * r_4 * h5_p2 + e_3 * fs_500_143_14 * r_6 * h3_0 - e_3 * fs_50_143_210 * r_6 * h3_p2 + e_4 * fs_36_221_14 * r_2 * h9_0 - e_4 * fs_6_221_385 * r_2 * h9_p2 + e_4 * fs_12_2431_14 * r_4 * h7_0 - e_4 * fs_354_2431_3 * r_4 * h7_p2 - e_4 * fs_2_39_14 * r_6 * h5_0 - e_4 * fs_1_39_30 * r_6 * h5_p2 - e_4 * fs_50_429_14 * r_8 * h3_0 + e_4 * fs_5_429_210 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_54 : simd::cache_line_size())
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

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_54[k] = e_0 * fs_225_16_210 * h3_p1 - e_0 * fs_225_16_14 * h3_p3 + e_1 * fs_85_4_21 * h5_p1 - e_1 * fs_85_1_2 * h5_p3 - e_1 * fs_25_2_210 * r_2 * h3_p1 + e_1 * fs_25_2_14 * r_2 * h3_p3 - e_2 * fs_2793_286_5 * h7_p1 + e_2 * fs_285_286_15 * h7_p3 - e_2 * fs_255_26_21 * r_2 * h5_p1 + e_2 * fs_510_13_2 * r_2 * h5_p3 + e_2 * fs_75_22_210 * r_4 * h3_p1 - e_2 * fs_75_22_14 * r_4 * h3_p3 - e_3 * fs_252_221_7 * h9_p1 + e_3 * fs_126_221_66 * h9_p3 + e_3 * fs_5586_2431_5 * r_2 * h7_p1 - e_3 * fs_570_2431_15 * r_2 * h7_p3 + e_3 * fs_17_13_21 * r_4 * h5_p1 - e_3 * fs_68_13_2 * r_4 * h5_p3 - e_3 * fs_50_143_210 * r_6 * h3_p1 + e_3 * fs_50_143_14 * r_6 * h3_p3 + e_4 * fs_24_221_7 * r_2 * h9_p1 - e_4 * fs_12_221_66 * r_2 * h9_p3 - e_4 * fs_294_2431_5 * r_4 * h7_p1 + e_4 * fs_30_2431_15 * r_4 * h7_p3 - e_4 * fs_2_39_21 * r_6 * h5_p1 + e_4 * fs_8_39_2 * r_6 * h5_p3 + e_4 * fs_5_429_210 * r_8 * h3_p1 - e_4 * fs_5_429_14 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_55 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_55[k] = - e_0 * fs_225_16_14 * h3_p2 - e_1 * fs_595_8_2 * h5_p2 + e_1 * fs_85_4_6 * h5_p4 + e_1 * fs_25_2_14 * r_2 * h3_p2 - e_2 * fs_3591_286_5 * h7_p2 + e_2 * fs_855_286_110 * h7_p4 + e_2 * fs_1785_52_2 * r_2 * h5_p2 - e_2 * fs_255_26_6 * r_2 * h5_p4 - e_2 * fs_75_22_14 * r_4 * h3_p2 - e_3 * fs_21_221_231 * h9_p2 + e_3 * fs_21_221_858 * h9_p4 + e_3 * fs_7182_2431_5 * r_2 * h7_p2 - e_3 * fs_1710_2431_110 * r_2 * h7_p4 - e_3 * fs_119_26_2 * r_4 * h5_p2 + e_3 * fs_17_13_6 * r_4 * h5_p4 + e_3 * fs_50_143_14 * r_6 * h3_p2 + e_4 * fs_2_221_231 * r_2 * h9_p2 - e_4 * fs_2_221_858 * r_2 * h9_p4 - e_4 * fs_378_2431_5 * r_4 * h7_p2 + e_4 * fs_90_2431_110 * r_4 * h7_p4 + e_4 * fs_7_39_2 * r_6 * h5_p2 - e_4 * fs_2_39_6 * r_6 * h5_p4 - e_4 * fs_5_429_14 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m1, ph5_m5, ph5_m4, ph5_m1, ph7_m5, ph7_m4, ph7_m1, ph9_m5, ph9_m4, ph9_m1, ab_2, pc_56, pc_57 : simd::cache_line_size())
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

        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m1 = ph9_m1[k];

        pc_56[k] = e_0 * fs_225_8_14 * h3_m1 + e_1 * fs_85_8_6 * h5_m5 + e_1 * fs_85_4_35 * h5_m1 - e_1 * fs_25_1_14 * r_2 * h3_m1 + e_2 * fs_2565_286_11 * h7_m5 + e_2 * fs_1995_143_3 * h7_m1 - e_2 * fs_255_52_6 * r_2 * h5_m5 - e_2 * fs_255_26_35 * r_2 * h5_m1 + e_2 * fs_75_11_14 * r_4 * h3_m1 + e_3 * fs_21_442_6006 * h9_m5 + e_3 * fs_21_221_105 * h9_m1 - e_3 * fs_5130_2431_11 * r_2 * h7_m5 - e_3 * fs_7980_2431_3 * r_2 * h7_m1 + e_3 * fs_17_26_6 * r_4 * h5_m5 + e_3 * fs_17_13_35 * r_4 * h5_m1 - e_3 * fs_100_143_14 * r_6 * h3_m1 - e_4 * fs_1_221_6006 * r_2 * h9_m5 - e_4 * fs_2_221_105 * r_2 * h9_m1 + e_4 * fs_270_2431_11 * r_4 * h7_m5 + e_4 * fs_420_2431_3 * r_4 * h7_m1 - e_4 * fs_1_39_6 * r_6 * h5_m5 - e_4 * fs_2_39_35 * r_6 * h5_m1 + e_4 * fs_10_429_14 * r_8 * h3_m1;

        pc_57[k] = - e_1 * fs_255_8_10 * h5_m4 - e_2 * fs_285_286_66 * h7_m4 + e_2 * fs_765_52_10 * r_2 * h5_m4 + e_3 * fs_63_442_1430 * h9_m4 + e_3 * fs_570_2431_66 * r_2 * h7_m4 - e_3 * fs_51_26_10 * r_4 * h5_m4 - e_4 * fs_3_221_1430 * r_2 * h9_m4 - e_4 * fs_30_2431_66 * r_4 * h7_m4 + e_4 * fs_1_13_10 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m1, ph3_p2, ph5_m3, ph5_m1, ph5_p2, ph7_m3, ph7_m1, ph7_p2, ph9_m3, ph9_m1, ph9_p2, ab_2, pc_58, pc_59 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_p2 = ph9_p2[k];

        pc_58[k] = e_0 * fs_225_8_14 * h3_m3 + e_0 * fs_225_8_210 * h3_m1 + e_1 * fs_595_8_2 * h5_m3 - e_1 * fs_85_4_21 * h5_m1 - e_1 * fs_25_1_14 * r_2 * h3_m3 - e_1 * fs_25_1_210 * r_2 * h3_m1 - e_2 * fs_912_143_15 * h7_m3 + e_2 * fs_1311_286_5 * h7_m1 - e_2 * fs_1785_52_2 * r_2 * h5_m3 + e_2 * fs_255_26_21 * r_2 * h5_m1 + e_2 * fs_75_11_14 * r_4 * h3_m3 + e_2 * fs_75_11_210 * r_4 * h3_m1 + e_3 * fs_315_442_66 * h9_m3 + e_3 * fs_315_221_7 * h9_m1 + e_3 * fs_3648_2431_15 * r_2 * h7_m3 - e_3 * fs_2622_2431_5 * r_2 * h7_m1 + e_3 * fs_119_26_2 * r_4 * h5_m3 - e_3 * fs_17_13_21 * r_4 * h5_m1 - e_3 * fs_100_143_14 * r_6 * h3_m3 - e_3 * fs_100_143_210 * r_6 * h3_m1 - e_4 * fs_15_221_66 * r_2 * h9_m3 - e_4 * fs_30_221_7 * r_2 * h9_m1 - e_4 * fs_192_2431_15 * r_4 * h7_m3 + e_4 * fs_138_2431_5 * r_4 * h7_m1 - e_4 * fs_7_39_2 * r_6 * h5_m3 + e_4 * fs_2_39_21 * r_6 * h5_m1 + e_4 * fs_10_429_14 * r_8 * h3_m3 + e_4 * fs_10_429_210 * r_8 * h3_m1;

        pc_59[k] = - e_0 * fs_225_2_14 * h3_p2 + e_1 * fs_85_2_2 * h5_p2 + e_1 * fs_100_1_14 * r_2 * h3_p2 - e_2 * fs_1197_143_5 * h7_p2 - e_2 * fs_255_13_2 * r_2 * h5_p2 - e_2 * fs_300_11_14 * r_4 * h3_p2 + e_3 * fs_105_221_231 * h9_p2 + e_3 * fs_4788_2431_5 * r_2 * h7_p2 + e_3 * fs_34_13_2 * r_4 * h5_p2 + e_3 * fs_400_143_14 * r_6 * h3_p2 - e_4 * fs_10_221_231 * r_2 * h9_p2 - e_4 * fs_252_2431_5 * r_4 * h7_p2 - e_4 * fs_4_39_2 * r_6 * h5_p2 - e_4 * fs_40_429_14 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_60 : simd::cache_line_size())
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

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_60[k] = - e_0 * fs_225_8_210 * h3_p1 + e_0 * fs_225_8_14 * h3_p3 + e_1 * fs_85_4_21 * h5_p1 + e_1 * fs_595_8_2 * h5_p3 + e_1 * fs_25_1_210 * r_2 * h3_p1 - e_1 * fs_25_1_14 * r_2 * h3_p3 - e_2 * fs_1311_286_5 * h7_p1 - e_2 * fs_912_143_15 * h7_p3 - e_2 * fs_255_26_21 * r_2 * h5_p1 - e_2 * fs_1785_52_2 * r_2 * h5_p3 - e_2 * fs_75_11_210 * r_4 * h3_p1 + e_2 * fs_75_11_14 * r_4 * h3_p3 - e_3 * fs_315_221_7 * h9_p1 + e_3 * fs_315_442_66 * h9_p3 + e_3 * fs_2622_2431_5 * r_2 * h7_p1 + e_3 * fs_3648_2431_15 * r_2 * h7_p3 + e_3 * fs_17_13_21 * r_4 * h5_p1 + e_3 * fs_119_26_2 * r_4 * h5_p3 + e_3 * fs_100_143_210 * r_6 * h3_p1 - e_3 * fs_100_143_14 * r_6 * h3_p3 + e_4 * fs_30_221_7 * r_2 * h9_p1 - e_4 * fs_15_221_66 * r_2 * h9_p3 - e_4 * fs_138_2431_5 * r_4 * h7_p1 - e_4 * fs_192_2431_15 * r_4 * h7_p3 - e_4 * fs_2_39_21 * r_6 * h5_p1 - e_4 * fs_7_39_2 * r_6 * h5_p3 - e_4 * fs_10_429_210 * r_8 * h3_p1 + e_4 * fs_10_429_14 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2, pc_61 : simd::cache_line_size())
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

        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_61[k] = - e_0 * fs_225_2_14 * h3_0 - e_1 * fs_85_4_14 * h5_0 - e_1 * fs_255_8_10 * h5_p4 + e_1 * fs_100_1_14 * r_2 * h3_0 + e_2 * fs_1425_143_14 * h7_0 - e_2 * fs_285_286_66 * h7_p4 + e_2 * fs_255_26_14 * r_2 * h5_0 + e_2 * fs_765_52_10 * r_2 * h5_p4 - e_2 * fs_300_11_14 * r_4 * h3_0 + e_3 * fs_189_221_14 * h9_0 + e_3 * fs_63_442_1430 * h9_p4 - e_3 * fs_5700_2431_14 * r_2 * h7_0 + e_3 * fs_570_2431_66 * r_2 * h7_p4 - e_3 * fs_17_13_14 * r_4 * h5_0 - e_3 * fs_51_26_10 * r_4 * h5_p4 + e_3 * fs_400_143_14 * r_6 * h3_0 - e_4 * fs_18_221_14 * r_2 * h9_0 - e_4 * fs_3_221_1430 * r_2 * h9_p4 + e_4 * fs_300_2431_14 * r_4 * h7_0 - e_4 * fs_30_2431_66 * r_4 * h7_p4 + e_4 * fs_2_39_14 * r_6 * h5_0 + e_4 * fs_1_13_10 * r_6 * h5_p4 - e_4 * fs_40_429_14 * r_8 * h3_0;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p1, ph5_p1, ph5_p5, ph7_m6, ph7_p1, ph7_p5, ph9_m6, ph9_p1, ph9_p5, ab_2, pc_62, pc_63 : simd::cache_line_size())
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

        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_62[k] = e_0 * fs_225_8_14 * h3_p1 + e_1 * fs_85_4_35 * h5_p1 + e_1 * fs_85_8_6 * h5_p5 - e_1 * fs_25_1_14 * r_2 * h3_p1 + e_2 * fs_1995_143_3 * h7_p1 + e_2 * fs_2565_286_11 * h7_p5 - e_2 * fs_255_26_35 * r_2 * h5_p1 - e_2 * fs_255_52_6 * r_2 * h5_p5 + e_2 * fs_75_11_14 * r_4 * h3_p1 + e_3 * fs_21_221_105 * h9_p1 + e_3 * fs_21_442_6006 * h9_p5 - e_3 * fs_7980_2431_3 * r_2 * h7_p1 - e_3 * fs_5130_2431_11 * r_2 * h7_p5 + e_3 * fs_17_13_35 * r_4 * h5_p1 + e_3 * fs_17_26_6 * r_4 * h5_p5 - e_3 * fs_100_143_14 * r_6 * h3_p1 - e_4 * fs_2_221_105 * r_2 * h9_p1 - e_4 * fs_1_221_6006 * r_2 * h9_p5 + e_4 * fs_420_2431_3 * r_4 * h7_p1 + e_4 * fs_270_2431_11 * r_4 * h7_p5 - e_4 * fs_2_39_35 * r_6 * h5_p1 - e_4 * fs_1_39_6 * r_6 * h5_p5 + e_4 * fs_10_429_14 * r_8 * h3_p1;

        pc_63[k] = e_2 * fs_855_572_286 * h7_m6 + e_3 * fs_21_442_10010 * h9_m6 - e_3 * fs_855_2431_286 * r_2 * h7_m6 - e_4 * fs_1_221_10010 * r_2 * h9_m6 + e_4 * fs_45_2431_286 * r_4 * h7_m6;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m1, ph5_m5, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ab_2, pc_64 : simd::cache_line_size())
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

        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_64[k] = e_0 * fs_675_8_21 * h3_m1 - e_1 * f_255_4 * h5_m5 - e_1 * fs_75_1_21 * r_2 * h3_m1 - e_2 * fs_1425_572_66 * h7_m5 - e_2 * fs_855_44_2 * h7_m1 + e_2 * f_765_26 * r_2 * h5_m5 + e_2 * fs_225_11_21 * r_4 * h3_m1 + e_3 * fs_42_221_1001 * h9_m5 - e_3 * fs_42_221_70 * h9_m1 + e_3 * fs_1425_2431_66 * r_2 * h7_m5 + e_3 * fs_855_187_2 * r_2 * h7_m1 - e_3 * f_51_13 * r_4 * h5_m5 - e_3 * fs_300_143_21 * r_6 * h3_m1 - e_4 * fs_4_221_1001 * r_2 * h9_m5 + e_4 * fs_4_221_70 * r_2 * h9_m1 - e_4 * fs_75_2431_66 * r_4 * h7_m5 - e_4 * fs_45_187_2 * r_4 * h7_m1 + e_4 * f_2_13 * r_6 * h5_m5 + e_4 * fs_10_143_21 * r_8 * h3_m1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph3_p3, ph5_m4, ph5_m2, ph5_p3, ph7_m4, ph7_m2, ph7_p3, ph9_m4, ph9_m2, ph9_p3, ab_2, pc_65, pc_66 : simd::cache_line_size())
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

        pc_65[k] = e_0 * fs_675_8_21 * h3_m2 + e_1 * f_255_2 * h5_m4 - e_1 * fs_255_4_3 * h5_m2 - e_1 * fs_75_1_21 * r_2 * h3_m2 - e_2 * fs_513_286_165 * h7_m4 + e_2 * fs_2109_572_30 * h7_m2 - e_2 * f_765_13 * r_2 * h5_m4 + e_2 * fs_765_26_3 * r_2 * h5_m2 + e_2 * fs_225_11_21 * r_4 * h3_m2 + e_3 * fs_105_221_143 * h9_m4 + e_3 * fs_105_442_154 * h9_m2 + e_3 * fs_1026_2431_165 * r_2 * h7_m4 - e_3 * fs_2109_2431_30 * r_2 * h7_m2 + e_3 * f_102_13 * r_4 * h5_m4 - e_3 * fs_51_13_3 * r_4 * h5_m2 - e_3 * fs_300_143_21 * r_6 * h3_m2 - e_4 * fs_10_221_143 * r_2 * h9_m4 - e_4 * fs_5_221_154 * r_2 * h9_m2 - e_4 * fs_54_2431_165 * r_4 * h7_m4 + e_4 * fs_111_2431_30 * r_4 * h7_m2 - e_4 * f_4_13 * r_6 * h5_m4 + e_4 * fs_2_13_3 * r_6 * h5_m2 + e_4 * fs_10_143_21 * r_8 * h3_m2;

        pc_66[k] = - e_0 * fs_225_4_21 * h3_p3 - e_1 * fs_85_4_3 * h5_p3 + e_1 * fs_50_1_21 * r_2 * h3_p3 - e_2 * fs_171_286_10 * h7_p3 + e_2 * fs_255_26_3 * r_2 * h5_p3 - e_2 * fs_150_11_21 * r_4 * h3_p3 + e_3 * fs_420_221_11 * h9_p3 + e_3 * fs_342_2431_10 * r_2 * h7_p3 - e_3 * fs_17_13_3 * r_4 * h5_p3 + e_3 * fs_200_143_21 * r_6 * h3_p3 - e_4 * fs_40_221_11 * r_2 * h9_p3 - e_4 * fs_18_2431_10 * r_4 * h7_p3 + e_4 * fs_2_39_3 * r_6 * h5_p3 - e_4 * fs_20_429_21 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_67 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_67[k] = - e_0 * fs_675_8_21 * h3_p2 + e_1 * fs_255_4_3 * h5_p2 + e_1 * f_255_2 * h5_p4 + e_1 * fs_75_1_21 * r_2 * h3_p2 - e_2 * fs_2109_572_30 * h7_p2 - e_2 * fs_513_286_165 * h7_p4 - e_2 * fs_765_26_3 * r_2 * h5_p2 - e_2 * f_765_13 * r_2 * h5_p4 - e_2 * fs_225_11_21 * r_4 * h3_p2 - e_3 * fs_105_442_154 * h9_p2 + e_3 * fs_105_221_143 * h9_p4 + e_3 * fs_2109_2431_30 * r_2 * h7_p2 + e_3 * fs_1026_2431_165 * r_2 * h7_p4 + e_3 * fs_51_13_3 * r_4 * h5_p2 + e_3 * f_102_13 * r_4 * h5_p4 + e_3 * fs_300_143_21 * r_6 * h3_p2 + e_4 * fs_5_221_154 * r_2 * h9_p2 - e_4 * fs_10_221_143 * r_2 * h9_p4 - e_4 * fs_111_2431_30 * r_4 * h7_p2 - e_4 * fs_54_2431_165 * r_4 * h7_p4 - e_4 * fs_2_13_3 * r_6 * h5_p2 - e_4 * f_4_13 * r_6 * h5_p4 - e_4 * fs_10_143_21 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_0, ph3_p1, ph5_0, ph5_p5, ph7_0, ph7_p1, ph7_p5, ph7_p6, ph9_0, ph9_p1, ph9_p5, ph9_p6, ab_2, pc_68, pc_69 : simd::cache_line_size())
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

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_68[k] = - e_0 * fs_675_8_21 * h3_p1 - e_1 * f_255_4 * h5_p5 + e_1 * fs_75_1_21 * r_2 * h3_p1 + e_2 * fs_855_44_2 * h7_p1 - e_2 * fs_1425_572_66 * h7_p5 + e_2 * f_765_26 * r_2 * h5_p5 - e_2 * fs_225_11_21 * r_4 * h3_p1 + e_3 * fs_42_221_70 * h9_p1 + e_3 * fs_42_221_1001 * h9_p5 - e_3 * fs_855_187_2 * r_2 * h7_p1 + e_3 * fs_1425_2431_66 * r_2 * h7_p5 - e_3 * f_51_13 * r_4 * h5_p5 + e_3 * fs_300_143_21 * r_6 * h3_p1 - e_4 * fs_4_221_70 * r_2 * h9_p1 - e_4 * fs_4_221_1001 * r_2 * h9_p5 + e_4 * fs_45_187_2 * r_4 * h7_p1 - e_4 * fs_75_2431_66 * r_4 * h7_p5 + e_4 * f_2_13 * r_6 * h5_p5 - e_4 * fs_10_143_21 * r_8 * h3_p1;

        pc_69[k] = - e_0 * fs_225_4_21 * h3_0 - e_1 * fs_85_2_21 * h5_0 + e_1 * fs_50_1_21 * r_2 * h3_0 - e_2 * fs_855_143_21 * h7_0 + e_2 * fs_855_572_286 * h7_p6 + e_2 * fs_255_13_21 * r_2 * h5_0 - e_2 * fs_150_11_21 * r_4 * h3_0 - e_3 * fs_42_221_21 * h9_0 + e_3 * fs_21_442_10010 * h9_p6 + e_3 * fs_3420_2431_21 * r_2 * h7_0 - e_3 * fs_855_2431_286 * r_2 * h7_p6 - e_3 * fs_34_13_21 * r_4 * h5_0 + e_3 * fs_200_143_21 * r_6 * h3_0 + e_4 * fs_4_221_21 * r_2 * h9_0 - e_4 * fs_1_221_10010 * r_2 * h9_p6 - e_4 * fs_180_2431_21 * r_4 * h7_0 + e_4 * fs_45_2431_286 * r_4 * h7_p6 + e_4 * fs_4_39_21 * r_6 * h5_0 - e_4 * fs_20_429_21 * r_8 * h3_0;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m7, ph7_m6, ph7_m2, ph7_m1, ph9_m7, ph9_m6, ph9_m2, ph9_m1, ab_2, pc_70, pc_71 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_70[k] = e_0 * fs_225_8_105 * h3_m1 + e_1 * fs_85_4_42 * h5_m1 - e_1 * fs_25_1_105 * r_2 * h3_m1 + e_2 * fs_57_572_30030 * h7_m7 + e_2 * fs_2565_572_10 * h7_m1 - e_2 * fs_255_26_42 * r_2 * h5_m1 + e_2 * fs_75_11_105 * r_4 * h3_m1 + e_3 * fs_42_221_1001 * h9_m7 + e_3 * fs_21_221_14 * h9_m1 - e_3 * fs_57_2431_30030 * r_2 * h7_m7 - e_3 * fs_2565_2431_10 * r_2 * h7_m1 + e_3 * fs_17_13_42 * r_4 * h5_m1 - e_3 * fs_100_143_105 * r_6 * h3_m1 - e_4 * fs_4_221_1001 * r_2 * h9_m7 - e_4 * fs_2_221_14 * r_2 * h9_m1 + e_4 * fs_3_2431_30030 * r_4 * h7_m7 + e_4 * fs_135_2431_10 * r_4 * h7_m1 - e_4 * fs_2_39_42 * r_6 * h5_m1 + e_4 * fs_10_429_105 * r_8 * h3_m1;

        pc_71[k] = e_0 * fs_675_4_7 * h3_m2 - e_1 * f_255_4 * h5_m2 - e_1 * fs_150_1_7 * r_2 * h3_m2 - e_2 * fs_114_143_1430 * h7_m6 - e_2 * fs_1140_143_10 * h7_m2 + e_2 * f_765_26 * r_2 * h5_m2 + e_2 * fs_450_11_7 * r_4 * h3_m2 + e_3 * fs_63_442_2002 * h9_m6 - e_3 * fs_21_442_462 * h9_m2 + e_3 * fs_456_2431_1430 * r_2 * h7_m6 + e_3 * fs_4560_2431_10 * r_2 * h7_m2 - e_3 * f_51_13 * r_4 * h5_m2 - e_3 * fs_600_143_7 * r_6 * h3_m2 - e_4 * fs_3_221_2002 * r_2 * h9_m6 + e_4 * fs_1_221_462 * r_2 * h9_m2 - e_4 * fs_24_2431_1430 * r_4 * h7_m6 - e_4 * fs_240_2431_10 * r_4 * h7_m2 + e_4 * f_2_13 * r_6 * h5_m2 + e_4 * fs_20_143_7 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m3, ph5_p4, ph7_m5, ph7_m3, ph7_p4, ph9_m5, ph9_m3, ph9_p4, ab_2, pc_72, pc_73 : simd::cache_line_size())
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

        pc_72[k] = e_0 * fs_225_8_105 * h3_m3 + e_1 * fs_255_4_3 * h5_m5 - e_1 * fs_85_4_15 * h5_m3 - e_1 * fs_25_1_105 * r_2 * h3_m3 - e_2 * fs_1653_572_22 * h7_m5 + e_2 * fs_10887_572_2 * h7_m3 - e_2 * fs_765_26_3 * r_2 * h5_m5 + e_2 * fs_255_26_15 * r_2 * h5_m3 + e_2 * fs_75_11_105 * r_4 * h3_m3 + e_3 * fs_21_221_3003 * h9_m5 + e_3 * fs_63_221_55 * h9_m3 + e_3 * fs_1653_2431_22 * r_2 * h7_m5 - e_3 * fs_10887_2431_2 * r_2 * h7_m3 + e_3 * fs_51_13_3 * r_4 * h5_m5 - e_3 * fs_17_13_15 * r_4 * h5_m3 - e_3 * fs_100_143_105 * r_6 * h3_m3 - e_4 * fs_2_221_3003 * r_2 * h9_m5 - e_4 * fs_6_221_55 * r_2 * h9_m3 - e_4 * fs_87_2431_22 * r_4 * h7_m5 + e_4 * fs_573_2431_2 * r_4 * h7_m3 - e_4 * fs_2_13_3 * r_6 * h5_m5 + e_4 * fs_2_39_15 * r_6 * h5_m3 + e_4 * fs_10_429_105 * r_8 * h3_m3;

        pc_73[k] = - e_1 * fs_255_4_5 * h5_p4 + e_2 * fs_456_143_33 * h7_p4 + e_2 * fs_765_26_5 * r_2 * h5_p4 + e_3 * fs_42_221_715 * h9_p4 - e_3 * fs_1824_2431_33 * r_2 * h7_p4 - e_3 * fs_51_13_5 * r_4 * h5_p4 - e_4 * fs_4_221_715 * r_2 * h9_p4 + e_4 * fs_96_2431_33 * r_4 * h7_p4 + e_4 * fs_2_13_5 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ab_2, pc_74 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];

        pc_74[k] = - e_0 * fs_225_8_105 * h3_p3 + e_1 * fs_85_4_15 * h5_p3 + e_1 * fs_255_4_3 * h5_p5 + e_1 * fs_25_1_105 * r_2 * h3_p3 - e_2 * fs_10887_572_2 * h7_p3 - e_2 * fs_1653_572_22 * h7_p5 - e_2 * fs_255_26_15 * r_2 * h5_p3 - e_2 * fs_765_26_3 * r_2 * h5_p5 - e_2 * fs_75_11_105 * r_4 * h3_p3 - e_3 * fs_63_221_55 * h9_p3 + e_3 * fs_21_221_3003 * h9_p5 + e_3 * fs_10887_2431_2 * r_2 * h7_p3 + e_3 * fs_1653_2431_22 * r_2 * h7_p5 + e_3 * fs_17_13_15 * r_4 * h5_p3 + e_3 * fs_51_13_3 * r_4 * h5_p5 + e_3 * fs_100_143_105 * r_6 * h3_p3 + e_4 * fs_6_221_55 * r_2 * h9_p3 - e_4 * fs_2_221_3003 * r_2 * h9_p5 - e_4 * fs_573_2431_2 * r_4 * h7_p3 - e_4 * fs_87_2431_22 * r_4 * h7_p5 - e_4 * fs_2_39_15 * r_6 * h5_p3 - e_4 * fs_2_13_3 * r_6 * h5_p5 - e_4 * fs_10_429_105 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph7_p6, ph7_p7, ph9_p1, ph9_p2, ph9_p6, ph9_p7, ab_2, pc_75, pc_76 : simd::cache_line_size())
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

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_75[k] = - e_0 * fs_675_4_7 * h3_p2 + e_1 * f_255_4 * h5_p2 + e_1 * fs_150_1_7 * r_2 * h3_p2 + e_2 * fs_1140_143_10 * h7_p2 - e_2 * fs_114_143_1430 * h7_p6 - e_2 * f_765_26 * r_2 * h5_p2 - e_2 * fs_450_11_7 * r_4 * h3_p2 + e_3 * fs_21_442_462 * h9_p2 + e_3 * fs_63_442_2002 * h9_p6 - e_3 * fs_4560_2431_10 * r_2 * h7_p2 + e_3 * fs_456_2431_1430 * r_2 * h7_p6 + e_3 * f_51_13 * r_4 * h5_p2 + e_3 * fs_600_143_7 * r_6 * h3_p2 - e_4 * fs_1_221_462 * r_2 * h9_p2 - e_4 * fs_3_221_2002 * r_2 * h9_p6 + e_4 * fs_240_2431_10 * r_4 * h7_p2 - e_4 * fs_24_2431_1430 * r_4 * h7_p6 - e_4 * f_2_13 * r_6 * h5_p2 - e_4 * fs_20_143_7 * r_8 * h3_p2;

        pc_76[k] = - e_0 * fs_225_8_105 * h3_p1 - e_1 * fs_85_4_42 * h5_p1 + e_1 * fs_25_1_105 * r_2 * h3_p1 - e_2 * fs_2565_572_10 * h7_p1 + e_2 * fs_57_572_30030 * h7_p7 + e_2 * fs_255_26_42 * r_2 * h5_p1 - e_2 * fs_75_11_105 * r_4 * h3_p1 - e_3 * fs_21_221_14 * h9_p1 + e_3 * fs_42_221_1001 * h9_p7 + e_3 * fs_2565_2431_10 * r_2 * h7_p1 - e_3 * fs_57_2431_30030 * r_2 * h7_p7 - e_3 * fs_17_13_42 * r_4 * h5_p1 + e_3 * fs_100_143_105 * r_6 * h3_p1 + e_4 * fs_2_221_14 * r_2 * h9_p1 - e_4 * fs_4_221_1001 * r_2 * h9_p7 - e_4 * fs_135_2431_10 * r_4 * h7_p1 + e_4 * fs_3_2431_30030 * r_4 * h7_p7 + e_4 * fs_2_39_42 * r_6 * h5_p1 - e_4 * fs_10_429_105 * r_8 * h3_p1;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m7, ph7_m3, ph7_m2, ph9_m8, ph9_m7, ph9_m3, ph9_m2, ab_2, pc_77, pc_78 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_77[k] = e_0 * fs_225_8_231 * h3_m2 + e_1 * fs_85_4_33 * h5_m2 - e_1 * fs_25_1_231 * r_2 * h3_m2 + e_2 * fs_285_572_330 * h7_m2 - e_2 * fs_255_26_33 * r_2 * h5_m2 + e_2 * fs_75_11_231 * r_4 * h3_m2 + e_3 * fs_42_221_1547 * h9_m8 + e_3 * fs_21_442_14 * h9_m2 - e_3 * fs_285_2431_330 * r_2 * h7_m2 + e_3 * fs_17_13_33 * r_4 * h5_m2 - e_3 * fs_100_143_231 * r_6 * h3_m2 - e_4 * fs_4_221_1547 * r_2 * h9_m8 - e_4 * fs_1_221_14 * r_2 * h9_m2 + e_4 * fs_15_2431_330 * r_4 * h7_m2 - e_4 * fs_2_39_33 * r_6 * h5_m2 + e_4 * fs_10_429_231 * r_8 * h3_m2;

        pc_78[k] = e_0 * fs_225_8_231 * h3_m3 - e_1 * fs_85_4_33 * h5_m3 - e_1 * fs_25_1_231 * r_2 * h3_m3 - e_2 * fs_57_52_910 * h7_m7 - e_2 * fs_1083_572_110 * h7_m3 + e_2 * fs_255_26_33 * r_2 * h5_m3 + e_2 * fs_75_11_231 * r_4 * h3_m3 + e_3 * fs_84_221_273 * h9_m7 - e_3 * f_126_221 * h9_m3 + e_3 * fs_57_221_910 * r_2 * h7_m7 + e_3 * fs_1083_2431_110 * r_2 * h7_m3 - e_3 * fs_17_13_33 * r_4 * h5_m3 - e_3 * fs_100_143_231 * r_6 * h3_m3 - e_4 * fs_8_221_273 * r_2 * h9_m7 + e_4 * f_12_221 * r_2 * h9_m3 - e_4 * fs_3_221_910 * r_4 * h7_m7 - e_4 * fs_57_2431_110 * r_4 * h7_m3 + e_4 * fs_2_39_33 * r_6 * h5_m3 + e_4 * fs_10_429_231 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p5, ph7_m6, ph7_m4, ph7_p4, ph7_p5, ph7_p6, ph9_m6, ph9_m4, ph9_p4, ph9_p5, ph9_p6, ab_2, pc_79, pc_80, pc_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_79[k] = e_2 * fs_57_52_26 * h7_m6 + e_2 * f_57_2 * h7_m4 + e_3 * fs_63_442_910 * h9_m6 + e_3 * fs_21_221_195 * h9_m4 - e_3 * fs_57_221_26 * r_2 * h7_m6 - e_3 * f_114_17 * r_2 * h7_m4 - e_4 * fs_3_221_910 * r_2 * h9_m6 - e_4 * fs_2_221_195 * r_2 * h9_m4 + e_4 * fs_3_221_26 * r_4 * h7_m6 + e_4 * f_6_17 * r_4 * h7_m4;

        pc_80[k] = - e_1 * fs_255_4_11 * h5_p5 + e_2 * fs_399_26_6 * h7_p5 + e_2 * fs_765_26_11 * r_2 * h5_p5 + e_3 * fs_84_221_91 * h9_p5 - e_3 * fs_798_221_6 * r_2 * h7_p5 - e_3 * fs_51_13_11 * r_4 * h5_p5 - e_4 * fs_8_221_91 * r_2 * h9_p5 + e_4 * fs_42_221_6 * r_4 * h7_p5 + e_4 * fs_2_13_11 * r_6 * h5_p5;

        pc_81[k] = - e_2 * f_57_2 * h7_p4 + e_2 * fs_57_52_26 * h7_p6 - e_3 * fs_21_221_195 * h9_p4 + e_3 * fs_63_442_910 * h9_p6 + e_3 * f_114_17 * r_2 * h7_p4 - e_3 * fs_57_221_26 * r_2 * h7_p6 + e_4 * fs_2_221_195 * r_2 * h9_p4 - e_4 * fs_3_221_910 * r_2 * h9_p6 - e_4 * f_6_17 * r_4 * h7_p4 + e_4 * fs_3_221_26 * r_4 * h7_p6;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph7_p7, ph9_p2, ph9_p3, ph9_p7, ph9_p8, ab_2, pc_82, pc_83 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h9_p8 = ph9_p8[k];

        pc_82[k] = - e_0 * fs_225_8_231 * h3_p3 + e_1 * fs_85_4_33 * h5_p3 + e_1 * fs_25_1_231 * r_2 * h3_p3 + e_2 * fs_1083_572_110 * h7_p3 - e_2 * fs_57_52_910 * h7_p7 - e_2 * fs_255_26_33 * r_2 * h5_p3 - e_2 * fs_75_11_231 * r_4 * h3_p3 + e_3 * f_126_221 * h9_p3 + e_3 * fs_84_221_273 * h9_p7 - e_3 * fs_1083_2431_110 * r_2 * h7_p3 + e_3 * fs_57_221_910 * r_2 * h7_p7 + e_3 * fs_17_13_33 * r_4 * h5_p3 + e_3 * fs_100_143_231 * r_6 * h3_p3 - e_4 * f_12_221 * r_2 * h9_p3 - e_4 * fs_8_221_273 * r_2 * h9_p7 + e_4 * fs_57_2431_110 * r_4 * h7_p3 - e_4 * fs_3_221_910 * r_4 * h7_p7 - e_4 * fs_2_39_33 * r_6 * h5_p3 - e_4 * fs_10_429_231 * r_8 * h3_p3;

        pc_83[k] = - e_0 * fs_225_8_231 * h3_p2 - e_1 * fs_85_4_33 * h5_p2 + e_1 * fs_25_1_231 * r_2 * h3_p2 - e_2 * fs_285_572_330 * h7_p2 + e_2 * fs_255_26_33 * r_2 * h5_p2 - e_2 * fs_75_11_231 * r_4 * h3_p2 - e_3 * fs_21_442_14 * h9_p2 + e_3 * fs_42_221_1547 * h9_p8 + e_3 * fs_285_2431_330 * r_2 * h7_p2 - e_3 * fs_17_13_33 * r_4 * h5_p2 + e_3 * fs_100_143_231 * r_6 * h3_p2 + e_4 * fs_1_221_14 * r_2 * h9_p2 - e_4 * fs_4_221_1547 * r_2 * h9_p8 - e_4 * fs_15_2431_330 * r_4 * h7_p2 + e_4 * fs_2_39_33 * r_6 * h5_p2 - e_4 * fs_10_429_231 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m4, ph5_m3, ph7_m4, ph7_m3, ph9_m9, ph9_m8, ph9_m4, ph9_m3, ab_2, pc_84, pc_85 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_84[k] = e_0 * fs_225_8_462 * h3_m3 + e_1 * fs_85_8_66 * h5_m3 - e_1 * fs_25_1_462 * r_2 * h3_m3 + e_2 * fs_171_286_55 * h7_m3 - e_2 * fs_255_52_66 * r_2 * h5_m3 + e_2 * fs_75_11_462 * r_4 * h3_m3 + e_3 * fs_21_221_9282 * h9_m9 + e_3 * fs_21_442_2 * h9_m3 - e_3 * fs_342_2431_55 * r_2 * h7_m3 + e_3 * fs_17_26_66 * r_4 * h5_m3 - e_3 * fs_100_143_462 * r_6 * h3_m3 - e_4 * fs_2_221_9282 * r_2 * h9_m9 - e_4 * fs_1_221_2 * r_2 * h9_m3 + e_4 * fs_18_2431_55 * r_4 * h7_m3 - e_4 * fs_1_39_66 * r_6 * h5_m3 + e_4 * fs_10_429_462 * r_8 * h3_m3;

        pc_85[k] = - e_1 * fs_255_8_22 * h5_m4 - e_2 * fs_57_26_30 * h7_m4 + e_2 * fs_765_52_22 * r_2 * h5_m4 + e_3 * fs_21_221_3094 * h9_m8 - e_3 * fs_21_442_26 * h9_m4 + e_3 * fs_114_221_30 * r_2 * h7_m4 - e_3 * fs_51_26_22 * r_4 * h5_m4 - e_4 * fs_2_221_3094 * r_2 * h9_m8 + e_4 * fs_1_221_26 * r_2 * h9_m4 - e_4 * fs_6_221_30 * r_4 * h7_m4 + e_4 * fs_1_13_22 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m5, ph5_p5, ph7_m7, ph7_m5, ph7_p5, ph7_p6, ph7_p7, ph9_m7, ph9_m5, ph9_p5, ph9_p6, ph9_p7, ab_2, pc_86, pc_87, pc_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m5 = ph5_m5[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_86[k] = e_1 * fs_255_8_22 * h5_m5 + e_2 * fs_57_26_273 * h7_m7 + e_2 * fs_171_13_3 * h7_m5 - e_2 * fs_765_52_22 * r_2 * h5_m5 + e_3 * fs_21_221_910 * h9_m7 + e_3 * fs_21_442_182 * h9_m5 - e_3 * fs_114_221_273 * r_2 * h7_m7 - e_3 * fs_684_221_3 * r_2 * h7_m5 + e_3 * fs_51_26_22 * r_4 * h5_m5 - e_4 * fs_2_221_910 * r_2 * h9_m7 - e_4 * fs_1_221_182 * r_2 * h9_m5 + e_4 * fs_6_221_273 * r_4 * h7_m7 + e_4 * fs_36_221_3 * r_4 * h7_m5 - e_4 * fs_1_13_22 * r_6 * h5_m5;

        pc_87[k] = e_2 * fs_171_13_13 * h7_p6 + e_3 * fs_21_221_455 * h9_p6 - e_3 * fs_684_221_13 * r_2 * h7_p6 - e_4 * fs_2_221_455 * r_2 * h9_p6 + e_4 * fs_36_221_13 * r_4 * h7_p6;

        pc_88[k] = - e_1 * fs_255_8_22 * h5_p5 - e_2 * fs_171_13_3 * h7_p5 + e_2 * fs_57_26_273 * h7_p7 + e_2 * fs_765_52_22 * r_2 * h5_p5 - e_3 * fs_21_442_182 * h9_p5 + e_3 * fs_21_221_910 * h9_p7 + e_3 * fs_684_221_3 * r_2 * h7_p5 - e_3 * fs_114_221_273 * r_2 * h7_p7 - e_3 * fs_51_26_22 * r_4 * h5_p5 + e_4 * fs_1_221_182 * r_2 * h9_p5 - e_4 * fs_2_221_910 * r_2 * h9_p7 - e_4 * fs_36_221_3 * r_4 * h7_p5 + e_4 * fs_6_221_273 * r_4 * h7_p7 + e_4 * fs_1_13_22 * r_6 * h5_p5;
    }

    // NOTE: the angular components are formed in 50 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph9_p3, ph9_p4, ph9_p8, ph9_p9, ab_2, pc_89, pc_90 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h9_p9 = ph9_p9[k];

        pc_89[k] = e_1 * fs_255_8_22 * h5_p4 + e_2 * fs_57_26_30 * h7_p4 - e_2 * fs_765_52_22 * r_2 * h5_p4 + e_3 * fs_21_442_26 * h9_p4 + e_3 * fs_21_221_3094 * h9_p8 - e_3 * fs_114_221_30 * r_2 * h7_p4 + e_3 * fs_51_26_22 * r_4 * h5_p4 - e_4 * fs_1_221_26 * r_2 * h9_p4 - e_4 * fs_2_221_3094 * r_2 * h9_p8 + e_4 * fs_6_221_30 * r_4 * h7_p4 - e_4 * fs_1_13_22 * r_6 * h5_p4;

        pc_90[k] = - e_0 * fs_225_8_462 * h3_p3 - e_1 * fs_85_8_66 * h5_p3 + e_1 * fs_25_1_462 * r_2 * h3_p3 - e_2 * fs_171_286_55 * h7_p3 + e_2 * fs_255_52_66 * r_2 * h5_p3 - e_2 * fs_75_11_462 * r_4 * h3_p3 - e_3 * fs_21_442_2 * h9_p3 + e_3 * fs_21_221_9282 * h9_p9 + e_3 * fs_342_2431_55 * r_2 * h7_p3 - e_3 * fs_17_26_66 * r_4 * h5_p3 + e_3 * fs_100_143_462 * r_6 * h3_p3 + e_4 * fs_1_221_2 * r_2 * h9_p3 - e_4 * fs_2_221_9282 * r_2 * h9_p9 - e_4 * fs_18_2431_55 * r_4 * h7_p3 + e_4 * fs_1_39_66 * r_6 * h5_p3 - e_4 * fs_10_429_462 * r_8 * h3_p3;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[91] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90};

    for (size_t n = 0; n < 91; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
