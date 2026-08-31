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



#include "SimdKineticEnergyRecHF.hpp"

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
compute_hf_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hf_kinetic_energy: Basis functions must be of angular momenta five and three"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hf_kinetic_energy: Harmonics must reach angular momentum 8"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hf_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 77 * nvalues, 0.0);

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

            const auto ff_0 = fbase * bexp * bexp * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * bexp * bexp * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_100_231 = 100.0 / 231.0;
    const auto f_1064_143 = 1064.0 / 143.0;
    const auto f_112_143 = 112.0 / 143.0;
    const auto f_135_1 = 135.0;
    const auto f_135_2 = 135.0 / 2.0;
    const auto f_14_33 = 14.0 / 33.0;
    const auto f_1620_143 = 1620.0 / 143.0;
    const auto f_1950_7 = 1950.0 / 7.0;
    const auto f_238_33 = 238.0 / 33.0;
    const auto f_255_11 = 255.0 / 11.0;
    const auto f_255_22 = 255.0 / 22.0;
    const auto f_2600_231 = 2600.0 / 231.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_315_4 = 315.0 / 4.0;
    const auto f_34_11 = 34.0 / 11.0;
    const auto f_360_1001 = 360.0 / 1001.0;
    const auto f_36_143 = 36.0 / 143.0;
    const auto f_399_143 = 399.0 / 143.0;
    const auto f_4050_77 = 4050.0 / 77.0;
    const auto f_405_11 = 405.0 / 11.0;
    const auto f_42_143 = 42.0 / 143.0;
    const auto f_4_11 = 4.0 / 11.0;
    const auto f_595_22 = 595.0 / 22.0;
    const auto f_650_7 = 650.0 / 7.0;
    const auto f_675_7 = 675.0 / 7.0;
    const auto f_68_11 = 68.0 / 11.0;
    const auto f_72_143 = 72.0 / 143.0;
    const auto f_8100_1001 = 8100.0 / 1001.0;
    const auto f_810_11 = 810.0 / 11.0;
    const auto f_810_143 = 810.0 / 143.0;
    const auto f_945_143 = 945.0 / 143.0;
    const auto f_945_22 = 945.0 / 22.0;
    const auto f_975_4 = 975.0 / 4.0;
    const auto fs_1040_231_5 = std::sqrt(5408000.0 / 53361.0);
    const auto fs_1080_77_7 = std::sqrt(1166400.0 / 847.0);
    const auto fs_10_143_66 = std::sqrt(600.0 / 1859.0);
    const auto fs_10_231_105 = std::sqrt(500.0 / 2541.0);
    const auto fs_10_231_35 = std::sqrt(500.0 / 7623.0);
    const auto fs_10_231_42 = std::sqrt(200.0 / 2541.0);
    const auto fs_10_231_5 = std::sqrt(500.0 / 53361.0);
    const auto fs_10_77_10 = std::sqrt(1000.0 / 5929.0);
    const auto fs_10_77_7 = std::sqrt(100.0 / 847.0);
    const auto fs_1125_56_10 = std::sqrt(6328125.0 / 1568.0);
    const auto fs_114_1001_6 = std::sqrt(77976.0 / 1002001.0);
    const auto fs_114_143_42 = std::sqrt(545832.0 / 20449.0);
    const auto fs_119_33_2 = std::sqrt(28322.0 / 1089.0);
    const auto fs_1215_1001_21 = std::sqrt(4428675.0 / 143143.0);
    const auto fs_1215_1001_35 = std::sqrt(7381125.0 / 143143.0);
    const auto fs_1215_154_21 = std::sqrt(4428675.0 / 3388.0);
    const auto fs_1215_154_35 = std::sqrt(7381125.0 / 3388.0);
    const auto fs_129_1001_2 = std::sqrt(33282.0 / 1002001.0);
    const auto fs_12_1001_21 = std::sqrt(432.0 / 143143.0);
    const auto fs_12_143_15 = std::sqrt(2160.0 / 20449.0);
    const auto fs_12_143_42 = std::sqrt(6048.0 / 20449.0);
    const auto fs_12_143_5 = std::sqrt(720.0 / 20449.0);
    const auto fs_1300_231_2 = std::sqrt(3380000.0 / 53361.0);
    const auto fs_130_231_14 = std::sqrt(33800.0 / 7623.0);
    const auto fs_130_231_2 = std::sqrt(33800.0 / 53361.0);
    const auto fs_130_231_210 = std::sqrt(169000.0 / 2541.0);
    const auto fs_130_231_30 = std::sqrt(169000.0 / 17787.0);
    const auto fs_130_7_14 = std::sqrt(33800.0 / 7.0);
    const auto fs_130_7_3 = std::sqrt(50700.0 / 49.0);
    const auto fs_130_7_7 = std::sqrt(16900.0 / 7.0);
    const auto fs_1350_1001_7 = std::sqrt(1822500.0 / 143143.0);
    const auto fs_135_1001_105 = std::sqrt(273375.0 / 143143.0);
    const auto fs_135_1001_15 = std::sqrt(273375.0 / 1002001.0);
    const auto fs_135_1001_210 = std::sqrt(546750.0 / 143143.0);
    const auto fs_135_11_15 = std::sqrt(273375.0 / 121.0);
    const auto fs_135_11_5 = std::sqrt(91125.0 / 121.0);
    const auto fs_135_143_3 = std::sqrt(54675.0 / 20449.0);
    const auto fs_135_14_10 = std::sqrt(91125.0 / 98.0);
    const auto fs_135_154_105 = std::sqrt(273375.0 / 3388.0);
    const auto fs_135_154_15 = std::sqrt(273375.0 / 23716.0);
    const auto fs_135_154_210 = std::sqrt(273375.0 / 1694.0);
    const auto fs_135_182_70 = std::sqrt(91125.0 / 2366.0);
    const auto fs_135_22_3 = std::sqrt(54675.0 / 484.0);
    const auto fs_135_28_70 = std::sqrt(91125.0 / 56.0);
    const auto fs_135_77_21 = std::sqrt(54675.0 / 847.0);
    const auto fs_135_91_10 = std::sqrt(182250.0 / 8281.0);
    const auto fs_136_33_2 = std::sqrt(36992.0 / 1089.0);
    const auto fs_152_143_7 = std::sqrt(161728.0 / 20449.0);
    const auto fs_15_143_6 = std::sqrt(1350.0 / 20449.0);
    const auto fs_16_143_7 = std::sqrt(1792.0 / 20449.0);
    const auto fs_170_11_2 = std::sqrt(57800.0 / 121.0);
    const auto fs_17_11_10 = std::sqrt(2890.0 / 121.0);
    const auto fs_17_11_22 = std::sqrt(578.0 / 11.0);
    const auto fs_17_33_30 = std::sqrt(2890.0 / 363.0);
    const auto fs_17_33_6 = std::sqrt(578.0 / 363.0);
    const auto fs_17_33_66 = std::sqrt(578.0 / 33.0);
    const auto fs_180_1001_7 = std::sqrt(32400.0 / 143143.0);
    const auto fs_180_7_7 = std::sqrt(32400.0 / 7.0);
    const auto fs_1935_56_2 = std::sqrt(3744225.0 / 1568.0);
    const auto fs_195_14_14 = std::sqrt(38025.0 / 14.0);
    const auto fs_195_14_2 = std::sqrt(38025.0 / 98.0);
    const auto fs_195_14_210 = std::sqrt(570375.0 / 14.0);
    const auto fs_195_14_30 = std::sqrt(570375.0 / 98.0);
    const auto fs_195_16_14 = std::sqrt(266175.0 / 128.0);
    const auto fs_195_16_2 = std::sqrt(38025.0 / 128.0);
    const auto fs_195_16_210 = std::sqrt(3992625.0 / 128.0);
    const auto fs_195_16_30 = std::sqrt(570375.0 / 128.0);
    const auto fs_195_2_5 = std::sqrt(190125.0 / 4.0);
    const auto fs_195_4_14 = std::sqrt(266175.0 / 8.0);
    const auto fs_195_4_3 = std::sqrt(114075.0 / 16.0);
    const auto fs_195_4_7 = std::sqrt(266175.0 / 16.0);
    const auto fs_195_7_10 = std::sqrt(380250.0 / 49.0);
    const auto fs_195_7_105 = std::sqrt(570375.0 / 7.0);
    const auto fs_195_7_35 = std::sqrt(190125.0 / 7.0);
    const auto fs_195_7_42 = std::sqrt(228150.0 / 7.0);
    const auto fs_195_7_5 = std::sqrt(190125.0 / 49.0);
    const auto fs_195_7_7 = std::sqrt(38025.0 / 7.0);
    const auto fs_195_8_105 = std::sqrt(3992625.0 / 64.0);
    const auto fs_195_8_35 = std::sqrt(1330875.0 / 64.0);
    const auto fs_195_8_42 = std::sqrt(798525.0 / 32.0);
    const auto fs_195_8_5 = std::sqrt(190125.0 / 64.0);
    const auto fs_19_143_105 = std::sqrt(37905.0 / 20449.0);
    const auto fs_19_143_1155 = std::sqrt(37905.0 / 1859.0);
    const auto fs_19_143_2002 = std::sqrt(5054.0 / 143.0);
    const auto fs_19_143_286 = std::sqrt(722.0 / 143.0);
    const auto fs_19_143_30 = std::sqrt(10830.0 / 20449.0);
    const auto fs_19_143_33 = std::sqrt(1083.0 / 1859.0);
    const auto fs_19_143_42 = std::sqrt(15162.0 / 20449.0);
    const auto fs_19_143_462 = std::sqrt(15162.0 / 1859.0);
    const auto fs_19_143_858 = std::sqrt(2166.0 / 143.0);
    const auto fs_19_286_10010 = std::sqrt(12635.0 / 286.0);
    const auto fs_19_286_14 = std::sqrt(2527.0 / 40898.0);
    const auto fs_19_286_2 = std::sqrt(361.0 / 40898.0);
    const auto fs_19_286_2002 = std::sqrt(2527.0 / 286.0);
    const auto fs_19_286_22 = std::sqrt(361.0 / 3718.0);
    const auto fs_19_286_4290 = std::sqrt(5415.0 / 286.0);
    const auto fs_19_286_6006 = std::sqrt(7581.0 / 286.0);
    const auto fs_1_11_10 = std::sqrt(10.0 / 121.0);
    const auto fs_1_11_22 = std::sqrt(2.0 / 11.0);
    const auto fs_1_143_10010 = std::sqrt(70.0 / 143.0);
    const auto fs_1_143_14 = std::sqrt(14.0 / 20449.0);
    const auto fs_1_143_2 = std::sqrt(2.0 / 20449.0);
    const auto fs_1_143_2002 = std::sqrt(14.0 / 143.0);
    const auto fs_1_143_22 = std::sqrt(2.0 / 1859.0);
    const auto fs_1_143_4290 = std::sqrt(30.0 / 143.0);
    const auto fs_1_143_6006 = std::sqrt(42.0 / 143.0);
    const auto fs_1_33_30 = std::sqrt(10.0 / 363.0);
    const auto fs_1_33_6 = std::sqrt(2.0 / 363.0);
    const auto fs_1_33_66 = std::sqrt(2.0 / 33.0);
    const auto fs_2025_1001_7 = std::sqrt(4100625.0 / 143143.0);
    const auto fs_2025_154_7 = std::sqrt(4100625.0 / 3388.0);
    const auto fs_2025_2002_30 = std::sqrt(61509375.0 / 2004002.0);
    const auto fs_2025_308_30 = std::sqrt(61509375.0 / 47432.0);
    const auto fs_2025_77_7 = std::sqrt(4100625.0 / 847.0);
    const auto fs_20_231_14 = std::sqrt(800.0 / 7623.0);
    const auto fs_20_231_3 = std::sqrt(400.0 / 17787.0);
    const auto fs_20_231_7 = std::sqrt(400.0 / 7623.0);
    const auto fs_2160_1001_7 = std::sqrt(4665600.0 / 143143.0);
    const auto fs_225_14_7 = std::sqrt(50625.0 / 28.0);
    const auto fs_225_28_105 = std::sqrt(759375.0 / 112.0);
    const auto fs_225_4_3 = std::sqrt(151875.0 / 16.0);
    const auto fs_228_143_11 = std::sqrt(51984.0 / 1859.0);
    const auto fs_240_1001_3 = std::sqrt(172800.0 / 1002001.0);
    const auto fs_24_1001_105 = std::sqrt(8640.0 / 143143.0);
    const auto fs_24_143_11 = std::sqrt(576.0 / 1859.0);
    const auto fs_24_143_3 = std::sqrt(1728.0 / 20449.0);
    const auto fs_255_22_11 = std::sqrt(65025.0 / 44.0);
    const auto fs_255_22_3 = std::sqrt(195075.0 / 484.0);
    const auto fs_255_22_5 = std::sqrt(325125.0 / 484.0);
    const auto fs_255_22_7 = std::sqrt(455175.0 / 484.0);
    const auto fs_255_44_10 = std::sqrt(325125.0 / 968.0);
    const auto fs_255_44_22 = std::sqrt(65025.0 / 88.0);
    const auto fs_2565_1001_6 = std::sqrt(39475350.0 / 1002001.0);
    const auto fs_2565_154_6 = std::sqrt(19737675.0 / 11858.0);
    const auto fs_260_231_105 = std::sqrt(338000.0 / 2541.0);
    const auto fs_260_231_35 = std::sqrt(338000.0 / 7623.0);
    const auto fs_260_231_42 = std::sqrt(135200.0 / 2541.0);
    const auto fs_260_231_5 = std::sqrt(338000.0 / 53361.0);
    const auto fs_260_77_10 = std::sqrt(676000.0 / 5929.0);
    const auto fs_260_77_7 = std::sqrt(67600.0 / 847.0);
    const auto fs_260_7_5 = std::sqrt(338000.0 / 49.0);
    const auto fs_266_143_10 = std::sqrt(707560.0 / 20449.0);
    const auto fs_266_143_15 = std::sqrt(1061340.0 / 20449.0);
    const auto fs_2700_77_3 = std::sqrt(21870000.0 / 5929.0);
    const auto fs_270_1001_21 = std::sqrt(218700.0 / 143143.0);
    const auto fs_270_11_3 = std::sqrt(218700.0 / 121.0);
    const auto fs_270_143_15 = std::sqrt(1093500.0 / 20449.0);
    const auto fs_270_143_5 = std::sqrt(364500.0 / 20449.0);
    const auto fs_270_77_105 = std::sqrt(1093500.0 / 847.0);
    const auto fs_285_143_7 = std::sqrt(568575.0 / 20449.0);
    const auto fs_285_286_6 = std::sqrt(243675.0 / 40898.0);
    const auto fs_28_143_10 = std::sqrt(7840.0 / 20449.0);
    const auto fs_28_143_15 = std::sqrt(11760.0 / 20449.0);
    const auto fs_2_11_11 = std::sqrt(4.0 / 11.0);
    const auto fs_2_11_3 = std::sqrt(12.0 / 121.0);
    const auto fs_2_11_5 = std::sqrt(20.0 / 121.0);
    const auto fs_2_11_7 = std::sqrt(28.0 / 121.0);
    const auto fs_2_143_105 = std::sqrt(420.0 / 20449.0);
    const auto fs_2_143_1155 = std::sqrt(420.0 / 1859.0);
    const auto fs_2_143_2002 = std::sqrt(56.0 / 143.0);
    const auto fs_2_143_286 = std::sqrt(8.0 / 143.0);
    const auto fs_2_143_30 = std::sqrt(120.0 / 20449.0);
    const auto fs_2_143_33 = std::sqrt(12.0 / 1859.0);
    const auto fs_2_143_42 = std::sqrt(168.0 / 20449.0);
    const auto fs_2_143_462 = std::sqrt(168.0 / 1859.0);
    const auto fs_2_143_858 = std::sqrt(24.0 / 143.0);
    const auto fs_2_33_14 = std::sqrt(56.0 / 1089.0);
    const auto fs_2_33_15 = std::sqrt(20.0 / 363.0);
    const auto fs_2_33_21 = std::sqrt(28.0 / 363.0);
    const auto fs_2_33_3 = std::sqrt(4.0 / 363.0);
    const auto fs_2_33_33 = std::sqrt(4.0 / 33.0);
    const auto fs_2_33_35 = std::sqrt(140.0 / 1089.0);
    const auto fs_2_33_42 = std::sqrt(56.0 / 363.0);
    const auto fs_2_33_6 = std::sqrt(8.0 / 363.0);
    const auto fs_30_1001_105 = std::sqrt(13500.0 / 143143.0);
    const auto fs_30_143_3 = std::sqrt(2700.0 / 20449.0);
    const auto fs_30_143_7 = std::sqrt(6300.0 / 20449.0);
    const auto fs_325_7_2 = std::sqrt(211250.0 / 49.0);
    const auto fs_3375_2002_10 = std::sqrt(56953125.0 / 2004002.0);
    const auto fs_3375_308_10 = std::sqrt(56953125.0 / 47432.0);
    const auto fs_34_11_11 = std::sqrt(1156.0 / 11.0);
    const auto fs_34_11_3 = std::sqrt(3468.0 / 121.0);
    const auto fs_34_11_5 = std::sqrt(5780.0 / 121.0);
    const auto fs_34_11_7 = std::sqrt(8092.0 / 121.0);
    const auto fs_34_33_14 = std::sqrt(16184.0 / 1089.0);
    const auto fs_34_33_15 = std::sqrt(5780.0 / 363.0);
    const auto fs_34_33_21 = std::sqrt(8092.0 / 363.0);
    const auto fs_34_33_3 = std::sqrt(1156.0 / 363.0);
    const auto fs_34_33_33 = std::sqrt(1156.0 / 33.0);
    const auto fs_34_33_35 = std::sqrt(40460.0 / 1089.0);
    const auto fs_34_33_42 = std::sqrt(16184.0 / 363.0);
    const auto fs_34_33_6 = std::sqrt(2312.0 / 363.0);
    const auto fs_380_143_6 = std::sqrt(866400.0 / 20449.0);
    const auto fs_38_143_1001 = std::sqrt(10108.0 / 143.0);
    const auto fs_38_143_7 = std::sqrt(10108.0 / 20449.0);
    const auto fs_390_7_14 = std::sqrt(304200.0 / 7.0);
    const auto fs_390_7_3 = std::sqrt(456300.0 / 49.0);
    const auto fs_390_7_7 = std::sqrt(152100.0 / 7.0);
    const auto fs_399_143_6 = std::sqrt(955206.0 / 20449.0);
    const auto fs_3_143_42 = std::sqrt(378.0 / 20449.0);
    const auto fs_3_143_858 = std::sqrt(54.0 / 143.0);
    const auto fs_3_91_70 = std::sqrt(90.0 / 1183.0);
    const auto fs_4050_1001_7 = std::sqrt(16402500.0 / 143143.0);
    const auto fs_405_28_21 = std::sqrt(492075.0 / 112.0);
    const auto fs_405_28_35 = std::sqrt(820125.0 / 112.0);
    const auto fs_40_143_6 = std::sqrt(9600.0 / 20449.0);
    const auto fs_40_231_5 = std::sqrt(8000.0 / 53361.0);
    const auto fs_42_143_6 = std::sqrt(10584.0 / 20449.0);
    const auto fs_450_7_3 = std::sqrt(607500.0 / 49.0);
    const auto fs_45_1001_30 = std::sqrt(60750.0 / 1002001.0);
    const auto fs_45_14_21 = std::sqrt(6075.0 / 28.0);
    const auto fs_45_1_3 = std::sqrt(6075.0);
    const auto fs_45_28_105 = std::sqrt(30375.0 / 112.0);
    const auto fs_45_28_15 = std::sqrt(30375.0 / 784.0);
    const auto fs_45_28_210 = std::sqrt(30375.0 / 56.0);
    const auto fs_45_2_15 = std::sqrt(30375.0 / 4.0);
    const auto fs_45_2_5 = std::sqrt(10125.0 / 4.0);
    const auto fs_45_4_3 = std::sqrt(6075.0 / 16.0);
    const auto fs_45_7_105 = std::sqrt(30375.0 / 7.0);
    const auto fs_495_28_10 = std::sqrt(1225125.0 / 392.0);
    const auto fs_495_56_70 = std::sqrt(1225125.0 / 224.0);
    const auto fs_4_143_1001 = std::sqrt(112.0 / 143.0);
    const auto fs_4_143_7 = std::sqrt(112.0 / 20449.0);
    const auto fs_4_33_2 = std::sqrt(32.0 / 1089.0);
    const auto fs_4_33_21 = std::sqrt(112.0 / 363.0);
    const auto fs_4_33_7 = std::sqrt(112.0 / 1089.0);
    const auto fs_50_231_2 = std::sqrt(5000.0 / 53361.0);
    const auto fs_520_231_14 = std::sqrt(540800.0 / 7623.0);
    const auto fs_520_231_3 = std::sqrt(270400.0 / 17787.0);
    const auto fs_520_231_7 = std::sqrt(270400.0 / 7623.0);
    const auto fs_5400_1001_3 = std::sqrt(87480000.0 / 1002001.0);
    const auto fs_540_1001_105 = std::sqrt(4374000.0 / 143143.0);
    const auto fs_540_143_3 = std::sqrt(874800.0 / 20449.0);
    const auto fs_54_1001_21 = std::sqrt(8748.0 / 143143.0);
    const auto fs_54_1001_35 = std::sqrt(14580.0 / 143143.0);
    const auto fs_57_143_165 = std::sqrt(48735.0 / 1859.0);
    const auto fs_57_143_70 = std::sqrt(227430.0 / 20449.0);
    const auto fs_57_286_42 = std::sqrt(68229.0 / 40898.0);
    const auto fs_57_286_858 = std::sqrt(9747.0 / 286.0);
    const auto fs_5805_2002_2 = std::sqrt(33698025.0 / 2004002.0);
    const auto fs_5805_308_2 = std::sqrt(33698025.0 / 47432.0);
    const auto fs_585_7_10 = std::sqrt(3422250.0 / 49.0);
    const auto fs_585_7_7 = std::sqrt(342225.0 / 7.0);
    const auto fs_585_8_10 = std::sqrt(1711125.0 / 32.0);
    const auto fs_585_8_7 = std::sqrt(2395575.0 / 64.0);
    const auto fs_595_44_2 = std::sqrt(354025.0 / 968.0);
    const auto fs_5_143_22 = std::sqrt(50.0 / 1859.0);
    const auto fs_5_231_14 = std::sqrt(50.0 / 7623.0);
    const auto fs_5_231_2 = std::sqrt(50.0 / 53361.0);
    const auto fs_5_231_210 = std::sqrt(250.0 / 2541.0);
    const auto fs_5_231_30 = std::sqrt(250.0 / 17787.0);
    const auto fs_60_1001_7 = std::sqrt(3600.0 / 143143.0);
    const auto fs_65_14_14 = std::sqrt(4225.0 / 14.0);
    const auto fs_65_14_2 = std::sqrt(4225.0 / 98.0);
    const auto fs_65_14_210 = std::sqrt(63375.0 / 14.0);
    const auto fs_65_14_30 = std::sqrt(63375.0 / 98.0);
    const auto fs_65_7_105 = std::sqrt(63375.0 / 7.0);
    const auto fs_65_7_35 = std::sqrt(21125.0 / 7.0);
    const auto fs_65_7_42 = std::sqrt(25350.0 / 7.0);
    const auto fs_65_7_5 = std::sqrt(21125.0 / 49.0);
    const auto fs_675_1001_105 = std::sqrt(6834375.0 / 143143.0);
    const auto fs_675_143_3 = std::sqrt(1366875.0 / 20449.0);
    const auto fs_675_14_7 = std::sqrt(455625.0 / 28.0);
    const auto fs_675_154_105 = std::sqrt(6834375.0 / 3388.0);
    const auto fs_675_22_3 = std::sqrt(1366875.0 / 484.0);
    const auto fs_675_28_7 = std::sqrt(455625.0 / 112.0);
    const auto fs_675_56_30 = std::sqrt(6834375.0 / 1568.0);
    const auto fs_675_77_7 = std::sqrt(455625.0 / 847.0);
    const auto fs_68_33_2 = std::sqrt(9248.0 / 1089.0);
    const auto fs_68_33_21 = std::sqrt(32368.0 / 363.0);
    const auto fs_68_33_7 = std::sqrt(32368.0 / 1089.0);
    const auto fs_6_1001_105 = std::sqrt(540.0 / 143143.0);
    const auto fs_6_1001_15 = std::sqrt(540.0 / 1002001.0);
    const auto fs_6_1001_210 = std::sqrt(1080.0 / 143143.0);
    const auto fs_6_143_165 = std::sqrt(540.0 / 1859.0);
    const auto fs_6_143_3 = std::sqrt(108.0 / 20449.0);
    const auto fs_6_143_70 = std::sqrt(2520.0 / 20449.0);
    const auto fs_6_91_10 = std::sqrt(360.0 / 8281.0);
    const auto fs_75_1001_10 = std::sqrt(56250.0 / 1002001.0);
    const auto fs_76_143_55 = std::sqrt(28880.0 / 1859.0);
    const auto fs_780_7_5 = std::sqrt(3042000.0 / 49.0);
    const auto fs_7_33_2 = std::sqrt(98.0 / 1089.0);
    const auto fs_855_28_6 = std::sqrt(2193075.0 / 392.0);
    const auto fs_85_11_2 = std::sqrt(14450.0 / 121.0);
    const auto fs_85_11_21 = std::sqrt(151725.0 / 121.0);
    const auto fs_85_11_7 = std::sqrt(50575.0 / 121.0);
    const auto fs_85_22_14 = std::sqrt(50575.0 / 242.0);
    const auto fs_85_22_15 = std::sqrt(108375.0 / 484.0);
    const auto fs_85_22_21 = std::sqrt(151725.0 / 484.0);
    const auto fs_85_22_3 = std::sqrt(21675.0 / 484.0);
    const auto fs_85_22_33 = std::sqrt(21675.0 / 44.0);
    const auto fs_85_22_35 = std::sqrt(252875.0 / 484.0);
    const auto fs_85_22_42 = std::sqrt(151725.0 / 242.0);
    const auto fs_85_22_6 = std::sqrt(21675.0 / 242.0);
    const auto fs_85_44_30 = std::sqrt(108375.0 / 968.0);
    const auto fs_85_44_6 = std::sqrt(21675.0 / 968.0);
    const auto fs_85_44_66 = std::sqrt(21675.0 / 88.0);
    const auto fs_8_143_55 = std::sqrt(320.0 / 1859.0);
    const auto fs_8_33_2 = std::sqrt(128.0 / 1089.0);
    const auto fs_90_1001_7 = std::sqrt(8100.0 / 143143.0);
    const auto fs_95_143_66 = std::sqrt(54150.0 / 1859.0);
    const auto fs_95_286_22 = std::sqrt(9025.0 / 3718.0);
    const auto fs_96_1001_7 = std::sqrt(9216.0 / 143143.0);
    const auto fs_975_7_2 = std::sqrt(1901250.0 / 49.0);
    const auto fs_975_8_2 = std::sqrt(950625.0 / 32.0);

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph8_p7, ph8_p8, ab_2, pc_0, pc_1 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * fs_195_8_105 * h2_p2 + e_1 * fs_675_28_7 * h4_p2 - e_1 * fs_195_7_105 * r_2 * h2_p2 + e_2 * fs_85_44_6 * h6_p2 - e_2 * fs_2025_154_7 * r_2 * h4_p2 + e_2 * fs_65_7_105 * r_4 * h2_p2 + e_3 * fs_19_286_2 * h8_p2 + e_3 * fs_38_143_1001 * h8_p8 - e_3 * fs_17_33_6 * r_2 * h6_p2 + e_3 * fs_2025_1001_7 * r_4 * h4_p2 - e_3 * fs_260_231_105 * r_6 * h2_p2 - e_4 * fs_1_143_2 * r_2 * h8_p2 - e_4 * fs_4_143_1001 * r_2 * h8_p8 + e_4 * fs_1_33_6 * r_4 * h6_p2 - e_4 * fs_90_1001_7 * r_6 * h4_p2 + e_4 * fs_10_231_105 * r_8 * h2_p2;

        pc_1[k] = - e_1 * fs_225_4_3 * h4_p3 - e_2 * f_255_22 * h6_p3 + e_2 * fs_675_22_3 * r_2 * h4_p3 - e_3 * fs_19_286_22 * h8_p3 + e_3 * fs_19_286_6006 * h8_p7 + e_3 * f_34_11 * r_2 * h6_p3 - e_3 * fs_675_143_3 * r_4 * h4_p3 + e_4 * fs_1_143_22 * r_2 * h8_p3 - e_4 * fs_1_143_6006 * r_2 * h8_p7 - e_4 * f_2_11 * r_4 * h6_p3 + e_4 * fs_30_143_3 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_p4, ph6_m6, ph6_m5, ph6_m4, ph6_p4, ph6_p6, ph8_m6, ph8_m5, ph8_m4, ph8_p4, ph8_p6, ab_2, pc_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];

        pc_2[k] = e_1 * fs_45_2_15 * h4_p4 + e_2 * fs_255_22_3 * h6_p4 + e_2 * fs_255_44_22 * h6_p6 - e_2 * fs_135_11_15 * r_2 * h4_p4 + e_3 * fs_19_143_33 * h8_p4 + e_3 * fs_19_286_2002 * h8_p6 - e_3 * fs_34_11_3 * r_2 * h6_p4 - e_3 * fs_17_11_22 * r_2 * h6_p6 + e_3 * fs_270_143_15 * r_4 * h4_p4 - e_4 * fs_2_143_33 * r_2 * h8_p4 - e_4 * fs_1_143_2002 * r_2 * h8_p6 + e_4 * fs_2_11_3 * r_4 * h6_p4 + e_4 * fs_1_11_22 * r_4 * h6_p6 - e_4 * fs_12_143_15 * r_6 * h4_p4;

        pc_3[k] = - e_2 * fs_255_22_11 * h6_m5 - e_3 * fs_19_143_286 * h8_m5 + e_3 * fs_34_11_11 * r_2 * h6_m5 + e_4 * fs_2_143_286 * r_2 * h8_m5 - e_4 * fs_2_11_11 * r_4 * h6_m5;

        pc_4[k] = e_1 * fs_45_2_15 * h4_m4 - e_2 * fs_255_44_22 * h6_m6 + e_2 * fs_255_22_3 * h6_m4 - e_2 * fs_135_11_15 * r_2 * h4_m4 - e_3 * fs_19_286_2002 * h8_m6 + e_3 * fs_19_143_33 * h8_m4 + e_3 * fs_17_11_22 * r_2 * h6_m6 - e_3 * fs_34_11_3 * r_2 * h6_m4 + e_3 * fs_270_143_15 * r_4 * h4_m4 + e_4 * fs_1_143_2002 * r_2 * h8_m6 - e_4 * fs_2_143_33 * r_2 * h8_m4 - e_4 * fs_1_11_22 * r_4 * h6_m6 + e_4 * fs_2_11_3 * r_4 * h6_m4 - e_4 * fs_12_143_15 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m3, ph4_m2, ph6_m3, ph6_m2, ph8_m8, ph8_m7, ph8_m3, ph8_m2, ab_2, pc_5, pc_6 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_5[k] = - e_1 * fs_225_4_3 * h4_m3 - e_2 * f_255_22 * h6_m3 + e_2 * fs_675_22_3 * r_2 * h4_m3 - e_3 * fs_19_286_6006 * h8_m7 - e_3 * fs_19_286_22 * h8_m3 + e_3 * f_34_11 * r_2 * h6_m3 - e_3 * fs_675_143_3 * r_4 * h4_m3 + e_4 * fs_1_143_6006 * r_2 * h8_m7 + e_4 * fs_1_143_22 * r_2 * h8_m3 - e_4 * f_2_11 * r_4 * h6_m3 + e_4 * fs_30_143_3 * r_6 * h4_m3;

        pc_6[k] = e_0 * fs_195_8_105 * h2_m2 + e_1 * fs_675_28_7 * h4_m2 - e_1 * fs_195_7_105 * r_2 * h2_m2 + e_2 * fs_85_44_6 * h6_m2 - e_2 * fs_2025_154_7 * r_2 * h4_m2 + e_2 * fs_65_7_105 * r_4 * h2_m2 - e_3 * fs_38_143_1001 * h8_m8 + e_3 * fs_19_286_2 * h8_m2 - e_3 * fs_17_33_6 * r_2 * h6_m2 + e_3 * fs_2025_1001_7 * r_4 * h4_m2 - e_3 * fs_260_231_105 * r_6 * h2_m2 + e_4 * fs_4_143_1001 * r_2 * h8_m8 - e_4 * fs_1_143_2 * r_2 * h8_m2 + e_4 * fs_1_33_6 * r_4 * h6_m2 - e_4 * fs_90_1001_7 * r_6 * h4_m2 + e_4 * fs_10_231_105 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph6_p6, ph8_p1, ph8_p2, ph8_p6, ph8_p7, ab_2, pc_7, pc_8 : simd::cache_line_size())
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
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_7[k] = e_0 * fs_195_8_42 * h2_p1 + e_1 * fs_405_28_35 * h4_p1 - e_1 * fs_195_7_42 * r_2 * h2_p1 + e_2 * fs_85_22_6 * h6_p1 - e_2 * fs_1215_154_35 * r_2 * h4_p1 + e_2 * fs_65_7_42 * r_4 * h2_p1 + e_3 * fs_19_286_14 * h8_p1 + e_3 * fs_19_286_10010 * h8_p7 - e_3 * fs_34_33_6 * r_2 * h6_p1 + e_3 * fs_1215_1001_35 * r_4 * h4_p1 - e_3 * fs_260_231_42 * r_6 * h2_p1 - e_4 * fs_1_143_14 * r_2 * h8_p1 - e_4 * fs_1_143_10010 * r_2 * h8_p7 + e_4 * fs_2_33_6 * r_4 * h6_p1 - e_4 * fs_54_1001_35 * r_6 * h4_p1 + e_4 * fs_10_231_42 * r_8 * h2_p1;

        pc_8[k] = e_0 * fs_585_8_7 * h2_p2 - e_1 * fs_45_7_105 * h4_p2 - e_1 * fs_585_7_7 * r_2 * h2_p2 - e_2 * fs_255_44_10 * h6_p2 - e_2 * fs_255_44_22 * h6_p6 + e_2 * fs_270_77_105 * r_2 * h4_p2 + e_2 * fs_195_7_7 * r_4 * h2_p2 - e_3 * fs_19_143_30 * h8_p2 + e_3 * fs_19_143_2002 * h8_p6 + e_3 * fs_17_11_10 * r_2 * h6_p2 + e_3 * fs_17_11_22 * r_2 * h6_p6 - e_3 * fs_540_1001_105 * r_4 * h4_p2 - e_3 * fs_260_77_7 * r_6 * h2_p2 + e_4 * fs_2_143_30 * r_2 * h8_p2 - e_4 * fs_2_143_2002 * r_2 * h8_p6 - e_4 * fs_1_11_10 * r_4 * h6_p2 - e_4 * fs_1_11_22 * r_4 * h6_p6 + e_4 * fs_24_1001_105 * r_6 * h4_p2 + e_4 * fs_10_77_7 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_m3, ph4_p3, ph6_m4, ph6_m3, ph6_p3, ph8_m5, ph8_m4, ph8_m3, ph8_p3, ph8_p5, ab_2, pc_9, pc_10, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];

        pc_9[k] = - e_1 * fs_45_4_3 * h4_p3 + e_2 * f_255_11 * h6_p3 + e_2 * fs_135_22_3 * r_2 * h4_p3 + e_3 * fs_95_286_22 * h8_p3 + e_3 * fs_19_286_4290 * h8_p5 - e_3 * f_68_11 * r_2 * h6_p3 - e_3 * fs_135_143_3 * r_4 * h4_p3 - e_4 * fs_5_143_22 * r_2 * h8_p3 - e_4 * fs_1_143_4290 * r_2 * h8_p5 + e_4 * f_4_11 * r_4 * h6_p3 + e_4 * fs_6_143_3 * r_6 * h4_p3;

        pc_10[k] = e_1 * f_135_1 * h4_m4 - e_2 * fs_255_22_5 * h6_m4 - e_2 * f_810_11 * r_2 * h4_m4 - e_3 * fs_76_143_55 * h8_m4 + e_3 * fs_34_11_5 * r_2 * h6_m4 + e_3 * f_1620_143 * r_4 * h4_m4 + e_4 * fs_8_143_55 * r_2 * h8_m4 - e_4 * fs_2_11_5 * r_4 * h6_m4 - e_4 * f_72_143 * r_6 * h4_m4;

        pc_11[k] = - e_1 * fs_45_4_3 * h4_m3 + e_2 * f_255_11 * h6_m3 + e_2 * fs_135_22_3 * r_2 * h4_m3 - e_3 * fs_19_286_4290 * h8_m5 + e_3 * fs_95_286_22 * h8_m3 - e_3 * f_68_11 * r_2 * h6_m3 - e_3 * fs_135_143_3 * r_4 * h4_m3 + e_4 * fs_1_143_4290 * r_2 * h8_m5 - e_4 * fs_5_143_22 * r_2 * h8_m3 + e_4 * f_4_11 * r_4 * h6_m3 + e_4 * fs_6_143_3 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2, pc_12, pc_13 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_12[k] = e_0 * fs_585_8_7 * h2_m2 - e_1 * fs_45_7_105 * h4_m2 - e_1 * fs_585_7_7 * r_2 * h2_m2 + e_2 * fs_255_44_22 * h6_m6 - e_2 * fs_255_44_10 * h6_m2 + e_2 * fs_270_77_105 * r_2 * h4_m2 + e_2 * fs_195_7_7 * r_4 * h2_m2 - e_3 * fs_19_143_2002 * h8_m6 - e_3 * fs_19_143_30 * h8_m2 - e_3 * fs_17_11_22 * r_2 * h6_m6 + e_3 * fs_17_11_10 * r_2 * h6_m2 - e_3 * fs_540_1001_105 * r_4 * h4_m2 - e_3 * fs_260_77_7 * r_6 * h2_m2 + e_4 * fs_2_143_2002 * r_2 * h8_m6 + e_4 * fs_2_143_30 * r_2 * h8_m2 + e_4 * fs_1_11_22 * r_4 * h6_m6 - e_4 * fs_1_11_10 * r_4 * h6_m2 + e_4 * fs_24_1001_105 * r_6 * h4_m2 + e_4 * fs_10_77_7 * r_8 * h2_m2;

        pc_13[k] = e_0 * fs_195_8_42 * h2_m1 + e_1 * fs_405_28_35 * h4_m1 - e_1 * fs_195_7_42 * r_2 * h2_m1 + e_2 * fs_85_22_6 * h6_m1 - e_2 * fs_1215_154_35 * r_2 * h4_m1 + e_2 * fs_65_7_42 * r_4 * h2_m1 - e_3 * fs_19_286_10010 * h8_m7 + e_3 * fs_19_286_14 * h8_m1 - e_3 * fs_34_33_6 * r_2 * h6_m1 + e_3 * fs_1215_1001_35 * r_4 * h4_m1 - e_3 * fs_260_231_42 * r_6 * h2_m1 + e_4 * fs_1_143_10010 * r_2 * h8_m7 - e_4 * fs_1_143_14 * r_2 * h8_m1 + e_4 * fs_2_33_6 * r_4 * h6_m1 - e_4 * fs_54_1001_35 * r_6 * h4_m1 + e_4 * fs_10_231_42 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ph8_0, ph8_p1, ph8_p5, ph8_p6, ab_2, pc_14, pc_15 : simd::cache_line_size())
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
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_14[k] = e_0 * fs_195_4_7 * h2_0 + e_1 * fs_675_14_7 * h4_0 - e_1 * fs_390_7_7 * r_2 * h2_0 + e_2 * fs_85_11_7 * h6_0 + e_2 * fs_85_44_66 * h6_p6 - e_2 * fs_2025_77_7 * r_2 * h4_0 + e_2 * fs_130_7_7 * r_4 * h2_0 + e_3 * fs_38_143_7 * h8_0 + e_3 * fs_19_286_6006 * h8_p6 - e_3 * fs_68_33_7 * r_2 * h6_0 - e_3 * fs_17_33_66 * r_2 * h6_p6 + e_3 * fs_4050_1001_7 * r_4 * h4_0 - e_3 * fs_520_231_7 * r_6 * h2_0 - e_4 * fs_4_143_7 * r_2 * h8_0 - e_4 * fs_1_143_6006 * r_2 * h8_p6 + e_4 * fs_4_33_7 * r_4 * h6_0 + e_4 * fs_1_33_66 * r_4 * h6_p6 - e_4 * fs_180_1001_7 * r_6 * h4_0 + e_4 * fs_20_231_7 * r_8 * h2_0;

        pc_15[k] = e_0 * fs_195_4_14 * h2_p1 - e_1 * fs_45_28_105 * h4_p1 - e_1 * fs_390_7_14 * r_2 * h2_p1 - e_2 * fs_170_11_2 * h6_p1 - e_2 * fs_85_22_33 * h6_p5 + e_2 * fs_135_154_105 * r_2 * h4_p1 + e_2 * fs_130_7_14 * r_4 * h2_p1 - e_3 * fs_57_286_42 * h8_p1 + e_3 * fs_57_286_858 * h8_p5 + e_3 * fs_136_33_2 * r_2 * h6_p1 + e_3 * fs_34_33_33 * r_2 * h6_p5 - e_3 * fs_135_1001_105 * r_4 * h4_p1 - e_3 * fs_520_231_14 * r_6 * h2_p1 + e_4 * fs_3_143_42 * r_2 * h8_p1 - e_4 * fs_3_143_858 * r_2 * h8_p5 - e_4 * fs_8_33_2 * r_4 * h6_p1 - e_4 * fs_2_33_33 * r_4 * h6_p5 + e_4 * fs_6_1001_105 * r_6 * h4_p1 + e_4 * fs_20_231_14 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_m3, ph4_p2, ph4_p4, ph6_m3, ph6_p2, ph6_p4, ph8_m3, ph8_p2, ph8_p4, ab_2, pc_16, pc_17 : simd::cache_line_size())
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

        pc_16[k] = e_0 * fs_195_8_35 * h2_p2 - e_1 * fs_405_28_21 * h4_p2 + e_1 * fs_45_1_3 * h4_p4 - e_1 * fs_195_7_35 * r_2 * h2_p2 + e_2 * fs_595_44_2 * h6_p2 - e_2 * fs_85_22_15 * h6_p4 + e_2 * fs_1215_154_21 * r_2 * h4_p2 - e_2 * fs_270_11_3 * r_2 * h4_p4 + e_2 * fs_65_7_35 * r_4 * h2_p2 + e_3 * fs_285_286_6 * h8_p2 + e_3 * fs_57_143_165 * h8_p4 - e_3 * fs_119_33_2 * r_2 * h6_p2 + e_3 * fs_34_33_15 * r_2 * h6_p4 - e_3 * fs_1215_1001_21 * r_4 * h4_p2 + e_3 * fs_540_143_3 * r_4 * h4_p4 - e_3 * fs_260_231_35 * r_6 * h2_p2 - e_4 * fs_15_143_6 * r_2 * h8_p2 - e_4 * fs_6_143_165 * r_2 * h8_p4 + e_4 * fs_7_33_2 * r_4 * h6_p2 - e_4 * fs_2_33_15 * r_4 * h6_p4 + e_4 * fs_54_1001_21 * r_6 * h4_p2 - e_4 * fs_24_143_3 * r_6 * h4_p4 + e_4 * fs_10_231_35 * r_8 * h2_p2;

        pc_17[k] = e_1 * f_135_2 * h4_m3 - e_2 * fs_85_22_3 * h6_m3 - e_2 * f_405_11 * r_2 * h4_m3 - e_3 * fs_95_143_66 * h8_m3 + e_3 * fs_34_33_3 * r_2 * h6_m3 + e_3 * f_810_143 * r_4 * h4_m3 + e_4 * fs_10_143_66 * r_2 * h8_m3 - e_4 * fs_2_33_3 * r_4 * h6_m3 - e_4 * f_36_143 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ab_2, pc_18 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];

        pc_18[k] = e_0 * fs_195_8_35 * h2_m2 - e_1 * fs_45_1_3 * h4_m4 - e_1 * fs_405_28_21 * h4_m2 - e_1 * fs_195_7_35 * r_2 * h2_m2 + e_2 * fs_85_22_15 * h6_m4 + e_2 * fs_595_44_2 * h6_m2 + e_2 * fs_270_11_3 * r_2 * h4_m4 + e_2 * fs_1215_154_21 * r_2 * h4_m2 + e_2 * fs_65_7_35 * r_4 * h2_m2 - e_3 * fs_57_143_165 * h8_m4 + e_3 * fs_285_286_6 * h8_m2 - e_3 * fs_34_33_15 * r_2 * h6_m4 - e_3 * fs_119_33_2 * r_2 * h6_m2 - e_3 * fs_540_143_3 * r_4 * h4_m4 - e_3 * fs_1215_1001_21 * r_4 * h4_m2 - e_3 * fs_260_231_35 * r_6 * h2_m2 + e_4 * fs_6_143_165 * r_2 * h8_m4 - e_4 * fs_15_143_6 * r_2 * h8_m2 + e_4 * fs_2_33_15 * r_4 * h6_m4 + e_4 * fs_7_33_2 * r_4 * h6_m2 + e_4 * fs_24_143_3 * r_6 * h4_m4 + e_4 * fs_54_1001_21 * r_6 * h4_m2 + e_4 * fs_10_231_35 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m6, ph6_m5, ph6_m1, ph8_m6, ph8_m5, ph8_m1, ab_2, pc_19, pc_20 : simd::cache_line_size())
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
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];

        pc_19[k] = e_0 * fs_195_4_14 * h2_m1 - e_1 * fs_45_28_105 * h4_m1 - e_1 * fs_390_7_14 * r_2 * h2_m1 + e_2 * fs_85_22_33 * h6_m5 - e_2 * fs_170_11_2 * h6_m1 + e_2 * fs_135_154_105 * r_2 * h4_m1 + e_2 * fs_130_7_14 * r_4 * h2_m1 - e_3 * fs_57_286_858 * h8_m5 - e_3 * fs_57_286_42 * h8_m1 - e_3 * fs_34_33_33 * r_2 * h6_m5 + e_3 * fs_136_33_2 * r_2 * h6_m1 - e_3 * fs_135_1001_105 * r_4 * h4_m1 - e_3 * fs_520_231_14 * r_6 * h2_m1 + e_4 * fs_3_143_858 * r_2 * h8_m5 + e_4 * fs_3_143_42 * r_2 * h8_m1 + e_4 * fs_2_33_33 * r_4 * h6_m5 - e_4 * fs_8_33_2 * r_4 * h6_m1 + e_4 * fs_6_1001_105 * r_6 * h4_m1 + e_4 * fs_20_231_14 * r_8 * h2_m1;

        pc_20[k] = - e_2 * fs_85_44_66 * h6_m6 - e_3 * fs_19_286_6006 * h8_m6 + e_3 * fs_17_33_66 * r_2 * h6_m6 + e_4 * fs_1_143_6006 * r_2 * h8_m6 - e_4 * fs_1_33_66 * r_4 * h6_m6;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ab_2, pc_21 : simd::cache_line_size())
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
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];

        pc_21[k] = - e_0 * fs_195_16_14 * h2_p1 - e_1 * fs_225_28_105 * h4_p1 + e_1 * fs_195_14_14 * r_2 * h2_p1 - e_2 * fs_595_44_2 * h6_p1 + e_2 * fs_85_22_33 * h6_p5 + e_2 * fs_675_154_105 * r_2 * h4_p1 - e_2 * fs_65_14_14 * r_4 * h2_p1 - e_3 * fs_19_143_42 * h8_p1 + e_3 * fs_19_143_858 * h8_p5 + e_3 * fs_119_33_2 * r_2 * h6_p1 - e_3 * fs_34_33_33 * r_2 * h6_p5 - e_3 * fs_675_1001_105 * r_4 * h4_p1 + e_3 * fs_130_231_14 * r_6 * h2_p1 + e_4 * fs_2_143_42 * r_2 * h8_p1 - e_4 * fs_2_143_858 * r_2 * h8_p5 - e_4 * fs_7_33_2 * r_4 * h6_p1 + e_4 * fs_2_33_33 * r_4 * h6_p5 + e_4 * fs_30_1001_105 * r_6 * h4_p1 - e_4 * fs_5_231_14 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ab_2, pc_22 : simd::cache_line_size())
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
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_22[k] = e_0 * fs_585_8_7 * h2_0 + e_1 * fs_225_14_7 * h4_0 - e_1 * fs_45_2_5 * h4_p4 - e_1 * fs_585_7_7 * r_2 * h2_0 - e_2 * fs_255_22_7 * h6_0 - e_2 * f_255_22 * h6_p4 - e_2 * fs_675_77_7 * r_2 * h4_0 + e_2 * fs_135_11_5 * r_2 * h4_p4 + e_2 * fs_195_7_7 * r_4 * h2_0 - e_3 * fs_152_143_7 * h8_0 + e_3 * fs_228_143_11 * h8_p4 + e_3 * fs_34_11_7 * r_2 * h6_0 + e_3 * f_34_11 * r_2 * h6_p4 + e_3 * fs_1350_1001_7 * r_4 * h4_0 - e_3 * fs_270_143_5 * r_4 * h4_p4 - e_3 * fs_260_77_7 * r_6 * h2_0 + e_4 * fs_16_143_7 * r_2 * h8_0 - e_4 * fs_24_143_11 * r_2 * h8_p4 - e_4 * fs_2_11_7 * r_4 * h6_0 - e_4 * f_2_11 * r_4 * h6_p4 - e_4 * fs_60_1001_7 * r_6 * h4_0 + e_4 * fs_12_143_5 * r_6 * h4_p4 + e_4 * fs_10_77_7 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph4_p3, ph6_m2, ph6_p1, ph6_p3, ph8_m2, ph8_p1, ph8_p3, ab_2, pc_23, pc_24 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_23[k] = e_0 * fs_195_16_210 * h2_p1 - e_1 * fs_180_7_7 * h4_p1 + e_1 * f_315_4 * h4_p3 - e_1 * fs_195_14_210 * r_2 * h2_p1 + e_2 * fs_85_44_30 * h6_p1 - e_2 * fs_255_22_3 * h6_p3 + e_2 * fs_1080_77_7 * r_2 * h4_p1 - e_2 * f_945_22 * r_2 * h4_p3 + e_2 * fs_65_14_210 * r_4 * h2_p1 + e_3 * fs_57_143_70 * h8_p1 + e_3 * fs_95_143_66 * h8_p3 - e_3 * fs_17_33_30 * r_2 * h6_p1 + e_3 * fs_34_11_3 * r_2 * h6_p3 - e_3 * fs_2160_1001_7 * r_4 * h4_p1 + e_3 * f_945_143 * r_4 * h4_p3 - e_3 * fs_130_231_210 * r_6 * h2_p1 - e_4 * fs_6_143_70 * r_2 * h8_p1 - e_4 * fs_10_143_66 * r_2 * h8_p3 + e_4 * fs_1_33_30 * r_4 * h6_p1 - e_4 * fs_2_11_3 * r_4 * h6_p3 + e_4 * fs_96_1001_7 * r_6 * h4_p1 - e_4 * f_42_143 * r_6 * h4_p3 + e_4 * fs_5_231_210 * r_8 * h2_p1;

        pc_24[k] = e_0 * fs_195_8_35 * h2_m2 - e_1 * fs_45_14_21 * h4_m2 - e_1 * fs_195_7_35 * r_2 * h2_m2 + e_2 * fs_85_11_2 * h6_m2 + e_2 * fs_135_77_21 * r_2 * h4_m2 + e_2 * fs_65_7_35 * r_4 * h2_m2 - e_3 * fs_380_143_6 * h8_m2 - e_3 * fs_68_33_2 * r_2 * h6_m2 - e_3 * fs_270_1001_21 * r_4 * h4_m2 - e_3 * fs_260_231_35 * r_6 * h2_m2 + e_4 * fs_40_143_6 * r_2 * h8_m2 + e_4 * fs_4_33_2 * r_4 * h6_m2 + e_4 * fs_12_1001_21 * r_6 * h4_m2 + e_4 * fs_10_231_35 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_25, pc_26, pc_27 : simd::cache_line_size())
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

        pc_25[k] = e_0 * fs_195_16_210 * h2_m1 - e_1 * f_315_4 * h4_m3 - e_1 * fs_180_7_7 * h4_m1 - e_1 * fs_195_14_210 * r_2 * h2_m1 + e_2 * fs_255_22_3 * h6_m3 + e_2 * fs_85_44_30 * h6_m1 + e_2 * f_945_22 * r_2 * h4_m3 + e_2 * fs_1080_77_7 * r_2 * h4_m1 + e_2 * fs_65_14_210 * r_4 * h2_m1 - e_3 * fs_95_143_66 * h8_m3 + e_3 * fs_57_143_70 * h8_m1 - e_3 * fs_34_11_3 * r_2 * h6_m3 - e_3 * fs_17_33_30 * r_2 * h6_m1 - e_3 * f_945_143 * r_4 * h4_m3 - e_3 * fs_2160_1001_7 * r_4 * h4_m1 - e_3 * fs_130_231_210 * r_6 * h2_m1 + e_4 * fs_10_143_66 * r_2 * h8_m3 - e_4 * fs_6_143_70 * r_2 * h8_m1 + e_4 * fs_2_11_3 * r_4 * h6_m3 + e_4 * fs_1_33_30 * r_4 * h6_m1 + e_4 * f_42_143 * r_6 * h4_m3 + e_4 * fs_96_1001_7 * r_6 * h4_m1 + e_4 * fs_5_231_210 * r_8 * h2_m1;

        pc_26[k] = e_1 * fs_45_2_5 * h4_m4 + e_2 * f_255_22 * h6_m4 - e_2 * fs_135_11_5 * r_2 * h4_m4 - e_3 * fs_228_143_11 * h8_m4 - e_3 * f_34_11 * r_2 * h6_m4 + e_3 * fs_270_143_5 * r_4 * h4_m4 + e_4 * fs_24_143_11 * r_2 * h8_m4 + e_4 * f_2_11 * r_4 * h6_m4 - e_4 * fs_12_143_5 * r_6 * h4_m4;

        pc_27[k] = e_0 * fs_195_16_14 * h2_m1 + e_1 * fs_225_28_105 * h4_m1 - e_1 * fs_195_14_14 * r_2 * h2_m1 - e_2 * fs_85_22_33 * h6_m5 + e_2 * fs_595_44_2 * h6_m1 - e_2 * fs_675_154_105 * r_2 * h4_m1 + e_2 * fs_65_14_14 * r_4 * h2_m1 - e_3 * fs_19_143_858 * h8_m5 + e_3 * fs_19_143_42 * h8_m1 + e_3 * fs_34_33_33 * r_2 * h6_m5 - e_3 * fs_119_33_2 * r_2 * h6_m1 + e_3 * fs_675_1001_105 * r_4 * h4_m1 - e_3 * fs_130_231_14 * r_6 * h2_m1 + e_4 * fs_2_143_858 * r_2 * h8_m5 - e_4 * fs_2_143_42 * r_2 * h8_m1 - e_4 * fs_2_33_33 * r_4 * h6_m5 + e_4 * fs_7_33_2 * r_4 * h6_m1 - e_4 * fs_30_1001_105 * r_6 * h4_m1 + e_4 * fs_5_231_14 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ab_2, pc_28 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];

        pc_28[k] = e_0 * fs_195_16_2 * h2_p2 + e_1 * fs_675_56_30 * h4_p2 + e_1 * fs_45_28_210 * h4_p4 - e_1 * fs_195_14_2 * r_2 * h2_p2 + e_2 * fs_85_22_35 * h6_p2 + e_2 * fs_85_22_42 * h6_p4 - e_2 * fs_2025_308_30 * r_2 * h4_p2 - e_2 * fs_135_154_210 * r_2 * h4_p4 + e_2 * fs_65_14_2 * r_4 * h2_p2 + e_3 * fs_19_143_105 * h8_p2 + e_3 * fs_19_143_462 * h8_p4 - e_3 * fs_34_33_35 * r_2 * h6_p2 - e_3 * fs_34_33_42 * r_2 * h6_p4 + e_3 * fs_2025_2002_30 * r_4 * h4_p2 + e_3 * fs_135_1001_210 * r_4 * h4_p4 - e_3 * fs_130_231_2 * r_6 * h2_p2 - e_4 * fs_2_143_105 * r_2 * h8_p2 - e_4 * fs_2_143_462 * r_2 * h8_p4 + e_4 * fs_2_33_35 * r_4 * h6_p2 + e_4 * fs_2_33_42 * r_4 * h6_p4 - e_4 * fs_45_1001_30 * r_6 * h4_p2 - e_4 * fs_6_1001_210 * r_6 * h4_p4 + e_4 * fs_5_231_2 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph8_p1, ph8_p3, ab_2, pc_29 : simd::cache_line_size())
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
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_29[k] = - e_0 * fs_195_4_3 * h2_p1 - e_1 * fs_1125_56_10 * h4_p1 - e_1 * fs_495_56_70 * h4_p3 + e_1 * fs_390_7_3 * r_2 * h2_p1 + e_2 * fs_85_22_21 * h6_p1 + e_2 * fs_3375_308_10 * r_2 * h4_p1 + e_2 * fs_135_28_70 * r_2 * h4_p3 - e_2 * fs_130_7_3 * r_4 * h2_p1 + e_3 * f_399_143 * h8_p1 + e_3 * fs_19_143_1155 * h8_p3 - e_3 * fs_34_33_21 * r_2 * h6_p1 - e_3 * fs_3375_2002_10 * r_4 * h4_p1 - e_3 * fs_135_182_70 * r_4 * h4_p3 + e_3 * fs_520_231_3 * r_6 * h2_p1 - e_4 * f_42_143 * r_2 * h8_p1 - e_4 * fs_2_143_1155 * r_2 * h8_p3 + e_4 * fs_2_33_21 * r_4 * h6_p1 + e_4 * fs_75_1001_10 * r_6 * h4_p1 + e_4 * fs_3_91_70 * r_6 * h4_p3 - e_4 * fs_20_231_3 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph2_0, ph2_p2, ph4_m1, ph4_0, ph4_p2, ph6_m1, ph6_p2, ph8_m1, ph8_0, ph8_p2, ab_2, pc_30, pc_31 : simd::cache_line_size())
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
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];

        pc_30[k] = e_0 * fs_585_8_10 * h2_0 + e_0 * fs_195_16_30 * h2_p2 - e_1 * fs_495_28_10 * h4_0 + e_1 * fs_1935_56_2 * h4_p2 - e_1 * fs_585_7_10 * r_2 * h2_0 - e_1 * fs_195_14_30 * r_2 * h2_p2 - e_2 * fs_85_22_21 * h6_p2 + e_2 * fs_135_14_10 * r_2 * h4_0 - e_2 * fs_5805_308_2 * r_2 * h4_p2 + e_2 * fs_195_7_10 * r_4 * h2_0 + e_2 * fs_65_14_30 * r_4 * h2_p2 + e_3 * fs_266_143_10 * h8_0 + e_3 * fs_285_143_7 * h8_p2 + e_3 * fs_34_33_21 * r_2 * h6_p2 - e_3 * fs_135_91_10 * r_4 * h4_0 + e_3 * fs_5805_2002_2 * r_4 * h4_p2 - e_3 * fs_260_77_10 * r_6 * h2_0 - e_3 * fs_130_231_30 * r_6 * h2_p2 - e_4 * fs_28_143_10 * r_2 * h8_0 - e_4 * fs_30_143_7 * r_2 * h8_p2 - e_4 * fs_2_33_21 * r_4 * h6_p2 + e_4 * fs_6_91_10 * r_6 * h4_0 - e_4 * fs_129_1001_2 * r_6 * h4_p2 + e_4 * fs_10_77_10 * r_8 * h2_0 + e_4 * fs_5_231_30 * r_8 * h2_p2;

        pc_31[k] = e_0 * fs_195_2_5 * h2_m1 - e_1 * fs_855_28_6 * h4_m1 - e_1 * fs_780_7_5 * r_2 * h2_m1 + e_2 * fs_85_22_35 * h6_m1 + e_2 * fs_2565_154_6 * r_2 * h4_m1 + e_2 * fs_260_7_5 * r_4 * h2_m1 - e_3 * fs_266_143_15 * h8_m1 - e_3 * fs_34_33_35 * r_2 * h6_m1 - e_3 * fs_2565_1001_6 * r_4 * h4_m1 - e_3 * fs_1040_231_5 * r_6 * h2_m1 + e_4 * fs_28_143_15 * r_2 * h8_m1 + e_4 * fs_2_33_35 * r_4 * h6_m1 + e_4 * fs_114_1001_6 * r_6 * h4_m1 + e_4 * fs_40_231_5 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_32, pc_33 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_32[k] = - e_0 * fs_195_16_30 * h2_m2 - e_1 * fs_1935_56_2 * h4_m2 + e_1 * fs_195_14_30 * r_2 * h2_m2 + e_2 * fs_85_22_21 * h6_m2 + e_2 * fs_5805_308_2 * r_2 * h4_m2 - e_2 * fs_65_14_30 * r_4 * h2_m2 - e_3 * fs_285_143_7 * h8_m2 - e_3 * fs_34_33_21 * r_2 * h6_m2 - e_3 * fs_5805_2002_2 * r_4 * h4_m2 + e_3 * fs_130_231_30 * r_6 * h2_m2 + e_4 * fs_30_143_7 * r_2 * h8_m2 + e_4 * fs_2_33_21 * r_4 * h6_m2 + e_4 * fs_129_1001_2 * r_6 * h4_m2 - e_4 * fs_5_231_30 * r_8 * h2_m2;

        pc_33[k] = e_0 * fs_195_4_3 * h2_m1 + e_1 * fs_495_56_70 * h4_m3 + e_1 * fs_1125_56_10 * h4_m1 - e_1 * fs_390_7_3 * r_2 * h2_m1 - e_2 * fs_85_22_21 * h6_m1 - e_2 * fs_135_28_70 * r_2 * h4_m3 - e_2 * fs_3375_308_10 * r_2 * h4_m1 + e_2 * fs_130_7_3 * r_4 * h2_m1 - e_3 * fs_19_143_1155 * h8_m3 - e_3 * f_399_143 * h8_m1 + e_3 * fs_34_33_21 * r_2 * h6_m1 + e_3 * fs_135_182_70 * r_4 * h4_m3 + e_3 * fs_3375_2002_10 * r_4 * h4_m1 - e_3 * fs_520_231_3 * r_6 * h2_m1 + e_4 * fs_2_143_1155 * r_2 * h8_m3 + e_4 * f_42_143 * r_2 * h8_m1 - e_4 * fs_2_33_21 * r_4 * h6_m1 - e_4 * fs_3_91_70 * r_6 * h4_m3 - e_4 * fs_75_1001_10 * r_6 * h4_m1 + e_4 * fs_20_231_3 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m4, ph6_m3, ph6_m2, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_34, pc_35, pc_36 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_34[k] = - e_0 * fs_195_16_2 * h2_m2 - e_1 * fs_45_28_210 * h4_m4 - e_1 * fs_675_56_30 * h4_m2 + e_1 * fs_195_14_2 * r_2 * h2_m2 - e_2 * fs_85_22_42 * h6_m4 - e_2 * fs_85_22_35 * h6_m2 + e_2 * fs_135_154_210 * r_2 * h4_m4 + e_2 * fs_2025_308_30 * r_2 * h4_m2 - e_2 * fs_65_14_2 * r_4 * h2_m2 - e_3 * fs_19_143_462 * h8_m4 - e_3 * fs_19_143_105 * h8_m2 + e_3 * fs_34_33_42 * r_2 * h6_m4 + e_3 * fs_34_33_35 * r_2 * h6_m2 - e_3 * fs_135_1001_210 * r_4 * h4_m4 - e_3 * fs_2025_2002_30 * r_4 * h4_m2 + e_3 * fs_130_231_2 * r_6 * h2_m2 + e_4 * fs_2_143_462 * r_2 * h8_m4 + e_4 * fs_2_143_105 * r_2 * h8_m2 - e_4 * fs_2_33_42 * r_4 * h6_m4 - e_4 * fs_2_33_35 * r_4 * h6_m2 + e_4 * fs_6_1001_210 * r_6 * h4_m4 + e_4 * fs_45_1001_30 * r_6 * h4_m2 - e_4 * fs_5_231_2 * r_8 * h2_m2;

        pc_35[k] = - e_1 * fs_675_28_7 * h4_m3 - e_2 * fs_85_11_21 * h6_m3 + e_2 * fs_2025_154_7 * r_2 * h4_m3 - e_3 * fs_19_143_462 * h8_m3 + e_3 * fs_68_33_21 * r_2 * h6_m3 - e_3 * fs_2025_1001_7 * r_4 * h4_m3 + e_4 * fs_2_143_462 * r_2 * h8_m3 - e_4 * fs_4_33_21 * r_4 * h6_m3 + e_4 * fs_90_1001_7 * r_6 * h4_m3;

        pc_36[k] = e_0 * fs_195_8_5 * h2_m2 + e_1 * fs_450_7_3 * h4_m2 - e_1 * fs_195_7_5 * r_2 * h2_m2 - e_2 * fs_85_22_14 * h6_m2 - e_2 * fs_2700_77_3 * r_2 * h4_m2 + e_2 * fs_65_7_5 * r_4 * h2_m2 - e_3 * fs_114_143_42 * h8_m2 + e_3 * fs_34_33_14 * r_2 * h6_m2 + e_3 * fs_5400_1001_3 * r_4 * h4_m2 - e_3 * fs_260_231_5 * r_6 * h2_m2 + e_4 * fs_12_143_42 * r_2 * h8_m2 - e_4 * fs_2_33_14 * r_4 * h6_m2 - e_4 * fs_240_1001_3 * r_6 * h4_m2 + e_4 * fs_10_231_5 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph2_0, ph2_p1, ph4_m1, ph4_0, ph4_p1, ph6_m1, ph6_0, ph6_p1, ph8_m1, ph8_0, ph8_p1, ab_2, pc_37, pc_38, pc_39 : simd::cache_line_size())
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
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];

        pc_37[k] = - e_0 * fs_975_8_2 * h2_m1 - e_1 * fs_45_28_15 * h4_m1 + e_1 * fs_975_7_2 * r_2 * h2_m1 + e_2 * fs_85_22_14 * h6_m1 + e_2 * fs_135_154_15 * r_2 * h4_m1 - e_2 * fs_325_7_2 * r_4 * h2_m1 - e_3 * fs_399_143_6 * h8_m1 - e_3 * fs_34_33_14 * r_2 * h6_m1 - e_3 * fs_135_1001_15 * r_4 * h4_m1 + e_3 * fs_1300_231_2 * r_6 * h2_m1 + e_4 * fs_42_143_6 * r_2 * h8_m1 + e_4 * fs_2_33_14 * r_4 * h6_m1 + e_4 * fs_6_1001_15 * r_6 * h4_m1 - e_4 * fs_50_231_2 * r_8 * h2_m1;

        pc_38[k] = e_0 * f_975_4 * h2_0 - e_1 * f_675_7 * h4_0 - e_1 * f_1950_7 * r_2 * h2_0 + e_2 * f_595_22 * h6_0 + e_2 * f_4050_77 * r_2 * h4_0 + e_2 * f_650_7 * r_4 * h2_0 - e_3 * f_1064_143 * h8_0 - e_3 * f_238_33 * r_2 * h6_0 - e_3 * f_8100_1001 * r_4 * h4_0 - e_3 * f_2600_231 * r_6 * h2_0 + e_4 * f_112_143 * r_2 * h8_0 + e_4 * f_14_33 * r_4 * h6_0 + e_4 * f_360_1001 * r_6 * h4_0 + e_4 * f_100_231 * r_8 * h2_0;

        pc_39[k] = - e_0 * fs_975_8_2 * h2_p1 - e_1 * fs_45_28_15 * h4_p1 + e_1 * fs_975_7_2 * r_2 * h2_p1 + e_2 * fs_85_22_14 * h6_p1 + e_2 * fs_135_154_15 * r_2 * h4_p1 - e_2 * fs_325_7_2 * r_4 * h2_p1 - e_3 * fs_399_143_6 * h8_p1 - e_3 * fs_34_33_14 * r_2 * h6_p1 - e_3 * fs_135_1001_15 * r_4 * h4_p1 + e_3 * fs_1300_231_2 * r_6 * h2_p1 + e_4 * fs_42_143_6 * r_2 * h8_p1 + e_4 * fs_2_33_14 * r_4 * h6_p1 + e_4 * fs_6_1001_15 * r_6 * h4_p1 - e_4 * fs_50_231_2 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ab_2, pc_40, pc_41 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_40[k] = e_0 * fs_195_8_5 * h2_p2 + e_1 * fs_450_7_3 * h4_p2 - e_1 * fs_195_7_5 * r_2 * h2_p2 - e_2 * fs_85_22_14 * h6_p2 - e_2 * fs_2700_77_3 * r_2 * h4_p2 + e_2 * fs_65_7_5 * r_4 * h2_p2 - e_3 * fs_114_143_42 * h8_p2 + e_3 * fs_34_33_14 * r_2 * h6_p2 + e_3 * fs_5400_1001_3 * r_4 * h4_p2 - e_3 * fs_260_231_5 * r_6 * h2_p2 + e_4 * fs_12_143_42 * r_2 * h8_p2 - e_4 * fs_2_33_14 * r_4 * h6_p2 - e_4 * fs_240_1001_3 * r_6 * h4_p2 + e_4 * fs_10_231_5 * r_8 * h2_p2;

        pc_41[k] = - e_1 * fs_675_28_7 * h4_p3 - e_2 * fs_85_11_21 * h6_p3 + e_2 * fs_2025_154_7 * r_2 * h4_p3 - e_3 * fs_19_143_462 * h8_p3 + e_3 * fs_68_33_21 * r_2 * h6_p3 - e_3 * fs_2025_1001_7 * r_4 * h4_p3 + e_4 * fs_2_143_462 * r_2 * h8_p3 - e_4 * fs_4_33_21 * r_4 * h6_p3 + e_4 * fs_90_1001_7 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ab_2, pc_42 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];

        pc_42[k] = e_0 * fs_195_16_2 * h2_m2 - e_1 * fs_45_28_210 * h4_m4 + e_1 * fs_675_56_30 * h4_m2 - e_1 * fs_195_14_2 * r_2 * h2_m2 - e_2 * fs_85_22_42 * h6_m4 + e_2 * fs_85_22_35 * h6_m2 + e_2 * fs_135_154_210 * r_2 * h4_m4 - e_2 * fs_2025_308_30 * r_2 * h4_m2 + e_2 * fs_65_14_2 * r_4 * h2_m2 - e_3 * fs_19_143_462 * h8_m4 + e_3 * fs_19_143_105 * h8_m2 + e_3 * fs_34_33_42 * r_2 * h6_m4 - e_3 * fs_34_33_35 * r_2 * h6_m2 - e_3 * fs_135_1001_210 * r_4 * h4_m4 + e_3 * fs_2025_2002_30 * r_4 * h4_m2 - e_3 * fs_130_231_2 * r_6 * h2_m2 + e_4 * fs_2_143_462 * r_2 * h8_m4 - e_4 * fs_2_143_105 * r_2 * h8_m2 - e_4 * fs_2_33_42 * r_4 * h6_m4 + e_4 * fs_2_33_35 * r_4 * h6_m2 + e_4 * fs_6_1001_210 * r_6 * h4_m4 - e_4 * fs_45_1001_30 * r_6 * h4_m2 + e_4 * fs_5_231_2 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_43, pc_44 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_43[k] = - e_0 * fs_195_4_3 * h2_m1 + e_1 * fs_495_56_70 * h4_m3 - e_1 * fs_1125_56_10 * h4_m1 + e_1 * fs_390_7_3 * r_2 * h2_m1 + e_2 * fs_85_22_21 * h6_m1 - e_2 * fs_135_28_70 * r_2 * h4_m3 + e_2 * fs_3375_308_10 * r_2 * h4_m1 - e_2 * fs_130_7_3 * r_4 * h2_m1 - e_3 * fs_19_143_1155 * h8_m3 + e_3 * f_399_143 * h8_m1 - e_3 * fs_34_33_21 * r_2 * h6_m1 + e_3 * fs_135_182_70 * r_4 * h4_m3 - e_3 * fs_3375_2002_10 * r_4 * h4_m1 + e_3 * fs_520_231_3 * r_6 * h2_m1 + e_4 * fs_2_143_1155 * r_2 * h8_m3 - e_4 * f_42_143 * r_2 * h8_m1 + e_4 * fs_2_33_21 * r_4 * h6_m1 - e_4 * fs_3_91_70 * r_6 * h4_m3 + e_4 * fs_75_1001_10 * r_6 * h4_m1 - e_4 * fs_20_231_3 * r_8 * h2_m1;

        pc_44[k] = - e_0 * fs_195_16_30 * h2_m2 - e_1 * fs_1935_56_2 * h4_m2 + e_1 * fs_195_14_30 * r_2 * h2_m2 + e_2 * fs_85_22_21 * h6_m2 + e_2 * fs_5805_308_2 * r_2 * h4_m2 - e_2 * fs_65_14_30 * r_4 * h2_m2 - e_3 * fs_285_143_7 * h8_m2 - e_3 * fs_34_33_21 * r_2 * h6_m2 - e_3 * fs_5805_2002_2 * r_4 * h4_m2 + e_3 * fs_130_231_30 * r_6 * h2_m2 + e_4 * fs_30_143_7 * r_2 * h8_m2 + e_4 * fs_2_33_21 * r_4 * h6_m2 + e_4 * fs_129_1001_2 * r_6 * h4_m2 - e_4 * fs_5_231_30 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph8_0, ph8_p1, ph8_p2, ab_2, pc_45, pc_46 : simd::cache_line_size())
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
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_45[k] = e_0 * fs_195_2_5 * h2_p1 - e_1 * fs_855_28_6 * h4_p1 - e_1 * fs_780_7_5 * r_2 * h2_p1 + e_2 * fs_85_22_35 * h6_p1 + e_2 * fs_2565_154_6 * r_2 * h4_p1 + e_2 * fs_260_7_5 * r_4 * h2_p1 - e_3 * fs_266_143_15 * h8_p1 - e_3 * fs_34_33_35 * r_2 * h6_p1 - e_3 * fs_2565_1001_6 * r_4 * h4_p1 - e_3 * fs_1040_231_5 * r_6 * h2_p1 + e_4 * fs_28_143_15 * r_2 * h8_p1 + e_4 * fs_2_33_35 * r_4 * h6_p1 + e_4 * fs_114_1001_6 * r_6 * h4_p1 + e_4 * fs_40_231_5 * r_8 * h2_p1;

        pc_46[k] = e_0 * fs_585_8_10 * h2_0 - e_0 * fs_195_16_30 * h2_p2 - e_1 * fs_495_28_10 * h4_0 - e_1 * fs_1935_56_2 * h4_p2 - e_1 * fs_585_7_10 * r_2 * h2_0 + e_1 * fs_195_14_30 * r_2 * h2_p2 + e_2 * fs_85_22_21 * h6_p2 + e_2 * fs_135_14_10 * r_2 * h4_0 + e_2 * fs_5805_308_2 * r_2 * h4_p2 + e_2 * fs_195_7_10 * r_4 * h2_0 - e_2 * fs_65_14_30 * r_4 * h2_p2 + e_3 * fs_266_143_10 * h8_0 - e_3 * fs_285_143_7 * h8_p2 - e_3 * fs_34_33_21 * r_2 * h6_p2 - e_3 * fs_135_91_10 * r_4 * h4_0 - e_3 * fs_5805_2002_2 * r_4 * h4_p2 - e_3 * fs_260_77_10 * r_6 * h2_0 + e_3 * fs_130_231_30 * r_6 * h2_p2 - e_4 * fs_28_143_10 * r_2 * h8_0 + e_4 * fs_30_143_7 * r_2 * h8_p2 + e_4 * fs_2_33_21 * r_4 * h6_p2 + e_4 * fs_6_91_10 * r_6 * h4_0 + e_4 * fs_129_1001_2 * r_6 * h4_p2 + e_4 * fs_10_77_10 * r_8 * h2_0 - e_4 * fs_5_231_30 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph8_p1, ph8_p3, ab_2, pc_47 : simd::cache_line_size())
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
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_47[k] = - e_0 * fs_195_4_3 * h2_p1 - e_1 * fs_1125_56_10 * h4_p1 + e_1 * fs_495_56_70 * h4_p3 + e_1 * fs_390_7_3 * r_2 * h2_p1 + e_2 * fs_85_22_21 * h6_p1 + e_2 * fs_3375_308_10 * r_2 * h4_p1 - e_2 * fs_135_28_70 * r_2 * h4_p3 - e_2 * fs_130_7_3 * r_4 * h2_p1 + e_3 * f_399_143 * h8_p1 - e_3 * fs_19_143_1155 * h8_p3 - e_3 * fs_34_33_21 * r_2 * h6_p1 - e_3 * fs_3375_2002_10 * r_4 * h4_p1 + e_3 * fs_135_182_70 * r_4 * h4_p3 + e_3 * fs_520_231_3 * r_6 * h2_p1 - e_4 * f_42_143 * r_2 * h8_p1 + e_4 * fs_2_143_1155 * r_2 * h8_p3 + e_4 * fs_2_33_21 * r_4 * h6_p1 + e_4 * fs_75_1001_10 * r_6 * h4_p1 - e_4 * fs_3_91_70 * r_6 * h4_p3 - e_4 * fs_20_231_3 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ab_2, pc_48 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];

        pc_48[k] = e_0 * fs_195_16_2 * h2_p2 + e_1 * fs_675_56_30 * h4_p2 - e_1 * fs_45_28_210 * h4_p4 - e_1 * fs_195_14_2 * r_2 * h2_p2 + e_2 * fs_85_22_35 * h6_p2 - e_2 * fs_85_22_42 * h6_p4 - e_2 * fs_2025_308_30 * r_2 * h4_p2 + e_2 * fs_135_154_210 * r_2 * h4_p4 + e_2 * fs_65_14_2 * r_4 * h2_p2 + e_3 * fs_19_143_105 * h8_p2 - e_3 * fs_19_143_462 * h8_p4 - e_3 * fs_34_33_35 * r_2 * h6_p2 + e_3 * fs_34_33_42 * r_2 * h6_p4 + e_3 * fs_2025_2002_30 * r_4 * h4_p2 - e_3 * fs_135_1001_210 * r_4 * h4_p4 - e_3 * fs_130_231_2 * r_6 * h2_p2 - e_4 * fs_2_143_105 * r_2 * h8_p2 + e_4 * fs_2_143_462 * r_2 * h8_p4 + e_4 * fs_2_33_35 * r_4 * h6_p2 - e_4 * fs_2_33_42 * r_4 * h6_p4 - e_4 * fs_45_1001_30 * r_6 * h4_p2 + e_4 * fs_6_1001_210 * r_6 * h4_p4 + e_4 * fs_5_231_2 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_49, pc_50, pc_51 : simd::cache_line_size())
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

        pc_49[k] = - e_0 * fs_195_16_14 * h2_m1 - e_1 * fs_225_28_105 * h4_m1 + e_1 * fs_195_14_14 * r_2 * h2_m1 - e_2 * fs_85_22_33 * h6_m5 - e_2 * fs_595_44_2 * h6_m1 + e_2 * fs_675_154_105 * r_2 * h4_m1 - e_2 * fs_65_14_14 * r_4 * h2_m1 - e_3 * fs_19_143_858 * h8_m5 - e_3 * fs_19_143_42 * h8_m1 + e_3 * fs_34_33_33 * r_2 * h6_m5 + e_3 * fs_119_33_2 * r_2 * h6_m1 - e_3 * fs_675_1001_105 * r_4 * h4_m1 + e_3 * fs_130_231_14 * r_6 * h2_m1 + e_4 * fs_2_143_858 * r_2 * h8_m5 + e_4 * fs_2_143_42 * r_2 * h8_m1 - e_4 * fs_2_33_33 * r_4 * h6_m5 - e_4 * fs_7_33_2 * r_4 * h6_m1 + e_4 * fs_30_1001_105 * r_6 * h4_m1 - e_4 * fs_5_231_14 * r_8 * h2_m1;

        pc_50[k] = e_1 * fs_45_2_5 * h4_m4 + e_2 * f_255_22 * h6_m4 - e_2 * fs_135_11_5 * r_2 * h4_m4 - e_3 * fs_228_143_11 * h8_m4 - e_3 * f_34_11 * r_2 * h6_m4 + e_3 * fs_270_143_5 * r_4 * h4_m4 + e_4 * fs_24_143_11 * r_2 * h8_m4 + e_4 * f_2_11 * r_4 * h6_m4 - e_4 * fs_12_143_5 * r_6 * h4_m4;

        pc_51[k] = - e_0 * fs_195_16_210 * h2_m1 - e_1 * f_315_4 * h4_m3 + e_1 * fs_180_7_7 * h4_m1 + e_1 * fs_195_14_210 * r_2 * h2_m1 + e_2 * fs_255_22_3 * h6_m3 - e_2 * fs_85_44_30 * h6_m1 + e_2 * f_945_22 * r_2 * h4_m3 - e_2 * fs_1080_77_7 * r_2 * h4_m1 - e_2 * fs_65_14_210 * r_4 * h2_m1 - e_3 * fs_95_143_66 * h8_m3 - e_3 * fs_57_143_70 * h8_m1 - e_3 * fs_34_11_3 * r_2 * h6_m3 + e_3 * fs_17_33_30 * r_2 * h6_m1 - e_3 * f_945_143 * r_4 * h4_m3 + e_3 * fs_2160_1001_7 * r_4 * h4_m1 + e_3 * fs_130_231_210 * r_6 * h2_m1 + e_4 * fs_10_143_66 * r_2 * h8_m3 + e_4 * fs_6_143_70 * r_2 * h8_m1 + e_4 * fs_2_11_3 * r_4 * h6_m3 - e_4 * fs_1_33_30 * r_4 * h6_m1 + e_4 * f_42_143 * r_6 * h4_m3 - e_4 * fs_96_1001_7 * r_6 * h4_m1 - e_4 * fs_5_231_210 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph4_p3, ph6_p1, ph6_p2, ph6_p3, ph8_p1, ph8_p2, ph8_p3, ab_2, pc_52, pc_53 : simd::cache_line_size())
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
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_52[k] = e_0 * fs_195_8_35 * h2_p2 - e_1 * fs_45_14_21 * h4_p2 - e_1 * fs_195_7_35 * r_2 * h2_p2 + e_2 * fs_85_11_2 * h6_p2 + e_2 * fs_135_77_21 * r_2 * h4_p2 + e_2 * fs_65_7_35 * r_4 * h2_p2 - e_3 * fs_380_143_6 * h8_p2 - e_3 * fs_68_33_2 * r_2 * h6_p2 - e_3 * fs_270_1001_21 * r_4 * h4_p2 - e_3 * fs_260_231_35 * r_6 * h2_p2 + e_4 * fs_40_143_6 * r_2 * h8_p2 + e_4 * fs_4_33_2 * r_4 * h6_p2 + e_4 * fs_12_1001_21 * r_6 * h4_p2 + e_4 * fs_10_231_35 * r_8 * h2_p2;

        pc_53[k] = e_0 * fs_195_16_210 * h2_p1 - e_1 * fs_180_7_7 * h4_p1 - e_1 * f_315_4 * h4_p3 - e_1 * fs_195_14_210 * r_2 * h2_p1 + e_2 * fs_85_44_30 * h6_p1 + e_2 * fs_255_22_3 * h6_p3 + e_2 * fs_1080_77_7 * r_2 * h4_p1 + e_2 * f_945_22 * r_2 * h4_p3 + e_2 * fs_65_14_210 * r_4 * h2_p1 + e_3 * fs_57_143_70 * h8_p1 - e_3 * fs_95_143_66 * h8_p3 - e_3 * fs_17_33_30 * r_2 * h6_p1 - e_3 * fs_34_11_3 * r_2 * h6_p3 - e_3 * fs_2160_1001_7 * r_4 * h4_p1 - e_3 * f_945_143 * r_4 * h4_p3 - e_3 * fs_130_231_210 * r_6 * h2_p1 - e_4 * fs_6_143_70 * r_2 * h8_p1 + e_4 * fs_10_143_66 * r_2 * h8_p3 + e_4 * fs_1_33_30 * r_4 * h6_p1 + e_4 * fs_2_11_3 * r_4 * h6_p3 + e_4 * fs_96_1001_7 * r_6 * h4_p1 + e_4 * f_42_143 * r_6 * h4_p3 + e_4 * fs_5_231_210 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ab_2, pc_54 : simd::cache_line_size())
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
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_54[k] = e_0 * fs_585_8_7 * h2_0 + e_1 * fs_225_14_7 * h4_0 + e_1 * fs_45_2_5 * h4_p4 - e_1 * fs_585_7_7 * r_2 * h2_0 - e_2 * fs_255_22_7 * h6_0 + e_2 * f_255_22 * h6_p4 - e_2 * fs_675_77_7 * r_2 * h4_0 - e_2 * fs_135_11_5 * r_2 * h4_p4 + e_2 * fs_195_7_7 * r_4 * h2_0 - e_3 * fs_152_143_7 * h8_0 - e_3 * fs_228_143_11 * h8_p4 + e_3 * fs_34_11_7 * r_2 * h6_0 - e_3 * f_34_11 * r_2 * h6_p4 + e_3 * fs_1350_1001_7 * r_4 * h4_0 + e_3 * fs_270_143_5 * r_4 * h4_p4 - e_3 * fs_260_77_7 * r_6 * h2_0 + e_4 * fs_16_143_7 * r_2 * h8_0 + e_4 * fs_24_143_11 * r_2 * h8_p4 - e_4 * fs_2_11_7 * r_4 * h6_0 + e_4 * f_2_11 * r_4 * h6_p4 - e_4 * fs_60_1001_7 * r_6 * h4_0 - e_4 * fs_12_143_5 * r_6 * h4_p4 + e_4 * fs_10_77_7 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_m6, ph6_p1, ph6_p5, ph8_m6, ph8_p1, ph8_p5, ab_2, pc_55, pc_56 : simd::cache_line_size())
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
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];

        pc_55[k] = - e_0 * fs_195_16_14 * h2_p1 - e_1 * fs_225_28_105 * h4_p1 + e_1 * fs_195_14_14 * r_2 * h2_p1 - e_2 * fs_595_44_2 * h6_p1 - e_2 * fs_85_22_33 * h6_p5 + e_2 * fs_675_154_105 * r_2 * h4_p1 - e_2 * fs_65_14_14 * r_4 * h2_p1 - e_3 * fs_19_143_42 * h8_p1 - e_3 * fs_19_143_858 * h8_p5 + e_3 * fs_119_33_2 * r_2 * h6_p1 + e_3 * fs_34_33_33 * r_2 * h6_p5 - e_3 * fs_675_1001_105 * r_4 * h4_p1 + e_3 * fs_130_231_14 * r_6 * h2_p1 + e_4 * fs_2_143_42 * r_2 * h8_p1 + e_4 * fs_2_143_858 * r_2 * h8_p5 - e_4 * fs_7_33_2 * r_4 * h6_p1 - e_4 * fs_2_33_33 * r_4 * h6_p5 + e_4 * fs_30_1001_105 * r_6 * h4_p1 - e_4 * fs_5_231_14 * r_8 * h2_p1;

        pc_56[k] = - e_2 * fs_85_44_66 * h6_m6 - e_3 * fs_19_286_6006 * h8_m6 + e_3 * fs_17_33_66 * r_2 * h6_m6 + e_4 * fs_1_143_6006 * r_2 * h8_m6 - e_4 * fs_1_33_66 * r_4 * h6_m6;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ab_2, pc_57 : simd::cache_line_size())
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
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];

        pc_57[k] = - e_0 * fs_195_4_14 * h2_m1 + e_1 * fs_45_28_105 * h4_m1 + e_1 * fs_390_7_14 * r_2 * h2_m1 + e_2 * fs_85_22_33 * h6_m5 + e_2 * fs_170_11_2 * h6_m1 - e_2 * fs_135_154_105 * r_2 * h4_m1 - e_2 * fs_130_7_14 * r_4 * h2_m1 - e_3 * fs_57_286_858 * h8_m5 + e_3 * fs_57_286_42 * h8_m1 - e_3 * fs_34_33_33 * r_2 * h6_m5 - e_3 * fs_136_33_2 * r_2 * h6_m1 + e_3 * fs_135_1001_105 * r_4 * h4_m1 + e_3 * fs_520_231_14 * r_6 * h2_m1 + e_4 * fs_3_143_858 * r_2 * h8_m5 - e_4 * fs_3_143_42 * r_2 * h8_m1 + e_4 * fs_2_33_33 * r_4 * h6_m5 + e_4 * fs_8_33_2 * r_4 * h6_m1 - e_4 * fs_6_1001_105 * r_6 * h4_m1 - e_4 * fs_20_231_14 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph4_p3, ph6_m4, ph6_m2, ph6_p3, ph8_m4, ph8_m2, ph8_p3, ab_2, pc_58, pc_59 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_58[k] = - e_0 * fs_195_8_35 * h2_m2 - e_1 * fs_45_1_3 * h4_m4 + e_1 * fs_405_28_21 * h4_m2 + e_1 * fs_195_7_35 * r_2 * h2_m2 + e_2 * fs_85_22_15 * h6_m4 - e_2 * fs_595_44_2 * h6_m2 + e_2 * fs_270_11_3 * r_2 * h4_m4 - e_2 * fs_1215_154_21 * r_2 * h4_m2 - e_2 * fs_65_7_35 * r_4 * h2_m2 - e_3 * fs_57_143_165 * h8_m4 - e_3 * fs_285_286_6 * h8_m2 - e_3 * fs_34_33_15 * r_2 * h6_m4 + e_3 * fs_119_33_2 * r_2 * h6_m2 - e_3 * fs_540_143_3 * r_4 * h4_m4 + e_3 * fs_1215_1001_21 * r_4 * h4_m2 + e_3 * fs_260_231_35 * r_6 * h2_m2 + e_4 * fs_6_143_165 * r_2 * h8_m4 + e_4 * fs_15_143_6 * r_2 * h8_m2 + e_4 * fs_2_33_15 * r_4 * h6_m4 - e_4 * fs_7_33_2 * r_4 * h6_m2 + e_4 * fs_24_143_3 * r_6 * h4_m4 - e_4 * fs_54_1001_21 * r_6 * h4_m2 - e_4 * fs_10_231_35 * r_8 * h2_m2;

        pc_59[k] = e_1 * f_135_2 * h4_p3 - e_2 * fs_85_22_3 * h6_p3 - e_2 * f_405_11 * r_2 * h4_p3 - e_3 * fs_95_143_66 * h8_p3 + e_3 * fs_34_33_3 * r_2 * h6_p3 + e_3 * f_810_143 * r_4 * h4_p3 + e_4 * fs_10_143_66 * r_2 * h8_p3 - e_4 * fs_2_33_3 * r_4 * h6_p3 - e_4 * f_36_143 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ab_2, pc_60 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];

        pc_60[k] = e_0 * fs_195_8_35 * h2_p2 - e_1 * fs_405_28_21 * h4_p2 - e_1 * fs_45_1_3 * h4_p4 - e_1 * fs_195_7_35 * r_2 * h2_p2 + e_2 * fs_595_44_2 * h6_p2 + e_2 * fs_85_22_15 * h6_p4 + e_2 * fs_1215_154_21 * r_2 * h4_p2 + e_2 * fs_270_11_3 * r_2 * h4_p4 + e_2 * fs_65_7_35 * r_4 * h2_p2 + e_3 * fs_285_286_6 * h8_p2 - e_3 * fs_57_143_165 * h8_p4 - e_3 * fs_119_33_2 * r_2 * h6_p2 - e_3 * fs_34_33_15 * r_2 * h6_p4 - e_3 * fs_1215_1001_21 * r_4 * h4_p2 - e_3 * fs_540_143_3 * r_4 * h4_p4 - e_3 * fs_260_231_35 * r_6 * h2_p2 - e_4 * fs_15_143_6 * r_2 * h8_p2 + e_4 * fs_6_143_165 * r_2 * h8_p4 + e_4 * fs_7_33_2 * r_4 * h6_p2 + e_4 * fs_2_33_15 * r_4 * h6_p4 + e_4 * fs_54_1001_21 * r_6 * h4_p2 + e_4 * fs_24_143_3 * r_6 * h4_p4 + e_4 * fs_10_231_35 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ph8_0, ph8_p1, ph8_p5, ph8_p6, ab_2, pc_61, pc_62 : simd::cache_line_size())
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
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_61[k] = e_0 * fs_195_4_14 * h2_p1 - e_1 * fs_45_28_105 * h4_p1 - e_1 * fs_390_7_14 * r_2 * h2_p1 - e_2 * fs_170_11_2 * h6_p1 + e_2 * fs_85_22_33 * h6_p5 + e_2 * fs_135_154_105 * r_2 * h4_p1 + e_2 * fs_130_7_14 * r_4 * h2_p1 - e_3 * fs_57_286_42 * h8_p1 - e_3 * fs_57_286_858 * h8_p5 + e_3 * fs_136_33_2 * r_2 * h6_p1 - e_3 * fs_34_33_33 * r_2 * h6_p5 - e_3 * fs_135_1001_105 * r_4 * h4_p1 - e_3 * fs_520_231_14 * r_6 * h2_p1 + e_4 * fs_3_143_42 * r_2 * h8_p1 + e_4 * fs_3_143_858 * r_2 * h8_p5 - e_4 * fs_8_33_2 * r_4 * h6_p1 + e_4 * fs_2_33_33 * r_4 * h6_p5 + e_4 * fs_6_1001_105 * r_6 * h4_p1 + e_4 * fs_20_231_14 * r_8 * h2_p1;

        pc_62[k] = e_0 * fs_195_4_7 * h2_0 + e_1 * fs_675_14_7 * h4_0 - e_1 * fs_390_7_7 * r_2 * h2_0 + e_2 * fs_85_11_7 * h6_0 - e_2 * fs_85_44_66 * h6_p6 - e_2 * fs_2025_77_7 * r_2 * h4_0 + e_2 * fs_130_7_7 * r_4 * h2_0 + e_3 * fs_38_143_7 * h8_0 - e_3 * fs_19_286_6006 * h8_p6 - e_3 * fs_68_33_7 * r_2 * h6_0 + e_3 * fs_17_33_66 * r_2 * h6_p6 + e_3 * fs_4050_1001_7 * r_4 * h4_0 - e_3 * fs_520_231_7 * r_6 * h2_0 - e_4 * fs_4_143_7 * r_2 * h8_0 + e_4 * fs_1_143_6006 * r_2 * h8_p6 + e_4 * fs_4_33_7 * r_4 * h6_0 - e_4 * fs_1_33_66 * r_4 * h6_p6 - e_4 * fs_180_1001_7 * r_6 * h4_0 + e_4 * fs_20_231_7 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2, pc_63, pc_64 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_63[k] = - e_0 * fs_195_8_42 * h2_m1 - e_1 * fs_405_28_35 * h4_m1 + e_1 * fs_195_7_42 * r_2 * h2_m1 - e_2 * fs_85_22_6 * h6_m1 + e_2 * fs_1215_154_35 * r_2 * h4_m1 - e_2 * fs_65_7_42 * r_4 * h2_m1 - e_3 * fs_19_286_10010 * h8_m7 - e_3 * fs_19_286_14 * h8_m1 + e_3 * fs_34_33_6 * r_2 * h6_m1 - e_3 * fs_1215_1001_35 * r_4 * h4_m1 + e_3 * fs_260_231_42 * r_6 * h2_m1 + e_4 * fs_1_143_10010 * r_2 * h8_m7 + e_4 * fs_1_143_14 * r_2 * h8_m1 - e_4 * fs_2_33_6 * r_4 * h6_m1 + e_4 * fs_54_1001_35 * r_6 * h4_m1 - e_4 * fs_10_231_42 * r_8 * h2_m1;

        pc_64[k] = - e_0 * fs_585_8_7 * h2_m2 + e_1 * fs_45_7_105 * h4_m2 + e_1 * fs_585_7_7 * r_2 * h2_m2 + e_2 * fs_255_44_22 * h6_m6 + e_2 * fs_255_44_10 * h6_m2 - e_2 * fs_270_77_105 * r_2 * h4_m2 - e_2 * fs_195_7_7 * r_4 * h2_m2 - e_3 * fs_19_143_2002 * h8_m6 + e_3 * fs_19_143_30 * h8_m2 - e_3 * fs_17_11_22 * r_2 * h6_m6 - e_3 * fs_17_11_10 * r_2 * h6_m2 + e_3 * fs_540_1001_105 * r_4 * h4_m2 + e_3 * fs_260_77_7 * r_6 * h2_m2 + e_4 * fs_2_143_2002 * r_2 * h8_m6 - e_4 * fs_2_143_30 * r_2 * h8_m2 + e_4 * fs_1_11_22 * r_4 * h6_m6 + e_4 * fs_1_11_10 * r_4 * h6_m2 - e_4 * fs_24_1001_105 * r_6 * h4_m2 - e_4 * fs_10_77_7 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph4_p3, ph4_p4, ph6_m3, ph6_p3, ph6_p4, ph8_m5, ph8_m3, ph8_p3, ph8_p4, ph8_p5, ab_2, pc_65, pc_66, pc_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_65[k] = e_1 * fs_45_4_3 * h4_m3 - e_2 * f_255_11 * h6_m3 - e_2 * fs_135_22_3 * r_2 * h4_m3 - e_3 * fs_19_286_4290 * h8_m5 - e_3 * fs_95_286_22 * h8_m3 + e_3 * f_68_11 * r_2 * h6_m3 + e_3 * fs_135_143_3 * r_4 * h4_m3 + e_4 * fs_1_143_4290 * r_2 * h8_m5 + e_4 * fs_5_143_22 * r_2 * h8_m3 - e_4 * f_4_11 * r_4 * h6_m3 - e_4 * fs_6_143_3 * r_6 * h4_m3;

        pc_66[k] = e_1 * f_135_1 * h4_p4 - e_2 * fs_255_22_5 * h6_p4 - e_2 * f_810_11 * r_2 * h4_p4 - e_3 * fs_76_143_55 * h8_p4 + e_3 * fs_34_11_5 * r_2 * h6_p4 + e_3 * f_1620_143 * r_4 * h4_p4 + e_4 * fs_8_143_55 * r_2 * h8_p4 - e_4 * fs_2_11_5 * r_4 * h6_p4 - e_4 * f_72_143 * r_6 * h4_p4;

        pc_67[k] = - e_1 * fs_45_4_3 * h4_p3 + e_2 * f_255_11 * h6_p3 + e_2 * fs_135_22_3 * r_2 * h4_p3 + e_3 * fs_95_286_22 * h8_p3 - e_3 * fs_19_286_4290 * h8_p5 - e_3 * f_68_11 * r_2 * h6_p3 - e_3 * fs_135_143_3 * r_4 * h4_p3 - e_4 * fs_5_143_22 * r_2 * h8_p3 + e_4 * fs_1_143_4290 * r_2 * h8_p5 + e_4 * f_4_11 * r_4 * h6_p3 + e_4 * fs_6_143_3 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph6_p6, ph8_p1, ph8_p2, ph8_p6, ph8_p7, ab_2, pc_68, pc_69 : simd::cache_line_size())
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
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_68[k] = e_0 * fs_585_8_7 * h2_p2 - e_1 * fs_45_7_105 * h4_p2 - e_1 * fs_585_7_7 * r_2 * h2_p2 - e_2 * fs_255_44_10 * h6_p2 + e_2 * fs_255_44_22 * h6_p6 + e_2 * fs_270_77_105 * r_2 * h4_p2 + e_2 * fs_195_7_7 * r_4 * h2_p2 - e_3 * fs_19_143_30 * h8_p2 - e_3 * fs_19_143_2002 * h8_p6 + e_3 * fs_17_11_10 * r_2 * h6_p2 - e_3 * fs_17_11_22 * r_2 * h6_p6 - e_3 * fs_540_1001_105 * r_4 * h4_p2 - e_3 * fs_260_77_7 * r_6 * h2_p2 + e_4 * fs_2_143_30 * r_2 * h8_p2 + e_4 * fs_2_143_2002 * r_2 * h8_p6 - e_4 * fs_1_11_10 * r_4 * h6_p2 + e_4 * fs_1_11_22 * r_4 * h6_p6 + e_4 * fs_24_1001_105 * r_6 * h4_p2 + e_4 * fs_10_77_7 * r_8 * h2_p2;

        pc_69[k] = e_0 * fs_195_8_42 * h2_p1 + e_1 * fs_405_28_35 * h4_p1 - e_1 * fs_195_7_42 * r_2 * h2_p1 + e_2 * fs_85_22_6 * h6_p1 - e_2 * fs_1215_154_35 * r_2 * h4_p1 + e_2 * fs_65_7_42 * r_4 * h2_p1 + e_3 * fs_19_286_14 * h8_p1 - e_3 * fs_19_286_10010 * h8_p7 - e_3 * fs_34_33_6 * r_2 * h6_p1 + e_3 * fs_1215_1001_35 * r_4 * h4_p1 - e_3 * fs_260_231_42 * r_6 * h2_p1 - e_4 * fs_1_143_14 * r_2 * h8_p1 + e_4 * fs_1_143_10010 * r_2 * h8_p7 + e_4 * fs_2_33_6 * r_4 * h6_p1 - e_4 * fs_54_1001_35 * r_6 * h4_p1 + e_4 * fs_10_231_42 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m3, ph4_m2, ph6_m3, ph6_m2, ph8_m8, ph8_m7, ph8_m3, ph8_m2, ab_2, pc_70, pc_71 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_70[k] = - e_0 * fs_195_8_105 * h2_m2 - e_1 * fs_675_28_7 * h4_m2 + e_1 * fs_195_7_105 * r_2 * h2_m2 - e_2 * fs_85_44_6 * h6_m2 + e_2 * fs_2025_154_7 * r_2 * h4_m2 - e_2 * fs_65_7_105 * r_4 * h2_m2 - e_3 * fs_38_143_1001 * h8_m8 - e_3 * fs_19_286_2 * h8_m2 + e_3 * fs_17_33_6 * r_2 * h6_m2 - e_3 * fs_2025_1001_7 * r_4 * h4_m2 + e_3 * fs_260_231_105 * r_6 * h2_m2 + e_4 * fs_4_143_1001 * r_2 * h8_m8 + e_4 * fs_1_143_2 * r_2 * h8_m2 - e_4 * fs_1_33_6 * r_4 * h6_m2 + e_4 * fs_90_1001_7 * r_6 * h4_m2 - e_4 * fs_10_231_105 * r_8 * h2_m2;

        pc_71[k] = e_1 * fs_225_4_3 * h4_m3 + e_2 * f_255_22 * h6_m3 - e_2 * fs_675_22_3 * r_2 * h4_m3 - e_3 * fs_19_286_6006 * h8_m7 + e_3 * fs_19_286_22 * h8_m3 - e_3 * f_34_11 * r_2 * h6_m3 + e_3 * fs_675_143_3 * r_4 * h4_m3 + e_4 * fs_1_143_6006 * r_2 * h8_m7 - e_4 * fs_1_143_22 * r_2 * h8_m3 + e_4 * f_2_11 * r_4 * h6_m3 - e_4 * fs_30_143_3 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_p4, ph6_m6, ph6_m4, ph6_p4, ph6_p5, ph6_p6, ph8_m6, ph8_m4, ph8_p4, ph8_p5, ph8_p6, ab_2, pc_72, pc_73, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_72[k] = - e_1 * fs_45_2_15 * h4_m4 - e_2 * fs_255_44_22 * h6_m6 - e_2 * fs_255_22_3 * h6_m4 + e_2 * fs_135_11_15 * r_2 * h4_m4 - e_3 * fs_19_286_2002 * h8_m6 - e_3 * fs_19_143_33 * h8_m4 + e_3 * fs_17_11_22 * r_2 * h6_m6 + e_3 * fs_34_11_3 * r_2 * h6_m4 - e_3 * fs_270_143_15 * r_4 * h4_m4 + e_4 * fs_1_143_2002 * r_2 * h8_m6 + e_4 * fs_2_143_33 * r_2 * h8_m4 - e_4 * fs_1_11_22 * r_4 * h6_m6 - e_4 * fs_2_11_3 * r_4 * h6_m4 + e_4 * fs_12_143_15 * r_6 * h4_m4;

        pc_73[k] = - e_2 * fs_255_22_11 * h6_p5 - e_3 * fs_19_143_286 * h8_p5 + e_3 * fs_34_11_11 * r_2 * h6_p5 + e_4 * fs_2_143_286 * r_2 * h8_p5 - e_4 * fs_2_11_11 * r_4 * h6_p5;

        pc_74[k] = e_1 * fs_45_2_15 * h4_p4 + e_2 * fs_255_22_3 * h6_p4 - e_2 * fs_255_44_22 * h6_p6 - e_2 * fs_135_11_15 * r_2 * h4_p4 + e_3 * fs_19_143_33 * h8_p4 - e_3 * fs_19_286_2002 * h8_p6 - e_3 * fs_34_11_3 * r_2 * h6_p4 + e_3 * fs_17_11_22 * r_2 * h6_p6 + e_3 * fs_270_143_15 * r_4 * h4_p4 - e_4 * fs_2_143_33 * r_2 * h8_p4 + e_4 * fs_1_143_2002 * r_2 * h8_p6 + e_4 * fs_2_11_3 * r_4 * h6_p4 - e_4 * fs_1_11_22 * r_4 * h6_p6 - e_4 * fs_12_143_15 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 40 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph8_p7, ph8_p8, ab_2, pc_75, pc_76 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_75[k] = - e_1 * fs_225_4_3 * h4_p3 - e_2 * f_255_22 * h6_p3 + e_2 * fs_675_22_3 * r_2 * h4_p3 - e_3 * fs_19_286_22 * h8_p3 - e_3 * fs_19_286_6006 * h8_p7 + e_3 * f_34_11 * r_2 * h6_p3 - e_3 * fs_675_143_3 * r_4 * h4_p3 + e_4 * fs_1_143_22 * r_2 * h8_p3 + e_4 * fs_1_143_6006 * r_2 * h8_p7 - e_4 * f_2_11 * r_4 * h6_p3 + e_4 * fs_30_143_3 * r_6 * h4_p3;

        pc_76[k] = e_0 * fs_195_8_105 * h2_p2 + e_1 * fs_675_28_7 * h4_p2 - e_1 * fs_195_7_105 * r_2 * h2_p2 + e_2 * fs_85_44_6 * h6_p2 - e_2 * fs_2025_154_7 * r_2 * h4_p2 + e_2 * fs_65_7_105 * r_4 * h2_p2 + e_3 * fs_19_286_2 * h8_p2 - e_3 * fs_38_143_1001 * h8_p8 - e_3 * fs_17_33_6 * r_2 * h6_p2 + e_3 * fs_2025_1001_7 * r_4 * h4_p2 - e_3 * fs_260_231_105 * r_6 * h2_p2 - e_4 * fs_1_143_2 * r_2 * h8_p2 + e_4 * fs_4_143_1001 * r_2 * h8_p8 + e_4 * fs_1_33_6 * r_4 * h6_p2 - e_4 * fs_90_1001_7 * r_6 * h4_p2 + e_4 * fs_10_231_105 * r_8 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[77] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76};

    for (size_t n = 0; n < 77; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
