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



#include "SimdKineticEnergyRecID.hpp"

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
compute_id_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 2))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_id_kinetic_energy: Basis functions must be of angular momenta six and two"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_id_kinetic_energy: Harmonics must reach angular momentum 8"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_id_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 65 * nvalues, 0.0);

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

            const auto ff_0 = fbase * bexp * bexp * bexp * bexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * bexp * bexp * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * bexp * bexp * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * bexp * bexp * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_2025_143 = 2025.0 / 143.0;
    const auto f_2025_22 = 2025.0 / 22.0;
    const auto f_24_143 = 24.0 / 143.0;
    const auto f_255_11 = 255.0 / 11.0;
    const auto f_255_22 = 255.0 / 22.0;
    const auto f_26_55 = 26.0 / 55.0;
    const auto f_270_11 = 270.0 / 11.0;
    const auto f_28_55 = 28.0 / 55.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_2_5 = 2.0 / 5.0;
    const auto f_34_11 = 34.0 / 11.0;
    const auto f_34_5 = 34.0 / 5.0;
    const auto f_357_11 = 357.0 / 11.0;
    const auto f_38_13 = 38.0 / 13.0;
    const auto f_442_55 = 442.0 / 55.0;
    const auto f_45_1 = 45.0;
    const auto f_476_55 = 476.0 / 55.0;
    const auto f_4_11 = 4.0 / 11.0;
    const auto f_4_13 = 4.0 / 13.0;
    const auto f_4_5 = 4.0 / 5.0;
    const auto f_4_55 = 4.0 / 55.0;
    const auto f_51_1 = 51.0;
    const auto f_51_11 = 51.0 / 11.0;
    const auto f_51_2 = 51.0 / 2.0;
    const auto f_532_65 = 532.0 / 65.0;
    const auto f_540_143 = 540.0 / 143.0;
    const auto f_56_65 = 56.0 / 65.0;
    const auto f_663_22 = 663.0 / 22.0;
    const auto f_675_4 = 675.0 / 4.0;
    const auto f_68_11 = 68.0 / 11.0;
    const auto f_68_5 = 68.0 / 5.0;
    const auto f_68_55 = 68.0 / 55.0;
    const auto f_90_143 = 90.0 / 143.0;
    const auto fs_102_55_10 = std::sqrt(20808.0 / 605.0);
    const auto fs_114_65_10 = std::sqrt(25992.0 / 845.0);
    const auto fs_114_65_11 = std::sqrt(142956.0 / 4225.0);
    const auto fs_114_65_21 = std::sqrt(272916.0 / 4225.0);
    const auto fs_119_55_10 = std::sqrt(28322.0 / 605.0);
    const auto fs_12_143_30 = std::sqrt(4320.0 / 20449.0);
    const auto fs_12_143_42 = std::sqrt(6048.0 / 20449.0);
    const auto fs_12_65_10 = std::sqrt(288.0 / 845.0);
    const auto fs_12_65_11 = std::sqrt(1584.0 / 4225.0);
    const auto fs_12_65_21 = std::sqrt(3024.0 / 4225.0);
    const auto fs_135_11_30 = std::sqrt(546750.0 / 121.0);
    const auto fs_135_11_42 = std::sqrt(765450.0 / 121.0);
    const auto fs_135_143_105 = std::sqrt(1913625.0 / 20449.0);
    const auto fs_135_143_15 = std::sqrt(273375.0 / 20449.0);
    const auto fs_135_143_165 = std::sqrt(273375.0 / 1859.0);
    const auto fs_135_143_210 = std::sqrt(3827250.0 / 20449.0);
    const auto fs_135_143_35 = std::sqrt(637875.0 / 20449.0);
    const auto fs_135_143_70 = std::sqrt(1275750.0 / 20449.0);
    const auto fs_135_22_105 = std::sqrt(1913625.0 / 484.0);
    const auto fs_135_22_15 = std::sqrt(273375.0 / 484.0);
    const auto fs_135_22_165 = std::sqrt(273375.0 / 44.0);
    const auto fs_135_22_210 = std::sqrt(1913625.0 / 242.0);
    const auto fs_135_22_35 = std::sqrt(637875.0 / 484.0);
    const auto fs_135_22_70 = std::sqrt(637875.0 / 242.0);
    const auto fs_135_286_10 = std::sqrt(91125.0 / 40898.0);
    const auto fs_135_286_2 = std::sqrt(18225.0 / 40898.0);
    const auto fs_135_286_330 = std::sqrt(273375.0 / 3718.0);
    const auto fs_135_286_70 = std::sqrt(637875.0 / 40898.0);
    const auto fs_135_2_3 = std::sqrt(54675.0 / 4.0);
    const auto fs_135_44_10 = std::sqrt(91125.0 / 968.0);
    const auto fs_135_44_2 = std::sqrt(18225.0 / 968.0);
    const auto fs_135_44_330 = std::sqrt(273375.0 / 88.0);
    const auto fs_135_44_70 = std::sqrt(637875.0 / 968.0);
    const auto fs_135_4_14 = std::sqrt(127575.0 / 8.0);
    const auto fs_135_4_5 = std::sqrt(91125.0 / 16.0);
    const auto fs_135_4_7 = std::sqrt(127575.0 / 16.0);
    const auto fs_135_8_110 = std::sqrt(1002375.0 / 32.0);
    const auto fs_135_8_2 = std::sqrt(18225.0 / 32.0);
    const auto fs_14_55_3 = std::sqrt(588.0 / 3025.0);
    const auto fs_152_65_7 = std::sqrt(161728.0 / 4225.0);
    const auto fs_153_22_10 = std::sqrt(117045.0 / 242.0);
    const auto fs_153_44_66 = std::sqrt(70227.0 / 88.0);
    const auto fs_16_65_7 = std::sqrt(1792.0 / 4225.0);
    const auto fs_18_143_14 = std::sqrt(4536.0 / 20449.0);
    const auto fs_18_143_5 = std::sqrt(1620.0 / 20449.0);
    const auto fs_18_143_7 = std::sqrt(2268.0 / 20449.0);
    const auto fs_19_130_10 = std::sqrt(361.0 / 1690.0);
    const auto fs_19_130_1430 = std::sqrt(3971.0 / 130.0);
    const auto fs_19_130_2 = std::sqrt(361.0 / 8450.0);
    const auto fs_19_130_2002 = std::sqrt(27797.0 / 650.0);
    const auto fs_19_130_26 = std::sqrt(361.0 / 650.0);
    const auto fs_19_130_2730 = std::sqrt(7581.0 / 130.0);
    const auto fs_19_130_30 = std::sqrt(1083.0 / 1690.0);
    const auto fs_19_130_70 = std::sqrt(2527.0 / 1690.0);
    const auto fs_19_130_910 = std::sqrt(2527.0 / 130.0);
    const auto fs_19_13_22 = std::sqrt(7942.0 / 169.0);
    const auto fs_19_26_66 = std::sqrt(11913.0 / 338.0);
    const auto fs_19_65_165 = std::sqrt(11913.0 / 845.0);
    const auto fs_19_65_210 = std::sqrt(15162.0 / 845.0);
    const auto fs_19_65_429 = std::sqrt(11913.0 / 325.0);
    const auto fs_19_65_55 = std::sqrt(3971.0 / 845.0);
    const auto fs_19_65_70 = std::sqrt(5054.0 / 845.0);
    const auto fs_19_65_91 = std::sqrt(2527.0 / 325.0);
    const auto fs_19_65_910 = std::sqrt(5054.0 / 65.0);
    const auto fs_1_13_66 = std::sqrt(66.0 / 169.0);
    const auto fs_1_65_10 = std::sqrt(2.0 / 845.0);
    const auto fs_1_65_1430 = std::sqrt(22.0 / 65.0);
    const auto fs_1_65_2 = std::sqrt(2.0 / 4225.0);
    const auto fs_1_65_2002 = std::sqrt(154.0 / 325.0);
    const auto fs_1_65_26 = std::sqrt(2.0 / 325.0);
    const auto fs_1_65_2730 = std::sqrt(42.0 / 65.0);
    const auto fs_1_65_30 = std::sqrt(6.0 / 845.0);
    const auto fs_1_65_70 = std::sqrt(14.0 / 845.0);
    const auto fs_1_65_910 = std::sqrt(14.0 / 65.0);
    const auto fs_225_4_7 = std::sqrt(354375.0 / 16.0);
    const auto fs_238_55_3 = std::sqrt(169932.0 / 3025.0);
    const auto fs_24_143_7 = std::sqrt(4032.0 / 20449.0);
    const auto fs_255_22_3 = std::sqrt(195075.0 / 484.0);
    const auto fs_266_65_3 = std::sqrt(212268.0 / 4225.0);
    const auto fs_270_11_7 = std::sqrt(510300.0 / 121.0);
    const auto fs_270_143_30 = std::sqrt(2187000.0 / 20449.0);
    const auto fs_270_143_42 = std::sqrt(3061800.0 / 20449.0);
    const auto fs_28_65_3 = std::sqrt(2352.0 / 4225.0);
    const auto fs_2_11_3 = std::sqrt(12.0 / 121.0);
    const auto fs_2_13_22 = std::sqrt(88.0 / 169.0);
    const auto fs_2_55_22 = std::sqrt(8.0 / 275.0);
    const auto fs_2_55_30 = std::sqrt(24.0 / 605.0);
    const auto fs_2_55_55 = std::sqrt(4.0 / 55.0);
    const auto fs_2_55_7 = std::sqrt(28.0 / 3025.0);
    const auto fs_2_65_165 = std::sqrt(132.0 / 845.0);
    const auto fs_2_65_210 = std::sqrt(168.0 / 845.0);
    const auto fs_2_65_429 = std::sqrt(132.0 / 325.0);
    const auto fs_2_65_55 = std::sqrt(44.0 / 845.0);
    const auto fs_2_65_70 = std::sqrt(56.0 / 845.0);
    const auto fs_2_65_91 = std::sqrt(28.0 / 325.0);
    const auto fs_2_65_910 = std::sqrt(56.0 / 65.0);
    const auto fs_30_143_7 = std::sqrt(6300.0 / 20449.0);
    const auto fs_34_11_3 = std::sqrt(3468.0 / 121.0);
    const auto fs_34_55_22 = std::sqrt(2312.0 / 275.0);
    const auto fs_34_55_30 = std::sqrt(6936.0 / 605.0);
    const auto fs_34_55_55 = std::sqrt(1156.0 / 55.0);
    const auto fs_34_55_7 = std::sqrt(8092.0 / 3025.0);
    const auto fs_357_22_3 = std::sqrt(382347.0 / 484.0);
    const auto fs_357_44_10 = std::sqrt(637245.0 / 968.0);
    const auto fs_36_143_3 = std::sqrt(3888.0 / 20449.0);
    const auto fs_38_65_110 = std::sqrt(31768.0 / 845.0);
    const auto fs_38_65_6 = std::sqrt(8664.0 / 4225.0);
    const auto fs_38_65_91 = std::sqrt(10108.0 / 325.0);
    const auto fs_3_143_10 = std::sqrt(90.0 / 20449.0);
    const auto fs_3_143_2 = std::sqrt(18.0 / 20449.0);
    const auto fs_3_143_330 = std::sqrt(270.0 / 1859.0);
    const auto fs_3_143_70 = std::sqrt(630.0 / 20449.0);
    const auto fs_3_55_66 = std::sqrt(54.0 / 275.0);
    const auto fs_3_65_110 = std::sqrt(198.0 / 845.0);
    const auto fs_3_65_70 = std::sqrt(126.0 / 845.0);
    const auto fs_405_11_3 = std::sqrt(492075.0 / 121.0);
    const auto fs_405_143_14 = std::sqrt(2296350.0 / 20449.0);
    const auto fs_405_143_5 = std::sqrt(820125.0 / 20449.0);
    const auto fs_405_143_7 = std::sqrt(1148175.0 / 20449.0);
    const auto fs_405_22_14 = std::sqrt(1148175.0 / 242.0);
    const auto fs_405_22_5 = std::sqrt(820125.0 / 484.0);
    const auto fs_405_22_7 = std::sqrt(1148175.0 / 484.0);
    const auto fs_405_286_110 = std::sqrt(820125.0 / 3718.0);
    const auto fs_405_286_2 = std::sqrt(164025.0 / 40898.0);
    const auto fs_405_44_110 = std::sqrt(820125.0 / 88.0);
    const auto fs_405_44_2 = std::sqrt(164025.0 / 968.0);
    const auto fs_45_1_7 = std::sqrt(14175.0);
    const auto fs_45_2_30 = std::sqrt(30375.0 / 2.0);
    const auto fs_45_2_42 = std::sqrt(42525.0 / 2.0);
    const auto fs_45_4_105 = std::sqrt(212625.0 / 16.0);
    const auto fs_45_4_15 = std::sqrt(30375.0 / 16.0);
    const auto fs_45_4_165 = std::sqrt(334125.0 / 16.0);
    const auto fs_45_4_210 = std::sqrt(212625.0 / 8.0);
    const auto fs_45_4_35 = std::sqrt(70875.0 / 16.0);
    const auto fs_45_4_70 = std::sqrt(70875.0 / 8.0);
    const auto fs_45_8_10 = std::sqrt(10125.0 / 32.0);
    const auto fs_45_8_2 = std::sqrt(2025.0 / 32.0);
    const auto fs_45_8_330 = std::sqrt(334125.0 / 32.0);
    const auto fs_45_8_70 = std::sqrt(70875.0 / 32.0);
    const auto fs_4_55_30 = std::sqrt(96.0 / 605.0);
    const auto fs_4_55_70 = std::sqrt(224.0 / 605.0);
    const auto fs_4_65_110 = std::sqrt(352.0 / 845.0);
    const auto fs_4_65_6 = std::sqrt(96.0 / 4225.0);
    const auto fs_4_65_91 = std::sqrt(112.0 / 325.0);
    const auto fs_51_11_30 = std::sqrt(78030.0 / 121.0);
    const auto fs_51_11_70 = std::sqrt(182070.0 / 121.0);
    const auto fs_51_22_22 = std::sqrt(2601.0 / 22.0);
    const auto fs_51_22_30 = std::sqrt(39015.0 / 242.0);
    const auto fs_51_22_55 = std::sqrt(13005.0 / 44.0);
    const auto fs_51_22_7 = std::sqrt(18207.0 / 484.0);
    const auto fs_51_55_66 = std::sqrt(15606.0 / 275.0);
    const auto fs_540_143_7 = std::sqrt(2041200.0 / 20449.0);
    const auto fs_57_130_110 = std::sqrt(35739.0 / 1690.0);
    const auto fs_57_130_70 = std::sqrt(22743.0 / 1690.0);
    const auto fs_57_13_3 = std::sqrt(9747.0 / 169.0);
    const auto fs_57_65_26 = std::sqrt(6498.0 / 325.0);
    const auto fs_57_65_7 = std::sqrt(22743.0 / 4225.0);
    const auto fs_675_143_7 = std::sqrt(3189375.0 / 20449.0);
    const auto fs_675_22_7 = std::sqrt(3189375.0 / 484.0);
    const auto fs_68_55_30 = std::sqrt(27744.0 / 605.0);
    const auto fs_68_55_70 = std::sqrt(64736.0 / 605.0);
    const auto fs_6_13_3 = std::sqrt(108.0 / 169.0);
    const auto fs_6_143_105 = std::sqrt(3780.0 / 20449.0);
    const auto fs_6_143_15 = std::sqrt(540.0 / 20449.0);
    const auto fs_6_143_165 = std::sqrt(540.0 / 1859.0);
    const auto fs_6_143_210 = std::sqrt(7560.0 / 20449.0);
    const auto fs_6_143_35 = std::sqrt(1260.0 / 20449.0);
    const auto fs_6_143_70 = std::sqrt(2520.0 / 20449.0);
    const auto fs_6_55_10 = std::sqrt(72.0 / 605.0);
    const auto fs_6_65_26 = std::sqrt(72.0 / 325.0);
    const auto fs_6_65_7 = std::sqrt(252.0 / 4225.0);
    const auto fs_7_55_10 = std::sqrt(98.0 / 605.0);
    const auto fs_810_143_3 = std::sqrt(1968300.0 / 20449.0);
    const auto fs_9_143_110 = std::sqrt(810.0 / 1859.0);
    const auto fs_9_143_2 = std::sqrt(162.0 / 20449.0);

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p4, ph6_m6, ph6_m5, ph6_p4, ph6_p5, ph8_m7, ph8_m6, ph8_m5, ph8_p4, ph8_p5, ph8_p7, ph8_p8, ab_2, pc_0, pc_1, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * fs_135_8_110 * h4_p4 + e_1 * fs_51_22_22 * h6_p4 - e_1 * fs_405_44_110 * r_2 * h4_p4 + e_2 * fs_19_130_2 * h8_p4 - e_2 * fs_19_65_910 * h8_p8 - e_2 * fs_34_55_22 * r_2 * h6_p4 + e_2 * fs_405_286_110 * r_4 * h4_p4 - e_3 * fs_1_65_2 * r_2 * h8_p4 + e_3 * fs_2_65_910 * r_2 * h8_p8 + e_3 * fs_2_55_22 * r_4 * h6_p4 - e_3 * fs_9_143_110 * r_6 * h4_p4;

        pc_1[k] = - e_1 * f_51_2 * h6_p5 - e_2 * fs_19_130_26 * h8_p5 - e_2 * fs_19_130_910 * h8_p7 + e_2 * f_34_5 * r_2 * h6_p5 + e_3 * fs_1_65_26 * r_2 * h8_p5 + e_3 * fs_1_65_910 * r_2 * h8_p7 - e_3 * f_2_5 * r_4 * h6_p5;

        pc_2[k] = e_1 * f_51_1 * h6_m6 + e_2 * fs_19_65_91 * h8_m6 - e_2 * f_68_5 * r_2 * h6_m6 - e_3 * fs_2_65_91 * r_2 * h8_m6 + e_3 * f_4_5 * r_4 * h6_m6;

        pc_3[k] = - e_1 * f_51_2 * h6_m5 + e_2 * fs_19_130_910 * h8_m7 - e_2 * fs_19_130_26 * h8_m5 + e_2 * f_34_5 * r_2 * h6_m5 - e_3 * fs_1_65_910 * r_2 * h8_m7 + e_3 * fs_1_65_26 * r_2 * h8_m5 - e_3 * f_2_5 * r_4 * h6_m5;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_p3, ph4_p4, ph6_m4, ph6_p3, ph6_p4, ph6_p6, ph8_m8, ph8_m4, ph8_p3, ph8_p4, ph8_p6, ph8_p7, ab_2, pc_4, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_4[k] = e_0 * fs_135_8_110 * h4_m4 + e_1 * fs_51_22_22 * h6_m4 - e_1 * fs_405_44_110 * r_2 * h4_m4 + e_2 * fs_19_65_910 * h8_m8 + e_2 * fs_19_130_2 * h8_m4 - e_2 * fs_34_55_22 * r_2 * h6_m4 + e_2 * fs_405_286_110 * r_4 * h4_m4 - e_3 * fs_2_65_910 * r_2 * h8_m8 - e_3 * fs_1_65_2 * r_2 * h8_m4 + e_3 * fs_2_55_22 * r_4 * h6_m4 - e_3 * fs_9_143_110 * r_6 * h4_m4;

        pc_5[k] = e_0 * fs_45_4_165 * h4_p3 + e_1 * fs_51_22_55 * h6_p3 - e_1 * fs_135_22_165 * r_2 * h4_p3 + e_2 * fs_19_130_10 * h8_p3 - e_2 * fs_19_130_2730 * h8_p7 - e_2 * fs_34_55_55 * r_2 * h6_p3 + e_2 * fs_135_143_165 * r_4 * h4_p3 - e_3 * fs_1_65_10 * r_2 * h8_p3 + e_3 * fs_1_65_2730 * r_2 * h8_p7 + e_3 * fs_2_55_55 * r_4 * h6_p3 - e_3 * fs_6_143_165 * r_6 * h4_p3;

        pc_6[k] = e_0 * fs_45_8_330 * h4_p4 - e_1 * fs_153_44_66 * h6_p4 + e_1 * f_51_2 * h6_p6 - e_1 * fs_135_44_330 * r_2 * h4_p4 - e_2 * fs_38_65_6 * h8_p4 - e_2 * fs_38_65_91 * h8_p6 + e_2 * fs_51_55_66 * r_2 * h6_p4 - e_2 * f_34_5 * r_2 * h6_p6 + e_2 * fs_135_286_330 * r_4 * h4_p4 + e_3 * fs_4_65_6 * r_2 * h8_p4 + e_3 * fs_4_65_91 * r_2 * h8_p6 - e_3 * fs_3_55_66 * r_4 * h6_p4 + e_3 * f_2_5 * r_4 * h6_p6 - e_3 * fs_3_143_330 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_m3, ph6_m6, ph6_m5, ph6_m4, ph6_m3, ph8_m7, ph8_m6, ph8_m5, ph8_m4, ph8_m3, ab_2, pc_7, pc_8, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];

        pc_7[k] = e_1 * f_51_2 * h6_m5 + e_2 * fs_57_65_26 * h8_m5 - e_2 * f_34_5 * r_2 * h6_m5 - e_3 * fs_6_65_26 * r_2 * h8_m5 + e_3 * f_2_5 * r_4 * h6_m5;

        pc_8[k] = e_0 * fs_45_8_330 * h4_m4 - e_1 * f_51_2 * h6_m6 - e_1 * fs_153_44_66 * h6_m4 - e_1 * fs_135_44_330 * r_2 * h4_m4 + e_2 * fs_38_65_91 * h8_m6 - e_2 * fs_38_65_6 * h8_m4 + e_2 * f_34_5 * r_2 * h6_m6 + e_2 * fs_51_55_66 * r_2 * h6_m4 + e_2 * fs_135_286_330 * r_4 * h4_m4 - e_3 * fs_4_65_91 * r_2 * h8_m6 + e_3 * fs_4_65_6 * r_2 * h8_m4 - e_3 * f_2_5 * r_4 * h6_m6 - e_3 * fs_3_55_66 * r_4 * h6_m4 - e_3 * fs_3_143_330 * r_6 * h4_m4;

        pc_9[k] = e_0 * fs_45_4_165 * h4_m3 + e_1 * fs_51_22_55 * h6_m3 - e_1 * fs_135_22_165 * r_2 * h4_m3 + e_2 * fs_19_130_2730 * h8_m7 + e_2 * fs_19_130_10 * h8_m3 - e_2 * fs_34_55_55 * r_2 * h6_m3 + e_2 * fs_135_143_165 * r_4 * h4_m3 - e_3 * fs_1_65_2730 * r_2 * h8_m7 - e_3 * fs_1_65_10 * r_2 * h8_m3 + e_3 * fs_2_55_55 * r_4 * h6_m3 - e_3 * fs_6_143_165 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_p2, ph4_p3, ph6_m4, ph6_p2, ph6_p3, ph6_p5, ph6_p6, ph8_m4, ph8_p2, ph8_p3, ph8_p5, ph8_p6, ab_2, pc_10, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_10[k] = e_0 * fs_45_4_105 * h4_p2 + e_1 * fs_153_22_10 * h6_p2 - e_1 * fs_51_22_22 * h6_p6 - e_1 * fs_135_22_105 * r_2 * h4_p2 + e_2 * fs_19_130_30 * h8_p2 - e_2 * fs_19_130_2002 * h8_p6 - e_2 * fs_102_55_10 * r_2 * h6_p2 + e_2 * fs_34_55_22 * r_2 * h6_p6 + e_2 * fs_135_143_105 * r_4 * h4_p2 - e_3 * fs_1_65_30 * r_2 * h8_p2 + e_3 * fs_1_65_2002 * r_2 * h8_p6 + e_3 * fs_6_55_10 * r_4 * h6_p2 - e_3 * fs_2_55_22 * r_4 * h6_p6 - e_3 * fs_6_143_105 * r_6 * h4_p2;

        pc_11[k] = e_0 * fs_45_2_30 * h4_p3 - e_1 * fs_357_44_10 * h6_p3 + e_1 * fs_153_44_66 * h6_p5 - e_1 * fs_135_11_30 * r_2 * h4_p3 - e_2 * fs_19_65_55 * h8_p3 - e_2 * fs_19_65_429 * h8_p5 + e_2 * fs_119_55_10 * r_2 * h6_p3 - e_2 * fs_51_55_66 * r_2 * h6_p5 + e_2 * fs_270_143_30 * r_4 * h4_p3 + e_3 * fs_2_65_55 * r_2 * h8_p3 + e_3 * fs_2_65_429 * r_2 * h8_p5 - e_3 * fs_7_55_10 * r_4 * h6_p3 + e_3 * fs_3_55_66 * r_4 * h6_p5 - e_3 * fs_12_143_30 * r_6 * h4_p3;

        pc_12[k] = e_0 * fs_135_4_5 * h4_m4 + e_1 * f_51_11 * h6_m4 - e_1 * fs_405_22_5 * r_2 * h4_m4 + e_2 * fs_114_65_11 * h8_m4 - e_2 * f_68_55 * r_2 * h6_m4 + e_2 * fs_405_143_5 * r_4 * h4_m4 - e_3 * fs_12_65_11 * r_2 * h8_m4 + e_3 * f_4_55 * r_4 * h6_m4 - e_3 * fs_18_143_5 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph4_m2, ph6_m6, ph6_m5, ph6_m3, ph6_m2, ph8_m6, ph8_m5, ph8_m3, ph8_m2, ab_2, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_13[k] = e_0 * fs_45_2_30 * h4_m3 - e_1 * fs_153_44_66 * h6_m5 - e_1 * fs_357_44_10 * h6_m3 - e_1 * fs_135_11_30 * r_2 * h4_m3 + e_2 * fs_19_65_429 * h8_m5 - e_2 * fs_19_65_55 * h8_m3 + e_2 * fs_51_55_66 * r_2 * h6_m5 + e_2 * fs_119_55_10 * r_2 * h6_m3 + e_2 * fs_270_143_30 * r_4 * h4_m3 - e_3 * fs_2_65_429 * r_2 * h8_m5 + e_3 * fs_2_65_55 * r_2 * h8_m3 - e_3 * fs_3_55_66 * r_4 * h6_m5 - e_3 * fs_7_55_10 * r_4 * h6_m3 - e_3 * fs_12_143_30 * r_6 * h4_m3;

        pc_14[k] = e_0 * fs_45_4_105 * h4_m2 + e_1 * fs_51_22_22 * h6_m6 + e_1 * fs_153_22_10 * h6_m2 - e_1 * fs_135_22_105 * r_2 * h4_m2 + e_2 * fs_19_130_2002 * h8_m6 + e_2 * fs_19_130_30 * h8_m2 - e_2 * fs_34_55_22 * r_2 * h6_m6 - e_2 * fs_102_55_10 * r_2 * h6_m2 + e_2 * fs_135_143_105 * r_4 * h4_m2 - e_3 * fs_1_65_2002 * r_2 * h8_m6 - e_3 * fs_1_65_30 * r_2 * h8_m2 + e_3 * fs_2_55_22 * r_4 * h6_m6 + e_3 * fs_6_55_10 * r_4 * h6_m2 - e_3 * fs_6_143_105 * r_6 * h4_m2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p1, ph4_p2, ph4_p4, ph6_p1, ph6_p2, ph6_p4, ph6_p5, ph8_p1, ph8_p2, ph8_p4, ph8_p5, ab_2, pc_15, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_15[k] = e_0 * fs_135_4_7 * h4_p1 + e_1 * fs_51_11_30 * h6_p1 - e_1 * fs_51_22_55 * h6_p5 - e_1 * fs_405_22_7 * r_2 * h4_p1 + e_2 * fs_19_130_70 * h8_p1 - e_2 * fs_19_130_1430 * h8_p5 - e_2 * fs_68_55_30 * r_2 * h6_p1 + e_2 * fs_34_55_55 * r_2 * h6_p5 + e_2 * fs_405_143_7 * r_4 * h4_p1 - e_3 * fs_1_65_70 * r_2 * h8_p1 + e_3 * fs_1_65_1430 * r_2 * h8_p5 + e_3 * fs_4_55_30 * r_4 * h6_p1 - e_3 * fs_2_55_55 * r_4 * h6_p5 - e_3 * fs_18_143_7 * r_6 * h4_p1;

        pc_16[k] = e_0 * fs_135_4_14 * h4_p2 + e_0 * fs_135_8_2 * h4_p4 - e_1 * fs_255_22_3 * h6_p2 + e_1 * fs_357_44_10 * h6_p4 - e_1 * fs_405_22_14 * r_2 * h4_p2 - e_1 * fs_405_44_2 * r_2 * h4_p4 - e_2 * f_38_13 * h8_p2 - e_2 * fs_38_65_110 * h8_p4 + e_2 * fs_34_11_3 * r_2 * h6_p2 - e_2 * fs_119_55_10 * r_2 * h6_p4 + e_2 * fs_405_143_14 * r_4 * h4_p2 + e_2 * fs_405_286_2 * r_4 * h4_p4 + e_3 * f_4_13 * r_2 * h8_p2 + e_3 * fs_4_65_110 * r_2 * h8_p4 - e_3 * fs_2_11_3 * r_4 * h6_p2 + e_3 * fs_7_55_10 * r_4 * h6_p4 - e_3 * fs_18_143_14 * r_6 * h4_p2 - e_3 * fs_9_143_2 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_m3, ph4_m2, ph6_m4, ph6_m3, ph6_m2, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_17, pc_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_17[k] = e_0 * fs_135_2_3 * h4_m3 - e_1 * f_255_22 * h6_m3 - e_1 * fs_405_11_3 * r_2 * h4_m3 + e_2 * fs_19_13_22 * h8_m3 + e_2 * f_34_11 * r_2 * h6_m3 + e_2 * fs_810_143_3 * r_4 * h4_m3 - e_3 * fs_2_13_22 * r_2 * h8_m3 - e_3 * f_2_11 * r_4 * h6_m3 - e_3 * fs_36_143_3 * r_6 * h4_m3;

        pc_18[k] = - e_0 * fs_135_8_2 * h4_m4 + e_0 * fs_135_4_14 * h4_m2 - e_1 * fs_357_44_10 * h6_m4 - e_1 * fs_255_22_3 * h6_m2 + e_1 * fs_405_44_2 * r_2 * h4_m4 - e_1 * fs_405_22_14 * r_2 * h4_m2 + e_2 * fs_38_65_110 * h8_m4 - e_2 * f_38_13 * h8_m2 + e_2 * fs_119_55_10 * r_2 * h6_m4 + e_2 * fs_34_11_3 * r_2 * h6_m2 - e_2 * fs_405_286_2 * r_4 * h4_m4 + e_2 * fs_405_143_14 * r_4 * h4_m2 - e_3 * fs_4_65_110 * r_2 * h8_m4 + e_3 * f_4_13 * r_2 * h8_m2 - e_3 * fs_7_55_10 * r_4 * h6_m4 - e_3 * fs_2_11_3 * r_4 * h6_m2 + e_3 * fs_9_143_2 * r_6 * h4_m4 - e_3 * fs_18_143_14 * r_6 * h4_m2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m1, ph4_0, ph4_p4, ph6_m5, ph6_m1, ph6_0, ph6_p4, ph8_m5, ph8_m1, ph8_0, ph8_p4, ab_2, pc_19, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_19[k] = e_0 * fs_135_4_7 * h4_m1 + e_1 * fs_51_22_55 * h6_m5 + e_1 * fs_51_11_30 * h6_m1 - e_1 * fs_405_22_7 * r_2 * h4_m1 + e_2 * fs_19_130_1430 * h8_m5 + e_2 * fs_19_130_70 * h8_m1 - e_2 * fs_34_55_55 * r_2 * h6_m5 - e_2 * fs_68_55_30 * r_2 * h6_m1 + e_2 * fs_405_143_7 * r_4 * h4_m1 - e_3 * fs_1_65_1430 * r_2 * h8_m5 - e_3 * fs_1_65_70 * r_2 * h8_m1 + e_3 * fs_2_55_55 * r_4 * h6_m5 + e_3 * fs_4_55_30 * r_4 * h6_m1 - e_3 * fs_18_143_7 * r_6 * h4_m1;

        pc_20[k] = e_0 * fs_45_4_70 * h4_0 - e_0 * fs_45_8_2 * h4_p4 + e_1 * fs_51_11_70 * h6_0 - e_1 * fs_153_22_10 * h6_p4 - e_1 * fs_135_22_70 * r_2 * h4_0 + e_1 * fs_135_44_2 * r_2 * h4_p4 + e_2 * fs_19_65_70 * h8_0 - e_2 * fs_57_130_110 * h8_p4 - e_2 * fs_68_55_70 * r_2 * h6_0 + e_2 * fs_102_55_10 * r_2 * h6_p4 + e_2 * fs_135_143_70 * r_4 * h4_0 - e_2 * fs_135_286_2 * r_4 * h4_p4 - e_3 * fs_2_65_70 * r_2 * h8_0 + e_3 * fs_3_65_110 * r_2 * h8_p4 + e_3 * fs_4_55_70 * r_4 * h6_0 - e_3 * fs_6_55_10 * r_4 * h6_p4 - e_3 * fs_6_143_70 * r_6 * h4_0 + e_3 * fs_3_143_2 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m2, ph4_p1, ph4_p3, ph6_m2, ph6_p1, ph6_p3, ph8_m2, ph8_p1, ph8_p3, ab_2, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_21[k] = e_0 * fs_45_1_7 * h4_p1 + e_0 * f_45_1 * h4_p3 - e_1 * fs_51_22_30 * h6_p1 + e_1 * fs_255_22_3 * h6_p3 - e_1 * fs_270_11_7 * r_2 * h4_p1 - e_1 * f_270_11 * r_2 * h4_p3 - e_2 * fs_57_130_70 * h8_p1 - e_2 * fs_19_26_66 * h8_p3 + e_2 * fs_34_55_30 * r_2 * h6_p1 - e_2 * fs_34_11_3 * r_2 * h6_p3 + e_2 * fs_540_143_7 * r_4 * h4_p1 + e_2 * f_540_143 * r_4 * h4_p3 + e_3 * fs_3_65_70 * r_2 * h8_p1 + e_3 * fs_1_13_66 * r_2 * h8_p3 - e_3 * fs_2_55_30 * r_4 * h6_p1 + e_3 * fs_2_11_3 * r_4 * h6_p3 - e_3 * fs_24_143_7 * r_6 * h4_p1 - e_3 * f_24_143 * r_6 * h4_p3;

        pc_22[k] = e_0 * fs_45_2_42 * h4_m2 - e_1 * f_255_11 * h6_m2 - e_1 * fs_135_11_42 * r_2 * h4_m2 + e_2 * fs_57_13_3 * h8_m2 + e_2 * f_68_11 * r_2 * h6_m2 + e_2 * fs_270_143_42 * r_4 * h4_m2 - e_3 * fs_6_13_3 * r_2 * h8_m2 - e_3 * f_4_11 * r_4 * h6_m2 - e_3 * fs_12_143_42 * r_6 * h4_m2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_m3, ph4_m1, ph6_m4, ph6_m3, ph6_m1, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_23, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];

        pc_23[k] = - e_0 * f_45_1 * h4_m3 + e_0 * fs_45_1_7 * h4_m1 - e_1 * fs_255_22_3 * h6_m3 - e_1 * fs_51_22_30 * h6_m1 + e_1 * f_270_11 * r_2 * h4_m3 - e_1 * fs_270_11_7 * r_2 * h4_m1 + e_2 * fs_19_26_66 * h8_m3 - e_2 * fs_57_130_70 * h8_m1 + e_2 * fs_34_11_3 * r_2 * h6_m3 + e_2 * fs_34_55_30 * r_2 * h6_m1 - e_2 * f_540_143 * r_4 * h4_m3 + e_2 * fs_540_143_7 * r_4 * h4_m1 - e_3 * fs_1_13_66 * r_2 * h8_m3 + e_3 * fs_3_65_70 * r_2 * h8_m1 - e_3 * fs_2_11_3 * r_4 * h6_m3 - e_3 * fs_2_55_30 * r_4 * h6_m1 + e_3 * f_24_143 * r_6 * h4_m3 - e_3 * fs_24_143_7 * r_6 * h4_m1;

        pc_24[k] = e_0 * fs_45_8_2 * h4_m4 + e_1 * fs_153_22_10 * h6_m4 - e_1 * fs_135_44_2 * r_2 * h4_m4 + e_2 * fs_57_130_110 * h8_m4 - e_2 * fs_102_55_10 * r_2 * h6_m4 + e_2 * fs_135_286_2 * r_4 * h4_m4 - e_3 * fs_3_65_110 * r_2 * h8_m4 + e_3 * fs_6_55_10 * r_4 * h6_m4 - e_3 * fs_3_143_2 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph8_0, ph8_p1, ph8_p2, ph8_p3, ab_2, pc_25, pc_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_25[k] = - e_0 * fs_45_8_70 * h4_p1 - e_0 * fs_45_8_10 * h4_p3 - e_1 * fs_357_22_3 * h6_p1 - e_1 * fs_51_11_30 * h6_p3 + e_1 * fs_135_44_70 * r_2 * h4_p1 + e_1 * fs_135_44_10 * r_2 * h4_p3 - e_2 * fs_57_65_7 * h8_p1 - e_2 * fs_19_65_165 * h8_p3 + e_2 * fs_238_55_3 * r_2 * h6_p1 + e_2 * fs_68_55_30 * r_2 * h6_p3 - e_2 * fs_135_286_70 * r_4 * h4_p1 - e_2 * fs_135_286_10 * r_4 * h4_p3 + e_3 * fs_6_65_7 * r_2 * h8_p1 + e_3 * fs_2_65_165 * r_2 * h8_p3 - e_3 * fs_14_55_3 * r_4 * h6_p1 - e_3 * fs_4_55_30 * r_4 * h6_p3 + e_3 * fs_3_143_70 * r_6 * h4_p1 + e_3 * fs_3_143_10 * r_6 * h4_p3;

        pc_26[k] = e_0 * fs_225_4_7 * h4_0 + e_0 * fs_45_4_35 * h4_p2 - e_1 * fs_51_22_7 * h6_0 + e_1 * fs_51_22_30 * h6_p2 - e_1 * fs_675_22_7 * r_2 * h4_0 - e_1 * fs_135_22_35 * r_2 * h4_p2 - e_2 * fs_152_65_7 * h8_0 - e_2 * fs_114_65_10 * h8_p2 + e_2 * fs_34_55_7 * r_2 * h6_0 - e_2 * fs_34_55_30 * r_2 * h6_p2 + e_2 * fs_675_143_7 * r_4 * h4_0 + e_2 * fs_135_143_35 * r_4 * h4_p2 + e_3 * fs_16_65_7 * r_2 * h8_0 + e_3 * fs_12_65_10 * r_2 * h8_p2 - e_3 * fs_2_55_7 * r_4 * h6_0 + e_3 * fs_2_55_30 * r_4 * h6_p2 - e_3 * fs_30_143_7 * r_6 * h4_0 - e_3 * fs_6_143_35 * r_6 * h4_p2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph4_m2, ph4_m1, ph4_0, ph6_m3, ph6_m2, ph6_m1, ph6_0, ph8_m3, ph8_m2, ph8_m1, ph8_0, ab_2, pc_27, pc_28, pc_29, pc_30, pc_31, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];

        pc_27[k] = e_0 * fs_45_4_210 * h4_m1 - e_1 * f_663_22 * h6_m1 - e_1 * fs_135_22_210 * r_2 * h4_m1 + e_2 * fs_114_65_21 * h8_m1 + e_2 * f_442_55 * r_2 * h6_m1 + e_2 * fs_135_143_210 * r_4 * h4_m1 - e_3 * fs_12_65_21 * r_2 * h8_m1 - e_3 * f_26_55 * r_4 * h6_m1 - e_3 * fs_6_143_210 * r_6 * h4_m1;

        pc_28[k] = - e_0 * fs_45_4_35 * h4_m2 - e_1 * fs_51_22_30 * h6_m2 + e_1 * fs_135_22_35 * r_2 * h4_m2 + e_2 * fs_114_65_10 * h8_m2 + e_2 * fs_34_55_30 * r_2 * h6_m2 - e_2 * fs_135_143_35 * r_4 * h4_m2 - e_3 * fs_12_65_10 * r_2 * h8_m2 - e_3 * fs_2_55_30 * r_4 * h6_m2 + e_3 * fs_6_143_35 * r_6 * h4_m2;

        pc_29[k] = e_0 * fs_45_8_10 * h4_m3 + e_0 * fs_45_8_70 * h4_m1 + e_1 * fs_51_11_30 * h6_m3 + e_1 * fs_357_22_3 * h6_m1 - e_1 * fs_135_44_10 * r_2 * h4_m3 - e_1 * fs_135_44_70 * r_2 * h4_m1 + e_2 * fs_19_65_165 * h8_m3 + e_2 * fs_57_65_7 * h8_m1 - e_2 * fs_68_55_30 * r_2 * h6_m3 - e_2 * fs_238_55_3 * r_2 * h6_m1 + e_2 * fs_135_286_10 * r_4 * h4_m3 + e_2 * fs_135_286_70 * r_4 * h4_m1 - e_3 * fs_2_65_165 * r_2 * h8_m3 - e_3 * fs_6_65_7 * r_2 * h8_m1 + e_3 * fs_4_55_30 * r_4 * h6_m3 + e_3 * fs_14_55_3 * r_4 * h6_m1 - e_3 * fs_3_143_10 * r_6 * h4_m3 - e_3 * fs_3_143_70 * r_6 * h4_m1;

        pc_30[k] = e_0 * fs_45_4_15 * h4_m2 + e_1 * fs_51_11_70 * h6_m2 - e_1 * fs_135_22_15 * r_2 * h4_m2 + e_2 * fs_19_65_210 * h8_m2 - e_2 * fs_68_55_70 * r_2 * h6_m2 + e_2 * fs_135_143_15 * r_4 * h4_m2 - e_3 * fs_2_65_210 * r_2 * h8_m2 + e_3 * fs_4_55_70 * r_4 * h6_m2 - e_3 * fs_6_143_15 * r_6 * h4_m2;

        pc_31[k] = - e_0 * fs_45_2_30 * h4_m1 - e_1 * fs_51_22_7 * h6_m1 + e_1 * fs_135_11_30 * r_2 * h4_m1 + e_2 * fs_266_65_3 * h8_m1 + e_2 * fs_34_55_7 * r_2 * h6_m1 - e_2 * fs_270_143_30 * r_4 * h4_m1 - e_3 * fs_28_65_3 * r_2 * h8_m1 - e_3 * fs_2_55_7 * r_4 * h6_m1 + e_3 * fs_12_143_30 * r_6 * h4_m1;

        pc_32[k] = e_0 * f_675_4 * h4_0 - e_1 * f_357_11 * h6_0 - e_1 * f_2025_22 * r_2 * h4_0 + e_2 * f_532_65 * h8_0 + e_2 * f_476_55 * r_2 * h6_0 + e_2 * f_2025_143 * r_4 * h4_0 - e_3 * f_56_65 * r_2 * h8_0 - e_3 * f_28_55 * r_4 * h6_0 - e_3 * f_90_143 * r_6 * h4_0;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph4_m1, ph4_p1, ph4_p2, ph6_m3, ph6_m1, ph6_p1, ph6_p2, ph8_m3, ph8_m1, ph8_p1, ph8_p2, ab_2, pc_33, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_33[k] = - e_0 * fs_45_2_30 * h4_p1 - e_1 * fs_51_22_7 * h6_p1 + e_1 * fs_135_11_30 * r_2 * h4_p1 + e_2 * fs_266_65_3 * h8_p1 + e_2 * fs_34_55_7 * r_2 * h6_p1 - e_2 * fs_270_143_30 * r_4 * h4_p1 - e_3 * fs_28_65_3 * r_2 * h8_p1 - e_3 * fs_2_55_7 * r_4 * h6_p1 + e_3 * fs_12_143_30 * r_6 * h4_p1;

        pc_34[k] = e_0 * fs_45_4_15 * h4_p2 + e_1 * fs_51_11_70 * h6_p2 - e_1 * fs_135_22_15 * r_2 * h4_p2 + e_2 * fs_19_65_210 * h8_p2 - e_2 * fs_68_55_70 * r_2 * h6_p2 + e_2 * fs_135_143_15 * r_4 * h4_p2 - e_3 * fs_2_65_210 * r_2 * h8_p2 + e_3 * fs_4_55_70 * r_4 * h6_p2 - e_3 * fs_6_143_15 * r_6 * h4_p2;

        pc_35[k] = e_0 * fs_45_8_10 * h4_m3 - e_0 * fs_45_8_70 * h4_m1 + e_1 * fs_51_11_30 * h6_m3 - e_1 * fs_357_22_3 * h6_m1 - e_1 * fs_135_44_10 * r_2 * h4_m3 + e_1 * fs_135_44_70 * r_2 * h4_m1 + e_2 * fs_19_65_165 * h8_m3 - e_2 * fs_57_65_7 * h8_m1 - e_2 * fs_68_55_30 * r_2 * h6_m3 + e_2 * fs_238_55_3 * r_2 * h6_m1 + e_2 * fs_135_286_10 * r_4 * h4_m3 - e_2 * fs_135_286_70 * r_4 * h4_m1 - e_3 * fs_2_65_165 * r_2 * h8_m3 + e_3 * fs_6_65_7 * r_2 * h8_m1 + e_3 * fs_4_55_30 * r_4 * h6_m3 - e_3 * fs_14_55_3 * r_4 * h6_m1 - e_3 * fs_3_143_10 * r_6 * h4_m3 + e_3 * fs_3_143_70 * r_6 * h4_m1;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m2, ph4_0, ph4_p1, ph4_p2, ph6_m2, ph6_0, ph6_p1, ph6_p2, ph8_m2, ph8_0, ph8_p1, ph8_p2, ab_2, pc_36, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m2 = ph4_m2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_36[k] = - e_0 * fs_45_4_35 * h4_m2 - e_1 * fs_51_22_30 * h6_m2 + e_1 * fs_135_22_35 * r_2 * h4_m2 + e_2 * fs_114_65_10 * h8_m2 + e_2 * fs_34_55_30 * r_2 * h6_m2 - e_2 * fs_135_143_35 * r_4 * h4_m2 - e_3 * fs_12_65_10 * r_2 * h8_m2 - e_3 * fs_2_55_30 * r_4 * h6_m2 + e_3 * fs_6_143_35 * r_6 * h4_m2;

        pc_37[k] = e_0 * fs_45_4_210 * h4_p1 - e_1 * f_663_22 * h6_p1 - e_1 * fs_135_22_210 * r_2 * h4_p1 + e_2 * fs_114_65_21 * h8_p1 + e_2 * f_442_55 * r_2 * h6_p1 + e_2 * fs_135_143_210 * r_4 * h4_p1 - e_3 * fs_12_65_21 * r_2 * h8_p1 - e_3 * f_26_55 * r_4 * h6_p1 - e_3 * fs_6_143_210 * r_6 * h4_p1;

        pc_38[k] = e_0 * fs_225_4_7 * h4_0 - e_0 * fs_45_4_35 * h4_p2 - e_1 * fs_51_22_7 * h6_0 - e_1 * fs_51_22_30 * h6_p2 - e_1 * fs_675_22_7 * r_2 * h4_0 + e_1 * fs_135_22_35 * r_2 * h4_p2 - e_2 * fs_152_65_7 * h8_0 + e_2 * fs_114_65_10 * h8_p2 + e_2 * fs_34_55_7 * r_2 * h6_0 + e_2 * fs_34_55_30 * r_2 * h6_p2 + e_2 * fs_675_143_7 * r_4 * h4_0 - e_2 * fs_135_143_35 * r_4 * h4_p2 + e_3 * fs_16_65_7 * r_2 * h8_0 - e_3 * fs_12_65_10 * r_2 * h8_p2 - e_3 * fs_2_55_7 * r_4 * h6_0 - e_3 * fs_2_55_30 * r_4 * h6_p2 - e_3 * fs_30_143_7 * r_6 * h4_0 + e_3 * fs_6_143_35 * r_6 * h4_p2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_p1, ph4_p3, ph6_m4, ph6_p1, ph6_p3, ph8_m4, ph8_p1, ph8_p3, ab_2, pc_39, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_39[k] = - e_0 * fs_45_8_70 * h4_p1 + e_0 * fs_45_8_10 * h4_p3 - e_1 * fs_357_22_3 * h6_p1 + e_1 * fs_51_11_30 * h6_p3 + e_1 * fs_135_44_70 * r_2 * h4_p1 - e_1 * fs_135_44_10 * r_2 * h4_p3 - e_2 * fs_57_65_7 * h8_p1 + e_2 * fs_19_65_165 * h8_p3 + e_2 * fs_238_55_3 * r_2 * h6_p1 - e_2 * fs_68_55_30 * r_2 * h6_p3 - e_2 * fs_135_286_70 * r_4 * h4_p1 + e_2 * fs_135_286_10 * r_4 * h4_p3 + e_3 * fs_6_65_7 * r_2 * h8_p1 - e_3 * fs_2_65_165 * r_2 * h8_p3 - e_3 * fs_14_55_3 * r_4 * h6_p1 + e_3 * fs_4_55_30 * r_4 * h6_p3 + e_3 * fs_3_143_70 * r_6 * h4_p1 - e_3 * fs_3_143_10 * r_6 * h4_p3;

        pc_40[k] = e_0 * fs_45_8_2 * h4_m4 + e_1 * fs_153_22_10 * h6_m4 - e_1 * fs_135_44_2 * r_2 * h4_m4 + e_2 * fs_57_130_110 * h8_m4 - e_2 * fs_102_55_10 * r_2 * h6_m4 + e_2 * fs_135_286_2 * r_4 * h4_m4 - e_3 * fs_3_65_110 * r_2 * h8_m4 + e_3 * fs_6_55_10 * r_4 * h6_m4 - e_3 * fs_3_143_2 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph4_m1, ph4_p2, ph6_m3, ph6_m1, ph6_p2, ph8_m3, ph8_m1, ph8_p2, ab_2, pc_41, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_41[k] = - e_0 * f_45_1 * h4_m3 - e_0 * fs_45_1_7 * h4_m1 - e_1 * fs_255_22_3 * h6_m3 + e_1 * fs_51_22_30 * h6_m1 + e_1 * f_270_11 * r_2 * h4_m3 + e_1 * fs_270_11_7 * r_2 * h4_m1 + e_2 * fs_19_26_66 * h8_m3 + e_2 * fs_57_130_70 * h8_m1 + e_2 * fs_34_11_3 * r_2 * h6_m3 - e_2 * fs_34_55_30 * r_2 * h6_m1 - e_2 * f_540_143 * r_4 * h4_m3 - e_2 * fs_540_143_7 * r_4 * h4_m1 - e_3 * fs_1_13_66 * r_2 * h8_m3 - e_3 * fs_3_65_70 * r_2 * h8_m1 - e_3 * fs_2_11_3 * r_4 * h6_m3 + e_3 * fs_2_55_30 * r_4 * h6_m1 + e_3 * f_24_143 * r_6 * h4_m3 + e_3 * fs_24_143_7 * r_6 * h4_m1;

        pc_42[k] = e_0 * fs_45_2_42 * h4_p2 - e_1 * f_255_11 * h6_p2 - e_1 * fs_135_11_42 * r_2 * h4_p2 + e_2 * fs_57_13_3 * h8_p2 + e_2 * f_68_11 * r_2 * h6_p2 + e_2 * fs_270_143_42 * r_4 * h4_p2 - e_3 * fs_6_13_3 * r_2 * h8_p2 - e_3 * f_4_11 * r_4 * h6_p2 - e_3 * fs_12_143_42 * r_6 * h4_p2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_0, ph4_p1, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p3, ph6_p4, ph8_0, ph8_p1, ph8_p3, ph8_p4, ab_2, pc_43, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_43[k] = e_0 * fs_45_1_7 * h4_p1 - e_0 * f_45_1 * h4_p3 - e_1 * fs_51_22_30 * h6_p1 - e_1 * fs_255_22_3 * h6_p3 - e_1 * fs_270_11_7 * r_2 * h4_p1 + e_1 * f_270_11 * r_2 * h4_p3 - e_2 * fs_57_130_70 * h8_p1 + e_2 * fs_19_26_66 * h8_p3 + e_2 * fs_34_55_30 * r_2 * h6_p1 + e_2 * fs_34_11_3 * r_2 * h6_p3 + e_2 * fs_540_143_7 * r_4 * h4_p1 - e_2 * f_540_143 * r_4 * h4_p3 + e_3 * fs_3_65_70 * r_2 * h8_p1 - e_3 * fs_1_13_66 * r_2 * h8_p3 - e_3 * fs_2_55_30 * r_4 * h6_p1 - e_3 * fs_2_11_3 * r_4 * h6_p3 - e_3 * fs_24_143_7 * r_6 * h4_p1 + e_3 * f_24_143 * r_6 * h4_p3;

        pc_44[k] = e_0 * fs_45_4_70 * h4_0 + e_0 * fs_45_8_2 * h4_p4 + e_1 * fs_51_11_70 * h6_0 + e_1 * fs_153_22_10 * h6_p4 - e_1 * fs_135_22_70 * r_2 * h4_0 - e_1 * fs_135_44_2 * r_2 * h4_p4 + e_2 * fs_19_65_70 * h8_0 + e_2 * fs_57_130_110 * h8_p4 - e_2 * fs_68_55_70 * r_2 * h6_0 - e_2 * fs_102_55_10 * r_2 * h6_p4 + e_2 * fs_135_143_70 * r_4 * h4_0 + e_2 * fs_135_286_2 * r_4 * h4_p4 - e_3 * fs_2_65_70 * r_2 * h8_0 - e_3 * fs_3_65_110 * r_2 * h8_p4 + e_3 * fs_4_55_70 * r_4 * h6_0 + e_3 * fs_6_55_10 * r_4 * h6_p4 - e_3 * fs_6_143_70 * r_6 * h4_0 - e_3 * fs_3_143_2 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m2, ph8_m1, ab_2, pc_45, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_45[k] = - e_0 * fs_135_4_7 * h4_m1 + e_1 * fs_51_22_55 * h6_m5 - e_1 * fs_51_11_30 * h6_m1 + e_1 * fs_405_22_7 * r_2 * h4_m1 + e_2 * fs_19_130_1430 * h8_m5 - e_2 * fs_19_130_70 * h8_m1 - e_2 * fs_34_55_55 * r_2 * h6_m5 + e_2 * fs_68_55_30 * r_2 * h6_m1 - e_2 * fs_405_143_7 * r_4 * h4_m1 - e_3 * fs_1_65_1430 * r_2 * h8_m5 + e_3 * fs_1_65_70 * r_2 * h8_m1 + e_3 * fs_2_55_55 * r_4 * h6_m5 - e_3 * fs_4_55_30 * r_4 * h6_m1 + e_3 * fs_18_143_7 * r_6 * h4_m1;

        pc_46[k] = - e_0 * fs_135_8_2 * h4_m4 - e_0 * fs_135_4_14 * h4_m2 - e_1 * fs_357_44_10 * h6_m4 + e_1 * fs_255_22_3 * h6_m2 + e_1 * fs_405_44_2 * r_2 * h4_m4 + e_1 * fs_405_22_14 * r_2 * h4_m2 + e_2 * fs_38_65_110 * h8_m4 + e_2 * f_38_13 * h8_m2 + e_2 * fs_119_55_10 * r_2 * h6_m4 - e_2 * fs_34_11_3 * r_2 * h6_m2 - e_2 * fs_405_286_2 * r_4 * h4_m4 - e_2 * fs_405_143_14 * r_4 * h4_m2 - e_3 * fs_4_65_110 * r_2 * h8_m4 - e_3 * f_4_13 * r_2 * h8_m2 - e_3 * fs_7_55_10 * r_4 * h6_m4 + e_3 * fs_2_11_3 * r_4 * h6_m2 + e_3 * fs_9_143_2 * r_6 * h4_m4 + e_3 * fs_18_143_14 * r_6 * h4_m2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph8_p2, ph8_p3, ph8_p4, ab_2, pc_47, pc_48 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_47[k] = e_0 * fs_135_2_3 * h4_p3 - e_1 * f_255_22 * h6_p3 - e_1 * fs_405_11_3 * r_2 * h4_p3 + e_2 * fs_19_13_22 * h8_p3 + e_2 * f_34_11 * r_2 * h6_p3 + e_2 * fs_810_143_3 * r_4 * h4_p3 - e_3 * fs_2_13_22 * r_2 * h8_p3 - e_3 * f_2_11 * r_4 * h6_p3 - e_3 * fs_36_143_3 * r_6 * h4_p3;

        pc_48[k] = e_0 * fs_135_4_14 * h4_p2 - e_0 * fs_135_8_2 * h4_p4 - e_1 * fs_255_22_3 * h6_p2 - e_1 * fs_357_44_10 * h6_p4 - e_1 * fs_405_22_14 * r_2 * h4_p2 + e_1 * fs_405_44_2 * r_2 * h4_p4 - e_2 * f_38_13 * h8_p2 + e_2 * fs_38_65_110 * h8_p4 + e_2 * fs_34_11_3 * r_2 * h6_p2 + e_2 * fs_119_55_10 * r_2 * h6_p4 + e_2 * fs_405_143_14 * r_4 * h4_p2 - e_2 * fs_405_286_2 * r_4 * h4_p4 + e_3 * f_4_13 * r_2 * h8_p2 - e_3 * fs_4_65_110 * r_2 * h8_p4 - e_3 * fs_2_11_3 * r_4 * h6_p2 - e_3 * fs_7_55_10 * r_4 * h6_p4 - e_3 * fs_18_143_14 * r_6 * h4_p2 + e_3 * fs_9_143_2 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m2, ph4_p1, ph6_m6, ph6_m2, ph6_p1, ph6_p5, ph8_m6, ph8_m2, ph8_p1, ph8_p5, ab_2, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];

        pc_49[k] = e_0 * fs_135_4_7 * h4_p1 + e_1 * fs_51_11_30 * h6_p1 + e_1 * fs_51_22_55 * h6_p5 - e_1 * fs_405_22_7 * r_2 * h4_p1 + e_2 * fs_19_130_70 * h8_p1 + e_2 * fs_19_130_1430 * h8_p5 - e_2 * fs_68_55_30 * r_2 * h6_p1 - e_2 * fs_34_55_55 * r_2 * h6_p5 + e_2 * fs_405_143_7 * r_4 * h4_p1 - e_3 * fs_1_65_70 * r_2 * h8_p1 - e_3 * fs_1_65_1430 * r_2 * h8_p5 + e_3 * fs_4_55_30 * r_4 * h6_p1 + e_3 * fs_2_55_55 * r_4 * h6_p5 - e_3 * fs_18_143_7 * r_6 * h4_p1;

        pc_50[k] = - e_0 * fs_45_4_105 * h4_m2 + e_1 * fs_51_22_22 * h6_m6 - e_1 * fs_153_22_10 * h6_m2 + e_1 * fs_135_22_105 * r_2 * h4_m2 + e_2 * fs_19_130_2002 * h8_m6 - e_2 * fs_19_130_30 * h8_m2 - e_2 * fs_34_55_22 * r_2 * h6_m6 + e_2 * fs_102_55_10 * r_2 * h6_m2 - e_2 * fs_135_143_105 * r_4 * h4_m2 - e_3 * fs_1_65_2002 * r_2 * h8_m6 + e_3 * fs_1_65_30 * r_2 * h8_m2 + e_3 * fs_2_55_22 * r_4 * h6_m6 - e_3 * fs_6_55_10 * r_4 * h6_m2 + e_3 * fs_6_143_105 * r_6 * h4_m2;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph4_p3, ph4_p4, ph6_m5, ph6_m3, ph6_p3, ph6_p4, ph6_p5, ph8_m5, ph8_m3, ph8_p3, ph8_p4, ph8_p5, ab_2, pc_51, pc_52, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_51[k] = - e_0 * fs_45_2_30 * h4_m3 - e_1 * fs_153_44_66 * h6_m5 + e_1 * fs_357_44_10 * h6_m3 + e_1 * fs_135_11_30 * r_2 * h4_m3 + e_2 * fs_19_65_429 * h8_m5 + e_2 * fs_19_65_55 * h8_m3 + e_2 * fs_51_55_66 * r_2 * h6_m5 - e_2 * fs_119_55_10 * r_2 * h6_m3 - e_2 * fs_270_143_30 * r_4 * h4_m3 - e_3 * fs_2_65_429 * r_2 * h8_m5 - e_3 * fs_2_65_55 * r_2 * h8_m3 - e_3 * fs_3_55_66 * r_4 * h6_m5 + e_3 * fs_7_55_10 * r_4 * h6_m3 + e_3 * fs_12_143_30 * r_6 * h4_m3;

        pc_52[k] = e_0 * fs_135_4_5 * h4_p4 + e_1 * f_51_11 * h6_p4 - e_1 * fs_405_22_5 * r_2 * h4_p4 + e_2 * fs_114_65_11 * h8_p4 - e_2 * f_68_55 * r_2 * h6_p4 + e_2 * fs_405_143_5 * r_4 * h4_p4 - e_3 * fs_12_65_11 * r_2 * h8_p4 + e_3 * f_4_55 * r_4 * h6_p4 - e_3 * fs_18_143_5 * r_6 * h4_p4;

        pc_53[k] = e_0 * fs_45_2_30 * h4_p3 - e_1 * fs_357_44_10 * h6_p3 - e_1 * fs_153_44_66 * h6_p5 - e_1 * fs_135_11_30 * r_2 * h4_p3 - e_2 * fs_19_65_55 * h8_p3 + e_2 * fs_19_65_429 * h8_p5 + e_2 * fs_119_55_10 * r_2 * h6_p3 + e_2 * fs_51_55_66 * r_2 * h6_p5 + e_2 * fs_270_143_30 * r_4 * h4_p3 + e_3 * fs_2_65_55 * r_2 * h8_p3 - e_3 * fs_2_65_429 * r_2 * h8_p5 - e_3 * fs_7_55_10 * r_4 * h6_p3 - e_3 * fs_3_55_66 * r_4 * h6_p5 - e_3 * fs_12_143_30 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m3, ph4_p2, ph6_m3, ph6_p2, ph6_p6, ph8_m7, ph8_m3, ph8_p2, ph8_p6, ab_2, pc_54, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];

        pc_54[k] = e_0 * fs_45_4_105 * h4_p2 + e_1 * fs_153_22_10 * h6_p2 + e_1 * fs_51_22_22 * h6_p6 - e_1 * fs_135_22_105 * r_2 * h4_p2 + e_2 * fs_19_130_30 * h8_p2 + e_2 * fs_19_130_2002 * h8_p6 - e_2 * fs_102_55_10 * r_2 * h6_p2 - e_2 * fs_34_55_22 * r_2 * h6_p6 + e_2 * fs_135_143_105 * r_4 * h4_p2 - e_3 * fs_1_65_30 * r_2 * h8_p2 - e_3 * fs_1_65_2002 * r_2 * h8_p6 + e_3 * fs_6_55_10 * r_4 * h6_p2 + e_3 * fs_2_55_22 * r_4 * h6_p6 - e_3 * fs_6_143_105 * r_6 * h4_p2;

        pc_55[k] = - e_0 * fs_45_4_165 * h4_m3 - e_1 * fs_51_22_55 * h6_m3 + e_1 * fs_135_22_165 * r_2 * h4_m3 + e_2 * fs_19_130_2730 * h8_m7 - e_2 * fs_19_130_10 * h8_m3 + e_2 * fs_34_55_55 * r_2 * h6_m3 - e_2 * fs_135_143_165 * r_4 * h4_m3 - e_3 * fs_1_65_2730 * r_2 * h8_m7 + e_3 * fs_1_65_10 * r_2 * h8_m3 - e_3 * fs_2_55_55 * r_4 * h6_m3 + e_3 * fs_6_143_165 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_p4, ph6_m6, ph6_m4, ph6_p4, ph6_p5, ph6_p6, ph8_m6, ph8_m4, ph8_p4, ph8_p5, ph8_p6, ab_2, pc_56, pc_57, pc_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

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

        pc_56[k] = - e_0 * fs_45_8_330 * h4_m4 - e_1 * f_51_2 * h6_m6 + e_1 * fs_153_44_66 * h6_m4 + e_1 * fs_135_44_330 * r_2 * h4_m4 + e_2 * fs_38_65_91 * h8_m6 + e_2 * fs_38_65_6 * h8_m4 + e_2 * f_34_5 * r_2 * h6_m6 - e_2 * fs_51_55_66 * r_2 * h6_m4 - e_2 * fs_135_286_330 * r_4 * h4_m4 - e_3 * fs_4_65_91 * r_2 * h8_m6 - e_3 * fs_4_65_6 * r_2 * h8_m4 - e_3 * f_2_5 * r_4 * h6_m6 + e_3 * fs_3_55_66 * r_4 * h6_m4 + e_3 * fs_3_143_330 * r_6 * h4_m4;

        pc_57[k] = e_1 * f_51_2 * h6_p5 + e_2 * fs_57_65_26 * h8_p5 - e_2 * f_34_5 * r_2 * h6_p5 - e_3 * fs_6_65_26 * r_2 * h8_p5 + e_3 * f_2_5 * r_4 * h6_p5;

        pc_58[k] = e_0 * fs_45_8_330 * h4_p4 - e_1 * fs_153_44_66 * h6_p4 - e_1 * f_51_2 * h6_p6 - e_1 * fs_135_44_330 * r_2 * h4_p4 - e_2 * fs_38_65_6 * h8_p4 + e_2 * fs_38_65_91 * h8_p6 + e_2 * fs_51_55_66 * r_2 * h6_p4 + e_2 * f_34_5 * r_2 * h6_p6 + e_2 * fs_135_286_330 * r_4 * h4_p4 + e_3 * fs_4_65_6 * r_2 * h8_p4 - e_3 * fs_4_65_91 * r_2 * h8_p6 - e_3 * fs_3_55_66 * r_4 * h6_p4 - e_3 * f_2_5 * r_4 * h6_p6 - e_3 * fs_3_143_330 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_m4, ph4_p3, ph6_m5, ph6_m4, ph6_p3, ph6_p6, ph8_m8, ph8_m7, ph8_m5, ph8_m4, ph8_p3, ph8_p6, ph8_p7, ab_2, pc_59, pc_60, pc_61, pc_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_59[k] = e_0 * fs_45_4_165 * h4_p3 + e_1 * fs_51_22_55 * h6_p3 - e_1 * fs_135_22_165 * r_2 * h4_p3 + e_2 * fs_19_130_10 * h8_p3 + e_2 * fs_19_130_2730 * h8_p7 - e_2 * fs_34_55_55 * r_2 * h6_p3 + e_2 * fs_135_143_165 * r_4 * h4_p3 - e_3 * fs_1_65_10 * r_2 * h8_p3 - e_3 * fs_1_65_2730 * r_2 * h8_p7 + e_3 * fs_2_55_55 * r_4 * h6_p3 - e_3 * fs_6_143_165 * r_6 * h4_p3;

        pc_60[k] = - e_0 * fs_135_8_110 * h4_m4 - e_1 * fs_51_22_22 * h6_m4 + e_1 * fs_405_44_110 * r_2 * h4_m4 + e_2 * fs_19_65_910 * h8_m8 - e_2 * fs_19_130_2 * h8_m4 + e_2 * fs_34_55_22 * r_2 * h6_m4 - e_2 * fs_405_286_110 * r_4 * h4_m4 - e_3 * fs_2_65_910 * r_2 * h8_m8 + e_3 * fs_1_65_2 * r_2 * h8_m4 - e_3 * fs_2_55_22 * r_4 * h6_m4 + e_3 * fs_9_143_110 * r_6 * h4_m4;

        pc_61[k] = e_1 * f_51_2 * h6_m5 + e_2 * fs_19_130_910 * h8_m7 + e_2 * fs_19_130_26 * h8_m5 - e_2 * f_34_5 * r_2 * h6_m5 - e_3 * fs_1_65_910 * r_2 * h8_m7 - e_3 * fs_1_65_26 * r_2 * h8_m5 + e_3 * f_2_5 * r_4 * h6_m5;

        pc_62[k] = e_1 * f_51_1 * h6_p6 + e_2 * fs_19_65_91 * h8_p6 - e_2 * f_68_5 * r_2 * h6_p6 - e_3 * fs_2_65_91 * r_2 * h8_p6 + e_3 * f_4_5 * r_4 * h6_p6;
    }

    // NOTE: the angular components are formed in 25 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph4_p4, ph6_p4, ph6_p5, ph8_p4, ph8_p5, ph8_p7, ph8_p8, ab_2, pc_63, pc_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_63[k] = - e_1 * f_51_2 * h6_p5 - e_2 * fs_19_130_26 * h8_p5 + e_2 * fs_19_130_910 * h8_p7 + e_2 * f_34_5 * r_2 * h6_p5 + e_3 * fs_1_65_26 * r_2 * h8_p5 - e_3 * fs_1_65_910 * r_2 * h8_p7 - e_3 * f_2_5 * r_4 * h6_p5;

        pc_64[k] = e_0 * fs_135_8_110 * h4_p4 + e_1 * fs_51_22_22 * h6_p4 - e_1 * fs_405_44_110 * r_2 * h4_p4 + e_2 * fs_19_130_2 * h8_p4 + e_2 * fs_19_65_910 * h8_p8 - e_2 * fs_34_55_22 * r_2 * h6_p4 + e_2 * fs_405_286_110 * r_4 * h4_p4 - e_3 * fs_1_65_2 * r_2 * h8_p4 - e_3 * fs_2_65_910 * r_2 * h8_p8 + e_3 * fs_2_55_22 * r_4 * h6_p4 - e_3 * fs_9_143_110 * r_6 * h4_p4;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[65] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64};

    for (size_t n = 0; n < 65; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
