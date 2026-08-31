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



#include "SimdKineticEnergyRecII.hpp"

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
compute_ii_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ii_kinetic_energy: Basis functions must be of angular momenta six and six"));
    }

    if (harmonics.size() < 12)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ii_kinetic_energy: Harmonics must reach angular momentum 12"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ii_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 169 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time.

    auto buffer = CSimdMatrix(8, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);
    auto *pe_6 = buffer.data(6);
    auto *pe_7 = buffer.data(7);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);
    std::fill(pe_4, pe_4 + nmax, 0.0);
    std::fill(pe_5, pe_5 + nmax, 0.0);
    std::fill(pe_6, pe_6 + nmax, 0.0);
    std::fill(pe_7, pe_7 + nmax, 0.0);

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

            const auto ff_0 = fbase * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_6 = fbase * fmu * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_7 = fbase * fmu * fmu * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ab_2 : simd::cache_line_size())
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
                pe_7[k] += ff_7 * fterm;
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
    const auto *ph12_m12 = harmonics[11].data(0);
    const auto *ph12_m11 = harmonics[11].data(1);
    const auto *ph12_m10 = harmonics[11].data(2);
    const auto *ph12_m9 = harmonics[11].data(3);
    const auto *ph12_m8 = harmonics[11].data(4);
    const auto *ph12_m7 = harmonics[11].data(5);
    const auto *ph12_m6 = harmonics[11].data(6);
    const auto *ph12_m5 = harmonics[11].data(7);
    const auto *ph12_m4 = harmonics[11].data(8);
    const auto *ph12_m3 = harmonics[11].data(9);
    const auto *ph12_m2 = harmonics[11].data(10);
    const auto *ph12_m1 = harmonics[11].data(11);
    const auto *ph12_0 = harmonics[11].data(12);
    const auto *ph12_p1 = harmonics[11].data(13);
    const auto *ph12_p2 = harmonics[11].data(14);
    const auto *ph12_p3 = harmonics[11].data(15);
    const auto *ph12_p4 = harmonics[11].data(16);
    const auto *ph12_p5 = harmonics[11].data(17);
    const auto *ph12_p6 = harmonics[11].data(18);
    const auto *ph12_p7 = harmonics[11].data(19);
    const auto *ph12_p8 = harmonics[11].data(20);
    const auto *ph12_p9 = harmonics[11].data(21);
    const auto *ph12_p10 = harmonics[11].data(22);
    const auto *ph12_p11 = harmonics[11].data(23);
    const auto *ph12_p12 = harmonics[11].data(24);

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
    auto *pc_13 = values + 14 * nvalues;
    auto *pc_14 = values + 15 * nvalues;
    auto *pc_15 = values + 16 * nvalues;
    auto *pc_16 = values + 17 * nvalues;
    auto *pc_17 = values + 18 * nvalues;
    auto *pc_18 = values + 19 * nvalues;
    auto *pc_19 = values + 20 * nvalues;
    auto *pc_20 = values + 21 * nvalues;
    auto *pc_21 = values + 22 * nvalues;
    auto *pc_22 = values + 23 * nvalues;
    auto *pc_23 = values + 24 * nvalues;
    auto *pc_24 = values + 25 * nvalues;
    auto *pc_25 = values + 28 * nvalues;
    auto *pc_26 = values + 29 * nvalues;
    auto *pc_27 = values + 30 * nvalues;
    auto *pc_28 = values + 31 * nvalues;
    auto *pc_29 = values + 32 * nvalues;
    auto *pc_30 = values + 33 * nvalues;
    auto *pc_31 = values + 34 * nvalues;
    auto *pc_32 = values + 35 * nvalues;
    auto *pc_33 = values + 36 * nvalues;
    auto *pc_34 = values + 37 * nvalues;
    auto *pc_35 = values + 38 * nvalues;
    auto *pc_36 = values + 42 * nvalues;
    auto *pc_37 = values + 43 * nvalues;
    auto *pc_38 = values + 44 * nvalues;
    auto *pc_39 = values + 45 * nvalues;
    auto *pc_40 = values + 46 * nvalues;
    auto *pc_41 = values + 47 * nvalues;
    auto *pc_42 = values + 48 * nvalues;
    auto *pc_43 = values + 49 * nvalues;
    auto *pc_44 = values + 50 * nvalues;
    auto *pc_45 = values + 51 * nvalues;
    auto *pc_46 = values + 56 * nvalues;
    auto *pc_47 = values + 57 * nvalues;
    auto *pc_48 = values + 58 * nvalues;
    auto *pc_49 = values + 59 * nvalues;
    auto *pc_50 = values + 60 * nvalues;
    auto *pc_51 = values + 61 * nvalues;
    auto *pc_52 = values + 62 * nvalues;
    auto *pc_53 = values + 63 * nvalues;
    auto *pc_54 = values + 64 * nvalues;
    auto *pc_55 = values + 70 * nvalues;
    auto *pc_56 = values + 71 * nvalues;
    auto *pc_57 = values + 72 * nvalues;
    auto *pc_58 = values + 73 * nvalues;
    auto *pc_59 = values + 74 * nvalues;
    auto *pc_60 = values + 75 * nvalues;
    auto *pc_61 = values + 76 * nvalues;
    auto *pc_62 = values + 77 * nvalues;
    auto *pc_63 = values + 84 * nvalues;
    auto *pc_64 = values + 85 * nvalues;
    auto *pc_65 = values + 86 * nvalues;
    auto *pc_66 = values + 87 * nvalues;
    auto *pc_67 = values + 88 * nvalues;
    auto *pc_68 = values + 89 * nvalues;
    auto *pc_69 = values + 90 * nvalues;
    auto *pc_70 = values + 98 * nvalues;
    auto *pc_71 = values + 99 * nvalues;
    auto *pc_72 = values + 100 * nvalues;
    auto *pc_73 = values + 101 * nvalues;
    auto *pc_74 = values + 102 * nvalues;
    auto *pc_75 = values + 103 * nvalues;
    auto *pc_76 = values + 112 * nvalues;
    auto *pc_77 = values + 113 * nvalues;
    auto *pc_78 = values + 114 * nvalues;
    auto *pc_79 = values + 115 * nvalues;
    auto *pc_80 = values + 116 * nvalues;
    auto *pc_81 = values + 126 * nvalues;
    auto *pc_82 = values + 127 * nvalues;
    auto *pc_83 = values + 128 * nvalues;
    auto *pc_84 = values + 129 * nvalues;
    auto *pc_85 = values + 140 * nvalues;
    auto *pc_86 = values + 141 * nvalues;
    auto *pc_87 = values + 142 * nvalues;
    auto *pc_88 = values + 154 * nvalues;
    auto *pc_89 = values + 155 * nvalues;
    auto *pc_90 = values + 168 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_100_2717 = 100.0 / 2717.0;
    const auto f_100_323 = 100.0 / 323.0;
    const auto f_1020_143 = 1020.0 / 143.0;
    const auto f_102_11 = 102.0 / 11.0;
    const auto f_102_13 = 102.0 / 13.0;
    const auto f_1050_1 = 1050.0;
    const auto f_105_13 = 105.0 / 13.0;
    const auto f_10692_96577 = 10692.0 / 96577.0;
    const auto f_10_143 = 10.0 / 143.0;
    const auto f_10_247 = 10.0 / 247.0;
    const auto f_112455_16 = 112455.0 / 16.0;
    const auto f_115425_143 = 115425.0 / 143.0;
    const auto f_115425_44 = 115425.0 / 44.0;
    const auto f_12075_286 = 12075.0 / 286.0;
    const auto f_12075_52 = 12075.0 / 52.0;
    const auto f_12600_187 = 12600.0 / 187.0;
    const auto f_1260_17 = 1260.0 / 17.0;
    const auto f_12750_11 = 12750.0 / 11.0;
    const auto f_1275_1 = 1275.0;
    const auto f_12825_13 = 12825.0 / 13.0;
    const auto f_12825_4 = 12825.0 / 4.0;
    const auto f_130725_4199 = 130725.0 / 4199.0;
    const auto f_136800_143 = 136800.0 / 143.0;
    const auto f_13986_96577 = 13986.0 / 96577.0;
    const auto f_1428_143 = 1428.0 / 143.0;
    const auto f_14520_96577 = 14520.0 / 96577.0;
    const auto f_15105_143 = 15105.0 / 143.0;
    const auto f_15105_2431 = 15105.0 / 2431.0;
    const auto f_15390_143 = 15390.0 / 143.0;
    const auto f_15390_2431 = 15390.0 / 2431.0;
    const auto f_154350_96577 = 154350.0 / 96577.0;
    const auto f_155925_64 = 155925.0 / 64.0;
    const auto f_15750_11 = 15750.0 / 11.0;
    const auto f_1575_2 = 1575.0 / 2.0;
    const auto f_16065_16 = 16065.0 / 16.0;
    const auto f_160_3553 = 160.0 / 3553.0;
    const auto f_165375_4199 = 165375.0 / 4199.0;
    const auto f_16575_11 = 16575.0 / 11.0;
    const auto f_16800_3553 = 16800.0 / 3553.0;
    const auto f_1680_11 = 1680.0 / 11.0;
    const auto f_1680_323 = 1680.0 / 323.0;
    const auto f_1710_13 = 1710.0 / 13.0;
    const auto f_1710_221 = 1710.0 / 221.0;
    const auto f_171465_572 = 171465.0 / 572.0;
    const auto f_1725_247 = 1725.0 / 247.0;
    const auto f_174825_8398 = 174825.0 / 8398.0;
    const auto f_176715_16 = 176715.0 / 16.0;
    const auto f_176715_32 = 176715.0 / 32.0;
    const auto f_17850_11 = 17850.0 / 11.0;
    const auto f_179550_143 = 179550.0 / 143.0;
    const auto f_17955_4 = 17955.0 / 4.0;
    const auto f_18240_143 = 18240.0 / 143.0;
    const auto f_18240_2431 = 18240.0 / 2431.0;
    const auto f_18900_7429 = 18900.0 / 7429.0;
    const auto f_18900_96577 = 18900.0 / 96577.0;
    const auto f_1890_5083 = 1890.0 / 5083.0;
    const auto f_19125_143 = 19125.0 / 143.0;
    const auto f_19125_286 = 19125.0 / 286.0;
    const auto f_19125_4 = 19125.0 / 4.0;
    const auto f_19125_8 = 19125.0 / 8.0;
    const auto f_196020_96577 = 196020.0 / 96577.0;
    const auto f_198450_4199 = 198450.0 / 4199.0;
    const auto f_204_13 = 204.0 / 13.0;
    const auto f_204_143 = 204.0 / 143.0;
    const auto f_205200_143 = 205200.0 / 143.0;
    const auto f_208845_32 = 208845.0 / 32.0;
    const auto f_20916_96577 = 20916.0 / 96577.0;
    const auto f_20_143 = 20.0 / 143.0;
    const auto f_20_209 = 20.0 / 209.0;
    const auto f_214935_572 = 214935.0 / 572.0;
    const auto f_218295_16 = 218295.0 / 16.0;
    const auto f_226575_286 = 226575.0 / 286.0;
    const auto f_226575_88 = 226575.0 / 88.0;
    const auto f_23085_8 = 23085.0 / 8.0;
    const auto f_23625_442 = 23625.0 / 442.0;
    const auto f_23940_143 = 23940.0 / 143.0;
    const auto f_23940_2431 = 23940.0 / 2431.0;
    const auto f_24150_2717 = 24150.0 / 2717.0;
    const auto f_2415_22 = 2415.0 / 22.0;
    const auto f_2415_52 = 2415.0 / 52.0;
    const auto f_24495_2717 = 24495.0 / 2717.0;
    const auto f_25200_187 = 25200.0 / 187.0;
    const auto f_253575_2717 = 253575.0 / 2717.0;
    const auto f_2550_1 = 2550.0;
    const auto f_2550_11 = 2550.0 / 11.0;
    const auto f_2565_13 = 2565.0 / 13.0;
    const auto f_2565_221 = 2565.0 / 221.0;
    const auto f_26460_96577 = 26460.0 / 96577.0;
    const auto f_26775_143 = 26775.0 / 143.0;
    const auto f_26775_4 = 26775.0 / 4.0;
    const auto f_27090_187 = 27090.0 / 187.0;
    const auto f_27360_143 = 27360.0 / 143.0;
    const auto f_27360_2431 = 27360.0 / 2431.0;
    const auto f_28215_8 = 28215.0 / 8.0;
    const auto f_285_13 = 285.0 / 13.0;
    const auto f_285_221 = 285.0 / 221.0;
    const auto f_28_143 = 28.0 / 143.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_30705_2717 = 30705.0 / 2717.0;
    const auto f_3150_11 = 3150.0 / 11.0;
    const auto f_3150_17 = 3150.0 / 17.0;
    const auto f_315_2 = 315.0 / 2.0;
    const auto f_31752_96577 = 31752.0 / 96577.0;
    const auto f_318_2431 = 318.0 / 2431.0;
    const auto f_324_2431 = 324.0 / 2431.0;
    const auto f_32670_96577 = 32670.0 / 96577.0;
    const auto f_33600_3553 = 33600.0 / 3553.0;
    const auto f_34200_11 = 34200.0 / 11.0;
    const auto f_3420_1 = 3420.0;
    const auto f_34425_4 = 34425.0 / 4.0;
    const auto f_34425_8 = 34425.0 / 8.0;
    const auto f_3450_2717 = 3450.0 / 2717.0;
    const auto f_345_143 = 345.0 / 143.0;
    const auto f_345_247 = 345.0 / 247.0;
    const auto f_349650_96577 = 349650.0 / 96577.0;
    const auto f_36120_3553 = 36120.0 / 3553.0;
    const auto f_36225_2717 = 36225.0 / 2717.0;
    const auto f_36225_494 = 36225.0 / 494.0;
    const auto f_363825_32 = 363825.0 / 32.0;
    const auto f_36_221 = 36.0 / 221.0;
    const auto f_3825_13 = 3825.0 / 13.0;
    const auto f_3825_143 = 3825.0 / 143.0;
    const auto f_3825_22 = 3825.0 / 22.0;
    const auto f_3825_26 = 3825.0 / 26.0;
    const auto f_3825_4 = 3825.0 / 4.0;
    const auto f_38475_26 = 38475.0 / 26.0;
    const auto f_38475_8 = 38475.0 / 8.0;
    const auto f_384_2431 = 384.0 / 2431.0;
    const auto f_400_3553 = 400.0 / 3553.0;
    const auto f_40_323 = 40.0 / 323.0;
    const auto f_4200_11 = 4200.0 / 11.0;
    const auto f_4200_323 = 4200.0 / 323.0;
    const auto f_42075_4 = 42075.0 / 4.0;
    const auto f_42075_8 = 42075.0 / 8.0;
    const auto f_420_1 = 420.0;
    const auto f_4275_26 = 4275.0 / 26.0;
    const auto f_4275_8 = 4275.0 / 8.0;
    const auto f_4356_96577 = 4356.0 / 96577.0;
    const auto f_441045_96577 = 441045.0 / 96577.0;
    const auto f_45315_16 = 45315.0 / 16.0;
    const auto f_45885_572 = 45885.0 / 572.0;
    const auto f_47250_5083 = 47250.0 / 5083.0;
    const auto f_4725_323 = 4725.0 / 323.0;
    const auto f_4725_4199 = 4725.0 / 4199.0;
    const auto f_48195_4 = 48195.0 / 4.0;
    const auto f_49725_8 = 49725.0 / 8.0;
    const auto f_4_13 = 4.0 / 13.0;
    const auto f_4_143 = 4.0 / 143.0;
    const auto f_5040_187 = 5040.0 / 187.0;
    const auto f_504_2431 = 504.0 / 2431.0;
    const auto f_50_247 = 50.0 / 247.0;
    const auto f_510_143 = 510.0 / 143.0;
    const auto f_51300_11 = 51300.0 / 11.0;
    const auto f_5130_1 = 5130.0;
    const auto f_514395_5434 = 514395.0 / 5434.0;
    const auto f_51975_8 = 51975.0 / 8.0;
    const auto f_52272_96577 = 52272.0 / 96577.0;
    const auto f_522900_96577 = 522900.0 / 96577.0;
    const auto f_54_221 = 54.0 / 221.0;
    const auto f_576_2431 = 576.0 / 2431.0;
    const auto f_5775_4 = 5775.0 / 4.0;
    const auto f_58806_96577 = 58806.0 / 96577.0;
    const auto f_60984_96577 = 60984.0 / 96577.0;
    const auto f_6174_96577 = 6174.0 / 96577.0;
    const auto f_6375_11 = 6375.0 / 11.0;
    const auto f_644805_5434 = 644805.0 / 5434.0;
    const auto f_661500_96577 = 661500.0 / 96577.0;
    const auto f_66_96577 = 66.0 / 96577.0;
    const auto f_6720_3553 = 6720.0 / 3553.0;
    const auto f_67725_44 = 67725.0 / 44.0;
    const auto f_6885_4 = 6885.0 / 4.0;
    const auto f_690_209 = 690.0 / 209.0;
    const auto f_6_221 = 6.0 / 221.0;
    const auto f_700_2717 = 700.0 / 2717.0;
    const auto f_705672_96577 = 705672.0 / 96577.0;
    const auto f_710_2717 = 710.0 / 2717.0;
    const auto f_7245_209 = 7245.0 / 209.0;
    const auto f_7245_286 = 7245.0 / 286.0;
    const auto f_7245_494 = 7245.0 / 494.0;
    const auto f_756_7429 = 756.0 / 7429.0;
    const auto f_756_96577 = 756.0 / 96577.0;
    const auto f_75735_4 = 75735.0 / 4.0;
    const auto f_75735_8 = 75735.0 / 8.0;
    const auto f_77175_8398 = 77175.0 / 8398.0;
    const auto f_7875_11 = 7875.0 / 11.0;
    const auto f_7875_4 = 7875.0 / 4.0;
    const auto f_792_96577 = 792.0 / 96577.0;
    const auto f_793800_96577 = 793800.0 / 96577.0;
    const auto f_800_3553 = 800.0 / 3553.0;
    const auto f_80325_16 = 80325.0 / 16.0;
    const auto f_80325_32 = 80325.0 / 32.0;
    const auto f_823284_96577 = 823284.0 / 96577.0;
    const auto f_8400_11 = 8400.0 / 11.0;
    const auto f_84525_286 = 84525.0 / 286.0;
    const auto f_84645_16 = 84645.0 / 16.0;
    const auto f_860_3553 = 860.0 / 3553.0;
    const auto f_890_2717 = 890.0 / 2717.0;
    const auto f_891_96577 = 891.0 / 96577.0;
    const auto f_89505_8 = 89505.0 / 8.0;
    const auto f_89775_22 = 89775.0 / 22.0;
    const auto f_9030_11 = 9030.0 / 11.0;
    const auto f_9405_16 = 9405.0 / 16.0;
    const auto fs_10080_3553_7 = std::sqrt(711244800.0 / 12623809.0);
    const auto fs_100_2717_21 = std::sqrt(210000.0 / 7382089.0);
    const auto fs_102_143_22 = std::sqrt(20808.0 / 1859.0);
    const auto fs_102_143_30 = std::sqrt(312120.0 / 20449.0);
    const auto fs_102_143_55 = std::sqrt(52020.0 / 1859.0);
    const auto fs_102_143_7 = std::sqrt(72828.0 / 20449.0);
    const auto fs_1035_2717_462 = std::sqrt(44991450.0 / 671099.0);
    const auto fs_1035_2717_715 = std::sqrt(5356125.0 / 51623.0);
    const auto fs_1035_5434_462 = std::sqrt(22495725.0 / 1342198.0);
    const auto fs_10692_96577_3003 = std::sqrt(26407657584.0 / 717470533.0);
    const auto fs_10692_96577_663 = std::sqrt(342956592.0 / 42204149.0);
    const auto fs_10_143_3 = std::sqrt(300.0 / 20449.0);
    const auto fs_10_143_7 = std::sqrt(700.0 / 20449.0);
    const auto fs_10_247_13 = std::sqrt(100.0 / 4693.0);
    const auto fs_10_247_21 = std::sqrt(2100.0 / 61009.0);
    const auto fs_10_247_3 = std::sqrt(300.0 / 61009.0);
    const auto fs_10_247_35 = std::sqrt(3500.0 / 61009.0);
    const auto fs_10_247_39 = std::sqrt(300.0 / 4693.0);
    const auto fs_10_247_91 = std::sqrt(700.0 / 4693.0);
    const auto fs_10_2717_1001 = std::sqrt(700.0 / 51623.0);
    const auto fs_10_2717_1155 = std::sqrt(10500.0 / 671099.0);
    const auto fs_10_2717_3003 = std::sqrt(2100.0 / 51623.0);
    const auto fs_10_2717_385 = std::sqrt(3500.0 / 671099.0);
    const auto fs_10_2717_715 = std::sqrt(500.0 / 51623.0);
    const auto fs_10_2717_858 = std::sqrt(600.0 / 51623.0);
    const auto fs_11025_4199_165 = std::sqrt(20055853125.0 / 17631601.0);
    const auto fs_11025_4199_39 = std::sqrt(364651875.0 / 1356277.0);
    const auto fs_11025_8398_429 = std::sqrt(4011170625.0 / 5425108.0);
    const auto fs_11025_8398_442 = std::sqrt(121550625.0 / 159562.0);
    const auto fs_112455_32_3 = std::sqrt(37938381075.0 / 1024.0);
    const auto fs_112455_64_10 = std::sqrt(63230635125.0 / 2048.0);
    const auto fs_1134_96577_10 = std::sqrt(12859560.0 / 9327116929.0);
    const auto fs_11475_16_66 = std::sqrt(4345295625.0 / 128.0);
    const auto fs_11475_286_10 = std::sqrt(658378125.0 / 40898.0);
    const auto fs_11475_572_66 = std::sqrt(395026875.0 / 14872.0);
    const auto fs_11475_8_10 = std::sqrt(658378125.0 / 32.0);
    const auto fs_11970_143_3 = std::sqrt(429842700.0 / 20449.0);
    const auto fs_11970_143_5 = std::sqrt(716404500.0 / 20449.0);
    const auto fs_11970_2431_3 = std::sqrt(429842700.0 / 5909761.0);
    const auto fs_11970_2431_5 = std::sqrt(716404500.0 / 5909761.0);
    const auto fs_12075_286_21 = std::sqrt(3061918125.0 / 81796.0);
    const auto fs_120_2717_22 = std::sqrt(28800.0 / 671099.0);
    const auto fs_12474_96577_2145 = std::sqrt(25674111540.0 / 717470533.0);
    const auto fs_12474_96577_442 = std::sqrt(311201352.0 / 42204149.0);
    const auto fs_12600_96577_2145 = std::sqrt(26195400000.0 / 717470533.0);
    const auto fs_12600_96577_3315 = std::sqrt(2381400000.0 / 42204149.0);
    const auto fs_1260_187_210 = std::sqrt(333396000.0 / 34969.0);
    const auto fs_1260_187_385 = std::sqrt(55566000.0 / 3179.0);
    const auto fs_1260_187_462 = std::sqrt(66679200.0 / 3179.0);
    const auto fs_1260_3553_462 = std::sqrt(66679200.0 / 1147619.0);
    const auto fs_126_2431_14 = std::sqrt(222264.0 / 5909761.0);
    const auto fs_126_2431_3 = std::sqrt(47628.0 / 5909761.0);
    const auto fs_126_5083_273 = std::sqrt(333396.0 / 1987453.0);
    const auto fs_126_5083_66 = std::sqrt(1047816.0 / 25836889.0);
    const auto fs_126_96577_165 = std::sqrt(2619540.0 / 9327116929.0);
    const auto fs_126_96577_20995 = std::sqrt(79380.0 / 2221271.0);
    const auto fs_126_96577_663 = std::sqrt(47628.0 / 42204149.0);
    const auto fs_126_96577_92378 = std::sqrt(349272.0 / 2221271.0);
    const auto fs_1275_11_22 = std::sqrt(3251250.0 / 11.0);
    const auto fs_1275_11_30 = std::sqrt(48768750.0 / 121.0);
    const auto fs_1275_11_55 = std::sqrt(8128125.0 / 11.0);
    const auto fs_1275_11_7 = std::sqrt(11379375.0 / 121.0);
    const auto fs_12825_16_30 = std::sqrt(2467209375.0 / 128.0);
    const auto fs_12825_22_2 = std::sqrt(164480625.0 / 242.0);
    const auto fs_12825_22_33 = std::sqrt(493441875.0 / 44.0);
    const auto fs_12825_286_55 = std::sqrt(822403125.0 / 7436.0);
    const auto fs_12825_286_77 = std::sqrt(1151364375.0 / 7436.0);
    const auto fs_12825_44_10 = std::sqrt(822403125.0 / 968.0);
    const auto fs_12825_52_30 = std::sqrt(2467209375.0 / 1352.0);
    const auto fs_12825_88_55 = std::sqrt(822403125.0 / 704.0);
    const auto fs_12825_88_77 = std::sqrt(1151364375.0 / 704.0);
    const auto fs_12_2431_210 = std::sqrt(30240.0 / 5909761.0);
    const auto fs_132_96577_138567 = std::sqrt(574992.0 / 2221271.0);
    const auto fs_132_96577_176358 = std::sqrt(731808.0 / 2221271.0);
    const auto fs_132_96577_30030 = std::sqrt(40249440.0 / 717470533.0);
    const auto fs_132_96577_455 = std::sqrt(609840.0 / 717470533.0);
    const auto fs_132_96577_4641 = std::sqrt(365904.0 / 42204149.0);
    const auto fs_132_96577_51051 = std::sqrt(4024944.0 / 42204149.0);
    const auto fs_13365_96577_2002 = std::sqrt(27507976650.0 / 717470533.0);
    const auto fs_135450_96577_22 = std::sqrt(403627455000.0 / 9327116929.0);
    const auto fs_135_2717_10 = std::sqrt(182250.0 / 7382089.0);
    const auto fs_1380_2717_21 = std::sqrt(39992400.0 / 7382089.0);
    const auto fs_1386_96577_1430 = std::sqrt(211309560.0 / 717470533.0);
    const auto fs_140_2717_22 = std::sqrt(39200.0 / 671099.0);
    const auto fs_14175_4199_221 = std::sqrt(200930625.0 / 79781.0);
    const auto fs_14175_8398_10 = std::sqrt(1004653125.0 / 35263202.0);
    const auto fs_1425_143_210 = std::sqrt(426431250.0 / 20449.0);
    const auto fs_1425_2431_210 = std::sqrt(426431250.0 / 5909761.0);
    const auto fs_14490_2717_21 = std::sqrt(4409162100.0 / 7382089.0);
    const auto fs_14_143_3 = std::sqrt(588.0 / 20449.0);
    const auto fs_1512_96577_22 = std::sqrt(50295168.0 / 9327116929.0);
    const auto fs_1512_96577_455 = std::sqrt(80015040.0 / 717470533.0);
    const auto fs_152145_5434_11 = std::sqrt(23148101025.0 / 2684396.0);
    const auto fs_153_143_66 = std::sqrt(140454.0 / 1859.0);
    const auto fs_15750_96577_858 = std::sqrt(16372125000.0 / 717470533.0);
    const auto fs_1575_187_70 = std::sqrt(173643750.0 / 34969.0);
    const auto fs_1575_22_210 = std::sqrt(260465625.0 / 242.0);
    const auto fs_1575_22_385 = std::sqrt(86821875.0 / 44.0);
    const auto fs_1575_22_462 = std::sqrt(52093125.0 / 22.0);
    const auto fs_1575_323_65 = std::sqrt(161240625.0 / 104329.0);
    const auto fs_1575_4199_12155 = std::sqrt(136434375.0 / 79781.0);
    const auto fs_1575_4199_12597 = std::sqrt(7441875.0 / 4199.0);
    const auto fs_1575_4199_143 = std::sqrt(27286875.0 / 1356277.0);
    const auto fs_1575_4199_15015 = std::sqrt(2865121875.0 / 1356277.0);
    const auto fs_1575_442_273 = std::sqrt(52093125.0 / 15028.0);
    const auto fs_1575_442_66 = std::sqrt(81860625.0 / 97682.0);
    const auto fs_1575_44_154 = std::sqrt(17364375.0 / 88.0);
    const auto fs_1575_4_7 = std::sqrt(17364375.0 / 16.0);
    const auto fs_1575_8398_165 = std::sqrt(409303125.0 / 70526404.0);
    const auto fs_1575_8398_20995 = std::sqrt(12403125.0 / 16796.0);
    const auto fs_1575_8398_663 = std::sqrt(7441875.0 / 319124.0);
    const auto fs_1575_8398_92378 = std::sqrt(27286875.0 / 8398.0);
    const auto fs_15_2717_462 = std::sqrt(9450.0 / 671099.0);
    const auto fs_16065_16_30 = std::sqrt(3871263375.0 / 128.0);
    const auto fs_16065_16_70 = std::sqrt(9032947875.0 / 128.0);
    const auto fs_16065_32_22 = std::sqrt(2838926475.0 / 512.0);
    const auto fs_16065_32_30 = std::sqrt(3871263375.0 / 512.0);
    const auto fs_16065_32_55 = std::sqrt(14194632375.0 / 1024.0);
    const auto fs_16065_32_7 = std::sqrt(1806589575.0 / 1024.0);
    const auto fs_166725_176_10 = std::sqrt(138986128125.0 / 15488.0);
    const auto fs_1680_3553_210 = std::sqrt(592704000.0 / 12623809.0);
    const auto fs_1680_3553_385 = std::sqrt(98784000.0 / 1147619.0);
    const auto fs_1680_3553_462 = std::sqrt(118540800.0 / 1147619.0);
    const auto fs_168_2431_5 = std::sqrt(141120.0 / 5909761.0);
    const auto fs_16905_286_22 = std::sqrt(285779025.0 / 3718.0);
    const auto fs_16905_572_3 = std::sqrt(857337075.0 / 327184.0);
    const auto fs_17325_16796_1430 = std::sqrt(16508559375.0 / 10850216.0);
    const auto fs_17325_16796_286 = std::sqrt(3301711875.0 / 10850216.0);
    const auto fs_17325_96577_1430 = std::sqrt(33017118750.0 / 717470533.0);
    const auto fs_17325_96577_286 = std::sqrt(6603423750.0 / 717470533.0);
    const auto fs_1764_96577_165 = std::sqrt(513429840.0 / 9327116929.0);
    const auto fs_1764_96577_39 = std::sqrt(9335088.0 / 717470533.0);
    const auto fs_1782_96577_138567 = std::sqrt(104792292.0 / 2221271.0);
    const auto fs_1782_96577_176358 = std::sqrt(133372008.0 / 2221271.0);
    const auto fs_1782_96577_30030 = std::sqrt(7335460440.0 / 717470533.0);
    const auto fs_1782_96577_455 = std::sqrt(111143340.0 / 717470533.0);
    const auto fs_1782_96577_4641 = std::sqrt(66686004.0 / 42204149.0);
    const auto fs_1782_96577_51051 = std::sqrt(733546044.0 / 42204149.0);
    const auto fs_17955_16_14 = std::sqrt(2256674175.0 / 128.0);
    const auto fs_17955_16_3 = std::sqrt(967146075.0 / 256.0);
    const auto fs_17955_8_3 = std::sqrt(967146075.0 / 64.0);
    const auto fs_17955_8_5 = std::sqrt(1611910125.0 / 64.0);
    const auto fs_1815_96577_78 = std::sqrt(19765350.0 / 717470533.0);
    const auto fs_1848_96577_286 = std::sqrt(75132288.0 / 717470533.0);
    const auto fs_18711_96577_1430 = std::sqrt(38511167310.0 / 717470533.0);
    const auto fs_18900_96577_1155 = std::sqrt(412577550000.0 / 9327116929.0);
    const auto fs_18900_96577_154 = std::sqrt(55010340000.0 / 9327116929.0);
    const auto fs_18900_96577_78 = std::sqrt(2143260000.0 / 717470533.0);
    const auto fs_1890_187_21 = std::sqrt(75014100.0 / 34969.0);
    const auto fs_189_96577_858 = std::sqrt(2357586.0 / 717470533.0);
    const auto fs_18_2431_55 = std::sqrt(1620.0 / 537251.0);
    const auto fs_18_2431_77 = std::sqrt(2268.0 / 537251.0);
    const auto fs_19125_286_3 = std::sqrt(1097296875.0 / 81796.0);
    const auto fs_19125_8_3 = std::sqrt(1097296875.0 / 64.0);
    const auto fs_195615_10868_10 = std::sqrt(191326141125.0 / 59056712.0);
    const auto fs_19665_32_6 = std::sqrt(1160136675.0 / 512.0);
    const auto fs_198_96577_20995 = std::sqrt(196020.0 / 2221271.0);
    const auto fs_198_96577_36465 = std::sqrt(6468660.0 / 42204149.0);
    const auto fs_198_96577_58786 = std::sqrt(548856.0 / 2221271.0);
    const auto fs_198_96577_6006 = std::sqrt(18112248.0 / 717470533.0);
    const auto fs_1995_143_7 = std::sqrt(27860175.0 / 20449.0);
    const auto fs_1995_2431_7 = std::sqrt(27860175.0 / 5909761.0);
    const auto fs_1995_286_462 = std::sqrt(83580525.0 / 3718.0);
    const auto fs_1995_4862_462 = std::sqrt(83580525.0 / 1074502.0);
    const auto fs_204_143_30 = std::sqrt(1248480.0 / 20449.0);
    const auto fs_204_143_70 = std::sqrt(2913120.0 / 20449.0);
    const auto fs_20655_16_66 = std::sqrt(14078757825.0 / 128.0);
    const auto fs_20655_8_10 = std::sqrt(2133145125.0 / 32.0);
    const auto fs_2070_2717_70 = std::sqrt(299943000.0 / 7382089.0);
    const auto fs_2070_2717_77 = std::sqrt(29994300.0 / 671099.0);
    const auto fs_20_247_7 = std::sqrt(2800.0 / 61009.0);
    const auto fs_20_2717_858 = std::sqrt(2400.0 / 51623.0);
    const auto fs_20_323_7 = std::sqrt(2800.0 / 104329.0);
    const auto fs_20_3553_154 = std::sqrt(5600.0 / 1147619.0);
    const auto fs_2100_3553_70 = std::sqrt(308700000.0 / 12623809.0);
    const auto fs_210_11_154 = std::sqrt(617400.0 / 11.0);
    const auto fs_210_1_7 = std::sqrt(308700.0);
    const auto fs_210_2717_11 = std::sqrt(44100.0 / 671099.0);
    const auto fs_21375_286_210 = std::sqrt(47973515625.0 / 40898.0);
    const auto fs_21375_88_210 = std::sqrt(47973515625.0 / 3872.0);
    const auto fs_21735_10868_462 = std::sqrt(9920614725.0 / 5368792.0);
    const auto fs_21735_2717_70 = std::sqrt(33068715750.0 / 7382089.0);
    const auto fs_21735_2717_77 = std::sqrt(3306871575.0 / 671099.0);
    const auto fs_21735_5434_462 = std::sqrt(9920614725.0 / 1342198.0);
    const auto fs_21735_5434_715 = std::sqrt(2362051125.0 / 206492.0);
    const auto fs_21735_572_30 = std::sqrt(7086153375.0 / 163592.0);
    const auto fs_2178_96577_195 = std::sqrt(71155260.0 / 717470533.0);
    const auto fs_21_2431_462 = std::sqrt(18522.0 / 537251.0);
    const auto fs_22050_96577_429 = std::sqrt(16044682500.0 / 717470533.0);
    const auto fs_22050_96577_442 = std::sqrt(972405000.0 / 42204149.0);
    const auto fs_2205_96577_130 = std::sqrt(48620250.0 / 717470533.0);
    const auto fs_2268_96577_221 = std::sqrt(5143824.0 / 42204149.0);
    const auto fs_240_3553_7 = std::sqrt(403200.0 / 12623809.0);
    const auto fs_2415_104_130 = std::sqrt(29161125.0 / 416.0);
    const auto fs_2415_104_182 = std::sqrt(40825575.0 / 416.0);
    const auto fs_2415_104_42 = std::sqrt(122476725.0 / 5408.0);
    const auto fs_2415_1144_770 = std::sqrt(204127875.0 / 59488.0);
    const auto fs_2415_143_21 = std::sqrt(122476725.0 / 20449.0);
    const auto fs_2415_26_7 = std::sqrt(40825575.0 / 676.0);
    const auto fs_2415_2717_3 = std::sqrt(17496675.0 / 7382089.0);
    const auto fs_2415_286_858 = std::sqrt(17496675.0 / 286.0);
    const auto fs_2415_52_13 = std::sqrt(5832225.0 / 208.0);
    const auto fs_2415_52_21 = std::sqrt(122476725.0 / 2704.0);
    const auto fs_2415_52_3 = std::sqrt(17496675.0 / 2704.0);
    const auto fs_2415_52_35 = std::sqrt(204127875.0 / 2704.0);
    const auto fs_2415_52_39 = std::sqrt(17496675.0 / 208.0);
    const auto fs_2415_52_91 = std::sqrt(40825575.0 / 208.0);
    const auto fs_2415_572_1001 = std::sqrt(40825575.0 / 2288.0);
    const auto fs_2415_572_1155 = std::sqrt(612383625.0 / 29744.0);
    const auto fs_2415_572_3003 = std::sqrt(122476725.0 / 2288.0);
    const auto fs_2415_572_385 = std::sqrt(204127875.0 / 29744.0);
    const auto fs_2415_572_715 = std::sqrt(29161125.0 / 2288.0);
    const auto fs_2415_572_858 = std::sqrt(17496675.0 / 1144.0);
    const auto fs_24948_96577_286 = std::sqrt(13692859488.0 / 717470533.0);
    const auto fs_2520_11_7 = std::sqrt(44452800.0 / 121.0);
    const auto fs_2520_3553_21 = std::sqrt(133358400.0 / 12623809.0);
    const auto fs_252_2431_3 = std::sqrt(190512.0 / 5909761.0);
    const auto fs_252_2431_5 = std::sqrt(317520.0 / 5909761.0);
    const auto fs_252_7429_65 = std::sqrt(4127760.0 / 55190041.0);
    const auto fs_252_96577_12155 = std::sqrt(3492720.0 / 42204149.0);
    const auto fs_252_96577_12597 = std::sqrt(190512.0 / 2221271.0);
    const auto fs_252_96577_143 = std::sqrt(698544.0 / 717470533.0);
    const auto fs_252_96577_15015 = std::sqrt(73347120.0 / 717470533.0);
    const auto fs_2550_11_30 = std::sqrt(195075000.0 / 121.0);
    const auto fs_2550_11_70 = std::sqrt(455175000.0 / 121.0);
    const auto fs_25650_143_2 = std::sqrt(1315845000.0 / 20449.0);
    const auto fs_25650_143_33 = std::sqrt(1973767500.0 / 1859.0);
    const auto fs_2565_16_55 = std::sqrt(361857375.0 / 256.0);
    const auto fs_2565_16_77 = std::sqrt(506600325.0 / 256.0);
    const auto fs_2565_286_154 = std::sqrt(46054575.0 / 3718.0);
    const auto fs_2565_286_330 = std::sqrt(98688375.0 / 3718.0);
    const auto fs_2565_4862_154 = std::sqrt(46054575.0 / 1074502.0);
    const auto fs_2565_4862_330 = std::sqrt(98688375.0 / 1074502.0);
    const auto fs_2565_4_2 = std::sqrt(6579225.0 / 8.0);
    const auto fs_2565_4_33 = std::sqrt(217114425.0 / 16.0);
    const auto fs_2583_96577_66 = std::sqrt(440344674.0 / 9327116929.0);
    const auto fs_2646_96577_39 = std::sqrt(21003948.0 / 717470533.0);
    const auto fs_264_96577_20995 = std::sqrt(348480.0 / 2221271.0);
    const auto fs_2673_193154_10010 = std::sqrt(2750797665.0 / 1434941066.0);
    const auto fs_2673_193154_2730 = std::sqrt(750217545.0 / 1434941066.0);
    const auto fs_2673_96577_20995 = std::sqrt(35724645.0 / 2221271.0);
    const auto fs_2673_96577_36465 = std::sqrt(1178913285.0 / 42204149.0);
    const auto fs_2673_96577_58786 = std::sqrt(100029006.0 / 2221271.0);
    const auto fs_2673_96577_6006 = std::sqrt(3300957198.0 / 717470533.0);
    const auto fs_26775_16_10 = std::sqrt(3584503125.0 / 128.0);
    const auto fs_26775_286_3 = std::sqrt(2150701875.0 / 81796.0);
    const auto fs_26775_572_10 = std::sqrt(3584503125.0 / 163592.0);
    const auto fs_26775_8_3 = std::sqrt(2150701875.0 / 64.0);
    const auto fs_2709_96577_10 = std::sqrt(73386810.0 / 9327116929.0);
    const auto fs_27_2431_154 = std::sqrt(10206.0 / 537251.0);
    const auto fs_27_2431_330 = std::sqrt(21870.0 / 537251.0);
    const auto fs_28215_32_30 = std::sqrt(11941293375.0 / 512.0);
    const auto fs_28350_96577_10 = std::sqrt(8037225000.0 / 9327116929.0);
    const auto fs_285_13_42 = std::sqrt(3411450.0 / 169.0);
    const auto fs_285_143_2310 = std::sqrt(17057250.0 / 1859.0);
    const auto fs_285_221_42 = std::sqrt(3411450.0 / 48841.0);
    const auto fs_285_2431_2310 = std::sqrt(17057250.0 / 537251.0);
    const auto fs_29403_96577_195 = std::sqrt(12968046135.0 / 717470533.0);
    const auto fs_29925_143_35 = std::sqrt(31342696875.0 / 20449.0);
    const auto fs_29925_176_462 = std::sqrt(18805618125.0 / 1408.0);
    const auto fs_29925_22_5 = std::sqrt(4477528125.0 / 484.0);
    const auto fs_29925_286_7 = std::sqrt(6268539375.0 / 81796.0);
    const auto fs_29925_44_35 = std::sqrt(31342696875.0 / 1936.0);
    const auto fs_29925_572_462 = std::sqrt(18805618125.0 / 14872.0);
    const auto fs_29925_88_7 = std::sqrt(6268539375.0 / 7744.0);
    const auto fs_2_143_22 = std::sqrt(8.0 / 1859.0);
    const auto fs_2_143_30 = std::sqrt(120.0 / 20449.0);
    const auto fs_2_143_55 = std::sqrt(20.0 / 1859.0);
    const auto fs_2_143_7 = std::sqrt(28.0 / 20449.0);
    const auto fs_306_143_10 = std::sqrt(936360.0 / 20449.0);
    const auto fs_30_2431_210 = std::sqrt(189000.0 / 5909761.0);
    const auto fs_30_2717_462 = std::sqrt(37800.0 / 671099.0);
    const auto fs_30_2717_715 = std::sqrt(4500.0 / 51623.0);
    const auto fs_30_3553_462 = std::sqrt(37800.0 / 1147619.0);
    const auto fs_3105_2717_30 = std::sqrt(289230750.0 / 7382089.0);
    const auto fs_3150_4199_2145 = std::sqrt(1637212500.0 / 1356277.0);
    const auto fs_3150_4199_3315 = std::sqrt(148837500.0 / 79781.0);
    const auto fs_3150_5083_273 = std::sqrt(208372500.0 / 1987453.0);
    const auto fs_3150_5083_66 = std::sqrt(654885000.0 / 25836889.0);
    const auto fs_3150_96577_165 = std::sqrt(1637212500.0 / 9327116929.0);
    const auto fs_3150_96577_20995 = std::sqrt(49612500.0 / 2221271.0);
    const auto fs_3150_96577_663 = std::sqrt(29767500.0 / 42204149.0);
    const auto fs_3150_96577_92378 = std::sqrt(218295000.0 / 2221271.0);
    const auto fs_315_11_462 = std::sqrt(4167450.0 / 11.0);
    const auto fs_33075_8398_39 = std::sqrt(3281866875.0 / 5425108.0);
    const auto fs_330_96577_17017 = std::sqrt(8385300.0 / 42204149.0);
    const auto fs_33345_32_10 = std::sqrt(5559445125.0 / 512.0);
    const auto fs_33_96577_182 = std::sqrt(15246.0 / 717470533.0);
    const auto fs_33_96577_26 = std::sqrt(2178.0 / 717470533.0);
    const auto fs_33_96577_910 = std::sqrt(76230.0 / 717470533.0);
    const auto fs_3402_96577_55 = std::sqrt(636548220.0 / 9327116929.0);
    const auto fs_3420_143_2 = std::sqrt(23392800.0 / 20449.0);
    const auto fs_3420_143_33 = std::sqrt(35089200.0 / 1859.0);
    const auto fs_3420_2431_2 = std::sqrt(23392800.0 / 5909761.0);
    const auto fs_3420_2431_33 = std::sqrt(35089200.0 / 537251.0);
    const auto fs_34425_8_3 = std::sqrt(3555241875.0 / 64.0);
    const auto fs_3450_2717_21 = std::sqrt(249952500.0 / 7382089.0);
    const auto fs_345_143_7 = std::sqrt(833175.0 / 20449.0);
    const auto fs_345_247_13 = std::sqrt(119025.0 / 4693.0);
    const auto fs_345_247_21 = std::sqrt(2499525.0 / 61009.0);
    const auto fs_345_247_3 = std::sqrt(357075.0 / 61009.0);
    const auto fs_345_247_35 = std::sqrt(4165875.0 / 61009.0);
    const auto fs_345_247_39 = std::sqrt(357075.0 / 4693.0);
    const auto fs_345_247_91 = std::sqrt(833175.0 / 4693.0);
    const auto fs_345_2717_1001 = std::sqrt(833175.0 / 51623.0);
    const auto fs_345_2717_1155 = std::sqrt(12497625.0 / 671099.0);
    const auto fs_345_2717_3003 = std::sqrt(2499525.0 / 51623.0);
    const auto fs_345_2717_385 = std::sqrt(4165875.0 / 671099.0);
    const auto fs_345_2717_715 = std::sqrt(595125.0 / 51623.0);
    const auto fs_345_2717_858 = std::sqrt(714150.0 / 51623.0);
    const auto fs_345_494_130 = std::sqrt(595125.0 / 9386.0);
    const auto fs_345_494_182 = std::sqrt(833175.0 / 9386.0);
    const auto fs_345_494_42 = std::sqrt(2499525.0 / 122018.0);
    const auto fs_345_5434_770 = std::sqrt(4165875.0 / 1342198.0);
    const auto fs_3564_96577_20995 = std::sqrt(63510480.0 / 2221271.0);
    const auto fs_357_143_10 = std::sqrt(1274490.0 / 20449.0);
    const auto fs_36225_2717_21 = std::sqrt(27557263125.0 / 7382089.0);
    const auto fs_37800_96577_22 = std::sqrt(31434480000.0 / 9327116929.0);
    const auto fs_37800_96577_455 = std::sqrt(50009400000.0 / 717470533.0);
    const auto fs_378_96577_1001 = std::sqrt(11002068.0 / 717470533.0);
    const auto fs_378_96577_2431 = std::sqrt(1571724.0 / 42204149.0);
    const auto fs_378_96577_4290 = std::sqrt(47151720.0 / 717470533.0);
    const auto fs_378_96577_8398 = std::sqrt(285768.0 / 2221271.0);
    const auto fs_3825_11_10 = std::sqrt(146306250.0 / 121.0);
    const auto fs_3825_143_30 = std::sqrt(438918750.0 / 20449.0);
    const auto fs_3825_143_70 = std::sqrt(1024143750.0 / 20449.0);
    const auto fs_3825_22_66 = std::sqrt(43891875.0 / 22.0);
    const auto fs_3825_286_22 = std::sqrt(14630625.0 / 3718.0);
    const auto fs_3825_286_30 = std::sqrt(219459375.0 / 40898.0);
    const auto fs_3825_286_55 = std::sqrt(73153125.0 / 7436.0);
    const auto fs_3825_286_7 = std::sqrt(102414375.0 / 81796.0);
    const auto fs_3825_4_30 = std::sqrt(219459375.0 / 8.0);
    const auto fs_3825_4_70 = std::sqrt(512071875.0 / 8.0);
    const auto fs_3825_8_22 = std::sqrt(160936875.0 / 32.0);
    const auto fs_3825_8_30 = std::sqrt(219459375.0 / 32.0);
    const auto fs_3825_8_55 = std::sqrt(804684375.0 / 64.0);
    const auto fs_3825_8_7 = std::sqrt(102414375.0 / 64.0);
    const auto fs_38475_176_154 = std::sqrt(10362279375.0 / 1408.0);
    const auto fs_38475_176_330 = std::sqrt(22204884375.0 / 1408.0);
    const auto fs_38475_572_154 = std::sqrt(10362279375.0 / 14872.0);
    const auto fs_38475_572_330 = std::sqrt(22204884375.0 / 14872.0);
    const auto fs_3906_96577_26 = std::sqrt(30513672.0 / 717470533.0);
    const auto fs_3990_143_35 = std::sqrt(557203500.0 / 20449.0);
    const auto fs_3990_2431_35 = std::sqrt(557203500.0 / 5909761.0);
    const auto fs_3_143_66 = std::sqrt(54.0 / 1859.0);
    const auto fs_40_2717_21 = std::sqrt(33600.0 / 7382089.0);
    const auto fs_40_3553_210 = std::sqrt(336000.0 / 12623809.0);
    const auto fs_40_3553_385 = std::sqrt(56000.0 / 1147619.0);
    const auto fs_40_3553_462 = std::sqrt(67200.0 / 1147619.0);
    const auto fs_4140_2717_22 = std::sqrt(34279200.0 / 671099.0);
    const auto fs_420_11_210 = std::sqrt(37044000.0 / 121.0);
    const auto fs_420_11_385 = std::sqrt(6174000.0 / 11.0);
    const auto fs_420_11_462 = std::sqrt(7408800.0 / 11.0);
    const auto fs_42525_8398_55 = std::sqrt(99460659375.0 / 70526404.0);
    const auto fs_4275_143_210 = std::sqrt(3837881250.0 / 20449.0);
    const auto fs_4275_16_210 = std::sqrt(1918940625.0 / 128.0);
    const auto fs_4275_26_42 = std::sqrt(383788125.0 / 338.0);
    const auto fs_4275_286_2310 = std::sqrt(1918940625.0 / 3718.0);
    const auto fs_4275_44_210 = std::sqrt(1918940625.0 / 968.0);
    const auto fs_4275_88_2310 = std::sqrt(1918940625.0 / 352.0);
    const auto fs_4275_8_42 = std::sqrt(383788125.0 / 32.0);
    const auto fs_42_2431_7 = std::sqrt(12348.0 / 5909761.0);
    const auto fs_43470_2717_22 = std::sqrt(3779281800.0 / 671099.0);
    const auto fs_4356_96577_182 = std::sqrt(265646304.0 / 717470533.0);
    const auto fs_44100_96577_165 = std::sqrt(320893650000.0 / 9327116929.0);
    const auto fs_44100_96577_39 = std::sqrt(5834430000.0 / 717470533.0);
    const auto fs_4455_96577_17017 = std::sqrt(1528220925.0 / 42204149.0);
    const auto fs_45885_572_7 = std::sqrt(14738032575.0 / 327184.0);
    const auto fs_4725_11_7 = std::sqrt(156279375.0 / 121.0);
    const auto fs_4725_16796_858 = std::sqrt(736745625.0 / 10850216.0);
    const auto fs_4725_4199_1155 = std::sqrt(25786096875.0 / 17631601.0);
    const auto fs_4725_4199_154 = std::sqrt(3438146250.0 / 17631601.0);
    const auto fs_4725_4199_78 = std::sqrt(133953750.0 / 1356277.0);
    const auto fs_4725_44_21 = std::sqrt(468838125.0 / 1936.0);
    const auto fs_4725_8398_1001 = std::sqrt(1719073125.0 / 5425108.0);
    const auto fs_4725_8398_2431 = std::sqrt(245581875.0 / 319124.0);
    const auto fs_4725_8398_4290 = std::sqrt(3683728125.0 / 2712554.0);
    const auto fs_4725_8398_8398 = std::sqrt(22325625.0 / 8398.0);
    const auto fs_4725_88_462 = std::sqrt(468838125.0 / 352.0);
    const auto fs_4725_96577_858 = std::sqrt(1473491250.0 / 717470533.0);
    const auto fs_48195_16_10 = std::sqrt(11613790125.0 / 128.0);
    const auto fs_48195_32_10 = std::sqrt(11613790125.0 / 512.0);
    const auto fs_48195_64_66 = std::sqrt(76651014825.0 / 2048.0);
    const auto fs_48195_8_3 = std::sqrt(6968274075.0 / 64.0);
    const auto fs_4830_2717_22 = std::sqrt(46657800.0 / 671099.0);
    const auto fs_48825_8398_26 = std::sqrt(2383880625.0 / 2712554.0);
    const auto fs_49005_193154_78 = std::sqrt(7204470075.0 / 1434941066.0);
    const auto fs_4_143_30 = std::sqrt(480.0 / 20449.0);
    const auto fs_4_143_70 = std::sqrt(1120.0 / 20449.0);
    const auto fs_504_96577_2145 = std::sqrt(41912640.0 / 717470533.0);
    const auto fs_504_96577_3315 = std::sqrt(3810240.0 / 42204149.0);
    const auto fs_50715_2717_22 = std::sqrt(5144022450.0 / 671099.0);
    const auto fs_50715_5434_3 = std::sqrt(7716033675.0 / 29528356.0);
    const auto fs_50715_572_11 = std::sqrt(2572011225.0 / 29744.0);
    const auto fs_50_3553_70 = std::sqrt(175000.0 / 12623809.0);
    const auto fs_510_143_3 = std::sqrt(780300.0 / 20449.0);
    const auto fs_525_11_70 = std::sqrt(19293750.0 / 121.0);
    const auto fs_528_96577_273 = std::sqrt(5854464.0 / 717470533.0);
    const auto fs_528_96577_5005 = std::sqrt(107331840.0 / 717470533.0);
    const auto fs_5418_96577_22 = std::sqrt(645803928.0 / 9327116929.0);
    const auto fs_55125_16796_130 = std::sqrt(15193828125.0 / 10850216.0);
    const auto fs_55125_96577_130 = std::sqrt(30387656250.0 / 717470533.0);
    const auto fs_56700_96577_221 = std::sqrt(3214890000.0 / 42204149.0);
    const auto fs_570_143_210 = std::sqrt(68229000.0 / 20449.0);
    const auto fs_570_2431_210 = std::sqrt(68229000.0 / 5909761.0);
    const auto fs_58806_96577_182 = std::sqrt(48414038904.0 / 717470533.0);
    const auto fs_594_96577_5005 = std::sqrt(135841860.0 / 717470533.0);
    const auto fs_59850_143_5 = std::sqrt(17910112500.0 / 20449.0);
    const auto fs_5985_143_14 = std::sqrt(501483150.0 / 20449.0);
    const auto fs_5985_143_3 = std::sqrt(107460675.0 / 20449.0);
    const auto fs_5985_16_7 = std::sqrt(250741575.0 / 256.0);
    const auto fs_5985_2431_14 = std::sqrt(501483150.0 / 5909761.0);
    const auto fs_5985_2431_3 = std::sqrt(107460675.0 / 5909761.0);
    const auto fs_5985_32_462 = std::sqrt(8274471975.0 / 512.0);
    const auto fs_5985_4_5 = std::sqrt(179101125.0 / 16.0);
    const auto fs_5985_8_35 = std::sqrt(1253707875.0 / 64.0);
    const auto fs_5_247_130 = std::sqrt(250.0 / 4693.0);
    const auto fs_5_247_182 = std::sqrt(350.0 / 4693.0);
    const auto fs_5_247_42 = std::sqrt(1050.0 / 61009.0);
    const auto fs_5_2717_770 = std::sqrt(1750.0 / 671099.0);
    const auto fs_60_2717_70 = std::sqrt(252000.0 / 7382089.0);
    const auto fs_60_2717_77 = std::sqrt(25200.0 / 671099.0);
    const auto fs_60_3553_21 = std::sqrt(75600.0 / 12623809.0);
    const auto fs_6300_7429_65 = std::sqrt(2579850000.0 / 55190041.0);
    const auto fs_6300_96577_12155 = std::sqrt(2182950000.0 / 42204149.0);
    const auto fs_6300_96577_12597 = std::sqrt(119070000.0 / 2221271.0);
    const auto fs_6300_96577_143 = std::sqrt(436590000.0 / 717470533.0);
    const auto fs_6300_96577_15015 = std::sqrt(45841950000.0 / 717470533.0);
    const auto fs_630_11_21 = std::sqrt(8334900.0 / 121.0);
    const auto fs_630_17_7 = std::sqrt(2778300.0 / 289.0);
    const auto fs_630_187_154 = std::sqrt(5556600.0 / 3179.0);
    const auto fs_630_96577_858 = std::sqrt(26195400.0 / 717470533.0);
    const auto fs_6375_11_3 = std::sqrt(121921875.0 / 121.0);
    const auto fs_64575_16796_66 = std::sqrt(137607710625.0 / 141052808.0);
    const auto fs_64575_96577_66 = std::sqrt(275215421250.0 / 9327116929.0);
    const auto fs_65205_1144_10 = std::sqrt(21258460125.0 / 654368.0);
    const auto fs_65205_5434_30 = std::sqrt(63775380375.0 / 14764178.0);
    const auto fs_6555_286_6 = std::sqrt(128904075.0 / 40898.0);
    const auto fs_6555_4862_6 = std::sqrt(128904075.0 / 11819522.0);
    const auto fs_660_96577_4862 = std::sqrt(9583200.0 / 42204149.0);
    const auto fs_66150_96577_39 = std::sqrt(13127467500.0 / 717470533.0);
    const auto fs_66_96577_1352078 = std::sqrt(60984.0 / 96577.0);
    const auto fs_66_96577_146965 = std::sqrt(152460.0 / 2221271.0);
    const auto fs_66_96577_25194 = std::sqrt(26136.0 / 2221271.0);
    const auto fs_66_96577_3094 = std::sqrt(60984.0 / 42204149.0);
    const auto fs_66_96577_323323 = std::sqrt(335412.0 / 2221271.0);
    const auto fs_66_96577_429 = std::sqrt(143748.0 / 717470533.0);
    const auto fs_66_96577_461890 = std::sqrt(479160.0 / 2221271.0);
    const auto fs_66_96577_62985 = std::sqrt(65340.0 / 2221271.0);
    const auto fs_66_96577_676039 = std::sqrt(30492.0 / 96577.0);
    const auto fs_66_96577_910 = std::sqrt(304920.0 / 717470533.0);
    const auto fs_67725_16796_10 = std::sqrt(22933378125.0 / 141052808.0);
    const auto fs_67725_8398_22 = std::sqrt(50453431875.0 / 35263202.0);
    const auto fs_67725_96577_10 = std::sqrt(45866756250.0 / 9327116929.0);
    const auto fs_6885_4_30 = std::sqrt(711048375.0 / 8.0);
    const auto fs_6885_4_70 = std::sqrt(1659112875.0 / 8.0);
    const auto fs_6885_8_22 = std::sqrt(521435475.0 / 32.0);
    const auto fs_6885_8_30 = std::sqrt(711048375.0 / 32.0);
    const auto fs_6885_8_55 = std::sqrt(2607177375.0 / 64.0);
    const auto fs_6885_8_7 = std::sqrt(331822575.0 / 64.0);
    const auto fs_690_247_7 = std::sqrt(3332700.0 / 61009.0);
    const auto fs_690_2717_858 = std::sqrt(2856600.0 / 51623.0);
    const auto fs_693_96577_1430 = std::sqrt(52827390.0 / 717470533.0);
    const auto fs_693_96577_286 = std::sqrt(10565478.0 / 717470533.0);
    const auto fs_69_2431_6 = std::sqrt(28566.0 / 5909761.0);
    const auto fs_6_143_10 = std::sqrt(360.0 / 20449.0);
    const auto fs_6_221_42 = std::sqrt(1512.0 / 48841.0);
    const auto fs_6_2431_2310 = std::sqrt(7560.0 / 537251.0);
    const auto fs_70_2717_3 = std::sqrt(14700.0 / 7382089.0);
    const auto fs_7128_96577_273 = std::sqrt(1066976064.0 / 717470533.0);
    const auto fs_7128_96577_5005 = std::sqrt(19561227840.0 / 717470533.0);
    const auto fs_714_143_3 = std::sqrt(1529388.0 / 20449.0);
    const auto fs_7245_10868_770 = std::sqrt(1837150875.0 / 5368792.0);
    const auto fs_7245_1144_462 = std::sqrt(1102290525.0 / 59488.0);
    const auto fs_7245_143_22 = std::sqrt(104980050.0 / 1859.0);
    const auto fs_7245_247_7 = std::sqrt(367430175.0 / 61009.0);
    const auto fs_7245_2717_11 = std::sqrt(52490025.0 / 671099.0);
    const auto fs_7245_2717_858 = std::sqrt(314940150.0 / 51623.0);
    const auto fs_7245_286_7 = std::sqrt(367430175.0 / 81796.0);
    const auto fs_7245_286_70 = std::sqrt(1837150875.0 / 40898.0);
    const auto fs_7245_286_77 = std::sqrt(367430175.0 / 7436.0);
    const auto fs_7245_494_13 = std::sqrt(52490025.0 / 18772.0);
    const auto fs_7245_494_21 = std::sqrt(1102290525.0 / 244036.0);
    const auto fs_7245_494_3 = std::sqrt(157470075.0 / 244036.0);
    const auto fs_7245_494_35 = std::sqrt(1837150875.0 / 244036.0);
    const auto fs_7245_494_39 = std::sqrt(157470075.0 / 18772.0);
    const auto fs_7245_494_91 = std::sqrt(367430175.0 / 18772.0);
    const auto fs_7245_5434_1001 = std::sqrt(367430175.0 / 206492.0);
    const auto fs_7245_5434_1155 = std::sqrt(5511452625.0 / 2684396.0);
    const auto fs_7245_5434_3003 = std::sqrt(1102290525.0 / 206492.0);
    const auto fs_7245_5434_385 = std::sqrt(1837150875.0 / 2684396.0);
    const auto fs_7245_5434_715 = std::sqrt(262450125.0 / 206492.0);
    const auto fs_7245_5434_858 = std::sqrt(157470075.0 / 103246.0);
    const auto fs_7245_572_462 = std::sqrt(1102290525.0 / 14872.0);
    const auto fs_7245_572_715 = std::sqrt(262450125.0 / 2288.0);
    const auto fs_7245_988_130 = std::sqrt(262450125.0 / 37544.0);
    const auto fs_7245_988_182 = std::sqrt(367430175.0 / 37544.0);
    const auto fs_7245_988_42 = std::sqrt(1102290525.0 / 488072.0);
    const auto fs_726_96577_65 = std::sqrt(2635380.0 / 717470533.0);
    const auto fs_72_2431_2 = std::sqrt(10368.0 / 5909761.0);
    const auto fs_72_2431_33 = std::sqrt(15552.0 / 537251.0);
    const auto fs_7560_187_7 = std::sqrt(400075200.0 / 34969.0);
    const auto fs_756_96577_1155 = std::sqrt(660124080.0 / 9327116929.0);
    const auto fs_756_96577_154 = std::sqrt(88016544.0 / 9327116929.0);
    const auto fs_756_96577_78 = std::sqrt(3429216.0 / 717470533.0);
    const auto fs_7695_32_154 = std::sqrt(4559402925.0 / 512.0);
    const auto fs_7695_32_330 = std::sqrt(9770149125.0 / 512.0);
    const auto fs_7875_8398_858 = std::sqrt(2046515625.0 / 2712554.0);
    const auto fs_7875_88_70 = std::sqrt(2170546875.0 / 3872.0);
    const auto fs_792_96577_3003 = std::sqrt(144897984.0 / 717470533.0);
    const auto fs_792_96577_663 = std::sqrt(1881792.0 / 42204149.0);
    const auto fs_7980_143_5 = std::sqrt(318402000.0 / 20449.0);
    const auto fs_7980_2431_5 = std::sqrt(318402000.0 / 5909761.0);
    const auto fs_7_143_10 = std::sqrt(490.0 / 20449.0);
    const auto fs_8019_96577_5005 = std::sqrt(24757178985.0 / 717470533.0);
    const auto fs_80325_32_3 = std::sqrt(19356316875.0 / 1024.0);
    const auto fs_840_323_7 = std::sqrt(4939200.0 / 104329.0);
    const auto fs_840_3553_154 = std::sqrt(9878400.0 / 1147619.0);
    const auto fs_84_2431_35 = std::sqrt(246960.0 / 5909761.0);
    const auto fs_85050_96577_55 = std::sqrt(397842637500.0 / 9327116929.0);
    const auto fs_855_143_55 = std::sqrt(3655125.0 / 1859.0);
    const auto fs_855_143_77 = std::sqrt(5117175.0 / 1859.0);
    const auto fs_855_16_2310 = std::sqrt(844333875.0 / 128.0);
    const auto fs_855_22_10 = std::sqrt(3655125.0 / 242.0);
    const auto fs_855_2431_55 = std::sqrt(3655125.0 / 537251.0);
    const auto fs_855_2431_77 = std::sqrt(5117175.0 / 537251.0);
    const auto fs_855_26_30 = std::sqrt(10965375.0 / 338.0);
    const auto fs_855_374_10 = std::sqrt(3655125.0 / 69938.0);
    const auto fs_855_442_30 = std::sqrt(10965375.0 / 97682.0);
    const auto fs_855_8_210 = std::sqrt(76757625.0 / 32.0);
    const auto fs_882_96577_429 = std::sqrt(25671492.0 / 717470533.0);
    const auto fs_882_96577_442 = std::sqrt(1555848.0 / 42204149.0);
    const auto fs_8910_96577_4862 = std::sqrt(1746538200.0 / 42204149.0);
    const auto fs_891_193154_182 = std::sqrt(5557167.0 / 1434941066.0);
    const auto fs_891_193154_26 = std::sqrt(793881.0 / 1434941066.0);
    const auto fs_891_193154_910 = std::sqrt(27785835.0 / 1434941066.0);
    const auto fs_891_96577_1352078 = std::sqrt(11114334.0 / 96577.0);
    const auto fs_891_96577_146965 = std::sqrt(27785835.0 / 2221271.0);
    const auto fs_891_96577_25194 = std::sqrt(4763286.0 / 2221271.0);
    const auto fs_891_96577_3094 = std::sqrt(11114334.0 / 42204149.0);
    const auto fs_891_96577_323323 = std::sqrt(61128837.0 / 2221271.0);
    const auto fs_891_96577_429 = std::sqrt(26198073.0 / 717470533.0);
    const auto fs_891_96577_461890 = std::sqrt(87326910.0 / 2221271.0);
    const auto fs_891_96577_62985 = std::sqrt(11908215.0 / 2221271.0);
    const auto fs_891_96577_676039 = std::sqrt(5557167.0 / 96577.0);
    const auto fs_891_96577_910 = std::sqrt(55571670.0 / 717470533.0);
    const auto fs_8925_11_3 = std::sqrt(238966875.0 / 121.0);
    const auto fs_8925_22_10 = std::sqrt(398278125.0 / 242.0);
    const auto fs_89775_143_3 = std::sqrt(24178651875.0 / 20449.0);
    const auto fs_89775_143_5 = std::sqrt(40297753125.0 / 20449.0);
    const auto fs_89775_286_14 = std::sqrt(56416854375.0 / 40898.0);
    const auto fs_89775_286_3 = std::sqrt(24178651875.0 / 81796.0);
    const auto fs_89775_44_3 = std::sqrt(24178651875.0 / 1936.0);
    const auto fs_89775_44_5 = std::sqrt(40297753125.0 / 1936.0);
    const auto fs_89775_88_14 = std::sqrt(56416854375.0 / 3872.0);
    const auto fs_89775_88_3 = std::sqrt(24178651875.0 / 7744.0);
    const auto fs_90_2717_30 = std::sqrt(243000.0 / 7382089.0);
    const auto fs_924_96577_2145 = std::sqrt(140873040.0 / 717470533.0);
    const auto fs_924_96577_442 = std::sqrt(1707552.0 / 42204149.0);
    const auto fs_9315_5434_10 = std::sqrt(433846125.0 / 14764178.0);
    const auto fs_9405_16_42 = std::sqrt(1857534525.0 / 128.0);
    const auto fs_9450_4199_22 = std::sqrt(1964655000.0 / 17631601.0);
    const auto fs_9450_4199_455 = std::sqrt(3125587500.0 / 1356277.0);
    const auto fs_9450_96577_1001 = std::sqrt(6876292500.0 / 717470533.0);
    const auto fs_9450_96577_2431 = std::sqrt(982327500.0 / 42204149.0);
    const auto fs_9450_96577_4290 = std::sqrt(29469825000.0 / 717470533.0);
    const auto fs_9450_96577_8398 = std::sqrt(178605000.0 / 2221271.0);
    const auto fs_945_187_462 = std::sqrt(37507050.0 / 3179.0);
    const auto fs_97650_96577_26 = std::sqrt(19071045000.0 / 717470533.0);
    const auto fs_9801_96577_65 = std::sqrt(480298005.0 / 717470533.0);
    const auto fs_98325_176_6 = std::sqrt(29003416875.0 / 15488.0);
    const auto fs_98325_572_6 = std::sqrt(29003416875.0 / 163592.0);
    const auto fs_990_96577_2002 = std::sqrt(150935400.0 / 717470533.0);
    const auto fs_99_96577_10010 = std::sqrt(7546770.0 / 717470533.0);
    const auto fs_99_96577_2730 = std::sqrt(2058210.0 / 717470533.0);
    const auto fs_9_187_10 = std::sqrt(810.0 / 34969.0);
    const auto fs_9_221_30 = std::sqrt(2430.0 / 48841.0);

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ph12_p12, ab_2, pc_0 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p12 = ph12_p12[k];

        pc_0[k] = e_0 * f_155925_64 + e_1 * f_176715_16 * h2_0 - e_1 * f_363825_32 * r_2 + e_2 * f_84645_16 * h4_0 - e_2 * f_75735_4 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 + e_3 * f_1575_2 * h6_0 - e_3 * f_38475_8 * r_2 * h4_0 + e_3 * f_42075_4 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 + e_4 * f_2415_52 * h8_0 - e_4 * f_420_1 * r_2 * h6_0 + e_4 * f_38475_26 * r_4 * h4_0 - e_4 * f_2550_1 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 + e_5 * f_4725_4199 * h10_0 - e_5 * f_7245_494 * r_2 * h8_0 + e_5 * f_1260_17 * r_4 * h6_0 - e_5 * f_2565_13 * r_6 * h4_0 + e_5 * f_3825_13 * r_8 * h2_0 - e_5 * f_315_2 * r_10 + e_6 * f_891_96577 * h12_0 - e_6 * fs_891_96577_1352078 * h12_p12 - e_6 * f_18900_96577 * r_2 * h10_0 + e_6 * f_345_247 * r_4 * h8_0 - e_6 * f_1680_323 * r_6 * h6_0 + e_6 * f_2565_221 * r_8 * h4_0 - e_6 * f_204_13 * r_10 * h2_0 + e_6 * f_105_13 * r_12 - e_7 * f_66_96577 * r_2 * h12_0 + e_7 * fs_66_96577_1352078 * r_2 * h12_p12 + e_7 * f_756_96577 * r_4 * h10_0 - e_7 * f_10_247 * r_6 * h8_0 + e_7 * f_40_323 * r_8 * h6_0 - e_7 * f_54_221 * r_10 * h4_0 + e_7 * f_4_13 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ph12_p11, ab_2, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p11 = ph12_p11[k];

        pc_1[k] = - e_1 * f_176715_32 * h2_p1 - e_2 * fs_28215_32_30 * h4_p1 + e_2 * f_75735_8 * r_2 * h2_p1 - e_3 * fs_1575_4_7 * h6_p1 + e_3 * fs_12825_16_30 * r_2 * h4_p1 - e_3 * f_42075_8 * r_4 * h2_p1 - e_4 * fs_2415_52_3 * h8_p1 + e_4 * fs_210_1_7 * r_2 * h6_p1 - e_4 * fs_12825_52_30 * r_4 * h4_p1 + e_4 * f_1275_1 * r_6 * h2_p1 - e_5 * fs_1575_8398_165 * h10_p1 + e_5 * fs_7245_494_3 * r_2 * h8_p1 - e_5 * fs_630_17_7 * r_4 * h6_p1 + e_5 * fs_855_26_30 * r_6 * h4_p1 - e_5 * f_3825_26 * r_8 * h2_p1 - e_6 * fs_891_193154_26 * h12_p1 - e_6 * fs_891_96577_676039 * h12_p11 + e_6 * fs_3150_96577_165 * r_2 * h10_p1 - e_6 * fs_345_247_3 * r_4 * h8_p1 + e_6 * fs_840_323_7 * r_6 * h6_p1 - e_6 * fs_855_442_30 * r_8 * h4_p1 + e_6 * f_102_13 * r_10 * h2_p1 + e_7 * fs_33_96577_26 * r_2 * h12_p1 + e_7 * fs_66_96577_676039 * r_2 * h12_p11 - e_7 * fs_126_96577_165 * r_4 * h10_p1 + e_7 * fs_10_247_3 * r_6 * h8_p1 - e_7 * fs_20_323_7 * r_8 * h6_p1 + e_7 * fs_9_221_30 * r_10 * h4_p1 - e_7 * f_2_13 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph10_p10, ph12_p2, ph12_p10, ab_2, pc_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p10 = ph12_p10[k];

        pc_2[k] = e_1 * fs_16065_32_22 * h2_p2 + e_2 * fs_7695_32_330 * h4_p2 - e_2 * fs_6885_8_22 * r_2 * h2_p2 + e_3 * fs_1575_22_385 * h6_p2 - e_3 * fs_38475_176_330 * r_2 * h4_p2 + e_3 * fs_3825_8_22 * r_4 * h2_p2 + e_4 * fs_2415_572_1155 * h8_p2 - e_4 * fs_420_11_385 * r_2 * h6_p2 + e_4 * fs_38475_572_330 * r_4 * h4_p2 - e_4 * fs_1275_11_22 * r_6 * h2_p2 + e_5 * fs_14175_8398_10 * h10_p2 - e_5 * fs_1575_4199_12597 * h10_p10 - e_5 * fs_7245_5434_1155 * r_2 * h8_p2 + e_5 * fs_1260_187_385 * r_4 * h6_p2 - e_5 * fs_2565_286_330 * r_6 * h4_p2 + e_5 * fs_3825_286_22 * r_8 * h2_p2 + e_6 * fs_891_193154_182 * h12_p2 - e_6 * fs_891_96577_323323 * h12_p10 - e_6 * fs_28350_96577_10 * r_2 * h10_p2 + e_6 * fs_6300_96577_12597 * r_2 * h10_p10 + e_6 * fs_345_2717_1155 * r_4 * h8_p2 - e_6 * fs_1680_3553_385 * r_6 * h6_p2 + e_6 * fs_2565_4862_330 * r_8 * h4_p2 - e_6 * fs_102_143_22 * r_10 * h2_p2 - e_7 * fs_33_96577_182 * r_2 * h12_p2 + e_7 * fs_66_96577_323323 * r_2 * h12_p10 + e_7 * fs_1134_96577_10 * r_4 * h10_p2 - e_7 * fs_252_96577_12597 * r_4 * h10_p10 - e_7 * fs_10_2717_1155 * r_6 * h8_p2 + e_7 * fs_40_3553_385 * r_8 * h6_p2 - e_7 * fs_27_2431_330 * r_10 * h4_p2 + e_7 * fs_2_143_22 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ph10_p9, ph12_p3, ph12_p9, ab_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p9 = ph12_p9[k];

        pc_3[k] = - e_2 * fs_7695_32_154 * h4_p3 - e_3 * fs_1575_22_462 * h6_p3 + e_3 * fs_38475_176_154 * r_2 * h4_p3 - e_4 * fs_2415_52_21 * h8_p3 + e_4 * fs_420_11_462 * r_2 * h6_p3 - e_4 * fs_38475_572_154 * r_4 * h4_p3 - e_5 * fs_4725_4199_78 * h10_p3 - e_5 * fs_4725_8398_8398 * h10_p9 + e_5 * fs_7245_494_21 * r_2 * h8_p3 - e_5 * fs_1260_187_462 * r_4 * h6_p3 + e_5 * fs_2565_286_154 * r_6 * h4_p3 - e_6 * fs_891_193154_910 * h12_p3 - e_6 * fs_891_96577_146965 * h12_p9 + e_6 * fs_18900_96577_78 * r_2 * h10_p3 + e_6 * fs_9450_96577_8398 * r_2 * h10_p9 - e_6 * fs_345_247_21 * r_4 * h8_p3 + e_6 * fs_1680_3553_462 * r_6 * h6_p3 - e_6 * fs_2565_4862_154 * r_8 * h4_p3 + e_7 * fs_33_96577_910 * r_2 * h12_p3 + e_7 * fs_66_96577_146965 * r_2 * h12_p9 - e_7 * fs_756_96577_78 * r_4 * h10_p3 - e_7 * fs_378_96577_8398 * r_4 * h10_p9 + e_7 * fs_10_247_21 * r_6 * h8_p3 - e_7 * fs_40_3553_462 * r_8 * h6_p3 + e_7 * fs_27_2431_154 * r_10 * h4_p3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p4, ph6_p4, ph8_p4, ph8_p8, ph10_p4, ph10_p8, ph12_p4, ph12_p8, ab_2, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p8 = ph12_p8[k];

        pc_4[k] = e_2 * fs_2565_16_77 * h4_p4 + e_3 * fs_1575_22_385 * h6_p4 - e_3 * fs_12825_88_77 * r_2 * h4_p4 + e_4 * fs_2415_52_35 * h8_p4 - e_4 * fs_2415_52_13 * h8_p8 - e_4 * fs_420_11_385 * r_2 * h6_p4 + e_4 * fs_12825_286_77 * r_4 * h4_p4 + e_5 * fs_11025_4199_39 * h10_p4 - e_5 * fs_14175_4199_221 * h10_p8 - e_5 * fs_7245_494_35 * r_2 * h8_p4 + e_5 * fs_7245_494_13 * r_2 * h8_p8 + e_5 * fs_1260_187_385 * r_4 * h6_p4 - e_5 * fs_855_143_77 * r_6 * h4_p4 + e_6 * fs_891_96577_910 * h12_p4 - e_6 * fs_891_96577_62985 * h12_p8 - e_6 * fs_44100_96577_39 * r_2 * h10_p4 + e_6 * fs_56700_96577_221 * r_2 * h10_p8 + e_6 * fs_345_247_35 * r_4 * h8_p4 - e_6 * fs_345_247_13 * r_4 * h8_p8 - e_6 * fs_1680_3553_385 * r_6 * h6_p4 + e_6 * fs_855_2431_77 * r_8 * h4_p4 - e_7 * fs_66_96577_910 * r_2 * h12_p4 + e_7 * fs_66_96577_62985 * r_2 * h12_p8 + e_7 * fs_1764_96577_39 * r_4 * h10_p4 - e_7 * fs_2268_96577_221 * r_4 * h10_p8 - e_7 * fs_10_247_35 * r_6 * h8_p4 + e_7 * fs_10_247_13 * r_6 * h8_p8 + e_7 * fs_40_3553_385 * r_8 * h6_p4 - e_7 * fs_18_2431_77 * r_10 * h4_p4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, pe_7, ph6_m6, ph6_p5, ph8_m6, ph8_p5, ph8_p7, ph10_m6, ph10_p5, ph10_p7, ph12_m6, ph12_p5, ph12_p7, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h6_m6 = ph6_m6[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_p5 = ph12_p5[k];
        const auto h12_p7 = ph12_p7[k];

        pc_5[k] = - e_3 * fs_1575_4_7 * h6_p5 - e_4 * fs_2415_104_182 * h8_p5 - e_4 * fs_2415_104_130 * h8_p7 + e_4 * fs_210_1_7 * r_2 * h6_p5 - e_5 * fs_33075_8398_39 * h10_p5 - e_5 * fs_3150_4199_3315 * h10_p7 + e_5 * fs_7245_988_182 * r_2 * h8_p5 + e_5 * fs_7245_988_130 * r_2 * h8_p7 - e_5 * fs_630_17_7 * r_4 * h6_p5 - e_6 * fs_891_96577_3094 * h12_p5 - e_6 * fs_891_96577_25194 * h12_p7 + e_6 * fs_66150_96577_39 * r_2 * h10_p5 + e_6 * fs_12600_96577_3315 * r_2 * h10_p7 - e_6 * fs_345_494_182 * r_4 * h8_p5 - e_6 * fs_345_494_130 * r_4 * h8_p7 + e_6 * fs_840_323_7 * r_6 * h6_p5 + e_7 * fs_66_96577_3094 * r_2 * h12_p5 + e_7 * fs_66_96577_25194 * r_2 * h12_p7 - e_7 * fs_2646_96577_39 * r_4 * h10_p5 - e_7 * fs_504_96577_3315 * r_4 * h10_p7 + e_7 * fs_5_247_182 * r_6 * h8_p5 + e_7 * fs_5_247_130 * r_6 * h8_p7 - e_7 * fs_20_323_7 * r_8 * h6_p5;

        pc_6[k] = e_3 * f_1575_2 * h6_m6 + e_4 * fs_2415_52_91 * h8_m6 - e_4 * f_420_1 * r_2 * h6_m6 + e_5 * fs_9450_4199_455 * h10_m6 - e_5 * fs_7245_494_91 * r_2 * h8_m6 + e_5 * f_1260_17 * r_4 * h6_m6 + e_6 * fs_1782_96577_4641 * h12_m6 - e_6 * fs_37800_96577_455 * r_2 * h10_m6 + e_6 * fs_345_247_91 * r_4 * h8_m6 - e_6 * f_1680_323 * r_6 * h6_m6 - e_7 * fs_132_96577_4641 * r_2 * h12_m6 + e_7 * fs_1512_96577_455 * r_4 * h10_m6 - e_7 * fs_10_247_91 * r_6 * h8_m6 + e_7 * f_40_323 * r_8 * h6_m6;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, pe_7, ph6_m5, ph8_m7, ph8_m5, ph10_m7, ph10_m5, ph12_m7, ph12_m5, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h6_m5 = ph6_m5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m5 = ph12_m5[k];

        pc_7[k] = - e_3 * fs_1575_4_7 * h6_m5 + e_4 * fs_2415_104_130 * h8_m7 - e_4 * fs_2415_104_182 * h8_m5 + e_4 * fs_210_1_7 * r_2 * h6_m5 + e_5 * fs_3150_4199_3315 * h10_m7 - e_5 * fs_33075_8398_39 * h10_m5 - e_5 * fs_7245_988_130 * r_2 * h8_m7 + e_5 * fs_7245_988_182 * r_2 * h8_m5 - e_5 * fs_630_17_7 * r_4 * h6_m5 + e_6 * fs_891_96577_25194 * h12_m7 - e_6 * fs_891_96577_3094 * h12_m5 - e_6 * fs_12600_96577_3315 * r_2 * h10_m7 + e_6 * fs_66150_96577_39 * r_2 * h10_m5 + e_6 * fs_345_494_130 * r_4 * h8_m7 - e_6 * fs_345_494_182 * r_4 * h8_m5 + e_6 * fs_840_323_7 * r_6 * h6_m5 - e_7 * fs_66_96577_25194 * r_2 * h12_m7 + e_7 * fs_66_96577_3094 * r_2 * h12_m5 + e_7 * fs_504_96577_3315 * r_4 * h10_m7 - e_7 * fs_2646_96577_39 * r_4 * h10_m5 - e_7 * fs_5_247_130 * r_6 * h8_m7 + e_7 * fs_5_247_182 * r_6 * h8_m5 - e_7 * fs_20_323_7 * r_8 * h6_m5;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m4, ph6_m4, ph8_m8, ph8_m4, ph10_m8, ph10_m4, ph12_m8, ph12_m4, ab_2, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m4 = ph12_m4[k];

        pc_8[k] = e_2 * fs_2565_16_77 * h4_m4 + e_3 * fs_1575_22_385 * h6_m4 - e_3 * fs_12825_88_77 * r_2 * h4_m4 + e_4 * fs_2415_52_13 * h8_m8 + e_4 * fs_2415_52_35 * h8_m4 - e_4 * fs_420_11_385 * r_2 * h6_m4 + e_4 * fs_12825_286_77 * r_4 * h4_m4 + e_5 * fs_14175_4199_221 * h10_m8 + e_5 * fs_11025_4199_39 * h10_m4 - e_5 * fs_7245_494_13 * r_2 * h8_m8 - e_5 * fs_7245_494_35 * r_2 * h8_m4 + e_5 * fs_1260_187_385 * r_4 * h6_m4 - e_5 * fs_855_143_77 * r_6 * h4_m4 + e_6 * fs_891_96577_62985 * h12_m8 + e_6 * fs_891_96577_910 * h12_m4 - e_6 * fs_56700_96577_221 * r_2 * h10_m8 - e_6 * fs_44100_96577_39 * r_2 * h10_m4 + e_6 * fs_345_247_13 * r_4 * h8_m8 + e_6 * fs_345_247_35 * r_4 * h8_m4 - e_6 * fs_1680_3553_385 * r_6 * h6_m4 + e_6 * fs_855_2431_77 * r_8 * h4_m4 - e_7 * fs_66_96577_62985 * r_2 * h12_m8 - e_7 * fs_66_96577_910 * r_2 * h12_m4 + e_7 * fs_2268_96577_221 * r_4 * h10_m8 + e_7 * fs_1764_96577_39 * r_4 * h10_m4 - e_7 * fs_10_247_13 * r_6 * h8_m8 - e_7 * fs_10_247_35 * r_6 * h8_m4 + e_7 * fs_40_3553_385 * r_8 * h6_m4 - e_7 * fs_18_2431_77 * r_10 * h4_m4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m3, ph6_m3, ph8_m3, ph10_m9, ph10_m3, ph12_m9, ph12_m3, ab_2, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m3 = ph12_m3[k];

        pc_9[k] = - e_2 * fs_7695_32_154 * h4_m3 - e_3 * fs_1575_22_462 * h6_m3 + e_3 * fs_38475_176_154 * r_2 * h4_m3 - e_4 * fs_2415_52_21 * h8_m3 + e_4 * fs_420_11_462 * r_2 * h6_m3 - e_4 * fs_38475_572_154 * r_4 * h4_m3 + e_5 * fs_4725_8398_8398 * h10_m9 - e_5 * fs_4725_4199_78 * h10_m3 + e_5 * fs_7245_494_21 * r_2 * h8_m3 - e_5 * fs_1260_187_462 * r_4 * h6_m3 + e_5 * fs_2565_286_154 * r_6 * h4_m3 + e_6 * fs_891_96577_146965 * h12_m9 - e_6 * fs_891_193154_910 * h12_m3 - e_6 * fs_9450_96577_8398 * r_2 * h10_m9 + e_6 * fs_18900_96577_78 * r_2 * h10_m3 - e_6 * fs_345_247_21 * r_4 * h8_m3 + e_6 * fs_1680_3553_462 * r_6 * h6_m3 - e_6 * fs_2565_4862_154 * r_8 * h4_m3 - e_7 * fs_66_96577_146965 * r_2 * h12_m9 + e_7 * fs_33_96577_910 * r_2 * h12_m3 + e_7 * fs_378_96577_8398 * r_4 * h10_m9 - e_7 * fs_756_96577_78 * r_4 * h10_m3 + e_7 * fs_10_247_21 * r_6 * h8_m3 - e_7 * fs_40_3553_462 * r_8 * h6_m3 + e_7 * fs_27_2431_154 * r_10 * h4_m3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m10, ph10_m2, ph12_m10, ph12_m2, ab_2, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m2 = ph12_m2[k];

        pc_10[k] = e_1 * fs_16065_32_22 * h2_m2 + e_2 * fs_7695_32_330 * h4_m2 - e_2 * fs_6885_8_22 * r_2 * h2_m2 + e_3 * fs_1575_22_385 * h6_m2 - e_3 * fs_38475_176_330 * r_2 * h4_m2 + e_3 * fs_3825_8_22 * r_4 * h2_m2 + e_4 * fs_2415_572_1155 * h8_m2 - e_4 * fs_420_11_385 * r_2 * h6_m2 + e_4 * fs_38475_572_330 * r_4 * h4_m2 - e_4 * fs_1275_11_22 * r_6 * h2_m2 + e_5 * fs_1575_4199_12597 * h10_m10 + e_5 * fs_14175_8398_10 * h10_m2 - e_5 * fs_7245_5434_1155 * r_2 * h8_m2 + e_5 * fs_1260_187_385 * r_4 * h6_m2 - e_5 * fs_2565_286_330 * r_6 * h4_m2 + e_5 * fs_3825_286_22 * r_8 * h2_m2 + e_6 * fs_891_96577_323323 * h12_m10 + e_6 * fs_891_193154_182 * h12_m2 - e_6 * fs_6300_96577_12597 * r_2 * h10_m10 - e_6 * fs_28350_96577_10 * r_2 * h10_m2 + e_6 * fs_345_2717_1155 * r_4 * h8_m2 - e_6 * fs_1680_3553_385 * r_6 * h6_m2 + e_6 * fs_2565_4862_330 * r_8 * h4_m2 - e_6 * fs_102_143_22 * r_10 * h2_m2 - e_7 * fs_66_96577_323323 * r_2 * h12_m10 - e_7 * fs_33_96577_182 * r_2 * h12_m2 + e_7 * fs_252_96577_12597 * r_4 * h10_m10 + e_7 * fs_1134_96577_10 * r_4 * h10_m2 - e_7 * fs_10_2717_1155 * r_6 * h8_m2 + e_7 * fs_40_3553_385 * r_8 * h6_m2 - e_7 * fs_27_2431_330 * r_10 * h4_m2 + e_7 * fs_2_143_22 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m1, ph12_m12, ph12_m11, ph12_m1, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m12 = ph12_m12[k];
        const auto h12_m11 = ph12_m11[k];
        const auto h12_m1 = ph12_m1[k];

        pc_11[k] = - e_1 * f_176715_32 * h2_m1 - e_2 * fs_28215_32_30 * h4_m1 + e_2 * f_75735_8 * r_2 * h2_m1 - e_3 * fs_1575_4_7 * h6_m1 + e_3 * fs_12825_16_30 * r_2 * h4_m1 - e_3 * f_42075_8 * r_4 * h2_m1 - e_4 * fs_2415_52_3 * h8_m1 + e_4 * fs_210_1_7 * r_2 * h6_m1 - e_4 * fs_12825_52_30 * r_4 * h4_m1 + e_4 * f_1275_1 * r_6 * h2_m1 - e_5 * fs_1575_8398_165 * h10_m1 + e_5 * fs_7245_494_3 * r_2 * h8_m1 - e_5 * fs_630_17_7 * r_4 * h6_m1 + e_5 * fs_855_26_30 * r_6 * h4_m1 - e_5 * f_3825_26 * r_8 * h2_m1 + e_6 * fs_891_96577_676039 * h12_m11 - e_6 * fs_891_193154_26 * h12_m1 + e_6 * fs_3150_96577_165 * r_2 * h10_m1 - e_6 * fs_345_247_3 * r_4 * h8_m1 + e_6 * fs_840_323_7 * r_6 * h6_m1 - e_6 * fs_855_442_30 * r_8 * h4_m1 + e_6 * f_102_13 * r_10 * h2_m1 - e_7 * fs_66_96577_676039 * r_2 * h12_m11 + e_7 * fs_33_96577_26 * r_2 * h12_m1 - e_7 * fs_126_96577_165 * r_4 * h10_m1 + e_7 * fs_10_247_3 * r_6 * h8_m1 - e_7 * fs_20_323_7 * r_8 * h6_m1 + e_7 * fs_9_221_30 * r_10 * h4_m1 - e_7 * f_2_13 * r_12 * h2_m1;

        pc_12[k] = e_6 * fs_891_96577_1352078 * h12_m12 - e_7 * fs_66_96577_1352078 * r_2 * h12_m12;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ph12_0, ph12_p10, ab_2, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p10 = ph12_p10[k];

        pc_13[k] = e_0 * f_155925_64 + e_1 * f_176715_32 * h2_0 - e_1 * f_363825_32 * r_2 - e_2 * f_28215_8 * h4_0 - e_2 * f_75735_8 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 - e_3 * f_7875_4 * h6_0 + e_3 * f_12825_4 * r_2 * h4_0 + e_3 * f_42075_8 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 - e_4 * f_12075_52 * h8_0 + e_4 * f_1050_1 * r_2 * h6_0 - e_4 * f_12825_13 * r_4 * h4_0 - e_4 * f_1275_1 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 - e_5 * f_77175_8398 * h10_0 + e_5 * fs_1575_8398_92378 * h10_p10 + e_5 * f_36225_494 * r_2 * h8_0 - e_5 * f_3150_17 * r_4 * h6_0 + e_5 * f_1710_13 * r_6 * h4_0 + e_5 * f_3825_26 * r_8 * h2_0 - e_5 * f_315_2 * r_10 - e_6 * f_10692_96577 * h12_0 - e_6 * fs_1782_96577_176358 * h12_p10 + e_6 * f_154350_96577 * r_2 * h10_0 - e_6 * fs_3150_96577_92378 * r_2 * h10_p10 - e_6 * f_1725_247 * r_4 * h8_0 + e_6 * f_4200_323 * r_6 * h6_0 - e_6 * f_1710_221 * r_8 * h4_0 - e_6 * f_102_13 * r_10 * h2_0 + e_6 * f_105_13 * r_12 + e_7 * f_792_96577 * r_2 * h12_0 + e_7 * fs_132_96577_176358 * r_2 * h12_p10 - e_7 * f_6174_96577 * r_4 * h10_0 + e_7 * fs_126_96577_92378 * r_4 * h10_p10 + e_7 * f_50_247 * r_6 * h8_0 - e_7 * f_100_323 * r_8 * h6_0 + e_7 * f_36_221 * r_10 * h4_0 + e_7 * f_2_13 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ph12_p1, ph12_p9, ab_2, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p9 = ph12_p9[k];

        pc_14[k] = - e_1 * fs_48195_64_66 * h2_p1 - e_2 * fs_2565_16_55 * h4_p1 + e_2 * fs_20655_16_66 * r_2 * h2_p1 + e_3 * fs_4725_88_462 * h6_p1 + e_3 * fs_12825_88_55 * r_2 * h4_p1 - e_3 * fs_11475_16_66 * r_4 * h2_p1 + e_4 * fs_7245_143_22 * h8_p1 - e_4 * fs_315_11_462 * r_2 * h6_p1 - e_4 * fs_12825_286_55 * r_4 * h4_p1 + e_4 * fs_3825_22_66 * r_6 * h2_p1 + e_5 * fs_67725_16796_10 * h10_p1 + e_5 * fs_1575_8398_20995 * h10_p9 - e_5 * fs_43470_2717_22 * r_2 * h8_p1 + e_5 * fs_945_187_462 * r_4 * h6_p1 + e_5 * fs_855_143_55 * r_6 * h4_p1 - e_5 * fs_11475_572_66 * r_8 * h2_p1 + e_6 * fs_891_96577_429 * h12_p1 - e_6 * fs_2673_96577_58786 * h12_p9 - e_6 * fs_67725_96577_10 * r_2 * h10_p1 - e_6 * fs_3150_96577_20995 * r_2 * h10_p9 + e_6 * fs_4140_2717_22 * r_4 * h8_p1 - e_6 * fs_1260_3553_462 * r_6 * h6_p1 - e_6 * fs_855_2431_55 * r_8 * h4_p1 + e_6 * fs_153_143_66 * r_10 * h2_p1 - e_7 * fs_66_96577_429 * r_2 * h12_p1 + e_7 * fs_198_96577_58786 * r_2 * h12_p9 + e_7 * fs_2709_96577_10 * r_4 * h10_p1 + e_7 * fs_126_96577_20995 * r_4 * h10_p9 - e_7 * fs_120_2717_22 * r_6 * h8_p1 + e_7 * fs_30_3553_462 * r_8 * h6_p1 + e_7 * fs_18_2431_55 * r_10 * h4_p1 - e_7 * fs_3_143_66 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ph12_p2, ph12_p8, ab_2, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p8 = ph12_p8[k];

        pc_15[k] = e_1 * fs_16065_32_55 * h2_p2 + e_2 * fs_2565_4_33 * h4_p2 - e_2 * fs_6885_8_55 * r_2 * h2_p2 - e_3 * fs_1575_44_154 * h6_p2 - e_3 * fs_12825_22_33 * r_2 * h4_p2 + e_3 * fs_3825_8_55 * r_4 * h2_p2 - e_4 * fs_7245_572_462 * h8_p2 + e_4 * fs_2415_52_39 * h8_p8 + e_4 * fs_210_11_154 * r_2 * h6_p2 + e_4 * fs_25650_143_33 * r_4 * h4_p2 - e_4 * fs_1275_11_55 * r_6 * h2_p2 - e_5 * f_174825_8398 * h10_p2 - e_5 * fs_1575_8398_663 * h10_p8 + e_5 * fs_21735_5434_462 * r_2 * h8_p2 - e_5 * fs_7245_494_39 * r_2 * h8_p8 - e_5 * fs_630_187_154 * r_4 * h6_p2 - e_5 * fs_3420_143_33 * r_6 * h4_p2 + e_5 * fs_3825_286_55 * r_8 * h2_p2 - e_6 * fs_1782_96577_455 * h12_p2 - e_6 * fs_3564_96577_20995 * h12_p8 + e_6 * f_349650_96577 * r_2 * h10_p2 + e_6 * fs_3150_96577_663 * r_2 * h10_p8 - e_6 * fs_1035_2717_462 * r_4 * h8_p2 + e_6 * fs_345_247_39 * r_4 * h8_p8 + e_6 * fs_840_3553_154 * r_6 * h6_p2 + e_6 * fs_3420_2431_33 * r_8 * h4_p2 - e_6 * fs_102_143_55 * r_10 * h2_p2 + e_7 * fs_132_96577_455 * r_2 * h12_p2 + e_7 * fs_264_96577_20995 * r_2 * h12_p8 - e_7 * f_13986_96577 * r_4 * h10_p2 - e_7 * fs_126_96577_663 * r_4 * h10_p8 + e_7 * fs_30_2717_462 * r_6 * h8_p2 - e_7 * fs_10_247_39 * r_6 * h8_p8 - e_7 * fs_20_3553_154 * r_8 * h6_p2 - e_7 * fs_72_2431_33 * r_10 * h4_p2 + e_7 * fs_2_143_55 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ph12_p3, ph12_p7, ab_2, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p7 = ph12_p7[k];

        pc_16[k] = - e_2 * fs_5985_32_462 * h4_p3 - e_3 * fs_1575_44_154 * h6_p3 + e_3 * fs_29925_176_462 * r_2 * h4_p3 + e_4 * fs_2415_26_7 * h8_p3 + e_4 * fs_2415_52_39 * h8_p7 + e_4 * fs_210_11_154 * r_2 * h6_p3 - e_4 * fs_29925_572_462 * r_4 * h4_p3 + e_5 * fs_48825_8398_26 * h10_p3 - e_5 * fs_11025_8398_442 * h10_p7 - e_5 * fs_7245_247_7 * r_2 * h8_p3 - e_5 * fs_7245_494_39 * r_2 * h8_p7 - e_5 * fs_630_187_154 * r_4 * h6_p3 + e_5 * fs_1995_286_462 * r_6 * h4_p3 + e_6 * fs_2673_193154_2730 * h12_p3 - e_6 * fs_2673_96577_20995 * h12_p7 - e_6 * fs_97650_96577_26 * r_2 * h10_p3 + e_6 * fs_22050_96577_442 * r_2 * h10_p7 + e_6 * fs_690_247_7 * r_4 * h8_p3 + e_6 * fs_345_247_39 * r_4 * h8_p7 + e_6 * fs_840_3553_154 * r_6 * h6_p3 - e_6 * fs_1995_4862_462 * r_8 * h4_p3 - e_7 * fs_99_96577_2730 * r_2 * h12_p3 + e_7 * fs_198_96577_20995 * r_2 * h12_p7 + e_7 * fs_3906_96577_26 * r_4 * h10_p3 - e_7 * fs_882_96577_442 * r_4 * h10_p7 - e_7 * fs_20_247_7 * r_6 * h8_p3 - e_7 * fs_10_247_39 * r_6 * h8_p7 - e_7 * fs_20_3553_154 * r_8 * h6_p3 + e_7 * fs_21_2431_462 * r_10 * h4_p3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p4, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p6, ph12_p4, ph12_p6, ab_2, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p6 = ph12_p6[k];

        pc_17[k] = e_2 * fs_855_16_2310 * h4_p4 + e_3 * fs_4725_88_462 * h6_p4 + e_3 * fs_1575_4_7 * h6_p6 - e_3 * fs_4275_88_2310 * r_2 * h4_p4 - e_4 * fs_2415_104_42 * h8_p4 + e_4 * fs_2415_52_13 * h8_p6 - e_4 * fs_315_11_462 * r_2 * h6_p4 - e_4 * fs_210_1_7 * r_2 * h6_p6 + e_4 * fs_4275_286_2310 * r_4 * h4_p4 - e_5 * fs_55125_16796_130 * h10_p4 - e_5 * fs_1575_323_65 * h10_p6 + e_5 * fs_7245_988_42 * r_2 * h8_p4 - e_5 * fs_7245_494_13 * r_2 * h8_p6 + e_5 * fs_945_187_462 * r_4 * h6_p4 + e_5 * fs_630_17_7 * r_4 * h6_p6 - e_5 * fs_285_143_2310 * r_6 * h4_p4 - e_6 * fs_7128_96577_273 * h12_p4 - e_6 * fs_10692_96577_663 * h12_p6 + e_6 * fs_55125_96577_130 * r_2 * h10_p4 + e_6 * fs_6300_7429_65 * r_2 * h10_p6 - e_6 * fs_345_494_42 * r_4 * h8_p4 + e_6 * fs_345_247_13 * r_4 * h8_p6 - e_6 * fs_1260_3553_462 * r_6 * h6_p4 - e_6 * fs_840_323_7 * r_6 * h6_p6 + e_6 * fs_285_2431_2310 * r_8 * h4_p4 + e_7 * fs_528_96577_273 * r_2 * h12_p4 + e_7 * fs_792_96577_663 * r_2 * h12_p6 - e_7 * fs_2205_96577_130 * r_4 * h10_p4 - e_7 * fs_252_7429_65 * r_4 * h10_p6 + e_7 * fs_5_247_42 * r_6 * h8_p4 - e_7 * fs_10_247_13 * r_6 * h8_p6 + e_7 * fs_30_3553_462 * r_8 * h6_p4 + e_7 * fs_20_323_7 * r_8 * h6_p6 - e_7 * fs_6_2431_2310 * r_10 * h4_p4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, pe_7, ph6_m5, ph10_m5, ph12_m5, ab_2, pc_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h6_m5 = ph6_m5[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h12_m5 = ph12_m5[k];

        pc_18[k] = - e_3 * f_7875_4 * h6_m5 + e_4 * f_1050_1 * r_2 * h6_m5 + e_5 * fs_1575_442_273 * h10_m5 - e_5 * f_3150_17 * r_4 * h6_m5 + e_6 * fs_12474_96577_442 * h12_m5 - e_6 * fs_3150_5083_273 * r_2 * h10_m5 + e_6 * f_4200_323 * r_6 * h6_m5 - e_7 * fs_924_96577_442 * r_2 * h12_m5 + e_7 * fs_126_5083_273 * r_4 * h10_m5 - e_7 * f_100_323 * r_8 * h6_m5;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ph12_m6, ph12_m4, ab_2, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m4 = ph12_m4[k];

        pc_19[k] = e_2 * fs_855_16_2310 * h4_m4 - e_3 * fs_1575_4_7 * h6_m6 + e_3 * fs_4725_88_462 * h6_m4 - e_3 * fs_4275_88_2310 * r_2 * h4_m4 - e_4 * fs_2415_52_13 * h8_m6 - e_4 * fs_2415_104_42 * h8_m4 + e_4 * fs_210_1_7 * r_2 * h6_m6 - e_4 * fs_315_11_462 * r_2 * h6_m4 + e_4 * fs_4275_286_2310 * r_4 * h4_m4 + e_5 * fs_1575_323_65 * h10_m6 - e_5 * fs_55125_16796_130 * h10_m4 + e_5 * fs_7245_494_13 * r_2 * h8_m6 + e_5 * fs_7245_988_42 * r_2 * h8_m4 - e_5 * fs_630_17_7 * r_4 * h6_m6 + e_5 * fs_945_187_462 * r_4 * h6_m4 - e_5 * fs_285_143_2310 * r_6 * h4_m4 + e_6 * fs_10692_96577_663 * h12_m6 - e_6 * fs_7128_96577_273 * h12_m4 - e_6 * fs_6300_7429_65 * r_2 * h10_m6 + e_6 * fs_55125_96577_130 * r_2 * h10_m4 - e_6 * fs_345_247_13 * r_4 * h8_m6 - e_6 * fs_345_494_42 * r_4 * h8_m4 + e_6 * fs_840_323_7 * r_6 * h6_m6 - e_6 * fs_1260_3553_462 * r_6 * h6_m4 + e_6 * fs_285_2431_2310 * r_8 * h4_m4 - e_7 * fs_792_96577_663 * r_2 * h12_m6 + e_7 * fs_528_96577_273 * r_2 * h12_m4 + e_7 * fs_252_7429_65 * r_4 * h10_m6 - e_7 * fs_2205_96577_130 * r_4 * h10_m4 + e_7 * fs_10_247_13 * r_6 * h8_m6 + e_7 * fs_5_247_42 * r_6 * h8_m4 - e_7 * fs_20_323_7 * r_8 * h6_m6 + e_7 * fs_30_3553_462 * r_8 * h6_m4 - e_7 * fs_6_2431_2310 * r_10 * h4_m4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ph12_m7, ph12_m3, ab_2, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m3 = ph12_m3[k];

        pc_20[k] = - e_2 * fs_5985_32_462 * h4_m3 - e_3 * fs_1575_44_154 * h6_m3 + e_3 * fs_29925_176_462 * r_2 * h4_m3 - e_4 * fs_2415_52_39 * h8_m7 + e_4 * fs_2415_26_7 * h8_m3 + e_4 * fs_210_11_154 * r_2 * h6_m3 - e_4 * fs_29925_572_462 * r_4 * h4_m3 + e_5 * fs_11025_8398_442 * h10_m7 + e_5 * fs_48825_8398_26 * h10_m3 + e_5 * fs_7245_494_39 * r_2 * h8_m7 - e_5 * fs_7245_247_7 * r_2 * h8_m3 - e_5 * fs_630_187_154 * r_4 * h6_m3 + e_5 * fs_1995_286_462 * r_6 * h4_m3 + e_6 * fs_2673_96577_20995 * h12_m7 + e_6 * fs_2673_193154_2730 * h12_m3 - e_6 * fs_22050_96577_442 * r_2 * h10_m7 - e_6 * fs_97650_96577_26 * r_2 * h10_m3 - e_6 * fs_345_247_39 * r_4 * h8_m7 + e_6 * fs_690_247_7 * r_4 * h8_m3 + e_6 * fs_840_3553_154 * r_6 * h6_m3 - e_6 * fs_1995_4862_462 * r_8 * h4_m3 - e_7 * fs_198_96577_20995 * r_2 * h12_m7 - e_7 * fs_99_96577_2730 * r_2 * h12_m3 + e_7 * fs_882_96577_442 * r_4 * h10_m7 + e_7 * fs_3906_96577_26 * r_4 * h10_m3 + e_7 * fs_10_247_39 * r_6 * h8_m7 - e_7 * fs_20_247_7 * r_6 * h8_m3 - e_7 * fs_20_3553_154 * r_8 * h6_m3 + e_7 * fs_21_2431_462 * r_10 * h4_m3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ph12_m8, ph12_m2, ab_2, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m2 = ph12_m2[k];

        pc_21[k] = e_1 * fs_16065_32_55 * h2_m2 + e_2 * fs_2565_4_33 * h4_m2 - e_2 * fs_6885_8_55 * r_2 * h2_m2 - e_3 * fs_1575_44_154 * h6_m2 - e_3 * fs_12825_22_33 * r_2 * h4_m2 + e_3 * fs_3825_8_55 * r_4 * h2_m2 - e_4 * fs_2415_52_39 * h8_m8 - e_4 * fs_7245_572_462 * h8_m2 + e_4 * fs_210_11_154 * r_2 * h6_m2 + e_4 * fs_25650_143_33 * r_4 * h4_m2 - e_4 * fs_1275_11_55 * r_6 * h2_m2 + e_5 * fs_1575_8398_663 * h10_m8 - e_5 * f_174825_8398 * h10_m2 + e_5 * fs_7245_494_39 * r_2 * h8_m8 + e_5 * fs_21735_5434_462 * r_2 * h8_m2 - e_5 * fs_630_187_154 * r_4 * h6_m2 - e_5 * fs_3420_143_33 * r_6 * h4_m2 + e_5 * fs_3825_286_55 * r_8 * h2_m2 + e_6 * fs_3564_96577_20995 * h12_m8 - e_6 * fs_1782_96577_455 * h12_m2 - e_6 * fs_3150_96577_663 * r_2 * h10_m8 + e_6 * f_349650_96577 * r_2 * h10_m2 - e_6 * fs_345_247_39 * r_4 * h8_m8 - e_6 * fs_1035_2717_462 * r_4 * h8_m2 + e_6 * fs_840_3553_154 * r_6 * h6_m2 + e_6 * fs_3420_2431_33 * r_8 * h4_m2 - e_6 * fs_102_143_55 * r_10 * h2_m2 - e_7 * fs_264_96577_20995 * r_2 * h12_m8 + e_7 * fs_132_96577_455 * r_2 * h12_m2 + e_7 * fs_126_96577_663 * r_4 * h10_m8 - e_7 * f_13986_96577 * r_4 * h10_m2 + e_7 * fs_10_247_39 * r_6 * h8_m8 + e_7 * fs_30_2717_462 * r_6 * h8_m2 - e_7 * fs_20_3553_154 * r_8 * h6_m2 - e_7 * fs_72_2431_33 * r_10 * h4_m2 + e_7 * fs_2_143_55 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m10, ph10_m9, ph10_m1, ph12_m10, ph12_m9, ph12_m1, ab_2, pc_22, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m1 = ph12_m1[k];

        pc_22[k] = - e_1 * fs_48195_64_66 * h2_m1 - e_2 * fs_2565_16_55 * h4_m1 + e_2 * fs_20655_16_66 * r_2 * h2_m1 + e_3 * fs_4725_88_462 * h6_m1 + e_3 * fs_12825_88_55 * r_2 * h4_m1 - e_3 * fs_11475_16_66 * r_4 * h2_m1 + e_4 * fs_7245_143_22 * h8_m1 - e_4 * fs_315_11_462 * r_2 * h6_m1 - e_4 * fs_12825_286_55 * r_4 * h4_m1 + e_4 * fs_3825_22_66 * r_6 * h2_m1 - e_5 * fs_1575_8398_20995 * h10_m9 + e_5 * fs_67725_16796_10 * h10_m1 - e_5 * fs_43470_2717_22 * r_2 * h8_m1 + e_5 * fs_945_187_462 * r_4 * h6_m1 + e_5 * fs_855_143_55 * r_6 * h4_m1 - e_5 * fs_11475_572_66 * r_8 * h2_m1 + e_6 * fs_2673_96577_58786 * h12_m9 + e_6 * fs_891_96577_429 * h12_m1 + e_6 * fs_3150_96577_20995 * r_2 * h10_m9 - e_6 * fs_67725_96577_10 * r_2 * h10_m1 + e_6 * fs_4140_2717_22 * r_4 * h8_m1 - e_6 * fs_1260_3553_462 * r_6 * h6_m1 - e_6 * fs_855_2431_55 * r_8 * h4_m1 + e_6 * fs_153_143_66 * r_10 * h2_m1 - e_7 * fs_198_96577_58786 * r_2 * h12_m9 - e_7 * fs_66_96577_429 * r_2 * h12_m1 - e_7 * fs_126_96577_20995 * r_4 * h10_m9 + e_7 * fs_2709_96577_10 * r_4 * h10_m1 - e_7 * fs_120_2717_22 * r_6 * h8_m1 + e_7 * fs_30_3553_462 * r_8 * h6_m1 + e_7 * fs_18_2431_55 * r_10 * h4_m1 - e_7 * fs_3_143_66 * r_12 * h2_m1;

        pc_23[k] = - e_5 * fs_1575_8398_92378 * h10_m10 + e_6 * fs_1782_96577_176358 * h12_m10 + e_6 * fs_3150_96577_92378 * r_2 * h10_m10 - e_7 * fs_132_96577_176358 * r_2 * h12_m10 - e_7 * fs_126_96577_92378 * r_4 * h10_m10;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m1, ph12_m11, ph12_m1, ab_2, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m11 = ph12_m11[k];
        const auto h12_m1 = ph12_m1[k];

        pc_24[k] = e_1 * f_176715_32 * h2_m1 + e_2 * fs_28215_32_30 * h4_m1 - e_2 * f_75735_8 * r_2 * h2_m1 + e_3 * fs_1575_4_7 * h6_m1 - e_3 * fs_12825_16_30 * r_2 * h4_m1 + e_3 * f_42075_8 * r_4 * h2_m1 + e_4 * fs_2415_52_3 * h8_m1 - e_4 * fs_210_1_7 * r_2 * h6_m1 + e_4 * fs_12825_52_30 * r_4 * h4_m1 - e_4 * f_1275_1 * r_6 * h2_m1 + e_5 * fs_1575_8398_165 * h10_m1 - e_5 * fs_7245_494_3 * r_2 * h8_m1 + e_5 * fs_630_17_7 * r_4 * h6_m1 - e_5 * fs_855_26_30 * r_6 * h4_m1 + e_5 * f_3825_26 * r_8 * h2_m1 + e_6 * fs_891_96577_676039 * h12_m11 + e_6 * fs_891_193154_26 * h12_m1 - e_6 * fs_3150_96577_165 * r_2 * h10_m1 + e_6 * fs_345_247_3 * r_4 * h8_m1 - e_6 * fs_840_323_7 * r_6 * h6_m1 + e_6 * fs_855_442_30 * r_8 * h4_m1 - e_6 * f_102_13 * r_10 * h2_m1 - e_7 * fs_66_96577_676039 * r_2 * h12_m11 - e_7 * fs_33_96577_26 * r_2 * h12_m1 + e_7 * fs_126_96577_165 * r_4 * h10_m1 - e_7 * fs_10_247_3 * r_6 * h8_m1 + e_7 * fs_20_323_7 * r_8 * h6_m1 - e_7 * fs_9_221_30 * r_10 * h4_m1 + e_7 * f_2_13 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ph12_0, ph12_p8, ab_2, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p8 = ph12_p8[k];

        pc_25[k] = e_0 * f_155925_64 + e_1 * f_16065_16 * h2_0 - e_1 * f_363825_32 * r_2 - e_2 * f_5130_1 * h4_0 - e_2 * f_6885_4 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 + e_3 * f_3150_11 * h6_0 + e_3 * f_51300_11 * r_2 * h4_0 + e_3 * f_3825_4 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 + e_4 * f_214935_572 * h8_0 - e_4 * fs_7245_572_715 * h8_p8 - e_4 * f_1680_11 * r_2 * h6_0 - e_4 * f_205200_143 * r_4 * h4_0 - e_4 * f_2550_11 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 + e_5 * f_130725_4199 * h10_0 + e_5 * fs_1575_4199_12155 * h10_p8 - e_5 * f_644805_5434 * r_2 * h8_0 + e_5 * fs_21735_5434_715 * r_2 * h8_p8 + e_5 * f_5040_187 * r_4 * h6_0 + e_5 * f_27360_143 * r_6 * h4_0 + e_5 * f_3825_143 * r_8 * h2_0 - e_5 * f_315_2 * r_10 + e_6 * f_58806_96577 * h12_0 - e_6 * fs_1782_96577_138567 * h12_p8 - e_6 * f_522900_96577 * r_2 * h10_0 - e_6 * fs_6300_96577_12155 * r_2 * h10_p8 + e_6 * f_30705_2717 * r_4 * h8_0 - e_6 * fs_1035_2717_715 * r_4 * h8_p8 - e_6 * f_6720_3553 * r_6 * h6_0 - e_6 * f_27360_2431 * r_8 * h4_0 - e_6 * f_204_143 * r_10 * h2_0 + e_6 * f_105_13 * r_12 - e_7 * f_4356_96577 * r_2 * h12_0 + e_7 * fs_132_96577_138567 * r_2 * h12_p8 + e_7 * f_20916_96577 * r_4 * h10_0 + e_7 * fs_252_96577_12155 * r_4 * h10_p8 - e_7 * f_890_2717 * r_6 * h8_0 + e_7 * fs_30_2717_715 * r_6 * h8_p8 + e_7 * f_160_3553 * r_8 * h6_0 + e_7 * f_576_2431 * r_10 * h4_0 + e_7 * f_4_143 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ph12_p1, ph12_p7, ab_2, pc_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p7 = ph12_p7[k];

        pc_26[k] = - e_1 * fs_112455_64_10 * h2_p1 + e_2 * fs_17955_16_3 * h4_p1 + e_2 * fs_48195_16_10 * r_2 * h2_p1 + e_3 * fs_7875_88_70 * h6_p1 - e_3 * fs_89775_88_3 * r_2 * h4_p1 - e_3 * fs_26775_16_10 * r_4 * h2_p1 - e_4 * fs_21735_572_30 * h8_p1 - e_4 * fs_2415_572_858 * h8_p7 - e_4 * fs_525_11_70 * r_2 * h6_p1 + e_4 * fs_89775_286_3 * r_4 * h4_p1 + e_4 * fs_8925_22_10 * r_6 * h2_p1 - e_5 * fs_64575_16796_66 * h10_p1 + e_5 * fs_4725_8398_2431 * h10_p7 + e_5 * fs_65205_5434_30 * r_2 * h8_p1 + e_5 * fs_7245_5434_858 * r_2 * h8_p7 + e_5 * fs_1575_187_70 * r_4 * h6_p1 - e_5 * fs_5985_143_3 * r_6 * h4_p1 - e_5 * fs_26775_572_10 * r_8 * h2_p1 - e_6 * fs_9801_96577_65 * h12_p1 - e_6 * fs_891_96577_461890 * h12_p7 + e_6 * fs_64575_96577_66 * r_2 * h10_p1 - e_6 * fs_9450_96577_2431 * r_2 * h10_p7 - e_6 * fs_3105_2717_30 * r_4 * h8_p1 - e_6 * fs_345_2717_858 * r_4 * h8_p7 - e_6 * fs_2100_3553_70 * r_6 * h6_p1 + e_6 * fs_5985_2431_3 * r_8 * h4_p1 + e_6 * fs_357_143_10 * r_10 * h2_p1 + e_7 * fs_726_96577_65 * r_2 * h12_p1 + e_7 * fs_66_96577_461890 * r_2 * h12_p7 - e_7 * fs_2583_96577_66 * r_4 * h10_p1 + e_7 * fs_378_96577_2431 * r_4 * h10_p7 + e_7 * fs_90_2717_30 * r_6 * h8_p1 + e_7 * fs_10_2717_858 * r_6 * h8_p7 + e_7 * fs_50_3553_70 * r_8 * h6_p1 - e_7 * fs_126_2431_3 * r_10 * h4_p1 - e_7 * fs_7_143_10 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ph12_p2, ph12_p6, ab_2, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p6 = ph12_p6[k];

        pc_27[k] = e_1 * fs_48195_32_10 * h2_p2 + e_2 * fs_19665_32_6 * h4_p2 - e_2 * fs_20655_8_10 * r_2 * h2_p2 - e_3 * fs_4725_11_7 * h6_p2 - e_3 * fs_1575_22_385 * h6_p6 - e_3 * fs_98325_176_6 * r_2 * h4_p2 + e_3 * fs_11475_8_10 * r_4 * h2_p2 + e_4 * fs_2415_143_21 * h8_p2 + e_4 * fs_2415_572_715 * h8_p6 + e_4 * fs_2520_11_7 * r_2 * h6_p2 + e_4 * fs_420_11_385 * r_2 * h6_p6 + e_4 * fs_98325_572_6 * r_4 * h4_p2 - e_4 * fs_3825_11_10 * r_6 * h2_p2 + e_5 * fs_67725_8398_22 * h10_p2 + e_5 * fs_1575_4199_143 * h10_p6 - e_5 * fs_14490_2717_21 * r_2 * h8_p2 - e_5 * fs_7245_5434_715 * r_2 * h8_p6 - e_5 * fs_7560_187_7 * r_4 * h6_p2 - e_5 * fs_1260_187_385 * r_4 * h6_p6 - e_5 * fs_6555_286_6 * r_6 * h4_p2 + e_5 * fs_11475_286_10 * r_8 * h2_p2 + e_6 * fs_2673_193154_10010 * h12_p2 - e_6 * fs_2673_96577_36465 * h12_p6 - e_6 * fs_135450_96577_22 * r_2 * h10_p2 - e_6 * fs_6300_96577_143 * r_2 * h10_p6 + e_6 * fs_1380_2717_21 * r_4 * h8_p2 + e_6 * fs_345_2717_715 * r_4 * h8_p6 + e_6 * fs_10080_3553_7 * r_6 * h6_p2 + e_6 * fs_1680_3553_385 * r_6 * h6_p6 + e_6 * fs_6555_4862_6 * r_8 * h4_p2 - e_6 * fs_306_143_10 * r_10 * h2_p2 - e_7 * fs_99_96577_10010 * r_2 * h12_p2 + e_7 * fs_198_96577_36465 * r_2 * h12_p6 + e_7 * fs_5418_96577_22 * r_4 * h10_p2 + e_7 * fs_252_96577_143 * r_4 * h10_p6 - e_7 * fs_40_2717_21 * r_6 * h8_p2 - e_7 * fs_10_2717_715 * r_6 * h8_p6 - e_7 * fs_240_3553_7 * r_8 * h6_p2 - e_7 * fs_40_3553_385 * r_8 * h6_p6 - e_7 * fs_69_2431_6 * r_10 * h4_p2 + e_7 * fs_6_143_10 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ph12_p3, ph12_p5, ab_2, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p5 = ph12_p5[k];

        pc_28[k] = - e_2 * fs_4275_16_210 * h4_p3 + e_3 * fs_7875_88_70 * h6_p3 - e_3 * fs_4725_88_462 * h6_p5 + e_3 * fs_21375_88_210 * r_2 * h4_p3 + e_4 * fs_2415_572_385 * h8_p3 + e_4 * fs_2415_572_3003 * h8_p5 - e_4 * fs_525_11_70 * r_2 * h6_p3 + e_4 * fs_315_11_462 * r_2 * h6_p5 - e_4 * fs_21375_286_210 * r_4 * h4_p3 - e_5 * fs_17325_16796_1430 * h10_p3 - e_5 * fs_17325_16796_286 * h10_p5 - e_5 * fs_7245_5434_385 * r_2 * h8_p3 - e_5 * fs_7245_5434_3003 * r_2 * h8_p5 + e_5 * fs_1575_187_70 * r_4 * h6_p3 - e_5 * fs_945_187_462 * r_4 * h6_p5 + e_5 * fs_1425_143_210 * r_6 * h4_p3 - e_6 * fs_2673_96577_6006 * h12_p3 - e_6 * fs_1782_96577_51051 * h12_p5 + e_6 * fs_17325_96577_1430 * r_2 * h10_p3 + e_6 * fs_17325_96577_286 * r_2 * h10_p5 + e_6 * fs_345_2717_385 * r_4 * h8_p3 + e_6 * fs_345_2717_3003 * r_4 * h8_p5 - e_6 * fs_2100_3553_70 * r_6 * h6_p3 + e_6 * fs_1260_3553_462 * r_6 * h6_p5 - e_6 * fs_1425_2431_210 * r_8 * h4_p3 + e_7 * fs_198_96577_6006 * r_2 * h12_p3 + e_7 * fs_132_96577_51051 * r_2 * h12_p5 - e_7 * fs_693_96577_1430 * r_4 * h10_p3 - e_7 * fs_693_96577_286 * r_4 * h10_p5 - e_7 * fs_10_2717_385 * r_6 * h8_p3 - e_7 * fs_10_2717_3003 * r_6 * h8_p5 + e_7 * fs_50_3553_70 * r_8 * h6_p3 - e_7 * fs_30_3553_462 * r_8 * h6_p5 + e_7 * fs_30_2431_210 * r_10 * h4_p3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ph12_m4, ab_2, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m4 = ph12_m4[k];

        pc_29[k] = e_2 * fs_17955_8_5 * h4_m4 + e_3 * f_3150_11 * h6_m4 - e_3 * fs_89775_44_5 * r_2 * h4_m4 - e_4 * fs_50715_572_11 * h8_m4 - e_4 * f_1680_11 * r_2 * h6_m4 + e_4 * fs_89775_143_5 * r_4 * h4_m4 + e_5 * fs_1575_4199_15015 * h10_m4 + e_5 * fs_152145_5434_11 * r_2 * h8_m4 + e_5 * f_5040_187 * r_4 * h6_m4 - e_5 * fs_11970_143_5 * r_6 * h4_m4 + e_6 * fs_24948_96577_286 * h12_m4 - e_6 * fs_6300_96577_15015 * r_2 * h10_m4 - e_6 * fs_7245_2717_11 * r_4 * h8_m4 - e_6 * f_6720_3553 * r_6 * h6_m4 + e_6 * fs_11970_2431_5 * r_8 * h4_m4 - e_7 * fs_1848_96577_286 * r_2 * h12_m4 + e_7 * fs_252_96577_15015 * r_4 * h10_m4 + e_7 * fs_210_2717_11 * r_6 * h8_m4 + e_7 * f_160_3553 * r_8 * h6_m4 - e_7 * fs_252_2431_5 * r_10 * h4_m4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ph12_m5, ph12_m3, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m3 = ph12_m3[k];

        pc_30[k] = - e_2 * fs_4275_16_210 * h4_m3 + e_3 * fs_4725_88_462 * h6_m5 + e_3 * fs_7875_88_70 * h6_m3 + e_3 * fs_21375_88_210 * r_2 * h4_m3 - e_4 * fs_2415_572_3003 * h8_m5 + e_4 * fs_2415_572_385 * h8_m3 - e_4 * fs_315_11_462 * r_2 * h6_m5 - e_4 * fs_525_11_70 * r_2 * h6_m3 - e_4 * fs_21375_286_210 * r_4 * h4_m3 + e_5 * fs_17325_16796_286 * h10_m5 - e_5 * fs_17325_16796_1430 * h10_m3 + e_5 * fs_7245_5434_3003 * r_2 * h8_m5 - e_5 * fs_7245_5434_385 * r_2 * h8_m3 + e_5 * fs_945_187_462 * r_4 * h6_m5 + e_5 * fs_1575_187_70 * r_4 * h6_m3 + e_5 * fs_1425_143_210 * r_6 * h4_m3 + e_6 * fs_1782_96577_51051 * h12_m5 - e_6 * fs_2673_96577_6006 * h12_m3 - e_6 * fs_17325_96577_286 * r_2 * h10_m5 + e_6 * fs_17325_96577_1430 * r_2 * h10_m3 - e_6 * fs_345_2717_3003 * r_4 * h8_m5 + e_6 * fs_345_2717_385 * r_4 * h8_m3 - e_6 * fs_1260_3553_462 * r_6 * h6_m5 - e_6 * fs_2100_3553_70 * r_6 * h6_m3 - e_6 * fs_1425_2431_210 * r_8 * h4_m3 - e_7 * fs_132_96577_51051 * r_2 * h12_m5 + e_7 * fs_198_96577_6006 * r_2 * h12_m3 + e_7 * fs_693_96577_286 * r_4 * h10_m5 - e_7 * fs_693_96577_1430 * r_4 * h10_m3 + e_7 * fs_10_2717_3003 * r_6 * h8_m5 - e_7 * fs_10_2717_385 * r_6 * h8_m3 + e_7 * fs_30_3553_462 * r_8 * h6_m5 + e_7 * fs_50_3553_70 * r_8 * h6_m3 + e_7 * fs_30_2431_210 * r_10 * h4_m3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ph12_m6, ph12_m2, ab_2, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m2 = ph12_m2[k];

        pc_31[k] = e_1 * fs_48195_32_10 * h2_m2 + e_2 * fs_19665_32_6 * h4_m2 - e_2 * fs_20655_8_10 * r_2 * h2_m2 + e_3 * fs_1575_22_385 * h6_m6 - e_3 * fs_4725_11_7 * h6_m2 - e_3 * fs_98325_176_6 * r_2 * h4_m2 + e_3 * fs_11475_8_10 * r_4 * h2_m2 - e_4 * fs_2415_572_715 * h8_m6 + e_4 * fs_2415_143_21 * h8_m2 - e_4 * fs_420_11_385 * r_2 * h6_m6 + e_4 * fs_2520_11_7 * r_2 * h6_m2 + e_4 * fs_98325_572_6 * r_4 * h4_m2 - e_4 * fs_3825_11_10 * r_6 * h2_m2 - e_5 * fs_1575_4199_143 * h10_m6 + e_5 * fs_67725_8398_22 * h10_m2 + e_5 * fs_7245_5434_715 * r_2 * h8_m6 - e_5 * fs_14490_2717_21 * r_2 * h8_m2 + e_5 * fs_1260_187_385 * r_4 * h6_m6 - e_5 * fs_7560_187_7 * r_4 * h6_m2 - e_5 * fs_6555_286_6 * r_6 * h4_m2 + e_5 * fs_11475_286_10 * r_8 * h2_m2 + e_6 * fs_2673_96577_36465 * h12_m6 + e_6 * fs_2673_193154_10010 * h12_m2 + e_6 * fs_6300_96577_143 * r_2 * h10_m6 - e_6 * fs_135450_96577_22 * r_2 * h10_m2 - e_6 * fs_345_2717_715 * r_4 * h8_m6 + e_6 * fs_1380_2717_21 * r_4 * h8_m2 - e_6 * fs_1680_3553_385 * r_6 * h6_m6 + e_6 * fs_10080_3553_7 * r_6 * h6_m2 + e_6 * fs_6555_4862_6 * r_8 * h4_m2 - e_6 * fs_306_143_10 * r_10 * h2_m2 - e_7 * fs_198_96577_36465 * r_2 * h12_m6 - e_7 * fs_99_96577_10010 * r_2 * h12_m2 - e_7 * fs_252_96577_143 * r_4 * h10_m6 + e_7 * fs_5418_96577_22 * r_4 * h10_m2 + e_7 * fs_10_2717_715 * r_6 * h8_m6 - e_7 * fs_40_2717_21 * r_6 * h8_m2 + e_7 * fs_40_3553_385 * r_8 * h6_m6 - e_7 * fs_240_3553_7 * r_8 * h6_m2 - e_7 * fs_69_2431_6 * r_10 * h4_m2 + e_7 * fs_6_143_10 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ph12_m7, ph12_m1, ab_2, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m1 = ph12_m1[k];

        pc_32[k] = - e_1 * fs_112455_64_10 * h2_m1 + e_2 * fs_17955_16_3 * h4_m1 + e_2 * fs_48195_16_10 * r_2 * h2_m1 + e_3 * fs_7875_88_70 * h6_m1 - e_3 * fs_89775_88_3 * r_2 * h4_m1 - e_3 * fs_26775_16_10 * r_4 * h2_m1 + e_4 * fs_2415_572_858 * h8_m7 - e_4 * fs_21735_572_30 * h8_m1 - e_4 * fs_525_11_70 * r_2 * h6_m1 + e_4 * fs_89775_286_3 * r_4 * h4_m1 + e_4 * fs_8925_22_10 * r_6 * h2_m1 - e_5 * fs_4725_8398_2431 * h10_m7 - e_5 * fs_64575_16796_66 * h10_m1 - e_5 * fs_7245_5434_858 * r_2 * h8_m7 + e_5 * fs_65205_5434_30 * r_2 * h8_m1 + e_5 * fs_1575_187_70 * r_4 * h6_m1 - e_5 * fs_5985_143_3 * r_6 * h4_m1 - e_5 * fs_26775_572_10 * r_8 * h2_m1 + e_6 * fs_891_96577_461890 * h12_m7 - e_6 * fs_9801_96577_65 * h12_m1 + e_6 * fs_9450_96577_2431 * r_2 * h10_m7 + e_6 * fs_64575_96577_66 * r_2 * h10_m1 + e_6 * fs_345_2717_858 * r_4 * h8_m7 - e_6 * fs_3105_2717_30 * r_4 * h8_m1 - e_6 * fs_2100_3553_70 * r_6 * h6_m1 + e_6 * fs_5985_2431_3 * r_8 * h4_m1 + e_6 * fs_357_143_10 * r_10 * h2_m1 - e_7 * fs_66_96577_461890 * r_2 * h12_m7 + e_7 * fs_726_96577_65 * r_2 * h12_m1 - e_7 * fs_378_96577_2431 * r_4 * h10_m7 - e_7 * fs_2583_96577_66 * r_4 * h10_m1 - e_7 * fs_10_2717_858 * r_6 * h8_m7 + e_7 * fs_90_2717_30 * r_6 * h8_m1 + e_7 * fs_50_3553_70 * r_8 * h6_m1 - e_7 * fs_126_2431_3 * r_10 * h4_m1 - e_7 * fs_7_143_10 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_4, pe_5, pe_6, pe_7, ph8_m8, ph10_m8, ph12_m8, ab_2, pc_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h8_m8 = ph8_m8[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h12_m8 = ph12_m8[k];

        pc_33[k] = e_4 * fs_7245_572_715 * h8_m8 - e_5 * fs_1575_4199_12155 * h10_m8 - e_5 * fs_21735_5434_715 * r_2 * h8_m8 + e_6 * fs_1782_96577_138567 * h12_m8 + e_6 * fs_6300_96577_12155 * r_2 * h10_m8 + e_6 * fs_1035_2717_715 * r_4 * h8_m8 - e_7 * fs_132_96577_138567 * r_2 * h12_m8 - e_7 * fs_252_96577_12155 * r_4 * h10_m8 - e_7 * fs_30_2717_715 * r_6 * h8_m8;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m9, ph10_m1, ph12_m9, ph12_m1, ab_2, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m1 = ph12_m1[k];

        pc_34[k] = e_1 * fs_48195_64_66 * h2_m1 + e_2 * fs_2565_16_55 * h4_m1 - e_2 * fs_20655_16_66 * r_2 * h2_m1 - e_3 * fs_4725_88_462 * h6_m1 - e_3 * fs_12825_88_55 * r_2 * h4_m1 + e_3 * fs_11475_16_66 * r_4 * h2_m1 - e_4 * fs_7245_143_22 * h8_m1 + e_4 * fs_315_11_462 * r_2 * h6_m1 + e_4 * fs_12825_286_55 * r_4 * h4_m1 - e_4 * fs_3825_22_66 * r_6 * h2_m1 - e_5 * fs_1575_8398_20995 * h10_m9 - e_5 * fs_67725_16796_10 * h10_m1 + e_5 * fs_43470_2717_22 * r_2 * h8_m1 - e_5 * fs_945_187_462 * r_4 * h6_m1 - e_5 * fs_855_143_55 * r_6 * h4_m1 + e_5 * fs_11475_572_66 * r_8 * h2_m1 + e_6 * fs_2673_96577_58786 * h12_m9 - e_6 * fs_891_96577_429 * h12_m1 + e_6 * fs_3150_96577_20995 * r_2 * h10_m9 + e_6 * fs_67725_96577_10 * r_2 * h10_m1 - e_6 * fs_4140_2717_22 * r_4 * h8_m1 + e_6 * fs_1260_3553_462 * r_6 * h6_m1 + e_6 * fs_855_2431_55 * r_8 * h4_m1 - e_6 * fs_153_143_66 * r_10 * h2_m1 - e_7 * fs_198_96577_58786 * r_2 * h12_m9 + e_7 * fs_66_96577_429 * r_2 * h12_m1 - e_7 * fs_126_96577_20995 * r_4 * h10_m9 - e_7 * fs_2709_96577_10 * r_4 * h10_m1 + e_7 * fs_120_2717_22 * r_6 * h8_m1 - e_7 * fs_30_3553_462 * r_8 * h6_m1 - e_7 * fs_18_2431_55 * r_10 * h4_m1 + e_7 * fs_3_143_66 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m10, ph10_m2, ph12_m10, ph12_m2, ab_2, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m2 = ph12_m2[k];

        pc_35[k] = - e_1 * fs_16065_32_22 * h2_m2 - e_2 * fs_7695_32_330 * h4_m2 + e_2 * fs_6885_8_22 * r_2 * h2_m2 - e_3 * fs_1575_22_385 * h6_m2 + e_3 * fs_38475_176_330 * r_2 * h4_m2 - e_3 * fs_3825_8_22 * r_4 * h2_m2 - e_4 * fs_2415_572_1155 * h8_m2 + e_4 * fs_420_11_385 * r_2 * h6_m2 - e_4 * fs_38475_572_330 * r_4 * h4_m2 + e_4 * fs_1275_11_22 * r_6 * h2_m2 + e_5 * fs_1575_4199_12597 * h10_m10 - e_5 * fs_14175_8398_10 * h10_m2 + e_5 * fs_7245_5434_1155 * r_2 * h8_m2 - e_5 * fs_1260_187_385 * r_4 * h6_m2 + e_5 * fs_2565_286_330 * r_6 * h4_m2 - e_5 * fs_3825_286_22 * r_8 * h2_m2 + e_6 * fs_891_96577_323323 * h12_m10 - e_6 * fs_891_193154_182 * h12_m2 - e_6 * fs_6300_96577_12597 * r_2 * h10_m10 + e_6 * fs_28350_96577_10 * r_2 * h10_m2 - e_6 * fs_345_2717_1155 * r_4 * h8_m2 + e_6 * fs_1680_3553_385 * r_6 * h6_m2 - e_6 * fs_2565_4862_330 * r_8 * h4_m2 + e_6 * fs_102_143_22 * r_10 * h2_m2 - e_7 * fs_66_96577_323323 * r_2 * h12_m10 + e_7 * fs_33_96577_182 * r_2 * h12_m2 + e_7 * fs_252_96577_12597 * r_4 * h10_m10 - e_7 * fs_1134_96577_10 * r_4 * h10_m2 + e_7 * fs_10_2717_1155 * r_6 * h8_m2 - e_7 * fs_40_3553_385 * r_8 * h6_m2 + e_7 * fs_27_2431_330 * r_10 * h4_m2 - e_7 * fs_2_143_22 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ph12_0, ph12_p6, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p6 = ph12_p6[k];

        pc_36[k] = e_0 * f_155925_64 - e_1 * f_80325_32 * h2_0 - e_1 * f_363825_32 * r_2 - e_2 * f_23085_8 * h4_0 + e_2 * f_34425_8 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 + e_3 * f_67725_44 * h6_0 + e_3 * fs_1575_22_462 * h6_p6 + e_3 * f_115425_44 * r_2 * h4_0 - e_3 * f_19125_8 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 - e_4 * f_45885_572 * h8_0 - e_4 * fs_2415_286_858 * h8_p6 - e_4 * f_9030_11 * r_2 * h6_0 - e_4 * fs_420_11_462 * r_2 * h6_p6 - e_4 * f_115425_143 * r_4 * h4_0 + e_4 * f_6375_11 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 - e_5 * f_23625_442 * h10_0 + e_5 * fs_4725_8398_4290 * h10_p6 + e_5 * f_7245_286 * r_2 * h8_0 + e_5 * fs_7245_2717_858 * r_2 * h8_p6 + e_5 * f_27090_187 * r_4 * h6_0 + e_5 * fs_1260_187_462 * r_4 * h6_p6 + e_5 * f_15390_143 * r_6 * h4_0 - e_5 * f_19125_286 * r_8 * h2_0 - e_5 * f_315_2 * r_10 - e_6 * f_196020_96577 * h12_0 - e_6 * fs_8910_96577_4862 * h12_p6 + e_6 * f_47250_5083 * r_2 * h10_0 - e_6 * fs_9450_96577_4290 * r_2 * h10_p6 - e_6 * f_345_143 * r_4 * h8_0 - e_6 * fs_690_2717_858 * r_4 * h8_p6 - e_6 * f_36120_3553 * r_6 * h6_0 - e_6 * fs_1680_3553_462 * r_6 * h6_p6 - e_6 * f_15390_2431 * r_8 * h4_0 + e_6 * f_510_143 * r_10 * h2_0 + e_6 * f_105_13 * r_12 + e_7 * f_14520_96577 * r_2 * h12_0 + e_7 * fs_660_96577_4862 * r_2 * h12_p6 - e_7 * f_1890_5083 * r_4 * h10_0 + e_7 * fs_378_96577_4290 * r_4 * h10_p6 + e_7 * f_10_143 * r_6 * h8_0 + e_7 * fs_20_2717_858 * r_6 * h8_p6 + e_7 * f_860_3553 * r_8 * h6_0 + e_7 * fs_40_3553_462 * r_8 * h6_p6 + e_7 * f_324_2431 * r_10 * h4_0 - e_7 * f_10_143 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ph12_p1, ph12_p5, ab_2, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p5 = ph12_p5[k];

        pc_37[k] = - e_1 * fs_80325_32_3 * h2_p1 + e_2 * fs_33345_32_10 * h4_p1 + e_2 * fs_34425_8_3 * r_2 * h2_p1 - e_3 * fs_4725_44_21 * h6_p1 + e_3 * fs_1575_44_154 * h6_p5 - e_3 * fs_166725_176_10 * r_2 * h4_p1 - e_3 * fs_19125_8_3 * r_4 * h2_p1 - e_4 * f_2415_22 * h8_p1 - e_4 * fs_2415_572_1001 * h8_p5 + e_4 * fs_630_11_21 * r_2 * h6_p1 - e_4 * fs_210_11_154 * r_2 * h6_p5 + e_4 * fs_12825_44_10 * r_4 * h4_p1 + e_4 * fs_6375_11_3 * r_6 * h2_p1 + e_5 * fs_42525_8398_55 * h10_p1 + e_5 * fs_7875_8398_858 * h10_p5 + e_5 * f_7245_209 * r_2 * h8_p1 + e_5 * fs_7245_5434_1001 * r_2 * h8_p5 - e_5 * fs_1890_187_21 * r_4 * h6_p1 + e_5 * fs_630_187_154 * r_4 * h6_p5 - e_5 * fs_855_22_10 * r_6 * h4_p1 - e_5 * fs_19125_286_3 * r_8 * h2_p1 + e_6 * fs_49005_193154_78 * h12_p1 - e_6 * fs_4455_96577_17017 * h12_p5 - e_6 * fs_85050_96577_55 * r_2 * h10_p1 - e_6 * fs_15750_96577_858 * r_2 * h10_p5 - e_6 * f_690_209 * r_4 * h8_p1 - e_6 * fs_345_2717_1001 * r_4 * h8_p5 + e_6 * fs_2520_3553_21 * r_6 * h6_p1 - e_6 * fs_840_3553_154 * r_6 * h6_p5 + e_6 * fs_855_374_10 * r_8 * h4_p1 + e_6 * fs_510_143_3 * r_10 * h2_p1 - e_7 * fs_1815_96577_78 * r_2 * h12_p1 + e_7 * fs_330_96577_17017 * r_2 * h12_p5 + e_7 * fs_3402_96577_55 * r_4 * h10_p1 + e_7 * fs_630_96577_858 * r_4 * h10_p5 + e_7 * f_20_209 * r_6 * h8_p1 + e_7 * fs_10_2717_1001 * r_6 * h8_p5 - e_7 * fs_60_3553_21 * r_8 * h6_p1 + e_7 * fs_20_3553_154 * r_8 * h6_p5 - e_7 * fs_9_187_10 * r_10 * h4_p1 - e_7 * fs_10_143_3 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ph12_p2, ph12_p4, ab_2, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p4 = ph12_p4[k];

        pc_38[k] = e_1 * fs_16065_16_30 * h2_p2 - e_2 * fs_2565_4_2 * h4_p2 + e_2 * fs_17955_16_14 * h4_p4 - e_2 * fs_6885_4_30 * r_2 * h2_p2 - e_3 * fs_4725_44_21 * h6_p2 - e_3 * fs_7875_88_70 * h6_p4 + e_3 * fs_12825_22_2 * r_2 * h4_p2 - e_3 * fs_89775_88_14 * r_2 * h4_p4 + e_3 * fs_3825_4_30 * r_4 * h2_p2 + e_4 * fs_45885_572_7 * h8_p2 + e_4 * fs_2415_1144_770 * h8_p4 + e_4 * fs_630_11_21 * r_2 * h6_p2 + e_4 * fs_525_11_70 * r_2 * h6_p4 - e_4 * fs_25650_143_2 * r_4 * h4_p2 + e_4 * fs_89775_286_14 * r_4 * h4_p4 - e_4 * fs_2550_11_30 * r_6 * h2_p2 - e_5 * fs_1575_442_66 * h10_p2 + e_5 * fs_4725_16796_858 * h10_p4 - e_5 * fs_7245_286_7 * r_2 * h8_p2 - e_5 * fs_7245_10868_770 * r_2 * h8_p4 - e_5 * fs_1890_187_21 * r_4 * h6_p2 - e_5 * fs_1575_187_70 * r_4 * h6_p4 + e_5 * fs_3420_143_2 * r_6 * h4_p2 - e_5 * fs_5985_143_14 * r_6 * h4_p4 + e_5 * fs_3825_143_30 * r_8 * h2_p2 - e_6 * fs_1782_96577_30030 * h12_p2 - e_6 * fs_7128_96577_5005 * h12_p4 + e_6 * fs_3150_5083_66 * r_2 * h10_p2 - e_6 * fs_4725_96577_858 * r_2 * h10_p4 + e_6 * fs_345_143_7 * r_4 * h8_p2 + e_6 * fs_345_5434_770 * r_4 * h8_p4 + e_6 * fs_2520_3553_21 * r_6 * h6_p2 + e_6 * fs_2100_3553_70 * r_6 * h6_p4 - e_6 * fs_3420_2431_2 * r_8 * h4_p2 + e_6 * fs_5985_2431_14 * r_8 * h4_p4 - e_6 * fs_204_143_30 * r_10 * h2_p2 + e_7 * fs_132_96577_30030 * r_2 * h12_p2 + e_7 * fs_528_96577_5005 * r_2 * h12_p4 - e_7 * fs_126_5083_66 * r_4 * h10_p2 + e_7 * fs_189_96577_858 * r_4 * h10_p4 - e_7 * fs_10_143_7 * r_6 * h8_p2 - e_7 * fs_5_2717_770 * r_6 * h8_p4 - e_7 * fs_60_3553_21 * r_8 * h6_p2 - e_7 * fs_50_3553_70 * r_8 * h6_p4 + e_7 * fs_72_2431_2 * r_10 * h4_p2 - e_7 * fs_126_2431_14 * r_10 * h4_p4 + e_7 * fs_4_143_30 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m3, ph6_m3, ph8_m3, ph10_m3, ph12_m3, ab_2, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m3 = ph12_m3[k];

        pc_39[k] = - e_2 * fs_17955_8_3 * h4_m3 + e_3 * f_67725_44 * h6_m3 + e_3 * fs_89775_44_3 * r_2 * h4_m3 - e_4 * fs_16905_286_22 * h8_m3 - e_4 * f_9030_11 * r_2 * h6_m3 - e_4 * fs_89775_143_3 * r_4 * h4_m3 + e_5 * fs_4725_8398_1001 * h10_m3 + e_5 * fs_50715_2717_22 * r_2 * h8_m3 + e_5 * f_27090_187 * r_4 * h6_m3 + e_5 * fs_11970_143_3 * r_6 * h4_m3 + e_6 * fs_12474_96577_2145 * h12_m3 - e_6 * fs_9450_96577_1001 * r_2 * h10_m3 - e_6 * fs_4830_2717_22 * r_4 * h8_m3 - e_6 * f_36120_3553 * r_6 * h6_m3 - e_6 * fs_11970_2431_3 * r_8 * h4_m3 - e_7 * fs_924_96577_2145 * r_2 * h12_m3 + e_7 * fs_378_96577_1001 * r_4 * h10_m3 + e_7 * fs_140_2717_22 * r_6 * h8_m3 + e_7 * f_860_3553 * r_8 * h6_m3 + e_7 * fs_252_2431_3 * r_10 * h4_m3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ph12_m4, ph12_m2, ab_2, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m4 = ph12_m4[k];
        const auto h12_m2 = ph12_m2[k];

        pc_40[k] = e_1 * fs_16065_16_30 * h2_m2 - e_2 * fs_17955_16_14 * h4_m4 - e_2 * fs_2565_4_2 * h4_m2 - e_2 * fs_6885_4_30 * r_2 * h2_m2 + e_3 * fs_7875_88_70 * h6_m4 - e_3 * fs_4725_44_21 * h6_m2 + e_3 * fs_89775_88_14 * r_2 * h4_m4 + e_3 * fs_12825_22_2 * r_2 * h4_m2 + e_3 * fs_3825_4_30 * r_4 * h2_m2 - e_4 * fs_2415_1144_770 * h8_m4 + e_4 * fs_45885_572_7 * h8_m2 - e_4 * fs_525_11_70 * r_2 * h6_m4 + e_4 * fs_630_11_21 * r_2 * h6_m2 - e_4 * fs_89775_286_14 * r_4 * h4_m4 - e_4 * fs_25650_143_2 * r_4 * h4_m2 - e_4 * fs_2550_11_30 * r_6 * h2_m2 - e_5 * fs_4725_16796_858 * h10_m4 - e_5 * fs_1575_442_66 * h10_m2 + e_5 * fs_7245_10868_770 * r_2 * h8_m4 - e_5 * fs_7245_286_7 * r_2 * h8_m2 + e_5 * fs_1575_187_70 * r_4 * h6_m4 - e_5 * fs_1890_187_21 * r_4 * h6_m2 + e_5 * fs_5985_143_14 * r_6 * h4_m4 + e_5 * fs_3420_143_2 * r_6 * h4_m2 + e_5 * fs_3825_143_30 * r_8 * h2_m2 + e_6 * fs_7128_96577_5005 * h12_m4 - e_6 * fs_1782_96577_30030 * h12_m2 + e_6 * fs_4725_96577_858 * r_2 * h10_m4 + e_6 * fs_3150_5083_66 * r_2 * h10_m2 - e_6 * fs_345_5434_770 * r_4 * h8_m4 + e_6 * fs_345_143_7 * r_4 * h8_m2 - e_6 * fs_2100_3553_70 * r_6 * h6_m4 + e_6 * fs_2520_3553_21 * r_6 * h6_m2 - e_6 * fs_5985_2431_14 * r_8 * h4_m4 - e_6 * fs_3420_2431_2 * r_8 * h4_m2 - e_6 * fs_204_143_30 * r_10 * h2_m2 - e_7 * fs_528_96577_5005 * r_2 * h12_m4 + e_7 * fs_132_96577_30030 * r_2 * h12_m2 - e_7 * fs_189_96577_858 * r_4 * h10_m4 - e_7 * fs_126_5083_66 * r_4 * h10_m2 + e_7 * fs_5_2717_770 * r_6 * h8_m4 - e_7 * fs_10_143_7 * r_6 * h8_m2 + e_7 * fs_50_3553_70 * r_8 * h6_m4 - e_7 * fs_60_3553_21 * r_8 * h6_m2 + e_7 * fs_126_2431_14 * r_10 * h4_m4 + e_7 * fs_72_2431_2 * r_10 * h4_m2 + e_7 * fs_4_143_30 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ph12_m5, ph12_m1, ab_2, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m1 = ph12_m1[k];

        pc_41[k] = - e_1 * fs_80325_32_3 * h2_m1 + e_2 * fs_33345_32_10 * h4_m1 + e_2 * fs_34425_8_3 * r_2 * h2_m1 - e_3 * fs_1575_44_154 * h6_m5 - e_3 * fs_4725_44_21 * h6_m1 - e_3 * fs_166725_176_10 * r_2 * h4_m1 - e_3 * fs_19125_8_3 * r_4 * h2_m1 + e_4 * fs_2415_572_1001 * h8_m5 - e_4 * f_2415_22 * h8_m1 + e_4 * fs_210_11_154 * r_2 * h6_m5 + e_4 * fs_630_11_21 * r_2 * h6_m1 + e_4 * fs_12825_44_10 * r_4 * h4_m1 + e_4 * fs_6375_11_3 * r_6 * h2_m1 - e_5 * fs_7875_8398_858 * h10_m5 + e_5 * fs_42525_8398_55 * h10_m1 - e_5 * fs_7245_5434_1001 * r_2 * h8_m5 + e_5 * f_7245_209 * r_2 * h8_m1 - e_5 * fs_630_187_154 * r_4 * h6_m5 - e_5 * fs_1890_187_21 * r_4 * h6_m1 - e_5 * fs_855_22_10 * r_6 * h4_m1 - e_5 * fs_19125_286_3 * r_8 * h2_m1 + e_6 * fs_4455_96577_17017 * h12_m5 + e_6 * fs_49005_193154_78 * h12_m1 + e_6 * fs_15750_96577_858 * r_2 * h10_m5 - e_6 * fs_85050_96577_55 * r_2 * h10_m1 + e_6 * fs_345_2717_1001 * r_4 * h8_m5 - e_6 * f_690_209 * r_4 * h8_m1 + e_6 * fs_840_3553_154 * r_6 * h6_m5 + e_6 * fs_2520_3553_21 * r_6 * h6_m1 + e_6 * fs_855_374_10 * r_8 * h4_m1 + e_6 * fs_510_143_3 * r_10 * h2_m1 - e_7 * fs_330_96577_17017 * r_2 * h12_m5 - e_7 * fs_1815_96577_78 * r_2 * h12_m1 - e_7 * fs_630_96577_858 * r_4 * h10_m5 + e_7 * fs_3402_96577_55 * r_4 * h10_m1 - e_7 * fs_10_2717_1001 * r_6 * h8_m5 + e_7 * f_20_209 * r_6 * h8_m1 - e_7 * fs_20_3553_154 * r_8 * h6_m5 - e_7 * fs_60_3553_21 * r_8 * h6_m1 - e_7 * fs_9_187_10 * r_10 * h4_m1 - e_7 * fs_10_143_3 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, pe_7, ph6_m6, ph8_m6, ph10_m6, ph12_m6, ab_2, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h6_m6 = ph6_m6[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h12_m6 = ph12_m6[k];

        pc_42[k] = - e_3 * fs_1575_22_462 * h6_m6 + e_4 * fs_2415_286_858 * h8_m6 + e_4 * fs_420_11_462 * r_2 * h6_m6 - e_5 * fs_4725_8398_4290 * h10_m6 - e_5 * fs_7245_2717_858 * r_2 * h8_m6 - e_5 * fs_1260_187_462 * r_4 * h6_m6 + e_6 * fs_8910_96577_4862 * h12_m6 + e_6 * fs_9450_96577_4290 * r_2 * h10_m6 + e_6 * fs_690_2717_858 * r_4 * h8_m6 + e_6 * fs_1680_3553_462 * r_6 * h6_m6 - e_7 * fs_660_96577_4862 * r_2 * h12_m6 - e_7 * fs_378_96577_4290 * r_4 * h10_m6 - e_7 * fs_20_2717_858 * r_6 * h8_m6 - e_7 * fs_40_3553_462 * r_8 * h6_m6;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ph12_m7, ph12_m1, ab_2, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m1 = ph12_m1[k];

        pc_43[k] = e_1 * fs_112455_64_10 * h2_m1 - e_2 * fs_17955_16_3 * h4_m1 - e_2 * fs_48195_16_10 * r_2 * h2_m1 - e_3 * fs_7875_88_70 * h6_m1 + e_3 * fs_89775_88_3 * r_2 * h4_m1 + e_3 * fs_26775_16_10 * r_4 * h2_m1 + e_4 * fs_2415_572_858 * h8_m7 + e_4 * fs_21735_572_30 * h8_m1 + e_4 * fs_525_11_70 * r_2 * h6_m1 - e_4 * fs_89775_286_3 * r_4 * h4_m1 - e_4 * fs_8925_22_10 * r_6 * h2_m1 - e_5 * fs_4725_8398_2431 * h10_m7 + e_5 * fs_64575_16796_66 * h10_m1 - e_5 * fs_7245_5434_858 * r_2 * h8_m7 - e_5 * fs_65205_5434_30 * r_2 * h8_m1 - e_5 * fs_1575_187_70 * r_4 * h6_m1 + e_5 * fs_5985_143_3 * r_6 * h4_m1 + e_5 * fs_26775_572_10 * r_8 * h2_m1 + e_6 * fs_891_96577_461890 * h12_m7 + e_6 * fs_9801_96577_65 * h12_m1 + e_6 * fs_9450_96577_2431 * r_2 * h10_m7 - e_6 * fs_64575_96577_66 * r_2 * h10_m1 + e_6 * fs_345_2717_858 * r_4 * h8_m7 + e_6 * fs_3105_2717_30 * r_4 * h8_m1 + e_6 * fs_2100_3553_70 * r_6 * h6_m1 - e_6 * fs_5985_2431_3 * r_8 * h4_m1 - e_6 * fs_357_143_10 * r_10 * h2_m1 - e_7 * fs_66_96577_461890 * r_2 * h12_m7 - e_7 * fs_726_96577_65 * r_2 * h12_m1 - e_7 * fs_378_96577_2431 * r_4 * h10_m7 + e_7 * fs_2583_96577_66 * r_4 * h10_m1 - e_7 * fs_10_2717_858 * r_6 * h8_m7 - e_7 * fs_90_2717_30 * r_6 * h8_m1 - e_7 * fs_50_3553_70 * r_8 * h6_m1 + e_7 * fs_126_2431_3 * r_10 * h4_m1 + e_7 * fs_7_143_10 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ph12_m8, ph12_m2, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m2 = ph12_m2[k];

        pc_44[k] = - e_1 * fs_16065_32_55 * h2_m2 - e_2 * fs_2565_4_33 * h4_m2 + e_2 * fs_6885_8_55 * r_2 * h2_m2 + e_3 * fs_1575_44_154 * h6_m2 + e_3 * fs_12825_22_33 * r_2 * h4_m2 - e_3 * fs_3825_8_55 * r_4 * h2_m2 - e_4 * fs_2415_52_39 * h8_m8 + e_4 * fs_7245_572_462 * h8_m2 - e_4 * fs_210_11_154 * r_2 * h6_m2 - e_4 * fs_25650_143_33 * r_4 * h4_m2 + e_4 * fs_1275_11_55 * r_6 * h2_m2 + e_5 * fs_1575_8398_663 * h10_m8 + e_5 * f_174825_8398 * h10_m2 + e_5 * fs_7245_494_39 * r_2 * h8_m8 - e_5 * fs_21735_5434_462 * r_2 * h8_m2 + e_5 * fs_630_187_154 * r_4 * h6_m2 + e_5 * fs_3420_143_33 * r_6 * h4_m2 - e_5 * fs_3825_286_55 * r_8 * h2_m2 + e_6 * fs_3564_96577_20995 * h12_m8 + e_6 * fs_1782_96577_455 * h12_m2 - e_6 * fs_3150_96577_663 * r_2 * h10_m8 - e_6 * f_349650_96577 * r_2 * h10_m2 - e_6 * fs_345_247_39 * r_4 * h8_m8 + e_6 * fs_1035_2717_462 * r_4 * h8_m2 - e_6 * fs_840_3553_154 * r_6 * h6_m2 - e_6 * fs_3420_2431_33 * r_8 * h4_m2 + e_6 * fs_102_143_55 * r_10 * h2_m2 - e_7 * fs_264_96577_20995 * r_2 * h12_m8 - e_7 * fs_132_96577_455 * r_2 * h12_m2 + e_7 * fs_126_96577_663 * r_4 * h10_m8 + e_7 * f_13986_96577 * r_4 * h10_m2 + e_7 * fs_10_247_39 * r_6 * h8_m8 - e_7 * fs_30_2717_462 * r_6 * h8_m2 + e_7 * fs_20_3553_154 * r_8 * h6_m2 + e_7 * fs_72_2431_33 * r_10 * h4_m2 - e_7 * fs_2_143_55 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m3, ph6_m3, ph8_m3, ph10_m9, ph10_m3, ph12_m9, ph12_m3, ab_2, pc_45 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m3 = ph12_m3[k];

        pc_45[k] = e_2 * fs_7695_32_154 * h4_m3 + e_3 * fs_1575_22_462 * h6_m3 - e_3 * fs_38475_176_154 * r_2 * h4_m3 + e_4 * fs_2415_52_21 * h8_m3 - e_4 * fs_420_11_462 * r_2 * h6_m3 + e_4 * fs_38475_572_154 * r_4 * h4_m3 + e_5 * fs_4725_8398_8398 * h10_m9 + e_5 * fs_4725_4199_78 * h10_m3 - e_5 * fs_7245_494_21 * r_2 * h8_m3 + e_5 * fs_1260_187_462 * r_4 * h6_m3 - e_5 * fs_2565_286_154 * r_6 * h4_m3 + e_6 * fs_891_96577_146965 * h12_m9 + e_6 * fs_891_193154_910 * h12_m3 - e_6 * fs_9450_96577_8398 * r_2 * h10_m9 - e_6 * fs_18900_96577_78 * r_2 * h10_m3 + e_6 * fs_345_247_21 * r_4 * h8_m3 - e_6 * fs_1680_3553_462 * r_6 * h6_m3 + e_6 * fs_2565_4862_154 * r_8 * h4_m3 - e_7 * fs_66_96577_146965 * r_2 * h12_m9 - e_7 * fs_33_96577_910 * r_2 * h12_m3 + e_7 * fs_378_96577_8398 * r_4 * h10_m9 + e_7 * fs_756_96577_78 * r_4 * h10_m3 - e_7 * fs_10_247_21 * r_6 * h8_m3 + e_7 * fs_40_3553_462 * r_8 * h6_m3 - e_7 * fs_27_2431_154 * r_10 * h4_m3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ph12_0, ph12_p4, ab_2, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p4 = ph12_p4[k];

        pc_46[k] = e_0 * f_155925_64 - e_1 * f_80325_16 * h2_0 - e_1 * f_363825_32 * r_2 + e_2 * f_9405_16 * h4_0 - e_2 * fs_5985_8_35 * h4_p4 + e_2 * f_34425_4 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 + e_3 * f_1575_2 * h6_0 + e_3 * fs_4725_11_7 * h6_p4 - e_3 * f_4275_8 * r_2 * h4_0 + e_3 * fs_29925_44_35 * r_2 * h4_p4 - e_3 * f_19125_4 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 - e_4 * f_171465_572 * h8_0 - e_4 * fs_7245_286_77 * h8_p4 - e_4 * f_420_1 * r_2 * h6_0 - e_4 * fs_2520_11_7 * r_2 * h6_p4 + e_4 * f_4275_26 * r_4 * h4_0 - e_4 * fs_29925_143_35 * r_4 * h4_p4 + e_4 * f_12750_11 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 + e_5 * f_165375_4199 * h10_0 + e_5 * fs_3150_4199_2145 * h10_p4 + e_5 * f_514395_5434 * r_2 * h8_0 + e_5 * fs_21735_2717_77 * r_2 * h8_p4 + e_5 * f_1260_17 * r_4 * h6_0 + e_5 * fs_7560_187_7 * r_4 * h6_p4 - e_5 * f_285_13 * r_6 * h4_0 + e_5 * fs_3990_143_35 * r_6 * h4_p4 - e_5 * f_19125_143 * r_8 * h2_0 - e_5 * f_315_2 * r_10 + e_6 * f_441045_96577 * h12_0 - e_6 * fs_13365_96577_2002 * h12_p4 - e_6 * f_661500_96577 * r_2 * h10_0 - e_6 * fs_12600_96577_2145 * r_2 * h10_p4 - e_6 * f_24495_2717 * r_4 * h8_0 - e_6 * fs_2070_2717_77 * r_4 * h8_p4 - e_6 * f_1680_323 * r_6 * h6_0 - e_6 * fs_10080_3553_7 * r_6 * h6_p4 + e_6 * f_285_221 * r_8 * h4_0 - e_6 * fs_3990_2431_35 * r_8 * h4_p4 + e_6 * f_1020_143 * r_10 * h2_0 + e_6 * f_105_13 * r_12 - e_7 * f_32670_96577 * r_2 * h12_0 + e_7 * fs_990_96577_2002 * r_2 * h12_p4 + e_7 * f_26460_96577 * r_4 * h10_0 + e_7 * fs_504_96577_2145 * r_4 * h10_p4 + e_7 * f_710_2717 * r_6 * h8_0 + e_7 * fs_60_2717_77 * r_6 * h8_p4 + e_7 * f_40_323 * r_8 * h6_0 + e_7 * fs_240_3553_7 * r_8 * h6_p4 - e_7 * f_6_221 * r_10 * h4_0 + e_7 * fs_84_2431_35 * r_10 * h4_p4 - e_7 * f_20_143 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ph12_p1, ph12_p3, ab_2, pc_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p3 = ph12_p3[k];

        pc_47[k] = - e_1 * fs_16065_32_30 * h2_p1 + e_2 * f_45315_16 * h4_p1 - e_2 * fs_5985_16_7 * h4_p3 + e_2 * fs_6885_8_30 * r_2 * h2_p1 - e_3 * fs_1575_22_210 * h6_p1 + e_3 * fs_4725_44_21 * h6_p3 - e_3 * f_226575_88 * r_2 * h4_p1 + e_3 * fs_29925_88_7 * r_2 * h4_p3 - e_3 * fs_3825_8_30 * r_4 * h2_p1 + e_4 * fs_65205_1144_10 * h8_p1 - e_4 * fs_7245_1144_462 * h8_p3 + e_4 * fs_420_11_210 * r_2 * h6_p1 - e_4 * fs_630_11_21 * r_2 * h6_p3 + e_4 * f_226575_286 * r_4 * h4_p1 - e_4 * fs_29925_286_7 * r_4 * h4_p3 + e_4 * fs_1275_11_30 * r_6 * h2_p1 - e_5 * fs_9450_4199_22 * h10_p1 + e_5 * fs_11025_8398_429 * h10_p3 - e_5 * fs_195615_10868_10 * r_2 * h8_p1 + e_5 * fs_21735_10868_462 * r_2 * h8_p3 - e_5 * fs_1260_187_210 * r_4 * h6_p1 + e_5 * fs_1890_187_21 * r_4 * h6_p3 - e_5 * f_15105_143 * r_6 * h4_p1 + e_5 * fs_1995_143_7 * r_6 * h4_p3 - e_5 * fs_3825_286_30 * r_8 * h2_p1 - e_6 * fs_29403_96577_195 * h12_p1 - e_6 * fs_8019_96577_5005 * h12_p3 + e_6 * fs_37800_96577_22 * r_2 * h10_p1 - e_6 * fs_22050_96577_429 * r_2 * h10_p3 + e_6 * fs_9315_5434_10 * r_4 * h8_p1 - e_6 * fs_1035_5434_462 * r_4 * h8_p3 + e_6 * fs_1680_3553_210 * r_6 * h6_p1 - e_6 * fs_2520_3553_21 * r_6 * h6_p3 + e_6 * f_15105_2431 * r_8 * h4_p1 - e_6 * fs_1995_2431_7 * r_8 * h4_p3 + e_6 * fs_102_143_30 * r_10 * h2_p1 + e_7 * fs_2178_96577_195 * r_2 * h12_p1 + e_7 * fs_594_96577_5005 * r_2 * h12_p3 - e_7 * fs_1512_96577_22 * r_4 * h10_p1 + e_7 * fs_882_96577_429 * r_4 * h10_p3 - e_7 * fs_135_2717_10 * r_6 * h8_p1 + e_7 * fs_15_2717_462 * r_6 * h8_p3 - e_7 * fs_40_3553_210 * r_8 * h6_p1 + e_7 * fs_60_3553_21 * r_8 * h6_p3 - e_7 * f_318_2431 * r_10 * h4_p1 + e_7 * fs_42_2431_7 * r_10 * h4_p3 - e_7 * fs_2_143_30 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ph12_m2, ab_2, pc_48 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m2 = ph12_m2[k];

        pc_48[k] = e_1 * fs_16065_16_70 * h2_m2 - e_2 * fs_9405_16_42 * h4_m2 - e_2 * fs_6885_4_70 * r_2 * h2_m2 + e_3 * f_1575_2 * h6_m2 + e_3 * fs_4275_8_42 * r_2 * h4_m2 + e_3 * fs_3825_4_70 * r_4 * h2_m2 - e_4 * fs_16905_572_3 * h8_m2 - e_4 * f_420_1 * r_2 * h6_m2 - e_4 * fs_4275_26_42 * r_4 * h4_m2 - e_4 * fs_2550_11_70 * r_6 * h2_m2 - e_5 * fs_4725_4199_154 * h10_m2 + e_5 * fs_50715_5434_3 * r_2 * h8_m2 + e_5 * f_1260_17 * r_4 * h6_m2 + e_5 * fs_285_13_42 * r_6 * h4_m2 + e_5 * fs_3825_143_70 * r_8 * h2_m2 + e_6 * fs_18711_96577_1430 * h12_m2 + e_6 * fs_18900_96577_154 * r_2 * h10_m2 - e_6 * fs_2415_2717_3 * r_4 * h8_m2 - e_6 * f_1680_323 * r_6 * h6_m2 - e_6 * fs_285_221_42 * r_8 * h4_m2 - e_6 * fs_204_143_70 * r_10 * h2_m2 - e_7 * fs_1386_96577_1430 * r_2 * h12_m2 - e_7 * fs_756_96577_154 * r_4 * h10_m2 + e_7 * fs_70_2717_3 * r_6 * h8_m2 + e_7 * f_40_323 * r_8 * h6_m2 + e_7 * fs_6_221_42 * r_10 * h4_m2 + e_7 * fs_4_143_70 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ph12_m3, ph12_m1, ab_2, pc_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m3 = ph12_m3[k];
        const auto h12_m1 = ph12_m1[k];

        pc_49[k] = - e_1 * fs_16065_32_30 * h2_m1 + e_2 * fs_5985_16_7 * h4_m3 + e_2 * f_45315_16 * h4_m1 + e_2 * fs_6885_8_30 * r_2 * h2_m1 - e_3 * fs_4725_44_21 * h6_m3 - e_3 * fs_1575_22_210 * h6_m1 - e_3 * fs_29925_88_7 * r_2 * h4_m3 - e_3 * f_226575_88 * r_2 * h4_m1 - e_3 * fs_3825_8_30 * r_4 * h2_m1 + e_4 * fs_7245_1144_462 * h8_m3 + e_4 * fs_65205_1144_10 * h8_m1 + e_4 * fs_630_11_21 * r_2 * h6_m3 + e_4 * fs_420_11_210 * r_2 * h6_m1 + e_4 * fs_29925_286_7 * r_4 * h4_m3 + e_4 * f_226575_286 * r_4 * h4_m1 + e_4 * fs_1275_11_30 * r_6 * h2_m1 - e_5 * fs_11025_8398_429 * h10_m3 - e_5 * fs_9450_4199_22 * h10_m1 - e_5 * fs_21735_10868_462 * r_2 * h8_m3 - e_5 * fs_195615_10868_10 * r_2 * h8_m1 - e_5 * fs_1890_187_21 * r_4 * h6_m3 - e_5 * fs_1260_187_210 * r_4 * h6_m1 - e_5 * fs_1995_143_7 * r_6 * h4_m3 - e_5 * f_15105_143 * r_6 * h4_m1 - e_5 * fs_3825_286_30 * r_8 * h2_m1 + e_6 * fs_8019_96577_5005 * h12_m3 - e_6 * fs_29403_96577_195 * h12_m1 + e_6 * fs_22050_96577_429 * r_2 * h10_m3 + e_6 * fs_37800_96577_22 * r_2 * h10_m1 + e_6 * fs_1035_5434_462 * r_4 * h8_m3 + e_6 * fs_9315_5434_10 * r_4 * h8_m1 + e_6 * fs_2520_3553_21 * r_6 * h6_m3 + e_6 * fs_1680_3553_210 * r_6 * h6_m1 + e_6 * fs_1995_2431_7 * r_8 * h4_m3 + e_6 * f_15105_2431 * r_8 * h4_m1 + e_6 * fs_102_143_30 * r_10 * h2_m1 - e_7 * fs_594_96577_5005 * r_2 * h12_m3 + e_7 * fs_2178_96577_195 * r_2 * h12_m1 - e_7 * fs_882_96577_429 * r_4 * h10_m3 - e_7 * fs_1512_96577_22 * r_4 * h10_m1 - e_7 * fs_15_2717_462 * r_6 * h8_m3 - e_7 * fs_135_2717_10 * r_6 * h8_m1 - e_7 * fs_60_3553_21 * r_8 * h6_m3 - e_7 * fs_40_3553_210 * r_8 * h6_m1 - e_7 * fs_42_2431_7 * r_10 * h4_m3 - e_7 * f_318_2431 * r_10 * h4_m1 - e_7 * fs_2_143_30 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ph12_m4, ab_2, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m4 = ph12_m4[k];

        pc_50[k] = e_2 * fs_5985_8_35 * h4_m4 - e_3 * fs_4725_11_7 * h6_m4 - e_3 * fs_29925_44_35 * r_2 * h4_m4 + e_4 * fs_7245_286_77 * h8_m4 + e_4 * fs_2520_11_7 * r_2 * h6_m4 + e_4 * fs_29925_143_35 * r_4 * h4_m4 - e_5 * fs_3150_4199_2145 * h10_m4 - e_5 * fs_21735_2717_77 * r_2 * h8_m4 - e_5 * fs_7560_187_7 * r_4 * h6_m4 - e_5 * fs_3990_143_35 * r_6 * h4_m4 + e_6 * fs_13365_96577_2002 * h12_m4 + e_6 * fs_12600_96577_2145 * r_2 * h10_m4 + e_6 * fs_2070_2717_77 * r_4 * h8_m4 + e_6 * fs_10080_3553_7 * r_6 * h6_m4 + e_6 * fs_3990_2431_35 * r_8 * h4_m4 - e_7 * fs_990_96577_2002 * r_2 * h12_m4 - e_7 * fs_504_96577_2145 * r_4 * h10_m4 - e_7 * fs_60_2717_77 * r_6 * h8_m4 - e_7 * fs_240_3553_7 * r_8 * h6_m4 - e_7 * fs_84_2431_35 * r_10 * h4_m4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ph12_m5, ph12_m1, ab_2, pc_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m1 = ph12_m1[k];

        pc_51[k] = e_1 * fs_80325_32_3 * h2_m1 - e_2 * fs_33345_32_10 * h4_m1 - e_2 * fs_34425_8_3 * r_2 * h2_m1 - e_3 * fs_1575_44_154 * h6_m5 + e_3 * fs_4725_44_21 * h6_m1 + e_3 * fs_166725_176_10 * r_2 * h4_m1 + e_3 * fs_19125_8_3 * r_4 * h2_m1 + e_4 * fs_2415_572_1001 * h8_m5 + e_4 * f_2415_22 * h8_m1 + e_4 * fs_210_11_154 * r_2 * h6_m5 - e_4 * fs_630_11_21 * r_2 * h6_m1 - e_4 * fs_12825_44_10 * r_4 * h4_m1 - e_4 * fs_6375_11_3 * r_6 * h2_m1 - e_5 * fs_7875_8398_858 * h10_m5 - e_5 * fs_42525_8398_55 * h10_m1 - e_5 * fs_7245_5434_1001 * r_2 * h8_m5 - e_5 * f_7245_209 * r_2 * h8_m1 - e_5 * fs_630_187_154 * r_4 * h6_m5 + e_5 * fs_1890_187_21 * r_4 * h6_m1 + e_5 * fs_855_22_10 * r_6 * h4_m1 + e_5 * fs_19125_286_3 * r_8 * h2_m1 + e_6 * fs_4455_96577_17017 * h12_m5 - e_6 * fs_49005_193154_78 * h12_m1 + e_6 * fs_15750_96577_858 * r_2 * h10_m5 + e_6 * fs_85050_96577_55 * r_2 * h10_m1 + e_6 * fs_345_2717_1001 * r_4 * h8_m5 + e_6 * f_690_209 * r_4 * h8_m1 + e_6 * fs_840_3553_154 * r_6 * h6_m5 - e_6 * fs_2520_3553_21 * r_6 * h6_m1 - e_6 * fs_855_374_10 * r_8 * h4_m1 - e_6 * fs_510_143_3 * r_10 * h2_m1 - e_7 * fs_330_96577_17017 * r_2 * h12_m5 + e_7 * fs_1815_96577_78 * r_2 * h12_m1 - e_7 * fs_630_96577_858 * r_4 * h10_m5 - e_7 * fs_3402_96577_55 * r_4 * h10_m1 - e_7 * fs_10_2717_1001 * r_6 * h8_m5 - e_7 * f_20_209 * r_6 * h8_m1 - e_7 * fs_20_3553_154 * r_8 * h6_m5 + e_7 * fs_60_3553_21 * r_8 * h6_m1 + e_7 * fs_9_187_10 * r_10 * h4_m1 + e_7 * fs_10_143_3 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ph12_m6, ph12_m2, ab_2, pc_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m2 = ph12_m2[k];

        pc_52[k] = - e_1 * fs_48195_32_10 * h2_m2 - e_2 * fs_19665_32_6 * h4_m2 + e_2 * fs_20655_8_10 * r_2 * h2_m2 + e_3 * fs_1575_22_385 * h6_m6 + e_3 * fs_4725_11_7 * h6_m2 + e_3 * fs_98325_176_6 * r_2 * h4_m2 - e_3 * fs_11475_8_10 * r_4 * h2_m2 - e_4 * fs_2415_572_715 * h8_m6 - e_4 * fs_2415_143_21 * h8_m2 - e_4 * fs_420_11_385 * r_2 * h6_m6 - e_4 * fs_2520_11_7 * r_2 * h6_m2 - e_4 * fs_98325_572_6 * r_4 * h4_m2 + e_4 * fs_3825_11_10 * r_6 * h2_m2 - e_5 * fs_1575_4199_143 * h10_m6 - e_5 * fs_67725_8398_22 * h10_m2 + e_5 * fs_7245_5434_715 * r_2 * h8_m6 + e_5 * fs_14490_2717_21 * r_2 * h8_m2 + e_5 * fs_1260_187_385 * r_4 * h6_m6 + e_5 * fs_7560_187_7 * r_4 * h6_m2 + e_5 * fs_6555_286_6 * r_6 * h4_m2 - e_5 * fs_11475_286_10 * r_8 * h2_m2 + e_6 * fs_2673_96577_36465 * h12_m6 - e_6 * fs_2673_193154_10010 * h12_m2 + e_6 * fs_6300_96577_143 * r_2 * h10_m6 + e_6 * fs_135450_96577_22 * r_2 * h10_m2 - e_6 * fs_345_2717_715 * r_4 * h8_m6 - e_6 * fs_1380_2717_21 * r_4 * h8_m2 - e_6 * fs_1680_3553_385 * r_6 * h6_m6 - e_6 * fs_10080_3553_7 * r_6 * h6_m2 - e_6 * fs_6555_4862_6 * r_8 * h4_m2 + e_6 * fs_306_143_10 * r_10 * h2_m2 - e_7 * fs_198_96577_36465 * r_2 * h12_m6 + e_7 * fs_99_96577_10010 * r_2 * h12_m2 - e_7 * fs_252_96577_143 * r_4 * h10_m6 - e_7 * fs_5418_96577_22 * r_4 * h10_m2 + e_7 * fs_10_2717_715 * r_6 * h8_m6 + e_7 * fs_40_2717_21 * r_6 * h8_m2 + e_7 * fs_40_3553_385 * r_8 * h6_m6 + e_7 * fs_240_3553_7 * r_8 * h6_m2 + e_7 * fs_69_2431_6 * r_10 * h4_m2 - e_7 * fs_6_143_10 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ph12_m7, ph12_m3, ab_2, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m3 = ph12_m3[k];

        pc_53[k] = e_2 * fs_5985_32_462 * h4_m3 + e_3 * fs_1575_44_154 * h6_m3 - e_3 * fs_29925_176_462 * r_2 * h4_m3 - e_4 * fs_2415_52_39 * h8_m7 - e_4 * fs_2415_26_7 * h8_m3 - e_4 * fs_210_11_154 * r_2 * h6_m3 + e_4 * fs_29925_572_462 * r_4 * h4_m3 + e_5 * fs_11025_8398_442 * h10_m7 - e_5 * fs_48825_8398_26 * h10_m3 + e_5 * fs_7245_494_39 * r_2 * h8_m7 + e_5 * fs_7245_247_7 * r_2 * h8_m3 + e_5 * fs_630_187_154 * r_4 * h6_m3 - e_5 * fs_1995_286_462 * r_6 * h4_m3 + e_6 * fs_2673_96577_20995 * h12_m7 - e_6 * fs_2673_193154_2730 * h12_m3 - e_6 * fs_22050_96577_442 * r_2 * h10_m7 + e_6 * fs_97650_96577_26 * r_2 * h10_m3 - e_6 * fs_345_247_39 * r_4 * h8_m7 - e_6 * fs_690_247_7 * r_4 * h8_m3 - e_6 * fs_840_3553_154 * r_6 * h6_m3 + e_6 * fs_1995_4862_462 * r_8 * h4_m3 - e_7 * fs_198_96577_20995 * r_2 * h12_m7 + e_7 * fs_99_96577_2730 * r_2 * h12_m3 + e_7 * fs_882_96577_442 * r_4 * h10_m7 - e_7 * fs_3906_96577_26 * r_4 * h10_m3 + e_7 * fs_10_247_39 * r_6 * h8_m7 + e_7 * fs_20_247_7 * r_6 * h8_m3 + e_7 * fs_20_3553_154 * r_8 * h6_m3 - e_7 * fs_21_2431_462 * r_10 * h4_m3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m4, ph6_m4, ph8_m8, ph8_m4, ph10_m8, ph10_m4, ph12_m8, ph12_m4, ab_2, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m4 = ph12_m4[k];

        pc_54[k] = - e_2 * fs_2565_16_77 * h4_m4 - e_3 * fs_1575_22_385 * h6_m4 + e_3 * fs_12825_88_77 * r_2 * h4_m4 + e_4 * fs_2415_52_13 * h8_m8 - e_4 * fs_2415_52_35 * h8_m4 + e_4 * fs_420_11_385 * r_2 * h6_m4 - e_4 * fs_12825_286_77 * r_4 * h4_m4 + e_5 * fs_14175_4199_221 * h10_m8 - e_5 * fs_11025_4199_39 * h10_m4 - e_5 * fs_7245_494_13 * r_2 * h8_m8 + e_5 * fs_7245_494_35 * r_2 * h8_m4 - e_5 * fs_1260_187_385 * r_4 * h6_m4 + e_5 * fs_855_143_77 * r_6 * h4_m4 + e_6 * fs_891_96577_62985 * h12_m8 - e_6 * fs_891_96577_910 * h12_m4 - e_6 * fs_56700_96577_221 * r_2 * h10_m8 + e_6 * fs_44100_96577_39 * r_2 * h10_m4 + e_6 * fs_345_247_13 * r_4 * h8_m8 - e_6 * fs_345_247_35 * r_4 * h8_m4 + e_6 * fs_1680_3553_385 * r_6 * h6_m4 - e_6 * fs_855_2431_77 * r_8 * h4_m4 - e_7 * fs_66_96577_62985 * r_2 * h12_m8 + e_7 * fs_66_96577_910 * r_2 * h12_m4 + e_7 * fs_2268_96577_221 * r_4 * h10_m8 - e_7 * fs_1764_96577_39 * r_4 * h10_m4 - e_7 * fs_10_247_13 * r_6 * h8_m8 + e_7 * fs_10_247_35 * r_6 * h8_m4 - e_7 * fs_40_3553_385 * r_8 * h6_m4 + e_7 * fs_18_2431_77 * r_10 * h4_m4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ph12_0, ph12_p2, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

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
        const auto h12_0 = ph12_0[k];
        const auto h12_p2 = ph12_p2[k];

        pc_55[k] = e_0 * f_155925_64 - e_1 * f_208845_32 * h2_0 + e_1 * fs_112455_32_3 * h2_p2 - e_1 * f_363825_32 * r_2 + e_2 * f_3420_1 * h4_0 - e_2 * fs_5985_4_5 * h4_p2 + e_2 * f_89505_8 * r_2 * h2_0 - e_2 * fs_48195_8_3 * r_2 * h2_p2 + e_2 * f_218295_16 * r_4 - e_3 * f_7875_11 * h6_0 + e_3 * fs_1575_22_210 * h6_p2 - e_3 * f_34200_11 * r_2 * h4_0 + e_3 * fs_29925_22_5 * r_2 * h4_p2 - e_3 * f_49725_8 * r_4 * h2_0 + e_3 * fs_26775_8_3 * r_4 * h2_p2 - e_3 * f_51975_8 * r_6 + e_4 * f_12075_286 * h8_0 - e_4 * fs_7245_286_70 * h8_p2 + e_4 * f_4200_11 * r_2 * h6_0 - e_4 * fs_420_11_210 * r_2 * h6_p2 + e_4 * f_136800_143 * r_4 * h4_0 - e_4 * fs_59850_143_5 * r_4 * h4_p2 + e_4 * f_16575_11 * r_6 * h2_0 - e_4 * fs_8925_11_3 * r_6 * h2_p2 + e_4 * f_5775_4 * r_8 + e_5 * f_4725_323 * h10_0 + e_5 * fs_11025_4199_165 * h10_p2 - e_5 * f_36225_2717 * r_2 * h8_0 + e_5 * fs_21735_2717_70 * r_2 * h8_p2 - e_5 * f_12600_187 * r_4 * h6_0 + e_5 * fs_1260_187_210 * r_4 * h6_p2 - e_5 * f_18240_143 * r_6 * h4_0 + e_5 * fs_7980_143_5 * r_6 * h4_p2 - e_5 * f_3825_22 * r_8 * h2_0 + e_5 * fs_26775_286_3 * r_8 * h2_p2 - e_5 * f_315_2 * r_10 - e_6 * f_705672_96577 * h12_0 - e_6 * fs_10692_96577_3003 * h12_p2 - e_6 * f_18900_7429 * r_2 * h10_0 - e_6 * fs_44100_96577_165 * r_2 * h10_p2 + e_6 * f_3450_2717 * r_4 * h8_0 - e_6 * fs_2070_2717_70 * r_4 * h8_p2 + e_6 * f_16800_3553 * r_6 * h6_0 - e_6 * fs_1680_3553_210 * r_6 * h6_p2 + e_6 * f_18240_2431 * r_8 * h4_0 - e_6 * fs_7980_2431_5 * r_8 * h4_p2 + e_6 * f_102_11 * r_10 * h2_0 - e_6 * fs_714_143_3 * r_10 * h2_p2 + e_6 * f_105_13 * r_12 + e_7 * f_52272_96577 * r_2 * h12_0 + e_7 * fs_792_96577_3003 * r_2 * h12_p2 + e_7 * f_756_7429 * r_4 * h10_0 + e_7 * fs_1764_96577_165 * r_4 * h10_p2 - e_7 * f_100_2717 * r_6 * h8_0 + e_7 * fs_60_2717_70 * r_6 * h8_p2 - e_7 * f_400_3553 * r_8 * h6_0 + e_7 * fs_40_3553_210 * r_8 * h6_p2 - e_7 * f_384_2431 * r_10 * h4_0 + e_7 * fs_168_2431_5 * r_10 * h4_p2 - e_7 * f_2_11 * r_12 * h2_0 + e_7 * fs_14_143_3 * r_12 * h2_p2 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m1, ph12_m1, ab_2, pc_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m1 = ph12_m1[k];

        pc_56[k] = - e_1 * fs_16065_32_7 * h2_m1 + e_2 * fs_855_8_210 * h4_m1 + e_2 * fs_6885_8_7 * r_2 * h2_m1 - e_3 * f_7875_11 * h6_m1 - e_3 * fs_4275_44_210 * r_2 * h4_m1 - e_3 * fs_3825_8_7 * r_4 * h2_m1 + e_4 * fs_12075_286_21 * h8_m1 + e_4 * f_4200_11 * r_2 * h6_m1 + e_4 * fs_4275_143_210 * r_4 * h4_m1 + e_4 * fs_1275_11_7 * r_6 * h2_m1 - e_5 * fs_4725_4199_1155 * h10_m1 - e_5 * fs_36225_2717_21 * r_2 * h8_m1 - e_5 * f_12600_187 * r_4 * h6_m1 - e_5 * fs_570_143_210 * r_6 * h4_m1 - e_5 * fs_3825_286_7 * r_8 * h2_m1 + e_6 * fs_58806_96577_182 * h12_m1 + e_6 * fs_18900_96577_1155 * r_2 * h10_m1 + e_6 * fs_3450_2717_21 * r_4 * h8_m1 + e_6 * f_16800_3553 * r_6 * h6_m1 + e_6 * fs_570_2431_210 * r_8 * h4_m1 + e_6 * fs_102_143_7 * r_10 * h2_m1 - e_7 * fs_4356_96577_182 * r_2 * h12_m1 - e_7 * fs_756_96577_1155 * r_4 * h10_m1 - e_7 * fs_100_2717_21 * r_6 * h8_m1 - e_7 * f_400_3553 * r_8 * h6_m1 - e_7 * fs_12_2431_210 * r_10 * h4_m1 - e_7 * fs_2_143_7 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ph12_m2, ab_2, pc_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m2 = ph12_m2[k];

        pc_57[k] = - e_1 * fs_112455_32_3 * h2_m2 + e_2 * fs_5985_4_5 * h4_m2 + e_2 * fs_48195_8_3 * r_2 * h2_m2 - e_3 * fs_1575_22_210 * h6_m2 - e_3 * fs_29925_22_5 * r_2 * h4_m2 - e_3 * fs_26775_8_3 * r_4 * h2_m2 + e_4 * fs_7245_286_70 * h8_m2 + e_4 * fs_420_11_210 * r_2 * h6_m2 + e_4 * fs_59850_143_5 * r_4 * h4_m2 + e_4 * fs_8925_11_3 * r_6 * h2_m2 - e_5 * fs_11025_4199_165 * h10_m2 - e_5 * fs_21735_2717_70 * r_2 * h8_m2 - e_5 * fs_1260_187_210 * r_4 * h6_m2 - e_5 * fs_7980_143_5 * r_6 * h4_m2 - e_5 * fs_26775_286_3 * r_8 * h2_m2 + e_6 * fs_10692_96577_3003 * h12_m2 + e_6 * fs_44100_96577_165 * r_2 * h10_m2 + e_6 * fs_2070_2717_70 * r_4 * h8_m2 + e_6 * fs_1680_3553_210 * r_6 * h6_m2 + e_6 * fs_7980_2431_5 * r_8 * h4_m2 + e_6 * fs_714_143_3 * r_10 * h2_m2 - e_7 * fs_792_96577_3003 * r_2 * h12_m2 - e_7 * fs_1764_96577_165 * r_4 * h10_m2 - e_7 * fs_60_2717_70 * r_6 * h8_m2 - e_7 * fs_40_3553_210 * r_8 * h6_m2 - e_7 * fs_168_2431_5 * r_10 * h4_m2 - e_7 * fs_14_143_3 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ph12_m3, ph12_m1, ab_2, pc_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m3 = ph12_m3[k];
        const auto h12_m1 = ph12_m1[k];

        pc_58[k] = e_1 * fs_16065_32_30 * h2_m1 + e_2 * fs_5985_16_7 * h4_m3 - e_2 * f_45315_16 * h4_m1 - e_2 * fs_6885_8_30 * r_2 * h2_m1 - e_3 * fs_4725_44_21 * h6_m3 + e_3 * fs_1575_22_210 * h6_m1 - e_3 * fs_29925_88_7 * r_2 * h4_m3 + e_3 * f_226575_88 * r_2 * h4_m1 + e_3 * fs_3825_8_30 * r_4 * h2_m1 + e_4 * fs_7245_1144_462 * h8_m3 - e_4 * fs_65205_1144_10 * h8_m1 + e_4 * fs_630_11_21 * r_2 * h6_m3 - e_4 * fs_420_11_210 * r_2 * h6_m1 + e_4 * fs_29925_286_7 * r_4 * h4_m3 - e_4 * f_226575_286 * r_4 * h4_m1 - e_4 * fs_1275_11_30 * r_6 * h2_m1 - e_5 * fs_11025_8398_429 * h10_m3 + e_5 * fs_9450_4199_22 * h10_m1 - e_5 * fs_21735_10868_462 * r_2 * h8_m3 + e_5 * fs_195615_10868_10 * r_2 * h8_m1 - e_5 * fs_1890_187_21 * r_4 * h6_m3 + e_5 * fs_1260_187_210 * r_4 * h6_m1 - e_5 * fs_1995_143_7 * r_6 * h4_m3 + e_5 * f_15105_143 * r_6 * h4_m1 + e_5 * fs_3825_286_30 * r_8 * h2_m1 + e_6 * fs_8019_96577_5005 * h12_m3 + e_6 * fs_29403_96577_195 * h12_m1 + e_6 * fs_22050_96577_429 * r_2 * h10_m3 - e_6 * fs_37800_96577_22 * r_2 * h10_m1 + e_6 * fs_1035_5434_462 * r_4 * h8_m3 - e_6 * fs_9315_5434_10 * r_4 * h8_m1 + e_6 * fs_2520_3553_21 * r_6 * h6_m3 - e_6 * fs_1680_3553_210 * r_6 * h6_m1 + e_6 * fs_1995_2431_7 * r_8 * h4_m3 - e_6 * f_15105_2431 * r_8 * h4_m1 - e_6 * fs_102_143_30 * r_10 * h2_m1 - e_7 * fs_594_96577_5005 * r_2 * h12_m3 - e_7 * fs_2178_96577_195 * r_2 * h12_m1 - e_7 * fs_882_96577_429 * r_4 * h10_m3 + e_7 * fs_1512_96577_22 * r_4 * h10_m1 - e_7 * fs_15_2717_462 * r_6 * h8_m3 + e_7 * fs_135_2717_10 * r_6 * h8_m1 - e_7 * fs_60_3553_21 * r_8 * h6_m3 + e_7 * fs_40_3553_210 * r_8 * h6_m1 - e_7 * fs_42_2431_7 * r_10 * h4_m3 + e_7 * f_318_2431 * r_10 * h4_m1 + e_7 * fs_2_143_30 * r_12 * h2_m1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ph12_m4, ph12_m2, ab_2, pc_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m4 = ph12_m4[k];
        const auto h12_m2 = ph12_m2[k];

        pc_59[k] = - e_1 * fs_16065_16_30 * h2_m2 - e_2 * fs_17955_16_14 * h4_m4 + e_2 * fs_2565_4_2 * h4_m2 + e_2 * fs_6885_4_30 * r_2 * h2_m2 + e_3 * fs_7875_88_70 * h6_m4 + e_3 * fs_4725_44_21 * h6_m2 + e_3 * fs_89775_88_14 * r_2 * h4_m4 - e_3 * fs_12825_22_2 * r_2 * h4_m2 - e_3 * fs_3825_4_30 * r_4 * h2_m2 - e_4 * fs_2415_1144_770 * h8_m4 - e_4 * fs_45885_572_7 * h8_m2 - e_4 * fs_525_11_70 * r_2 * h6_m4 - e_4 * fs_630_11_21 * r_2 * h6_m2 - e_4 * fs_89775_286_14 * r_4 * h4_m4 + e_4 * fs_25650_143_2 * r_4 * h4_m2 + e_4 * fs_2550_11_30 * r_6 * h2_m2 - e_5 * fs_4725_16796_858 * h10_m4 + e_5 * fs_1575_442_66 * h10_m2 + e_5 * fs_7245_10868_770 * r_2 * h8_m4 + e_5 * fs_7245_286_7 * r_2 * h8_m2 + e_5 * fs_1575_187_70 * r_4 * h6_m4 + e_5 * fs_1890_187_21 * r_4 * h6_m2 + e_5 * fs_5985_143_14 * r_6 * h4_m4 - e_5 * fs_3420_143_2 * r_6 * h4_m2 - e_5 * fs_3825_143_30 * r_8 * h2_m2 + e_6 * fs_7128_96577_5005 * h12_m4 + e_6 * fs_1782_96577_30030 * h12_m2 + e_6 * fs_4725_96577_858 * r_2 * h10_m4 - e_6 * fs_3150_5083_66 * r_2 * h10_m2 - e_6 * fs_345_5434_770 * r_4 * h8_m4 - e_6 * fs_345_143_7 * r_4 * h8_m2 - e_6 * fs_2100_3553_70 * r_6 * h6_m4 - e_6 * fs_2520_3553_21 * r_6 * h6_m2 - e_6 * fs_5985_2431_14 * r_8 * h4_m4 + e_6 * fs_3420_2431_2 * r_8 * h4_m2 + e_6 * fs_204_143_30 * r_10 * h2_m2 - e_7 * fs_528_96577_5005 * r_2 * h12_m4 - e_7 * fs_132_96577_30030 * r_2 * h12_m2 - e_7 * fs_189_96577_858 * r_4 * h10_m4 + e_7 * fs_126_5083_66 * r_4 * h10_m2 + e_7 * fs_5_2717_770 * r_6 * h8_m4 + e_7 * fs_10_143_7 * r_6 * h8_m2 + e_7 * fs_50_3553_70 * r_8 * h6_m4 + e_7 * fs_60_3553_21 * r_8 * h6_m2 + e_7 * fs_126_2431_14 * r_10 * h4_m4 - e_7 * fs_72_2431_2 * r_10 * h4_m2 - e_7 * fs_4_143_30 * r_12 * h2_m2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ph12_m5, ph12_m3, ab_2, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m3 = ph12_m3[k];

        pc_60[k] = e_2 * fs_4275_16_210 * h4_m3 + e_3 * fs_4725_88_462 * h6_m5 - e_3 * fs_7875_88_70 * h6_m3 - e_3 * fs_21375_88_210 * r_2 * h4_m3 - e_4 * fs_2415_572_3003 * h8_m5 - e_4 * fs_2415_572_385 * h8_m3 - e_4 * fs_315_11_462 * r_2 * h6_m5 + e_4 * fs_525_11_70 * r_2 * h6_m3 + e_4 * fs_21375_286_210 * r_4 * h4_m3 + e_5 * fs_17325_16796_286 * h10_m5 + e_5 * fs_17325_16796_1430 * h10_m3 + e_5 * fs_7245_5434_3003 * r_2 * h8_m5 + e_5 * fs_7245_5434_385 * r_2 * h8_m3 + e_5 * fs_945_187_462 * r_4 * h6_m5 - e_5 * fs_1575_187_70 * r_4 * h6_m3 - e_5 * fs_1425_143_210 * r_6 * h4_m3 + e_6 * fs_1782_96577_51051 * h12_m5 + e_6 * fs_2673_96577_6006 * h12_m3 - e_6 * fs_17325_96577_286 * r_2 * h10_m5 - e_6 * fs_17325_96577_1430 * r_2 * h10_m3 - e_6 * fs_345_2717_3003 * r_4 * h8_m5 - e_6 * fs_345_2717_385 * r_4 * h8_m3 - e_6 * fs_1260_3553_462 * r_6 * h6_m5 + e_6 * fs_2100_3553_70 * r_6 * h6_m3 + e_6 * fs_1425_2431_210 * r_8 * h4_m3 - e_7 * fs_132_96577_51051 * r_2 * h12_m5 - e_7 * fs_198_96577_6006 * r_2 * h12_m3 + e_7 * fs_693_96577_286 * r_4 * h10_m5 + e_7 * fs_693_96577_1430 * r_4 * h10_m3 + e_7 * fs_10_2717_3003 * r_6 * h8_m5 + e_7 * fs_10_2717_385 * r_6 * h8_m3 + e_7 * fs_30_3553_462 * r_8 * h6_m5 - e_7 * fs_50_3553_70 * r_8 * h6_m3 - e_7 * fs_30_2431_210 * r_10 * h4_m3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ph12_m6, ph12_m4, ab_2, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m4 = ph12_m4[k];

        pc_61[k] = - e_2 * fs_855_16_2310 * h4_m4 - e_3 * fs_1575_4_7 * h6_m6 - e_3 * fs_4725_88_462 * h6_m4 + e_3 * fs_4275_88_2310 * r_2 * h4_m4 - e_4 * fs_2415_52_13 * h8_m6 + e_4 * fs_2415_104_42 * h8_m4 + e_4 * fs_210_1_7 * r_2 * h6_m6 + e_4 * fs_315_11_462 * r_2 * h6_m4 - e_4 * fs_4275_286_2310 * r_4 * h4_m4 + e_5 * fs_1575_323_65 * h10_m6 + e_5 * fs_55125_16796_130 * h10_m4 + e_5 * fs_7245_494_13 * r_2 * h8_m6 - e_5 * fs_7245_988_42 * r_2 * h8_m4 - e_5 * fs_630_17_7 * r_4 * h6_m6 - e_5 * fs_945_187_462 * r_4 * h6_m4 + e_5 * fs_285_143_2310 * r_6 * h4_m4 + e_6 * fs_10692_96577_663 * h12_m6 + e_6 * fs_7128_96577_273 * h12_m4 - e_6 * fs_6300_7429_65 * r_2 * h10_m6 - e_6 * fs_55125_96577_130 * r_2 * h10_m4 - e_6 * fs_345_247_13 * r_4 * h8_m6 + e_6 * fs_345_494_42 * r_4 * h8_m4 + e_6 * fs_840_323_7 * r_6 * h6_m6 + e_6 * fs_1260_3553_462 * r_6 * h6_m4 - e_6 * fs_285_2431_2310 * r_8 * h4_m4 - e_7 * fs_792_96577_663 * r_2 * h12_m6 - e_7 * fs_528_96577_273 * r_2 * h12_m4 + e_7 * fs_252_7429_65 * r_4 * h10_m6 + e_7 * fs_2205_96577_130 * r_4 * h10_m4 + e_7 * fs_10_247_13 * r_6 * h8_m6 - e_7 * fs_5_247_42 * r_6 * h8_m4 - e_7 * fs_20_323_7 * r_8 * h6_m6 - e_7 * fs_30_3553_462 * r_8 * h6_m4 + e_7 * fs_6_2431_2310 * r_10 * h4_m4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, pe_7, ph6_m5, ph8_m7, ph8_m5, ph10_m7, ph10_m5, ph12_m7, ph12_m5, ab_2, pc_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h6_m5 = ph6_m5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m5 = ph12_m5[k];

        pc_62[k] = e_3 * fs_1575_4_7 * h6_m5 + e_4 * fs_2415_104_130 * h8_m7 + e_4 * fs_2415_104_182 * h8_m5 - e_4 * fs_210_1_7 * r_2 * h6_m5 + e_5 * fs_3150_4199_3315 * h10_m7 + e_5 * fs_33075_8398_39 * h10_m5 - e_5 * fs_7245_988_130 * r_2 * h8_m7 - e_5 * fs_7245_988_182 * r_2 * h8_m5 + e_5 * fs_630_17_7 * r_4 * h6_m5 + e_6 * fs_891_96577_25194 * h12_m7 + e_6 * fs_891_96577_3094 * h12_m5 - e_6 * fs_12600_96577_3315 * r_2 * h10_m7 - e_6 * fs_66150_96577_39 * r_2 * h10_m5 + e_6 * fs_345_494_130 * r_4 * h8_m7 + e_6 * fs_345_494_182 * r_4 * h8_m5 - e_6 * fs_840_323_7 * r_6 * h6_m5 - e_7 * fs_66_96577_25194 * r_2 * h12_m7 - e_7 * fs_66_96577_3094 * r_2 * h12_m5 + e_7 * fs_504_96577_3315 * r_4 * h10_m7 + e_7 * fs_2646_96577_39 * r_4 * h10_m5 - e_7 * fs_5_247_130 * r_6 * h8_m7 - e_7 * fs_5_247_182 * r_6 * h8_m5 + e_7 * fs_20_323_7 * r_8 * h6_m5;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ab_2, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];

        pc_63[k] = e_0 * f_155925_64 - e_1 * f_112455_16 * h2_0 - e_1 * f_363825_32 * r_2 + e_2 * f_17955_4 * h4_0 + e_2 * f_48195_4 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 - e_3 * f_15750_11 * h6_0 - e_3 * f_89775_22 * r_2 * h4_0 - e_3 * f_26775_4 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 + e_4 * f_84525_286 * h8_0 + e_4 * f_8400_11 * r_2 * h6_0 + e_4 * f_179550_143 * r_4 * h4_0 + e_4 * f_17850_11 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 - e_5 * f_198450_4199 * h10_0 - e_5 * f_253575_2717 * r_2 * h8_0 - e_5 * f_25200_187 * r_4 * h6_0 - e_5 * f_23940_143 * r_6 * h4_0 - e_5 * f_26775_143 * r_8 * h2_0 - e_5 * f_315_2 * r_10 + e_6 * f_823284_96577 * h12_0 + e_6 * f_793800_96577 * r_2 * h10_0 + e_6 * f_24150_2717 * r_4 * h8_0 + e_6 * f_33600_3553 * r_6 * h6_0 + e_6 * f_23940_2431 * r_8 * h4_0 + e_6 * f_1428_143 * r_10 * h2_0 + e_6 * f_105_13 * r_12 - e_7 * f_60984_96577 * r_2 * h12_0 - e_7 * f_31752_96577 * r_4 * h10_0 - e_7 * f_700_2717 * r_6 * h8_0 - e_7 * f_800_3553 * r_8 * h6_0 - e_7 * f_504_2431 * r_10 * h4_0 - e_7 * f_28_143 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ab_2, pc_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h12_p1 = ph12_p1[k];

        pc_64[k] = - e_1 * fs_16065_32_7 * h2_p1 + e_2 * fs_855_8_210 * h4_p1 + e_2 * fs_6885_8_7 * r_2 * h2_p1 - e_3 * f_7875_11 * h6_p1 - e_3 * fs_4275_44_210 * r_2 * h4_p1 - e_3 * fs_3825_8_7 * r_4 * h2_p1 + e_4 * fs_12075_286_21 * h8_p1 + e_4 * f_4200_11 * r_2 * h6_p1 + e_4 * fs_4275_143_210 * r_4 * h4_p1 + e_4 * fs_1275_11_7 * r_6 * h2_p1 - e_5 * fs_4725_4199_1155 * h10_p1 - e_5 * fs_36225_2717_21 * r_2 * h8_p1 - e_5 * f_12600_187 * r_4 * h6_p1 - e_5 * fs_570_143_210 * r_6 * h4_p1 - e_5 * fs_3825_286_7 * r_8 * h2_p1 + e_6 * fs_58806_96577_182 * h12_p1 + e_6 * fs_18900_96577_1155 * r_2 * h10_p1 + e_6 * fs_3450_2717_21 * r_4 * h8_p1 + e_6 * f_16800_3553 * r_6 * h6_p1 + e_6 * fs_570_2431_210 * r_8 * h4_p1 + e_6 * fs_102_143_7 * r_10 * h2_p1 - e_7 * fs_4356_96577_182 * r_2 * h12_p1 - e_7 * fs_756_96577_1155 * r_4 * h10_p1 - e_7 * fs_100_2717_21 * r_6 * h8_p1 - e_7 * f_400_3553 * r_8 * h6_p1 - e_7 * fs_12_2431_210 * r_10 * h4_p1 - e_7 * fs_2_143_7 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph12_p2, ab_2, pc_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h12_p2 = ph12_p2[k];

        pc_65[k] = e_1 * fs_16065_16_70 * h2_p2 - e_2 * fs_9405_16_42 * h4_p2 - e_2 * fs_6885_4_70 * r_2 * h2_p2 + e_3 * f_1575_2 * h6_p2 + e_3 * fs_4275_8_42 * r_2 * h4_p2 + e_3 * fs_3825_4_70 * r_4 * h2_p2 - e_4 * fs_16905_572_3 * h8_p2 - e_4 * f_420_1 * r_2 * h6_p2 - e_4 * fs_4275_26_42 * r_4 * h4_p2 - e_4 * fs_2550_11_70 * r_6 * h2_p2 - e_5 * fs_4725_4199_154 * h10_p2 + e_5 * fs_50715_5434_3 * r_2 * h8_p2 + e_5 * f_1260_17 * r_4 * h6_p2 + e_5 * fs_285_13_42 * r_6 * h4_p2 + e_5 * fs_3825_143_70 * r_8 * h2_p2 + e_6 * fs_18711_96577_1430 * h12_p2 + e_6 * fs_18900_96577_154 * r_2 * h10_p2 - e_6 * fs_2415_2717_3 * r_4 * h8_p2 - e_6 * f_1680_323 * r_6 * h6_p2 - e_6 * fs_285_221_42 * r_8 * h4_p2 - e_6 * fs_204_143_70 * r_10 * h2_p2 - e_7 * fs_1386_96577_1430 * r_2 * h12_p2 - e_7 * fs_756_96577_154 * r_4 * h10_p2 + e_7 * fs_70_2717_3 * r_6 * h8_p2 + e_7 * f_40_323 * r_8 * h6_p2 + e_7 * fs_6_221_42 * r_10 * h4_p2 + e_7 * fs_4_143_70 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ph10_p3, ph10_p4, ph12_p3, ph12_p4, ab_2, pc_66, pc_67 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p4 = ph12_p4[k];

        pc_66[k] = - e_2 * fs_17955_8_3 * h4_p3 + e_3 * f_67725_44 * h6_p3 + e_3 * fs_89775_44_3 * r_2 * h4_p3 - e_4 * fs_16905_286_22 * h8_p3 - e_4 * f_9030_11 * r_2 * h6_p3 - e_4 * fs_89775_143_3 * r_4 * h4_p3 + e_5 * fs_4725_8398_1001 * h10_p3 + e_5 * fs_50715_2717_22 * r_2 * h8_p3 + e_5 * f_27090_187 * r_4 * h6_p3 + e_5 * fs_11970_143_3 * r_6 * h4_p3 + e_6 * fs_12474_96577_2145 * h12_p3 - e_6 * fs_9450_96577_1001 * r_2 * h10_p3 - e_6 * fs_4830_2717_22 * r_4 * h8_p3 - e_6 * f_36120_3553 * r_6 * h6_p3 - e_6 * fs_11970_2431_3 * r_8 * h4_p3 - e_7 * fs_924_96577_2145 * r_2 * h12_p3 + e_7 * fs_378_96577_1001 * r_4 * h10_p3 + e_7 * fs_140_2717_22 * r_6 * h8_p3 + e_7 * f_860_3553 * r_8 * h6_p3 + e_7 * fs_252_2431_3 * r_10 * h4_p3;

        pc_67[k] = e_2 * fs_17955_8_5 * h4_p4 + e_3 * f_3150_11 * h6_p4 - e_3 * fs_89775_44_5 * r_2 * h4_p4 - e_4 * fs_50715_572_11 * h8_p4 - e_4 * f_1680_11 * r_2 * h6_p4 + e_4 * fs_89775_143_5 * r_4 * h4_p4 + e_5 * fs_1575_4199_15015 * h10_p4 + e_5 * fs_152145_5434_11 * r_2 * h8_p4 + e_5 * f_5040_187 * r_4 * h6_p4 - e_5 * fs_11970_143_5 * r_6 * h4_p4 + e_6 * fs_24948_96577_286 * h12_p4 - e_6 * fs_6300_96577_15015 * r_2 * h10_p4 - e_6 * fs_7245_2717_11 * r_4 * h8_p4 - e_6 * f_6720_3553 * r_6 * h6_p4 + e_6 * fs_11970_2431_5 * r_8 * h4_p4 - e_7 * fs_1848_96577_286 * r_2 * h12_p4 + e_7 * fs_252_96577_15015 * r_4 * h10_p4 + e_7 * fs_210_2717_11 * r_6 * h8_p4 + e_7 * f_160_3553 * r_8 * h6_p4 - e_7 * fs_252_2431_5 * r_10 * h4_p4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, pe_7, ph6_p5, ph6_p6, ph8_p6, ph10_p5, ph10_p6, ph12_p5, ph12_p6, ab_2, pc_68, pc_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p5 = ph12_p5[k];
        const auto h12_p6 = ph12_p6[k];

        pc_68[k] = - e_3 * f_7875_4 * h6_p5 + e_4 * f_1050_1 * r_2 * h6_p5 + e_5 * fs_1575_442_273 * h10_p5 - e_5 * f_3150_17 * r_4 * h6_p5 + e_6 * fs_12474_96577_442 * h12_p5 - e_6 * fs_3150_5083_273 * r_2 * h10_p5 + e_6 * f_4200_323 * r_6 * h6_p5 - e_7 * fs_924_96577_442 * r_2 * h12_p5 + e_7 * fs_126_5083_273 * r_4 * h10_p5 - e_7 * f_100_323 * r_8 * h6_p5;

        pc_69[k] = e_3 * f_1575_2 * h6_p6 + e_4 * fs_2415_52_91 * h8_p6 - e_4 * f_420_1 * r_2 * h6_p6 + e_5 * fs_9450_4199_455 * h10_p6 - e_5 * fs_7245_494_91 * r_2 * h8_p6 + e_5 * f_1260_17 * r_4 * h6_p6 + e_6 * fs_1782_96577_4641 * h12_p6 - e_6 * fs_37800_96577_455 * r_2 * h10_p6 + e_6 * fs_345_247_91 * r_4 * h8_p6 - e_6 * f_1680_323 * r_6 * h6_p6 - e_7 * fs_132_96577_4641 * r_2 * h12_p6 + e_7 * fs_1512_96577_455 * r_4 * h10_p6 - e_7 * fs_10_247_91 * r_6 * h8_p6 + e_7 * f_40_323 * r_8 * h6_p6;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ph12_0, ph12_p2, ab_2, pc_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

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
        const auto h12_0 = ph12_0[k];
        const auto h12_p2 = ph12_p2[k];

        pc_70[k] = e_0 * f_155925_64 - e_1 * f_208845_32 * h2_0 - e_1 * fs_112455_32_3 * h2_p2 - e_1 * f_363825_32 * r_2 + e_2 * f_3420_1 * h4_0 + e_2 * fs_5985_4_5 * h4_p2 + e_2 * f_89505_8 * r_2 * h2_0 + e_2 * fs_48195_8_3 * r_2 * h2_p2 + e_2 * f_218295_16 * r_4 - e_3 * f_7875_11 * h6_0 - e_3 * fs_1575_22_210 * h6_p2 - e_3 * f_34200_11 * r_2 * h4_0 - e_3 * fs_29925_22_5 * r_2 * h4_p2 - e_3 * f_49725_8 * r_4 * h2_0 - e_3 * fs_26775_8_3 * r_4 * h2_p2 - e_3 * f_51975_8 * r_6 + e_4 * f_12075_286 * h8_0 + e_4 * fs_7245_286_70 * h8_p2 + e_4 * f_4200_11 * r_2 * h6_0 + e_4 * fs_420_11_210 * r_2 * h6_p2 + e_4 * f_136800_143 * r_4 * h4_0 + e_4 * fs_59850_143_5 * r_4 * h4_p2 + e_4 * f_16575_11 * r_6 * h2_0 + e_4 * fs_8925_11_3 * r_6 * h2_p2 + e_4 * f_5775_4 * r_8 + e_5 * f_4725_323 * h10_0 - e_5 * fs_11025_4199_165 * h10_p2 - e_5 * f_36225_2717 * r_2 * h8_0 - e_5 * fs_21735_2717_70 * r_2 * h8_p2 - e_5 * f_12600_187 * r_4 * h6_0 - e_5 * fs_1260_187_210 * r_4 * h6_p2 - e_5 * f_18240_143 * r_6 * h4_0 - e_5 * fs_7980_143_5 * r_6 * h4_p2 - e_5 * f_3825_22 * r_8 * h2_0 - e_5 * fs_26775_286_3 * r_8 * h2_p2 - e_5 * f_315_2 * r_10 - e_6 * f_705672_96577 * h12_0 + e_6 * fs_10692_96577_3003 * h12_p2 - e_6 * f_18900_7429 * r_2 * h10_0 + e_6 * fs_44100_96577_165 * r_2 * h10_p2 + e_6 * f_3450_2717 * r_4 * h8_0 + e_6 * fs_2070_2717_70 * r_4 * h8_p2 + e_6 * f_16800_3553 * r_6 * h6_0 + e_6 * fs_1680_3553_210 * r_6 * h6_p2 + e_6 * f_18240_2431 * r_8 * h4_0 + e_6 * fs_7980_2431_5 * r_8 * h4_p2 + e_6 * f_102_11 * r_10 * h2_0 + e_6 * fs_714_143_3 * r_10 * h2_p2 + e_6 * f_105_13 * r_12 + e_7 * f_52272_96577 * r_2 * h12_0 - e_7 * fs_792_96577_3003 * r_2 * h12_p2 + e_7 * f_756_7429 * r_4 * h10_0 - e_7 * fs_1764_96577_165 * r_4 * h10_p2 - e_7 * f_100_2717 * r_6 * h8_0 - e_7 * fs_60_2717_70 * r_6 * h8_p2 - e_7 * f_400_3553 * r_8 * h6_0 - e_7 * fs_40_3553_210 * r_8 * h6_p2 - e_7 * f_384_2431 * r_10 * h4_0 - e_7 * fs_168_2431_5 * r_10 * h4_p2 - e_7 * f_2_11 * r_12 * h2_0 - e_7 * fs_14_143_3 * r_12 * h2_p2 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ph12_p1, ph12_p3, ab_2, pc_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p3 = ph12_p3[k];

        pc_71[k] = - e_1 * fs_16065_32_30 * h2_p1 + e_2 * f_45315_16 * h4_p1 + e_2 * fs_5985_16_7 * h4_p3 + e_2 * fs_6885_8_30 * r_2 * h2_p1 - e_3 * fs_1575_22_210 * h6_p1 - e_3 * fs_4725_44_21 * h6_p3 - e_3 * f_226575_88 * r_2 * h4_p1 - e_3 * fs_29925_88_7 * r_2 * h4_p3 - e_3 * fs_3825_8_30 * r_4 * h2_p1 + e_4 * fs_65205_1144_10 * h8_p1 + e_4 * fs_7245_1144_462 * h8_p3 + e_4 * fs_420_11_210 * r_2 * h6_p1 + e_4 * fs_630_11_21 * r_2 * h6_p3 + e_4 * f_226575_286 * r_4 * h4_p1 + e_4 * fs_29925_286_7 * r_4 * h4_p3 + e_4 * fs_1275_11_30 * r_6 * h2_p1 - e_5 * fs_9450_4199_22 * h10_p1 - e_5 * fs_11025_8398_429 * h10_p3 - e_5 * fs_195615_10868_10 * r_2 * h8_p1 - e_5 * fs_21735_10868_462 * r_2 * h8_p3 - e_5 * fs_1260_187_210 * r_4 * h6_p1 - e_5 * fs_1890_187_21 * r_4 * h6_p3 - e_5 * f_15105_143 * r_6 * h4_p1 - e_5 * fs_1995_143_7 * r_6 * h4_p3 - e_5 * fs_3825_286_30 * r_8 * h2_p1 - e_6 * fs_29403_96577_195 * h12_p1 + e_6 * fs_8019_96577_5005 * h12_p3 + e_6 * fs_37800_96577_22 * r_2 * h10_p1 + e_6 * fs_22050_96577_429 * r_2 * h10_p3 + e_6 * fs_9315_5434_10 * r_4 * h8_p1 + e_6 * fs_1035_5434_462 * r_4 * h8_p3 + e_6 * fs_1680_3553_210 * r_6 * h6_p1 + e_6 * fs_2520_3553_21 * r_6 * h6_p3 + e_6 * f_15105_2431 * r_8 * h4_p1 + e_6 * fs_1995_2431_7 * r_8 * h4_p3 + e_6 * fs_102_143_30 * r_10 * h2_p1 + e_7 * fs_2178_96577_195 * r_2 * h12_p1 - e_7 * fs_594_96577_5005 * r_2 * h12_p3 - e_7 * fs_1512_96577_22 * r_4 * h10_p1 - e_7 * fs_882_96577_429 * r_4 * h10_p3 - e_7 * fs_135_2717_10 * r_6 * h8_p1 - e_7 * fs_15_2717_462 * r_6 * h8_p3 - e_7 * fs_40_3553_210 * r_8 * h6_p1 - e_7 * fs_60_3553_21 * r_8 * h6_p3 - e_7 * f_318_2431 * r_10 * h4_p1 - e_7 * fs_42_2431_7 * r_10 * h4_p3 - e_7 * fs_2_143_30 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ph12_p2, ph12_p4, ab_2, pc_72 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p4 = ph12_p4[k];

        pc_72[k] = e_1 * fs_16065_16_30 * h2_p2 - e_2 * fs_2565_4_2 * h4_p2 - e_2 * fs_17955_16_14 * h4_p4 - e_2 * fs_6885_4_30 * r_2 * h2_p2 - e_3 * fs_4725_44_21 * h6_p2 + e_3 * fs_7875_88_70 * h6_p4 + e_3 * fs_12825_22_2 * r_2 * h4_p2 + e_3 * fs_89775_88_14 * r_2 * h4_p4 + e_3 * fs_3825_4_30 * r_4 * h2_p2 + e_4 * fs_45885_572_7 * h8_p2 - e_4 * fs_2415_1144_770 * h8_p4 + e_4 * fs_630_11_21 * r_2 * h6_p2 - e_4 * fs_525_11_70 * r_2 * h6_p4 - e_4 * fs_25650_143_2 * r_4 * h4_p2 - e_4 * fs_89775_286_14 * r_4 * h4_p4 - e_4 * fs_2550_11_30 * r_6 * h2_p2 - e_5 * fs_1575_442_66 * h10_p2 - e_5 * fs_4725_16796_858 * h10_p4 - e_5 * fs_7245_286_7 * r_2 * h8_p2 + e_5 * fs_7245_10868_770 * r_2 * h8_p4 - e_5 * fs_1890_187_21 * r_4 * h6_p2 + e_5 * fs_1575_187_70 * r_4 * h6_p4 + e_5 * fs_3420_143_2 * r_6 * h4_p2 + e_5 * fs_5985_143_14 * r_6 * h4_p4 + e_5 * fs_3825_143_30 * r_8 * h2_p2 - e_6 * fs_1782_96577_30030 * h12_p2 + e_6 * fs_7128_96577_5005 * h12_p4 + e_6 * fs_3150_5083_66 * r_2 * h10_p2 + e_6 * fs_4725_96577_858 * r_2 * h10_p4 + e_6 * fs_345_143_7 * r_4 * h8_p2 - e_6 * fs_345_5434_770 * r_4 * h8_p4 + e_6 * fs_2520_3553_21 * r_6 * h6_p2 - e_6 * fs_2100_3553_70 * r_6 * h6_p4 - e_6 * fs_3420_2431_2 * r_8 * h4_p2 - e_6 * fs_5985_2431_14 * r_8 * h4_p4 - e_6 * fs_204_143_30 * r_10 * h2_p2 + e_7 * fs_132_96577_30030 * r_2 * h12_p2 - e_7 * fs_528_96577_5005 * r_2 * h12_p4 - e_7 * fs_126_5083_66 * r_4 * h10_p2 - e_7 * fs_189_96577_858 * r_4 * h10_p4 - e_7 * fs_10_143_7 * r_6 * h8_p2 + e_7 * fs_5_2717_770 * r_6 * h8_p4 - e_7 * fs_60_3553_21 * r_8 * h6_p2 + e_7 * fs_50_3553_70 * r_8 * h6_p4 + e_7 * fs_72_2431_2 * r_10 * h4_p2 + e_7 * fs_126_2431_14 * r_10 * h4_p4 + e_7 * fs_4_143_30 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ph12_p3, ph12_p5, ab_2, pc_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p5 = ph12_p5[k];

        pc_73[k] = - e_2 * fs_4275_16_210 * h4_p3 + e_3 * fs_7875_88_70 * h6_p3 + e_3 * fs_4725_88_462 * h6_p5 + e_3 * fs_21375_88_210 * r_2 * h4_p3 + e_4 * fs_2415_572_385 * h8_p3 - e_4 * fs_2415_572_3003 * h8_p5 - e_4 * fs_525_11_70 * r_2 * h6_p3 - e_4 * fs_315_11_462 * r_2 * h6_p5 - e_4 * fs_21375_286_210 * r_4 * h4_p3 - e_5 * fs_17325_16796_1430 * h10_p3 + e_5 * fs_17325_16796_286 * h10_p5 - e_5 * fs_7245_5434_385 * r_2 * h8_p3 + e_5 * fs_7245_5434_3003 * r_2 * h8_p5 + e_5 * fs_1575_187_70 * r_4 * h6_p3 + e_5 * fs_945_187_462 * r_4 * h6_p5 + e_5 * fs_1425_143_210 * r_6 * h4_p3 - e_6 * fs_2673_96577_6006 * h12_p3 + e_6 * fs_1782_96577_51051 * h12_p5 + e_6 * fs_17325_96577_1430 * r_2 * h10_p3 - e_6 * fs_17325_96577_286 * r_2 * h10_p5 + e_6 * fs_345_2717_385 * r_4 * h8_p3 - e_6 * fs_345_2717_3003 * r_4 * h8_p5 - e_6 * fs_2100_3553_70 * r_6 * h6_p3 - e_6 * fs_1260_3553_462 * r_6 * h6_p5 - e_6 * fs_1425_2431_210 * r_8 * h4_p3 + e_7 * fs_198_96577_6006 * r_2 * h12_p3 - e_7 * fs_132_96577_51051 * r_2 * h12_p5 - e_7 * fs_693_96577_1430 * r_4 * h10_p3 + e_7 * fs_693_96577_286 * r_4 * h10_p5 - e_7 * fs_10_2717_385 * r_6 * h8_p3 + e_7 * fs_10_2717_3003 * r_6 * h8_p5 + e_7 * fs_50_3553_70 * r_8 * h6_p3 + e_7 * fs_30_3553_462 * r_8 * h6_p5 + e_7 * fs_30_2431_210 * r_10 * h4_p3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p4, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p6, ph12_p4, ph12_p6, ab_2, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p6 = ph12_p6[k];

        pc_74[k] = e_2 * fs_855_16_2310 * h4_p4 + e_3 * fs_4725_88_462 * h6_p4 - e_3 * fs_1575_4_7 * h6_p6 - e_3 * fs_4275_88_2310 * r_2 * h4_p4 - e_4 * fs_2415_104_42 * h8_p4 - e_4 * fs_2415_52_13 * h8_p6 - e_4 * fs_315_11_462 * r_2 * h6_p4 + e_4 * fs_210_1_7 * r_2 * h6_p6 + e_4 * fs_4275_286_2310 * r_4 * h4_p4 - e_5 * fs_55125_16796_130 * h10_p4 + e_5 * fs_1575_323_65 * h10_p6 + e_5 * fs_7245_988_42 * r_2 * h8_p4 + e_5 * fs_7245_494_13 * r_2 * h8_p6 + e_5 * fs_945_187_462 * r_4 * h6_p4 - e_5 * fs_630_17_7 * r_4 * h6_p6 - e_5 * fs_285_143_2310 * r_6 * h4_p4 - e_6 * fs_7128_96577_273 * h12_p4 + e_6 * fs_10692_96577_663 * h12_p6 + e_6 * fs_55125_96577_130 * r_2 * h10_p4 - e_6 * fs_6300_7429_65 * r_2 * h10_p6 - e_6 * fs_345_494_42 * r_4 * h8_p4 - e_6 * fs_345_247_13 * r_4 * h8_p6 - e_6 * fs_1260_3553_462 * r_6 * h6_p4 + e_6 * fs_840_323_7 * r_6 * h6_p6 + e_6 * fs_285_2431_2310 * r_8 * h4_p4 + e_7 * fs_528_96577_273 * r_2 * h12_p4 - e_7 * fs_792_96577_663 * r_2 * h12_p6 - e_7 * fs_2205_96577_130 * r_4 * h10_p4 + e_7 * fs_252_7429_65 * r_4 * h10_p6 + e_7 * fs_5_247_42 * r_6 * h8_p4 + e_7 * fs_10_247_13 * r_6 * h8_p6 + e_7 * fs_30_3553_462 * r_8 * h6_p4 - e_7 * fs_20_323_7 * r_8 * h6_p6 - e_7 * fs_6_2431_2310 * r_10 * h4_p4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, pe_7, ph6_p5, ph8_p5, ph8_p7, ph10_p5, ph10_p7, ph12_p5, ph12_p7, ab_2, pc_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h6_p5 = ph6_p5[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p5 = ph12_p5[k];
        const auto h12_p7 = ph12_p7[k];

        pc_75[k] = - e_3 * fs_1575_4_7 * h6_p5 - e_4 * fs_2415_104_182 * h8_p5 + e_4 * fs_2415_104_130 * h8_p7 + e_4 * fs_210_1_7 * r_2 * h6_p5 - e_5 * fs_33075_8398_39 * h10_p5 + e_5 * fs_3150_4199_3315 * h10_p7 + e_5 * fs_7245_988_182 * r_2 * h8_p5 - e_5 * fs_7245_988_130 * r_2 * h8_p7 - e_5 * fs_630_17_7 * r_4 * h6_p5 - e_6 * fs_891_96577_3094 * h12_p5 + e_6 * fs_891_96577_25194 * h12_p7 + e_6 * fs_66150_96577_39 * r_2 * h10_p5 - e_6 * fs_12600_96577_3315 * r_2 * h10_p7 - e_6 * fs_345_494_182 * r_4 * h8_p5 + e_6 * fs_345_494_130 * r_4 * h8_p7 + e_6 * fs_840_323_7 * r_6 * h6_p5 + e_7 * fs_66_96577_3094 * r_2 * h12_p5 - e_7 * fs_66_96577_25194 * r_2 * h12_p7 - e_7 * fs_2646_96577_39 * r_4 * h10_p5 + e_7 * fs_504_96577_3315 * r_4 * h10_p7 + e_7 * fs_5_247_182 * r_6 * h8_p5 - e_7 * fs_5_247_130 * r_6 * h8_p7 - e_7 * fs_20_323_7 * r_8 * h6_p5;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ph12_0, ph12_p4, ab_2, pc_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p4 = ph12_p4[k];

        pc_76[k] = e_0 * f_155925_64 - e_1 * f_80325_16 * h2_0 - e_1 * f_363825_32 * r_2 + e_2 * f_9405_16 * h4_0 + e_2 * fs_5985_8_35 * h4_p4 + e_2 * f_34425_4 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 + e_3 * f_1575_2 * h6_0 - e_3 * fs_4725_11_7 * h6_p4 - e_3 * f_4275_8 * r_2 * h4_0 - e_3 * fs_29925_44_35 * r_2 * h4_p4 - e_3 * f_19125_4 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 - e_4 * f_171465_572 * h8_0 + e_4 * fs_7245_286_77 * h8_p4 - e_4 * f_420_1 * r_2 * h6_0 + e_4 * fs_2520_11_7 * r_2 * h6_p4 + e_4 * f_4275_26 * r_4 * h4_0 + e_4 * fs_29925_143_35 * r_4 * h4_p4 + e_4 * f_12750_11 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 + e_5 * f_165375_4199 * h10_0 - e_5 * fs_3150_4199_2145 * h10_p4 + e_5 * f_514395_5434 * r_2 * h8_0 - e_5 * fs_21735_2717_77 * r_2 * h8_p4 + e_5 * f_1260_17 * r_4 * h6_0 - e_5 * fs_7560_187_7 * r_4 * h6_p4 - e_5 * f_285_13 * r_6 * h4_0 - e_5 * fs_3990_143_35 * r_6 * h4_p4 - e_5 * f_19125_143 * r_8 * h2_0 - e_5 * f_315_2 * r_10 + e_6 * f_441045_96577 * h12_0 + e_6 * fs_13365_96577_2002 * h12_p4 - e_6 * f_661500_96577 * r_2 * h10_0 + e_6 * fs_12600_96577_2145 * r_2 * h10_p4 - e_6 * f_24495_2717 * r_4 * h8_0 + e_6 * fs_2070_2717_77 * r_4 * h8_p4 - e_6 * f_1680_323 * r_6 * h6_0 + e_6 * fs_10080_3553_7 * r_6 * h6_p4 + e_6 * f_285_221 * r_8 * h4_0 + e_6 * fs_3990_2431_35 * r_8 * h4_p4 + e_6 * f_1020_143 * r_10 * h2_0 + e_6 * f_105_13 * r_12 - e_7 * f_32670_96577 * r_2 * h12_0 - e_7 * fs_990_96577_2002 * r_2 * h12_p4 + e_7 * f_26460_96577 * r_4 * h10_0 - e_7 * fs_504_96577_2145 * r_4 * h10_p4 + e_7 * f_710_2717 * r_6 * h8_0 - e_7 * fs_60_2717_77 * r_6 * h8_p4 + e_7 * f_40_323 * r_8 * h6_0 - e_7 * fs_240_3553_7 * r_8 * h6_p4 - e_7 * f_6_221 * r_10 * h4_0 - e_7 * fs_84_2431_35 * r_10 * h4_p4 - e_7 * f_20_143 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ph12_p1, ph12_p5, ab_2, pc_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p5 = ph12_p5[k];

        pc_77[k] = - e_1 * fs_80325_32_3 * h2_p1 + e_2 * fs_33345_32_10 * h4_p1 + e_2 * fs_34425_8_3 * r_2 * h2_p1 - e_3 * fs_4725_44_21 * h6_p1 - e_3 * fs_1575_44_154 * h6_p5 - e_3 * fs_166725_176_10 * r_2 * h4_p1 - e_3 * fs_19125_8_3 * r_4 * h2_p1 - e_4 * f_2415_22 * h8_p1 + e_4 * fs_2415_572_1001 * h8_p5 + e_4 * fs_630_11_21 * r_2 * h6_p1 + e_4 * fs_210_11_154 * r_2 * h6_p5 + e_4 * fs_12825_44_10 * r_4 * h4_p1 + e_4 * fs_6375_11_3 * r_6 * h2_p1 + e_5 * fs_42525_8398_55 * h10_p1 - e_5 * fs_7875_8398_858 * h10_p5 + e_5 * f_7245_209 * r_2 * h8_p1 - e_5 * fs_7245_5434_1001 * r_2 * h8_p5 - e_5 * fs_1890_187_21 * r_4 * h6_p1 - e_5 * fs_630_187_154 * r_4 * h6_p5 - e_5 * fs_855_22_10 * r_6 * h4_p1 - e_5 * fs_19125_286_3 * r_8 * h2_p1 + e_6 * fs_49005_193154_78 * h12_p1 + e_6 * fs_4455_96577_17017 * h12_p5 - e_6 * fs_85050_96577_55 * r_2 * h10_p1 + e_6 * fs_15750_96577_858 * r_2 * h10_p5 - e_6 * f_690_209 * r_4 * h8_p1 + e_6 * fs_345_2717_1001 * r_4 * h8_p5 + e_6 * fs_2520_3553_21 * r_6 * h6_p1 + e_6 * fs_840_3553_154 * r_6 * h6_p5 + e_6 * fs_855_374_10 * r_8 * h4_p1 + e_6 * fs_510_143_3 * r_10 * h2_p1 - e_7 * fs_1815_96577_78 * r_2 * h12_p1 - e_7 * fs_330_96577_17017 * r_2 * h12_p5 + e_7 * fs_3402_96577_55 * r_4 * h10_p1 - e_7 * fs_630_96577_858 * r_4 * h10_p5 + e_7 * f_20_209 * r_6 * h8_p1 - e_7 * fs_10_2717_1001 * r_6 * h8_p5 - e_7 * fs_60_3553_21 * r_8 * h6_p1 - e_7 * fs_20_3553_154 * r_8 * h6_p5 - e_7 * fs_9_187_10 * r_10 * h4_p1 - e_7 * fs_10_143_3 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ph12_p2, ph12_p6, ab_2, pc_78 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p6 = ph12_p6[k];

        pc_78[k] = e_1 * fs_48195_32_10 * h2_p2 + e_2 * fs_19665_32_6 * h4_p2 - e_2 * fs_20655_8_10 * r_2 * h2_p2 - e_3 * fs_4725_11_7 * h6_p2 + e_3 * fs_1575_22_385 * h6_p6 - e_3 * fs_98325_176_6 * r_2 * h4_p2 + e_3 * fs_11475_8_10 * r_4 * h2_p2 + e_4 * fs_2415_143_21 * h8_p2 - e_4 * fs_2415_572_715 * h8_p6 + e_4 * fs_2520_11_7 * r_2 * h6_p2 - e_4 * fs_420_11_385 * r_2 * h6_p6 + e_4 * fs_98325_572_6 * r_4 * h4_p2 - e_4 * fs_3825_11_10 * r_6 * h2_p2 + e_5 * fs_67725_8398_22 * h10_p2 - e_5 * fs_1575_4199_143 * h10_p6 - e_5 * fs_14490_2717_21 * r_2 * h8_p2 + e_5 * fs_7245_5434_715 * r_2 * h8_p6 - e_5 * fs_7560_187_7 * r_4 * h6_p2 + e_5 * fs_1260_187_385 * r_4 * h6_p6 - e_5 * fs_6555_286_6 * r_6 * h4_p2 + e_5 * fs_11475_286_10 * r_8 * h2_p2 + e_6 * fs_2673_193154_10010 * h12_p2 + e_6 * fs_2673_96577_36465 * h12_p6 - e_6 * fs_135450_96577_22 * r_2 * h10_p2 + e_6 * fs_6300_96577_143 * r_2 * h10_p6 + e_6 * fs_1380_2717_21 * r_4 * h8_p2 - e_6 * fs_345_2717_715 * r_4 * h8_p6 + e_6 * fs_10080_3553_7 * r_6 * h6_p2 - e_6 * fs_1680_3553_385 * r_6 * h6_p6 + e_6 * fs_6555_4862_6 * r_8 * h4_p2 - e_6 * fs_306_143_10 * r_10 * h2_p2 - e_7 * fs_99_96577_10010 * r_2 * h12_p2 - e_7 * fs_198_96577_36465 * r_2 * h12_p6 + e_7 * fs_5418_96577_22 * r_4 * h10_p2 - e_7 * fs_252_96577_143 * r_4 * h10_p6 - e_7 * fs_40_2717_21 * r_6 * h8_p2 + e_7 * fs_10_2717_715 * r_6 * h8_p6 - e_7 * fs_240_3553_7 * r_8 * h6_p2 + e_7 * fs_40_3553_385 * r_8 * h6_p6 - e_7 * fs_69_2431_6 * r_10 * h4_p2 + e_7 * fs_6_143_10 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ph12_p3, ph12_p7, ab_2, pc_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p7 = ph12_p7[k];

        pc_79[k] = - e_2 * fs_5985_32_462 * h4_p3 - e_3 * fs_1575_44_154 * h6_p3 + e_3 * fs_29925_176_462 * r_2 * h4_p3 + e_4 * fs_2415_26_7 * h8_p3 - e_4 * fs_2415_52_39 * h8_p7 + e_4 * fs_210_11_154 * r_2 * h6_p3 - e_4 * fs_29925_572_462 * r_4 * h4_p3 + e_5 * fs_48825_8398_26 * h10_p3 + e_5 * fs_11025_8398_442 * h10_p7 - e_5 * fs_7245_247_7 * r_2 * h8_p3 + e_5 * fs_7245_494_39 * r_2 * h8_p7 - e_5 * fs_630_187_154 * r_4 * h6_p3 + e_5 * fs_1995_286_462 * r_6 * h4_p3 + e_6 * fs_2673_193154_2730 * h12_p3 + e_6 * fs_2673_96577_20995 * h12_p7 - e_6 * fs_97650_96577_26 * r_2 * h10_p3 - e_6 * fs_22050_96577_442 * r_2 * h10_p7 + e_6 * fs_690_247_7 * r_4 * h8_p3 - e_6 * fs_345_247_39 * r_4 * h8_p7 + e_6 * fs_840_3553_154 * r_6 * h6_p3 - e_6 * fs_1995_4862_462 * r_8 * h4_p3 - e_7 * fs_99_96577_2730 * r_2 * h12_p3 - e_7 * fs_198_96577_20995 * r_2 * h12_p7 + e_7 * fs_3906_96577_26 * r_4 * h10_p3 + e_7 * fs_882_96577_442 * r_4 * h10_p7 - e_7 * fs_20_247_7 * r_6 * h8_p3 + e_7 * fs_10_247_39 * r_6 * h8_p7 - e_7 * fs_20_3553_154 * r_8 * h6_p3 + e_7 * fs_21_2431_462 * r_10 * h4_p3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p4, ph6_p4, ph8_p4, ph8_p8, ph10_p4, ph10_p8, ph12_p4, ph12_p8, ab_2, pc_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p8 = ph12_p8[k];

        pc_80[k] = e_2 * fs_2565_16_77 * h4_p4 + e_3 * fs_1575_22_385 * h6_p4 - e_3 * fs_12825_88_77 * r_2 * h4_p4 + e_4 * fs_2415_52_35 * h8_p4 + e_4 * fs_2415_52_13 * h8_p8 - e_4 * fs_420_11_385 * r_2 * h6_p4 + e_4 * fs_12825_286_77 * r_4 * h4_p4 + e_5 * fs_11025_4199_39 * h10_p4 + e_5 * fs_14175_4199_221 * h10_p8 - e_5 * fs_7245_494_35 * r_2 * h8_p4 - e_5 * fs_7245_494_13 * r_2 * h8_p8 + e_5 * fs_1260_187_385 * r_4 * h6_p4 - e_5 * fs_855_143_77 * r_6 * h4_p4 + e_6 * fs_891_96577_910 * h12_p4 + e_6 * fs_891_96577_62985 * h12_p8 - e_6 * fs_44100_96577_39 * r_2 * h10_p4 - e_6 * fs_56700_96577_221 * r_2 * h10_p8 + e_6 * fs_345_247_35 * r_4 * h8_p4 + e_6 * fs_345_247_13 * r_4 * h8_p8 - e_6 * fs_1680_3553_385 * r_6 * h6_p4 + e_6 * fs_855_2431_77 * r_8 * h4_p4 - e_7 * fs_66_96577_910 * r_2 * h12_p4 - e_7 * fs_66_96577_62985 * r_2 * h12_p8 + e_7 * fs_1764_96577_39 * r_4 * h10_p4 + e_7 * fs_2268_96577_221 * r_4 * h10_p8 - e_7 * fs_10_247_35 * r_6 * h8_p4 - e_7 * fs_10_247_13 * r_6 * h8_p8 + e_7 * fs_40_3553_385 * r_8 * h6_p4 - e_7 * fs_18_2431_77 * r_10 * h4_p4;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ph12_0, ph12_p6, ab_2, pc_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p6 = ph12_p6[k];

        pc_81[k] = e_0 * f_155925_64 - e_1 * f_80325_32 * h2_0 - e_1 * f_363825_32 * r_2 - e_2 * f_23085_8 * h4_0 + e_2 * f_34425_8 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 + e_3 * f_67725_44 * h6_0 - e_3 * fs_1575_22_462 * h6_p6 + e_3 * f_115425_44 * r_2 * h4_0 - e_3 * f_19125_8 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 - e_4 * f_45885_572 * h8_0 + e_4 * fs_2415_286_858 * h8_p6 - e_4 * f_9030_11 * r_2 * h6_0 + e_4 * fs_420_11_462 * r_2 * h6_p6 - e_4 * f_115425_143 * r_4 * h4_0 + e_4 * f_6375_11 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 - e_5 * f_23625_442 * h10_0 - e_5 * fs_4725_8398_4290 * h10_p6 + e_5 * f_7245_286 * r_2 * h8_0 - e_5 * fs_7245_2717_858 * r_2 * h8_p6 + e_5 * f_27090_187 * r_4 * h6_0 - e_5 * fs_1260_187_462 * r_4 * h6_p6 + e_5 * f_15390_143 * r_6 * h4_0 - e_5 * f_19125_286 * r_8 * h2_0 - e_5 * f_315_2 * r_10 - e_6 * f_196020_96577 * h12_0 + e_6 * fs_8910_96577_4862 * h12_p6 + e_6 * f_47250_5083 * r_2 * h10_0 + e_6 * fs_9450_96577_4290 * r_2 * h10_p6 - e_6 * f_345_143 * r_4 * h8_0 + e_6 * fs_690_2717_858 * r_4 * h8_p6 - e_6 * f_36120_3553 * r_6 * h6_0 + e_6 * fs_1680_3553_462 * r_6 * h6_p6 - e_6 * f_15390_2431 * r_8 * h4_0 + e_6 * f_510_143 * r_10 * h2_0 + e_6 * f_105_13 * r_12 + e_7 * f_14520_96577 * r_2 * h12_0 - e_7 * fs_660_96577_4862 * r_2 * h12_p6 - e_7 * f_1890_5083 * r_4 * h10_0 - e_7 * fs_378_96577_4290 * r_4 * h10_p6 + e_7 * f_10_143 * r_6 * h8_0 - e_7 * fs_20_2717_858 * r_6 * h8_p6 + e_7 * f_860_3553 * r_8 * h6_0 - e_7 * fs_40_3553_462 * r_8 * h6_p6 + e_7 * f_324_2431 * r_10 * h4_0 - e_7 * f_10_143 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ph12_p1, ph12_p7, ab_2, pc_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p7 = ph12_p7[k];

        pc_82[k] = - e_1 * fs_112455_64_10 * h2_p1 + e_2 * fs_17955_16_3 * h4_p1 + e_2 * fs_48195_16_10 * r_2 * h2_p1 + e_3 * fs_7875_88_70 * h6_p1 - e_3 * fs_89775_88_3 * r_2 * h4_p1 - e_3 * fs_26775_16_10 * r_4 * h2_p1 - e_4 * fs_21735_572_30 * h8_p1 + e_4 * fs_2415_572_858 * h8_p7 - e_4 * fs_525_11_70 * r_2 * h6_p1 + e_4 * fs_89775_286_3 * r_4 * h4_p1 + e_4 * fs_8925_22_10 * r_6 * h2_p1 - e_5 * fs_64575_16796_66 * h10_p1 - e_5 * fs_4725_8398_2431 * h10_p7 + e_5 * fs_65205_5434_30 * r_2 * h8_p1 - e_5 * fs_7245_5434_858 * r_2 * h8_p7 + e_5 * fs_1575_187_70 * r_4 * h6_p1 - e_5 * fs_5985_143_3 * r_6 * h4_p1 - e_5 * fs_26775_572_10 * r_8 * h2_p1 - e_6 * fs_9801_96577_65 * h12_p1 + e_6 * fs_891_96577_461890 * h12_p7 + e_6 * fs_64575_96577_66 * r_2 * h10_p1 + e_6 * fs_9450_96577_2431 * r_2 * h10_p7 - e_6 * fs_3105_2717_30 * r_4 * h8_p1 + e_6 * fs_345_2717_858 * r_4 * h8_p7 - e_6 * fs_2100_3553_70 * r_6 * h6_p1 + e_6 * fs_5985_2431_3 * r_8 * h4_p1 + e_6 * fs_357_143_10 * r_10 * h2_p1 + e_7 * fs_726_96577_65 * r_2 * h12_p1 - e_7 * fs_66_96577_461890 * r_2 * h12_p7 - e_7 * fs_2583_96577_66 * r_4 * h10_p1 - e_7 * fs_378_96577_2431 * r_4 * h10_p7 + e_7 * fs_90_2717_30 * r_6 * h8_p1 - e_7 * fs_10_2717_858 * r_6 * h8_p7 + e_7 * fs_50_3553_70 * r_8 * h6_p1 - e_7 * fs_126_2431_3 * r_10 * h4_p1 - e_7 * fs_7_143_10 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ph12_p2, ph12_p8, ab_2, pc_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p8 = ph12_p8[k];

        pc_83[k] = e_1 * fs_16065_32_55 * h2_p2 + e_2 * fs_2565_4_33 * h4_p2 - e_2 * fs_6885_8_55 * r_2 * h2_p2 - e_3 * fs_1575_44_154 * h6_p2 - e_3 * fs_12825_22_33 * r_2 * h4_p2 + e_3 * fs_3825_8_55 * r_4 * h2_p2 - e_4 * fs_7245_572_462 * h8_p2 - e_4 * fs_2415_52_39 * h8_p8 + e_4 * fs_210_11_154 * r_2 * h6_p2 + e_4 * fs_25650_143_33 * r_4 * h4_p2 - e_4 * fs_1275_11_55 * r_6 * h2_p2 - e_5 * f_174825_8398 * h10_p2 + e_5 * fs_1575_8398_663 * h10_p8 + e_5 * fs_21735_5434_462 * r_2 * h8_p2 + e_5 * fs_7245_494_39 * r_2 * h8_p8 - e_5 * fs_630_187_154 * r_4 * h6_p2 - e_5 * fs_3420_143_33 * r_6 * h4_p2 + e_5 * fs_3825_286_55 * r_8 * h2_p2 - e_6 * fs_1782_96577_455 * h12_p2 + e_6 * fs_3564_96577_20995 * h12_p8 + e_6 * f_349650_96577 * r_2 * h10_p2 - e_6 * fs_3150_96577_663 * r_2 * h10_p8 - e_6 * fs_1035_2717_462 * r_4 * h8_p2 - e_6 * fs_345_247_39 * r_4 * h8_p8 + e_6 * fs_840_3553_154 * r_6 * h6_p2 + e_6 * fs_3420_2431_33 * r_8 * h4_p2 - e_6 * fs_102_143_55 * r_10 * h2_p2 + e_7 * fs_132_96577_455 * r_2 * h12_p2 - e_7 * fs_264_96577_20995 * r_2 * h12_p8 - e_7 * f_13986_96577 * r_4 * h10_p2 + e_7 * fs_126_96577_663 * r_4 * h10_p8 + e_7 * fs_30_2717_462 * r_6 * h8_p2 + e_7 * fs_10_247_39 * r_6 * h8_p8 - e_7 * fs_20_3553_154 * r_8 * h6_p2 - e_7 * fs_72_2431_33 * r_10 * h4_p2 + e_7 * fs_2_143_55 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ph10_p9, ph12_p3, ph12_p9, ab_2, pc_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p9 = ph12_p9[k];

        pc_84[k] = - e_2 * fs_7695_32_154 * h4_p3 - e_3 * fs_1575_22_462 * h6_p3 + e_3 * fs_38475_176_154 * r_2 * h4_p3 - e_4 * fs_2415_52_21 * h8_p3 + e_4 * fs_420_11_462 * r_2 * h6_p3 - e_4 * fs_38475_572_154 * r_4 * h4_p3 - e_5 * fs_4725_4199_78 * h10_p3 + e_5 * fs_4725_8398_8398 * h10_p9 + e_5 * fs_7245_494_21 * r_2 * h8_p3 - e_5 * fs_1260_187_462 * r_4 * h6_p3 + e_5 * fs_2565_286_154 * r_6 * h4_p3 - e_6 * fs_891_193154_910 * h12_p3 + e_6 * fs_891_96577_146965 * h12_p9 + e_6 * fs_18900_96577_78 * r_2 * h10_p3 - e_6 * fs_9450_96577_8398 * r_2 * h10_p9 - e_6 * fs_345_247_21 * r_4 * h8_p3 + e_6 * fs_1680_3553_462 * r_6 * h6_p3 - e_6 * fs_2565_4862_154 * r_8 * h4_p3 + e_7 * fs_33_96577_910 * r_2 * h12_p3 - e_7 * fs_66_96577_146965 * r_2 * h12_p9 - e_7 * fs_756_96577_78 * r_4 * h10_p3 + e_7 * fs_378_96577_8398 * r_4 * h10_p9 + e_7 * fs_10_247_21 * r_6 * h8_p3 - e_7 * fs_40_3553_462 * r_8 * h6_p3 + e_7 * fs_27_2431_154 * r_10 * h4_p3;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ph12_0, ph12_p8, ab_2, pc_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p8 = ph12_p8[k];

        pc_85[k] = e_0 * f_155925_64 + e_1 * f_16065_16 * h2_0 - e_1 * f_363825_32 * r_2 - e_2 * f_5130_1 * h4_0 - e_2 * f_6885_4 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 + e_3 * f_3150_11 * h6_0 + e_3 * f_51300_11 * r_2 * h4_0 + e_3 * f_3825_4 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 + e_4 * f_214935_572 * h8_0 + e_4 * fs_7245_572_715 * h8_p8 - e_4 * f_1680_11 * r_2 * h6_0 - e_4 * f_205200_143 * r_4 * h4_0 - e_4 * f_2550_11 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 + e_5 * f_130725_4199 * h10_0 - e_5 * fs_1575_4199_12155 * h10_p8 - e_5 * f_644805_5434 * r_2 * h8_0 - e_5 * fs_21735_5434_715 * r_2 * h8_p8 + e_5 * f_5040_187 * r_4 * h6_0 + e_5 * f_27360_143 * r_6 * h4_0 + e_5 * f_3825_143 * r_8 * h2_0 - e_5 * f_315_2 * r_10 + e_6 * f_58806_96577 * h12_0 + e_6 * fs_1782_96577_138567 * h12_p8 - e_6 * f_522900_96577 * r_2 * h10_0 + e_6 * fs_6300_96577_12155 * r_2 * h10_p8 + e_6 * f_30705_2717 * r_4 * h8_0 + e_6 * fs_1035_2717_715 * r_4 * h8_p8 - e_6 * f_6720_3553 * r_6 * h6_0 - e_6 * f_27360_2431 * r_8 * h4_0 - e_6 * f_204_143 * r_10 * h2_0 + e_6 * f_105_13 * r_12 - e_7 * f_4356_96577 * r_2 * h12_0 - e_7 * fs_132_96577_138567 * r_2 * h12_p8 + e_7 * f_20916_96577 * r_4 * h10_0 - e_7 * fs_252_96577_12155 * r_4 * h10_p8 - e_7 * f_890_2717 * r_6 * h8_0 - e_7 * fs_30_2717_715 * r_6 * h8_p8 + e_7 * f_160_3553 * r_8 * h6_0 + e_7 * f_576_2431 * r_10 * h4_0 + e_7 * f_4_143 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ph12_p1, ph12_p9, ab_2, pc_86 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p9 = ph12_p9[k];

        pc_86[k] = - e_1 * fs_48195_64_66 * h2_p1 - e_2 * fs_2565_16_55 * h4_p1 + e_2 * fs_20655_16_66 * r_2 * h2_p1 + e_3 * fs_4725_88_462 * h6_p1 + e_3 * fs_12825_88_55 * r_2 * h4_p1 - e_3 * fs_11475_16_66 * r_4 * h2_p1 + e_4 * fs_7245_143_22 * h8_p1 - e_4 * fs_315_11_462 * r_2 * h6_p1 - e_4 * fs_12825_286_55 * r_4 * h4_p1 + e_4 * fs_3825_22_66 * r_6 * h2_p1 + e_5 * fs_67725_16796_10 * h10_p1 - e_5 * fs_1575_8398_20995 * h10_p9 - e_5 * fs_43470_2717_22 * r_2 * h8_p1 + e_5 * fs_945_187_462 * r_4 * h6_p1 + e_5 * fs_855_143_55 * r_6 * h4_p1 - e_5 * fs_11475_572_66 * r_8 * h2_p1 + e_6 * fs_891_96577_429 * h12_p1 + e_6 * fs_2673_96577_58786 * h12_p9 - e_6 * fs_67725_96577_10 * r_2 * h10_p1 + e_6 * fs_3150_96577_20995 * r_2 * h10_p9 + e_6 * fs_4140_2717_22 * r_4 * h8_p1 - e_6 * fs_1260_3553_462 * r_6 * h6_p1 - e_6 * fs_855_2431_55 * r_8 * h4_p1 + e_6 * fs_153_143_66 * r_10 * h2_p1 - e_7 * fs_66_96577_429 * r_2 * h12_p1 - e_7 * fs_198_96577_58786 * r_2 * h12_p9 + e_7 * fs_2709_96577_10 * r_4 * h10_p1 - e_7 * fs_126_96577_20995 * r_4 * h10_p9 - e_7 * fs_120_2717_22 * r_6 * h8_p1 + e_7 * fs_30_3553_462 * r_8 * h6_p1 + e_7 * fs_18_2431_55 * r_10 * h4_p1 - e_7 * fs_3_143_66 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph10_p10, ph12_p2, ph12_p10, ab_2, pc_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p10 = ph12_p10[k];

        pc_87[k] = e_1 * fs_16065_32_22 * h2_p2 + e_2 * fs_7695_32_330 * h4_p2 - e_2 * fs_6885_8_22 * r_2 * h2_p2 + e_3 * fs_1575_22_385 * h6_p2 - e_3 * fs_38475_176_330 * r_2 * h4_p2 + e_3 * fs_3825_8_22 * r_4 * h2_p2 + e_4 * fs_2415_572_1155 * h8_p2 - e_4 * fs_420_11_385 * r_2 * h6_p2 + e_4 * fs_38475_572_330 * r_4 * h4_p2 - e_4 * fs_1275_11_22 * r_6 * h2_p2 + e_5 * fs_14175_8398_10 * h10_p2 + e_5 * fs_1575_4199_12597 * h10_p10 - e_5 * fs_7245_5434_1155 * r_2 * h8_p2 + e_5 * fs_1260_187_385 * r_4 * h6_p2 - e_5 * fs_2565_286_330 * r_6 * h4_p2 + e_5 * fs_3825_286_22 * r_8 * h2_p2 + e_6 * fs_891_193154_182 * h12_p2 + e_6 * fs_891_96577_323323 * h12_p10 - e_6 * fs_28350_96577_10 * r_2 * h10_p2 - e_6 * fs_6300_96577_12597 * r_2 * h10_p10 + e_6 * fs_345_2717_1155 * r_4 * h8_p2 - e_6 * fs_1680_3553_385 * r_6 * h6_p2 + e_6 * fs_2565_4862_330 * r_8 * h4_p2 - e_6 * fs_102_143_22 * r_10 * h2_p2 - e_7 * fs_33_96577_182 * r_2 * h12_p2 - e_7 * fs_66_96577_323323 * r_2 * h12_p10 + e_7 * fs_1134_96577_10 * r_4 * h10_p2 + e_7 * fs_252_96577_12597 * r_4 * h10_p10 - e_7 * fs_10_2717_1155 * r_6 * h8_p2 + e_7 * fs_40_3553_385 * r_8 * h6_p2 - e_7 * fs_27_2431_330 * r_10 * h4_p2 + e_7 * fs_2_143_22 * r_12 * h2_p2;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ph12_0, ph12_p10, ab_2, pc_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p10 = ph12_p10[k];

        pc_88[k] = e_0 * f_155925_64 + e_1 * f_176715_32 * h2_0 - e_1 * f_363825_32 * r_2 - e_2 * f_28215_8 * h4_0 - e_2 * f_75735_8 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 - e_3 * f_7875_4 * h6_0 + e_3 * f_12825_4 * r_2 * h4_0 + e_3 * f_42075_8 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 - e_4 * f_12075_52 * h8_0 + e_4 * f_1050_1 * r_2 * h6_0 - e_4 * f_12825_13 * r_4 * h4_0 - e_4 * f_1275_1 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 - e_5 * f_77175_8398 * h10_0 - e_5 * fs_1575_8398_92378 * h10_p10 + e_5 * f_36225_494 * r_2 * h8_0 - e_5 * f_3150_17 * r_4 * h6_0 + e_5 * f_1710_13 * r_6 * h4_0 + e_5 * f_3825_26 * r_8 * h2_0 - e_5 * f_315_2 * r_10 - e_6 * f_10692_96577 * h12_0 + e_6 * fs_1782_96577_176358 * h12_p10 + e_6 * f_154350_96577 * r_2 * h10_0 + e_6 * fs_3150_96577_92378 * r_2 * h10_p10 - e_6 * f_1725_247 * r_4 * h8_0 + e_6 * f_4200_323 * r_6 * h6_0 - e_6 * f_1710_221 * r_8 * h4_0 - e_6 * f_102_13 * r_10 * h2_0 + e_6 * f_105_13 * r_12 + e_7 * f_792_96577 * r_2 * h12_0 - e_7 * fs_132_96577_176358 * r_2 * h12_p10 - e_7 * f_6174_96577 * r_4 * h10_0 - e_7 * fs_126_96577_92378 * r_4 * h10_p10 + e_7 * f_50_247 * r_6 * h8_0 - e_7 * f_100_323 * r_8 * h6_0 + e_7 * f_36_221 * r_10 * h4_0 + e_7 * f_2_13 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ph12_p11, ab_2, pc_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p11 = ph12_p11[k];

        pc_89[k] = - e_1 * f_176715_32 * h2_p1 - e_2 * fs_28215_32_30 * h4_p1 + e_2 * f_75735_8 * r_2 * h2_p1 - e_3 * fs_1575_4_7 * h6_p1 + e_3 * fs_12825_16_30 * r_2 * h4_p1 - e_3 * f_42075_8 * r_4 * h2_p1 - e_4 * fs_2415_52_3 * h8_p1 + e_4 * fs_210_1_7 * r_2 * h6_p1 - e_4 * fs_12825_52_30 * r_4 * h4_p1 + e_4 * f_1275_1 * r_6 * h2_p1 - e_5 * fs_1575_8398_165 * h10_p1 + e_5 * fs_7245_494_3 * r_2 * h8_p1 - e_5 * fs_630_17_7 * r_4 * h6_p1 + e_5 * fs_855_26_30 * r_6 * h4_p1 - e_5 * f_3825_26 * r_8 * h2_p1 - e_6 * fs_891_193154_26 * h12_p1 + e_6 * fs_891_96577_676039 * h12_p11 + e_6 * fs_3150_96577_165 * r_2 * h10_p1 - e_6 * fs_345_247_3 * r_4 * h8_p1 + e_6 * fs_840_323_7 * r_6 * h6_p1 - e_6 * fs_855_442_30 * r_8 * h4_p1 + e_6 * f_102_13 * r_10 * h2_p1 + e_7 * fs_33_96577_26 * r_2 * h12_p1 - e_7 * fs_66_96577_676039 * r_2 * h12_p11 - e_7 * fs_126_96577_165 * r_4 * h10_p1 + e_7 * fs_10_247_3 * r_6 * h8_p1 - e_7 * fs_20_323_7 * r_8 * h6_p1 + e_7 * fs_9_221_30 * r_10 * h4_p1 - e_7 * f_2_13 * r_12 * h2_p1;
    }

    // NOTE: the angular components are formed in 86 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, pe_7, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ph12_p12, ab_2, pc_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];
        const auto e_7 = pe_7[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;
        const auto r_14 = r_2 * r_12;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p12 = ph12_p12[k];

        pc_90[k] = e_0 * f_155925_64 + e_1 * f_176715_16 * h2_0 - e_1 * f_363825_32 * r_2 + e_2 * f_84645_16 * h4_0 - e_2 * f_75735_4 * r_2 * h2_0 + e_2 * f_218295_16 * r_4 + e_3 * f_1575_2 * h6_0 - e_3 * f_38475_8 * r_2 * h4_0 + e_3 * f_42075_4 * r_4 * h2_0 - e_3 * f_51975_8 * r_6 + e_4 * f_2415_52 * h8_0 - e_4 * f_420_1 * r_2 * h6_0 + e_4 * f_38475_26 * r_4 * h4_0 - e_4 * f_2550_1 * r_6 * h2_0 + e_4 * f_5775_4 * r_8 + e_5 * f_4725_4199 * h10_0 - e_5 * f_7245_494 * r_2 * h8_0 + e_5 * f_1260_17 * r_4 * h6_0 - e_5 * f_2565_13 * r_6 * h4_0 + e_5 * f_3825_13 * r_8 * h2_0 - e_5 * f_315_2 * r_10 + e_6 * f_891_96577 * h12_0 + e_6 * fs_891_96577_1352078 * h12_p12 - e_6 * f_18900_96577 * r_2 * h10_0 + e_6 * f_345_247 * r_4 * h8_0 - e_6 * f_1680_323 * r_6 * h6_0 + e_6 * f_2565_221 * r_8 * h4_0 - e_6 * f_204_13 * r_10 * h2_0 + e_6 * f_105_13 * r_12 - e_7 * f_66_96577 * r_2 * h12_0 - e_7 * fs_66_96577_1352078 * r_2 * h12_p12 + e_7 * f_756_96577 * r_4 * h10_0 - e_7 * f_10_247 * r_6 * h8_0 + e_7 * f_40_323 * r_8 * h6_0 - e_7 * f_54_221 * r_10 * h4_0 + e_7 * f_4_13 * r_12 * h2_0 - e_7 * f_2_13 * r_14;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[169] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 2, 15, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 3, 16, 29, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 4, 17, 30, 43, 56, 57, 58, 59, 60, 61, 62, 63, 64, 5, 18, 31, 44, 57, 70, 71, 72, 73, 74, 75, 76, 77, 6, 19, 32, 45, 58, 71, 84, 85, 86, 87, 88, 89, 90, 7, 20, 33, 46, 59, 72, 85, 98, 99, 100, 101, 102, 103, 8, 21, 34, 47, 60, 73, 86, 99, 112, 113, 114, 115, 116, 9, 22, 35, 48, 61, 74, 87, 100, 113, 126, 127, 128, 129, 10, 23, 36, 49, 62, 75, 88, 101, 114, 127, 140, 141, 142, 11, 24, 37, 50, 63, 76, 89, 102, 115, 128, 141, 154, 155, 12, 25, 38, 51, 64, 77, 90, 103, 116, 129, 142, 155, 168};

    for (size_t n = 0; n < 169; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
