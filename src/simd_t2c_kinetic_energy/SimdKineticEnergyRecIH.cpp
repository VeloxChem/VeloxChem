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



#include "SimdKineticEnergyRecIH.hpp"

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
compute_ih_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ih_kinetic_energy: Basis functions must be of angular momenta six and five"));
    }

    if (harmonics.size() < 11)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ih_kinetic_energy: Harmonics must reach angular momentum 11"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ih_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 143 * nvalues, 0.0);

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

            const auto ff_0 = fbase * bexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * bexp * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * bexp * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_6 = fbase * bexp * fmu * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

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
    const auto *ph11_m11 = harmonics[10].data(0);
    const auto *ph11_m10 = harmonics[10].data(1);
    const auto *ph11_m9 = harmonics[10].data(2);
    const auto *ph11_m8 = harmonics[10].data(3);
    const auto *ph11_m7 = harmonics[10].data(4);
    const auto *ph11_m6 = harmonics[10].data(5);
    const auto *ph11_m5 = harmonics[10].data(6);
    const auto *ph11_m4 = harmonics[10].data(7);
    const auto *ph11_m3 = harmonics[10].data(8);
    const auto *ph11_m2 = harmonics[10].data(9);
    const auto *ph11_m1 = harmonics[10].data(10);
    const auto *ph11_0 = harmonics[10].data(11);
    const auto *ph11_p1 = harmonics[10].data(12);
    const auto *ph11_p2 = harmonics[10].data(13);
    const auto *ph11_p3 = harmonics[10].data(14);
    const auto *ph11_p4 = harmonics[10].data(15);
    const auto *ph11_p5 = harmonics[10].data(16);
    const auto *ph11_p6 = harmonics[10].data(17);
    const auto *ph11_p7 = harmonics[10].data(18);
    const auto *ph11_p8 = harmonics[10].data(19);
    const auto *ph11_p9 = harmonics[10].data(20);
    const auto *ph11_p10 = harmonics[10].data(21);
    const auto *ph11_p11 = harmonics[10].data(22);

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
    auto *pc_99 = values + 99 * nvalues;
    auto *pc_100 = values + 100 * nvalues;
    auto *pc_101 = values + 101 * nvalues;
    auto *pc_102 = values + 102 * nvalues;
    auto *pc_103 = values + 103 * nvalues;
    auto *pc_104 = values + 104 * nvalues;
    auto *pc_105 = values + 105 * nvalues;
    auto *pc_106 = values + 106 * nvalues;
    auto *pc_107 = values + 107 * nvalues;
    auto *pc_108 = values + 108 * nvalues;
    auto *pc_109 = values + 109 * nvalues;
    auto *pc_110 = values + 110 * nvalues;
    auto *pc_111 = values + 111 * nvalues;
    auto *pc_112 = values + 112 * nvalues;
    auto *pc_113 = values + 113 * nvalues;
    auto *pc_114 = values + 114 * nvalues;
    auto *pc_115 = values + 115 * nvalues;
    auto *pc_116 = values + 116 * nvalues;
    auto *pc_117 = values + 117 * nvalues;
    auto *pc_118 = values + 118 * nvalues;
    auto *pc_119 = values + 119 * nvalues;
    auto *pc_120 = values + 120 * nvalues;
    auto *pc_121 = values + 121 * nvalues;
    auto *pc_122 = values + 122 * nvalues;
    auto *pc_123 = values + 123 * nvalues;
    auto *pc_124 = values + 124 * nvalues;
    auto *pc_125 = values + 125 * nvalues;
    auto *pc_126 = values + 126 * nvalues;
    auto *pc_127 = values + 127 * nvalues;
    auto *pc_128 = values + 128 * nvalues;
    auto *pc_129 = values + 129 * nvalues;
    auto *pc_130 = values + 130 * nvalues;
    auto *pc_131 = values + 131 * nvalues;
    auto *pc_132 = values + 132 * nvalues;
    auto *pc_133 = values + 133 * nvalues;
    auto *pc_134 = values + 134 * nvalues;
    auto *pc_135 = values + 135 * nvalues;
    auto *pc_136 = values + 136 * nvalues;
    auto *pc_137 = values + 137 * nvalues;
    auto *pc_138 = values + 138 * nvalues;
    auto *pc_139 = values + 139 * nvalues;
    auto *pc_140 = values + 140 * nvalues;
    auto *pc_141 = values + 141 * nvalues;
    auto *pc_142 = values + 142 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_101430_2431 = 101430.0 / 2431.0;
    const auto f_1102_13 = 1102.0 / 13.0;
    const auto f_1140_1 = 1140.0;
    const auto f_116_663 = 116.0 / 663.0;
    const auto f_1190_143 = 1190.0 / 143.0;
    const auto f_12600_46189 = 12600.0 / 46189.0;
    const auto f_1350_1 = 1350.0;
    const auto f_1368_13 = 1368.0 / 13.0;
    const auto f_1425_4 = 1425.0 / 4.0;
    const auto f_14_143 = 14.0 / 143.0;
    const auto f_1520_13 = 1520.0 / 13.0;
    const auto f_16065_16 = 16065.0 / 16.0;
    const auto f_16065_8 = 16065.0 / 8.0;
    const auto f_160_663 = 160.0 / 663.0;
    const auto f_1620_143 = 1620.0 / 143.0;
    const auto f_16_221 = 16.0 / 221.0;
    const auto f_17850_143 = 17850.0 / 143.0;
    const auto f_18225_4 = 18225.0 / 4.0;
    const auto f_1824_13 = 1824.0 / 13.0;
    const auto f_1824_221 = 1824.0 / 221.0;
    const auto f_19320_2431 = 19320.0 / 2431.0;
    const auto f_198450_2431 = 198450.0 / 2431.0;
    const auto f_2025_11 = 2025.0 / 11.0;
    const auto f_20825_143 = 20825.0 / 143.0;
    const auto f_20825_22 = 20825.0 / 22.0;
    const auto f_20825_8 = 20825.0 / 8.0;
    const auto f_20_221 = 20.0 / 221.0;
    const auto f_2280_13 = 2280.0 / 13.0;
    const auto f_2432_221 = 2432.0 / 221.0;
    const auto f_25515_4 = 25515.0 / 4.0;
    const auto f_2755_4 = 2755.0 / 4.0;
    const auto f_276_2431 = 276.0 / 2431.0;
    const auto f_2772_4199 = 2772.0 / 4199.0;
    const auto f_2850_13 = 2850.0 / 13.0;
    const auto f_285_1 = 285.0;
    const auto f_28_143 = 28.0 / 143.0;
    const auto f_33075_143 = 33075.0 / 143.0;
    const auto f_33327_2431 = 33327.0 / 2431.0;
    const auto f_34650_4199 = 34650.0 / 4199.0;
    const auto f_36_143 = 36.0 / 143.0;
    const auto f_37485_16 = 37485.0 / 16.0;
    const auto f_396900_46189 = 396900.0 / 46189.0;
    const auto f_4165_429 = 4165.0 / 429.0;
    const auto f_42525_16 = 42525.0 / 16.0;
    const auto f_4408_663 = 4408.0 / 663.0;
    const auto f_456_13 = 456.0 / 13.0;
    const auto f_48_221 = 48.0 / 221.0;
    const auto f_5510_13 = 5510.0 / 13.0;
    const auto f_570_13 = 570.0 / 13.0;
    const auto f_595_143 = 595.0 / 143.0;
    const auto f_6080_663 = 6080.0 / 663.0;
    const auto f_608_221 = 608.0 / 221.0;
    const auto f_6348_2431 = 6348.0 / 2431.0;
    const auto f_64_221 = 64.0 / 221.0;
    const auto f_6840_13 = 6840.0 / 13.0;
    const auto f_7600_13 = 7600.0 / 13.0;
    const auto f_760_221 = 760.0 / 221.0;
    const auto f_840_2431 = 840.0 / 2431.0;
    const auto f_855_1 = 855.0;
    const auto f_8925_11 = 8925.0 / 11.0;
    const auto f_8925_143 = 8925.0 / 143.0;
    const auto f_8925_22 = 8925.0 / 22.0;
    const auto f_8925_4 = 8925.0 / 4.0;
    const auto f_8925_8 = 8925.0 / 8.0;
    const auto f_9120_13 = 9120.0 / 13.0;
    const auto f_950_1 = 950.0;
    const auto f_98_429 = 98.0 / 429.0;
    const auto fs_100_663_3 = std::sqrt(10000.0 / 146523.0);
    const auto fs_10143_2431_143 = std::sqrt(102880449.0 / 41327.0);
    const auto fs_10143_4862_165 = std::sqrt(1543206735.0 / 2149004.0);
    const auto fs_1035_2431_10 = std::sqrt(10712250.0 / 5909761.0);
    const auto fs_1035_46189_14 = std::sqrt(14997150.0 / 2133423721.0);
    const auto fs_1045_13_42 = std::sqrt(45865050.0 / 169.0);
    const auto fs_1045_4_5 = std::sqrt(5460125.0 / 16.0);
    const auto fs_1045_8_42 = std::sqrt(22932525.0 / 32.0);
    const auto fs_104895_46189_10 = std::sqrt(110029610250.0 / 2133423721.0);
    const auto fs_104895_4862_10 = std::sqrt(55014805125.0 / 11819522.0);
    const auto fs_1050_4199_429 = std::sqrt(36382500.0 / 1356277.0);
    const auto fs_1050_46189_154 = std::sqrt(15435000.0 / 193947611.0);
    const auto fs_105_4199_130 = std::sqrt(110250.0 / 1356277.0);
    const auto fs_105_4199_286 = std::sqrt(242550.0 / 1356277.0);
    const auto fs_106_2431_33 = std::sqrt(33708.0 / 537251.0);
    const auto fs_1080_143_2 = std::sqrt(2332800.0 / 20449.0);
    const auto fs_10_221_30 = std::sqrt(3000.0 / 48841.0);
    const auto fs_10_221_91 = std::sqrt(700.0 / 3757.0);
    const auto fs_10_2431_7293 = std::sqrt(300.0 / 2431.0);
    const auto fs_10_2431_858 = std::sqrt(600.0 / 41327.0);
    const auto fs_11025_286_11 = std::sqrt(121550625.0 / 7436.0);
    const auto fs_11025_286_14 = std::sqrt(850854375.0 / 40898.0);
    const auto fs_11025_286_35 = std::sqrt(4254271875.0 / 81796.0);
    const auto fs_11025_286_66 = std::sqrt(364651875.0 / 3718.0);
    const auto fs_11025_572_154 = std::sqrt(850854375.0 / 14872.0);
    const auto fs_1104_2431_154 = std::sqrt(17063424.0 / 537251.0);
    const auto fs_1104_2431_165 = std::sqrt(18282240.0 / 537251.0);
    const auto fs_1110_46189_55 = std::sqrt(6160500.0 / 193947611.0);
    const auto fs_1125_4199_231 = std::sqrt(292359375.0 / 17631601.0);
    const auto fs_1125_4199_429 = std::sqrt(41765625.0 / 1356277.0);
    const auto fs_1125_8398_2002 = std::sqrt(97453125.0 / 2712554.0);
    const auto fs_114_13_154 = std::sqrt(2001384.0 / 169.0);
    const auto fs_114_13_30 = std::sqrt(389880.0 / 169.0);
    const auto fs_114_13_55 = std::sqrt(714780.0 / 169.0);
    const auto fs_114_13_70 = std::sqrt(909720.0 / 169.0);
    const auto fs_1150_2431_143 = std::sqrt(1322500.0 / 41327.0);
    const auto fs_115605_1144_2 = std::sqrt(13364516025.0 / 654368.0);
    const auto fs_11655_572_55 = std::sqrt(679195125.0 / 29744.0);
    const auto fs_11900_143_2 = std::sqrt(283220000.0 / 20449.0);
    const auto fs_119070_46189_3 = std::sqrt(42532994700.0 / 2133423721.0);
    const auto fs_1190_429_2 = std::sqrt(2832200.0 / 184041.0);
    const auto fs_1190_429_21 = std::sqrt(9912700.0 / 61347.0);
    const auto fs_1190_429_7 = std::sqrt(9912700.0 / 184041.0);
    const auto fs_12075_4862_143 = std::sqrt(145805625.0 / 165308.0);
    const auto fs_120_4199_286 = std::sqrt(316800.0 / 1356277.0);
    const auto fs_120_4199_33 = std::sqrt(475200.0 / 17631601.0);
    const auto fs_120_4199_65 = std::sqrt(72000.0 / 1356277.0);
    const auto fs_12190_2431_3 = std::sqrt(445788300.0 / 5909761.0);
    const auto fs_1235_8_6 = std::sqrt(4575675.0 / 32.0);
    const auto fs_12375_4199_2 = std::sqrt(306281250.0 / 17631601.0);
    const auto fs_12420_2431_2 = std::sqrt(308512800.0 / 5909761.0);
    const auto fs_1260_13_3 = std::sqrt(4762800.0 / 169.0);
    const auto fs_1260_143_35 = std::sqrt(55566000.0 / 20449.0);
    const auto fs_1260_143_385 = std::sqrt(55566000.0 / 1859.0);
    const auto fs_1260_46189_35 = std::sqrt(55566000.0 / 2133423721.0);
    const auto fs_1260_46189_70 = std::sqrt(111132000.0 / 2133423721.0);
    const auto fs_126_4199_143 = std::sqrt(174636.0 / 1356277.0);
    const auto fs_127995_4862_3 = std::sqrt(49148160075.0 / 23639044.0);
    const auto fs_12_143_5 = std::sqrt(720.0 / 20449.0);
    const auto fs_12_2431_2431 = std::sqrt(144.0 / 2431.0);
    const auto fs_12_2431_715 = std::sqrt(720.0 / 41327.0);
    const auto fs_12_4199_15 = std::sqrt(2160.0 / 17631601.0);
    const auto fs_12_4199_15015 = std::sqrt(166320.0 / 1356277.0);
    const auto fs_12_4199_1547 = std::sqrt(1008.0 / 79781.0);
    const auto fs_12_4199_2002 = std::sqrt(22176.0 / 1356277.0);
    const auto fs_12_4199_2431 = std::sqrt(1584.0 / 79781.0);
    const auto fs_12_4199_24310 = std::sqrt(15840.0 / 79781.0);
    const auto fs_12_4199_41990 = std::sqrt(1440.0 / 4199.0);
    const auto fs_12_4199_7735 = std::sqrt(5040.0 / 79781.0);
    const auto fs_1311_2431_154 = std::sqrt(24062094.0 / 537251.0);
    const auto fs_13230_2717_5 = std::sqrt(875164500.0 / 7382089.0);
    const auto fs_13230_46189_165 = std::sqrt(2625493500.0 / 193947611.0);
    const auto fs_13230_46189_21 = std::sqrt(3675690900.0 / 2133423721.0);
    const auto fs_13230_46189_35 = std::sqrt(6126151500.0 / 2133423721.0);
    const auto fs_1334_2431_143 = std::sqrt(1779556.0 / 41327.0);
    const auto fs_1350_11_2 = std::sqrt(3645000.0 / 121.0);
    const auto fs_135_143_110 = std::sqrt(182250.0 / 1859.0);
    const auto fs_135_143_2 = std::sqrt(36450.0 / 20449.0);
    const auto fs_135_143_6 = std::sqrt(109350.0 / 20449.0);
    const auto fs_135_4199_22 = std::sqrt(400950.0 / 17631601.0);
    const auto fs_138_221_221 = std::sqrt(19044.0 / 221.0);
    const auto fs_138_221_91 = std::sqrt(133308.0 / 3757.0);
    const auto fs_138_2431_10010 = std::sqrt(1333080.0 / 41327.0);
    const auto fs_138_2431_12155 = std::sqrt(95220.0 / 2431.0);
    const auto fs_14007_4862_143 = std::sqrt(196196049.0 / 165308.0);
    const auto fs_14175_16_5 = std::sqrt(1004653125.0 / 256.0);
    const auto fs_14175_32_11 = std::sqrt(2210236875.0 / 1024.0);
    const auto fs_14175_32_14 = std::sqrt(1406514375.0 / 512.0);
    const auto fs_14175_32_15 = std::sqrt(3013959375.0 / 1024.0);
    const auto fs_14175_32_21 = std::sqrt(4219543125.0 / 1024.0);
    const auto fs_14175_32_3 = std::sqrt(602791875.0 / 1024.0);
    const auto fs_14175_32_33 = std::sqrt(6630710625.0 / 1024.0);
    const auto fs_14175_32_35 = std::sqrt(7032571875.0 / 1024.0);
    const auto fs_14175_32_5 = std::sqrt(1004653125.0 / 1024.0);
    const auto fs_14175_4199_13 = std::sqrt(200930625.0 / 1356277.0);
    const auto fs_14175_442_13 = std::sqrt(200930625.0 / 15028.0);
    const auto fs_14175_46189_385 = std::sqrt(7032571875.0 / 193947611.0);
    const auto fs_14175_4862_385 = std::sqrt(7032571875.0 / 2149004.0);
    const auto fs_14175_64_110 = std::sqrt(11051184375.0 / 2048.0);
    const auto fs_14175_64_2 = std::sqrt(200930625.0 / 2048.0);
    const auto fs_14175_64_6 = std::sqrt(602791875.0 / 2048.0);
    const auto fs_14175_8398_30 = std::sqrt(3013959375.0 / 35263202.0);
    const auto fs_14175_884_30 = std::sqrt(3013959375.0 / 390728.0);
    const auto fs_14175_8_2 = std::sqrt(200930625.0 / 32.0);
    const auto fs_1425_13_30 = std::sqrt(60918750.0 / 169.0);
    const auto fs_1425_4_11 = std::sqrt(22336875.0 / 16.0);
    const auto fs_1425_4_7 = std::sqrt(14214375.0 / 16.0);
    const auto fs_1425_8_30 = std::sqrt(30459375.0 / 32.0);
    const auto fs_1449_143_5 = std::sqrt(10498005.0 / 20449.0);
    const auto fs_1449_2431_2431 = std::sqrt(2099601.0 / 2431.0);
    const auto fs_1449_2431_715 = std::sqrt(10498005.0 / 41327.0);
    const auto fs_1449_442_221 = std::sqrt(2099601.0 / 884.0);
    const auto fs_1449_442_91 = std::sqrt(14697207.0 / 15028.0);
    const auto fs_1449_4862_10010 = std::sqrt(73486035.0 / 82654.0);
    const auto fs_1449_4862_12155 = std::sqrt(10498005.0 / 9724.0);
    const auto fs_1449_748_22 = std::sqrt(2099601.0 / 25432.0);
    const auto fs_1449_9724_30030 = std::sqrt(220458105.0 / 330616.0);
    const auto fs_1470_46189_35 = std::sqrt(75631500.0 / 2133423721.0);
    const auto fs_147_2431_22 = std::sqrt(43218.0 / 537251.0);
    const auto fs_14_143_11 = std::sqrt(196.0 / 1859.0);
    const auto fs_14_143_3 = std::sqrt(588.0 / 20449.0);
    const auto fs_14_143_5 = std::sqrt(980.0 / 20449.0);
    const auto fs_14_143_7 = std::sqrt(1372.0 / 20449.0);
    const auto fs_14_429_14 = std::sqrt(2744.0 / 184041.0);
    const auto fs_14_429_15 = std::sqrt(980.0 / 61347.0);
    const auto fs_14_429_21 = std::sqrt(1372.0 / 61347.0);
    const auto fs_14_429_3 = std::sqrt(196.0 / 61347.0);
    const auto fs_14_429_33 = std::sqrt(196.0 / 5577.0);
    const auto fs_14_429_35 = std::sqrt(6860.0 / 184041.0);
    const auto fs_14_429_42 = std::sqrt(2744.0 / 61347.0);
    const auto fs_14_429_6 = std::sqrt(392.0 / 61347.0);
    const auto fs_1500_4199_286 = std::sqrt(49500000.0 / 1356277.0);
    const auto fs_1500_4199_33 = std::sqrt(74250000.0 / 17631601.0);
    const auto fs_150_3553_5 = std::sqrt(112500.0 / 12623809.0);
    const auto fs_150_4199_10 = std::sqrt(225000.0 / 17631601.0);
    const auto fs_150_4199_15 = std::sqrt(337500.0 / 17631601.0);
    const auto fs_150_4199_15015 = std::sqrt(25987500.0 / 1356277.0);
    const auto fs_150_4199_1547 = std::sqrt(157500.0 / 79781.0);
    const auto fs_150_4199_2002 = std::sqrt(3465000.0 / 1356277.0);
    const auto fs_150_4199_2431 = std::sqrt(247500.0 / 79781.0);
    const auto fs_150_4199_24310 = std::sqrt(2475000.0 / 79781.0);
    const auto fs_150_4199_41990 = std::sqrt(225000.0 / 4199.0);
    const auto fs_150_4199_6 = std::sqrt(135000.0 / 17631601.0);
    const auto fs_150_4199_7735 = std::sqrt(787500.0 / 79781.0);
    const auto fs_150_46189_231 = std::sqrt(472500.0 / 193947611.0);
    const auto fs_15120_4199_3 = std::sqrt(685843200.0 / 17631601.0);
    const auto fs_15120_46189_35 = std::sqrt(8001504000.0 / 2133423721.0);
    const auto fs_15120_46189_385 = std::sqrt(8001504000.0 / 193947611.0);
    const auto fs_1520_663_2 = std::sqrt(4620800.0 / 439569.0);
    const auto fs_152_13_35 = std::sqrt(808640.0 / 169.0);
    const auto fs_152_221_154 = std::sqrt(3558016.0 / 48841.0);
    const auto fs_152_221_30 = std::sqrt(693120.0 / 48841.0);
    const auto fs_152_221_55 = std::sqrt(1270720.0 / 48841.0);
    const auto fs_152_221_70 = std::sqrt(1617280.0 / 48841.0);
    const auto fs_152_663_14 = std::sqrt(323456.0 / 439569.0);
    const auto fs_152_663_231 = std::sqrt(1779008.0 / 146523.0);
    const auto fs_15435_286_7 = std::sqrt(1667674575.0 / 81796.0);
    const auto fs_15435_572_35 = std::sqrt(8338372875.0 / 327184.0);
    const auto fs_1575_104_210 = std::sqrt(260465625.0 / 5408.0);
    const auto fs_1575_286_1001 = std::sqrt(17364375.0 / 572.0);
    const auto fs_1575_4199_143 = std::sqrt(27286875.0 / 1356277.0);
    const auto fs_1575_44_5 = std::sqrt(12403125.0 / 1936.0);
    const auto fs_1575_52_10 = std::sqrt(12403125.0 / 1352.0);
    const auto fs_1575_52_6 = std::sqrt(7441875.0 / 1352.0);
    const auto fs_1575_572_231 = std::sqrt(52093125.0 / 29744.0);
    const auto fs_15_4199_2730 = std::sqrt(47250.0 / 1356277.0);
    const auto fs_15_46189_10010 = std::sqrt(15750.0 / 14919047.0);
    const auto fs_16065_16_11 = std::sqrt(2838926475.0 / 256.0);
    const auto fs_16065_16_3 = std::sqrt(774252675.0 / 256.0);
    const auto fs_16065_16_5 = std::sqrt(1290421125.0 / 256.0);
    const auto fs_16065_16_7 = std::sqrt(1806589575.0 / 256.0);
    const auto fs_16065_32_10 = std::sqrt(1290421125.0 / 512.0);
    const auto fs_16065_32_22 = std::sqrt(2838926475.0 / 512.0);
    const auto fs_1610_2431_11 = std::sqrt(2592100.0 / 537251.0);
    const auto fs_161_221_78 = std::sqrt(155526.0 / 3757.0);
    const auto fs_161_2431_1430 = std::sqrt(259210.0 / 41327.0);
    const auto fs_1620_46189_30 = std::sqrt(78732000.0 / 2133423721.0);
    const auto fs_1672_663_5 = std::sqrt(13977920.0 / 439569.0);
    const auto fs_1680_46189_42 = std::sqrt(118540800.0 / 2133423721.0);
    const auto fs_168_4199_143 = std::sqrt(310464.0 / 1356277.0);
    const auto fs_16905_4862_11 = std::sqrt(285779025.0 / 2149004.0);
    const auto fs_16_221_5 = std::sqrt(1280.0 / 48841.0);
    const auto fs_16_221_7 = std::sqrt(1792.0 / 48841.0);
    const auto fs_16_2431_1430 = std::sqrt(2560.0 / 41327.0);
    const auto fs_16_663_35 = std::sqrt(8960.0 / 439569.0);
    const auto fs_171_13_66 = std::sqrt(1929906.0 / 169.0);
    const auto fs_177261_9724_2 = std::sqrt(31421462121.0 / 47278088.0);
    const auto fs_178_2431_15 = std::sqrt(475260.0 / 5909761.0);
    const auto fs_180_46189_2002 = std::sqrt(453600.0 / 14919047.0);
    const auto fs_18225_16_10 = std::sqrt(1660753125.0 / 128.0);
    const auto fs_18225_8_2 = std::sqrt(332150625.0 / 32.0);
    const auto fs_18225_8_3 = std::sqrt(996451875.0 / 64.0);
    const auto fs_184_221_91 = std::sqrt(236992.0 / 3757.0);
    const auto fs_1886_2431_66 = std::sqrt(21341976.0 / 537251.0);
    const auto fs_1890_221_65 = std::sqrt(17860500.0 / 3757.0);
    const auto fs_1890_4199_105 = std::sqrt(375070500.0 / 17631601.0);
    const auto fs_1890_46189_10 = std::sqrt(35721000.0 / 2133423721.0);
    const auto fs_18_143_2 = std::sqrt(648.0 / 20449.0);
    const auto fs_18_143_3 = std::sqrt(972.0 / 20449.0);
    const auto fs_18_221_7 = std::sqrt(2268.0 / 48841.0);
    const auto fs_18_2431_286 = std::sqrt(648.0 / 41327.0);
    const auto fs_18_2431_462 = std::sqrt(13608.0 / 537251.0);
    const auto fs_18_4199_12597 = std::sqrt(972.0 / 4199.0);
    const auto fs_18_4199_14586 = std::sqrt(21384.0 / 79781.0);
    const auto fs_18_4199_165 = std::sqrt(53460.0 / 17631601.0);
    const auto fs_18_4199_2431 = std::sqrt(3564.0 / 79781.0);
    const auto fs_18_4199_39 = std::sqrt(972.0 / 1356277.0);
    const auto fs_1900_13_2 = std::sqrt(7220000.0 / 169.0);
    const auto fs_190_13_14 = std::sqrt(505400.0 / 169.0);
    const auto fs_190_13_231 = std::sqrt(8339100.0 / 169.0);
    const auto fs_190_13_42 = std::sqrt(1516200.0 / 169.0);
    const auto fs_1932_2431_143 = std::sqrt(3732624.0 / 41327.0);
    const auto fs_1932_2431_1430 = std::sqrt(37326240.0 / 41327.0);
    const auto fs_195615_92378_6 = std::sqrt(114795684675.0 / 4266847442.0);
    const auto fs_195615_9724_6 = std::sqrt(114795684675.0 / 47278088.0);
    const auto fs_19803_4862_66 = std::sqrt(1176476427.0 / 1074502.0);
    const auto fs_19845_2431_35 = std::sqrt(13783840875.0 / 5909761.0);
    const auto fs_19845_2431_70 = std::sqrt(27567681750.0 / 5909761.0);
    const auto fs_19845_286_3 = std::sqrt(1181472075.0 / 81796.0);
    const auto fs_19845_46189_286 = std::sqrt(787648050.0 / 14919047.0);
    const auto fs_19845_4862_286 = std::sqrt(393824025.0 / 82654.0);
    const auto fs_19845_572_10 = std::sqrt(1969120125.0 / 163592.0);
    const auto fs_19_1_6 = std::sqrt(2166.0);
    const auto fs_1_2431_30030 = std::sqrt(210.0 / 41327.0);
    const auto fs_1_2431_6006 = std::sqrt(42.0 / 41327.0);
    const auto fs_2001_2431_30 = std::sqrt(120120030.0 / 5909761.0);
    const auto fs_2025_22_2 = std::sqrt(4100625.0 / 242.0);
    const auto fs_2025_22_3 = std::sqrt(12301875.0 / 484.0);
    const auto fs_2025_44_10 = std::sqrt(20503125.0 / 968.0);
    const auto fs_207_2431_4290 = std::sqrt(1285470.0 / 41327.0);
    const auto fs_20825_16_2 = std::sqrt(433680625.0 / 128.0);
    const auto fs_20825_286_2 = std::sqrt(433680625.0 / 40898.0);
    const auto fs_20825_44_2 = std::sqrt(433680625.0 / 968.0);
    const auto fs_2090_13_5 = std::sqrt(21840500.0 / 169.0);
    const auto fs_209_13_42 = std::sqrt(1834602.0 / 169.0);
    const auto fs_20_221_11 = std::sqrt(4400.0 / 48841.0);
    const auto fs_20_221_7 = std::sqrt(2800.0 / 48841.0);
    const auto fs_20_2431_231 = std::sqrt(8400.0 / 537251.0);
    const auto fs_20_663_42 = std::sqrt(5600.0 / 146523.0);
    const auto fs_2100_4199_143 = std::sqrt(48510000.0 / 1356277.0);
    const auto fs_2100_46189_11 = std::sqrt(4410000.0 / 193947611.0);
    const auto fs_2100_46189_14 = std::sqrt(61740000.0 / 2133423721.0);
    const auto fs_2100_46189_35 = std::sqrt(154350000.0 / 2133423721.0);
    const auto fs_2100_46189_66 = std::sqrt(26460000.0 / 193947611.0);
    const auto fs_210_46189_4290 = std::sqrt(1323000.0 / 14919047.0);
    const auto fs_210_46189_715 = std::sqrt(220500.0 / 14919047.0);
    const auto fs_214_2431_10 = std::sqrt(457960.0 / 5909761.0);
    const auto fs_216_4199_154 = std::sqrt(7185024.0 / 17631601.0);
    const auto fs_21735_1144_14 = std::sqrt(3306871575.0 / 654368.0);
    const auto fs_21735_9724_10 = std::sqrt(2362051125.0 / 47278088.0);
    const auto fs_21_4199_390 = std::sqrt(13230.0 / 1356277.0);
    const auto fs_2205_104_130 = std::sqrt(24310125.0 / 416.0);
    const auto fs_2205_286_165 = std::sqrt(72930375.0 / 7436.0);
    const auto fs_2205_286_21 = std::sqrt(102102525.0 / 81796.0);
    const auto fs_2205_286_35 = std::sqrt(170170875.0 / 81796.0);
    const auto fs_2205_572_4290 = std::sqrt(72930375.0 / 1144.0);
    const auto fs_2205_572_715 = std::sqrt(24310125.0 / 2288.0);
    const auto fs_2254_2431_33 = std::sqrt(15241548.0 / 537251.0);
    const auto fs_225_1_11 = std::sqrt(556875.0);
    const auto fs_225_1_14 = std::sqrt(708750.0);
    const auto fs_225_1_15 = std::sqrt(759375.0);
    const auto fs_225_1_21 = std::sqrt(1063125.0);
    const auto fs_225_1_3 = std::sqrt(151875.0);
    const auto fs_225_1_33 = std::sqrt(1670625.0);
    const auto fs_225_1_35 = std::sqrt(1771875.0);
    const auto fs_225_1_5 = std::sqrt(253125.0);
    const auto fs_225_2_110 = std::sqrt(2784375.0 / 2.0);
    const auto fs_225_2_2 = std::sqrt(50625.0 / 2.0);
    const auto fs_225_2_6 = std::sqrt(151875.0 / 2.0);
    const auto fs_225_4199_12597 = std::sqrt(151875.0 / 4199.0);
    const auto fs_225_4199_14586 = std::sqrt(3341250.0 / 79781.0);
    const auto fs_225_4199_165 = std::sqrt(8353125.0 / 17631601.0);
    const auto fs_225_4199_2431 = std::sqrt(556875.0 / 79781.0);
    const auto fs_225_4199_30 = std::sqrt(1518750.0 / 17631601.0);
    const auto fs_225_4199_39 = std::sqrt(151875.0 / 1356277.0);
    const auto fs_225_8398_1430 = std::sqrt(2784375.0 / 2712554.0);
    const auto fs_225_8398_30030 = std::sqrt(58471875.0 / 2712554.0);
    const auto fs_2280_13_5 = std::sqrt(25992000.0 / 169.0);
    const auto fs_2280_13_7 = std::sqrt(36388800.0 / 169.0);
    const auto fs_228_221_66 = std::sqrt(3430944.0 / 48841.0);
    const auto fs_22_221_2 = std::sqrt(968.0 / 48841.0);
    const auto fs_22_663_42 = std::sqrt(6776.0 / 146523.0);
    const auto fs_230_221_91 = std::sqrt(370300.0 / 3757.0);
    const auto fs_230_2431_7293 = std::sqrt(158700.0 / 2431.0);
    const auto fs_230_2431_858 = std::sqrt(317400.0 / 41327.0);
    const auto fs_23625_92378_770 = std::sqrt(19534921875.0 / 387895222.0);
    const auto fs_23625_9724_770 = std::sqrt(19534921875.0 / 4298008.0);
    const auto fs_23667_4862_33 = std::sqrt(1680380667.0 / 2149004.0);
    const auto fs_2375_4_3 = std::sqrt(16921875.0 / 16.0);
    const auto fs_2380_429_2 = std::sqrt(11328800.0 / 184041.0);
    const auto fs_23_2431_30030 = std::sqrt(111090.0 / 41327.0);
    const auto fs_23_2431_6006 = std::sqrt(22218.0 / 41327.0);
    const auto fs_2400_46189_33 = std::sqrt(17280000.0 / 193947611.0);
    const auto fs_240_46189_2310 = std::sqrt(12096000.0 / 193947611.0);
    const auto fs_2415_2431_231 = std::sqrt(122476725.0 / 537251.0);
    const auto fs_2415_442_91 = std::sqrt(40825575.0 / 15028.0);
    const auto fs_2415_4862_7293 = std::sqrt(17496675.0 / 9724.0);
    const auto fs_2415_4862_858 = std::sqrt(17496675.0 / 82654.0);
    const auto fs_2432_663_3 = std::sqrt(5914624.0 / 146523.0);
    const auto fs_2438_2431_33 = std::sqrt(17831532.0 / 537251.0);
    const auto fs_24_143_2 = std::sqrt(1152.0 / 20449.0);
    const auto fs_24_2431_70 = std::sqrt(40320.0 / 5909761.0);
    const auto fs_24_4199_2145 = std::sqrt(95040.0 / 1356277.0);
    const auto fs_24_4199_6006 = std::sqrt(266112.0 / 1356277.0);
    const auto fs_24_4199_91 = std::sqrt(4032.0 / 1356277.0);
    const auto fs_2520_46189_30 = std::sqrt(190512000.0 / 2133423721.0);
    const auto fs_252_4199_110 = std::sqrt(6985440.0 / 17631601.0);
    const auto fs_25515_16_10 = std::sqrt(3255076125.0 / 128.0);
    const auto fs_25515_2431_30 = std::sqrt(19530456750.0 / 5909761.0);
    const auto fs_25515_8_2 = std::sqrt(651015225.0 / 32.0);
    const auto fs_25515_8_3 = std::sqrt(1953045675.0 / 64.0);
    const auto fs_25515_92378_462 = std::sqrt(13671319725.0 / 387895222.0);
    const auto fs_25515_9724_462 = std::sqrt(13671319725.0 / 4298008.0);
    const auto fs_25599_4862_33 = std::sqrt(1965926403.0 / 2149004.0);
    const auto fs_2625_8398_286 = std::sqrt(75796875.0 / 2712554.0);
    const auto fs_26460_2431_42 = std::sqrt(29405527200.0 / 5909761.0);
    const auto fs_26775_286_3 = std::sqrt(2150701875.0 / 81796.0);
    const auto fs_2700_4199_154 = std::sqrt(1122660000.0 / 17631601.0);
    const auto fs_270_143_11 = std::sqrt(72900.0 / 1859.0);
    const auto fs_270_143_14 = std::sqrt(1020600.0 / 20449.0);
    const auto fs_270_143_15 = std::sqrt(1093500.0 / 20449.0);
    const auto fs_270_143_21 = std::sqrt(1530900.0 / 20449.0);
    const auto fs_270_143_3 = std::sqrt(218700.0 / 20449.0);
    const auto fs_270_143_33 = std::sqrt(218700.0 / 1859.0);
    const auto fs_270_143_35 = std::sqrt(2551500.0 / 20449.0);
    const auto fs_270_143_5 = std::sqrt(364500.0 / 20449.0);
    const auto fs_270_46189_105 = std::sqrt(7654500.0 / 2133423721.0);
    const auto fs_270_46189_385 = std::sqrt(2551500.0 / 193947611.0);
    const auto fs_270_46189_5 = std::sqrt(364500.0 / 2133423721.0);
    const auto fs_270_46189_770 = std::sqrt(5103000.0 / 193947611.0);
    const auto fs_271215_46189_2 = std::sqrt(147115152450.0 / 2133423721.0);
    const auto fs_271215_4862_2 = std::sqrt(73557576225.0 / 11819522.0);
    const auto fs_27405_92378_462 = std::sqrt(15771714525.0 / 387895222.0);
    const auto fs_27405_92378_70 = std::sqrt(26286190875.0 / 4266847442.0);
    const auto fs_27405_9724_462 = std::sqrt(15771714525.0 / 4298008.0);
    const auto fs_27405_9724_70 = std::sqrt(26286190875.0 / 47278088.0);
    const auto fs_27531_9724_154 = std::sqrt(5305691727.0 / 4298008.0);
    const auto fs_276_143_5 = std::sqrt(380880.0 / 20449.0);
    const auto fs_276_2431_2431 = std::sqrt(76176.0 / 2431.0);
    const auto fs_276_2431_715 = std::sqrt(380880.0 / 41327.0);
    const auto fs_2835_2431_2002 = std::sqrt(112521150.0 / 41327.0);
    const auto fs_2835_3553_70 = std::sqrt(562605750.0 / 12623809.0);
    const auto fs_2835_374_70 = std::sqrt(281302875.0 / 69938.0);
    const auto fs_2835_46189_66 = std::sqrt(48223350.0 / 193947611.0);
    const auto fs_2835_46189_77 = std::sqrt(56260575.0 / 193947611.0);
    const auto fs_2835_4862_66 = std::sqrt(24111675.0 / 1074502.0);
    const auto fs_2835_4862_77 = std::sqrt(56260575.0 / 2149004.0);
    const auto fs_2835_572_105 = std::sqrt(843908625.0 / 327184.0);
    const auto fs_2835_572_385 = std::sqrt(281302875.0 / 29744.0);
    const auto fs_2835_572_5 = std::sqrt(40186125.0 / 327184.0);
    const auto fs_2835_572_770 = std::sqrt(281302875.0 / 14872.0);
    const auto fs_2850_13_11 = std::sqrt(89347500.0 / 169.0);
    const auto fs_2850_13_7 = std::sqrt(56857500.0 / 169.0);
    const auto fs_285_13_10 = std::sqrt(812250.0 / 169.0);
    const auto fs_285_13_210 = std::sqrt(17057250.0 / 169.0);
    const auto fs_285_13_30 = std::sqrt(2436750.0 / 169.0);
    const auto fs_285_1_5 = std::sqrt(406125.0);
    const auto fs_285_1_7 = std::sqrt(568575.0);
    const auto fs_285_4_154 = std::sqrt(6254325.0 / 8.0);
    const auto fs_285_4_30 = std::sqrt(1218375.0 / 8.0);
    const auto fs_285_4_55 = std::sqrt(4467375.0 / 16.0);
    const auto fs_285_4_70 = std::sqrt(2842875.0 / 8.0);
    const auto fs_285_8_10 = std::sqrt(406125.0 / 32.0);
    const auto fs_285_8_210 = std::sqrt(8528625.0 / 32.0);
    const auto fs_2898_2431_70 = std::sqrt(587888280.0 / 5909761.0);
    const auto fs_28_187_3 = std::sqrt(2352.0 / 34969.0);
    const auto fs_28_2431_715 = std::sqrt(3920.0 / 41327.0);
    const auto fs_28_429_2 = std::sqrt(1568.0 / 184041.0);
    const auto fs_28_429_21 = std::sqrt(5488.0 / 61347.0);
    const auto fs_28_429_7 = std::sqrt(5488.0 / 184041.0);
    const auto fs_2940_46189_7 = std::sqrt(60505200.0 / 2133423721.0);
    const auto fs_2975_11_2 = std::sqrt(17701250.0 / 121.0);
    const auto fs_2975_11_21 = std::sqrt(185863125.0 / 121.0);
    const auto fs_2975_11_7 = std::sqrt(61954375.0 / 121.0);
    const auto fs_2975_143_14 = std::sqrt(123908750.0 / 20449.0);
    const auto fs_2975_143_15 = std::sqrt(132759375.0 / 20449.0);
    const auto fs_2975_143_21 = std::sqrt(185863125.0 / 20449.0);
    const auto fs_2975_143_3 = std::sqrt(26551875.0 / 20449.0);
    const auto fs_2975_143_33 = std::sqrt(26551875.0 / 1859.0);
    const auto fs_2975_143_35 = std::sqrt(309771875.0 / 20449.0);
    const auto fs_2975_143_42 = std::sqrt(371726250.0 / 20449.0);
    const auto fs_2975_143_6 = std::sqrt(53103750.0 / 20449.0);
    const auto fs_2975_16_30 = std::sqrt(132759375.0 / 128.0);
    const auto fs_2975_16_6 = std::sqrt(26551875.0 / 128.0);
    const auto fs_2975_16_66 = std::sqrt(292070625.0 / 128.0);
    const auto fs_2975_22_14 = std::sqrt(61954375.0 / 242.0);
    const auto fs_2975_22_15 = std::sqrt(132759375.0 / 484.0);
    const auto fs_2975_22_21 = std::sqrt(185863125.0 / 484.0);
    const auto fs_2975_22_3 = std::sqrt(26551875.0 / 484.0);
    const auto fs_2975_22_33 = std::sqrt(26551875.0 / 44.0);
    const auto fs_2975_22_35 = std::sqrt(309771875.0 / 484.0);
    const auto fs_2975_22_42 = std::sqrt(185863125.0 / 242.0);
    const auto fs_2975_22_6 = std::sqrt(26551875.0 / 242.0);
    const auto fs_2975_286_30 = std::sqrt(132759375.0 / 40898.0);
    const auto fs_2975_286_6 = std::sqrt(26551875.0 / 40898.0);
    const auto fs_2975_286_66 = std::sqrt(26551875.0 / 3718.0);
    const auto fs_2975_2_2 = std::sqrt(8850625.0 / 2.0);
    const auto fs_2975_44_30 = std::sqrt(132759375.0 / 968.0);
    const auto fs_2975_44_6 = std::sqrt(26551875.0 / 968.0);
    const auto fs_2975_44_66 = std::sqrt(26551875.0 / 88.0);
    const auto fs_2975_4_2 = std::sqrt(8850625.0 / 8.0);
    const auto fs_2975_4_21 = std::sqrt(185863125.0 / 16.0);
    const auto fs_2975_4_7 = std::sqrt(61954375.0 / 16.0);
    const auto fs_2975_8_14 = std::sqrt(61954375.0 / 32.0);
    const auto fs_2975_8_15 = std::sqrt(132759375.0 / 64.0);
    const auto fs_2975_8_21 = std::sqrt(185863125.0 / 64.0);
    const auto fs_2975_8_3 = std::sqrt(26551875.0 / 64.0);
    const auto fs_2975_8_33 = std::sqrt(292070625.0 / 64.0);
    const auto fs_2975_8_35 = std::sqrt(309771875.0 / 64.0);
    const auto fs_2975_8_42 = std::sqrt(185863125.0 / 32.0);
    const auto fs_2975_8_6 = std::sqrt(26551875.0 / 32.0);
    const auto fs_2_187_385 = std::sqrt(140.0 / 3179.0);
    const auto fs_2_221_10 = std::sqrt(40.0 / 48841.0);
    const auto fs_2_221_1326 = std::sqrt(24.0 / 221.0);
    const auto fs_2_221_182 = std::sqrt(56.0 / 3757.0);
    const auto fs_2_221_210 = std::sqrt(840.0 / 48841.0);
    const auto fs_2_221_221 = std::sqrt(4.0 / 221.0);
    const auto fs_2_221_455 = std::sqrt(140.0 / 3757.0);
    const auto fs_2_2431_165 = std::sqrt(60.0 / 537251.0);
    const auto fs_2_2431_2145 = std::sqrt(60.0 / 41327.0);
    const auto fs_2_51_6 = std::sqrt(8.0 / 867.0);
    const auto fs_300_2717_3 = std::sqrt(270000.0 / 7382089.0);
    const auto fs_300_4199_2145 = std::sqrt(14850000.0 / 1356277.0);
    const auto fs_300_4199_6006 = std::sqrt(41580000.0 / 1356277.0);
    const auto fs_300_4199_91 = std::sqrt(630000.0 / 1356277.0);
    const auto fs_300_46189_1001 = std::sqrt(630000.0 / 14919047.0);
    const auto fs_3040_13_3 = std::sqrt(27724800.0 / 169.0);
    const auto fs_304_663_210 = std::sqrt(6469120.0 / 146523.0);
    const auto fs_304_663_462 = std::sqrt(14232064.0 / 146523.0);
    const auto fs_30_2431_1001 = std::sqrt(6300.0 / 41327.0);
    const auto fs_30_2431_11 = std::sqrt(900.0 / 537251.0);
    const auto fs_30_2431_154 = std::sqrt(12600.0 / 537251.0);
    const auto fs_30_2431_66 = std::sqrt(5400.0 / 537251.0);
    const auto fs_30_2431_715 = std::sqrt(4500.0 / 41327.0);
    const auto fs_30_2431_858 = std::sqrt(5400.0 / 41327.0);
    const auto fs_30_3553_1001 = std::sqrt(81900.0 / 1147619.0);
    const auto fs_30_4199_2431 = std::sqrt(9900.0 / 79781.0);
    const auto fs_30_4199_4862 = std::sqrt(19800.0 / 79781.0);
    const auto fs_30_4199_770 = std::sqrt(693000.0 / 17631601.0);
    const auto fs_3105_46189_6 = std::sqrt(57846150.0 / 2133423721.0);
    const auto fs_3150_4199_110 = std::sqrt(1091475000.0 / 17631601.0);
    const auto fs_315_104_2730 = std::sqrt(10418625.0 / 416.0);
    const auto fs_315_1144_10010 = std::sqrt(3472875.0 / 4576.0);
    const auto fs_315_13_65 = std::sqrt(496125.0 / 13.0);
    const auto fs_315_26_105 = std::sqrt(10418625.0 / 676.0);
    const auto fs_315_44_1001 = std::sqrt(9029475.0 / 176.0);
    const auto fs_315_52_770 = std::sqrt(38201625.0 / 1352.0);
    const auto fs_33075_2431_11 = std::sqrt(1093955625.0 / 537251.0);
    const auto fs_33075_2431_14 = std::sqrt(15315378750.0 / 5909761.0);
    const auto fs_33075_2431_35 = std::sqrt(38288446875.0 / 5909761.0);
    const auto fs_33075_2431_66 = std::sqrt(6563733750.0 / 537251.0);
    const auto fs_33075_46189_154 = std::sqrt(15315378750.0 / 193947611.0);
    const auto fs_33075_4862_154 = std::sqrt(7657689375.0 / 1074502.0);
    const auto fs_330_4199_3 = std::sqrt(326700.0 / 17631601.0);
    const auto fs_3330_46189_10 = std::sqrt(110889000.0 / 2133423721.0);
    const auto fs_3375_8398_22 = std::sqrt(125296875.0 / 35263202.0);
    const auto fs_3381_187_3 = std::sqrt(34293483.0 / 34969.0);
    const auto fs_3381_2431_22 = std::sqrt(22862322.0 / 537251.0);
    const auto fs_3381_2431_715 = std::sqrt(57155805.0 / 41327.0);
    const auto fs_3381_884_78 = std::sqrt(34293483.0 / 30056.0);
    const auto fs_3381_9724_1430 = std::sqrt(57155805.0 / 330616.0);
    const auto fs_346815_92378_2 = std::sqrt(120280644225.0 / 4266847442.0);
    const auto fs_346815_9724_2 = std::sqrt(120280644225.0 / 47278088.0);
    const auto fs_34965_46189_55 = std::sqrt(6112756125.0 / 193947611.0);
    const auto fs_34965_4862_55 = std::sqrt(6112756125.0 / 2149004.0);
    const auto fs_34965_572_10 = std::sqrt(6112756125.0 / 163592.0);
    const auto fs_367_2431_2 = std::sqrt(269378.0 / 5909761.0);
    const auto fs_368_2431_1430 = std::sqrt(1354240.0 / 41327.0);
    const auto fs_36_2431_35 = std::sqrt(45360.0 / 5909761.0);
    const auto fs_36_4199_143 = std::sqrt(14256.0 / 1356277.0);
    const auto fs_36_4199_3003 = std::sqrt(299376.0 / 1356277.0);
    const auto fs_37485_286_5 = std::sqrt(7025626125.0 / 81796.0);
    const auto fs_37485_32_2 = std::sqrt(1405125225.0 / 512.0);
    const auto fs_375_4199_2431 = std::sqrt(1546875.0 / 79781.0);
    const auto fs_375_4199_4862 = std::sqrt(3093750.0 / 79781.0);
    const auto fs_375_46189_770 = std::sqrt(9843750.0 / 193947611.0);
    const auto fs_37800_2431_33 = std::sqrt(4286520000.0 / 537251.0);
    const auto fs_3780_2431_2310 = std::sqrt(3000564000.0 / 537251.0);
    const auto fs_3780_4199_65 = std::sqrt(71442000.0 / 1356277.0);
    const auto fs_3780_46189_3 = std::sqrt(42865200.0 / 2133423721.0);
    const auto fs_3800_663_3 = std::sqrt(14440000.0 / 146523.0);
    const auto fs_380_13_2 = std::sqrt(288800.0 / 169.0);
    const auto fs_380_13_210 = std::sqrt(30324000.0 / 169.0);
    const auto fs_380_13_462 = std::sqrt(66712800.0 / 169.0);
    const auto fs_380_1_3 = std::sqrt(433200.0);
    const auto fs_380_221_30 = std::sqrt(4332000.0 / 48841.0);
    const auto fs_38_13_14 = std::sqrt(20216.0 / 169.0);
    const auto fs_38_13_231 = std::sqrt(333564.0 / 169.0);
    const auto fs_39690_2431_30 = std::sqrt(47258883000.0 / 5909761.0);
    const auto fs_39690_46189_35 = std::sqrt(55135363500.0 / 2133423721.0);
    const auto fs_39690_46189_70 = std::sqrt(110270727000.0 / 2133423721.0);
    const auto fs_396_4199_35 = std::sqrt(5488560.0 / 17631601.0);
    const auto fs_3_143_110 = std::sqrt(90.0 / 1859.0);
    const auto fs_3_143_2 = std::sqrt(18.0 / 20449.0);
    const auto fs_3_143_6 = std::sqrt(54.0 / 20449.0);
    const auto fs_3_187_22 = std::sqrt(18.0 / 3179.0);
    const auto fs_3_2431_30030 = std::sqrt(1890.0 / 41327.0);
    const auto fs_3_4199_182 = std::sqrt(126.0 / 1356277.0);
    const auto fs_3_4199_2 = std::sqrt(18.0 / 17631601.0);
    const auto fs_3_4199_2002 = std::sqrt(1386.0 / 1356277.0);
    const auto fs_3_4199_26 = std::sqrt(18.0 / 1356277.0);
    const auto fs_3_4199_6006 = std::sqrt(4158.0 / 1356277.0);
    const auto fs_3_4199_910 = std::sqrt(630.0 / 1356277.0);
    const auto fs_405_143_10 = std::sqrt(1640250.0 / 20449.0);
    const auto fs_405_46189_462 = std::sqrt(6889050.0 / 193947611.0);
    const auto fs_4094_2431_15 = std::sqrt(251412540.0 / 5909761.0);
    const auto fs_40_663_2 = std::sqrt(3200.0 / 439569.0);
    const auto fs_4125_4199_3 = std::sqrt(51046875.0 / 17631601.0);
    const auto fs_414_221_7 = std::sqrt(1199772.0 / 48841.0);
    const auto fs_414_2431_286 = std::sqrt(342792.0 / 41327.0);
    const auto fs_414_2431_462 = std::sqrt(7198632.0 / 537251.0);
    const auto fs_4165_858_2 = std::sqrt(17347225.0 / 368082.0);
    const auto fs_418_13_5 = std::sqrt(873620.0 / 169.0);
    const auto fs_42021_9724_30 = std::sqrt(26486466615.0 / 47278088.0);
    const auto fs_420_2717_5 = std::sqrt(882000.0 / 7382089.0);
    const auto fs_420_46189_165 = std::sqrt(2646000.0 / 193947611.0);
    const auto fs_420_46189_21 = std::sqrt(3704400.0 / 2133423721.0);
    const auto fs_420_46189_35 = std::sqrt(6174000.0 / 2133423721.0);
    const auto fs_42525_32_2 = std::sqrt(1808375625.0 / 512.0);
    const auto fs_42525_32_3 = std::sqrt(5425126875.0 / 1024.0);
    const auto fs_42525_64_10 = std::sqrt(9041878125.0 / 2048.0);
    const auto fs_42987_4862_15 = std::sqrt(27718232535.0 / 23639044.0);
    const auto fs_42_2431_165 = std::sqrt(26460.0 / 537251.0);
    const auto fs_42_4199_2145 = std::sqrt(291060.0 / 1356277.0);
    const auto fs_42_4199_286 = std::sqrt(38808.0 / 1356277.0);
    const auto fs_42_4199_715 = std::sqrt(97020.0 / 1356277.0);
    const auto fs_4347_2431_35 = std::sqrt(661374315.0 / 5909761.0);
    const auto fs_4347_442_7 = std::sqrt(132274863.0 / 195364.0);
    const auto fs_4347_4862_286 = std::sqrt(18896409.0 / 82654.0);
    const auto fs_4347_4862_462 = std::sqrt(396824589.0 / 1074502.0);
    const auto fs_4347_9724_4290 = std::sqrt(283446135.0 / 330616.0);
    const auto fs_435_46189_462 = std::sqrt(7947450.0 / 193947611.0);
    const auto fs_435_46189_70 = std::sqrt(13245750.0 / 2133423721.0);
    const auto fs_4410_143_42 = std::sqrt(816820200.0 / 20449.0);
    const auto fs_44_663_5 = std::sqrt(9680.0 / 439569.0);
    const auto fs_450_1_5 = std::sqrt(1012500.0);
    const auto fs_450_4199_13 = std::sqrt(202500.0 / 1356277.0);
    const auto fs_450_4199_143 = std::sqrt(2227500.0 / 1356277.0);
    const auto fs_450_4199_3003 = std::sqrt(46777500.0 / 1356277.0);
    const auto fs_450_46189_385 = std::sqrt(7087500.0 / 193947611.0);
    const auto fs_456_13_5 = std::sqrt(1039680.0 / 169.0);
    const auto fs_456_13_7 = std::sqrt(1455552.0 / 169.0);
    const auto fs_45_2431_10 = std::sqrt(20250.0 / 5909761.0);
    const auto fs_45_4199_2002 = std::sqrt(311850.0 / 1356277.0);
    const auto fs_460_2431_231 = std::sqrt(4443600.0 / 537251.0);
    const auto fs_46305_2431_7 = std::sqrt(15009071175.0 / 5909761.0);
    const auto fs_46305_46189_35 = std::sqrt(75045355875.0 / 2133423721.0);
    const auto fs_46305_4862_35 = std::sqrt(75045355875.0 / 23639044.0);
    const auto fs_46_187_385 = std::sqrt(74060.0 / 3179.0);
    const auto fs_46_221_1326 = std::sqrt(12696.0 / 221.0);
    const auto fs_46_221_182 = std::sqrt(29624.0 / 3757.0);
    const auto fs_46_221_221 = std::sqrt(2116.0 / 221.0);
    const auto fs_46_221_455 = std::sqrt(74060.0 / 3757.0);
    const auto fs_46_2431_165 = std::sqrt(31740.0 / 537251.0);
    const auto fs_46_2431_2145 = std::sqrt(31740.0 / 41327.0);
    const auto fs_4725_104_30 = std::sqrt(334884375.0 / 5408.0);
    const auto fs_4725_143_3 = std::sqrt(66976875.0 / 20449.0);
    const auto fs_4725_2431_1001 = std::sqrt(156279375.0 / 41327.0);
    const auto fs_4725_3553_5 = std::sqrt(111628125.0 / 12623809.0);
    const auto fs_4725_374_5 = std::sqrt(111628125.0 / 139876.0);
    const auto fs_4725_4199_10 = std::sqrt(223256250.0 / 17631601.0);
    const auto fs_4725_4199_6 = std::sqrt(133953750.0 / 17631601.0);
    const auto fs_4725_442_10 = std::sqrt(111628125.0 / 97682.0);
    const auto fs_4725_442_6 = std::sqrt(66976875.0 / 97682.0);
    const auto fs_4725_46189_231 = std::sqrt(468838125.0 / 193947611.0);
    const auto fs_4725_4862_231 = std::sqrt(468838125.0 / 2149004.0);
    const auto fs_4725_52_13 = std::sqrt(22325625.0 / 208.0);
    const auto fs_4725_572_385 = std::sqrt(781396875.0 / 29744.0);
    const auto fs_4725_8398_210 = std::sqrt(2344190625.0 / 35263202.0);
    const auto fs_4725_884_210 = std::sqrt(2344190625.0 / 390728.0);
    const auto fs_4750_13_3 = std::sqrt(67687500.0 / 169.0);
    const auto fs_475_2_2 = std::sqrt(225625.0 / 2.0);
    const auto fs_475_4_42 = std::sqrt(4738125.0 / 8.0);
    const auto fs_480_4199_3 = std::sqrt(691200.0 / 17631601.0);
    const auto fs_480_46189_35 = std::sqrt(8064000.0 / 2133423721.0);
    const auto fs_480_46189_385 = std::sqrt(8064000.0 / 193947611.0);
    const auto fs_483_143_105 = std::sqrt(24495345.0 / 20449.0);
    const auto fs_483_221_14 = std::sqrt(3266046.0 / 48841.0);
    const auto fs_483_221_3 = std::sqrt(699867.0 / 48841.0);
    const auto fs_483_221_39 = std::sqrt(699867.0 / 3757.0);
    const auto fs_483_221_442 = std::sqrt(466578.0 / 221.0);
    const auto fs_483_221_455 = std::sqrt(8165115.0 / 3757.0);
    const auto fs_483_2431_462 = std::sqrt(9798138.0 / 537251.0);
    const auto fs_483_374_385 = std::sqrt(8165115.0 / 12716.0);
    const auto fs_483_442_1326 = std::sqrt(699867.0 / 442.0);
    const auto fs_483_442_182 = std::sqrt(1633023.0 / 7514.0);
    const auto fs_483_442_221 = std::sqrt(233289.0 / 884.0);
    const auto fs_483_442_455 = std::sqrt(8165115.0 / 15028.0);
    const auto fs_483_4862_165 = std::sqrt(3499335.0 / 2149004.0);
    const auto fs_483_4862_2145 = std::sqrt(3499335.0 / 165308.0);
    const auto fs_483_9724_30030 = std::sqrt(24495345.0 / 330616.0);
    const auto fs_483_9724_6006 = std::sqrt(4899069.0 / 330616.0);
    const auto fs_48_2431_154 = std::sqrt(32256.0 / 537251.0);
    const auto fs_48_2431_165 = std::sqrt(34560.0 / 537251.0);
    const auto fs_4922_2431_10 = std::sqrt(242260840.0 / 5909761.0);
    const auto fs_4950_4199_35 = std::sqrt(857587500.0 / 17631601.0);
    const auto fs_49_429_2 = std::sqrt(4802.0 / 184041.0);
    const auto fs_4_143_105 = std::sqrt(1680.0 / 20449.0);
    const auto fs_4_221_14 = std::sqrt(224.0 / 48841.0);
    const auto fs_4_221_154 = std::sqrt(2464.0 / 48841.0);
    const auto fs_4_221_3 = std::sqrt(48.0 / 48841.0);
    const auto fs_4_221_30 = std::sqrt(480.0 / 48841.0);
    const auto fs_4_221_39 = std::sqrt(48.0 / 3757.0);
    const auto fs_4_221_442 = std::sqrt(32.0 / 221.0);
    const auto fs_4_221_455 = std::sqrt(560.0 / 3757.0);
    const auto fs_4_221_55 = std::sqrt(880.0 / 48841.0);
    const auto fs_4_221_70 = std::sqrt(1120.0 / 48841.0);
    const auto fs_4_2431_462 = std::sqrt(672.0 / 537251.0);
    const auto fs_4_663_14 = std::sqrt(224.0 / 439569.0);
    const auto fs_4_663_231 = std::sqrt(1232.0 / 146523.0);
    const auto fs_506_221_2 = std::sqrt(512072.0 / 48841.0);
    const auto fs_50_2431_143 = std::sqrt(2500.0 / 41327.0);
    const auto fs_51030_46189_30 = std::sqrt(78121827000.0 / 2133423721.0);
    const auto fs_51681_4862_10 = std::sqrt(13354628805.0 / 11819522.0);
    const auto fs_525_4199_2145 = std::sqrt(45478125.0 / 1356277.0);
    const auto fs_525_4199_286 = std::sqrt(6063750.0 / 1356277.0);
    const auto fs_525_4199_715 = std::sqrt(15159375.0 / 1356277.0);
    const auto fs_525_8398_390 = std::sqrt(4134375.0 / 2712554.0);
    const auto fs_52920_46189_42 = std::sqrt(117622108800.0 / 2133423721.0);
    const auto fs_530_2431_3 = std::sqrt(842700.0 / 5909761.0);
    const auto fs_5313_442_2 = std::sqrt(28227969.0 / 97682.0);
    const auto fs_5355_16_14 = std::sqrt(200732175.0 / 128.0);
    const auto fs_5355_16_15 = std::sqrt(430140375.0 / 256.0);
    const auto fs_5355_16_21 = std::sqrt(602196525.0 / 256.0);
    const auto fs_5355_16_3 = std::sqrt(86028075.0 / 256.0);
    const auto fs_5355_16_33 = std::sqrt(946308825.0 / 256.0);
    const auto fs_5355_16_35 = std::sqrt(1003660875.0 / 256.0);
    const auto fs_5355_16_42 = std::sqrt(602196525.0 / 128.0);
    const auto fs_5355_16_6 = std::sqrt(86028075.0 / 128.0);
    const auto fs_5355_32_30 = std::sqrt(430140375.0 / 512.0);
    const auto fs_5355_32_6 = std::sqrt(86028075.0 / 512.0);
    const auto fs_5355_32_66 = std::sqrt(946308825.0 / 512.0);
    const auto fs_5355_4_2 = std::sqrt(28676025.0 / 8.0);
    const auto fs_5355_8_2 = std::sqrt(28676025.0 / 32.0);
    const auto fs_5355_8_21 = std::sqrt(602196525.0 / 64.0);
    const auto fs_5355_8_7 = std::sqrt(200732175.0 / 64.0);
    const auto fs_540_143_5 = std::sqrt(1458000.0 / 20449.0);
    const auto fs_540_2431_2 = std::sqrt(583200.0 / 5909761.0);
    const auto fs_54_4199_2002 = std::sqrt(449064.0 / 1356277.0);
    const auto fs_54_4199_286 = std::sqrt(64152.0 / 1356277.0);
    const auto fs_5505_46189_2 = std::sqrt(60610050.0 / 2133423721.0);
    const auto fs_552_2431_70 = std::sqrt(21329280.0 / 5909761.0);
    const auto fs_5670_46189_2002 = std::sqrt(450084600.0 / 14919047.0);
    const auto fs_56_429_2 = std::sqrt(6272.0 / 184041.0);
    const auto fs_570_13_11 = std::sqrt(3573900.0 / 169.0);
    const auto fs_570_13_154 = std::sqrt(50034600.0 / 169.0);
    const auto fs_570_13_30 = std::sqrt(9747000.0 / 169.0);
    const auto fs_570_13_55 = std::sqrt(17869500.0 / 169.0);
    const auto fs_570_13_7 = std::sqrt(2274300.0 / 169.0);
    const auto fs_570_13_70 = std::sqrt(22743000.0 / 169.0);
    const auto fs_5796_2431_154 = std::sqrt(470310624.0 / 537251.0);
    const auto fs_5796_2431_165 = std::sqrt(503904240.0 / 537251.0);
    const auto fs_57_13_10 = std::sqrt(32490.0 / 169.0);
    const auto fs_57_13_210 = std::sqrt(682290.0 / 169.0);
    const auto fs_57_2431_154 = std::sqrt(45486.0 / 537251.0);
    const auto fs_58_2431_143 = std::sqrt(3364.0 / 41327.0);
    const auto fs_5950_11_2 = std::sqrt(70805000.0 / 121.0);
    const auto fs_5950_143_2 = std::sqrt(70805000.0 / 20449.0);
    const auto fs_5950_143_21 = std::sqrt(743452500.0 / 20449.0);
    const auto fs_5950_143_7 = std::sqrt(247817500.0 / 20449.0);
    const auto fs_59535_2431_3 = std::sqrt(10633248675.0 / 5909761.0);
    const auto fs_59535_46189_10 = std::sqrt(35444162250.0 / 2133423721.0);
    const auto fs_59535_4862_10 = std::sqrt(17722081125.0 / 11819522.0);
    const auto fs_595_143_11 = std::sqrt(354025.0 / 1859.0);
    const auto fs_595_143_3 = std::sqrt(1062075.0 / 20449.0);
    const auto fs_595_143_5 = std::sqrt(1770125.0 / 20449.0);
    const auto fs_595_143_7 = std::sqrt(2478175.0 / 20449.0);
    const auto fs_595_286_10 = std::sqrt(1770125.0 / 40898.0);
    const auto fs_595_286_22 = std::sqrt(354025.0 / 3718.0);
    const auto fs_595_429_14 = std::sqrt(4956350.0 / 184041.0);
    const auto fs_595_429_15 = std::sqrt(1770125.0 / 61347.0);
    const auto fs_595_429_21 = std::sqrt(2478175.0 / 61347.0);
    const auto fs_595_429_3 = std::sqrt(354025.0 / 61347.0);
    const auto fs_595_429_33 = std::sqrt(354025.0 / 5577.0);
    const auto fs_595_429_35 = std::sqrt(12390875.0 / 184041.0);
    const auto fs_595_429_42 = std::sqrt(4956350.0 / 61347.0);
    const auto fs_595_429_6 = std::sqrt(708050.0 / 61347.0);
    const auto fs_595_858_30 = std::sqrt(1770125.0 / 122694.0);
    const auto fs_595_858_6 = std::sqrt(354025.0 / 122694.0);
    const auto fs_595_858_66 = std::sqrt(354025.0 / 11154.0);
    const auto fs_6075_16_110 = std::sqrt(2029809375.0 / 128.0);
    const auto fs_6075_16_2 = std::sqrt(36905625.0 / 128.0);
    const auto fs_6075_16_6 = std::sqrt(110716875.0 / 128.0);
    const auto fs_6075_2_2 = std::sqrt(36905625.0 / 2.0);
    const auto fs_6075_4_5 = std::sqrt(184528125.0 / 16.0);
    const auto fs_6075_8_11 = std::sqrt(405961875.0 / 64.0);
    const auto fs_6075_8_14 = std::sqrt(258339375.0 / 32.0);
    const auto fs_6075_8_15 = std::sqrt(553584375.0 / 64.0);
    const auto fs_6075_8_21 = std::sqrt(775018125.0 / 64.0);
    const auto fs_6075_8_3 = std::sqrt(110716875.0 / 64.0);
    const auto fs_6075_8_33 = std::sqrt(1217885625.0 / 64.0);
    const auto fs_6075_8_35 = std::sqrt(1291696875.0 / 64.0);
    const auto fs_6075_8_5 = std::sqrt(184528125.0 / 64.0);
    const auto fs_608_13_3 = std::sqrt(1108992.0 / 169.0);
    const auto fs_608_221_5 = std::sqrt(1848320.0 / 48841.0);
    const auto fs_608_221_7 = std::sqrt(2587648.0 / 48841.0);
    const auto fs_608_663_35 = std::sqrt(12938240.0 / 439569.0);
    const auto fs_60_4199_105 = std::sqrt(378000.0 / 17631601.0);
    const auto fs_60_4199_11 = std::sqrt(39600.0 / 17631601.0);
    const auto fs_6300_143_33 = std::sqrt(119070000.0 / 1859.0);
    const auto fs_630_143_2310 = std::sqrt(83349000.0 / 1859.0);
    const auto fs_630_46189_286 = std::sqrt(793800.0 / 14919047.0);
    const auto fs_644_187_3 = std::sqrt(1244208.0 / 34969.0);
    const auto fs_644_2431_715 = std::sqrt(2073680.0 / 41327.0);
    const auto fs_64_663_3 = std::sqrt(4096.0 / 146523.0);
    const auto fs_65205_1144_6 = std::sqrt(12755076075.0 / 654368.0);
    const auto fs_65205_2431_2 = std::sqrt(8503384050.0 / 5909761.0);
    const auto fs_65205_92378_14 = std::sqrt(29761844175.0 / 4266847442.0);
    const auto fs_65205_9724_14 = std::sqrt(29761844175.0 / 47278088.0);
    const auto fs_66150_46189_11 = std::sqrt(4375822500.0 / 193947611.0);
    const auto fs_66150_46189_14 = std::sqrt(61261515000.0 / 2133423721.0);
    const auto fs_66150_46189_35 = std::sqrt(153153787500.0 / 2133423721.0);
    const auto fs_66150_46189_66 = std::sqrt(26254935000.0 / 193947611.0);
    const auto fs_6615_143_30 = std::sqrt(1312746750.0 / 20449.0);
    const auto fs_6615_143_5 = std::sqrt(218791125.0 / 20449.0);
    const auto fs_6615_2431_165 = std::sqrt(656373375.0 / 537251.0);
    const auto fs_6615_2431_21 = std::sqrt(918922725.0 / 5909761.0);
    const auto fs_6615_2431_35 = std::sqrt(1531537875.0 / 5909761.0);
    const auto fs_6615_286_35 = std::sqrt(1531537875.0 / 81796.0);
    const auto fs_6615_286_70 = std::sqrt(1531537875.0 / 40898.0);
    const auto fs_6615_46189_4290 = std::sqrt(1312746750.0 / 14919047.0);
    const auto fs_6615_46189_715 = std::sqrt(218791125.0 / 14919047.0);
    const auto fs_6615_4862_4290 = std::sqrt(656373375.0 / 82654.0);
    const auto fs_6615_4862_715 = std::sqrt(218791125.0 / 165308.0);
    const auto fs_6615_572_286 = std::sqrt(43758225.0 / 1144.0);
    const auto fs_6615_8398_130 = std::sqrt(218791125.0 / 2712554.0);
    const auto fs_6615_884_130 = std::sqrt(218791125.0 / 30056.0);
    const auto fs_66_4199_5 = std::sqrt(21780.0 / 17631601.0);
    const auto fs_675_11_5 = std::sqrt(2278125.0 / 121.0);
    const auto fs_675_1_2 = std::sqrt(911250.0);
    const auto fs_675_1_3 = std::sqrt(1366875.0);
    const auto fs_675_22_11 = std::sqrt(455625.0 / 44.0);
    const auto fs_675_22_14 = std::sqrt(3189375.0 / 242.0);
    const auto fs_675_22_15 = std::sqrt(6834375.0 / 484.0);
    const auto fs_675_22_21 = std::sqrt(9568125.0 / 484.0);
    const auto fs_675_22_3 = std::sqrt(1366875.0 / 484.0);
    const auto fs_675_22_33 = std::sqrt(1366875.0 / 44.0);
    const auto fs_675_22_35 = std::sqrt(15946875.0 / 484.0);
    const auto fs_675_22_5 = std::sqrt(2278125.0 / 484.0);
    const auto fs_675_2_10 = std::sqrt(2278125.0 / 2.0);
    const auto fs_675_4199_2002 = std::sqrt(70166250.0 / 1356277.0);
    const auto fs_675_4199_286 = std::sqrt(10023750.0 / 1356277.0);
    const auto fs_675_44_110 = std::sqrt(2278125.0 / 88.0);
    const auto fs_675_44_2 = std::sqrt(455625.0 / 968.0);
    const auto fs_675_44_6 = std::sqrt(1366875.0 / 968.0);
    const auto fs_690_2431_1001 = std::sqrt(3332700.0 / 41327.0);
    const auto fs_690_2431_11 = std::sqrt(476100.0 / 537251.0);
    const auto fs_690_2431_154 = std::sqrt(6665400.0 / 537251.0);
    const auto fs_690_2431_66 = std::sqrt(2856600.0 / 537251.0);
    const auto fs_690_2431_715 = std::sqrt(2380500.0 / 41327.0);
    const auto fs_690_2431_858 = std::sqrt(2856600.0 / 41327.0);
    const auto fs_69_187_22 = std::sqrt(9522.0 / 3179.0);
    const auto fs_69_2431_30030 = std::sqrt(999810.0 / 41327.0);
    const auto fs_6_143_11 = std::sqrt(36.0 / 1859.0);
    const auto fs_6_143_14 = std::sqrt(504.0 / 20449.0);
    const auto fs_6_143_15 = std::sqrt(540.0 / 20449.0);
    const auto fs_6_143_21 = std::sqrt(756.0 / 20449.0);
    const auto fs_6_143_3 = std::sqrt(108.0 / 20449.0);
    const auto fs_6_143_33 = std::sqrt(108.0 / 1859.0);
    const auto fs_6_143_35 = std::sqrt(1260.0 / 20449.0);
    const auto fs_6_143_5 = std::sqrt(180.0 / 20449.0);
    const auto fs_6_221_221 = std::sqrt(36.0 / 221.0);
    const auto fs_6_221_66 = std::sqrt(2376.0 / 48841.0);
    const auto fs_6_221_91 = std::sqrt(252.0 / 3757.0);
    const auto fs_6_2431_10010 = std::sqrt(2520.0 / 41327.0);
    const auto fs_6_2431_12155 = std::sqrt(180.0 / 2431.0);
    const auto fs_6_4199_11 = std::sqrt(396.0 / 17631601.0);
    const auto fs_6_4199_138567 = std::sqrt(1188.0 / 4199.0);
    const auto fs_6_4199_143 = std::sqrt(396.0 / 1356277.0);
    const auto fs_6_4199_146965 = std::sqrt(1260.0 / 4199.0);
    const auto fs_6_4199_176358 = std::sqrt(1512.0 / 4199.0);
    const auto fs_6_4199_25194 = std::sqrt(216.0 / 4199.0);
    const auto fs_6_4199_30030 = std::sqrt(83160.0 / 1356277.0);
    const auto fs_6_4199_323323 = std::sqrt(2772.0 / 4199.0);
    const auto fs_6_4199_33 = std::sqrt(1188.0 / 17631601.0);
    const auto fs_6_4199_46189 = std::sqrt(396.0 / 4199.0);
    const auto fs_6_4199_62985 = std::sqrt(540.0 / 4199.0);
    const auto fs_6_4199_910 = std::sqrt(2520.0 / 1356277.0);
    const auto fs_6_4199_92378 = std::sqrt(792.0 / 4199.0);
    const auto fs_6_4199_9282 = std::sqrt(1512.0 / 79781.0);
    const auto fs_70_2431_11 = std::sqrt(4900.0 / 537251.0);
    const auto fs_71001_9724_22 = std::sqrt(5041142001.0 / 4298008.0);
    const auto fs_7245_4862_1001 = std::sqrt(367430175.0 / 165308.0);
    const auto fs_7245_4862_11 = std::sqrt(52490025.0 / 2149004.0);
    const auto fs_7245_4862_154 = std::sqrt(367430175.0 / 1074502.0);
    const auto fs_7245_4862_66 = std::sqrt(157470075.0 / 1074502.0);
    const auto fs_7245_4862_715 = std::sqrt(262450125.0 / 165308.0);
    const auto fs_7245_4862_858 = std::sqrt(157470075.0 / 82654.0);
    const auto fs_72_4199_1001 = std::sqrt(399168.0 / 1356277.0);
    const auto fs_72_4199_182 = std::sqrt(72576.0 / 1356277.0);
    const auto fs_72_4199_330 = std::sqrt(1710720.0 / 17631601.0);
    const auto fs_72_4199_442 = std::sqrt(10368.0 / 79781.0);
    const auto fs_750_4199_11 = std::sqrt(6187500.0 / 17631601.0);
    const auto fs_75600_46189_33 = std::sqrt(17146080000.0 / 193947611.0);
    const auto fs_7560_221_3 = std::sqrt(171460800.0 / 48841.0);
    const auto fs_7560_2431_35 = std::sqrt(2000376000.0 / 5909761.0);
    const auto fs_7560_2431_385 = std::sqrt(2000376000.0 / 537251.0);
    const auto fs_7560_46189_2310 = std::sqrt(12002256000.0 / 193947611.0);
    const auto fs_75_4199_11 = std::sqrt(61875.0 / 17631601.0);
    const auto fs_75_4199_138567 = std::sqrt(185625.0 / 4199.0);
    const auto fs_75_4199_143 = std::sqrt(61875.0 / 1356277.0);
    const auto fs_75_4199_146965 = std::sqrt(196875.0 / 4199.0);
    const auto fs_75_4199_176358 = std::sqrt(236250.0 / 4199.0);
    const auto fs_75_4199_210 = std::sqrt(1181250.0 / 17631601.0);
    const auto fs_75_4199_25194 = std::sqrt(33750.0 / 4199.0);
    const auto fs_75_4199_30030 = std::sqrt(12993750.0 / 1356277.0);
    const auto fs_75_4199_323323 = std::sqrt(433125.0 / 4199.0);
    const auto fs_75_4199_33 = std::sqrt(185625.0 / 17631601.0);
    const auto fs_75_4199_46189 = std::sqrt(61875.0 / 4199.0);
    const auto fs_75_4199_62985 = std::sqrt(84375.0 / 4199.0);
    const auto fs_75_4199_910 = std::sqrt(393750.0 / 1356277.0);
    const auto fs_75_4199_92378 = std::sqrt(123750.0 / 4199.0);
    const auto fs_75_4199_9282 = std::sqrt(236250.0 / 79781.0);
    const auto fs_75_8398_182 = std::sqrt(39375.0 / 2712554.0);
    const auto fs_75_8398_2 = std::sqrt(5625.0 / 35263202.0);
    const auto fs_75_8398_2002 = std::sqrt(433125.0 / 2712554.0);
    const auto fs_75_8398_26 = std::sqrt(5625.0 / 2712554.0);
    const auto fs_75_8398_6006 = std::sqrt(1299375.0 / 2712554.0);
    const auto fs_75_8398_910 = std::sqrt(196875.0 / 2712554.0);
    const auto fs_760_13_35 = std::sqrt(20216000.0 / 169.0);
    const auto fs_760_221_11 = std::sqrt(6353600.0 / 48841.0);
    const auto fs_760_221_7 = std::sqrt(4043200.0 / 48841.0);
    const auto fs_760_663_42 = std::sqrt(8086400.0 / 146523.0);
    const auto fs_76_13_210 = std::sqrt(1212960.0 / 169.0);
    const auto fs_76_13_462 = std::sqrt(2668512.0 / 169.0);
    const auto fs_76_221_10 = std::sqrt(57760.0 / 48841.0);
    const auto fs_76_221_210 = std::sqrt(1212960.0 / 48841.0);
    const auto fs_76_51_6 = std::sqrt(11552.0 / 867.0);
    const auto fs_7875_1144_770 = std::sqrt(2170546875.0 / 59488.0);
    const auto fs_79380_46189_30 = std::sqrt(189035532000.0 / 2133423721.0);
    const auto fs_7_143_10 = std::sqrt(490.0 / 20449.0);
    const auto fs_7_143_22 = std::sqrt(98.0 / 1859.0);
    const auto fs_7_221_78 = std::sqrt(294.0 / 3757.0);
    const auto fs_7_2431_1430 = std::sqrt(490.0 / 41327.0);
    const auto fs_7_429_30 = std::sqrt(490.0 / 61347.0);
    const auto fs_7_429_6 = std::sqrt(98.0 / 61347.0);
    const auto fs_7_429_66 = std::sqrt(98.0 / 5577.0);
    const auto fs_810_143_2 = std::sqrt(1312200.0 / 20449.0);
    const auto fs_810_143_3 = std::sqrt(1968300.0 / 20449.0);
    const auto fs_825_4199_5 = std::sqrt(3403125.0 / 17631601.0);
    const auto fs_828_2431_35 = std::sqrt(23995440.0 / 5909761.0);
    const auto fs_82_2431_66 = std::sqrt(40344.0 / 537251.0);
    const auto fs_836_663_42 = std::sqrt(9784544.0 / 146523.0);
    const auto fs_8441_2431_2 = std::sqrt(142500962.0 / 5909761.0);
    const auto fs_84_2431_143 = std::sqrt(7056.0 / 41327.0);
    const auto fs_84_4199_429 = std::sqrt(232848.0 / 1356277.0);
    const auto fs_8505_1144_462 = std::sqrt(1519035525.0 / 59488.0);
    const auto fs_8505_16_110 = std::sqrt(3978426375.0 / 128.0);
    const auto fs_8505_16_2 = std::sqrt(72335025.0 / 128.0);
    const auto fs_8505_16_6 = std::sqrt(217005075.0 / 128.0);
    const auto fs_8505_286_30 = std::sqrt(1085025375.0 / 40898.0);
    const auto fs_8505_2_2 = std::sqrt(72335025.0 / 2.0);
    const auto fs_8505_46189_105 = std::sqrt(7595177625.0 / 2133423721.0);
    const auto fs_8505_46189_385 = std::sqrt(2531725875.0 / 193947611.0);
    const auto fs_8505_46189_5 = std::sqrt(361675125.0 / 2133423721.0);
    const auto fs_8505_46189_770 = std::sqrt(5063451750.0 / 193947611.0);
    const auto fs_8505_4862_105 = std::sqrt(7595177625.0 / 23639044.0);
    const auto fs_8505_4862_385 = std::sqrt(2531725875.0 / 2149004.0);
    const auto fs_8505_4862_5 = std::sqrt(361675125.0 / 23639044.0);
    const auto fs_8505_4862_770 = std::sqrt(2531725875.0 / 1074502.0);
    const auto fs_8505_4_5 = std::sqrt(361675125.0 / 16.0);
    const auto fs_8505_8_11 = std::sqrt(795685275.0 / 64.0);
    const auto fs_8505_8_14 = std::sqrt(506345175.0 / 32.0);
    const auto fs_8505_8_15 = std::sqrt(1085025375.0 / 64.0);
    const auto fs_8505_8_21 = std::sqrt(1519035525.0 / 64.0);
    const auto fs_8505_8_3 = std::sqrt(217005075.0 / 64.0);
    const auto fs_8505_8_33 = std::sqrt(2387055825.0 / 64.0);
    const auto fs_8505_8_35 = std::sqrt(2531725875.0 / 64.0);
    const auto fs_8505_8_5 = std::sqrt(361675125.0 / 64.0);
    const auto fs_855_13_66 = std::sqrt(48247650.0 / 169.0);
    const auto fs_855_8_66 = std::sqrt(24123825.0 / 32.0);
    const auto fs_8610_46189_2 = std::sqrt(148264200.0 / 2133423721.0);
    const auto fs_87_2431_30 = std::sqrt(227070.0 / 5909761.0);
    const auto fs_8925_143_11 = std::sqrt(79655625.0 / 1859.0);
    const auto fs_8925_143_3 = std::sqrt(238966875.0 / 20449.0);
    const auto fs_8925_143_5 = std::sqrt(398278125.0 / 20449.0);
    const auto fs_8925_143_7 = std::sqrt(557589375.0 / 20449.0);
    const auto fs_8925_16_10 = std::sqrt(398278125.0 / 128.0);
    const auto fs_8925_16_22 = std::sqrt(876211875.0 / 128.0);
    const auto fs_8925_22_11 = std::sqrt(79655625.0 / 44.0);
    const auto fs_8925_22_3 = std::sqrt(238966875.0 / 484.0);
    const auto fs_8925_22_5 = std::sqrt(398278125.0 / 484.0);
    const auto fs_8925_22_7 = std::sqrt(557589375.0 / 484.0);
    const auto fs_8925_286_10 = std::sqrt(398278125.0 / 40898.0);
    const auto fs_8925_286_22 = std::sqrt(79655625.0 / 3718.0);
    const auto fs_8925_44_10 = std::sqrt(398278125.0 / 968.0);
    const auto fs_8925_44_22 = std::sqrt(79655625.0 / 88.0);
    const auto fs_8925_8_11 = std::sqrt(876211875.0 / 64.0);
    const auto fs_8925_8_3 = std::sqrt(238966875.0 / 64.0);
    const auto fs_8925_8_5 = std::sqrt(398278125.0 / 64.0);
    const auto fs_8925_8_7 = std::sqrt(557589375.0 / 64.0);
    const auto fs_8_221_91 = std::sqrt(448.0 / 3757.0);
    const auto fs_8_663_210 = std::sqrt(4480.0 / 146523.0);
    const auto fs_8_663_462 = std::sqrt(9856.0 / 146523.0);
    const auto fs_900_1_2 = std::sqrt(1620000.0);
    const auto fs_900_4199_1001 = std::sqrt(62370000.0 / 1356277.0);
    const auto fs_900_4199_182 = std::sqrt(11340000.0 / 1356277.0);
    const auto fs_900_4199_330 = std::sqrt(267300000.0 / 17631601.0);
    const auto fs_900_4199_442 = std::sqrt(1620000.0 / 79781.0);
    const auto fs_90405_572_2 = std::sqrt(8173064025.0 / 163592.0);
    const auto fs_90_3553_70 = std::sqrt(567000.0 / 12623809.0);
    const auto fs_90_4199_231 = std::sqrt(1871100.0 / 17631601.0);
    const auto fs_90_4199_429 = std::sqrt(267300.0 / 1356277.0);
    const auto fs_90_46189_66 = std::sqrt(48600.0 / 193947611.0);
    const auto fs_90_46189_77 = std::sqrt(56700.0 / 193947611.0);
    const auto fs_9135_1144_462 = std::sqrt(1752412725.0 / 59488.0);
    const auto fs_9135_1144_70 = std::sqrt(2920687875.0 / 654368.0);
    const auto fs_92610_46189_7 = std::sqrt(60036284700.0 / 2133423721.0);
    const auto fs_92_143_105 = std::sqrt(888720.0 / 20449.0);
    const auto fs_92_221_14 = std::sqrt(118496.0 / 48841.0);
    const auto fs_92_221_3 = std::sqrt(25392.0 / 48841.0);
    const auto fs_92_221_39 = std::sqrt(25392.0 / 3757.0);
    const auto fs_92_221_442 = std::sqrt(16928.0 / 221.0);
    const auto fs_92_221_455 = std::sqrt(296240.0 / 3757.0);
    const auto fs_92_2431_462 = std::sqrt(355488.0 / 537251.0);
    const auto fs_9450_2717_3 = std::sqrt(267907500.0 / 7382089.0);
    const auto fs_9450_46189_1001 = std::sqrt(625117500.0 / 14919047.0);
    const auto fs_945_221_105 = std::sqrt(93767625.0 / 48841.0);
    const auto fs_945_286_2002 = std::sqrt(6251175.0 / 286.0);
    const auto fs_945_3553_1001 = std::sqrt(81265275.0 / 1147619.0);
    const auto fs_945_374_1001 = std::sqrt(81265275.0 / 12716.0);
    const auto fs_945_4199_770 = std::sqrt(687629250.0 / 17631601.0);
    const auto fs_945_442_770 = std::sqrt(343814625.0 / 97682.0);
    const auto fs_945_44_70 = std::sqrt(31255875.0 / 968.0);
    const auto fs_945_572_66 = std::sqrt(2679075.0 / 14872.0);
    const auto fs_945_572_77 = std::sqrt(6251175.0 / 29744.0);
    const auto fs_945_8398_2730 = std::sqrt(93767625.0 / 2712554.0);
    const auto fs_945_884_2730 = std::sqrt(93767625.0 / 30056.0);
    const auto fs_945_92378_10010 = std::sqrt(31255875.0 / 29838094.0);
    const auto fs_945_9724_10010 = std::sqrt(31255875.0 / 330616.0);
    const auto fs_950_13_3 = std::sqrt(2707500.0 / 169.0);
    const auto fs_950_13_42 = std::sqrt(37905000.0 / 169.0);
    const auto fs_95_1_35 = std::sqrt(315875.0);
    const auto fs_95_1_6 = std::sqrt(54150.0);
    const auto fs_95_2_210 = std::sqrt(947625.0 / 2.0);
    const auto fs_95_2_462 = std::sqrt(2084775.0 / 2.0);
    const auto fs_95_4_14 = std::sqrt(63175.0 / 8.0);
    const auto fs_95_4_231 = std::sqrt(2084775.0 / 16.0);
    const auto fs_966_221_91 = std::sqrt(6532092.0 / 3757.0);
    const auto fs_966_2431_165 = std::sqrt(13997340.0 / 537251.0);
    const auto fs_98_2431_33 = std::sqrt(28812.0 / 537251.0);
    const auto fs_990_4199_2 = std::sqrt(1960200.0 / 17631601.0);
    const auto fs_9_143_10 = std::sqrt(810.0 / 20449.0);
    const auto fs_9_2431_4290 = std::sqrt(2430.0 / 41327.0);
    const auto fs_9_4199_1430 = std::sqrt(8910.0 / 1356277.0);
    const auto fs_9_4199_30030 = std::sqrt(187110.0 / 1356277.0);

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph11_p1, ph11_p11, ab_2, pc_0 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p11 = ph11_p11[k];

        pc_0[k] = - e_0 * fs_14175_32_33 * h1_p1 - e_1 * fs_16065_32_22 * h3_p1 + e_1 * fs_8505_8_33 * r_2 * h1_p1 - e_2 * fs_285_4_55 * h5_p1 + e_2 * fs_8925_16_22 * r_2 * h3_p1 - e_2 * fs_6075_8_33 * r_4 * h1_p1 - e_3 * fs_1575_572_231 * h7_p1 + e_3 * fs_570_13_55 * r_2 * h5_p1 - e_3 * fs_8925_44_22 * r_4 * h3_p1 + e_3 * fs_225_1_33 * r_6 * h1_p1 - e_4 * fs_483_4862_165 * h9_p1 + e_4 * fs_4725_4862_231 * r_2 * h7_p1 - e_4 * fs_114_13_55 * r_4 * h5_p1 + e_4 * fs_8925_286_22 * r_6 * h3_p1 - e_4 * fs_675_22_33 * r_8 * h1_p1 - e_5 * fs_75_8398_2 * h11_p1 - e_5 * fs_75_4199_323323 * h11_p11 + e_5 * fs_46_2431_165 * r_2 * h9_p1 - e_5 * fs_4725_46189_231 * r_4 * h7_p1 + e_5 * fs_152_221_55 * r_6 * h5_p1 - e_5 * fs_595_286_22 * r_8 * h3_p1 + e_5 * fs_270_143_33 * r_10 * h1_p1 + e_6 * fs_3_4199_2 * r_2 * h11_p1 + e_6 * fs_6_4199_323323 * r_2 * h11_p11 - e_6 * fs_2_2431_165 * r_4 * h9_p1 + e_6 * fs_150_46189_231 * r_6 * h7_p1 - e_6 * fs_4_221_55 * r_8 * h5_p1 + e_6 * fs_7_143_22 * r_10 * h3_p1 - e_6 * fs_6_143_33 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ph11_p10, ab_2, pc_1 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p10 = ph11_p10[k];

        pc_1[k] = e_1 * fs_16065_32_22 * h3_p2 + e_2 * fs_285_4_154 * h5_p2 - e_2 * fs_8925_16_22 * r_2 * h3_p2 + e_3 * fs_2835_572_385 * h7_p2 - e_3 * fs_570_13_154 * r_2 * h5_p2 + e_3 * fs_8925_44_22 * r_4 * h3_p2 + e_4 * fs_483_221_3 * h9_p2 - e_4 * fs_8505_4862_385 * r_2 * h7_p2 + e_4 * fs_114_13_154 * r_4 * h5_p2 - e_4 * fs_8925_286_22 * r_6 * h3_p2 + e_5 * fs_75_8398_26 * h11_p2 - e_5 * fs_75_4199_146965 * h11_p10 - e_5 * fs_92_221_3 * r_2 * h9_p2 + e_5 * fs_8505_46189_385 * r_4 * h7_p2 - e_5 * fs_152_221_154 * r_6 * h5_p2 + e_5 * fs_595_286_22 * r_8 * h3_p2 - e_6 * fs_3_4199_26 * r_2 * h11_p2 + e_6 * fs_6_4199_146965 * r_2 * h11_p10 + e_6 * fs_4_221_3 * r_4 * h9_p2 - e_6 * fs_270_46189_385 * r_6 * h7_p2 + e_6 * fs_4_221_154 * r_8 * h5_p2 - e_6 * fs_7_143_22 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph9_p9, ph11_p3, ph11_p9, ab_2, pc_2 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p9 = ph11_p9[k];

        pc_2[k] = - e_1 * fs_5355_32_66 * h3_p3 - e_2 * fs_95_2_462 * h5_p3 + e_2 * fs_2975_16_66 * r_2 * h3_p3 - e_3 * fs_4725_572_385 * h7_p3 + e_3 * fs_380_13_462 * r_2 * h5_p3 - e_3 * fs_2975_44_66 * r_4 * h3_p3 - e_4 * fs_483_221_14 * h9_p3 - e_4 * fs_483_442_1326 * h9_p9 + e_4 * fs_14175_4862_385 * r_2 * h7_p3 - e_4 * fs_76_13_462 * r_4 * h5_p3 + e_4 * fs_2975_286_66 * r_6 * h3_p3 - e_5 * fs_75_8398_182 * h11_p3 - e_5 * fs_75_4199_62985 * h11_p9 + e_5 * fs_92_221_14 * r_2 * h9_p3 + e_5 * fs_46_221_1326 * r_2 * h9_p9 - e_5 * fs_14175_46189_385 * r_4 * h7_p3 + e_5 * fs_304_663_462 * r_6 * h5_p3 - e_5 * fs_595_858_66 * r_8 * h3_p3 + e_6 * fs_3_4199_182 * r_2 * h11_p3 + e_6 * fs_6_4199_62985 * r_2 * h11_p9 - e_6 * fs_4_221_14 * r_4 * h9_p3 - e_6 * fs_2_221_1326 * r_4 * h9_p9 + e_6 * fs_450_46189_385 * r_6 * h7_p3 - e_6 * fs_8_663_462 * r_8 * h5_p3 + e_6 * fs_7_429_66 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p4, ph7_p4, ph9_p4, ph9_p8, ph11_p4, ph11_p8, ab_2, pc_3 : simd::cache_line_size())
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

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p8 = ph11_p8[k];

        pc_3[k] = e_2 * fs_285_4_154 * h5_p4 + e_3 * fs_1575_104_210 * h7_p4 - e_3 * fs_570_13_154 * r_2 * h5_p4 + e_4 * fs_483_442_182 * h9_p4 - e_4 * fs_483_221_442 * h9_p8 - e_4 * fs_4725_884_210 * r_2 * h7_p4 + e_4 * fs_114_13_154 * r_4 * h5_p4 + e_5 * fs_75_8398_910 * h11_p4 - e_5 * fs_75_4199_25194 * h11_p8 - e_5 * fs_46_221_182 * r_2 * h9_p4 + e_5 * fs_92_221_442 * r_2 * h9_p8 + e_5 * fs_4725_8398_210 * r_4 * h7_p4 - e_5 * fs_152_221_154 * r_6 * h5_p4 - e_6 * fs_3_4199_910 * r_2 * h11_p4 + e_6 * fs_6_4199_25194 * r_2 * h11_p8 + e_6 * fs_2_221_182 * r_4 * h9_p4 - e_6 * fs_4_221_442 * r_4 * h9_p8 - e_6 * fs_75_4199_210 * r_6 * h7_p4 + e_6 * fs_4_221_154 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p5, ph7_m6, ph7_p5, ph7_p7, ph9_m6, ph9_p5, ph9_p7, ph11_m6, ph11_p5, ph11_p7, ab_2, pc_4, pc_5 : simd::cache_line_size())
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

        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_p5 = ph11_p5[k];
        const auto h11_p7 = ph11_p7[k];

        pc_4[k] = - e_2 * fs_285_4_55 * h5_p5 - e_3 * fs_4725_104_30 * h7_p5 - e_3 * fs_315_104_2730 * h7_p7 + e_3 * fs_570_13_55 * r_2 * h5_p5 - e_4 * fs_483_442_455 * h9_p5 - e_4 * fs_966_221_91 * h9_p7 + e_4 * fs_14175_884_30 * r_2 * h7_p5 + e_4 * fs_945_884_2730 * r_2 * h7_p7 - e_4 * fs_114_13_55 * r_4 * h5_p5 - e_5 * fs_75_4199_910 * h11_p5 - e_5 * fs_75_4199_9282 * h11_p7 + e_5 * fs_46_221_455 * r_2 * h9_p5 + e_5 * fs_184_221_91 * r_2 * h9_p7 - e_5 * fs_14175_8398_30 * r_4 * h7_p5 - e_5 * fs_945_8398_2730 * r_4 * h7_p7 + e_5 * fs_152_221_55 * r_6 * h5_p5 + e_6 * fs_6_4199_910 * r_2 * h11_p5 + e_6 * fs_6_4199_9282 * r_2 * h11_p7 - e_6 * fs_2_221_455 * r_4 * h9_p5 - e_6 * fs_8_221_91 * r_4 * h9_p7 + e_6 * fs_225_4199_30 * r_6 * h7_p5 + e_6 * fs_15_4199_2730 * r_6 * h7_p7 - e_6 * fs_4_221_55 * r_8 * h5_p5;

        pc_5[k] = e_3 * fs_4725_52_13 * h7_m6 + e_4 * fs_483_221_455 * h9_m6 - e_4 * fs_14175_442_13 * r_2 * h7_m6 + e_5 * fs_150_4199_1547 * h11_m6 - e_5 * fs_92_221_455 * r_2 * h9_m6 + e_5 * fs_14175_4199_13 * r_4 * h7_m6 - e_6 * fs_12_4199_1547 * r_2 * h11_m6 + e_6 * fs_4_221_455 * r_4 * h9_m6 - e_6 * fs_450_4199_13 * r_6 * h7_m6;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m5, ph7_m7, ph7_m5, ph9_m7, ph9_m5, ph11_m7, ph11_m5, ab_2, pc_6 : simd::cache_line_size())
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

        const auto h5_m5 = ph5_m5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m5 = ph11_m5[k];

        pc_6[k] = - e_2 * fs_285_4_55 * h5_m5 + e_3 * fs_315_104_2730 * h7_m7 - e_3 * fs_4725_104_30 * h7_m5 + e_3 * fs_570_13_55 * r_2 * h5_m5 + e_4 * fs_966_221_91 * h9_m7 - e_4 * fs_483_442_455 * h9_m5 - e_4 * fs_945_884_2730 * r_2 * h7_m7 + e_4 * fs_14175_884_30 * r_2 * h7_m5 - e_4 * fs_114_13_55 * r_4 * h5_m5 + e_5 * fs_75_4199_9282 * h11_m7 - e_5 * fs_75_4199_910 * h11_m5 - e_5 * fs_184_221_91 * r_2 * h9_m7 + e_5 * fs_46_221_455 * r_2 * h9_m5 + e_5 * fs_945_8398_2730 * r_4 * h7_m7 - e_5 * fs_14175_8398_30 * r_4 * h7_m5 + e_5 * fs_152_221_55 * r_6 * h5_m5 - e_6 * fs_6_4199_9282 * r_2 * h11_m7 + e_6 * fs_6_4199_910 * r_2 * h11_m5 + e_6 * fs_8_221_91 * r_4 * h9_m7 - e_6 * fs_2_221_455 * r_4 * h9_m5 - e_6 * fs_15_4199_2730 * r_6 * h7_m7 + e_6 * fs_225_4199_30 * r_6 * h7_m5 - e_6 * fs_4_221_55 * r_8 * h5_m5;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph7_m4, ph9_m8, ph9_m4, ph11_m8, ph11_m4, ab_2, pc_7 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m4 = ph11_m4[k];

        pc_7[k] = e_2 * fs_285_4_154 * h5_m4 + e_3 * fs_1575_104_210 * h7_m4 - e_3 * fs_570_13_154 * r_2 * h5_m4 + e_4 * fs_483_221_442 * h9_m8 + e_4 * fs_483_442_182 * h9_m4 - e_4 * fs_4725_884_210 * r_2 * h7_m4 + e_4 * fs_114_13_154 * r_4 * h5_m4 + e_5 * fs_75_4199_25194 * h11_m8 + e_5 * fs_75_8398_910 * h11_m4 - e_5 * fs_92_221_442 * r_2 * h9_m8 - e_5 * fs_46_221_182 * r_2 * h9_m4 + e_5 * fs_4725_8398_210 * r_4 * h7_m4 - e_5 * fs_152_221_154 * r_6 * h5_m4 - e_6 * fs_6_4199_25194 * r_2 * h11_m8 - e_6 * fs_3_4199_910 * r_2 * h11_m4 + e_6 * fs_4_221_442 * r_4 * h9_m8 + e_6 * fs_2_221_182 * r_4 * h9_m4 - e_6 * fs_75_4199_210 * r_6 * h7_m4 + e_6 * fs_4_221_154 * r_8 * h5_m4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m3, ph7_m3, ph9_m9, ph9_m3, ph11_m9, ph11_m3, ab_2, pc_8 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m3 = ph11_m3[k];

        pc_8[k] = - e_1 * fs_5355_32_66 * h3_m3 - e_2 * fs_95_2_462 * h5_m3 + e_2 * fs_2975_16_66 * r_2 * h3_m3 - e_3 * fs_4725_572_385 * h7_m3 + e_3 * fs_380_13_462 * r_2 * h5_m3 - e_3 * fs_2975_44_66 * r_4 * h3_m3 + e_4 * fs_483_442_1326 * h9_m9 - e_4 * fs_483_221_14 * h9_m3 + e_4 * fs_14175_4862_385 * r_2 * h7_m3 - e_4 * fs_76_13_462 * r_4 * h5_m3 + e_4 * fs_2975_286_66 * r_6 * h3_m3 + e_5 * fs_75_4199_62985 * h11_m9 - e_5 * fs_75_8398_182 * h11_m3 - e_5 * fs_46_221_1326 * r_2 * h9_m9 + e_5 * fs_92_221_14 * r_2 * h9_m3 - e_5 * fs_14175_46189_385 * r_4 * h7_m3 + e_5 * fs_304_663_462 * r_6 * h5_m3 - e_5 * fs_595_858_66 * r_8 * h3_m3 - e_6 * fs_6_4199_62985 * r_2 * h11_m9 + e_6 * fs_3_4199_182 * r_2 * h11_m3 + e_6 * fs_2_221_1326 * r_4 * h9_m9 - e_6 * fs_4_221_14 * r_4 * h9_m3 + e_6 * fs_450_46189_385 * r_6 * h7_m3 - e_6 * fs_8_663_462 * r_8 * h5_m3 + e_6 * fs_7_429_66 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m10, ph11_m2, ab_2, pc_9 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m2 = ph11_m2[k];

        pc_9[k] = e_1 * fs_16065_32_22 * h3_m2 + e_2 * fs_285_4_154 * h5_m2 - e_2 * fs_8925_16_22 * r_2 * h3_m2 + e_3 * fs_2835_572_385 * h7_m2 - e_3 * fs_570_13_154 * r_2 * h5_m2 + e_3 * fs_8925_44_22 * r_4 * h3_m2 + e_4 * fs_483_221_3 * h9_m2 - e_4 * fs_8505_4862_385 * r_2 * h7_m2 + e_4 * fs_114_13_154 * r_4 * h5_m2 - e_4 * fs_8925_286_22 * r_6 * h3_m2 + e_5 * fs_75_4199_146965 * h11_m10 + e_5 * fs_75_8398_26 * h11_m2 - e_5 * fs_92_221_3 * r_2 * h9_m2 + e_5 * fs_8505_46189_385 * r_4 * h7_m2 - e_5 * fs_152_221_154 * r_6 * h5_m2 + e_5 * fs_595_286_22 * r_8 * h3_m2 - e_6 * fs_6_4199_146965 * r_2 * h11_m10 - e_6 * fs_3_4199_26 * r_2 * h11_m2 + e_6 * fs_4_221_3 * r_4 * h9_m2 - e_6 * fs_270_46189_385 * r_6 * h7_m2 + e_6 * fs_4_221_154 * r_8 * h5_m2 - e_6 * fs_7_143_22 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m1, ph11_m11, ph11_m1, ab_2, pc_10 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m11 = ph11_m11[k];
        const auto h11_m1 = ph11_m1[k];

        pc_10[k] = - e_0 * fs_14175_32_33 * h1_m1 - e_1 * fs_16065_32_22 * h3_m1 + e_1 * fs_8505_8_33 * r_2 * h1_m1 - e_2 * fs_285_4_55 * h5_m1 + e_2 * fs_8925_16_22 * r_2 * h3_m1 - e_2 * fs_6075_8_33 * r_4 * h1_m1 - e_3 * fs_1575_572_231 * h7_m1 + e_3 * fs_570_13_55 * r_2 * h5_m1 - e_3 * fs_8925_44_22 * r_4 * h3_m1 + e_3 * fs_225_1_33 * r_6 * h1_m1 - e_4 * fs_483_4862_165 * h9_m1 + e_4 * fs_4725_4862_231 * r_2 * h7_m1 - e_4 * fs_114_13_55 * r_4 * h5_m1 + e_4 * fs_8925_286_22 * r_6 * h3_m1 - e_4 * fs_675_22_33 * r_8 * h1_m1 + e_5 * fs_75_4199_323323 * h11_m11 - e_5 * fs_75_8398_2 * h11_m1 + e_5 * fs_46_2431_165 * r_2 * h9_m1 - e_5 * fs_4725_46189_231 * r_4 * h7_m1 + e_5 * fs_152_221_55 * r_6 * h5_m1 - e_5 * fs_595_286_22 * r_8 * h3_m1 + e_5 * fs_270_143_33 * r_10 * h1_m1 - e_6 * fs_6_4199_323323 * r_2 * h11_m11 + e_6 * fs_3_4199_2 * r_2 * h11_m1 - e_6 * fs_2_2431_165 * r_4 * h9_m1 + e_6 * fs_150_46189_231 * r_6 * h7_m1 - e_6 * fs_4_221_55 * r_8 * h5_m1 + e_6 * fs_7_143_22 * r_10 * h3_m1 - e_6 * fs_6_143_33 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph11_0, ph11_p10, ab_2, pc_11 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p10 = ph11_p10[k];

        pc_11[k] = - e_0 * fs_14175_32_11 * h1_0 - e_1 * fs_16065_16_11 * h3_0 + e_1 * fs_8505_8_11 * r_2 * h1_0 - e_2 * fs_1425_4_11 * h5_0 + e_2 * fs_8925_8_11 * r_2 * h3_0 - e_2 * fs_6075_8_11 * r_4 * h1_0 - e_3 * fs_11025_286_11 * h7_0 + e_3 * fs_2850_13_11 * r_2 * h5_0 - e_3 * fs_8925_22_11 * r_4 * h3_0 + e_3 * fs_225_1_11 * r_6 * h1_0 - e_4 * fs_7245_4862_11 * h9_0 + e_4 * fs_33075_2431_11 * r_2 * h7_0 - e_4 * fs_570_13_11 * r_4 * h5_0 + e_4 * fs_8925_143_11 * r_6 * h3_0 - e_4 * fs_675_22_11 * r_8 * h1_0 - e_5 * fs_75_4199_11 * h11_0 - e_5 * fs_75_4199_176358 * h11_p10 + e_5 * fs_690_2431_11 * r_2 * h9_0 - e_5 * fs_66150_46189_11 * r_4 * h7_0 + e_5 * fs_760_221_11 * r_6 * h5_0 - e_5 * fs_595_143_11 * r_8 * h3_0 + e_5 * fs_270_143_11 * r_10 * h1_0 + e_6 * fs_6_4199_11 * r_2 * h11_0 + e_6 * fs_6_4199_176358 * r_2 * h11_p10 - e_6 * fs_30_2431_11 * r_4 * h9_0 + e_6 * fs_2100_46189_11 * r_6 * h7_0 - e_6 * fs_20_221_11 * r_8 * h5_0 + e_6 * fs_14_143_11 * r_10 * h3_0 - e_6 * fs_6_143_11 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2, pc_12 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_12[k] = - e_0 * fs_14175_64_110 * h1_p1 + e_1 * fs_8505_16_110 * r_2 * h1_p1 + e_2 * fs_855_8_66 * h5_p1 - e_2 * fs_6075_16_110 * r_4 * h1_p1 + e_3 * fs_315_52_770 * h7_p1 - e_3 * fs_855_13_66 * r_2 * h5_p1 + e_3 * fs_225_2_110 * r_6 * h1_p1 + e_4 * fs_1449_748_22 * h9_p1 + e_4 * fs_1449_442_221 * h9_p9 - e_4 * fs_945_442_770 * r_2 * h7_p1 + e_4 * fs_171_13_66 * r_4 * h5_p1 - e_4 * fs_675_44_110 * r_8 * h1_p1 + e_5 * fs_150_4199_15 * h11_p1 - e_5 * fs_150_4199_41990 * h11_p9 - e_5 * fs_69_187_22 * r_2 * h9_p1 - e_5 * fs_138_221_221 * r_2 * h9_p9 + e_5 * fs_945_4199_770 * r_4 * h7_p1 - e_5 * fs_228_221_66 * r_6 * h5_p1 + e_5 * fs_135_143_110 * r_10 * h1_p1 - e_6 * fs_12_4199_15 * r_2 * h11_p1 + e_6 * fs_12_4199_41990 * r_2 * h11_p9 + e_6 * fs_3_187_22 * r_4 * h9_p1 + e_6 * fs_6_221_221 * r_4 * h9_p9 - e_6 * fs_30_4199_770 * r_6 * h7_p1 + e_6 * fs_6_221_66 * r_8 * h5_p1 - e_6 * fs_3_143_110 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2, pc_13 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_13[k] = e_1 * fs_5355_16_33 * h3_p2 - e_2 * fs_95_4_231 * h5_p2 - e_2 * fs_2975_8_33 * r_2 * h3_p2 - e_3 * fs_630_143_2310 * h7_p2 + e_3 * fs_190_13_231 * r_2 * h5_p2 + e_3 * fs_2975_22_33 * r_4 * h3_p2 - e_4 * fs_5313_442_2 * h9_p2 + e_4 * fs_483_442_221 * h9_p8 + e_4 * fs_3780_2431_2310 * r_2 * h7_p2 - e_4 * fs_38_13_231 * r_4 * h5_p2 - e_4 * fs_2975_143_33 * r_6 * h3_p2 - e_5 * fs_225_4199_39 * h11_p2 - e_5 * fs_225_4199_12597 * h11_p8 + e_5 * fs_506_221_2 * r_2 * h9_p2 - e_5 * fs_46_221_221 * r_2 * h9_p8 - e_5 * fs_7560_46189_2310 * r_4 * h7_p2 + e_5 * fs_152_663_231 * r_6 * h5_p2 + e_5 * fs_595_429_33 * r_8 * h3_p2 + e_6 * fs_18_4199_39 * r_2 * h11_p2 + e_6 * fs_18_4199_12597 * r_2 * h11_p8 - e_6 * fs_22_221_2 * r_4 * h9_p2 + e_6 * fs_2_221_221 * r_4 * h9_p8 + e_6 * fs_240_46189_2310 * r_6 * h7_p2 - e_6 * fs_4_663_231 * r_8 * h5_p2 - e_6 * fs_14_429_33 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2, pc_14 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_14[k] = - e_1 * fs_5355_16_33 * h3_p3 - e_2 * fs_95_4_231 * h5_p3 + e_2 * fs_2975_8_33 * r_2 * h3_p3 + e_3 * fs_7875_1144_770 * h7_p3 + e_3 * fs_2205_104_130 * h7_p7 + e_3 * fs_190_13_231 * r_2 * h5_p3 - e_3 * fs_2975_22_33 * r_4 * h3_p3 + e_4 * fs_4347_442_7 * h9_p3 - e_4 * fs_483_221_39 * h9_p7 - e_4 * fs_23625_9724_770 * r_2 * h7_p3 - e_4 * fs_6615_884_130 * r_2 * h7_p7 - e_4 * fs_38_13_231 * r_4 * h5_p3 + e_4 * fs_2975_143_33 * r_6 * h3_p3 + e_5 * fs_300_4199_91 * h11_p3 - e_5 * fs_900_4199_442 * h11_p7 - e_5 * fs_414_221_7 * r_2 * h9_p3 + e_5 * fs_92_221_39 * r_2 * h9_p7 + e_5 * fs_23625_92378_770 * r_4 * h7_p3 + e_5 * fs_6615_8398_130 * r_4 * h7_p7 + e_5 * fs_152_663_231 * r_6 * h5_p3 - e_5 * fs_595_429_33 * r_8 * h3_p3 - e_6 * fs_24_4199_91 * r_2 * h11_p3 + e_6 * fs_72_4199_442 * r_2 * h11_p7 + e_6 * fs_18_221_7 * r_4 * h9_p3 - e_6 * fs_4_221_39 * r_4 * h9_p7 - e_6 * fs_375_46189_770 * r_6 * h7_p3 - e_6 * fs_105_4199_130 * r_6 * h7_p7 - e_6 * fs_4_663_231 * r_8 * h5_p3 + e_6 * fs_14_429_33 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m5, ph5_p4, ph7_m5, ph7_p4, ph7_p6, ph9_m5, ph9_p4, ph9_p6, ph11_m5, ph11_p4, ph11_p6, ab_2, pc_15, pc_16 : simd::cache_line_size())
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

        const auto h5_m5 = ph5_m5[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_15[k] = e_2 * fs_855_8_66 * h5_p4 - e_3 * fs_1575_52_10 * h7_p4 + e_3 * fs_315_13_65 * h7_p6 - e_3 * fs_855_13_66 * r_2 * h5_p4 - e_4 * fs_3381_884_78 * h9_p4 - e_4 * fs_1449_442_91 * h9_p6 + e_4 * fs_4725_442_10 * r_2 * h7_p4 - e_4 * fs_1890_221_65 * r_2 * h7_p6 + e_4 * fs_171_13_66 * r_4 * h5_p4 - e_5 * fs_525_8398_390 * h11_p4 - e_5 * fs_150_4199_7735 * h11_p6 + e_5 * fs_161_221_78 * r_2 * h9_p4 + e_5 * fs_138_221_91 * r_2 * h9_p6 - e_5 * fs_4725_4199_10 * r_4 * h7_p4 + e_5 * fs_3780_4199_65 * r_4 * h7_p6 - e_5 * fs_228_221_66 * r_6 * h5_p4 + e_6 * fs_21_4199_390 * r_2 * h11_p4 + e_6 * fs_12_4199_7735 * r_2 * h11_p6 - e_6 * fs_7_221_78 * r_4 * h9_p4 - e_6 * fs_6_221_91 * r_4 * h9_p6 + e_6 * fs_150_4199_10 * r_6 * h7_p4 - e_6 * fs_120_4199_65 * r_6 * h7_p6 + e_6 * fs_6_221_66 * r_8 * h5_p4;

        pc_16[k] = - e_2 * fs_1425_4_11 * h5_m5 - e_3 * fs_1575_52_6 * h7_m5 + e_3 * fs_2850_13_11 * r_2 * h5_m5 + e_4 * fs_2415_442_91 * h9_m5 + e_4 * fs_4725_442_6 * r_2 * h7_m5 - e_4 * fs_570_13_11 * r_4 * h5_m5 + e_5 * fs_900_4199_182 * h11_m5 - e_5 * fs_230_221_91 * r_2 * h9_m5 - e_5 * fs_4725_4199_6 * r_4 * h7_m5 + e_5 * fs_760_221_11 * r_6 * h5_m5 - e_6 * fs_72_4199_182 * r_2 * h11_m5 + e_6 * fs_10_221_91 * r_4 * h9_m5 + e_6 * fs_150_4199_6 * r_6 * h7_m5 - e_6 * fs_20_221_11 * r_8 * h5_m5;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph7_m6, ph7_m4, ph9_m6, ph9_m4, ph11_m6, ph11_m4, ab_2, pc_17 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m4 = ph11_m4[k];

        pc_17[k] = e_2 * fs_855_8_66 * h5_m4 - e_3 * fs_315_13_65 * h7_m6 - e_3 * fs_1575_52_10 * h7_m4 - e_3 * fs_855_13_66 * r_2 * h5_m4 + e_4 * fs_1449_442_91 * h9_m6 - e_4 * fs_3381_884_78 * h9_m4 + e_4 * fs_1890_221_65 * r_2 * h7_m6 + e_4 * fs_4725_442_10 * r_2 * h7_m4 + e_4 * fs_171_13_66 * r_4 * h5_m4 + e_5 * fs_150_4199_7735 * h11_m6 - e_5 * fs_525_8398_390 * h11_m4 - e_5 * fs_138_221_91 * r_2 * h9_m6 + e_5 * fs_161_221_78 * r_2 * h9_m4 - e_5 * fs_3780_4199_65 * r_4 * h7_m6 - e_5 * fs_4725_4199_10 * r_4 * h7_m4 - e_5 * fs_228_221_66 * r_6 * h5_m4 - e_6 * fs_12_4199_7735 * r_2 * h11_m6 + e_6 * fs_21_4199_390 * r_2 * h11_m4 + e_6 * fs_6_221_91 * r_4 * h9_m6 - e_6 * fs_7_221_78 * r_4 * h9_m4 + e_6 * fs_120_4199_65 * r_6 * h7_m6 + e_6 * fs_150_4199_10 * r_6 * h7_m4 + e_6 * fs_6_221_66 * r_8 * h5_m4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2, pc_18 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_18[k] = - e_1 * fs_5355_16_33 * h3_m3 - e_2 * fs_95_4_231 * h5_m3 + e_2 * fs_2975_8_33 * r_2 * h3_m3 - e_3 * fs_2205_104_130 * h7_m7 + e_3 * fs_7875_1144_770 * h7_m3 + e_3 * fs_190_13_231 * r_2 * h5_m3 - e_3 * fs_2975_22_33 * r_4 * h3_m3 + e_4 * fs_483_221_39 * h9_m7 + e_4 * fs_4347_442_7 * h9_m3 + e_4 * fs_6615_884_130 * r_2 * h7_m7 - e_4 * fs_23625_9724_770 * r_2 * h7_m3 - e_4 * fs_38_13_231 * r_4 * h5_m3 + e_4 * fs_2975_143_33 * r_6 * h3_m3 + e_5 * fs_900_4199_442 * h11_m7 + e_5 * fs_300_4199_91 * h11_m3 - e_5 * fs_92_221_39 * r_2 * h9_m7 - e_5 * fs_414_221_7 * r_2 * h9_m3 - e_5 * fs_6615_8398_130 * r_4 * h7_m7 + e_5 * fs_23625_92378_770 * r_4 * h7_m3 + e_5 * fs_152_663_231 * r_6 * h5_m3 - e_5 * fs_595_429_33 * r_8 * h3_m3 - e_6 * fs_72_4199_442 * r_2 * h11_m7 - e_6 * fs_24_4199_91 * r_2 * h11_m3 + e_6 * fs_4_221_39 * r_4 * h9_m7 + e_6 * fs_18_221_7 * r_4 * h9_m3 + e_6 * fs_105_4199_130 * r_6 * h7_m7 - e_6 * fs_375_46189_770 * r_6 * h7_m3 - e_6 * fs_4_663_231 * r_8 * h5_m3 + e_6 * fs_14_429_33 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2, pc_19 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_19[k] = e_1 * fs_5355_16_33 * h3_m2 - e_2 * fs_95_4_231 * h5_m2 - e_2 * fs_2975_8_33 * r_2 * h3_m2 - e_3 * fs_630_143_2310 * h7_m2 + e_3 * fs_190_13_231 * r_2 * h5_m2 + e_3 * fs_2975_22_33 * r_4 * h3_m2 - e_4 * fs_483_442_221 * h9_m8 - e_4 * fs_5313_442_2 * h9_m2 + e_4 * fs_3780_2431_2310 * r_2 * h7_m2 - e_4 * fs_38_13_231 * r_4 * h5_m2 - e_4 * fs_2975_143_33 * r_6 * h3_m2 + e_5 * fs_225_4199_12597 * h11_m8 - e_5 * fs_225_4199_39 * h11_m2 + e_5 * fs_46_221_221 * r_2 * h9_m8 + e_5 * fs_506_221_2 * r_2 * h9_m2 - e_5 * fs_7560_46189_2310 * r_4 * h7_m2 + e_5 * fs_152_663_231 * r_6 * h5_m2 + e_5 * fs_595_429_33 * r_8 * h3_m2 - e_6 * fs_18_4199_12597 * r_2 * h11_m8 + e_6 * fs_18_4199_39 * r_2 * h11_m2 - e_6 * fs_2_221_221 * r_4 * h9_m8 - e_6 * fs_22_221_2 * r_4 * h9_m2 + e_6 * fs_240_46189_2310 * r_6 * h7_m2 - e_6 * fs_4_663_231 * r_8 * h5_m2 - e_6 * fs_14_429_33 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m1, ph11_m10, ph11_m9, ph11_m1, ab_2, pc_20, pc_21 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m1 = ph11_m1[k];

        pc_20[k] = - e_0 * fs_14175_64_110 * h1_m1 + e_1 * fs_8505_16_110 * r_2 * h1_m1 + e_2 * fs_855_8_66 * h5_m1 - e_2 * fs_6075_16_110 * r_4 * h1_m1 + e_3 * fs_315_52_770 * h7_m1 - e_3 * fs_855_13_66 * r_2 * h5_m1 + e_3 * fs_225_2_110 * r_6 * h1_m1 - e_4 * fs_1449_442_221 * h9_m9 + e_4 * fs_1449_748_22 * h9_m1 - e_4 * fs_945_442_770 * r_2 * h7_m1 + e_4 * fs_171_13_66 * r_4 * h5_m1 - e_4 * fs_675_44_110 * r_8 * h1_m1 + e_5 * fs_150_4199_41990 * h11_m9 + e_5 * fs_150_4199_15 * h11_m1 + e_5 * fs_138_221_221 * r_2 * h9_m9 - e_5 * fs_69_187_22 * r_2 * h9_m1 + e_5 * fs_945_4199_770 * r_4 * h7_m1 - e_5 * fs_228_221_66 * r_6 * h5_m1 + e_5 * fs_135_143_110 * r_10 * h1_m1 - e_6 * fs_12_4199_41990 * r_2 * h11_m9 - e_6 * fs_12_4199_15 * r_2 * h11_m1 - e_6 * fs_6_221_221 * r_4 * h9_m9 + e_6 * fs_3_187_22 * r_4 * h9_m1 - e_6 * fs_30_4199_770 * r_6 * h7_m1 + e_6 * fs_6_221_66 * r_8 * h5_m1 - e_6 * fs_3_143_110 * r_12 * h1_m1;

        pc_21[k] = e_5 * fs_75_4199_176358 * h11_m10 - e_6 * fs_6_4199_176358 * r_2 * h11_m10;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2, pc_22 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_22[k] = e_0 * fs_14175_64_2 * h1_p1 + e_1 * fs_16065_16_3 * h3_p1 - e_1 * fs_8505_16_2 * r_2 * h1_p1 + e_2 * fs_1425_8_30 * h5_p1 - e_2 * fs_8925_8_3 * r_2 * h3_p1 + e_2 * fs_6075_16_2 * r_4 * h1_p1 + e_3 * fs_11025_286_14 * h7_p1 - e_3 * fs_1425_13_30 * r_2 * h5_p1 + e_3 * fs_8925_22_3 * r_4 * h3_p1 - e_3 * fs_225_2_2 * r_6 * h1_p1 + e_4 * fs_21735_9724_10 * h9_p1 - e_4 * fs_1449_4862_12155 * h9_p9 - e_4 * fs_33075_2431_14 * r_2 * h7_p1 + e_4 * fs_285_13_30 * r_4 * h5_p1 - e_4 * fs_8925_143_3 * r_6 * h3_p1 + e_4 * fs_675_44_2 * r_8 * h1_p1 + e_5 * fs_75_4199_33 * h11_p1 - e_5 * fs_75_4199_92378 * h11_p9 - e_5 * fs_1035_2431_10 * r_2 * h9_p1 + e_5 * fs_138_2431_12155 * r_2 * h9_p9 + e_5 * fs_66150_46189_14 * r_4 * h7_p1 - e_5 * fs_380_221_30 * r_6 * h5_p1 + e_5 * fs_595_143_3 * r_8 * h3_p1 - e_5 * fs_135_143_2 * r_10 * h1_p1 - e_6 * fs_6_4199_33 * r_2 * h11_p1 + e_6 * fs_6_4199_92378 * r_2 * h11_p9 + e_6 * fs_45_2431_10 * r_4 * h9_p1 - e_6 * fs_6_2431_12155 * r_4 * h9_p9 - e_6 * fs_2100_46189_14 * r_6 * h7_p1 + e_6 * fs_10_221_30 * r_8 * h5_p1 - e_6 * fs_14_143_3 * r_10 * h3_p1 + e_6 * fs_3_143_2 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph9_p8, ph11_0, ph11_p8, ab_2, pc_23 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p8 = ph11_p8[k];

        pc_23[k] = - e_0 * fs_14175_16_5 * h1_0 - e_1 * fs_16065_16_5 * h3_0 + e_1 * fs_8505_4_5 * r_2 * h1_0 + e_2 * fs_285_1_5 * h5_0 + e_2 * fs_8925_8_5 * r_2 * h3_0 - e_2 * fs_6075_4_5 * r_4 * h1_0 + e_3 * fs_37485_286_5 * h7_0 - e_3 * fs_2280_13_5 * r_2 * h5_0 - e_3 * fs_8925_22_5 * r_4 * h3_0 + e_3 * fs_450_1_5 * r_6 * h1_0 + e_4 * fs_1449_143_5 * h9_0 + e_4 * fs_1449_2431_2431 * h9_p8 - e_4 * fs_6615_143_5 * r_2 * h7_0 + e_4 * fs_456_13_5 * r_4 * h5_0 + e_4 * fs_8925_143_5 * r_6 * h3_0 - e_4 * fs_675_11_5 * r_8 * h1_0 + e_5 * fs_825_4199_5 * h11_0 - e_5 * fs_75_4199_138567 * h11_p8 - e_5 * fs_276_143_5 * r_2 * h9_0 - e_5 * fs_276_2431_2431 * r_2 * h9_p8 + e_5 * fs_13230_2717_5 * r_4 * h7_0 - e_5 * fs_608_221_5 * r_6 * h5_0 - e_5 * fs_595_143_5 * r_8 * h3_0 + e_5 * fs_540_143_5 * r_10 * h1_0 - e_6 * fs_66_4199_5 * r_2 * h11_0 + e_6 * fs_6_4199_138567 * r_2 * h11_p8 + e_6 * fs_12_143_5 * r_4 * h9_0 + e_6 * fs_12_2431_2431 * r_4 * h9_p8 - e_6 * fs_420_2717_5 * r_6 * h7_0 + e_6 * fs_16_221_5 * r_8 * h5_0 + e_6 * fs_14_143_5 * r_10 * h3_0 - e_6 * fs_12_143_5 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2, pc_24 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_24[k] = - e_0 * fs_42525_64_10 * h1_p1 + e_1 * fs_5355_16_15 * h3_p1 + e_1 * fs_25515_16_10 * r_2 * h1_p1 + e_2 * fs_1235_8_6 * h5_p1 - e_2 * fs_2975_8_15 * r_2 * h3_p1 - e_2 * fs_18225_16_10 * r_4 * h1_p1 - e_3 * fs_945_44_70 * h7_p1 - e_3 * fs_2205_572_4290 * h7_p7 - e_3 * fs_95_1_6 * r_2 * h5_p1 + e_3 * fs_2975_22_15 * r_4 * h3_p1 + e_3 * fs_675_2_10 * r_6 * h1_p1 - e_4 * fs_177261_9724_2 * h9_p1 + e_4 * fs_14007_4862_143 * h9_p7 + e_4 * fs_2835_374_70 * r_2 * h7_p1 + e_4 * fs_6615_4862_4290 * r_2 * h7_p7 + e_4 * fs_19_1_6 * r_4 * h5_p1 - e_4 * fs_2975_143_15 * r_6 * h3_p1 - e_4 * fs_2025_44_10 * r_8 * h1_p1 - e_5 * fs_225_4199_165 * h11_p1 - e_5 * fs_225_4199_14586 * h11_p7 + e_5 * fs_8441_2431_2 * r_2 * h9_p1 - e_5 * fs_1334_2431_143 * r_2 * h9_p7 - e_5 * fs_2835_3553_70 * r_4 * h7_p1 - e_5 * fs_6615_46189_4290 * r_4 * h7_p7 - e_5 * fs_76_51_6 * r_6 * h5_p1 + e_5 * fs_595_429_15 * r_8 * h3_p1 + e_5 * fs_405_143_10 * r_10 * h1_p1 + e_6 * fs_18_4199_165 * r_2 * h11_p1 + e_6 * fs_18_4199_14586 * r_2 * h11_p7 - e_6 * fs_367_2431_2 * r_4 * h9_p1 + e_6 * fs_58_2431_143 * r_4 * h9_p7 + e_6 * fs_90_3553_70 * r_6 * h7_p1 + e_6 * fs_210_46189_4290 * r_6 * h7_p7 + e_6 * fs_2_51_6 * r_8 * h5_p1 - e_6 * fs_14_429_15 * r_10 * h3_p1 - e_6 * fs_9_143_10 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2, pc_25 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_25[k] = e_1 * f_16065_16 * h3_p2 - e_2 * fs_285_1_7 * h5_p2 - e_2 * f_8925_8 * r_2 * h3_p2 + e_3 * fs_9135_1144_70 * h7_p2 - e_3 * fs_315_1144_10010 * h7_p6 + e_3 * fs_2280_13_7 * r_2 * h5_p2 + e_3 * f_8925_22 * r_4 * h3_p2 + e_4 * fs_19803_4862_66 * h9_p2 + e_4 * fs_4347_4862_286 * h9_p6 - e_4 * fs_27405_9724_70 * r_2 * h7_p2 + e_4 * fs_945_9724_10010 * r_2 * h7_p6 - e_4 * fs_456_13_7 * r_4 * h5_p2 - e_4 * f_8925_143 * r_6 * h3_p2 + e_5 * fs_450_4199_143 * h11_p2 - e_5 * fs_150_4199_24310 * h11_p6 - e_5 * fs_1886_2431_66 * r_2 * h9_p2 - e_5 * fs_414_2431_286 * r_2 * h9_p6 + e_5 * fs_27405_92378_70 * r_4 * h7_p2 - e_5 * fs_945_92378_10010 * r_4 * h7_p6 + e_5 * fs_608_221_7 * r_6 * h5_p2 + e_5 * f_595_143 * r_8 * h3_p2 - e_6 * fs_36_4199_143 * r_2 * h11_p2 + e_6 * fs_12_4199_24310 * r_2 * h11_p6 + e_6 * fs_82_2431_66 * r_4 * h9_p2 + e_6 * fs_18_2431_286 * r_4 * h9_p6 - e_6 * fs_435_46189_70 * r_6 * h7_p2 + e_6 * fs_15_46189_10010 * r_6 * h7_p6 - e_6 * fs_16_221_7 * r_8 * h5_p2 - e_6 * f_14_143 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2, pc_26 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_26[k] = - e_1 * fs_5355_16_42 * h3_p3 + e_2 * fs_1235_8_6 * h5_p3 - e_2 * fs_1425_8_30 * h5_p5 + e_2 * fs_2975_8_42 * r_2 * h3_p3 + e_3 * fs_1575_44_5 * h7_p3 + e_3 * fs_11655_572_55 * h7_p5 - e_3 * fs_95_1_6 * r_2 * h5_p3 + e_3 * fs_1425_13_30 * r_2 * h5_p5 - e_3 * fs_2975_22_42 * r_4 * h3_p3 - e_4 * fs_71001_9724_22 * h9_p3 - e_4 * fs_483_9724_30030 * h9_p5 - e_4 * fs_4725_374_5 * r_2 * h7_p3 - e_4 * fs_34965_4862_55 * r_2 * h7_p5 + e_4 * fs_19_1_6 * r_4 * h5_p3 - e_4 * fs_285_13_30 * r_4 * h5_p5 + e_4 * fs_2975_143_42 * r_6 * h3_p3 - e_5 * fs_525_4199_286 * h11_p3 - e_5 * fs_150_4199_15015 * h11_p5 + e_5 * fs_3381_2431_22 * r_2 * h9_p3 + e_5 * fs_23_2431_30030 * r_2 * h9_p5 + e_5 * fs_4725_3553_5 * r_4 * h7_p3 + e_5 * fs_34965_46189_55 * r_4 * h7_p5 - e_5 * fs_76_51_6 * r_6 * h5_p3 + e_5 * fs_380_221_30 * r_6 * h5_p5 - e_5 * fs_595_429_42 * r_8 * h3_p3 + e_6 * fs_42_4199_286 * r_2 * h11_p3 + e_6 * fs_12_4199_15015 * r_2 * h11_p5 - e_6 * fs_147_2431_22 * r_4 * h9_p3 - e_6 * fs_1_2431_30030 * r_4 * h9_p5 - e_6 * fs_150_3553_5 * r_6 * h7_p3 - e_6 * fs_1110_46189_55 * r_6 * h7_p5 + e_6 * fs_2_51_6 * r_8 * h5_p3 - e_6 * fs_10_221_30 * r_8 * h5_p5 + e_6 * fs_14_429_42 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2, pc_27 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_27[k] = e_2 * fs_285_1_5 * h5_m4 - e_3 * fs_6300_143_33 * h7_m4 - e_3 * fs_2280_13_5 * r_2 * h5_m4 + e_4 * fs_3381_2431_715 * h9_m4 + e_4 * fs_37800_2431_33 * r_2 * h7_m4 + e_4 * fs_456_13_5 * r_4 * h5_m4 + e_5 * fs_1575_4199_143 * h11_m4 - e_5 * fs_644_2431_715 * r_2 * h9_m4 - e_5 * fs_75600_46189_33 * r_4 * h7_m4 - e_5 * fs_608_221_5 * r_6 * h5_m4 - e_6 * fs_126_4199_143 * r_2 * h11_m4 + e_6 * fs_28_2431_715 * r_4 * h9_m4 + e_6 * fs_2400_46189_33 * r_6 * h7_m4 + e_6 * fs_16_221_5 * r_8 * h5_m4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2, pc_28 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_28[k] = - e_1 * fs_5355_16_42 * h3_m3 + e_2 * fs_1425_8_30 * h5_m5 + e_2 * fs_1235_8_6 * h5_m3 + e_2 * fs_2975_8_42 * r_2 * h3_m3 - e_3 * fs_11655_572_55 * h7_m5 + e_3 * fs_1575_44_5 * h7_m3 - e_3 * fs_1425_13_30 * r_2 * h5_m5 - e_3 * fs_95_1_6 * r_2 * h5_m3 - e_3 * fs_2975_22_42 * r_4 * h3_m3 + e_4 * fs_483_9724_30030 * h9_m5 - e_4 * fs_71001_9724_22 * h9_m3 + e_4 * fs_34965_4862_55 * r_2 * h7_m5 - e_4 * fs_4725_374_5 * r_2 * h7_m3 + e_4 * fs_285_13_30 * r_4 * h5_m5 + e_4 * fs_19_1_6 * r_4 * h5_m3 + e_4 * fs_2975_143_42 * r_6 * h3_m3 + e_5 * fs_150_4199_15015 * h11_m5 - e_5 * fs_525_4199_286 * h11_m3 - e_5 * fs_23_2431_30030 * r_2 * h9_m5 + e_5 * fs_3381_2431_22 * r_2 * h9_m3 - e_5 * fs_34965_46189_55 * r_4 * h7_m5 + e_5 * fs_4725_3553_5 * r_4 * h7_m3 - e_5 * fs_380_221_30 * r_6 * h5_m5 - e_5 * fs_76_51_6 * r_6 * h5_m3 - e_5 * fs_595_429_42 * r_8 * h3_m3 - e_6 * fs_12_4199_15015 * r_2 * h11_m5 + e_6 * fs_42_4199_286 * r_2 * h11_m3 + e_6 * fs_1_2431_30030 * r_4 * h9_m5 - e_6 * fs_147_2431_22 * r_4 * h9_m3 + e_6 * fs_1110_46189_55 * r_6 * h7_m5 - e_6 * fs_150_3553_5 * r_6 * h7_m3 + e_6 * fs_10_221_30 * r_8 * h5_m5 + e_6 * fs_2_51_6 * r_8 * h5_m3 + e_6 * fs_14_429_42 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2, pc_29 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_29[k] = e_1 * f_16065_16 * h3_m2 - e_2 * fs_285_1_7 * h5_m2 - e_2 * f_8925_8 * r_2 * h3_m2 + e_3 * fs_315_1144_10010 * h7_m6 + e_3 * fs_9135_1144_70 * h7_m2 + e_3 * fs_2280_13_7 * r_2 * h5_m2 + e_3 * f_8925_22 * r_4 * h3_m2 - e_4 * fs_4347_4862_286 * h9_m6 + e_4 * fs_19803_4862_66 * h9_m2 - e_4 * fs_945_9724_10010 * r_2 * h7_m6 - e_4 * fs_27405_9724_70 * r_2 * h7_m2 - e_4 * fs_456_13_7 * r_4 * h5_m2 - e_4 * f_8925_143 * r_6 * h3_m2 + e_5 * fs_150_4199_24310 * h11_m6 + e_5 * fs_450_4199_143 * h11_m2 + e_5 * fs_414_2431_286 * r_2 * h9_m6 - e_5 * fs_1886_2431_66 * r_2 * h9_m2 + e_5 * fs_945_92378_10010 * r_4 * h7_m6 + e_5 * fs_27405_92378_70 * r_4 * h7_m2 + e_5 * fs_608_221_7 * r_6 * h5_m2 + e_5 * f_595_143 * r_8 * h3_m2 - e_6 * fs_12_4199_24310 * r_2 * h11_m6 - e_6 * fs_36_4199_143 * r_2 * h11_m2 - e_6 * fs_18_2431_286 * r_4 * h9_m6 + e_6 * fs_82_2431_66 * r_4 * h9_m2 - e_6 * fs_15_46189_10010 * r_6 * h7_m6 - e_6 * fs_435_46189_70 * r_6 * h7_m2 - e_6 * fs_16_221_7 * r_8 * h5_m2 - e_6 * f_14_143 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2, pc_30 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_30[k] = - e_0 * fs_42525_64_10 * h1_m1 + e_1 * fs_5355_16_15 * h3_m1 + e_1 * fs_25515_16_10 * r_2 * h1_m1 + e_2 * fs_1235_8_6 * h5_m1 - e_2 * fs_2975_8_15 * r_2 * h3_m1 - e_2 * fs_18225_16_10 * r_4 * h1_m1 + e_3 * fs_2205_572_4290 * h7_m7 - e_3 * fs_945_44_70 * h7_m1 - e_3 * fs_95_1_6 * r_2 * h5_m1 + e_3 * fs_2975_22_15 * r_4 * h3_m1 + e_3 * fs_675_2_10 * r_6 * h1_m1 - e_4 * fs_14007_4862_143 * h9_m7 - e_4 * fs_177261_9724_2 * h9_m1 - e_4 * fs_6615_4862_4290 * r_2 * h7_m7 + e_4 * fs_2835_374_70 * r_2 * h7_m1 + e_4 * fs_19_1_6 * r_4 * h5_m1 - e_4 * fs_2975_143_15 * r_6 * h3_m1 - e_4 * fs_2025_44_10 * r_8 * h1_m1 + e_5 * fs_225_4199_14586 * h11_m7 - e_5 * fs_225_4199_165 * h11_m1 + e_5 * fs_1334_2431_143 * r_2 * h9_m7 + e_5 * fs_8441_2431_2 * r_2 * h9_m1 + e_5 * fs_6615_46189_4290 * r_4 * h7_m7 - e_5 * fs_2835_3553_70 * r_4 * h7_m1 - e_5 * fs_76_51_6 * r_6 * h5_m1 + e_5 * fs_595_429_15 * r_8 * h3_m1 + e_5 * fs_405_143_10 * r_10 * h1_m1 - e_6 * fs_18_4199_14586 * r_2 * h11_m7 + e_6 * fs_18_4199_165 * r_2 * h11_m1 - e_6 * fs_58_2431_143 * r_4 * h9_m7 - e_6 * fs_367_2431_2 * r_4 * h9_m1 - e_6 * fs_210_46189_4290 * r_6 * h7_m7 + e_6 * fs_90_3553_70 * r_6 * h7_m1 + e_6 * fs_2_51_6 * r_8 * h5_m1 - e_6 * fs_14_429_15 * r_10 * h3_m1 - e_6 * fs_9_143_10 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m8, ph9_m1, ph11_m9, ph11_m8, ph11_m1, ab_2, pc_31, pc_32 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m1 = ph11_m1[k];

        pc_31[k] = - e_4 * fs_1449_2431_2431 * h9_m8 + e_5 * fs_75_4199_138567 * h11_m8 + e_5 * fs_276_2431_2431 * r_2 * h9_m8 - e_6 * fs_6_4199_138567 * r_2 * h11_m8 - e_6 * fs_12_2431_2431 * r_4 * h9_m8;

        pc_32[k] = - e_0 * fs_14175_64_2 * h1_m1 - e_1 * fs_16065_16_3 * h3_m1 + e_1 * fs_8505_16_2 * r_2 * h1_m1 - e_2 * fs_1425_8_30 * h5_m1 + e_2 * fs_8925_8_3 * r_2 * h3_m1 - e_2 * fs_6075_16_2 * r_4 * h1_m1 - e_3 * fs_11025_286_14 * h7_m1 + e_3 * fs_1425_13_30 * r_2 * h5_m1 - e_3 * fs_8925_22_3 * r_4 * h3_m1 + e_3 * fs_225_2_2 * r_6 * h1_m1 + e_4 * fs_1449_4862_12155 * h9_m9 - e_4 * fs_21735_9724_10 * h9_m1 + e_4 * fs_33075_2431_14 * r_2 * h7_m1 - e_4 * fs_285_13_30 * r_4 * h5_m1 + e_4 * fs_8925_143_3 * r_6 * h3_m1 - e_4 * fs_675_44_2 * r_8 * h1_m1 + e_5 * fs_75_4199_92378 * h11_m9 - e_5 * fs_75_4199_33 * h11_m1 - e_5 * fs_138_2431_12155 * r_2 * h9_m9 + e_5 * fs_1035_2431_10 * r_2 * h9_m1 - e_5 * fs_66150_46189_14 * r_4 * h7_m1 + e_5 * fs_380_221_30 * r_6 * h5_m1 - e_5 * fs_595_143_3 * r_8 * h3_m1 + e_5 * fs_135_143_2 * r_10 * h1_m1 - e_6 * fs_6_4199_92378 * r_2 * h11_m9 + e_6 * fs_6_4199_33 * r_2 * h11_m1 + e_6 * fs_6_2431_12155 * r_4 * h9_m9 - e_6 * fs_45_2431_10 * r_4 * h9_m1 + e_6 * fs_2100_46189_14 * r_6 * h7_m1 - e_6 * fs_10_221_30 * r_8 * h5_m1 + e_6 * fs_14_143_3 * r_10 * h3_m1 - e_6 * fs_3_143_2 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2, pc_33 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_33[k] = - e_1 * f_16065_16 * h3_p2 - e_2 * fs_1425_4_7 * h5_p2 + e_2 * f_8925_8 * r_2 * h3_p2 - e_3 * fs_6615_286_70 * h7_p2 + e_3 * fs_2850_13_7 * r_2 * h5_p2 - e_3 * f_8925_22 * r_4 * h3_p2 - e_4 * fs_7245_4862_66 * h9_p2 - e_4 * fs_2415_4862_7293 * h9_p8 + e_4 * fs_19845_2431_70 * r_2 * h7_p2 - e_4 * fs_570_13_7 * r_4 * h5_p2 + e_4 * f_8925_143 * r_6 * h3_p2 - e_5 * fs_75_4199_143 * h11_p2 - e_5 * fs_75_4199_46189 * h11_p8 + e_5 * fs_690_2431_66 * r_2 * h9_p2 + e_5 * fs_230_2431_7293 * r_2 * h9_p8 - e_5 * fs_39690_46189_70 * r_4 * h7_p2 + e_5 * fs_760_221_7 * r_6 * h5_p2 - e_5 * f_595_143 * r_8 * h3_p2 + e_6 * fs_6_4199_143 * r_2 * h11_p2 + e_6 * fs_6_4199_46189 * r_2 * h11_p8 - e_6 * fs_30_2431_66 * r_4 * h9_p2 - e_6 * fs_10_2431_7293 * r_4 * h9_p8 + e_6 * fs_1260_46189_70 * r_6 * h7_p2 - e_6 * fs_20_221_7 * r_8 * h5_p2 + e_6 * f_14_143 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2, pc_34 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_34[k] = e_0 * fs_14175_64_6 * h1_p1 + e_1 * f_16065_8 * h3_p1 - e_1 * fs_8505_16_6 * r_2 * h1_p1 + e_2 * fs_285_8_10 * h5_p1 - e_2 * f_8925_4 * r_2 * h3_p1 + e_2 * fs_6075_16_6 * r_4 * h1_p1 - e_3 * fs_4410_143_42 * h7_p1 + e_3 * fs_6615_572_286 * h7_p7 - e_3 * fs_285_13_10 * r_2 * h5_p1 + e_3 * f_8925_11 * r_4 * h3_p1 - e_3 * fs_225_2_6 * r_6 * h1_p1 - e_4 * fs_42021_9724_30 * h9_p1 + e_4 * fs_483_4862_2145 * h9_p7 + e_4 * fs_26460_2431_42 * r_2 * h7_p1 - e_4 * fs_19845_4862_286 * r_2 * h7_p7 + e_4 * fs_57_13_10 * r_4 * h5_p1 - e_4 * f_17850_143 * r_6 * h3_p1 + e_4 * fs_675_44_6 * r_8 * h1_p1 - e_5 * fs_750_4199_11 * h11_p1 - e_5 * fs_150_4199_24310 * h11_p7 + e_5 * fs_2001_2431_30 * r_2 * h9_p1 - e_5 * fs_46_2431_2145 * r_2 * h9_p7 - e_5 * fs_52920_46189_42 * r_4 * h7_p1 + e_5 * fs_19845_46189_286 * r_4 * h7_p7 - e_5 * fs_76_221_10 * r_6 * h5_p1 + e_5 * f_1190_143 * r_8 * h3_p1 - e_5 * fs_135_143_6 * r_10 * h1_p1 + e_6 * fs_60_4199_11 * r_2 * h11_p1 + e_6 * fs_12_4199_24310 * r_2 * h11_p7 - e_6 * fs_87_2431_30 * r_4 * h9_p1 + e_6 * fs_2_2431_2145 * r_4 * h9_p7 + e_6 * fs_1680_46189_42 * r_6 * h7_p1 - e_6 * fs_630_46189_286 * r_6 * h7_p7 + e_6 * fs_2_221_10 * r_8 * h5_p1 - e_6 * f_28_143 * r_10 * h3_p1 + e_6 * fs_3_143_6 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ph11_0, ph11_p6, ab_2, pc_35 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p6 = ph11_p6[k];

        pc_35[k] = - e_0 * fs_42525_32_3 * h1_0 - e_1 * fs_5355_16_3 * h3_0 + e_1 * fs_25515_8_3 * r_2 * h1_0 + e_2 * fs_2375_4_3 * h5_0 + e_2 * fs_2975_8_3 * r_2 * h3_0 - e_2 * fs_18225_8_3 * r_4 * h1_0 - e_3 * fs_19845_286_3 * h7_0 - e_3 * fs_945_286_2002 * h7_p6 - e_3 * fs_4750_13_3 * r_2 * h5_0 - e_3 * fs_2975_22_3 * r_4 * h3_0 + e_3 * fs_675_1_3 * r_6 * h1_0 - e_4 * fs_127995_4862_3 * h9_0 + e_4 * fs_1932_2431_1430 * h9_p6 + e_4 * fs_59535_2431_3 * r_2 * h7_0 + e_4 * fs_2835_2431_2002 * r_2 * h7_p6 + e_4 * fs_950_13_3 * r_4 * h5_0 + e_4 * fs_2975_143_3 * r_6 * h3_0 - e_4 * fs_2025_22_3 * r_8 * h1_0 - e_5 * fs_4125_4199_3 * h11_0 - e_5 * fs_375_4199_4862 * h11_p6 + e_5 * fs_12190_2431_3 * r_2 * h9_0 - e_5 * fs_368_2431_1430 * r_2 * h9_p6 - e_5 * fs_119070_46189_3 * r_4 * h7_0 - e_5 * fs_5670_46189_2002 * r_4 * h7_p6 - e_5 * fs_3800_663_3 * r_6 * h5_0 - e_5 * fs_595_429_3 * r_8 * h3_0 + e_5 * fs_810_143_3 * r_10 * h1_0 + e_6 * fs_330_4199_3 * r_2 * h11_0 + e_6 * fs_30_4199_4862 * r_2 * h11_p6 - e_6 * fs_530_2431_3 * r_4 * h9_0 + e_6 * fs_16_2431_1430 * r_4 * h9_p6 + e_6 * fs_3780_46189_3 * r_6 * h7_0 + e_6 * fs_180_46189_2002 * r_6 * h7_p6 + e_6 * fs_100_663_3 * r_8 * h5_0 + e_6 * fs_14_429_3 * r_10 * h3_0 - e_6 * fs_18_143_3 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2, pc_36 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p5 = ph11_p5[k];

        pc_36[k] = - e_0 * fs_42525_32_2 * h1_p1 + e_1 * fs_16065_16_3 * h3_p1 + e_1 * fs_25515_8_2 * r_2 * h1_p1 - e_2 * fs_285_4_30 * h5_p1 + e_2 * fs_1425_4_7 * h5_p5 - e_2 * fs_8925_8_3 * r_2 * h3_p1 - e_2 * fs_18225_8_2 * r_4 * h1_p1 - e_3 * fs_21735_1144_14 * h7_p1 - e_3 * fs_9135_1144_462 * h7_p5 + e_3 * fs_570_13_30 * r_2 * h5_p1 - e_3 * fs_2850_13_7 * r_2 * h5_p5 + e_3 * fs_8925_22_3 * r_4 * h3_p1 + e_3 * fs_675_1_2 * r_6 * h1_p1 + e_4 * fs_51681_4862_10 * h9_p1 + e_4 * fs_12075_4862_143 * h9_p5 + e_4 * fs_65205_9724_14 * r_2 * h7_p1 + e_4 * fs_27405_9724_462 * r_2 * h7_p5 - e_4 * fs_114_13_30 * r_4 * h5_p1 + e_4 * fs_570_13_7 * r_4 * h5_p5 - e_4 * fs_8925_143_3 * r_6 * h3_p1 - e_4 * fs_2025_22_2 * r_8 * h1_p1 + e_5 * fs_1500_4199_33 * h11_p1 - e_5 * fs_1500_4199_286 * h11_p5 - e_5 * fs_4922_2431_10 * r_2 * h9_p1 - e_5 * fs_1150_2431_143 * r_2 * h9_p5 - e_5 * fs_65205_92378_14 * r_4 * h7_p1 - e_5 * fs_27405_92378_462 * r_4 * h7_p5 + e_5 * fs_152_221_30 * r_6 * h5_p1 - e_5 * fs_760_221_7 * r_6 * h5_p5 + e_5 * fs_595_143_3 * r_8 * h3_p1 + e_5 * fs_810_143_2 * r_10 * h1_p1 - e_6 * fs_120_4199_33 * r_2 * h11_p1 + e_6 * fs_120_4199_286 * r_2 * h11_p5 + e_6 * fs_214_2431_10 * r_4 * h9_p1 + e_6 * fs_50_2431_143 * r_4 * h9_p5 + e_6 * fs_1035_46189_14 * r_6 * h7_p1 + e_6 * fs_435_46189_462 * r_6 * h7_p5 - e_6 * fs_4_221_30 * r_8 * h5_p1 + e_6 * fs_20_221_7 * r_8 * h5_p5 - e_6 * fs_14_143_3 * r_10 * h3_p1 - e_6 * fs_18_143_2 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2, pc_37 : simd::cache_line_size())
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

        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_37[k] = - e_2 * fs_285_4_30 * h5_p2 - e_2 * fs_285_8_10 * h5_p4 + e_3 * fs_1260_13_3 * h7_p2 - e_3 * fs_945_572_66 * h7_p4 + e_3 * fs_570_13_30 * r_2 * h5_p2 + e_3 * fs_285_13_10 * r_2 * h5_p4 - e_4 * fs_483_374_385 * h9_p2 + e_4 * fs_3381_9724_1430 * h9_p4 - e_4 * fs_7560_221_3 * r_2 * h7_p2 + e_4 * fs_2835_4862_66 * r_2 * h7_p4 - e_4 * fs_114_13_30 * r_4 * h5_p2 - e_4 * fs_57_13_10 * r_4 * h5_p4 - e_5 * fs_75_4199_30030 * h11_p2 - e_5 * fs_2625_8398_286 * h11_p4 + e_5 * fs_46_187_385 * r_2 * h9_p2 - e_5 * fs_161_2431_1430 * r_2 * h9_p4 + e_5 * fs_15120_4199_3 * r_4 * h7_p2 - e_5 * fs_2835_46189_66 * r_4 * h7_p4 + e_5 * fs_152_221_30 * r_6 * h5_p2 + e_5 * fs_76_221_10 * r_6 * h5_p4 + e_6 * fs_6_4199_30030 * r_2 * h11_p2 + e_6 * fs_105_4199_286 * r_2 * h11_p4 - e_6 * fs_2_187_385 * r_4 * h9_p2 + e_6 * fs_7_2431_1430 * r_4 * h9_p4 - e_6 * fs_480_4199_3 * r_6 * h7_p2 + e_6 * fs_90_46189_66 * r_6 * h7_p4 - e_6 * fs_4_221_30 * r_8 * h5_p2 - e_6 * fs_2_221_10 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m3, ph7_m3, ph9_m3, ph11_m3, ab_2, pc_38 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m3 = ph11_m3[k];

        pc_38[k] = - e_1 * fs_5355_8_21 * h3_m3 + e_2 * fs_2375_4_3 * h5_m3 + e_2 * fs_2975_4_21 * r_2 * h3_m3 - e_3 * fs_34965_572_10 * h7_m3 - e_3 * fs_4750_13_3 * r_2 * h5_m3 - e_3 * fs_2975_11_21 * r_4 * h3_m3 + e_4 * fs_16905_4862_11 * h9_m3 + e_4 * fs_104895_4862_10 * r_2 * h7_m3 + e_4 * fs_950_13_3 * r_4 * h5_m3 + e_4 * fs_5950_143_21 * r_6 * h3_m3 + e_5 * fs_2100_4199_143 * h11_m3 - e_5 * fs_1610_2431_11 * r_2 * h9_m3 - e_5 * fs_104895_46189_10 * r_4 * h7_m3 - e_5 * fs_3800_663_3 * r_6 * h5_m3 - e_5 * fs_1190_429_21 * r_8 * h3_m3 - e_6 * fs_168_4199_143 * r_2 * h11_m3 + e_6 * fs_70_2431_11 * r_4 * h9_m3 + e_6 * fs_3330_46189_10 * r_6 * h7_m3 + e_6 * fs_100_663_3 * r_8 * h5_m3 + e_6 * fs_28_429_21 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2, pc_39 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_39[k] = e_2 * fs_285_8_10 * h5_m4 - e_2 * fs_285_4_30 * h5_m2 + e_3 * fs_945_572_66 * h7_m4 + e_3 * fs_1260_13_3 * h7_m2 - e_3 * fs_285_13_10 * r_2 * h5_m4 + e_3 * fs_570_13_30 * r_2 * h5_m2 - e_4 * fs_3381_9724_1430 * h9_m4 - e_4 * fs_483_374_385 * h9_m2 - e_4 * fs_2835_4862_66 * r_2 * h7_m4 - e_4 * fs_7560_221_3 * r_2 * h7_m2 + e_4 * fs_57_13_10 * r_4 * h5_m4 - e_4 * fs_114_13_30 * r_4 * h5_m2 + e_5 * fs_2625_8398_286 * h11_m4 - e_5 * fs_75_4199_30030 * h11_m2 + e_5 * fs_161_2431_1430 * r_2 * h9_m4 + e_5 * fs_46_187_385 * r_2 * h9_m2 + e_5 * fs_2835_46189_66 * r_4 * h7_m4 + e_5 * fs_15120_4199_3 * r_4 * h7_m2 - e_5 * fs_76_221_10 * r_6 * h5_m4 + e_5 * fs_152_221_30 * r_6 * h5_m2 - e_6 * fs_105_4199_286 * r_2 * h11_m4 + e_6 * fs_6_4199_30030 * r_2 * h11_m2 - e_6 * fs_7_2431_1430 * r_4 * h9_m4 - e_6 * fs_2_187_385 * r_4 * h9_m2 - e_6 * fs_90_46189_66 * r_6 * h7_m4 - e_6 * fs_480_4199_3 * r_6 * h7_m2 + e_6 * fs_2_221_10 * r_8 * h5_m4 - e_6 * fs_4_221_30 * r_8 * h5_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2, pc_40 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m1 = ph11_m1[k];

        pc_40[k] = - e_0 * fs_42525_32_2 * h1_m1 + e_1 * fs_16065_16_3 * h3_m1 + e_1 * fs_25515_8_2 * r_2 * h1_m1 - e_2 * fs_1425_4_7 * h5_m5 - e_2 * fs_285_4_30 * h5_m1 - e_2 * fs_8925_8_3 * r_2 * h3_m1 - e_2 * fs_18225_8_2 * r_4 * h1_m1 + e_3 * fs_9135_1144_462 * h7_m5 - e_3 * fs_21735_1144_14 * h7_m1 + e_3 * fs_2850_13_7 * r_2 * h5_m5 + e_3 * fs_570_13_30 * r_2 * h5_m1 + e_3 * fs_8925_22_3 * r_4 * h3_m1 + e_3 * fs_675_1_2 * r_6 * h1_m1 - e_4 * fs_12075_4862_143 * h9_m5 + e_4 * fs_51681_4862_10 * h9_m1 - e_4 * fs_27405_9724_462 * r_2 * h7_m5 + e_4 * fs_65205_9724_14 * r_2 * h7_m1 - e_4 * fs_570_13_7 * r_4 * h5_m5 - e_4 * fs_114_13_30 * r_4 * h5_m1 - e_4 * fs_8925_143_3 * r_6 * h3_m1 - e_4 * fs_2025_22_2 * r_8 * h1_m1 + e_5 * fs_1500_4199_286 * h11_m5 + e_5 * fs_1500_4199_33 * h11_m1 + e_5 * fs_1150_2431_143 * r_2 * h9_m5 - e_5 * fs_4922_2431_10 * r_2 * h9_m1 + e_5 * fs_27405_92378_462 * r_4 * h7_m5 - e_5 * fs_65205_92378_14 * r_4 * h7_m1 + e_5 * fs_760_221_7 * r_6 * h5_m5 + e_5 * fs_152_221_30 * r_6 * h5_m1 + e_5 * fs_595_143_3 * r_8 * h3_m1 + e_5 * fs_810_143_2 * r_10 * h1_m1 - e_6 * fs_120_4199_286 * r_2 * h11_m5 - e_6 * fs_120_4199_33 * r_2 * h11_m1 - e_6 * fs_50_2431_143 * r_4 * h9_m5 + e_6 * fs_214_2431_10 * r_4 * h9_m1 - e_6 * fs_435_46189_462 * r_6 * h7_m5 + e_6 * fs_1035_46189_14 * r_6 * h7_m1 - e_6 * fs_20_221_7 * r_8 * h5_m5 - e_6 * fs_4_221_30 * r_8 * h5_m1 - e_6 * fs_14_143_3 * r_10 * h3_m1 - e_6 * fs_18_143_2 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, ph7_m6, ph9_m6, ph11_m6, ab_2, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h7_m6 = ph7_m6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h11_m6 = ph11_m6[k];

        pc_41[k] = e_3 * fs_945_286_2002 * h7_m6 - e_4 * fs_1932_2431_1430 * h9_m6 - e_4 * fs_2835_2431_2002 * r_2 * h7_m6 + e_5 * fs_375_4199_4862 * h11_m6 + e_5 * fs_368_2431_1430 * r_2 * h9_m6 + e_5 * fs_5670_46189_2002 * r_4 * h7_m6 - e_6 * fs_30_4199_4862 * r_2 * h11_m6 - e_6 * fs_16_2431_1430 * r_4 * h9_m6 - e_6 * fs_180_46189_2002 * r_6 * h7_m6;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2, pc_42 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_42[k] = - e_0 * fs_14175_64_6 * h1_m1 - e_1 * f_16065_8 * h3_m1 + e_1 * fs_8505_16_6 * r_2 * h1_m1 - e_2 * fs_285_8_10 * h5_m1 + e_2 * f_8925_4 * r_2 * h3_m1 - e_2 * fs_6075_16_6 * r_4 * h1_m1 - e_3 * fs_6615_572_286 * h7_m7 + e_3 * fs_4410_143_42 * h7_m1 + e_3 * fs_285_13_10 * r_2 * h5_m1 - e_3 * f_8925_11 * r_4 * h3_m1 + e_3 * fs_225_2_6 * r_6 * h1_m1 - e_4 * fs_483_4862_2145 * h9_m7 + e_4 * fs_42021_9724_30 * h9_m1 + e_4 * fs_19845_4862_286 * r_2 * h7_m7 - e_4 * fs_26460_2431_42 * r_2 * h7_m1 - e_4 * fs_57_13_10 * r_4 * h5_m1 + e_4 * f_17850_143 * r_6 * h3_m1 - e_4 * fs_675_44_6 * r_8 * h1_m1 + e_5 * fs_150_4199_24310 * h11_m7 + e_5 * fs_750_4199_11 * h11_m1 + e_5 * fs_46_2431_2145 * r_2 * h9_m7 - e_5 * fs_2001_2431_30 * r_2 * h9_m1 - e_5 * fs_19845_46189_286 * r_4 * h7_m7 + e_5 * fs_52920_46189_42 * r_4 * h7_m1 + e_5 * fs_76_221_10 * r_6 * h5_m1 - e_5 * f_1190_143 * r_8 * h3_m1 + e_5 * fs_135_143_6 * r_10 * h1_m1 - e_6 * fs_12_4199_24310 * r_2 * h11_m7 - e_6 * fs_60_4199_11 * r_2 * h11_m1 - e_6 * fs_2_2431_2145 * r_4 * h9_m7 + e_6 * fs_87_2431_30 * r_4 * h9_m1 + e_6 * fs_630_46189_286 * r_6 * h7_m7 - e_6 * fs_1680_46189_42 * r_6 * h7_m1 - e_6 * fs_2_221_10 * r_8 * h5_m1 + e_6 * f_28_143 * r_10 * h3_m1 - e_6 * fs_3_143_6 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2, pc_43 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_43[k] = e_1 * f_16065_16 * h3_m2 + e_2 * fs_1425_4_7 * h5_m2 - e_2 * f_8925_8 * r_2 * h3_m2 + e_3 * fs_6615_286_70 * h7_m2 - e_3 * fs_2850_13_7 * r_2 * h5_m2 + e_3 * f_8925_22 * r_4 * h3_m2 + e_4 * fs_2415_4862_7293 * h9_m8 + e_4 * fs_7245_4862_66 * h9_m2 - e_4 * fs_19845_2431_70 * r_2 * h7_m2 + e_4 * fs_570_13_7 * r_4 * h5_m2 - e_4 * f_8925_143 * r_6 * h3_m2 + e_5 * fs_75_4199_46189 * h11_m8 + e_5 * fs_75_4199_143 * h11_m2 - e_5 * fs_230_2431_7293 * r_2 * h9_m8 - e_5 * fs_690_2431_66 * r_2 * h9_m2 + e_5 * fs_39690_46189_70 * r_4 * h7_m2 - e_5 * fs_760_221_7 * r_6 * h5_m2 + e_5 * f_595_143 * r_8 * h3_m2 - e_6 * fs_6_4199_46189 * r_2 * h11_m8 - e_6 * fs_6_4199_143 * r_2 * h11_m2 + e_6 * fs_10_2431_7293 * r_4 * h9_m8 + e_6 * fs_30_2431_66 * r_4 * h9_m2 - e_6 * fs_1260_46189_70 * r_6 * h7_m2 + e_6 * fs_20_221_7 * r_8 * h5_m2 - e_6 * f_14_143 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2, pc_44 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_44[k] = e_1 * fs_5355_32_6 * h3_p3 + e_2 * fs_475_4_42 * h5_p3 - e_2 * fs_2975_16_6 * r_2 * h3_p3 + e_3 * fs_11025_286_35 * h7_p3 - e_3 * fs_2205_572_715 * h7_p7 - e_3 * fs_950_13_42 * r_2 * h5_p3 + e_3 * fs_2975_44_6 * r_4 * h3_p3 + e_4 * fs_7245_4862_154 * h9_p3 - e_4 * fs_7245_4862_858 * h9_p7 - e_4 * fs_33075_2431_35 * r_2 * h7_p3 + e_4 * fs_6615_4862_715 * r_2 * h7_p7 + e_4 * fs_190_13_42 * r_4 * h5_p3 - e_4 * fs_2975_286_6 * r_6 * h3_p3 + e_5 * fs_75_8398_2002 * h11_p3 - e_5 * fs_225_4199_2431 * h11_p7 - e_5 * fs_690_2431_154 * r_2 * h9_p3 + e_5 * fs_690_2431_858 * r_2 * h9_p7 + e_5 * fs_66150_46189_35 * r_4 * h7_p3 - e_5 * fs_6615_46189_715 * r_4 * h7_p7 - e_5 * fs_760_663_42 * r_6 * h5_p3 + e_5 * fs_595_858_6 * r_8 * h3_p3 - e_6 * fs_3_4199_2002 * r_2 * h11_p3 + e_6 * fs_18_4199_2431 * r_2 * h11_p7 + e_6 * fs_30_2431_154 * r_4 * h9_p3 - e_6 * fs_30_2431_858 * r_4 * h9_p7 - e_6 * fs_2100_46189_35 * r_6 * h7_p3 + e_6 * fs_210_46189_715 * r_6 * h7_p7 + e_6 * fs_20_663_42 * r_8 * h5_p3 - e_6 * fs_7_429_6 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2, pc_45 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_45[k] = - e_1 * fs_16065_32_10 * h3_p2 - e_2 * fs_285_4_70 * h5_p2 + e_2 * fs_8925_16_10 * r_2 * h3_p2 + e_3 * fs_15435_286_7 * h7_p2 + e_3 * fs_315_44_1001 * h7_p6 + e_3 * fs_570_13_70 * r_2 * h5_p2 - e_3 * fs_8925_44_10 * r_4 * h3_p2 + e_4 * fs_5796_2431_165 * h9_p2 - e_4 * fs_1449_2431_715 * h9_p6 - e_4 * fs_46305_2431_7 * r_2 * h7_p2 - e_4 * fs_945_374_1001 * r_2 * h7_p6 - e_4 * fs_114_13_70 * r_4 * h5_p2 + e_4 * fs_8925_286_10 * r_6 * h3_p2 + e_5 * fs_225_8398_1430 * h11_p2 - e_5 * fs_375_4199_2431 * h11_p6 - e_5 * fs_1104_2431_165 * r_2 * h9_p2 + e_5 * fs_276_2431_715 * r_2 * h9_p6 + e_5 * fs_92610_46189_7 * r_4 * h7_p2 + e_5 * fs_945_3553_1001 * r_4 * h7_p6 + e_5 * fs_152_221_70 * r_6 * h5_p2 - e_5 * fs_595_286_10 * r_8 * h3_p2 - e_6 * fs_9_4199_1430 * r_2 * h11_p2 + e_6 * fs_30_4199_2431 * r_2 * h11_p6 + e_6 * fs_48_2431_165 * r_4 * h9_p2 - e_6 * fs_12_2431_715 * r_4 * h9_p6 - e_6 * fs_2940_46189_7 * r_6 * h7_p2 - e_6 * fs_30_3553_1001 * r_6 * h7_p6 - e_6 * fs_4_221_70 * r_8 * h5_p2 + e_6 * fs_7_143_10 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2, pc_46 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p5 = ph11_p5[k];

        pc_46[k] = e_0 * fs_14175_32_3 * h1_p1 + e_1 * fs_37485_32_2 * h3_p1 - e_1 * fs_8505_8_3 * r_2 * h1_p1 - e_2 * fs_1045_4_5 * h5_p1 - e_2 * fs_475_4_42 * h5_p5 - e_2 * fs_20825_16_2 * r_2 * h3_p1 + e_2 * fs_6075_8_3 * r_4 * h1_p1 - e_3 * fs_2205_286_21 * h7_p1 + e_3 * fs_945_572_77 * h7_p5 + e_3 * fs_2090_13_5 * r_2 * h5_p1 + e_3 * fs_950_13_42 * r_2 * h5_p5 + e_3 * fs_20825_44_2 * r_4 * h3_p1 - e_3 * fs_225_1_3 * r_6 * h1_p1 + e_4 * fs_42987_4862_15 * h9_p1 + e_4 * fs_2415_4862_858 * h9_p5 + e_4 * fs_6615_2431_21 * r_2 * h7_p1 - e_4 * fs_2835_4862_77 * r_2 * h7_p5 - e_4 * fs_418_13_5 * r_4 * h5_p1 - e_4 * fs_190_13_42 * r_4 * h5_p5 - e_4 * fs_20825_286_2 * r_6 * h3_p1 + e_4 * fs_675_22_3 * r_8 * h1_p1 + e_5 * fs_3375_8398_22 * h11_p1 - e_5 * fs_1125_4199_429 * h11_p5 - e_5 * fs_4094_2431_15 * r_2 * h9_p1 - e_5 * fs_230_2431_858 * r_2 * h9_p5 - e_5 * fs_13230_46189_21 * r_4 * h7_p1 + e_5 * fs_2835_46189_77 * r_4 * h7_p5 + e_5 * fs_1672_663_5 * r_6 * h5_p1 + e_5 * fs_760_663_42 * r_6 * h5_p5 + e_5 * fs_4165_858_2 * r_8 * h3_p1 - e_5 * fs_270_143_3 * r_10 * h1_p1 - e_6 * fs_135_4199_22 * r_2 * h11_p1 + e_6 * fs_90_4199_429 * r_2 * h11_p5 + e_6 * fs_178_2431_15 * r_4 * h9_p1 + e_6 * fs_10_2431_858 * r_4 * h9_p5 + e_6 * fs_420_46189_21 * r_6 * h7_p1 - e_6 * fs_90_46189_77 * r_6 * h7_p5 - e_6 * fs_44_663_5 * r_8 * h5_p1 - e_6 * fs_20_663_42 * r_8 * h5_p5 - e_6 * fs_49_429_2 * r_10 * h3_p1 + e_6 * fs_6_143_3 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ph11_0, ph11_p4, ab_2, pc_47 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p4 = ph11_p4[k];

        pc_47[k] = - e_0 * fs_14175_8_2 * h1_0 + e_1 * fs_5355_8_2 * h3_0 + e_1 * fs_8505_2_2 * r_2 * h1_0 + e_2 * fs_475_2_2 * h5_0 + e_2 * fs_285_4_70 * h5_p4 - e_2 * fs_2975_4_2 * r_2 * h3_0 - e_2 * fs_6075_2_2 * r_4 * h1_0 - e_3 * fs_90405_572_2 * h7_0 - e_3 * fs_8505_1144_462 * h7_p4 - e_3 * fs_1900_13_2 * r_2 * h5_0 - e_3 * fs_570_13_70 * r_2 * h5_p4 + e_3 * fs_2975_11_2 * r_4 * h3_0 + e_3 * fs_900_1_2 * r_6 * h1_0 + e_4 * fs_65205_2431_2 * h9_0 + e_4 * fs_1449_4862_10010 * h9_p4 + e_4 * fs_271215_4862_2 * r_2 * h7_0 + e_4 * fs_25515_9724_462 * r_2 * h7_p4 + e_4 * fs_380_13_2 * r_4 * h5_0 + e_4 * fs_114_13_70 * r_4 * h5_p4 - e_4 * fs_5950_143_2 * r_6 * h3_0 - e_4 * fs_1350_11_2 * r_8 * h1_0 + e_5 * fs_12375_4199_2 * h11_0 - e_5 * fs_1125_8398_2002 * h11_p4 - e_5 * fs_12420_2431_2 * r_2 * h9_0 - e_5 * fs_138_2431_10010 * r_2 * h9_p4 - e_5 * fs_271215_46189_2 * r_4 * h7_0 - e_5 * fs_25515_92378_462 * r_4 * h7_p4 - e_5 * fs_1520_663_2 * r_6 * h5_0 - e_5 * fs_152_221_70 * r_6 * h5_p4 + e_5 * fs_1190_429_2 * r_8 * h3_0 + e_5 * fs_1080_143_2 * r_10 * h1_0 - e_6 * fs_990_4199_2 * r_2 * h11_0 + e_6 * fs_45_4199_2002 * r_2 * h11_p4 + e_6 * fs_540_2431_2 * r_4 * h9_0 + e_6 * fs_6_2431_10010 * r_4 * h9_p4 + e_6 * fs_8610_46189_2 * r_6 * h7_0 + e_6 * fs_405_46189_462 * r_6 * h7_p4 + e_6 * fs_40_663_2 * r_8 * h5_0 + e_6 * fs_4_221_70 * r_8 * h5_p4 - e_6 * fs_28_429_2 * r_10 * h3_0 - e_6 * fs_24_143_2 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2, pc_48 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p3 = ph11_p3[k];

        pc_48[k] = - e_0 * fs_14175_32_14 * h1_p1 + e_1 * fs_5355_16_21 * h3_p1 - e_1 * fs_5355_16_35 * h3_p3 + e_1 * fs_8505_8_14 * r_2 * h1_p1 - e_2 * fs_95_2_210 * h5_p1 + e_2 * fs_1045_4_5 * h5_p3 - e_2 * fs_2975_8_21 * r_2 * h3_p1 + e_2 * fs_2975_8_35 * r_2 * h3_p3 - e_2 * fs_6075_8_14 * r_4 * h1_p1 + e_3 * fs_115605_1144_2 * h7_p1 - e_3 * fs_65205_1144_6 * h7_p3 + e_3 * fs_380_13_210 * r_2 * h5_p1 - e_3 * fs_2090_13_5 * r_2 * h5_p3 + e_3 * fs_2975_22_21 * r_4 * h3_p1 - e_3 * fs_2975_22_35 * r_4 * h3_p3 + e_3 * fs_225_1_14 * r_6 * h1_p1 - e_4 * fs_2898_2431_70 * h9_p1 + e_4 * fs_10143_4862_165 * h9_p3 - e_4 * fs_346815_9724_2 * r_2 * h7_p1 + e_4 * fs_195615_9724_6 * r_2 * h7_p3 - e_4 * fs_76_13_210 * r_4 * h5_p1 + e_4 * fs_418_13_5 * r_4 * h5_p3 - e_4 * fs_2975_143_21 * r_6 * h3_p1 + e_4 * fs_2975_143_35 * r_6 * h3_p3 - e_4 * fs_675_22_14 * r_8 * h1_p1 - e_5 * fs_1125_4199_231 * h11_p1 - e_5 * fs_525_4199_2145 * h11_p3 + e_5 * fs_552_2431_70 * r_2 * h9_p1 - e_5 * fs_966_2431_165 * r_2 * h9_p3 + e_5 * fs_346815_92378_2 * r_4 * h7_p1 - e_5 * fs_195615_92378_6 * r_4 * h7_p3 + e_5 * fs_304_663_210 * r_6 * h5_p1 - e_5 * fs_1672_663_5 * r_6 * h5_p3 + e_5 * fs_595_429_21 * r_8 * h3_p1 - e_5 * fs_595_429_35 * r_8 * h3_p3 + e_5 * fs_270_143_14 * r_10 * h1_p1 + e_6 * fs_90_4199_231 * r_2 * h11_p1 + e_6 * fs_42_4199_2145 * r_2 * h11_p3 - e_6 * fs_24_2431_70 * r_4 * h9_p1 + e_6 * fs_42_2431_165 * r_4 * h9_p3 - e_6 * fs_5505_46189_2 * r_6 * h7_p1 + e_6 * fs_3105_46189_6 * r_6 * h7_p3 - e_6 * fs_8_663_210 * r_8 * h5_p1 + e_6 * fs_44_663_5 * r_8 * h5_p3 - e_6 * fs_14_429_21 * r_10 * h3_p1 + e_6 * fs_14_429_35 * r_10 * h3_p3 - e_6 * fs_6_143_14 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m2, ab_2, pc_49 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m2 = ph11_m2[k];

        pc_49[k] = - e_1 * fs_5355_16_14 * h3_m2 + e_2 * fs_475_2_2 * h5_m2 + e_2 * fs_2975_8_14 * r_2 * h3_m2 - e_3 * fs_2835_572_5 * h7_m2 - e_3 * fs_1900_13_2 * r_2 * h5_m2 - e_3 * fs_2975_22_14 * r_4 * h3_m2 - e_4 * fs_2415_2431_231 * h9_m2 + e_4 * fs_8505_4862_5 * r_2 * h7_m2 + e_4 * fs_380_13_2 * r_4 * h5_m2 + e_4 * fs_2975_143_14 * r_6 * h3_m2 + e_5 * fs_675_4199_2002 * h11_m2 + e_5 * fs_460_2431_231 * r_2 * h9_m2 - e_5 * fs_8505_46189_5 * r_4 * h7_m2 - e_5 * fs_1520_663_2 * r_6 * h5_m2 - e_5 * fs_595_429_14 * r_8 * h3_m2 - e_6 * fs_54_4199_2002 * r_2 * h11_m2 - e_6 * fs_20_2431_231 * r_4 * h9_m2 + e_6 * fs_270_46189_5 * r_6 * h7_m2 + e_6 * fs_40_663_2 * r_8 * h5_m2 + e_6 * fs_14_429_14 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2, pc_50 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m1 = ph11_m1[k];

        pc_50[k] = - e_0 * fs_14175_32_14 * h1_m1 + e_1 * fs_5355_16_35 * h3_m3 + e_1 * fs_5355_16_21 * h3_m1 + e_1 * fs_8505_8_14 * r_2 * h1_m1 - e_2 * fs_1045_4_5 * h5_m3 - e_2 * fs_95_2_210 * h5_m1 - e_2 * fs_2975_8_35 * r_2 * h3_m3 - e_2 * fs_2975_8_21 * r_2 * h3_m1 - e_2 * fs_6075_8_14 * r_4 * h1_m1 + e_3 * fs_65205_1144_6 * h7_m3 + e_3 * fs_115605_1144_2 * h7_m1 + e_3 * fs_2090_13_5 * r_2 * h5_m3 + e_3 * fs_380_13_210 * r_2 * h5_m1 + e_3 * fs_2975_22_35 * r_4 * h3_m3 + e_3 * fs_2975_22_21 * r_4 * h3_m1 + e_3 * fs_225_1_14 * r_6 * h1_m1 - e_4 * fs_10143_4862_165 * h9_m3 - e_4 * fs_2898_2431_70 * h9_m1 - e_4 * fs_195615_9724_6 * r_2 * h7_m3 - e_4 * fs_346815_9724_2 * r_2 * h7_m1 - e_4 * fs_418_13_5 * r_4 * h5_m3 - e_4 * fs_76_13_210 * r_4 * h5_m1 - e_4 * fs_2975_143_35 * r_6 * h3_m3 - e_4 * fs_2975_143_21 * r_6 * h3_m1 - e_4 * fs_675_22_14 * r_8 * h1_m1 + e_5 * fs_525_4199_2145 * h11_m3 - e_5 * fs_1125_4199_231 * h11_m1 + e_5 * fs_966_2431_165 * r_2 * h9_m3 + e_5 * fs_552_2431_70 * r_2 * h9_m1 + e_5 * fs_195615_92378_6 * r_4 * h7_m3 + e_5 * fs_346815_92378_2 * r_4 * h7_m1 + e_5 * fs_1672_663_5 * r_6 * h5_m3 + e_5 * fs_304_663_210 * r_6 * h5_m1 + e_5 * fs_595_429_35 * r_8 * h3_m3 + e_5 * fs_595_429_21 * r_8 * h3_m1 + e_5 * fs_270_143_14 * r_10 * h1_m1 - e_6 * fs_42_4199_2145 * r_2 * h11_m3 + e_6 * fs_90_4199_231 * r_2 * h11_m1 - e_6 * fs_42_2431_165 * r_4 * h9_m3 - e_6 * fs_24_2431_70 * r_4 * h9_m1 - e_6 * fs_3105_46189_6 * r_6 * h7_m3 - e_6 * fs_5505_46189_2 * r_6 * h7_m1 - e_6 * fs_44_663_5 * r_8 * h5_m3 - e_6 * fs_8_663_210 * r_8 * h5_m1 - e_6 * fs_14_429_35 * r_10 * h3_m3 - e_6 * fs_14_429_21 * r_10 * h3_m1 - e_6 * fs_6_143_14 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2, pc_51 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_51[k] = - e_2 * fs_285_4_70 * h5_m4 + e_3 * fs_8505_1144_462 * h7_m4 + e_3 * fs_570_13_70 * r_2 * h5_m4 - e_4 * fs_1449_4862_10010 * h9_m4 - e_4 * fs_25515_9724_462 * r_2 * h7_m4 - e_4 * fs_114_13_70 * r_4 * h5_m4 + e_5 * fs_1125_8398_2002 * h11_m4 + e_5 * fs_138_2431_10010 * r_2 * h9_m4 + e_5 * fs_25515_92378_462 * r_4 * h7_m4 + e_5 * fs_152_221_70 * r_6 * h5_m4 - e_6 * fs_45_4199_2002 * r_2 * h11_m4 - e_6 * fs_6_2431_10010 * r_4 * h9_m4 - e_6 * fs_405_46189_462 * r_6 * h7_m4 - e_6 * fs_4_221_70 * r_8 * h5_m4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2, pc_52 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m1 = ph11_m1[k];

        pc_52[k] = - e_0 * fs_14175_32_3 * h1_m1 - e_1 * fs_37485_32_2 * h3_m1 + e_1 * fs_8505_8_3 * r_2 * h1_m1 + e_2 * fs_475_4_42 * h5_m5 + e_2 * fs_1045_4_5 * h5_m1 + e_2 * fs_20825_16_2 * r_2 * h3_m1 - e_2 * fs_6075_8_3 * r_4 * h1_m1 - e_3 * fs_945_572_77 * h7_m5 + e_3 * fs_2205_286_21 * h7_m1 - e_3 * fs_950_13_42 * r_2 * h5_m5 - e_3 * fs_2090_13_5 * r_2 * h5_m1 - e_3 * fs_20825_44_2 * r_4 * h3_m1 + e_3 * fs_225_1_3 * r_6 * h1_m1 - e_4 * fs_2415_4862_858 * h9_m5 - e_4 * fs_42987_4862_15 * h9_m1 + e_4 * fs_2835_4862_77 * r_2 * h7_m5 - e_4 * fs_6615_2431_21 * r_2 * h7_m1 + e_4 * fs_190_13_42 * r_4 * h5_m5 + e_4 * fs_418_13_5 * r_4 * h5_m1 + e_4 * fs_20825_286_2 * r_6 * h3_m1 - e_4 * fs_675_22_3 * r_8 * h1_m1 + e_5 * fs_1125_4199_429 * h11_m5 - e_5 * fs_3375_8398_22 * h11_m1 + e_5 * fs_230_2431_858 * r_2 * h9_m5 + e_5 * fs_4094_2431_15 * r_2 * h9_m1 - e_5 * fs_2835_46189_77 * r_4 * h7_m5 + e_5 * fs_13230_46189_21 * r_4 * h7_m1 - e_5 * fs_760_663_42 * r_6 * h5_m5 - e_5 * fs_1672_663_5 * r_6 * h5_m1 - e_5 * fs_4165_858_2 * r_8 * h3_m1 + e_5 * fs_270_143_3 * r_10 * h1_m1 - e_6 * fs_90_4199_429 * r_2 * h11_m5 + e_6 * fs_135_4199_22 * r_2 * h11_m1 - e_6 * fs_10_2431_858 * r_4 * h9_m5 - e_6 * fs_178_2431_15 * r_4 * h9_m1 + e_6 * fs_90_46189_77 * r_6 * h7_m5 - e_6 * fs_420_46189_21 * r_6 * h7_m1 + e_6 * fs_20_663_42 * r_8 * h5_m5 + e_6 * fs_44_663_5 * r_8 * h5_m1 + e_6 * fs_49_429_2 * r_10 * h3_m1 - e_6 * fs_6_143_3 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2, pc_53 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_53[k] = e_1 * fs_16065_32_10 * h3_m2 + e_2 * fs_285_4_70 * h5_m2 - e_2 * fs_8925_16_10 * r_2 * h3_m2 - e_3 * fs_315_44_1001 * h7_m6 - e_3 * fs_15435_286_7 * h7_m2 - e_3 * fs_570_13_70 * r_2 * h5_m2 + e_3 * fs_8925_44_10 * r_4 * h3_m2 + e_4 * fs_1449_2431_715 * h9_m6 - e_4 * fs_5796_2431_165 * h9_m2 + e_4 * fs_945_374_1001 * r_2 * h7_m6 + e_4 * fs_46305_2431_7 * r_2 * h7_m2 + e_4 * fs_114_13_70 * r_4 * h5_m2 - e_4 * fs_8925_286_10 * r_6 * h3_m2 + e_5 * fs_375_4199_2431 * h11_m6 - e_5 * fs_225_8398_1430 * h11_m2 - e_5 * fs_276_2431_715 * r_2 * h9_m6 + e_5 * fs_1104_2431_165 * r_2 * h9_m2 - e_5 * fs_945_3553_1001 * r_4 * h7_m6 - e_5 * fs_92610_46189_7 * r_4 * h7_m2 - e_5 * fs_152_221_70 * r_6 * h5_m2 + e_5 * fs_595_286_10 * r_8 * h3_m2 - e_6 * fs_30_4199_2431 * r_2 * h11_m6 + e_6 * fs_9_4199_1430 * r_2 * h11_m2 + e_6 * fs_12_2431_715 * r_4 * h9_m6 - e_6 * fs_48_2431_165 * r_4 * h9_m2 + e_6 * fs_30_3553_1001 * r_6 * h7_m6 + e_6 * fs_2940_46189_7 * r_6 * h7_m2 + e_6 * fs_4_221_70 * r_8 * h5_m2 - e_6 * fs_7_143_10 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2, pc_54 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_54[k] = - e_1 * fs_5355_32_6 * h3_m3 - e_2 * fs_475_4_42 * h5_m3 + e_2 * fs_2975_16_6 * r_2 * h3_m3 + e_3 * fs_2205_572_715 * h7_m7 - e_3 * fs_11025_286_35 * h7_m3 + e_3 * fs_950_13_42 * r_2 * h5_m3 - e_3 * fs_2975_44_6 * r_4 * h3_m3 + e_4 * fs_7245_4862_858 * h9_m7 - e_4 * fs_7245_4862_154 * h9_m3 - e_4 * fs_6615_4862_715 * r_2 * h7_m7 + e_4 * fs_33075_2431_35 * r_2 * h7_m3 - e_4 * fs_190_13_42 * r_4 * h5_m3 + e_4 * fs_2975_286_6 * r_6 * h3_m3 + e_5 * fs_225_4199_2431 * h11_m7 - e_5 * fs_75_8398_2002 * h11_m3 - e_5 * fs_690_2431_858 * r_2 * h9_m7 + e_5 * fs_690_2431_154 * r_2 * h9_m3 + e_5 * fs_6615_46189_715 * r_4 * h7_m7 - e_5 * fs_66150_46189_35 * r_4 * h7_m3 + e_5 * fs_760_663_42 * r_6 * h5_m3 - e_5 * fs_595_858_6 * r_8 * h3_m3 - e_6 * fs_18_4199_2431 * r_2 * h11_m7 + e_6 * fs_3_4199_2002 * r_2 * h11_m3 + e_6 * fs_30_2431_858 * r_4 * h9_m7 - e_6 * fs_30_2431_154 * r_4 * h9_m3 - e_6 * fs_210_46189_715 * r_6 * h7_m7 + e_6 * fs_2100_46189_35 * r_6 * h7_m3 - e_6 * fs_20_663_42 * r_8 * h5_m3 + e_6 * fs_7_429_6 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2, pc_55 : simd::cache_line_size())
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

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_55[k] = - e_2 * fs_285_8_210 * h5_p4 - e_3 * fs_11025_572_154 * h7_p4 - e_3 * fs_1575_286_1001 * h7_p6 + e_3 * fs_285_13_210 * r_2 * h5_p4 - e_4 * fs_1449_9724_30030 * h9_p4 - e_4 * fs_7245_4862_715 * h9_p6 + e_4 * fs_33075_4862_154 * r_2 * h7_p4 + e_4 * fs_4725_2431_1001 * r_2 * h7_p6 - e_4 * fs_57_13_210 * r_4 * h5_p4 - e_5 * fs_75_8398_6006 * h11_p4 - e_5 * fs_150_4199_2431 * h11_p6 + e_5 * fs_69_2431_30030 * r_2 * h9_p4 + e_5 * fs_690_2431_715 * r_2 * h9_p6 - e_5 * fs_33075_46189_154 * r_4 * h7_p4 - e_5 * fs_9450_46189_1001 * r_4 * h7_p6 + e_5 * fs_76_221_210 * r_6 * h5_p4 + e_6 * fs_3_4199_6006 * r_2 * h11_p4 + e_6 * fs_12_4199_2431 * r_2 * h11_p6 - e_6 * fs_3_2431_30030 * r_4 * h9_p4 - e_6 * fs_30_2431_715 * r_4 * h9_p6 + e_6 * fs_1050_46189_154 * r_6 * h7_p4 + e_6 * fs_300_46189_1001 * r_6 * h7_p6 - e_6 * fs_2_221_210 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2, pc_56 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_56[k] = e_1 * fs_5355_16_6 * h3_p3 + e_2 * fs_1045_8_42 * h5_p3 + e_2 * fs_285_8_210 * h5_p5 - e_2 * fs_2975_8_6 * r_2 * h3_p3 - e_3 * fs_2205_286_35 * h7_p3 + e_3 * fs_1260_143_385 * h7_p5 - e_3 * fs_1045_13_42 * r_2 * h5_p3 - e_3 * fs_285_13_210 * r_2 * h5_p5 + e_3 * fs_2975_22_6 * r_4 * h3_p3 - e_4 * fs_27531_9724_154 * h9_p3 - e_4 * fs_4347_9724_4290 * h9_p5 + e_4 * fs_6615_2431_35 * r_2 * h7_p3 - e_4 * fs_7560_2431_385 * r_2 * h7_p5 + e_4 * fs_209_13_42 * r_4 * h5_p3 + e_4 * fs_57_13_210 * r_4 * h5_p5 - e_4 * fs_2975_143_6 * r_6 * h3_p3 - e_5 * fs_150_4199_2002 * h11_p3 - e_5 * fs_300_4199_2145 * h11_p5 + e_5 * fs_1311_2431_154 * r_2 * h9_p3 + e_5 * fs_207_2431_4290 * r_2 * h9_p5 - e_5 * fs_13230_46189_35 * r_4 * h7_p3 + e_5 * fs_15120_46189_385 * r_4 * h7_p5 - e_5 * fs_836_663_42 * r_6 * h5_p3 - e_5 * fs_76_221_210 * r_6 * h5_p5 + e_5 * fs_595_429_6 * r_8 * h3_p3 + e_6 * fs_12_4199_2002 * r_2 * h11_p3 + e_6 * fs_24_4199_2145 * r_2 * h11_p5 - e_6 * fs_57_2431_154 * r_4 * h9_p3 - e_6 * fs_9_2431_4290 * r_4 * h9_p5 + e_6 * fs_420_46189_35 * r_6 * h7_p3 - e_6 * fs_480_46189_385 * r_6 * h7_p5 + e_6 * fs_22_663_42 * r_8 * h5_p3 + e_6 * fs_2_221_210 * r_8 * h5_p5 - e_6 * fs_14_429_6 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2, pc_57 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_57[k] = - e_1 * fs_5355_4_2 * h3_p2 + e_2 * fs_95_4_14 * h5_p2 - e_2 * fs_1045_8_42 * h5_p4 + e_2 * fs_2975_2_2 * r_2 * h3_p2 + e_3 * fs_6615_286_35 * h7_p2 + e_3 * fs_2835_572_770 * h7_p4 - e_3 * fs_190_13_14 * r_2 * h5_p2 + e_3 * fs_1045_13_42 * r_2 * h5_p4 - e_3 * fs_5950_11_2 * r_4 * h3_p2 - e_4 * fs_25599_4862_33 * h9_p2 - e_4 * fs_483_9724_6006 * h9_p4 - e_4 * fs_19845_2431_35 * r_2 * h7_p2 - e_4 * fs_8505_4862_770 * r_2 * h7_p4 + e_4 * fs_38_13_14 * r_4 * h5_p2 - e_4 * fs_209_13_42 * r_4 * h5_p4 + e_4 * fs_11900_143_2 * r_6 * h3_p2 - e_5 * fs_675_4199_286 * h11_p2 - e_5 * fs_225_8398_30030 * h11_p4 + e_5 * fs_2438_2431_33 * r_2 * h9_p2 + e_5 * fs_23_2431_6006 * r_2 * h9_p4 + e_5 * fs_39690_46189_35 * r_4 * h7_p2 + e_5 * fs_8505_46189_770 * r_4 * h7_p4 - e_5 * fs_152_663_14 * r_6 * h5_p2 + e_5 * fs_836_663_42 * r_6 * h5_p4 - e_5 * fs_2380_429_2 * r_8 * h3_p2 + e_6 * fs_54_4199_286 * r_2 * h11_p2 + e_6 * fs_9_4199_30030 * r_2 * h11_p4 - e_6 * fs_106_2431_33 * r_4 * h9_p2 - e_6 * fs_1_2431_6006 * r_4 * h9_p4 - e_6 * fs_1260_46189_35 * r_6 * h7_p2 - e_6 * fs_270_46189_770 * r_6 * h7_p4 + e_6 * fs_4_663_14 * r_8 * h5_p2 - e_6 * fs_22_663_42 * r_8 * h5_p4 + e_6 * fs_56_429_2 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2, pc_58 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p3 = ph11_p3[k];

        pc_58[k] = e_0 * fs_14175_32_5 * h1_p1 + e_1 * fs_5355_32_30 * h3_p1 + e_1 * fs_37485_32_2 * h3_p3 - e_1 * fs_8505_8_5 * r_2 * h1_p1 - e_2 * fs_380_1_3 * h5_p1 - e_2 * fs_95_4_14 * h5_p3 - e_2 * fs_2975_16_30 * r_2 * h3_p1 - e_2 * fs_20825_16_2 * r_2 * h3_p3 + e_2 * fs_6075_8_5 * r_4 * h1_p1 + e_3 * fs_15435_572_35 * h7_p1 - e_3 * fs_2835_572_105 * h7_p3 + e_3 * fs_3040_13_3 * r_2 * h5_p1 + e_3 * fs_190_13_14 * r_2 * h5_p3 + e_3 * fs_2975_44_30 * r_4 * h3_p1 + e_3 * fs_20825_44_2 * r_4 * h3_p3 - e_3 * fs_225_1_5 * r_6 * h1_p1 - e_4 * f_33327_2431 * h9_p1 + e_4 * fs_4347_4862_462 * h9_p3 - e_4 * fs_46305_4862_35 * r_2 * h7_p1 + e_4 * fs_8505_4862_105 * r_2 * h7_p3 - e_4 * fs_608_13_3 * r_4 * h5_p1 - e_4 * fs_38_13_14 * r_4 * h5_p3 - e_4 * fs_2975_286_30 * r_6 * h3_p1 - e_4 * fs_20825_286_2 * r_6 * h3_p3 + e_4 * fs_675_22_5 * r_8 * h1_p1 - e_5 * fs_900_4199_330 * h11_p1 - e_5 * fs_300_4199_6006 * h11_p3 + e_5 * f_6348_2431 * r_2 * h9_p1 - e_5 * fs_414_2431_462 * r_2 * h9_p3 + e_5 * fs_46305_46189_35 * r_4 * h7_p1 - e_5 * fs_8505_46189_105 * r_4 * h7_p3 + e_5 * fs_2432_663_3 * r_6 * h5_p1 + e_5 * fs_152_663_14 * r_6 * h5_p3 + e_5 * fs_595_858_30 * r_8 * h3_p1 + e_5 * fs_4165_858_2 * r_8 * h3_p3 - e_5 * fs_270_143_5 * r_10 * h1_p1 + e_6 * fs_72_4199_330 * r_2 * h11_p1 + e_6 * fs_24_4199_6006 * r_2 * h11_p3 - e_6 * f_276_2431 * r_4 * h9_p1 + e_6 * fs_18_2431_462 * r_4 * h9_p3 - e_6 * fs_1470_46189_35 * r_6 * h7_p1 + e_6 * fs_270_46189_105 * r_6 * h7_p3 - e_6 * fs_64_663_3 * r_8 * h5_p1 - e_6 * fs_4_663_14 * r_8 * h5_p3 - e_6 * fs_7_429_30 * r_10 * h3_p1 - e_6 * fs_49_429_2 * r_10 * h3_p3 + e_6 * fs_6_143_5 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ph11_0, ph11_p2, ab_2, pc_59 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p2 = ph11_p2[k];

        pc_59[k] = - e_0 * fs_14175_32_35 * h1_0 + e_1 * fs_5355_16_35 * h3_0 - e_1 * fs_5355_16_21 * h3_p2 + e_1 * fs_8505_8_35 * r_2 * h1_0 - e_2 * fs_95_1_35 * h5_0 + e_2 * fs_380_1_3 * h5_p2 - e_2 * fs_2975_8_35 * r_2 * h3_0 + e_2 * fs_2975_8_21 * r_2 * h3_p2 - e_2 * fs_6075_8_35 * r_4 * h1_0 + e_3 * fs_1260_143_35 * h7_0 - e_3 * fs_8505_286_30 * h7_p2 + e_3 * fs_760_13_35 * r_2 * h5_0 - e_3 * fs_3040_13_3 * r_2 * h5_p2 + e_3 * fs_2975_22_35 * r_4 * h3_0 - e_3 * fs_2975_22_21 * r_4 * h3_p2 + e_3 * fs_225_1_35 * r_6 * h1_0 + e_4 * fs_4347_2431_35 * h9_0 + e_4 * fs_5796_2431_154 * h9_p2 - e_4 * fs_7560_2431_35 * r_2 * h7_0 + e_4 * fs_25515_2431_30 * r_2 * h7_p2 - e_4 * fs_152_13_35 * r_4 * h5_0 + e_4 * fs_608_13_3 * r_4 * h5_p2 - e_4 * fs_2975_143_35 * r_6 * h3_0 + e_4 * fs_2975_143_21 * r_6 * h3_p2 - e_4 * fs_675_22_35 * r_8 * h1_0 - e_5 * fs_4950_4199_35 * h11_0 - e_5 * fs_450_4199_3003 * h11_p2 - e_5 * fs_828_2431_35 * r_2 * h9_0 - e_5 * fs_1104_2431_154 * r_2 * h9_p2 + e_5 * fs_15120_46189_35 * r_4 * h7_0 - e_5 * fs_51030_46189_30 * r_4 * h7_p2 + e_5 * fs_608_663_35 * r_6 * h5_0 - e_5 * fs_2432_663_3 * r_6 * h5_p2 + e_5 * fs_595_429_35 * r_8 * h3_0 - e_5 * fs_595_429_21 * r_8 * h3_p2 + e_5 * fs_270_143_35 * r_10 * h1_0 + e_6 * fs_396_4199_35 * r_2 * h11_0 + e_6 * fs_36_4199_3003 * r_2 * h11_p2 + e_6 * fs_36_2431_35 * r_4 * h9_0 + e_6 * fs_48_2431_154 * r_4 * h9_p2 - e_6 * fs_480_46189_35 * r_6 * h7_0 + e_6 * fs_1620_46189_30 * r_6 * h7_p2 - e_6 * fs_16_663_35 * r_8 * h5_0 + e_6 * fs_64_663_3 * r_8 * h5_p2 - e_6 * fs_14_429_35 * r_10 * h3_0 + e_6 * fs_14_429_21 * r_10 * h3_p2 - e_6 * fs_6_143_35 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m1, ph11_m1, ab_2, pc_60 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m1 = ph11_m1[k];

        pc_60[k] = - e_0 * fs_14175_32_21 * h1_m1 + e_1 * fs_5355_16_14 * h3_m1 + e_1 * fs_8505_8_21 * r_2 * h1_m1 - e_2 * fs_95_1_35 * h5_m1 - e_2 * fs_2975_8_14 * r_2 * h3_m1 - e_2 * fs_6075_8_21 * r_4 * h1_m1 + e_3 * fs_26775_286_3 * h7_m1 + e_3 * fs_760_13_35 * r_2 * h5_m1 + e_3 * fs_2975_22_14 * r_4 * h3_m1 + e_3 * fs_225_1_21 * r_6 * h1_m1 - e_4 * fs_483_143_105 * h9_m1 - e_4 * fs_4725_143_3 * r_2 * h7_m1 - e_4 * fs_152_13_35 * r_4 * h5_m1 - e_4 * fs_2975_143_14 * r_6 * h3_m1 - e_4 * fs_675_22_21 * r_8 * h1_m1 + e_5 * fs_2700_4199_154 * h11_m1 + e_5 * fs_92_143_105 * r_2 * h9_m1 + e_5 * fs_9450_2717_3 * r_4 * h7_m1 + e_5 * fs_608_663_35 * r_6 * h5_m1 + e_5 * fs_595_429_14 * r_8 * h3_m1 + e_5 * fs_270_143_21 * r_10 * h1_m1 - e_6 * fs_216_4199_154 * r_2 * h11_m1 - e_6 * fs_4_143_105 * r_4 * h9_m1 - e_6 * fs_300_2717_3 * r_6 * h7_m1 - e_6 * fs_16_663_35 * r_8 * h5_m1 - e_6 * fs_14_429_14 * r_10 * h3_m1 - e_6 * fs_6_143_21 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m2, ab_2, pc_61 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m2 = ph11_m2[k];

        pc_61[k] = e_1 * fs_5355_16_21 * h3_m2 - e_2 * fs_380_1_3 * h5_m2 - e_2 * fs_2975_8_21 * r_2 * h3_m2 + e_3 * fs_8505_286_30 * h7_m2 + e_3 * fs_3040_13_3 * r_2 * h5_m2 + e_3 * fs_2975_22_21 * r_4 * h3_m2 - e_4 * fs_5796_2431_154 * h9_m2 - e_4 * fs_25515_2431_30 * r_2 * h7_m2 - e_4 * fs_608_13_3 * r_4 * h5_m2 - e_4 * fs_2975_143_21 * r_6 * h3_m2 + e_5 * fs_450_4199_3003 * h11_m2 + e_5 * fs_1104_2431_154 * r_2 * h9_m2 + e_5 * fs_51030_46189_30 * r_4 * h7_m2 + e_5 * fs_2432_663_3 * r_6 * h5_m2 + e_5 * fs_595_429_21 * r_8 * h3_m2 - e_6 * fs_36_4199_3003 * r_2 * h11_m2 - e_6 * fs_48_2431_154 * r_4 * h9_m2 - e_6 * fs_1620_46189_30 * r_6 * h7_m2 - e_6 * fs_64_663_3 * r_8 * h5_m2 - e_6 * fs_14_429_21 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2, pc_62 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m1 = ph11_m1[k];

        pc_62[k] = - e_0 * fs_14175_32_5 * h1_m1 - e_1 * fs_37485_32_2 * h3_m3 - e_1 * fs_5355_32_30 * h3_m1 + e_1 * fs_8505_8_5 * r_2 * h1_m1 + e_2 * fs_95_4_14 * h5_m3 + e_2 * fs_380_1_3 * h5_m1 + e_2 * fs_20825_16_2 * r_2 * h3_m3 + e_2 * fs_2975_16_30 * r_2 * h3_m1 - e_2 * fs_6075_8_5 * r_4 * h1_m1 + e_3 * fs_2835_572_105 * h7_m3 - e_3 * fs_15435_572_35 * h7_m1 - e_3 * fs_190_13_14 * r_2 * h5_m3 - e_3 * fs_3040_13_3 * r_2 * h5_m1 - e_3 * fs_20825_44_2 * r_4 * h3_m3 - e_3 * fs_2975_44_30 * r_4 * h3_m1 + e_3 * fs_225_1_5 * r_6 * h1_m1 - e_4 * fs_4347_4862_462 * h9_m3 + e_4 * f_33327_2431 * h9_m1 - e_4 * fs_8505_4862_105 * r_2 * h7_m3 + e_4 * fs_46305_4862_35 * r_2 * h7_m1 + e_4 * fs_38_13_14 * r_4 * h5_m3 + e_4 * fs_608_13_3 * r_4 * h5_m1 + e_4 * fs_20825_286_2 * r_6 * h3_m3 + e_4 * fs_2975_286_30 * r_6 * h3_m1 - e_4 * fs_675_22_5 * r_8 * h1_m1 + e_5 * fs_300_4199_6006 * h11_m3 + e_5 * fs_900_4199_330 * h11_m1 + e_5 * fs_414_2431_462 * r_2 * h9_m3 - e_5 * f_6348_2431 * r_2 * h9_m1 + e_5 * fs_8505_46189_105 * r_4 * h7_m3 - e_5 * fs_46305_46189_35 * r_4 * h7_m1 - e_5 * fs_152_663_14 * r_6 * h5_m3 - e_5 * fs_2432_663_3 * r_6 * h5_m1 - e_5 * fs_4165_858_2 * r_8 * h3_m3 - e_5 * fs_595_858_30 * r_8 * h3_m1 + e_5 * fs_270_143_5 * r_10 * h1_m1 - e_6 * fs_24_4199_6006 * r_2 * h11_m3 - e_6 * fs_72_4199_330 * r_2 * h11_m1 - e_6 * fs_18_2431_462 * r_4 * h9_m3 + e_6 * f_276_2431 * r_4 * h9_m1 - e_6 * fs_270_46189_105 * r_6 * h7_m3 + e_6 * fs_1470_46189_35 * r_6 * h7_m1 + e_6 * fs_4_663_14 * r_8 * h5_m3 + e_6 * fs_64_663_3 * r_8 * h5_m1 + e_6 * fs_49_429_2 * r_10 * h3_m3 + e_6 * fs_7_429_30 * r_10 * h3_m1 - e_6 * fs_6_143_5 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2, pc_63 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_63[k] = e_1 * fs_5355_4_2 * h3_m2 + e_2 * fs_1045_8_42 * h5_m4 - e_2 * fs_95_4_14 * h5_m2 - e_2 * fs_2975_2_2 * r_2 * h3_m2 - e_3 * fs_2835_572_770 * h7_m4 - e_3 * fs_6615_286_35 * h7_m2 - e_3 * fs_1045_13_42 * r_2 * h5_m4 + e_3 * fs_190_13_14 * r_2 * h5_m2 + e_3 * fs_5950_11_2 * r_4 * h3_m2 + e_4 * fs_483_9724_6006 * h9_m4 + e_4 * fs_25599_4862_33 * h9_m2 + e_4 * fs_8505_4862_770 * r_2 * h7_m4 + e_4 * fs_19845_2431_35 * r_2 * h7_m2 + e_4 * fs_209_13_42 * r_4 * h5_m4 - e_4 * fs_38_13_14 * r_4 * h5_m2 - e_4 * fs_11900_143_2 * r_6 * h3_m2 + e_5 * fs_225_8398_30030 * h11_m4 + e_5 * fs_675_4199_286 * h11_m2 - e_5 * fs_23_2431_6006 * r_2 * h9_m4 - e_5 * fs_2438_2431_33 * r_2 * h9_m2 - e_5 * fs_8505_46189_770 * r_4 * h7_m4 - e_5 * fs_39690_46189_35 * r_4 * h7_m2 - e_5 * fs_836_663_42 * r_6 * h5_m4 + e_5 * fs_152_663_14 * r_6 * h5_m2 + e_5 * fs_2380_429_2 * r_8 * h3_m2 - e_6 * fs_9_4199_30030 * r_2 * h11_m4 - e_6 * fs_54_4199_286 * r_2 * h11_m2 + e_6 * fs_1_2431_6006 * r_4 * h9_m4 + e_6 * fs_106_2431_33 * r_4 * h9_m2 + e_6 * fs_270_46189_770 * r_6 * h7_m4 + e_6 * fs_1260_46189_35 * r_6 * h7_m2 + e_6 * fs_22_663_42 * r_8 * h5_m4 - e_6 * fs_4_663_14 * r_8 * h5_m2 - e_6 * fs_56_429_2 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2, pc_64 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_64[k] = - e_1 * fs_5355_16_6 * h3_m3 - e_2 * fs_285_8_210 * h5_m5 - e_2 * fs_1045_8_42 * h5_m3 + e_2 * fs_2975_8_6 * r_2 * h3_m3 - e_3 * fs_1260_143_385 * h7_m5 + e_3 * fs_2205_286_35 * h7_m3 + e_3 * fs_285_13_210 * r_2 * h5_m5 + e_3 * fs_1045_13_42 * r_2 * h5_m3 - e_3 * fs_2975_22_6 * r_4 * h3_m3 + e_4 * fs_4347_9724_4290 * h9_m5 + e_4 * fs_27531_9724_154 * h9_m3 + e_4 * fs_7560_2431_385 * r_2 * h7_m5 - e_4 * fs_6615_2431_35 * r_2 * h7_m3 - e_4 * fs_57_13_210 * r_4 * h5_m5 - e_4 * fs_209_13_42 * r_4 * h5_m3 + e_4 * fs_2975_143_6 * r_6 * h3_m3 + e_5 * fs_300_4199_2145 * h11_m5 + e_5 * fs_150_4199_2002 * h11_m3 - e_5 * fs_207_2431_4290 * r_2 * h9_m5 - e_5 * fs_1311_2431_154 * r_2 * h9_m3 - e_5 * fs_15120_46189_385 * r_4 * h7_m5 + e_5 * fs_13230_46189_35 * r_4 * h7_m3 + e_5 * fs_76_221_210 * r_6 * h5_m5 + e_5 * fs_836_663_42 * r_6 * h5_m3 - e_5 * fs_595_429_6 * r_8 * h3_m3 - e_6 * fs_24_4199_2145 * r_2 * h11_m5 - e_6 * fs_12_4199_2002 * r_2 * h11_m3 + e_6 * fs_9_2431_4290 * r_4 * h9_m5 + e_6 * fs_57_2431_154 * r_4 * h9_m3 + e_6 * fs_480_46189_385 * r_6 * h7_m5 - e_6 * fs_420_46189_35 * r_6 * h7_m3 - e_6 * fs_2_221_210 * r_8 * h5_m5 - e_6 * fs_22_663_42 * r_8 * h5_m3 + e_6 * fs_14_429_6 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m5, ph5_m4, ph7_m6, ph7_m5, ph7_m4, ph9_m6, ph9_m5, ph9_m4, ph11_m6, ph11_m5, ph11_m4, ab_2, pc_65, pc_66, pc_67 : simd::cache_line_size())
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

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m4 = ph11_m4[k];

        pc_65[k] = e_2 * fs_285_8_210 * h5_m4 + e_3 * fs_1575_286_1001 * h7_m6 + e_3 * fs_11025_572_154 * h7_m4 - e_3 * fs_285_13_210 * r_2 * h5_m4 + e_4 * fs_7245_4862_715 * h9_m6 + e_4 * fs_1449_9724_30030 * h9_m4 - e_4 * fs_4725_2431_1001 * r_2 * h7_m6 - e_4 * fs_33075_4862_154 * r_2 * h7_m4 + e_4 * fs_57_13_210 * r_4 * h5_m4 + e_5 * fs_150_4199_2431 * h11_m6 + e_5 * fs_75_8398_6006 * h11_m4 - e_5 * fs_690_2431_715 * r_2 * h9_m6 - e_5 * fs_69_2431_30030 * r_2 * h9_m4 + e_5 * fs_9450_46189_1001 * r_4 * h7_m6 + e_5 * fs_33075_46189_154 * r_4 * h7_m4 - e_5 * fs_76_221_210 * r_6 * h5_m4 - e_6 * fs_12_4199_2431 * r_2 * h11_m6 - e_6 * fs_3_4199_6006 * r_2 * h11_m4 + e_6 * fs_30_2431_715 * r_4 * h9_m6 + e_6 * fs_3_2431_30030 * r_4 * h9_m4 - e_6 * fs_300_46189_1001 * r_6 * h7_m6 - e_6 * fs_1050_46189_154 * r_6 * h7_m4 + e_6 * fs_2_221_210 * r_8 * h5_m4;

        pc_66[k] = e_2 * f_1425_4 * h5_m5 + e_3 * fs_11025_286_66 * h7_m5 - e_3 * f_2850_13 * r_2 * h5_m5 + e_4 * fs_7245_4862_1001 * h9_m5 - e_4 * fs_33075_2431_66 * r_2 * h7_m5 + e_4 * f_570_13 * r_4 * h5_m5 + e_5 * fs_150_4199_2002 * h11_m5 - e_5 * fs_690_2431_1001 * r_2 * h9_m5 + e_5 * fs_66150_46189_66 * r_4 * h7_m5 - e_5 * f_760_221 * r_6 * h5_m5 - e_6 * fs_12_4199_2002 * r_2 * h11_m5 + e_6 * fs_30_2431_1001 * r_4 * h9_m5 - e_6 * fs_2100_46189_66 * r_6 * h7_m5 + e_6 * f_20_221 * r_8 * h5_m5;

        pc_67[k] = - e_2 * f_1140_1 * h5_m4 - e_3 * fs_2205_286_165 * h7_m4 + e_3 * f_9120_13 * r_2 * h5_m4 + e_4 * fs_10143_2431_143 * h9_m4 + e_4 * fs_6615_2431_165 * r_2 * h7_m4 - e_4 * f_1824_13 * r_4 * h5_m4 + e_5 * fs_525_4199_715 * h11_m4 - e_5 * fs_1932_2431_143 * r_2 * h9_m4 - e_5 * fs_13230_46189_165 * r_4 * h7_m4 + e_5 * f_2432_221 * r_6 * h5_m4 - e_6 * fs_42_4199_715 * r_2 * h11_m4 + e_6 * fs_84_2431_143 * r_4 * h9_m4 + e_6 * fs_420_46189_165 * r_6 * h7_m4 - e_6 * f_64_221 * r_8 * h5_m4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m3, ph7_m2, ph9_m3, ph9_m2, ph11_m3, ph11_m2, ab_2, pc_68, pc_69 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m2 = ph11_m2[k];

        pc_68[k] = e_1 * fs_5355_8_7 * h3_m3 + e_2 * f_2755_4 * h5_m3 - e_2 * fs_2975_4_7 * r_2 * h3_m3 - e_3 * fs_6615_143_30 * h7_m3 - e_3 * f_5510_13 * r_2 * h5_m3 + e_3 * fs_2975_11_7 * r_4 * h3_m3 + e_4 * fs_23667_4862_33 * h9_m3 + e_4 * fs_39690_2431_30 * r_2 * h7_m3 + e_4 * f_1102_13 * r_4 * h5_m3 - e_4 * fs_5950_143_7 * r_6 * h3_m3 + e_5 * fs_1050_4199_429 * h11_m3 - e_5 * fs_2254_2431_33 * r_2 * h9_m3 - e_5 * fs_79380_46189_30 * r_4 * h7_m3 - e_5 * f_4408_663 * r_6 * h5_m3 + e_5 * fs_1190_429_7 * r_8 * h3_m3 - e_6 * fs_84_4199_429 * r_2 * h11_m3 + e_6 * fs_98_2431_33 * r_4 * h9_m3 + e_6 * fs_2520_46189_30 * r_6 * h7_m3 + e_6 * f_116_663 * r_8 * h5_m3 - e_6 * fs_28_429_7 * r_10 * h3_m3;

        pc_69[k] = - e_1 * fs_16065_16_7 * h3_m2 + e_2 * f_855_1 * h5_m2 + e_2 * fs_8925_8_7 * r_2 * h3_m2 - e_3 * fs_19845_572_10 * h7_m2 - e_3 * f_6840_13 * r_2 * h5_m2 - e_3 * fs_8925_22_7 * r_4 * h3_m2 - e_4 * fs_483_2431_462 * h9_m2 + e_4 * fs_59535_4862_10 * r_2 * h7_m2 + e_4 * f_1368_13 * r_4 * h5_m2 + e_4 * fs_8925_143_7 * r_6 * h3_m2 + e_5 * fs_900_4199_1001 * h11_m2 + e_5 * fs_92_2431_462 * r_2 * h9_m2 - e_5 * fs_59535_46189_10 * r_4 * h7_m2 - e_5 * f_1824_221 * r_6 * h5_m2 - e_5 * fs_595_143_7 * r_8 * h3_m2 - e_6 * fs_72_4199_1001 * r_2 * h11_m2 - e_6 * fs_4_2431_462 * r_4 * h9_m2 + e_6 * fs_1890_46189_10 * r_6 * h7_m2 + e_6 * f_48_221 * r_8 * h5_m2 + e_6 * fs_14_143_7 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph5_m1, ph7_m1, ph9_m1, ph11_m1, ab_2, pc_70 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m1 = ph11_m1[k];

        pc_70[k] = e_0 * fs_14175_32_15 * h1_m1 - e_1 * fs_8505_8_15 * r_2 * h1_m1 - e_2 * f_285_1 * h5_m1 + e_2 * fs_6075_8_15 * r_4 * h1_m1 + e_3 * fs_315_26_105 * h7_m1 + e_3 * f_2280_13 * r_2 * h5_m1 - e_3 * fs_225_1_15 * r_6 * h1_m1 - e_4 * fs_3381_187_3 * h9_m1 - e_4 * fs_945_221_105 * r_2 * h7_m1 - e_4 * f_456_13 * r_4 * h5_m1 + e_4 * fs_675_22_15 * r_8 * h1_m1 + e_5 * fs_3150_4199_110 * h11_m1 + e_5 * fs_644_187_3 * r_2 * h9_m1 + e_5 * fs_1890_4199_105 * r_4 * h7_m1 + e_5 * f_608_221 * r_6 * h5_m1 - e_5 * fs_270_143_15 * r_10 * h1_m1 - e_6 * fs_252_4199_110 * r_2 * h11_m1 - e_6 * fs_28_187_3 * r_4 * h9_m1 - e_6 * fs_60_4199_105 * r_6 * h7_m1 - e_6 * f_16_221 * r_8 * h5_m1 + e_6 * fs_6_143_15 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph11_0, ab_2, pc_71 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_0 = ph11_0[k];

        pc_71[k] = - e_0 * f_42525_16 * h1_0 + e_1 * f_37485_16 * h3_0 + e_1 * f_25515_4 * r_2 * h1_0 - e_2 * f_950_1 * h5_0 - e_2 * f_20825_8 * r_2 * h3_0 - e_2 * f_18225_4 * r_4 * h1_0 + e_3 * f_33075_143 * h7_0 + e_3 * f_7600_13 * r_2 * h5_0 + e_3 * f_20825_22 * r_4 * h3_0 + e_3 * f_1350_1 * r_6 * h1_0 - e_4 * f_101430_2431 * h9_0 - e_4 * f_198450_2431 * r_2 * h7_0 - e_4 * f_1520_13 * r_4 * h5_0 - e_4 * f_20825_143 * r_6 * h3_0 - e_4 * f_2025_11 * r_8 * h1_0 + e_5 * f_34650_4199 * h11_0 + e_5 * f_19320_2431 * r_2 * h9_0 + e_5 * f_396900_46189 * r_4 * h7_0 + e_5 * f_6080_663 * r_6 * h5_0 + e_5 * f_4165_429 * r_8 * h3_0 + e_5 * f_1620_143 * r_10 * h1_0 - e_6 * f_2772_4199 * r_2 * h11_0 - e_6 * f_840_2431 * r_4 * h9_0 - e_6 * f_12600_46189 * r_6 * h7_0 - e_6 * f_160_663 * r_8 * h5_0 - e_6 * f_98_429 * r_10 * h3_0 - e_6 * f_36_143 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph9_p1, ph9_p2, ph11_p1, ph11_p2, ab_2, pc_72, pc_73 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p2 = ph11_p2[k];

        pc_72[k] = e_0 * fs_14175_32_15 * h1_p1 - e_1 * fs_8505_8_15 * r_2 * h1_p1 - e_2 * f_285_1 * h5_p1 + e_2 * fs_6075_8_15 * r_4 * h1_p1 + e_3 * fs_315_26_105 * h7_p1 + e_3 * f_2280_13 * r_2 * h5_p1 - e_3 * fs_225_1_15 * r_6 * h1_p1 - e_4 * fs_3381_187_3 * h9_p1 - e_4 * fs_945_221_105 * r_2 * h7_p1 - e_4 * f_456_13 * r_4 * h5_p1 + e_4 * fs_675_22_15 * r_8 * h1_p1 + e_5 * fs_3150_4199_110 * h11_p1 + e_5 * fs_644_187_3 * r_2 * h9_p1 + e_5 * fs_1890_4199_105 * r_4 * h7_p1 + e_5 * f_608_221 * r_6 * h5_p1 - e_5 * fs_270_143_15 * r_10 * h1_p1 - e_6 * fs_252_4199_110 * r_2 * h11_p1 - e_6 * fs_28_187_3 * r_4 * h9_p1 - e_6 * fs_60_4199_105 * r_6 * h7_p1 - e_6 * f_16_221 * r_8 * h5_p1 + e_6 * fs_6_143_15 * r_12 * h1_p1;

        pc_73[k] = - e_1 * fs_16065_16_7 * h3_p2 + e_2 * f_855_1 * h5_p2 + e_2 * fs_8925_8_7 * r_2 * h3_p2 - e_3 * fs_19845_572_10 * h7_p2 - e_3 * f_6840_13 * r_2 * h5_p2 - e_3 * fs_8925_22_7 * r_4 * h3_p2 - e_4 * fs_483_2431_462 * h9_p2 + e_4 * fs_59535_4862_10 * r_2 * h7_p2 + e_4 * f_1368_13 * r_4 * h5_p2 + e_4 * fs_8925_143_7 * r_6 * h3_p2 + e_5 * fs_900_4199_1001 * h11_p2 + e_5 * fs_92_2431_462 * r_2 * h9_p2 - e_5 * fs_59535_46189_10 * r_4 * h7_p2 - e_5 * f_1824_221 * r_6 * h5_p2 - e_5 * fs_595_143_7 * r_8 * h3_p2 - e_6 * fs_72_4199_1001 * r_2 * h11_p2 - e_6 * fs_4_2431_462 * r_4 * h9_p2 + e_6 * fs_1890_46189_10 * r_6 * h7_p2 + e_6 * f_48_221 * r_8 * h5_p2 + e_6 * fs_14_143_7 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph9_p3, ph9_p4, ph11_p3, ph11_p4, ab_2, pc_74, pc_75 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p4 = ph11_p4[k];

        pc_74[k] = e_1 * fs_5355_8_7 * h3_p3 + e_2 * f_2755_4 * h5_p3 - e_2 * fs_2975_4_7 * r_2 * h3_p3 - e_3 * fs_6615_143_30 * h7_p3 - e_3 * f_5510_13 * r_2 * h5_p3 + e_3 * fs_2975_11_7 * r_4 * h3_p3 + e_4 * fs_23667_4862_33 * h9_p3 + e_4 * fs_39690_2431_30 * r_2 * h7_p3 + e_4 * f_1102_13 * r_4 * h5_p3 - e_4 * fs_5950_143_7 * r_6 * h3_p3 + e_5 * fs_1050_4199_429 * h11_p3 - e_5 * fs_2254_2431_33 * r_2 * h9_p3 - e_5 * fs_79380_46189_30 * r_4 * h7_p3 - e_5 * f_4408_663 * r_6 * h5_p3 + e_5 * fs_1190_429_7 * r_8 * h3_p3 - e_6 * fs_84_4199_429 * r_2 * h11_p3 + e_6 * fs_98_2431_33 * r_4 * h9_p3 + e_6 * fs_2520_46189_30 * r_6 * h7_p3 + e_6 * f_116_663 * r_8 * h5_p3 - e_6 * fs_28_429_7 * r_10 * h3_p3;

        pc_75[k] = - e_2 * f_1140_1 * h5_p4 - e_3 * fs_2205_286_165 * h7_p4 + e_3 * f_9120_13 * r_2 * h5_p4 + e_4 * fs_10143_2431_143 * h9_p4 + e_4 * fs_6615_2431_165 * r_2 * h7_p4 - e_4 * f_1824_13 * r_4 * h5_p4 + e_5 * fs_525_4199_715 * h11_p4 - e_5 * fs_1932_2431_143 * r_2 * h9_p4 - e_5 * fs_13230_46189_165 * r_4 * h7_p4 + e_5 * f_2432_221 * r_6 * h5_p4 - e_6 * fs_42_4199_715 * r_2 * h11_p4 + e_6 * fs_84_2431_143 * r_4 * h9_p4 + e_6 * fs_420_46189_165 * r_6 * h7_p4 - e_6 * f_64_221 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph5_p5, ph7_m6, ph7_m4, ph7_p5, ph9_m6, ph9_m4, ph9_p5, ph11_m6, ph11_m4, ph11_p5, ab_2, pc_76, pc_77 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_p5 = ph11_p5[k];

        pc_76[k] = e_2 * f_1425_4 * h5_p5 + e_3 * fs_11025_286_66 * h7_p5 - e_3 * f_2850_13 * r_2 * h5_p5 + e_4 * fs_7245_4862_1001 * h9_p5 - e_4 * fs_33075_2431_66 * r_2 * h7_p5 + e_4 * f_570_13 * r_4 * h5_p5 + e_5 * fs_150_4199_2002 * h11_p5 - e_5 * fs_690_2431_1001 * r_2 * h9_p5 + e_5 * fs_66150_46189_66 * r_4 * h7_p5 - e_5 * f_760_221 * r_6 * h5_p5 - e_6 * fs_12_4199_2002 * r_2 * h11_p5 + e_6 * fs_30_2431_1001 * r_4 * h9_p5 - e_6 * fs_2100_46189_66 * r_6 * h7_p5 + e_6 * f_20_221 * r_8 * h5_p5;

        pc_77[k] = - e_2 * fs_285_8_210 * h5_m4 + e_3 * fs_1575_286_1001 * h7_m6 - e_3 * fs_11025_572_154 * h7_m4 + e_3 * fs_285_13_210 * r_2 * h5_m4 + e_4 * fs_7245_4862_715 * h9_m6 - e_4 * fs_1449_9724_30030 * h9_m4 - e_4 * fs_4725_2431_1001 * r_2 * h7_m6 + e_4 * fs_33075_4862_154 * r_2 * h7_m4 - e_4 * fs_57_13_210 * r_4 * h5_m4 + e_5 * fs_150_4199_2431 * h11_m6 - e_5 * fs_75_8398_6006 * h11_m4 - e_5 * fs_690_2431_715 * r_2 * h9_m6 + e_5 * fs_69_2431_30030 * r_2 * h9_m4 + e_5 * fs_9450_46189_1001 * r_4 * h7_m6 - e_5 * fs_33075_46189_154 * r_4 * h7_m4 + e_5 * fs_76_221_210 * r_6 * h5_m4 - e_6 * fs_12_4199_2431 * r_2 * h11_m6 + e_6 * fs_3_4199_6006 * r_2 * h11_m4 + e_6 * fs_30_2431_715 * r_4 * h9_m6 - e_6 * fs_3_2431_30030 * r_4 * h9_m4 - e_6 * fs_300_46189_1001 * r_6 * h7_m6 + e_6 * fs_1050_46189_154 * r_6 * h7_m4 - e_6 * fs_2_221_210 * r_8 * h5_m4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2, pc_78 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_78[k] = e_1 * fs_5355_16_6 * h3_m3 - e_2 * fs_285_8_210 * h5_m5 + e_2 * fs_1045_8_42 * h5_m3 - e_2 * fs_2975_8_6 * r_2 * h3_m3 - e_3 * fs_1260_143_385 * h7_m5 - e_3 * fs_2205_286_35 * h7_m3 + e_3 * fs_285_13_210 * r_2 * h5_m5 - e_3 * fs_1045_13_42 * r_2 * h5_m3 + e_3 * fs_2975_22_6 * r_4 * h3_m3 + e_4 * fs_4347_9724_4290 * h9_m5 - e_4 * fs_27531_9724_154 * h9_m3 + e_4 * fs_7560_2431_385 * r_2 * h7_m5 + e_4 * fs_6615_2431_35 * r_2 * h7_m3 - e_4 * fs_57_13_210 * r_4 * h5_m5 + e_4 * fs_209_13_42 * r_4 * h5_m3 - e_4 * fs_2975_143_6 * r_6 * h3_m3 + e_5 * fs_300_4199_2145 * h11_m5 - e_5 * fs_150_4199_2002 * h11_m3 - e_5 * fs_207_2431_4290 * r_2 * h9_m5 + e_5 * fs_1311_2431_154 * r_2 * h9_m3 - e_5 * fs_15120_46189_385 * r_4 * h7_m5 - e_5 * fs_13230_46189_35 * r_4 * h7_m3 + e_5 * fs_76_221_210 * r_6 * h5_m5 - e_5 * fs_836_663_42 * r_6 * h5_m3 + e_5 * fs_595_429_6 * r_8 * h3_m3 - e_6 * fs_24_4199_2145 * r_2 * h11_m5 + e_6 * fs_12_4199_2002 * r_2 * h11_m3 + e_6 * fs_9_2431_4290 * r_4 * h9_m5 - e_6 * fs_57_2431_154 * r_4 * h9_m3 + e_6 * fs_480_46189_385 * r_6 * h7_m5 + e_6 * fs_420_46189_35 * r_6 * h7_m3 - e_6 * fs_2_221_210 * r_8 * h5_m5 + e_6 * fs_22_663_42 * r_8 * h5_m3 - e_6 * fs_14_429_6 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2, pc_79 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_79[k] = - e_1 * fs_5355_4_2 * h3_m2 + e_2 * fs_1045_8_42 * h5_m4 + e_2 * fs_95_4_14 * h5_m2 + e_2 * fs_2975_2_2 * r_2 * h3_m2 - e_3 * fs_2835_572_770 * h7_m4 + e_3 * fs_6615_286_35 * h7_m2 - e_3 * fs_1045_13_42 * r_2 * h5_m4 - e_3 * fs_190_13_14 * r_2 * h5_m2 - e_3 * fs_5950_11_2 * r_4 * h3_m2 + e_4 * fs_483_9724_6006 * h9_m4 - e_4 * fs_25599_4862_33 * h9_m2 + e_4 * fs_8505_4862_770 * r_2 * h7_m4 - e_4 * fs_19845_2431_35 * r_2 * h7_m2 + e_4 * fs_209_13_42 * r_4 * h5_m4 + e_4 * fs_38_13_14 * r_4 * h5_m2 + e_4 * fs_11900_143_2 * r_6 * h3_m2 + e_5 * fs_225_8398_30030 * h11_m4 - e_5 * fs_675_4199_286 * h11_m2 - e_5 * fs_23_2431_6006 * r_2 * h9_m4 + e_5 * fs_2438_2431_33 * r_2 * h9_m2 - e_5 * fs_8505_46189_770 * r_4 * h7_m4 + e_5 * fs_39690_46189_35 * r_4 * h7_m2 - e_5 * fs_836_663_42 * r_6 * h5_m4 - e_5 * fs_152_663_14 * r_6 * h5_m2 - e_5 * fs_2380_429_2 * r_8 * h3_m2 - e_6 * fs_9_4199_30030 * r_2 * h11_m4 + e_6 * fs_54_4199_286 * r_2 * h11_m2 + e_6 * fs_1_2431_6006 * r_4 * h9_m4 - e_6 * fs_106_2431_33 * r_4 * h9_m2 + e_6 * fs_270_46189_770 * r_6 * h7_m4 - e_6 * fs_1260_46189_35 * r_6 * h7_m2 + e_6 * fs_22_663_42 * r_8 * h5_m4 + e_6 * fs_4_663_14 * r_8 * h5_m2 + e_6 * fs_56_429_2 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2, pc_80 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m1 = ph11_m1[k];

        pc_80[k] = e_0 * fs_14175_32_5 * h1_m1 - e_1 * fs_37485_32_2 * h3_m3 + e_1 * fs_5355_32_30 * h3_m1 - e_1 * fs_8505_8_5 * r_2 * h1_m1 + e_2 * fs_95_4_14 * h5_m3 - e_2 * fs_380_1_3 * h5_m1 + e_2 * fs_20825_16_2 * r_2 * h3_m3 - e_2 * fs_2975_16_30 * r_2 * h3_m1 + e_2 * fs_6075_8_5 * r_4 * h1_m1 + e_3 * fs_2835_572_105 * h7_m3 + e_3 * fs_15435_572_35 * h7_m1 - e_3 * fs_190_13_14 * r_2 * h5_m3 + e_3 * fs_3040_13_3 * r_2 * h5_m1 - e_3 * fs_20825_44_2 * r_4 * h3_m3 + e_3 * fs_2975_44_30 * r_4 * h3_m1 - e_3 * fs_225_1_5 * r_6 * h1_m1 - e_4 * fs_4347_4862_462 * h9_m3 - e_4 * f_33327_2431 * h9_m1 - e_4 * fs_8505_4862_105 * r_2 * h7_m3 - e_4 * fs_46305_4862_35 * r_2 * h7_m1 + e_4 * fs_38_13_14 * r_4 * h5_m3 - e_4 * fs_608_13_3 * r_4 * h5_m1 + e_4 * fs_20825_286_2 * r_6 * h3_m3 - e_4 * fs_2975_286_30 * r_6 * h3_m1 + e_4 * fs_675_22_5 * r_8 * h1_m1 + e_5 * fs_300_4199_6006 * h11_m3 - e_5 * fs_900_4199_330 * h11_m1 + e_5 * fs_414_2431_462 * r_2 * h9_m3 + e_5 * f_6348_2431 * r_2 * h9_m1 + e_5 * fs_8505_46189_105 * r_4 * h7_m3 + e_5 * fs_46305_46189_35 * r_4 * h7_m1 - e_5 * fs_152_663_14 * r_6 * h5_m3 + e_5 * fs_2432_663_3 * r_6 * h5_m1 - e_5 * fs_4165_858_2 * r_8 * h3_m3 + e_5 * fs_595_858_30 * r_8 * h3_m1 - e_5 * fs_270_143_5 * r_10 * h1_m1 - e_6 * fs_24_4199_6006 * r_2 * h11_m3 + e_6 * fs_72_4199_330 * r_2 * h11_m1 - e_6 * fs_18_2431_462 * r_4 * h9_m3 - e_6 * f_276_2431 * r_4 * h9_m1 - e_6 * fs_270_46189_105 * r_6 * h7_m3 - e_6 * fs_1470_46189_35 * r_6 * h7_m1 + e_6 * fs_4_663_14 * r_8 * h5_m3 - e_6 * fs_64_663_3 * r_8 * h5_m1 + e_6 * fs_49_429_2 * r_10 * h3_m3 - e_6 * fs_7_429_30 * r_10 * h3_m1 + e_6 * fs_6_143_5 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m2, ab_2, pc_81 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m2 = ph11_m2[k];

        pc_81[k] = e_1 * fs_5355_16_21 * h3_m2 - e_2 * fs_380_1_3 * h5_m2 - e_2 * fs_2975_8_21 * r_2 * h3_m2 + e_3 * fs_8505_286_30 * h7_m2 + e_3 * fs_3040_13_3 * r_2 * h5_m2 + e_3 * fs_2975_22_21 * r_4 * h3_m2 - e_4 * fs_5796_2431_154 * h9_m2 - e_4 * fs_25515_2431_30 * r_2 * h7_m2 - e_4 * fs_608_13_3 * r_4 * h5_m2 - e_4 * fs_2975_143_21 * r_6 * h3_m2 + e_5 * fs_450_4199_3003 * h11_m2 + e_5 * fs_1104_2431_154 * r_2 * h9_m2 + e_5 * fs_51030_46189_30 * r_4 * h7_m2 + e_5 * fs_2432_663_3 * r_6 * h5_m2 + e_5 * fs_595_429_21 * r_8 * h3_m2 - e_6 * fs_36_4199_3003 * r_2 * h11_m2 - e_6 * fs_48_2431_154 * r_4 * h9_m2 - e_6 * fs_1620_46189_30 * r_6 * h7_m2 - e_6 * fs_64_663_3 * r_8 * h5_m2 - e_6 * fs_14_429_21 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph11_p1, ab_2, pc_82 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_p1 = ph11_p1[k];

        pc_82[k] = - e_0 * fs_14175_32_21 * h1_p1 + e_1 * fs_5355_16_14 * h3_p1 + e_1 * fs_8505_8_21 * r_2 * h1_p1 - e_2 * fs_95_1_35 * h5_p1 - e_2 * fs_2975_8_14 * r_2 * h3_p1 - e_2 * fs_6075_8_21 * r_4 * h1_p1 + e_3 * fs_26775_286_3 * h7_p1 + e_3 * fs_760_13_35 * r_2 * h5_p1 + e_3 * fs_2975_22_14 * r_4 * h3_p1 + e_3 * fs_225_1_21 * r_6 * h1_p1 - e_4 * fs_483_143_105 * h9_p1 - e_4 * fs_4725_143_3 * r_2 * h7_p1 - e_4 * fs_152_13_35 * r_4 * h5_p1 - e_4 * fs_2975_143_14 * r_6 * h3_p1 - e_4 * fs_675_22_21 * r_8 * h1_p1 + e_5 * fs_2700_4199_154 * h11_p1 + e_5 * fs_92_143_105 * r_2 * h9_p1 + e_5 * fs_9450_2717_3 * r_4 * h7_p1 + e_5 * fs_608_663_35 * r_6 * h5_p1 + e_5 * fs_595_429_14 * r_8 * h3_p1 + e_5 * fs_270_143_21 * r_10 * h1_p1 - e_6 * fs_216_4199_154 * r_2 * h11_p1 - e_6 * fs_4_143_105 * r_4 * h9_p1 - e_6 * fs_300_2717_3 * r_6 * h7_p1 - e_6 * fs_16_663_35 * r_8 * h5_p1 - e_6 * fs_14_429_14 * r_10 * h3_p1 - e_6 * fs_6_143_21 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ph11_0, ph11_p2, ab_2, pc_83 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p2 = ph11_p2[k];

        pc_83[k] = - e_0 * fs_14175_32_35 * h1_0 + e_1 * fs_5355_16_35 * h3_0 + e_1 * fs_5355_16_21 * h3_p2 + e_1 * fs_8505_8_35 * r_2 * h1_0 - e_2 * fs_95_1_35 * h5_0 - e_2 * fs_380_1_3 * h5_p2 - e_2 * fs_2975_8_35 * r_2 * h3_0 - e_2 * fs_2975_8_21 * r_2 * h3_p2 - e_2 * fs_6075_8_35 * r_4 * h1_0 + e_3 * fs_1260_143_35 * h7_0 + e_3 * fs_8505_286_30 * h7_p2 + e_3 * fs_760_13_35 * r_2 * h5_0 + e_3 * fs_3040_13_3 * r_2 * h5_p2 + e_3 * fs_2975_22_35 * r_4 * h3_0 + e_3 * fs_2975_22_21 * r_4 * h3_p2 + e_3 * fs_225_1_35 * r_6 * h1_0 + e_4 * fs_4347_2431_35 * h9_0 - e_4 * fs_5796_2431_154 * h9_p2 - e_4 * fs_7560_2431_35 * r_2 * h7_0 - e_4 * fs_25515_2431_30 * r_2 * h7_p2 - e_4 * fs_152_13_35 * r_4 * h5_0 - e_4 * fs_608_13_3 * r_4 * h5_p2 - e_4 * fs_2975_143_35 * r_6 * h3_0 - e_4 * fs_2975_143_21 * r_6 * h3_p2 - e_4 * fs_675_22_35 * r_8 * h1_0 - e_5 * fs_4950_4199_35 * h11_0 + e_5 * fs_450_4199_3003 * h11_p2 - e_5 * fs_828_2431_35 * r_2 * h9_0 + e_5 * fs_1104_2431_154 * r_2 * h9_p2 + e_5 * fs_15120_46189_35 * r_4 * h7_0 + e_5 * fs_51030_46189_30 * r_4 * h7_p2 + e_5 * fs_608_663_35 * r_6 * h5_0 + e_5 * fs_2432_663_3 * r_6 * h5_p2 + e_5 * fs_595_429_35 * r_8 * h3_0 + e_5 * fs_595_429_21 * r_8 * h3_p2 + e_5 * fs_270_143_35 * r_10 * h1_0 + e_6 * fs_396_4199_35 * r_2 * h11_0 - e_6 * fs_36_4199_3003 * r_2 * h11_p2 + e_6 * fs_36_2431_35 * r_4 * h9_0 - e_6 * fs_48_2431_154 * r_4 * h9_p2 - e_6 * fs_480_46189_35 * r_6 * h7_0 - e_6 * fs_1620_46189_30 * r_6 * h7_p2 - e_6 * fs_16_663_35 * r_8 * h5_0 - e_6 * fs_64_663_3 * r_8 * h5_p2 - e_6 * fs_14_429_35 * r_10 * h3_0 - e_6 * fs_14_429_21 * r_10 * h3_p2 - e_6 * fs_6_143_35 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2, pc_84 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p3 = ph11_p3[k];

        pc_84[k] = e_0 * fs_14175_32_5 * h1_p1 + e_1 * fs_5355_32_30 * h3_p1 - e_1 * fs_37485_32_2 * h3_p3 - e_1 * fs_8505_8_5 * r_2 * h1_p1 - e_2 * fs_380_1_3 * h5_p1 + e_2 * fs_95_4_14 * h5_p3 - e_2 * fs_2975_16_30 * r_2 * h3_p1 + e_2 * fs_20825_16_2 * r_2 * h3_p3 + e_2 * fs_6075_8_5 * r_4 * h1_p1 + e_3 * fs_15435_572_35 * h7_p1 + e_3 * fs_2835_572_105 * h7_p3 + e_3 * fs_3040_13_3 * r_2 * h5_p1 - e_3 * fs_190_13_14 * r_2 * h5_p3 + e_3 * fs_2975_44_30 * r_4 * h3_p1 - e_3 * fs_20825_44_2 * r_4 * h3_p3 - e_3 * fs_225_1_5 * r_6 * h1_p1 - e_4 * f_33327_2431 * h9_p1 - e_4 * fs_4347_4862_462 * h9_p3 - e_4 * fs_46305_4862_35 * r_2 * h7_p1 - e_4 * fs_8505_4862_105 * r_2 * h7_p3 - e_4 * fs_608_13_3 * r_4 * h5_p1 + e_4 * fs_38_13_14 * r_4 * h5_p3 - e_4 * fs_2975_286_30 * r_6 * h3_p1 + e_4 * fs_20825_286_2 * r_6 * h3_p3 + e_4 * fs_675_22_5 * r_8 * h1_p1 - e_5 * fs_900_4199_330 * h11_p1 + e_5 * fs_300_4199_6006 * h11_p3 + e_5 * f_6348_2431 * r_2 * h9_p1 + e_5 * fs_414_2431_462 * r_2 * h9_p3 + e_5 * fs_46305_46189_35 * r_4 * h7_p1 + e_5 * fs_8505_46189_105 * r_4 * h7_p3 + e_5 * fs_2432_663_3 * r_6 * h5_p1 - e_5 * fs_152_663_14 * r_6 * h5_p3 + e_5 * fs_595_858_30 * r_8 * h3_p1 - e_5 * fs_4165_858_2 * r_8 * h3_p3 - e_5 * fs_270_143_5 * r_10 * h1_p1 + e_6 * fs_72_4199_330 * r_2 * h11_p1 - e_6 * fs_24_4199_6006 * r_2 * h11_p3 - e_6 * f_276_2431 * r_4 * h9_p1 - e_6 * fs_18_2431_462 * r_4 * h9_p3 - e_6 * fs_1470_46189_35 * r_6 * h7_p1 - e_6 * fs_270_46189_105 * r_6 * h7_p3 - e_6 * fs_64_663_3 * r_8 * h5_p1 + e_6 * fs_4_663_14 * r_8 * h5_p3 - e_6 * fs_7_429_30 * r_10 * h3_p1 + e_6 * fs_49_429_2 * r_10 * h3_p3 + e_6 * fs_6_143_5 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2, pc_85 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_85[k] = - e_1 * fs_5355_4_2 * h3_p2 + e_2 * fs_95_4_14 * h5_p2 + e_2 * fs_1045_8_42 * h5_p4 + e_2 * fs_2975_2_2 * r_2 * h3_p2 + e_3 * fs_6615_286_35 * h7_p2 - e_3 * fs_2835_572_770 * h7_p4 - e_3 * fs_190_13_14 * r_2 * h5_p2 - e_3 * fs_1045_13_42 * r_2 * h5_p4 - e_3 * fs_5950_11_2 * r_4 * h3_p2 - e_4 * fs_25599_4862_33 * h9_p2 + e_4 * fs_483_9724_6006 * h9_p4 - e_4 * fs_19845_2431_35 * r_2 * h7_p2 + e_4 * fs_8505_4862_770 * r_2 * h7_p4 + e_4 * fs_38_13_14 * r_4 * h5_p2 + e_4 * fs_209_13_42 * r_4 * h5_p4 + e_4 * fs_11900_143_2 * r_6 * h3_p2 - e_5 * fs_675_4199_286 * h11_p2 + e_5 * fs_225_8398_30030 * h11_p4 + e_5 * fs_2438_2431_33 * r_2 * h9_p2 - e_5 * fs_23_2431_6006 * r_2 * h9_p4 + e_5 * fs_39690_46189_35 * r_4 * h7_p2 - e_5 * fs_8505_46189_770 * r_4 * h7_p4 - e_5 * fs_152_663_14 * r_6 * h5_p2 - e_5 * fs_836_663_42 * r_6 * h5_p4 - e_5 * fs_2380_429_2 * r_8 * h3_p2 + e_6 * fs_54_4199_286 * r_2 * h11_p2 - e_6 * fs_9_4199_30030 * r_2 * h11_p4 - e_6 * fs_106_2431_33 * r_4 * h9_p2 + e_6 * fs_1_2431_6006 * r_4 * h9_p4 - e_6 * fs_1260_46189_35 * r_6 * h7_p2 + e_6 * fs_270_46189_770 * r_6 * h7_p4 + e_6 * fs_4_663_14 * r_8 * h5_p2 + e_6 * fs_22_663_42 * r_8 * h5_p4 + e_6 * fs_56_429_2 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2, pc_86 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_86[k] = e_1 * fs_5355_16_6 * h3_p3 + e_2 * fs_1045_8_42 * h5_p3 - e_2 * fs_285_8_210 * h5_p5 - e_2 * fs_2975_8_6 * r_2 * h3_p3 - e_3 * fs_2205_286_35 * h7_p3 - e_3 * fs_1260_143_385 * h7_p5 - e_3 * fs_1045_13_42 * r_2 * h5_p3 + e_3 * fs_285_13_210 * r_2 * h5_p5 + e_3 * fs_2975_22_6 * r_4 * h3_p3 - e_4 * fs_27531_9724_154 * h9_p3 + e_4 * fs_4347_9724_4290 * h9_p5 + e_4 * fs_6615_2431_35 * r_2 * h7_p3 + e_4 * fs_7560_2431_385 * r_2 * h7_p5 + e_4 * fs_209_13_42 * r_4 * h5_p3 - e_4 * fs_57_13_210 * r_4 * h5_p5 - e_4 * fs_2975_143_6 * r_6 * h3_p3 - e_5 * fs_150_4199_2002 * h11_p3 + e_5 * fs_300_4199_2145 * h11_p5 + e_5 * fs_1311_2431_154 * r_2 * h9_p3 - e_5 * fs_207_2431_4290 * r_2 * h9_p5 - e_5 * fs_13230_46189_35 * r_4 * h7_p3 - e_5 * fs_15120_46189_385 * r_4 * h7_p5 - e_5 * fs_836_663_42 * r_6 * h5_p3 + e_5 * fs_76_221_210 * r_6 * h5_p5 + e_5 * fs_595_429_6 * r_8 * h3_p3 + e_6 * fs_12_4199_2002 * r_2 * h11_p3 - e_6 * fs_24_4199_2145 * r_2 * h11_p5 - e_6 * fs_57_2431_154 * r_4 * h9_p3 + e_6 * fs_9_2431_4290 * r_4 * h9_p5 + e_6 * fs_420_46189_35 * r_6 * h7_p3 + e_6 * fs_480_46189_385 * r_6 * h7_p5 + e_6 * fs_22_663_42 * r_8 * h5_p3 - e_6 * fs_2_221_210 * r_8 * h5_p5 - e_6 * fs_14_429_6 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2, pc_87 : simd::cache_line_size())
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

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_87[k] = - e_2 * fs_285_8_210 * h5_p4 - e_3 * fs_11025_572_154 * h7_p4 + e_3 * fs_1575_286_1001 * h7_p6 + e_3 * fs_285_13_210 * r_2 * h5_p4 - e_4 * fs_1449_9724_30030 * h9_p4 + e_4 * fs_7245_4862_715 * h9_p6 + e_4 * fs_33075_4862_154 * r_2 * h7_p4 - e_4 * fs_4725_2431_1001 * r_2 * h7_p6 - e_4 * fs_57_13_210 * r_4 * h5_p4 - e_5 * fs_75_8398_6006 * h11_p4 + e_5 * fs_150_4199_2431 * h11_p6 + e_5 * fs_69_2431_30030 * r_2 * h9_p4 - e_5 * fs_690_2431_715 * r_2 * h9_p6 - e_5 * fs_33075_46189_154 * r_4 * h7_p4 + e_5 * fs_9450_46189_1001 * r_4 * h7_p6 + e_5 * fs_76_221_210 * r_6 * h5_p4 + e_6 * fs_3_4199_6006 * r_2 * h11_p4 - e_6 * fs_12_4199_2431 * r_2 * h11_p6 - e_6 * fs_3_2431_30030 * r_4 * h9_p4 + e_6 * fs_30_2431_715 * r_4 * h9_p6 + e_6 * fs_1050_46189_154 * r_6 * h7_p4 - e_6 * fs_300_46189_1001 * r_6 * h7_p6 - e_6 * fs_2_221_210 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2, pc_88 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_88[k] = e_1 * fs_5355_32_6 * h3_m3 + e_2 * fs_475_4_42 * h5_m3 - e_2 * fs_2975_16_6 * r_2 * h3_m3 + e_3 * fs_2205_572_715 * h7_m7 + e_3 * fs_11025_286_35 * h7_m3 - e_3 * fs_950_13_42 * r_2 * h5_m3 + e_3 * fs_2975_44_6 * r_4 * h3_m3 + e_4 * fs_7245_4862_858 * h9_m7 + e_4 * fs_7245_4862_154 * h9_m3 - e_4 * fs_6615_4862_715 * r_2 * h7_m7 - e_4 * fs_33075_2431_35 * r_2 * h7_m3 + e_4 * fs_190_13_42 * r_4 * h5_m3 - e_4 * fs_2975_286_6 * r_6 * h3_m3 + e_5 * fs_225_4199_2431 * h11_m7 + e_5 * fs_75_8398_2002 * h11_m3 - e_5 * fs_690_2431_858 * r_2 * h9_m7 - e_5 * fs_690_2431_154 * r_2 * h9_m3 + e_5 * fs_6615_46189_715 * r_4 * h7_m7 + e_5 * fs_66150_46189_35 * r_4 * h7_m3 - e_5 * fs_760_663_42 * r_6 * h5_m3 + e_5 * fs_595_858_6 * r_8 * h3_m3 - e_6 * fs_18_4199_2431 * r_2 * h11_m7 - e_6 * fs_3_4199_2002 * r_2 * h11_m3 + e_6 * fs_30_2431_858 * r_4 * h9_m7 + e_6 * fs_30_2431_154 * r_4 * h9_m3 - e_6 * fs_210_46189_715 * r_6 * h7_m7 - e_6 * fs_2100_46189_35 * r_6 * h7_m3 + e_6 * fs_20_663_42 * r_8 * h5_m3 - e_6 * fs_7_429_6 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2, pc_89 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_89[k] = - e_1 * fs_16065_32_10 * h3_m2 - e_2 * fs_285_4_70 * h5_m2 + e_2 * fs_8925_16_10 * r_2 * h3_m2 - e_3 * fs_315_44_1001 * h7_m6 + e_3 * fs_15435_286_7 * h7_m2 + e_3 * fs_570_13_70 * r_2 * h5_m2 - e_3 * fs_8925_44_10 * r_4 * h3_m2 + e_4 * fs_1449_2431_715 * h9_m6 + e_4 * fs_5796_2431_165 * h9_m2 + e_4 * fs_945_374_1001 * r_2 * h7_m6 - e_4 * fs_46305_2431_7 * r_2 * h7_m2 - e_4 * fs_114_13_70 * r_4 * h5_m2 + e_4 * fs_8925_286_10 * r_6 * h3_m2 + e_5 * fs_375_4199_2431 * h11_m6 + e_5 * fs_225_8398_1430 * h11_m2 - e_5 * fs_276_2431_715 * r_2 * h9_m6 - e_5 * fs_1104_2431_165 * r_2 * h9_m2 - e_5 * fs_945_3553_1001 * r_4 * h7_m6 + e_5 * fs_92610_46189_7 * r_4 * h7_m2 + e_5 * fs_152_221_70 * r_6 * h5_m2 - e_5 * fs_595_286_10 * r_8 * h3_m2 - e_6 * fs_30_4199_2431 * r_2 * h11_m6 - e_6 * fs_9_4199_1430 * r_2 * h11_m2 + e_6 * fs_12_2431_715 * r_4 * h9_m6 + e_6 * fs_48_2431_165 * r_4 * h9_m2 + e_6 * fs_30_3553_1001 * r_6 * h7_m6 - e_6 * fs_2940_46189_7 * r_6 * h7_m2 - e_6 * fs_4_221_70 * r_8 * h5_m2 + e_6 * fs_7_143_10 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2, pc_90 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m1 = ph11_m1[k];

        pc_90[k] = e_0 * fs_14175_32_3 * h1_m1 + e_1 * fs_37485_32_2 * h3_m1 - e_1 * fs_8505_8_3 * r_2 * h1_m1 + e_2 * fs_475_4_42 * h5_m5 - e_2 * fs_1045_4_5 * h5_m1 - e_2 * fs_20825_16_2 * r_2 * h3_m1 + e_2 * fs_6075_8_3 * r_4 * h1_m1 - e_3 * fs_945_572_77 * h7_m5 - e_3 * fs_2205_286_21 * h7_m1 - e_3 * fs_950_13_42 * r_2 * h5_m5 + e_3 * fs_2090_13_5 * r_2 * h5_m1 + e_3 * fs_20825_44_2 * r_4 * h3_m1 - e_3 * fs_225_1_3 * r_6 * h1_m1 - e_4 * fs_2415_4862_858 * h9_m5 + e_4 * fs_42987_4862_15 * h9_m1 + e_4 * fs_2835_4862_77 * r_2 * h7_m5 + e_4 * fs_6615_2431_21 * r_2 * h7_m1 + e_4 * fs_190_13_42 * r_4 * h5_m5 - e_4 * fs_418_13_5 * r_4 * h5_m1 - e_4 * fs_20825_286_2 * r_6 * h3_m1 + e_4 * fs_675_22_3 * r_8 * h1_m1 + e_5 * fs_1125_4199_429 * h11_m5 + e_5 * fs_3375_8398_22 * h11_m1 + e_5 * fs_230_2431_858 * r_2 * h9_m5 - e_5 * fs_4094_2431_15 * r_2 * h9_m1 - e_5 * fs_2835_46189_77 * r_4 * h7_m5 - e_5 * fs_13230_46189_21 * r_4 * h7_m1 - e_5 * fs_760_663_42 * r_6 * h5_m5 + e_5 * fs_1672_663_5 * r_6 * h5_m1 + e_5 * fs_4165_858_2 * r_8 * h3_m1 - e_5 * fs_270_143_3 * r_10 * h1_m1 - e_6 * fs_90_4199_429 * r_2 * h11_m5 - e_6 * fs_135_4199_22 * r_2 * h11_m1 - e_6 * fs_10_2431_858 * r_4 * h9_m5 + e_6 * fs_178_2431_15 * r_4 * h9_m1 + e_6 * fs_90_46189_77 * r_6 * h7_m5 + e_6 * fs_420_46189_21 * r_6 * h7_m1 + e_6 * fs_20_663_42 * r_8 * h5_m5 - e_6 * fs_44_663_5 * r_8 * h5_m1 - e_6 * fs_49_429_2 * r_10 * h3_m1 + e_6 * fs_6_143_3 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2, pc_91 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_91[k] = - e_2 * fs_285_4_70 * h5_m4 + e_3 * fs_8505_1144_462 * h7_m4 + e_3 * fs_570_13_70 * r_2 * h5_m4 - e_4 * fs_1449_4862_10010 * h9_m4 - e_4 * fs_25515_9724_462 * r_2 * h7_m4 - e_4 * fs_114_13_70 * r_4 * h5_m4 + e_5 * fs_1125_8398_2002 * h11_m4 + e_5 * fs_138_2431_10010 * r_2 * h9_m4 + e_5 * fs_25515_92378_462 * r_4 * h7_m4 + e_5 * fs_152_221_70 * r_6 * h5_m4 - e_6 * fs_45_4199_2002 * r_2 * h11_m4 - e_6 * fs_6_2431_10010 * r_4 * h9_m4 - e_6 * fs_405_46189_462 * r_6 * h7_m4 - e_6 * fs_4_221_70 * r_8 * h5_m4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2, pc_92 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m1 = ph11_m1[k];

        pc_92[k] = e_0 * fs_14175_32_14 * h1_m1 + e_1 * fs_5355_16_35 * h3_m3 - e_1 * fs_5355_16_21 * h3_m1 - e_1 * fs_8505_8_14 * r_2 * h1_m1 - e_2 * fs_1045_4_5 * h5_m3 + e_2 * fs_95_2_210 * h5_m1 - e_2 * fs_2975_8_35 * r_2 * h3_m3 + e_2 * fs_2975_8_21 * r_2 * h3_m1 + e_2 * fs_6075_8_14 * r_4 * h1_m1 + e_3 * fs_65205_1144_6 * h7_m3 - e_3 * fs_115605_1144_2 * h7_m1 + e_3 * fs_2090_13_5 * r_2 * h5_m3 - e_3 * fs_380_13_210 * r_2 * h5_m1 + e_3 * fs_2975_22_35 * r_4 * h3_m3 - e_3 * fs_2975_22_21 * r_4 * h3_m1 - e_3 * fs_225_1_14 * r_6 * h1_m1 - e_4 * fs_10143_4862_165 * h9_m3 + e_4 * fs_2898_2431_70 * h9_m1 - e_4 * fs_195615_9724_6 * r_2 * h7_m3 + e_4 * fs_346815_9724_2 * r_2 * h7_m1 - e_4 * fs_418_13_5 * r_4 * h5_m3 + e_4 * fs_76_13_210 * r_4 * h5_m1 - e_4 * fs_2975_143_35 * r_6 * h3_m3 + e_4 * fs_2975_143_21 * r_6 * h3_m1 + e_4 * fs_675_22_14 * r_8 * h1_m1 + e_5 * fs_525_4199_2145 * h11_m3 + e_5 * fs_1125_4199_231 * h11_m1 + e_5 * fs_966_2431_165 * r_2 * h9_m3 - e_5 * fs_552_2431_70 * r_2 * h9_m1 + e_5 * fs_195615_92378_6 * r_4 * h7_m3 - e_5 * fs_346815_92378_2 * r_4 * h7_m1 + e_5 * fs_1672_663_5 * r_6 * h5_m3 - e_5 * fs_304_663_210 * r_6 * h5_m1 + e_5 * fs_595_429_35 * r_8 * h3_m3 - e_5 * fs_595_429_21 * r_8 * h3_m1 - e_5 * fs_270_143_14 * r_10 * h1_m1 - e_6 * fs_42_4199_2145 * r_2 * h11_m3 - e_6 * fs_90_4199_231 * r_2 * h11_m1 - e_6 * fs_42_2431_165 * r_4 * h9_m3 + e_6 * fs_24_2431_70 * r_4 * h9_m1 - e_6 * fs_3105_46189_6 * r_6 * h7_m3 + e_6 * fs_5505_46189_2 * r_6 * h7_m1 - e_6 * fs_44_663_5 * r_8 * h5_m3 + e_6 * fs_8_663_210 * r_8 * h5_m1 - e_6 * fs_14_429_35 * r_10 * h3_m3 + e_6 * fs_14_429_21 * r_10 * h3_m1 + e_6 * fs_6_143_14 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ab_2, pc_93 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];

        pc_93[k] = - e_1 * fs_5355_16_14 * h3_p2 + e_2 * fs_475_2_2 * h5_p2 + e_2 * fs_2975_8_14 * r_2 * h3_p2 - e_3 * fs_2835_572_5 * h7_p2 - e_3 * fs_1900_13_2 * r_2 * h5_p2 - e_3 * fs_2975_22_14 * r_4 * h3_p2 - e_4 * fs_2415_2431_231 * h9_p2 + e_4 * fs_8505_4862_5 * r_2 * h7_p2 + e_4 * fs_380_13_2 * r_4 * h5_p2 + e_4 * fs_2975_143_14 * r_6 * h3_p2 + e_5 * fs_675_4199_2002 * h11_p2 + e_5 * fs_460_2431_231 * r_2 * h9_p2 - e_5 * fs_8505_46189_5 * r_4 * h7_p2 - e_5 * fs_1520_663_2 * r_6 * h5_p2 - e_5 * fs_595_429_14 * r_8 * h3_p2 - e_6 * fs_54_4199_2002 * r_2 * h11_p2 - e_6 * fs_20_2431_231 * r_4 * h9_p2 + e_6 * fs_270_46189_5 * r_6 * h7_p2 + e_6 * fs_40_663_2 * r_8 * h5_p2 + e_6 * fs_14_429_14 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2, pc_94 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p3 = ph11_p3[k];

        pc_94[k] = - e_0 * fs_14175_32_14 * h1_p1 + e_1 * fs_5355_16_21 * h3_p1 + e_1 * fs_5355_16_35 * h3_p3 + e_1 * fs_8505_8_14 * r_2 * h1_p1 - e_2 * fs_95_2_210 * h5_p1 - e_2 * fs_1045_4_5 * h5_p3 - e_2 * fs_2975_8_21 * r_2 * h3_p1 - e_2 * fs_2975_8_35 * r_2 * h3_p3 - e_2 * fs_6075_8_14 * r_4 * h1_p1 + e_3 * fs_115605_1144_2 * h7_p1 + e_3 * fs_65205_1144_6 * h7_p3 + e_3 * fs_380_13_210 * r_2 * h5_p1 + e_3 * fs_2090_13_5 * r_2 * h5_p3 + e_3 * fs_2975_22_21 * r_4 * h3_p1 + e_3 * fs_2975_22_35 * r_4 * h3_p3 + e_3 * fs_225_1_14 * r_6 * h1_p1 - e_4 * fs_2898_2431_70 * h9_p1 - e_4 * fs_10143_4862_165 * h9_p3 - e_4 * fs_346815_9724_2 * r_2 * h7_p1 - e_4 * fs_195615_9724_6 * r_2 * h7_p3 - e_4 * fs_76_13_210 * r_4 * h5_p1 - e_4 * fs_418_13_5 * r_4 * h5_p3 - e_4 * fs_2975_143_21 * r_6 * h3_p1 - e_4 * fs_2975_143_35 * r_6 * h3_p3 - e_4 * fs_675_22_14 * r_8 * h1_p1 - e_5 * fs_1125_4199_231 * h11_p1 + e_5 * fs_525_4199_2145 * h11_p3 + e_5 * fs_552_2431_70 * r_2 * h9_p1 + e_5 * fs_966_2431_165 * r_2 * h9_p3 + e_5 * fs_346815_92378_2 * r_4 * h7_p1 + e_5 * fs_195615_92378_6 * r_4 * h7_p3 + e_5 * fs_304_663_210 * r_6 * h5_p1 + e_5 * fs_1672_663_5 * r_6 * h5_p3 + e_5 * fs_595_429_21 * r_8 * h3_p1 + e_5 * fs_595_429_35 * r_8 * h3_p3 + e_5 * fs_270_143_14 * r_10 * h1_p1 + e_6 * fs_90_4199_231 * r_2 * h11_p1 - e_6 * fs_42_4199_2145 * r_2 * h11_p3 - e_6 * fs_24_2431_70 * r_4 * h9_p1 - e_6 * fs_42_2431_165 * r_4 * h9_p3 - e_6 * fs_5505_46189_2 * r_6 * h7_p1 - e_6 * fs_3105_46189_6 * r_6 * h7_p3 - e_6 * fs_8_663_210 * r_8 * h5_p1 - e_6 * fs_44_663_5 * r_8 * h5_p3 - e_6 * fs_14_429_21 * r_10 * h3_p1 - e_6 * fs_14_429_35 * r_10 * h3_p3 - e_6 * fs_6_143_14 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ph11_0, ph11_p4, ab_2, pc_95 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p4 = ph11_p4[k];

        pc_95[k] = - e_0 * fs_14175_8_2 * h1_0 + e_1 * fs_5355_8_2 * h3_0 + e_1 * fs_8505_2_2 * r_2 * h1_0 + e_2 * fs_475_2_2 * h5_0 - e_2 * fs_285_4_70 * h5_p4 - e_2 * fs_2975_4_2 * r_2 * h3_0 - e_2 * fs_6075_2_2 * r_4 * h1_0 - e_3 * fs_90405_572_2 * h7_0 + e_3 * fs_8505_1144_462 * h7_p4 - e_3 * fs_1900_13_2 * r_2 * h5_0 + e_3 * fs_570_13_70 * r_2 * h5_p4 + e_3 * fs_2975_11_2 * r_4 * h3_0 + e_3 * fs_900_1_2 * r_6 * h1_0 + e_4 * fs_65205_2431_2 * h9_0 - e_4 * fs_1449_4862_10010 * h9_p4 + e_4 * fs_271215_4862_2 * r_2 * h7_0 - e_4 * fs_25515_9724_462 * r_2 * h7_p4 + e_4 * fs_380_13_2 * r_4 * h5_0 - e_4 * fs_114_13_70 * r_4 * h5_p4 - e_4 * fs_5950_143_2 * r_6 * h3_0 - e_4 * fs_1350_11_2 * r_8 * h1_0 + e_5 * fs_12375_4199_2 * h11_0 + e_5 * fs_1125_8398_2002 * h11_p4 - e_5 * fs_12420_2431_2 * r_2 * h9_0 + e_5 * fs_138_2431_10010 * r_2 * h9_p4 - e_5 * fs_271215_46189_2 * r_4 * h7_0 + e_5 * fs_25515_92378_462 * r_4 * h7_p4 - e_5 * fs_1520_663_2 * r_6 * h5_0 + e_5 * fs_152_221_70 * r_6 * h5_p4 + e_5 * fs_1190_429_2 * r_8 * h3_0 + e_5 * fs_1080_143_2 * r_10 * h1_0 - e_6 * fs_990_4199_2 * r_2 * h11_0 - e_6 * fs_45_4199_2002 * r_2 * h11_p4 + e_6 * fs_540_2431_2 * r_4 * h9_0 - e_6 * fs_6_2431_10010 * r_4 * h9_p4 + e_6 * fs_8610_46189_2 * r_6 * h7_0 - e_6 * fs_405_46189_462 * r_6 * h7_p4 + e_6 * fs_40_663_2 * r_8 * h5_0 - e_6 * fs_4_221_70 * r_8 * h5_p4 - e_6 * fs_28_429_2 * r_10 * h3_0 - e_6 * fs_24_143_2 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2, pc_96 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p5 = ph11_p5[k];

        pc_96[k] = e_0 * fs_14175_32_3 * h1_p1 + e_1 * fs_37485_32_2 * h3_p1 - e_1 * fs_8505_8_3 * r_2 * h1_p1 - e_2 * fs_1045_4_5 * h5_p1 + e_2 * fs_475_4_42 * h5_p5 - e_2 * fs_20825_16_2 * r_2 * h3_p1 + e_2 * fs_6075_8_3 * r_4 * h1_p1 - e_3 * fs_2205_286_21 * h7_p1 - e_3 * fs_945_572_77 * h7_p5 + e_3 * fs_2090_13_5 * r_2 * h5_p1 - e_3 * fs_950_13_42 * r_2 * h5_p5 + e_3 * fs_20825_44_2 * r_4 * h3_p1 - e_3 * fs_225_1_3 * r_6 * h1_p1 + e_4 * fs_42987_4862_15 * h9_p1 - e_4 * fs_2415_4862_858 * h9_p5 + e_4 * fs_6615_2431_21 * r_2 * h7_p1 + e_4 * fs_2835_4862_77 * r_2 * h7_p5 - e_4 * fs_418_13_5 * r_4 * h5_p1 + e_4 * fs_190_13_42 * r_4 * h5_p5 - e_4 * fs_20825_286_2 * r_6 * h3_p1 + e_4 * fs_675_22_3 * r_8 * h1_p1 + e_5 * fs_3375_8398_22 * h11_p1 + e_5 * fs_1125_4199_429 * h11_p5 - e_5 * fs_4094_2431_15 * r_2 * h9_p1 + e_5 * fs_230_2431_858 * r_2 * h9_p5 - e_5 * fs_13230_46189_21 * r_4 * h7_p1 - e_5 * fs_2835_46189_77 * r_4 * h7_p5 + e_5 * fs_1672_663_5 * r_6 * h5_p1 - e_5 * fs_760_663_42 * r_6 * h5_p5 + e_5 * fs_4165_858_2 * r_8 * h3_p1 - e_5 * fs_270_143_3 * r_10 * h1_p1 - e_6 * fs_135_4199_22 * r_2 * h11_p1 - e_6 * fs_90_4199_429 * r_2 * h11_p5 + e_6 * fs_178_2431_15 * r_4 * h9_p1 - e_6 * fs_10_2431_858 * r_4 * h9_p5 + e_6 * fs_420_46189_21 * r_6 * h7_p1 + e_6 * fs_90_46189_77 * r_6 * h7_p5 - e_6 * fs_44_663_5 * r_8 * h5_p1 + e_6 * fs_20_663_42 * r_8 * h5_p5 - e_6 * fs_49_429_2 * r_10 * h3_p1 + e_6 * fs_6_143_3 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2, pc_97 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_97[k] = - e_1 * fs_16065_32_10 * h3_p2 - e_2 * fs_285_4_70 * h5_p2 + e_2 * fs_8925_16_10 * r_2 * h3_p2 + e_3 * fs_15435_286_7 * h7_p2 - e_3 * fs_315_44_1001 * h7_p6 + e_3 * fs_570_13_70 * r_2 * h5_p2 - e_3 * fs_8925_44_10 * r_4 * h3_p2 + e_4 * fs_5796_2431_165 * h9_p2 + e_4 * fs_1449_2431_715 * h9_p6 - e_4 * fs_46305_2431_7 * r_2 * h7_p2 + e_4 * fs_945_374_1001 * r_2 * h7_p6 - e_4 * fs_114_13_70 * r_4 * h5_p2 + e_4 * fs_8925_286_10 * r_6 * h3_p2 + e_5 * fs_225_8398_1430 * h11_p2 + e_5 * fs_375_4199_2431 * h11_p6 - e_5 * fs_1104_2431_165 * r_2 * h9_p2 - e_5 * fs_276_2431_715 * r_2 * h9_p6 + e_5 * fs_92610_46189_7 * r_4 * h7_p2 - e_5 * fs_945_3553_1001 * r_4 * h7_p6 + e_5 * fs_152_221_70 * r_6 * h5_p2 - e_5 * fs_595_286_10 * r_8 * h3_p2 - e_6 * fs_9_4199_1430 * r_2 * h11_p2 - e_6 * fs_30_4199_2431 * r_2 * h11_p6 + e_6 * fs_48_2431_165 * r_4 * h9_p2 + e_6 * fs_12_2431_715 * r_4 * h9_p6 - e_6 * fs_2940_46189_7 * r_6 * h7_p2 + e_6 * fs_30_3553_1001 * r_6 * h7_p6 - e_6 * fs_4_221_70 * r_8 * h5_p2 + e_6 * fs_7_143_10 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2, pc_98 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_98[k] = e_1 * fs_5355_32_6 * h3_p3 + e_2 * fs_475_4_42 * h5_p3 - e_2 * fs_2975_16_6 * r_2 * h3_p3 + e_3 * fs_11025_286_35 * h7_p3 + e_3 * fs_2205_572_715 * h7_p7 - e_3 * fs_950_13_42 * r_2 * h5_p3 + e_3 * fs_2975_44_6 * r_4 * h3_p3 + e_4 * fs_7245_4862_154 * h9_p3 + e_4 * fs_7245_4862_858 * h9_p7 - e_4 * fs_33075_2431_35 * r_2 * h7_p3 - e_4 * fs_6615_4862_715 * r_2 * h7_p7 + e_4 * fs_190_13_42 * r_4 * h5_p3 - e_4 * fs_2975_286_6 * r_6 * h3_p3 + e_5 * fs_75_8398_2002 * h11_p3 + e_5 * fs_225_4199_2431 * h11_p7 - e_5 * fs_690_2431_154 * r_2 * h9_p3 - e_5 * fs_690_2431_858 * r_2 * h9_p7 + e_5 * fs_66150_46189_35 * r_4 * h7_p3 + e_5 * fs_6615_46189_715 * r_4 * h7_p7 - e_5 * fs_760_663_42 * r_6 * h5_p3 + e_5 * fs_595_858_6 * r_8 * h3_p3 - e_6 * fs_3_4199_2002 * r_2 * h11_p3 - e_6 * fs_18_4199_2431 * r_2 * h11_p7 + e_6 * fs_30_2431_154 * r_4 * h9_p3 + e_6 * fs_30_2431_858 * r_4 * h9_p7 - e_6 * fs_2100_46189_35 * r_6 * h7_p3 - e_6 * fs_210_46189_715 * r_6 * h7_p7 + e_6 * fs_20_663_42 * r_8 * h5_p3 - e_6 * fs_7_429_6 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2, pc_99 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_99[k] = - e_1 * f_16065_16 * h3_m2 - e_2 * fs_1425_4_7 * h5_m2 + e_2 * f_8925_8 * r_2 * h3_m2 - e_3 * fs_6615_286_70 * h7_m2 + e_3 * fs_2850_13_7 * r_2 * h5_m2 - e_3 * f_8925_22 * r_4 * h3_m2 + e_4 * fs_2415_4862_7293 * h9_m8 - e_4 * fs_7245_4862_66 * h9_m2 + e_4 * fs_19845_2431_70 * r_2 * h7_m2 - e_4 * fs_570_13_7 * r_4 * h5_m2 + e_4 * f_8925_143 * r_6 * h3_m2 + e_5 * fs_75_4199_46189 * h11_m8 - e_5 * fs_75_4199_143 * h11_m2 - e_5 * fs_230_2431_7293 * r_2 * h9_m8 + e_5 * fs_690_2431_66 * r_2 * h9_m2 - e_5 * fs_39690_46189_70 * r_4 * h7_m2 + e_5 * fs_760_221_7 * r_6 * h5_m2 - e_5 * f_595_143 * r_8 * h3_m2 - e_6 * fs_6_4199_46189 * r_2 * h11_m8 + e_6 * fs_6_4199_143 * r_2 * h11_m2 + e_6 * fs_10_2431_7293 * r_4 * h9_m8 - e_6 * fs_30_2431_66 * r_4 * h9_m2 + e_6 * fs_1260_46189_70 * r_6 * h7_m2 - e_6 * fs_20_221_7 * r_8 * h5_m2 + e_6 * f_14_143 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2, pc_100 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_100[k] = e_0 * fs_14175_64_6 * h1_m1 + e_1 * f_16065_8 * h3_m1 - e_1 * fs_8505_16_6 * r_2 * h1_m1 + e_2 * fs_285_8_10 * h5_m1 - e_2 * f_8925_4 * r_2 * h3_m1 + e_2 * fs_6075_16_6 * r_4 * h1_m1 - e_3 * fs_6615_572_286 * h7_m7 - e_3 * fs_4410_143_42 * h7_m1 - e_3 * fs_285_13_10 * r_2 * h5_m1 + e_3 * f_8925_11 * r_4 * h3_m1 - e_3 * fs_225_2_6 * r_6 * h1_m1 - e_4 * fs_483_4862_2145 * h9_m7 - e_4 * fs_42021_9724_30 * h9_m1 + e_4 * fs_19845_4862_286 * r_2 * h7_m7 + e_4 * fs_26460_2431_42 * r_2 * h7_m1 + e_4 * fs_57_13_10 * r_4 * h5_m1 - e_4 * f_17850_143 * r_6 * h3_m1 + e_4 * fs_675_44_6 * r_8 * h1_m1 + e_5 * fs_150_4199_24310 * h11_m7 - e_5 * fs_750_4199_11 * h11_m1 + e_5 * fs_46_2431_2145 * r_2 * h9_m7 + e_5 * fs_2001_2431_30 * r_2 * h9_m1 - e_5 * fs_19845_46189_286 * r_4 * h7_m7 - e_5 * fs_52920_46189_42 * r_4 * h7_m1 - e_5 * fs_76_221_10 * r_6 * h5_m1 + e_5 * f_1190_143 * r_8 * h3_m1 - e_5 * fs_135_143_6 * r_10 * h1_m1 - e_6 * fs_12_4199_24310 * r_2 * h11_m7 + e_6 * fs_60_4199_11 * r_2 * h11_m1 - e_6 * fs_2_2431_2145 * r_4 * h9_m7 - e_6 * fs_87_2431_30 * r_4 * h9_m1 + e_6 * fs_630_46189_286 * r_6 * h7_m7 + e_6 * fs_1680_46189_42 * r_6 * h7_m1 + e_6 * fs_2_221_10 * r_8 * h5_m1 - e_6 * f_28_143 * r_10 * h3_m1 + e_6 * fs_3_143_6 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, ph7_m6, ph9_m6, ph11_m6, ab_2, pc_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h7_m6 = ph7_m6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h11_m6 = ph11_m6[k];

        pc_101[k] = e_3 * fs_945_286_2002 * h7_m6 - e_4 * fs_1932_2431_1430 * h9_m6 - e_4 * fs_2835_2431_2002 * r_2 * h7_m6 + e_5 * fs_375_4199_4862 * h11_m6 + e_5 * fs_368_2431_1430 * r_2 * h9_m6 + e_5 * fs_5670_46189_2002 * r_4 * h7_m6 - e_6 * fs_30_4199_4862 * r_2 * h11_m6 - e_6 * fs_16_2431_1430 * r_4 * h9_m6 - e_6 * fs_180_46189_2002 * r_6 * h7_m6;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2, pc_102 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m1 = ph11_m1[k];

        pc_102[k] = e_0 * fs_42525_32_2 * h1_m1 - e_1 * fs_16065_16_3 * h3_m1 - e_1 * fs_25515_8_2 * r_2 * h1_m1 - e_2 * fs_1425_4_7 * h5_m5 + e_2 * fs_285_4_30 * h5_m1 + e_2 * fs_8925_8_3 * r_2 * h3_m1 + e_2 * fs_18225_8_2 * r_4 * h1_m1 + e_3 * fs_9135_1144_462 * h7_m5 + e_3 * fs_21735_1144_14 * h7_m1 + e_3 * fs_2850_13_7 * r_2 * h5_m5 - e_3 * fs_570_13_30 * r_2 * h5_m1 - e_3 * fs_8925_22_3 * r_4 * h3_m1 - e_3 * fs_675_1_2 * r_6 * h1_m1 - e_4 * fs_12075_4862_143 * h9_m5 - e_4 * fs_51681_4862_10 * h9_m1 - e_4 * fs_27405_9724_462 * r_2 * h7_m5 - e_4 * fs_65205_9724_14 * r_2 * h7_m1 - e_4 * fs_570_13_7 * r_4 * h5_m5 + e_4 * fs_114_13_30 * r_4 * h5_m1 + e_4 * fs_8925_143_3 * r_6 * h3_m1 + e_4 * fs_2025_22_2 * r_8 * h1_m1 + e_5 * fs_1500_4199_286 * h11_m5 - e_5 * fs_1500_4199_33 * h11_m1 + e_5 * fs_1150_2431_143 * r_2 * h9_m5 + e_5 * fs_4922_2431_10 * r_2 * h9_m1 + e_5 * fs_27405_92378_462 * r_4 * h7_m5 + e_5 * fs_65205_92378_14 * r_4 * h7_m1 + e_5 * fs_760_221_7 * r_6 * h5_m5 - e_5 * fs_152_221_30 * r_6 * h5_m1 - e_5 * fs_595_143_3 * r_8 * h3_m1 - e_5 * fs_810_143_2 * r_10 * h1_m1 - e_6 * fs_120_4199_286 * r_2 * h11_m5 + e_6 * fs_120_4199_33 * r_2 * h11_m1 - e_6 * fs_50_2431_143 * r_4 * h9_m5 - e_6 * fs_214_2431_10 * r_4 * h9_m1 - e_6 * fs_435_46189_462 * r_6 * h7_m5 - e_6 * fs_1035_46189_14 * r_6 * h7_m1 - e_6 * fs_20_221_7 * r_8 * h5_m5 + e_6 * fs_4_221_30 * r_8 * h5_m1 + e_6 * fs_14_143_3 * r_10 * h3_m1 + e_6 * fs_18_143_2 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2, pc_103 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_103[k] = e_2 * fs_285_8_10 * h5_m4 + e_2 * fs_285_4_30 * h5_m2 + e_3 * fs_945_572_66 * h7_m4 - e_3 * fs_1260_13_3 * h7_m2 - e_3 * fs_285_13_10 * r_2 * h5_m4 - e_3 * fs_570_13_30 * r_2 * h5_m2 - e_4 * fs_3381_9724_1430 * h9_m4 + e_4 * fs_483_374_385 * h9_m2 - e_4 * fs_2835_4862_66 * r_2 * h7_m4 + e_4 * fs_7560_221_3 * r_2 * h7_m2 + e_4 * fs_57_13_10 * r_4 * h5_m4 + e_4 * fs_114_13_30 * r_4 * h5_m2 + e_5 * fs_2625_8398_286 * h11_m4 + e_5 * fs_75_4199_30030 * h11_m2 + e_5 * fs_161_2431_1430 * r_2 * h9_m4 - e_5 * fs_46_187_385 * r_2 * h9_m2 + e_5 * fs_2835_46189_66 * r_4 * h7_m4 - e_5 * fs_15120_4199_3 * r_4 * h7_m2 - e_5 * fs_76_221_10 * r_6 * h5_m4 - e_5 * fs_152_221_30 * r_6 * h5_m2 - e_6 * fs_105_4199_286 * r_2 * h11_m4 - e_6 * fs_6_4199_30030 * r_2 * h11_m2 - e_6 * fs_7_2431_1430 * r_4 * h9_m4 + e_6 * fs_2_187_385 * r_4 * h9_m2 - e_6 * fs_90_46189_66 * r_6 * h7_m4 + e_6 * fs_480_4199_3 * r_6 * h7_m2 + e_6 * fs_2_221_10 * r_8 * h5_m4 + e_6 * fs_4_221_30 * r_8 * h5_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph11_p3, ab_2, pc_104 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p3 = ph11_p3[k];

        pc_104[k] = - e_1 * fs_5355_8_21 * h3_p3 + e_2 * fs_2375_4_3 * h5_p3 + e_2 * fs_2975_4_21 * r_2 * h3_p3 - e_3 * fs_34965_572_10 * h7_p3 - e_3 * fs_4750_13_3 * r_2 * h5_p3 - e_3 * fs_2975_11_21 * r_4 * h3_p3 + e_4 * fs_16905_4862_11 * h9_p3 + e_4 * fs_104895_4862_10 * r_2 * h7_p3 + e_4 * fs_950_13_3 * r_4 * h5_p3 + e_4 * fs_5950_143_21 * r_6 * h3_p3 + e_5 * fs_2100_4199_143 * h11_p3 - e_5 * fs_1610_2431_11 * r_2 * h9_p3 - e_5 * fs_104895_46189_10 * r_4 * h7_p3 - e_5 * fs_3800_663_3 * r_6 * h5_p3 - e_5 * fs_1190_429_21 * r_8 * h3_p3 - e_6 * fs_168_4199_143 * r_2 * h11_p3 + e_6 * fs_70_2431_11 * r_4 * h9_p3 + e_6 * fs_3330_46189_10 * r_6 * h7_p3 + e_6 * fs_100_663_3 * r_8 * h5_p3 + e_6 * fs_28_429_21 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2, pc_105 : simd::cache_line_size())
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

        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_105[k] = - e_2 * fs_285_4_30 * h5_p2 + e_2 * fs_285_8_10 * h5_p4 + e_3 * fs_1260_13_3 * h7_p2 + e_3 * fs_945_572_66 * h7_p4 + e_3 * fs_570_13_30 * r_2 * h5_p2 - e_3 * fs_285_13_10 * r_2 * h5_p4 - e_4 * fs_483_374_385 * h9_p2 - e_4 * fs_3381_9724_1430 * h9_p4 - e_4 * fs_7560_221_3 * r_2 * h7_p2 - e_4 * fs_2835_4862_66 * r_2 * h7_p4 - e_4 * fs_114_13_30 * r_4 * h5_p2 + e_4 * fs_57_13_10 * r_4 * h5_p4 - e_5 * fs_75_4199_30030 * h11_p2 + e_5 * fs_2625_8398_286 * h11_p4 + e_5 * fs_46_187_385 * r_2 * h9_p2 + e_5 * fs_161_2431_1430 * r_2 * h9_p4 + e_5 * fs_15120_4199_3 * r_4 * h7_p2 + e_5 * fs_2835_46189_66 * r_4 * h7_p4 + e_5 * fs_152_221_30 * r_6 * h5_p2 - e_5 * fs_76_221_10 * r_6 * h5_p4 + e_6 * fs_6_4199_30030 * r_2 * h11_p2 - e_6 * fs_105_4199_286 * r_2 * h11_p4 - e_6 * fs_2_187_385 * r_4 * h9_p2 - e_6 * fs_7_2431_1430 * r_4 * h9_p4 - e_6 * fs_480_4199_3 * r_6 * h7_p2 - e_6 * fs_90_46189_66 * r_6 * h7_p4 - e_6 * fs_4_221_30 * r_8 * h5_p2 + e_6 * fs_2_221_10 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2, pc_106 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p5 = ph11_p5[k];

        pc_106[k] = - e_0 * fs_42525_32_2 * h1_p1 + e_1 * fs_16065_16_3 * h3_p1 + e_1 * fs_25515_8_2 * r_2 * h1_p1 - e_2 * fs_285_4_30 * h5_p1 - e_2 * fs_1425_4_7 * h5_p5 - e_2 * fs_8925_8_3 * r_2 * h3_p1 - e_2 * fs_18225_8_2 * r_4 * h1_p1 - e_3 * fs_21735_1144_14 * h7_p1 + e_3 * fs_9135_1144_462 * h7_p5 + e_3 * fs_570_13_30 * r_2 * h5_p1 + e_3 * fs_2850_13_7 * r_2 * h5_p5 + e_3 * fs_8925_22_3 * r_4 * h3_p1 + e_3 * fs_675_1_2 * r_6 * h1_p1 + e_4 * fs_51681_4862_10 * h9_p1 - e_4 * fs_12075_4862_143 * h9_p5 + e_4 * fs_65205_9724_14 * r_2 * h7_p1 - e_4 * fs_27405_9724_462 * r_2 * h7_p5 - e_4 * fs_114_13_30 * r_4 * h5_p1 - e_4 * fs_570_13_7 * r_4 * h5_p5 - e_4 * fs_8925_143_3 * r_6 * h3_p1 - e_4 * fs_2025_22_2 * r_8 * h1_p1 + e_5 * fs_1500_4199_33 * h11_p1 + e_5 * fs_1500_4199_286 * h11_p5 - e_5 * fs_4922_2431_10 * r_2 * h9_p1 + e_5 * fs_1150_2431_143 * r_2 * h9_p5 - e_5 * fs_65205_92378_14 * r_4 * h7_p1 + e_5 * fs_27405_92378_462 * r_4 * h7_p5 + e_5 * fs_152_221_30 * r_6 * h5_p1 + e_5 * fs_760_221_7 * r_6 * h5_p5 + e_5 * fs_595_143_3 * r_8 * h3_p1 + e_5 * fs_810_143_2 * r_10 * h1_p1 - e_6 * fs_120_4199_33 * r_2 * h11_p1 - e_6 * fs_120_4199_286 * r_2 * h11_p5 + e_6 * fs_214_2431_10 * r_4 * h9_p1 - e_6 * fs_50_2431_143 * r_4 * h9_p5 + e_6 * fs_1035_46189_14 * r_6 * h7_p1 - e_6 * fs_435_46189_462 * r_6 * h7_p5 - e_6 * fs_4_221_30 * r_8 * h5_p1 - e_6 * fs_20_221_7 * r_8 * h5_p5 - e_6 * fs_14_143_3 * r_10 * h3_p1 - e_6 * fs_18_143_2 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ph11_0, ph11_p6, ab_2, pc_107 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p6 = ph11_p6[k];

        pc_107[k] = - e_0 * fs_42525_32_3 * h1_0 - e_1 * fs_5355_16_3 * h3_0 + e_1 * fs_25515_8_3 * r_2 * h1_0 + e_2 * fs_2375_4_3 * h5_0 + e_2 * fs_2975_8_3 * r_2 * h3_0 - e_2 * fs_18225_8_3 * r_4 * h1_0 - e_3 * fs_19845_286_3 * h7_0 + e_3 * fs_945_286_2002 * h7_p6 - e_3 * fs_4750_13_3 * r_2 * h5_0 - e_3 * fs_2975_22_3 * r_4 * h3_0 + e_3 * fs_675_1_3 * r_6 * h1_0 - e_4 * fs_127995_4862_3 * h9_0 - e_4 * fs_1932_2431_1430 * h9_p6 + e_4 * fs_59535_2431_3 * r_2 * h7_0 - e_4 * fs_2835_2431_2002 * r_2 * h7_p6 + e_4 * fs_950_13_3 * r_4 * h5_0 + e_4 * fs_2975_143_3 * r_6 * h3_0 - e_4 * fs_2025_22_3 * r_8 * h1_0 - e_5 * fs_4125_4199_3 * h11_0 + e_5 * fs_375_4199_4862 * h11_p6 + e_5 * fs_12190_2431_3 * r_2 * h9_0 + e_5 * fs_368_2431_1430 * r_2 * h9_p6 - e_5 * fs_119070_46189_3 * r_4 * h7_0 + e_5 * fs_5670_46189_2002 * r_4 * h7_p6 - e_5 * fs_3800_663_3 * r_6 * h5_0 - e_5 * fs_595_429_3 * r_8 * h3_0 + e_5 * fs_810_143_3 * r_10 * h1_0 + e_6 * fs_330_4199_3 * r_2 * h11_0 - e_6 * fs_30_4199_4862 * r_2 * h11_p6 - e_6 * fs_530_2431_3 * r_4 * h9_0 - e_6 * fs_16_2431_1430 * r_4 * h9_p6 + e_6 * fs_3780_46189_3 * r_6 * h7_0 - e_6 * fs_180_46189_2002 * r_6 * h7_p6 + e_6 * fs_100_663_3 * r_8 * h5_0 + e_6 * fs_14_429_3 * r_10 * h3_0 - e_6 * fs_18_143_3 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2, pc_108 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_108[k] = e_0 * fs_14175_64_6 * h1_p1 + e_1 * f_16065_8 * h3_p1 - e_1 * fs_8505_16_6 * r_2 * h1_p1 + e_2 * fs_285_8_10 * h5_p1 - e_2 * f_8925_4 * r_2 * h3_p1 + e_2 * fs_6075_16_6 * r_4 * h1_p1 - e_3 * fs_4410_143_42 * h7_p1 - e_3 * fs_6615_572_286 * h7_p7 - e_3 * fs_285_13_10 * r_2 * h5_p1 + e_3 * f_8925_11 * r_4 * h3_p1 - e_3 * fs_225_2_6 * r_6 * h1_p1 - e_4 * fs_42021_9724_30 * h9_p1 - e_4 * fs_483_4862_2145 * h9_p7 + e_4 * fs_26460_2431_42 * r_2 * h7_p1 + e_4 * fs_19845_4862_286 * r_2 * h7_p7 + e_4 * fs_57_13_10 * r_4 * h5_p1 - e_4 * f_17850_143 * r_6 * h3_p1 + e_4 * fs_675_44_6 * r_8 * h1_p1 - e_5 * fs_750_4199_11 * h11_p1 + e_5 * fs_150_4199_24310 * h11_p7 + e_5 * fs_2001_2431_30 * r_2 * h9_p1 + e_5 * fs_46_2431_2145 * r_2 * h9_p7 - e_5 * fs_52920_46189_42 * r_4 * h7_p1 - e_5 * fs_19845_46189_286 * r_4 * h7_p7 - e_5 * fs_76_221_10 * r_6 * h5_p1 + e_5 * f_1190_143 * r_8 * h3_p1 - e_5 * fs_135_143_6 * r_10 * h1_p1 + e_6 * fs_60_4199_11 * r_2 * h11_p1 - e_6 * fs_12_4199_24310 * r_2 * h11_p7 - e_6 * fs_87_2431_30 * r_4 * h9_p1 - e_6 * fs_2_2431_2145 * r_4 * h9_p7 + e_6 * fs_1680_46189_42 * r_6 * h7_p1 + e_6 * fs_630_46189_286 * r_6 * h7_p7 + e_6 * fs_2_221_10 * r_8 * h5_p1 - e_6 * f_28_143 * r_10 * h3_p1 + e_6 * fs_3_143_6 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2, pc_109 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_109[k] = - e_1 * f_16065_16 * h3_p2 - e_2 * fs_1425_4_7 * h5_p2 + e_2 * f_8925_8 * r_2 * h3_p2 - e_3 * fs_6615_286_70 * h7_p2 + e_3 * fs_2850_13_7 * r_2 * h5_p2 - e_3 * f_8925_22 * r_4 * h3_p2 - e_4 * fs_7245_4862_66 * h9_p2 + e_4 * fs_2415_4862_7293 * h9_p8 + e_4 * fs_19845_2431_70 * r_2 * h7_p2 - e_4 * fs_570_13_7 * r_4 * h5_p2 + e_4 * f_8925_143 * r_6 * h3_p2 - e_5 * fs_75_4199_143 * h11_p2 + e_5 * fs_75_4199_46189 * h11_p8 + e_5 * fs_690_2431_66 * r_2 * h9_p2 - e_5 * fs_230_2431_7293 * r_2 * h9_p8 - e_5 * fs_39690_46189_70 * r_4 * h7_p2 + e_5 * fs_760_221_7 * r_6 * h5_p2 - e_5 * f_595_143 * r_8 * h3_p2 + e_6 * fs_6_4199_143 * r_2 * h11_p2 - e_6 * fs_6_4199_46189 * r_2 * h11_p8 - e_6 * fs_30_2431_66 * r_4 * h9_p2 + e_6 * fs_10_2431_7293 * r_4 * h9_p8 + e_6 * fs_1260_46189_70 * r_6 * h7_p2 - e_6 * fs_20_221_7 * r_8 * h5_p2 + e_6 * f_14_143 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m8, ph9_m1, ph11_m9, ph11_m8, ph11_m1, ab_2, pc_110, pc_111 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m1 = ph11_m1[k];

        pc_110[k] = e_0 * fs_14175_64_2 * h1_m1 + e_1 * fs_16065_16_3 * h3_m1 - e_1 * fs_8505_16_2 * r_2 * h1_m1 + e_2 * fs_1425_8_30 * h5_m1 - e_2 * fs_8925_8_3 * r_2 * h3_m1 + e_2 * fs_6075_16_2 * r_4 * h1_m1 + e_3 * fs_11025_286_14 * h7_m1 - e_3 * fs_1425_13_30 * r_2 * h5_m1 + e_3 * fs_8925_22_3 * r_4 * h3_m1 - e_3 * fs_225_2_2 * r_6 * h1_m1 + e_4 * fs_1449_4862_12155 * h9_m9 + e_4 * fs_21735_9724_10 * h9_m1 - e_4 * fs_33075_2431_14 * r_2 * h7_m1 + e_4 * fs_285_13_30 * r_4 * h5_m1 - e_4 * fs_8925_143_3 * r_6 * h3_m1 + e_4 * fs_675_44_2 * r_8 * h1_m1 + e_5 * fs_75_4199_92378 * h11_m9 + e_5 * fs_75_4199_33 * h11_m1 - e_5 * fs_138_2431_12155 * r_2 * h9_m9 - e_5 * fs_1035_2431_10 * r_2 * h9_m1 + e_5 * fs_66150_46189_14 * r_4 * h7_m1 - e_5 * fs_380_221_30 * r_6 * h5_m1 + e_5 * fs_595_143_3 * r_8 * h3_m1 - e_5 * fs_135_143_2 * r_10 * h1_m1 - e_6 * fs_6_4199_92378 * r_2 * h11_m9 - e_6 * fs_6_4199_33 * r_2 * h11_m1 + e_6 * fs_6_2431_12155 * r_4 * h9_m9 + e_6 * fs_45_2431_10 * r_4 * h9_m1 - e_6 * fs_2100_46189_14 * r_6 * h7_m1 + e_6 * fs_10_221_30 * r_8 * h5_m1 - e_6 * fs_14_143_3 * r_10 * h3_m1 + e_6 * fs_3_143_2 * r_12 * h1_m1;

        pc_111[k] = - e_4 * fs_1449_2431_2431 * h9_m8 + e_5 * fs_75_4199_138567 * h11_m8 + e_5 * fs_276_2431_2431 * r_2 * h9_m8 - e_6 * fs_6_4199_138567 * r_2 * h11_m8 - e_6 * fs_12_2431_2431 * r_4 * h9_m8;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2, pc_112 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_112[k] = e_0 * fs_42525_64_10 * h1_m1 - e_1 * fs_5355_16_15 * h3_m1 - e_1 * fs_25515_16_10 * r_2 * h1_m1 - e_2 * fs_1235_8_6 * h5_m1 + e_2 * fs_2975_8_15 * r_2 * h3_m1 + e_2 * fs_18225_16_10 * r_4 * h1_m1 + e_3 * fs_2205_572_4290 * h7_m7 + e_3 * fs_945_44_70 * h7_m1 + e_3 * fs_95_1_6 * r_2 * h5_m1 - e_3 * fs_2975_22_15 * r_4 * h3_m1 - e_3 * fs_675_2_10 * r_6 * h1_m1 - e_4 * fs_14007_4862_143 * h9_m7 + e_4 * fs_177261_9724_2 * h9_m1 - e_4 * fs_6615_4862_4290 * r_2 * h7_m7 - e_4 * fs_2835_374_70 * r_2 * h7_m1 - e_4 * fs_19_1_6 * r_4 * h5_m1 + e_4 * fs_2975_143_15 * r_6 * h3_m1 + e_4 * fs_2025_44_10 * r_8 * h1_m1 + e_5 * fs_225_4199_14586 * h11_m7 + e_5 * fs_225_4199_165 * h11_m1 + e_5 * fs_1334_2431_143 * r_2 * h9_m7 - e_5 * fs_8441_2431_2 * r_2 * h9_m1 + e_5 * fs_6615_46189_4290 * r_4 * h7_m7 + e_5 * fs_2835_3553_70 * r_4 * h7_m1 + e_5 * fs_76_51_6 * r_6 * h5_m1 - e_5 * fs_595_429_15 * r_8 * h3_m1 - e_5 * fs_405_143_10 * r_10 * h1_m1 - e_6 * fs_18_4199_14586 * r_2 * h11_m7 - e_6 * fs_18_4199_165 * r_2 * h11_m1 - e_6 * fs_58_2431_143 * r_4 * h9_m7 + e_6 * fs_367_2431_2 * r_4 * h9_m1 - e_6 * fs_210_46189_4290 * r_6 * h7_m7 - e_6 * fs_90_3553_70 * r_6 * h7_m1 - e_6 * fs_2_51_6 * r_8 * h5_m1 + e_6 * fs_14_429_15 * r_10 * h3_m1 + e_6 * fs_9_143_10 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2, pc_113 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_113[k] = - e_1 * f_16065_16 * h3_m2 + e_2 * fs_285_1_7 * h5_m2 + e_2 * f_8925_8 * r_2 * h3_m2 + e_3 * fs_315_1144_10010 * h7_m6 - e_3 * fs_9135_1144_70 * h7_m2 - e_3 * fs_2280_13_7 * r_2 * h5_m2 - e_3 * f_8925_22 * r_4 * h3_m2 - e_4 * fs_4347_4862_286 * h9_m6 - e_4 * fs_19803_4862_66 * h9_m2 - e_4 * fs_945_9724_10010 * r_2 * h7_m6 + e_4 * fs_27405_9724_70 * r_2 * h7_m2 + e_4 * fs_456_13_7 * r_4 * h5_m2 + e_4 * f_8925_143 * r_6 * h3_m2 + e_5 * fs_150_4199_24310 * h11_m6 - e_5 * fs_450_4199_143 * h11_m2 + e_5 * fs_414_2431_286 * r_2 * h9_m6 + e_5 * fs_1886_2431_66 * r_2 * h9_m2 + e_5 * fs_945_92378_10010 * r_4 * h7_m6 - e_5 * fs_27405_92378_70 * r_4 * h7_m2 - e_5 * fs_608_221_7 * r_6 * h5_m2 - e_5 * f_595_143 * r_8 * h3_m2 - e_6 * fs_12_4199_24310 * r_2 * h11_m6 + e_6 * fs_36_4199_143 * r_2 * h11_m2 - e_6 * fs_18_2431_286 * r_4 * h9_m6 - e_6 * fs_82_2431_66 * r_4 * h9_m2 - e_6 * fs_15_46189_10010 * r_6 * h7_m6 + e_6 * fs_435_46189_70 * r_6 * h7_m2 + e_6 * fs_16_221_7 * r_8 * h5_m2 + e_6 * f_14_143 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2, pc_114 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_114[k] = e_1 * fs_5355_16_42 * h3_m3 + e_2 * fs_1425_8_30 * h5_m5 - e_2 * fs_1235_8_6 * h5_m3 - e_2 * fs_2975_8_42 * r_2 * h3_m3 - e_3 * fs_11655_572_55 * h7_m5 - e_3 * fs_1575_44_5 * h7_m3 - e_3 * fs_1425_13_30 * r_2 * h5_m5 + e_3 * fs_95_1_6 * r_2 * h5_m3 + e_3 * fs_2975_22_42 * r_4 * h3_m3 + e_4 * fs_483_9724_30030 * h9_m5 + e_4 * fs_71001_9724_22 * h9_m3 + e_4 * fs_34965_4862_55 * r_2 * h7_m5 + e_4 * fs_4725_374_5 * r_2 * h7_m3 + e_4 * fs_285_13_30 * r_4 * h5_m5 - e_4 * fs_19_1_6 * r_4 * h5_m3 - e_4 * fs_2975_143_42 * r_6 * h3_m3 + e_5 * fs_150_4199_15015 * h11_m5 + e_5 * fs_525_4199_286 * h11_m3 - e_5 * fs_23_2431_30030 * r_2 * h9_m5 - e_5 * fs_3381_2431_22 * r_2 * h9_m3 - e_5 * fs_34965_46189_55 * r_4 * h7_m5 - e_5 * fs_4725_3553_5 * r_4 * h7_m3 - e_5 * fs_380_221_30 * r_6 * h5_m5 + e_5 * fs_76_51_6 * r_6 * h5_m3 + e_5 * fs_595_429_42 * r_8 * h3_m3 - e_6 * fs_12_4199_15015 * r_2 * h11_m5 - e_6 * fs_42_4199_286 * r_2 * h11_m3 + e_6 * fs_1_2431_30030 * r_4 * h9_m5 + e_6 * fs_147_2431_22 * r_4 * h9_m3 + e_6 * fs_1110_46189_55 * r_6 * h7_m5 + e_6 * fs_150_3553_5 * r_6 * h7_m3 + e_6 * fs_10_221_30 * r_8 * h5_m5 - e_6 * fs_2_51_6 * r_8 * h5_m3 - e_6 * fs_14_429_42 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p4, ph7_p4, ph9_p4, ph11_p4, ab_2, pc_115 : simd::cache_line_size())
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

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p4 = ph11_p4[k];

        pc_115[k] = e_2 * fs_285_1_5 * h5_p4 - e_3 * fs_6300_143_33 * h7_p4 - e_3 * fs_2280_13_5 * r_2 * h5_p4 + e_4 * fs_3381_2431_715 * h9_p4 + e_4 * fs_37800_2431_33 * r_2 * h7_p4 + e_4 * fs_456_13_5 * r_4 * h5_p4 + e_5 * fs_1575_4199_143 * h11_p4 - e_5 * fs_644_2431_715 * r_2 * h9_p4 - e_5 * fs_75600_46189_33 * r_4 * h7_p4 - e_5 * fs_608_221_5 * r_6 * h5_p4 - e_6 * fs_126_4199_143 * r_2 * h11_p4 + e_6 * fs_28_2431_715 * r_4 * h9_p4 + e_6 * fs_2400_46189_33 * r_6 * h7_p4 + e_6 * fs_16_221_5 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2, pc_116 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_116[k] = - e_1 * fs_5355_16_42 * h3_p3 + e_2 * fs_1235_8_6 * h5_p3 + e_2 * fs_1425_8_30 * h5_p5 + e_2 * fs_2975_8_42 * r_2 * h3_p3 + e_3 * fs_1575_44_5 * h7_p3 - e_3 * fs_11655_572_55 * h7_p5 - e_3 * fs_95_1_6 * r_2 * h5_p3 - e_3 * fs_1425_13_30 * r_2 * h5_p5 - e_3 * fs_2975_22_42 * r_4 * h3_p3 - e_4 * fs_71001_9724_22 * h9_p3 + e_4 * fs_483_9724_30030 * h9_p5 - e_4 * fs_4725_374_5 * r_2 * h7_p3 + e_4 * fs_34965_4862_55 * r_2 * h7_p5 + e_4 * fs_19_1_6 * r_4 * h5_p3 + e_4 * fs_285_13_30 * r_4 * h5_p5 + e_4 * fs_2975_143_42 * r_6 * h3_p3 - e_5 * fs_525_4199_286 * h11_p3 + e_5 * fs_150_4199_15015 * h11_p5 + e_5 * fs_3381_2431_22 * r_2 * h9_p3 - e_5 * fs_23_2431_30030 * r_2 * h9_p5 + e_5 * fs_4725_3553_5 * r_4 * h7_p3 - e_5 * fs_34965_46189_55 * r_4 * h7_p5 - e_5 * fs_76_51_6 * r_6 * h5_p3 - e_5 * fs_380_221_30 * r_6 * h5_p5 - e_5 * fs_595_429_42 * r_8 * h3_p3 + e_6 * fs_42_4199_286 * r_2 * h11_p3 - e_6 * fs_12_4199_15015 * r_2 * h11_p5 - e_6 * fs_147_2431_22 * r_4 * h9_p3 + e_6 * fs_1_2431_30030 * r_4 * h9_p5 - e_6 * fs_150_3553_5 * r_6 * h7_p3 + e_6 * fs_1110_46189_55 * r_6 * h7_p5 + e_6 * fs_2_51_6 * r_8 * h5_p3 + e_6 * fs_10_221_30 * r_8 * h5_p5 + e_6 * fs_14_429_42 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2, pc_117 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_117[k] = e_1 * f_16065_16 * h3_p2 - e_2 * fs_285_1_7 * h5_p2 - e_2 * f_8925_8 * r_2 * h3_p2 + e_3 * fs_9135_1144_70 * h7_p2 + e_3 * fs_315_1144_10010 * h7_p6 + e_3 * fs_2280_13_7 * r_2 * h5_p2 + e_3 * f_8925_22 * r_4 * h3_p2 + e_4 * fs_19803_4862_66 * h9_p2 - e_4 * fs_4347_4862_286 * h9_p6 - e_4 * fs_27405_9724_70 * r_2 * h7_p2 - e_4 * fs_945_9724_10010 * r_2 * h7_p6 - e_4 * fs_456_13_7 * r_4 * h5_p2 - e_4 * f_8925_143 * r_6 * h3_p2 + e_5 * fs_450_4199_143 * h11_p2 + e_5 * fs_150_4199_24310 * h11_p6 - e_5 * fs_1886_2431_66 * r_2 * h9_p2 + e_5 * fs_414_2431_286 * r_2 * h9_p6 + e_5 * fs_27405_92378_70 * r_4 * h7_p2 + e_5 * fs_945_92378_10010 * r_4 * h7_p6 + e_5 * fs_608_221_7 * r_6 * h5_p2 + e_5 * f_595_143 * r_8 * h3_p2 - e_6 * fs_36_4199_143 * r_2 * h11_p2 - e_6 * fs_12_4199_24310 * r_2 * h11_p6 + e_6 * fs_82_2431_66 * r_4 * h9_p2 - e_6 * fs_18_2431_286 * r_4 * h9_p6 - e_6 * fs_435_46189_70 * r_6 * h7_p2 - e_6 * fs_15_46189_10010 * r_6 * h7_p6 - e_6 * fs_16_221_7 * r_8 * h5_p2 - e_6 * f_14_143 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2, pc_118 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_118[k] = - e_0 * fs_42525_64_10 * h1_p1 + e_1 * fs_5355_16_15 * h3_p1 + e_1 * fs_25515_16_10 * r_2 * h1_p1 + e_2 * fs_1235_8_6 * h5_p1 - e_2 * fs_2975_8_15 * r_2 * h3_p1 - e_2 * fs_18225_16_10 * r_4 * h1_p1 - e_3 * fs_945_44_70 * h7_p1 + e_3 * fs_2205_572_4290 * h7_p7 - e_3 * fs_95_1_6 * r_2 * h5_p1 + e_3 * fs_2975_22_15 * r_4 * h3_p1 + e_3 * fs_675_2_10 * r_6 * h1_p1 - e_4 * fs_177261_9724_2 * h9_p1 - e_4 * fs_14007_4862_143 * h9_p7 + e_4 * fs_2835_374_70 * r_2 * h7_p1 - e_4 * fs_6615_4862_4290 * r_2 * h7_p7 + e_4 * fs_19_1_6 * r_4 * h5_p1 - e_4 * fs_2975_143_15 * r_6 * h3_p1 - e_4 * fs_2025_44_10 * r_8 * h1_p1 - e_5 * fs_225_4199_165 * h11_p1 + e_5 * fs_225_4199_14586 * h11_p7 + e_5 * fs_8441_2431_2 * r_2 * h9_p1 + e_5 * fs_1334_2431_143 * r_2 * h9_p7 - e_5 * fs_2835_3553_70 * r_4 * h7_p1 + e_5 * fs_6615_46189_4290 * r_4 * h7_p7 - e_5 * fs_76_51_6 * r_6 * h5_p1 + e_5 * fs_595_429_15 * r_8 * h3_p1 + e_5 * fs_405_143_10 * r_10 * h1_p1 + e_6 * fs_18_4199_165 * r_2 * h11_p1 - e_6 * fs_18_4199_14586 * r_2 * h11_p7 - e_6 * fs_367_2431_2 * r_4 * h9_p1 - e_6 * fs_58_2431_143 * r_4 * h9_p7 + e_6 * fs_90_3553_70 * r_6 * h7_p1 - e_6 * fs_210_46189_4290 * r_6 * h7_p7 + e_6 * fs_2_51_6 * r_8 * h5_p1 - e_6 * fs_14_429_15 * r_10 * h3_p1 - e_6 * fs_9_143_10 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph9_p8, ph11_0, ph11_p8, ab_2, pc_119 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p8 = ph11_p8[k];

        pc_119[k] = - e_0 * fs_14175_16_5 * h1_0 - e_1 * fs_16065_16_5 * h3_0 + e_1 * fs_8505_4_5 * r_2 * h1_0 + e_2 * fs_285_1_5 * h5_0 + e_2 * fs_8925_8_5 * r_2 * h3_0 - e_2 * fs_6075_4_5 * r_4 * h1_0 + e_3 * fs_37485_286_5 * h7_0 - e_3 * fs_2280_13_5 * r_2 * h5_0 - e_3 * fs_8925_22_5 * r_4 * h3_0 + e_3 * fs_450_1_5 * r_6 * h1_0 + e_4 * fs_1449_143_5 * h9_0 - e_4 * fs_1449_2431_2431 * h9_p8 - e_4 * fs_6615_143_5 * r_2 * h7_0 + e_4 * fs_456_13_5 * r_4 * h5_0 + e_4 * fs_8925_143_5 * r_6 * h3_0 - e_4 * fs_675_11_5 * r_8 * h1_0 + e_5 * fs_825_4199_5 * h11_0 + e_5 * fs_75_4199_138567 * h11_p8 - e_5 * fs_276_143_5 * r_2 * h9_0 + e_5 * fs_276_2431_2431 * r_2 * h9_p8 + e_5 * fs_13230_2717_5 * r_4 * h7_0 - e_5 * fs_608_221_5 * r_6 * h5_0 - e_5 * fs_595_143_5 * r_8 * h3_0 + e_5 * fs_540_143_5 * r_10 * h1_0 - e_6 * fs_66_4199_5 * r_2 * h11_0 - e_6 * fs_6_4199_138567 * r_2 * h11_p8 + e_6 * fs_12_143_5 * r_4 * h9_0 - e_6 * fs_12_2431_2431 * r_4 * h9_p8 - e_6 * fs_420_2717_5 * r_6 * h7_0 + e_6 * fs_16_221_5 * r_8 * h5_0 + e_6 * fs_14_143_5 * r_10 * h3_0 - e_6 * fs_12_143_5 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_m10, ph11_p1, ph11_p9, ab_2, pc_120, pc_121 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_120[k] = e_0 * fs_14175_64_2 * h1_p1 + e_1 * fs_16065_16_3 * h3_p1 - e_1 * fs_8505_16_2 * r_2 * h1_p1 + e_2 * fs_1425_8_30 * h5_p1 - e_2 * fs_8925_8_3 * r_2 * h3_p1 + e_2 * fs_6075_16_2 * r_4 * h1_p1 + e_3 * fs_11025_286_14 * h7_p1 - e_3 * fs_1425_13_30 * r_2 * h5_p1 + e_3 * fs_8925_22_3 * r_4 * h3_p1 - e_3 * fs_225_2_2 * r_6 * h1_p1 + e_4 * fs_21735_9724_10 * h9_p1 + e_4 * fs_1449_4862_12155 * h9_p9 - e_4 * fs_33075_2431_14 * r_2 * h7_p1 + e_4 * fs_285_13_30 * r_4 * h5_p1 - e_4 * fs_8925_143_3 * r_6 * h3_p1 + e_4 * fs_675_44_2 * r_8 * h1_p1 + e_5 * fs_75_4199_33 * h11_p1 + e_5 * fs_75_4199_92378 * h11_p9 - e_5 * fs_1035_2431_10 * r_2 * h9_p1 - e_5 * fs_138_2431_12155 * r_2 * h9_p9 + e_5 * fs_66150_46189_14 * r_4 * h7_p1 - e_5 * fs_380_221_30 * r_6 * h5_p1 + e_5 * fs_595_143_3 * r_8 * h3_p1 - e_5 * fs_135_143_2 * r_10 * h1_p1 - e_6 * fs_6_4199_33 * r_2 * h11_p1 - e_6 * fs_6_4199_92378 * r_2 * h11_p9 + e_6 * fs_45_2431_10 * r_4 * h9_p1 + e_6 * fs_6_2431_12155 * r_4 * h9_p9 - e_6 * fs_2100_46189_14 * r_6 * h7_p1 + e_6 * fs_10_221_30 * r_8 * h5_p1 - e_6 * fs_14_143_3 * r_10 * h3_p1 + e_6 * fs_3_143_2 * r_12 * h1_p1;

        pc_121[k] = e_5 * fs_75_4199_176358 * h11_m10 - e_6 * fs_6_4199_176358 * r_2 * h11_m10;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m1, ph11_m9, ph11_m1, ab_2, pc_122 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m1 = ph11_m1[k];

        pc_122[k] = e_0 * fs_14175_64_110 * h1_m1 - e_1 * fs_8505_16_110 * r_2 * h1_m1 - e_2 * fs_855_8_66 * h5_m1 + e_2 * fs_6075_16_110 * r_4 * h1_m1 - e_3 * fs_315_52_770 * h7_m1 + e_3 * fs_855_13_66 * r_2 * h5_m1 - e_3 * fs_225_2_110 * r_6 * h1_m1 - e_4 * fs_1449_442_221 * h9_m9 - e_4 * fs_1449_748_22 * h9_m1 + e_4 * fs_945_442_770 * r_2 * h7_m1 - e_4 * fs_171_13_66 * r_4 * h5_m1 + e_4 * fs_675_44_110 * r_8 * h1_m1 + e_5 * fs_150_4199_41990 * h11_m9 - e_5 * fs_150_4199_15 * h11_m1 + e_5 * fs_138_221_221 * r_2 * h9_m9 + e_5 * fs_69_187_22 * r_2 * h9_m1 - e_5 * fs_945_4199_770 * r_4 * h7_m1 + e_5 * fs_228_221_66 * r_6 * h5_m1 - e_5 * fs_135_143_110 * r_10 * h1_m1 - e_6 * fs_12_4199_41990 * r_2 * h11_m9 + e_6 * fs_12_4199_15 * r_2 * h11_m1 - e_6 * fs_6_221_221 * r_4 * h9_m9 - e_6 * fs_3_187_22 * r_4 * h9_m1 + e_6 * fs_30_4199_770 * r_6 * h7_m1 - e_6 * fs_6_221_66 * r_8 * h5_m1 + e_6 * fs_3_143_110 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2, pc_123 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_123[k] = - e_1 * fs_5355_16_33 * h3_m2 + e_2 * fs_95_4_231 * h5_m2 + e_2 * fs_2975_8_33 * r_2 * h3_m2 + e_3 * fs_630_143_2310 * h7_m2 - e_3 * fs_190_13_231 * r_2 * h5_m2 - e_3 * fs_2975_22_33 * r_4 * h3_m2 - e_4 * fs_483_442_221 * h9_m8 + e_4 * fs_5313_442_2 * h9_m2 - e_4 * fs_3780_2431_2310 * r_2 * h7_m2 + e_4 * fs_38_13_231 * r_4 * h5_m2 + e_4 * fs_2975_143_33 * r_6 * h3_m2 + e_5 * fs_225_4199_12597 * h11_m8 + e_5 * fs_225_4199_39 * h11_m2 + e_5 * fs_46_221_221 * r_2 * h9_m8 - e_5 * fs_506_221_2 * r_2 * h9_m2 + e_5 * fs_7560_46189_2310 * r_4 * h7_m2 - e_5 * fs_152_663_231 * r_6 * h5_m2 - e_5 * fs_595_429_33 * r_8 * h3_m2 - e_6 * fs_18_4199_12597 * r_2 * h11_m8 - e_6 * fs_18_4199_39 * r_2 * h11_m2 - e_6 * fs_2_221_221 * r_4 * h9_m8 + e_6 * fs_22_221_2 * r_4 * h9_m2 - e_6 * fs_240_46189_2310 * r_6 * h7_m2 + e_6 * fs_4_663_231 * r_8 * h5_m2 + e_6 * fs_14_429_33 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2, pc_124 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_124[k] = e_1 * fs_5355_16_33 * h3_m3 + e_2 * fs_95_4_231 * h5_m3 - e_2 * fs_2975_8_33 * r_2 * h3_m3 - e_3 * fs_2205_104_130 * h7_m7 - e_3 * fs_7875_1144_770 * h7_m3 - e_3 * fs_190_13_231 * r_2 * h5_m3 + e_3 * fs_2975_22_33 * r_4 * h3_m3 + e_4 * fs_483_221_39 * h9_m7 - e_4 * fs_4347_442_7 * h9_m3 + e_4 * fs_6615_884_130 * r_2 * h7_m7 + e_4 * fs_23625_9724_770 * r_2 * h7_m3 + e_4 * fs_38_13_231 * r_4 * h5_m3 - e_4 * fs_2975_143_33 * r_6 * h3_m3 + e_5 * fs_900_4199_442 * h11_m7 - e_5 * fs_300_4199_91 * h11_m3 - e_5 * fs_92_221_39 * r_2 * h9_m7 + e_5 * fs_414_221_7 * r_2 * h9_m3 - e_5 * fs_6615_8398_130 * r_4 * h7_m7 - e_5 * fs_23625_92378_770 * r_4 * h7_m3 - e_5 * fs_152_663_231 * r_6 * h5_m3 + e_5 * fs_595_429_33 * r_8 * h3_m3 - e_6 * fs_72_4199_442 * r_2 * h11_m7 + e_6 * fs_24_4199_91 * r_2 * h11_m3 + e_6 * fs_4_221_39 * r_4 * h9_m7 - e_6 * fs_18_221_7 * r_4 * h9_m3 + e_6 * fs_105_4199_130 * r_6 * h7_m7 + e_6 * fs_375_46189_770 * r_6 * h7_m3 + e_6 * fs_4_663_231 * r_8 * h5_m3 - e_6 * fs_14_429_33 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph5_p5, ph7_m6, ph7_m4, ph7_p5, ph9_m6, ph9_m4, ph9_p5, ph11_m6, ph11_m4, ph11_p5, ab_2, pc_125, pc_126 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_p5 = ph11_p5[k];

        pc_125[k] = - e_2 * fs_855_8_66 * h5_m4 - e_3 * fs_315_13_65 * h7_m6 + e_3 * fs_1575_52_10 * h7_m4 + e_3 * fs_855_13_66 * r_2 * h5_m4 + e_4 * fs_1449_442_91 * h9_m6 + e_4 * fs_3381_884_78 * h9_m4 + e_4 * fs_1890_221_65 * r_2 * h7_m6 - e_4 * fs_4725_442_10 * r_2 * h7_m4 - e_4 * fs_171_13_66 * r_4 * h5_m4 + e_5 * fs_150_4199_7735 * h11_m6 + e_5 * fs_525_8398_390 * h11_m4 - e_5 * fs_138_221_91 * r_2 * h9_m6 - e_5 * fs_161_221_78 * r_2 * h9_m4 - e_5 * fs_3780_4199_65 * r_4 * h7_m6 + e_5 * fs_4725_4199_10 * r_4 * h7_m4 + e_5 * fs_228_221_66 * r_6 * h5_m4 - e_6 * fs_12_4199_7735 * r_2 * h11_m6 - e_6 * fs_21_4199_390 * r_2 * h11_m4 + e_6 * fs_6_221_91 * r_4 * h9_m6 + e_6 * fs_7_221_78 * r_4 * h9_m4 + e_6 * fs_120_4199_65 * r_6 * h7_m6 - e_6 * fs_150_4199_10 * r_6 * h7_m4 - e_6 * fs_6_221_66 * r_8 * h5_m4;

        pc_126[k] = - e_2 * fs_1425_4_11 * h5_p5 - e_3 * fs_1575_52_6 * h7_p5 + e_3 * fs_2850_13_11 * r_2 * h5_p5 + e_4 * fs_2415_442_91 * h9_p5 + e_4 * fs_4725_442_6 * r_2 * h7_p5 - e_4 * fs_570_13_11 * r_4 * h5_p5 + e_5 * fs_900_4199_182 * h11_p5 - e_5 * fs_230_221_91 * r_2 * h9_p5 - e_5 * fs_4725_4199_6 * r_4 * h7_p5 + e_5 * fs_760_221_11 * r_6 * h5_p5 - e_6 * fs_72_4199_182 * r_2 * h11_p5 + e_6 * fs_10_221_91 * r_4 * h9_p5 + e_6 * fs_150_4199_6 * r_6 * h7_p5 - e_6 * fs_20_221_11 * r_8 * h5_p5;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2, pc_127 : simd::cache_line_size())
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

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_127[k] = e_2 * fs_855_8_66 * h5_p4 - e_3 * fs_1575_52_10 * h7_p4 - e_3 * fs_315_13_65 * h7_p6 - e_3 * fs_855_13_66 * r_2 * h5_p4 - e_4 * fs_3381_884_78 * h9_p4 + e_4 * fs_1449_442_91 * h9_p6 + e_4 * fs_4725_442_10 * r_2 * h7_p4 + e_4 * fs_1890_221_65 * r_2 * h7_p6 + e_4 * fs_171_13_66 * r_4 * h5_p4 - e_5 * fs_525_8398_390 * h11_p4 + e_5 * fs_150_4199_7735 * h11_p6 + e_5 * fs_161_221_78 * r_2 * h9_p4 - e_5 * fs_138_221_91 * r_2 * h9_p6 - e_5 * fs_4725_4199_10 * r_4 * h7_p4 - e_5 * fs_3780_4199_65 * r_4 * h7_p6 - e_5 * fs_228_221_66 * r_6 * h5_p4 + e_6 * fs_21_4199_390 * r_2 * h11_p4 - e_6 * fs_12_4199_7735 * r_2 * h11_p6 - e_6 * fs_7_221_78 * r_4 * h9_p4 + e_6 * fs_6_221_91 * r_4 * h9_p6 + e_6 * fs_150_4199_10 * r_6 * h7_p4 + e_6 * fs_120_4199_65 * r_6 * h7_p6 + e_6 * fs_6_221_66 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2, pc_128 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_128[k] = - e_1 * fs_5355_16_33 * h3_p3 - e_2 * fs_95_4_231 * h5_p3 + e_2 * fs_2975_8_33 * r_2 * h3_p3 + e_3 * fs_7875_1144_770 * h7_p3 - e_3 * fs_2205_104_130 * h7_p7 + e_3 * fs_190_13_231 * r_2 * h5_p3 - e_3 * fs_2975_22_33 * r_4 * h3_p3 + e_4 * fs_4347_442_7 * h9_p3 + e_4 * fs_483_221_39 * h9_p7 - e_4 * fs_23625_9724_770 * r_2 * h7_p3 + e_4 * fs_6615_884_130 * r_2 * h7_p7 - e_4 * fs_38_13_231 * r_4 * h5_p3 + e_4 * fs_2975_143_33 * r_6 * h3_p3 + e_5 * fs_300_4199_91 * h11_p3 + e_5 * fs_900_4199_442 * h11_p7 - e_5 * fs_414_221_7 * r_2 * h9_p3 - e_5 * fs_92_221_39 * r_2 * h9_p7 + e_5 * fs_23625_92378_770 * r_4 * h7_p3 - e_5 * fs_6615_8398_130 * r_4 * h7_p7 + e_5 * fs_152_663_231 * r_6 * h5_p3 - e_5 * fs_595_429_33 * r_8 * h3_p3 - e_6 * fs_24_4199_91 * r_2 * h11_p3 - e_6 * fs_72_4199_442 * r_2 * h11_p7 + e_6 * fs_18_221_7 * r_4 * h9_p3 + e_6 * fs_4_221_39 * r_4 * h9_p7 - e_6 * fs_375_46189_770 * r_6 * h7_p3 + e_6 * fs_105_4199_130 * r_6 * h7_p7 - e_6 * fs_4_663_231 * r_8 * h5_p3 + e_6 * fs_14_429_33 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2, pc_129 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_129[k] = e_1 * fs_5355_16_33 * h3_p2 - e_2 * fs_95_4_231 * h5_p2 - e_2 * fs_2975_8_33 * r_2 * h3_p2 - e_3 * fs_630_143_2310 * h7_p2 + e_3 * fs_190_13_231 * r_2 * h5_p2 + e_3 * fs_2975_22_33 * r_4 * h3_p2 - e_4 * fs_5313_442_2 * h9_p2 - e_4 * fs_483_442_221 * h9_p8 + e_4 * fs_3780_2431_2310 * r_2 * h7_p2 - e_4 * fs_38_13_231 * r_4 * h5_p2 - e_4 * fs_2975_143_33 * r_6 * h3_p2 - e_5 * fs_225_4199_39 * h11_p2 + e_5 * fs_225_4199_12597 * h11_p8 + e_5 * fs_506_221_2 * r_2 * h9_p2 + e_5 * fs_46_221_221 * r_2 * h9_p8 - e_5 * fs_7560_46189_2310 * r_4 * h7_p2 + e_5 * fs_152_663_231 * r_6 * h5_p2 + e_5 * fs_595_429_33 * r_8 * h3_p2 + e_6 * fs_18_4199_39 * r_2 * h11_p2 - e_6 * fs_18_4199_12597 * r_2 * h11_p8 - e_6 * fs_22_221_2 * r_4 * h9_p2 - e_6 * fs_2_221_221 * r_4 * h9_p8 + e_6 * fs_240_46189_2310 * r_6 * h7_p2 - e_6 * fs_4_663_231 * r_8 * h5_p2 - e_6 * fs_14_429_33 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2, pc_130 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_130[k] = - e_0 * fs_14175_64_110 * h1_p1 + e_1 * fs_8505_16_110 * r_2 * h1_p1 + e_2 * fs_855_8_66 * h5_p1 - e_2 * fs_6075_16_110 * r_4 * h1_p1 + e_3 * fs_315_52_770 * h7_p1 - e_3 * fs_855_13_66 * r_2 * h5_p1 + e_3 * fs_225_2_110 * r_6 * h1_p1 + e_4 * fs_1449_748_22 * h9_p1 - e_4 * fs_1449_442_221 * h9_p9 - e_4 * fs_945_442_770 * r_2 * h7_p1 + e_4 * fs_171_13_66 * r_4 * h5_p1 - e_4 * fs_675_44_110 * r_8 * h1_p1 + e_5 * fs_150_4199_15 * h11_p1 + e_5 * fs_150_4199_41990 * h11_p9 - e_5 * fs_69_187_22 * r_2 * h9_p1 + e_5 * fs_138_221_221 * r_2 * h9_p9 + e_5 * fs_945_4199_770 * r_4 * h7_p1 - e_5 * fs_228_221_66 * r_6 * h5_p1 + e_5 * fs_135_143_110 * r_10 * h1_p1 - e_6 * fs_12_4199_15 * r_2 * h11_p1 - e_6 * fs_12_4199_41990 * r_2 * h11_p9 + e_6 * fs_3_187_22 * r_4 * h9_p1 - e_6 * fs_6_221_221 * r_4 * h9_p9 - e_6 * fs_30_4199_770 * r_6 * h7_p1 + e_6 * fs_6_221_66 * r_8 * h5_p1 - e_6 * fs_3_143_110 * r_12 * h1_p1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph11_0, ph11_p10, ab_2, pc_131 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p10 = ph11_p10[k];

        pc_131[k] = - e_0 * fs_14175_32_11 * h1_0 - e_1 * fs_16065_16_11 * h3_0 + e_1 * fs_8505_8_11 * r_2 * h1_0 - e_2 * fs_1425_4_11 * h5_0 + e_2 * fs_8925_8_11 * r_2 * h3_0 - e_2 * fs_6075_8_11 * r_4 * h1_0 - e_3 * fs_11025_286_11 * h7_0 + e_3 * fs_2850_13_11 * r_2 * h5_0 - e_3 * fs_8925_22_11 * r_4 * h3_0 + e_3 * fs_225_1_11 * r_6 * h1_0 - e_4 * fs_7245_4862_11 * h9_0 + e_4 * fs_33075_2431_11 * r_2 * h7_0 - e_4 * fs_570_13_11 * r_4 * h5_0 + e_4 * fs_8925_143_11 * r_6 * h3_0 - e_4 * fs_675_22_11 * r_8 * h1_0 - e_5 * fs_75_4199_11 * h11_0 + e_5 * fs_75_4199_176358 * h11_p10 + e_5 * fs_690_2431_11 * r_2 * h9_0 - e_5 * fs_66150_46189_11 * r_4 * h7_0 + e_5 * fs_760_221_11 * r_6 * h5_0 - e_5 * fs_595_143_11 * r_8 * h3_0 + e_5 * fs_270_143_11 * r_10 * h1_0 + e_6 * fs_6_4199_11 * r_2 * h11_0 - e_6 * fs_6_4199_176358 * r_2 * h11_p10 - e_6 * fs_30_2431_11 * r_4 * h9_0 + e_6 * fs_2100_46189_11 * r_6 * h7_0 - e_6 * fs_20_221_11 * r_8 * h5_0 + e_6 * fs_14_143_11 * r_10 * h3_0 - e_6 * fs_6_143_11 * r_12 * h1_0;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m1, ph11_m11, ph11_m1, ab_2, pc_132 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m11 = ph11_m11[k];
        const auto h11_m1 = ph11_m1[k];

        pc_132[k] = e_0 * fs_14175_32_33 * h1_m1 + e_1 * fs_16065_32_22 * h3_m1 - e_1 * fs_8505_8_33 * r_2 * h1_m1 + e_2 * fs_285_4_55 * h5_m1 - e_2 * fs_8925_16_22 * r_2 * h3_m1 + e_2 * fs_6075_8_33 * r_4 * h1_m1 + e_3 * fs_1575_572_231 * h7_m1 - e_3 * fs_570_13_55 * r_2 * h5_m1 + e_3 * fs_8925_44_22 * r_4 * h3_m1 - e_3 * fs_225_1_33 * r_6 * h1_m1 + e_4 * fs_483_4862_165 * h9_m1 - e_4 * fs_4725_4862_231 * r_2 * h7_m1 + e_4 * fs_114_13_55 * r_4 * h5_m1 - e_4 * fs_8925_286_22 * r_6 * h3_m1 + e_4 * fs_675_22_33 * r_8 * h1_m1 + e_5 * fs_75_4199_323323 * h11_m11 + e_5 * fs_75_8398_2 * h11_m1 - e_5 * fs_46_2431_165 * r_2 * h9_m1 + e_5 * fs_4725_46189_231 * r_4 * h7_m1 - e_5 * fs_152_221_55 * r_6 * h5_m1 + e_5 * fs_595_286_22 * r_8 * h3_m1 - e_5 * fs_270_143_33 * r_10 * h1_m1 - e_6 * fs_6_4199_323323 * r_2 * h11_m11 - e_6 * fs_3_4199_2 * r_2 * h11_m1 + e_6 * fs_2_2431_165 * r_4 * h9_m1 - e_6 * fs_150_46189_231 * r_6 * h7_m1 + e_6 * fs_4_221_55 * r_8 * h5_m1 - e_6 * fs_7_143_22 * r_10 * h3_m1 + e_6 * fs_6_143_33 * r_12 * h1_m1;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m10, ph11_m2, ab_2, pc_133 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m2 = ph11_m2[k];

        pc_133[k] = - e_1 * fs_16065_32_22 * h3_m2 - e_2 * fs_285_4_154 * h5_m2 + e_2 * fs_8925_16_22 * r_2 * h3_m2 - e_3 * fs_2835_572_385 * h7_m2 + e_3 * fs_570_13_154 * r_2 * h5_m2 - e_3 * fs_8925_44_22 * r_4 * h3_m2 - e_4 * fs_483_221_3 * h9_m2 + e_4 * fs_8505_4862_385 * r_2 * h7_m2 - e_4 * fs_114_13_154 * r_4 * h5_m2 + e_4 * fs_8925_286_22 * r_6 * h3_m2 + e_5 * fs_75_4199_146965 * h11_m10 - e_5 * fs_75_8398_26 * h11_m2 + e_5 * fs_92_221_3 * r_2 * h9_m2 - e_5 * fs_8505_46189_385 * r_4 * h7_m2 + e_5 * fs_152_221_154 * r_6 * h5_m2 - e_5 * fs_595_286_22 * r_8 * h3_m2 - e_6 * fs_6_4199_146965 * r_2 * h11_m10 + e_6 * fs_3_4199_26 * r_2 * h11_m2 - e_6 * fs_4_221_3 * r_4 * h9_m2 + e_6 * fs_270_46189_385 * r_6 * h7_m2 - e_6 * fs_4_221_154 * r_8 * h5_m2 + e_6 * fs_7_143_22 * r_10 * h3_m2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_m3, ph5_m3, ph7_m3, ph9_m9, ph9_m3, ph11_m9, ph11_m3, ab_2, pc_134 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m3 = ph11_m3[k];

        pc_134[k] = e_1 * fs_5355_32_66 * h3_m3 + e_2 * fs_95_2_462 * h5_m3 - e_2 * fs_2975_16_66 * r_2 * h3_m3 + e_3 * fs_4725_572_385 * h7_m3 - e_3 * fs_380_13_462 * r_2 * h5_m3 + e_3 * fs_2975_44_66 * r_4 * h3_m3 + e_4 * fs_483_442_1326 * h9_m9 + e_4 * fs_483_221_14 * h9_m3 - e_4 * fs_14175_4862_385 * r_2 * h7_m3 + e_4 * fs_76_13_462 * r_4 * h5_m3 - e_4 * fs_2975_286_66 * r_6 * h3_m3 + e_5 * fs_75_4199_62985 * h11_m9 + e_5 * fs_75_8398_182 * h11_m3 - e_5 * fs_46_221_1326 * r_2 * h9_m9 - e_5 * fs_92_221_14 * r_2 * h9_m3 + e_5 * fs_14175_46189_385 * r_4 * h7_m3 - e_5 * fs_304_663_462 * r_6 * h5_m3 + e_5 * fs_595_858_66 * r_8 * h3_m3 - e_6 * fs_6_4199_62985 * r_2 * h11_m9 - e_6 * fs_3_4199_182 * r_2 * h11_m3 + e_6 * fs_2_221_1326 * r_4 * h9_m9 + e_6 * fs_4_221_14 * r_4 * h9_m3 - e_6 * fs_450_46189_385 * r_6 * h7_m3 + e_6 * fs_8_663_462 * r_8 * h5_m3 - e_6 * fs_7_429_66 * r_10 * h3_m3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m4, ph7_m4, ph9_m8, ph9_m4, ph11_m8, ph11_m4, ab_2, pc_135 : simd::cache_line_size())
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

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m4 = ph11_m4[k];

        pc_135[k] = - e_2 * fs_285_4_154 * h5_m4 - e_3 * fs_1575_104_210 * h7_m4 + e_3 * fs_570_13_154 * r_2 * h5_m4 + e_4 * fs_483_221_442 * h9_m8 - e_4 * fs_483_442_182 * h9_m4 + e_4 * fs_4725_884_210 * r_2 * h7_m4 - e_4 * fs_114_13_154 * r_4 * h5_m4 + e_5 * fs_75_4199_25194 * h11_m8 - e_5 * fs_75_8398_910 * h11_m4 - e_5 * fs_92_221_442 * r_2 * h9_m8 + e_5 * fs_46_221_182 * r_2 * h9_m4 - e_5 * fs_4725_8398_210 * r_4 * h7_m4 + e_5 * fs_152_221_154 * r_6 * h5_m4 - e_6 * fs_6_4199_25194 * r_2 * h11_m8 + e_6 * fs_3_4199_910 * r_2 * h11_m4 + e_6 * fs_4_221_442 * r_4 * h9_m8 - e_6 * fs_2_221_182 * r_4 * h9_m4 + e_6 * fs_75_4199_210 * r_6 * h7_m4 - e_6 * fs_4_221_154 * r_8 * h5_m4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_m5, ph7_m7, ph7_m5, ph7_p6, ph9_m7, ph9_m5, ph9_p6, ph11_m7, ph11_m5, ph11_p6, ab_2, pc_136, pc_137 : simd::cache_line_size())
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

        const auto h5_m5 = ph5_m5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_p6 = ph11_p6[k];

        pc_136[k] = e_2 * fs_285_4_55 * h5_m5 + e_3 * fs_315_104_2730 * h7_m7 + e_3 * fs_4725_104_30 * h7_m5 - e_3 * fs_570_13_55 * r_2 * h5_m5 + e_4 * fs_966_221_91 * h9_m7 + e_4 * fs_483_442_455 * h9_m5 - e_4 * fs_945_884_2730 * r_2 * h7_m7 - e_4 * fs_14175_884_30 * r_2 * h7_m5 + e_4 * fs_114_13_55 * r_4 * h5_m5 + e_5 * fs_75_4199_9282 * h11_m7 + e_5 * fs_75_4199_910 * h11_m5 - e_5 * fs_184_221_91 * r_2 * h9_m7 - e_5 * fs_46_221_455 * r_2 * h9_m5 + e_5 * fs_945_8398_2730 * r_4 * h7_m7 + e_5 * fs_14175_8398_30 * r_4 * h7_m5 - e_5 * fs_152_221_55 * r_6 * h5_m5 - e_6 * fs_6_4199_9282 * r_2 * h11_m7 - e_6 * fs_6_4199_910 * r_2 * h11_m5 + e_6 * fs_8_221_91 * r_4 * h9_m7 + e_6 * fs_2_221_455 * r_4 * h9_m5 - e_6 * fs_15_4199_2730 * r_6 * h7_m7 - e_6 * fs_225_4199_30 * r_6 * h7_m5 + e_6 * fs_4_221_55 * r_8 * h5_m5;

        pc_137[k] = e_3 * fs_4725_52_13 * h7_p6 + e_4 * fs_483_221_455 * h9_p6 - e_4 * fs_14175_442_13 * r_2 * h7_p6 + e_5 * fs_150_4199_1547 * h11_p6 - e_5 * fs_92_221_455 * r_2 * h9_p6 + e_5 * fs_14175_4199_13 * r_4 * h7_p6 - e_6 * fs_12_4199_1547 * r_2 * h11_p6 + e_6 * fs_4_221_455 * r_4 * h9_p6 - e_6 * fs_450_4199_13 * r_6 * h7_p6;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p5, ph7_p5, ph7_p7, ph9_p5, ph9_p7, ph11_p5, ph11_p7, ab_2, pc_138 : simd::cache_line_size())
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

        const auto h5_p5 = ph5_p5[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p5 = ph11_p5[k];
        const auto h11_p7 = ph11_p7[k];

        pc_138[k] = - e_2 * fs_285_4_55 * h5_p5 - e_3 * fs_4725_104_30 * h7_p5 + e_3 * fs_315_104_2730 * h7_p7 + e_3 * fs_570_13_55 * r_2 * h5_p5 - e_4 * fs_483_442_455 * h9_p5 + e_4 * fs_966_221_91 * h9_p7 + e_4 * fs_14175_884_30 * r_2 * h7_p5 - e_4 * fs_945_884_2730 * r_2 * h7_p7 - e_4 * fs_114_13_55 * r_4 * h5_p5 - e_5 * fs_75_4199_910 * h11_p5 + e_5 * fs_75_4199_9282 * h11_p7 + e_5 * fs_46_221_455 * r_2 * h9_p5 - e_5 * fs_184_221_91 * r_2 * h9_p7 - e_5 * fs_14175_8398_30 * r_4 * h7_p5 + e_5 * fs_945_8398_2730 * r_4 * h7_p7 + e_5 * fs_152_221_55 * r_6 * h5_p5 + e_6 * fs_6_4199_910 * r_2 * h11_p5 - e_6 * fs_6_4199_9282 * r_2 * h11_p7 - e_6 * fs_2_221_455 * r_4 * h9_p5 + e_6 * fs_8_221_91 * r_4 * h9_p7 + e_6 * fs_225_4199_30 * r_6 * h7_p5 - e_6 * fs_15_4199_2730 * r_6 * h7_p7 - e_6 * fs_4_221_55 * r_8 * h5_p5;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph5_p4, ph7_p4, ph9_p4, ph9_p8, ph11_p4, ph11_p8, ab_2, pc_139 : simd::cache_line_size())
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

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p8 = ph11_p8[k];

        pc_139[k] = e_2 * fs_285_4_154 * h5_p4 + e_3 * fs_1575_104_210 * h7_p4 - e_3 * fs_570_13_154 * r_2 * h5_p4 + e_4 * fs_483_442_182 * h9_p4 + e_4 * fs_483_221_442 * h9_p8 - e_4 * fs_4725_884_210 * r_2 * h7_p4 + e_4 * fs_114_13_154 * r_4 * h5_p4 + e_5 * fs_75_8398_910 * h11_p4 + e_5 * fs_75_4199_25194 * h11_p8 - e_5 * fs_46_221_182 * r_2 * h9_p4 - e_5 * fs_92_221_442 * r_2 * h9_p8 + e_5 * fs_4725_8398_210 * r_4 * h7_p4 - e_5 * fs_152_221_154 * r_6 * h5_p4 - e_6 * fs_3_4199_910 * r_2 * h11_p4 - e_6 * fs_6_4199_25194 * r_2 * h11_p8 + e_6 * fs_2_221_182 * r_4 * h9_p4 + e_6 * fs_4_221_442 * r_4 * h9_p8 - e_6 * fs_75_4199_210 * r_6 * h7_p4 + e_6 * fs_4_221_154 * r_8 * h5_p4;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph9_p9, ph11_p3, ph11_p9, ab_2, pc_140 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p9 = ph11_p9[k];

        pc_140[k] = - e_1 * fs_5355_32_66 * h3_p3 - e_2 * fs_95_2_462 * h5_p3 + e_2 * fs_2975_16_66 * r_2 * h3_p3 - e_3 * fs_4725_572_385 * h7_p3 + e_3 * fs_380_13_462 * r_2 * h5_p3 - e_3 * fs_2975_44_66 * r_4 * h3_p3 - e_4 * fs_483_221_14 * h9_p3 + e_4 * fs_483_442_1326 * h9_p9 + e_4 * fs_14175_4862_385 * r_2 * h7_p3 - e_4 * fs_76_13_462 * r_4 * h5_p3 + e_4 * fs_2975_286_66 * r_6 * h3_p3 - e_5 * fs_75_8398_182 * h11_p3 + e_5 * fs_75_4199_62985 * h11_p9 + e_5 * fs_92_221_14 * r_2 * h9_p3 - e_5 * fs_46_221_1326 * r_2 * h9_p9 - e_5 * fs_14175_46189_385 * r_4 * h7_p3 + e_5 * fs_304_663_462 * r_6 * h5_p3 - e_5 * fs_595_858_66 * r_8 * h3_p3 + e_6 * fs_3_4199_182 * r_2 * h11_p3 - e_6 * fs_6_4199_62985 * r_2 * h11_p9 - e_6 * fs_4_221_14 * r_4 * h9_p3 + e_6 * fs_2_221_1326 * r_4 * h9_p9 + e_6 * fs_450_46189_385 * r_6 * h7_p3 - e_6 * fs_8_663_462 * r_8 * h5_p3 + e_6 * fs_7_429_66 * r_10 * h3_p3;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ph11_p10, ab_2, pc_141 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p10 = ph11_p10[k];

        pc_141[k] = e_1 * fs_16065_32_22 * h3_p2 + e_2 * fs_285_4_154 * h5_p2 - e_2 * fs_8925_16_22 * r_2 * h3_p2 + e_3 * fs_2835_572_385 * h7_p2 - e_3 * fs_570_13_154 * r_2 * h5_p2 + e_3 * fs_8925_44_22 * r_4 * h3_p2 + e_4 * fs_483_221_3 * h9_p2 - e_4 * fs_8505_4862_385 * r_2 * h7_p2 + e_4 * fs_114_13_154 * r_4 * h5_p2 - e_4 * fs_8925_286_22 * r_6 * h3_p2 + e_5 * fs_75_8398_26 * h11_p2 + e_5 * fs_75_4199_146965 * h11_p10 - e_5 * fs_92_221_3 * r_2 * h9_p2 + e_5 * fs_8505_46189_385 * r_4 * h7_p2 - e_5 * fs_152_221_154 * r_6 * h5_p2 + e_5 * fs_595_286_22 * r_8 * h3_p2 - e_6 * fs_3_4199_26 * r_2 * h11_p2 - e_6 * fs_6_4199_146965 * r_2 * h11_p10 + e_6 * fs_4_221_3 * r_4 * h9_p2 - e_6 * fs_270_46189_385 * r_6 * h7_p2 + e_6 * fs_4_221_154 * r_8 * h5_p2 - e_6 * fs_7_143_22 * r_10 * h3_p2;
    }

    // NOTE: the angular components are formed in 129 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph11_p1, ph11_p11, ab_2, pc_142 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p11 = ph11_p11[k];

        pc_142[k] = - e_0 * fs_14175_32_33 * h1_p1 - e_1 * fs_16065_32_22 * h3_p1 + e_1 * fs_8505_8_33 * r_2 * h1_p1 - e_2 * fs_285_4_55 * h5_p1 + e_2 * fs_8925_16_22 * r_2 * h3_p1 - e_2 * fs_6075_8_33 * r_4 * h1_p1 - e_3 * fs_1575_572_231 * h7_p1 + e_3 * fs_570_13_55 * r_2 * h5_p1 - e_3 * fs_8925_44_22 * r_4 * h3_p1 + e_3 * fs_225_1_33 * r_6 * h1_p1 - e_4 * fs_483_4862_165 * h9_p1 + e_4 * fs_4725_4862_231 * r_2 * h7_p1 - e_4 * fs_114_13_55 * r_4 * h5_p1 + e_4 * fs_8925_286_22 * r_6 * h3_p1 - e_4 * fs_675_22_33 * r_8 * h1_p1 - e_5 * fs_75_8398_2 * h11_p1 + e_5 * fs_75_4199_323323 * h11_p11 + e_5 * fs_46_2431_165 * r_2 * h9_p1 - e_5 * fs_4725_46189_231 * r_4 * h7_p1 + e_5 * fs_152_221_55 * r_6 * h5_p1 - e_5 * fs_595_286_22 * r_8 * h3_p1 + e_5 * fs_270_143_33 * r_10 * h1_p1 + e_6 * fs_3_4199_2 * r_2 * h11_p1 - e_6 * fs_6_4199_323323 * r_2 * h11_p11 - e_6 * fs_2_2431_165 * r_4 * h9_p1 + e_6 * fs_150_46189_231 * r_6 * h7_p1 - e_6 * fs_4_221_55 * r_8 * h5_p1 + e_6 * fs_7_143_22 * r_10 * h3_p1 - e_6 * fs_6_143_33 * r_12 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[143] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142};

    for (size_t n = 0; n < 143; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
