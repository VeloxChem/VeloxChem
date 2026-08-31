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



#include "SimdOverlapRecIH.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ih_overlap(double                         *values,
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
            false, std::string("SimdOverlapRecIH.compute_ih_overlap: Basis functions must be of angular momenta six and five"));
    }

    if (harmonics.size() < 11)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIH.compute_ih_overlap: Harmonics must reach angular momentum eleven"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIH.compute_ih_overlap: Number of values exceeds number of atom pairs"));
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
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

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

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

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

    // accumulate the prefactor of each term over the pairs of primitives

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

            const auto f_0 = fbase * bexp * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * bexp * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * bexp * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * bexp * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_4 = fbase * bexp * fmu * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_5 = fbase * bexp / fexp / fexp / fexp / fexp / fexp / fexp;

            // NOTE: the exponential depends on the pair of primitives alone, so it is
            // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fss = std::exp(-fmu * ab_2[k]);

                pe_0[k] += f_0 * fss;
                pe_1[k] += f_1 * fss;
                pe_2[k] += f_2 * fss;
                pe_3[k] += f_3 * fss;
                pe_4[k] += f_4 * fss;
                pe_5[k] += f_5 * fss;
            }
        }
    }

    // NOTE: the geometry of a term is a solid harmonic of the vector between the
    // atoms times a power of their squared distance.

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

    // NOTE: the rows of the values are not aligned, as they start at the offset
    // of this combination of basis functions in the values block, so they are kept
    // out of the aligned clauses below.

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

    const auto fs_9823275_512 = std::sqrt(19186.083984375);
    const auto fs_29469825_256 = std::sqrt(115116.50390625);
    const auto fs_12375_16 = std::sqrt(773.4375);
    const auto fs_121275_8 = std::sqrt(15159.375);
    const auto fs_601425_16 = std::sqrt(37589.0625);
    const auto fs_118125_29744 = std::sqrt(118125.0 / 29744.0);
    const auto fs_111375_676 = std::sqrt(111375.0 / 676.0);
    const auto fs_99225_88 = std::sqrt(99225.0 / 88.0);
    const auto fs_7425_4 = std::sqrt(1856.25);
    const auto fs_6615_2149004 = std::sqrt(6615.0 / 2149004.0);
    const auto fs_118125_537251 = std::sqrt(118125.0 / 537251.0);
    const auto fs_495_169 = std::sqrt(495.0 / 169.0);
    const auto fs_22050_1859 = std::sqrt(22050.0 / 1859.0);
    const auto fs_675_44 = std::sqrt(675.0 / 44.0);
    const auto fs_9_35263202 = std::sqrt(9.0 / 35263202.0);
    const auto fs_693_4199 = std::sqrt(693.0 / 4199.0);
    const auto fs_15_537251 = std::sqrt(15.0 / 537251.0);
    const auto fs_118125_193947611 = std::sqrt(118125.0 / 193947611.0);
    const auto fs_220_48841 = std::sqrt(220.0 / 48841.0);
    const auto fs_49_3718 = std::sqrt(49.0 / 3718.0);
    const auto fs_27_1859 = std::sqrt(27.0 / 1859.0);
    const auto fs_29469825_1024 = std::sqrt(28779.1259765625);
    const auto fs_17325_8 = std::sqrt(2165.625);
    const auto fs_637875_29744 = std::sqrt(637875.0 / 29744.0);
    const auto fs_155925_338 = std::sqrt(155925.0 / 338.0);
    const auto fs_1323_48841 = std::sqrt(1323.0 / 48841.0);
    const auto fs_637875_537251 = std::sqrt(637875.0 / 537251.0);
    const auto fs_1386_169 = std::sqrt(1386.0 / 169.0);
    const auto fs_9_2712554 = std::sqrt(9.0 / 2712554.0);
    const auto fs_315_4199 = std::sqrt(315.0 / 4199.0);
    const auto fs_12_48841 = std::sqrt(12.0 / 48841.0);
    const auto fs_637875_193947611 = std::sqrt(637875.0 / 193947611.0);
    const auto fs_616_48841 = std::sqrt(616.0 / 48841.0);
    const auto fs_3274425_512 = std::sqrt(6395.361328125);
    const auto fs_5775_2 = std::sqrt(2887.5);
    const auto fs_40425_8 = std::sqrt(5053.125);
    const auto fs_1771875_29744 = std::sqrt(1771875.0 / 29744.0);
    const auto fs_103950_169 = std::sqrt(103950.0 / 169.0);
    const auto fs_33075_88 = std::sqrt(33075.0 / 88.0);
    const auto fs_6174_48841 = std::sqrt(6174.0 / 48841.0);
    const auto fs_1323_442 = std::sqrt(1323.0 / 442.0);
    const auto fs_1771875_537251 = std::sqrt(1771875.0 / 537251.0);
    const auto fs_1848_169 = std::sqrt(1848.0 / 169.0);
    const auto fs_7350_1859 = std::sqrt(7350.0 / 1859.0);
    const auto fs_63_2712554 = std::sqrt(63.0 / 2712554.0);
    const auto fs_135_4199 = std::sqrt(135.0 / 4199.0);
    const auto fs_56_48841 = std::sqrt(56.0 / 48841.0);
    const auto fs_6_221 = std::sqrt(6.0 / 221.0);
    const auto fs_1771875_193947611 = std::sqrt(1771875.0 / 193947611.0);
    const auto fs_2464_146523 = std::sqrt(2464.0 / 146523.0);
    const auto fs_49_11154 = std::sqrt(49.0 / 11154.0);
    const auto fs_590625_5408 = std::sqrt(590625.0 / 5408.0);
    const auto fs_3087_7514 = std::sqrt(3087.0 / 7514.0);
    const auto fs_882_221 = std::sqrt(882.0 / 221.0);
    const auto fs_590625_97682 = std::sqrt(590625.0 / 97682.0);
    const auto fs_315_2712554 = std::sqrt(315.0 / 2712554.0);
    const auto fs_54_4199 = std::sqrt(54.0 / 4199.0);
    const auto fs_14_3757 = std::sqrt(14.0 / 3757.0);
    const auto fs_8_221 = std::sqrt(8.0 / 221.0);
    const auto fs_590625_35263202 = std::sqrt(590625.0 / 35263202.0);
    const auto fs_759375_5408 = std::sqrt(759375.0 / 5408.0);
    const auto fs_23625_416 = std::sqrt(23625.0 / 416.0);
    const auto fs_15435_15028 = std::sqrt(15435.0 / 15028.0);
    const auto fs_12348_3757 = std::sqrt(12348.0 / 3757.0);
    const auto fs_759375_97682 = std::sqrt(759375.0 / 97682.0);
    const auto fs_23625_7514 = std::sqrt(23625.0 / 7514.0);
    const auto fs_630_1356277 = std::sqrt(630.0 / 1356277.0);
    const auto fs_378_79781 = std::sqrt(378.0 / 79781.0);
    const auto fs_35_3757 = std::sqrt(35.0 / 3757.0);
    const auto fs_112_3757 = std::sqrt(112.0 / 3757.0);
    const auto fs_759375_35263202 = std::sqrt(759375.0 / 35263202.0);
    const auto fs_23625_2712554 = std::sqrt(23625.0 / 2712554.0);
    const auto fs_50625_208 = std::sqrt(50625.0 / 208.0);
    const auto fs_15435_3757 = std::sqrt(15435.0 / 3757.0);
    const auto fs_50625_3757 = std::sqrt(50625.0 / 3757.0);
    const auto fs_252_79781 = std::sqrt(252.0 / 79781.0);
    const auto fs_140_3757 = std::sqrt(140.0 / 3757.0);
    const auto fs_50625_1356277 = std::sqrt(50625.0 / 1356277.0);
    const auto fs_9823275_256 = std::sqrt(38372.16796875);
    const auto fs_61875_16 = std::sqrt(3867.1875);
    const auto fs_121275_4 = std::sqrt(30318.75);
    const auto fs_200475_16 = std::sqrt(12529.6875);
    const auto fs_275625_7436 = std::sqrt(275625.0 / 7436.0);
    const auto fs_556875_676 = std::sqrt(556875.0 / 676.0);
    const auto fs_99225_44 = std::sqrt(99225.0 / 44.0);
    const auto fs_2475_4 = std::sqrt(618.75);
    const auto fs_99225_2149004 = std::sqrt(99225.0 / 2149004.0);
    const auto fs_1102500_537251 = std::sqrt(1102500.0 / 537251.0);
    const auto fs_2475_169 = std::sqrt(2475.0 / 169.0);
    const auto fs_44100_1859 = std::sqrt(44100.0 / 1859.0);
    const auto fs_225_44 = std::sqrt(225.0 / 44.0);
    const auto fs_99_17631601 = std::sqrt(99.0 / 17631601.0);
    const auto fs_378_4199 = std::sqrt(378.0 / 4199.0);
    const auto fs_225_537251 = std::sqrt(225.0 / 537251.0);
    const auto fs_1102500_193947611 = std::sqrt(1102500.0 / 193947611.0);
    const auto fs_1100_48841 = std::sqrt(1100.0 / 48841.0);
    const auto fs_49_1859 = std::sqrt(49.0 / 1859.0);
    const auto fs_9_1859 = std::sqrt(9.0 / 1859.0);
    const auto fs_9823275_1024 = std::sqrt(9593.0419921875);
    const auto fs_49116375_512 = std::sqrt(95930.419921875);
    const auto fs_66825_32 = std::sqrt(2088.28125);
    const auto fs_1002375_32 = std::sqrt(31324.21875);
    const auto fs_86625_1352 = std::sqrt(86625.0 / 1352.0);
    const auto fs_601425_1352 = std::sqrt(601425.0 / 1352.0);
    const auto fs_12375_8 = std::sqrt(1546.875);
    const auto fs_3969_25432 = std::sqrt(3969.0 / 25432.0);
    const auto fs_3969_884 = std::sqrt(3969.0 / 884.0);
    const auto fs_173250_48841 = std::sqrt(173250.0 / 48841.0);
    const auto fs_2673_338 = std::sqrt(2673.0 / 338.0);
    const auto fs_1125_88 = std::sqrt(1125.0 / 88.0);
    const auto fs_540_17631601 = std::sqrt(540.0 / 17631601.0);
    const auto fs_360_4199 = std::sqrt(360.0 / 4199.0);
    const auto fs_9_6358 = std::sqrt(9.0 / 6358.0);
    const auto fs_9_221 = std::sqrt(9.0 / 221.0);
    const auto fs_173250_17631601 = std::sqrt(173250.0 / 17631601.0);
    const auto fs_594_48841 = std::sqrt(594.0 / 48841.0);
    const auto fs_45_3718 = std::sqrt(45.0 / 3718.0);
    const auto fs_49116375_2048 = std::sqrt(23982.60498046875);
    const auto fs_3274425_256 = std::sqrt(12790.72265625);
    const auto fs_5775_16 = std::sqrt(360.9375);
    const auto fs_40425_4 = std::sqrt(10106.25);
    const auto fs_189000_1859 = std::sqrt(189000.0 / 1859.0);
    const auto fs_51975_676 = std::sqrt(51975.0 / 676.0);
    const auto fs_33075_44 = std::sqrt(33075.0 / 44.0);
    const auto fs_53361_97682 = std::sqrt(53361.0 / 97682.0);
    const auto fs_441_884 = std::sqrt(441.0 / 884.0);
    const auto fs_3024000_537251 = std::sqrt(3024000.0 / 537251.0);
    const auto fs_231_169 = std::sqrt(231.0 / 169.0);
    const auto fs_14700_1859 = std::sqrt(14700.0 / 1859.0);
    const auto fs_243_1356277 = std::sqrt(243.0 / 1356277.0);
    const auto fs_243_4199 = std::sqrt(243.0 / 4199.0);
    const auto fs_242_48841 = std::sqrt(242.0 / 48841.0);
    const auto fs_1_221 = std::sqrt(1.0 / 221.0);
    const auto fs_3024000_193947611 = std::sqrt(3024000.0 / 193947611.0);
    const auto fs_308_146523 = std::sqrt(308.0 / 146523.0);
    const auto fs_49_5577 = std::sqrt(49.0 / 5577.0);
    const auto fs_4921875_59488 = std::sqrt(4921875.0 / 59488.0);
    const auto fs_55125_416 = std::sqrt(55125.0 / 416.0);
    const auto fs_250047_195364 = std::sqrt(250047.0 / 195364.0);
    const auto fs_1323_3757 = std::sqrt(1323.0 / 3757.0);
    const auto fs_4921875_1074502 = std::sqrt(4921875.0 / 1074502.0);
    const auto fs_55125_7514 = std::sqrt(55125.0 / 7514.0);
    const auto fs_1008_1356277 = std::sqrt(1008.0 / 1356277.0);
    const auto fs_2592_79781 = std::sqrt(2592.0 / 79781.0);
    const auto fs_567_48841 = std::sqrt(567.0 / 48841.0);
    const auto fs_12_3757 = std::sqrt(12.0 / 3757.0);
    const auto fs_4921875_387895222 = std::sqrt(4921875.0 / 387895222.0);
    const auto fs_55125_2712554 = std::sqrt(55125.0 / 2712554.0);
    const auto fs_28125_1352 = std::sqrt(28125.0 / 1352.0);
    const auto fs_1125_13 = std::sqrt(1125.0 / 13.0);
    const auto fs_64827_30056 = std::sqrt(64827.0 / 30056.0);
    const auto fs_27783_15028 = std::sqrt(27783.0 / 15028.0);
    const auto fs_56250_48841 = std::sqrt(56250.0 / 48841.0);
    const auto fs_18000_3757 = std::sqrt(18000.0 / 3757.0);
    const auto fs_6615_2712554 = std::sqrt(6615.0 / 2712554.0);
    const auto fs_1260_79781 = std::sqrt(1260.0 / 79781.0);
    const auto fs_147_7514 = std::sqrt(147.0 / 7514.0);
    const auto fs_63_3757 = std::sqrt(63.0 / 3757.0);
    const auto fs_56250_17631601 = std::sqrt(56250.0 / 17631601.0);
    const auto fs_18000_1356277 = std::sqrt(18000.0 / 1356277.0);
    const auto fs_16875_1352 = std::sqrt(16875.0 / 1352.0);
    const auto fs_77175_15028 = std::sqrt(77175.0 / 15028.0);
    const auto fs_33750_48841 = std::sqrt(33750.0 / 48841.0);
    const auto fs_18144_1356277 = std::sqrt(18144.0 / 1356277.0);
    const auto fs_175_3757 = std::sqrt(175.0 / 3757.0);
    const auto fs_33750_17631601 = std::sqrt(33750.0 / 17631601.0);
    const auto fs_2679075_256 = std::sqrt(10465.13671875);
    const auto fs_893025_512 = std::sqrt(1744.189453125);
    const auto fs_84375_32 = std::sqrt(2636.71875);
    const auto fs_33075_4 = std::sqrt(8268.75);
    const auto fs_18225_32 = std::sqrt(569.53125);
    const auto fs_1929375_40898 = std::sqrt(1929375.0 / 40898.0);
    const auto fs_759375_1352 = std::sqrt(759375.0 / 1352.0);
    const auto fs_297675_484 = std::sqrt(297675.0 / 484.0);
    const auto fs_225_8 = std::sqrt(28.125);
    const auto fs_4465125_47278088 = std::sqrt(4465125.0 / 47278088.0);
    const auto fs_19845_9724 = std::sqrt(19845.0 / 9724.0);
    const auto fs_15435000_5909761 = std::sqrt(15435000.0 / 5909761.0);
    const auto fs_3375_338 = std::sqrt(3375.0 / 338.0);
    const auto fs_132300_20449 = std::sqrt(132300.0 / 20449.0);
    const auto fs_225_968 = std::sqrt(225.0 / 968.0);
    const auto fs_297_17631601 = std::sqrt(297.0 / 17631601.0);
    const auto fs_198_4199 = std::sqrt(198.0 / 4199.0);
    const auto fs_10125_11819522 = std::sqrt(10125.0 / 11819522.0);
    const auto fs_45_2431 = std::sqrt(45.0 / 2431.0);
    const auto fs_15435000_2133423721 = std::sqrt(15435000.0 / 2133423721.0);
    const auto fs_750_48841 = std::sqrt(750.0 / 48841.0);
    const auto fs_147_20449 = std::sqrt(147.0 / 20449.0);
    const auto fs_9_40898 = std::sqrt(9.0 / 40898.0);
    const auto fs_893025_2048 = std::sqrt(436.04736328125);
    const auto fs_4465125_256 = std::sqrt(17441.89453125);
    const auto fs_4465125_64 = std::sqrt(69767.578125);
    const auto fs_1125 = std::sqrt(1125.0);
    const auto fs_55125_4 = std::sqrt(13781.25);
    const auto fs_91125_4 = std::sqrt(22781.25);
    const auto fs_15931125_81796 = std::sqrt(15931125.0 / 81796.0);
    const auto fs_40500_169 = std::sqrt(40500.0 / 169.0);
    const auto fs_496125_484 = std::sqrt(496125.0 / 484.0);
    const auto fs_19845_20449 = std::sqrt(19845.0 / 20449.0);
    const auto fs_3969_2431 = std::sqrt(3969.0 / 2431.0);
    const auto fs_220500_20449 = std::sqrt(220500.0 / 20449.0);
    const auto fs_720_169 = std::sqrt(720.0 / 169.0);
    const auto fs_1125_121 = std::sqrt(1125.0 / 121.0);
    const auto fs_5445_17631601 = std::sqrt(5445.0 / 17631601.0);
    const auto fs_297_4199 = std::sqrt(297.0 / 4199.0);
    const auto fs_180_20449 = std::sqrt(180.0 / 20449.0);
    const auto fs_36_2431 = std::sqrt(36.0 / 2431.0);
    const auto fs_220500_7382089 = std::sqrt(220500.0 / 7382089.0);
    const auto fs_320_48841 = std::sqrt(320.0 / 48841.0);
    const auto fs_245_20449 = std::sqrt(245.0 / 20449.0);
    const auto fs_1488375_256 = std::sqrt(5813.96484375);
    const auto fs_40186125_512 = std::sqrt(78488.525390625);
    const auto fs_12675_32 = std::sqrt(396.09375);
    const auto fs_18375_4 = std::sqrt(4593.75);
    const auto fs_820125_32 = std::sqrt(25628.90625);
    const auto fs_70875_968 = std::sqrt(70875.0 / 968.0);
    const auto fs_165375_1144 = std::sqrt(165375.0 / 1144.0);
    const auto fs_675_8 = std::sqrt(84.375);
    const auto fs_165375_484 = std::sqrt(165375.0 / 484.0);
    const auto fs_10125_8 = std::sqrt(1265.625);
    const auto fs_59397849_47278088 = std::sqrt(59397849.0 / 47278088.0);
    const auto fs_370881_165308 = std::sqrt(370881.0 / 165308.0);
    const auto fs_141750_34969 = std::sqrt(141750.0 / 34969.0);
    const auto fs_330750_41327 = std::sqrt(330750.0 / 41327.0);
    const auto fs_3_2 = std::sqrt(1.5);
    const auto fs_73500_20449 = std::sqrt(73500.0 / 20449.0);
    const auto fs_10125_968 = std::sqrt(10125.0 / 968.0);
    const auto fs_13365_17631601 = std::sqrt(13365.0 / 17631601.0);
    const auto fs_5346_79781 = std::sqrt(5346.0 / 79781.0);
    const auto fs_134689_11819522 = std::sqrt(134689.0 / 11819522.0);
    const auto fs_841_41327 = std::sqrt(841.0 / 41327.0);
    const auto fs_141750_12623809 = std::sqrt(141750.0 / 12623809.0);
    const auto fs_330750_14919047 = std::sqrt(330750.0 / 14919047.0);
    const auto fs_2_867 = std::sqrt(2.0 / 867.0);
    const auto fs_245_61347 = std::sqrt(245.0 / 61347.0);
    const auto fs_405_40898 = std::sqrt(405.0 / 40898.0);
    const auto fs_40186125_2048 = std::sqrt(19622.13134765625);
    const auto f_945_16 = 59.0625;
    const auto fs_1575 = std::sqrt(1575.0);
    const auto f_105_2 = 52.5;
    const auto fs_6622875_654368 = std::sqrt(6622875.0 / 654368.0);
    const auto fs_7875_4576 = std::sqrt(7875.0 / 4576.0);
    const auto fs_56700_169 = std::sqrt(56700.0 / 169.0);
    const auto f_315_22 = 315.0 / 22.0;
    const auto fs_2223963_1074502 = std::sqrt(2223963.0 / 1074502.0);
    const auto fs_35721_82654 = std::sqrt(35721.0 / 82654.0);
    const auto fs_6622875_11819522 = std::sqrt(6622875.0 / 11819522.0);
    const auto fs_7875_82654 = std::sqrt(7875.0 / 82654.0);
    const auto fs_1008_169 = std::sqrt(1008.0 / 169.0);
    const auto f_210_143 = 210.0 / 143.0;
    const auto fs_3564_1356277 = std::sqrt(3564.0 / 1356277.0);
    const auto fs_3960_79781 = std::sqrt(3960.0 / 79781.0);
    const auto fs_10086_537251 = std::sqrt(10086.0 / 537251.0);
    const auto fs_162_41327 = std::sqrt(162.0 / 41327.0);
    const auto fs_6622875_4266847442 = std::sqrt(6622875.0 / 4266847442.0);
    const auto fs_7875_29838094 = std::sqrt(7875.0 / 29838094.0);
    const auto fs_448_48841 = std::sqrt(448.0 / 48841.0);
    const auto f_7_143 = 7.0 / 143.0;
    const auto fs_2083725_128 = std::sqrt(16279.1015625);
    const auto fs_25725_2 = std::sqrt(12862.5);
    const auto fs_28125_1936 = std::sqrt(28125.0 / 1936.0);
    const auto fs_1540125_29744 = std::sqrt(1540125.0 / 29744.0);
    const auto fs_231525_242 = std::sqrt(231525.0 / 242.0);
    const auto fs_9529569_4298008 = std::sqrt(9529569.0 / 4298008.0);
    const auto fs_46305_330616 = std::sqrt(46305.0 / 330616.0);
    const auto fs_28125_34969 = std::sqrt(28125.0 / 34969.0);
    const auto fs_1540125_537251 = std::sqrt(1540125.0 / 537251.0);
    const auto fs_205800_20449 = std::sqrt(205800.0 / 20449.0);
    const auto fs_9702_1356277 = std::sqrt(9702.0 / 1356277.0);
    const auto fs_41580_1356277 = std::sqrt(41580.0 / 1356277.0);
    const auto fs_21609_1074502 = std::sqrt(21609.0 / 1074502.0);
    const auto fs_105_82654 = std::sqrt(105.0 / 82654.0);
    const auto fs_28125_12623809 = std::sqrt(28125.0 / 12623809.0);
    const auto fs_1540125_193947611 = std::sqrt(1540125.0 / 193947611.0);
    const auto fs_686_61347 = std::sqrt(686.0 / 61347.0);
    const auto fs_270000_1859 = std::sqrt(270000.0 / 1859.0);
    const auto fs_108045_41327 = std::sqrt(108045.0 / 41327.0);
    const auto fs_4320000_537251 = std::sqrt(4320000.0 / 537251.0);
    const auto fs_43659_1356277 = std::sqrt(43659.0 / 1356277.0);
    const auto fs_980_41327 = std::sqrt(980.0 / 41327.0);
    const auto fs_4320000_193947611 = std::sqrt(4320000.0 / 193947611.0);
    const auto fs_39375_16 = std::sqrt(2460.9375);
    const auto fs_3472875_40898 = std::sqrt(3472875.0 / 40898.0);
    const auto fs_354375_676 = std::sqrt(354375.0 / 676.0);
    const auto fs_297675_1074502 = std::sqrt(297675.0 / 1074502.0);
    const auto fs_33075_9724 = std::sqrt(33075.0 / 9724.0);
    const auto fs_27783000_5909761 = std::sqrt(27783000.0 / 5909761.0);
    const auto fs_1575_169 = std::sqrt(1575.0 / 169.0);
    const auto fs_99_1356277 = std::sqrt(99.0 / 1356277.0);
    const auto fs_99_4199 = std::sqrt(99.0 / 4199.0);
    const auto fs_1350_537251 = std::sqrt(1350.0 / 537251.0);
    const auto fs_75_2431 = std::sqrt(75.0 / 2431.0);
    const auto fs_27783000_2133423721 = std::sqrt(27783000.0 / 2133423721.0);
    const auto fs_700_48841 = std::sqrt(700.0 / 48841.0);
    const auto f_945_8 = 118.125;
    const auto fs_2679075_512 = std::sqrt(5232.568359375);
    const auto fs_1125_32 = std::sqrt(35.15625);
    const auto f_105 = 105.0;
    const auto fs_54675_32 = std::sqrt(1708.59375);
    const auto fs_1852200_20449 = std::sqrt(1852200.0 / 20449.0);
    const auto fs_99225_1144 = std::sqrt(99225.0 / 1144.0);
    const auto fs_10125_1352 = std::sqrt(10125.0 / 1352.0);
    const auto f_315_11 = 315.0 / 11.0;
    const auto fs_50068935_47278088 = std::sqrt(50068935.0 / 47278088.0);
    const auto fs_6615_165308 = std::sqrt(6615.0 / 165308.0);
    const auto fs_29635200_5909761 = std::sqrt(29635200.0 / 5909761.0);
    const auto fs_198450_41327 = std::sqrt(198450.0 / 41327.0);
    const auto fs_45_338 = std::sqrt(45.0 / 338.0);
    const auto f_420_143 = 420.0 / 143.0;
    const auto fs_675_968 = std::sqrt(675.0 / 968.0);
    const auto fs_9900_17631601 = std::sqrt(9900.0 / 17631601.0);
    const auto fs_113535_11819522 = std::sqrt(113535.0 / 11819522.0);
    const auto fs_15_41327 = std::sqrt(15.0 / 41327.0);
    const auto fs_29635200_2133423721 = std::sqrt(29635200.0 / 2133423721.0);
    const auto fs_198450_14919047 = std::sqrt(198450.0 / 14919047.0);
    const auto fs_10_48841 = std::sqrt(10.0 / 48841.0);
    const auto f_14_143 = 14.0 / 143.0;
    const auto fs_27_40898 = std::sqrt(27.0 / 40898.0);
    const auto fs_2679075_2048 = std::sqrt(1308.14208984375);
    const auto fs_297675_256 = std::sqrt(1162.79296875);
    const auto fs_24111675_256 = std::sqrt(94186.23046875);
    const auto fs_46875_16 = std::sqrt(2929.6875);
    const auto fs_3675_4 = std::sqrt(918.75);
    const auto fs_492075_16 = std::sqrt(30754.6875);
    const auto fs_2679075_81796 = std::sqrt(2679075.0 / 81796.0);
    const auto fs_14175_286 = std::sqrt(14175.0 / 286.0);
    const auto fs_421875_676 = std::sqrt(421875.0 / 676.0);
    const auto fs_33075_484 = std::sqrt(33075.0 / 484.0);
    const auto fs_6075_4 = std::sqrt(1518.75);
    const auto fs_92907675_23639044 = std::sqrt(92907675.0 / 23639044.0);
    const auto fs_70560_41327 = std::sqrt(70560.0 / 41327.0);
    const auto fs_10716300_5909761 = std::sqrt(10716300.0 / 5909761.0);
    const auto fs_113400_41327 = std::sqrt(113400.0 / 41327.0);
    const auto fs_1875_169 = std::sqrt(1875.0 / 169.0);
    const auto fs_14700_20449 = std::sqrt(14700.0 / 20449.0);
    const auto fs_6075_484 = std::sqrt(6075.0 / 484.0);
    const auto fs_81675_17631601 = std::sqrt(81675.0 / 17631601.0);
    const auto fs_4950_79781 = std::sqrt(4950.0 / 79781.0);
    const auto fs_210675_5909761 = std::sqrt(210675.0 / 5909761.0);
    const auto fs_640_41327 = std::sqrt(640.0 / 41327.0);
    const auto fs_10716300_2133423721 = std::sqrt(10716300.0 / 2133423721.0);
    const auto fs_113400_14919047 = std::sqrt(113400.0 / 14919047.0);
    const auto fs_2500_146523 = std::sqrt(2500.0 / 146523.0);
    const auto fs_49_61347 = std::sqrt(49.0 / 61347.0);
    const auto fs_243_20449 = std::sqrt(243.0 / 20449.0);
    const auto fs_24111675_1024 = std::sqrt(23546.5576171875);
    const auto fs_8037225_128 = std::sqrt(62790.8203125);
    const auto fs_3375_8 = std::sqrt(421.875);
    const auto fs_164025_8 = std::sqrt(20503.125);
    const auto fs_7498575_654368 = std::sqrt(7498575.0 / 654368.0);
    const auto fs_3973725_59488 = std::sqrt(3973725.0 / 59488.0);
    const auto fs_30375_338 = std::sqrt(30375.0 / 338.0);
    const auto fs_2025_2 = std::sqrt(1012.5);
    const auto fs_25245045_11819522 = std::sqrt(25245045.0 / 11819522.0);
    const auto fs_275625_165308 = std::sqrt(275625.0 / 165308.0);
    const auto fs_7498575_11819522 = std::sqrt(7498575.0 / 11819522.0);
    const auto fs_3973725_1074502 = std::sqrt(3973725.0 / 1074502.0);
    const auto fs_270_169 = std::sqrt(270.0 / 169.0);
    const auto fs_2025_242 = std::sqrt(2025.0 / 242.0);
    const auto fs_118800_17631601 = std::sqrt(118800.0 / 17631601.0);
    const auto fs_79200_1356277 = std::sqrt(79200.0 / 1356277.0);
    const auto fs_114490_5909761 = std::sqrt(114490.0 / 5909761.0);
    const auto fs_625_41327 = std::sqrt(625.0 / 41327.0);
    const auto fs_7498575_4266847442 = std::sqrt(7498575.0 / 4266847442.0);
    const auto fs_3973725_387895222 = std::sqrt(3973725.0 / 387895222.0);
    const auto fs_120_48841 = std::sqrt(120.0 / 48841.0);
    const auto fs_162_20449 = std::sqrt(162.0 / 20449.0);
    const auto fs_8037225_512 = std::sqrt(15697.705078125);
    const auto fs_10800_169 = std::sqrt(10800.0 / 169.0);
    const auto fs_6075_14872 = std::sqrt(6075.0 / 14872.0);
    const auto fs_15435_12716 = std::sqrt(15435.0 / 12716.0);
    const auto fs_108045_330616 = std::sqrt(108045.0 / 330616.0);
    const auto fs_172800_48841 = std::sqrt(172800.0 / 48841.0);
    const auto fs_12150_537251 = std::sqrt(12150.0 / 537251.0);
    const auto fs_20790_1356277 = std::sqrt(20790.0 / 1356277.0);
    const auto fs_121275_2712554 = std::sqrt(121275.0 / 2712554.0);
    const auto fs_35_3179 = std::sqrt(35.0 / 3179.0);
    const auto fs_245_82654 = std::sqrt(245.0 / 82654.0);
    const auto fs_172800_17631601 = std::sqrt(172800.0 / 17631601.0);
    const auto fs_12150_193947611 = std::sqrt(12150.0 / 193947611.0);
    const auto fs_2083725_64 = std::sqrt(32558.203125);
    const auto fs_25725 = std::sqrt(25725.0);
    const auto fs_13861125_163592 = std::sqrt(13861125.0 / 163592.0);
    const auto fs_231525_121 = std::sqrt(231525.0 / 121.0);
    const auto fs_540225_2149004 = std::sqrt(540225.0 / 2149004.0);
    const auto fs_27722250_5909761 = std::sqrt(27722250.0 / 5909761.0);
    const auto fs_411600_20449 = std::sqrt(411600.0 / 20449.0);
    const auto fs_77616_1356277 = std::sqrt(77616.0 / 1356277.0);
    const auto fs_1225_537251 = std::sqrt(1225.0 / 537251.0);
    const auto fs_27722250_2133423721 = std::sqrt(27722250.0 / 2133423721.0);
    const auto fs_1372_61347 = std::sqrt(1372.0 / 61347.0);
    const auto fs_297675_512 = std::sqrt(581.396484375);
    const auto fs_13125_8 = std::sqrt(1640.625);
    const auto fs_3675_8 = std::sqrt(459.375);
    const auto fs_9646875_81796 = std::sqrt(9646875.0 / 81796.0);
    const auto fs_55125_2288 = std::sqrt(55125.0 / 2288.0);
    const auto fs_118125_338 = std::sqrt(118125.0 / 338.0);
    const auto fs_33075_968 = std::sqrt(33075.0 / 968.0);
    const auto fs_694575_1074502 = std::sqrt(694575.0 / 1074502.0);
    const auto fs_297675_82654 = std::sqrt(297675.0 / 82654.0);
    const auto fs_38587500_5909761 = std::sqrt(38587500.0 / 5909761.0);
    const auto fs_55125_41327 = std::sqrt(55125.0 / 41327.0);
    const auto fs_1050_169 = std::sqrt(1050.0 / 169.0);
    const auto fs_7350_20449 = std::sqrt(7350.0 / 20449.0);
    const auto fs_693_2712554 = std::sqrt(693.0 / 2712554.0);
    const auto fs_891_79781 = std::sqrt(891.0 / 79781.0);
    const auto fs_3150_537251 = std::sqrt(3150.0 / 537251.0);
    const auto fs_1350_41327 = std::sqrt(1350.0 / 41327.0);
    const auto fs_38587500_2133423721 = std::sqrt(38587500.0 / 2133423721.0);
    const auto fs_55125_14919047 = std::sqrt(55125.0 / 14919047.0);
    const auto fs_1400_146523 = std::sqrt(1400.0 / 146523.0);
    const auto fs_49_122694 = std::sqrt(49.0 / 122694.0);
    const auto fs_4465125_512 = std::sqrt(8720.947265625);
    const auto fs_7875_8 = std::sqrt(984.375);
    const auto fs_55125_8 = std::sqrt(6890.625);
    const auto fs_3781575_81796 = std::sqrt(3781575.0 / 81796.0);
    const auto fs_20475_176 = std::sqrt(20475.0 / 176.0);
    const auto fs_70875_338 = std::sqrt(70875.0 / 338.0);
    const auto fs_496125_968 = std::sqrt(496125.0 / 968.0);
    const auto fs_952560_537251 = std::sqrt(952560.0 / 537251.0);
    const auto fs_19845_41327 = std::sqrt(19845.0 / 41327.0);
    const auto fs_15126300_5909761 = std::sqrt(15126300.0 / 5909761.0);
    const auto fs_20475_3179 = std::sqrt(20475.0 / 3179.0);
    const auto fs_630_169 = std::sqrt(630.0 / 169.0);
    const auto fs_110250_20449 = std::sqrt(110250.0 / 20449.0);
    const auto fs_4455_2712554 = std::sqrt(4455.0 / 2712554.0);
    const auto fs_2475_79781 = std::sqrt(2475.0 / 79781.0);
    const auto fs_8640_537251 = std::sqrt(8640.0 / 537251.0);
    const auto fs_180_41327 = std::sqrt(180.0 / 41327.0);
    const auto fs_15126300_2133423721 = std::sqrt(15126300.0 / 2133423721.0);
    const auto fs_20475_1147619 = std::sqrt(20475.0 / 1147619.0);
    const auto fs_280_48841 = std::sqrt(280.0 / 48841.0);
    const auto fs_245_40898 = std::sqrt(245.0 / 40898.0);
    const auto fs_4862025_512 = std::sqrt(9496.142578125);
    const auto fs_15125_16 = std::sqrt(945.3125);
    const auto fs_60025_8 = std::sqrt(7503.125);
    const auto fs_54675_16 = std::sqrt(3417.1875);
    const auto fs_231525_81796 = std::sqrt(231525.0 / 81796.0);
    const auto fs_14175_29744 = std::sqrt(14175.0 / 29744.0);
    const auto fs_136125_676 = std::sqrt(136125.0 / 676.0);
    const auto fs_540225_968 = std::sqrt(540225.0 / 968.0);
    const auto fs_675_4 = std::sqrt(168.75);
    const auto fs_52397415_23639044 = std::sqrt(52397415.0 / 23639044.0);
    const auto fs_33075_82654 = std::sqrt(33075.0 / 82654.0);
    const auto fs_926100_5909761 = std::sqrt(926100.0 / 5909761.0);
    const auto fs_14175_537251 = std::sqrt(14175.0 / 537251.0);
    const auto fs_605_169 = std::sqrt(605.0 / 169.0);
    const auto fs_120050_20449 = std::sqrt(120050.0 / 20449.0);
    const auto fs_675_484 = std::sqrt(675.0 / 484.0);
    const auto fs_200475_35263202 = std::sqrt(200475.0 / 35263202.0);
    const auto fs_66825_1356277 = std::sqrt(66825.0 / 1356277.0);
    const auto fs_118815_5909761 = std::sqrt(118815.0 / 5909761.0);
    const auto fs_150_41327 = std::sqrt(150.0 / 41327.0);
    const auto fs_926100_2133423721 = std::sqrt(926100.0 / 2133423721.0);
    const auto fs_14175_193947611 = std::sqrt(14175.0 / 193947611.0);
    const auto fs_2420_439569 = std::sqrt(2420.0 / 439569.0);
    const auto fs_2401_368082 = std::sqrt(2401.0 / 368082.0);
    const auto fs_27_20449 = std::sqrt(27.0 / 20449.0);
    const auto fs_2679075_1024 = std::sqrt(2616.2841796875);
    const auto fs_99225_32 = std::sqrt(3100.78125);
    const auto fs_893025_8 = std::sqrt(111628.125);
    const auto fs_625_2 = std::sqrt(312.5);
    const auto fs_2450 = std::sqrt(2450.0);
    const auto fs_36450 = std::sqrt(36450.0);
    const auto fs_18533025_163592 = std::sqrt(18533025.0 / 163592.0);
    const auto fs_3444525_59488 = std::sqrt(3444525.0 / 59488.0);
    const auto fs_11250_169 = std::sqrt(11250.0 / 169.0);
    const auto fs_22050_121 = std::sqrt(22050.0 / 121.0);
    const auto fs_1800 = std::sqrt(1800.0);
    const auto fs_16074450_5909761 = std::sqrt(16074450.0 / 5909761.0);
    const auto fs_138915_82654 = std::sqrt(138915.0 / 82654.0);
    const auto fs_37066050_5909761 = std::sqrt(37066050.0 / 5909761.0);
    const auto fs_3444525_1074502 = std::sqrt(3444525.0 / 1074502.0);
    const auto fs_200_169 = std::sqrt(200.0 / 169.0);
    const auto fs_39200_20449 = std::sqrt(39200.0 / 20449.0);
    const auto fs_1800_121 = std::sqrt(1800.0 / 121.0);
    const auto fs_490050_17631601 = std::sqrt(490050.0 / 17631601.0);
    const auto fs_155925_2712554 = std::sqrt(155925.0 / 2712554.0);
    const auto fs_145800_5909761 = std::sqrt(145800.0 / 5909761.0);
    const auto fs_630_41327 = std::sqrt(630.0 / 41327.0);
    const auto fs_37066050_2133423721 = std::sqrt(37066050.0 / 2133423721.0);
    const auto fs_3444525_387895222 = std::sqrt(3444525.0 / 387895222.0);
    const auto fs_800_439569 = std::sqrt(800.0 / 439569.0);
    const auto fs_392_184041 = std::sqrt(392.0 / 184041.0);
    const auto fs_288_20449 = std::sqrt(288.0 / 20449.0);
    const auto fs_893025_32 = std::sqrt(27907.03125);
    const auto fs_2083725_256 = std::sqrt(8139.55078125);
    const auto fs_3472875_256 = std::sqrt(13565.91796875);
    const auto fs_6251175_128 = std::sqrt(48837.3046875);
    const auto fs_2625_2 = std::sqrt(1312.5);
    const auto fs_25725_4 = std::sqrt(6431.25);
    const auto fs_42875_4 = std::sqrt(10718.75);
    const auto fs_127575_8 = std::sqrt(15946.875);
    const auto fs_30305025_654368 = std::sqrt(30305025.0 / 654368.0);
    const auto fs_28923075_654368 = std::sqrt(28923075.0 / 654368.0);
    const auto fs_47250_169 = std::sqrt(47250.0 / 169.0);
    const auto fs_231525_484 = std::sqrt(231525.0 / 484.0);
    const auto fs_385875_484 = std::sqrt(385875.0 / 484.0);
    const auto fs_1575_2 = std::sqrt(787.5);
    const auto fs_1111320_5909761 = std::sqrt(1111320.0 / 5909761.0);
    const auto fs_2917215_2149004 = std::sqrt(2917215.0 / 2149004.0);
    const auto fs_30305025_11819522 = std::sqrt(30305025.0 / 11819522.0);
    const auto fs_28923075_11819522 = std::sqrt(28923075.0 / 11819522.0);
    const auto fs_840_169 = std::sqrt(840.0 / 169.0);
    const auto fs_102900_20449 = std::sqrt(102900.0 / 20449.0);
    const auto fs_171500_20449 = std::sqrt(171500.0 / 20449.0);
    const auto fs_1575_242 = std::sqrt(1575.0 / 242.0);
    const auto fs_467775_17631601 = std::sqrt(467775.0 / 17631601.0);
    const auto fs_72765_1356277 = std::sqrt(72765.0 / 1356277.0);
    const auto fs_10080_5909761 = std::sqrt(10080.0 / 5909761.0);
    const auto fs_6615_537251 = std::sqrt(6615.0 / 537251.0);
    const auto fs_30305025_4266847442 = std::sqrt(30305025.0 / 4266847442.0);
    const auto fs_28923075_4266847442 = std::sqrt(28923075.0 / 4266847442.0);
    const auto fs_1120_146523 = std::sqrt(1120.0 / 146523.0);
    const auto fs_343_61347 = std::sqrt(343.0 / 61347.0);
    const auto fs_1715_184041 = std::sqrt(1715.0 / 184041.0);
    const auto fs_126_20449 = std::sqrt(126.0 / 20449.0);
    const auto fs_6251175_512 = std::sqrt(12209.326171875);
    const auto fs_694575_128 = std::sqrt(5426.3671875);
    const auto fs_8575_2 = std::sqrt(4287.5);
    const auto fs_91125_327184 = std::sqrt(91125.0 / 327184.0);
    const auto fs_77175_242 = std::sqrt(77175.0 / 242.0);
    const auto fs_231525_537251 = std::sqrt(231525.0 / 537251.0);
    const auto fs_91125_5909761 = std::sqrt(91125.0 / 5909761.0);
    const auto fs_68600_20449 = std::sqrt(68600.0 / 20449.0);
    const auto fs_112266_1356277 = std::sqrt(112266.0 / 1356277.0);
    const auto fs_2100_537251 = std::sqrt(2100.0 / 537251.0);
    const auto fs_91125_2133423721 = std::sqrt(91125.0 / 2133423721.0);
    const auto fs_686_184041 = std::sqrt(686.0 / 184041.0);
    const auto fs_23625_32 = std::sqrt(738.28125);
    const auto fs_1929375_14872 = std::sqrt(1929375.0 / 14872.0);
    const auto fs_39375_572 = std::sqrt(39375.0 / 572.0);
    const auto fs_212625_1352 = std::sqrt(212625.0 / 1352.0);
    const auto fs_416745_330616 = std::sqrt(416745.0 / 330616.0);
    const auto fs_496125_165308 = std::sqrt(496125.0 / 165308.0);
    const auto fs_3858750_537251 = std::sqrt(3858750.0 / 537251.0);
    const auto fs_157500_41327 = std::sqrt(157500.0 / 41327.0);
    const auto fs_945_338 = std::sqrt(945.0 / 338.0);
    const auto fs_2079_2712554 = std::sqrt(2079.0 / 2712554.0);
    const auto fs_396_79781 = std::sqrt(396.0 / 79781.0);
    const auto fs_945_82654 = std::sqrt(945.0 / 82654.0);
    const auto fs_1125_41327 = std::sqrt(1125.0 / 41327.0);
    const auto fs_3858750_193947611 = std::sqrt(3858750.0 / 193947611.0);
    const auto fs_157500_14919047 = std::sqrt(157500.0 / 14919047.0);
    const auto fs_210_48841 = std::sqrt(210.0 / 48841.0);
    const auto fs_297675_128 = std::sqrt(2325.5859375);
    const auto fs_63525_32 = std::sqrt(1985.15625);
    const auto fs_3675_2 = std::sqrt(1837.5);
    const auto fs_385875_81796 = std::sqrt(385875.0 / 81796.0);
    const auto fs_126000_1859 = std::sqrt(126000.0 / 1859.0);
    const auto fs_571725_1352 = std::sqrt(571725.0 / 1352.0);
    const auto fs_33075_242 = std::sqrt(33075.0 / 242.0);
    const auto fs_10029663_4298008 = std::sqrt(10029663.0 / 4298008.0);
    const auto fs_535815_330616 = std::sqrt(535815.0 / 330616.0);
    const auto fs_1543500_5909761 = std::sqrt(1543500.0 / 5909761.0);
    const auto fs_2016000_537251 = std::sqrt(2016000.0 / 537251.0);
    const auto fs_2541_338 = std::sqrt(2541.0 / 338.0);
    const auto fs_29400_20449 = std::sqrt(29400.0 / 20449.0);
    const auto fs_5544_1356277 = std::sqrt(5544.0 / 1356277.0);
    const auto fs_23760_1356277 = std::sqrt(23760.0 / 1356277.0);
    const auto fs_22743_1074502 = std::sqrt(22743.0 / 1074502.0);
    const auto fs_1215_82654 = std::sqrt(1215.0 / 82654.0);
    const auto fs_1543500_2133423721 = std::sqrt(1543500.0 / 2133423721.0);
    const auto fs_2016000_193947611 = std::sqrt(2016000.0 / 193947611.0);
    const auto fs_1694_146523 = std::sqrt(1694.0 / 146523.0);
    const auto fs_98_61347 = std::sqrt(98.0 / 61347.0);
    const auto fs_99225_8 = std::sqrt(12403.125);
    const auto fs_175_8 = std::sqrt(21.875);
    const auto fs_9800 = std::sqrt(9800.0);
    const auto fs_3472875_81796 = std::sqrt(3472875.0 / 81796.0);
    const auto fs_637875_14872 = std::sqrt(637875.0 / 14872.0);
    const auto fs_1575_338 = std::sqrt(1575.0 / 338.0);
    const auto fs_88200_121 = std::sqrt(88200.0 / 121.0);
    const auto fs_3716307_2149004 = std::sqrt(3716307.0 / 2149004.0);
    const auto fs_9261_330616 = std::sqrt(9261.0 / 330616.0);
    const auto fs_13891500_5909761 = std::sqrt(13891500.0 / 5909761.0);
    const auto fs_1275750_537251 = std::sqrt(1275750.0 / 537251.0);
    const auto fs_14_169 = std::sqrt(14.0 / 169.0);
    const auto fs_156800_20449 = std::sqrt(156800.0 / 20449.0);
    const auto fs_16038_1356277 = std::sqrt(16038.0 / 1356277.0);
    const auto fs_93555_2712554 = std::sqrt(93555.0 / 2712554.0);
    const auto fs_8427_537251 = std::sqrt(8427.0 / 537251.0);
    const auto fs_21_82654 = std::sqrt(21.0 / 82654.0);
    const auto fs_13891500_2133423721 = std::sqrt(13891500.0 / 2133423721.0);
    const auto fs_1275750_193947611 = std::sqrt(1275750.0 / 193947611.0);
    const auto fs_56_439569 = std::sqrt(56.0 / 439569.0);
    const auto fs_1568_184041 = std::sqrt(1568.0 / 184041.0);
    const auto fs_1488375_512 = std::sqrt(2906.982421875);
    const auto fs_1200 = std::sqrt(1200.0);
    const auto fs_18375_8 = std::sqrt(2296.875);
    const auto fs_91125_16 = std::sqrt(5695.3125);
    const auto fs_18907875_327184 = std::sqrt(18907875.0 / 327184.0);
    const auto fs_1913625_327184 = std::sqrt(1913625.0 / 327184.0);
    const auto fs_43200_169 = std::sqrt(43200.0 / 169.0);
    const auto fs_165375_968 = std::sqrt(165375.0 / 968.0);
    const auto fs_1125_4 = std::sqrt(281.25);
    const auto f_1449_2431 = 1449.0 / 2431.0;
    const auto fs_750141_1074502 = std::sqrt(750141.0 / 1074502.0);
    const auto fs_18907875_5909761 = std::sqrt(18907875.0 / 5909761.0);
    const auto fs_1913625_5909761 = std::sqrt(1913625.0 / 5909761.0);
    const auto fs_768_169 = std::sqrt(768.0 / 169.0);
    const auto fs_36750_20449 = std::sqrt(36750.0 / 20449.0);
    const auto fs_1125_484 = std::sqrt(1125.0 / 484.0);
    const auto fs_427680_17631601 = std::sqrt(427680.0 / 17631601.0);
    const auto fs_66528_1356277 = std::sqrt(66528.0 / 1356277.0);
    const auto f_138_2431 = 138.0 / 2431.0;
    const auto fs_3402_537251 = std::sqrt(3402.0 / 537251.0);
    const auto fs_18907875_2133423721 = std::sqrt(18907875.0 / 2133423721.0);
    const auto fs_1913625_2133423721 = std::sqrt(1913625.0 / 2133423721.0);
    const auto fs_1024_146523 = std::sqrt(1024.0 / 146523.0);
    const auto fs_245_122694 = std::sqrt(245.0 / 122694.0);
    const auto fs_45_20449 = std::sqrt(45.0 / 20449.0);
    const auto fs_4465125_1024 = std::sqrt(4360.4736328125);
    const auto fs_31255875_256 = std::sqrt(122093.26171875);
    const auto fs_875 = std::sqrt(875.0);
    const auto fs_637875_16 = std::sqrt(39867.1875);
    const auto fs_126000_20449 = std::sqrt(126000.0 / 20449.0);
    const auto fs_2460375_40898 = std::sqrt(2460375.0 / 40898.0);
    const auto fs_31500_169 = std::sqrt(31500.0 / 169.0);
    const auto fs_7875_4 = std::sqrt(1968.75);
    const auto fs_1250235_5909761 = std::sqrt(1250235.0 / 5909761.0);
    const auto fs_889056_537251 = std::sqrt(889056.0 / 537251.0);
    const auto fs_2016000_5909761 = std::sqrt(2016000.0 / 5909761.0);
    const auto fs_19683000_5909761 = std::sqrt(19683000.0 / 5909761.0);
    const auto fs_560_169 = std::sqrt(560.0 / 169.0);
    const auto fs_7875_484 = std::sqrt(7875.0 / 484.0);
    const auto fs_1372140_17631601 = std::sqrt(1372140.0 / 17631601.0);
    const auto fs_74844_1356277 = std::sqrt(74844.0 / 1356277.0);
    const auto fs_11340_5909761 = std::sqrt(11340.0 / 5909761.0);
    const auto fs_8064_537251 = std::sqrt(8064.0 / 537251.0);
    const auto fs_2016000_2133423721 = std::sqrt(2016000.0 / 2133423721.0);
    const auto fs_19683000_2133423721 = std::sqrt(19683000.0 / 2133423721.0);
    const auto fs_2240_439569 = std::sqrt(2240.0 / 439569.0);
    const auto fs_315_20449 = std::sqrt(315.0 / 20449.0);
    const auto fs_31255875_1024 = std::sqrt(30523.3154296875);
    const auto fs_18753525_256 = std::sqrt(73255.95703125);
    const auto fs_382725_16 = std::sqrt(23920.3125);
    const auto fs_4876875_81796 = std::sqrt(4876875.0 / 81796.0);
    const auto fs_4725_4 = std::sqrt(1181.25);
    const auto fs_46305_20449 = std::sqrt(46305.0 / 20449.0);
    const auto fs_67500_20449 = std::sqrt(67500.0 / 20449.0);
    const auto fs_4725_484 = std::sqrt(4725.0 / 484.0);
    const auto fs_1796256_17631601 = std::sqrt(1796256.0 / 17631601.0);
    const auto fs_420_20449 = std::sqrt(420.0 / 20449.0);
    const auto fs_67500_7382089 = std::sqrt(67500.0 / 7382089.0);
    const auto fs_189_20449 = std::sqrt(189.0 / 20449.0);
    const auto fs_18753525_1024 = std::sqrt(18313.9892578125);
    const auto f_75_4 = 18.75;
    const auto fs_826875_3718 = std::sqrt(826875.0 / 3718.0);
    const auto f_225_26 = 225.0 / 26.0;
    const auto fs_694575_165308 = std::sqrt(694575.0 / 165308.0);
    const auto fs_6615000_537251 = std::sqrt(6615000.0 / 537251.0);
    const auto f_15_13 = 15.0 / 13.0;
    const auto fs_1575_41327 = std::sqrt(1575.0 / 41327.0);
    const auto fs_6615000_193947611 = std::sqrt(6615000.0 / 193947611.0);
    const auto f_10_221 = 10.0 / 221.0;
    const auto f_60 = 60.0;
    const auto fs_165375_7436 = std::sqrt(165375.0 / 7436.0);
    const auto f_360_13 = 360.0 / 13.0;
    const auto fs_194481_41327 = std::sqrt(194481.0 / 41327.0);
    const auto fs_661500_537251 = std::sqrt(661500.0 / 537251.0);
    const auto f_48_13 = 48.0 / 13.0;
    const auto fs_24255_1356277 = std::sqrt(24255.0 / 1356277.0);
    const auto fs_1764_41327 = std::sqrt(1764.0 / 41327.0);
    const auto fs_661500_193947611 = std::sqrt(661500.0 / 193947611.0);
    const auto f_32_221 = 32.0 / 221.0;
    const auto fs_694575_64 = std::sqrt(10852.734375);
    const auto f_145_4 = 36.25;
    const auto fs_8575 = std::sqrt(8575.0);
    const auto fs_2976750_20449 = std::sqrt(2976750.0 / 20449.0);
    const auto f_435_26 = 435.0 / 26.0;
    const auto fs_77175_121 = std::sqrt(77175.0 / 121.0);
    const auto fs_3176523_2149004 = std::sqrt(3176523.0 / 2149004.0);
    const auto fs_47628000_5909761 = std::sqrt(47628000.0 / 5909761.0);
    const auto f_29_13 = 29.0 / 13.0;
    const auto fs_137200_20449 = std::sqrt(137200.0 / 20449.0);
    const auto fs_58212_1356277 = std::sqrt(58212.0 / 1356277.0);
    const auto fs_7203_537251 = std::sqrt(7203.0 / 537251.0);
    const auto fs_47628000_2133423721 = std::sqrt(47628000.0 / 2133423721.0);
    const auto f_58_663 = 58.0 / 663.0;
    const auto fs_1372_184041 = std::sqrt(1372.0 / 184041.0);
    const auto fs_6251175_256 = std::sqrt(24418.65234375);
    const auto f_45 = 45.0;
    const auto fs_77175_4 = std::sqrt(19293.75);
    const auto fs_4465125_163592 = std::sqrt(4465125.0 / 163592.0);
    const auto f_270_13 = 270.0 / 13.0;
    const auto fs_694575_484 = std::sqrt(694575.0 / 484.0);
    const auto fs_18522_537251 = std::sqrt(18522.0 / 537251.0);
    const auto fs_8930250_5909761 = std::sqrt(8930250.0 / 5909761.0);
    const auto f_36_13 = 36.0 / 13.0;
    const auto fs_308700_20449 = std::sqrt(308700.0 / 20449.0);
    const auto fs_99792_1356277 = std::sqrt(99792.0 / 1356277.0);
    const auto fs_168_537251 = std::sqrt(168.0 / 537251.0);
    const auto fs_8930250_2133423721 = std::sqrt(8930250.0 / 2133423721.0);
    const auto f_24_221 = 24.0 / 221.0;
    const auto fs_343_20449 = std::sqrt(343.0 / 20449.0);
    const auto fs_13395375_256 = std::sqrt(52325.68359375);
    const auto f_15 = 15.0;
    const auto fs_273375_16 = std::sqrt(17085.9375);
    const auto fs_23625_676 = std::sqrt(23625.0 / 676.0);
    const auto f_90_13 = 90.0 / 13.0;
    const auto fs_3375_4 = std::sqrt(843.75);
    const auto fs_64827_34969 = std::sqrt(64827.0 / 34969.0);
    const auto fs_94500_48841 = std::sqrt(94500.0 / 48841.0);
    const auto f_12_13 = 12.0 / 13.0;
    const auto fs_3375_484 = std::sqrt(3375.0 / 484.0);
    const auto fs_1746360_17631601 = std::sqrt(1746360.0 / 17631601.0);
    const auto fs_588_34969 = std::sqrt(588.0 / 34969.0);
    const auto fs_94500_17631601 = std::sqrt(94500.0 / 17631601.0);
    const auto f_8_221 = 8.0 / 221.0;
    const auto fs_135_20449 = std::sqrt(135.0 / 20449.0);
    const auto fs_13395375_1024 = std::sqrt(13081.4208984375);
    const auto f_2205_16 = 137.8125;
    const auto f_2835_8 = 354.375;
    const auto f_50 = 50.0;
    const auto f_245_2 = 122.5;
    const auto f_405_2 = 202.5;
    const auto f_1575_143 = 1575.0 / 143.0;
    const auto f_300_13 = 300.0 / 13.0;
    const auto f_735_22 = 735.0 / 22.0;
    const auto f_4410_2431 = 4410.0 / 2431.0;
    const auto f_6300_2431 = 6300.0 / 2431.0;
    const auto f_40_13 = 40.0 / 13.0;
    const auto f_490_143 = 490.0 / 143.0;
    const auto f_45_11 = 45.0 / 11.0;
    const auto f_1386_4199 = 1386.0 / 4199.0;
    const auto f_420_2431 = 420.0 / 2431.0;
    const auto f_6300_46189 = 6300.0 / 46189.0;
    const auto f_80_663 = 80.0 / 663.0;
    const auto f_49_429 = 49.0 / 429.0;
    const auto f_18_143 = 18.0 / 143.0;
    const auto f_2835_16 = 177.1875;

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph11_p1, ph11_p11, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p11 = ph11_p11[k];

        pc_0[k] = e_0 * (-fs_9823275_512 * h3_p1 + fs_29469825_256 * r_2 * h1_p1) + e_1 * (-fs_12375_16 * h5_p1 + fs_121275_8 * r_2 * h3_p1 - fs_601425_16 * r_4 * h1_p1) + e_2 * (-fs_118125_29744 * h7_p1 + fs_111375_676 * r_2 * h5_p1 - fs_99225_88 * r_4 * h3_p1 + fs_7425_4 * r_6 * h1_p1) + e_3 * (-fs_6615_2149004 * h9_p1 + fs_118125_537251 * r_2 * h7_p1 - fs_495_169 * r_4 * h5_p1 + fs_22050_1859 * r_6 * h3_p1 - fs_675_44 * r_8 * h1_p1) + e_4 * (-fs_9_35263202 * h11_p1 - fs_693_4199 * h11_p11 + fs_15_537251 * r_2 * h9_p1 - fs_118125_193947611 * r_4 * h7_p1 + fs_220_48841 * r_6 * h5_p1 - fs_49_3718 * r_8 * h3_p1 + fs_27_1859 * r_10 * h1_p1) - fs_29469825_1024 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ph11_p10, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p10 = ph11_p10[k];

        pc_1[k] = fs_9823275_512 * e_0 * h3_p2 + e_1 * (fs_17325_8 * h5_p2 - fs_121275_8 * r_2 * h3_p2) + e_2 * (fs_637875_29744 * h7_p2 - fs_155925_338 * r_2 * h5_p2 + fs_99225_88 * r_4 * h3_p2) + e_3 * (fs_1323_48841 * h9_p2 - fs_637875_537251 * r_2 * h7_p2 + fs_1386_169 * r_4 * h5_p2 - fs_22050_1859 * r_6 * h3_p2) + e_4 * (fs_9_2712554 * h11_p2 - fs_315_4199 * h11_p10 - fs_12_48841 * r_2 * h9_p2 + fs_637875_193947611 * r_4 * h7_p2 - fs_616_48841 * r_6 * h5_p2 + fs_49_3718 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph9_p9, ph11_p3, ph11_p9, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p9 = ph11_p9[k];

        pc_2[k] = -fs_3274425_512 * e_0 * h3_p3 + e_1 * (-fs_5775_2 * h5_p3 + fs_40425_8 * r_2 * h3_p3) + e_2 * (-fs_1771875_29744 * h7_p3 + fs_103950_169 * r_2 * h5_p3 - fs_33075_88 * r_4 * h3_p3) + e_3 * (-fs_6174_48841 * h9_p3 - fs_1323_442 * h9_p9 + fs_1771875_537251 * r_2 * h7_p3 - fs_1848_169 * r_4 * h5_p3 + fs_7350_1859 * r_6 * h3_p3) + e_4 * (-fs_63_2712554 * h11_p3 - fs_135_4199 * h11_p9 + fs_56_48841 * r_2 * h9_p3 + fs_6_221 * r_2 * h9_p9 - fs_1771875_193947611 * r_4 * h7_p3 + fs_2464_146523 * r_6 * h5_p3 - fs_49_11154 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p4, ph5_p5, ph7_p4, ph7_p5, ph7_p7, ph9_p4, ph9_p5, ph9_p7, ph9_p8, ph11_p4, ph11_p5, ph11_p7, ph11_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p5 = ph11_p5[k];
        const auto h11_p7 = ph11_p7[k];
        const auto h11_p8 = ph11_p8[k];

        pc_3[k] = fs_17325_8 * e_1 * h5_p4 + e_2 * (fs_590625_5408 * h7_p4 - fs_155925_338 * r_2 * h5_p4) + e_3 * (fs_3087_7514 * h9_p4 - fs_882_221 * h9_p8 - fs_590625_97682 * r_2 * h7_p4 + fs_1386_169 * r_4 * h5_p4) + e_4 * (fs_315_2712554 * h11_p4 - fs_54_4199 * h11_p8 - fs_14_3757 * r_2 * h9_p4 + fs_8_221 * r_2 * h9_p8 + fs_590625_35263202 * r_4 * h7_p4 - fs_616_48841 * r_6 * h5_p4);

        pc_4[k] = -fs_12375_16 * e_1 * h5_p5 + e_2 * (-fs_759375_5408 * h7_p5 - fs_23625_416 * h7_p7 + fs_111375_676 * r_2 * h5_p5) + e_3 * (-fs_15435_15028 * h9_p5 - fs_12348_3757 * h9_p7 + fs_759375_97682 * r_2 * h7_p5 + fs_23625_7514 * r_2 * h7_p7 - fs_495_169 * r_4 * h5_p5) + e_4 * (-fs_630_1356277 * h11_p5 - fs_378_79781 * h11_p7 + fs_35_3757 * r_2 * h9_p5 + fs_112_3757 * r_2 * h9_p7 - fs_759375_35263202 * r_4 * h7_p5 - fs_23625_2712554 * r_4 * h7_p7 + fs_220_48841 * r_6 * h5_p5);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m5, ph7_m7, ph7_m6, ph7_m5, ph9_m7, ph9_m6, ph9_m5, ph11_m7, ph11_m6, ph11_m5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m5 = ph5_m5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m5 = ph11_m5[k];

        pc_5[k] = fs_50625_208 * e_2 * h7_m6 + e_3 * (fs_15435_3757 * h9_m6 - fs_50625_3757 * r_2 * h7_m6) + e_4 * (fs_252_79781 * h11_m6 - fs_140_3757 * r_2 * h9_m6 + fs_50625_1356277 * r_4 * h7_m6);

        pc_6[k] = -fs_12375_16 * e_1 * h5_m5 + e_2 * (fs_23625_416 * h7_m7 - fs_759375_5408 * h7_m5 + fs_111375_676 * r_2 * h5_m5) + e_3 * (fs_12348_3757 * h9_m7 - fs_15435_15028 * h9_m5 - fs_23625_7514 * r_2 * h7_m7 + fs_759375_97682 * r_2 * h7_m5 - fs_495_169 * r_4 * h5_m5) + e_4 * (fs_378_79781 * h11_m7 - fs_630_1356277 * h11_m5 - fs_112_3757 * r_2 * h9_m7 + fs_35_3757 * r_2 * h9_m5 + fs_23625_2712554 * r_4 * h7_m7 - fs_759375_35263202 * r_4 * h7_m5 + fs_220_48841 * r_6 * h5_m5);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph7_m4, ph9_m8, ph9_m4, ph11_m8, ph11_m4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m4 = ph11_m4[k];

        pc_7[k] = fs_17325_8 * e_1 * h5_m4 + e_2 * (fs_590625_5408 * h7_m4 - fs_155925_338 * r_2 * h5_m4) + e_3 * (fs_882_221 * h9_m8 + fs_3087_7514 * h9_m4 - fs_590625_97682 * r_2 * h7_m4 + fs_1386_169 * r_4 * h5_m4) + e_4 * (fs_54_4199 * h11_m8 + fs_315_2712554 * h11_m4 - fs_8_221 * r_2 * h9_m8 - fs_14_3757 * r_2 * h9_m4 + fs_590625_35263202 * r_4 * h7_m4 - fs_616_48841 * r_6 * h5_m4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m3, ph7_m3, ph9_m9, ph9_m3, ph11_m9, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m3 = ph11_m3[k];

        pc_8[k] = -fs_3274425_512 * e_0 * h3_m3 + e_1 * (-fs_5775_2 * h5_m3 + fs_40425_8 * r_2 * h3_m3) + e_2 * (-fs_1771875_29744 * h7_m3 + fs_103950_169 * r_2 * h5_m3 - fs_33075_88 * r_4 * h3_m3) + e_3 * (fs_1323_442 * h9_m9 - fs_6174_48841 * h9_m3 + fs_1771875_537251 * r_2 * h7_m3 - fs_1848_169 * r_4 * h5_m3 + fs_7350_1859 * r_6 * h3_m3) + e_4 * (fs_135_4199 * h11_m9 - fs_63_2712554 * h11_m3 - fs_6_221 * r_2 * h9_m9 + fs_56_48841 * r_2 * h9_m3 - fs_1771875_193947611 * r_4 * h7_m3 + fs_2464_146523 * r_6 * h5_m3 - fs_49_11154 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m10, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m2 = ph11_m2[k];

        pc_9[k] = fs_9823275_512 * e_0 * h3_m2 + e_1 * (fs_17325_8 * h5_m2 - fs_121275_8 * r_2 * h3_m2) + e_2 * (fs_637875_29744 * h7_m2 - fs_155925_338 * r_2 * h5_m2 + fs_99225_88 * r_4 * h3_m2) + e_3 * (fs_1323_48841 * h9_m2 - fs_637875_537251 * r_2 * h7_m2 + fs_1386_169 * r_4 * h5_m2 - fs_22050_1859 * r_6 * h3_m2) + e_4 * (fs_315_4199 * h11_m10 + fs_9_2712554 * h11_m2 - fs_12_48841 * r_2 * h9_m2 + fs_637875_193947611 * r_4 * h7_m2 - fs_616_48841 * r_6 * h5_m2 + fs_49_3718 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m1, ph11_m11, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m11 = ph11_m11[k];
        const auto h11_m1 = ph11_m1[k];

        pc_10[k] = e_0 * (-fs_9823275_512 * h3_m1 + fs_29469825_256 * r_2 * h1_m1) + e_1 * (-fs_12375_16 * h5_m1 + fs_121275_8 * r_2 * h3_m1 - fs_601425_16 * r_4 * h1_m1) + e_2 * (-fs_118125_29744 * h7_m1 + fs_111375_676 * r_2 * h5_m1 - fs_99225_88 * r_4 * h3_m1 + fs_7425_4 * r_6 * h1_m1) + e_3 * (-fs_6615_2149004 * h9_m1 + fs_118125_537251 * r_2 * h7_m1 - fs_495_169 * r_4 * h5_m1 + fs_22050_1859 * r_6 * h3_m1 - fs_675_44 * r_8 * h1_m1) + e_4 * (fs_693_4199 * h11_m11 - fs_9_35263202 * h11_m1 + fs_15_537251 * r_2 * h9_m1 - fs_118125_193947611 * r_4 * h7_m1 + fs_220_48841 * r_6 * h5_m1 - fs_49_3718 * r_8 * h3_m1 + fs_27_1859 * r_10 * h1_m1) - fs_29469825_1024 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph11_0, ph11_p10, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p10 = ph11_p10[k];

        pc_11[k] = e_0 * (-fs_9823275_256 * h3_0 + fs_9823275_256 * r_2 * h1_0) + e_1 * (-fs_61875_16 * h5_0 + fs_121275_4 * r_2 * h3_0 - fs_200475_16 * r_4 * h1_0) + e_2 * (-fs_275625_7436 * h7_0 + fs_556875_676 * r_2 * h5_0 - fs_99225_44 * r_4 * h3_0 + fs_2475_4 * r_6 * h1_0) + e_3 * (-fs_99225_2149004 * h9_0 + fs_1102500_537251 * r_2 * h7_0 - fs_2475_169 * r_4 * h5_0 + fs_44100_1859 * r_6 * h3_0 - fs_225_44 * r_8 * h1_0) + e_4 * (-fs_99_17631601 * h11_0 - fs_378_4199 * h11_p10 + fs_225_537251 * r_2 * h9_0 - fs_1102500_193947611 * r_4 * h7_0 + fs_1100_48841 * r_6 * h5_0 - fs_49_1859 * r_8 * h3_0 + fs_9_1859 * r_10 * h1_0) - fs_9823275_1024 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_12[k] = fs_49116375_512 * e_0 * r_2 * h1_p1 + e_1 * (fs_66825_32 * h5_p1 - fs_1002375_32 * r_4 * h1_p1) + e_2 * (fs_86625_1352 * h7_p1 - fs_601425_1352 * r_2 * h5_p1 + fs_12375_8 * r_6 * h1_p1) + e_3 * (fs_3969_25432 * h9_p1 + fs_3969_884 * h9_p9 - fs_173250_48841 * r_2 * h7_p1 + fs_2673_338 * r_4 * h5_p1 - fs_1125_88 * r_8 * h1_p1) + e_4 * (fs_540_17631601 * h11_p1 - fs_360_4199 * h11_p9 - fs_9_6358 * r_2 * h9_p1 - fs_9_221 * r_2 * h9_p9 + fs_173250_17631601 * r_4 * h7_p1 - fs_594_48841 * r_6 * h5_p1 + fs_45_3718 * r_10 * h1_p1) - fs_49116375_2048 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_13[k] = fs_3274425_256 * e_0 * h3_p2 + e_1 * (-fs_5775_16 * h5_p2 - fs_40425_4 * r_2 * h3_p2) + e_2 * (-fs_189000_1859 * h7_p2 + fs_51975_676 * r_2 * h5_p2 + fs_33075_44 * r_4 * h3_p2) + e_3 * (-fs_53361_97682 * h9_p2 + fs_441_884 * h9_p8 + fs_3024000_537251 * r_2 * h7_p2 - fs_231_169 * r_4 * h5_p2 - fs_14700_1859 * r_6 * h3_p2) + e_4 * (-fs_243_1356277 * h11_p2 - fs_243_4199 * h11_p8 + fs_242_48841 * r_2 * h9_p2 - fs_1_221 * r_2 * h9_p8 - fs_3024000_193947611 * r_4 * h7_p2 + fs_308_146523 * r_6 * h5_p2 + fs_49_5577 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_14[k] = -fs_3274425_256 * e_0 * h3_p3 + e_1 * (-fs_5775_16 * h5_p3 + fs_40425_4 * r_2 * h3_p3) + e_2 * (fs_4921875_59488 * h7_p3 + fs_55125_416 * h7_p7 + fs_51975_676 * r_2 * h5_p3 - fs_33075_44 * r_4 * h3_p3) + e_3 * (fs_250047_195364 * h9_p3 - fs_1323_3757 * h9_p7 - fs_4921875_1074502 * r_2 * h7_p3 - fs_55125_7514 * r_2 * h7_p7 - fs_231_169 * r_4 * h5_p3 + fs_14700_1859 * r_6 * h3_p3) + e_4 * (fs_1008_1356277 * h11_p3 - fs_2592_79781 * h11_p7 - fs_567_48841 * r_2 * h9_p3 + fs_12_3757 * r_2 * h9_p7 + fs_4921875_387895222 * r_4 * h7_p3 + fs_55125_2712554 * r_4 * h7_p7 + fs_308_146523 * r_6 * h5_p3 - fs_49_5577 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m5, ph5_p4, ph7_m5, ph7_p4, ph7_p6, ph9_m5, ph9_p4, ph9_p6, ph11_m5, ph11_p4, ph11_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_15[k] = fs_66825_32 * e_1 * h5_p4 + e_2 * (-fs_28125_1352 * h7_p4 + fs_1125_13 * h7_p6 - fs_601425_1352 * r_2 * h5_p4) + e_3 * (-fs_64827_30056 * h9_p4 - fs_27783_15028 * h9_p6 + fs_56250_48841 * r_2 * h7_p4 - fs_18000_3757 * r_2 * h7_p6 + fs_2673_338 * r_4 * h5_p4) + e_4 * (-fs_6615_2712554 * h11_p4 - fs_1260_79781 * h11_p6 + fs_147_7514 * r_2 * h9_p4 + fs_63_3757 * r_2 * h9_p6 - fs_56250_17631601 * r_4 * h7_p4 + fs_18000_1356277 * r_4 * h7_p6 - fs_594_48841 * r_6 * h5_p4);

        pc_16[k] = -fs_61875_16 * e_1 * h5_m5 + e_2 * (-fs_16875_1352 * h7_m5 + fs_556875_676 * r_2 * h5_m5) + e_3 * (fs_77175_15028 * h9_m5 + fs_33750_48841 * r_2 * h7_m5 - fs_2475_169 * r_4 * h5_m5) + e_4 * (fs_18144_1356277 * h11_m5 - fs_175_3757 * r_2 * h9_m5 - fs_33750_17631601 * r_4 * h7_m5 + fs_1100_48841 * r_6 * h5_m5);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph7_m6, ph7_m4, ph9_m6, ph9_m4, ph11_m6, ph11_m4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m4 = ph11_m4[k];

        pc_17[k] = fs_66825_32 * e_1 * h5_m4 + e_2 * (-fs_1125_13 * h7_m6 - fs_28125_1352 * h7_m4 - fs_601425_1352 * r_2 * h5_m4) + e_3 * (fs_27783_15028 * h9_m6 - fs_64827_30056 * h9_m4 + fs_18000_3757 * r_2 * h7_m6 + fs_56250_48841 * r_2 * h7_m4 + fs_2673_338 * r_4 * h5_m4) + e_4 * (fs_1260_79781 * h11_m6 - fs_6615_2712554 * h11_m4 - fs_63_3757 * r_2 * h9_m6 + fs_147_7514 * r_2 * h9_m4 - fs_18000_1356277 * r_4 * h7_m6 - fs_56250_17631601 * r_4 * h7_m4 - fs_594_48841 * r_6 * h5_m4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_18[k] = -fs_3274425_256 * e_0 * h3_m3 + e_1 * (-fs_5775_16 * h5_m3 + fs_40425_4 * r_2 * h3_m3) + e_2 * (-fs_55125_416 * h7_m7 + fs_4921875_59488 * h7_m3 + fs_51975_676 * r_2 * h5_m3 - fs_33075_44 * r_4 * h3_m3) + e_3 * (fs_1323_3757 * h9_m7 + fs_250047_195364 * h9_m3 + fs_55125_7514 * r_2 * h7_m7 - fs_4921875_1074502 * r_2 * h7_m3 - fs_231_169 * r_4 * h5_m3 + fs_14700_1859 * r_6 * h3_m3) + e_4 * (fs_2592_79781 * h11_m7 + fs_1008_1356277 * h11_m3 - fs_12_3757 * r_2 * h9_m7 - fs_567_48841 * r_2 * h9_m3 - fs_55125_2712554 * r_4 * h7_m7 + fs_4921875_387895222 * r_4 * h7_m3 + fs_308_146523 * r_6 * h5_m3 - fs_49_5577 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_19[k] = fs_3274425_256 * e_0 * h3_m2 + e_1 * (-fs_5775_16 * h5_m2 - fs_40425_4 * r_2 * h3_m2) + e_2 * (-fs_189000_1859 * h7_m2 + fs_51975_676 * r_2 * h5_m2 + fs_33075_44 * r_4 * h3_m2) + e_3 * (-fs_441_884 * h9_m8 - fs_53361_97682 * h9_m2 + fs_3024000_537251 * r_2 * h7_m2 - fs_231_169 * r_4 * h5_m2 - fs_14700_1859 * r_6 * h3_m2) + e_4 * (fs_243_4199 * h11_m8 - fs_243_1356277 * h11_m2 + fs_1_221 * r_2 * h9_m8 + fs_242_48841 * r_2 * h9_m2 - fs_3024000_193947611 * r_4 * h7_m2 + fs_308_146523 * r_6 * h5_m2 + fs_49_5577 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m1, ph11_m10, ph11_m9, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m1 = ph11_m1[k];

        pc_20[k] = fs_49116375_512 * e_0 * r_2 * h1_m1 + e_1 * (fs_66825_32 * h5_m1 - fs_1002375_32 * r_4 * h1_m1) + e_2 * (fs_86625_1352 * h7_m1 - fs_601425_1352 * r_2 * h5_m1 + fs_12375_8 * r_6 * h1_m1) + e_3 * (-fs_3969_884 * h9_m9 + fs_3969_25432 * h9_m1 - fs_173250_48841 * r_2 * h7_m1 + fs_2673_338 * r_4 * h5_m1 - fs_1125_88 * r_8 * h1_m1) + e_4 * (fs_360_4199 * h11_m9 + fs_540_17631601 * h11_m1 + fs_9_221 * r_2 * h9_m9 - fs_9_6358 * r_2 * h9_m1 + fs_173250_17631601 * r_4 * h7_m1 - fs_594_48841 * r_6 * h5_m1 + fs_45_3718 * r_10 * h1_m1) - fs_49116375_2048 * e_5 * h1_m1;

        pc_21[k] = fs_378_4199 * e_4 * h11_m10;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_22[k] = e_0 * (fs_2679075_256 * h3_p1 - fs_893025_512 * r_2 * h1_p1) + e_1 * (fs_84375_32 * h5_p1 - fs_33075_4 * r_2 * h3_p1 + fs_18225_32 * r_4 * h1_p1) + e_2 * (fs_1929375_40898 * h7_p1 - fs_759375_1352 * r_2 * h5_p1 + fs_297675_484 * r_4 * h3_p1 - fs_225_8 * r_6 * h1_p1) + e_3 * (fs_4465125_47278088 * h9_p1 - fs_19845_9724 * h9_p9 - fs_15435000_5909761 * r_2 * h7_p1 + fs_3375_338 * r_4 * h5_p1 - fs_132300_20449 * r_6 * h3_p1 + fs_225_968 * r_8 * h1_p1) + e_4 * (fs_297_17631601 * h11_p1 - fs_198_4199 * h11_p9 - fs_10125_11819522 * r_2 * h9_p1 + fs_45_2431 * r_2 * h9_p9 + fs_15435000_2133423721 * r_4 * h7_p1 - fs_750_48841 * r_6 * h5_p1 + fs_147_20449 * r_8 * h3_p1 - fs_9_40898 * r_10 * h1_p1) + fs_893025_2048 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph9_p8, ph11_0, ph11_p8, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p8 = ph11_p8[k];

        pc_23[k] = e_0 * (-fs_4465125_256 * h3_0 + fs_4465125_64 * r_2 * h1_0) + e_1 * (fs_1125 * h5_0 + fs_55125_4 * r_2 * h3_0 - fs_91125_4 * r_4 * h1_0) + e_2 * (fs_15931125_81796 * h7_0 - fs_40500_169 * r_2 * h5_0 - fs_496125_484 * r_4 * h3_0 + fs_1125 * r_6 * h1_0) + e_3 * (fs_19845_20449 * h9_0 + fs_3969_2431 * h9_p8 - fs_220500_20449 * r_2 * h7_0 + fs_720_169 * r_4 * h5_0 + fs_220500_20449 * r_6 * h3_0 - fs_1125_121 * r_8 * h1_0) + e_4 * (fs_5445_17631601 * h11_0 - fs_297_4199 * h11_p8 - fs_180_20449 * r_2 * h9_0 - fs_36_2431 * r_2 * h9_p8 + fs_220500_7382089 * r_4 * h7_0 - fs_320_48841 * r_6 * h5_0 - fs_245_20449 * r_8 * h3_0 + fs_180_20449 * r_10 * h1_0) - fs_4465125_256 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_24[k] = e_0 * (fs_1488375_256 * h3_p1 + fs_40186125_512 * r_2 * h1_p1) + e_1 * (fs_12675_32 * h5_p1 - fs_18375_4 * r_2 * h3_p1 - fs_820125_32 * r_4 * h1_p1) + e_2 * (-fs_70875_968 * h7_p1 - fs_165375_1144 * h7_p7 - fs_675_8 * r_2 * h5_p1 + fs_165375_484 * r_4 * h3_p1 + fs_10125_8 * r_6 * h1_p1) + e_3 * (-fs_59397849_47278088 * h9_p1 + fs_370881_165308 * h9_p7 + fs_141750_34969 * r_2 * h7_p1 + fs_330750_41327 * r_2 * h7_p7 + fs_3_2 * r_4 * h5_p1 - fs_73500_20449 * r_6 * h3_p1 - fs_10125_968 * r_8 * h1_p1) + e_4 * (-fs_13365_17631601 * h11_p1 - fs_5346_79781 * h11_p7 + fs_134689_11819522 * r_2 * h9_p1 - fs_841_41327 * r_2 * h9_p7 - fs_141750_12623809 * r_4 * h7_p1 - fs_330750_14919047 * r_4 * h7_p7 - fs_2_867 * r_6 * h5_p1 + fs_245_61347 * r_8 * h3_p1 + fs_405_40898 * r_10 * h1_p1) - fs_40186125_2048 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_25[k] = f_945_16 * e_0 * h3_p2 + e_1 * (-fs_1575 * h5_p2 - f_105_2 * r_2 * h3_p2) + e_2 * (fs_6622875_654368 * h7_p2 - fs_7875_4576 * h7_p6 + fs_56700_169 * r_2 * h5_p2 + f_315_22 * r_4 * h3_p2) + e_3 * (fs_2223963_1074502 * h9_p2 + fs_35721_82654 * h9_p6 - fs_6622875_11819522 * r_2 * h7_p2 + fs_7875_82654 * r_2 * h7_p6 - fs_1008_169 * r_4 * h5_p2 - f_210_143 * r_6 * h3_p2) + e_4 * (fs_3564_1356277 * h11_p2 - fs_3960_79781 * h11_p6 - fs_10086_537251 * r_2 * h9_p2 - fs_162_41327 * r_2 * h9_p6 + fs_6622875_4266847442 * r_4 * h7_p2 - fs_7875_29838094 * r_4 * h7_p6 + fs_448_48841 * r_6 * h5_p2 + f_7_143 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_26[k] = -fs_2083725_128 * e_0 * h3_p3 + e_1 * (fs_12675_32 * h5_p3 - fs_84375_32 * h5_p5 + fs_25725_2 * r_2 * h3_p3) + e_2 * (fs_28125_1936 * h7_p3 + fs_1540125_29744 * h7_p5 - fs_675_8 * r_2 * h5_p3 + fs_759375_1352 * r_2 * h5_p5 - fs_231525_242 * r_4 * h3_p3) + e_3 * (-fs_9529569_4298008 * h9_p3 - fs_46305_330616 * h9_p5 - fs_28125_34969 * r_2 * h7_p3 - fs_1540125_537251 * r_2 * h7_p5 + fs_3_2 * r_4 * h5_p3 - fs_3375_338 * r_4 * h5_p5 + fs_205800_20449 * r_6 * h3_p3) + e_4 * (-fs_9702_1356277 * h11_p3 - fs_41580_1356277 * h11_p5 + fs_21609_1074502 * r_2 * h9_p3 + fs_105_82654 * r_2 * h9_p5 + fs_28125_12623809 * r_4 * h7_p3 + fs_1540125_193947611 * r_4 * h7_p5 - fs_2_867 * r_6 * h5_p3 + fs_750_48841 * r_6 * h5_p5 - fs_686_61347 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_27[k] = fs_1125 * e_1 * h5_m4 + e_2 * (-fs_270000_1859 * h7_m4 - fs_40500_169 * r_2 * h5_m4) + e_3 * (fs_108045_41327 * h9_m4 + fs_4320000_537251 * r_2 * h7_m4 + fs_720_169 * r_4 * h5_m4) + e_4 * (fs_43659_1356277 * h11_m4 - fs_980_41327 * r_2 * h9_m4 - fs_4320000_193947611 * r_4 * h7_m4 - fs_320_48841 * r_6 * h5_m4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_28[k] = -fs_2083725_128 * e_0 * h3_m3 + e_1 * (fs_84375_32 * h5_m5 + fs_12675_32 * h5_m3 + fs_25725_2 * r_2 * h3_m3) + e_2 * (-fs_1540125_29744 * h7_m5 + fs_28125_1936 * h7_m3 - fs_759375_1352 * r_2 * h5_m5 - fs_675_8 * r_2 * h5_m3 - fs_231525_242 * r_4 * h3_m3) + e_3 * (fs_46305_330616 * h9_m5 - fs_9529569_4298008 * h9_m3 + fs_1540125_537251 * r_2 * h7_m5 - fs_28125_34969 * r_2 * h7_m3 + fs_3375_338 * r_4 * h5_m5 + fs_3_2 * r_4 * h5_m3 + fs_205800_20449 * r_6 * h3_m3) + e_4 * (fs_41580_1356277 * h11_m5 - fs_9702_1356277 * h11_m3 - fs_105_82654 * r_2 * h9_m5 + fs_21609_1074502 * r_2 * h9_m3 - fs_1540125_193947611 * r_4 * h7_m5 + fs_28125_12623809 * r_4 * h7_m3 - fs_750_48841 * r_6 * h5_m5 - fs_2_867 * r_6 * h5_m3 - fs_686_61347 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_29[k] = f_945_16 * e_0 * h3_m2 + e_1 * (-fs_1575 * h5_m2 - f_105_2 * r_2 * h3_m2) + e_2 * (fs_7875_4576 * h7_m6 + fs_6622875_654368 * h7_m2 + fs_56700_169 * r_2 * h5_m2 + f_315_22 * r_4 * h3_m2) + e_3 * (-fs_35721_82654 * h9_m6 + fs_2223963_1074502 * h9_m2 - fs_7875_82654 * r_2 * h7_m6 - fs_6622875_11819522 * r_2 * h7_m2 - fs_1008_169 * r_4 * h5_m2 - f_210_143 * r_6 * h3_m2) + e_4 * (fs_3960_79781 * h11_m6 + fs_3564_1356277 * h11_m2 + fs_162_41327 * r_2 * h9_m6 - fs_10086_537251 * r_2 * h9_m2 + fs_7875_29838094 * r_4 * h7_m6 + fs_6622875_4266847442 * r_4 * h7_m2 + fs_448_48841 * r_6 * h5_m2 + f_7_143 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m8, ph9_m7, ph9_m1, ph11_m8, ph11_m7, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_30[k] = e_0 * (fs_1488375_256 * h3_m1 + fs_40186125_512 * r_2 * h1_m1) + e_1 * (fs_12675_32 * h5_m1 - fs_18375_4 * r_2 * h3_m1 - fs_820125_32 * r_4 * h1_m1) + e_2 * (fs_165375_1144 * h7_m7 - fs_70875_968 * h7_m1 - fs_675_8 * r_2 * h5_m1 + fs_165375_484 * r_4 * h3_m1 + fs_10125_8 * r_6 * h1_m1) + e_3 * (-fs_370881_165308 * h9_m7 - fs_59397849_47278088 * h9_m1 - fs_330750_41327 * r_2 * h7_m7 + fs_141750_34969 * r_2 * h7_m1 + fs_3_2 * r_4 * h5_m1 - fs_73500_20449 * r_6 * h3_m1 - fs_10125_968 * r_8 * h1_m1) + e_4 * (fs_5346_79781 * h11_m7 - fs_13365_17631601 * h11_m1 + fs_841_41327 * r_2 * h9_m7 + fs_134689_11819522 * r_2 * h9_m1 + fs_330750_14919047 * r_4 * h7_m7 - fs_141750_12623809 * r_4 * h7_m1 - fs_2_867 * r_6 * h5_m1 + fs_245_61347 * r_8 * h3_m1 + fs_405_40898 * r_10 * h1_m1) - fs_40186125_2048 * e_5 * h1_m1;

        pc_31[k] = -fs_3969_2431 * e_3 * h9_m8 + e_4 * (fs_297_4199 * h11_m8 + fs_36_2431 * r_2 * h9_m8);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m1, ph11_m9, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m1 = ph11_m1[k];

        pc_32[k] = e_0 * (-fs_2679075_256 * h3_m1 + fs_893025_512 * r_2 * h1_m1) + e_1 * (-fs_84375_32 * h5_m1 + fs_33075_4 * r_2 * h3_m1 - fs_18225_32 * r_4 * h1_m1) + e_2 * (-fs_1929375_40898 * h7_m1 + fs_759375_1352 * r_2 * h5_m1 - fs_297675_484 * r_4 * h3_m1 + fs_225_8 * r_6 * h1_m1) + e_3 * (fs_19845_9724 * h9_m9 - fs_4465125_47278088 * h9_m1 + fs_15435000_5909761 * r_2 * h7_m1 - fs_3375_338 * r_4 * h5_m1 + fs_132300_20449 * r_6 * h3_m1 - fs_225_968 * r_8 * h1_m1) + e_4 * (fs_198_4199 * h11_m9 - fs_297_17631601 * h11_m1 - fs_45_2431 * r_2 * h9_m9 + fs_10125_11819522 * r_2 * h9_m1 - fs_15435000_2133423721 * r_4 * h7_m1 + fs_750_48841 * r_6 * h5_m1 - fs_147_20449 * r_8 * h3_m1 + fs_9_40898 * r_10 * h1_m1) - fs_893025_2048 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_33[k] = -f_945_16 * e_0 * h3_p2 + e_1 * (-fs_39375_16 * h5_p2 + f_105_2 * r_2 * h3_p2) + e_2 * (-fs_3472875_40898 * h7_p2 + fs_354375_676 * r_2 * h5_p2 - f_315_22 * r_4 * h3_p2) + e_3 * (-fs_297675_1074502 * h9_p2 - fs_33075_9724 * h9_p8 + fs_27783000_5909761 * r_2 * h7_p2 - fs_1575_169 * r_4 * h5_p2 + f_210_143 * r_6 * h3_p2) + e_4 * (-fs_99_1356277 * h11_p2 - fs_99_4199 * h11_p8 + fs_1350_537251 * r_2 * h9_p2 + fs_75_2431 * r_2 * h9_p8 - fs_27783000_2133423721 * r_4 * h7_p2 + fs_700_48841 * r_6 * h5_p2 - f_7_143 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_34[k] = e_0 * (f_945_8 * h3_p1 - fs_2679075_512 * r_2 * h1_p1) + e_1 * (fs_1125_32 * h5_p1 - f_105 * r_2 * h3_p1 + fs_54675_32 * r_4 * h1_p1) + e_2 * (-fs_1852200_20449 * h7_p1 + fs_99225_1144 * h7_p7 - fs_10125_1352 * r_2 * h5_p1 + f_315_11 * r_4 * h3_p1 - fs_675_8 * r_6 * h1_p1) + e_3 * (-fs_50068935_47278088 * h9_p1 + fs_6615_165308 * h9_p7 + fs_29635200_5909761 * r_2 * h7_p1 - fs_198450_41327 * r_2 * h7_p7 + fs_45_338 * r_4 * h5_p1 - f_420_143 * r_6 * h3_p1 + fs_675_968 * r_8 * h1_p1) + e_4 * (-fs_9900_17631601 * h11_p1 - fs_3960_79781 * h11_p7 + fs_113535_11819522 * r_2 * h9_p1 - fs_15_41327 * r_2 * h9_p7 - fs_29635200_2133423721 * r_4 * h7_p1 + fs_198450_14919047 * r_4 * h7_p7 - fs_10_48841 * r_6 * h5_p1 + f_14_143 * r_8 * h3_p1 - fs_27_40898 * r_10 * h1_p1) + fs_2679075_2048 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ph11_0, ph11_p6, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p6 = ph11_p6[k];

        pc_35[k] = e_0 * (-fs_297675_256 * h3_0 + fs_24111675_256 * r_2 * h1_0) + e_1 * (fs_46875_16 * h5_0 + fs_3675_4 * r_2 * h3_0 - fs_492075_16 * r_4 * h1_0) + e_2 * (-fs_2679075_81796 * h7_0 - fs_14175_286 * h7_p6 - fs_421875_676 * r_2 * h5_0 - fs_33075_484 * r_4 * h3_0 + fs_6075_4 * r_6 * h1_0) + e_3 * (-fs_92907675_23639044 * h9_0 + fs_70560_41327 * h9_p6 + fs_10716300_5909761 * r_2 * h7_0 + fs_113400_41327 * r_2 * h7_p6 + fs_1875_169 * r_4 * h5_0 + fs_14700_20449 * r_6 * h3_0 - fs_6075_484 * r_8 * h1_0) + e_4 * (-fs_81675_17631601 * h11_0 - fs_4950_79781 * h11_p6 + fs_210675_5909761 * r_2 * h9_0 - fs_640_41327 * r_2 * h9_p6 - fs_10716300_2133423721 * r_4 * h7_0 - fs_113400_14919047 * r_4 * h7_p6 - fs_2500_146523 * r_6 * h5_0 - fs_49_61347 * r_8 * h3_0 + fs_243_20449 * r_10 * h1_0) - fs_24111675_1024 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_36[k] = e_0 * (fs_2679075_256 * h3_p1 + fs_8037225_128 * r_2 * h1_p1) + e_1 * (-fs_3375_8 * h5_p1 + fs_39375_16 * h5_p5 - fs_33075_4 * r_2 * h3_p1 - fs_164025_8 * r_4 * h1_p1) + e_2 * (-fs_7498575_654368 * h7_p1 - fs_3973725_59488 * h7_p5 + fs_30375_338 * r_2 * h5_p1 - fs_354375_676 * r_2 * h5_p5 + fs_297675_484 * r_4 * h3_p1 + fs_2025_2 * r_6 * h1_p1) + e_3 * (fs_25245045_11819522 * h9_p1 + fs_275625_165308 * h9_p5 + fs_7498575_11819522 * r_2 * h7_p1 + fs_3973725_1074502 * r_2 * h7_p5 - fs_270_169 * r_4 * h5_p1 + fs_1575_169 * r_4 * h5_p5 - fs_132300_20449 * r_6 * h3_p1 - fs_2025_242 * r_8 * h1_p1) + e_4 * (fs_118800_17631601 * h11_p1 - fs_79200_1356277 * h11_p5 - fs_114490_5909761 * r_2 * h9_p1 - fs_625_41327 * r_2 * h9_p5 - fs_7498575_4266847442 * r_4 * h7_p1 - fs_3973725_387895222 * r_4 * h7_p5 + fs_120_48841 * r_6 * h5_p1 - fs_700_48841 * r_6 * h5_p5 + fs_147_20449 * r_8 * h3_p1 + fs_162_20449 * r_10 * h1_p1) - fs_8037225_512 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_37[k] = e_1 * (-fs_3375_8 * h5_p2 - fs_1125_32 * h5_p4) + e_2 * (fs_10800_169 * h7_p2 - fs_6075_14872 * h7_p4 + fs_30375_338 * r_2 * h5_p2 + fs_10125_1352 * r_2 * h5_p4) + e_3 * (-fs_15435_12716 * h9_p2 + fs_108045_330616 * h9_p4 - fs_172800_48841 * r_2 * h7_p2 + fs_12150_537251 * r_2 * h7_p4 - fs_270_169 * r_4 * h5_p2 - fs_45_338 * r_4 * h5_p4) + e_4 * (-fs_20790_1356277 * h11_p2 - fs_121275_2712554 * h11_p4 + fs_35_3179 * r_2 * h9_p2 - fs_245_82654 * r_2 * h9_p4 + fs_172800_17631601 * r_4 * h7_p2 - fs_12150_193947611 * r_4 * h7_p4 + fs_120_48841 * r_6 * h5_p2 + fs_10_48841 * r_6 * h5_p4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m3, ph7_m3, ph9_m3, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m3 = ph11_m3[k];

        pc_38[k] = -fs_2083725_64 * e_0 * h3_m3 + e_1 * (fs_46875_16 * h5_m3 + fs_25725 * r_2 * h3_m3) + e_2 * (-fs_13861125_163592 * h7_m3 - fs_421875_676 * r_2 * h5_m3 - fs_231525_121 * r_4 * h3_m3) + e_3 * (fs_540225_2149004 * h9_m3 + fs_27722250_5909761 * r_2 * h7_m3 + fs_1875_169 * r_4 * h5_m3 + fs_411600_20449 * r_6 * h3_m3) + e_4 * (fs_77616_1356277 * h11_m3 - fs_1225_537251 * r_2 * h9_m3 - fs_27722250_2133423721 * r_4 * h7_m3 - fs_2500_146523 * r_6 * h5_m3 - fs_1372_61347 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_39[k] = e_1 * (fs_1125_32 * h5_m4 - fs_3375_8 * h5_m2) + e_2 * (fs_6075_14872 * h7_m4 + fs_10800_169 * h7_m2 - fs_10125_1352 * r_2 * h5_m4 + fs_30375_338 * r_2 * h5_m2) + e_3 * (-fs_108045_330616 * h9_m4 - fs_15435_12716 * h9_m2 - fs_12150_537251 * r_2 * h7_m4 - fs_172800_48841 * r_2 * h7_m2 + fs_45_338 * r_4 * h5_m4 - fs_270_169 * r_4 * h5_m2) + e_4 * (fs_121275_2712554 * h11_m4 - fs_20790_1356277 * h11_m2 + fs_245_82654 * r_2 * h9_m4 + fs_35_3179 * r_2 * h9_m2 + fs_12150_193947611 * r_4 * h7_m4 + fs_172800_17631601 * r_4 * h7_m2 - fs_10_48841 * r_6 * h5_m4 + fs_120_48841 * r_6 * h5_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_40[k] = e_0 * (fs_2679075_256 * h3_m1 + fs_8037225_128 * r_2 * h1_m1) + e_1 * (-fs_39375_16 * h5_m5 - fs_3375_8 * h5_m1 - fs_33075_4 * r_2 * h3_m1 - fs_164025_8 * r_4 * h1_m1) + e_2 * (fs_3973725_59488 * h7_m5 - fs_7498575_654368 * h7_m1 + fs_354375_676 * r_2 * h5_m5 + fs_30375_338 * r_2 * h5_m1 + fs_297675_484 * r_4 * h3_m1 + fs_2025_2 * r_6 * h1_m1) + e_3 * (-fs_275625_165308 * h9_m5 + fs_25245045_11819522 * h9_m1 - fs_3973725_1074502 * r_2 * h7_m5 + fs_7498575_11819522 * r_2 * h7_m1 - fs_1575_169 * r_4 * h5_m5 - fs_270_169 * r_4 * h5_m1 - fs_132300_20449 * r_6 * h3_m1 - fs_2025_242 * r_8 * h1_m1) + e_4 * (fs_79200_1356277 * h11_m5 + fs_118800_17631601 * h11_m1 + fs_625_41327 * r_2 * h9_m5 - fs_114490_5909761 * r_2 * h9_m1 + fs_3973725_387895222 * r_4 * h7_m5 - fs_7498575_4266847442 * r_4 * h7_m1 + fs_700_48841 * r_6 * h5_m5 + fs_120_48841 * r_6 * h5_m1 + fs_147_20449 * r_8 * h3_m1 + fs_162_20449 * r_10 * h1_m1) - fs_8037225_512 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_2, pe_3, pe_4, ph7_m6, ph9_m6, ph11_m6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h7_m6 = ph7_m6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h11_m6 = ph11_m6[k];

        pc_41[k] = fs_14175_286 * e_2 * h7_m6 + e_3 * (-fs_70560_41327 * h9_m6 - fs_113400_41327 * r_2 * h7_m6) + e_4 * (fs_4950_79781 * h11_m6 + fs_640_41327 * r_2 * h9_m6 + fs_113400_14919047 * r_4 * h7_m6);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_42[k] = e_0 * (-f_945_8 * h3_m1 + fs_2679075_512 * r_2 * h1_m1) + e_1 * (-fs_1125_32 * h5_m1 + f_105 * r_2 * h3_m1 - fs_54675_32 * r_4 * h1_m1) + e_2 * (-fs_99225_1144 * h7_m7 + fs_1852200_20449 * h7_m1 + fs_10125_1352 * r_2 * h5_m1 - f_315_11 * r_4 * h3_m1 + fs_675_8 * r_6 * h1_m1) + e_3 * (-fs_6615_165308 * h9_m7 + fs_50068935_47278088 * h9_m1 + fs_198450_41327 * r_2 * h7_m7 - fs_29635200_5909761 * r_2 * h7_m1 - fs_45_338 * r_4 * h5_m1 + f_420_143 * r_6 * h3_m1 - fs_675_968 * r_8 * h1_m1) + e_4 * (fs_3960_79781 * h11_m7 + fs_9900_17631601 * h11_m1 + fs_15_41327 * r_2 * h9_m7 - fs_113535_11819522 * r_2 * h9_m1 - fs_198450_14919047 * r_4 * h7_m7 + fs_29635200_2133423721 * r_4 * h7_m1 + fs_10_48841 * r_6 * h5_m1 - f_14_143 * r_8 * h3_m1 + fs_27_40898 * r_10 * h1_m1) - fs_2679075_2048 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_43[k] = f_945_16 * e_0 * h3_m2 + e_1 * (fs_39375_16 * h5_m2 - f_105_2 * r_2 * h3_m2) + e_2 * (fs_3472875_40898 * h7_m2 - fs_354375_676 * r_2 * h5_m2 + f_315_22 * r_4 * h3_m2) + e_3 * (fs_33075_9724 * h9_m8 + fs_297675_1074502 * h9_m2 - fs_27783000_5909761 * r_2 * h7_m2 + fs_1575_169 * r_4 * h5_m2 - f_210_143 * r_6 * h3_m2) + e_4 * (fs_99_4199 * h11_m8 + fs_99_1356277 * h11_m2 - fs_75_2431 * r_2 * h9_m8 - fs_1350_537251 * r_2 * h9_m2 + fs_27783000_2133423721 * r_4 * h7_m2 - fs_700_48841 * r_6 * h5_m2 + f_7_143 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_44[k] = fs_297675_512 * e_0 * h3_p3 + e_1 * (fs_13125_8 * h5_p3 - fs_3675_8 * r_2 * h3_p3) + e_2 * (fs_9646875_81796 * h7_p3 - fs_55125_2288 * h7_p7 - fs_118125_338 * r_2 * h5_p3 + fs_33075_968 * r_4 * h3_p3) + e_3 * (fs_694575_1074502 * h9_p3 - fs_297675_82654 * h9_p7 - fs_38587500_5909761 * r_2 * h7_p3 + fs_55125_41327 * r_2 * h7_p7 + fs_1050_169 * r_4 * h5_p3 - fs_7350_20449 * r_6 * h3_p3) + e_4 * (fs_693_2712554 * h11_p3 - fs_891_79781 * h11_p7 - fs_3150_537251 * r_2 * h9_p3 + fs_1350_41327 * r_2 * h9_p7 + fs_38587500_2133423721 * r_4 * h7_p3 - fs_55125_14919047 * r_4 * h7_p7 - fs_1400_146523 * r_6 * h5_p3 + fs_49_122694 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_45[k] = -fs_4465125_512 * e_0 * h3_p2 + e_1 * (-fs_7875_8 * h5_p2 + fs_55125_8 * r_2 * h3_p2) + e_2 * (fs_3781575_81796 * h7_p2 + fs_20475_176 * h7_p6 + fs_70875_338 * r_2 * h5_p2 - fs_496125_968 * r_4 * h3_p2) + e_3 * (fs_952560_537251 * h9_p2 - fs_19845_41327 * h9_p6 - fs_15126300_5909761 * r_2 * h7_p2 - fs_20475_3179 * r_2 * h7_p6 - fs_630_169 * r_4 * h5_p2 + fs_110250_20449 * r_6 * h3_p2) + e_4 * (fs_4455_2712554 * h11_p2 - fs_2475_79781 * h11_p6 - fs_8640_537251 * r_2 * h9_p2 + fs_180_41327 * r_2 * h9_p6 + fs_15126300_2133423721 * r_4 * h7_p2 + fs_20475_1147619 * r_4 * h7_p6 + fs_280_48841 * r_6 * h5_p2 - fs_245_40898 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_46[k] = e_0 * (fs_4862025_512 * h3_p1 - fs_2679075_256 * r_2 * h1_p1) + e_1 * (-fs_15125_16 * h5_p1 - fs_13125_8 * h5_p5 - fs_60025_8 * r_2 * h3_p1 + fs_54675_16 * r_4 * h1_p1) + e_2 * (-fs_231525_81796 * h7_p1 + fs_14175_29744 * h7_p5 + fs_136125_676 * r_2 * h5_p1 + fs_118125_338 * r_2 * h5_p5 + fs_540225_968 * r_4 * h3_p1 - fs_675_4 * r_6 * h1_p1) + e_3 * (fs_52397415_23639044 * h9_p1 + fs_33075_82654 * h9_p5 + fs_926100_5909761 * r_2 * h7_p1 - fs_14175_537251 * r_2 * h7_p5 - fs_605_169 * r_4 * h5_p1 - fs_1050_169 * r_4 * h5_p5 - fs_120050_20449 * r_6 * h3_p1 + fs_675_484 * r_8 * h1_p1) + e_4 * (fs_200475_35263202 * h11_p1 - fs_66825_1356277 * h11_p5 - fs_118815_5909761 * r_2 * h9_p1 - fs_150_41327 * r_2 * h9_p5 - fs_926100_2133423721 * r_4 * h7_p1 + fs_14175_193947611 * r_4 * h7_p5 + fs_2420_439569 * r_6 * h5_p1 + fs_1400_146523 * r_6 * h5_p5 + fs_2401_368082 * r_8 * h3_p1 - fs_27_20449 * r_10 * h1_p1) + fs_2679075_1024 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ph11_0, ph11_p4, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_47[k] = e_0 * (fs_99225_32 * h3_0 + fs_893025_8 * r_2 * h1_0) + e_1 * (fs_625_2 * h5_0 + fs_7875_8 * h5_p4 - fs_2450 * r_2 * h3_0 - fs_36450 * r_4 * h1_0) + e_2 * (-fs_18533025_163592 * h7_0 - fs_3444525_59488 * h7_p4 - fs_11250_169 * r_2 * h5_0 - fs_70875_338 * r_2 * h5_p4 + fs_22050_121 * r_4 * h3_0 + fs_1800 * r_6 * h1_0) + e_3 * (fs_16074450_5909761 * h9_0 + fs_138915_82654 * h9_p4 + fs_37066050_5909761 * r_2 * h7_0 + fs_3444525_1074502 * r_2 * h7_p4 + fs_200_169 * r_4 * h5_0 + fs_630_169 * r_4 * h5_p4 - fs_39200_20449 * r_6 * h3_0 - fs_1800_121 * r_8 * h1_0) + e_4 * (fs_490050_17631601 * h11_0 - fs_155925_2712554 * h11_p4 - fs_145800_5909761 * r_2 * h9_0 - fs_630_41327 * r_2 * h9_p4 - fs_37066050_2133423721 * r_4 * h7_0 - fs_3444525_387895222 * r_4 * h7_p4 - fs_800_439569 * r_6 * h5_0 - fs_280_48841 * r_6 * h5_p4 + fs_392_184041 * r_8 * h3_0 + fs_288_20449 * r_10 * h1_0) - fs_893025_32 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_48[k] = e_0 * (fs_2083725_256 * h3_p1 - fs_3472875_256 * h3_p3 + fs_6251175_128 * r_2 * h1_p1) + e_1 * (-fs_2625_2 * h5_p1 + fs_15125_16 * h5_p3 - fs_25725_4 * r_2 * h3_p1 + fs_42875_4 * r_2 * h3_p3 - fs_127575_8 * r_4 * h1_p1) + e_2 * (fs_30305025_654368 * h7_p1 - fs_28923075_654368 * h7_p3 + fs_47250_169 * r_2 * h5_p1 - fs_136125_676 * r_2 * h5_p3 + fs_231525_484 * r_4 * h3_p1 - fs_385875_484 * r_4 * h3_p3 + fs_1575_2 * r_6 * h1_p1) + e_3 * (-fs_1111320_5909761 * h9_p1 + fs_2917215_2149004 * h9_p3 - fs_30305025_11819522 * r_2 * h7_p1 + fs_28923075_11819522 * r_2 * h7_p3 - fs_840_169 * r_4 * h5_p1 + fs_605_169 * r_4 * h5_p3 - fs_102900_20449 * r_6 * h3_p1 + fs_171500_20449 * r_6 * h3_p3 - fs_1575_242 * r_8 * h1_p1) + e_4 * (-fs_467775_17631601 * h11_p1 - fs_72765_1356277 * h11_p3 + fs_10080_5909761 * r_2 * h9_p1 - fs_6615_537251 * r_2 * h9_p3 + fs_30305025_4266847442 * r_4 * h7_p1 - fs_28923075_4266847442 * r_4 * h7_p3 + fs_1120_146523 * r_6 * h5_p1 - fs_2420_439569 * r_6 * h5_p3 + fs_343_61347 * r_8 * h3_p1 - fs_1715_184041 * r_8 * h3_p3 + fs_126_20449 * r_10 * h1_p1) - fs_6251175_512 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m2 = ph11_m2[k];

        pc_49[k] = -fs_694575_128 * e_0 * h3_m2 + e_1 * (fs_625_2 * h5_m2 + fs_8575_2 * r_2 * h3_m2) + e_2 * (-fs_91125_327184 * h7_m2 - fs_11250_169 * r_2 * h5_m2 - fs_77175_242 * r_4 * h3_m2) + e_3 * (-fs_231525_537251 * h9_m2 + fs_91125_5909761 * r_2 * h7_m2 + fs_200_169 * r_4 * h5_m2 + fs_68600_20449 * r_6 * h3_m2) + e_4 * (fs_112266_1356277 * h11_m2 + fs_2100_537251 * r_2 * h9_m2 - fs_91125_2133423721 * r_4 * h7_m2 - fs_800_439569 * r_6 * h5_m2 - fs_686_184041 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_50[k] = e_0 * (fs_3472875_256 * h3_m3 + fs_2083725_256 * h3_m1 + fs_6251175_128 * r_2 * h1_m1) + e_1 * (-fs_15125_16 * h5_m3 - fs_2625_2 * h5_m1 - fs_42875_4 * r_2 * h3_m3 - fs_25725_4 * r_2 * h3_m1 - fs_127575_8 * r_4 * h1_m1) + e_2 * (fs_28923075_654368 * h7_m3 + fs_30305025_654368 * h7_m1 + fs_136125_676 * r_2 * h5_m3 + fs_47250_169 * r_2 * h5_m1 + fs_385875_484 * r_4 * h3_m3 + fs_231525_484 * r_4 * h3_m1 + fs_1575_2 * r_6 * h1_m1) + e_3 * (-fs_2917215_2149004 * h9_m3 - fs_1111320_5909761 * h9_m1 - fs_28923075_11819522 * r_2 * h7_m3 - fs_30305025_11819522 * r_2 * h7_m1 - fs_605_169 * r_4 * h5_m3 - fs_840_169 * r_4 * h5_m1 - fs_171500_20449 * r_6 * h3_m3 - fs_102900_20449 * r_6 * h3_m1 - fs_1575_242 * r_8 * h1_m1) + e_4 * (fs_72765_1356277 * h11_m3 - fs_467775_17631601 * h11_m1 + fs_6615_537251 * r_2 * h9_m3 + fs_10080_5909761 * r_2 * h9_m1 + fs_28923075_4266847442 * r_4 * h7_m3 + fs_30305025_4266847442 * r_4 * h7_m1 + fs_2420_439569 * r_6 * h5_m3 + fs_1120_146523 * r_6 * h5_m1 + fs_1715_184041 * r_8 * h3_m3 + fs_343_61347 * r_8 * h3_m1 + fs_126_20449 * r_10 * h1_m1) - fs_6251175_512 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_51[k] = -fs_7875_8 * e_1 * h5_m4 + e_2 * (fs_3444525_59488 * h7_m4 + fs_70875_338 * r_2 * h5_m4) + e_3 * (-fs_138915_82654 * h9_m4 - fs_3444525_1074502 * r_2 * h7_m4 - fs_630_169 * r_4 * h5_m4) + e_4 * (fs_155925_2712554 * h11_m4 + fs_630_41327 * r_2 * h9_m4 + fs_3444525_387895222 * r_4 * h7_m4 + fs_280_48841 * r_6 * h5_m4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_52[k] = e_0 * (-fs_4862025_512 * h3_m1 + fs_2679075_256 * r_2 * h1_m1) + e_1 * (fs_13125_8 * h5_m5 + fs_15125_16 * h5_m1 + fs_60025_8 * r_2 * h3_m1 - fs_54675_16 * r_4 * h1_m1) + e_2 * (-fs_14175_29744 * h7_m5 + fs_231525_81796 * h7_m1 - fs_118125_338 * r_2 * h5_m5 - fs_136125_676 * r_2 * h5_m1 - fs_540225_968 * r_4 * h3_m1 + fs_675_4 * r_6 * h1_m1) + e_3 * (-fs_33075_82654 * h9_m5 - fs_52397415_23639044 * h9_m1 + fs_14175_537251 * r_2 * h7_m5 - fs_926100_5909761 * r_2 * h7_m1 + fs_1050_169 * r_4 * h5_m5 + fs_605_169 * r_4 * h5_m1 + fs_120050_20449 * r_6 * h3_m1 - fs_675_484 * r_8 * h1_m1) + e_4 * (fs_66825_1356277 * h11_m5 - fs_200475_35263202 * h11_m1 + fs_150_41327 * r_2 * h9_m5 + fs_118815_5909761 * r_2 * h9_m1 - fs_14175_193947611 * r_4 * h7_m5 + fs_926100_2133423721 * r_4 * h7_m1 - fs_1400_146523 * r_6 * h5_m5 - fs_2420_439569 * r_6 * h5_m1 - fs_2401_368082 * r_8 * h3_m1 + fs_27_20449 * r_10 * h1_m1) - fs_2679075_1024 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_53[k] = fs_4465125_512 * e_0 * h3_m2 + e_1 * (fs_7875_8 * h5_m2 - fs_55125_8 * r_2 * h3_m2) + e_2 * (-fs_20475_176 * h7_m6 - fs_3781575_81796 * h7_m2 - fs_70875_338 * r_2 * h5_m2 + fs_496125_968 * r_4 * h3_m2) + e_3 * (fs_19845_41327 * h9_m6 - fs_952560_537251 * h9_m2 + fs_20475_3179 * r_2 * h7_m6 + fs_15126300_5909761 * r_2 * h7_m2 + fs_630_169 * r_4 * h5_m2 - fs_110250_20449 * r_6 * h3_m2) + e_4 * (fs_2475_79781 * h11_m6 - fs_4455_2712554 * h11_m2 - fs_180_41327 * r_2 * h9_m6 + fs_8640_537251 * r_2 * h9_m2 - fs_20475_1147619 * r_4 * h7_m6 - fs_15126300_2133423721 * r_4 * h7_m2 - fs_280_48841 * r_6 * h5_m2 + fs_245_40898 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_54[k] = -fs_297675_512 * e_0 * h3_m3 + e_1 * (-fs_13125_8 * h5_m3 + fs_3675_8 * r_2 * h3_m3) + e_2 * (fs_55125_2288 * h7_m7 - fs_9646875_81796 * h7_m3 + fs_118125_338 * r_2 * h5_m3 - fs_33075_968 * r_4 * h3_m3) + e_3 * (fs_297675_82654 * h9_m7 - fs_694575_1074502 * h9_m3 - fs_55125_41327 * r_2 * h7_m7 + fs_38587500_5909761 * r_2 * h7_m3 - fs_1050_169 * r_4 * h5_m3 + fs_7350_20449 * r_6 * h3_m3) + e_4 * (fs_891_79781 * h11_m7 - fs_693_2712554 * h11_m3 - fs_1350_41327 * r_2 * h9_m7 + fs_3150_537251 * r_2 * h9_m3 + fs_55125_14919047 * r_4 * h7_m7 - fs_38587500_2133423721 * r_4 * h7_m3 + fs_1400_146523 * r_6 * h5_m3 - fs_49_122694 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_55[k] = -fs_23625_32 * e_1 * h5_p4 + e_2 * (-fs_1929375_14872 * h7_p4 - fs_39375_572 * h7_p6 + fs_212625_1352 * r_2 * h5_p4) + e_3 * (-fs_416745_330616 * h9_p4 - fs_496125_165308 * h9_p6 + fs_3858750_537251 * r_2 * h7_p4 + fs_157500_41327 * r_2 * h7_p6 - fs_945_338 * r_4 * h5_p4) + e_4 * (-fs_2079_2712554 * h11_p4 - fs_396_79781 * h11_p6 + fs_945_82654 * r_2 * h9_p4 + fs_1125_41327 * r_2 * h9_p6 - fs_3858750_193947611 * r_4 * h7_p4 - fs_157500_14919047 * r_4 * h7_p6 + fs_210_48841 * r_6 * h5_p4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_56[k] = fs_297675_128 * e_0 * h3_p3 + e_1 * (fs_63525_32 * h5_p3 + fs_23625_32 * h5_p5 - fs_3675_2 * r_2 * h3_p3) + e_2 * (-fs_385875_81796 * h7_p3 + fs_126000_1859 * h7_p5 - fs_571725_1352 * r_2 * h5_p3 - fs_212625_1352 * r_2 * h5_p5 + fs_33075_242 * r_4 * h3_p3) + e_3 * (-fs_10029663_4298008 * h9_p3 - fs_535815_330616 * h9_p5 + fs_1543500_5909761 * r_2 * h7_p3 - fs_2016000_537251 * r_2 * h7_p5 + fs_2541_338 * r_4 * h5_p3 + fs_945_338 * r_4 * h5_p5 - fs_29400_20449 * r_6 * h3_p3) + e_4 * (-fs_5544_1356277 * h11_p3 - fs_23760_1356277 * h11_p5 + fs_22743_1074502 * r_2 * h9_p3 + fs_1215_82654 * r_2 * h9_p5 - fs_1543500_2133423721 * r_4 * h7_p3 + fs_2016000_193947611 * r_4 * h7_p5 - fs_1694_146523 * r_6 * h5_p3 - fs_210_48841 * r_6 * h5_p5 + fs_98_61347 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_57[k] = -fs_99225_8 * e_0 * h3_p2 + e_1 * (fs_175_8 * h5_p2 - fs_63525_32 * h5_p4 + fs_9800 * r_2 * h3_p2) + e_2 * (fs_3472875_81796 * h7_p2 + fs_637875_14872 * h7_p4 - fs_1575_338 * r_2 * h5_p2 + fs_571725_1352 * r_2 * h5_p4 - fs_88200_121 * r_4 * h3_p2) + e_3 * (-fs_3716307_2149004 * h9_p2 - fs_9261_330616 * h9_p4 - fs_13891500_5909761 * r_2 * h7_p2 - fs_1275750_537251 * r_2 * h7_p4 + fs_14_169 * r_4 * h5_p2 - fs_2541_338 * r_4 * h5_p4 + fs_156800_20449 * r_6 * h3_p2) + e_4 * (-fs_16038_1356277 * h11_p2 - fs_93555_2712554 * h11_p4 + fs_8427_537251 * r_2 * h9_p2 + fs_21_82654 * r_2 * h9_p4 + fs_13891500_2133423721 * r_4 * h7_p2 + fs_1275750_193947611 * r_4 * h7_p4 - fs_56_439569 * r_6 * h5_p2 + fs_1694_146523 * r_6 * h5_p4 - fs_1568_184041 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_58[k] = e_0 * (fs_1488375_512 * h3_p1 + fs_4862025_512 * h3_p3 - fs_4465125_256 * r_2 * h1_p1) + e_1 * (-fs_1200 * h5_p1 - fs_175_8 * h5_p3 - fs_18375_8 * r_2 * h3_p1 - fs_60025_8 * r_2 * h3_p3 + fs_91125_16 * r_4 * h1_p1) + e_2 * (fs_18907875_327184 * h7_p1 - fs_1913625_327184 * h7_p3 + fs_43200_169 * r_2 * h5_p1 + fs_1575_338 * r_2 * h5_p3 + fs_165375_968 * r_4 * h3_p1 + fs_540225_968 * r_4 * h3_p3 - fs_1125_4 * r_6 * h1_p1) + e_3 * (-f_1449_2431 * h9_p1 + fs_750141_1074502 * h9_p3 - fs_18907875_5909761 * r_2 * h7_p1 + fs_1913625_5909761 * r_2 * h7_p3 - fs_768_169 * r_4 * h5_p1 - fs_14_169 * r_4 * h5_p3 - fs_36750_20449 * r_6 * h3_p1 - fs_120050_20449 * r_6 * h3_p3 + fs_1125_484 * r_8 * h1_p1) + e_4 * (-fs_427680_17631601 * h11_p1 - fs_66528_1356277 * h11_p3 + f_138_2431 * r_2 * h9_p1 - fs_3402_537251 * r_2 * h9_p3 + fs_18907875_2133423721 * r_4 * h7_p1 - fs_1913625_2133423721 * r_4 * h7_p3 + fs_1024_146523 * r_6 * h5_p1 + fs_56_439569 * r_6 * h5_p3 + fs_245_122694 * r_8 * h3_p1 + fs_2401_368082 * r_8 * h3_p3 - fs_45_20449 * r_10 * h1_p1) + fs_4465125_1024 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ph11_0, ph11_p2, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_59[k] = e_0 * (fs_3472875_256 * h3_0 - fs_2083725_256 * h3_p2 + fs_31255875_256 * r_2 * h1_0) + e_1 * (-fs_875 * h5_0 + fs_1200 * h5_p2 - fs_42875_4 * r_2 * h3_0 + fs_25725_4 * r_2 * h3_p2 - fs_637875_16 * r_4 * h1_0) + e_2 * (fs_126000_20449 * h7_0 - fs_2460375_40898 * h7_p2 + fs_31500_169 * r_2 * h5_0 - fs_43200_169 * r_2 * h5_p2 + fs_385875_484 * r_4 * h3_0 - fs_231525_484 * r_4 * h3_p2 + fs_7875_4 * r_6 * h1_0) + e_3 * (fs_1250235_5909761 * h9_0 + fs_889056_537251 * h9_p2 - fs_2016000_5909761 * r_2 * h7_0 + fs_19683000_5909761 * r_2 * h7_p2 - fs_560_169 * r_4 * h5_0 + fs_768_169 * r_4 * h5_p2 - fs_171500_20449 * r_6 * h3_0 + fs_102900_20449 * r_6 * h3_p2 - fs_7875_484 * r_8 * h1_0) + e_4 * (-fs_1372140_17631601 * h11_0 - fs_74844_1356277 * h11_p2 - fs_11340_5909761 * r_2 * h9_0 - fs_8064_537251 * r_2 * h9_p2 + fs_2016000_2133423721 * r_4 * h7_0 - fs_19683000_2133423721 * r_4 * h7_p2 + fs_2240_439569 * r_6 * h5_0 - fs_1024_146523 * r_6 * h5_p2 + fs_1715_184041 * r_8 * h3_0 - fs_343_61347 * r_8 * h3_p2 + fs_315_20449 * r_10 * h1_0) - fs_31255875_1024 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ph11_m2, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m2 = ph11_m2[k];
        const auto h11_m1 = ph11_m1[k];

        pc_60[k] = e_0 * (fs_694575_128 * h3_m1 + fs_18753525_256 * r_2 * h1_m1) + e_1 * (-fs_875 * h5_m1 - fs_8575_2 * r_2 * h3_m1 - fs_382725_16 * r_4 * h1_m1) + e_2 * (fs_4876875_81796 * h7_m1 + fs_31500_169 * r_2 * h5_m1 + fs_77175_242 * r_4 * h3_m1 + fs_4725_4 * r_6 * h1_m1) + e_3 * (-fs_46305_20449 * h9_m1 - fs_67500_20449 * r_2 * h7_m1 - fs_560_169 * r_4 * h5_m1 - fs_68600_20449 * r_6 * h3_m1 - fs_4725_484 * r_8 * h1_m1) + e_4 * (fs_1796256_17631601 * h11_m1 + fs_420_20449 * r_2 * h9_m1 + fs_67500_7382089 * r_4 * h7_m1 + fs_2240_439569 * r_6 * h5_m1 + fs_686_184041 * r_8 * h3_m1 + fs_189_20449 * r_10 * h1_m1) - fs_18753525_1024 * e_5 * h1_m1;

        pc_61[k] = fs_2083725_256 * e_0 * h3_m2 + e_1 * (-fs_1200 * h5_m2 - fs_25725_4 * r_2 * h3_m2) + e_2 * (fs_2460375_40898 * h7_m2 + fs_43200_169 * r_2 * h5_m2 + fs_231525_484 * r_4 * h3_m2) + e_3 * (-fs_889056_537251 * h9_m2 - fs_19683000_5909761 * r_2 * h7_m2 - fs_768_169 * r_4 * h5_m2 - fs_102900_20449 * r_6 * h3_m2) + e_4 * (fs_74844_1356277 * h11_m2 + fs_8064_537251 * r_2 * h9_m2 + fs_19683000_2133423721 * r_4 * h7_m2 + fs_1024_146523 * r_6 * h5_m2 + fs_343_61347 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_62[k] = e_0 * (-fs_4862025_512 * h3_m3 - fs_1488375_512 * h3_m1 + fs_4465125_256 * r_2 * h1_m1) + e_1 * (fs_175_8 * h5_m3 + fs_1200 * h5_m1 + fs_60025_8 * r_2 * h3_m3 + fs_18375_8 * r_2 * h3_m1 - fs_91125_16 * r_4 * h1_m1) + e_2 * (fs_1913625_327184 * h7_m3 - fs_18907875_327184 * h7_m1 - fs_1575_338 * r_2 * h5_m3 - fs_43200_169 * r_2 * h5_m1 - fs_540225_968 * r_4 * h3_m3 - fs_165375_968 * r_4 * h3_m1 + fs_1125_4 * r_6 * h1_m1) + e_3 * (-fs_750141_1074502 * h9_m3 + f_1449_2431 * h9_m1 - fs_1913625_5909761 * r_2 * h7_m3 + fs_18907875_5909761 * r_2 * h7_m1 + fs_14_169 * r_4 * h5_m3 + fs_768_169 * r_4 * h5_m1 + fs_120050_20449 * r_6 * h3_m3 + fs_36750_20449 * r_6 * h3_m1 - fs_1125_484 * r_8 * h1_m1) + e_4 * (fs_66528_1356277 * h11_m3 + fs_427680_17631601 * h11_m1 + fs_3402_537251 * r_2 * h9_m3 - f_138_2431 * r_2 * h9_m1 + fs_1913625_2133423721 * r_4 * h7_m3 - fs_18907875_2133423721 * r_4 * h7_m1 - fs_56_439569 * r_6 * h5_m3 - fs_1024_146523 * r_6 * h5_m1 - fs_2401_368082 * r_8 * h3_m3 - fs_245_122694 * r_8 * h3_m1 + fs_45_20449 * r_10 * h1_m1) - fs_4465125_1024 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_63[k] = fs_99225_8 * e_0 * h3_m2 + e_1 * (fs_63525_32 * h5_m4 - fs_175_8 * h5_m2 - fs_9800 * r_2 * h3_m2) + e_2 * (-fs_637875_14872 * h7_m4 - fs_3472875_81796 * h7_m2 - fs_571725_1352 * r_2 * h5_m4 + fs_1575_338 * r_2 * h5_m2 + fs_88200_121 * r_4 * h3_m2) + e_3 * (fs_9261_330616 * h9_m4 + fs_3716307_2149004 * h9_m2 + fs_1275750_537251 * r_2 * h7_m4 + fs_13891500_5909761 * r_2 * h7_m2 + fs_2541_338 * r_4 * h5_m4 - fs_14_169 * r_4 * h5_m2 - fs_156800_20449 * r_6 * h3_m2) + e_4 * (fs_93555_2712554 * h11_m4 + fs_16038_1356277 * h11_m2 - fs_21_82654 * r_2 * h9_m4 - fs_8427_537251 * r_2 * h9_m2 - fs_1275750_193947611 * r_4 * h7_m4 - fs_13891500_2133423721 * r_4 * h7_m2 - fs_1694_146523 * r_6 * h5_m4 + fs_56_439569 * r_6 * h5_m2 + fs_1568_184041 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_64[k] = -fs_297675_128 * e_0 * h3_m3 + e_1 * (-fs_23625_32 * h5_m5 - fs_63525_32 * h5_m3 + fs_3675_2 * r_2 * h3_m3) + e_2 * (-fs_126000_1859 * h7_m5 + fs_385875_81796 * h7_m3 + fs_212625_1352 * r_2 * h5_m5 + fs_571725_1352 * r_2 * h5_m3 - fs_33075_242 * r_4 * h3_m3) + e_3 * (fs_535815_330616 * h9_m5 + fs_10029663_4298008 * h9_m3 + fs_2016000_537251 * r_2 * h7_m5 - fs_1543500_5909761 * r_2 * h7_m3 - fs_945_338 * r_4 * h5_m5 - fs_2541_338 * r_4 * h5_m3 + fs_29400_20449 * r_6 * h3_m3) + e_4 * (fs_23760_1356277 * h11_m5 + fs_5544_1356277 * h11_m3 - fs_1215_82654 * r_2 * h9_m5 - fs_22743_1074502 * r_2 * h9_m3 - fs_2016000_193947611 * r_4 * h7_m5 + fs_1543500_2133423721 * r_4 * h7_m3 + fs_210_48841 * r_6 * h5_m5 + fs_1694_146523 * r_6 * h5_m3 - fs_98_61347 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m5, ph5_m4, ph7_m6, ph7_m5, ph7_m4, ph9_m6, ph9_m5, ph9_m4, ph11_m6, ph11_m5, ph11_m4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_65[k] = fs_23625_32 * e_1 * h5_m4 + e_2 * (fs_39375_572 * h7_m6 + fs_1929375_14872 * h7_m4 - fs_212625_1352 * r_2 * h5_m4) + e_3 * (fs_496125_165308 * h9_m6 + fs_416745_330616 * h9_m4 - fs_157500_41327 * r_2 * h7_m6 - fs_3858750_537251 * r_2 * h7_m4 + fs_945_338 * r_4 * h5_m4) + e_4 * (fs_396_79781 * h11_m6 + fs_2079_2712554 * h11_m4 - fs_1125_41327 * r_2 * h9_m6 - fs_945_82654 * r_2 * h9_m4 + fs_157500_14919047 * r_4 * h7_m6 + fs_3858750_193947611 * r_4 * h7_m4 - fs_210_48841 * r_6 * h5_m4);

        pc_66[k] = f_75_4 * e_1 * h5_m5 + e_2 * (fs_826875_3718 * h7_m5 - f_225_26 * r_2 * h5_m5) + e_3 * (fs_694575_165308 * h9_m5 - fs_6615000_537251 * r_2 * h7_m5 + f_15_13 * r_4 * h5_m5) + e_4 * (fs_5544_1356277 * h11_m5 - fs_1575_41327 * r_2 * h9_m5 + fs_6615000_193947611 * r_4 * h7_m5 - f_10_221 * r_6 * h5_m5);

        pc_67[k] = -f_60 * e_1 * h5_m4 + e_2 * (-fs_165375_7436 * h7_m4 + f_360_13 * r_2 * h5_m4) + e_3 * (fs_194481_41327 * h9_m4 + fs_661500_537251 * r_2 * h7_m4 - f_48_13 * r_4 * h5_m4) + e_4 * (fs_24255_1356277 * h11_m4 - fs_1764_41327 * r_2 * h9_m4 - fs_661500_193947611 * r_4 * h7_m4 + f_32_221 * r_6 * h5_m4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m3, ph7_m2, ph9_m3, ph9_m2, ph11_m3, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

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

        pc_68[k] = fs_694575_64 * e_0 * h3_m3 + e_1 * (f_145_4 * h5_m3 - fs_8575 * r_2 * h3_m3) + e_2 * (-fs_2976750_20449 * h7_m3 - f_435_26 * r_2 * h5_m3 + fs_77175_121 * r_4 * h3_m3) + e_3 * (fs_3176523_2149004 * h9_m3 + fs_47628000_5909761 * r_2 * h7_m3 + f_29_13 * r_4 * h5_m3 - fs_137200_20449 * r_6 * h3_m3) + e_4 * (fs_58212_1356277 * h11_m3 - fs_7203_537251 * r_2 * h9_m3 - fs_47628000_2133423721 * r_4 * h7_m3 - f_58_663 * r_6 * h5_m3 + fs_1372_184041 * r_8 * h3_m3);

        pc_69[k] = -fs_6251175_256 * e_0 * h3_m2 + e_1 * (f_45 * h5_m2 + fs_77175_4 * r_2 * h3_m2) + e_2 * (-fs_4465125_163592 * h7_m2 - f_270_13 * r_2 * h5_m2 - fs_694575_484 * r_4 * h3_m2) + e_3 * (-fs_18522_537251 * h9_m2 + fs_8930250_5909761 * r_2 * h7_m2 + f_36_13 * r_4 * h5_m2 + fs_308700_20449 * r_6 * h3_m2) + e_4 * (fs_99792_1356277 * h11_m2 + fs_168_537251 * r_2 * h9_m2 - fs_8930250_2133423721 * r_4 * h7_m2 - f_24_221 * r_6 * h5_m2 - fs_343_20449 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph1_0, ph3_0, ph5_m1, ph5_0, ph7_m1, ph7_0, ph9_m1, ph9_0, ph11_m1, ph11_0, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_m1 = ph11_m1[k];
        const auto h11_0 = ph11_0[k];

        pc_70[k] = -fs_13395375_256 * e_0 * r_2 * h1_m1 + e_1 * (-f_15 * h5_m1 + fs_273375_16 * r_4 * h1_m1) + e_2 * (fs_23625_676 * h7_m1 + f_90_13 * r_2 * h5_m1 - fs_3375_4 * r_6 * h1_m1) + e_3 * (-fs_64827_34969 * h9_m1 - fs_94500_48841 * r_2 * h7_m1 - f_12_13 * r_4 * h5_m1 + fs_3375_484 * r_8 * h1_m1) + e_4 * (fs_1746360_17631601 * h11_m1 + fs_588_34969 * r_2 * h9_m1 + fs_94500_17631601 * r_4 * h7_m1 + f_8_221 * r_6 * h5_m1 - fs_135_20449 * r_10 * h1_m1) + fs_13395375_1024 * e_5 * h1_m1;

        pc_71[k] = e_0 * (f_2205_16 * h3_0 + f_2835_8 * r_2 * h1_0) + e_1 * (-f_50 * h5_0 - f_245_2 * r_2 * h3_0 - f_405_2 * r_4 * h1_0) + e_2 * (f_1575_143 * h7_0 + f_300_13 * r_2 * h5_0 + f_735_22 * r_4 * h3_0 + f_45 * r_6 * h1_0) + e_3 * (-f_4410_2431 * h9_0 - f_6300_2431 * r_2 * h7_0 - f_40_13 * r_4 * h5_0 - f_490_143 * r_6 * h3_0 - f_45_11 * r_8 * h1_0) + e_4 * (f_1386_4199 * h11_0 + f_420_2431 * r_2 * h9_0 + f_6300_46189 * r_4 * h7_0 + f_80_663 * r_6 * h5_0 + f_49_429 * r_8 * h3_0 + f_18_143 * r_10 * h1_0) - f_2835_16 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph9_p1, ph9_p2, ph11_p1, ph11_p2, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_72[k] = -fs_13395375_256 * e_0 * r_2 * h1_p1 + e_1 * (-f_15 * h5_p1 + fs_273375_16 * r_4 * h1_p1) + e_2 * (fs_23625_676 * h7_p1 + f_90_13 * r_2 * h5_p1 - fs_3375_4 * r_6 * h1_p1) + e_3 * (-fs_64827_34969 * h9_p1 - fs_94500_48841 * r_2 * h7_p1 - f_12_13 * r_4 * h5_p1 + fs_3375_484 * r_8 * h1_p1) + e_4 * (fs_1746360_17631601 * h11_p1 + fs_588_34969 * r_2 * h9_p1 + fs_94500_17631601 * r_4 * h7_p1 + f_8_221 * r_6 * h5_p1 - fs_135_20449 * r_10 * h1_p1) + fs_13395375_1024 * e_5 * h1_p1;

        pc_73[k] = -fs_6251175_256 * e_0 * h3_p2 + e_1 * (f_45 * h5_p2 + fs_77175_4 * r_2 * h3_p2) + e_2 * (-fs_4465125_163592 * h7_p2 - f_270_13 * r_2 * h5_p2 - fs_694575_484 * r_4 * h3_p2) + e_3 * (-fs_18522_537251 * h9_p2 + fs_8930250_5909761 * r_2 * h7_p2 + f_36_13 * r_4 * h5_p2 + fs_308700_20449 * r_6 * h3_p2) + e_4 * (fs_99792_1356277 * h11_p2 + fs_168_537251 * r_2 * h9_p2 - fs_8930250_2133423721 * r_4 * h7_p2 - f_24_221 * r_6 * h5_p2 - fs_343_20449 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph9_p3, ph9_p4, ph11_p3, ph11_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p4 = ph11_p4[k];

        pc_74[k] = fs_694575_64 * e_0 * h3_p3 + e_1 * (f_145_4 * h5_p3 - fs_8575 * r_2 * h3_p3) + e_2 * (-fs_2976750_20449 * h7_p3 - f_435_26 * r_2 * h5_p3 + fs_77175_121 * r_4 * h3_p3) + e_3 * (fs_3176523_2149004 * h9_p3 + fs_47628000_5909761 * r_2 * h7_p3 + f_29_13 * r_4 * h5_p3 - fs_137200_20449 * r_6 * h3_p3) + e_4 * (fs_58212_1356277 * h11_p3 - fs_7203_537251 * r_2 * h9_p3 - fs_47628000_2133423721 * r_4 * h7_p3 - f_58_663 * r_6 * h5_p3 + fs_1372_184041 * r_8 * h3_p3);

        pc_75[k] = -f_60 * e_1 * h5_p4 + e_2 * (-fs_165375_7436 * h7_p4 + f_360_13 * r_2 * h5_p4) + e_3 * (fs_194481_41327 * h9_p4 + fs_661500_537251 * r_2 * h7_p4 - f_48_13 * r_4 * h5_p4) + e_4 * (fs_24255_1356277 * h11_p4 - fs_1764_41327 * r_2 * h9_p4 - fs_661500_193947611 * r_4 * h7_p4 + f_32_221 * r_6 * h5_p4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph5_p5, ph7_m6, ph7_m4, ph7_p5, ph9_m6, ph9_m4, ph9_p5, ph11_m6, ph11_m4, ph11_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_76[k] = f_75_4 * e_1 * h5_p5 + e_2 * (fs_826875_3718 * h7_p5 - f_225_26 * r_2 * h5_p5) + e_3 * (fs_694575_165308 * h9_p5 - fs_6615000_537251 * r_2 * h7_p5 + f_15_13 * r_4 * h5_p5) + e_4 * (fs_5544_1356277 * h11_p5 - fs_1575_41327 * r_2 * h9_p5 + fs_6615000_193947611 * r_4 * h7_p5 - f_10_221 * r_6 * h5_p5);

        pc_77[k] = -fs_23625_32 * e_1 * h5_m4 + e_2 * (fs_39375_572 * h7_m6 - fs_1929375_14872 * h7_m4 + fs_212625_1352 * r_2 * h5_m4) + e_3 * (fs_496125_165308 * h9_m6 - fs_416745_330616 * h9_m4 - fs_157500_41327 * r_2 * h7_m6 + fs_3858750_537251 * r_2 * h7_m4 - fs_945_338 * r_4 * h5_m4) + e_4 * (fs_396_79781 * h11_m6 - fs_2079_2712554 * h11_m4 - fs_1125_41327 * r_2 * h9_m6 + fs_945_82654 * r_2 * h9_m4 + fs_157500_14919047 * r_4 * h7_m6 - fs_3858750_193947611 * r_4 * h7_m4 + fs_210_48841 * r_6 * h5_m4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_78[k] = fs_297675_128 * e_0 * h3_m3 + e_1 * (-fs_23625_32 * h5_m5 + fs_63525_32 * h5_m3 - fs_3675_2 * r_2 * h3_m3) + e_2 * (-fs_126000_1859 * h7_m5 - fs_385875_81796 * h7_m3 + fs_212625_1352 * r_2 * h5_m5 - fs_571725_1352 * r_2 * h5_m3 + fs_33075_242 * r_4 * h3_m3) + e_3 * (fs_535815_330616 * h9_m5 - fs_10029663_4298008 * h9_m3 + fs_2016000_537251 * r_2 * h7_m5 + fs_1543500_5909761 * r_2 * h7_m3 - fs_945_338 * r_4 * h5_m5 + fs_2541_338 * r_4 * h5_m3 - fs_29400_20449 * r_6 * h3_m3) + e_4 * (fs_23760_1356277 * h11_m5 - fs_5544_1356277 * h11_m3 - fs_1215_82654 * r_2 * h9_m5 + fs_22743_1074502 * r_2 * h9_m3 - fs_2016000_193947611 * r_4 * h7_m5 - fs_1543500_2133423721 * r_4 * h7_m3 + fs_210_48841 * r_6 * h5_m5 - fs_1694_146523 * r_6 * h5_m3 + fs_98_61347 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_79[k] = -fs_99225_8 * e_0 * h3_m2 + e_1 * (fs_63525_32 * h5_m4 + fs_175_8 * h5_m2 + fs_9800 * r_2 * h3_m2) + e_2 * (-fs_637875_14872 * h7_m4 + fs_3472875_81796 * h7_m2 - fs_571725_1352 * r_2 * h5_m4 - fs_1575_338 * r_2 * h5_m2 - fs_88200_121 * r_4 * h3_m2) + e_3 * (fs_9261_330616 * h9_m4 - fs_3716307_2149004 * h9_m2 + fs_1275750_537251 * r_2 * h7_m4 - fs_13891500_5909761 * r_2 * h7_m2 + fs_2541_338 * r_4 * h5_m4 + fs_14_169 * r_4 * h5_m2 + fs_156800_20449 * r_6 * h3_m2) + e_4 * (fs_93555_2712554 * h11_m4 - fs_16038_1356277 * h11_m2 - fs_21_82654 * r_2 * h9_m4 + fs_8427_537251 * r_2 * h9_m2 - fs_1275750_193947611 * r_4 * h7_m4 + fs_13891500_2133423721 * r_4 * h7_m2 - fs_1694_146523 * r_6 * h5_m4 - fs_56_439569 * r_6 * h5_m2 - fs_1568_184041 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_80[k] = e_0 * (-fs_4862025_512 * h3_m3 + fs_1488375_512 * h3_m1 - fs_4465125_256 * r_2 * h1_m1) + e_1 * (fs_175_8 * h5_m3 - fs_1200 * h5_m1 + fs_60025_8 * r_2 * h3_m3 - fs_18375_8 * r_2 * h3_m1 + fs_91125_16 * r_4 * h1_m1) + e_2 * (fs_1913625_327184 * h7_m3 + fs_18907875_327184 * h7_m1 - fs_1575_338 * r_2 * h5_m3 + fs_43200_169 * r_2 * h5_m1 - fs_540225_968 * r_4 * h3_m3 + fs_165375_968 * r_4 * h3_m1 - fs_1125_4 * r_6 * h1_m1) + e_3 * (-fs_750141_1074502 * h9_m3 - f_1449_2431 * h9_m1 - fs_1913625_5909761 * r_2 * h7_m3 - fs_18907875_5909761 * r_2 * h7_m1 + fs_14_169 * r_4 * h5_m3 - fs_768_169 * r_4 * h5_m1 + fs_120050_20449 * r_6 * h3_m3 - fs_36750_20449 * r_6 * h3_m1 + fs_1125_484 * r_8 * h1_m1) + e_4 * (fs_66528_1356277 * h11_m3 - fs_427680_17631601 * h11_m1 + fs_3402_537251 * r_2 * h9_m3 + f_138_2431 * r_2 * h9_m1 + fs_1913625_2133423721 * r_4 * h7_m3 + fs_18907875_2133423721 * r_4 * h7_m1 - fs_56_439569 * r_6 * h5_m3 + fs_1024_146523 * r_6 * h5_m1 - fs_2401_368082 * r_8 * h3_m3 + fs_245_122694 * r_8 * h3_m1 - fs_45_20449 * r_10 * h1_m1) + fs_4465125_1024 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_m2, ph3_p1, ph5_m2, ph5_p1, ph7_m2, ph7_p1, ph9_m2, ph9_p1, ph11_m2, ph11_p1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_m2 = ph11_m2[k];
        const auto h11_p1 = ph11_p1[k];

        pc_81[k] = fs_2083725_256 * e_0 * h3_m2 + e_1 * (-fs_1200 * h5_m2 - fs_25725_4 * r_2 * h3_m2) + e_2 * (fs_2460375_40898 * h7_m2 + fs_43200_169 * r_2 * h5_m2 + fs_231525_484 * r_4 * h3_m2) + e_3 * (-fs_889056_537251 * h9_m2 - fs_19683000_5909761 * r_2 * h7_m2 - fs_768_169 * r_4 * h5_m2 - fs_102900_20449 * r_6 * h3_m2) + e_4 * (fs_74844_1356277 * h11_m2 + fs_8064_537251 * r_2 * h9_m2 + fs_19683000_2133423721 * r_4 * h7_m2 + fs_1024_146523 * r_6 * h5_m2 + fs_343_61347 * r_8 * h3_m2);

        pc_82[k] = e_0 * (fs_694575_128 * h3_p1 + fs_18753525_256 * r_2 * h1_p1) + e_1 * (-fs_875 * h5_p1 - fs_8575_2 * r_2 * h3_p1 - fs_382725_16 * r_4 * h1_p1) + e_2 * (fs_4876875_81796 * h7_p1 + fs_31500_169 * r_2 * h5_p1 + fs_77175_242 * r_4 * h3_p1 + fs_4725_4 * r_6 * h1_p1) + e_3 * (-fs_46305_20449 * h9_p1 - fs_67500_20449 * r_2 * h7_p1 - fs_560_169 * r_4 * h5_p1 - fs_68600_20449 * r_6 * h3_p1 - fs_4725_484 * r_8 * h1_p1) + e_4 * (fs_1796256_17631601 * h11_p1 + fs_420_20449 * r_2 * h9_p1 + fs_67500_7382089 * r_4 * h7_p1 + fs_2240_439569 * r_6 * h5_p1 + fs_686_184041 * r_8 * h3_p1 + fs_189_20449 * r_10 * h1_p1) - fs_18753525_1024 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ph11_0, ph11_p2, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_83[k] = e_0 * (fs_3472875_256 * h3_0 + fs_2083725_256 * h3_p2 + fs_31255875_256 * r_2 * h1_0) + e_1 * (-fs_875 * h5_0 - fs_1200 * h5_p2 - fs_42875_4 * r_2 * h3_0 - fs_25725_4 * r_2 * h3_p2 - fs_637875_16 * r_4 * h1_0) + e_2 * (fs_126000_20449 * h7_0 + fs_2460375_40898 * h7_p2 + fs_31500_169 * r_2 * h5_0 + fs_43200_169 * r_2 * h5_p2 + fs_385875_484 * r_4 * h3_0 + fs_231525_484 * r_4 * h3_p2 + fs_7875_4 * r_6 * h1_0) + e_3 * (fs_1250235_5909761 * h9_0 - fs_889056_537251 * h9_p2 - fs_2016000_5909761 * r_2 * h7_0 - fs_19683000_5909761 * r_2 * h7_p2 - fs_560_169 * r_4 * h5_0 - fs_768_169 * r_4 * h5_p2 - fs_171500_20449 * r_6 * h3_0 - fs_102900_20449 * r_6 * h3_p2 - fs_7875_484 * r_8 * h1_0) + e_4 * (-fs_1372140_17631601 * h11_0 + fs_74844_1356277 * h11_p2 - fs_11340_5909761 * r_2 * h9_0 + fs_8064_537251 * r_2 * h9_p2 + fs_2016000_2133423721 * r_4 * h7_0 + fs_19683000_2133423721 * r_4 * h7_p2 + fs_2240_439569 * r_6 * h5_0 + fs_1024_146523 * r_6 * h5_p2 + fs_1715_184041 * r_8 * h3_0 + fs_343_61347 * r_8 * h3_p2 + fs_315_20449 * r_10 * h1_0) - fs_31255875_1024 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_84[k] = e_0 * (fs_1488375_512 * h3_p1 - fs_4862025_512 * h3_p3 - fs_4465125_256 * r_2 * h1_p1) + e_1 * (-fs_1200 * h5_p1 + fs_175_8 * h5_p3 - fs_18375_8 * r_2 * h3_p1 + fs_60025_8 * r_2 * h3_p3 + fs_91125_16 * r_4 * h1_p1) + e_2 * (fs_18907875_327184 * h7_p1 + fs_1913625_327184 * h7_p3 + fs_43200_169 * r_2 * h5_p1 - fs_1575_338 * r_2 * h5_p3 + fs_165375_968 * r_4 * h3_p1 - fs_540225_968 * r_4 * h3_p3 - fs_1125_4 * r_6 * h1_p1) + e_3 * (-f_1449_2431 * h9_p1 - fs_750141_1074502 * h9_p3 - fs_18907875_5909761 * r_2 * h7_p1 - fs_1913625_5909761 * r_2 * h7_p3 - fs_768_169 * r_4 * h5_p1 + fs_14_169 * r_4 * h5_p3 - fs_36750_20449 * r_6 * h3_p1 + fs_120050_20449 * r_6 * h3_p3 + fs_1125_484 * r_8 * h1_p1) + e_4 * (-fs_427680_17631601 * h11_p1 + fs_66528_1356277 * h11_p3 + f_138_2431 * r_2 * h9_p1 + fs_3402_537251 * r_2 * h9_p3 + fs_18907875_2133423721 * r_4 * h7_p1 + fs_1913625_2133423721 * r_4 * h7_p3 + fs_1024_146523 * r_6 * h5_p1 - fs_56_439569 * r_6 * h5_p3 + fs_245_122694 * r_8 * h3_p1 - fs_2401_368082 * r_8 * h3_p3 - fs_45_20449 * r_10 * h1_p1) + fs_4465125_1024 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_85[k] = -fs_99225_8 * e_0 * h3_p2 + e_1 * (fs_175_8 * h5_p2 + fs_63525_32 * h5_p4 + fs_9800 * r_2 * h3_p2) + e_2 * (fs_3472875_81796 * h7_p2 - fs_637875_14872 * h7_p4 - fs_1575_338 * r_2 * h5_p2 - fs_571725_1352 * r_2 * h5_p4 - fs_88200_121 * r_4 * h3_p2) + e_3 * (-fs_3716307_2149004 * h9_p2 + fs_9261_330616 * h9_p4 - fs_13891500_5909761 * r_2 * h7_p2 + fs_1275750_537251 * r_2 * h7_p4 + fs_14_169 * r_4 * h5_p2 + fs_2541_338 * r_4 * h5_p4 + fs_156800_20449 * r_6 * h3_p2) + e_4 * (-fs_16038_1356277 * h11_p2 + fs_93555_2712554 * h11_p4 + fs_8427_537251 * r_2 * h9_p2 - fs_21_82654 * r_2 * h9_p4 + fs_13891500_2133423721 * r_4 * h7_p2 - fs_1275750_193947611 * r_4 * h7_p4 - fs_56_439569 * r_6 * h5_p2 - fs_1694_146523 * r_6 * h5_p4 - fs_1568_184041 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_86[k] = fs_297675_128 * e_0 * h3_p3 + e_1 * (fs_63525_32 * h5_p3 - fs_23625_32 * h5_p5 - fs_3675_2 * r_2 * h3_p3) + e_2 * (-fs_385875_81796 * h7_p3 - fs_126000_1859 * h7_p5 - fs_571725_1352 * r_2 * h5_p3 + fs_212625_1352 * r_2 * h5_p5 + fs_33075_242 * r_4 * h3_p3) + e_3 * (-fs_10029663_4298008 * h9_p3 + fs_535815_330616 * h9_p5 + fs_1543500_5909761 * r_2 * h7_p3 + fs_2016000_537251 * r_2 * h7_p5 + fs_2541_338 * r_4 * h5_p3 - fs_945_338 * r_4 * h5_p5 - fs_29400_20449 * r_6 * h3_p3) + e_4 * (-fs_5544_1356277 * h11_p3 + fs_23760_1356277 * h11_p5 + fs_22743_1074502 * r_2 * h9_p3 - fs_1215_82654 * r_2 * h9_p5 - fs_1543500_2133423721 * r_4 * h7_p3 - fs_2016000_193947611 * r_4 * h7_p5 - fs_1694_146523 * r_6 * h5_p3 + fs_210_48841 * r_6 * h5_p5 + fs_98_61347 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_87[k] = -fs_23625_32 * e_1 * h5_p4 + e_2 * (-fs_1929375_14872 * h7_p4 + fs_39375_572 * h7_p6 + fs_212625_1352 * r_2 * h5_p4) + e_3 * (-fs_416745_330616 * h9_p4 + fs_496125_165308 * h9_p6 + fs_3858750_537251 * r_2 * h7_p4 - fs_157500_41327 * r_2 * h7_p6 - fs_945_338 * r_4 * h5_p4) + e_4 * (-fs_2079_2712554 * h11_p4 + fs_396_79781 * h11_p6 + fs_945_82654 * r_2 * h9_p4 - fs_1125_41327 * r_2 * h9_p6 - fs_3858750_193947611 * r_4 * h7_p4 + fs_157500_14919047 * r_4 * h7_p6 + fs_210_48841 * r_6 * h5_p4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_88[k] = fs_297675_512 * e_0 * h3_m3 + e_1 * (fs_13125_8 * h5_m3 - fs_3675_8 * r_2 * h3_m3) + e_2 * (fs_55125_2288 * h7_m7 + fs_9646875_81796 * h7_m3 - fs_118125_338 * r_2 * h5_m3 + fs_33075_968 * r_4 * h3_m3) + e_3 * (fs_297675_82654 * h9_m7 + fs_694575_1074502 * h9_m3 - fs_55125_41327 * r_2 * h7_m7 - fs_38587500_5909761 * r_2 * h7_m3 + fs_1050_169 * r_4 * h5_m3 - fs_7350_20449 * r_6 * h3_m3) + e_4 * (fs_891_79781 * h11_m7 + fs_693_2712554 * h11_m3 - fs_1350_41327 * r_2 * h9_m7 - fs_3150_537251 * r_2 * h9_m3 + fs_55125_14919047 * r_4 * h7_m7 + fs_38587500_2133423721 * r_4 * h7_m3 - fs_1400_146523 * r_6 * h5_m3 + fs_49_122694 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_89[k] = -fs_4465125_512 * e_0 * h3_m2 + e_1 * (-fs_7875_8 * h5_m2 + fs_55125_8 * r_2 * h3_m2) + e_2 * (-fs_20475_176 * h7_m6 + fs_3781575_81796 * h7_m2 + fs_70875_338 * r_2 * h5_m2 - fs_496125_968 * r_4 * h3_m2) + e_3 * (fs_19845_41327 * h9_m6 + fs_952560_537251 * h9_m2 + fs_20475_3179 * r_2 * h7_m6 - fs_15126300_5909761 * r_2 * h7_m2 - fs_630_169 * r_4 * h5_m2 + fs_110250_20449 * r_6 * h3_m2) + e_4 * (fs_2475_79781 * h11_m6 + fs_4455_2712554 * h11_m2 - fs_180_41327 * r_2 * h9_m6 - fs_8640_537251 * r_2 * h9_m2 - fs_20475_1147619 * r_4 * h7_m6 + fs_15126300_2133423721 * r_4 * h7_m2 + fs_280_48841 * r_6 * h5_m2 - fs_245_40898 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_90[k] = e_0 * (fs_4862025_512 * h3_m1 - fs_2679075_256 * r_2 * h1_m1) + e_1 * (fs_13125_8 * h5_m5 - fs_15125_16 * h5_m1 - fs_60025_8 * r_2 * h3_m1 + fs_54675_16 * r_4 * h1_m1) + e_2 * (-fs_14175_29744 * h7_m5 - fs_231525_81796 * h7_m1 - fs_118125_338 * r_2 * h5_m5 + fs_136125_676 * r_2 * h5_m1 + fs_540225_968 * r_4 * h3_m1 - fs_675_4 * r_6 * h1_m1) + e_3 * (-fs_33075_82654 * h9_m5 + fs_52397415_23639044 * h9_m1 + fs_14175_537251 * r_2 * h7_m5 + fs_926100_5909761 * r_2 * h7_m1 + fs_1050_169 * r_4 * h5_m5 - fs_605_169 * r_4 * h5_m1 - fs_120050_20449 * r_6 * h3_m1 + fs_675_484 * r_8 * h1_m1) + e_4 * (fs_66825_1356277 * h11_m5 + fs_200475_35263202 * h11_m1 + fs_150_41327 * r_2 * h9_m5 - fs_118815_5909761 * r_2 * h9_m1 - fs_14175_193947611 * r_4 * h7_m5 - fs_926100_2133423721 * r_4 * h7_m1 - fs_1400_146523 * r_6 * h5_m5 + fs_2420_439569 * r_6 * h5_m1 + fs_2401_368082 * r_8 * h3_m1 - fs_27_20449 * r_10 * h1_m1) + fs_2679075_1024 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_91[k] = -fs_7875_8 * e_1 * h5_m4 + e_2 * (fs_3444525_59488 * h7_m4 + fs_70875_338 * r_2 * h5_m4) + e_3 * (-fs_138915_82654 * h9_m4 - fs_3444525_1074502 * r_2 * h7_m4 - fs_630_169 * r_4 * h5_m4) + e_4 * (fs_155925_2712554 * h11_m4 + fs_630_41327 * r_2 * h9_m4 + fs_3444525_387895222 * r_4 * h7_m4 + fs_280_48841 * r_6 * h5_m4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_92[k] = e_0 * (fs_3472875_256 * h3_m3 - fs_2083725_256 * h3_m1 - fs_6251175_128 * r_2 * h1_m1) + e_1 * (-fs_15125_16 * h5_m3 + fs_2625_2 * h5_m1 - fs_42875_4 * r_2 * h3_m3 + fs_25725_4 * r_2 * h3_m1 + fs_127575_8 * r_4 * h1_m1) + e_2 * (fs_28923075_654368 * h7_m3 - fs_30305025_654368 * h7_m1 + fs_136125_676 * r_2 * h5_m3 - fs_47250_169 * r_2 * h5_m1 + fs_385875_484 * r_4 * h3_m3 - fs_231525_484 * r_4 * h3_m1 - fs_1575_2 * r_6 * h1_m1) + e_3 * (-fs_2917215_2149004 * h9_m3 + fs_1111320_5909761 * h9_m1 - fs_28923075_11819522 * r_2 * h7_m3 + fs_30305025_11819522 * r_2 * h7_m1 - fs_605_169 * r_4 * h5_m3 + fs_840_169 * r_4 * h5_m1 - fs_171500_20449 * r_6 * h3_m3 + fs_102900_20449 * r_6 * h3_m1 + fs_1575_242 * r_8 * h1_m1) + e_4 * (fs_72765_1356277 * h11_m3 + fs_467775_17631601 * h11_m1 + fs_6615_537251 * r_2 * h9_m3 - fs_10080_5909761 * r_2 * h9_m1 + fs_28923075_4266847442 * r_4 * h7_m3 - fs_30305025_4266847442 * r_4 * h7_m1 + fs_2420_439569 * r_6 * h5_m3 - fs_1120_146523 * r_6 * h5_m1 + fs_1715_184041 * r_8 * h3_m3 - fs_343_61347 * r_8 * h3_m1 - fs_126_20449 * r_10 * h1_m1) + fs_6251175_512 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];

        pc_93[k] = -fs_694575_128 * e_0 * h3_p2 + e_1 * (fs_625_2 * h5_p2 + fs_8575_2 * r_2 * h3_p2) + e_2 * (-fs_91125_327184 * h7_p2 - fs_11250_169 * r_2 * h5_p2 - fs_77175_242 * r_4 * h3_p2) + e_3 * (-fs_231525_537251 * h9_p2 + fs_91125_5909761 * r_2 * h7_p2 + fs_200_169 * r_4 * h5_p2 + fs_68600_20449 * r_6 * h3_p2) + e_4 * (fs_112266_1356277 * h11_p2 + fs_2100_537251 * r_2 * h9_p2 - fs_91125_2133423721 * r_4 * h7_p2 - fs_800_439569 * r_6 * h5_p2 - fs_686_184041 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_94[k] = e_0 * (fs_2083725_256 * h3_p1 + fs_3472875_256 * h3_p3 + fs_6251175_128 * r_2 * h1_p1) + e_1 * (-fs_2625_2 * h5_p1 - fs_15125_16 * h5_p3 - fs_25725_4 * r_2 * h3_p1 - fs_42875_4 * r_2 * h3_p3 - fs_127575_8 * r_4 * h1_p1) + e_2 * (fs_30305025_654368 * h7_p1 + fs_28923075_654368 * h7_p3 + fs_47250_169 * r_2 * h5_p1 + fs_136125_676 * r_2 * h5_p3 + fs_231525_484 * r_4 * h3_p1 + fs_385875_484 * r_4 * h3_p3 + fs_1575_2 * r_6 * h1_p1) + e_3 * (-fs_1111320_5909761 * h9_p1 - fs_2917215_2149004 * h9_p3 - fs_30305025_11819522 * r_2 * h7_p1 - fs_28923075_11819522 * r_2 * h7_p3 - fs_840_169 * r_4 * h5_p1 - fs_605_169 * r_4 * h5_p3 - fs_102900_20449 * r_6 * h3_p1 - fs_171500_20449 * r_6 * h3_p3 - fs_1575_242 * r_8 * h1_p1) + e_4 * (-fs_467775_17631601 * h11_p1 + fs_72765_1356277 * h11_p3 + fs_10080_5909761 * r_2 * h9_p1 + fs_6615_537251 * r_2 * h9_p3 + fs_30305025_4266847442 * r_4 * h7_p1 + fs_28923075_4266847442 * r_4 * h7_p3 + fs_1120_146523 * r_6 * h5_p1 + fs_2420_439569 * r_6 * h5_p3 + fs_343_61347 * r_8 * h3_p1 + fs_1715_184041 * r_8 * h3_p3 + fs_126_20449 * r_10 * h1_p1) - fs_6251175_512 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ph11_0, ph11_p4, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_95[k] = e_0 * (fs_99225_32 * h3_0 + fs_893025_8 * r_2 * h1_0) + e_1 * (fs_625_2 * h5_0 - fs_7875_8 * h5_p4 - fs_2450 * r_2 * h3_0 - fs_36450 * r_4 * h1_0) + e_2 * (-fs_18533025_163592 * h7_0 + fs_3444525_59488 * h7_p4 - fs_11250_169 * r_2 * h5_0 + fs_70875_338 * r_2 * h5_p4 + fs_22050_121 * r_4 * h3_0 + fs_1800 * r_6 * h1_0) + e_3 * (fs_16074450_5909761 * h9_0 - fs_138915_82654 * h9_p4 + fs_37066050_5909761 * r_2 * h7_0 - fs_3444525_1074502 * r_2 * h7_p4 + fs_200_169 * r_4 * h5_0 - fs_630_169 * r_4 * h5_p4 - fs_39200_20449 * r_6 * h3_0 - fs_1800_121 * r_8 * h1_0) + e_4 * (fs_490050_17631601 * h11_0 + fs_155925_2712554 * h11_p4 - fs_145800_5909761 * r_2 * h9_0 + fs_630_41327 * r_2 * h9_p4 - fs_37066050_2133423721 * r_4 * h7_0 + fs_3444525_387895222 * r_4 * h7_p4 - fs_800_439569 * r_6 * h5_0 + fs_280_48841 * r_6 * h5_p4 + fs_392_184041 * r_8 * h3_0 + fs_288_20449 * r_10 * h1_0) - fs_893025_32 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_96[k] = e_0 * (fs_4862025_512 * h3_p1 - fs_2679075_256 * r_2 * h1_p1) + e_1 * (-fs_15125_16 * h5_p1 + fs_13125_8 * h5_p5 - fs_60025_8 * r_2 * h3_p1 + fs_54675_16 * r_4 * h1_p1) + e_2 * (-fs_231525_81796 * h7_p1 - fs_14175_29744 * h7_p5 + fs_136125_676 * r_2 * h5_p1 - fs_118125_338 * r_2 * h5_p5 + fs_540225_968 * r_4 * h3_p1 - fs_675_4 * r_6 * h1_p1) + e_3 * (fs_52397415_23639044 * h9_p1 - fs_33075_82654 * h9_p5 + fs_926100_5909761 * r_2 * h7_p1 + fs_14175_537251 * r_2 * h7_p5 - fs_605_169 * r_4 * h5_p1 + fs_1050_169 * r_4 * h5_p5 - fs_120050_20449 * r_6 * h3_p1 + fs_675_484 * r_8 * h1_p1) + e_4 * (fs_200475_35263202 * h11_p1 + fs_66825_1356277 * h11_p5 - fs_118815_5909761 * r_2 * h9_p1 + fs_150_41327 * r_2 * h9_p5 - fs_926100_2133423721 * r_4 * h7_p1 - fs_14175_193947611 * r_4 * h7_p5 + fs_2420_439569 * r_6 * h5_p1 - fs_1400_146523 * r_6 * h5_p5 + fs_2401_368082 * r_8 * h3_p1 - fs_27_20449 * r_10 * h1_p1) + fs_2679075_1024 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_97[k] = -fs_4465125_512 * e_0 * h3_p2 + e_1 * (-fs_7875_8 * h5_p2 + fs_55125_8 * r_2 * h3_p2) + e_2 * (fs_3781575_81796 * h7_p2 - fs_20475_176 * h7_p6 + fs_70875_338 * r_2 * h5_p2 - fs_496125_968 * r_4 * h3_p2) + e_3 * (fs_952560_537251 * h9_p2 + fs_19845_41327 * h9_p6 - fs_15126300_5909761 * r_2 * h7_p2 + fs_20475_3179 * r_2 * h7_p6 - fs_630_169 * r_4 * h5_p2 + fs_110250_20449 * r_6 * h3_p2) + e_4 * (fs_4455_2712554 * h11_p2 + fs_2475_79781 * h11_p6 - fs_8640_537251 * r_2 * h9_p2 - fs_180_41327 * r_2 * h9_p6 + fs_15126300_2133423721 * r_4 * h7_p2 - fs_20475_1147619 * r_4 * h7_p6 + fs_280_48841 * r_6 * h5_p2 - fs_245_40898 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_98[k] = fs_297675_512 * e_0 * h3_p3 + e_1 * (fs_13125_8 * h5_p3 - fs_3675_8 * r_2 * h3_p3) + e_2 * (fs_9646875_81796 * h7_p3 + fs_55125_2288 * h7_p7 - fs_118125_338 * r_2 * h5_p3 + fs_33075_968 * r_4 * h3_p3) + e_3 * (fs_694575_1074502 * h9_p3 + fs_297675_82654 * h9_p7 - fs_38587500_5909761 * r_2 * h7_p3 - fs_55125_41327 * r_2 * h7_p7 + fs_1050_169 * r_4 * h5_p3 - fs_7350_20449 * r_6 * h3_p3) + e_4 * (fs_693_2712554 * h11_p3 + fs_891_79781 * h11_p7 - fs_3150_537251 * r_2 * h9_p3 - fs_1350_41327 * r_2 * h9_p7 + fs_38587500_2133423721 * r_4 * h7_p3 + fs_55125_14919047 * r_4 * h7_p7 - fs_1400_146523 * r_6 * h5_p3 + fs_49_122694 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_99[k] = -f_945_16 * e_0 * h3_m2 + e_1 * (-fs_39375_16 * h5_m2 + f_105_2 * r_2 * h3_m2) + e_2 * (-fs_3472875_40898 * h7_m2 + fs_354375_676 * r_2 * h5_m2 - f_315_22 * r_4 * h3_m2) + e_3 * (fs_33075_9724 * h9_m8 - fs_297675_1074502 * h9_m2 + fs_27783000_5909761 * r_2 * h7_m2 - fs_1575_169 * r_4 * h5_m2 + f_210_143 * r_6 * h3_m2) + e_4 * (fs_99_4199 * h11_m8 - fs_99_1356277 * h11_m2 - fs_75_2431 * r_2 * h9_m8 + fs_1350_537251 * r_2 * h9_m2 - fs_27783000_2133423721 * r_4 * h7_m2 + fs_700_48841 * r_6 * h5_m2 - f_7_143 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_100[k] = e_0 * (f_945_8 * h3_m1 - fs_2679075_512 * r_2 * h1_m1) + e_1 * (fs_1125_32 * h5_m1 - f_105 * r_2 * h3_m1 + fs_54675_32 * r_4 * h1_m1) + e_2 * (-fs_99225_1144 * h7_m7 - fs_1852200_20449 * h7_m1 - fs_10125_1352 * r_2 * h5_m1 + f_315_11 * r_4 * h3_m1 - fs_675_8 * r_6 * h1_m1) + e_3 * (-fs_6615_165308 * h9_m7 - fs_50068935_47278088 * h9_m1 + fs_198450_41327 * r_2 * h7_m7 + fs_29635200_5909761 * r_2 * h7_m1 + fs_45_338 * r_4 * h5_m1 - f_420_143 * r_6 * h3_m1 + fs_675_968 * r_8 * h1_m1) + e_4 * (fs_3960_79781 * h11_m7 - fs_9900_17631601 * h11_m1 + fs_15_41327 * r_2 * h9_m7 + fs_113535_11819522 * r_2 * h9_m1 - fs_198450_14919047 * r_4 * h7_m7 - fs_29635200_2133423721 * r_4 * h7_m1 - fs_10_48841 * r_6 * h5_m1 + f_14_143 * r_8 * h3_m1 - fs_27_40898 * r_10 * h1_m1) + fs_2679075_2048 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_2, pe_3, pe_4, ph7_m6, ph9_m6, ph11_m6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h7_m6 = ph7_m6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h11_m6 = ph11_m6[k];

        pc_101[k] = fs_14175_286 * e_2 * h7_m6 + e_3 * (-fs_70560_41327 * h9_m6 - fs_113400_41327 * r_2 * h7_m6) + e_4 * (fs_4950_79781 * h11_m6 + fs_640_41327 * r_2 * h9_m6 + fs_113400_14919047 * r_4 * h7_m6);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_102[k] = e_0 * (-fs_2679075_256 * h3_m1 - fs_8037225_128 * r_2 * h1_m1) + e_1 * (-fs_39375_16 * h5_m5 + fs_3375_8 * h5_m1 + fs_33075_4 * r_2 * h3_m1 + fs_164025_8 * r_4 * h1_m1) + e_2 * (fs_3973725_59488 * h7_m5 + fs_7498575_654368 * h7_m1 + fs_354375_676 * r_2 * h5_m5 - fs_30375_338 * r_2 * h5_m1 - fs_297675_484 * r_4 * h3_m1 - fs_2025_2 * r_6 * h1_m1) + e_3 * (-fs_275625_165308 * h9_m5 - fs_25245045_11819522 * h9_m1 - fs_3973725_1074502 * r_2 * h7_m5 - fs_7498575_11819522 * r_2 * h7_m1 - fs_1575_169 * r_4 * h5_m5 + fs_270_169 * r_4 * h5_m1 + fs_132300_20449 * r_6 * h3_m1 + fs_2025_242 * r_8 * h1_m1) + e_4 * (fs_79200_1356277 * h11_m5 - fs_118800_17631601 * h11_m1 + fs_625_41327 * r_2 * h9_m5 + fs_114490_5909761 * r_2 * h9_m1 + fs_3973725_387895222 * r_4 * h7_m5 + fs_7498575_4266847442 * r_4 * h7_m1 + fs_700_48841 * r_6 * h5_m5 - fs_120_48841 * r_6 * h5_m1 - fs_147_20449 * r_8 * h3_m1 - fs_162_20449 * r_10 * h1_m1) + fs_8037225_512 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_103[k] = e_1 * (fs_1125_32 * h5_m4 + fs_3375_8 * h5_m2) + e_2 * (fs_6075_14872 * h7_m4 - fs_10800_169 * h7_m2 - fs_10125_1352 * r_2 * h5_m4 - fs_30375_338 * r_2 * h5_m2) + e_3 * (-fs_108045_330616 * h9_m4 + fs_15435_12716 * h9_m2 - fs_12150_537251 * r_2 * h7_m4 + fs_172800_48841 * r_2 * h7_m2 + fs_45_338 * r_4 * h5_m4 + fs_270_169 * r_4 * h5_m2) + e_4 * (fs_121275_2712554 * h11_m4 + fs_20790_1356277 * h11_m2 + fs_245_82654 * r_2 * h9_m4 - fs_35_3179 * r_2 * h9_m2 + fs_12150_193947611 * r_4 * h7_m4 - fs_172800_17631601 * r_4 * h7_m2 - fs_10_48841 * r_6 * h5_m4 - fs_120_48841 * r_6 * h5_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph11_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p3 = ph11_p3[k];

        pc_104[k] = -fs_2083725_64 * e_0 * h3_p3 + e_1 * (fs_46875_16 * h5_p3 + fs_25725 * r_2 * h3_p3) + e_2 * (-fs_13861125_163592 * h7_p3 - fs_421875_676 * r_2 * h5_p3 - fs_231525_121 * r_4 * h3_p3) + e_3 * (fs_540225_2149004 * h9_p3 + fs_27722250_5909761 * r_2 * h7_p3 + fs_1875_169 * r_4 * h5_p3 + fs_411600_20449 * r_6 * h3_p3) + e_4 * (fs_77616_1356277 * h11_p3 - fs_1225_537251 * r_2 * h9_p3 - fs_27722250_2133423721 * r_4 * h7_p3 - fs_2500_146523 * r_6 * h5_p3 - fs_1372_61347 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_105[k] = e_1 * (-fs_3375_8 * h5_p2 + fs_1125_32 * h5_p4) + e_2 * (fs_10800_169 * h7_p2 + fs_6075_14872 * h7_p4 + fs_30375_338 * r_2 * h5_p2 - fs_10125_1352 * r_2 * h5_p4) + e_3 * (-fs_15435_12716 * h9_p2 - fs_108045_330616 * h9_p4 - fs_172800_48841 * r_2 * h7_p2 - fs_12150_537251 * r_2 * h7_p4 - fs_270_169 * r_4 * h5_p2 + fs_45_338 * r_4 * h5_p4) + e_4 * (-fs_20790_1356277 * h11_p2 + fs_121275_2712554 * h11_p4 + fs_35_3179 * r_2 * h9_p2 + fs_245_82654 * r_2 * h9_p4 + fs_172800_17631601 * r_4 * h7_p2 + fs_12150_193947611 * r_4 * h7_p4 + fs_120_48841 * r_6 * h5_p2 - fs_10_48841 * r_6 * h5_p4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_106[k] = e_0 * (fs_2679075_256 * h3_p1 + fs_8037225_128 * r_2 * h1_p1) + e_1 * (-fs_3375_8 * h5_p1 - fs_39375_16 * h5_p5 - fs_33075_4 * r_2 * h3_p1 - fs_164025_8 * r_4 * h1_p1) + e_2 * (-fs_7498575_654368 * h7_p1 + fs_3973725_59488 * h7_p5 + fs_30375_338 * r_2 * h5_p1 + fs_354375_676 * r_2 * h5_p5 + fs_297675_484 * r_4 * h3_p1 + fs_2025_2 * r_6 * h1_p1) + e_3 * (fs_25245045_11819522 * h9_p1 - fs_275625_165308 * h9_p5 + fs_7498575_11819522 * r_2 * h7_p1 - fs_3973725_1074502 * r_2 * h7_p5 - fs_270_169 * r_4 * h5_p1 - fs_1575_169 * r_4 * h5_p5 - fs_132300_20449 * r_6 * h3_p1 - fs_2025_242 * r_8 * h1_p1) + e_4 * (fs_118800_17631601 * h11_p1 + fs_79200_1356277 * h11_p5 - fs_114490_5909761 * r_2 * h9_p1 + fs_625_41327 * r_2 * h9_p5 - fs_7498575_4266847442 * r_4 * h7_p1 + fs_3973725_387895222 * r_4 * h7_p5 + fs_120_48841 * r_6 * h5_p1 + fs_700_48841 * r_6 * h5_p5 + fs_147_20449 * r_8 * h3_p1 + fs_162_20449 * r_10 * h1_p1) - fs_8037225_512 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ph11_0, ph11_p6, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p6 = ph11_p6[k];

        pc_107[k] = e_0 * (-fs_297675_256 * h3_0 + fs_24111675_256 * r_2 * h1_0) + e_1 * (fs_46875_16 * h5_0 + fs_3675_4 * r_2 * h3_0 - fs_492075_16 * r_4 * h1_0) + e_2 * (-fs_2679075_81796 * h7_0 + fs_14175_286 * h7_p6 - fs_421875_676 * r_2 * h5_0 - fs_33075_484 * r_4 * h3_0 + fs_6075_4 * r_6 * h1_0) + e_3 * (-fs_92907675_23639044 * h9_0 - fs_70560_41327 * h9_p6 + fs_10716300_5909761 * r_2 * h7_0 - fs_113400_41327 * r_2 * h7_p6 + fs_1875_169 * r_4 * h5_0 + fs_14700_20449 * r_6 * h3_0 - fs_6075_484 * r_8 * h1_0) + e_4 * (-fs_81675_17631601 * h11_0 + fs_4950_79781 * h11_p6 + fs_210675_5909761 * r_2 * h9_0 + fs_640_41327 * r_2 * h9_p6 - fs_10716300_2133423721 * r_4 * h7_0 + fs_113400_14919047 * r_4 * h7_p6 - fs_2500_146523 * r_6 * h5_0 - fs_49_61347 * r_8 * h3_0 + fs_243_20449 * r_10 * h1_0) - fs_24111675_1024 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_108[k] = e_0 * (f_945_8 * h3_p1 - fs_2679075_512 * r_2 * h1_p1) + e_1 * (fs_1125_32 * h5_p1 - f_105 * r_2 * h3_p1 + fs_54675_32 * r_4 * h1_p1) + e_2 * (-fs_1852200_20449 * h7_p1 - fs_99225_1144 * h7_p7 - fs_10125_1352 * r_2 * h5_p1 + f_315_11 * r_4 * h3_p1 - fs_675_8 * r_6 * h1_p1) + e_3 * (-fs_50068935_47278088 * h9_p1 - fs_6615_165308 * h9_p7 + fs_29635200_5909761 * r_2 * h7_p1 + fs_198450_41327 * r_2 * h7_p7 + fs_45_338 * r_4 * h5_p1 - f_420_143 * r_6 * h3_p1 + fs_675_968 * r_8 * h1_p1) + e_4 * (-fs_9900_17631601 * h11_p1 + fs_3960_79781 * h11_p7 + fs_113535_11819522 * r_2 * h9_p1 + fs_15_41327 * r_2 * h9_p7 - fs_29635200_2133423721 * r_4 * h7_p1 - fs_198450_14919047 * r_4 * h7_p7 - fs_10_48841 * r_6 * h5_p1 + f_14_143 * r_8 * h3_p1 - fs_27_40898 * r_10 * h1_p1) + fs_2679075_2048 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_109[k] = -f_945_16 * e_0 * h3_p2 + e_1 * (-fs_39375_16 * h5_p2 + f_105_2 * r_2 * h3_p2) + e_2 * (-fs_3472875_40898 * h7_p2 + fs_354375_676 * r_2 * h5_p2 - f_315_22 * r_4 * h3_p2) + e_3 * (-fs_297675_1074502 * h9_p2 + fs_33075_9724 * h9_p8 + fs_27783000_5909761 * r_2 * h7_p2 - fs_1575_169 * r_4 * h5_p2 + f_210_143 * r_6 * h3_p2) + e_4 * (-fs_99_1356277 * h11_p2 + fs_99_4199 * h11_p8 + fs_1350_537251 * r_2 * h9_p2 - fs_75_2431 * r_2 * h9_p8 - fs_27783000_2133423721 * r_4 * h7_p2 + fs_700_48841 * r_6 * h5_p2 - f_7_143 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m8, ph9_m1, ph11_m9, ph11_m8, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

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

        pc_110[k] = e_0 * (fs_2679075_256 * h3_m1 - fs_893025_512 * r_2 * h1_m1) + e_1 * (fs_84375_32 * h5_m1 - fs_33075_4 * r_2 * h3_m1 + fs_18225_32 * r_4 * h1_m1) + e_2 * (fs_1929375_40898 * h7_m1 - fs_759375_1352 * r_2 * h5_m1 + fs_297675_484 * r_4 * h3_m1 - fs_225_8 * r_6 * h1_m1) + e_3 * (fs_19845_9724 * h9_m9 + fs_4465125_47278088 * h9_m1 - fs_15435000_5909761 * r_2 * h7_m1 + fs_3375_338 * r_4 * h5_m1 - fs_132300_20449 * r_6 * h3_m1 + fs_225_968 * r_8 * h1_m1) + e_4 * (fs_198_4199 * h11_m9 + fs_297_17631601 * h11_m1 - fs_45_2431 * r_2 * h9_m9 - fs_10125_11819522 * r_2 * h9_m1 + fs_15435000_2133423721 * r_4 * h7_m1 - fs_750_48841 * r_6 * h5_m1 + fs_147_20449 * r_8 * h3_m1 - fs_9_40898 * r_10 * h1_m1) + fs_893025_2048 * e_5 * h1_m1;

        pc_111[k] = -fs_3969_2431 * e_3 * h9_m8 + e_4 * (fs_297_4199 * h11_m8 + fs_36_2431 * r_2 * h9_m8);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_112[k] = e_0 * (-fs_1488375_256 * h3_m1 - fs_40186125_512 * r_2 * h1_m1) + e_1 * (-fs_12675_32 * h5_m1 + fs_18375_4 * r_2 * h3_m1 + fs_820125_32 * r_4 * h1_m1) + e_2 * (fs_165375_1144 * h7_m7 + fs_70875_968 * h7_m1 + fs_675_8 * r_2 * h5_m1 - fs_165375_484 * r_4 * h3_m1 - fs_10125_8 * r_6 * h1_m1) + e_3 * (-fs_370881_165308 * h9_m7 + fs_59397849_47278088 * h9_m1 - fs_330750_41327 * r_2 * h7_m7 - fs_141750_34969 * r_2 * h7_m1 - fs_3_2 * r_4 * h5_m1 + fs_73500_20449 * r_6 * h3_m1 + fs_10125_968 * r_8 * h1_m1) + e_4 * (fs_5346_79781 * h11_m7 + fs_13365_17631601 * h11_m1 + fs_841_41327 * r_2 * h9_m7 - fs_134689_11819522 * r_2 * h9_m1 + fs_330750_14919047 * r_4 * h7_m7 + fs_141750_12623809 * r_4 * h7_m1 + fs_2_867 * r_6 * h5_m1 - fs_245_61347 * r_8 * h3_m1 - fs_405_40898 * r_10 * h1_m1) + fs_40186125_2048 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_113[k] = -f_945_16 * e_0 * h3_m2 + e_1 * (fs_1575 * h5_m2 + f_105_2 * r_2 * h3_m2) + e_2 * (fs_7875_4576 * h7_m6 - fs_6622875_654368 * h7_m2 - fs_56700_169 * r_2 * h5_m2 - f_315_22 * r_4 * h3_m2) + e_3 * (-fs_35721_82654 * h9_m6 - fs_2223963_1074502 * h9_m2 - fs_7875_82654 * r_2 * h7_m6 + fs_6622875_11819522 * r_2 * h7_m2 + fs_1008_169 * r_4 * h5_m2 + f_210_143 * r_6 * h3_m2) + e_4 * (fs_3960_79781 * h11_m6 - fs_3564_1356277 * h11_m2 + fs_162_41327 * r_2 * h9_m6 + fs_10086_537251 * r_2 * h9_m2 + fs_7875_29838094 * r_4 * h7_m6 - fs_6622875_4266847442 * r_4 * h7_m2 - fs_448_48841 * r_6 * h5_m2 - f_7_143 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_114[k] = fs_2083725_128 * e_0 * h3_m3 + e_1 * (fs_84375_32 * h5_m5 - fs_12675_32 * h5_m3 - fs_25725_2 * r_2 * h3_m3) + e_2 * (-fs_1540125_29744 * h7_m5 - fs_28125_1936 * h7_m3 - fs_759375_1352 * r_2 * h5_m5 + fs_675_8 * r_2 * h5_m3 + fs_231525_242 * r_4 * h3_m3) + e_3 * (fs_46305_330616 * h9_m5 + fs_9529569_4298008 * h9_m3 + fs_1540125_537251 * r_2 * h7_m5 + fs_28125_34969 * r_2 * h7_m3 + fs_3375_338 * r_4 * h5_m5 - fs_3_2 * r_4 * h5_m3 - fs_205800_20449 * r_6 * h3_m3) + e_4 * (fs_41580_1356277 * h11_m5 + fs_9702_1356277 * h11_m3 - fs_105_82654 * r_2 * h9_m5 - fs_21609_1074502 * r_2 * h9_m3 - fs_1540125_193947611 * r_4 * h7_m5 - fs_28125_12623809 * r_4 * h7_m3 - fs_750_48841 * r_6 * h5_m5 + fs_2_867 * r_6 * h5_m3 + fs_686_61347 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p4, ph7_p4, ph9_p4, ph11_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p4 = ph11_p4[k];

        pc_115[k] = fs_1125 * e_1 * h5_p4 + e_2 * (-fs_270000_1859 * h7_p4 - fs_40500_169 * r_2 * h5_p4) + e_3 * (fs_108045_41327 * h9_p4 + fs_4320000_537251 * r_2 * h7_p4 + fs_720_169 * r_4 * h5_p4) + e_4 * (fs_43659_1356277 * h11_p4 - fs_980_41327 * r_2 * h9_p4 - fs_4320000_193947611 * r_4 * h7_p4 - fs_320_48841 * r_6 * h5_p4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_116[k] = -fs_2083725_128 * e_0 * h3_p3 + e_1 * (fs_12675_32 * h5_p3 + fs_84375_32 * h5_p5 + fs_25725_2 * r_2 * h3_p3) + e_2 * (fs_28125_1936 * h7_p3 - fs_1540125_29744 * h7_p5 - fs_675_8 * r_2 * h5_p3 - fs_759375_1352 * r_2 * h5_p5 - fs_231525_242 * r_4 * h3_p3) + e_3 * (-fs_9529569_4298008 * h9_p3 + fs_46305_330616 * h9_p5 - fs_28125_34969 * r_2 * h7_p3 + fs_1540125_537251 * r_2 * h7_p5 + fs_3_2 * r_4 * h5_p3 + fs_3375_338 * r_4 * h5_p5 + fs_205800_20449 * r_6 * h3_p3) + e_4 * (-fs_9702_1356277 * h11_p3 + fs_41580_1356277 * h11_p5 + fs_21609_1074502 * r_2 * h9_p3 - fs_105_82654 * r_2 * h9_p5 + fs_28125_12623809 * r_4 * h7_p3 - fs_1540125_193947611 * r_4 * h7_p5 - fs_2_867 * r_6 * h5_p3 - fs_750_48841 * r_6 * h5_p5 - fs_686_61347 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_117[k] = f_945_16 * e_0 * h3_p2 + e_1 * (-fs_1575 * h5_p2 - f_105_2 * r_2 * h3_p2) + e_2 * (fs_6622875_654368 * h7_p2 + fs_7875_4576 * h7_p6 + fs_56700_169 * r_2 * h5_p2 + f_315_22 * r_4 * h3_p2) + e_3 * (fs_2223963_1074502 * h9_p2 - fs_35721_82654 * h9_p6 - fs_6622875_11819522 * r_2 * h7_p2 - fs_7875_82654 * r_2 * h7_p6 - fs_1008_169 * r_4 * h5_p2 - f_210_143 * r_6 * h3_p2) + e_4 * (fs_3564_1356277 * h11_p2 + fs_3960_79781 * h11_p6 - fs_10086_537251 * r_2 * h9_p2 + fs_162_41327 * r_2 * h9_p6 + fs_6622875_4266847442 * r_4 * h7_p2 + fs_7875_29838094 * r_4 * h7_p6 + fs_448_48841 * r_6 * h5_p2 + f_7_143 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_118[k] = e_0 * (fs_1488375_256 * h3_p1 + fs_40186125_512 * r_2 * h1_p1) + e_1 * (fs_12675_32 * h5_p1 - fs_18375_4 * r_2 * h3_p1 - fs_820125_32 * r_4 * h1_p1) + e_2 * (-fs_70875_968 * h7_p1 + fs_165375_1144 * h7_p7 - fs_675_8 * r_2 * h5_p1 + fs_165375_484 * r_4 * h3_p1 + fs_10125_8 * r_6 * h1_p1) + e_3 * (-fs_59397849_47278088 * h9_p1 - fs_370881_165308 * h9_p7 + fs_141750_34969 * r_2 * h7_p1 - fs_330750_41327 * r_2 * h7_p7 + fs_3_2 * r_4 * h5_p1 - fs_73500_20449 * r_6 * h3_p1 - fs_10125_968 * r_8 * h1_p1) + e_4 * (-fs_13365_17631601 * h11_p1 + fs_5346_79781 * h11_p7 + fs_134689_11819522 * r_2 * h9_p1 + fs_841_41327 * r_2 * h9_p7 - fs_141750_12623809 * r_4 * h7_p1 + fs_330750_14919047 * r_4 * h7_p7 - fs_2_867 * r_6 * h5_p1 + fs_245_61347 * r_8 * h3_p1 + fs_405_40898 * r_10 * h1_p1) - fs_40186125_2048 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph9_p8, ph11_0, ph11_p8, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p8 = ph11_p8[k];

        pc_119[k] = e_0 * (-fs_4465125_256 * h3_0 + fs_4465125_64 * r_2 * h1_0) + e_1 * (fs_1125 * h5_0 + fs_55125_4 * r_2 * h3_0 - fs_91125_4 * r_4 * h1_0) + e_2 * (fs_15931125_81796 * h7_0 - fs_40500_169 * r_2 * h5_0 - fs_496125_484 * r_4 * h3_0 + fs_1125 * r_6 * h1_0) + e_3 * (fs_19845_20449 * h9_0 - fs_3969_2431 * h9_p8 - fs_220500_20449 * r_2 * h7_0 + fs_720_169 * r_4 * h5_0 + fs_220500_20449 * r_6 * h3_0 - fs_1125_121 * r_8 * h1_0) + e_4 * (fs_5445_17631601 * h11_0 + fs_297_4199 * h11_p8 - fs_180_20449 * r_2 * h9_0 + fs_36_2431 * r_2 * h9_p8 + fs_220500_7382089 * r_4 * h7_0 - fs_320_48841 * r_6 * h5_0 - fs_245_20449 * r_8 * h3_0 + fs_180_20449 * r_10 * h1_0) - fs_4465125_256 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_m10, ph11_p1, ph11_p9, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_120[k] = e_0 * (fs_2679075_256 * h3_p1 - fs_893025_512 * r_2 * h1_p1) + e_1 * (fs_84375_32 * h5_p1 - fs_33075_4 * r_2 * h3_p1 + fs_18225_32 * r_4 * h1_p1) + e_2 * (fs_1929375_40898 * h7_p1 - fs_759375_1352 * r_2 * h5_p1 + fs_297675_484 * r_4 * h3_p1 - fs_225_8 * r_6 * h1_p1) + e_3 * (fs_4465125_47278088 * h9_p1 + fs_19845_9724 * h9_p9 - fs_15435000_5909761 * r_2 * h7_p1 + fs_3375_338 * r_4 * h5_p1 - fs_132300_20449 * r_6 * h3_p1 + fs_225_968 * r_8 * h1_p1) + e_4 * (fs_297_17631601 * h11_p1 + fs_198_4199 * h11_p9 - fs_10125_11819522 * r_2 * h9_p1 - fs_45_2431 * r_2 * h9_p9 + fs_15435000_2133423721 * r_4 * h7_p1 - fs_750_48841 * r_6 * h5_p1 + fs_147_20449 * r_8 * h3_p1 - fs_9_40898 * r_10 * h1_p1) + fs_893025_2048 * e_5 * h1_p1;

        pc_121[k] = fs_378_4199 * e_4 * h11_m10;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m1, ph11_m9, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m1 = ph11_m1[k];

        pc_122[k] = -fs_49116375_512 * e_0 * r_2 * h1_m1 + e_1 * (-fs_66825_32 * h5_m1 + fs_1002375_32 * r_4 * h1_m1) + e_2 * (-fs_86625_1352 * h7_m1 + fs_601425_1352 * r_2 * h5_m1 - fs_12375_8 * r_6 * h1_m1) + e_3 * (-fs_3969_884 * h9_m9 - fs_3969_25432 * h9_m1 + fs_173250_48841 * r_2 * h7_m1 - fs_2673_338 * r_4 * h5_m1 + fs_1125_88 * r_8 * h1_m1) + e_4 * (fs_360_4199 * h11_m9 - fs_540_17631601 * h11_m1 + fs_9_221 * r_2 * h9_m9 + fs_9_6358 * r_2 * h9_m1 - fs_173250_17631601 * r_4 * h7_m1 + fs_594_48841 * r_6 * h5_m1 - fs_45_3718 * r_10 * h1_m1) + fs_49116375_2048 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_123[k] = -fs_3274425_256 * e_0 * h3_m2 + e_1 * (fs_5775_16 * h5_m2 + fs_40425_4 * r_2 * h3_m2) + e_2 * (fs_189000_1859 * h7_m2 - fs_51975_676 * r_2 * h5_m2 - fs_33075_44 * r_4 * h3_m2) + e_3 * (-fs_441_884 * h9_m8 + fs_53361_97682 * h9_m2 - fs_3024000_537251 * r_2 * h7_m2 + fs_231_169 * r_4 * h5_m2 + fs_14700_1859 * r_6 * h3_m2) + e_4 * (fs_243_4199 * h11_m8 + fs_243_1356277 * h11_m2 + fs_1_221 * r_2 * h9_m8 - fs_242_48841 * r_2 * h9_m2 + fs_3024000_193947611 * r_4 * h7_m2 - fs_308_146523 * r_6 * h5_m2 - fs_49_5577 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_124[k] = fs_3274425_256 * e_0 * h3_m3 + e_1 * (fs_5775_16 * h5_m3 - fs_40425_4 * r_2 * h3_m3) + e_2 * (-fs_55125_416 * h7_m7 - fs_4921875_59488 * h7_m3 - fs_51975_676 * r_2 * h5_m3 + fs_33075_44 * r_4 * h3_m3) + e_3 * (fs_1323_3757 * h9_m7 - fs_250047_195364 * h9_m3 + fs_55125_7514 * r_2 * h7_m7 + fs_4921875_1074502 * r_2 * h7_m3 + fs_231_169 * r_4 * h5_m3 - fs_14700_1859 * r_6 * h3_m3) + e_4 * (fs_2592_79781 * h11_m7 - fs_1008_1356277 * h11_m3 - fs_12_3757 * r_2 * h9_m7 + fs_567_48841 * r_2 * h9_m3 - fs_55125_2712554 * r_4 * h7_m7 - fs_4921875_387895222 * r_4 * h7_m3 - fs_308_146523 * r_6 * h5_m3 + fs_49_5577 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m4, ph5_p5, ph7_m6, ph7_m4, ph7_p5, ph9_m6, ph9_m4, ph9_p5, ph11_m6, ph11_m4, ph11_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_125[k] = -fs_66825_32 * e_1 * h5_m4 + e_2 * (-fs_1125_13 * h7_m6 + fs_28125_1352 * h7_m4 + fs_601425_1352 * r_2 * h5_m4) + e_3 * (fs_27783_15028 * h9_m6 + fs_64827_30056 * h9_m4 + fs_18000_3757 * r_2 * h7_m6 - fs_56250_48841 * r_2 * h7_m4 - fs_2673_338 * r_4 * h5_m4) + e_4 * (fs_1260_79781 * h11_m6 + fs_6615_2712554 * h11_m4 - fs_63_3757 * r_2 * h9_m6 - fs_147_7514 * r_2 * h9_m4 - fs_18000_1356277 * r_4 * h7_m6 + fs_56250_17631601 * r_4 * h7_m4 + fs_594_48841 * r_6 * h5_m4);

        pc_126[k] = -fs_61875_16 * e_1 * h5_p5 + e_2 * (-fs_16875_1352 * h7_p5 + fs_556875_676 * r_2 * h5_p5) + e_3 * (fs_77175_15028 * h9_p5 + fs_33750_48841 * r_2 * h7_p5 - fs_2475_169 * r_4 * h5_p5) + e_4 * (fs_18144_1356277 * h11_p5 - fs_175_3757 * r_2 * h9_p5 - fs_33750_17631601 * r_4 * h7_p5 + fs_1100_48841 * r_6 * h5_p5);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_127[k] = fs_66825_32 * e_1 * h5_p4 + e_2 * (-fs_28125_1352 * h7_p4 - fs_1125_13 * h7_p6 - fs_601425_1352 * r_2 * h5_p4) + e_3 * (-fs_64827_30056 * h9_p4 + fs_27783_15028 * h9_p6 + fs_56250_48841 * r_2 * h7_p4 + fs_18000_3757 * r_2 * h7_p6 + fs_2673_338 * r_4 * h5_p4) + e_4 * (-fs_6615_2712554 * h11_p4 + fs_1260_79781 * h11_p6 + fs_147_7514 * r_2 * h9_p4 - fs_63_3757 * r_2 * h9_p6 - fs_56250_17631601 * r_4 * h7_p4 - fs_18000_1356277 * r_4 * h7_p6 - fs_594_48841 * r_6 * h5_p4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_128[k] = -fs_3274425_256 * e_0 * h3_p3 + e_1 * (-fs_5775_16 * h5_p3 + fs_40425_4 * r_2 * h3_p3) + e_2 * (fs_4921875_59488 * h7_p3 - fs_55125_416 * h7_p7 + fs_51975_676 * r_2 * h5_p3 - fs_33075_44 * r_4 * h3_p3) + e_3 * (fs_250047_195364 * h9_p3 + fs_1323_3757 * h9_p7 - fs_4921875_1074502 * r_2 * h7_p3 + fs_55125_7514 * r_2 * h7_p7 - fs_231_169 * r_4 * h5_p3 + fs_14700_1859 * r_6 * h3_p3) + e_4 * (fs_1008_1356277 * h11_p3 + fs_2592_79781 * h11_p7 - fs_567_48841 * r_2 * h9_p3 - fs_12_3757 * r_2 * h9_p7 + fs_4921875_387895222 * r_4 * h7_p3 - fs_55125_2712554 * r_4 * h7_p7 + fs_308_146523 * r_6 * h5_p3 - fs_49_5577 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_129[k] = fs_3274425_256 * e_0 * h3_p2 + e_1 * (-fs_5775_16 * h5_p2 - fs_40425_4 * r_2 * h3_p2) + e_2 * (-fs_189000_1859 * h7_p2 + fs_51975_676 * r_2 * h5_p2 + fs_33075_44 * r_4 * h3_p2) + e_3 * (-fs_53361_97682 * h9_p2 - fs_441_884 * h9_p8 + fs_3024000_537251 * r_2 * h7_p2 - fs_231_169 * r_4 * h5_p2 - fs_14700_1859 * r_6 * h3_p2) + e_4 * (-fs_243_1356277 * h11_p2 + fs_243_4199 * h11_p8 + fs_242_48841 * r_2 * h9_p2 + fs_1_221 * r_2 * h9_p8 - fs_3024000_193947611 * r_4 * h7_p2 + fs_308_146523 * r_6 * h5_p2 + fs_49_5577 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_130[k] = fs_49116375_512 * e_0 * r_2 * h1_p1 + e_1 * (fs_66825_32 * h5_p1 - fs_1002375_32 * r_4 * h1_p1) + e_2 * (fs_86625_1352 * h7_p1 - fs_601425_1352 * r_2 * h5_p1 + fs_12375_8 * r_6 * h1_p1) + e_3 * (fs_3969_25432 * h9_p1 - fs_3969_884 * h9_p9 - fs_173250_48841 * r_2 * h7_p1 + fs_2673_338 * r_4 * h5_p1 - fs_1125_88 * r_8 * h1_p1) + e_4 * (fs_540_17631601 * h11_p1 + fs_360_4199 * h11_p9 - fs_9_6358 * r_2 * h9_p1 + fs_9_221 * r_2 * h9_p9 + fs_173250_17631601 * r_4 * h7_p1 - fs_594_48841 * r_6 * h5_p1 + fs_45_3718 * r_10 * h1_p1) - fs_49116375_2048 * e_5 * h1_p1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph11_0, ph11_p10, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p10 = ph11_p10[k];

        pc_131[k] = e_0 * (-fs_9823275_256 * h3_0 + fs_9823275_256 * r_2 * h1_0) + e_1 * (-fs_61875_16 * h5_0 + fs_121275_4 * r_2 * h3_0 - fs_200475_16 * r_4 * h1_0) + e_2 * (-fs_275625_7436 * h7_0 + fs_556875_676 * r_2 * h5_0 - fs_99225_44 * r_4 * h3_0 + fs_2475_4 * r_6 * h1_0) + e_3 * (-fs_99225_2149004 * h9_0 + fs_1102500_537251 * r_2 * h7_0 - fs_2475_169 * r_4 * h5_0 + fs_44100_1859 * r_6 * h3_0 - fs_225_44 * r_8 * h1_0) + e_4 * (-fs_99_17631601 * h11_0 + fs_378_4199 * h11_p10 + fs_225_537251 * r_2 * h9_0 - fs_1102500_193947611 * r_4 * h7_0 + fs_1100_48841 * r_6 * h5_0 - fs_49_1859 * r_8 * h3_0 + fs_9_1859 * r_10 * h1_0) - fs_9823275_1024 * e_5 * h1_0;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m1, ph11_m11, ph11_m1, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m11 = ph11_m11[k];
        const auto h11_m1 = ph11_m1[k];

        pc_132[k] = e_0 * (fs_9823275_512 * h3_m1 - fs_29469825_256 * r_2 * h1_m1) + e_1 * (fs_12375_16 * h5_m1 - fs_121275_8 * r_2 * h3_m1 + fs_601425_16 * r_4 * h1_m1) + e_2 * (fs_118125_29744 * h7_m1 - fs_111375_676 * r_2 * h5_m1 + fs_99225_88 * r_4 * h3_m1 - fs_7425_4 * r_6 * h1_m1) + e_3 * (fs_6615_2149004 * h9_m1 - fs_118125_537251 * r_2 * h7_m1 + fs_495_169 * r_4 * h5_m1 - fs_22050_1859 * r_6 * h3_m1 + fs_675_44 * r_8 * h1_m1) + e_4 * (fs_693_4199 * h11_m11 + fs_9_35263202 * h11_m1 - fs_15_537251 * r_2 * h9_m1 + fs_118125_193947611 * r_4 * h7_m1 - fs_220_48841 * r_6 * h5_m1 + fs_49_3718 * r_8 * h3_m1 - fs_27_1859 * r_10 * h1_m1) + fs_29469825_1024 * e_5 * h1_m1;
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m10, ph11_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m2 = ph11_m2[k];

        pc_133[k] = -fs_9823275_512 * e_0 * h3_m2 + e_1 * (-fs_17325_8 * h5_m2 + fs_121275_8 * r_2 * h3_m2) + e_2 * (-fs_637875_29744 * h7_m2 + fs_155925_338 * r_2 * h5_m2 - fs_99225_88 * r_4 * h3_m2) + e_3 * (-fs_1323_48841 * h9_m2 + fs_637875_537251 * r_2 * h7_m2 - fs_1386_169 * r_4 * h5_m2 + fs_22050_1859 * r_6 * h3_m2) + e_4 * (fs_315_4199 * h11_m10 - fs_9_2712554 * h11_m2 + fs_12_48841 * r_2 * h9_m2 - fs_637875_193947611 * r_4 * h7_m2 + fs_616_48841 * r_6 * h5_m2 - fs_49_3718 * r_8 * h3_m2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m3, ph7_m3, ph9_m9, ph9_m3, ph11_m9, ph11_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m3 = ph11_m3[k];

        pc_134[k] = fs_3274425_512 * e_0 * h3_m3 + e_1 * (fs_5775_2 * h5_m3 - fs_40425_8 * r_2 * h3_m3) + e_2 * (fs_1771875_29744 * h7_m3 - fs_103950_169 * r_2 * h5_m3 + fs_33075_88 * r_4 * h3_m3) + e_3 * (fs_1323_442 * h9_m9 + fs_6174_48841 * h9_m3 - fs_1771875_537251 * r_2 * h7_m3 + fs_1848_169 * r_4 * h5_m3 - fs_7350_1859 * r_6 * h3_m3) + e_4 * (fs_135_4199 * h11_m9 + fs_63_2712554 * h11_m3 - fs_6_221 * r_2 * h9_m9 - fs_56_48841 * r_2 * h9_m3 + fs_1771875_193947611 * r_4 * h7_m3 - fs_2464_146523 * r_6 * h5_m3 + fs_49_11154 * r_8 * h3_m3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_m5, ph5_m4, ph7_m7, ph7_m5, ph7_m4, ph9_m8, ph9_m7, ph9_m5, ph9_m4, ph11_m8, ph11_m7, ph11_m5, ph11_m4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m4 = ph11_m4[k];

        pc_135[k] = -fs_17325_8 * e_1 * h5_m4 + e_2 * (-fs_590625_5408 * h7_m4 + fs_155925_338 * r_2 * h5_m4) + e_3 * (fs_882_221 * h9_m8 - fs_3087_7514 * h9_m4 + fs_590625_97682 * r_2 * h7_m4 - fs_1386_169 * r_4 * h5_m4) + e_4 * (fs_54_4199 * h11_m8 - fs_315_2712554 * h11_m4 - fs_8_221 * r_2 * h9_m8 + fs_14_3757 * r_2 * h9_m4 - fs_590625_35263202 * r_4 * h7_m4 + fs_616_48841 * r_6 * h5_m4);

        pc_136[k] = fs_12375_16 * e_1 * h5_m5 + e_2 * (fs_23625_416 * h7_m7 + fs_759375_5408 * h7_m5 - fs_111375_676 * r_2 * h5_m5) + e_3 * (fs_12348_3757 * h9_m7 + fs_15435_15028 * h9_m5 - fs_23625_7514 * r_2 * h7_m7 - fs_759375_97682 * r_2 * h7_m5 + fs_495_169 * r_4 * h5_m5) + e_4 * (fs_378_79781 * h11_m7 + fs_630_1356277 * h11_m5 - fs_112_3757 * r_2 * h9_m7 - fs_35_3757 * r_2 * h9_m5 + fs_23625_2712554 * r_4 * h7_m7 + fs_759375_35263202 * r_4 * h7_m5 - fs_220_48841 * r_6 * h5_m5);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p5, ph7_p5, ph7_p6, ph7_p7, ph9_p5, ph9_p6, ph9_p7, ph11_p5, ph11_p6, ph11_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p5 = ph5_p5[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p5 = ph11_p5[k];
        const auto h11_p6 = ph11_p6[k];
        const auto h11_p7 = ph11_p7[k];

        pc_137[k] = fs_50625_208 * e_2 * h7_p6 + e_3 * (fs_15435_3757 * h9_p6 - fs_50625_3757 * r_2 * h7_p6) + e_4 * (fs_252_79781 * h11_p6 - fs_140_3757 * r_2 * h9_p6 + fs_50625_1356277 * r_4 * h7_p6);

        pc_138[k] = -fs_12375_16 * e_1 * h5_p5 + e_2 * (-fs_759375_5408 * h7_p5 + fs_23625_416 * h7_p7 + fs_111375_676 * r_2 * h5_p5) + e_3 * (-fs_15435_15028 * h9_p5 + fs_12348_3757 * h9_p7 + fs_759375_97682 * r_2 * h7_p5 - fs_23625_7514 * r_2 * h7_p7 - fs_495_169 * r_4 * h5_p5) + e_4 * (-fs_630_1356277 * h11_p5 + fs_378_79781 * h11_p7 + fs_35_3757 * r_2 * h9_p5 - fs_112_3757 * r_2 * h9_p7 - fs_759375_35263202 * r_4 * h7_p5 + fs_23625_2712554 * r_4 * h7_p7 + fs_220_48841 * r_6 * h5_p5);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph5_p4, ph7_p4, ph9_p4, ph9_p8, ph11_p4, ph11_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p8 = ph11_p8[k];

        pc_139[k] = fs_17325_8 * e_1 * h5_p4 + e_2 * (fs_590625_5408 * h7_p4 - fs_155925_338 * r_2 * h5_p4) + e_3 * (fs_3087_7514 * h9_p4 + fs_882_221 * h9_p8 - fs_590625_97682 * r_2 * h7_p4 + fs_1386_169 * r_4 * h5_p4) + e_4 * (fs_315_2712554 * h11_p4 + fs_54_4199 * h11_p8 - fs_14_3757 * r_2 * h9_p4 - fs_8_221 * r_2 * h9_p8 + fs_590625_35263202 * r_4 * h7_p4 - fs_616_48841 * r_6 * h5_p4);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph9_p9, ph11_p3, ph11_p9, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p9 = ph11_p9[k];

        pc_140[k] = -fs_3274425_512 * e_0 * h3_p3 + e_1 * (-fs_5775_2 * h5_p3 + fs_40425_8 * r_2 * h3_p3) + e_2 * (-fs_1771875_29744 * h7_p3 + fs_103950_169 * r_2 * h5_p3 - fs_33075_88 * r_4 * h3_p3) + e_3 * (-fs_6174_48841 * h9_p3 + fs_1323_442 * h9_p9 + fs_1771875_537251 * r_2 * h7_p3 - fs_1848_169 * r_4 * h5_p3 + fs_7350_1859 * r_6 * h3_p3) + e_4 * (-fs_63_2712554 * h11_p3 + fs_135_4199 * h11_p9 + fs_56_48841 * r_2 * h9_p3 - fs_6_221 * r_2 * h9_p9 - fs_1771875_193947611 * r_4 * h7_p3 + fs_2464_146523 * r_6 * h5_p3 - fs_49_11154 * r_8 * h3_p3);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ph11_p10, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p10 = ph11_p10[k];

        pc_141[k] = fs_9823275_512 * e_0 * h3_p2 + e_1 * (fs_17325_8 * h5_p2 - fs_121275_8 * r_2 * h3_p2) + e_2 * (fs_637875_29744 * h7_p2 - fs_155925_338 * r_2 * h5_p2 + fs_99225_88 * r_4 * h3_p2) + e_3 * (fs_1323_48841 * h9_p2 - fs_637875_537251 * r_2 * h7_p2 + fs_1386_169 * r_4 * h5_p2 - fs_22050_1859 * r_6 * h3_p2) + e_4 * (fs_9_2712554 * h11_p2 + fs_315_4199 * h11_p10 - fs_12_48841 * r_2 * h9_p2 + fs_637875_193947611 * r_4 * h7_p2 - fs_616_48841 * r_6 * h5_p2 + fs_49_3718 * r_8 * h3_p2);
    }

    // NOTE: the rows are formed in 124 loops, as the vectorizer runs out of
    // registers with all 143 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph11_p1, ph11_p11, ab_2 : simd::cache_line_size())
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
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p11 = ph11_p11[k];

        pc_142[k] = e_0 * (-fs_9823275_512 * h3_p1 + fs_29469825_256 * r_2 * h1_p1) + e_1 * (-fs_12375_16 * h5_p1 + fs_121275_8 * r_2 * h3_p1 - fs_601425_16 * r_4 * h1_p1) + e_2 * (-fs_118125_29744 * h7_p1 + fs_111375_676 * r_2 * h5_p1 - fs_99225_88 * r_4 * h3_p1 + fs_7425_4 * r_6 * h1_p1) + e_3 * (-fs_6615_2149004 * h9_p1 + fs_118125_537251 * r_2 * h7_p1 - fs_495_169 * r_4 * h5_p1 + fs_22050_1859 * r_6 * h3_p1 - fs_675_44 * r_8 * h1_p1) + e_4 * (-fs_9_35263202 * h11_p1 + fs_693_4199 * h11_p11 + fs_15_537251 * r_2 * h9_p1 - fs_118125_193947611 * r_4 * h7_p1 + fs_220_48841 * r_6 * h5_p1 - fs_49_3718 * r_8 * h3_p1 + fs_27_1859 * r_10 * h1_p1) - fs_29469825_1024 * e_5 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[143] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142};

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
