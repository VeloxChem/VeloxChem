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



#include "SimdOverlapRecIF.hpp"

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
compute_if_overlap(double                         *values,
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
            false, std::string("SimdOverlapRecIF.compute_if_overlap: Basis functions must be of angular momenta six and three"));
    }

    if (harmonics.size() < 9)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIF.compute_if_overlap: Harmonics must reach angular momentum nine"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIF.compute_if_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 91 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

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

            const auto f_0 = fbase * bexp * bexp * bexp * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * bexp * bexp * bexp * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * bexp * bexp * bexp * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * bexp * bexp * bexp / fexp / fexp / fexp / fexp / fexp / fexp;

            // NOTE: the exponential depends on the pair of primitives alone, so it is
            // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fss = std::exp(-fmu * ab_2[k]);

                pe_0[k] += f_0 * fss;
                pe_1[k] += f_1 * fss;
                pe_2[k] += f_2 * fss;
                pe_3[k] += f_3 * fss;
            }
        }
    }

    // NOTE: the geometry of a term is a solid harmonic of the vector between the
    // atoms times a power of their squared distance.

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_825_32 = std::sqrt(25.78125);
    const auto fs_5775_8 = std::sqrt(721.875);
    const auto fs_405_7436 = std::sqrt(405.0 / 7436.0);
    const auto fs_825_338 = std::sqrt(825.0 / 338.0);
    const auto fs_525_22 = std::sqrt(525.0 / 22.0);
    const auto fs_1_97682 = std::sqrt(1.0 / 97682.0);
    const auto fs_42_221 = std::sqrt(42.0 / 221.0);
    const auto fs_405_537251 = std::sqrt(405.0 / 537251.0);
    const auto fs_11_1014 = std::sqrt(11.0 / 1014.0);
    const auto fs_350_5577 = std::sqrt(350.0 / 5577.0);
    const auto fs_51975_32 = std::sqrt(1624.21875);
    const auto fs_2475_32 = std::sqrt(77.34375);
    const auto fs_135_338 = std::sqrt(135.0 / 338.0);
    const auto fs_2475_338 = std::sqrt(2475.0 / 338.0);
    const auto fs_1_7514 = std::sqrt(1.0 / 7514.0);
    const auto fs_14_221 = std::sqrt(14.0 / 221.0);
    const auto fs_270_48841 = std::sqrt(270.0 / 48841.0);
    const auto fs_11_338 = std::sqrt(11.0 / 338.0);
    const auto fs_243_169 = std::sqrt(243.0 / 169.0);
    const auto fs_189_52 = std::sqrt(189.0 / 52.0);
    const auto fs_7_7514 = std::sqrt(7.0 / 7514.0);
    const auto fs_70_3757 = std::sqrt(70.0 / 3757.0);
    const auto fs_972_48841 = std::sqrt(972.0 / 48841.0);
    const auto fs_189_3757 = std::sqrt(189.0 / 3757.0);
    const auto fs_81_13 = std::sqrt(81.0 / 13.0);
    const auto fs_35_3757 = std::sqrt(35.0 / 3757.0);
    const auto fs_324_3757 = std::sqrt(324.0 / 3757.0);
    const auto fs_825_16 = std::sqrt(51.5625);
    const auto fs_5775_16 = std::sqrt(360.9375);
    const auto fs_3375_14872 = std::sqrt(3375.0 / 14872.0);
    const auto fs_825_169 = std::sqrt(825.0 / 169.0);
    const auto fs_525_44 = std::sqrt(525.0 / 44.0);
    const auto fs_7_97682 = std::sqrt(7.0 / 97682.0);
    const auto fs_28_221 = std::sqrt(28.0 / 221.0);
    const auto fs_3375_1074502 = std::sqrt(3375.0 / 1074502.0);
    const auto fs_11_507 = std::sqrt(11.0 / 507.0);
    const auto fs_175_5577 = std::sqrt(175.0 / 5577.0);
    const auto fs_51975_64 = std::sqrt(812.109375);
    const auto fs_16245_14872 = std::sqrt(16245.0 / 14872.0);
    const auto fs_315_104 = std::sqrt(315.0 / 104.0);
    const auto f_6_221 = 6.0 / 221.0;
    const auto fs_336_3757 = std::sqrt(336.0 / 3757.0);
    const auto fs_16245_1074502 = std::sqrt(16245.0 / 1074502.0);
    const auto fs_315_7514 = std::sqrt(315.0 / 7514.0);
    const auto f_3_2 = 1.5;
    const auto fs_9_104 = std::sqrt(9.0 / 104.0);
    const auto fs_15_3757 = std::sqrt(15.0 / 3757.0);
    const auto f_3_17 = 3.0 / 17.0;
    const auto fs_9_7514 = std::sqrt(9.0 / 7514.0);
    const auto fs_2475_16 = std::sqrt(154.6875);
    const auto fs_1323_338 = std::sqrt(1323.0 / 338.0);
    const auto fs_2475_169 = std::sqrt(2475.0 / 169.0);
    const auto fs_112_3757 = std::sqrt(112.0 / 3757.0);
    const auto fs_2646_48841 = std::sqrt(2646.0 / 48841.0);
    const auto fs_11_169 = std::sqrt(11.0 / 169.0);
    const auto fs_525_8 = std::sqrt(65.625);
    const auto fs_2625_16 = std::sqrt(164.0625);
    const auto fs_91125_163592 = std::sqrt(91125.0 / 163592.0);
    const auto fs_945_1144 = std::sqrt(945.0 / 1144.0);
    const auto fs_1050_169 = std::sqrt(1050.0 / 169.0);
    const auto fs_2625_484 = std::sqrt(2625.0 / 484.0);
    const auto fs_14_48841 = std::sqrt(14.0 / 48841.0);
    const auto fs_308_3757 = std::sqrt(308.0 / 3757.0);
    const auto fs_91125_11819522 = std::sqrt(91125.0 / 11819522.0);
    const auto fs_945_82654 = std::sqrt(945.0 / 82654.0);
    const auto fs_14_507 = std::sqrt(14.0 / 507.0);
    const auto fs_875_61347 = std::sqrt(875.0 / 61347.0);
    const auto fs_23625_64 = std::sqrt(369.140625);
    const auto f_15_4 = 3.75;
    const auto fs_1575_4 = std::sqrt(393.75);
    const auto fs_36000_20449 = std::sqrt(36000.0 / 20449.0);
    const auto fs_360_143 = std::sqrt(360.0 / 143.0);
    const auto f_15_13 = 15.0 / 13.0;
    const auto fs_1575_121 = std::sqrt(1575.0 / 121.0);
    const auto fs_231_97682 = std::sqrt(231.0 / 97682.0);
    const auto fs_693_7514 = std::sqrt(693.0 / 7514.0);
    const auto fs_144000_5909761 = std::sqrt(144000.0 / 5909761.0);
    const auto fs_1440_41327 = std::sqrt(1440.0 / 41327.0);
    const auto f_1_13 = 1.0 / 13.0;
    const auto fs_700_20449 = std::sqrt(700.0 / 20449.0);
    const auto fs_14175_16 = std::sqrt(885.9375);
    const auto fs_375_16 = std::sqrt(23.4375);
    const auto fs_675_16 = std::sqrt(42.1875);
    const auto fs_328329_163592 = std::sqrt(328329.0 / 163592.0);
    const auto fs_7569_14872 = std::sqrt(7569.0 / 14872.0);
    const auto fs_375_169 = std::sqrt(375.0 / 169.0);
    const auto fs_675_169 = std::sqrt(675.0 / 169.0);
    const auto fs_495_48841 = std::sqrt(495.0 / 48841.0);
    const auto fs_231_3757 = std::sqrt(231.0 / 3757.0);
    const auto fs_328329_11819522 = std::sqrt(328329.0 / 11819522.0);
    const auto fs_7569_1074502 = std::sqrt(7569.0 / 1074502.0);
    const auto fs_5_507 = std::sqrt(5.0 / 507.0);
    const auto fs_3_169 = std::sqrt(3.0 / 169.0);
    const auto fs_1125_16 = std::sqrt(70.3125);
    const auto fs_1728_1859 = std::sqrt(1728.0 / 1859.0);
    const auto fs_1125_169 = std::sqrt(1125.0 / 169.0);
    const auto fs_220_3757 = std::sqrt(220.0 / 3757.0);
    const auto fs_6912_537251 = std::sqrt(6912.0 / 537251.0);
    const auto fs_5_169 = std::sqrt(5.0 / 169.0);
    const auto fs_525_4 = std::sqrt(131.25);
    const auto fs_42525_20449 = std::sqrt(42525.0 / 20449.0);
    const auto fs_2025_1144 = std::sqrt(2025.0 / 1144.0);
    const auto fs_2100_169 = std::sqrt(2100.0 / 169.0);
    const auto fs_525_121 = std::sqrt(525.0 / 121.0);
    const auto fs_84_48841 = std::sqrt(84.0 / 48841.0);
    const auto fs_385_7514 = std::sqrt(385.0 / 7514.0);
    const auto fs_170100_5909761 = std::sqrt(170100.0 / 5909761.0);
    const auto fs_2025_82654 = std::sqrt(2025.0 / 82654.0);
    const auto fs_28_507 = std::sqrt(28.0 / 507.0);
    const auto fs_700_61347 = std::sqrt(700.0 / 61347.0);
    const auto fs_4725_16 = std::sqrt(295.3125);
    const auto fs_2025_968 = std::sqrt(2025.0 / 968.0);
    const auto fs_16875_14872 = std::sqrt(16875.0 / 14872.0);
    const auto fs_4725_484 = std::sqrt(4725.0 / 484.0);
    const auto fs_280_48841 = std::sqrt(280.0 / 48841.0);
    const auto fs_2025_69938 = std::sqrt(2025.0 / 69938.0);
    const auto fs_16875_1074502 = std::sqrt(16875.0 / 1074502.0);
    const auto fs_525_20449 = std::sqrt(525.0 / 20449.0);
    const auto fs_42525_64 = std::sqrt(664.453125);
    const auto f_15_2 = 7.5;
    const auto fs_184815_163592 = std::sqrt(184815.0 / 163592.0);
    const auto fs_10935_7436 = std::sqrt(10935.0 / 7436.0);
    const auto f_30_13 = 30.0 / 13.0;
    const auto fs_1925_97682 = std::sqrt(1925.0 / 97682.0);
    const auto fs_275_3757 = std::sqrt(275.0 / 3757.0);
    const auto fs_184815_11819522 = std::sqrt(184815.0 / 11819522.0);
    const auto fs_10935_537251 = std::sqrt(10935.0 / 537251.0);
    const auto f_2_13 = 2.0 / 13.0;
    const auto fs_75_16 = std::sqrt(4.6875);
    const auto fs_405_40898 = std::sqrt(405.0 / 40898.0);
    const auto fs_75_169 = std::sqrt(75.0 / 169.0);
    const auto fs_4400_48841 = std::sqrt(4400.0 / 48841.0);
    const auto fs_810_5909761 = std::sqrt(810.0 / 5909761.0);
    const auto fs_1_507 = std::sqrt(1.0 / 507.0);
    const auto fs_875_16 = std::sqrt(54.6875);
    const auto fs_75_32 = std::sqrt(2.34375);
    const auto fs_175_8 = std::sqrt(21.875);
    const auto fs_33075_20449 = std::sqrt(33075.0 / 20449.0);
    const auto fs_18225_7436 = std::sqrt(18225.0 / 7436.0);
    const auto fs_875_169 = std::sqrt(875.0 / 169.0);
    const auto fs_75_338 = std::sqrt(75.0 / 338.0);
    const auto fs_175_242 = std::sqrt(175.0 / 242.0);
    const auto fs_105_48841 = std::sqrt(105.0 / 48841.0);
    const auto fs_231_7514 = std::sqrt(231.0 / 7514.0);
    const auto fs_132300_5909761 = std::sqrt(132300.0 / 5909761.0);
    const auto fs_18225_537251 = std::sqrt(18225.0 / 537251.0);
    const auto fs_35_1521 = std::sqrt(35.0 / 1521.0);
    const auto fs_1_1014 = std::sqrt(1.0 / 1014.0);
    const auto fs_350_184041 = std::sqrt(350.0 / 184041.0);
    const auto fs_1575_32 = std::sqrt(49.21875);
    const auto fs_1125_32 = std::sqrt(35.15625);
    const auto fs_350 = std::sqrt(350.0);
    const auto fs_78750_20449 = std::sqrt(78750.0 / 20449.0);
    const auto fs_675_3718 = std::sqrt(675.0 / 3718.0);
    const auto fs_350_169 = std::sqrt(350.0 / 169.0);
    const auto fs_1125_338 = std::sqrt(1125.0 / 338.0);
    const auto fs_1400_121 = std::sqrt(1400.0 / 121.0);
    const auto fs_1134_48841 = std::sqrt(1134.0 / 48841.0);
    const auto fs_495_7514 = std::sqrt(495.0 / 7514.0);
    const auto fs_315000_5909761 = std::sqrt(315000.0 / 5909761.0);
    const auto fs_1350_537251 = std::sqrt(1350.0 / 537251.0);
    const auto fs_14_1521 = std::sqrt(14.0 / 1521.0);
    const auto fs_5_338 = std::sqrt(5.0 / 338.0);
    const auto fs_5600_184041 = std::sqrt(5600.0 / 184041.0);
    const auto fs_1575_2 = std::sqrt(787.5);
    const auto fs_525_16 = std::sqrt(32.8125);
    const auto fs_1225_32 = std::sqrt(38.28125);
    const auto fs_2625_8 = std::sqrt(328.125);
    const auto fs_23805_81796 = std::sqrt(23805.0 / 81796.0);
    const auto fs_34560_20449 = std::sqrt(34560.0 / 20449.0);
    const auto fs_525_169 = std::sqrt(525.0 / 169.0);
    const auto fs_1225_338 = std::sqrt(1225.0 / 338.0);
    const auto fs_2625_242 = std::sqrt(2625.0 / 242.0);
    const auto fs_1575_48841 = std::sqrt(1575.0 / 48841.0);
    const auto fs_7425_97682 = std::sqrt(7425.0 / 97682.0);
    const auto fs_23805_5909761 = std::sqrt(23805.0 / 5909761.0);
    const auto fs_138240_5909761 = std::sqrt(138240.0 / 5909761.0);
    const auto fs_7_507 = std::sqrt(7.0 / 507.0);
    const auto fs_49_3042 = std::sqrt(49.0 / 3042.0);
    const auto fs_1750_61347 = std::sqrt(1750.0 / 61347.0);
    const auto fs_23625_32 = std::sqrt(738.28125);
    const auto fs_25_2 = std::sqrt(12.5);
    const auto fs_19845_20449 = std::sqrt(19845.0 / 20449.0);
    const auto fs_200_169 = std::sqrt(200.0 / 169.0);
    const auto fs_5775_48841 = std::sqrt(5775.0 / 48841.0);
    const auto fs_79380_5909761 = std::sqrt(79380.0 / 5909761.0);
    const auto fs_8_1521 = std::sqrt(8.0 / 1521.0);
    const auto fs_75_8 = std::sqrt(9.375);
    const auto fs_175_32 = std::sqrt(5.46875);
    const auto fs_178605_81796 = std::sqrt(178605.0 / 81796.0);
    const auto fs_10125_3718 = std::sqrt(10125.0 / 3718.0);
    const auto fs_150_169 = std::sqrt(150.0 / 169.0);
    const auto fs_175_968 = std::sqrt(175.0 / 968.0);
    const auto fs_231_48841 = std::sqrt(231.0 / 48841.0);
    const auto fs_66_3757 = std::sqrt(66.0 / 3757.0);
    const auto fs_178605_5909761 = std::sqrt(178605.0 / 5909761.0);
    const auto fs_20250_537251 = std::sqrt(20250.0 / 537251.0);
    const auto fs_2_507 = std::sqrt(2.0 / 507.0);
    const auto fs_175_368082 = std::sqrt(175.0 / 368082.0);
    const auto fs_1575_128 = std::sqrt(12.3046875);
    const auto fs_50 = std::sqrt(50.0);
    const auto fs_2625_32 = std::sqrt(82.03125);
    const auto fs_108045_81796 = std::sqrt(108045.0 / 81796.0);
    const auto fs_3375_81796 = std::sqrt(3375.0 / 81796.0);
    const auto fs_800_169 = std::sqrt(800.0 / 169.0);
    const auto fs_2625_968 = std::sqrt(2625.0 / 968.0);
    const auto fs_1008_48841 = std::sqrt(1008.0 / 48841.0);
    const auto fs_2376_48841 = std::sqrt(2376.0 / 48841.0);
    const auto fs_108045_5909761 = std::sqrt(108045.0 / 5909761.0);
    const auto fs_3375_5909761 = std::sqrt(3375.0 / 5909761.0);
    const auto fs_32_1521 = std::sqrt(32.0 / 1521.0);
    const auto fs_875_122694 = std::sqrt(875.0 / 122694.0);
    const auto fs_23625_128 = std::sqrt(184.5703125);
    const auto fs_375_32 = std::sqrt(11.71875);
    const auto fs_4375_8 = std::sqrt(546.875);
    const auto fs_126_20449 = std::sqrt(126.0 / 20449.0);
    const auto fs_93987_81796 = std::sqrt(93987.0 / 81796.0);
    const auto fs_375_338 = std::sqrt(375.0 / 338.0);
    const auto fs_4375_242 = std::sqrt(4375.0 / 242.0);
    const auto fs_4536_48841 = std::sqrt(4536.0 / 48841.0);
    const auto fs_3465_48841 = std::sqrt(3465.0 / 48841.0);
    const auto fs_504_5909761 = std::sqrt(504.0 / 5909761.0);
    const auto fs_93987_5909761 = std::sqrt(93987.0 / 5909761.0);
    const auto fs_5_1014 = std::sqrt(5.0 / 1014.0);
    const auto fs_8750_184041 = std::sqrt(8750.0 / 184041.0);
    const auto fs_39375_32 = std::sqrt(1230.46875);
    const auto fs_49923_20449 = std::sqrt(49923.0 / 20449.0);
    const auto fs_6720_48841 = std::sqrt(6720.0 / 48841.0);
    const auto fs_199692_5909761 = std::sqrt(199692.0 / 5909761.0);
    const auto fs_175_4 = std::sqrt(43.75);
    const auto f_5_4 = 1.25;
    const auto fs_212625_40898 = std::sqrt(212625.0 / 40898.0);
    const auto fs_700_169 = std::sqrt(700.0 / 169.0);
    const auto f_5_22 = 5.0 / 22.0;
    const auto fs_924_48841 = std::sqrt(924.0 / 48841.0);
    const auto fs_425250_5909761 = std::sqrt(425250.0 / 5909761.0);
    const auto fs_28_1521 = std::sqrt(28.0 / 1521.0);
    const auto f_5_429 = 5.0 / 429.0;
    const auto f_15_8 = 1.875;
    const auto fs_1575_16 = std::sqrt(98.4375);
    const auto fs_22680_20449 = std::sqrt(22680.0 / 20449.0);
    const auto fs_1575_169 = std::sqrt(1575.0 / 169.0);
    const auto f_15_11 = 15.0 / 11.0;
    const auto fs_3234_48841 = std::sqrt(3234.0 / 48841.0);
    const auto fs_90720_5909761 = std::sqrt(90720.0 / 5909761.0);
    const auto fs_7_169 = std::sqrt(7.0 / 169.0);
    const auto f_10_143 = 10.0 / 143.0;
    const auto f_45_4 = 11.25;
    const auto f_75_4 = 18.75;
    const auto fs_189_242 = std::sqrt(189.0 / 242.0);
    const auto f_75_22 = 75.0 / 22.0;
    const auto fs_5880_48841 = std::sqrt(5880.0 / 48841.0);
    const auto fs_378_34969 = std::sqrt(378.0 / 34969.0);
    const auto f_25_143 = 25.0 / 143.0;
    const auto f_225_8 = 28.125;
    const auto f_35_4 = 8.75;
    const auto f_25 = 25.0;
    const auto f_252_143 = 252.0 / 143.0;
    const auto f_35_13 = 35.0 / 13.0;
    const auto f_50_11 = 50.0 / 11.0;
    const auto f_84_221 = 84.0 / 221.0;
    const auto f_504_2431 = 504.0 / 2431.0;
    const auto f_7_39 = 7.0 / 39.0;
    const auto f_100_429 = 100.0 / 429.0;
    const auto f_75_2 = 37.5;

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph9_p3, ph9_p4, ph9_p8, ph9_p9, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h9_p9 = ph9_p9[k];

        pc_0[k] = e_0 * (-fs_825_32 * h5_p3 + fs_5775_8 * r_2 * h3_p3) + e_1 * (-fs_405_7436 * h7_p3 + fs_825_338 * r_2 * h5_p3 - fs_525_22 * r_4 * h3_p3) + e_2 * (-fs_1_97682 * h9_p3 - fs_42_221 * h9_p9 + fs_405_537251 * r_2 * h7_p3 - fs_11_1014 * r_4 * h5_p3 + fs_350_5577 * r_6 * h3_p3) - fs_51975_32 * e_3 * h3_p3;

        pc_1[k] = fs_2475_32 * e_0 * h5_p4 + e_1 * (fs_135_338 * h7_p4 - fs_2475_338 * r_2 * h5_p4) + e_2 * (fs_1_7514 * h9_p4 - fs_14_221 * h9_p8 - fs_270_48841 * r_2 * h7_p4 + fs_11_338 * r_4 * h5_p4);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_m5, ph5_p5, ph7_m7, ph7_m6, ph7_m5, ph7_p5, ph7_p7, ph9_m7, ph9_m6, ph9_m5, ph9_p5, ph9_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_2[k] = -fs_2475_32 * e_0 * h5_p5 + e_1 * (-fs_243_169 * h7_p5 - fs_189_52 * h7_p7 + fs_2475_338 * r_2 * h5_p5) + e_2 * (-fs_7_7514 * h9_p5 - fs_70_3757 * h9_p7 + fs_972_48841 * r_2 * h7_p5 + fs_189_3757 * r_2 * h7_p7 - fs_11_338 * r_4 * h5_p5);

        pc_3[k] = fs_81_13 * e_1 * h7_m6 + e_2 * (fs_35_3757 * h9_m6 - fs_324_3757 * r_2 * h7_m6);

        pc_4[k] = -fs_2475_32 * e_0 * h5_m5 + e_1 * (fs_189_52 * h7_m7 - fs_243_169 * h7_m5 + fs_2475_338 * r_2 * h5_m5) + e_2 * (fs_70_3757 * h9_m7 - fs_7_7514 * h9_m5 - fs_189_3757 * r_2 * h7_m7 + fs_972_48841 * r_2 * h7_m5 - fs_11_338 * r_4 * h5_m5);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m4, ph5_m3, ph7_m4, ph7_m3, ph9_m9, ph9_m8, ph9_m4, ph9_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_5[k] = fs_2475_32 * e_0 * h5_m4 + e_1 * (fs_135_338 * h7_m4 - fs_2475_338 * r_2 * h5_m4) + e_2 * (fs_14_221 * h9_m8 + fs_1_7514 * h9_m4 - fs_270_48841 * r_2 * h7_m4 + fs_11_338 * r_4 * h5_m4);

        pc_6[k] = e_0 * (-fs_825_32 * h5_m3 + fs_5775_8 * r_2 * h3_m3) + e_1 * (-fs_405_7436 * h7_m3 + fs_825_338 * r_2 * h5_m3 - fs_525_22 * r_4 * h3_m3) + e_2 * (fs_42_221 * h9_m9 - fs_1_97682 * h9_m3 + fs_405_537251 * r_2 * h7_m3 - fs_11_1014 * r_4 * h5_m3 + fs_350_5577 * r_6 * h3_m3) - fs_51975_32 * e_3 * h3_m3;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph7_p7, ph9_p2, ph9_p3, ph9_p7, ph9_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_7[k] = e_0 * (-fs_825_16 * h5_p2 + fs_5775_16 * r_2 * h3_p2) + e_1 * (-fs_3375_14872 * h7_p2 + fs_825_169 * r_2 * h5_p2 - fs_525_44 * r_4 * h3_p2) + e_2 * (-fs_7_97682 * h9_p2 - fs_28_221 * h9_p8 + fs_3375_1074502 * r_2 * h7_p2 - fs_11_507 * r_4 * h5_p2 + fs_175_5577 * r_6 * h3_p2) - fs_51975_64 * e_3 * h3_p2;

        pc_8[k] = e_0 * (fs_825_16 * h5_p3 + fs_5775_16 * r_2 * h3_p3) + e_1 * (fs_16245_14872 * h7_p3 + fs_315_104 * h7_p7 - fs_825_169 * r_2 * h5_p3 - fs_525_44 * r_4 * h3_p3) + e_2 * (f_6_221 * h9_p3 - fs_336_3757 * h9_p7 - fs_16245_1074502 * r_2 * h7_p3 - fs_315_7514 * r_2 * h7_p7 + fs_11_507 * r_4 * h5_p3 + fs_175_5577 * r_6 * h3_p3) - fs_51975_64 * e_3 * h3_p3;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_m5, ph7_m6, ph7_m5, ph7_m4, ph7_p4, ph7_p6, ph9_m6, ph9_m5, ph9_m4, ph9_p4, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_9[k] = e_1 * (-f_3_2 * h7_p4 - fs_9_104 * h7_p6) + e_2 * (-fs_15_3757 * h9_p4 - fs_315_7514 * h9_p6 + f_3_17 * r_2 * h7_p4 + fs_9_7514 * r_2 * h7_p6);

        pc_10[k] = -fs_2475_16 * e_0 * h5_m5 + e_1 * (fs_1323_338 * h7_m5 + fs_2475_169 * r_2 * h5_m5) + e_2 * (fs_112_3757 * h9_m5 - fs_2646_48841 * r_2 * h7_m5 - fs_11_169 * r_4 * h5_m5);

        pc_11[k] = e_1 * (fs_9_104 * h7_m6 - f_3_2 * h7_m4) + e_2 * (fs_315_7514 * h9_m6 - fs_15_3757 * h9_m4 - fs_9_7514 * r_2 * h7_m6 + f_3_17 * r_2 * h7_m4);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m7, ph7_m3, ph7_m2, ph9_m8, ph9_m7, ph9_m3, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_12[k] = e_0 * (fs_825_16 * h5_m3 + fs_5775_16 * r_2 * h3_m3) + e_1 * (-fs_315_104 * h7_m7 + fs_16245_14872 * h7_m3 - fs_825_169 * r_2 * h5_m3 - fs_525_44 * r_4 * h3_m3) + e_2 * (fs_336_3757 * h9_m7 + f_6_221 * h9_m3 + fs_315_7514 * r_2 * h7_m7 - fs_16245_1074502 * r_2 * h7_m3 + fs_11_507 * r_4 * h5_m3 + fs_175_5577 * r_6 * h3_m3) - fs_51975_64 * e_3 * h3_m3;

        pc_13[k] = e_0 * (-fs_825_16 * h5_m2 + fs_5775_16 * r_2 * h3_m2) + e_1 * (-fs_3375_14872 * h7_m2 + fs_825_169 * r_2 * h5_m2 - fs_525_44 * r_4 * h3_m2) + e_2 * (fs_28_221 * h9_m8 - fs_7_97682 * h9_m2 + fs_3375_1074502 * r_2 * h7_m2 - fs_11_507 * r_4 * h5_m2 + fs_175_5577 * r_6 * h3_m2) - fs_51975_64 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph7_p6, ph7_p7, ph9_p1, ph9_p2, ph9_p6, ph9_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_14[k] = e_0 * (-fs_525_8 * h5_p1 + fs_2625_16 * r_2 * h3_p1) + e_1 * (-fs_91125_163592 * h7_p1 - fs_945_1144 * h7_p7 + fs_1050_169 * r_2 * h5_p1 - fs_2625_484 * r_4 * h3_p1) + e_2 * (-fs_14_48841 * h9_p1 - fs_308_3757 * h9_p7 + fs_91125_11819522 * r_2 * h7_p1 + fs_945_82654 * r_2 * h7_p7 - fs_14_507 * r_4 * h5_p1 + fs_875_61347 * r_6 * h3_p1) - fs_23625_64 * e_3 * h3_p1;

        pc_15[k] = e_0 * (f_15_4 * h5_p2 + fs_1575_4 * r_2 * h3_p2) + e_1 * (fs_36000_20449 * h7_p2 + fs_360_143 * h7_p6 - f_15_13 * r_2 * h5_p2 - fs_1575_121 * r_4 * h3_p2) + e_2 * (fs_231_97682 * h9_p2 - fs_693_7514 * h9_p6 - fs_144000_5909761 * r_2 * h7_p2 - fs_1440_41327 * r_2 * h7_p6 + f_1_13 * r_4 * h5_p2 + fs_700_20449 * r_6 * h3_p2) - fs_14175_16 * e_3 * h3_p2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_m4, ph5_p3, ph5_p5, ph7_m4, ph7_p3, ph7_p5, ph9_m4, ph9_p3, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_16[k] = e_0 * (fs_375_16 * h5_p3 - fs_675_16 * h5_p5 + fs_2625_16 * r_2 * h3_p3) + e_1 * (-fs_328329_163592 * h7_p3 + fs_7569_14872 * h7_p5 - fs_375_169 * r_2 * h5_p3 + fs_675_169 * r_2 * h5_p5 - fs_2625_484 * r_4 * h3_p3) + e_2 * (-fs_495_48841 * h9_p3 - fs_231_3757 * h9_p5 + fs_328329_11819522 * r_2 * h7_p3 - fs_7569_1074502 * r_2 * h7_p5 + fs_5_507 * r_4 * h5_p3 - fs_3_169 * r_4 * h5_p5 + fs_875_61347 * r_6 * h3_p3) - fs_23625_64 * e_3 * h3_p3;

        pc_17[k] = -fs_1125_16 * e_0 * h5_m4 + e_1 * (fs_1728_1859 * h7_m4 + fs_1125_169 * r_2 * h5_m4) + e_2 * (fs_220_3757 * h9_m4 - fs_6912_537251 * r_2 * h7_m4 - fs_5_169 * r_4 * h5_m4);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m5, ph5_m3, ph5_m2, ph7_m6, ph7_m5, ph7_m3, ph7_m2, ph9_m6, ph9_m5, ph9_m3, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_18[k] = e_0 * (fs_675_16 * h5_m5 + fs_375_16 * h5_m3 + fs_2625_16 * r_2 * h3_m3) + e_1 * (-fs_7569_14872 * h7_m5 - fs_328329_163592 * h7_m3 - fs_675_169 * r_2 * h5_m5 - fs_375_169 * r_2 * h5_m3 - fs_2625_484 * r_4 * h3_m3) + e_2 * (fs_231_3757 * h9_m5 - fs_495_48841 * h9_m3 + fs_7569_1074502 * r_2 * h7_m5 + fs_328329_11819522 * r_2 * h7_m3 + fs_3_169 * r_4 * h5_m5 + fs_5_507 * r_4 * h5_m3 + fs_875_61347 * r_6 * h3_m3) - fs_23625_64 * e_3 * h3_m3;

        pc_19[k] = e_0 * (f_15_4 * h5_m2 + fs_1575_4 * r_2 * h3_m2) + e_1 * (-fs_360_143 * h7_m6 + fs_36000_20449 * h7_m2 - f_15_13 * r_2 * h5_m2 - fs_1575_121 * r_4 * h3_m2) + e_2 * (fs_693_7514 * h9_m6 + fs_231_97682 * h9_m2 + fs_1440_41327 * r_2 * h7_m6 - fs_144000_5909761 * r_2 * h7_m2 + f_1_13 * r_4 * h5_m2 + fs_700_20449 * r_6 * h3_m2) - fs_14175_16 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph3_0, ph5_m1, ph5_0, ph7_m7, ph7_m1, ph7_0, ph7_p6, ph9_m7, ph9_m1, ph9_0, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];

        pc_20[k] = e_0 * (-fs_525_8 * h5_m1 + fs_2625_16 * r_2 * h3_m1) + e_1 * (fs_945_1144 * h7_m7 - fs_91125_163592 * h7_m1 + fs_1050_169 * r_2 * h5_m1 - fs_2625_484 * r_4 * h3_m1) + e_2 * (fs_308_3757 * h9_m7 - fs_14_48841 * h9_m1 - fs_945_82654 * r_2 * h7_m7 + fs_91125_11819522 * r_2 * h7_m1 - fs_14_507 * r_4 * h5_m1 + fs_875_61347 * r_6 * h3_m1) - fs_23625_64 * e_3 * h3_m1;

        pc_21[k] = e_0 * (-fs_525_4 * h5_0 + fs_525_4 * r_2 * h3_0) + e_1 * (-fs_42525_20449 * h7_0 - fs_2025_1144 * h7_p6 + fs_2100_169 * r_2 * h5_0 - fs_525_121 * r_4 * h3_0) + e_2 * (-fs_84_48841 * h9_0 - fs_385_7514 * h9_p6 + fs_170100_5909761 * r_2 * h7_0 + fs_2025_82654 * r_2 * h7_p6 - fs_28_507 * r_4 * h5_0 + fs_700_61347 * r_6 * h3_0) - fs_4725_16 * e_3 * h3_0;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p2, ph5_p2, ph5_p4, ph5_p5, ph7_p1, ph7_p2, ph7_p4, ph7_p5, ph9_p1, ph9_p2, ph9_p4, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];

        pc_22[k] = e_0 * (f_15_4 * h5_p5 + fs_4725_16 * r_2 * h3_p1) + e_1 * (fs_2025_968 * h7_p1 + fs_16875_14872 * h7_p5 - f_15_13 * r_2 * h5_p5 - fs_4725_484 * r_4 * h3_p1) + e_2 * (fs_280_48841 * h9_p1 - fs_308_3757 * h9_p5 - fs_2025_69938 * r_2 * h7_p1 - fs_16875_1074502 * r_2 * h7_p5 + f_1_13 * r_4 * h5_p5 + fs_525_20449 * r_6 * h3_p1) - fs_42525_64 * e_3 * h3_p1;

        pc_23[k] = e_0 * (fs_675_16 * h5_p2 - f_15_2 * h5_p4 + fs_4725_16 * r_2 * h3_p2) + e_1 * (-fs_184815_163592 * h7_p2 + fs_10935_7436 * h7_p4 - fs_675_169 * r_2 * h5_p2 + f_30_13 * r_2 * h5_p4 - fs_4725_484 * r_4 * h3_p2) + e_2 * (-fs_1925_97682 * h9_p2 - fs_275_3757 * h9_p4 + fs_184815_11819522 * r_2 * h7_p2 - fs_10935_537251 * r_2 * h7_p4 + fs_3_169 * r_4 * h5_p2 - f_2_13 * r_4 * h5_p4 + fs_525_20449 * r_6 * h3_p2) - fs_42525_64 * e_3 * h3_p2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m4, ph5_m3, ph5_m2, ph7_m4, ph7_m3, ph7_m2, ph9_m4, ph9_m3, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_24[k] = e_0 * (-fs_75_16 * h5_m3 + fs_525_4 * r_2 * h3_m3) + e_1 * (-fs_405_40898 * h7_m3 + fs_75_169 * r_2 * h5_m3 - fs_525_121 * r_4 * h3_m3) + e_2 * (fs_4400_48841 * h9_m3 + fs_810_5909761 * r_2 * h7_m3 - fs_1_507 * r_4 * h5_m3 + fs_700_61347 * r_6 * h3_m3) - fs_4725_16 * e_3 * h3_m3;

        pc_25[k] = e_0 * (f_15_2 * h5_m4 + fs_675_16 * h5_m2 + fs_4725_16 * r_2 * h3_m2) + e_1 * (-fs_10935_7436 * h7_m4 - fs_184815_163592 * h7_m2 - f_30_13 * r_2 * h5_m4 - fs_675_169 * r_2 * h5_m2 - fs_4725_484 * r_4 * h3_m2) + e_2 * (fs_275_3757 * h9_m4 - fs_1925_97682 * h9_m2 + fs_10935_537251 * r_2 * h7_m4 + fs_184815_11819522 * r_2 * h7_m2 + f_2_13 * r_4 * h5_m4 + fs_3_169 * r_4 * h5_m2 + fs_525_20449 * r_6 * h3_m2) - fs_42525_64 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph5_m5, ph7_m6, ph7_m5, ph7_m1, ph9_m6, ph9_m5, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_26[k] = e_0 * (-f_15_4 * h5_m5 + fs_4725_16 * r_2 * h3_m1) + e_1 * (-fs_16875_14872 * h7_m5 + fs_2025_968 * h7_m1 + f_15_13 * r_2 * h5_m5 - fs_4725_484 * r_4 * h3_m1) + e_2 * (fs_308_3757 * h9_m5 + fs_280_48841 * h9_m1 + fs_16875_1074502 * r_2 * h7_m5 - fs_2025_69938 * r_2 * h7_m1 - f_1_13 * r_4 * h5_m5 + fs_525_20449 * r_6 * h3_m1) - fs_42525_64 * e_3 * h3_m1;

        pc_27[k] = fs_2025_1144 * e_1 * h7_m6 + e_2 * (fs_385_7514 * h9_m6 - fs_2025_82654 * r_2 * h7_m6);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_28[k] = e_0 * (fs_875_16 * h5_p1 - fs_75_32 * h5_p5 - fs_175_8 * r_2 * h3_p1) + e_1 * (fs_33075_20449 * h7_p1 - fs_18225_7436 * h7_p5 - fs_875_169 * r_2 * h5_p1 + fs_75_338 * r_2 * h5_p5 + fs_175_242 * r_4 * h3_p1) + e_2 * (fs_105_48841 * h9_p1 - fs_231_7514 * h9_p5 - fs_132300_5909761 * r_2 * h7_p1 + fs_18225_537251 * r_2 * h7_p5 + fs_35_1521 * r_4 * h5_p1 - fs_1_1014 * r_4 * h5_p5 - fs_350_184041 * r_6 * h3_p1) + fs_1575_32 * e_3 * h3_p1;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_29[k] = e_0 * (-fs_175_8 * h5_0 + fs_1125_32 * h5_p4 + fs_350 * r_2 * h3_0) + e_1 * (fs_78750_20449 * h7_0 + fs_675_3718 * h7_p4 + fs_350_169 * r_2 * h5_0 - fs_1125_338 * r_2 * h5_p4 - fs_1400_121 * r_4 * h3_0) + e_2 * (fs_1134_48841 * h9_0 - fs_495_7514 * h9_p4 - fs_315000_5909761 * r_2 * h7_0 - fs_1350_537251 * r_2 * h7_p4 - fs_14_1521 * r_4 * h5_0 + fs_5_338 * r_4 * h5_p4 + fs_5600_184041 * r_6 * h3_0) - fs_1575_2 * e_3 * h3_0;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_p1, ph3_p3, ph5_m2, ph5_p1, ph5_p3, ph7_m2, ph7_p1, ph7_p3, ph9_m2, ph9_p1, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_30[k] = e_0 * (fs_525_16 * h5_p1 - fs_1225_32 * h5_p3 + fs_2625_8 * r_2 * h3_p1 + fs_175_8 * r_2 * h3_p3) + e_1 * (-fs_23805_81796 * h7_p1 + fs_34560_20449 * h7_p3 - fs_525_169 * r_2 * h5_p1 + fs_1225_338 * r_2 * h5_p3 - fs_2625_242 * r_4 * h3_p1 - fs_175_242 * r_4 * h3_p3) + e_2 * (-fs_1575_48841 * h9_p1 - fs_7425_97682 * h9_p3 + fs_23805_5909761 * r_2 * h7_p1 - fs_138240_5909761 * r_2 * h7_p3 + fs_7_507 * r_4 * h5_p1 - fs_49_3042 * r_4 * h5_p3 + fs_1750_61347 * r_6 * h3_p1 + fs_350_184041 * r_6 * h3_p3) + e_3 * (-fs_23625_32 * h3_p1 - fs_1575_32 * h3_p3);

        pc_31[k] = e_0 * (fs_25_2 * h5_m2 + fs_350 * r_2 * h3_m2) + e_1 * (-fs_19845_20449 * h7_m2 - fs_200_169 * r_2 * h5_m2 - fs_1400_121 * r_4 * h3_m2) + e_2 * (fs_5775_48841 * h9_m2 + fs_79380_5909761 * r_2 * h7_m2 + fs_8_1521 * r_4 * h5_m2 + fs_5600_184041 * r_6 * h3_m2) - fs_1575_2 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m1, ph5_m4, ph5_m3, ph5_m1, ph7_m4, ph7_m3, ph7_m1, ph9_m4, ph9_m3, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_32[k] = e_0 * (fs_1225_32 * h5_m3 + fs_525_16 * h5_m1 - fs_175_8 * r_2 * h3_m3 + fs_2625_8 * r_2 * h3_m1) + e_1 * (-fs_34560_20449 * h7_m3 - fs_23805_81796 * h7_m1 - fs_1225_338 * r_2 * h5_m3 - fs_525_169 * r_2 * h5_m1 + fs_175_242 * r_4 * h3_m3 - fs_2625_242 * r_4 * h3_m1) + e_2 * (fs_7425_97682 * h9_m3 - fs_1575_48841 * h9_m1 + fs_138240_5909761 * r_2 * h7_m3 + fs_23805_5909761 * r_2 * h7_m1 + fs_49_3042 * r_4 * h5_m3 + fs_7_507 * r_4 * h5_m1 - fs_350_184041 * r_6 * h3_m3 + fs_1750_61347 * r_6 * h3_m1) + e_3 * (fs_1575_32 * h3_m3 - fs_23625_32 * h3_m1);

        pc_33[k] = -fs_1125_32 * e_0 * h5_m4 + e_1 * (-fs_675_3718 * h7_m4 + fs_1125_338 * r_2 * h5_m4) + e_2 * (fs_495_7514 * h9_m4 + fs_1350_537251 * r_2 * h7_m4 - fs_5_338 * r_4 * h5_m4);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_34[k] = e_0 * (fs_75_32 * h5_m5 - fs_875_16 * h5_m1 + fs_175_8 * r_2 * h3_m1) + e_1 * (fs_18225_7436 * h7_m5 - fs_33075_20449 * h7_m1 - fs_75_338 * r_2 * h5_m5 + fs_875_169 * r_2 * h5_m1 - fs_175_242 * r_4 * h3_m1) + e_2 * (fs_231_7514 * h9_m5 - fs_105_48841 * h9_m1 - fs_18225_537251 * r_2 * h7_m5 + fs_132300_5909761 * r_2 * h7_m1 + fs_1_1014 * r_4 * h5_m5 - fs_35_1521 * r_4 * h5_m1 + fs_350_184041 * r_6 * h3_m1) - fs_1575_32 * e_3 * h3_m1;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_35[k] = e_0 * (-fs_1225_32 * h5_p2 - fs_75_8 * h5_p4 + fs_175_32 * r_2 * h3_p2) + e_1 * (-fs_178605_81796 * h7_p2 - fs_10125_3718 * h7_p4 + fs_1225_338 * r_2 * h5_p2 + fs_150_169 * r_2 * h5_p4 - fs_175_968 * r_4 * h3_p2) + e_2 * (-fs_231_48841 * h9_p2 - fs_66_3757 * h9_p4 + fs_178605_5909761 * r_2 * h7_p2 + fs_20250_537251 * r_2 * h7_p4 - fs_49_3042 * r_4 * h5_p2 - fs_2_507 * r_4 * h5_p4 + fs_175_368082 * r_6 * h3_p2) - fs_1575_128 * e_3 * h3_p2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_36[k] = e_0 * (fs_525_16 * h5_p1 + fs_50 * h5_p3 - fs_2625_32 * r_2 * h3_p1 - fs_175_32 * r_2 * h3_p3) + e_1 * (-fs_108045_81796 * h7_p1 - fs_3375_81796 * h7_p3 - fs_525_169 * r_2 * h5_p1 - fs_800_169 * r_2 * h5_p3 + fs_2625_968 * r_4 * h3_p1 + fs_175_968 * r_4 * h3_p3) + e_2 * (-fs_1008_48841 * h9_p1 - fs_2376_48841 * h9_p3 + fs_108045_5909761 * r_2 * h7_p1 + fs_3375_5909761 * r_2 * h7_p3 + fs_7_507 * r_4 * h5_p1 + fs_32_1521 * r_4 * h5_p3 - fs_875_122694 * r_6 * h3_p1 - fs_175_368082 * r_6 * h3_p3) + e_3 * (fs_23625_128 * h3_p1 + fs_1575_128 * h3_p3);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph3_0, ph3_p2, ph5_m1, ph5_0, ph5_p2, ph7_m1, ph7_0, ph7_p2, ph9_m1, ph9_0, ph9_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_37[k] = e_0 * (fs_175_8 * h5_0 - fs_375_32 * h5_p2 + fs_4375_8 * r_2 * h3_0 + fs_2625_32 * r_2 * h3_p2) + e_1 * (fs_126_20449 * h7_0 + fs_93987_81796 * h7_p2 - fs_350_169 * r_2 * h5_0 + fs_375_338 * r_2 * h5_p2 - fs_4375_242 * r_4 * h3_0 - fs_2625_968 * r_4 * h3_p2) + e_2 * (-fs_4536_48841 * h9_0 - fs_3465_48841 * h9_p2 - fs_504_5909761 * r_2 * h7_0 - fs_93987_5909761 * r_2 * h7_p2 + fs_14_1521 * r_4 * h5_0 - fs_5_1014 * r_4 * h5_p2 + fs_8750_184041 * r_6 * h3_0 + fs_875_122694 * r_6 * h3_p2) + e_3 * (-fs_39375_32 * h3_0 - fs_23625_128 * h3_p2);

        pc_38[k] = e_0 * (fs_875_16 * h5_m1 + fs_4375_8 * r_2 * h3_m1) + e_1 * (-fs_49923_20449 * h7_m1 - fs_875_169 * r_2 * h5_m1 - fs_4375_242 * r_4 * h3_m1) + e_2 * (fs_6720_48841 * h9_m1 + fs_199692_5909761 * r_2 * h7_m1 + fs_35_1521 * r_4 * h5_m1 + fs_8750_184041 * r_6 * h3_m1) - fs_39375_32 * e_3 * h3_m1;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph3_m1, ph5_m3, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ph9_m3, ph9_m2, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_39[k] = e_0 * (fs_375_32 * h5_m2 - fs_2625_32 * r_2 * h3_m2) + e_1 * (-fs_93987_81796 * h7_m2 - fs_375_338 * r_2 * h5_m2 + fs_2625_968 * r_4 * h3_m2) + e_2 * (fs_3465_48841 * h9_m2 + fs_93987_5909761 * r_2 * h7_m2 + fs_5_1014 * r_4 * h5_m2 - fs_875_122694 * r_6 * h3_m2) + fs_23625_128 * e_3 * h3_m2;

        pc_40[k] = e_0 * (-fs_50 * h5_m3 - fs_525_16 * h5_m1 + fs_175_32 * r_2 * h3_m3 + fs_2625_32 * r_2 * h3_m1) + e_1 * (fs_3375_81796 * h7_m3 + fs_108045_81796 * h7_m1 + fs_800_169 * r_2 * h5_m3 + fs_525_169 * r_2 * h5_m1 - fs_175_968 * r_4 * h3_m3 - fs_2625_968 * r_4 * h3_m1) + e_2 * (fs_2376_48841 * h9_m3 + fs_1008_48841 * h9_m1 - fs_3375_5909761 * r_2 * h7_m3 - fs_108045_5909761 * r_2 * h7_m1 - fs_32_1521 * r_4 * h5_m3 - fs_7_507 * r_4 * h5_m1 + fs_175_368082 * r_6 * h3_m3 + fs_875_122694 * r_6 * h3_m1) + e_3 * (-fs_1575_128 * h3_m3 - fs_23625_128 * h3_m1);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m4, ph5_m3, ph5_m2, ph7_m4, ph7_m3, ph7_m2, ph9_m4, ph9_m3, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_41[k] = e_0 * (fs_75_8 * h5_m4 + fs_1225_32 * h5_m2 - fs_175_32 * r_2 * h3_m2) + e_1 * (fs_10125_3718 * h7_m4 + fs_178605_81796 * h7_m2 - fs_150_169 * r_2 * h5_m4 - fs_1225_338 * r_2 * h5_m2 + fs_175_968 * r_4 * h3_m2) + e_2 * (fs_66_3757 * h9_m4 + fs_231_48841 * h9_m2 - fs_20250_537251 * r_2 * h7_m4 - fs_178605_5909761 * r_2 * h7_m2 + fs_2_507 * r_4 * h5_m4 + fs_49_3042 * r_4 * h5_m2 - fs_175_368082 * r_6 * h3_m2) + fs_1575_128 * e_3 * h3_m2;

        pc_42[k] = e_0 * (fs_175_4 * h5_m3 - f_5_4 * r_2 * h3_m3) + e_1 * (fs_212625_40898 * h7_m3 - fs_700_169 * r_2 * h5_m3 + f_5_22 * r_4 * h3_m3) + e_2 * (fs_924_48841 * h9_m3 - fs_425250_5909761 * r_2 * h7_m3 + fs_28_1521 * r_4 * h5_m3 - f_5_429 * r_6 * h3_m3) + f_15_8 * e_3 * h3_m3;

        pc_43[k] = e_0 * (-fs_1575_16 * h5_m2 + f_15_2 * r_2 * h3_m2) + e_1 * (fs_22680_20449 * h7_m2 + fs_1575_169 * r_2 * h5_m2 - f_15_11 * r_4 * h3_m2) + e_2 * (fs_3234_48841 * h9_m2 - fs_90720_5909761 * r_2 * h7_m2 - fs_7_169 * r_4 * h5_m2 + f_10_143 * r_6 * h3_m2) - f_45_4 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph3_0, ph3_p1, ph5_0, ph7_m1, ph7_0, ph7_p1, ph9_m1, ph9_0, ph9_p1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_44[k] = -f_75_4 * e_0 * r_2 * h3_m1 + e_1 * (-fs_189_242 * h7_m1 + f_75_22 * r_4 * h3_m1) + e_2 * (fs_5880_48841 * h9_m1 + fs_378_34969 * r_2 * h7_m1 - f_25_143 * r_6 * h3_m1) + f_225_8 * e_3 * h3_m1;

        pc_45[k] = e_0 * (f_35_4 * h5_0 + f_25 * r_2 * h3_0) + e_1 * (-f_252_143 * h7_0 - f_35_13 * r_2 * h5_0 - f_50_11 * r_4 * h3_0) + e_2 * (f_84_221 * h9_0 + f_504_2431 * r_2 * h7_0 + f_7_39 * r_4 * h5_0 + f_100_429 * r_6 * h3_0) - f_75_2 * e_3 * h3_0;

        pc_46[k] = -f_75_4 * e_0 * r_2 * h3_p1 + e_1 * (-fs_189_242 * h7_p1 + f_75_22 * r_4 * h3_p1) + e_2 * (fs_5880_48841 * h9_p1 + fs_378_34969 * r_2 * h7_p1 - f_25_143 * r_6 * h3_p1) + f_225_8 * e_3 * h3_p1;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph9_p2, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];

        pc_47[k] = e_0 * (-fs_1575_16 * h5_p2 + f_15_2 * r_2 * h3_p2) + e_1 * (fs_22680_20449 * h7_p2 + fs_1575_169 * r_2 * h5_p2 - f_15_11 * r_4 * h3_p2) + e_2 * (fs_3234_48841 * h9_p2 - fs_90720_5909761 * r_2 * h7_p2 - fs_7_169 * r_4 * h5_p2 + f_10_143 * r_6 * h3_p2) - f_45_4 * e_3 * h3_p2;

        pc_48[k] = e_0 * (fs_175_4 * h5_p3 - f_5_4 * r_2 * h3_p3) + e_1 * (fs_212625_40898 * h7_p3 - fs_700_169 * r_2 * h5_p3 + f_5_22 * r_4 * h3_p3) + e_2 * (fs_924_48841 * h9_p3 - fs_425250_5909761 * r_2 * h7_p3 + fs_28_1521 * r_4 * h5_p3 - f_5_429 * r_6 * h3_p3) + f_15_8 * e_3 * h3_p3;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];

        pc_49[k] = e_0 * (fs_75_8 * h5_m4 - fs_1225_32 * h5_m2 + fs_175_32 * r_2 * h3_m2) + e_1 * (fs_10125_3718 * h7_m4 - fs_178605_81796 * h7_m2 - fs_150_169 * r_2 * h5_m4 + fs_1225_338 * r_2 * h5_m2 - fs_175_968 * r_4 * h3_m2) + e_2 * (fs_66_3757 * h9_m4 - fs_231_48841 * h9_m2 - fs_20250_537251 * r_2 * h7_m4 + fs_178605_5909761 * r_2 * h7_m2 + fs_2_507 * r_4 * h5_m4 - fs_49_3042 * r_4 * h5_m2 + fs_175_368082 * r_6 * h3_m2) - fs_1575_128 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph3_m1, ph5_m3, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ph9_m3, ph9_m2, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_50[k] = e_0 * (-fs_50 * h5_m3 + fs_525_16 * h5_m1 + fs_175_32 * r_2 * h3_m3 - fs_2625_32 * r_2 * h3_m1) + e_1 * (fs_3375_81796 * h7_m3 - fs_108045_81796 * h7_m1 + fs_800_169 * r_2 * h5_m3 - fs_525_169 * r_2 * h5_m1 - fs_175_968 * r_4 * h3_m3 + fs_2625_968 * r_4 * h3_m1) + e_2 * (fs_2376_48841 * h9_m3 - fs_1008_48841 * h9_m1 - fs_3375_5909761 * r_2 * h7_m3 + fs_108045_5909761 * r_2 * h7_m1 - fs_32_1521 * r_4 * h5_m3 + fs_7_507 * r_4 * h5_m1 + fs_175_368082 * r_6 * h3_m3 - fs_875_122694 * r_6 * h3_m1) + e_3 * (-fs_1575_128 * h3_m3 + fs_23625_128 * h3_m1);

        pc_51[k] = e_0 * (fs_375_32 * h5_m2 - fs_2625_32 * r_2 * h3_m2) + e_1 * (-fs_93987_81796 * h7_m2 - fs_375_338 * r_2 * h5_m2 + fs_2625_968 * r_4 * h3_m2) + e_2 * (fs_3465_48841 * h9_m2 + fs_93987_5909761 * r_2 * h7_m2 + fs_5_1014 * r_4 * h5_m2 - fs_875_122694 * r_6 * h3_m2) + fs_23625_128 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph3_p2, ph5_0, ph5_p1, ph5_p2, ph7_0, ph7_p1, ph7_p2, ph9_0, ph9_p1, ph9_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_52[k] = e_0 * (fs_875_16 * h5_p1 + fs_4375_8 * r_2 * h3_p1) + e_1 * (-fs_49923_20449 * h7_p1 - fs_875_169 * r_2 * h5_p1 - fs_4375_242 * r_4 * h3_p1) + e_2 * (fs_6720_48841 * h9_p1 + fs_199692_5909761 * r_2 * h7_p1 + fs_35_1521 * r_4 * h5_p1 + fs_8750_184041 * r_6 * h3_p1) - fs_39375_32 * e_3 * h3_p1;

        pc_53[k] = e_0 * (fs_175_8 * h5_0 + fs_375_32 * h5_p2 + fs_4375_8 * r_2 * h3_0 - fs_2625_32 * r_2 * h3_p2) + e_1 * (fs_126_20449 * h7_0 - fs_93987_81796 * h7_p2 - fs_350_169 * r_2 * h5_0 - fs_375_338 * r_2 * h5_p2 - fs_4375_242 * r_4 * h3_0 + fs_2625_968 * r_4 * h3_p2) + e_2 * (-fs_4536_48841 * h9_0 + fs_3465_48841 * h9_p2 - fs_504_5909761 * r_2 * h7_0 + fs_93987_5909761 * r_2 * h7_p2 + fs_14_1521 * r_4 * h5_0 + fs_5_1014 * r_4 * h5_p2 + fs_8750_184041 * r_6 * h3_0 - fs_875_122694 * r_6 * h3_p2) + e_3 * (-fs_39375_32 * h3_0 + fs_23625_128 * h3_p2);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_54[k] = e_0 * (fs_525_16 * h5_p1 - fs_50 * h5_p3 - fs_2625_32 * r_2 * h3_p1 + fs_175_32 * r_2 * h3_p3) + e_1 * (-fs_108045_81796 * h7_p1 + fs_3375_81796 * h7_p3 - fs_525_169 * r_2 * h5_p1 + fs_800_169 * r_2 * h5_p3 + fs_2625_968 * r_4 * h3_p1 - fs_175_968 * r_4 * h3_p3) + e_2 * (-fs_1008_48841 * h9_p1 + fs_2376_48841 * h9_p3 + fs_108045_5909761 * r_2 * h7_p1 - fs_3375_5909761 * r_2 * h7_p3 + fs_7_507 * r_4 * h5_p1 - fs_32_1521 * r_4 * h5_p3 - fs_875_122694 * r_6 * h3_p1 + fs_175_368082 * r_6 * h3_p3) + e_3 * (fs_23625_128 * h3_p1 - fs_1575_128 * h3_p3);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_55[k] = e_0 * (-fs_1225_32 * h5_p2 + fs_75_8 * h5_p4 + fs_175_32 * r_2 * h3_p2) + e_1 * (-fs_178605_81796 * h7_p2 + fs_10125_3718 * h7_p4 + fs_1225_338 * r_2 * h5_p2 - fs_150_169 * r_2 * h5_p4 - fs_175_968 * r_4 * h3_p2) + e_2 * (-fs_231_48841 * h9_p2 + fs_66_3757 * h9_p4 + fs_178605_5909761 * r_2 * h7_p2 - fs_20250_537251 * r_2 * h7_p4 - fs_49_3042 * r_4 * h5_p2 + fs_2_507 * r_4 * h5_p4 + fs_175_368082 * r_6 * h3_p2) - fs_1575_128 * e_3 * h3_p2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph5_m5, ph5_m4, ph5_m1, ph7_m5, ph7_m4, ph7_m1, ph9_m5, ph9_m4, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_56[k] = e_0 * (fs_75_32 * h5_m5 + fs_875_16 * h5_m1 - fs_175_8 * r_2 * h3_m1) + e_1 * (fs_18225_7436 * h7_m5 + fs_33075_20449 * h7_m1 - fs_75_338 * r_2 * h5_m5 - fs_875_169 * r_2 * h5_m1 + fs_175_242 * r_4 * h3_m1) + e_2 * (fs_231_7514 * h9_m5 + fs_105_48841 * h9_m1 - fs_18225_537251 * r_2 * h7_m5 - fs_132300_5909761 * r_2 * h7_m1 + fs_1_1014 * r_4 * h5_m5 + fs_35_1521 * r_4 * h5_m1 - fs_350_184041 * r_6 * h3_m1) + fs_1575_32 * e_3 * h3_m1;

        pc_57[k] = -fs_1125_32 * e_0 * h5_m4 + e_1 * (-fs_675_3718 * h7_m4 + fs_1125_338 * r_2 * h5_m4) + e_2 * (fs_495_7514 * h9_m4 + fs_1350_537251 * r_2 * h7_m4 - fs_5_338 * r_4 * h5_m4);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m1, ph3_p2, ph5_m3, ph5_m1, ph5_p2, ph7_m3, ph7_m1, ph7_p2, ph9_m3, ph9_m1, ph9_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_58[k] = e_0 * (fs_1225_32 * h5_m3 - fs_525_16 * h5_m1 - fs_175_8 * r_2 * h3_m3 - fs_2625_8 * r_2 * h3_m1) + e_1 * (-fs_34560_20449 * h7_m3 + fs_23805_81796 * h7_m1 - fs_1225_338 * r_2 * h5_m3 + fs_525_169 * r_2 * h5_m1 + fs_175_242 * r_4 * h3_m3 + fs_2625_242 * r_4 * h3_m1) + e_2 * (fs_7425_97682 * h9_m3 + fs_1575_48841 * h9_m1 + fs_138240_5909761 * r_2 * h7_m3 - fs_23805_5909761 * r_2 * h7_m1 + fs_49_3042 * r_4 * h5_m3 - fs_7_507 * r_4 * h5_m1 - fs_350_184041 * r_6 * h3_m3 - fs_1750_61347 * r_6 * h3_m1) + e_3 * (fs_1575_32 * h3_m3 + fs_23625_32 * h3_m1);

        pc_59[k] = e_0 * (fs_25_2 * h5_p2 + fs_350 * r_2 * h3_p2) + e_1 * (-fs_19845_20449 * h7_p2 - fs_200_169 * r_2 * h5_p2 - fs_1400_121 * r_4 * h3_p2) + e_2 * (fs_5775_48841 * h9_p2 + fs_79380_5909761 * r_2 * h7_p2 + fs_8_1521 * r_4 * h5_p2 + fs_5600_184041 * r_6 * h3_p2) - fs_1575_2 * e_3 * h3_p2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_60[k] = e_0 * (fs_525_16 * h5_p1 + fs_1225_32 * h5_p3 + fs_2625_8 * r_2 * h3_p1 - fs_175_8 * r_2 * h3_p3) + e_1 * (-fs_23805_81796 * h7_p1 - fs_34560_20449 * h7_p3 - fs_525_169 * r_2 * h5_p1 - fs_1225_338 * r_2 * h5_p3 - fs_2625_242 * r_4 * h3_p1 + fs_175_242 * r_4 * h3_p3) + e_2 * (-fs_1575_48841 * h9_p1 + fs_7425_97682 * h9_p3 + fs_23805_5909761 * r_2 * h7_p1 + fs_138240_5909761 * r_2 * h7_p3 + fs_7_507 * r_4 * h5_p1 + fs_49_3042 * r_4 * h5_p3 + fs_1750_61347 * r_6 * h3_p1 - fs_350_184041 * r_6 * h3_p3) + e_3 * (-fs_23625_32 * h3_p1 + fs_1575_32 * h3_p3);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_61[k] = e_0 * (-fs_175_8 * h5_0 - fs_1125_32 * h5_p4 + fs_350 * r_2 * h3_0) + e_1 * (fs_78750_20449 * h7_0 - fs_675_3718 * h7_p4 + fs_350_169 * r_2 * h5_0 + fs_1125_338 * r_2 * h5_p4 - fs_1400_121 * r_4 * h3_0) + e_2 * (fs_1134_48841 * h9_0 + fs_495_7514 * h9_p4 - fs_315000_5909761 * r_2 * h7_0 + fs_1350_537251 * r_2 * h7_p4 - fs_14_1521 * r_4 * h5_0 - fs_5_338 * r_4 * h5_p4 + fs_5600_184041 * r_6 * h3_0) - fs_1575_2 * e_3 * h3_0;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph5_p1, ph5_p5, ph7_m6, ph7_p1, ph7_p5, ph9_m6, ph9_p1, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_62[k] = e_0 * (fs_875_16 * h5_p1 + fs_75_32 * h5_p5 - fs_175_8 * r_2 * h3_p1) + e_1 * (fs_33075_20449 * h7_p1 + fs_18225_7436 * h7_p5 - fs_875_169 * r_2 * h5_p1 - fs_75_338 * r_2 * h5_p5 + fs_175_242 * r_4 * h3_p1) + e_2 * (fs_105_48841 * h9_p1 + fs_231_7514 * h9_p5 - fs_132300_5909761 * r_2 * h7_p1 - fs_18225_537251 * r_2 * h7_p5 + fs_35_1521 * r_4 * h5_p1 + fs_1_1014 * r_4 * h5_p5 - fs_350_184041 * r_6 * h3_p1) + fs_1575_32 * e_3 * h3_p1;

        pc_63[k] = fs_2025_1144 * e_1 * h7_m6 + e_2 * (fs_385_7514 * h9_m6 - fs_2025_82654 * r_2 * h7_m6);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_m1, ph5_m5, ph5_m4, ph5_m2, ph7_m5, ph7_m4, ph7_m2, ph7_m1, ph9_m5, ph9_m4, ph9_m2, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_64[k] = e_0 * (-f_15_4 * h5_m5 - fs_4725_16 * r_2 * h3_m1) + e_1 * (-fs_16875_14872 * h7_m5 - fs_2025_968 * h7_m1 + f_15_13 * r_2 * h5_m5 + fs_4725_484 * r_4 * h3_m1) + e_2 * (fs_308_3757 * h9_m5 - fs_280_48841 * h9_m1 + fs_16875_1074502 * r_2 * h7_m5 + fs_2025_69938 * r_2 * h7_m1 - f_1_13 * r_4 * h5_m5 - fs_525_20449 * r_6 * h3_m1) + fs_42525_64 * e_3 * h3_m1;

        pc_65[k] = e_0 * (f_15_2 * h5_m4 - fs_675_16 * h5_m2 - fs_4725_16 * r_2 * h3_m2) + e_1 * (-fs_10935_7436 * h7_m4 + fs_184815_163592 * h7_m2 - f_30_13 * r_2 * h5_m4 + fs_675_169 * r_2 * h5_m2 + fs_4725_484 * r_4 * h3_m2) + e_2 * (fs_275_3757 * h9_m4 + fs_1925_97682 * h9_m2 + fs_10935_537251 * r_2 * h7_m4 - fs_184815_11819522 * r_2 * h7_m2 + f_2_13 * r_4 * h5_m4 - fs_3_169 * r_4 * h5_m2 - fs_525_20449 * r_6 * h3_m2) + fs_42525_64 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p4, ph7_p2, ph7_p3, ph7_p4, ph9_p2, ph9_p3, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_66[k] = e_0 * (-fs_75_16 * h5_p3 + fs_525_4 * r_2 * h3_p3) + e_1 * (-fs_405_40898 * h7_p3 + fs_75_169 * r_2 * h5_p3 - fs_525_121 * r_4 * h3_p3) + e_2 * (fs_4400_48841 * h9_p3 + fs_810_5909761 * r_2 * h7_p3 - fs_1_507 * r_4 * h5_p3 + fs_700_61347 * r_6 * h3_p3) - fs_4725_16 * e_3 * h3_p3;

        pc_67[k] = e_0 * (fs_675_16 * h5_p2 + f_15_2 * h5_p4 + fs_4725_16 * r_2 * h3_p2) + e_1 * (-fs_184815_163592 * h7_p2 - fs_10935_7436 * h7_p4 - fs_675_169 * r_2 * h5_p2 - f_30_13 * r_2 * h5_p4 - fs_4725_484 * r_4 * h3_p2) + e_2 * (-fs_1925_97682 * h9_p2 + fs_275_3757 * h9_p4 + fs_184815_11819522 * r_2 * h7_p2 + fs_10935_537251 * r_2 * h7_p4 + fs_3_169 * r_4 * h5_p2 + f_2_13 * r_4 * h5_p4 + fs_525_20449 * r_6 * h3_p2) - fs_42525_64 * e_3 * h3_p2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph5_0, ph5_p5, ph7_0, ph7_p1, ph7_p5, ph7_p6, ph9_0, ph9_p1, ph9_p5, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_68[k] = e_0 * (-f_15_4 * h5_p5 + fs_4725_16 * r_2 * h3_p1) + e_1 * (fs_2025_968 * h7_p1 - fs_16875_14872 * h7_p5 + f_15_13 * r_2 * h5_p5 - fs_4725_484 * r_4 * h3_p1) + e_2 * (fs_280_48841 * h9_p1 + fs_308_3757 * h9_p5 - fs_2025_69938 * r_2 * h7_p1 + fs_16875_1074502 * r_2 * h7_p5 - f_1_13 * r_4 * h5_p5 + fs_525_20449 * r_6 * h3_p1) - fs_42525_64 * e_3 * h3_p1;

        pc_69[k] = e_0 * (-fs_525_4 * h5_0 + fs_525_4 * r_2 * h3_0) + e_1 * (-fs_42525_20449 * h7_0 + fs_2025_1144 * h7_p6 + fs_2100_169 * r_2 * h5_0 - fs_525_121 * r_4 * h3_0) + e_2 * (-fs_84_48841 * h9_0 + fs_385_7514 * h9_p6 + fs_170100_5909761 * r_2 * h7_0 - fs_2025_82654 * r_2 * h7_p6 - fs_28_507 * r_4 * h5_0 + fs_700_61347 * r_6 * h3_0) - fs_4725_16 * e_3 * h3_0;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m7, ph7_m6, ph7_m2, ph7_m1, ph9_m7, ph9_m6, ph9_m2, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_70[k] = e_0 * (fs_525_8 * h5_m1 - fs_2625_16 * r_2 * h3_m1) + e_1 * (fs_945_1144 * h7_m7 + fs_91125_163592 * h7_m1 - fs_1050_169 * r_2 * h5_m1 + fs_2625_484 * r_4 * h3_m1) + e_2 * (fs_308_3757 * h9_m7 + fs_14_48841 * h9_m1 - fs_945_82654 * r_2 * h7_m7 - fs_91125_11819522 * r_2 * h7_m1 + fs_14_507 * r_4 * h5_m1 - fs_875_61347 * r_6 * h3_m1) + fs_23625_64 * e_3 * h3_m1;

        pc_71[k] = e_0 * (-f_15_4 * h5_m2 - fs_1575_4 * r_2 * h3_m2) + e_1 * (-fs_360_143 * h7_m6 - fs_36000_20449 * h7_m2 + f_15_13 * r_2 * h5_m2 + fs_1575_121 * r_4 * h3_m2) + e_2 * (fs_693_7514 * h9_m6 - fs_231_97682 * h9_m2 + fs_1440_41327 * r_2 * h7_m6 + fs_144000_5909761 * r_2 * h7_m2 - f_1_13 * r_4 * h5_m2 - fs_700_20449 * r_6 * h3_m2) + fs_14175_16 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m5, ph5_m3, ph5_p4, ph7_m5, ph7_m3, ph7_p4, ph9_m5, ph9_m3, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_72[k] = e_0 * (fs_675_16 * h5_m5 - fs_375_16 * h5_m3 - fs_2625_16 * r_2 * h3_m3) + e_1 * (-fs_7569_14872 * h7_m5 + fs_328329_163592 * h7_m3 - fs_675_169 * r_2 * h5_m5 + fs_375_169 * r_2 * h5_m3 + fs_2625_484 * r_4 * h3_m3) + e_2 * (fs_231_3757 * h9_m5 + fs_495_48841 * h9_m3 + fs_7569_1074502 * r_2 * h7_m5 - fs_328329_11819522 * r_2 * h7_m3 + fs_3_169 * r_4 * h5_m5 - fs_5_507 * r_4 * h5_m3 - fs_875_61347 * r_6 * h3_m3) + fs_23625_64 * e_3 * h3_m3;

        pc_73[k] = -fs_1125_16 * e_0 * h5_p4 + e_1 * (fs_1728_1859 * h7_p4 + fs_1125_169 * r_2 * h5_p4) + e_2 * (fs_220_3757 * h9_p4 - fs_6912_537251 * r_2 * h7_p4 - fs_5_169 * r_4 * h5_p4);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p5, ph7_p2, ph7_p3, ph7_p5, ph7_p6, ph9_p2, ph9_p3, ph9_p5, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_74[k] = e_0 * (fs_375_16 * h5_p3 + fs_675_16 * h5_p5 + fs_2625_16 * r_2 * h3_p3) + e_1 * (-fs_328329_163592 * h7_p3 - fs_7569_14872 * h7_p5 - fs_375_169 * r_2 * h5_p3 - fs_675_169 * r_2 * h5_p5 - fs_2625_484 * r_4 * h3_p3) + e_2 * (-fs_495_48841 * h9_p3 + fs_231_3757 * h9_p5 + fs_328329_11819522 * r_2 * h7_p3 + fs_7569_1074502 * r_2 * h7_p5 + fs_5_507 * r_4 * h5_p3 + fs_3_169 * r_4 * h5_p5 + fs_875_61347 * r_6 * h3_p3) - fs_23625_64 * e_3 * h3_p3;

        pc_75[k] = e_0 * (f_15_4 * h5_p2 + fs_1575_4 * r_2 * h3_p2) + e_1 * (fs_36000_20449 * h7_p2 - fs_360_143 * h7_p6 - f_15_13 * r_2 * h5_p2 - fs_1575_121 * r_4 * h3_p2) + e_2 * (fs_231_97682 * h9_p2 + fs_693_7514 * h9_p6 - fs_144000_5909761 * r_2 * h7_p2 + fs_1440_41327 * r_2 * h7_p6 + f_1_13 * r_4 * h5_p2 + fs_700_20449 * r_6 * h3_p2) - fs_14175_16 * e_3 * h3_p2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_p1, ph5_m2, ph5_p1, ph7_m2, ph7_p1, ph7_p7, ph9_m8, ph9_m2, ph9_p1, ph9_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_76[k] = e_0 * (-fs_525_8 * h5_p1 + fs_2625_16 * r_2 * h3_p1) + e_1 * (-fs_91125_163592 * h7_p1 + fs_945_1144 * h7_p7 + fs_1050_169 * r_2 * h5_p1 - fs_2625_484 * r_4 * h3_p1) + e_2 * (-fs_14_48841 * h9_p1 + fs_308_3757 * h9_p7 + fs_91125_11819522 * r_2 * h7_p1 - fs_945_82654 * r_2 * h7_p7 - fs_14_507 * r_4 * h5_p1 + fs_875_61347 * r_6 * h3_p1) - fs_23625_64 * e_3 * h3_p1;

        pc_77[k] = e_0 * (fs_825_16 * h5_m2 - fs_5775_16 * r_2 * h3_m2) + e_1 * (fs_3375_14872 * h7_m2 - fs_825_169 * r_2 * h5_m2 + fs_525_44 * r_4 * h3_m2) + e_2 * (fs_28_221 * h9_m8 + fs_7_97682 * h9_m2 - fs_3375_1074502 * r_2 * h7_m2 + fs_11_507 * r_4 * h5_m2 - fs_175_5577 * r_6 * h3_m2) + fs_51975_64 * e_3 * h3_m2;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m3, ph5_p5, ph7_m7, ph7_m6, ph7_m4, ph7_m3, ph7_p5, ph9_m7, ph9_m6, ph9_m4, ph9_m3, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p5 = ph9_p5[k];

        pc_78[k] = e_0 * (-fs_825_16 * h5_m3 - fs_5775_16 * r_2 * h3_m3) + e_1 * (-fs_315_104 * h7_m7 - fs_16245_14872 * h7_m3 + fs_825_169 * r_2 * h5_m3 + fs_525_44 * r_4 * h3_m3) + e_2 * (fs_336_3757 * h9_m7 - f_6_221 * h9_m3 + fs_315_7514 * r_2 * h7_m7 + fs_16245_1074502 * r_2 * h7_m3 - fs_11_507 * r_4 * h5_m3 - fs_175_5577 * r_6 * h3_m3) + fs_51975_64 * e_3 * h3_m3;

        pc_79[k] = e_1 * (fs_9_104 * h7_m6 + f_3_2 * h7_m4) + e_2 * (fs_315_7514 * h9_m6 + fs_15_3757 * h9_m4 - fs_9_7514 * r_2 * h7_m6 - f_3_17 * r_2 * h7_m4);

        pc_80[k] = -fs_2475_16 * e_0 * h5_p5 + e_1 * (fs_1323_338 * h7_p5 + fs_2475_169 * r_2 * h5_p5) + e_2 * (fs_112_3757 * h9_p5 - fs_2646_48841 * r_2 * h7_p5 - fs_11_169 * r_4 * h5_p5);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph7_p3, ph7_p4, ph7_p6, ph7_p7, ph9_p3, ph9_p4, ph9_p6, ph9_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_81[k] = e_1 * (-f_3_2 * h7_p4 + fs_9_104 * h7_p6) + e_2 * (-fs_15_3757 * h9_p4 + fs_315_7514 * h9_p6 + f_3_17 * r_2 * h7_p4 - fs_9_7514 * r_2 * h7_p6);

        pc_82[k] = e_0 * (fs_825_16 * h5_p3 + fs_5775_16 * r_2 * h3_p3) + e_1 * (fs_16245_14872 * h7_p3 - fs_315_104 * h7_p7 - fs_825_169 * r_2 * h5_p3 - fs_525_44 * r_4 * h3_p3) + e_2 * (f_6_221 * h9_p3 + fs_336_3757 * h9_p7 - fs_16245_1074502 * r_2 * h7_p3 + fs_315_7514 * r_2 * h7_p7 + fs_11_507 * r_4 * h5_p3 + fs_175_5577 * r_6 * h3_p3) - fs_51975_64 * e_3 * h3_p3;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_p2, ph5_m3, ph5_p2, ph7_m3, ph7_p2, ph9_m9, ph9_m3, ph9_p2, ph9_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];

        pc_83[k] = e_0 * (-fs_825_16 * h5_p2 + fs_5775_16 * r_2 * h3_p2) + e_1 * (-fs_3375_14872 * h7_p2 + fs_825_169 * r_2 * h5_p2 - fs_525_44 * r_4 * h3_p2) + e_2 * (-fs_7_97682 * h9_p2 + fs_28_221 * h9_p8 + fs_3375_1074502 * r_2 * h7_p2 - fs_11_507 * r_4 * h5_p2 + fs_175_5577 * r_6 * h3_p2) - fs_51975_64 * e_3 * h3_p2;

        pc_84[k] = e_0 * (fs_825_32 * h5_m3 - fs_5775_8 * r_2 * h3_m3) + e_1 * (fs_405_7436 * h7_m3 - fs_825_338 * r_2 * h5_m3 + fs_525_22 * r_4 * h3_m3) + e_2 * (fs_42_221 * h9_m9 + fs_1_97682 * h9_m3 - fs_405_537251 * r_2 * h7_m3 + fs_11_1014 * r_4 * h5_m3 - fs_350_5577 * r_6 * h3_m3) + fs_51975_32 * e_3 * h3_m3;
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_m5, ph5_m4, ph7_m7, ph7_m5, ph7_m4, ph7_p6, ph9_m8, ph9_m7, ph9_m5, ph9_m4, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p6 = ph9_p6[k];

        pc_85[k] = -fs_2475_32 * e_0 * h5_m4 + e_1 * (-fs_135_338 * h7_m4 + fs_2475_338 * r_2 * h5_m4) + e_2 * (fs_14_221 * h9_m8 - fs_1_7514 * h9_m4 + fs_270_48841 * r_2 * h7_m4 - fs_11_338 * r_4 * h5_m4);

        pc_86[k] = fs_2475_32 * e_0 * h5_m5 + e_1 * (fs_189_52 * h7_m7 + fs_243_169 * h7_m5 - fs_2475_338 * r_2 * h5_m5) + e_2 * (fs_70_3757 * h9_m7 + fs_7_7514 * h9_m5 - fs_189_3757 * r_2 * h7_m7 - fs_972_48841 * r_2 * h7_m5 + fs_11_338 * r_4 * h5_m5);

        pc_87[k] = fs_81_13 * e_1 * h7_p6 + e_2 * (fs_35_3757 * h9_p6 - fs_324_3757 * r_2 * h7_p6);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_p4, ph5_p5, ph7_p4, ph7_p5, ph7_p7, ph9_p4, ph9_p5, ph9_p7, ph9_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h9_p8 = ph9_p8[k];

        pc_88[k] = -fs_2475_32 * e_0 * h5_p5 + e_1 * (-fs_243_169 * h7_p5 + fs_189_52 * h7_p7 + fs_2475_338 * r_2 * h5_p5) + e_2 * (-fs_7_7514 * h9_p5 + fs_70_3757 * h9_p7 + fs_972_48841 * r_2 * h7_p5 - fs_189_3757 * r_2 * h7_p7 - fs_11_338 * r_4 * h5_p5);

        pc_89[k] = fs_2475_32 * e_0 * h5_p4 + e_1 * (fs_135_338 * h7_p4 - fs_2475_338 * r_2 * h5_p4) + e_2 * (fs_1_7514 * h9_p4 + fs_14_221 * h9_p8 - fs_270_48841 * r_2 * h7_p4 + fs_11_338 * r_4 * h5_p4);
    }

    // NOTE: the rows are formed in 48 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph9_p9, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p9 = ph9_p9[k];

        pc_90[k] = e_0 * (-fs_825_32 * h5_p3 + fs_5775_8 * r_2 * h3_p3) + e_1 * (-fs_405_7436 * h7_p3 + fs_825_338 * r_2 * h5_p3 - fs_525_22 * r_4 * h3_p3) + e_2 * (-fs_1_97682 * h9_p3 + fs_42_221 * h9_p9 + fs_405_537251 * r_2 * h7_p3 - fs_11_1014 * r_4 * h5_p3 + fs_350_5577 * r_6 * h3_p3) - fs_51975_32 * e_3 * h3_p3;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[91] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90};

    for (size_t m = 0; m < 91; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
