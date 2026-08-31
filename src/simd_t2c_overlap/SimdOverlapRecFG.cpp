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



#include "SimdOverlapRecFG.hpp"

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
compute_fg_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 3) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFG.compute_fg_overlap: Basis functions must be of angular momenta three and four"));
    }

    if (harmonics.size() < 7)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFG.compute_fg_overlap: Harmonics must reach angular momentum seven"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFG.compute_fg_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 63 * nvalues, 0.0);

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

            const auto f_0 = fbase * aexp * fmu / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * aexp * fmu * fmu / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * aexp * fmu * fmu * fmu / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * aexp / fexp / fexp / fexp / fexp;

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_189_16 = std::sqrt(11.8125);
    const auto fs_567_8 = std::sqrt(70.875);
    const auto fs_5_42 = std::sqrt(5.0 / 42.0);
    const auto fs_7_3 = std::sqrt(7.0 / 3.0);
    const auto fs_81_14 = std::sqrt(81.0 / 14.0);
    const auto fs_25_368082 = std::sqrt(25.0 / 368082.0);
    const auto fs_175_858 = std::sqrt(175.0 / 858.0);
    const auto fs_10_3549 = std::sqrt(10.0 / 3549.0);
    const auto fs_7_363 = std::sqrt(7.0 / 363.0);
    const auto fs_2_63 = std::sqrt(2.0 / 63.0);
    const auto fs_1575_32 = std::sqrt(49.21875);
    const auto fs_567_16 = std::sqrt(35.4375);
    const auto fs_25_28 = std::sqrt(25.0 / 28.0);
    const auto fs_7 = std::sqrt(7.0);
    const auto fs_81_28 = std::sqrt(81.0 / 28.0);
    const auto fs_175_184041 = std::sqrt(175.0 / 184041.0);
    const auto fs_50_429 = std::sqrt(50.0 / 429.0);
    const auto fs_25_1183 = std::sqrt(25.0 / 1183.0);
    const auto fs_7_121 = std::sqrt(7.0 / 121.0);
    const auto fs_1_63 = std::sqrt(1.0 / 63.0);
    const auto fs_1575_64 = std::sqrt(24.609375);
    const auto fs_243_16 = std::sqrt(15.1875);
    const auto fs_81_32 = std::sqrt(2.53125);
    const auto fs_375_392 = std::sqrt(375.0 / 392.0);
    const auto fs_3 = std::sqrt(3.0);
    const auto fs_81_392 = std::sqrt(81.0 / 392.0);
    const auto fs_350_184041 = std::sqrt(350.0 / 184041.0);
    const auto fs_350_5577 = std::sqrt(350.0 / 5577.0);
    const auto fs_375_16562 = std::sqrt(375.0 / 16562.0);
    const auto fs_3_121 = std::sqrt(3.0 / 121.0);
    const auto fs_1_882 = std::sqrt(1.0 / 882.0);
    const auto fs_225_128 = std::sqrt(1.7578125);
    const auto fs_135_16 = std::sqrt(8.4375);
    const auto fs_125_84 = std::sqrt(125.0 / 84.0);
    const auto fs_45_28 = std::sqrt(45.0 / 28.0);
    const auto fs_5_3 = std::sqrt(5.0 / 3.0);
    const auto fs_350_61347 = std::sqrt(350.0 / 61347.0);
    const auto fs_175_5577 = std::sqrt(175.0 / 5577.0);
    const auto fs_125_3549 = std::sqrt(125.0 / 3549.0);
    const auto fs_45_1183 = std::sqrt(45.0 / 1183.0);
    const auto fs_5_363 = std::sqrt(5.0 / 363.0);
    const auto f_9_4 = 2.25;
    const auto fs_25_7 = std::sqrt(25.0 / 7.0);
    const auto f_1 = 1.0;
    const auto fs_1750_61347 = std::sqrt(1750.0 / 61347.0);
    const auto fs_100_1183 = std::sqrt(100.0 / 1183.0);
    const auto f_1_11 = 1.0 / 11.0;
    const auto fs_315_16 = std::sqrt(19.6875);
    const auto fs_5_9 = std::sqrt(5.0 / 9.0);
    const auto fs_35_9 = std::sqrt(35.0 / 9.0);
    const auto fs_25_40898 = std::sqrt(25.0 / 40898.0);
    const auto fs_25_286 = std::sqrt(25.0 / 286.0);
    const auto fs_20_1521 = std::sqrt(20.0 / 1521.0);
    const auto fs_35_1089 = std::sqrt(35.0 / 1089.0);
    const auto fs_63_16 = std::sqrt(3.9375);
    const auto fs_1701_32 = std::sqrt(53.15625);
    const auto fs_605_504 = std::sqrt(605.0 / 504.0);
    const auto fs_25_12 = std::sqrt(25.0 / 12.0);
    const auto fs_7_9 = std::sqrt(7.0 / 9.0);
    const auto fs_243_56 = std::sqrt(243.0 / 56.0);
    const auto fs_200_61347 = std::sqrt(200.0 / 61347.0);
    const auto fs_200_1859 = std::sqrt(200.0 / 1859.0);
    const auto fs_605_21294 = std::sqrt(605.0 / 21294.0);
    const auto fs_25_507 = std::sqrt(25.0 / 507.0);
    const auto fs_7_1089 = std::sqrt(7.0 / 1089.0);
    const auto fs_1_42 = std::sqrt(1.0 / 42.0);
    const auto fs_4725_128 = std::sqrt(36.9140625);
    const auto fs_27_16 = std::sqrt(1.6875);
    const auto fs_243_4 = std::sqrt(60.75);
    const auto fs_400_147 = std::sqrt(400.0 / 147.0);
    const auto fs_20_21 = std::sqrt(20.0 / 21.0);
    const auto fs_1_3 = std::sqrt(1.0 / 3.0);
    const auto fs_243_49 = std::sqrt(243.0 / 49.0);
    const auto fs_1225_61347 = std::sqrt(1225.0 / 61347.0);
    const auto fs_175_1859 = std::sqrt(175.0 / 1859.0);
    const auto fs_1600_24843 = std::sqrt(1600.0 / 24843.0);
    const auto fs_80_3549 = std::sqrt(80.0 / 3549.0);
    const auto fs_1_363 = std::sqrt(1.0 / 363.0);
    const auto fs_4_147 = std::sqrt(4.0 / 147.0);
    const auto fs_675_16 = std::sqrt(42.1875);
    const auto f_3 = 3.0;
    const auto fs_243_32 = std::sqrt(7.59375);
    const auto fs_3125_3528 = std::sqrt(3125.0 / 3528.0);
    const auto fs_5_84 = std::sqrt(5.0 / 84.0);
    const auto f_4_3 = 4.0 / 3.0;
    const auto fs_243_392 = std::sqrt(243.0 / 392.0);
    const auto fs_1400_61347 = std::sqrt(1400.0 / 61347.0);
    const auto fs_1400_20449 = std::sqrt(1400.0 / 20449.0);
    const auto fs_3125_149058 = std::sqrt(3125.0 / 149058.0);
    const auto fs_5_3549 = std::sqrt(5.0 / 3549.0);
    const auto f_4_33 = 4.0 / 33.0;
    const auto fs_1_294 = std::sqrt(1.0 / 294.0);
    const auto fs_675_128 = std::sqrt(5.2734375);
    const auto f_21_4 = 5.25;
    const auto fs_25_63 = std::sqrt(25.0 / 63.0);
    const auto f_7_3 = 7.0 / 3.0;
    const auto fs_1750_20449 = std::sqrt(1750.0 / 20449.0);
    const auto fs_100_10647 = std::sqrt(100.0 / 10647.0);
    const auto f_7_33 = 7.0 / 33.0;
    const auto fs_4_3 = std::sqrt(4.0 / 3.0);
    const auto fs_125_40898 = std::sqrt(125.0 / 40898.0);
    const auto fs_125_3718 = std::sqrt(125.0 / 3718.0);
    const auto fs_16_507 = std::sqrt(16.0 / 507.0);
    const auto fs_20_507 = std::sqrt(20.0 / 507.0);
    const auto f_7_6 = 7.0 / 6.0;
    const auto fs_1_12 = std::sqrt(1.0 / 12.0);
    const auto fs_250_20449 = std::sqrt(250.0 / 20449.0);
    const auto fs_125_1859 = std::sqrt(125.0 / 1859.0);
    const auto f_7_39 = 7.0 / 39.0;
    const auto fs_1_507 = std::sqrt(1.0 / 507.0);
    const auto fs_45_4 = std::sqrt(11.25);
    const auto fs_1215_32 = std::sqrt(37.96875);
    const auto fs_1849_3528 = std::sqrt(1849.0 / 3528.0);
    const auto fs_27_28 = std::sqrt(27.0 / 28.0);
    const auto fs_20_9 = std::sqrt(20.0 / 9.0);
    const auto fs_1215_392 = std::sqrt(1215.0 / 392.0);
    const auto fs_1849_149058 = std::sqrt(1849.0 / 149058.0);
    const auto fs_27_1183 = std::sqrt(27.0 / 1183.0);
    const auto fs_20_1089 = std::sqrt(20.0 / 1089.0);
    const auto fs_5_294 = std::sqrt(5.0 / 294.0);
    const auto fs_3375_128 = std::sqrt(26.3671875);
    const auto fs_1215_16 = std::sqrt(75.9375);
    const auto fs_5_588 = std::sqrt(5.0 / 588.0);
    const auto fs_64_63 = std::sqrt(64.0 / 63.0);
    const auto fs_1215_196 = std::sqrt(1215.0 / 196.0);
    const auto fs_6125_61347 = std::sqrt(6125.0 / 61347.0);
    const auto fs_5_24843 = std::sqrt(5.0 / 24843.0);
    const auto fs_256_10647 = std::sqrt(256.0 / 10647.0);
    const auto fs_5_147 = std::sqrt(5.0 / 147.0);
    const auto fs_3375_64 = std::sqrt(52.734375);
    const auto f_3_4 = 0.75;
    const auto fs_243_8 = std::sqrt(30.375);
    const auto fs_605_882 = std::sqrt(605.0 / 882.0);
    const auto f_1_3 = 1.0 / 3.0;
    const auto fs_243_98 = std::sqrt(243.0 / 98.0);
    const auto fs_8750_61347 = std::sqrt(8750.0 / 61347.0);
    const auto fs_1210_74529 = std::sqrt(1210.0 / 74529.0);
    const auto f_1_33 = 1.0 / 33.0;
    const auto fs_2_147 = std::sqrt(2.0 / 147.0);
    const auto fs_675_32 = std::sqrt(21.09375);
    const auto f_2 = 2.0;
    const auto fs_125_5577 = std::sqrt(125.0 / 5577.0);
    const auto f_4_13 = 4.0 / 13.0;
    const auto fs_4000_61347 = std::sqrt(4000.0 / 61347.0);
    const auto f_2_13 = 2.0 / 13.0;
    const auto fs_1_21 = std::sqrt(1.0 / 21.0);
    const auto fs_7000_61347 = std::sqrt(7000.0 / 61347.0);
    const auto fs_4_3549 = std::sqrt(4.0 / 3549.0);
    const auto fs_405_8 = std::sqrt(50.625);
    const auto fs_361_294 = std::sqrt(361.0 / 294.0);
    const auto fs_405_98 = std::sqrt(405.0 / 98.0);
    const auto fs_28000_184041 = std::sqrt(28000.0 / 184041.0);
    const auto fs_722_24843 = std::sqrt(722.0 / 24843.0);
    const auto fs_10_441 = std::sqrt(10.0 / 441.0);
    const auto fs_1125_32 = std::sqrt(35.15625);
    const auto f_9_2 = 4.5;
    const auto f_9 = 9.0;
    const auto f_10_7 = 10.0 / 7.0;
    const auto f_18_7 = 18.0 / 7.0;
    const auto f_175_429 = 175.0 / 429.0;
    const auto f_20_91 = 20.0 / 91.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_4_21 = 4.0 / 21.0;
    const auto f_15_2 = 7.5;

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p5, ph7_0, ph7_p1, ph7_p5, ph7_p6, ph7_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_0[k] = e_0 * (fs_189_16 * h3_p1 - fs_567_8 * r_2 * h1_p1) + e_1 * (fs_5_42 * h5_p1 - fs_7_3 * r_2 * h3_p1 + fs_81_14 * r_4 * h1_p1) + e_2 * (fs_25_368082 * h7_p1 + fs_175_858 * h7_p7 - fs_10_3549 * r_2 * h5_p1 + fs_7_363 * r_4 * h3_p1 - fs_2_63 * r_6 * h1_p1) + fs_1575_32 * e_3 * h1_p1;

        pc_1[k] = e_0 * (fs_567_16 * h3_0 - fs_567_16 * r_2 * h1_0) + e_1 * (fs_25_28 * h5_0 - fs_7 * r_2 * h3_0 + fs_81_28 * r_4 * h1_0) + e_2 * (fs_175_184041 * h7_0 + fs_50_429 * h7_p6 - fs_25_1183 * r_2 * h5_0 + fs_7_121 * r_4 * h3_0 - fs_1_63 * r_6 * h1_0) + fs_1575_64 * e_3 * h1_0;

        pc_2[k] = e_0 * (-fs_243_16 * h3_p1 + fs_81_32 * r_2 * h1_p1) + e_1 * (-fs_375_392 * h5_p1 + fs_25_28 * h5_p5 + fs_3 * r_2 * h3_p1 - fs_81_392 * r_4 * h1_p1) + e_2 * (-fs_350_184041 * h7_p1 + fs_350_5577 * h7_p5 + fs_375_16562 * r_2 * h5_p1 - fs_25_1183 * r_2 * h5_p5 - fs_3_121 * r_4 * h3_p1 + fs_1_882 * r_6 * h1_p1) - fs_225_128 * e_3 * h1_p1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph3_p2, ph5_m4, ph5_m3, ph5_m2, ph5_p2, ph5_p4, ph7_m4, ph7_m3, ph7_m2, ph7_p2, ph7_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];

        pc_3[k] = fs_135_16 * e_0 * h3_p2 + e_1 * (fs_125_84 * h5_p2 + fs_45_28 * h5_p4 - fs_5_3 * r_2 * h3_p2) + e_2 * (fs_350_61347 * h7_p2 + fs_175_5577 * h7_p4 - fs_125_3549 * r_2 * h5_p2 - fs_45_1183 * r_2 * h5_p4 + fs_5_363 * r_4 * h3_p2);

        pc_4[k] = -f_9_4 * e_0 * h3_m3 + e_1 * (-fs_25_7 * h5_m3 + f_1 * r_2 * h3_m3) + e_2 * (-fs_1750_61347 * h7_m3 + fs_100_1183 * r_2 * h5_m3 - f_1_11 * r_4 * h3_m3);

        pc_5[k] = fs_135_16 * e_0 * h3_m2 + e_1 * (-fs_45_28 * h5_m4 + fs_125_84 * h5_m2 - fs_5_3 * r_2 * h3_m2) + e_2 * (-fs_175_5577 * h7_m4 + fs_350_61347 * h7_m2 + fs_45_1183 * r_2 * h5_m4 - fs_125_3549 * r_2 * h5_m2 + fs_5_363 * r_4 * h3_m2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph3_m1, ph3_p2, ph5_m5, ph5_m1, ph5_p2, ph7_m7, ph7_m6, ph7_m5, ph7_m1, ph7_p2, ph7_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];

        pc_6[k] = e_0 * (-fs_243_16 * h3_m1 + fs_81_32 * r_2 * h1_m1) + e_1 * (-fs_25_28 * h5_m5 - fs_375_392 * h5_m1 + fs_3 * r_2 * h3_m1 - fs_81_392 * r_4 * h1_m1) + e_2 * (-fs_350_5577 * h7_m5 - fs_350_184041 * h7_m1 + fs_25_1183 * r_2 * h5_m5 + fs_375_16562 * r_2 * h5_m1 - fs_3_121 * r_4 * h3_m1 + fs_1_882 * r_6 * h1_m1) - fs_225_128 * e_3 * h1_m1;

        pc_7[k] = -fs_50_429 * e_2 * h7_m6;

        pc_8[k] = e_0 * (-fs_189_16 * h3_m1 + fs_567_8 * r_2 * h1_m1) + e_1 * (-fs_5_42 * h5_m1 + fs_7_3 * r_2 * h3_m1 - fs_81_14 * r_4 * h1_m1) + e_2 * (-fs_175_858 * h7_m7 - fs_25_368082 * h7_m1 + fs_10_3549 * r_2 * h5_m1 - fs_7_363 * r_4 * h3_m1 + fs_2_63 * r_6 * h1_m1) - fs_1575_32 * e_3 * h1_m1;

        pc_9[k] = -fs_315_16 * e_0 * h3_p2 + e_1 * (-fs_5_9 * h5_p2 + fs_35_9 * r_2 * h3_p2) + e_2 * (-fs_25_40898 * h7_p2 + fs_25_286 * h7_p6 + fs_20_1521 * r_2 * h5_p2 - fs_35_1089 * r_4 * h3_p2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p4, ph5_p5, ph7_0, ph7_p1, ph7_p4, ph7_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_10[k] = e_0 * (-fs_63_16 * h3_p1 - fs_1701_32 * r_2 * h1_p1) + e_1 * (-fs_605_504 * h5_p1 - fs_25_12 * h5_p5 + fs_7_9 * r_2 * h3_p1 + fs_243_56 * r_4 * h1_p1) + e_2 * (-fs_200_61347 * h7_p1 + fs_200_1859 * h7_p5 + fs_605_21294 * r_2 * h5_p1 + fs_25_507 * r_2 * h5_p5 - fs_7_1089 * r_4 * h3_p1 - fs_1_42 * r_6 * h1_p1) + fs_4725_128 * e_3 * h1_p1;

        pc_11[k] = e_0 * (fs_27_16 * h3_0 - fs_243_4 * r_2 * h1_0) + e_1 * (-fs_400_147 * h5_0 - fs_20_21 * h5_p4 - fs_1_3 * r_2 * h3_0 + fs_243_49 * r_4 * h1_0) + e_2 * (-fs_1225_61347 * h7_0 + fs_175_1859 * h7_p4 + fs_1600_24843 * r_2 * h5_0 + fs_80_3549 * r_2 * h5_p4 + fs_1_363 * r_4 * h3_0 - fs_4_147 * r_6 * h1_0) + fs_675_16 * e_3 * h1_0;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_p1, ph3_m2, ph3_p1, ph3_p3, ph5_m2, ph5_p1, ph5_p3, ph7_m2, ph7_p1, ph7_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];

        pc_12[k] = e_0 * (-f_3 * h3_p1 - fs_135_16 * h3_p3 + fs_243_32 * r_2 * h1_p1) + e_1 * (fs_3125_3528 * h5_p1 - fs_5_84 * h5_p3 + f_4_3 * r_2 * h3_p1 + fs_5_3 * r_2 * h3_p3 - fs_243_392 * r_4 * h1_p1) + e_2 * (fs_1400_61347 * h7_p1 + fs_1400_20449 * h7_p3 - fs_3125_149058 * r_2 * h5_p1 + fs_5_3549 * r_2 * h5_p3 - f_4_33 * r_4 * h3_p1 - fs_5_363 * r_4 * h3_p3 + fs_1_294 * r_6 * h1_p1) - fs_675_128 * e_3 * h1_p1;

        pc_13[k] = f_21_4 * e_0 * h3_m2 + e_1 * (-fs_25_63 * h5_m2 - f_7_3 * r_2 * h3_m2) + e_2 * (-fs_1750_20449 * h7_m2 + fs_100_10647 * r_2 * h5_m2 + f_7_33 * r_4 * h3_m2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph3_m3, ph3_m1, ph5_m5, ph5_m4, ph5_m3, ph5_m1, ph7_m5, ph7_m4, ph7_m3, ph7_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];

        pc_14[k] = e_0 * (fs_135_16 * h3_m3 - f_3 * h3_m1 + fs_243_32 * r_2 * h1_m1) + e_1 * (fs_5_84 * h5_m3 + fs_3125_3528 * h5_m1 - fs_5_3 * r_2 * h3_m3 + f_4_3 * r_2 * h3_m1 - fs_243_392 * r_4 * h1_m1) + e_2 * (-fs_1400_20449 * h7_m3 + fs_1400_61347 * h7_m1 - fs_5_3549 * r_2 * h5_m3 - fs_3125_149058 * r_2 * h5_m1 + fs_5_363 * r_4 * h3_m3 - f_4_33 * r_4 * h3_m1 + fs_1_294 * r_6 * h1_m1) - fs_675_128 * e_3 * h1_m1;

        pc_15[k] = fs_20_21 * e_1 * h5_m4 + e_2 * (-fs_175_1859 * h7_m4 - fs_80_3549 * r_2 * h5_m4);

        pc_16[k] = e_0 * (fs_63_16 * h3_m1 + fs_1701_32 * r_2 * h1_m1) + e_1 * (fs_25_12 * h5_m5 + fs_605_504 * h5_m1 - fs_7_9 * r_2 * h3_m1 - fs_243_56 * r_4 * h1_m1) + e_2 * (-fs_200_1859 * h7_m5 + fs_200_61347 * h7_m1 - fs_25_507 * r_2 * h5_m5 - fs_605_21294 * r_2 * h5_m1 + fs_7_1089 * r_4 * h3_m1 + fs_1_42 * r_6 * h1_m1) - fs_4725_128 * e_3 * h1_m1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_p2, ph3_p3, ph5_m2, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ph7_m6, ph7_m2, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_17[k] = fs_315_16 * e_0 * h3_m2 + e_1 * (fs_5_9 * h5_m2 - fs_35_9 * r_2 * h3_m2) + e_2 * (-fs_25_286 * h7_m6 + fs_25_40898 * h7_m2 - fs_20_1521 * r_2 * h5_m2 + fs_35_1089 * r_4 * h3_m2);

        pc_18[k] = fs_189_16 * e_0 * h3_p3 + e_1 * (fs_4_3 * h5_p3 + fs_5_3 * h5_p5 - fs_7_3 * r_2 * h3_p3) + e_2 * (fs_125_40898 * h7_p3 + fs_125_3718 * h7_p5 - fs_16_507 * r_2 * h5_p3 - fs_20_507 * r_2 * h5_p5 + fs_7_363 * r_4 * h3_p3);

        pc_19[k] = -fs_63_16 * e_0 * h3_p2 + e_1 * (f_7_6 * h5_p2 - fs_1_12 * h5_p4 + fs_7_9 * r_2 * h3_p2) + e_2 * (fs_250_20449 * h7_p2 + fs_125_1859 * h7_p4 - f_7_39 * r_2 * h5_p2 + fs_1_507 * r_2 * h5_p4 - fs_7_1089 * r_4 * h3_p2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];

        pc_20[k] = e_0 * (-fs_45_4 * h3_p1 + fs_243_16 * h3_p3 - fs_1215_32 * r_2 * h1_p1) + e_1 * (fs_1849_3528 * h5_p1 - fs_27_28 * h5_p3 + fs_20_9 * r_2 * h3_p1 - fs_3 * r_2 * h3_p3 + fs_1215_392 * r_4 * h1_p1) + e_2 * (fs_1750_61347 * h7_p1 + fs_1750_20449 * h7_p3 - fs_1849_149058 * r_2 * h5_p1 + fs_27_1183 * r_2 * h5_p3 - fs_20_1089 * r_4 * h3_p1 + fs_3_121 * r_4 * h3_p3 - fs_5_294 * r_6 * h1_p1) + fs_3375_128 * e_3 * h1_p1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph1_0, ph3_m1, ph3_0, ph3_p2, ph5_m1, ph5_0, ph5_p2, ph7_m1, ph7_0, ph7_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];

        pc_21[k] = e_0 * (-fs_135_16 * h3_0 + f_3 * h3_p2 - fs_1215_16 * r_2 * h1_0) + e_1 * (fs_5_588 * h5_0 - fs_64_63 * h5_p2 + fs_5_3 * r_2 * h3_0 - f_4_3 * r_2 * h3_p2 + fs_1215_196 * r_4 * h1_0) + e_2 * (fs_6125_61347 * h7_0 + fs_1750_20449 * h7_p2 - fs_5_24843 * r_2 * h5_0 + fs_256_10647 * r_2 * h5_p2 - fs_5_363 * r_4 * h3_0 + f_4_33 * r_4 * h3_p2 - fs_5_147 * r_6 * h1_0) + fs_3375_64 * e_3 * h1_0;

        pc_22[k] = e_0 * (-f_3_4 * h3_m1 + fs_243_8 * r_2 * h1_m1) + e_1 * (fs_605_882 * h5_m1 + f_1_3 * r_2 * h3_m1 - fs_243_98 * r_4 * h1_m1) + e_2 * (-fs_8750_61347 * h7_m1 - fs_1210_74529 * r_2 * h5_m1 - f_1_33 * r_4 * h3_m1 + fs_2_147 * r_6 * h1_m1) - fs_675_32 * e_3 * h1_m1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph3_m3, ph3_m2, ph3_m1, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_23[k] = -f_3 * e_0 * h3_m2 + e_1 * (fs_64_63 * h5_m2 + f_4_3 * r_2 * h3_m2) + e_2 * (-fs_1750_20449 * h7_m2 - fs_256_10647 * r_2 * h5_m2 - f_4_33 * r_4 * h3_m2);

        pc_24[k] = e_0 * (-fs_243_16 * h3_m3 + fs_45_4 * h3_m1 + fs_1215_32 * r_2 * h1_m1) + e_1 * (fs_27_28 * h5_m3 - fs_1849_3528 * h5_m1 + fs_3 * r_2 * h3_m3 - fs_20_9 * r_2 * h3_m1 - fs_1215_392 * r_4 * h1_m1) + e_2 * (-fs_1750_20449 * h7_m3 - fs_1750_61347 * h7_m1 - fs_27_1183 * r_2 * h5_m3 + fs_1849_149058 * r_2 * h5_m1 - fs_3_121 * r_4 * h3_m3 + fs_20_1089 * r_4 * h3_m1 + fs_5_294 * r_6 * h1_m1) - fs_3375_128 * e_3 * h1_m1;

        pc_25[k] = fs_63_16 * e_0 * h3_m2 + e_1 * (fs_1_12 * h5_m4 - f_7_6 * h5_m2 - fs_7_9 * r_2 * h3_m2) + e_2 * (-fs_125_1859 * h7_m4 - fs_250_20449 * h7_m2 - fs_1_507 * r_2 * h5_m4 + f_7_39 * r_2 * h5_m2 + fs_7_1089 * r_4 * h3_m2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];

        pc_26[k] = -fs_189_16 * e_0 * h3_m3 + e_1 * (-fs_5_3 * h5_m5 - fs_4_3 * h5_m3 + fs_7_3 * r_2 * h3_m3) + e_2 * (-fs_125_3718 * h7_m5 - fs_125_40898 * h7_m3 + fs_20_507 * r_2 * h5_m5 + fs_16_507 * r_2 * h5_m3 - fs_7_363 * r_4 * h3_m3);

        pc_27[k] = -f_2 * e_1 * h5_m4 + e_2 * (-fs_125_5577 * h7_m4 + f_4_13 * r_2 * h5_m4);

        pc_28[k] = fs_567_16 * e_0 * h3_m3 + e_1 * (-f_1 * h5_m3 - fs_7 * r_2 * h3_m3) + e_2 * (-fs_4000_61347 * h7_m3 + f_2_13 * r_2 * h5_m3 + fs_7_121 * r_4 * h3_m3);

        pc_29[k] = fs_27_16 * e_0 * h3_m2 + e_1 * (fs_1_21 * h5_m2 - fs_1_3 * r_2 * h3_m2) + e_2 * (-fs_7000_61347 * h7_m2 - fs_4_3549 * r_2 * h5_m2 + fs_1_363 * r_4 * h3_m2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph1_0, ph1_p1, ph3_m1, ph3_0, ph3_p1, ph5_m1, ph5_0, ph5_p1, ph7_m1, ph7_0, ph7_p1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];

        pc_30[k] = e_0 * (-fs_135_16 * h3_m1 - fs_405_8 * r_2 * h1_m1) + e_1 * (fs_361_294 * h5_m1 + fs_5_3 * r_2 * h3_m1 + fs_405_98 * r_4 * h1_m1) + e_2 * (-fs_28000_184041 * h7_m1 - fs_722_24843 * r_2 * h5_m1 - fs_5_363 * r_4 * h3_m1 - fs_10_441 * r_6 * h1_m1) + fs_1125_32 * e_3 * h1_m1;

        pc_31[k] = e_0 * (-f_9_2 * h3_0 - f_9 * r_2 * h1_0) + e_1 * (f_10_7 * h5_0 + f_2 * r_2 * h3_0 + f_18_7 * r_4 * h1_0) + e_2 * (-f_175_429 * h7_0 - f_20_91 * r_2 * h5_0 - f_2_11 * r_4 * h3_0 - f_4_21 * r_6 * h1_0) + f_15_2 * e_3 * h1_0;

        pc_32[k] = e_0 * (-fs_135_16 * h3_p1 - fs_405_8 * r_2 * h1_p1) + e_1 * (fs_361_294 * h5_p1 + fs_5_3 * r_2 * h3_p1 + fs_405_98 * r_4 * h1_p1) + e_2 * (-fs_28000_184041 * h7_p1 - fs_722_24843 * r_2 * h5_p1 - fs_5_363 * r_4 * h3_p1 - fs_10_441 * r_6 * h1_p1) + fs_1125_32 * e_3 * h1_p1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p2, ph3_p3, ph5_m5, ph5_m3, ph5_p2, ph5_p3, ph5_p4, ph7_m5, ph7_m3, ph7_p2, ph7_p3, ph7_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_33[k] = fs_27_16 * e_0 * h3_p2 + e_1 * (fs_1_21 * h5_p2 - fs_1_3 * r_2 * h3_p2) + e_2 * (-fs_7000_61347 * h7_p2 - fs_4_3549 * r_2 * h5_p2 + fs_1_363 * r_4 * h3_p2);

        pc_34[k] = fs_567_16 * e_0 * h3_p3 + e_1 * (-f_1 * h5_p3 - fs_7 * r_2 * h3_p3) + e_2 * (-fs_4000_61347 * h7_p3 + f_2_13 * r_2 * h5_p3 + fs_7_121 * r_4 * h3_p3);

        pc_35[k] = -f_2 * e_1 * h5_p4 + e_2 * (-fs_125_5577 * h7_p4 + f_4_13 * r_2 * h5_p4);

        pc_36[k] = fs_189_16 * e_0 * h3_m3 + e_1 * (-fs_5_3 * h5_m5 + fs_4_3 * h5_m3 - fs_7_3 * r_2 * h3_m3) + e_2 * (-fs_125_3718 * h7_m5 + fs_125_40898 * h7_m3 + fs_20_507 * r_2 * h5_m5 - fs_16_507 * r_2 * h5_m3 + fs_7_363 * r_4 * h3_m3);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph3_m3, ph3_m2, ph3_m1, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_37[k] = -fs_63_16 * e_0 * h3_m2 + e_1 * (fs_1_12 * h5_m4 + f_7_6 * h5_m2 + fs_7_9 * r_2 * h3_m2) + e_2 * (-fs_125_1859 * h7_m4 + fs_250_20449 * h7_m2 - fs_1_507 * r_2 * h5_m4 - f_7_39 * r_2 * h5_m2 - fs_7_1089 * r_4 * h3_m2);

        pc_38[k] = e_0 * (-fs_243_16 * h3_m3 - fs_45_4 * h3_m1 - fs_1215_32 * r_2 * h1_m1) + e_1 * (fs_27_28 * h5_m3 + fs_1849_3528 * h5_m1 + fs_3 * r_2 * h3_m3 + fs_20_9 * r_2 * h3_m1 + fs_1215_392 * r_4 * h1_m1) + e_2 * (-fs_1750_20449 * h7_m3 + fs_1750_61347 * h7_m1 - fs_27_1183 * r_2 * h5_m3 - fs_1849_149058 * r_2 * h5_m1 - fs_3_121 * r_4 * h3_m3 - fs_20_1089 * r_4 * h3_m1 - fs_5_294 * r_6 * h1_m1) + fs_3375_128 * e_3 * h1_m1;

        pc_39[k] = -f_3 * e_0 * h3_m2 + e_1 * (fs_64_63 * h5_m2 + f_4_3 * r_2 * h3_m2) + e_2 * (-fs_1750_20449 * h7_m2 - fs_256_10647 * r_2 * h5_m2 - f_4_33 * r_4 * h3_m2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_0, ph1_p1, ph3_0, ph3_p1, ph3_p2, ph5_0, ph5_p1, ph5_p2, ph7_0, ph7_p1, ph7_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_40[k] = e_0 * (-f_3_4 * h3_p1 + fs_243_8 * r_2 * h1_p1) + e_1 * (fs_605_882 * h5_p1 + f_1_3 * r_2 * h3_p1 - fs_243_98 * r_4 * h1_p1) + e_2 * (-fs_8750_61347 * h7_p1 - fs_1210_74529 * r_2 * h5_p1 - f_1_33 * r_4 * h3_p1 + fs_2_147 * r_6 * h1_p1) - fs_675_32 * e_3 * h1_p1;

        pc_41[k] = e_0 * (-fs_135_16 * h3_0 - f_3 * h3_p2 - fs_1215_16 * r_2 * h1_0) + e_1 * (fs_5_588 * h5_0 + fs_64_63 * h5_p2 + fs_5_3 * r_2 * h3_0 + f_4_3 * r_2 * h3_p2 + fs_1215_196 * r_4 * h1_0) + e_2 * (fs_6125_61347 * h7_0 - fs_1750_20449 * h7_p2 - fs_5_24843 * r_2 * h5_0 - fs_256_10647 * r_2 * h5_p2 - fs_5_363 * r_4 * h3_0 - f_4_33 * r_4 * h3_p2 - fs_5_147 * r_6 * h1_0) + fs_3375_64 * e_3 * h1_0;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_p1, ph3_p1, ph3_p2, ph3_p3, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph7_p1, ph7_p2, ph7_p3, ph7_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_42[k] = e_0 * (-fs_45_4 * h3_p1 - fs_243_16 * h3_p3 - fs_1215_32 * r_2 * h1_p1) + e_1 * (fs_1849_3528 * h5_p1 + fs_27_28 * h5_p3 + fs_20_9 * r_2 * h3_p1 + fs_3 * r_2 * h3_p3 + fs_1215_392 * r_4 * h1_p1) + e_2 * (fs_1750_61347 * h7_p1 - fs_1750_20449 * h7_p3 - fs_1849_149058 * r_2 * h5_p1 - fs_27_1183 * r_2 * h5_p3 - fs_20_1089 * r_4 * h3_p1 - fs_3_121 * r_4 * h3_p3 - fs_5_294 * r_6 * h1_p1) + fs_3375_128 * e_3 * h1_p1;

        pc_43[k] = -fs_63_16 * e_0 * h3_p2 + e_1 * (f_7_6 * h5_p2 + fs_1_12 * h5_p4 + fs_7_9 * r_2 * h3_p2) + e_2 * (fs_250_20449 * h7_p2 - fs_125_1859 * h7_p4 - f_7_39 * r_2 * h5_p2 - fs_1_507 * r_2 * h5_p4 - fs_7_1089 * r_4 * h3_p2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_p3, ph5_m2, ph5_p3, ph5_p5, ph7_m6, ph7_m2, ph7_p3, ph7_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];

        pc_44[k] = fs_189_16 * e_0 * h3_p3 + e_1 * (fs_4_3 * h5_p3 - fs_5_3 * h5_p5 - fs_7_3 * r_2 * h3_p3) + e_2 * (fs_125_40898 * h7_p3 - fs_125_3718 * h7_p5 - fs_16_507 * r_2 * h5_p3 + fs_20_507 * r_2 * h5_p5 + fs_7_363 * r_4 * h3_p3);

        pc_45[k] = -fs_315_16 * e_0 * h3_m2 + e_1 * (-fs_5_9 * h5_m2 + fs_35_9 * r_2 * h3_m2) + e_2 * (-fs_25_286 * h7_m6 - fs_25_40898 * h7_m2 + fs_20_1521 * r_2 * h5_m2 - fs_35_1089 * r_4 * h3_m2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph3_m3, ph3_m1, ph5_m5, ph5_m4, ph5_m3, ph5_m1, ph7_m5, ph7_m4, ph7_m3, ph7_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];

        pc_46[k] = e_0 * (-fs_63_16 * h3_m1 - fs_1701_32 * r_2 * h1_m1) + e_1 * (fs_25_12 * h5_m5 - fs_605_504 * h5_m1 + fs_7_9 * r_2 * h3_m1 + fs_243_56 * r_4 * h1_m1) + e_2 * (-fs_200_1859 * h7_m5 - fs_200_61347 * h7_m1 - fs_25_507 * r_2 * h5_m5 + fs_605_21294 * r_2 * h5_m1 - fs_7_1089 * r_4 * h3_m1 - fs_1_42 * r_6 * h1_m1) + fs_4725_128 * e_3 * h1_m1;

        pc_47[k] = fs_20_21 * e_1 * h5_m4 + e_2 * (-fs_175_1859 * h7_m4 - fs_80_3549 * r_2 * h5_m4);

        pc_48[k] = e_0 * (fs_135_16 * h3_m3 + f_3 * h3_m1 - fs_243_32 * r_2 * h1_m1) + e_1 * (fs_5_84 * h5_m3 - fs_3125_3528 * h5_m1 - fs_5_3 * r_2 * h3_m3 - f_4_3 * r_2 * h3_m1 + fs_243_392 * r_4 * h1_m1) + e_2 * (-fs_1400_20449 * h7_m3 - fs_1400_61347 * h7_m1 - fs_5_3549 * r_2 * h5_m3 + fs_3125_149058 * r_2 * h5_m1 + fs_5_363 * r_4 * h3_m3 + f_4_33 * r_4 * h3_m1 - fs_1_294 * r_6 * h1_m1) + fs_675_128 * e_3 * h1_m1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_p1, ph3_p1, ph3_p2, ph3_p3, ph5_p1, ph5_p2, ph5_p3, ph7_p1, ph7_p2, ph7_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];

        pc_49[k] = f_21_4 * e_0 * h3_p2 + e_1 * (-fs_25_63 * h5_p2 - f_7_3 * r_2 * h3_p2) + e_2 * (-fs_1750_20449 * h7_p2 + fs_100_10647 * r_2 * h5_p2 + f_7_33 * r_4 * h3_p2);

        pc_50[k] = e_0 * (-f_3 * h3_p1 + fs_135_16 * h3_p3 + fs_243_32 * r_2 * h1_p1) + e_1 * (fs_3125_3528 * h5_p1 + fs_5_84 * h5_p3 + f_4_3 * r_2 * h3_p1 - fs_5_3 * r_2 * h3_p3 - fs_243_392 * r_4 * h1_p1) + e_2 * (fs_1400_61347 * h7_p1 - fs_1400_20449 * h7_p3 - fs_3125_149058 * r_2 * h5_p1 - fs_5_3549 * r_2 * h5_p3 - f_4_33 * r_4 * h3_p1 + fs_5_363 * r_4 * h3_p3 + fs_1_294 * r_6 * h1_p1) - fs_675_128 * e_3 * h1_p1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p4, ph5_p5, ph7_0, ph7_p1, ph7_p4, ph7_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_51[k] = e_0 * (fs_27_16 * h3_0 - fs_243_4 * r_2 * h1_0) + e_1 * (-fs_400_147 * h5_0 + fs_20_21 * h5_p4 - fs_1_3 * r_2 * h3_0 + fs_243_49 * r_4 * h1_0) + e_2 * (-fs_1225_61347 * h7_0 - fs_175_1859 * h7_p4 + fs_1600_24843 * r_2 * h5_0 - fs_80_3549 * r_2 * h5_p4 + fs_1_363 * r_4 * h3_0 - fs_4_147 * r_6 * h1_0) + fs_675_16 * e_3 * h1_0;

        pc_52[k] = e_0 * (-fs_63_16 * h3_p1 - fs_1701_32 * r_2 * h1_p1) + e_1 * (-fs_605_504 * h5_p1 + fs_25_12 * h5_p5 + fs_7_9 * r_2 * h3_p1 + fs_243_56 * r_4 * h1_p1) + e_2 * (-fs_200_61347 * h7_p1 - fs_200_1859 * h7_p5 + fs_605_21294 * r_2 * h5_p1 - fs_25_507 * r_2 * h5_p5 - fs_7_1089 * r_4 * h3_p1 - fs_1_42 * r_6 * h1_p1) + fs_4725_128 * e_3 * h1_p1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_m1, ph3_m1, ph3_p2, ph5_m5, ph5_m1, ph5_p2, ph7_m7, ph7_m6, ph7_m5, ph7_m1, ph7_p2, ph7_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];

        pc_53[k] = -fs_315_16 * e_0 * h3_p2 + e_1 * (-fs_5_9 * h5_p2 + fs_35_9 * r_2 * h3_p2) + e_2 * (-fs_25_40898 * h7_p2 - fs_25_286 * h7_p6 + fs_20_1521 * r_2 * h5_p2 - fs_35_1089 * r_4 * h3_p2);

        pc_54[k] = e_0 * (fs_189_16 * h3_m1 - fs_567_8 * r_2 * h1_m1) + e_1 * (fs_5_42 * h5_m1 - fs_7_3 * r_2 * h3_m1 + fs_81_14 * r_4 * h1_m1) + e_2 * (-fs_175_858 * h7_m7 + fs_25_368082 * h7_m1 - fs_10_3549 * r_2 * h5_m1 + fs_7_363 * r_4 * h3_m1 - fs_2_63 * r_6 * h1_m1) + fs_1575_32 * e_3 * h1_m1;

        pc_55[k] = -fs_50_429 * e_2 * h7_m6;

        pc_56[k] = e_0 * (fs_243_16 * h3_m1 - fs_81_32 * r_2 * h1_m1) + e_1 * (-fs_25_28 * h5_m5 + fs_375_392 * h5_m1 - fs_3 * r_2 * h3_m1 + fs_81_392 * r_4 * h1_m1) + e_2 * (-fs_350_5577 * h7_m5 + fs_350_184041 * h7_m1 + fs_25_1183 * r_2 * h5_m5 - fs_375_16562 * r_2 * h5_m1 + fs_3_121 * r_4 * h3_m1 - fs_1_882 * r_6 * h1_m1) + fs_225_128 * e_3 * h1_m1;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_p2, ph3_p3, ph5_m4, ph5_m2, ph5_p2, ph5_p3, ph5_p4, ph7_m4, ph7_m2, ph7_p2, ph7_p3, ph7_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_57[k] = -fs_135_16 * e_0 * h3_m2 + e_1 * (-fs_45_28 * h5_m4 - fs_125_84 * h5_m2 + fs_5_3 * r_2 * h3_m2) + e_2 * (-fs_175_5577 * h7_m4 - fs_350_61347 * h7_m2 + fs_45_1183 * r_2 * h5_m4 + fs_125_3549 * r_2 * h5_m2 - fs_5_363 * r_4 * h3_m2);

        pc_58[k] = -f_9_4 * e_0 * h3_p3 + e_1 * (-fs_25_7 * h5_p3 + f_1 * r_2 * h3_p3) + e_2 * (-fs_1750_61347 * h7_p3 + fs_100_1183 * r_2 * h5_p3 - f_1_11 * r_4 * h3_p3);

        pc_59[k] = fs_135_16 * e_0 * h3_p2 + e_1 * (fs_125_84 * h5_p2 - fs_45_28 * h5_p4 - fs_5_3 * r_2 * h3_p2) + e_2 * (fs_350_61347 * h7_p2 - fs_175_5577 * h7_p4 - fs_125_3549 * r_2 * h5_p2 + fs_45_1183 * r_2 * h5_p4 + fs_5_363 * r_4 * h3_p2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 63 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p5, ph7_0, ph7_p1, ph7_p5, ph7_p6, ph7_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_60[k] = e_0 * (-fs_243_16 * h3_p1 + fs_81_32 * r_2 * h1_p1) + e_1 * (-fs_375_392 * h5_p1 - fs_25_28 * h5_p5 + fs_3 * r_2 * h3_p1 - fs_81_392 * r_4 * h1_p1) + e_2 * (-fs_350_184041 * h7_p1 - fs_350_5577 * h7_p5 + fs_375_16562 * r_2 * h5_p1 + fs_25_1183 * r_2 * h5_p5 - fs_3_121 * r_4 * h3_p1 + fs_1_882 * r_6 * h1_p1) - fs_225_128 * e_3 * h1_p1;

        pc_61[k] = e_0 * (fs_567_16 * h3_0 - fs_567_16 * r_2 * h1_0) + e_1 * (fs_25_28 * h5_0 - fs_7 * r_2 * h3_0 + fs_81_28 * r_4 * h1_0) + e_2 * (fs_175_184041 * h7_0 - fs_50_429 * h7_p6 - fs_25_1183 * r_2 * h5_0 + fs_7_121 * r_4 * h3_0 - fs_1_63 * r_6 * h1_0) + fs_1575_64 * e_3 * h1_0;

        pc_62[k] = e_0 * (fs_189_16 * h3_p1 - fs_567_8 * r_2 * h1_p1) + e_1 * (fs_5_42 * h5_p1 - fs_7_3 * r_2 * h3_p1 + fs_81_14 * r_4 * h1_p1) + e_2 * (fs_25_368082 * h7_p1 - fs_175_858 * h7_p7 - fs_10_3549 * r_2 * h5_p1 + fs_7_363 * r_4 * h3_p1 - fs_2_63 * r_6 * h1_p1) + fs_1575_32 * e_3 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[63] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62};

    for (size_t m = 0; m < 63; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
