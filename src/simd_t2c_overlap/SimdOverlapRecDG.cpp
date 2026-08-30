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



#include "SimdOverlapRecDG.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_dg_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDG.compute_dg_overlap: Basis functions must be of angular momenta two and four"));
    }

    if (harmonics.size() < 6)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDG.compute_dg_overlap: Harmonics must reach angular momentum six"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDG.compute_dg_overlap: Number of values exceeds number of atom pairs"));
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

    const auto nmax = dimensions.back();

    if (nmax == 0)
    {
        std::fill(values, values + 45 * nvalues, 0.0);

        return;
    }

    // NOTE: the first three rows accumulate the contracted prefactors of the terms,
    // and the remaining 45 rows hold the integrals of the combinations of angular
    // components.

    auto buffer = CSimdMatrix(48, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);

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

            const auto f_0 = fbase * aexp * aexp * fmu / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * aexp * aexp * fmu * fmu / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * aexp * aexp / fexp / fexp / fexp / fexp;

            // NOTE: the exponential depends on the pair of primitives alone, so it is
            // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fss = std::exp(-fmu * ab_2[k]);

                pe_0[k] += f_0 * fss;
                pe_1[k] += f_1 * fss;
                pe_2[k] += f_2 * fss;
            }
        }
    }

    // NOTE: the geometry of a term is a solid harmonic of the vector between the
    // atoms times a power of their squared distance.

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

    auto *pc_0 = buffer.data(3);
    auto *pc_1 = buffer.data(4);
    auto *pc_2 = buffer.data(5);
    auto *pc_3 = buffer.data(6);
    auto *pc_4 = buffer.data(7);
    auto *pc_5 = buffer.data(8);
    auto *pc_6 = buffer.data(9);
    auto *pc_7 = buffer.data(10);
    auto *pc_8 = buffer.data(11);
    auto *pc_9 = buffer.data(12);
    auto *pc_10 = buffer.data(13);
    auto *pc_11 = buffer.data(14);
    auto *pc_12 = buffer.data(15);
    auto *pc_13 = buffer.data(16);
    auto *pc_14 = buffer.data(17);
    auto *pc_15 = buffer.data(18);
    auto *pc_16 = buffer.data(19);
    auto *pc_17 = buffer.data(20);
    auto *pc_18 = buffer.data(21);
    auto *pc_19 = buffer.data(22);
    auto *pc_20 = buffer.data(23);
    auto *pc_21 = buffer.data(24);
    auto *pc_22 = buffer.data(25);
    auto *pc_23 = buffer.data(26);
    auto *pc_24 = buffer.data(27);
    auto *pc_25 = buffer.data(28);
    auto *pc_26 = buffer.data(29);
    auto *pc_27 = buffer.data(30);
    auto *pc_28 = buffer.data(31);
    auto *pc_29 = buffer.data(32);
    auto *pc_30 = buffer.data(33);
    auto *pc_31 = buffer.data(34);
    auto *pc_32 = buffer.data(35);
    auto *pc_33 = buffer.data(36);
    auto *pc_34 = buffer.data(37);
    auto *pc_35 = buffer.data(38);
    auto *pc_36 = buffer.data(39);
    auto *pc_37 = buffer.data(40);
    auto *pc_38 = buffer.data(41);
    auto *pc_39 = buffer.data(42);
    auto *pc_40 = buffer.data(43);
    auto *pc_41 = buffer.data(44);
    auto *pc_42 = buffer.data(45);
    auto *pc_43 = buffer.data(46);
    auto *pc_44 = buffer.data(47);

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_3_7 = std::sqrt(3.0 / 7.0);
    const auto fs_45_7 = std::sqrt(45.0 / 7.0);
    const auto fs_1_2178 = std::sqrt(1.0 / 2178.0);
    const auto fs_5_22 = std::sqrt(5.0 / 22.0);
    const auto fs_12_847 = std::sqrt(12.0 / 847.0);
    const auto fs_5_63 = std::sqrt(5.0 / 63.0);
    const auto fs_315_16 = std::sqrt(19.6875);
    const auto fs_27_28 = std::sqrt(27.0 / 28.0);
    const auto fs_45_14 = std::sqrt(45.0 / 14.0);
    const auto fs_5_2178 = std::sqrt(5.0 / 2178.0);
    const auto fs_5_33 = std::sqrt(5.0 / 33.0);
    const auto fs_27_847 = std::sqrt(27.0 / 847.0);
    const auto fs_5_126 = std::sqrt(5.0 / 126.0);
    const auto fs_315_32 = std::sqrt(9.84375);
    const auto fs_135_49 = std::sqrt(135.0 / 49.0);
    const auto fs_5_363 = std::sqrt(5.0 / 363.0);
    const auto fs_35_363 = std::sqrt(35.0 / 363.0);
    const auto fs_540_5929 = std::sqrt(540.0 / 5929.0);
    const auto fs_5_147 = std::sqrt(5.0 / 147.0);
    const auto fs_135_16 = std::sqrt(8.4375);
    const auto fs_75_49 = std::sqrt(75.0 / 49.0);
    const auto fs_45_98 = std::sqrt(45.0 / 98.0);
    const auto fs_35_2178 = std::sqrt(35.0 / 2178.0);
    const auto fs_7_121 = std::sqrt(7.0 / 121.0);
    const auto fs_300_5929 = std::sqrt(300.0 / 5929.0);
    const auto fs_5_882 = std::sqrt(5.0 / 882.0);
    const auto fs_45_32 = std::sqrt(1.40625);
    const auto f_3_7 = 3.0 / 7.0;
    const auto fs_70_1089 = std::sqrt(70.0 / 1089.0);
    const auto f_1_21 = 1.0 / 21.0;
    const auto f_3_4 = 0.75;
    const auto fs_3_2 = std::sqrt(1.5);
    const auto fs_1_242 = std::sqrt(1.0 / 242.0);
    const auto fs_5_66 = std::sqrt(5.0 / 66.0);
    const auto fs_6_121 = std::sqrt(6.0 / 121.0);
    const auto fs_75_56 = std::sqrt(75.0 / 56.0);
    const auto f_4_33 = 4.0 / 33.0;
    const auto fs_40_363 = std::sqrt(40.0 / 363.0);
    const auto fs_75_1694 = std::sqrt(75.0 / 1694.0);
    const auto fs_243_392 = std::sqrt(243.0 / 392.0);
    const auto fs_180_49 = std::sqrt(180.0 / 49.0);
    const auto fs_35_1089 = std::sqrt(35.0 / 1089.0);
    const auto fs_14_121 = std::sqrt(14.0 / 121.0);
    const auto fs_243_11858 = std::sqrt(243.0 / 11858.0);
    const auto fs_20_441 = std::sqrt(20.0 / 441.0);
    const auto fs_45_4 = std::sqrt(11.25);
    const auto fs_15_98 = std::sqrt(15.0 / 98.0);
    const auto fs_270_49 = std::sqrt(270.0 / 49.0);
    const auto fs_112_1089 = std::sqrt(112.0 / 1089.0);
    const auto fs_30_5929 = std::sqrt(30.0 / 5929.0);
    const auto fs_10_147 = std::sqrt(10.0 / 147.0);
    const auto fs_135_8 = std::sqrt(16.875);
    const auto f_12_7 = 12.0 / 7.0;
    const auto fs_175_1089 = std::sqrt(175.0 / 1089.0);
    const auto f_4_21 = 4.0 / 21.0;
    const auto f_3 = 3.0;
    const auto f_2 = 2.0;
    const auto fs_5_121 = std::sqrt(5.0 / 121.0);
    const auto f_4_11 = 4.0 / 11.0;
    const auto f_1_2 = 0.5;
    const auto fs_12_121 = std::sqrt(12.0 / 121.0);
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_4_7 = 4.0 / 7.0;
    const auto fs_56_363 = std::sqrt(56.0 / 363.0);
    const auto f_8_77 = 8.0 / 77.0;
    const auto f_17_14 = 17.0 / 14.0;
    const auto fs_70_363 = std::sqrt(70.0 / 363.0);
    const auto f_17_77 = 17.0 / 77.0;
    const auto f_10_7 = 10.0 / 7.0;
    const auto f_18_7 = 18.0 / 7.0;
    const auto f_5_11 = 5.0 / 11.0;
    const auto f_20_77 = 20.0 / 77.0;
    const auto f_2_7 = 2.0 / 7.0;
    const auto f_9_2 = 4.5;

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p4, ph6_0, ph6_p1, ph6_p2, ph6_p4, ph6_p5, ph6_p6, ab_2, pc_0, pc_1, pc_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];

        pc_0[k] = e_0 * (fs_3_7 * h4_p2 - fs_45_7 * r_2 * h2_p2) + e_1 * (fs_1_2178 * h6_p2 - fs_5_22 * h6_p6 - fs_12_847 * r_2 * h4_p2 + fs_5_63 * r_4 * h2_p2) + fs_315_16 * e_2 * h2_p2;

        pc_1[k] = e_0 * (fs_27_28 * h4_p1 - fs_45_14 * r_2 * h2_p1) + e_1 * (fs_5_2178 * h6_p1 - fs_5_33 * h6_p5 - fs_27_847 * r_2 * h4_p1 + fs_5_126 * r_4 * h2_p1) + fs_315_32 * e_2 * h2_p1;

        pc_2[k] = e_0 * (fs_135_49 * h4_0 - fs_3_7 * h4_p4 - fs_135_49 * r_2 * h2_0) + e_1 * (fs_5_363 * h6_0 - fs_35_363 * h6_p4 - fs_540_5929 * r_2 * h4_0 + fs_12_847 * r_2 * h4_p4 + fs_5_147 * r_4 * h2_0) + fs_135_16 * e_2 * h2_0;
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph2_p1, ph4_m3, ph4_m2, ph4_m1, ph4_p1, ph4_p3, ph6_m3, ph6_m2, ph6_m1, ph6_p1, ph6_p3, ab_2, pc_3, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];

        pc_3[k] = e_0 * (-fs_75_49 * h4_p1 - fs_27_28 * h4_p3 + fs_45_98 * r_2 * h2_p1) + e_1 * (-fs_35_2178 * h6_p1 - fs_7_121 * h6_p3 + fs_300_5929 * r_2 * h4_p1 + fs_27_847 * r_2 * h4_p3 - fs_5_882 * r_4 * h2_p1) - fs_45_32 * e_2 * h2_p1;

        pc_4[k] = e_0 * (fs_135_49 * h4_m2 - f_3_7 * r_2 * h2_m2) + e_1 * (fs_70_1089 * h6_m2 - fs_540_5929 * r_2 * h4_m2 + f_1_21 * r_4 * h2_m2) + f_3_4 * e_2 * h2_m2;

        pc_5[k] = e_0 * (fs_27_28 * h4_m3 - fs_75_49 * h4_m1 + fs_45_98 * r_2 * h2_m1) + e_1 * (fs_7_121 * h6_m3 - fs_35_2178 * h6_m1 - fs_27_847 * r_2 * h4_m3 + fs_300_5929 * r_2 * h4_m1 - fs_5_882 * r_4 * h2_m1) - fs_45_32 * e_2 * h2_m1;
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph4_m4, ph4_m2, ph4_m1, ph4_p3, ph6_m6, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph6_p3, ph6_p5, ab_2, pc_6, pc_7, pc_8, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];

        pc_6[k] = fs_3_7 * e_0 * h4_m4 + e_1 * (fs_35_363 * h6_m4 - fs_12_847 * r_2 * h4_m4);

        pc_7[k] = e_0 * (-fs_27_28 * h4_m1 + fs_45_14 * r_2 * h2_m1) + e_1 * (fs_5_33 * h6_m5 - fs_5_2178 * h6_m1 + fs_27_847 * r_2 * h4_m1 - fs_5_126 * r_4 * h2_m1) - fs_315_32 * e_2 * h2_m1;

        pc_8[k] = e_0 * (-fs_3_7 * h4_m2 + fs_45_7 * r_2 * h2_m2) + e_1 * (fs_5_22 * h6_m6 - fs_1_2178 * h6_m2 + fs_12_847 * r_2 * h4_m2 - fs_5_63 * r_4 * h2_m2) - fs_315_16 * e_2 * h2_m2;

        pc_9[k] = -fs_3_2 * e_0 * h4_p3 + e_1 * (-fs_1_242 * h6_p3 - fs_5_66 * h6_p5 + fs_6_121 * r_2 * h4_p3);
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph6_p4, ab_2, pc_10, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];

        pc_10[k] = e_0 * (-fs_75_56 * h4_p2 + fs_3_2 * h4_p4 - fs_45_14 * r_2 * h2_p2) + e_1 * (-f_4_33 * h6_p2 - fs_40_363 * h6_p4 + fs_75_1694 * r_2 * h4_p2 - fs_6_121 * r_2 * h4_p4 + fs_5_126 * r_4 * h2_p2) + fs_315_32 * e_2 * h2_p2;

        pc_11[k] = e_0 * (-fs_243_392 * h4_p1 + fs_75_56 * h4_p3 - fs_180_49 * r_2 * h2_p1) + e_1 * (-fs_35_1089 * h6_p1 - fs_14_121 * h6_p3 + fs_243_11858 * r_2 * h4_p1 - fs_75_1694 * r_2 * h4_p3 + fs_20_441 * r_4 * h2_p1) + fs_45_4 * e_2 * h2_p1;

        pc_12[k] = e_0 * (-fs_15_98 * h4_0 + fs_243_392 * h4_p2 - fs_270_49 * r_2 * h2_0 - fs_45_98 * r_2 * h2_p2) + e_1 * (-fs_40_363 * h6_0 - fs_112_1089 * h6_p2 + fs_30_5929 * r_2 * h4_0 - fs_243_11858 * r_2 * h4_p2 + fs_10_147 * r_4 * h2_0 + fs_5_882 * r_4 * h2_p2) + e_2 * (fs_135_8 * h2_0 + fs_45_32 * h2_p2);
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph2_0, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph4_0, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ph6_0, ab_2, pc_13, pc_14, pc_15, pc_16, pc_17, pc_18, pc_19, pc_20, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];

        pc_13[k] = e_0 * (-fs_15_98 * h4_m1 + f_12_7 * r_2 * h2_m1) + e_1 * (fs_175_1089 * h6_m1 + fs_30_5929 * r_2 * h4_m1 - f_4_21 * r_4 * h2_m1) - f_3 * e_2 * h2_m1;

        pc_14[k] = e_0 * (-fs_243_392 * h4_m2 + fs_45_98 * r_2 * h2_m2) + e_1 * (fs_112_1089 * h6_m2 + fs_243_11858 * r_2 * h4_m2 - fs_5_882 * r_4 * h2_m2) - fs_45_32 * e_2 * h2_m2;

        pc_15[k] = e_0 * (-fs_75_56 * h4_m3 + fs_243_392 * h4_m1 + fs_180_49 * r_2 * h2_m1) + e_1 * (fs_14_121 * h6_m3 + fs_35_1089 * h6_m1 + fs_75_1694 * r_2 * h4_m3 - fs_243_11858 * r_2 * h4_m1 - fs_20_441 * r_4 * h2_m1) - fs_45_4 * e_2 * h2_m1;

        pc_16[k] = e_0 * (-fs_3_2 * h4_m4 + fs_75_56 * h4_m2 + fs_45_14 * r_2 * h2_m2) + e_1 * (fs_40_363 * h6_m4 + f_4_33 * h6_m2 + fs_6_121 * r_2 * h4_m4 - fs_75_1694 * r_2 * h4_m2 - fs_5_126 * r_4 * h2_m2) - fs_315_32 * e_2 * h2_m2;

        pc_17[k] = fs_3_2 * e_0 * h4_m3 + e_1 * (fs_5_66 * h6_m5 + fs_1_242 * h6_m3 - fs_6_121 * r_2 * h4_m3);

        pc_18[k] = f_2 * e_0 * h4_m4 + e_1 * (fs_5_121 * h6_m4 - f_4_11 * r_2 * h4_m4);

        pc_19[k] = f_1_2 * e_0 * h4_m3 + e_1 * (fs_12_121 * h6_m3 - f_1_11 * r_2 * h4_m3);

        pc_20[k] = e_0 * (-f_4_7 * h4_m2 - fs_135_49 * r_2 * h2_m2) + e_1 * (fs_56_363 * h6_m2 + f_8_77 * r_2 * h4_m2 + fs_5_147 * r_4 * h2_m2) + fs_135_16 * e_2 * h2_m2;

        pc_21[k] = e_0 * (-f_17_14 * h4_m1 - fs_270_49 * r_2 * h2_m1) + e_1 * (fs_70_363 * h6_m1 + f_17_77 * r_2 * h4_m1 + fs_10_147 * r_4 * h2_m1) + fs_135_8 * e_2 * h2_m1;

        pc_22[k] = e_0 * (-f_10_7 * h4_0 - f_18_7 * r_2 * h2_0) + e_1 * (f_5_11 * h6_0 + f_20_77 * r_2 * h4_0 + f_2_7 * r_4 * h2_0) + f_9_2 * e_2 * h2_0;
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_p1, ph2_p2, ph4_m3, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ph6_m5, ph6_m3, ph6_p1, ph6_p2, ph6_p3, ph6_p4, ab_2, pc_23, pc_24, pc_25, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];

        pc_23[k] = e_0 * (-f_17_14 * h4_p1 - fs_270_49 * r_2 * h2_p1) + e_1 * (fs_70_363 * h6_p1 + f_17_77 * r_2 * h4_p1 + fs_10_147 * r_4 * h2_p1) + fs_135_8 * e_2 * h2_p1;

        pc_24[k] = e_0 * (-f_4_7 * h4_p2 - fs_135_49 * r_2 * h2_p2) + e_1 * (fs_56_363 * h6_p2 + f_8_77 * r_2 * h4_p2 + fs_5_147 * r_4 * h2_p2) + fs_135_16 * e_2 * h2_p2;

        pc_25[k] = f_1_2 * e_0 * h4_p3 + e_1 * (fs_12_121 * h6_p3 - f_1_11 * r_2 * h4_p3);

        pc_26[k] = f_2 * e_0 * h4_p4 + e_1 * (fs_5_121 * h6_p4 - f_4_11 * r_2 * h4_p4);

        pc_27[k] = -fs_3_2 * e_0 * h4_m3 + e_1 * (fs_5_66 * h6_m5 - fs_1_242 * h6_m3 + fs_6_121 * r_2 * h4_m3);
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph2_p1, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph4_p1, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ph6_p1, ab_2, pc_28, pc_29, pc_30, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];

        pc_28[k] = e_0 * (-fs_3_2 * h4_m4 - fs_75_56 * h4_m2 - fs_45_14 * r_2 * h2_m2) + e_1 * (fs_40_363 * h6_m4 - f_4_33 * h6_m2 + fs_6_121 * r_2 * h4_m4 + fs_75_1694 * r_2 * h4_m2 + fs_5_126 * r_4 * h2_m2) + fs_315_32 * e_2 * h2_m2;

        pc_29[k] = e_0 * (-fs_75_56 * h4_m3 - fs_243_392 * h4_m1 - fs_180_49 * r_2 * h2_m1) + e_1 * (fs_14_121 * h6_m3 - fs_35_1089 * h6_m1 + fs_75_1694 * r_2 * h4_m3 + fs_243_11858 * r_2 * h4_m1 + fs_20_441 * r_4 * h2_m1) + fs_45_4 * e_2 * h2_m1;

        pc_30[k] = e_0 * (-fs_243_392 * h4_m2 + fs_45_98 * r_2 * h2_m2) + e_1 * (fs_112_1089 * h6_m2 + fs_243_11858 * r_2 * h4_m2 - fs_5_882 * r_4 * h2_m2) - fs_45_32 * e_2 * h2_m2;

        pc_31[k] = e_0 * (-fs_15_98 * h4_p1 + f_12_7 * r_2 * h2_p1) + e_1 * (fs_175_1089 * h6_p1 + fs_30_5929 * r_2 * h4_p1 - f_4_21 * r_4 * h2_p1) - f_3 * e_2 * h2_p1;
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ab_2, pc_32, pc_33, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];

        pc_32[k] = e_0 * (-fs_15_98 * h4_0 - fs_243_392 * h4_p2 - fs_270_49 * r_2 * h2_0 + fs_45_98 * r_2 * h2_p2) + e_1 * (-fs_40_363 * h6_0 + fs_112_1089 * h6_p2 + fs_30_5929 * r_2 * h4_0 + fs_243_11858 * r_2 * h4_p2 + fs_10_147 * r_4 * h2_0 - fs_5_882 * r_4 * h2_p2) + e_2 * (fs_135_8 * h2_0 - fs_45_32 * h2_p2);

        pc_33[k] = e_0 * (-fs_243_392 * h4_p1 - fs_75_56 * h4_p3 - fs_180_49 * r_2 * h2_p1) + e_1 * (-fs_35_1089 * h6_p1 + fs_14_121 * h6_p3 + fs_243_11858 * r_2 * h4_p1 + fs_75_1694 * r_2 * h4_p3 + fs_20_441 * r_4 * h2_p1) + fs_45_4 * e_2 * h2_p1;

        pc_34[k] = e_0 * (-fs_75_56 * h4_p2 - fs_3_2 * h4_p4 - fs_45_14 * r_2 * h2_p2) + e_1 * (-f_4_33 * h6_p2 + fs_40_363 * h6_p4 + fs_75_1694 * r_2 * h4_p2 + fs_6_121 * r_2 * h4_p4 + fs_5_126 * r_4 * h2_p2) + fs_315_32 * e_2 * h2_p2;

        pc_35[k] = -fs_3_2 * e_0 * h4_p3 + e_1 * (-fs_1_242 * h6_p3 + fs_5_66 * h6_p5 + fs_6_121 * r_2 * h4_p3);
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph6_m6, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ab_2, pc_36, pc_37, pc_38, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];

        pc_36[k] = e_0 * (fs_3_7 * h4_m2 - fs_45_7 * r_2 * h2_m2) + e_1 * (fs_5_22 * h6_m6 + fs_1_2178 * h6_m2 - fs_12_847 * r_2 * h4_m2 + fs_5_63 * r_4 * h2_m2) + fs_315_16 * e_2 * h2_m2;

        pc_37[k] = e_0 * (fs_27_28 * h4_m1 - fs_45_14 * r_2 * h2_m1) + e_1 * (fs_5_33 * h6_m5 + fs_5_2178 * h6_m1 - fs_27_847 * r_2 * h4_m1 + fs_5_126 * r_4 * h2_m1) + fs_315_32 * e_2 * h2_m1;

        pc_38[k] = fs_3_7 * e_0 * h4_m4 + e_1 * (fs_35_363 * h6_m4 - fs_12_847 * r_2 * h4_m4);

        pc_39[k] = e_0 * (fs_27_28 * h4_m3 + fs_75_49 * h4_m1 - fs_45_98 * r_2 * h2_m1) + e_1 * (fs_7_121 * h6_m3 + fs_35_2178 * h6_m1 - fs_27_847 * r_2 * h4_m3 - fs_300_5929 * r_2 * h4_m1 + fs_5_882 * r_4 * h2_m1) + fs_45_32 * e_2 * h2_m1;
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ab_2, pc_40, pc_41, pc_42, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];

        pc_40[k] = e_0 * (fs_135_49 * h4_p2 - f_3_7 * r_2 * h2_p2) + e_1 * (fs_70_1089 * h6_p2 - fs_540_5929 * r_2 * h4_p2 + f_1_21 * r_4 * h2_p2) + f_3_4 * e_2 * h2_p2;

        pc_41[k] = e_0 * (-fs_75_49 * h4_p1 + fs_27_28 * h4_p3 + fs_45_98 * r_2 * h2_p1) + e_1 * (-fs_35_2178 * h6_p1 + fs_7_121 * h6_p3 + fs_300_5929 * r_2 * h4_p1 - fs_27_847 * r_2 * h4_p3 - fs_5_882 * r_4 * h2_p1) - fs_45_32 * e_2 * h2_p1;

        pc_42[k] = e_0 * (fs_135_49 * h4_0 + fs_3_7 * h4_p4 - fs_135_49 * r_2 * h2_0) + e_1 * (fs_5_363 * h6_0 + fs_35_363 * h6_p4 - fs_540_5929 * r_2 * h4_0 - fs_12_847 * r_2 * h4_p4 + fs_5_147 * r_4 * h2_0) + fs_135_16 * e_2 * h2_0;

        pc_43[k] = e_0 * (fs_27_28 * h4_p1 - fs_45_14 * r_2 * h2_p1) + e_1 * (fs_5_2178 * h6_p1 + fs_5_33 * h6_p5 - fs_27_847 * r_2 * h4_p1 + fs_5_126 * r_4 * h2_p1) + fs_315_32 * e_2 * h2_p1;
    }

    // NOTE: the rows are formed in 11 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];

        pc_44[k] = e_0 * (fs_3_7 * h4_p2 - fs_45_7 * r_2 * h2_p2) + e_1 * (fs_1_2178 * h6_p2 + fs_5_22 * h6_p6 - fs_12_847 * r_2 * h4_p2 + fs_5_63 * r_4 * h2_p2) + fs_315_16 * e_2 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest, and
    // the atom pairs beyond the reach of every pair of primitives are set to zero.

    const size_t sources[45] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44};

    for (size_t m = 0; m < 45; m++)
    {
        const auto *pc = buffer.data(3 + sources[m]);

        std::copy(pc, pc + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdovl
