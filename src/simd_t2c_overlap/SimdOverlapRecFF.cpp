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



#include "SimdOverlapRecFF.hpp"

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
compute_ff_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 3) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFF.compute_ff_overlap: Basis functions must be of angular momenta three and three"));
    }

    if (harmonics.size() < 6)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFF.compute_ff_overlap: Harmonics must reach angular momentum six"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFF.compute_ff_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 49 * nvalues, 0.0);

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

            const auto f_0 = fbase * fmu / fexp / fexp / fexp;

            const auto f_1 = fbase * fmu * fmu / fexp / fexp / fexp;

            const auto f_2 = fbase * fmu * fmu * fmu / fexp / fexp / fexp;

            const auto f_3 = fbase / fexp / fexp / fexp;

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
    auto *pc_7 = values + 8 * nvalues;
    auto *pc_8 = values + 9 * nvalues;
    auto *pc_9 = values + 10 * nvalues;
    auto *pc_10 = values + 11 * nvalues;
    auto *pc_11 = values + 12 * nvalues;
    auto *pc_12 = values + 13 * nvalues;
    auto *pc_13 = values + 16 * nvalues;
    auto *pc_14 = values + 17 * nvalues;
    auto *pc_15 = values + 18 * nvalues;
    auto *pc_16 = values + 19 * nvalues;
    auto *pc_17 = values + 20 * nvalues;
    auto *pc_18 = values + 24 * nvalues;
    auto *pc_19 = values + 25 * nvalues;
    auto *pc_20 = values + 26 * nvalues;
    auto *pc_21 = values + 27 * nvalues;
    auto *pc_22 = values + 32 * nvalues;
    auto *pc_23 = values + 33 * nvalues;
    auto *pc_24 = values + 34 * nvalues;
    auto *pc_25 = values + 40 * nvalues;
    auto *pc_26 = values + 41 * nvalues;
    auto *pc_27 = values + 48 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_15_4 = 3.75;
    const auto f_9_14 = 9.0 / 14.0;
    const auto f_15_7 = 15.0 / 7.0;
    const auto f_3_2 = 1.5;
    const auto f_5_231 = 5.0 / 231.0;
    const auto fs_50_231 = std::sqrt(50.0 / 231.0);
    const auto f_9_77 = 9.0 / 77.0;
    const auto f_5_21 = 5.0 / 21.0;
    const auto f_1_7 = 1.0 / 7.0;
    const auto f_15_8 = 1.875;
    const auto fs_225_32 = std::sqrt(7.03125);
    const auto fs_135_196 = std::sqrt(135.0 / 196.0);
    const auto fs_225_98 = std::sqrt(225.0 / 98.0);
    const auto fs_25_15246 = std::sqrt(25.0 / 15246.0);
    const auto fs_25_231 = std::sqrt(25.0 / 231.0);
    const auto fs_135_5929 = std::sqrt(135.0 / 5929.0);
    const auto fs_25_882 = std::sqrt(25.0 / 882.0);
    const auto fs_45_16 = std::sqrt(2.8125);
    const auto fs_243_196 = std::sqrt(243.0 / 196.0);
    const auto fs_27_28 = std::sqrt(27.0 / 28.0);
    const auto fs_45_49 = std::sqrt(45.0 / 49.0);
    const auto fs_50_7623 = std::sqrt(50.0 / 7623.0);
    const auto fs_125_2541 = std::sqrt(125.0 / 2541.0);
    const auto fs_243_5929 = std::sqrt(243.0 / 5929.0);
    const auto fs_27_847 = std::sqrt(27.0 / 847.0);
    const auto fs_5_441 = std::sqrt(5.0 / 441.0);
    const auto fs_81_28 = std::sqrt(81.0 / 28.0);
    const auto fs_100_2541 = std::sqrt(100.0 / 2541.0);
    const auto fs_81_847 = std::sqrt(81.0 / 847.0);
    const auto fs_45_28 = std::sqrt(45.0 / 28.0);
    const auto f_10_77 = 10.0 / 77.0;
    const auto fs_100_847 = std::sqrt(100.0 / 847.0);
    const auto f_3_11 = 3.0 / 11.0;
    const auto fs_45_847 = std::sqrt(45.0 / 847.0);
    const auto fs_135_32 = std::sqrt(4.21875);
    const auto f_6_7 = 6.0 / 7.0;
    const auto fs_9_28 = std::sqrt(9.0 / 28.0);
    const auto fs_135_98 = std::sqrt(135.0 / 98.0);
    const auto fs_125_5082 = std::sqrt(125.0 / 5082.0);
    const auto fs_75_847 = std::sqrt(75.0 / 847.0);
    const auto f_12_77 = 12.0 / 77.0;
    const auto fs_9_847 = std::sqrt(9.0 / 847.0);
    const auto fs_5_294 = std::sqrt(5.0 / 294.0);
    const auto fs_45_4 = std::sqrt(11.25);
    const auto fs_27_196 = std::sqrt(27.0 / 196.0);
    const auto fs_180_49 = std::sqrt(180.0 / 49.0);
    const auto fs_800_7623 = std::sqrt(800.0 / 7623.0);
    const auto fs_27_5929 = std::sqrt(27.0 / 5929.0);
    const auto fs_20_441 = std::sqrt(20.0 / 441.0);
    const auto f_9_4 = 2.25;
    const auto fs_27_4 = std::sqrt(6.75);
    const auto f_3_14 = 3.0 / 14.0;
    const auto f_9_7 = 9.0 / 7.0;
    const auto fs_108_49 = std::sqrt(108.0 / 49.0);
    const auto f_25_77 = 25.0 / 77.0;
    const auto fs_250_2541 = std::sqrt(250.0 / 2541.0);
    const auto f_3_77 = 3.0 / 77.0;
    const auto fs_180_5929 = std::sqrt(180.0 / 5929.0);
    const auto fs_4_147 = std::sqrt(4.0 / 147.0);
    const auto fs_9_8 = std::sqrt(1.125);
    const auto fs_18_49 = std::sqrt(18.0 / 49.0);
    const auto fs_1250_7623 = std::sqrt(1250.0 / 7623.0);
    const auto fs_2_441 = std::sqrt(2.0 / 441.0);
    const auto f_3 = 3.0;
    const auto f_12_7 = 12.0 / 7.0;
    const auto f_100_231 = 100.0 / 231.0;
    const auto f_18_77 = 18.0 / 77.0;
    const auto f_4_21 = 4.0 / 21.0;

    // NOTE: the rows are formed in 8 loops, as the vectorizer runs out of
    // registers with all 28 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p4, ph6_0, ph6_p1, ph6_p2, ph6_p4, ph6_p5, ph6_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_0[k] = e_0 * (f_15_4 * h2_0 - f_15_4 * r_2) + e_1 * (f_9_14 * h4_0 - f_15_7 * r_2 * h2_0 + f_3_2 * r_4) + e_2 * (f_5_231 * h6_0 + fs_50_231 * h6_p6 - f_9_77 * r_2 * h4_0 + f_5_21 * r_4 * h2_0 - f_1_7 * r_6) + f_15_8 * e_3;

        pc_1[k] = -fs_225_32 * e_0 * h2_p1 + e_1 * (-fs_135_196 * h4_p1 + fs_225_98 * r_2 * h2_p1) + e_2 * (-fs_25_15246 * h6_p1 + fs_25_231 * h6_p5 + fs_135_5929 * r_2 * h4_p1 - fs_25_882 * r_4 * h2_p1);

        pc_2[k] = fs_45_16 * e_0 * h2_p2 + e_1 * (fs_243_196 * h4_p2 + fs_27_28 * h4_p4 - fs_45_49 * r_2 * h2_p2) + e_2 * (fs_50_7623 * h6_p2 + fs_125_2541 * h6_p4 - fs_243_5929 * r_2 * h4_p2 - fs_27_847 * r_2 * h4_p4 + fs_5_441 * r_4 * h2_p2);
    }

    // NOTE: the rows are formed in 8 loops, as the vectorizer runs out of
    // registers with all 28 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph6_m6, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ab_2 : simd::cache_line_size())
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

        pc_3[k] = -fs_81_28 * e_1 * h4_m3 + e_2 * (-fs_100_2541 * h6_m3 + fs_81_847 * r_2 * h4_m3);

        pc_4[k] = fs_45_16 * e_0 * h2_m2 + e_1 * (-fs_27_28 * h4_m4 + fs_243_196 * h4_m2 - fs_45_49 * r_2 * h2_m2) + e_2 * (-fs_125_2541 * h6_m4 + fs_50_7623 * h6_m2 + fs_27_847 * r_2 * h4_m4 - fs_243_5929 * r_2 * h4_m2 + fs_5_441 * r_4 * h2_m2);

        pc_5[k] = -fs_225_32 * e_0 * h2_m1 + e_1 * (-fs_135_196 * h4_m1 + fs_225_98 * r_2 * h2_m1) + e_2 * (-fs_25_231 * h6_m5 - fs_25_15246 * h6_m1 + fs_135_5929 * r_2 * h4_m1 - fs_25_882 * r_4 * h2_m1);

        pc_6[k] = -fs_50_231 * e_2 * h6_m6;
    }

    // NOTE: the rows are formed in 8 loops, as the vectorizer runs out of
    // registers with all 28 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_p1, ph4_m2, ph4_0, ph4_p1, ph4_p3, ph4_p4, ph6_m2, ph6_0, ph6_p1, ph6_p3, ph6_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];

        pc_7[k] = -f_15_4 * e_0 * r_2 + e_1 * (-f_3_2 * h4_0 - fs_45_28 * h4_p4 + f_3_2 * r_4) + e_2 * (-f_10_77 * h6_0 + fs_100_847 * h6_p4 + f_3_11 * r_2 * h4_0 + fs_45_847 * r_2 * h4_p4 - f_1_7 * r_6) + f_15_8 * e_3;

        pc_8[k] = -fs_135_32 * e_0 * h2_p1 + e_1 * (f_6_7 * h4_p1 - fs_9_28 * h4_p3 + fs_135_98 * r_2 * h2_p1) + e_2 * (fs_125_5082 * h6_p1 + fs_75_847 * h6_p3 - f_12_77 * r_2 * h4_p1 + fs_9_847 * r_2 * h4_p3 - fs_5_294 * r_4 * h2_p1);

        pc_9[k] = fs_45_4 * e_0 * h2_m2 + e_1 * (-fs_27_196 * h4_m2 - fs_180_49 * r_2 * h2_m2) + e_2 * (-fs_800_7623 * h6_m2 + fs_27_5929 * r_2 * h4_m2 + fs_20_441 * r_4 * h2_m2);
    }

    // NOTE: the rows are formed in 8 loops, as the vectorizer runs out of
    // registers with all 28 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];

        pc_10[k] = -fs_135_32 * e_0 * h2_m1 + e_1 * (fs_9_28 * h4_m3 + f_6_7 * h4_m1 + fs_135_98 * r_2 * h2_m1) + e_2 * (-fs_75_847 * h6_m3 + fs_125_5082 * h6_m1 - fs_9_847 * r_2 * h4_m3 - f_12_77 * r_2 * h4_m1 - fs_5_294 * r_4 * h2_m1);

        pc_11[k] = fs_45_28 * e_1 * h4_m4 + e_2 * (-fs_100_847 * h6_m4 - fs_45_847 * r_2 * h4_m4);

        pc_12[k] = fs_225_32 * e_0 * h2_m1 + e_1 * (fs_135_196 * h4_m1 - fs_225_98 * r_2 * h2_m1) + e_2 * (-fs_25_231 * h6_m5 + fs_25_15246 * h6_m1 - fs_135_5929 * r_2 * h4_m1 + fs_25_882 * r_4 * h2_m1);
    }

    // NOTE: the rows are formed in 8 loops, as the vectorizer runs out of
    // registers with all 28 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph2_0, ph2_p2, ph4_m2, ph4_m1, ph4_0, ph4_p2, ph6_m2, ph6_m1, ph6_0, ph6_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];

        pc_13[k] = e_0 * (-f_9_4 * h2_0 + fs_27_4 * h2_p2 - f_15_4 * r_2) + e_1 * (f_3_14 * h4_0 - fs_45_49 * h4_p2 + f_9_7 * r_2 * h2_0 - fs_108_49 * r_2 * h2_p2 + f_3_2 * r_4) + e_2 * (f_25_77 * h6_0 + fs_250_2541 * h6_p2 - f_3_77 * r_2 * h4_0 + fs_180_5929 * r_2 * h4_p2 - f_1_7 * r_4 * h2_0 + fs_4_147 * r_4 * h2_p2 - f_1_7 * r_6) + f_15_8 * e_3;

        pc_14[k] = -fs_9_8 * e_0 * h2_m1 + e_1 * (fs_135_196 * h4_m1 + fs_18_49 * r_2 * h2_m1) + e_2 * (-fs_1250_7623 * h6_m1 - fs_135_5929 * r_2 * h4_m1 - fs_2_441 * r_4 * h2_m1);

        pc_15[k] = -fs_27_4 * e_0 * h2_m2 + e_1 * (fs_45_49 * h4_m2 + fs_108_49 * r_2 * h2_m2) + e_2 * (-fs_250_2541 * h6_m2 - fs_180_5929 * r_2 * h4_m2 - fs_4_147 * r_4 * h2_m2);
    }

    // NOTE: the rows are formed in 8 loops, as the vectorizer runs out of
    // registers with all 28 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph2_0, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph4_0, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ph6_0, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];

        pc_16[k] = fs_135_32 * e_0 * h2_m1 + e_1 * (fs_9_28 * h4_m3 - f_6_7 * h4_m1 - fs_135_98 * r_2 * h2_m1) + e_2 * (-fs_75_847 * h6_m3 - fs_125_5082 * h6_m1 - fs_9_847 * r_2 * h4_m3 + f_12_77 * r_2 * h4_m1 + fs_5_294 * r_4 * h2_m1);

        pc_17[k] = -fs_45_16 * e_0 * h2_m2 + e_1 * (-fs_27_28 * h4_m4 - fs_243_196 * h4_m2 + fs_45_49 * r_2 * h2_m2) + e_2 * (-fs_125_2541 * h6_m4 - fs_50_7623 * h6_m2 + fs_27_847 * r_2 * h4_m4 + fs_243_5929 * r_2 * h4_m2 - fs_5_441 * r_4 * h2_m2);

        pc_18[k] = e_0 * (-f_3 * h2_0 - f_15_4 * r_2) + e_1 * (f_9_7 * h4_0 + f_12_7 * r_2 * h2_0 + f_3_2 * r_4) + e_2 * (-f_100_231 * h6_0 - f_18_77 * r_2 * h4_0 - f_4_21 * r_4 * h2_0 - f_1_7 * r_6) + f_15_8 * e_3;
    }

    // NOTE: the rows are formed in 8 loops, as the vectorizer runs out of
    // registers with all 28 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph6_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

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

        pc_19[k] = -fs_9_8 * e_0 * h2_p1 + e_1 * (fs_135_196 * h4_p1 + fs_18_49 * r_2 * h2_p1) + e_2 * (-fs_1250_7623 * h6_p1 - fs_135_5929 * r_2 * h4_p1 - fs_2_441 * r_4 * h2_p1);

        pc_20[k] = fs_45_4 * e_0 * h2_p2 + e_1 * (-fs_27_196 * h4_p2 - fs_180_49 * r_2 * h2_p2) + e_2 * (-fs_800_7623 * h6_p2 + fs_27_5929 * r_2 * h4_p2 + fs_20_441 * r_4 * h2_p2);

        pc_21[k] = -fs_81_28 * e_1 * h4_p3 + e_2 * (-fs_100_2541 * h6_p3 + fs_81_847 * r_2 * h4_p3);

        pc_22[k] = e_0 * (-f_9_4 * h2_0 - fs_27_4 * h2_p2 - f_15_4 * r_2) + e_1 * (f_3_14 * h4_0 + fs_45_49 * h4_p2 + f_9_7 * r_2 * h2_0 + fs_108_49 * r_2 * h2_p2 + f_3_2 * r_4) + e_2 * (f_25_77 * h6_0 - fs_250_2541 * h6_p2 - f_3_77 * r_2 * h4_0 - fs_180_5929 * r_2 * h4_p2 - f_1_7 * r_4 * h2_0 - fs_4_147 * r_4 * h2_p2 - f_1_7 * r_6) + f_15_8 * e_3;

        pc_23[k] = -fs_135_32 * e_0 * h2_p1 + e_1 * (f_6_7 * h4_p1 + fs_9_28 * h4_p3 + fs_135_98 * r_2 * h2_p1) + e_2 * (fs_125_5082 * h6_p1 - fs_75_847 * h6_p3 - f_12_77 * r_2 * h4_p1 - fs_9_847 * r_2 * h4_p3 - fs_5_294 * r_4 * h2_p1);

        pc_24[k] = fs_45_16 * e_0 * h2_p2 + e_1 * (fs_243_196 * h4_p2 - fs_27_28 * h4_p4 - fs_45_49 * r_2 * h2_p2) + e_2 * (fs_50_7623 * h6_p2 - fs_125_2541 * h6_p4 - fs_243_5929 * r_2 * h4_p2 + fs_27_847 * r_2 * h4_p4 + fs_5_441 * r_4 * h2_p2);

        pc_25[k] = -f_15_4 * e_0 * r_2 + e_1 * (-f_3_2 * h4_0 + fs_45_28 * h4_p4 + f_3_2 * r_4) + e_2 * (-f_10_77 * h6_0 - fs_100_847 * h6_p4 + f_3_11 * r_2 * h4_0 - fs_45_847 * r_2 * h4_p4 - f_1_7 * r_6) + f_15_8 * e_3;
    }

    // NOTE: the rows are formed in 8 loops, as the vectorizer runs out of
    // registers with all 28 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];

        pc_26[k] = -fs_225_32 * e_0 * h2_p1 + e_1 * (-fs_135_196 * h4_p1 + fs_225_98 * r_2 * h2_p1) + e_2 * (-fs_25_15246 * h6_p1 - fs_25_231 * h6_p5 + fs_135_5929 * r_2 * h4_p1 - fs_25_882 * r_4 * h2_p1);

        pc_27[k] = e_0 * (f_15_4 * h2_0 - f_15_4 * r_2) + e_1 * (f_9_14 * h4_0 - f_15_7 * r_2 * h2_0 + f_3_2 * r_4) + e_2 * (f_5_231 * h6_0 - fs_50_231 * h6_p6 - f_9_77 * r_2 * h4_0 + f_5_21 * r_4 * h2_0 - f_1_7 * r_6) + f_15_8 * e_3;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[49] = {0, 1, 2, 3, 4, 5, 6, 1, 8, 9, 10, 11, 12, 13, 2, 9, 16, 17, 18, 19, 20, 3, 10, 17, 24, 25, 26, 27, 4, 11, 18, 25, 32, 33, 34, 5, 12, 19, 26, 33, 40, 41, 6, 13, 20, 27, 34, 41, 48};

    for (size_t m = 0; m < 49; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
