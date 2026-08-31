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



#include "SimdOverlapRecFD.hpp"

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
compute_fd_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 3) || (ket.get_angular_momentum() != 2))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFD.compute_fd_overlap: Basis functions must be of angular momenta three and two"));
    }

    if (harmonics.size() < 5)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFD.compute_fd_overlap: Harmonics must reach angular momentum five"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecFD.compute_fd_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 35 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

    auto buffer = CSimdMatrix(3, nmax);

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

            const auto f_0 = fbase * bexp * fmu / fexp / fexp / fexp;

            const auto f_1 = fbase * bexp * fmu * fmu / fexp / fexp / fexp;

            const auto f_2 = fbase * bexp / fexp / fexp / fexp;

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_9_20 = std::sqrt(9.0 / 20.0);
    const auto fs_27_10 = std::sqrt(27.0 / 10.0);
    const auto fs_1_882 = std::sqrt(1.0 / 882.0);
    const auto fs_5_21 = std::sqrt(5.0 / 21.0);
    const auto fs_1_45 = std::sqrt(1.0 / 45.0);
    const auto fs_27_490 = std::sqrt(27.0 / 490.0);
    const auto fs_135_32 = std::sqrt(4.21875);
    const auto fs_9_8 = std::sqrt(1.125);
    const auto fs_1_126 = std::sqrt(1.0 / 126.0);
    const auto fs_2_21 = std::sqrt(2.0 / 21.0);
    const auto fs_1_18 = std::sqrt(1.0 / 18.0);
    const auto f_3_2 = 1.5;
    const auto fs_4_63 = std::sqrt(4.0 / 63.0);
    const auto f_1_3 = 1.0 / 3.0;
    const auto fs_9_5 = std::sqrt(9.0 / 5.0);
    const auto fs_5_441 = std::sqrt(5.0 / 441.0);
    const auto fs_1_7 = std::sqrt(1.0 / 7.0);
    const auto fs_4_45 = std::sqrt(4.0 / 45.0);
    const auto fs_9_245 = std::sqrt(9.0 / 245.0);
    const auto fs_45_16 = std::sqrt(2.8125);
    const auto fs_27_40 = std::sqrt(27.0 / 40.0);
    const auto fs_4_147 = std::sqrt(4.0 / 147.0);
    const auto fs_8_63 = std::sqrt(8.0 / 63.0);
    const auto fs_1_30 = std::sqrt(1.0 / 30.0);
    const auto fs_27_25 = std::sqrt(27.0 / 25.0);
    const auto fs_9_50 = std::sqrt(9.0 / 50.0);
    const auto fs_5_294 = std::sqrt(5.0 / 294.0);
    const auto fs_5_63 = std::sqrt(5.0 / 63.0);
    const auto fs_4_75 = std::sqrt(4.0 / 75.0);
    const auto fs_9_2450 = std::sqrt(9.0 / 2450.0);
    const auto fs_9_32 = std::sqrt(0.28125);
    const auto fs_72_25 = std::sqrt(72.0 / 25.0);
    const auto fs_50_441 = std::sqrt(50.0 / 441.0);
    const auto fs_5_42 = std::sqrt(5.0 / 42.0);
    const auto fs_2_225 = std::sqrt(2.0 / 225.0);
    const auto fs_72_1225 = std::sqrt(72.0 / 1225.0);
    const auto fs_9_2 = std::sqrt(4.5);
    const auto f_9_10 = 9.0 / 10.0;
    const auto fs_54_25 = std::sqrt(54.0 / 25.0);
    const auto fs_10_49 = std::sqrt(10.0 / 49.0);
    const auto f_1_5 = 1.0 / 5.0;
    const auto fs_54_1225 = std::sqrt(54.0 / 1225.0);
    const auto fs_27_8 = std::sqrt(3.375);
    const auto fs_80_441 = std::sqrt(80.0 / 441.0);
    const auto fs_27_1225 = std::sqrt(27.0 / 1225.0);
    const auto fs_27_16 = std::sqrt(1.6875);
    const auto f_6_5 = 6.0 / 5.0;
    const auto f_9_5 = 9.0 / 5.0;
    const auto f_10_21 = 10.0 / 21.0;
    const auto f_4_15 = 4.0 / 15.0;
    const auto f_9_35 = 9.0 / 35.0;
    const auto f_9_4 = 2.25;

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_p1, ph3_m3, ph3_m2, ph3_p1, ph3_p2, ph5_m4, ph5_m3, ph5_m2, ph5_p1, ph5_p2, ph5_p4, ph5_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];

        pc_0[k] = e_0 * (-fs_9_20 * h3_p1 + fs_27_10 * r_2 * h1_p1) + e_1 * (-fs_1_882 * h5_p1 + fs_5_21 * h5_p5 + fs_1_45 * r_2 * h3_p1 - fs_27_490 * r_4 * h1_p1) - fs_135_32 * e_2 * h1_p1;

        pc_1[k] = fs_9_8 * e_0 * h3_p2 + e_1 * (fs_1_126 * h5_p2 + fs_2_21 * h5_p4 - fs_1_18 * r_2 * h3_p2);

        pc_2[k] = -f_3_2 * e_0 * h3_m3 + e_1 * (-fs_4_63 * h5_m3 + f_1_3 * r_2 * h3_m3);

        pc_3[k] = fs_9_8 * e_0 * h3_m2 + e_1 * (-fs_2_21 * h5_m4 + fs_1_126 * h5_m2 - fs_1_18 * r_2 * h3_m2);
    }

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_m1, ph1_0, ph1_p1, ph3_m1, ph3_0, ph3_p1, ph3_p3, ph5_m5, ph5_m2, ph5_m1, ph5_0, ph5_p1, ph5_p3, ph5_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];

        pc_4[k] = e_0 * (-fs_9_20 * h3_m1 + fs_27_10 * r_2 * h1_m1) + e_1 * (-fs_5_21 * h5_m5 - fs_1_882 * h5_m1 + fs_1_45 * r_2 * h3_m1 - fs_27_490 * r_4 * h1_m1) - fs_135_32 * e_2 * h1_m1;

        pc_5[k] = e_0 * (-fs_9_5 * h3_0 + fs_9_5 * r_2 * h1_0) + e_1 * (-fs_5_441 * h5_0 + fs_1_7 * h5_p4 + fs_4_45 * r_2 * h3_0 - fs_9_245 * r_4 * h1_0) - fs_45_16 * e_2 * h1_0;

        pc_6[k] = e_0 * (fs_27_40 * h3_p1 - fs_9_8 * h3_p3 + fs_9_5 * r_2 * h1_p1) + e_1 * (fs_4_147 * h5_p1 + fs_8_63 * h5_p3 - fs_1_30 * r_2 * h3_p1 + fs_1_18 * r_2 * h3_p3 - fs_9_245 * r_4 * h1_p1) - fs_45_16 * e_2 * h1_p1;

        pc_7[k] = -fs_1_7 * e_1 * h5_m2;
    }

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_m1, ph1_p1, ph3_m3, ph3_m1, ph3_p1, ph3_p3, ph5_m4, ph5_m3, ph5_m1, ph5_p1, ph5_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];

        pc_8[k] = e_0 * (fs_9_8 * h3_m3 + fs_27_40 * h3_m1 + fs_9_5 * r_2 * h1_m1) + e_1 * (-fs_8_63 * h5_m3 + fs_4_147 * h5_m1 - fs_1_18 * r_2 * h3_m3 - fs_1_30 * r_2 * h3_m1 - fs_9_245 * r_4 * h1_m1) - fs_45_16 * e_2 * h1_m1;

        pc_9[k] = -fs_1_7 * e_1 * h5_m4;

        pc_10[k] = e_0 * (fs_27_25 * h3_p1 + fs_9_20 * h3_p3 - fs_9_50 * r_2 * h1_p1) + e_1 * (fs_5_294 * h5_p1 + fs_5_63 * h5_p3 - fs_4_75 * r_2 * h3_p1 - fs_1_45 * r_2 * h3_p3 + fs_9_2450 * r_4 * h1_p1) + fs_9_32 * e_2 * h1_p1;
    }

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_m1, ph1_0, ph3_m3, ph3_m2, ph3_m1, ph3_0, ph3_p2, ph5_m3, ph5_m2, ph5_m1, ph5_0, ph5_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];

        pc_11[k] = e_0 * (fs_9_50 * h3_0 - fs_27_40 * h3_p2 + fs_72_25 * r_2 * h1_0) + e_1 * (fs_50_441 * h5_0 + fs_5_42 * h5_p2 - fs_2_225 * r_2 * h3_0 + fs_1_30 * r_2 * h3_p2 - fs_72_1225 * r_4 * h1_0) - fs_9_2 * e_2 * h1_0;

        pc_12[k] = e_0 * (f_9_10 * h3_m1 + fs_54_25 * r_2 * h1_m1) + e_1 * (-fs_10_49 * h5_m1 - f_1_5 * r_2 * h3_m1 - fs_54_1225 * r_4 * h1_m1) - fs_27_8 * e_2 * h1_m1;

        pc_13[k] = fs_27_40 * e_0 * h3_m2 + e_1 * (-fs_5_42 * h5_m2 - fs_1_30 * r_2 * h3_m2);

        pc_14[k] = e_0 * (-fs_9_20 * h3_m3 - fs_27_25 * h3_m1 + fs_9_50 * r_2 * h1_m1) + e_1 * (-fs_5_63 * h5_m3 - fs_5_294 * h5_m1 + fs_1_45 * r_2 * h3_m3 + fs_4_75 * r_2 * h3_m1 - fs_9_2450 * r_4 * h1_m1) - fs_9_32 * e_2 * h1_m1;

        pc_15[k] = -fs_9_5 * e_0 * h3_m2 + e_1 * (-fs_5_63 * h5_m2 + fs_4_45 * r_2 * h3_m2);

        pc_16[k] = e_0 * (fs_9_50 * h3_m1 - fs_27_25 * r_2 * h1_m1) + e_1 * (-fs_80_441 * h5_m1 - fs_2_225 * r_2 * h3_m1 + fs_27_1225 * r_4 * h1_m1) + fs_27_16 * e_2 * h1_m1;

        pc_17[k] = e_0 * (f_6_5 * h3_0 + f_9_5 * r_2 * h1_0) + e_1 * (-f_10_21 * h5_0 - f_4_15 * r_2 * h3_0 - f_9_35 * r_4 * h1_0) - f_9_4 * e_2 * h1_0;
    }

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_m1, ph1_p1, ph3_m3, ph3_m2, ph3_m1, ph3_p1, ph3_p2, ph5_m3, ph5_m2, ph5_m1, ph5_p1, ph5_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];

        pc_18[k] = e_0 * (fs_9_50 * h3_p1 - fs_27_25 * r_2 * h1_p1) + e_1 * (-fs_80_441 * h5_p1 - fs_2_225 * r_2 * h3_p1 + fs_27_1225 * r_4 * h1_p1) + fs_27_16 * e_2 * h1_p1;

        pc_19[k] = -fs_9_5 * e_0 * h3_p2 + e_1 * (-fs_5_63 * h5_p2 + fs_4_45 * r_2 * h3_p2);

        pc_20[k] = e_0 * (-fs_9_20 * h3_m3 + fs_27_25 * h3_m1 - fs_9_50 * r_2 * h1_m1) + e_1 * (-fs_5_63 * h5_m3 + fs_5_294 * h5_m1 + fs_1_45 * r_2 * h3_m3 - fs_4_75 * r_2 * h3_m1 + fs_9_2450 * r_4 * h1_m1) + fs_9_32 * e_2 * h1_m1;

        pc_21[k] = fs_27_40 * e_0 * h3_m2 + e_1 * (-fs_5_42 * h5_m2 - fs_1_30 * r_2 * h3_m2);

        pc_22[k] = e_0 * (f_9_10 * h3_p1 + fs_54_25 * r_2 * h1_p1) + e_1 * (-fs_10_49 * h5_p1 - f_1_5 * r_2 * h3_p1 - fs_54_1225 * r_4 * h1_p1) - fs_27_8 * e_2 * h1_p1;
    }

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_0, ph1_p1, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_m4, ph5_0, ph5_p1, ph5_p2, ph5_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];

        pc_23[k] = e_0 * (fs_9_50 * h3_0 + fs_27_40 * h3_p2 + fs_72_25 * r_2 * h1_0) + e_1 * (fs_50_441 * h5_0 - fs_5_42 * h5_p2 - fs_2_225 * r_2 * h3_0 - fs_1_30 * r_2 * h3_p2 - fs_72_1225 * r_4 * h1_0) - fs_9_2 * e_2 * h1_0;

        pc_24[k] = e_0 * (fs_27_25 * h3_p1 - fs_9_20 * h3_p3 - fs_9_50 * r_2 * h1_p1) + e_1 * (fs_5_294 * h5_p1 - fs_5_63 * h5_p3 - fs_4_75 * r_2 * h3_p1 + fs_1_45 * r_2 * h3_p3 + fs_9_2450 * r_4 * h1_p1) + fs_9_32 * e_2 * h1_p1;

        pc_25[k] = -fs_1_7 * e_1 * h5_m4;
    }

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_m1, ph1_p1, ph3_m3, ph3_m1, ph3_p1, ph3_p3, ph5_m3, ph5_m1, ph5_p1, ph5_p2, ph5_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];

        pc_26[k] = e_0 * (fs_9_8 * h3_m3 - fs_27_40 * h3_m1 - fs_9_5 * r_2 * h1_m1) + e_1 * (-fs_8_63 * h5_m3 - fs_4_147 * h5_m1 - fs_1_18 * r_2 * h3_m3 + fs_1_30 * r_2 * h3_m1 + fs_9_245 * r_4 * h1_m1) + fs_45_16 * e_2 * h1_m1;

        pc_27[k] = -fs_1_7 * e_1 * h5_p2;

        pc_28[k] = e_0 * (fs_27_40 * h3_p1 + fs_9_8 * h3_p3 + fs_9_5 * r_2 * h1_p1) + e_1 * (fs_4_147 * h5_p1 - fs_8_63 * h5_p3 - fs_1_30 * r_2 * h3_p1 - fs_1_18 * r_2 * h3_p3 - fs_9_245 * r_4 * h1_p1) - fs_45_16 * e_2 * h1_p1;
    }

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_m1, ph1_0, ph3_m2, ph3_m1, ph3_0, ph3_p3, ph5_m5, ph5_m4, ph5_m2, ph5_m1, ph5_0, ph5_p3, ph5_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];

        pc_29[k] = e_0 * (-fs_9_5 * h3_0 + fs_9_5 * r_2 * h1_0) + e_1 * (-fs_5_441 * h5_0 - fs_1_7 * h5_p4 + fs_4_45 * r_2 * h3_0 - fs_9_245 * r_4 * h1_0) - fs_45_16 * e_2 * h1_0;

        pc_30[k] = e_0 * (fs_9_20 * h3_m1 - fs_27_10 * r_2 * h1_m1) + e_1 * (-fs_5_21 * h5_m5 + fs_1_882 * h5_m1 - fs_1_45 * r_2 * h3_m1 + fs_27_490 * r_4 * h1_m1) + fs_135_32 * e_2 * h1_m1;

        pc_31[k] = -fs_9_8 * e_0 * h3_m2 + e_1 * (-fs_2_21 * h5_m4 - fs_1_126 * h5_m2 + fs_1_18 * r_2 * h3_m2);

        pc_32[k] = -f_3_2 * e_0 * h3_p3 + e_1 * (-fs_4_63 * h5_p3 + f_1_3 * r_2 * h3_p3);
    }

    // NOTE: the rows are formed in 9 loops, as the vectorizer runs out of
    // registers with all 35 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph1_p1, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph5_p4, ph5_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];

        pc_33[k] = fs_9_8 * e_0 * h3_p2 + e_1 * (fs_1_126 * h5_p2 - fs_2_21 * h5_p4 - fs_1_18 * r_2 * h3_p2);

        pc_34[k] = e_0 * (-fs_9_20 * h3_p1 + fs_27_10 * r_2 * h1_p1) + e_1 * (-fs_1_882 * h5_p1 - fs_5_21 * h5_p5 + fs_1_45 * r_2 * h3_p1 - fs_27_490 * r_4 * h1_p1) - fs_135_32 * e_2 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[35] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};

    for (size_t m = 0; m < 35; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
