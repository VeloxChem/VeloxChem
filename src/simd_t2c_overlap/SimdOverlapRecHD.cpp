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



#include "SimdOverlapRecHD.hpp"

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
compute_hd_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 2))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHD.compute_hd_overlap: Basis functions must be of angular momenta five and two"));
    }

    if (harmonics.size() < 7)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHD.compute_hd_overlap: Harmonics must reach angular momentum seven"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHD.compute_hd_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 55 * nvalues, 0.0);

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

            const auto f_0 = fbase * bexp * bexp * bexp * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * bexp * bexp * bexp * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * bexp * bexp * bexp / fexp / fexp / fexp / fexp / fexp;

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_5_12 = std::sqrt(5.0 / 12.0);
    const auto fs_35_3 = std::sqrt(35.0 / 3.0);
    const auto fs_9_40898 = std::sqrt(9.0 / 40898.0);
    const auto fs_63_286 = std::sqrt(63.0 / 286.0);
    const auto fs_5_507 = std::sqrt(5.0 / 507.0);
    const auto fs_35_363 = std::sqrt(35.0 / 363.0);
    const auto fs_945_16 = std::sqrt(59.0625);
    const auto fs_15_8 = std::sqrt(1.875);
    const auto fs_9_3718 = std::sqrt(9.0 / 3718.0);
    const auto fs_9_143 = std::sqrt(9.0 / 143.0);
    const auto fs_15_338 = std::sqrt(15.0 / 338.0);
    const auto f_5_2 = 2.5;
    const auto fs_54_1859 = std::sqrt(54.0 / 1859.0);
    const auto f_5_13 = 5.0 / 13.0;
    const auto f_1 = 1.0;
    const auto fs_7 = std::sqrt(7.0);
    const auto fs_45_40898 = std::sqrt(45.0 / 40898.0);
    const auto fs_45_286 = std::sqrt(45.0 / 286.0);
    const auto f_2_13 = 2.0 / 13.0;
    const auto fs_7_121 = std::sqrt(7.0 / 121.0);
    const auto fs_567_16 = std::sqrt(35.4375);
    const auto fs_49_24 = std::sqrt(49.0 / 24.0);
    const auto fs_14_3 = std::sqrt(14.0 / 3.0);
    const auto fs_180_20449 = std::sqrt(180.0 / 20449.0);
    const auto fs_180_1859 = std::sqrt(180.0 / 1859.0);
    const auto fs_49_1014 = std::sqrt(49.0 / 1014.0);
    const auto fs_14_363 = std::sqrt(14.0 / 363.0);
    const auto fs_189_8 = std::sqrt(23.625);
    const auto fs_135_1859 = std::sqrt(135.0 / 1859.0);
    const auto fs_14_9 = std::sqrt(14.0 / 9.0);
    const auto fs_35_9 = std::sqrt(35.0 / 9.0);
    const auto fs_135_40898 = std::sqrt(135.0 / 40898.0);
    const auto fs_405_3718 = std::sqrt(405.0 / 3718.0);
    const auto fs_56_1521 = std::sqrt(56.0 / 1521.0);
    const auto fs_35_1089 = std::sqrt(35.0 / 1089.0);
    const auto fs_315_16 = std::sqrt(19.6875);
    const auto fs_25_18 = std::sqrt(25.0 / 18.0);
    const auto fs_56_9 = std::sqrt(56.0 / 9.0);
    const auto fs_405_20449 = std::sqrt(405.0 / 20449.0);
    const auto fs_50_1521 = std::sqrt(50.0 / 1521.0);
    const auto fs_56_1089 = std::sqrt(56.0 / 1089.0);
    const auto fs_63_2 = std::sqrt(31.5);
    const auto f_1_6 = 1.0 / 6.0;
    const auto fs_28_9 = std::sqrt(28.0 / 9.0);
    const auto fs_2430_20449 = std::sqrt(2430.0 / 20449.0);
    const auto f_1_39 = 1.0 / 39.0;
    const auto fs_28_1089 = std::sqrt(28.0 / 1089.0);
    const auto fs_63_4 = std::sqrt(15.75);
    const auto fs_315_20449 = std::sqrt(315.0 / 20449.0);
    const auto fs_140_1521 = std::sqrt(140.0 / 1521.0);
    const auto fs_7_12 = std::sqrt(7.0 / 12.0);
    const auto fs_35_6 = std::sqrt(35.0 / 6.0);
    const auto fs_7_18 = std::sqrt(7.0 / 18.0);
    const auto fs_720_20449 = std::sqrt(720.0 / 20449.0);
    const auto fs_2160_20449 = std::sqrt(2160.0 / 20449.0);
    const auto fs_7_507 = std::sqrt(7.0 / 507.0);
    const auto fs_35_726 = std::sqrt(35.0 / 726.0);
    const auto fs_7_2178 = std::sqrt(7.0 / 2178.0);
    const auto fs_945_32 = std::sqrt(29.53125);
    const auto fs_63_32 = std::sqrt(1.96875);
    const auto fs_3240_20449 = std::sqrt(3240.0 / 20449.0);
    const auto fs_25_12 = std::sqrt(25.0 / 12.0);
    const auto fs_5_6 = std::sqrt(5.0 / 6.0);
    const auto fs_1_18 = std::sqrt(1.0 / 18.0);
    const auto fs_945_20449 = std::sqrt(945.0 / 20449.0);
    const auto fs_25_507 = std::sqrt(25.0 / 507.0);
    const auto fs_5_726 = std::sqrt(5.0 / 726.0);
    const auto fs_1_2178 = std::sqrt(1.0 / 2178.0);
    const auto fs_135_32 = std::sqrt(4.21875);
    const auto fs_9_32 = std::sqrt(0.28125);
    const auto fs_5_36 = std::sqrt(5.0 / 36.0);
    const auto fs_80_9 = std::sqrt(80.0 / 9.0);
    const auto fs_4_3 = std::sqrt(4.0 / 3.0);
    const auto fs_2205_20449 = std::sqrt(2205.0 / 20449.0);
    const auto fs_1890_20449 = std::sqrt(1890.0 / 20449.0);
    const auto fs_5_1521 = std::sqrt(5.0 / 1521.0);
    const auto fs_80_1089 = std::sqrt(80.0 / 1089.0);
    const auto fs_4_363 = std::sqrt(4.0 / 363.0);
    const auto fs_45 = std::sqrt(45.0);
    const auto fs_27_4 = std::sqrt(6.75);
    const auto f_3_2 = 1.5;
    const auto fs_10 = std::sqrt(10.0);
    const auto fs_3780_20449 = std::sqrt(3780.0 / 20449.0);
    const auto f_3_13 = 3.0 / 13.0;
    const auto fs_10_121 = std::sqrt(10.0 / 121.0);
    const auto fs_405_8 = std::sqrt(50.625);
    const auto fs_5_9 = std::sqrt(5.0 / 9.0);
    const auto fs_1134_20449 = std::sqrt(1134.0 / 20449.0);
    const auto fs_5_1089 = std::sqrt(5.0 / 1089.0);
    const auto fs_45_16 = std::sqrt(2.8125);
    const auto fs_50_9 = std::sqrt(50.0 / 9.0);
    const auto fs_3024_20449 = std::sqrt(3024.0 / 20449.0);
    const auto fs_50_1089 = std::sqrt(50.0 / 1089.0);
    const auto fs_225_8 = std::sqrt(28.125);
    const auto f_5_3 = 5.0 / 3.0;
    const auto f_10_3 = 10.0 / 3.0;
    const auto f_63_143 = 63.0 / 143.0;
    const auto f_10_39 = 10.0 / 39.0;
    const auto f_10_33 = 10.0 / 33.0;
    const auto f_15_2 = 7.5;

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p3, ph5_m5, ph5_m4, ph5_p3, ph5_p4, ph7_m6, ph7_m5, ph7_m4, ph7_p3, ph7_p4, ph7_p6, ph7_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_0[k] = e_0 * (-fs_5_12 * h5_p3 + fs_35_3 * r_2 * h3_p3) + e_1 * (-fs_9_40898 * h7_p3 + fs_63_286 * h7_p7 + fs_5_507 * r_2 * h5_p3 - fs_35_363 * r_4 * h3_p3) - fs_945_16 * e_2 * h3_p3;

        pc_1[k] = fs_15_8 * e_0 * h5_p4 + e_1 * (fs_9_3718 * h7_p4 + fs_9_143 * h7_p6 - fs_15_338 * r_2 * h5_p4);

        pc_2[k] = -f_5_2 * e_0 * h5_m5 + e_1 * (-fs_54_1859 * h7_m5 + f_5_13 * r_2 * h5_m5);

        pc_3[k] = fs_15_8 * e_0 * h5_m4 + e_1 * (-fs_9_143 * h7_m6 + fs_9_3718 * h7_m4 - fs_15_338 * r_2 * h5_m4);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p2, ph3_p3, ph5_m3, ph5_p2, ph5_p3, ph5_p5, ph7_m7, ph7_m3, ph7_p2, ph7_p3, ph7_p5, ph7_p6, ab_2 : simd::cache_line_size())
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
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_4[k] = e_0 * (-fs_5_12 * h5_m3 + fs_35_3 * r_2 * h3_m3) + e_1 * (-fs_63_286 * h7_m7 - fs_9_40898 * h7_m3 + fs_5_507 * r_2 * h5_m3 - fs_35_363 * r_4 * h3_m3) - fs_945_16 * e_2 * h3_m3;

        pc_5[k] = e_0 * (-f_1 * h5_p2 + fs_7 * r_2 * h3_p2) + e_1 * (-fs_45_40898 * h7_p2 + fs_45_286 * h7_p6 + f_2_13 * r_2 * h5_p2 - fs_7_121 * r_4 * h3_p2) - fs_567_16 * e_2 * h3_p2;

        pc_6[k] = e_0 * (fs_49_24 * h5_p3 - fs_15_8 * h5_p5 + fs_14_3 * r_2 * h3_p3) + e_1 * (fs_180_20449 * h7_p3 + fs_180_1859 * h7_p5 - fs_49_1014 * r_2 * h5_p3 + fs_15_338 * r_2 * h5_p5 - fs_14_363 * r_4 * h3_p3) - fs_189_8 * e_2 * h3_p3;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ab_2 : simd::cache_line_size())
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
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];

        pc_7[k] = -f_1 * e_0 * h5_m4 + e_1 * (-fs_135_1859 * h7_m4 + f_2_13 * r_2 * h5_m4);

        pc_8[k] = e_0 * (fs_15_8 * h5_m5 + fs_49_24 * h5_m3 + fs_14_3 * r_2 * h3_m3) + e_1 * (-fs_180_1859 * h7_m5 + fs_180_20449 * h7_m3 - fs_15_338 * r_2 * h5_m5 - fs_49_1014 * r_2 * h5_m3 - fs_14_363 * r_4 * h3_m3) - fs_189_8 * e_2 * h3_m3;

        pc_9[k] = e_0 * (-f_1 * h5_m2 + fs_7 * r_2 * h3_m2) + e_1 * (-fs_45_286 * h7_m6 - fs_45_40898 * h7_m2 + f_2_13 * r_2 * h5_m2 - fs_7_121 * r_4 * h3_m2) - fs_567_16 * e_2 * h3_m2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p1, ph3_p2, ph5_m3, ph5_p1, ph5_p2, ph5_p4, ph5_p5, ph7_m3, ph7_p1, ph7_p2, ph7_p4, ph7_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_10[k] = e_0 * (-fs_14_9 * h5_p1 + fs_5_12 * h5_p5 + fs_35_9 * r_2 * h3_p1) + e_1 * (-fs_135_40898 * h7_p1 + fs_405_3718 * h7_p5 + fs_56_1521 * r_2 * h5_p1 - fs_5_507 * r_2 * h5_p5 - fs_35_1089 * r_4 * h3_p1) - fs_315_16 * e_2 * h3_p1;

        pc_11[k] = e_0 * (fs_25_18 * h5_p2 - fs_49_24 * h5_p4 + fs_56_9 * r_2 * h3_p2) + e_1 * (fs_405_20449 * h7_p2 + fs_405_3718 * h7_p4 - fs_50_1521 * r_2 * h5_p2 + fs_49_1014 * r_2 * h5_p4 - fs_56_1089 * r_4 * h3_p2) - fs_63_2 * e_2 * h3_p2;

        pc_12[k] = e_0 * (f_1_6 * h5_m3 + fs_28_9 * r_2 * h3_m3) + e_1 * (-fs_2430_20449 * h7_m3 - f_1_39 * r_2 * h5_m3 - fs_28_1089 * r_4 * h3_m3) - fs_63_4 * e_2 * h3_m3;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_m1, ph5_m5, ph5_m4, ph5_m2, ph5_m1, ph7_m5, ph7_m4, ph7_m2, ph7_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_13[k] = e_0 * (fs_49_24 * h5_m4 + fs_25_18 * h5_m2 + fs_56_9 * r_2 * h3_m2) + e_1 * (-fs_405_3718 * h7_m4 + fs_405_20449 * h7_m2 - fs_49_1014 * r_2 * h5_m4 - fs_50_1521 * r_2 * h5_m2 - fs_56_1089 * r_4 * h3_m2) - fs_63_2 * e_2 * h3_m2;

        pc_14[k] = e_0 * (-fs_5_12 * h5_m5 - fs_14_9 * h5_m1 + fs_35_9 * r_2 * h3_m1) + e_1 * (-fs_405_3718 * h7_m5 - fs_135_40898 * h7_m1 + fs_5_507 * r_2 * h5_m5 + fs_56_1521 * r_2 * h5_m1 - fs_35_1089 * r_4 * h3_m1) - fs_315_16 * e_2 * h3_m1;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_0, ph3_p1, ph3_p3, ph5_m2, ph5_0, ph5_p1, ph5_p3, ph5_p4, ph7_m2, ph7_0, ph7_p1, ph7_p3, ph7_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_15[k] = e_0 * (-fs_35_9 * h5_0 + f_1 * h5_p4 + fs_35_9 * r_2 * h3_0) + e_1 * (-fs_315_20449 * h7_0 + fs_135_1859 * h7_p4 + fs_140_1521 * r_2 * h5_0 - f_2_13 * r_2 * h5_p4 - fs_35_1089 * r_4 * h3_0) - fs_315_16 * e_2 * h3_0;

        pc_16[k] = e_0 * (fs_7_12 * h5_p1 - fs_25_18 * h5_p3 + fs_35_6 * r_2 * h3_p1 + fs_7_18 * r_2 * h3_p3) + e_1 * (fs_720_20449 * h7_p1 + fs_2160_20449 * h7_p3 - fs_7_507 * r_2 * h5_p1 + fs_50_1521 * r_2 * h5_p3 - fs_35_726 * r_4 * h3_p1 - fs_7_2178 * r_4 * h3_p3) + e_2 * (-fs_945_32 * h3_p1 - fs_63_32 * h3_p3);

        pc_17[k] = e_0 * (f_1 * h5_m2 + fs_7 * r_2 * h3_m2) + e_1 * (-fs_3240_20449 * h7_m2 - f_2_13 * r_2 * h5_m2 - fs_7_121 * r_4 * h3_m2) - fs_567_16 * e_2 * h3_m2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m1, ph3_p1, ph3_p3, ph5_m4, ph5_m3, ph5_m1, ph5_p1, ph5_p3, ph7_m4, ph7_m3, ph7_m1, ph7_p1, ph7_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];

        pc_18[k] = e_0 * (fs_25_18 * h5_m3 + fs_7_12 * h5_m1 - fs_7_18 * r_2 * h3_m3 + fs_35_6 * r_2 * h3_m1) + e_1 * (-fs_2160_20449 * h7_m3 + fs_720_20449 * h7_m1 - fs_50_1521 * r_2 * h5_m3 - fs_7_507 * r_2 * h5_m1 + fs_7_2178 * r_4 * h3_m3 - fs_35_726 * r_4 * h3_m1) + e_2 * (fs_63_32 * h3_m3 - fs_945_32 * h3_m1);

        pc_19[k] = -f_1 * e_0 * h5_m4 + e_1 * (-fs_135_1859 * h7_m4 + f_2_13 * r_2 * h5_m4);

        pc_20[k] = e_0 * (fs_25_12 * h5_p1 + fs_14_9 * h5_p3 - fs_5_6 * r_2 * h3_p1 - fs_1_18 * r_2 * h3_p3) + e_1 * (fs_315_20449 * h7_p1 + fs_945_20449 * h7_p3 - fs_25_507 * r_2 * h5_p1 - fs_56_1521 * r_2 * h5_p3 + fs_5_726 * r_4 * h3_p1 + fs_1_2178 * r_4 * h3_p3) + e_2 * (fs_135_32 * h3_p1 + fs_9_32 * h3_p3);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_m1, ph3_0, ph3_p2, ph5_m2, ph5_m1, ph5_0, ph5_p2, ph7_m2, ph7_m1, ph7_0, ph7_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];

        pc_21[k] = e_0 * (fs_5_36 * h5_0 - fs_7_12 * h5_p2 + fs_80_9 * r_2 * h3_0 + fs_4_3 * r_2 * h3_p2) + e_1 * (fs_2205_20449 * h7_0 + fs_1890_20449 * h7_p2 - fs_5_1521 * r_2 * h5_0 + fs_7_507 * r_2 * h5_p2 - fs_80_1089 * r_4 * h3_0 - fs_4_363 * r_4 * h3_p2) + e_2 * (-fs_45 * h3_0 - fs_27_4 * h3_p2);

        pc_22[k] = e_0 * (f_3_2 * h5_m1 + fs_10 * r_2 * h3_m1) + e_1 * (-fs_3780_20449 * h7_m1 - f_3_13 * r_2 * h5_m1 - fs_10_121 * r_4 * h3_m1) - fs_405_8 * e_2 * h3_m1;

        pc_23[k] = e_0 * (fs_7_12 * h5_m2 - fs_4_3 * r_2 * h3_m2) + e_1 * (-fs_1890_20449 * h7_m2 - fs_7_507 * r_2 * h5_m2 + fs_4_363 * r_4 * h3_m2) + fs_27_4 * e_2 * h3_m2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph3_m1, ph3_0, ph5_m3, ph5_m2, ph5_m1, ph5_0, ph7_m3, ph7_m2, ph7_m1, ph7_0, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];

        pc_24[k] = e_0 * (-fs_14_9 * h5_m3 - fs_25_12 * h5_m1 + fs_1_18 * r_2 * h3_m3 + fs_5_6 * r_2 * h3_m1) + e_1 * (-fs_945_20449 * h7_m3 - fs_315_20449 * h7_m1 + fs_56_1521 * r_2 * h5_m3 + fs_25_507 * r_2 * h5_m1 - fs_1_2178 * r_4 * h3_m3 - fs_5_726 * r_4 * h3_m1) + e_2 * (-fs_9_32 * h3_m3 - fs_135_32 * h3_m1);

        pc_25[k] = e_0 * (-fs_35_9 * h5_m2 + fs_5_9 * r_2 * h3_m2) + e_1 * (-fs_1134_20449 * h7_m2 + fs_140_1521 * r_2 * h5_m2 - fs_5_1089 * r_4 * h3_m2) - fs_45_16 * e_2 * h3_m2;

        pc_26[k] = e_0 * (fs_5_36 * h5_m1 - fs_50_9 * r_2 * h3_m1) + e_1 * (-fs_3024_20449 * h7_m1 - fs_5_1521 * r_2 * h5_m1 + fs_50_1089 * r_4 * h3_m1) + fs_225_8 * e_2 * h3_m1;

        pc_27[k] = e_0 * (f_5_3 * h5_0 + f_10_3 * r_2 * h3_0) + e_1 * (-f_63_143 * h7_0 - f_10_39 * r_2 * h5_0 - f_10_33 * r_4 * h3_0) - f_15_2 * e_2 * h3_0;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m1, ph3_p1, ph3_p2, ph5_m3, ph5_m1, ph5_p1, ph5_p2, ph7_m3, ph7_m1, ph7_p1, ph7_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_28[k] = e_0 * (fs_5_36 * h5_p1 - fs_50_9 * r_2 * h3_p1) + e_1 * (-fs_3024_20449 * h7_p1 - fs_5_1521 * r_2 * h5_p1 + fs_50_1089 * r_4 * h3_p1) + fs_225_8 * e_2 * h3_p1;

        pc_29[k] = e_0 * (-fs_35_9 * h5_p2 + fs_5_9 * r_2 * h3_p2) + e_1 * (-fs_1134_20449 * h7_p2 + fs_140_1521 * r_2 * h5_p2 - fs_5_1089 * r_4 * h3_p2) - fs_45_16 * e_2 * h3_p2;

        pc_30[k] = e_0 * (-fs_14_9 * h5_m3 + fs_25_12 * h5_m1 + fs_1_18 * r_2 * h3_m3 - fs_5_6 * r_2 * h3_m1) + e_1 * (-fs_945_20449 * h7_m3 + fs_315_20449 * h7_m1 + fs_56_1521 * r_2 * h5_m3 - fs_25_507 * r_2 * h5_m1 - fs_1_2178 * r_4 * h3_m3 + fs_5_726 * r_4 * h3_m1) + e_2 * (-fs_9_32 * h3_m3 + fs_135_32 * h3_m1);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_0, ph3_p1, ph3_p2, ph5_m2, ph5_0, ph5_p1, ph5_p2, ph7_m2, ph7_0, ph7_p1, ph7_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_31[k] = e_0 * (fs_7_12 * h5_m2 - fs_4_3 * r_2 * h3_m2) + e_1 * (-fs_1890_20449 * h7_m2 - fs_7_507 * r_2 * h5_m2 + fs_4_363 * r_4 * h3_m2) + fs_27_4 * e_2 * h3_m2;

        pc_32[k] = e_0 * (f_3_2 * h5_p1 + fs_10 * r_2 * h3_p1) + e_1 * (-fs_3780_20449 * h7_p1 - f_3_13 * r_2 * h5_p1 - fs_10_121 * r_4 * h3_p1) - fs_405_8 * e_2 * h3_p1;

        pc_33[k] = e_0 * (fs_5_36 * h5_0 + fs_7_12 * h5_p2 + fs_80_9 * r_2 * h3_0 - fs_4_3 * r_2 * h3_p2) + e_1 * (fs_2205_20449 * h7_0 - fs_1890_20449 * h7_p2 - fs_5_1521 * r_2 * h5_0 - fs_7_507 * r_2 * h5_p2 - fs_80_1089 * r_4 * h3_0 + fs_4_363 * r_4 * h3_p2) + e_2 * (-fs_45 * h3_0 + fs_27_4 * h3_p2);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m1, ph3_p1, ph3_p3, ph5_m4, ph5_m3, ph5_m1, ph5_p1, ph5_p3, ph7_m4, ph7_m3, ph7_m1, ph7_p1, ph7_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];

        pc_34[k] = e_0 * (fs_25_12 * h5_p1 - fs_14_9 * h5_p3 - fs_5_6 * r_2 * h3_p1 + fs_1_18 * r_2 * h3_p3) + e_1 * (fs_315_20449 * h7_p1 - fs_945_20449 * h7_p3 - fs_25_507 * r_2 * h5_p1 + fs_56_1521 * r_2 * h5_p3 + fs_5_726 * r_4 * h3_p1 - fs_1_2178 * r_4 * h3_p3) + e_2 * (fs_135_32 * h3_p1 - fs_9_32 * h3_p3);

        pc_35[k] = -f_1 * e_0 * h5_m4 + e_1 * (-fs_135_1859 * h7_m4 + f_2_13 * r_2 * h5_m4);

        pc_36[k] = e_0 * (fs_25_18 * h5_m3 - fs_7_12 * h5_m1 - fs_7_18 * r_2 * h3_m3 - fs_35_6 * r_2 * h3_m1) + e_1 * (-fs_2160_20449 * h7_m3 - fs_720_20449 * h7_m1 - fs_50_1521 * r_2 * h5_m3 + fs_7_507 * r_2 * h5_m1 + fs_7_2178 * r_4 * h3_m3 + fs_35_726 * r_4 * h3_m1) + e_2 * (fs_63_32 * h3_m3 + fs_945_32 * h3_m1);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph7_0, ph7_p1, ph7_p2, ph7_p3, ph7_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_37[k] = e_0 * (f_1 * h5_p2 + fs_7 * r_2 * h3_p2) + e_1 * (-fs_3240_20449 * h7_p2 - f_2_13 * r_2 * h5_p2 - fs_7_121 * r_4 * h3_p2) - fs_567_16 * e_2 * h3_p2;

        pc_38[k] = e_0 * (fs_7_12 * h5_p1 + fs_25_18 * h5_p3 + fs_35_6 * r_2 * h3_p1 - fs_7_18 * r_2 * h3_p3) + e_1 * (fs_720_20449 * h7_p1 - fs_2160_20449 * h7_p3 - fs_7_507 * r_2 * h5_p1 - fs_50_1521 * r_2 * h5_p3 - fs_35_726 * r_4 * h3_p1 + fs_7_2178 * r_4 * h3_p3) + e_2 * (-fs_945_32 * h3_p1 + fs_63_32 * h3_p3);

        pc_39[k] = e_0 * (-fs_35_9 * h5_0 - f_1 * h5_p4 + fs_35_9 * r_2 * h3_0) + e_1 * (-fs_315_20449 * h7_0 - fs_135_1859 * h7_p4 + fs_140_1521 * r_2 * h5_0 + f_2_13 * r_2 * h5_p4 - fs_35_1089 * r_4 * h3_0) - fs_315_16 * e_2 * h3_0;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_m1, ph3_p3, ph5_m5, ph5_m4, ph5_m2, ph5_m1, ph5_p3, ph7_m5, ph7_m4, ph7_m2, ph7_m1, ph7_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p3 = ph7_p3[k];

        pc_40[k] = e_0 * (-fs_5_12 * h5_m5 + fs_14_9 * h5_m1 - fs_35_9 * r_2 * h3_m1) + e_1 * (-fs_405_3718 * h7_m5 + fs_135_40898 * h7_m1 + fs_5_507 * r_2 * h5_m5 - fs_56_1521 * r_2 * h5_m1 + fs_35_1089 * r_4 * h3_m1) + fs_315_16 * e_2 * h3_m1;

        pc_41[k] = e_0 * (fs_49_24 * h5_m4 - fs_25_18 * h5_m2 - fs_56_9 * r_2 * h3_m2) + e_1 * (-fs_405_3718 * h7_m4 - fs_405_20449 * h7_m2 - fs_49_1014 * r_2 * h5_m4 + fs_50_1521 * r_2 * h5_m2 + fs_56_1089 * r_4 * h3_m2) + fs_63_2 * e_2 * h3_m2;

        pc_42[k] = e_0 * (f_1_6 * h5_p3 + fs_28_9 * r_2 * h3_p3) + e_1 * (-fs_2430_20449 * h7_p3 - f_1_39 * r_2 * h5_p3 - fs_28_1089 * r_4 * h3_p3) - fs_63_4 * e_2 * h3_p3;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_p1, ph3_p2, ph5_m2, ph5_p1, ph5_p2, ph5_p4, ph5_p5, ph7_m6, ph7_m2, ph7_p1, ph7_p2, ph7_p4, ph7_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_43[k] = e_0 * (fs_25_18 * h5_p2 + fs_49_24 * h5_p4 + fs_56_9 * r_2 * h3_p2) + e_1 * (fs_405_20449 * h7_p2 - fs_405_3718 * h7_p4 - fs_50_1521 * r_2 * h5_p2 - fs_49_1014 * r_2 * h5_p4 - fs_56_1089 * r_4 * h3_p2) - fs_63_2 * e_2 * h3_p2;

        pc_44[k] = e_0 * (-fs_14_9 * h5_p1 - fs_5_12 * h5_p5 + fs_35_9 * r_2 * h3_p1) + e_1 * (-fs_135_40898 * h7_p1 - fs_405_3718 * h7_p5 + fs_56_1521 * r_2 * h5_p1 + fs_5_507 * r_2 * h5_p5 - fs_35_1089 * r_4 * h3_p1) - fs_315_16 * e_2 * h3_p1;

        pc_45[k] = e_0 * (f_1 * h5_m2 - fs_7 * r_2 * h3_m2) + e_1 * (-fs_45_286 * h7_m6 + fs_45_40898 * h7_m2 - f_2_13 * r_2 * h5_m2 + fs_7_121 * r_4 * h3_m2) + fs_567_16 * e_2 * h3_m2;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p3, ph5_m5, ph5_m3, ph5_p3, ph5_p4, ph5_p5, ph7_m5, ph7_m3, ph7_p3, ph7_p4, ph7_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_46[k] = e_0 * (fs_15_8 * h5_m5 - fs_49_24 * h5_m3 - fs_14_3 * r_2 * h3_m3) + e_1 * (-fs_180_1859 * h7_m5 - fs_180_20449 * h7_m3 - fs_15_338 * r_2 * h5_m5 + fs_49_1014 * r_2 * h5_m3 + fs_14_363 * r_4 * h3_m3) + fs_189_8 * e_2 * h3_m3;

        pc_47[k] = -f_1 * e_0 * h5_p4 + e_1 * (-fs_135_1859 * h7_p4 + f_2_13 * r_2 * h5_p4);

        pc_48[k] = e_0 * (fs_49_24 * h5_p3 + fs_15_8 * h5_p5 + fs_14_3 * r_2 * h3_p3) + e_1 * (fs_180_20449 * h7_p3 - fs_180_1859 * h7_p5 - fs_49_1014 * r_2 * h5_p3 - fs_15_338 * r_2 * h5_p5 - fs_14_363 * r_4 * h3_p3) - fs_189_8 * e_2 * h3_p3;
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p2, ph5_m4, ph5_m3, ph5_p2, ph5_p5, ph7_m7, ph7_m6, ph7_m4, ph7_m3, ph7_p2, ph7_p5, ph7_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_49[k] = e_0 * (-f_1 * h5_p2 + fs_7 * r_2 * h3_p2) + e_1 * (-fs_45_40898 * h7_p2 - fs_45_286 * h7_p6 + f_2_13 * r_2 * h5_p2 - fs_7_121 * r_4 * h3_p2) - fs_567_16 * e_2 * h3_p2;

        pc_50[k] = e_0 * (fs_5_12 * h5_m3 - fs_35_3 * r_2 * h3_m3) + e_1 * (-fs_63_286 * h7_m7 + fs_9_40898 * h7_m3 - fs_5_507 * r_2 * h5_m3 + fs_35_363 * r_4 * h3_m3) + fs_945_16 * e_2 * h3_m3;

        pc_51[k] = -fs_15_8 * e_0 * h5_m4 + e_1 * (-fs_9_143 * h7_m6 - fs_9_3718 * h7_m4 + fs_15_338 * r_2 * h5_m4);

        pc_52[k] = -f_5_2 * e_0 * h5_p5 + e_1 * (-fs_54_1859 * h7_p5 + f_5_13 * r_2 * h5_p5);
    }

    // NOTE: the rows are formed in 18 loops, as the vectorizer runs out of
    // registers with all 55 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph7_p6, ph7_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_53[k] = fs_15_8 * e_0 * h5_p4 + e_1 * (fs_9_3718 * h7_p4 - fs_9_143 * h7_p6 - fs_15_338 * r_2 * h5_p4);

        pc_54[k] = e_0 * (-fs_5_12 * h5_p3 + fs_35_3 * r_2 * h3_p3) + e_1 * (-fs_9_40898 * h7_p3 - fs_63_286 * h7_p7 + fs_5_507 * r_2 * h5_p3 - fs_35_363 * r_4 * h3_p3) - fs_945_16 * e_2 * h3_p3;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[55] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54};

    for (size_t m = 0; m < 55; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
