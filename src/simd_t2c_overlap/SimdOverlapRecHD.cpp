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

    const auto nmax = dimensions.back();

    if (nmax == 0)
    {
        std::fill(values, values + 55 * nvalues, 0.0);

        return;
    }

    // NOTE: the first three rows accumulate the contracted prefactors of the terms,
    // and the remaining 55 rows hold the integrals of the combinations of angular
    // components.

    auto buffer = CSimdMatrix(58, nmax);

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
    auto *pc_45 = buffer.data(48);
    auto *pc_46 = buffer.data(49);
    auto *pc_47 = buffer.data(50);
    auto *pc_48 = buffer.data(51);
    auto *pc_49 = buffer.data(52);
    auto *pc_50 = buffer.data(53);
    auto *pc_51 = buffer.data(54);
    auto *pc_52 = buffer.data(55);
    auto *pc_53 = buffer.data(56);
    auto *pc_54 = buffer.data(57);

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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p3, ph5_m5, ph5_m4, ph5_p3, ph5_p4, ph7_m6, ph7_m5, ph7_m4, ph7_p3, ph7_p4, ph7_p6, ph7_p7, ab_2, pc_0, pc_1, pc_2, pc_3 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p2, ph3_p3, ph5_m3, ph5_p2, ph5_p3, ph5_p5, ph7_m7, ph7_m3, ph7_p2, ph7_p3, ph7_p5, ph7_p6, ab_2, pc_4, pc_5, pc_6 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ab_2, pc_7, pc_8, pc_9 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p1, ph3_p2, ph5_m3, ph5_p1, ph5_p2, ph5_p4, ph5_p5, ph7_m3, ph7_p1, ph7_p2, ph7_p4, ph7_p5, ab_2, pc_10, pc_11, pc_12 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_m1, ph5_m5, ph5_m4, ph5_m2, ph5_m1, ph7_m5, ph7_m4, ph7_m2, ph7_m1, ab_2, pc_13, pc_14 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_0, ph3_p1, ph3_p3, ph5_m2, ph5_0, ph5_p1, ph5_p3, ph5_p4, ph7_m2, ph7_0, ph7_p1, ph7_p3, ph7_p4, ab_2, pc_15, pc_16, pc_17 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m1, ph3_p1, ph3_p3, ph5_m4, ph5_m3, ph5_m1, ph5_p1, ph5_p3, ph7_m4, ph7_m3, ph7_m1, ph7_p1, ph7_p3, ab_2, pc_18, pc_19, pc_20 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_m1, ph3_0, ph3_p2, ph5_m2, ph5_m1, ph5_0, ph5_p2, ph7_m2, ph7_m1, ph7_0, ph7_p2, ab_2, pc_21, pc_22, pc_23 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph3_m1, ph3_0, ph5_m3, ph5_m2, ph5_m1, ph5_0, ph7_m3, ph7_m2, ph7_m1, ph7_0, ab_2, pc_24, pc_25, pc_26, pc_27 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m1, ph3_p1, ph3_p2, ph5_m3, ph5_m1, ph5_p1, ph5_p2, ph7_m3, ph7_m1, ph7_p1, ph7_p2, ab_2, pc_28, pc_29, pc_30 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_0, ph3_p1, ph3_p2, ph5_m2, ph5_0, ph5_p1, ph5_p2, ph7_m2, ph7_0, ph7_p1, ph7_p2, ab_2, pc_31, pc_32, pc_33 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m1, ph3_p1, ph3_p3, ph5_m4, ph5_m3, ph5_m1, ph5_p1, ph5_p3, ph7_m4, ph7_m3, ph7_m1, ph7_p1, ph7_p3, ab_2, pc_34, pc_35, pc_36 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph7_0, ph7_p1, ph7_p2, ph7_p3, ph7_p4, ab_2, pc_37, pc_38, pc_39 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_m1, ph3_p3, ph5_m5, ph5_m4, ph5_m2, ph5_m1, ph5_p3, ph7_m5, ph7_m4, ph7_m2, ph7_m1, ph7_p3, ab_2, pc_40, pc_41, pc_42 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_p1, ph3_p2, ph5_m2, ph5_p1, ph5_p2, ph5_p4, ph5_p5, ph7_m6, ph7_m2, ph7_p1, ph7_p2, ph7_p4, ph7_p5, ab_2, pc_43, pc_44, pc_45 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p3, ph5_m5, ph5_m3, ph5_p3, ph5_p4, ph5_p5, ph7_m5, ph7_m3, ph7_p3, ph7_p4, ph7_p5, ab_2, pc_46, pc_47, pc_48 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_p2, ph5_m4, ph5_m3, ph5_p2, ph5_p5, ph7_m7, ph7_m6, ph7_m4, ph7_m3, ph7_p2, ph7_p5, ph7_p6, ab_2, pc_49, pc_50, pc_51, pc_52 : simd::cache_line_size())
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

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph7_p6, ph7_p7, ab_2, pc_53, pc_54 : simd::cache_line_size())
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
    // row of nvalues columns, with the component on bra side running slowest, and
    // the atom pairs beyond the reach of every pair of primitives are set to zero.

    const size_t sources[55] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54};

    for (size_t m = 0; m < 55; m++)
    {
        const auto *pc = buffer.data(3 + sources[m]);

        std::copy(pc, pc + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdovl
