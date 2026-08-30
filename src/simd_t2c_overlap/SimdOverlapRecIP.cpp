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



#include "SimdOverlapRecIP.hpp"

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
compute_ip_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 1))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIP.compute_ip_overlap: Basis functions must be of angular momenta six and one"));
    }

    if (harmonics.size() < 7)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIP.compute_ip_overlap: Harmonics must reach angular momentum seven"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIP.compute_ip_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 39 * nvalues, 0.0);

        return;
    }

    // NOTE: the first two rows accumulate the contracted prefactors of the terms,
    // and the remaining 39 rows hold the integrals of the combinations of angular
    // components.

    auto buffer = CSimdMatrix(41, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);

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

            const auto f_0 = fbase * bexp * bexp * bexp * bexp * bexp * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * bexp * bexp * bexp * bexp * bexp / fexp / fexp / fexp / fexp / fexp / fexp;

            // NOTE: the exponential depends on the pair of primitives alone, so it is
            // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fss = std::exp(-fmu * ab_2[k]);

                pe_0[k] += f_0 * fss;
                pe_1[k] += f_1 * fss;
            }
        }
    }

    // NOTE: the geometry of a term is a solid harmonic of the vector between the
    // atoms times a power of their squared distance.

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

    auto *pc_0 = buffer.data(2);
    auto *pc_1 = buffer.data(3);
    auto *pc_2 = buffer.data(4);
    auto *pc_3 = buffer.data(5);
    auto *pc_4 = buffer.data(6);
    auto *pc_5 = buffer.data(7);
    auto *pc_6 = buffer.data(8);
    auto *pc_7 = buffer.data(9);
    auto *pc_8 = buffer.data(10);
    auto *pc_9 = buffer.data(11);
    auto *pc_10 = buffer.data(12);
    auto *pc_11 = buffer.data(13);
    auto *pc_12 = buffer.data(14);
    auto *pc_13 = buffer.data(15);
    auto *pc_14 = buffer.data(16);
    auto *pc_15 = buffer.data(17);
    auto *pc_16 = buffer.data(18);
    auto *pc_17 = buffer.data(19);
    auto *pc_18 = buffer.data(20);
    auto *pc_19 = buffer.data(21);
    auto *pc_20 = buffer.data(22);
    auto *pc_21 = buffer.data(23);
    auto *pc_22 = buffer.data(24);
    auto *pc_23 = buffer.data(25);
    auto *pc_24 = buffer.data(26);
    auto *pc_25 = buffer.data(27);
    auto *pc_26 = buffer.data(28);
    auto *pc_27 = buffer.data(29);
    auto *pc_28 = buffer.data(30);
    auto *pc_29 = buffer.data(31);
    auto *pc_30 = buffer.data(32);
    auto *pc_31 = buffer.data(33);
    auto *pc_32 = buffer.data(34);
    auto *pc_33 = buffer.data(35);
    auto *pc_34 = buffer.data(36);
    auto *pc_35 = buffer.data(37);
    auto *pc_36 = buffer.data(38);
    auto *pc_37 = buffer.data(39);
    auto *pc_38 = buffer.data(40);

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_1_338 = std::sqrt(1.0 / 338.0);
    const auto fs_7_26 = std::sqrt(7.0 / 26.0);
    const auto fs_33_169 = std::sqrt(33.0 / 169.0);
    const auto fs_33_4 = std::sqrt(8.25);
    const auto fs_1_13 = std::sqrt(1.0 / 13.0);
    const auto fs_3_338 = std::sqrt(3.0 / 338.0);
    const auto fs_3_13 = std::sqrt(3.0 / 13.0);
    const auto fs_55_338 = std::sqrt(55.0 / 338.0);
    const auto fs_55_8 = std::sqrt(6.875);
    const auto fs_24_169 = std::sqrt(24.0 / 169.0);
    const auto fs_11_169 = std::sqrt(11.0 / 169.0);
    const auto fs_11_4 = std::sqrt(2.75);
    const auto fs_3_169 = std::sqrt(3.0 / 169.0);
    const auto fs_45_338 = std::sqrt(45.0 / 338.0);
    const auto fs_45_8 = std::sqrt(5.625);
    const auto fs_1_8 = std::sqrt(0.125);
    const auto fs_20_169 = std::sqrt(20.0 / 169.0);
    const auto fs_5 = std::sqrt(5.0);
    const auto fs_5_169 = std::sqrt(5.0 / 169.0);
    const auto fs_18_169 = std::sqrt(18.0 / 169.0);
    const auto fs_9_2 = std::sqrt(4.5);
    const auto fs_3_8 = std::sqrt(0.375);
    const auto fs_40_169 = std::sqrt(40.0 / 169.0);
    const auto fs_27_169 = std::sqrt(27.0 / 169.0);
    const auto fs_27_4 = std::sqrt(6.75);
    const auto fs_15_338 = std::sqrt(15.0 / 338.0);
    const auto fs_14_169 = std::sqrt(14.0 / 169.0);
    const auto fs_7_2 = std::sqrt(3.5);
    const auto fs_3_4 = std::sqrt(0.75);
    const auto fs_45_169 = std::sqrt(45.0 / 169.0);
    const auto fs_32_169 = std::sqrt(32.0 / 169.0);
    const auto fs_8 = std::sqrt(8.0);
    const auto fs_21_169 = std::sqrt(21.0 / 169.0);
    const auto fs_21_4 = std::sqrt(5.25);
    const auto fs_5_4 = std::sqrt(1.25);
    const auto fs_48_169 = std::sqrt(48.0 / 169.0);
    const auto fs_35_169 = std::sqrt(35.0 / 169.0);
    const auto fs_35_4 = std::sqrt(8.75);
    const auto fs_28_169 = std::sqrt(28.0 / 169.0);
    const auto fs_15_169 = std::sqrt(15.0 / 169.0);
    const auto fs_15_4 = std::sqrt(3.75);
    const auto f_7_13 = 7.0 / 13.0;
    const auto f_6_13 = 6.0 / 13.0;
    const auto f_3 = 3.0;

    // NOTE: the rows are formed in 6 loops, as the vectorizer runs out of
    // registers with all 39 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph5_m5, ph5_m4, ph5_p3, ph5_p4, ph5_p5, ph7_m7, ph7_m6, ph7_m5, ph7_m4, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ph7_p7, ab_2, pc_0, pc_1, pc_2, pc_3, pc_4, pc_5, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_0[k] = e_0 * (-fs_1_338 * h7_p5 - fs_7_26 * h7_p7 + fs_33_169 * r_2 * h5_p5) - fs_33_4 * e_1 * h5_p5;

        pc_1[k] = fs_1_13 * e_0 * h7_m6;

        pc_2[k] = e_0 * (fs_7_26 * h7_m7 - fs_1_338 * h7_m5 + fs_33_169 * r_2 * h5_m5) - fs_33_4 * e_1 * h5_m5;

        pc_3[k] = e_0 * (-fs_3_338 * h7_p4 - fs_3_13 * h7_p6 + fs_55_338 * r_2 * h5_p4) - fs_55_8 * e_1 * h5_p4;

        pc_4[k] = e_0 * (fs_24_169 * h7_m5 + fs_11_169 * r_2 * h5_m5) - fs_11_4 * e_1 * h5_m5;

        pc_5[k] = e_0 * (fs_3_13 * h7_m6 - fs_3_338 * h7_m4 + fs_55_338 * r_2 * h5_m4) - fs_55_8 * e_1 * h5_m4;

        pc_6[k] = e_0 * (-fs_3_169 * h7_p3 - fs_33_169 * h7_p5 + fs_45_338 * r_2 * h5_p3 + fs_1_338 * r_2 * h5_p5) + e_1 * (-fs_45_8 * h5_p3 - fs_1_8 * h5_p5);

        pc_7[k] = e_0 * (fs_33_169 * h7_m4 + fs_20_169 * r_2 * h5_m4) - fs_5 * e_1 * h5_m4;
    }

    // NOTE: the rows are formed in 6 loops, as the vectorizer runs out of
    // registers with all 39 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph5_p2, ph5_p4, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ph7_p2, ph7_p4, ab_2, pc_8, pc_9, pc_10, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];

        pc_8[k] = e_0 * (fs_33_169 * h7_m5 - fs_3_169 * h7_m3 - fs_1_338 * r_2 * h5_m5 + fs_45_338 * r_2 * h5_m3) + e_1 * (fs_1_8 * h5_m5 - fs_45_8 * h5_m3);

        pc_9[k] = e_0 * (-fs_5_169 * h7_p2 - fs_55_338 * h7_p4 + fs_18_169 * r_2 * h5_p2 + fs_3_338 * r_2 * h5_p4) + e_1 * (-fs_9_2 * h5_p2 - fs_3_8 * h5_p4);

        pc_10[k] = e_0 * (fs_40_169 * h7_m3 + fs_27_169 * r_2 * h5_m3) - fs_27_4 * e_1 * h5_m3;

        pc_11[k] = e_0 * (fs_55_338 * h7_m4 - fs_5_169 * h7_m2 - fs_3_338 * r_2 * h5_m4 + fs_18_169 * r_2 * h5_m2) + e_1 * (fs_3_8 * h5_m4 - fs_9_2 * h5_m2);
    }

    // NOTE: the rows are formed in 6 loops, as the vectorizer runs out of
    // registers with all 39 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph5_m3, ph5_m2, ph5_m1, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph7_m3, ph7_m2, ph7_m1, ph7_0, ph7_p1, ph7_p2, ph7_p3, ab_2, pc_12, pc_13, pc_14, pc_15, pc_16, pc_17, pc_18, pc_19, pc_20, pc_21, pc_22, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];

        pc_12[k] = e_0 * (-fs_15_338 * h7_p1 - fs_45_338 * h7_p3 + fs_14_169 * r_2 * h5_p1 + fs_3_169 * r_2 * h5_p3) + e_1 * (-fs_7_2 * h5_p1 - fs_3_4 * h5_p3);

        pc_13[k] = e_0 * (fs_45_169 * h7_m2 + fs_32_169 * r_2 * h5_m2) - fs_8 * e_1 * h5_m2;

        pc_14[k] = e_0 * (fs_45_338 * h7_m3 - fs_15_338 * h7_m1 - fs_3_169 * r_2 * h5_m3 + fs_14_169 * r_2 * h5_m1) + e_1 * (fs_3_4 * h5_m3 - fs_7_2 * h5_m1);

        pc_15[k] = e_0 * (-fs_21_169 * h7_0 - fs_18_169 * h7_p2 + fs_21_169 * r_2 * h5_0 + fs_5_169 * r_2 * h5_p2) + e_1 * (-fs_21_4 * h5_0 - fs_5_4 * h5_p2);

        pc_16[k] = e_0 * (fs_48_169 * h7_m1 + fs_35_169 * r_2 * h5_m1) - fs_35_4 * e_1 * h5_m1;

        pc_17[k] = e_0 * (fs_18_169 * h7_m2 - fs_5_169 * r_2 * h5_m2) + fs_5_4 * e_1 * h5_m2;

        pc_18[k] = e_0 * (fs_28_169 * h7_m1 - fs_15_169 * r_2 * h5_m1) + fs_15_4 * e_1 * h5_m1;

        pc_19[k] = e_0 * (f_7_13 * h7_0 + f_6_13 * r_2 * h5_0) - f_3 * e_1 * h5_0;

        pc_20[k] = e_0 * (fs_28_169 * h7_p1 - fs_15_169 * r_2 * h5_p1) + fs_15_4 * e_1 * h5_p1;

        pc_21[k] = e_0 * (fs_18_169 * h7_m2 - fs_5_169 * r_2 * h5_m2) + fs_5_4 * e_1 * h5_m2;

        pc_22[k] = e_0 * (fs_48_169 * h7_p1 + fs_35_169 * r_2 * h5_p1) - fs_35_4 * e_1 * h5_p1;

        pc_23[k] = e_0 * (-fs_21_169 * h7_0 + fs_18_169 * h7_p2 + fs_21_169 * r_2 * h5_0 - fs_5_169 * r_2 * h5_p2) + e_1 * (-fs_21_4 * h5_0 + fs_5_4 * h5_p2);
    }

    // NOTE: the rows are formed in 6 loops, as the vectorizer runs out of
    // registers with all 39 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph5_p1, ph5_p2, ph5_p3, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ph7_p1, ph7_p2, ph7_p3, ab_2, pc_24, pc_25, pc_26, pc_27, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];

        pc_24[k] = e_0 * (fs_45_338 * h7_m3 + fs_15_338 * h7_m1 - fs_3_169 * r_2 * h5_m3 - fs_14_169 * r_2 * h5_m1) + e_1 * (fs_3_4 * h5_m3 + fs_7_2 * h5_m1);

        pc_25[k] = e_0 * (fs_45_169 * h7_p2 + fs_32_169 * r_2 * h5_p2) - fs_8 * e_1 * h5_p2;

        pc_26[k] = e_0 * (-fs_15_338 * h7_p1 + fs_45_338 * h7_p3 + fs_14_169 * r_2 * h5_p1 - fs_3_169 * r_2 * h5_p3) + e_1 * (-fs_7_2 * h5_p1 + fs_3_4 * h5_p3);

        pc_27[k] = e_0 * (fs_55_338 * h7_m4 + fs_5_169 * h7_m2 - fs_3_338 * r_2 * h5_m4 - fs_18_169 * r_2 * h5_m2) + e_1 * (fs_3_8 * h5_m4 + fs_9_2 * h5_m2);

        pc_28[k] = e_0 * (fs_40_169 * h7_p3 + fs_27_169 * r_2 * h5_p3) - fs_27_4 * e_1 * h5_p3;
    }

    // NOTE: the rows are formed in 6 loops, as the vectorizer runs out of
    // registers with all 39 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph5_m5, ph5_m4, ph5_m3, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ab_2, pc_29, pc_30, pc_31, pc_32, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_29[k] = e_0 * (-fs_5_169 * h7_p2 + fs_55_338 * h7_p4 + fs_18_169 * r_2 * h5_p2 - fs_3_338 * r_2 * h5_p4) + e_1 * (-fs_9_2 * h5_p2 + fs_3_8 * h5_p4);

        pc_30[k] = e_0 * (fs_33_169 * h7_m5 + fs_3_169 * h7_m3 - fs_1_338 * r_2 * h5_m5 - fs_45_338 * r_2 * h5_m3) + e_1 * (fs_1_8 * h5_m5 + fs_45_8 * h5_m3);

        pc_31[k] = e_0 * (fs_33_169 * h7_p4 + fs_20_169 * r_2 * h5_p4) - fs_5 * e_1 * h5_p4;

        pc_32[k] = e_0 * (-fs_3_169 * h7_p3 + fs_33_169 * h7_p5 + fs_45_338 * r_2 * h5_p3 - fs_1_338 * r_2 * h5_p5) + e_1 * (-fs_45_8 * h5_p3 + fs_1_8 * h5_p5);

        pc_33[k] = e_0 * (fs_3_13 * h7_m6 + fs_3_338 * h7_m4 - fs_55_338 * r_2 * h5_m4) + fs_55_8 * e_1 * h5_m4;

        pc_34[k] = e_0 * (fs_24_169 * h7_p5 + fs_11_169 * r_2 * h5_p5) - fs_11_4 * e_1 * h5_p5;
    }

    // NOTE: the rows are formed in 6 loops, as the vectorizer runs out of
    // registers with all 39 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph5_m5, ph5_p4, ph5_p5, ph7_m7, ph7_m5, ph7_p4, ph7_p5, ph7_p6, ph7_p7, ab_2, pc_35, pc_36, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h5_m5 = ph5_m5[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_35[k] = e_0 * (-fs_3_338 * h7_p4 + fs_3_13 * h7_p6 + fs_55_338 * r_2 * h5_p4) - fs_55_8 * e_1 * h5_p4;

        pc_36[k] = e_0 * (fs_7_26 * h7_m7 + fs_1_338 * h7_m5 - fs_33_169 * r_2 * h5_m5) + fs_33_4 * e_1 * h5_m5;

        pc_37[k] = fs_1_13 * e_0 * h7_p6;

        pc_38[k] = e_0 * (-fs_1_338 * h7_p5 + fs_7_26 * h7_p7 + fs_33_169 * r_2 * h5_p5) - fs_33_4 * e_1 * h5_p5;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest, and
    // the atom pairs beyond the reach of every pair of primitives are set to zero.

    const size_t sources[39] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};

    for (size_t m = 0; m < 39; m++)
    {
        const auto *pc = buffer.data(2 + sources[m]);

        std::copy(pc, pc + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdovl
