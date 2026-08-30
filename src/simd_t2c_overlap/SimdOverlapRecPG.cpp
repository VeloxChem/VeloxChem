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



#include "SimdOverlapRecPG.hpp"

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
compute_pg_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 1) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecPG.compute_pg_overlap: Basis functions must be of angular momenta one and four"));
    }

    if (harmonics.size() < 5)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecPG.compute_pg_overlap: Harmonics must reach angular momentum five"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecPG.compute_pg_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 27 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

    auto buffer = CSimdMatrix(2, nmax);

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

            const auto f_0 = fbase * aexp * aexp * aexp * fmu / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * aexp * aexp * aexp / fexp / fexp / fexp / fexp;

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_1_162 = std::sqrt(1.0 / 162.0);
    const auto fs_5_18 = std::sqrt(5.0 / 18.0);
    const auto fs_14_81 = std::sqrt(14.0 / 81.0);
    const auto fs_7_2 = std::sqrt(3.5);
    const auto fs_1_54 = std::sqrt(1.0 / 54.0);
    const auto fs_2_9 = std::sqrt(2.0 / 9.0);
    const auto fs_7_54 = std::sqrt(7.0 / 54.0);
    const auto fs_21_8 = std::sqrt(2.625);
    const auto fs_1_27 = std::sqrt(1.0 / 27.0);
    const auto fs_5_54 = std::sqrt(5.0 / 54.0);
    const auto fs_15_8 = std::sqrt(1.875);
    const auto fs_1_8 = std::sqrt(0.125);
    const auto fs_10_81 = std::sqrt(10.0 / 81.0);
    const auto fs_5_2 = std::sqrt(2.5);
    const auto fs_3_8 = std::sqrt(0.375);
    const auto fs_5_27 = std::sqrt(5.0 / 27.0);
    const auto fs_2_27 = std::sqrt(2.0 / 27.0);
    const auto fs_3_2 = std::sqrt(1.5);
    const auto f_1_3 = 1.0 / 3.0;
    const auto f_4_9 = 4.0 / 9.0;
    const auto fs_7_81 = std::sqrt(7.0 / 81.0);
    const auto fs_7_4 = std::sqrt(1.75);
    const auto fs_7_27 = std::sqrt(7.0 / 27.0);
    const auto fs_4_27 = std::sqrt(4.0 / 27.0);
    const auto fs_3 = std::sqrt(3.0);
    const auto fs_8_27 = std::sqrt(8.0 / 27.0);
    const auto fs_15_4 = std::sqrt(3.75);
    const auto f_5_9 = 5.0 / 9.0;
    const auto f_2 = 2.0;

    // NOTE: the rows are formed in 4 loops, as the vectorizer runs out of
    // registers with all 27 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph3_m2, ph3_m1, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_m2, ph5_m1, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];

        pc_0[k] = e_0 * (fs_1_162 * h5_p3 + fs_5_18 * h5_p5 - fs_14_81 * r_2 * h3_p3) + fs_7_2 * e_1 * h3_p3;

        pc_1[k] = e_0 * (fs_1_54 * h5_p2 + fs_2_9 * h5_p4 - fs_7_54 * r_2 * h3_p2) + fs_21_8 * e_1 * h3_p2;

        pc_2[k] = e_0 * (fs_1_27 * h5_p1 + fs_14_81 * h5_p3 - fs_5_54 * r_2 * h3_p1 - fs_1_162 * r_2 * h3_p3) + e_1 * (fs_15_8 * h3_p1 + fs_1_8 * h3_p3);

        pc_3[k] = e_0 * (fs_10_81 * h5_0 + fs_7_54 * h5_p2 - fs_10_81 * r_2 * h3_0 - fs_1_54 * r_2 * h3_p2) + e_1 * (fs_5_2 * h3_0 + fs_3_8 * h3_p2);

        pc_4[k] = e_0 * (-fs_5_27 * h5_m1 + fs_2_27 * r_2 * h3_m1) - fs_3_2 * e_1 * h3_m1;

        pc_5[k] = e_0 * (-fs_7_54 * h5_m2 + fs_1_54 * r_2 * h3_m2) - fs_3_8 * e_1 * h3_m2;
    }

    // NOTE: the rows are formed in 4 loops, as the vectorizer runs out of
    // registers with all 27 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph3_m3, ph3_m2, ph3_m1, ph3_0, ph3_p1, ph3_p2, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph5_0, ph5_p1, ph5_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];

        pc_6[k] = e_0 * (-fs_14_81 * h5_m3 - fs_1_27 * h5_m1 + fs_1_162 * r_2 * h3_m3 + fs_5_54 * r_2 * h3_m1) + e_1 * (-fs_1_8 * h3_m3 - fs_15_8 * h3_m1);

        pc_7[k] = e_0 * (-fs_2_9 * h5_m4 - fs_1_54 * h5_m2 + fs_7_54 * r_2 * h3_m2) - fs_21_8 * e_1 * h3_m2;

        pc_8[k] = e_0 * (-fs_5_18 * h5_m5 - fs_1_162 * h5_m3 + fs_14_81 * r_2 * h3_m3) - fs_7_2 * e_1 * h3_m3;

        pc_9[k] = -f_1_3 * e_0 * h5_m4;

        pc_10[k] = e_0 * (-f_4_9 * h5_m3 - fs_7_81 * r_2 * h3_m3) + fs_7_4 * e_1 * h3_m3;

        pc_11[k] = e_0 * (-fs_7_27 * h5_m2 - fs_4_27 * r_2 * h3_m2) + fs_3 * e_1 * h3_m2;

        pc_12[k] = e_0 * (-fs_8_27 * h5_m1 - fs_5_27 * r_2 * h3_m1) + fs_15_4 * e_1 * h3_m1;

        pc_13[k] = e_0 * (-f_5_9 * h5_0 - f_4_9 * r_2 * h3_0) + f_2 * e_1 * h3_0;

        pc_14[k] = e_0 * (-fs_8_27 * h5_p1 - fs_5_27 * r_2 * h3_p1) + fs_15_4 * e_1 * h3_p1;

        pc_15[k] = e_0 * (-fs_7_27 * h5_p2 - fs_4_27 * r_2 * h3_p2) + fs_3 * e_1 * h3_p2;
    }

    // NOTE: the rows are formed in 4 loops, as the vectorizer runs out of
    // registers with all 27 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph3_m3, ph3_m2, ph3_m1, ph3_p1, ph3_p3, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph5_p1, ph5_p3, ph5_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];

        pc_16[k] = e_0 * (-f_4_9 * h5_p3 - fs_7_81 * r_2 * h3_p3) + fs_7_4 * e_1 * h3_p3;

        pc_17[k] = -f_1_3 * e_0 * h5_p4;

        pc_18[k] = e_0 * (-fs_5_18 * h5_m5 + fs_1_162 * h5_m3 - fs_14_81 * r_2 * h3_m3) + fs_7_2 * e_1 * h3_m3;

        pc_19[k] = e_0 * (-fs_2_9 * h5_m4 + fs_1_54 * h5_m2 - fs_7_54 * r_2 * h3_m2) + fs_21_8 * e_1 * h3_m2;

        pc_20[k] = e_0 * (-fs_14_81 * h5_m3 + fs_1_27 * h5_m1 + fs_1_162 * r_2 * h3_m3 - fs_5_54 * r_2 * h3_m1) + e_1 * (-fs_1_8 * h3_m3 + fs_15_8 * h3_m1);

        pc_21[k] = e_0 * (-fs_7_54 * h5_m2 + fs_1_54 * r_2 * h3_m2) - fs_3_8 * e_1 * h3_m2;

        pc_22[k] = e_0 * (-fs_5_27 * h5_p1 + fs_2_27 * r_2 * h3_p1) - fs_3_2 * e_1 * h3_p1;
    }

    // NOTE: the rows are formed in 4 loops, as the vectorizer runs out of
    // registers with all 27 of them in one.

#pragma omp simd aligned(pe_0, pe_1, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        const auto r_2 = ab_2[k];

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];

        pc_23[k] = e_0 * (fs_10_81 * h5_0 - fs_7_54 * h5_p2 - fs_10_81 * r_2 * h3_0 + fs_1_54 * r_2 * h3_p2) + e_1 * (fs_5_2 * h3_0 - fs_3_8 * h3_p2);

        pc_24[k] = e_0 * (fs_1_27 * h5_p1 - fs_14_81 * h5_p3 - fs_5_54 * r_2 * h3_p1 + fs_1_162 * r_2 * h3_p3) + e_1 * (fs_15_8 * h3_p1 - fs_1_8 * h3_p3);

        pc_25[k] = e_0 * (fs_1_54 * h5_p2 - fs_2_9 * h5_p4 - fs_7_54 * r_2 * h3_p2) + fs_21_8 * e_1 * h3_p2;

        pc_26[k] = e_0 * (fs_1_162 * h5_p3 - fs_5_18 * h5_p5 - fs_14_81 * r_2 * h3_p3) + fs_7_2 * e_1 * h3_p3;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};

    for (size_t m = 0; m < 27; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
