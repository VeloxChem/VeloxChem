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



#include "SimdOverlapRecDD.hpp"

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
compute_dd_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 2))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDD.compute_dd_overlap: Basis functions must be of angular momenta two and two"));
    }

    if (harmonics.size() < 4)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDD.compute_dd_overlap: Harmonics must reach angular momentum four"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDD.compute_dd_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 25 * nvalues, 0.0);

        return;
    }

    // NOTE: the first three rows accumulate the contracted prefactors of the terms,
    // and the remaining 15 rows hold the integrals of the combinations of angular
    // components which are not related by symmetry.

    auto buffer = CSimdMatrix(18, nmax);

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

            const auto f_0 = fbase * fmu / fexp / fexp;

            const auto f_1 = fbase * fmu * fmu / fexp / fexp;

            const auto f_2 = fbase / fexp / fexp;

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_1 = 1.0;
    const auto f_3_35 = 3.0 / 35.0;
    const auto fs_9_35 = std::sqrt(9.0 / 35.0);
    const auto f_2_7 = 2.0 / 7.0;
    const auto f_1_5 = 1.0 / 5.0;
    const auto f_3_4 = 0.75;
    const auto fs_3_4 = std::sqrt(0.75);
    const auto fs_9_490 = std::sqrt(9.0 / 490.0);
    const auto fs_9_70 = std::sqrt(9.0 / 70.0);
    const auto fs_3_49 = std::sqrt(3.0 / 49.0);
    const auto fs_27_245 = std::sqrt(27.0 / 245.0);
    const auto f_1_2 = 0.5;
    const auto f_12_35 = 12.0 / 35.0;
    const auto fs_36_245 = std::sqrt(36.0 / 245.0);
    const auto f_1_7 = 1.0 / 7.0;
    const auto fs_54_245 = std::sqrt(54.0 / 245.0);
    const auto f_18_35 = 18.0 / 35.0;

    // NOTE: the rows are formed in 2 loops, as the vectorizer runs out of
    // registers with all 15 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph2_0, ph2_p1, ph2_p2, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ab_2, pc_0, pc_1, pc_2, pc_3, pc_4, pc_5, pc_6, pc_7, pc_8, pc_9, pc_10, pc_11 : simd::cache_line_size())
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
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];

        pc_0[k] = e_0 * (f_1 * h2_0 - f_1 * r_2) + e_1 * (f_3_35 * h4_0 - fs_9_35 * h4_p4 - f_2_7 * r_2 * h2_0 + f_1_5 * r_4) + f_3_4 * e_2;

        pc_1[k] = -fs_3_4 * e_0 * h2_p1 + e_1 * (-fs_9_490 * h4_p1 - fs_9_70 * h4_p3 + fs_3_49 * r_2 * h2_p1);

        pc_2[k] = f_1 * e_0 * h2_m2 + e_1 * (fs_27_245 * h4_m2 - f_2_7 * r_2 * h2_m2);

        pc_3[k] = -fs_3_4 * e_0 * h2_m1 + e_1 * (fs_9_70 * h4_m3 - fs_9_490 * h4_m1 + fs_3_49 * r_2 * h2_m1);

        pc_4[k] = fs_9_35 * e_1 * h4_m4;

        pc_5[k] = e_0 * (-f_1_2 * h2_0 + fs_3_4 * h2_p2 - f_1 * r_2) + e_1 * (-f_12_35 * h4_0 - fs_36_245 * h4_p2 + f_1_7 * r_2 * h2_0 - fs_3_49 * r_2 * h2_p2 + f_1_5 * r_4) + f_3_4 * e_2;

        pc_6[k] = -f_1_2 * e_0 * h2_m1 + e_1 * (fs_54_245 * h4_m1 + f_1_7 * r_2 * h2_m1);

        pc_7[k] = -fs_3_4 * e_0 * h2_m2 + e_1 * (fs_36_245 * h4_m2 + fs_3_49 * r_2 * h2_m2);

        pc_8[k] = fs_3_4 * e_0 * h2_m1 + e_1 * (fs_9_70 * h4_m3 + fs_9_490 * h4_m1 - fs_3_49 * r_2 * h2_m1);

        pc_9[k] = e_0 * (-f_1 * h2_0 - f_1 * r_2) + e_1 * (f_18_35 * h4_0 + f_2_7 * r_2 * h2_0 + f_1_5 * r_4) + f_3_4 * e_2;

        pc_10[k] = -f_1_2 * e_0 * h2_p1 + e_1 * (fs_54_245 * h4_p1 + f_1_7 * r_2 * h2_p1);

        pc_11[k] = f_1 * e_0 * h2_p2 + e_1 * (fs_27_245 * h4_p2 - f_2_7 * r_2 * h2_p2);
    }

    // NOTE: the rows are formed in 2 loops, as the vectorizer runs out of
    // registers with all 15 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ab_2, pc_12, pc_13, pc_14 : simd::cache_line_size())
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

        pc_12[k] = e_0 * (-f_1_2 * h2_0 - fs_3_4 * h2_p2 - f_1 * r_2) + e_1 * (-f_12_35 * h4_0 + fs_36_245 * h4_p2 + f_1_7 * r_2 * h2_0 + fs_3_49 * r_2 * h2_p2 + f_1_5 * r_4) + f_3_4 * e_2;

        pc_13[k] = -fs_3_4 * e_0 * h2_p1 + e_1 * (-fs_9_490 * h4_p1 + fs_9_70 * h4_p3 + fs_3_49 * r_2 * h2_p1);

        pc_14[k] = e_0 * (f_1 * h2_0 - f_1 * r_2) + e_1 * (f_3_35 * h4_0 + fs_9_35 * h4_p4 - f_2_7 * r_2 * h2_0 + f_1_5 * r_4) + f_3_4 * e_2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest, and
    // the atom pairs beyond the reach of every pair of primitives are set to zero.

    const size_t sources[25] = {0, 1, 2, 3, 4, 1, 5, 6, 7, 8, 2, 6, 9, 10, 11, 3, 7, 10, 12, 13, 4, 8, 11, 13, 14};

    for (size_t m = 0; m < 25; m++)
    {
        const auto *pc = buffer.data(3 + sources[m]);

        std::copy(pc, pc + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdovl
