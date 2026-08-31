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



#include "SimdKineticEnergyRecPF.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_pf_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 1) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_pf_kinetic_energy: Basis functions must be of angular momenta one and three"));
    }

    if (harmonics.size() < 4)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_pf_kinetic_energy: Harmonics must reach angular momentum four"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_pf_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        bra, ket, nvalues, coordinates, screenfunc::two_center_kinetic_energy_primitive_bound,
        threshold / static_cast<double>(nprims));

    // NOTE: the buffer spans the atom pairs reached by the pair of primitives
    // reaching furthest, which is searched for rather than assumed. The
    // primitives are sorted by descending exponent, but the bound of a pair of
    // primitives carries their prefactor as well as their decay, so a tighter
    // pair with a larger prefactor reaches further than a more diffuse pair with
    // a smaller one, and the last pair is not always the furthest reaching.

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 21 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time.

    auto buffer = CSimdMatrix(3, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);

    const auto *ab_2 = coordinates.data(6);

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the prefactor of each exponent factor over the pairs of primitives

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

            const auto ff_0 = fbase * aexp * aexp * fmu / (fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * aexp * fmu * fmu / (fexp * fexp * fexp);

            const auto ff_2 = fbase * aexp * aexp * fmu * fmu * fmu / (fexp * fexp * fexp);

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fterm = std::exp(-fmu * ab_2[k]);

                pe_0[k] += ff_0 * fterm;
                pe_1[k] += ff_1 * fterm;
                pe_2[k] += ff_2 * fterm;
            }
        }
    }

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_27_2 = 27.0 / 2.0;
    const auto f_44_7 = 44.0 / 7.0;
    const auto f_54_7 = 54.0 / 7.0;
    const auto f_6_7 = 6.0 / 7.0;
    const auto f_8_7 = 8.0 / 7.0;
    const auto fs_11_14_2 = std::sqrt(121.0 / 98.0);
    const auto fs_11_14_30 = std::sqrt(1815.0 / 98.0);
    const auto fs_11_14_42 = std::sqrt(363.0 / 14.0);
    const auto fs_11_14_6 = std::sqrt(363.0 / 98.0);
    const auto fs_11_7_10 = std::sqrt(1210.0 / 49.0);
    const auto fs_11_7_14 = std::sqrt(242.0 / 7.0);
    const auto fs_11_7_15 = std::sqrt(1815.0 / 49.0);
    const auto fs_11_7_6 = std::sqrt(726.0 / 49.0);
    const auto fs_11_7_7 = std::sqrt(121.0 / 7.0);
    const auto fs_18_7_3 = std::sqrt(972.0 / 49.0);
    const auto fs_18_7_5 = std::sqrt(1620.0 / 49.0);
    const auto fs_18_7_6 = std::sqrt(1944.0 / 49.0);
    const auto fs_1_7_2 = std::sqrt(2.0 / 49.0);
    const auto fs_1_7_30 = std::sqrt(30.0 / 49.0);
    const auto fs_1_7_42 = std::sqrt(6.0 / 7.0);
    const auto fs_1_7_6 = std::sqrt(6.0 / 49.0);
    const auto fs_22_7_3 = std::sqrt(1452.0 / 49.0);
    const auto fs_2_7_10 = std::sqrt(40.0 / 49.0);
    const auto fs_2_7_14 = std::sqrt(8.0 / 7.0);
    const auto fs_2_7_15 = std::sqrt(60.0 / 49.0);
    const auto fs_2_7_3 = std::sqrt(12.0 / 49.0);
    const auto fs_2_7_5 = std::sqrt(20.0 / 49.0);
    const auto fs_2_7_6 = std::sqrt(24.0 / 49.0);
    const auto fs_2_7_7 = std::sqrt(4.0 / 7.0);
    const auto fs_36_7_2 = std::sqrt(2592.0 / 49.0);
    const auto fs_4_7_2 = std::sqrt(32.0 / 49.0);
    const auto fs_4_7_3 = std::sqrt(48.0 / 49.0);
    const auto fs_9_1_2 = std::sqrt(162.0);
    const auto fs_9_2_3 = std::sqrt(243.0 / 4.0);
    const auto fs_9_2_5 = std::sqrt(405.0 / 4.0);
    const auto fs_9_2_6 = std::sqrt(243.0 / 2.0);
    const auto fs_9_4_2 = std::sqrt(81.0 / 8.0);
    const auto fs_9_4_30 = std::sqrt(1215.0 / 8.0);
    const auto fs_9_7_2 = std::sqrt(162.0 / 49.0);
    const auto fs_9_7_30 = std::sqrt(2430.0 / 49.0);

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph2_m2, ph2_m1, ph2_0, ph2_p1, ph2_p2, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ab_2, pc_0, pc_1, pc_2, pc_3, pc_4, pc_5, pc_6, pc_7, pc_8, pc_9, pc_10, pc_11, pc_12, pc_13, pc_14, pc_15, pc_16, pc_17, pc_18, pc_19, pc_20 : simd::cache_line_size())
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

        pc_0[k] = e_0 * fs_9_4_30 * h2_p2 + e_1 * fs_11_14_2 * h4_p2 + e_1 * fs_11_7_14 * h4_p4 - e_1 * fs_9_7_30 * r_2 * h2_p2 - e_2 * fs_1_7_2 * r_2 * h4_p2 - e_2 * fs_2_7_14 * r_2 * h4_p4 + e_2 * fs_1_7_30 * r_4 * h2_p2;

        pc_1[k] = e_0 * fs_9_2_5 * h2_p1 + e_1 * fs_11_14_6 * h4_p1 + e_1 * fs_11_14_42 * h4_p3 - e_1 * fs_18_7_5 * r_2 * h2_p1 - e_2 * fs_1_7_6 * r_2 * h4_p1 - e_2 * fs_1_7_42 * r_2 * h4_p3 + e_2 * fs_2_7_5 * r_4 * h2_p1;

        pc_2[k] = e_0 * fs_9_2_6 * h2_0 + e_0 * fs_9_4_2 * h2_p2 + e_1 * fs_11_7_6 * h4_0 + e_1 * fs_11_14_30 * h4_p2 - e_1 * fs_18_7_6 * r_2 * h2_0 - e_1 * fs_9_7_2 * r_2 * h2_p2 - e_2 * fs_2_7_6 * r_2 * h4_0 - e_2 * fs_1_7_30 * r_2 * h4_p2 + e_2 * fs_2_7_6 * r_4 * h2_0 + e_2 * fs_1_7_2 * r_4 * h2_p2;

        pc_3[k] = - e_0 * fs_9_2_3 * h2_m1 - e_1 * fs_11_7_10 * h4_m1 + e_1 * fs_18_7_3 * r_2 * h2_m1 + e_2 * fs_2_7_10 * r_2 * h4_m1 - e_2 * fs_2_7_3 * r_4 * h2_m1;

        pc_4[k] = - e_0 * fs_9_4_2 * h2_m2 - e_1 * fs_11_14_30 * h4_m2 + e_1 * fs_9_7_2 * r_2 * h2_m2 + e_2 * fs_1_7_30 * r_2 * h4_m2 - e_2 * fs_1_7_2 * r_4 * h2_m2;

        pc_5[k] = - e_0 * fs_9_2_5 * h2_m1 - e_1 * fs_11_14_42 * h4_m3 - e_1 * fs_11_14_6 * h4_m1 + e_1 * fs_18_7_5 * r_2 * h2_m1 + e_2 * fs_1_7_42 * r_2 * h4_m3 + e_2 * fs_1_7_6 * r_2 * h4_m1 - e_2 * fs_2_7_5 * r_4 * h2_m1;

        pc_6[k] = - e_0 * fs_9_4_30 * h2_m2 - e_1 * fs_11_7_14 * h4_m4 - e_1 * fs_11_14_2 * h4_m2 + e_1 * fs_9_7_30 * r_2 * h2_m2 + e_2 * fs_2_7_14 * r_2 * h4_m4 + e_2 * fs_1_7_2 * r_2 * h4_m2 - e_2 * fs_1_7_30 * r_4 * h2_m2;

        pc_7[k] = - e_1 * fs_11_7_7 * h4_m3 + e_2 * fs_2_7_7 * r_2 * h4_m3;

        pc_8[k] = e_0 * fs_9_2_5 * h2_m2 - e_1 * fs_22_7_3 * h4_m2 - e_1 * fs_18_7_5 * r_2 * h2_m2 + e_2 * fs_4_7_3 * r_2 * h4_m2 + e_2 * fs_2_7_5 * r_4 * h2_m2;

        pc_9[k] = e_0 * fs_9_1_2 * h2_m1 - e_1 * fs_11_7_15 * h4_m1 - e_1 * fs_36_7_2 * r_2 * h2_m1 + e_2 * fs_2_7_15 * r_2 * h4_m1 + e_2 * fs_4_7_2 * r_4 * h2_m1;

        pc_10[k] = e_0 * f_27_2 * h2_0 - e_1 * f_44_7 * h4_0 - e_1 * f_54_7 * r_2 * h2_0 + e_2 * f_8_7 * r_2 * h4_0 + e_2 * f_6_7 * r_4 * h2_0;

        pc_11[k] = e_0 * fs_9_1_2 * h2_p1 - e_1 * fs_11_7_15 * h4_p1 - e_1 * fs_36_7_2 * r_2 * h2_p1 + e_2 * fs_2_7_15 * r_2 * h4_p1 + e_2 * fs_4_7_2 * r_4 * h2_p1;

        pc_12[k] = e_0 * fs_9_2_5 * h2_p2 - e_1 * fs_22_7_3 * h4_p2 - e_1 * fs_18_7_5 * r_2 * h2_p2 + e_2 * fs_4_7_3 * r_2 * h4_p2 + e_2 * fs_2_7_5 * r_4 * h2_p2;

        pc_13[k] = - e_1 * fs_11_7_7 * h4_p3 + e_2 * fs_2_7_7 * r_2 * h4_p3;

        pc_14[k] = e_0 * fs_9_4_30 * h2_m2 - e_1 * fs_11_7_14 * h4_m4 + e_1 * fs_11_14_2 * h4_m2 - e_1 * fs_9_7_30 * r_2 * h2_m2 + e_2 * fs_2_7_14 * r_2 * h4_m4 - e_2 * fs_1_7_2 * r_2 * h4_m2 + e_2 * fs_1_7_30 * r_4 * h2_m2;

        pc_15[k] = e_0 * fs_9_2_5 * h2_m1 - e_1 * fs_11_14_42 * h4_m3 + e_1 * fs_11_14_6 * h4_m1 - e_1 * fs_18_7_5 * r_2 * h2_m1 + e_2 * fs_1_7_42 * r_2 * h4_m3 - e_2 * fs_1_7_6 * r_2 * h4_m1 + e_2 * fs_2_7_5 * r_4 * h2_m1;

        pc_16[k] = - e_0 * fs_9_4_2 * h2_m2 - e_1 * fs_11_14_30 * h4_m2 + e_1 * fs_9_7_2 * r_2 * h2_m2 + e_2 * fs_1_7_30 * r_2 * h4_m2 - e_2 * fs_1_7_2 * r_4 * h2_m2;

        pc_17[k] = - e_0 * fs_9_2_3 * h2_p1 - e_1 * fs_11_7_10 * h4_p1 + e_1 * fs_18_7_3 * r_2 * h2_p1 + e_2 * fs_2_7_10 * r_2 * h4_p1 - e_2 * fs_2_7_3 * r_4 * h2_p1;

        pc_18[k] = e_0 * fs_9_2_6 * h2_0 - e_0 * fs_9_4_2 * h2_p2 + e_1 * fs_11_7_6 * h4_0 - e_1 * fs_11_14_30 * h4_p2 - e_1 * fs_18_7_6 * r_2 * h2_0 + e_1 * fs_9_7_2 * r_2 * h2_p2 - e_2 * fs_2_7_6 * r_2 * h4_0 + e_2 * fs_1_7_30 * r_2 * h4_p2 + e_2 * fs_2_7_6 * r_4 * h2_0 - e_2 * fs_1_7_2 * r_4 * h2_p2;

        pc_19[k] = e_0 * fs_9_2_5 * h2_p1 + e_1 * fs_11_14_6 * h4_p1 - e_1 * fs_11_14_42 * h4_p3 - e_1 * fs_18_7_5 * r_2 * h2_p1 - e_2 * fs_1_7_6 * r_2 * h4_p1 + e_2 * fs_1_7_42 * r_2 * h4_p3 + e_2 * fs_2_7_5 * r_4 * h2_p1;

        pc_20[k] = e_0 * fs_9_4_30 * h2_p2 + e_1 * fs_11_14_2 * h4_p2 - e_1 * fs_11_7_14 * h4_p4 - e_1 * fs_9_7_30 * r_2 * h2_p2 - e_2 * fs_1_7_2 * r_2 * h4_p2 + e_2 * fs_2_7_14 * r_2 * h4_p4 + e_2 * fs_1_7_30 * r_4 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[21] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

    for (size_t n = 0; n < 21; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
