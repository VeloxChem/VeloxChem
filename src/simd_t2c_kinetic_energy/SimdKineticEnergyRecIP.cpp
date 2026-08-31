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



#include "SimdKineticEnergyRecIP.hpp"

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
compute_ip_kinetic_energy(double                         *values,
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
            false, std::string("SimdKineticEnergyFunc.compute_ip_kinetic_energy: Basis functions must be of angular momenta six and one"));
    }

    if (harmonics.size() < 7)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ip_kinetic_energy: Harmonics must reach angular momentum 7"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ip_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 39 * nvalues, 0.0);

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

            const auto ff_0 = fbase * bexp * bexp * bexp * bexp * bexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * bexp * bexp * bexp * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * bexp * bexp * bexp * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_119_13 = 119.0 / 13.0;
    const auto f_12_13 = 12.0 / 13.0;
    const auto f_14_13 = 14.0 / 13.0;
    const auto f_180_13 = 180.0 / 13.0;
    const auto f_45_1 = 45.0;
    const auto fs_120_13_2 = std::sqrt(28800.0 / 169.0);
    const auto fs_15_13_110 = std::sqrt(24750.0 / 169.0);
    const auto fs_15_13_2 = std::sqrt(450.0 / 169.0);
    const auto fs_15_13_6 = std::sqrt(1350.0 / 169.0);
    const auto fs_15_1_5 = std::sqrt(1125.0);
    const auto fs_15_2_11 = std::sqrt(2475.0 / 4.0);
    const auto fs_15_2_14 = std::sqrt(1575.0 / 2.0);
    const auto fs_15_2_15 = std::sqrt(3375.0 / 4.0);
    const auto fs_15_2_21 = std::sqrt(4725.0 / 4.0);
    const auto fs_15_2_3 = std::sqrt(675.0 / 4.0);
    const auto fs_15_2_33 = std::sqrt(7425.0 / 4.0);
    const auto fs_15_2_35 = std::sqrt(7875.0 / 4.0);
    const auto fs_15_2_5 = std::sqrt(1125.0 / 4.0);
    const auto fs_15_4_110 = std::sqrt(12375.0 / 8.0);
    const auto fs_15_4_2 = std::sqrt(225.0 / 8.0);
    const auto fs_15_4_6 = std::sqrt(675.0 / 8.0);
    const auto fs_17_13_13 = std::sqrt(289.0 / 13.0);
    const auto fs_17_13_21 = std::sqrt(6069.0 / 169.0);
    const auto fs_17_13_3 = std::sqrt(867.0 / 169.0);
    const auto fs_17_13_33 = std::sqrt(9537.0 / 169.0);
    const auto fs_17_13_39 = std::sqrt(867.0 / 13.0);
    const auto fs_17_13_5 = std::sqrt(1445.0 / 169.0);
    const auto fs_17_26_110 = std::sqrt(15895.0 / 338.0);
    const auto fs_17_26_182 = std::sqrt(2023.0 / 26.0);
    const auto fs_17_26_2 = std::sqrt(289.0 / 338.0);
    const auto fs_17_26_30 = std::sqrt(4335.0 / 338.0);
    const auto fs_17_26_6 = std::sqrt(867.0 / 338.0);
    const auto fs_1_13_110 = std::sqrt(110.0 / 169.0);
    const auto fs_1_13_182 = std::sqrt(14.0 / 13.0);
    const auto fs_1_13_2 = std::sqrt(2.0 / 169.0);
    const auto fs_1_13_30 = std::sqrt(30.0 / 169.0);
    const auto fs_1_13_6 = std::sqrt(6.0 / 169.0);
    const auto fs_2_13_11 = std::sqrt(44.0 / 169.0);
    const auto fs_2_13_13 = std::sqrt(4.0 / 13.0);
    const auto fs_2_13_14 = std::sqrt(56.0 / 169.0);
    const auto fs_2_13_15 = std::sqrt(60.0 / 169.0);
    const auto fs_2_13_21 = std::sqrt(84.0 / 169.0);
    const auto fs_2_13_3 = std::sqrt(12.0 / 169.0);
    const auto fs_2_13_33 = std::sqrt(132.0 / 169.0);
    const auto fs_2_13_35 = std::sqrt(140.0 / 169.0);
    const auto fs_2_13_39 = std::sqrt(12.0 / 13.0);
    const auto fs_2_13_5 = std::sqrt(20.0 / 169.0);
    const auto fs_30_13_11 = std::sqrt(9900.0 / 169.0);
    const auto fs_30_13_14 = std::sqrt(12600.0 / 169.0);
    const auto fs_30_13_15 = std::sqrt(13500.0 / 169.0);
    const auto fs_30_13_21 = std::sqrt(18900.0 / 169.0);
    const auto fs_30_13_3 = std::sqrt(2700.0 / 169.0);
    const auto fs_30_13_33 = std::sqrt(29700.0 / 169.0);
    const auto fs_30_13_35 = std::sqrt(31500.0 / 169.0);
    const auto fs_30_13_5 = std::sqrt(4500.0 / 169.0);
    const auto fs_30_1_2 = std::sqrt(1800.0);
    const auto fs_34_13_10 = std::sqrt(11560.0 / 169.0);
    const auto fs_34_13_6 = std::sqrt(6936.0 / 169.0);
    const auto fs_34_13_7 = std::sqrt(8092.0 / 169.0);
    const auto fs_3_13_10 = std::sqrt(90.0 / 169.0);
    const auto fs_45_13_10 = std::sqrt(20250.0 / 169.0);
    const auto fs_45_2_2 = std::sqrt(2025.0 / 2.0);
    const auto fs_45_2_3 = std::sqrt(6075.0 / 4.0);
    const auto fs_45_4_10 = std::sqrt(10125.0 / 8.0);
    const auto fs_4_13_10 = std::sqrt(160.0 / 169.0);
    const auto fs_4_13_5 = std::sqrt(80.0 / 169.0);
    const auto fs_4_13_6 = std::sqrt(96.0 / 169.0);
    const auto fs_4_13_7 = std::sqrt(112.0 / 169.0);
    const auto fs_51_13_2 = std::sqrt(5202.0 / 169.0);
    const auto fs_51_13_5 = std::sqrt(13005.0 / 169.0);
    const auto fs_51_26_10 = std::sqrt(13005.0 / 338.0);
    const auto fs_60_13_5 = std::sqrt(18000.0 / 169.0);
    const auto fs_68_13_3 = std::sqrt(13872.0 / 169.0);
    const auto fs_6_13_2 = std::sqrt(72.0 / 169.0);
    const auto fs_6_13_3 = std::sqrt(108.0 / 169.0);
    const auto fs_6_13_5 = std::sqrt(180.0 / 169.0);
    const auto fs_8_13_2 = std::sqrt(128.0 / 169.0);
    const auto fs_8_13_3 = std::sqrt(192.0 / 169.0);
    const auto fs_90_13_2 = std::sqrt(16200.0 / 169.0);
    const auto fs_90_13_3 = std::sqrt(24300.0 / 169.0);

    // NOTE: the angular components are formed in 5 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_m5, ph5_m4, ph5_p3, ph5_p4, ph5_p5, ph7_m7, ph7_m6, ph7_m5, ph7_m4, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ph7_p7, ab_2, pc_0, pc_1, pc_2, pc_3, pc_4, pc_5, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_0[k] = - e_0 * fs_15_2_33 * h5_p5 - e_1 * fs_17_26_2 * h7_p5 - e_1 * fs_17_26_182 * h7_p7 + e_1 * fs_30_13_33 * r_2 * h5_p5 + e_2 * fs_1_13_2 * r_2 * h7_p5 + e_2 * fs_1_13_182 * r_2 * h7_p7 - e_2 * fs_2_13_33 * r_4 * h5_p5;

        pc_1[k] = e_1 * fs_17_13_13 * h7_m6 - e_2 * fs_2_13_13 * r_2 * h7_m6;

        pc_2[k] = - e_0 * fs_15_2_33 * h5_m5 + e_1 * fs_17_26_182 * h7_m7 - e_1 * fs_17_26_2 * h7_m5 + e_1 * fs_30_13_33 * r_2 * h5_m5 - e_2 * fs_1_13_182 * r_2 * h7_m7 + e_2 * fs_1_13_2 * r_2 * h7_m5 - e_2 * fs_2_13_33 * r_4 * h5_m5;

        pc_3[k] = - e_0 * fs_15_4_110 * h5_p4 - e_1 * fs_17_26_6 * h7_p4 - e_1 * fs_17_13_39 * h7_p6 + e_1 * fs_15_13_110 * r_2 * h5_p4 + e_2 * fs_1_13_6 * r_2 * h7_p4 + e_2 * fs_2_13_39 * r_2 * h7_p6 - e_2 * fs_1_13_110 * r_4 * h5_p4;

        pc_4[k] = - e_0 * fs_15_2_11 * h5_m5 + e_1 * fs_34_13_6 * h7_m5 + e_1 * fs_30_13_11 * r_2 * h5_m5 - e_2 * fs_4_13_6 * r_2 * h7_m5 - e_2 * fs_2_13_11 * r_4 * h5_m5;

        pc_5[k] = - e_0 * fs_15_4_110 * h5_m4 + e_1 * fs_17_13_39 * h7_m6 - e_1 * fs_17_26_6 * h7_m4 + e_1 * fs_15_13_110 * r_2 * h5_m4 - e_2 * fs_2_13_39 * r_2 * h7_m6 + e_2 * fs_1_13_6 * r_2 * h7_m4 - e_2 * fs_1_13_110 * r_4 * h5_m4;

        pc_6[k] = - e_0 * fs_45_4_10 * h5_p3 - e_0 * fs_15_4_2 * h5_p5 - e_1 * fs_17_13_3 * h7_p3 - e_1 * fs_17_13_33 * h7_p5 + e_1 * fs_45_13_10 * r_2 * h5_p3 + e_1 * fs_15_13_2 * r_2 * h5_p5 + e_2 * fs_2_13_3 * r_2 * h7_p3 + e_2 * fs_2_13_33 * r_2 * h7_p5 - e_2 * fs_3_13_10 * r_4 * h5_p3 - e_2 * fs_1_13_2 * r_4 * h5_p5;

        pc_7[k] = - e_0 * fs_15_1_5 * h5_m4 + e_1 * fs_17_13_33 * h7_m4 + e_1 * fs_60_13_5 * r_2 * h5_m4 - e_2 * fs_2_13_33 * r_2 * h7_m4 - e_2 * fs_4_13_5 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 5 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph5_p2, ph5_p4, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ph7_p2, ph7_p4, ab_2, pc_8, pc_9, pc_10, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_8[k] = e_0 * fs_15_4_2 * h5_m5 - e_0 * fs_45_4_10 * h5_m3 + e_1 * fs_17_13_33 * h7_m5 - e_1 * fs_17_13_3 * h7_m3 - e_1 * fs_15_13_2 * r_2 * h5_m5 + e_1 * fs_45_13_10 * r_2 * h5_m3 - e_2 * fs_2_13_33 * r_2 * h7_m5 + e_2 * fs_2_13_3 * r_2 * h7_m3 + e_2 * fs_1_13_2 * r_4 * h5_m5 - e_2 * fs_3_13_10 * r_4 * h5_m3;

        pc_9[k] = - e_0 * fs_45_2_2 * h5_p2 - e_0 * fs_15_4_6 * h5_p4 - e_1 * fs_17_13_5 * h7_p2 - e_1 * fs_17_26_110 * h7_p4 + e_1 * fs_90_13_2 * r_2 * h5_p2 + e_1 * fs_15_13_6 * r_2 * h5_p4 + e_2 * fs_2_13_5 * r_2 * h7_p2 + e_2 * fs_1_13_110 * r_2 * h7_p4 - e_2 * fs_6_13_2 * r_4 * h5_p2 - e_2 * fs_1_13_6 * r_4 * h5_p4;

        pc_10[k] = - e_0 * fs_45_2_3 * h5_m3 + e_1 * fs_34_13_10 * h7_m3 + e_1 * fs_90_13_3 * r_2 * h5_m3 - e_2 * fs_4_13_10 * r_2 * h7_m3 - e_2 * fs_6_13_3 * r_4 * h5_m3;

        pc_11[k] = e_0 * fs_15_4_6 * h5_m4 - e_0 * fs_45_2_2 * h5_m2 + e_1 * fs_17_26_110 * h7_m4 - e_1 * fs_17_13_5 * h7_m2 - e_1 * fs_15_13_6 * r_2 * h5_m4 + e_1 * fs_90_13_2 * r_2 * h5_m2 - e_2 * fs_1_13_110 * r_2 * h7_m4 + e_2 * fs_2_13_5 * r_2 * h7_m2 + e_2 * fs_1_13_6 * r_4 * h5_m4 - e_2 * fs_6_13_2 * r_4 * h5_m2;
    }

    // NOTE: the angular components are formed in 5 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_m3, ph5_m2, ph5_m1, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph7_m3, ph7_m2, ph7_m1, ph7_0, ph7_p1, ph7_p2, ph7_p3, ab_2, pc_12, pc_13, pc_14, pc_15, pc_16, pc_17, pc_18, pc_19, pc_20, pc_21, pc_22, pc_23, pc_24, pc_25, pc_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_12[k] = - e_0 * fs_15_2_14 * h5_p1 - e_0 * fs_15_2_3 * h5_p3 - e_1 * fs_17_26_30 * h7_p1 - e_1 * fs_51_26_10 * h7_p3 + e_1 * fs_30_13_14 * r_2 * h5_p1 + e_1 * fs_30_13_3 * r_2 * h5_p3 + e_2 * fs_1_13_30 * r_2 * h7_p1 + e_2 * fs_3_13_10 * r_2 * h7_p3 - e_2 * fs_2_13_14 * r_4 * h5_p1 - e_2 * fs_2_13_3 * r_4 * h5_p3;

        pc_13[k] = - e_0 * fs_30_1_2 * h5_m2 + e_1 * fs_51_13_5 * h7_m2 + e_1 * fs_120_13_2 * r_2 * h5_m2 - e_2 * fs_6_13_5 * r_2 * h7_m2 - e_2 * fs_8_13_2 * r_4 * h5_m2;

        pc_14[k] = e_0 * fs_15_2_3 * h5_m3 - e_0 * fs_15_2_14 * h5_m1 + e_1 * fs_51_26_10 * h7_m3 - e_1 * fs_17_26_30 * h7_m1 - e_1 * fs_30_13_3 * r_2 * h5_m3 + e_1 * fs_30_13_14 * r_2 * h5_m1 - e_2 * fs_3_13_10 * r_2 * h7_m3 + e_2 * fs_1_13_30 * r_2 * h7_m1 + e_2 * fs_2_13_3 * r_4 * h5_m3 - e_2 * fs_2_13_14 * r_4 * h5_m1;

        pc_15[k] = - e_0 * fs_15_2_21 * h5_0 - e_0 * fs_15_2_5 * h5_p2 - e_1 * fs_17_13_21 * h7_0 - e_1 * fs_51_13_2 * h7_p2 + e_1 * fs_30_13_21 * r_2 * h5_0 + e_1 * fs_30_13_5 * r_2 * h5_p2 + e_2 * fs_2_13_21 * r_2 * h7_0 + e_2 * fs_6_13_2 * r_2 * h7_p2 - e_2 * fs_2_13_21 * r_4 * h5_0 - e_2 * fs_2_13_5 * r_4 * h5_p2;

        pc_16[k] = - e_0 * fs_15_2_35 * h5_m1 + e_1 * fs_68_13_3 * h7_m1 + e_1 * fs_30_13_35 * r_2 * h5_m1 - e_2 * fs_8_13_3 * r_2 * h7_m1 - e_2 * fs_2_13_35 * r_4 * h5_m1;

        pc_17[k] = e_0 * fs_15_2_5 * h5_m2 + e_1 * fs_51_13_2 * h7_m2 - e_1 * fs_30_13_5 * r_2 * h5_m2 - e_2 * fs_6_13_2 * r_2 * h7_m2 + e_2 * fs_2_13_5 * r_4 * h5_m2;

        pc_18[k] = e_0 * fs_15_2_15 * h5_m1 + e_1 * fs_34_13_7 * h7_m1 - e_1 * fs_30_13_15 * r_2 * h5_m1 - e_2 * fs_4_13_7 * r_2 * h7_m1 + e_2 * fs_2_13_15 * r_4 * h5_m1;

        pc_19[k] = - e_0 * f_45_1 * h5_0 + e_1 * f_119_13 * h7_0 + e_1 * f_180_13 * r_2 * h5_0 - e_2 * f_14_13 * r_2 * h7_0 - e_2 * f_12_13 * r_4 * h5_0;

        pc_20[k] = e_0 * fs_15_2_15 * h5_p1 + e_1 * fs_34_13_7 * h7_p1 - e_1 * fs_30_13_15 * r_2 * h5_p1 - e_2 * fs_4_13_7 * r_2 * h7_p1 + e_2 * fs_2_13_15 * r_4 * h5_p1;

        pc_21[k] = e_0 * fs_15_2_5 * h5_m2 + e_1 * fs_51_13_2 * h7_m2 - e_1 * fs_30_13_5 * r_2 * h5_m2 - e_2 * fs_6_13_2 * r_2 * h7_m2 + e_2 * fs_2_13_5 * r_4 * h5_m2;

        pc_22[k] = - e_0 * fs_15_2_35 * h5_p1 + e_1 * fs_68_13_3 * h7_p1 + e_1 * fs_30_13_35 * r_2 * h5_p1 - e_2 * fs_8_13_3 * r_2 * h7_p1 - e_2 * fs_2_13_35 * r_4 * h5_p1;

        pc_23[k] = - e_0 * fs_15_2_21 * h5_0 + e_0 * fs_15_2_5 * h5_p2 - e_1 * fs_17_13_21 * h7_0 + e_1 * fs_51_13_2 * h7_p2 + e_1 * fs_30_13_21 * r_2 * h5_0 - e_1 * fs_30_13_5 * r_2 * h5_p2 + e_2 * fs_2_13_21 * r_2 * h7_0 - e_2 * fs_6_13_2 * r_2 * h7_p2 - e_2 * fs_2_13_21 * r_4 * h5_0 + e_2 * fs_2_13_5 * r_4 * h5_p2;

        pc_24[k] = e_0 * fs_15_2_3 * h5_m3 + e_0 * fs_15_2_14 * h5_m1 + e_1 * fs_51_26_10 * h7_m3 + e_1 * fs_17_26_30 * h7_m1 - e_1 * fs_30_13_3 * r_2 * h5_m3 - e_1 * fs_30_13_14 * r_2 * h5_m1 - e_2 * fs_3_13_10 * r_2 * h7_m3 - e_2 * fs_1_13_30 * r_2 * h7_m1 + e_2 * fs_2_13_3 * r_4 * h5_m3 + e_2 * fs_2_13_14 * r_4 * h5_m1;

        pc_25[k] = - e_0 * fs_30_1_2 * h5_p2 + e_1 * fs_51_13_5 * h7_p2 + e_1 * fs_120_13_2 * r_2 * h5_p2 - e_2 * fs_6_13_5 * r_2 * h7_p2 - e_2 * fs_8_13_2 * r_4 * h5_p2;

        pc_26[k] = - e_0 * fs_15_2_14 * h5_p1 + e_0 * fs_15_2_3 * h5_p3 - e_1 * fs_17_26_30 * h7_p1 + e_1 * fs_51_26_10 * h7_p3 + e_1 * fs_30_13_14 * r_2 * h5_p1 - e_1 * fs_30_13_3 * r_2 * h5_p3 + e_2 * fs_1_13_30 * r_2 * h7_p1 - e_2 * fs_3_13_10 * r_2 * h7_p3 - e_2 * fs_2_13_14 * r_4 * h5_p1 + e_2 * fs_2_13_3 * r_4 * h5_p3;
    }

    // NOTE: the angular components are formed in 5 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph5_p2, ph5_p3, ph5_p4, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ph7_p2, ph7_p3, ph7_p4, ab_2, pc_27, pc_28, pc_29, pc_30, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_27[k] = e_0 * fs_15_4_6 * h5_m4 + e_0 * fs_45_2_2 * h5_m2 + e_1 * fs_17_26_110 * h7_m4 + e_1 * fs_17_13_5 * h7_m2 - e_1 * fs_15_13_6 * r_2 * h5_m4 - e_1 * fs_90_13_2 * r_2 * h5_m2 - e_2 * fs_1_13_110 * r_2 * h7_m4 - e_2 * fs_2_13_5 * r_2 * h7_m2 + e_2 * fs_1_13_6 * r_4 * h5_m4 + e_2 * fs_6_13_2 * r_4 * h5_m2;

        pc_28[k] = - e_0 * fs_45_2_3 * h5_p3 + e_1 * fs_34_13_10 * h7_p3 + e_1 * fs_90_13_3 * r_2 * h5_p3 - e_2 * fs_4_13_10 * r_2 * h7_p3 - e_2 * fs_6_13_3 * r_4 * h5_p3;

        pc_29[k] = - e_0 * fs_45_2_2 * h5_p2 + e_0 * fs_15_4_6 * h5_p4 - e_1 * fs_17_13_5 * h7_p2 + e_1 * fs_17_26_110 * h7_p4 + e_1 * fs_90_13_2 * r_2 * h5_p2 - e_1 * fs_15_13_6 * r_2 * h5_p4 + e_2 * fs_2_13_5 * r_2 * h7_p2 - e_2 * fs_1_13_110 * r_2 * h7_p4 - e_2 * fs_6_13_2 * r_4 * h5_p2 + e_2 * fs_1_13_6 * r_4 * h5_p4;

        pc_30[k] = e_0 * fs_15_4_2 * h5_m5 + e_0 * fs_45_4_10 * h5_m3 + e_1 * fs_17_13_33 * h7_m5 + e_1 * fs_17_13_3 * h7_m3 - e_1 * fs_15_13_2 * r_2 * h5_m5 - e_1 * fs_45_13_10 * r_2 * h5_m3 - e_2 * fs_2_13_33 * r_2 * h7_m5 - e_2 * fs_2_13_3 * r_2 * h7_m3 + e_2 * fs_1_13_2 * r_4 * h5_m5 + e_2 * fs_3_13_10 * r_4 * h5_m3;

        pc_31[k] = - e_0 * fs_15_1_5 * h5_p4 + e_1 * fs_17_13_33 * h7_p4 + e_1 * fs_60_13_5 * r_2 * h5_p4 - e_2 * fs_2_13_33 * r_2 * h7_p4 - e_2 * fs_4_13_5 * r_4 * h5_p4;
    }

    // NOTE: the angular components are formed in 5 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph5_m5, ph5_m4, ph5_p3, ph5_p4, ph5_p5, ph7_m7, ph7_m6, ph7_m5, ph7_m4, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ph7_p7, ab_2, pc_32, pc_33, pc_34, pc_35, pc_36, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_32[k] = - e_0 * fs_45_4_10 * h5_p3 + e_0 * fs_15_4_2 * h5_p5 - e_1 * fs_17_13_3 * h7_p3 + e_1 * fs_17_13_33 * h7_p5 + e_1 * fs_45_13_10 * r_2 * h5_p3 - e_1 * fs_15_13_2 * r_2 * h5_p5 + e_2 * fs_2_13_3 * r_2 * h7_p3 - e_2 * fs_2_13_33 * r_2 * h7_p5 - e_2 * fs_3_13_10 * r_4 * h5_p3 + e_2 * fs_1_13_2 * r_4 * h5_p5;

        pc_33[k] = e_0 * fs_15_4_110 * h5_m4 + e_1 * fs_17_13_39 * h7_m6 + e_1 * fs_17_26_6 * h7_m4 - e_1 * fs_15_13_110 * r_2 * h5_m4 - e_2 * fs_2_13_39 * r_2 * h7_m6 - e_2 * fs_1_13_6 * r_2 * h7_m4 + e_2 * fs_1_13_110 * r_4 * h5_m4;

        pc_34[k] = - e_0 * fs_15_2_11 * h5_p5 + e_1 * fs_34_13_6 * h7_p5 + e_1 * fs_30_13_11 * r_2 * h5_p5 - e_2 * fs_4_13_6 * r_2 * h7_p5 - e_2 * fs_2_13_11 * r_4 * h5_p5;

        pc_35[k] = - e_0 * fs_15_4_110 * h5_p4 - e_1 * fs_17_26_6 * h7_p4 + e_1 * fs_17_13_39 * h7_p6 + e_1 * fs_15_13_110 * r_2 * h5_p4 + e_2 * fs_1_13_6 * r_2 * h7_p4 - e_2 * fs_2_13_39 * r_2 * h7_p6 - e_2 * fs_1_13_110 * r_4 * h5_p4;

        pc_36[k] = e_0 * fs_15_2_33 * h5_m5 + e_1 * fs_17_26_182 * h7_m7 + e_1 * fs_17_26_2 * h7_m5 - e_1 * fs_30_13_33 * r_2 * h5_m5 - e_2 * fs_1_13_182 * r_2 * h7_m7 - e_2 * fs_1_13_2 * r_2 * h7_m5 + e_2 * fs_2_13_33 * r_4 * h5_m5;

        pc_37[k] = e_1 * fs_17_13_13 * h7_p6 - e_2 * fs_2_13_13 * r_2 * h7_p6;

        pc_38[k] = - e_0 * fs_15_2_33 * h5_p5 - e_1 * fs_17_26_2 * h7_p5 + e_1 * fs_17_26_182 * h7_p7 + e_1 * fs_30_13_33 * r_2 * h5_p5 + e_2 * fs_1_13_2 * r_2 * h7_p5 - e_2 * fs_1_13_182 * r_2 * h7_p7 - e_2 * fs_2_13_33 * r_4 * h5_p5;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[39] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};

    for (size_t n = 0; n < 39; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
