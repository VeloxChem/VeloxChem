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



#include "SimdThreeCenterElectronRepulsionRecGSS.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdBoysFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdVariableMatrix.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gss_electron_repulsion(double                         *values,
                               const size_t                    npairs,
                               const size_t                    natoms,
                               const size_t                    iatom,
                               const CBasisFunction           &a_function,
                               const CBasisFunction           &b_function,
                               const CBasisFunction           &c_function,
                               const std::vector<CSimdMatrix> &ab_harmonics,
                               const std::vector<CSimdMatrix> &bc_harmonics,
                               const CSimdMatrix              &ab_coordinates,
                               const CSimdMatrix              &bc_coordinates,
                               const double                    threshold) -> void
{
    if ((a_function.get_angular_momentum() != 4) || (b_function.get_angular_momentum() != 0) ||
        (c_function.get_angular_momentum() != 0))
    {
        errors::assertMsgCritical(
            false,
            std::string("SimdThreeCenterElectronRepulsionRecGSS.compute_gss_electron_repulsion: Basis functions must be of angular momenta four, zero and zero"));
    }

    if ((ab_harmonics.size() < 4) || (bc_harmonics.size() < 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdThreeCenterElectronRepulsionRecGSS.compute_gss_electron_repulsion: Harmonics must reach angular momentum four"));
    }

    if (npairs > ab_coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdThreeCenterElectronRepulsionRecGSS.compute_gss_electron_repulsion: Number of atom pairs exceeds coordinates"));
    }

    if (iatom >= natoms)
    {
        errors::assertMsgCritical(
            false, std::string("SimdThreeCenterElectronRepulsionRecGSS.compute_gss_electron_repulsion: Index of atom on c side is out of range"));
    }

    if (npairs == 0) return;

    const auto &a_exps = a_function.exponents();

    const auto &b_exps = b_function.exponents();

    const auto &c_exps = c_function.exponents();

    const auto &a_norms = a_function.normalization_factors();

    const auto &b_norms = b_function.normalization_factors();

    const auto &c_norms = c_function.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprim_c = c_exps.size();

    const auto nprims = nprim_a * nprim_b * nprim_c;

    // NOTE: the triples of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(a_function,
                                                             b_function,
                                                             c_function,
                                                             npairs,
                                                             ab_coordinates,
                                                             screenfunc::three_center_electron_repulsion_primitive_bound,
                                                             threshold / static_cast<double>(nprims));

    const auto nmax = *std::ranges::max_element(dimensions);

    // NOTE: the values of one atom on c side are contiguous over the atom pairs,
    // the atoms are npairs apart and the angular components of the atoms on c
    // side are npairs times the number of those atoms apart, as the sparsity
    // pattern lays them out.

    const auto stride = natoms * npairs;

    double *slices[9];

    for (size_t m = 0; m < 9; m++) slices[m] = values + m * stride + iatom * npairs;

    if (nmax == 0)
    {
        for (size_t m = 0; m < 9; m++) std::fill(slices[m], slices[m] + npairs, 0.0);

        return;
    }

    const auto *a_x = ab_coordinates.data(0);
    const auto *a_y = ab_coordinates.data(1);
    const auto *a_z = ab_coordinates.data(2);

    const auto *ab_2 = ab_coordinates.data(6);

    const auto *b_x = bc_coordinates.data(0);
    const auto *b_y = bc_coordinates.data(1);
    const auto *b_z = bc_coordinates.data(2);

    const auto *c_x = bc_coordinates.data(3);
    const auto *c_y = bc_coordinates.data(4);
    const auto *c_z = bc_coordinates.data(5);

    auto factors = CSimdMatrix(2, nmax);

    auto *e_ab = factors.data(0);

    auto *pc_2 = factors.data(1);

    // NOTE: one row accumulates for each bidegree of the addition theorem, as
    // they carry different powers of the exponents and cannot share an
    // accumulator. The row of index l1 multiplies the harmonics of degree l1 of
    // the atom pairs and of degree 4 less l1 of the atoms on c side.

    auto buffer = CSimdMatrix(5, nmax);

    buffer.zero();

    auto *acc_0 = buffer.data(0);

    auto *acc_1 = buffer.data(1);

    auto *acc_2 = buffer.data(2);

    auto *acc_3 = buffer.data(3);

    auto *acc_4 = buffer.data(4);

    constexpr auto fpi = mathconst::pi_value();

    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto bexp = b_exps[j];

            // NOTE: the widest of the triples of this pair of primitives is
            // searched for rather than assumed to be the last, as the bound of a
            // triple carries its prefactor as well as its decay.

            const auto first = dimensions.begin() + static_cast<long>((i * nprim_b + j) * nprim_c);

            const auto npair_max = *std::ranges::max_element(first, first + static_cast<long>(nprim_c));

            if (npair_max == 0) continue;

            const auto pexp = aexp + bexp;

            const auto fmu = aexp * bexp / pexp;

            const auto frp = 1.0 / pexp;

#pragma omp simd aligned(e_ab, pc_2, ab_2, a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z : simd::cache_line_size())
            for (size_t k = 0; k < npair_max; k++)
            {
                e_ab[k] = std::exp(-fmu * ab_2[k]);

                const auto p_x = frp * (aexp * (a_x[k] - c_x[k]) + bexp * (b_x[k] - c_x[k]));

                const auto p_y = frp * (aexp * (a_y[k] - c_y[k]) + bexp * (b_y[k] - c_y[k]));

                const auto p_z = frp * (aexp * (a_z[k] - c_z[k]) + bexp * (b_z[k] - c_z[k]));

                pc_2[k] = p_x * p_x + p_y * p_y + p_z * p_z;
            }

            // NOTE: the Boys function of every primitive on c side of this pair
            // is computed by one call, which fills the orders zero to four of
            // every row. The integrals need the order four alone, and the lower
            // orders are formed on the way to it by the recursion.

            auto boys = CSimdVariableMatrix(std::vector<size_t>(first, first + static_cast<long>(nprim_c)), 6);

            for (size_t k = 0; k < nprim_c; k++)
            {
                const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                if (ncols == 0) continue;

                const auto frho = pexp * c_exps[k] / (pexp + c_exps[k]);

                auto *bargs = boys.data(0, k);

#pragma omp simd aligned(bargs, pc_2 : simd::cache_line_size())
                for (size_t l = 0; l < ncols; l++)
                {
                    bargs[l] = frho * pc_2[l];
                }
            }

            simdfunc::compute_boys_function(boys);

            for (size_t k = 0; k < nprim_c; k++)
            {
                const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                if (ncols == 0) continue;

                const auto cexp = c_exps[k];

                const auto qexp = pexp + cexp;

                const auto fbase = fcoul * anorm * b_norms[j] * c_norms[k] / (pexp * cexp * std::sqrt(qexp));

                const auto frq = 1.0 / qexp;

                const auto frat = bexp * frp;

                const auto fag = aexp * cexp * frp * frq;

                const auto fgq = cexp * frq;

                const auto w_0_0 = fgq * fgq * fgq * fgq;

                const auto w_1_0 = frat * fgq * fgq * fgq;

                const auto w_1_1 = fag * fgq * fgq * fgq;

                const auto w_2_0 = frat * frat * fgq * fgq;

                const auto w_2_1 = 2.0 * frat * fag * fgq * fgq;

                const auto w_2_2 = fag * fag * fgq * fgq;

                const auto w_3_0 = frat * frat * frat * fgq;

                const auto w_3_1 = 3.0 * frat * frat * fag * fgq;

                const auto w_3_2 = 3.0 * frat * fag * fag * fgq;

                const auto w_3_3 = fag * fag * fag * fgq;

                const auto w_4_0 = frat * frat * frat * frat;

                const auto w_4_1 = 4.0 * frat * frat * frat * fag;

                const auto w_4_2 = 6.0 * frat * frat * fag * fag;

                const auto w_4_3 = 4.0 * frat * fag * fag * fag;

                const auto w_4_4 = fag * fag * fag * fag;

                const auto *bv_0 = boys.data(1, k);

                const auto *bv_1 = boys.data(2, k);

                const auto *bv_2 = boys.data(3, k);

                const auto *bv_3 = boys.data(4, k);

                const auto *bv_4 = boys.data(5, k);

#pragma omp simd aligned(acc_0, acc_1, acc_2, acc_3, acc_4, e_ab, bv_0, bv_1, bv_2, bv_3, bv_4 : simd::cache_line_size())
                for (size_t l = 0; l < ncols; l++)
                {
                    const auto fval = fbase * e_ab[l];

                    acc_0[l] += fval * (w_0_0 * bv_4[l]);

                    acc_1[l] += fval * (w_1_0 * bv_3[l] + w_1_1 * bv_4[l]);

                    acc_2[l] += fval * (w_2_0 * bv_2[l] + w_2_1 * bv_3[l] + w_2_2 * bv_4[l]);

                    acc_3[l] += fval * (w_3_0 * bv_1[l] + w_3_1 * bv_2[l] + w_3_2 * bv_3[l] + w_3_3 * bv_4[l]);

                    acc_4[l] += fval * (w_4_0 * bv_0[l] + w_4_1 * bv_1[l] + w_4_2 * bv_2[l] + w_4_3 * bv_3[l] + w_4_4 * bv_4[l]);
                }
            }
        }
    }

    // NOTE: the bidegrees are accumulated into the angular components one at a
    // time, so that no loop holds the harmonics of every degree at once and the
    // vectorizer keeps its registers.

    auto components = CSimdMatrix(9, nmax);

    components.zero();

    auto *out_m4 = components.data(0);
    auto *out_m3 = components.data(1);
    auto *out_m2 = components.data(2);
    auto *out_m1 = components.data(3);
    auto *out_0 = components.data(4);
    auto *out_p1 = components.data(5);
    auto *out_p2 = components.data(6);
    auto *out_p3 = components.data(7);
    auto *out_p4 = components.data(8);

    // the bidegree of degree zero on the atom pairs and four on the atoms
    // on c side, whose coefficients are one: the harmonic of the other side is
    // of degree zero and is one for every atom pair

    {
        const auto *h_m4 = bc_harmonics[3].data(0);
        const auto *h_m3 = bc_harmonics[3].data(1);
        const auto *h_m2 = bc_harmonics[3].data(2);
        const auto *h_m1 = bc_harmonics[3].data(3);
        const auto *h_0 = bc_harmonics[3].data(4);
        const auto *h_p1 = bc_harmonics[3].data(5);
        const auto *h_p2 = bc_harmonics[3].data(6);
        const auto *h_p3 = bc_harmonics[3].data(7);
        const auto *h_p4 = bc_harmonics[3].data(8);

#pragma omp simd aligned(out_m4, out_m3, out_m2, out_m1, out_0, out_p1, out_p2, out_p3, out_p4, acc_0, h_m4, h_m3, h_m2, h_m1, h_0, h_p1, h_p2, h_p3, h_p4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            const auto f = acc_0[k];

            out_m4[k] += f * h_m4[k];
            out_m3[k] += f * h_m3[k];
            out_m2[k] += f * h_m2[k];
            out_m1[k] += f * h_m1[k];
            out_0[k] += f * h_0[k];
            out_p1[k] += f * h_p1[k];
            out_p2[k] += f * h_p2[k];
            out_p3[k] += f * h_p3[k];
            out_p4[k] += f * h_p4[k];
        }
    }

    // the bidegree of degree one on the atom pairs and three on the atoms
    // on c side, whose 21 products of harmonics are formed once and read by
    // the angular components which carry them

    {
        auto products = CSimdMatrix(21, nmax);

        const auto *p_m1 = ab_harmonics[0].data(0);
        const auto *p_0 = ab_harmonics[0].data(1);
        const auto *p_p1 = ab_harmonics[0].data(2);

        const auto *q_m3 = bc_harmonics[2].data(0);
        const auto *q_m2 = bc_harmonics[2].data(1);
        const auto *q_m1 = bc_harmonics[2].data(2);
        const auto *q_0 = bc_harmonics[2].data(3);
        const auto *q_p1 = bc_harmonics[2].data(4);
        const auto *q_p2 = bc_harmonics[2].data(5);
        const auto *q_p3 = bc_harmonics[2].data(6);

        auto *r_m1_m3 = products.data(0);
        auto *r_m1_m2 = products.data(1);
        auto *r_m1_m1 = products.data(2);
        auto *r_m1_0 = products.data(3);
        auto *r_m1_p1 = products.data(4);
        auto *r_m1_p2 = products.data(5);
        auto *r_m1_p3 = products.data(6);
        auto *r_0_m3 = products.data(7);
        auto *r_0_m2 = products.data(8);
        auto *r_0_m1 = products.data(9);
        auto *r_0_0 = products.data(10);
        auto *r_0_p1 = products.data(11);
        auto *r_0_p2 = products.data(12);
        auto *r_0_p3 = products.data(13);
        auto *r_p1_m3 = products.data(14);
        auto *r_p1_m2 = products.data(15);
        auto *r_p1_m1 = products.data(16);
        auto *r_p1_0 = products.data(17);
        auto *r_p1_p1 = products.data(18);
        auto *r_p1_p2 = products.data(19);
        auto *r_p1_p3 = products.data(20);

#pragma omp simd aligned(r_m1_m3, r_m1_m2, r_m1_m1, r_m1_0, r_m1_p1, r_m1_p2, r_m1_p3, r_0_m3, r_0_m2, r_0_m1, r_0_0, r_0_p1, r_0_p2, r_0_p3, r_p1_m3, r_p1_m2, r_p1_m1, r_p1_0, r_p1_p1, r_p1_p2, r_p1_p3, p_m1, p_0, p_p1, q_m3, q_m2, q_m1, q_0, q_p1, q_p2, q_p3 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            r_m1_m3[k] = p_m1[k] * q_m3[k];
            r_m1_m2[k] = p_m1[k] * q_m2[k];
            r_m1_m1[k] = p_m1[k] * q_m1[k];
            r_m1_0[k] = p_m1[k] * q_0[k];
            r_m1_p1[k] = p_m1[k] * q_p1[k];
            r_m1_p2[k] = p_m1[k] * q_p2[k];
            r_m1_p3[k] = p_m1[k] * q_p3[k];
            r_0_m3[k] = p_0[k] * q_m3[k];
            r_0_m2[k] = p_0[k] * q_m2[k];
            r_0_m1[k] = p_0[k] * q_m1[k];
            r_0_0[k] = p_0[k] * q_0[k];
            r_0_p1[k] = p_0[k] * q_p1[k];
            r_0_p2[k] = p_0[k] * q_p2[k];
            r_0_p3[k] = p_0[k] * q_p3[k];
            r_p1_m3[k] = p_p1[k] * q_m3[k];
            r_p1_m2[k] = p_p1[k] * q_m2[k];
            r_p1_m1[k] = p_p1[k] * q_m1[k];
            r_p1_0[k] = p_p1[k] * q_0[k];
            r_p1_p1[k] = p_p1[k] * q_p1[k];
            r_p1_p2[k] = p_p1[k] * q_p2[k];
            r_p1_p3[k] = p_p1[k] * q_p3[k];
        }

        // NOTE: an angular component is accumulated by a loop of its own, as
        // it reads only the products which carry it and the vectorizer would
        // otherwise hold every product of the bidegree at once.

#pragma omp simd aligned(out_m4, acc_1, r_m1_p3, r_p1_m3 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m4[k] += acc_1[k] * (std::sqrt(14.0) * r_m1_p3[k] + std::sqrt(14.0) * r_p1_m3[k]);
        }

#pragma omp simd aligned(out_m3, acc_1, r_m1_p2, r_0_m3, r_p1_m2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m3[k] += acc_1[k] * (std::sqrt(10.5) * r_m1_p2[k] + std::sqrt(7.0) * r_0_m3[k] + std::sqrt(10.5) * r_p1_m2[k]);
        }

#pragma omp simd aligned(out_m2, acc_1, r_m1_p1, r_m1_p3, r_0_m2, r_p1_m3, r_p1_m1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m2[k] += acc_1[k] * (std::sqrt(7.5) * r_m1_p1[k] + std::sqrt(0.5) * r_m1_p3[k] + std::sqrt(12.0) * r_0_m2[k] - std::sqrt(0.5) * r_p1_m3[k] + std::sqrt(7.5) * r_p1_m1[k]);
        }

#pragma omp simd aligned(out_m1, acc_1, r_m1_0, r_m1_p2, r_0_m1, r_p1_m2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m1[k] += acc_1[k] * (std::sqrt(10.0) * r_m1_0[k] + std::sqrt(1.5) * r_m1_p2[k] + std::sqrt(15.0) * r_0_m1[k] - std::sqrt(1.5) * r_p1_m2[k]);
        }

#pragma omp simd aligned(out_0, acc_1, r_m1_m1, r_0_0, r_p1_p1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_0[k] += acc_1[k] * (-std::sqrt(6.0) * r_m1_m1[k] + 4.0 * r_0_0[k] - std::sqrt(6.0) * r_p1_p1[k]);
        }

#pragma omp simd aligned(out_p1, acc_1, r_m1_m2, r_0_p1, r_p1_0, r_p1_p2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p1[k] += acc_1[k] * (-std::sqrt(1.5) * r_m1_m2[k] + std::sqrt(15.0) * r_0_p1[k] + std::sqrt(10.0) * r_p1_0[k] - std::sqrt(1.5) * r_p1_p2[k]);
        }

#pragma omp simd aligned(out_p2, acc_1, r_m1_m3, r_m1_m1, r_0_p2, r_p1_p1, r_p1_p3 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p2[k] += acc_1[k] * (-std::sqrt(0.5) * r_m1_m3[k] - std::sqrt(7.5) * r_m1_m1[k] + std::sqrt(12.0) * r_0_p2[k] + std::sqrt(7.5) * r_p1_p1[k] - std::sqrt(0.5) * r_p1_p3[k]);
        }

#pragma omp simd aligned(out_p3, acc_1, r_m1_m2, r_0_p3, r_p1_p2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p3[k] += acc_1[k] * (-std::sqrt(10.5) * r_m1_m2[k] + std::sqrt(7.0) * r_0_p3[k] + std::sqrt(10.5) * r_p1_p2[k]);
        }

#pragma omp simd aligned(out_p4, acc_1, r_m1_m3, r_p1_p3 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p4[k] += acc_1[k] * (-std::sqrt(14.0) * r_m1_m3[k] + std::sqrt(14.0) * r_p1_p3[k]);
        }
    }

    // the bidegree of degree two on the atom pairs and two on the atoms
    // on c side, whose 25 products of harmonics are formed once and read by
    // the angular components which carry them

    {
        auto products = CSimdMatrix(25, nmax);

        const auto *p_m2 = ab_harmonics[1].data(0);
        const auto *p_m1 = ab_harmonics[1].data(1);
        const auto *p_0 = ab_harmonics[1].data(2);
        const auto *p_p1 = ab_harmonics[1].data(3);
        const auto *p_p2 = ab_harmonics[1].data(4);

        const auto *q_m2 = bc_harmonics[1].data(0);
        const auto *q_m1 = bc_harmonics[1].data(1);
        const auto *q_0 = bc_harmonics[1].data(2);
        const auto *q_p1 = bc_harmonics[1].data(3);
        const auto *q_p2 = bc_harmonics[1].data(4);

        auto *r_m2_m2 = products.data(0);
        auto *r_m2_m1 = products.data(1);
        auto *r_m2_0 = products.data(2);
        auto *r_m2_p1 = products.data(3);
        auto *r_m2_p2 = products.data(4);
        auto *r_m1_m2 = products.data(5);
        auto *r_m1_m1 = products.data(6);
        auto *r_m1_0 = products.data(7);
        auto *r_m1_p1 = products.data(8);
        auto *r_m1_p2 = products.data(9);
        auto *r_0_m2 = products.data(10);
        auto *r_0_m1 = products.data(11);
        auto *r_0_0 = products.data(12);
        auto *r_0_p1 = products.data(13);
        auto *r_0_p2 = products.data(14);
        auto *r_p1_m2 = products.data(15);
        auto *r_p1_m1 = products.data(16);
        auto *r_p1_0 = products.data(17);
        auto *r_p1_p1 = products.data(18);
        auto *r_p1_p2 = products.data(19);
        auto *r_p2_m2 = products.data(20);
        auto *r_p2_m1 = products.data(21);
        auto *r_p2_0 = products.data(22);
        auto *r_p2_p1 = products.data(23);
        auto *r_p2_p2 = products.data(24);

#pragma omp simd aligned(r_m2_m2, r_m2_m1, r_m2_0, r_m2_p1, r_m2_p2, r_m1_m2, r_m1_m1, r_m1_0, r_m1_p1, r_m1_p2, r_0_m2, r_0_m1, r_0_0, r_0_p1, r_0_p2, r_p1_m2, r_p1_m1, r_p1_0, r_p1_p1, r_p1_p2, r_p2_m2, r_p2_m1, r_p2_0, r_p2_p1, r_p2_p2, p_m2, p_m1, p_0, p_p1, p_p2, q_m2, q_m1, q_0, q_p1, q_p2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            r_m2_m2[k] = p_m2[k] * q_m2[k];
            r_m2_m1[k] = p_m2[k] * q_m1[k];
            r_m2_0[k] = p_m2[k] * q_0[k];
            r_m2_p1[k] = p_m2[k] * q_p1[k];
            r_m2_p2[k] = p_m2[k] * q_p2[k];
            r_m1_m2[k] = p_m1[k] * q_m2[k];
            r_m1_m1[k] = p_m1[k] * q_m1[k];
            r_m1_0[k] = p_m1[k] * q_0[k];
            r_m1_p1[k] = p_m1[k] * q_p1[k];
            r_m1_p2[k] = p_m1[k] * q_p2[k];
            r_0_m2[k] = p_0[k] * q_m2[k];
            r_0_m1[k] = p_0[k] * q_m1[k];
            r_0_0[k] = p_0[k] * q_0[k];
            r_0_p1[k] = p_0[k] * q_p1[k];
            r_0_p2[k] = p_0[k] * q_p2[k];
            r_p1_m2[k] = p_p1[k] * q_m2[k];
            r_p1_m1[k] = p_p1[k] * q_m1[k];
            r_p1_0[k] = p_p1[k] * q_0[k];
            r_p1_p1[k] = p_p1[k] * q_p1[k];
            r_p1_p2[k] = p_p1[k] * q_p2[k];
            r_p2_m2[k] = p_p2[k] * q_m2[k];
            r_p2_m1[k] = p_p2[k] * q_m1[k];
            r_p2_0[k] = p_p2[k] * q_0[k];
            r_p2_p1[k] = p_p2[k] * q_p1[k];
            r_p2_p2[k] = p_p2[k] * q_p2[k];
        }

        // NOTE: an angular component is accumulated by a loop of its own, as
        // it reads only the products which carry it and the vectorizer would
        // otherwise hold every product of the bidegree at once.

#pragma omp simd aligned(out_m4, acc_2, r_m2_p2, r_p2_m2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m4[k] += acc_2[k] * (std::sqrt(35.0) * r_m2_p2[k] + std::sqrt(35.0) * r_p2_m2[k]);
        }

#pragma omp simd aligned(out_m3, acc_2, r_m2_p1, r_m1_p2, r_p1_m2, r_p2_m1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m3[k] += acc_2[k] * (std::sqrt(17.5) * r_m2_p1[k] + std::sqrt(17.5) * r_m1_p2[k] + std::sqrt(17.5) * r_p1_m2[k] + std::sqrt(17.5) * r_p2_m1[k]);
        }

#pragma omp simd aligned(out_m2, acc_2, r_m2_0, r_m1_p1, r_0_m2, r_p1_m1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m2[k] += acc_2[k] * (std::sqrt(15.0) * r_m2_0[k] + std::sqrt(20.0) * r_m1_p1[k] + std::sqrt(15.0) * r_0_m2[k] + std::sqrt(20.0) * r_p1_m1[k]);
        }

#pragma omp simd aligned(out_m1, acc_2, r_m2_p1, r_m1_0, r_m1_p2, r_0_m1, r_p1_m2, r_p2_m1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m1[k] += acc_2[k] * (-std::sqrt(2.5) * r_m2_p1[k] + std::sqrt(30.0) * r_m1_0[k] + std::sqrt(2.5) * r_m1_p2[k] + std::sqrt(30.0) * r_0_m1[k] - std::sqrt(2.5) * r_p1_m2[k] + std::sqrt(2.5) * r_p2_m1[k]);
        }

#pragma omp simd aligned(out_0, acc_2, r_m2_m2, r_m1_m1, r_0_0, r_p1_p1, r_p2_p2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_0[k] += acc_2[k] * (r_m2_m2[k] - 4.0 * r_m1_m1[k] + 6.0 * r_0_0[k] - 4.0 * r_p1_p1[k] + r_p2_p2[k]);
        }

#pragma omp simd aligned(out_p1, acc_2, r_m2_m1, r_m1_m2, r_0_p1, r_p1_0, r_p1_p2, r_p2_p1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p1[k] += acc_2[k] * (-std::sqrt(2.5) * r_m2_m1[k] - std::sqrt(2.5) * r_m1_m2[k] + std::sqrt(30.0) * r_0_p1[k] + std::sqrt(30.0) * r_p1_0[k] - std::sqrt(2.5) * r_p1_p2[k] - std::sqrt(2.5) * r_p2_p1[k]);
        }

#pragma omp simd aligned(out_p2, acc_2, r_m1_m1, r_0_p2, r_p1_p1, r_p2_0 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p2[k] += acc_2[k] * (-std::sqrt(20.0) * r_m1_m1[k] + std::sqrt(15.0) * r_0_p2[k] + std::sqrt(20.0) * r_p1_p1[k] + std::sqrt(15.0) * r_p2_0[k]);
        }

#pragma omp simd aligned(out_p3, acc_2, r_m2_m1, r_m1_m2, r_p1_p2, r_p2_p1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p3[k] += acc_2[k] * (-std::sqrt(17.5) * r_m2_m1[k] - std::sqrt(17.5) * r_m1_m2[k] + std::sqrt(17.5) * r_p1_p2[k] + std::sqrt(17.5) * r_p2_p1[k]);
        }

#pragma omp simd aligned(out_p4, acc_2, r_m2_m2, r_p2_p2 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p4[k] += acc_2[k] * (-std::sqrt(35.0) * r_m2_m2[k] + std::sqrt(35.0) * r_p2_p2[k]);
        }
    }

    // the bidegree of degree three on the atom pairs and one on the atoms
    // on c side, whose 21 products of harmonics are formed once and read by
    // the angular components which carry them

    {
        auto products = CSimdMatrix(21, nmax);

        const auto *p_m3 = ab_harmonics[2].data(0);
        const auto *p_m2 = ab_harmonics[2].data(1);
        const auto *p_m1 = ab_harmonics[2].data(2);
        const auto *p_0 = ab_harmonics[2].data(3);
        const auto *p_p1 = ab_harmonics[2].data(4);
        const auto *p_p2 = ab_harmonics[2].data(5);
        const auto *p_p3 = ab_harmonics[2].data(6);

        const auto *q_m1 = bc_harmonics[0].data(0);
        const auto *q_0 = bc_harmonics[0].data(1);
        const auto *q_p1 = bc_harmonics[0].data(2);

        auto *r_m3_m1 = products.data(0);
        auto *r_m3_0 = products.data(1);
        auto *r_m3_p1 = products.data(2);
        auto *r_m2_m1 = products.data(3);
        auto *r_m2_0 = products.data(4);
        auto *r_m2_p1 = products.data(5);
        auto *r_m1_m1 = products.data(6);
        auto *r_m1_0 = products.data(7);
        auto *r_m1_p1 = products.data(8);
        auto *r_0_m1 = products.data(9);
        auto *r_0_0 = products.data(10);
        auto *r_0_p1 = products.data(11);
        auto *r_p1_m1 = products.data(12);
        auto *r_p1_0 = products.data(13);
        auto *r_p1_p1 = products.data(14);
        auto *r_p2_m1 = products.data(15);
        auto *r_p2_0 = products.data(16);
        auto *r_p2_p1 = products.data(17);
        auto *r_p3_m1 = products.data(18);
        auto *r_p3_0 = products.data(19);
        auto *r_p3_p1 = products.data(20);

#pragma omp simd aligned(r_m3_m1, r_m3_0, r_m3_p1, r_m2_m1, r_m2_0, r_m2_p1, r_m1_m1, r_m1_0, r_m1_p1, r_0_m1, r_0_0, r_0_p1, r_p1_m1, r_p1_0, r_p1_p1, r_p2_m1, r_p2_0, r_p2_p1, r_p3_m1, r_p3_0, r_p3_p1, p_m3, p_m2, p_m1, p_0, p_p1, p_p2, p_p3, q_m1, q_0, q_p1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            r_m3_m1[k] = p_m3[k] * q_m1[k];
            r_m3_0[k] = p_m3[k] * q_0[k];
            r_m3_p1[k] = p_m3[k] * q_p1[k];
            r_m2_m1[k] = p_m2[k] * q_m1[k];
            r_m2_0[k] = p_m2[k] * q_0[k];
            r_m2_p1[k] = p_m2[k] * q_p1[k];
            r_m1_m1[k] = p_m1[k] * q_m1[k];
            r_m1_0[k] = p_m1[k] * q_0[k];
            r_m1_p1[k] = p_m1[k] * q_p1[k];
            r_0_m1[k] = p_0[k] * q_m1[k];
            r_0_0[k] = p_0[k] * q_0[k];
            r_0_p1[k] = p_0[k] * q_p1[k];
            r_p1_m1[k] = p_p1[k] * q_m1[k];
            r_p1_0[k] = p_p1[k] * q_0[k];
            r_p1_p1[k] = p_p1[k] * q_p1[k];
            r_p2_m1[k] = p_p2[k] * q_m1[k];
            r_p2_0[k] = p_p2[k] * q_0[k];
            r_p2_p1[k] = p_p2[k] * q_p1[k];
            r_p3_m1[k] = p_p3[k] * q_m1[k];
            r_p3_0[k] = p_p3[k] * q_0[k];
            r_p3_p1[k] = p_p3[k] * q_p1[k];
        }

        // NOTE: an angular component is accumulated by a loop of its own, as
        // it reads only the products which carry it and the vectorizer would
        // otherwise hold every product of the bidegree at once.

#pragma omp simd aligned(out_m4, acc_3, r_m3_p1, r_p3_m1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m4[k] += acc_3[k] * (std::sqrt(14.0) * r_m3_p1[k] + std::sqrt(14.0) * r_p3_m1[k]);
        }

#pragma omp simd aligned(out_m3, acc_3, r_m3_0, r_m2_p1, r_p2_m1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m3[k] += acc_3[k] * (std::sqrt(7.0) * r_m3_0[k] + std::sqrt(10.5) * r_m2_p1[k] + std::sqrt(10.5) * r_p2_m1[k]);
        }

#pragma omp simd aligned(out_m2, acc_3, r_m3_p1, r_m2_0, r_m1_p1, r_p1_m1, r_p3_m1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m2[k] += acc_3[k] * (-std::sqrt(0.5) * r_m3_p1[k] + std::sqrt(12.0) * r_m2_0[k] + std::sqrt(7.5) * r_m1_p1[k] + std::sqrt(7.5) * r_p1_m1[k] + std::sqrt(0.5) * r_p3_m1[k]);
        }

#pragma omp simd aligned(out_m1, acc_3, r_m2_p1, r_m1_0, r_0_m1, r_p2_m1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_m1[k] += acc_3[k] * (-std::sqrt(1.5) * r_m2_p1[k] + std::sqrt(15.0) * r_m1_0[k] + std::sqrt(10.0) * r_0_m1[k] + std::sqrt(1.5) * r_p2_m1[k]);
        }

#pragma omp simd aligned(out_0, acc_3, r_m1_m1, r_0_0, r_p1_p1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_0[k] += acc_3[k] * (-std::sqrt(6.0) * r_m1_m1[k] + 4.0 * r_0_0[k] - std::sqrt(6.0) * r_p1_p1[k]);
        }

#pragma omp simd aligned(out_p1, acc_3, r_m2_m1, r_0_p1, r_p1_0, r_p2_p1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p1[k] += acc_3[k] * (-std::sqrt(1.5) * r_m2_m1[k] + std::sqrt(10.0) * r_0_p1[k] + std::sqrt(15.0) * r_p1_0[k] - std::sqrt(1.5) * r_p2_p1[k]);
        }

#pragma omp simd aligned(out_p2, acc_3, r_m3_m1, r_m1_m1, r_p1_p1, r_p2_0, r_p3_p1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p2[k] += acc_3[k] * (-std::sqrt(0.5) * r_m3_m1[k] - std::sqrt(7.5) * r_m1_m1[k] + std::sqrt(7.5) * r_p1_p1[k] + std::sqrt(12.0) * r_p2_0[k] - std::sqrt(0.5) * r_p3_p1[k]);
        }

#pragma omp simd aligned(out_p3, acc_3, r_m2_m1, r_p2_p1, r_p3_0 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p3[k] += acc_3[k] * (-std::sqrt(10.5) * r_m2_m1[k] + std::sqrt(10.5) * r_p2_p1[k] + std::sqrt(7.0) * r_p3_0[k]);
        }

#pragma omp simd aligned(out_p4, acc_3, r_m3_m1, r_p3_p1 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            out_p4[k] += acc_3[k] * (-std::sqrt(14.0) * r_m3_m1[k] + std::sqrt(14.0) * r_p3_p1[k]);
        }
    }

    // the bidegree of degree four on the atom pairs and zero on the atoms
    // on c side, whose coefficients are one: the harmonic of the other side is
    // of degree zero and is one for every atom pair

    {
        const auto *h_m4 = ab_harmonics[3].data(0);
        const auto *h_m3 = ab_harmonics[3].data(1);
        const auto *h_m2 = ab_harmonics[3].data(2);
        const auto *h_m1 = ab_harmonics[3].data(3);
        const auto *h_0 = ab_harmonics[3].data(4);
        const auto *h_p1 = ab_harmonics[3].data(5);
        const auto *h_p2 = ab_harmonics[3].data(6);
        const auto *h_p3 = ab_harmonics[3].data(7);
        const auto *h_p4 = ab_harmonics[3].data(8);

#pragma omp simd aligned(out_m4, out_m3, out_m2, out_m1, out_0, out_p1, out_p2, out_p3, out_p4, acc_4, h_m4, h_m3, h_m2, h_m1, h_0, h_p1, h_p2, h_p3, h_p4 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            const auto f = acc_4[k];

            out_m4[k] += f * h_m4[k];
            out_m3[k] += f * h_m3[k];
            out_m2[k] += f * h_m2[k];
            out_m1[k] += f * h_m1[k];
            out_0[k] += f * h_0[k];
            out_p1[k] += f * h_p1[k];
            out_p2[k] += f * h_p2[k];
            out_p3[k] += f * h_p3[k];
            out_p4[k] += f * h_p4[k];
        }
    }

    // NOTE: the atom pairs beyond the reach of every triple of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 9; m++)
    {
        const auto *row = components.data(m);

        std::copy(row, row + nmax, slices[m]);

        std::fill(slices[m] + nmax, slices[m] + npairs, 0.0);
    }
}

}  // namespace simdt3ceri
