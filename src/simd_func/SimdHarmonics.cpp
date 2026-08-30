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



#include "SimdHarmonics.hpp"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

namespace simdfunc {  // simdfunc namespace

auto
make_solid_harmonics(const CSimdMatrix &coordinates, const int lmax) -> std::vector<CSimdMatrix>
{
    errors::assertMsgCritical(lmax <= 12,
                              std::string("SimdHarmonics.make_solid_harmonics: Solid harmonics above angular momentum twelve are not implemented"));

    std::vector<CSimdMatrix> harmonics;

    if (lmax < 1) return harmonics;

    // NOTE: the storage is reserved before the harmonics are created, so that the
    // matrices of the lower angular momenta are not moved while they are read by
    // the recursion.

    harmonics.reserve(static_cast<size_t>(lmax));

    harmonics.push_back(make_p_solid_harmonics(coordinates));

    if (lmax > 1) harmonics.push_back(make_d_solid_harmonics(harmonics[0], coordinates));

    if (lmax > 2) harmonics.push_back(make_f_solid_harmonics(harmonics[1], harmonics[0], coordinates));

    if (lmax > 3) harmonics.push_back(make_g_solid_harmonics(harmonics[2], harmonics[1], harmonics[0], coordinates));

    if (lmax > 4) harmonics.push_back(make_h_solid_harmonics(harmonics[3], harmonics[2], harmonics[0], coordinates));

    if (lmax > 5) harmonics.push_back(make_i_solid_harmonics(harmonics[4], harmonics[3], harmonics[0], coordinates));

    if (lmax > 6) harmonics.push_back(make_k_solid_harmonics(harmonics[5], harmonics[4], harmonics[0], coordinates));

    if (lmax > 7) harmonics.push_back(make_l_solid_harmonics(harmonics[6], harmonics[5], harmonics[0], coordinates));

    if (lmax > 8) harmonics.push_back(make_m_solid_harmonics(harmonics[7], harmonics[6], harmonics[0], coordinates));

    if (lmax > 9) harmonics.push_back(make_n_solid_harmonics(harmonics[8], harmonics[7], harmonics[0], coordinates));

    if (lmax > 10) harmonics.push_back(make_o_solid_harmonics(harmonics[9], harmonics[8], harmonics[0], coordinates));

    if (lmax > 11) harmonics.push_back(make_q_solid_harmonics(harmonics[10], harmonics[9], harmonics[0], coordinates));

    return harmonics;
}

auto
make_p_solid_harmonics(const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_p_solid_harmonics: Coordinates must have seven rows"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(3, ncols);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *a_x = coordinates.data(0);
    const auto *a_y = coordinates.data(1);
    const auto *a_z = coordinates.data(2);
    const auto *b_x = coordinates.data(3);
    const auto *b_y = coordinates.data(4);
    const auto *b_z = coordinates.data(5);

    auto *s_m1 = matrix.data(0);
    auto *s_z0 = matrix.data(1);
    auto *s_p1 = matrix.data(2);

    // NOTE: the real solid harmonics of angular momentum one are the components
    // of the vector between the atoms, taken in the order of the spherical
    // components, and need no normalization factor.

    // NOTE: the rows of the coordinates and of the harmonics start at a cache
    // line boundary, so the loop is vectorized with aligned loads and stores.

#pragma omp simd aligned(s_m1, s_z0, s_p1, a_x, a_y, a_z, b_x, b_y, b_z : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        s_m1[i] = a_y[i] - b_y[i];

        s_z0[i] = a_z[i] - b_z[i];

        s_p1[i] = a_x[i] - b_x[i];
    }

    return matrix;
}

auto
make_d_solid_harmonics(const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_d_solid_harmonics: Harmonics of angular momentum one must have 3 rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_d_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(p_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_d_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(5, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = 0.5 * std::sqrt(3.0);
    const auto fact_xy = std::sqrt(3.0);
    const auto fz_m1 = std::sqrt(3.0);
    const auto fz_0 = 1.5;
    const auto fr_0 = 0.5;
    const auto fz_p1 = std::sqrt(3.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    auto *ps_m2 = matrix.data(0);
    auto *ps_m1 = matrix.data(1);
    auto *ps_0 = matrix.data(2);
    auto *ps_p1 = matrix.data(3);
    auto *ps_p2 = matrix.data(4);

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, ps_m2, ps_m1, ps_0, ps_p1, ps_p2 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        ps_p2[i] = fact * (x * x - y * y);

        ps_m2[i] = fact_xy * x * y;

        ps_m1[i] = fz_m1 * z * y;

        ps_0[i] = fz_0 * z * z - fr_0 * r_2;

        ps_p1[i] = fz_p1 * z * x;
    }

    return matrix;
}

auto
make_f_solid_harmonics(const CSimdMatrix &d_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(d_harmonics.number_of_rows() == 5,
                              std::string("SimdHarmonics.make_f_solid_harmonics: Harmonics of angular momentum two must have 5 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_f_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_f_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(d_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_f_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(7, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = (1.0 / 6.0) * std::sqrt(30.0);
    const auto fz_m2 = std::sqrt(5.0);
    const auto fz_m1 = 1.25 * std::sqrt(2.0);
    const auto fr_m1 = 0.25 * std::sqrt(6.0);
    const auto fz_0 = 5.0 / 3.0;
    const auto fr_0 = 2.0 / 3.0;
    const auto fz_p1 = 1.25 * std::sqrt(2.0);
    const auto fr_p1 = 0.25 * std::sqrt(6.0);
    const auto fz_p2 = std::sqrt(5.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m2 = d_harmonics.data(0);
    const auto *pt_m1 = d_harmonics.data(1);
    const auto *pt_0 = d_harmonics.data(2);
    const auto *pt_p1 = d_harmonics.data(3);
    const auto *pt_p2 = d_harmonics.data(4);

    auto *ps_m3 = matrix.data(0);
    auto *ps_m2 = matrix.data(1);
    auto *ps_m1 = matrix.data(2);
    auto *ps_0 = matrix.data(3);
    auto *ps_p1 = matrix.data(4);
    auto *ps_p2 = matrix.data(5);
    auto *ps_p3 = matrix.data(6);

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];

        ps_p3[i] = fact * (x * t_p2 - y * t_m2);

        ps_m3[i] = fact * (y * t_p2 + x * t_m2);

        ps_m2[i] = fz_m2 * z * t_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * y;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * z;

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * x;

        ps_p2[i] = fz_p2 * z * t_p2;
    }

    return matrix;
}

auto
make_g_solid_harmonics(const CSimdMatrix &f_harmonics, const CSimdMatrix &d_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(f_harmonics.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_g_solid_harmonics: Harmonics of angular momentum three must have 7 rows"));

    errors::assertMsgCritical(d_harmonics.number_of_rows() == 5,
                              std::string("SimdHarmonics.make_g_solid_harmonics: Harmonics of angular momentum two must have 5 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_g_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_g_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(f_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_g_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(9, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = 0.25 * std::sqrt(14.0);
    const auto fz_m3 = std::sqrt(7.0);
    const auto fz_m2 = (7.0 / 6.0) * std::sqrt(3.0);
    const auto fr_m2 = (1.0 / 6.0) * std::sqrt(15.0);
    const auto fz_m1 = (7.0 / 15.0) * std::sqrt(15.0);
    const auto fr_m1 = (2.0 / 15.0) * std::sqrt(30.0);
    const auto fz_0 = 1.75;
    const auto fr_0 = 0.75;
    const auto fz_p1 = (7.0 / 15.0) * std::sqrt(15.0);
    const auto fr_p1 = (2.0 / 15.0) * std::sqrt(30.0);
    const auto fz_p2 = (7.0 / 6.0) * std::sqrt(3.0);
    const auto fr_p2 = (1.0 / 6.0) * std::sqrt(15.0);
    const auto fz_p3 = std::sqrt(7.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m3 = f_harmonics.data(0);
    const auto *pt_m2 = f_harmonics.data(1);
    const auto *pt_m1 = f_harmonics.data(2);
    const auto *pt_0 = f_harmonics.data(3);
    const auto *pt_p1 = f_harmonics.data(4);
    const auto *pt_p2 = f_harmonics.data(5);
    const auto *pt_p3 = f_harmonics.data(6);

    const auto *pu_m2 = d_harmonics.data(0);
    const auto *pu_m1 = d_harmonics.data(1);
    const auto *pu_0 = d_harmonics.data(2);
    const auto *pu_p1 = d_harmonics.data(3);
    const auto *pu_p2 = d_harmonics.data(4);

    auto *ps_m4 = matrix.data(0);
    auto *ps_m3 = matrix.data(1);
    auto *ps_m2 = matrix.data(2);
    auto *ps_m1 = matrix.data(3);
    auto *ps_0 = matrix.data(4);
    auto *ps_p1 = matrix.data(5);
    auto *ps_p2 = matrix.data(6);
    auto *ps_p3 = matrix.data(7);
    auto *ps_p4 = matrix.data(8);

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m3, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, pt_p3, pu_m2, pu_m1, pu_0, pu_p1, pu_p2, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3, ps_p4 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];

        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];
        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];

        ps_p4[i] = fact * (x * t_p3 - y * t_m3);

        ps_m4[i] = fact * (y * t_p3 + x * t_m3);

        ps_m3[i] = fz_m3 * z * t_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3;
    }

    return matrix;
}

auto
make_h_solid_harmonics(const CSimdMatrix &g_harmonics, const CSimdMatrix &f_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(g_harmonics.number_of_rows() == 9,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics of angular momentum four must have 9 rows"));

    errors::assertMsgCritical(f_harmonics.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics of angular momentum three must have 7 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(g_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(11, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = (3.0 / 10.0) * std::sqrt(10.0);
    const auto fz_m4 = 3.0;
    const auto fz_m3 = 2.25;
    const auto fr_m3 = 0.25 * std::sqrt(7.0);
    const auto fz_m2 = (3.0 / 7.0) * std::sqrt(21.0);
    const auto fr_m2 = (2.0 / 7.0) * std::sqrt(7.0);
    const auto fz_m1 = 0.75 * std::sqrt(6.0);
    const auto fr_m1 = 0.25 * std::sqrt(10.0);
    const auto fz_0 = 9.0 / 5.0;
    const auto fr_0 = 4.0 / 5.0;
    const auto fz_p1 = 0.75 * std::sqrt(6.0);
    const auto fr_p1 = 0.25 * std::sqrt(10.0);
    const auto fz_p2 = (3.0 / 7.0) * std::sqrt(21.0);
    const auto fr_p2 = (2.0 / 7.0) * std::sqrt(7.0);
    const auto fz_p3 = 2.25;
    const auto fr_p3 = 0.25 * std::sqrt(7.0);
    const auto fz_p4 = 3.0;

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m4 = g_harmonics.data(0);
    const auto *pt_m3 = g_harmonics.data(1);
    const auto *pt_m2 = g_harmonics.data(2);
    const auto *pt_m1 = g_harmonics.data(3);
    const auto *pt_0 = g_harmonics.data(4);
    const auto *pt_p1 = g_harmonics.data(5);
    const auto *pt_p2 = g_harmonics.data(6);
    const auto *pt_p3 = g_harmonics.data(7);
    const auto *pt_p4 = g_harmonics.data(8);

    const auto *pu_m3 = f_harmonics.data(0);
    const auto *pu_m2 = f_harmonics.data(1);
    const auto *pu_m1 = f_harmonics.data(2);
    const auto *pu_0 = f_harmonics.data(3);
    const auto *pu_p1 = f_harmonics.data(4);
    const auto *pu_p2 = f_harmonics.data(5);
    const auto *pu_p3 = f_harmonics.data(6);

    auto *ps_m5 = matrix.data(0);
    auto *ps_m4 = matrix.data(1);
    auto *ps_m3 = matrix.data(2);
    auto *ps_m2 = matrix.data(3);
    auto *ps_m1 = matrix.data(4);
    auto *ps_0 = matrix.data(5);
    auto *ps_p1 = matrix.data(6);
    auto *ps_p2 = matrix.data(7);
    auto *ps_p3 = matrix.data(8);
    auto *ps_p4 = matrix.data(9);
    auto *ps_p5 = matrix.data(10);

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m4, pt_m3, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, pt_p3, pt_p4, pu_m3, pu_m2, pu_m1, pu_0, pu_p1, pu_p2, pu_p3, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m4 = pt_m4[i];
        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];

        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];
        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];

        ps_p5[i] = fact * (x * t_p4 - y * t_m4);

        ps_m5[i] = fact * (y * t_p4 + x * t_m4);

        ps_m4[i] = fz_m4 * z * t_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4;
    }

    return matrix;
}

auto
make_i_solid_harmonics(const CSimdMatrix &h_harmonics, const CSimdMatrix &g_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(h_harmonics.number_of_rows() == 11,
                              std::string("SimdHarmonics.make_i_solid_harmonics: Harmonics of angular momentum five must have 11 rows"));

    errors::assertMsgCritical(g_harmonics.number_of_rows() == 9,
                              std::string("SimdHarmonics.make_i_solid_harmonics: Harmonics of angular momentum four must have 9 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_i_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_i_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(h_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_i_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(13, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = (1.0 / 6.0) * std::sqrt(33.0);
    const auto fz_m5 = std::sqrt(11.0);
    const auto fz_m4 = (11.0 / 10.0) * std::sqrt(5.0);
    const auto fr_m4 = (3.0 / 10.0) * std::sqrt(5.0);
    const auto fz_m3 = (11.0 / 9.0) * std::sqrt(3.0);
    const auto fr_m3 = (4.0 / 9.0) * std::sqrt(3.0);
    const auto fz_m2 = 1.375 * std::sqrt(2.0);
    const auto fr_m2 = 0.125 * std::sqrt(42.0);
    const auto fz_m1 = (11.0 / 35.0) * std::sqrt(35.0);
    const auto fr_m1 = (2.0 / 35.0) * std::sqrt(210.0);
    const auto fz_0 = 11.0 / 6.0;
    const auto fr_0 = 5.0 / 6.0;
    const auto fz_p1 = (11.0 / 35.0) * std::sqrt(35.0);
    const auto fr_p1 = (2.0 / 35.0) * std::sqrt(210.0);
    const auto fz_p2 = 1.375 * std::sqrt(2.0);
    const auto fr_p2 = 0.125 * std::sqrt(42.0);
    const auto fz_p3 = (11.0 / 9.0) * std::sqrt(3.0);
    const auto fr_p3 = (4.0 / 9.0) * std::sqrt(3.0);
    const auto fz_p4 = (11.0 / 10.0) * std::sqrt(5.0);
    const auto fr_p4 = (3.0 / 10.0) * std::sqrt(5.0);
    const auto fz_p5 = std::sqrt(11.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m5 = h_harmonics.data(0);
    const auto *pt_m4 = h_harmonics.data(1);
    const auto *pt_m3 = h_harmonics.data(2);
    const auto *pt_m2 = h_harmonics.data(3);
    const auto *pt_m1 = h_harmonics.data(4);
    const auto *pt_0 = h_harmonics.data(5);
    const auto *pt_p1 = h_harmonics.data(6);
    const auto *pt_p2 = h_harmonics.data(7);
    const auto *pt_p3 = h_harmonics.data(8);
    const auto *pt_p4 = h_harmonics.data(9);
    const auto *pt_p5 = h_harmonics.data(10);

    const auto *pu_m4 = g_harmonics.data(0);
    const auto *pu_m3 = g_harmonics.data(1);
    const auto *pu_m2 = g_harmonics.data(2);
    const auto *pu_m1 = g_harmonics.data(3);
    const auto *pu_0 = g_harmonics.data(4);
    const auto *pu_p1 = g_harmonics.data(5);
    const auto *pu_p2 = g_harmonics.data(6);
    const auto *pu_p3 = g_harmonics.data(7);
    const auto *pu_p4 = g_harmonics.data(8);

    auto *ps_m6 = matrix.data(0);
    auto *ps_m5 = matrix.data(1);
    auto *ps_m4 = matrix.data(2);
    auto *ps_m3 = matrix.data(3);
    auto *ps_m2 = matrix.data(4);
    auto *ps_m1 = matrix.data(5);
    auto *ps_0 = matrix.data(6);
    auto *ps_p1 = matrix.data(7);
    auto *ps_p2 = matrix.data(8);
    auto *ps_p3 = matrix.data(9);
    auto *ps_p4 = matrix.data(10);
    auto *ps_p5 = matrix.data(11);
    auto *ps_p6 = matrix.data(12);

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m5, pt_m4, pt_m3, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, pt_p3, pt_p4, pt_p5, pu_m4, pu_m3, pu_m2, pu_m1, pu_0, pu_p1, pu_p2, pu_p3, pu_p4, ps_m6, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5, ps_p6 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m5 = pt_m5[i];
        const auto t_m4 = pt_m4[i];
        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];
        const auto t_p5 = pt_p5[i];

        const auto u_m4 = pu_m4[i];
        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];
        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];
        const auto u_p4 = pu_p4[i];

        ps_p6[i] = fact * (x * t_p5 - y * t_m5);

        ps_m6[i] = fact * (y * t_p5 + x * t_m5);

        ps_m5[i] = fz_m5 * z * t_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4 - fr_p4 * r_2 * u_p4;

        ps_p5[i] = fz_p5 * z * t_p5;
    }

    return matrix;
}

auto
make_k_solid_harmonics(const CSimdMatrix &i_harmonics, const CSimdMatrix &h_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(i_harmonics.number_of_rows() == 13,
                              std::string("SimdHarmonics.make_k_solid_harmonics: Harmonics of angular momentum six must have 13 rows"));

    errors::assertMsgCritical(h_harmonics.number_of_rows() == 11,
                              std::string("SimdHarmonics.make_k_solid_harmonics: Harmonics of angular momentum five must have 11 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_k_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_k_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(i_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_k_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(15, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = (1.0 / 14.0) * std::sqrt(182.0);
    const auto fz_m6 = std::sqrt(13.0);
    const auto fz_m5 = (13.0 / 12.0) * std::sqrt(6.0);
    const auto fr_m5 = (1.0 / 12.0) * std::sqrt(66.0);
    const auto fz_m4 = (13.0 / 33.0) * std::sqrt(33.0);
    const auto fr_m4 = (2.0 / 33.0) * std::sqrt(165.0);
    const auto fz_m3 = (13.0 / 20.0) * std::sqrt(10.0);
    const auto fr_m3 = (3.0 / 20.0) * std::sqrt(30.0);
    const auto fz_m2 = (13.0 / 15.0) * std::sqrt(5.0);
    const auto fr_m2 = (4.0 / 15.0) * std::sqrt(10.0);
    const auto fz_m1 = (13.0 / 12.0) * std::sqrt(3.0);
    const auto fr_m1 = (1.0 / 12.0) * std::sqrt(105.0);
    const auto fz_0 = 13.0 / 7.0;
    const auto fr_0 = 6.0 / 7.0;
    const auto fz_p1 = (13.0 / 12.0) * std::sqrt(3.0);
    const auto fr_p1 = (1.0 / 12.0) * std::sqrt(105.0);
    const auto fz_p2 = (13.0 / 15.0) * std::sqrt(5.0);
    const auto fr_p2 = (4.0 / 15.0) * std::sqrt(10.0);
    const auto fz_p3 = (13.0 / 20.0) * std::sqrt(10.0);
    const auto fr_p3 = (3.0 / 20.0) * std::sqrt(30.0);
    const auto fz_p4 = (13.0 / 33.0) * std::sqrt(33.0);
    const auto fr_p4 = (2.0 / 33.0) * std::sqrt(165.0);
    const auto fz_p5 = (13.0 / 12.0) * std::sqrt(6.0);
    const auto fr_p5 = (1.0 / 12.0) * std::sqrt(66.0);
    const auto fz_p6 = std::sqrt(13.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m6 = i_harmonics.data(0);
    const auto *pt_m5 = i_harmonics.data(1);
    const auto *pt_m4 = i_harmonics.data(2);
    const auto *pt_m3 = i_harmonics.data(3);
    const auto *pt_m2 = i_harmonics.data(4);
    const auto *pt_m1 = i_harmonics.data(5);
    const auto *pt_0 = i_harmonics.data(6);
    const auto *pt_p1 = i_harmonics.data(7);
    const auto *pt_p2 = i_harmonics.data(8);
    const auto *pt_p3 = i_harmonics.data(9);
    const auto *pt_p4 = i_harmonics.data(10);
    const auto *pt_p5 = i_harmonics.data(11);
    const auto *pt_p6 = i_harmonics.data(12);

    const auto *pu_m5 = h_harmonics.data(0);
    const auto *pu_m4 = h_harmonics.data(1);
    const auto *pu_m3 = h_harmonics.data(2);
    const auto *pu_m2 = h_harmonics.data(3);
    const auto *pu_m1 = h_harmonics.data(4);
    const auto *pu_0 = h_harmonics.data(5);
    const auto *pu_p1 = h_harmonics.data(6);
    const auto *pu_p2 = h_harmonics.data(7);
    const auto *pu_p3 = h_harmonics.data(8);
    const auto *pu_p4 = h_harmonics.data(9);
    const auto *pu_p5 = h_harmonics.data(10);

    auto *ps_m7 = matrix.data(0);
    auto *ps_m6 = matrix.data(1);
    auto *ps_m5 = matrix.data(2);
    auto *ps_m4 = matrix.data(3);
    auto *ps_m3 = matrix.data(4);
    auto *ps_m2 = matrix.data(5);
    auto *ps_m1 = matrix.data(6);
    auto *ps_0 = matrix.data(7);
    auto *ps_p1 = matrix.data(8);
    auto *ps_p2 = matrix.data(9);
    auto *ps_p3 = matrix.data(10);
    auto *ps_p4 = matrix.data(11);
    auto *ps_p5 = matrix.data(12);
    auto *ps_p6 = matrix.data(13);
    auto *ps_p7 = matrix.data(14);

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m6, pt_m5, pt_m4, pt_m3, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, pt_p3, pt_p4, pt_p5, pt_p6, pu_m5, pu_m4, pu_m3, pu_m2, pu_m1, pu_0, pu_p1, pu_p2, pu_p3, pu_p4, pu_p5, ps_m7, ps_m6, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5, ps_p6, ps_p7 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m6 = pt_m6[i];
        const auto t_m5 = pt_m5[i];
        const auto t_m4 = pt_m4[i];
        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];
        const auto t_p5 = pt_p5[i];
        const auto t_p6 = pt_p6[i];

        const auto u_m5 = pu_m5[i];
        const auto u_m4 = pu_m4[i];
        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];
        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];
        const auto u_p4 = pu_p4[i];
        const auto u_p5 = pu_p5[i];

        ps_p7[i] = fact * (x * t_p6 - y * t_m6);

        ps_m7[i] = fact * (y * t_p6 + x * t_m6);

        ps_m6[i] = fz_m6 * z * t_m6;

        ps_m5[i] = fz_m5 * z * t_m5 - fr_m5 * r_2 * u_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4 - fr_p4 * r_2 * u_p4;

        ps_p5[i] = fz_p5 * z * t_p5 - fr_p5 * r_2 * u_p5;

        ps_p6[i] = fz_p6 * z * t_p6;
    }

    return matrix;
}

auto
make_l_solid_harmonics(const CSimdMatrix &k_harmonics, const CSimdMatrix &i_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(k_harmonics.number_of_rows() == 15,
                              std::string("SimdHarmonics.make_l_solid_harmonics: Harmonics of angular momentum seven must have 15 rows"));

    errors::assertMsgCritical(i_harmonics.number_of_rows() == 13,
                              std::string("SimdHarmonics.make_l_solid_harmonics: Harmonics of angular momentum six must have 13 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_l_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_l_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(k_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_l_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(17, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = 0.25 * std::sqrt(15.0);
    const auto fz_m7 = std::sqrt(15.0);
    const auto fz_m6 = (15.0 / 14.0) * std::sqrt(7.0);
    const auto fr_m6 = (1.0 / 14.0) * std::sqrt(91.0);
    const auto fz_m5 = (5.0 / 13.0) * std::sqrt(39.0);
    const auto fr_m5 = (2.0 / 13.0) * std::sqrt(26.0);
    const auto fz_m4 = 1.25 * std::sqrt(3.0);
    const auto fr_m4 = 0.25 * std::sqrt(11.0);
    const auto fz_m3 = (3.0 / 11.0) * std::sqrt(55.0);
    const auto fr_m3 = (2.0 / 11.0) * std::sqrt(22.0);
    const auto fz_m2 = 0.5 * std::sqrt(15.0);
    const auto fr_m2 = 0.5 * std::sqrt(3.0);
    const auto fz_m1 = (5.0 / 7.0) * std::sqrt(7.0);
    const auto fr_m1 = (4.0 / 21.0) * std::sqrt(21.0);
    const auto fz_0 = 1.875;
    const auto fr_0 = 0.875;
    const auto fz_p1 = (5.0 / 7.0) * std::sqrt(7.0);
    const auto fr_p1 = (4.0 / 21.0) * std::sqrt(21.0);
    const auto fz_p2 = 0.5 * std::sqrt(15.0);
    const auto fr_p2 = 0.5 * std::sqrt(3.0);
    const auto fz_p3 = (3.0 / 11.0) * std::sqrt(55.0);
    const auto fr_p3 = (2.0 / 11.0) * std::sqrt(22.0);
    const auto fz_p4 = 1.25 * std::sqrt(3.0);
    const auto fr_p4 = 0.25 * std::sqrt(11.0);
    const auto fz_p5 = (5.0 / 13.0) * std::sqrt(39.0);
    const auto fr_p5 = (2.0 / 13.0) * std::sqrt(26.0);
    const auto fz_p6 = (15.0 / 14.0) * std::sqrt(7.0);
    const auto fr_p6 = (1.0 / 14.0) * std::sqrt(91.0);
    const auto fz_p7 = std::sqrt(15.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m7 = k_harmonics.data(0);
    const auto *pt_m6 = k_harmonics.data(1);
    const auto *pt_m5 = k_harmonics.data(2);
    const auto *pt_m4 = k_harmonics.data(3);
    const auto *pt_m3 = k_harmonics.data(4);
    const auto *pt_m2 = k_harmonics.data(5);
    const auto *pt_m1 = k_harmonics.data(6);
    const auto *pt_0 = k_harmonics.data(7);
    const auto *pt_p1 = k_harmonics.data(8);
    const auto *pt_p2 = k_harmonics.data(9);
    const auto *pt_p3 = k_harmonics.data(10);
    const auto *pt_p4 = k_harmonics.data(11);
    const auto *pt_p5 = k_harmonics.data(12);
    const auto *pt_p6 = k_harmonics.data(13);
    const auto *pt_p7 = k_harmonics.data(14);

    const auto *pu_m6 = i_harmonics.data(0);
    const auto *pu_m5 = i_harmonics.data(1);
    const auto *pu_m4 = i_harmonics.data(2);
    const auto *pu_m3 = i_harmonics.data(3);
    const auto *pu_m2 = i_harmonics.data(4);
    const auto *pu_m1 = i_harmonics.data(5);
    const auto *pu_0 = i_harmonics.data(6);
    const auto *pu_p1 = i_harmonics.data(7);
    const auto *pu_p2 = i_harmonics.data(8);
    const auto *pu_p3 = i_harmonics.data(9);
    const auto *pu_p4 = i_harmonics.data(10);
    const auto *pu_p5 = i_harmonics.data(11);
    const auto *pu_p6 = i_harmonics.data(12);

    auto *ps_m8 = matrix.data(0);
    auto *ps_m7 = matrix.data(1);
    auto *ps_m6 = matrix.data(2);
    auto *ps_m5 = matrix.data(3);
    auto *ps_m4 = matrix.data(4);
    auto *ps_m3 = matrix.data(5);
    auto *ps_m2 = matrix.data(6);
    auto *ps_m1 = matrix.data(7);
    auto *ps_0 = matrix.data(8);
    auto *ps_p1 = matrix.data(9);
    auto *ps_p2 = matrix.data(10);
    auto *ps_p3 = matrix.data(11);
    auto *ps_p4 = matrix.data(12);
    auto *ps_p5 = matrix.data(13);
    auto *ps_p6 = matrix.data(14);
    auto *ps_p7 = matrix.data(15);
    auto *ps_p8 = matrix.data(16);

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 17 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m7, pt_m6, pt_m5, pt_m4, pt_m3, pt_m2, pt_m1, pt_0, pt_p7, pu_m6, pu_m5, pu_m4, pu_m3, pu_m2, pu_m1, pu_0, ps_m8, ps_m7, ps_m6, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p8 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m7 = pt_m7[i];
        const auto t_m6 = pt_m6[i];
        const auto t_m5 = pt_m5[i];
        const auto t_m4 = pt_m4[i];
        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p7 = pt_p7[i];

        const auto u_m6 = pu_m6[i];
        const auto u_m5 = pu_m5[i];
        const auto u_m4 = pu_m4[i];
        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];

        ps_p8[i] = fact * (x * t_p7 - y * t_m7);

        ps_m8[i] = fact * (y * t_p7 + x * t_m7);

        ps_m7[i] = fz_m7 * z * t_m7;

        ps_m6[i] = fz_m6 * z * t_m6 - fr_m6 * r_2 * u_m6;

        ps_m5[i] = fz_m5 * z * t_m5 - fr_m5 * r_2 * u_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;
    }

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 17 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_p1, pt_p2, pt_p3, pt_p4, pt_p5, pt_p6, pt_p7, pu_p1, pu_p2, pu_p3, pu_p4, pu_p5, pu_p6, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5, ps_p6, ps_p7 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];
        const auto t_p5 = pt_p5[i];
        const auto t_p6 = pt_p6[i];
        const auto t_p7 = pt_p7[i];

        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];
        const auto u_p4 = pu_p4[i];
        const auto u_p5 = pu_p5[i];
        const auto u_p6 = pu_p6[i];

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4 - fr_p4 * r_2 * u_p4;

        ps_p5[i] = fz_p5 * z * t_p5 - fr_p5 * r_2 * u_p5;

        ps_p6[i] = fz_p6 * z * t_p6 - fr_p6 * r_2 * u_p6;

        ps_p7[i] = fz_p7 * z * t_p7;
    }

    return matrix;
}

auto
make_m_solid_harmonics(const CSimdMatrix &l_harmonics, const CSimdMatrix &k_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(l_harmonics.number_of_rows() == 17,
                              std::string("SimdHarmonics.make_m_solid_harmonics: Harmonics of angular momentum eight must have 17 rows"));

    errors::assertMsgCritical(k_harmonics.number_of_rows() == 15,
                              std::string("SimdHarmonics.make_m_solid_harmonics: Harmonics of angular momentum seven must have 15 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_m_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_m_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(l_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_m_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(19, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = (1.0 / 6.0) * std::sqrt(34.0);
    const auto fz_m8 = std::sqrt(17.0);
    const auto fz_m7 = 2.125 * std::sqrt(2.0);
    const auto fr_m7 = 0.125 * std::sqrt(30.0);
    const auto fz_m6 = (17.0 / 15.0) * std::sqrt(5.0);
    const auto fr_m6 = (2.0 / 15.0) * std::sqrt(35.0);
    const auto fz_m5 = (17.0 / 28.0) * std::sqrt(14.0);
    const auto fr_m5 = (1.0 / 28.0) * std::sqrt(546.0);
    const auto fz_m4 = (17.0 / 65.0) * std::sqrt(65.0);
    const auto fr_m4 = (4.0 / 65.0) * std::sqrt(195.0);
    const auto fz_m3 = (17.0 / 12.0) * std::sqrt(2.0);
    const auto fr_m3 = (1.0 / 12.0) * std::sqrt(110.0);
    const auto fz_m2 = (17.0 / 77.0) * std::sqrt(77.0);
    const auto fr_m2 = (2.0 / 77.0) * std::sqrt(1155.0);
    const auto fz_m1 = (17.0 / 20.0) * std::sqrt(5.0);
    const auto fr_m1 = (3.0 / 20.0) * std::sqrt(35.0);
    const auto fz_0 = 17.0 / 9.0;
    const auto fr_0 = 8.0 / 9.0;
    const auto fz_p1 = (17.0 / 20.0) * std::sqrt(5.0);
    const auto fr_p1 = (3.0 / 20.0) * std::sqrt(35.0);
    const auto fz_p2 = (17.0 / 77.0) * std::sqrt(77.0);
    const auto fr_p2 = (2.0 / 77.0) * std::sqrt(1155.0);
    const auto fz_p3 = (17.0 / 12.0) * std::sqrt(2.0);
    const auto fr_p3 = (1.0 / 12.0) * std::sqrt(110.0);
    const auto fz_p4 = (17.0 / 65.0) * std::sqrt(65.0);
    const auto fr_p4 = (4.0 / 65.0) * std::sqrt(195.0);
    const auto fz_p5 = (17.0 / 28.0) * std::sqrt(14.0);
    const auto fr_p5 = (1.0 / 28.0) * std::sqrt(546.0);
    const auto fz_p6 = (17.0 / 15.0) * std::sqrt(5.0);
    const auto fr_p6 = (2.0 / 15.0) * std::sqrt(35.0);
    const auto fz_p7 = 2.125 * std::sqrt(2.0);
    const auto fr_p7 = 0.125 * std::sqrt(30.0);
    const auto fz_p8 = std::sqrt(17.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m8 = l_harmonics.data(0);
    const auto *pt_m7 = l_harmonics.data(1);
    const auto *pt_m6 = l_harmonics.data(2);
    const auto *pt_m5 = l_harmonics.data(3);
    const auto *pt_m4 = l_harmonics.data(4);
    const auto *pt_m3 = l_harmonics.data(5);
    const auto *pt_m2 = l_harmonics.data(6);
    const auto *pt_m1 = l_harmonics.data(7);
    const auto *pt_0 = l_harmonics.data(8);
    const auto *pt_p1 = l_harmonics.data(9);
    const auto *pt_p2 = l_harmonics.data(10);
    const auto *pt_p3 = l_harmonics.data(11);
    const auto *pt_p4 = l_harmonics.data(12);
    const auto *pt_p5 = l_harmonics.data(13);
    const auto *pt_p6 = l_harmonics.data(14);
    const auto *pt_p7 = l_harmonics.data(15);
    const auto *pt_p8 = l_harmonics.data(16);

    const auto *pu_m7 = k_harmonics.data(0);
    const auto *pu_m6 = k_harmonics.data(1);
    const auto *pu_m5 = k_harmonics.data(2);
    const auto *pu_m4 = k_harmonics.data(3);
    const auto *pu_m3 = k_harmonics.data(4);
    const auto *pu_m2 = k_harmonics.data(5);
    const auto *pu_m1 = k_harmonics.data(6);
    const auto *pu_0 = k_harmonics.data(7);
    const auto *pu_p1 = k_harmonics.data(8);
    const auto *pu_p2 = k_harmonics.data(9);
    const auto *pu_p3 = k_harmonics.data(10);
    const auto *pu_p4 = k_harmonics.data(11);
    const auto *pu_p5 = k_harmonics.data(12);
    const auto *pu_p6 = k_harmonics.data(13);
    const auto *pu_p7 = k_harmonics.data(14);

    auto *ps_m9 = matrix.data(0);
    auto *ps_m8 = matrix.data(1);
    auto *ps_m7 = matrix.data(2);
    auto *ps_m6 = matrix.data(3);
    auto *ps_m5 = matrix.data(4);
    auto *ps_m4 = matrix.data(5);
    auto *ps_m3 = matrix.data(6);
    auto *ps_m2 = matrix.data(7);
    auto *ps_m1 = matrix.data(8);
    auto *ps_0 = matrix.data(9);
    auto *ps_p1 = matrix.data(10);
    auto *ps_p2 = matrix.data(11);
    auto *ps_p3 = matrix.data(12);
    auto *ps_p4 = matrix.data(13);
    auto *ps_p5 = matrix.data(14);
    auto *ps_p6 = matrix.data(15);
    auto *ps_p7 = matrix.data(16);
    auto *ps_p8 = matrix.data(17);
    auto *ps_p9 = matrix.data(18);

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 19 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m8, pt_m7, pt_m6, pt_m5, pt_m4, pt_m3, pt_m2, pt_m1, pt_0, pt_p8, pu_m7, pu_m6, pu_m5, pu_m4, pu_m3, pu_m2, pu_m1, pu_0, ps_m9, ps_m8, ps_m7, ps_m6, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p9 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m8 = pt_m8[i];
        const auto t_m7 = pt_m7[i];
        const auto t_m6 = pt_m6[i];
        const auto t_m5 = pt_m5[i];
        const auto t_m4 = pt_m4[i];
        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p8 = pt_p8[i];

        const auto u_m7 = pu_m7[i];
        const auto u_m6 = pu_m6[i];
        const auto u_m5 = pu_m5[i];
        const auto u_m4 = pu_m4[i];
        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];

        ps_p9[i] = fact * (x * t_p8 - y * t_m8);

        ps_m9[i] = fact * (y * t_p8 + x * t_m8);

        ps_m8[i] = fz_m8 * z * t_m8;

        ps_m7[i] = fz_m7 * z * t_m7 - fr_m7 * r_2 * u_m7;

        ps_m6[i] = fz_m6 * z * t_m6 - fr_m6 * r_2 * u_m6;

        ps_m5[i] = fz_m5 * z * t_m5 - fr_m5 * r_2 * u_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;
    }

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 19 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_p1, pt_p2, pt_p3, pt_p4, pt_p5, pt_p6, pt_p7, pt_p8, pu_p1, pu_p2, pu_p3, pu_p4, pu_p5, pu_p6, pu_p7, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5, ps_p6, ps_p7, ps_p8 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];
        const auto t_p5 = pt_p5[i];
        const auto t_p6 = pt_p6[i];
        const auto t_p7 = pt_p7[i];
        const auto t_p8 = pt_p8[i];

        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];
        const auto u_p4 = pu_p4[i];
        const auto u_p5 = pu_p5[i];
        const auto u_p6 = pu_p6[i];
        const auto u_p7 = pu_p7[i];

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4 - fr_p4 * r_2 * u_p4;

        ps_p5[i] = fz_p5 * z * t_p5 - fr_p5 * r_2 * u_p5;

        ps_p6[i] = fz_p6 * z * t_p6 - fr_p6 * r_2 * u_p6;

        ps_p7[i] = fz_p7 * z * t_p7 - fr_p7 * r_2 * u_p7;

        ps_p8[i] = fz_p8 * z * t_p8;
    }

    return matrix;
}

auto
make_n_solid_harmonics(const CSimdMatrix &m_harmonics, const CSimdMatrix &l_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(m_harmonics.number_of_rows() == 19,
                              std::string("SimdHarmonics.make_n_solid_harmonics: Harmonics of angular momentum nine must have 19 rows"));

    errors::assertMsgCritical(l_harmonics.number_of_rows() == 17,
                              std::string("SimdHarmonics.make_n_solid_harmonics: Harmonics of angular momentum eight must have 17 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_n_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_n_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(m_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_n_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(21, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = (1.0 / 10.0) * std::sqrt(95.0);
    const auto fz_m9 = std::sqrt(19.0);
    const auto fz_m8 = 19.0 / 6.0;
    const auto fr_m8 = (1.0 / 6.0) * std::sqrt(17.0);
    const auto fz_m7 = (19.0 / 51.0) * std::sqrt(51.0);
    const auto fr_m7 = (4.0 / 51.0) * std::sqrt(102.0);
    const auto fz_m6 = 2.375;
    const auto fr_m6 = 0.375 * std::sqrt(5.0);
    const auto fz_m5 = (19.0 / 15.0) * std::sqrt(3.0);
    const auto fr_m5 = (2.0 / 15.0) * std::sqrt(42.0);
    const auto fz_m4 = (19.0 / 42.0) * std::sqrt(21.0);
    const auto fr_m4 = (1.0 / 42.0) * std::sqrt(1365.0);
    const auto fz_m3 = (19.0 / 91.0) * std::sqrt(91.0);
    const auto fr_m3 = (6.0 / 91.0) * std::sqrt(182.0);
    const auto fz_m2 = (19.0 / 24.0) * std::sqrt(6.0);
    const auto fr_m2 = (1.0 / 24.0) * std::sqrt(462.0);
    const auto fz_m1 = (19.0 / 33.0) * std::sqrt(11.0);
    const auto fr_m1 = (4.0 / 33.0) * std::sqrt(55.0);
    const auto fz_0 = 19.0 / 10.0;
    const auto fr_0 = 9.0 / 10.0;
    const auto fz_p1 = (19.0 / 33.0) * std::sqrt(11.0);
    const auto fr_p1 = (4.0 / 33.0) * std::sqrt(55.0);
    const auto fz_p2 = (19.0 / 24.0) * std::sqrt(6.0);
    const auto fr_p2 = (1.0 / 24.0) * std::sqrt(462.0);
    const auto fz_p3 = (19.0 / 91.0) * std::sqrt(91.0);
    const auto fr_p3 = (6.0 / 91.0) * std::sqrt(182.0);
    const auto fz_p4 = (19.0 / 42.0) * std::sqrt(21.0);
    const auto fr_p4 = (1.0 / 42.0) * std::sqrt(1365.0);
    const auto fz_p5 = (19.0 / 15.0) * std::sqrt(3.0);
    const auto fr_p5 = (2.0 / 15.0) * std::sqrt(42.0);
    const auto fz_p6 = 2.375;
    const auto fr_p6 = 0.375 * std::sqrt(5.0);
    const auto fz_p7 = (19.0 / 51.0) * std::sqrt(51.0);
    const auto fr_p7 = (4.0 / 51.0) * std::sqrt(102.0);
    const auto fz_p8 = 19.0 / 6.0;
    const auto fr_p8 = (1.0 / 6.0) * std::sqrt(17.0);
    const auto fz_p9 = std::sqrt(19.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m9 = m_harmonics.data(0);
    const auto *pt_m8 = m_harmonics.data(1);
    const auto *pt_m7 = m_harmonics.data(2);
    const auto *pt_m6 = m_harmonics.data(3);
    const auto *pt_m5 = m_harmonics.data(4);
    const auto *pt_m4 = m_harmonics.data(5);
    const auto *pt_m3 = m_harmonics.data(6);
    const auto *pt_m2 = m_harmonics.data(7);
    const auto *pt_m1 = m_harmonics.data(8);
    const auto *pt_0 = m_harmonics.data(9);
    const auto *pt_p1 = m_harmonics.data(10);
    const auto *pt_p2 = m_harmonics.data(11);
    const auto *pt_p3 = m_harmonics.data(12);
    const auto *pt_p4 = m_harmonics.data(13);
    const auto *pt_p5 = m_harmonics.data(14);
    const auto *pt_p6 = m_harmonics.data(15);
    const auto *pt_p7 = m_harmonics.data(16);
    const auto *pt_p8 = m_harmonics.data(17);
    const auto *pt_p9 = m_harmonics.data(18);

    const auto *pu_m8 = l_harmonics.data(0);
    const auto *pu_m7 = l_harmonics.data(1);
    const auto *pu_m6 = l_harmonics.data(2);
    const auto *pu_m5 = l_harmonics.data(3);
    const auto *pu_m4 = l_harmonics.data(4);
    const auto *pu_m3 = l_harmonics.data(5);
    const auto *pu_m2 = l_harmonics.data(6);
    const auto *pu_m1 = l_harmonics.data(7);
    const auto *pu_0 = l_harmonics.data(8);
    const auto *pu_p1 = l_harmonics.data(9);
    const auto *pu_p2 = l_harmonics.data(10);
    const auto *pu_p3 = l_harmonics.data(11);
    const auto *pu_p4 = l_harmonics.data(12);
    const auto *pu_p5 = l_harmonics.data(13);
    const auto *pu_p6 = l_harmonics.data(14);
    const auto *pu_p7 = l_harmonics.data(15);
    const auto *pu_p8 = l_harmonics.data(16);

    auto *ps_m10 = matrix.data(0);
    auto *ps_m9 = matrix.data(1);
    auto *ps_m8 = matrix.data(2);
    auto *ps_m7 = matrix.data(3);
    auto *ps_m6 = matrix.data(4);
    auto *ps_m5 = matrix.data(5);
    auto *ps_m4 = matrix.data(6);
    auto *ps_m3 = matrix.data(7);
    auto *ps_m2 = matrix.data(8);
    auto *ps_m1 = matrix.data(9);
    auto *ps_0 = matrix.data(10);
    auto *ps_p1 = matrix.data(11);
    auto *ps_p2 = matrix.data(12);
    auto *ps_p3 = matrix.data(13);
    auto *ps_p4 = matrix.data(14);
    auto *ps_p5 = matrix.data(15);
    auto *ps_p6 = matrix.data(16);
    auto *ps_p7 = matrix.data(17);
    auto *ps_p8 = matrix.data(18);
    auto *ps_p9 = matrix.data(19);
    auto *ps_p10 = matrix.data(20);

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 21 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m9, pt_m8, pt_m7, pt_m6, pt_m5, pt_m4, pt_m3, pt_m2, pt_m1, pt_0, pt_p9, pu_m8, pu_m7, pu_m6, pu_m5, pu_m4, pu_m3, pu_m2, pu_m1, pu_0, ps_m10, ps_m9, ps_m8, ps_m7, ps_m6, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p10 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m9 = pt_m9[i];
        const auto t_m8 = pt_m8[i];
        const auto t_m7 = pt_m7[i];
        const auto t_m6 = pt_m6[i];
        const auto t_m5 = pt_m5[i];
        const auto t_m4 = pt_m4[i];
        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p9 = pt_p9[i];

        const auto u_m8 = pu_m8[i];
        const auto u_m7 = pu_m7[i];
        const auto u_m6 = pu_m6[i];
        const auto u_m5 = pu_m5[i];
        const auto u_m4 = pu_m4[i];
        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];

        ps_p10[i] = fact * (x * t_p9 - y * t_m9);

        ps_m10[i] = fact * (y * t_p9 + x * t_m9);

        ps_m9[i] = fz_m9 * z * t_m9;

        ps_m8[i] = fz_m8 * z * t_m8 - fr_m8 * r_2 * u_m8;

        ps_m7[i] = fz_m7 * z * t_m7 - fr_m7 * r_2 * u_m7;

        ps_m6[i] = fz_m6 * z * t_m6 - fr_m6 * r_2 * u_m6;

        ps_m5[i] = fz_m5 * z * t_m5 - fr_m5 * r_2 * u_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;
    }

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 21 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_p1, pt_p2, pt_p3, pt_p4, pt_p5, pt_p6, pt_p7, pt_p8, pt_p9, pu_p1, pu_p2, pu_p3, pu_p4, pu_p5, pu_p6, pu_p7, pu_p8, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5, ps_p6, ps_p7, ps_p8, ps_p9 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];
        const auto t_p5 = pt_p5[i];
        const auto t_p6 = pt_p6[i];
        const auto t_p7 = pt_p7[i];
        const auto t_p8 = pt_p8[i];
        const auto t_p9 = pt_p9[i];

        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];
        const auto u_p4 = pu_p4[i];
        const auto u_p5 = pu_p5[i];
        const auto u_p6 = pu_p6[i];
        const auto u_p7 = pu_p7[i];
        const auto u_p8 = pu_p8[i];

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4 - fr_p4 * r_2 * u_p4;

        ps_p5[i] = fz_p5 * z * t_p5 - fr_p5 * r_2 * u_p5;

        ps_p6[i] = fz_p6 * z * t_p6 - fr_p6 * r_2 * u_p6;

        ps_p7[i] = fz_p7 * z * t_p7 - fr_p7 * r_2 * u_p7;

        ps_p8[i] = fz_p8 * z * t_p8 - fr_p8 * r_2 * u_p8;

        ps_p9[i] = fz_p9 * z * t_p9;
    }

    return matrix;
}

auto
make_o_solid_harmonics(const CSimdMatrix &n_harmonics, const CSimdMatrix &m_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(n_harmonics.number_of_rows() == 21,
                              std::string("SimdHarmonics.make_o_solid_harmonics: Harmonics of angular momentum ten must have 21 rows"));

    errors::assertMsgCritical(m_harmonics.number_of_rows() == 19,
                              std::string("SimdHarmonics.make_o_solid_harmonics: Harmonics of angular momentum nine must have 19 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_o_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_o_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(n_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_o_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(23, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = (1.0 / 22.0) * std::sqrt(462.0);
    const auto fz_m10 = std::sqrt(21.0);
    const auto fz_m9 = (21.0 / 20.0) * std::sqrt(10.0);
    const auto fr_m9 = (1.0 / 20.0) * std::sqrt(190.0);
    const auto fz_m8 = (7.0 / 19.0) * std::sqrt(57.0);
    const auto fr_m8 = (2.0 / 19.0) * std::sqrt(57.0);
    const auto fz_m7 = 1.75 * std::sqrt(2.0);
    const auto fr_m7 = (1.0 / 12.0) * std::sqrt(102.0);
    const auto fz_m6 = (21.0 / 85.0) * std::sqrt(85.0);
    const auto fr_m6 = (8.0 / 85.0) * std::sqrt(85.0);
    const auto fz_m5 = 0.875 * std::sqrt(6.0);
    const auto fr_m5 = 0.625 * std::sqrt(2.0);
    const auto fz_m4 = (1.0 / 5.0) * std::sqrt(105.0);
    const auto fr_m4 = (2.0 / 5.0) * std::sqrt(5.0);
    const auto fz_m3 = 0.75 * std::sqrt(7.0);
    const auto fr_m3 = 0.25 * std::sqrt(13.0);
    const auto fz_m2 = (7.0 / 13.0) * std::sqrt(13.0);
    const auto fr_m2 = (4.0 / 39.0) * std::sqrt(78.0);
    const auto fz_m1 = (7.0 / 20.0) * std::sqrt(30.0);
    const auto fr_m1 = (1.0 / 20.0) * std::sqrt(330.0);
    const auto fz_0 = 21.0 / 11.0;
    const auto fr_0 = 10.0 / 11.0;
    const auto fz_p1 = (7.0 / 20.0) * std::sqrt(30.0);
    const auto fr_p1 = (1.0 / 20.0) * std::sqrt(330.0);
    const auto fz_p2 = (7.0 / 13.0) * std::sqrt(13.0);
    const auto fr_p2 = (4.0 / 39.0) * std::sqrt(78.0);
    const auto fz_p3 = 0.75 * std::sqrt(7.0);
    const auto fr_p3 = 0.25 * std::sqrt(13.0);
    const auto fz_p4 = (1.0 / 5.0) * std::sqrt(105.0);
    const auto fr_p4 = (2.0 / 5.0) * std::sqrt(5.0);
    const auto fz_p5 = 0.875 * std::sqrt(6.0);
    const auto fr_p5 = 0.625 * std::sqrt(2.0);
    const auto fz_p6 = (21.0 / 85.0) * std::sqrt(85.0);
    const auto fr_p6 = (8.0 / 85.0) * std::sqrt(85.0);
    const auto fz_p7 = 1.75 * std::sqrt(2.0);
    const auto fr_p7 = (1.0 / 12.0) * std::sqrt(102.0);
    const auto fz_p8 = (7.0 / 19.0) * std::sqrt(57.0);
    const auto fr_p8 = (2.0 / 19.0) * std::sqrt(57.0);
    const auto fz_p9 = (21.0 / 20.0) * std::sqrt(10.0);
    const auto fr_p9 = (1.0 / 20.0) * std::sqrt(190.0);
    const auto fz_p10 = std::sqrt(21.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m10 = n_harmonics.data(0);
    const auto *pt_m9 = n_harmonics.data(1);
    const auto *pt_m8 = n_harmonics.data(2);
    const auto *pt_m7 = n_harmonics.data(3);
    const auto *pt_m6 = n_harmonics.data(4);
    const auto *pt_m5 = n_harmonics.data(5);
    const auto *pt_m4 = n_harmonics.data(6);
    const auto *pt_m3 = n_harmonics.data(7);
    const auto *pt_m2 = n_harmonics.data(8);
    const auto *pt_m1 = n_harmonics.data(9);
    const auto *pt_0 = n_harmonics.data(10);
    const auto *pt_p1 = n_harmonics.data(11);
    const auto *pt_p2 = n_harmonics.data(12);
    const auto *pt_p3 = n_harmonics.data(13);
    const auto *pt_p4 = n_harmonics.data(14);
    const auto *pt_p5 = n_harmonics.data(15);
    const auto *pt_p6 = n_harmonics.data(16);
    const auto *pt_p7 = n_harmonics.data(17);
    const auto *pt_p8 = n_harmonics.data(18);
    const auto *pt_p9 = n_harmonics.data(19);
    const auto *pt_p10 = n_harmonics.data(20);

    const auto *pu_m9 = m_harmonics.data(0);
    const auto *pu_m8 = m_harmonics.data(1);
    const auto *pu_m7 = m_harmonics.data(2);
    const auto *pu_m6 = m_harmonics.data(3);
    const auto *pu_m5 = m_harmonics.data(4);
    const auto *pu_m4 = m_harmonics.data(5);
    const auto *pu_m3 = m_harmonics.data(6);
    const auto *pu_m2 = m_harmonics.data(7);
    const auto *pu_m1 = m_harmonics.data(8);
    const auto *pu_0 = m_harmonics.data(9);
    const auto *pu_p1 = m_harmonics.data(10);
    const auto *pu_p2 = m_harmonics.data(11);
    const auto *pu_p3 = m_harmonics.data(12);
    const auto *pu_p4 = m_harmonics.data(13);
    const auto *pu_p5 = m_harmonics.data(14);
    const auto *pu_p6 = m_harmonics.data(15);
    const auto *pu_p7 = m_harmonics.data(16);
    const auto *pu_p8 = m_harmonics.data(17);
    const auto *pu_p9 = m_harmonics.data(18);

    auto *ps_m11 = matrix.data(0);
    auto *ps_m10 = matrix.data(1);
    auto *ps_m9 = matrix.data(2);
    auto *ps_m8 = matrix.data(3);
    auto *ps_m7 = matrix.data(4);
    auto *ps_m6 = matrix.data(5);
    auto *ps_m5 = matrix.data(6);
    auto *ps_m4 = matrix.data(7);
    auto *ps_m3 = matrix.data(8);
    auto *ps_m2 = matrix.data(9);
    auto *ps_m1 = matrix.data(10);
    auto *ps_0 = matrix.data(11);
    auto *ps_p1 = matrix.data(12);
    auto *ps_p2 = matrix.data(13);
    auto *ps_p3 = matrix.data(14);
    auto *ps_p4 = matrix.data(15);
    auto *ps_p5 = matrix.data(16);
    auto *ps_p6 = matrix.data(17);
    auto *ps_p7 = matrix.data(18);
    auto *ps_p8 = matrix.data(19);
    auto *ps_p9 = matrix.data(20);
    auto *ps_p10 = matrix.data(21);
    auto *ps_p11 = matrix.data(22);

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 23 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m10, pt_m9, pt_m8, pt_m7, pt_m6, pt_m5, pt_m4, pt_m3, pt_m2, pt_m1, pt_0, pt_p10, pu_m9, pu_m8, pu_m7, pu_m6, pu_m5, pu_m4, pu_m3, pu_m2, pu_m1, pu_0, ps_m11, ps_m10, ps_m9, ps_m8, ps_m7, ps_m6, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p11 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m10 = pt_m10[i];
        const auto t_m9 = pt_m9[i];
        const auto t_m8 = pt_m8[i];
        const auto t_m7 = pt_m7[i];
        const auto t_m6 = pt_m6[i];
        const auto t_m5 = pt_m5[i];
        const auto t_m4 = pt_m4[i];
        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p10 = pt_p10[i];

        const auto u_m9 = pu_m9[i];
        const auto u_m8 = pu_m8[i];
        const auto u_m7 = pu_m7[i];
        const auto u_m6 = pu_m6[i];
        const auto u_m5 = pu_m5[i];
        const auto u_m4 = pu_m4[i];
        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];

        ps_p11[i] = fact * (x * t_p10 - y * t_m10);

        ps_m11[i] = fact * (y * t_p10 + x * t_m10);

        ps_m10[i] = fz_m10 * z * t_m10;

        ps_m9[i] = fz_m9 * z * t_m9 - fr_m9 * r_2 * u_m9;

        ps_m8[i] = fz_m8 * z * t_m8 - fr_m8 * r_2 * u_m8;

        ps_m7[i] = fz_m7 * z * t_m7 - fr_m7 * r_2 * u_m7;

        ps_m6[i] = fz_m6 * z * t_m6 - fr_m6 * r_2 * u_m6;

        ps_m5[i] = fz_m5 * z * t_m5 - fr_m5 * r_2 * u_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;
    }

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 23 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_p1, pt_p2, pt_p3, pt_p4, pt_p5, pt_p6, pt_p7, pt_p8, pt_p9, pt_p10, pu_p1, pu_p2, pu_p3, pu_p4, pu_p5, pu_p6, pu_p7, pu_p8, pu_p9, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5, ps_p6, ps_p7, ps_p8, ps_p9, ps_p10 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];
        const auto t_p5 = pt_p5[i];
        const auto t_p6 = pt_p6[i];
        const auto t_p7 = pt_p7[i];
        const auto t_p8 = pt_p8[i];
        const auto t_p9 = pt_p9[i];
        const auto t_p10 = pt_p10[i];

        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];
        const auto u_p4 = pu_p4[i];
        const auto u_p5 = pu_p5[i];
        const auto u_p6 = pu_p6[i];
        const auto u_p7 = pu_p7[i];
        const auto u_p8 = pu_p8[i];
        const auto u_p9 = pu_p9[i];

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4 - fr_p4 * r_2 * u_p4;

        ps_p5[i] = fz_p5 * z * t_p5 - fr_p5 * r_2 * u_p5;

        ps_p6[i] = fz_p6 * z * t_p6 - fr_p6 * r_2 * u_p6;

        ps_p7[i] = fz_p7 * z * t_p7 - fr_p7 * r_2 * u_p7;

        ps_p8[i] = fz_p8 * z * t_p8 - fr_p8 * r_2 * u_p8;

        ps_p9[i] = fz_p9 * z * t_p9 - fr_p9 * r_2 * u_p9;

        ps_p10[i] = fz_p10 * z * t_p10;
    }

    return matrix;
}

auto
make_q_solid_harmonics(const CSimdMatrix &o_harmonics, const CSimdMatrix &n_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(o_harmonics.number_of_rows() == 23,
                              std::string("SimdHarmonics.make_q_solid_harmonics: Harmonics of angular momentum eleven must have 23 rows"));

    errors::assertMsgCritical(n_harmonics.number_of_rows() == 21,
                              std::string("SimdHarmonics.make_q_solid_harmonics: Harmonics of angular momentum ten must have 21 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_q_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_q_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(o_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_q_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(25, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = (1.0 / 12.0) * std::sqrt(138.0);
    const auto fz_m11 = std::sqrt(23.0);
    const auto fz_m10 = (23.0 / 22.0) * std::sqrt(11.0);
    const auto fr_m10 = (1.0 / 22.0) * std::sqrt(231.0);
    const auto fz_m9 = (23.0 / 21.0) * std::sqrt(7.0);
    const auto fr_m9 = (2.0 / 21.0) * std::sqrt(70.0);
    const auto fz_m8 = (23.0 / 20.0) * std::sqrt(5.0);
    const auto fr_m8 = (1.0 / 20.0) * std::sqrt(285.0);
    const auto fz_m7 = (23.0 / 95.0) * std::sqrt(95.0);
    const auto fr_m7 = (6.0 / 95.0) * std::sqrt(190.0);
    const auto fz_m6 = (23.0 / 18.0) * std::sqrt(3.0);
    const auto fr_m6 = (1.0 / 18.0) * std::sqrt(255.0);
    const auto fz_m5 = (23.0 / 119.0) * std::sqrt(119.0);
    const auto fr_m5 = (4.0 / 119.0) * std::sqrt(714.0);
    const auto fz_m4 = 1.4375 * std::sqrt(2.0);
    const auto fr_m4 = 0.0625 * std::sqrt(210.0);
    const auto fz_m3 = (23.0 / 45.0) * std::sqrt(15.0);
    const auto fr_m3 = (4.0 / 45.0) * std::sqrt(105.0);
    const auto fz_m2 = (23.0 / 70.0) * std::sqrt(35.0);
    const auto fr_m2 = (3.0 / 70.0) * std::sqrt(455.0);
    const auto fz_m1 = (23.0 / 143.0) * std::sqrt(143.0);
    const auto fr_m1 = (2.0 / 143.0) * std::sqrt(4290.0);
    const auto fz_0 = 23.0 / 12.0;
    const auto fr_0 = 11.0 / 12.0;
    const auto fz_p1 = (23.0 / 143.0) * std::sqrt(143.0);
    const auto fr_p1 = (2.0 / 143.0) * std::sqrt(4290.0);
    const auto fz_p2 = (23.0 / 70.0) * std::sqrt(35.0);
    const auto fr_p2 = (3.0 / 70.0) * std::sqrt(455.0);
    const auto fz_p3 = (23.0 / 45.0) * std::sqrt(15.0);
    const auto fr_p3 = (4.0 / 45.0) * std::sqrt(105.0);
    const auto fz_p4 = 1.4375 * std::sqrt(2.0);
    const auto fr_p4 = 0.0625 * std::sqrt(210.0);
    const auto fz_p5 = (23.0 / 119.0) * std::sqrt(119.0);
    const auto fr_p5 = (4.0 / 119.0) * std::sqrt(714.0);
    const auto fz_p6 = (23.0 / 18.0) * std::sqrt(3.0);
    const auto fr_p6 = (1.0 / 18.0) * std::sqrt(255.0);
    const auto fz_p7 = (23.0 / 95.0) * std::sqrt(95.0);
    const auto fr_p7 = (6.0 / 95.0) * std::sqrt(190.0);
    const auto fz_p8 = (23.0 / 20.0) * std::sqrt(5.0);
    const auto fr_p8 = (1.0 / 20.0) * std::sqrt(285.0);
    const auto fz_p9 = (23.0 / 21.0) * std::sqrt(7.0);
    const auto fr_p9 = (2.0 / 21.0) * std::sqrt(70.0);
    const auto fz_p10 = (23.0 / 22.0) * std::sqrt(11.0);
    const auto fr_p10 = (1.0 / 22.0) * std::sqrt(231.0);
    const auto fz_p11 = std::sqrt(23.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m11 = o_harmonics.data(0);
    const auto *pt_m10 = o_harmonics.data(1);
    const auto *pt_m9 = o_harmonics.data(2);
    const auto *pt_m8 = o_harmonics.data(3);
    const auto *pt_m7 = o_harmonics.data(4);
    const auto *pt_m6 = o_harmonics.data(5);
    const auto *pt_m5 = o_harmonics.data(6);
    const auto *pt_m4 = o_harmonics.data(7);
    const auto *pt_m3 = o_harmonics.data(8);
    const auto *pt_m2 = o_harmonics.data(9);
    const auto *pt_m1 = o_harmonics.data(10);
    const auto *pt_0 = o_harmonics.data(11);
    const auto *pt_p1 = o_harmonics.data(12);
    const auto *pt_p2 = o_harmonics.data(13);
    const auto *pt_p3 = o_harmonics.data(14);
    const auto *pt_p4 = o_harmonics.data(15);
    const auto *pt_p5 = o_harmonics.data(16);
    const auto *pt_p6 = o_harmonics.data(17);
    const auto *pt_p7 = o_harmonics.data(18);
    const auto *pt_p8 = o_harmonics.data(19);
    const auto *pt_p9 = o_harmonics.data(20);
    const auto *pt_p10 = o_harmonics.data(21);
    const auto *pt_p11 = o_harmonics.data(22);

    const auto *pu_m10 = n_harmonics.data(0);
    const auto *pu_m9 = n_harmonics.data(1);
    const auto *pu_m8 = n_harmonics.data(2);
    const auto *pu_m7 = n_harmonics.data(3);
    const auto *pu_m6 = n_harmonics.data(4);
    const auto *pu_m5 = n_harmonics.data(5);
    const auto *pu_m4 = n_harmonics.data(6);
    const auto *pu_m3 = n_harmonics.data(7);
    const auto *pu_m2 = n_harmonics.data(8);
    const auto *pu_m1 = n_harmonics.data(9);
    const auto *pu_0 = n_harmonics.data(10);
    const auto *pu_p1 = n_harmonics.data(11);
    const auto *pu_p2 = n_harmonics.data(12);
    const auto *pu_p3 = n_harmonics.data(13);
    const auto *pu_p4 = n_harmonics.data(14);
    const auto *pu_p5 = n_harmonics.data(15);
    const auto *pu_p6 = n_harmonics.data(16);
    const auto *pu_p7 = n_harmonics.data(17);
    const auto *pu_p8 = n_harmonics.data(18);
    const auto *pu_p9 = n_harmonics.data(19);
    const auto *pu_p10 = n_harmonics.data(20);

    auto *ps_m12 = matrix.data(0);
    auto *ps_m11 = matrix.data(1);
    auto *ps_m10 = matrix.data(2);
    auto *ps_m9 = matrix.data(3);
    auto *ps_m8 = matrix.data(4);
    auto *ps_m7 = matrix.data(5);
    auto *ps_m6 = matrix.data(6);
    auto *ps_m5 = matrix.data(7);
    auto *ps_m4 = matrix.data(8);
    auto *ps_m3 = matrix.data(9);
    auto *ps_m2 = matrix.data(10);
    auto *ps_m1 = matrix.data(11);
    auto *ps_0 = matrix.data(12);
    auto *ps_p1 = matrix.data(13);
    auto *ps_p2 = matrix.data(14);
    auto *ps_p3 = matrix.data(15);
    auto *ps_p4 = matrix.data(16);
    auto *ps_p5 = matrix.data(17);
    auto *ps_p6 = matrix.data(18);
    auto *ps_p7 = matrix.data(19);
    auto *ps_p8 = matrix.data(20);
    auto *ps_p9 = matrix.data(21);
    auto *ps_p10 = matrix.data(22);
    auto *ps_p11 = matrix.data(23);
    auto *ps_p12 = matrix.data(24);

    // NOTE: the rows are formed in three loops, as the vectorizer runs out of
    // registers with all 25 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m11, pt_m10, pt_m9, pt_m8, pt_m7, pt_m6, pt_m5, pt_m4, pt_p11, pu_m10, pu_m9, pu_m8, pu_m7, pu_m6, pu_m5, pu_m4, ps_m12, ps_m11, ps_m10, ps_m9, ps_m8, ps_m7, ps_m6, ps_m5, ps_m4, ps_p12 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m11 = pt_m11[i];
        const auto t_m10 = pt_m10[i];
        const auto t_m9 = pt_m9[i];
        const auto t_m8 = pt_m8[i];
        const auto t_m7 = pt_m7[i];
        const auto t_m6 = pt_m6[i];
        const auto t_m5 = pt_m5[i];
        const auto t_m4 = pt_m4[i];
        const auto t_p11 = pt_p11[i];

        const auto u_m10 = pu_m10[i];
        const auto u_m9 = pu_m9[i];
        const auto u_m8 = pu_m8[i];
        const auto u_m7 = pu_m7[i];
        const auto u_m6 = pu_m6[i];
        const auto u_m5 = pu_m5[i];
        const auto u_m4 = pu_m4[i];

        ps_p12[i] = fact * (x * t_p11 - y * t_m11);

        ps_m12[i] = fact * (y * t_p11 + x * t_m11);

        ps_m11[i] = fz_m11 * z * t_m11;

        ps_m10[i] = fz_m10 * z * t_m10 - fr_m10 * r_2 * u_m10;

        ps_m9[i] = fz_m9 * z * t_m9 - fr_m9 * r_2 * u_m9;

        ps_m8[i] = fz_m8 * z * t_m8 - fr_m8 * r_2 * u_m8;

        ps_m7[i] = fz_m7 * z * t_m7 - fr_m7 * r_2 * u_m7;

        ps_m6[i] = fz_m6 * z * t_m6 - fr_m6 * r_2 * u_m6;

        ps_m5[i] = fz_m5 * z * t_m5 - fr_m5 * r_2 * u_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;
    }

    // NOTE: the rows are formed in three loops, as the vectorizer runs out of
    // registers with all 25 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_m3, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, pt_p3, pt_p4, pu_m3, pu_m2, pu_m1, pu_0, pu_p1, pu_p2, pu_p3, pu_p4, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3, ps_p4 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];

        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];
        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];
        const auto u_p4 = pu_p4[i];

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4 - fr_p4 * r_2 * u_p4;
    }

    // NOTE: the rows are formed in three loops, as the vectorizer runs out of
    // registers with all 25 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_p5, pt_p6, pt_p7, pt_p8, pt_p9, pt_p10, pt_p11, pu_p5, pu_p6, pu_p7, pu_p8, pu_p9, pu_p10, ps_p5, ps_p6, ps_p7, ps_p8, ps_p9, ps_p10, ps_p11 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_p5 = pt_p5[i];
        const auto t_p6 = pt_p6[i];
        const auto t_p7 = pt_p7[i];
        const auto t_p8 = pt_p8[i];
        const auto t_p9 = pt_p9[i];
        const auto t_p10 = pt_p10[i];
        const auto t_p11 = pt_p11[i];

        const auto u_p5 = pu_p5[i];
        const auto u_p6 = pu_p6[i];
        const auto u_p7 = pu_p7[i];
        const auto u_p8 = pu_p8[i];
        const auto u_p9 = pu_p9[i];
        const auto u_p10 = pu_p10[i];

        ps_p5[i] = fz_p5 * z * t_p5 - fr_p5 * r_2 * u_p5;

        ps_p6[i] = fz_p6 * z * t_p6 - fr_p6 * r_2 * u_p6;

        ps_p7[i] = fz_p7 * z * t_p7 - fr_p7 * r_2 * u_p7;

        ps_p8[i] = fz_p8 * z * t_p8 - fr_p8 * r_2 * u_p8;

        ps_p9[i] = fz_p9 * z * t_p9 - fr_p9 * r_2 * u_p9;

        ps_p10[i] = fz_p10 * z * t_p10 - fr_p10 * r_2 * u_p10;

        ps_p11[i] = fz_p11 * z * t_p11;
    }

    return matrix;
}

}  // namespace simdfunc
