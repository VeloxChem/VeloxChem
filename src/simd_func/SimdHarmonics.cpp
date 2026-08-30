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

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

namespace simdfunc {  // simdfunc namespace

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
                              std::string("SimdHarmonics.make_d_solid_harmonics: Harmonics of angular momentum 1 must have 3 rows"));

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
                              std::string("SimdHarmonics.make_f_solid_harmonics: Harmonics of angular momentum 2 must have 5 rows"));

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
    const auto fz_0 = 1.6666666666666667;
    const auto fr_0 = 0.6666666666666666;
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
                              std::string("SimdHarmonics.make_g_solid_harmonics: Harmonics of angular momentum 3 must have 7 rows"));

    errors::assertMsgCritical(d_harmonics.number_of_rows() == 5,
                              std::string("SimdHarmonics.make_g_solid_harmonics: Harmonics of angular momentum 2 must have 5 rows"));

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
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics of angular momentum 4 must have 9 rows"));

    errors::assertMsgCritical(f_harmonics.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics of angular momentum 3 must have 7 rows"));

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
    const auto fz_0 = 1.8;
    const auto fr_0 = 0.8;
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
                              std::string("SimdHarmonics.make_i_solid_harmonics: Harmonics of angular momentum 5 must have 11 rows"));

    errors::assertMsgCritical(g_harmonics.number_of_rows() == 9,
                              std::string("SimdHarmonics.make_i_solid_harmonics: Harmonics of angular momentum 4 must have 9 rows"));

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
    const auto fz_0 = 1.8333333333333333;
    const auto fr_0 = 0.8333333333333334;
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
                              std::string("SimdHarmonics.make_k_solid_harmonics: Harmonics of angular momentum 6 must have 13 rows"));

    errors::assertMsgCritical(h_harmonics.number_of_rows() == 11,
                              std::string("SimdHarmonics.make_k_solid_harmonics: Harmonics of angular momentum 5 must have 11 rows"));

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
    const auto fz_0 = 1.8571428571428572;
    const auto fr_0 = 0.8571428571428571;
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
                              std::string("SimdHarmonics.make_l_solid_harmonics: Harmonics of angular momentum 7 must have 15 rows"));

    errors::assertMsgCritical(i_harmonics.number_of_rows() == 13,
                              std::string("SimdHarmonics.make_l_solid_harmonics: Harmonics of angular momentum 6 must have 13 rows"));

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
    // distance are loaded by both loops, every other value once.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m7, pt_m6, pt_m5, pt_m4, pt_m3, pt_m2, pt_m1, pt_p7, pu_m6, pu_m5, pu_m4, pu_m3, pu_m2, pu_m1, ps_m8, ps_m7, ps_m6, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_p8 : simd::cache_line_size())
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
        const auto t_p7 = pt_p7[i];

        const auto u_m6 = pu_m6[i];
        const auto u_m5 = pu_m5[i];
        const auto u_m4 = pu_m4[i];
        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];

        ps_p8[i] = fact * (x * t_p7 - y * t_m7);

        ps_m8[i] = fact * (y * t_p7 + x * t_m7);

        ps_m7[i] = fz_m7 * z * t_m7;

        ps_m6[i] = fz_m6 * z * t_m6 - fr_m6 * r_2 * u_m6;

        ps_m5[i] = fz_m5 * z * t_m5 - fr_m5 * r_2 * u_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;
    }

    // NOTE: the rows are formed in two loops, as the vectorizer runs out of
    // registers with all 17 of them in one. Only the components and the squared
    // distance are loaded by both loops, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_0, pt_p1, pt_p2, pt_p3, pt_p4, pt_p5, pt_p6, pt_p7, pu_0, pu_p1, pu_p2, pu_p3, pu_p4, pu_p5, pu_p6, ps_0, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5, ps_p6, ps_p7 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_0 = pt_0[i];
        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];
        const auto t_p5 = pt_p5[i];
        const auto t_p6 = pt_p6[i];
        const auto t_p7 = pt_p7[i];

        const auto u_0 = pu_0[i];
        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];
        const auto u_p4 = pu_p4[i];
        const auto u_p5 = pu_p5[i];
        const auto u_p6 = pu_p6[i];

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;

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

}  // namespace simdfunc
