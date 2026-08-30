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



#include "SimdHarmonicsRecF.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>

#include "ErrorHandler.hpp"
#include "OpenMPFunc.hpp"
#include "SimdAlign.hpp"

namespace simdfunc {  // simdfunc namespace

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

    const auto fact = std::sqrt(5.0 / 6.0);
    const auto fz_m2 = std::sqrt(5.0);
    const auto fz_m1 = std::sqrt(3.125);
    const auto fr_m1 = std::sqrt(0.375);
    const auto fz_0 = 5.0 / 3.0;
    const auto fr_0 = 2.0 / 3.0;
    const auto fz_p1 = std::sqrt(3.125);
    const auto fr_p1 = std::sqrt(0.375);
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

    // NOTE: the rows are filled by the threads only when they hold enough work
    // to repay the fixed cost of forking and joining them. The two loops are
    // written out, as a loop under a pragma with an if clause does not vectorize
    // as well as one without it when the clause is false.

    const auto npnts = static_cast<int64_t>(7 * ncols);

    // NOTE: the work needed to repay the threads is taken per thread, as forking
    // and joining them costs more the more of them there are. The constant is
    // measured on fourteen threads, where the threshold it gives is the smallest
    // one which does not lose on the blocks too small to fill them.

    const auto nthreshold = int64_t{12000} * omp::get_number_of_threads();

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

    if (npnts > nthreshold)
    {
#pragma omp parallel for simd schedule(static) aligned(pr_x, pr_y, pr_z, pr_2, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3 : simd::cache_line_size())
        for (int64_t i = 0; i < static_cast<int64_t>(ncols); i++)
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
    }
    else
    {
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
    }

    return matrix;
}

}  // namespace simdfunc
