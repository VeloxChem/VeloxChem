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


#include "SimdOverlapCtrVrrFS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ctr_fs_overlap_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                         const size_t pa, const size_t ps, const size_t ds, const size_t ncols,
                         const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.25 * std::sqrt(10.0) / p;
    const auto f_1 = 0.75 * std::sqrt(10.0);
    const auto f_2 = 0.25 * std::sqrt(10.0);
    const auto f_3 = std::sqrt(15.0);
    const auto f_4 = 0.25 * std::sqrt(6.0) / p;
    const auto f_5 = 0.25 * std::sqrt(6.0);
    const auto f_6 = std::sqrt(6.0);
    const auto f_7 = 1.0 / p;
    const auto f_8 = 0.5 * std::sqrt(15.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

    const auto *ds_0 = buffer.data(ds + 0);
    const auto *ds_1 = buffer.data(ds + 1);
    const auto *ds_2 = buffer.data(ds + 2);
    const auto *ds_3 = buffer.data(ds + 3);

#pragma omp simd aligned(pa_x, pa_y, pa_z, ps_0, ps_1, ps_2, ds_0, ds_1, ds_2, \
                         ds_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += -f_0 * ps_1[k]
                  + f_1 * pa_y[k] * ds_0[k]
                  - f_2 * pa_y[k] * ds_1[k];

        g_1[k] += f_3 * pa_x[k] * ds_2[k];

        g_2[k] += -f_4 * ps_1[k]
                  - f_5 * pa_y[k] * ds_0[k]
                  - f_5 * pa_y[k] * ds_1[k]
                  + f_6 * pa_y[k] * ds_3[k];

        g_3[k] += f_7 * ps_2[k]
                  - 1.5 * pa_z[k] * ds_0[k]
                  - 1.5 * pa_z[k] * ds_1[k]
                  + pa_z[k] * ds_3[k];

        g_4[k] += -f_4 * ps_0[k]
                  - f_5 * pa_x[k] * ds_0[k]
                  - f_5 * pa_x[k] * ds_1[k]
                  + f_6 * pa_x[k] * ds_3[k];

        g_5[k] += f_8 * pa_z[k] * ds_0[k]
                  - f_8 * pa_z[k] * ds_1[k];
    }

#pragma omp simd aligned(pa_x, ps_0, ds_0, ds_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_6[k] += f_0 * ps_0[k]
                  + f_2 * pa_x[k] * ds_0[k]
                  - f_1 * pa_x[k] * ds_1[k];
    }
}

}  // namespace simdovl
