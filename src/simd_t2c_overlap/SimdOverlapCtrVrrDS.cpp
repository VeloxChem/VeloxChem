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


#include "SimdOverlapCtrVrrDS.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ctr_ds_overlap_0(double *values, const size_t nvalues, CSimdMatrix &buffer,
                         const size_t pa, const size_t ps, const size_t ncols) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = std::sqrt(3.0);
    const auto f_1 = 0.5 * std::sqrt(3.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ps_0 = buffer.data(ps + 0);
    const auto *ps_1 = buffer.data(ps + 1);
    const auto *ps_2 = buffer.data(ps + 2);

#pragma omp simd aligned(pa_x, pa_y, pa_z, ps_0, ps_1, ps_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        g_0[k] += f_0 * pa_x[k] * ps_1[k];

        g_1[k] += f_0 * pa_y[k] * ps_2[k];

        g_2[k] += -0.5 * pa_x[k] * ps_0[k]
                  - 0.5 * pa_y[k] * ps_1[k]
                  + pa_z[k] * ps_2[k];

        g_3[k] += f_0 * pa_x[k] * ps_2[k];

        g_4[k] += f_1 * pa_x[k] * ps_0[k]
                  - f_1 * pa_y[k] * ps_1[k];
    }
}

}  // namespace simdovl
