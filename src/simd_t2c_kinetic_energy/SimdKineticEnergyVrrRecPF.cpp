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


#include "SimdKineticEnergyVrrRecPF.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_pf_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t sd, const size_t sf,
                                 const size_t pf_s, const size_t pd, const size_t ncols,
                                 const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sd_0 = buffer.data(sd + 0);
    const auto *sd_3 = buffer.data(sd + 3);
    const auto *sd_5 = buffer.data(sd + 5);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_9 = buffer.data(sf + 9);

    const auto *pf_s_0 = buffer.data(pf_s + 0);
    const auto *pf_s_1 = buffer.data(pf_s + 1);
    const auto *pf_s_2 = buffer.data(pf_s + 2);
    const auto *pf_s_3 = buffer.data(pf_s + 3);
    const auto *pf_s_4 = buffer.data(pf_s + 4);
    const auto *pf_s_5 = buffer.data(pf_s + 5);
    const auto *pf_s_6 = buffer.data(pf_s + 6);
    const auto *pf_s_7 = buffer.data(pf_s + 7);
    const auto *pf_s_8 = buffer.data(pf_s + 8);
    const auto *pf_s_9 = buffer.data(pf_s + 9);
    const auto *pf_s_10 = buffer.data(pf_s + 10);
    const auto *pf_s_11 = buffer.data(pf_s + 11);
    const auto *pf_s_12 = buffer.data(pf_s + 12);
    const auto *pf_s_13 = buffer.data(pf_s + 13);
    const auto *pf_s_14 = buffer.data(pf_s + 14);
    const auto *pf_s_15 = buffer.data(pf_s + 15);
    const auto *pf_s_16 = buffer.data(pf_s + 16);
    const auto *pf_s_17 = buffer.data(pf_s + 17);
    const auto *pf_s_18 = buffer.data(pf_s + 18);
    const auto *pf_s_19 = buffer.data(pf_s + 19);
    const auto *pf_s_20 = buffer.data(pf_s + 20);
    const auto *pf_s_21 = buffer.data(pf_s + 21);
    const auto *pf_s_22 = buffer.data(pf_s + 22);
    const auto *pf_s_23 = buffer.data(pf_s + 23);
    const auto *pf_s_24 = buffer.data(pf_s + 24);
    const auto *pf_s_25 = buffer.data(pf_s + 25);
    const auto *pf_s_26 = buffer.data(pf_s + 26);
    const auto *pf_s_27 = buffer.data(pf_s + 27);
    const auto *pf_s_28 = buffer.data(pf_s + 28);
    const auto *pf_s_29 = buffer.data(pf_s + 29);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_10 = buffer.data(pd + 10);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_16 = buffer.data(pd + 16);
    const auto *pd_17 = buffer.data(pd + 17);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pb_y, pb_z, sd_0, sf_0, pf_s_0, pf_s_1, pf_s_2, \
                         pd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sd_0[k]
                 + pa_x[k] * sf_0[k]
                 + f_1 * pf_s_0[k];

        t_1[k] = f_1 * pf_s_1[k]
                 + pb_y[k] * pd_0[k];

        t_2[k] = f_1 * pf_s_2[k]
                 + pb_z[k] * pd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, sd_3, sd_5, pf_s_3, pf_s_4, pf_s_5, pd_2, \
                         pd_3, pd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_2 * sd_3[k]
                 + f_1 * pf_s_3[k]
                 + pb_x[k] * pd_3[k];

        t_4[k] = f_1 * pf_s_4[k]
                 + pb_y[k] * pd_2[k];

        t_5[k] = f_2 * sd_5[k]
                 + f_1 * pf_s_5[k]
                 + pb_x[k] * pd_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pb_y, pb_z, sf_6, sf_9, pf_s_6, pf_s_7, \
                         pf_s_8, pf_s_9, pd_3, pd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_x[k] * sf_6[k]
                 + f_1 * pf_s_6[k];

        t_7[k] = f_1 * pf_s_7[k]
                 + pb_z[k] * pd_3[k];

        t_8[k] = f_1 * pf_s_8[k]
                 + pb_y[k] * pd_5[k];

        t_9[k] = pa_x[k] * sf_9[k]
                 + f_1 * pf_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_x, pb_z, sd_0, sf_0, sf_1, pf_s_10, \
                         pf_s_11, pf_s_12, pf_s_13, pd_6, pd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * sf_0[k]
                  + f_1 * pf_s_10[k];

        t_11[k] = f_2 * sd_0[k]
                  + pa_y[k] * sf_1[k]
                  + f_1 * pf_s_11[k];

        t_12[k] = f_1 * pf_s_12[k]
                  + pb_z[k] * pd_6[k];

        t_13[k] = f_1 * pf_s_13[k]
                  + pb_x[k] * pd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pb_x, pb_z, sd_3, sf_6, pf_s_14, \
                         pf_s_15, pf_s_16, pf_s_17, pd_9, pd_10, \
                         pd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * pf_s_14[k]
                  + pb_x[k] * pd_10[k];

        t_15[k] = f_1 * pf_s_15[k]
                  + pb_x[k] * pd_11[k];

        t_16[k] = f_0 * sd_3[k]
                  + pa_y[k] * sf_6[k]
                  + f_1 * pf_s_16[k];

        t_17[k] = f_1 * pf_s_17[k]
                  + pb_z[k] * pd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_y, pa_z, pb_y, sd_5, sf_0, sf_9, pf_s_18, \
                         pf_s_19, pf_s_20, pf_s_21, pd_11, pd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * sd_5[k]
                  + f_1 * pf_s_18[k]
                  + pb_y[k] * pd_11[k];

        t_19[k] = pa_y[k] * sf_9[k]
                  + f_1 * pf_s_19[k];

        t_20[k] = pa_z[k] * sf_0[k]
                  + f_1 * pf_s_20[k];

        t_21[k] = f_1 * pf_s_21[k]
                  + pb_y[k] * pd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_z, pb_x, sd_0, sf_2, pf_s_22, pf_s_23, \
                         pf_s_24, pf_s_25, pd_15, pd_16, pd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * sd_0[k]
                  + pa_z[k] * sf_2[k]
                  + f_1 * pf_s_22[k];

        t_23[k] = f_1 * pf_s_23[k]
                  + pb_x[k] * pd_15[k];

        t_24[k] = f_1 * pf_s_24[k]
                  + pb_x[k] * pd_16[k];

        t_25[k] = f_1 * pf_s_25[k]
                  + pb_x[k] * pd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_y, sd_3, sd_5, sf_6, sf_7, sf_9, \
                         pf_s_26, pf_s_27, pf_s_28, pf_s_29, pd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * sf_6[k]
                  + f_1 * pf_s_26[k];

        t_27[k] = f_2 * sd_3[k]
                  + pa_z[k] * sf_7[k]
                  + f_1 * pf_s_27[k];

        t_28[k] = f_1 * pf_s_28[k]
                  + pb_y[k] * pd_17[k];

        t_29[k] = f_0 * sd_5[k]
                  + pa_z[k] * sf_9[k]
                  + f_1 * pf_s_29[k];
    }
}

}  // namespace simdkin
