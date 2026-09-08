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


#include "SimdElectronRepulsionVrrRecPG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_pg_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t sf, const size_t sg,
                                     const size_t pf, const size_t ncols,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 1.0 / p;
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
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_2 = buffer.data(sf + 2);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_6 = buffer.data(sf + 6);
    const auto *sf_7 = buffer.data(sf + 7);
    const auto *sf_8 = buffer.data(sf + 8);
    const auto *sf_9 = buffer.data(sf + 9);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_y, pb_z, sf_0, sf_3, sg_0, sg_3, \
                         pf_0, pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * pf_0[k];

        t_2[k] = pb_z[k] * pf_0[k];

        t_3[k] = f_1 * sf_3[k]
                 + pa_x[k] * sg_3[k];

        t_4[k] = pb_y[k] * pf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_x, pb_x, pb_y, pb_z, sf_5, sf_6, sg_5, pf_3, \
                         pf_5, pf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sf_5[k]
                 + pa_x[k] * sg_5[k];

        t_6[k] = f_2 * sf_6[k]
                 + pb_x[k] * pf_6[k];

        t_7[k] = pb_z[k] * pf_3[k];

        t_8[k] = pb_y[k] * pf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, t_14, pa_x, pb_x, pb_y, pb_z, sf_9, \
                         sg_10, sg_12, sg_14, pf_6, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * sf_9[k]
                 + pb_x[k] * pf_9[k];

        t_10[k] = pa_x[k] * sg_10[k];

        t_11[k] = pb_z[k] * pf_6[k];

        t_12[k] = pa_x[k] * sg_12[k];

        t_13[k] = pb_y[k] * pf_9[k];

        t_14[k] = pa_x[k] * sg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, sf_0, sf_1, \
                         sg_0, sg_3, sg_5, pf_10, pf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * sg_0[k];

        t_16[k] = f_2 * sf_0[k]
                  + pb_y[k] * pf_10[k];

        t_17[k] = pb_z[k] * pf_10[k];

        t_18[k] = f_1 * sf_1[k]
                  + pa_y[k] * sg_3[k];

        t_19[k] = pb_z[k] * pf_11[k];

        t_20[k] = pa_y[k] * sg_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pa_y, pb_x, pb_z, sf_6, sg_10, \
                         pf_16, pf_17, pf_18, pf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_x[k] * pf_16[k];

        t_22[k] = pb_x[k] * pf_17[k];

        t_23[k] = pb_x[k] * pf_18[k];

        t_24[k] = pb_x[k] * pf_19[k];

        t_25[k] = f_0 * sf_6[k]
                  + pa_y[k] * sg_10[k];

        t_26[k] = pb_z[k] * pf_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, sf_8, sf_9, sg_0, \
                         sg_12, sg_14, pf_19, pf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * sf_8[k]
                  + pa_y[k] * sg_12[k];

        t_28[k] = f_2 * sf_9[k]
                  + pb_y[k] * pf_19[k];

        t_29[k] = pa_y[k] * sg_14[k];

        t_30[k] = pa_z[k] * sg_0[k];

        t_31[k] = pb_y[k] * pf_20[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_z, pb_x, pb_y, pb_z, sf_0, sf_2, \
                         sg_3, sg_5, pf_20, pf_22, pf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * sf_0[k]
                  + pb_z[k] * pf_20[k];

        t_33[k] = pa_z[k] * sg_3[k];

        t_34[k] = pb_y[k] * pf_22[k];

        t_35[k] = f_1 * sf_2[k]
                  + pa_z[k] * sg_5[k];

        t_36[k] = pb_x[k] * pf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_z, pb_x, pb_z, sf_6, sg_10, pf_26, \
                         pf_27, pf_28, pf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_x[k] * pf_27[k];

        t_38[k] = pb_x[k] * pf_28[k];

        t_39[k] = pb_x[k] * pf_29[k];

        t_40[k] = pa_z[k] * sg_10[k];

        t_41[k] = f_2 * sf_6[k]
                  + pb_z[k] * pf_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_z, pb_y, sf_7, sf_9, sg_12, sg_14, \
                         pf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * sf_7[k]
                  + pa_z[k] * sg_12[k];

        t_43[k] = pb_y[k] * pf_29[k];

        t_44[k] = f_0 * sf_9[k]
                  + pa_z[k] * sg_14[k];
    }
}

}  // namespace simdt2ceri
