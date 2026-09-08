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


#include "SimdOverlapVrrRecPG.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_pg_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_19 = buffer.data(pf + 19);
    const auto *pf_20 = buffer.data(pf + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pb_y, pb_z, sf_0, sg_0, pd_0, \
                         pf_0, pf_1, pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * pf_0[k];

        t_2[k] = pb_z[k] * pf_0[k];

        t_3[k] = f_1 * pd_0[k]
                 + pb_y[k] * pf_1[k];

        t_4[k] = pb_y[k] * pf_2[k];

        t_5[k] = f_1 * pd_0[k]
                 + pb_z[k] * pf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pb_x, pb_y, pb_z, sf_5, sf_8, sg_5, \
                         pf_3, pf_4, pf_5, pf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_1 * sf_5[k]
                 + pb_x[k] * pf_5[k];

        t_7[k] = pb_z[k] * pf_3[k];

        t_8[k] = pb_y[k] * pf_4[k];

        t_9[k] = f_1 * sf_8[k]
                 + pb_x[k] * pf_6[k];

        t_10[k] = pa_x[k] * sg_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, pb_y, pb_z, sg_0, sg_7, \
                         sg_8, pf_5, pf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * pf_5[k];

        t_12[k] = pa_x[k] * sg_7[k];

        t_13[k] = pb_y[k] * pf_6[k];

        t_14[k] = pa_x[k] * sg_8[k];

        t_15[k] = pa_y[k] * sg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, pa_y, pb_x, pb_z, sg_4, pd_3, \
                         pd_4, pf_7, pf_8, pf_9, pf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * pd_3[k]
                  + pb_x[k] * pf_8[k];

        t_17[k] = pb_z[k] * pf_7[k];

        t_18[k] = f_1 * pd_4[k]
                  + pb_x[k] * pf_9[k];

        t_19[k] = pb_z[k] * pf_8[k];

        t_20[k] = pa_y[k] * sg_4[k];

        t_21[k] = pb_x[k] * pf_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pa_y, pb_x, pb_z, sf_5, sg_5, \
                         pd_4, pf_10, pf_11, pf_12, pf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_x[k] * pf_11[k];

        t_23[k] = pb_x[k] * pf_12[k];

        t_24[k] = pb_x[k] * pf_13[k];

        t_25[k] = f_0 * sf_5[k]
                  + pa_y[k] * sg_5[k];

        t_26[k] = pb_z[k] * pf_10[k];

        t_27[k] = f_1 * pd_4[k]
                  + pb_z[k] * pf_11[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_y, pa_z, pb_x, pb_y, sf_8, sg_0, \
                         sg_8, pd_6, pf_13, pf_14, pf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * sf_8[k]
                  + pb_y[k] * pf_13[k];

        t_29[k] = pa_y[k] * sg_8[k];

        t_30[k] = pa_z[k] * sg_0[k];

        t_31[k] = pb_y[k] * pf_14[k];

        t_32[k] = f_2 * pd_6[k]
                  + pb_x[k] * pf_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_z, pb_x, pb_y, sg_3, pd_8, \
                         pf_15, pf_16, pf_17, pf_18, pf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * sg_3[k];

        t_34[k] = pb_y[k] * pf_15[k];

        t_35[k] = f_1 * pd_8[k]
                  + pb_x[k] * pf_16[k];

        t_36[k] = pb_x[k] * pf_17[k];

        t_37[k] = pb_x[k] * pf_18[k];

        t_38[k] = pb_x[k] * pf_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_z, pb_x, pb_y, sg_5, pd_7, pd_8, \
                         pf_18, pf_19, pf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_x[k] * pf_20[k];

        t_40[k] = pa_z[k] * sg_5[k];

        t_41[k] = f_2 * pd_7[k]
                  + pb_y[k] * pf_18[k];

        t_42[k] = f_1 * pd_8[k]
                  + pb_y[k] * pf_19[k];

        t_43[k] = pb_y[k] * pf_20[k];
    }

#pragma omp simd aligned(t_44, pa_z, sf_8, sg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * sf_8[k]
                  + pa_z[k] * sg_8[k];
    }
}

auto
compute_prim_pg_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_y, pb_z, sf_0, sg_0, pd_0, pf_0, \
                         pf_1, pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * pf_0[k];

        t_2[k] = pb_z[k] * pf_0[k];

        t_3[k] = f_1 * pd_0[k]
                 + pb_y[k] * pf_1[k];

        t_4[k] = f_1 * pd_0[k]
                 + pb_z[k] * pf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_x, pb_x, sf_5, sf_8, sg_7, sg_9, sg_11, \
                         pf_3, pf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sf_5[k]
                 + pb_x[k] * pf_3[k];

        t_6[k] = f_1 * sf_8[k]
                 + pb_x[k] * pf_4[k];

        t_7[k] = pa_x[k] * sg_7[k];

        t_8[k] = pa_x[k] * sg_9[k];

        t_9[k] = pa_x[k] * sg_11[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, sg_0, pd_1, pd_2, pf_5, \
                         pf_6, pf_7, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * sg_0[k];

        t_11[k] = f_2 * pd_1[k]
                  + pb_x[k] * pf_5[k];

        t_12[k] = f_1 * pd_2[k]
                  + pb_x[k] * pf_6[k];

        t_13[k] = pb_x[k] * pf_7[k];

        t_14[k] = pb_x[k] * pf_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, sf_5, sf_8, sg_7, \
                         sg_11, pd_2, pf_7, pf_8, pf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * sf_5[k]
                  + pa_y[k] * sg_7[k];

        t_16[k] = pb_z[k] * pf_7[k];

        t_17[k] = f_1 * pd_2[k]
                  + pb_z[k] * pf_8[k];

        t_18[k] = f_1 * sf_8[k]
                  + pb_y[k] * pf_10[k];

        t_19[k] = pa_y[k] * sg_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pa_z, pb_x, sg_0, sg_7, pd_4, \
                         pd_6, pf_11, pf_12, pf_13, pf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * sg_0[k];

        t_21[k] = f_2 * pd_4[k]
                  + pb_x[k] * pf_11[k];

        t_22[k] = f_1 * pd_6[k]
                  + pb_x[k] * pf_12[k];

        t_23[k] = pb_x[k] * pf_13[k];

        t_24[k] = pb_x[k] * pf_15[k];

        t_25[k] = pa_z[k] * sg_7[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_y, sf_8, sg_11, pd_5, pd_6, pf_13, \
                         pf_14, pf_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * pd_5[k]
                  + pb_y[k] * pf_13[k];

        t_27[k] = f_1 * pd_6[k]
                  + pb_y[k] * pf_14[k];

        t_28[k] = pb_y[k] * pf_15[k];

        t_29[k] = f_0 * sf_8[k]
                  + pa_z[k] * sg_11[k];
    }
}

auto
compute_prim_pg_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_3 = buffer.data(sf + 3);
    const auto *sf_5 = buffer.data(sf + 5);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pb_y, pb_z, sf_0, sg_0, sg_7, \
                         pd_0, pf_0, pf_1, pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * pf_0[k];

        t_2[k] = pb_z[k] * pf_0[k];

        t_3[k] = f_1 * pd_0[k]
                 + pb_y[k] * pf_1[k];

        t_4[k] = f_1 * pd_0[k]
                 + pb_z[k] * pf_2[k];

        t_5[k] = pa_x[k] * sg_7[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pb_x, sg_11, pd_1, pd_2, pf_3, pf_4, \
                         pf_5, pf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_x[k] * sg_11[k];

        t_7[k] = f_2 * pd_1[k]
                 + pb_x[k] * pf_3[k];

        t_8[k] = f_1 * pd_2[k]
                 + pb_x[k] * pf_4[k];

        t_9[k] = pb_x[k] * pf_5[k];

        t_10[k] = pb_x[k] * pf_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, sf_3, sf_5, sg_7, \
                         sg_11, pd_2, pf_5, pf_6, pf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * sf_3[k]
                  + pa_y[k] * sg_7[k];

        t_12[k] = pb_z[k] * pf_5[k];

        t_13[k] = f_1 * pd_2[k]
                  + pb_z[k] * pf_6[k];

        t_14[k] = f_1 * sf_5[k]
                  + pb_y[k] * pf_8[k];

        t_15[k] = pa_y[k] * sg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, pd_4, pd_5, pd_6, \
                         pf_9, pf_10, pf_11, pf_12, pf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * pd_4[k]
                  + pb_x[k] * pf_9[k];

        t_17[k] = f_1 * pd_6[k]
                  + pb_x[k] * pf_10[k];

        t_18[k] = pb_x[k] * pf_11[k];

        t_19[k] = pb_x[k] * pf_13[k];

        t_20[k] = f_2 * pd_5[k]
                  + pb_y[k] * pf_11[k];

        t_21[k] = f_1 * pd_6[k]
                  + pb_y[k] * pf_12[k];
    }

#pragma omp simd aligned(t_22, t_23, pa_z, pb_y, sf_5, sg_11, pf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_y[k] * pf_13[k];

        t_23[k] = f_0 * sf_5[k]
                  + pa_z[k] * sg_11[k];
    }
}

auto
compute_prim_pg_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_8 = buffer.data(sg + 8);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_x, pb_y, pb_z, sf_0, sf_5, sg_0, pd_0, \
                         pf_1, pf_2, pf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = f_1 * pd_0[k]
                 + pb_y[k] * pf_1[k];

        t_2[k] = f_1 * pd_0[k]
                 + pb_z[k] * pf_2[k];

        t_3[k] = f_1 * sf_5[k]
                 + pb_x[k] * pf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, sf_5, sf_8, sg_5, pd_3, pd_4, pf_4, \
                         pf_5, pf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * sf_8[k]
                 + pb_x[k] * pf_4[k];

        t_5[k] = f_2 * pd_3[k]
                 + pb_x[k] * pf_5[k];

        t_6[k] = f_1 * pd_4[k]
                 + pb_x[k] * pf_6[k];

        t_7[k] = f_0 * sf_5[k]
                 + pa_y[k] * sg_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pa_z, pb_y, pb_z, sf_8, sg_0, sg_8, pd_4, \
                         pf_8, pf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * pd_4[k]
                 + pb_z[k] * pf_8[k];

        t_9[k] = f_1 * sf_8[k]
                 + pb_y[k] * pf_10[k];

        t_10[k] = pa_y[k] * sg_8[k];

        t_11[k] = pa_z[k] * sg_0[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_z, pb_x, pb_y, sg_5, pd_6, pd_7, \
                         pd_8, pf_11, pf_12, pf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * pd_6[k]
                  + pb_x[k] * pf_11[k];

        t_13[k] = pb_y[k] * pf_11[k];

        t_14[k] = f_1 * pd_8[k]
                  + pb_x[k] * pf_12[k];

        t_15[k] = pa_z[k] * sg_5[k];

        t_16[k] = f_2 * pd_7[k]
                  + pb_y[k] * pf_13[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_z, pb_y, sf_8, sg_8, pd_8, \
                         pf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * pd_8[k]
                  + pb_y[k] * pf_14[k];

        t_18[k] = f_0 * sf_8[k]
                  + pa_z[k] * sg_8[k];
    }
}

auto
compute_prim_pg_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pb_y, pb_z, sf_0, sg_0, sg_7, \
                         pd_0, pf_0, pf_1, pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * pf_0[k];

        t_2[k] = pb_z[k] * pf_0[k];

        t_3[k] = f_1 * pd_0[k]
                 + pb_y[k] * pf_1[k];

        t_4[k] = f_1 * pd_0[k]
                 + pb_z[k] * pf_2[k];

        t_5[k] = pa_x[k] * sg_7[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pb_x, sg_9, sg_11, pd_1, pd_2, pf_5, \
                         pf_6, pf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_x[k] * sg_9[k];

        t_7[k] = pa_x[k] * sg_11[k];

        t_8[k] = f_2 * pd_1[k]
                 + pb_x[k] * pf_5[k];

        t_9[k] = f_1 * pd_2[k]
                 + pb_x[k] * pf_6[k];

        t_10[k] = pb_x[k] * pf_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, sf_5, sf_8, sg_7, \
                         sg_11, pd_2, pf_7, pf_8, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * sf_5[k]
                  + pa_y[k] * sg_7[k];

        t_12[k] = pb_z[k] * pf_7[k];

        t_13[k] = f_1 * pd_2[k]
                  + pb_z[k] * pf_8[k];

        t_14[k] = f_1 * sf_8[k]
                  + pb_y[k] * pf_9[k];

        t_15[k] = pa_y[k] * sg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, pd_4, pd_5, pd_6, \
                         pf_10, pf_11, pf_12, pf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * pd_4[k]
                  + pb_x[k] * pf_10[k];

        t_17[k] = pb_y[k] * pf_10[k];

        t_18[k] = f_1 * pd_6[k]
                  + pb_x[k] * pf_11[k];

        t_19[k] = pb_x[k] * pf_12[k];

        t_20[k] = pb_x[k] * pf_14[k];

        t_21[k] = f_2 * pd_5[k]
                  + pb_y[k] * pf_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_z, pb_y, sf_8, sg_11, pd_6, pf_13, \
                         pf_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * pd_6[k]
                  + pb_y[k] * pf_13[k];

        t_23[k] = pb_y[k] * pf_14[k];

        t_24[k] = f_0 * sf_8[k]
                  + pa_z[k] * sg_11[k];
    }
}

auto
compute_prim_pg_overlap_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pb_y, pb_z, sf_0, sg_0, sg_7, \
                         pd_0, pf_0, pf_1, pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = pb_y[k] * pf_0[k];

        t_2[k] = pb_z[k] * pf_0[k];

        t_3[k] = f_1 * pd_0[k]
                 + pb_y[k] * pf_1[k];

        t_4[k] = f_1 * pd_0[k]
                 + pb_z[k] * pf_2[k];

        t_5[k] = pa_x[k] * sg_7[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pb_x, sf_5, sg_7, sg_11, pd_1, \
                         pd_2, pf_3, pf_4, pf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_x[k] * sg_11[k];

        t_7[k] = f_2 * pd_1[k]
                 + pb_x[k] * pf_3[k];

        t_8[k] = f_1 * pd_2[k]
                 + pb_x[k] * pf_4[k];

        t_9[k] = pb_x[k] * pf_5[k];

        t_10[k] = f_0 * sf_5[k]
                  + pa_y[k] * sg_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_x, pb_z, sg_11, pd_2, pd_4, \
                         pd_6, pf_5, pf_6, pf_8, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * pf_5[k];

        t_12[k] = f_1 * pd_2[k]
                  + pb_z[k] * pf_6[k];

        t_13[k] = pa_y[k] * sg_11[k];

        t_14[k] = f_2 * pd_4[k]
                  + pb_x[k] * pf_8[k];

        t_15[k] = f_1 * pd_6[k]
                  + pb_x[k] * pf_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, pa_z, pb_x, pb_y, sf_8, sg_11, \
                         pd_5, pd_6, pf_10, pf_11, pf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_x[k] * pf_10[k];

        t_17[k] = pb_x[k] * pf_12[k];

        t_18[k] = f_2 * pd_5[k]
                  + pb_y[k] * pf_10[k];

        t_19[k] = f_1 * pd_6[k]
                  + pb_y[k] * pf_11[k];

        t_20[k] = pb_y[k] * pf_12[k];

        t_21[k] = f_0 * sf_8[k]
                  + pa_z[k] * sg_11[k];
    }
}

auto
compute_prim_pg_overlap_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_5 = buffer.data(sg + 5);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_4 = buffer.data(pd + 4);

    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_9 = buffer.data(pf + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pb_x, sf_0, sf_5, sg_0, sg_3, pd_1, \
                         pd_4, pf_3, pf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = f_1 * pd_1[k]
                 + pb_x[k] * pf_3[k];

        t_2[k] = f_0 * sf_5[k]
                 + pa_y[k] * sg_3[k];

        t_3[k] = f_1 * pd_4[k]
                 + pb_x[k] * pf_7[k];
    }

#pragma omp simd aligned(t_4, t_5, pa_z, pb_y, sf_8, sg_5, pd_4, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * pd_4[k]
                 + pb_y[k] * pf_9[k];

        t_5[k] = f_0 * sf_8[k]
                 + pa_z[k] * sg_5[k];
    }
}

auto
compute_prim_pg_overlap_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_x, pb_z, sf_0, sg_0, sg_7, sg_11, \
                         pd_0, pd_1, pf_1, pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = f_1 * pd_0[k]
                 + pb_z[k] * pf_1[k];

        t_2[k] = pa_x[k] * sg_7[k];

        t_3[k] = pa_x[k] * sg_11[k];

        t_4[k] = f_1 * pd_1[k]
                 + pb_x[k] * pf_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_y, pb_z, sf_5, sf_8, sg_7, pd_1, \
                         pf_3, pf_4, pf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_x[k] * pf_3[k];

        t_6[k] = f_0 * sf_5[k]
                 + pa_y[k] * sg_7[k];

        t_7[k] = f_1 * pd_1[k]
                 + pb_z[k] * pf_4[k];

        t_8[k] = f_1 * sf_8[k]
                 + pb_y[k] * pf_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, sg_11, pd_3, \
                         pd_4, pf_6, pf_7, pf_8, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = pa_y[k] * sg_11[k];

        t_10[k] = f_1 * pd_4[k]
                  + pb_x[k] * pf_6[k];

        t_11[k] = pb_x[k] * pf_9[k];

        t_12[k] = f_2 * pd_3[k]
                  + pb_y[k] * pf_7[k];

        t_13[k] = f_1 * pd_4[k]
                  + pb_y[k] * pf_8[k];

        t_14[k] = pb_y[k] * pf_9[k];
    }

#pragma omp simd aligned(t_15, pa_z, sf_8, sg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * sf_8[k]
                  + pa_z[k] * sg_11[k];
    }
}

auto
compute_prim_pg_overlap_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_6 = buffer.data(pf + 6);
    const auto *pf_7 = buffer.data(pf + 7);
    const auto *pf_8 = buffer.data(pf + 8);
    const auto *pf_9 = buffer.data(pf + 9);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pb_x, pb_z, sf_0, sg_0, pd_0, pd_1, \
                         pf_0, pf_1, pf_2, pf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = pb_z[k] * pf_0[k];

        t_2[k] = f_1 * pd_0[k]
                 + pb_z[k] * pf_1[k];

        t_3[k] = f_1 * pd_1[k]
                 + pb_x[k] * pf_2[k];

        t_4[k] = pb_x[k] * pf_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_x, pb_z, sf_5, sg_7, sg_11, pd_1, \
                         pd_4, pf_3, pf_4, pf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * sf_5[k]
                 + pa_y[k] * sg_7[k];

        t_6[k] = pb_z[k] * pf_3[k];

        t_7[k] = f_1 * pd_1[k]
                 + pb_z[k] * pf_4[k];

        t_8[k] = pa_y[k] * sg_11[k];

        t_9[k] = f_1 * pd_4[k]
                 + pb_x[k] * pf_6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_z, pb_x, pb_y, sf_8, sg_11, pd_3, \
                         pd_4, pf_7, pf_8, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_x[k] * pf_9[k];

        t_11[k] = f_2 * pd_3[k]
                  + pb_y[k] * pf_7[k];

        t_12[k] = f_1 * pd_4[k]
                  + pb_y[k] * pf_8[k];

        t_13[k] = pb_y[k] * pf_9[k];

        t_14[k] = f_0 * sf_8[k]
                  + pa_z[k] * sg_11[k];
    }
}

auto
compute_prim_pg_overlap_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                          const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_4 = buffer.data(sf + 4);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_2 = buffer.data(sg + 2);

    const auto *pd_3 = buffer.data(pd + 3);

    const auto *pf_2 = buffer.data(pf + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pa_y, pa_z, pb_y, sf_1, sf_4, sg_0, sg_2, pd_3, \
                         pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_1[k]
                 + pa_y[k] * sg_0[k];

        t_1[k] = f_1 * pd_3[k]
                 + pb_y[k] * pf_2[k];

        t_2[k] = f_0 * sf_4[k]
                 + pa_z[k] * sg_2[k];
    }
}

auto
compute_prim_pg_overlap_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                           const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_7 = buffer.data(sg + 7);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_6 = buffer.data(pf + 6);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pb_x, sf_0, sf_5, sg_0, sg_3, pd_1, \
                         pd_2, pf_1, pf_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = f_1 * pd_1[k]
                 + pb_x[k] * pf_1[k];

        t_2[k] = f_0 * sf_5[k]
                 + pa_y[k] * sg_3[k];

        t_3[k] = f_1 * pd_2[k]
                 + pb_x[k] * pf_4[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_z, pb_y, sf_8, sg_7, pd_2, pf_5, \
                         pf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * pd_2[k]
                 + pb_y[k] * pf_5[k];

        t_5[k] = pb_y[k] * pf_6[k];

        t_6[k] = f_0 * sf_8[k]
                 + pa_z[k] * sg_7[k];
    }
}

auto
compute_prim_pg_overlap_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                           const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_5 = buffer.data(sf + 5);
    const auto *sf_8 = buffer.data(sf + 8);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_11 = buffer.data(sg + 11);

    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pa_x, pa_y, pb_x, sf_0, sf_5, sg_0, sg_7, \
                         sg_11, pd_1, pf_1, pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = f_1 * pd_1[k]
                 + pb_x[k] * pf_1[k];

        t_2[k] = pb_x[k] * pf_2[k];

        t_3[k] = f_0 * sf_5[k]
                 + pa_y[k] * sg_7[k];

        t_4[k] = pa_y[k] * sg_11[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_z, pb_x, pb_y, sf_8, sg_11, pd_2, pf_3, \
                         pf_4, pf_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * pd_2[k]
                 + pb_x[k] * pf_3[k];

        t_6[k] = pb_x[k] * pf_5[k];

        t_7[k] = f_1 * pd_2[k]
                 + pb_y[k] * pf_4[k];

        t_8[k] = pb_y[k] * pf_5[k];

        t_9[k] = f_0 * sf_8[k]
                 + pa_z[k] * sg_11[k];
    }
}

auto
compute_prim_pg_overlap_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                           const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_4 = buffer.data(sf + 4);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_2 = buffer.data(sg + 2);

    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_2 = buffer.data(pf + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pa_y, pa_z, pb_y, sf_1, sf_4, sg_0, sg_2, pd_2, \
                         pf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_1[k]
                 + pa_y[k] * sg_0[k];

        t_1[k] = f_1 * pd_2[k]
                 + pb_y[k] * pf_2[k];

        t_2[k] = f_0 * sf_4[k]
                 + pa_z[k] * sg_2[k];
    }
}

auto
compute_prim_pg_overlap_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                           const size_t pb, const size_t sf, const size_t sg, const size_t pd,
                           const size_t pf, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *sf_0 = buffer.data(sf + 0);
    const auto *sf_1 = buffer.data(sf + 1);
    const auto *sf_4 = buffer.data(sf + 4);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_5 = buffer.data(sg + 5);

    const auto *pd_2 = buffer.data(pd + 2);

    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pb_y, sf_0, sf_1, sg_0, sg_1, pd_2, \
                         pf_2, pf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k];

        t_1[k] = f_0 * sf_1[k]
                 + pa_y[k] * sg_1[k];

        t_2[k] = f_1 * pd_2[k]
                 + pb_y[k] * pf_2[k];

        t_3[k] = pb_y[k] * pf_3[k];
    }

#pragma omp simd aligned(t_4, pa_z, sf_4, sg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * sf_4[k]
                 + pa_z[k] * sg_5[k];
    }
}

}  // namespace simdovl
