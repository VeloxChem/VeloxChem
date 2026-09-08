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


#include "SimdKineticEnergyVrrRecPG.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_pg_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t sf, const size_t sg,
                                 const size_t pg_s, const size_t pf, const size_t ncols,
                                 const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

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
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *pg_s_0 = buffer.data(pg_s + 0);
    const auto *pg_s_1 = buffer.data(pg_s + 1);
    const auto *pg_s_2 = buffer.data(pg_s + 2);
    const auto *pg_s_3 = buffer.data(pg_s + 3);
    const auto *pg_s_4 = buffer.data(pg_s + 4);
    const auto *pg_s_5 = buffer.data(pg_s + 5);
    const auto *pg_s_6 = buffer.data(pg_s + 6);
    const auto *pg_s_7 = buffer.data(pg_s + 7);
    const auto *pg_s_8 = buffer.data(pg_s + 8);
    const auto *pg_s_9 = buffer.data(pg_s + 9);
    const auto *pg_s_10 = buffer.data(pg_s + 10);
    const auto *pg_s_11 = buffer.data(pg_s + 11);
    const auto *pg_s_12 = buffer.data(pg_s + 12);
    const auto *pg_s_13 = buffer.data(pg_s + 13);
    const auto *pg_s_14 = buffer.data(pg_s + 14);
    const auto *pg_s_15 = buffer.data(pg_s + 15);
    const auto *pg_s_16 = buffer.data(pg_s + 16);
    const auto *pg_s_17 = buffer.data(pg_s + 17);
    const auto *pg_s_18 = buffer.data(pg_s + 18);
    const auto *pg_s_19 = buffer.data(pg_s + 19);
    const auto *pg_s_20 = buffer.data(pg_s + 20);
    const auto *pg_s_21 = buffer.data(pg_s + 21);
    const auto *pg_s_22 = buffer.data(pg_s + 22);
    const auto *pg_s_23 = buffer.data(pg_s + 23);
    const auto *pg_s_24 = buffer.data(pg_s + 24);
    const auto *pg_s_25 = buffer.data(pg_s + 25);
    const auto *pg_s_26 = buffer.data(pg_s + 26);
    const auto *pg_s_27 = buffer.data(pg_s + 27);
    const auto *pg_s_28 = buffer.data(pg_s + 28);
    const auto *pg_s_29 = buffer.data(pg_s + 29);
    const auto *pg_s_30 = buffer.data(pg_s + 30);
    const auto *pg_s_31 = buffer.data(pg_s + 31);
    const auto *pg_s_32 = buffer.data(pg_s + 32);
    const auto *pg_s_33 = buffer.data(pg_s + 33);
    const auto *pg_s_34 = buffer.data(pg_s + 34);
    const auto *pg_s_35 = buffer.data(pg_s + 35);
    const auto *pg_s_36 = buffer.data(pg_s + 36);
    const auto *pg_s_37 = buffer.data(pg_s + 37);
    const auto *pg_s_38 = buffer.data(pg_s + 38);
    const auto *pg_s_39 = buffer.data(pg_s + 39);
    const auto *pg_s_40 = buffer.data(pg_s + 40);
    const auto *pg_s_41 = buffer.data(pg_s + 41);
    const auto *pg_s_42 = buffer.data(pg_s + 42);
    const auto *pg_s_43 = buffer.data(pg_s + 43);
    const auto *pg_s_44 = buffer.data(pg_s + 44);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_y, pb_z, sf_0, sf_3, sg_0, sg_3, pg_s_0, \
                         pg_s_1, pg_s_2, pg_s_3, pf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sf_0[k]
                 + pa_x[k] * sg_0[k]
                 + f_1 * pg_s_0[k];

        t_1[k] = f_1 * pg_s_1[k]
                 + pb_y[k] * pf_0[k];

        t_2[k] = f_1 * pg_s_2[k]
                 + pb_z[k] * pf_0[k];

        t_3[k] = f_2 * sf_3[k]
                 + pa_x[k] * sg_3[k]
                 + f_1 * pg_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_x, pb_y, sf_5, sf_6, sg_5, pg_s_4, pg_s_5, \
                         pg_s_6, pf_2, pf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * pg_s_4[k]
                 + pb_y[k] * pf_2[k];

        t_5[k] = f_2 * sf_5[k]
                 + pa_x[k] * sg_5[k]
                 + f_1 * pg_s_5[k];

        t_6[k] = f_3 * sf_6[k]
                 + f_1 * pg_s_6[k]
                 + pb_x[k] * pf_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, sf_9, pg_s_7, pg_s_8, pg_s_9, pf_3, \
                         pf_5, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pg_s_7[k]
                 + pb_z[k] * pf_3[k];

        t_8[k] = f_1 * pg_s_8[k]
                 + pb_y[k] * pf_5[k];

        t_9[k] = f_3 * sf_9[k]
                 + f_1 * pg_s_9[k]
                 + pb_x[k] * pf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_y, pb_z, sg_10, sg_12, pg_s_10, \
                         pg_s_11, pg_s_12, pg_s_13, pf_6, pf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * sg_10[k]
                  + f_1 * pg_s_10[k];

        t_11[k] = f_1 * pg_s_11[k]
                  + pb_z[k] * pf_6[k];

        t_12[k] = pa_x[k] * sg_12[k]
                  + f_1 * pg_s_12[k];

        t_13[k] = f_1 * pg_s_13[k]
                  + pb_y[k] * pf_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_y, pb_z, sf_0, sg_0, sg_1, sg_14, \
                         pg_s_14, pg_s_15, pg_s_16, pg_s_17, pf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_x[k] * sg_14[k]
                  + f_1 * pg_s_14[k];

        t_15[k] = pa_y[k] * sg_0[k]
                  + f_1 * pg_s_15[k];

        t_16[k] = f_3 * sf_0[k]
                  + pa_y[k] * sg_1[k]
                  + f_1 * pg_s_16[k];

        t_17[k] = f_1 * pg_s_17[k]
                  + pb_z[k] * pf_10[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_y, pb_x, pb_z, sf_1, sg_3, sg_5, pg_s_18, \
                         pg_s_19, pg_s_20, pg_s_21, pf_11, pf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * sf_1[k]
                  + pa_y[k] * sg_3[k]
                  + f_1 * pg_s_18[k];

        t_19[k] = f_1 * pg_s_19[k]
                  + pb_z[k] * pf_11[k];

        t_20[k] = pa_y[k] * sg_5[k]
                  + f_1 * pg_s_20[k];

        t_21[k] = f_1 * pg_s_21[k]
                  + pb_x[k] * pf_16[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_y, pb_x, sf_6, sg_10, pg_s_22, pg_s_23, \
                         pg_s_24, pg_s_25, pf_17, pf_18, pf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * pg_s_22[k]
                  + pb_x[k] * pf_17[k];

        t_23[k] = f_1 * pg_s_23[k]
                  + pb_x[k] * pf_18[k];

        t_24[k] = f_1 * pg_s_24[k]
                  + pb_x[k] * pf_19[k];

        t_25[k] = f_0 * sf_6[k]
                  + pa_y[k] * sg_10[k]
                  + f_1 * pg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_y, pb_z, sf_8, sf_9, sg_12, pg_s_26, \
                         pg_s_27, pg_s_28, pf_16, pf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * pg_s_26[k]
                  + pb_z[k] * pf_16[k];

        t_27[k] = f_2 * sf_8[k]
                  + pa_y[k] * sg_12[k]
                  + f_1 * pg_s_27[k];

        t_28[k] = f_3 * sf_9[k]
                  + f_1 * pg_s_28[k]
                  + pb_y[k] * pf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_y, sf_0, sg_0, sg_2, sg_14, \
                         pg_s_29, pg_s_30, pg_s_31, pg_s_32, pf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * sg_14[k]
                  + f_1 * pg_s_29[k];

        t_30[k] = pa_z[k] * sg_0[k]
                  + f_1 * pg_s_30[k];

        t_31[k] = f_1 * pg_s_31[k]
                  + pb_y[k] * pf_20[k];

        t_32[k] = f_3 * sf_0[k]
                  + pa_z[k] * sg_2[k]
                  + f_1 * pg_s_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_x, pb_y, sf_2, sg_3, sg_5, pg_s_33, \
                         pg_s_34, pg_s_35, pg_s_36, pf_22, pf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * sg_3[k]
                  + f_1 * pg_s_33[k];

        t_34[k] = f_1 * pg_s_34[k]
                  + pb_y[k] * pf_22[k];

        t_35[k] = f_2 * sf_2[k]
                  + pa_z[k] * sg_5[k]
                  + f_1 * pg_s_35[k];

        t_36[k] = f_1 * pg_s_36[k]
                  + pb_x[k] * pf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, sg_10, pg_s_37, pg_s_38, pg_s_39, \
                         pg_s_40, pf_27, pf_28, pf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_1 * pg_s_37[k]
                  + pb_x[k] * pf_27[k];

        t_38[k] = f_1 * pg_s_38[k]
                  + pb_x[k] * pf_28[k];

        t_39[k] = f_1 * pg_s_39[k]
                  + pb_x[k] * pf_29[k];

        t_40[k] = pa_z[k] * sg_10[k]
                  + f_1 * pg_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_z, pb_y, sf_6, sf_7, sg_11, sg_12, pg_s_41, \
                         pg_s_42, pg_s_43, pf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_3 * sf_6[k]
                  + pa_z[k] * sg_11[k]
                  + f_1 * pg_s_41[k];

        t_42[k] = f_2 * sf_7[k]
                  + pa_z[k] * sg_12[k]
                  + f_1 * pg_s_42[k];

        t_43[k] = f_1 * pg_s_43[k]
                  + pb_y[k] * pf_29[k];
    }

#pragma omp simd aligned(t_44, pa_z, sf_9, sg_14, pg_s_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * sf_9[k]
                  + pa_z[k] * sg_14[k]
                  + f_1 * pg_s_44[k];
    }
}

}  // namespace simdkin
