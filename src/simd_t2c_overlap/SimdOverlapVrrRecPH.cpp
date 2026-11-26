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


#include "SimdOverlapVrrRecPH.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_ph_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sg, const size_t sh, const size_t pf,
                          const size_t pg, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;

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
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_20 = buffer.data(sh + 20);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_21 = buffer.data(pg + 21);
    const auto *pg_23 = buffer.data(pg + 23);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_26 = buffer.data(pg + 26);
    const auto *pg_27 = buffer.data(pg + 27);
    const auto *pg_28 = buffer.data(pg + 28);
    const auto *pg_29 = buffer.data(pg + 29);
    const auto *pg_30 = buffer.data(pg + 30);
    const auto *pg_32 = buffer.data(pg + 32);
    const auto *pg_35 = buffer.data(pg + 35);
    const auto *pg_37 = buffer.data(pg + 37);
    const auto *pg_39 = buffer.data(pg + 39);
    const auto *pg_40 = buffer.data(pg + 40);
    const auto *pg_41 = buffer.data(pg + 41);
    const auto *pg_42 = buffer.data(pg + 42);
    const auto *pg_43 = buffer.data(pg + 43);
    const auto *pg_44 = buffer.data(pg + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pb_y, pb_z, sg_0, sh_0, pf_0, \
                         pg_0, pg_1, pg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sg_0[k]
                 + pa_x[k] * sh_0[k];

        t_1[k] = pb_y[k] * pg_0[k];

        t_2[k] = pb_z[k] * pg_0[k];

        t_3[k] = f_1 * pf_0[k]
                 + pb_y[k] * pg_1[k];

        t_4[k] = pb_y[k] * pg_2[k];

        t_5[k] = f_1 * pf_0[k]
                 + pb_z[k] * pg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, sg_10, pf_1, pf_2, \
                         pg_3, pg_5, pg_6, pg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * pf_1[k]
                 + pb_y[k] * pg_3[k];

        t_7[k] = pb_z[k] * pg_3[k];

        t_8[k] = pb_y[k] * pg_5[k];

        t_9[k] = f_2 * pf_2[k]
                 + pb_z[k] * pg_5[k];

        t_10[k] = f_1 * sg_10[k]
                  + pb_x[k] * pg_10[k];

        t_11[k] = pb_z[k] * pg_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_x, pb_x, pb_y, pb_z, sg_12, sg_14, \
                         sh_15, pg_9, pg_10, pg_12, pg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * sg_12[k]
                  + pb_x[k] * pg_12[k];

        t_13[k] = pb_y[k] * pg_9[k];

        t_14[k] = f_1 * sg_14[k]
                  + pb_x[k] * pg_14[k];

        t_15[k] = pa_x[k] * sh_15[k];

        t_16[k] = pb_z[k] * pg_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, sh_0, sh_17, sh_18, \
                         sh_20, pg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_x[k] * sh_17[k];

        t_18[k] = pa_x[k] * sh_18[k];

        t_19[k] = pb_y[k] * pg_14[k];

        t_20[k] = pa_x[k] * sh_20[k];

        t_21[k] = pa_y[k] * sh_0[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_y, pb_x, pb_z, sh_5, pf_11, pf_13, \
                         pg_15, pg_16, pg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pf_11[k]
                  + pb_x[k] * pg_16[k];

        t_23[k] = pb_z[k] * pg_15[k];

        t_24[k] = f_2 * pf_13[k]
                  + pb_x[k] * pg_18[k];

        t_25[k] = pb_z[k] * pg_16[k];

        t_26[k] = pa_y[k] * sh_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pb_x, pb_z, sh_9, pf_16, pf_18, \
                         pg_18, pg_21, pg_23, pg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * pf_16[k]
                  + pb_x[k] * pg_21[k];

        t_28[k] = pb_z[k] * pg_18[k];

        t_29[k] = f_1 * pf_18[k]
                  + pb_x[k] * pg_23[k];

        t_30[k] = pa_y[k] * sh_9[k];

        t_31[k] = pb_x[k] * pg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_z, sg_10, sh_15, \
                         pg_25, pg_26, pg_27, pg_28, pg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = pb_x[k] * pg_26[k];

        t_33[k] = pb_x[k] * pg_27[k];

        t_34[k] = pb_x[k] * pg_28[k];

        t_35[k] = pb_x[k] * pg_29[k];

        t_36[k] = f_0 * sg_10[k]
                  + pa_y[k] * sh_15[k];

        t_37[k] = pb_z[k] * pg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pb_y, pb_z, sg_14, sh_20, pf_16, pf_17, \
                         pg_26, pg_27, pg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_1 * pf_16[k]
                  + pb_z[k] * pg_26[k];

        t_39[k] = f_2 * pf_17[k]
                  + pb_z[k] * pg_27[k];

        t_40[k] = f_1 * sg_14[k]
                  + pb_y[k] * pg_29[k];

        t_41[k] = pa_y[k] * sh_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_z, pb_x, pb_y, sh_0, sh_3, \
                         pf_22, pf_25, pg_30, pg_32, pg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * sh_0[k];

        t_43[k] = pb_y[k] * pg_30[k];

        t_44[k] = f_3 * pf_22[k]
                  + pb_x[k] * pg_32[k];

        t_45[k] = pa_z[k] * sh_3[k];

        t_46[k] = pb_y[k] * pg_32[k];

        t_47[k] = f_2 * pf_25[k]
                  + pb_x[k] * pg_35[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_z, pb_x, pb_y, sh_6, pf_27, pf_29, \
                         pg_35, pg_37, pg_39, pg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * sh_6[k];

        t_49[k] = f_1 * pf_27[k]
                  + pb_x[k] * pg_37[k];

        t_50[k] = pb_y[k] * pg_35[k];

        t_51[k] = f_1 * pf_29[k]
                  + pb_x[k] * pg_39[k];

        t_52[k] = pb_x[k] * pg_40[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, pa_z, pb_x, pb_y, sh_15, pf_27, \
                         pg_41, pg_42, pg_43, pg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_x[k] * pg_41[k];

        t_54[k] = pb_x[k] * pg_42[k];

        t_55[k] = pb_x[k] * pg_43[k];

        t_56[k] = pb_x[k] * pg_44[k];

        t_57[k] = pa_z[k] * sh_15[k];

        t_58[k] = f_3 * pf_27[k]
                  + pb_y[k] * pg_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_y, sg_14, sh_20, pf_28, pf_29, \
                         pg_42, pg_43, pg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_2 * pf_28[k]
                  + pb_y[k] * pg_42[k];

        t_60[k] = f_1 * pf_29[k]
                  + pb_y[k] * pg_43[k];

        t_61[k] = pb_y[k] * pg_44[k];

        t_62[k] = f_0 * sg_14[k]
                  + pa_z[k] * sh_20[k];
    }
}

}  // namespace simdovl
