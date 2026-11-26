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


#include "SimdOverlapVrrRecSN.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_sn_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                          const size_t sl, const size_t sm, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 3.5 / p;
    const auto f_2 = 3.0 / p;
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 0.5 / p;

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
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sl_0 = buffer.data(sl + 0);
    const auto *sl_3 = buffer.data(sl + 3);
    const auto *sl_5 = buffer.data(sl + 5);
    const auto *sl_6 = buffer.data(sl + 6);
    const auto *sl_9 = buffer.data(sl + 9);
    const auto *sl_10 = buffer.data(sl + 10);
    const auto *sl_12 = buffer.data(sl + 12);
    const auto *sl_14 = buffer.data(sl + 14);
    const auto *sl_15 = buffer.data(sl + 15);
    const auto *sl_17 = buffer.data(sl + 17);
    const auto *sl_18 = buffer.data(sl + 18);
    const auto *sl_20 = buffer.data(sl + 20);
    const auto *sl_21 = buffer.data(sl + 21);
    const auto *sl_23 = buffer.data(sl + 23);
    const auto *sl_24 = buffer.data(sl + 24);
    const auto *sl_25 = buffer.data(sl + 25);
    const auto *sl_27 = buffer.data(sl + 27);
    const auto *sl_28 = buffer.data(sl + 28);
    const auto *sl_30 = buffer.data(sl + 30);
    const auto *sl_31 = buffer.data(sl + 31);
    const auto *sl_32 = buffer.data(sl + 32);
    const auto *sl_33 = buffer.data(sl + 33);
    const auto *sl_35 = buffer.data(sl + 35);
    const auto *sl_36 = buffer.data(sl + 36);
    const auto *sl_38 = buffer.data(sl + 38);
    const auto *sl_39 = buffer.data(sl + 39);
    const auto *sl_40 = buffer.data(sl + 40);
    const auto *sl_41 = buffer.data(sl + 41);
    const auto *sl_42 = buffer.data(sl + 42);
    const auto *sl_43 = buffer.data(sl + 43);
    const auto *sl_44 = buffer.data(sl + 44);

    const auto *sm_0 = buffer.data(sm + 0);
    const auto *sm_2 = buffer.data(sm + 2);
    const auto *sm_3 = buffer.data(sm + 3);
    const auto *sm_5 = buffer.data(sm + 5);
    const auto *sm_6 = buffer.data(sm + 6);
    const auto *sm_9 = buffer.data(sm + 9);
    const auto *sm_10 = buffer.data(sm + 10);
    const auto *sm_12 = buffer.data(sm + 12);
    const auto *sm_14 = buffer.data(sm + 14);
    const auto *sm_15 = buffer.data(sm + 15);
    const auto *sm_17 = buffer.data(sm + 17);
    const auto *sm_18 = buffer.data(sm + 18);
    const auto *sm_20 = buffer.data(sm + 20);
    const auto *sm_21 = buffer.data(sm + 21);
    const auto *sm_23 = buffer.data(sm + 23);
    const auto *sm_24 = buffer.data(sm + 24);
    const auto *sm_25 = buffer.data(sm + 25);
    const auto *sm_27 = buffer.data(sm + 27);
    const auto *sm_28 = buffer.data(sm + 28);
    const auto *sm_30 = buffer.data(sm + 30);
    const auto *sm_31 = buffer.data(sm + 31);
    const auto *sm_32 = buffer.data(sm + 32);
    const auto *sm_33 = buffer.data(sm + 33);
    const auto *sm_35 = buffer.data(sm + 35);
    const auto *sm_36 = buffer.data(sm + 36);
    const auto *sm_38 = buffer.data(sm + 38);
    const auto *sm_39 = buffer.data(sm + 39);
    const auto *sm_40 = buffer.data(sm + 40);
    const auto *sm_41 = buffer.data(sm + 41);
    const auto *sm_42 = buffer.data(sm + 42);
    const auto *sm_44 = buffer.data(sm + 44);
    const auto *sm_45 = buffer.data(sm + 45);
    const auto *sm_46 = buffer.data(sm + 46);
    const auto *sm_47 = buffer.data(sm + 47);
    const auto *sm_48 = buffer.data(sm + 48);
    const auto *sm_49 = buffer.data(sm + 49);
    const auto *sm_50 = buffer.data(sm + 50);
    const auto *sm_51 = buffer.data(sm + 51);
    const auto *sm_52 = buffer.data(sm + 52);
    const auto *sm_53 = buffer.data(sm + 53);
    const auto *sm_54 = buffer.data(sm + 54);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, sl_0, sl_3, sl_5, \
                         sm_0, sm_2, sm_3, sm_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sl_0[k]
                 + pb_x[k] * sm_0[k];

        t_1[k] = pb_y[k] * sm_0[k];

        t_2[k] = pb_z[k] * sm_0[k];

        t_3[k] = f_1 * sl_3[k]
                 + pb_x[k] * sm_3[k];

        t_4[k] = pb_y[k] * sm_2[k];

        t_5[k] = f_1 * sl_5[k]
                 + pb_x[k] * sm_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, sl_6, sl_9, sl_10, sm_3, \
                         sm_5, sm_6, sm_9, sm_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * sl_6[k]
                 + pb_x[k] * sm_6[k];

        t_7[k] = pb_z[k] * sm_3[k];

        t_8[k] = pb_y[k] * sm_5[k];

        t_9[k] = f_2 * sl_9[k]
                 + pb_x[k] * sm_9[k];

        t_10[k] = f_3 * sl_10[k]
                  + pb_x[k] * sm_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, sl_12, sl_14, sl_15, \
                         sm_6, sm_9, sm_12, sm_14, sm_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * sm_6[k];

        t_12[k] = f_3 * sl_12[k]
                  + pb_x[k] * sm_12[k];

        t_13[k] = pb_y[k] * sm_9[k];

        t_14[k] = f_3 * sl_14[k]
                  + pb_x[k] * sm_14[k];

        t_15[k] = f_4 * sl_15[k]
                  + pb_x[k] * sm_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pb_y, pb_z, sl_17, sl_18, sl_20, \
                         sm_10, sm_14, sm_17, sm_18, sm_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * sm_10[k];

        t_17[k] = f_4 * sl_17[k]
                  + pb_x[k] * sm_17[k];

        t_18[k] = f_4 * sl_18[k]
                  + pb_x[k] * sm_18[k];

        t_19[k] = pb_y[k] * sm_14[k];

        t_20[k] = f_4 * sl_20[k]
                  + pb_x[k] * sm_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, sl_21, sl_23, sl_24, sl_25, \
                         sm_15, sm_21, sm_23, sm_24, sm_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * sl_21[k]
                  + pb_x[k] * sm_21[k];

        t_22[k] = pb_z[k] * sm_15[k];

        t_23[k] = f_5 * sl_23[k]
                  + pb_x[k] * sm_23[k];

        t_24[k] = f_5 * sl_24[k]
                  + pb_x[k] * sm_24[k];

        t_25[k] = f_5 * sl_25[k]
                  + pb_x[k] * sm_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pb_x, pb_y, pb_z, sl_27, sl_28, sl_30, \
                         sm_20, sm_21, sm_27, sm_28, sm_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * sm_20[k];

        t_27[k] = f_5 * sl_27[k]
                  + pb_x[k] * sm_27[k];

        t_28[k] = f_6 * sl_28[k]
                  + pb_x[k] * sm_28[k];

        t_29[k] = pb_z[k] * sm_21[k];

        t_30[k] = f_6 * sl_30[k]
                  + pb_x[k] * sm_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, sl_31, sl_32, sl_33, sl_35, \
                         sm_27, sm_31, sm_32, sm_33, sm_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * sl_31[k]
                  + pb_x[k] * sm_31[k];

        t_32[k] = f_6 * sl_32[k]
                  + pb_x[k] * sm_32[k];

        t_33[k] = f_6 * sl_33[k]
                  + pb_x[k] * sm_33[k];

        t_34[k] = pb_y[k] * sm_27[k];

        t_35[k] = f_6 * sl_35[k]
                  + pb_x[k] * sm_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pb_x, pb_z, sl_36, sl_38, sl_39, sl_40, \
                         sm_28, sm_36, sm_38, sm_39, sm_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_7 * sl_36[k]
                  + pb_x[k] * sm_36[k];

        t_37[k] = pb_z[k] * sm_28[k];

        t_38[k] = f_7 * sl_38[k]
                  + pb_x[k] * sm_38[k];

        t_39[k] = f_7 * sl_39[k]
                  + pb_x[k] * sm_39[k];

        t_40[k] = f_7 * sl_40[k]
                  + pb_x[k] * sm_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pb_x, pb_y, sl_41, sl_42, sl_44, sm_35, \
                         sm_41, sm_42, sm_44, sm_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * sl_41[k]
                  + pb_x[k] * sm_41[k];

        t_42[k] = f_7 * sl_42[k]
                  + pb_x[k] * sm_42[k];

        t_43[k] = pb_y[k] * sm_35[k];

        t_44[k] = f_7 * sl_44[k]
                  + pb_x[k] * sm_44[k];

        t_45[k] = pb_x[k] * sm_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, t_51, t_52, pb_x, sm_46, sm_47, sm_48, \
                         sm_49, sm_50, sm_51, sm_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_x[k] * sm_46[k];

        t_47[k] = pb_x[k] * sm_47[k];

        t_48[k] = pb_x[k] * sm_48[k];

        t_49[k] = pb_x[k] * sm_49[k];

        t_50[k] = pb_x[k] * sm_50[k];

        t_51[k] = pb_x[k] * sm_51[k];

        t_52[k] = pb_x[k] * sm_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pb_x, pb_y, pb_z, sl_36, sl_38, sm_45, \
                         sm_47, sm_53, sm_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_x[k] * sm_53[k];

        t_54[k] = pb_x[k] * sm_54[k];

        t_55[k] = f_0 * sl_36[k]
                  + pb_y[k] * sm_45[k];

        t_56[k] = pb_z[k] * sm_45[k];

        t_57[k] = f_1 * sl_38[k]
                  + pb_y[k] * sm_47[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_y, sl_39, sl_40, sl_41, sl_42, \
                         sl_43, sm_48, sm_49, sm_50, sm_51, sm_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * sl_39[k]
                  + pb_y[k] * sm_48[k];

        t_59[k] = f_3 * sl_40[k]
                  + pb_y[k] * sm_49[k];

        t_60[k] = f_4 * sl_41[k]
                  + pb_y[k] * sm_50[k];

        t_61[k] = f_5 * sl_42[k]
                  + pb_y[k] * sm_51[k];

        t_62[k] = f_6 * sl_43[k]
                  + pb_y[k] * sm_52[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, pb_z, sl_44, sm_53, \
                         sm_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_7 * sl_44[k]
                  + pb_y[k] * sm_53[k];

        t_64[k] = pb_y[k] * sm_54[k];

        t_65[k] = f_0 * sl_44[k]
                  + pb_z[k] * sm_54[k];
    }
}

}  // namespace simdovl
