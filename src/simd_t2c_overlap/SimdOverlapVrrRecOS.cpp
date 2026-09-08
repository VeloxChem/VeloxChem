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


#include "SimdOverlapVrrRecOS.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_os_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t ms, const size_t ns, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / p;
    const auto f_1 = 4.0 / p;
    const auto f_2 = 3.5 / p;
    const auto f_3 = 3.0 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 2.0 / p;
    const auto f_6 = 1.5 / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 0.5 / p;

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
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *ms_0 = buffer.data(ms + 0);
    const auto *ms_3 = buffer.data(ms + 3);
    const auto *ms_5 = buffer.data(ms + 5);
    const auto *ms_6 = buffer.data(ms + 6);
    const auto *ms_9 = buffer.data(ms + 9);
    const auto *ms_10 = buffer.data(ms + 10);
    const auto *ms_12 = buffer.data(ms + 12);
    const auto *ms_14 = buffer.data(ms + 14);
    const auto *ms_15 = buffer.data(ms + 15);
    const auto *ms_17 = buffer.data(ms + 17);
    const auto *ms_18 = buffer.data(ms + 18);
    const auto *ms_20 = buffer.data(ms + 20);
    const auto *ms_21 = buffer.data(ms + 21);
    const auto *ms_23 = buffer.data(ms + 23);
    const auto *ms_24 = buffer.data(ms + 24);
    const auto *ms_25 = buffer.data(ms + 25);
    const auto *ms_27 = buffer.data(ms + 27);
    const auto *ms_28 = buffer.data(ms + 28);
    const auto *ms_30 = buffer.data(ms + 30);
    const auto *ms_31 = buffer.data(ms + 31);
    const auto *ms_32 = buffer.data(ms + 32);
    const auto *ms_33 = buffer.data(ms + 33);
    const auto *ms_35 = buffer.data(ms + 35);
    const auto *ms_36 = buffer.data(ms + 36);
    const auto *ms_38 = buffer.data(ms + 38);
    const auto *ms_39 = buffer.data(ms + 39);
    const auto *ms_40 = buffer.data(ms + 40);
    const auto *ms_41 = buffer.data(ms + 41);
    const auto *ms_42 = buffer.data(ms + 42);
    const auto *ms_44 = buffer.data(ms + 44);
    const auto *ms_45 = buffer.data(ms + 45);
    const auto *ms_47 = buffer.data(ms + 47);
    const auto *ms_48 = buffer.data(ms + 48);
    const auto *ms_49 = buffer.data(ms + 49);
    const auto *ms_50 = buffer.data(ms + 50);
    const auto *ms_51 = buffer.data(ms + 51);
    const auto *ms_52 = buffer.data(ms + 52);
    const auto *ms_53 = buffer.data(ms + 53);
    const auto *ms_54 = buffer.data(ms + 54);

    const auto *ns_0 = buffer.data(ns + 0);
    const auto *ns_2 = buffer.data(ns + 2);
    const auto *ns_3 = buffer.data(ns + 3);
    const auto *ns_5 = buffer.data(ns + 5);
    const auto *ns_6 = buffer.data(ns + 6);
    const auto *ns_9 = buffer.data(ns + 9);
    const auto *ns_10 = buffer.data(ns + 10);
    const auto *ns_12 = buffer.data(ns + 12);
    const auto *ns_14 = buffer.data(ns + 14);
    const auto *ns_15 = buffer.data(ns + 15);
    const auto *ns_17 = buffer.data(ns + 17);
    const auto *ns_18 = buffer.data(ns + 18);
    const auto *ns_20 = buffer.data(ns + 20);
    const auto *ns_21 = buffer.data(ns + 21);
    const auto *ns_23 = buffer.data(ns + 23);
    const auto *ns_24 = buffer.data(ns + 24);
    const auto *ns_25 = buffer.data(ns + 25);
    const auto *ns_27 = buffer.data(ns + 27);
    const auto *ns_28 = buffer.data(ns + 28);
    const auto *ns_30 = buffer.data(ns + 30);
    const auto *ns_31 = buffer.data(ns + 31);
    const auto *ns_32 = buffer.data(ns + 32);
    const auto *ns_33 = buffer.data(ns + 33);
    const auto *ns_35 = buffer.data(ns + 35);
    const auto *ns_36 = buffer.data(ns + 36);
    const auto *ns_38 = buffer.data(ns + 38);
    const auto *ns_39 = buffer.data(ns + 39);
    const auto *ns_40 = buffer.data(ns + 40);
    const auto *ns_41 = buffer.data(ns + 41);
    const auto *ns_42 = buffer.data(ns + 42);
    const auto *ns_44 = buffer.data(ns + 44);
    const auto *ns_45 = buffer.data(ns + 45);
    const auto *ns_47 = buffer.data(ns + 47);
    const auto *ns_48 = buffer.data(ns + 48);
    const auto *ns_49 = buffer.data(ns + 49);
    const auto *ns_50 = buffer.data(ns + 50);
    const auto *ns_51 = buffer.data(ns + 51);
    const auto *ns_52 = buffer.data(ns + 52);
    const auto *ns_54 = buffer.data(ns + 54);
    const auto *ns_55 = buffer.data(ns + 55);
    const auto *ns_56 = buffer.data(ns + 56);
    const auto *ns_57 = buffer.data(ns + 57);
    const auto *ns_58 = buffer.data(ns + 58);
    const auto *ns_59 = buffer.data(ns + 59);
    const auto *ns_60 = buffer.data(ns + 60);
    const auto *ns_61 = buffer.data(ns + 61);
    const auto *ns_62 = buffer.data(ns + 62);
    const auto *ns_63 = buffer.data(ns + 63);
    const auto *ns_64 = buffer.data(ns + 64);
    const auto *ns_65 = buffer.data(ns + 65);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pa_y, pa_z, ms_0, ms_3, ms_5, \
                         ns_0, ns_2, ns_3, ns_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ms_0[k]
                 + pa_x[k] * ns_0[k];

        t_1[k] = pa_y[k] * ns_0[k];

        t_2[k] = pa_z[k] * ns_0[k];

        t_3[k] = f_1 * ms_3[k]
                 + pa_x[k] * ns_3[k];

        t_4[k] = pa_y[k] * ns_2[k];

        t_5[k] = f_1 * ms_5[k]
                 + pa_x[k] * ns_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_x, pa_y, pa_z, ms_6, ms_9, ms_10, ns_3, \
                         ns_5, ns_6, ns_9, ns_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * ms_6[k]
                 + pa_x[k] * ns_6[k];

        t_7[k] = pa_z[k] * ns_3[k];

        t_8[k] = pa_y[k] * ns_5[k];

        t_9[k] = f_2 * ms_9[k]
                 + pa_x[k] * ns_9[k];

        t_10[k] = f_3 * ms_10[k]
                  + pa_x[k] * ns_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, ms_12, ms_14, ms_15, \
                         ns_6, ns_9, ns_12, ns_14, ns_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * ns_6[k];

        t_12[k] = f_3 * ms_12[k]
                  + pa_x[k] * ns_12[k];

        t_13[k] = pa_y[k] * ns_9[k];

        t_14[k] = f_3 * ms_14[k]
                  + pa_x[k] * ns_14[k];

        t_15[k] = f_4 * ms_15[k]
                  + pa_x[k] * ns_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pa_y, pa_z, ms_17, ms_18, ms_20, \
                         ns_10, ns_14, ns_17, ns_18, ns_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pa_z[k] * ns_10[k];

        t_17[k] = f_4 * ms_17[k]
                  + pa_x[k] * ns_17[k];

        t_18[k] = f_4 * ms_18[k]
                  + pa_x[k] * ns_18[k];

        t_19[k] = pa_y[k] * ns_14[k];

        t_20[k] = f_4 * ms_20[k]
                  + pa_x[k] * ns_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pa_z, ms_21, ms_23, ms_24, ms_25, \
                         ns_15, ns_21, ns_23, ns_24, ns_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * ms_21[k]
                  + pa_x[k] * ns_21[k];

        t_22[k] = pa_z[k] * ns_15[k];

        t_23[k] = f_5 * ms_23[k]
                  + pa_x[k] * ns_23[k];

        t_24[k] = f_5 * ms_24[k]
                  + pa_x[k] * ns_24[k];

        t_25[k] = f_5 * ms_25[k]
                  + pa_x[k] * ns_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_x, pa_y, pa_z, ms_27, ms_28, ms_30, \
                         ns_20, ns_21, ns_27, ns_28, ns_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_y[k] * ns_20[k];

        t_27[k] = f_5 * ms_27[k]
                  + pa_x[k] * ns_27[k];

        t_28[k] = f_6 * ms_28[k]
                  + pa_x[k] * ns_28[k];

        t_29[k] = pa_z[k] * ns_21[k];

        t_30[k] = f_6 * ms_30[k]
                  + pa_x[k] * ns_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pa_y, ms_31, ms_32, ms_33, ms_35, \
                         ns_27, ns_31, ns_32, ns_33, ns_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * ms_31[k]
                  + pa_x[k] * ns_31[k];

        t_32[k] = f_6 * ms_32[k]
                  + pa_x[k] * ns_32[k];

        t_33[k] = f_6 * ms_33[k]
                  + pa_x[k] * ns_33[k];

        t_34[k] = pa_y[k] * ns_27[k];

        t_35[k] = f_6 * ms_35[k]
                  + pa_x[k] * ns_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_x, pa_z, ms_36, ms_38, ms_39, ms_40, \
                         ns_28, ns_36, ns_38, ns_39, ns_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_7 * ms_36[k]
                  + pa_x[k] * ns_36[k];

        t_37[k] = pa_z[k] * ns_28[k];

        t_38[k] = f_7 * ms_38[k]
                  + pa_x[k] * ns_38[k];

        t_39[k] = f_7 * ms_39[k]
                  + pa_x[k] * ns_39[k];

        t_40[k] = f_7 * ms_40[k]
                  + pa_x[k] * ns_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_x, pa_y, ms_41, ms_42, ms_44, ms_45, \
                         ns_35, ns_41, ns_42, ns_44, ns_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * ms_41[k]
                  + pa_x[k] * ns_41[k];

        t_42[k] = f_7 * ms_42[k]
                  + pa_x[k] * ns_42[k];

        t_43[k] = pa_y[k] * ns_35[k];

        t_44[k] = f_7 * ms_44[k]
                  + pa_x[k] * ns_44[k];

        t_45[k] = f_8 * ms_45[k]
                  + pa_x[k] * ns_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pa_x, pa_z, ms_47, ms_48, ms_49, ms_50, \
                         ns_36, ns_47, ns_48, ns_49, ns_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * ns_36[k];

        t_47[k] = f_8 * ms_47[k]
                  + pa_x[k] * ns_47[k];

        t_48[k] = f_8 * ms_48[k]
                  + pa_x[k] * ns_48[k];

        t_49[k] = f_8 * ms_49[k]
                  + pa_x[k] * ns_49[k];

        t_50[k] = f_8 * ms_50[k]
                  + pa_x[k] * ns_50[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_x, pa_y, ms_51, ms_52, ms_54, ns_44, \
                         ns_51, ns_52, ns_54, ns_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * ms_51[k]
                  + pa_x[k] * ns_51[k];

        t_52[k] = f_8 * ms_52[k]
                  + pa_x[k] * ns_52[k];

        t_53[k] = pa_y[k] * ns_44[k];

        t_54[k] = f_8 * ms_54[k]
                  + pa_x[k] * ns_54[k];

        t_55[k] = pa_x[k] * ns_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, t_62, pa_x, ns_56, ns_57, ns_58, \
                         ns_59, ns_60, ns_61, ns_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_x[k] * ns_56[k];

        t_57[k] = pa_x[k] * ns_57[k];

        t_58[k] = pa_x[k] * ns_58[k];

        t_59[k] = pa_x[k] * ns_59[k];

        t_60[k] = pa_x[k] * ns_60[k];

        t_61[k] = pa_x[k] * ns_61[k];

        t_62[k] = pa_x[k] * ns_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, t_68, pa_x, pa_y, pa_z, ms_45, ms_47, \
                         ns_55, ns_57, ns_63, ns_64, ns_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_x[k] * ns_63[k];

        t_64[k] = pa_x[k] * ns_64[k];

        t_65[k] = pa_x[k] * ns_65[k];

        t_66[k] = f_0 * ms_45[k]
                  + pa_y[k] * ns_55[k];

        t_67[k] = pa_z[k] * ns_55[k];

        t_68[k] = f_1 * ms_47[k]
                  + pa_y[k] * ns_57[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_y, ms_48, ms_49, ms_50, ms_51, \
                         ms_52, ns_58, ns_59, ns_60, ns_61, ns_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_2 * ms_48[k]
                  + pa_y[k] * ns_58[k];

        t_70[k] = f_3 * ms_49[k]
                  + pa_y[k] * ns_59[k];

        t_71[k] = f_4 * ms_50[k]
                  + pa_y[k] * ns_60[k];

        t_72[k] = f_5 * ms_51[k]
                  + pa_y[k] * ns_61[k];

        t_73[k] = f_6 * ms_52[k]
                  + pa_y[k] * ns_62[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_y, pa_z, ms_53, ms_54, ns_63, ns_64, \
                         ns_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_7 * ms_53[k]
                  + pa_y[k] * ns_63[k];

        t_75[k] = f_8 * ms_54[k]
                  + pa_y[k] * ns_64[k];

        t_76[k] = pa_y[k] * ns_65[k];

        t_77[k] = f_0 * ms_54[k]
                  + pa_z[k] * ns_65[k];
    }
}

}  // namespace simdovl
