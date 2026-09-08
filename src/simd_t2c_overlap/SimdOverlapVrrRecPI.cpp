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


#include "SimdOverlapVrrRecPI.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_pi_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sh, const size_t si, const size_t pg,
                          const size_t ph, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_17 = buffer.data(sh + 17);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);
    const auto *pg_20 = buffer.data(pg + 20);
    const auto *pg_21 = buffer.data(pg + 21);
    const auto *pg_22 = buffer.data(pg + 22);
    const auto *pg_23 = buffer.data(pg + 23);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_16 = buffer.data(ph + 16);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_22 = buffer.data(ph + 22);
    const auto *ph_23 = buffer.data(ph + 23);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_25 = buffer.data(ph + 25);
    const auto *ph_26 = buffer.data(ph + 26);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_28 = buffer.data(ph + 28);
    const auto *ph_29 = buffer.data(ph + 29);
    const auto *ph_30 = buffer.data(ph + 30);
    const auto *ph_31 = buffer.data(ph + 31);
    const auto *ph_32 = buffer.data(ph + 32);
    const auto *ph_33 = buffer.data(ph + 33);
    const auto *ph_34 = buffer.data(ph + 34);
    const auto *ph_35 = buffer.data(ph + 35);
    const auto *ph_36 = buffer.data(ph + 36);
    const auto *ph_37 = buffer.data(ph + 37);
    const auto *ph_38 = buffer.data(ph + 38);
    const auto *ph_39 = buffer.data(ph + 39);
    const auto *ph_40 = buffer.data(ph + 40);
    const auto *ph_41 = buffer.data(ph + 41);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pa_x, pb_y, pb_z, sh_0, si_0, pg_0, \
                         ph_0, ph_1, ph_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_0[k]
                 + pa_x[k] * si_0[k];

        t_1[k] = pb_y[k] * ph_0[k];

        t_2[k] = pb_z[k] * ph_0[k];

        t_3[k] = f_1 * pg_0[k]
                 + pb_y[k] * ph_1[k];

        t_4[k] = pb_y[k] * ph_2[k];

        t_5[k] = f_1 * pg_0[k]
                 + pb_z[k] * ph_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_y, pb_z, pg_1, pg_2, pg_3, ph_3, \
                         ph_4, ph_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * pg_1[k]
                 + pb_y[k] * ph_3[k];

        t_7[k] = pb_z[k] * ph_3[k];

        t_8[k] = pb_y[k] * ph_4[k];

        t_9[k] = f_2 * pg_2[k]
                 + pb_z[k] * ph_4[k];

        t_10[k] = f_3 * pg_3[k]
                  + pb_y[k] * ph_5[k];

        t_11[k] = pb_z[k] * ph_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, sh_12, pg_4, ph_6, \
                         ph_7, ph_8, ph_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * pg_4[k]
                  + pb_y[k] * ph_6[k];

        t_13[k] = pb_y[k] * ph_7[k];

        t_14[k] = f_3 * pg_4[k]
                  + pb_z[k] * ph_7[k];

        t_15[k] = f_1 * sh_12[k]
                  + pb_x[k] * ph_10[k];

        t_16[k] = pb_z[k] * ph_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_x, pb_x, pb_y, sh_14, sh_15, sh_17, \
                         si_14, ph_9, ph_11, ph_12, ph_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_1 * sh_14[k]
                  + pb_x[k] * ph_11[k];

        t_18[k] = f_1 * sh_15[k]
                  + pb_x[k] * ph_12[k];

        t_19[k] = pb_y[k] * ph_9[k];

        t_20[k] = f_1 * sh_17[k]
                  + pb_x[k] * ph_13[k];

        t_21[k] = pa_x[k] * si_14[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pa_x, pb_y, pb_z, si_16, si_17, \
                         si_18, si_19, ph_10, ph_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * ph_10[k];

        t_23[k] = pa_x[k] * si_16[k];

        t_24[k] = pa_x[k] * si_17[k];

        t_25[k] = pa_x[k] * si_18[k];

        t_26[k] = pb_y[k] * ph_13[k];

        t_27[k] = pa_x[k] * si_19[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, t_33, pa_y, pb_x, pb_z, si_0, si_4, \
                         pg_7, pg_8, ph_14, ph_15, ph_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * si_0[k];

        t_29[k] = f_4 * pg_7[k]
                  + pb_x[k] * ph_15[k];

        t_30[k] = pb_z[k] * ph_14[k];

        t_31[k] = f_3 * pg_8[k]
                  + pb_x[k] * ph_16[k];

        t_32[k] = pb_z[k] * ph_15[k];

        t_33[k] = pa_y[k] * si_4[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pb_x, pb_z, si_8, pg_9, pg_10, \
                         pg_11, ph_16, ph_17, ph_18, ph_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * pg_9[k]
                  + pb_x[k] * ph_17[k];

        t_35[k] = pb_z[k] * ph_16[k];

        t_36[k] = f_2 * pg_10[k]
                  + pb_x[k] * ph_18[k];

        t_37[k] = pa_y[k] * si_8[k];

        t_38[k] = f_1 * pg_11[k]
                  + pb_x[k] * ph_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_x, pb_z, si_13, pg_13, pg_14, \
                         ph_17, ph_20, ph_21, ph_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_z[k] * ph_17[k];

        t_40[k] = f_1 * pg_13[k]
                  + pb_x[k] * ph_20[k];

        t_41[k] = f_1 * pg_14[k]
                  + pb_x[k] * ph_21[k];

        t_42[k] = pa_y[k] * si_13[k];

        t_43[k] = pb_x[k] * ph_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, pa_y, pb_x, sh_12, si_14, ph_23, \
                         ph_24, ph_25, ph_26, ph_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_x[k] * ph_23[k];

        t_45[k] = pb_x[k] * ph_24[k];

        t_46[k] = pb_x[k] * ph_25[k];

        t_47[k] = pb_x[k] * ph_26[k];

        t_48[k] = pb_x[k] * ph_27[k];

        t_49[k] = f_0 * sh_12[k]
                  + pa_y[k] * si_14[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pb_y, pb_z, sh_17, pg_11, pg_12, pg_13, \
                         ph_22, ph_23, ph_24, ph_25, ph_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_z[k] * ph_22[k];

        t_51[k] = f_1 * pg_11[k]
                  + pb_z[k] * ph_23[k];

        t_52[k] = f_2 * pg_12[k]
                  + pb_z[k] * ph_24[k];

        t_53[k] = f_3 * pg_13[k]
                  + pb_z[k] * ph_25[k];

        t_54[k] = f_1 * sh_17[k]
                  + pb_y[k] * ph_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, t_60, pa_y, pa_z, pb_x, pb_y, si_0, \
                         si_3, si_19, pg_16, ph_28, ph_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * si_19[k];

        t_56[k] = pa_z[k] * si_0[k];

        t_57[k] = pb_y[k] * ph_28[k];

        t_58[k] = f_4 * pg_16[k]
                  + pb_x[k] * ph_29[k];

        t_59[k] = pa_z[k] * si_3[k];

        t_60[k] = pb_y[k] * ph_29[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pa_z, pb_x, pb_y, si_5, pg_17, pg_18, \
                         pg_19, ph_30, ph_31, ph_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pg_17[k]
                  + pb_x[k] * ph_30[k];

        t_62[k] = pa_z[k] * si_5[k];

        t_63[k] = f_2 * pg_18[k]
                  + pb_x[k] * ph_31[k];

        t_64[k] = pb_y[k] * ph_30[k];

        t_65[k] = f_2 * pg_19[k]
                  + pb_x[k] * ph_32[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_z, pb_x, pb_y, si_9, pg_20, pg_21, \
                         pg_23, ph_32, ph_33, ph_34, ph_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * si_9[k];

        t_67[k] = f_1 * pg_20[k]
                  + pb_x[k] * ph_33[k];

        t_68[k] = f_1 * pg_21[k]
                  + pb_x[k] * ph_34[k];

        t_69[k] = pb_y[k] * ph_32[k];

        t_70[k] = f_1 * pg_23[k]
                  + pb_x[k] * ph_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, si_14, ph_36, \
                         ph_37, ph_38, ph_39, ph_40, ph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pb_x[k] * ph_36[k];

        t_72[k] = pb_x[k] * ph_37[k];

        t_73[k] = pb_x[k] * ph_38[k];

        t_74[k] = pb_x[k] * ph_39[k];

        t_75[k] = pb_x[k] * ph_40[k];

        t_76[k] = pb_x[k] * ph_41[k];

        t_77[k] = pa_z[k] * si_14[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pb_y, pg_20, pg_21, pg_22, pg_23, \
                         ph_37, ph_38, ph_39, ph_40, ph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * pg_20[k]
                  + pb_y[k] * ph_37[k];

        t_79[k] = f_3 * pg_21[k]
                  + pb_y[k] * ph_38[k];

        t_80[k] = f_2 * pg_22[k]
                  + pb_y[k] * ph_39[k];

        t_81[k] = f_1 * pg_23[k]
                  + pb_y[k] * ph_40[k];

        t_82[k] = pb_y[k] * ph_41[k];
    }

#pragma omp simd aligned(t_83, pa_z, sh_17, si_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * sh_17[k]
                  + pa_z[k] * si_19[k];
    }
}

auto
compute_prim_pi_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sh, const size_t si, const size_t pg,
                          const size_t ph, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_17 = buffer.data(sh + 17);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_17 = buffer.data(si + 17);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);
    const auto *pg_20 = buffer.data(pg + 20);
    const auto *pg_21 = buffer.data(pg + 21);

    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_22 = buffer.data(ph + 22);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_25 = buffer.data(ph + 25);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_28 = buffer.data(ph + 28);
    const auto *ph_29 = buffer.data(ph + 29);
    const auto *ph_30 = buffer.data(ph + 30);
    const auto *ph_31 = buffer.data(ph + 31);
    const auto *ph_32 = buffer.data(ph + 32);
    const auto *ph_33 = buffer.data(ph + 33);
    const auto *ph_34 = buffer.data(ph + 34);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_y, pb_z, sh_0, si_0, pg_0, pg_1, ph_1, \
                         ph_2, ph_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_0[k]
                 + pa_x[k] * si_0[k];

        t_1[k] = f_1 * pg_0[k]
                 + pb_y[k] * ph_1[k];

        t_2[k] = f_1 * pg_0[k]
                 + pb_z[k] * ph_2[k];

        t_3[k] = f_2 * pg_1[k]
                 + pb_y[k] * ph_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, sh_12, pg_2, pg_3, pg_4, ph_4, \
                         ph_5, ph_6, ph_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pg_2[k]
                 + pb_z[k] * ph_4[k];

        t_5[k] = f_3 * pg_3[k]
                 + pb_y[k] * ph_5[k];

        t_6[k] = f_3 * pg_4[k]
                 + pb_z[k] * ph_6[k];

        t_7[k] = f_1 * sh_12[k]
                 + pb_x[k] * ph_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, pb_x, sh_17, pg_7, pg_8, pg_9, pg_10, \
                         ph_10, ph_11, ph_12, ph_13, ph_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_1 * sh_17[k]
                 + pb_x[k] * ph_10[k];

        t_9[k] = f_4 * pg_7[k]
                 + pb_x[k] * ph_11[k];

        t_10[k] = f_3 * pg_8[k]
                  + pb_x[k] * ph_12[k];

        t_11[k] = f_2 * pg_9[k]
                  + pb_x[k] * ph_13[k];

        t_12[k] = f_1 * pg_10[k]
                  + pb_x[k] * ph_14[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_y, pb_x, pb_z, sh_12, si_12, pg_10, pg_11, \
                         pg_12, ph_15, ph_17, ph_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * pg_12[k]
                  + pb_x[k] * ph_15[k];

        t_14[k] = f_0 * sh_12[k]
                  + pa_y[k] * si_12[k];

        t_15[k] = f_1 * pg_10[k]
                  + pb_z[k] * ph_17[k];

        t_16[k] = f_2 * pg_11[k]
                  + pb_z[k] * ph_18[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, pb_y, pb_z, sh_17, si_0, si_17, \
                         pg_12, ph_19, ph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pg_12[k]
                  + pb_z[k] * ph_19[k];

        t_18[k] = f_1 * sh_17[k]
                  + pb_y[k] * ph_21[k];

        t_19[k] = pa_y[k] * si_17[k];

        t_20[k] = pa_z[k] * si_0[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, pb_x, pb_y, pg_14, pg_15, pg_16, \
                         pg_17, ph_22, ph_24, ph_25, ph_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_4 * pg_14[k]
                  + pb_x[k] * ph_22[k];

        t_22[k] = pb_y[k] * ph_22[k];

        t_23[k] = f_3 * pg_15[k]
                  + pb_x[k] * ph_24[k];

        t_24[k] = f_2 * pg_16[k]
                  + pb_x[k] * ph_25[k];

        t_25[k] = pb_y[k] * ph_24[k];

        t_26[k] = f_2 * pg_17[k]
                  + pb_x[k] * ph_27[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_z, pb_x, pb_y, si_12, pg_18, pg_19, \
                         pg_21, ph_27, ph_28, ph_29, ph_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * pg_18[k]
                  + pb_x[k] * ph_28[k];

        t_28[k] = f_1 * pg_19[k]
                  + pb_x[k] * ph_29[k];

        t_29[k] = pb_y[k] * ph_27[k];

        t_30[k] = f_1 * pg_21[k]
                  + pb_x[k] * ph_30[k];

        t_31[k] = pa_z[k] * si_12[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pg_18, pg_19, pg_20, pg_21, ph_31, \
                         ph_32, ph_33, ph_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * pg_18[k]
                  + pb_y[k] * ph_31[k];

        t_33[k] = f_3 * pg_19[k]
                  + pb_y[k] * ph_32[k];

        t_34[k] = f_2 * pg_20[k]
                  + pb_y[k] * ph_33[k];

        t_35[k] = f_1 * pg_21[k]
                  + pb_y[k] * ph_34[k];
    }

#pragma omp simd aligned(t_36, pa_z, sh_17, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * sh_17[k]
                  + pa_z[k] * si_17[k];
    }
}

auto
compute_prim_pi_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sh, const size_t si, const size_t pg,
                          const size_t ph, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_16 = buffer.data(sh + 16);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_12 = buffer.data(si + 12);

    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);

    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_16 = buffer.data(ph + 16);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_x, sh_0, si_0, pg_3, pg_4, pg_5, ph_5, \
                         ph_6, ph_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_0[k]
                 + pa_x[k] * si_0[k];

        t_1[k] = f_1 * pg_3[k]
                 + pb_x[k] * ph_5[k];

        t_2[k] = f_2 * pg_4[k]
                 + pb_x[k] * ph_6[k];

        t_3[k] = f_3 * pg_5[k]
                 + pb_x[k] * ph_7[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, sh_11, si_8, pg_9, pg_10, pg_12, \
                         ph_13, ph_15, ph_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * sh_11[k]
                 + pa_y[k] * si_8[k];

        t_5[k] = f_1 * pg_9[k]
                 + pb_x[k] * ph_13[k];

        t_6[k] = f_2 * pg_10[k]
                 + pb_x[k] * ph_15[k];

        t_7[k] = f_3 * pg_12[k]
                 + pb_x[k] * ph_16[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_y, pg_12, pg_13, pg_14, ph_17, ph_19, \
                         ph_20, ph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_3 * pg_14[k]
                 + pb_x[k] * ph_17[k];

        t_9[k] = f_1 * pg_12[k]
                 + pb_y[k] * ph_19[k];

        t_10[k] = f_2 * pg_13[k]
                  + pb_y[k] * ph_20[k];

        t_11[k] = f_3 * pg_14[k]
                  + pb_y[k] * ph_21[k];
    }

#pragma omp simd aligned(t_12, pa_z, sh_16, si_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sh_16[k]
                  + pa_z[k] * si_12[k];
    }
}

auto
compute_prim_pi_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sh, const size_t si, const size_t pg,
                          const size_t ph, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_10 = buffer.data(sh + 10);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_4 = buffer.data(si + 4);

    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);

    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_y, sh_5, si_0, pg_7, pg_8, pg_9, ph_6, \
                         ph_7, ph_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_5[k]
                 + pa_y[k] * si_0[k];

        t_1[k] = f_1 * pg_7[k]
                 + pb_y[k] * ph_6[k];

        t_2[k] = f_2 * pg_8[k]
                 + pb_y[k] * ph_7[k];

        t_3[k] = f_3 * pg_9[k]
                 + pb_y[k] * ph_8[k];
    }

#pragma omp simd aligned(t_4, pa_z, sh_10, si_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * sh_10[k]
                 + pa_z[k] * si_4[k];
    }
}

auto
compute_prim_pi_overlap_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t sh, const size_t si, const size_t pg,
                          const size_t ph, const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 0.5 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_6 = buffer.data(sh + 6);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_4 = buffer.data(si + 4);

    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_y, pb_y, sh_1, si_0, pg_2, pg_3, pg_4, ph_2, \
                         ph_3, ph_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_1[k]
                 + pa_y[k] * si_0[k];

        t_1[k] = f_1 * pg_2[k]
                 + pb_y[k] * ph_2[k];

        t_2[k] = f_2 * pg_3[k]
                 + pb_y[k] * ph_3[k];

        t_3[k] = f_3 * pg_4[k]
                 + pb_y[k] * ph_4[k];
    }

#pragma omp simd aligned(t_4, pa_z, sh_6, si_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_0 * sh_6[k]
                 + pa_z[k] * si_4[k];
    }
}

}  // namespace simdovl
