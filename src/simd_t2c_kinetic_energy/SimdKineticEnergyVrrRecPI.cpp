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


#include "SimdKineticEnergyVrrRecPI.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_pi_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t sh, const size_t si,
                                 const size_t pi_s, const size_t ph, const size_t ncols,
                                 const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);

    const auto *pi_s_0 = buffer.data(pi_s + 0);
    const auto *pi_s_1 = buffer.data(pi_s + 1);
    const auto *pi_s_2 = buffer.data(pi_s + 2);
    const auto *pi_s_3 = buffer.data(pi_s + 3);
    const auto *pi_s_4 = buffer.data(pi_s + 4);
    const auto *pi_s_5 = buffer.data(pi_s + 5);
    const auto *pi_s_6 = buffer.data(pi_s + 6);
    const auto *pi_s_7 = buffer.data(pi_s + 7);
    const auto *pi_s_8 = buffer.data(pi_s + 8);
    const auto *pi_s_9 = buffer.data(pi_s + 9);
    const auto *pi_s_10 = buffer.data(pi_s + 10);
    const auto *pi_s_11 = buffer.data(pi_s + 11);
    const auto *pi_s_12 = buffer.data(pi_s + 12);
    const auto *pi_s_13 = buffer.data(pi_s + 13);
    const auto *pi_s_14 = buffer.data(pi_s + 14);
    const auto *pi_s_15 = buffer.data(pi_s + 15);
    const auto *pi_s_16 = buffer.data(pi_s + 16);
    const auto *pi_s_17 = buffer.data(pi_s + 17);
    const auto *pi_s_18 = buffer.data(pi_s + 18);
    const auto *pi_s_19 = buffer.data(pi_s + 19);
    const auto *pi_s_20 = buffer.data(pi_s + 20);
    const auto *pi_s_21 = buffer.data(pi_s + 21);
    const auto *pi_s_22 = buffer.data(pi_s + 22);
    const auto *pi_s_23 = buffer.data(pi_s + 23);
    const auto *pi_s_24 = buffer.data(pi_s + 24);
    const auto *pi_s_25 = buffer.data(pi_s + 25);
    const auto *pi_s_26 = buffer.data(pi_s + 26);
    const auto *pi_s_27 = buffer.data(pi_s + 27);
    const auto *pi_s_28 = buffer.data(pi_s + 28);
    const auto *pi_s_29 = buffer.data(pi_s + 29);
    const auto *pi_s_30 = buffer.data(pi_s + 30);
    const auto *pi_s_31 = buffer.data(pi_s + 31);
    const auto *pi_s_32 = buffer.data(pi_s + 32);
    const auto *pi_s_33 = buffer.data(pi_s + 33);
    const auto *pi_s_34 = buffer.data(pi_s + 34);
    const auto *pi_s_35 = buffer.data(pi_s + 35);
    const auto *pi_s_36 = buffer.data(pi_s + 36);
    const auto *pi_s_37 = buffer.data(pi_s + 37);
    const auto *pi_s_38 = buffer.data(pi_s + 38);
    const auto *pi_s_39 = buffer.data(pi_s + 39);
    const auto *pi_s_40 = buffer.data(pi_s + 40);
    const auto *pi_s_41 = buffer.data(pi_s + 41);
    const auto *pi_s_42 = buffer.data(pi_s + 42);
    const auto *pi_s_43 = buffer.data(pi_s + 43);
    const auto *pi_s_44 = buffer.data(pi_s + 44);
    const auto *pi_s_45 = buffer.data(pi_s + 45);
    const auto *pi_s_46 = buffer.data(pi_s + 46);
    const auto *pi_s_47 = buffer.data(pi_s + 47);
    const auto *pi_s_48 = buffer.data(pi_s + 48);
    const auto *pi_s_49 = buffer.data(pi_s + 49);
    const auto *pi_s_50 = buffer.data(pi_s + 50);
    const auto *pi_s_51 = buffer.data(pi_s + 51);
    const auto *pi_s_52 = buffer.data(pi_s + 52);
    const auto *pi_s_53 = buffer.data(pi_s + 53);
    const auto *pi_s_54 = buffer.data(pi_s + 54);
    const auto *pi_s_55 = buffer.data(pi_s + 55);
    const auto *pi_s_56 = buffer.data(pi_s + 56);
    const auto *pi_s_57 = buffer.data(pi_s + 57);
    const auto *pi_s_58 = buffer.data(pi_s + 58);
    const auto *pi_s_59 = buffer.data(pi_s + 59);
    const auto *pi_s_60 = buffer.data(pi_s + 60);
    const auto *pi_s_61 = buffer.data(pi_s + 61);
    const auto *pi_s_62 = buffer.data(pi_s + 62);
    const auto *pi_s_63 = buffer.data(pi_s + 63);
    const auto *pi_s_64 = buffer.data(pi_s + 64);
    const auto *pi_s_65 = buffer.data(pi_s + 65);
    const auto *pi_s_66 = buffer.data(pi_s + 66);
    const auto *pi_s_67 = buffer.data(pi_s + 67);
    const auto *pi_s_68 = buffer.data(pi_s + 68);
    const auto *pi_s_69 = buffer.data(pi_s + 69);
    const auto *pi_s_70 = buffer.data(pi_s + 70);
    const auto *pi_s_71 = buffer.data(pi_s + 71);
    const auto *pi_s_72 = buffer.data(pi_s + 72);
    const auto *pi_s_73 = buffer.data(pi_s + 73);
    const auto *pi_s_74 = buffer.data(pi_s + 74);
    const auto *pi_s_75 = buffer.data(pi_s + 75);
    const auto *pi_s_76 = buffer.data(pi_s + 76);
    const auto *pi_s_77 = buffer.data(pi_s + 77);
    const auto *pi_s_78 = buffer.data(pi_s + 78);
    const auto *pi_s_79 = buffer.data(pi_s + 79);
    const auto *pi_s_80 = buffer.data(pi_s + 80);
    const auto *pi_s_81 = buffer.data(pi_s + 81);
    const auto *pi_s_82 = buffer.data(pi_s + 82);
    const auto *pi_s_83 = buffer.data(pi_s + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pb_y, pb_z, sh_0, sh_3, si_0, si_3, pi_s_0, \
                         pi_s_1, pi_s_2, pi_s_3, ph_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_0[k]
                 + pa_x[k] * si_0[k]
                 + f_1 * pi_s_0[k];

        t_1[k] = f_1 * pi_s_1[k]
                 + pb_y[k] * ph_0[k];

        t_2[k] = f_1 * pi_s_2[k]
                 + pb_z[k] * ph_0[k];

        t_3[k] = f_2 * sh_3[k]
                 + pa_x[k] * si_3[k]
                 + f_1 * pi_s_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, pb_y, sh_4, sh_5, si_4, si_5, pi_s_4, pi_s_5, \
                         pi_s_6, ph_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * pi_s_4[k]
                 + pb_y[k] * ph_1[k];

        t_5[k] = f_2 * sh_4[k]
                 + pa_x[k] * si_4[k]
                 + f_1 * pi_s_5[k];

        t_6[k] = f_3 * sh_5[k]
                 + pa_x[k] * si_5[k]
                 + f_1 * pi_s_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, pb_y, pb_z, sh_8, si_8, pi_s_7, pi_s_8, pi_s_9, \
                         ph_2, ph_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_1 * pi_s_7[k]
                 + pb_z[k] * ph_2[k];

        t_8[k] = f_1 * pi_s_8[k]
                 + pb_y[k] * ph_3[k];

        t_9[k] = f_3 * sh_8[k]
                 + pa_x[k] * si_8[k]
                 + f_1 * pi_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_z, sh_9, sh_10, si_9, si_11, pi_s_10, \
                         pi_s_11, pi_s_12, ph_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * sh_9[k]
                  + pa_x[k] * si_9[k]
                  + f_1 * pi_s_10[k];

        t_11[k] = f_1 * pi_s_11[k]
                  + pb_z[k] * ph_4[k];

        t_12[k] = f_4 * sh_10[k]
                  + pa_x[k] * si_11[k]
                  + f_1 * pi_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_x, pb_y, sh_11, sh_12, si_13, pi_s_13, \
                         pi_s_14, pi_s_15, ph_5, ph_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_1 * pi_s_13[k]
                  + pb_y[k] * ph_5[k];

        t_14[k] = f_4 * sh_11[k]
                  + pa_x[k] * si_13[k]
                  + f_1 * pi_s_14[k];

        t_15[k] = f_5 * sh_12[k]
                  + f_1 * pi_s_15[k]
                  + pb_x[k] * ph_8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_z, sh_14, sh_15, pi_s_16, pi_s_17, \
                         pi_s_18, ph_6, ph_9, ph_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_1 * pi_s_16[k]
                  + pb_z[k] * ph_6[k];

        t_17[k] = f_5 * sh_14[k]
                  + f_1 * pi_s_17[k]
                  + pb_x[k] * ph_9[k];

        t_18[k] = f_5 * sh_15[k]
                  + f_1 * pi_s_18[k]
                  + pb_x[k] * ph_10[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pb_x, pb_y, sh_17, si_14, pi_s_19, pi_s_20, \
                         pi_s_21, ph_7, ph_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * pi_s_19[k]
                  + pb_y[k] * ph_7[k];

        t_20[k] = f_5 * sh_17[k]
                  + f_1 * pi_s_20[k]
                  + pb_x[k] * ph_11[k];

        t_21[k] = pa_x[k] * si_14[k]
                  + f_1 * pi_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_z, si_16, si_17, si_18, pi_s_22, \
                         pi_s_23, pi_s_24, pi_s_25, ph_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_1 * pi_s_22[k]
                  + pb_z[k] * ph_8[k];

        t_23[k] = pa_x[k] * si_16[k]
                  + f_1 * pi_s_23[k];

        t_24[k] = pa_x[k] * si_17[k]
                  + f_1 * pi_s_24[k];

        t_25[k] = pa_x[k] * si_18[k]
                  + f_1 * pi_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pa_y, pb_y, sh_0, si_0, si_1, si_19, \
                         pi_s_26, pi_s_27, pi_s_28, pi_s_29, ph_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * pi_s_26[k]
                  + pb_y[k] * ph_11[k];

        t_27[k] = pa_x[k] * si_19[k]
                  + f_1 * pi_s_27[k];

        t_28[k] = pa_y[k] * si_0[k]
                  + f_1 * pi_s_28[k];

        t_29[k] = f_5 * sh_0[k]
                  + pa_y[k] * si_1[k]
                  + f_1 * pi_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_z, sh_1, si_3, si_4, pi_s_30, \
                         pi_s_31, pi_s_32, pi_s_33, ph_12, ph_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_1 * pi_s_30[k]
                  + pb_z[k] * ph_12[k];

        t_31[k] = f_4 * sh_1[k]
                  + pa_y[k] * si_3[k]
                  + f_1 * pi_s_31[k];

        t_32[k] = f_1 * pi_s_32[k]
                  + pb_z[k] * ph_13[k];

        t_33[k] = pa_y[k] * si_4[k]
                  + f_1 * pi_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pb_z, sh_3, sh_4, si_5, si_7, si_8, \
                         pi_s_34, pi_s_35, pi_s_36, pi_s_37, ph_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * sh_3[k]
                  + pa_y[k] * si_5[k]
                  + f_1 * pi_s_34[k];

        t_35[k] = f_1 * pi_s_35[k]
                  + pb_z[k] * ph_14[k];

        t_36[k] = f_5 * sh_4[k]
                  + pa_y[k] * si_7[k]
                  + f_1 * pi_s_36[k];

        t_37[k] = pa_y[k] * si_8[k]
                  + f_1 * pi_s_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pb_z, sh_5, sh_7, si_9, si_11, pi_s_38, \
                         pi_s_39, pi_s_40, ph_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_2 * sh_5[k]
                  + pa_y[k] * si_9[k]
                  + f_1 * pi_s_38[k];

        t_39[k] = f_1 * pi_s_39[k]
                  + pb_z[k] * ph_15[k];

        t_40[k] = f_4 * sh_7[k]
                  + pa_y[k] * si_11[k]
                  + f_1 * pi_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pb_x, sh_8, si_12, si_13, pi_s_41, \
                         pi_s_42, pi_s_43, pi_s_44, ph_16, ph_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * sh_8[k]
                  + pa_y[k] * si_12[k]
                  + f_1 * pi_s_41[k];

        t_42[k] = pa_y[k] * si_13[k]
                  + f_1 * pi_s_42[k];

        t_43[k] = f_1 * pi_s_43[k]
                  + pb_x[k] * ph_16[k];

        t_44[k] = f_1 * pi_s_44[k]
                  + pb_x[k] * ph_17[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_x, pi_s_45, pi_s_46, pi_s_47, pi_s_48, \
                         ph_18, ph_19, ph_20, ph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_1 * pi_s_45[k]
                  + pb_x[k] * ph_18[k];

        t_46[k] = f_1 * pi_s_46[k]
                  + pb_x[k] * ph_19[k];

        t_47[k] = f_1 * pi_s_47[k]
                  + pb_x[k] * ph_20[k];

        t_48[k] = f_1 * pi_s_48[k]
                  + pb_x[k] * ph_21[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pb_z, sh_12, sh_14, si_14, si_16, pi_s_49, \
                         pi_s_50, pi_s_51, ph_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_0 * sh_12[k]
                  + pa_y[k] * si_14[k]
                  + f_1 * pi_s_49[k];

        t_50[k] = f_1 * pi_s_50[k]
                  + pb_z[k] * ph_16[k];

        t_51[k] = f_2 * sh_14[k]
                  + pa_y[k] * si_16[k]
                  + f_1 * pi_s_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_y, pb_y, sh_15, sh_16, sh_17, si_17, si_18, \
                         pi_s_52, pi_s_53, pi_s_54, ph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * sh_15[k]
                  + pa_y[k] * si_17[k]
                  + f_1 * pi_s_52[k];

        t_53[k] = f_4 * sh_16[k]
                  + pa_y[k] * si_18[k]
                  + f_1 * pi_s_53[k];

        t_54[k] = f_5 * sh_17[k]
                  + f_1 * pi_s_54[k]
                  + pb_y[k] * ph_21[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_y, pa_z, pb_y, sh_0, si_0, si_2, si_19, \
                         pi_s_55, pi_s_56, pi_s_57, pi_s_58, ph_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * si_19[k]
                  + f_1 * pi_s_55[k];

        t_56[k] = pa_z[k] * si_0[k]
                  + f_1 * pi_s_56[k];

        t_57[k] = f_1 * pi_s_57[k]
                  + pb_y[k] * ph_22[k];

        t_58[k] = f_5 * sh_0[k]
                  + pa_z[k] * si_2[k]
                  + f_1 * pi_s_58[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_y, sh_2, si_3, si_4, si_5, pi_s_59, \
                         pi_s_60, pi_s_61, pi_s_62, ph_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_z[k] * si_3[k]
                  + f_1 * pi_s_59[k];

        t_60[k] = f_1 * pi_s_60[k]
                  + pb_y[k] * ph_23[k];

        t_61[k] = f_4 * sh_2[k]
                  + pa_z[k] * si_4[k]
                  + f_1 * pi_s_61[k];

        t_62[k] = pa_z[k] * si_5[k]
                  + f_1 * pi_s_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_y, sh_3, sh_4, si_6, si_8, si_9, \
                         pi_s_63, pi_s_64, pi_s_65, pi_s_66, ph_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_5 * sh_3[k]
                  + pa_z[k] * si_6[k]
                  + f_1 * pi_s_63[k];

        t_64[k] = f_1 * pi_s_64[k]
                  + pb_y[k] * ph_24[k];

        t_65[k] = f_3 * sh_4[k]
                  + pa_z[k] * si_8[k]
                  + f_1 * pi_s_65[k];

        t_66[k] = pa_z[k] * si_9[k]
                  + f_1 * pi_s_66[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pa_z, pb_y, sh_5, sh_6, si_10, si_11, pi_s_67, \
                         pi_s_68, pi_s_69, ph_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_5 * sh_5[k]
                  + pa_z[k] * si_10[k]
                  + f_1 * pi_s_67[k];

        t_68[k] = f_4 * sh_6[k]
                  + pa_z[k] * si_11[k]
                  + f_1 * pi_s_68[k];

        t_69[k] = f_1 * pi_s_69[k]
                  + pb_y[k] * ph_25[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pb_x, sh_8, si_13, pi_s_70, pi_s_71, \
                         pi_s_72, pi_s_73, ph_26, ph_27, ph_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_2 * sh_8[k]
                  + pa_z[k] * si_13[k]
                  + f_1 * pi_s_70[k];

        t_71[k] = f_1 * pi_s_71[k]
                  + pb_x[k] * ph_26[k];

        t_72[k] = f_1 * pi_s_72[k]
                  + pb_x[k] * ph_27[k];

        t_73[k] = f_1 * pi_s_73[k]
                  + pb_x[k] * ph_28[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pb_x, si_14, pi_s_74, pi_s_75, pi_s_76, \
                         pi_s_77, ph_29, ph_30, ph_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_1 * pi_s_74[k]
                  + pb_x[k] * ph_29[k];

        t_75[k] = f_1 * pi_s_75[k]
                  + pb_x[k] * ph_30[k];

        t_76[k] = f_1 * pi_s_76[k]
                  + pb_x[k] * ph_31[k];

        t_77[k] = pa_z[k] * si_14[k]
                  + f_1 * pi_s_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_z, sh_12, sh_13, sh_14, si_15, si_16, si_17, \
                         pi_s_78, pi_s_79, pi_s_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_5 * sh_12[k]
                  + pa_z[k] * si_15[k]
                  + f_1 * pi_s_78[k];

        t_79[k] = f_4 * sh_13[k]
                  + pa_z[k] * si_16[k]
                  + f_1 * pi_s_79[k];

        t_80[k] = f_3 * sh_14[k]
                  + pa_z[k] * si_17[k]
                  + f_1 * pi_s_80[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_y, sh_15, sh_17, si_18, si_19, pi_s_81, \
                         pi_s_82, pi_s_83, ph_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_2 * sh_15[k]
                  + pa_z[k] * si_18[k]
                  + f_1 * pi_s_81[k];

        t_82[k] = f_1 * pi_s_82[k]
                  + pb_y[k] * ph_31[k];

        t_83[k] = f_0 * sh_17[k]
                  + pa_z[k] * si_19[k]
                  + f_1 * pi_s_83[k];
    }
}

auto
compute_prim_pi_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t sh, const size_t si,
                                 const size_t pi_s, const size_t ph, const size_t ncols,
                                 const double alpha, const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / p;

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

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);

    const auto *pi_s_0 = buffer.data(pi_s + 0);
    const auto *pi_s_1 = buffer.data(pi_s + 1);
    const auto *pi_s_2 = buffer.data(pi_s + 2);
    const auto *pi_s_3 = buffer.data(pi_s + 3);
    const auto *pi_s_4 = buffer.data(pi_s + 4);
    const auto *pi_s_5 = buffer.data(pi_s + 5);
    const auto *pi_s_6 = buffer.data(pi_s + 6);
    const auto *pi_s_7 = buffer.data(pi_s + 7);
    const auto *pi_s_8 = buffer.data(pi_s + 8);
    const auto *pi_s_9 = buffer.data(pi_s + 9);
    const auto *pi_s_10 = buffer.data(pi_s + 10);
    const auto *pi_s_11 = buffer.data(pi_s + 11);
    const auto *pi_s_12 = buffer.data(pi_s + 12);
    const auto *pi_s_13 = buffer.data(pi_s + 13);
    const auto *pi_s_14 = buffer.data(pi_s + 14);
    const auto *pi_s_15 = buffer.data(pi_s + 15);
    const auto *pi_s_16 = buffer.data(pi_s + 16);
    const auto *pi_s_17 = buffer.data(pi_s + 17);
    const auto *pi_s_18 = buffer.data(pi_s + 18);
    const auto *pi_s_19 = buffer.data(pi_s + 19);
    const auto *pi_s_20 = buffer.data(pi_s + 20);
    const auto *pi_s_21 = buffer.data(pi_s + 21);
    const auto *pi_s_22 = buffer.data(pi_s + 22);
    const auto *pi_s_23 = buffer.data(pi_s + 23);
    const auto *pi_s_24 = buffer.data(pi_s + 24);
    const auto *pi_s_25 = buffer.data(pi_s + 25);
    const auto *pi_s_26 = buffer.data(pi_s + 26);
    const auto *pi_s_27 = buffer.data(pi_s + 27);
    const auto *pi_s_28 = buffer.data(pi_s + 28);
    const auto *pi_s_29 = buffer.data(pi_s + 29);
    const auto *pi_s_30 = buffer.data(pi_s + 30);
    const auto *pi_s_31 = buffer.data(pi_s + 31);
    const auto *pi_s_32 = buffer.data(pi_s + 32);
    const auto *pi_s_33 = buffer.data(pi_s + 33);
    const auto *pi_s_34 = buffer.data(pi_s + 34);
    const auto *pi_s_35 = buffer.data(pi_s + 35);
    const auto *pi_s_36 = buffer.data(pi_s + 36);

    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_16 = buffer.data(ph + 16);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_21 = buffer.data(ph + 21);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, sh_0, sh_3, sh_4, si_0, si_3, si_4, pi_s_0, \
                         pi_s_1, pi_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_0[k]
                 + pa_x[k] * si_0[k]
                 + f_1 * pi_s_0[k];

        t_1[k] = f_2 * sh_3[k]
                 + pa_x[k] * si_3[k]
                 + f_1 * pi_s_1[k];

        t_2[k] = f_2 * sh_4[k]
                 + pa_x[k] * si_4[k]
                 + f_1 * pi_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, sh_5, sh_8, sh_9, si_5, si_7, si_8, pi_s_3, \
                         pi_s_4, pi_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * sh_5[k]
                 + pa_x[k] * si_5[k]
                 + f_1 * pi_s_3[k];

        t_4[k] = f_3 * sh_8[k]
                 + pa_x[k] * si_7[k]
                 + f_1 * pi_s_4[k];

        t_5[k] = f_4 * sh_9[k]
                 + pa_x[k] * si_8[k]
                 + f_1 * pi_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pb_x, sh_11, sh_12, sh_17, si_11, pi_s_6, \
                         pi_s_7, pi_s_8, ph_3, ph_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_4 * sh_11[k]
                 + pa_x[k] * si_11[k]
                 + f_1 * pi_s_6[k];

        t_7[k] = f_5 * sh_12[k]
                 + f_1 * pi_s_7[k]
                 + pb_x[k] * ph_3[k];

        t_8[k] = f_5 * sh_17[k]
                 + f_1 * pi_s_8[k]
                 + pb_x[k] * ph_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, sh_0, sh_1, sh_3, si_1, si_3, si_5, pi_s_9, \
                         pi_s_10, pi_s_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * sh_0[k]
                 + pa_y[k] * si_1[k]
                 + f_1 * pi_s_9[k];

        t_10[k] = f_4 * sh_1[k]
                  + pa_y[k] * si_3[k]
                  + f_1 * pi_s_10[k];

        t_11[k] = f_3 * sh_3[k]
                  + pa_y[k] * si_5[k]
                  + f_1 * pi_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, sh_5, sh_7, sh_12, si_8, si_10, si_12, \
                         pi_s_12, pi_s_13, pi_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * sh_5[k]
                  + pa_y[k] * si_8[k]
                  + f_1 * pi_s_12[k];

        t_13[k] = f_4 * sh_7[k]
                  + pa_y[k] * si_10[k]
                  + f_1 * pi_s_13[k];

        t_14[k] = f_0 * sh_12[k]
                  + pa_y[k] * si_12[k]
                  + f_1 * pi_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_y, sh_14, sh_15, sh_16, si_14, si_15, si_16, \
                         pi_s_15, pi_s_16, pi_s_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * sh_14[k]
                  + pa_y[k] * si_14[k]
                  + f_1 * pi_s_15[k];

        t_16[k] = f_3 * sh_15[k]
                  + pa_y[k] * si_15[k]
                  + f_1 * pi_s_16[k];

        t_17[k] = f_4 * sh_16[k]
                  + pa_y[k] * si_16[k]
                  + f_1 * pi_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pb_y, sh_17, si_0, si_17, pi_s_18, \
                         pi_s_19, pi_s_20, ph_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * sh_17[k]
                  + f_1 * pi_s_18[k]
                  + pb_y[k] * ph_15[k];

        t_19[k] = pa_y[k] * si_17[k]
                  + f_1 * pi_s_19[k];

        t_20[k] = pa_z[k] * si_0[k]
                  + f_1 * pi_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_z, pb_y, sh_0, sh_2, si_2, si_4, pi_s_21, \
                         pi_s_22, pi_s_23, ph_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_5 * sh_0[k]
                  + pa_z[k] * si_2[k]
                  + f_1 * pi_s_21[k];

        t_22[k] = f_1 * pi_s_22[k]
                  + pb_y[k] * ph_16[k];

        t_23[k] = f_4 * sh_2[k]
                  + pa_z[k] * si_4[k]
                  + f_1 * pi_s_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_z, pb_y, sh_3, sh_4, si_6, si_7, pi_s_24, \
                         pi_s_25, pi_s_26, ph_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * sh_3[k]
                  + pa_z[k] * si_6[k]
                  + f_1 * pi_s_24[k];

        t_25[k] = f_1 * pi_s_25[k]
                  + pb_y[k] * ph_18[k];

        t_26[k] = f_3 * sh_4[k]
                  + pa_z[k] * si_7[k]
                  + f_1 * pi_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_z, pb_y, sh_5, sh_6, si_9, si_10, pi_s_27, \
                         pi_s_28, pi_s_29, ph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_5 * sh_5[k]
                  + pa_z[k] * si_9[k]
                  + f_1 * pi_s_27[k];

        t_28[k] = f_4 * sh_6[k]
                  + pa_z[k] * si_10[k]
                  + f_1 * pi_s_28[k];

        t_29[k] = f_1 * pi_s_29[k]
                  + pb_y[k] * ph_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, sh_8, sh_12, sh_13, si_11, si_12, \
                         si_13, si_14, pi_s_30, pi_s_31, pi_s_32, \
                         pi_s_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * sh_8[k]
                  + pa_z[k] * si_11[k]
                  + f_1 * pi_s_30[k];

        t_31[k] = pa_z[k] * si_12[k]
                  + f_1 * pi_s_31[k];

        t_32[k] = f_5 * sh_12[k]
                  + pa_z[k] * si_13[k]
                  + f_1 * pi_s_32[k];

        t_33[k] = f_4 * sh_13[k]
                  + pa_z[k] * si_14[k]
                  + f_1 * pi_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_z, sh_14, sh_15, sh_17, si_15, si_16, si_17, \
                         pi_s_34, pi_s_35, pi_s_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * sh_14[k]
                  + pa_z[k] * si_15[k]
                  + f_1 * pi_s_34[k];

        t_35[k] = f_2 * sh_15[k]
                  + pa_z[k] * si_16[k]
                  + f_1 * pi_s_35[k];

        t_36[k] = f_0 * sh_17[k]
                  + pa_z[k] * si_17[k]
                  + f_1 * pi_s_36[k];
    }
}

auto
compute_prim_pi_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sh, const size_t si, const size_t pi_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_16 = buffer.data(sh + 16);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);

    const auto *pi_s_0 = buffer.data(pi_s + 0);
    const auto *pi_s_1 = buffer.data(pi_s + 1);
    const auto *pi_s_2 = buffer.data(pi_s + 2);
    const auto *pi_s_3 = buffer.data(pi_s + 3);
    const auto *pi_s_4 = buffer.data(pi_s + 4);
    const auto *pi_s_5 = buffer.data(pi_s + 5);
    const auto *pi_s_6 = buffer.data(pi_s + 6);
    const auto *pi_s_7 = buffer.data(pi_s + 7);
    const auto *pi_s_8 = buffer.data(pi_s + 8);
    const auto *pi_s_9 = buffer.data(pi_s + 9);
    const auto *pi_s_10 = buffer.data(pi_s + 10);
    const auto *pi_s_11 = buffer.data(pi_s + 11);
    const auto *pi_s_12 = buffer.data(pi_s + 12);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, sh_0, sh_1, sh_3, si_0, si_1, si_3, \
                         pi_s_0, pi_s_1, pi_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_0[k]
                 + pa_x[k] * si_0[k]
                 + f_1 * pi_s_0[k];

        t_1[k] = f_2 * sh_1[k]
                 + pa_y[k] * si_1[k]
                 + f_1 * pi_s_1[k];

        t_2[k] = f_3 * sh_3[k]
                 + pa_y[k] * si_3[k]
                 + f_1 * pi_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_y, pa_z, sh_2, sh_5, sh_11, si_2, si_5, si_8, \
                         pi_s_3, pi_s_4, pi_s_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh_5[k]
                 + pa_y[k] * si_5[k]
                 + f_1 * pi_s_3[k];

        t_4[k] = f_0 * sh_11[k]
                 + pa_y[k] * si_8[k]
                 + f_1 * pi_s_4[k];

        t_5[k] = f_2 * sh_2[k]
                 + pa_z[k] * si_2[k]
                 + f_1 * pi_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_z, sh_4, sh_6, sh_7, si_4, si_6, si_7, pi_s_6, \
                         pi_s_7, pi_s_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * sh_4[k]
                 + pa_z[k] * si_4[k]
                 + f_1 * pi_s_6[k];

        t_7[k] = f_2 * sh_6[k]
                 + pa_z[k] * si_6[k]
                 + f_1 * pi_s_7[k];

        t_8[k] = f_4 * sh_7[k]
                 + pa_z[k] * si_7[k]
                 + f_1 * pi_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_z, sh_12, sh_13, sh_14, si_9, si_10, si_11, \
                         pi_s_9, pi_s_10, pi_s_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * sh_12[k]
                 + pa_z[k] * si_9[k]
                 + f_1 * pi_s_9[k];

        t_10[k] = f_3 * sh_13[k]
                  + pa_z[k] * si_10[k]
                  + f_1 * pi_s_10[k];

        t_11[k] = f_4 * sh_14[k]
                  + pa_z[k] * si_11[k]
                  + f_1 * pi_s_11[k];
    }

#pragma omp simd aligned(t_12, pa_z, sh_16, si_12, pi_s_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sh_16[k]
                  + pa_z[k] * si_12[k]
                  + f_1 * pi_s_12[k];
    }
}

auto
compute_prim_pi_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sh, const size_t si, const size_t pi_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_10 = buffer.data(sh + 10);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);

    const auto *pi_s_0 = buffer.data(pi_s + 0);
    const auto *pi_s_1 = buffer.data(pi_s + 1);
    const auto *pi_s_2 = buffer.data(pi_s + 2);
    const auto *pi_s_3 = buffer.data(pi_s + 3);
    const auto *pi_s_4 = buffer.data(pi_s + 4);

#pragma omp simd aligned(t_0, t_1, t_2, pa_y, pa_z, sh_5, sh_6, sh_7, si_0, si_1, si_2, \
                         pi_s_0, pi_s_1, pi_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_5[k]
                 + pa_y[k] * si_0[k]
                 + f_1 * pi_s_0[k];

        t_1[k] = f_2 * sh_6[k]
                 + pa_z[k] * si_1[k]
                 + f_1 * pi_s_1[k];

        t_2[k] = f_3 * sh_7[k]
                 + pa_z[k] * si_2[k]
                 + f_1 * pi_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_z, sh_8, sh_10, si_3, si_4, pi_s_3, \
                         pi_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh_8[k]
                 + pa_z[k] * si_3[k]
                 + f_1 * pi_s_3[k];

        t_4[k] = f_0 * sh_10[k]
                 + pa_z[k] * si_4[k]
                 + f_1 * pi_s_4[k];
    }
}

auto
compute_prim_pi_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t sh, const size_t si, const size_t pi_s,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_6 = buffer.data(sh + 6);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);

    const auto *pi_s_0 = buffer.data(pi_s + 0);
    const auto *pi_s_1 = buffer.data(pi_s + 1);
    const auto *pi_s_2 = buffer.data(pi_s + 2);
    const auto *pi_s_3 = buffer.data(pi_s + 3);
    const auto *pi_s_4 = buffer.data(pi_s + 4);

#pragma omp simd aligned(t_0, t_1, t_2, pa_y, pa_z, sh_1, sh_2, sh_3, si_0, si_1, si_2, \
                         pi_s_0, pi_s_1, pi_s_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh_1[k]
                 + pa_y[k] * si_0[k]
                 + f_1 * pi_s_0[k];

        t_1[k] = f_2 * sh_2[k]
                 + pa_z[k] * si_1[k]
                 + f_1 * pi_s_1[k];

        t_2[k] = f_3 * sh_3[k]
                 + pa_z[k] * si_2[k]
                 + f_1 * pi_s_2[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_z, sh_4, sh_6, si_3, si_4, pi_s_3, \
                         pi_s_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh_4[k]
                 + pa_z[k] * si_3[k]
                 + f_1 * pi_s_3[k];

        t_4[k] = f_0 * sh_6[k]
                 + pa_z[k] * si_4[k]
                 + f_1 * pi_s_4[k];
    }
}

}  // namespace simdkin
