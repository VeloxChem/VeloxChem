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


#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_psi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ssi0, const size_t ssh,
                                                   const size_t ssi1, const size_t psg0,
                                                   const size_t psg1, const size_t psh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 0.5 / gamma;
    const auto f_4 = 0.5 * p / (gamma * q);
    const auto f_5 = 1.0 / gamma;
    const auto f_6 = p / (gamma * q);
    const auto f_7 = 1.5 / gamma;
    const auto f_8 = 1.5 * p / (gamma * q);
    const auto f_9 = 0.5 / q;
    const auto f_10 = 2.0 / gamma;
    const auto f_11 = 2.0 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssi0_0 = buffer.data(ssi0 + 0);
    const auto *ssi0_3 = buffer.data(ssi0 + 3);
    const auto *ssi0_5 = buffer.data(ssi0 + 5);
    const auto *ssi0_6 = buffer.data(ssi0 + 6);
    const auto *ssi0_9 = buffer.data(ssi0 + 9);
    const auto *ssi0_10 = buffer.data(ssi0 + 10);
    const auto *ssi0_14 = buffer.data(ssi0 + 14);
    const auto *ssi0_21 = buffer.data(ssi0 + 21);
    const auto *ssi0_23 = buffer.data(ssi0 + 23);
    const auto *ssi0_24 = buffer.data(ssi0 + 24);
    const auto *ssi0_25 = buffer.data(ssi0 + 25);
    const auto *ssi0_27 = buffer.data(ssi0 + 27);

    const auto *ssh_0 = buffer.data(ssh + 0);
    const auto *ssh_15 = buffer.data(ssh + 15);
    const auto *ssh_17 = buffer.data(ssh + 17);
    const auto *ssh_18 = buffer.data(ssh + 18);
    const auto *ssh_20 = buffer.data(ssh + 20);

    const auto *ssi1_0 = buffer.data(ssi1 + 0);
    const auto *ssi1_3 = buffer.data(ssi1 + 3);
    const auto *ssi1_5 = buffer.data(ssi1 + 5);
    const auto *ssi1_6 = buffer.data(ssi1 + 6);
    const auto *ssi1_9 = buffer.data(ssi1 + 9);
    const auto *ssi1_10 = buffer.data(ssi1 + 10);
    const auto *ssi1_14 = buffer.data(ssi1 + 14);
    const auto *ssi1_21 = buffer.data(ssi1 + 21);
    const auto *ssi1_23 = buffer.data(ssi1 + 23);
    const auto *ssi1_24 = buffer.data(ssi1 + 24);
    const auto *ssi1_25 = buffer.data(ssi1 + 25);
    const auto *ssi1_27 = buffer.data(ssi1 + 27);

    const auto *psg0_0 = buffer.data(psg0 + 0);
    const auto *psg0_1 = buffer.data(psg0 + 1);
    const auto *psg0_2 = buffer.data(psg0 + 2);
    const auto *psg0_3 = buffer.data(psg0 + 3);
    const auto *psg0_5 = buffer.data(psg0 + 5);
    const auto *psg0_16 = buffer.data(psg0 + 16);
    const auto *psg0_18 = buffer.data(psg0 + 18);
    const auto *psg0_21 = buffer.data(psg0 + 21);
    const auto *psg0_23 = buffer.data(psg0 + 23);
    const auto *psg0_25 = buffer.data(psg0 + 25);
    const auto *psg0_26 = buffer.data(psg0 + 26);
    const auto *psg0_27 = buffer.data(psg0 + 27);
    const auto *psg0_28 = buffer.data(psg0 + 28);
    const auto *psg0_32 = buffer.data(psg0 + 32);
    const auto *psg0_35 = buffer.data(psg0 + 35);
    const auto *psg0_37 = buffer.data(psg0 + 37);
    const auto *psg0_39 = buffer.data(psg0 + 39);
    const auto *psg0_41 = buffer.data(psg0 + 41);
    const auto *psg0_42 = buffer.data(psg0 + 42);
    const auto *psg0_43 = buffer.data(psg0 + 43);
    const auto *psg0_44 = buffer.data(psg0 + 44);

    const auto *psg1_0 = buffer.data(psg1 + 0);
    const auto *psg1_1 = buffer.data(psg1 + 1);
    const auto *psg1_2 = buffer.data(psg1 + 2);
    const auto *psg1_3 = buffer.data(psg1 + 3);
    const auto *psg1_5 = buffer.data(psg1 + 5);
    const auto *psg1_16 = buffer.data(psg1 + 16);
    const auto *psg1_18 = buffer.data(psg1 + 18);
    const auto *psg1_21 = buffer.data(psg1 + 21);
    const auto *psg1_23 = buffer.data(psg1 + 23);
    const auto *psg1_25 = buffer.data(psg1 + 25);
    const auto *psg1_26 = buffer.data(psg1 + 26);
    const auto *psg1_27 = buffer.data(psg1 + 27);
    const auto *psg1_28 = buffer.data(psg1 + 28);
    const auto *psg1_32 = buffer.data(psg1 + 32);
    const auto *psg1_35 = buffer.data(psg1 + 35);
    const auto *psg1_37 = buffer.data(psg1 + 37);
    const auto *psg1_39 = buffer.data(psg1 + 39);
    const auto *psg1_41 = buffer.data(psg1 + 41);
    const auto *psg1_42 = buffer.data(psg1 + 42);
    const auto *psg1_43 = buffer.data(psg1 + 43);
    const auto *psg1_44 = buffer.data(psg1 + 44);

    const auto *psh_0 = buffer.data(psh + 0);
    const auto *psh_1 = buffer.data(psh + 1);
    const auto *psh_2 = buffer.data(psh + 2);
    const auto *psh_3 = buffer.data(psh + 3);
    const auto *psh_5 = buffer.data(psh + 5);
    const auto *psh_6 = buffer.data(psh + 6);
    const auto *psh_8 = buffer.data(psh + 8);
    const auto *psh_9 = buffer.data(psh + 9);
    const auto *psh_10 = buffer.data(psh + 10);
    const auto *psh_14 = buffer.data(psh + 14);
    const auto *psh_15 = buffer.data(psh + 15);
    const auto *psh_17 = buffer.data(psh + 17);
    const auto *psh_18 = buffer.data(psh + 18);
    const auto *psh_20 = buffer.data(psh + 20);
    const auto *psh_21 = buffer.data(psh + 21);
    const auto *psh_22 = buffer.data(psh + 22);
    const auto *psh_24 = buffer.data(psh + 24);
    const auto *psh_27 = buffer.data(psh + 27);
    const auto *psh_29 = buffer.data(psh + 29);
    const auto *psh_31 = buffer.data(psh + 31);
    const auto *psh_33 = buffer.data(psh + 33);
    const auto *psh_34 = buffer.data(psh + 34);
    const auto *psh_36 = buffer.data(psh + 36);
    const auto *psh_37 = buffer.data(psh + 37);
    const auto *psh_38 = buffer.data(psh + 38);
    const auto *psh_39 = buffer.data(psh + 39);
    const auto *psh_40 = buffer.data(psh + 40);
    const auto *psh_41 = buffer.data(psh + 41);
    const auto *psh_42 = buffer.data(psh + 42);
    const auto *psh_44 = buffer.data(psh + 44);
    const auto *psh_47 = buffer.data(psh + 47);
    const auto *psh_49 = buffer.data(psh + 49);
    const auto *psh_51 = buffer.data(psh + 51);
    const auto *psh_53 = buffer.data(psh + 53);
    const auto *psh_54 = buffer.data(psh + 54);
    const auto *psh_56 = buffer.data(psh + 56);
    const auto *psh_57 = buffer.data(psh + 57);
    const auto *psh_58 = buffer.data(psh + 58);
    const auto *psh_59 = buffer.data(psh + 59);
    const auto *psh_60 = buffer.data(psh + 60);
    const auto *psh_61 = buffer.data(psh + 61);
    const auto *psh_62 = buffer.data(psh + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pc_x, pc_y, pc_z, ssi0_0, ssh_0, ssi1_0, \
                         psg0_0, psg1_0, psh_0, psh_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pa_x[k] * ssi0_0[k]
                 + f_0 * ssh_0[k]
                 - f_1 * pc_x[k] * ssi1_0[k];

        t_1[k] = f_2 * pc_y[k] * psh_0[k];

        t_2[k] = f_2 * pc_z[k] * psh_0[k];

        t_3[k] = f_3 * psg0_0[k]
                 - f_4 * psg1_0[k]
                 + f_2 * pc_y[k] * psh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_y, pc_z, psg0_0, psg0_1, psg1_0, psg1_1, \
                         psh_2, psh_3, psh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * psh_2[k];

        t_5[k] = f_3 * psg0_0[k]
                 - f_4 * psg1_0[k]
                 + f_2 * pc_z[k] * psh_2[k];

        t_6[k] = f_5 * psg0_1[k]
                 - f_6 * psg1_1[k]
                 + f_2 * pc_y[k] * psh_3[k];

        t_7[k] = f_2 * pc_z[k] * psh_3[k];

        t_8[k] = f_2 * pc_y[k] * psh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pc_y, pc_z, psg0_2, psg0_3, psg0_5, psg1_2, \
                         psg1_3, psg1_5, psh_5, psh_6, psh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * psg0_2[k]
                 - f_6 * psg1_2[k]
                 + f_2 * pc_z[k] * psh_5[k];

        t_10[k] = f_7 * psg0_3[k]
                  - f_8 * psg1_3[k]
                  + f_2 * pc_y[k] * psh_6[k];

        t_11[k] = f_2 * pc_z[k] * psh_6[k];

        t_12[k] = f_3 * psg0_5[k]
                  - f_4 * psg1_5[k]
                  + f_2 * pc_y[k] * psh_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, ssh_15, ssh_17, \
                         psg0_5, psg1_5, psh_9, psh_10, psh_15, \
                         psh_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * pc_y[k] * psh_9[k];

        t_14[k] = f_7 * psg0_5[k]
                  - f_8 * psg1_5[k]
                  + f_2 * pc_z[k] * psh_9[k];

        t_15[k] = f_9 * ssh_15[k]
                  + f_2 * pc_x[k] * psh_15[k];

        t_16[k] = f_2 * pc_z[k] * psh_10[k];

        t_17[k] = f_9 * ssh_17[k]
                  + f_2 * pc_x[k] * psh_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pc_x, pc_y, ssi0_21, ssh_18, ssh_20, \
                         ssi1_21, psh_14, psh_18, psh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * ssh_18[k]
                  + f_2 * pc_x[k] * psh_18[k];

        t_19[k] = f_2 * pc_y[k] * psh_14[k];

        t_20[k] = f_9 * ssh_20[k]
                  + f_2 * pc_x[k] * psh_20[k];

        t_21[k] = pa_x[k] * ssi0_21[k]
                  - f_1 * pc_x[k] * ssi1_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pc_x, pc_z, ssi0_23, ssi0_24, ssi0_25, \
                         ssi1_23, ssi1_24, ssi1_25, psh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * pc_z[k] * psh_15[k];

        t_23[k] = pa_x[k] * ssi0_23[k]
                  - f_1 * pc_x[k] * ssi1_23[k];

        t_24[k] = pa_x[k] * ssi0_24[k]
                  - f_1 * pc_x[k] * ssi1_24[k];

        t_25[k] = pa_x[k] * ssi0_25[k]
                  - f_1 * pc_x[k] * ssi1_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_x, pa_y, pc_x, pc_y, ssi0_0, ssi0_27, \
                         ssi1_0, ssi1_27, psg0_16, psg1_16, psh_20, \
                         psh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * pc_y[k] * psh_20[k];

        t_27[k] = pa_x[k] * ssi0_27[k]
                  - f_1 * pc_x[k] * ssi1_27[k];

        t_28[k] = pa_y[k] * ssi0_0[k]
                  - f_1 * pc_y[k] * ssi1_0[k];

        t_29[k] = f_10 * psg0_16[k]
                  - f_11 * psg1_16[k]
                  + f_2 * pc_x[k] * psh_22[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_x, pc_y, pc_z, ssi0_5, ssi1_5, \
                         psg0_18, psg1_18, psh_21, psh_22, psh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * pc_z[k] * psh_21[k];

        t_31[k] = f_7 * psg0_18[k]
                  - f_8 * psg1_18[k]
                  + f_2 * pc_x[k] * psh_24[k];

        t_32[k] = f_2 * pc_z[k] * psh_22[k];

        t_33[k] = pa_y[k] * ssi0_5[k]
                  - f_1 * pc_y[k] * ssi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pc_x, pc_z, psg0_21, psg0_23, psg1_21, psg1_23, \
                         psh_24, psh_27, psh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * psg0_21[k]
                  - f_6 * psg1_21[k]
                  + f_2 * pc_x[k] * psh_27[k];

        t_35[k] = f_2 * pc_z[k] * psh_24[k];

        t_36[k] = f_5 * psg0_23[k]
                  - f_6 * psg1_23[k]
                  + f_2 * pc_x[k] * psh_29[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pc_x, pc_y, pc_z, ssi0_9, ssi1_9, psg0_25, \
                         psg1_25, psh_27, psh_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_y[k] * ssi0_9[k]
                  - f_1 * pc_y[k] * ssi1_9[k];

        t_38[k] = f_3 * psg0_25[k]
                  - f_4 * psg1_25[k]
                  + f_2 * pc_x[k] * psh_31[k];

        t_39[k] = f_2 * pc_z[k] * psh_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, ssi0_14, ssi1_14, psg0_27, \
                         psg0_28, psg1_27, psg1_28, psh_33, psh_34, \
                         psh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_3 * psg0_27[k]
                  - f_4 * psg1_27[k]
                  + f_2 * pc_x[k] * psh_33[k];

        t_41[k] = f_3 * psg0_28[k]
                  - f_4 * psg1_28[k]
                  + f_2 * pc_x[k] * psh_34[k];

        t_42[k] = pa_y[k] * ssi0_14[k]
                  - f_1 * pc_y[k] * ssi1_14[k];

        t_43[k] = f_2 * pc_x[k] * psh_36[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pc_x, psh_37, psh_38, psh_39, psh_40, \
                         psh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * pc_x[k] * psh_37[k];

        t_45[k] = f_2 * pc_x[k] * psh_38[k];

        t_46[k] = f_2 * pc_x[k] * psh_39[k];

        t_47[k] = f_2 * pc_x[k] * psh_40[k];

        t_48[k] = f_2 * pc_x[k] * psh_41[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pc_y, pc_z, ssi0_21, ssh_15, ssi1_21, \
                         psg0_25, psg1_25, psh_36, psh_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_y[k] * ssi0_21[k]
                  + f_0 * ssh_15[k]
                  - f_1 * pc_y[k] * ssi1_21[k];

        t_50[k] = f_2 * pc_z[k] * psh_36[k];

        t_51[k] = f_3 * psg0_25[k]
                  - f_4 * psg1_25[k]
                  + f_2 * pc_z[k] * psh_37[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pc_y, pc_z, ssh_20, psg0_26, psg0_27, psg1_26, \
                         psg1_27, psh_38, psh_39, psh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_5 * psg0_26[k]
                  - f_6 * psg1_26[k]
                  + f_2 * pc_z[k] * psh_38[k];

        t_53[k] = f_7 * psg0_27[k]
                  - f_8 * psg1_27[k]
                  + f_2 * pc_z[k] * psh_39[k];

        t_54[k] = f_9 * ssh_20[k]
                  + f_2 * pc_y[k] * psh_41[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, ssi0_0, ssi0_27, ssi1_0, \
                         ssi1_27, psh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * ssi0_27[k]
                  - f_1 * pc_y[k] * ssi1_27[k];

        t_56[k] = pa_z[k] * ssi0_0[k]
                  - f_1 * pc_z[k] * ssi1_0[k];

        t_57[k] = f_2 * pc_y[k] * psh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_x, pc_y, pc_z, ssi0_3, ssi1_3, \
                         psg0_32, psg0_35, psg1_32, psg1_35, psh_44, \
                         psh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_10 * psg0_32[k]
                  - f_11 * psg1_32[k]
                  + f_2 * pc_x[k] * psh_44[k];

        t_59[k] = pa_z[k] * ssi0_3[k]
                  - f_1 * pc_z[k] * ssi1_3[k];

        t_60[k] = f_2 * pc_y[k] * psh_44[k];

        t_61[k] = f_7 * psg0_35[k]
                  - f_8 * psg1_35[k]
                  + f_2 * pc_x[k] * psh_47[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_z, pc_x, pc_y, pc_z, ssi0_6, ssi1_6, psg0_37, \
                         psg1_37, psh_47, psh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * ssi0_6[k]
                  - f_1 * pc_z[k] * ssi1_6[k];

        t_63[k] = f_5 * psg0_37[k]
                  - f_6 * psg1_37[k]
                  + f_2 * pc_x[k] * psh_49[k];

        t_64[k] = f_2 * pc_y[k] * psh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_z, pc_x, pc_z, ssi0_10, ssi1_10, psg0_39, \
                         psg0_41, psg1_39, psg1_41, psh_51, psh_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_5 * psg0_39[k]
                  - f_6 * psg1_39[k]
                  + f_2 * pc_x[k] * psh_51[k];

        t_66[k] = pa_z[k] * ssi0_10[k]
                  - f_1 * pc_z[k] * ssi1_10[k];

        t_67[k] = f_3 * psg0_41[k]
                  - f_4 * psg1_41[k]
                  + f_2 * pc_x[k] * psh_53[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pc_x, pc_y, psg0_42, psg0_44, psg1_42, \
                         psg1_44, psh_51, psh_54, psh_56, psh_57, \
                         psh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * psg0_42[k]
                  - f_4 * psg1_42[k]
                  + f_2 * pc_x[k] * psh_54[k];

        t_69[k] = f_2 * pc_y[k] * psh_51[k];

        t_70[k] = f_3 * psg0_44[k]
                  - f_4 * psg1_44[k]
                  + f_2 * pc_x[k] * psh_56[k];

        t_71[k] = f_2 * pc_x[k] * psh_57[k];

        t_72[k] = f_2 * pc_x[k] * psh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pc_x, pc_z, ssi0_21, ssi1_21, \
                         psh_59, psh_60, psh_61, psh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_2 * pc_x[k] * psh_59[k];

        t_74[k] = f_2 * pc_x[k] * psh_60[k];

        t_75[k] = f_2 * pc_x[k] * psh_61[k];

        t_76[k] = f_2 * pc_x[k] * psh_62[k];

        t_77[k] = pa_z[k] * ssi0_21[k]
                  - f_1 * pc_z[k] * ssi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, psg0_41, psg0_42, psg0_43, psg1_41, psg1_42, \
                         psg1_43, psh_58, psh_59, psh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * psg0_41[k]
                  - f_11 * psg1_41[k]
                  + f_2 * pc_y[k] * psh_58[k];

        t_79[k] = f_7 * psg0_42[k]
                  - f_8 * psg1_42[k]
                  + f_2 * pc_y[k] * psh_59[k];

        t_80[k] = f_5 * psg0_43[k]
                  - f_6 * psg1_43[k]
                  + f_2 * pc_y[k] * psh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pc_y, pc_z, ssi0_27, ssh_20, ssi1_27, \
                         psg0_44, psg1_44, psh_61, psh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * psg0_44[k]
                  - f_4 * psg1_44[k]
                  + f_2 * pc_y[k] * psh_61[k];

        t_82[k] = f_2 * pc_y[k] * psh_62[k];

        t_83[k] = pa_z[k] * ssi0_27[k]
                  + f_0 * ssh_20[k]
                  - f_1 * pc_z[k] * ssi1_27[k];
    }
}

}  // namespace simdt3ceri
