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


#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sdg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t spg0, const size_t spf,
                                                   const size_t spg1, const size_t sdd0,
                                                   const size_t sdd1, const size_t sdf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;

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
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spg0_0 = buffer.data(spg0 + 0);
    const auto *spg0_3 = buffer.data(spg0 + 3);
    const auto *spg0_5 = buffer.data(spg0 + 5);
    const auto *spg0_18 = buffer.data(spg0 + 18);
    const auto *spg0_25 = buffer.data(spg0 + 25);
    const auto *spg0_27 = buffer.data(spg0 + 27);
    const auto *spg0_29 = buffer.data(spg0 + 29);
    const auto *spg0_30 = buffer.data(spg0 + 30);
    const auto *spg0_35 = buffer.data(spg0 + 35);
    const auto *spg0_40 = buffer.data(spg0 + 40);
    const auto *spg0_42 = buffer.data(spg0 + 42);
    const auto *spg0_44 = buffer.data(spg0 + 44);

    const auto *spf_0 = buffer.data(spf + 0);
    const auto *spf_2 = buffer.data(spf + 2);
    const auto *spf_3 = buffer.data(spf + 3);
    const auto *spf_5 = buffer.data(spf + 5);
    const auto *spf_6 = buffer.data(spf + 6);
    const auto *spf_7 = buffer.data(spf + 7);
    const auto *spf_8 = buffer.data(spf + 8);
    const auto *spf_9 = buffer.data(spf + 9);
    const auto *spf_10 = buffer.data(spf + 10);
    const auto *spf_12 = buffer.data(spf + 12);
    const auto *spf_13 = buffer.data(spf + 13);
    const auto *spf_16 = buffer.data(spf + 16);
    const auto *spf_17 = buffer.data(spf + 17);
    const auto *spf_18 = buffer.data(spf + 18);
    const auto *spf_19 = buffer.data(spf + 19);
    const auto *spf_20 = buffer.data(spf + 20);
    const auto *spf_22 = buffer.data(spf + 22);
    const auto *spf_25 = buffer.data(spf + 25);
    const auto *spf_26 = buffer.data(spf + 26);
    const auto *spf_27 = buffer.data(spf + 27);
    const auto *spf_28 = buffer.data(spf + 28);
    const auto *spf_29 = buffer.data(spf + 29);

    const auto *spg1_0 = buffer.data(spg1 + 0);
    const auto *spg1_3 = buffer.data(spg1 + 3);
    const auto *spg1_5 = buffer.data(spg1 + 5);
    const auto *spg1_18 = buffer.data(spg1 + 18);
    const auto *spg1_25 = buffer.data(spg1 + 25);
    const auto *spg1_27 = buffer.data(spg1 + 27);
    const auto *spg1_29 = buffer.data(spg1 + 29);
    const auto *spg1_30 = buffer.data(spg1 + 30);
    const auto *spg1_35 = buffer.data(spg1 + 35);
    const auto *spg1_40 = buffer.data(spg1 + 40);
    const auto *spg1_42 = buffer.data(spg1 + 42);
    const auto *spg1_44 = buffer.data(spg1 + 44);

    const auto *sdd0_0 = buffer.data(sdd0 + 0);
    const auto *sdd0_3 = buffer.data(sdd0 + 3);
    const auto *sdd0_5 = buffer.data(sdd0 + 5);
    const auto *sdd0_18 = buffer.data(sdd0 + 18);
    const auto *sdd0_21 = buffer.data(sdd0 + 21);
    const auto *sdd0_23 = buffer.data(sdd0 + 23);
    const auto *sdd0_30 = buffer.data(sdd0 + 30);
    const auto *sdd0_33 = buffer.data(sdd0 + 33);
    const auto *sdd0_35 = buffer.data(sdd0 + 35);

    const auto *sdd1_0 = buffer.data(sdd1 + 0);
    const auto *sdd1_3 = buffer.data(sdd1 + 3);
    const auto *sdd1_5 = buffer.data(sdd1 + 5);
    const auto *sdd1_18 = buffer.data(sdd1 + 18);
    const auto *sdd1_21 = buffer.data(sdd1 + 21);
    const auto *sdd1_23 = buffer.data(sdd1 + 23);
    const auto *sdd1_30 = buffer.data(sdd1 + 30);
    const auto *sdd1_33 = buffer.data(sdd1 + 33);
    const auto *sdd1_35 = buffer.data(sdd1 + 35);

    const auto *sdf_0 = buffer.data(sdf + 0);
    const auto *sdf_2 = buffer.data(sdf + 2);
    const auto *sdf_3 = buffer.data(sdf + 3);
    const auto *sdf_5 = buffer.data(sdf + 5);
    const auto *sdf_6 = buffer.data(sdf + 6);
    const auto *sdf_7 = buffer.data(sdf + 7);
    const auto *sdf_8 = buffer.data(sdf + 8);
    const auto *sdf_9 = buffer.data(sdf + 9);
    const auto *sdf_10 = buffer.data(sdf + 10);
    const auto *sdf_12 = buffer.data(sdf + 12);
    const auto *sdf_16 = buffer.data(sdf + 16);
    const auto *sdf_17 = buffer.data(sdf + 17);
    const auto *sdf_18 = buffer.data(sdf + 18);
    const auto *sdf_19 = buffer.data(sdf + 19);
    const auto *sdf_20 = buffer.data(sdf + 20);
    const auto *sdf_22 = buffer.data(sdf + 22);
    const auto *sdf_26 = buffer.data(sdf + 26);
    const auto *sdf_27 = buffer.data(sdf + 27);
    const auto *sdf_28 = buffer.data(sdf + 28);
    const auto *sdf_29 = buffer.data(sdf + 29);
    const auto *sdf_30 = buffer.data(sdf + 30);
    const auto *sdf_32 = buffer.data(sdf + 32);
    const auto *sdf_33 = buffer.data(sdf + 33);
    const auto *sdf_35 = buffer.data(sdf + 35);
    const auto *sdf_36 = buffer.data(sdf + 36);
    const auto *sdf_37 = buffer.data(sdf + 37);
    const auto *sdf_38 = buffer.data(sdf + 38);
    const auto *sdf_39 = buffer.data(sdf + 39);
    const auto *sdf_40 = buffer.data(sdf + 40);
    const auto *sdf_42 = buffer.data(sdf + 42);
    const auto *sdf_46 = buffer.data(sdf + 46);
    const auto *sdf_47 = buffer.data(sdf + 47);
    const auto *sdf_48 = buffer.data(sdf + 48);
    const auto *sdf_49 = buffer.data(sdf + 49);
    const auto *sdf_50 = buffer.data(sdf + 50);
    const auto *sdf_52 = buffer.data(sdf + 52);
    const auto *sdf_53 = buffer.data(sdf + 53);
    const auto *sdf_55 = buffer.data(sdf + 55);
    const auto *sdf_56 = buffer.data(sdf + 56);
    const auto *sdf_57 = buffer.data(sdf + 57);
    const auto *sdf_58 = buffer.data(sdf + 58);
    const auto *sdf_59 = buffer.data(sdf + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, spf_0, spf_3, sdd0_0, sdd0_3, \
                         sdd1_0, sdd1_3, sdf_0, sdf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spf_0[k]
                 + f_1 * sdd0_0[k]
                 - f_2 * sdd1_0[k]
                 + f_3 * pc_x[k] * sdf_0[k];

        t_1[k] = f_3 * pc_y[k] * sdf_0[k];

        t_2[k] = f_3 * pc_z[k] * sdf_0[k];

        t_3[k] = f_0 * spf_3[k]
                 + f_4 * sdd0_3[k]
                 - f_5 * sdd1_3[k]
                 + f_3 * pc_x[k] * sdf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, spf_5, spf_6, spf_7, sdd0_5, sdd1_5, \
                         sdf_2, sdf_5, sdf_6, sdf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sdf_2[k];

        t_5[k] = f_0 * spf_5[k]
                 + f_4 * sdd0_5[k]
                 - f_5 * sdd1_5[k]
                 + f_3 * pc_x[k] * sdf_5[k];

        t_6[k] = f_0 * spf_6[k]
                 + f_3 * pc_x[k] * sdf_6[k];

        t_7[k] = f_0 * spf_7[k]
                 + f_3 * pc_x[k] * sdf_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, spf_8, spf_9, sdd0_3, sdd1_3, \
                         sdf_6, sdf_8, sdf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * spf_8[k]
                 + f_3 * pc_x[k] * sdf_8[k];

        t_9[k] = f_0 * spf_9[k]
                 + f_3 * pc_x[k] * sdf_9[k];

        t_10[k] = f_1 * sdd0_3[k]
                  - f_2 * sdd1_3[k]
                  + f_3 * pc_y[k] * sdf_6[k];

        t_11[k] = f_3 * pc_z[k] * sdf_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, spg0_0, spf_0, \
                         spg1_0, sdd0_5, sdd1_5, sdf_8, sdf_9, sdf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * sdd0_5[k]
                  - f_5 * sdd1_5[k]
                  + f_3 * pc_y[k] * sdf_8[k];

        t_13[k] = f_3 * pc_y[k] * sdf_9[k];

        t_14[k] = f_1 * sdd0_5[k]
                  - f_2 * sdd1_5[k]
                  + f_3 * pc_z[k] * sdf_9[k];

        t_15[k] = pb_y[k] * spg0_0[k]
                  - f_6 * pc_y[k] * spg1_0[k];

        t_16[k] = f_7 * spf_0[k]
                  + f_3 * pc_y[k] * sdf_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_x, pc_x, pc_y, pc_z, spg0_18, spf_2, spf_13, \
                         spg1_18, sdf_10, sdf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * sdf_10[k];

        t_18[k] = pb_x[k] * spg0_18[k]
                  + f_0 * spf_13[k]
                  - f_6 * pc_x[k] * spg1_18[k];

        t_19[k] = f_7 * spf_2[k]
                  + f_3 * pc_y[k] * sdf_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_y, pc_x, pc_y, spg0_5, spf_16, spf_17, \
                         spf_18, spg1_5, sdf_16, sdf_17, sdf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_y[k] * spg0_5[k]
                  - f_6 * pc_y[k] * spg1_5[k];

        t_21[k] = f_7 * spf_16[k]
                  + f_3 * pc_x[k] * sdf_16[k];

        t_22[k] = f_7 * spf_17[k]
                  + f_3 * pc_x[k] * sdf_17[k];

        t_23[k] = f_7 * spf_18[k]
                  + f_3 * pc_x[k] * sdf_18[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, pc_x, pc_z, spg0_25, spg0_27, spf_19, \
                         spg1_25, spg1_27, sdf_16, sdf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * spf_19[k]
                  + f_3 * pc_x[k] * sdf_19[k];

        t_25[k] = pb_x[k] * spg0_25[k]
                  - f_6 * pc_x[k] * spg1_25[k];

        t_26[k] = f_3 * pc_z[k] * sdf_16[k];

        t_27[k] = pb_x[k] * spg0_27[k]
                  - f_6 * pc_x[k] * spg1_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_x, pb_z, pc_x, pc_y, pc_z, spg0_0, \
                         spg0_29, spf_9, spg1_0, spg1_29, sdf_19, \
                         sdf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_7 * spf_9[k]
                  + f_3 * pc_y[k] * sdf_19[k];

        t_29[k] = pb_x[k] * spg0_29[k]
                  - f_6 * pc_x[k] * spg1_29[k];

        t_30[k] = pb_z[k] * spg0_0[k]
                  - f_6 * pc_z[k] * spg1_0[k];

        t_31[k] = f_3 * pc_y[k] * sdf_20[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_z, pc_y, pc_z, spg0_3, spf_0, spg1_3, sdf_20, \
                         sdf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_7 * spf_0[k]
                  + f_3 * pc_z[k] * sdf_20[k];

        t_33[k] = pb_z[k] * spg0_3[k]
                  - f_6 * pc_z[k] * spg1_3[k];

        t_34[k] = f_3 * pc_y[k] * sdf_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pc_x, spg0_35, spf_25, spf_26, spf_27, \
                         spf_28, spg1_35, sdf_26, sdf_27, sdf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_x[k] * spg0_35[k]
                  + f_0 * spf_25[k]
                  - f_6 * pc_x[k] * spg1_35[k];

        t_36[k] = f_7 * spf_26[k]
                  + f_3 * pc_x[k] * sdf_26[k];

        t_37[k] = f_7 * spf_27[k]
                  + f_3 * pc_x[k] * sdf_27[k];

        t_38[k] = f_7 * spf_28[k]
                  + f_3 * pc_x[k] * sdf_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pb_x, pc_x, pc_z, spg0_40, spg0_42, spf_6, \
                         spf_29, spg1_40, spg1_42, sdf_26, sdf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_7 * spf_29[k]
                  + f_3 * pc_x[k] * sdf_29[k];

        t_40[k] = pb_x[k] * spg0_40[k]
                  - f_6 * pc_x[k] * spg1_40[k];

        t_41[k] = f_7 * spf_6[k]
                  + f_3 * pc_z[k] * sdf_26[k];

        t_42[k] = pb_x[k] * spg0_42[k]
                  - f_6 * pc_x[k] * spg1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pc_x, pc_y, pc_z, spg0_44, \
                         spf_10, spg1_44, sdd0_18, sdd1_18, sdf_29, \
                         sdf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * pc_y[k] * sdf_29[k];

        t_44[k] = pb_x[k] * spg0_44[k]
                  - f_6 * pc_x[k] * spg1_44[k];

        t_45[k] = f_1 * sdd0_18[k]
                  - f_2 * sdd1_18[k]
                  + f_3 * pc_x[k] * sdf_30[k];

        t_46[k] = f_0 * spf_10[k]
                  + f_3 * pc_y[k] * sdf_30[k];

        t_47[k] = f_3 * pc_z[k] * sdf_30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pc_x, pc_y, spf_12, sdd0_21, sdd0_23, \
                         sdd1_21, sdd1_23, sdf_32, sdf_33, sdf_35, \
                         sdf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_4 * sdd0_21[k]
                  - f_5 * sdd1_21[k]
                  + f_3 * pc_x[k] * sdf_33[k];

        t_49[k] = f_0 * spf_12[k]
                  + f_3 * pc_y[k] * sdf_32[k];

        t_50[k] = f_4 * sdd0_23[k]
                  - f_5 * sdd1_23[k]
                  + f_3 * pc_x[k] * sdf_35[k];

        t_51[k] = f_3 * pc_x[k] * sdf_36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, spf_16, sdd0_21, \
                         sdd1_21, sdf_36, sdf_37, sdf_38, sdf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_x[k] * sdf_37[k];

        t_53[k] = f_3 * pc_x[k] * sdf_38[k];

        t_54[k] = f_3 * pc_x[k] * sdf_39[k];

        t_55[k] = f_0 * spf_16[k]
                  + f_1 * sdd0_21[k]
                  - f_2 * sdd1_21[k]
                  + f_3 * pc_y[k] * sdf_36[k];

        t_56[k] = f_3 * pc_z[k] * sdf_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, spg0_30, spf_18, spf_19, \
                         spg1_30, sdd0_23, sdd1_23, sdf_38, sdf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_0 * spf_18[k]
                  + f_4 * sdd0_23[k]
                  - f_5 * sdd1_23[k]
                  + f_3 * pc_y[k] * sdf_38[k];

        t_58[k] = f_0 * spf_19[k]
                  + f_3 * pc_y[k] * sdf_39[k];

        t_59[k] = f_1 * sdd0_23[k]
                  - f_2 * sdd1_23[k]
                  + f_3 * pc_z[k] * sdf_39[k];

        t_60[k] = pb_y[k] * spg0_30[k]
                  - f_6 * pc_y[k] * spg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, spg0_18, spf_10, spf_20, \
                         spf_22, spg1_18, sdf_40, sdf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * spf_20[k]
                  + f_3 * pc_y[k] * sdf_40[k];

        t_62[k] = f_7 * spf_10[k]
                  + f_3 * pc_z[k] * sdf_40[k];

        t_63[k] = pb_z[k] * spg0_18[k]
                  - f_6 * pc_z[k] * spg1_18[k];

        t_64[k] = f_7 * spf_22[k]
                  + f_3 * pc_y[k] * sdf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pb_y, pc_x, pc_y, spg0_35, spg1_35, \
                         sdf_46, sdf_47, sdf_48, sdf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * spg0_35[k]
                  - f_6 * pc_y[k] * spg1_35[k];

        t_66[k] = f_3 * pc_x[k] * sdf_46[k];

        t_67[k] = f_3 * pc_x[k] * sdf_47[k];

        t_68[k] = f_3 * pc_x[k] * sdf_48[k];

        t_69[k] = f_3 * pc_x[k] * sdf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_y, pb_z, pc_y, pc_z, spg0_25, spg0_42, spf_16, \
                         spf_28, spg1_25, spg1_42, sdf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * spg0_25[k]
                  - f_6 * pc_z[k] * spg1_25[k];

        t_71[k] = f_7 * spf_16[k]
                  + f_3 * pc_z[k] * sdf_46[k];

        t_72[k] = pb_y[k] * spg0_42[k]
                  + f_0 * spf_28[k]
                  - f_6 * pc_y[k] * spg1_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pb_y, pc_x, pc_y, spg0_44, spf_29, spg1_44, \
                         sdd0_30, sdd1_30, sdf_49, sdf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * spf_29[k]
                  + f_3 * pc_y[k] * sdf_49[k];

        t_74[k] = pb_y[k] * spg0_44[k]
                  - f_6 * pc_y[k] * spg1_44[k];

        t_75[k] = f_1 * sdd0_30[k]
                  - f_2 * sdd1_30[k]
                  + f_3 * pc_x[k] * sdf_50[k];

        t_76[k] = f_3 * pc_y[k] * sdf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pc_x, pc_y, pc_z, spf_20, sdd0_33, sdd0_35, \
                         sdd1_33, sdd1_35, sdf_50, sdf_52, sdf_53, \
                         sdf_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * spf_20[k]
                  + f_3 * pc_z[k] * sdf_50[k];

        t_78[k] = f_4 * sdd0_33[k]
                  - f_5 * sdd1_33[k]
                  + f_3 * pc_x[k] * sdf_53[k];

        t_79[k] = f_3 * pc_y[k] * sdf_52[k];

        t_80[k] = f_4 * sdd0_35[k]
                  - f_5 * sdd1_35[k]
                  + f_3 * pc_x[k] * sdf_55[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, spf_26, \
                         sdd0_33, sdd1_33, sdf_56, sdf_57, sdf_58, \
                         sdf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * pc_x[k] * sdf_56[k];

        t_82[k] = f_3 * pc_x[k] * sdf_57[k];

        t_83[k] = f_3 * pc_x[k] * sdf_58[k];

        t_84[k] = f_3 * pc_x[k] * sdf_59[k];

        t_85[k] = f_1 * sdd0_33[k]
                  - f_2 * sdd1_33[k]
                  + f_3 * pc_y[k] * sdf_56[k];

        t_86[k] = f_0 * spf_26[k]
                  + f_3 * pc_z[k] * sdf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pc_y, pc_z, spf_29, sdd0_35, sdd1_35, sdf_58, \
                         sdf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * sdd0_35[k]
                  - f_5 * sdd1_35[k]
                  + f_3 * pc_y[k] * sdf_58[k];

        t_88[k] = f_3 * pc_y[k] * sdf_59[k];

        t_89[k] = f_0 * spf_29[k]
                  + f_1 * sdd0_35[k]
                  - f_2 * sdd1_35[k]
                  + f_3 * pc_z[k] * sdf_59[k];
    }
}

}  // namespace simdt3ceri
