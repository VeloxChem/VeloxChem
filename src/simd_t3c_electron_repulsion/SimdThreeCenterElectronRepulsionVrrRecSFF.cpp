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


#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sff_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sdf0, const size_t sdd,
                                                   const size_t sdf1, const size_t sfp0,
                                                   const size_t sfp1, const size_t sfd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 1.0 / q;

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
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdf0_0 = buffer.data(sdf0 + 0);
    const auto *sdf0_6 = buffer.data(sdf0 + 6);
    const auto *sdf0_9 = buffer.data(sdf0 + 9);
    const auto *sdf0_20 = buffer.data(sdf0 + 20);
    const auto *sdf0_30 = buffer.data(sdf0 + 30);
    const auto *sdf0_36 = buffer.data(sdf0 + 36);
    const auto *sdf0_39 = buffer.data(sdf0 + 39);
    const auto *sdf0_46 = buffer.data(sdf0 + 46);
    const auto *sdf0_49 = buffer.data(sdf0 + 49);
    const auto *sdf0_50 = buffer.data(sdf0 + 50);
    const auto *sdf0_56 = buffer.data(sdf0 + 56);
    const auto *sdf0_59 = buffer.data(sdf0 + 59);

    const auto *sdd_0 = buffer.data(sdd + 0);
    const auto *sdd_3 = buffer.data(sdd + 3);
    const auto *sdd_4 = buffer.data(sdd + 4);
    const auto *sdd_5 = buffer.data(sdd + 5);
    const auto *sdd_6 = buffer.data(sdd + 6);
    const auto *sdd_9 = buffer.data(sdd + 9);
    const auto *sdd_10 = buffer.data(sdd + 10);
    const auto *sdd_11 = buffer.data(sdd + 11);
    const auto *sdd_12 = buffer.data(sdd + 12);
    const auto *sdd_15 = buffer.data(sdd + 15);
    const auto *sdd_16 = buffer.data(sdd + 16);
    const auto *sdd_17 = buffer.data(sdd + 17);
    const auto *sdd_18 = buffer.data(sdd + 18);
    const auto *sdd_21 = buffer.data(sdd + 21);
    const auto *sdd_22 = buffer.data(sdd + 22);
    const auto *sdd_23 = buffer.data(sdd + 23);
    const auto *sdd_24 = buffer.data(sdd + 24);
    const auto *sdd_27 = buffer.data(sdd + 27);
    const auto *sdd_28 = buffer.data(sdd + 28);
    const auto *sdd_29 = buffer.data(sdd + 29);
    const auto *sdd_30 = buffer.data(sdd + 30);
    const auto *sdd_33 = buffer.data(sdd + 33);
    const auto *sdd_34 = buffer.data(sdd + 34);
    const auto *sdd_35 = buffer.data(sdd + 35);

    const auto *sdf1_0 = buffer.data(sdf1 + 0);
    const auto *sdf1_6 = buffer.data(sdf1 + 6);
    const auto *sdf1_9 = buffer.data(sdf1 + 9);
    const auto *sdf1_20 = buffer.data(sdf1 + 20);
    const auto *sdf1_30 = buffer.data(sdf1 + 30);
    const auto *sdf1_36 = buffer.data(sdf1 + 36);
    const auto *sdf1_39 = buffer.data(sdf1 + 39);
    const auto *sdf1_46 = buffer.data(sdf1 + 46);
    const auto *sdf1_49 = buffer.data(sdf1 + 49);
    const auto *sdf1_50 = buffer.data(sdf1 + 50);
    const auto *sdf1_56 = buffer.data(sdf1 + 56);
    const auto *sdf1_59 = buffer.data(sdf1 + 59);

    const auto *sfp0_0 = buffer.data(sfp0 + 0);
    const auto *sfp0_1 = buffer.data(sfp0 + 1);
    const auto *sfp0_2 = buffer.data(sfp0 + 2);
    const auto *sfp0_4 = buffer.data(sfp0 + 4);
    const auto *sfp0_8 = buffer.data(sfp0 + 8);
    const auto *sfp0_18 = buffer.data(sfp0 + 18);
    const auto *sfp0_19 = buffer.data(sfp0 + 19);
    const auto *sfp0_20 = buffer.data(sfp0 + 20);
    const auto *sfp0_23 = buffer.data(sfp0 + 23);
    const auto *sfp0_27 = buffer.data(sfp0 + 27);
    const auto *sfp0_28 = buffer.data(sfp0 + 28);
    const auto *sfp0_29 = buffer.data(sfp0 + 29);

    const auto *sfp1_0 = buffer.data(sfp1 + 0);
    const auto *sfp1_1 = buffer.data(sfp1 + 1);
    const auto *sfp1_2 = buffer.data(sfp1 + 2);
    const auto *sfp1_4 = buffer.data(sfp1 + 4);
    const auto *sfp1_8 = buffer.data(sfp1 + 8);
    const auto *sfp1_18 = buffer.data(sfp1 + 18);
    const auto *sfp1_19 = buffer.data(sfp1 + 19);
    const auto *sfp1_20 = buffer.data(sfp1 + 20);
    const auto *sfp1_23 = buffer.data(sfp1 + 23);
    const auto *sfp1_27 = buffer.data(sfp1 + 27);
    const auto *sfp1_28 = buffer.data(sfp1 + 28);
    const auto *sfp1_29 = buffer.data(sfp1 + 29);

    const auto *sfd_0 = buffer.data(sfd + 0);
    const auto *sfd_3 = buffer.data(sfd + 3);
    const auto *sfd_4 = buffer.data(sfd + 4);
    const auto *sfd_5 = buffer.data(sfd + 5);
    const auto *sfd_6 = buffer.data(sfd + 6);
    const auto *sfd_9 = buffer.data(sfd + 9);
    const auto *sfd_10 = buffer.data(sfd + 10);
    const auto *sfd_11 = buffer.data(sfd + 11);
    const auto *sfd_12 = buffer.data(sfd + 12);
    const auto *sfd_15 = buffer.data(sfd + 15);
    const auto *sfd_16 = buffer.data(sfd + 16);
    const auto *sfd_17 = buffer.data(sfd + 17);
    const auto *sfd_18 = buffer.data(sfd + 18);
    const auto *sfd_21 = buffer.data(sfd + 21);
    const auto *sfd_22 = buffer.data(sfd + 22);
    const auto *sfd_23 = buffer.data(sfd + 23);
    const auto *sfd_24 = buffer.data(sfd + 24);
    const auto *sfd_27 = buffer.data(sfd + 27);
    const auto *sfd_28 = buffer.data(sfd + 28);
    const auto *sfd_29 = buffer.data(sfd + 29);
    const auto *sfd_30 = buffer.data(sfd + 30);
    const auto *sfd_33 = buffer.data(sfd + 33);
    const auto *sfd_34 = buffer.data(sfd + 34);
    const auto *sfd_35 = buffer.data(sfd + 35);
    const auto *sfd_36 = buffer.data(sfd + 36);
    const auto *sfd_39 = buffer.data(sfd + 39);
    const auto *sfd_40 = buffer.data(sfd + 40);
    const auto *sfd_41 = buffer.data(sfd + 41);
    const auto *sfd_42 = buffer.data(sfd + 42);
    const auto *sfd_45 = buffer.data(sfd + 45);
    const auto *sfd_46 = buffer.data(sfd + 46);
    const auto *sfd_47 = buffer.data(sfd + 47);
    const auto *sfd_48 = buffer.data(sfd + 48);
    const auto *sfd_51 = buffer.data(sfd + 51);
    const auto *sfd_52 = buffer.data(sfd + 52);
    const auto *sfd_53 = buffer.data(sfd + 53);
    const auto *sfd_54 = buffer.data(sfd + 54);
    const auto *sfd_57 = buffer.data(sfd + 57);
    const auto *sfd_58 = buffer.data(sfd + 58);
    const auto *sfd_59 = buffer.data(sfd + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, sdd_0, sdd_3, sdd_4, \
                         sfp0_0, sfp1_0, sfd_0, sfd_3, sfd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdd_0[k]
                 + f_1 * sfp0_0[k]
                 - f_2 * sfp1_0[k]
                 + f_3 * pc_x[k] * sfd_0[k];

        t_1[k] = f_3 * pc_y[k] * sfd_0[k];

        t_2[k] = f_3 * pc_z[k] * sfd_0[k];

        t_3[k] = f_0 * sdd_3[k]
                 + f_3 * pc_x[k] * sfd_3[k];

        t_4[k] = f_0 * sdd_4[k]
                 + f_3 * pc_x[k] * sfd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, sdd_5, sfp0_1, sfp0_2, \
                         sfp1_1, sfp1_2, sfd_3, sfd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * sdd_5[k]
                 + f_3 * pc_x[k] * sfd_5[k];

        t_6[k] = f_1 * sfp0_1[k]
                 - f_2 * sfp1_1[k]
                 + f_3 * pc_y[k] * sfd_3[k];

        t_7[k] = f_3 * pc_z[k] * sfd_3[k];

        t_8[k] = f_3 * pc_y[k] * sfd_5[k];

        t_9[k] = f_1 * sfp0_2[k]
                 - f_2 * sfp1_2[k]
                 + f_3 * pc_z[k] * sfd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, sdf0_0, sdd_0, sdd_9, \
                         sdf1_0, sfd_6, sfd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * sdf0_0[k]
                  - f_4 * pc_y[k] * sdf1_0[k];

        t_11[k] = f_5 * sdd_0[k]
                  + f_3 * pc_y[k] * sfd_6[k];

        t_12[k] = f_3 * pc_z[k] * sfd_6[k];

        t_13[k] = f_6 * sdd_9[k]
                  + f_3 * pc_x[k] * sfd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sdd_3, sdd_10, sdd_11, \
                         sfp0_4, sfp1_4, sfd_9, sfd_10, sfd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * sdd_10[k]
                  + f_3 * pc_x[k] * sfd_10[k];

        t_15[k] = f_6 * sdd_11[k]
                  + f_3 * pc_x[k] * sfd_11[k];

        t_16[k] = f_5 * sdd_3[k]
                  + f_1 * sfp0_4[k]
                  - f_2 * sfp1_4[k]
                  + f_3 * pc_y[k] * sfd_9[k];

        t_17[k] = f_3 * pc_z[k] * sfd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, sdf0_0, sdf0_9, \
                         sdd_5, sdf1_0, sdf1_9, sfd_11, sfd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * sdd_5[k]
                  + f_3 * pc_y[k] * sfd_11[k];

        t_19[k] = pb_y[k] * sdf0_9[k]
                  - f_4 * pc_y[k] * sdf1_9[k];

        t_20[k] = pb_z[k] * sdf0_0[k]
                  - f_4 * pc_z[k] * sdf1_0[k];

        t_21[k] = f_3 * pc_y[k] * sfd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, sdd_0, sdd_15, sdd_16, sdd_17, \
                         sfd_12, sfd_15, sfd_16, sfd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * sdd_0[k]
                  + f_3 * pc_z[k] * sfd_12[k];

        t_23[k] = f_6 * sdd_15[k]
                  + f_3 * pc_x[k] * sfd_15[k];

        t_24[k] = f_6 * sdd_16[k]
                  + f_3 * pc_x[k] * sfd_16[k];

        t_25[k] = f_6 * sdd_17[k]
                  + f_3 * pc_x[k] * sfd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, sdf0_6, sdd_3, sdd_5, \
                         sdf1_6, sfp0_8, sfp1_8, sfd_15, sfd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * sdf0_6[k]
                  - f_4 * pc_z[k] * sdf1_6[k];

        t_27[k] = f_5 * sdd_3[k]
                  + f_3 * pc_z[k] * sfd_15[k];

        t_28[k] = f_3 * pc_y[k] * sfd_17[k];

        t_29[k] = f_5 * sdd_5[k]
                  + f_1 * sfp0_8[k]
                  - f_2 * sfp1_8[k]
                  + f_3 * pc_z[k] * sfd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_x, pc_x, pc_y, pc_z, sdf0_30, sdd_6, \
                         sdd_18, sdd_21, sdf1_30, sfd_18, sfd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_x[k] * sdf0_30[k]
                  + f_0 * sdd_18[k]
                  - f_4 * pc_x[k] * sdf1_30[k];

        t_31[k] = f_6 * sdd_6[k]
                  + f_3 * pc_y[k] * sfd_18[k];

        t_32[k] = f_3 * pc_z[k] * sfd_18[k];

        t_33[k] = f_5 * sdd_21[k]
                  + f_3 * pc_x[k] * sfd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, pc_x, pc_z, sdf0_36, sdd_22, sdd_23, \
                         sdf1_36, sfd_21, sfd_22, sfd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * sdd_22[k]
                  + f_3 * pc_x[k] * sfd_22[k];

        t_35[k] = f_5 * sdd_23[k]
                  + f_3 * pc_x[k] * sfd_23[k];

        t_36[k] = pb_x[k] * sdf0_36[k]
                  - f_4 * pc_x[k] * sdf1_36[k];

        t_37[k] = f_3 * pc_z[k] * sfd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, pc_x, pc_y, sdf0_20, sdf0_39, \
                         sdd_11, sdd_12, sdf1_20, sdf1_39, sfd_23, \
                         sfd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * sdd_11[k]
                  + f_3 * pc_y[k] * sfd_23[k];

        t_39[k] = pb_x[k] * sdf0_39[k]
                  - f_4 * pc_x[k] * sdf1_39[k];

        t_40[k] = pb_y[k] * sdf0_20[k]
                  - f_4 * pc_y[k] * sdf1_20[k];

        t_41[k] = f_5 * sdd_12[k]
                  + f_3 * pc_y[k] * sfd_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, sdd_6, sdd_27, sdd_28, sdd_29, \
                         sfd_24, sfd_27, sfd_28, sfd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * sdd_6[k]
                  + f_3 * pc_z[k] * sfd_24[k];

        t_43[k] = f_5 * sdd_27[k]
                  + f_3 * pc_x[k] * sfd_27[k];

        t_44[k] = f_5 * sdd_28[k]
                  + f_3 * pc_x[k] * sfd_28[k];

        t_45[k] = f_5 * sdd_29[k]
                  + f_3 * pc_x[k] * sfd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, pc_x, pc_y, pc_z, sdf0_46, sdf0_49, \
                         sdd_9, sdd_17, sdf1_46, sdf1_49, sfd_27, \
                         sfd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_x[k] * sdf0_46[k]
                  - f_4 * pc_x[k] * sdf1_46[k];

        t_47[k] = f_5 * sdd_9[k]
                  + f_3 * pc_z[k] * sfd_27[k];

        t_48[k] = f_5 * sdd_17[k]
                  + f_3 * pc_y[k] * sfd_29[k];

        t_49[k] = pb_x[k] * sdf0_49[k]
                  - f_4 * pc_x[k] * sdf1_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, pc_x, pc_y, pc_z, sdf0_50, sdd_12, \
                         sdd_30, sdd_33, sdf1_50, sfd_30, sfd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_x[k] * sdf0_50[k]
                  + f_0 * sdd_30[k]
                  - f_4 * pc_x[k] * sdf1_50[k];

        t_51[k] = f_3 * pc_y[k] * sfd_30[k];

        t_52[k] = f_6 * sdd_12[k]
                  + f_3 * pc_z[k] * sfd_30[k];

        t_53[k] = f_5 * sdd_33[k]
                  + f_3 * pc_x[k] * sfd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_x, pc_x, pc_z, sdf0_56, sdd_15, sdd_34, \
                         sdd_35, sdf1_56, sfd_33, sfd_34, sfd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_5 * sdd_34[k]
                  + f_3 * pc_x[k] * sfd_34[k];

        t_55[k] = f_5 * sdd_35[k]
                  + f_3 * pc_x[k] * sfd_35[k];

        t_56[k] = pb_x[k] * sdf0_56[k]
                  - f_4 * pc_x[k] * sdf1_56[k];

        t_57[k] = f_6 * sdd_15[k]
                  + f_3 * pc_z[k] * sfd_33[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pc_x, pc_y, pc_z, sdf0_59, \
                         sdd_18, sdf1_59, sfp0_18, sfp1_18, sfd_35, \
                         sfd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * pc_y[k] * sfd_35[k];

        t_59[k] = pb_x[k] * sdf0_59[k]
                  - f_4 * pc_x[k] * sdf1_59[k];

        t_60[k] = f_1 * sfp0_18[k]
                  - f_2 * sfp1_18[k]
                  + f_3 * pc_x[k] * sfd_36[k];

        t_61[k] = f_0 * sdd_18[k]
                  + f_3 * pc_y[k] * sfd_36[k];

        t_62[k] = f_3 * pc_z[k] * sfd_36[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, sdd_21, sdd_23, \
                         sfp0_19, sfp1_19, sfd_39, sfd_40, sfd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * pc_x[k] * sfd_39[k];

        t_64[k] = f_3 * pc_x[k] * sfd_40[k];

        t_65[k] = f_3 * pc_x[k] * sfd_41[k];

        t_66[k] = f_0 * sdd_21[k]
                  + f_1 * sfp0_19[k]
                  - f_2 * sfp1_19[k]
                  + f_3 * pc_y[k] * sfd_39[k];

        t_67[k] = f_3 * pc_z[k] * sfd_39[k];

        t_68[k] = f_0 * sdd_23[k]
                  + f_3 * pc_y[k] * sfd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_y, pc_z, sdf0_30, sdd_18, sdd_24, \
                         sdf1_30, sfp0_20, sfp1_20, sfd_41, sfd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * sfp0_20[k]
                  - f_2 * sfp1_20[k]
                  + f_3 * pc_z[k] * sfd_41[k];

        t_70[k] = pb_z[k] * sdf0_30[k]
                  - f_4 * pc_z[k] * sdf1_30[k];

        t_71[k] = f_6 * sdd_24[k]
                  + f_3 * pc_y[k] * sfd_42[k];

        t_72[k] = f_5 * sdd_18[k]
                  + f_3 * pc_z[k] * sfd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_z, pc_x, pc_z, sdf0_36, sdd_21, \
                         sdf1_36, sfd_45, sfd_46, sfd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_3 * pc_x[k] * sfd_45[k];

        t_74[k] = f_3 * pc_x[k] * sfd_46[k];

        t_75[k] = f_3 * pc_x[k] * sfd_47[k];

        t_76[k] = pb_z[k] * sdf0_36[k]
                  - f_4 * pc_z[k] * sdf1_36[k];

        t_77[k] = f_5 * sdd_21[k]
                  + f_3 * pc_z[k] * sfd_45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_y, pc_y, pc_z, sdf0_50, sdd_23, sdd_29, \
                         sdd_30, sdf1_50, sfp0_23, sfp1_23, sfd_47, \
                         sfd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_6 * sdd_29[k]
                  + f_3 * pc_y[k] * sfd_47[k];

        t_79[k] = f_5 * sdd_23[k]
                  + f_1 * sfp0_23[k]
                  - f_2 * sfp1_23[k]
                  + f_3 * pc_z[k] * sfd_47[k];

        t_80[k] = pb_y[k] * sdf0_50[k]
                  - f_4 * pc_y[k] * sdf1_50[k];

        t_81[k] = f_5 * sdd_30[k]
                  + f_3 * pc_y[k] * sfd_48[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pc_x, pc_z, sdd_24, sfd_48, sfd_51, sfd_52, \
                         sfd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * sdd_24[k]
                  + f_3 * pc_z[k] * sfd_48[k];

        t_83[k] = f_3 * pc_x[k] * sfd_51[k];

        t_84[k] = f_3 * pc_x[k] * sfd_52[k];

        t_85[k] = f_3 * pc_x[k] * sfd_53[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_y, pc_y, pc_z, sdf0_56, sdf0_59, sdd_27, \
                         sdd_33, sdd_35, sdf1_56, sdf1_59, sfd_51, \
                         sfd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_y[k] * sdf0_56[k]
                  + f_0 * sdd_33[k]
                  - f_4 * pc_y[k] * sdf1_56[k];

        t_87[k] = f_6 * sdd_27[k]
                  + f_3 * pc_z[k] * sfd_51[k];

        t_88[k] = f_5 * sdd_35[k]
                  + f_3 * pc_y[k] * sfd_53[k];

        t_89[k] = pb_y[k] * sdf0_59[k]
                  - f_4 * pc_y[k] * sdf1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, sdd_30, \
                         sfp0_27, sfp1_27, sfd_54, sfd_57, sfd_58, \
                         sfd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * sfp0_27[k]
                  - f_2 * sfp1_27[k]
                  + f_3 * pc_x[k] * sfd_54[k];

        t_91[k] = f_3 * pc_y[k] * sfd_54[k];

        t_92[k] = f_0 * sdd_30[k]
                  + f_3 * pc_z[k] * sfd_54[k];

        t_93[k] = f_3 * pc_x[k] * sfd_57[k];

        t_94[k] = f_3 * pc_x[k] * sfd_58[k];

        t_95[k] = f_3 * pc_x[k] * sfd_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_y, pc_z, sdd_33, sdd_35, sfp0_28, sfp0_29, \
                         sfp1_28, sfp1_29, sfd_57, sfd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * sfp0_28[k]
                  - f_2 * sfp1_28[k]
                  + f_3 * pc_y[k] * sfd_57[k];

        t_97[k] = f_0 * sdd_33[k]
                  + f_3 * pc_z[k] * sfd_57[k];

        t_98[k] = f_3 * pc_y[k] * sfd_59[k];

        t_99[k] = f_0 * sdd_35[k]
                  + f_1 * sfp0_29[k]
                  - f_2 * sfp1_29[k]
                  + f_3 * pc_z[k] * sfd_59[k];
    }
}

}  // namespace simdt3ceri
