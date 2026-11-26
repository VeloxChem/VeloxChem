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


#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sgd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sfd0, const size_t sfp,
                                                   const size_t sfd1, const size_t sgs0,
                                                   const size_t sgs1, const size_t sgp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 1.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 1.0 / q;

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

    const auto *sfd0_0 = buffer.data(sfd0 + 0);
    const auto *sfd0_3 = buffer.data(sfd0 + 3);
    const auto *sfd0_5 = buffer.data(sfd0 + 5);
    const auto *sfd0_9 = buffer.data(sfd0 + 9);
    const auto *sfd0_12 = buffer.data(sfd0 + 12);
    const auto *sfd0_17 = buffer.data(sfd0 + 17);
    const auto *sfd0_18 = buffer.data(sfd0 + 18);
    const auto *sfd0_30 = buffer.data(sfd0 + 30);
    const auto *sfd0_36 = buffer.data(sfd0 + 36);
    const auto *sfd0_39 = buffer.data(sfd0 + 39);
    const auto *sfd0_41 = buffer.data(sfd0 + 41);
    const auto *sfd0_45 = buffer.data(sfd0 + 45);
    const auto *sfd0_47 = buffer.data(sfd0 + 47);
    const auto *sfd0_51 = buffer.data(sfd0 + 51);
    const auto *sfd0_53 = buffer.data(sfd0 + 53);
    const auto *sfd0_54 = buffer.data(sfd0 + 54);
    const auto *sfd0_57 = buffer.data(sfd0 + 57);
    const auto *sfd0_59 = buffer.data(sfd0 + 59);

    const auto *sfp_0 = buffer.data(sfp + 0);
    const auto *sfp_1 = buffer.data(sfp + 1);
    const auto *sfp_2 = buffer.data(sfp + 2);
    const auto *sfp_4 = buffer.data(sfp + 4);
    const auto *sfp_5 = buffer.data(sfp + 5);
    const auto *sfp_7 = buffer.data(sfp + 7);
    const auto *sfp_8 = buffer.data(sfp + 8);
    const auto *sfp_9 = buffer.data(sfp + 9);
    const auto *sfp_10 = buffer.data(sfp + 10);
    const auto *sfp_11 = buffer.data(sfp + 11);
    const auto *sfp_13 = buffer.data(sfp + 13);
    const auto *sfp_14 = buffer.data(sfp + 14);
    const auto *sfp_15 = buffer.data(sfp + 15);
    const auto *sfp_16 = buffer.data(sfp + 16);
    const auto *sfp_17 = buffer.data(sfp + 17);
    const auto *sfp_18 = buffer.data(sfp + 18);
    const auto *sfp_19 = buffer.data(sfp + 19);
    const auto *sfp_20 = buffer.data(sfp + 20);
    const auto *sfp_22 = buffer.data(sfp + 22);
    const auto *sfp_23 = buffer.data(sfp + 23);
    const auto *sfp_25 = buffer.data(sfp + 25);
    const auto *sfp_26 = buffer.data(sfp + 26);
    const auto *sfp_27 = buffer.data(sfp + 27);
    const auto *sfp_28 = buffer.data(sfp + 28);
    const auto *sfp_29 = buffer.data(sfp + 29);

    const auto *sfd1_0 = buffer.data(sfd1 + 0);
    const auto *sfd1_3 = buffer.data(sfd1 + 3);
    const auto *sfd1_5 = buffer.data(sfd1 + 5);
    const auto *sfd1_9 = buffer.data(sfd1 + 9);
    const auto *sfd1_12 = buffer.data(sfd1 + 12);
    const auto *sfd1_17 = buffer.data(sfd1 + 17);
    const auto *sfd1_18 = buffer.data(sfd1 + 18);
    const auto *sfd1_30 = buffer.data(sfd1 + 30);
    const auto *sfd1_36 = buffer.data(sfd1 + 36);
    const auto *sfd1_39 = buffer.data(sfd1 + 39);
    const auto *sfd1_41 = buffer.data(sfd1 + 41);
    const auto *sfd1_45 = buffer.data(sfd1 + 45);
    const auto *sfd1_47 = buffer.data(sfd1 + 47);
    const auto *sfd1_51 = buffer.data(sfd1 + 51);
    const auto *sfd1_53 = buffer.data(sfd1 + 53);
    const auto *sfd1_54 = buffer.data(sfd1 + 54);
    const auto *sfd1_57 = buffer.data(sfd1 + 57);
    const auto *sfd1_59 = buffer.data(sfd1 + 59);

    const auto *sgs0_0 = buffer.data(sgs0 + 0);
    const auto *sgs0_1 = buffer.data(sgs0 + 1);
    const auto *sgs0_2 = buffer.data(sgs0 + 2);
    const auto *sgs0_3 = buffer.data(sgs0 + 3);
    const auto *sgs0_5 = buffer.data(sgs0 + 5);
    const auto *sgs0_10 = buffer.data(sgs0 + 10);
    const auto *sgs0_11 = buffer.data(sgs0 + 11);
    const auto *sgs0_12 = buffer.data(sgs0 + 12);
    const auto *sgs0_14 = buffer.data(sgs0 + 14);

    const auto *sgs1_0 = buffer.data(sgs1 + 0);
    const auto *sgs1_1 = buffer.data(sgs1 + 1);
    const auto *sgs1_2 = buffer.data(sgs1 + 2);
    const auto *sgs1_3 = buffer.data(sgs1 + 3);
    const auto *sgs1_5 = buffer.data(sgs1 + 5);
    const auto *sgs1_10 = buffer.data(sgs1 + 10);
    const auto *sgs1_11 = buffer.data(sgs1 + 11);
    const auto *sgs1_12 = buffer.data(sgs1 + 12);
    const auto *sgs1_14 = buffer.data(sgs1 + 14);

    const auto *sgp_0 = buffer.data(sgp + 0);
    const auto *sgp_1 = buffer.data(sgp + 1);
    const auto *sgp_2 = buffer.data(sgp + 2);
    const auto *sgp_4 = buffer.data(sgp + 4);
    const auto *sgp_5 = buffer.data(sgp + 5);
    const auto *sgp_7 = buffer.data(sgp + 7);
    const auto *sgp_8 = buffer.data(sgp + 8);
    const auto *sgp_9 = buffer.data(sgp + 9);
    const auto *sgp_10 = buffer.data(sgp + 10);
    const auto *sgp_11 = buffer.data(sgp + 11);
    const auto *sgp_13 = buffer.data(sgp + 13);
    const auto *sgp_14 = buffer.data(sgp + 14);
    const auto *sgp_15 = buffer.data(sgp + 15);
    const auto *sgp_16 = buffer.data(sgp + 16);
    const auto *sgp_17 = buffer.data(sgp + 17);
    const auto *sgp_19 = buffer.data(sgp + 19);
    const auto *sgp_20 = buffer.data(sgp + 20);
    const auto *sgp_22 = buffer.data(sgp + 22);
    const auto *sgp_23 = buffer.data(sgp + 23);
    const auto *sgp_25 = buffer.data(sgp + 25);
    const auto *sgp_26 = buffer.data(sgp + 26);
    const auto *sgp_28 = buffer.data(sgp + 28);
    const auto *sgp_29 = buffer.data(sgp + 29);
    const auto *sgp_30 = buffer.data(sgp + 30);
    const auto *sgp_31 = buffer.data(sgp + 31);
    const auto *sgp_32 = buffer.data(sgp + 32);
    const auto *sgp_34 = buffer.data(sgp + 34);
    const auto *sgp_35 = buffer.data(sgp + 35);
    const auto *sgp_36 = buffer.data(sgp + 36);
    const auto *sgp_37 = buffer.data(sgp + 37);
    const auto *sgp_38 = buffer.data(sgp + 38);
    const auto *sgp_40 = buffer.data(sgp + 40);
    const auto *sgp_41 = buffer.data(sgp + 41);
    const auto *sgp_42 = buffer.data(sgp + 42);
    const auto *sgp_43 = buffer.data(sgp + 43);
    const auto *sgp_44 = buffer.data(sgp + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, sfp_0, sfp_1, sfp_2, sgs0_0, \
                         sgs1_0, sgp_0, sgp_1, sgp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfp_0[k]
                 + f_1 * sgs0_0[k]
                 - f_2 * sgs1_0[k]
                 + f_3 * pc_x[k] * sgp_0[k];

        t_1[k] = f_0 * sfp_1[k]
                 + f_3 * pc_x[k] * sgp_1[k];

        t_2[k] = f_0 * sfp_2[k]
                 + f_3 * pc_x[k] * sgp_2[k];

        t_3[k] = f_1 * sgs0_0[k]
                 - f_2 * sgs1_0[k]
                 + f_3 * pc_y[k] * sgp_1[k];

        t_4[k] = f_3 * pc_y[k] * sgp_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, sfd0_0, sfp_4, sfd1_0, sgs0_0, \
                         sgs1_0, sgp_2, sgp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sgs0_0[k]
                 - f_2 * sgs1_0[k]
                 + f_3 * pc_z[k] * sgp_2[k];

        t_6[k] = pb_y[k] * sfd0_0[k]
                 - f_4 * pc_y[k] * sfd1_0[k];

        t_7[k] = f_5 * sfp_4[k]
                 + f_3 * pc_x[k] * sgp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, sfd0_5, sfp_1, sfp_2, sfp_5, \
                         sfd1_5, sgs0_1, sgs1_1, sgp_4, sgp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * sfp_5[k]
                 + f_3 * pc_x[k] * sgp_5[k];

        t_9[k] = f_6 * sfp_1[k]
                 + f_1 * sgs0_1[k]
                 - f_2 * sgs1_1[k]
                 + f_3 * pc_y[k] * sgp_4[k];

        t_10[k] = f_6 * sfp_2[k]
                  + f_3 * pc_y[k] * sgp_5[k];

        t_11[k] = pb_y[k] * sfd0_5[k]
                  - f_4 * pc_y[k] * sfd1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, sfd0_0, sfd0_3, sfp_7, \
                         sfp_8, sfd1_0, sfd1_3, sgp_7, sgp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * sfd0_0[k]
                  - f_4 * pc_z[k] * sfd1_0[k];

        t_13[k] = f_5 * sfp_7[k]
                  + f_3 * pc_x[k] * sgp_7[k];

        t_14[k] = f_5 * sfp_8[k]
                  + f_3 * pc_x[k] * sgp_8[k];

        t_15[k] = pb_z[k] * sfd0_3[k]
                  - f_4 * pc_z[k] * sfd1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_y, pc_z, sfp_2, sfp_9, sgs0_2, sgs0_3, \
                         sgs1_2, sgs1_3, sgp_8, sgp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * sgp_8[k];

        t_17[k] = f_6 * sfp_2[k]
                  + f_1 * sgs0_2[k]
                  - f_2 * sgs1_2[k]
                  + f_3 * pc_z[k] * sgp_8[k];

        t_18[k] = f_7 * sfp_9[k]
                  + f_1 * sgs0_3[k]
                  - f_2 * sgs1_3[k]
                  + f_3 * pc_x[k] * sgp_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, sfp_4, sfp_5, sfp_10, \
                         sfp_11, sgs0_3, sgs1_3, sgp_10, sgp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * sfp_10[k]
                  + f_3 * pc_x[k] * sgp_10[k];

        t_20[k] = f_7 * sfp_11[k]
                  + f_3 * pc_x[k] * sgp_11[k];

        t_21[k] = f_7 * sfp_4[k]
                  + f_1 * sgs0_3[k]
                  - f_2 * sgs1_3[k]
                  + f_3 * pc_y[k] * sgp_10[k];

        t_22[k] = f_7 * sfp_5[k]
                  + f_3 * pc_y[k] * sgp_11[k];

        t_23[k] = f_1 * sgs0_3[k]
                  - f_2 * sgs1_3[k]
                  + f_3 * pc_z[k] * sgp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, sfd0_12, sfp_13, sfp_14, sfd1_12, \
                         sgp_13, sgp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * sfd0_12[k]
                  - f_4 * pc_y[k] * sfd1_12[k];

        t_25[k] = f_7 * sfp_13[k]
                  + f_3 * pc_x[k] * sgp_13[k];

        t_26[k] = f_7 * sfp_14[k]
                  + f_3 * pc_x[k] * sgp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pc_y, pc_z, sfd0_9, sfd0_17, sfp_8, \
                         sfd1_9, sfd1_17, sgp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_z[k] * sfd0_9[k]
                  - f_4 * pc_z[k] * sfd1_9[k];

        t_28[k] = f_6 * sfp_8[k]
                  + f_3 * pc_y[k] * sgp_14[k];

        t_29[k] = pb_y[k] * sfd0_17[k]
                  - f_4 * pc_y[k] * sfd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, sfp_15, sfp_16, sfp_17, \
                         sgs0_5, sgs1_5, sgp_15, sgp_16, sgp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * sfp_15[k]
                  + f_1 * sgs0_5[k]
                  - f_2 * sgs1_5[k]
                  + f_3 * pc_x[k] * sgp_15[k];

        t_31[k] = f_7 * sfp_16[k]
                  + f_3 * pc_x[k] * sgp_16[k];

        t_32[k] = f_7 * sfp_17[k]
                  + f_3 * pc_x[k] * sgp_17[k];

        t_33[k] = f_1 * sgs0_5[k]
                  - f_2 * sgs1_5[k]
                  + f_3 * pc_y[k] * sgp_16[k];

        t_34[k] = f_3 * pc_y[k] * sgp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, pc_x, pc_z, sfd0_36, sfp_8, sfp_18, sfp_19, \
                         sfd1_36, sgs0_5, sgs1_5, sgp_17, sgp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_7 * sfp_8[k]
                  + f_1 * sgs0_5[k]
                  - f_2 * sgs1_5[k]
                  + f_3 * pc_z[k] * sgp_17[k];

        t_36[k] = pb_x[k] * sfd0_36[k]
                  + f_7 * sfp_18[k]
                  - f_4 * pc_x[k] * sfd1_36[k];

        t_37[k] = f_6 * sfp_19[k]
                  + f_3 * pc_x[k] * sgp_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pc_x, pc_y, sfd0_39, sfd0_41, sfp_11, \
                         sfp_20, sfd1_39, sfd1_41, sgp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * sfp_20[k]
                  + f_3 * pc_x[k] * sgp_20[k];

        t_39[k] = pb_x[k] * sfd0_39[k]
                  - f_4 * pc_x[k] * sfd1_39[k];

        t_40[k] = f_5 * sfp_11[k]
                  + f_3 * pc_y[k] * sgp_20[k];

        t_41[k] = pb_x[k] * sfd0_41[k]
                  - f_4 * pc_x[k] * sfd1_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_x, pb_z, pc_x, pc_z, sfd0_18, sfd0_45, \
                         sfp_22, sfp_23, sfd1_18, sfd1_45, sgp_22, \
                         sgp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sfd0_18[k]
                  - f_4 * pc_z[k] * sfd1_18[k];

        t_43[k] = f_6 * sfp_22[k]
                  + f_3 * pc_x[k] * sgp_22[k];

        t_44[k] = f_6 * sfp_23[k]
                  + f_3 * pc_x[k] * sgp_23[k];

        t_45[k] = pb_x[k] * sfd0_45[k]
                  - f_4 * pc_x[k] * sfd1_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, pb_y, pc_x, pc_y, sfd0_30, sfd0_47, \
                         sfp_14, sfp_25, sfd1_30, sfd1_47, sgp_23, \
                         sgp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_7 * sfp_14[k]
                  + f_3 * pc_y[k] * sgp_23[k];

        t_47[k] = pb_x[k] * sfd0_47[k]
                  - f_4 * pc_x[k] * sfd1_47[k];

        t_48[k] = pb_y[k] * sfd0_30[k]
                  - f_4 * pc_y[k] * sfd1_30[k];

        t_49[k] = f_6 * sfp_25[k]
                  + f_3 * pc_x[k] * sgp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, pc_x, pc_y, sfd0_51, sfd0_53, sfp_17, \
                         sfp_26, sfd1_51, sfd1_53, sgp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * sfp_26[k]
                  + f_3 * pc_x[k] * sgp_26[k];

        t_51[k] = pb_x[k] * sfd0_51[k]
                  - f_4 * pc_x[k] * sfd1_51[k];

        t_52[k] = f_6 * sfp_17[k]
                  + f_3 * pc_y[k] * sgp_26[k];

        t_53[k] = pb_x[k] * sfd0_53[k]
                  - f_4 * pc_x[k] * sfd1_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_x, pc_x, sfd0_54, sfd0_57, sfp_27, sfp_28, \
                         sfp_29, sfd1_54, sfd1_57, sgp_28, sgp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_x[k] * sfd0_54[k]
                  + f_7 * sfp_27[k]
                  - f_4 * pc_x[k] * sfd1_54[k];

        t_55[k] = f_6 * sfp_28[k]
                  + f_3 * pc_x[k] * sgp_28[k];

        t_56[k] = f_6 * sfp_29[k]
                  + f_3 * pc_x[k] * sgp_29[k];

        t_57[k] = pb_x[k] * sfd0_57[k]
                  - f_4 * pc_x[k] * sfd1_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pc_x, pc_y, sfd0_59, sfd1_59, \
                         sgs0_10, sgs1_10, sgp_29, sgp_30, sgp_31, \
                         sgp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * pc_y[k] * sgp_29[k];

        t_59[k] = pb_x[k] * sfd0_59[k]
                  - f_4 * pc_x[k] * sfd1_59[k];

        t_60[k] = f_1 * sgs0_10[k]
                  - f_2 * sgs1_10[k]
                  + f_3 * pc_x[k] * sgp_30[k];

        t_61[k] = f_3 * pc_x[k] * sgp_31[k];

        t_62[k] = f_3 * pc_x[k] * sgp_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pb_z, pc_y, pc_z, sfd0_36, sfp_19, sfp_20, \
                         sfd1_36, sgs0_10, sgs1_10, sgp_31, sgp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_0 * sfp_19[k]
                  + f_1 * sgs0_10[k]
                  - f_2 * sgs1_10[k]
                  + f_3 * pc_y[k] * sgp_31[k];

        t_64[k] = f_0 * sfp_20[k]
                  + f_3 * pc_y[k] * sgp_32[k];

        t_65[k] = f_1 * sgs0_10[k]
                  - f_2 * sgs1_10[k]
                  + f_3 * pc_z[k] * sgp_32[k];

        t_66[k] = pb_z[k] * sfd0_36[k]
                  - f_4 * pc_z[k] * sfd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_z, pc_x, pc_y, pc_z, sfd0_39, sfp_23, \
                         sfd1_39, sgp_34, sgp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_x[k] * sgp_34[k];

        t_68[k] = f_3 * pc_x[k] * sgp_35[k];

        t_69[k] = pb_z[k] * sfd0_39[k]
                  - f_4 * pc_z[k] * sfd1_39[k];

        t_70[k] = f_5 * sfp_23[k]
                  + f_3 * pc_y[k] * sgp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_z, sfp_20, sgs0_11, sgs0_12, \
                         sgs1_11, sgs1_12, sgp_35, sgp_36, sgp_37, \
                         sgp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * sfp_20[k]
                  + f_1 * sgs0_11[k]
                  - f_2 * sgs1_11[k]
                  + f_3 * pc_z[k] * sgp_35[k];

        t_72[k] = f_1 * sgs0_12[k]
                  - f_2 * sgs1_12[k]
                  + f_3 * pc_x[k] * sgp_36[k];

        t_73[k] = f_3 * pc_x[k] * sgp_37[k];

        t_74[k] = f_3 * pc_x[k] * sgp_38[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_y, pc_y, pc_z, sfd0_54, sfp_23, sfp_25, \
                         sfp_26, sfd1_54, sgs0_12, sgs1_12, sgp_37, \
                         sgp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_7 * sfp_25[k]
                  + f_1 * sgs0_12[k]
                  - f_2 * sgs1_12[k]
                  + f_3 * pc_y[k] * sgp_37[k];

        t_76[k] = f_7 * sfp_26[k]
                  + f_3 * pc_y[k] * sgp_38[k];

        t_77[k] = f_7 * sfp_23[k]
                  + f_1 * sgs0_12[k]
                  - f_2 * sgs1_12[k]
                  + f_3 * pc_z[k] * sgp_38[k];

        t_78[k] = pb_y[k] * sfd0_54[k]
                  - f_4 * pc_y[k] * sfd1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_y, pc_x, pc_y, sfd0_57, sfd0_59, \
                         sfp_28, sfp_29, sfd1_57, sfd1_59, sgp_40, \
                         sgp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_x[k] * sgp_40[k];

        t_80[k] = f_3 * pc_x[k] * sgp_41[k];

        t_81[k] = pb_y[k] * sfd0_57[k]
                  + f_7 * sfp_28[k]
                  - f_4 * pc_y[k] * sfd1_57[k];

        t_82[k] = f_6 * sfp_29[k]
                  + f_3 * pc_y[k] * sgp_41[k];

        t_83[k] = pb_y[k] * sfd0_59[k]
                  - f_4 * pc_y[k] * sfd1_59[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, sfp_29, \
                         sgs0_14, sgs1_14, sgp_42, sgp_43, sgp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * sgs0_14[k]
                  - f_2 * sgs1_14[k]
                  + f_3 * pc_x[k] * sgp_42[k];

        t_85[k] = f_3 * pc_x[k] * sgp_43[k];

        t_86[k] = f_3 * pc_x[k] * sgp_44[k];

        t_87[k] = f_1 * sgs0_14[k]
                  - f_2 * sgs1_14[k]
                  + f_3 * pc_y[k] * sgp_43[k];

        t_88[k] = f_3 * pc_y[k] * sgp_44[k];

        t_89[k] = f_0 * sfp_29[k]
                  + f_1 * sgs0_14[k]
                  - f_2 * sgs1_14[k]
                  + f_3 * pc_z[k] * sgp_44[k];
    }
}

}  // namespace simdt3ceri
