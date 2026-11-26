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


#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sfg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdg0,
                                                          const size_t sdf, const size_t sdg1,
                                                          const size_t sfd0, const size_t sfd1,
                                                          const size_t sff, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;

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
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdg0_0 = buffer.data(sdg0 + 0);
    const auto *sdg0_3 = buffer.data(sdg0 + 3);
    const auto *sdg0_5 = buffer.data(sdg0 + 5);
    const auto *sdg0_10 = buffer.data(sdg0 + 10);
    const auto *sdg0_14 = buffer.data(sdg0 + 14);
    const auto *sdg0_18 = buffer.data(sdg0 + 18);
    const auto *sdg0_30 = buffer.data(sdg0 + 30);
    const auto *sdg0_35 = buffer.data(sdg0 + 35);
    const auto *sdg0_45 = buffer.data(sdg0 + 45);
    const auto *sdg0_48 = buffer.data(sdg0 + 48);
    const auto *sdg0_50 = buffer.data(sdg0 + 50);
    const auto *sdg0_55 = buffer.data(sdg0 + 55);
    const auto *sdg0_57 = buffer.data(sdg0 + 57);
    const auto *sdg0_59 = buffer.data(sdg0 + 59);
    const auto *sdg0_70 = buffer.data(sdg0 + 70);
    const auto *sdg0_72 = buffer.data(sdg0 + 72);
    const auto *sdg0_74 = buffer.data(sdg0 + 74);
    const auto *sdg0_75 = buffer.data(sdg0 + 75);
    const auto *sdg0_78 = buffer.data(sdg0 + 78);
    const auto *sdg0_80 = buffer.data(sdg0 + 80);
    const auto *sdg0_85 = buffer.data(sdg0 + 85);
    const auto *sdg0_87 = buffer.data(sdg0 + 87);
    const auto *sdg0_89 = buffer.data(sdg0 + 89);

    const auto *sdf_0 = buffer.data(sdf + 0);
    const auto *sdf_1 = buffer.data(sdf + 1);
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

    const auto *sdg1_0 = buffer.data(sdg1 + 0);
    const auto *sdg1_3 = buffer.data(sdg1 + 3);
    const auto *sdg1_5 = buffer.data(sdg1 + 5);
    const auto *sdg1_10 = buffer.data(sdg1 + 10);
    const auto *sdg1_14 = buffer.data(sdg1 + 14);
    const auto *sdg1_18 = buffer.data(sdg1 + 18);
    const auto *sdg1_30 = buffer.data(sdg1 + 30);
    const auto *sdg1_35 = buffer.data(sdg1 + 35);
    const auto *sdg1_45 = buffer.data(sdg1 + 45);
    const auto *sdg1_48 = buffer.data(sdg1 + 48);
    const auto *sdg1_50 = buffer.data(sdg1 + 50);
    const auto *sdg1_55 = buffer.data(sdg1 + 55);
    const auto *sdg1_57 = buffer.data(sdg1 + 57);
    const auto *sdg1_59 = buffer.data(sdg1 + 59);
    const auto *sdg1_70 = buffer.data(sdg1 + 70);
    const auto *sdg1_72 = buffer.data(sdg1 + 72);
    const auto *sdg1_74 = buffer.data(sdg1 + 74);
    const auto *sdg1_75 = buffer.data(sdg1 + 75);
    const auto *sdg1_78 = buffer.data(sdg1 + 78);
    const auto *sdg1_80 = buffer.data(sdg1 + 80);
    const auto *sdg1_85 = buffer.data(sdg1 + 85);
    const auto *sdg1_87 = buffer.data(sdg1 + 87);
    const auto *sdg1_89 = buffer.data(sdg1 + 89);

    const auto *sfd0_0 = buffer.data(sfd0 + 0);
    const auto *sfd0_3 = buffer.data(sfd0 + 3);
    const auto *sfd0_5 = buffer.data(sfd0 + 5);
    const auto *sfd0_9 = buffer.data(sfd0 + 9);
    const auto *sfd0_11 = buffer.data(sfd0 + 11);
    const auto *sfd0_17 = buffer.data(sfd0 + 17);
    const auto *sfd0_36 = buffer.data(sfd0 + 36);
    const auto *sfd0_39 = buffer.data(sfd0 + 39);
    const auto *sfd0_41 = buffer.data(sfd0 + 41);
    const auto *sfd0_47 = buffer.data(sfd0 + 47);
    const auto *sfd0_51 = buffer.data(sfd0 + 51);

    const auto *sfd1_0 = buffer.data(sfd1 + 0);
    const auto *sfd1_3 = buffer.data(sfd1 + 3);
    const auto *sfd1_5 = buffer.data(sfd1 + 5);
    const auto *sfd1_9 = buffer.data(sfd1 + 9);
    const auto *sfd1_11 = buffer.data(sfd1 + 11);
    const auto *sfd1_17 = buffer.data(sfd1 + 17);
    const auto *sfd1_36 = buffer.data(sfd1 + 36);
    const auto *sfd1_39 = buffer.data(sfd1 + 39);
    const auto *sfd1_41 = buffer.data(sfd1 + 41);
    const auto *sfd1_47 = buffer.data(sfd1 + 47);
    const auto *sfd1_51 = buffer.data(sfd1 + 51);

    const auto *sff_0 = buffer.data(sff + 0);
    const auto *sff_2 = buffer.data(sff + 2);
    const auto *sff_3 = buffer.data(sff + 3);
    const auto *sff_5 = buffer.data(sff + 5);
    const auto *sff_6 = buffer.data(sff + 6);
    const auto *sff_7 = buffer.data(sff + 7);
    const auto *sff_8 = buffer.data(sff + 8);
    const auto *sff_9 = buffer.data(sff + 9);
    const auto *sff_10 = buffer.data(sff + 10);
    const auto *sff_12 = buffer.data(sff + 12);
    const auto *sff_16 = buffer.data(sff + 16);
    const auto *sff_17 = buffer.data(sff + 17);
    const auto *sff_18 = buffer.data(sff + 18);
    const auto *sff_19 = buffer.data(sff + 19);
    const auto *sff_20 = buffer.data(sff + 20);
    const auto *sff_22 = buffer.data(sff + 22);
    const auto *sff_26 = buffer.data(sff + 26);
    const auto *sff_27 = buffer.data(sff + 27);
    const auto *sff_28 = buffer.data(sff + 28);
    const auto *sff_29 = buffer.data(sff + 29);
    const auto *sff_30 = buffer.data(sff + 30);
    const auto *sff_32 = buffer.data(sff + 32);
    const auto *sff_36 = buffer.data(sff + 36);
    const auto *sff_37 = buffer.data(sff + 37);
    const auto *sff_38 = buffer.data(sff + 38);
    const auto *sff_39 = buffer.data(sff + 39);
    const auto *sff_40 = buffer.data(sff + 40);
    const auto *sff_42 = buffer.data(sff + 42);
    const auto *sff_46 = buffer.data(sff + 46);
    const auto *sff_47 = buffer.data(sff + 47);
    const auto *sff_48 = buffer.data(sff + 48);
    const auto *sff_49 = buffer.data(sff + 49);
    const auto *sff_50 = buffer.data(sff + 50);
    const auto *sff_52 = buffer.data(sff + 52);
    const auto *sff_56 = buffer.data(sff + 56);
    const auto *sff_57 = buffer.data(sff + 57);
    const auto *sff_58 = buffer.data(sff + 58);
    const auto *sff_59 = buffer.data(sff + 59);
    const auto *sff_60 = buffer.data(sff + 60);
    const auto *sff_62 = buffer.data(sff + 62);
    const auto *sff_63 = buffer.data(sff + 63);
    const auto *sff_65 = buffer.data(sff + 65);
    const auto *sff_66 = buffer.data(sff + 66);
    const auto *sff_67 = buffer.data(sff + 67);
    const auto *sff_68 = buffer.data(sff + 68);
    const auto *sff_69 = buffer.data(sff + 69);
    const auto *sff_70 = buffer.data(sff + 70);
    const auto *sff_72 = buffer.data(sff + 72);
    const auto *sff_75 = buffer.data(sff + 75);
    const auto *sff_76 = buffer.data(sff + 76);
    const auto *sff_77 = buffer.data(sff + 77);
    const auto *sff_78 = buffer.data(sff + 78);
    const auto *sff_79 = buffer.data(sff + 79);
    const auto *sff_80 = buffer.data(sff + 80);
    const auto *sff_82 = buffer.data(sff + 82);
    const auto *sff_83 = buffer.data(sff + 83);
    const auto *sff_86 = buffer.data(sff + 86);
    const auto *sff_87 = buffer.data(sff + 87);
    const auto *sff_88 = buffer.data(sff + 88);
    const auto *sff_89 = buffer.data(sff + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sdf_0, sdf_3, sfd0_0, sfd0_3, \
                         sfd1_0, sfd1_3, sff_0, sff_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdf_0[k]
                 + f_1 * sfd0_0[k]
                 - f_2 * sfd1_0[k]
                 + f_3 * pc_x[k] * sff_0[k];

        t_1[k] = f_3 * pc_y[k] * sff_0[k];

        t_2[k] = f_3 * pc_z[k] * sff_0[k];

        t_3[k] = f_0 * sdf_3[k]
                 + f_4 * sfd0_3[k]
                 - f_5 * sfd1_3[k]
                 + f_3 * pc_x[k] * sff_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, sdf_5, sdf_6, sdf_7, sfd0_5, sfd1_5, \
                         sff_2, sff_5, sff_6, sff_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sff_2[k];

        t_5[k] = f_0 * sdf_5[k]
                 + f_4 * sfd0_5[k]
                 - f_5 * sfd1_5[k]
                 + f_3 * pc_x[k] * sff_5[k];

        t_6[k] = f_0 * sdf_6[k]
                 + f_3 * pc_x[k] * sff_6[k];

        t_7[k] = f_0 * sdf_7[k]
                 + f_3 * pc_x[k] * sff_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, sdf_8, sdf_9, sfd0_3, sfd1_3, \
                         sff_6, sff_8, sff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sdf_8[k]
                 + f_3 * pc_x[k] * sff_8[k];

        t_9[k] = f_0 * sdf_9[k]
                 + f_3 * pc_x[k] * sff_9[k];

        t_10[k] = f_1 * sfd0_3[k]
                  - f_2 * sfd1_3[k]
                  + f_3 * pc_y[k] * sff_6[k];

        t_11[k] = f_3 * pc_z[k] * sff_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, sdg0_0, sdf_0, \
                         sdg1_0, sfd0_5, sfd1_5, sff_8, sff_9, sff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * sfd0_5[k]
                  - f_5 * sfd1_5[k]
                  + f_3 * pc_y[k] * sff_8[k];

        t_13[k] = f_3 * pc_y[k] * sff_9[k];

        t_14[k] = f_1 * sfd0_5[k]
                  - f_2 * sfd1_5[k]
                  + f_3 * pc_z[k] * sff_9[k];

        t_15[k] = pb_y[k] * sdg0_0[k]
                  - f_6 * pc_y[k] * sdg1_0[k];

        t_16[k] = f_7 * sdf_0[k]
                  + f_3 * pc_y[k] * sff_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, sdg0_3, sdg0_5, sdf_1, \
                         sdf_2, sdg1_3, sdg1_5, sff_10, sff_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * sff_10[k];

        t_18[k] = pb_y[k] * sdg0_3[k]
                  + f_8 * sdf_1[k]
                  - f_6 * pc_y[k] * sdg1_3[k];

        t_19[k] = f_7 * sdf_2[k]
                  + f_3 * pc_y[k] * sff_12[k];

        t_20[k] = pb_y[k] * sdg0_5[k]
                  - f_6 * pc_y[k] * sdg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, sdf_16, sdf_17, sdf_18, sdf_19, sff_16, \
                         sff_17, sff_18, sff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * sdf_16[k]
                  + f_3 * pc_x[k] * sff_16[k];

        t_22[k] = f_8 * sdf_17[k]
                  + f_3 * pc_x[k] * sff_17[k];

        t_23[k] = f_8 * sdf_18[k]
                  + f_3 * pc_x[k] * sff_18[k];

        t_24[k] = f_8 * sdf_19[k]
                  + f_3 * pc_x[k] * sff_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, sdf_6, sdf_8, sdf_9, sfd0_9, \
                         sfd0_11, sfd1_9, sfd1_11, sff_16, sff_18, \
                         sff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * sdf_6[k]
                  + f_1 * sfd0_9[k]
                  - f_2 * sfd1_9[k]
                  + f_3 * pc_y[k] * sff_16[k];

        t_26[k] = f_3 * pc_z[k] * sff_16[k];

        t_27[k] = f_7 * sdf_8[k]
                  + f_4 * sfd0_11[k]
                  - f_5 * sfd1_11[k]
                  + f_3 * pc_y[k] * sff_18[k];

        t_28[k] = f_7 * sdf_9[k]
                  + f_3 * pc_y[k] * sff_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, sdg0_0, sdg0_14, \
                         sdf_0, sdg1_0, sdg1_14, sff_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * sdg0_14[k]
                  - f_6 * pc_y[k] * sdg1_14[k];

        t_30[k] = pb_z[k] * sdg0_0[k]
                  - f_6 * pc_z[k] * sdg1_0[k];

        t_31[k] = f_3 * pc_y[k] * sff_20[k];

        t_32[k] = f_7 * sdf_0[k]
                  + f_3 * pc_z[k] * sff_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, sdg0_3, sdg0_5, \
                         sdf_2, sdf_26, sdg1_3, sdg1_5, sff_22, \
                         sff_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * sdg0_3[k]
                  - f_6 * pc_z[k] * sdg1_3[k];

        t_34[k] = f_3 * pc_y[k] * sff_22[k];

        t_35[k] = pb_z[k] * sdg0_5[k]
                  + f_8 * sdf_2[k]
                  - f_6 * pc_z[k] * sdg1_5[k];

        t_36[k] = f_8 * sdf_26[k]
                  + f_3 * pc_x[k] * sff_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, sdg0_10, sdf_27, sdf_28, \
                         sdf_29, sdg1_10, sff_27, sff_28, sff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_8 * sdf_27[k]
                  + f_3 * pc_x[k] * sff_27[k];

        t_38[k] = f_8 * sdf_28[k]
                  + f_3 * pc_x[k] * sff_28[k];

        t_39[k] = f_8 * sdf_29[k]
                  + f_3 * pc_x[k] * sff_29[k];

        t_40[k] = pb_z[k] * sdg0_10[k]
                  - f_6 * pc_z[k] * sdg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sdf_6, sdf_9, sfd0_17, sfd1_17, \
                         sff_26, sff_28, sff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * sdf_6[k]
                  + f_3 * pc_z[k] * sff_26[k];

        t_42[k] = f_4 * sfd0_17[k]
                  - f_5 * sfd1_17[k]
                  + f_3 * pc_y[k] * sff_28[k];

        t_43[k] = f_3 * pc_y[k] * sff_29[k];

        t_44[k] = f_7 * sdf_9[k]
                  + f_1 * sfd0_17[k]
                  - f_2 * sfd1_17[k]
                  + f_3 * pc_z[k] * sff_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_x, pc_x, pc_y, pc_z, sdg0_45, sdg0_48, \
                         sdf_10, sdf_30, sdf_33, sdg1_45, sdg1_48, \
                         sff_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_x[k] * sdg0_45[k]
                  + f_9 * sdf_30[k]
                  - f_6 * pc_x[k] * sdg1_45[k];

        t_46[k] = f_8 * sdf_10[k]
                  + f_3 * pc_y[k] * sff_30[k];

        t_47[k] = f_3 * pc_z[k] * sff_30[k];

        t_48[k] = pb_x[k] * sdg0_48[k]
                  + f_8 * sdf_33[k]
                  - f_6 * pc_x[k] * sdg1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_x, pc_x, pc_y, sdg0_50, sdf_12, sdf_35, \
                         sdf_36, sdf_37, sdg1_50, sff_32, sff_36, \
                         sff_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * sdf_12[k]
                  + f_3 * pc_y[k] * sff_32[k];

        t_50[k] = pb_x[k] * sdg0_50[k]
                  + f_8 * sdf_35[k]
                  - f_6 * pc_x[k] * sdg1_50[k];

        t_51[k] = f_7 * sdf_36[k]
                  + f_3 * pc_x[k] * sff_36[k];

        t_52[k] = f_7 * sdf_37[k]
                  + f_3 * pc_x[k] * sff_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, pc_x, pc_z, sdg0_55, sdf_38, sdf_39, \
                         sdg1_55, sff_36, sff_38, sff_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * sdf_38[k]
                  + f_3 * pc_x[k] * sff_38[k];

        t_54[k] = f_7 * sdf_39[k]
                  + f_3 * pc_x[k] * sff_39[k];

        t_55[k] = pb_x[k] * sdg0_55[k]
                  - f_6 * pc_x[k] * sdg1_55[k];

        t_56[k] = f_3 * pc_z[k] * sff_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_x, pb_y, pc_x, pc_y, sdg0_30, sdg0_57, \
                         sdg0_59, sdf_19, sdg1_30, sdg1_57, sdg1_59, \
                         sff_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_x[k] * sdg0_57[k]
                  - f_6 * pc_x[k] * sdg1_57[k];

        t_58[k] = f_8 * sdf_19[k]
                  + f_3 * pc_y[k] * sff_39[k];

        t_59[k] = pb_x[k] * sdg0_59[k]
                  - f_6 * pc_x[k] * sdg1_59[k];

        t_60[k] = pb_y[k] * sdg0_30[k]
                  - f_6 * pc_y[k] * sdg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sdg0_18, sdf_10, sdf_20, \
                         sdf_22, sdg1_18, sff_40, sff_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * sdf_20[k]
                  + f_3 * pc_y[k] * sff_40[k];

        t_62[k] = f_7 * sdf_10[k]
                  + f_3 * pc_z[k] * sff_40[k];

        t_63[k] = pb_z[k] * sdg0_18[k]
                  - f_6 * pc_z[k] * sdg1_18[k];

        t_64[k] = f_7 * sdf_22[k]
                  + f_3 * pc_y[k] * sff_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, sdg0_35, sdf_46, sdf_47, \
                         sdf_48, sdg1_35, sff_46, sff_47, sff_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sdg0_35[k]
                  - f_6 * pc_y[k] * sdg1_35[k];

        t_66[k] = f_7 * sdf_46[k]
                  + f_3 * pc_x[k] * sff_46[k];

        t_67[k] = f_7 * sdf_47[k]
                  + f_3 * pc_x[k] * sff_47[k];

        t_68[k] = f_7 * sdf_48[k]
                  + f_3 * pc_x[k] * sff_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pc_x, pc_z, sdg0_70, sdg0_72, sdf_16, \
                         sdf_49, sdg1_70, sdg1_72, sff_46, sff_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_7 * sdf_49[k]
                  + f_3 * pc_x[k] * sff_49[k];

        t_70[k] = pb_x[k] * sdg0_70[k]
                  - f_6 * pc_x[k] * sdg1_70[k];

        t_71[k] = f_7 * sdf_16[k]
                  + f_3 * pc_z[k] * sff_46[k];

        t_72[k] = pb_x[k] * sdg0_72[k]
                  - f_6 * pc_x[k] * sdg1_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pb_x, pc_x, pc_y, sdg0_74, sdg0_75, sdf_29, \
                         sdf_50, sdg1_74, sdg1_75, sff_49, sff_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * sdf_29[k]
                  + f_3 * pc_y[k] * sff_49[k];

        t_74[k] = pb_x[k] * sdg0_74[k]
                  - f_6 * pc_x[k] * sdg1_74[k];

        t_75[k] = pb_x[k] * sdg0_75[k]
                  + f_9 * sdf_50[k]
                  - f_6 * pc_x[k] * sdg1_75[k];

        t_76[k] = f_3 * pc_y[k] * sff_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_x, pc_x, pc_y, pc_z, sdg0_78, sdf_20, sdf_53, \
                         sdg1_78, sff_50, sff_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * sdf_20[k]
                  + f_3 * pc_z[k] * sff_50[k];

        t_78[k] = pb_x[k] * sdg0_78[k]
                  + f_8 * sdf_53[k]
                  - f_6 * pc_x[k] * sdg1_78[k];

        t_79[k] = f_3 * pc_y[k] * sff_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pb_x, pc_x, sdg0_80, sdf_55, sdf_56, sdf_57, \
                         sdf_58, sdg1_80, sff_56, sff_57, sff_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * sdg0_80[k]
                  + f_8 * sdf_55[k]
                  - f_6 * pc_x[k] * sdg1_80[k];

        t_81[k] = f_7 * sdf_56[k]
                  + f_3 * pc_x[k] * sff_56[k];

        t_82[k] = f_7 * sdf_57[k]
                  + f_3 * pc_x[k] * sff_57[k];

        t_83[k] = f_7 * sdf_58[k]
                  + f_3 * pc_x[k] * sff_58[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pc_x, pc_z, sdg0_85, sdg0_87, sdf_26, \
                         sdf_59, sdg1_85, sdg1_87, sff_56, sff_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_7 * sdf_59[k]
                  + f_3 * pc_x[k] * sff_59[k];

        t_85[k] = pb_x[k] * sdg0_85[k]
                  - f_6 * pc_x[k] * sdg1_85[k];

        t_86[k] = f_8 * sdf_26[k]
                  + f_3 * pc_z[k] * sff_56[k];

        t_87[k] = pb_x[k] * sdg0_87[k]
                  - f_6 * pc_x[k] * sdg1_87[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pb_x, pc_x, pc_y, pc_z, sdg0_89, \
                         sdf_30, sdg1_89, sfd0_36, sfd1_36, sff_59, \
                         sff_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * pc_y[k] * sff_59[k];

        t_89[k] = pb_x[k] * sdg0_89[k]
                  - f_6 * pc_x[k] * sdg1_89[k];

        t_90[k] = f_1 * sfd0_36[k]
                  - f_2 * sfd1_36[k]
                  + f_3 * pc_x[k] * sff_60[k];

        t_91[k] = f_0 * sdf_30[k]
                  + f_3 * pc_y[k] * sff_60[k];

        t_92[k] = f_3 * pc_z[k] * sff_60[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, sdf_32, sfd0_39, sfd0_41, \
                         sfd1_39, sfd1_41, sff_62, sff_63, sff_65, \
                         sff_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_4 * sfd0_39[k]
                  - f_5 * sfd1_39[k]
                  + f_3 * pc_x[k] * sff_63[k];

        t_94[k] = f_0 * sdf_32[k]
                  + f_3 * pc_y[k] * sff_62[k];

        t_95[k] = f_4 * sfd0_41[k]
                  - f_5 * sfd1_41[k]
                  + f_3 * pc_x[k] * sff_65[k];

        t_96[k] = f_3 * pc_x[k] * sff_66[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pc_x, pc_y, pc_z, sdf_36, sfd0_39, \
                         sfd1_39, sff_66, sff_67, sff_68, sff_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * pc_x[k] * sff_67[k];

        t_98[k] = f_3 * pc_x[k] * sff_68[k];

        t_99[k] = f_3 * pc_x[k] * sff_69[k];

        t_100[k] = f_0 * sdf_36[k]
                   + f_1 * sfd0_39[k]
                   - f_2 * sfd1_39[k]
                   + f_3 * pc_y[k] * sff_66[k];

        t_101[k] = f_3 * pc_z[k] * sff_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, sdg0_45, sdf_38, \
                         sdf_39, sdg1_45, sfd0_41, sfd1_41, sff_68, \
                         sff_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_0 * sdf_38[k]
                   + f_4 * sfd0_41[k]
                   - f_5 * sfd1_41[k]
                   + f_3 * pc_y[k] * sff_68[k];

        t_103[k] = f_0 * sdf_39[k]
                   + f_3 * pc_y[k] * sff_69[k];

        t_104[k] = f_1 * sfd0_41[k]
                   - f_2 * sfd1_41[k]
                   + f_3 * pc_z[k] * sff_69[k];

        t_105[k] = pb_z[k] * sdg0_45[k]
                   - f_6 * pc_z[k] * sdg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, sdg0_48, sdf_30, \
                         sdf_40, sdf_42, sdg1_48, sff_70, sff_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * sdf_40[k]
                   + f_3 * pc_y[k] * sff_70[k];

        t_107[k] = f_7 * sdf_30[k]
                   + f_3 * pc_z[k] * sff_70[k];

        t_108[k] = pb_z[k] * sdg0_48[k]
                   - f_6 * pc_z[k] * sdg1_48[k];

        t_109[k] = f_8 * sdf_42[k]
                   + f_3 * pc_y[k] * sff_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pc_x, sfd0_47, sfd1_47, sff_75, \
                         sff_76, sff_77, sff_78, sff_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_4 * sfd0_47[k]
                   - f_5 * sfd1_47[k]
                   + f_3 * pc_x[k] * sff_75[k];

        t_111[k] = f_3 * pc_x[k] * sff_76[k];

        t_112[k] = f_3 * pc_x[k] * sff_77[k];

        t_113[k] = f_3 * pc_x[k] * sff_78[k];

        t_114[k] = f_3 * pc_x[k] * sff_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_y, pc_z, sdg0_55, sdg0_57, \
                         sdf_36, sdf_37, sdf_49, sdg1_55, sdg1_57, sff_76, \
                         sff_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * sdg0_55[k]
                   - f_6 * pc_z[k] * sdg1_55[k];

        t_116[k] = f_7 * sdf_36[k]
                   + f_3 * pc_z[k] * sff_76[k];

        t_117[k] = pb_z[k] * sdg0_57[k]
                   + f_8 * sdf_37[k]
                   - f_6 * pc_z[k] * sdg1_57[k];

        t_118[k] = f_8 * sdf_49[k]
                   + f_3 * pc_y[k] * sff_79[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_y, pc_y, pc_z, sdg0_75, sdf_39, \
                         sdf_40, sdf_50, sdg1_75, sfd0_47, sfd1_47, sff_79, \
                         sff_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_7 * sdf_39[k]
                   + f_1 * sfd0_47[k]
                   - f_2 * sfd1_47[k]
                   + f_3 * pc_z[k] * sff_79[k];

        t_120[k] = pb_y[k] * sdg0_75[k]
                   - f_6 * pc_y[k] * sdg1_75[k];

        t_121[k] = f_7 * sdf_50[k]
                   + f_3 * pc_y[k] * sff_80[k];

        t_122[k] = f_8 * sdf_40[k]
                   + f_3 * pc_z[k] * sff_80[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_y, pc_x, pc_y, sdg0_80, sdf_52, \
                         sdg1_80, sfd0_51, sfd1_51, sff_82, sff_83, \
                         sff_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_4 * sfd0_51[k]
                   - f_5 * sfd1_51[k]
                   + f_3 * pc_x[k] * sff_83[k];

        t_124[k] = f_7 * sdf_52[k]
                   + f_3 * pc_y[k] * sff_82[k];

        t_125[k] = pb_y[k] * sdg0_80[k]
                   - f_6 * pc_y[k] * sdg1_80[k];

        t_126[k] = f_3 * pc_x[k] * sff_86[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pb_y, pc_x, pc_y, sdg0_85, sdf_56, \
                         sdg1_85, sff_87, sff_88, sff_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_3 * pc_x[k] * sff_87[k];

        t_128[k] = f_3 * pc_x[k] * sff_88[k];

        t_129[k] = f_3 * pc_x[k] * sff_89[k];

        t_130[k] = pb_y[k] * sdg0_85[k]
                   + f_9 * sdf_56[k]
                   - f_6 * pc_y[k] * sdg1_85[k];
    }
}

static auto
compute_prim_sfg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdg0,
                                                          const size_t sdf, const size_t sdg1,
                                                          const size_t sfd0, const size_t sfd1,
                                                          const size_t sff, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;

    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdg0_87 = buffer.data(sdg0 + 87);
    const auto *sdg0_89 = buffer.data(sdg0 + 89);

    const auto *sdf_46 = buffer.data(sdf + 46);
    const auto *sdf_50 = buffer.data(sdf + 50);
    const auto *sdf_56 = buffer.data(sdf + 56);
    const auto *sdf_58 = buffer.data(sdf + 58);
    const auto *sdf_59 = buffer.data(sdf + 59);

    const auto *sdg1_87 = buffer.data(sdg1 + 87);
    const auto *sdg1_89 = buffer.data(sdg1 + 89);

    const auto *sfd0_54 = buffer.data(sfd0 + 54);
    const auto *sfd0_57 = buffer.data(sfd0 + 57);
    const auto *sfd0_59 = buffer.data(sfd0 + 59);

    const auto *sfd1_54 = buffer.data(sfd1 + 54);
    const auto *sfd1_57 = buffer.data(sfd1 + 57);
    const auto *sfd1_59 = buffer.data(sfd1 + 59);

    const auto *sff_86 = buffer.data(sff + 86);
    const auto *sff_89 = buffer.data(sff + 89);
    const auto *sff_90 = buffer.data(sff + 90);
    const auto *sff_92 = buffer.data(sff + 92);
    const auto *sff_93 = buffer.data(sff + 93);
    const auto *sff_95 = buffer.data(sff + 95);
    const auto *sff_96 = buffer.data(sff + 96);
    const auto *sff_97 = buffer.data(sff + 97);
    const auto *sff_98 = buffer.data(sff + 98);
    const auto *sff_99 = buffer.data(sff + 99);

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pb_y, pc_y, pc_z, sdg0_87, sdg0_89, \
                         sdf_46, sdf_58, sdf_59, sdg1_87, sdg1_89, sff_86, \
                         sff_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_8 * sdf_46[k]
                   + f_3 * pc_z[k] * sff_86[k];

        t_132[k] = pb_y[k] * sdg0_87[k]
                   + f_8 * sdf_58[k]
                   - f_6 * pc_y[k] * sdg1_87[k];

        t_133[k] = f_7 * sdf_59[k]
                   + f_3 * pc_y[k] * sff_89[k];

        t_134[k] = pb_y[k] * sdg0_89[k]
                   - f_6 * pc_y[k] * sdg1_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, sdf_50, sfd0_54, \
                         sfd0_57, sfd1_54, sfd1_57, sff_90, sff_92, \
                         sff_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_1 * sfd0_54[k]
                   - f_2 * sfd1_54[k]
                   + f_3 * pc_x[k] * sff_90[k];

        t_136[k] = f_3 * pc_y[k] * sff_90[k];

        t_137[k] = f_0 * sdf_50[k]
                   + f_3 * pc_z[k] * sff_90[k];

        t_138[k] = f_4 * sfd0_57[k]
                   - f_5 * sfd1_57[k]
                   + f_3 * pc_x[k] * sff_93[k];

        t_139[k] = f_3 * pc_y[k] * sff_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, sfd0_59, sfd1_59, sff_95, \
                         sff_96, sff_97, sff_98, sff_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * sfd0_59[k]
                   - f_5 * sfd1_59[k]
                   + f_3 * pc_x[k] * sff_95[k];

        t_141[k] = f_3 * pc_x[k] * sff_96[k];

        t_142[k] = f_3 * pc_x[k] * sff_97[k];

        t_143[k] = f_3 * pc_x[k] * sff_98[k];

        t_144[k] = f_3 * pc_x[k] * sff_99[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pc_y, pc_z, sdf_56, sdf_59, \
                         sfd0_57, sfd0_59, sfd1_57, sfd1_59, sff_96, sff_98, \
                         sff_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_1 * sfd0_57[k]
                   - f_2 * sfd1_57[k]
                   + f_3 * pc_y[k] * sff_96[k];

        t_146[k] = f_0 * sdf_56[k]
                   + f_3 * pc_z[k] * sff_96[k];

        t_147[k] = f_4 * sfd0_59[k]
                   - f_5 * sfd1_59[k]
                   + f_3 * pc_y[k] * sff_98[k];

        t_148[k] = f_3 * pc_y[k] * sff_99[k];

        t_149[k] = f_0 * sdf_59[k]
                   + f_1 * sfd0_59[k]
                   - f_2 * sfd1_59[k]
                   + f_3 * pc_z[k] * sff_99[k];
    }
}

auto
compute_prim_sfg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sdg0, const size_t sdf,
                                                   const size_t sdg1, const size_t sfd0,
                                                   const size_t sfd1, const size_t sff,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sfg_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sdg0, sdf,
                                                              sdg1, sfd0, sfd1, sff, ncols,
                                                              gamma, p, q);

    compute_prim_sfg_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sdg0, sdf,
                                                              sdg1, sfd0, sfd1, sff, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
