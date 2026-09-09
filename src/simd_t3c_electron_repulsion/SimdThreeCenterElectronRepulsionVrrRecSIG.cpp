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


#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sig_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shg0,
                                                          const size_t shf, const size_t shg1,
                                                          const size_t sid0, const size_t sid1,
                                                          const size_t sif, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 2.0 / q;
    const auto f_11 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shg0_0 = buffer.data(shg0 + 0);
    const auto *shg0_3 = buffer.data(shg0 + 3);
    const auto *shg0_5 = buffer.data(shg0 + 5);
    const auto *shg0_10 = buffer.data(shg0 + 10);
    const auto *shg0_14 = buffer.data(shg0 + 14);
    const auto *shg0_18 = buffer.data(shg0 + 18);
    const auto *shg0_25 = buffer.data(shg0 + 25);
    const auto *shg0_30 = buffer.data(shg0 + 30);
    const auto *shg0_35 = buffer.data(shg0 + 35);
    const auto *shg0_44 = buffer.data(shg0 + 44);
    const auto *shg0_45 = buffer.data(shg0 + 45);
    const auto *shg0_48 = buffer.data(shg0 + 48);
    const auto *shg0_55 = buffer.data(shg0 + 55);
    const auto *shg0_75 = buffer.data(shg0 + 75);
    const auto *shg0_78 = buffer.data(shg0 + 78);

    const auto *shf_0 = buffer.data(shf + 0);
    const auto *shf_1 = buffer.data(shf + 1);
    const auto *shf_2 = buffer.data(shf + 2);
    const auto *shf_3 = buffer.data(shf + 3);
    const auto *shf_5 = buffer.data(shf + 5);
    const auto *shf_6 = buffer.data(shf + 6);
    const auto *shf_7 = buffer.data(shf + 7);
    const auto *shf_8 = buffer.data(shf + 8);
    const auto *shf_9 = buffer.data(shf + 9);
    const auto *shf_10 = buffer.data(shf + 10);
    const auto *shf_12 = buffer.data(shf + 12);
    const auto *shf_16 = buffer.data(shf + 16);
    const auto *shf_17 = buffer.data(shf + 17);
    const auto *shf_18 = buffer.data(shf + 18);
    const auto *shf_19 = buffer.data(shf + 19);
    const auto *shf_20 = buffer.data(shf + 20);
    const auto *shf_22 = buffer.data(shf + 22);
    const auto *shf_26 = buffer.data(shf + 26);
    const auto *shf_27 = buffer.data(shf + 27);
    const auto *shf_28 = buffer.data(shf + 28);
    const auto *shf_29 = buffer.data(shf + 29);
    const auto *shf_30 = buffer.data(shf + 30);
    const auto *shf_32 = buffer.data(shf + 32);
    const auto *shf_33 = buffer.data(shf + 33);
    const auto *shf_35 = buffer.data(shf + 35);
    const auto *shf_36 = buffer.data(shf + 36);
    const auto *shf_37 = buffer.data(shf + 37);
    const auto *shf_38 = buffer.data(shf + 38);
    const auto *shf_39 = buffer.data(shf + 39);
    const auto *shf_40 = buffer.data(shf + 40);
    const auto *shf_42 = buffer.data(shf + 42);
    const auto *shf_46 = buffer.data(shf + 46);
    const auto *shf_47 = buffer.data(shf + 47);
    const auto *shf_48 = buffer.data(shf + 48);
    const auto *shf_49 = buffer.data(shf + 49);
    const auto *shf_50 = buffer.data(shf + 50);
    const auto *shf_51 = buffer.data(shf + 51);
    const auto *shf_52 = buffer.data(shf + 52);
    const auto *shf_53 = buffer.data(shf + 53);
    const auto *shf_55 = buffer.data(shf + 55);
    const auto *shf_56 = buffer.data(shf + 56);
    const auto *shf_57 = buffer.data(shf + 57);
    const auto *shf_58 = buffer.data(shf + 58);
    const auto *shf_59 = buffer.data(shf + 59);
    const auto *shf_60 = buffer.data(shf + 60);
    const auto *shf_63 = buffer.data(shf + 63);
    const auto *shf_65 = buffer.data(shf + 65);
    const auto *shf_66 = buffer.data(shf + 66);
    const auto *shf_67 = buffer.data(shf + 67);
    const auto *shf_68 = buffer.data(shf + 68);
    const auto *shf_69 = buffer.data(shf + 69);
    const auto *shf_75 = buffer.data(shf + 75);
    const auto *shf_76 = buffer.data(shf + 76);
    const auto *shf_77 = buffer.data(shf + 77);
    const auto *shf_78 = buffer.data(shf + 78);
    const auto *shf_79 = buffer.data(shf + 79);

    const auto *shg1_0 = buffer.data(shg1 + 0);
    const auto *shg1_3 = buffer.data(shg1 + 3);
    const auto *shg1_5 = buffer.data(shg1 + 5);
    const auto *shg1_10 = buffer.data(shg1 + 10);
    const auto *shg1_14 = buffer.data(shg1 + 14);
    const auto *shg1_18 = buffer.data(shg1 + 18);
    const auto *shg1_25 = buffer.data(shg1 + 25);
    const auto *shg1_30 = buffer.data(shg1 + 30);
    const auto *shg1_35 = buffer.data(shg1 + 35);
    const auto *shg1_44 = buffer.data(shg1 + 44);
    const auto *shg1_45 = buffer.data(shg1 + 45);
    const auto *shg1_48 = buffer.data(shg1 + 48);
    const auto *shg1_55 = buffer.data(shg1 + 55);
    const auto *shg1_75 = buffer.data(shg1 + 75);
    const auto *shg1_78 = buffer.data(shg1 + 78);

    const auto *sid0_0 = buffer.data(sid0 + 0);
    const auto *sid0_3 = buffer.data(sid0 + 3);
    const auto *sid0_5 = buffer.data(sid0 + 5);
    const auto *sid0_9 = buffer.data(sid0 + 9);
    const auto *sid0_11 = buffer.data(sid0 + 11);
    const auto *sid0_17 = buffer.data(sid0 + 17);
    const auto *sid0_18 = buffer.data(sid0 + 18);
    const auto *sid0_21 = buffer.data(sid0 + 21);
    const auto *sid0_23 = buffer.data(sid0 + 23);
    const auto *sid0_29 = buffer.data(sid0 + 29);
    const auto *sid0_30 = buffer.data(sid0 + 30);
    const auto *sid0_33 = buffer.data(sid0 + 33);
    const auto *sid0_35 = buffer.data(sid0 + 35);
    const auto *sid0_36 = buffer.data(sid0 + 36);
    const auto *sid0_39 = buffer.data(sid0 + 39);
    const auto *sid0_41 = buffer.data(sid0 + 41);
    const auto *sid0_47 = buffer.data(sid0 + 47);

    const auto *sid1_0 = buffer.data(sid1 + 0);
    const auto *sid1_3 = buffer.data(sid1 + 3);
    const auto *sid1_5 = buffer.data(sid1 + 5);
    const auto *sid1_9 = buffer.data(sid1 + 9);
    const auto *sid1_11 = buffer.data(sid1 + 11);
    const auto *sid1_17 = buffer.data(sid1 + 17);
    const auto *sid1_18 = buffer.data(sid1 + 18);
    const auto *sid1_21 = buffer.data(sid1 + 21);
    const auto *sid1_23 = buffer.data(sid1 + 23);
    const auto *sid1_29 = buffer.data(sid1 + 29);
    const auto *sid1_30 = buffer.data(sid1 + 30);
    const auto *sid1_33 = buffer.data(sid1 + 33);
    const auto *sid1_35 = buffer.data(sid1 + 35);
    const auto *sid1_36 = buffer.data(sid1 + 36);
    const auto *sid1_39 = buffer.data(sid1 + 39);
    const auto *sid1_41 = buffer.data(sid1 + 41);
    const auto *sid1_47 = buffer.data(sid1 + 47);

    const auto *sif_0 = buffer.data(sif + 0);
    const auto *sif_2 = buffer.data(sif + 2);
    const auto *sif_3 = buffer.data(sif + 3);
    const auto *sif_5 = buffer.data(sif + 5);
    const auto *sif_6 = buffer.data(sif + 6);
    const auto *sif_7 = buffer.data(sif + 7);
    const auto *sif_8 = buffer.data(sif + 8);
    const auto *sif_9 = buffer.data(sif + 9);
    const auto *sif_10 = buffer.data(sif + 10);
    const auto *sif_12 = buffer.data(sif + 12);
    const auto *sif_16 = buffer.data(sif + 16);
    const auto *sif_17 = buffer.data(sif + 17);
    const auto *sif_18 = buffer.data(sif + 18);
    const auto *sif_19 = buffer.data(sif + 19);
    const auto *sif_20 = buffer.data(sif + 20);
    const auto *sif_22 = buffer.data(sif + 22);
    const auto *sif_26 = buffer.data(sif + 26);
    const auto *sif_27 = buffer.data(sif + 27);
    const auto *sif_28 = buffer.data(sif + 28);
    const auto *sif_29 = buffer.data(sif + 29);
    const auto *sif_30 = buffer.data(sif + 30);
    const auto *sif_32 = buffer.data(sif + 32);
    const auto *sif_33 = buffer.data(sif + 33);
    const auto *sif_35 = buffer.data(sif + 35);
    const auto *sif_36 = buffer.data(sif + 36);
    const auto *sif_37 = buffer.data(sif + 37);
    const auto *sif_38 = buffer.data(sif + 38);
    const auto *sif_39 = buffer.data(sif + 39);
    const auto *sif_40 = buffer.data(sif + 40);
    const auto *sif_42 = buffer.data(sif + 42);
    const auto *sif_46 = buffer.data(sif + 46);
    const auto *sif_47 = buffer.data(sif + 47);
    const auto *sif_48 = buffer.data(sif + 48);
    const auto *sif_49 = buffer.data(sif + 49);
    const auto *sif_50 = buffer.data(sif + 50);
    const auto *sif_52 = buffer.data(sif + 52);
    const auto *sif_53 = buffer.data(sif + 53);
    const auto *sif_55 = buffer.data(sif + 55);
    const auto *sif_56 = buffer.data(sif + 56);
    const auto *sif_57 = buffer.data(sif + 57);
    const auto *sif_58 = buffer.data(sif + 58);
    const auto *sif_59 = buffer.data(sif + 59);
    const auto *sif_60 = buffer.data(sif + 60);
    const auto *sif_62 = buffer.data(sif + 62);
    const auto *sif_63 = buffer.data(sif + 63);
    const auto *sif_65 = buffer.data(sif + 65);
    const auto *sif_66 = buffer.data(sif + 66);
    const auto *sif_67 = buffer.data(sif + 67);
    const auto *sif_68 = buffer.data(sif + 68);
    const auto *sif_69 = buffer.data(sif + 69);
    const auto *sif_70 = buffer.data(sif + 70);
    const auto *sif_72 = buffer.data(sif + 72);
    const auto *sif_75 = buffer.data(sif + 75);
    const auto *sif_76 = buffer.data(sif + 76);
    const auto *sif_77 = buffer.data(sif + 77);
    const auto *sif_78 = buffer.data(sif + 78);
    const auto *sif_79 = buffer.data(sif + 79);
    const auto *sif_80 = buffer.data(sif + 80);
    const auto *sif_82 = buffer.data(sif + 82);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, shf_0, shf_3, sid0_0, sid0_3, \
                         sid1_0, sid1_3, sif_0, sif_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * shf_0[k]
                 + f_1 * sid0_0[k]
                 - f_2 * sid1_0[k]
                 + f_3 * pc_x[k] * sif_0[k];

        t_1[k] = f_3 * pc_y[k] * sif_0[k];

        t_2[k] = f_3 * pc_z[k] * sif_0[k];

        t_3[k] = f_0 * shf_3[k]
                 + f_4 * sid0_3[k]
                 - f_5 * sid1_3[k]
                 + f_3 * pc_x[k] * sif_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, shf_5, shf_6, shf_7, sid0_5, sid1_5, \
                         sif_2, sif_5, sif_6, sif_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sif_2[k];

        t_5[k] = f_0 * shf_5[k]
                 + f_4 * sid0_5[k]
                 - f_5 * sid1_5[k]
                 + f_3 * pc_x[k] * sif_5[k];

        t_6[k] = f_0 * shf_6[k]
                 + f_3 * pc_x[k] * sif_6[k];

        t_7[k] = f_0 * shf_7[k]
                 + f_3 * pc_x[k] * sif_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, shf_8, shf_9, sid0_3, sid1_3, \
                         sif_6, sif_8, sif_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * shf_8[k]
                 + f_3 * pc_x[k] * sif_8[k];

        t_9[k] = f_0 * shf_9[k]
                 + f_3 * pc_x[k] * sif_9[k];

        t_10[k] = f_1 * sid0_3[k]
                  - f_2 * sid1_3[k]
                  + f_3 * pc_y[k] * sif_6[k];

        t_11[k] = f_3 * pc_z[k] * sif_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, shg0_0, shf_0, \
                         shg1_0, sid0_5, sid1_5, sif_8, sif_9, sif_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * sid0_5[k]
                  - f_5 * sid1_5[k]
                  + f_3 * pc_y[k] * sif_8[k];

        t_13[k] = f_3 * pc_y[k] * sif_9[k];

        t_14[k] = f_1 * sid0_5[k]
                  - f_2 * sid1_5[k]
                  + f_3 * pc_z[k] * sif_9[k];

        t_15[k] = pb_y[k] * shg0_0[k]
                  - f_6 * pc_y[k] * shg1_0[k];

        t_16[k] = f_7 * shf_0[k]
                  + f_3 * pc_y[k] * sif_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, shg0_3, shg0_5, shf_1, \
                         shf_2, shg1_3, shg1_5, sif_10, sif_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * sif_10[k];

        t_18[k] = pb_y[k] * shg0_3[k]
                  + f_8 * shf_1[k]
                  - f_6 * pc_y[k] * shg1_3[k];

        t_19[k] = f_7 * shf_2[k]
                  + f_3 * pc_y[k] * sif_12[k];

        t_20[k] = pb_y[k] * shg0_5[k]
                  - f_6 * pc_y[k] * shg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, shf_16, shf_17, shf_18, shf_19, sif_16, \
                         sif_17, sif_18, sif_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * shf_16[k]
                  + f_3 * pc_x[k] * sif_16[k];

        t_22[k] = f_9 * shf_17[k]
                  + f_3 * pc_x[k] * sif_17[k];

        t_23[k] = f_9 * shf_18[k]
                  + f_3 * pc_x[k] * sif_18[k];

        t_24[k] = f_9 * shf_19[k]
                  + f_3 * pc_x[k] * sif_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, shf_6, shf_8, shf_9, sid0_9, \
                         sid0_11, sid1_9, sid1_11, sif_16, sif_18, \
                         sif_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * shf_6[k]
                  + f_1 * sid0_9[k]
                  - f_2 * sid1_9[k]
                  + f_3 * pc_y[k] * sif_16[k];

        t_26[k] = f_3 * pc_z[k] * sif_16[k];

        t_27[k] = f_7 * shf_8[k]
                  + f_4 * sid0_11[k]
                  - f_5 * sid1_11[k]
                  + f_3 * pc_y[k] * sif_18[k];

        t_28[k] = f_7 * shf_9[k]
                  + f_3 * pc_y[k] * sif_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, shg0_0, shg0_14, \
                         shf_0, shg1_0, shg1_14, sif_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * shg0_14[k]
                  - f_6 * pc_y[k] * shg1_14[k];

        t_30[k] = pb_z[k] * shg0_0[k]
                  - f_6 * pc_z[k] * shg1_0[k];

        t_31[k] = f_3 * pc_y[k] * sif_20[k];

        t_32[k] = f_7 * shf_0[k]
                  + f_3 * pc_z[k] * sif_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, shg0_3, shg0_5, \
                         shf_2, shf_26, shg1_3, shg1_5, sif_22, \
                         sif_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * shg0_3[k]
                  - f_6 * pc_z[k] * shg1_3[k];

        t_34[k] = f_3 * pc_y[k] * sif_22[k];

        t_35[k] = pb_z[k] * shg0_5[k]
                  + f_8 * shf_2[k]
                  - f_6 * pc_z[k] * shg1_5[k];

        t_36[k] = f_9 * shf_26[k]
                  + f_3 * pc_x[k] * sif_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, shg0_10, shf_27, shf_28, \
                         shf_29, shg1_10, sif_27, sif_28, sif_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * shf_27[k]
                  + f_3 * pc_x[k] * sif_27[k];

        t_38[k] = f_9 * shf_28[k]
                  + f_3 * pc_x[k] * sif_28[k];

        t_39[k] = f_9 * shf_29[k]
                  + f_3 * pc_x[k] * sif_29[k];

        t_40[k] = pb_z[k] * shg0_10[k]
                  - f_6 * pc_z[k] * shg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, shf_6, shf_9, sid0_17, sid1_17, \
                         sif_26, sif_28, sif_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * shf_6[k]
                  + f_3 * pc_z[k] * sif_26[k];

        t_42[k] = f_4 * sid0_17[k]
                  - f_5 * sid1_17[k]
                  + f_3 * pc_y[k] * sif_28[k];

        t_43[k] = f_3 * pc_y[k] * sif_29[k];

        t_44[k] = f_7 * shf_9[k]
                  + f_1 * sid0_17[k]
                  - f_2 * sid1_17[k]
                  + f_3 * pc_z[k] * sif_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, shf_10, shf_30, shf_33, \
                         sid0_18, sid0_21, sid1_18, sid1_21, sif_30, \
                         sif_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_10 * shf_30[k]
                  + f_1 * sid0_18[k]
                  - f_2 * sid1_18[k]
                  + f_3 * pc_x[k] * sif_30[k];

        t_46[k] = f_8 * shf_10[k]
                  + f_3 * pc_y[k] * sif_30[k];

        t_47[k] = f_3 * pc_z[k] * sif_30[k];

        t_48[k] = f_10 * shf_33[k]
                  + f_4 * sid0_21[k]
                  - f_5 * sid1_21[k]
                  + f_3 * pc_x[k] * sif_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, shf_12, shf_35, shf_36, shf_37, \
                         sid0_23, sid1_23, sif_32, sif_35, sif_36, \
                         sif_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * shf_12[k]
                  + f_3 * pc_y[k] * sif_32[k];

        t_50[k] = f_10 * shf_35[k]
                  + f_4 * sid0_23[k]
                  - f_5 * sid1_23[k]
                  + f_3 * pc_x[k] * sif_35[k];

        t_51[k] = f_10 * shf_36[k]
                  + f_3 * pc_x[k] * sif_36[k];

        t_52[k] = f_10 * shf_37[k]
                  + f_3 * pc_x[k] * sif_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, shf_16, shf_38, shf_39, \
                         sid0_21, sid1_21, sif_36, sif_38, sif_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * shf_38[k]
                  + f_3 * pc_x[k] * sif_38[k];

        t_54[k] = f_10 * shf_39[k]
                  + f_3 * pc_x[k] * sif_39[k];

        t_55[k] = f_8 * shf_16[k]
                  + f_1 * sid0_21[k]
                  - f_2 * sid1_21[k]
                  + f_3 * pc_y[k] * sif_36[k];

        t_56[k] = f_3 * pc_z[k] * sif_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, shg0_30, shf_18, shf_19, \
                         shg1_30, sid0_23, sid1_23, sif_38, sif_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * shf_18[k]
                  + f_4 * sid0_23[k]
                  - f_5 * sid1_23[k]
                  + f_3 * pc_y[k] * sif_38[k];

        t_58[k] = f_8 * shf_19[k]
                  + f_3 * pc_y[k] * sif_39[k];

        t_59[k] = f_1 * sid0_23[k]
                  - f_2 * sid1_23[k]
                  + f_3 * pc_z[k] * sif_39[k];

        t_60[k] = pb_y[k] * shg0_30[k]
                  - f_6 * pc_y[k] * shg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, shg0_18, shf_10, shf_20, \
                         shf_22, shg1_18, sif_40, sif_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * shf_20[k]
                  + f_3 * pc_y[k] * sif_40[k];

        t_62[k] = f_7 * shf_10[k]
                  + f_3 * pc_z[k] * sif_40[k];

        t_63[k] = pb_z[k] * shg0_18[k]
                  - f_6 * pc_z[k] * shg1_18[k];

        t_64[k] = f_7 * shf_22[k]
                  + f_3 * pc_y[k] * sif_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, shg0_35, shf_46, shf_47, \
                         shf_48, shg1_35, sif_46, sif_47, sif_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * shg0_35[k]
                  - f_6 * pc_y[k] * shg1_35[k];

        t_66[k] = f_10 * shf_46[k]
                  + f_3 * pc_x[k] * sif_46[k];

        t_67[k] = f_10 * shf_47[k]
                  + f_3 * pc_x[k] * sif_47[k];

        t_68[k] = f_10 * shf_48[k]
                  + f_3 * pc_x[k] * sif_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, pc_x, pc_z, shg0_25, shf_16, shf_49, shg1_25, \
                         sif_46, sif_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * shf_49[k]
                  + f_3 * pc_x[k] * sif_49[k];

        t_70[k] = pb_z[k] * shg0_25[k]
                  - f_6 * pc_z[k] * shg1_25[k];

        t_71[k] = f_7 * shf_16[k]
                  + f_3 * pc_z[k] * sif_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, pc_y, shg0_44, shf_28, shf_29, shg1_44, \
                         sid0_29, sid1_29, sif_48, sif_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * shf_28[k]
                  + f_4 * sid0_29[k]
                  - f_5 * sid1_29[k]
                  + f_3 * pc_y[k] * sif_48[k];

        t_73[k] = f_7 * shf_29[k]
                  + f_3 * pc_y[k] * sif_49[k];

        t_74[k] = pb_y[k] * shg0_44[k]
                  - f_6 * pc_y[k] * shg1_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, shf_20, shf_50, shf_53, \
                         sid0_30, sid0_33, sid1_30, sid1_33, sif_50, \
                         sif_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * shf_50[k]
                  + f_1 * sid0_30[k]
                  - f_2 * sid1_30[k]
                  + f_3 * pc_x[k] * sif_50[k];

        t_76[k] = f_3 * pc_y[k] * sif_50[k];

        t_77[k] = f_8 * shf_20[k]
                  + f_3 * pc_z[k] * sif_50[k];

        t_78[k] = f_10 * shf_53[k]
                  + f_4 * sid0_33[k]
                  - f_5 * sid1_33[k]
                  + f_3 * pc_x[k] * sif_53[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, shf_55, shf_56, shf_57, sid0_35, \
                         sid1_35, sif_52, sif_55, sif_56, sif_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_y[k] * sif_52[k];

        t_80[k] = f_10 * shf_55[k]
                  + f_4 * sid0_35[k]
                  - f_5 * sid1_35[k]
                  + f_3 * pc_x[k] * sif_55[k];

        t_81[k] = f_10 * shf_56[k]
                  + f_3 * pc_x[k] * sif_56[k];

        t_82[k] = f_10 * shf_57[k]
                  + f_3 * pc_x[k] * sif_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, shf_26, shf_58, shf_59, \
                         sid0_33, sid1_33, sif_56, sif_58, sif_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_10 * shf_58[k]
                  + f_3 * pc_x[k] * sif_58[k];

        t_84[k] = f_10 * shf_59[k]
                  + f_3 * pc_x[k] * sif_59[k];

        t_85[k] = f_1 * sid0_33[k]
                  - f_2 * sid1_33[k]
                  + f_3 * pc_y[k] * sif_56[k];

        t_86[k] = f_8 * shf_26[k]
                  + f_3 * pc_z[k] * sif_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, shf_29, shf_60, sid0_35, \
                         sid0_36, sid1_35, sid1_36, sif_58, sif_59, \
                         sif_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * sid0_35[k]
                  - f_5 * sid1_35[k]
                  + f_3 * pc_y[k] * sif_58[k];

        t_88[k] = f_3 * pc_y[k] * sif_59[k];

        t_89[k] = f_8 * shf_29[k]
                  + f_1 * sid0_35[k]
                  - f_2 * sid1_35[k]
                  + f_3 * pc_z[k] * sif_59[k];

        t_90[k] = f_11 * shf_60[k]
                  + f_1 * sid0_36[k]
                  - f_2 * sid1_36[k]
                  + f_3 * pc_x[k] * sif_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, shf_30, shf_32, shf_63, \
                         sid0_39, sid1_39, sif_60, sif_62, sif_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_11 * shf_30[k]
                  + f_3 * pc_y[k] * sif_60[k];

        t_92[k] = f_3 * pc_z[k] * sif_60[k];

        t_93[k] = f_11 * shf_63[k]
                  + f_4 * sid0_39[k]
                  - f_5 * sid1_39[k]
                  + f_3 * pc_x[k] * sif_63[k];

        t_94[k] = f_11 * shf_32[k]
                  + f_3 * pc_y[k] * sif_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, shf_65, shf_66, shf_67, shf_68, \
                         sid0_41, sid1_41, sif_65, sif_66, sif_67, \
                         sif_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * shf_65[k]
                  + f_4 * sid0_41[k]
                  - f_5 * sid1_41[k]
                  + f_3 * pc_x[k] * sif_65[k];

        t_96[k] = f_11 * shf_66[k]
                  + f_3 * pc_x[k] * sif_66[k];

        t_97[k] = f_11 * shf_67[k]
                  + f_3 * pc_x[k] * sif_67[k];

        t_98[k] = f_11 * shf_68[k]
                  + f_3 * pc_x[k] * sif_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_x, pc_y, pc_z, shf_36, shf_69, sid0_39, \
                         sid1_39, sif_66, sif_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_11 * shf_69[k]
                  + f_3 * pc_x[k] * sif_69[k];

        t_100[k] = f_11 * shf_36[k]
                   + f_1 * sid0_39[k]
                   - f_2 * sid1_39[k]
                   + f_3 * pc_y[k] * sif_66[k];

        t_101[k] = f_3 * pc_z[k] * sif_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, shg0_45, shf_38, \
                         shf_39, shg1_45, sid0_41, sid1_41, sif_68, \
                         sif_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_11 * shf_38[k]
                   + f_4 * sid0_41[k]
                   - f_5 * sid1_41[k]
                   + f_3 * pc_y[k] * sif_68[k];

        t_103[k] = f_11 * shf_39[k]
                   + f_3 * pc_y[k] * sif_69[k];

        t_104[k] = f_1 * sid0_41[k]
                   - f_2 * sid1_41[k]
                   + f_3 * pc_z[k] * sif_69[k];

        t_105[k] = pb_z[k] * shg0_45[k]
                   - f_6 * pc_z[k] * shg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, shg0_48, shf_30, \
                         shf_40, shf_42, shg1_48, sif_70, sif_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * shf_40[k]
                   + f_3 * pc_y[k] * sif_70[k];

        t_107[k] = f_7 * shf_30[k]
                   + f_3 * pc_z[k] * sif_70[k];

        t_108[k] = pb_z[k] * shg0_48[k]
                   - f_6 * pc_z[k] * shg1_48[k];

        t_109[k] = f_8 * shf_42[k]
                   + f_3 * pc_y[k] * sif_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, shf_75, shf_76, shf_77, shf_78, \
                         sid0_47, sid1_47, sif_75, sif_76, sif_77, \
                         sif_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_11 * shf_75[k]
                   + f_4 * sid0_47[k]
                   - f_5 * sid1_47[k]
                   + f_3 * pc_x[k] * sif_75[k];

        t_111[k] = f_11 * shf_76[k]
                   + f_3 * pc_x[k] * sif_76[k];

        t_112[k] = f_11 * shf_77[k]
                   + f_3 * pc_x[k] * sif_77[k];

        t_113[k] = f_11 * shf_78[k]
                   + f_3 * pc_x[k] * sif_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_z, pc_x, pc_z, shg0_55, shf_36, shf_79, \
                         shg1_55, sif_76, sif_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_11 * shf_79[k]
                   + f_3 * pc_x[k] * sif_79[k];

        t_115[k] = pb_z[k] * shg0_55[k]
                   - f_6 * pc_z[k] * shg1_55[k];

        t_116[k] = f_7 * shf_36[k]
                   + f_3 * pc_z[k] * sif_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pc_y, pc_z, shg0_75, shf_39, \
                         shf_48, shf_49, shg1_75, sid0_47, sid1_47, sif_78, \
                         sif_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * shf_48[k]
                   + f_4 * sid0_47[k]
                   - f_5 * sid1_47[k]
                   + f_3 * pc_y[k] * sif_78[k];

        t_118[k] = f_8 * shf_49[k]
                   + f_3 * pc_y[k] * sif_79[k];

        t_119[k] = f_7 * shf_39[k]
                   + f_1 * sid0_47[k]
                   - f_2 * sid1_47[k]
                   + f_3 * pc_z[k] * sif_79[k];

        t_120[k] = pb_y[k] * shg0_75[k]
                   - f_6 * pc_y[k] * shg1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pc_y, pc_z, shg0_78, shf_40, \
                         shf_50, shf_51, shf_52, shg1_78, sif_80, \
                         sif_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * shf_50[k]
                   + f_3 * pc_y[k] * sif_80[k];

        t_122[k] = f_8 * shf_40[k]
                   + f_3 * pc_z[k] * sif_80[k];

        t_123[k] = pb_y[k] * shg0_78[k]
                   + f_8 * shf_51[k]
                   - f_6 * pc_y[k] * shg1_78[k];

        t_124[k] = f_7 * shf_52[k]
                   + f_3 * pc_y[k] * sif_82[k];
    }
}

static auto
compute_prim_sig_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shg0,
                                                          const size_t shf, const size_t shg1,
                                                          const size_t sid0, const size_t sid1,
                                                          const size_t sif, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 2.0 / q;
    const auto f_11 = 1.5 / q;

    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
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
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shg0_80 = buffer.data(shg0 + 80);
    const auto *shg0_89 = buffer.data(shg0 + 89);
    const auto *shg0_90 = buffer.data(shg0 + 90);
    const auto *shg0_93 = buffer.data(shg0 + 93);
    const auto *shg0_100 = buffer.data(shg0 + 100);
    const auto *shg0_135 = buffer.data(shg0 + 135);
    const auto *shg0_138 = buffer.data(shg0 + 138);
    const auto *shg0_140 = buffer.data(shg0 + 140);
    const auto *shg0_149 = buffer.data(shg0 + 149);
    const auto *shg0_150 = buffer.data(shg0 + 150);
    const auto *shg0_153 = buffer.data(shg0 + 153);
    const auto *shg0_225 = buffer.data(shg0 + 225);
    const auto *shg0_228 = buffer.data(shg0 + 228);
    const auto *shg0_230 = buffer.data(shg0 + 230);
    const auto *shg0_235 = buffer.data(shg0 + 235);
    const auto *shg0_237 = buffer.data(shg0 + 237);
    const auto *shg0_239 = buffer.data(shg0 + 239);

    const auto *shf_46 = buffer.data(shf + 46);
    const auto *shf_50 = buffer.data(shf + 50);
    const auto *shf_56 = buffer.data(shf + 56);
    const auto *shf_58 = buffer.data(shf + 58);
    const auto *shf_59 = buffer.data(shf + 59);
    const auto *shf_60 = buffer.data(shf + 60);
    const auto *shf_62 = buffer.data(shf + 62);
    const auto *shf_66 = buffer.data(shf + 66);
    const auto *shf_68 = buffer.data(shf + 68);
    const auto *shf_69 = buffer.data(shf + 69);
    const auto *shf_70 = buffer.data(shf + 70);
    const auto *shf_72 = buffer.data(shf + 72);
    const auto *shf_76 = buffer.data(shf + 76);
    const auto *shf_78 = buffer.data(shf + 78);
    const auto *shf_79 = buffer.data(shf + 79);
    const auto *shf_80 = buffer.data(shf + 80);
    const auto *shf_82 = buffer.data(shf + 82);
    const auto *shf_86 = buffer.data(shf + 86);
    const auto *shf_87 = buffer.data(shf + 87);
    const auto *shf_88 = buffer.data(shf + 88);
    const auto *shf_89 = buffer.data(shf + 89);
    const auto *shf_90 = buffer.data(shf + 90);
    const auto *shf_91 = buffer.data(shf + 91);
    const auto *shf_92 = buffer.data(shf + 92);
    const auto *shf_93 = buffer.data(shf + 93);
    const auto *shf_95 = buffer.data(shf + 95);
    const auto *shf_96 = buffer.data(shf + 96);
    const auto *shf_97 = buffer.data(shf + 97);
    const auto *shf_98 = buffer.data(shf + 98);
    const auto *shf_99 = buffer.data(shf + 99);
    const auto *shf_100 = buffer.data(shf + 100);
    const auto *shf_102 = buffer.data(shf + 102);
    const auto *shf_103 = buffer.data(shf + 103);
    const auto *shf_105 = buffer.data(shf + 105);
    const auto *shf_106 = buffer.data(shf + 106);
    const auto *shf_107 = buffer.data(shf + 107);
    const auto *shf_108 = buffer.data(shf + 108);
    const auto *shf_109 = buffer.data(shf + 109);
    const auto *shf_110 = buffer.data(shf + 110);
    const auto *shf_112 = buffer.data(shf + 112);
    const auto *shf_115 = buffer.data(shf + 115);
    const auto *shf_116 = buffer.data(shf + 116);
    const auto *shf_117 = buffer.data(shf + 117);
    const auto *shf_118 = buffer.data(shf + 118);
    const auto *shf_119 = buffer.data(shf + 119);
    const auto *shf_120 = buffer.data(shf + 120);
    const auto *shf_123 = buffer.data(shf + 123);
    const auto *shf_125 = buffer.data(shf + 125);
    const auto *shf_126 = buffer.data(shf + 126);
    const auto *shf_127 = buffer.data(shf + 127);
    const auto *shf_128 = buffer.data(shf + 128);
    const auto *shf_129 = buffer.data(shf + 129);
    const auto *shf_136 = buffer.data(shf + 136);
    const auto *shf_137 = buffer.data(shf + 137);
    const auto *shf_138 = buffer.data(shf + 138);
    const auto *shf_139 = buffer.data(shf + 139);
    const auto *shf_140 = buffer.data(shf + 140);
    const auto *shf_143 = buffer.data(shf + 143);
    const auto *shf_145 = buffer.data(shf + 145);
    const auto *shf_146 = buffer.data(shf + 146);
    const auto *shf_147 = buffer.data(shf + 147);
    const auto *shf_148 = buffer.data(shf + 148);
    const auto *shf_149 = buffer.data(shf + 149);
    const auto *shf_150 = buffer.data(shf + 150);
    const auto *shf_153 = buffer.data(shf + 153);
    const auto *shf_155 = buffer.data(shf + 155);
    const auto *shf_156 = buffer.data(shf + 156);
    const auto *shf_157 = buffer.data(shf + 157);
    const auto *shf_158 = buffer.data(shf + 158);
    const auto *shf_159 = buffer.data(shf + 159);

    const auto *shg1_80 = buffer.data(shg1 + 80);
    const auto *shg1_89 = buffer.data(shg1 + 89);
    const auto *shg1_90 = buffer.data(shg1 + 90);
    const auto *shg1_93 = buffer.data(shg1 + 93);
    const auto *shg1_100 = buffer.data(shg1 + 100);
    const auto *shg1_135 = buffer.data(shg1 + 135);
    const auto *shg1_138 = buffer.data(shg1 + 138);
    const auto *shg1_140 = buffer.data(shg1 + 140);
    const auto *shg1_149 = buffer.data(shg1 + 149);
    const auto *shg1_150 = buffer.data(shg1 + 150);
    const auto *shg1_153 = buffer.data(shg1 + 153);
    const auto *shg1_225 = buffer.data(shg1 + 225);
    const auto *shg1_228 = buffer.data(shg1 + 228);
    const auto *shg1_230 = buffer.data(shg1 + 230);
    const auto *shg1_235 = buffer.data(shg1 + 235);
    const auto *shg1_237 = buffer.data(shg1 + 237);
    const auto *shg1_239 = buffer.data(shg1 + 239);

    const auto *sid0_51 = buffer.data(sid0 + 51);
    const auto *sid0_53 = buffer.data(sid0 + 53);
    const auto *sid0_54 = buffer.data(sid0 + 54);
    const auto *sid0_57 = buffer.data(sid0 + 57);
    const auto *sid0_59 = buffer.data(sid0 + 59);
    const auto *sid0_60 = buffer.data(sid0 + 60);
    const auto *sid0_63 = buffer.data(sid0 + 63);
    const auto *sid0_65 = buffer.data(sid0 + 65);
    const auto *sid0_71 = buffer.data(sid0 + 71);
    const auto *sid0_72 = buffer.data(sid0 + 72);
    const auto *sid0_75 = buffer.data(sid0 + 75);
    const auto *sid0_77 = buffer.data(sid0 + 77);
    const auto *sid0_81 = buffer.data(sid0 + 81);
    const auto *sid0_83 = buffer.data(sid0 + 83);
    const auto *sid0_84 = buffer.data(sid0 + 84);
    const auto *sid0_87 = buffer.data(sid0 + 87);
    const auto *sid0_89 = buffer.data(sid0 + 89);

    const auto *sid1_51 = buffer.data(sid1 + 51);
    const auto *sid1_53 = buffer.data(sid1 + 53);
    const auto *sid1_54 = buffer.data(sid1 + 54);
    const auto *sid1_57 = buffer.data(sid1 + 57);
    const auto *sid1_59 = buffer.data(sid1 + 59);
    const auto *sid1_60 = buffer.data(sid1 + 60);
    const auto *sid1_63 = buffer.data(sid1 + 63);
    const auto *sid1_65 = buffer.data(sid1 + 65);
    const auto *sid1_71 = buffer.data(sid1 + 71);
    const auto *sid1_72 = buffer.data(sid1 + 72);
    const auto *sid1_75 = buffer.data(sid1 + 75);
    const auto *sid1_77 = buffer.data(sid1 + 77);
    const auto *sid1_81 = buffer.data(sid1 + 81);
    const auto *sid1_83 = buffer.data(sid1 + 83);
    const auto *sid1_84 = buffer.data(sid1 + 84);
    const auto *sid1_87 = buffer.data(sid1 + 87);
    const auto *sid1_89 = buffer.data(sid1 + 89);

    const auto *sif_86 = buffer.data(sif + 86);
    const auto *sif_87 = buffer.data(sif + 87);
    const auto *sif_88 = buffer.data(sif + 88);
    const auto *sif_89 = buffer.data(sif + 89);
    const auto *sif_90 = buffer.data(sif + 90);
    const auto *sif_92 = buffer.data(sif + 92);
    const auto *sif_93 = buffer.data(sif + 93);
    const auto *sif_95 = buffer.data(sif + 95);
    const auto *sif_96 = buffer.data(sif + 96);
    const auto *sif_97 = buffer.data(sif + 97);
    const auto *sif_98 = buffer.data(sif + 98);
    const auto *sif_99 = buffer.data(sif + 99);
    const auto *sif_100 = buffer.data(sif + 100);
    const auto *sif_102 = buffer.data(sif + 102);
    const auto *sif_103 = buffer.data(sif + 103);
    const auto *sif_105 = buffer.data(sif + 105);
    const auto *sif_106 = buffer.data(sif + 106);
    const auto *sif_107 = buffer.data(sif + 107);
    const auto *sif_108 = buffer.data(sif + 108);
    const auto *sif_109 = buffer.data(sif + 109);
    const auto *sif_110 = buffer.data(sif + 110);
    const auto *sif_112 = buffer.data(sif + 112);
    const auto *sif_115 = buffer.data(sif + 115);
    const auto *sif_116 = buffer.data(sif + 116);
    const auto *sif_117 = buffer.data(sif + 117);
    const auto *sif_118 = buffer.data(sif + 118);
    const auto *sif_119 = buffer.data(sif + 119);
    const auto *sif_120 = buffer.data(sif + 120);
    const auto *sif_122 = buffer.data(sif + 122);
    const auto *sif_123 = buffer.data(sif + 123);
    const auto *sif_125 = buffer.data(sif + 125);
    const auto *sif_126 = buffer.data(sif + 126);
    const auto *sif_127 = buffer.data(sif + 127);
    const auto *sif_128 = buffer.data(sif + 128);
    const auto *sif_129 = buffer.data(sif + 129);
    const auto *sif_130 = buffer.data(sif + 130);
    const auto *sif_132 = buffer.data(sif + 132);
    const auto *sif_136 = buffer.data(sif + 136);
    const auto *sif_137 = buffer.data(sif + 137);
    const auto *sif_138 = buffer.data(sif + 138);
    const auto *sif_139 = buffer.data(sif + 139);
    const auto *sif_140 = buffer.data(sif + 140);
    const auto *sif_142 = buffer.data(sif + 142);
    const auto *sif_143 = buffer.data(sif + 143);
    const auto *sif_145 = buffer.data(sif + 145);
    const auto *sif_146 = buffer.data(sif + 146);
    const auto *sif_147 = buffer.data(sif + 147);
    const auto *sif_148 = buffer.data(sif + 148);
    const auto *sif_149 = buffer.data(sif + 149);
    const auto *sif_150 = buffer.data(sif + 150);
    const auto *sif_152 = buffer.data(sif + 152);
    const auto *sif_156 = buffer.data(sif + 156);
    const auto *sif_157 = buffer.data(sif + 157);
    const auto *sif_158 = buffer.data(sif + 158);
    const auto *sif_159 = buffer.data(sif + 159);
    const auto *sif_160 = buffer.data(sif + 160);
    const auto *sif_162 = buffer.data(sif + 162);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, shg0_80, shf_86, \
                         shf_87, shf_88, shg1_80, sif_86, sif_87, \
                         sif_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * shg0_80[k]
                   - f_6 * pc_y[k] * shg1_80[k];

        t_126[k] = f_11 * shf_86[k]
                   + f_3 * pc_x[k] * sif_86[k];

        t_127[k] = f_11 * shf_87[k]
                   + f_3 * pc_x[k] * sif_87[k];

        t_128[k] = f_11 * shf_88[k]
                   + f_3 * pc_x[k] * sif_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, pc_z, shf_46, shf_56, shf_89, \
                         sid0_51, sid1_51, sif_86, sif_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * shf_89[k]
                   + f_3 * pc_x[k] * sif_89[k];

        t_130[k] = f_7 * shf_56[k]
                   + f_1 * sid0_51[k]
                   - f_2 * sid1_51[k]
                   + f_3 * pc_y[k] * sif_86[k];

        t_131[k] = f_8 * shf_46[k]
                   + f_3 * pc_z[k] * sif_86[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pc_y, shg0_89, shf_58, shf_59, shg1_89, \
                         sid0_53, sid1_53, sif_88, sif_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_7 * shf_58[k]
                   + f_4 * sid0_53[k]
                   - f_5 * sid1_53[k]
                   + f_3 * pc_y[k] * sif_88[k];

        t_133[k] = f_7 * shf_59[k]
                   + f_3 * pc_y[k] * sif_89[k];

        t_134[k] = pb_y[k] * shg0_89[k]
                   - f_6 * pc_y[k] * shg1_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, shf_50, shf_90, shf_93, \
                         sid0_54, sid0_57, sid1_54, sid1_57, sif_90, \
                         sif_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_11 * shf_90[k]
                   + f_1 * sid0_54[k]
                   - f_2 * sid1_54[k]
                   + f_3 * pc_x[k] * sif_90[k];

        t_136[k] = f_3 * pc_y[k] * sif_90[k];

        t_137[k] = f_11 * shf_50[k]
                   + f_3 * pc_z[k] * sif_90[k];

        t_138[k] = f_11 * shf_93[k]
                   + f_4 * sid0_57[k]
                   - f_5 * sid1_57[k]
                   + f_3 * pc_x[k] * sif_93[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, shf_95, shf_96, shf_97, \
                         sid0_59, sid1_59, sif_92, sif_95, sif_96, \
                         sif_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * pc_y[k] * sif_92[k];

        t_140[k] = f_11 * shf_95[k]
                   + f_4 * sid0_59[k]
                   - f_5 * sid1_59[k]
                   + f_3 * pc_x[k] * sif_95[k];

        t_141[k] = f_11 * shf_96[k]
                   + f_3 * pc_x[k] * sif_96[k];

        t_142[k] = f_11 * shf_97[k]
                   + f_3 * pc_x[k] * sif_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, shf_56, shf_98, shf_99, \
                         sid0_57, sid1_57, sif_96, sif_98, sif_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * shf_98[k]
                   + f_3 * pc_x[k] * sif_98[k];

        t_144[k] = f_11 * shf_99[k]
                   + f_3 * pc_x[k] * sif_99[k];

        t_145[k] = f_1 * sid0_57[k]
                   - f_2 * sid1_57[k]
                   + f_3 * pc_y[k] * sif_96[k];

        t_146[k] = f_11 * shf_56[k]
                   + f_3 * pc_z[k] * sif_96[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, shf_59, shf_100, \
                         sid0_59, sid0_60, sid1_59, sid1_60, sif_98, sif_99, \
                         sif_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * sid0_59[k]
                   - f_5 * sid1_59[k]
                   + f_3 * pc_y[k] * sif_98[k];

        t_148[k] = f_3 * pc_y[k] * sif_99[k];

        t_149[k] = f_11 * shf_59[k]
                   + f_1 * sid0_59[k]
                   - f_2 * sid1_59[k]
                   + f_3 * pc_z[k] * sif_99[k];

        t_150[k] = f_8 * shf_100[k]
                   + f_1 * sid0_60[k]
                   - f_2 * sid1_60[k]
                   + f_3 * pc_x[k] * sif_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, shf_60, shf_62, \
                         shf_103, sid0_63, sid1_63, sif_100, sif_102, \
                         sif_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_10 * shf_60[k]
                   + f_3 * pc_y[k] * sif_100[k];

        t_152[k] = f_3 * pc_z[k] * sif_100[k];

        t_153[k] = f_8 * shf_103[k]
                   + f_4 * sid0_63[k]
                   - f_5 * sid1_63[k]
                   + f_3 * pc_x[k] * sif_103[k];

        t_154[k] = f_10 * shf_62[k]
                   + f_3 * pc_y[k] * sif_102[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, shf_105, shf_106, shf_107, shf_108, \
                         sid0_65, sid1_65, sif_105, sif_106, sif_107, \
                         sif_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_8 * shf_105[k]
                   + f_4 * sid0_65[k]
                   - f_5 * sid1_65[k]
                   + f_3 * pc_x[k] * sif_105[k];

        t_156[k] = f_8 * shf_106[k]
                   + f_3 * pc_x[k] * sif_106[k];

        t_157[k] = f_8 * shf_107[k]
                   + f_3 * pc_x[k] * sif_107[k];

        t_158[k] = f_8 * shf_108[k]
                   + f_3 * pc_x[k] * sif_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, pc_z, shf_66, shf_109, sid0_63, \
                         sid1_63, sif_106, sif_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_8 * shf_109[k]
                   + f_3 * pc_x[k] * sif_109[k];

        t_160[k] = f_10 * shf_66[k]
                   + f_1 * sid0_63[k]
                   - f_2 * sid1_63[k]
                   + f_3 * pc_y[k] * sif_106[k];

        t_161[k] = f_3 * pc_z[k] * sif_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_z, pc_y, pc_z, shg0_90, shf_68, \
                         shf_69, shg1_90, sid0_65, sid1_65, sif_108, \
                         sif_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_10 * shf_68[k]
                   + f_4 * sid0_65[k]
                   - f_5 * sid1_65[k]
                   + f_3 * pc_y[k] * sif_108[k];

        t_163[k] = f_10 * shf_69[k]
                   + f_3 * pc_y[k] * sif_109[k];

        t_164[k] = f_1 * sid0_65[k]
                   - f_2 * sid1_65[k]
                   + f_3 * pc_z[k] * sif_109[k];

        t_165[k] = pb_z[k] * shg0_90[k]
                   - f_6 * pc_z[k] * shg1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, shg0_93, shf_60, \
                         shf_70, shf_72, shg1_93, sif_110, sif_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_11 * shf_70[k]
                   + f_3 * pc_y[k] * sif_110[k];

        t_167[k] = f_7 * shf_60[k]
                   + f_3 * pc_z[k] * sif_110[k];

        t_168[k] = pb_z[k] * shg0_93[k]
                   - f_6 * pc_z[k] * shg1_93[k];

        t_169[k] = f_11 * shf_72[k]
                   + f_3 * pc_y[k] * sif_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, shf_115, shf_116, shf_117, shf_118, \
                         sid0_71, sid1_71, sif_115, sif_116, sif_117, \
                         sif_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_8 * shf_115[k]
                   + f_4 * sid0_71[k]
                   - f_5 * sid1_71[k]
                   + f_3 * pc_x[k] * sif_115[k];

        t_171[k] = f_8 * shf_116[k]
                   + f_3 * pc_x[k] * sif_116[k];

        t_172[k] = f_8 * shf_117[k]
                   + f_3 * pc_x[k] * sif_117[k];

        t_173[k] = f_8 * shf_118[k]
                   + f_3 * pc_x[k] * sif_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_z, pc_x, pc_z, shg0_100, shf_66, shf_119, \
                         shg1_100, sif_116, sif_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_8 * shf_119[k]
                   + f_3 * pc_x[k] * sif_119[k];

        t_175[k] = pb_z[k] * shg0_100[k]
                   - f_6 * pc_z[k] * shg1_100[k];

        t_176[k] = f_7 * shf_66[k]
                   + f_3 * pc_z[k] * sif_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pc_y, pc_z, shf_69, shf_78, shf_79, sid0_71, \
                         sid1_71, sif_118, sif_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_11 * shf_78[k]
                   + f_4 * sid0_71[k]
                   - f_5 * sid1_71[k]
                   + f_3 * pc_y[k] * sif_118[k];

        t_178[k] = f_11 * shf_79[k]
                   + f_3 * pc_y[k] * sif_119[k];

        t_179[k] = f_7 * shf_69[k]
                   + f_1 * sid0_71[k]
                   - f_2 * sid1_71[k]
                   + f_3 * pc_z[k] * sif_119[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pc_x, pc_y, pc_z, shf_70, shf_80, shf_120, \
                         sid0_72, sid1_72, sif_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_8 * shf_120[k]
                   + f_1 * sid0_72[k]
                   - f_2 * sid1_72[k]
                   + f_3 * pc_x[k] * sif_120[k];

        t_181[k] = f_8 * shf_80[k]
                   + f_3 * pc_y[k] * sif_120[k];

        t_182[k] = f_8 * shf_70[k]
                   + f_3 * pc_z[k] * sif_120[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_x, pc_y, shf_82, shf_123, shf_125, sid0_75, \
                         sid0_77, sid1_75, sid1_77, sif_122, sif_123, \
                         sif_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_8 * shf_123[k]
                   + f_4 * sid0_75[k]
                   - f_5 * sid1_75[k]
                   + f_3 * pc_x[k] * sif_123[k];

        t_184[k] = f_8 * shf_82[k]
                   + f_3 * pc_y[k] * sif_122[k];

        t_185[k] = f_8 * shf_125[k]
                   + f_4 * sid0_77[k]
                   - f_5 * sid1_77[k]
                   + f_3 * pc_x[k] * sif_125[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, shf_126, shf_127, shf_128, shf_129, \
                         sif_126, sif_127, sif_128, sif_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_8 * shf_126[k]
                   + f_3 * pc_x[k] * sif_126[k];

        t_187[k] = f_8 * shf_127[k]
                   + f_3 * pc_x[k] * sif_127[k];

        t_188[k] = f_8 * shf_128[k]
                   + f_3 * pc_x[k] * sif_128[k];

        t_189[k] = f_8 * shf_129[k]
                   + f_3 * pc_x[k] * sif_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_y, pc_z, shf_76, shf_86, shf_88, sid0_75, \
                         sid0_77, sid1_75, sid1_77, sif_126, sif_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * shf_86[k]
                   + f_1 * sid0_75[k]
                   - f_2 * sid1_75[k]
                   + f_3 * pc_y[k] * sif_126[k];

        t_191[k] = f_8 * shf_76[k]
                   + f_3 * pc_z[k] * sif_126[k];

        t_192[k] = f_8 * shf_88[k]
                   + f_4 * sid0_77[k]
                   - f_5 * sid1_77[k]
                   + f_3 * pc_y[k] * sif_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_y, pc_y, pc_z, shg0_135, shf_79, \
                         shf_89, shf_90, shg1_135, sid0_77, sid1_77, sif_129, \
                         sif_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * shf_89[k]
                   + f_3 * pc_y[k] * sif_129[k];

        t_194[k] = f_8 * shf_79[k]
                   + f_1 * sid0_77[k]
                   - f_2 * sid1_77[k]
                   + f_3 * pc_z[k] * sif_129[k];

        t_195[k] = pb_y[k] * shg0_135[k]
                   - f_6 * pc_y[k] * shg1_135[k];

        t_196[k] = f_7 * shf_90[k]
                   + f_3 * pc_y[k] * sif_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_y, pc_y, pc_z, shg0_138, shg0_140, \
                         shf_80, shf_91, shf_92, shg1_138, shg1_140, sif_130, \
                         sif_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_11 * shf_80[k]
                   + f_3 * pc_z[k] * sif_130[k];

        t_198[k] = pb_y[k] * shg0_138[k]
                   + f_8 * shf_91[k]
                   - f_6 * pc_y[k] * shg1_138[k];

        t_199[k] = f_7 * shf_92[k]
                   + f_3 * pc_y[k] * sif_132[k];

        t_200[k] = pb_y[k] * shg0_140[k]
                   - f_6 * pc_y[k] * shg1_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, shf_136, shf_137, shf_138, shf_139, \
                         sif_136, sif_137, sif_138, sif_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_8 * shf_136[k]
                   + f_3 * pc_x[k] * sif_136[k];

        t_202[k] = f_8 * shf_137[k]
                   + f_3 * pc_x[k] * sif_137[k];

        t_203[k] = f_8 * shf_138[k]
                   + f_3 * pc_x[k] * sif_138[k];

        t_204[k] = f_8 * shf_139[k]
                   + f_3 * pc_x[k] * sif_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pc_y, pc_z, shf_86, shf_96, shf_98, sid0_81, \
                         sid0_83, sid1_81, sid1_83, sif_136, sif_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_7 * shf_96[k]
                   + f_1 * sid0_81[k]
                   - f_2 * sid1_81[k]
                   + f_3 * pc_y[k] * sif_136[k];

        t_206[k] = f_11 * shf_86[k]
                   + f_3 * pc_z[k] * sif_136[k];

        t_207[k] = f_7 * shf_98[k]
                   + f_4 * sid0_83[k]
                   - f_5 * sid1_83[k]
                   + f_3 * pc_y[k] * sif_138[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pc_x, pc_y, shg0_149, shf_99, \
                         shf_140, shg1_149, sid0_84, sid1_84, sif_139, \
                         sif_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * shf_99[k]
                   + f_3 * pc_y[k] * sif_139[k];

        t_209[k] = pb_y[k] * shg0_149[k]
                   - f_6 * pc_y[k] * shg1_149[k];

        t_210[k] = f_8 * shf_140[k]
                   + f_1 * sid0_84[k]
                   - f_2 * sid1_84[k]
                   + f_3 * pc_x[k] * sif_140[k];

        t_211[k] = f_3 * pc_y[k] * sif_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, shf_90, shf_143, sid0_87, \
                         sid1_87, sif_140, sif_142, sif_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_10 * shf_90[k]
                   + f_3 * pc_z[k] * sif_140[k];

        t_213[k] = f_8 * shf_143[k]
                   + f_4 * sid0_87[k]
                   - f_5 * sid1_87[k]
                   + f_3 * pc_x[k] * sif_143[k];

        t_214[k] = f_3 * pc_y[k] * sif_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, shf_145, shf_146, shf_147, shf_148, \
                         sid0_89, sid1_89, sif_145, sif_146, sif_147, \
                         sif_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_8 * shf_145[k]
                   + f_4 * sid0_89[k]
                   - f_5 * sid1_89[k]
                   + f_3 * pc_x[k] * sif_145[k];

        t_216[k] = f_8 * shf_146[k]
                   + f_3 * pc_x[k] * sif_146[k];

        t_217[k] = f_8 * shf_147[k]
                   + f_3 * pc_x[k] * sif_147[k];

        t_218[k] = f_8 * shf_148[k]
                   + f_3 * pc_x[k] * sif_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pc_x, pc_y, pc_z, shf_96, shf_149, \
                         sid0_87, sid0_89, sid1_87, sid1_89, sif_146, sif_148, \
                         sif_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_8 * shf_149[k]
                   + f_3 * pc_x[k] * sif_149[k];

        t_220[k] = f_1 * sid0_87[k]
                   - f_2 * sid1_87[k]
                   + f_3 * pc_y[k] * sif_146[k];

        t_221[k] = f_10 * shf_96[k]
                   + f_3 * pc_z[k] * sif_146[k];

        t_222[k] = f_4 * sid0_89[k]
                   - f_5 * sid1_89[k]
                   + f_3 * pc_y[k] * sif_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pb_x, pc_x, pc_y, pc_z, shg0_225, shf_99, \
                         shf_150, shg1_225, sid0_89, sid1_89, sif_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_y[k] * sif_149[k];

        t_224[k] = f_10 * shf_99[k]
                   + f_1 * sid0_89[k]
                   - f_2 * sid1_89[k]
                   + f_3 * pc_z[k] * sif_149[k];

        t_225[k] = pb_x[k] * shg0_225[k]
                   + f_10 * shf_150[k]
                   - f_6 * pc_x[k] * shg1_225[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pb_x, pc_x, pc_y, pc_z, shg0_228, \
                         shf_100, shf_102, shf_153, shg1_228, sif_150, \
                         sif_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_9 * shf_100[k]
                   + f_3 * pc_y[k] * sif_150[k];

        t_227[k] = f_3 * pc_z[k] * sif_150[k];

        t_228[k] = pb_x[k] * shg0_228[k]
                   + f_8 * shf_153[k]
                   - f_6 * pc_x[k] * shg1_228[k];

        t_229[k] = f_9 * shf_102[k]
                   + f_3 * pc_y[k] * sif_152[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pb_x, pc_x, shg0_230, shf_155, shf_156, \
                         shf_157, shf_158, shg1_230, sif_156, sif_157, \
                         sif_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = pb_x[k] * shg0_230[k]
                   + f_8 * shf_155[k]
                   - f_6 * pc_x[k] * shg1_230[k];

        t_231[k] = f_7 * shf_156[k]
                   + f_3 * pc_x[k] * sif_156[k];

        t_232[k] = f_7 * shf_157[k]
                   + f_3 * pc_x[k] * sif_157[k];

        t_233[k] = f_7 * shf_158[k]
                   + f_3 * pc_x[k] * sif_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pb_x, pc_x, pc_z, shg0_235, shg0_237, \
                         shf_159, shg1_235, shg1_237, sif_156, \
                         sif_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_7 * shf_159[k]
                   + f_3 * pc_x[k] * sif_159[k];

        t_235[k] = pb_x[k] * shg0_235[k]
                   - f_6 * pc_x[k] * shg1_235[k];

        t_236[k] = f_3 * pc_z[k] * sif_156[k];

        t_237[k] = pb_x[k] * shg0_237[k]
                   - f_6 * pc_x[k] * shg1_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pb_x, pb_z, pc_x, pc_y, pc_z, shg0_150, \
                         shg0_239, shf_109, shg1_150, shg1_239, \
                         sif_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_9 * shf_109[k]
                   + f_3 * pc_y[k] * sif_159[k];

        t_239[k] = pb_x[k] * shg0_239[k]
                   - f_6 * pc_x[k] * shg1_239[k];

        t_240[k] = pb_z[k] * shg0_150[k]
                   - f_6 * pc_z[k] * shg1_150[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_z, pc_y, pc_z, shg0_153, shf_100, \
                         shf_110, shf_112, shg1_153, sif_160, sif_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_10 * shf_110[k]
                   + f_3 * pc_y[k] * sif_160[k];

        t_242[k] = f_7 * shf_100[k]
                   + f_3 * pc_z[k] * sif_160[k];

        t_243[k] = pb_z[k] * shg0_153[k]
                   - f_6 * pc_z[k] * shg1_153[k];

        t_244[k] = f_10 * shf_112[k]
                   + f_3 * pc_y[k] * sif_162[k];
    }
}

static auto
compute_prim_sig_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shg0,
                                                          const size_t shf, const size_t shg1,
                                                          const size_t sid0, const size_t sid1,
                                                          const size_t sif, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 2.0 / q;
    const auto f_11 = 1.5 / q;

    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);
    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shg0_210 = buffer.data(shg0 + 210);
    const auto *shg0_215 = buffer.data(shg0 + 215);
    const auto *shg0_225 = buffer.data(shg0 + 225);
    const auto *shg0_228 = buffer.data(shg0 + 228);
    const auto *shg0_235 = buffer.data(shg0 + 235);
    const auto *shg0_237 = buffer.data(shg0 + 237);
    const auto *shg0_245 = buffer.data(shg0 + 245);
    const auto *shg0_250 = buffer.data(shg0 + 250);
    const auto *shg0_252 = buffer.data(shg0 + 252);
    const auto *shg0_254 = buffer.data(shg0 + 254);
    const auto *shg0_255 = buffer.data(shg0 + 255);
    const auto *shg0_258 = buffer.data(shg0 + 258);
    const auto *shg0_260 = buffer.data(shg0 + 260);
    const auto *shg0_265 = buffer.data(shg0 + 265);
    const auto *shg0_267 = buffer.data(shg0 + 267);
    const auto *shg0_269 = buffer.data(shg0 + 269);
    const auto *shg0_270 = buffer.data(shg0 + 270);
    const auto *shg0_273 = buffer.data(shg0 + 273);
    const auto *shg0_275 = buffer.data(shg0 + 275);
    const auto *shg0_280 = buffer.data(shg0 + 280);
    const auto *shg0_282 = buffer.data(shg0 + 282);
    const auto *shg0_284 = buffer.data(shg0 + 284);
    const auto *shg0_288 = buffer.data(shg0 + 288);
    const auto *shg0_295 = buffer.data(shg0 + 295);
    const auto *shg0_297 = buffer.data(shg0 + 297);
    const auto *shg0_299 = buffer.data(shg0 + 299);
    const auto *shg0_300 = buffer.data(shg0 + 300);
    const auto *shg0_303 = buffer.data(shg0 + 303);
    const auto *shg0_305 = buffer.data(shg0 + 305);
    const auto *shg0_310 = buffer.data(shg0 + 310);
    const auto *shg0_312 = buffer.data(shg0 + 312);
    const auto *shg0_314 = buffer.data(shg0 + 314);

    const auto *shf_106 = buffer.data(shf + 106);
    const auto *shf_110 = buffer.data(shf + 110);
    const auto *shf_116 = buffer.data(shf + 116);
    const auto *shf_119 = buffer.data(shf + 119);
    const auto *shf_120 = buffer.data(shf + 120);
    const auto *shf_122 = buffer.data(shf + 122);
    const auto *shf_126 = buffer.data(shf + 126);
    const auto *shf_129 = buffer.data(shf + 129);
    const auto *shf_130 = buffer.data(shf + 130);
    const auto *shf_132 = buffer.data(shf + 132);
    const auto *shf_136 = buffer.data(shf + 136);
    const auto *shf_139 = buffer.data(shf + 139);
    const auto *shf_140 = buffer.data(shf + 140);
    const auto *shf_142 = buffer.data(shf + 142);
    const auto *shf_146 = buffer.data(shf + 146);
    const auto *shf_149 = buffer.data(shf + 149);
    const auto *shf_150 = buffer.data(shf + 150);
    const auto *shf_152 = buffer.data(shf + 152);
    const auto *shf_156 = buffer.data(shf + 156);
    const auto *shf_157 = buffer.data(shf + 157);
    const auto *shf_158 = buffer.data(shf + 158);
    const auto *shf_159 = buffer.data(shf + 159);
    const auto *shf_160 = buffer.data(shf + 160);
    const auto *shf_162 = buffer.data(shf + 162);
    const auto *shf_165 = buffer.data(shf + 165);
    const auto *shf_166 = buffer.data(shf + 166);
    const auto *shf_167 = buffer.data(shf + 167);
    const auto *shf_168 = buffer.data(shf + 168);
    const auto *shf_169 = buffer.data(shf + 169);
    const auto *shf_170 = buffer.data(shf + 170);
    const auto *shf_172 = buffer.data(shf + 172);
    const auto *shf_173 = buffer.data(shf + 173);
    const auto *shf_175 = buffer.data(shf + 175);
    const auto *shf_176 = buffer.data(shf + 176);
    const auto *shf_177 = buffer.data(shf + 177);
    const auto *shf_178 = buffer.data(shf + 178);
    const auto *shf_179 = buffer.data(shf + 179);
    const auto *shf_180 = buffer.data(shf + 180);
    const auto *shf_182 = buffer.data(shf + 182);
    const auto *shf_183 = buffer.data(shf + 183);
    const auto *shf_185 = buffer.data(shf + 185);
    const auto *shf_186 = buffer.data(shf + 186);
    const auto *shf_187 = buffer.data(shf + 187);
    const auto *shf_188 = buffer.data(shf + 188);
    const auto *shf_189 = buffer.data(shf + 189);
    const auto *shf_193 = buffer.data(shf + 193);
    const auto *shf_196 = buffer.data(shf + 196);
    const auto *shf_197 = buffer.data(shf + 197);
    const auto *shf_198 = buffer.data(shf + 198);
    const auto *shf_199 = buffer.data(shf + 199);
    const auto *shf_200 = buffer.data(shf + 200);
    const auto *shf_203 = buffer.data(shf + 203);
    const auto *shf_205 = buffer.data(shf + 205);
    const auto *shf_206 = buffer.data(shf + 206);
    const auto *shf_207 = buffer.data(shf + 207);
    const auto *shf_208 = buffer.data(shf + 208);
    const auto *shf_209 = buffer.data(shf + 209);

    const auto *shg1_210 = buffer.data(shg1 + 210);
    const auto *shg1_215 = buffer.data(shg1 + 215);
    const auto *shg1_225 = buffer.data(shg1 + 225);
    const auto *shg1_228 = buffer.data(shg1 + 228);
    const auto *shg1_235 = buffer.data(shg1 + 235);
    const auto *shg1_237 = buffer.data(shg1 + 237);
    const auto *shg1_245 = buffer.data(shg1 + 245);
    const auto *shg1_250 = buffer.data(shg1 + 250);
    const auto *shg1_252 = buffer.data(shg1 + 252);
    const auto *shg1_254 = buffer.data(shg1 + 254);
    const auto *shg1_255 = buffer.data(shg1 + 255);
    const auto *shg1_258 = buffer.data(shg1 + 258);
    const auto *shg1_260 = buffer.data(shg1 + 260);
    const auto *shg1_265 = buffer.data(shg1 + 265);
    const auto *shg1_267 = buffer.data(shg1 + 267);
    const auto *shg1_269 = buffer.data(shg1 + 269);
    const auto *shg1_270 = buffer.data(shg1 + 270);
    const auto *shg1_273 = buffer.data(shg1 + 273);
    const auto *shg1_275 = buffer.data(shg1 + 275);
    const auto *shg1_280 = buffer.data(shg1 + 280);
    const auto *shg1_282 = buffer.data(shg1 + 282);
    const auto *shg1_284 = buffer.data(shg1 + 284);
    const auto *shg1_288 = buffer.data(shg1 + 288);
    const auto *shg1_295 = buffer.data(shg1 + 295);
    const auto *shg1_297 = buffer.data(shg1 + 297);
    const auto *shg1_299 = buffer.data(shg1 + 299);
    const auto *shg1_300 = buffer.data(shg1 + 300);
    const auto *shg1_303 = buffer.data(shg1 + 303);
    const auto *shg1_305 = buffer.data(shg1 + 305);
    const auto *shg1_310 = buffer.data(shg1 + 310);
    const auto *shg1_312 = buffer.data(shg1 + 312);
    const auto *shg1_314 = buffer.data(shg1 + 314);

    const auto *sid0_126 = buffer.data(sid0 + 126);
    const auto *sid0_129 = buffer.data(sid0 + 129);
    const auto *sid0_131 = buffer.data(sid0 + 131);
    const auto *sid0_137 = buffer.data(sid0 + 137);
    const auto *sid0_138 = buffer.data(sid0 + 138);
    const auto *sid0_141 = buffer.data(sid0 + 141);
    const auto *sid0_143 = buffer.data(sid0 + 143);
    const auto *sid0_144 = buffer.data(sid0 + 144);
    const auto *sid0_147 = buffer.data(sid0 + 147);
    const auto *sid0_149 = buffer.data(sid0 + 149);

    const auto *sid1_126 = buffer.data(sid1 + 126);
    const auto *sid1_129 = buffer.data(sid1 + 129);
    const auto *sid1_131 = buffer.data(sid1 + 131);
    const auto *sid1_137 = buffer.data(sid1 + 137);
    const auto *sid1_138 = buffer.data(sid1 + 138);
    const auto *sid1_141 = buffer.data(sid1 + 141);
    const auto *sid1_143 = buffer.data(sid1 + 143);
    const auto *sid1_144 = buffer.data(sid1 + 144);
    const auto *sid1_147 = buffer.data(sid1 + 147);
    const auto *sid1_149 = buffer.data(sid1 + 149);

    const auto *sif_166 = buffer.data(sif + 166);
    const auto *sif_167 = buffer.data(sif + 167);
    const auto *sif_168 = buffer.data(sif + 168);
    const auto *sif_169 = buffer.data(sif + 169);
    const auto *sif_170 = buffer.data(sif + 170);
    const auto *sif_172 = buffer.data(sif + 172);
    const auto *sif_176 = buffer.data(sif + 176);
    const auto *sif_177 = buffer.data(sif + 177);
    const auto *sif_178 = buffer.data(sif + 178);
    const auto *sif_179 = buffer.data(sif + 179);
    const auto *sif_180 = buffer.data(sif + 180);
    const auto *sif_182 = buffer.data(sif + 182);
    const auto *sif_186 = buffer.data(sif + 186);
    const auto *sif_187 = buffer.data(sif + 187);
    const auto *sif_188 = buffer.data(sif + 188);
    const auto *sif_189 = buffer.data(sif + 189);
    const auto *sif_190 = buffer.data(sif + 190);
    const auto *sif_192 = buffer.data(sif + 192);
    const auto *sif_196 = buffer.data(sif + 196);
    const auto *sif_197 = buffer.data(sif + 197);
    const auto *sif_198 = buffer.data(sif + 198);
    const auto *sif_199 = buffer.data(sif + 199);
    const auto *sif_200 = buffer.data(sif + 200);
    const auto *sif_202 = buffer.data(sif + 202);
    const auto *sif_206 = buffer.data(sif + 206);
    const auto *sif_207 = buffer.data(sif + 207);
    const auto *sif_208 = buffer.data(sif + 208);
    const auto *sif_209 = buffer.data(sif + 209);
    const auto *sif_210 = buffer.data(sif + 210);
    const auto *sif_212 = buffer.data(sif + 212);
    const auto *sif_213 = buffer.data(sif + 213);
    const auto *sif_215 = buffer.data(sif + 215);
    const auto *sif_216 = buffer.data(sif + 216);
    const auto *sif_217 = buffer.data(sif + 217);
    const auto *sif_218 = buffer.data(sif + 218);
    const auto *sif_219 = buffer.data(sif + 219);
    const auto *sif_220 = buffer.data(sif + 220);
    const auto *sif_222 = buffer.data(sif + 222);
    const auto *sif_225 = buffer.data(sif + 225);
    const auto *sif_226 = buffer.data(sif + 226);
    const auto *sif_227 = buffer.data(sif + 227);
    const auto *sif_228 = buffer.data(sif + 228);
    const auto *sif_229 = buffer.data(sif + 229);
    const auto *sif_230 = buffer.data(sif + 230);
    const auto *sif_232 = buffer.data(sif + 232);
    const auto *sif_233 = buffer.data(sif + 233);
    const auto *sif_235 = buffer.data(sif + 235);
    const auto *sif_236 = buffer.data(sif + 236);
    const auto *sif_237 = buffer.data(sif + 237);
    const auto *sif_238 = buffer.data(sif + 238);
    const auto *sif_239 = buffer.data(sif + 239);
    const auto *sif_240 = buffer.data(sif + 240);
    const auto *sif_242 = buffer.data(sif + 242);
    const auto *sif_243 = buffer.data(sif + 243);
    const auto *sif_245 = buffer.data(sif + 245);
    const auto *sif_246 = buffer.data(sif + 246);
    const auto *sif_247 = buffer.data(sif + 247);
    const auto *sif_248 = buffer.data(sif + 248);
    const auto *sif_249 = buffer.data(sif + 249);

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pb_x, pc_x, shg0_245, shf_165, shf_166, \
                         shf_167, shf_168, shg1_245, sif_166, sif_167, \
                         sif_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = pb_x[k] * shg0_245[k]
                   + f_8 * shf_165[k]
                   - f_6 * pc_x[k] * shg1_245[k];

        t_246[k] = f_7 * shf_166[k]
                   + f_3 * pc_x[k] * sif_166[k];

        t_247[k] = f_7 * shf_167[k]
                   + f_3 * pc_x[k] * sif_167[k];

        t_248[k] = f_7 * shf_168[k]
                   + f_3 * pc_x[k] * sif_168[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pb_x, pc_x, pc_z, shg0_250, shg0_252, \
                         shf_106, shf_169, shg1_250, shg1_252, sif_166, \
                         sif_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_7 * shf_169[k]
                   + f_3 * pc_x[k] * sif_169[k];

        t_250[k] = pb_x[k] * shg0_250[k]
                   - f_6 * pc_x[k] * shg1_250[k];

        t_251[k] = f_7 * shf_106[k]
                   + f_3 * pc_z[k] * sif_166[k];

        t_252[k] = pb_x[k] * shg0_252[k]
                   - f_6 * pc_x[k] * shg1_252[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pb_x, pc_x, pc_y, shg0_254, shg0_255, \
                         shf_119, shf_120, shf_170, shg1_254, shg1_255, sif_169, \
                         sif_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_10 * shf_119[k]
                   + f_3 * pc_y[k] * sif_169[k];

        t_254[k] = pb_x[k] * shg0_254[k]
                   - f_6 * pc_x[k] * shg1_254[k];

        t_255[k] = pb_x[k] * shg0_255[k]
                   + f_10 * shf_170[k]
                   - f_6 * pc_x[k] * shg1_255[k];

        t_256[k] = f_11 * shf_120[k]
                   + f_3 * pc_y[k] * sif_170[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_x, pc_x, pc_y, pc_z, shg0_258, shf_110, \
                         shf_122, shf_173, shg1_258, sif_170, sif_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_8 * shf_110[k]
                   + f_3 * pc_z[k] * sif_170[k];

        t_258[k] = pb_x[k] * shg0_258[k]
                   + f_8 * shf_173[k]
                   - f_6 * pc_x[k] * shg1_258[k];

        t_259[k] = f_11 * shf_122[k]
                   + f_3 * pc_y[k] * sif_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pb_x, pc_x, shg0_260, shf_175, shf_176, \
                         shf_177, shf_178, shg1_260, sif_176, sif_177, \
                         sif_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pb_x[k] * shg0_260[k]
                   + f_8 * shf_175[k]
                   - f_6 * pc_x[k] * shg1_260[k];

        t_261[k] = f_7 * shf_176[k]
                   + f_3 * pc_x[k] * sif_176[k];

        t_262[k] = f_7 * shf_177[k]
                   + f_3 * pc_x[k] * sif_177[k];

        t_263[k] = f_7 * shf_178[k]
                   + f_3 * pc_x[k] * sif_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pb_x, pc_x, pc_z, shg0_265, shg0_267, \
                         shf_116, shf_179, shg1_265, shg1_267, sif_176, \
                         sif_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_7 * shf_179[k]
                   + f_3 * pc_x[k] * sif_179[k];

        t_265[k] = pb_x[k] * shg0_265[k]
                   - f_6 * pc_x[k] * shg1_265[k];

        t_266[k] = f_8 * shf_116[k]
                   + f_3 * pc_z[k] * sif_176[k];

        t_267[k] = pb_x[k] * shg0_267[k]
                   - f_6 * pc_x[k] * shg1_267[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pb_x, pc_x, pc_y, shg0_269, shg0_270, \
                         shf_129, shf_130, shf_180, shg1_269, shg1_270, sif_179, \
                         sif_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_11 * shf_129[k]
                   + f_3 * pc_y[k] * sif_179[k];

        t_269[k] = pb_x[k] * shg0_269[k]
                   - f_6 * pc_x[k] * shg1_269[k];

        t_270[k] = pb_x[k] * shg0_270[k]
                   + f_10 * shf_180[k]
                   - f_6 * pc_x[k] * shg1_270[k];

        t_271[k] = f_8 * shf_130[k]
                   + f_3 * pc_y[k] * sif_180[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pb_x, pc_x, pc_y, pc_z, shg0_273, shf_120, \
                         shf_132, shf_183, shg1_273, sif_180, sif_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_11 * shf_120[k]
                   + f_3 * pc_z[k] * sif_180[k];

        t_273[k] = pb_x[k] * shg0_273[k]
                   + f_8 * shf_183[k]
                   - f_6 * pc_x[k] * shg1_273[k];

        t_274[k] = f_8 * shf_132[k]
                   + f_3 * pc_y[k] * sif_182[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_x, pc_x, shg0_275, shf_185, shf_186, \
                         shf_187, shf_188, shg1_275, sif_186, sif_187, \
                         sif_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = pb_x[k] * shg0_275[k]
                   + f_8 * shf_185[k]
                   - f_6 * pc_x[k] * shg1_275[k];

        t_276[k] = f_7 * shf_186[k]
                   + f_3 * pc_x[k] * sif_186[k];

        t_277[k] = f_7 * shf_187[k]
                   + f_3 * pc_x[k] * sif_187[k];

        t_278[k] = f_7 * shf_188[k]
                   + f_3 * pc_x[k] * sif_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_x, pc_x, pc_z, shg0_280, shg0_282, \
                         shf_126, shf_189, shg1_280, shg1_282, sif_186, \
                         sif_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_7 * shf_189[k]
                   + f_3 * pc_x[k] * sif_189[k];

        t_280[k] = pb_x[k] * shg0_280[k]
                   - f_6 * pc_x[k] * shg1_280[k];

        t_281[k] = f_11 * shf_126[k]
                   + f_3 * pc_z[k] * sif_186[k];

        t_282[k] = pb_x[k] * shg0_282[k]
                   - f_6 * pc_x[k] * shg1_282[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pb_x, pb_y, pc_x, pc_y, shg0_210, \
                         shg0_284, shf_139, shf_140, shg1_210, shg1_284, sif_189, \
                         sif_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * shf_139[k]
                   + f_3 * pc_y[k] * sif_189[k];

        t_284[k] = pb_x[k] * shg0_284[k]
                   - f_6 * pc_x[k] * shg1_284[k];

        t_285[k] = pb_y[k] * shg0_210[k]
                   - f_6 * pc_y[k] * shg1_210[k];

        t_286[k] = f_7 * shf_140[k]
                   + f_3 * pc_y[k] * sif_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pb_x, pc_x, pc_y, pc_z, shg0_288, shf_130, \
                         shf_142, shf_193, shg1_288, sif_190, sif_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_10 * shf_130[k]
                   + f_3 * pc_z[k] * sif_190[k];

        t_288[k] = pb_x[k] * shg0_288[k]
                   + f_8 * shf_193[k]
                   - f_6 * pc_x[k] * shg1_288[k];

        t_289[k] = f_7 * shf_142[k]
                   + f_3 * pc_y[k] * sif_192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_x, pc_y, shg0_215, shf_196, \
                         shf_197, shf_198, shg1_215, sif_196, sif_197, \
                         sif_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_y[k] * shg0_215[k]
                   - f_6 * pc_y[k] * shg1_215[k];

        t_291[k] = f_7 * shf_196[k]
                   + f_3 * pc_x[k] * sif_196[k];

        t_292[k] = f_7 * shf_197[k]
                   + f_3 * pc_x[k] * sif_197[k];

        t_293[k] = f_7 * shf_198[k]
                   + f_3 * pc_x[k] * sif_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_x, pc_x, pc_z, shg0_295, shg0_297, \
                         shf_136, shf_199, shg1_295, shg1_297, sif_196, \
                         sif_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_7 * shf_199[k]
                   + f_3 * pc_x[k] * sif_199[k];

        t_295[k] = pb_x[k] * shg0_295[k]
                   - f_6 * pc_x[k] * shg1_295[k];

        t_296[k] = f_10 * shf_136[k]
                   + f_3 * pc_z[k] * sif_196[k];

        t_297[k] = pb_x[k] * shg0_297[k]
                   - f_6 * pc_x[k] * shg1_297[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pb_x, pc_x, pc_y, shg0_299, shg0_300, \
                         shf_149, shf_200, shg1_299, shg1_300, sif_199, \
                         sif_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * shf_149[k]
                   + f_3 * pc_y[k] * sif_199[k];

        t_299[k] = pb_x[k] * shg0_299[k]
                   - f_6 * pc_x[k] * shg1_299[k];

        t_300[k] = pb_x[k] * shg0_300[k]
                   + f_10 * shf_200[k]
                   - f_6 * pc_x[k] * shg1_300[k];

        t_301[k] = f_3 * pc_y[k] * sif_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pb_x, pc_x, pc_y, pc_z, shg0_303, shf_140, \
                         shf_203, shg1_303, sif_200, sif_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_9 * shf_140[k]
                   + f_3 * pc_z[k] * sif_200[k];

        t_303[k] = pb_x[k] * shg0_303[k]
                   + f_8 * shf_203[k]
                   - f_6 * pc_x[k] * shg1_303[k];

        t_304[k] = f_3 * pc_y[k] * sif_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_x, pc_x, shg0_305, shf_205, shf_206, \
                         shf_207, shf_208, shg1_305, sif_206, sif_207, \
                         sif_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_x[k] * shg0_305[k]
                   + f_8 * shf_205[k]
                   - f_6 * pc_x[k] * shg1_305[k];

        t_306[k] = f_7 * shf_206[k]
                   + f_3 * pc_x[k] * sif_206[k];

        t_307[k] = f_7 * shf_207[k]
                   + f_3 * pc_x[k] * sif_207[k];

        t_308[k] = f_7 * shf_208[k]
                   + f_3 * pc_x[k] * sif_208[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_x, pc_x, pc_z, shg0_310, shg0_312, \
                         shf_146, shf_209, shg1_310, shg1_312, sif_206, \
                         sif_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_7 * shf_209[k]
                   + f_3 * pc_x[k] * sif_209[k];

        t_310[k] = pb_x[k] * shg0_310[k]
                   - f_6 * pc_x[k] * shg1_310[k];

        t_311[k] = f_9 * shf_146[k]
                   + f_3 * pc_z[k] * sif_206[k];

        t_312[k] = pb_x[k] * shg0_312[k]
                   - f_6 * pc_x[k] * shg1_312[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, pb_x, pc_x, pc_y, pc_z, shg0_314, \
                         shf_150, shg1_314, sid0_126, sid1_126, sif_209, \
                         sif_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * sif_209[k];

        t_314[k] = pb_x[k] * shg0_314[k]
                   - f_6 * pc_x[k] * shg1_314[k];

        t_315[k] = f_1 * sid0_126[k]
                   - f_2 * sid1_126[k]
                   + f_3 * pc_x[k] * sif_210[k];

        t_316[k] = f_0 * shf_150[k]
                   + f_3 * pc_y[k] * sif_210[k];

        t_317[k] = f_3 * pc_z[k] * sif_210[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pc_x, pc_y, shf_152, sid0_129, sid0_131, \
                         sid1_129, sid1_131, sif_212, sif_213, sif_215, \
                         sif_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_4 * sid0_129[k]
                   - f_5 * sid1_129[k]
                   + f_3 * pc_x[k] * sif_213[k];

        t_319[k] = f_0 * shf_152[k]
                   + f_3 * pc_y[k] * sif_212[k];

        t_320[k] = f_4 * sid0_131[k]
                   - f_5 * sid1_131[k]
                   + f_3 * pc_x[k] * sif_215[k];

        t_321[k] = f_3 * pc_x[k] * sif_216[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, pc_x, pc_y, pc_z, shf_156, \
                         sid0_129, sid1_129, sif_216, sif_217, sif_218, \
                         sif_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_3 * pc_x[k] * sif_217[k];

        t_323[k] = f_3 * pc_x[k] * sif_218[k];

        t_324[k] = f_3 * pc_x[k] * sif_219[k];

        t_325[k] = f_0 * shf_156[k]
                   + f_1 * sid0_129[k]
                   - f_2 * sid1_129[k]
                   + f_3 * pc_y[k] * sif_216[k];

        t_326[k] = f_3 * pc_z[k] * sif_216[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pb_z, pc_y, pc_z, shg0_225, shf_158, \
                         shf_159, shg1_225, sid0_131, sid1_131, sif_218, \
                         sif_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_0 * shf_158[k]
                   + f_4 * sid0_131[k]
                   - f_5 * sid1_131[k]
                   + f_3 * pc_y[k] * sif_218[k];

        t_328[k] = f_0 * shf_159[k]
                   + f_3 * pc_y[k] * sif_219[k];

        t_329[k] = f_1 * sid0_131[k]
                   - f_2 * sid1_131[k]
                   + f_3 * pc_z[k] * sif_219[k];

        t_330[k] = pb_z[k] * shg0_225[k]
                   - f_6 * pc_z[k] * shg1_225[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_z, pc_y, pc_z, shg0_228, shf_150, \
                         shf_160, shf_162, shg1_228, sif_220, sif_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_9 * shf_160[k]
                   + f_3 * pc_y[k] * sif_220[k];

        t_332[k] = f_7 * shf_150[k]
                   + f_3 * pc_z[k] * sif_220[k];

        t_333[k] = pb_z[k] * shg0_228[k]
                   - f_6 * pc_z[k] * shg1_228[k];

        t_334[k] = f_9 * shf_162[k]
                   + f_3 * pc_y[k] * sif_222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, pc_x, sid0_137, sid1_137, sif_225, \
                         sif_226, sif_227, sif_228, sif_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_4 * sid0_137[k]
                   - f_5 * sid1_137[k]
                   + f_3 * pc_x[k] * sif_225[k];

        t_336[k] = f_3 * pc_x[k] * sif_226[k];

        t_337[k] = f_3 * pc_x[k] * sif_227[k];

        t_338[k] = f_3 * pc_x[k] * sif_228[k];

        t_339[k] = f_3 * pc_x[k] * sif_229[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pb_z, pc_y, pc_z, shg0_235, shg0_237, \
                         shf_156, shf_157, shf_169, shg1_235, shg1_237, sif_226, \
                         sif_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pb_z[k] * shg0_235[k]
                   - f_6 * pc_z[k] * shg1_235[k];

        t_341[k] = f_7 * shf_156[k]
                   + f_3 * pc_z[k] * sif_226[k];

        t_342[k] = pb_z[k] * shg0_237[k]
                   + f_8 * shf_157[k]
                   - f_6 * pc_z[k] * shg1_237[k];

        t_343[k] = f_9 * shf_169[k]
                   + f_3 * pc_y[k] * sif_229[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, pc_z, shf_159, shf_160, \
                         shf_170, sid0_137, sid0_138, sid1_137, sid1_138, sif_229, \
                         sif_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_7 * shf_159[k]
                   + f_1 * sid0_137[k]
                   - f_2 * sid1_137[k]
                   + f_3 * pc_z[k] * sif_229[k];

        t_345[k] = f_1 * sid0_138[k]
                   - f_2 * sid1_138[k]
                   + f_3 * pc_x[k] * sif_230[k];

        t_346[k] = f_10 * shf_170[k]
                   + f_3 * pc_y[k] * sif_230[k];

        t_347[k] = f_8 * shf_160[k]
                   + f_3 * pc_z[k] * sif_230[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pc_x, pc_y, shf_172, sid0_141, sid0_143, \
                         sid1_141, sid1_143, sif_232, sif_233, sif_235, \
                         sif_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_4 * sid0_141[k]
                   - f_5 * sid1_141[k]
                   + f_3 * pc_x[k] * sif_233[k];

        t_349[k] = f_10 * shf_172[k]
                   + f_3 * pc_y[k] * sif_232[k];

        t_350[k] = f_4 * sid0_143[k]
                   - f_5 * sid1_143[k]
                   + f_3 * pc_x[k] * sif_235[k];

        t_351[k] = f_3 * pc_x[k] * sif_236[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, pc_x, pc_y, pc_z, shf_166, \
                         shf_176, sid0_141, sid1_141, sif_236, sif_237, sif_238, \
                         sif_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_3 * pc_x[k] * sif_237[k];

        t_353[k] = f_3 * pc_x[k] * sif_238[k];

        t_354[k] = f_3 * pc_x[k] * sif_239[k];

        t_355[k] = f_10 * shf_176[k]
                   + f_1 * sid0_141[k]
                   - f_2 * sid1_141[k]
                   + f_3 * pc_y[k] * sif_236[k];

        t_356[k] = f_8 * shf_166[k]
                   + f_3 * pc_z[k] * sif_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, shf_169, shf_178, shf_179, sid0_143, \
                         sid1_143, sif_238, sif_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_10 * shf_178[k]
                   + f_4 * sid0_143[k]
                   - f_5 * sid1_143[k]
                   + f_3 * pc_y[k] * sif_238[k];

        t_358[k] = f_10 * shf_179[k]
                   + f_3 * pc_y[k] * sif_239[k];

        t_359[k] = f_8 * shf_169[k]
                   + f_1 * sid0_143[k]
                   - f_2 * sid1_143[k]
                   + f_3 * pc_z[k] * sif_239[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, shf_170, shf_180, \
                         sid0_144, sid0_147, sid1_144, sid1_147, sif_240, \
                         sif_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_1 * sid0_144[k]
                   - f_2 * sid1_144[k]
                   + f_3 * pc_x[k] * sif_240[k];

        t_361[k] = f_11 * shf_180[k]
                   + f_3 * pc_y[k] * sif_240[k];

        t_362[k] = f_11 * shf_170[k]
                   + f_3 * pc_z[k] * sif_240[k];

        t_363[k] = f_4 * sid0_147[k]
                   - f_5 * sid1_147[k]
                   + f_3 * pc_x[k] * sif_243[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, pc_x, pc_y, shf_182, sid0_149, \
                         sid1_149, sif_242, sif_245, sif_246, sif_247, \
                         sif_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_11 * shf_182[k]
                   + f_3 * pc_y[k] * sif_242[k];

        t_365[k] = f_4 * sid0_149[k]
                   - f_5 * sid1_149[k]
                   + f_3 * pc_x[k] * sif_245[k];

        t_366[k] = f_3 * pc_x[k] * sif_246[k];

        t_367[k] = f_3 * pc_x[k] * sif_247[k];

        t_368[k] = f_3 * pc_x[k] * sif_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pc_x, pc_y, pc_z, shf_176, shf_186, sid0_147, \
                         sid1_147, sif_246, sif_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_3 * pc_x[k] * sif_249[k];

        t_370[k] = f_11 * shf_186[k]
                   + f_1 * sid0_147[k]
                   - f_2 * sid1_147[k]
                   + f_3 * pc_y[k] * sif_246[k];

        t_371[k] = f_11 * shf_176[k]
                   + f_3 * pc_z[k] * sif_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, shf_179, shf_188, shf_189, sid0_149, \
                         sid1_149, sif_248, sif_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_11 * shf_188[k]
                   + f_4 * sid0_149[k]
                   - f_5 * sid1_149[k]
                   + f_3 * pc_y[k] * sif_248[k];

        t_373[k] = f_11 * shf_189[k]
                   + f_3 * pc_y[k] * sif_249[k];

        t_374[k] = f_11 * shf_179[k]
                   + f_1 * sid0_149[k]
                   - f_2 * sid1_149[k]
                   + f_3 * pc_z[k] * sif_249[k];
    }
}

static auto
compute_prim_sig_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shg0,
                                                          const size_t shf, const size_t shg1,
                                                          const size_t sid0, const size_t sid1,
                                                          const size_t sif, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 2.0 / q;

    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shg0_300 = buffer.data(shg0 + 300);
    const auto *shg0_305 = buffer.data(shg0 + 305);
    const auto *shg0_310 = buffer.data(shg0 + 310);
    const auto *shg0_312 = buffer.data(shg0 + 312);
    const auto *shg0_314 = buffer.data(shg0 + 314);

    const auto *shf_180 = buffer.data(shf + 180);
    const auto *shf_186 = buffer.data(shf + 186);
    const auto *shf_189 = buffer.data(shf + 189);
    const auto *shf_190 = buffer.data(shf + 190);
    const auto *shf_192 = buffer.data(shf + 192);
    const auto *shf_196 = buffer.data(shf + 196);
    const auto *shf_198 = buffer.data(shf + 198);
    const auto *shf_199 = buffer.data(shf + 199);
    const auto *shf_200 = buffer.data(shf + 200);
    const auto *shf_202 = buffer.data(shf + 202);
    const auto *shf_206 = buffer.data(shf + 206);
    const auto *shf_208 = buffer.data(shf + 208);
    const auto *shf_209 = buffer.data(shf + 209);

    const auto *shg1_300 = buffer.data(shg1 + 300);
    const auto *shg1_305 = buffer.data(shg1 + 305);
    const auto *shg1_310 = buffer.data(shg1 + 310);
    const auto *shg1_312 = buffer.data(shg1 + 312);
    const auto *shg1_314 = buffer.data(shg1 + 314);

    const auto *sid0_150 = buffer.data(sid0 + 150);
    const auto *sid0_153 = buffer.data(sid0 + 153);
    const auto *sid0_155 = buffer.data(sid0 + 155);
    const auto *sid0_159 = buffer.data(sid0 + 159);
    const auto *sid0_162 = buffer.data(sid0 + 162);
    const auto *sid0_165 = buffer.data(sid0 + 165);
    const auto *sid0_167 = buffer.data(sid0 + 167);

    const auto *sid1_150 = buffer.data(sid1 + 150);
    const auto *sid1_153 = buffer.data(sid1 + 153);
    const auto *sid1_155 = buffer.data(sid1 + 155);
    const auto *sid1_159 = buffer.data(sid1 + 159);
    const auto *sid1_162 = buffer.data(sid1 + 162);
    const auto *sid1_165 = buffer.data(sid1 + 165);
    const auto *sid1_167 = buffer.data(sid1 + 167);

    const auto *sif_250 = buffer.data(sif + 250);
    const auto *sif_252 = buffer.data(sif + 252);
    const auto *sif_253 = buffer.data(sif + 253);
    const auto *sif_255 = buffer.data(sif + 255);
    const auto *sif_256 = buffer.data(sif + 256);
    const auto *sif_257 = buffer.data(sif + 257);
    const auto *sif_258 = buffer.data(sif + 258);
    const auto *sif_259 = buffer.data(sif + 259);
    const auto *sif_260 = buffer.data(sif + 260);
    const auto *sif_262 = buffer.data(sif + 262);
    const auto *sif_263 = buffer.data(sif + 263);
    const auto *sif_266 = buffer.data(sif + 266);
    const auto *sif_267 = buffer.data(sif + 267);
    const auto *sif_268 = buffer.data(sif + 268);
    const auto *sif_269 = buffer.data(sif + 269);
    const auto *sif_270 = buffer.data(sif + 270);
    const auto *sif_272 = buffer.data(sif + 272);
    const auto *sif_273 = buffer.data(sif + 273);
    const auto *sif_275 = buffer.data(sif + 275);
    const auto *sif_276 = buffer.data(sif + 276);
    const auto *sif_277 = buffer.data(sif + 277);
    const auto *sif_278 = buffer.data(sif + 278);
    const auto *sif_279 = buffer.data(sif + 279);

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pc_x, pc_y, pc_z, shf_180, shf_190, \
                         sid0_150, sid0_153, sid1_150, sid1_153, sif_250, \
                         sif_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_1 * sid0_150[k]
                   - f_2 * sid1_150[k]
                   + f_3 * pc_x[k] * sif_250[k];

        t_376[k] = f_8 * shf_190[k]
                   + f_3 * pc_y[k] * sif_250[k];

        t_377[k] = f_10 * shf_180[k]
                   + f_3 * pc_z[k] * sif_250[k];

        t_378[k] = f_4 * sid0_153[k]
                   - f_5 * sid1_153[k]
                   + f_3 * pc_x[k] * sif_253[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, pc_x, pc_y, shf_192, sid0_155, \
                         sid1_155, sif_252, sif_255, sif_256, sif_257, \
                         sif_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_8 * shf_192[k]
                   + f_3 * pc_y[k] * sif_252[k];

        t_380[k] = f_4 * sid0_155[k]
                   - f_5 * sid1_155[k]
                   + f_3 * pc_x[k] * sif_255[k];

        t_381[k] = f_3 * pc_x[k] * sif_256[k];

        t_382[k] = f_3 * pc_x[k] * sif_257[k];

        t_383[k] = f_3 * pc_x[k] * sif_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_y, pc_z, shf_186, shf_196, sid0_153, \
                         sid1_153, sif_256, sif_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_3 * pc_x[k] * sif_259[k];

        t_385[k] = f_8 * shf_196[k]
                   + f_1 * sid0_153[k]
                   - f_2 * sid1_153[k]
                   + f_3 * pc_y[k] * sif_256[k];

        t_386[k] = f_10 * shf_186[k]
                   + f_3 * pc_z[k] * sif_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pb_y, pc_y, pc_z, shg0_300, shf_189, \
                         shf_198, shf_199, shg1_300, sid0_155, sid1_155, sif_258, \
                         sif_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * shf_198[k]
                   + f_4 * sid0_155[k]
                   - f_5 * sid1_155[k]
                   + f_3 * pc_y[k] * sif_258[k];

        t_388[k] = f_8 * shf_199[k]
                   + f_3 * pc_y[k] * sif_259[k];

        t_389[k] = f_10 * shf_189[k]
                   + f_1 * sid0_155[k]
                   - f_2 * sid1_155[k]
                   + f_3 * pc_z[k] * sif_259[k];

        t_390[k] = pb_y[k] * shg0_300[k]
                   - f_6 * pc_y[k] * shg1_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, shf_190, shf_200, \
                         shf_202, sid0_159, sid1_159, sif_260, sif_262, \
                         sif_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_7 * shf_200[k]
                   + f_3 * pc_y[k] * sif_260[k];

        t_392[k] = f_9 * shf_190[k]
                   + f_3 * pc_z[k] * sif_260[k];

        t_393[k] = f_4 * sid0_159[k]
                   - f_5 * sid1_159[k]
                   + f_3 * pc_x[k] * sif_263[k];

        t_394[k] = f_7 * shf_202[k]
                   + f_3 * pc_y[k] * sif_262[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, pb_y, pc_x, pc_y, shg0_305, \
                         shg1_305, sif_266, sif_267, sif_268, sif_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pb_y[k] * shg0_305[k]
                   - f_6 * pc_y[k] * shg1_305[k];

        t_396[k] = f_3 * pc_x[k] * sif_266[k];

        t_397[k] = f_3 * pc_x[k] * sif_267[k];

        t_398[k] = f_3 * pc_x[k] * sif_268[k];

        t_399[k] = f_3 * pc_x[k] * sif_269[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pb_y, pc_y, pc_z, shg0_310, shg0_312, shf_196, \
                         shf_206, shf_208, shg1_310, shg1_312, \
                         sif_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = pb_y[k] * shg0_310[k]
                   + f_10 * shf_206[k]
                   - f_6 * pc_y[k] * shg1_310[k];

        t_401[k] = f_9 * shf_196[k]
                   + f_3 * pc_z[k] * sif_266[k];

        t_402[k] = pb_y[k] * shg0_312[k]
                   + f_8 * shf_208[k]
                   - f_6 * pc_y[k] * shg1_312[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pb_y, pc_x, pc_y, shg0_314, shf_209, \
                         shg1_314, sid0_162, sid1_162, sif_269, \
                         sif_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_7 * shf_209[k]
                   + f_3 * pc_y[k] * sif_269[k];

        t_404[k] = pb_y[k] * shg0_314[k]
                   - f_6 * pc_y[k] * shg1_314[k];

        t_405[k] = f_1 * sid0_162[k]
                   - f_2 * sid1_162[k]
                   + f_3 * pc_x[k] * sif_270[k];

        t_406[k] = f_3 * pc_y[k] * sif_270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, pc_x, pc_y, pc_z, shf_200, sid0_165, \
                         sid0_167, sid1_165, sid1_167, sif_270, sif_272, sif_273, \
                         sif_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_0 * shf_200[k]
                   + f_3 * pc_z[k] * sif_270[k];

        t_408[k] = f_4 * sid0_165[k]
                   - f_5 * sid1_165[k]
                   + f_3 * pc_x[k] * sif_273[k];

        t_409[k] = f_3 * pc_y[k] * sif_272[k];

        t_410[k] = f_4 * sid0_167[k]
                   - f_5 * sid1_167[k]
                   + f_3 * pc_x[k] * sif_275[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, t_416, pc_x, pc_y, pc_z, shf_206, \
                         sid0_165, sid1_165, sif_276, sif_277, sif_278, \
                         sif_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_3 * pc_x[k] * sif_276[k];

        t_412[k] = f_3 * pc_x[k] * sif_277[k];

        t_413[k] = f_3 * pc_x[k] * sif_278[k];

        t_414[k] = f_3 * pc_x[k] * sif_279[k];

        t_415[k] = f_1 * sid0_165[k]
                   - f_2 * sid1_165[k]
                   + f_3 * pc_y[k] * sif_276[k];

        t_416[k] = f_0 * shf_206[k]
                   + f_3 * pc_z[k] * sif_276[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_y, pc_z, shf_209, sid0_167, sid1_167, \
                         sif_278, sif_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * sid0_167[k]
                   - f_5 * sid1_167[k]
                   + f_3 * pc_y[k] * sif_278[k];

        t_418[k] = f_3 * pc_y[k] * sif_279[k];

        t_419[k] = f_0 * shf_209[k]
                   + f_1 * sid0_167[k]
                   - f_2 * sid1_167[k]
                   + f_3 * pc_z[k] * sif_279[k];
    }
}

auto
compute_prim_sig_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t shg0, const size_t shf,
                                                   const size_t shg1, const size_t sid0,
                                                   const size_t sid1, const size_t sif,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sig_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, shg0, shf,
                                                              shg1, sid0, sid1, sif, ncols,
                                                              gamma, p, q);

    compute_prim_sig_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, shg0, shf,
                                                              shg1, sid0, sid1, sif, ncols,
                                                              gamma, p, q);

    compute_prim_sig_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, shg0, shf,
                                                              shg1, sid0, sid1, sif, ncols,
                                                              gamma, p, q);

    compute_prim_sig_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, shg0, shf,
                                                              shg1, sid0, sid1, sif, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
