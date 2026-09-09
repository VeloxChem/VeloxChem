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


#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_0 = buffer.data(sli0 + 0);
    const auto *sli0_3 = buffer.data(sli0 + 3);
    const auto *sli0_5 = buffer.data(sli0 + 5);
    const auto *sli0_6 = buffer.data(sli0 + 6);
    const auto *sli0_9 = buffer.data(sli0 + 9);
    const auto *sli0_10 = buffer.data(sli0 + 10);
    const auto *sli0_12 = buffer.data(sli0 + 12);
    const auto *sli0_14 = buffer.data(sli0 + 14);
    const auto *sli0_21 = buffer.data(sli0 + 21);
    const auto *sli0_27 = buffer.data(sli0 + 27);
    const auto *sli0_31 = buffer.data(sli0 + 31);
    const auto *sli0_34 = buffer.data(sli0 + 34);
    const auto *sli0_38 = buffer.data(sli0 + 38);
    const auto *sli0_56 = buffer.data(sli0 + 56);
    const auto *sli0_61 = buffer.data(sli0 + 61);
    const auto *sli0_65 = buffer.data(sli0 + 65);

    const auto *slh_0 = buffer.data(slh + 0);
    const auto *slh_1 = buffer.data(slh + 1);
    const auto *slh_2 = buffer.data(slh + 2);
    const auto *slh_3 = buffer.data(slh + 3);
    const auto *slh_5 = buffer.data(slh + 5);
    const auto *slh_6 = buffer.data(slh + 6);
    const auto *slh_7 = buffer.data(slh + 7);
    const auto *slh_8 = buffer.data(slh + 8);
    const auto *slh_9 = buffer.data(slh + 9);
    const auto *slh_10 = buffer.data(slh + 10);
    const auto *slh_12 = buffer.data(slh + 12);
    const auto *slh_14 = buffer.data(slh + 14);
    const auto *slh_15 = buffer.data(slh + 15);
    const auto *slh_16 = buffer.data(slh + 16);
    const auto *slh_17 = buffer.data(slh + 17);
    const auto *slh_18 = buffer.data(slh + 18);
    const auto *slh_19 = buffer.data(slh + 19);
    const auto *slh_20 = buffer.data(slh + 20);
    const auto *slh_21 = buffer.data(slh + 21);
    const auto *slh_23 = buffer.data(slh + 23);
    const auto *slh_24 = buffer.data(slh + 24);
    const auto *slh_26 = buffer.data(slh + 26);
    const auto *slh_27 = buffer.data(slh + 27);
    const auto *slh_30 = buffer.data(slh + 30);
    const auto *slh_36 = buffer.data(slh + 36);
    const auto *slh_37 = buffer.data(slh + 37);
    const auto *slh_38 = buffer.data(slh + 38);
    const auto *slh_39 = buffer.data(slh + 39);
    const auto *slh_40 = buffer.data(slh + 40);
    const auto *slh_41 = buffer.data(slh + 41);
    const auto *slh_42 = buffer.data(slh + 42);
    const auto *slh_44 = buffer.data(slh + 44);
    const auto *slh_47 = buffer.data(slh + 47);
    const auto *slh_57 = buffer.data(slh + 57);
    const auto *slh_58 = buffer.data(slh + 58);
    const auto *slh_59 = buffer.data(slh + 59);
    const auto *slh_60 = buffer.data(slh + 60);
    const auto *slh_61 = buffer.data(slh + 61);
    const auto *slh_62 = buffer.data(slh + 62);
    const auto *slh_63 = buffer.data(slh + 63);
    const auto *slh_66 = buffer.data(slh + 66);
    const auto *slh_68 = buffer.data(slh + 68);
    const auto *slh_69 = buffer.data(slh + 69);
    const auto *slh_72 = buffer.data(slh + 72);
    const auto *slh_73 = buffer.data(slh + 73);
    const auto *slh_75 = buffer.data(slh + 75);
    const auto *slh_77 = buffer.data(slh + 77);
    const auto *slh_78 = buffer.data(slh + 78);
    const auto *slh_79 = buffer.data(slh + 79);
    const auto *slh_80 = buffer.data(slh + 80);
    const auto *slh_81 = buffer.data(slh + 81);
    const auto *slh_82 = buffer.data(slh + 82);
    const auto *slh_83 = buffer.data(slh + 83);

    const auto *sli1_0 = buffer.data(sli1 + 0);
    const auto *sli1_3 = buffer.data(sli1 + 3);
    const auto *sli1_5 = buffer.data(sli1 + 5);
    const auto *sli1_6 = buffer.data(sli1 + 6);
    const auto *sli1_9 = buffer.data(sli1 + 9);
    const auto *sli1_10 = buffer.data(sli1 + 10);
    const auto *sli1_12 = buffer.data(sli1 + 12);
    const auto *sli1_14 = buffer.data(sli1 + 14);
    const auto *sli1_21 = buffer.data(sli1 + 21);
    const auto *sli1_27 = buffer.data(sli1 + 27);
    const auto *sli1_31 = buffer.data(sli1 + 31);
    const auto *sli1_34 = buffer.data(sli1 + 34);
    const auto *sli1_38 = buffer.data(sli1 + 38);
    const auto *sli1_56 = buffer.data(sli1 + 56);
    const auto *sli1_61 = buffer.data(sli1 + 61);
    const auto *sli1_65 = buffer.data(sli1 + 65);

    const auto *smg0_0 = buffer.data(smg0 + 0);
    const auto *smg0_3 = buffer.data(smg0 + 3);
    const auto *smg0_5 = buffer.data(smg0 + 5);
    const auto *smg0_6 = buffer.data(smg0 + 6);
    const auto *smg0_9 = buffer.data(smg0 + 9);
    const auto *smg0_10 = buffer.data(smg0 + 10);
    const auto *smg0_12 = buffer.data(smg0 + 12);
    const auto *smg0_13 = buffer.data(smg0 + 13);
    const auto *smg0_14 = buffer.data(smg0 + 14);
    const auto *smg0_25 = buffer.data(smg0 + 25);
    const auto *smg0_27 = buffer.data(smg0 + 27);
    const auto *smg0_28 = buffer.data(smg0 + 28);
    const auto *smg0_29 = buffer.data(smg0 + 29);
    const auto *smg0_42 = buffer.data(smg0 + 42);
    const auto *smg0_43 = buffer.data(smg0 + 43);
    const auto *smg0_44 = buffer.data(smg0 + 44);
    const auto *smg0_45 = buffer.data(smg0 + 45);
    const auto *smg0_48 = buffer.data(smg0 + 48);
    const auto *smg0_50 = buffer.data(smg0 + 50);
    const auto *smg0_51 = buffer.data(smg0 + 51);
    const auto *smg0_54 = buffer.data(smg0 + 54);
    const auto *smg0_55 = buffer.data(smg0 + 55);
    const auto *smg0_57 = buffer.data(smg0 + 57);
    const auto *smg0_58 = buffer.data(smg0 + 58);
    const auto *smg0_59 = buffer.data(smg0 + 59);

    const auto *smg1_0 = buffer.data(smg1 + 0);
    const auto *smg1_3 = buffer.data(smg1 + 3);
    const auto *smg1_5 = buffer.data(smg1 + 5);
    const auto *smg1_6 = buffer.data(smg1 + 6);
    const auto *smg1_9 = buffer.data(smg1 + 9);
    const auto *smg1_10 = buffer.data(smg1 + 10);
    const auto *smg1_12 = buffer.data(smg1 + 12);
    const auto *smg1_13 = buffer.data(smg1 + 13);
    const auto *smg1_14 = buffer.data(smg1 + 14);
    const auto *smg1_25 = buffer.data(smg1 + 25);
    const auto *smg1_27 = buffer.data(smg1 + 27);
    const auto *smg1_28 = buffer.data(smg1 + 28);
    const auto *smg1_29 = buffer.data(smg1 + 29);
    const auto *smg1_42 = buffer.data(smg1 + 42);
    const auto *smg1_43 = buffer.data(smg1 + 43);
    const auto *smg1_44 = buffer.data(smg1 + 44);
    const auto *smg1_45 = buffer.data(smg1 + 45);
    const auto *smg1_48 = buffer.data(smg1 + 48);
    const auto *smg1_50 = buffer.data(smg1 + 50);
    const auto *smg1_51 = buffer.data(smg1 + 51);
    const auto *smg1_54 = buffer.data(smg1 + 54);
    const auto *smg1_55 = buffer.data(smg1 + 55);
    const auto *smg1_57 = buffer.data(smg1 + 57);
    const auto *smg1_58 = buffer.data(smg1 + 58);
    const auto *smg1_59 = buffer.data(smg1 + 59);

    const auto *smh_0 = buffer.data(smh + 0);
    const auto *smh_2 = buffer.data(smh + 2);
    const auto *smh_3 = buffer.data(smh + 3);
    const auto *smh_5 = buffer.data(smh + 5);
    const auto *smh_6 = buffer.data(smh + 6);
    const auto *smh_9 = buffer.data(smh + 9);
    const auto *smh_10 = buffer.data(smh + 10);
    const auto *smh_12 = buffer.data(smh + 12);
    const auto *smh_14 = buffer.data(smh + 14);
    const auto *smh_15 = buffer.data(smh + 15);
    const auto *smh_16 = buffer.data(smh + 16);
    const auto *smh_17 = buffer.data(smh + 17);
    const auto *smh_18 = buffer.data(smh + 18);
    const auto *smh_19 = buffer.data(smh + 19);
    const auto *smh_20 = buffer.data(smh + 20);
    const auto *smh_21 = buffer.data(smh + 21);
    const auto *smh_23 = buffer.data(smh + 23);
    const auto *smh_24 = buffer.data(smh + 24);
    const auto *smh_26 = buffer.data(smh + 26);
    const auto *smh_27 = buffer.data(smh + 27);
    const auto *smh_30 = buffer.data(smh + 30);
    const auto *smh_36 = buffer.data(smh + 36);
    const auto *smh_37 = buffer.data(smh + 37);
    const auto *smh_38 = buffer.data(smh + 38);
    const auto *smh_39 = buffer.data(smh + 39);
    const auto *smh_40 = buffer.data(smh + 40);
    const auto *smh_41 = buffer.data(smh + 41);
    const auto *smh_42 = buffer.data(smh + 42);
    const auto *smh_44 = buffer.data(smh + 44);
    const auto *smh_45 = buffer.data(smh + 45);
    const auto *smh_47 = buffer.data(smh + 47);
    const auto *smh_48 = buffer.data(smh + 48);
    const auto *smh_51 = buffer.data(smh + 51);
    const auto *smh_57 = buffer.data(smh + 57);
    const auto *smh_58 = buffer.data(smh + 58);
    const auto *smh_59 = buffer.data(smh + 59);
    const auto *smh_60 = buffer.data(smh + 60);
    const auto *smh_61 = buffer.data(smh + 61);
    const auto *smh_62 = buffer.data(smh + 62);
    const auto *smh_63 = buffer.data(smh + 63);
    const auto *smh_65 = buffer.data(smh + 65);
    const auto *smh_66 = buffer.data(smh + 66);
    const auto *smh_68 = buffer.data(smh + 68);
    const auto *smh_69 = buffer.data(smh + 69);
    const auto *smh_72 = buffer.data(smh + 72);
    const auto *smh_73 = buffer.data(smh + 73);
    const auto *smh_75 = buffer.data(smh + 75);
    const auto *smh_77 = buffer.data(smh + 77);
    const auto *smh_78 = buffer.data(smh + 78);
    const auto *smh_79 = buffer.data(smh + 79);
    const auto *smh_80 = buffer.data(smh + 80);
    const auto *smh_81 = buffer.data(smh + 81);
    const auto *smh_82 = buffer.data(smh + 82);
    const auto *smh_83 = buffer.data(smh + 83);
    const auto *smh_84 = buffer.data(smh + 84);
    const auto *smh_86 = buffer.data(smh + 86);
    const auto *smh_87 = buffer.data(smh + 87);
    const auto *smh_89 = buffer.data(smh + 89);
    const auto *smh_90 = buffer.data(smh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, slh_0, slh_3, smg0_0, smg0_3, \
                         smg1_0, smg1_3, smh_0, smh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * slh_0[k]
                 + f_1 * smg0_0[k]
                 - f_2 * smg1_0[k]
                 + f_3 * pc_x[k] * smh_0[k];

        t_1[k] = f_3 * pc_y[k] * smh_0[k];

        t_2[k] = f_3 * pc_z[k] * smh_0[k];

        t_3[k] = f_0 * slh_3[k]
                 + f_4 * smg0_3[k]
                 - f_5 * smg1_3[k]
                 + f_3 * pc_x[k] * smh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, slh_5, slh_6, smg0_5, smg0_6, smg1_5, \
                         smg1_6, smh_2, smh_5, smh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * smh_2[k];

        t_5[k] = f_0 * slh_5[k]
                 + f_4 * smg0_5[k]
                 - f_5 * smg1_5[k]
                 + f_3 * pc_x[k] * smh_5[k];

        t_6[k] = f_0 * slh_6[k]
                 + f_6 * smg0_6[k]
                 - f_7 * smg1_6[k]
                 + f_3 * pc_x[k] * smh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, slh_9, smg0_9, smg1_9, smh_3, smh_5, \
                         smh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * smh_3[k];

        t_8[k] = f_3 * pc_y[k] * smh_5[k];

        t_9[k] = f_0 * slh_9[k]
                 + f_6 * smg0_9[k]
                 - f_7 * smg1_9[k]
                 + f_3 * pc_x[k] * smh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, slh_10, slh_12, smg0_10, smg0_12, \
                         smg1_10, smg1_12, smh_6, smh_10, smh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * slh_10[k]
                  + f_8 * smg0_10[k]
                  - f_9 * smg1_10[k]
                  + f_3 * pc_x[k] * smh_10[k];

        t_11[k] = f_3 * pc_z[k] * smh_6[k];

        t_12[k] = f_0 * slh_12[k]
                  + f_8 * smg0_12[k]
                  - f_9 * smg1_12[k]
                  + f_3 * pc_x[k] * smh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, slh_14, slh_15, slh_16, smg0_14, \
                         smg1_14, smh_9, smh_14, smh_15, smh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * smh_9[k];

        t_14[k] = f_0 * slh_14[k]
                  + f_8 * smg0_14[k]
                  - f_9 * smg1_14[k]
                  + f_3 * pc_x[k] * smh_14[k];

        t_15[k] = f_0 * slh_15[k]
                  + f_3 * pc_x[k] * smh_15[k];

        t_16[k] = f_0 * slh_16[k]
                  + f_3 * pc_x[k] * smh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, slh_17, slh_18, slh_19, slh_20, smh_17, \
                         smh_18, smh_19, smh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * slh_17[k]
                  + f_3 * pc_x[k] * smh_17[k];

        t_18[k] = f_0 * slh_18[k]
                  + f_3 * pc_x[k] * smh_18[k];

        t_19[k] = f_0 * slh_19[k]
                  + f_3 * pc_x[k] * smh_19[k];

        t_20[k] = f_0 * slh_20[k]
                  + f_3 * pc_x[k] * smh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, smg0_10, smg0_12, smg0_13, \
                         smg1_10, smg1_12, smg1_13, smh_15, smh_17, \
                         smh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * smg0_10[k]
                  - f_2 * smg1_10[k]
                  + f_3 * pc_y[k] * smh_15[k];

        t_22[k] = f_3 * pc_z[k] * smh_15[k];

        t_23[k] = f_4 * smg0_12[k]
                  - f_5 * smg1_12[k]
                  + f_3 * pc_y[k] * smh_17[k];

        t_24[k] = f_6 * smg0_13[k]
                  - f_7 * smg1_13[k]
                  + f_3 * pc_y[k] * smh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, sli0_0, slh_0, \
                         sli1_0, smg0_14, smg1_14, smh_19, smh_20, \
                         smh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * smg0_14[k]
                  - f_9 * smg1_14[k]
                  + f_3 * pc_y[k] * smh_19[k];

        t_26[k] = f_3 * pc_y[k] * smh_20[k];

        t_27[k] = f_1 * smg0_14[k]
                  - f_2 * smg1_14[k]
                  + f_3 * pc_z[k] * smh_20[k];

        t_28[k] = pb_y[k] * sli0_0[k]
                  - f_10 * pc_y[k] * sli1_0[k];

        t_29[k] = f_11 * slh_0[k]
                  + f_3 * pc_y[k] * smh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, sli0_3, sli0_5, slh_1, \
                         slh_2, sli1_3, sli1_5, smh_21, smh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * smh_21[k];

        t_31[k] = pb_y[k] * sli0_3[k]
                  + f_12 * slh_1[k]
                  - f_10 * pc_y[k] * sli1_3[k];

        t_32[k] = f_11 * slh_2[k]
                  + f_3 * pc_y[k] * smh_23[k];

        t_33[k] = pb_y[k] * sli0_5[k]
                  - f_10 * pc_y[k] * sli1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, sli0_6, sli0_9, slh_3, \
                         slh_5, sli1_6, sli1_9, smh_24, smh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * sli0_6[k]
                  + f_13 * slh_3[k]
                  - f_10 * pc_y[k] * sli1_6[k];

        t_35[k] = f_3 * pc_z[k] * smh_24[k];

        t_36[k] = f_11 * slh_5[k]
                  + f_3 * pc_y[k] * smh_26[k];

        t_37[k] = pb_y[k] * sli0_9[k]
                  - f_10 * pc_y[k] * sli1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sli0_10, sli0_12, slh_6, \
                         slh_8, slh_9, sli1_10, sli1_12, smh_27, \
                         smh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * sli0_10[k]
                  + f_14 * slh_6[k]
                  - f_10 * pc_y[k] * sli1_10[k];

        t_39[k] = f_3 * pc_z[k] * smh_27[k];

        t_40[k] = pb_y[k] * sli0_12[k]
                  + f_12 * slh_8[k]
                  - f_10 * pc_y[k] * sli1_12[k];

        t_41[k] = f_11 * slh_9[k]
                  + f_3 * pc_y[k] * smh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, sli0_14, slh_36, slh_37, \
                         slh_38, sli1_14, smh_36, smh_37, smh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * sli0_14[k]
                  - f_10 * pc_y[k] * sli1_14[k];

        t_43[k] = f_15 * slh_36[k]
                  + f_3 * pc_x[k] * smh_36[k];

        t_44[k] = f_15 * slh_37[k]
                  + f_3 * pc_x[k] * smh_37[k];

        t_45[k] = f_15 * slh_38[k]
                  + f_3 * pc_x[k] * smh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, slh_15, slh_39, slh_40, slh_41, \
                         smg0_25, smg1_25, smh_36, smh_39, smh_40, \
                         smh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * slh_39[k]
                  + f_3 * pc_x[k] * smh_39[k];

        t_47[k] = f_15 * slh_40[k]
                  + f_3 * pc_x[k] * smh_40[k];

        t_48[k] = f_15 * slh_41[k]
                  + f_3 * pc_x[k] * smh_41[k];

        t_49[k] = f_11 * slh_15[k]
                  + f_1 * smg0_25[k]
                  - f_2 * smg1_25[k]
                  + f_3 * pc_y[k] * smh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, slh_17, slh_18, smg0_27, smg0_28, \
                         smg1_27, smg1_28, smh_36, smh_38, smh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * smh_36[k];

        t_51[k] = f_11 * slh_17[k]
                  + f_4 * smg0_27[k]
                  - f_5 * smg1_27[k]
                  + f_3 * pc_y[k] * smh_38[k];

        t_52[k] = f_11 * slh_18[k]
                  + f_6 * smg0_28[k]
                  - f_7 * smg1_28[k]
                  + f_3 * pc_y[k] * smh_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, sli0_27, slh_19, slh_20, sli1_27, \
                         smg0_29, smg1_29, smh_40, smh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * slh_19[k]
                  + f_8 * smg0_29[k]
                  - f_9 * smg1_29[k]
                  + f_3 * pc_y[k] * smh_40[k];

        t_54[k] = f_11 * slh_20[k]
                  + f_3 * pc_y[k] * smh_41[k];

        t_55[k] = pb_y[k] * sli0_27[k]
                  - f_10 * pc_y[k] * sli1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, sli0_0, sli0_3, \
                         slh_0, sli1_0, sli1_3, smh_42, smh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * sli0_0[k]
                  - f_10 * pc_z[k] * sli1_0[k];

        t_57[k] = f_3 * pc_y[k] * smh_42[k];

        t_58[k] = f_11 * slh_0[k]
                  + f_3 * pc_z[k] * smh_42[k];

        t_59[k] = pb_z[k] * sli0_3[k]
                  - f_10 * pc_z[k] * sli1_3[k];

        t_60[k] = f_3 * pc_y[k] * smh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sli0_5, sli0_6, slh_2, \
                         slh_3, sli1_5, sli1_6, smh_45, smh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * sli0_5[k]
                  + f_12 * slh_2[k]
                  - f_10 * pc_z[k] * sli1_5[k];

        t_62[k] = pb_z[k] * sli0_6[k]
                  - f_10 * pc_z[k] * sli1_6[k];

        t_63[k] = f_11 * slh_3[k]
                  + f_3 * pc_z[k] * smh_45[k];

        t_64[k] = f_3 * pc_y[k] * smh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, sli0_9, sli0_10, sli0_12, slh_5, \
                         slh_6, slh_7, sli1_9, sli1_10, sli1_12, \
                         smh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * sli0_9[k]
                  + f_13 * slh_5[k]
                  - f_10 * pc_z[k] * sli1_9[k];

        t_66[k] = pb_z[k] * sli0_10[k]
                  - f_10 * pc_z[k] * sli1_10[k];

        t_67[k] = f_11 * slh_6[k]
                  + f_3 * pc_z[k] * smh_48[k];

        t_68[k] = pb_z[k] * sli0_12[k]
                  + f_12 * slh_7[k]
                  - f_10 * pc_z[k] * sli1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, sli0_14, slh_9, \
                         slh_57, slh_58, sli1_14, smh_51, smh_57, \
                         smh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * smh_51[k];

        t_70[k] = pb_z[k] * sli0_14[k]
                  + f_14 * slh_9[k]
                  - f_10 * pc_z[k] * sli1_14[k];

        t_71[k] = f_15 * slh_57[k]
                  + f_3 * pc_x[k] * smh_57[k];

        t_72[k] = f_15 * slh_58[k]
                  + f_3 * pc_x[k] * smh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, slh_59, slh_60, slh_61, slh_62, smh_59, \
                         smh_60, smh_61, smh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * slh_59[k]
                  + f_3 * pc_x[k] * smh_59[k];

        t_74[k] = f_15 * slh_60[k]
                  + f_3 * pc_x[k] * smh_60[k];

        t_75[k] = f_15 * slh_61[k]
                  + f_3 * pc_x[k] * smh_61[k];

        t_76[k] = f_15 * slh_62[k]
                  + f_3 * pc_x[k] * smh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, sli0_21, slh_15, sli1_21, \
                         smg0_42, smg1_42, smh_57, smh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * sli0_21[k]
                  - f_10 * pc_z[k] * sli1_21[k];

        t_78[k] = f_11 * slh_15[k]
                  + f_3 * pc_z[k] * smh_57[k];

        t_79[k] = f_4 * smg0_42[k]
                  - f_5 * smg1_42[k]
                  + f_3 * pc_y[k] * smh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, slh_20, smg0_43, smg0_44, \
                         smg1_43, smg1_44, smh_60, smh_61, smh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * smg0_43[k]
                  - f_7 * smg1_43[k]
                  + f_3 * pc_y[k] * smh_60[k];

        t_81[k] = f_8 * smg0_44[k]
                  - f_9 * smg1_44[k]
                  + f_3 * pc_y[k] * smh_61[k];

        t_82[k] = f_3 * pc_y[k] * smh_62[k];

        t_83[k] = f_11 * slh_20[k]
                  + f_1 * smg0_44[k]
                  - f_2 * smg1_44[k]
                  + f_3 * pc_z[k] * smh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, slh_21, slh_63, slh_66, \
                         smg0_45, smg0_48, smg1_45, smg1_48, smh_63, \
                         smh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_16 * slh_63[k]
                  + f_1 * smg0_45[k]
                  - f_2 * smg1_45[k]
                  + f_3 * pc_x[k] * smh_63[k];

        t_85[k] = f_12 * slh_21[k]
                  + f_3 * pc_y[k] * smh_63[k];

        t_86[k] = f_3 * pc_z[k] * smh_63[k];

        t_87[k] = f_16 * slh_66[k]
                  + f_4 * smg0_48[k]
                  - f_5 * smg1_48[k]
                  + f_3 * pc_x[k] * smh_66[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pc_x, pc_y, slh_23, slh_68, slh_69, smg0_50, \
                         smg0_51, smg1_50, smg1_51, smh_65, smh_68, \
                         smh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * slh_23[k]
                  + f_3 * pc_y[k] * smh_65[k];

        t_89[k] = f_16 * slh_68[k]
                  + f_4 * smg0_50[k]
                  - f_5 * smg1_50[k]
                  + f_3 * pc_x[k] * smh_68[k];

        t_90[k] = f_16 * slh_69[k]
                  + f_6 * smg0_51[k]
                  - f_7 * smg1_51[k]
                  + f_3 * pc_x[k] * smh_69[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pc_x, pc_y, pc_z, slh_26, slh_72, smg0_54, smg1_54, \
                         smh_66, smh_68, smh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * smh_66[k];

        t_92[k] = f_12 * slh_26[k]
                  + f_3 * pc_y[k] * smh_68[k];

        t_93[k] = f_16 * slh_72[k]
                  + f_6 * smg0_54[k]
                  - f_7 * smg1_54[k]
                  + f_3 * pc_x[k] * smh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pc_x, pc_z, slh_73, slh_75, smg0_55, smg0_57, \
                         smg1_55, smg1_57, smh_69, smh_73, smh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_16 * slh_73[k]
                  + f_8 * smg0_55[k]
                  - f_9 * smg1_55[k]
                  + f_3 * pc_x[k] * smh_73[k];

        t_95[k] = f_3 * pc_z[k] * smh_69[k];

        t_96[k] = f_16 * slh_75[k]
                  + f_8 * smg0_57[k]
                  - f_9 * smg1_57[k]
                  + f_3 * pc_x[k] * smh_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, slh_30, slh_77, slh_78, slh_79, \
                         smg0_59, smg1_59, smh_72, smh_77, smh_78, \
                         smh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * slh_30[k]
                  + f_3 * pc_y[k] * smh_72[k];

        t_98[k] = f_16 * slh_77[k]
                  + f_8 * smg0_59[k]
                  - f_9 * smg1_59[k]
                  + f_3 * pc_x[k] * smh_77[k];

        t_99[k] = f_16 * slh_78[k]
                  + f_3 * pc_x[k] * smh_78[k];

        t_100[k] = f_16 * slh_79[k]
                   + f_3 * pc_x[k] * smh_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, slh_80, slh_81, slh_82, slh_83, \
                         smh_80, smh_81, smh_82, smh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_16 * slh_80[k]
                   + f_3 * pc_x[k] * smh_80[k];

        t_102[k] = f_16 * slh_81[k]
                   + f_3 * pc_x[k] * smh_81[k];

        t_103[k] = f_16 * slh_82[k]
                   + f_3 * pc_x[k] * smh_82[k];

        t_104[k] = f_16 * slh_83[k]
                   + f_3 * pc_x[k] * smh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, slh_36, slh_38, smg0_55, smg0_57, \
                         smg1_55, smg1_57, smh_78, smh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * slh_36[k]
                   + f_1 * smg0_55[k]
                   - f_2 * smg1_55[k]
                   + f_3 * pc_y[k] * smh_78[k];

        t_106[k] = f_3 * pc_z[k] * smh_78[k];

        t_107[k] = f_12 * slh_38[k]
                   + f_4 * smg0_57[k]
                   - f_5 * smg1_57[k]
                   + f_3 * pc_y[k] * smh_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, slh_39, slh_40, slh_41, \
                         smg0_58, smg0_59, smg1_58, smg1_59, smh_81, smh_82, \
                         smh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_12 * slh_39[k]
                   + f_6 * smg0_58[k]
                   - f_7 * smg1_58[k]
                   + f_3 * pc_y[k] * smh_81[k];

        t_109[k] = f_12 * slh_40[k]
                   + f_8 * smg0_59[k]
                   - f_9 * smg1_59[k]
                   + f_3 * pc_y[k] * smh_82[k];

        t_110[k] = f_12 * slh_41[k]
                   + f_3 * pc_y[k] * smh_83[k];

        t_111[k] = f_1 * smg0_59[k]
                   - f_2 * smg1_59[k]
                   + f_3 * pc_z[k] * smh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, sli0_31, sli0_56, \
                         slh_21, slh_42, sli1_31, sli1_56, smh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * sli0_56[k]
                   - f_10 * pc_y[k] * sli1_56[k];

        t_113[k] = f_11 * slh_42[k]
                   + f_3 * pc_y[k] * smh_84[k];

        t_114[k] = f_11 * slh_21[k]
                   + f_3 * pc_z[k] * smh_84[k];

        t_115[k] = pb_z[k] * sli0_31[k]
                   - f_10 * pc_z[k] * sli1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, sli0_34, sli0_61, \
                         slh_24, slh_44, sli1_34, sli1_61, smh_86, \
                         smh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * slh_44[k]
                   + f_3 * pc_y[k] * smh_86[k];

        t_117[k] = pb_y[k] * sli0_61[k]
                   - f_10 * pc_y[k] * sli1_61[k];

        t_118[k] = pb_z[k] * sli0_34[k]
                   - f_10 * pc_z[k] * sli1_34[k];

        t_119[k] = f_11 * slh_24[k]
                   + f_3 * pc_z[k] * smh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, sli0_38, sli0_65, \
                         slh_27, slh_47, sli1_38, sli1_65, smh_89, \
                         smh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * slh_47[k]
                   + f_3 * pc_y[k] * smh_89[k];

        t_121[k] = pb_y[k] * sli0_65[k]
                   - f_10 * pc_y[k] * sli1_65[k];

        t_122[k] = pb_z[k] * sli0_38[k]
                   - f_10 * pc_z[k] * sli1_38[k];

        t_123[k] = f_11 * slh_27[k]
                   + f_3 * pc_z[k] * smh_90[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;

    auto *t_124 = buffer.data(target + 124);
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_49 = buffer.data(sli0 + 49);
    const auto *sli0_68 = buffer.data(sli0 + 68);
    const auto *sli0_70 = buffer.data(sli0 + 70);
    const auto *sli0_83 = buffer.data(sli0 + 83);
    const auto *sli0_84 = buffer.data(sli0 + 84);
    const auto *sli0_87 = buffer.data(sli0 + 87);
    const auto *sli0_90 = buffer.data(sli0 + 90);
    const auto *sli0_94 = buffer.data(sli0 + 94);
    const auto *sli0_96 = buffer.data(sli0 + 96);
    const auto *sli0_105 = buffer.data(sli0 + 105);
    const auto *sli0_140 = buffer.data(sli0 + 140);
    const auto *sli0_143 = buffer.data(sli0 + 143);
    const auto *sli0_145 = buffer.data(sli0 + 145);
    const auto *sli0_146 = buffer.data(sli0 + 146);
    const auto *sli0_149 = buffer.data(sli0 + 149);
    const auto *sli0_150 = buffer.data(sli0 + 150);
    const auto *sli0_152 = buffer.data(sli0 + 152);
    const auto *sli0_154 = buffer.data(sli0 + 154);

    const auto *slh_36 = buffer.data(slh + 36);
    const auto *slh_42 = buffer.data(slh + 42);
    const auto *slh_45 = buffer.data(slh + 45);
    const auto *slh_48 = buffer.data(slh + 48);
    const auto *slh_50 = buffer.data(slh + 50);
    const auto *slh_51 = buffer.data(slh + 51);
    const auto *slh_57 = buffer.data(slh + 57);
    const auto *slh_59 = buffer.data(slh + 59);
    const auto *slh_60 = buffer.data(slh + 60);
    const auto *slh_61 = buffer.data(slh + 61);
    const auto *slh_62 = buffer.data(slh + 62);
    const auto *slh_63 = buffer.data(slh + 63);
    const auto *slh_65 = buffer.data(slh + 65);
    const auto *slh_66 = buffer.data(slh + 66);
    const auto *slh_68 = buffer.data(slh + 68);
    const auto *slh_69 = buffer.data(slh + 69);
    const auto *slh_70 = buffer.data(slh + 70);
    const auto *slh_72 = buffer.data(slh + 72);
    const auto *slh_78 = buffer.data(slh + 78);
    const auto *slh_80 = buffer.data(slh + 80);
    const auto *slh_81 = buffer.data(slh + 81);
    const auto *slh_82 = buffer.data(slh + 82);
    const auto *slh_83 = buffer.data(slh + 83);
    const auto *slh_84 = buffer.data(slh + 84);
    const auto *slh_86 = buffer.data(slh + 86);
    const auto *slh_87 = buffer.data(slh + 87);
    const auto *slh_89 = buffer.data(slh + 89);
    const auto *slh_90 = buffer.data(slh + 90);
    const auto *slh_93 = buffer.data(slh + 93);
    const auto *slh_99 = buffer.data(slh + 99);
    const auto *slh_100 = buffer.data(slh + 100);
    const auto *slh_101 = buffer.data(slh + 101);
    const auto *slh_102 = buffer.data(slh + 102);
    const auto *slh_103 = buffer.data(slh + 103);
    const auto *slh_104 = buffer.data(slh + 104);
    const auto *slh_105 = buffer.data(slh + 105);
    const auto *slh_106 = buffer.data(slh + 106);
    const auto *slh_107 = buffer.data(slh + 107);
    const auto *slh_108 = buffer.data(slh + 108);
    const auto *slh_110 = buffer.data(slh + 110);
    const auto *slh_111 = buffer.data(slh + 111);
    const auto *slh_113 = buffer.data(slh + 113);
    const auto *slh_114 = buffer.data(slh + 114);
    const auto *slh_115 = buffer.data(slh + 115);
    const auto *slh_117 = buffer.data(slh + 117);
    const auto *slh_119 = buffer.data(slh + 119);
    const auto *slh_120 = buffer.data(slh + 120);
    const auto *slh_121 = buffer.data(slh + 121);
    const auto *slh_122 = buffer.data(slh + 122);
    const auto *slh_123 = buffer.data(slh + 123);
    const auto *slh_124 = buffer.data(slh + 124);
    const auto *slh_125 = buffer.data(slh + 125);
    const auto *slh_126 = buffer.data(slh + 126);
    const auto *slh_129 = buffer.data(slh + 129);
    const auto *slh_131 = buffer.data(slh + 131);
    const auto *slh_132 = buffer.data(slh + 132);
    const auto *slh_135 = buffer.data(slh + 135);
    const auto *slh_136 = buffer.data(slh + 136);
    const auto *slh_138 = buffer.data(slh + 138);
    const auto *slh_140 = buffer.data(slh + 140);
    const auto *slh_141 = buffer.data(slh + 141);
    const auto *slh_142 = buffer.data(slh + 142);
    const auto *slh_143 = buffer.data(slh + 143);
    const auto *slh_144 = buffer.data(slh + 144);
    const auto *slh_145 = buffer.data(slh + 145);
    const auto *slh_146 = buffer.data(slh + 146);
    const auto *slh_152 = buffer.data(slh + 152);
    const auto *slh_156 = buffer.data(slh + 156);
    const auto *slh_161 = buffer.data(slh + 161);
    const auto *slh_162 = buffer.data(slh + 162);
    const auto *slh_163 = buffer.data(slh + 163);
    const auto *slh_164 = buffer.data(slh + 164);
    const auto *slh_165 = buffer.data(slh + 165);
    const auto *slh_166 = buffer.data(slh + 166);
    const auto *slh_167 = buffer.data(slh + 167);
    const auto *slh_183 = buffer.data(slh + 183);

    const auto *sli1_49 = buffer.data(sli1 + 49);
    const auto *sli1_68 = buffer.data(sli1 + 68);
    const auto *sli1_70 = buffer.data(sli1 + 70);
    const auto *sli1_83 = buffer.data(sli1 + 83);
    const auto *sli1_84 = buffer.data(sli1 + 84);
    const auto *sli1_87 = buffer.data(sli1 + 87);
    const auto *sli1_90 = buffer.data(sli1 + 90);
    const auto *sli1_94 = buffer.data(sli1 + 94);
    const auto *sli1_96 = buffer.data(sli1 + 96);
    const auto *sli1_105 = buffer.data(sli1 + 105);
    const auto *sli1_140 = buffer.data(sli1 + 140);
    const auto *sli1_143 = buffer.data(sli1 + 143);
    const auto *sli1_145 = buffer.data(sli1 + 145);
    const auto *sli1_146 = buffer.data(sli1 + 146);
    const auto *sli1_149 = buffer.data(sli1 + 149);
    const auto *sli1_150 = buffer.data(sli1 + 150);
    const auto *sli1_152 = buffer.data(sli1 + 152);
    const auto *sli1_154 = buffer.data(sli1 + 154);

    const auto *smg0_72 = buffer.data(smg0 + 72);
    const auto *smg0_73 = buffer.data(smg0 + 73);
    const auto *smg0_74 = buffer.data(smg0 + 74);
    const auto *smg0_75 = buffer.data(smg0 + 75);
    const auto *smg0_78 = buffer.data(smg0 + 78);
    const auto *smg0_80 = buffer.data(smg0 + 80);
    const auto *smg0_81 = buffer.data(smg0 + 81);
    const auto *smg0_84 = buffer.data(smg0 + 84);
    const auto *smg0_85 = buffer.data(smg0 + 85);
    const auto *smg0_87 = buffer.data(smg0 + 87);
    const auto *smg0_88 = buffer.data(smg0 + 88);
    const auto *smg0_89 = buffer.data(smg0 + 89);
    const auto *smg0_90 = buffer.data(smg0 + 90);
    const auto *smg0_93 = buffer.data(smg0 + 93);
    const auto *smg0_95 = buffer.data(smg0 + 95);
    const auto *smg0_96 = buffer.data(smg0 + 96);
    const auto *smg0_99 = buffer.data(smg0 + 99);
    const auto *smg0_100 = buffer.data(smg0 + 100);
    const auto *smg0_102 = buffer.data(smg0 + 102);
    const auto *smg0_103 = buffer.data(smg0 + 103);
    const auto *smg0_104 = buffer.data(smg0 + 104);
    const auto *smg0_110 = buffer.data(smg0 + 110);
    const auto *smg0_114 = buffer.data(smg0 + 114);
    const auto *smg0_117 = buffer.data(smg0 + 117);
    const auto *smg0_118 = buffer.data(smg0 + 118);
    const auto *smg0_119 = buffer.data(smg0 + 119);

    const auto *smg1_72 = buffer.data(smg1 + 72);
    const auto *smg1_73 = buffer.data(smg1 + 73);
    const auto *smg1_74 = buffer.data(smg1 + 74);
    const auto *smg1_75 = buffer.data(smg1 + 75);
    const auto *smg1_78 = buffer.data(smg1 + 78);
    const auto *smg1_80 = buffer.data(smg1 + 80);
    const auto *smg1_81 = buffer.data(smg1 + 81);
    const auto *smg1_84 = buffer.data(smg1 + 84);
    const auto *smg1_85 = buffer.data(smg1 + 85);
    const auto *smg1_87 = buffer.data(smg1 + 87);
    const auto *smg1_88 = buffer.data(smg1 + 88);
    const auto *smg1_89 = buffer.data(smg1 + 89);
    const auto *smg1_90 = buffer.data(smg1 + 90);
    const auto *smg1_93 = buffer.data(smg1 + 93);
    const auto *smg1_95 = buffer.data(smg1 + 95);
    const auto *smg1_96 = buffer.data(smg1 + 96);
    const auto *smg1_99 = buffer.data(smg1 + 99);
    const auto *smg1_100 = buffer.data(smg1 + 100);
    const auto *smg1_102 = buffer.data(smg1 + 102);
    const auto *smg1_103 = buffer.data(smg1 + 103);
    const auto *smg1_104 = buffer.data(smg1 + 104);
    const auto *smg1_110 = buffer.data(smg1 + 110);
    const auto *smg1_114 = buffer.data(smg1 + 114);
    const auto *smg1_117 = buffer.data(smg1 + 117);
    const auto *smg1_118 = buffer.data(smg1 + 118);
    const auto *smg1_119 = buffer.data(smg1 + 119);

    const auto *smh_93 = buffer.data(smh + 93);
    const auto *smh_99 = buffer.data(smh + 99);
    const auto *smh_100 = buffer.data(smh + 100);
    const auto *smh_101 = buffer.data(smh + 101);
    const auto *smh_102 = buffer.data(smh + 102);
    const auto *smh_103 = buffer.data(smh + 103);
    const auto *smh_104 = buffer.data(smh + 104);
    const auto *smh_105 = buffer.data(smh + 105);
    const auto *smh_107 = buffer.data(smh + 107);
    const auto *smh_108 = buffer.data(smh + 108);
    const auto *smh_110 = buffer.data(smh + 110);
    const auto *smh_111 = buffer.data(smh + 111);
    const auto *smh_114 = buffer.data(smh + 114);
    const auto *smh_115 = buffer.data(smh + 115);
    const auto *smh_117 = buffer.data(smh + 117);
    const auto *smh_119 = buffer.data(smh + 119);
    const auto *smh_120 = buffer.data(smh + 120);
    const auto *smh_121 = buffer.data(smh + 121);
    const auto *smh_122 = buffer.data(smh + 122);
    const auto *smh_123 = buffer.data(smh + 123);
    const auto *smh_124 = buffer.data(smh + 124);
    const auto *smh_125 = buffer.data(smh + 125);
    const auto *smh_126 = buffer.data(smh + 126);
    const auto *smh_128 = buffer.data(smh + 128);
    const auto *smh_129 = buffer.data(smh + 129);
    const auto *smh_131 = buffer.data(smh + 131);
    const auto *smh_132 = buffer.data(smh + 132);
    const auto *smh_135 = buffer.data(smh + 135);
    const auto *smh_136 = buffer.data(smh + 136);
    const auto *smh_138 = buffer.data(smh + 138);
    const auto *smh_140 = buffer.data(smh + 140);
    const auto *smh_141 = buffer.data(smh + 141);
    const auto *smh_142 = buffer.data(smh + 142);
    const auto *smh_143 = buffer.data(smh + 143);
    const auto *smh_144 = buffer.data(smh + 144);
    const auto *smh_145 = buffer.data(smh + 145);
    const auto *smh_146 = buffer.data(smh + 146);
    const auto *smh_147 = buffer.data(smh + 147);
    const auto *smh_149 = buffer.data(smh + 149);
    const auto *smh_150 = buffer.data(smh + 150);
    const auto *smh_152 = buffer.data(smh + 152);
    const auto *smh_153 = buffer.data(smh + 153);
    const auto *smh_156 = buffer.data(smh + 156);
    const auto *smh_161 = buffer.data(smh + 161);
    const auto *smh_162 = buffer.data(smh + 162);
    const auto *smh_163 = buffer.data(smh + 163);
    const auto *smh_164 = buffer.data(smh + 164);
    const auto *smh_165 = buffer.data(smh + 165);
    const auto *smh_166 = buffer.data(smh + 166);
    const auto *smh_167 = buffer.data(smh + 167);
    const auto *smh_168 = buffer.data(smh + 168);
    const auto *smh_170 = buffer.data(smh + 170);
    const auto *smh_171 = buffer.data(smh + 171);
    const auto *smh_173 = buffer.data(smh + 173);
    const auto *smh_174 = buffer.data(smh + 174);
    const auto *smh_177 = buffer.data(smh + 177);
    const auto *smh_183 = buffer.data(smh + 183);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, sli0_68, sli0_70, \
                         slh_50, slh_51, slh_99, sli1_68, sli1_70, smh_93, \
                         smh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * sli0_68[k]
                   + f_12 * slh_50[k]
                   - f_10 * pc_y[k] * sli1_68[k];

        t_125[k] = f_11 * slh_51[k]
                   + f_3 * pc_y[k] * smh_93[k];

        t_126[k] = pb_y[k] * sli0_70[k]
                   - f_10 * pc_y[k] * sli1_70[k];

        t_127[k] = f_16 * slh_99[k]
                   + f_3 * pc_x[k] * smh_99[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, slh_100, slh_101, slh_102, \
                         slh_103, slh_104, smh_100, smh_101, smh_102, smh_103, \
                         smh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_16 * slh_100[k]
                   + f_3 * pc_x[k] * smh_100[k];

        t_129[k] = f_16 * slh_101[k]
                   + f_3 * pc_x[k] * smh_101[k];

        t_130[k] = f_16 * slh_102[k]
                   + f_3 * pc_x[k] * smh_102[k];

        t_131[k] = f_16 * slh_103[k]
                   + f_3 * pc_x[k] * smh_103[k];

        t_132[k] = f_16 * slh_104[k]
                   + f_3 * pc_x[k] * smh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, pc_y, pc_z, sli0_49, slh_36, slh_59, \
                         sli1_49, smg0_72, smg1_72, smh_99, smh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_z[k] * sli0_49[k]
                   - f_10 * pc_z[k] * sli1_49[k];

        t_134[k] = f_11 * slh_36[k]
                   + f_3 * pc_z[k] * smh_99[k];

        t_135[k] = f_11 * slh_59[k]
                   + f_4 * smg0_72[k]
                   - f_5 * smg1_72[k]
                   + f_3 * pc_y[k] * smh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, slh_60, slh_61, slh_62, smg0_73, smg0_74, \
                         smg1_73, smg1_74, smh_102, smh_103, smh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * slh_60[k]
                   + f_6 * smg0_73[k]
                   - f_7 * smg1_73[k]
                   + f_3 * pc_y[k] * smh_102[k];

        t_137[k] = f_11 * slh_61[k]
                   + f_8 * smg0_74[k]
                   - f_9 * smg1_74[k]
                   + f_3 * pc_y[k] * smh_103[k];

        t_138[k] = f_11 * slh_62[k]
                   + f_3 * pc_y[k] * smh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_y, pc_x, pc_y, pc_z, sli0_83, slh_42, \
                         slh_105, sli1_83, smg0_75, smg1_75, smh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * sli0_83[k]
                   - f_10 * pc_y[k] * sli1_83[k];

        t_140[k] = f_16 * slh_105[k]
                   + f_1 * smg0_75[k]
                   - f_2 * smg1_75[k]
                   + f_3 * pc_x[k] * smh_105[k];

        t_141[k] = f_3 * pc_y[k] * smh_105[k];

        t_142[k] = f_12 * slh_42[k]
                   + f_3 * pc_z[k] * smh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, slh_108, slh_110, smg0_78, smg0_80, \
                         smg1_78, smg1_80, smh_107, smh_108, smh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * slh_108[k]
                   + f_4 * smg0_78[k]
                   - f_5 * smg1_78[k]
                   + f_3 * pc_x[k] * smh_108[k];

        t_144[k] = f_3 * pc_y[k] * smh_107[k];

        t_145[k] = f_16 * slh_110[k]
                   + f_4 * smg0_80[k]
                   - f_5 * smg1_80[k]
                   + f_3 * pc_x[k] * smh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_x, pc_y, pc_z, slh_45, slh_111, smg0_81, \
                         smg1_81, smh_108, smh_110, smh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_16 * slh_111[k]
                   + f_6 * smg0_81[k]
                   - f_7 * smg1_81[k]
                   + f_3 * pc_x[k] * smh_111[k];

        t_147[k] = f_12 * slh_45[k]
                   + f_3 * pc_z[k] * smh_108[k];

        t_148[k] = f_3 * pc_y[k] * smh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, slh_48, slh_114, slh_115, smg0_84, \
                         smg0_85, smg1_84, smg1_85, smh_111, smh_114, \
                         smh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_16 * slh_114[k]
                   + f_6 * smg0_84[k]
                   - f_7 * smg1_84[k]
                   + f_3 * pc_x[k] * smh_114[k];

        t_150[k] = f_16 * slh_115[k]
                   + f_8 * smg0_85[k]
                   - f_9 * smg1_85[k]
                   + f_3 * pc_x[k] * smh_115[k];

        t_151[k] = f_12 * slh_48[k]
                   + f_3 * pc_z[k] * smh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pc_x, pc_y, slh_117, slh_119, smg0_87, smg0_89, \
                         smg1_87, smg1_89, smh_114, smh_117, smh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_16 * slh_117[k]
                   + f_8 * smg0_87[k]
                   - f_9 * smg1_87[k]
                   + f_3 * pc_x[k] * smh_117[k];

        t_153[k] = f_3 * pc_y[k] * smh_114[k];

        t_154[k] = f_16 * slh_119[k]
                   + f_8 * smg0_89[k]
                   - f_9 * smg1_89[k]
                   + f_3 * pc_x[k] * smh_119[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, slh_120, slh_121, slh_122, \
                         slh_123, slh_124, smh_120, smh_121, smh_122, smh_123, \
                         smh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_16 * slh_120[k]
                   + f_3 * pc_x[k] * smh_120[k];

        t_156[k] = f_16 * slh_121[k]
                   + f_3 * pc_x[k] * smh_121[k];

        t_157[k] = f_16 * slh_122[k]
                   + f_3 * pc_x[k] * smh_122[k];

        t_158[k] = f_16 * slh_123[k]
                   + f_3 * pc_x[k] * smh_123[k];

        t_159[k] = f_16 * slh_124[k]
                   + f_3 * pc_x[k] * smh_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, slh_57, slh_125, \
                         smg0_85, smg0_87, smg1_85, smg1_87, smh_120, smh_122, \
                         smh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * slh_125[k]
                   + f_3 * pc_x[k] * smh_125[k];

        t_161[k] = f_1 * smg0_85[k]
                   - f_2 * smg1_85[k]
                   + f_3 * pc_y[k] * smh_120[k];

        t_162[k] = f_12 * slh_57[k]
                   + f_3 * pc_z[k] * smh_120[k];

        t_163[k] = f_4 * smg0_87[k]
                   - f_5 * smg1_87[k]
                   + f_3 * pc_y[k] * smh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, slh_62, smg0_88, smg0_89, \
                         smg1_88, smg1_89, smh_123, smh_124, smh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * smg0_88[k]
                   - f_7 * smg1_88[k]
                   + f_3 * pc_y[k] * smh_123[k];

        t_165[k] = f_8 * smg0_89[k]
                   - f_9 * smg1_89[k]
                   + f_3 * pc_y[k] * smh_124[k];

        t_166[k] = f_3 * pc_y[k] * smh_125[k];

        t_167[k] = f_12 * slh_62[k]
                   + f_1 * smg0_89[k]
                   - f_2 * smg1_89[k]
                   + f_3 * pc_z[k] * smh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, slh_63, slh_126, \
                         slh_129, smg0_90, smg0_93, smg1_90, smg1_93, smh_126, \
                         smh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_17 * slh_126[k]
                   + f_1 * smg0_90[k]
                   - f_2 * smg1_90[k]
                   + f_3 * pc_x[k] * smh_126[k];

        t_169[k] = f_13 * slh_63[k]
                   + f_3 * pc_y[k] * smh_126[k];

        t_170[k] = f_3 * pc_z[k] * smh_126[k];

        t_171[k] = f_17 * slh_129[k]
                   + f_4 * smg0_93[k]
                   - f_5 * smg1_93[k]
                   + f_3 * pc_x[k] * smh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pc_x, pc_y, slh_65, slh_131, slh_132, smg0_95, \
                         smg0_96, smg1_95, smg1_96, smh_128, smh_131, \
                         smh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * slh_65[k]
                   + f_3 * pc_y[k] * smh_128[k];

        t_173[k] = f_17 * slh_131[k]
                   + f_4 * smg0_95[k]
                   - f_5 * smg1_95[k]
                   + f_3 * pc_x[k] * smh_131[k];

        t_174[k] = f_17 * slh_132[k]
                   + f_6 * smg0_96[k]
                   - f_7 * smg1_96[k]
                   + f_3 * pc_x[k] * smh_132[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pc_x, pc_y, pc_z, slh_68, slh_135, smg0_99, \
                         smg1_99, smh_129, smh_131, smh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * pc_z[k] * smh_129[k];

        t_176[k] = f_13 * slh_68[k]
                   + f_3 * pc_y[k] * smh_131[k];

        t_177[k] = f_17 * slh_135[k]
                   + f_6 * smg0_99[k]
                   - f_7 * smg1_99[k]
                   + f_3 * pc_x[k] * smh_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_z, slh_136, slh_138, smg0_100, \
                         smg0_102, smg1_100, smg1_102, smh_132, smh_136, \
                         smh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_17 * slh_136[k]
                   + f_8 * smg0_100[k]
                   - f_9 * smg1_100[k]
                   + f_3 * pc_x[k] * smh_136[k];

        t_179[k] = f_3 * pc_z[k] * smh_132[k];

        t_180[k] = f_17 * slh_138[k]
                   + f_8 * smg0_102[k]
                   - f_9 * smg1_102[k]
                   + f_3 * pc_x[k] * smh_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, slh_72, slh_140, slh_141, \
                         slh_142, smg0_104, smg1_104, smh_135, smh_140, smh_141, \
                         smh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * slh_72[k]
                   + f_3 * pc_y[k] * smh_135[k];

        t_182[k] = f_17 * slh_140[k]
                   + f_8 * smg0_104[k]
                   - f_9 * smg1_104[k]
                   + f_3 * pc_x[k] * smh_140[k];

        t_183[k] = f_17 * slh_141[k]
                   + f_3 * pc_x[k] * smh_141[k];

        t_184[k] = f_17 * slh_142[k]
                   + f_3 * pc_x[k] * smh_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, slh_143, slh_144, slh_145, slh_146, \
                         smh_143, smh_144, smh_145, smh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_17 * slh_143[k]
                   + f_3 * pc_x[k] * smh_143[k];

        t_186[k] = f_17 * slh_144[k]
                   + f_3 * pc_x[k] * smh_144[k];

        t_187[k] = f_17 * slh_145[k]
                   + f_3 * pc_x[k] * smh_145[k];

        t_188[k] = f_17 * slh_146[k]
                   + f_3 * pc_x[k] * smh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_y, pc_z, slh_78, slh_80, smg0_100, smg0_102, \
                         smg1_100, smg1_102, smh_141, smh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * slh_78[k]
                   + f_1 * smg0_100[k]
                   - f_2 * smg1_100[k]
                   + f_3 * pc_y[k] * smh_141[k];

        t_190[k] = f_3 * pc_z[k] * smh_141[k];

        t_191[k] = f_13 * slh_80[k]
                   + f_4 * smg0_102[k]
                   - f_5 * smg1_102[k]
                   + f_3 * pc_y[k] * smh_143[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, slh_81, slh_82, slh_83, \
                         smg0_103, smg0_104, smg1_103, smg1_104, smh_144, smh_145, \
                         smh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_13 * slh_81[k]
                   + f_6 * smg0_103[k]
                   - f_7 * smg1_103[k]
                   + f_3 * pc_y[k] * smh_144[k];

        t_193[k] = f_13 * slh_82[k]
                   + f_8 * smg0_104[k]
                   - f_9 * smg1_104[k]
                   + f_3 * pc_y[k] * smh_145[k];

        t_194[k] = f_13 * slh_83[k]
                   + f_3 * pc_y[k] * smh_146[k];

        t_195[k] = f_1 * smg0_104[k]
                   - f_2 * smg1_104[k]
                   + f_3 * pc_z[k] * smh_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_z, pc_y, pc_z, sli0_84, sli0_87, \
                         slh_63, slh_84, sli1_84, sli1_87, smh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * sli0_84[k]
                   - f_10 * pc_z[k] * sli1_84[k];

        t_197[k] = f_12 * slh_84[k]
                   + f_3 * pc_y[k] * smh_147[k];

        t_198[k] = f_11 * slh_63[k]
                   + f_3 * pc_z[k] * smh_147[k];

        t_199[k] = pb_z[k] * sli0_87[k]
                   - f_10 * pc_z[k] * sli1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_z, pc_x, pc_y, pc_z, sli0_90, slh_86, \
                         slh_152, sli1_90, smg0_110, smg1_110, smh_149, \
                         smh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * slh_86[k]
                   + f_3 * pc_y[k] * smh_149[k];

        t_201[k] = f_17 * slh_152[k]
                   + f_4 * smg0_110[k]
                   - f_5 * smg1_110[k]
                   + f_3 * pc_x[k] * smh_152[k];

        t_202[k] = pb_z[k] * sli0_90[k]
                   - f_10 * pc_z[k] * sli1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, pc_y, pc_z, slh_66, slh_89, slh_156, \
                         smg0_114, smg1_114, smh_150, smh_152, \
                         smh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * slh_66[k]
                   + f_3 * pc_z[k] * smh_150[k];

        t_204[k] = f_12 * slh_89[k]
                   + f_3 * pc_y[k] * smh_152[k];

        t_205[k] = f_17 * slh_156[k]
                   + f_6 * smg0_114[k]
                   - f_7 * smg1_114[k]
                   + f_3 * pc_x[k] * smh_156[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_z, pc_y, pc_z, sli0_94, sli0_96, \
                         slh_69, slh_70, slh_93, sli1_94, sli1_96, smh_153, \
                         smh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * sli0_94[k]
                   - f_10 * pc_z[k] * sli1_94[k];

        t_207[k] = f_11 * slh_69[k]
                   + f_3 * pc_z[k] * smh_153[k];

        t_208[k] = pb_z[k] * sli0_96[k]
                   + f_12 * slh_70[k]
                   - f_10 * pc_z[k] * sli1_96[k];

        t_209[k] = f_12 * slh_93[k]
                   + f_3 * pc_y[k] * smh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, slh_161, slh_162, slh_163, slh_164, \
                         smg0_119, smg1_119, smh_161, smh_162, smh_163, \
                         smh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * slh_161[k]
                   + f_8 * smg0_119[k]
                   - f_9 * smg1_119[k]
                   + f_3 * pc_x[k] * smh_161[k];

        t_211[k] = f_17 * slh_162[k]
                   + f_3 * pc_x[k] * smh_162[k];

        t_212[k] = f_17 * slh_163[k]
                   + f_3 * pc_x[k] * smh_163[k];

        t_213[k] = f_17 * slh_164[k]
                   + f_3 * pc_x[k] * smh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, sli0_105, slh_165, \
                         slh_166, slh_167, sli1_105, smh_165, smh_166, \
                         smh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_17 * slh_165[k]
                   + f_3 * pc_x[k] * smh_165[k];

        t_215[k] = f_17 * slh_166[k]
                   + f_3 * pc_x[k] * smh_166[k];

        t_216[k] = f_17 * slh_167[k]
                   + f_3 * pc_x[k] * smh_167[k];

        t_217[k] = pb_z[k] * sli0_105[k]
                   - f_10 * pc_z[k] * sli1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, slh_78, slh_101, slh_102, smg0_117, \
                         smg0_118, smg1_117, smg1_118, smh_162, smh_164, \
                         smh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * slh_78[k]
                   + f_3 * pc_z[k] * smh_162[k];

        t_219[k] = f_12 * slh_101[k]
                   + f_4 * smg0_117[k]
                   - f_5 * smg1_117[k]
                   + f_3 * pc_y[k] * smh_164[k];

        t_220[k] = f_12 * slh_102[k]
                   + f_6 * smg0_118[k]
                   - f_7 * smg1_118[k]
                   + f_3 * pc_y[k] * smh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_y, pc_y, pc_z, sli0_140, slh_83, \
                         slh_103, slh_104, sli1_140, smg0_119, smg1_119, smh_166, \
                         smh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * slh_103[k]
                   + f_8 * smg0_119[k]
                   - f_9 * smg1_119[k]
                   + f_3 * pc_y[k] * smh_166[k];

        t_222[k] = f_12 * slh_104[k]
                   + f_3 * pc_y[k] * smh_167[k];

        t_223[k] = f_11 * slh_83[k]
                   + f_1 * smg0_119[k]
                   - f_2 * smg1_119[k]
                   + f_3 * pc_z[k] * smh_167[k];

        t_224[k] = pb_y[k] * sli0_140[k]
                   - f_10 * pc_y[k] * sli1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_y, pc_y, pc_z, sli0_143, slh_84, \
                         slh_105, slh_106, slh_107, sli1_143, smh_168, \
                         smh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * slh_105[k]
                   + f_3 * pc_y[k] * smh_168[k];

        t_226[k] = f_12 * slh_84[k]
                   + f_3 * pc_z[k] * smh_168[k];

        t_227[k] = pb_y[k] * sli0_143[k]
                   + f_12 * slh_106[k]
                   - f_10 * pc_y[k] * sli1_143[k];

        t_228[k] = f_11 * slh_107[k]
                   + f_3 * pc_y[k] * smh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_y, pc_y, pc_z, sli0_145, sli0_146, \
                         slh_87, slh_108, slh_110, sli1_145, sli1_146, smh_171, \
                         smh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * sli0_145[k]
                   - f_10 * pc_y[k] * sli1_145[k];

        t_230[k] = pb_y[k] * sli0_146[k]
                   + f_13 * slh_108[k]
                   - f_10 * pc_y[k] * sli1_146[k];

        t_231[k] = f_12 * slh_87[k]
                   + f_3 * pc_z[k] * smh_171[k];

        t_232[k] = f_11 * slh_110[k]
                   + f_3 * pc_y[k] * smh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pc_y, pc_z, sli0_149, sli0_150, slh_90, \
                         slh_111, sli1_149, sli1_150, smh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * sli0_149[k]
                   - f_10 * pc_y[k] * sli1_149[k];

        t_234[k] = pb_y[k] * sli0_150[k]
                   + f_14 * slh_111[k]
                   - f_10 * pc_y[k] * sli1_150[k];

        t_235[k] = f_12 * slh_90[k]
                   + f_3 * pc_z[k] * smh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_y, pc_x, pc_y, sli0_152, sli0_154, \
                         slh_113, slh_114, slh_183, sli1_152, sli1_154, smh_177, \
                         smh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_y[k] * sli0_152[k]
                   + f_12 * slh_113[k]
                   - f_10 * pc_y[k] * sli1_152[k];

        t_237[k] = f_11 * slh_114[k]
                   + f_3 * pc_y[k] * smh_177[k];

        t_238[k] = pb_y[k] * sli0_154[k]
                   - f_10 * pc_y[k] * sli1_154[k];

        t_239[k] = f_17 * slh_183[k]
                   + f_3 * pc_x[k] * smh_183[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_167 = buffer.data(sli0 + 167);
    const auto *sli0_168 = buffer.data(sli0 + 168);
    const auto *sli0_171 = buffer.data(sli0 + 171);
    const auto *sli0_174 = buffer.data(sli0 + 174);
    const auto *sli0_178 = buffer.data(sli0 + 178);
    const auto *sli0_180 = buffer.data(sli0 + 180);
    const auto *sli0_189 = buffer.data(sli0 + 189);

    const auto *slh_99 = buffer.data(slh + 99);
    const auto *slh_105 = buffer.data(slh + 105);
    const auto *slh_108 = buffer.data(slh + 108);
    const auto *slh_111 = buffer.data(slh + 111);
    const auto *slh_120 = buffer.data(slh + 120);
    const auto *slh_122 = buffer.data(slh + 122);
    const auto *slh_123 = buffer.data(slh + 123);
    const auto *slh_124 = buffer.data(slh + 124);
    const auto *slh_125 = buffer.data(slh + 125);
    const auto *slh_126 = buffer.data(slh + 126);
    const auto *slh_128 = buffer.data(slh + 128);
    const auto *slh_129 = buffer.data(slh + 129);
    const auto *slh_131 = buffer.data(slh + 131);
    const auto *slh_132 = buffer.data(slh + 132);
    const auto *slh_133 = buffer.data(slh + 133);
    const auto *slh_135 = buffer.data(slh + 135);
    const auto *slh_141 = buffer.data(slh + 141);
    const auto *slh_143 = buffer.data(slh + 143);
    const auto *slh_144 = buffer.data(slh + 144);
    const auto *slh_145 = buffer.data(slh + 145);
    const auto *slh_146 = buffer.data(slh + 146);
    const auto *slh_147 = buffer.data(slh + 147);
    const auto *slh_149 = buffer.data(slh + 149);
    const auto *slh_150 = buffer.data(slh + 150);
    const auto *slh_152 = buffer.data(slh + 152);
    const auto *slh_153 = buffer.data(slh + 153);
    const auto *slh_156 = buffer.data(slh + 156);
    const auto *slh_164 = buffer.data(slh + 164);
    const auto *slh_165 = buffer.data(slh + 165);
    const auto *slh_166 = buffer.data(slh + 166);
    const auto *slh_167 = buffer.data(slh + 167);
    const auto *slh_168 = buffer.data(slh + 168);
    const auto *slh_170 = buffer.data(slh + 170);
    const auto *slh_173 = buffer.data(slh + 173);
    const auto *slh_177 = buffer.data(slh + 177);
    const auto *slh_184 = buffer.data(slh + 184);
    const auto *slh_185 = buffer.data(slh + 185);
    const auto *slh_186 = buffer.data(slh + 186);
    const auto *slh_187 = buffer.data(slh + 187);
    const auto *slh_188 = buffer.data(slh + 188);
    const auto *slh_189 = buffer.data(slh + 189);
    const auto *slh_192 = buffer.data(slh + 192);
    const auto *slh_194 = buffer.data(slh + 194);
    const auto *slh_195 = buffer.data(slh + 195);
    const auto *slh_198 = buffer.data(slh + 198);
    const auto *slh_199 = buffer.data(slh + 199);
    const auto *slh_201 = buffer.data(slh + 201);
    const auto *slh_203 = buffer.data(slh + 203);
    const auto *slh_204 = buffer.data(slh + 204);
    const auto *slh_205 = buffer.data(slh + 205);
    const auto *slh_206 = buffer.data(slh + 206);
    const auto *slh_207 = buffer.data(slh + 207);
    const auto *slh_208 = buffer.data(slh + 208);
    const auto *slh_209 = buffer.data(slh + 209);
    const auto *slh_210 = buffer.data(slh + 210);
    const auto *slh_213 = buffer.data(slh + 213);
    const auto *slh_215 = buffer.data(slh + 215);
    const auto *slh_216 = buffer.data(slh + 216);
    const auto *slh_219 = buffer.data(slh + 219);
    const auto *slh_220 = buffer.data(slh + 220);
    const auto *slh_222 = buffer.data(slh + 222);
    const auto *slh_224 = buffer.data(slh + 224);
    const auto *slh_225 = buffer.data(slh + 225);
    const auto *slh_226 = buffer.data(slh + 226);
    const auto *slh_227 = buffer.data(slh + 227);
    const auto *slh_228 = buffer.data(slh + 228);
    const auto *slh_229 = buffer.data(slh + 229);
    const auto *slh_230 = buffer.data(slh + 230);
    const auto *slh_236 = buffer.data(slh + 236);
    const auto *slh_240 = buffer.data(slh + 240);
    const auto *slh_245 = buffer.data(slh + 245);
    const auto *slh_246 = buffer.data(slh + 246);
    const auto *slh_247 = buffer.data(slh + 247);
    const auto *slh_248 = buffer.data(slh + 248);
    const auto *slh_249 = buffer.data(slh + 249);
    const auto *slh_250 = buffer.data(slh + 250);
    const auto *slh_251 = buffer.data(slh + 251);
    const auto *slh_252 = buffer.data(slh + 252);
    const auto *slh_255 = buffer.data(slh + 255);
    const auto *slh_257 = buffer.data(slh + 257);
    const auto *slh_258 = buffer.data(slh + 258);
    const auto *slh_261 = buffer.data(slh + 261);
    const auto *slh_262 = buffer.data(slh + 262);
    const auto *slh_264 = buffer.data(slh + 264);
    const auto *slh_266 = buffer.data(slh + 266);

    const auto *sli1_167 = buffer.data(sli1 + 167);
    const auto *sli1_168 = buffer.data(sli1 + 168);
    const auto *sli1_171 = buffer.data(sli1 + 171);
    const auto *sli1_174 = buffer.data(sli1 + 174);
    const auto *sli1_178 = buffer.data(sli1 + 178);
    const auto *sli1_180 = buffer.data(sli1 + 180);
    const auto *sli1_189 = buffer.data(sli1 + 189);

    const auto *smg0_130 = buffer.data(smg0 + 130);
    const auto *smg0_132 = buffer.data(smg0 + 132);
    const auto *smg0_133 = buffer.data(smg0 + 133);
    const auto *smg0_134 = buffer.data(smg0 + 134);
    const auto *smg0_135 = buffer.data(smg0 + 135);
    const auto *smg0_138 = buffer.data(smg0 + 138);
    const auto *smg0_140 = buffer.data(smg0 + 140);
    const auto *smg0_141 = buffer.data(smg0 + 141);
    const auto *smg0_144 = buffer.data(smg0 + 144);
    const auto *smg0_145 = buffer.data(smg0 + 145);
    const auto *smg0_147 = buffer.data(smg0 + 147);
    const auto *smg0_148 = buffer.data(smg0 + 148);
    const auto *smg0_149 = buffer.data(smg0 + 149);
    const auto *smg0_150 = buffer.data(smg0 + 150);
    const auto *smg0_153 = buffer.data(smg0 + 153);
    const auto *smg0_155 = buffer.data(smg0 + 155);
    const auto *smg0_156 = buffer.data(smg0 + 156);
    const auto *smg0_159 = buffer.data(smg0 + 159);
    const auto *smg0_160 = buffer.data(smg0 + 160);
    const auto *smg0_162 = buffer.data(smg0 + 162);
    const auto *smg0_163 = buffer.data(smg0 + 163);
    const auto *smg0_164 = buffer.data(smg0 + 164);
    const auto *smg0_170 = buffer.data(smg0 + 170);
    const auto *smg0_174 = buffer.data(smg0 + 174);
    const auto *smg0_177 = buffer.data(smg0 + 177);
    const auto *smg0_178 = buffer.data(smg0 + 178);
    const auto *smg0_179 = buffer.data(smg0 + 179);
    const auto *smg0_180 = buffer.data(smg0 + 180);
    const auto *smg0_183 = buffer.data(smg0 + 183);
    const auto *smg0_185 = buffer.data(smg0 + 185);
    const auto *smg0_186 = buffer.data(smg0 + 186);
    const auto *smg0_189 = buffer.data(smg0 + 189);
    const auto *smg0_190 = buffer.data(smg0 + 190);
    const auto *smg0_192 = buffer.data(smg0 + 192);
    const auto *smg0_194 = buffer.data(smg0 + 194);

    const auto *smg1_130 = buffer.data(smg1 + 130);
    const auto *smg1_132 = buffer.data(smg1 + 132);
    const auto *smg1_133 = buffer.data(smg1 + 133);
    const auto *smg1_134 = buffer.data(smg1 + 134);
    const auto *smg1_135 = buffer.data(smg1 + 135);
    const auto *smg1_138 = buffer.data(smg1 + 138);
    const auto *smg1_140 = buffer.data(smg1 + 140);
    const auto *smg1_141 = buffer.data(smg1 + 141);
    const auto *smg1_144 = buffer.data(smg1 + 144);
    const auto *smg1_145 = buffer.data(smg1 + 145);
    const auto *smg1_147 = buffer.data(smg1 + 147);
    const auto *smg1_148 = buffer.data(smg1 + 148);
    const auto *smg1_149 = buffer.data(smg1 + 149);
    const auto *smg1_150 = buffer.data(smg1 + 150);
    const auto *smg1_153 = buffer.data(smg1 + 153);
    const auto *smg1_155 = buffer.data(smg1 + 155);
    const auto *smg1_156 = buffer.data(smg1 + 156);
    const auto *smg1_159 = buffer.data(smg1 + 159);
    const auto *smg1_160 = buffer.data(smg1 + 160);
    const auto *smg1_162 = buffer.data(smg1 + 162);
    const auto *smg1_163 = buffer.data(smg1 + 163);
    const auto *smg1_164 = buffer.data(smg1 + 164);
    const auto *smg1_170 = buffer.data(smg1 + 170);
    const auto *smg1_174 = buffer.data(smg1 + 174);
    const auto *smg1_177 = buffer.data(smg1 + 177);
    const auto *smg1_178 = buffer.data(smg1 + 178);
    const auto *smg1_179 = buffer.data(smg1 + 179);
    const auto *smg1_180 = buffer.data(smg1 + 180);
    const auto *smg1_183 = buffer.data(smg1 + 183);
    const auto *smg1_185 = buffer.data(smg1 + 185);
    const auto *smg1_186 = buffer.data(smg1 + 186);
    const auto *smg1_189 = buffer.data(smg1 + 189);
    const auto *smg1_190 = buffer.data(smg1 + 190);
    const auto *smg1_192 = buffer.data(smg1 + 192);
    const auto *smg1_194 = buffer.data(smg1 + 194);

    const auto *smh_183 = buffer.data(smh + 183);
    const auto *smh_184 = buffer.data(smh + 184);
    const auto *smh_185 = buffer.data(smh + 185);
    const auto *smh_186 = buffer.data(smh + 186);
    const auto *smh_187 = buffer.data(smh + 187);
    const auto *smh_188 = buffer.data(smh + 188);
    const auto *smh_189 = buffer.data(smh + 189);
    const auto *smh_191 = buffer.data(smh + 191);
    const auto *smh_192 = buffer.data(smh + 192);
    const auto *smh_194 = buffer.data(smh + 194);
    const auto *smh_195 = buffer.data(smh + 195);
    const auto *smh_198 = buffer.data(smh + 198);
    const auto *smh_199 = buffer.data(smh + 199);
    const auto *smh_201 = buffer.data(smh + 201);
    const auto *smh_203 = buffer.data(smh + 203);
    const auto *smh_204 = buffer.data(smh + 204);
    const auto *smh_205 = buffer.data(smh + 205);
    const auto *smh_206 = buffer.data(smh + 206);
    const auto *smh_207 = buffer.data(smh + 207);
    const auto *smh_208 = buffer.data(smh + 208);
    const auto *smh_209 = buffer.data(smh + 209);
    const auto *smh_210 = buffer.data(smh + 210);
    const auto *smh_212 = buffer.data(smh + 212);
    const auto *smh_213 = buffer.data(smh + 213);
    const auto *smh_215 = buffer.data(smh + 215);
    const auto *smh_216 = buffer.data(smh + 216);
    const auto *smh_219 = buffer.data(smh + 219);
    const auto *smh_220 = buffer.data(smh + 220);
    const auto *smh_222 = buffer.data(smh + 222);
    const auto *smh_224 = buffer.data(smh + 224);
    const auto *smh_225 = buffer.data(smh + 225);
    const auto *smh_226 = buffer.data(smh + 226);
    const auto *smh_227 = buffer.data(smh + 227);
    const auto *smh_228 = buffer.data(smh + 228);
    const auto *smh_229 = buffer.data(smh + 229);
    const auto *smh_230 = buffer.data(smh + 230);
    const auto *smh_231 = buffer.data(smh + 231);
    const auto *smh_233 = buffer.data(smh + 233);
    const auto *smh_234 = buffer.data(smh + 234);
    const auto *smh_236 = buffer.data(smh + 236);
    const auto *smh_237 = buffer.data(smh + 237);
    const auto *smh_240 = buffer.data(smh + 240);
    const auto *smh_245 = buffer.data(smh + 245);
    const auto *smh_246 = buffer.data(smh + 246);
    const auto *smh_247 = buffer.data(smh + 247);
    const auto *smh_248 = buffer.data(smh + 248);
    const auto *smh_249 = buffer.data(smh + 249);
    const auto *smh_250 = buffer.data(smh + 250);
    const auto *smh_251 = buffer.data(smh + 251);
    const auto *smh_252 = buffer.data(smh + 252);
    const auto *smh_254 = buffer.data(smh + 254);
    const auto *smh_255 = buffer.data(smh + 255);
    const auto *smh_257 = buffer.data(smh + 257);
    const auto *smh_258 = buffer.data(smh + 258);
    const auto *smh_261 = buffer.data(smh + 261);
    const auto *smh_262 = buffer.data(smh + 262);
    const auto *smh_264 = buffer.data(smh + 264);
    const auto *smh_266 = buffer.data(smh + 266);

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, slh_184, slh_185, slh_186, \
                         slh_187, slh_188, smh_184, smh_185, smh_186, smh_187, \
                         smh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * slh_184[k]
                   + f_3 * pc_x[k] * smh_184[k];

        t_241[k] = f_17 * slh_185[k]
                   + f_3 * pc_x[k] * smh_185[k];

        t_242[k] = f_17 * slh_186[k]
                   + f_3 * pc_x[k] * smh_186[k];

        t_243[k] = f_17 * slh_187[k]
                   + f_3 * pc_x[k] * smh_187[k];

        t_244[k] = f_17 * slh_188[k]
                   + f_3 * pc_x[k] * smh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, slh_99, slh_120, slh_122, smg0_130, \
                         smg0_132, smg1_130, smg1_132, smh_183, \
                         smh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * slh_120[k]
                   + f_1 * smg0_130[k]
                   - f_2 * smg1_130[k]
                   + f_3 * pc_y[k] * smh_183[k];

        t_246[k] = f_12 * slh_99[k]
                   + f_3 * pc_z[k] * smh_183[k];

        t_247[k] = f_11 * slh_122[k]
                   + f_4 * smg0_132[k]
                   - f_5 * smg1_132[k]
                   + f_3 * pc_y[k] * smh_185[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, slh_123, slh_124, slh_125, smg0_133, \
                         smg0_134, smg1_133, smg1_134, smh_186, smh_187, \
                         smh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * slh_123[k]
                   + f_6 * smg0_133[k]
                   - f_7 * smg1_133[k]
                   + f_3 * pc_y[k] * smh_186[k];

        t_249[k] = f_11 * slh_124[k]
                   + f_8 * smg0_134[k]
                   - f_9 * smg1_134[k]
                   + f_3 * pc_y[k] * smh_187[k];

        t_250[k] = f_11 * slh_125[k]
                   + f_3 * pc_y[k] * smh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_x, pc_y, pc_z, sli0_167, \
                         slh_105, slh_189, sli1_167, smg0_135, smg1_135, \
                         smh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * sli0_167[k]
                   - f_10 * pc_y[k] * sli1_167[k];

        t_252[k] = f_17 * slh_189[k]
                   + f_1 * smg0_135[k]
                   - f_2 * smg1_135[k]
                   + f_3 * pc_x[k] * smh_189[k];

        t_253[k] = f_3 * pc_y[k] * smh_189[k];

        t_254[k] = f_13 * slh_105[k]
                   + f_3 * pc_z[k] * smh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, slh_192, slh_194, smg0_138, \
                         smg0_140, smg1_138, smg1_140, smh_191, smh_192, \
                         smh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_17 * slh_192[k]
                   + f_4 * smg0_138[k]
                   - f_5 * smg1_138[k]
                   + f_3 * pc_x[k] * smh_192[k];

        t_256[k] = f_3 * pc_y[k] * smh_191[k];

        t_257[k] = f_17 * slh_194[k]
                   + f_4 * smg0_140[k]
                   - f_5 * smg1_140[k]
                   + f_3 * pc_x[k] * smh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, pc_z, slh_108, slh_195, smg0_141, \
                         smg1_141, smh_192, smh_194, smh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_17 * slh_195[k]
                   + f_6 * smg0_141[k]
                   - f_7 * smg1_141[k]
                   + f_3 * pc_x[k] * smh_195[k];

        t_259[k] = f_13 * slh_108[k]
                   + f_3 * pc_z[k] * smh_192[k];

        t_260[k] = f_3 * pc_y[k] * smh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, slh_111, slh_198, slh_199, smg0_144, \
                         smg0_145, smg1_144, smg1_145, smh_195, smh_198, \
                         smh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_17 * slh_198[k]
                   + f_6 * smg0_144[k]
                   - f_7 * smg1_144[k]
                   + f_3 * pc_x[k] * smh_198[k];

        t_262[k] = f_17 * slh_199[k]
                   + f_8 * smg0_145[k]
                   - f_9 * smg1_145[k]
                   + f_3 * pc_x[k] * smh_199[k];

        t_263[k] = f_13 * slh_111[k]
                   + f_3 * pc_z[k] * smh_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, slh_201, slh_203, smg0_147, \
                         smg0_149, smg1_147, smg1_149, smh_198, smh_201, \
                         smh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * slh_201[k]
                   + f_8 * smg0_147[k]
                   - f_9 * smg1_147[k]
                   + f_3 * pc_x[k] * smh_201[k];

        t_265[k] = f_3 * pc_y[k] * smh_198[k];

        t_266[k] = f_17 * slh_203[k]
                   + f_8 * smg0_149[k]
                   - f_9 * smg1_149[k]
                   + f_3 * pc_x[k] * smh_203[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pc_x, slh_204, slh_205, slh_206, \
                         slh_207, slh_208, smh_204, smh_205, smh_206, smh_207, \
                         smh_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_17 * slh_204[k]
                   + f_3 * pc_x[k] * smh_204[k];

        t_268[k] = f_17 * slh_205[k]
                   + f_3 * pc_x[k] * smh_205[k];

        t_269[k] = f_17 * slh_206[k]
                   + f_3 * pc_x[k] * smh_206[k];

        t_270[k] = f_17 * slh_207[k]
                   + f_3 * pc_x[k] * smh_207[k];

        t_271[k] = f_17 * slh_208[k]
                   + f_3 * pc_x[k] * smh_208[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, slh_120, slh_209, \
                         smg0_145, smg0_147, smg1_145, smg1_147, smh_204, smh_206, \
                         smh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_17 * slh_209[k]
                   + f_3 * pc_x[k] * smh_209[k];

        t_273[k] = f_1 * smg0_145[k]
                   - f_2 * smg1_145[k]
                   + f_3 * pc_y[k] * smh_204[k];

        t_274[k] = f_13 * slh_120[k]
                   + f_3 * pc_z[k] * smh_204[k];

        t_275[k] = f_4 * smg0_147[k]
                   - f_5 * smg1_147[k]
                   + f_3 * pc_y[k] * smh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, slh_125, smg0_148, smg0_149, \
                         smg1_148, smg1_149, smh_207, smh_208, \
                         smh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * smg0_148[k]
                   - f_7 * smg1_148[k]
                   + f_3 * pc_y[k] * smh_207[k];

        t_277[k] = f_8 * smg0_149[k]
                   - f_9 * smg1_149[k]
                   + f_3 * pc_y[k] * smh_208[k];

        t_278[k] = f_3 * pc_y[k] * smh_209[k];

        t_279[k] = f_13 * slh_125[k]
                   + f_1 * smg0_149[k]
                   - f_2 * smg1_149[k]
                   + f_3 * pc_z[k] * smh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, slh_126, slh_210, \
                         slh_213, smg0_150, smg0_153, smg1_150, smg1_153, smh_210, \
                         smh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_18 * slh_210[k]
                   + f_1 * smg0_150[k]
                   - f_2 * smg1_150[k]
                   + f_3 * pc_x[k] * smh_210[k];

        t_281[k] = f_14 * slh_126[k]
                   + f_3 * pc_y[k] * smh_210[k];

        t_282[k] = f_3 * pc_z[k] * smh_210[k];

        t_283[k] = f_18 * slh_213[k]
                   + f_4 * smg0_153[k]
                   - f_5 * smg1_153[k]
                   + f_3 * pc_x[k] * smh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_y, slh_128, slh_215, slh_216, smg0_155, \
                         smg0_156, smg1_155, smg1_156, smh_212, smh_215, \
                         smh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_14 * slh_128[k]
                   + f_3 * pc_y[k] * smh_212[k];

        t_285[k] = f_18 * slh_215[k]
                   + f_4 * smg0_155[k]
                   - f_5 * smg1_155[k]
                   + f_3 * pc_x[k] * smh_215[k];

        t_286[k] = f_18 * slh_216[k]
                   + f_6 * smg0_156[k]
                   - f_7 * smg1_156[k]
                   + f_3 * pc_x[k] * smh_216[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, pc_z, slh_131, slh_219, smg0_159, \
                         smg1_159, smh_213, smh_215, smh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * smh_213[k];

        t_288[k] = f_14 * slh_131[k]
                   + f_3 * pc_y[k] * smh_215[k];

        t_289[k] = f_18 * slh_219[k]
                   + f_6 * smg0_159[k]
                   - f_7 * smg1_159[k]
                   + f_3 * pc_x[k] * smh_219[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, slh_220, slh_222, smg0_160, \
                         smg0_162, smg1_160, smg1_162, smh_216, smh_220, \
                         smh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_18 * slh_220[k]
                   + f_8 * smg0_160[k]
                   - f_9 * smg1_160[k]
                   + f_3 * pc_x[k] * smh_220[k];

        t_291[k] = f_3 * pc_z[k] * smh_216[k];

        t_292[k] = f_18 * slh_222[k]
                   + f_8 * smg0_162[k]
                   - f_9 * smg1_162[k]
                   + f_3 * pc_x[k] * smh_222[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pc_x, pc_y, slh_135, slh_224, slh_225, \
                         slh_226, smg0_164, smg1_164, smh_219, smh_224, smh_225, \
                         smh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_14 * slh_135[k]
                   + f_3 * pc_y[k] * smh_219[k];

        t_294[k] = f_18 * slh_224[k]
                   + f_8 * smg0_164[k]
                   - f_9 * smg1_164[k]
                   + f_3 * pc_x[k] * smh_224[k];

        t_295[k] = f_18 * slh_225[k]
                   + f_3 * pc_x[k] * smh_225[k];

        t_296[k] = f_18 * slh_226[k]
                   + f_3 * pc_x[k] * smh_226[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pc_x, slh_227, slh_228, slh_229, slh_230, \
                         smh_227, smh_228, smh_229, smh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_18 * slh_227[k]
                   + f_3 * pc_x[k] * smh_227[k];

        t_298[k] = f_18 * slh_228[k]
                   + f_3 * pc_x[k] * smh_228[k];

        t_299[k] = f_18 * slh_229[k]
                   + f_3 * pc_x[k] * smh_229[k];

        t_300[k] = f_18 * slh_230[k]
                   + f_3 * pc_x[k] * smh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, pc_y, pc_z, slh_141, slh_143, smg0_160, \
                         smg0_162, smg1_160, smg1_162, smh_225, \
                         smh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * slh_141[k]
                   + f_1 * smg0_160[k]
                   - f_2 * smg1_160[k]
                   + f_3 * pc_y[k] * smh_225[k];

        t_302[k] = f_3 * pc_z[k] * smh_225[k];

        t_303[k] = f_14 * slh_143[k]
                   + f_4 * smg0_162[k]
                   - f_5 * smg1_162[k]
                   + f_3 * pc_y[k] * smh_227[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_y, pc_z, slh_144, slh_145, slh_146, \
                         smg0_163, smg0_164, smg1_163, smg1_164, smh_228, smh_229, \
                         smh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_14 * slh_144[k]
                   + f_6 * smg0_163[k]
                   - f_7 * smg1_163[k]
                   + f_3 * pc_y[k] * smh_228[k];

        t_305[k] = f_14 * slh_145[k]
                   + f_8 * smg0_164[k]
                   - f_9 * smg1_164[k]
                   + f_3 * pc_y[k] * smh_229[k];

        t_306[k] = f_14 * slh_146[k]
                   + f_3 * pc_y[k] * smh_230[k];

        t_307[k] = f_1 * smg0_164[k]
                   - f_2 * smg1_164[k]
                   + f_3 * pc_z[k] * smh_230[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, sli0_168, sli0_171, \
                         slh_126, slh_147, sli1_168, sli1_171, \
                         smh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * sli0_168[k]
                   - f_10 * pc_z[k] * sli1_168[k];

        t_309[k] = f_13 * slh_147[k]
                   + f_3 * pc_y[k] * smh_231[k];

        t_310[k] = f_11 * slh_126[k]
                   + f_3 * pc_z[k] * smh_231[k];

        t_311[k] = pb_z[k] * sli0_171[k]
                   - f_10 * pc_z[k] * sli1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_z, pc_x, pc_y, pc_z, sli0_174, slh_149, \
                         slh_236, sli1_174, smg0_170, smg1_170, smh_233, \
                         smh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * slh_149[k]
                   + f_3 * pc_y[k] * smh_233[k];

        t_313[k] = f_18 * slh_236[k]
                   + f_4 * smg0_170[k]
                   - f_5 * smg1_170[k]
                   + f_3 * pc_x[k] * smh_236[k];

        t_314[k] = pb_z[k] * sli0_174[k]
                   - f_10 * pc_z[k] * sli1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pc_x, pc_y, pc_z, slh_129, slh_152, slh_240, \
                         smg0_174, smg1_174, smh_234, smh_236, \
                         smh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * slh_129[k]
                   + f_3 * pc_z[k] * smh_234[k];

        t_316[k] = f_13 * slh_152[k]
                   + f_3 * pc_y[k] * smh_236[k];

        t_317[k] = f_18 * slh_240[k]
                   + f_6 * smg0_174[k]
                   - f_7 * smg1_174[k]
                   + f_3 * pc_x[k] * smh_240[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_z, pc_y, pc_z, sli0_178, sli0_180, \
                         slh_132, slh_133, slh_156, sli1_178, sli1_180, smh_237, \
                         smh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * sli0_178[k]
                   - f_10 * pc_z[k] * sli1_178[k];

        t_319[k] = f_11 * slh_132[k]
                   + f_3 * pc_z[k] * smh_237[k];

        t_320[k] = pb_z[k] * sli0_180[k]
                   + f_12 * slh_133[k]
                   - f_10 * pc_z[k] * sli1_180[k];

        t_321[k] = f_13 * slh_156[k]
                   + f_3 * pc_y[k] * smh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, slh_245, slh_246, slh_247, slh_248, \
                         smg0_179, smg1_179, smh_245, smh_246, smh_247, \
                         smh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_18 * slh_245[k]
                   + f_8 * smg0_179[k]
                   - f_9 * smg1_179[k]
                   + f_3 * pc_x[k] * smh_245[k];

        t_323[k] = f_18 * slh_246[k]
                   + f_3 * pc_x[k] * smh_246[k];

        t_324[k] = f_18 * slh_247[k]
                   + f_3 * pc_x[k] * smh_247[k];

        t_325[k] = f_18 * slh_248[k]
                   + f_3 * pc_x[k] * smh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pb_z, pc_x, pc_z, sli0_189, slh_249, \
                         slh_250, slh_251, sli1_189, smh_249, smh_250, \
                         smh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_18 * slh_249[k]
                   + f_3 * pc_x[k] * smh_249[k];

        t_327[k] = f_18 * slh_250[k]
                   + f_3 * pc_x[k] * smh_250[k];

        t_328[k] = f_18 * slh_251[k]
                   + f_3 * pc_x[k] * smh_251[k];

        t_329[k] = pb_z[k] * sli0_189[k]
                   - f_10 * pc_z[k] * sli1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, slh_141, slh_164, slh_165, smg0_177, \
                         smg0_178, smg1_177, smg1_178, smh_246, smh_248, \
                         smh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * slh_141[k]
                   + f_3 * pc_z[k] * smh_246[k];

        t_331[k] = f_13 * slh_164[k]
                   + f_4 * smg0_177[k]
                   - f_5 * smg1_177[k]
                   + f_3 * pc_y[k] * smh_248[k];

        t_332[k] = f_13 * slh_165[k]
                   + f_6 * smg0_178[k]
                   - f_7 * smg1_178[k]
                   + f_3 * pc_y[k] * smh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, slh_146, slh_166, slh_167, smg0_179, \
                         smg1_179, smh_250, smh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * slh_166[k]
                   + f_8 * smg0_179[k]
                   - f_9 * smg1_179[k]
                   + f_3 * pc_y[k] * smh_250[k];

        t_334[k] = f_13 * slh_167[k]
                   + f_3 * pc_y[k] * smh_251[k];

        t_335[k] = f_11 * slh_146[k]
                   + f_1 * smg0_179[k]
                   - f_2 * smg1_179[k]
                   + f_3 * pc_z[k] * smh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, slh_147, slh_168, slh_252, \
                         smg0_180, smg1_180, smh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_18 * slh_252[k]
                   + f_1 * smg0_180[k]
                   - f_2 * smg1_180[k]
                   + f_3 * pc_x[k] * smh_252[k];

        t_337[k] = f_12 * slh_168[k]
                   + f_3 * pc_y[k] * smh_252[k];

        t_338[k] = f_12 * slh_147[k]
                   + f_3 * pc_z[k] * smh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, slh_170, slh_255, slh_257, smg0_183, \
                         smg0_185, smg1_183, smg1_185, smh_254, smh_255, \
                         smh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_18 * slh_255[k]
                   + f_4 * smg0_183[k]
                   - f_5 * smg1_183[k]
                   + f_3 * pc_x[k] * smh_255[k];

        t_340[k] = f_12 * slh_170[k]
                   + f_3 * pc_y[k] * smh_254[k];

        t_341[k] = f_18 * slh_257[k]
                   + f_4 * smg0_185[k]
                   - f_5 * smg1_185[k]
                   + f_3 * pc_x[k] * smh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, slh_150, slh_173, slh_258, \
                         smg0_186, smg1_186, smh_255, smh_257, \
                         smh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_18 * slh_258[k]
                   + f_6 * smg0_186[k]
                   - f_7 * smg1_186[k]
                   + f_3 * pc_x[k] * smh_258[k];

        t_343[k] = f_12 * slh_150[k]
                   + f_3 * pc_z[k] * smh_255[k];

        t_344[k] = f_12 * slh_173[k]
                   + f_3 * pc_y[k] * smh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, slh_153, slh_261, slh_262, smg0_189, \
                         smg0_190, smg1_189, smg1_190, smh_258, smh_261, \
                         smh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_18 * slh_261[k]
                   + f_6 * smg0_189[k]
                   - f_7 * smg1_189[k]
                   + f_3 * pc_x[k] * smh_261[k];

        t_346[k] = f_18 * slh_262[k]
                   + f_8 * smg0_190[k]
                   - f_9 * smg1_190[k]
                   + f_3 * pc_x[k] * smh_262[k];

        t_347[k] = f_12 * slh_153[k]
                   + f_3 * pc_z[k] * smh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, slh_177, slh_264, slh_266, smg0_192, \
                         smg0_194, smg1_192, smg1_194, smh_261, smh_264, \
                         smh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_18 * slh_264[k]
                   + f_8 * smg0_192[k]
                   - f_9 * smg1_192[k]
                   + f_3 * pc_x[k] * smh_264[k];

        t_349[k] = f_12 * slh_177[k]
                   + f_3 * pc_y[k] * smh_261[k];

        t_350[k] = f_18 * slh_266[k]
                   + f_8 * smg0_194[k]
                   - f_9 * smg1_194[k]
                   + f_3 * pc_x[k] * smh_266[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_252 = buffer.data(sli0 + 252);
    const auto *sli0_255 = buffer.data(sli0 + 255);
    const auto *sli0_257 = buffer.data(sli0 + 257);
    const auto *sli0_258 = buffer.data(sli0 + 258);
    const auto *sli0_261 = buffer.data(sli0 + 261);
    const auto *sli0_262 = buffer.data(sli0 + 262);
    const auto *sli0_264 = buffer.data(sli0 + 264);
    const auto *sli0_266 = buffer.data(sli0 + 266);
    const auto *sli0_279 = buffer.data(sli0 + 279);
    const auto *sli0_280 = buffer.data(sli0 + 280);
    const auto *sli0_283 = buffer.data(sli0 + 283);
    const auto *sli0_286 = buffer.data(sli0 + 286);
    const auto *sli0_290 = buffer.data(sli0 + 290);
    const auto *sli0_292 = buffer.data(sli0 + 292);

    const auto *slh_162 = buffer.data(slh + 162);
    const auto *slh_167 = buffer.data(slh + 167);
    const auto *slh_168 = buffer.data(slh + 168);
    const auto *slh_171 = buffer.data(slh + 171);
    const auto *slh_174 = buffer.data(slh + 174);
    const auto *slh_183 = buffer.data(slh + 183);
    const auto *slh_185 = buffer.data(slh + 185);
    const auto *slh_186 = buffer.data(slh + 186);
    const auto *slh_187 = buffer.data(slh + 187);
    const auto *slh_188 = buffer.data(slh + 188);
    const auto *slh_189 = buffer.data(slh + 189);
    const auto *slh_190 = buffer.data(slh + 190);
    const auto *slh_191 = buffer.data(slh + 191);
    const auto *slh_192 = buffer.data(slh + 192);
    const auto *slh_194 = buffer.data(slh + 194);
    const auto *slh_195 = buffer.data(slh + 195);
    const auto *slh_197 = buffer.data(slh + 197);
    const auto *slh_198 = buffer.data(slh + 198);
    const auto *slh_204 = buffer.data(slh + 204);
    const auto *slh_206 = buffer.data(slh + 206);
    const auto *slh_207 = buffer.data(slh + 207);
    const auto *slh_208 = buffer.data(slh + 208);
    const auto *slh_209 = buffer.data(slh + 209);
    const auto *slh_210 = buffer.data(slh + 210);
    const auto *slh_212 = buffer.data(slh + 212);
    const auto *slh_213 = buffer.data(slh + 213);
    const auto *slh_215 = buffer.data(slh + 215);
    const auto *slh_216 = buffer.data(slh + 216);
    const auto *slh_217 = buffer.data(slh + 217);
    const auto *slh_219 = buffer.data(slh + 219);
    const auto *slh_225 = buffer.data(slh + 225);
    const auto *slh_227 = buffer.data(slh + 227);
    const auto *slh_228 = buffer.data(slh + 228);
    const auto *slh_229 = buffer.data(slh + 229);
    const auto *slh_230 = buffer.data(slh + 230);
    const auto *slh_231 = buffer.data(slh + 231);
    const auto *slh_233 = buffer.data(slh + 233);
    const auto *slh_236 = buffer.data(slh + 236);
    const auto *slh_240 = buffer.data(slh + 240);
    const auto *slh_267 = buffer.data(slh + 267);
    const auto *slh_268 = buffer.data(slh + 268);
    const auto *slh_269 = buffer.data(slh + 269);
    const auto *slh_270 = buffer.data(slh + 270);
    const auto *slh_271 = buffer.data(slh + 271);
    const auto *slh_272 = buffer.data(slh + 272);
    const auto *slh_288 = buffer.data(slh + 288);
    const auto *slh_289 = buffer.data(slh + 289);
    const auto *slh_290 = buffer.data(slh + 290);
    const auto *slh_291 = buffer.data(slh + 291);
    const auto *slh_292 = buffer.data(slh + 292);
    const auto *slh_293 = buffer.data(slh + 293);
    const auto *slh_294 = buffer.data(slh + 294);
    const auto *slh_297 = buffer.data(slh + 297);
    const auto *slh_299 = buffer.data(slh + 299);
    const auto *slh_300 = buffer.data(slh + 300);
    const auto *slh_303 = buffer.data(slh + 303);
    const auto *slh_304 = buffer.data(slh + 304);
    const auto *slh_306 = buffer.data(slh + 306);
    const auto *slh_308 = buffer.data(slh + 308);
    const auto *slh_309 = buffer.data(slh + 309);
    const auto *slh_310 = buffer.data(slh + 310);
    const auto *slh_311 = buffer.data(slh + 311);
    const auto *slh_312 = buffer.data(slh + 312);
    const auto *slh_313 = buffer.data(slh + 313);
    const auto *slh_314 = buffer.data(slh + 314);
    const auto *slh_315 = buffer.data(slh + 315);
    const auto *slh_318 = buffer.data(slh + 318);
    const auto *slh_320 = buffer.data(slh + 320);
    const auto *slh_321 = buffer.data(slh + 321);
    const auto *slh_324 = buffer.data(slh + 324);
    const auto *slh_325 = buffer.data(slh + 325);
    const auto *slh_327 = buffer.data(slh + 327);
    const auto *slh_329 = buffer.data(slh + 329);
    const auto *slh_330 = buffer.data(slh + 330);
    const auto *slh_331 = buffer.data(slh + 331);
    const auto *slh_332 = buffer.data(slh + 332);
    const auto *slh_333 = buffer.data(slh + 333);
    const auto *slh_334 = buffer.data(slh + 334);
    const auto *slh_335 = buffer.data(slh + 335);
    const auto *slh_341 = buffer.data(slh + 341);
    const auto *slh_345 = buffer.data(slh + 345);
    const auto *slh_350 = buffer.data(slh + 350);
    const auto *slh_351 = buffer.data(slh + 351);
    const auto *slh_352 = buffer.data(slh + 352);
    const auto *slh_353 = buffer.data(slh + 353);

    const auto *sli1_252 = buffer.data(sli1 + 252);
    const auto *sli1_255 = buffer.data(sli1 + 255);
    const auto *sli1_257 = buffer.data(sli1 + 257);
    const auto *sli1_258 = buffer.data(sli1 + 258);
    const auto *sli1_261 = buffer.data(sli1 + 261);
    const auto *sli1_262 = buffer.data(sli1 + 262);
    const auto *sli1_264 = buffer.data(sli1 + 264);
    const auto *sli1_266 = buffer.data(sli1 + 266);
    const auto *sli1_279 = buffer.data(sli1 + 279);
    const auto *sli1_280 = buffer.data(sli1 + 280);
    const auto *sli1_283 = buffer.data(sli1 + 283);
    const auto *sli1_286 = buffer.data(sli1 + 286);
    const auto *sli1_290 = buffer.data(sli1 + 290);
    const auto *sli1_292 = buffer.data(sli1 + 292);

    const auto *smg0_190 = buffer.data(smg0 + 190);
    const auto *smg0_192 = buffer.data(smg0 + 192);
    const auto *smg0_193 = buffer.data(smg0 + 193);
    const auto *smg0_194 = buffer.data(smg0 + 194);
    const auto *smg0_205 = buffer.data(smg0 + 205);
    const auto *smg0_207 = buffer.data(smg0 + 207);
    const auto *smg0_208 = buffer.data(smg0 + 208);
    const auto *smg0_209 = buffer.data(smg0 + 209);
    const auto *smg0_210 = buffer.data(smg0 + 210);
    const auto *smg0_213 = buffer.data(smg0 + 213);
    const auto *smg0_215 = buffer.data(smg0 + 215);
    const auto *smg0_216 = buffer.data(smg0 + 216);
    const auto *smg0_219 = buffer.data(smg0 + 219);
    const auto *smg0_220 = buffer.data(smg0 + 220);
    const auto *smg0_222 = buffer.data(smg0 + 222);
    const auto *smg0_223 = buffer.data(smg0 + 223);
    const auto *smg0_224 = buffer.data(smg0 + 224);
    const auto *smg0_225 = buffer.data(smg0 + 225);
    const auto *smg0_228 = buffer.data(smg0 + 228);
    const auto *smg0_230 = buffer.data(smg0 + 230);
    const auto *smg0_231 = buffer.data(smg0 + 231);
    const auto *smg0_234 = buffer.data(smg0 + 234);
    const auto *smg0_235 = buffer.data(smg0 + 235);
    const auto *smg0_237 = buffer.data(smg0 + 237);
    const auto *smg0_238 = buffer.data(smg0 + 238);
    const auto *smg0_239 = buffer.data(smg0 + 239);
    const auto *smg0_245 = buffer.data(smg0 + 245);
    const auto *smg0_249 = buffer.data(smg0 + 249);
    const auto *smg0_254 = buffer.data(smg0 + 254);

    const auto *smg1_190 = buffer.data(smg1 + 190);
    const auto *smg1_192 = buffer.data(smg1 + 192);
    const auto *smg1_193 = buffer.data(smg1 + 193);
    const auto *smg1_194 = buffer.data(smg1 + 194);
    const auto *smg1_205 = buffer.data(smg1 + 205);
    const auto *smg1_207 = buffer.data(smg1 + 207);
    const auto *smg1_208 = buffer.data(smg1 + 208);
    const auto *smg1_209 = buffer.data(smg1 + 209);
    const auto *smg1_210 = buffer.data(smg1 + 210);
    const auto *smg1_213 = buffer.data(smg1 + 213);
    const auto *smg1_215 = buffer.data(smg1 + 215);
    const auto *smg1_216 = buffer.data(smg1 + 216);
    const auto *smg1_219 = buffer.data(smg1 + 219);
    const auto *smg1_220 = buffer.data(smg1 + 220);
    const auto *smg1_222 = buffer.data(smg1 + 222);
    const auto *smg1_223 = buffer.data(smg1 + 223);
    const auto *smg1_224 = buffer.data(smg1 + 224);
    const auto *smg1_225 = buffer.data(smg1 + 225);
    const auto *smg1_228 = buffer.data(smg1 + 228);
    const auto *smg1_230 = buffer.data(smg1 + 230);
    const auto *smg1_231 = buffer.data(smg1 + 231);
    const auto *smg1_234 = buffer.data(smg1 + 234);
    const auto *smg1_235 = buffer.data(smg1 + 235);
    const auto *smg1_237 = buffer.data(smg1 + 237);
    const auto *smg1_238 = buffer.data(smg1 + 238);
    const auto *smg1_239 = buffer.data(smg1 + 239);
    const auto *smg1_245 = buffer.data(smg1 + 245);
    const auto *smg1_249 = buffer.data(smg1 + 249);
    const auto *smg1_254 = buffer.data(smg1 + 254);

    const auto *smh_267 = buffer.data(smh + 267);
    const auto *smh_268 = buffer.data(smh + 268);
    const auto *smh_269 = buffer.data(smh + 269);
    const auto *smh_270 = buffer.data(smh + 270);
    const auto *smh_271 = buffer.data(smh + 271);
    const auto *smh_272 = buffer.data(smh + 272);
    const auto *smh_273 = buffer.data(smh + 273);
    const auto *smh_275 = buffer.data(smh + 275);
    const auto *smh_276 = buffer.data(smh + 276);
    const auto *smh_278 = buffer.data(smh + 278);
    const auto *smh_279 = buffer.data(smh + 279);
    const auto *smh_282 = buffer.data(smh + 282);
    const auto *smh_288 = buffer.data(smh + 288);
    const auto *smh_289 = buffer.data(smh + 289);
    const auto *smh_290 = buffer.data(smh + 290);
    const auto *smh_291 = buffer.data(smh + 291);
    const auto *smh_292 = buffer.data(smh + 292);
    const auto *smh_293 = buffer.data(smh + 293);
    const auto *smh_294 = buffer.data(smh + 294);
    const auto *smh_296 = buffer.data(smh + 296);
    const auto *smh_297 = buffer.data(smh + 297);
    const auto *smh_299 = buffer.data(smh + 299);
    const auto *smh_300 = buffer.data(smh + 300);
    const auto *smh_303 = buffer.data(smh + 303);
    const auto *smh_304 = buffer.data(smh + 304);
    const auto *smh_306 = buffer.data(smh + 306);
    const auto *smh_308 = buffer.data(smh + 308);
    const auto *smh_309 = buffer.data(smh + 309);
    const auto *smh_310 = buffer.data(smh + 310);
    const auto *smh_311 = buffer.data(smh + 311);
    const auto *smh_312 = buffer.data(smh + 312);
    const auto *smh_313 = buffer.data(smh + 313);
    const auto *smh_314 = buffer.data(smh + 314);
    const auto *smh_315 = buffer.data(smh + 315);
    const auto *smh_317 = buffer.data(smh + 317);
    const auto *smh_318 = buffer.data(smh + 318);
    const auto *smh_320 = buffer.data(smh + 320);
    const auto *smh_321 = buffer.data(smh + 321);
    const auto *smh_324 = buffer.data(smh + 324);
    const auto *smh_325 = buffer.data(smh + 325);
    const auto *smh_327 = buffer.data(smh + 327);
    const auto *smh_329 = buffer.data(smh + 329);
    const auto *smh_330 = buffer.data(smh + 330);
    const auto *smh_331 = buffer.data(smh + 331);
    const auto *smh_332 = buffer.data(smh + 332);
    const auto *smh_333 = buffer.data(smh + 333);
    const auto *smh_334 = buffer.data(smh + 334);
    const auto *smh_335 = buffer.data(smh + 335);
    const auto *smh_336 = buffer.data(smh + 336);
    const auto *smh_338 = buffer.data(smh + 338);
    const auto *smh_339 = buffer.data(smh + 339);
    const auto *smh_341 = buffer.data(smh + 341);
    const auto *smh_342 = buffer.data(smh + 342);
    const auto *smh_345 = buffer.data(smh + 345);
    const auto *smh_350 = buffer.data(smh + 350);
    const auto *smh_351 = buffer.data(smh + 351);
    const auto *smh_352 = buffer.data(smh + 352);
    const auto *smh_353 = buffer.data(smh + 353);

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, slh_267, slh_268, slh_269, \
                         slh_270, slh_271, smh_267, smh_268, smh_269, smh_270, \
                         smh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_18 * slh_267[k]
                   + f_3 * pc_x[k] * smh_267[k];

        t_352[k] = f_18 * slh_268[k]
                   + f_3 * pc_x[k] * smh_268[k];

        t_353[k] = f_18 * slh_269[k]
                   + f_3 * pc_x[k] * smh_269[k];

        t_354[k] = f_18 * slh_270[k]
                   + f_3 * pc_x[k] * smh_270[k];

        t_355[k] = f_18 * slh_271[k]
                   + f_3 * pc_x[k] * smh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, slh_162, slh_183, slh_272, \
                         smg0_190, smg1_190, smh_267, smh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_18 * slh_272[k]
                   + f_3 * pc_x[k] * smh_272[k];

        t_357[k] = f_12 * slh_183[k]
                   + f_1 * smg0_190[k]
                   - f_2 * smg1_190[k]
                   + f_3 * pc_y[k] * smh_267[k];

        t_358[k] = f_12 * slh_162[k]
                   + f_3 * pc_z[k] * smh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, slh_185, slh_186, slh_187, smg0_192, \
                         smg0_193, smg0_194, smg1_192, smg1_193, smg1_194, smh_269, smh_270, \
                         smh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * slh_185[k]
                   + f_4 * smg0_192[k]
                   - f_5 * smg1_192[k]
                   + f_3 * pc_y[k] * smh_269[k];

        t_360[k] = f_12 * slh_186[k]
                   + f_6 * smg0_193[k]
                   - f_7 * smg1_193[k]
                   + f_3 * pc_y[k] * smh_270[k];

        t_361[k] = f_12 * slh_187[k]
                   + f_8 * smg0_194[k]
                   - f_9 * smg1_194[k]
                   + f_3 * pc_y[k] * smh_271[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pb_y, pc_y, pc_z, sli0_252, slh_167, \
                         slh_188, slh_189, sli1_252, smg0_194, smg1_194, smh_272, \
                         smh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * slh_188[k]
                   + f_3 * pc_y[k] * smh_272[k];

        t_363[k] = f_12 * slh_167[k]
                   + f_1 * smg0_194[k]
                   - f_2 * smg1_194[k]
                   + f_3 * pc_z[k] * smh_272[k];

        t_364[k] = pb_y[k] * sli0_252[k]
                   - f_10 * pc_y[k] * sli1_252[k];

        t_365[k] = f_11 * slh_189[k]
                   + f_3 * pc_y[k] * smh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_y, pc_y, pc_z, sli0_255, sli0_257, \
                         slh_168, slh_190, slh_191, sli1_255, sli1_257, smh_273, \
                         smh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * slh_168[k]
                   + f_3 * pc_z[k] * smh_273[k];

        t_367[k] = pb_y[k] * sli0_255[k]
                   + f_12 * slh_190[k]
                   - f_10 * pc_y[k] * sli1_255[k];

        t_368[k] = f_11 * slh_191[k]
                   + f_3 * pc_y[k] * smh_275[k];

        t_369[k] = pb_y[k] * sli0_257[k]
                   - f_10 * pc_y[k] * sli1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_y, pc_y, pc_z, sli0_258, sli0_261, \
                         slh_171, slh_192, slh_194, sli1_258, sli1_261, smh_276, \
                         smh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_y[k] * sli0_258[k]
                   + f_13 * slh_192[k]
                   - f_10 * pc_y[k] * sli1_258[k];

        t_371[k] = f_13 * slh_171[k]
                   + f_3 * pc_z[k] * smh_276[k];

        t_372[k] = f_11 * slh_194[k]
                   + f_3 * pc_y[k] * smh_278[k];

        t_373[k] = pb_y[k] * sli0_261[k]
                   - f_10 * pc_y[k] * sli1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pb_y, pc_y, pc_z, sli0_262, sli0_264, slh_174, \
                         slh_195, slh_197, sli1_262, sli1_264, \
                         smh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_y[k] * sli0_262[k]
                   + f_14 * slh_195[k]
                   - f_10 * pc_y[k] * sli1_262[k];

        t_375[k] = f_13 * slh_174[k]
                   + f_3 * pc_z[k] * smh_279[k];

        t_376[k] = pb_y[k] * sli0_264[k]
                   + f_12 * slh_197[k]
                   - f_10 * pc_y[k] * sli1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_y, pc_x, pc_y, sli0_266, slh_198, \
                         slh_288, slh_289, sli1_266, smh_282, smh_288, \
                         smh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * slh_198[k]
                   + f_3 * pc_y[k] * smh_282[k];

        t_378[k] = pb_y[k] * sli0_266[k]
                   - f_10 * pc_y[k] * sli1_266[k];

        t_379[k] = f_18 * slh_288[k]
                   + f_3 * pc_x[k] * smh_288[k];

        t_380[k] = f_18 * slh_289[k]
                   + f_3 * pc_x[k] * smh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, slh_290, slh_291, slh_292, slh_293, \
                         smh_290, smh_291, smh_292, smh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_18 * slh_290[k]
                   + f_3 * pc_x[k] * smh_290[k];

        t_382[k] = f_18 * slh_291[k]
                   + f_3 * pc_x[k] * smh_291[k];

        t_383[k] = f_18 * slh_292[k]
                   + f_3 * pc_x[k] * smh_292[k];

        t_384[k] = f_18 * slh_293[k]
                   + f_3 * pc_x[k] * smh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, slh_183, slh_204, slh_206, smg0_205, \
                         smg0_207, smg1_205, smg1_207, smh_288, \
                         smh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * slh_204[k]
                   + f_1 * smg0_205[k]
                   - f_2 * smg1_205[k]
                   + f_3 * pc_y[k] * smh_288[k];

        t_386[k] = f_13 * slh_183[k]
                   + f_3 * pc_z[k] * smh_288[k];

        t_387[k] = f_11 * slh_206[k]
                   + f_4 * smg0_207[k]
                   - f_5 * smg1_207[k]
                   + f_3 * pc_y[k] * smh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, slh_207, slh_208, slh_209, smg0_208, \
                         smg0_209, smg1_208, smg1_209, smh_291, smh_292, \
                         smh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * slh_207[k]
                   + f_6 * smg0_208[k]
                   - f_7 * smg1_208[k]
                   + f_3 * pc_y[k] * smh_291[k];

        t_389[k] = f_11 * slh_208[k]
                   + f_8 * smg0_209[k]
                   - f_9 * smg1_209[k]
                   + f_3 * pc_y[k] * smh_292[k];

        t_390[k] = f_11 * slh_209[k]
                   + f_3 * pc_y[k] * smh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_y, pc_x, pc_y, pc_z, sli0_279, \
                         slh_189, slh_294, sli1_279, smg0_210, smg1_210, \
                         smh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pb_y[k] * sli0_279[k]
                   - f_10 * pc_y[k] * sli1_279[k];

        t_392[k] = f_18 * slh_294[k]
                   + f_1 * smg0_210[k]
                   - f_2 * smg1_210[k]
                   + f_3 * pc_x[k] * smh_294[k];

        t_393[k] = f_3 * pc_y[k] * smh_294[k];

        t_394[k] = f_14 * slh_189[k]
                   + f_3 * pc_z[k] * smh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, slh_297, slh_299, smg0_213, \
                         smg0_215, smg1_213, smg1_215, smh_296, smh_297, \
                         smh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_18 * slh_297[k]
                   + f_4 * smg0_213[k]
                   - f_5 * smg1_213[k]
                   + f_3 * pc_x[k] * smh_297[k];

        t_396[k] = f_3 * pc_y[k] * smh_296[k];

        t_397[k] = f_18 * slh_299[k]
                   + f_4 * smg0_215[k]
                   - f_5 * smg1_215[k]
                   + f_3 * pc_x[k] * smh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_x, pc_y, pc_z, slh_192, slh_300, smg0_216, \
                         smg1_216, smh_297, smh_299, smh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_18 * slh_300[k]
                   + f_6 * smg0_216[k]
                   - f_7 * smg1_216[k]
                   + f_3 * pc_x[k] * smh_300[k];

        t_399[k] = f_14 * slh_192[k]
                   + f_3 * pc_z[k] * smh_297[k];

        t_400[k] = f_3 * pc_y[k] * smh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_z, slh_195, slh_303, slh_304, smg0_219, \
                         smg0_220, smg1_219, smg1_220, smh_300, smh_303, \
                         smh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_18 * slh_303[k]
                   + f_6 * smg0_219[k]
                   - f_7 * smg1_219[k]
                   + f_3 * pc_x[k] * smh_303[k];

        t_402[k] = f_18 * slh_304[k]
                   + f_8 * smg0_220[k]
                   - f_9 * smg1_220[k]
                   + f_3 * pc_x[k] * smh_304[k];

        t_403[k] = f_14 * slh_195[k]
                   + f_3 * pc_z[k] * smh_300[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pc_x, pc_y, slh_306, slh_308, smg0_222, \
                         smg0_224, smg1_222, smg1_224, smh_303, smh_306, \
                         smh_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_18 * slh_306[k]
                   + f_8 * smg0_222[k]
                   - f_9 * smg1_222[k]
                   + f_3 * pc_x[k] * smh_306[k];

        t_405[k] = f_3 * pc_y[k] * smh_303[k];

        t_406[k] = f_18 * slh_308[k]
                   + f_8 * smg0_224[k]
                   - f_9 * smg1_224[k]
                   + f_3 * pc_x[k] * smh_308[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pc_x, slh_309, slh_310, slh_311, \
                         slh_312, slh_313, smh_309, smh_310, smh_311, smh_312, \
                         smh_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_18 * slh_309[k]
                   + f_3 * pc_x[k] * smh_309[k];

        t_408[k] = f_18 * slh_310[k]
                   + f_3 * pc_x[k] * smh_310[k];

        t_409[k] = f_18 * slh_311[k]
                   + f_3 * pc_x[k] * smh_311[k];

        t_410[k] = f_18 * slh_312[k]
                   + f_3 * pc_x[k] * smh_312[k];

        t_411[k] = f_18 * slh_313[k]
                   + f_3 * pc_x[k] * smh_313[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pc_x, pc_y, pc_z, slh_204, slh_314, \
                         smg0_220, smg0_222, smg1_220, smg1_222, smh_309, smh_311, \
                         smh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_18 * slh_314[k]
                   + f_3 * pc_x[k] * smh_314[k];

        t_413[k] = f_1 * smg0_220[k]
                   - f_2 * smg1_220[k]
                   + f_3 * pc_y[k] * smh_309[k];

        t_414[k] = f_14 * slh_204[k]
                   + f_3 * pc_z[k] * smh_309[k];

        t_415[k] = f_4 * smg0_222[k]
                   - f_5 * smg1_222[k]
                   + f_3 * pc_y[k] * smh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, slh_209, smg0_223, smg0_224, \
                         smg1_223, smg1_224, smh_312, smh_313, \
                         smh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * smg0_223[k]
                   - f_7 * smg1_223[k]
                   + f_3 * pc_y[k] * smh_312[k];

        t_417[k] = f_8 * smg0_224[k]
                   - f_9 * smg1_224[k]
                   + f_3 * pc_y[k] * smh_313[k];

        t_418[k] = f_3 * pc_y[k] * smh_314[k];

        t_419[k] = f_14 * slh_209[k]
                   + f_1 * smg0_224[k]
                   - f_2 * smg1_224[k]
                   + f_3 * pc_z[k] * smh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, slh_210, slh_315, \
                         slh_318, smg0_225, smg0_228, smg1_225, smg1_228, smh_315, \
                         smh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_14 * slh_315[k]
                   + f_1 * smg0_225[k]
                   - f_2 * smg1_225[k]
                   + f_3 * pc_x[k] * smh_315[k];

        t_421[k] = f_18 * slh_210[k]
                   + f_3 * pc_y[k] * smh_315[k];

        t_422[k] = f_3 * pc_z[k] * smh_315[k];

        t_423[k] = f_14 * slh_318[k]
                   + f_4 * smg0_228[k]
                   - f_5 * smg1_228[k]
                   + f_3 * pc_x[k] * smh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, slh_212, slh_320, slh_321, smg0_230, \
                         smg0_231, smg1_230, smg1_231, smh_317, smh_320, \
                         smh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_18 * slh_212[k]
                   + f_3 * pc_y[k] * smh_317[k];

        t_425[k] = f_14 * slh_320[k]
                   + f_4 * smg0_230[k]
                   - f_5 * smg1_230[k]
                   + f_3 * pc_x[k] * smh_320[k];

        t_426[k] = f_14 * slh_321[k]
                   + f_6 * smg0_231[k]
                   - f_7 * smg1_231[k]
                   + f_3 * pc_x[k] * smh_321[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_x, pc_y, pc_z, slh_215, slh_324, smg0_234, \
                         smg1_234, smh_318, smh_320, smh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_3 * pc_z[k] * smh_318[k];

        t_428[k] = f_18 * slh_215[k]
                   + f_3 * pc_y[k] * smh_320[k];

        t_429[k] = f_14 * slh_324[k]
                   + f_6 * smg0_234[k]
                   - f_7 * smg1_234[k]
                   + f_3 * pc_x[k] * smh_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_z, slh_325, slh_327, smg0_235, \
                         smg0_237, smg1_235, smg1_237, smh_321, smh_325, \
                         smh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_14 * slh_325[k]
                   + f_8 * smg0_235[k]
                   - f_9 * smg1_235[k]
                   + f_3 * pc_x[k] * smh_325[k];

        t_431[k] = f_3 * pc_z[k] * smh_321[k];

        t_432[k] = f_14 * slh_327[k]
                   + f_8 * smg0_237[k]
                   - f_9 * smg1_237[k]
                   + f_3 * pc_x[k] * smh_327[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, slh_219, slh_329, slh_330, \
                         slh_331, smg0_239, smg1_239, smh_324, smh_329, smh_330, \
                         smh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_18 * slh_219[k]
                   + f_3 * pc_y[k] * smh_324[k];

        t_434[k] = f_14 * slh_329[k]
                   + f_8 * smg0_239[k]
                   - f_9 * smg1_239[k]
                   + f_3 * pc_x[k] * smh_329[k];

        t_435[k] = f_14 * slh_330[k]
                   + f_3 * pc_x[k] * smh_330[k];

        t_436[k] = f_14 * slh_331[k]
                   + f_3 * pc_x[k] * smh_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, slh_332, slh_333, slh_334, slh_335, \
                         smh_332, smh_333, smh_334, smh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_14 * slh_332[k]
                   + f_3 * pc_x[k] * smh_332[k];

        t_438[k] = f_14 * slh_333[k]
                   + f_3 * pc_x[k] * smh_333[k];

        t_439[k] = f_14 * slh_334[k]
                   + f_3 * pc_x[k] * smh_334[k];

        t_440[k] = f_14 * slh_335[k]
                   + f_3 * pc_x[k] * smh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, pc_z, slh_225, slh_227, smg0_235, \
                         smg0_237, smg1_235, smg1_237, smh_330, \
                         smh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_18 * slh_225[k]
                   + f_1 * smg0_235[k]
                   - f_2 * smg1_235[k]
                   + f_3 * pc_y[k] * smh_330[k];

        t_442[k] = f_3 * pc_z[k] * smh_330[k];

        t_443[k] = f_18 * slh_227[k]
                   + f_4 * smg0_237[k]
                   - f_5 * smg1_237[k]
                   + f_3 * pc_y[k] * smh_332[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, slh_228, slh_229, slh_230, \
                         smg0_238, smg0_239, smg1_238, smg1_239, smh_333, smh_334, \
                         smh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_18 * slh_228[k]
                   + f_6 * smg0_238[k]
                   - f_7 * smg1_238[k]
                   + f_3 * pc_y[k] * smh_333[k];

        t_445[k] = f_18 * slh_229[k]
                   + f_8 * smg0_239[k]
                   - f_9 * smg1_239[k]
                   + f_3 * pc_y[k] * smh_334[k];

        t_446[k] = f_18 * slh_230[k]
                   + f_3 * pc_y[k] * smh_335[k];

        t_447[k] = f_1 * smg0_239[k]
                   - f_2 * smg1_239[k]
                   + f_3 * pc_z[k] * smh_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_z, pc_y, pc_z, sli0_280, sli0_283, \
                         slh_210, slh_231, sli1_280, sli1_283, \
                         smh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_z[k] * sli0_280[k]
                   - f_10 * pc_z[k] * sli1_280[k];

        t_449[k] = f_14 * slh_231[k]
                   + f_3 * pc_y[k] * smh_336[k];

        t_450[k] = f_11 * slh_210[k]
                   + f_3 * pc_z[k] * smh_336[k];

        t_451[k] = pb_z[k] * sli0_283[k]
                   - f_10 * pc_z[k] * sli1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_z, pc_x, pc_y, pc_z, sli0_286, slh_233, \
                         slh_341, sli1_286, smg0_245, smg1_245, smh_338, \
                         smh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * slh_233[k]
                   + f_3 * pc_y[k] * smh_338[k];

        t_453[k] = f_14 * slh_341[k]
                   + f_4 * smg0_245[k]
                   - f_5 * smg1_245[k]
                   + f_3 * pc_x[k] * smh_341[k];

        t_454[k] = pb_z[k] * sli0_286[k]
                   - f_10 * pc_z[k] * sli1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_y, pc_z, slh_213, slh_236, slh_345, \
                         smg0_249, smg1_249, smh_339, smh_341, \
                         smh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * slh_213[k]
                   + f_3 * pc_z[k] * smh_339[k];

        t_456[k] = f_14 * slh_236[k]
                   + f_3 * pc_y[k] * smh_341[k];

        t_457[k] = f_14 * slh_345[k]
                   + f_6 * smg0_249[k]
                   - f_7 * smg1_249[k]
                   + f_3 * pc_x[k] * smh_345[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pb_z, pc_y, pc_z, sli0_290, sli0_292, \
                         slh_216, slh_217, slh_240, sli1_290, sli1_292, smh_342, \
                         smh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pb_z[k] * sli0_290[k]
                   - f_10 * pc_z[k] * sli1_290[k];

        t_459[k] = f_11 * slh_216[k]
                   + f_3 * pc_z[k] * smh_342[k];

        t_460[k] = pb_z[k] * sli0_292[k]
                   + f_12 * slh_217[k]
                   - f_10 * pc_z[k] * sli1_292[k];

        t_461[k] = f_14 * slh_240[k]
                   + f_3 * pc_y[k] * smh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, slh_350, slh_351, slh_352, slh_353, \
                         smg0_254, smg1_254, smh_350, smh_351, smh_352, \
                         smh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_14 * slh_350[k]
                   + f_8 * smg0_254[k]
                   - f_9 * smg1_254[k]
                   + f_3 * pc_x[k] * smh_350[k];

        t_463[k] = f_14 * slh_351[k]
                   + f_3 * pc_x[k] * smh_351[k];

        t_464[k] = f_14 * slh_352[k]
                   + f_3 * pc_x[k] * smh_352[k];

        t_465[k] = f_14 * slh_353[k]
                   + f_3 * pc_x[k] * smh_353[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_301 = buffer.data(sli0 + 301);
    const auto *sli0_392 = buffer.data(sli0 + 392);
    const auto *sli0_395 = buffer.data(sli0 + 395);
    const auto *sli0_397 = buffer.data(sli0 + 397);
    const auto *sli0_398 = buffer.data(sli0 + 398);
    const auto *sli0_401 = buffer.data(sli0 + 401);
    const auto *sli0_402 = buffer.data(sli0 + 402);
    const auto *sli0_404 = buffer.data(sli0 + 404);
    const auto *sli0_406 = buffer.data(sli0 + 406);
    const auto *sli0_419 = buffer.data(sli0 + 419);

    const auto *slh_225 = buffer.data(slh + 225);
    const auto *slh_230 = buffer.data(slh + 230);
    const auto *slh_231 = buffer.data(slh + 231);
    const auto *slh_234 = buffer.data(slh + 234);
    const auto *slh_237 = buffer.data(slh + 237);
    const auto *slh_246 = buffer.data(slh + 246);
    const auto *slh_248 = buffer.data(slh + 248);
    const auto *slh_249 = buffer.data(slh + 249);
    const auto *slh_250 = buffer.data(slh + 250);
    const auto *slh_251 = buffer.data(slh + 251);
    const auto *slh_252 = buffer.data(slh + 252);
    const auto *slh_254 = buffer.data(slh + 254);
    const auto *slh_255 = buffer.data(slh + 255);
    const auto *slh_257 = buffer.data(slh + 257);
    const auto *slh_258 = buffer.data(slh + 258);
    const auto *slh_261 = buffer.data(slh + 261);
    const auto *slh_267 = buffer.data(slh + 267);
    const auto *slh_269 = buffer.data(slh + 269);
    const auto *slh_270 = buffer.data(slh + 270);
    const auto *slh_271 = buffer.data(slh + 271);
    const auto *slh_272 = buffer.data(slh + 272);
    const auto *slh_273 = buffer.data(slh + 273);
    const auto *slh_275 = buffer.data(slh + 275);
    const auto *slh_276 = buffer.data(slh + 276);
    const auto *slh_278 = buffer.data(slh + 278);
    const auto *slh_279 = buffer.data(slh + 279);
    const auto *slh_282 = buffer.data(slh + 282);
    const auto *slh_288 = buffer.data(slh + 288);
    const auto *slh_290 = buffer.data(slh + 290);
    const auto *slh_291 = buffer.data(slh + 291);
    const auto *slh_292 = buffer.data(slh + 292);
    const auto *slh_293 = buffer.data(slh + 293);
    const auto *slh_294 = buffer.data(slh + 294);
    const auto *slh_295 = buffer.data(slh + 295);
    const auto *slh_296 = buffer.data(slh + 296);
    const auto *slh_297 = buffer.data(slh + 297);
    const auto *slh_299 = buffer.data(slh + 299);
    const auto *slh_300 = buffer.data(slh + 300);
    const auto *slh_302 = buffer.data(slh + 302);
    const auto *slh_303 = buffer.data(slh + 303);
    const auto *slh_309 = buffer.data(slh + 309);
    const auto *slh_311 = buffer.data(slh + 311);
    const auto *slh_312 = buffer.data(slh + 312);
    const auto *slh_313 = buffer.data(slh + 313);
    const auto *slh_314 = buffer.data(slh + 314);
    const auto *slh_354 = buffer.data(slh + 354);
    const auto *slh_355 = buffer.data(slh + 355);
    const auto *slh_356 = buffer.data(slh + 356);
    const auto *slh_357 = buffer.data(slh + 357);
    const auto *slh_360 = buffer.data(slh + 360);
    const auto *slh_362 = buffer.data(slh + 362);
    const auto *slh_363 = buffer.data(slh + 363);
    const auto *slh_366 = buffer.data(slh + 366);
    const auto *slh_367 = buffer.data(slh + 367);
    const auto *slh_369 = buffer.data(slh + 369);
    const auto *slh_371 = buffer.data(slh + 371);
    const auto *slh_372 = buffer.data(slh + 372);
    const auto *slh_373 = buffer.data(slh + 373);
    const auto *slh_374 = buffer.data(slh + 374);
    const auto *slh_375 = buffer.data(slh + 375);
    const auto *slh_376 = buffer.data(slh + 376);
    const auto *slh_377 = buffer.data(slh + 377);
    const auto *slh_378 = buffer.data(slh + 378);
    const auto *slh_381 = buffer.data(slh + 381);
    const auto *slh_383 = buffer.data(slh + 383);
    const auto *slh_384 = buffer.data(slh + 384);
    const auto *slh_387 = buffer.data(slh + 387);
    const auto *slh_388 = buffer.data(slh + 388);
    const auto *slh_390 = buffer.data(slh + 390);
    const auto *slh_392 = buffer.data(slh + 392);
    const auto *slh_393 = buffer.data(slh + 393);
    const auto *slh_394 = buffer.data(slh + 394);
    const auto *slh_395 = buffer.data(slh + 395);
    const auto *slh_396 = buffer.data(slh + 396);
    const auto *slh_397 = buffer.data(slh + 397);
    const auto *slh_398 = buffer.data(slh + 398);
    const auto *slh_414 = buffer.data(slh + 414);
    const auto *slh_415 = buffer.data(slh + 415);
    const auto *slh_416 = buffer.data(slh + 416);
    const auto *slh_417 = buffer.data(slh + 417);
    const auto *slh_418 = buffer.data(slh + 418);
    const auto *slh_419 = buffer.data(slh + 419);
    const auto *slh_420 = buffer.data(slh + 420);
    const auto *slh_423 = buffer.data(slh + 423);
    const auto *slh_425 = buffer.data(slh + 425);
    const auto *slh_426 = buffer.data(slh + 426);
    const auto *slh_429 = buffer.data(slh + 429);
    const auto *slh_430 = buffer.data(slh + 430);
    const auto *slh_432 = buffer.data(slh + 432);
    const auto *slh_434 = buffer.data(slh + 434);

    const auto *sli1_301 = buffer.data(sli1 + 301);
    const auto *sli1_392 = buffer.data(sli1 + 392);
    const auto *sli1_395 = buffer.data(sli1 + 395);
    const auto *sli1_397 = buffer.data(sli1 + 397);
    const auto *sli1_398 = buffer.data(sli1 + 398);
    const auto *sli1_401 = buffer.data(sli1 + 401);
    const auto *sli1_402 = buffer.data(sli1 + 402);
    const auto *sli1_404 = buffer.data(sli1 + 404);
    const auto *sli1_406 = buffer.data(sli1 + 406);
    const auto *sli1_419 = buffer.data(sli1 + 419);

    const auto *smg0_252 = buffer.data(smg0 + 252);
    const auto *smg0_253 = buffer.data(smg0 + 253);
    const auto *smg0_254 = buffer.data(smg0 + 254);
    const auto *smg0_255 = buffer.data(smg0 + 255);
    const auto *smg0_258 = buffer.data(smg0 + 258);
    const auto *smg0_260 = buffer.data(smg0 + 260);
    const auto *smg0_261 = buffer.data(smg0 + 261);
    const auto *smg0_264 = buffer.data(smg0 + 264);
    const auto *smg0_265 = buffer.data(smg0 + 265);
    const auto *smg0_267 = buffer.data(smg0 + 267);
    const auto *smg0_268 = buffer.data(smg0 + 268);
    const auto *smg0_269 = buffer.data(smg0 + 269);
    const auto *smg0_270 = buffer.data(smg0 + 270);
    const auto *smg0_273 = buffer.data(smg0 + 273);
    const auto *smg0_275 = buffer.data(smg0 + 275);
    const auto *smg0_276 = buffer.data(smg0 + 276);
    const auto *smg0_279 = buffer.data(smg0 + 279);
    const auto *smg0_280 = buffer.data(smg0 + 280);
    const auto *smg0_282 = buffer.data(smg0 + 282);
    const auto *smg0_283 = buffer.data(smg0 + 283);
    const auto *smg0_284 = buffer.data(smg0 + 284);
    const auto *smg0_295 = buffer.data(smg0 + 295);
    const auto *smg0_297 = buffer.data(smg0 + 297);
    const auto *smg0_298 = buffer.data(smg0 + 298);
    const auto *smg0_299 = buffer.data(smg0 + 299);
    const auto *smg0_300 = buffer.data(smg0 + 300);
    const auto *smg0_303 = buffer.data(smg0 + 303);
    const auto *smg0_305 = buffer.data(smg0 + 305);
    const auto *smg0_306 = buffer.data(smg0 + 306);
    const auto *smg0_309 = buffer.data(smg0 + 309);
    const auto *smg0_310 = buffer.data(smg0 + 310);
    const auto *smg0_312 = buffer.data(smg0 + 312);
    const auto *smg0_314 = buffer.data(smg0 + 314);

    const auto *smg1_252 = buffer.data(smg1 + 252);
    const auto *smg1_253 = buffer.data(smg1 + 253);
    const auto *smg1_254 = buffer.data(smg1 + 254);
    const auto *smg1_255 = buffer.data(smg1 + 255);
    const auto *smg1_258 = buffer.data(smg1 + 258);
    const auto *smg1_260 = buffer.data(smg1 + 260);
    const auto *smg1_261 = buffer.data(smg1 + 261);
    const auto *smg1_264 = buffer.data(smg1 + 264);
    const auto *smg1_265 = buffer.data(smg1 + 265);
    const auto *smg1_267 = buffer.data(smg1 + 267);
    const auto *smg1_268 = buffer.data(smg1 + 268);
    const auto *smg1_269 = buffer.data(smg1 + 269);
    const auto *smg1_270 = buffer.data(smg1 + 270);
    const auto *smg1_273 = buffer.data(smg1 + 273);
    const auto *smg1_275 = buffer.data(smg1 + 275);
    const auto *smg1_276 = buffer.data(smg1 + 276);
    const auto *smg1_279 = buffer.data(smg1 + 279);
    const auto *smg1_280 = buffer.data(smg1 + 280);
    const auto *smg1_282 = buffer.data(smg1 + 282);
    const auto *smg1_283 = buffer.data(smg1 + 283);
    const auto *smg1_284 = buffer.data(smg1 + 284);
    const auto *smg1_295 = buffer.data(smg1 + 295);
    const auto *smg1_297 = buffer.data(smg1 + 297);
    const auto *smg1_298 = buffer.data(smg1 + 298);
    const auto *smg1_299 = buffer.data(smg1 + 299);
    const auto *smg1_300 = buffer.data(smg1 + 300);
    const auto *smg1_303 = buffer.data(smg1 + 303);
    const auto *smg1_305 = buffer.data(smg1 + 305);
    const auto *smg1_306 = buffer.data(smg1 + 306);
    const auto *smg1_309 = buffer.data(smg1 + 309);
    const auto *smg1_310 = buffer.data(smg1 + 310);
    const auto *smg1_312 = buffer.data(smg1 + 312);
    const auto *smg1_314 = buffer.data(smg1 + 314);

    const auto *smh_351 = buffer.data(smh + 351);
    const auto *smh_353 = buffer.data(smh + 353);
    const auto *smh_354 = buffer.data(smh + 354);
    const auto *smh_355 = buffer.data(smh + 355);
    const auto *smh_356 = buffer.data(smh + 356);
    const auto *smh_357 = buffer.data(smh + 357);
    const auto *smh_359 = buffer.data(smh + 359);
    const auto *smh_360 = buffer.data(smh + 360);
    const auto *smh_362 = buffer.data(smh + 362);
    const auto *smh_363 = buffer.data(smh + 363);
    const auto *smh_366 = buffer.data(smh + 366);
    const auto *smh_367 = buffer.data(smh + 367);
    const auto *smh_369 = buffer.data(smh + 369);
    const auto *smh_371 = buffer.data(smh + 371);
    const auto *smh_372 = buffer.data(smh + 372);
    const auto *smh_373 = buffer.data(smh + 373);
    const auto *smh_374 = buffer.data(smh + 374);
    const auto *smh_375 = buffer.data(smh + 375);
    const auto *smh_376 = buffer.data(smh + 376);
    const auto *smh_377 = buffer.data(smh + 377);
    const auto *smh_378 = buffer.data(smh + 378);
    const auto *smh_380 = buffer.data(smh + 380);
    const auto *smh_381 = buffer.data(smh + 381);
    const auto *smh_383 = buffer.data(smh + 383);
    const auto *smh_384 = buffer.data(smh + 384);
    const auto *smh_387 = buffer.data(smh + 387);
    const auto *smh_388 = buffer.data(smh + 388);
    const auto *smh_390 = buffer.data(smh + 390);
    const auto *smh_392 = buffer.data(smh + 392);
    const auto *smh_393 = buffer.data(smh + 393);
    const auto *smh_394 = buffer.data(smh + 394);
    const auto *smh_395 = buffer.data(smh + 395);
    const auto *smh_396 = buffer.data(smh + 396);
    const auto *smh_397 = buffer.data(smh + 397);
    const auto *smh_398 = buffer.data(smh + 398);
    const auto *smh_399 = buffer.data(smh + 399);
    const auto *smh_401 = buffer.data(smh + 401);
    const auto *smh_402 = buffer.data(smh + 402);
    const auto *smh_404 = buffer.data(smh + 404);
    const auto *smh_405 = buffer.data(smh + 405);
    const auto *smh_408 = buffer.data(smh + 408);
    const auto *smh_414 = buffer.data(smh + 414);
    const auto *smh_415 = buffer.data(smh + 415);
    const auto *smh_416 = buffer.data(smh + 416);
    const auto *smh_417 = buffer.data(smh + 417);
    const auto *smh_418 = buffer.data(smh + 418);
    const auto *smh_419 = buffer.data(smh + 419);
    const auto *smh_420 = buffer.data(smh + 420);
    const auto *smh_422 = buffer.data(smh + 422);
    const auto *smh_423 = buffer.data(smh + 423);
    const auto *smh_425 = buffer.data(smh + 425);
    const auto *smh_426 = buffer.data(smh + 426);
    const auto *smh_429 = buffer.data(smh + 429);
    const auto *smh_430 = buffer.data(smh + 430);
    const auto *smh_432 = buffer.data(smh + 432);
    const auto *smh_434 = buffer.data(smh + 434);

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_z, pc_x, pc_z, sli0_301, slh_354, \
                         slh_355, slh_356, sli1_301, smh_354, smh_355, \
                         smh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * slh_354[k]
                   + f_3 * pc_x[k] * smh_354[k];

        t_467[k] = f_14 * slh_355[k]
                   + f_3 * pc_x[k] * smh_355[k];

        t_468[k] = f_14 * slh_356[k]
                   + f_3 * pc_x[k] * smh_356[k];

        t_469[k] = pb_z[k] * sli0_301[k]
                   - f_10 * pc_z[k] * sli1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, slh_225, slh_248, slh_249, smg0_252, \
                         smg0_253, smg1_252, smg1_253, smh_351, smh_353, \
                         smh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * slh_225[k]
                   + f_3 * pc_z[k] * smh_351[k];

        t_471[k] = f_14 * slh_248[k]
                   + f_4 * smg0_252[k]
                   - f_5 * smg1_252[k]
                   + f_3 * pc_y[k] * smh_353[k];

        t_472[k] = f_14 * slh_249[k]
                   + f_6 * smg0_253[k]
                   - f_7 * smg1_253[k]
                   + f_3 * pc_y[k] * smh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, slh_230, slh_250, slh_251, smg0_254, \
                         smg1_254, smh_355, smh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * slh_250[k]
                   + f_8 * smg0_254[k]
                   - f_9 * smg1_254[k]
                   + f_3 * pc_y[k] * smh_355[k];

        t_474[k] = f_14 * slh_251[k]
                   + f_3 * pc_y[k] * smh_356[k];

        t_475[k] = f_11 * slh_230[k]
                   + f_1 * smg0_254[k]
                   - f_2 * smg1_254[k]
                   + f_3 * pc_z[k] * smh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, slh_231, slh_252, slh_357, \
                         smg0_255, smg1_255, smh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_14 * slh_357[k]
                   + f_1 * smg0_255[k]
                   - f_2 * smg1_255[k]
                   + f_3 * pc_x[k] * smh_357[k];

        t_477[k] = f_13 * slh_252[k]
                   + f_3 * pc_y[k] * smh_357[k];

        t_478[k] = f_12 * slh_231[k]
                   + f_3 * pc_z[k] * smh_357[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, slh_254, slh_360, slh_362, smg0_258, \
                         smg0_260, smg1_258, smg1_260, smh_359, smh_360, \
                         smh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_14 * slh_360[k]
                   + f_4 * smg0_258[k]
                   - f_5 * smg1_258[k]
                   + f_3 * pc_x[k] * smh_360[k];

        t_480[k] = f_13 * slh_254[k]
                   + f_3 * pc_y[k] * smh_359[k];

        t_481[k] = f_14 * slh_362[k]
                   + f_4 * smg0_260[k]
                   - f_5 * smg1_260[k]
                   + f_3 * pc_x[k] * smh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, slh_234, slh_257, slh_363, \
                         smg0_261, smg1_261, smh_360, smh_362, \
                         smh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_14 * slh_363[k]
                   + f_6 * smg0_261[k]
                   - f_7 * smg1_261[k]
                   + f_3 * pc_x[k] * smh_363[k];

        t_483[k] = f_12 * slh_234[k]
                   + f_3 * pc_z[k] * smh_360[k];

        t_484[k] = f_13 * slh_257[k]
                   + f_3 * pc_y[k] * smh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, slh_237, slh_366, slh_367, smg0_264, \
                         smg0_265, smg1_264, smg1_265, smh_363, smh_366, \
                         smh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_14 * slh_366[k]
                   + f_6 * smg0_264[k]
                   - f_7 * smg1_264[k]
                   + f_3 * pc_x[k] * smh_366[k];

        t_486[k] = f_14 * slh_367[k]
                   + f_8 * smg0_265[k]
                   - f_9 * smg1_265[k]
                   + f_3 * pc_x[k] * smh_367[k];

        t_487[k] = f_12 * slh_237[k]
                   + f_3 * pc_z[k] * smh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, slh_261, slh_369, slh_371, smg0_267, \
                         smg0_269, smg1_267, smg1_269, smh_366, smh_369, \
                         smh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_14 * slh_369[k]
                   + f_8 * smg0_267[k]
                   - f_9 * smg1_267[k]
                   + f_3 * pc_x[k] * smh_369[k];

        t_489[k] = f_13 * slh_261[k]
                   + f_3 * pc_y[k] * smh_366[k];

        t_490[k] = f_14 * slh_371[k]
                   + f_8 * smg0_269[k]
                   - f_9 * smg1_269[k]
                   + f_3 * pc_x[k] * smh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, slh_372, slh_373, slh_374, \
                         slh_375, slh_376, smh_372, smh_373, smh_374, smh_375, \
                         smh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_14 * slh_372[k]
                   + f_3 * pc_x[k] * smh_372[k];

        t_492[k] = f_14 * slh_373[k]
                   + f_3 * pc_x[k] * smh_373[k];

        t_493[k] = f_14 * slh_374[k]
                   + f_3 * pc_x[k] * smh_374[k];

        t_494[k] = f_14 * slh_375[k]
                   + f_3 * pc_x[k] * smh_375[k];

        t_495[k] = f_14 * slh_376[k]
                   + f_3 * pc_x[k] * smh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, slh_246, slh_267, slh_377, \
                         smg0_265, smg1_265, smh_372, smh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_14 * slh_377[k]
                   + f_3 * pc_x[k] * smh_377[k];

        t_497[k] = f_13 * slh_267[k]
                   + f_1 * smg0_265[k]
                   - f_2 * smg1_265[k]
                   + f_3 * pc_y[k] * smh_372[k];

        t_498[k] = f_12 * slh_246[k]
                   + f_3 * pc_z[k] * smh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, slh_269, slh_270, slh_271, smg0_267, \
                         smg0_268, smg0_269, smg1_267, smg1_268, smg1_269, smh_374, smh_375, \
                         smh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * slh_269[k]
                   + f_4 * smg0_267[k]
                   - f_5 * smg1_267[k]
                   + f_3 * pc_y[k] * smh_374[k];

        t_500[k] = f_13 * slh_270[k]
                   + f_6 * smg0_268[k]
                   - f_7 * smg1_268[k]
                   + f_3 * pc_y[k] * smh_375[k];

        t_501[k] = f_13 * slh_271[k]
                   + f_8 * smg0_269[k]
                   - f_9 * smg1_269[k]
                   + f_3 * pc_y[k] * smh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, slh_251, slh_272, slh_378, \
                         smg0_269, smg0_270, smg1_269, smg1_270, smh_377, \
                         smh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * slh_272[k]
                   + f_3 * pc_y[k] * smh_377[k];

        t_503[k] = f_12 * slh_251[k]
                   + f_1 * smg0_269[k]
                   - f_2 * smg1_269[k]
                   + f_3 * pc_z[k] * smh_377[k];

        t_504[k] = f_14 * slh_378[k]
                   + f_1 * smg0_270[k]
                   - f_2 * smg1_270[k]
                   + f_3 * pc_x[k] * smh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, slh_252, slh_273, \
                         slh_275, slh_381, smg0_273, smg1_273, smh_378, smh_380, \
                         smh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * slh_273[k]
                   + f_3 * pc_y[k] * smh_378[k];

        t_506[k] = f_13 * slh_252[k]
                   + f_3 * pc_z[k] * smh_378[k];

        t_507[k] = f_14 * slh_381[k]
                   + f_4 * smg0_273[k]
                   - f_5 * smg1_273[k]
                   + f_3 * pc_x[k] * smh_381[k];

        t_508[k] = f_12 * slh_275[k]
                   + f_3 * pc_y[k] * smh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, slh_255, slh_383, slh_384, smg0_275, \
                         smg0_276, smg1_275, smg1_276, smh_381, smh_383, \
                         smh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_14 * slh_383[k]
                   + f_4 * smg0_275[k]
                   - f_5 * smg1_275[k]
                   + f_3 * pc_x[k] * smh_383[k];

        t_510[k] = f_14 * slh_384[k]
                   + f_6 * smg0_276[k]
                   - f_7 * smg1_276[k]
                   + f_3 * pc_x[k] * smh_384[k];

        t_511[k] = f_13 * slh_255[k]
                   + f_3 * pc_z[k] * smh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, slh_278, slh_387, slh_388, smg0_279, \
                         smg0_280, smg1_279, smg1_280, smh_383, smh_387, \
                         smh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * slh_278[k]
                   + f_3 * pc_y[k] * smh_383[k];

        t_513[k] = f_14 * slh_387[k]
                   + f_6 * smg0_279[k]
                   - f_7 * smg1_279[k]
                   + f_3 * pc_x[k] * smh_387[k];

        t_514[k] = f_14 * slh_388[k]
                   + f_8 * smg0_280[k]
                   - f_9 * smg1_280[k]
                   + f_3 * pc_x[k] * smh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, slh_258, slh_282, slh_390, \
                         smg0_282, smg1_282, smh_384, smh_387, \
                         smh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * slh_258[k]
                   + f_3 * pc_z[k] * smh_384[k];

        t_516[k] = f_14 * slh_390[k]
                   + f_8 * smg0_282[k]
                   - f_9 * smg1_282[k]
                   + f_3 * pc_x[k] * smh_390[k];

        t_517[k] = f_12 * slh_282[k]
                   + f_3 * pc_y[k] * smh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, slh_392, slh_393, slh_394, slh_395, \
                         smg0_284, smg1_284, smh_392, smh_393, smh_394, \
                         smh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_14 * slh_392[k]
                   + f_8 * smg0_284[k]
                   - f_9 * smg1_284[k]
                   + f_3 * pc_x[k] * smh_392[k];

        t_519[k] = f_14 * slh_393[k]
                   + f_3 * pc_x[k] * smh_393[k];

        t_520[k] = f_14 * slh_394[k]
                   + f_3 * pc_x[k] * smh_394[k];

        t_521[k] = f_14 * slh_395[k]
                   + f_3 * pc_x[k] * smh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, slh_288, slh_396, slh_397, \
                         slh_398, smg0_280, smg1_280, smh_393, smh_396, smh_397, \
                         smh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_14 * slh_396[k]
                   + f_3 * pc_x[k] * smh_396[k];

        t_523[k] = f_14 * slh_397[k]
                   + f_3 * pc_x[k] * smh_397[k];

        t_524[k] = f_14 * slh_398[k]
                   + f_3 * pc_x[k] * smh_398[k];

        t_525[k] = f_12 * slh_288[k]
                   + f_1 * smg0_280[k]
                   - f_2 * smg1_280[k]
                   + f_3 * pc_y[k] * smh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, slh_267, slh_290, slh_291, smg0_282, \
                         smg0_283, smg1_282, smg1_283, smh_393, smh_395, \
                         smh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * slh_267[k]
                   + f_3 * pc_z[k] * smh_393[k];

        t_527[k] = f_12 * slh_290[k]
                   + f_4 * smg0_282[k]
                   - f_5 * smg1_282[k]
                   + f_3 * pc_y[k] * smh_395[k];

        t_528[k] = f_12 * slh_291[k]
                   + f_6 * smg0_283[k]
                   - f_7 * smg1_283[k]
                   + f_3 * pc_y[k] * smh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pb_y, pc_y, pc_z, sli0_392, slh_272, \
                         slh_292, slh_293, sli1_392, smg0_284, smg1_284, smh_397, \
                         smh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * slh_292[k]
                   + f_8 * smg0_284[k]
                   - f_9 * smg1_284[k]
                   + f_3 * pc_y[k] * smh_397[k];

        t_530[k] = f_12 * slh_293[k]
                   + f_3 * pc_y[k] * smh_398[k];

        t_531[k] = f_13 * slh_272[k]
                   + f_1 * smg0_284[k]
                   - f_2 * smg1_284[k]
                   + f_3 * pc_z[k] * smh_398[k];

        t_532[k] = pb_y[k] * sli0_392[k]
                   - f_10 * pc_y[k] * sli1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pb_y, pc_y, pc_z, sli0_395, slh_273, \
                         slh_294, slh_295, slh_296, sli1_395, smh_399, \
                         smh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * slh_294[k]
                   + f_3 * pc_y[k] * smh_399[k];

        t_534[k] = f_14 * slh_273[k]
                   + f_3 * pc_z[k] * smh_399[k];

        t_535[k] = pb_y[k] * sli0_395[k]
                   + f_12 * slh_295[k]
                   - f_10 * pc_y[k] * sli1_395[k];

        t_536[k] = f_11 * slh_296[k]
                   + f_3 * pc_y[k] * smh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pb_y, pc_y, pc_z, sli0_397, sli0_398, \
                         slh_276, slh_297, slh_299, sli1_397, sli1_398, smh_402, \
                         smh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pb_y[k] * sli0_397[k]
                   - f_10 * pc_y[k] * sli1_397[k];

        t_538[k] = pb_y[k] * sli0_398[k]
                   + f_13 * slh_297[k]
                   - f_10 * pc_y[k] * sli1_398[k];

        t_539[k] = f_14 * slh_276[k]
                   + f_3 * pc_z[k] * smh_402[k];

        t_540[k] = f_11 * slh_299[k]
                   + f_3 * pc_y[k] * smh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pb_y, pc_y, pc_z, sli0_401, sli0_402, slh_279, \
                         slh_300, sli1_401, sli1_402, smh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pb_y[k] * sli0_401[k]
                   - f_10 * pc_y[k] * sli1_401[k];

        t_542[k] = pb_y[k] * sli0_402[k]
                   + f_14 * slh_300[k]
                   - f_10 * pc_y[k] * sli1_402[k];

        t_543[k] = f_14 * slh_279[k]
                   + f_3 * pc_z[k] * smh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_x, pc_y, sli0_404, sli0_406, \
                         slh_302, slh_303, slh_414, sli1_404, sli1_406, smh_408, \
                         smh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pb_y[k] * sli0_404[k]
                   + f_12 * slh_302[k]
                   - f_10 * pc_y[k] * sli1_404[k];

        t_545[k] = f_11 * slh_303[k]
                   + f_3 * pc_y[k] * smh_408[k];

        t_546[k] = pb_y[k] * sli0_406[k]
                   - f_10 * pc_y[k] * sli1_406[k];

        t_547[k] = f_14 * slh_414[k]
                   + f_3 * pc_x[k] * smh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, slh_415, slh_416, slh_417, \
                         slh_418, slh_419, smh_415, smh_416, smh_417, smh_418, \
                         smh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_14 * slh_415[k]
                   + f_3 * pc_x[k] * smh_415[k];

        t_549[k] = f_14 * slh_416[k]
                   + f_3 * pc_x[k] * smh_416[k];

        t_550[k] = f_14 * slh_417[k]
                   + f_3 * pc_x[k] * smh_417[k];

        t_551[k] = f_14 * slh_418[k]
                   + f_3 * pc_x[k] * smh_418[k];

        t_552[k] = f_14 * slh_419[k]
                   + f_3 * pc_x[k] * smh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, slh_288, slh_309, slh_311, smg0_295, \
                         smg0_297, smg1_295, smg1_297, smh_414, \
                         smh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * slh_309[k]
                   + f_1 * smg0_295[k]
                   - f_2 * smg1_295[k]
                   + f_3 * pc_y[k] * smh_414[k];

        t_554[k] = f_14 * slh_288[k]
                   + f_3 * pc_z[k] * smh_414[k];

        t_555[k] = f_11 * slh_311[k]
                   + f_4 * smg0_297[k]
                   - f_5 * smg1_297[k]
                   + f_3 * pc_y[k] * smh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, slh_312, slh_313, slh_314, smg0_298, \
                         smg0_299, smg1_298, smg1_299, smh_417, smh_418, \
                         smh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * slh_312[k]
                   + f_6 * smg0_298[k]
                   - f_7 * smg1_298[k]
                   + f_3 * pc_y[k] * smh_417[k];

        t_557[k] = f_11 * slh_313[k]
                   + f_8 * smg0_299[k]
                   - f_9 * smg1_299[k]
                   + f_3 * pc_y[k] * smh_418[k];

        t_558[k] = f_11 * slh_314[k]
                   + f_3 * pc_y[k] * smh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_y, pc_x, pc_y, pc_z, sli0_419, \
                         slh_294, slh_420, sli1_419, smg0_300, smg1_300, \
                         smh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pb_y[k] * sli0_419[k]
                   - f_10 * pc_y[k] * sli1_419[k];

        t_560[k] = f_14 * slh_420[k]
                   + f_1 * smg0_300[k]
                   - f_2 * smg1_300[k]
                   + f_3 * pc_x[k] * smh_420[k];

        t_561[k] = f_3 * pc_y[k] * smh_420[k];

        t_562[k] = f_18 * slh_294[k]
                   + f_3 * pc_z[k] * smh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, slh_423, slh_425, smg0_303, \
                         smg0_305, smg1_303, smg1_305, smh_422, smh_423, \
                         smh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_14 * slh_423[k]
                   + f_4 * smg0_303[k]
                   - f_5 * smg1_303[k]
                   + f_3 * pc_x[k] * smh_423[k];

        t_564[k] = f_3 * pc_y[k] * smh_422[k];

        t_565[k] = f_14 * slh_425[k]
                   + f_4 * smg0_305[k]
                   - f_5 * smg1_305[k]
                   + f_3 * pc_x[k] * smh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_x, pc_y, pc_z, slh_297, slh_426, smg0_306, \
                         smg1_306, smh_423, smh_425, smh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_14 * slh_426[k]
                   + f_6 * smg0_306[k]
                   - f_7 * smg1_306[k]
                   + f_3 * pc_x[k] * smh_426[k];

        t_567[k] = f_18 * slh_297[k]
                   + f_3 * pc_z[k] * smh_423[k];

        t_568[k] = f_3 * pc_y[k] * smh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_z, slh_300, slh_429, slh_430, smg0_309, \
                         smg0_310, smg1_309, smg1_310, smh_426, smh_429, \
                         smh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_14 * slh_429[k]
                   + f_6 * smg0_309[k]
                   - f_7 * smg1_309[k]
                   + f_3 * pc_x[k] * smh_429[k];

        t_570[k] = f_14 * slh_430[k]
                   + f_8 * smg0_310[k]
                   - f_9 * smg1_310[k]
                   + f_3 * pc_x[k] * smh_430[k];

        t_571[k] = f_18 * slh_300[k]
                   + f_3 * pc_z[k] * smh_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pc_x, pc_y, slh_432, slh_434, smg0_312, \
                         smg0_314, smg1_312, smg1_314, smh_429, smh_432, \
                         smh_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_14 * slh_432[k]
                   + f_8 * smg0_312[k]
                   - f_9 * smg1_312[k]
                   + f_3 * pc_x[k] * smh_432[k];

        t_573[k] = f_3 * pc_y[k] * smh_429[k];

        t_574[k] = f_14 * slh_434[k]
                   + f_8 * smg0_314[k]
                   - f_9 * smg1_314[k]
                   + f_3 * pc_x[k] * smh_434[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_420 = buffer.data(sli0 + 420);
    const auto *sli0_423 = buffer.data(sli0 + 423);
    const auto *sli0_426 = buffer.data(sli0 + 426);
    const auto *sli0_430 = buffer.data(sli0 + 430);
    const auto *sli0_432 = buffer.data(sli0 + 432);
    const auto *sli0_441 = buffer.data(sli0 + 441);

    const auto *slh_309 = buffer.data(slh + 309);
    const auto *slh_314 = buffer.data(slh + 314);
    const auto *slh_315 = buffer.data(slh + 315);
    const auto *slh_317 = buffer.data(slh + 317);
    const auto *slh_318 = buffer.data(slh + 318);
    const auto *slh_320 = buffer.data(slh + 320);
    const auto *slh_321 = buffer.data(slh + 321);
    const auto *slh_322 = buffer.data(slh + 322);
    const auto *slh_324 = buffer.data(slh + 324);
    const auto *slh_330 = buffer.data(slh + 330);
    const auto *slh_332 = buffer.data(slh + 332);
    const auto *slh_333 = buffer.data(slh + 333);
    const auto *slh_334 = buffer.data(slh + 334);
    const auto *slh_335 = buffer.data(slh + 335);
    const auto *slh_336 = buffer.data(slh + 336);
    const auto *slh_338 = buffer.data(slh + 338);
    const auto *slh_339 = buffer.data(slh + 339);
    const auto *slh_341 = buffer.data(slh + 341);
    const auto *slh_342 = buffer.data(slh + 342);
    const auto *slh_345 = buffer.data(slh + 345);
    const auto *slh_351 = buffer.data(slh + 351);
    const auto *slh_353 = buffer.data(slh + 353);
    const auto *slh_354 = buffer.data(slh + 354);
    const auto *slh_355 = buffer.data(slh + 355);
    const auto *slh_356 = buffer.data(slh + 356);
    const auto *slh_357 = buffer.data(slh + 357);
    const auto *slh_359 = buffer.data(slh + 359);
    const auto *slh_360 = buffer.data(slh + 360);
    const auto *slh_362 = buffer.data(slh + 362);
    const auto *slh_363 = buffer.data(slh + 363);
    const auto *slh_366 = buffer.data(slh + 366);
    const auto *slh_372 = buffer.data(slh + 372);
    const auto *slh_374 = buffer.data(slh + 374);
    const auto *slh_375 = buffer.data(slh + 375);
    const auto *slh_376 = buffer.data(slh + 376);
    const auto *slh_377 = buffer.data(slh + 377);
    const auto *slh_378 = buffer.data(slh + 378);
    const auto *slh_380 = buffer.data(slh + 380);
    const auto *slh_383 = buffer.data(slh + 383);
    const auto *slh_387 = buffer.data(slh + 387);
    const auto *slh_435 = buffer.data(slh + 435);
    const auto *slh_436 = buffer.data(slh + 436);
    const auto *slh_437 = buffer.data(slh + 437);
    const auto *slh_438 = buffer.data(slh + 438);
    const auto *slh_439 = buffer.data(slh + 439);
    const auto *slh_440 = buffer.data(slh + 440);
    const auto *slh_441 = buffer.data(slh + 441);
    const auto *slh_444 = buffer.data(slh + 444);
    const auto *slh_446 = buffer.data(slh + 446);
    const auto *slh_447 = buffer.data(slh + 447);
    const auto *slh_450 = buffer.data(slh + 450);
    const auto *slh_451 = buffer.data(slh + 451);
    const auto *slh_453 = buffer.data(slh + 453);
    const auto *slh_455 = buffer.data(slh + 455);
    const auto *slh_456 = buffer.data(slh + 456);
    const auto *slh_457 = buffer.data(slh + 457);
    const auto *slh_458 = buffer.data(slh + 458);
    const auto *slh_459 = buffer.data(slh + 459);
    const auto *slh_460 = buffer.data(slh + 460);
    const auto *slh_461 = buffer.data(slh + 461);
    const auto *slh_467 = buffer.data(slh + 467);
    const auto *slh_471 = buffer.data(slh + 471);
    const auto *slh_476 = buffer.data(slh + 476);
    const auto *slh_477 = buffer.data(slh + 477);
    const auto *slh_478 = buffer.data(slh + 478);
    const auto *slh_479 = buffer.data(slh + 479);
    const auto *slh_480 = buffer.data(slh + 480);
    const auto *slh_481 = buffer.data(slh + 481);
    const auto *slh_482 = buffer.data(slh + 482);
    const auto *slh_483 = buffer.data(slh + 483);
    const auto *slh_486 = buffer.data(slh + 486);
    const auto *slh_488 = buffer.data(slh + 488);
    const auto *slh_489 = buffer.data(slh + 489);
    const auto *slh_492 = buffer.data(slh + 492);
    const auto *slh_493 = buffer.data(slh + 493);
    const auto *slh_495 = buffer.data(slh + 495);
    const auto *slh_497 = buffer.data(slh + 497);
    const auto *slh_498 = buffer.data(slh + 498);
    const auto *slh_499 = buffer.data(slh + 499);
    const auto *slh_500 = buffer.data(slh + 500);
    const auto *slh_501 = buffer.data(slh + 501);
    const auto *slh_502 = buffer.data(slh + 502);
    const auto *slh_503 = buffer.data(slh + 503);
    const auto *slh_504 = buffer.data(slh + 504);
    const auto *slh_507 = buffer.data(slh + 507);
    const auto *slh_509 = buffer.data(slh + 509);
    const auto *slh_510 = buffer.data(slh + 510);
    const auto *slh_513 = buffer.data(slh + 513);
    const auto *slh_514 = buffer.data(slh + 514);
    const auto *slh_516 = buffer.data(slh + 516);

    const auto *sli1_420 = buffer.data(sli1 + 420);
    const auto *sli1_423 = buffer.data(sli1 + 423);
    const auto *sli1_426 = buffer.data(sli1 + 426);
    const auto *sli1_430 = buffer.data(sli1 + 430);
    const auto *sli1_432 = buffer.data(sli1 + 432);
    const auto *sli1_441 = buffer.data(sli1 + 441);

    const auto *smg0_310 = buffer.data(smg0 + 310);
    const auto *smg0_312 = buffer.data(smg0 + 312);
    const auto *smg0_313 = buffer.data(smg0 + 313);
    const auto *smg0_314 = buffer.data(smg0 + 314);
    const auto *smg0_315 = buffer.data(smg0 + 315);
    const auto *smg0_318 = buffer.data(smg0 + 318);
    const auto *smg0_320 = buffer.data(smg0 + 320);
    const auto *smg0_321 = buffer.data(smg0 + 321);
    const auto *smg0_324 = buffer.data(smg0 + 324);
    const auto *smg0_325 = buffer.data(smg0 + 325);
    const auto *smg0_327 = buffer.data(smg0 + 327);
    const auto *smg0_328 = buffer.data(smg0 + 328);
    const auto *smg0_329 = buffer.data(smg0 + 329);
    const auto *smg0_335 = buffer.data(smg0 + 335);
    const auto *smg0_339 = buffer.data(smg0 + 339);
    const auto *smg0_342 = buffer.data(smg0 + 342);
    const auto *smg0_343 = buffer.data(smg0 + 343);
    const auto *smg0_344 = buffer.data(smg0 + 344);
    const auto *smg0_345 = buffer.data(smg0 + 345);
    const auto *smg0_348 = buffer.data(smg0 + 348);
    const auto *smg0_350 = buffer.data(smg0 + 350);
    const auto *smg0_351 = buffer.data(smg0 + 351);
    const auto *smg0_354 = buffer.data(smg0 + 354);
    const auto *smg0_355 = buffer.data(smg0 + 355);
    const auto *smg0_357 = buffer.data(smg0 + 357);
    const auto *smg0_358 = buffer.data(smg0 + 358);
    const auto *smg0_359 = buffer.data(smg0 + 359);
    const auto *smg0_360 = buffer.data(smg0 + 360);
    const auto *smg0_363 = buffer.data(smg0 + 363);
    const auto *smg0_365 = buffer.data(smg0 + 365);
    const auto *smg0_366 = buffer.data(smg0 + 366);
    const auto *smg0_369 = buffer.data(smg0 + 369);
    const auto *smg0_370 = buffer.data(smg0 + 370);
    const auto *smg0_372 = buffer.data(smg0 + 372);

    const auto *smg1_310 = buffer.data(smg1 + 310);
    const auto *smg1_312 = buffer.data(smg1 + 312);
    const auto *smg1_313 = buffer.data(smg1 + 313);
    const auto *smg1_314 = buffer.data(smg1 + 314);
    const auto *smg1_315 = buffer.data(smg1 + 315);
    const auto *smg1_318 = buffer.data(smg1 + 318);
    const auto *smg1_320 = buffer.data(smg1 + 320);
    const auto *smg1_321 = buffer.data(smg1 + 321);
    const auto *smg1_324 = buffer.data(smg1 + 324);
    const auto *smg1_325 = buffer.data(smg1 + 325);
    const auto *smg1_327 = buffer.data(smg1 + 327);
    const auto *smg1_328 = buffer.data(smg1 + 328);
    const auto *smg1_329 = buffer.data(smg1 + 329);
    const auto *smg1_335 = buffer.data(smg1 + 335);
    const auto *smg1_339 = buffer.data(smg1 + 339);
    const auto *smg1_342 = buffer.data(smg1 + 342);
    const auto *smg1_343 = buffer.data(smg1 + 343);
    const auto *smg1_344 = buffer.data(smg1 + 344);
    const auto *smg1_345 = buffer.data(smg1 + 345);
    const auto *smg1_348 = buffer.data(smg1 + 348);
    const auto *smg1_350 = buffer.data(smg1 + 350);
    const auto *smg1_351 = buffer.data(smg1 + 351);
    const auto *smg1_354 = buffer.data(smg1 + 354);
    const auto *smg1_355 = buffer.data(smg1 + 355);
    const auto *smg1_357 = buffer.data(smg1 + 357);
    const auto *smg1_358 = buffer.data(smg1 + 358);
    const auto *smg1_359 = buffer.data(smg1 + 359);
    const auto *smg1_360 = buffer.data(smg1 + 360);
    const auto *smg1_363 = buffer.data(smg1 + 363);
    const auto *smg1_365 = buffer.data(smg1 + 365);
    const auto *smg1_366 = buffer.data(smg1 + 366);
    const auto *smg1_369 = buffer.data(smg1 + 369);
    const auto *smg1_370 = buffer.data(smg1 + 370);
    const auto *smg1_372 = buffer.data(smg1 + 372);

    const auto *smh_435 = buffer.data(smh + 435);
    const auto *smh_436 = buffer.data(smh + 436);
    const auto *smh_437 = buffer.data(smh + 437);
    const auto *smh_438 = buffer.data(smh + 438);
    const auto *smh_439 = buffer.data(smh + 439);
    const auto *smh_440 = buffer.data(smh + 440);
    const auto *smh_441 = buffer.data(smh + 441);
    const auto *smh_443 = buffer.data(smh + 443);
    const auto *smh_444 = buffer.data(smh + 444);
    const auto *smh_446 = buffer.data(smh + 446);
    const auto *smh_447 = buffer.data(smh + 447);
    const auto *smh_450 = buffer.data(smh + 450);
    const auto *smh_451 = buffer.data(smh + 451);
    const auto *smh_453 = buffer.data(smh + 453);
    const auto *smh_455 = buffer.data(smh + 455);
    const auto *smh_456 = buffer.data(smh + 456);
    const auto *smh_457 = buffer.data(smh + 457);
    const auto *smh_458 = buffer.data(smh + 458);
    const auto *smh_459 = buffer.data(smh + 459);
    const auto *smh_460 = buffer.data(smh + 460);
    const auto *smh_461 = buffer.data(smh + 461);
    const auto *smh_462 = buffer.data(smh + 462);
    const auto *smh_464 = buffer.data(smh + 464);
    const auto *smh_465 = buffer.data(smh + 465);
    const auto *smh_467 = buffer.data(smh + 467);
    const auto *smh_468 = buffer.data(smh + 468);
    const auto *smh_471 = buffer.data(smh + 471);
    const auto *smh_476 = buffer.data(smh + 476);
    const auto *smh_477 = buffer.data(smh + 477);
    const auto *smh_478 = buffer.data(smh + 478);
    const auto *smh_479 = buffer.data(smh + 479);
    const auto *smh_480 = buffer.data(smh + 480);
    const auto *smh_481 = buffer.data(smh + 481);
    const auto *smh_482 = buffer.data(smh + 482);
    const auto *smh_483 = buffer.data(smh + 483);
    const auto *smh_485 = buffer.data(smh + 485);
    const auto *smh_486 = buffer.data(smh + 486);
    const auto *smh_488 = buffer.data(smh + 488);
    const auto *smh_489 = buffer.data(smh + 489);
    const auto *smh_492 = buffer.data(smh + 492);
    const auto *smh_493 = buffer.data(smh + 493);
    const auto *smh_495 = buffer.data(smh + 495);
    const auto *smh_497 = buffer.data(smh + 497);
    const auto *smh_498 = buffer.data(smh + 498);
    const auto *smh_499 = buffer.data(smh + 499);
    const auto *smh_500 = buffer.data(smh + 500);
    const auto *smh_501 = buffer.data(smh + 501);
    const auto *smh_502 = buffer.data(smh + 502);
    const auto *smh_503 = buffer.data(smh + 503);
    const auto *smh_504 = buffer.data(smh + 504);
    const auto *smh_506 = buffer.data(smh + 506);
    const auto *smh_507 = buffer.data(smh + 507);
    const auto *smh_509 = buffer.data(smh + 509);
    const auto *smh_510 = buffer.data(smh + 510);
    const auto *smh_513 = buffer.data(smh + 513);
    const auto *smh_514 = buffer.data(smh + 514);
    const auto *smh_516 = buffer.data(smh + 516);

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, pc_x, slh_435, slh_436, slh_437, \
                         slh_438, slh_439, smh_435, smh_436, smh_437, smh_438, \
                         smh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_14 * slh_435[k]
                   + f_3 * pc_x[k] * smh_435[k];

        t_576[k] = f_14 * slh_436[k]
                   + f_3 * pc_x[k] * smh_436[k];

        t_577[k] = f_14 * slh_437[k]
                   + f_3 * pc_x[k] * smh_437[k];

        t_578[k] = f_14 * slh_438[k]
                   + f_3 * pc_x[k] * smh_438[k];

        t_579[k] = f_14 * slh_439[k]
                   + f_3 * pc_x[k] * smh_439[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, pc_z, slh_309, slh_440, \
                         smg0_310, smg0_312, smg1_310, smg1_312, smh_435, smh_437, \
                         smh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_14 * slh_440[k]
                   + f_3 * pc_x[k] * smh_440[k];

        t_581[k] = f_1 * smg0_310[k]
                   - f_2 * smg1_310[k]
                   + f_3 * pc_y[k] * smh_435[k];

        t_582[k] = f_18 * slh_309[k]
                   + f_3 * pc_z[k] * smh_435[k];

        t_583[k] = f_4 * smg0_312[k]
                   - f_5 * smg1_312[k]
                   + f_3 * pc_y[k] * smh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, slh_314, smg0_313, smg0_314, \
                         smg1_313, smg1_314, smh_438, smh_439, \
                         smh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * smg0_313[k]
                   - f_7 * smg1_313[k]
                   + f_3 * pc_y[k] * smh_438[k];

        t_585[k] = f_8 * smg0_314[k]
                   - f_9 * smg1_314[k]
                   + f_3 * pc_y[k] * smh_439[k];

        t_586[k] = f_3 * pc_y[k] * smh_440[k];

        t_587[k] = f_18 * slh_314[k]
                   + f_1 * smg0_314[k]
                   - f_2 * smg1_314[k]
                   + f_3 * pc_z[k] * smh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, slh_315, slh_441, \
                         slh_444, smg0_315, smg0_318, smg1_315, smg1_318, smh_441, \
                         smh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_13 * slh_441[k]
                   + f_1 * smg0_315[k]
                   - f_2 * smg1_315[k]
                   + f_3 * pc_x[k] * smh_441[k];

        t_589[k] = f_17 * slh_315[k]
                   + f_3 * pc_y[k] * smh_441[k];

        t_590[k] = f_3 * pc_z[k] * smh_441[k];

        t_591[k] = f_13 * slh_444[k]
                   + f_4 * smg0_318[k]
                   - f_5 * smg1_318[k]
                   + f_3 * pc_x[k] * smh_444[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, pc_x, pc_y, slh_317, slh_446, slh_447, smg0_320, \
                         smg0_321, smg1_320, smg1_321, smh_443, smh_446, \
                         smh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * slh_317[k]
                   + f_3 * pc_y[k] * smh_443[k];

        t_593[k] = f_13 * slh_446[k]
                   + f_4 * smg0_320[k]
                   - f_5 * smg1_320[k]
                   + f_3 * pc_x[k] * smh_446[k];

        t_594[k] = f_13 * slh_447[k]
                   + f_6 * smg0_321[k]
                   - f_7 * smg1_321[k]
                   + f_3 * pc_x[k] * smh_447[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_x, pc_y, pc_z, slh_320, slh_450, smg0_324, \
                         smg1_324, smh_444, smh_446, smh_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_3 * pc_z[k] * smh_444[k];

        t_596[k] = f_17 * slh_320[k]
                   + f_3 * pc_y[k] * smh_446[k];

        t_597[k] = f_13 * slh_450[k]
                   + f_6 * smg0_324[k]
                   - f_7 * smg1_324[k]
                   + f_3 * pc_x[k] * smh_450[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pc_x, pc_z, slh_451, slh_453, smg0_325, \
                         smg0_327, smg1_325, smg1_327, smh_447, smh_451, \
                         smh_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_13 * slh_451[k]
                   + f_8 * smg0_325[k]
                   - f_9 * smg1_325[k]
                   + f_3 * pc_x[k] * smh_451[k];

        t_599[k] = f_3 * pc_z[k] * smh_447[k];

        t_600[k] = f_13 * slh_453[k]
                   + f_8 * smg0_327[k]
                   - f_9 * smg1_327[k]
                   + f_3 * pc_x[k] * smh_453[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pc_x, pc_y, slh_324, slh_455, slh_456, \
                         slh_457, smg0_329, smg1_329, smh_450, smh_455, smh_456, \
                         smh_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_17 * slh_324[k]
                   + f_3 * pc_y[k] * smh_450[k];

        t_602[k] = f_13 * slh_455[k]
                   + f_8 * smg0_329[k]
                   - f_9 * smg1_329[k]
                   + f_3 * pc_x[k] * smh_455[k];

        t_603[k] = f_13 * slh_456[k]
                   + f_3 * pc_x[k] * smh_456[k];

        t_604[k] = f_13 * slh_457[k]
                   + f_3 * pc_x[k] * smh_457[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pc_x, slh_458, slh_459, slh_460, slh_461, \
                         smh_458, smh_459, smh_460, smh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_13 * slh_458[k]
                   + f_3 * pc_x[k] * smh_458[k];

        t_606[k] = f_13 * slh_459[k]
                   + f_3 * pc_x[k] * smh_459[k];

        t_607[k] = f_13 * slh_460[k]
                   + f_3 * pc_x[k] * smh_460[k];

        t_608[k] = f_13 * slh_461[k]
                   + f_3 * pc_x[k] * smh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_y, pc_z, slh_330, slh_332, smg0_325, \
                         smg0_327, smg1_325, smg1_327, smh_456, \
                         smh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_17 * slh_330[k]
                   + f_1 * smg0_325[k]
                   - f_2 * smg1_325[k]
                   + f_3 * pc_y[k] * smh_456[k];

        t_610[k] = f_3 * pc_z[k] * smh_456[k];

        t_611[k] = f_17 * slh_332[k]
                   + f_4 * smg0_327[k]
                   - f_5 * smg1_327[k]
                   + f_3 * pc_y[k] * smh_458[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pc_y, pc_z, slh_333, slh_334, slh_335, \
                         smg0_328, smg0_329, smg1_328, smg1_329, smh_459, smh_460, \
                         smh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_17 * slh_333[k]
                   + f_6 * smg0_328[k]
                   - f_7 * smg1_328[k]
                   + f_3 * pc_y[k] * smh_459[k];

        t_613[k] = f_17 * slh_334[k]
                   + f_8 * smg0_329[k]
                   - f_9 * smg1_329[k]
                   + f_3 * pc_y[k] * smh_460[k];

        t_614[k] = f_17 * slh_335[k]
                   + f_3 * pc_y[k] * smh_461[k];

        t_615[k] = f_1 * smg0_329[k]
                   - f_2 * smg1_329[k]
                   + f_3 * pc_z[k] * smh_461[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pb_z, pc_y, pc_z, sli0_420, sli0_423, \
                         slh_315, slh_336, sli1_420, sli1_423, \
                         smh_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pb_z[k] * sli0_420[k]
                   - f_10 * pc_z[k] * sli1_420[k];

        t_617[k] = f_18 * slh_336[k]
                   + f_3 * pc_y[k] * smh_462[k];

        t_618[k] = f_11 * slh_315[k]
                   + f_3 * pc_z[k] * smh_462[k];

        t_619[k] = pb_z[k] * sli0_423[k]
                   - f_10 * pc_z[k] * sli1_423[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pb_z, pc_x, pc_y, pc_z, sli0_426, slh_338, \
                         slh_467, sli1_426, smg0_335, smg1_335, smh_464, \
                         smh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_18 * slh_338[k]
                   + f_3 * pc_y[k] * smh_464[k];

        t_621[k] = f_13 * slh_467[k]
                   + f_4 * smg0_335[k]
                   - f_5 * smg1_335[k]
                   + f_3 * pc_x[k] * smh_467[k];

        t_622[k] = pb_z[k] * sli0_426[k]
                   - f_10 * pc_z[k] * sli1_426[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, slh_318, slh_341, slh_471, \
                         smg0_339, smg1_339, smh_465, smh_467, \
                         smh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_11 * slh_318[k]
                   + f_3 * pc_z[k] * smh_465[k];

        t_624[k] = f_18 * slh_341[k]
                   + f_3 * pc_y[k] * smh_467[k];

        t_625[k] = f_13 * slh_471[k]
                   + f_6 * smg0_339[k]
                   - f_7 * smg1_339[k]
                   + f_3 * pc_x[k] * smh_471[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pb_z, pc_y, pc_z, sli0_430, sli0_432, \
                         slh_321, slh_322, slh_345, sli1_430, sli1_432, smh_468, \
                         smh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pb_z[k] * sli0_430[k]
                   - f_10 * pc_z[k] * sli1_430[k];

        t_627[k] = f_11 * slh_321[k]
                   + f_3 * pc_z[k] * smh_468[k];

        t_628[k] = pb_z[k] * sli0_432[k]
                   + f_12 * slh_322[k]
                   - f_10 * pc_z[k] * sli1_432[k];

        t_629[k] = f_18 * slh_345[k]
                   + f_3 * pc_y[k] * smh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, slh_476, slh_477, slh_478, slh_479, \
                         smg0_344, smg1_344, smh_476, smh_477, smh_478, \
                         smh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_13 * slh_476[k]
                   + f_8 * smg0_344[k]
                   - f_9 * smg1_344[k]
                   + f_3 * pc_x[k] * smh_476[k];

        t_631[k] = f_13 * slh_477[k]
                   + f_3 * pc_x[k] * smh_477[k];

        t_632[k] = f_13 * slh_478[k]
                   + f_3 * pc_x[k] * smh_478[k];

        t_633[k] = f_13 * slh_479[k]
                   + f_3 * pc_x[k] * smh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pb_z, pc_x, pc_z, sli0_441, slh_480, \
                         slh_481, slh_482, sli1_441, smh_480, smh_481, \
                         smh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_13 * slh_480[k]
                   + f_3 * pc_x[k] * smh_480[k];

        t_635[k] = f_13 * slh_481[k]
                   + f_3 * pc_x[k] * smh_481[k];

        t_636[k] = f_13 * slh_482[k]
                   + f_3 * pc_x[k] * smh_482[k];

        t_637[k] = pb_z[k] * sli0_441[k]
                   - f_10 * pc_z[k] * sli1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, slh_330, slh_353, slh_354, smg0_342, \
                         smg0_343, smg1_342, smg1_343, smh_477, smh_479, \
                         smh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * slh_330[k]
                   + f_3 * pc_z[k] * smh_477[k];

        t_639[k] = f_18 * slh_353[k]
                   + f_4 * smg0_342[k]
                   - f_5 * smg1_342[k]
                   + f_3 * pc_y[k] * smh_479[k];

        t_640[k] = f_18 * slh_354[k]
                   + f_6 * smg0_343[k]
                   - f_7 * smg1_343[k]
                   + f_3 * pc_y[k] * smh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, slh_335, slh_355, slh_356, smg0_344, \
                         smg1_344, smh_481, smh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_18 * slh_355[k]
                   + f_8 * smg0_344[k]
                   - f_9 * smg1_344[k]
                   + f_3 * pc_y[k] * smh_481[k];

        t_642[k] = f_18 * slh_356[k]
                   + f_3 * pc_y[k] * smh_482[k];

        t_643[k] = f_11 * slh_335[k]
                   + f_1 * smg0_344[k]
                   - f_2 * smg1_344[k]
                   + f_3 * pc_z[k] * smh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, slh_336, slh_357, slh_483, \
                         smg0_345, smg1_345, smh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_13 * slh_483[k]
                   + f_1 * smg0_345[k]
                   - f_2 * smg1_345[k]
                   + f_3 * pc_x[k] * smh_483[k];

        t_645[k] = f_14 * slh_357[k]
                   + f_3 * pc_y[k] * smh_483[k];

        t_646[k] = f_12 * slh_336[k]
                   + f_3 * pc_z[k] * smh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, slh_359, slh_486, slh_488, smg0_348, \
                         smg0_350, smg1_348, smg1_350, smh_485, smh_486, \
                         smh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_13 * slh_486[k]
                   + f_4 * smg0_348[k]
                   - f_5 * smg1_348[k]
                   + f_3 * pc_x[k] * smh_486[k];

        t_648[k] = f_14 * slh_359[k]
                   + f_3 * pc_y[k] * smh_485[k];

        t_649[k] = f_13 * slh_488[k]
                   + f_4 * smg0_350[k]
                   - f_5 * smg1_350[k]
                   + f_3 * pc_x[k] * smh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, slh_339, slh_362, slh_489, \
                         smg0_351, smg1_351, smh_486, smh_488, \
                         smh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_13 * slh_489[k]
                   + f_6 * smg0_351[k]
                   - f_7 * smg1_351[k]
                   + f_3 * pc_x[k] * smh_489[k];

        t_651[k] = f_12 * slh_339[k]
                   + f_3 * pc_z[k] * smh_486[k];

        t_652[k] = f_14 * slh_362[k]
                   + f_3 * pc_y[k] * smh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, slh_342, slh_492, slh_493, smg0_354, \
                         smg0_355, smg1_354, smg1_355, smh_489, smh_492, \
                         smh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_13 * slh_492[k]
                   + f_6 * smg0_354[k]
                   - f_7 * smg1_354[k]
                   + f_3 * pc_x[k] * smh_492[k];

        t_654[k] = f_13 * slh_493[k]
                   + f_8 * smg0_355[k]
                   - f_9 * smg1_355[k]
                   + f_3 * pc_x[k] * smh_493[k];

        t_655[k] = f_12 * slh_342[k]
                   + f_3 * pc_z[k] * smh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, slh_366, slh_495, slh_497, smg0_357, \
                         smg0_359, smg1_357, smg1_359, smh_492, smh_495, \
                         smh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_13 * slh_495[k]
                   + f_8 * smg0_357[k]
                   - f_9 * smg1_357[k]
                   + f_3 * pc_x[k] * smh_495[k];

        t_657[k] = f_14 * slh_366[k]
                   + f_3 * pc_y[k] * smh_492[k];

        t_658[k] = f_13 * slh_497[k]
                   + f_8 * smg0_359[k]
                   - f_9 * smg1_359[k]
                   + f_3 * pc_x[k] * smh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, slh_498, slh_499, slh_500, \
                         slh_501, slh_502, smh_498, smh_499, smh_500, smh_501, \
                         smh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_13 * slh_498[k]
                   + f_3 * pc_x[k] * smh_498[k];

        t_660[k] = f_13 * slh_499[k]
                   + f_3 * pc_x[k] * smh_499[k];

        t_661[k] = f_13 * slh_500[k]
                   + f_3 * pc_x[k] * smh_500[k];

        t_662[k] = f_13 * slh_501[k]
                   + f_3 * pc_x[k] * smh_501[k];

        t_663[k] = f_13 * slh_502[k]
                   + f_3 * pc_x[k] * smh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, slh_351, slh_372, slh_503, \
                         smg0_355, smg1_355, smh_498, smh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_13 * slh_503[k]
                   + f_3 * pc_x[k] * smh_503[k];

        t_665[k] = f_14 * slh_372[k]
                   + f_1 * smg0_355[k]
                   - f_2 * smg1_355[k]
                   + f_3 * pc_y[k] * smh_498[k];

        t_666[k] = f_12 * slh_351[k]
                   + f_3 * pc_z[k] * smh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, slh_374, slh_375, slh_376, smg0_357, \
                         smg0_358, smg0_359, smg1_357, smg1_358, smg1_359, smh_500, smh_501, \
                         smh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * slh_374[k]
                   + f_4 * smg0_357[k]
                   - f_5 * smg1_357[k]
                   + f_3 * pc_y[k] * smh_500[k];

        t_668[k] = f_14 * slh_375[k]
                   + f_6 * smg0_358[k]
                   - f_7 * smg1_358[k]
                   + f_3 * pc_y[k] * smh_501[k];

        t_669[k] = f_14 * slh_376[k]
                   + f_8 * smg0_359[k]
                   - f_9 * smg1_359[k]
                   + f_3 * pc_y[k] * smh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, slh_356, slh_377, slh_504, \
                         smg0_359, smg0_360, smg1_359, smg1_360, smh_503, \
                         smh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * slh_377[k]
                   + f_3 * pc_y[k] * smh_503[k];

        t_671[k] = f_12 * slh_356[k]
                   + f_1 * smg0_359[k]
                   - f_2 * smg1_359[k]
                   + f_3 * pc_z[k] * smh_503[k];

        t_672[k] = f_13 * slh_504[k]
                   + f_1 * smg0_360[k]
                   - f_2 * smg1_360[k]
                   + f_3 * pc_x[k] * smh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, slh_357, slh_378, \
                         slh_380, slh_507, smg0_363, smg1_363, smh_504, smh_506, \
                         smh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * slh_378[k]
                   + f_3 * pc_y[k] * smh_504[k];

        t_674[k] = f_13 * slh_357[k]
                   + f_3 * pc_z[k] * smh_504[k];

        t_675[k] = f_13 * slh_507[k]
                   + f_4 * smg0_363[k]
                   - f_5 * smg1_363[k]
                   + f_3 * pc_x[k] * smh_507[k];

        t_676[k] = f_13 * slh_380[k]
                   + f_3 * pc_y[k] * smh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, slh_360, slh_509, slh_510, smg0_365, \
                         smg0_366, smg1_365, smg1_366, smh_507, smh_509, \
                         smh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_13 * slh_509[k]
                   + f_4 * smg0_365[k]
                   - f_5 * smg1_365[k]
                   + f_3 * pc_x[k] * smh_509[k];

        t_678[k] = f_13 * slh_510[k]
                   + f_6 * smg0_366[k]
                   - f_7 * smg1_366[k]
                   + f_3 * pc_x[k] * smh_510[k];

        t_679[k] = f_13 * slh_360[k]
                   + f_3 * pc_z[k] * smh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, slh_383, slh_513, slh_514, smg0_369, \
                         smg0_370, smg1_369, smg1_370, smh_509, smh_513, \
                         smh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * slh_383[k]
                   + f_3 * pc_y[k] * smh_509[k];

        t_681[k] = f_13 * slh_513[k]
                   + f_6 * smg0_369[k]
                   - f_7 * smg1_369[k]
                   + f_3 * pc_x[k] * smh_513[k];

        t_682[k] = f_13 * slh_514[k]
                   + f_8 * smg0_370[k]
                   - f_9 * smg1_370[k]
                   + f_3 * pc_x[k] * smh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, slh_363, slh_387, slh_516, \
                         smg0_372, smg1_372, smh_510, smh_513, \
                         smh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * slh_363[k]
                   + f_3 * pc_z[k] * smh_510[k];

        t_684[k] = f_13 * slh_516[k]
                   + f_8 * smg0_372[k]
                   - f_9 * smg1_372[k]
                   + f_3 * pc_x[k] * smh_516[k];

        t_685[k] = f_13 * slh_387[k]
                   + f_3 * pc_y[k] * smh_513[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_560 = buffer.data(sli0 + 560);
    const auto *sli0_563 = buffer.data(sli0 + 563);
    const auto *sli0_565 = buffer.data(sli0 + 565);
    const auto *sli0_566 = buffer.data(sli0 + 566);
    const auto *sli0_569 = buffer.data(sli0 + 569);
    const auto *sli0_570 = buffer.data(sli0 + 570);
    const auto *sli0_572 = buffer.data(sli0 + 572);
    const auto *sli0_574 = buffer.data(sli0 + 574);
    const auto *sli0_587 = buffer.data(sli0 + 587);

    const auto *slh_372 = buffer.data(slh + 372);
    const auto *slh_377 = buffer.data(slh + 377);
    const auto *slh_378 = buffer.data(slh + 378);
    const auto *slh_381 = buffer.data(slh + 381);
    const auto *slh_384 = buffer.data(slh + 384);
    const auto *slh_393 = buffer.data(slh + 393);
    const auto *slh_395 = buffer.data(slh + 395);
    const auto *slh_396 = buffer.data(slh + 396);
    const auto *slh_397 = buffer.data(slh + 397);
    const auto *slh_398 = buffer.data(slh + 398);
    const auto *slh_399 = buffer.data(slh + 399);
    const auto *slh_401 = buffer.data(slh + 401);
    const auto *slh_402 = buffer.data(slh + 402);
    const auto *slh_404 = buffer.data(slh + 404);
    const auto *slh_405 = buffer.data(slh + 405);
    const auto *slh_408 = buffer.data(slh + 408);
    const auto *slh_414 = buffer.data(slh + 414);
    const auto *slh_416 = buffer.data(slh + 416);
    const auto *slh_417 = buffer.data(slh + 417);
    const auto *slh_418 = buffer.data(slh + 418);
    const auto *slh_419 = buffer.data(slh + 419);
    const auto *slh_420 = buffer.data(slh + 420);
    const auto *slh_421 = buffer.data(slh + 421);
    const auto *slh_422 = buffer.data(slh + 422);
    const auto *slh_423 = buffer.data(slh + 423);
    const auto *slh_425 = buffer.data(slh + 425);
    const auto *slh_426 = buffer.data(slh + 426);
    const auto *slh_428 = buffer.data(slh + 428);
    const auto *slh_429 = buffer.data(slh + 429);
    const auto *slh_435 = buffer.data(slh + 435);
    const auto *slh_437 = buffer.data(slh + 437);
    const auto *slh_438 = buffer.data(slh + 438);
    const auto *slh_439 = buffer.data(slh + 439);
    const auto *slh_440 = buffer.data(slh + 440);
    const auto *slh_441 = buffer.data(slh + 441);
    const auto *slh_443 = buffer.data(slh + 443);
    const auto *slh_446 = buffer.data(slh + 446);
    const auto *slh_518 = buffer.data(slh + 518);
    const auto *slh_519 = buffer.data(slh + 519);
    const auto *slh_520 = buffer.data(slh + 520);
    const auto *slh_521 = buffer.data(slh + 521);
    const auto *slh_522 = buffer.data(slh + 522);
    const auto *slh_523 = buffer.data(slh + 523);
    const auto *slh_524 = buffer.data(slh + 524);
    const auto *slh_525 = buffer.data(slh + 525);
    const auto *slh_528 = buffer.data(slh + 528);
    const auto *slh_530 = buffer.data(slh + 530);
    const auto *slh_531 = buffer.data(slh + 531);
    const auto *slh_534 = buffer.data(slh + 534);
    const auto *slh_535 = buffer.data(slh + 535);
    const auto *slh_537 = buffer.data(slh + 537);
    const auto *slh_539 = buffer.data(slh + 539);
    const auto *slh_540 = buffer.data(slh + 540);
    const auto *slh_541 = buffer.data(slh + 541);
    const auto *slh_542 = buffer.data(slh + 542);
    const auto *slh_543 = buffer.data(slh + 543);
    const auto *slh_544 = buffer.data(slh + 544);
    const auto *slh_545 = buffer.data(slh + 545);
    const auto *slh_561 = buffer.data(slh + 561);
    const auto *slh_562 = buffer.data(slh + 562);
    const auto *slh_563 = buffer.data(slh + 563);
    const auto *slh_564 = buffer.data(slh + 564);
    const auto *slh_565 = buffer.data(slh + 565);
    const auto *slh_566 = buffer.data(slh + 566);
    const auto *slh_567 = buffer.data(slh + 567);
    const auto *slh_570 = buffer.data(slh + 570);
    const auto *slh_572 = buffer.data(slh + 572);
    const auto *slh_573 = buffer.data(slh + 573);
    const auto *slh_576 = buffer.data(slh + 576);
    const auto *slh_577 = buffer.data(slh + 577);
    const auto *slh_579 = buffer.data(slh + 579);
    const auto *slh_581 = buffer.data(slh + 581);
    const auto *slh_582 = buffer.data(slh + 582);
    const auto *slh_583 = buffer.data(slh + 583);
    const auto *slh_584 = buffer.data(slh + 584);
    const auto *slh_585 = buffer.data(slh + 585);
    const auto *slh_586 = buffer.data(slh + 586);
    const auto *slh_587 = buffer.data(slh + 587);
    const auto *slh_588 = buffer.data(slh + 588);
    const auto *slh_591 = buffer.data(slh + 591);
    const auto *slh_593 = buffer.data(slh + 593);
    const auto *slh_594 = buffer.data(slh + 594);
    const auto *slh_597 = buffer.data(slh + 597);
    const auto *slh_598 = buffer.data(slh + 598);
    const auto *slh_600 = buffer.data(slh + 600);

    const auto *sli1_560 = buffer.data(sli1 + 560);
    const auto *sli1_563 = buffer.data(sli1 + 563);
    const auto *sli1_565 = buffer.data(sli1 + 565);
    const auto *sli1_566 = buffer.data(sli1 + 566);
    const auto *sli1_569 = buffer.data(sli1 + 569);
    const auto *sli1_570 = buffer.data(sli1 + 570);
    const auto *sli1_572 = buffer.data(sli1 + 572);
    const auto *sli1_574 = buffer.data(sli1 + 574);
    const auto *sli1_587 = buffer.data(sli1 + 587);

    const auto *smg0_370 = buffer.data(smg0 + 370);
    const auto *smg0_372 = buffer.data(smg0 + 372);
    const auto *smg0_373 = buffer.data(smg0 + 373);
    const auto *smg0_374 = buffer.data(smg0 + 374);
    const auto *smg0_375 = buffer.data(smg0 + 375);
    const auto *smg0_378 = buffer.data(smg0 + 378);
    const auto *smg0_380 = buffer.data(smg0 + 380);
    const auto *smg0_381 = buffer.data(smg0 + 381);
    const auto *smg0_384 = buffer.data(smg0 + 384);
    const auto *smg0_385 = buffer.data(smg0 + 385);
    const auto *smg0_387 = buffer.data(smg0 + 387);
    const auto *smg0_388 = buffer.data(smg0 + 388);
    const auto *smg0_389 = buffer.data(smg0 + 389);
    const auto *smg0_400 = buffer.data(smg0 + 400);
    const auto *smg0_402 = buffer.data(smg0 + 402);
    const auto *smg0_403 = buffer.data(smg0 + 403);
    const auto *smg0_404 = buffer.data(smg0 + 404);
    const auto *smg0_405 = buffer.data(smg0 + 405);
    const auto *smg0_408 = buffer.data(smg0 + 408);
    const auto *smg0_410 = buffer.data(smg0 + 410);
    const auto *smg0_411 = buffer.data(smg0 + 411);
    const auto *smg0_414 = buffer.data(smg0 + 414);
    const auto *smg0_415 = buffer.data(smg0 + 415);
    const auto *smg0_417 = buffer.data(smg0 + 417);
    const auto *smg0_418 = buffer.data(smg0 + 418);
    const auto *smg0_419 = buffer.data(smg0 + 419);
    const auto *smg0_420 = buffer.data(smg0 + 420);
    const auto *smg0_423 = buffer.data(smg0 + 423);
    const auto *smg0_425 = buffer.data(smg0 + 425);
    const auto *smg0_426 = buffer.data(smg0 + 426);
    const auto *smg0_429 = buffer.data(smg0 + 429);
    const auto *smg0_430 = buffer.data(smg0 + 430);
    const auto *smg0_432 = buffer.data(smg0 + 432);

    const auto *smg1_370 = buffer.data(smg1 + 370);
    const auto *smg1_372 = buffer.data(smg1 + 372);
    const auto *smg1_373 = buffer.data(smg1 + 373);
    const auto *smg1_374 = buffer.data(smg1 + 374);
    const auto *smg1_375 = buffer.data(smg1 + 375);
    const auto *smg1_378 = buffer.data(smg1 + 378);
    const auto *smg1_380 = buffer.data(smg1 + 380);
    const auto *smg1_381 = buffer.data(smg1 + 381);
    const auto *smg1_384 = buffer.data(smg1 + 384);
    const auto *smg1_385 = buffer.data(smg1 + 385);
    const auto *smg1_387 = buffer.data(smg1 + 387);
    const auto *smg1_388 = buffer.data(smg1 + 388);
    const auto *smg1_389 = buffer.data(smg1 + 389);
    const auto *smg1_400 = buffer.data(smg1 + 400);
    const auto *smg1_402 = buffer.data(smg1 + 402);
    const auto *smg1_403 = buffer.data(smg1 + 403);
    const auto *smg1_404 = buffer.data(smg1 + 404);
    const auto *smg1_405 = buffer.data(smg1 + 405);
    const auto *smg1_408 = buffer.data(smg1 + 408);
    const auto *smg1_410 = buffer.data(smg1 + 410);
    const auto *smg1_411 = buffer.data(smg1 + 411);
    const auto *smg1_414 = buffer.data(smg1 + 414);
    const auto *smg1_415 = buffer.data(smg1 + 415);
    const auto *smg1_417 = buffer.data(smg1 + 417);
    const auto *smg1_418 = buffer.data(smg1 + 418);
    const auto *smg1_419 = buffer.data(smg1 + 419);
    const auto *smg1_420 = buffer.data(smg1 + 420);
    const auto *smg1_423 = buffer.data(smg1 + 423);
    const auto *smg1_425 = buffer.data(smg1 + 425);
    const auto *smg1_426 = buffer.data(smg1 + 426);
    const auto *smg1_429 = buffer.data(smg1 + 429);
    const auto *smg1_430 = buffer.data(smg1 + 430);
    const auto *smg1_432 = buffer.data(smg1 + 432);

    const auto *smh_518 = buffer.data(smh + 518);
    const auto *smh_519 = buffer.data(smh + 519);
    const auto *smh_520 = buffer.data(smh + 520);
    const auto *smh_521 = buffer.data(smh + 521);
    const auto *smh_522 = buffer.data(smh + 522);
    const auto *smh_523 = buffer.data(smh + 523);
    const auto *smh_524 = buffer.data(smh + 524);
    const auto *smh_525 = buffer.data(smh + 525);
    const auto *smh_527 = buffer.data(smh + 527);
    const auto *smh_528 = buffer.data(smh + 528);
    const auto *smh_530 = buffer.data(smh + 530);
    const auto *smh_531 = buffer.data(smh + 531);
    const auto *smh_534 = buffer.data(smh + 534);
    const auto *smh_535 = buffer.data(smh + 535);
    const auto *smh_537 = buffer.data(smh + 537);
    const auto *smh_539 = buffer.data(smh + 539);
    const auto *smh_540 = buffer.data(smh + 540);
    const auto *smh_541 = buffer.data(smh + 541);
    const auto *smh_542 = buffer.data(smh + 542);
    const auto *smh_543 = buffer.data(smh + 543);
    const auto *smh_544 = buffer.data(smh + 544);
    const auto *smh_545 = buffer.data(smh + 545);
    const auto *smh_546 = buffer.data(smh + 546);
    const auto *smh_548 = buffer.data(smh + 548);
    const auto *smh_549 = buffer.data(smh + 549);
    const auto *smh_551 = buffer.data(smh + 551);
    const auto *smh_552 = buffer.data(smh + 552);
    const auto *smh_555 = buffer.data(smh + 555);
    const auto *smh_561 = buffer.data(smh + 561);
    const auto *smh_562 = buffer.data(smh + 562);
    const auto *smh_563 = buffer.data(smh + 563);
    const auto *smh_564 = buffer.data(smh + 564);
    const auto *smh_565 = buffer.data(smh + 565);
    const auto *smh_566 = buffer.data(smh + 566);
    const auto *smh_567 = buffer.data(smh + 567);
    const auto *smh_569 = buffer.data(smh + 569);
    const auto *smh_570 = buffer.data(smh + 570);
    const auto *smh_572 = buffer.data(smh + 572);
    const auto *smh_573 = buffer.data(smh + 573);
    const auto *smh_576 = buffer.data(smh + 576);
    const auto *smh_577 = buffer.data(smh + 577);
    const auto *smh_579 = buffer.data(smh + 579);
    const auto *smh_581 = buffer.data(smh + 581);
    const auto *smh_582 = buffer.data(smh + 582);
    const auto *smh_583 = buffer.data(smh + 583);
    const auto *smh_584 = buffer.data(smh + 584);
    const auto *smh_585 = buffer.data(smh + 585);
    const auto *smh_586 = buffer.data(smh + 586);
    const auto *smh_587 = buffer.data(smh + 587);
    const auto *smh_588 = buffer.data(smh + 588);
    const auto *smh_590 = buffer.data(smh + 590);
    const auto *smh_591 = buffer.data(smh + 591);
    const auto *smh_593 = buffer.data(smh + 593);
    const auto *smh_594 = buffer.data(smh + 594);
    const auto *smh_597 = buffer.data(smh + 597);
    const auto *smh_598 = buffer.data(smh + 598);
    const auto *smh_600 = buffer.data(smh + 600);

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, slh_518, slh_519, slh_520, slh_521, \
                         smg0_374, smg1_374, smh_518, smh_519, smh_520, \
                         smh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_13 * slh_518[k]
                   + f_8 * smg0_374[k]
                   - f_9 * smg1_374[k]
                   + f_3 * pc_x[k] * smh_518[k];

        t_687[k] = f_13 * slh_519[k]
                   + f_3 * pc_x[k] * smh_519[k];

        t_688[k] = f_13 * slh_520[k]
                   + f_3 * pc_x[k] * smh_520[k];

        t_689[k] = f_13 * slh_521[k]
                   + f_3 * pc_x[k] * smh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, slh_393, slh_522, slh_523, \
                         slh_524, smg0_370, smg1_370, smh_519, smh_522, smh_523, \
                         smh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_13 * slh_522[k]
                   + f_3 * pc_x[k] * smh_522[k];

        t_691[k] = f_13 * slh_523[k]
                   + f_3 * pc_x[k] * smh_523[k];

        t_692[k] = f_13 * slh_524[k]
                   + f_3 * pc_x[k] * smh_524[k];

        t_693[k] = f_13 * slh_393[k]
                   + f_1 * smg0_370[k]
                   - f_2 * smg1_370[k]
                   + f_3 * pc_y[k] * smh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, slh_372, slh_395, slh_396, smg0_372, \
                         smg0_373, smg1_372, smg1_373, smh_519, smh_521, \
                         smh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * slh_372[k]
                   + f_3 * pc_z[k] * smh_519[k];

        t_695[k] = f_13 * slh_395[k]
                   + f_4 * smg0_372[k]
                   - f_5 * smg1_372[k]
                   + f_3 * pc_y[k] * smh_521[k];

        t_696[k] = f_13 * slh_396[k]
                   + f_6 * smg0_373[k]
                   - f_7 * smg1_373[k]
                   + f_3 * pc_y[k] * smh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, slh_377, slh_397, slh_398, smg0_374, \
                         smg1_374, smh_523, smh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * slh_397[k]
                   + f_8 * smg0_374[k]
                   - f_9 * smg1_374[k]
                   + f_3 * pc_y[k] * smh_523[k];

        t_698[k] = f_13 * slh_398[k]
                   + f_3 * pc_y[k] * smh_524[k];

        t_699[k] = f_13 * slh_377[k]
                   + f_1 * smg0_374[k]
                   - f_2 * smg1_374[k]
                   + f_3 * pc_z[k] * smh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, slh_378, slh_399, slh_525, \
                         smg0_375, smg1_375, smh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_13 * slh_525[k]
                   + f_1 * smg0_375[k]
                   - f_2 * smg1_375[k]
                   + f_3 * pc_x[k] * smh_525[k];

        t_701[k] = f_12 * slh_399[k]
                   + f_3 * pc_y[k] * smh_525[k];

        t_702[k] = f_14 * slh_378[k]
                   + f_3 * pc_z[k] * smh_525[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, slh_401, slh_528, slh_530, smg0_378, \
                         smg0_380, smg1_378, smg1_380, smh_527, smh_528, \
                         smh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_13 * slh_528[k]
                   + f_4 * smg0_378[k]
                   - f_5 * smg1_378[k]
                   + f_3 * pc_x[k] * smh_528[k];

        t_704[k] = f_12 * slh_401[k]
                   + f_3 * pc_y[k] * smh_527[k];

        t_705[k] = f_13 * slh_530[k]
                   + f_4 * smg0_380[k]
                   - f_5 * smg1_380[k]
                   + f_3 * pc_x[k] * smh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, slh_381, slh_404, slh_531, \
                         smg0_381, smg1_381, smh_528, smh_530, \
                         smh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_13 * slh_531[k]
                   + f_6 * smg0_381[k]
                   - f_7 * smg1_381[k]
                   + f_3 * pc_x[k] * smh_531[k];

        t_707[k] = f_14 * slh_381[k]
                   + f_3 * pc_z[k] * smh_528[k];

        t_708[k] = f_12 * slh_404[k]
                   + f_3 * pc_y[k] * smh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, slh_384, slh_534, slh_535, smg0_384, \
                         smg0_385, smg1_384, smg1_385, smh_531, smh_534, \
                         smh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_13 * slh_534[k]
                   + f_6 * smg0_384[k]
                   - f_7 * smg1_384[k]
                   + f_3 * pc_x[k] * smh_534[k];

        t_710[k] = f_13 * slh_535[k]
                   + f_8 * smg0_385[k]
                   - f_9 * smg1_385[k]
                   + f_3 * pc_x[k] * smh_535[k];

        t_711[k] = f_14 * slh_384[k]
                   + f_3 * pc_z[k] * smh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, slh_408, slh_537, slh_539, smg0_387, \
                         smg0_389, smg1_387, smg1_389, smh_534, smh_537, \
                         smh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_13 * slh_537[k]
                   + f_8 * smg0_387[k]
                   - f_9 * smg1_387[k]
                   + f_3 * pc_x[k] * smh_537[k];

        t_713[k] = f_12 * slh_408[k]
                   + f_3 * pc_y[k] * smh_534[k];

        t_714[k] = f_13 * slh_539[k]
                   + f_8 * smg0_389[k]
                   - f_9 * smg1_389[k]
                   + f_3 * pc_x[k] * smh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, slh_540, slh_541, slh_542, \
                         slh_543, slh_544, smh_540, smh_541, smh_542, smh_543, \
                         smh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_13 * slh_540[k]
                   + f_3 * pc_x[k] * smh_540[k];

        t_716[k] = f_13 * slh_541[k]
                   + f_3 * pc_x[k] * smh_541[k];

        t_717[k] = f_13 * slh_542[k]
                   + f_3 * pc_x[k] * smh_542[k];

        t_718[k] = f_13 * slh_543[k]
                   + f_3 * pc_x[k] * smh_543[k];

        t_719[k] = f_13 * slh_544[k]
                   + f_3 * pc_x[k] * smh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, slh_393, slh_414, slh_545, \
                         smg0_385, smg1_385, smh_540, smh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_13 * slh_545[k]
                   + f_3 * pc_x[k] * smh_545[k];

        t_721[k] = f_12 * slh_414[k]
                   + f_1 * smg0_385[k]
                   - f_2 * smg1_385[k]
                   + f_3 * pc_y[k] * smh_540[k];

        t_722[k] = f_14 * slh_393[k]
                   + f_3 * pc_z[k] * smh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, slh_416, slh_417, slh_418, smg0_387, \
                         smg0_388, smg0_389, smg1_387, smg1_388, smg1_389, smh_542, smh_543, \
                         smh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * slh_416[k]
                   + f_4 * smg0_387[k]
                   - f_5 * smg1_387[k]
                   + f_3 * pc_y[k] * smh_542[k];

        t_724[k] = f_12 * slh_417[k]
                   + f_6 * smg0_388[k]
                   - f_7 * smg1_388[k]
                   + f_3 * pc_y[k] * smh_543[k];

        t_725[k] = f_12 * slh_418[k]
                   + f_8 * smg0_389[k]
                   - f_9 * smg1_389[k]
                   + f_3 * pc_y[k] * smh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pb_y, pc_y, pc_z, sli0_560, slh_398, \
                         slh_419, slh_420, sli1_560, smg0_389, smg1_389, smh_545, \
                         smh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * slh_419[k]
                   + f_3 * pc_y[k] * smh_545[k];

        t_727[k] = f_14 * slh_398[k]
                   + f_1 * smg0_389[k]
                   - f_2 * smg1_389[k]
                   + f_3 * pc_z[k] * smh_545[k];

        t_728[k] = pb_y[k] * sli0_560[k]
                   - f_10 * pc_y[k] * sli1_560[k];

        t_729[k] = f_11 * slh_420[k]
                   + f_3 * pc_y[k] * smh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pb_y, pc_y, pc_z, sli0_563, sli0_565, \
                         slh_399, slh_421, slh_422, sli1_563, sli1_565, smh_546, \
                         smh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_18 * slh_399[k]
                   + f_3 * pc_z[k] * smh_546[k];

        t_731[k] = pb_y[k] * sli0_563[k]
                   + f_12 * slh_421[k]
                   - f_10 * pc_y[k] * sli1_563[k];

        t_732[k] = f_11 * slh_422[k]
                   + f_3 * pc_y[k] * smh_548[k];

        t_733[k] = pb_y[k] * sli0_565[k]
                   - f_10 * pc_y[k] * sli1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pb_y, pc_y, pc_z, sli0_566, sli0_569, \
                         slh_402, slh_423, slh_425, sli1_566, sli1_569, smh_549, \
                         smh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pb_y[k] * sli0_566[k]
                   + f_13 * slh_423[k]
                   - f_10 * pc_y[k] * sli1_566[k];

        t_735[k] = f_18 * slh_402[k]
                   + f_3 * pc_z[k] * smh_549[k];

        t_736[k] = f_11 * slh_425[k]
                   + f_3 * pc_y[k] * smh_551[k];

        t_737[k] = pb_y[k] * sli0_569[k]
                   - f_10 * pc_y[k] * sli1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pb_y, pc_y, pc_z, sli0_570, sli0_572, slh_405, \
                         slh_426, slh_428, sli1_570, sli1_572, \
                         smh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pb_y[k] * sli0_570[k]
                   + f_14 * slh_426[k]
                   - f_10 * pc_y[k] * sli1_570[k];

        t_739[k] = f_18 * slh_405[k]
                   + f_3 * pc_z[k] * smh_552[k];

        t_740[k] = pb_y[k] * sli0_572[k]
                   + f_12 * slh_428[k]
                   - f_10 * pc_y[k] * sli1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pb_y, pc_x, pc_y, sli0_574, slh_429, \
                         slh_561, slh_562, sli1_574, smh_555, smh_561, \
                         smh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * slh_429[k]
                   + f_3 * pc_y[k] * smh_555[k];

        t_742[k] = pb_y[k] * sli0_574[k]
                   - f_10 * pc_y[k] * sli1_574[k];

        t_743[k] = f_13 * slh_561[k]
                   + f_3 * pc_x[k] * smh_561[k];

        t_744[k] = f_13 * slh_562[k]
                   + f_3 * pc_x[k] * smh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, slh_563, slh_564, slh_565, slh_566, \
                         smh_563, smh_564, smh_565, smh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_13 * slh_563[k]
                   + f_3 * pc_x[k] * smh_563[k];

        t_746[k] = f_13 * slh_564[k]
                   + f_3 * pc_x[k] * smh_564[k];

        t_747[k] = f_13 * slh_565[k]
                   + f_3 * pc_x[k] * smh_565[k];

        t_748[k] = f_13 * slh_566[k]
                   + f_3 * pc_x[k] * smh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, slh_414, slh_435, slh_437, smg0_400, \
                         smg0_402, smg1_400, smg1_402, smh_561, \
                         smh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * slh_435[k]
                   + f_1 * smg0_400[k]
                   - f_2 * smg1_400[k]
                   + f_3 * pc_y[k] * smh_561[k];

        t_750[k] = f_18 * slh_414[k]
                   + f_3 * pc_z[k] * smh_561[k];

        t_751[k] = f_11 * slh_437[k]
                   + f_4 * smg0_402[k]
                   - f_5 * smg1_402[k]
                   + f_3 * pc_y[k] * smh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, slh_438, slh_439, slh_440, smg0_403, \
                         smg0_404, smg1_403, smg1_404, smh_564, smh_565, \
                         smh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * slh_438[k]
                   + f_6 * smg0_403[k]
                   - f_7 * smg1_403[k]
                   + f_3 * pc_y[k] * smh_564[k];

        t_753[k] = f_11 * slh_439[k]
                   + f_8 * smg0_404[k]
                   - f_9 * smg1_404[k]
                   + f_3 * pc_y[k] * smh_565[k];

        t_754[k] = f_11 * slh_440[k]
                   + f_3 * pc_y[k] * smh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pb_y, pc_x, pc_y, pc_z, sli0_587, \
                         slh_420, slh_567, sli1_587, smg0_405, smg1_405, \
                         smh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pb_y[k] * sli0_587[k]
                   - f_10 * pc_y[k] * sli1_587[k];

        t_756[k] = f_13 * slh_567[k]
                   + f_1 * smg0_405[k]
                   - f_2 * smg1_405[k]
                   + f_3 * pc_x[k] * smh_567[k];

        t_757[k] = f_3 * pc_y[k] * smh_567[k];

        t_758[k] = f_17 * slh_420[k]
                   + f_3 * pc_z[k] * smh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, slh_570, slh_572, smg0_408, \
                         smg0_410, smg1_408, smg1_410, smh_569, smh_570, \
                         smh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_13 * slh_570[k]
                   + f_4 * smg0_408[k]
                   - f_5 * smg1_408[k]
                   + f_3 * pc_x[k] * smh_570[k];

        t_760[k] = f_3 * pc_y[k] * smh_569[k];

        t_761[k] = f_13 * slh_572[k]
                   + f_4 * smg0_410[k]
                   - f_5 * smg1_410[k]
                   + f_3 * pc_x[k] * smh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_x, pc_y, pc_z, slh_423, slh_573, smg0_411, \
                         smg1_411, smh_570, smh_572, smh_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_13 * slh_573[k]
                   + f_6 * smg0_411[k]
                   - f_7 * smg1_411[k]
                   + f_3 * pc_x[k] * smh_573[k];

        t_763[k] = f_17 * slh_423[k]
                   + f_3 * pc_z[k] * smh_570[k];

        t_764[k] = f_3 * pc_y[k] * smh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_z, slh_426, slh_576, slh_577, smg0_414, \
                         smg0_415, smg1_414, smg1_415, smh_573, smh_576, \
                         smh_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_13 * slh_576[k]
                   + f_6 * smg0_414[k]
                   - f_7 * smg1_414[k]
                   + f_3 * pc_x[k] * smh_576[k];

        t_766[k] = f_13 * slh_577[k]
                   + f_8 * smg0_415[k]
                   - f_9 * smg1_415[k]
                   + f_3 * pc_x[k] * smh_577[k];

        t_767[k] = f_17 * slh_426[k]
                   + f_3 * pc_z[k] * smh_573[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_x, pc_y, slh_579, slh_581, smg0_417, \
                         smg0_419, smg1_417, smg1_419, smh_576, smh_579, \
                         smh_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_13 * slh_579[k]
                   + f_8 * smg0_417[k]
                   - f_9 * smg1_417[k]
                   + f_3 * pc_x[k] * smh_579[k];

        t_769[k] = f_3 * pc_y[k] * smh_576[k];

        t_770[k] = f_13 * slh_581[k]
                   + f_8 * smg0_419[k]
                   - f_9 * smg1_419[k]
                   + f_3 * pc_x[k] * smh_581[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, t_775, pc_x, slh_582, slh_583, slh_584, \
                         slh_585, slh_586, smh_582, smh_583, smh_584, smh_585, \
                         smh_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_13 * slh_582[k]
                   + f_3 * pc_x[k] * smh_582[k];

        t_772[k] = f_13 * slh_583[k]
                   + f_3 * pc_x[k] * smh_583[k];

        t_773[k] = f_13 * slh_584[k]
                   + f_3 * pc_x[k] * smh_584[k];

        t_774[k] = f_13 * slh_585[k]
                   + f_3 * pc_x[k] * smh_585[k];

        t_775[k] = f_13 * slh_586[k]
                   + f_3 * pc_x[k] * smh_586[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pc_x, pc_y, pc_z, slh_435, slh_587, \
                         smg0_415, smg0_417, smg1_415, smg1_417, smh_582, smh_584, \
                         smh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_13 * slh_587[k]
                   + f_3 * pc_x[k] * smh_587[k];

        t_777[k] = f_1 * smg0_415[k]
                   - f_2 * smg1_415[k]
                   + f_3 * pc_y[k] * smh_582[k];

        t_778[k] = f_17 * slh_435[k]
                   + f_3 * pc_z[k] * smh_582[k];

        t_779[k] = f_4 * smg0_417[k]
                   - f_5 * smg1_417[k]
                   + f_3 * pc_y[k] * smh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, slh_440, smg0_418, smg0_419, \
                         smg1_418, smg1_419, smh_585, smh_586, \
                         smh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * smg0_418[k]
                   - f_7 * smg1_418[k]
                   + f_3 * pc_y[k] * smh_585[k];

        t_781[k] = f_8 * smg0_419[k]
                   - f_9 * smg1_419[k]
                   + f_3 * pc_y[k] * smh_586[k];

        t_782[k] = f_3 * pc_y[k] * smh_587[k];

        t_783[k] = f_17 * slh_440[k]
                   + f_1 * smg0_419[k]
                   - f_2 * smg1_419[k]
                   + f_3 * pc_z[k] * smh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pc_x, pc_y, pc_z, slh_441, slh_588, \
                         slh_591, smg0_420, smg0_423, smg1_420, smg1_423, smh_588, \
                         smh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_12 * slh_588[k]
                   + f_1 * smg0_420[k]
                   - f_2 * smg1_420[k]
                   + f_3 * pc_x[k] * smh_588[k];

        t_785[k] = f_16 * slh_441[k]
                   + f_3 * pc_y[k] * smh_588[k];

        t_786[k] = f_3 * pc_z[k] * smh_588[k];

        t_787[k] = f_12 * slh_591[k]
                   + f_4 * smg0_423[k]
                   - f_5 * smg1_423[k]
                   + f_3 * pc_x[k] * smh_591[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, pc_y, slh_443, slh_593, slh_594, smg0_425, \
                         smg0_426, smg1_425, smg1_426, smh_590, smh_593, \
                         smh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_16 * slh_443[k]
                   + f_3 * pc_y[k] * smh_590[k];

        t_789[k] = f_12 * slh_593[k]
                   + f_4 * smg0_425[k]
                   - f_5 * smg1_425[k]
                   + f_3 * pc_x[k] * smh_593[k];

        t_790[k] = f_12 * slh_594[k]
                   + f_6 * smg0_426[k]
                   - f_7 * smg1_426[k]
                   + f_3 * pc_x[k] * smh_594[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, pc_x, pc_y, pc_z, slh_446, slh_597, smg0_429, \
                         smg1_429, smh_591, smh_593, smh_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_3 * pc_z[k] * smh_591[k];

        t_792[k] = f_16 * slh_446[k]
                   + f_3 * pc_y[k] * smh_593[k];

        t_793[k] = f_12 * slh_597[k]
                   + f_6 * smg0_429[k]
                   - f_7 * smg1_429[k]
                   + f_3 * pc_x[k] * smh_597[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pc_x, pc_z, slh_598, slh_600, smg0_430, \
                         smg0_432, smg1_430, smg1_432, smh_594, smh_598, \
                         smh_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_12 * slh_598[k]
                   + f_8 * smg0_430[k]
                   - f_9 * smg1_430[k]
                   + f_3 * pc_x[k] * smh_598[k];

        t_795[k] = f_3 * pc_z[k] * smh_594[k];

        t_796[k] = f_12 * slh_600[k]
                   + f_8 * smg0_432[k]
                   - f_9 * smg1_432[k]
                   + f_3 * pc_x[k] * smh_600[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_588 = buffer.data(sli0 + 588);
    const auto *sli0_591 = buffer.data(sli0 + 591);
    const auto *sli0_594 = buffer.data(sli0 + 594);
    const auto *sli0_598 = buffer.data(sli0 + 598);
    const auto *sli0_600 = buffer.data(sli0 + 600);
    const auto *sli0_609 = buffer.data(sli0 + 609);

    const auto *slh_441 = buffer.data(slh + 441);
    const auto *slh_444 = buffer.data(slh + 444);
    const auto *slh_447 = buffer.data(slh + 447);
    const auto *slh_448 = buffer.data(slh + 448);
    const auto *slh_450 = buffer.data(slh + 450);
    const auto *slh_456 = buffer.data(slh + 456);
    const auto *slh_458 = buffer.data(slh + 458);
    const auto *slh_459 = buffer.data(slh + 459);
    const auto *slh_460 = buffer.data(slh + 460);
    const auto *slh_461 = buffer.data(slh + 461);
    const auto *slh_462 = buffer.data(slh + 462);
    const auto *slh_464 = buffer.data(slh + 464);
    const auto *slh_465 = buffer.data(slh + 465);
    const auto *slh_467 = buffer.data(slh + 467);
    const auto *slh_468 = buffer.data(slh + 468);
    const auto *slh_471 = buffer.data(slh + 471);
    const auto *slh_477 = buffer.data(slh + 477);
    const auto *slh_479 = buffer.data(slh + 479);
    const auto *slh_480 = buffer.data(slh + 480);
    const auto *slh_481 = buffer.data(slh + 481);
    const auto *slh_482 = buffer.data(slh + 482);
    const auto *slh_483 = buffer.data(slh + 483);
    const auto *slh_485 = buffer.data(slh + 485);
    const auto *slh_486 = buffer.data(slh + 486);
    const auto *slh_488 = buffer.data(slh + 488);
    const auto *slh_489 = buffer.data(slh + 489);
    const auto *slh_492 = buffer.data(slh + 492);
    const auto *slh_498 = buffer.data(slh + 498);
    const auto *slh_500 = buffer.data(slh + 500);
    const auto *slh_501 = buffer.data(slh + 501);
    const auto *slh_502 = buffer.data(slh + 502);
    const auto *slh_503 = buffer.data(slh + 503);
    const auto *slh_504 = buffer.data(slh + 504);
    const auto *slh_506 = buffer.data(slh + 506);
    const auto *slh_507 = buffer.data(slh + 507);
    const auto *slh_509 = buffer.data(slh + 509);
    const auto *slh_510 = buffer.data(slh + 510);
    const auto *slh_513 = buffer.data(slh + 513);
    const auto *slh_519 = buffer.data(slh + 519);
    const auto *slh_521 = buffer.data(slh + 521);
    const auto *slh_522 = buffer.data(slh + 522);
    const auto *slh_523 = buffer.data(slh + 523);
    const auto *slh_524 = buffer.data(slh + 524);
    const auto *slh_525 = buffer.data(slh + 525);
    const auto *slh_527 = buffer.data(slh + 527);
    const auto *slh_530 = buffer.data(slh + 530);
    const auto *slh_602 = buffer.data(slh + 602);
    const auto *slh_603 = buffer.data(slh + 603);
    const auto *slh_604 = buffer.data(slh + 604);
    const auto *slh_605 = buffer.data(slh + 605);
    const auto *slh_606 = buffer.data(slh + 606);
    const auto *slh_607 = buffer.data(slh + 607);
    const auto *slh_608 = buffer.data(slh + 608);
    const auto *slh_614 = buffer.data(slh + 614);
    const auto *slh_618 = buffer.data(slh + 618);
    const auto *slh_623 = buffer.data(slh + 623);
    const auto *slh_624 = buffer.data(slh + 624);
    const auto *slh_625 = buffer.data(slh + 625);
    const auto *slh_626 = buffer.data(slh + 626);
    const auto *slh_627 = buffer.data(slh + 627);
    const auto *slh_628 = buffer.data(slh + 628);
    const auto *slh_629 = buffer.data(slh + 629);
    const auto *slh_630 = buffer.data(slh + 630);
    const auto *slh_633 = buffer.data(slh + 633);
    const auto *slh_635 = buffer.data(slh + 635);
    const auto *slh_636 = buffer.data(slh + 636);
    const auto *slh_639 = buffer.data(slh + 639);
    const auto *slh_640 = buffer.data(slh + 640);
    const auto *slh_642 = buffer.data(slh + 642);
    const auto *slh_644 = buffer.data(slh + 644);
    const auto *slh_645 = buffer.data(slh + 645);
    const auto *slh_646 = buffer.data(slh + 646);
    const auto *slh_647 = buffer.data(slh + 647);
    const auto *slh_648 = buffer.data(slh + 648);
    const auto *slh_649 = buffer.data(slh + 649);
    const auto *slh_650 = buffer.data(slh + 650);
    const auto *slh_651 = buffer.data(slh + 651);
    const auto *slh_654 = buffer.data(slh + 654);
    const auto *slh_656 = buffer.data(slh + 656);
    const auto *slh_657 = buffer.data(slh + 657);
    const auto *slh_660 = buffer.data(slh + 660);
    const auto *slh_661 = buffer.data(slh + 661);
    const auto *slh_663 = buffer.data(slh + 663);
    const auto *slh_665 = buffer.data(slh + 665);
    const auto *slh_666 = buffer.data(slh + 666);
    const auto *slh_667 = buffer.data(slh + 667);
    const auto *slh_668 = buffer.data(slh + 668);
    const auto *slh_669 = buffer.data(slh + 669);
    const auto *slh_670 = buffer.data(slh + 670);
    const auto *slh_671 = buffer.data(slh + 671);
    const auto *slh_672 = buffer.data(slh + 672);
    const auto *slh_675 = buffer.data(slh + 675);
    const auto *slh_677 = buffer.data(slh + 677);
    const auto *slh_678 = buffer.data(slh + 678);
    const auto *slh_681 = buffer.data(slh + 681);
    const auto *slh_682 = buffer.data(slh + 682);

    const auto *sli1_588 = buffer.data(sli1 + 588);
    const auto *sli1_591 = buffer.data(sli1 + 591);
    const auto *sli1_594 = buffer.data(sli1 + 594);
    const auto *sli1_598 = buffer.data(sli1 + 598);
    const auto *sli1_600 = buffer.data(sli1 + 600);
    const auto *sli1_609 = buffer.data(sli1 + 609);

    const auto *smg0_430 = buffer.data(smg0 + 430);
    const auto *smg0_432 = buffer.data(smg0 + 432);
    const auto *smg0_433 = buffer.data(smg0 + 433);
    const auto *smg0_434 = buffer.data(smg0 + 434);
    const auto *smg0_440 = buffer.data(smg0 + 440);
    const auto *smg0_444 = buffer.data(smg0 + 444);
    const auto *smg0_447 = buffer.data(smg0 + 447);
    const auto *smg0_448 = buffer.data(smg0 + 448);
    const auto *smg0_449 = buffer.data(smg0 + 449);
    const auto *smg0_450 = buffer.data(smg0 + 450);
    const auto *smg0_453 = buffer.data(smg0 + 453);
    const auto *smg0_455 = buffer.data(smg0 + 455);
    const auto *smg0_456 = buffer.data(smg0 + 456);
    const auto *smg0_459 = buffer.data(smg0 + 459);
    const auto *smg0_460 = buffer.data(smg0 + 460);
    const auto *smg0_462 = buffer.data(smg0 + 462);
    const auto *smg0_463 = buffer.data(smg0 + 463);
    const auto *smg0_464 = buffer.data(smg0 + 464);
    const auto *smg0_465 = buffer.data(smg0 + 465);
    const auto *smg0_468 = buffer.data(smg0 + 468);
    const auto *smg0_470 = buffer.data(smg0 + 470);
    const auto *smg0_471 = buffer.data(smg0 + 471);
    const auto *smg0_474 = buffer.data(smg0 + 474);
    const auto *smg0_475 = buffer.data(smg0 + 475);
    const auto *smg0_477 = buffer.data(smg0 + 477);
    const auto *smg0_478 = buffer.data(smg0 + 478);
    const auto *smg0_479 = buffer.data(smg0 + 479);
    const auto *smg0_480 = buffer.data(smg0 + 480);
    const auto *smg0_483 = buffer.data(smg0 + 483);
    const auto *smg0_485 = buffer.data(smg0 + 485);
    const auto *smg0_486 = buffer.data(smg0 + 486);
    const auto *smg0_489 = buffer.data(smg0 + 489);
    const auto *smg0_490 = buffer.data(smg0 + 490);

    const auto *smg1_430 = buffer.data(smg1 + 430);
    const auto *smg1_432 = buffer.data(smg1 + 432);
    const auto *smg1_433 = buffer.data(smg1 + 433);
    const auto *smg1_434 = buffer.data(smg1 + 434);
    const auto *smg1_440 = buffer.data(smg1 + 440);
    const auto *smg1_444 = buffer.data(smg1 + 444);
    const auto *smg1_447 = buffer.data(smg1 + 447);
    const auto *smg1_448 = buffer.data(smg1 + 448);
    const auto *smg1_449 = buffer.data(smg1 + 449);
    const auto *smg1_450 = buffer.data(smg1 + 450);
    const auto *smg1_453 = buffer.data(smg1 + 453);
    const auto *smg1_455 = buffer.data(smg1 + 455);
    const auto *smg1_456 = buffer.data(smg1 + 456);
    const auto *smg1_459 = buffer.data(smg1 + 459);
    const auto *smg1_460 = buffer.data(smg1 + 460);
    const auto *smg1_462 = buffer.data(smg1 + 462);
    const auto *smg1_463 = buffer.data(smg1 + 463);
    const auto *smg1_464 = buffer.data(smg1 + 464);
    const auto *smg1_465 = buffer.data(smg1 + 465);
    const auto *smg1_468 = buffer.data(smg1 + 468);
    const auto *smg1_470 = buffer.data(smg1 + 470);
    const auto *smg1_471 = buffer.data(smg1 + 471);
    const auto *smg1_474 = buffer.data(smg1 + 474);
    const auto *smg1_475 = buffer.data(smg1 + 475);
    const auto *smg1_477 = buffer.data(smg1 + 477);
    const auto *smg1_478 = buffer.data(smg1 + 478);
    const auto *smg1_479 = buffer.data(smg1 + 479);
    const auto *smg1_480 = buffer.data(smg1 + 480);
    const auto *smg1_483 = buffer.data(smg1 + 483);
    const auto *smg1_485 = buffer.data(smg1 + 485);
    const auto *smg1_486 = buffer.data(smg1 + 486);
    const auto *smg1_489 = buffer.data(smg1 + 489);
    const auto *smg1_490 = buffer.data(smg1 + 490);

    const auto *smh_597 = buffer.data(smh + 597);
    const auto *smh_602 = buffer.data(smh + 602);
    const auto *smh_603 = buffer.data(smh + 603);
    const auto *smh_604 = buffer.data(smh + 604);
    const auto *smh_605 = buffer.data(smh + 605);
    const auto *smh_606 = buffer.data(smh + 606);
    const auto *smh_607 = buffer.data(smh + 607);
    const auto *smh_608 = buffer.data(smh + 608);
    const auto *smh_609 = buffer.data(smh + 609);
    const auto *smh_611 = buffer.data(smh + 611);
    const auto *smh_612 = buffer.data(smh + 612);
    const auto *smh_614 = buffer.data(smh + 614);
    const auto *smh_615 = buffer.data(smh + 615);
    const auto *smh_618 = buffer.data(smh + 618);
    const auto *smh_623 = buffer.data(smh + 623);
    const auto *smh_624 = buffer.data(smh + 624);
    const auto *smh_625 = buffer.data(smh + 625);
    const auto *smh_626 = buffer.data(smh + 626);
    const auto *smh_627 = buffer.data(smh + 627);
    const auto *smh_628 = buffer.data(smh + 628);
    const auto *smh_629 = buffer.data(smh + 629);
    const auto *smh_630 = buffer.data(smh + 630);
    const auto *smh_632 = buffer.data(smh + 632);
    const auto *smh_633 = buffer.data(smh + 633);
    const auto *smh_635 = buffer.data(smh + 635);
    const auto *smh_636 = buffer.data(smh + 636);
    const auto *smh_639 = buffer.data(smh + 639);
    const auto *smh_640 = buffer.data(smh + 640);
    const auto *smh_642 = buffer.data(smh + 642);
    const auto *smh_644 = buffer.data(smh + 644);
    const auto *smh_645 = buffer.data(smh + 645);
    const auto *smh_646 = buffer.data(smh + 646);
    const auto *smh_647 = buffer.data(smh + 647);
    const auto *smh_648 = buffer.data(smh + 648);
    const auto *smh_649 = buffer.data(smh + 649);
    const auto *smh_650 = buffer.data(smh + 650);
    const auto *smh_651 = buffer.data(smh + 651);
    const auto *smh_653 = buffer.data(smh + 653);
    const auto *smh_654 = buffer.data(smh + 654);
    const auto *smh_656 = buffer.data(smh + 656);
    const auto *smh_657 = buffer.data(smh + 657);
    const auto *smh_660 = buffer.data(smh + 660);
    const auto *smh_661 = buffer.data(smh + 661);
    const auto *smh_663 = buffer.data(smh + 663);
    const auto *smh_665 = buffer.data(smh + 665);
    const auto *smh_666 = buffer.data(smh + 666);
    const auto *smh_667 = buffer.data(smh + 667);
    const auto *smh_668 = buffer.data(smh + 668);
    const auto *smh_669 = buffer.data(smh + 669);
    const auto *smh_670 = buffer.data(smh + 670);
    const auto *smh_671 = buffer.data(smh + 671);
    const auto *smh_672 = buffer.data(smh + 672);
    const auto *smh_674 = buffer.data(smh + 674);
    const auto *smh_675 = buffer.data(smh + 675);
    const auto *smh_677 = buffer.data(smh + 677);
    const auto *smh_678 = buffer.data(smh + 678);
    const auto *smh_681 = buffer.data(smh + 681);
    const auto *smh_682 = buffer.data(smh + 682);

#pragma omp simd aligned(t_797, t_798, t_799, t_800, pc_x, pc_y, slh_450, slh_602, slh_603, \
                         slh_604, smg0_434, smg1_434, smh_597, smh_602, smh_603, \
                         smh_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_16 * slh_450[k]
                   + f_3 * pc_y[k] * smh_597[k];

        t_798[k] = f_12 * slh_602[k]
                   + f_8 * smg0_434[k]
                   - f_9 * smg1_434[k]
                   + f_3 * pc_x[k] * smh_602[k];

        t_799[k] = f_12 * slh_603[k]
                   + f_3 * pc_x[k] * smh_603[k];

        t_800[k] = f_12 * slh_604[k]
                   + f_3 * pc_x[k] * smh_604[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, t_804, pc_x, slh_605, slh_606, slh_607, slh_608, \
                         smh_605, smh_606, smh_607, smh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_12 * slh_605[k]
                   + f_3 * pc_x[k] * smh_605[k];

        t_802[k] = f_12 * slh_606[k]
                   + f_3 * pc_x[k] * smh_606[k];

        t_803[k] = f_12 * slh_607[k]
                   + f_3 * pc_x[k] * smh_607[k];

        t_804[k] = f_12 * slh_608[k]
                   + f_3 * pc_x[k] * smh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pc_y, pc_z, slh_456, slh_458, smg0_430, \
                         smg0_432, smg1_430, smg1_432, smh_603, \
                         smh_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_16 * slh_456[k]
                   + f_1 * smg0_430[k]
                   - f_2 * smg1_430[k]
                   + f_3 * pc_y[k] * smh_603[k];

        t_806[k] = f_3 * pc_z[k] * smh_603[k];

        t_807[k] = f_16 * slh_458[k]
                   + f_4 * smg0_432[k]
                   - f_5 * smg1_432[k]
                   + f_3 * pc_y[k] * smh_605[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pc_y, pc_z, slh_459, slh_460, slh_461, \
                         smg0_433, smg0_434, smg1_433, smg1_434, smh_606, smh_607, \
                         smh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_16 * slh_459[k]
                   + f_6 * smg0_433[k]
                   - f_7 * smg1_433[k]
                   + f_3 * pc_y[k] * smh_606[k];

        t_809[k] = f_16 * slh_460[k]
                   + f_8 * smg0_434[k]
                   - f_9 * smg1_434[k]
                   + f_3 * pc_y[k] * smh_607[k];

        t_810[k] = f_16 * slh_461[k]
                   + f_3 * pc_y[k] * smh_608[k];

        t_811[k] = f_1 * smg0_434[k]
                   - f_2 * smg1_434[k]
                   + f_3 * pc_z[k] * smh_608[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pb_z, pc_y, pc_z, sli0_588, sli0_591, \
                         slh_441, slh_462, sli1_588, sli1_591, \
                         smh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = pb_z[k] * sli0_588[k]
                   - f_10 * pc_z[k] * sli1_588[k];

        t_813[k] = f_17 * slh_462[k]
                   + f_3 * pc_y[k] * smh_609[k];

        t_814[k] = f_11 * slh_441[k]
                   + f_3 * pc_z[k] * smh_609[k];

        t_815[k] = pb_z[k] * sli0_591[k]
                   - f_10 * pc_z[k] * sli1_591[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pb_z, pc_x, pc_y, pc_z, sli0_594, slh_464, \
                         slh_614, sli1_594, smg0_440, smg1_440, smh_611, \
                         smh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_17 * slh_464[k]
                   + f_3 * pc_y[k] * smh_611[k];

        t_817[k] = f_12 * slh_614[k]
                   + f_4 * smg0_440[k]
                   - f_5 * smg1_440[k]
                   + f_3 * pc_x[k] * smh_614[k];

        t_818[k] = pb_z[k] * sli0_594[k]
                   - f_10 * pc_z[k] * sli1_594[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, pc_y, pc_z, slh_444, slh_467, slh_618, \
                         smg0_444, smg1_444, smh_612, smh_614, \
                         smh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_11 * slh_444[k]
                   + f_3 * pc_z[k] * smh_612[k];

        t_820[k] = f_17 * slh_467[k]
                   + f_3 * pc_y[k] * smh_614[k];

        t_821[k] = f_12 * slh_618[k]
                   + f_6 * smg0_444[k]
                   - f_7 * smg1_444[k]
                   + f_3 * pc_x[k] * smh_618[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pb_z, pc_y, pc_z, sli0_598, sli0_600, \
                         slh_447, slh_448, slh_471, sli1_598, sli1_600, smh_615, \
                         smh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pb_z[k] * sli0_598[k]
                   - f_10 * pc_z[k] * sli1_598[k];

        t_823[k] = f_11 * slh_447[k]
                   + f_3 * pc_z[k] * smh_615[k];

        t_824[k] = pb_z[k] * sli0_600[k]
                   + f_12 * slh_448[k]
                   - f_10 * pc_z[k] * sli1_600[k];

        t_825[k] = f_17 * slh_471[k]
                   + f_3 * pc_y[k] * smh_618[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, slh_623, slh_624, slh_625, slh_626, \
                         smg0_449, smg1_449, smh_623, smh_624, smh_625, \
                         smh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_12 * slh_623[k]
                   + f_8 * smg0_449[k]
                   - f_9 * smg1_449[k]
                   + f_3 * pc_x[k] * smh_623[k];

        t_827[k] = f_12 * slh_624[k]
                   + f_3 * pc_x[k] * smh_624[k];

        t_828[k] = f_12 * slh_625[k]
                   + f_3 * pc_x[k] * smh_625[k];

        t_829[k] = f_12 * slh_626[k]
                   + f_3 * pc_x[k] * smh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pb_z, pc_x, pc_z, sli0_609, slh_627, \
                         slh_628, slh_629, sli1_609, smh_627, smh_628, \
                         smh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_12 * slh_627[k]
                   + f_3 * pc_x[k] * smh_627[k];

        t_831[k] = f_12 * slh_628[k]
                   + f_3 * pc_x[k] * smh_628[k];

        t_832[k] = f_12 * slh_629[k]
                   + f_3 * pc_x[k] * smh_629[k];

        t_833[k] = pb_z[k] * sli0_609[k]
                   - f_10 * pc_z[k] * sli1_609[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_y, pc_z, slh_456, slh_479, slh_480, smg0_447, \
                         smg0_448, smg1_447, smg1_448, smh_624, smh_626, \
                         smh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_11 * slh_456[k]
                   + f_3 * pc_z[k] * smh_624[k];

        t_835[k] = f_17 * slh_479[k]
                   + f_4 * smg0_447[k]
                   - f_5 * smg1_447[k]
                   + f_3 * pc_y[k] * smh_626[k];

        t_836[k] = f_17 * slh_480[k]
                   + f_6 * smg0_448[k]
                   - f_7 * smg1_448[k]
                   + f_3 * pc_y[k] * smh_627[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, slh_461, slh_481, slh_482, smg0_449, \
                         smg1_449, smh_628, smh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_17 * slh_481[k]
                   + f_8 * smg0_449[k]
                   - f_9 * smg1_449[k]
                   + f_3 * pc_y[k] * smh_628[k];

        t_838[k] = f_17 * slh_482[k]
                   + f_3 * pc_y[k] * smh_629[k];

        t_839[k] = f_11 * slh_461[k]
                   + f_1 * smg0_449[k]
                   - f_2 * smg1_449[k]
                   + f_3 * pc_z[k] * smh_629[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, pc_y, pc_z, slh_462, slh_483, slh_630, \
                         smg0_450, smg1_450, smh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_12 * slh_630[k]
                   + f_1 * smg0_450[k]
                   - f_2 * smg1_450[k]
                   + f_3 * pc_x[k] * smh_630[k];

        t_841[k] = f_18 * slh_483[k]
                   + f_3 * pc_y[k] * smh_630[k];

        t_842[k] = f_12 * slh_462[k]
                   + f_3 * pc_z[k] * smh_630[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, pc_y, slh_485, slh_633, slh_635, smg0_453, \
                         smg0_455, smg1_453, smg1_455, smh_632, smh_633, \
                         smh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_12 * slh_633[k]
                   + f_4 * smg0_453[k]
                   - f_5 * smg1_453[k]
                   + f_3 * pc_x[k] * smh_633[k];

        t_844[k] = f_18 * slh_485[k]
                   + f_3 * pc_y[k] * smh_632[k];

        t_845[k] = f_12 * slh_635[k]
                   + f_4 * smg0_455[k]
                   - f_5 * smg1_455[k]
                   + f_3 * pc_x[k] * smh_635[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pc_x, pc_y, pc_z, slh_465, slh_488, slh_636, \
                         smg0_456, smg1_456, smh_633, smh_635, \
                         smh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_12 * slh_636[k]
                   + f_6 * smg0_456[k]
                   - f_7 * smg1_456[k]
                   + f_3 * pc_x[k] * smh_636[k];

        t_847[k] = f_12 * slh_465[k]
                   + f_3 * pc_z[k] * smh_633[k];

        t_848[k] = f_18 * slh_488[k]
                   + f_3 * pc_y[k] * smh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, slh_468, slh_639, slh_640, smg0_459, \
                         smg0_460, smg1_459, smg1_460, smh_636, smh_639, \
                         smh_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_12 * slh_639[k]
                   + f_6 * smg0_459[k]
                   - f_7 * smg1_459[k]
                   + f_3 * pc_x[k] * smh_639[k];

        t_850[k] = f_12 * slh_640[k]
                   + f_8 * smg0_460[k]
                   - f_9 * smg1_460[k]
                   + f_3 * pc_x[k] * smh_640[k];

        t_851[k] = f_12 * slh_468[k]
                   + f_3 * pc_z[k] * smh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, pc_x, pc_y, slh_492, slh_642, slh_644, smg0_462, \
                         smg0_464, smg1_462, smg1_464, smh_639, smh_642, \
                         smh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_12 * slh_642[k]
                   + f_8 * smg0_462[k]
                   - f_9 * smg1_462[k]
                   + f_3 * pc_x[k] * smh_642[k];

        t_853[k] = f_18 * slh_492[k]
                   + f_3 * pc_y[k] * smh_639[k];

        t_854[k] = f_12 * slh_644[k]
                   + f_8 * smg0_464[k]
                   - f_9 * smg1_464[k]
                   + f_3 * pc_x[k] * smh_644[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, pc_x, slh_645, slh_646, slh_647, \
                         slh_648, slh_649, smh_645, smh_646, smh_647, smh_648, \
                         smh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_12 * slh_645[k]
                   + f_3 * pc_x[k] * smh_645[k];

        t_856[k] = f_12 * slh_646[k]
                   + f_3 * pc_x[k] * smh_646[k];

        t_857[k] = f_12 * slh_647[k]
                   + f_3 * pc_x[k] * smh_647[k];

        t_858[k] = f_12 * slh_648[k]
                   + f_3 * pc_x[k] * smh_648[k];

        t_859[k] = f_12 * slh_649[k]
                   + f_3 * pc_x[k] * smh_649[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pc_x, pc_y, pc_z, slh_477, slh_498, slh_650, \
                         smg0_460, smg1_460, smh_645, smh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_12 * slh_650[k]
                   + f_3 * pc_x[k] * smh_650[k];

        t_861[k] = f_18 * slh_498[k]
                   + f_1 * smg0_460[k]
                   - f_2 * smg1_460[k]
                   + f_3 * pc_y[k] * smh_645[k];

        t_862[k] = f_12 * slh_477[k]
                   + f_3 * pc_z[k] * smh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_y, slh_500, slh_501, slh_502, smg0_462, \
                         smg0_463, smg0_464, smg1_462, smg1_463, smg1_464, smh_647, smh_648, \
                         smh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_18 * slh_500[k]
                   + f_4 * smg0_462[k]
                   - f_5 * smg1_462[k]
                   + f_3 * pc_y[k] * smh_647[k];

        t_864[k] = f_18 * slh_501[k]
                   + f_6 * smg0_463[k]
                   - f_7 * smg1_463[k]
                   + f_3 * pc_y[k] * smh_648[k];

        t_865[k] = f_18 * slh_502[k]
                   + f_8 * smg0_464[k]
                   - f_9 * smg1_464[k]
                   + f_3 * pc_y[k] * smh_649[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, slh_482, slh_503, slh_651, \
                         smg0_464, smg0_465, smg1_464, smg1_465, smh_650, \
                         smh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * slh_503[k]
                   + f_3 * pc_y[k] * smh_650[k];

        t_867[k] = f_12 * slh_482[k]
                   + f_1 * smg0_464[k]
                   - f_2 * smg1_464[k]
                   + f_3 * pc_z[k] * smh_650[k];

        t_868[k] = f_12 * slh_651[k]
                   + f_1 * smg0_465[k]
                   - f_2 * smg1_465[k]
                   + f_3 * pc_x[k] * smh_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, pc_z, slh_483, slh_504, \
                         slh_506, slh_654, smg0_468, smg1_468, smh_651, smh_653, \
                         smh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * slh_504[k]
                   + f_3 * pc_y[k] * smh_651[k];

        t_870[k] = f_13 * slh_483[k]
                   + f_3 * pc_z[k] * smh_651[k];

        t_871[k] = f_12 * slh_654[k]
                   + f_4 * smg0_468[k]
                   - f_5 * smg1_468[k]
                   + f_3 * pc_x[k] * smh_654[k];

        t_872[k] = f_14 * slh_506[k]
                   + f_3 * pc_y[k] * smh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_z, slh_486, slh_656, slh_657, smg0_470, \
                         smg0_471, smg1_470, smg1_471, smh_654, smh_656, \
                         smh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_12 * slh_656[k]
                   + f_4 * smg0_470[k]
                   - f_5 * smg1_470[k]
                   + f_3 * pc_x[k] * smh_656[k];

        t_874[k] = f_12 * slh_657[k]
                   + f_6 * smg0_471[k]
                   - f_7 * smg1_471[k]
                   + f_3 * pc_x[k] * smh_657[k];

        t_875[k] = f_13 * slh_486[k]
                   + f_3 * pc_z[k] * smh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_y, slh_509, slh_660, slh_661, smg0_474, \
                         smg0_475, smg1_474, smg1_475, smh_656, smh_660, \
                         smh_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * slh_509[k]
                   + f_3 * pc_y[k] * smh_656[k];

        t_877[k] = f_12 * slh_660[k]
                   + f_6 * smg0_474[k]
                   - f_7 * smg1_474[k]
                   + f_3 * pc_x[k] * smh_660[k];

        t_878[k] = f_12 * slh_661[k]
                   + f_8 * smg0_475[k]
                   - f_9 * smg1_475[k]
                   + f_3 * pc_x[k] * smh_661[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, slh_489, slh_513, slh_663, \
                         smg0_477, smg1_477, smh_657, smh_660, \
                         smh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * slh_489[k]
                   + f_3 * pc_z[k] * smh_657[k];

        t_880[k] = f_12 * slh_663[k]
                   + f_8 * smg0_477[k]
                   - f_9 * smg1_477[k]
                   + f_3 * pc_x[k] * smh_663[k];

        t_881[k] = f_14 * slh_513[k]
                   + f_3 * pc_y[k] * smh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pc_x, slh_665, slh_666, slh_667, slh_668, \
                         smg0_479, smg1_479, smh_665, smh_666, smh_667, \
                         smh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_12 * slh_665[k]
                   + f_8 * smg0_479[k]
                   - f_9 * smg1_479[k]
                   + f_3 * pc_x[k] * smh_665[k];

        t_883[k] = f_12 * slh_666[k]
                   + f_3 * pc_x[k] * smh_666[k];

        t_884[k] = f_12 * slh_667[k]
                   + f_3 * pc_x[k] * smh_667[k];

        t_885[k] = f_12 * slh_668[k]
                   + f_3 * pc_x[k] * smh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pc_x, pc_y, slh_519, slh_669, slh_670, \
                         slh_671, smg0_475, smg1_475, smh_666, smh_669, smh_670, \
                         smh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_12 * slh_669[k]
                   + f_3 * pc_x[k] * smh_669[k];

        t_887[k] = f_12 * slh_670[k]
                   + f_3 * pc_x[k] * smh_670[k];

        t_888[k] = f_12 * slh_671[k]
                   + f_3 * pc_x[k] * smh_671[k];

        t_889[k] = f_14 * slh_519[k]
                   + f_1 * smg0_475[k]
                   - f_2 * smg1_475[k]
                   + f_3 * pc_y[k] * smh_666[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pc_y, pc_z, slh_498, slh_521, slh_522, smg0_477, \
                         smg0_478, smg1_477, smg1_478, smh_666, smh_668, \
                         smh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * slh_498[k]
                   + f_3 * pc_z[k] * smh_666[k];

        t_891[k] = f_14 * slh_521[k]
                   + f_4 * smg0_477[k]
                   - f_5 * smg1_477[k]
                   + f_3 * pc_y[k] * smh_668[k];

        t_892[k] = f_14 * slh_522[k]
                   + f_6 * smg0_478[k]
                   - f_7 * smg1_478[k]
                   + f_3 * pc_y[k] * smh_669[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, pc_z, slh_503, slh_523, slh_524, smg0_479, \
                         smg1_479, smh_670, smh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_14 * slh_523[k]
                   + f_8 * smg0_479[k]
                   - f_9 * smg1_479[k]
                   + f_3 * pc_y[k] * smh_670[k];

        t_894[k] = f_14 * slh_524[k]
                   + f_3 * pc_y[k] * smh_671[k];

        t_895[k] = f_13 * slh_503[k]
                   + f_1 * smg0_479[k]
                   - f_2 * smg1_479[k]
                   + f_3 * pc_z[k] * smh_671[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_x, pc_y, pc_z, slh_504, slh_525, slh_672, \
                         smg0_480, smg1_480, smh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_12 * slh_672[k]
                   + f_1 * smg0_480[k]
                   - f_2 * smg1_480[k]
                   + f_3 * pc_x[k] * smh_672[k];

        t_897[k] = f_13 * slh_525[k]
                   + f_3 * pc_y[k] * smh_672[k];

        t_898[k] = f_14 * slh_504[k]
                   + f_3 * pc_z[k] * smh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, slh_527, slh_675, slh_677, smg0_483, \
                         smg0_485, smg1_483, smg1_485, smh_674, smh_675, \
                         smh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_12 * slh_675[k]
                   + f_4 * smg0_483[k]
                   - f_5 * smg1_483[k]
                   + f_3 * pc_x[k] * smh_675[k];

        t_900[k] = f_13 * slh_527[k]
                   + f_3 * pc_y[k] * smh_674[k];

        t_901[k] = f_12 * slh_677[k]
                   + f_4 * smg0_485[k]
                   - f_5 * smg1_485[k]
                   + f_3 * pc_x[k] * smh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, slh_507, slh_530, slh_678, \
                         smg0_486, smg1_486, smh_675, smh_677, \
                         smh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_12 * slh_678[k]
                   + f_6 * smg0_486[k]
                   - f_7 * smg1_486[k]
                   + f_3 * pc_x[k] * smh_678[k];

        t_903[k] = f_14 * slh_507[k]
                   + f_3 * pc_z[k] * smh_675[k];

        t_904[k] = f_13 * slh_530[k]
                   + f_3 * pc_y[k] * smh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, slh_510, slh_681, slh_682, smg0_489, \
                         smg0_490, smg1_489, smg1_490, smh_678, smh_681, \
                         smh_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_12 * slh_681[k]
                   + f_6 * smg0_489[k]
                   - f_7 * smg1_489[k]
                   + f_3 * pc_x[k] * smh_681[k];

        t_906[k] = f_12 * slh_682[k]
                   + f_8 * smg0_490[k]
                   - f_9 * smg1_490[k]
                   + f_3 * pc_x[k] * smh_682[k];

        t_907[k] = f_14 * slh_510[k]
                   + f_3 * pc_z[k] * smh_678[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smg0, const size_t smg1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_756 = buffer.data(sli0 + 756);
    const auto *sli0_759 = buffer.data(sli0 + 759);
    const auto *sli0_761 = buffer.data(sli0 + 761);
    const auto *sli0_762 = buffer.data(sli0 + 762);
    const auto *sli0_765 = buffer.data(sli0 + 765);
    const auto *sli0_766 = buffer.data(sli0 + 766);
    const auto *sli0_768 = buffer.data(sli0 + 768);
    const auto *sli0_770 = buffer.data(sli0 + 770);
    const auto *sli0_783 = buffer.data(sli0 + 783);
    const auto *sli0_1008 = buffer.data(sli0 + 1008);
    const auto *sli0_1011 = buffer.data(sli0 + 1011);
    const auto *sli0_1013 = buffer.data(sli0 + 1013);
    const auto *sli0_1014 = buffer.data(sli0 + 1014);
    const auto *sli0_1017 = buffer.data(sli0 + 1017);
    const auto *sli0_1018 = buffer.data(sli0 + 1018);
    const auto *sli0_1020 = buffer.data(sli0 + 1020);

    const auto *slh_519 = buffer.data(slh + 519);
    const auto *slh_524 = buffer.data(slh + 524);
    const auto *slh_525 = buffer.data(slh + 525);
    const auto *slh_528 = buffer.data(slh + 528);
    const auto *slh_531 = buffer.data(slh + 531);
    const auto *slh_534 = buffer.data(slh + 534);
    const auto *slh_540 = buffer.data(slh + 540);
    const auto *slh_542 = buffer.data(slh + 542);
    const auto *slh_543 = buffer.data(slh + 543);
    const auto *slh_544 = buffer.data(slh + 544);
    const auto *slh_545 = buffer.data(slh + 545);
    const auto *slh_546 = buffer.data(slh + 546);
    const auto *slh_548 = buffer.data(slh + 548);
    const auto *slh_549 = buffer.data(slh + 549);
    const auto *slh_551 = buffer.data(slh + 551);
    const auto *slh_552 = buffer.data(slh + 552);
    const auto *slh_555 = buffer.data(slh + 555);
    const auto *slh_561 = buffer.data(slh + 561);
    const auto *slh_563 = buffer.data(slh + 563);
    const auto *slh_564 = buffer.data(slh + 564);
    const auto *slh_565 = buffer.data(slh + 565);
    const auto *slh_566 = buffer.data(slh + 566);
    const auto *slh_567 = buffer.data(slh + 567);
    const auto *slh_568 = buffer.data(slh + 568);
    const auto *slh_569 = buffer.data(slh + 569);
    const auto *slh_570 = buffer.data(slh + 570);
    const auto *slh_572 = buffer.data(slh + 572);
    const auto *slh_573 = buffer.data(slh + 573);
    const auto *slh_575 = buffer.data(slh + 575);
    const auto *slh_576 = buffer.data(slh + 576);
    const auto *slh_582 = buffer.data(slh + 582);
    const auto *slh_584 = buffer.data(slh + 584);
    const auto *slh_585 = buffer.data(slh + 585);
    const auto *slh_586 = buffer.data(slh + 586);
    const auto *slh_587 = buffer.data(slh + 587);
    const auto *slh_588 = buffer.data(slh + 588);
    const auto *slh_590 = buffer.data(slh + 590);
    const auto *slh_593 = buffer.data(slh + 593);
    const auto *slh_684 = buffer.data(slh + 684);
    const auto *slh_686 = buffer.data(slh + 686);
    const auto *slh_687 = buffer.data(slh + 687);
    const auto *slh_688 = buffer.data(slh + 688);
    const auto *slh_689 = buffer.data(slh + 689);
    const auto *slh_690 = buffer.data(slh + 690);
    const auto *slh_691 = buffer.data(slh + 691);
    const auto *slh_692 = buffer.data(slh + 692);
    const auto *slh_693 = buffer.data(slh + 693);
    const auto *slh_696 = buffer.data(slh + 696);
    const auto *slh_698 = buffer.data(slh + 698);
    const auto *slh_699 = buffer.data(slh + 699);
    const auto *slh_702 = buffer.data(slh + 702);
    const auto *slh_703 = buffer.data(slh + 703);
    const auto *slh_705 = buffer.data(slh + 705);
    const auto *slh_707 = buffer.data(slh + 707);
    const auto *slh_708 = buffer.data(slh + 708);
    const auto *slh_709 = buffer.data(slh + 709);
    const auto *slh_710 = buffer.data(slh + 710);
    const auto *slh_711 = buffer.data(slh + 711);
    const auto *slh_712 = buffer.data(slh + 712);
    const auto *slh_713 = buffer.data(slh + 713);
    const auto *slh_729 = buffer.data(slh + 729);
    const auto *slh_730 = buffer.data(slh + 730);
    const auto *slh_731 = buffer.data(slh + 731);
    const auto *slh_732 = buffer.data(slh + 732);
    const auto *slh_733 = buffer.data(slh + 733);
    const auto *slh_734 = buffer.data(slh + 734);
    const auto *slh_735 = buffer.data(slh + 735);
    const auto *slh_738 = buffer.data(slh + 738);
    const auto *slh_740 = buffer.data(slh + 740);
    const auto *slh_741 = buffer.data(slh + 741);
    const auto *slh_744 = buffer.data(slh + 744);
    const auto *slh_745 = buffer.data(slh + 745);
    const auto *slh_747 = buffer.data(slh + 747);
    const auto *slh_749 = buffer.data(slh + 749);
    const auto *slh_750 = buffer.data(slh + 750);
    const auto *slh_751 = buffer.data(slh + 751);
    const auto *slh_752 = buffer.data(slh + 752);
    const auto *slh_753 = buffer.data(slh + 753);
    const auto *slh_754 = buffer.data(slh + 754);
    const auto *slh_755 = buffer.data(slh + 755);
    const auto *slh_756 = buffer.data(slh + 756);
    const auto *slh_759 = buffer.data(slh + 759);
    const auto *slh_761 = buffer.data(slh + 761);
    const auto *slh_762 = buffer.data(slh + 762);
    const auto *slh_765 = buffer.data(slh + 765);
    const auto *slh_766 = buffer.data(slh + 766);
    const auto *slh_768 = buffer.data(slh + 768);

    const auto *sli1_756 = buffer.data(sli1 + 756);
    const auto *sli1_759 = buffer.data(sli1 + 759);
    const auto *sli1_761 = buffer.data(sli1 + 761);
    const auto *sli1_762 = buffer.data(sli1 + 762);
    const auto *sli1_765 = buffer.data(sli1 + 765);
    const auto *sli1_766 = buffer.data(sli1 + 766);
    const auto *sli1_768 = buffer.data(sli1 + 768);
    const auto *sli1_770 = buffer.data(sli1 + 770);
    const auto *sli1_783 = buffer.data(sli1 + 783);
    const auto *sli1_1008 = buffer.data(sli1 + 1008);
    const auto *sli1_1011 = buffer.data(sli1 + 1011);
    const auto *sli1_1013 = buffer.data(sli1 + 1013);
    const auto *sli1_1014 = buffer.data(sli1 + 1014);
    const auto *sli1_1017 = buffer.data(sli1 + 1017);
    const auto *sli1_1018 = buffer.data(sli1 + 1018);
    const auto *sli1_1020 = buffer.data(sli1 + 1020);

    const auto *smg0_490 = buffer.data(smg0 + 490);
    const auto *smg0_492 = buffer.data(smg0 + 492);
    const auto *smg0_493 = buffer.data(smg0 + 493);
    const auto *smg0_494 = buffer.data(smg0 + 494);
    const auto *smg0_495 = buffer.data(smg0 + 495);
    const auto *smg0_498 = buffer.data(smg0 + 498);
    const auto *smg0_500 = buffer.data(smg0 + 500);
    const auto *smg0_501 = buffer.data(smg0 + 501);
    const auto *smg0_504 = buffer.data(smg0 + 504);
    const auto *smg0_505 = buffer.data(smg0 + 505);
    const auto *smg0_507 = buffer.data(smg0 + 507);
    const auto *smg0_508 = buffer.data(smg0 + 508);
    const auto *smg0_509 = buffer.data(smg0 + 509);
    const auto *smg0_520 = buffer.data(smg0 + 520);
    const auto *smg0_522 = buffer.data(smg0 + 522);
    const auto *smg0_523 = buffer.data(smg0 + 523);
    const auto *smg0_524 = buffer.data(smg0 + 524);
    const auto *smg0_525 = buffer.data(smg0 + 525);
    const auto *smg0_528 = buffer.data(smg0 + 528);
    const auto *smg0_530 = buffer.data(smg0 + 530);
    const auto *smg0_531 = buffer.data(smg0 + 531);
    const auto *smg0_534 = buffer.data(smg0 + 534);
    const auto *smg0_535 = buffer.data(smg0 + 535);
    const auto *smg0_537 = buffer.data(smg0 + 537);
    const auto *smg0_538 = buffer.data(smg0 + 538);
    const auto *smg0_539 = buffer.data(smg0 + 539);

    const auto *smg1_490 = buffer.data(smg1 + 490);
    const auto *smg1_492 = buffer.data(smg1 + 492);
    const auto *smg1_493 = buffer.data(smg1 + 493);
    const auto *smg1_494 = buffer.data(smg1 + 494);
    const auto *smg1_495 = buffer.data(smg1 + 495);
    const auto *smg1_498 = buffer.data(smg1 + 498);
    const auto *smg1_500 = buffer.data(smg1 + 500);
    const auto *smg1_501 = buffer.data(smg1 + 501);
    const auto *smg1_504 = buffer.data(smg1 + 504);
    const auto *smg1_505 = buffer.data(smg1 + 505);
    const auto *smg1_507 = buffer.data(smg1 + 507);
    const auto *smg1_508 = buffer.data(smg1 + 508);
    const auto *smg1_509 = buffer.data(smg1 + 509);
    const auto *smg1_520 = buffer.data(smg1 + 520);
    const auto *smg1_522 = buffer.data(smg1 + 522);
    const auto *smg1_523 = buffer.data(smg1 + 523);
    const auto *smg1_524 = buffer.data(smg1 + 524);
    const auto *smg1_525 = buffer.data(smg1 + 525);
    const auto *smg1_528 = buffer.data(smg1 + 528);
    const auto *smg1_530 = buffer.data(smg1 + 530);
    const auto *smg1_531 = buffer.data(smg1 + 531);
    const auto *smg1_534 = buffer.data(smg1 + 534);
    const auto *smg1_535 = buffer.data(smg1 + 535);
    const auto *smg1_537 = buffer.data(smg1 + 537);
    const auto *smg1_538 = buffer.data(smg1 + 538);
    const auto *smg1_539 = buffer.data(smg1 + 539);

    const auto *smh_681 = buffer.data(smh + 681);
    const auto *smh_684 = buffer.data(smh + 684);
    const auto *smh_686 = buffer.data(smh + 686);
    const auto *smh_687 = buffer.data(smh + 687);
    const auto *smh_688 = buffer.data(smh + 688);
    const auto *smh_689 = buffer.data(smh + 689);
    const auto *smh_690 = buffer.data(smh + 690);
    const auto *smh_691 = buffer.data(smh + 691);
    const auto *smh_692 = buffer.data(smh + 692);
    const auto *smh_693 = buffer.data(smh + 693);
    const auto *smh_695 = buffer.data(smh + 695);
    const auto *smh_696 = buffer.data(smh + 696);
    const auto *smh_698 = buffer.data(smh + 698);
    const auto *smh_699 = buffer.data(smh + 699);
    const auto *smh_702 = buffer.data(smh + 702);
    const auto *smh_703 = buffer.data(smh + 703);
    const auto *smh_705 = buffer.data(smh + 705);
    const auto *smh_707 = buffer.data(smh + 707);
    const auto *smh_708 = buffer.data(smh + 708);
    const auto *smh_709 = buffer.data(smh + 709);
    const auto *smh_710 = buffer.data(smh + 710);
    const auto *smh_711 = buffer.data(smh + 711);
    const auto *smh_712 = buffer.data(smh + 712);
    const auto *smh_713 = buffer.data(smh + 713);
    const auto *smh_714 = buffer.data(smh + 714);
    const auto *smh_716 = buffer.data(smh + 716);
    const auto *smh_717 = buffer.data(smh + 717);
    const auto *smh_719 = buffer.data(smh + 719);
    const auto *smh_720 = buffer.data(smh + 720);
    const auto *smh_723 = buffer.data(smh + 723);
    const auto *smh_729 = buffer.data(smh + 729);
    const auto *smh_730 = buffer.data(smh + 730);
    const auto *smh_731 = buffer.data(smh + 731);
    const auto *smh_732 = buffer.data(smh + 732);
    const auto *smh_733 = buffer.data(smh + 733);
    const auto *smh_734 = buffer.data(smh + 734);
    const auto *smh_735 = buffer.data(smh + 735);
    const auto *smh_737 = buffer.data(smh + 737);
    const auto *smh_738 = buffer.data(smh + 738);
    const auto *smh_740 = buffer.data(smh + 740);
    const auto *smh_741 = buffer.data(smh + 741);
    const auto *smh_744 = buffer.data(smh + 744);
    const auto *smh_745 = buffer.data(smh + 745);
    const auto *smh_747 = buffer.data(smh + 747);
    const auto *smh_749 = buffer.data(smh + 749);
    const auto *smh_750 = buffer.data(smh + 750);
    const auto *smh_751 = buffer.data(smh + 751);
    const auto *smh_752 = buffer.data(smh + 752);
    const auto *smh_753 = buffer.data(smh + 753);
    const auto *smh_754 = buffer.data(smh + 754);
    const auto *smh_755 = buffer.data(smh + 755);
    const auto *smh_756 = buffer.data(smh + 756);
    const auto *smh_758 = buffer.data(smh + 758);
    const auto *smh_759 = buffer.data(smh + 759);
    const auto *smh_761 = buffer.data(smh + 761);
    const auto *smh_762 = buffer.data(smh + 762);

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, slh_534, slh_684, slh_686, smg0_492, \
                         smg0_494, smg1_492, smg1_494, smh_681, smh_684, \
                         smh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_12 * slh_684[k]
                   + f_8 * smg0_492[k]
                   - f_9 * smg1_492[k]
                   + f_3 * pc_x[k] * smh_684[k];

        t_909[k] = f_13 * slh_534[k]
                   + f_3 * pc_y[k] * smh_681[k];

        t_910[k] = f_12 * slh_686[k]
                   + f_8 * smg0_494[k]
                   - f_9 * smg1_494[k]
                   + f_3 * pc_x[k] * smh_686[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pc_x, slh_687, slh_688, slh_689, \
                         slh_690, slh_691, smh_687, smh_688, smh_689, smh_690, \
                         smh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_12 * slh_687[k]
                   + f_3 * pc_x[k] * smh_687[k];

        t_912[k] = f_12 * slh_688[k]
                   + f_3 * pc_x[k] * smh_688[k];

        t_913[k] = f_12 * slh_689[k]
                   + f_3 * pc_x[k] * smh_689[k];

        t_914[k] = f_12 * slh_690[k]
                   + f_3 * pc_x[k] * smh_690[k];

        t_915[k] = f_12 * slh_691[k]
                   + f_3 * pc_x[k] * smh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pc_x, pc_y, pc_z, slh_519, slh_540, slh_692, \
                         smg0_490, smg1_490, smh_687, smh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_12 * slh_692[k]
                   + f_3 * pc_x[k] * smh_692[k];

        t_917[k] = f_13 * slh_540[k]
                   + f_1 * smg0_490[k]
                   - f_2 * smg1_490[k]
                   + f_3 * pc_y[k] * smh_687[k];

        t_918[k] = f_14 * slh_519[k]
                   + f_3 * pc_z[k] * smh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pc_y, slh_542, slh_543, slh_544, smg0_492, \
                         smg0_493, smg0_494, smg1_492, smg1_493, smg1_494, smh_689, smh_690, \
                         smh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_13 * slh_542[k]
                   + f_4 * smg0_492[k]
                   - f_5 * smg1_492[k]
                   + f_3 * pc_y[k] * smh_689[k];

        t_920[k] = f_13 * slh_543[k]
                   + f_6 * smg0_493[k]
                   - f_7 * smg1_493[k]
                   + f_3 * pc_y[k] * smh_690[k];

        t_921[k] = f_13 * slh_544[k]
                   + f_8 * smg0_494[k]
                   - f_9 * smg1_494[k]
                   + f_3 * pc_y[k] * smh_691[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pc_x, pc_y, pc_z, slh_524, slh_545, slh_693, \
                         smg0_494, smg0_495, smg1_494, smg1_495, smh_692, \
                         smh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * slh_545[k]
                   + f_3 * pc_y[k] * smh_692[k];

        t_923[k] = f_14 * slh_524[k]
                   + f_1 * smg0_494[k]
                   - f_2 * smg1_494[k]
                   + f_3 * pc_z[k] * smh_692[k];

        t_924[k] = f_12 * slh_693[k]
                   + f_1 * smg0_495[k]
                   - f_2 * smg1_495[k]
                   + f_3 * pc_x[k] * smh_693[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, slh_525, slh_546, \
                         slh_548, slh_696, smg0_498, smg1_498, smh_693, smh_695, \
                         smh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_12 * slh_546[k]
                   + f_3 * pc_y[k] * smh_693[k];

        t_926[k] = f_18 * slh_525[k]
                   + f_3 * pc_z[k] * smh_693[k];

        t_927[k] = f_12 * slh_696[k]
                   + f_4 * smg0_498[k]
                   - f_5 * smg1_498[k]
                   + f_3 * pc_x[k] * smh_696[k];

        t_928[k] = f_12 * slh_548[k]
                   + f_3 * pc_y[k] * smh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pc_x, pc_z, slh_528, slh_698, slh_699, smg0_500, \
                         smg0_501, smg1_500, smg1_501, smh_696, smh_698, \
                         smh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_12 * slh_698[k]
                   + f_4 * smg0_500[k]
                   - f_5 * smg1_500[k]
                   + f_3 * pc_x[k] * smh_698[k];

        t_930[k] = f_12 * slh_699[k]
                   + f_6 * smg0_501[k]
                   - f_7 * smg1_501[k]
                   + f_3 * pc_x[k] * smh_699[k];

        t_931[k] = f_18 * slh_528[k]
                   + f_3 * pc_z[k] * smh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pc_x, pc_y, slh_551, slh_702, slh_703, smg0_504, \
                         smg0_505, smg1_504, smg1_505, smh_698, smh_702, \
                         smh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * slh_551[k]
                   + f_3 * pc_y[k] * smh_698[k];

        t_933[k] = f_12 * slh_702[k]
                   + f_6 * smg0_504[k]
                   - f_7 * smg1_504[k]
                   + f_3 * pc_x[k] * smh_702[k];

        t_934[k] = f_12 * slh_703[k]
                   + f_8 * smg0_505[k]
                   - f_9 * smg1_505[k]
                   + f_3 * pc_x[k] * smh_703[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pc_x, pc_y, pc_z, slh_531, slh_555, slh_705, \
                         smg0_507, smg1_507, smh_699, smh_702, \
                         smh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_18 * slh_531[k]
                   + f_3 * pc_z[k] * smh_699[k];

        t_936[k] = f_12 * slh_705[k]
                   + f_8 * smg0_507[k]
                   - f_9 * smg1_507[k]
                   + f_3 * pc_x[k] * smh_705[k];

        t_937[k] = f_12 * slh_555[k]
                   + f_3 * pc_y[k] * smh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pc_x, slh_707, slh_708, slh_709, slh_710, \
                         smg0_509, smg1_509, smh_707, smh_708, smh_709, \
                         smh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_12 * slh_707[k]
                   + f_8 * smg0_509[k]
                   - f_9 * smg1_509[k]
                   + f_3 * pc_x[k] * smh_707[k];

        t_939[k] = f_12 * slh_708[k]
                   + f_3 * pc_x[k] * smh_708[k];

        t_940[k] = f_12 * slh_709[k]
                   + f_3 * pc_x[k] * smh_709[k];

        t_941[k] = f_12 * slh_710[k]
                   + f_3 * pc_x[k] * smh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, slh_561, slh_711, slh_712, \
                         slh_713, smg0_505, smg1_505, smh_708, smh_711, smh_712, \
                         smh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_12 * slh_711[k]
                   + f_3 * pc_x[k] * smh_711[k];

        t_943[k] = f_12 * slh_712[k]
                   + f_3 * pc_x[k] * smh_712[k];

        t_944[k] = f_12 * slh_713[k]
                   + f_3 * pc_x[k] * smh_713[k];

        t_945[k] = f_12 * slh_561[k]
                   + f_1 * smg0_505[k]
                   - f_2 * smg1_505[k]
                   + f_3 * pc_y[k] * smh_708[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pc_y, pc_z, slh_540, slh_563, slh_564, smg0_507, \
                         smg0_508, smg1_507, smg1_508, smh_708, smh_710, \
                         smh_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_18 * slh_540[k]
                   + f_3 * pc_z[k] * smh_708[k];

        t_947[k] = f_12 * slh_563[k]
                   + f_4 * smg0_507[k]
                   - f_5 * smg1_507[k]
                   + f_3 * pc_y[k] * smh_710[k];

        t_948[k] = f_12 * slh_564[k]
                   + f_6 * smg0_508[k]
                   - f_7 * smg1_508[k]
                   + f_3 * pc_y[k] * smh_711[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pb_y, pc_y, pc_z, sli0_756, slh_545, \
                         slh_565, slh_566, sli1_756, smg0_509, smg1_509, smh_712, \
                         smh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_12 * slh_565[k]
                   + f_8 * smg0_509[k]
                   - f_9 * smg1_509[k]
                   + f_3 * pc_y[k] * smh_712[k];

        t_950[k] = f_12 * slh_566[k]
                   + f_3 * pc_y[k] * smh_713[k];

        t_951[k] = f_18 * slh_545[k]
                   + f_1 * smg0_509[k]
                   - f_2 * smg1_509[k]
                   + f_3 * pc_z[k] * smh_713[k];

        t_952[k] = pb_y[k] * sli0_756[k]
                   - f_10 * pc_y[k] * sli1_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pb_y, pc_y, pc_z, sli0_759, slh_546, \
                         slh_567, slh_568, slh_569, sli1_759, smh_714, \
                         smh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_11 * slh_567[k]
                   + f_3 * pc_y[k] * smh_714[k];

        t_954[k] = f_17 * slh_546[k]
                   + f_3 * pc_z[k] * smh_714[k];

        t_955[k] = pb_y[k] * sli0_759[k]
                   + f_12 * slh_568[k]
                   - f_10 * pc_y[k] * sli1_759[k];

        t_956[k] = f_11 * slh_569[k]
                   + f_3 * pc_y[k] * smh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pb_y, pc_y, pc_z, sli0_761, sli0_762, \
                         slh_549, slh_570, slh_572, sli1_761, sli1_762, smh_717, \
                         smh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pb_y[k] * sli0_761[k]
                   - f_10 * pc_y[k] * sli1_761[k];

        t_958[k] = pb_y[k] * sli0_762[k]
                   + f_13 * slh_570[k]
                   - f_10 * pc_y[k] * sli1_762[k];

        t_959[k] = f_17 * slh_549[k]
                   + f_3 * pc_z[k] * smh_717[k];

        t_960[k] = f_11 * slh_572[k]
                   + f_3 * pc_y[k] * smh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pb_y, pc_y, pc_z, sli0_765, sli0_766, slh_552, \
                         slh_573, sli1_765, sli1_766, smh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pb_y[k] * sli0_765[k]
                   - f_10 * pc_y[k] * sli1_765[k];

        t_962[k] = pb_y[k] * sli0_766[k]
                   + f_14 * slh_573[k]
                   - f_10 * pc_y[k] * sli1_766[k];

        t_963[k] = f_17 * slh_552[k]
                   + f_3 * pc_z[k] * smh_720[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pb_y, pc_x, pc_y, sli0_768, sli0_770, \
                         slh_575, slh_576, slh_729, sli1_768, sli1_770, smh_723, \
                         smh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pb_y[k] * sli0_768[k]
                   + f_12 * slh_575[k]
                   - f_10 * pc_y[k] * sli1_768[k];

        t_965[k] = f_11 * slh_576[k]
                   + f_3 * pc_y[k] * smh_723[k];

        t_966[k] = pb_y[k] * sli0_770[k]
                   - f_10 * pc_y[k] * sli1_770[k];

        t_967[k] = f_12 * slh_729[k]
                   + f_3 * pc_x[k] * smh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, slh_730, slh_731, slh_732, \
                         slh_733, slh_734, smh_730, smh_731, smh_732, smh_733, \
                         smh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_12 * slh_730[k]
                   + f_3 * pc_x[k] * smh_730[k];

        t_969[k] = f_12 * slh_731[k]
                   + f_3 * pc_x[k] * smh_731[k];

        t_970[k] = f_12 * slh_732[k]
                   + f_3 * pc_x[k] * smh_732[k];

        t_971[k] = f_12 * slh_733[k]
                   + f_3 * pc_x[k] * smh_733[k];

        t_972[k] = f_12 * slh_734[k]
                   + f_3 * pc_x[k] * smh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pc_y, pc_z, slh_561, slh_582, slh_584, smg0_520, \
                         smg0_522, smg1_520, smg1_522, smh_729, \
                         smh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_11 * slh_582[k]
                   + f_1 * smg0_520[k]
                   - f_2 * smg1_520[k]
                   + f_3 * pc_y[k] * smh_729[k];

        t_974[k] = f_17 * slh_561[k]
                   + f_3 * pc_z[k] * smh_729[k];

        t_975[k] = f_11 * slh_584[k]
                   + f_4 * smg0_522[k]
                   - f_5 * smg1_522[k]
                   + f_3 * pc_y[k] * smh_731[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_y, slh_585, slh_586, slh_587, smg0_523, \
                         smg0_524, smg1_523, smg1_524, smh_732, smh_733, \
                         smh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_11 * slh_585[k]
                   + f_6 * smg0_523[k]
                   - f_7 * smg1_523[k]
                   + f_3 * pc_y[k] * smh_732[k];

        t_977[k] = f_11 * slh_586[k]
                   + f_8 * smg0_524[k]
                   - f_9 * smg1_524[k]
                   + f_3 * pc_y[k] * smh_733[k];

        t_978[k] = f_11 * slh_587[k]
                   + f_3 * pc_y[k] * smh_734[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pb_y, pc_x, pc_y, pc_z, sli0_783, \
                         slh_567, slh_735, sli1_783, smg0_525, smg1_525, \
                         smh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pb_y[k] * sli0_783[k]
                   - f_10 * pc_y[k] * sli1_783[k];

        t_980[k] = f_12 * slh_735[k]
                   + f_1 * smg0_525[k]
                   - f_2 * smg1_525[k]
                   + f_3 * pc_x[k] * smh_735[k];

        t_981[k] = f_3 * pc_y[k] * smh_735[k];

        t_982[k] = f_16 * slh_567[k]
                   + f_3 * pc_z[k] * smh_735[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_x, pc_y, slh_738, slh_740, smg0_528, \
                         smg0_530, smg1_528, smg1_530, smh_737, smh_738, \
                         smh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_12 * slh_738[k]
                   + f_4 * smg0_528[k]
                   - f_5 * smg1_528[k]
                   + f_3 * pc_x[k] * smh_738[k];

        t_984[k] = f_3 * pc_y[k] * smh_737[k];

        t_985[k] = f_12 * slh_740[k]
                   + f_4 * smg0_530[k]
                   - f_5 * smg1_530[k]
                   + f_3 * pc_x[k] * smh_740[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_x, pc_y, pc_z, slh_570, slh_741, smg0_531, \
                         smg1_531, smh_738, smh_740, smh_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_12 * slh_741[k]
                   + f_6 * smg0_531[k]
                   - f_7 * smg1_531[k]
                   + f_3 * pc_x[k] * smh_741[k];

        t_987[k] = f_16 * slh_570[k]
                   + f_3 * pc_z[k] * smh_738[k];

        t_988[k] = f_3 * pc_y[k] * smh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_z, slh_573, slh_744, slh_745, smg0_534, \
                         smg0_535, smg1_534, smg1_535, smh_741, smh_744, \
                         smh_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_12 * slh_744[k]
                   + f_6 * smg0_534[k]
                   - f_7 * smg1_534[k]
                   + f_3 * pc_x[k] * smh_744[k];

        t_990[k] = f_12 * slh_745[k]
                   + f_8 * smg0_535[k]
                   - f_9 * smg1_535[k]
                   + f_3 * pc_x[k] * smh_745[k];

        t_991[k] = f_16 * slh_573[k]
                   + f_3 * pc_z[k] * smh_741[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, pc_x, pc_y, slh_747, slh_749, smg0_537, \
                         smg0_539, smg1_537, smg1_539, smh_744, smh_747, \
                         smh_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_12 * slh_747[k]
                   + f_8 * smg0_537[k]
                   - f_9 * smg1_537[k]
                   + f_3 * pc_x[k] * smh_747[k];

        t_993[k] = f_3 * pc_y[k] * smh_744[k];

        t_994[k] = f_12 * slh_749[k]
                   + f_8 * smg0_539[k]
                   - f_9 * smg1_539[k]
                   + f_3 * pc_x[k] * smh_749[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pc_x, slh_750, slh_751, slh_752, \
                         slh_753, slh_754, smh_750, smh_751, smh_752, smh_753, \
                         smh_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_12 * slh_750[k]
                   + f_3 * pc_x[k] * smh_750[k];

        t_996[k] = f_12 * slh_751[k]
                   + f_3 * pc_x[k] * smh_751[k];

        t_997[k] = f_12 * slh_752[k]
                   + f_3 * pc_x[k] * smh_752[k];

        t_998[k] = f_12 * slh_753[k]
                   + f_3 * pc_x[k] * smh_753[k];

        t_999[k] = f_12 * slh_754[k]
                   + f_3 * pc_x[k] * smh_754[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pc_x, pc_y, pc_z, slh_582, slh_755, \
                         smg0_535, smg0_537, smg1_535, smg1_537, smh_750, smh_752, \
                         smh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_12 * slh_755[k]
                    + f_3 * pc_x[k] * smh_755[k];

        t_1001[k] = f_1 * smg0_535[k]
                    - f_2 * smg1_535[k]
                    + f_3 * pc_y[k] * smh_750[k];

        t_1002[k] = f_16 * slh_582[k]
                    + f_3 * pc_z[k] * smh_750[k];

        t_1003[k] = f_4 * smg0_537[k]
                    - f_5 * smg1_537[k]
                    + f_3 * pc_y[k] * smh_752[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, slh_587, smg0_538, \
                         smg0_539, smg1_538, smg1_539, smh_753, smh_754, \
                         smh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_6 * smg0_538[k]
                    - f_7 * smg1_538[k]
                    + f_3 * pc_y[k] * smh_753[k];

        t_1005[k] = f_8 * smg0_539[k]
                    - f_9 * smg1_539[k]
                    + f_3 * pc_y[k] * smh_754[k];

        t_1006[k] = f_3 * pc_y[k] * smh_755[k];

        t_1007[k] = f_16 * slh_587[k]
                    + f_1 * smg0_539[k]
                    - f_2 * smg1_539[k]
                    + f_3 * pc_z[k] * smh_755[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pb_x, pc_x, pc_y, pc_z, sli0_1008, \
                         sli0_1011, slh_588, slh_756, slh_759, sli1_1008, sli1_1011, \
                         smh_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = pb_x[k] * sli0_1008[k]
                    + f_17 * slh_756[k]
                    - f_10 * pc_x[k] * sli1_1008[k];

        t_1009[k] = f_15 * slh_588[k]
                    + f_3 * pc_y[k] * smh_756[k];

        t_1010[k] = f_3 * pc_z[k] * smh_756[k];

        t_1011[k] = pb_x[k] * sli0_1011[k]
                    + f_14 * slh_759[k]
                    - f_10 * pc_x[k] * sli1_1011[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pb_x, pc_x, pc_y, sli0_1013, sli0_1014, \
                         slh_590, slh_761, slh_762, sli1_1013, sli1_1014, \
                         smh_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_15 * slh_590[k]
                    + f_3 * pc_y[k] * smh_758[k];

        t_1013[k] = pb_x[k] * sli0_1013[k]
                    + f_14 * slh_761[k]
                    - f_10 * pc_x[k] * sli1_1013[k];

        t_1014[k] = pb_x[k] * sli0_1014[k]
                    + f_13 * slh_762[k]
                    - f_10 * pc_x[k] * sli1_1014[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, pb_x, pc_x, pc_y, pc_z, sli0_1017, slh_593, \
                         slh_765, sli1_1017, smh_759, smh_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = f_3 * pc_z[k] * smh_759[k];

        t_1016[k] = f_15 * slh_593[k]
                    + f_3 * pc_y[k] * smh_761[k];

        t_1017[k] = pb_x[k] * sli0_1017[k]
                    + f_13 * slh_765[k]
                    - f_10 * pc_x[k] * sli1_1017[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, pb_x, pc_x, pc_z, sli0_1018, sli0_1020, \
                         slh_766, slh_768, sli1_1018, sli1_1020, \
                         smh_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = pb_x[k] * sli0_1018[k]
                    + f_12 * slh_766[k]
                    - f_10 * pc_x[k] * sli1_1018[k];

        t_1019[k] = f_3 * pc_z[k] * smh_762[k];

        t_1020[k] = pb_x[k] * sli0_1020[k]
                    + f_12 * slh_768[k]
                    - f_10 * pc_x[k] * sli1_1020[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sli0,
                                                          const size_t slh, const size_t sli1,
                                                          const size_t smh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_784 = buffer.data(sli0 + 784);
    const auto *sli0_787 = buffer.data(sli0 + 787);
    const auto *sli0_790 = buffer.data(sli0 + 790);
    const auto *sli0_794 = buffer.data(sli0 + 794);
    const auto *sli0_1022 = buffer.data(sli0 + 1022);
    const auto *sli0_1029 = buffer.data(sli0 + 1029);
    const auto *sli0_1031 = buffer.data(sli0 + 1031);
    const auto *sli0_1032 = buffer.data(sli0 + 1032);
    const auto *sli0_1033 = buffer.data(sli0 + 1033);
    const auto *sli0_1035 = buffer.data(sli0 + 1035);
    const auto *sli0_1041 = buffer.data(sli0 + 1041);
    const auto *sli0_1045 = buffer.data(sli0 + 1045);
    const auto *sli0_1048 = buffer.data(sli0 + 1048);
    const auto *sli0_1050 = buffer.data(sli0 + 1050);
    const auto *sli0_1057 = buffer.data(sli0 + 1057);
    const auto *sli0_1059 = buffer.data(sli0 + 1059);
    const auto *sli0_1060 = buffer.data(sli0 + 1060);
    const auto *sli0_1061 = buffer.data(sli0 + 1061);
    const auto *sli0_1063 = buffer.data(sli0 + 1063);
    const auto *sli0_1064 = buffer.data(sli0 + 1064);
    const auto *sli0_1067 = buffer.data(sli0 + 1067);
    const auto *sli0_1069 = buffer.data(sli0 + 1069);
    const auto *sli0_1070 = buffer.data(sli0 + 1070);
    const auto *sli0_1073 = buffer.data(sli0 + 1073);
    const auto *sli0_1074 = buffer.data(sli0 + 1074);
    const auto *sli0_1076 = buffer.data(sli0 + 1076);
    const auto *sli0_1078 = buffer.data(sli0 + 1078);
    const auto *sli0_1085 = buffer.data(sli0 + 1085);
    const auto *sli0_1087 = buffer.data(sli0 + 1087);
    const auto *sli0_1088 = buffer.data(sli0 + 1088);
    const auto *sli0_1089 = buffer.data(sli0 + 1089);
    const auto *sli0_1091 = buffer.data(sli0 + 1091);
    const auto *sli0_1092 = buffer.data(sli0 + 1092);
    const auto *sli0_1095 = buffer.data(sli0 + 1095);
    const auto *sli0_1097 = buffer.data(sli0 + 1097);
    const auto *sli0_1098 = buffer.data(sli0 + 1098);
    const auto *sli0_1101 = buffer.data(sli0 + 1101);
    const auto *sli0_1102 = buffer.data(sli0 + 1102);
    const auto *sli0_1104 = buffer.data(sli0 + 1104);
    const auto *sli0_1106 = buffer.data(sli0 + 1106);
    const auto *sli0_1113 = buffer.data(sli0 + 1113);
    const auto *sli0_1115 = buffer.data(sli0 + 1115);
    const auto *sli0_1116 = buffer.data(sli0 + 1116);
    const auto *sli0_1117 = buffer.data(sli0 + 1117);
    const auto *sli0_1119 = buffer.data(sli0 + 1119);
    const auto *sli0_1120 = buffer.data(sli0 + 1120);
    const auto *sli0_1123 = buffer.data(sli0 + 1123);
    const auto *sli0_1125 = buffer.data(sli0 + 1125);
    const auto *sli0_1126 = buffer.data(sli0 + 1126);
    const auto *sli0_1129 = buffer.data(sli0 + 1129);
    const auto *sli0_1130 = buffer.data(sli0 + 1130);
    const auto *sli0_1132 = buffer.data(sli0 + 1132);
    const auto *sli0_1134 = buffer.data(sli0 + 1134);
    const auto *sli0_1141 = buffer.data(sli0 + 1141);

    const auto *slh_588 = buffer.data(slh + 588);
    const auto *slh_591 = buffer.data(slh + 591);
    const auto *slh_594 = buffer.data(slh + 594);
    const auto *slh_597 = buffer.data(slh + 597);
    const auto *slh_603 = buffer.data(slh + 603);
    const auto *slh_608 = buffer.data(slh + 608);
    const auto *slh_609 = buffer.data(slh + 609);
    const auto *slh_611 = buffer.data(slh + 611);
    const auto *slh_612 = buffer.data(slh + 612);
    const auto *slh_614 = buffer.data(slh + 614);
    const auto *slh_615 = buffer.data(slh + 615);
    const auto *slh_618 = buffer.data(slh + 618);
    const auto *slh_624 = buffer.data(slh + 624);
    const auto *slh_629 = buffer.data(slh + 629);
    const auto *slh_630 = buffer.data(slh + 630);
    const auto *slh_632 = buffer.data(slh + 632);
    const auto *slh_633 = buffer.data(slh + 633);
    const auto *slh_635 = buffer.data(slh + 635);
    const auto *slh_636 = buffer.data(slh + 636);
    const auto *slh_639 = buffer.data(slh + 639);
    const auto *slh_645 = buffer.data(slh + 645);
    const auto *slh_650 = buffer.data(slh + 650);
    const auto *slh_651 = buffer.data(slh + 651);
    const auto *slh_653 = buffer.data(slh + 653);
    const auto *slh_654 = buffer.data(slh + 654);
    const auto *slh_656 = buffer.data(slh + 656);
    const auto *slh_657 = buffer.data(slh + 657);
    const auto *slh_660 = buffer.data(slh + 660);
    const auto *slh_671 = buffer.data(slh + 671);
    const auto *slh_672 = buffer.data(slh + 672);
    const auto *slh_674 = buffer.data(slh + 674);
    const auto *slh_677 = buffer.data(slh + 677);
    const auto *slh_681 = buffer.data(slh + 681);
    const auto *slh_770 = buffer.data(slh + 770);
    const auto *slh_771 = buffer.data(slh + 771);
    const auto *slh_772 = buffer.data(slh + 772);
    const auto *slh_773 = buffer.data(slh + 773);
    const auto *slh_774 = buffer.data(slh + 774);
    const auto *slh_775 = buffer.data(slh + 775);
    const auto *slh_776 = buffer.data(slh + 776);
    const auto *slh_782 = buffer.data(slh + 782);
    const auto *slh_786 = buffer.data(slh + 786);
    const auto *slh_789 = buffer.data(slh + 789);
    const auto *slh_791 = buffer.data(slh + 791);
    const auto *slh_792 = buffer.data(slh + 792);
    const auto *slh_793 = buffer.data(slh + 793);
    const auto *slh_794 = buffer.data(slh + 794);
    const auto *slh_795 = buffer.data(slh + 795);
    const auto *slh_796 = buffer.data(slh + 796);
    const auto *slh_797 = buffer.data(slh + 797);
    const auto *slh_798 = buffer.data(slh + 798);
    const auto *slh_801 = buffer.data(slh + 801);
    const auto *slh_803 = buffer.data(slh + 803);
    const auto *slh_804 = buffer.data(slh + 804);
    const auto *slh_807 = buffer.data(slh + 807);
    const auto *slh_808 = buffer.data(slh + 808);
    const auto *slh_810 = buffer.data(slh + 810);
    const auto *slh_812 = buffer.data(slh + 812);
    const auto *slh_813 = buffer.data(slh + 813);
    const auto *slh_814 = buffer.data(slh + 814);
    const auto *slh_815 = buffer.data(slh + 815);
    const auto *slh_816 = buffer.data(slh + 816);
    const auto *slh_817 = buffer.data(slh + 817);
    const auto *slh_818 = buffer.data(slh + 818);
    const auto *slh_819 = buffer.data(slh + 819);
    const auto *slh_822 = buffer.data(slh + 822);
    const auto *slh_824 = buffer.data(slh + 824);
    const auto *slh_825 = buffer.data(slh + 825);
    const auto *slh_828 = buffer.data(slh + 828);
    const auto *slh_829 = buffer.data(slh + 829);
    const auto *slh_831 = buffer.data(slh + 831);
    const auto *slh_833 = buffer.data(slh + 833);
    const auto *slh_834 = buffer.data(slh + 834);
    const auto *slh_835 = buffer.data(slh + 835);
    const auto *slh_836 = buffer.data(slh + 836);
    const auto *slh_837 = buffer.data(slh + 837);
    const auto *slh_838 = buffer.data(slh + 838);
    const auto *slh_839 = buffer.data(slh + 839);
    const auto *slh_840 = buffer.data(slh + 840);
    const auto *slh_843 = buffer.data(slh + 843);
    const auto *slh_845 = buffer.data(slh + 845);
    const auto *slh_846 = buffer.data(slh + 846);
    const auto *slh_849 = buffer.data(slh + 849);
    const auto *slh_850 = buffer.data(slh + 850);
    const auto *slh_852 = buffer.data(slh + 852);
    const auto *slh_854 = buffer.data(slh + 854);
    const auto *slh_855 = buffer.data(slh + 855);
    const auto *slh_856 = buffer.data(slh + 856);
    const auto *slh_857 = buffer.data(slh + 857);
    const auto *slh_858 = buffer.data(slh + 858);
    const auto *slh_859 = buffer.data(slh + 859);
    const auto *slh_860 = buffer.data(slh + 860);

    const auto *sli1_784 = buffer.data(sli1 + 784);
    const auto *sli1_787 = buffer.data(sli1 + 787);
    const auto *sli1_790 = buffer.data(sli1 + 790);
    const auto *sli1_794 = buffer.data(sli1 + 794);
    const auto *sli1_1022 = buffer.data(sli1 + 1022);
    const auto *sli1_1029 = buffer.data(sli1 + 1029);
    const auto *sli1_1031 = buffer.data(sli1 + 1031);
    const auto *sli1_1032 = buffer.data(sli1 + 1032);
    const auto *sli1_1033 = buffer.data(sli1 + 1033);
    const auto *sli1_1035 = buffer.data(sli1 + 1035);
    const auto *sli1_1041 = buffer.data(sli1 + 1041);
    const auto *sli1_1045 = buffer.data(sli1 + 1045);
    const auto *sli1_1048 = buffer.data(sli1 + 1048);
    const auto *sli1_1050 = buffer.data(sli1 + 1050);
    const auto *sli1_1057 = buffer.data(sli1 + 1057);
    const auto *sli1_1059 = buffer.data(sli1 + 1059);
    const auto *sli1_1060 = buffer.data(sli1 + 1060);
    const auto *sli1_1061 = buffer.data(sli1 + 1061);
    const auto *sli1_1063 = buffer.data(sli1 + 1063);
    const auto *sli1_1064 = buffer.data(sli1 + 1064);
    const auto *sli1_1067 = buffer.data(sli1 + 1067);
    const auto *sli1_1069 = buffer.data(sli1 + 1069);
    const auto *sli1_1070 = buffer.data(sli1 + 1070);
    const auto *sli1_1073 = buffer.data(sli1 + 1073);
    const auto *sli1_1074 = buffer.data(sli1 + 1074);
    const auto *sli1_1076 = buffer.data(sli1 + 1076);
    const auto *sli1_1078 = buffer.data(sli1 + 1078);
    const auto *sli1_1085 = buffer.data(sli1 + 1085);
    const auto *sli1_1087 = buffer.data(sli1 + 1087);
    const auto *sli1_1088 = buffer.data(sli1 + 1088);
    const auto *sli1_1089 = buffer.data(sli1 + 1089);
    const auto *sli1_1091 = buffer.data(sli1 + 1091);
    const auto *sli1_1092 = buffer.data(sli1 + 1092);
    const auto *sli1_1095 = buffer.data(sli1 + 1095);
    const auto *sli1_1097 = buffer.data(sli1 + 1097);
    const auto *sli1_1098 = buffer.data(sli1 + 1098);
    const auto *sli1_1101 = buffer.data(sli1 + 1101);
    const auto *sli1_1102 = buffer.data(sli1 + 1102);
    const auto *sli1_1104 = buffer.data(sli1 + 1104);
    const auto *sli1_1106 = buffer.data(sli1 + 1106);
    const auto *sli1_1113 = buffer.data(sli1 + 1113);
    const auto *sli1_1115 = buffer.data(sli1 + 1115);
    const auto *sli1_1116 = buffer.data(sli1 + 1116);
    const auto *sli1_1117 = buffer.data(sli1 + 1117);
    const auto *sli1_1119 = buffer.data(sli1 + 1119);
    const auto *sli1_1120 = buffer.data(sli1 + 1120);
    const auto *sli1_1123 = buffer.data(sli1 + 1123);
    const auto *sli1_1125 = buffer.data(sli1 + 1125);
    const auto *sli1_1126 = buffer.data(sli1 + 1126);
    const auto *sli1_1129 = buffer.data(sli1 + 1129);
    const auto *sli1_1130 = buffer.data(sli1 + 1130);
    const auto *sli1_1132 = buffer.data(sli1 + 1132);
    const auto *sli1_1134 = buffer.data(sli1 + 1134);
    const auto *sli1_1141 = buffer.data(sli1 + 1141);

    const auto *smh_765 = buffer.data(smh + 765);
    const auto *smh_771 = buffer.data(smh + 771);
    const auto *smh_772 = buffer.data(smh + 772);
    const auto *smh_773 = buffer.data(smh + 773);
    const auto *smh_774 = buffer.data(smh + 774);
    const auto *smh_775 = buffer.data(smh + 775);
    const auto *smh_776 = buffer.data(smh + 776);
    const auto *smh_777 = buffer.data(smh + 777);
    const auto *smh_779 = buffer.data(smh + 779);
    const auto *smh_780 = buffer.data(smh + 780);
    const auto *smh_782 = buffer.data(smh + 782);
    const auto *smh_783 = buffer.data(smh + 783);
    const auto *smh_786 = buffer.data(smh + 786);
    const auto *smh_792 = buffer.data(smh + 792);
    const auto *smh_793 = buffer.data(smh + 793);
    const auto *smh_794 = buffer.data(smh + 794);
    const auto *smh_795 = buffer.data(smh + 795);
    const auto *smh_796 = buffer.data(smh + 796);
    const auto *smh_797 = buffer.data(smh + 797);
    const auto *smh_798 = buffer.data(smh + 798);
    const auto *smh_800 = buffer.data(smh + 800);
    const auto *smh_801 = buffer.data(smh + 801);
    const auto *smh_803 = buffer.data(smh + 803);
    const auto *smh_804 = buffer.data(smh + 804);
    const auto *smh_807 = buffer.data(smh + 807);
    const auto *smh_813 = buffer.data(smh + 813);
    const auto *smh_814 = buffer.data(smh + 814);
    const auto *smh_815 = buffer.data(smh + 815);
    const auto *smh_816 = buffer.data(smh + 816);
    const auto *smh_817 = buffer.data(smh + 817);
    const auto *smh_818 = buffer.data(smh + 818);
    const auto *smh_819 = buffer.data(smh + 819);
    const auto *smh_821 = buffer.data(smh + 821);
    const auto *smh_822 = buffer.data(smh + 822);
    const auto *smh_824 = buffer.data(smh + 824);
    const auto *smh_825 = buffer.data(smh + 825);
    const auto *smh_828 = buffer.data(smh + 828);
    const auto *smh_834 = buffer.data(smh + 834);
    const auto *smh_835 = buffer.data(smh + 835);
    const auto *smh_836 = buffer.data(smh + 836);
    const auto *smh_837 = buffer.data(smh + 837);
    const auto *smh_838 = buffer.data(smh + 838);
    const auto *smh_839 = buffer.data(smh + 839);
    const auto *smh_840 = buffer.data(smh + 840);
    const auto *smh_842 = buffer.data(smh + 842);
    const auto *smh_843 = buffer.data(smh + 843);
    const auto *smh_845 = buffer.data(smh + 845);
    const auto *smh_846 = buffer.data(smh + 846);
    const auto *smh_849 = buffer.data(smh + 849);
    const auto *smh_855 = buffer.data(smh + 855);
    const auto *smh_856 = buffer.data(smh + 856);
    const auto *smh_857 = buffer.data(smh + 857);
    const auto *smh_858 = buffer.data(smh + 858);
    const auto *smh_859 = buffer.data(smh + 859);
    const auto *smh_860 = buffer.data(smh + 860);

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pb_x, pc_x, pc_y, sli0_1022, slh_597, \
                         slh_770, slh_771, slh_772, sli1_1022, smh_765, smh_771, \
                         smh_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_15 * slh_597[k]
                    + f_3 * pc_y[k] * smh_765[k];

        t_1022[k] = pb_x[k] * sli0_1022[k]
                    + f_12 * slh_770[k]
                    - f_10 * pc_x[k] * sli1_1022[k];

        t_1023[k] = f_11 * slh_771[k]
                    + f_3 * pc_x[k] * smh_771[k];

        t_1024[k] = f_11 * slh_772[k]
                    + f_3 * pc_x[k] * smh_772[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pc_x, slh_773, slh_774, slh_775, \
                         slh_776, smh_773, smh_774, smh_775, smh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_11 * slh_773[k]
                    + f_3 * pc_x[k] * smh_773[k];

        t_1026[k] = f_11 * slh_774[k]
                    + f_3 * pc_x[k] * smh_774[k];

        t_1027[k] = f_11 * slh_775[k]
                    + f_3 * pc_x[k] * smh_775[k];

        t_1028[k] = f_11 * slh_776[k]
                    + f_3 * pc_x[k] * smh_776[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pb_x, pc_x, pc_z, sli0_1029, \
                         sli0_1031, sli0_1032, sli1_1029, sli1_1031, sli1_1032, \
                         smh_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = pb_x[k] * sli0_1029[k]
                    - f_10 * pc_x[k] * sli1_1029[k];

        t_1030[k] = f_3 * pc_z[k] * smh_771[k];

        t_1031[k] = pb_x[k] * sli0_1031[k]
                    - f_10 * pc_x[k] * sli1_1031[k];

        t_1032[k] = pb_x[k] * sli0_1032[k]
                    - f_10 * pc_x[k] * sli1_1032[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, pb_x, pc_x, pc_y, sli0_1033, sli0_1035, \
                         slh_608, sli1_1033, sli1_1035, smh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = pb_x[k] * sli0_1033[k]
                    - f_10 * pc_x[k] * sli1_1033[k];

        t_1034[k] = f_15 * slh_608[k]
                    + f_3 * pc_y[k] * smh_776[k];

        t_1035[k] = pb_x[k] * sli0_1035[k]
                    - f_10 * pc_x[k] * sli1_1035[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, t_1039, pb_z, pc_y, pc_z, sli0_784, sli0_787, \
                         slh_588, slh_609, sli1_784, sli1_787, \
                         smh_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = pb_z[k] * sli0_784[k]
                    - f_10 * pc_z[k] * sli1_784[k];

        t_1037[k] = f_16 * slh_609[k]
                    + f_3 * pc_y[k] * smh_777[k];

        t_1038[k] = f_11 * slh_588[k]
                    + f_3 * pc_z[k] * smh_777[k];

        t_1039[k] = pb_z[k] * sli0_787[k]
                    - f_10 * pc_z[k] * sli1_787[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pb_x, pb_z, pc_x, pc_y, pc_z, sli0_790, \
                         sli0_1041, slh_611, slh_782, sli1_790, sli1_1041, \
                         smh_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_16 * slh_611[k]
                    + f_3 * pc_y[k] * smh_779[k];

        t_1041[k] = pb_x[k] * sli0_1041[k]
                    + f_14 * slh_782[k]
                    - f_10 * pc_x[k] * sli1_1041[k];

        t_1042[k] = pb_z[k] * sli0_790[k]
                    - f_10 * pc_z[k] * sli1_790[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pb_x, pc_x, pc_y, pc_z, sli0_1045, slh_591, \
                         slh_614, slh_786, sli1_1045, smh_780, \
                         smh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_11 * slh_591[k]
                    + f_3 * pc_z[k] * smh_780[k];

        t_1044[k] = f_16 * slh_614[k]
                    + f_3 * pc_y[k] * smh_782[k];

        t_1045[k] = pb_x[k] * sli0_1045[k]
                    + f_13 * slh_786[k]
                    - f_10 * pc_x[k] * sli1_1045[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pb_x, pb_z, pc_x, pc_z, sli0_794, sli0_1048, \
                         slh_594, slh_789, sli1_794, sli1_1048, \
                         smh_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pb_z[k] * sli0_794[k]
                    - f_10 * pc_z[k] * sli1_794[k];

        t_1047[k] = f_11 * slh_594[k]
                    + f_3 * pc_z[k] * smh_783[k];

        t_1048[k] = pb_x[k] * sli0_1048[k]
                    + f_12 * slh_789[k]
                    - f_10 * pc_x[k] * sli1_1048[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, t_1052, pb_x, pc_x, pc_y, sli0_1050, slh_618, \
                         slh_791, slh_792, slh_793, sli1_1050, smh_786, smh_792, \
                         smh_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_16 * slh_618[k]
                    + f_3 * pc_y[k] * smh_786[k];

        t_1050[k] = pb_x[k] * sli0_1050[k]
                    + f_12 * slh_791[k]
                    - f_10 * pc_x[k] * sli1_1050[k];

        t_1051[k] = f_11 * slh_792[k]
                    + f_3 * pc_x[k] * smh_792[k];

        t_1052[k] = f_11 * slh_793[k]
                    + f_3 * pc_x[k] * smh_793[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, t_1056, pc_x, slh_794, slh_795, slh_796, \
                         slh_797, smh_794, smh_795, smh_796, smh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = f_11 * slh_794[k]
                    + f_3 * pc_x[k] * smh_794[k];

        t_1054[k] = f_11 * slh_795[k]
                    + f_3 * pc_x[k] * smh_795[k];

        t_1055[k] = f_11 * slh_796[k]
                    + f_3 * pc_x[k] * smh_796[k];

        t_1056[k] = f_11 * slh_797[k]
                    + f_3 * pc_x[k] * smh_797[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, t_1060, pb_x, pc_x, pc_z, sli0_1057, \
                         sli0_1059, sli0_1060, slh_603, sli1_1057, sli1_1059, sli1_1060, \
                         smh_792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = pb_x[k] * sli0_1057[k]
                    - f_10 * pc_x[k] * sli1_1057[k];

        t_1058[k] = f_11 * slh_603[k]
                    + f_3 * pc_z[k] * smh_792[k];

        t_1059[k] = pb_x[k] * sli0_1059[k]
                    - f_10 * pc_x[k] * sli1_1059[k];

        t_1060[k] = pb_x[k] * sli0_1060[k]
                    - f_10 * pc_x[k] * sli1_1060[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, pb_x, pc_x, pc_y, sli0_1061, \
                         sli0_1063, sli0_1064, slh_629, slh_798, sli1_1061, sli1_1063, \
                         sli1_1064, smh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = pb_x[k] * sli0_1061[k]
                    - f_10 * pc_x[k] * sli1_1061[k];

        t_1062[k] = f_16 * slh_629[k]
                    + f_3 * pc_y[k] * smh_797[k];

        t_1063[k] = pb_x[k] * sli0_1063[k]
                    - f_10 * pc_x[k] * sli1_1063[k];

        t_1064[k] = pb_x[k] * sli0_1064[k]
                    + f_17 * slh_798[k]
                    - f_10 * pc_x[k] * sli1_1064[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, pb_x, pc_x, pc_y, pc_z, sli0_1067, \
                         slh_609, slh_630, slh_632, slh_801, sli1_1067, smh_798, \
                         smh_800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_17 * slh_630[k]
                    + f_3 * pc_y[k] * smh_798[k];

        t_1066[k] = f_12 * slh_609[k]
                    + f_3 * pc_z[k] * smh_798[k];

        t_1067[k] = pb_x[k] * sli0_1067[k]
                    + f_14 * slh_801[k]
                    - f_10 * pc_x[k] * sli1_1067[k];

        t_1068[k] = f_17 * slh_632[k]
                    + f_3 * pc_y[k] * smh_800[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, pb_x, pc_x, pc_z, sli0_1069, sli0_1070, \
                         slh_612, slh_803, slh_804, sli1_1069, sli1_1070, \
                         smh_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = pb_x[k] * sli0_1069[k]
                    + f_14 * slh_803[k]
                    - f_10 * pc_x[k] * sli1_1069[k];

        t_1070[k] = pb_x[k] * sli0_1070[k]
                    + f_13 * slh_804[k]
                    - f_10 * pc_x[k] * sli1_1070[k];

        t_1071[k] = f_12 * slh_612[k]
                    + f_3 * pc_z[k] * smh_801[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, pb_x, pc_x, pc_y, sli0_1073, sli0_1074, \
                         slh_635, slh_807, slh_808, sli1_1073, sli1_1074, \
                         smh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = f_17 * slh_635[k]
                    + f_3 * pc_y[k] * smh_803[k];

        t_1073[k] = pb_x[k] * sli0_1073[k]
                    + f_13 * slh_807[k]
                    - f_10 * pc_x[k] * sli1_1073[k];

        t_1074[k] = pb_x[k] * sli0_1074[k]
                    + f_12 * slh_808[k]
                    - f_10 * pc_x[k] * sli1_1074[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, pb_x, pc_x, pc_y, pc_z, sli0_1076, slh_615, \
                         slh_639, slh_810, sli1_1076, smh_804, \
                         smh_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_12 * slh_615[k]
                    + f_3 * pc_z[k] * smh_804[k];

        t_1076[k] = pb_x[k] * sli0_1076[k]
                    + f_12 * slh_810[k]
                    - f_10 * pc_x[k] * sli1_1076[k];

        t_1077[k] = f_17 * slh_639[k]
                    + f_3 * pc_y[k] * smh_807[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, t_1081, pb_x, pc_x, sli0_1078, slh_812, \
                         slh_813, slh_814, slh_815, sli1_1078, smh_813, smh_814, \
                         smh_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = pb_x[k] * sli0_1078[k]
                    + f_12 * slh_812[k]
                    - f_10 * pc_x[k] * sli1_1078[k];

        t_1079[k] = f_11 * slh_813[k]
                    + f_3 * pc_x[k] * smh_813[k];

        t_1080[k] = f_11 * slh_814[k]
                    + f_3 * pc_x[k] * smh_814[k];

        t_1081[k] = f_11 * slh_815[k]
                    + f_3 * pc_x[k] * smh_815[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, t_1085, pb_x, pc_x, sli0_1085, slh_816, \
                         slh_817, slh_818, sli1_1085, smh_816, smh_817, \
                         smh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_11 * slh_816[k]
                    + f_3 * pc_x[k] * smh_816[k];

        t_1083[k] = f_11 * slh_817[k]
                    + f_3 * pc_x[k] * smh_817[k];

        t_1084[k] = f_11 * slh_818[k]
                    + f_3 * pc_x[k] * smh_818[k];

        t_1085[k] = pb_x[k] * sli0_1085[k]
                    - f_10 * pc_x[k] * sli1_1085[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, t_1089, pb_x, pc_x, pc_z, sli0_1087, \
                         sli0_1088, sli0_1089, slh_624, sli1_1087, sli1_1088, sli1_1089, \
                         smh_813 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_12 * slh_624[k]
                    + f_3 * pc_z[k] * smh_813[k];

        t_1087[k] = pb_x[k] * sli0_1087[k]
                    - f_10 * pc_x[k] * sli1_1087[k];

        t_1088[k] = pb_x[k] * sli0_1088[k]
                    - f_10 * pc_x[k] * sli1_1088[k];

        t_1089[k] = pb_x[k] * sli0_1089[k]
                    - f_10 * pc_x[k] * sli1_1089[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, t_1093, pb_x, pc_x, pc_y, sli0_1091, \
                         sli0_1092, slh_650, slh_651, slh_819, sli1_1091, sli1_1092, smh_818, \
                         smh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_17 * slh_650[k]
                    + f_3 * pc_y[k] * smh_818[k];

        t_1091[k] = pb_x[k] * sli0_1091[k]
                    - f_10 * pc_x[k] * sli1_1091[k];

        t_1092[k] = pb_x[k] * sli0_1092[k]
                    + f_17 * slh_819[k]
                    - f_10 * pc_x[k] * sli1_1092[k];

        t_1093[k] = f_18 * slh_651[k]
                    + f_3 * pc_y[k] * smh_819[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pb_x, pc_x, pc_y, pc_z, sli0_1095, slh_630, \
                         slh_653, slh_822, sli1_1095, smh_819, \
                         smh_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_13 * slh_630[k]
                    + f_3 * pc_z[k] * smh_819[k];

        t_1095[k] = pb_x[k] * sli0_1095[k]
                    + f_14 * slh_822[k]
                    - f_10 * pc_x[k] * sli1_1095[k];

        t_1096[k] = f_18 * slh_653[k]
                    + f_3 * pc_y[k] * smh_821[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pb_x, pc_x, pc_z, sli0_1097, sli0_1098, \
                         slh_633, slh_824, slh_825, sli1_1097, sli1_1098, \
                         smh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = pb_x[k] * sli0_1097[k]
                    + f_14 * slh_824[k]
                    - f_10 * pc_x[k] * sli1_1097[k];

        t_1098[k] = pb_x[k] * sli0_1098[k]
                    + f_13 * slh_825[k]
                    - f_10 * pc_x[k] * sli1_1098[k];

        t_1099[k] = f_13 * slh_633[k]
                    + f_3 * pc_z[k] * smh_822[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pb_x, pc_x, pc_y, sli0_1101, sli0_1102, \
                         slh_656, slh_828, slh_829, sli1_1101, sli1_1102, \
                         smh_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_18 * slh_656[k]
                    + f_3 * pc_y[k] * smh_824[k];

        t_1101[k] = pb_x[k] * sli0_1101[k]
                    + f_13 * slh_828[k]
                    - f_10 * pc_x[k] * sli1_1101[k];

        t_1102[k] = pb_x[k] * sli0_1102[k]
                    + f_12 * slh_829[k]
                    - f_10 * pc_x[k] * sli1_1102[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pb_x, pc_x, pc_y, pc_z, sli0_1104, slh_636, \
                         slh_660, slh_831, sli1_1104, smh_825, \
                         smh_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_13 * slh_636[k]
                    + f_3 * pc_z[k] * smh_825[k];

        t_1104[k] = pb_x[k] * sli0_1104[k]
                    + f_12 * slh_831[k]
                    - f_10 * pc_x[k] * sli1_1104[k];

        t_1105[k] = f_18 * slh_660[k]
                    + f_3 * pc_y[k] * smh_828[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pb_x, pc_x, sli0_1106, slh_833, \
                         slh_834, slh_835, slh_836, sli1_1106, smh_834, smh_835, \
                         smh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = pb_x[k] * sli0_1106[k]
                    + f_12 * slh_833[k]
                    - f_10 * pc_x[k] * sli1_1106[k];

        t_1107[k] = f_11 * slh_834[k]
                    + f_3 * pc_x[k] * smh_834[k];

        t_1108[k] = f_11 * slh_835[k]
                    + f_3 * pc_x[k] * smh_835[k];

        t_1109[k] = f_11 * slh_836[k]
                    + f_3 * pc_x[k] * smh_836[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pb_x, pc_x, sli0_1113, slh_837, \
                         slh_838, slh_839, sli1_1113, smh_837, smh_838, \
                         smh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_11 * slh_837[k]
                    + f_3 * pc_x[k] * smh_837[k];

        t_1111[k] = f_11 * slh_838[k]
                    + f_3 * pc_x[k] * smh_838[k];

        t_1112[k] = f_11 * slh_839[k]
                    + f_3 * pc_x[k] * smh_839[k];

        t_1113[k] = pb_x[k] * sli0_1113[k]
                    - f_10 * pc_x[k] * sli1_1113[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, t_1117, pb_x, pc_x, pc_z, sli0_1115, \
                         sli0_1116, sli0_1117, slh_645, sli1_1115, sli1_1116, sli1_1117, \
                         smh_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_13 * slh_645[k]
                    + f_3 * pc_z[k] * smh_834[k];

        t_1115[k] = pb_x[k] * sli0_1115[k]
                    - f_10 * pc_x[k] * sli1_1115[k];

        t_1116[k] = pb_x[k] * sli0_1116[k]
                    - f_10 * pc_x[k] * sli1_1116[k];

        t_1117[k] = pb_x[k] * sli0_1117[k]
                    - f_10 * pc_x[k] * sli1_1117[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, t_1121, pb_x, pc_x, pc_y, sli0_1119, \
                         sli0_1120, slh_671, slh_672, slh_840, sli1_1119, sli1_1120, smh_839, \
                         smh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_18 * slh_671[k]
                    + f_3 * pc_y[k] * smh_839[k];

        t_1119[k] = pb_x[k] * sli0_1119[k]
                    - f_10 * pc_x[k] * sli1_1119[k];

        t_1120[k] = pb_x[k] * sli0_1120[k]
                    + f_17 * slh_840[k]
                    - f_10 * pc_x[k] * sli1_1120[k];

        t_1121[k] = f_14 * slh_672[k]
                    + f_3 * pc_y[k] * smh_840[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, pb_x, pc_x, pc_y, pc_z, sli0_1123, slh_651, \
                         slh_674, slh_843, sli1_1123, smh_840, \
                         smh_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = f_14 * slh_651[k]
                    + f_3 * pc_z[k] * smh_840[k];

        t_1123[k] = pb_x[k] * sli0_1123[k]
                    + f_14 * slh_843[k]
                    - f_10 * pc_x[k] * sli1_1123[k];

        t_1124[k] = f_14 * slh_674[k]
                    + f_3 * pc_y[k] * smh_842[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, pb_x, pc_x, pc_z, sli0_1125, sli0_1126, \
                         slh_654, slh_845, slh_846, sli1_1125, sli1_1126, \
                         smh_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = pb_x[k] * sli0_1125[k]
                    + f_14 * slh_845[k]
                    - f_10 * pc_x[k] * sli1_1125[k];

        t_1126[k] = pb_x[k] * sli0_1126[k]
                    + f_13 * slh_846[k]
                    - f_10 * pc_x[k] * sli1_1126[k];

        t_1127[k] = f_14 * slh_654[k]
                    + f_3 * pc_z[k] * smh_843[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pb_x, pc_x, pc_y, sli0_1129, sli0_1130, \
                         slh_677, slh_849, slh_850, sli1_1129, sli1_1130, \
                         smh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_14 * slh_677[k]
                    + f_3 * pc_y[k] * smh_845[k];

        t_1129[k] = pb_x[k] * sli0_1129[k]
                    + f_13 * slh_849[k]
                    - f_10 * pc_x[k] * sli1_1129[k];

        t_1130[k] = pb_x[k] * sli0_1130[k]
                    + f_12 * slh_850[k]
                    - f_10 * pc_x[k] * sli1_1130[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pb_x, pc_x, pc_y, pc_z, sli0_1132, slh_657, \
                         slh_681, slh_852, sli1_1132, smh_846, \
                         smh_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_14 * slh_657[k]
                    + f_3 * pc_z[k] * smh_846[k];

        t_1132[k] = pb_x[k] * sli0_1132[k]
                    + f_12 * slh_852[k]
                    - f_10 * pc_x[k] * sli1_1132[k];

        t_1133[k] = f_14 * slh_681[k]
                    + f_3 * pc_y[k] * smh_849[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, t_1137, pb_x, pc_x, sli0_1134, slh_854, \
                         slh_855, slh_856, slh_857, sli1_1134, smh_855, smh_856, \
                         smh_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = pb_x[k] * sli0_1134[k]
                    + f_12 * slh_854[k]
                    - f_10 * pc_x[k] * sli1_1134[k];

        t_1135[k] = f_11 * slh_855[k]
                    + f_3 * pc_x[k] * smh_855[k];

        t_1136[k] = f_11 * slh_856[k]
                    + f_3 * pc_x[k] * smh_856[k];

        t_1137[k] = f_11 * slh_857[k]
                    + f_3 * pc_x[k] * smh_857[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, t_1141, pb_x, pc_x, sli0_1141, slh_858, \
                         slh_859, slh_860, sli1_1141, smh_858, smh_859, \
                         smh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_11 * slh_858[k]
                    + f_3 * pc_x[k] * smh_858[k];

        t_1139[k] = f_11 * slh_859[k]
                    + f_3 * pc_x[k] * smh_859[k];

        t_1140[k] = f_11 * slh_860[k]
                    + f_3 * pc_x[k] * smh_860[k];

        t_1141[k] = pb_x[k] * sli0_1141[k]
                    - f_10 * pc_x[k] * sli1_1141[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sli0,
                                                           const size_t slh, const size_t sli1,
                                                           const size_t smg0, const size_t smg1,
                                                           const size_t smh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);
    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);
    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_980 = buffer.data(sli0 + 980);
    const auto *sli0_985 = buffer.data(sli0 + 985);
    const auto *sli0_989 = buffer.data(sli0 + 989);
    const auto *sli0_994 = buffer.data(sli0 + 994);
    const auto *sli0_1143 = buffer.data(sli0 + 1143);
    const auto *sli0_1144 = buffer.data(sli0 + 1144);
    const auto *sli0_1145 = buffer.data(sli0 + 1145);
    const auto *sli0_1147 = buffer.data(sli0 + 1147);
    const auto *sli0_1148 = buffer.data(sli0 + 1148);
    const auto *sli0_1151 = buffer.data(sli0 + 1151);
    const auto *sli0_1153 = buffer.data(sli0 + 1153);
    const auto *sli0_1154 = buffer.data(sli0 + 1154);
    const auto *sli0_1157 = buffer.data(sli0 + 1157);
    const auto *sli0_1158 = buffer.data(sli0 + 1158);
    const auto *sli0_1160 = buffer.data(sli0 + 1160);
    const auto *sli0_1162 = buffer.data(sli0 + 1162);
    const auto *sli0_1169 = buffer.data(sli0 + 1169);
    const auto *sli0_1171 = buffer.data(sli0 + 1171);
    const auto *sli0_1172 = buffer.data(sli0 + 1172);
    const auto *sli0_1173 = buffer.data(sli0 + 1173);
    const auto *sli0_1175 = buffer.data(sli0 + 1175);
    const auto *sli0_1176 = buffer.data(sli0 + 1176);
    const auto *sli0_1179 = buffer.data(sli0 + 1179);
    const auto *sli0_1181 = buffer.data(sli0 + 1181);
    const auto *sli0_1182 = buffer.data(sli0 + 1182);
    const auto *sli0_1185 = buffer.data(sli0 + 1185);
    const auto *sli0_1186 = buffer.data(sli0 + 1186);
    const auto *sli0_1188 = buffer.data(sli0 + 1188);
    const auto *sli0_1190 = buffer.data(sli0 + 1190);
    const auto *sli0_1197 = buffer.data(sli0 + 1197);
    const auto *sli0_1199 = buffer.data(sli0 + 1199);
    const auto *sli0_1200 = buffer.data(sli0 + 1200);
    const auto *sli0_1201 = buffer.data(sli0 + 1201);
    const auto *sli0_1203 = buffer.data(sli0 + 1203);
    const auto *sli0_1207 = buffer.data(sli0 + 1207);
    const auto *sli0_1210 = buffer.data(sli0 + 1210);
    const auto *sli0_1214 = buffer.data(sli0 + 1214);
    const auto *sli0_1216 = buffer.data(sli0 + 1216);
    const auto *sli0_1225 = buffer.data(sli0 + 1225);
    const auto *sli0_1227 = buffer.data(sli0 + 1227);
    const auto *sli0_1228 = buffer.data(sli0 + 1228);
    const auto *sli0_1229 = buffer.data(sli0 + 1229);
    const auto *sli0_1231 = buffer.data(sli0 + 1231);
    const auto *sli0_1232 = buffer.data(sli0 + 1232);
    const auto *sli0_1235 = buffer.data(sli0 + 1235);
    const auto *sli0_1237 = buffer.data(sli0 + 1237);
    const auto *sli0_1238 = buffer.data(sli0 + 1238);
    const auto *sli0_1241 = buffer.data(sli0 + 1241);
    const auto *sli0_1242 = buffer.data(sli0 + 1242);
    const auto *sli0_1244 = buffer.data(sli0 + 1244);
    const auto *sli0_1246 = buffer.data(sli0 + 1246);
    const auto *sli0_1253 = buffer.data(sli0 + 1253);
    const auto *sli0_1255 = buffer.data(sli0 + 1255);
    const auto *sli0_1256 = buffer.data(sli0 + 1256);
    const auto *sli0_1257 = buffer.data(sli0 + 1257);
    const auto *sli0_1259 = buffer.data(sli0 + 1259);

    const auto *slh_666 = buffer.data(slh + 666);
    const auto *slh_672 = buffer.data(slh + 672);
    const auto *slh_675 = buffer.data(slh + 675);
    const auto *slh_678 = buffer.data(slh + 678);
    const auto *slh_687 = buffer.data(slh + 687);
    const auto *slh_692 = buffer.data(slh + 692);
    const auto *slh_693 = buffer.data(slh + 693);
    const auto *slh_695 = buffer.data(slh + 695);
    const auto *slh_696 = buffer.data(slh + 696);
    const auto *slh_698 = buffer.data(slh + 698);
    const auto *slh_699 = buffer.data(slh + 699);
    const auto *slh_702 = buffer.data(slh + 702);
    const auto *slh_708 = buffer.data(slh + 708);
    const auto *slh_713 = buffer.data(slh + 713);
    const auto *slh_714 = buffer.data(slh + 714);
    const auto *slh_716 = buffer.data(slh + 716);
    const auto *slh_717 = buffer.data(slh + 717);
    const auto *slh_719 = buffer.data(slh + 719);
    const auto *slh_720 = buffer.data(slh + 720);
    const auto *slh_723 = buffer.data(slh + 723);
    const auto *slh_729 = buffer.data(slh + 729);
    const auto *slh_734 = buffer.data(slh + 734);
    const auto *slh_735 = buffer.data(slh + 735);
    const auto *slh_737 = buffer.data(slh + 737);
    const auto *slh_738 = buffer.data(slh + 738);
    const auto *slh_740 = buffer.data(slh + 740);
    const auto *slh_741 = buffer.data(slh + 741);
    const auto *slh_744 = buffer.data(slh + 744);
    const auto *slh_750 = buffer.data(slh + 750);
    const auto *slh_755 = buffer.data(slh + 755);
    const auto *slh_756 = buffer.data(slh + 756);
    const auto *slh_861 = buffer.data(slh + 861);
    const auto *slh_864 = buffer.data(slh + 864);
    const auto *slh_866 = buffer.data(slh + 866);
    const auto *slh_867 = buffer.data(slh + 867);
    const auto *slh_870 = buffer.data(slh + 870);
    const auto *slh_871 = buffer.data(slh + 871);
    const auto *slh_873 = buffer.data(slh + 873);
    const auto *slh_875 = buffer.data(slh + 875);
    const auto *slh_876 = buffer.data(slh + 876);
    const auto *slh_877 = buffer.data(slh + 877);
    const auto *slh_878 = buffer.data(slh + 878);
    const auto *slh_879 = buffer.data(slh + 879);
    const auto *slh_880 = buffer.data(slh + 880);
    const auto *slh_881 = buffer.data(slh + 881);
    const auto *slh_882 = buffer.data(slh + 882);
    const auto *slh_885 = buffer.data(slh + 885);
    const auto *slh_887 = buffer.data(slh + 887);
    const auto *slh_888 = buffer.data(slh + 888);
    const auto *slh_891 = buffer.data(slh + 891);
    const auto *slh_892 = buffer.data(slh + 892);
    const auto *slh_894 = buffer.data(slh + 894);
    const auto *slh_896 = buffer.data(slh + 896);
    const auto *slh_897 = buffer.data(slh + 897);
    const auto *slh_898 = buffer.data(slh + 898);
    const auto *slh_899 = buffer.data(slh + 899);
    const auto *slh_900 = buffer.data(slh + 900);
    const auto *slh_901 = buffer.data(slh + 901);
    const auto *slh_902 = buffer.data(slh + 902);
    const auto *slh_906 = buffer.data(slh + 906);
    const auto *slh_909 = buffer.data(slh + 909);
    const auto *slh_913 = buffer.data(slh + 913);
    const auto *slh_915 = buffer.data(slh + 915);
    const auto *slh_918 = buffer.data(slh + 918);
    const auto *slh_919 = buffer.data(slh + 919);
    const auto *slh_920 = buffer.data(slh + 920);
    const auto *slh_921 = buffer.data(slh + 921);
    const auto *slh_922 = buffer.data(slh + 922);
    const auto *slh_923 = buffer.data(slh + 923);
    const auto *slh_924 = buffer.data(slh + 924);
    const auto *slh_927 = buffer.data(slh + 927);
    const auto *slh_929 = buffer.data(slh + 929);
    const auto *slh_930 = buffer.data(slh + 930);
    const auto *slh_933 = buffer.data(slh + 933);
    const auto *slh_934 = buffer.data(slh + 934);
    const auto *slh_936 = buffer.data(slh + 936);
    const auto *slh_938 = buffer.data(slh + 938);
    const auto *slh_939 = buffer.data(slh + 939);
    const auto *slh_940 = buffer.data(slh + 940);
    const auto *slh_941 = buffer.data(slh + 941);
    const auto *slh_942 = buffer.data(slh + 942);
    const auto *slh_943 = buffer.data(slh + 943);
    const auto *slh_944 = buffer.data(slh + 944);

    const auto *sli1_980 = buffer.data(sli1 + 980);
    const auto *sli1_985 = buffer.data(sli1 + 985);
    const auto *sli1_989 = buffer.data(sli1 + 989);
    const auto *sli1_994 = buffer.data(sli1 + 994);
    const auto *sli1_1143 = buffer.data(sli1 + 1143);
    const auto *sli1_1144 = buffer.data(sli1 + 1144);
    const auto *sli1_1145 = buffer.data(sli1 + 1145);
    const auto *sli1_1147 = buffer.data(sli1 + 1147);
    const auto *sli1_1148 = buffer.data(sli1 + 1148);
    const auto *sli1_1151 = buffer.data(sli1 + 1151);
    const auto *sli1_1153 = buffer.data(sli1 + 1153);
    const auto *sli1_1154 = buffer.data(sli1 + 1154);
    const auto *sli1_1157 = buffer.data(sli1 + 1157);
    const auto *sli1_1158 = buffer.data(sli1 + 1158);
    const auto *sli1_1160 = buffer.data(sli1 + 1160);
    const auto *sli1_1162 = buffer.data(sli1 + 1162);
    const auto *sli1_1169 = buffer.data(sli1 + 1169);
    const auto *sli1_1171 = buffer.data(sli1 + 1171);
    const auto *sli1_1172 = buffer.data(sli1 + 1172);
    const auto *sli1_1173 = buffer.data(sli1 + 1173);
    const auto *sli1_1175 = buffer.data(sli1 + 1175);
    const auto *sli1_1176 = buffer.data(sli1 + 1176);
    const auto *sli1_1179 = buffer.data(sli1 + 1179);
    const auto *sli1_1181 = buffer.data(sli1 + 1181);
    const auto *sli1_1182 = buffer.data(sli1 + 1182);
    const auto *sli1_1185 = buffer.data(sli1 + 1185);
    const auto *sli1_1186 = buffer.data(sli1 + 1186);
    const auto *sli1_1188 = buffer.data(sli1 + 1188);
    const auto *sli1_1190 = buffer.data(sli1 + 1190);
    const auto *sli1_1197 = buffer.data(sli1 + 1197);
    const auto *sli1_1199 = buffer.data(sli1 + 1199);
    const auto *sli1_1200 = buffer.data(sli1 + 1200);
    const auto *sli1_1201 = buffer.data(sli1 + 1201);
    const auto *sli1_1203 = buffer.data(sli1 + 1203);
    const auto *sli1_1207 = buffer.data(sli1 + 1207);
    const auto *sli1_1210 = buffer.data(sli1 + 1210);
    const auto *sli1_1214 = buffer.data(sli1 + 1214);
    const auto *sli1_1216 = buffer.data(sli1 + 1216);
    const auto *sli1_1225 = buffer.data(sli1 + 1225);
    const auto *sli1_1227 = buffer.data(sli1 + 1227);
    const auto *sli1_1228 = buffer.data(sli1 + 1228);
    const auto *sli1_1229 = buffer.data(sli1 + 1229);
    const auto *sli1_1231 = buffer.data(sli1 + 1231);
    const auto *sli1_1232 = buffer.data(sli1 + 1232);
    const auto *sli1_1235 = buffer.data(sli1 + 1235);
    const auto *sli1_1237 = buffer.data(sli1 + 1237);
    const auto *sli1_1238 = buffer.data(sli1 + 1238);
    const auto *sli1_1241 = buffer.data(sli1 + 1241);
    const auto *sli1_1242 = buffer.data(sli1 + 1242);
    const auto *sli1_1244 = buffer.data(sli1 + 1244);
    const auto *sli1_1246 = buffer.data(sli1 + 1246);
    const auto *sli1_1253 = buffer.data(sli1 + 1253);
    const auto *sli1_1255 = buffer.data(sli1 + 1255);
    const auto *sli1_1256 = buffer.data(sli1 + 1256);
    const auto *sli1_1257 = buffer.data(sli1 + 1257);
    const auto *sli1_1259 = buffer.data(sli1 + 1259);

    const auto *smg0_675 = buffer.data(smg0 + 675);

    const auto *smg1_675 = buffer.data(smg1 + 675);

    const auto *smh_855 = buffer.data(smh + 855);
    const auto *smh_860 = buffer.data(smh + 860);
    const auto *smh_861 = buffer.data(smh + 861);
    const auto *smh_863 = buffer.data(smh + 863);
    const auto *smh_864 = buffer.data(smh + 864);
    const auto *smh_866 = buffer.data(smh + 866);
    const auto *smh_867 = buffer.data(smh + 867);
    const auto *smh_870 = buffer.data(smh + 870);
    const auto *smh_876 = buffer.data(smh + 876);
    const auto *smh_877 = buffer.data(smh + 877);
    const auto *smh_878 = buffer.data(smh + 878);
    const auto *smh_879 = buffer.data(smh + 879);
    const auto *smh_880 = buffer.data(smh + 880);
    const auto *smh_881 = buffer.data(smh + 881);
    const auto *smh_882 = buffer.data(smh + 882);
    const auto *smh_884 = buffer.data(smh + 884);
    const auto *smh_885 = buffer.data(smh + 885);
    const auto *smh_887 = buffer.data(smh + 887);
    const auto *smh_888 = buffer.data(smh + 888);
    const auto *smh_891 = buffer.data(smh + 891);
    const auto *smh_897 = buffer.data(smh + 897);
    const auto *smh_898 = buffer.data(smh + 898);
    const auto *smh_899 = buffer.data(smh + 899);
    const auto *smh_900 = buffer.data(smh + 900);
    const auto *smh_901 = buffer.data(smh + 901);
    const auto *smh_902 = buffer.data(smh + 902);
    const auto *smh_903 = buffer.data(smh + 903);
    const auto *smh_905 = buffer.data(smh + 905);
    const auto *smh_906 = buffer.data(smh + 906);
    const auto *smh_908 = buffer.data(smh + 908);
    const auto *smh_909 = buffer.data(smh + 909);
    const auto *smh_912 = buffer.data(smh + 912);
    const auto *smh_918 = buffer.data(smh + 918);
    const auto *smh_919 = buffer.data(smh + 919);
    const auto *smh_920 = buffer.data(smh + 920);
    const auto *smh_921 = buffer.data(smh + 921);
    const auto *smh_922 = buffer.data(smh + 922);
    const auto *smh_923 = buffer.data(smh + 923);
    const auto *smh_924 = buffer.data(smh + 924);
    const auto *smh_926 = buffer.data(smh + 926);
    const auto *smh_927 = buffer.data(smh + 927);
    const auto *smh_929 = buffer.data(smh + 929);
    const auto *smh_930 = buffer.data(smh + 930);
    const auto *smh_933 = buffer.data(smh + 933);
    const auto *smh_939 = buffer.data(smh + 939);
    const auto *smh_940 = buffer.data(smh + 940);
    const auto *smh_941 = buffer.data(smh + 941);
    const auto *smh_942 = buffer.data(smh + 942);
    const auto *smh_943 = buffer.data(smh + 943);
    const auto *smh_944 = buffer.data(smh + 944);
    const auto *smh_945 = buffer.data(smh + 945);

#pragma omp simd aligned(t_1142, t_1143, t_1144, t_1145, pb_x, pc_x, pc_z, sli0_1143, \
                         sli0_1144, sli0_1145, slh_666, sli1_1143, sli1_1144, sli1_1145, \
                         smh_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_14 * slh_666[k]
                    + f_3 * pc_z[k] * smh_855[k];

        t_1143[k] = pb_x[k] * sli0_1143[k]
                    - f_10 * pc_x[k] * sli1_1143[k];

        t_1144[k] = pb_x[k] * sli0_1144[k]
                    - f_10 * pc_x[k] * sli1_1144[k];

        t_1145[k] = pb_x[k] * sli0_1145[k]
                    - f_10 * pc_x[k] * sli1_1145[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, t_1149, pb_x, pc_x, pc_y, sli0_1147, \
                         sli0_1148, slh_692, slh_693, slh_861, sli1_1147, sli1_1148, smh_860, \
                         smh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_14 * slh_692[k]
                    + f_3 * pc_y[k] * smh_860[k];

        t_1147[k] = pb_x[k] * sli0_1147[k]
                    - f_10 * pc_x[k] * sli1_1147[k];

        t_1148[k] = pb_x[k] * sli0_1148[k]
                    + f_17 * slh_861[k]
                    - f_10 * pc_x[k] * sli1_1148[k];

        t_1149[k] = f_13 * slh_693[k]
                    + f_3 * pc_y[k] * smh_861[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, pb_x, pc_x, pc_y, pc_z, sli0_1151, slh_672, \
                         slh_695, slh_864, sli1_1151, smh_861, \
                         smh_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_18 * slh_672[k]
                    + f_3 * pc_z[k] * smh_861[k];

        t_1151[k] = pb_x[k] * sli0_1151[k]
                    + f_14 * slh_864[k]
                    - f_10 * pc_x[k] * sli1_1151[k];

        t_1152[k] = f_13 * slh_695[k]
                    + f_3 * pc_y[k] * smh_863[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, pb_x, pc_x, pc_z, sli0_1153, sli0_1154, \
                         slh_675, slh_866, slh_867, sli1_1153, sli1_1154, \
                         smh_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = pb_x[k] * sli0_1153[k]
                    + f_14 * slh_866[k]
                    - f_10 * pc_x[k] * sli1_1153[k];

        t_1154[k] = pb_x[k] * sli0_1154[k]
                    + f_13 * slh_867[k]
                    - f_10 * pc_x[k] * sli1_1154[k];

        t_1155[k] = f_18 * slh_675[k]
                    + f_3 * pc_z[k] * smh_864[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, pb_x, pc_x, pc_y, sli0_1157, sli0_1158, \
                         slh_698, slh_870, slh_871, sli1_1157, sli1_1158, \
                         smh_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_13 * slh_698[k]
                    + f_3 * pc_y[k] * smh_866[k];

        t_1157[k] = pb_x[k] * sli0_1157[k]
                    + f_13 * slh_870[k]
                    - f_10 * pc_x[k] * sli1_1157[k];

        t_1158[k] = pb_x[k] * sli0_1158[k]
                    + f_12 * slh_871[k]
                    - f_10 * pc_x[k] * sli1_1158[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, pb_x, pc_x, pc_y, pc_z, sli0_1160, slh_678, \
                         slh_702, slh_873, sli1_1160, smh_867, \
                         smh_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_18 * slh_678[k]
                    + f_3 * pc_z[k] * smh_867[k];

        t_1160[k] = pb_x[k] * sli0_1160[k]
                    + f_12 * slh_873[k]
                    - f_10 * pc_x[k] * sli1_1160[k];

        t_1161[k] = f_13 * slh_702[k]
                    + f_3 * pc_y[k] * smh_870[k];
    }

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, pb_x, pc_x, sli0_1162, slh_875, \
                         slh_876, slh_877, slh_878, sli1_1162, smh_876, smh_877, \
                         smh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = pb_x[k] * sli0_1162[k]
                    + f_12 * slh_875[k]
                    - f_10 * pc_x[k] * sli1_1162[k];

        t_1163[k] = f_11 * slh_876[k]
                    + f_3 * pc_x[k] * smh_876[k];

        t_1164[k] = f_11 * slh_877[k]
                    + f_3 * pc_x[k] * smh_877[k];

        t_1165[k] = f_11 * slh_878[k]
                    + f_3 * pc_x[k] * smh_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, pb_x, pc_x, sli0_1169, slh_879, \
                         slh_880, slh_881, sli1_1169, smh_879, smh_880, \
                         smh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_11 * slh_879[k]
                    + f_3 * pc_x[k] * smh_879[k];

        t_1167[k] = f_11 * slh_880[k]
                    + f_3 * pc_x[k] * smh_880[k];

        t_1168[k] = f_11 * slh_881[k]
                    + f_3 * pc_x[k] * smh_881[k];

        t_1169[k] = pb_x[k] * sli0_1169[k]
                    - f_10 * pc_x[k] * sli1_1169[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, t_1173, pb_x, pc_x, pc_z, sli0_1171, \
                         sli0_1172, sli0_1173, slh_687, sli1_1171, sli1_1172, sli1_1173, \
                         smh_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_18 * slh_687[k]
                    + f_3 * pc_z[k] * smh_876[k];

        t_1171[k] = pb_x[k] * sli0_1171[k]
                    - f_10 * pc_x[k] * sli1_1171[k];

        t_1172[k] = pb_x[k] * sli0_1172[k]
                    - f_10 * pc_x[k] * sli1_1172[k];

        t_1173[k] = pb_x[k] * sli0_1173[k]
                    - f_10 * pc_x[k] * sli1_1173[k];
    }

#pragma omp simd aligned(t_1174, t_1175, t_1176, t_1177, pb_x, pc_x, pc_y, sli0_1175, \
                         sli0_1176, slh_713, slh_714, slh_882, sli1_1175, sli1_1176, smh_881, \
                         smh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1174[k] = f_13 * slh_713[k]
                    + f_3 * pc_y[k] * smh_881[k];

        t_1175[k] = pb_x[k] * sli0_1175[k]
                    - f_10 * pc_x[k] * sli1_1175[k];

        t_1176[k] = pb_x[k] * sli0_1176[k]
                    + f_17 * slh_882[k]
                    - f_10 * pc_x[k] * sli1_1176[k];

        t_1177[k] = f_12 * slh_714[k]
                    + f_3 * pc_y[k] * smh_882[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, pb_x, pc_x, pc_y, pc_z, sli0_1179, slh_693, \
                         slh_716, slh_885, sli1_1179, smh_882, \
                         smh_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_17 * slh_693[k]
                    + f_3 * pc_z[k] * smh_882[k];

        t_1179[k] = pb_x[k] * sli0_1179[k]
                    + f_14 * slh_885[k]
                    - f_10 * pc_x[k] * sli1_1179[k];

        t_1180[k] = f_12 * slh_716[k]
                    + f_3 * pc_y[k] * smh_884[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, pb_x, pc_x, pc_z, sli0_1181, sli0_1182, \
                         slh_696, slh_887, slh_888, sli1_1181, sli1_1182, \
                         smh_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = pb_x[k] * sli0_1181[k]
                    + f_14 * slh_887[k]
                    - f_10 * pc_x[k] * sli1_1181[k];

        t_1182[k] = pb_x[k] * sli0_1182[k]
                    + f_13 * slh_888[k]
                    - f_10 * pc_x[k] * sli1_1182[k];

        t_1183[k] = f_17 * slh_696[k]
                    + f_3 * pc_z[k] * smh_885[k];
    }

#pragma omp simd aligned(t_1184, t_1185, t_1186, pb_x, pc_x, pc_y, sli0_1185, sli0_1186, \
                         slh_719, slh_891, slh_892, sli1_1185, sli1_1186, \
                         smh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1184[k] = f_12 * slh_719[k]
                    + f_3 * pc_y[k] * smh_887[k];

        t_1185[k] = pb_x[k] * sli0_1185[k]
                    + f_13 * slh_891[k]
                    - f_10 * pc_x[k] * sli1_1185[k];

        t_1186[k] = pb_x[k] * sli0_1186[k]
                    + f_12 * slh_892[k]
                    - f_10 * pc_x[k] * sli1_1186[k];
    }

#pragma omp simd aligned(t_1187, t_1188, t_1189, pb_x, pc_x, pc_y, pc_z, sli0_1188, slh_699, \
                         slh_723, slh_894, sli1_1188, smh_888, \
                         smh_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1187[k] = f_17 * slh_699[k]
                    + f_3 * pc_z[k] * smh_888[k];

        t_1188[k] = pb_x[k] * sli0_1188[k]
                    + f_12 * slh_894[k]
                    - f_10 * pc_x[k] * sli1_1188[k];

        t_1189[k] = f_12 * slh_723[k]
                    + f_3 * pc_y[k] * smh_891[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, t_1193, pb_x, pc_x, sli0_1190, slh_896, \
                         slh_897, slh_898, slh_899, sli1_1190, smh_897, smh_898, \
                         smh_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = pb_x[k] * sli0_1190[k]
                    + f_12 * slh_896[k]
                    - f_10 * pc_x[k] * sli1_1190[k];

        t_1191[k] = f_11 * slh_897[k]
                    + f_3 * pc_x[k] * smh_897[k];

        t_1192[k] = f_11 * slh_898[k]
                    + f_3 * pc_x[k] * smh_898[k];

        t_1193[k] = f_11 * slh_899[k]
                    + f_3 * pc_x[k] * smh_899[k];
    }

#pragma omp simd aligned(t_1194, t_1195, t_1196, t_1197, pb_x, pc_x, sli0_1197, slh_900, \
                         slh_901, slh_902, sli1_1197, smh_900, smh_901, \
                         smh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1194[k] = f_11 * slh_900[k]
                    + f_3 * pc_x[k] * smh_900[k];

        t_1195[k] = f_11 * slh_901[k]
                    + f_3 * pc_x[k] * smh_901[k];

        t_1196[k] = f_11 * slh_902[k]
                    + f_3 * pc_x[k] * smh_902[k];

        t_1197[k] = pb_x[k] * sli0_1197[k]
                    - f_10 * pc_x[k] * sli1_1197[k];
    }

#pragma omp simd aligned(t_1198, t_1199, t_1200, t_1201, pb_x, pc_x, pc_z, sli0_1199, \
                         sli0_1200, sli0_1201, slh_708, sli1_1199, sli1_1200, sli1_1201, \
                         smh_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1198[k] = f_17 * slh_708[k]
                    + f_3 * pc_z[k] * smh_897[k];

        t_1199[k] = pb_x[k] * sli0_1199[k]
                    - f_10 * pc_x[k] * sli1_1199[k];

        t_1200[k] = pb_x[k] * sli0_1200[k]
                    - f_10 * pc_x[k] * sli1_1200[k];

        t_1201[k] = pb_x[k] * sli0_1201[k]
                    - f_10 * pc_x[k] * sli1_1201[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, pb_x, pb_y, pc_x, pc_y, sli0_980, \
                         sli0_1203, slh_734, slh_735, sli1_980, sli1_1203, smh_902, \
                         smh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_12 * slh_734[k]
                    + f_3 * pc_y[k] * smh_902[k];

        t_1203[k] = pb_x[k] * sli0_1203[k]
                    - f_10 * pc_x[k] * sli1_1203[k];

        t_1204[k] = pb_y[k] * sli0_980[k]
                    - f_10 * pc_y[k] * sli1_980[k];

        t_1205[k] = f_11 * slh_735[k]
                    + f_3 * pc_y[k] * smh_903[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, pb_x, pc_x, pc_y, pc_z, sli0_1207, slh_714, \
                         slh_737, slh_906, sli1_1207, smh_903, \
                         smh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_16 * slh_714[k]
                    + f_3 * pc_z[k] * smh_903[k];

        t_1207[k] = pb_x[k] * sli0_1207[k]
                    + f_14 * slh_906[k]
                    - f_10 * pc_x[k] * sli1_1207[k];

        t_1208[k] = f_11 * slh_737[k]
                    + f_3 * pc_y[k] * smh_905[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, pb_x, pb_y, pc_x, pc_y, pc_z, sli0_985, \
                         sli0_1210, slh_717, slh_909, sli1_985, sli1_1210, \
                         smh_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = pb_y[k] * sli0_985[k]
                    - f_10 * pc_y[k] * sli1_985[k];

        t_1210[k] = pb_x[k] * sli0_1210[k]
                    + f_13 * slh_909[k]
                    - f_10 * pc_x[k] * sli1_1210[k];

        t_1211[k] = f_16 * slh_717[k]
                    + f_3 * pc_z[k] * smh_906[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, pb_x, pb_y, pc_x, pc_y, sli0_989, sli0_1214, \
                         slh_740, slh_913, sli1_989, sli1_1214, \
                         smh_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_11 * slh_740[k]
                    + f_3 * pc_y[k] * smh_908[k];

        t_1213[k] = pb_y[k] * sli0_989[k]
                    - f_10 * pc_y[k] * sli1_989[k];

        t_1214[k] = pb_x[k] * sli0_1214[k]
                    + f_12 * slh_913[k]
                    - f_10 * pc_x[k] * sli1_1214[k];
    }

#pragma omp simd aligned(t_1215, t_1216, t_1217, pb_x, pc_x, pc_y, pc_z, sli0_1216, slh_720, \
                         slh_744, slh_915, sli1_1216, smh_909, \
                         smh_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1215[k] = f_16 * slh_720[k]
                    + f_3 * pc_z[k] * smh_909[k];

        t_1216[k] = pb_x[k] * sli0_1216[k]
                    + f_12 * slh_915[k]
                    - f_10 * pc_x[k] * sli1_1216[k];

        t_1217[k] = f_11 * slh_744[k]
                    + f_3 * pc_y[k] * smh_912[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, t_1221, pb_y, pc_x, pc_y, sli0_994, slh_918, \
                         slh_919, slh_920, sli1_994, smh_918, smh_919, \
                         smh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = pb_y[k] * sli0_994[k]
                    - f_10 * pc_y[k] * sli1_994[k];

        t_1219[k] = f_11 * slh_918[k]
                    + f_3 * pc_x[k] * smh_918[k];

        t_1220[k] = f_11 * slh_919[k]
                    + f_3 * pc_x[k] * smh_919[k];

        t_1221[k] = f_11 * slh_920[k]
                    + f_3 * pc_x[k] * smh_920[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pb_x, pc_x, sli0_1225, slh_921, \
                         slh_922, slh_923, sli1_1225, smh_921, smh_922, \
                         smh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_11 * slh_921[k]
                    + f_3 * pc_x[k] * smh_921[k];

        t_1223[k] = f_11 * slh_922[k]
                    + f_3 * pc_x[k] * smh_922[k];

        t_1224[k] = f_11 * slh_923[k]
                    + f_3 * pc_x[k] * smh_923[k];

        t_1225[k] = pb_x[k] * sli0_1225[k]
                    - f_10 * pc_x[k] * sli1_1225[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, pb_x, pc_x, pc_z, sli0_1227, \
                         sli0_1228, sli0_1229, slh_729, sli1_1227, sli1_1228, sli1_1229, \
                         smh_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_16 * slh_729[k]
                    + f_3 * pc_z[k] * smh_918[k];

        t_1227[k] = pb_x[k] * sli0_1227[k]
                    - f_10 * pc_x[k] * sli1_1227[k];

        t_1228[k] = pb_x[k] * sli0_1228[k]
                    - f_10 * pc_x[k] * sli1_1228[k];

        t_1229[k] = pb_x[k] * sli0_1229[k]
                    - f_10 * pc_x[k] * sli1_1229[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pb_x, pc_x, pc_y, sli0_1231, \
                         sli0_1232, slh_755, slh_924, sli1_1231, sli1_1232, smh_923, \
                         smh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_11 * slh_755[k]
                    + f_3 * pc_y[k] * smh_923[k];

        t_1231[k] = pb_x[k] * sli0_1231[k]
                    - f_10 * pc_x[k] * sli1_1231[k];

        t_1232[k] = pb_x[k] * sli0_1232[k]
                    + f_17 * slh_924[k]
                    - f_10 * pc_x[k] * sli1_1232[k];

        t_1233[k] = f_3 * pc_y[k] * smh_924[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pb_x, pc_x, pc_y, pc_z, sli0_1235, slh_735, \
                         slh_927, sli1_1235, smh_924, smh_926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = f_15 * slh_735[k]
                    + f_3 * pc_z[k] * smh_924[k];

        t_1235[k] = pb_x[k] * sli0_1235[k]
                    + f_14 * slh_927[k]
                    - f_10 * pc_x[k] * sli1_1235[k];

        t_1236[k] = f_3 * pc_y[k] * smh_926[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, pb_x, pc_x, pc_z, sli0_1237, sli0_1238, \
                         slh_738, slh_929, slh_930, sli1_1237, sli1_1238, \
                         smh_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = pb_x[k] * sli0_1237[k]
                    + f_14 * slh_929[k]
                    - f_10 * pc_x[k] * sli1_1237[k];

        t_1238[k] = pb_x[k] * sli0_1238[k]
                    + f_13 * slh_930[k]
                    - f_10 * pc_x[k] * sli1_1238[k];

        t_1239[k] = f_15 * slh_738[k]
                    + f_3 * pc_z[k] * smh_927[k];
    }

#pragma omp simd aligned(t_1240, t_1241, t_1242, pb_x, pc_x, pc_y, sli0_1241, sli0_1242, \
                         slh_933, slh_934, sli1_1241, sli1_1242, \
                         smh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1240[k] = f_3 * pc_y[k] * smh_929[k];

        t_1241[k] = pb_x[k] * sli0_1241[k]
                    + f_13 * slh_933[k]
                    - f_10 * pc_x[k] * sli1_1241[k];

        t_1242[k] = pb_x[k] * sli0_1242[k]
                    + f_12 * slh_934[k]
                    - f_10 * pc_x[k] * sli1_1242[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, pb_x, pc_x, pc_y, pc_z, sli0_1244, slh_741, \
                         slh_936, sli1_1244, smh_930, smh_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_15 * slh_741[k]
                    + f_3 * pc_z[k] * smh_930[k];

        t_1244[k] = pb_x[k] * sli0_1244[k]
                    + f_12 * slh_936[k]
                    - f_10 * pc_x[k] * sli1_1244[k];

        t_1245[k] = f_3 * pc_y[k] * smh_933[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, pb_x, pc_x, sli0_1246, slh_938, \
                         slh_939, slh_940, slh_941, sli1_1246, smh_939, smh_940, \
                         smh_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = pb_x[k] * sli0_1246[k]
                    + f_12 * slh_938[k]
                    - f_10 * pc_x[k] * sli1_1246[k];

        t_1247[k] = f_11 * slh_939[k]
                    + f_3 * pc_x[k] * smh_939[k];

        t_1248[k] = f_11 * slh_940[k]
                    + f_3 * pc_x[k] * smh_940[k];

        t_1249[k] = f_11 * slh_941[k]
                    + f_3 * pc_x[k] * smh_941[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pb_x, pc_x, sli0_1253, slh_942, \
                         slh_943, slh_944, sli1_1253, smh_942, smh_943, \
                         smh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_11 * slh_942[k]
                    + f_3 * pc_x[k] * smh_942[k];

        t_1251[k] = f_11 * slh_943[k]
                    + f_3 * pc_x[k] * smh_943[k];

        t_1252[k] = f_11 * slh_944[k]
                    + f_3 * pc_x[k] * smh_944[k];

        t_1253[k] = pb_x[k] * sli0_1253[k]
                    - f_10 * pc_x[k] * sli1_1253[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, t_1257, pb_x, pc_x, pc_z, sli0_1255, \
                         sli0_1256, sli0_1257, slh_750, sli1_1255, sli1_1256, sli1_1257, \
                         smh_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_15 * slh_750[k]
                    + f_3 * pc_z[k] * smh_939[k];

        t_1255[k] = pb_x[k] * sli0_1255[k]
                    - f_10 * pc_x[k] * sli1_1255[k];

        t_1256[k] = pb_x[k] * sli0_1256[k]
                    - f_10 * pc_x[k] * sli1_1256[k];

        t_1257[k] = pb_x[k] * sli0_1257[k]
                    - f_10 * pc_x[k] * sli1_1257[k];
    }

#pragma omp simd aligned(t_1258, t_1259, t_1260, t_1261, t_1262, pb_x, pc_x, pc_y, pc_z, \
                         sli0_1259, slh_756, sli1_1259, smg0_675, smg1_675, smh_944, \
                         smh_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1258[k] = f_3 * pc_y[k] * smh_944[k];

        t_1259[k] = pb_x[k] * sli0_1259[k]
                    - f_10 * pc_x[k] * sli1_1259[k];

        t_1260[k] = f_1 * smg0_675[k]
                    - f_2 * smg1_675[k]
                    + f_3 * pc_x[k] * smh_945[k];

        t_1261[k] = f_0 * slh_756[k]
                    + f_3 * pc_y[k] * smh_945[k];

        t_1262[k] = f_3 * pc_z[k] * smh_945[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sli0,
                                                           const size_t slh, const size_t sli1,
                                                           const size_t smg0, const size_t smg1,
                                                           const size_t smh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);
    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);
    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);
    auto *t_1296 = buffer.data(target + 1296);
    auto *t_1297 = buffer.data(target + 1297);
    auto *t_1298 = buffer.data(target + 1298);
    auto *t_1299 = buffer.data(target + 1299);
    auto *t_1300 = buffer.data(target + 1300);
    auto *t_1301 = buffer.data(target + 1301);
    auto *t_1302 = buffer.data(target + 1302);
    auto *t_1303 = buffer.data(target + 1303);
    auto *t_1304 = buffer.data(target + 1304);
    auto *t_1305 = buffer.data(target + 1305);
    auto *t_1306 = buffer.data(target + 1306);
    auto *t_1307 = buffer.data(target + 1307);
    auto *t_1308 = buffer.data(target + 1308);
    auto *t_1309 = buffer.data(target + 1309);
    auto *t_1310 = buffer.data(target + 1310);
    auto *t_1311 = buffer.data(target + 1311);
    auto *t_1312 = buffer.data(target + 1312);
    auto *t_1313 = buffer.data(target + 1313);
    auto *t_1314 = buffer.data(target + 1314);
    auto *t_1315 = buffer.data(target + 1315);
    auto *t_1316 = buffer.data(target + 1316);
    auto *t_1317 = buffer.data(target + 1317);
    auto *t_1318 = buffer.data(target + 1318);
    auto *t_1319 = buffer.data(target + 1319);
    auto *t_1320 = buffer.data(target + 1320);
    auto *t_1321 = buffer.data(target + 1321);
    auto *t_1322 = buffer.data(target + 1322);
    auto *t_1323 = buffer.data(target + 1323);
    auto *t_1324 = buffer.data(target + 1324);
    auto *t_1325 = buffer.data(target + 1325);
    auto *t_1326 = buffer.data(target + 1326);
    auto *t_1327 = buffer.data(target + 1327);
    auto *t_1328 = buffer.data(target + 1328);
    auto *t_1329 = buffer.data(target + 1329);
    auto *t_1330 = buffer.data(target + 1330);
    auto *t_1331 = buffer.data(target + 1331);
    auto *t_1332 = buffer.data(target + 1332);
    auto *t_1333 = buffer.data(target + 1333);
    auto *t_1334 = buffer.data(target + 1334);
    auto *t_1335 = buffer.data(target + 1335);
    auto *t_1336 = buffer.data(target + 1336);
    auto *t_1337 = buffer.data(target + 1337);
    auto *t_1338 = buffer.data(target + 1338);
    auto *t_1339 = buffer.data(target + 1339);
    auto *t_1340 = buffer.data(target + 1340);
    auto *t_1341 = buffer.data(target + 1341);
    auto *t_1342 = buffer.data(target + 1342);
    auto *t_1343 = buffer.data(target + 1343);
    auto *t_1344 = buffer.data(target + 1344);
    auto *t_1345 = buffer.data(target + 1345);
    auto *t_1346 = buffer.data(target + 1346);
    auto *t_1347 = buffer.data(target + 1347);
    auto *t_1348 = buffer.data(target + 1348);
    auto *t_1349 = buffer.data(target + 1349);
    auto *t_1350 = buffer.data(target + 1350);
    auto *t_1351 = buffer.data(target + 1351);
    auto *t_1352 = buffer.data(target + 1352);
    auto *t_1353 = buffer.data(target + 1353);
    auto *t_1354 = buffer.data(target + 1354);
    auto *t_1355 = buffer.data(target + 1355);
    auto *t_1356 = buffer.data(target + 1356);
    auto *t_1357 = buffer.data(target + 1357);
    auto *t_1358 = buffer.data(target + 1358);
    auto *t_1359 = buffer.data(target + 1359);
    auto *t_1360 = buffer.data(target + 1360);
    auto *t_1361 = buffer.data(target + 1361);
    auto *t_1362 = buffer.data(target + 1362);
    auto *t_1363 = buffer.data(target + 1363);
    auto *t_1364 = buffer.data(target + 1364);
    auto *t_1365 = buffer.data(target + 1365);
    auto *t_1366 = buffer.data(target + 1366);
    auto *t_1367 = buffer.data(target + 1367);
    auto *t_1368 = buffer.data(target + 1368);
    auto *t_1369 = buffer.data(target + 1369);
    auto *t_1370 = buffer.data(target + 1370);
    auto *t_1371 = buffer.data(target + 1371);
    auto *t_1372 = buffer.data(target + 1372);
    auto *t_1373 = buffer.data(target + 1373);
    auto *t_1374 = buffer.data(target + 1374);
    auto *t_1375 = buffer.data(target + 1375);
    auto *t_1376 = buffer.data(target + 1376);
    auto *t_1377 = buffer.data(target + 1377);
    auto *t_1378 = buffer.data(target + 1378);
    auto *t_1379 = buffer.data(target + 1379);
    auto *t_1380 = buffer.data(target + 1380);
    auto *t_1381 = buffer.data(target + 1381);
    auto *t_1382 = buffer.data(target + 1382);
    auto *t_1383 = buffer.data(target + 1383);
    auto *t_1384 = buffer.data(target + 1384);
    auto *t_1385 = buffer.data(target + 1385);
    auto *t_1386 = buffer.data(target + 1386);
    auto *t_1387 = buffer.data(target + 1387);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_1008 = buffer.data(sli0 + 1008);
    const auto *sli0_1011 = buffer.data(sli0 + 1011);
    const auto *sli0_1014 = buffer.data(sli0 + 1014);
    const auto *sli0_1018 = buffer.data(sli0 + 1018);
    const auto *sli0_1029 = buffer.data(sli0 + 1029);
    const auto *sli0_1031 = buffer.data(sli0 + 1031);
    const auto *sli0_1032 = buffer.data(sli0 + 1032);
    const auto *sli0_1033 = buffer.data(sli0 + 1033);

    const auto *slh_756 = buffer.data(slh + 756);
    const auto *slh_758 = buffer.data(slh + 758);
    const auto *slh_759 = buffer.data(slh + 759);
    const auto *slh_761 = buffer.data(slh + 761);
    const auto *slh_762 = buffer.data(slh + 762);
    const auto *slh_765 = buffer.data(slh + 765);
    const auto *slh_771 = buffer.data(slh + 771);
    const auto *slh_772 = buffer.data(slh + 772);
    const auto *slh_773 = buffer.data(slh + 773);
    const auto *slh_774 = buffer.data(slh + 774);
    const auto *slh_775 = buffer.data(slh + 775);
    const auto *slh_776 = buffer.data(slh + 776);
    const auto *slh_777 = buffer.data(slh + 777);
    const auto *slh_779 = buffer.data(slh + 779);
    const auto *slh_780 = buffer.data(slh + 780);
    const auto *slh_782 = buffer.data(slh + 782);
    const auto *slh_783 = buffer.data(slh + 783);
    const auto *slh_786 = buffer.data(slh + 786);
    const auto *slh_792 = buffer.data(slh + 792);
    const auto *slh_797 = buffer.data(slh + 797);
    const auto *slh_798 = buffer.data(slh + 798);
    const auto *slh_800 = buffer.data(slh + 800);
    const auto *slh_801 = buffer.data(slh + 801);
    const auto *slh_803 = buffer.data(slh + 803);
    const auto *slh_804 = buffer.data(slh + 804);
    const auto *slh_807 = buffer.data(slh + 807);
    const auto *slh_813 = buffer.data(slh + 813);
    const auto *slh_815 = buffer.data(slh + 815);
    const auto *slh_816 = buffer.data(slh + 816);
    const auto *slh_817 = buffer.data(slh + 817);
    const auto *slh_818 = buffer.data(slh + 818);
    const auto *slh_819 = buffer.data(slh + 819);
    const auto *slh_821 = buffer.data(slh + 821);
    const auto *slh_822 = buffer.data(slh + 822);
    const auto *slh_824 = buffer.data(slh + 824);
    const auto *slh_825 = buffer.data(slh + 825);
    const auto *slh_828 = buffer.data(slh + 828);
    const auto *slh_834 = buffer.data(slh + 834);
    const auto *slh_836 = buffer.data(slh + 836);
    const auto *slh_837 = buffer.data(slh + 837);
    const auto *slh_838 = buffer.data(slh + 838);
    const auto *slh_839 = buffer.data(slh + 839);
    const auto *slh_840 = buffer.data(slh + 840);
    const auto *slh_842 = buffer.data(slh + 842);
    const auto *slh_845 = buffer.data(slh + 845);
    const auto *slh_849 = buffer.data(slh + 849);

    const auto *sli1_1008 = buffer.data(sli1 + 1008);
    const auto *sli1_1011 = buffer.data(sli1 + 1011);
    const auto *sli1_1014 = buffer.data(sli1 + 1014);
    const auto *sli1_1018 = buffer.data(sli1 + 1018);
    const auto *sli1_1029 = buffer.data(sli1 + 1029);
    const auto *sli1_1031 = buffer.data(sli1 + 1031);
    const auto *sli1_1032 = buffer.data(sli1 + 1032);
    const auto *sli1_1033 = buffer.data(sli1 + 1033);

    const auto *smg0_678 = buffer.data(smg0 + 678);
    const auto *smg0_680 = buffer.data(smg0 + 680);
    const auto *smg0_681 = buffer.data(smg0 + 681);
    const auto *smg0_684 = buffer.data(smg0 + 684);
    const auto *smg0_685 = buffer.data(smg0 + 685);
    const auto *smg0_687 = buffer.data(smg0 + 687);
    const auto *smg0_688 = buffer.data(smg0 + 688);
    const auto *smg0_689 = buffer.data(smg0 + 689);
    const auto *smg0_695 = buffer.data(smg0 + 695);
    const auto *smg0_699 = buffer.data(smg0 + 699);
    const auto *smg0_702 = buffer.data(smg0 + 702);
    const auto *smg0_704 = buffer.data(smg0 + 704);
    const auto *smg0_705 = buffer.data(smg0 + 705);
    const auto *smg0_708 = buffer.data(smg0 + 708);
    const auto *smg0_710 = buffer.data(smg0 + 710);
    const auto *smg0_711 = buffer.data(smg0 + 711);
    const auto *smg0_714 = buffer.data(smg0 + 714);
    const auto *smg0_715 = buffer.data(smg0 + 715);
    const auto *smg0_717 = buffer.data(smg0 + 717);
    const auto *smg0_718 = buffer.data(smg0 + 718);
    const auto *smg0_719 = buffer.data(smg0 + 719);
    const auto *smg0_720 = buffer.data(smg0 + 720);
    const auto *smg0_723 = buffer.data(smg0 + 723);
    const auto *smg0_725 = buffer.data(smg0 + 725);
    const auto *smg0_726 = buffer.data(smg0 + 726);
    const auto *smg0_729 = buffer.data(smg0 + 729);
    const auto *smg0_730 = buffer.data(smg0 + 730);
    const auto *smg0_732 = buffer.data(smg0 + 732);
    const auto *smg0_733 = buffer.data(smg0 + 733);
    const auto *smg0_734 = buffer.data(smg0 + 734);
    const auto *smg0_735 = buffer.data(smg0 + 735);
    const auto *smg0_738 = buffer.data(smg0 + 738);
    const auto *smg0_740 = buffer.data(smg0 + 740);
    const auto *smg0_741 = buffer.data(smg0 + 741);
    const auto *smg0_744 = buffer.data(smg0 + 744);
    const auto *smg0_745 = buffer.data(smg0 + 745);
    const auto *smg0_747 = buffer.data(smg0 + 747);
    const auto *smg0_749 = buffer.data(smg0 + 749);

    const auto *smg1_678 = buffer.data(smg1 + 678);
    const auto *smg1_680 = buffer.data(smg1 + 680);
    const auto *smg1_681 = buffer.data(smg1 + 681);
    const auto *smg1_684 = buffer.data(smg1 + 684);
    const auto *smg1_685 = buffer.data(smg1 + 685);
    const auto *smg1_687 = buffer.data(smg1 + 687);
    const auto *smg1_688 = buffer.data(smg1 + 688);
    const auto *smg1_689 = buffer.data(smg1 + 689);
    const auto *smg1_695 = buffer.data(smg1 + 695);
    const auto *smg1_699 = buffer.data(smg1 + 699);
    const auto *smg1_702 = buffer.data(smg1 + 702);
    const auto *smg1_704 = buffer.data(smg1 + 704);
    const auto *smg1_705 = buffer.data(smg1 + 705);
    const auto *smg1_708 = buffer.data(smg1 + 708);
    const auto *smg1_710 = buffer.data(smg1 + 710);
    const auto *smg1_711 = buffer.data(smg1 + 711);
    const auto *smg1_714 = buffer.data(smg1 + 714);
    const auto *smg1_715 = buffer.data(smg1 + 715);
    const auto *smg1_717 = buffer.data(smg1 + 717);
    const auto *smg1_718 = buffer.data(smg1 + 718);
    const auto *smg1_719 = buffer.data(smg1 + 719);
    const auto *smg1_720 = buffer.data(smg1 + 720);
    const auto *smg1_723 = buffer.data(smg1 + 723);
    const auto *smg1_725 = buffer.data(smg1 + 725);
    const auto *smg1_726 = buffer.data(smg1 + 726);
    const auto *smg1_729 = buffer.data(smg1 + 729);
    const auto *smg1_730 = buffer.data(smg1 + 730);
    const auto *smg1_732 = buffer.data(smg1 + 732);
    const auto *smg1_733 = buffer.data(smg1 + 733);
    const auto *smg1_734 = buffer.data(smg1 + 734);
    const auto *smg1_735 = buffer.data(smg1 + 735);
    const auto *smg1_738 = buffer.data(smg1 + 738);
    const auto *smg1_740 = buffer.data(smg1 + 740);
    const auto *smg1_741 = buffer.data(smg1 + 741);
    const auto *smg1_744 = buffer.data(smg1 + 744);
    const auto *smg1_745 = buffer.data(smg1 + 745);
    const auto *smg1_747 = buffer.data(smg1 + 747);
    const auto *smg1_749 = buffer.data(smg1 + 749);

    const auto *smh_947 = buffer.data(smh + 947);
    const auto *smh_948 = buffer.data(smh + 948);
    const auto *smh_950 = buffer.data(smh + 950);
    const auto *smh_951 = buffer.data(smh + 951);
    const auto *smh_954 = buffer.data(smh + 954);
    const auto *smh_955 = buffer.data(smh + 955);
    const auto *smh_957 = buffer.data(smh + 957);
    const auto *smh_959 = buffer.data(smh + 959);
    const auto *smh_960 = buffer.data(smh + 960);
    const auto *smh_961 = buffer.data(smh + 961);
    const auto *smh_962 = buffer.data(smh + 962);
    const auto *smh_963 = buffer.data(smh + 963);
    const auto *smh_964 = buffer.data(smh + 964);
    const auto *smh_965 = buffer.data(smh + 965);
    const auto *smh_966 = buffer.data(smh + 966);
    const auto *smh_968 = buffer.data(smh + 968);
    const auto *smh_969 = buffer.data(smh + 969);
    const auto *smh_971 = buffer.data(smh + 971);
    const auto *smh_972 = buffer.data(smh + 972);
    const auto *smh_975 = buffer.data(smh + 975);
    const auto *smh_978 = buffer.data(smh + 978);
    const auto *smh_980 = buffer.data(smh + 980);
    const auto *smh_981 = buffer.data(smh + 981);
    const auto *smh_982 = buffer.data(smh + 982);
    const auto *smh_983 = buffer.data(smh + 983);
    const auto *smh_984 = buffer.data(smh + 984);
    const auto *smh_985 = buffer.data(smh + 985);
    const auto *smh_986 = buffer.data(smh + 986);
    const auto *smh_987 = buffer.data(smh + 987);
    const auto *smh_989 = buffer.data(smh + 989);
    const auto *smh_990 = buffer.data(smh + 990);
    const auto *smh_992 = buffer.data(smh + 992);
    const auto *smh_993 = buffer.data(smh + 993);
    const auto *smh_996 = buffer.data(smh + 996);
    const auto *smh_997 = buffer.data(smh + 997);
    const auto *smh_999 = buffer.data(smh + 999);
    const auto *smh_1001 = buffer.data(smh + 1001);
    const auto *smh_1002 = buffer.data(smh + 1002);
    const auto *smh_1003 = buffer.data(smh + 1003);
    const auto *smh_1004 = buffer.data(smh + 1004);
    const auto *smh_1005 = buffer.data(smh + 1005);
    const auto *smh_1006 = buffer.data(smh + 1006);
    const auto *smh_1007 = buffer.data(smh + 1007);
    const auto *smh_1008 = buffer.data(smh + 1008);
    const auto *smh_1010 = buffer.data(smh + 1010);
    const auto *smh_1011 = buffer.data(smh + 1011);
    const auto *smh_1013 = buffer.data(smh + 1013);
    const auto *smh_1014 = buffer.data(smh + 1014);
    const auto *smh_1017 = buffer.data(smh + 1017);
    const auto *smh_1018 = buffer.data(smh + 1018);
    const auto *smh_1020 = buffer.data(smh + 1020);
    const auto *smh_1022 = buffer.data(smh + 1022);
    const auto *smh_1023 = buffer.data(smh + 1023);
    const auto *smh_1024 = buffer.data(smh + 1024);
    const auto *smh_1025 = buffer.data(smh + 1025);
    const auto *smh_1026 = buffer.data(smh + 1026);
    const auto *smh_1027 = buffer.data(smh + 1027);
    const auto *smh_1028 = buffer.data(smh + 1028);
    const auto *smh_1029 = buffer.data(smh + 1029);
    const auto *smh_1031 = buffer.data(smh + 1031);
    const auto *smh_1032 = buffer.data(smh + 1032);
    const auto *smh_1034 = buffer.data(smh + 1034);
    const auto *smh_1035 = buffer.data(smh + 1035);
    const auto *smh_1038 = buffer.data(smh + 1038);
    const auto *smh_1039 = buffer.data(smh + 1039);
    const auto *smh_1041 = buffer.data(smh + 1041);
    const auto *smh_1043 = buffer.data(smh + 1043);
    const auto *smh_1044 = buffer.data(smh + 1044);

#pragma omp simd aligned(t_1263, t_1264, t_1265, pc_x, pc_y, slh_758, smg0_678, smg0_680, \
                         smg1_678, smg1_680, smh_947, smh_948, \
                         smh_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = f_4 * smg0_678[k]
                    - f_5 * smg1_678[k]
                    + f_3 * pc_x[k] * smh_948[k];

        t_1264[k] = f_0 * slh_758[k]
                    + f_3 * pc_y[k] * smh_947[k];

        t_1265[k] = f_4 * smg0_680[k]
                    - f_5 * smg1_680[k]
                    + f_3 * pc_x[k] * smh_950[k];
    }

#pragma omp simd aligned(t_1266, t_1267, t_1268, t_1269, pc_x, pc_y, pc_z, slh_761, smg0_681, \
                         smg0_684, smg1_681, smg1_684, smh_948, smh_950, smh_951, \
                         smh_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1266[k] = f_6 * smg0_681[k]
                    - f_7 * smg1_681[k]
                    + f_3 * pc_x[k] * smh_951[k];

        t_1267[k] = f_3 * pc_z[k] * smh_948[k];

        t_1268[k] = f_0 * slh_761[k]
                    + f_3 * pc_y[k] * smh_950[k];

        t_1269[k] = f_6 * smg0_684[k]
                    - f_7 * smg1_684[k]
                    + f_3 * pc_x[k] * smh_954[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, t_1273, pc_x, pc_y, pc_z, slh_765, smg0_685, \
                         smg0_687, smg1_685, smg1_687, smh_951, smh_954, smh_955, \
                         smh_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = f_8 * smg0_685[k]
                    - f_9 * smg1_685[k]
                    + f_3 * pc_x[k] * smh_955[k];

        t_1271[k] = f_3 * pc_z[k] * smh_951[k];

        t_1272[k] = f_8 * smg0_687[k]
                    - f_9 * smg1_687[k]
                    + f_3 * pc_x[k] * smh_957[k];

        t_1273[k] = f_0 * slh_765[k]
                    + f_3 * pc_y[k] * smh_954[k];
    }

#pragma omp simd aligned(t_1274, t_1275, t_1276, t_1277, t_1278, t_1279, pc_x, smg0_689, \
                         smg1_689, smh_959, smh_960, smh_961, smh_962, smh_963, \
                         smh_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1274[k] = f_8 * smg0_689[k]
                    - f_9 * smg1_689[k]
                    + f_3 * pc_x[k] * smh_959[k];

        t_1275[k] = f_3 * pc_x[k] * smh_960[k];

        t_1276[k] = f_3 * pc_x[k] * smh_961[k];

        t_1277[k] = f_3 * pc_x[k] * smh_962[k];

        t_1278[k] = f_3 * pc_x[k] * smh_963[k];

        t_1279[k] = f_3 * pc_x[k] * smh_964[k];
    }

#pragma omp simd aligned(t_1280, t_1281, t_1282, t_1283, pc_x, pc_y, pc_z, slh_771, slh_773, \
                         smg0_685, smg0_687, smg1_685, smg1_687, smh_960, smh_962, \
                         smh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1280[k] = f_3 * pc_x[k] * smh_965[k];

        t_1281[k] = f_0 * slh_771[k]
                    + f_1 * smg0_685[k]
                    - f_2 * smg1_685[k]
                    + f_3 * pc_y[k] * smh_960[k];

        t_1282[k] = f_3 * pc_z[k] * smh_960[k];

        t_1283[k] = f_0 * slh_773[k]
                    + f_4 * smg0_687[k]
                    - f_5 * smg1_687[k]
                    + f_3 * pc_y[k] * smh_962[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, t_1287, pc_y, pc_z, slh_774, slh_775, \
                         slh_776, smg0_688, smg0_689, smg1_688, smg1_689, smh_963, smh_964, \
                         smh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_0 * slh_774[k]
                    + f_6 * smg0_688[k]
                    - f_7 * smg1_688[k]
                    + f_3 * pc_y[k] * smh_963[k];

        t_1285[k] = f_0 * slh_775[k]
                    + f_8 * smg0_689[k]
                    - f_9 * smg1_689[k]
                    + f_3 * pc_y[k] * smh_964[k];

        t_1286[k] = f_0 * slh_776[k]
                    + f_3 * pc_y[k] * smh_965[k];

        t_1287[k] = f_1 * smg0_689[k]
                    - f_2 * smg1_689[k]
                    + f_3 * pc_z[k] * smh_965[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pb_z, pc_y, pc_z, sli0_1008, \
                         sli0_1011, slh_756, slh_777, sli1_1008, sli1_1011, \
                         smh_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = pb_z[k] * sli0_1008[k]
                    - f_10 * pc_z[k] * sli1_1008[k];

        t_1289[k] = f_15 * slh_777[k]
                    + f_3 * pc_y[k] * smh_966[k];

        t_1290[k] = f_11 * slh_756[k]
                    + f_3 * pc_z[k] * smh_966[k];

        t_1291[k] = pb_z[k] * sli0_1011[k]
                    - f_10 * pc_z[k] * sli1_1011[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, pb_z, pc_x, pc_y, pc_z, sli0_1014, slh_779, \
                         sli1_1014, smg0_695, smg1_695, smh_968, \
                         smh_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_15 * slh_779[k]
                    + f_3 * pc_y[k] * smh_968[k];

        t_1293[k] = f_4 * smg0_695[k]
                    - f_5 * smg1_695[k]
                    + f_3 * pc_x[k] * smh_971[k];

        t_1294[k] = pb_z[k] * sli0_1014[k]
                    - f_10 * pc_z[k] * sli1_1014[k];
    }

#pragma omp simd aligned(t_1295, t_1296, t_1297, pc_x, pc_y, pc_z, slh_759, slh_782, smg0_699, \
                         smg1_699, smh_969, smh_971, smh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1295[k] = f_11 * slh_759[k]
                    + f_3 * pc_z[k] * smh_969[k];

        t_1296[k] = f_15 * slh_782[k]
                    + f_3 * pc_y[k] * smh_971[k];

        t_1297[k] = f_6 * smg0_699[k]
                    - f_7 * smg1_699[k]
                    + f_3 * pc_x[k] * smh_975[k];
    }

#pragma omp simd aligned(t_1298, t_1299, t_1300, pb_z, pc_x, pc_z, sli0_1018, slh_762, \
                         sli1_1018, smg0_702, smg1_702, smh_972, \
                         smh_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1298[k] = pb_z[k] * sli0_1018[k]
                    - f_10 * pc_z[k] * sli1_1018[k];

        t_1299[k] = f_11 * slh_762[k]
                    + f_3 * pc_z[k] * smh_972[k];

        t_1300[k] = f_8 * smg0_702[k]
                    - f_9 * smg1_702[k]
                    + f_3 * pc_x[k] * smh_978[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, t_1305, pc_x, pc_y, slh_786, \
                         smg0_704, smg1_704, smh_975, smh_980, smh_981, smh_982, \
                         smh_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_15 * slh_786[k]
                    + f_3 * pc_y[k] * smh_975[k];

        t_1302[k] = f_8 * smg0_704[k]
                    - f_9 * smg1_704[k]
                    + f_3 * pc_x[k] * smh_980[k];

        t_1303[k] = f_3 * pc_x[k] * smh_981[k];

        t_1304[k] = f_3 * pc_x[k] * smh_982[k];

        t_1305[k] = f_3 * pc_x[k] * smh_983[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, t_1309, t_1310, pb_z, pc_x, pc_z, sli0_1029, \
                         slh_771, sli1_1029, smh_981, smh_984, smh_985, \
                         smh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_3 * pc_x[k] * smh_984[k];

        t_1307[k] = f_3 * pc_x[k] * smh_985[k];

        t_1308[k] = f_3 * pc_x[k] * smh_986[k];

        t_1309[k] = pb_z[k] * sli0_1029[k]
                    - f_10 * pc_z[k] * sli1_1029[k];

        t_1310[k] = f_11 * slh_771[k]
                    + f_3 * pc_z[k] * smh_981[k];
    }

#pragma omp simd aligned(t_1311, t_1312, t_1313, pb_z, pc_z, sli0_1031, sli0_1032, sli0_1033, \
                         slh_772, slh_773, slh_774, sli1_1031, sli1_1032, \
                         sli1_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1311[k] = pb_z[k] * sli0_1031[k]
                    + f_12 * slh_772[k]
                    - f_10 * pc_z[k] * sli1_1031[k];

        t_1312[k] = pb_z[k] * sli0_1032[k]
                    + f_13 * slh_773[k]
                    - f_10 * pc_z[k] * sli1_1032[k];

        t_1313[k] = pb_z[k] * sli0_1033[k]
                    + f_14 * slh_774[k]
                    - f_10 * pc_z[k] * sli1_1033[k];
    }

#pragma omp simd aligned(t_1314, t_1315, t_1316, t_1317, pc_x, pc_y, pc_z, slh_776, slh_797, \
                         slh_798, smg0_704, smg0_705, smg1_704, smg1_705, smh_986, \
                         smh_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1314[k] = f_15 * slh_797[k]
                    + f_3 * pc_y[k] * smh_986[k];

        t_1315[k] = f_11 * slh_776[k]
                    + f_1 * smg0_704[k]
                    - f_2 * smg1_704[k]
                    + f_3 * pc_z[k] * smh_986[k];

        t_1316[k] = f_1 * smg0_705[k]
                    - f_2 * smg1_705[k]
                    + f_3 * pc_x[k] * smh_987[k];

        t_1317[k] = f_16 * slh_798[k]
                    + f_3 * pc_y[k] * smh_987[k];
    }

#pragma omp simd aligned(t_1318, t_1319, t_1320, pc_x, pc_y, pc_z, slh_777, slh_800, smg0_708, \
                         smg1_708, smh_987, smh_989, smh_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1318[k] = f_12 * slh_777[k]
                    + f_3 * pc_z[k] * smh_987[k];

        t_1319[k] = f_4 * smg0_708[k]
                    - f_5 * smg1_708[k]
                    + f_3 * pc_x[k] * smh_990[k];

        t_1320[k] = f_16 * slh_800[k]
                    + f_3 * pc_y[k] * smh_989[k];
    }

#pragma omp simd aligned(t_1321, t_1322, t_1323, t_1324, pc_x, pc_y, pc_z, slh_780, slh_803, \
                         smg0_710, smg0_711, smg1_710, smg1_711, smh_990, smh_992, \
                         smh_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1321[k] = f_4 * smg0_710[k]
                    - f_5 * smg1_710[k]
                    + f_3 * pc_x[k] * smh_992[k];

        t_1322[k] = f_6 * smg0_711[k]
                    - f_7 * smg1_711[k]
                    + f_3 * pc_x[k] * smh_993[k];

        t_1323[k] = f_12 * slh_780[k]
                    + f_3 * pc_z[k] * smh_990[k];

        t_1324[k] = f_16 * slh_803[k]
                    + f_3 * pc_y[k] * smh_992[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pc_x, pc_z, slh_783, smg0_714, smg0_715, \
                         smg1_714, smg1_715, smh_993, smh_996, \
                         smh_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_6 * smg0_714[k]
                    - f_7 * smg1_714[k]
                    + f_3 * pc_x[k] * smh_996[k];

        t_1326[k] = f_8 * smg0_715[k]
                    - f_9 * smg1_715[k]
                    + f_3 * pc_x[k] * smh_997[k];

        t_1327[k] = f_12 * slh_783[k]
                    + f_3 * pc_z[k] * smh_993[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, t_1331, pc_x, pc_y, slh_807, smg0_717, \
                         smg0_719, smg1_717, smg1_719, smh_996, smh_999, smh_1001, \
                         smh_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_8 * smg0_717[k]
                    - f_9 * smg1_717[k]
                    + f_3 * pc_x[k] * smh_999[k];

        t_1329[k] = f_16 * slh_807[k]
                    + f_3 * pc_y[k] * smh_996[k];

        t_1330[k] = f_8 * smg0_719[k]
                    - f_9 * smg1_719[k]
                    + f_3 * pc_x[k] * smh_1001[k];

        t_1331[k] = f_3 * pc_x[k] * smh_1002[k];
    }

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, t_1336, pc_x, smh_1003, smh_1004, \
                         smh_1005, smh_1006, smh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = f_3 * pc_x[k] * smh_1003[k];

        t_1333[k] = f_3 * pc_x[k] * smh_1004[k];

        t_1334[k] = f_3 * pc_x[k] * smh_1005[k];

        t_1335[k] = f_3 * pc_x[k] * smh_1006[k];

        t_1336[k] = f_3 * pc_x[k] * smh_1007[k];
    }

#pragma omp simd aligned(t_1337, t_1338, t_1339, pc_y, pc_z, slh_792, slh_813, slh_815, \
                         smg0_715, smg0_717, smg1_715, smg1_717, smh_1002, \
                         smh_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1337[k] = f_16 * slh_813[k]
                    + f_1 * smg0_715[k]
                    - f_2 * smg1_715[k]
                    + f_3 * pc_y[k] * smh_1002[k];

        t_1338[k] = f_12 * slh_792[k]
                    + f_3 * pc_z[k] * smh_1002[k];

        t_1339[k] = f_16 * slh_815[k]
                    + f_4 * smg0_717[k]
                    - f_5 * smg1_717[k]
                    + f_3 * pc_y[k] * smh_1004[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, pc_y, slh_816, slh_817, slh_818, smg0_718, \
                         smg0_719, smg1_718, smg1_719, smh_1005, smh_1006, \
                         smh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_16 * slh_816[k]
                    + f_6 * smg0_718[k]
                    - f_7 * smg1_718[k]
                    + f_3 * pc_y[k] * smh_1005[k];

        t_1341[k] = f_16 * slh_817[k]
                    + f_8 * smg0_719[k]
                    - f_9 * smg1_719[k]
                    + f_3 * pc_y[k] * smh_1006[k];

        t_1342[k] = f_16 * slh_818[k]
                    + f_3 * pc_y[k] * smh_1007[k];
    }

#pragma omp simd aligned(t_1343, t_1344, t_1345, t_1346, pc_x, pc_y, pc_z, slh_797, slh_798, \
                         slh_819, smg0_719, smg0_720, smg1_719, smg1_720, smh_1007, \
                         smh_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1343[k] = f_12 * slh_797[k]
                    + f_1 * smg0_719[k]
                    - f_2 * smg1_719[k]
                    + f_3 * pc_z[k] * smh_1007[k];

        t_1344[k] = f_1 * smg0_720[k]
                    - f_2 * smg1_720[k]
                    + f_3 * pc_x[k] * smh_1008[k];

        t_1345[k] = f_17 * slh_819[k]
                    + f_3 * pc_y[k] * smh_1008[k];

        t_1346[k] = f_13 * slh_798[k]
                    + f_3 * pc_z[k] * smh_1008[k];
    }

#pragma omp simd aligned(t_1347, t_1348, t_1349, pc_x, pc_y, slh_821, smg0_723, smg0_725, \
                         smg1_723, smg1_725, smh_1010, smh_1011, \
                         smh_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1347[k] = f_4 * smg0_723[k]
                    - f_5 * smg1_723[k]
                    + f_3 * pc_x[k] * smh_1011[k];

        t_1348[k] = f_17 * slh_821[k]
                    + f_3 * pc_y[k] * smh_1010[k];

        t_1349[k] = f_4 * smg0_725[k]
                    - f_5 * smg1_725[k]
                    + f_3 * pc_x[k] * smh_1013[k];
    }

#pragma omp simd aligned(t_1350, t_1351, t_1352, pc_x, pc_y, pc_z, slh_801, slh_824, smg0_726, \
                         smg1_726, smh_1011, smh_1013, smh_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1350[k] = f_6 * smg0_726[k]
                    - f_7 * smg1_726[k]
                    + f_3 * pc_x[k] * smh_1014[k];

        t_1351[k] = f_13 * slh_801[k]
                    + f_3 * pc_z[k] * smh_1011[k];

        t_1352[k] = f_17 * slh_824[k]
                    + f_3 * pc_y[k] * smh_1013[k];
    }

#pragma omp simd aligned(t_1353, t_1354, t_1355, pc_x, pc_z, slh_804, smg0_729, smg0_730, \
                         smg1_729, smg1_730, smh_1014, smh_1017, \
                         smh_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1353[k] = f_6 * smg0_729[k]
                    - f_7 * smg1_729[k]
                    + f_3 * pc_x[k] * smh_1017[k];

        t_1354[k] = f_8 * smg0_730[k]
                    - f_9 * smg1_730[k]
                    + f_3 * pc_x[k] * smh_1018[k];

        t_1355[k] = f_13 * slh_804[k]
                    + f_3 * pc_z[k] * smh_1014[k];
    }

#pragma omp simd aligned(t_1356, t_1357, t_1358, t_1359, pc_x, pc_y, slh_828, smg0_732, \
                         smg0_734, smg1_732, smg1_734, smh_1017, smh_1020, smh_1022, \
                         smh_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1356[k] = f_8 * smg0_732[k]
                    - f_9 * smg1_732[k]
                    + f_3 * pc_x[k] * smh_1020[k];

        t_1357[k] = f_17 * slh_828[k]
                    + f_3 * pc_y[k] * smh_1017[k];

        t_1358[k] = f_8 * smg0_734[k]
                    - f_9 * smg1_734[k]
                    + f_3 * pc_x[k] * smh_1022[k];

        t_1359[k] = f_3 * pc_x[k] * smh_1023[k];
    }

#pragma omp simd aligned(t_1360, t_1361, t_1362, t_1363, t_1364, pc_x, smh_1024, smh_1025, \
                         smh_1026, smh_1027, smh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1360[k] = f_3 * pc_x[k] * smh_1024[k];

        t_1361[k] = f_3 * pc_x[k] * smh_1025[k];

        t_1362[k] = f_3 * pc_x[k] * smh_1026[k];

        t_1363[k] = f_3 * pc_x[k] * smh_1027[k];

        t_1364[k] = f_3 * pc_x[k] * smh_1028[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, pc_y, pc_z, slh_813, slh_834, slh_836, \
                         smg0_730, smg0_732, smg1_730, smg1_732, smh_1023, \
                         smh_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = f_17 * slh_834[k]
                    + f_1 * smg0_730[k]
                    - f_2 * smg1_730[k]
                    + f_3 * pc_y[k] * smh_1023[k];

        t_1366[k] = f_13 * slh_813[k]
                    + f_3 * pc_z[k] * smh_1023[k];

        t_1367[k] = f_17 * slh_836[k]
                    + f_4 * smg0_732[k]
                    - f_5 * smg1_732[k]
                    + f_3 * pc_y[k] * smh_1025[k];
    }

#pragma omp simd aligned(t_1368, t_1369, t_1370, pc_y, slh_837, slh_838, slh_839, smg0_733, \
                         smg0_734, smg1_733, smg1_734, smh_1026, smh_1027, \
                         smh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1368[k] = f_17 * slh_837[k]
                    + f_6 * smg0_733[k]
                    - f_7 * smg1_733[k]
                    + f_3 * pc_y[k] * smh_1026[k];

        t_1369[k] = f_17 * slh_838[k]
                    + f_8 * smg0_734[k]
                    - f_9 * smg1_734[k]
                    + f_3 * pc_y[k] * smh_1027[k];

        t_1370[k] = f_17 * slh_839[k]
                    + f_3 * pc_y[k] * smh_1028[k];
    }

#pragma omp simd aligned(t_1371, t_1372, t_1373, t_1374, pc_x, pc_y, pc_z, slh_818, slh_819, \
                         slh_840, smg0_734, smg0_735, smg1_734, smg1_735, smh_1028, \
                         smh_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1371[k] = f_13 * slh_818[k]
                    + f_1 * smg0_734[k]
                    - f_2 * smg1_734[k]
                    + f_3 * pc_z[k] * smh_1028[k];

        t_1372[k] = f_1 * smg0_735[k]
                    - f_2 * smg1_735[k]
                    + f_3 * pc_x[k] * smh_1029[k];

        t_1373[k] = f_18 * slh_840[k]
                    + f_3 * pc_y[k] * smh_1029[k];

        t_1374[k] = f_14 * slh_819[k]
                    + f_3 * pc_z[k] * smh_1029[k];
    }

#pragma omp simd aligned(t_1375, t_1376, t_1377, pc_x, pc_y, slh_842, smg0_738, smg0_740, \
                         smg1_738, smg1_740, smh_1031, smh_1032, \
                         smh_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1375[k] = f_4 * smg0_738[k]
                    - f_5 * smg1_738[k]
                    + f_3 * pc_x[k] * smh_1032[k];

        t_1376[k] = f_18 * slh_842[k]
                    + f_3 * pc_y[k] * smh_1031[k];

        t_1377[k] = f_4 * smg0_740[k]
                    - f_5 * smg1_740[k]
                    + f_3 * pc_x[k] * smh_1034[k];
    }

#pragma omp simd aligned(t_1378, t_1379, t_1380, pc_x, pc_y, pc_z, slh_822, slh_845, smg0_741, \
                         smg1_741, smh_1032, smh_1034, smh_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1378[k] = f_6 * smg0_741[k]
                    - f_7 * smg1_741[k]
                    + f_3 * pc_x[k] * smh_1035[k];

        t_1379[k] = f_14 * slh_822[k]
                    + f_3 * pc_z[k] * smh_1032[k];

        t_1380[k] = f_18 * slh_845[k]
                    + f_3 * pc_y[k] * smh_1034[k];
    }

#pragma omp simd aligned(t_1381, t_1382, t_1383, pc_x, pc_z, slh_825, smg0_744, smg0_745, \
                         smg1_744, smg1_745, smh_1035, smh_1038, \
                         smh_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1381[k] = f_6 * smg0_744[k]
                    - f_7 * smg1_744[k]
                    + f_3 * pc_x[k] * smh_1038[k];

        t_1382[k] = f_8 * smg0_745[k]
                    - f_9 * smg1_745[k]
                    + f_3 * pc_x[k] * smh_1039[k];

        t_1383[k] = f_14 * slh_825[k]
                    + f_3 * pc_z[k] * smh_1035[k];
    }

#pragma omp simd aligned(t_1384, t_1385, t_1386, t_1387, pc_x, pc_y, slh_849, smg0_747, \
                         smg0_749, smg1_747, smg1_749, smh_1038, smh_1041, smh_1043, \
                         smh_1044 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1384[k] = f_8 * smg0_747[k]
                    - f_9 * smg1_747[k]
                    + f_3 * pc_x[k] * smh_1041[k];

        t_1385[k] = f_18 * slh_849[k]
                    + f_3 * pc_y[k] * smh_1038[k];

        t_1386[k] = f_8 * smg0_749[k]
                    - f_9 * smg1_749[k]
                    + f_3 * pc_x[k] * smh_1043[k];

        t_1387[k] = f_3 * pc_x[k] * smh_1044[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sli0,
                                                           const size_t slh, const size_t sli1,
                                                           const size_t smg0, const size_t smg1,
                                                           const size_t smh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.0 / q;
    const auto f_16 = 3.5 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1388 = buffer.data(target + 1388);
    auto *t_1389 = buffer.data(target + 1389);
    auto *t_1390 = buffer.data(target + 1390);
    auto *t_1391 = buffer.data(target + 1391);
    auto *t_1392 = buffer.data(target + 1392);
    auto *t_1393 = buffer.data(target + 1393);
    auto *t_1394 = buffer.data(target + 1394);
    auto *t_1395 = buffer.data(target + 1395);
    auto *t_1396 = buffer.data(target + 1396);
    auto *t_1397 = buffer.data(target + 1397);
    auto *t_1398 = buffer.data(target + 1398);
    auto *t_1399 = buffer.data(target + 1399);
    auto *t_1400 = buffer.data(target + 1400);
    auto *t_1401 = buffer.data(target + 1401);
    auto *t_1402 = buffer.data(target + 1402);
    auto *t_1403 = buffer.data(target + 1403);
    auto *t_1404 = buffer.data(target + 1404);
    auto *t_1405 = buffer.data(target + 1405);
    auto *t_1406 = buffer.data(target + 1406);
    auto *t_1407 = buffer.data(target + 1407);
    auto *t_1408 = buffer.data(target + 1408);
    auto *t_1409 = buffer.data(target + 1409);
    auto *t_1410 = buffer.data(target + 1410);
    auto *t_1411 = buffer.data(target + 1411);
    auto *t_1412 = buffer.data(target + 1412);
    auto *t_1413 = buffer.data(target + 1413);
    auto *t_1414 = buffer.data(target + 1414);
    auto *t_1415 = buffer.data(target + 1415);
    auto *t_1416 = buffer.data(target + 1416);
    auto *t_1417 = buffer.data(target + 1417);
    auto *t_1418 = buffer.data(target + 1418);
    auto *t_1419 = buffer.data(target + 1419);
    auto *t_1420 = buffer.data(target + 1420);
    auto *t_1421 = buffer.data(target + 1421);
    auto *t_1422 = buffer.data(target + 1422);
    auto *t_1423 = buffer.data(target + 1423);
    auto *t_1424 = buffer.data(target + 1424);
    auto *t_1425 = buffer.data(target + 1425);
    auto *t_1426 = buffer.data(target + 1426);
    auto *t_1427 = buffer.data(target + 1427);
    auto *t_1428 = buffer.data(target + 1428);
    auto *t_1429 = buffer.data(target + 1429);
    auto *t_1430 = buffer.data(target + 1430);
    auto *t_1431 = buffer.data(target + 1431);
    auto *t_1432 = buffer.data(target + 1432);
    auto *t_1433 = buffer.data(target + 1433);
    auto *t_1434 = buffer.data(target + 1434);
    auto *t_1435 = buffer.data(target + 1435);
    auto *t_1436 = buffer.data(target + 1436);
    auto *t_1437 = buffer.data(target + 1437);
    auto *t_1438 = buffer.data(target + 1438);
    auto *t_1439 = buffer.data(target + 1439);
    auto *t_1440 = buffer.data(target + 1440);
    auto *t_1441 = buffer.data(target + 1441);
    auto *t_1442 = buffer.data(target + 1442);
    auto *t_1443 = buffer.data(target + 1443);
    auto *t_1444 = buffer.data(target + 1444);
    auto *t_1445 = buffer.data(target + 1445);
    auto *t_1446 = buffer.data(target + 1446);
    auto *t_1447 = buffer.data(target + 1447);
    auto *t_1448 = buffer.data(target + 1448);
    auto *t_1449 = buffer.data(target + 1449);
    auto *t_1450 = buffer.data(target + 1450);
    auto *t_1451 = buffer.data(target + 1451);
    auto *t_1452 = buffer.data(target + 1452);
    auto *t_1453 = buffer.data(target + 1453);
    auto *t_1454 = buffer.data(target + 1454);
    auto *t_1455 = buffer.data(target + 1455);
    auto *t_1456 = buffer.data(target + 1456);
    auto *t_1457 = buffer.data(target + 1457);
    auto *t_1458 = buffer.data(target + 1458);
    auto *t_1459 = buffer.data(target + 1459);
    auto *t_1460 = buffer.data(target + 1460);
    auto *t_1461 = buffer.data(target + 1461);
    auto *t_1462 = buffer.data(target + 1462);
    auto *t_1463 = buffer.data(target + 1463);
    auto *t_1464 = buffer.data(target + 1464);
    auto *t_1465 = buffer.data(target + 1465);
    auto *t_1466 = buffer.data(target + 1466);
    auto *t_1467 = buffer.data(target + 1467);
    auto *t_1468 = buffer.data(target + 1468);
    auto *t_1469 = buffer.data(target + 1469);
    auto *t_1470 = buffer.data(target + 1470);
    auto *t_1471 = buffer.data(target + 1471);
    auto *t_1472 = buffer.data(target + 1472);
    auto *t_1473 = buffer.data(target + 1473);
    auto *t_1474 = buffer.data(target + 1474);
    auto *t_1475 = buffer.data(target + 1475);
    auto *t_1476 = buffer.data(target + 1476);
    auto *t_1477 = buffer.data(target + 1477);
    auto *t_1478 = buffer.data(target + 1478);
    auto *t_1479 = buffer.data(target + 1479);
    auto *t_1480 = buffer.data(target + 1480);
    auto *t_1481 = buffer.data(target + 1481);
    auto *t_1482 = buffer.data(target + 1482);
    auto *t_1483 = buffer.data(target + 1483);
    auto *t_1484 = buffer.data(target + 1484);
    auto *t_1485 = buffer.data(target + 1485);
    auto *t_1486 = buffer.data(target + 1486);
    auto *t_1487 = buffer.data(target + 1487);
    auto *t_1488 = buffer.data(target + 1488);
    auto *t_1489 = buffer.data(target + 1489);
    auto *t_1490 = buffer.data(target + 1490);
    auto *t_1491 = buffer.data(target + 1491);
    auto *t_1492 = buffer.data(target + 1492);
    auto *t_1493 = buffer.data(target + 1493);
    auto *t_1494 = buffer.data(target + 1494);
    auto *t_1495 = buffer.data(target + 1495);
    auto *t_1496 = buffer.data(target + 1496);
    auto *t_1497 = buffer.data(target + 1497);
    auto *t_1498 = buffer.data(target + 1498);
    auto *t_1499 = buffer.data(target + 1499);
    auto *t_1500 = buffer.data(target + 1500);
    auto *t_1501 = buffer.data(target + 1501);
    auto *t_1502 = buffer.data(target + 1502);
    auto *t_1503 = buffer.data(target + 1503);
    auto *t_1504 = buffer.data(target + 1504);
    auto *t_1505 = buffer.data(target + 1505);
    auto *t_1506 = buffer.data(target + 1506);
    auto *t_1507 = buffer.data(target + 1507);
    auto *t_1508 = buffer.data(target + 1508);
    auto *t_1509 = buffer.data(target + 1509);
    auto *t_1510 = buffer.data(target + 1510);
    auto *t_1511 = buffer.data(target + 1511);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli0_1232 = buffer.data(sli0 + 1232);
    const auto *sli0_1237 = buffer.data(sli0 + 1237);
    const auto *sli0_1241 = buffer.data(sli0 + 1241);
    const auto *sli0_1246 = buffer.data(sli0 + 1246);
    const auto *sli0_1253 = buffer.data(sli0 + 1253);
    const auto *sli0_1255 = buffer.data(sli0 + 1255);
    const auto *sli0_1256 = buffer.data(sli0 + 1256);
    const auto *sli0_1257 = buffer.data(sli0 + 1257);
    const auto *sli0_1259 = buffer.data(sli0 + 1259);

    const auto *slh_834 = buffer.data(slh + 834);
    const auto *slh_839 = buffer.data(slh + 839);
    const auto *slh_840 = buffer.data(slh + 840);
    const auto *slh_843 = buffer.data(slh + 843);
    const auto *slh_846 = buffer.data(slh + 846);
    const auto *slh_855 = buffer.data(slh + 855);
    const auto *slh_857 = buffer.data(slh + 857);
    const auto *slh_858 = buffer.data(slh + 858);
    const auto *slh_859 = buffer.data(slh + 859);
    const auto *slh_860 = buffer.data(slh + 860);
    const auto *slh_861 = buffer.data(slh + 861);
    const auto *slh_863 = buffer.data(slh + 863);
    const auto *slh_864 = buffer.data(slh + 864);
    const auto *slh_866 = buffer.data(slh + 866);
    const auto *slh_867 = buffer.data(slh + 867);
    const auto *slh_870 = buffer.data(slh + 870);
    const auto *slh_876 = buffer.data(slh + 876);
    const auto *slh_878 = buffer.data(slh + 878);
    const auto *slh_879 = buffer.data(slh + 879);
    const auto *slh_880 = buffer.data(slh + 880);
    const auto *slh_881 = buffer.data(slh + 881);
    const auto *slh_882 = buffer.data(slh + 882);
    const auto *slh_884 = buffer.data(slh + 884);
    const auto *slh_885 = buffer.data(slh + 885);
    const auto *slh_887 = buffer.data(slh + 887);
    const auto *slh_888 = buffer.data(slh + 888);
    const auto *slh_891 = buffer.data(slh + 891);
    const auto *slh_897 = buffer.data(slh + 897);
    const auto *slh_899 = buffer.data(slh + 899);
    const auto *slh_900 = buffer.data(slh + 900);
    const auto *slh_901 = buffer.data(slh + 901);
    const auto *slh_902 = buffer.data(slh + 902);
    const auto *slh_903 = buffer.data(slh + 903);
    const auto *slh_905 = buffer.data(slh + 905);
    const auto *slh_906 = buffer.data(slh + 906);
    const auto *slh_908 = buffer.data(slh + 908);
    const auto *slh_909 = buffer.data(slh + 909);
    const auto *slh_912 = buffer.data(slh + 912);
    const auto *slh_918 = buffer.data(slh + 918);
    const auto *slh_920 = buffer.data(slh + 920);
    const auto *slh_921 = buffer.data(slh + 921);
    const auto *slh_922 = buffer.data(slh + 922);
    const auto *slh_923 = buffer.data(slh + 923);
    const auto *slh_924 = buffer.data(slh + 924);
    const auto *slh_926 = buffer.data(slh + 926);
    const auto *slh_929 = buffer.data(slh + 929);
    const auto *slh_933 = buffer.data(slh + 933);
    const auto *slh_939 = buffer.data(slh + 939);
    const auto *slh_941 = buffer.data(slh + 941);
    const auto *slh_942 = buffer.data(slh + 942);
    const auto *slh_943 = buffer.data(slh + 943);
    const auto *slh_944 = buffer.data(slh + 944);

    const auto *sli1_1232 = buffer.data(sli1 + 1232);
    const auto *sli1_1237 = buffer.data(sli1 + 1237);
    const auto *sli1_1241 = buffer.data(sli1 + 1241);
    const auto *sli1_1246 = buffer.data(sli1 + 1246);
    const auto *sli1_1253 = buffer.data(sli1 + 1253);
    const auto *sli1_1255 = buffer.data(sli1 + 1255);
    const auto *sli1_1256 = buffer.data(sli1 + 1256);
    const auto *sli1_1257 = buffer.data(sli1 + 1257);
    const auto *sli1_1259 = buffer.data(sli1 + 1259);

    const auto *smg0_745 = buffer.data(smg0 + 745);
    const auto *smg0_747 = buffer.data(smg0 + 747);
    const auto *smg0_748 = buffer.data(smg0 + 748);
    const auto *smg0_749 = buffer.data(smg0 + 749);
    const auto *smg0_750 = buffer.data(smg0 + 750);
    const auto *smg0_753 = buffer.data(smg0 + 753);
    const auto *smg0_755 = buffer.data(smg0 + 755);
    const auto *smg0_756 = buffer.data(smg0 + 756);
    const auto *smg0_759 = buffer.data(smg0 + 759);
    const auto *smg0_760 = buffer.data(smg0 + 760);
    const auto *smg0_762 = buffer.data(smg0 + 762);
    const auto *smg0_763 = buffer.data(smg0 + 763);
    const auto *smg0_764 = buffer.data(smg0 + 764);
    const auto *smg0_765 = buffer.data(smg0 + 765);
    const auto *smg0_768 = buffer.data(smg0 + 768);
    const auto *smg0_770 = buffer.data(smg0 + 770);
    const auto *smg0_771 = buffer.data(smg0 + 771);
    const auto *smg0_774 = buffer.data(smg0 + 774);
    const auto *smg0_775 = buffer.data(smg0 + 775);
    const auto *smg0_777 = buffer.data(smg0 + 777);
    const auto *smg0_778 = buffer.data(smg0 + 778);
    const auto *smg0_779 = buffer.data(smg0 + 779);
    const auto *smg0_780 = buffer.data(smg0 + 780);
    const auto *smg0_783 = buffer.data(smg0 + 783);
    const auto *smg0_785 = buffer.data(smg0 + 785);
    const auto *smg0_786 = buffer.data(smg0 + 786);
    const auto *smg0_789 = buffer.data(smg0 + 789);
    const auto *smg0_790 = buffer.data(smg0 + 790);
    const auto *smg0_792 = buffer.data(smg0 + 792);
    const auto *smg0_793 = buffer.data(smg0 + 793);
    const auto *smg0_794 = buffer.data(smg0 + 794);
    const auto *smg0_798 = buffer.data(smg0 + 798);
    const auto *smg0_801 = buffer.data(smg0 + 801);
    const auto *smg0_805 = buffer.data(smg0 + 805);
    const auto *smg0_807 = buffer.data(smg0 + 807);

    const auto *smg1_745 = buffer.data(smg1 + 745);
    const auto *smg1_747 = buffer.data(smg1 + 747);
    const auto *smg1_748 = buffer.data(smg1 + 748);
    const auto *smg1_749 = buffer.data(smg1 + 749);
    const auto *smg1_750 = buffer.data(smg1 + 750);
    const auto *smg1_753 = buffer.data(smg1 + 753);
    const auto *smg1_755 = buffer.data(smg1 + 755);
    const auto *smg1_756 = buffer.data(smg1 + 756);
    const auto *smg1_759 = buffer.data(smg1 + 759);
    const auto *smg1_760 = buffer.data(smg1 + 760);
    const auto *smg1_762 = buffer.data(smg1 + 762);
    const auto *smg1_763 = buffer.data(smg1 + 763);
    const auto *smg1_764 = buffer.data(smg1 + 764);
    const auto *smg1_765 = buffer.data(smg1 + 765);
    const auto *smg1_768 = buffer.data(smg1 + 768);
    const auto *smg1_770 = buffer.data(smg1 + 770);
    const auto *smg1_771 = buffer.data(smg1 + 771);
    const auto *smg1_774 = buffer.data(smg1 + 774);
    const auto *smg1_775 = buffer.data(smg1 + 775);
    const auto *smg1_777 = buffer.data(smg1 + 777);
    const auto *smg1_778 = buffer.data(smg1 + 778);
    const auto *smg1_779 = buffer.data(smg1 + 779);
    const auto *smg1_780 = buffer.data(smg1 + 780);
    const auto *smg1_783 = buffer.data(smg1 + 783);
    const auto *smg1_785 = buffer.data(smg1 + 785);
    const auto *smg1_786 = buffer.data(smg1 + 786);
    const auto *smg1_789 = buffer.data(smg1 + 789);
    const auto *smg1_790 = buffer.data(smg1 + 790);
    const auto *smg1_792 = buffer.data(smg1 + 792);
    const auto *smg1_793 = buffer.data(smg1 + 793);
    const auto *smg1_794 = buffer.data(smg1 + 794);
    const auto *smg1_798 = buffer.data(smg1 + 798);
    const auto *smg1_801 = buffer.data(smg1 + 801);
    const auto *smg1_805 = buffer.data(smg1 + 805);
    const auto *smg1_807 = buffer.data(smg1 + 807);

    const auto *smh_1044 = buffer.data(smh + 1044);
    const auto *smh_1045 = buffer.data(smh + 1045);
    const auto *smh_1046 = buffer.data(smh + 1046);
    const auto *smh_1047 = buffer.data(smh + 1047);
    const auto *smh_1048 = buffer.data(smh + 1048);
    const auto *smh_1049 = buffer.data(smh + 1049);
    const auto *smh_1050 = buffer.data(smh + 1050);
    const auto *smh_1052 = buffer.data(smh + 1052);
    const auto *smh_1053 = buffer.data(smh + 1053);
    const auto *smh_1055 = buffer.data(smh + 1055);
    const auto *smh_1056 = buffer.data(smh + 1056);
    const auto *smh_1059 = buffer.data(smh + 1059);
    const auto *smh_1060 = buffer.data(smh + 1060);
    const auto *smh_1062 = buffer.data(smh + 1062);
    const auto *smh_1064 = buffer.data(smh + 1064);
    const auto *smh_1065 = buffer.data(smh + 1065);
    const auto *smh_1066 = buffer.data(smh + 1066);
    const auto *smh_1067 = buffer.data(smh + 1067);
    const auto *smh_1068 = buffer.data(smh + 1068);
    const auto *smh_1069 = buffer.data(smh + 1069);
    const auto *smh_1070 = buffer.data(smh + 1070);
    const auto *smh_1071 = buffer.data(smh + 1071);
    const auto *smh_1073 = buffer.data(smh + 1073);
    const auto *smh_1074 = buffer.data(smh + 1074);
    const auto *smh_1076 = buffer.data(smh + 1076);
    const auto *smh_1077 = buffer.data(smh + 1077);
    const auto *smh_1080 = buffer.data(smh + 1080);
    const auto *smh_1081 = buffer.data(smh + 1081);
    const auto *smh_1083 = buffer.data(smh + 1083);
    const auto *smh_1085 = buffer.data(smh + 1085);
    const auto *smh_1086 = buffer.data(smh + 1086);
    const auto *smh_1087 = buffer.data(smh + 1087);
    const auto *smh_1088 = buffer.data(smh + 1088);
    const auto *smh_1089 = buffer.data(smh + 1089);
    const auto *smh_1090 = buffer.data(smh + 1090);
    const auto *smh_1091 = buffer.data(smh + 1091);
    const auto *smh_1092 = buffer.data(smh + 1092);
    const auto *smh_1094 = buffer.data(smh + 1094);
    const auto *smh_1095 = buffer.data(smh + 1095);
    const auto *smh_1097 = buffer.data(smh + 1097);
    const auto *smh_1098 = buffer.data(smh + 1098);
    const auto *smh_1101 = buffer.data(smh + 1101);
    const auto *smh_1102 = buffer.data(smh + 1102);
    const auto *smh_1104 = buffer.data(smh + 1104);
    const auto *smh_1106 = buffer.data(smh + 1106);
    const auto *smh_1107 = buffer.data(smh + 1107);
    const auto *smh_1108 = buffer.data(smh + 1108);
    const auto *smh_1109 = buffer.data(smh + 1109);
    const auto *smh_1110 = buffer.data(smh + 1110);
    const auto *smh_1111 = buffer.data(smh + 1111);
    const auto *smh_1112 = buffer.data(smh + 1112);
    const auto *smh_1113 = buffer.data(smh + 1113);
    const auto *smh_1115 = buffer.data(smh + 1115);
    const auto *smh_1116 = buffer.data(smh + 1116);
    const auto *smh_1118 = buffer.data(smh + 1118);
    const auto *smh_1119 = buffer.data(smh + 1119);
    const auto *smh_1122 = buffer.data(smh + 1122);
    const auto *smh_1123 = buffer.data(smh + 1123);
    const auto *smh_1125 = buffer.data(smh + 1125);
    const auto *smh_1128 = buffer.data(smh + 1128);
    const auto *smh_1129 = buffer.data(smh + 1129);
    const auto *smh_1130 = buffer.data(smh + 1130);
    const auto *smh_1131 = buffer.data(smh + 1131);
    const auto *smh_1132 = buffer.data(smh + 1132);
    const auto *smh_1133 = buffer.data(smh + 1133);

#pragma omp simd aligned(t_1388, t_1389, t_1390, t_1391, t_1392, pc_x, smh_1045, smh_1046, \
                         smh_1047, smh_1048, smh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = f_3 * pc_x[k] * smh_1045[k];

        t_1389[k] = f_3 * pc_x[k] * smh_1046[k];

        t_1390[k] = f_3 * pc_x[k] * smh_1047[k];

        t_1391[k] = f_3 * pc_x[k] * smh_1048[k];

        t_1392[k] = f_3 * pc_x[k] * smh_1049[k];
    }

#pragma omp simd aligned(t_1393, t_1394, t_1395, pc_y, pc_z, slh_834, slh_855, slh_857, \
                         smg0_745, smg0_747, smg1_745, smg1_747, smh_1044, \
                         smh_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1393[k] = f_18 * slh_855[k]
                    + f_1 * smg0_745[k]
                    - f_2 * smg1_745[k]
                    + f_3 * pc_y[k] * smh_1044[k];

        t_1394[k] = f_14 * slh_834[k]
                    + f_3 * pc_z[k] * smh_1044[k];

        t_1395[k] = f_18 * slh_857[k]
                    + f_4 * smg0_747[k]
                    - f_5 * smg1_747[k]
                    + f_3 * pc_y[k] * smh_1046[k];
    }

#pragma omp simd aligned(t_1396, t_1397, t_1398, pc_y, slh_858, slh_859, slh_860, smg0_748, \
                         smg0_749, smg1_748, smg1_749, smh_1047, smh_1048, \
                         smh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1396[k] = f_18 * slh_858[k]
                    + f_6 * smg0_748[k]
                    - f_7 * smg1_748[k]
                    + f_3 * pc_y[k] * smh_1047[k];

        t_1397[k] = f_18 * slh_859[k]
                    + f_8 * smg0_749[k]
                    - f_9 * smg1_749[k]
                    + f_3 * pc_y[k] * smh_1048[k];

        t_1398[k] = f_18 * slh_860[k]
                    + f_3 * pc_y[k] * smh_1049[k];
    }

#pragma omp simd aligned(t_1399, t_1400, t_1401, t_1402, pc_x, pc_y, pc_z, slh_839, slh_840, \
                         slh_861, smg0_749, smg0_750, smg1_749, smg1_750, smh_1049, \
                         smh_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1399[k] = f_14 * slh_839[k]
                    + f_1 * smg0_749[k]
                    - f_2 * smg1_749[k]
                    + f_3 * pc_z[k] * smh_1049[k];

        t_1400[k] = f_1 * smg0_750[k]
                    - f_2 * smg1_750[k]
                    + f_3 * pc_x[k] * smh_1050[k];

        t_1401[k] = f_14 * slh_861[k]
                    + f_3 * pc_y[k] * smh_1050[k];

        t_1402[k] = f_18 * slh_840[k]
                    + f_3 * pc_z[k] * smh_1050[k];
    }

#pragma omp simd aligned(t_1403, t_1404, t_1405, pc_x, pc_y, slh_863, smg0_753, smg0_755, \
                         smg1_753, smg1_755, smh_1052, smh_1053, \
                         smh_1055 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1403[k] = f_4 * smg0_753[k]
                    - f_5 * smg1_753[k]
                    + f_3 * pc_x[k] * smh_1053[k];

        t_1404[k] = f_14 * slh_863[k]
                    + f_3 * pc_y[k] * smh_1052[k];

        t_1405[k] = f_4 * smg0_755[k]
                    - f_5 * smg1_755[k]
                    + f_3 * pc_x[k] * smh_1055[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, pc_x, pc_y, pc_z, slh_843, slh_866, smg0_756, \
                         smg1_756, smh_1053, smh_1055, smh_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_6 * smg0_756[k]
                    - f_7 * smg1_756[k]
                    + f_3 * pc_x[k] * smh_1056[k];

        t_1407[k] = f_18 * slh_843[k]
                    + f_3 * pc_z[k] * smh_1053[k];

        t_1408[k] = f_14 * slh_866[k]
                    + f_3 * pc_y[k] * smh_1055[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pc_x, pc_z, slh_846, smg0_759, smg0_760, \
                         smg1_759, smg1_760, smh_1056, smh_1059, \
                         smh_1060 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_6 * smg0_759[k]
                    - f_7 * smg1_759[k]
                    + f_3 * pc_x[k] * smh_1059[k];

        t_1410[k] = f_8 * smg0_760[k]
                    - f_9 * smg1_760[k]
                    + f_3 * pc_x[k] * smh_1060[k];

        t_1411[k] = f_18 * slh_846[k]
                    + f_3 * pc_z[k] * smh_1056[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, t_1415, pc_x, pc_y, slh_870, smg0_762, \
                         smg0_764, smg1_762, smg1_764, smh_1059, smh_1062, smh_1064, \
                         smh_1065 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_8 * smg0_762[k]
                    - f_9 * smg1_762[k]
                    + f_3 * pc_x[k] * smh_1062[k];

        t_1413[k] = f_14 * slh_870[k]
                    + f_3 * pc_y[k] * smh_1059[k];

        t_1414[k] = f_8 * smg0_764[k]
                    - f_9 * smg1_764[k]
                    + f_3 * pc_x[k] * smh_1064[k];

        t_1415[k] = f_3 * pc_x[k] * smh_1065[k];
    }

#pragma omp simd aligned(t_1416, t_1417, t_1418, t_1419, t_1420, pc_x, smh_1066, smh_1067, \
                         smh_1068, smh_1069, smh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1416[k] = f_3 * pc_x[k] * smh_1066[k];

        t_1417[k] = f_3 * pc_x[k] * smh_1067[k];

        t_1418[k] = f_3 * pc_x[k] * smh_1068[k];

        t_1419[k] = f_3 * pc_x[k] * smh_1069[k];

        t_1420[k] = f_3 * pc_x[k] * smh_1070[k];
    }

#pragma omp simd aligned(t_1421, t_1422, t_1423, pc_y, pc_z, slh_855, slh_876, slh_878, \
                         smg0_760, smg0_762, smg1_760, smg1_762, smh_1065, \
                         smh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1421[k] = f_14 * slh_876[k]
                    + f_1 * smg0_760[k]
                    - f_2 * smg1_760[k]
                    + f_3 * pc_y[k] * smh_1065[k];

        t_1422[k] = f_18 * slh_855[k]
                    + f_3 * pc_z[k] * smh_1065[k];

        t_1423[k] = f_14 * slh_878[k]
                    + f_4 * smg0_762[k]
                    - f_5 * smg1_762[k]
                    + f_3 * pc_y[k] * smh_1067[k];
    }

#pragma omp simd aligned(t_1424, t_1425, t_1426, pc_y, slh_879, slh_880, slh_881, smg0_763, \
                         smg0_764, smg1_763, smg1_764, smh_1068, smh_1069, \
                         smh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1424[k] = f_14 * slh_879[k]
                    + f_6 * smg0_763[k]
                    - f_7 * smg1_763[k]
                    + f_3 * pc_y[k] * smh_1068[k];

        t_1425[k] = f_14 * slh_880[k]
                    + f_8 * smg0_764[k]
                    - f_9 * smg1_764[k]
                    + f_3 * pc_y[k] * smh_1069[k];

        t_1426[k] = f_14 * slh_881[k]
                    + f_3 * pc_y[k] * smh_1070[k];
    }

#pragma omp simd aligned(t_1427, t_1428, t_1429, t_1430, pc_x, pc_y, pc_z, slh_860, slh_861, \
                         slh_882, smg0_764, smg0_765, smg1_764, smg1_765, smh_1070, \
                         smh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1427[k] = f_18 * slh_860[k]
                    + f_1 * smg0_764[k]
                    - f_2 * smg1_764[k]
                    + f_3 * pc_z[k] * smh_1070[k];

        t_1428[k] = f_1 * smg0_765[k]
                    - f_2 * smg1_765[k]
                    + f_3 * pc_x[k] * smh_1071[k];

        t_1429[k] = f_13 * slh_882[k]
                    + f_3 * pc_y[k] * smh_1071[k];

        t_1430[k] = f_17 * slh_861[k]
                    + f_3 * pc_z[k] * smh_1071[k];
    }

#pragma omp simd aligned(t_1431, t_1432, t_1433, pc_x, pc_y, slh_884, smg0_768, smg0_770, \
                         smg1_768, smg1_770, smh_1073, smh_1074, \
                         smh_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1431[k] = f_4 * smg0_768[k]
                    - f_5 * smg1_768[k]
                    + f_3 * pc_x[k] * smh_1074[k];

        t_1432[k] = f_13 * slh_884[k]
                    + f_3 * pc_y[k] * smh_1073[k];

        t_1433[k] = f_4 * smg0_770[k]
                    - f_5 * smg1_770[k]
                    + f_3 * pc_x[k] * smh_1076[k];
    }

#pragma omp simd aligned(t_1434, t_1435, t_1436, pc_x, pc_y, pc_z, slh_864, slh_887, smg0_771, \
                         smg1_771, smh_1074, smh_1076, smh_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1434[k] = f_6 * smg0_771[k]
                    - f_7 * smg1_771[k]
                    + f_3 * pc_x[k] * smh_1077[k];

        t_1435[k] = f_17 * slh_864[k]
                    + f_3 * pc_z[k] * smh_1074[k];

        t_1436[k] = f_13 * slh_887[k]
                    + f_3 * pc_y[k] * smh_1076[k];
    }

#pragma omp simd aligned(t_1437, t_1438, t_1439, pc_x, pc_z, slh_867, smg0_774, smg0_775, \
                         smg1_774, smg1_775, smh_1077, smh_1080, \
                         smh_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1437[k] = f_6 * smg0_774[k]
                    - f_7 * smg1_774[k]
                    + f_3 * pc_x[k] * smh_1080[k];

        t_1438[k] = f_8 * smg0_775[k]
                    - f_9 * smg1_775[k]
                    + f_3 * pc_x[k] * smh_1081[k];

        t_1439[k] = f_17 * slh_867[k]
                    + f_3 * pc_z[k] * smh_1077[k];
    }

#pragma omp simd aligned(t_1440, t_1441, t_1442, t_1443, pc_x, pc_y, slh_891, smg0_777, \
                         smg0_779, smg1_777, smg1_779, smh_1080, smh_1083, smh_1085, \
                         smh_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1440[k] = f_8 * smg0_777[k]
                    - f_9 * smg1_777[k]
                    + f_3 * pc_x[k] * smh_1083[k];

        t_1441[k] = f_13 * slh_891[k]
                    + f_3 * pc_y[k] * smh_1080[k];

        t_1442[k] = f_8 * smg0_779[k]
                    - f_9 * smg1_779[k]
                    + f_3 * pc_x[k] * smh_1085[k];

        t_1443[k] = f_3 * pc_x[k] * smh_1086[k];
    }

#pragma omp simd aligned(t_1444, t_1445, t_1446, t_1447, t_1448, pc_x, smh_1087, smh_1088, \
                         smh_1089, smh_1090, smh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1444[k] = f_3 * pc_x[k] * smh_1087[k];

        t_1445[k] = f_3 * pc_x[k] * smh_1088[k];

        t_1446[k] = f_3 * pc_x[k] * smh_1089[k];

        t_1447[k] = f_3 * pc_x[k] * smh_1090[k];

        t_1448[k] = f_3 * pc_x[k] * smh_1091[k];
    }

#pragma omp simd aligned(t_1449, t_1450, t_1451, pc_y, pc_z, slh_876, slh_897, slh_899, \
                         smg0_775, smg0_777, smg1_775, smg1_777, smh_1086, \
                         smh_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1449[k] = f_13 * slh_897[k]
                    + f_1 * smg0_775[k]
                    - f_2 * smg1_775[k]
                    + f_3 * pc_y[k] * smh_1086[k];

        t_1450[k] = f_17 * slh_876[k]
                    + f_3 * pc_z[k] * smh_1086[k];

        t_1451[k] = f_13 * slh_899[k]
                    + f_4 * smg0_777[k]
                    - f_5 * smg1_777[k]
                    + f_3 * pc_y[k] * smh_1088[k];
    }

#pragma omp simd aligned(t_1452, t_1453, t_1454, pc_y, slh_900, slh_901, slh_902, smg0_778, \
                         smg0_779, smg1_778, smg1_779, smh_1089, smh_1090, \
                         smh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1452[k] = f_13 * slh_900[k]
                    + f_6 * smg0_778[k]
                    - f_7 * smg1_778[k]
                    + f_3 * pc_y[k] * smh_1089[k];

        t_1453[k] = f_13 * slh_901[k]
                    + f_8 * smg0_779[k]
                    - f_9 * smg1_779[k]
                    + f_3 * pc_y[k] * smh_1090[k];

        t_1454[k] = f_13 * slh_902[k]
                    + f_3 * pc_y[k] * smh_1091[k];
    }

#pragma omp simd aligned(t_1455, t_1456, t_1457, t_1458, pc_x, pc_y, pc_z, slh_881, slh_882, \
                         slh_903, smg0_779, smg0_780, smg1_779, smg1_780, smh_1091, \
                         smh_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1455[k] = f_17 * slh_881[k]
                    + f_1 * smg0_779[k]
                    - f_2 * smg1_779[k]
                    + f_3 * pc_z[k] * smh_1091[k];

        t_1456[k] = f_1 * smg0_780[k]
                    - f_2 * smg1_780[k]
                    + f_3 * pc_x[k] * smh_1092[k];

        t_1457[k] = f_12 * slh_903[k]
                    + f_3 * pc_y[k] * smh_1092[k];

        t_1458[k] = f_16 * slh_882[k]
                    + f_3 * pc_z[k] * smh_1092[k];
    }

#pragma omp simd aligned(t_1459, t_1460, t_1461, pc_x, pc_y, slh_905, smg0_783, smg0_785, \
                         smg1_783, smg1_785, smh_1094, smh_1095, \
                         smh_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1459[k] = f_4 * smg0_783[k]
                    - f_5 * smg1_783[k]
                    + f_3 * pc_x[k] * smh_1095[k];

        t_1460[k] = f_12 * slh_905[k]
                    + f_3 * pc_y[k] * smh_1094[k];

        t_1461[k] = f_4 * smg0_785[k]
                    - f_5 * smg1_785[k]
                    + f_3 * pc_x[k] * smh_1097[k];
    }

#pragma omp simd aligned(t_1462, t_1463, t_1464, pc_x, pc_y, pc_z, slh_885, slh_908, smg0_786, \
                         smg1_786, smh_1095, smh_1097, smh_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1462[k] = f_6 * smg0_786[k]
                    - f_7 * smg1_786[k]
                    + f_3 * pc_x[k] * smh_1098[k];

        t_1463[k] = f_16 * slh_885[k]
                    + f_3 * pc_z[k] * smh_1095[k];

        t_1464[k] = f_12 * slh_908[k]
                    + f_3 * pc_y[k] * smh_1097[k];
    }

#pragma omp simd aligned(t_1465, t_1466, t_1467, pc_x, pc_z, slh_888, smg0_789, smg0_790, \
                         smg1_789, smg1_790, smh_1098, smh_1101, \
                         smh_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1465[k] = f_6 * smg0_789[k]
                    - f_7 * smg1_789[k]
                    + f_3 * pc_x[k] * smh_1101[k];

        t_1466[k] = f_8 * smg0_790[k]
                    - f_9 * smg1_790[k]
                    + f_3 * pc_x[k] * smh_1102[k];

        t_1467[k] = f_16 * slh_888[k]
                    + f_3 * pc_z[k] * smh_1098[k];
    }

#pragma omp simd aligned(t_1468, t_1469, t_1470, t_1471, pc_x, pc_y, slh_912, smg0_792, \
                         smg0_794, smg1_792, smg1_794, smh_1101, smh_1104, smh_1106, \
                         smh_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1468[k] = f_8 * smg0_792[k]
                    - f_9 * smg1_792[k]
                    + f_3 * pc_x[k] * smh_1104[k];

        t_1469[k] = f_12 * slh_912[k]
                    + f_3 * pc_y[k] * smh_1101[k];

        t_1470[k] = f_8 * smg0_794[k]
                    - f_9 * smg1_794[k]
                    + f_3 * pc_x[k] * smh_1106[k];

        t_1471[k] = f_3 * pc_x[k] * smh_1107[k];
    }

#pragma omp simd aligned(t_1472, t_1473, t_1474, t_1475, t_1476, pc_x, smh_1108, smh_1109, \
                         smh_1110, smh_1111, smh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1472[k] = f_3 * pc_x[k] * smh_1108[k];

        t_1473[k] = f_3 * pc_x[k] * smh_1109[k];

        t_1474[k] = f_3 * pc_x[k] * smh_1110[k];

        t_1475[k] = f_3 * pc_x[k] * smh_1111[k];

        t_1476[k] = f_3 * pc_x[k] * smh_1112[k];
    }

#pragma omp simd aligned(t_1477, t_1478, t_1479, pc_y, pc_z, slh_897, slh_918, slh_920, \
                         smg0_790, smg0_792, smg1_790, smg1_792, smh_1107, \
                         smh_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1477[k] = f_12 * slh_918[k]
                    + f_1 * smg0_790[k]
                    - f_2 * smg1_790[k]
                    + f_3 * pc_y[k] * smh_1107[k];

        t_1478[k] = f_16 * slh_897[k]
                    + f_3 * pc_z[k] * smh_1107[k];

        t_1479[k] = f_12 * slh_920[k]
                    + f_4 * smg0_792[k]
                    - f_5 * smg1_792[k]
                    + f_3 * pc_y[k] * smh_1109[k];
    }

#pragma omp simd aligned(t_1480, t_1481, t_1482, pc_y, slh_921, slh_922, slh_923, smg0_793, \
                         smg0_794, smg1_793, smg1_794, smh_1110, smh_1111, \
                         smh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1480[k] = f_12 * slh_921[k]
                    + f_6 * smg0_793[k]
                    - f_7 * smg1_793[k]
                    + f_3 * pc_y[k] * smh_1110[k];

        t_1481[k] = f_12 * slh_922[k]
                    + f_8 * smg0_794[k]
                    - f_9 * smg1_794[k]
                    + f_3 * pc_y[k] * smh_1111[k];

        t_1482[k] = f_12 * slh_923[k]
                    + f_3 * pc_y[k] * smh_1112[k];
    }

#pragma omp simd aligned(t_1483, t_1484, t_1485, t_1486, pb_y, pc_y, pc_z, sli0_1232, slh_902, \
                         slh_903, slh_924, sli1_1232, smg0_794, smg1_794, smh_1112, \
                         smh_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1483[k] = f_16 * slh_902[k]
                    + f_1 * smg0_794[k]
                    - f_2 * smg1_794[k]
                    + f_3 * pc_z[k] * smh_1112[k];

        t_1484[k] = pb_y[k] * sli0_1232[k]
                    - f_10 * pc_y[k] * sli1_1232[k];

        t_1485[k] = f_11 * slh_924[k]
                    + f_3 * pc_y[k] * smh_1113[k];

        t_1486[k] = f_15 * slh_903[k]
                    + f_3 * pc_z[k] * smh_1113[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pb_y, pc_x, pc_y, sli0_1237, slh_926, \
                         sli1_1237, smg0_798, smg1_798, smh_1115, \
                         smh_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_4 * smg0_798[k]
                    - f_5 * smg1_798[k]
                    + f_3 * pc_x[k] * smh_1116[k];

        t_1488[k] = f_11 * slh_926[k]
                    + f_3 * pc_y[k] * smh_1115[k];

        t_1489[k] = pb_y[k] * sli0_1237[k]
                    - f_10 * pc_y[k] * sli1_1237[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pc_x, pc_y, pc_z, slh_906, slh_929, smg0_801, \
                         smg1_801, smh_1116, smh_1118, smh_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = f_6 * smg0_801[k]
                    - f_7 * smg1_801[k]
                    + f_3 * pc_x[k] * smh_1119[k];

        t_1491[k] = f_15 * slh_906[k]
                    + f_3 * pc_z[k] * smh_1116[k];

        t_1492[k] = f_11 * slh_929[k]
                    + f_3 * pc_y[k] * smh_1118[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pb_y, pc_x, pc_y, pc_z, sli0_1241, slh_909, \
                         sli1_1241, smg0_805, smg1_805, smh_1119, \
                         smh_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = pb_y[k] * sli0_1241[k]
                    - f_10 * pc_y[k] * sli1_1241[k];

        t_1494[k] = f_8 * smg0_805[k]
                    - f_9 * smg1_805[k]
                    + f_3 * pc_x[k] * smh_1123[k];

        t_1495[k] = f_15 * slh_909[k]
                    + f_3 * pc_z[k] * smh_1119[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, t_1499, pb_y, pc_x, pc_y, sli0_1246, slh_933, \
                         sli1_1246, smg0_807, smg1_807, smh_1122, smh_1125, \
                         smh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = f_8 * smg0_807[k]
                    - f_9 * smg1_807[k]
                    + f_3 * pc_x[k] * smh_1125[k];

        t_1497[k] = f_11 * slh_933[k]
                    + f_3 * pc_y[k] * smh_1122[k];

        t_1498[k] = pb_y[k] * sli0_1246[k]
                    - f_10 * pc_y[k] * sli1_1246[k];

        t_1499[k] = f_3 * pc_x[k] * smh_1128[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, t_1503, t_1504, pc_x, smh_1129, smh_1130, \
                         smh_1131, smh_1132, smh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_3 * pc_x[k] * smh_1129[k];

        t_1501[k] = f_3 * pc_x[k] * smh_1130[k];

        t_1502[k] = f_3 * pc_x[k] * smh_1131[k];

        t_1503[k] = f_3 * pc_x[k] * smh_1132[k];

        t_1504[k] = f_3 * pc_x[k] * smh_1133[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, pb_y, pc_y, pc_z, sli0_1253, sli0_1255, \
                         slh_918, slh_939, slh_941, sli1_1253, sli1_1255, \
                         smh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = pb_y[k] * sli0_1253[k]
                    + f_17 * slh_939[k]
                    - f_10 * pc_y[k] * sli1_1253[k];

        t_1506[k] = f_15 * slh_918[k]
                    + f_3 * pc_z[k] * smh_1128[k];

        t_1507[k] = pb_y[k] * sli0_1255[k]
                    + f_14 * slh_941[k]
                    - f_10 * pc_y[k] * sli1_1255[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, t_1511, pb_y, pc_y, sli0_1256, sli0_1257, \
                         sli0_1259, slh_942, slh_943, slh_944, sli1_1256, sli1_1257, \
                         sli1_1259, smh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = pb_y[k] * sli0_1256[k]
                    + f_13 * slh_942[k]
                    - f_10 * pc_y[k] * sli1_1256[k];

        t_1509[k] = pb_y[k] * sli0_1257[k]
                    + f_12 * slh_943[k]
                    - f_10 * pc_y[k] * sli1_1257[k];

        t_1510[k] = f_11 * slh_944[k]
                    + f_3 * pc_y[k] * smh_1133[k];

        t_1511[k] = pb_y[k] * sli0_1259[k]
                    - f_10 * pc_y[k] * sli1_1259[k];
    }
}

static auto
compute_prim_smi_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t slh, const size_t smg0,
                                                           const size_t smg1, const size_t smh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);

    auto *t_1512 = buffer.data(target + 1512);
    auto *t_1513 = buffer.data(target + 1513);
    auto *t_1514 = buffer.data(target + 1514);
    auto *t_1515 = buffer.data(target + 1515);
    auto *t_1516 = buffer.data(target + 1516);
    auto *t_1517 = buffer.data(target + 1517);
    auto *t_1518 = buffer.data(target + 1518);
    auto *t_1519 = buffer.data(target + 1519);
    auto *t_1520 = buffer.data(target + 1520);
    auto *t_1521 = buffer.data(target + 1521);
    auto *t_1522 = buffer.data(target + 1522);
    auto *t_1523 = buffer.data(target + 1523);
    auto *t_1524 = buffer.data(target + 1524);
    auto *t_1525 = buffer.data(target + 1525);
    auto *t_1526 = buffer.data(target + 1526);
    auto *t_1527 = buffer.data(target + 1527);
    auto *t_1528 = buffer.data(target + 1528);
    auto *t_1529 = buffer.data(target + 1529);
    auto *t_1530 = buffer.data(target + 1530);
    auto *t_1531 = buffer.data(target + 1531);
    auto *t_1532 = buffer.data(target + 1532);
    auto *t_1533 = buffer.data(target + 1533);
    auto *t_1534 = buffer.data(target + 1534);
    auto *t_1535 = buffer.data(target + 1535);
    auto *t_1536 = buffer.data(target + 1536);
    auto *t_1537 = buffer.data(target + 1537);
    auto *t_1538 = buffer.data(target + 1538);
    auto *t_1539 = buffer.data(target + 1539);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slh_924 = buffer.data(slh + 924);
    const auto *slh_927 = buffer.data(slh + 927);
    const auto *slh_930 = buffer.data(slh + 930);
    const auto *slh_939 = buffer.data(slh + 939);
    const auto *slh_944 = buffer.data(slh + 944);

    const auto *smg0_810 = buffer.data(smg0 + 810);
    const auto *smg0_813 = buffer.data(smg0 + 813);
    const auto *smg0_815 = buffer.data(smg0 + 815);
    const auto *smg0_816 = buffer.data(smg0 + 816);
    const auto *smg0_819 = buffer.data(smg0 + 819);
    const auto *smg0_820 = buffer.data(smg0 + 820);
    const auto *smg0_822 = buffer.data(smg0 + 822);
    const auto *smg0_823 = buffer.data(smg0 + 823);
    const auto *smg0_824 = buffer.data(smg0 + 824);

    const auto *smg1_810 = buffer.data(smg1 + 810);
    const auto *smg1_813 = buffer.data(smg1 + 813);
    const auto *smg1_815 = buffer.data(smg1 + 815);
    const auto *smg1_816 = buffer.data(smg1 + 816);
    const auto *smg1_819 = buffer.data(smg1 + 819);
    const auto *smg1_820 = buffer.data(smg1 + 820);
    const auto *smg1_822 = buffer.data(smg1 + 822);
    const auto *smg1_823 = buffer.data(smg1 + 823);
    const auto *smg1_824 = buffer.data(smg1 + 824);

    const auto *smh_1134 = buffer.data(smh + 1134);
    const auto *smh_1136 = buffer.data(smh + 1136);
    const auto *smh_1137 = buffer.data(smh + 1137);
    const auto *smh_1139 = buffer.data(smh + 1139);
    const auto *smh_1140 = buffer.data(smh + 1140);
    const auto *smh_1143 = buffer.data(smh + 1143);
    const auto *smh_1144 = buffer.data(smh + 1144);
    const auto *smh_1146 = buffer.data(smh + 1146);
    const auto *smh_1148 = buffer.data(smh + 1148);
    const auto *smh_1149 = buffer.data(smh + 1149);
    const auto *smh_1150 = buffer.data(smh + 1150);
    const auto *smh_1151 = buffer.data(smh + 1151);
    const auto *smh_1152 = buffer.data(smh + 1152);
    const auto *smh_1153 = buffer.data(smh + 1153);
    const auto *smh_1154 = buffer.data(smh + 1154);

#pragma omp simd aligned(t_1512, t_1513, t_1514, t_1515, t_1516, pc_x, pc_y, pc_z, slh_924, \
                         smg0_810, smg0_813, smg1_810, smg1_813, smh_1134, smh_1136, \
                         smh_1137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1512[k] = f_1 * smg0_810[k]
                    - f_2 * smg1_810[k]
                    + f_3 * pc_x[k] * smh_1134[k];

        t_1513[k] = f_3 * pc_y[k] * smh_1134[k];

        t_1514[k] = f_0 * slh_924[k]
                    + f_3 * pc_z[k] * smh_1134[k];

        t_1515[k] = f_4 * smg0_813[k]
                    - f_5 * smg1_813[k]
                    + f_3 * pc_x[k] * smh_1137[k];

        t_1516[k] = f_3 * pc_y[k] * smh_1136[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, t_1520, pc_x, pc_y, pc_z, slh_927, smg0_815, \
                         smg0_816, smg1_815, smg1_816, smh_1137, smh_1139, \
                         smh_1140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = f_4 * smg0_815[k]
                    - f_5 * smg1_815[k]
                    + f_3 * pc_x[k] * smh_1139[k];

        t_1518[k] = f_6 * smg0_816[k]
                    - f_7 * smg1_816[k]
                    + f_3 * pc_x[k] * smh_1140[k];

        t_1519[k] = f_0 * slh_927[k]
                    + f_3 * pc_z[k] * smh_1137[k];

        t_1520[k] = f_3 * pc_y[k] * smh_1139[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, pc_x, pc_z, slh_930, smg0_819, smg0_820, \
                         smg1_819, smg1_820, smh_1140, smh_1143, \
                         smh_1144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = f_6 * smg0_819[k]
                    - f_7 * smg1_819[k]
                    + f_3 * pc_x[k] * smh_1143[k];

        t_1522[k] = f_8 * smg0_820[k]
                    - f_9 * smg1_820[k]
                    + f_3 * pc_x[k] * smh_1144[k];

        t_1523[k] = f_0 * slh_930[k]
                    + f_3 * pc_z[k] * smh_1140[k];
    }

#pragma omp simd aligned(t_1524, t_1525, t_1526, t_1527, t_1528, pc_x, pc_y, smg0_822, \
                         smg0_824, smg1_822, smg1_824, smh_1143, smh_1146, smh_1148, smh_1149, \
                         smh_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1524[k] = f_8 * smg0_822[k]
                    - f_9 * smg1_822[k]
                    + f_3 * pc_x[k] * smh_1146[k];

        t_1525[k] = f_3 * pc_y[k] * smh_1143[k];

        t_1526[k] = f_8 * smg0_824[k]
                    - f_9 * smg1_824[k]
                    + f_3 * pc_x[k] * smh_1148[k];

        t_1527[k] = f_3 * pc_x[k] * smh_1149[k];

        t_1528[k] = f_3 * pc_x[k] * smh_1150[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, t_1532, t_1533, pc_x, pc_y, smg0_820, \
                         smg1_820, smh_1149, smh_1151, smh_1152, smh_1153, \
                         smh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_3 * pc_x[k] * smh_1151[k];

        t_1530[k] = f_3 * pc_x[k] * smh_1152[k];

        t_1531[k] = f_3 * pc_x[k] * smh_1153[k];

        t_1532[k] = f_3 * pc_x[k] * smh_1154[k];

        t_1533[k] = f_1 * smg0_820[k]
                    - f_2 * smg1_820[k]
                    + f_3 * pc_y[k] * smh_1149[k];
    }

#pragma omp simd aligned(t_1534, t_1535, t_1536, pc_y, pc_z, slh_939, smg0_822, smg0_823, \
                         smg1_822, smg1_823, smh_1149, smh_1151, \
                         smh_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1534[k] = f_0 * slh_939[k]
                    + f_3 * pc_z[k] * smh_1149[k];

        t_1535[k] = f_4 * smg0_822[k]
                    - f_5 * smg1_822[k]
                    + f_3 * pc_y[k] * smh_1151[k];

        t_1536[k] = f_6 * smg0_823[k]
                    - f_7 * smg1_823[k]
                    + f_3 * pc_y[k] * smh_1152[k];
    }

#pragma omp simd aligned(t_1537, t_1538, t_1539, pc_y, pc_z, slh_944, smg0_824, smg1_824, \
                         smh_1153, smh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1537[k] = f_8 * smg0_824[k]
                    - f_9 * smg1_824[k]
                    + f_3 * pc_y[k] * smh_1153[k];

        t_1538[k] = f_3 * pc_y[k] * smh_1154[k];

        t_1539[k] = f_0 * slh_944[k]
                    + f_1 * smg0_824[k]
                    - f_2 * smg1_824[k]
                    + f_3 * pc_z[k] * smh_1154[k];
    }
}

auto
compute_prim_smi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sli0, const size_t slh,
                                                   const size_t sli1, const size_t smg0,
                                                   const size_t smg1, const size_t smh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_smi_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smg0, smg1, smh, ncols,
                                                              gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, sli0, slh,
                                                              sli1, smh, ncols, gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, sli0,
                                                               slh, sli1, smg0, smg1, smh,
                                                               ncols, gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, sli0,
                                                               slh, sli1, smg0, smg1, smh,
                                                               ncols, gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece12(buffer, target, pb, pc, sli0,
                                                               slh, sli1, smg0, smg1, smh,
                                                               ncols, gamma, p, q);

    compute_prim_smi_three_center_electron_repulsion_0_piece13(buffer, target, pc, slh, smg0,
                                                               smg1, smh, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
