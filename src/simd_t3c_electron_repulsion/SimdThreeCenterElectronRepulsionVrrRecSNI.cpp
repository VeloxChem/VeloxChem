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


#include "SimdThreeCenterElectronRepulsionVrrRecSNI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;

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

    const auto *smi0_0 = buffer.data(smi0 + 0);
    const auto *smi0_3 = buffer.data(smi0 + 3);
    const auto *smi0_5 = buffer.data(smi0 + 5);
    const auto *smi0_6 = buffer.data(smi0 + 6);
    const auto *smi0_9 = buffer.data(smi0 + 9);
    const auto *smi0_10 = buffer.data(smi0 + 10);
    const auto *smi0_12 = buffer.data(smi0 + 12);
    const auto *smi0_14 = buffer.data(smi0 + 14);
    const auto *smi0_21 = buffer.data(smi0 + 21);
    const auto *smi0_27 = buffer.data(smi0 + 27);
    const auto *smi0_31 = buffer.data(smi0 + 31);
    const auto *smi0_34 = buffer.data(smi0 + 34);
    const auto *smi0_38 = buffer.data(smi0 + 38);
    const auto *smi0_56 = buffer.data(smi0 + 56);
    const auto *smi0_61 = buffer.data(smi0 + 61);
    const auto *smi0_65 = buffer.data(smi0 + 65);

    const auto *smh_0 = buffer.data(smh + 0);
    const auto *smh_1 = buffer.data(smh + 1);
    const auto *smh_2 = buffer.data(smh + 2);
    const auto *smh_3 = buffer.data(smh + 3);
    const auto *smh_5 = buffer.data(smh + 5);
    const auto *smh_6 = buffer.data(smh + 6);
    const auto *smh_7 = buffer.data(smh + 7);
    const auto *smh_8 = buffer.data(smh + 8);
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
    const auto *smh_47 = buffer.data(smh + 47);
    const auto *smh_57 = buffer.data(smh + 57);
    const auto *smh_58 = buffer.data(smh + 58);
    const auto *smh_59 = buffer.data(smh + 59);
    const auto *smh_60 = buffer.data(smh + 60);
    const auto *smh_61 = buffer.data(smh + 61);
    const auto *smh_62 = buffer.data(smh + 62);
    const auto *smh_63 = buffer.data(smh + 63);
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

    const auto *smi1_0 = buffer.data(smi1 + 0);
    const auto *smi1_3 = buffer.data(smi1 + 3);
    const auto *smi1_5 = buffer.data(smi1 + 5);
    const auto *smi1_6 = buffer.data(smi1 + 6);
    const auto *smi1_9 = buffer.data(smi1 + 9);
    const auto *smi1_10 = buffer.data(smi1 + 10);
    const auto *smi1_12 = buffer.data(smi1 + 12);
    const auto *smi1_14 = buffer.data(smi1 + 14);
    const auto *smi1_21 = buffer.data(smi1 + 21);
    const auto *smi1_27 = buffer.data(smi1 + 27);
    const auto *smi1_31 = buffer.data(smi1 + 31);
    const auto *smi1_34 = buffer.data(smi1 + 34);
    const auto *smi1_38 = buffer.data(smi1 + 38);
    const auto *smi1_56 = buffer.data(smi1 + 56);
    const auto *smi1_61 = buffer.data(smi1 + 61);
    const auto *smi1_65 = buffer.data(smi1 + 65);

    const auto *sng0_0 = buffer.data(sng0 + 0);
    const auto *sng0_3 = buffer.data(sng0 + 3);
    const auto *sng0_5 = buffer.data(sng0 + 5);
    const auto *sng0_6 = buffer.data(sng0 + 6);
    const auto *sng0_9 = buffer.data(sng0 + 9);
    const auto *sng0_10 = buffer.data(sng0 + 10);
    const auto *sng0_12 = buffer.data(sng0 + 12);
    const auto *sng0_13 = buffer.data(sng0 + 13);
    const auto *sng0_14 = buffer.data(sng0 + 14);
    const auto *sng0_25 = buffer.data(sng0 + 25);
    const auto *sng0_27 = buffer.data(sng0 + 27);
    const auto *sng0_28 = buffer.data(sng0 + 28);
    const auto *sng0_29 = buffer.data(sng0 + 29);
    const auto *sng0_42 = buffer.data(sng0 + 42);
    const auto *sng0_43 = buffer.data(sng0 + 43);
    const auto *sng0_44 = buffer.data(sng0 + 44);
    const auto *sng0_45 = buffer.data(sng0 + 45);
    const auto *sng0_48 = buffer.data(sng0 + 48);
    const auto *sng0_50 = buffer.data(sng0 + 50);
    const auto *sng0_51 = buffer.data(sng0 + 51);
    const auto *sng0_54 = buffer.data(sng0 + 54);
    const auto *sng0_55 = buffer.data(sng0 + 55);
    const auto *sng0_57 = buffer.data(sng0 + 57);
    const auto *sng0_58 = buffer.data(sng0 + 58);
    const auto *sng0_59 = buffer.data(sng0 + 59);

    const auto *sng1_0 = buffer.data(sng1 + 0);
    const auto *sng1_3 = buffer.data(sng1 + 3);
    const auto *sng1_5 = buffer.data(sng1 + 5);
    const auto *sng1_6 = buffer.data(sng1 + 6);
    const auto *sng1_9 = buffer.data(sng1 + 9);
    const auto *sng1_10 = buffer.data(sng1 + 10);
    const auto *sng1_12 = buffer.data(sng1 + 12);
    const auto *sng1_13 = buffer.data(sng1 + 13);
    const auto *sng1_14 = buffer.data(sng1 + 14);
    const auto *sng1_25 = buffer.data(sng1 + 25);
    const auto *sng1_27 = buffer.data(sng1 + 27);
    const auto *sng1_28 = buffer.data(sng1 + 28);
    const auto *sng1_29 = buffer.data(sng1 + 29);
    const auto *sng1_42 = buffer.data(sng1 + 42);
    const auto *sng1_43 = buffer.data(sng1 + 43);
    const auto *sng1_44 = buffer.data(sng1 + 44);
    const auto *sng1_45 = buffer.data(sng1 + 45);
    const auto *sng1_48 = buffer.data(sng1 + 48);
    const auto *sng1_50 = buffer.data(sng1 + 50);
    const auto *sng1_51 = buffer.data(sng1 + 51);
    const auto *sng1_54 = buffer.data(sng1 + 54);
    const auto *sng1_55 = buffer.data(sng1 + 55);
    const auto *sng1_57 = buffer.data(sng1 + 57);
    const auto *sng1_58 = buffer.data(sng1 + 58);
    const auto *sng1_59 = buffer.data(sng1 + 59);

    const auto *snh_0 = buffer.data(snh + 0);
    const auto *snh_2 = buffer.data(snh + 2);
    const auto *snh_3 = buffer.data(snh + 3);
    const auto *snh_5 = buffer.data(snh + 5);
    const auto *snh_6 = buffer.data(snh + 6);
    const auto *snh_9 = buffer.data(snh + 9);
    const auto *snh_10 = buffer.data(snh + 10);
    const auto *snh_12 = buffer.data(snh + 12);
    const auto *snh_14 = buffer.data(snh + 14);
    const auto *snh_15 = buffer.data(snh + 15);
    const auto *snh_16 = buffer.data(snh + 16);
    const auto *snh_17 = buffer.data(snh + 17);
    const auto *snh_18 = buffer.data(snh + 18);
    const auto *snh_19 = buffer.data(snh + 19);
    const auto *snh_20 = buffer.data(snh + 20);
    const auto *snh_21 = buffer.data(snh + 21);
    const auto *snh_23 = buffer.data(snh + 23);
    const auto *snh_24 = buffer.data(snh + 24);
    const auto *snh_26 = buffer.data(snh + 26);
    const auto *snh_27 = buffer.data(snh + 27);
    const auto *snh_30 = buffer.data(snh + 30);
    const auto *snh_36 = buffer.data(snh + 36);
    const auto *snh_37 = buffer.data(snh + 37);
    const auto *snh_38 = buffer.data(snh + 38);
    const auto *snh_39 = buffer.data(snh + 39);
    const auto *snh_40 = buffer.data(snh + 40);
    const auto *snh_41 = buffer.data(snh + 41);
    const auto *snh_42 = buffer.data(snh + 42);
    const auto *snh_44 = buffer.data(snh + 44);
    const auto *snh_45 = buffer.data(snh + 45);
    const auto *snh_47 = buffer.data(snh + 47);
    const auto *snh_48 = buffer.data(snh + 48);
    const auto *snh_51 = buffer.data(snh + 51);
    const auto *snh_57 = buffer.data(snh + 57);
    const auto *snh_58 = buffer.data(snh + 58);
    const auto *snh_59 = buffer.data(snh + 59);
    const auto *snh_60 = buffer.data(snh + 60);
    const auto *snh_61 = buffer.data(snh + 61);
    const auto *snh_62 = buffer.data(snh + 62);
    const auto *snh_63 = buffer.data(snh + 63);
    const auto *snh_65 = buffer.data(snh + 65);
    const auto *snh_66 = buffer.data(snh + 66);
    const auto *snh_68 = buffer.data(snh + 68);
    const auto *snh_69 = buffer.data(snh + 69);
    const auto *snh_72 = buffer.data(snh + 72);
    const auto *snh_73 = buffer.data(snh + 73);
    const auto *snh_75 = buffer.data(snh + 75);
    const auto *snh_77 = buffer.data(snh + 77);
    const auto *snh_78 = buffer.data(snh + 78);
    const auto *snh_79 = buffer.data(snh + 79);
    const auto *snh_80 = buffer.data(snh + 80);
    const auto *snh_81 = buffer.data(snh + 81);
    const auto *snh_82 = buffer.data(snh + 82);
    const auto *snh_83 = buffer.data(snh + 83);
    const auto *snh_84 = buffer.data(snh + 84);
    const auto *snh_86 = buffer.data(snh + 86);
    const auto *snh_87 = buffer.data(snh + 87);
    const auto *snh_89 = buffer.data(snh + 89);
    const auto *snh_90 = buffer.data(snh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, smh_0, smh_3, sng0_0, sng0_3, \
                         sng1_0, sng1_3, snh_0, snh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * smh_0[k]
                 + f_1 * sng0_0[k]
                 - f_2 * sng1_0[k]
                 + f_3 * pc_x[k] * snh_0[k];

        t_1[k] = f_3 * pc_y[k] * snh_0[k];

        t_2[k] = f_3 * pc_z[k] * snh_0[k];

        t_3[k] = f_0 * smh_3[k]
                 + f_4 * sng0_3[k]
                 - f_5 * sng1_3[k]
                 + f_3 * pc_x[k] * snh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, smh_5, smh_6, sng0_5, sng0_6, sng1_5, \
                         sng1_6, snh_2, snh_5, snh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * snh_2[k];

        t_5[k] = f_0 * smh_5[k]
                 + f_4 * sng0_5[k]
                 - f_5 * sng1_5[k]
                 + f_3 * pc_x[k] * snh_5[k];

        t_6[k] = f_0 * smh_6[k]
                 + f_6 * sng0_6[k]
                 - f_7 * sng1_6[k]
                 + f_3 * pc_x[k] * snh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, smh_9, sng0_9, sng1_9, snh_3, snh_5, \
                         snh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * snh_3[k];

        t_8[k] = f_3 * pc_y[k] * snh_5[k];

        t_9[k] = f_0 * smh_9[k]
                 + f_6 * sng0_9[k]
                 - f_7 * sng1_9[k]
                 + f_3 * pc_x[k] * snh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, smh_10, smh_12, sng0_10, sng0_12, \
                         sng1_10, sng1_12, snh_6, snh_10, snh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * smh_10[k]
                  + f_8 * sng0_10[k]
                  - f_9 * sng1_10[k]
                  + f_3 * pc_x[k] * snh_10[k];

        t_11[k] = f_3 * pc_z[k] * snh_6[k];

        t_12[k] = f_0 * smh_12[k]
                  + f_8 * sng0_12[k]
                  - f_9 * sng1_12[k]
                  + f_3 * pc_x[k] * snh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, smh_14, smh_15, smh_16, sng0_14, \
                         sng1_14, snh_9, snh_14, snh_15, snh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * snh_9[k];

        t_14[k] = f_0 * smh_14[k]
                  + f_8 * sng0_14[k]
                  - f_9 * sng1_14[k]
                  + f_3 * pc_x[k] * snh_14[k];

        t_15[k] = f_0 * smh_15[k]
                  + f_3 * pc_x[k] * snh_15[k];

        t_16[k] = f_0 * smh_16[k]
                  + f_3 * pc_x[k] * snh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, smh_17, smh_18, smh_19, smh_20, snh_17, \
                         snh_18, snh_19, snh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * smh_17[k]
                  + f_3 * pc_x[k] * snh_17[k];

        t_18[k] = f_0 * smh_18[k]
                  + f_3 * pc_x[k] * snh_18[k];

        t_19[k] = f_0 * smh_19[k]
                  + f_3 * pc_x[k] * snh_19[k];

        t_20[k] = f_0 * smh_20[k]
                  + f_3 * pc_x[k] * snh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, sng0_10, sng0_12, sng0_13, \
                         sng1_10, sng1_12, sng1_13, snh_15, snh_17, \
                         snh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * sng0_10[k]
                  - f_2 * sng1_10[k]
                  + f_3 * pc_y[k] * snh_15[k];

        t_22[k] = f_3 * pc_z[k] * snh_15[k];

        t_23[k] = f_4 * sng0_12[k]
                  - f_5 * sng1_12[k]
                  + f_3 * pc_y[k] * snh_17[k];

        t_24[k] = f_6 * sng0_13[k]
                  - f_7 * sng1_13[k]
                  + f_3 * pc_y[k] * snh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, smi0_0, smh_0, \
                         smi1_0, sng0_14, sng1_14, snh_19, snh_20, \
                         snh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * sng0_14[k]
                  - f_9 * sng1_14[k]
                  + f_3 * pc_y[k] * snh_19[k];

        t_26[k] = f_3 * pc_y[k] * snh_20[k];

        t_27[k] = f_1 * sng0_14[k]
                  - f_2 * sng1_14[k]
                  + f_3 * pc_z[k] * snh_20[k];

        t_28[k] = pb_y[k] * smi0_0[k]
                  - f_10 * pc_y[k] * smi1_0[k];

        t_29[k] = f_11 * smh_0[k]
                  + f_3 * pc_y[k] * snh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, smi0_3, smi0_5, smh_1, \
                         smh_2, smi1_3, smi1_5, snh_21, snh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * snh_21[k];

        t_31[k] = pb_y[k] * smi0_3[k]
                  + f_12 * smh_1[k]
                  - f_10 * pc_y[k] * smi1_3[k];

        t_32[k] = f_11 * smh_2[k]
                  + f_3 * pc_y[k] * snh_23[k];

        t_33[k] = pb_y[k] * smi0_5[k]
                  - f_10 * pc_y[k] * smi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, smi0_6, smi0_9, smh_3, \
                         smh_5, smi1_6, smi1_9, snh_24, snh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * smi0_6[k]
                  + f_13 * smh_3[k]
                  - f_10 * pc_y[k] * smi1_6[k];

        t_35[k] = f_3 * pc_z[k] * snh_24[k];

        t_36[k] = f_11 * smh_5[k]
                  + f_3 * pc_y[k] * snh_26[k];

        t_37[k] = pb_y[k] * smi0_9[k]
                  - f_10 * pc_y[k] * smi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, smi0_10, smi0_12, smh_6, \
                         smh_8, smh_9, smi1_10, smi1_12, snh_27, \
                         snh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * smi0_10[k]
                  + f_14 * smh_6[k]
                  - f_10 * pc_y[k] * smi1_10[k];

        t_39[k] = f_3 * pc_z[k] * snh_27[k];

        t_40[k] = pb_y[k] * smi0_12[k]
                  + f_12 * smh_8[k]
                  - f_10 * pc_y[k] * smi1_12[k];

        t_41[k] = f_11 * smh_9[k]
                  + f_3 * pc_y[k] * snh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, smi0_14, smh_36, smh_37, \
                         smh_38, smi1_14, snh_36, snh_37, snh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * smi0_14[k]
                  - f_10 * pc_y[k] * smi1_14[k];

        t_43[k] = f_15 * smh_36[k]
                  + f_3 * pc_x[k] * snh_36[k];

        t_44[k] = f_15 * smh_37[k]
                  + f_3 * pc_x[k] * snh_37[k];

        t_45[k] = f_15 * smh_38[k]
                  + f_3 * pc_x[k] * snh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, smh_15, smh_39, smh_40, smh_41, \
                         sng0_25, sng1_25, snh_36, snh_39, snh_40, \
                         snh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * smh_39[k]
                  + f_3 * pc_x[k] * snh_39[k];

        t_47[k] = f_15 * smh_40[k]
                  + f_3 * pc_x[k] * snh_40[k];

        t_48[k] = f_15 * smh_41[k]
                  + f_3 * pc_x[k] * snh_41[k];

        t_49[k] = f_11 * smh_15[k]
                  + f_1 * sng0_25[k]
                  - f_2 * sng1_25[k]
                  + f_3 * pc_y[k] * snh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, smh_17, smh_18, sng0_27, sng0_28, \
                         sng1_27, sng1_28, snh_36, snh_38, snh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * snh_36[k];

        t_51[k] = f_11 * smh_17[k]
                  + f_4 * sng0_27[k]
                  - f_5 * sng1_27[k]
                  + f_3 * pc_y[k] * snh_38[k];

        t_52[k] = f_11 * smh_18[k]
                  + f_6 * sng0_28[k]
                  - f_7 * sng1_28[k]
                  + f_3 * pc_y[k] * snh_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, smi0_27, smh_19, smh_20, smi1_27, \
                         sng0_29, sng1_29, snh_40, snh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * smh_19[k]
                  + f_8 * sng0_29[k]
                  - f_9 * sng1_29[k]
                  + f_3 * pc_y[k] * snh_40[k];

        t_54[k] = f_11 * smh_20[k]
                  + f_3 * pc_y[k] * snh_41[k];

        t_55[k] = pb_y[k] * smi0_27[k]
                  - f_10 * pc_y[k] * smi1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, smi0_0, smi0_3, \
                         smh_0, smi1_0, smi1_3, snh_42, snh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * smi0_0[k]
                  - f_10 * pc_z[k] * smi1_0[k];

        t_57[k] = f_3 * pc_y[k] * snh_42[k];

        t_58[k] = f_11 * smh_0[k]
                  + f_3 * pc_z[k] * snh_42[k];

        t_59[k] = pb_z[k] * smi0_3[k]
                  - f_10 * pc_z[k] * smi1_3[k];

        t_60[k] = f_3 * pc_y[k] * snh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, smi0_5, smi0_6, smh_2, \
                         smh_3, smi1_5, smi1_6, snh_45, snh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * smi0_5[k]
                  + f_12 * smh_2[k]
                  - f_10 * pc_z[k] * smi1_5[k];

        t_62[k] = pb_z[k] * smi0_6[k]
                  - f_10 * pc_z[k] * smi1_6[k];

        t_63[k] = f_11 * smh_3[k]
                  + f_3 * pc_z[k] * snh_45[k];

        t_64[k] = f_3 * pc_y[k] * snh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, smi0_9, smi0_10, smi0_12, smh_5, \
                         smh_6, smh_7, smi1_9, smi1_10, smi1_12, \
                         snh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * smi0_9[k]
                  + f_13 * smh_5[k]
                  - f_10 * pc_z[k] * smi1_9[k];

        t_66[k] = pb_z[k] * smi0_10[k]
                  - f_10 * pc_z[k] * smi1_10[k];

        t_67[k] = f_11 * smh_6[k]
                  + f_3 * pc_z[k] * snh_48[k];

        t_68[k] = pb_z[k] * smi0_12[k]
                  + f_12 * smh_7[k]
                  - f_10 * pc_z[k] * smi1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, smi0_14, smh_9, \
                         smh_57, smh_58, smi1_14, snh_51, snh_57, \
                         snh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * snh_51[k];

        t_70[k] = pb_z[k] * smi0_14[k]
                  + f_14 * smh_9[k]
                  - f_10 * pc_z[k] * smi1_14[k];

        t_71[k] = f_15 * smh_57[k]
                  + f_3 * pc_x[k] * snh_57[k];

        t_72[k] = f_15 * smh_58[k]
                  + f_3 * pc_x[k] * snh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, smh_59, smh_60, smh_61, smh_62, snh_59, \
                         snh_60, snh_61, snh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * smh_59[k]
                  + f_3 * pc_x[k] * snh_59[k];

        t_74[k] = f_15 * smh_60[k]
                  + f_3 * pc_x[k] * snh_60[k];

        t_75[k] = f_15 * smh_61[k]
                  + f_3 * pc_x[k] * snh_61[k];

        t_76[k] = f_15 * smh_62[k]
                  + f_3 * pc_x[k] * snh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, smi0_21, smh_15, smi1_21, \
                         sng0_42, sng1_42, snh_57, snh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * smi0_21[k]
                  - f_10 * pc_z[k] * smi1_21[k];

        t_78[k] = f_11 * smh_15[k]
                  + f_3 * pc_z[k] * snh_57[k];

        t_79[k] = f_4 * sng0_42[k]
                  - f_5 * sng1_42[k]
                  + f_3 * pc_y[k] * snh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, smh_20, sng0_43, sng0_44, \
                         sng1_43, sng1_44, snh_60, snh_61, snh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * sng0_43[k]
                  - f_7 * sng1_43[k]
                  + f_3 * pc_y[k] * snh_60[k];

        t_81[k] = f_8 * sng0_44[k]
                  - f_9 * sng1_44[k]
                  + f_3 * pc_y[k] * snh_61[k];

        t_82[k] = f_3 * pc_y[k] * snh_62[k];

        t_83[k] = f_11 * smh_20[k]
                  + f_1 * sng0_44[k]
                  - f_2 * sng1_44[k]
                  + f_3 * pc_z[k] * snh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, smh_21, smh_63, smh_66, \
                         sng0_45, sng0_48, sng1_45, sng1_48, snh_63, \
                         snh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_16 * smh_63[k]
                  + f_1 * sng0_45[k]
                  - f_2 * sng1_45[k]
                  + f_3 * pc_x[k] * snh_63[k];

        t_85[k] = f_12 * smh_21[k]
                  + f_3 * pc_y[k] * snh_63[k];

        t_86[k] = f_3 * pc_z[k] * snh_63[k];

        t_87[k] = f_16 * smh_66[k]
                  + f_4 * sng0_48[k]
                  - f_5 * sng1_48[k]
                  + f_3 * pc_x[k] * snh_66[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pc_x, pc_y, smh_23, smh_68, smh_69, sng0_50, \
                         sng0_51, sng1_50, sng1_51, snh_65, snh_68, \
                         snh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * smh_23[k]
                  + f_3 * pc_y[k] * snh_65[k];

        t_89[k] = f_16 * smh_68[k]
                  + f_4 * sng0_50[k]
                  - f_5 * sng1_50[k]
                  + f_3 * pc_x[k] * snh_68[k];

        t_90[k] = f_16 * smh_69[k]
                  + f_6 * sng0_51[k]
                  - f_7 * sng1_51[k]
                  + f_3 * pc_x[k] * snh_69[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pc_x, pc_y, pc_z, smh_26, smh_72, sng0_54, sng1_54, \
                         snh_66, snh_68, snh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * snh_66[k];

        t_92[k] = f_12 * smh_26[k]
                  + f_3 * pc_y[k] * snh_68[k];

        t_93[k] = f_16 * smh_72[k]
                  + f_6 * sng0_54[k]
                  - f_7 * sng1_54[k]
                  + f_3 * pc_x[k] * snh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pc_x, pc_z, smh_73, smh_75, sng0_55, sng0_57, \
                         sng1_55, sng1_57, snh_69, snh_73, snh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_16 * smh_73[k]
                  + f_8 * sng0_55[k]
                  - f_9 * sng1_55[k]
                  + f_3 * pc_x[k] * snh_73[k];

        t_95[k] = f_3 * pc_z[k] * snh_69[k];

        t_96[k] = f_16 * smh_75[k]
                  + f_8 * sng0_57[k]
                  - f_9 * sng1_57[k]
                  + f_3 * pc_x[k] * snh_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, smh_30, smh_77, smh_78, smh_79, \
                         sng0_59, sng1_59, snh_72, snh_77, snh_78, \
                         snh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * smh_30[k]
                  + f_3 * pc_y[k] * snh_72[k];

        t_98[k] = f_16 * smh_77[k]
                  + f_8 * sng0_59[k]
                  - f_9 * sng1_59[k]
                  + f_3 * pc_x[k] * snh_77[k];

        t_99[k] = f_16 * smh_78[k]
                  + f_3 * pc_x[k] * snh_78[k];

        t_100[k] = f_16 * smh_79[k]
                   + f_3 * pc_x[k] * snh_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, smh_80, smh_81, smh_82, smh_83, \
                         snh_80, snh_81, snh_82, snh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_16 * smh_80[k]
                   + f_3 * pc_x[k] * snh_80[k];

        t_102[k] = f_16 * smh_81[k]
                   + f_3 * pc_x[k] * snh_81[k];

        t_103[k] = f_16 * smh_82[k]
                   + f_3 * pc_x[k] * snh_82[k];

        t_104[k] = f_16 * smh_83[k]
                   + f_3 * pc_x[k] * snh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, smh_36, smh_38, sng0_55, sng0_57, \
                         sng1_55, sng1_57, snh_78, snh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * smh_36[k]
                   + f_1 * sng0_55[k]
                   - f_2 * sng1_55[k]
                   + f_3 * pc_y[k] * snh_78[k];

        t_106[k] = f_3 * pc_z[k] * snh_78[k];

        t_107[k] = f_12 * smh_38[k]
                   + f_4 * sng0_57[k]
                   - f_5 * sng1_57[k]
                   + f_3 * pc_y[k] * snh_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, smh_39, smh_40, smh_41, \
                         sng0_58, sng0_59, sng1_58, sng1_59, snh_81, snh_82, \
                         snh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_12 * smh_39[k]
                   + f_6 * sng0_58[k]
                   - f_7 * sng1_58[k]
                   + f_3 * pc_y[k] * snh_81[k];

        t_109[k] = f_12 * smh_40[k]
                   + f_8 * sng0_59[k]
                   - f_9 * sng1_59[k]
                   + f_3 * pc_y[k] * snh_82[k];

        t_110[k] = f_12 * smh_41[k]
                   + f_3 * pc_y[k] * snh_83[k];

        t_111[k] = f_1 * sng0_59[k]
                   - f_2 * sng1_59[k]
                   + f_3 * pc_z[k] * snh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, smi0_31, smi0_56, \
                         smh_21, smh_42, smi1_31, smi1_56, snh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * smi0_56[k]
                   - f_10 * pc_y[k] * smi1_56[k];

        t_113[k] = f_11 * smh_42[k]
                   + f_3 * pc_y[k] * snh_84[k];

        t_114[k] = f_11 * smh_21[k]
                   + f_3 * pc_z[k] * snh_84[k];

        t_115[k] = pb_z[k] * smi0_31[k]
                   - f_10 * pc_z[k] * smi1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, smi0_34, smi0_61, \
                         smh_24, smh_44, smi1_34, smi1_61, snh_86, \
                         snh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * smh_44[k]
                   + f_3 * pc_y[k] * snh_86[k];

        t_117[k] = pb_y[k] * smi0_61[k]
                   - f_10 * pc_y[k] * smi1_61[k];

        t_118[k] = pb_z[k] * smi0_34[k]
                   - f_10 * pc_z[k] * smi1_34[k];

        t_119[k] = f_11 * smh_24[k]
                   + f_3 * pc_z[k] * snh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, smi0_38, smi0_65, \
                         smh_27, smh_47, smi1_38, smi1_65, snh_89, \
                         snh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * smh_47[k]
                   + f_3 * pc_y[k] * snh_89[k];

        t_121[k] = pb_y[k] * smi0_65[k]
                   - f_10 * pc_y[k] * smi1_65[k];

        t_122[k] = pb_z[k] * smi0_38[k]
                   - f_10 * pc_z[k] * smi1_38[k];

        t_123[k] = f_11 * smh_27[k]
                   + f_3 * pc_z[k] * snh_90[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;

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

    const auto *smi0_49 = buffer.data(smi0 + 49);
    const auto *smi0_68 = buffer.data(smi0 + 68);
    const auto *smi0_70 = buffer.data(smi0 + 70);
    const auto *smi0_83 = buffer.data(smi0 + 83);
    const auto *smi0_84 = buffer.data(smi0 + 84);
    const auto *smi0_87 = buffer.data(smi0 + 87);
    const auto *smi0_90 = buffer.data(smi0 + 90);
    const auto *smi0_94 = buffer.data(smi0 + 94);
    const auto *smi0_96 = buffer.data(smi0 + 96);
    const auto *smi0_105 = buffer.data(smi0 + 105);
    const auto *smi0_140 = buffer.data(smi0 + 140);
    const auto *smi0_143 = buffer.data(smi0 + 143);
    const auto *smi0_145 = buffer.data(smi0 + 145);
    const auto *smi0_146 = buffer.data(smi0 + 146);
    const auto *smi0_149 = buffer.data(smi0 + 149);
    const auto *smi0_150 = buffer.data(smi0 + 150);
    const auto *smi0_152 = buffer.data(smi0 + 152);
    const auto *smi0_154 = buffer.data(smi0 + 154);

    const auto *smh_36 = buffer.data(smh + 36);
    const auto *smh_42 = buffer.data(smh + 42);
    const auto *smh_45 = buffer.data(smh + 45);
    const auto *smh_48 = buffer.data(smh + 48);
    const auto *smh_50 = buffer.data(smh + 50);
    const auto *smh_51 = buffer.data(smh + 51);
    const auto *smh_57 = buffer.data(smh + 57);
    const auto *smh_59 = buffer.data(smh + 59);
    const auto *smh_60 = buffer.data(smh + 60);
    const auto *smh_61 = buffer.data(smh + 61);
    const auto *smh_62 = buffer.data(smh + 62);
    const auto *smh_63 = buffer.data(smh + 63);
    const auto *smh_65 = buffer.data(smh + 65);
    const auto *smh_66 = buffer.data(smh + 66);
    const auto *smh_68 = buffer.data(smh + 68);
    const auto *smh_69 = buffer.data(smh + 69);
    const auto *smh_70 = buffer.data(smh + 70);
    const auto *smh_72 = buffer.data(smh + 72);
    const auto *smh_78 = buffer.data(smh + 78);
    const auto *smh_80 = buffer.data(smh + 80);
    const auto *smh_81 = buffer.data(smh + 81);
    const auto *smh_82 = buffer.data(smh + 82);
    const auto *smh_83 = buffer.data(smh + 83);
    const auto *smh_84 = buffer.data(smh + 84);
    const auto *smh_86 = buffer.data(smh + 86);
    const auto *smh_87 = buffer.data(smh + 87);
    const auto *smh_89 = buffer.data(smh + 89);
    const auto *smh_90 = buffer.data(smh + 90);
    const auto *smh_93 = buffer.data(smh + 93);
    const auto *smh_99 = buffer.data(smh + 99);
    const auto *smh_100 = buffer.data(smh + 100);
    const auto *smh_101 = buffer.data(smh + 101);
    const auto *smh_102 = buffer.data(smh + 102);
    const auto *smh_103 = buffer.data(smh + 103);
    const auto *smh_104 = buffer.data(smh + 104);
    const auto *smh_105 = buffer.data(smh + 105);
    const auto *smh_106 = buffer.data(smh + 106);
    const auto *smh_107 = buffer.data(smh + 107);
    const auto *smh_108 = buffer.data(smh + 108);
    const auto *smh_110 = buffer.data(smh + 110);
    const auto *smh_111 = buffer.data(smh + 111);
    const auto *smh_113 = buffer.data(smh + 113);
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
    const auto *smh_152 = buffer.data(smh + 152);
    const auto *smh_156 = buffer.data(smh + 156);
    const auto *smh_161 = buffer.data(smh + 161);
    const auto *smh_162 = buffer.data(smh + 162);
    const auto *smh_163 = buffer.data(smh + 163);
    const auto *smh_164 = buffer.data(smh + 164);
    const auto *smh_165 = buffer.data(smh + 165);
    const auto *smh_166 = buffer.data(smh + 166);
    const auto *smh_167 = buffer.data(smh + 167);
    const auto *smh_183 = buffer.data(smh + 183);

    const auto *smi1_49 = buffer.data(smi1 + 49);
    const auto *smi1_68 = buffer.data(smi1 + 68);
    const auto *smi1_70 = buffer.data(smi1 + 70);
    const auto *smi1_83 = buffer.data(smi1 + 83);
    const auto *smi1_84 = buffer.data(smi1 + 84);
    const auto *smi1_87 = buffer.data(smi1 + 87);
    const auto *smi1_90 = buffer.data(smi1 + 90);
    const auto *smi1_94 = buffer.data(smi1 + 94);
    const auto *smi1_96 = buffer.data(smi1 + 96);
    const auto *smi1_105 = buffer.data(smi1 + 105);
    const auto *smi1_140 = buffer.data(smi1 + 140);
    const auto *smi1_143 = buffer.data(smi1 + 143);
    const auto *smi1_145 = buffer.data(smi1 + 145);
    const auto *smi1_146 = buffer.data(smi1 + 146);
    const auto *smi1_149 = buffer.data(smi1 + 149);
    const auto *smi1_150 = buffer.data(smi1 + 150);
    const auto *smi1_152 = buffer.data(smi1 + 152);
    const auto *smi1_154 = buffer.data(smi1 + 154);

    const auto *sng0_72 = buffer.data(sng0 + 72);
    const auto *sng0_73 = buffer.data(sng0 + 73);
    const auto *sng0_74 = buffer.data(sng0 + 74);
    const auto *sng0_75 = buffer.data(sng0 + 75);
    const auto *sng0_78 = buffer.data(sng0 + 78);
    const auto *sng0_80 = buffer.data(sng0 + 80);
    const auto *sng0_81 = buffer.data(sng0 + 81);
    const auto *sng0_84 = buffer.data(sng0 + 84);
    const auto *sng0_85 = buffer.data(sng0 + 85);
    const auto *sng0_87 = buffer.data(sng0 + 87);
    const auto *sng0_88 = buffer.data(sng0 + 88);
    const auto *sng0_89 = buffer.data(sng0 + 89);
    const auto *sng0_90 = buffer.data(sng0 + 90);
    const auto *sng0_93 = buffer.data(sng0 + 93);
    const auto *sng0_95 = buffer.data(sng0 + 95);
    const auto *sng0_96 = buffer.data(sng0 + 96);
    const auto *sng0_99 = buffer.data(sng0 + 99);
    const auto *sng0_100 = buffer.data(sng0 + 100);
    const auto *sng0_102 = buffer.data(sng0 + 102);
    const auto *sng0_103 = buffer.data(sng0 + 103);
    const auto *sng0_104 = buffer.data(sng0 + 104);
    const auto *sng0_110 = buffer.data(sng0 + 110);
    const auto *sng0_114 = buffer.data(sng0 + 114);
    const auto *sng0_117 = buffer.data(sng0 + 117);
    const auto *sng0_118 = buffer.data(sng0 + 118);
    const auto *sng0_119 = buffer.data(sng0 + 119);

    const auto *sng1_72 = buffer.data(sng1 + 72);
    const auto *sng1_73 = buffer.data(sng1 + 73);
    const auto *sng1_74 = buffer.data(sng1 + 74);
    const auto *sng1_75 = buffer.data(sng1 + 75);
    const auto *sng1_78 = buffer.data(sng1 + 78);
    const auto *sng1_80 = buffer.data(sng1 + 80);
    const auto *sng1_81 = buffer.data(sng1 + 81);
    const auto *sng1_84 = buffer.data(sng1 + 84);
    const auto *sng1_85 = buffer.data(sng1 + 85);
    const auto *sng1_87 = buffer.data(sng1 + 87);
    const auto *sng1_88 = buffer.data(sng1 + 88);
    const auto *sng1_89 = buffer.data(sng1 + 89);
    const auto *sng1_90 = buffer.data(sng1 + 90);
    const auto *sng1_93 = buffer.data(sng1 + 93);
    const auto *sng1_95 = buffer.data(sng1 + 95);
    const auto *sng1_96 = buffer.data(sng1 + 96);
    const auto *sng1_99 = buffer.data(sng1 + 99);
    const auto *sng1_100 = buffer.data(sng1 + 100);
    const auto *sng1_102 = buffer.data(sng1 + 102);
    const auto *sng1_103 = buffer.data(sng1 + 103);
    const auto *sng1_104 = buffer.data(sng1 + 104);
    const auto *sng1_110 = buffer.data(sng1 + 110);
    const auto *sng1_114 = buffer.data(sng1 + 114);
    const auto *sng1_117 = buffer.data(sng1 + 117);
    const auto *sng1_118 = buffer.data(sng1 + 118);
    const auto *sng1_119 = buffer.data(sng1 + 119);

    const auto *snh_93 = buffer.data(snh + 93);
    const auto *snh_99 = buffer.data(snh + 99);
    const auto *snh_100 = buffer.data(snh + 100);
    const auto *snh_101 = buffer.data(snh + 101);
    const auto *snh_102 = buffer.data(snh + 102);
    const auto *snh_103 = buffer.data(snh + 103);
    const auto *snh_104 = buffer.data(snh + 104);
    const auto *snh_105 = buffer.data(snh + 105);
    const auto *snh_107 = buffer.data(snh + 107);
    const auto *snh_108 = buffer.data(snh + 108);
    const auto *snh_110 = buffer.data(snh + 110);
    const auto *snh_111 = buffer.data(snh + 111);
    const auto *snh_114 = buffer.data(snh + 114);
    const auto *snh_115 = buffer.data(snh + 115);
    const auto *snh_117 = buffer.data(snh + 117);
    const auto *snh_119 = buffer.data(snh + 119);
    const auto *snh_120 = buffer.data(snh + 120);
    const auto *snh_121 = buffer.data(snh + 121);
    const auto *snh_122 = buffer.data(snh + 122);
    const auto *snh_123 = buffer.data(snh + 123);
    const auto *snh_124 = buffer.data(snh + 124);
    const auto *snh_125 = buffer.data(snh + 125);
    const auto *snh_126 = buffer.data(snh + 126);
    const auto *snh_128 = buffer.data(snh + 128);
    const auto *snh_129 = buffer.data(snh + 129);
    const auto *snh_131 = buffer.data(snh + 131);
    const auto *snh_132 = buffer.data(snh + 132);
    const auto *snh_135 = buffer.data(snh + 135);
    const auto *snh_136 = buffer.data(snh + 136);
    const auto *snh_138 = buffer.data(snh + 138);
    const auto *snh_140 = buffer.data(snh + 140);
    const auto *snh_141 = buffer.data(snh + 141);
    const auto *snh_142 = buffer.data(snh + 142);
    const auto *snh_143 = buffer.data(snh + 143);
    const auto *snh_144 = buffer.data(snh + 144);
    const auto *snh_145 = buffer.data(snh + 145);
    const auto *snh_146 = buffer.data(snh + 146);
    const auto *snh_147 = buffer.data(snh + 147);
    const auto *snh_149 = buffer.data(snh + 149);
    const auto *snh_150 = buffer.data(snh + 150);
    const auto *snh_152 = buffer.data(snh + 152);
    const auto *snh_153 = buffer.data(snh + 153);
    const auto *snh_156 = buffer.data(snh + 156);
    const auto *snh_161 = buffer.data(snh + 161);
    const auto *snh_162 = buffer.data(snh + 162);
    const auto *snh_163 = buffer.data(snh + 163);
    const auto *snh_164 = buffer.data(snh + 164);
    const auto *snh_165 = buffer.data(snh + 165);
    const auto *snh_166 = buffer.data(snh + 166);
    const auto *snh_167 = buffer.data(snh + 167);
    const auto *snh_168 = buffer.data(snh + 168);
    const auto *snh_170 = buffer.data(snh + 170);
    const auto *snh_171 = buffer.data(snh + 171);
    const auto *snh_173 = buffer.data(snh + 173);
    const auto *snh_174 = buffer.data(snh + 174);
    const auto *snh_177 = buffer.data(snh + 177);
    const auto *snh_183 = buffer.data(snh + 183);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, smi0_68, smi0_70, \
                         smh_50, smh_51, smh_99, smi1_68, smi1_70, snh_93, \
                         snh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * smi0_68[k]
                   + f_12 * smh_50[k]
                   - f_10 * pc_y[k] * smi1_68[k];

        t_125[k] = f_11 * smh_51[k]
                   + f_3 * pc_y[k] * snh_93[k];

        t_126[k] = pb_y[k] * smi0_70[k]
                   - f_10 * pc_y[k] * smi1_70[k];

        t_127[k] = f_16 * smh_99[k]
                   + f_3 * pc_x[k] * snh_99[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, smh_100, smh_101, smh_102, \
                         smh_103, smh_104, snh_100, snh_101, snh_102, snh_103, \
                         snh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_16 * smh_100[k]
                   + f_3 * pc_x[k] * snh_100[k];

        t_129[k] = f_16 * smh_101[k]
                   + f_3 * pc_x[k] * snh_101[k];

        t_130[k] = f_16 * smh_102[k]
                   + f_3 * pc_x[k] * snh_102[k];

        t_131[k] = f_16 * smh_103[k]
                   + f_3 * pc_x[k] * snh_103[k];

        t_132[k] = f_16 * smh_104[k]
                   + f_3 * pc_x[k] * snh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, pc_y, pc_z, smi0_49, smh_36, smh_59, \
                         smi1_49, sng0_72, sng1_72, snh_99, snh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_z[k] * smi0_49[k]
                   - f_10 * pc_z[k] * smi1_49[k];

        t_134[k] = f_11 * smh_36[k]
                   + f_3 * pc_z[k] * snh_99[k];

        t_135[k] = f_11 * smh_59[k]
                   + f_4 * sng0_72[k]
                   - f_5 * sng1_72[k]
                   + f_3 * pc_y[k] * snh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, smh_60, smh_61, smh_62, sng0_73, sng0_74, \
                         sng1_73, sng1_74, snh_102, snh_103, snh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * smh_60[k]
                   + f_6 * sng0_73[k]
                   - f_7 * sng1_73[k]
                   + f_3 * pc_y[k] * snh_102[k];

        t_137[k] = f_11 * smh_61[k]
                   + f_8 * sng0_74[k]
                   - f_9 * sng1_74[k]
                   + f_3 * pc_y[k] * snh_103[k];

        t_138[k] = f_11 * smh_62[k]
                   + f_3 * pc_y[k] * snh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_y, pc_x, pc_y, pc_z, smi0_83, smh_42, \
                         smh_105, smi1_83, sng0_75, sng1_75, snh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * smi0_83[k]
                   - f_10 * pc_y[k] * smi1_83[k];

        t_140[k] = f_16 * smh_105[k]
                   + f_1 * sng0_75[k]
                   - f_2 * sng1_75[k]
                   + f_3 * pc_x[k] * snh_105[k];

        t_141[k] = f_3 * pc_y[k] * snh_105[k];

        t_142[k] = f_12 * smh_42[k]
                   + f_3 * pc_z[k] * snh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, smh_108, smh_110, sng0_78, sng0_80, \
                         sng1_78, sng1_80, snh_107, snh_108, snh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * smh_108[k]
                   + f_4 * sng0_78[k]
                   - f_5 * sng1_78[k]
                   + f_3 * pc_x[k] * snh_108[k];

        t_144[k] = f_3 * pc_y[k] * snh_107[k];

        t_145[k] = f_16 * smh_110[k]
                   + f_4 * sng0_80[k]
                   - f_5 * sng1_80[k]
                   + f_3 * pc_x[k] * snh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_x, pc_y, pc_z, smh_45, smh_111, sng0_81, \
                         sng1_81, snh_108, snh_110, snh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_16 * smh_111[k]
                   + f_6 * sng0_81[k]
                   - f_7 * sng1_81[k]
                   + f_3 * pc_x[k] * snh_111[k];

        t_147[k] = f_12 * smh_45[k]
                   + f_3 * pc_z[k] * snh_108[k];

        t_148[k] = f_3 * pc_y[k] * snh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, smh_48, smh_114, smh_115, sng0_84, \
                         sng0_85, sng1_84, sng1_85, snh_111, snh_114, \
                         snh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_16 * smh_114[k]
                   + f_6 * sng0_84[k]
                   - f_7 * sng1_84[k]
                   + f_3 * pc_x[k] * snh_114[k];

        t_150[k] = f_16 * smh_115[k]
                   + f_8 * sng0_85[k]
                   - f_9 * sng1_85[k]
                   + f_3 * pc_x[k] * snh_115[k];

        t_151[k] = f_12 * smh_48[k]
                   + f_3 * pc_z[k] * snh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pc_x, pc_y, smh_117, smh_119, sng0_87, sng0_89, \
                         sng1_87, sng1_89, snh_114, snh_117, snh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_16 * smh_117[k]
                   + f_8 * sng0_87[k]
                   - f_9 * sng1_87[k]
                   + f_3 * pc_x[k] * snh_117[k];

        t_153[k] = f_3 * pc_y[k] * snh_114[k];

        t_154[k] = f_16 * smh_119[k]
                   + f_8 * sng0_89[k]
                   - f_9 * sng1_89[k]
                   + f_3 * pc_x[k] * snh_119[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, smh_120, smh_121, smh_122, \
                         smh_123, smh_124, snh_120, snh_121, snh_122, snh_123, \
                         snh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_16 * smh_120[k]
                   + f_3 * pc_x[k] * snh_120[k];

        t_156[k] = f_16 * smh_121[k]
                   + f_3 * pc_x[k] * snh_121[k];

        t_157[k] = f_16 * smh_122[k]
                   + f_3 * pc_x[k] * snh_122[k];

        t_158[k] = f_16 * smh_123[k]
                   + f_3 * pc_x[k] * snh_123[k];

        t_159[k] = f_16 * smh_124[k]
                   + f_3 * pc_x[k] * snh_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, smh_57, smh_125, \
                         sng0_85, sng0_87, sng1_85, sng1_87, snh_120, snh_122, \
                         snh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * smh_125[k]
                   + f_3 * pc_x[k] * snh_125[k];

        t_161[k] = f_1 * sng0_85[k]
                   - f_2 * sng1_85[k]
                   + f_3 * pc_y[k] * snh_120[k];

        t_162[k] = f_12 * smh_57[k]
                   + f_3 * pc_z[k] * snh_120[k];

        t_163[k] = f_4 * sng0_87[k]
                   - f_5 * sng1_87[k]
                   + f_3 * pc_y[k] * snh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, smh_62, sng0_88, sng0_89, \
                         sng1_88, sng1_89, snh_123, snh_124, snh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * sng0_88[k]
                   - f_7 * sng1_88[k]
                   + f_3 * pc_y[k] * snh_123[k];

        t_165[k] = f_8 * sng0_89[k]
                   - f_9 * sng1_89[k]
                   + f_3 * pc_y[k] * snh_124[k];

        t_166[k] = f_3 * pc_y[k] * snh_125[k];

        t_167[k] = f_12 * smh_62[k]
                   + f_1 * sng0_89[k]
                   - f_2 * sng1_89[k]
                   + f_3 * pc_z[k] * snh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, smh_63, smh_126, \
                         smh_129, sng0_90, sng0_93, sng1_90, sng1_93, snh_126, \
                         snh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_17 * smh_126[k]
                   + f_1 * sng0_90[k]
                   - f_2 * sng1_90[k]
                   + f_3 * pc_x[k] * snh_126[k];

        t_169[k] = f_13 * smh_63[k]
                   + f_3 * pc_y[k] * snh_126[k];

        t_170[k] = f_3 * pc_z[k] * snh_126[k];

        t_171[k] = f_17 * smh_129[k]
                   + f_4 * sng0_93[k]
                   - f_5 * sng1_93[k]
                   + f_3 * pc_x[k] * snh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pc_x, pc_y, smh_65, smh_131, smh_132, sng0_95, \
                         sng0_96, sng1_95, sng1_96, snh_128, snh_131, \
                         snh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * smh_65[k]
                   + f_3 * pc_y[k] * snh_128[k];

        t_173[k] = f_17 * smh_131[k]
                   + f_4 * sng0_95[k]
                   - f_5 * sng1_95[k]
                   + f_3 * pc_x[k] * snh_131[k];

        t_174[k] = f_17 * smh_132[k]
                   + f_6 * sng0_96[k]
                   - f_7 * sng1_96[k]
                   + f_3 * pc_x[k] * snh_132[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pc_x, pc_y, pc_z, smh_68, smh_135, sng0_99, \
                         sng1_99, snh_129, snh_131, snh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * pc_z[k] * snh_129[k];

        t_176[k] = f_13 * smh_68[k]
                   + f_3 * pc_y[k] * snh_131[k];

        t_177[k] = f_17 * smh_135[k]
                   + f_6 * sng0_99[k]
                   - f_7 * sng1_99[k]
                   + f_3 * pc_x[k] * snh_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_z, smh_136, smh_138, sng0_100, \
                         sng0_102, sng1_100, sng1_102, snh_132, snh_136, \
                         snh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_17 * smh_136[k]
                   + f_8 * sng0_100[k]
                   - f_9 * sng1_100[k]
                   + f_3 * pc_x[k] * snh_136[k];

        t_179[k] = f_3 * pc_z[k] * snh_132[k];

        t_180[k] = f_17 * smh_138[k]
                   + f_8 * sng0_102[k]
                   - f_9 * sng1_102[k]
                   + f_3 * pc_x[k] * snh_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, smh_72, smh_140, smh_141, \
                         smh_142, sng0_104, sng1_104, snh_135, snh_140, snh_141, \
                         snh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * smh_72[k]
                   + f_3 * pc_y[k] * snh_135[k];

        t_182[k] = f_17 * smh_140[k]
                   + f_8 * sng0_104[k]
                   - f_9 * sng1_104[k]
                   + f_3 * pc_x[k] * snh_140[k];

        t_183[k] = f_17 * smh_141[k]
                   + f_3 * pc_x[k] * snh_141[k];

        t_184[k] = f_17 * smh_142[k]
                   + f_3 * pc_x[k] * snh_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, smh_143, smh_144, smh_145, smh_146, \
                         snh_143, snh_144, snh_145, snh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_17 * smh_143[k]
                   + f_3 * pc_x[k] * snh_143[k];

        t_186[k] = f_17 * smh_144[k]
                   + f_3 * pc_x[k] * snh_144[k];

        t_187[k] = f_17 * smh_145[k]
                   + f_3 * pc_x[k] * snh_145[k];

        t_188[k] = f_17 * smh_146[k]
                   + f_3 * pc_x[k] * snh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_y, pc_z, smh_78, smh_80, sng0_100, sng0_102, \
                         sng1_100, sng1_102, snh_141, snh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * smh_78[k]
                   + f_1 * sng0_100[k]
                   - f_2 * sng1_100[k]
                   + f_3 * pc_y[k] * snh_141[k];

        t_190[k] = f_3 * pc_z[k] * snh_141[k];

        t_191[k] = f_13 * smh_80[k]
                   + f_4 * sng0_102[k]
                   - f_5 * sng1_102[k]
                   + f_3 * pc_y[k] * snh_143[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, smh_81, smh_82, smh_83, \
                         sng0_103, sng0_104, sng1_103, sng1_104, snh_144, snh_145, \
                         snh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_13 * smh_81[k]
                   + f_6 * sng0_103[k]
                   - f_7 * sng1_103[k]
                   + f_3 * pc_y[k] * snh_144[k];

        t_193[k] = f_13 * smh_82[k]
                   + f_8 * sng0_104[k]
                   - f_9 * sng1_104[k]
                   + f_3 * pc_y[k] * snh_145[k];

        t_194[k] = f_13 * smh_83[k]
                   + f_3 * pc_y[k] * snh_146[k];

        t_195[k] = f_1 * sng0_104[k]
                   - f_2 * sng1_104[k]
                   + f_3 * pc_z[k] * snh_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_z, pc_y, pc_z, smi0_84, smi0_87, \
                         smh_63, smh_84, smi1_84, smi1_87, snh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * smi0_84[k]
                   - f_10 * pc_z[k] * smi1_84[k];

        t_197[k] = f_12 * smh_84[k]
                   + f_3 * pc_y[k] * snh_147[k];

        t_198[k] = f_11 * smh_63[k]
                   + f_3 * pc_z[k] * snh_147[k];

        t_199[k] = pb_z[k] * smi0_87[k]
                   - f_10 * pc_z[k] * smi1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_z, pc_x, pc_y, pc_z, smi0_90, smh_86, \
                         smh_152, smi1_90, sng0_110, sng1_110, snh_149, \
                         snh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * smh_86[k]
                   + f_3 * pc_y[k] * snh_149[k];

        t_201[k] = f_17 * smh_152[k]
                   + f_4 * sng0_110[k]
                   - f_5 * sng1_110[k]
                   + f_3 * pc_x[k] * snh_152[k];

        t_202[k] = pb_z[k] * smi0_90[k]
                   - f_10 * pc_z[k] * smi1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, pc_y, pc_z, smh_66, smh_89, smh_156, \
                         sng0_114, sng1_114, snh_150, snh_152, \
                         snh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * smh_66[k]
                   + f_3 * pc_z[k] * snh_150[k];

        t_204[k] = f_12 * smh_89[k]
                   + f_3 * pc_y[k] * snh_152[k];

        t_205[k] = f_17 * smh_156[k]
                   + f_6 * sng0_114[k]
                   - f_7 * sng1_114[k]
                   + f_3 * pc_x[k] * snh_156[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_z, pc_y, pc_z, smi0_94, smi0_96, \
                         smh_69, smh_70, smh_93, smi1_94, smi1_96, snh_153, \
                         snh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * smi0_94[k]
                   - f_10 * pc_z[k] * smi1_94[k];

        t_207[k] = f_11 * smh_69[k]
                   + f_3 * pc_z[k] * snh_153[k];

        t_208[k] = pb_z[k] * smi0_96[k]
                   + f_12 * smh_70[k]
                   - f_10 * pc_z[k] * smi1_96[k];

        t_209[k] = f_12 * smh_93[k]
                   + f_3 * pc_y[k] * snh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, smh_161, smh_162, smh_163, smh_164, \
                         sng0_119, sng1_119, snh_161, snh_162, snh_163, \
                         snh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * smh_161[k]
                   + f_8 * sng0_119[k]
                   - f_9 * sng1_119[k]
                   + f_3 * pc_x[k] * snh_161[k];

        t_211[k] = f_17 * smh_162[k]
                   + f_3 * pc_x[k] * snh_162[k];

        t_212[k] = f_17 * smh_163[k]
                   + f_3 * pc_x[k] * snh_163[k];

        t_213[k] = f_17 * smh_164[k]
                   + f_3 * pc_x[k] * snh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, smi0_105, smh_165, \
                         smh_166, smh_167, smi1_105, snh_165, snh_166, \
                         snh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_17 * smh_165[k]
                   + f_3 * pc_x[k] * snh_165[k];

        t_215[k] = f_17 * smh_166[k]
                   + f_3 * pc_x[k] * snh_166[k];

        t_216[k] = f_17 * smh_167[k]
                   + f_3 * pc_x[k] * snh_167[k];

        t_217[k] = pb_z[k] * smi0_105[k]
                   - f_10 * pc_z[k] * smi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, smh_78, smh_101, smh_102, sng0_117, \
                         sng0_118, sng1_117, sng1_118, snh_162, snh_164, \
                         snh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * smh_78[k]
                   + f_3 * pc_z[k] * snh_162[k];

        t_219[k] = f_12 * smh_101[k]
                   + f_4 * sng0_117[k]
                   - f_5 * sng1_117[k]
                   + f_3 * pc_y[k] * snh_164[k];

        t_220[k] = f_12 * smh_102[k]
                   + f_6 * sng0_118[k]
                   - f_7 * sng1_118[k]
                   + f_3 * pc_y[k] * snh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_y, pc_y, pc_z, smi0_140, smh_83, \
                         smh_103, smh_104, smi1_140, sng0_119, sng1_119, snh_166, \
                         snh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * smh_103[k]
                   + f_8 * sng0_119[k]
                   - f_9 * sng1_119[k]
                   + f_3 * pc_y[k] * snh_166[k];

        t_222[k] = f_12 * smh_104[k]
                   + f_3 * pc_y[k] * snh_167[k];

        t_223[k] = f_11 * smh_83[k]
                   + f_1 * sng0_119[k]
                   - f_2 * sng1_119[k]
                   + f_3 * pc_z[k] * snh_167[k];

        t_224[k] = pb_y[k] * smi0_140[k]
                   - f_10 * pc_y[k] * smi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_y, pc_y, pc_z, smi0_143, smh_84, \
                         smh_105, smh_106, smh_107, smi1_143, snh_168, \
                         snh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * smh_105[k]
                   + f_3 * pc_y[k] * snh_168[k];

        t_226[k] = f_12 * smh_84[k]
                   + f_3 * pc_z[k] * snh_168[k];

        t_227[k] = pb_y[k] * smi0_143[k]
                   + f_12 * smh_106[k]
                   - f_10 * pc_y[k] * smi1_143[k];

        t_228[k] = f_11 * smh_107[k]
                   + f_3 * pc_y[k] * snh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_y, pc_y, pc_z, smi0_145, smi0_146, \
                         smh_87, smh_108, smh_110, smi1_145, smi1_146, snh_171, \
                         snh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * smi0_145[k]
                   - f_10 * pc_y[k] * smi1_145[k];

        t_230[k] = pb_y[k] * smi0_146[k]
                   + f_13 * smh_108[k]
                   - f_10 * pc_y[k] * smi1_146[k];

        t_231[k] = f_12 * smh_87[k]
                   + f_3 * pc_z[k] * snh_171[k];

        t_232[k] = f_11 * smh_110[k]
                   + f_3 * pc_y[k] * snh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pc_y, pc_z, smi0_149, smi0_150, smh_90, \
                         smh_111, smi1_149, smi1_150, snh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * smi0_149[k]
                   - f_10 * pc_y[k] * smi1_149[k];

        t_234[k] = pb_y[k] * smi0_150[k]
                   + f_14 * smh_111[k]
                   - f_10 * pc_y[k] * smi1_150[k];

        t_235[k] = f_12 * smh_90[k]
                   + f_3 * pc_z[k] * snh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_y, pc_x, pc_y, smi0_152, smi0_154, \
                         smh_113, smh_114, smh_183, smi1_152, smi1_154, snh_177, \
                         snh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_y[k] * smi0_152[k]
                   + f_12 * smh_113[k]
                   - f_10 * pc_y[k] * smi1_152[k];

        t_237[k] = f_11 * smh_114[k]
                   + f_3 * pc_y[k] * snh_177[k];

        t_238[k] = pb_y[k] * smi0_154[k]
                   - f_10 * pc_y[k] * smi1_154[k];

        t_239[k] = f_17 * smh_183[k]
                   + f_3 * pc_x[k] * snh_183[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;

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

    const auto *smi0_167 = buffer.data(smi0 + 167);
    const auto *smi0_168 = buffer.data(smi0 + 168);
    const auto *smi0_171 = buffer.data(smi0 + 171);
    const auto *smi0_174 = buffer.data(smi0 + 174);
    const auto *smi0_178 = buffer.data(smi0 + 178);
    const auto *smi0_180 = buffer.data(smi0 + 180);
    const auto *smi0_189 = buffer.data(smi0 + 189);

    const auto *smh_99 = buffer.data(smh + 99);
    const auto *smh_105 = buffer.data(smh + 105);
    const auto *smh_108 = buffer.data(smh + 108);
    const auto *smh_111 = buffer.data(smh + 111);
    const auto *smh_120 = buffer.data(smh + 120);
    const auto *smh_122 = buffer.data(smh + 122);
    const auto *smh_123 = buffer.data(smh + 123);
    const auto *smh_124 = buffer.data(smh + 124);
    const auto *smh_125 = buffer.data(smh + 125);
    const auto *smh_126 = buffer.data(smh + 126);
    const auto *smh_128 = buffer.data(smh + 128);
    const auto *smh_129 = buffer.data(smh + 129);
    const auto *smh_131 = buffer.data(smh + 131);
    const auto *smh_132 = buffer.data(smh + 132);
    const auto *smh_133 = buffer.data(smh + 133);
    const auto *smh_135 = buffer.data(smh + 135);
    const auto *smh_141 = buffer.data(smh + 141);
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
    const auto *smh_164 = buffer.data(smh + 164);
    const auto *smh_165 = buffer.data(smh + 165);
    const auto *smh_166 = buffer.data(smh + 166);
    const auto *smh_167 = buffer.data(smh + 167);
    const auto *smh_168 = buffer.data(smh + 168);
    const auto *smh_170 = buffer.data(smh + 170);
    const auto *smh_173 = buffer.data(smh + 173);
    const auto *smh_177 = buffer.data(smh + 177);
    const auto *smh_184 = buffer.data(smh + 184);
    const auto *smh_185 = buffer.data(smh + 185);
    const auto *smh_186 = buffer.data(smh + 186);
    const auto *smh_187 = buffer.data(smh + 187);
    const auto *smh_188 = buffer.data(smh + 188);
    const auto *smh_189 = buffer.data(smh + 189);
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
    const auto *smh_236 = buffer.data(smh + 236);
    const auto *smh_240 = buffer.data(smh + 240);
    const auto *smh_245 = buffer.data(smh + 245);
    const auto *smh_246 = buffer.data(smh + 246);
    const auto *smh_247 = buffer.data(smh + 247);
    const auto *smh_248 = buffer.data(smh + 248);
    const auto *smh_249 = buffer.data(smh + 249);
    const auto *smh_250 = buffer.data(smh + 250);
    const auto *smh_251 = buffer.data(smh + 251);
    const auto *smh_252 = buffer.data(smh + 252);
    const auto *smh_255 = buffer.data(smh + 255);
    const auto *smh_257 = buffer.data(smh + 257);
    const auto *smh_258 = buffer.data(smh + 258);
    const auto *smh_261 = buffer.data(smh + 261);
    const auto *smh_262 = buffer.data(smh + 262);
    const auto *smh_264 = buffer.data(smh + 264);
    const auto *smh_266 = buffer.data(smh + 266);

    const auto *smi1_167 = buffer.data(smi1 + 167);
    const auto *smi1_168 = buffer.data(smi1 + 168);
    const auto *smi1_171 = buffer.data(smi1 + 171);
    const auto *smi1_174 = buffer.data(smi1 + 174);
    const auto *smi1_178 = buffer.data(smi1 + 178);
    const auto *smi1_180 = buffer.data(smi1 + 180);
    const auto *smi1_189 = buffer.data(smi1 + 189);

    const auto *sng0_130 = buffer.data(sng0 + 130);
    const auto *sng0_132 = buffer.data(sng0 + 132);
    const auto *sng0_133 = buffer.data(sng0 + 133);
    const auto *sng0_134 = buffer.data(sng0 + 134);
    const auto *sng0_135 = buffer.data(sng0 + 135);
    const auto *sng0_138 = buffer.data(sng0 + 138);
    const auto *sng0_140 = buffer.data(sng0 + 140);
    const auto *sng0_141 = buffer.data(sng0 + 141);
    const auto *sng0_144 = buffer.data(sng0 + 144);
    const auto *sng0_145 = buffer.data(sng0 + 145);
    const auto *sng0_147 = buffer.data(sng0 + 147);
    const auto *sng0_148 = buffer.data(sng0 + 148);
    const auto *sng0_149 = buffer.data(sng0 + 149);
    const auto *sng0_150 = buffer.data(sng0 + 150);
    const auto *sng0_153 = buffer.data(sng0 + 153);
    const auto *sng0_155 = buffer.data(sng0 + 155);
    const auto *sng0_156 = buffer.data(sng0 + 156);
    const auto *sng0_159 = buffer.data(sng0 + 159);
    const auto *sng0_160 = buffer.data(sng0 + 160);
    const auto *sng0_162 = buffer.data(sng0 + 162);
    const auto *sng0_163 = buffer.data(sng0 + 163);
    const auto *sng0_164 = buffer.data(sng0 + 164);
    const auto *sng0_170 = buffer.data(sng0 + 170);
    const auto *sng0_174 = buffer.data(sng0 + 174);
    const auto *sng0_177 = buffer.data(sng0 + 177);
    const auto *sng0_178 = buffer.data(sng0 + 178);
    const auto *sng0_179 = buffer.data(sng0 + 179);
    const auto *sng0_180 = buffer.data(sng0 + 180);
    const auto *sng0_183 = buffer.data(sng0 + 183);
    const auto *sng0_185 = buffer.data(sng0 + 185);
    const auto *sng0_186 = buffer.data(sng0 + 186);
    const auto *sng0_189 = buffer.data(sng0 + 189);
    const auto *sng0_190 = buffer.data(sng0 + 190);
    const auto *sng0_192 = buffer.data(sng0 + 192);
    const auto *sng0_194 = buffer.data(sng0 + 194);

    const auto *sng1_130 = buffer.data(sng1 + 130);
    const auto *sng1_132 = buffer.data(sng1 + 132);
    const auto *sng1_133 = buffer.data(sng1 + 133);
    const auto *sng1_134 = buffer.data(sng1 + 134);
    const auto *sng1_135 = buffer.data(sng1 + 135);
    const auto *sng1_138 = buffer.data(sng1 + 138);
    const auto *sng1_140 = buffer.data(sng1 + 140);
    const auto *sng1_141 = buffer.data(sng1 + 141);
    const auto *sng1_144 = buffer.data(sng1 + 144);
    const auto *sng1_145 = buffer.data(sng1 + 145);
    const auto *sng1_147 = buffer.data(sng1 + 147);
    const auto *sng1_148 = buffer.data(sng1 + 148);
    const auto *sng1_149 = buffer.data(sng1 + 149);
    const auto *sng1_150 = buffer.data(sng1 + 150);
    const auto *sng1_153 = buffer.data(sng1 + 153);
    const auto *sng1_155 = buffer.data(sng1 + 155);
    const auto *sng1_156 = buffer.data(sng1 + 156);
    const auto *sng1_159 = buffer.data(sng1 + 159);
    const auto *sng1_160 = buffer.data(sng1 + 160);
    const auto *sng1_162 = buffer.data(sng1 + 162);
    const auto *sng1_163 = buffer.data(sng1 + 163);
    const auto *sng1_164 = buffer.data(sng1 + 164);
    const auto *sng1_170 = buffer.data(sng1 + 170);
    const auto *sng1_174 = buffer.data(sng1 + 174);
    const auto *sng1_177 = buffer.data(sng1 + 177);
    const auto *sng1_178 = buffer.data(sng1 + 178);
    const auto *sng1_179 = buffer.data(sng1 + 179);
    const auto *sng1_180 = buffer.data(sng1 + 180);
    const auto *sng1_183 = buffer.data(sng1 + 183);
    const auto *sng1_185 = buffer.data(sng1 + 185);
    const auto *sng1_186 = buffer.data(sng1 + 186);
    const auto *sng1_189 = buffer.data(sng1 + 189);
    const auto *sng1_190 = buffer.data(sng1 + 190);
    const auto *sng1_192 = buffer.data(sng1 + 192);
    const auto *sng1_194 = buffer.data(sng1 + 194);

    const auto *snh_183 = buffer.data(snh + 183);
    const auto *snh_184 = buffer.data(snh + 184);
    const auto *snh_185 = buffer.data(snh + 185);
    const auto *snh_186 = buffer.data(snh + 186);
    const auto *snh_187 = buffer.data(snh + 187);
    const auto *snh_188 = buffer.data(snh + 188);
    const auto *snh_189 = buffer.data(snh + 189);
    const auto *snh_191 = buffer.data(snh + 191);
    const auto *snh_192 = buffer.data(snh + 192);
    const auto *snh_194 = buffer.data(snh + 194);
    const auto *snh_195 = buffer.data(snh + 195);
    const auto *snh_198 = buffer.data(snh + 198);
    const auto *snh_199 = buffer.data(snh + 199);
    const auto *snh_201 = buffer.data(snh + 201);
    const auto *snh_203 = buffer.data(snh + 203);
    const auto *snh_204 = buffer.data(snh + 204);
    const auto *snh_205 = buffer.data(snh + 205);
    const auto *snh_206 = buffer.data(snh + 206);
    const auto *snh_207 = buffer.data(snh + 207);
    const auto *snh_208 = buffer.data(snh + 208);
    const auto *snh_209 = buffer.data(snh + 209);
    const auto *snh_210 = buffer.data(snh + 210);
    const auto *snh_212 = buffer.data(snh + 212);
    const auto *snh_213 = buffer.data(snh + 213);
    const auto *snh_215 = buffer.data(snh + 215);
    const auto *snh_216 = buffer.data(snh + 216);
    const auto *snh_219 = buffer.data(snh + 219);
    const auto *snh_220 = buffer.data(snh + 220);
    const auto *snh_222 = buffer.data(snh + 222);
    const auto *snh_224 = buffer.data(snh + 224);
    const auto *snh_225 = buffer.data(snh + 225);
    const auto *snh_226 = buffer.data(snh + 226);
    const auto *snh_227 = buffer.data(snh + 227);
    const auto *snh_228 = buffer.data(snh + 228);
    const auto *snh_229 = buffer.data(snh + 229);
    const auto *snh_230 = buffer.data(snh + 230);
    const auto *snh_231 = buffer.data(snh + 231);
    const auto *snh_233 = buffer.data(snh + 233);
    const auto *snh_234 = buffer.data(snh + 234);
    const auto *snh_236 = buffer.data(snh + 236);
    const auto *snh_237 = buffer.data(snh + 237);
    const auto *snh_240 = buffer.data(snh + 240);
    const auto *snh_245 = buffer.data(snh + 245);
    const auto *snh_246 = buffer.data(snh + 246);
    const auto *snh_247 = buffer.data(snh + 247);
    const auto *snh_248 = buffer.data(snh + 248);
    const auto *snh_249 = buffer.data(snh + 249);
    const auto *snh_250 = buffer.data(snh + 250);
    const auto *snh_251 = buffer.data(snh + 251);
    const auto *snh_252 = buffer.data(snh + 252);
    const auto *snh_254 = buffer.data(snh + 254);
    const auto *snh_255 = buffer.data(snh + 255);
    const auto *snh_257 = buffer.data(snh + 257);
    const auto *snh_258 = buffer.data(snh + 258);
    const auto *snh_261 = buffer.data(snh + 261);
    const auto *snh_262 = buffer.data(snh + 262);
    const auto *snh_264 = buffer.data(snh + 264);
    const auto *snh_266 = buffer.data(snh + 266);

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, smh_184, smh_185, smh_186, \
                         smh_187, smh_188, snh_184, snh_185, snh_186, snh_187, \
                         snh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * smh_184[k]
                   + f_3 * pc_x[k] * snh_184[k];

        t_241[k] = f_17 * smh_185[k]
                   + f_3 * pc_x[k] * snh_185[k];

        t_242[k] = f_17 * smh_186[k]
                   + f_3 * pc_x[k] * snh_186[k];

        t_243[k] = f_17 * smh_187[k]
                   + f_3 * pc_x[k] * snh_187[k];

        t_244[k] = f_17 * smh_188[k]
                   + f_3 * pc_x[k] * snh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, smh_99, smh_120, smh_122, sng0_130, \
                         sng0_132, sng1_130, sng1_132, snh_183, \
                         snh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * smh_120[k]
                   + f_1 * sng0_130[k]
                   - f_2 * sng1_130[k]
                   + f_3 * pc_y[k] * snh_183[k];

        t_246[k] = f_12 * smh_99[k]
                   + f_3 * pc_z[k] * snh_183[k];

        t_247[k] = f_11 * smh_122[k]
                   + f_4 * sng0_132[k]
                   - f_5 * sng1_132[k]
                   + f_3 * pc_y[k] * snh_185[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, smh_123, smh_124, smh_125, sng0_133, \
                         sng0_134, sng1_133, sng1_134, snh_186, snh_187, \
                         snh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * smh_123[k]
                   + f_6 * sng0_133[k]
                   - f_7 * sng1_133[k]
                   + f_3 * pc_y[k] * snh_186[k];

        t_249[k] = f_11 * smh_124[k]
                   + f_8 * sng0_134[k]
                   - f_9 * sng1_134[k]
                   + f_3 * pc_y[k] * snh_187[k];

        t_250[k] = f_11 * smh_125[k]
                   + f_3 * pc_y[k] * snh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_x, pc_y, pc_z, smi0_167, \
                         smh_105, smh_189, smi1_167, sng0_135, sng1_135, \
                         snh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * smi0_167[k]
                   - f_10 * pc_y[k] * smi1_167[k];

        t_252[k] = f_17 * smh_189[k]
                   + f_1 * sng0_135[k]
                   - f_2 * sng1_135[k]
                   + f_3 * pc_x[k] * snh_189[k];

        t_253[k] = f_3 * pc_y[k] * snh_189[k];

        t_254[k] = f_13 * smh_105[k]
                   + f_3 * pc_z[k] * snh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, smh_192, smh_194, sng0_138, \
                         sng0_140, sng1_138, sng1_140, snh_191, snh_192, \
                         snh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_17 * smh_192[k]
                   + f_4 * sng0_138[k]
                   - f_5 * sng1_138[k]
                   + f_3 * pc_x[k] * snh_192[k];

        t_256[k] = f_3 * pc_y[k] * snh_191[k];

        t_257[k] = f_17 * smh_194[k]
                   + f_4 * sng0_140[k]
                   - f_5 * sng1_140[k]
                   + f_3 * pc_x[k] * snh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, pc_z, smh_108, smh_195, sng0_141, \
                         sng1_141, snh_192, snh_194, snh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_17 * smh_195[k]
                   + f_6 * sng0_141[k]
                   - f_7 * sng1_141[k]
                   + f_3 * pc_x[k] * snh_195[k];

        t_259[k] = f_13 * smh_108[k]
                   + f_3 * pc_z[k] * snh_192[k];

        t_260[k] = f_3 * pc_y[k] * snh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, smh_111, smh_198, smh_199, sng0_144, \
                         sng0_145, sng1_144, sng1_145, snh_195, snh_198, \
                         snh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_17 * smh_198[k]
                   + f_6 * sng0_144[k]
                   - f_7 * sng1_144[k]
                   + f_3 * pc_x[k] * snh_198[k];

        t_262[k] = f_17 * smh_199[k]
                   + f_8 * sng0_145[k]
                   - f_9 * sng1_145[k]
                   + f_3 * pc_x[k] * snh_199[k];

        t_263[k] = f_13 * smh_111[k]
                   + f_3 * pc_z[k] * snh_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, smh_201, smh_203, sng0_147, \
                         sng0_149, sng1_147, sng1_149, snh_198, snh_201, \
                         snh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * smh_201[k]
                   + f_8 * sng0_147[k]
                   - f_9 * sng1_147[k]
                   + f_3 * pc_x[k] * snh_201[k];

        t_265[k] = f_3 * pc_y[k] * snh_198[k];

        t_266[k] = f_17 * smh_203[k]
                   + f_8 * sng0_149[k]
                   - f_9 * sng1_149[k]
                   + f_3 * pc_x[k] * snh_203[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pc_x, smh_204, smh_205, smh_206, \
                         smh_207, smh_208, snh_204, snh_205, snh_206, snh_207, \
                         snh_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_17 * smh_204[k]
                   + f_3 * pc_x[k] * snh_204[k];

        t_268[k] = f_17 * smh_205[k]
                   + f_3 * pc_x[k] * snh_205[k];

        t_269[k] = f_17 * smh_206[k]
                   + f_3 * pc_x[k] * snh_206[k];

        t_270[k] = f_17 * smh_207[k]
                   + f_3 * pc_x[k] * snh_207[k];

        t_271[k] = f_17 * smh_208[k]
                   + f_3 * pc_x[k] * snh_208[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, smh_120, smh_209, \
                         sng0_145, sng0_147, sng1_145, sng1_147, snh_204, snh_206, \
                         snh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_17 * smh_209[k]
                   + f_3 * pc_x[k] * snh_209[k];

        t_273[k] = f_1 * sng0_145[k]
                   - f_2 * sng1_145[k]
                   + f_3 * pc_y[k] * snh_204[k];

        t_274[k] = f_13 * smh_120[k]
                   + f_3 * pc_z[k] * snh_204[k];

        t_275[k] = f_4 * sng0_147[k]
                   - f_5 * sng1_147[k]
                   + f_3 * pc_y[k] * snh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, smh_125, sng0_148, sng0_149, \
                         sng1_148, sng1_149, snh_207, snh_208, \
                         snh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * sng0_148[k]
                   - f_7 * sng1_148[k]
                   + f_3 * pc_y[k] * snh_207[k];

        t_277[k] = f_8 * sng0_149[k]
                   - f_9 * sng1_149[k]
                   + f_3 * pc_y[k] * snh_208[k];

        t_278[k] = f_3 * pc_y[k] * snh_209[k];

        t_279[k] = f_13 * smh_125[k]
                   + f_1 * sng0_149[k]
                   - f_2 * sng1_149[k]
                   + f_3 * pc_z[k] * snh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, smh_126, smh_210, \
                         smh_213, sng0_150, sng0_153, sng1_150, sng1_153, snh_210, \
                         snh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_18 * smh_210[k]
                   + f_1 * sng0_150[k]
                   - f_2 * sng1_150[k]
                   + f_3 * pc_x[k] * snh_210[k];

        t_281[k] = f_14 * smh_126[k]
                   + f_3 * pc_y[k] * snh_210[k];

        t_282[k] = f_3 * pc_z[k] * snh_210[k];

        t_283[k] = f_18 * smh_213[k]
                   + f_4 * sng0_153[k]
                   - f_5 * sng1_153[k]
                   + f_3 * pc_x[k] * snh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_y, smh_128, smh_215, smh_216, sng0_155, \
                         sng0_156, sng1_155, sng1_156, snh_212, snh_215, \
                         snh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_14 * smh_128[k]
                   + f_3 * pc_y[k] * snh_212[k];

        t_285[k] = f_18 * smh_215[k]
                   + f_4 * sng0_155[k]
                   - f_5 * sng1_155[k]
                   + f_3 * pc_x[k] * snh_215[k];

        t_286[k] = f_18 * smh_216[k]
                   + f_6 * sng0_156[k]
                   - f_7 * sng1_156[k]
                   + f_3 * pc_x[k] * snh_216[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, pc_z, smh_131, smh_219, sng0_159, \
                         sng1_159, snh_213, snh_215, snh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * snh_213[k];

        t_288[k] = f_14 * smh_131[k]
                   + f_3 * pc_y[k] * snh_215[k];

        t_289[k] = f_18 * smh_219[k]
                   + f_6 * sng0_159[k]
                   - f_7 * sng1_159[k]
                   + f_3 * pc_x[k] * snh_219[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, smh_220, smh_222, sng0_160, \
                         sng0_162, sng1_160, sng1_162, snh_216, snh_220, \
                         snh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_18 * smh_220[k]
                   + f_8 * sng0_160[k]
                   - f_9 * sng1_160[k]
                   + f_3 * pc_x[k] * snh_220[k];

        t_291[k] = f_3 * pc_z[k] * snh_216[k];

        t_292[k] = f_18 * smh_222[k]
                   + f_8 * sng0_162[k]
                   - f_9 * sng1_162[k]
                   + f_3 * pc_x[k] * snh_222[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pc_x, pc_y, smh_135, smh_224, smh_225, \
                         smh_226, sng0_164, sng1_164, snh_219, snh_224, snh_225, \
                         snh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_14 * smh_135[k]
                   + f_3 * pc_y[k] * snh_219[k];

        t_294[k] = f_18 * smh_224[k]
                   + f_8 * sng0_164[k]
                   - f_9 * sng1_164[k]
                   + f_3 * pc_x[k] * snh_224[k];

        t_295[k] = f_18 * smh_225[k]
                   + f_3 * pc_x[k] * snh_225[k];

        t_296[k] = f_18 * smh_226[k]
                   + f_3 * pc_x[k] * snh_226[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pc_x, smh_227, smh_228, smh_229, smh_230, \
                         snh_227, snh_228, snh_229, snh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_18 * smh_227[k]
                   + f_3 * pc_x[k] * snh_227[k];

        t_298[k] = f_18 * smh_228[k]
                   + f_3 * pc_x[k] * snh_228[k];

        t_299[k] = f_18 * smh_229[k]
                   + f_3 * pc_x[k] * snh_229[k];

        t_300[k] = f_18 * smh_230[k]
                   + f_3 * pc_x[k] * snh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, pc_y, pc_z, smh_141, smh_143, sng0_160, \
                         sng0_162, sng1_160, sng1_162, snh_225, \
                         snh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * smh_141[k]
                   + f_1 * sng0_160[k]
                   - f_2 * sng1_160[k]
                   + f_3 * pc_y[k] * snh_225[k];

        t_302[k] = f_3 * pc_z[k] * snh_225[k];

        t_303[k] = f_14 * smh_143[k]
                   + f_4 * sng0_162[k]
                   - f_5 * sng1_162[k]
                   + f_3 * pc_y[k] * snh_227[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_y, pc_z, smh_144, smh_145, smh_146, \
                         sng0_163, sng0_164, sng1_163, sng1_164, snh_228, snh_229, \
                         snh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_14 * smh_144[k]
                   + f_6 * sng0_163[k]
                   - f_7 * sng1_163[k]
                   + f_3 * pc_y[k] * snh_228[k];

        t_305[k] = f_14 * smh_145[k]
                   + f_8 * sng0_164[k]
                   - f_9 * sng1_164[k]
                   + f_3 * pc_y[k] * snh_229[k];

        t_306[k] = f_14 * smh_146[k]
                   + f_3 * pc_y[k] * snh_230[k];

        t_307[k] = f_1 * sng0_164[k]
                   - f_2 * sng1_164[k]
                   + f_3 * pc_z[k] * snh_230[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, smi0_168, smi0_171, \
                         smh_126, smh_147, smi1_168, smi1_171, \
                         snh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * smi0_168[k]
                   - f_10 * pc_z[k] * smi1_168[k];

        t_309[k] = f_13 * smh_147[k]
                   + f_3 * pc_y[k] * snh_231[k];

        t_310[k] = f_11 * smh_126[k]
                   + f_3 * pc_z[k] * snh_231[k];

        t_311[k] = pb_z[k] * smi0_171[k]
                   - f_10 * pc_z[k] * smi1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_z, pc_x, pc_y, pc_z, smi0_174, smh_149, \
                         smh_236, smi1_174, sng0_170, sng1_170, snh_233, \
                         snh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * smh_149[k]
                   + f_3 * pc_y[k] * snh_233[k];

        t_313[k] = f_18 * smh_236[k]
                   + f_4 * sng0_170[k]
                   - f_5 * sng1_170[k]
                   + f_3 * pc_x[k] * snh_236[k];

        t_314[k] = pb_z[k] * smi0_174[k]
                   - f_10 * pc_z[k] * smi1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pc_x, pc_y, pc_z, smh_129, smh_152, smh_240, \
                         sng0_174, sng1_174, snh_234, snh_236, \
                         snh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * smh_129[k]
                   + f_3 * pc_z[k] * snh_234[k];

        t_316[k] = f_13 * smh_152[k]
                   + f_3 * pc_y[k] * snh_236[k];

        t_317[k] = f_18 * smh_240[k]
                   + f_6 * sng0_174[k]
                   - f_7 * sng1_174[k]
                   + f_3 * pc_x[k] * snh_240[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_z, pc_y, pc_z, smi0_178, smi0_180, \
                         smh_132, smh_133, smh_156, smi1_178, smi1_180, snh_237, \
                         snh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * smi0_178[k]
                   - f_10 * pc_z[k] * smi1_178[k];

        t_319[k] = f_11 * smh_132[k]
                   + f_3 * pc_z[k] * snh_237[k];

        t_320[k] = pb_z[k] * smi0_180[k]
                   + f_12 * smh_133[k]
                   - f_10 * pc_z[k] * smi1_180[k];

        t_321[k] = f_13 * smh_156[k]
                   + f_3 * pc_y[k] * snh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, smh_245, smh_246, smh_247, smh_248, \
                         sng0_179, sng1_179, snh_245, snh_246, snh_247, \
                         snh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_18 * smh_245[k]
                   + f_8 * sng0_179[k]
                   - f_9 * sng1_179[k]
                   + f_3 * pc_x[k] * snh_245[k];

        t_323[k] = f_18 * smh_246[k]
                   + f_3 * pc_x[k] * snh_246[k];

        t_324[k] = f_18 * smh_247[k]
                   + f_3 * pc_x[k] * snh_247[k];

        t_325[k] = f_18 * smh_248[k]
                   + f_3 * pc_x[k] * snh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pb_z, pc_x, pc_z, smi0_189, smh_249, \
                         smh_250, smh_251, smi1_189, snh_249, snh_250, \
                         snh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_18 * smh_249[k]
                   + f_3 * pc_x[k] * snh_249[k];

        t_327[k] = f_18 * smh_250[k]
                   + f_3 * pc_x[k] * snh_250[k];

        t_328[k] = f_18 * smh_251[k]
                   + f_3 * pc_x[k] * snh_251[k];

        t_329[k] = pb_z[k] * smi0_189[k]
                   - f_10 * pc_z[k] * smi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, smh_141, smh_164, smh_165, sng0_177, \
                         sng0_178, sng1_177, sng1_178, snh_246, snh_248, \
                         snh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * smh_141[k]
                   + f_3 * pc_z[k] * snh_246[k];

        t_331[k] = f_13 * smh_164[k]
                   + f_4 * sng0_177[k]
                   - f_5 * sng1_177[k]
                   + f_3 * pc_y[k] * snh_248[k];

        t_332[k] = f_13 * smh_165[k]
                   + f_6 * sng0_178[k]
                   - f_7 * sng1_178[k]
                   + f_3 * pc_y[k] * snh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, smh_146, smh_166, smh_167, sng0_179, \
                         sng1_179, snh_250, snh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * smh_166[k]
                   + f_8 * sng0_179[k]
                   - f_9 * sng1_179[k]
                   + f_3 * pc_y[k] * snh_250[k];

        t_334[k] = f_13 * smh_167[k]
                   + f_3 * pc_y[k] * snh_251[k];

        t_335[k] = f_11 * smh_146[k]
                   + f_1 * sng0_179[k]
                   - f_2 * sng1_179[k]
                   + f_3 * pc_z[k] * snh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, smh_147, smh_168, smh_252, \
                         sng0_180, sng1_180, snh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_18 * smh_252[k]
                   + f_1 * sng0_180[k]
                   - f_2 * sng1_180[k]
                   + f_3 * pc_x[k] * snh_252[k];

        t_337[k] = f_12 * smh_168[k]
                   + f_3 * pc_y[k] * snh_252[k];

        t_338[k] = f_12 * smh_147[k]
                   + f_3 * pc_z[k] * snh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, smh_170, smh_255, smh_257, sng0_183, \
                         sng0_185, sng1_183, sng1_185, snh_254, snh_255, \
                         snh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_18 * smh_255[k]
                   + f_4 * sng0_183[k]
                   - f_5 * sng1_183[k]
                   + f_3 * pc_x[k] * snh_255[k];

        t_340[k] = f_12 * smh_170[k]
                   + f_3 * pc_y[k] * snh_254[k];

        t_341[k] = f_18 * smh_257[k]
                   + f_4 * sng0_185[k]
                   - f_5 * sng1_185[k]
                   + f_3 * pc_x[k] * snh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, smh_150, smh_173, smh_258, \
                         sng0_186, sng1_186, snh_255, snh_257, \
                         snh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_18 * smh_258[k]
                   + f_6 * sng0_186[k]
                   - f_7 * sng1_186[k]
                   + f_3 * pc_x[k] * snh_258[k];

        t_343[k] = f_12 * smh_150[k]
                   + f_3 * pc_z[k] * snh_255[k];

        t_344[k] = f_12 * smh_173[k]
                   + f_3 * pc_y[k] * snh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, smh_153, smh_261, smh_262, sng0_189, \
                         sng0_190, sng1_189, sng1_190, snh_258, snh_261, \
                         snh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_18 * smh_261[k]
                   + f_6 * sng0_189[k]
                   - f_7 * sng1_189[k]
                   + f_3 * pc_x[k] * snh_261[k];

        t_346[k] = f_18 * smh_262[k]
                   + f_8 * sng0_190[k]
                   - f_9 * sng1_190[k]
                   + f_3 * pc_x[k] * snh_262[k];

        t_347[k] = f_12 * smh_153[k]
                   + f_3 * pc_z[k] * snh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, smh_177, smh_264, smh_266, sng0_192, \
                         sng0_194, sng1_192, sng1_194, snh_261, snh_264, \
                         snh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_18 * smh_264[k]
                   + f_8 * sng0_192[k]
                   - f_9 * sng1_192[k]
                   + f_3 * pc_x[k] * snh_264[k];

        t_349[k] = f_12 * smh_177[k]
                   + f_3 * pc_y[k] * snh_261[k];

        t_350[k] = f_18 * smh_266[k]
                   + f_8 * sng0_194[k]
                   - f_9 * sng1_194[k]
                   + f_3 * pc_x[k] * snh_266[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *smi0_252 = buffer.data(smi0 + 252);
    const auto *smi0_255 = buffer.data(smi0 + 255);
    const auto *smi0_257 = buffer.data(smi0 + 257);
    const auto *smi0_258 = buffer.data(smi0 + 258);
    const auto *smi0_261 = buffer.data(smi0 + 261);
    const auto *smi0_262 = buffer.data(smi0 + 262);
    const auto *smi0_264 = buffer.data(smi0 + 264);
    const auto *smi0_266 = buffer.data(smi0 + 266);
    const auto *smi0_279 = buffer.data(smi0 + 279);
    const auto *smi0_280 = buffer.data(smi0 + 280);
    const auto *smi0_283 = buffer.data(smi0 + 283);
    const auto *smi0_286 = buffer.data(smi0 + 286);
    const auto *smi0_290 = buffer.data(smi0 + 290);
    const auto *smi0_292 = buffer.data(smi0 + 292);

    const auto *smh_162 = buffer.data(smh + 162);
    const auto *smh_167 = buffer.data(smh + 167);
    const auto *smh_168 = buffer.data(smh + 168);
    const auto *smh_171 = buffer.data(smh + 171);
    const auto *smh_174 = buffer.data(smh + 174);
    const auto *smh_183 = buffer.data(smh + 183);
    const auto *smh_185 = buffer.data(smh + 185);
    const auto *smh_186 = buffer.data(smh + 186);
    const auto *smh_187 = buffer.data(smh + 187);
    const auto *smh_188 = buffer.data(smh + 188);
    const auto *smh_189 = buffer.data(smh + 189);
    const auto *smh_190 = buffer.data(smh + 190);
    const auto *smh_191 = buffer.data(smh + 191);
    const auto *smh_192 = buffer.data(smh + 192);
    const auto *smh_194 = buffer.data(smh + 194);
    const auto *smh_195 = buffer.data(smh + 195);
    const auto *smh_197 = buffer.data(smh + 197);
    const auto *smh_198 = buffer.data(smh + 198);
    const auto *smh_204 = buffer.data(smh + 204);
    const auto *smh_206 = buffer.data(smh + 206);
    const auto *smh_207 = buffer.data(smh + 207);
    const auto *smh_208 = buffer.data(smh + 208);
    const auto *smh_209 = buffer.data(smh + 209);
    const auto *smh_210 = buffer.data(smh + 210);
    const auto *smh_212 = buffer.data(smh + 212);
    const auto *smh_213 = buffer.data(smh + 213);
    const auto *smh_215 = buffer.data(smh + 215);
    const auto *smh_216 = buffer.data(smh + 216);
    const auto *smh_217 = buffer.data(smh + 217);
    const auto *smh_219 = buffer.data(smh + 219);
    const auto *smh_225 = buffer.data(smh + 225);
    const auto *smh_227 = buffer.data(smh + 227);
    const auto *smh_228 = buffer.data(smh + 228);
    const auto *smh_229 = buffer.data(smh + 229);
    const auto *smh_230 = buffer.data(smh + 230);
    const auto *smh_231 = buffer.data(smh + 231);
    const auto *smh_233 = buffer.data(smh + 233);
    const auto *smh_236 = buffer.data(smh + 236);
    const auto *smh_240 = buffer.data(smh + 240);
    const auto *smh_267 = buffer.data(smh + 267);
    const auto *smh_268 = buffer.data(smh + 268);
    const auto *smh_269 = buffer.data(smh + 269);
    const auto *smh_270 = buffer.data(smh + 270);
    const auto *smh_271 = buffer.data(smh + 271);
    const auto *smh_272 = buffer.data(smh + 272);
    const auto *smh_288 = buffer.data(smh + 288);
    const auto *smh_289 = buffer.data(smh + 289);
    const auto *smh_290 = buffer.data(smh + 290);
    const auto *smh_291 = buffer.data(smh + 291);
    const auto *smh_292 = buffer.data(smh + 292);
    const auto *smh_293 = buffer.data(smh + 293);
    const auto *smh_294 = buffer.data(smh + 294);
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
    const auto *smh_341 = buffer.data(smh + 341);
    const auto *smh_345 = buffer.data(smh + 345);
    const auto *smh_350 = buffer.data(smh + 350);
    const auto *smh_351 = buffer.data(smh + 351);
    const auto *smh_352 = buffer.data(smh + 352);
    const auto *smh_353 = buffer.data(smh + 353);

    const auto *smi1_252 = buffer.data(smi1 + 252);
    const auto *smi1_255 = buffer.data(smi1 + 255);
    const auto *smi1_257 = buffer.data(smi1 + 257);
    const auto *smi1_258 = buffer.data(smi1 + 258);
    const auto *smi1_261 = buffer.data(smi1 + 261);
    const auto *smi1_262 = buffer.data(smi1 + 262);
    const auto *smi1_264 = buffer.data(smi1 + 264);
    const auto *smi1_266 = buffer.data(smi1 + 266);
    const auto *smi1_279 = buffer.data(smi1 + 279);
    const auto *smi1_280 = buffer.data(smi1 + 280);
    const auto *smi1_283 = buffer.data(smi1 + 283);
    const auto *smi1_286 = buffer.data(smi1 + 286);
    const auto *smi1_290 = buffer.data(smi1 + 290);
    const auto *smi1_292 = buffer.data(smi1 + 292);

    const auto *sng0_190 = buffer.data(sng0 + 190);
    const auto *sng0_192 = buffer.data(sng0 + 192);
    const auto *sng0_193 = buffer.data(sng0 + 193);
    const auto *sng0_194 = buffer.data(sng0 + 194);
    const auto *sng0_205 = buffer.data(sng0 + 205);
    const auto *sng0_207 = buffer.data(sng0 + 207);
    const auto *sng0_208 = buffer.data(sng0 + 208);
    const auto *sng0_209 = buffer.data(sng0 + 209);
    const auto *sng0_210 = buffer.data(sng0 + 210);
    const auto *sng0_213 = buffer.data(sng0 + 213);
    const auto *sng0_215 = buffer.data(sng0 + 215);
    const auto *sng0_216 = buffer.data(sng0 + 216);
    const auto *sng0_219 = buffer.data(sng0 + 219);
    const auto *sng0_220 = buffer.data(sng0 + 220);
    const auto *sng0_222 = buffer.data(sng0 + 222);
    const auto *sng0_223 = buffer.data(sng0 + 223);
    const auto *sng0_224 = buffer.data(sng0 + 224);
    const auto *sng0_225 = buffer.data(sng0 + 225);
    const auto *sng0_228 = buffer.data(sng0 + 228);
    const auto *sng0_230 = buffer.data(sng0 + 230);
    const auto *sng0_231 = buffer.data(sng0 + 231);
    const auto *sng0_234 = buffer.data(sng0 + 234);
    const auto *sng0_235 = buffer.data(sng0 + 235);
    const auto *sng0_237 = buffer.data(sng0 + 237);
    const auto *sng0_238 = buffer.data(sng0 + 238);
    const auto *sng0_239 = buffer.data(sng0 + 239);
    const auto *sng0_245 = buffer.data(sng0 + 245);
    const auto *sng0_249 = buffer.data(sng0 + 249);
    const auto *sng0_254 = buffer.data(sng0 + 254);

    const auto *sng1_190 = buffer.data(sng1 + 190);
    const auto *sng1_192 = buffer.data(sng1 + 192);
    const auto *sng1_193 = buffer.data(sng1 + 193);
    const auto *sng1_194 = buffer.data(sng1 + 194);
    const auto *sng1_205 = buffer.data(sng1 + 205);
    const auto *sng1_207 = buffer.data(sng1 + 207);
    const auto *sng1_208 = buffer.data(sng1 + 208);
    const auto *sng1_209 = buffer.data(sng1 + 209);
    const auto *sng1_210 = buffer.data(sng1 + 210);
    const auto *sng1_213 = buffer.data(sng1 + 213);
    const auto *sng1_215 = buffer.data(sng1 + 215);
    const auto *sng1_216 = buffer.data(sng1 + 216);
    const auto *sng1_219 = buffer.data(sng1 + 219);
    const auto *sng1_220 = buffer.data(sng1 + 220);
    const auto *sng1_222 = buffer.data(sng1 + 222);
    const auto *sng1_223 = buffer.data(sng1 + 223);
    const auto *sng1_224 = buffer.data(sng1 + 224);
    const auto *sng1_225 = buffer.data(sng1 + 225);
    const auto *sng1_228 = buffer.data(sng1 + 228);
    const auto *sng1_230 = buffer.data(sng1 + 230);
    const auto *sng1_231 = buffer.data(sng1 + 231);
    const auto *sng1_234 = buffer.data(sng1 + 234);
    const auto *sng1_235 = buffer.data(sng1 + 235);
    const auto *sng1_237 = buffer.data(sng1 + 237);
    const auto *sng1_238 = buffer.data(sng1 + 238);
    const auto *sng1_239 = buffer.data(sng1 + 239);
    const auto *sng1_245 = buffer.data(sng1 + 245);
    const auto *sng1_249 = buffer.data(sng1 + 249);
    const auto *sng1_254 = buffer.data(sng1 + 254);

    const auto *snh_267 = buffer.data(snh + 267);
    const auto *snh_268 = buffer.data(snh + 268);
    const auto *snh_269 = buffer.data(snh + 269);
    const auto *snh_270 = buffer.data(snh + 270);
    const auto *snh_271 = buffer.data(snh + 271);
    const auto *snh_272 = buffer.data(snh + 272);
    const auto *snh_273 = buffer.data(snh + 273);
    const auto *snh_275 = buffer.data(snh + 275);
    const auto *snh_276 = buffer.data(snh + 276);
    const auto *snh_278 = buffer.data(snh + 278);
    const auto *snh_279 = buffer.data(snh + 279);
    const auto *snh_282 = buffer.data(snh + 282);
    const auto *snh_288 = buffer.data(snh + 288);
    const auto *snh_289 = buffer.data(snh + 289);
    const auto *snh_290 = buffer.data(snh + 290);
    const auto *snh_291 = buffer.data(snh + 291);
    const auto *snh_292 = buffer.data(snh + 292);
    const auto *snh_293 = buffer.data(snh + 293);
    const auto *snh_294 = buffer.data(snh + 294);
    const auto *snh_296 = buffer.data(snh + 296);
    const auto *snh_297 = buffer.data(snh + 297);
    const auto *snh_299 = buffer.data(snh + 299);
    const auto *snh_300 = buffer.data(snh + 300);
    const auto *snh_303 = buffer.data(snh + 303);
    const auto *snh_304 = buffer.data(snh + 304);
    const auto *snh_306 = buffer.data(snh + 306);
    const auto *snh_308 = buffer.data(snh + 308);
    const auto *snh_309 = buffer.data(snh + 309);
    const auto *snh_310 = buffer.data(snh + 310);
    const auto *snh_311 = buffer.data(snh + 311);
    const auto *snh_312 = buffer.data(snh + 312);
    const auto *snh_313 = buffer.data(snh + 313);
    const auto *snh_314 = buffer.data(snh + 314);
    const auto *snh_315 = buffer.data(snh + 315);
    const auto *snh_317 = buffer.data(snh + 317);
    const auto *snh_318 = buffer.data(snh + 318);
    const auto *snh_320 = buffer.data(snh + 320);
    const auto *snh_321 = buffer.data(snh + 321);
    const auto *snh_324 = buffer.data(snh + 324);
    const auto *snh_325 = buffer.data(snh + 325);
    const auto *snh_327 = buffer.data(snh + 327);
    const auto *snh_329 = buffer.data(snh + 329);
    const auto *snh_330 = buffer.data(snh + 330);
    const auto *snh_331 = buffer.data(snh + 331);
    const auto *snh_332 = buffer.data(snh + 332);
    const auto *snh_333 = buffer.data(snh + 333);
    const auto *snh_334 = buffer.data(snh + 334);
    const auto *snh_335 = buffer.data(snh + 335);
    const auto *snh_336 = buffer.data(snh + 336);
    const auto *snh_338 = buffer.data(snh + 338);
    const auto *snh_339 = buffer.data(snh + 339);
    const auto *snh_341 = buffer.data(snh + 341);
    const auto *snh_342 = buffer.data(snh + 342);
    const auto *snh_345 = buffer.data(snh + 345);
    const auto *snh_350 = buffer.data(snh + 350);
    const auto *snh_351 = buffer.data(snh + 351);
    const auto *snh_352 = buffer.data(snh + 352);
    const auto *snh_353 = buffer.data(snh + 353);

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, smh_267, smh_268, smh_269, \
                         smh_270, smh_271, snh_267, snh_268, snh_269, snh_270, \
                         snh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_18 * smh_267[k]
                   + f_3 * pc_x[k] * snh_267[k];

        t_352[k] = f_18 * smh_268[k]
                   + f_3 * pc_x[k] * snh_268[k];

        t_353[k] = f_18 * smh_269[k]
                   + f_3 * pc_x[k] * snh_269[k];

        t_354[k] = f_18 * smh_270[k]
                   + f_3 * pc_x[k] * snh_270[k];

        t_355[k] = f_18 * smh_271[k]
                   + f_3 * pc_x[k] * snh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, smh_162, smh_183, smh_272, \
                         sng0_190, sng1_190, snh_267, snh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_18 * smh_272[k]
                   + f_3 * pc_x[k] * snh_272[k];

        t_357[k] = f_12 * smh_183[k]
                   + f_1 * sng0_190[k]
                   - f_2 * sng1_190[k]
                   + f_3 * pc_y[k] * snh_267[k];

        t_358[k] = f_12 * smh_162[k]
                   + f_3 * pc_z[k] * snh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, smh_185, smh_186, smh_187, sng0_192, \
                         sng0_193, sng0_194, sng1_192, sng1_193, sng1_194, snh_269, snh_270, \
                         snh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * smh_185[k]
                   + f_4 * sng0_192[k]
                   - f_5 * sng1_192[k]
                   + f_3 * pc_y[k] * snh_269[k];

        t_360[k] = f_12 * smh_186[k]
                   + f_6 * sng0_193[k]
                   - f_7 * sng1_193[k]
                   + f_3 * pc_y[k] * snh_270[k];

        t_361[k] = f_12 * smh_187[k]
                   + f_8 * sng0_194[k]
                   - f_9 * sng1_194[k]
                   + f_3 * pc_y[k] * snh_271[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pb_y, pc_y, pc_z, smi0_252, smh_167, \
                         smh_188, smh_189, smi1_252, sng0_194, sng1_194, snh_272, \
                         snh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * smh_188[k]
                   + f_3 * pc_y[k] * snh_272[k];

        t_363[k] = f_12 * smh_167[k]
                   + f_1 * sng0_194[k]
                   - f_2 * sng1_194[k]
                   + f_3 * pc_z[k] * snh_272[k];

        t_364[k] = pb_y[k] * smi0_252[k]
                   - f_10 * pc_y[k] * smi1_252[k];

        t_365[k] = f_11 * smh_189[k]
                   + f_3 * pc_y[k] * snh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_y, pc_y, pc_z, smi0_255, smi0_257, \
                         smh_168, smh_190, smh_191, smi1_255, smi1_257, snh_273, \
                         snh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * smh_168[k]
                   + f_3 * pc_z[k] * snh_273[k];

        t_367[k] = pb_y[k] * smi0_255[k]
                   + f_12 * smh_190[k]
                   - f_10 * pc_y[k] * smi1_255[k];

        t_368[k] = f_11 * smh_191[k]
                   + f_3 * pc_y[k] * snh_275[k];

        t_369[k] = pb_y[k] * smi0_257[k]
                   - f_10 * pc_y[k] * smi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_y, pc_y, pc_z, smi0_258, smi0_261, \
                         smh_171, smh_192, smh_194, smi1_258, smi1_261, snh_276, \
                         snh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_y[k] * smi0_258[k]
                   + f_13 * smh_192[k]
                   - f_10 * pc_y[k] * smi1_258[k];

        t_371[k] = f_13 * smh_171[k]
                   + f_3 * pc_z[k] * snh_276[k];

        t_372[k] = f_11 * smh_194[k]
                   + f_3 * pc_y[k] * snh_278[k];

        t_373[k] = pb_y[k] * smi0_261[k]
                   - f_10 * pc_y[k] * smi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pb_y, pc_y, pc_z, smi0_262, smi0_264, smh_174, \
                         smh_195, smh_197, smi1_262, smi1_264, \
                         snh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_y[k] * smi0_262[k]
                   + f_14 * smh_195[k]
                   - f_10 * pc_y[k] * smi1_262[k];

        t_375[k] = f_13 * smh_174[k]
                   + f_3 * pc_z[k] * snh_279[k];

        t_376[k] = pb_y[k] * smi0_264[k]
                   + f_12 * smh_197[k]
                   - f_10 * pc_y[k] * smi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_y, pc_x, pc_y, smi0_266, smh_198, \
                         smh_288, smh_289, smi1_266, snh_282, snh_288, \
                         snh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * smh_198[k]
                   + f_3 * pc_y[k] * snh_282[k];

        t_378[k] = pb_y[k] * smi0_266[k]
                   - f_10 * pc_y[k] * smi1_266[k];

        t_379[k] = f_18 * smh_288[k]
                   + f_3 * pc_x[k] * snh_288[k];

        t_380[k] = f_18 * smh_289[k]
                   + f_3 * pc_x[k] * snh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, smh_290, smh_291, smh_292, smh_293, \
                         snh_290, snh_291, snh_292, snh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_18 * smh_290[k]
                   + f_3 * pc_x[k] * snh_290[k];

        t_382[k] = f_18 * smh_291[k]
                   + f_3 * pc_x[k] * snh_291[k];

        t_383[k] = f_18 * smh_292[k]
                   + f_3 * pc_x[k] * snh_292[k];

        t_384[k] = f_18 * smh_293[k]
                   + f_3 * pc_x[k] * snh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, smh_183, smh_204, smh_206, sng0_205, \
                         sng0_207, sng1_205, sng1_207, snh_288, \
                         snh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * smh_204[k]
                   + f_1 * sng0_205[k]
                   - f_2 * sng1_205[k]
                   + f_3 * pc_y[k] * snh_288[k];

        t_386[k] = f_13 * smh_183[k]
                   + f_3 * pc_z[k] * snh_288[k];

        t_387[k] = f_11 * smh_206[k]
                   + f_4 * sng0_207[k]
                   - f_5 * sng1_207[k]
                   + f_3 * pc_y[k] * snh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, smh_207, smh_208, smh_209, sng0_208, \
                         sng0_209, sng1_208, sng1_209, snh_291, snh_292, \
                         snh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * smh_207[k]
                   + f_6 * sng0_208[k]
                   - f_7 * sng1_208[k]
                   + f_3 * pc_y[k] * snh_291[k];

        t_389[k] = f_11 * smh_208[k]
                   + f_8 * sng0_209[k]
                   - f_9 * sng1_209[k]
                   + f_3 * pc_y[k] * snh_292[k];

        t_390[k] = f_11 * smh_209[k]
                   + f_3 * pc_y[k] * snh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_y, pc_x, pc_y, pc_z, smi0_279, \
                         smh_189, smh_294, smi1_279, sng0_210, sng1_210, \
                         snh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pb_y[k] * smi0_279[k]
                   - f_10 * pc_y[k] * smi1_279[k];

        t_392[k] = f_18 * smh_294[k]
                   + f_1 * sng0_210[k]
                   - f_2 * sng1_210[k]
                   + f_3 * pc_x[k] * snh_294[k];

        t_393[k] = f_3 * pc_y[k] * snh_294[k];

        t_394[k] = f_14 * smh_189[k]
                   + f_3 * pc_z[k] * snh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, smh_297, smh_299, sng0_213, \
                         sng0_215, sng1_213, sng1_215, snh_296, snh_297, \
                         snh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_18 * smh_297[k]
                   + f_4 * sng0_213[k]
                   - f_5 * sng1_213[k]
                   + f_3 * pc_x[k] * snh_297[k];

        t_396[k] = f_3 * pc_y[k] * snh_296[k];

        t_397[k] = f_18 * smh_299[k]
                   + f_4 * sng0_215[k]
                   - f_5 * sng1_215[k]
                   + f_3 * pc_x[k] * snh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_x, pc_y, pc_z, smh_192, smh_300, sng0_216, \
                         sng1_216, snh_297, snh_299, snh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_18 * smh_300[k]
                   + f_6 * sng0_216[k]
                   - f_7 * sng1_216[k]
                   + f_3 * pc_x[k] * snh_300[k];

        t_399[k] = f_14 * smh_192[k]
                   + f_3 * pc_z[k] * snh_297[k];

        t_400[k] = f_3 * pc_y[k] * snh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_z, smh_195, smh_303, smh_304, sng0_219, \
                         sng0_220, sng1_219, sng1_220, snh_300, snh_303, \
                         snh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_18 * smh_303[k]
                   + f_6 * sng0_219[k]
                   - f_7 * sng1_219[k]
                   + f_3 * pc_x[k] * snh_303[k];

        t_402[k] = f_18 * smh_304[k]
                   + f_8 * sng0_220[k]
                   - f_9 * sng1_220[k]
                   + f_3 * pc_x[k] * snh_304[k];

        t_403[k] = f_14 * smh_195[k]
                   + f_3 * pc_z[k] * snh_300[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pc_x, pc_y, smh_306, smh_308, sng0_222, \
                         sng0_224, sng1_222, sng1_224, snh_303, snh_306, \
                         snh_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_18 * smh_306[k]
                   + f_8 * sng0_222[k]
                   - f_9 * sng1_222[k]
                   + f_3 * pc_x[k] * snh_306[k];

        t_405[k] = f_3 * pc_y[k] * snh_303[k];

        t_406[k] = f_18 * smh_308[k]
                   + f_8 * sng0_224[k]
                   - f_9 * sng1_224[k]
                   + f_3 * pc_x[k] * snh_308[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pc_x, smh_309, smh_310, smh_311, \
                         smh_312, smh_313, snh_309, snh_310, snh_311, snh_312, \
                         snh_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_18 * smh_309[k]
                   + f_3 * pc_x[k] * snh_309[k];

        t_408[k] = f_18 * smh_310[k]
                   + f_3 * pc_x[k] * snh_310[k];

        t_409[k] = f_18 * smh_311[k]
                   + f_3 * pc_x[k] * snh_311[k];

        t_410[k] = f_18 * smh_312[k]
                   + f_3 * pc_x[k] * snh_312[k];

        t_411[k] = f_18 * smh_313[k]
                   + f_3 * pc_x[k] * snh_313[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pc_x, pc_y, pc_z, smh_204, smh_314, \
                         sng0_220, sng0_222, sng1_220, sng1_222, snh_309, snh_311, \
                         snh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_18 * smh_314[k]
                   + f_3 * pc_x[k] * snh_314[k];

        t_413[k] = f_1 * sng0_220[k]
                   - f_2 * sng1_220[k]
                   + f_3 * pc_y[k] * snh_309[k];

        t_414[k] = f_14 * smh_204[k]
                   + f_3 * pc_z[k] * snh_309[k];

        t_415[k] = f_4 * sng0_222[k]
                   - f_5 * sng1_222[k]
                   + f_3 * pc_y[k] * snh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, smh_209, sng0_223, sng0_224, \
                         sng1_223, sng1_224, snh_312, snh_313, \
                         snh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * sng0_223[k]
                   - f_7 * sng1_223[k]
                   + f_3 * pc_y[k] * snh_312[k];

        t_417[k] = f_8 * sng0_224[k]
                   - f_9 * sng1_224[k]
                   + f_3 * pc_y[k] * snh_313[k];

        t_418[k] = f_3 * pc_y[k] * snh_314[k];

        t_419[k] = f_14 * smh_209[k]
                   + f_1 * sng0_224[k]
                   - f_2 * sng1_224[k]
                   + f_3 * pc_z[k] * snh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, smh_210, smh_315, \
                         smh_318, sng0_225, sng0_228, sng1_225, sng1_228, snh_315, \
                         snh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_19 * smh_315[k]
                   + f_1 * sng0_225[k]
                   - f_2 * sng1_225[k]
                   + f_3 * pc_x[k] * snh_315[k];

        t_421[k] = f_19 * smh_210[k]
                   + f_3 * pc_y[k] * snh_315[k];

        t_422[k] = f_3 * pc_z[k] * snh_315[k];

        t_423[k] = f_19 * smh_318[k]
                   + f_4 * sng0_228[k]
                   - f_5 * sng1_228[k]
                   + f_3 * pc_x[k] * snh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, smh_212, smh_320, smh_321, sng0_230, \
                         sng0_231, sng1_230, sng1_231, snh_317, snh_320, \
                         snh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_19 * smh_212[k]
                   + f_3 * pc_y[k] * snh_317[k];

        t_425[k] = f_19 * smh_320[k]
                   + f_4 * sng0_230[k]
                   - f_5 * sng1_230[k]
                   + f_3 * pc_x[k] * snh_320[k];

        t_426[k] = f_19 * smh_321[k]
                   + f_6 * sng0_231[k]
                   - f_7 * sng1_231[k]
                   + f_3 * pc_x[k] * snh_321[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_x, pc_y, pc_z, smh_215, smh_324, sng0_234, \
                         sng1_234, snh_318, snh_320, snh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_3 * pc_z[k] * snh_318[k];

        t_428[k] = f_19 * smh_215[k]
                   + f_3 * pc_y[k] * snh_320[k];

        t_429[k] = f_19 * smh_324[k]
                   + f_6 * sng0_234[k]
                   - f_7 * sng1_234[k]
                   + f_3 * pc_x[k] * snh_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_z, smh_325, smh_327, sng0_235, \
                         sng0_237, sng1_235, sng1_237, snh_321, snh_325, \
                         snh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_19 * smh_325[k]
                   + f_8 * sng0_235[k]
                   - f_9 * sng1_235[k]
                   + f_3 * pc_x[k] * snh_325[k];

        t_431[k] = f_3 * pc_z[k] * snh_321[k];

        t_432[k] = f_19 * smh_327[k]
                   + f_8 * sng0_237[k]
                   - f_9 * sng1_237[k]
                   + f_3 * pc_x[k] * snh_327[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, smh_219, smh_329, smh_330, \
                         smh_331, sng0_239, sng1_239, snh_324, snh_329, snh_330, \
                         snh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_19 * smh_219[k]
                   + f_3 * pc_y[k] * snh_324[k];

        t_434[k] = f_19 * smh_329[k]
                   + f_8 * sng0_239[k]
                   - f_9 * sng1_239[k]
                   + f_3 * pc_x[k] * snh_329[k];

        t_435[k] = f_19 * smh_330[k]
                   + f_3 * pc_x[k] * snh_330[k];

        t_436[k] = f_19 * smh_331[k]
                   + f_3 * pc_x[k] * snh_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, smh_332, smh_333, smh_334, smh_335, \
                         snh_332, snh_333, snh_334, snh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_19 * smh_332[k]
                   + f_3 * pc_x[k] * snh_332[k];

        t_438[k] = f_19 * smh_333[k]
                   + f_3 * pc_x[k] * snh_333[k];

        t_439[k] = f_19 * smh_334[k]
                   + f_3 * pc_x[k] * snh_334[k];

        t_440[k] = f_19 * smh_335[k]
                   + f_3 * pc_x[k] * snh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, pc_z, smh_225, smh_227, sng0_235, \
                         sng0_237, sng1_235, sng1_237, snh_330, \
                         snh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_19 * smh_225[k]
                   + f_1 * sng0_235[k]
                   - f_2 * sng1_235[k]
                   + f_3 * pc_y[k] * snh_330[k];

        t_442[k] = f_3 * pc_z[k] * snh_330[k];

        t_443[k] = f_19 * smh_227[k]
                   + f_4 * sng0_237[k]
                   - f_5 * sng1_237[k]
                   + f_3 * pc_y[k] * snh_332[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, smh_228, smh_229, smh_230, \
                         sng0_238, sng0_239, sng1_238, sng1_239, snh_333, snh_334, \
                         snh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_19 * smh_228[k]
                   + f_6 * sng0_238[k]
                   - f_7 * sng1_238[k]
                   + f_3 * pc_y[k] * snh_333[k];

        t_445[k] = f_19 * smh_229[k]
                   + f_8 * sng0_239[k]
                   - f_9 * sng1_239[k]
                   + f_3 * pc_y[k] * snh_334[k];

        t_446[k] = f_19 * smh_230[k]
                   + f_3 * pc_y[k] * snh_335[k];

        t_447[k] = f_1 * sng0_239[k]
                   - f_2 * sng1_239[k]
                   + f_3 * pc_z[k] * snh_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_z, pc_y, pc_z, smi0_280, smi0_283, \
                         smh_210, smh_231, smi1_280, smi1_283, \
                         snh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_z[k] * smi0_280[k]
                   - f_10 * pc_z[k] * smi1_280[k];

        t_449[k] = f_14 * smh_231[k]
                   + f_3 * pc_y[k] * snh_336[k];

        t_450[k] = f_11 * smh_210[k]
                   + f_3 * pc_z[k] * snh_336[k];

        t_451[k] = pb_z[k] * smi0_283[k]
                   - f_10 * pc_z[k] * smi1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_z, pc_x, pc_y, pc_z, smi0_286, smh_233, \
                         smh_341, smi1_286, sng0_245, sng1_245, snh_338, \
                         snh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * smh_233[k]
                   + f_3 * pc_y[k] * snh_338[k];

        t_453[k] = f_19 * smh_341[k]
                   + f_4 * sng0_245[k]
                   - f_5 * sng1_245[k]
                   + f_3 * pc_x[k] * snh_341[k];

        t_454[k] = pb_z[k] * smi0_286[k]
                   - f_10 * pc_z[k] * smi1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_y, pc_z, smh_213, smh_236, smh_345, \
                         sng0_249, sng1_249, snh_339, snh_341, \
                         snh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * smh_213[k]
                   + f_3 * pc_z[k] * snh_339[k];

        t_456[k] = f_14 * smh_236[k]
                   + f_3 * pc_y[k] * snh_341[k];

        t_457[k] = f_19 * smh_345[k]
                   + f_6 * sng0_249[k]
                   - f_7 * sng1_249[k]
                   + f_3 * pc_x[k] * snh_345[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pb_z, pc_y, pc_z, smi0_290, smi0_292, \
                         smh_216, smh_217, smh_240, smi1_290, smi1_292, snh_342, \
                         snh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pb_z[k] * smi0_290[k]
                   - f_10 * pc_z[k] * smi1_290[k];

        t_459[k] = f_11 * smh_216[k]
                   + f_3 * pc_z[k] * snh_342[k];

        t_460[k] = pb_z[k] * smi0_292[k]
                   + f_12 * smh_217[k]
                   - f_10 * pc_z[k] * smi1_292[k];

        t_461[k] = f_14 * smh_240[k]
                   + f_3 * pc_y[k] * snh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, smh_350, smh_351, smh_352, smh_353, \
                         sng0_254, sng1_254, snh_350, snh_351, snh_352, \
                         snh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_19 * smh_350[k]
                   + f_8 * sng0_254[k]
                   - f_9 * sng1_254[k]
                   + f_3 * pc_x[k] * snh_350[k];

        t_463[k] = f_19 * smh_351[k]
                   + f_3 * pc_x[k] * snh_351[k];

        t_464[k] = f_19 * smh_352[k]
                   + f_3 * pc_x[k] * snh_352[k];

        t_465[k] = f_19 * smh_353[k]
                   + f_3 * pc_x[k] * snh_353[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_19 = 2.5 / q;

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

    const auto *smi0_301 = buffer.data(smi0 + 301);
    const auto *smi0_392 = buffer.data(smi0 + 392);
    const auto *smi0_395 = buffer.data(smi0 + 395);
    const auto *smi0_397 = buffer.data(smi0 + 397);
    const auto *smi0_398 = buffer.data(smi0 + 398);
    const auto *smi0_401 = buffer.data(smi0 + 401);
    const auto *smi0_402 = buffer.data(smi0 + 402);
    const auto *smi0_404 = buffer.data(smi0 + 404);
    const auto *smi0_406 = buffer.data(smi0 + 406);
    const auto *smi0_419 = buffer.data(smi0 + 419);

    const auto *smh_225 = buffer.data(smh + 225);
    const auto *smh_230 = buffer.data(smh + 230);
    const auto *smh_231 = buffer.data(smh + 231);
    const auto *smh_234 = buffer.data(smh + 234);
    const auto *smh_237 = buffer.data(smh + 237);
    const auto *smh_246 = buffer.data(smh + 246);
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
    const auto *smh_267 = buffer.data(smh + 267);
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
    const auto *smh_290 = buffer.data(smh + 290);
    const auto *smh_291 = buffer.data(smh + 291);
    const auto *smh_292 = buffer.data(smh + 292);
    const auto *smh_293 = buffer.data(smh + 293);
    const auto *smh_294 = buffer.data(smh + 294);
    const auto *smh_295 = buffer.data(smh + 295);
    const auto *smh_296 = buffer.data(smh + 296);
    const auto *smh_297 = buffer.data(smh + 297);
    const auto *smh_299 = buffer.data(smh + 299);
    const auto *smh_300 = buffer.data(smh + 300);
    const auto *smh_302 = buffer.data(smh + 302);
    const auto *smh_303 = buffer.data(smh + 303);
    const auto *smh_309 = buffer.data(smh + 309);
    const auto *smh_311 = buffer.data(smh + 311);
    const auto *smh_312 = buffer.data(smh + 312);
    const auto *smh_313 = buffer.data(smh + 313);
    const auto *smh_314 = buffer.data(smh + 314);
    const auto *smh_354 = buffer.data(smh + 354);
    const auto *smh_355 = buffer.data(smh + 355);
    const auto *smh_356 = buffer.data(smh + 356);
    const auto *smh_357 = buffer.data(smh + 357);
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
    const auto *smh_414 = buffer.data(smh + 414);
    const auto *smh_415 = buffer.data(smh + 415);
    const auto *smh_416 = buffer.data(smh + 416);
    const auto *smh_417 = buffer.data(smh + 417);
    const auto *smh_418 = buffer.data(smh + 418);
    const auto *smh_419 = buffer.data(smh + 419);
    const auto *smh_420 = buffer.data(smh + 420);
    const auto *smh_423 = buffer.data(smh + 423);
    const auto *smh_425 = buffer.data(smh + 425);
    const auto *smh_426 = buffer.data(smh + 426);
    const auto *smh_429 = buffer.data(smh + 429);
    const auto *smh_430 = buffer.data(smh + 430);
    const auto *smh_432 = buffer.data(smh + 432);
    const auto *smh_434 = buffer.data(smh + 434);

    const auto *smi1_301 = buffer.data(smi1 + 301);
    const auto *smi1_392 = buffer.data(smi1 + 392);
    const auto *smi1_395 = buffer.data(smi1 + 395);
    const auto *smi1_397 = buffer.data(smi1 + 397);
    const auto *smi1_398 = buffer.data(smi1 + 398);
    const auto *smi1_401 = buffer.data(smi1 + 401);
    const auto *smi1_402 = buffer.data(smi1 + 402);
    const auto *smi1_404 = buffer.data(smi1 + 404);
    const auto *smi1_406 = buffer.data(smi1 + 406);
    const auto *smi1_419 = buffer.data(smi1 + 419);

    const auto *sng0_252 = buffer.data(sng0 + 252);
    const auto *sng0_253 = buffer.data(sng0 + 253);
    const auto *sng0_254 = buffer.data(sng0 + 254);
    const auto *sng0_255 = buffer.data(sng0 + 255);
    const auto *sng0_258 = buffer.data(sng0 + 258);
    const auto *sng0_260 = buffer.data(sng0 + 260);
    const auto *sng0_261 = buffer.data(sng0 + 261);
    const auto *sng0_264 = buffer.data(sng0 + 264);
    const auto *sng0_265 = buffer.data(sng0 + 265);
    const auto *sng0_267 = buffer.data(sng0 + 267);
    const auto *sng0_268 = buffer.data(sng0 + 268);
    const auto *sng0_269 = buffer.data(sng0 + 269);
    const auto *sng0_270 = buffer.data(sng0 + 270);
    const auto *sng0_273 = buffer.data(sng0 + 273);
    const auto *sng0_275 = buffer.data(sng0 + 275);
    const auto *sng0_276 = buffer.data(sng0 + 276);
    const auto *sng0_279 = buffer.data(sng0 + 279);
    const auto *sng0_280 = buffer.data(sng0 + 280);
    const auto *sng0_282 = buffer.data(sng0 + 282);
    const auto *sng0_283 = buffer.data(sng0 + 283);
    const auto *sng0_284 = buffer.data(sng0 + 284);
    const auto *sng0_295 = buffer.data(sng0 + 295);
    const auto *sng0_297 = buffer.data(sng0 + 297);
    const auto *sng0_298 = buffer.data(sng0 + 298);
    const auto *sng0_299 = buffer.data(sng0 + 299);
    const auto *sng0_300 = buffer.data(sng0 + 300);
    const auto *sng0_303 = buffer.data(sng0 + 303);
    const auto *sng0_305 = buffer.data(sng0 + 305);
    const auto *sng0_306 = buffer.data(sng0 + 306);
    const auto *sng0_309 = buffer.data(sng0 + 309);
    const auto *sng0_310 = buffer.data(sng0 + 310);
    const auto *sng0_312 = buffer.data(sng0 + 312);
    const auto *sng0_314 = buffer.data(sng0 + 314);

    const auto *sng1_252 = buffer.data(sng1 + 252);
    const auto *sng1_253 = buffer.data(sng1 + 253);
    const auto *sng1_254 = buffer.data(sng1 + 254);
    const auto *sng1_255 = buffer.data(sng1 + 255);
    const auto *sng1_258 = buffer.data(sng1 + 258);
    const auto *sng1_260 = buffer.data(sng1 + 260);
    const auto *sng1_261 = buffer.data(sng1 + 261);
    const auto *sng1_264 = buffer.data(sng1 + 264);
    const auto *sng1_265 = buffer.data(sng1 + 265);
    const auto *sng1_267 = buffer.data(sng1 + 267);
    const auto *sng1_268 = buffer.data(sng1 + 268);
    const auto *sng1_269 = buffer.data(sng1 + 269);
    const auto *sng1_270 = buffer.data(sng1 + 270);
    const auto *sng1_273 = buffer.data(sng1 + 273);
    const auto *sng1_275 = buffer.data(sng1 + 275);
    const auto *sng1_276 = buffer.data(sng1 + 276);
    const auto *sng1_279 = buffer.data(sng1 + 279);
    const auto *sng1_280 = buffer.data(sng1 + 280);
    const auto *sng1_282 = buffer.data(sng1 + 282);
    const auto *sng1_283 = buffer.data(sng1 + 283);
    const auto *sng1_284 = buffer.data(sng1 + 284);
    const auto *sng1_295 = buffer.data(sng1 + 295);
    const auto *sng1_297 = buffer.data(sng1 + 297);
    const auto *sng1_298 = buffer.data(sng1 + 298);
    const auto *sng1_299 = buffer.data(sng1 + 299);
    const auto *sng1_300 = buffer.data(sng1 + 300);
    const auto *sng1_303 = buffer.data(sng1 + 303);
    const auto *sng1_305 = buffer.data(sng1 + 305);
    const auto *sng1_306 = buffer.data(sng1 + 306);
    const auto *sng1_309 = buffer.data(sng1 + 309);
    const auto *sng1_310 = buffer.data(sng1 + 310);
    const auto *sng1_312 = buffer.data(sng1 + 312);
    const auto *sng1_314 = buffer.data(sng1 + 314);

    const auto *snh_351 = buffer.data(snh + 351);
    const auto *snh_353 = buffer.data(snh + 353);
    const auto *snh_354 = buffer.data(snh + 354);
    const auto *snh_355 = buffer.data(snh + 355);
    const auto *snh_356 = buffer.data(snh + 356);
    const auto *snh_357 = buffer.data(snh + 357);
    const auto *snh_359 = buffer.data(snh + 359);
    const auto *snh_360 = buffer.data(snh + 360);
    const auto *snh_362 = buffer.data(snh + 362);
    const auto *snh_363 = buffer.data(snh + 363);
    const auto *snh_366 = buffer.data(snh + 366);
    const auto *snh_367 = buffer.data(snh + 367);
    const auto *snh_369 = buffer.data(snh + 369);
    const auto *snh_371 = buffer.data(snh + 371);
    const auto *snh_372 = buffer.data(snh + 372);
    const auto *snh_373 = buffer.data(snh + 373);
    const auto *snh_374 = buffer.data(snh + 374);
    const auto *snh_375 = buffer.data(snh + 375);
    const auto *snh_376 = buffer.data(snh + 376);
    const auto *snh_377 = buffer.data(snh + 377);
    const auto *snh_378 = buffer.data(snh + 378);
    const auto *snh_380 = buffer.data(snh + 380);
    const auto *snh_381 = buffer.data(snh + 381);
    const auto *snh_383 = buffer.data(snh + 383);
    const auto *snh_384 = buffer.data(snh + 384);
    const auto *snh_387 = buffer.data(snh + 387);
    const auto *snh_388 = buffer.data(snh + 388);
    const auto *snh_390 = buffer.data(snh + 390);
    const auto *snh_392 = buffer.data(snh + 392);
    const auto *snh_393 = buffer.data(snh + 393);
    const auto *snh_394 = buffer.data(snh + 394);
    const auto *snh_395 = buffer.data(snh + 395);
    const auto *snh_396 = buffer.data(snh + 396);
    const auto *snh_397 = buffer.data(snh + 397);
    const auto *snh_398 = buffer.data(snh + 398);
    const auto *snh_399 = buffer.data(snh + 399);
    const auto *snh_401 = buffer.data(snh + 401);
    const auto *snh_402 = buffer.data(snh + 402);
    const auto *snh_404 = buffer.data(snh + 404);
    const auto *snh_405 = buffer.data(snh + 405);
    const auto *snh_408 = buffer.data(snh + 408);
    const auto *snh_414 = buffer.data(snh + 414);
    const auto *snh_415 = buffer.data(snh + 415);
    const auto *snh_416 = buffer.data(snh + 416);
    const auto *snh_417 = buffer.data(snh + 417);
    const auto *snh_418 = buffer.data(snh + 418);
    const auto *snh_419 = buffer.data(snh + 419);
    const auto *snh_420 = buffer.data(snh + 420);
    const auto *snh_422 = buffer.data(snh + 422);
    const auto *snh_423 = buffer.data(snh + 423);
    const auto *snh_425 = buffer.data(snh + 425);
    const auto *snh_426 = buffer.data(snh + 426);
    const auto *snh_429 = buffer.data(snh + 429);
    const auto *snh_430 = buffer.data(snh + 430);
    const auto *snh_432 = buffer.data(snh + 432);
    const auto *snh_434 = buffer.data(snh + 434);

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_z, pc_x, pc_z, smi0_301, smh_354, \
                         smh_355, smh_356, smi1_301, snh_354, snh_355, \
                         snh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_19 * smh_354[k]
                   + f_3 * pc_x[k] * snh_354[k];

        t_467[k] = f_19 * smh_355[k]
                   + f_3 * pc_x[k] * snh_355[k];

        t_468[k] = f_19 * smh_356[k]
                   + f_3 * pc_x[k] * snh_356[k];

        t_469[k] = pb_z[k] * smi0_301[k]
                   - f_10 * pc_z[k] * smi1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, smh_225, smh_248, smh_249, sng0_252, \
                         sng0_253, sng1_252, sng1_253, snh_351, snh_353, \
                         snh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * smh_225[k]
                   + f_3 * pc_z[k] * snh_351[k];

        t_471[k] = f_14 * smh_248[k]
                   + f_4 * sng0_252[k]
                   - f_5 * sng1_252[k]
                   + f_3 * pc_y[k] * snh_353[k];

        t_472[k] = f_14 * smh_249[k]
                   + f_6 * sng0_253[k]
                   - f_7 * sng1_253[k]
                   + f_3 * pc_y[k] * snh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, smh_230, smh_250, smh_251, sng0_254, \
                         sng1_254, snh_355, snh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * smh_250[k]
                   + f_8 * sng0_254[k]
                   - f_9 * sng1_254[k]
                   + f_3 * pc_y[k] * snh_355[k];

        t_474[k] = f_14 * smh_251[k]
                   + f_3 * pc_y[k] * snh_356[k];

        t_475[k] = f_11 * smh_230[k]
                   + f_1 * sng0_254[k]
                   - f_2 * sng1_254[k]
                   + f_3 * pc_z[k] * snh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, smh_231, smh_252, smh_357, \
                         sng0_255, sng1_255, snh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_19 * smh_357[k]
                   + f_1 * sng0_255[k]
                   - f_2 * sng1_255[k]
                   + f_3 * pc_x[k] * snh_357[k];

        t_477[k] = f_13 * smh_252[k]
                   + f_3 * pc_y[k] * snh_357[k];

        t_478[k] = f_12 * smh_231[k]
                   + f_3 * pc_z[k] * snh_357[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, smh_254, smh_360, smh_362, sng0_258, \
                         sng0_260, sng1_258, sng1_260, snh_359, snh_360, \
                         snh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_19 * smh_360[k]
                   + f_4 * sng0_258[k]
                   - f_5 * sng1_258[k]
                   + f_3 * pc_x[k] * snh_360[k];

        t_480[k] = f_13 * smh_254[k]
                   + f_3 * pc_y[k] * snh_359[k];

        t_481[k] = f_19 * smh_362[k]
                   + f_4 * sng0_260[k]
                   - f_5 * sng1_260[k]
                   + f_3 * pc_x[k] * snh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, smh_234, smh_257, smh_363, \
                         sng0_261, sng1_261, snh_360, snh_362, \
                         snh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_19 * smh_363[k]
                   + f_6 * sng0_261[k]
                   - f_7 * sng1_261[k]
                   + f_3 * pc_x[k] * snh_363[k];

        t_483[k] = f_12 * smh_234[k]
                   + f_3 * pc_z[k] * snh_360[k];

        t_484[k] = f_13 * smh_257[k]
                   + f_3 * pc_y[k] * snh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, smh_237, smh_366, smh_367, sng0_264, \
                         sng0_265, sng1_264, sng1_265, snh_363, snh_366, \
                         snh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_19 * smh_366[k]
                   + f_6 * sng0_264[k]
                   - f_7 * sng1_264[k]
                   + f_3 * pc_x[k] * snh_366[k];

        t_486[k] = f_19 * smh_367[k]
                   + f_8 * sng0_265[k]
                   - f_9 * sng1_265[k]
                   + f_3 * pc_x[k] * snh_367[k];

        t_487[k] = f_12 * smh_237[k]
                   + f_3 * pc_z[k] * snh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, smh_261, smh_369, smh_371, sng0_267, \
                         sng0_269, sng1_267, sng1_269, snh_366, snh_369, \
                         snh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_19 * smh_369[k]
                   + f_8 * sng0_267[k]
                   - f_9 * sng1_267[k]
                   + f_3 * pc_x[k] * snh_369[k];

        t_489[k] = f_13 * smh_261[k]
                   + f_3 * pc_y[k] * snh_366[k];

        t_490[k] = f_19 * smh_371[k]
                   + f_8 * sng0_269[k]
                   - f_9 * sng1_269[k]
                   + f_3 * pc_x[k] * snh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, smh_372, smh_373, smh_374, \
                         smh_375, smh_376, snh_372, snh_373, snh_374, snh_375, \
                         snh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_19 * smh_372[k]
                   + f_3 * pc_x[k] * snh_372[k];

        t_492[k] = f_19 * smh_373[k]
                   + f_3 * pc_x[k] * snh_373[k];

        t_493[k] = f_19 * smh_374[k]
                   + f_3 * pc_x[k] * snh_374[k];

        t_494[k] = f_19 * smh_375[k]
                   + f_3 * pc_x[k] * snh_375[k];

        t_495[k] = f_19 * smh_376[k]
                   + f_3 * pc_x[k] * snh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, smh_246, smh_267, smh_377, \
                         sng0_265, sng1_265, snh_372, snh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_19 * smh_377[k]
                   + f_3 * pc_x[k] * snh_377[k];

        t_497[k] = f_13 * smh_267[k]
                   + f_1 * sng0_265[k]
                   - f_2 * sng1_265[k]
                   + f_3 * pc_y[k] * snh_372[k];

        t_498[k] = f_12 * smh_246[k]
                   + f_3 * pc_z[k] * snh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, smh_269, smh_270, smh_271, sng0_267, \
                         sng0_268, sng0_269, sng1_267, sng1_268, sng1_269, snh_374, snh_375, \
                         snh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * smh_269[k]
                   + f_4 * sng0_267[k]
                   - f_5 * sng1_267[k]
                   + f_3 * pc_y[k] * snh_374[k];

        t_500[k] = f_13 * smh_270[k]
                   + f_6 * sng0_268[k]
                   - f_7 * sng1_268[k]
                   + f_3 * pc_y[k] * snh_375[k];

        t_501[k] = f_13 * smh_271[k]
                   + f_8 * sng0_269[k]
                   - f_9 * sng1_269[k]
                   + f_3 * pc_y[k] * snh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, smh_251, smh_272, smh_378, \
                         sng0_269, sng0_270, sng1_269, sng1_270, snh_377, \
                         snh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * smh_272[k]
                   + f_3 * pc_y[k] * snh_377[k];

        t_503[k] = f_12 * smh_251[k]
                   + f_1 * sng0_269[k]
                   - f_2 * sng1_269[k]
                   + f_3 * pc_z[k] * snh_377[k];

        t_504[k] = f_19 * smh_378[k]
                   + f_1 * sng0_270[k]
                   - f_2 * sng1_270[k]
                   + f_3 * pc_x[k] * snh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, smh_252, smh_273, \
                         smh_275, smh_381, sng0_273, sng1_273, snh_378, snh_380, \
                         snh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * smh_273[k]
                   + f_3 * pc_y[k] * snh_378[k];

        t_506[k] = f_13 * smh_252[k]
                   + f_3 * pc_z[k] * snh_378[k];

        t_507[k] = f_19 * smh_381[k]
                   + f_4 * sng0_273[k]
                   - f_5 * sng1_273[k]
                   + f_3 * pc_x[k] * snh_381[k];

        t_508[k] = f_12 * smh_275[k]
                   + f_3 * pc_y[k] * snh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, smh_255, smh_383, smh_384, sng0_275, \
                         sng0_276, sng1_275, sng1_276, snh_381, snh_383, \
                         snh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_19 * smh_383[k]
                   + f_4 * sng0_275[k]
                   - f_5 * sng1_275[k]
                   + f_3 * pc_x[k] * snh_383[k];

        t_510[k] = f_19 * smh_384[k]
                   + f_6 * sng0_276[k]
                   - f_7 * sng1_276[k]
                   + f_3 * pc_x[k] * snh_384[k];

        t_511[k] = f_13 * smh_255[k]
                   + f_3 * pc_z[k] * snh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, smh_278, smh_387, smh_388, sng0_279, \
                         sng0_280, sng1_279, sng1_280, snh_383, snh_387, \
                         snh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * smh_278[k]
                   + f_3 * pc_y[k] * snh_383[k];

        t_513[k] = f_19 * smh_387[k]
                   + f_6 * sng0_279[k]
                   - f_7 * sng1_279[k]
                   + f_3 * pc_x[k] * snh_387[k];

        t_514[k] = f_19 * smh_388[k]
                   + f_8 * sng0_280[k]
                   - f_9 * sng1_280[k]
                   + f_3 * pc_x[k] * snh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, smh_258, smh_282, smh_390, \
                         sng0_282, sng1_282, snh_384, snh_387, \
                         snh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * smh_258[k]
                   + f_3 * pc_z[k] * snh_384[k];

        t_516[k] = f_19 * smh_390[k]
                   + f_8 * sng0_282[k]
                   - f_9 * sng1_282[k]
                   + f_3 * pc_x[k] * snh_390[k];

        t_517[k] = f_12 * smh_282[k]
                   + f_3 * pc_y[k] * snh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, smh_392, smh_393, smh_394, smh_395, \
                         sng0_284, sng1_284, snh_392, snh_393, snh_394, \
                         snh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_19 * smh_392[k]
                   + f_8 * sng0_284[k]
                   - f_9 * sng1_284[k]
                   + f_3 * pc_x[k] * snh_392[k];

        t_519[k] = f_19 * smh_393[k]
                   + f_3 * pc_x[k] * snh_393[k];

        t_520[k] = f_19 * smh_394[k]
                   + f_3 * pc_x[k] * snh_394[k];

        t_521[k] = f_19 * smh_395[k]
                   + f_3 * pc_x[k] * snh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, smh_288, smh_396, smh_397, \
                         smh_398, sng0_280, sng1_280, snh_393, snh_396, snh_397, \
                         snh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_19 * smh_396[k]
                   + f_3 * pc_x[k] * snh_396[k];

        t_523[k] = f_19 * smh_397[k]
                   + f_3 * pc_x[k] * snh_397[k];

        t_524[k] = f_19 * smh_398[k]
                   + f_3 * pc_x[k] * snh_398[k];

        t_525[k] = f_12 * smh_288[k]
                   + f_1 * sng0_280[k]
                   - f_2 * sng1_280[k]
                   + f_3 * pc_y[k] * snh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, smh_267, smh_290, smh_291, sng0_282, \
                         sng0_283, sng1_282, sng1_283, snh_393, snh_395, \
                         snh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * smh_267[k]
                   + f_3 * pc_z[k] * snh_393[k];

        t_527[k] = f_12 * smh_290[k]
                   + f_4 * sng0_282[k]
                   - f_5 * sng1_282[k]
                   + f_3 * pc_y[k] * snh_395[k];

        t_528[k] = f_12 * smh_291[k]
                   + f_6 * sng0_283[k]
                   - f_7 * sng1_283[k]
                   + f_3 * pc_y[k] * snh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pb_y, pc_y, pc_z, smi0_392, smh_272, \
                         smh_292, smh_293, smi1_392, sng0_284, sng1_284, snh_397, \
                         snh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * smh_292[k]
                   + f_8 * sng0_284[k]
                   - f_9 * sng1_284[k]
                   + f_3 * pc_y[k] * snh_397[k];

        t_530[k] = f_12 * smh_293[k]
                   + f_3 * pc_y[k] * snh_398[k];

        t_531[k] = f_13 * smh_272[k]
                   + f_1 * sng0_284[k]
                   - f_2 * sng1_284[k]
                   + f_3 * pc_z[k] * snh_398[k];

        t_532[k] = pb_y[k] * smi0_392[k]
                   - f_10 * pc_y[k] * smi1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pb_y, pc_y, pc_z, smi0_395, smh_273, \
                         smh_294, smh_295, smh_296, smi1_395, snh_399, \
                         snh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * smh_294[k]
                   + f_3 * pc_y[k] * snh_399[k];

        t_534[k] = f_14 * smh_273[k]
                   + f_3 * pc_z[k] * snh_399[k];

        t_535[k] = pb_y[k] * smi0_395[k]
                   + f_12 * smh_295[k]
                   - f_10 * pc_y[k] * smi1_395[k];

        t_536[k] = f_11 * smh_296[k]
                   + f_3 * pc_y[k] * snh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pb_y, pc_y, pc_z, smi0_397, smi0_398, \
                         smh_276, smh_297, smh_299, smi1_397, smi1_398, snh_402, \
                         snh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pb_y[k] * smi0_397[k]
                   - f_10 * pc_y[k] * smi1_397[k];

        t_538[k] = pb_y[k] * smi0_398[k]
                   + f_13 * smh_297[k]
                   - f_10 * pc_y[k] * smi1_398[k];

        t_539[k] = f_14 * smh_276[k]
                   + f_3 * pc_z[k] * snh_402[k];

        t_540[k] = f_11 * smh_299[k]
                   + f_3 * pc_y[k] * snh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pb_y, pc_y, pc_z, smi0_401, smi0_402, smh_279, \
                         smh_300, smi1_401, smi1_402, snh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pb_y[k] * smi0_401[k]
                   - f_10 * pc_y[k] * smi1_401[k];

        t_542[k] = pb_y[k] * smi0_402[k]
                   + f_14 * smh_300[k]
                   - f_10 * pc_y[k] * smi1_402[k];

        t_543[k] = f_14 * smh_279[k]
                   + f_3 * pc_z[k] * snh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_x, pc_y, smi0_404, smi0_406, \
                         smh_302, smh_303, smh_414, smi1_404, smi1_406, snh_408, \
                         snh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pb_y[k] * smi0_404[k]
                   + f_12 * smh_302[k]
                   - f_10 * pc_y[k] * smi1_404[k];

        t_545[k] = f_11 * smh_303[k]
                   + f_3 * pc_y[k] * snh_408[k];

        t_546[k] = pb_y[k] * smi0_406[k]
                   - f_10 * pc_y[k] * smi1_406[k];

        t_547[k] = f_19 * smh_414[k]
                   + f_3 * pc_x[k] * snh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, smh_415, smh_416, smh_417, \
                         smh_418, smh_419, snh_415, snh_416, snh_417, snh_418, \
                         snh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_19 * smh_415[k]
                   + f_3 * pc_x[k] * snh_415[k];

        t_549[k] = f_19 * smh_416[k]
                   + f_3 * pc_x[k] * snh_416[k];

        t_550[k] = f_19 * smh_417[k]
                   + f_3 * pc_x[k] * snh_417[k];

        t_551[k] = f_19 * smh_418[k]
                   + f_3 * pc_x[k] * snh_418[k];

        t_552[k] = f_19 * smh_419[k]
                   + f_3 * pc_x[k] * snh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, smh_288, smh_309, smh_311, sng0_295, \
                         sng0_297, sng1_295, sng1_297, snh_414, \
                         snh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * smh_309[k]
                   + f_1 * sng0_295[k]
                   - f_2 * sng1_295[k]
                   + f_3 * pc_y[k] * snh_414[k];

        t_554[k] = f_14 * smh_288[k]
                   + f_3 * pc_z[k] * snh_414[k];

        t_555[k] = f_11 * smh_311[k]
                   + f_4 * sng0_297[k]
                   - f_5 * sng1_297[k]
                   + f_3 * pc_y[k] * snh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, smh_312, smh_313, smh_314, sng0_298, \
                         sng0_299, sng1_298, sng1_299, snh_417, snh_418, \
                         snh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * smh_312[k]
                   + f_6 * sng0_298[k]
                   - f_7 * sng1_298[k]
                   + f_3 * pc_y[k] * snh_417[k];

        t_557[k] = f_11 * smh_313[k]
                   + f_8 * sng0_299[k]
                   - f_9 * sng1_299[k]
                   + f_3 * pc_y[k] * snh_418[k];

        t_558[k] = f_11 * smh_314[k]
                   + f_3 * pc_y[k] * snh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_y, pc_x, pc_y, pc_z, smi0_419, \
                         smh_294, smh_420, smi1_419, sng0_300, sng1_300, \
                         snh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pb_y[k] * smi0_419[k]
                   - f_10 * pc_y[k] * smi1_419[k];

        t_560[k] = f_19 * smh_420[k]
                   + f_1 * sng0_300[k]
                   - f_2 * sng1_300[k]
                   + f_3 * pc_x[k] * snh_420[k];

        t_561[k] = f_3 * pc_y[k] * snh_420[k];

        t_562[k] = f_19 * smh_294[k]
                   + f_3 * pc_z[k] * snh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, smh_423, smh_425, sng0_303, \
                         sng0_305, sng1_303, sng1_305, snh_422, snh_423, \
                         snh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_19 * smh_423[k]
                   + f_4 * sng0_303[k]
                   - f_5 * sng1_303[k]
                   + f_3 * pc_x[k] * snh_423[k];

        t_564[k] = f_3 * pc_y[k] * snh_422[k];

        t_565[k] = f_19 * smh_425[k]
                   + f_4 * sng0_305[k]
                   - f_5 * sng1_305[k]
                   + f_3 * pc_x[k] * snh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_x, pc_y, pc_z, smh_297, smh_426, sng0_306, \
                         sng1_306, snh_423, snh_425, snh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_19 * smh_426[k]
                   + f_6 * sng0_306[k]
                   - f_7 * sng1_306[k]
                   + f_3 * pc_x[k] * snh_426[k];

        t_567[k] = f_19 * smh_297[k]
                   + f_3 * pc_z[k] * snh_423[k];

        t_568[k] = f_3 * pc_y[k] * snh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_z, smh_300, smh_429, smh_430, sng0_309, \
                         sng0_310, sng1_309, sng1_310, snh_426, snh_429, \
                         snh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_19 * smh_429[k]
                   + f_6 * sng0_309[k]
                   - f_7 * sng1_309[k]
                   + f_3 * pc_x[k] * snh_429[k];

        t_570[k] = f_19 * smh_430[k]
                   + f_8 * sng0_310[k]
                   - f_9 * sng1_310[k]
                   + f_3 * pc_x[k] * snh_430[k];

        t_571[k] = f_19 * smh_300[k]
                   + f_3 * pc_z[k] * snh_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pc_x, pc_y, smh_432, smh_434, sng0_312, \
                         sng0_314, sng1_312, sng1_314, snh_429, snh_432, \
                         snh_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_19 * smh_432[k]
                   + f_8 * sng0_312[k]
                   - f_9 * sng1_312[k]
                   + f_3 * pc_x[k] * snh_432[k];

        t_573[k] = f_3 * pc_y[k] * snh_429[k];

        t_574[k] = f_19 * smh_434[k]
                   + f_8 * sng0_314[k]
                   - f_9 * sng1_314[k]
                   + f_3 * pc_x[k] * snh_434[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *smi0_420 = buffer.data(smi0 + 420);
    const auto *smi0_423 = buffer.data(smi0 + 423);
    const auto *smi0_426 = buffer.data(smi0 + 426);
    const auto *smi0_430 = buffer.data(smi0 + 430);
    const auto *smi0_432 = buffer.data(smi0 + 432);
    const auto *smi0_441 = buffer.data(smi0 + 441);

    const auto *smh_309 = buffer.data(smh + 309);
    const auto *smh_314 = buffer.data(smh + 314);
    const auto *smh_315 = buffer.data(smh + 315);
    const auto *smh_317 = buffer.data(smh + 317);
    const auto *smh_318 = buffer.data(smh + 318);
    const auto *smh_320 = buffer.data(smh + 320);
    const auto *smh_321 = buffer.data(smh + 321);
    const auto *smh_322 = buffer.data(smh + 322);
    const auto *smh_324 = buffer.data(smh + 324);
    const auto *smh_330 = buffer.data(smh + 330);
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
    const auto *smh_372 = buffer.data(smh + 372);
    const auto *smh_374 = buffer.data(smh + 374);
    const auto *smh_375 = buffer.data(smh + 375);
    const auto *smh_376 = buffer.data(smh + 376);
    const auto *smh_377 = buffer.data(smh + 377);
    const auto *smh_378 = buffer.data(smh + 378);
    const auto *smh_380 = buffer.data(smh + 380);
    const auto *smh_383 = buffer.data(smh + 383);
    const auto *smh_387 = buffer.data(smh + 387);
    const auto *smh_435 = buffer.data(smh + 435);
    const auto *smh_436 = buffer.data(smh + 436);
    const auto *smh_437 = buffer.data(smh + 437);
    const auto *smh_438 = buffer.data(smh + 438);
    const auto *smh_439 = buffer.data(smh + 439);
    const auto *smh_440 = buffer.data(smh + 440);
    const auto *smh_441 = buffer.data(smh + 441);
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
    const auto *smh_467 = buffer.data(smh + 467);
    const auto *smh_471 = buffer.data(smh + 471);
    const auto *smh_476 = buffer.data(smh + 476);
    const auto *smh_477 = buffer.data(smh + 477);
    const auto *smh_478 = buffer.data(smh + 478);
    const auto *smh_479 = buffer.data(smh + 479);
    const auto *smh_480 = buffer.data(smh + 480);
    const auto *smh_481 = buffer.data(smh + 481);
    const auto *smh_482 = buffer.data(smh + 482);
    const auto *smh_483 = buffer.data(smh + 483);
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
    const auto *smh_507 = buffer.data(smh + 507);
    const auto *smh_509 = buffer.data(smh + 509);
    const auto *smh_510 = buffer.data(smh + 510);
    const auto *smh_513 = buffer.data(smh + 513);
    const auto *smh_514 = buffer.data(smh + 514);
    const auto *smh_516 = buffer.data(smh + 516);

    const auto *smi1_420 = buffer.data(smi1 + 420);
    const auto *smi1_423 = buffer.data(smi1 + 423);
    const auto *smi1_426 = buffer.data(smi1 + 426);
    const auto *smi1_430 = buffer.data(smi1 + 430);
    const auto *smi1_432 = buffer.data(smi1 + 432);
    const auto *smi1_441 = buffer.data(smi1 + 441);

    const auto *sng0_310 = buffer.data(sng0 + 310);
    const auto *sng0_312 = buffer.data(sng0 + 312);
    const auto *sng0_313 = buffer.data(sng0 + 313);
    const auto *sng0_314 = buffer.data(sng0 + 314);
    const auto *sng0_315 = buffer.data(sng0 + 315);
    const auto *sng0_318 = buffer.data(sng0 + 318);
    const auto *sng0_320 = buffer.data(sng0 + 320);
    const auto *sng0_321 = buffer.data(sng0 + 321);
    const auto *sng0_324 = buffer.data(sng0 + 324);
    const auto *sng0_325 = buffer.data(sng0 + 325);
    const auto *sng0_327 = buffer.data(sng0 + 327);
    const auto *sng0_328 = buffer.data(sng0 + 328);
    const auto *sng0_329 = buffer.data(sng0 + 329);
    const auto *sng0_335 = buffer.data(sng0 + 335);
    const auto *sng0_339 = buffer.data(sng0 + 339);
    const auto *sng0_342 = buffer.data(sng0 + 342);
    const auto *sng0_343 = buffer.data(sng0 + 343);
    const auto *sng0_344 = buffer.data(sng0 + 344);
    const auto *sng0_345 = buffer.data(sng0 + 345);
    const auto *sng0_348 = buffer.data(sng0 + 348);
    const auto *sng0_350 = buffer.data(sng0 + 350);
    const auto *sng0_351 = buffer.data(sng0 + 351);
    const auto *sng0_354 = buffer.data(sng0 + 354);
    const auto *sng0_355 = buffer.data(sng0 + 355);
    const auto *sng0_357 = buffer.data(sng0 + 357);
    const auto *sng0_358 = buffer.data(sng0 + 358);
    const auto *sng0_359 = buffer.data(sng0 + 359);
    const auto *sng0_360 = buffer.data(sng0 + 360);
    const auto *sng0_363 = buffer.data(sng0 + 363);
    const auto *sng0_365 = buffer.data(sng0 + 365);
    const auto *sng0_366 = buffer.data(sng0 + 366);
    const auto *sng0_369 = buffer.data(sng0 + 369);
    const auto *sng0_370 = buffer.data(sng0 + 370);
    const auto *sng0_372 = buffer.data(sng0 + 372);

    const auto *sng1_310 = buffer.data(sng1 + 310);
    const auto *sng1_312 = buffer.data(sng1 + 312);
    const auto *sng1_313 = buffer.data(sng1 + 313);
    const auto *sng1_314 = buffer.data(sng1 + 314);
    const auto *sng1_315 = buffer.data(sng1 + 315);
    const auto *sng1_318 = buffer.data(sng1 + 318);
    const auto *sng1_320 = buffer.data(sng1 + 320);
    const auto *sng1_321 = buffer.data(sng1 + 321);
    const auto *sng1_324 = buffer.data(sng1 + 324);
    const auto *sng1_325 = buffer.data(sng1 + 325);
    const auto *sng1_327 = buffer.data(sng1 + 327);
    const auto *sng1_328 = buffer.data(sng1 + 328);
    const auto *sng1_329 = buffer.data(sng1 + 329);
    const auto *sng1_335 = buffer.data(sng1 + 335);
    const auto *sng1_339 = buffer.data(sng1 + 339);
    const auto *sng1_342 = buffer.data(sng1 + 342);
    const auto *sng1_343 = buffer.data(sng1 + 343);
    const auto *sng1_344 = buffer.data(sng1 + 344);
    const auto *sng1_345 = buffer.data(sng1 + 345);
    const auto *sng1_348 = buffer.data(sng1 + 348);
    const auto *sng1_350 = buffer.data(sng1 + 350);
    const auto *sng1_351 = buffer.data(sng1 + 351);
    const auto *sng1_354 = buffer.data(sng1 + 354);
    const auto *sng1_355 = buffer.data(sng1 + 355);
    const auto *sng1_357 = buffer.data(sng1 + 357);
    const auto *sng1_358 = buffer.data(sng1 + 358);
    const auto *sng1_359 = buffer.data(sng1 + 359);
    const auto *sng1_360 = buffer.data(sng1 + 360);
    const auto *sng1_363 = buffer.data(sng1 + 363);
    const auto *sng1_365 = buffer.data(sng1 + 365);
    const auto *sng1_366 = buffer.data(sng1 + 366);
    const auto *sng1_369 = buffer.data(sng1 + 369);
    const auto *sng1_370 = buffer.data(sng1 + 370);
    const auto *sng1_372 = buffer.data(sng1 + 372);

    const auto *snh_435 = buffer.data(snh + 435);
    const auto *snh_436 = buffer.data(snh + 436);
    const auto *snh_437 = buffer.data(snh + 437);
    const auto *snh_438 = buffer.data(snh + 438);
    const auto *snh_439 = buffer.data(snh + 439);
    const auto *snh_440 = buffer.data(snh + 440);
    const auto *snh_441 = buffer.data(snh + 441);
    const auto *snh_443 = buffer.data(snh + 443);
    const auto *snh_444 = buffer.data(snh + 444);
    const auto *snh_446 = buffer.data(snh + 446);
    const auto *snh_447 = buffer.data(snh + 447);
    const auto *snh_450 = buffer.data(snh + 450);
    const auto *snh_451 = buffer.data(snh + 451);
    const auto *snh_453 = buffer.data(snh + 453);
    const auto *snh_455 = buffer.data(snh + 455);
    const auto *snh_456 = buffer.data(snh + 456);
    const auto *snh_457 = buffer.data(snh + 457);
    const auto *snh_458 = buffer.data(snh + 458);
    const auto *snh_459 = buffer.data(snh + 459);
    const auto *snh_460 = buffer.data(snh + 460);
    const auto *snh_461 = buffer.data(snh + 461);
    const auto *snh_462 = buffer.data(snh + 462);
    const auto *snh_464 = buffer.data(snh + 464);
    const auto *snh_465 = buffer.data(snh + 465);
    const auto *snh_467 = buffer.data(snh + 467);
    const auto *snh_468 = buffer.data(snh + 468);
    const auto *snh_471 = buffer.data(snh + 471);
    const auto *snh_476 = buffer.data(snh + 476);
    const auto *snh_477 = buffer.data(snh + 477);
    const auto *snh_478 = buffer.data(snh + 478);
    const auto *snh_479 = buffer.data(snh + 479);
    const auto *snh_480 = buffer.data(snh + 480);
    const auto *snh_481 = buffer.data(snh + 481);
    const auto *snh_482 = buffer.data(snh + 482);
    const auto *snh_483 = buffer.data(snh + 483);
    const auto *snh_485 = buffer.data(snh + 485);
    const auto *snh_486 = buffer.data(snh + 486);
    const auto *snh_488 = buffer.data(snh + 488);
    const auto *snh_489 = buffer.data(snh + 489);
    const auto *snh_492 = buffer.data(snh + 492);
    const auto *snh_493 = buffer.data(snh + 493);
    const auto *snh_495 = buffer.data(snh + 495);
    const auto *snh_497 = buffer.data(snh + 497);
    const auto *snh_498 = buffer.data(snh + 498);
    const auto *snh_499 = buffer.data(snh + 499);
    const auto *snh_500 = buffer.data(snh + 500);
    const auto *snh_501 = buffer.data(snh + 501);
    const auto *snh_502 = buffer.data(snh + 502);
    const auto *snh_503 = buffer.data(snh + 503);
    const auto *snh_504 = buffer.data(snh + 504);
    const auto *snh_506 = buffer.data(snh + 506);
    const auto *snh_507 = buffer.data(snh + 507);
    const auto *snh_509 = buffer.data(snh + 509);
    const auto *snh_510 = buffer.data(snh + 510);
    const auto *snh_513 = buffer.data(snh + 513);
    const auto *snh_514 = buffer.data(snh + 514);
    const auto *snh_516 = buffer.data(snh + 516);

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, pc_x, smh_435, smh_436, smh_437, \
                         smh_438, smh_439, snh_435, snh_436, snh_437, snh_438, \
                         snh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_19 * smh_435[k]
                   + f_3 * pc_x[k] * snh_435[k];

        t_576[k] = f_19 * smh_436[k]
                   + f_3 * pc_x[k] * snh_436[k];

        t_577[k] = f_19 * smh_437[k]
                   + f_3 * pc_x[k] * snh_437[k];

        t_578[k] = f_19 * smh_438[k]
                   + f_3 * pc_x[k] * snh_438[k];

        t_579[k] = f_19 * smh_439[k]
                   + f_3 * pc_x[k] * snh_439[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, pc_z, smh_309, smh_440, \
                         sng0_310, sng0_312, sng1_310, sng1_312, snh_435, snh_437, \
                         snh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_19 * smh_440[k]
                   + f_3 * pc_x[k] * snh_440[k];

        t_581[k] = f_1 * sng0_310[k]
                   - f_2 * sng1_310[k]
                   + f_3 * pc_y[k] * snh_435[k];

        t_582[k] = f_19 * smh_309[k]
                   + f_3 * pc_z[k] * snh_435[k];

        t_583[k] = f_4 * sng0_312[k]
                   - f_5 * sng1_312[k]
                   + f_3 * pc_y[k] * snh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, smh_314, sng0_313, sng0_314, \
                         sng1_313, sng1_314, snh_438, snh_439, \
                         snh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * sng0_313[k]
                   - f_7 * sng1_313[k]
                   + f_3 * pc_y[k] * snh_438[k];

        t_585[k] = f_8 * sng0_314[k]
                   - f_9 * sng1_314[k]
                   + f_3 * pc_y[k] * snh_439[k];

        t_586[k] = f_3 * pc_y[k] * snh_440[k];

        t_587[k] = f_19 * smh_314[k]
                   + f_1 * sng0_314[k]
                   - f_2 * sng1_314[k]
                   + f_3 * pc_z[k] * snh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, smh_315, smh_441, \
                         smh_444, sng0_315, sng0_318, sng1_315, sng1_318, snh_441, \
                         snh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_14 * smh_441[k]
                   + f_1 * sng0_315[k]
                   - f_2 * sng1_315[k]
                   + f_3 * pc_x[k] * snh_441[k];

        t_589[k] = f_18 * smh_315[k]
                   + f_3 * pc_y[k] * snh_441[k];

        t_590[k] = f_3 * pc_z[k] * snh_441[k];

        t_591[k] = f_14 * smh_444[k]
                   + f_4 * sng0_318[k]
                   - f_5 * sng1_318[k]
                   + f_3 * pc_x[k] * snh_444[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, pc_x, pc_y, smh_317, smh_446, smh_447, sng0_320, \
                         sng0_321, sng1_320, sng1_321, snh_443, snh_446, \
                         snh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_18 * smh_317[k]
                   + f_3 * pc_y[k] * snh_443[k];

        t_593[k] = f_14 * smh_446[k]
                   + f_4 * sng0_320[k]
                   - f_5 * sng1_320[k]
                   + f_3 * pc_x[k] * snh_446[k];

        t_594[k] = f_14 * smh_447[k]
                   + f_6 * sng0_321[k]
                   - f_7 * sng1_321[k]
                   + f_3 * pc_x[k] * snh_447[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_x, pc_y, pc_z, smh_320, smh_450, sng0_324, \
                         sng1_324, snh_444, snh_446, snh_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_3 * pc_z[k] * snh_444[k];

        t_596[k] = f_18 * smh_320[k]
                   + f_3 * pc_y[k] * snh_446[k];

        t_597[k] = f_14 * smh_450[k]
                   + f_6 * sng0_324[k]
                   - f_7 * sng1_324[k]
                   + f_3 * pc_x[k] * snh_450[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pc_x, pc_z, smh_451, smh_453, sng0_325, \
                         sng0_327, sng1_325, sng1_327, snh_447, snh_451, \
                         snh_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_14 * smh_451[k]
                   + f_8 * sng0_325[k]
                   - f_9 * sng1_325[k]
                   + f_3 * pc_x[k] * snh_451[k];

        t_599[k] = f_3 * pc_z[k] * snh_447[k];

        t_600[k] = f_14 * smh_453[k]
                   + f_8 * sng0_327[k]
                   - f_9 * sng1_327[k]
                   + f_3 * pc_x[k] * snh_453[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pc_x, pc_y, smh_324, smh_455, smh_456, \
                         smh_457, sng0_329, sng1_329, snh_450, snh_455, snh_456, \
                         snh_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_18 * smh_324[k]
                   + f_3 * pc_y[k] * snh_450[k];

        t_602[k] = f_14 * smh_455[k]
                   + f_8 * sng0_329[k]
                   - f_9 * sng1_329[k]
                   + f_3 * pc_x[k] * snh_455[k];

        t_603[k] = f_14 * smh_456[k]
                   + f_3 * pc_x[k] * snh_456[k];

        t_604[k] = f_14 * smh_457[k]
                   + f_3 * pc_x[k] * snh_457[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pc_x, smh_458, smh_459, smh_460, smh_461, \
                         snh_458, snh_459, snh_460, snh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_14 * smh_458[k]
                   + f_3 * pc_x[k] * snh_458[k];

        t_606[k] = f_14 * smh_459[k]
                   + f_3 * pc_x[k] * snh_459[k];

        t_607[k] = f_14 * smh_460[k]
                   + f_3 * pc_x[k] * snh_460[k];

        t_608[k] = f_14 * smh_461[k]
                   + f_3 * pc_x[k] * snh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_y, pc_z, smh_330, smh_332, sng0_325, \
                         sng0_327, sng1_325, sng1_327, snh_456, \
                         snh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_18 * smh_330[k]
                   + f_1 * sng0_325[k]
                   - f_2 * sng1_325[k]
                   + f_3 * pc_y[k] * snh_456[k];

        t_610[k] = f_3 * pc_z[k] * snh_456[k];

        t_611[k] = f_18 * smh_332[k]
                   + f_4 * sng0_327[k]
                   - f_5 * sng1_327[k]
                   + f_3 * pc_y[k] * snh_458[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pc_y, pc_z, smh_333, smh_334, smh_335, \
                         sng0_328, sng0_329, sng1_328, sng1_329, snh_459, snh_460, \
                         snh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_18 * smh_333[k]
                   + f_6 * sng0_328[k]
                   - f_7 * sng1_328[k]
                   + f_3 * pc_y[k] * snh_459[k];

        t_613[k] = f_18 * smh_334[k]
                   + f_8 * sng0_329[k]
                   - f_9 * sng1_329[k]
                   + f_3 * pc_y[k] * snh_460[k];

        t_614[k] = f_18 * smh_335[k]
                   + f_3 * pc_y[k] * snh_461[k];

        t_615[k] = f_1 * sng0_329[k]
                   - f_2 * sng1_329[k]
                   + f_3 * pc_z[k] * snh_461[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pb_z, pc_y, pc_z, smi0_420, smi0_423, \
                         smh_315, smh_336, smi1_420, smi1_423, \
                         snh_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pb_z[k] * smi0_420[k]
                   - f_10 * pc_z[k] * smi1_420[k];

        t_617[k] = f_19 * smh_336[k]
                   + f_3 * pc_y[k] * snh_462[k];

        t_618[k] = f_11 * smh_315[k]
                   + f_3 * pc_z[k] * snh_462[k];

        t_619[k] = pb_z[k] * smi0_423[k]
                   - f_10 * pc_z[k] * smi1_423[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pb_z, pc_x, pc_y, pc_z, smi0_426, smh_338, \
                         smh_467, smi1_426, sng0_335, sng1_335, snh_464, \
                         snh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_19 * smh_338[k]
                   + f_3 * pc_y[k] * snh_464[k];

        t_621[k] = f_14 * smh_467[k]
                   + f_4 * sng0_335[k]
                   - f_5 * sng1_335[k]
                   + f_3 * pc_x[k] * snh_467[k];

        t_622[k] = pb_z[k] * smi0_426[k]
                   - f_10 * pc_z[k] * smi1_426[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, smh_318, smh_341, smh_471, \
                         sng0_339, sng1_339, snh_465, snh_467, \
                         snh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_11 * smh_318[k]
                   + f_3 * pc_z[k] * snh_465[k];

        t_624[k] = f_19 * smh_341[k]
                   + f_3 * pc_y[k] * snh_467[k];

        t_625[k] = f_14 * smh_471[k]
                   + f_6 * sng0_339[k]
                   - f_7 * sng1_339[k]
                   + f_3 * pc_x[k] * snh_471[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pb_z, pc_y, pc_z, smi0_430, smi0_432, \
                         smh_321, smh_322, smh_345, smi1_430, smi1_432, snh_468, \
                         snh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pb_z[k] * smi0_430[k]
                   - f_10 * pc_z[k] * smi1_430[k];

        t_627[k] = f_11 * smh_321[k]
                   + f_3 * pc_z[k] * snh_468[k];

        t_628[k] = pb_z[k] * smi0_432[k]
                   + f_12 * smh_322[k]
                   - f_10 * pc_z[k] * smi1_432[k];

        t_629[k] = f_19 * smh_345[k]
                   + f_3 * pc_y[k] * snh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, smh_476, smh_477, smh_478, smh_479, \
                         sng0_344, sng1_344, snh_476, snh_477, snh_478, \
                         snh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_14 * smh_476[k]
                   + f_8 * sng0_344[k]
                   - f_9 * sng1_344[k]
                   + f_3 * pc_x[k] * snh_476[k];

        t_631[k] = f_14 * smh_477[k]
                   + f_3 * pc_x[k] * snh_477[k];

        t_632[k] = f_14 * smh_478[k]
                   + f_3 * pc_x[k] * snh_478[k];

        t_633[k] = f_14 * smh_479[k]
                   + f_3 * pc_x[k] * snh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pb_z, pc_x, pc_z, smi0_441, smh_480, \
                         smh_481, smh_482, smi1_441, snh_480, snh_481, \
                         snh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_14 * smh_480[k]
                   + f_3 * pc_x[k] * snh_480[k];

        t_635[k] = f_14 * smh_481[k]
                   + f_3 * pc_x[k] * snh_481[k];

        t_636[k] = f_14 * smh_482[k]
                   + f_3 * pc_x[k] * snh_482[k];

        t_637[k] = pb_z[k] * smi0_441[k]
                   - f_10 * pc_z[k] * smi1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, smh_330, smh_353, smh_354, sng0_342, \
                         sng0_343, sng1_342, sng1_343, snh_477, snh_479, \
                         snh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * smh_330[k]
                   + f_3 * pc_z[k] * snh_477[k];

        t_639[k] = f_19 * smh_353[k]
                   + f_4 * sng0_342[k]
                   - f_5 * sng1_342[k]
                   + f_3 * pc_y[k] * snh_479[k];

        t_640[k] = f_19 * smh_354[k]
                   + f_6 * sng0_343[k]
                   - f_7 * sng1_343[k]
                   + f_3 * pc_y[k] * snh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, smh_335, smh_355, smh_356, sng0_344, \
                         sng1_344, snh_481, snh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_19 * smh_355[k]
                   + f_8 * sng0_344[k]
                   - f_9 * sng1_344[k]
                   + f_3 * pc_y[k] * snh_481[k];

        t_642[k] = f_19 * smh_356[k]
                   + f_3 * pc_y[k] * snh_482[k];

        t_643[k] = f_11 * smh_335[k]
                   + f_1 * sng0_344[k]
                   - f_2 * sng1_344[k]
                   + f_3 * pc_z[k] * snh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, smh_336, smh_357, smh_483, \
                         sng0_345, sng1_345, snh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_14 * smh_483[k]
                   + f_1 * sng0_345[k]
                   - f_2 * sng1_345[k]
                   + f_3 * pc_x[k] * snh_483[k];

        t_645[k] = f_14 * smh_357[k]
                   + f_3 * pc_y[k] * snh_483[k];

        t_646[k] = f_12 * smh_336[k]
                   + f_3 * pc_z[k] * snh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, smh_359, smh_486, smh_488, sng0_348, \
                         sng0_350, sng1_348, sng1_350, snh_485, snh_486, \
                         snh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_14 * smh_486[k]
                   + f_4 * sng0_348[k]
                   - f_5 * sng1_348[k]
                   + f_3 * pc_x[k] * snh_486[k];

        t_648[k] = f_14 * smh_359[k]
                   + f_3 * pc_y[k] * snh_485[k];

        t_649[k] = f_14 * smh_488[k]
                   + f_4 * sng0_350[k]
                   - f_5 * sng1_350[k]
                   + f_3 * pc_x[k] * snh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, smh_339, smh_362, smh_489, \
                         sng0_351, sng1_351, snh_486, snh_488, \
                         snh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_14 * smh_489[k]
                   + f_6 * sng0_351[k]
                   - f_7 * sng1_351[k]
                   + f_3 * pc_x[k] * snh_489[k];

        t_651[k] = f_12 * smh_339[k]
                   + f_3 * pc_z[k] * snh_486[k];

        t_652[k] = f_14 * smh_362[k]
                   + f_3 * pc_y[k] * snh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, smh_342, smh_492, smh_493, sng0_354, \
                         sng0_355, sng1_354, sng1_355, snh_489, snh_492, \
                         snh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_14 * smh_492[k]
                   + f_6 * sng0_354[k]
                   - f_7 * sng1_354[k]
                   + f_3 * pc_x[k] * snh_492[k];

        t_654[k] = f_14 * smh_493[k]
                   + f_8 * sng0_355[k]
                   - f_9 * sng1_355[k]
                   + f_3 * pc_x[k] * snh_493[k];

        t_655[k] = f_12 * smh_342[k]
                   + f_3 * pc_z[k] * snh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, smh_366, smh_495, smh_497, sng0_357, \
                         sng0_359, sng1_357, sng1_359, snh_492, snh_495, \
                         snh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * smh_495[k]
                   + f_8 * sng0_357[k]
                   - f_9 * sng1_357[k]
                   + f_3 * pc_x[k] * snh_495[k];

        t_657[k] = f_14 * smh_366[k]
                   + f_3 * pc_y[k] * snh_492[k];

        t_658[k] = f_14 * smh_497[k]
                   + f_8 * sng0_359[k]
                   - f_9 * sng1_359[k]
                   + f_3 * pc_x[k] * snh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, smh_498, smh_499, smh_500, \
                         smh_501, smh_502, snh_498, snh_499, snh_500, snh_501, \
                         snh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_14 * smh_498[k]
                   + f_3 * pc_x[k] * snh_498[k];

        t_660[k] = f_14 * smh_499[k]
                   + f_3 * pc_x[k] * snh_499[k];

        t_661[k] = f_14 * smh_500[k]
                   + f_3 * pc_x[k] * snh_500[k];

        t_662[k] = f_14 * smh_501[k]
                   + f_3 * pc_x[k] * snh_501[k];

        t_663[k] = f_14 * smh_502[k]
                   + f_3 * pc_x[k] * snh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, smh_351, smh_372, smh_503, \
                         sng0_355, sng1_355, snh_498, snh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_14 * smh_503[k]
                   + f_3 * pc_x[k] * snh_503[k];

        t_665[k] = f_14 * smh_372[k]
                   + f_1 * sng0_355[k]
                   - f_2 * sng1_355[k]
                   + f_3 * pc_y[k] * snh_498[k];

        t_666[k] = f_12 * smh_351[k]
                   + f_3 * pc_z[k] * snh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, smh_374, smh_375, smh_376, sng0_357, \
                         sng0_358, sng0_359, sng1_357, sng1_358, sng1_359, snh_500, snh_501, \
                         snh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * smh_374[k]
                   + f_4 * sng0_357[k]
                   - f_5 * sng1_357[k]
                   + f_3 * pc_y[k] * snh_500[k];

        t_668[k] = f_14 * smh_375[k]
                   + f_6 * sng0_358[k]
                   - f_7 * sng1_358[k]
                   + f_3 * pc_y[k] * snh_501[k];

        t_669[k] = f_14 * smh_376[k]
                   + f_8 * sng0_359[k]
                   - f_9 * sng1_359[k]
                   + f_3 * pc_y[k] * snh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, smh_356, smh_377, smh_504, \
                         sng0_359, sng0_360, sng1_359, sng1_360, snh_503, \
                         snh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * smh_377[k]
                   + f_3 * pc_y[k] * snh_503[k];

        t_671[k] = f_12 * smh_356[k]
                   + f_1 * sng0_359[k]
                   - f_2 * sng1_359[k]
                   + f_3 * pc_z[k] * snh_503[k];

        t_672[k] = f_14 * smh_504[k]
                   + f_1 * sng0_360[k]
                   - f_2 * sng1_360[k]
                   + f_3 * pc_x[k] * snh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, smh_357, smh_378, \
                         smh_380, smh_507, sng0_363, sng1_363, snh_504, snh_506, \
                         snh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * smh_378[k]
                   + f_3 * pc_y[k] * snh_504[k];

        t_674[k] = f_13 * smh_357[k]
                   + f_3 * pc_z[k] * snh_504[k];

        t_675[k] = f_14 * smh_507[k]
                   + f_4 * sng0_363[k]
                   - f_5 * sng1_363[k]
                   + f_3 * pc_x[k] * snh_507[k];

        t_676[k] = f_13 * smh_380[k]
                   + f_3 * pc_y[k] * snh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, smh_360, smh_509, smh_510, sng0_365, \
                         sng0_366, sng1_365, sng1_366, snh_507, snh_509, \
                         snh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_14 * smh_509[k]
                   + f_4 * sng0_365[k]
                   - f_5 * sng1_365[k]
                   + f_3 * pc_x[k] * snh_509[k];

        t_678[k] = f_14 * smh_510[k]
                   + f_6 * sng0_366[k]
                   - f_7 * sng1_366[k]
                   + f_3 * pc_x[k] * snh_510[k];

        t_679[k] = f_13 * smh_360[k]
                   + f_3 * pc_z[k] * snh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, smh_383, smh_513, smh_514, sng0_369, \
                         sng0_370, sng1_369, sng1_370, snh_509, snh_513, \
                         snh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * smh_383[k]
                   + f_3 * pc_y[k] * snh_509[k];

        t_681[k] = f_14 * smh_513[k]
                   + f_6 * sng0_369[k]
                   - f_7 * sng1_369[k]
                   + f_3 * pc_x[k] * snh_513[k];

        t_682[k] = f_14 * smh_514[k]
                   + f_8 * sng0_370[k]
                   - f_9 * sng1_370[k]
                   + f_3 * pc_x[k] * snh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, smh_363, smh_387, smh_516, \
                         sng0_372, sng1_372, snh_510, snh_513, \
                         snh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * smh_363[k]
                   + f_3 * pc_z[k] * snh_510[k];

        t_684[k] = f_14 * smh_516[k]
                   + f_8 * sng0_372[k]
                   - f_9 * sng1_372[k]
                   + f_3 * pc_x[k] * snh_516[k];

        t_685[k] = f_13 * smh_387[k]
                   + f_3 * pc_y[k] * snh_513[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *smi0_560 = buffer.data(smi0 + 560);
    const auto *smi0_563 = buffer.data(smi0 + 563);
    const auto *smi0_565 = buffer.data(smi0 + 565);
    const auto *smi0_566 = buffer.data(smi0 + 566);
    const auto *smi0_569 = buffer.data(smi0 + 569);
    const auto *smi0_570 = buffer.data(smi0 + 570);
    const auto *smi0_572 = buffer.data(smi0 + 572);
    const auto *smi0_574 = buffer.data(smi0 + 574);
    const auto *smi0_587 = buffer.data(smi0 + 587);

    const auto *smh_372 = buffer.data(smh + 372);
    const auto *smh_377 = buffer.data(smh + 377);
    const auto *smh_378 = buffer.data(smh + 378);
    const auto *smh_381 = buffer.data(smh + 381);
    const auto *smh_384 = buffer.data(smh + 384);
    const auto *smh_393 = buffer.data(smh + 393);
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
    const auto *smh_416 = buffer.data(smh + 416);
    const auto *smh_417 = buffer.data(smh + 417);
    const auto *smh_418 = buffer.data(smh + 418);
    const auto *smh_419 = buffer.data(smh + 419);
    const auto *smh_420 = buffer.data(smh + 420);
    const auto *smh_421 = buffer.data(smh + 421);
    const auto *smh_422 = buffer.data(smh + 422);
    const auto *smh_423 = buffer.data(smh + 423);
    const auto *smh_425 = buffer.data(smh + 425);
    const auto *smh_426 = buffer.data(smh + 426);
    const auto *smh_428 = buffer.data(smh + 428);
    const auto *smh_429 = buffer.data(smh + 429);
    const auto *smh_435 = buffer.data(smh + 435);
    const auto *smh_437 = buffer.data(smh + 437);
    const auto *smh_438 = buffer.data(smh + 438);
    const auto *smh_439 = buffer.data(smh + 439);
    const auto *smh_440 = buffer.data(smh + 440);
    const auto *smh_441 = buffer.data(smh + 441);
    const auto *smh_443 = buffer.data(smh + 443);
    const auto *smh_446 = buffer.data(smh + 446);
    const auto *smh_518 = buffer.data(smh + 518);
    const auto *smh_519 = buffer.data(smh + 519);
    const auto *smh_520 = buffer.data(smh + 520);
    const auto *smh_521 = buffer.data(smh + 521);
    const auto *smh_522 = buffer.data(smh + 522);
    const auto *smh_523 = buffer.data(smh + 523);
    const auto *smh_524 = buffer.data(smh + 524);
    const auto *smh_525 = buffer.data(smh + 525);
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
    const auto *smh_561 = buffer.data(smh + 561);
    const auto *smh_562 = buffer.data(smh + 562);
    const auto *smh_563 = buffer.data(smh + 563);
    const auto *smh_564 = buffer.data(smh + 564);
    const auto *smh_565 = buffer.data(smh + 565);
    const auto *smh_566 = buffer.data(smh + 566);
    const auto *smh_567 = buffer.data(smh + 567);
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
    const auto *smh_591 = buffer.data(smh + 591);
    const auto *smh_593 = buffer.data(smh + 593);
    const auto *smh_594 = buffer.data(smh + 594);
    const auto *smh_597 = buffer.data(smh + 597);
    const auto *smh_598 = buffer.data(smh + 598);
    const auto *smh_600 = buffer.data(smh + 600);

    const auto *smi1_560 = buffer.data(smi1 + 560);
    const auto *smi1_563 = buffer.data(smi1 + 563);
    const auto *smi1_565 = buffer.data(smi1 + 565);
    const auto *smi1_566 = buffer.data(smi1 + 566);
    const auto *smi1_569 = buffer.data(smi1 + 569);
    const auto *smi1_570 = buffer.data(smi1 + 570);
    const auto *smi1_572 = buffer.data(smi1 + 572);
    const auto *smi1_574 = buffer.data(smi1 + 574);
    const auto *smi1_587 = buffer.data(smi1 + 587);

    const auto *sng0_370 = buffer.data(sng0 + 370);
    const auto *sng0_372 = buffer.data(sng0 + 372);
    const auto *sng0_373 = buffer.data(sng0 + 373);
    const auto *sng0_374 = buffer.data(sng0 + 374);
    const auto *sng0_375 = buffer.data(sng0 + 375);
    const auto *sng0_378 = buffer.data(sng0 + 378);
    const auto *sng0_380 = buffer.data(sng0 + 380);
    const auto *sng0_381 = buffer.data(sng0 + 381);
    const auto *sng0_384 = buffer.data(sng0 + 384);
    const auto *sng0_385 = buffer.data(sng0 + 385);
    const auto *sng0_387 = buffer.data(sng0 + 387);
    const auto *sng0_388 = buffer.data(sng0 + 388);
    const auto *sng0_389 = buffer.data(sng0 + 389);
    const auto *sng0_400 = buffer.data(sng0 + 400);
    const auto *sng0_402 = buffer.data(sng0 + 402);
    const auto *sng0_403 = buffer.data(sng0 + 403);
    const auto *sng0_404 = buffer.data(sng0 + 404);
    const auto *sng0_405 = buffer.data(sng0 + 405);
    const auto *sng0_408 = buffer.data(sng0 + 408);
    const auto *sng0_410 = buffer.data(sng0 + 410);
    const auto *sng0_411 = buffer.data(sng0 + 411);
    const auto *sng0_414 = buffer.data(sng0 + 414);
    const auto *sng0_415 = buffer.data(sng0 + 415);
    const auto *sng0_417 = buffer.data(sng0 + 417);
    const auto *sng0_418 = buffer.data(sng0 + 418);
    const auto *sng0_419 = buffer.data(sng0 + 419);
    const auto *sng0_420 = buffer.data(sng0 + 420);
    const auto *sng0_423 = buffer.data(sng0 + 423);
    const auto *sng0_425 = buffer.data(sng0 + 425);
    const auto *sng0_426 = buffer.data(sng0 + 426);
    const auto *sng0_429 = buffer.data(sng0 + 429);
    const auto *sng0_430 = buffer.data(sng0 + 430);
    const auto *sng0_432 = buffer.data(sng0 + 432);

    const auto *sng1_370 = buffer.data(sng1 + 370);
    const auto *sng1_372 = buffer.data(sng1 + 372);
    const auto *sng1_373 = buffer.data(sng1 + 373);
    const auto *sng1_374 = buffer.data(sng1 + 374);
    const auto *sng1_375 = buffer.data(sng1 + 375);
    const auto *sng1_378 = buffer.data(sng1 + 378);
    const auto *sng1_380 = buffer.data(sng1 + 380);
    const auto *sng1_381 = buffer.data(sng1 + 381);
    const auto *sng1_384 = buffer.data(sng1 + 384);
    const auto *sng1_385 = buffer.data(sng1 + 385);
    const auto *sng1_387 = buffer.data(sng1 + 387);
    const auto *sng1_388 = buffer.data(sng1 + 388);
    const auto *sng1_389 = buffer.data(sng1 + 389);
    const auto *sng1_400 = buffer.data(sng1 + 400);
    const auto *sng1_402 = buffer.data(sng1 + 402);
    const auto *sng1_403 = buffer.data(sng1 + 403);
    const auto *sng1_404 = buffer.data(sng1 + 404);
    const auto *sng1_405 = buffer.data(sng1 + 405);
    const auto *sng1_408 = buffer.data(sng1 + 408);
    const auto *sng1_410 = buffer.data(sng1 + 410);
    const auto *sng1_411 = buffer.data(sng1 + 411);
    const auto *sng1_414 = buffer.data(sng1 + 414);
    const auto *sng1_415 = buffer.data(sng1 + 415);
    const auto *sng1_417 = buffer.data(sng1 + 417);
    const auto *sng1_418 = buffer.data(sng1 + 418);
    const auto *sng1_419 = buffer.data(sng1 + 419);
    const auto *sng1_420 = buffer.data(sng1 + 420);
    const auto *sng1_423 = buffer.data(sng1 + 423);
    const auto *sng1_425 = buffer.data(sng1 + 425);
    const auto *sng1_426 = buffer.data(sng1 + 426);
    const auto *sng1_429 = buffer.data(sng1 + 429);
    const auto *sng1_430 = buffer.data(sng1 + 430);
    const auto *sng1_432 = buffer.data(sng1 + 432);

    const auto *snh_518 = buffer.data(snh + 518);
    const auto *snh_519 = buffer.data(snh + 519);
    const auto *snh_520 = buffer.data(snh + 520);
    const auto *snh_521 = buffer.data(snh + 521);
    const auto *snh_522 = buffer.data(snh + 522);
    const auto *snh_523 = buffer.data(snh + 523);
    const auto *snh_524 = buffer.data(snh + 524);
    const auto *snh_525 = buffer.data(snh + 525);
    const auto *snh_527 = buffer.data(snh + 527);
    const auto *snh_528 = buffer.data(snh + 528);
    const auto *snh_530 = buffer.data(snh + 530);
    const auto *snh_531 = buffer.data(snh + 531);
    const auto *snh_534 = buffer.data(snh + 534);
    const auto *snh_535 = buffer.data(snh + 535);
    const auto *snh_537 = buffer.data(snh + 537);
    const auto *snh_539 = buffer.data(snh + 539);
    const auto *snh_540 = buffer.data(snh + 540);
    const auto *snh_541 = buffer.data(snh + 541);
    const auto *snh_542 = buffer.data(snh + 542);
    const auto *snh_543 = buffer.data(snh + 543);
    const auto *snh_544 = buffer.data(snh + 544);
    const auto *snh_545 = buffer.data(snh + 545);
    const auto *snh_546 = buffer.data(snh + 546);
    const auto *snh_548 = buffer.data(snh + 548);
    const auto *snh_549 = buffer.data(snh + 549);
    const auto *snh_551 = buffer.data(snh + 551);
    const auto *snh_552 = buffer.data(snh + 552);
    const auto *snh_555 = buffer.data(snh + 555);
    const auto *snh_561 = buffer.data(snh + 561);
    const auto *snh_562 = buffer.data(snh + 562);
    const auto *snh_563 = buffer.data(snh + 563);
    const auto *snh_564 = buffer.data(snh + 564);
    const auto *snh_565 = buffer.data(snh + 565);
    const auto *snh_566 = buffer.data(snh + 566);
    const auto *snh_567 = buffer.data(snh + 567);
    const auto *snh_569 = buffer.data(snh + 569);
    const auto *snh_570 = buffer.data(snh + 570);
    const auto *snh_572 = buffer.data(snh + 572);
    const auto *snh_573 = buffer.data(snh + 573);
    const auto *snh_576 = buffer.data(snh + 576);
    const auto *snh_577 = buffer.data(snh + 577);
    const auto *snh_579 = buffer.data(snh + 579);
    const auto *snh_581 = buffer.data(snh + 581);
    const auto *snh_582 = buffer.data(snh + 582);
    const auto *snh_583 = buffer.data(snh + 583);
    const auto *snh_584 = buffer.data(snh + 584);
    const auto *snh_585 = buffer.data(snh + 585);
    const auto *snh_586 = buffer.data(snh + 586);
    const auto *snh_587 = buffer.data(snh + 587);
    const auto *snh_588 = buffer.data(snh + 588);
    const auto *snh_590 = buffer.data(snh + 590);
    const auto *snh_591 = buffer.data(snh + 591);
    const auto *snh_593 = buffer.data(snh + 593);
    const auto *snh_594 = buffer.data(snh + 594);
    const auto *snh_597 = buffer.data(snh + 597);
    const auto *snh_598 = buffer.data(snh + 598);
    const auto *snh_600 = buffer.data(snh + 600);

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, smh_518, smh_519, smh_520, smh_521, \
                         sng0_374, sng1_374, snh_518, snh_519, snh_520, \
                         snh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_14 * smh_518[k]
                   + f_8 * sng0_374[k]
                   - f_9 * sng1_374[k]
                   + f_3 * pc_x[k] * snh_518[k];

        t_687[k] = f_14 * smh_519[k]
                   + f_3 * pc_x[k] * snh_519[k];

        t_688[k] = f_14 * smh_520[k]
                   + f_3 * pc_x[k] * snh_520[k];

        t_689[k] = f_14 * smh_521[k]
                   + f_3 * pc_x[k] * snh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, smh_393, smh_522, smh_523, \
                         smh_524, sng0_370, sng1_370, snh_519, snh_522, snh_523, \
                         snh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_14 * smh_522[k]
                   + f_3 * pc_x[k] * snh_522[k];

        t_691[k] = f_14 * smh_523[k]
                   + f_3 * pc_x[k] * snh_523[k];

        t_692[k] = f_14 * smh_524[k]
                   + f_3 * pc_x[k] * snh_524[k];

        t_693[k] = f_13 * smh_393[k]
                   + f_1 * sng0_370[k]
                   - f_2 * sng1_370[k]
                   + f_3 * pc_y[k] * snh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, smh_372, smh_395, smh_396, sng0_372, \
                         sng0_373, sng1_372, sng1_373, snh_519, snh_521, \
                         snh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * smh_372[k]
                   + f_3 * pc_z[k] * snh_519[k];

        t_695[k] = f_13 * smh_395[k]
                   + f_4 * sng0_372[k]
                   - f_5 * sng1_372[k]
                   + f_3 * pc_y[k] * snh_521[k];

        t_696[k] = f_13 * smh_396[k]
                   + f_6 * sng0_373[k]
                   - f_7 * sng1_373[k]
                   + f_3 * pc_y[k] * snh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, smh_377, smh_397, smh_398, sng0_374, \
                         sng1_374, snh_523, snh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * smh_397[k]
                   + f_8 * sng0_374[k]
                   - f_9 * sng1_374[k]
                   + f_3 * pc_y[k] * snh_523[k];

        t_698[k] = f_13 * smh_398[k]
                   + f_3 * pc_y[k] * snh_524[k];

        t_699[k] = f_13 * smh_377[k]
                   + f_1 * sng0_374[k]
                   - f_2 * sng1_374[k]
                   + f_3 * pc_z[k] * snh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, smh_378, smh_399, smh_525, \
                         sng0_375, sng1_375, snh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_14 * smh_525[k]
                   + f_1 * sng0_375[k]
                   - f_2 * sng1_375[k]
                   + f_3 * pc_x[k] * snh_525[k];

        t_701[k] = f_12 * smh_399[k]
                   + f_3 * pc_y[k] * snh_525[k];

        t_702[k] = f_14 * smh_378[k]
                   + f_3 * pc_z[k] * snh_525[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, smh_401, smh_528, smh_530, sng0_378, \
                         sng0_380, sng1_378, sng1_380, snh_527, snh_528, \
                         snh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_14 * smh_528[k]
                   + f_4 * sng0_378[k]
                   - f_5 * sng1_378[k]
                   + f_3 * pc_x[k] * snh_528[k];

        t_704[k] = f_12 * smh_401[k]
                   + f_3 * pc_y[k] * snh_527[k];

        t_705[k] = f_14 * smh_530[k]
                   + f_4 * sng0_380[k]
                   - f_5 * sng1_380[k]
                   + f_3 * pc_x[k] * snh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, smh_381, smh_404, smh_531, \
                         sng0_381, sng1_381, snh_528, snh_530, \
                         snh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_14 * smh_531[k]
                   + f_6 * sng0_381[k]
                   - f_7 * sng1_381[k]
                   + f_3 * pc_x[k] * snh_531[k];

        t_707[k] = f_14 * smh_381[k]
                   + f_3 * pc_z[k] * snh_528[k];

        t_708[k] = f_12 * smh_404[k]
                   + f_3 * pc_y[k] * snh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, smh_384, smh_534, smh_535, sng0_384, \
                         sng0_385, sng1_384, sng1_385, snh_531, snh_534, \
                         snh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_14 * smh_534[k]
                   + f_6 * sng0_384[k]
                   - f_7 * sng1_384[k]
                   + f_3 * pc_x[k] * snh_534[k];

        t_710[k] = f_14 * smh_535[k]
                   + f_8 * sng0_385[k]
                   - f_9 * sng1_385[k]
                   + f_3 * pc_x[k] * snh_535[k];

        t_711[k] = f_14 * smh_384[k]
                   + f_3 * pc_z[k] * snh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, smh_408, smh_537, smh_539, sng0_387, \
                         sng0_389, sng1_387, sng1_389, snh_534, snh_537, \
                         snh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_14 * smh_537[k]
                   + f_8 * sng0_387[k]
                   - f_9 * sng1_387[k]
                   + f_3 * pc_x[k] * snh_537[k];

        t_713[k] = f_12 * smh_408[k]
                   + f_3 * pc_y[k] * snh_534[k];

        t_714[k] = f_14 * smh_539[k]
                   + f_8 * sng0_389[k]
                   - f_9 * sng1_389[k]
                   + f_3 * pc_x[k] * snh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, smh_540, smh_541, smh_542, \
                         smh_543, smh_544, snh_540, snh_541, snh_542, snh_543, \
                         snh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_14 * smh_540[k]
                   + f_3 * pc_x[k] * snh_540[k];

        t_716[k] = f_14 * smh_541[k]
                   + f_3 * pc_x[k] * snh_541[k];

        t_717[k] = f_14 * smh_542[k]
                   + f_3 * pc_x[k] * snh_542[k];

        t_718[k] = f_14 * smh_543[k]
                   + f_3 * pc_x[k] * snh_543[k];

        t_719[k] = f_14 * smh_544[k]
                   + f_3 * pc_x[k] * snh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, smh_393, smh_414, smh_545, \
                         sng0_385, sng1_385, snh_540, snh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_14 * smh_545[k]
                   + f_3 * pc_x[k] * snh_545[k];

        t_721[k] = f_12 * smh_414[k]
                   + f_1 * sng0_385[k]
                   - f_2 * sng1_385[k]
                   + f_3 * pc_y[k] * snh_540[k];

        t_722[k] = f_14 * smh_393[k]
                   + f_3 * pc_z[k] * snh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, smh_416, smh_417, smh_418, sng0_387, \
                         sng0_388, sng0_389, sng1_387, sng1_388, sng1_389, snh_542, snh_543, \
                         snh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * smh_416[k]
                   + f_4 * sng0_387[k]
                   - f_5 * sng1_387[k]
                   + f_3 * pc_y[k] * snh_542[k];

        t_724[k] = f_12 * smh_417[k]
                   + f_6 * sng0_388[k]
                   - f_7 * sng1_388[k]
                   + f_3 * pc_y[k] * snh_543[k];

        t_725[k] = f_12 * smh_418[k]
                   + f_8 * sng0_389[k]
                   - f_9 * sng1_389[k]
                   + f_3 * pc_y[k] * snh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pb_y, pc_y, pc_z, smi0_560, smh_398, \
                         smh_419, smh_420, smi1_560, sng0_389, sng1_389, snh_545, \
                         snh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * smh_419[k]
                   + f_3 * pc_y[k] * snh_545[k];

        t_727[k] = f_14 * smh_398[k]
                   + f_1 * sng0_389[k]
                   - f_2 * sng1_389[k]
                   + f_3 * pc_z[k] * snh_545[k];

        t_728[k] = pb_y[k] * smi0_560[k]
                   - f_10 * pc_y[k] * smi1_560[k];

        t_729[k] = f_11 * smh_420[k]
                   + f_3 * pc_y[k] * snh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pb_y, pc_y, pc_z, smi0_563, smi0_565, \
                         smh_399, smh_421, smh_422, smi1_563, smi1_565, snh_546, \
                         snh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_19 * smh_399[k]
                   + f_3 * pc_z[k] * snh_546[k];

        t_731[k] = pb_y[k] * smi0_563[k]
                   + f_12 * smh_421[k]
                   - f_10 * pc_y[k] * smi1_563[k];

        t_732[k] = f_11 * smh_422[k]
                   + f_3 * pc_y[k] * snh_548[k];

        t_733[k] = pb_y[k] * smi0_565[k]
                   - f_10 * pc_y[k] * smi1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pb_y, pc_y, pc_z, smi0_566, smi0_569, \
                         smh_402, smh_423, smh_425, smi1_566, smi1_569, snh_549, \
                         snh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pb_y[k] * smi0_566[k]
                   + f_13 * smh_423[k]
                   - f_10 * pc_y[k] * smi1_566[k];

        t_735[k] = f_19 * smh_402[k]
                   + f_3 * pc_z[k] * snh_549[k];

        t_736[k] = f_11 * smh_425[k]
                   + f_3 * pc_y[k] * snh_551[k];

        t_737[k] = pb_y[k] * smi0_569[k]
                   - f_10 * pc_y[k] * smi1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pb_y, pc_y, pc_z, smi0_570, smi0_572, smh_405, \
                         smh_426, smh_428, smi1_570, smi1_572, \
                         snh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pb_y[k] * smi0_570[k]
                   + f_14 * smh_426[k]
                   - f_10 * pc_y[k] * smi1_570[k];

        t_739[k] = f_19 * smh_405[k]
                   + f_3 * pc_z[k] * snh_552[k];

        t_740[k] = pb_y[k] * smi0_572[k]
                   + f_12 * smh_428[k]
                   - f_10 * pc_y[k] * smi1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pb_y, pc_x, pc_y, smi0_574, smh_429, \
                         smh_561, smh_562, smi1_574, snh_555, snh_561, \
                         snh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * smh_429[k]
                   + f_3 * pc_y[k] * snh_555[k];

        t_742[k] = pb_y[k] * smi0_574[k]
                   - f_10 * pc_y[k] * smi1_574[k];

        t_743[k] = f_14 * smh_561[k]
                   + f_3 * pc_x[k] * snh_561[k];

        t_744[k] = f_14 * smh_562[k]
                   + f_3 * pc_x[k] * snh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, smh_563, smh_564, smh_565, smh_566, \
                         snh_563, snh_564, snh_565, snh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_14 * smh_563[k]
                   + f_3 * pc_x[k] * snh_563[k];

        t_746[k] = f_14 * smh_564[k]
                   + f_3 * pc_x[k] * snh_564[k];

        t_747[k] = f_14 * smh_565[k]
                   + f_3 * pc_x[k] * snh_565[k];

        t_748[k] = f_14 * smh_566[k]
                   + f_3 * pc_x[k] * snh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, smh_414, smh_435, smh_437, sng0_400, \
                         sng0_402, sng1_400, sng1_402, snh_561, \
                         snh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * smh_435[k]
                   + f_1 * sng0_400[k]
                   - f_2 * sng1_400[k]
                   + f_3 * pc_y[k] * snh_561[k];

        t_750[k] = f_19 * smh_414[k]
                   + f_3 * pc_z[k] * snh_561[k];

        t_751[k] = f_11 * smh_437[k]
                   + f_4 * sng0_402[k]
                   - f_5 * sng1_402[k]
                   + f_3 * pc_y[k] * snh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, smh_438, smh_439, smh_440, sng0_403, \
                         sng0_404, sng1_403, sng1_404, snh_564, snh_565, \
                         snh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * smh_438[k]
                   + f_6 * sng0_403[k]
                   - f_7 * sng1_403[k]
                   + f_3 * pc_y[k] * snh_564[k];

        t_753[k] = f_11 * smh_439[k]
                   + f_8 * sng0_404[k]
                   - f_9 * sng1_404[k]
                   + f_3 * pc_y[k] * snh_565[k];

        t_754[k] = f_11 * smh_440[k]
                   + f_3 * pc_y[k] * snh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pb_y, pc_x, pc_y, pc_z, smi0_587, \
                         smh_420, smh_567, smi1_587, sng0_405, sng1_405, \
                         snh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pb_y[k] * smi0_587[k]
                   - f_10 * pc_y[k] * smi1_587[k];

        t_756[k] = f_14 * smh_567[k]
                   + f_1 * sng0_405[k]
                   - f_2 * sng1_405[k]
                   + f_3 * pc_x[k] * snh_567[k];

        t_757[k] = f_3 * pc_y[k] * snh_567[k];

        t_758[k] = f_18 * smh_420[k]
                   + f_3 * pc_z[k] * snh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, smh_570, smh_572, sng0_408, \
                         sng0_410, sng1_408, sng1_410, snh_569, snh_570, \
                         snh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_14 * smh_570[k]
                   + f_4 * sng0_408[k]
                   - f_5 * sng1_408[k]
                   + f_3 * pc_x[k] * snh_570[k];

        t_760[k] = f_3 * pc_y[k] * snh_569[k];

        t_761[k] = f_14 * smh_572[k]
                   + f_4 * sng0_410[k]
                   - f_5 * sng1_410[k]
                   + f_3 * pc_x[k] * snh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_x, pc_y, pc_z, smh_423, smh_573, sng0_411, \
                         sng1_411, snh_570, snh_572, snh_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_14 * smh_573[k]
                   + f_6 * sng0_411[k]
                   - f_7 * sng1_411[k]
                   + f_3 * pc_x[k] * snh_573[k];

        t_763[k] = f_18 * smh_423[k]
                   + f_3 * pc_z[k] * snh_570[k];

        t_764[k] = f_3 * pc_y[k] * snh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_z, smh_426, smh_576, smh_577, sng0_414, \
                         sng0_415, sng1_414, sng1_415, snh_573, snh_576, \
                         snh_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_14 * smh_576[k]
                   + f_6 * sng0_414[k]
                   - f_7 * sng1_414[k]
                   + f_3 * pc_x[k] * snh_576[k];

        t_766[k] = f_14 * smh_577[k]
                   + f_8 * sng0_415[k]
                   - f_9 * sng1_415[k]
                   + f_3 * pc_x[k] * snh_577[k];

        t_767[k] = f_18 * smh_426[k]
                   + f_3 * pc_z[k] * snh_573[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_x, pc_y, smh_579, smh_581, sng0_417, \
                         sng0_419, sng1_417, sng1_419, snh_576, snh_579, \
                         snh_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_14 * smh_579[k]
                   + f_8 * sng0_417[k]
                   - f_9 * sng1_417[k]
                   + f_3 * pc_x[k] * snh_579[k];

        t_769[k] = f_3 * pc_y[k] * snh_576[k];

        t_770[k] = f_14 * smh_581[k]
                   + f_8 * sng0_419[k]
                   - f_9 * sng1_419[k]
                   + f_3 * pc_x[k] * snh_581[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, t_775, pc_x, smh_582, smh_583, smh_584, \
                         smh_585, smh_586, snh_582, snh_583, snh_584, snh_585, \
                         snh_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_14 * smh_582[k]
                   + f_3 * pc_x[k] * snh_582[k];

        t_772[k] = f_14 * smh_583[k]
                   + f_3 * pc_x[k] * snh_583[k];

        t_773[k] = f_14 * smh_584[k]
                   + f_3 * pc_x[k] * snh_584[k];

        t_774[k] = f_14 * smh_585[k]
                   + f_3 * pc_x[k] * snh_585[k];

        t_775[k] = f_14 * smh_586[k]
                   + f_3 * pc_x[k] * snh_586[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pc_x, pc_y, pc_z, smh_435, smh_587, \
                         sng0_415, sng0_417, sng1_415, sng1_417, snh_582, snh_584, \
                         snh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_14 * smh_587[k]
                   + f_3 * pc_x[k] * snh_587[k];

        t_777[k] = f_1 * sng0_415[k]
                   - f_2 * sng1_415[k]
                   + f_3 * pc_y[k] * snh_582[k];

        t_778[k] = f_18 * smh_435[k]
                   + f_3 * pc_z[k] * snh_582[k];

        t_779[k] = f_4 * sng0_417[k]
                   - f_5 * sng1_417[k]
                   + f_3 * pc_y[k] * snh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, smh_440, sng0_418, sng0_419, \
                         sng1_418, sng1_419, snh_585, snh_586, \
                         snh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * sng0_418[k]
                   - f_7 * sng1_418[k]
                   + f_3 * pc_y[k] * snh_585[k];

        t_781[k] = f_8 * sng0_419[k]
                   - f_9 * sng1_419[k]
                   + f_3 * pc_y[k] * snh_586[k];

        t_782[k] = f_3 * pc_y[k] * snh_587[k];

        t_783[k] = f_18 * smh_440[k]
                   + f_1 * sng0_419[k]
                   - f_2 * sng1_419[k]
                   + f_3 * pc_z[k] * snh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pc_x, pc_y, pc_z, smh_441, smh_588, \
                         smh_591, sng0_420, sng0_423, sng1_420, sng1_423, snh_588, \
                         snh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_13 * smh_588[k]
                   + f_1 * sng0_420[k]
                   - f_2 * sng1_420[k]
                   + f_3 * pc_x[k] * snh_588[k];

        t_785[k] = f_17 * smh_441[k]
                   + f_3 * pc_y[k] * snh_588[k];

        t_786[k] = f_3 * pc_z[k] * snh_588[k];

        t_787[k] = f_13 * smh_591[k]
                   + f_4 * sng0_423[k]
                   - f_5 * sng1_423[k]
                   + f_3 * pc_x[k] * snh_591[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, pc_y, smh_443, smh_593, smh_594, sng0_425, \
                         sng0_426, sng1_425, sng1_426, snh_590, snh_593, \
                         snh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_17 * smh_443[k]
                   + f_3 * pc_y[k] * snh_590[k];

        t_789[k] = f_13 * smh_593[k]
                   + f_4 * sng0_425[k]
                   - f_5 * sng1_425[k]
                   + f_3 * pc_x[k] * snh_593[k];

        t_790[k] = f_13 * smh_594[k]
                   + f_6 * sng0_426[k]
                   - f_7 * sng1_426[k]
                   + f_3 * pc_x[k] * snh_594[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, pc_x, pc_y, pc_z, smh_446, smh_597, sng0_429, \
                         sng1_429, snh_591, snh_593, snh_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_3 * pc_z[k] * snh_591[k];

        t_792[k] = f_17 * smh_446[k]
                   + f_3 * pc_y[k] * snh_593[k];

        t_793[k] = f_13 * smh_597[k]
                   + f_6 * sng0_429[k]
                   - f_7 * sng1_429[k]
                   + f_3 * pc_x[k] * snh_597[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pc_x, pc_z, smh_598, smh_600, sng0_430, \
                         sng0_432, sng1_430, sng1_432, snh_594, snh_598, \
                         snh_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_13 * smh_598[k]
                   + f_8 * sng0_430[k]
                   - f_9 * sng1_430[k]
                   + f_3 * pc_x[k] * snh_598[k];

        t_795[k] = f_3 * pc_z[k] * snh_594[k];

        t_796[k] = f_13 * smh_600[k]
                   + f_8 * sng0_432[k]
                   - f_9 * sng1_432[k]
                   + f_3 * pc_x[k] * snh_600[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *smi0_588 = buffer.data(smi0 + 588);
    const auto *smi0_591 = buffer.data(smi0 + 591);
    const auto *smi0_594 = buffer.data(smi0 + 594);
    const auto *smi0_598 = buffer.data(smi0 + 598);
    const auto *smi0_600 = buffer.data(smi0 + 600);
    const auto *smi0_609 = buffer.data(smi0 + 609);

    const auto *smh_441 = buffer.data(smh + 441);
    const auto *smh_444 = buffer.data(smh + 444);
    const auto *smh_447 = buffer.data(smh + 447);
    const auto *smh_448 = buffer.data(smh + 448);
    const auto *smh_450 = buffer.data(smh + 450);
    const auto *smh_456 = buffer.data(smh + 456);
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
    const auto *smh_477 = buffer.data(smh + 477);
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
    const auto *smh_498 = buffer.data(smh + 498);
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
    const auto *smh_519 = buffer.data(smh + 519);
    const auto *smh_521 = buffer.data(smh + 521);
    const auto *smh_522 = buffer.data(smh + 522);
    const auto *smh_523 = buffer.data(smh + 523);
    const auto *smh_524 = buffer.data(smh + 524);
    const auto *smh_525 = buffer.data(smh + 525);
    const auto *smh_527 = buffer.data(smh + 527);
    const auto *smh_530 = buffer.data(smh + 530);
    const auto *smh_602 = buffer.data(smh + 602);
    const auto *smh_603 = buffer.data(smh + 603);
    const auto *smh_604 = buffer.data(smh + 604);
    const auto *smh_605 = buffer.data(smh + 605);
    const auto *smh_606 = buffer.data(smh + 606);
    const auto *smh_607 = buffer.data(smh + 607);
    const auto *smh_608 = buffer.data(smh + 608);
    const auto *smh_614 = buffer.data(smh + 614);
    const auto *smh_618 = buffer.data(smh + 618);
    const auto *smh_623 = buffer.data(smh + 623);
    const auto *smh_624 = buffer.data(smh + 624);
    const auto *smh_625 = buffer.data(smh + 625);
    const auto *smh_626 = buffer.data(smh + 626);
    const auto *smh_627 = buffer.data(smh + 627);
    const auto *smh_628 = buffer.data(smh + 628);
    const auto *smh_629 = buffer.data(smh + 629);
    const auto *smh_630 = buffer.data(smh + 630);
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
    const auto *smh_675 = buffer.data(smh + 675);
    const auto *smh_677 = buffer.data(smh + 677);
    const auto *smh_678 = buffer.data(smh + 678);
    const auto *smh_681 = buffer.data(smh + 681);
    const auto *smh_682 = buffer.data(smh + 682);

    const auto *smi1_588 = buffer.data(smi1 + 588);
    const auto *smi1_591 = buffer.data(smi1 + 591);
    const auto *smi1_594 = buffer.data(smi1 + 594);
    const auto *smi1_598 = buffer.data(smi1 + 598);
    const auto *smi1_600 = buffer.data(smi1 + 600);
    const auto *smi1_609 = buffer.data(smi1 + 609);

    const auto *sng0_430 = buffer.data(sng0 + 430);
    const auto *sng0_432 = buffer.data(sng0 + 432);
    const auto *sng0_433 = buffer.data(sng0 + 433);
    const auto *sng0_434 = buffer.data(sng0 + 434);
    const auto *sng0_440 = buffer.data(sng0 + 440);
    const auto *sng0_444 = buffer.data(sng0 + 444);
    const auto *sng0_447 = buffer.data(sng0 + 447);
    const auto *sng0_448 = buffer.data(sng0 + 448);
    const auto *sng0_449 = buffer.data(sng0 + 449);
    const auto *sng0_450 = buffer.data(sng0 + 450);
    const auto *sng0_453 = buffer.data(sng0 + 453);
    const auto *sng0_455 = buffer.data(sng0 + 455);
    const auto *sng0_456 = buffer.data(sng0 + 456);
    const auto *sng0_459 = buffer.data(sng0 + 459);
    const auto *sng0_460 = buffer.data(sng0 + 460);
    const auto *sng0_462 = buffer.data(sng0 + 462);
    const auto *sng0_463 = buffer.data(sng0 + 463);
    const auto *sng0_464 = buffer.data(sng0 + 464);
    const auto *sng0_465 = buffer.data(sng0 + 465);
    const auto *sng0_468 = buffer.data(sng0 + 468);
    const auto *sng0_470 = buffer.data(sng0 + 470);
    const auto *sng0_471 = buffer.data(sng0 + 471);
    const auto *sng0_474 = buffer.data(sng0 + 474);
    const auto *sng0_475 = buffer.data(sng0 + 475);
    const auto *sng0_477 = buffer.data(sng0 + 477);
    const auto *sng0_478 = buffer.data(sng0 + 478);
    const auto *sng0_479 = buffer.data(sng0 + 479);
    const auto *sng0_480 = buffer.data(sng0 + 480);
    const auto *sng0_483 = buffer.data(sng0 + 483);
    const auto *sng0_485 = buffer.data(sng0 + 485);
    const auto *sng0_486 = buffer.data(sng0 + 486);
    const auto *sng0_489 = buffer.data(sng0 + 489);
    const auto *sng0_490 = buffer.data(sng0 + 490);

    const auto *sng1_430 = buffer.data(sng1 + 430);
    const auto *sng1_432 = buffer.data(sng1 + 432);
    const auto *sng1_433 = buffer.data(sng1 + 433);
    const auto *sng1_434 = buffer.data(sng1 + 434);
    const auto *sng1_440 = buffer.data(sng1 + 440);
    const auto *sng1_444 = buffer.data(sng1 + 444);
    const auto *sng1_447 = buffer.data(sng1 + 447);
    const auto *sng1_448 = buffer.data(sng1 + 448);
    const auto *sng1_449 = buffer.data(sng1 + 449);
    const auto *sng1_450 = buffer.data(sng1 + 450);
    const auto *sng1_453 = buffer.data(sng1 + 453);
    const auto *sng1_455 = buffer.data(sng1 + 455);
    const auto *sng1_456 = buffer.data(sng1 + 456);
    const auto *sng1_459 = buffer.data(sng1 + 459);
    const auto *sng1_460 = buffer.data(sng1 + 460);
    const auto *sng1_462 = buffer.data(sng1 + 462);
    const auto *sng1_463 = buffer.data(sng1 + 463);
    const auto *sng1_464 = buffer.data(sng1 + 464);
    const auto *sng1_465 = buffer.data(sng1 + 465);
    const auto *sng1_468 = buffer.data(sng1 + 468);
    const auto *sng1_470 = buffer.data(sng1 + 470);
    const auto *sng1_471 = buffer.data(sng1 + 471);
    const auto *sng1_474 = buffer.data(sng1 + 474);
    const auto *sng1_475 = buffer.data(sng1 + 475);
    const auto *sng1_477 = buffer.data(sng1 + 477);
    const auto *sng1_478 = buffer.data(sng1 + 478);
    const auto *sng1_479 = buffer.data(sng1 + 479);
    const auto *sng1_480 = buffer.data(sng1 + 480);
    const auto *sng1_483 = buffer.data(sng1 + 483);
    const auto *sng1_485 = buffer.data(sng1 + 485);
    const auto *sng1_486 = buffer.data(sng1 + 486);
    const auto *sng1_489 = buffer.data(sng1 + 489);
    const auto *sng1_490 = buffer.data(sng1 + 490);

    const auto *snh_597 = buffer.data(snh + 597);
    const auto *snh_602 = buffer.data(snh + 602);
    const auto *snh_603 = buffer.data(snh + 603);
    const auto *snh_604 = buffer.data(snh + 604);
    const auto *snh_605 = buffer.data(snh + 605);
    const auto *snh_606 = buffer.data(snh + 606);
    const auto *snh_607 = buffer.data(snh + 607);
    const auto *snh_608 = buffer.data(snh + 608);
    const auto *snh_609 = buffer.data(snh + 609);
    const auto *snh_611 = buffer.data(snh + 611);
    const auto *snh_612 = buffer.data(snh + 612);
    const auto *snh_614 = buffer.data(snh + 614);
    const auto *snh_615 = buffer.data(snh + 615);
    const auto *snh_618 = buffer.data(snh + 618);
    const auto *snh_623 = buffer.data(snh + 623);
    const auto *snh_624 = buffer.data(snh + 624);
    const auto *snh_625 = buffer.data(snh + 625);
    const auto *snh_626 = buffer.data(snh + 626);
    const auto *snh_627 = buffer.data(snh + 627);
    const auto *snh_628 = buffer.data(snh + 628);
    const auto *snh_629 = buffer.data(snh + 629);
    const auto *snh_630 = buffer.data(snh + 630);
    const auto *snh_632 = buffer.data(snh + 632);
    const auto *snh_633 = buffer.data(snh + 633);
    const auto *snh_635 = buffer.data(snh + 635);
    const auto *snh_636 = buffer.data(snh + 636);
    const auto *snh_639 = buffer.data(snh + 639);
    const auto *snh_640 = buffer.data(snh + 640);
    const auto *snh_642 = buffer.data(snh + 642);
    const auto *snh_644 = buffer.data(snh + 644);
    const auto *snh_645 = buffer.data(snh + 645);
    const auto *snh_646 = buffer.data(snh + 646);
    const auto *snh_647 = buffer.data(snh + 647);
    const auto *snh_648 = buffer.data(snh + 648);
    const auto *snh_649 = buffer.data(snh + 649);
    const auto *snh_650 = buffer.data(snh + 650);
    const auto *snh_651 = buffer.data(snh + 651);
    const auto *snh_653 = buffer.data(snh + 653);
    const auto *snh_654 = buffer.data(snh + 654);
    const auto *snh_656 = buffer.data(snh + 656);
    const auto *snh_657 = buffer.data(snh + 657);
    const auto *snh_660 = buffer.data(snh + 660);
    const auto *snh_661 = buffer.data(snh + 661);
    const auto *snh_663 = buffer.data(snh + 663);
    const auto *snh_665 = buffer.data(snh + 665);
    const auto *snh_666 = buffer.data(snh + 666);
    const auto *snh_667 = buffer.data(snh + 667);
    const auto *snh_668 = buffer.data(snh + 668);
    const auto *snh_669 = buffer.data(snh + 669);
    const auto *snh_670 = buffer.data(snh + 670);
    const auto *snh_671 = buffer.data(snh + 671);
    const auto *snh_672 = buffer.data(snh + 672);
    const auto *snh_674 = buffer.data(snh + 674);
    const auto *snh_675 = buffer.data(snh + 675);
    const auto *snh_677 = buffer.data(snh + 677);
    const auto *snh_678 = buffer.data(snh + 678);
    const auto *snh_681 = buffer.data(snh + 681);
    const auto *snh_682 = buffer.data(snh + 682);

#pragma omp simd aligned(t_797, t_798, t_799, t_800, pc_x, pc_y, smh_450, smh_602, smh_603, \
                         smh_604, sng0_434, sng1_434, snh_597, snh_602, snh_603, \
                         snh_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_17 * smh_450[k]
                   + f_3 * pc_y[k] * snh_597[k];

        t_798[k] = f_13 * smh_602[k]
                   + f_8 * sng0_434[k]
                   - f_9 * sng1_434[k]
                   + f_3 * pc_x[k] * snh_602[k];

        t_799[k] = f_13 * smh_603[k]
                   + f_3 * pc_x[k] * snh_603[k];

        t_800[k] = f_13 * smh_604[k]
                   + f_3 * pc_x[k] * snh_604[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, t_804, pc_x, smh_605, smh_606, smh_607, smh_608, \
                         snh_605, snh_606, snh_607, snh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_13 * smh_605[k]
                   + f_3 * pc_x[k] * snh_605[k];

        t_802[k] = f_13 * smh_606[k]
                   + f_3 * pc_x[k] * snh_606[k];

        t_803[k] = f_13 * smh_607[k]
                   + f_3 * pc_x[k] * snh_607[k];

        t_804[k] = f_13 * smh_608[k]
                   + f_3 * pc_x[k] * snh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pc_y, pc_z, smh_456, smh_458, sng0_430, \
                         sng0_432, sng1_430, sng1_432, snh_603, \
                         snh_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_17 * smh_456[k]
                   + f_1 * sng0_430[k]
                   - f_2 * sng1_430[k]
                   + f_3 * pc_y[k] * snh_603[k];

        t_806[k] = f_3 * pc_z[k] * snh_603[k];

        t_807[k] = f_17 * smh_458[k]
                   + f_4 * sng0_432[k]
                   - f_5 * sng1_432[k]
                   + f_3 * pc_y[k] * snh_605[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pc_y, pc_z, smh_459, smh_460, smh_461, \
                         sng0_433, sng0_434, sng1_433, sng1_434, snh_606, snh_607, \
                         snh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_17 * smh_459[k]
                   + f_6 * sng0_433[k]
                   - f_7 * sng1_433[k]
                   + f_3 * pc_y[k] * snh_606[k];

        t_809[k] = f_17 * smh_460[k]
                   + f_8 * sng0_434[k]
                   - f_9 * sng1_434[k]
                   + f_3 * pc_y[k] * snh_607[k];

        t_810[k] = f_17 * smh_461[k]
                   + f_3 * pc_y[k] * snh_608[k];

        t_811[k] = f_1 * sng0_434[k]
                   - f_2 * sng1_434[k]
                   + f_3 * pc_z[k] * snh_608[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pb_z, pc_y, pc_z, smi0_588, smi0_591, \
                         smh_441, smh_462, smi1_588, smi1_591, \
                         snh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = pb_z[k] * smi0_588[k]
                   - f_10 * pc_z[k] * smi1_588[k];

        t_813[k] = f_18 * smh_462[k]
                   + f_3 * pc_y[k] * snh_609[k];

        t_814[k] = f_11 * smh_441[k]
                   + f_3 * pc_z[k] * snh_609[k];

        t_815[k] = pb_z[k] * smi0_591[k]
                   - f_10 * pc_z[k] * smi1_591[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pb_z, pc_x, pc_y, pc_z, smi0_594, smh_464, \
                         smh_614, smi1_594, sng0_440, sng1_440, snh_611, \
                         snh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_18 * smh_464[k]
                   + f_3 * pc_y[k] * snh_611[k];

        t_817[k] = f_13 * smh_614[k]
                   + f_4 * sng0_440[k]
                   - f_5 * sng1_440[k]
                   + f_3 * pc_x[k] * snh_614[k];

        t_818[k] = pb_z[k] * smi0_594[k]
                   - f_10 * pc_z[k] * smi1_594[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, pc_y, pc_z, smh_444, smh_467, smh_618, \
                         sng0_444, sng1_444, snh_612, snh_614, \
                         snh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_11 * smh_444[k]
                   + f_3 * pc_z[k] * snh_612[k];

        t_820[k] = f_18 * smh_467[k]
                   + f_3 * pc_y[k] * snh_614[k];

        t_821[k] = f_13 * smh_618[k]
                   + f_6 * sng0_444[k]
                   - f_7 * sng1_444[k]
                   + f_3 * pc_x[k] * snh_618[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pb_z, pc_y, pc_z, smi0_598, smi0_600, \
                         smh_447, smh_448, smh_471, smi1_598, smi1_600, snh_615, \
                         snh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pb_z[k] * smi0_598[k]
                   - f_10 * pc_z[k] * smi1_598[k];

        t_823[k] = f_11 * smh_447[k]
                   + f_3 * pc_z[k] * snh_615[k];

        t_824[k] = pb_z[k] * smi0_600[k]
                   + f_12 * smh_448[k]
                   - f_10 * pc_z[k] * smi1_600[k];

        t_825[k] = f_18 * smh_471[k]
                   + f_3 * pc_y[k] * snh_618[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, smh_623, smh_624, smh_625, smh_626, \
                         sng0_449, sng1_449, snh_623, snh_624, snh_625, \
                         snh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_13 * smh_623[k]
                   + f_8 * sng0_449[k]
                   - f_9 * sng1_449[k]
                   + f_3 * pc_x[k] * snh_623[k];

        t_827[k] = f_13 * smh_624[k]
                   + f_3 * pc_x[k] * snh_624[k];

        t_828[k] = f_13 * smh_625[k]
                   + f_3 * pc_x[k] * snh_625[k];

        t_829[k] = f_13 * smh_626[k]
                   + f_3 * pc_x[k] * snh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pb_z, pc_x, pc_z, smi0_609, smh_627, \
                         smh_628, smh_629, smi1_609, snh_627, snh_628, \
                         snh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_13 * smh_627[k]
                   + f_3 * pc_x[k] * snh_627[k];

        t_831[k] = f_13 * smh_628[k]
                   + f_3 * pc_x[k] * snh_628[k];

        t_832[k] = f_13 * smh_629[k]
                   + f_3 * pc_x[k] * snh_629[k];

        t_833[k] = pb_z[k] * smi0_609[k]
                   - f_10 * pc_z[k] * smi1_609[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_y, pc_z, smh_456, smh_479, smh_480, sng0_447, \
                         sng0_448, sng1_447, sng1_448, snh_624, snh_626, \
                         snh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_11 * smh_456[k]
                   + f_3 * pc_z[k] * snh_624[k];

        t_835[k] = f_18 * smh_479[k]
                   + f_4 * sng0_447[k]
                   - f_5 * sng1_447[k]
                   + f_3 * pc_y[k] * snh_626[k];

        t_836[k] = f_18 * smh_480[k]
                   + f_6 * sng0_448[k]
                   - f_7 * sng1_448[k]
                   + f_3 * pc_y[k] * snh_627[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, smh_461, smh_481, smh_482, sng0_449, \
                         sng1_449, snh_628, snh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_18 * smh_481[k]
                   + f_8 * sng0_449[k]
                   - f_9 * sng1_449[k]
                   + f_3 * pc_y[k] * snh_628[k];

        t_838[k] = f_18 * smh_482[k]
                   + f_3 * pc_y[k] * snh_629[k];

        t_839[k] = f_11 * smh_461[k]
                   + f_1 * sng0_449[k]
                   - f_2 * sng1_449[k]
                   + f_3 * pc_z[k] * snh_629[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, pc_y, pc_z, smh_462, smh_483, smh_630, \
                         sng0_450, sng1_450, snh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_13 * smh_630[k]
                   + f_1 * sng0_450[k]
                   - f_2 * sng1_450[k]
                   + f_3 * pc_x[k] * snh_630[k];

        t_841[k] = f_19 * smh_483[k]
                   + f_3 * pc_y[k] * snh_630[k];

        t_842[k] = f_12 * smh_462[k]
                   + f_3 * pc_z[k] * snh_630[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, pc_y, smh_485, smh_633, smh_635, sng0_453, \
                         sng0_455, sng1_453, sng1_455, snh_632, snh_633, \
                         snh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_13 * smh_633[k]
                   + f_4 * sng0_453[k]
                   - f_5 * sng1_453[k]
                   + f_3 * pc_x[k] * snh_633[k];

        t_844[k] = f_19 * smh_485[k]
                   + f_3 * pc_y[k] * snh_632[k];

        t_845[k] = f_13 * smh_635[k]
                   + f_4 * sng0_455[k]
                   - f_5 * sng1_455[k]
                   + f_3 * pc_x[k] * snh_635[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pc_x, pc_y, pc_z, smh_465, smh_488, smh_636, \
                         sng0_456, sng1_456, snh_633, snh_635, \
                         snh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_13 * smh_636[k]
                   + f_6 * sng0_456[k]
                   - f_7 * sng1_456[k]
                   + f_3 * pc_x[k] * snh_636[k];

        t_847[k] = f_12 * smh_465[k]
                   + f_3 * pc_z[k] * snh_633[k];

        t_848[k] = f_19 * smh_488[k]
                   + f_3 * pc_y[k] * snh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, smh_468, smh_639, smh_640, sng0_459, \
                         sng0_460, sng1_459, sng1_460, snh_636, snh_639, \
                         snh_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_13 * smh_639[k]
                   + f_6 * sng0_459[k]
                   - f_7 * sng1_459[k]
                   + f_3 * pc_x[k] * snh_639[k];

        t_850[k] = f_13 * smh_640[k]
                   + f_8 * sng0_460[k]
                   - f_9 * sng1_460[k]
                   + f_3 * pc_x[k] * snh_640[k];

        t_851[k] = f_12 * smh_468[k]
                   + f_3 * pc_z[k] * snh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, pc_x, pc_y, smh_492, smh_642, smh_644, sng0_462, \
                         sng0_464, sng1_462, sng1_464, snh_639, snh_642, \
                         snh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_13 * smh_642[k]
                   + f_8 * sng0_462[k]
                   - f_9 * sng1_462[k]
                   + f_3 * pc_x[k] * snh_642[k];

        t_853[k] = f_19 * smh_492[k]
                   + f_3 * pc_y[k] * snh_639[k];

        t_854[k] = f_13 * smh_644[k]
                   + f_8 * sng0_464[k]
                   - f_9 * sng1_464[k]
                   + f_3 * pc_x[k] * snh_644[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, pc_x, smh_645, smh_646, smh_647, \
                         smh_648, smh_649, snh_645, snh_646, snh_647, snh_648, \
                         snh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_13 * smh_645[k]
                   + f_3 * pc_x[k] * snh_645[k];

        t_856[k] = f_13 * smh_646[k]
                   + f_3 * pc_x[k] * snh_646[k];

        t_857[k] = f_13 * smh_647[k]
                   + f_3 * pc_x[k] * snh_647[k];

        t_858[k] = f_13 * smh_648[k]
                   + f_3 * pc_x[k] * snh_648[k];

        t_859[k] = f_13 * smh_649[k]
                   + f_3 * pc_x[k] * snh_649[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pc_x, pc_y, pc_z, smh_477, smh_498, smh_650, \
                         sng0_460, sng1_460, snh_645, snh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_13 * smh_650[k]
                   + f_3 * pc_x[k] * snh_650[k];

        t_861[k] = f_19 * smh_498[k]
                   + f_1 * sng0_460[k]
                   - f_2 * sng1_460[k]
                   + f_3 * pc_y[k] * snh_645[k];

        t_862[k] = f_12 * smh_477[k]
                   + f_3 * pc_z[k] * snh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_y, smh_500, smh_501, smh_502, sng0_462, \
                         sng0_463, sng0_464, sng1_462, sng1_463, sng1_464, snh_647, snh_648, \
                         snh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_19 * smh_500[k]
                   + f_4 * sng0_462[k]
                   - f_5 * sng1_462[k]
                   + f_3 * pc_y[k] * snh_647[k];

        t_864[k] = f_19 * smh_501[k]
                   + f_6 * sng0_463[k]
                   - f_7 * sng1_463[k]
                   + f_3 * pc_y[k] * snh_648[k];

        t_865[k] = f_19 * smh_502[k]
                   + f_8 * sng0_464[k]
                   - f_9 * sng1_464[k]
                   + f_3 * pc_y[k] * snh_649[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, smh_482, smh_503, smh_651, \
                         sng0_464, sng0_465, sng1_464, sng1_465, snh_650, \
                         snh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_19 * smh_503[k]
                   + f_3 * pc_y[k] * snh_650[k];

        t_867[k] = f_12 * smh_482[k]
                   + f_1 * sng0_464[k]
                   - f_2 * sng1_464[k]
                   + f_3 * pc_z[k] * snh_650[k];

        t_868[k] = f_13 * smh_651[k]
                   + f_1 * sng0_465[k]
                   - f_2 * sng1_465[k]
                   + f_3 * pc_x[k] * snh_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, pc_z, smh_483, smh_504, \
                         smh_506, smh_654, sng0_468, sng1_468, snh_651, snh_653, \
                         snh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * smh_504[k]
                   + f_3 * pc_y[k] * snh_651[k];

        t_870[k] = f_13 * smh_483[k]
                   + f_3 * pc_z[k] * snh_651[k];

        t_871[k] = f_13 * smh_654[k]
                   + f_4 * sng0_468[k]
                   - f_5 * sng1_468[k]
                   + f_3 * pc_x[k] * snh_654[k];

        t_872[k] = f_14 * smh_506[k]
                   + f_3 * pc_y[k] * snh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_z, smh_486, smh_656, smh_657, sng0_470, \
                         sng0_471, sng1_470, sng1_471, snh_654, snh_656, \
                         snh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_13 * smh_656[k]
                   + f_4 * sng0_470[k]
                   - f_5 * sng1_470[k]
                   + f_3 * pc_x[k] * snh_656[k];

        t_874[k] = f_13 * smh_657[k]
                   + f_6 * sng0_471[k]
                   - f_7 * sng1_471[k]
                   + f_3 * pc_x[k] * snh_657[k];

        t_875[k] = f_13 * smh_486[k]
                   + f_3 * pc_z[k] * snh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_y, smh_509, smh_660, smh_661, sng0_474, \
                         sng0_475, sng1_474, sng1_475, snh_656, snh_660, \
                         snh_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * smh_509[k]
                   + f_3 * pc_y[k] * snh_656[k];

        t_877[k] = f_13 * smh_660[k]
                   + f_6 * sng0_474[k]
                   - f_7 * sng1_474[k]
                   + f_3 * pc_x[k] * snh_660[k];

        t_878[k] = f_13 * smh_661[k]
                   + f_8 * sng0_475[k]
                   - f_9 * sng1_475[k]
                   + f_3 * pc_x[k] * snh_661[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, smh_489, smh_513, smh_663, \
                         sng0_477, sng1_477, snh_657, snh_660, \
                         snh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * smh_489[k]
                   + f_3 * pc_z[k] * snh_657[k];

        t_880[k] = f_13 * smh_663[k]
                   + f_8 * sng0_477[k]
                   - f_9 * sng1_477[k]
                   + f_3 * pc_x[k] * snh_663[k];

        t_881[k] = f_14 * smh_513[k]
                   + f_3 * pc_y[k] * snh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pc_x, smh_665, smh_666, smh_667, smh_668, \
                         sng0_479, sng1_479, snh_665, snh_666, snh_667, \
                         snh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_13 * smh_665[k]
                   + f_8 * sng0_479[k]
                   - f_9 * sng1_479[k]
                   + f_3 * pc_x[k] * snh_665[k];

        t_883[k] = f_13 * smh_666[k]
                   + f_3 * pc_x[k] * snh_666[k];

        t_884[k] = f_13 * smh_667[k]
                   + f_3 * pc_x[k] * snh_667[k];

        t_885[k] = f_13 * smh_668[k]
                   + f_3 * pc_x[k] * snh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pc_x, pc_y, smh_519, smh_669, smh_670, \
                         smh_671, sng0_475, sng1_475, snh_666, snh_669, snh_670, \
                         snh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_13 * smh_669[k]
                   + f_3 * pc_x[k] * snh_669[k];

        t_887[k] = f_13 * smh_670[k]
                   + f_3 * pc_x[k] * snh_670[k];

        t_888[k] = f_13 * smh_671[k]
                   + f_3 * pc_x[k] * snh_671[k];

        t_889[k] = f_14 * smh_519[k]
                   + f_1 * sng0_475[k]
                   - f_2 * sng1_475[k]
                   + f_3 * pc_y[k] * snh_666[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pc_y, pc_z, smh_498, smh_521, smh_522, sng0_477, \
                         sng0_478, sng1_477, sng1_478, snh_666, snh_668, \
                         snh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * smh_498[k]
                   + f_3 * pc_z[k] * snh_666[k];

        t_891[k] = f_14 * smh_521[k]
                   + f_4 * sng0_477[k]
                   - f_5 * sng1_477[k]
                   + f_3 * pc_y[k] * snh_668[k];

        t_892[k] = f_14 * smh_522[k]
                   + f_6 * sng0_478[k]
                   - f_7 * sng1_478[k]
                   + f_3 * pc_y[k] * snh_669[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, pc_z, smh_503, smh_523, smh_524, sng0_479, \
                         sng1_479, snh_670, snh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_14 * smh_523[k]
                   + f_8 * sng0_479[k]
                   - f_9 * sng1_479[k]
                   + f_3 * pc_y[k] * snh_670[k];

        t_894[k] = f_14 * smh_524[k]
                   + f_3 * pc_y[k] * snh_671[k];

        t_895[k] = f_13 * smh_503[k]
                   + f_1 * sng0_479[k]
                   - f_2 * sng1_479[k]
                   + f_3 * pc_z[k] * snh_671[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_x, pc_y, pc_z, smh_504, smh_525, smh_672, \
                         sng0_480, sng1_480, snh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_13 * smh_672[k]
                   + f_1 * sng0_480[k]
                   - f_2 * sng1_480[k]
                   + f_3 * pc_x[k] * snh_672[k];

        t_897[k] = f_13 * smh_525[k]
                   + f_3 * pc_y[k] * snh_672[k];

        t_898[k] = f_14 * smh_504[k]
                   + f_3 * pc_z[k] * snh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, smh_527, smh_675, smh_677, sng0_483, \
                         sng0_485, sng1_483, sng1_485, snh_674, snh_675, \
                         snh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_13 * smh_675[k]
                   + f_4 * sng0_483[k]
                   - f_5 * sng1_483[k]
                   + f_3 * pc_x[k] * snh_675[k];

        t_900[k] = f_13 * smh_527[k]
                   + f_3 * pc_y[k] * snh_674[k];

        t_901[k] = f_13 * smh_677[k]
                   + f_4 * sng0_485[k]
                   - f_5 * sng1_485[k]
                   + f_3 * pc_x[k] * snh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, smh_507, smh_530, smh_678, \
                         sng0_486, sng1_486, snh_675, snh_677, \
                         snh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_13 * smh_678[k]
                   + f_6 * sng0_486[k]
                   - f_7 * sng1_486[k]
                   + f_3 * pc_x[k] * snh_678[k];

        t_903[k] = f_14 * smh_507[k]
                   + f_3 * pc_z[k] * snh_675[k];

        t_904[k] = f_13 * smh_530[k]
                   + f_3 * pc_y[k] * snh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, smh_510, smh_681, smh_682, sng0_489, \
                         sng0_490, sng1_489, sng1_490, snh_678, snh_681, \
                         snh_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_13 * smh_681[k]
                   + f_6 * sng0_489[k]
                   - f_7 * sng1_489[k]
                   + f_3 * pc_x[k] * snh_681[k];

        t_906[k] = f_13 * smh_682[k]
                   + f_8 * sng0_490[k]
                   - f_9 * sng1_490[k]
                   + f_3 * pc_x[k] * snh_682[k];

        t_907[k] = f_14 * smh_510[k]
                   + f_3 * pc_z[k] * snh_678[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi0_756 = buffer.data(smi0 + 756);
    const auto *smi0_759 = buffer.data(smi0 + 759);
    const auto *smi0_761 = buffer.data(smi0 + 761);
    const auto *smi0_762 = buffer.data(smi0 + 762);
    const auto *smi0_765 = buffer.data(smi0 + 765);
    const auto *smi0_766 = buffer.data(smi0 + 766);
    const auto *smi0_768 = buffer.data(smi0 + 768);
    const auto *smi0_770 = buffer.data(smi0 + 770);
    const auto *smi0_783 = buffer.data(smi0 + 783);

    const auto *smh_519 = buffer.data(smh + 519);
    const auto *smh_524 = buffer.data(smh + 524);
    const auto *smh_525 = buffer.data(smh + 525);
    const auto *smh_528 = buffer.data(smh + 528);
    const auto *smh_531 = buffer.data(smh + 531);
    const auto *smh_534 = buffer.data(smh + 534);
    const auto *smh_540 = buffer.data(smh + 540);
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
    const auto *smh_563 = buffer.data(smh + 563);
    const auto *smh_564 = buffer.data(smh + 564);
    const auto *smh_565 = buffer.data(smh + 565);
    const auto *smh_566 = buffer.data(smh + 566);
    const auto *smh_567 = buffer.data(smh + 567);
    const auto *smh_568 = buffer.data(smh + 568);
    const auto *smh_569 = buffer.data(smh + 569);
    const auto *smh_570 = buffer.data(smh + 570);
    const auto *smh_572 = buffer.data(smh + 572);
    const auto *smh_573 = buffer.data(smh + 573);
    const auto *smh_575 = buffer.data(smh + 575);
    const auto *smh_576 = buffer.data(smh + 576);
    const auto *smh_582 = buffer.data(smh + 582);
    const auto *smh_584 = buffer.data(smh + 584);
    const auto *smh_585 = buffer.data(smh + 585);
    const auto *smh_586 = buffer.data(smh + 586);
    const auto *smh_587 = buffer.data(smh + 587);
    const auto *smh_588 = buffer.data(smh + 588);
    const auto *smh_590 = buffer.data(smh + 590);
    const auto *smh_593 = buffer.data(smh + 593);
    const auto *smh_684 = buffer.data(smh + 684);
    const auto *smh_686 = buffer.data(smh + 686);
    const auto *smh_687 = buffer.data(smh + 687);
    const auto *smh_688 = buffer.data(smh + 688);
    const auto *smh_689 = buffer.data(smh + 689);
    const auto *smh_690 = buffer.data(smh + 690);
    const auto *smh_691 = buffer.data(smh + 691);
    const auto *smh_692 = buffer.data(smh + 692);
    const auto *smh_693 = buffer.data(smh + 693);
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
    const auto *smh_729 = buffer.data(smh + 729);
    const auto *smh_730 = buffer.data(smh + 730);
    const auto *smh_731 = buffer.data(smh + 731);
    const auto *smh_732 = buffer.data(smh + 732);
    const auto *smh_733 = buffer.data(smh + 733);
    const auto *smh_734 = buffer.data(smh + 734);
    const auto *smh_735 = buffer.data(smh + 735);
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
    const auto *smh_759 = buffer.data(smh + 759);
    const auto *smh_761 = buffer.data(smh + 761);
    const auto *smh_762 = buffer.data(smh + 762);
    const auto *smh_765 = buffer.data(smh + 765);
    const auto *smh_766 = buffer.data(smh + 766);
    const auto *smh_768 = buffer.data(smh + 768);

    const auto *smi1_756 = buffer.data(smi1 + 756);
    const auto *smi1_759 = buffer.data(smi1 + 759);
    const auto *smi1_761 = buffer.data(smi1 + 761);
    const auto *smi1_762 = buffer.data(smi1 + 762);
    const auto *smi1_765 = buffer.data(smi1 + 765);
    const auto *smi1_766 = buffer.data(smi1 + 766);
    const auto *smi1_768 = buffer.data(smi1 + 768);
    const auto *smi1_770 = buffer.data(smi1 + 770);
    const auto *smi1_783 = buffer.data(smi1 + 783);

    const auto *sng0_490 = buffer.data(sng0 + 490);
    const auto *sng0_492 = buffer.data(sng0 + 492);
    const auto *sng0_493 = buffer.data(sng0 + 493);
    const auto *sng0_494 = buffer.data(sng0 + 494);
    const auto *sng0_495 = buffer.data(sng0 + 495);
    const auto *sng0_498 = buffer.data(sng0 + 498);
    const auto *sng0_500 = buffer.data(sng0 + 500);
    const auto *sng0_501 = buffer.data(sng0 + 501);
    const auto *sng0_504 = buffer.data(sng0 + 504);
    const auto *sng0_505 = buffer.data(sng0 + 505);
    const auto *sng0_507 = buffer.data(sng0 + 507);
    const auto *sng0_508 = buffer.data(sng0 + 508);
    const auto *sng0_509 = buffer.data(sng0 + 509);
    const auto *sng0_520 = buffer.data(sng0 + 520);
    const auto *sng0_522 = buffer.data(sng0 + 522);
    const auto *sng0_523 = buffer.data(sng0 + 523);
    const auto *sng0_524 = buffer.data(sng0 + 524);
    const auto *sng0_525 = buffer.data(sng0 + 525);
    const auto *sng0_528 = buffer.data(sng0 + 528);
    const auto *sng0_530 = buffer.data(sng0 + 530);
    const auto *sng0_531 = buffer.data(sng0 + 531);
    const auto *sng0_534 = buffer.data(sng0 + 534);
    const auto *sng0_535 = buffer.data(sng0 + 535);
    const auto *sng0_537 = buffer.data(sng0 + 537);
    const auto *sng0_538 = buffer.data(sng0 + 538);
    const auto *sng0_539 = buffer.data(sng0 + 539);
    const auto *sng0_540 = buffer.data(sng0 + 540);
    const auto *sng0_543 = buffer.data(sng0 + 543);
    const auto *sng0_545 = buffer.data(sng0 + 545);
    const auto *sng0_546 = buffer.data(sng0 + 546);
    const auto *sng0_549 = buffer.data(sng0 + 549);
    const auto *sng0_550 = buffer.data(sng0 + 550);
    const auto *sng0_552 = buffer.data(sng0 + 552);

    const auto *sng1_490 = buffer.data(sng1 + 490);
    const auto *sng1_492 = buffer.data(sng1 + 492);
    const auto *sng1_493 = buffer.data(sng1 + 493);
    const auto *sng1_494 = buffer.data(sng1 + 494);
    const auto *sng1_495 = buffer.data(sng1 + 495);
    const auto *sng1_498 = buffer.data(sng1 + 498);
    const auto *sng1_500 = buffer.data(sng1 + 500);
    const auto *sng1_501 = buffer.data(sng1 + 501);
    const auto *sng1_504 = buffer.data(sng1 + 504);
    const auto *sng1_505 = buffer.data(sng1 + 505);
    const auto *sng1_507 = buffer.data(sng1 + 507);
    const auto *sng1_508 = buffer.data(sng1 + 508);
    const auto *sng1_509 = buffer.data(sng1 + 509);
    const auto *sng1_520 = buffer.data(sng1 + 520);
    const auto *sng1_522 = buffer.data(sng1 + 522);
    const auto *sng1_523 = buffer.data(sng1 + 523);
    const auto *sng1_524 = buffer.data(sng1 + 524);
    const auto *sng1_525 = buffer.data(sng1 + 525);
    const auto *sng1_528 = buffer.data(sng1 + 528);
    const auto *sng1_530 = buffer.data(sng1 + 530);
    const auto *sng1_531 = buffer.data(sng1 + 531);
    const auto *sng1_534 = buffer.data(sng1 + 534);
    const auto *sng1_535 = buffer.data(sng1 + 535);
    const auto *sng1_537 = buffer.data(sng1 + 537);
    const auto *sng1_538 = buffer.data(sng1 + 538);
    const auto *sng1_539 = buffer.data(sng1 + 539);
    const auto *sng1_540 = buffer.data(sng1 + 540);
    const auto *sng1_543 = buffer.data(sng1 + 543);
    const auto *sng1_545 = buffer.data(sng1 + 545);
    const auto *sng1_546 = buffer.data(sng1 + 546);
    const auto *sng1_549 = buffer.data(sng1 + 549);
    const auto *sng1_550 = buffer.data(sng1 + 550);
    const auto *sng1_552 = buffer.data(sng1 + 552);

    const auto *snh_681 = buffer.data(snh + 681);
    const auto *snh_684 = buffer.data(snh + 684);
    const auto *snh_686 = buffer.data(snh + 686);
    const auto *snh_687 = buffer.data(snh + 687);
    const auto *snh_688 = buffer.data(snh + 688);
    const auto *snh_689 = buffer.data(snh + 689);
    const auto *snh_690 = buffer.data(snh + 690);
    const auto *snh_691 = buffer.data(snh + 691);
    const auto *snh_692 = buffer.data(snh + 692);
    const auto *snh_693 = buffer.data(snh + 693);
    const auto *snh_695 = buffer.data(snh + 695);
    const auto *snh_696 = buffer.data(snh + 696);
    const auto *snh_698 = buffer.data(snh + 698);
    const auto *snh_699 = buffer.data(snh + 699);
    const auto *snh_702 = buffer.data(snh + 702);
    const auto *snh_703 = buffer.data(snh + 703);
    const auto *snh_705 = buffer.data(snh + 705);
    const auto *snh_707 = buffer.data(snh + 707);
    const auto *snh_708 = buffer.data(snh + 708);
    const auto *snh_709 = buffer.data(snh + 709);
    const auto *snh_710 = buffer.data(snh + 710);
    const auto *snh_711 = buffer.data(snh + 711);
    const auto *snh_712 = buffer.data(snh + 712);
    const auto *snh_713 = buffer.data(snh + 713);
    const auto *snh_714 = buffer.data(snh + 714);
    const auto *snh_716 = buffer.data(snh + 716);
    const auto *snh_717 = buffer.data(snh + 717);
    const auto *snh_719 = buffer.data(snh + 719);
    const auto *snh_720 = buffer.data(snh + 720);
    const auto *snh_723 = buffer.data(snh + 723);
    const auto *snh_729 = buffer.data(snh + 729);
    const auto *snh_730 = buffer.data(snh + 730);
    const auto *snh_731 = buffer.data(snh + 731);
    const auto *snh_732 = buffer.data(snh + 732);
    const auto *snh_733 = buffer.data(snh + 733);
    const auto *snh_734 = buffer.data(snh + 734);
    const auto *snh_735 = buffer.data(snh + 735);
    const auto *snh_737 = buffer.data(snh + 737);
    const auto *snh_738 = buffer.data(snh + 738);
    const auto *snh_740 = buffer.data(snh + 740);
    const auto *snh_741 = buffer.data(snh + 741);
    const auto *snh_744 = buffer.data(snh + 744);
    const auto *snh_745 = buffer.data(snh + 745);
    const auto *snh_747 = buffer.data(snh + 747);
    const auto *snh_749 = buffer.data(snh + 749);
    const auto *snh_750 = buffer.data(snh + 750);
    const auto *snh_751 = buffer.data(snh + 751);
    const auto *snh_752 = buffer.data(snh + 752);
    const auto *snh_753 = buffer.data(snh + 753);
    const auto *snh_754 = buffer.data(snh + 754);
    const auto *snh_755 = buffer.data(snh + 755);
    const auto *snh_756 = buffer.data(snh + 756);
    const auto *snh_758 = buffer.data(snh + 758);
    const auto *snh_759 = buffer.data(snh + 759);
    const auto *snh_761 = buffer.data(snh + 761);
    const auto *snh_762 = buffer.data(snh + 762);
    const auto *snh_765 = buffer.data(snh + 765);
    const auto *snh_766 = buffer.data(snh + 766);
    const auto *snh_768 = buffer.data(snh + 768);

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, smh_534, smh_684, smh_686, sng0_492, \
                         sng0_494, sng1_492, sng1_494, snh_681, snh_684, \
                         snh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_13 * smh_684[k]
                   + f_8 * sng0_492[k]
                   - f_9 * sng1_492[k]
                   + f_3 * pc_x[k] * snh_684[k];

        t_909[k] = f_13 * smh_534[k]
                   + f_3 * pc_y[k] * snh_681[k];

        t_910[k] = f_13 * smh_686[k]
                   + f_8 * sng0_494[k]
                   - f_9 * sng1_494[k]
                   + f_3 * pc_x[k] * snh_686[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pc_x, smh_687, smh_688, smh_689, \
                         smh_690, smh_691, snh_687, snh_688, snh_689, snh_690, \
                         snh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_13 * smh_687[k]
                   + f_3 * pc_x[k] * snh_687[k];

        t_912[k] = f_13 * smh_688[k]
                   + f_3 * pc_x[k] * snh_688[k];

        t_913[k] = f_13 * smh_689[k]
                   + f_3 * pc_x[k] * snh_689[k];

        t_914[k] = f_13 * smh_690[k]
                   + f_3 * pc_x[k] * snh_690[k];

        t_915[k] = f_13 * smh_691[k]
                   + f_3 * pc_x[k] * snh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pc_x, pc_y, pc_z, smh_519, smh_540, smh_692, \
                         sng0_490, sng1_490, snh_687, snh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_13 * smh_692[k]
                   + f_3 * pc_x[k] * snh_692[k];

        t_917[k] = f_13 * smh_540[k]
                   + f_1 * sng0_490[k]
                   - f_2 * sng1_490[k]
                   + f_3 * pc_y[k] * snh_687[k];

        t_918[k] = f_14 * smh_519[k]
                   + f_3 * pc_z[k] * snh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pc_y, smh_542, smh_543, smh_544, sng0_492, \
                         sng0_493, sng0_494, sng1_492, sng1_493, sng1_494, snh_689, snh_690, \
                         snh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_13 * smh_542[k]
                   + f_4 * sng0_492[k]
                   - f_5 * sng1_492[k]
                   + f_3 * pc_y[k] * snh_689[k];

        t_920[k] = f_13 * smh_543[k]
                   + f_6 * sng0_493[k]
                   - f_7 * sng1_493[k]
                   + f_3 * pc_y[k] * snh_690[k];

        t_921[k] = f_13 * smh_544[k]
                   + f_8 * sng0_494[k]
                   - f_9 * sng1_494[k]
                   + f_3 * pc_y[k] * snh_691[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pc_x, pc_y, pc_z, smh_524, smh_545, smh_693, \
                         sng0_494, sng0_495, sng1_494, sng1_495, snh_692, \
                         snh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * smh_545[k]
                   + f_3 * pc_y[k] * snh_692[k];

        t_923[k] = f_14 * smh_524[k]
                   + f_1 * sng0_494[k]
                   - f_2 * sng1_494[k]
                   + f_3 * pc_z[k] * snh_692[k];

        t_924[k] = f_13 * smh_693[k]
                   + f_1 * sng0_495[k]
                   - f_2 * sng1_495[k]
                   + f_3 * pc_x[k] * snh_693[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, smh_525, smh_546, \
                         smh_548, smh_696, sng0_498, sng1_498, snh_693, snh_695, \
                         snh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_12 * smh_546[k]
                   + f_3 * pc_y[k] * snh_693[k];

        t_926[k] = f_19 * smh_525[k]
                   + f_3 * pc_z[k] * snh_693[k];

        t_927[k] = f_13 * smh_696[k]
                   + f_4 * sng0_498[k]
                   - f_5 * sng1_498[k]
                   + f_3 * pc_x[k] * snh_696[k];

        t_928[k] = f_12 * smh_548[k]
                   + f_3 * pc_y[k] * snh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pc_x, pc_z, smh_528, smh_698, smh_699, sng0_500, \
                         sng0_501, sng1_500, sng1_501, snh_696, snh_698, \
                         snh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_13 * smh_698[k]
                   + f_4 * sng0_500[k]
                   - f_5 * sng1_500[k]
                   + f_3 * pc_x[k] * snh_698[k];

        t_930[k] = f_13 * smh_699[k]
                   + f_6 * sng0_501[k]
                   - f_7 * sng1_501[k]
                   + f_3 * pc_x[k] * snh_699[k];

        t_931[k] = f_19 * smh_528[k]
                   + f_3 * pc_z[k] * snh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pc_x, pc_y, smh_551, smh_702, smh_703, sng0_504, \
                         sng0_505, sng1_504, sng1_505, snh_698, snh_702, \
                         snh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * smh_551[k]
                   + f_3 * pc_y[k] * snh_698[k];

        t_933[k] = f_13 * smh_702[k]
                   + f_6 * sng0_504[k]
                   - f_7 * sng1_504[k]
                   + f_3 * pc_x[k] * snh_702[k];

        t_934[k] = f_13 * smh_703[k]
                   + f_8 * sng0_505[k]
                   - f_9 * sng1_505[k]
                   + f_3 * pc_x[k] * snh_703[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pc_x, pc_y, pc_z, smh_531, smh_555, smh_705, \
                         sng0_507, sng1_507, snh_699, snh_702, \
                         snh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_19 * smh_531[k]
                   + f_3 * pc_z[k] * snh_699[k];

        t_936[k] = f_13 * smh_705[k]
                   + f_8 * sng0_507[k]
                   - f_9 * sng1_507[k]
                   + f_3 * pc_x[k] * snh_705[k];

        t_937[k] = f_12 * smh_555[k]
                   + f_3 * pc_y[k] * snh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pc_x, smh_707, smh_708, smh_709, smh_710, \
                         sng0_509, sng1_509, snh_707, snh_708, snh_709, \
                         snh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_13 * smh_707[k]
                   + f_8 * sng0_509[k]
                   - f_9 * sng1_509[k]
                   + f_3 * pc_x[k] * snh_707[k];

        t_939[k] = f_13 * smh_708[k]
                   + f_3 * pc_x[k] * snh_708[k];

        t_940[k] = f_13 * smh_709[k]
                   + f_3 * pc_x[k] * snh_709[k];

        t_941[k] = f_13 * smh_710[k]
                   + f_3 * pc_x[k] * snh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, smh_561, smh_711, smh_712, \
                         smh_713, sng0_505, sng1_505, snh_708, snh_711, snh_712, \
                         snh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_13 * smh_711[k]
                   + f_3 * pc_x[k] * snh_711[k];

        t_943[k] = f_13 * smh_712[k]
                   + f_3 * pc_x[k] * snh_712[k];

        t_944[k] = f_13 * smh_713[k]
                   + f_3 * pc_x[k] * snh_713[k];

        t_945[k] = f_12 * smh_561[k]
                   + f_1 * sng0_505[k]
                   - f_2 * sng1_505[k]
                   + f_3 * pc_y[k] * snh_708[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pc_y, pc_z, smh_540, smh_563, smh_564, sng0_507, \
                         sng0_508, sng1_507, sng1_508, snh_708, snh_710, \
                         snh_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_19 * smh_540[k]
                   + f_3 * pc_z[k] * snh_708[k];

        t_947[k] = f_12 * smh_563[k]
                   + f_4 * sng0_507[k]
                   - f_5 * sng1_507[k]
                   + f_3 * pc_y[k] * snh_710[k];

        t_948[k] = f_12 * smh_564[k]
                   + f_6 * sng0_508[k]
                   - f_7 * sng1_508[k]
                   + f_3 * pc_y[k] * snh_711[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pb_y, pc_y, pc_z, smi0_756, smh_545, \
                         smh_565, smh_566, smi1_756, sng0_509, sng1_509, snh_712, \
                         snh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_12 * smh_565[k]
                   + f_8 * sng0_509[k]
                   - f_9 * sng1_509[k]
                   + f_3 * pc_y[k] * snh_712[k];

        t_950[k] = f_12 * smh_566[k]
                   + f_3 * pc_y[k] * snh_713[k];

        t_951[k] = f_19 * smh_545[k]
                   + f_1 * sng0_509[k]
                   - f_2 * sng1_509[k]
                   + f_3 * pc_z[k] * snh_713[k];

        t_952[k] = pb_y[k] * smi0_756[k]
                   - f_10 * pc_y[k] * smi1_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pb_y, pc_y, pc_z, smi0_759, smh_546, \
                         smh_567, smh_568, smh_569, smi1_759, snh_714, \
                         snh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_11 * smh_567[k]
                   + f_3 * pc_y[k] * snh_714[k];

        t_954[k] = f_18 * smh_546[k]
                   + f_3 * pc_z[k] * snh_714[k];

        t_955[k] = pb_y[k] * smi0_759[k]
                   + f_12 * smh_568[k]
                   - f_10 * pc_y[k] * smi1_759[k];

        t_956[k] = f_11 * smh_569[k]
                   + f_3 * pc_y[k] * snh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pb_y, pc_y, pc_z, smi0_761, smi0_762, \
                         smh_549, smh_570, smh_572, smi1_761, smi1_762, snh_717, \
                         snh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pb_y[k] * smi0_761[k]
                   - f_10 * pc_y[k] * smi1_761[k];

        t_958[k] = pb_y[k] * smi0_762[k]
                   + f_13 * smh_570[k]
                   - f_10 * pc_y[k] * smi1_762[k];

        t_959[k] = f_18 * smh_549[k]
                   + f_3 * pc_z[k] * snh_717[k];

        t_960[k] = f_11 * smh_572[k]
                   + f_3 * pc_y[k] * snh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pb_y, pc_y, pc_z, smi0_765, smi0_766, smh_552, \
                         smh_573, smi1_765, smi1_766, snh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pb_y[k] * smi0_765[k]
                   - f_10 * pc_y[k] * smi1_765[k];

        t_962[k] = pb_y[k] * smi0_766[k]
                   + f_14 * smh_573[k]
                   - f_10 * pc_y[k] * smi1_766[k];

        t_963[k] = f_18 * smh_552[k]
                   + f_3 * pc_z[k] * snh_720[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pb_y, pc_x, pc_y, smi0_768, smi0_770, \
                         smh_575, smh_576, smh_729, smi1_768, smi1_770, snh_723, \
                         snh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pb_y[k] * smi0_768[k]
                   + f_12 * smh_575[k]
                   - f_10 * pc_y[k] * smi1_768[k];

        t_965[k] = f_11 * smh_576[k]
                   + f_3 * pc_y[k] * snh_723[k];

        t_966[k] = pb_y[k] * smi0_770[k]
                   - f_10 * pc_y[k] * smi1_770[k];

        t_967[k] = f_13 * smh_729[k]
                   + f_3 * pc_x[k] * snh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, smh_730, smh_731, smh_732, \
                         smh_733, smh_734, snh_730, snh_731, snh_732, snh_733, \
                         snh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_13 * smh_730[k]
                   + f_3 * pc_x[k] * snh_730[k];

        t_969[k] = f_13 * smh_731[k]
                   + f_3 * pc_x[k] * snh_731[k];

        t_970[k] = f_13 * smh_732[k]
                   + f_3 * pc_x[k] * snh_732[k];

        t_971[k] = f_13 * smh_733[k]
                   + f_3 * pc_x[k] * snh_733[k];

        t_972[k] = f_13 * smh_734[k]
                   + f_3 * pc_x[k] * snh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pc_y, pc_z, smh_561, smh_582, smh_584, sng0_520, \
                         sng0_522, sng1_520, sng1_522, snh_729, \
                         snh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_11 * smh_582[k]
                   + f_1 * sng0_520[k]
                   - f_2 * sng1_520[k]
                   + f_3 * pc_y[k] * snh_729[k];

        t_974[k] = f_18 * smh_561[k]
                   + f_3 * pc_z[k] * snh_729[k];

        t_975[k] = f_11 * smh_584[k]
                   + f_4 * sng0_522[k]
                   - f_5 * sng1_522[k]
                   + f_3 * pc_y[k] * snh_731[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_y, smh_585, smh_586, smh_587, sng0_523, \
                         sng0_524, sng1_523, sng1_524, snh_732, snh_733, \
                         snh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_11 * smh_585[k]
                   + f_6 * sng0_523[k]
                   - f_7 * sng1_523[k]
                   + f_3 * pc_y[k] * snh_732[k];

        t_977[k] = f_11 * smh_586[k]
                   + f_8 * sng0_524[k]
                   - f_9 * sng1_524[k]
                   + f_3 * pc_y[k] * snh_733[k];

        t_978[k] = f_11 * smh_587[k]
                   + f_3 * pc_y[k] * snh_734[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pb_y, pc_x, pc_y, pc_z, smi0_783, \
                         smh_567, smh_735, smi1_783, sng0_525, sng1_525, \
                         snh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pb_y[k] * smi0_783[k]
                   - f_10 * pc_y[k] * smi1_783[k];

        t_980[k] = f_13 * smh_735[k]
                   + f_1 * sng0_525[k]
                   - f_2 * sng1_525[k]
                   + f_3 * pc_x[k] * snh_735[k];

        t_981[k] = f_3 * pc_y[k] * snh_735[k];

        t_982[k] = f_17 * smh_567[k]
                   + f_3 * pc_z[k] * snh_735[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_x, pc_y, smh_738, smh_740, sng0_528, \
                         sng0_530, sng1_528, sng1_530, snh_737, snh_738, \
                         snh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_13 * smh_738[k]
                   + f_4 * sng0_528[k]
                   - f_5 * sng1_528[k]
                   + f_3 * pc_x[k] * snh_738[k];

        t_984[k] = f_3 * pc_y[k] * snh_737[k];

        t_985[k] = f_13 * smh_740[k]
                   + f_4 * sng0_530[k]
                   - f_5 * sng1_530[k]
                   + f_3 * pc_x[k] * snh_740[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_x, pc_y, pc_z, smh_570, smh_741, sng0_531, \
                         sng1_531, snh_738, snh_740, snh_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_13 * smh_741[k]
                   + f_6 * sng0_531[k]
                   - f_7 * sng1_531[k]
                   + f_3 * pc_x[k] * snh_741[k];

        t_987[k] = f_17 * smh_570[k]
                   + f_3 * pc_z[k] * snh_738[k];

        t_988[k] = f_3 * pc_y[k] * snh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_z, smh_573, smh_744, smh_745, sng0_534, \
                         sng0_535, sng1_534, sng1_535, snh_741, snh_744, \
                         snh_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_13 * smh_744[k]
                   + f_6 * sng0_534[k]
                   - f_7 * sng1_534[k]
                   + f_3 * pc_x[k] * snh_744[k];

        t_990[k] = f_13 * smh_745[k]
                   + f_8 * sng0_535[k]
                   - f_9 * sng1_535[k]
                   + f_3 * pc_x[k] * snh_745[k];

        t_991[k] = f_17 * smh_573[k]
                   + f_3 * pc_z[k] * snh_741[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, pc_x, pc_y, smh_747, smh_749, sng0_537, \
                         sng0_539, sng1_537, sng1_539, snh_744, snh_747, \
                         snh_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_13 * smh_747[k]
                   + f_8 * sng0_537[k]
                   - f_9 * sng1_537[k]
                   + f_3 * pc_x[k] * snh_747[k];

        t_993[k] = f_3 * pc_y[k] * snh_744[k];

        t_994[k] = f_13 * smh_749[k]
                   + f_8 * sng0_539[k]
                   - f_9 * sng1_539[k]
                   + f_3 * pc_x[k] * snh_749[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pc_x, smh_750, smh_751, smh_752, \
                         smh_753, smh_754, snh_750, snh_751, snh_752, snh_753, \
                         snh_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_13 * smh_750[k]
                   + f_3 * pc_x[k] * snh_750[k];

        t_996[k] = f_13 * smh_751[k]
                   + f_3 * pc_x[k] * snh_751[k];

        t_997[k] = f_13 * smh_752[k]
                   + f_3 * pc_x[k] * snh_752[k];

        t_998[k] = f_13 * smh_753[k]
                   + f_3 * pc_x[k] * snh_753[k];

        t_999[k] = f_13 * smh_754[k]
                   + f_3 * pc_x[k] * snh_754[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pc_x, pc_y, pc_z, smh_582, smh_755, \
                         sng0_535, sng0_537, sng1_535, sng1_537, snh_750, snh_752, \
                         snh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_13 * smh_755[k]
                    + f_3 * pc_x[k] * snh_755[k];

        t_1001[k] = f_1 * sng0_535[k]
                    - f_2 * sng1_535[k]
                    + f_3 * pc_y[k] * snh_750[k];

        t_1002[k] = f_17 * smh_582[k]
                    + f_3 * pc_z[k] * snh_750[k];

        t_1003[k] = f_4 * sng0_537[k]
                    - f_5 * sng1_537[k]
                    + f_3 * pc_y[k] * snh_752[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, smh_587, sng0_538, \
                         sng0_539, sng1_538, sng1_539, snh_753, snh_754, \
                         snh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_6 * sng0_538[k]
                    - f_7 * sng1_538[k]
                    + f_3 * pc_y[k] * snh_753[k];

        t_1005[k] = f_8 * sng0_539[k]
                    - f_9 * sng1_539[k]
                    + f_3 * pc_y[k] * snh_754[k];

        t_1006[k] = f_3 * pc_y[k] * snh_755[k];

        t_1007[k] = f_17 * smh_587[k]
                    + f_1 * sng0_539[k]
                    - f_2 * sng1_539[k]
                    + f_3 * pc_z[k] * snh_755[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, smh_588, smh_756, \
                         smh_759, sng0_540, sng0_543, sng1_540, sng1_543, snh_756, \
                         snh_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_12 * smh_756[k]
                    + f_1 * sng0_540[k]
                    - f_2 * sng1_540[k]
                    + f_3 * pc_x[k] * snh_756[k];

        t_1009[k] = f_16 * smh_588[k]
                    + f_3 * pc_y[k] * snh_756[k];

        t_1010[k] = f_3 * pc_z[k] * snh_756[k];

        t_1011[k] = f_12 * smh_759[k]
                    + f_4 * sng0_543[k]
                    - f_5 * sng1_543[k]
                    + f_3 * pc_x[k] * snh_759[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pc_x, pc_y, smh_590, smh_761, smh_762, \
                         sng0_545, sng0_546, sng1_545, sng1_546, snh_758, snh_761, \
                         snh_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_16 * smh_590[k]
                    + f_3 * pc_y[k] * snh_758[k];

        t_1013[k] = f_12 * smh_761[k]
                    + f_4 * sng0_545[k]
                    - f_5 * sng1_545[k]
                    + f_3 * pc_x[k] * snh_761[k];

        t_1014[k] = f_12 * smh_762[k]
                    + f_6 * sng0_546[k]
                    - f_7 * sng1_546[k]
                    + f_3 * pc_x[k] * snh_762[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, pc_x, pc_y, pc_z, smh_593, smh_765, sng0_549, \
                         sng1_549, snh_759, snh_761, snh_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = f_3 * pc_z[k] * snh_759[k];

        t_1016[k] = f_16 * smh_593[k]
                    + f_3 * pc_y[k] * snh_761[k];

        t_1017[k] = f_12 * smh_765[k]
                    + f_6 * sng0_549[k]
                    - f_7 * sng1_549[k]
                    + f_3 * pc_x[k] * snh_765[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, pc_x, pc_z, smh_766, smh_768, sng0_550, \
                         sng0_552, sng1_550, sng1_552, snh_762, snh_766, \
                         snh_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_12 * smh_766[k]
                    + f_8 * sng0_550[k]
                    - f_9 * sng1_550[k]
                    + f_3 * pc_x[k] * snh_766[k];

        t_1019[k] = f_3 * pc_z[k] * snh_762[k];

        t_1020[k] = f_12 * smh_768[k]
                    + f_8 * sng0_552[k]
                    - f_9 * sng1_552[k]
                    + f_3 * pc_x[k] * snh_768[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smi0,
                                                          const size_t smh, const size_t smi1,
                                                          const size_t sng0, const size_t sng1,
                                                          const size_t snh, const size_t ncols,
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
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi0_784 = buffer.data(smi0 + 784);
    const auto *smi0_787 = buffer.data(smi0 + 787);
    const auto *smi0_790 = buffer.data(smi0 + 790);
    const auto *smi0_794 = buffer.data(smi0 + 794);
    const auto *smi0_796 = buffer.data(smi0 + 796);
    const auto *smi0_805 = buffer.data(smi0 + 805);

    const auto *smh_588 = buffer.data(smh + 588);
    const auto *smh_591 = buffer.data(smh + 591);
    const auto *smh_594 = buffer.data(smh + 594);
    const auto *smh_595 = buffer.data(smh + 595);
    const auto *smh_597 = buffer.data(smh + 597);
    const auto *smh_603 = buffer.data(smh + 603);
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
    const auto *smh_624 = buffer.data(smh + 624);
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
    const auto *smh_645 = buffer.data(smh + 645);
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
    const auto *smh_666 = buffer.data(smh + 666);
    const auto *smh_668 = buffer.data(smh + 668);
    const auto *smh_669 = buffer.data(smh + 669);
    const auto *smh_670 = buffer.data(smh + 670);
    const auto *smh_671 = buffer.data(smh + 671);
    const auto *smh_672 = buffer.data(smh + 672);
    const auto *smh_674 = buffer.data(smh + 674);
    const auto *smh_677 = buffer.data(smh + 677);
    const auto *smh_770 = buffer.data(smh + 770);
    const auto *smh_771 = buffer.data(smh + 771);
    const auto *smh_772 = buffer.data(smh + 772);
    const auto *smh_773 = buffer.data(smh + 773);
    const auto *smh_774 = buffer.data(smh + 774);
    const auto *smh_775 = buffer.data(smh + 775);
    const auto *smh_776 = buffer.data(smh + 776);
    const auto *smh_782 = buffer.data(smh + 782);
    const auto *smh_786 = buffer.data(smh + 786);
    const auto *smh_791 = buffer.data(smh + 791);
    const auto *smh_792 = buffer.data(smh + 792);
    const auto *smh_793 = buffer.data(smh + 793);
    const auto *smh_794 = buffer.data(smh + 794);
    const auto *smh_795 = buffer.data(smh + 795);
    const auto *smh_796 = buffer.data(smh + 796);
    const auto *smh_797 = buffer.data(smh + 797);
    const auto *smh_798 = buffer.data(smh + 798);
    const auto *smh_801 = buffer.data(smh + 801);
    const auto *smh_803 = buffer.data(smh + 803);
    const auto *smh_804 = buffer.data(smh + 804);
    const auto *smh_807 = buffer.data(smh + 807);
    const auto *smh_808 = buffer.data(smh + 808);
    const auto *smh_810 = buffer.data(smh + 810);
    const auto *smh_812 = buffer.data(smh + 812);
    const auto *smh_813 = buffer.data(smh + 813);
    const auto *smh_814 = buffer.data(smh + 814);
    const auto *smh_815 = buffer.data(smh + 815);
    const auto *smh_816 = buffer.data(smh + 816);
    const auto *smh_817 = buffer.data(smh + 817);
    const auto *smh_818 = buffer.data(smh + 818);
    const auto *smh_819 = buffer.data(smh + 819);
    const auto *smh_822 = buffer.data(smh + 822);
    const auto *smh_824 = buffer.data(smh + 824);
    const auto *smh_825 = buffer.data(smh + 825);
    const auto *smh_828 = buffer.data(smh + 828);
    const auto *smh_829 = buffer.data(smh + 829);
    const auto *smh_831 = buffer.data(smh + 831);
    const auto *smh_833 = buffer.data(smh + 833);
    const auto *smh_834 = buffer.data(smh + 834);
    const auto *smh_835 = buffer.data(smh + 835);
    const auto *smh_836 = buffer.data(smh + 836);
    const auto *smh_837 = buffer.data(smh + 837);
    const auto *smh_838 = buffer.data(smh + 838);
    const auto *smh_839 = buffer.data(smh + 839);
    const auto *smh_840 = buffer.data(smh + 840);
    const auto *smh_843 = buffer.data(smh + 843);
    const auto *smh_845 = buffer.data(smh + 845);
    const auto *smh_846 = buffer.data(smh + 846);
    const auto *smh_849 = buffer.data(smh + 849);
    const auto *smh_850 = buffer.data(smh + 850);

    const auto *smi1_784 = buffer.data(smi1 + 784);
    const auto *smi1_787 = buffer.data(smi1 + 787);
    const auto *smi1_790 = buffer.data(smi1 + 790);
    const auto *smi1_794 = buffer.data(smi1 + 794);
    const auto *smi1_796 = buffer.data(smi1 + 796);
    const auto *smi1_805 = buffer.data(smi1 + 805);

    const auto *sng0_550 = buffer.data(sng0 + 550);
    const auto *sng0_552 = buffer.data(sng0 + 552);
    const auto *sng0_553 = buffer.data(sng0 + 553);
    const auto *sng0_554 = buffer.data(sng0 + 554);
    const auto *sng0_560 = buffer.data(sng0 + 560);
    const auto *sng0_564 = buffer.data(sng0 + 564);
    const auto *sng0_567 = buffer.data(sng0 + 567);
    const auto *sng0_568 = buffer.data(sng0 + 568);
    const auto *sng0_569 = buffer.data(sng0 + 569);
    const auto *sng0_570 = buffer.data(sng0 + 570);
    const auto *sng0_573 = buffer.data(sng0 + 573);
    const auto *sng0_575 = buffer.data(sng0 + 575);
    const auto *sng0_576 = buffer.data(sng0 + 576);
    const auto *sng0_579 = buffer.data(sng0 + 579);
    const auto *sng0_580 = buffer.data(sng0 + 580);
    const auto *sng0_582 = buffer.data(sng0 + 582);
    const auto *sng0_583 = buffer.data(sng0 + 583);
    const auto *sng0_584 = buffer.data(sng0 + 584);
    const auto *sng0_585 = buffer.data(sng0 + 585);
    const auto *sng0_588 = buffer.data(sng0 + 588);
    const auto *sng0_590 = buffer.data(sng0 + 590);
    const auto *sng0_591 = buffer.data(sng0 + 591);
    const auto *sng0_594 = buffer.data(sng0 + 594);
    const auto *sng0_595 = buffer.data(sng0 + 595);
    const auto *sng0_597 = buffer.data(sng0 + 597);
    const auto *sng0_598 = buffer.data(sng0 + 598);
    const auto *sng0_599 = buffer.data(sng0 + 599);
    const auto *sng0_600 = buffer.data(sng0 + 600);
    const auto *sng0_603 = buffer.data(sng0 + 603);
    const auto *sng0_605 = buffer.data(sng0 + 605);
    const auto *sng0_606 = buffer.data(sng0 + 606);
    const auto *sng0_609 = buffer.data(sng0 + 609);
    const auto *sng0_610 = buffer.data(sng0 + 610);

    const auto *sng1_550 = buffer.data(sng1 + 550);
    const auto *sng1_552 = buffer.data(sng1 + 552);
    const auto *sng1_553 = buffer.data(sng1 + 553);
    const auto *sng1_554 = buffer.data(sng1 + 554);
    const auto *sng1_560 = buffer.data(sng1 + 560);
    const auto *sng1_564 = buffer.data(sng1 + 564);
    const auto *sng1_567 = buffer.data(sng1 + 567);
    const auto *sng1_568 = buffer.data(sng1 + 568);
    const auto *sng1_569 = buffer.data(sng1 + 569);
    const auto *sng1_570 = buffer.data(sng1 + 570);
    const auto *sng1_573 = buffer.data(sng1 + 573);
    const auto *sng1_575 = buffer.data(sng1 + 575);
    const auto *sng1_576 = buffer.data(sng1 + 576);
    const auto *sng1_579 = buffer.data(sng1 + 579);
    const auto *sng1_580 = buffer.data(sng1 + 580);
    const auto *sng1_582 = buffer.data(sng1 + 582);
    const auto *sng1_583 = buffer.data(sng1 + 583);
    const auto *sng1_584 = buffer.data(sng1 + 584);
    const auto *sng1_585 = buffer.data(sng1 + 585);
    const auto *sng1_588 = buffer.data(sng1 + 588);
    const auto *sng1_590 = buffer.data(sng1 + 590);
    const auto *sng1_591 = buffer.data(sng1 + 591);
    const auto *sng1_594 = buffer.data(sng1 + 594);
    const auto *sng1_595 = buffer.data(sng1 + 595);
    const auto *sng1_597 = buffer.data(sng1 + 597);
    const auto *sng1_598 = buffer.data(sng1 + 598);
    const auto *sng1_599 = buffer.data(sng1 + 599);
    const auto *sng1_600 = buffer.data(sng1 + 600);
    const auto *sng1_603 = buffer.data(sng1 + 603);
    const auto *sng1_605 = buffer.data(sng1 + 605);
    const auto *sng1_606 = buffer.data(sng1 + 606);
    const auto *sng1_609 = buffer.data(sng1 + 609);
    const auto *sng1_610 = buffer.data(sng1 + 610);

    const auto *snh_765 = buffer.data(snh + 765);
    const auto *snh_770 = buffer.data(snh + 770);
    const auto *snh_771 = buffer.data(snh + 771);
    const auto *snh_772 = buffer.data(snh + 772);
    const auto *snh_773 = buffer.data(snh + 773);
    const auto *snh_774 = buffer.data(snh + 774);
    const auto *snh_775 = buffer.data(snh + 775);
    const auto *snh_776 = buffer.data(snh + 776);
    const auto *snh_777 = buffer.data(snh + 777);
    const auto *snh_779 = buffer.data(snh + 779);
    const auto *snh_780 = buffer.data(snh + 780);
    const auto *snh_782 = buffer.data(snh + 782);
    const auto *snh_783 = buffer.data(snh + 783);
    const auto *snh_786 = buffer.data(snh + 786);
    const auto *snh_791 = buffer.data(snh + 791);
    const auto *snh_792 = buffer.data(snh + 792);
    const auto *snh_793 = buffer.data(snh + 793);
    const auto *snh_794 = buffer.data(snh + 794);
    const auto *snh_795 = buffer.data(snh + 795);
    const auto *snh_796 = buffer.data(snh + 796);
    const auto *snh_797 = buffer.data(snh + 797);
    const auto *snh_798 = buffer.data(snh + 798);
    const auto *snh_800 = buffer.data(snh + 800);
    const auto *snh_801 = buffer.data(snh + 801);
    const auto *snh_803 = buffer.data(snh + 803);
    const auto *snh_804 = buffer.data(snh + 804);
    const auto *snh_807 = buffer.data(snh + 807);
    const auto *snh_808 = buffer.data(snh + 808);
    const auto *snh_810 = buffer.data(snh + 810);
    const auto *snh_812 = buffer.data(snh + 812);
    const auto *snh_813 = buffer.data(snh + 813);
    const auto *snh_814 = buffer.data(snh + 814);
    const auto *snh_815 = buffer.data(snh + 815);
    const auto *snh_816 = buffer.data(snh + 816);
    const auto *snh_817 = buffer.data(snh + 817);
    const auto *snh_818 = buffer.data(snh + 818);
    const auto *snh_819 = buffer.data(snh + 819);
    const auto *snh_821 = buffer.data(snh + 821);
    const auto *snh_822 = buffer.data(snh + 822);
    const auto *snh_824 = buffer.data(snh + 824);
    const auto *snh_825 = buffer.data(snh + 825);
    const auto *snh_828 = buffer.data(snh + 828);
    const auto *snh_829 = buffer.data(snh + 829);
    const auto *snh_831 = buffer.data(snh + 831);
    const auto *snh_833 = buffer.data(snh + 833);
    const auto *snh_834 = buffer.data(snh + 834);
    const auto *snh_835 = buffer.data(snh + 835);
    const auto *snh_836 = buffer.data(snh + 836);
    const auto *snh_837 = buffer.data(snh + 837);
    const auto *snh_838 = buffer.data(snh + 838);
    const auto *snh_839 = buffer.data(snh + 839);
    const auto *snh_840 = buffer.data(snh + 840);
    const auto *snh_842 = buffer.data(snh + 842);
    const auto *snh_843 = buffer.data(snh + 843);
    const auto *snh_845 = buffer.data(snh + 845);
    const auto *snh_846 = buffer.data(snh + 846);
    const auto *snh_849 = buffer.data(snh + 849);
    const auto *snh_850 = buffer.data(snh + 850);

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pc_x, pc_y, smh_597, smh_770, \
                         smh_771, smh_772, sng0_554, sng1_554, snh_765, snh_770, snh_771, \
                         snh_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_16 * smh_597[k]
                    + f_3 * pc_y[k] * snh_765[k];

        t_1022[k] = f_12 * smh_770[k]
                    + f_8 * sng0_554[k]
                    - f_9 * sng1_554[k]
                    + f_3 * pc_x[k] * snh_770[k];

        t_1023[k] = f_12 * smh_771[k]
                    + f_3 * pc_x[k] * snh_771[k];

        t_1024[k] = f_12 * smh_772[k]
                    + f_3 * pc_x[k] * snh_772[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pc_x, smh_773, smh_774, smh_775, \
                         smh_776, snh_773, snh_774, snh_775, snh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_12 * smh_773[k]
                    + f_3 * pc_x[k] * snh_773[k];

        t_1026[k] = f_12 * smh_774[k]
                    + f_3 * pc_x[k] * snh_774[k];

        t_1027[k] = f_12 * smh_775[k]
                    + f_3 * pc_x[k] * snh_775[k];

        t_1028[k] = f_12 * smh_776[k]
                    + f_3 * pc_x[k] * snh_776[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, pc_y, pc_z, smh_603, smh_605, sng0_550, \
                         sng0_552, sng1_550, sng1_552, snh_771, \
                         snh_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = f_16 * smh_603[k]
                    + f_1 * sng0_550[k]
                    - f_2 * sng1_550[k]
                    + f_3 * pc_y[k] * snh_771[k];

        t_1030[k] = f_3 * pc_z[k] * snh_771[k];

        t_1031[k] = f_16 * smh_605[k]
                    + f_4 * sng0_552[k]
                    - f_5 * sng1_552[k]
                    + f_3 * pc_y[k] * snh_773[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, t_1035, pc_y, pc_z, smh_606, smh_607, \
                         smh_608, sng0_553, sng0_554, sng1_553, sng1_554, snh_774, snh_775, \
                         snh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = f_16 * smh_606[k]
                    + f_6 * sng0_553[k]
                    - f_7 * sng1_553[k]
                    + f_3 * pc_y[k] * snh_774[k];

        t_1033[k] = f_16 * smh_607[k]
                    + f_8 * sng0_554[k]
                    - f_9 * sng1_554[k]
                    + f_3 * pc_y[k] * snh_775[k];

        t_1034[k] = f_16 * smh_608[k]
                    + f_3 * pc_y[k] * snh_776[k];

        t_1035[k] = f_1 * sng0_554[k]
                    - f_2 * sng1_554[k]
                    + f_3 * pc_z[k] * snh_776[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, t_1039, pb_z, pc_y, pc_z, smi0_784, smi0_787, \
                         smh_588, smh_609, smi1_784, smi1_787, \
                         snh_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = pb_z[k] * smi0_784[k]
                    - f_10 * pc_z[k] * smi1_784[k];

        t_1037[k] = f_17 * smh_609[k]
                    + f_3 * pc_y[k] * snh_777[k];

        t_1038[k] = f_11 * smh_588[k]
                    + f_3 * pc_z[k] * snh_777[k];

        t_1039[k] = pb_z[k] * smi0_787[k]
                    - f_10 * pc_z[k] * smi1_787[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pb_z, pc_x, pc_y, pc_z, smi0_790, smh_611, \
                         smh_782, smi1_790, sng0_560, sng1_560, snh_779, \
                         snh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_17 * smh_611[k]
                    + f_3 * pc_y[k] * snh_779[k];

        t_1041[k] = f_12 * smh_782[k]
                    + f_4 * sng0_560[k]
                    - f_5 * sng1_560[k]
                    + f_3 * pc_x[k] * snh_782[k];

        t_1042[k] = pb_z[k] * smi0_790[k]
                    - f_10 * pc_z[k] * smi1_790[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pc_x, pc_y, pc_z, smh_591, smh_614, smh_786, \
                         sng0_564, sng1_564, snh_780, snh_782, \
                         snh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_11 * smh_591[k]
                    + f_3 * pc_z[k] * snh_780[k];

        t_1044[k] = f_17 * smh_614[k]
                    + f_3 * pc_y[k] * snh_782[k];

        t_1045[k] = f_12 * smh_786[k]
                    + f_6 * sng0_564[k]
                    - f_7 * sng1_564[k]
                    + f_3 * pc_x[k] * snh_786[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, t_1049, pb_z, pc_y, pc_z, smi0_794, smi0_796, \
                         smh_594, smh_595, smh_618, smi1_794, smi1_796, snh_783, \
                         snh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pb_z[k] * smi0_794[k]
                    - f_10 * pc_z[k] * smi1_794[k];

        t_1047[k] = f_11 * smh_594[k]
                    + f_3 * pc_z[k] * snh_783[k];

        t_1048[k] = pb_z[k] * smi0_796[k]
                    + f_12 * smh_595[k]
                    - f_10 * pc_z[k] * smi1_796[k];

        t_1049[k] = f_17 * smh_618[k]
                    + f_3 * pc_y[k] * snh_786[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pc_x, smh_791, smh_792, smh_793, \
                         smh_794, sng0_569, sng1_569, snh_791, snh_792, snh_793, \
                         snh_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_12 * smh_791[k]
                    + f_8 * sng0_569[k]
                    - f_9 * sng1_569[k]
                    + f_3 * pc_x[k] * snh_791[k];

        t_1051[k] = f_12 * smh_792[k]
                    + f_3 * pc_x[k] * snh_792[k];

        t_1052[k] = f_12 * smh_793[k]
                    + f_3 * pc_x[k] * snh_793[k];

        t_1053[k] = f_12 * smh_794[k]
                    + f_3 * pc_x[k] * snh_794[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, pb_z, pc_x, pc_z, smi0_805, smh_795, \
                         smh_796, smh_797, smi1_805, snh_795, snh_796, \
                         snh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_12 * smh_795[k]
                    + f_3 * pc_x[k] * snh_795[k];

        t_1055[k] = f_12 * smh_796[k]
                    + f_3 * pc_x[k] * snh_796[k];

        t_1056[k] = f_12 * smh_797[k]
                    + f_3 * pc_x[k] * snh_797[k];

        t_1057[k] = pb_z[k] * smi0_805[k]
                    - f_10 * pc_z[k] * smi1_805[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_y, pc_z, smh_603, smh_626, smh_627, \
                         sng0_567, sng0_568, sng1_567, sng1_568, snh_792, snh_794, \
                         snh_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_11 * smh_603[k]
                    + f_3 * pc_z[k] * snh_792[k];

        t_1059[k] = f_17 * smh_626[k]
                    + f_4 * sng0_567[k]
                    - f_5 * sng1_567[k]
                    + f_3 * pc_y[k] * snh_794[k];

        t_1060[k] = f_17 * smh_627[k]
                    + f_6 * sng0_568[k]
                    - f_7 * sng1_568[k]
                    + f_3 * pc_y[k] * snh_795[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pc_y, pc_z, smh_608, smh_628, smh_629, \
                         sng0_569, sng1_569, snh_796, snh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_17 * smh_628[k]
                    + f_8 * sng0_569[k]
                    - f_9 * sng1_569[k]
                    + f_3 * pc_y[k] * snh_796[k];

        t_1062[k] = f_17 * smh_629[k]
                    + f_3 * pc_y[k] * snh_797[k];

        t_1063[k] = f_11 * smh_608[k]
                    + f_1 * sng0_569[k]
                    - f_2 * sng1_569[k]
                    + f_3 * pc_z[k] * snh_797[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, pc_x, pc_y, pc_z, smh_609, smh_630, smh_798, \
                         sng0_570, sng1_570, snh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_12 * smh_798[k]
                    + f_1 * sng0_570[k]
                    - f_2 * sng1_570[k]
                    + f_3 * pc_x[k] * snh_798[k];

        t_1065[k] = f_18 * smh_630[k]
                    + f_3 * pc_y[k] * snh_798[k];

        t_1066[k] = f_12 * smh_609[k]
                    + f_3 * pc_z[k] * snh_798[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pc_x, pc_y, smh_632, smh_801, smh_803, \
                         sng0_573, sng0_575, sng1_573, sng1_575, snh_800, snh_801, \
                         snh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_12 * smh_801[k]
                    + f_4 * sng0_573[k]
                    - f_5 * sng1_573[k]
                    + f_3 * pc_x[k] * snh_801[k];

        t_1068[k] = f_18 * smh_632[k]
                    + f_3 * pc_y[k] * snh_800[k];

        t_1069[k] = f_12 * smh_803[k]
                    + f_4 * sng0_575[k]
                    - f_5 * sng1_575[k]
                    + f_3 * pc_x[k] * snh_803[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, smh_612, smh_635, smh_804, \
                         sng0_576, sng1_576, snh_801, snh_803, \
                         snh_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_12 * smh_804[k]
                    + f_6 * sng0_576[k]
                    - f_7 * sng1_576[k]
                    + f_3 * pc_x[k] * snh_804[k];

        t_1071[k] = f_12 * smh_612[k]
                    + f_3 * pc_z[k] * snh_801[k];

        t_1072[k] = f_18 * smh_635[k]
                    + f_3 * pc_y[k] * snh_803[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, pc_z, smh_615, smh_807, smh_808, \
                         sng0_579, sng0_580, sng1_579, sng1_580, snh_804, snh_807, \
                         snh_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_12 * smh_807[k]
                    + f_6 * sng0_579[k]
                    - f_7 * sng1_579[k]
                    + f_3 * pc_x[k] * snh_807[k];

        t_1074[k] = f_12 * smh_808[k]
                    + f_8 * sng0_580[k]
                    - f_9 * sng1_580[k]
                    + f_3 * pc_x[k] * snh_808[k];

        t_1075[k] = f_12 * smh_615[k]
                    + f_3 * pc_z[k] * snh_804[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, pc_y, smh_639, smh_810, smh_812, \
                         sng0_582, sng0_584, sng1_582, sng1_584, snh_807, snh_810, \
                         snh_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_12 * smh_810[k]
                    + f_8 * sng0_582[k]
                    - f_9 * sng1_582[k]
                    + f_3 * pc_x[k] * snh_810[k];

        t_1077[k] = f_18 * smh_639[k]
                    + f_3 * pc_y[k] * snh_807[k];

        t_1078[k] = f_12 * smh_812[k]
                    + f_8 * sng0_584[k]
                    - f_9 * sng1_584[k]
                    + f_3 * pc_x[k] * snh_812[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, t_1083, pc_x, smh_813, smh_814, \
                         smh_815, smh_816, smh_817, snh_813, snh_814, snh_815, snh_816, \
                         snh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_12 * smh_813[k]
                    + f_3 * pc_x[k] * snh_813[k];

        t_1080[k] = f_12 * smh_814[k]
                    + f_3 * pc_x[k] * snh_814[k];

        t_1081[k] = f_12 * smh_815[k]
                    + f_3 * pc_x[k] * snh_815[k];

        t_1082[k] = f_12 * smh_816[k]
                    + f_3 * pc_x[k] * snh_816[k];

        t_1083[k] = f_12 * smh_817[k]
                    + f_3 * pc_x[k] * snh_817[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, pc_x, pc_y, pc_z, smh_624, smh_645, smh_818, \
                         sng0_580, sng1_580, snh_813, snh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_12 * smh_818[k]
                    + f_3 * pc_x[k] * snh_818[k];

        t_1085[k] = f_18 * smh_645[k]
                    + f_1 * sng0_580[k]
                    - f_2 * sng1_580[k]
                    + f_3 * pc_y[k] * snh_813[k];

        t_1086[k] = f_12 * smh_624[k]
                    + f_3 * pc_z[k] * snh_813[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, smh_647, smh_648, smh_649, sng0_582, \
                         sng0_583, sng0_584, sng1_582, sng1_583, sng1_584, snh_815, snh_816, \
                         snh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_18 * smh_647[k]
                    + f_4 * sng0_582[k]
                    - f_5 * sng1_582[k]
                    + f_3 * pc_y[k] * snh_815[k];

        t_1088[k] = f_18 * smh_648[k]
                    + f_6 * sng0_583[k]
                    - f_7 * sng1_583[k]
                    + f_3 * pc_y[k] * snh_816[k];

        t_1089[k] = f_18 * smh_649[k]
                    + f_8 * sng0_584[k]
                    - f_9 * sng1_584[k]
                    + f_3 * pc_y[k] * snh_817[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, smh_629, smh_650, smh_819, \
                         sng0_584, sng0_585, sng1_584, sng1_585, snh_818, \
                         snh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_18 * smh_650[k]
                    + f_3 * pc_y[k] * snh_818[k];

        t_1091[k] = f_12 * smh_629[k]
                    + f_1 * sng0_584[k]
                    - f_2 * sng1_584[k]
                    + f_3 * pc_z[k] * snh_818[k];

        t_1092[k] = f_12 * smh_819[k]
                    + f_1 * sng0_585[k]
                    - f_2 * sng1_585[k]
                    + f_3 * pc_x[k] * snh_819[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pc_x, pc_y, pc_z, smh_630, smh_651, \
                         smh_653, smh_822, sng0_588, sng1_588, snh_819, snh_821, \
                         snh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_19 * smh_651[k]
                    + f_3 * pc_y[k] * snh_819[k];

        t_1094[k] = f_13 * smh_630[k]
                    + f_3 * pc_z[k] * snh_819[k];

        t_1095[k] = f_12 * smh_822[k]
                    + f_4 * sng0_588[k]
                    - f_5 * sng1_588[k]
                    + f_3 * pc_x[k] * snh_822[k];

        t_1096[k] = f_19 * smh_653[k]
                    + f_3 * pc_y[k] * snh_821[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_z, smh_633, smh_824, smh_825, \
                         sng0_590, sng0_591, sng1_590, sng1_591, snh_822, snh_824, \
                         snh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_12 * smh_824[k]
                    + f_4 * sng0_590[k]
                    - f_5 * sng1_590[k]
                    + f_3 * pc_x[k] * snh_824[k];

        t_1098[k] = f_12 * smh_825[k]
                    + f_6 * sng0_591[k]
                    - f_7 * sng1_591[k]
                    + f_3 * pc_x[k] * snh_825[k];

        t_1099[k] = f_13 * smh_633[k]
                    + f_3 * pc_z[k] * snh_822[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_y, smh_656, smh_828, smh_829, \
                         sng0_594, sng0_595, sng1_594, sng1_595, snh_824, snh_828, \
                         snh_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_19 * smh_656[k]
                    + f_3 * pc_y[k] * snh_824[k];

        t_1101[k] = f_12 * smh_828[k]
                    + f_6 * sng0_594[k]
                    - f_7 * sng1_594[k]
                    + f_3 * pc_x[k] * snh_828[k];

        t_1102[k] = f_12 * smh_829[k]
                    + f_8 * sng0_595[k]
                    - f_9 * sng1_595[k]
                    + f_3 * pc_x[k] * snh_829[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, pc_y, pc_z, smh_636, smh_660, smh_831, \
                         sng0_597, sng1_597, snh_825, snh_828, \
                         snh_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_13 * smh_636[k]
                    + f_3 * pc_z[k] * snh_825[k];

        t_1104[k] = f_12 * smh_831[k]
                    + f_8 * sng0_597[k]
                    - f_9 * sng1_597[k]
                    + f_3 * pc_x[k] * snh_831[k];

        t_1105[k] = f_19 * smh_660[k]
                    + f_3 * pc_y[k] * snh_828[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, smh_833, smh_834, smh_835, \
                         smh_836, sng0_599, sng1_599, snh_833, snh_834, snh_835, \
                         snh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_12 * smh_833[k]
                    + f_8 * sng0_599[k]
                    - f_9 * sng1_599[k]
                    + f_3 * pc_x[k] * snh_833[k];

        t_1107[k] = f_12 * smh_834[k]
                    + f_3 * pc_x[k] * snh_834[k];

        t_1108[k] = f_12 * smh_835[k]
                    + f_3 * pc_x[k] * snh_835[k];

        t_1109[k] = f_12 * smh_836[k]
                    + f_3 * pc_x[k] * snh_836[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pc_x, pc_y, smh_666, smh_837, \
                         smh_838, smh_839, sng0_595, sng1_595, snh_834, snh_837, snh_838, \
                         snh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_12 * smh_837[k]
                    + f_3 * pc_x[k] * snh_837[k];

        t_1111[k] = f_12 * smh_838[k]
                    + f_3 * pc_x[k] * snh_838[k];

        t_1112[k] = f_12 * smh_839[k]
                    + f_3 * pc_x[k] * snh_839[k];

        t_1113[k] = f_19 * smh_666[k]
                    + f_1 * sng0_595[k]
                    - f_2 * sng1_595[k]
                    + f_3 * pc_y[k] * snh_834[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_y, pc_z, smh_645, smh_668, smh_669, \
                         sng0_597, sng0_598, sng1_597, sng1_598, snh_834, snh_836, \
                         snh_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_13 * smh_645[k]
                    + f_3 * pc_z[k] * snh_834[k];

        t_1115[k] = f_19 * smh_668[k]
                    + f_4 * sng0_597[k]
                    - f_5 * sng1_597[k]
                    + f_3 * pc_y[k] * snh_836[k];

        t_1116[k] = f_19 * smh_669[k]
                    + f_6 * sng0_598[k]
                    - f_7 * sng1_598[k]
                    + f_3 * pc_y[k] * snh_837[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, pc_y, pc_z, smh_650, smh_670, smh_671, \
                         sng0_599, sng1_599, snh_838, snh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_19 * smh_670[k]
                    + f_8 * sng0_599[k]
                    - f_9 * sng1_599[k]
                    + f_3 * pc_y[k] * snh_838[k];

        t_1118[k] = f_19 * smh_671[k]
                    + f_3 * pc_y[k] * snh_839[k];

        t_1119[k] = f_13 * smh_650[k]
                    + f_1 * sng0_599[k]
                    - f_2 * sng1_599[k]
                    + f_3 * pc_z[k] * snh_839[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, pc_x, pc_y, pc_z, smh_651, smh_672, smh_840, \
                         sng0_600, sng1_600, snh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = f_12 * smh_840[k]
                    + f_1 * sng0_600[k]
                    - f_2 * sng1_600[k]
                    + f_3 * pc_x[k] * snh_840[k];

        t_1121[k] = f_14 * smh_672[k]
                    + f_3 * pc_y[k] * snh_840[k];

        t_1122[k] = f_14 * smh_651[k]
                    + f_3 * pc_z[k] * snh_840[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, pc_x, pc_y, smh_674, smh_843, smh_845, \
                         sng0_603, sng0_605, sng1_603, sng1_605, snh_842, snh_843, \
                         snh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_12 * smh_843[k]
                    + f_4 * sng0_603[k]
                    - f_5 * sng1_603[k]
                    + f_3 * pc_x[k] * snh_843[k];

        t_1124[k] = f_14 * smh_674[k]
                    + f_3 * pc_y[k] * snh_842[k];

        t_1125[k] = f_12 * smh_845[k]
                    + f_4 * sng0_605[k]
                    - f_5 * sng1_605[k]
                    + f_3 * pc_x[k] * snh_845[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, pc_x, pc_y, pc_z, smh_654, smh_677, smh_846, \
                         sng0_606, sng1_606, snh_843, snh_845, \
                         snh_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_12 * smh_846[k]
                    + f_6 * sng0_606[k]
                    - f_7 * sng1_606[k]
                    + f_3 * pc_x[k] * snh_846[k];

        t_1127[k] = f_14 * smh_654[k]
                    + f_3 * pc_z[k] * snh_843[k];

        t_1128[k] = f_14 * smh_677[k]
                    + f_3 * pc_y[k] * snh_845[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pc_x, pc_z, smh_657, smh_849, smh_850, \
                         sng0_609, sng0_610, sng1_609, sng1_610, snh_846, snh_849, \
                         snh_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_12 * smh_849[k]
                    + f_6 * sng0_609[k]
                    - f_7 * sng1_609[k]
                    + f_3 * pc_x[k] * snh_849[k];

        t_1130[k] = f_12 * smh_850[k]
                    + f_8 * sng0_610[k]
                    - f_9 * sng1_610[k]
                    + f_3 * pc_x[k] * snh_850[k];

        t_1131[k] = f_14 * smh_657[k]
                    + f_3 * pc_z[k] * snh_846[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smi0,
                                                           const size_t smh, const size_t smi1,
                                                           const size_t sng0, const size_t sng1,
                                                           const size_t snh, const size_t ncols,
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
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi0_980 = buffer.data(smi0 + 980);
    const auto *smi0_983 = buffer.data(smi0 + 983);
    const auto *smi0_985 = buffer.data(smi0 + 985);
    const auto *smi0_986 = buffer.data(smi0 + 986);
    const auto *smi0_989 = buffer.data(smi0 + 989);
    const auto *smi0_990 = buffer.data(smi0 + 990);
    const auto *smi0_992 = buffer.data(smi0 + 992);
    const auto *smi0_994 = buffer.data(smi0 + 994);
    const auto *smi0_1007 = buffer.data(smi0 + 1007);

    const auto *smh_666 = buffer.data(smh + 666);
    const auto *smh_671 = buffer.data(smh + 671);
    const auto *smh_672 = buffer.data(smh + 672);
    const auto *smh_675 = buffer.data(smh + 675);
    const auto *smh_678 = buffer.data(smh + 678);
    const auto *smh_681 = buffer.data(smh + 681);
    const auto *smh_687 = buffer.data(smh + 687);
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
    const auto *smh_708 = buffer.data(smh + 708);
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
    const auto *smh_731 = buffer.data(smh + 731);
    const auto *smh_732 = buffer.data(smh + 732);
    const auto *smh_733 = buffer.data(smh + 733);
    const auto *smh_734 = buffer.data(smh + 734);
    const auto *smh_735 = buffer.data(smh + 735);
    const auto *smh_736 = buffer.data(smh + 736);
    const auto *smh_737 = buffer.data(smh + 737);
    const auto *smh_738 = buffer.data(smh + 738);
    const auto *smh_740 = buffer.data(smh + 740);
    const auto *smh_741 = buffer.data(smh + 741);
    const auto *smh_743 = buffer.data(smh + 743);
    const auto *smh_744 = buffer.data(smh + 744);
    const auto *smh_750 = buffer.data(smh + 750);
    const auto *smh_752 = buffer.data(smh + 752);
    const auto *smh_753 = buffer.data(smh + 753);
    const auto *smh_754 = buffer.data(smh + 754);
    const auto *smh_755 = buffer.data(smh + 755);
    const auto *smh_852 = buffer.data(smh + 852);
    const auto *smh_854 = buffer.data(smh + 854);
    const auto *smh_855 = buffer.data(smh + 855);
    const auto *smh_856 = buffer.data(smh + 856);
    const auto *smh_857 = buffer.data(smh + 857);
    const auto *smh_858 = buffer.data(smh + 858);
    const auto *smh_859 = buffer.data(smh + 859);
    const auto *smh_860 = buffer.data(smh + 860);
    const auto *smh_861 = buffer.data(smh + 861);
    const auto *smh_864 = buffer.data(smh + 864);
    const auto *smh_866 = buffer.data(smh + 866);
    const auto *smh_867 = buffer.data(smh + 867);
    const auto *smh_870 = buffer.data(smh + 870);
    const auto *smh_871 = buffer.data(smh + 871);
    const auto *smh_873 = buffer.data(smh + 873);
    const auto *smh_875 = buffer.data(smh + 875);
    const auto *smh_876 = buffer.data(smh + 876);
    const auto *smh_877 = buffer.data(smh + 877);
    const auto *smh_878 = buffer.data(smh + 878);
    const auto *smh_879 = buffer.data(smh + 879);
    const auto *smh_880 = buffer.data(smh + 880);
    const auto *smh_881 = buffer.data(smh + 881);
    const auto *smh_882 = buffer.data(smh + 882);
    const auto *smh_885 = buffer.data(smh + 885);
    const auto *smh_887 = buffer.data(smh + 887);
    const auto *smh_888 = buffer.data(smh + 888);
    const auto *smh_891 = buffer.data(smh + 891);
    const auto *smh_892 = buffer.data(smh + 892);
    const auto *smh_894 = buffer.data(smh + 894);
    const auto *smh_896 = buffer.data(smh + 896);
    const auto *smh_897 = buffer.data(smh + 897);
    const auto *smh_898 = buffer.data(smh + 898);
    const auto *smh_899 = buffer.data(smh + 899);
    const auto *smh_900 = buffer.data(smh + 900);
    const auto *smh_901 = buffer.data(smh + 901);
    const auto *smh_902 = buffer.data(smh + 902);
    const auto *smh_918 = buffer.data(smh + 918);
    const auto *smh_919 = buffer.data(smh + 919);
    const auto *smh_920 = buffer.data(smh + 920);
    const auto *smh_921 = buffer.data(smh + 921);
    const auto *smh_922 = buffer.data(smh + 922);
    const auto *smh_923 = buffer.data(smh + 923);
    const auto *smh_924 = buffer.data(smh + 924);
    const auto *smh_927 = buffer.data(smh + 927);
    const auto *smh_929 = buffer.data(smh + 929);
    const auto *smh_930 = buffer.data(smh + 930);

    const auto *smi1_980 = buffer.data(smi1 + 980);
    const auto *smi1_983 = buffer.data(smi1 + 983);
    const auto *smi1_985 = buffer.data(smi1 + 985);
    const auto *smi1_986 = buffer.data(smi1 + 986);
    const auto *smi1_989 = buffer.data(smi1 + 989);
    const auto *smi1_990 = buffer.data(smi1 + 990);
    const auto *smi1_992 = buffer.data(smi1 + 992);
    const auto *smi1_994 = buffer.data(smi1 + 994);
    const auto *smi1_1007 = buffer.data(smi1 + 1007);

    const auto *sng0_610 = buffer.data(sng0 + 610);
    const auto *sng0_612 = buffer.data(sng0 + 612);
    const auto *sng0_613 = buffer.data(sng0 + 613);
    const auto *sng0_614 = buffer.data(sng0 + 614);
    const auto *sng0_615 = buffer.data(sng0 + 615);
    const auto *sng0_618 = buffer.data(sng0 + 618);
    const auto *sng0_620 = buffer.data(sng0 + 620);
    const auto *sng0_621 = buffer.data(sng0 + 621);
    const auto *sng0_624 = buffer.data(sng0 + 624);
    const auto *sng0_625 = buffer.data(sng0 + 625);
    const auto *sng0_627 = buffer.data(sng0 + 627);
    const auto *sng0_628 = buffer.data(sng0 + 628);
    const auto *sng0_629 = buffer.data(sng0 + 629);
    const auto *sng0_630 = buffer.data(sng0 + 630);
    const auto *sng0_633 = buffer.data(sng0 + 633);
    const auto *sng0_635 = buffer.data(sng0 + 635);
    const auto *sng0_636 = buffer.data(sng0 + 636);
    const auto *sng0_639 = buffer.data(sng0 + 639);
    const auto *sng0_640 = buffer.data(sng0 + 640);
    const auto *sng0_642 = buffer.data(sng0 + 642);
    const auto *sng0_643 = buffer.data(sng0 + 643);
    const auto *sng0_644 = buffer.data(sng0 + 644);
    const auto *sng0_655 = buffer.data(sng0 + 655);
    const auto *sng0_657 = buffer.data(sng0 + 657);
    const auto *sng0_658 = buffer.data(sng0 + 658);
    const auto *sng0_659 = buffer.data(sng0 + 659);
    const auto *sng0_660 = buffer.data(sng0 + 660);
    const auto *sng0_663 = buffer.data(sng0 + 663);
    const auto *sng0_665 = buffer.data(sng0 + 665);
    const auto *sng0_666 = buffer.data(sng0 + 666);

    const auto *sng1_610 = buffer.data(sng1 + 610);
    const auto *sng1_612 = buffer.data(sng1 + 612);
    const auto *sng1_613 = buffer.data(sng1 + 613);
    const auto *sng1_614 = buffer.data(sng1 + 614);
    const auto *sng1_615 = buffer.data(sng1 + 615);
    const auto *sng1_618 = buffer.data(sng1 + 618);
    const auto *sng1_620 = buffer.data(sng1 + 620);
    const auto *sng1_621 = buffer.data(sng1 + 621);
    const auto *sng1_624 = buffer.data(sng1 + 624);
    const auto *sng1_625 = buffer.data(sng1 + 625);
    const auto *sng1_627 = buffer.data(sng1 + 627);
    const auto *sng1_628 = buffer.data(sng1 + 628);
    const auto *sng1_629 = buffer.data(sng1 + 629);
    const auto *sng1_630 = buffer.data(sng1 + 630);
    const auto *sng1_633 = buffer.data(sng1 + 633);
    const auto *sng1_635 = buffer.data(sng1 + 635);
    const auto *sng1_636 = buffer.data(sng1 + 636);
    const auto *sng1_639 = buffer.data(sng1 + 639);
    const auto *sng1_640 = buffer.data(sng1 + 640);
    const auto *sng1_642 = buffer.data(sng1 + 642);
    const auto *sng1_643 = buffer.data(sng1 + 643);
    const auto *sng1_644 = buffer.data(sng1 + 644);
    const auto *sng1_655 = buffer.data(sng1 + 655);
    const auto *sng1_657 = buffer.data(sng1 + 657);
    const auto *sng1_658 = buffer.data(sng1 + 658);
    const auto *sng1_659 = buffer.data(sng1 + 659);
    const auto *sng1_660 = buffer.data(sng1 + 660);
    const auto *sng1_663 = buffer.data(sng1 + 663);
    const auto *sng1_665 = buffer.data(sng1 + 665);
    const auto *sng1_666 = buffer.data(sng1 + 666);

    const auto *snh_849 = buffer.data(snh + 849);
    const auto *snh_852 = buffer.data(snh + 852);
    const auto *snh_854 = buffer.data(snh + 854);
    const auto *snh_855 = buffer.data(snh + 855);
    const auto *snh_856 = buffer.data(snh + 856);
    const auto *snh_857 = buffer.data(snh + 857);
    const auto *snh_858 = buffer.data(snh + 858);
    const auto *snh_859 = buffer.data(snh + 859);
    const auto *snh_860 = buffer.data(snh + 860);
    const auto *snh_861 = buffer.data(snh + 861);
    const auto *snh_863 = buffer.data(snh + 863);
    const auto *snh_864 = buffer.data(snh + 864);
    const auto *snh_866 = buffer.data(snh + 866);
    const auto *snh_867 = buffer.data(snh + 867);
    const auto *snh_870 = buffer.data(snh + 870);
    const auto *snh_871 = buffer.data(snh + 871);
    const auto *snh_873 = buffer.data(snh + 873);
    const auto *snh_875 = buffer.data(snh + 875);
    const auto *snh_876 = buffer.data(snh + 876);
    const auto *snh_877 = buffer.data(snh + 877);
    const auto *snh_878 = buffer.data(snh + 878);
    const auto *snh_879 = buffer.data(snh + 879);
    const auto *snh_880 = buffer.data(snh + 880);
    const auto *snh_881 = buffer.data(snh + 881);
    const auto *snh_882 = buffer.data(snh + 882);
    const auto *snh_884 = buffer.data(snh + 884);
    const auto *snh_885 = buffer.data(snh + 885);
    const auto *snh_887 = buffer.data(snh + 887);
    const auto *snh_888 = buffer.data(snh + 888);
    const auto *snh_891 = buffer.data(snh + 891);
    const auto *snh_892 = buffer.data(snh + 892);
    const auto *snh_894 = buffer.data(snh + 894);
    const auto *snh_896 = buffer.data(snh + 896);
    const auto *snh_897 = buffer.data(snh + 897);
    const auto *snh_898 = buffer.data(snh + 898);
    const auto *snh_899 = buffer.data(snh + 899);
    const auto *snh_900 = buffer.data(snh + 900);
    const auto *snh_901 = buffer.data(snh + 901);
    const auto *snh_902 = buffer.data(snh + 902);
    const auto *snh_903 = buffer.data(snh + 903);
    const auto *snh_905 = buffer.data(snh + 905);
    const auto *snh_906 = buffer.data(snh + 906);
    const auto *snh_908 = buffer.data(snh + 908);
    const auto *snh_909 = buffer.data(snh + 909);
    const auto *snh_912 = buffer.data(snh + 912);
    const auto *snh_918 = buffer.data(snh + 918);
    const auto *snh_919 = buffer.data(snh + 919);
    const auto *snh_920 = buffer.data(snh + 920);
    const auto *snh_921 = buffer.data(snh + 921);
    const auto *snh_922 = buffer.data(snh + 922);
    const auto *snh_923 = buffer.data(snh + 923);
    const auto *snh_924 = buffer.data(snh + 924);
    const auto *snh_926 = buffer.data(snh + 926);
    const auto *snh_927 = buffer.data(snh + 927);
    const auto *snh_929 = buffer.data(snh + 929);
    const auto *snh_930 = buffer.data(snh + 930);

#pragma omp simd aligned(t_1132, t_1133, t_1134, pc_x, pc_y, smh_681, smh_852, smh_854, \
                         sng0_612, sng0_614, sng1_612, sng1_614, snh_849, snh_852, \
                         snh_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_12 * smh_852[k]
                    + f_8 * sng0_612[k]
                    - f_9 * sng1_612[k]
                    + f_3 * pc_x[k] * snh_852[k];

        t_1133[k] = f_14 * smh_681[k]
                    + f_3 * pc_y[k] * snh_849[k];

        t_1134[k] = f_12 * smh_854[k]
                    + f_8 * sng0_614[k]
                    - f_9 * sng1_614[k]
                    + f_3 * pc_x[k] * snh_854[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, t_1138, t_1139, pc_x, smh_855, smh_856, \
                         smh_857, smh_858, smh_859, snh_855, snh_856, snh_857, snh_858, \
                         snh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = f_12 * smh_855[k]
                    + f_3 * pc_x[k] * snh_855[k];

        t_1136[k] = f_12 * smh_856[k]
                    + f_3 * pc_x[k] * snh_856[k];

        t_1137[k] = f_12 * smh_857[k]
                    + f_3 * pc_x[k] * snh_857[k];

        t_1138[k] = f_12 * smh_858[k]
                    + f_3 * pc_x[k] * snh_858[k];

        t_1139[k] = f_12 * smh_859[k]
                    + f_3 * pc_x[k] * snh_859[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, pc_x, pc_y, pc_z, smh_666, smh_687, smh_860, \
                         sng0_610, sng1_610, snh_855, snh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_12 * smh_860[k]
                    + f_3 * pc_x[k] * snh_860[k];

        t_1141[k] = f_14 * smh_687[k]
                    + f_1 * sng0_610[k]
                    - f_2 * sng1_610[k]
                    + f_3 * pc_y[k] * snh_855[k];

        t_1142[k] = f_14 * smh_666[k]
                    + f_3 * pc_z[k] * snh_855[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, pc_y, smh_689, smh_690, smh_691, sng0_612, \
                         sng0_613, sng0_614, sng1_612, sng1_613, sng1_614, snh_857, snh_858, \
                         snh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_14 * smh_689[k]
                    + f_4 * sng0_612[k]
                    - f_5 * sng1_612[k]
                    + f_3 * pc_y[k] * snh_857[k];

        t_1144[k] = f_14 * smh_690[k]
                    + f_6 * sng0_613[k]
                    - f_7 * sng1_613[k]
                    + f_3 * pc_y[k] * snh_858[k];

        t_1145[k] = f_14 * smh_691[k]
                    + f_8 * sng0_614[k]
                    - f_9 * sng1_614[k]
                    + f_3 * pc_y[k] * snh_859[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pc_x, pc_y, pc_z, smh_671, smh_692, smh_861, \
                         sng0_614, sng0_615, sng1_614, sng1_615, snh_860, \
                         snh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_14 * smh_692[k]
                    + f_3 * pc_y[k] * snh_860[k];

        t_1147[k] = f_14 * smh_671[k]
                    + f_1 * sng0_614[k]
                    - f_2 * sng1_614[k]
                    + f_3 * pc_z[k] * snh_860[k];

        t_1148[k] = f_12 * smh_861[k]
                    + f_1 * sng0_615[k]
                    - f_2 * sng1_615[k]
                    + f_3 * pc_x[k] * snh_861[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, smh_672, smh_693, \
                         smh_695, smh_864, sng0_618, sng1_618, snh_861, snh_863, \
                         snh_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_13 * smh_693[k]
                    + f_3 * pc_y[k] * snh_861[k];

        t_1150[k] = f_19 * smh_672[k]
                    + f_3 * pc_z[k] * snh_861[k];

        t_1151[k] = f_12 * smh_864[k]
                    + f_4 * sng0_618[k]
                    - f_5 * sng1_618[k]
                    + f_3 * pc_x[k] * snh_864[k];

        t_1152[k] = f_13 * smh_695[k]
                    + f_3 * pc_y[k] * snh_863[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, pc_x, pc_z, smh_675, smh_866, smh_867, \
                         sng0_620, sng0_621, sng1_620, sng1_621, snh_864, snh_866, \
                         snh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_12 * smh_866[k]
                    + f_4 * sng0_620[k]
                    - f_5 * sng1_620[k]
                    + f_3 * pc_x[k] * snh_866[k];

        t_1154[k] = f_12 * smh_867[k]
                    + f_6 * sng0_621[k]
                    - f_7 * sng1_621[k]
                    + f_3 * pc_x[k] * snh_867[k];

        t_1155[k] = f_19 * smh_675[k]
                    + f_3 * pc_z[k] * snh_864[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, pc_x, pc_y, smh_698, smh_870, smh_871, \
                         sng0_624, sng0_625, sng1_624, sng1_625, snh_866, snh_870, \
                         snh_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_13 * smh_698[k]
                    + f_3 * pc_y[k] * snh_866[k];

        t_1157[k] = f_12 * smh_870[k]
                    + f_6 * sng0_624[k]
                    - f_7 * sng1_624[k]
                    + f_3 * pc_x[k] * snh_870[k];

        t_1158[k] = f_12 * smh_871[k]
                    + f_8 * sng0_625[k]
                    - f_9 * sng1_625[k]
                    + f_3 * pc_x[k] * snh_871[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, pc_x, pc_y, pc_z, smh_678, smh_702, smh_873, \
                         sng0_627, sng1_627, snh_867, snh_870, \
                         snh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_19 * smh_678[k]
                    + f_3 * pc_z[k] * snh_867[k];

        t_1160[k] = f_12 * smh_873[k]
                    + f_8 * sng0_627[k]
                    - f_9 * sng1_627[k]
                    + f_3 * pc_x[k] * snh_873[k];

        t_1161[k] = f_13 * smh_702[k]
                    + f_3 * pc_y[k] * snh_870[k];
    }

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, pc_x, smh_875, smh_876, smh_877, \
                         smh_878, sng0_629, sng1_629, snh_875, snh_876, snh_877, \
                         snh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = f_12 * smh_875[k]
                    + f_8 * sng0_629[k]
                    - f_9 * sng1_629[k]
                    + f_3 * pc_x[k] * snh_875[k];

        t_1163[k] = f_12 * smh_876[k]
                    + f_3 * pc_x[k] * snh_876[k];

        t_1164[k] = f_12 * smh_877[k]
                    + f_3 * pc_x[k] * snh_877[k];

        t_1165[k] = f_12 * smh_878[k]
                    + f_3 * pc_x[k] * snh_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, pc_x, pc_y, smh_708, smh_879, \
                         smh_880, smh_881, sng0_625, sng1_625, snh_876, snh_879, snh_880, \
                         snh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_12 * smh_879[k]
                    + f_3 * pc_x[k] * snh_879[k];

        t_1167[k] = f_12 * smh_880[k]
                    + f_3 * pc_x[k] * snh_880[k];

        t_1168[k] = f_12 * smh_881[k]
                    + f_3 * pc_x[k] * snh_881[k];

        t_1169[k] = f_13 * smh_708[k]
                    + f_1 * sng0_625[k]
                    - f_2 * sng1_625[k]
                    + f_3 * pc_y[k] * snh_876[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, pc_y, pc_z, smh_687, smh_710, smh_711, \
                         sng0_627, sng0_628, sng1_627, sng1_628, snh_876, snh_878, \
                         snh_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_19 * smh_687[k]
                    + f_3 * pc_z[k] * snh_876[k];

        t_1171[k] = f_13 * smh_710[k]
                    + f_4 * sng0_627[k]
                    - f_5 * sng1_627[k]
                    + f_3 * pc_y[k] * snh_878[k];

        t_1172[k] = f_13 * smh_711[k]
                    + f_6 * sng0_628[k]
                    - f_7 * sng1_628[k]
                    + f_3 * pc_y[k] * snh_879[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, pc_y, pc_z, smh_692, smh_712, smh_713, \
                         sng0_629, sng1_629, snh_880, snh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_13 * smh_712[k]
                    + f_8 * sng0_629[k]
                    - f_9 * sng1_629[k]
                    + f_3 * pc_y[k] * snh_880[k];

        t_1174[k] = f_13 * smh_713[k]
                    + f_3 * pc_y[k] * snh_881[k];

        t_1175[k] = f_19 * smh_692[k]
                    + f_1 * sng0_629[k]
                    - f_2 * sng1_629[k]
                    + f_3 * pc_z[k] * snh_881[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, pc_x, pc_y, pc_z, smh_693, smh_714, smh_882, \
                         sng0_630, sng1_630, snh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_12 * smh_882[k]
                    + f_1 * sng0_630[k]
                    - f_2 * sng1_630[k]
                    + f_3 * pc_x[k] * snh_882[k];

        t_1177[k] = f_12 * smh_714[k]
                    + f_3 * pc_y[k] * snh_882[k];

        t_1178[k] = f_18 * smh_693[k]
                    + f_3 * pc_z[k] * snh_882[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pc_x, pc_y, smh_716, smh_885, smh_887, \
                         sng0_633, sng0_635, sng1_633, sng1_635, snh_884, snh_885, \
                         snh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = f_12 * smh_885[k]
                    + f_4 * sng0_633[k]
                    - f_5 * sng1_633[k]
                    + f_3 * pc_x[k] * snh_885[k];

        t_1180[k] = f_12 * smh_716[k]
                    + f_3 * pc_y[k] * snh_884[k];

        t_1181[k] = f_12 * smh_887[k]
                    + f_4 * sng0_635[k]
                    - f_5 * sng1_635[k]
                    + f_3 * pc_x[k] * snh_887[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pc_x, pc_y, pc_z, smh_696, smh_719, smh_888, \
                         sng0_636, sng1_636, snh_885, snh_887, \
                         snh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_12 * smh_888[k]
                    + f_6 * sng0_636[k]
                    - f_7 * sng1_636[k]
                    + f_3 * pc_x[k] * snh_888[k];

        t_1183[k] = f_18 * smh_696[k]
                    + f_3 * pc_z[k] * snh_885[k];

        t_1184[k] = f_12 * smh_719[k]
                    + f_3 * pc_y[k] * snh_887[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_x, pc_z, smh_699, smh_891, smh_892, \
                         sng0_639, sng0_640, sng1_639, sng1_640, snh_888, snh_891, \
                         snh_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_12 * smh_891[k]
                    + f_6 * sng0_639[k]
                    - f_7 * sng1_639[k]
                    + f_3 * pc_x[k] * snh_891[k];

        t_1186[k] = f_12 * smh_892[k]
                    + f_8 * sng0_640[k]
                    - f_9 * sng1_640[k]
                    + f_3 * pc_x[k] * snh_892[k];

        t_1187[k] = f_18 * smh_699[k]
                    + f_3 * pc_z[k] * snh_888[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, pc_x, pc_y, smh_723, smh_894, smh_896, \
                         sng0_642, sng0_644, sng1_642, sng1_644, snh_891, snh_894, \
                         snh_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_12 * smh_894[k]
                    + f_8 * sng0_642[k]
                    - f_9 * sng1_642[k]
                    + f_3 * pc_x[k] * snh_894[k];

        t_1189[k] = f_12 * smh_723[k]
                    + f_3 * pc_y[k] * snh_891[k];

        t_1190[k] = f_12 * smh_896[k]
                    + f_8 * sng0_644[k]
                    - f_9 * sng1_644[k]
                    + f_3 * pc_x[k] * snh_896[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, t_1194, t_1195, pc_x, smh_897, smh_898, \
                         smh_899, smh_900, smh_901, snh_897, snh_898, snh_899, snh_900, \
                         snh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = f_12 * smh_897[k]
                    + f_3 * pc_x[k] * snh_897[k];

        t_1192[k] = f_12 * smh_898[k]
                    + f_3 * pc_x[k] * snh_898[k];

        t_1193[k] = f_12 * smh_899[k]
                    + f_3 * pc_x[k] * snh_899[k];

        t_1194[k] = f_12 * smh_900[k]
                    + f_3 * pc_x[k] * snh_900[k];

        t_1195[k] = f_12 * smh_901[k]
                    + f_3 * pc_x[k] * snh_901[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, pc_z, smh_708, smh_729, smh_902, \
                         sng0_640, sng1_640, snh_897, snh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_12 * smh_902[k]
                    + f_3 * pc_x[k] * snh_902[k];

        t_1197[k] = f_12 * smh_729[k]
                    + f_1 * sng0_640[k]
                    - f_2 * sng1_640[k]
                    + f_3 * pc_y[k] * snh_897[k];

        t_1198[k] = f_18 * smh_708[k]
                    + f_3 * pc_z[k] * snh_897[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_y, smh_731, smh_732, smh_733, sng0_642, \
                         sng0_643, sng0_644, sng1_642, sng1_643, sng1_644, snh_899, snh_900, \
                         snh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_12 * smh_731[k]
                    + f_4 * sng0_642[k]
                    - f_5 * sng1_642[k]
                    + f_3 * pc_y[k] * snh_899[k];

        t_1200[k] = f_12 * smh_732[k]
                    + f_6 * sng0_643[k]
                    - f_7 * sng1_643[k]
                    + f_3 * pc_y[k] * snh_900[k];

        t_1201[k] = f_12 * smh_733[k]
                    + f_8 * sng0_644[k]
                    - f_9 * sng1_644[k]
                    + f_3 * pc_y[k] * snh_901[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, pb_y, pc_y, pc_z, smi0_980, smh_713, \
                         smh_734, smh_735, smi1_980, sng0_644, sng1_644, snh_902, \
                         snh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_12 * smh_734[k]
                    + f_3 * pc_y[k] * snh_902[k];

        t_1203[k] = f_18 * smh_713[k]
                    + f_1 * sng0_644[k]
                    - f_2 * sng1_644[k]
                    + f_3 * pc_z[k] * snh_902[k];

        t_1204[k] = pb_y[k] * smi0_980[k]
                    - f_10 * pc_y[k] * smi1_980[k];

        t_1205[k] = f_11 * smh_735[k]
                    + f_3 * pc_y[k] * snh_903[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pb_y, pc_y, pc_z, smi0_983, smi0_985, \
                         smh_714, smh_736, smh_737, smi1_983, smi1_985, snh_903, \
                         snh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_17 * smh_714[k]
                    + f_3 * pc_z[k] * snh_903[k];

        t_1207[k] = pb_y[k] * smi0_983[k]
                    + f_12 * smh_736[k]
                    - f_10 * pc_y[k] * smi1_983[k];

        t_1208[k] = f_11 * smh_737[k]
                    + f_3 * pc_y[k] * snh_905[k];

        t_1209[k] = pb_y[k] * smi0_985[k]
                    - f_10 * pc_y[k] * smi1_985[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pb_y, pc_y, pc_z, smi0_986, smi0_989, \
                         smh_717, smh_738, smh_740, smi1_986, smi1_989, snh_906, \
                         snh_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = pb_y[k] * smi0_986[k]
                    + f_13 * smh_738[k]
                    - f_10 * pc_y[k] * smi1_986[k];

        t_1211[k] = f_17 * smh_717[k]
                    + f_3 * pc_z[k] * snh_906[k];

        t_1212[k] = f_11 * smh_740[k]
                    + f_3 * pc_y[k] * snh_908[k];

        t_1213[k] = pb_y[k] * smi0_989[k]
                    - f_10 * pc_y[k] * smi1_989[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, pb_y, pc_y, pc_z, smi0_990, smi0_992, \
                         smh_720, smh_741, smh_743, smi1_990, smi1_992, \
                         snh_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pb_y[k] * smi0_990[k]
                    + f_14 * smh_741[k]
                    - f_10 * pc_y[k] * smi1_990[k];

        t_1215[k] = f_17 * smh_720[k]
                    + f_3 * pc_z[k] * snh_909[k];

        t_1216[k] = pb_y[k] * smi0_992[k]
                    + f_12 * smh_743[k]
                    - f_10 * pc_y[k] * smi1_992[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, pb_y, pc_x, pc_y, smi0_994, smh_744, \
                         smh_918, smh_919, smi1_994, snh_912, snh_918, \
                         snh_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_11 * smh_744[k]
                    + f_3 * pc_y[k] * snh_912[k];

        t_1218[k] = pb_y[k] * smi0_994[k]
                    - f_10 * pc_y[k] * smi1_994[k];

        t_1219[k] = f_12 * smh_918[k]
                    + f_3 * pc_x[k] * snh_918[k];

        t_1220[k] = f_12 * smh_919[k]
                    + f_3 * pc_x[k] * snh_919[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, t_1224, pc_x, smh_920, smh_921, smh_922, \
                         smh_923, snh_920, snh_921, snh_922, snh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_12 * smh_920[k]
                    + f_3 * pc_x[k] * snh_920[k];

        t_1222[k] = f_12 * smh_921[k]
                    + f_3 * pc_x[k] * snh_921[k];

        t_1223[k] = f_12 * smh_922[k]
                    + f_3 * pc_x[k] * snh_922[k];

        t_1224[k] = f_12 * smh_923[k]
                    + f_3 * pc_x[k] * snh_923[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, pc_y, pc_z, smh_729, smh_750, smh_752, \
                         sng0_655, sng0_657, sng1_655, sng1_657, snh_918, \
                         snh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_11 * smh_750[k]
                    + f_1 * sng0_655[k]
                    - f_2 * sng1_655[k]
                    + f_3 * pc_y[k] * snh_918[k];

        t_1226[k] = f_17 * smh_729[k]
                    + f_3 * pc_z[k] * snh_918[k];

        t_1227[k] = f_11 * smh_752[k]
                    + f_4 * sng0_657[k]
                    - f_5 * sng1_657[k]
                    + f_3 * pc_y[k] * snh_920[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, pc_y, smh_753, smh_754, smh_755, sng0_658, \
                         sng0_659, sng1_658, sng1_659, snh_921, snh_922, \
                         snh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = f_11 * smh_753[k]
                    + f_6 * sng0_658[k]
                    - f_7 * sng1_658[k]
                    + f_3 * pc_y[k] * snh_921[k];

        t_1229[k] = f_11 * smh_754[k]
                    + f_8 * sng0_659[k]
                    - f_9 * sng1_659[k]
                    + f_3 * pc_y[k] * snh_922[k];

        t_1230[k] = f_11 * smh_755[k]
                    + f_3 * pc_y[k] * snh_923[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pb_y, pc_x, pc_y, pc_z, smi0_1007, \
                         smh_735, smh_924, smi1_1007, sng0_660, sng1_660, \
                         snh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = pb_y[k] * smi0_1007[k]
                    - f_10 * pc_y[k] * smi1_1007[k];

        t_1232[k] = f_12 * smh_924[k]
                    + f_1 * sng0_660[k]
                    - f_2 * sng1_660[k]
                    + f_3 * pc_x[k] * snh_924[k];

        t_1233[k] = f_3 * pc_y[k] * snh_924[k];

        t_1234[k] = f_16 * smh_735[k]
                    + f_3 * pc_z[k] * snh_924[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pc_x, pc_y, smh_927, smh_929, sng0_663, \
                         sng0_665, sng1_663, sng1_665, snh_926, snh_927, \
                         snh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_12 * smh_927[k]
                    + f_4 * sng0_663[k]
                    - f_5 * sng1_663[k]
                    + f_3 * pc_x[k] * snh_927[k];

        t_1236[k] = f_3 * pc_y[k] * snh_926[k];

        t_1237[k] = f_12 * smh_929[k]
                    + f_4 * sng0_665[k]
                    - f_5 * sng1_665[k]
                    + f_3 * pc_x[k] * snh_929[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, pc_x, pc_y, pc_z, smh_738, smh_930, sng0_666, \
                         sng1_666, snh_927, snh_929, snh_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_12 * smh_930[k]
                    + f_6 * sng0_666[k]
                    - f_7 * sng1_666[k]
                    + f_3 * pc_x[k] * snh_930[k];

        t_1239[k] = f_16 * smh_738[k]
                    + f_3 * pc_z[k] * snh_927[k];

        t_1240[k] = f_3 * pc_y[k] * snh_929[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smi0,
                                                           const size_t smh, const size_t smi1,
                                                           const size_t sng0, const size_t sng1,
                                                           const size_t snh, const size_t ncols,
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
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi0_1008 = buffer.data(smi0 + 1008);
    const auto *smi0_1011 = buffer.data(smi0 + 1011);
    const auto *smi0_1014 = buffer.data(smi0 + 1014);
    const auto *smi0_1018 = buffer.data(smi0 + 1018);
    const auto *smi0_1260 = buffer.data(smi0 + 1260);
    const auto *smi0_1263 = buffer.data(smi0 + 1263);
    const auto *smi0_1265 = buffer.data(smi0 + 1265);
    const auto *smi0_1266 = buffer.data(smi0 + 1266);
    const auto *smi0_1269 = buffer.data(smi0 + 1269);
    const auto *smi0_1270 = buffer.data(smi0 + 1270);
    const auto *smi0_1272 = buffer.data(smi0 + 1272);
    const auto *smi0_1274 = buffer.data(smi0 + 1274);
    const auto *smi0_1281 = buffer.data(smi0 + 1281);
    const auto *smi0_1283 = buffer.data(smi0 + 1283);
    const auto *smi0_1284 = buffer.data(smi0 + 1284);
    const auto *smi0_1285 = buffer.data(smi0 + 1285);
    const auto *smi0_1287 = buffer.data(smi0 + 1287);
    const auto *smi0_1293 = buffer.data(smi0 + 1293);
    const auto *smi0_1297 = buffer.data(smi0 + 1297);
    const auto *smi0_1300 = buffer.data(smi0 + 1300);
    const auto *smi0_1302 = buffer.data(smi0 + 1302);
    const auto *smi0_1309 = buffer.data(smi0 + 1309);
    const auto *smi0_1311 = buffer.data(smi0 + 1311);
    const auto *smi0_1312 = buffer.data(smi0 + 1312);
    const auto *smi0_1313 = buffer.data(smi0 + 1313);
    const auto *smi0_1315 = buffer.data(smi0 + 1315);
    const auto *smi0_1316 = buffer.data(smi0 + 1316);
    const auto *smi0_1319 = buffer.data(smi0 + 1319);
    const auto *smi0_1321 = buffer.data(smi0 + 1321);
    const auto *smi0_1322 = buffer.data(smi0 + 1322);
    const auto *smi0_1325 = buffer.data(smi0 + 1325);
    const auto *smi0_1326 = buffer.data(smi0 + 1326);
    const auto *smi0_1328 = buffer.data(smi0 + 1328);
    const auto *smi0_1330 = buffer.data(smi0 + 1330);
    const auto *smi0_1337 = buffer.data(smi0 + 1337);
    const auto *smi0_1339 = buffer.data(smi0 + 1339);
    const auto *smi0_1340 = buffer.data(smi0 + 1340);
    const auto *smi0_1341 = buffer.data(smi0 + 1341);
    const auto *smi0_1343 = buffer.data(smi0 + 1343);
    const auto *smi0_1344 = buffer.data(smi0 + 1344);
    const auto *smi0_1347 = buffer.data(smi0 + 1347);
    const auto *smi0_1349 = buffer.data(smi0 + 1349);
    const auto *smi0_1350 = buffer.data(smi0 + 1350);
    const auto *smi0_1353 = buffer.data(smi0 + 1353);
    const auto *smi0_1354 = buffer.data(smi0 + 1354);
    const auto *smi0_1356 = buffer.data(smi0 + 1356);
    const auto *smi0_1358 = buffer.data(smi0 + 1358);

    const auto *smh_741 = buffer.data(smh + 741);
    const auto *smh_750 = buffer.data(smh + 750);
    const auto *smh_755 = buffer.data(smh + 755);
    const auto *smh_756 = buffer.data(smh + 756);
    const auto *smh_758 = buffer.data(smh + 758);
    const auto *smh_759 = buffer.data(smh + 759);
    const auto *smh_761 = buffer.data(smh + 761);
    const auto *smh_762 = buffer.data(smh + 762);
    const auto *smh_765 = buffer.data(smh + 765);
    const auto *smh_771 = buffer.data(smh + 771);
    const auto *smh_776 = buffer.data(smh + 776);
    const auto *smh_777 = buffer.data(smh + 777);
    const auto *smh_779 = buffer.data(smh + 779);
    const auto *smh_780 = buffer.data(smh + 780);
    const auto *smh_782 = buffer.data(smh + 782);
    const auto *smh_783 = buffer.data(smh + 783);
    const auto *smh_786 = buffer.data(smh + 786);
    const auto *smh_792 = buffer.data(smh + 792);
    const auto *smh_797 = buffer.data(smh + 797);
    const auto *smh_798 = buffer.data(smh + 798);
    const auto *smh_800 = buffer.data(smh + 800);
    const auto *smh_801 = buffer.data(smh + 801);
    const auto *smh_803 = buffer.data(smh + 803);
    const auto *smh_804 = buffer.data(smh + 804);
    const auto *smh_807 = buffer.data(smh + 807);
    const auto *smh_818 = buffer.data(smh + 818);
    const auto *smh_819 = buffer.data(smh + 819);
    const auto *smh_821 = buffer.data(smh + 821);
    const auto *smh_824 = buffer.data(smh + 824);
    const auto *smh_828 = buffer.data(smh + 828);
    const auto *smh_933 = buffer.data(smh + 933);
    const auto *smh_934 = buffer.data(smh + 934);
    const auto *smh_936 = buffer.data(smh + 936);
    const auto *smh_938 = buffer.data(smh + 938);
    const auto *smh_939 = buffer.data(smh + 939);
    const auto *smh_940 = buffer.data(smh + 940);
    const auto *smh_941 = buffer.data(smh + 941);
    const auto *smh_942 = buffer.data(smh + 942);
    const auto *smh_943 = buffer.data(smh + 943);
    const auto *smh_944 = buffer.data(smh + 944);
    const auto *smh_945 = buffer.data(smh + 945);
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
    const auto *smh_971 = buffer.data(smh + 971);
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

    const auto *smi1_1008 = buffer.data(smi1 + 1008);
    const auto *smi1_1011 = buffer.data(smi1 + 1011);
    const auto *smi1_1014 = buffer.data(smi1 + 1014);
    const auto *smi1_1018 = buffer.data(smi1 + 1018);
    const auto *smi1_1260 = buffer.data(smi1 + 1260);
    const auto *smi1_1263 = buffer.data(smi1 + 1263);
    const auto *smi1_1265 = buffer.data(smi1 + 1265);
    const auto *smi1_1266 = buffer.data(smi1 + 1266);
    const auto *smi1_1269 = buffer.data(smi1 + 1269);
    const auto *smi1_1270 = buffer.data(smi1 + 1270);
    const auto *smi1_1272 = buffer.data(smi1 + 1272);
    const auto *smi1_1274 = buffer.data(smi1 + 1274);
    const auto *smi1_1281 = buffer.data(smi1 + 1281);
    const auto *smi1_1283 = buffer.data(smi1 + 1283);
    const auto *smi1_1284 = buffer.data(smi1 + 1284);
    const auto *smi1_1285 = buffer.data(smi1 + 1285);
    const auto *smi1_1287 = buffer.data(smi1 + 1287);
    const auto *smi1_1293 = buffer.data(smi1 + 1293);
    const auto *smi1_1297 = buffer.data(smi1 + 1297);
    const auto *smi1_1300 = buffer.data(smi1 + 1300);
    const auto *smi1_1302 = buffer.data(smi1 + 1302);
    const auto *smi1_1309 = buffer.data(smi1 + 1309);
    const auto *smi1_1311 = buffer.data(smi1 + 1311);
    const auto *smi1_1312 = buffer.data(smi1 + 1312);
    const auto *smi1_1313 = buffer.data(smi1 + 1313);
    const auto *smi1_1315 = buffer.data(smi1 + 1315);
    const auto *smi1_1316 = buffer.data(smi1 + 1316);
    const auto *smi1_1319 = buffer.data(smi1 + 1319);
    const auto *smi1_1321 = buffer.data(smi1 + 1321);
    const auto *smi1_1322 = buffer.data(smi1 + 1322);
    const auto *smi1_1325 = buffer.data(smi1 + 1325);
    const auto *smi1_1326 = buffer.data(smi1 + 1326);
    const auto *smi1_1328 = buffer.data(smi1 + 1328);
    const auto *smi1_1330 = buffer.data(smi1 + 1330);
    const auto *smi1_1337 = buffer.data(smi1 + 1337);
    const auto *smi1_1339 = buffer.data(smi1 + 1339);
    const auto *smi1_1340 = buffer.data(smi1 + 1340);
    const auto *smi1_1341 = buffer.data(smi1 + 1341);
    const auto *smi1_1343 = buffer.data(smi1 + 1343);
    const auto *smi1_1344 = buffer.data(smi1 + 1344);
    const auto *smi1_1347 = buffer.data(smi1 + 1347);
    const auto *smi1_1349 = buffer.data(smi1 + 1349);
    const auto *smi1_1350 = buffer.data(smi1 + 1350);
    const auto *smi1_1353 = buffer.data(smi1 + 1353);
    const auto *smi1_1354 = buffer.data(smi1 + 1354);
    const auto *smi1_1356 = buffer.data(smi1 + 1356);
    const auto *smi1_1358 = buffer.data(smi1 + 1358);

    const auto *sng0_669 = buffer.data(sng0 + 669);
    const auto *sng0_670 = buffer.data(sng0 + 670);
    const auto *sng0_672 = buffer.data(sng0 + 672);
    const auto *sng0_673 = buffer.data(sng0 + 673);
    const auto *sng0_674 = buffer.data(sng0 + 674);

    const auto *sng1_669 = buffer.data(sng1 + 669);
    const auto *sng1_670 = buffer.data(sng1 + 670);
    const auto *sng1_672 = buffer.data(sng1 + 672);
    const auto *sng1_673 = buffer.data(sng1 + 673);
    const auto *sng1_674 = buffer.data(sng1 + 674);

    const auto *snh_930 = buffer.data(snh + 930);
    const auto *snh_933 = buffer.data(snh + 933);
    const auto *snh_934 = buffer.data(snh + 934);
    const auto *snh_936 = buffer.data(snh + 936);
    const auto *snh_938 = buffer.data(snh + 938);
    const auto *snh_939 = buffer.data(snh + 939);
    const auto *snh_940 = buffer.data(snh + 940);
    const auto *snh_941 = buffer.data(snh + 941);
    const auto *snh_942 = buffer.data(snh + 942);
    const auto *snh_943 = buffer.data(snh + 943);
    const auto *snh_944 = buffer.data(snh + 944);
    const auto *snh_945 = buffer.data(snh + 945);
    const auto *snh_947 = buffer.data(snh + 947);
    const auto *snh_948 = buffer.data(snh + 948);
    const auto *snh_950 = buffer.data(snh + 950);
    const auto *snh_951 = buffer.data(snh + 951);
    const auto *snh_954 = buffer.data(snh + 954);
    const auto *snh_960 = buffer.data(snh + 960);
    const auto *snh_961 = buffer.data(snh + 961);
    const auto *snh_962 = buffer.data(snh + 962);
    const auto *snh_963 = buffer.data(snh + 963);
    const auto *snh_964 = buffer.data(snh + 964);
    const auto *snh_965 = buffer.data(snh + 965);
    const auto *snh_966 = buffer.data(snh + 966);
    const auto *snh_968 = buffer.data(snh + 968);
    const auto *snh_969 = buffer.data(snh + 969);
    const auto *snh_971 = buffer.data(snh + 971);
    const auto *snh_972 = buffer.data(snh + 972);
    const auto *snh_975 = buffer.data(snh + 975);
    const auto *snh_981 = buffer.data(snh + 981);
    const auto *snh_982 = buffer.data(snh + 982);
    const auto *snh_983 = buffer.data(snh + 983);
    const auto *snh_984 = buffer.data(snh + 984);
    const auto *snh_985 = buffer.data(snh + 985);
    const auto *snh_986 = buffer.data(snh + 986);
    const auto *snh_987 = buffer.data(snh + 987);
    const auto *snh_989 = buffer.data(snh + 989);
    const auto *snh_990 = buffer.data(snh + 990);
    const auto *snh_992 = buffer.data(snh + 992);
    const auto *snh_993 = buffer.data(snh + 993);
    const auto *snh_996 = buffer.data(snh + 996);
    const auto *snh_1002 = buffer.data(snh + 1002);
    const auto *snh_1003 = buffer.data(snh + 1003);
    const auto *snh_1004 = buffer.data(snh + 1004);
    const auto *snh_1005 = buffer.data(snh + 1005);
    const auto *snh_1006 = buffer.data(snh + 1006);
    const auto *snh_1007 = buffer.data(snh + 1007);
    const auto *snh_1008 = buffer.data(snh + 1008);
    const auto *snh_1010 = buffer.data(snh + 1010);
    const auto *snh_1011 = buffer.data(snh + 1011);
    const auto *snh_1013 = buffer.data(snh + 1013);
    const auto *snh_1014 = buffer.data(snh + 1014);
    const auto *snh_1017 = buffer.data(snh + 1017);
    const auto *snh_1023 = buffer.data(snh + 1023);
    const auto *snh_1024 = buffer.data(snh + 1024);
    const auto *snh_1025 = buffer.data(snh + 1025);

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pc_z, smh_741, smh_933, smh_934, \
                         sng0_669, sng0_670, sng1_669, sng1_670, snh_930, snh_933, \
                         snh_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_12 * smh_933[k]
                    + f_6 * sng0_669[k]
                    - f_7 * sng1_669[k]
                    + f_3 * pc_x[k] * snh_933[k];

        t_1242[k] = f_12 * smh_934[k]
                    + f_8 * sng0_670[k]
                    - f_9 * sng1_670[k]
                    + f_3 * pc_x[k] * snh_934[k];

        t_1243[k] = f_16 * smh_741[k]
                    + f_3 * pc_z[k] * snh_930[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, pc_x, pc_y, smh_936, smh_938, sng0_672, \
                         sng0_674, sng1_672, sng1_674, snh_933, snh_936, \
                         snh_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_12 * smh_936[k]
                    + f_8 * sng0_672[k]
                    - f_9 * sng1_672[k]
                    + f_3 * pc_x[k] * snh_936[k];

        t_1245[k] = f_3 * pc_y[k] * snh_933[k];

        t_1246[k] = f_12 * smh_938[k]
                    + f_8 * sng0_674[k]
                    - f_9 * sng1_674[k]
                    + f_3 * pc_x[k] * snh_938[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, t_1251, pc_x, smh_939, smh_940, \
                         smh_941, smh_942, smh_943, snh_939, snh_940, snh_941, snh_942, \
                         snh_943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = f_12 * smh_939[k]
                    + f_3 * pc_x[k] * snh_939[k];

        t_1248[k] = f_12 * smh_940[k]
                    + f_3 * pc_x[k] * snh_940[k];

        t_1249[k] = f_12 * smh_941[k]
                    + f_3 * pc_x[k] * snh_941[k];

        t_1250[k] = f_12 * smh_942[k]
                    + f_3 * pc_x[k] * snh_942[k];

        t_1251[k] = f_12 * smh_943[k]
                    + f_3 * pc_x[k] * snh_943[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, pc_x, pc_y, pc_z, smh_750, smh_944, \
                         sng0_670, sng0_672, sng1_670, sng1_672, snh_939, snh_941, \
                         snh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_12 * smh_944[k]
                    + f_3 * pc_x[k] * snh_944[k];

        t_1253[k] = f_1 * sng0_670[k]
                    - f_2 * sng1_670[k]
                    + f_3 * pc_y[k] * snh_939[k];

        t_1254[k] = f_16 * smh_750[k]
                    + f_3 * pc_z[k] * snh_939[k];

        t_1255[k] = f_4 * sng0_672[k]
                    - f_5 * sng1_672[k]
                    + f_3 * pc_y[k] * snh_941[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pc_y, pc_z, smh_755, sng0_673, \
                         sng0_674, sng1_673, sng1_674, snh_942, snh_943, \
                         snh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_6 * sng0_673[k]
                    - f_7 * sng1_673[k]
                    + f_3 * pc_y[k] * snh_942[k];

        t_1257[k] = f_8 * sng0_674[k]
                    - f_9 * sng1_674[k]
                    + f_3 * pc_y[k] * snh_943[k];

        t_1258[k] = f_3 * pc_y[k] * snh_944[k];

        t_1259[k] = f_16 * smh_755[k]
                    + f_1 * sng0_674[k]
                    - f_2 * sng1_674[k]
                    + f_3 * pc_z[k] * snh_944[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pb_x, pc_x, pc_y, pc_z, smi0_1260, \
                         smi0_1263, smh_756, smh_945, smh_948, smi1_1260, smi1_1263, \
                         snh_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = pb_x[k] * smi0_1260[k]
                    + f_18 * smh_945[k]
                    - f_10 * pc_x[k] * smi1_1260[k];

        t_1261[k] = f_15 * smh_756[k]
                    + f_3 * pc_y[k] * snh_945[k];

        t_1262[k] = f_3 * pc_z[k] * snh_945[k];

        t_1263[k] = pb_x[k] * smi0_1263[k]
                    + f_14 * smh_948[k]
                    - f_10 * pc_x[k] * smi1_1263[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, pb_x, pc_x, pc_y, smi0_1265, smi0_1266, \
                         smh_758, smh_950, smh_951, smi1_1265, smi1_1266, \
                         snh_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_15 * smh_758[k]
                    + f_3 * pc_y[k] * snh_947[k];

        t_1265[k] = pb_x[k] * smi0_1265[k]
                    + f_14 * smh_950[k]
                    - f_10 * pc_x[k] * smi1_1265[k];

        t_1266[k] = pb_x[k] * smi0_1266[k]
                    + f_13 * smh_951[k]
                    - f_10 * pc_x[k] * smi1_1266[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, pb_x, pc_x, pc_y, pc_z, smi0_1269, smh_761, \
                         smh_954, smi1_1269, snh_948, snh_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_3 * pc_z[k] * snh_948[k];

        t_1268[k] = f_15 * smh_761[k]
                    + f_3 * pc_y[k] * snh_950[k];

        t_1269[k] = pb_x[k] * smi0_1269[k]
                    + f_13 * smh_954[k]
                    - f_10 * pc_x[k] * smi1_1269[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, pb_x, pc_x, pc_z, smi0_1270, smi0_1272, \
                         smh_955, smh_957, smi1_1270, smi1_1272, \
                         snh_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = pb_x[k] * smi0_1270[k]
                    + f_12 * smh_955[k]
                    - f_10 * pc_x[k] * smi1_1270[k];

        t_1271[k] = f_3 * pc_z[k] * snh_951[k];

        t_1272[k] = pb_x[k] * smi0_1272[k]
                    + f_12 * smh_957[k]
                    - f_10 * pc_x[k] * smi1_1272[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, t_1276, pb_x, pc_x, pc_y, smi0_1274, smh_765, \
                         smh_959, smh_960, smh_961, smi1_1274, snh_954, snh_960, \
                         snh_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_15 * smh_765[k]
                    + f_3 * pc_y[k] * snh_954[k];

        t_1274[k] = pb_x[k] * smi0_1274[k]
                    + f_12 * smh_959[k]
                    - f_10 * pc_x[k] * smi1_1274[k];

        t_1275[k] = f_11 * smh_960[k]
                    + f_3 * pc_x[k] * snh_960[k];

        t_1276[k] = f_11 * smh_961[k]
                    + f_3 * pc_x[k] * snh_961[k];
    }

#pragma omp simd aligned(t_1277, t_1278, t_1279, t_1280, pc_x, smh_962, smh_963, smh_964, \
                         smh_965, snh_962, snh_963, snh_964, snh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = f_11 * smh_962[k]
                    + f_3 * pc_x[k] * snh_962[k];

        t_1278[k] = f_11 * smh_963[k]
                    + f_3 * pc_x[k] * snh_963[k];

        t_1279[k] = f_11 * smh_964[k]
                    + f_3 * pc_x[k] * snh_964[k];

        t_1280[k] = f_11 * smh_965[k]
                    + f_3 * pc_x[k] * snh_965[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, t_1284, pb_x, pc_x, pc_z, smi0_1281, \
                         smi0_1283, smi0_1284, smi1_1281, smi1_1283, smi1_1284, \
                         snh_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = pb_x[k] * smi0_1281[k]
                    - f_10 * pc_x[k] * smi1_1281[k];

        t_1282[k] = f_3 * pc_z[k] * snh_960[k];

        t_1283[k] = pb_x[k] * smi0_1283[k]
                    - f_10 * pc_x[k] * smi1_1283[k];

        t_1284[k] = pb_x[k] * smi0_1284[k]
                    - f_10 * pc_x[k] * smi1_1284[k];
    }

#pragma omp simd aligned(t_1285, t_1286, t_1287, pb_x, pc_x, pc_y, smi0_1285, smi0_1287, \
                         smh_776, smi1_1285, smi1_1287, snh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1285[k] = pb_x[k] * smi0_1285[k]
                    - f_10 * pc_x[k] * smi1_1285[k];

        t_1286[k] = f_15 * smh_776[k]
                    + f_3 * pc_y[k] * snh_965[k];

        t_1287[k] = pb_x[k] * smi0_1287[k]
                    - f_10 * pc_x[k] * smi1_1287[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pb_z, pc_y, pc_z, smi0_1008, \
                         smi0_1011, smh_756, smh_777, smi1_1008, smi1_1011, \
                         snh_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = pb_z[k] * smi0_1008[k]
                    - f_10 * pc_z[k] * smi1_1008[k];

        t_1289[k] = f_16 * smh_777[k]
                    + f_3 * pc_y[k] * snh_966[k];

        t_1290[k] = f_11 * smh_756[k]
                    + f_3 * pc_z[k] * snh_966[k];

        t_1291[k] = pb_z[k] * smi0_1011[k]
                    - f_10 * pc_z[k] * smi1_1011[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, pb_x, pb_z, pc_x, pc_y, pc_z, smi0_1014, \
                         smi0_1293, smh_779, smh_971, smi1_1014, smi1_1293, \
                         snh_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_16 * smh_779[k]
                    + f_3 * pc_y[k] * snh_968[k];

        t_1293[k] = pb_x[k] * smi0_1293[k]
                    + f_14 * smh_971[k]
                    - f_10 * pc_x[k] * smi1_1293[k];

        t_1294[k] = pb_z[k] * smi0_1014[k]
                    - f_10 * pc_z[k] * smi1_1014[k];
    }

#pragma omp simd aligned(t_1295, t_1296, t_1297, pb_x, pc_x, pc_y, pc_z, smi0_1297, smh_759, \
                         smh_782, smh_975, smi1_1297, snh_969, \
                         snh_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1295[k] = f_11 * smh_759[k]
                    + f_3 * pc_z[k] * snh_969[k];

        t_1296[k] = f_16 * smh_782[k]
                    + f_3 * pc_y[k] * snh_971[k];

        t_1297[k] = pb_x[k] * smi0_1297[k]
                    + f_13 * smh_975[k]
                    - f_10 * pc_x[k] * smi1_1297[k];
    }

#pragma omp simd aligned(t_1298, t_1299, t_1300, pb_x, pb_z, pc_x, pc_z, smi0_1018, smi0_1300, \
                         smh_762, smh_978, smi1_1018, smi1_1300, \
                         snh_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1298[k] = pb_z[k] * smi0_1018[k]
                    - f_10 * pc_z[k] * smi1_1018[k];

        t_1299[k] = f_11 * smh_762[k]
                    + f_3 * pc_z[k] * snh_972[k];

        t_1300[k] = pb_x[k] * smi0_1300[k]
                    + f_12 * smh_978[k]
                    - f_10 * pc_x[k] * smi1_1300[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, pb_x, pc_x, pc_y, smi0_1302, smh_786, \
                         smh_980, smh_981, smh_982, smi1_1302, snh_975, snh_981, \
                         snh_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_16 * smh_786[k]
                    + f_3 * pc_y[k] * snh_975[k];

        t_1302[k] = pb_x[k] * smi0_1302[k]
                    + f_12 * smh_980[k]
                    - f_10 * pc_x[k] * smi1_1302[k];

        t_1303[k] = f_11 * smh_981[k]
                    + f_3 * pc_x[k] * snh_981[k];

        t_1304[k] = f_11 * smh_982[k]
                    + f_3 * pc_x[k] * snh_982[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, t_1308, pc_x, smh_983, smh_984, smh_985, \
                         smh_986, snh_983, snh_984, snh_985, snh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = f_11 * smh_983[k]
                    + f_3 * pc_x[k] * snh_983[k];

        t_1306[k] = f_11 * smh_984[k]
                    + f_3 * pc_x[k] * snh_984[k];

        t_1307[k] = f_11 * smh_985[k]
                    + f_3 * pc_x[k] * snh_985[k];

        t_1308[k] = f_11 * smh_986[k]
                    + f_3 * pc_x[k] * snh_986[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, t_1312, pb_x, pc_x, pc_z, smi0_1309, \
                         smi0_1311, smi0_1312, smh_771, smi1_1309, smi1_1311, smi1_1312, \
                         snh_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = pb_x[k] * smi0_1309[k]
                    - f_10 * pc_x[k] * smi1_1309[k];

        t_1310[k] = f_11 * smh_771[k]
                    + f_3 * pc_z[k] * snh_981[k];

        t_1311[k] = pb_x[k] * smi0_1311[k]
                    - f_10 * pc_x[k] * smi1_1311[k];

        t_1312[k] = pb_x[k] * smi0_1312[k]
                    - f_10 * pc_x[k] * smi1_1312[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, t_1316, pb_x, pc_x, pc_y, smi0_1313, \
                         smi0_1315, smi0_1316, smh_797, smh_987, smi1_1313, smi1_1315, \
                         smi1_1316, snh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = pb_x[k] * smi0_1313[k]
                    - f_10 * pc_x[k] * smi1_1313[k];

        t_1314[k] = f_16 * smh_797[k]
                    + f_3 * pc_y[k] * snh_986[k];

        t_1315[k] = pb_x[k] * smi0_1315[k]
                    - f_10 * pc_x[k] * smi1_1315[k];

        t_1316[k] = pb_x[k] * smi0_1316[k]
                    + f_18 * smh_987[k]
                    - f_10 * pc_x[k] * smi1_1316[k];
    }

#pragma omp simd aligned(t_1317, t_1318, t_1319, t_1320, pb_x, pc_x, pc_y, pc_z, smi0_1319, \
                         smh_777, smh_798, smh_800, smh_990, smi1_1319, snh_987, \
                         snh_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1317[k] = f_17 * smh_798[k]
                    + f_3 * pc_y[k] * snh_987[k];

        t_1318[k] = f_12 * smh_777[k]
                    + f_3 * pc_z[k] * snh_987[k];

        t_1319[k] = pb_x[k] * smi0_1319[k]
                    + f_14 * smh_990[k]
                    - f_10 * pc_x[k] * smi1_1319[k];

        t_1320[k] = f_17 * smh_800[k]
                    + f_3 * pc_y[k] * snh_989[k];
    }

#pragma omp simd aligned(t_1321, t_1322, t_1323, pb_x, pc_x, pc_z, smi0_1321, smi0_1322, \
                         smh_780, smh_992, smh_993, smi1_1321, smi1_1322, \
                         snh_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1321[k] = pb_x[k] * smi0_1321[k]
                    + f_14 * smh_992[k]
                    - f_10 * pc_x[k] * smi1_1321[k];

        t_1322[k] = pb_x[k] * smi0_1322[k]
                    + f_13 * smh_993[k]
                    - f_10 * pc_x[k] * smi1_1322[k];

        t_1323[k] = f_12 * smh_780[k]
                    + f_3 * pc_z[k] * snh_990[k];
    }

#pragma omp simd aligned(t_1324, t_1325, t_1326, pb_x, pc_x, pc_y, smi0_1325, smi0_1326, \
                         smh_803, smh_996, smh_997, smi1_1325, smi1_1326, \
                         snh_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1324[k] = f_17 * smh_803[k]
                    + f_3 * pc_y[k] * snh_992[k];

        t_1325[k] = pb_x[k] * smi0_1325[k]
                    + f_13 * smh_996[k]
                    - f_10 * pc_x[k] * smi1_1325[k];

        t_1326[k] = pb_x[k] * smi0_1326[k]
                    + f_12 * smh_997[k]
                    - f_10 * pc_x[k] * smi1_1326[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, pb_x, pc_x, pc_y, pc_z, smi0_1328, smh_783, \
                         smh_807, smh_999, smi1_1328, snh_993, \
                         snh_996 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = f_12 * smh_783[k]
                    + f_3 * pc_z[k] * snh_993[k];

        t_1328[k] = pb_x[k] * smi0_1328[k]
                    + f_12 * smh_999[k]
                    - f_10 * pc_x[k] * smi1_1328[k];

        t_1329[k] = f_17 * smh_807[k]
                    + f_3 * pc_y[k] * snh_996[k];
    }

#pragma omp simd aligned(t_1330, t_1331, t_1332, t_1333, pb_x, pc_x, smi0_1330, smh_1001, \
                         smh_1002, smh_1003, smh_1004, smi1_1330, snh_1002, snh_1003, \
                         snh_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1330[k] = pb_x[k] * smi0_1330[k]
                    + f_12 * smh_1001[k]
                    - f_10 * pc_x[k] * smi1_1330[k];

        t_1331[k] = f_11 * smh_1002[k]
                    + f_3 * pc_x[k] * snh_1002[k];

        t_1332[k] = f_11 * smh_1003[k]
                    + f_3 * pc_x[k] * snh_1003[k];

        t_1333[k] = f_11 * smh_1004[k]
                    + f_3 * pc_x[k] * snh_1004[k];
    }

#pragma omp simd aligned(t_1334, t_1335, t_1336, t_1337, pb_x, pc_x, smi0_1337, smh_1005, \
                         smh_1006, smh_1007, smi1_1337, snh_1005, snh_1006, \
                         snh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1334[k] = f_11 * smh_1005[k]
                    + f_3 * pc_x[k] * snh_1005[k];

        t_1335[k] = f_11 * smh_1006[k]
                    + f_3 * pc_x[k] * snh_1006[k];

        t_1336[k] = f_11 * smh_1007[k]
                    + f_3 * pc_x[k] * snh_1007[k];

        t_1337[k] = pb_x[k] * smi0_1337[k]
                    - f_10 * pc_x[k] * smi1_1337[k];
    }

#pragma omp simd aligned(t_1338, t_1339, t_1340, t_1341, pb_x, pc_x, pc_z, smi0_1339, \
                         smi0_1340, smi0_1341, smh_792, smi1_1339, smi1_1340, smi1_1341, \
                         snh_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1338[k] = f_12 * smh_792[k]
                    + f_3 * pc_z[k] * snh_1002[k];

        t_1339[k] = pb_x[k] * smi0_1339[k]
                    - f_10 * pc_x[k] * smi1_1339[k];

        t_1340[k] = pb_x[k] * smi0_1340[k]
                    - f_10 * pc_x[k] * smi1_1340[k];

        t_1341[k] = pb_x[k] * smi0_1341[k]
                    - f_10 * pc_x[k] * smi1_1341[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, t_1345, pb_x, pc_x, pc_y, smi0_1343, \
                         smi0_1344, smh_818, smh_819, smh_1008, smi1_1343, smi1_1344, \
                         snh_1007, snh_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_17 * smh_818[k]
                    + f_3 * pc_y[k] * snh_1007[k];

        t_1343[k] = pb_x[k] * smi0_1343[k]
                    - f_10 * pc_x[k] * smi1_1343[k];

        t_1344[k] = pb_x[k] * smi0_1344[k]
                    + f_18 * smh_1008[k]
                    - f_10 * pc_x[k] * smi1_1344[k];

        t_1345[k] = f_18 * smh_819[k]
                    + f_3 * pc_y[k] * snh_1008[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pb_x, pc_x, pc_y, pc_z, smi0_1347, smh_798, \
                         smh_821, smh_1011, smi1_1347, snh_1008, \
                         snh_1010 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_13 * smh_798[k]
                    + f_3 * pc_z[k] * snh_1008[k];

        t_1347[k] = pb_x[k] * smi0_1347[k]
                    + f_14 * smh_1011[k]
                    - f_10 * pc_x[k] * smi1_1347[k];

        t_1348[k] = f_18 * smh_821[k]
                    + f_3 * pc_y[k] * snh_1010[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pb_x, pc_x, pc_z, smi0_1349, smi0_1350, \
                         smh_801, smh_1013, smh_1014, smi1_1349, smi1_1350, \
                         snh_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = pb_x[k] * smi0_1349[k]
                    + f_14 * smh_1013[k]
                    - f_10 * pc_x[k] * smi1_1349[k];

        t_1350[k] = pb_x[k] * smi0_1350[k]
                    + f_13 * smh_1014[k]
                    - f_10 * pc_x[k] * smi1_1350[k];

        t_1351[k] = f_13 * smh_801[k]
                    + f_3 * pc_z[k] * snh_1011[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, pb_x, pc_x, pc_y, smi0_1353, smi0_1354, \
                         smh_824, smh_1017, smh_1018, smi1_1353, smi1_1354, \
                         snh_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_18 * smh_824[k]
                    + f_3 * pc_y[k] * snh_1013[k];

        t_1353[k] = pb_x[k] * smi0_1353[k]
                    + f_13 * smh_1017[k]
                    - f_10 * pc_x[k] * smi1_1353[k];

        t_1354[k] = pb_x[k] * smi0_1354[k]
                    + f_12 * smh_1018[k]
                    - f_10 * pc_x[k] * smi1_1354[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pb_x, pc_x, pc_y, pc_z, smi0_1356, smh_804, \
                         smh_828, smh_1020, smi1_1356, snh_1014, \
                         snh_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_13 * smh_804[k]
                    + f_3 * pc_z[k] * snh_1014[k];

        t_1356[k] = pb_x[k] * smi0_1356[k]
                    + f_12 * smh_1020[k]
                    - f_10 * pc_x[k] * smi1_1356[k];

        t_1357[k] = f_18 * smh_828[k]
                    + f_3 * pc_y[k] * snh_1017[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, t_1361, pb_x, pc_x, smi0_1358, smh_1022, \
                         smh_1023, smh_1024, smh_1025, smi1_1358, snh_1023, snh_1024, \
                         snh_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = pb_x[k] * smi0_1358[k]
                    + f_12 * smh_1022[k]
                    - f_10 * pc_x[k] * smi1_1358[k];

        t_1359[k] = f_11 * smh_1023[k]
                    + f_3 * pc_x[k] * snh_1023[k];

        t_1360[k] = f_11 * smh_1024[k]
                    + f_3 * pc_x[k] * snh_1024[k];

        t_1361[k] = f_11 * smh_1025[k]
                    + f_3 * pc_x[k] * snh_1025[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smi0,
                                                           const size_t smh, const size_t smi1,
                                                           const size_t snh, const size_t ncols,
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
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi0_1365 = buffer.data(smi0 + 1365);
    const auto *smi0_1367 = buffer.data(smi0 + 1367);
    const auto *smi0_1368 = buffer.data(smi0 + 1368);
    const auto *smi0_1369 = buffer.data(smi0 + 1369);
    const auto *smi0_1371 = buffer.data(smi0 + 1371);
    const auto *smi0_1372 = buffer.data(smi0 + 1372);
    const auto *smi0_1375 = buffer.data(smi0 + 1375);
    const auto *smi0_1377 = buffer.data(smi0 + 1377);
    const auto *smi0_1378 = buffer.data(smi0 + 1378);
    const auto *smi0_1381 = buffer.data(smi0 + 1381);
    const auto *smi0_1382 = buffer.data(smi0 + 1382);
    const auto *smi0_1384 = buffer.data(smi0 + 1384);
    const auto *smi0_1386 = buffer.data(smi0 + 1386);
    const auto *smi0_1393 = buffer.data(smi0 + 1393);
    const auto *smi0_1395 = buffer.data(smi0 + 1395);
    const auto *smi0_1396 = buffer.data(smi0 + 1396);
    const auto *smi0_1397 = buffer.data(smi0 + 1397);
    const auto *smi0_1399 = buffer.data(smi0 + 1399);
    const auto *smi0_1400 = buffer.data(smi0 + 1400);
    const auto *smi0_1403 = buffer.data(smi0 + 1403);
    const auto *smi0_1405 = buffer.data(smi0 + 1405);
    const auto *smi0_1406 = buffer.data(smi0 + 1406);
    const auto *smi0_1409 = buffer.data(smi0 + 1409);
    const auto *smi0_1410 = buffer.data(smi0 + 1410);
    const auto *smi0_1412 = buffer.data(smi0 + 1412);
    const auto *smi0_1414 = buffer.data(smi0 + 1414);
    const auto *smi0_1421 = buffer.data(smi0 + 1421);
    const auto *smi0_1423 = buffer.data(smi0 + 1423);
    const auto *smi0_1424 = buffer.data(smi0 + 1424);
    const auto *smi0_1425 = buffer.data(smi0 + 1425);
    const auto *smi0_1427 = buffer.data(smi0 + 1427);
    const auto *smi0_1428 = buffer.data(smi0 + 1428);
    const auto *smi0_1431 = buffer.data(smi0 + 1431);
    const auto *smi0_1433 = buffer.data(smi0 + 1433);
    const auto *smi0_1434 = buffer.data(smi0 + 1434);
    const auto *smi0_1437 = buffer.data(smi0 + 1437);
    const auto *smi0_1438 = buffer.data(smi0 + 1438);
    const auto *smi0_1440 = buffer.data(smi0 + 1440);
    const auto *smi0_1442 = buffer.data(smi0 + 1442);
    const auto *smi0_1449 = buffer.data(smi0 + 1449);
    const auto *smi0_1451 = buffer.data(smi0 + 1451);
    const auto *smi0_1452 = buffer.data(smi0 + 1452);
    const auto *smi0_1453 = buffer.data(smi0 + 1453);
    const auto *smi0_1455 = buffer.data(smi0 + 1455);
    const auto *smi0_1456 = buffer.data(smi0 + 1456);
    const auto *smi0_1459 = buffer.data(smi0 + 1459);
    const auto *smi0_1461 = buffer.data(smi0 + 1461);
    const auto *smi0_1462 = buffer.data(smi0 + 1462);
    const auto *smi0_1465 = buffer.data(smi0 + 1465);
    const auto *smi0_1466 = buffer.data(smi0 + 1466);
    const auto *smi0_1468 = buffer.data(smi0 + 1468);
    const auto *smi0_1470 = buffer.data(smi0 + 1470);
    const auto *smi0_1477 = buffer.data(smi0 + 1477);
    const auto *smi0_1479 = buffer.data(smi0 + 1479);
    const auto *smi0_1480 = buffer.data(smi0 + 1480);
    const auto *smi0_1481 = buffer.data(smi0 + 1481);

    const auto *smh_813 = buffer.data(smh + 813);
    const auto *smh_819 = buffer.data(smh + 819);
    const auto *smh_822 = buffer.data(smh + 822);
    const auto *smh_825 = buffer.data(smh + 825);
    const auto *smh_834 = buffer.data(smh + 834);
    const auto *smh_839 = buffer.data(smh + 839);
    const auto *smh_840 = buffer.data(smh + 840);
    const auto *smh_842 = buffer.data(smh + 842);
    const auto *smh_843 = buffer.data(smh + 843);
    const auto *smh_845 = buffer.data(smh + 845);
    const auto *smh_846 = buffer.data(smh + 846);
    const auto *smh_849 = buffer.data(smh + 849);
    const auto *smh_855 = buffer.data(smh + 855);
    const auto *smh_860 = buffer.data(smh + 860);
    const auto *smh_861 = buffer.data(smh + 861);
    const auto *smh_863 = buffer.data(smh + 863);
    const auto *smh_864 = buffer.data(smh + 864);
    const auto *smh_866 = buffer.data(smh + 866);
    const auto *smh_867 = buffer.data(smh + 867);
    const auto *smh_870 = buffer.data(smh + 870);
    const auto *smh_876 = buffer.data(smh + 876);
    const auto *smh_881 = buffer.data(smh + 881);
    const auto *smh_882 = buffer.data(smh + 882);
    const auto *smh_884 = buffer.data(smh + 884);
    const auto *smh_885 = buffer.data(smh + 885);
    const auto *smh_887 = buffer.data(smh + 887);
    const auto *smh_888 = buffer.data(smh + 888);
    const auto *smh_891 = buffer.data(smh + 891);
    const auto *smh_897 = buffer.data(smh + 897);
    const auto *smh_902 = buffer.data(smh + 902);
    const auto *smh_903 = buffer.data(smh + 903);
    const auto *smh_905 = buffer.data(smh + 905);
    const auto *smh_908 = buffer.data(smh + 908);
    const auto *smh_912 = buffer.data(smh + 912);
    const auto *smh_1026 = buffer.data(smh + 1026);
    const auto *smh_1027 = buffer.data(smh + 1027);
    const auto *smh_1028 = buffer.data(smh + 1028);
    const auto *smh_1029 = buffer.data(smh + 1029);
    const auto *smh_1032 = buffer.data(smh + 1032);
    const auto *smh_1034 = buffer.data(smh + 1034);
    const auto *smh_1035 = buffer.data(smh + 1035);
    const auto *smh_1038 = buffer.data(smh + 1038);
    const auto *smh_1039 = buffer.data(smh + 1039);
    const auto *smh_1041 = buffer.data(smh + 1041);
    const auto *smh_1043 = buffer.data(smh + 1043);
    const auto *smh_1044 = buffer.data(smh + 1044);
    const auto *smh_1045 = buffer.data(smh + 1045);
    const auto *smh_1046 = buffer.data(smh + 1046);
    const auto *smh_1047 = buffer.data(smh + 1047);
    const auto *smh_1048 = buffer.data(smh + 1048);
    const auto *smh_1049 = buffer.data(smh + 1049);
    const auto *smh_1050 = buffer.data(smh + 1050);
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

    const auto *smi1_1365 = buffer.data(smi1 + 1365);
    const auto *smi1_1367 = buffer.data(smi1 + 1367);
    const auto *smi1_1368 = buffer.data(smi1 + 1368);
    const auto *smi1_1369 = buffer.data(smi1 + 1369);
    const auto *smi1_1371 = buffer.data(smi1 + 1371);
    const auto *smi1_1372 = buffer.data(smi1 + 1372);
    const auto *smi1_1375 = buffer.data(smi1 + 1375);
    const auto *smi1_1377 = buffer.data(smi1 + 1377);
    const auto *smi1_1378 = buffer.data(smi1 + 1378);
    const auto *smi1_1381 = buffer.data(smi1 + 1381);
    const auto *smi1_1382 = buffer.data(smi1 + 1382);
    const auto *smi1_1384 = buffer.data(smi1 + 1384);
    const auto *smi1_1386 = buffer.data(smi1 + 1386);
    const auto *smi1_1393 = buffer.data(smi1 + 1393);
    const auto *smi1_1395 = buffer.data(smi1 + 1395);
    const auto *smi1_1396 = buffer.data(smi1 + 1396);
    const auto *smi1_1397 = buffer.data(smi1 + 1397);
    const auto *smi1_1399 = buffer.data(smi1 + 1399);
    const auto *smi1_1400 = buffer.data(smi1 + 1400);
    const auto *smi1_1403 = buffer.data(smi1 + 1403);
    const auto *smi1_1405 = buffer.data(smi1 + 1405);
    const auto *smi1_1406 = buffer.data(smi1 + 1406);
    const auto *smi1_1409 = buffer.data(smi1 + 1409);
    const auto *smi1_1410 = buffer.data(smi1 + 1410);
    const auto *smi1_1412 = buffer.data(smi1 + 1412);
    const auto *smi1_1414 = buffer.data(smi1 + 1414);
    const auto *smi1_1421 = buffer.data(smi1 + 1421);
    const auto *smi1_1423 = buffer.data(smi1 + 1423);
    const auto *smi1_1424 = buffer.data(smi1 + 1424);
    const auto *smi1_1425 = buffer.data(smi1 + 1425);
    const auto *smi1_1427 = buffer.data(smi1 + 1427);
    const auto *smi1_1428 = buffer.data(smi1 + 1428);
    const auto *smi1_1431 = buffer.data(smi1 + 1431);
    const auto *smi1_1433 = buffer.data(smi1 + 1433);
    const auto *smi1_1434 = buffer.data(smi1 + 1434);
    const auto *smi1_1437 = buffer.data(smi1 + 1437);
    const auto *smi1_1438 = buffer.data(smi1 + 1438);
    const auto *smi1_1440 = buffer.data(smi1 + 1440);
    const auto *smi1_1442 = buffer.data(smi1 + 1442);
    const auto *smi1_1449 = buffer.data(smi1 + 1449);
    const auto *smi1_1451 = buffer.data(smi1 + 1451);
    const auto *smi1_1452 = buffer.data(smi1 + 1452);
    const auto *smi1_1453 = buffer.data(smi1 + 1453);
    const auto *smi1_1455 = buffer.data(smi1 + 1455);
    const auto *smi1_1456 = buffer.data(smi1 + 1456);
    const auto *smi1_1459 = buffer.data(smi1 + 1459);
    const auto *smi1_1461 = buffer.data(smi1 + 1461);
    const auto *smi1_1462 = buffer.data(smi1 + 1462);
    const auto *smi1_1465 = buffer.data(smi1 + 1465);
    const auto *smi1_1466 = buffer.data(smi1 + 1466);
    const auto *smi1_1468 = buffer.data(smi1 + 1468);
    const auto *smi1_1470 = buffer.data(smi1 + 1470);
    const auto *smi1_1477 = buffer.data(smi1 + 1477);
    const auto *smi1_1479 = buffer.data(smi1 + 1479);
    const auto *smi1_1480 = buffer.data(smi1 + 1480);
    const auto *smi1_1481 = buffer.data(smi1 + 1481);

    const auto *snh_1023 = buffer.data(snh + 1023);
    const auto *snh_1026 = buffer.data(snh + 1026);
    const auto *snh_1027 = buffer.data(snh + 1027);
    const auto *snh_1028 = buffer.data(snh + 1028);
    const auto *snh_1029 = buffer.data(snh + 1029);
    const auto *snh_1031 = buffer.data(snh + 1031);
    const auto *snh_1032 = buffer.data(snh + 1032);
    const auto *snh_1034 = buffer.data(snh + 1034);
    const auto *snh_1035 = buffer.data(snh + 1035);
    const auto *snh_1038 = buffer.data(snh + 1038);
    const auto *snh_1044 = buffer.data(snh + 1044);
    const auto *snh_1045 = buffer.data(snh + 1045);
    const auto *snh_1046 = buffer.data(snh + 1046);
    const auto *snh_1047 = buffer.data(snh + 1047);
    const auto *snh_1048 = buffer.data(snh + 1048);
    const auto *snh_1049 = buffer.data(snh + 1049);
    const auto *snh_1050 = buffer.data(snh + 1050);
    const auto *snh_1052 = buffer.data(snh + 1052);
    const auto *snh_1053 = buffer.data(snh + 1053);
    const auto *snh_1055 = buffer.data(snh + 1055);
    const auto *snh_1056 = buffer.data(snh + 1056);
    const auto *snh_1059 = buffer.data(snh + 1059);
    const auto *snh_1065 = buffer.data(snh + 1065);
    const auto *snh_1066 = buffer.data(snh + 1066);
    const auto *snh_1067 = buffer.data(snh + 1067);
    const auto *snh_1068 = buffer.data(snh + 1068);
    const auto *snh_1069 = buffer.data(snh + 1069);
    const auto *snh_1070 = buffer.data(snh + 1070);
    const auto *snh_1071 = buffer.data(snh + 1071);
    const auto *snh_1073 = buffer.data(snh + 1073);
    const auto *snh_1074 = buffer.data(snh + 1074);
    const auto *snh_1076 = buffer.data(snh + 1076);
    const auto *snh_1077 = buffer.data(snh + 1077);
    const auto *snh_1080 = buffer.data(snh + 1080);
    const auto *snh_1086 = buffer.data(snh + 1086);
    const auto *snh_1087 = buffer.data(snh + 1087);
    const auto *snh_1088 = buffer.data(snh + 1088);
    const auto *snh_1089 = buffer.data(snh + 1089);
    const auto *snh_1090 = buffer.data(snh + 1090);
    const auto *snh_1091 = buffer.data(snh + 1091);
    const auto *snh_1092 = buffer.data(snh + 1092);
    const auto *snh_1094 = buffer.data(snh + 1094);
    const auto *snh_1095 = buffer.data(snh + 1095);
    const auto *snh_1097 = buffer.data(snh + 1097);
    const auto *snh_1098 = buffer.data(snh + 1098);
    const auto *snh_1101 = buffer.data(snh + 1101);
    const auto *snh_1107 = buffer.data(snh + 1107);
    const auto *snh_1108 = buffer.data(snh + 1108);
    const auto *snh_1109 = buffer.data(snh + 1109);
    const auto *snh_1110 = buffer.data(snh + 1110);
    const auto *snh_1111 = buffer.data(snh + 1111);
    const auto *snh_1112 = buffer.data(snh + 1112);

#pragma omp simd aligned(t_1362, t_1363, t_1364, t_1365, pb_x, pc_x, smi0_1365, smh_1026, \
                         smh_1027, smh_1028, smi1_1365, snh_1026, snh_1027, \
                         snh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = f_11 * smh_1026[k]
                    + f_3 * pc_x[k] * snh_1026[k];

        t_1363[k] = f_11 * smh_1027[k]
                    + f_3 * pc_x[k] * snh_1027[k];

        t_1364[k] = f_11 * smh_1028[k]
                    + f_3 * pc_x[k] * snh_1028[k];

        t_1365[k] = pb_x[k] * smi0_1365[k]
                    - f_10 * pc_x[k] * smi1_1365[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, t_1369, pb_x, pc_x, pc_z, smi0_1367, \
                         smi0_1368, smi0_1369, smh_813, smi1_1367, smi1_1368, smi1_1369, \
                         snh_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_13 * smh_813[k]
                    + f_3 * pc_z[k] * snh_1023[k];

        t_1367[k] = pb_x[k] * smi0_1367[k]
                    - f_10 * pc_x[k] * smi1_1367[k];

        t_1368[k] = pb_x[k] * smi0_1368[k]
                    - f_10 * pc_x[k] * smi1_1368[k];

        t_1369[k] = pb_x[k] * smi0_1369[k]
                    - f_10 * pc_x[k] * smi1_1369[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, t_1373, pb_x, pc_x, pc_y, smi0_1371, \
                         smi0_1372, smh_839, smh_840, smh_1029, smi1_1371, smi1_1372, \
                         snh_1028, snh_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_18 * smh_839[k]
                    + f_3 * pc_y[k] * snh_1028[k];

        t_1371[k] = pb_x[k] * smi0_1371[k]
                    - f_10 * pc_x[k] * smi1_1371[k];

        t_1372[k] = pb_x[k] * smi0_1372[k]
                    + f_18 * smh_1029[k]
                    - f_10 * pc_x[k] * smi1_1372[k];

        t_1373[k] = f_19 * smh_840[k]
                    + f_3 * pc_y[k] * snh_1029[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, pb_x, pc_x, pc_y, pc_z, smi0_1375, smh_819, \
                         smh_842, smh_1032, smi1_1375, snh_1029, \
                         snh_1031 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_14 * smh_819[k]
                    + f_3 * pc_z[k] * snh_1029[k];

        t_1375[k] = pb_x[k] * smi0_1375[k]
                    + f_14 * smh_1032[k]
                    - f_10 * pc_x[k] * smi1_1375[k];

        t_1376[k] = f_19 * smh_842[k]
                    + f_3 * pc_y[k] * snh_1031[k];
    }

#pragma omp simd aligned(t_1377, t_1378, t_1379, pb_x, pc_x, pc_z, smi0_1377, smi0_1378, \
                         smh_822, smh_1034, smh_1035, smi1_1377, smi1_1378, \
                         snh_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1377[k] = pb_x[k] * smi0_1377[k]
                    + f_14 * smh_1034[k]
                    - f_10 * pc_x[k] * smi1_1377[k];

        t_1378[k] = pb_x[k] * smi0_1378[k]
                    + f_13 * smh_1035[k]
                    - f_10 * pc_x[k] * smi1_1378[k];

        t_1379[k] = f_14 * smh_822[k]
                    + f_3 * pc_z[k] * snh_1032[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, pb_x, pc_x, pc_y, smi0_1381, smi0_1382, \
                         smh_845, smh_1038, smh_1039, smi1_1381, smi1_1382, \
                         snh_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_19 * smh_845[k]
                    + f_3 * pc_y[k] * snh_1034[k];

        t_1381[k] = pb_x[k] * smi0_1381[k]
                    + f_13 * smh_1038[k]
                    - f_10 * pc_x[k] * smi1_1381[k];

        t_1382[k] = pb_x[k] * smi0_1382[k]
                    + f_12 * smh_1039[k]
                    - f_10 * pc_x[k] * smi1_1382[k];
    }

#pragma omp simd aligned(t_1383, t_1384, t_1385, pb_x, pc_x, pc_y, pc_z, smi0_1384, smh_825, \
                         smh_849, smh_1041, smi1_1384, snh_1035, \
                         snh_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1383[k] = f_14 * smh_825[k]
                    + f_3 * pc_z[k] * snh_1035[k];

        t_1384[k] = pb_x[k] * smi0_1384[k]
                    + f_12 * smh_1041[k]
                    - f_10 * pc_x[k] * smi1_1384[k];

        t_1385[k] = f_19 * smh_849[k]
                    + f_3 * pc_y[k] * snh_1038[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, t_1389, pb_x, pc_x, smi0_1386, smh_1043, \
                         smh_1044, smh_1045, smh_1046, smi1_1386, snh_1044, snh_1045, \
                         snh_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = pb_x[k] * smi0_1386[k]
                    + f_12 * smh_1043[k]
                    - f_10 * pc_x[k] * smi1_1386[k];

        t_1387[k] = f_11 * smh_1044[k]
                    + f_3 * pc_x[k] * snh_1044[k];

        t_1388[k] = f_11 * smh_1045[k]
                    + f_3 * pc_x[k] * snh_1045[k];

        t_1389[k] = f_11 * smh_1046[k]
                    + f_3 * pc_x[k] * snh_1046[k];
    }

#pragma omp simd aligned(t_1390, t_1391, t_1392, t_1393, pb_x, pc_x, smi0_1393, smh_1047, \
                         smh_1048, smh_1049, smi1_1393, snh_1047, snh_1048, \
                         snh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1390[k] = f_11 * smh_1047[k]
                    + f_3 * pc_x[k] * snh_1047[k];

        t_1391[k] = f_11 * smh_1048[k]
                    + f_3 * pc_x[k] * snh_1048[k];

        t_1392[k] = f_11 * smh_1049[k]
                    + f_3 * pc_x[k] * snh_1049[k];

        t_1393[k] = pb_x[k] * smi0_1393[k]
                    - f_10 * pc_x[k] * smi1_1393[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, t_1397, pb_x, pc_x, pc_z, smi0_1395, \
                         smi0_1396, smi0_1397, smh_834, smi1_1395, smi1_1396, smi1_1397, \
                         snh_1044 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_14 * smh_834[k]
                    + f_3 * pc_z[k] * snh_1044[k];

        t_1395[k] = pb_x[k] * smi0_1395[k]
                    - f_10 * pc_x[k] * smi1_1395[k];

        t_1396[k] = pb_x[k] * smi0_1396[k]
                    - f_10 * pc_x[k] * smi1_1396[k];

        t_1397[k] = pb_x[k] * smi0_1397[k]
                    - f_10 * pc_x[k] * smi1_1397[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, t_1401, pb_x, pc_x, pc_y, smi0_1399, \
                         smi0_1400, smh_860, smh_861, smh_1050, smi1_1399, smi1_1400, \
                         snh_1049, snh_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_19 * smh_860[k]
                    + f_3 * pc_y[k] * snh_1049[k];

        t_1399[k] = pb_x[k] * smi0_1399[k]
                    - f_10 * pc_x[k] * smi1_1399[k];

        t_1400[k] = pb_x[k] * smi0_1400[k]
                    + f_18 * smh_1050[k]
                    - f_10 * pc_x[k] * smi1_1400[k];

        t_1401[k] = f_14 * smh_861[k]
                    + f_3 * pc_y[k] * snh_1050[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, pb_x, pc_x, pc_y, pc_z, smi0_1403, smh_840, \
                         smh_863, smh_1053, smi1_1403, snh_1050, \
                         snh_1052 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_19 * smh_840[k]
                    + f_3 * pc_z[k] * snh_1050[k];

        t_1403[k] = pb_x[k] * smi0_1403[k]
                    + f_14 * smh_1053[k]
                    - f_10 * pc_x[k] * smi1_1403[k];

        t_1404[k] = f_14 * smh_863[k]
                    + f_3 * pc_y[k] * snh_1052[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, pb_x, pc_x, pc_z, smi0_1405, smi0_1406, \
                         smh_843, smh_1055, smh_1056, smi1_1405, smi1_1406, \
                         snh_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = pb_x[k] * smi0_1405[k]
                    + f_14 * smh_1055[k]
                    - f_10 * pc_x[k] * smi1_1405[k];

        t_1406[k] = pb_x[k] * smi0_1406[k]
                    + f_13 * smh_1056[k]
                    - f_10 * pc_x[k] * smi1_1406[k];

        t_1407[k] = f_19 * smh_843[k]
                    + f_3 * pc_z[k] * snh_1053[k];
    }

#pragma omp simd aligned(t_1408, t_1409, t_1410, pb_x, pc_x, pc_y, smi0_1409, smi0_1410, \
                         smh_866, smh_1059, smh_1060, smi1_1409, smi1_1410, \
                         snh_1055 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1408[k] = f_14 * smh_866[k]
                    + f_3 * pc_y[k] * snh_1055[k];

        t_1409[k] = pb_x[k] * smi0_1409[k]
                    + f_13 * smh_1059[k]
                    - f_10 * pc_x[k] * smi1_1409[k];

        t_1410[k] = pb_x[k] * smi0_1410[k]
                    + f_12 * smh_1060[k]
                    - f_10 * pc_x[k] * smi1_1410[k];
    }

#pragma omp simd aligned(t_1411, t_1412, t_1413, pb_x, pc_x, pc_y, pc_z, smi0_1412, smh_846, \
                         smh_870, smh_1062, smi1_1412, snh_1056, \
                         snh_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1411[k] = f_19 * smh_846[k]
                    + f_3 * pc_z[k] * snh_1056[k];

        t_1412[k] = pb_x[k] * smi0_1412[k]
                    + f_12 * smh_1062[k]
                    - f_10 * pc_x[k] * smi1_1412[k];

        t_1413[k] = f_14 * smh_870[k]
                    + f_3 * pc_y[k] * snh_1059[k];
    }

#pragma omp simd aligned(t_1414, t_1415, t_1416, t_1417, pb_x, pc_x, smi0_1414, smh_1064, \
                         smh_1065, smh_1066, smh_1067, smi1_1414, snh_1065, snh_1066, \
                         snh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1414[k] = pb_x[k] * smi0_1414[k]
                    + f_12 * smh_1064[k]
                    - f_10 * pc_x[k] * smi1_1414[k];

        t_1415[k] = f_11 * smh_1065[k]
                    + f_3 * pc_x[k] * snh_1065[k];

        t_1416[k] = f_11 * smh_1066[k]
                    + f_3 * pc_x[k] * snh_1066[k];

        t_1417[k] = f_11 * smh_1067[k]
                    + f_3 * pc_x[k] * snh_1067[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, t_1421, pb_x, pc_x, smi0_1421, smh_1068, \
                         smh_1069, smh_1070, smi1_1421, snh_1068, snh_1069, \
                         snh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_11 * smh_1068[k]
                    + f_3 * pc_x[k] * snh_1068[k];

        t_1419[k] = f_11 * smh_1069[k]
                    + f_3 * pc_x[k] * snh_1069[k];

        t_1420[k] = f_11 * smh_1070[k]
                    + f_3 * pc_x[k] * snh_1070[k];

        t_1421[k] = pb_x[k] * smi0_1421[k]
                    - f_10 * pc_x[k] * smi1_1421[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, t_1425, pb_x, pc_x, pc_z, smi0_1423, \
                         smi0_1424, smi0_1425, smh_855, smi1_1423, smi1_1424, smi1_1425, \
                         snh_1065 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_19 * smh_855[k]
                    + f_3 * pc_z[k] * snh_1065[k];

        t_1423[k] = pb_x[k] * smi0_1423[k]
                    - f_10 * pc_x[k] * smi1_1423[k];

        t_1424[k] = pb_x[k] * smi0_1424[k]
                    - f_10 * pc_x[k] * smi1_1424[k];

        t_1425[k] = pb_x[k] * smi0_1425[k]
                    - f_10 * pc_x[k] * smi1_1425[k];
    }

#pragma omp simd aligned(t_1426, t_1427, t_1428, t_1429, pb_x, pc_x, pc_y, smi0_1427, \
                         smi0_1428, smh_881, smh_882, smh_1071, smi1_1427, smi1_1428, \
                         snh_1070, snh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1426[k] = f_14 * smh_881[k]
                    + f_3 * pc_y[k] * snh_1070[k];

        t_1427[k] = pb_x[k] * smi0_1427[k]
                    - f_10 * pc_x[k] * smi1_1427[k];

        t_1428[k] = pb_x[k] * smi0_1428[k]
                    + f_18 * smh_1071[k]
                    - f_10 * pc_x[k] * smi1_1428[k];

        t_1429[k] = f_13 * smh_882[k]
                    + f_3 * pc_y[k] * snh_1071[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, pb_x, pc_x, pc_y, pc_z, smi0_1431, smh_861, \
                         smh_884, smh_1074, smi1_1431, snh_1071, \
                         snh_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_18 * smh_861[k]
                    + f_3 * pc_z[k] * snh_1071[k];

        t_1431[k] = pb_x[k] * smi0_1431[k]
                    + f_14 * smh_1074[k]
                    - f_10 * pc_x[k] * smi1_1431[k];

        t_1432[k] = f_13 * smh_884[k]
                    + f_3 * pc_y[k] * snh_1073[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pb_x, pc_x, pc_z, smi0_1433, smi0_1434, \
                         smh_864, smh_1076, smh_1077, smi1_1433, smi1_1434, \
                         snh_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = pb_x[k] * smi0_1433[k]
                    + f_14 * smh_1076[k]
                    - f_10 * pc_x[k] * smi1_1433[k];

        t_1434[k] = pb_x[k] * smi0_1434[k]
                    + f_13 * smh_1077[k]
                    - f_10 * pc_x[k] * smi1_1434[k];

        t_1435[k] = f_18 * smh_864[k]
                    + f_3 * pc_z[k] * snh_1074[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, pb_x, pc_x, pc_y, smi0_1437, smi0_1438, \
                         smh_887, smh_1080, smh_1081, smi1_1437, smi1_1438, \
                         snh_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_13 * smh_887[k]
                    + f_3 * pc_y[k] * snh_1076[k];

        t_1437[k] = pb_x[k] * smi0_1437[k]
                    + f_13 * smh_1080[k]
                    - f_10 * pc_x[k] * smi1_1437[k];

        t_1438[k] = pb_x[k] * smi0_1438[k]
                    + f_12 * smh_1081[k]
                    - f_10 * pc_x[k] * smi1_1438[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, pb_x, pc_x, pc_y, pc_z, smi0_1440, smh_867, \
                         smh_891, smh_1083, smi1_1440, snh_1077, \
                         snh_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = f_18 * smh_867[k]
                    + f_3 * pc_z[k] * snh_1077[k];

        t_1440[k] = pb_x[k] * smi0_1440[k]
                    + f_12 * smh_1083[k]
                    - f_10 * pc_x[k] * smi1_1440[k];

        t_1441[k] = f_13 * smh_891[k]
                    + f_3 * pc_y[k] * snh_1080[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, t_1445, pb_x, pc_x, smi0_1442, smh_1085, \
                         smh_1086, smh_1087, smh_1088, smi1_1442, snh_1086, snh_1087, \
                         snh_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = pb_x[k] * smi0_1442[k]
                    + f_12 * smh_1085[k]
                    - f_10 * pc_x[k] * smi1_1442[k];

        t_1443[k] = f_11 * smh_1086[k]
                    + f_3 * pc_x[k] * snh_1086[k];

        t_1444[k] = f_11 * smh_1087[k]
                    + f_3 * pc_x[k] * snh_1087[k];

        t_1445[k] = f_11 * smh_1088[k]
                    + f_3 * pc_x[k] * snh_1088[k];
    }

#pragma omp simd aligned(t_1446, t_1447, t_1448, t_1449, pb_x, pc_x, smi0_1449, smh_1089, \
                         smh_1090, smh_1091, smi1_1449, snh_1089, snh_1090, \
                         snh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1446[k] = f_11 * smh_1089[k]
                    + f_3 * pc_x[k] * snh_1089[k];

        t_1447[k] = f_11 * smh_1090[k]
                    + f_3 * pc_x[k] * snh_1090[k];

        t_1448[k] = f_11 * smh_1091[k]
                    + f_3 * pc_x[k] * snh_1091[k];

        t_1449[k] = pb_x[k] * smi0_1449[k]
                    - f_10 * pc_x[k] * smi1_1449[k];
    }

#pragma omp simd aligned(t_1450, t_1451, t_1452, t_1453, pb_x, pc_x, pc_z, smi0_1451, \
                         smi0_1452, smi0_1453, smh_876, smi1_1451, smi1_1452, smi1_1453, \
                         snh_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1450[k] = f_18 * smh_876[k]
                    + f_3 * pc_z[k] * snh_1086[k];

        t_1451[k] = pb_x[k] * smi0_1451[k]
                    - f_10 * pc_x[k] * smi1_1451[k];

        t_1452[k] = pb_x[k] * smi0_1452[k]
                    - f_10 * pc_x[k] * smi1_1452[k];

        t_1453[k] = pb_x[k] * smi0_1453[k]
                    - f_10 * pc_x[k] * smi1_1453[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, t_1457, pb_x, pc_x, pc_y, smi0_1455, \
                         smi0_1456, smh_902, smh_903, smh_1092, smi1_1455, smi1_1456, \
                         snh_1091, snh_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_13 * smh_902[k]
                    + f_3 * pc_y[k] * snh_1091[k];

        t_1455[k] = pb_x[k] * smi0_1455[k]
                    - f_10 * pc_x[k] * smi1_1455[k];

        t_1456[k] = pb_x[k] * smi0_1456[k]
                    + f_18 * smh_1092[k]
                    - f_10 * pc_x[k] * smi1_1456[k];

        t_1457[k] = f_12 * smh_903[k]
                    + f_3 * pc_y[k] * snh_1092[k];
    }

#pragma omp simd aligned(t_1458, t_1459, t_1460, pb_x, pc_x, pc_y, pc_z, smi0_1459, smh_882, \
                         smh_905, smh_1095, smi1_1459, snh_1092, \
                         snh_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1458[k] = f_17 * smh_882[k]
                    + f_3 * pc_z[k] * snh_1092[k];

        t_1459[k] = pb_x[k] * smi0_1459[k]
                    + f_14 * smh_1095[k]
                    - f_10 * pc_x[k] * smi1_1459[k];

        t_1460[k] = f_12 * smh_905[k]
                    + f_3 * pc_y[k] * snh_1094[k];
    }

#pragma omp simd aligned(t_1461, t_1462, t_1463, pb_x, pc_x, pc_z, smi0_1461, smi0_1462, \
                         smh_885, smh_1097, smh_1098, smi1_1461, smi1_1462, \
                         snh_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1461[k] = pb_x[k] * smi0_1461[k]
                    + f_14 * smh_1097[k]
                    - f_10 * pc_x[k] * smi1_1461[k];

        t_1462[k] = pb_x[k] * smi0_1462[k]
                    + f_13 * smh_1098[k]
                    - f_10 * pc_x[k] * smi1_1462[k];

        t_1463[k] = f_17 * smh_885[k]
                    + f_3 * pc_z[k] * snh_1095[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, pb_x, pc_x, pc_y, smi0_1465, smi0_1466, \
                         smh_908, smh_1101, smh_1102, smi1_1465, smi1_1466, \
                         snh_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_12 * smh_908[k]
                    + f_3 * pc_y[k] * snh_1097[k];

        t_1465[k] = pb_x[k] * smi0_1465[k]
                    + f_13 * smh_1101[k]
                    - f_10 * pc_x[k] * smi1_1465[k];

        t_1466[k] = pb_x[k] * smi0_1466[k]
                    + f_12 * smh_1102[k]
                    - f_10 * pc_x[k] * smi1_1466[k];
    }

#pragma omp simd aligned(t_1467, t_1468, t_1469, pb_x, pc_x, pc_y, pc_z, smi0_1468, smh_888, \
                         smh_912, smh_1104, smi1_1468, snh_1098, \
                         snh_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1467[k] = f_17 * smh_888[k]
                    + f_3 * pc_z[k] * snh_1098[k];

        t_1468[k] = pb_x[k] * smi0_1468[k]
                    + f_12 * smh_1104[k]
                    - f_10 * pc_x[k] * smi1_1468[k];

        t_1469[k] = f_12 * smh_912[k]
                    + f_3 * pc_y[k] * snh_1101[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, pb_x, pc_x, smi0_1470, smh_1106, \
                         smh_1107, smh_1108, smh_1109, smi1_1470, snh_1107, snh_1108, \
                         snh_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = pb_x[k] * smi0_1470[k]
                    + f_12 * smh_1106[k]
                    - f_10 * pc_x[k] * smi1_1470[k];

        t_1471[k] = f_11 * smh_1107[k]
                    + f_3 * pc_x[k] * snh_1107[k];

        t_1472[k] = f_11 * smh_1108[k]
                    + f_3 * pc_x[k] * snh_1108[k];

        t_1473[k] = f_11 * smh_1109[k]
                    + f_3 * pc_x[k] * snh_1109[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, pb_x, pc_x, smi0_1477, smh_1110, \
                         smh_1111, smh_1112, smi1_1477, snh_1110, snh_1111, \
                         snh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_11 * smh_1110[k]
                    + f_3 * pc_x[k] * snh_1110[k];

        t_1475[k] = f_11 * smh_1111[k]
                    + f_3 * pc_x[k] * snh_1111[k];

        t_1476[k] = f_11 * smh_1112[k]
                    + f_3 * pc_x[k] * snh_1112[k];

        t_1477[k] = pb_x[k] * smi0_1477[k]
                    - f_10 * pc_x[k] * smi1_1477[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, t_1481, pb_x, pc_x, pc_z, smi0_1479, \
                         smi0_1480, smi0_1481, smh_897, smi1_1479, smi1_1480, smi1_1481, \
                         snh_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_17 * smh_897[k]
                    + f_3 * pc_z[k] * snh_1107[k];

        t_1479[k] = pb_x[k] * smi0_1479[k]
                    - f_10 * pc_x[k] * smi1_1479[k];

        t_1480[k] = pb_x[k] * smi0_1480[k]
                    - f_10 * pc_x[k] * smi1_1480[k];

        t_1481[k] = pb_x[k] * smi0_1481[k]
                    - f_10 * pc_x[k] * smi1_1481[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smi0,
                                                           const size_t smh, const size_t smi1,
                                                           const size_t sng0, const size_t sng1,
                                                           const size_t snh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_18 = 3.0 / q;

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
    auto *t_1540 = buffer.data(target + 1540);
    auto *t_1541 = buffer.data(target + 1541);
    auto *t_1542 = buffer.data(target + 1542);
    auto *t_1543 = buffer.data(target + 1543);
    auto *t_1544 = buffer.data(target + 1544);
    auto *t_1545 = buffer.data(target + 1545);
    auto *t_1546 = buffer.data(target + 1546);
    auto *t_1547 = buffer.data(target + 1547);
    auto *t_1548 = buffer.data(target + 1548);
    auto *t_1549 = buffer.data(target + 1549);
    auto *t_1550 = buffer.data(target + 1550);
    auto *t_1551 = buffer.data(target + 1551);
    auto *t_1552 = buffer.data(target + 1552);
    auto *t_1553 = buffer.data(target + 1553);
    auto *t_1554 = buffer.data(target + 1554);
    auto *t_1555 = buffer.data(target + 1555);
    auto *t_1556 = buffer.data(target + 1556);
    auto *t_1557 = buffer.data(target + 1557);
    auto *t_1558 = buffer.data(target + 1558);
    auto *t_1559 = buffer.data(target + 1559);
    auto *t_1560 = buffer.data(target + 1560);
    auto *t_1561 = buffer.data(target + 1561);
    auto *t_1562 = buffer.data(target + 1562);
    auto *t_1563 = buffer.data(target + 1563);
    auto *t_1564 = buffer.data(target + 1564);
    auto *t_1565 = buffer.data(target + 1565);
    auto *t_1566 = buffer.data(target + 1566);
    auto *t_1567 = buffer.data(target + 1567);
    auto *t_1568 = buffer.data(target + 1568);
    auto *t_1569 = buffer.data(target + 1569);
    auto *t_1570 = buffer.data(target + 1570);
    auto *t_1571 = buffer.data(target + 1571);
    auto *t_1572 = buffer.data(target + 1572);
    auto *t_1573 = buffer.data(target + 1573);
    auto *t_1574 = buffer.data(target + 1574);
    auto *t_1575 = buffer.data(target + 1575);
    auto *t_1576 = buffer.data(target + 1576);
    auto *t_1577 = buffer.data(target + 1577);
    auto *t_1578 = buffer.data(target + 1578);
    auto *t_1579 = buffer.data(target + 1579);
    auto *t_1580 = buffer.data(target + 1580);
    auto *t_1581 = buffer.data(target + 1581);
    auto *t_1582 = buffer.data(target + 1582);
    auto *t_1583 = buffer.data(target + 1583);
    auto *t_1584 = buffer.data(target + 1584);
    auto *t_1585 = buffer.data(target + 1585);
    auto *t_1586 = buffer.data(target + 1586);
    auto *t_1587 = buffer.data(target + 1587);
    auto *t_1588 = buffer.data(target + 1588);
    auto *t_1589 = buffer.data(target + 1589);
    auto *t_1590 = buffer.data(target + 1590);
    auto *t_1591 = buffer.data(target + 1591);
    auto *t_1592 = buffer.data(target + 1592);
    auto *t_1593 = buffer.data(target + 1593);
    auto *t_1594 = buffer.data(target + 1594);
    auto *t_1595 = buffer.data(target + 1595);
    auto *t_1596 = buffer.data(target + 1596);
    auto *t_1597 = buffer.data(target + 1597);
    auto *t_1598 = buffer.data(target + 1598);
    auto *t_1599 = buffer.data(target + 1599);
    auto *t_1600 = buffer.data(target + 1600);
    auto *t_1601 = buffer.data(target + 1601);
    auto *t_1602 = buffer.data(target + 1602);
    auto *t_1603 = buffer.data(target + 1603);
    auto *t_1604 = buffer.data(target + 1604);
    auto *t_1605 = buffer.data(target + 1605);
    auto *t_1606 = buffer.data(target + 1606);
    auto *t_1607 = buffer.data(target + 1607);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi0_1232 = buffer.data(smi0 + 1232);
    const auto *smi0_1237 = buffer.data(smi0 + 1237);
    const auto *smi0_1241 = buffer.data(smi0 + 1241);
    const auto *smi0_1246 = buffer.data(smi0 + 1246);
    const auto *smi0_1260 = buffer.data(smi0 + 1260);
    const auto *smi0_1263 = buffer.data(smi0 + 1263);
    const auto *smi0_1266 = buffer.data(smi0 + 1266);
    const auto *smi0_1270 = buffer.data(smi0 + 1270);
    const auto *smi0_1281 = buffer.data(smi0 + 1281);
    const auto *smi0_1283 = buffer.data(smi0 + 1283);
    const auto *smi0_1284 = buffer.data(smi0 + 1284);
    const auto *smi0_1285 = buffer.data(smi0 + 1285);
    const auto *smi0_1483 = buffer.data(smi0 + 1483);
    const auto *smi0_1487 = buffer.data(smi0 + 1487);
    const auto *smi0_1490 = buffer.data(smi0 + 1490);
    const auto *smi0_1494 = buffer.data(smi0 + 1494);
    const auto *smi0_1496 = buffer.data(smi0 + 1496);
    const auto *smi0_1505 = buffer.data(smi0 + 1505);
    const auto *smi0_1507 = buffer.data(smi0 + 1507);
    const auto *smi0_1508 = buffer.data(smi0 + 1508);
    const auto *smi0_1509 = buffer.data(smi0 + 1509);
    const auto *smi0_1511 = buffer.data(smi0 + 1511);
    const auto *smi0_1512 = buffer.data(smi0 + 1512);
    const auto *smi0_1515 = buffer.data(smi0 + 1515);
    const auto *smi0_1517 = buffer.data(smi0 + 1517);
    const auto *smi0_1518 = buffer.data(smi0 + 1518);
    const auto *smi0_1521 = buffer.data(smi0 + 1521);
    const auto *smi0_1522 = buffer.data(smi0 + 1522);
    const auto *smi0_1524 = buffer.data(smi0 + 1524);
    const auto *smi0_1526 = buffer.data(smi0 + 1526);
    const auto *smi0_1533 = buffer.data(smi0 + 1533);
    const auto *smi0_1535 = buffer.data(smi0 + 1535);
    const auto *smi0_1536 = buffer.data(smi0 + 1536);
    const auto *smi0_1537 = buffer.data(smi0 + 1537);
    const auto *smi0_1539 = buffer.data(smi0 + 1539);

    const auto *smh_903 = buffer.data(smh + 903);
    const auto *smh_906 = buffer.data(smh + 906);
    const auto *smh_909 = buffer.data(smh + 909);
    const auto *smh_918 = buffer.data(smh + 918);
    const auto *smh_923 = buffer.data(smh + 923);
    const auto *smh_924 = buffer.data(smh + 924);
    const auto *smh_926 = buffer.data(smh + 926);
    const auto *smh_927 = buffer.data(smh + 927);
    const auto *smh_929 = buffer.data(smh + 929);
    const auto *smh_930 = buffer.data(smh + 930);
    const auto *smh_933 = buffer.data(smh + 933);
    const auto *smh_939 = buffer.data(smh + 939);
    const auto *smh_944 = buffer.data(smh + 944);
    const auto *smh_945 = buffer.data(smh + 945);
    const auto *smh_947 = buffer.data(smh + 947);
    const auto *smh_948 = buffer.data(smh + 948);
    const auto *smh_950 = buffer.data(smh + 950);
    const auto *smh_951 = buffer.data(smh + 951);
    const auto *smh_954 = buffer.data(smh + 954);
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
    const auto *smh_986 = buffer.data(smh + 986);
    const auto *smh_987 = buffer.data(smh + 987);
    const auto *smh_989 = buffer.data(smh + 989);
    const auto *smh_992 = buffer.data(smh + 992);
    const auto *smh_1116 = buffer.data(smh + 1116);
    const auto *smh_1119 = buffer.data(smh + 1119);
    const auto *smh_1123 = buffer.data(smh + 1123);
    const auto *smh_1125 = buffer.data(smh + 1125);
    const auto *smh_1128 = buffer.data(smh + 1128);
    const auto *smh_1129 = buffer.data(smh + 1129);
    const auto *smh_1130 = buffer.data(smh + 1130);
    const auto *smh_1131 = buffer.data(smh + 1131);
    const auto *smh_1132 = buffer.data(smh + 1132);
    const auto *smh_1133 = buffer.data(smh + 1133);
    const auto *smh_1134 = buffer.data(smh + 1134);
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

    const auto *smi1_1232 = buffer.data(smi1 + 1232);
    const auto *smi1_1237 = buffer.data(smi1 + 1237);
    const auto *smi1_1241 = buffer.data(smi1 + 1241);
    const auto *smi1_1246 = buffer.data(smi1 + 1246);
    const auto *smi1_1260 = buffer.data(smi1 + 1260);
    const auto *smi1_1263 = buffer.data(smi1 + 1263);
    const auto *smi1_1266 = buffer.data(smi1 + 1266);
    const auto *smi1_1270 = buffer.data(smi1 + 1270);
    const auto *smi1_1281 = buffer.data(smi1 + 1281);
    const auto *smi1_1283 = buffer.data(smi1 + 1283);
    const auto *smi1_1284 = buffer.data(smi1 + 1284);
    const auto *smi1_1285 = buffer.data(smi1 + 1285);
    const auto *smi1_1483 = buffer.data(smi1 + 1483);
    const auto *smi1_1487 = buffer.data(smi1 + 1487);
    const auto *smi1_1490 = buffer.data(smi1 + 1490);
    const auto *smi1_1494 = buffer.data(smi1 + 1494);
    const auto *smi1_1496 = buffer.data(smi1 + 1496);
    const auto *smi1_1505 = buffer.data(smi1 + 1505);
    const auto *smi1_1507 = buffer.data(smi1 + 1507);
    const auto *smi1_1508 = buffer.data(smi1 + 1508);
    const auto *smi1_1509 = buffer.data(smi1 + 1509);
    const auto *smi1_1511 = buffer.data(smi1 + 1511);
    const auto *smi1_1512 = buffer.data(smi1 + 1512);
    const auto *smi1_1515 = buffer.data(smi1 + 1515);
    const auto *smi1_1517 = buffer.data(smi1 + 1517);
    const auto *smi1_1518 = buffer.data(smi1 + 1518);
    const auto *smi1_1521 = buffer.data(smi1 + 1521);
    const auto *smi1_1522 = buffer.data(smi1 + 1522);
    const auto *smi1_1524 = buffer.data(smi1 + 1524);
    const auto *smi1_1526 = buffer.data(smi1 + 1526);
    const auto *smi1_1533 = buffer.data(smi1 + 1533);
    const auto *smi1_1535 = buffer.data(smi1 + 1535);
    const auto *smi1_1536 = buffer.data(smi1 + 1536);
    const auto *smi1_1537 = buffer.data(smi1 + 1537);
    const auto *smi1_1539 = buffer.data(smi1 + 1539);

    const auto *sng0_825 = buffer.data(sng0 + 825);
    const auto *sng0_828 = buffer.data(sng0 + 828);
    const auto *sng0_830 = buffer.data(sng0 + 830);
    const auto *sng0_831 = buffer.data(sng0 + 831);
    const auto *sng0_834 = buffer.data(sng0 + 834);
    const auto *sng0_835 = buffer.data(sng0 + 835);
    const auto *sng0_837 = buffer.data(sng0 + 837);
    const auto *sng0_838 = buffer.data(sng0 + 838);
    const auto *sng0_839 = buffer.data(sng0 + 839);
    const auto *sng0_845 = buffer.data(sng0 + 845);
    const auto *sng0_849 = buffer.data(sng0 + 849);
    const auto *sng0_852 = buffer.data(sng0 + 852);
    const auto *sng0_854 = buffer.data(sng0 + 854);
    const auto *sng0_855 = buffer.data(sng0 + 855);
    const auto *sng0_858 = buffer.data(sng0 + 858);
    const auto *sng0_860 = buffer.data(sng0 + 860);
    const auto *sng0_861 = buffer.data(sng0 + 861);
    const auto *sng0_864 = buffer.data(sng0 + 864);
    const auto *sng0_865 = buffer.data(sng0 + 865);

    const auto *sng1_825 = buffer.data(sng1 + 825);
    const auto *sng1_828 = buffer.data(sng1 + 828);
    const auto *sng1_830 = buffer.data(sng1 + 830);
    const auto *sng1_831 = buffer.data(sng1 + 831);
    const auto *sng1_834 = buffer.data(sng1 + 834);
    const auto *sng1_835 = buffer.data(sng1 + 835);
    const auto *sng1_837 = buffer.data(sng1 + 837);
    const auto *sng1_838 = buffer.data(sng1 + 838);
    const auto *sng1_839 = buffer.data(sng1 + 839);
    const auto *sng1_845 = buffer.data(sng1 + 845);
    const auto *sng1_849 = buffer.data(sng1 + 849);
    const auto *sng1_852 = buffer.data(sng1 + 852);
    const auto *sng1_854 = buffer.data(sng1 + 854);
    const auto *sng1_855 = buffer.data(sng1 + 855);
    const auto *sng1_858 = buffer.data(sng1 + 858);
    const auto *sng1_860 = buffer.data(sng1 + 860);
    const auto *sng1_861 = buffer.data(sng1 + 861);
    const auto *sng1_864 = buffer.data(sng1 + 864);
    const auto *sng1_865 = buffer.data(sng1 + 865);

    const auto *snh_1112 = buffer.data(snh + 1112);
    const auto *snh_1113 = buffer.data(snh + 1113);
    const auto *snh_1115 = buffer.data(snh + 1115);
    const auto *snh_1116 = buffer.data(snh + 1116);
    const auto *snh_1118 = buffer.data(snh + 1118);
    const auto *snh_1119 = buffer.data(snh + 1119);
    const auto *snh_1122 = buffer.data(snh + 1122);
    const auto *snh_1128 = buffer.data(snh + 1128);
    const auto *snh_1129 = buffer.data(snh + 1129);
    const auto *snh_1130 = buffer.data(snh + 1130);
    const auto *snh_1131 = buffer.data(snh + 1131);
    const auto *snh_1132 = buffer.data(snh + 1132);
    const auto *snh_1133 = buffer.data(snh + 1133);
    const auto *snh_1134 = buffer.data(snh + 1134);
    const auto *snh_1136 = buffer.data(snh + 1136);
    const auto *snh_1137 = buffer.data(snh + 1137);
    const auto *snh_1139 = buffer.data(snh + 1139);
    const auto *snh_1140 = buffer.data(snh + 1140);
    const auto *snh_1143 = buffer.data(snh + 1143);
    const auto *snh_1149 = buffer.data(snh + 1149);
    const auto *snh_1150 = buffer.data(snh + 1150);
    const auto *snh_1151 = buffer.data(snh + 1151);
    const auto *snh_1152 = buffer.data(snh + 1152);
    const auto *snh_1153 = buffer.data(snh + 1153);
    const auto *snh_1154 = buffer.data(snh + 1154);
    const auto *snh_1155 = buffer.data(snh + 1155);
    const auto *snh_1157 = buffer.data(snh + 1157);
    const auto *snh_1158 = buffer.data(snh + 1158);
    const auto *snh_1160 = buffer.data(snh + 1160);
    const auto *snh_1161 = buffer.data(snh + 1161);
    const auto *snh_1164 = buffer.data(snh + 1164);
    const auto *snh_1165 = buffer.data(snh + 1165);
    const auto *snh_1167 = buffer.data(snh + 1167);
    const auto *snh_1169 = buffer.data(snh + 1169);
    const auto *snh_1170 = buffer.data(snh + 1170);
    const auto *snh_1171 = buffer.data(snh + 1171);
    const auto *snh_1172 = buffer.data(snh + 1172);
    const auto *snh_1173 = buffer.data(snh + 1173);
    const auto *snh_1174 = buffer.data(snh + 1174);
    const auto *snh_1175 = buffer.data(snh + 1175);
    const auto *snh_1176 = buffer.data(snh + 1176);
    const auto *snh_1178 = buffer.data(snh + 1178);
    const auto *snh_1179 = buffer.data(snh + 1179);
    const auto *snh_1181 = buffer.data(snh + 1181);
    const auto *snh_1182 = buffer.data(snh + 1182);
    const auto *snh_1185 = buffer.data(snh + 1185);
    const auto *snh_1188 = buffer.data(snh + 1188);
    const auto *snh_1190 = buffer.data(snh + 1190);
    const auto *snh_1191 = buffer.data(snh + 1191);
    const auto *snh_1192 = buffer.data(snh + 1192);
    const auto *snh_1193 = buffer.data(snh + 1193);
    const auto *snh_1194 = buffer.data(snh + 1194);
    const auto *snh_1195 = buffer.data(snh + 1195);
    const auto *snh_1196 = buffer.data(snh + 1196);
    const auto *snh_1197 = buffer.data(snh + 1197);
    const auto *snh_1199 = buffer.data(snh + 1199);
    const auto *snh_1200 = buffer.data(snh + 1200);
    const auto *snh_1202 = buffer.data(snh + 1202);
    const auto *snh_1203 = buffer.data(snh + 1203);
    const auto *snh_1206 = buffer.data(snh + 1206);
    const auto *snh_1207 = buffer.data(snh + 1207);

#pragma omp simd aligned(t_1482, t_1483, t_1484, t_1485, pb_x, pb_y, pc_x, pc_y, smi0_1232, \
                         smi0_1483, smh_923, smh_924, smi1_1232, smi1_1483, snh_1112, \
                         snh_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1482[k] = f_12 * smh_923[k]
                    + f_3 * pc_y[k] * snh_1112[k];

        t_1483[k] = pb_x[k] * smi0_1483[k]
                    - f_10 * pc_x[k] * smi1_1483[k];

        t_1484[k] = pb_y[k] * smi0_1232[k]
                    - f_10 * pc_y[k] * smi1_1232[k];

        t_1485[k] = f_11 * smh_924[k]
                    + f_3 * pc_y[k] * snh_1113[k];
    }

#pragma omp simd aligned(t_1486, t_1487, t_1488, pb_x, pc_x, pc_y, pc_z, smi0_1487, smh_903, \
                         smh_926, smh_1116, smi1_1487, snh_1113, \
                         snh_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1486[k] = f_16 * smh_903[k]
                    + f_3 * pc_z[k] * snh_1113[k];

        t_1487[k] = pb_x[k] * smi0_1487[k]
                    + f_14 * smh_1116[k]
                    - f_10 * pc_x[k] * smi1_1487[k];

        t_1488[k] = f_11 * smh_926[k]
                    + f_3 * pc_y[k] * snh_1115[k];
    }

#pragma omp simd aligned(t_1489, t_1490, t_1491, pb_x, pb_y, pc_x, pc_y, pc_z, smi0_1237, \
                         smi0_1490, smh_906, smh_1119, smi1_1237, smi1_1490, \
                         snh_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1489[k] = pb_y[k] * smi0_1237[k]
                    - f_10 * pc_y[k] * smi1_1237[k];

        t_1490[k] = pb_x[k] * smi0_1490[k]
                    + f_13 * smh_1119[k]
                    - f_10 * pc_x[k] * smi1_1490[k];

        t_1491[k] = f_16 * smh_906[k]
                    + f_3 * pc_z[k] * snh_1116[k];
    }

#pragma omp simd aligned(t_1492, t_1493, t_1494, pb_x, pb_y, pc_x, pc_y, smi0_1241, smi0_1494, \
                         smh_929, smh_1123, smi1_1241, smi1_1494, \
                         snh_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1492[k] = f_11 * smh_929[k]
                    + f_3 * pc_y[k] * snh_1118[k];

        t_1493[k] = pb_y[k] * smi0_1241[k]
                    - f_10 * pc_y[k] * smi1_1241[k];

        t_1494[k] = pb_x[k] * smi0_1494[k]
                    + f_12 * smh_1123[k]
                    - f_10 * pc_x[k] * smi1_1494[k];
    }

#pragma omp simd aligned(t_1495, t_1496, t_1497, pb_x, pc_x, pc_y, pc_z, smi0_1496, smh_909, \
                         smh_933, smh_1125, smi1_1496, snh_1119, \
                         snh_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1495[k] = f_16 * smh_909[k]
                    + f_3 * pc_z[k] * snh_1119[k];

        t_1496[k] = pb_x[k] * smi0_1496[k]
                    + f_12 * smh_1125[k]
                    - f_10 * pc_x[k] * smi1_1496[k];

        t_1497[k] = f_11 * smh_933[k]
                    + f_3 * pc_y[k] * snh_1122[k];
    }

#pragma omp simd aligned(t_1498, t_1499, t_1500, t_1501, pb_y, pc_x, pc_y, smi0_1246, \
                         smh_1128, smh_1129, smh_1130, smi1_1246, snh_1128, snh_1129, \
                         snh_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1498[k] = pb_y[k] * smi0_1246[k]
                    - f_10 * pc_y[k] * smi1_1246[k];

        t_1499[k] = f_11 * smh_1128[k]
                    + f_3 * pc_x[k] * snh_1128[k];

        t_1500[k] = f_11 * smh_1129[k]
                    + f_3 * pc_x[k] * snh_1129[k];

        t_1501[k] = f_11 * smh_1130[k]
                    + f_3 * pc_x[k] * snh_1130[k];
    }

#pragma omp simd aligned(t_1502, t_1503, t_1504, t_1505, pb_x, pc_x, smi0_1505, smh_1131, \
                         smh_1132, smh_1133, smi1_1505, snh_1131, snh_1132, \
                         snh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = f_11 * smh_1131[k]
                    + f_3 * pc_x[k] * snh_1131[k];

        t_1503[k] = f_11 * smh_1132[k]
                    + f_3 * pc_x[k] * snh_1132[k];

        t_1504[k] = f_11 * smh_1133[k]
                    + f_3 * pc_x[k] * snh_1133[k];

        t_1505[k] = pb_x[k] * smi0_1505[k]
                    - f_10 * pc_x[k] * smi1_1505[k];
    }

#pragma omp simd aligned(t_1506, t_1507, t_1508, t_1509, pb_x, pc_x, pc_z, smi0_1507, \
                         smi0_1508, smi0_1509, smh_918, smi1_1507, smi1_1508, smi1_1509, \
                         snh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1506[k] = f_16 * smh_918[k]
                    + f_3 * pc_z[k] * snh_1128[k];

        t_1507[k] = pb_x[k] * smi0_1507[k]
                    - f_10 * pc_x[k] * smi1_1507[k];

        t_1508[k] = pb_x[k] * smi0_1508[k]
                    - f_10 * pc_x[k] * smi1_1508[k];

        t_1509[k] = pb_x[k] * smi0_1509[k]
                    - f_10 * pc_x[k] * smi1_1509[k];
    }

#pragma omp simd aligned(t_1510, t_1511, t_1512, t_1513, pb_x, pc_x, pc_y, smi0_1511, \
                         smi0_1512, smh_944, smh_1134, smi1_1511, smi1_1512, snh_1133, \
                         snh_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1510[k] = f_11 * smh_944[k]
                    + f_3 * pc_y[k] * snh_1133[k];

        t_1511[k] = pb_x[k] * smi0_1511[k]
                    - f_10 * pc_x[k] * smi1_1511[k];

        t_1512[k] = pb_x[k] * smi0_1512[k]
                    + f_18 * smh_1134[k]
                    - f_10 * pc_x[k] * smi1_1512[k];

        t_1513[k] = f_3 * pc_y[k] * snh_1134[k];
    }

#pragma omp simd aligned(t_1514, t_1515, t_1516, pb_x, pc_x, pc_y, pc_z, smi0_1515, smh_924, \
                         smh_1137, smi1_1515, snh_1134, snh_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1514[k] = f_15 * smh_924[k]
                    + f_3 * pc_z[k] * snh_1134[k];

        t_1515[k] = pb_x[k] * smi0_1515[k]
                    + f_14 * smh_1137[k]
                    - f_10 * pc_x[k] * smi1_1515[k];

        t_1516[k] = f_3 * pc_y[k] * snh_1136[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, pb_x, pc_x, pc_z, smi0_1517, smi0_1518, \
                         smh_927, smh_1139, smh_1140, smi1_1517, smi1_1518, \
                         snh_1137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = pb_x[k] * smi0_1517[k]
                    + f_14 * smh_1139[k]
                    - f_10 * pc_x[k] * smi1_1517[k];

        t_1518[k] = pb_x[k] * smi0_1518[k]
                    + f_13 * smh_1140[k]
                    - f_10 * pc_x[k] * smi1_1518[k];

        t_1519[k] = f_15 * smh_927[k]
                    + f_3 * pc_z[k] * snh_1137[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, pb_x, pc_x, pc_y, smi0_1521, smi0_1522, \
                         smh_1143, smh_1144, smi1_1521, smi1_1522, \
                         snh_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_3 * pc_y[k] * snh_1139[k];

        t_1521[k] = pb_x[k] * smi0_1521[k]
                    + f_13 * smh_1143[k]
                    - f_10 * pc_x[k] * smi1_1521[k];

        t_1522[k] = pb_x[k] * smi0_1522[k]
                    + f_12 * smh_1144[k]
                    - f_10 * pc_x[k] * smi1_1522[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pb_x, pc_x, pc_y, pc_z, smi0_1524, smh_930, \
                         smh_1146, smi1_1524, snh_1140, snh_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_15 * smh_930[k]
                    + f_3 * pc_z[k] * snh_1140[k];

        t_1524[k] = pb_x[k] * smi0_1524[k]
                    + f_12 * smh_1146[k]
                    - f_10 * pc_x[k] * smi1_1524[k];

        t_1525[k] = f_3 * pc_y[k] * snh_1143[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, t_1529, pb_x, pc_x, smi0_1526, smh_1148, \
                         smh_1149, smh_1150, smh_1151, smi1_1526, snh_1149, snh_1150, \
                         snh_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = pb_x[k] * smi0_1526[k]
                    + f_12 * smh_1148[k]
                    - f_10 * pc_x[k] * smi1_1526[k];

        t_1527[k] = f_11 * smh_1149[k]
                    + f_3 * pc_x[k] * snh_1149[k];

        t_1528[k] = f_11 * smh_1150[k]
                    + f_3 * pc_x[k] * snh_1150[k];

        t_1529[k] = f_11 * smh_1151[k]
                    + f_3 * pc_x[k] * snh_1151[k];
    }

#pragma omp simd aligned(t_1530, t_1531, t_1532, t_1533, pb_x, pc_x, smi0_1533, smh_1152, \
                         smh_1153, smh_1154, smi1_1533, snh_1152, snh_1153, \
                         snh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1530[k] = f_11 * smh_1152[k]
                    + f_3 * pc_x[k] * snh_1152[k];

        t_1531[k] = f_11 * smh_1153[k]
                    + f_3 * pc_x[k] * snh_1153[k];

        t_1532[k] = f_11 * smh_1154[k]
                    + f_3 * pc_x[k] * snh_1154[k];

        t_1533[k] = pb_x[k] * smi0_1533[k]
                    - f_10 * pc_x[k] * smi1_1533[k];
    }

#pragma omp simd aligned(t_1534, t_1535, t_1536, t_1537, pb_x, pc_x, pc_z, smi0_1535, \
                         smi0_1536, smi0_1537, smh_939, smi1_1535, smi1_1536, smi1_1537, \
                         snh_1149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1534[k] = f_15 * smh_939[k]
                    + f_3 * pc_z[k] * snh_1149[k];

        t_1535[k] = pb_x[k] * smi0_1535[k]
                    - f_10 * pc_x[k] * smi1_1535[k];

        t_1536[k] = pb_x[k] * smi0_1536[k]
                    - f_10 * pc_x[k] * smi1_1536[k];

        t_1537[k] = pb_x[k] * smi0_1537[k]
                    - f_10 * pc_x[k] * smi1_1537[k];
    }

#pragma omp simd aligned(t_1538, t_1539, t_1540, t_1541, t_1542, pb_x, pc_x, pc_y, pc_z, \
                         smi0_1539, smh_945, smi1_1539, sng0_825, sng1_825, snh_1154, \
                         snh_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1538[k] = f_3 * pc_y[k] * snh_1154[k];

        t_1539[k] = pb_x[k] * smi0_1539[k]
                    - f_10 * pc_x[k] * smi1_1539[k];

        t_1540[k] = f_1 * sng0_825[k]
                    - f_2 * sng1_825[k]
                    + f_3 * pc_x[k] * snh_1155[k];

        t_1541[k] = f_0 * smh_945[k]
                    + f_3 * pc_y[k] * snh_1155[k];

        t_1542[k] = f_3 * pc_z[k] * snh_1155[k];
    }

#pragma omp simd aligned(t_1543, t_1544, t_1545, pc_x, pc_y, smh_947, sng0_828, sng0_830, \
                         sng1_828, sng1_830, snh_1157, snh_1158, \
                         snh_1160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1543[k] = f_4 * sng0_828[k]
                    - f_5 * sng1_828[k]
                    + f_3 * pc_x[k] * snh_1158[k];

        t_1544[k] = f_0 * smh_947[k]
                    + f_3 * pc_y[k] * snh_1157[k];

        t_1545[k] = f_4 * sng0_830[k]
                    - f_5 * sng1_830[k]
                    + f_3 * pc_x[k] * snh_1160[k];
    }

#pragma omp simd aligned(t_1546, t_1547, t_1548, t_1549, pc_x, pc_y, pc_z, smh_950, sng0_831, \
                         sng0_834, sng1_831, sng1_834, snh_1158, snh_1160, snh_1161, \
                         snh_1164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1546[k] = f_6 * sng0_831[k]
                    - f_7 * sng1_831[k]
                    + f_3 * pc_x[k] * snh_1161[k];

        t_1547[k] = f_3 * pc_z[k] * snh_1158[k];

        t_1548[k] = f_0 * smh_950[k]
                    + f_3 * pc_y[k] * snh_1160[k];

        t_1549[k] = f_6 * sng0_834[k]
                    - f_7 * sng1_834[k]
                    + f_3 * pc_x[k] * snh_1164[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, t_1553, pc_x, pc_y, pc_z, smh_954, sng0_835, \
                         sng0_837, sng1_835, sng1_837, snh_1161, snh_1164, snh_1165, \
                         snh_1167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = f_8 * sng0_835[k]
                    - f_9 * sng1_835[k]
                    + f_3 * pc_x[k] * snh_1165[k];

        t_1551[k] = f_3 * pc_z[k] * snh_1161[k];

        t_1552[k] = f_8 * sng0_837[k]
                    - f_9 * sng1_837[k]
                    + f_3 * pc_x[k] * snh_1167[k];

        t_1553[k] = f_0 * smh_954[k]
                    + f_3 * pc_y[k] * snh_1164[k];
    }

#pragma omp simd aligned(t_1554, t_1555, t_1556, t_1557, t_1558, t_1559, pc_x, sng0_839, \
                         sng1_839, snh_1169, snh_1170, snh_1171, snh_1172, snh_1173, \
                         snh_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1554[k] = f_8 * sng0_839[k]
                    - f_9 * sng1_839[k]
                    + f_3 * pc_x[k] * snh_1169[k];

        t_1555[k] = f_3 * pc_x[k] * snh_1170[k];

        t_1556[k] = f_3 * pc_x[k] * snh_1171[k];

        t_1557[k] = f_3 * pc_x[k] * snh_1172[k];

        t_1558[k] = f_3 * pc_x[k] * snh_1173[k];

        t_1559[k] = f_3 * pc_x[k] * snh_1174[k];
    }

#pragma omp simd aligned(t_1560, t_1561, t_1562, t_1563, pc_x, pc_y, pc_z, smh_960, smh_962, \
                         sng0_835, sng0_837, sng1_835, sng1_837, snh_1170, snh_1172, \
                         snh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1560[k] = f_3 * pc_x[k] * snh_1175[k];

        t_1561[k] = f_0 * smh_960[k]
                    + f_1 * sng0_835[k]
                    - f_2 * sng1_835[k]
                    + f_3 * pc_y[k] * snh_1170[k];

        t_1562[k] = f_3 * pc_z[k] * snh_1170[k];

        t_1563[k] = f_0 * smh_962[k]
                    + f_4 * sng0_837[k]
                    - f_5 * sng1_837[k]
                    + f_3 * pc_y[k] * snh_1172[k];
    }

#pragma omp simd aligned(t_1564, t_1565, t_1566, t_1567, pc_y, pc_z, smh_963, smh_964, \
                         smh_965, sng0_838, sng0_839, sng1_838, sng1_839, snh_1173, snh_1174, \
                         snh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1564[k] = f_0 * smh_963[k]
                    + f_6 * sng0_838[k]
                    - f_7 * sng1_838[k]
                    + f_3 * pc_y[k] * snh_1173[k];

        t_1565[k] = f_0 * smh_964[k]
                    + f_8 * sng0_839[k]
                    - f_9 * sng1_839[k]
                    + f_3 * pc_y[k] * snh_1174[k];

        t_1566[k] = f_0 * smh_965[k]
                    + f_3 * pc_y[k] * snh_1175[k];

        t_1567[k] = f_1 * sng0_839[k]
                    - f_2 * sng1_839[k]
                    + f_3 * pc_z[k] * snh_1175[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, t_1571, pb_z, pc_y, pc_z, smi0_1260, \
                         smi0_1263, smh_945, smh_966, smi1_1260, smi1_1263, \
                         snh_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = pb_z[k] * smi0_1260[k]
                    - f_10 * pc_z[k] * smi1_1260[k];

        t_1569[k] = f_15 * smh_966[k]
                    + f_3 * pc_y[k] * snh_1176[k];

        t_1570[k] = f_11 * smh_945[k]
                    + f_3 * pc_z[k] * snh_1176[k];

        t_1571[k] = pb_z[k] * smi0_1263[k]
                    - f_10 * pc_z[k] * smi1_1263[k];
    }

#pragma omp simd aligned(t_1572, t_1573, t_1574, pb_z, pc_x, pc_y, pc_z, smi0_1266, smh_968, \
                         smi1_1266, sng0_845, sng1_845, snh_1178, \
                         snh_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1572[k] = f_15 * smh_968[k]
                    + f_3 * pc_y[k] * snh_1178[k];

        t_1573[k] = f_4 * sng0_845[k]
                    - f_5 * sng1_845[k]
                    + f_3 * pc_x[k] * snh_1181[k];

        t_1574[k] = pb_z[k] * smi0_1266[k]
                    - f_10 * pc_z[k] * smi1_1266[k];
    }

#pragma omp simd aligned(t_1575, t_1576, t_1577, pc_x, pc_y, pc_z, smh_948, smh_971, sng0_849, \
                         sng1_849, snh_1179, snh_1181, snh_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1575[k] = f_11 * smh_948[k]
                    + f_3 * pc_z[k] * snh_1179[k];

        t_1576[k] = f_15 * smh_971[k]
                    + f_3 * pc_y[k] * snh_1181[k];

        t_1577[k] = f_6 * sng0_849[k]
                    - f_7 * sng1_849[k]
                    + f_3 * pc_x[k] * snh_1185[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pb_z, pc_x, pc_z, smi0_1270, smh_951, \
                         smi1_1270, sng0_852, sng1_852, snh_1182, \
                         snh_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = pb_z[k] * smi0_1270[k]
                    - f_10 * pc_z[k] * smi1_1270[k];

        t_1579[k] = f_11 * smh_951[k]
                    + f_3 * pc_z[k] * snh_1182[k];

        t_1580[k] = f_8 * sng0_852[k]
                    - f_9 * sng1_852[k]
                    + f_3 * pc_x[k] * snh_1188[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, t_1584, t_1585, pc_x, pc_y, smh_975, \
                         sng0_854, sng1_854, snh_1185, snh_1190, snh_1191, snh_1192, \
                         snh_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_15 * smh_975[k]
                    + f_3 * pc_y[k] * snh_1185[k];

        t_1582[k] = f_8 * sng0_854[k]
                    - f_9 * sng1_854[k]
                    + f_3 * pc_x[k] * snh_1190[k];

        t_1583[k] = f_3 * pc_x[k] * snh_1191[k];

        t_1584[k] = f_3 * pc_x[k] * snh_1192[k];

        t_1585[k] = f_3 * pc_x[k] * snh_1193[k];
    }

#pragma omp simd aligned(t_1586, t_1587, t_1588, t_1589, t_1590, pb_z, pc_x, pc_z, smi0_1281, \
                         smh_960, smi1_1281, snh_1191, snh_1194, snh_1195, \
                         snh_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1586[k] = f_3 * pc_x[k] * snh_1194[k];

        t_1587[k] = f_3 * pc_x[k] * snh_1195[k];

        t_1588[k] = f_3 * pc_x[k] * snh_1196[k];

        t_1589[k] = pb_z[k] * smi0_1281[k]
                    - f_10 * pc_z[k] * smi1_1281[k];

        t_1590[k] = f_11 * smh_960[k]
                    + f_3 * pc_z[k] * snh_1191[k];
    }

#pragma omp simd aligned(t_1591, t_1592, t_1593, pb_z, pc_z, smi0_1283, smi0_1284, smi0_1285, \
                         smh_961, smh_962, smh_963, smi1_1283, smi1_1284, \
                         smi1_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1591[k] = pb_z[k] * smi0_1283[k]
                    + f_12 * smh_961[k]
                    - f_10 * pc_z[k] * smi1_1283[k];

        t_1592[k] = pb_z[k] * smi0_1284[k]
                    + f_13 * smh_962[k]
                    - f_10 * pc_z[k] * smi1_1284[k];

        t_1593[k] = pb_z[k] * smi0_1285[k]
                    + f_14 * smh_963[k]
                    - f_10 * pc_z[k] * smi1_1285[k];
    }

#pragma omp simd aligned(t_1594, t_1595, t_1596, t_1597, pc_x, pc_y, pc_z, smh_965, smh_986, \
                         smh_987, sng0_854, sng0_855, sng1_854, sng1_855, snh_1196, \
                         snh_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1594[k] = f_15 * smh_986[k]
                    + f_3 * pc_y[k] * snh_1196[k];

        t_1595[k] = f_11 * smh_965[k]
                    + f_1 * sng0_854[k]
                    - f_2 * sng1_854[k]
                    + f_3 * pc_z[k] * snh_1196[k];

        t_1596[k] = f_1 * sng0_855[k]
                    - f_2 * sng1_855[k]
                    + f_3 * pc_x[k] * snh_1197[k];

        t_1597[k] = f_16 * smh_987[k]
                    + f_3 * pc_y[k] * snh_1197[k];
    }

#pragma omp simd aligned(t_1598, t_1599, t_1600, pc_x, pc_y, pc_z, smh_966, smh_989, sng0_858, \
                         sng1_858, snh_1197, snh_1199, snh_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1598[k] = f_12 * smh_966[k]
                    + f_3 * pc_z[k] * snh_1197[k];

        t_1599[k] = f_4 * sng0_858[k]
                    - f_5 * sng1_858[k]
                    + f_3 * pc_x[k] * snh_1200[k];

        t_1600[k] = f_16 * smh_989[k]
                    + f_3 * pc_y[k] * snh_1199[k];
    }

#pragma omp simd aligned(t_1601, t_1602, t_1603, t_1604, pc_x, pc_y, pc_z, smh_969, smh_992, \
                         sng0_860, sng0_861, sng1_860, sng1_861, snh_1200, snh_1202, \
                         snh_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1601[k] = f_4 * sng0_860[k]
                    - f_5 * sng1_860[k]
                    + f_3 * pc_x[k] * snh_1202[k];

        t_1602[k] = f_6 * sng0_861[k]
                    - f_7 * sng1_861[k]
                    + f_3 * pc_x[k] * snh_1203[k];

        t_1603[k] = f_12 * smh_969[k]
                    + f_3 * pc_z[k] * snh_1200[k];

        t_1604[k] = f_16 * smh_992[k]
                    + f_3 * pc_y[k] * snh_1202[k];
    }

#pragma omp simd aligned(t_1605, t_1606, t_1607, pc_x, pc_z, smh_972, sng0_864, sng0_865, \
                         sng1_864, sng1_865, snh_1203, snh_1206, \
                         snh_1207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1605[k] = f_6 * sng0_864[k]
                    - f_7 * sng1_864[k]
                    + f_3 * pc_x[k] * snh_1206[k];

        t_1606[k] = f_8 * sng0_865[k]
                    - f_9 * sng1_865[k]
                    + f_3 * pc_x[k] * snh_1207[k];

        t_1607[k] = f_12 * smh_972[k]
                    + f_3 * pc_z[k] * snh_1203[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t smh, const size_t sng0,
                                                           const size_t sng1, const size_t snh,
                                                           const size_t ncols,
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_1608 = buffer.data(target + 1608);
    auto *t_1609 = buffer.data(target + 1609);
    auto *t_1610 = buffer.data(target + 1610);
    auto *t_1611 = buffer.data(target + 1611);
    auto *t_1612 = buffer.data(target + 1612);
    auto *t_1613 = buffer.data(target + 1613);
    auto *t_1614 = buffer.data(target + 1614);
    auto *t_1615 = buffer.data(target + 1615);
    auto *t_1616 = buffer.data(target + 1616);
    auto *t_1617 = buffer.data(target + 1617);
    auto *t_1618 = buffer.data(target + 1618);
    auto *t_1619 = buffer.data(target + 1619);
    auto *t_1620 = buffer.data(target + 1620);
    auto *t_1621 = buffer.data(target + 1621);
    auto *t_1622 = buffer.data(target + 1622);
    auto *t_1623 = buffer.data(target + 1623);
    auto *t_1624 = buffer.data(target + 1624);
    auto *t_1625 = buffer.data(target + 1625);
    auto *t_1626 = buffer.data(target + 1626);
    auto *t_1627 = buffer.data(target + 1627);
    auto *t_1628 = buffer.data(target + 1628);
    auto *t_1629 = buffer.data(target + 1629);
    auto *t_1630 = buffer.data(target + 1630);
    auto *t_1631 = buffer.data(target + 1631);
    auto *t_1632 = buffer.data(target + 1632);
    auto *t_1633 = buffer.data(target + 1633);
    auto *t_1634 = buffer.data(target + 1634);
    auto *t_1635 = buffer.data(target + 1635);
    auto *t_1636 = buffer.data(target + 1636);
    auto *t_1637 = buffer.data(target + 1637);
    auto *t_1638 = buffer.data(target + 1638);
    auto *t_1639 = buffer.data(target + 1639);
    auto *t_1640 = buffer.data(target + 1640);
    auto *t_1641 = buffer.data(target + 1641);
    auto *t_1642 = buffer.data(target + 1642);
    auto *t_1643 = buffer.data(target + 1643);
    auto *t_1644 = buffer.data(target + 1644);
    auto *t_1645 = buffer.data(target + 1645);
    auto *t_1646 = buffer.data(target + 1646);
    auto *t_1647 = buffer.data(target + 1647);
    auto *t_1648 = buffer.data(target + 1648);
    auto *t_1649 = buffer.data(target + 1649);
    auto *t_1650 = buffer.data(target + 1650);
    auto *t_1651 = buffer.data(target + 1651);
    auto *t_1652 = buffer.data(target + 1652);
    auto *t_1653 = buffer.data(target + 1653);
    auto *t_1654 = buffer.data(target + 1654);
    auto *t_1655 = buffer.data(target + 1655);
    auto *t_1656 = buffer.data(target + 1656);
    auto *t_1657 = buffer.data(target + 1657);
    auto *t_1658 = buffer.data(target + 1658);
    auto *t_1659 = buffer.data(target + 1659);
    auto *t_1660 = buffer.data(target + 1660);
    auto *t_1661 = buffer.data(target + 1661);
    auto *t_1662 = buffer.data(target + 1662);
    auto *t_1663 = buffer.data(target + 1663);
    auto *t_1664 = buffer.data(target + 1664);
    auto *t_1665 = buffer.data(target + 1665);
    auto *t_1666 = buffer.data(target + 1666);
    auto *t_1667 = buffer.data(target + 1667);
    auto *t_1668 = buffer.data(target + 1668);
    auto *t_1669 = buffer.data(target + 1669);
    auto *t_1670 = buffer.data(target + 1670);
    auto *t_1671 = buffer.data(target + 1671);
    auto *t_1672 = buffer.data(target + 1672);
    auto *t_1673 = buffer.data(target + 1673);
    auto *t_1674 = buffer.data(target + 1674);
    auto *t_1675 = buffer.data(target + 1675);
    auto *t_1676 = buffer.data(target + 1676);
    auto *t_1677 = buffer.data(target + 1677);
    auto *t_1678 = buffer.data(target + 1678);
    auto *t_1679 = buffer.data(target + 1679);
    auto *t_1680 = buffer.data(target + 1680);
    auto *t_1681 = buffer.data(target + 1681);
    auto *t_1682 = buffer.data(target + 1682);
    auto *t_1683 = buffer.data(target + 1683);
    auto *t_1684 = buffer.data(target + 1684);
    auto *t_1685 = buffer.data(target + 1685);
    auto *t_1686 = buffer.data(target + 1686);
    auto *t_1687 = buffer.data(target + 1687);
    auto *t_1688 = buffer.data(target + 1688);
    auto *t_1689 = buffer.data(target + 1689);
    auto *t_1690 = buffer.data(target + 1690);
    auto *t_1691 = buffer.data(target + 1691);
    auto *t_1692 = buffer.data(target + 1692);
    auto *t_1693 = buffer.data(target + 1693);
    auto *t_1694 = buffer.data(target + 1694);
    auto *t_1695 = buffer.data(target + 1695);
    auto *t_1696 = buffer.data(target + 1696);
    auto *t_1697 = buffer.data(target + 1697);
    auto *t_1698 = buffer.data(target + 1698);
    auto *t_1699 = buffer.data(target + 1699);
    auto *t_1700 = buffer.data(target + 1700);
    auto *t_1701 = buffer.data(target + 1701);
    auto *t_1702 = buffer.data(target + 1702);
    auto *t_1703 = buffer.data(target + 1703);
    auto *t_1704 = buffer.data(target + 1704);
    auto *t_1705 = buffer.data(target + 1705);
    auto *t_1706 = buffer.data(target + 1706);
    auto *t_1707 = buffer.data(target + 1707);
    auto *t_1708 = buffer.data(target + 1708);
    auto *t_1709 = buffer.data(target + 1709);
    auto *t_1710 = buffer.data(target + 1710);
    auto *t_1711 = buffer.data(target + 1711);
    auto *t_1712 = buffer.data(target + 1712);
    auto *t_1713 = buffer.data(target + 1713);
    auto *t_1714 = buffer.data(target + 1714);
    auto *t_1715 = buffer.data(target + 1715);
    auto *t_1716 = buffer.data(target + 1716);
    auto *t_1717 = buffer.data(target + 1717);
    auto *t_1718 = buffer.data(target + 1718);
    auto *t_1719 = buffer.data(target + 1719);
    auto *t_1720 = buffer.data(target + 1720);
    auto *t_1721 = buffer.data(target + 1721);
    auto *t_1722 = buffer.data(target + 1722);
    auto *t_1723 = buffer.data(target + 1723);
    auto *t_1724 = buffer.data(target + 1724);
    auto *t_1725 = buffer.data(target + 1725);
    auto *t_1726 = buffer.data(target + 1726);
    auto *t_1727 = buffer.data(target + 1727);
    auto *t_1728 = buffer.data(target + 1728);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smh_981 = buffer.data(smh + 981);
    const auto *smh_986 = buffer.data(smh + 986);
    const auto *smh_987 = buffer.data(smh + 987);
    const auto *smh_990 = buffer.data(smh + 990);
    const auto *smh_993 = buffer.data(smh + 993);
    const auto *smh_996 = buffer.data(smh + 996);
    const auto *smh_1002 = buffer.data(smh + 1002);
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
    const auto *smh_1023 = buffer.data(smh + 1023);
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
    const auto *smh_1044 = buffer.data(smh + 1044);
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
    const auto *smh_1065 = buffer.data(smh + 1065);
    const auto *smh_1067 = buffer.data(smh + 1067);
    const auto *smh_1068 = buffer.data(smh + 1068);
    const auto *smh_1069 = buffer.data(smh + 1069);
    const auto *smh_1070 = buffer.data(smh + 1070);
    const auto *smh_1071 = buffer.data(smh + 1071);
    const auto *smh_1073 = buffer.data(smh + 1073);
    const auto *smh_1076 = buffer.data(smh + 1076);
    const auto *smh_1080 = buffer.data(smh + 1080);

    const auto *sng0_865 = buffer.data(sng0 + 865);
    const auto *sng0_867 = buffer.data(sng0 + 867);
    const auto *sng0_868 = buffer.data(sng0 + 868);
    const auto *sng0_869 = buffer.data(sng0 + 869);
    const auto *sng0_870 = buffer.data(sng0 + 870);
    const auto *sng0_873 = buffer.data(sng0 + 873);
    const auto *sng0_875 = buffer.data(sng0 + 875);
    const auto *sng0_876 = buffer.data(sng0 + 876);
    const auto *sng0_879 = buffer.data(sng0 + 879);
    const auto *sng0_880 = buffer.data(sng0 + 880);
    const auto *sng0_882 = buffer.data(sng0 + 882);
    const auto *sng0_883 = buffer.data(sng0 + 883);
    const auto *sng0_884 = buffer.data(sng0 + 884);
    const auto *sng0_885 = buffer.data(sng0 + 885);
    const auto *sng0_888 = buffer.data(sng0 + 888);
    const auto *sng0_890 = buffer.data(sng0 + 890);
    const auto *sng0_891 = buffer.data(sng0 + 891);
    const auto *sng0_894 = buffer.data(sng0 + 894);
    const auto *sng0_895 = buffer.data(sng0 + 895);
    const auto *sng0_897 = buffer.data(sng0 + 897);
    const auto *sng0_898 = buffer.data(sng0 + 898);
    const auto *sng0_899 = buffer.data(sng0 + 899);
    const auto *sng0_900 = buffer.data(sng0 + 900);
    const auto *sng0_903 = buffer.data(sng0 + 903);
    const auto *sng0_905 = buffer.data(sng0 + 905);
    const auto *sng0_906 = buffer.data(sng0 + 906);
    const auto *sng0_909 = buffer.data(sng0 + 909);
    const auto *sng0_910 = buffer.data(sng0 + 910);
    const auto *sng0_912 = buffer.data(sng0 + 912);
    const auto *sng0_913 = buffer.data(sng0 + 913);
    const auto *sng0_914 = buffer.data(sng0 + 914);
    const auto *sng0_915 = buffer.data(sng0 + 915);
    const auto *sng0_918 = buffer.data(sng0 + 918);
    const auto *sng0_920 = buffer.data(sng0 + 920);
    const auto *sng0_921 = buffer.data(sng0 + 921);
    const auto *sng0_924 = buffer.data(sng0 + 924);
    const auto *sng0_925 = buffer.data(sng0 + 925);
    const auto *sng0_927 = buffer.data(sng0 + 927);
    const auto *sng0_929 = buffer.data(sng0 + 929);

    const auto *sng1_865 = buffer.data(sng1 + 865);
    const auto *sng1_867 = buffer.data(sng1 + 867);
    const auto *sng1_868 = buffer.data(sng1 + 868);
    const auto *sng1_869 = buffer.data(sng1 + 869);
    const auto *sng1_870 = buffer.data(sng1 + 870);
    const auto *sng1_873 = buffer.data(sng1 + 873);
    const auto *sng1_875 = buffer.data(sng1 + 875);
    const auto *sng1_876 = buffer.data(sng1 + 876);
    const auto *sng1_879 = buffer.data(sng1 + 879);
    const auto *sng1_880 = buffer.data(sng1 + 880);
    const auto *sng1_882 = buffer.data(sng1 + 882);
    const auto *sng1_883 = buffer.data(sng1 + 883);
    const auto *sng1_884 = buffer.data(sng1 + 884);
    const auto *sng1_885 = buffer.data(sng1 + 885);
    const auto *sng1_888 = buffer.data(sng1 + 888);
    const auto *sng1_890 = buffer.data(sng1 + 890);
    const auto *sng1_891 = buffer.data(sng1 + 891);
    const auto *sng1_894 = buffer.data(sng1 + 894);
    const auto *sng1_895 = buffer.data(sng1 + 895);
    const auto *sng1_897 = buffer.data(sng1 + 897);
    const auto *sng1_898 = buffer.data(sng1 + 898);
    const auto *sng1_899 = buffer.data(sng1 + 899);
    const auto *sng1_900 = buffer.data(sng1 + 900);
    const auto *sng1_903 = buffer.data(sng1 + 903);
    const auto *sng1_905 = buffer.data(sng1 + 905);
    const auto *sng1_906 = buffer.data(sng1 + 906);
    const auto *sng1_909 = buffer.data(sng1 + 909);
    const auto *sng1_910 = buffer.data(sng1 + 910);
    const auto *sng1_912 = buffer.data(sng1 + 912);
    const auto *sng1_913 = buffer.data(sng1 + 913);
    const auto *sng1_914 = buffer.data(sng1 + 914);
    const auto *sng1_915 = buffer.data(sng1 + 915);
    const auto *sng1_918 = buffer.data(sng1 + 918);
    const auto *sng1_920 = buffer.data(sng1 + 920);
    const auto *sng1_921 = buffer.data(sng1 + 921);
    const auto *sng1_924 = buffer.data(sng1 + 924);
    const auto *sng1_925 = buffer.data(sng1 + 925);
    const auto *sng1_927 = buffer.data(sng1 + 927);
    const auto *sng1_929 = buffer.data(sng1 + 929);

    const auto *snh_1206 = buffer.data(snh + 1206);
    const auto *snh_1209 = buffer.data(snh + 1209);
    const auto *snh_1211 = buffer.data(snh + 1211);
    const auto *snh_1212 = buffer.data(snh + 1212);
    const auto *snh_1213 = buffer.data(snh + 1213);
    const auto *snh_1214 = buffer.data(snh + 1214);
    const auto *snh_1215 = buffer.data(snh + 1215);
    const auto *snh_1216 = buffer.data(snh + 1216);
    const auto *snh_1217 = buffer.data(snh + 1217);
    const auto *snh_1218 = buffer.data(snh + 1218);
    const auto *snh_1220 = buffer.data(snh + 1220);
    const auto *snh_1221 = buffer.data(snh + 1221);
    const auto *snh_1223 = buffer.data(snh + 1223);
    const auto *snh_1224 = buffer.data(snh + 1224);
    const auto *snh_1227 = buffer.data(snh + 1227);
    const auto *snh_1228 = buffer.data(snh + 1228);
    const auto *snh_1230 = buffer.data(snh + 1230);
    const auto *snh_1232 = buffer.data(snh + 1232);
    const auto *snh_1233 = buffer.data(snh + 1233);
    const auto *snh_1234 = buffer.data(snh + 1234);
    const auto *snh_1235 = buffer.data(snh + 1235);
    const auto *snh_1236 = buffer.data(snh + 1236);
    const auto *snh_1237 = buffer.data(snh + 1237);
    const auto *snh_1238 = buffer.data(snh + 1238);
    const auto *snh_1239 = buffer.data(snh + 1239);
    const auto *snh_1241 = buffer.data(snh + 1241);
    const auto *snh_1242 = buffer.data(snh + 1242);
    const auto *snh_1244 = buffer.data(snh + 1244);
    const auto *snh_1245 = buffer.data(snh + 1245);
    const auto *snh_1248 = buffer.data(snh + 1248);
    const auto *snh_1249 = buffer.data(snh + 1249);
    const auto *snh_1251 = buffer.data(snh + 1251);
    const auto *snh_1253 = buffer.data(snh + 1253);
    const auto *snh_1254 = buffer.data(snh + 1254);
    const auto *snh_1255 = buffer.data(snh + 1255);
    const auto *snh_1256 = buffer.data(snh + 1256);
    const auto *snh_1257 = buffer.data(snh + 1257);
    const auto *snh_1258 = buffer.data(snh + 1258);
    const auto *snh_1259 = buffer.data(snh + 1259);
    const auto *snh_1260 = buffer.data(snh + 1260);
    const auto *snh_1262 = buffer.data(snh + 1262);
    const auto *snh_1263 = buffer.data(snh + 1263);
    const auto *snh_1265 = buffer.data(snh + 1265);
    const auto *snh_1266 = buffer.data(snh + 1266);
    const auto *snh_1269 = buffer.data(snh + 1269);
    const auto *snh_1270 = buffer.data(snh + 1270);
    const auto *snh_1272 = buffer.data(snh + 1272);
    const auto *snh_1274 = buffer.data(snh + 1274);
    const auto *snh_1275 = buffer.data(snh + 1275);
    const auto *snh_1276 = buffer.data(snh + 1276);
    const auto *snh_1277 = buffer.data(snh + 1277);
    const auto *snh_1278 = buffer.data(snh + 1278);
    const auto *snh_1279 = buffer.data(snh + 1279);
    const auto *snh_1280 = buffer.data(snh + 1280);
    const auto *snh_1281 = buffer.data(snh + 1281);
    const auto *snh_1283 = buffer.data(snh + 1283);
    const auto *snh_1284 = buffer.data(snh + 1284);
    const auto *snh_1286 = buffer.data(snh + 1286);
    const auto *snh_1287 = buffer.data(snh + 1287);
    const auto *snh_1290 = buffer.data(snh + 1290);
    const auto *snh_1291 = buffer.data(snh + 1291);
    const auto *snh_1293 = buffer.data(snh + 1293);
    const auto *snh_1295 = buffer.data(snh + 1295);
    const auto *snh_1296 = buffer.data(snh + 1296);
    const auto *snh_1297 = buffer.data(snh + 1297);
    const auto *snh_1298 = buffer.data(snh + 1298);
    const auto *snh_1299 = buffer.data(snh + 1299);
    const auto *snh_1300 = buffer.data(snh + 1300);
    const auto *snh_1301 = buffer.data(snh + 1301);

#pragma omp simd aligned(t_1608, t_1609, t_1610, t_1611, pc_x, pc_y, smh_996, sng0_867, \
                         sng0_869, sng1_867, sng1_869, snh_1206, snh_1209, snh_1211, \
                         snh_1212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1608[k] = f_8 * sng0_867[k]
                    - f_9 * sng1_867[k]
                    + f_3 * pc_x[k] * snh_1209[k];

        t_1609[k] = f_16 * smh_996[k]
                    + f_3 * pc_y[k] * snh_1206[k];

        t_1610[k] = f_8 * sng0_869[k]
                    - f_9 * sng1_869[k]
                    + f_3 * pc_x[k] * snh_1211[k];

        t_1611[k] = f_3 * pc_x[k] * snh_1212[k];
    }

#pragma omp simd aligned(t_1612, t_1613, t_1614, t_1615, t_1616, pc_x, snh_1213, snh_1214, \
                         snh_1215, snh_1216, snh_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1612[k] = f_3 * pc_x[k] * snh_1213[k];

        t_1613[k] = f_3 * pc_x[k] * snh_1214[k];

        t_1614[k] = f_3 * pc_x[k] * snh_1215[k];

        t_1615[k] = f_3 * pc_x[k] * snh_1216[k];

        t_1616[k] = f_3 * pc_x[k] * snh_1217[k];
    }

#pragma omp simd aligned(t_1617, t_1618, t_1619, pc_y, pc_z, smh_981, smh_1002, smh_1004, \
                         sng0_865, sng0_867, sng1_865, sng1_867, snh_1212, \
                         snh_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = f_16 * smh_1002[k]
                    + f_1 * sng0_865[k]
                    - f_2 * sng1_865[k]
                    + f_3 * pc_y[k] * snh_1212[k];

        t_1618[k] = f_12 * smh_981[k]
                    + f_3 * pc_z[k] * snh_1212[k];

        t_1619[k] = f_16 * smh_1004[k]
                    + f_4 * sng0_867[k]
                    - f_5 * sng1_867[k]
                    + f_3 * pc_y[k] * snh_1214[k];
    }

#pragma omp simd aligned(t_1620, t_1621, t_1622, pc_y, smh_1005, smh_1006, smh_1007, sng0_868, \
                         sng0_869, sng1_868, sng1_869, snh_1215, snh_1216, \
                         snh_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1620[k] = f_16 * smh_1005[k]
                    + f_6 * sng0_868[k]
                    - f_7 * sng1_868[k]
                    + f_3 * pc_y[k] * snh_1215[k];

        t_1621[k] = f_16 * smh_1006[k]
                    + f_8 * sng0_869[k]
                    - f_9 * sng1_869[k]
                    + f_3 * pc_y[k] * snh_1216[k];

        t_1622[k] = f_16 * smh_1007[k]
                    + f_3 * pc_y[k] * snh_1217[k];
    }

#pragma omp simd aligned(t_1623, t_1624, t_1625, t_1626, pc_x, pc_y, pc_z, smh_986, smh_987, \
                         smh_1008, sng0_869, sng0_870, sng1_869, sng1_870, snh_1217, \
                         snh_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1623[k] = f_12 * smh_986[k]
                    + f_1 * sng0_869[k]
                    - f_2 * sng1_869[k]
                    + f_3 * pc_z[k] * snh_1217[k];

        t_1624[k] = f_1 * sng0_870[k]
                    - f_2 * sng1_870[k]
                    + f_3 * pc_x[k] * snh_1218[k];

        t_1625[k] = f_17 * smh_1008[k]
                    + f_3 * pc_y[k] * snh_1218[k];

        t_1626[k] = f_13 * smh_987[k]
                    + f_3 * pc_z[k] * snh_1218[k];
    }

#pragma omp simd aligned(t_1627, t_1628, t_1629, pc_x, pc_y, smh_1010, sng0_873, sng0_875, \
                         sng1_873, sng1_875, snh_1220, snh_1221, \
                         snh_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1627[k] = f_4 * sng0_873[k]
                    - f_5 * sng1_873[k]
                    + f_3 * pc_x[k] * snh_1221[k];

        t_1628[k] = f_17 * smh_1010[k]
                    + f_3 * pc_y[k] * snh_1220[k];

        t_1629[k] = f_4 * sng0_875[k]
                    - f_5 * sng1_875[k]
                    + f_3 * pc_x[k] * snh_1223[k];
    }

#pragma omp simd aligned(t_1630, t_1631, t_1632, pc_x, pc_y, pc_z, smh_990, smh_1013, \
                         sng0_876, sng1_876, snh_1221, snh_1223, \
                         snh_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1630[k] = f_6 * sng0_876[k]
                    - f_7 * sng1_876[k]
                    + f_3 * pc_x[k] * snh_1224[k];

        t_1631[k] = f_13 * smh_990[k]
                    + f_3 * pc_z[k] * snh_1221[k];

        t_1632[k] = f_17 * smh_1013[k]
                    + f_3 * pc_y[k] * snh_1223[k];
    }

#pragma omp simd aligned(t_1633, t_1634, t_1635, pc_x, pc_z, smh_993, sng0_879, sng0_880, \
                         sng1_879, sng1_880, snh_1224, snh_1227, \
                         snh_1228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1633[k] = f_6 * sng0_879[k]
                    - f_7 * sng1_879[k]
                    + f_3 * pc_x[k] * snh_1227[k];

        t_1634[k] = f_8 * sng0_880[k]
                    - f_9 * sng1_880[k]
                    + f_3 * pc_x[k] * snh_1228[k];

        t_1635[k] = f_13 * smh_993[k]
                    + f_3 * pc_z[k] * snh_1224[k];
    }

#pragma omp simd aligned(t_1636, t_1637, t_1638, t_1639, pc_x, pc_y, smh_1017, sng0_882, \
                         sng0_884, sng1_882, sng1_884, snh_1227, snh_1230, snh_1232, \
                         snh_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1636[k] = f_8 * sng0_882[k]
                    - f_9 * sng1_882[k]
                    + f_3 * pc_x[k] * snh_1230[k];

        t_1637[k] = f_17 * smh_1017[k]
                    + f_3 * pc_y[k] * snh_1227[k];

        t_1638[k] = f_8 * sng0_884[k]
                    - f_9 * sng1_884[k]
                    + f_3 * pc_x[k] * snh_1232[k];

        t_1639[k] = f_3 * pc_x[k] * snh_1233[k];
    }

#pragma omp simd aligned(t_1640, t_1641, t_1642, t_1643, t_1644, pc_x, snh_1234, snh_1235, \
                         snh_1236, snh_1237, snh_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1640[k] = f_3 * pc_x[k] * snh_1234[k];

        t_1641[k] = f_3 * pc_x[k] * snh_1235[k];

        t_1642[k] = f_3 * pc_x[k] * snh_1236[k];

        t_1643[k] = f_3 * pc_x[k] * snh_1237[k];

        t_1644[k] = f_3 * pc_x[k] * snh_1238[k];
    }

#pragma omp simd aligned(t_1645, t_1646, t_1647, pc_y, pc_z, smh_1002, smh_1023, smh_1025, \
                         sng0_880, sng0_882, sng1_880, sng1_882, snh_1233, \
                         snh_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1645[k] = f_17 * smh_1023[k]
                    + f_1 * sng0_880[k]
                    - f_2 * sng1_880[k]
                    + f_3 * pc_y[k] * snh_1233[k];

        t_1646[k] = f_13 * smh_1002[k]
                    + f_3 * pc_z[k] * snh_1233[k];

        t_1647[k] = f_17 * smh_1025[k]
                    + f_4 * sng0_882[k]
                    - f_5 * sng1_882[k]
                    + f_3 * pc_y[k] * snh_1235[k];
    }

#pragma omp simd aligned(t_1648, t_1649, t_1650, pc_y, smh_1026, smh_1027, smh_1028, sng0_883, \
                         sng0_884, sng1_883, sng1_884, snh_1236, snh_1237, \
                         snh_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1648[k] = f_17 * smh_1026[k]
                    + f_6 * sng0_883[k]
                    - f_7 * sng1_883[k]
                    + f_3 * pc_y[k] * snh_1236[k];

        t_1649[k] = f_17 * smh_1027[k]
                    + f_8 * sng0_884[k]
                    - f_9 * sng1_884[k]
                    + f_3 * pc_y[k] * snh_1237[k];

        t_1650[k] = f_17 * smh_1028[k]
                    + f_3 * pc_y[k] * snh_1238[k];
    }

#pragma omp simd aligned(t_1651, t_1652, t_1653, t_1654, pc_x, pc_y, pc_z, smh_1007, smh_1008, \
                         smh_1029, sng0_884, sng0_885, sng1_884, sng1_885, snh_1238, \
                         snh_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1651[k] = f_13 * smh_1007[k]
                    + f_1 * sng0_884[k]
                    - f_2 * sng1_884[k]
                    + f_3 * pc_z[k] * snh_1238[k];

        t_1652[k] = f_1 * sng0_885[k]
                    - f_2 * sng1_885[k]
                    + f_3 * pc_x[k] * snh_1239[k];

        t_1653[k] = f_18 * smh_1029[k]
                    + f_3 * pc_y[k] * snh_1239[k];

        t_1654[k] = f_14 * smh_1008[k]
                    + f_3 * pc_z[k] * snh_1239[k];
    }

#pragma omp simd aligned(t_1655, t_1656, t_1657, pc_x, pc_y, smh_1031, sng0_888, sng0_890, \
                         sng1_888, sng1_890, snh_1241, snh_1242, \
                         snh_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1655[k] = f_4 * sng0_888[k]
                    - f_5 * sng1_888[k]
                    + f_3 * pc_x[k] * snh_1242[k];

        t_1656[k] = f_18 * smh_1031[k]
                    + f_3 * pc_y[k] * snh_1241[k];

        t_1657[k] = f_4 * sng0_890[k]
                    - f_5 * sng1_890[k]
                    + f_3 * pc_x[k] * snh_1244[k];
    }

#pragma omp simd aligned(t_1658, t_1659, t_1660, pc_x, pc_y, pc_z, smh_1011, smh_1034, \
                         sng0_891, sng1_891, snh_1242, snh_1244, \
                         snh_1245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1658[k] = f_6 * sng0_891[k]
                    - f_7 * sng1_891[k]
                    + f_3 * pc_x[k] * snh_1245[k];

        t_1659[k] = f_14 * smh_1011[k]
                    + f_3 * pc_z[k] * snh_1242[k];

        t_1660[k] = f_18 * smh_1034[k]
                    + f_3 * pc_y[k] * snh_1244[k];
    }

#pragma omp simd aligned(t_1661, t_1662, t_1663, pc_x, pc_z, smh_1014, sng0_894, sng0_895, \
                         sng1_894, sng1_895, snh_1245, snh_1248, \
                         snh_1249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1661[k] = f_6 * sng0_894[k]
                    - f_7 * sng1_894[k]
                    + f_3 * pc_x[k] * snh_1248[k];

        t_1662[k] = f_8 * sng0_895[k]
                    - f_9 * sng1_895[k]
                    + f_3 * pc_x[k] * snh_1249[k];

        t_1663[k] = f_14 * smh_1014[k]
                    + f_3 * pc_z[k] * snh_1245[k];
    }

#pragma omp simd aligned(t_1664, t_1665, t_1666, t_1667, pc_x, pc_y, smh_1038, sng0_897, \
                         sng0_899, sng1_897, sng1_899, snh_1248, snh_1251, snh_1253, \
                         snh_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1664[k] = f_8 * sng0_897[k]
                    - f_9 * sng1_897[k]
                    + f_3 * pc_x[k] * snh_1251[k];

        t_1665[k] = f_18 * smh_1038[k]
                    + f_3 * pc_y[k] * snh_1248[k];

        t_1666[k] = f_8 * sng0_899[k]
                    - f_9 * sng1_899[k]
                    + f_3 * pc_x[k] * snh_1253[k];

        t_1667[k] = f_3 * pc_x[k] * snh_1254[k];
    }

#pragma omp simd aligned(t_1668, t_1669, t_1670, t_1671, t_1672, pc_x, snh_1255, snh_1256, \
                         snh_1257, snh_1258, snh_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1668[k] = f_3 * pc_x[k] * snh_1255[k];

        t_1669[k] = f_3 * pc_x[k] * snh_1256[k];

        t_1670[k] = f_3 * pc_x[k] * snh_1257[k];

        t_1671[k] = f_3 * pc_x[k] * snh_1258[k];

        t_1672[k] = f_3 * pc_x[k] * snh_1259[k];
    }

#pragma omp simd aligned(t_1673, t_1674, t_1675, pc_y, pc_z, smh_1023, smh_1044, smh_1046, \
                         sng0_895, sng0_897, sng1_895, sng1_897, snh_1254, \
                         snh_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1673[k] = f_18 * smh_1044[k]
                    + f_1 * sng0_895[k]
                    - f_2 * sng1_895[k]
                    + f_3 * pc_y[k] * snh_1254[k];

        t_1674[k] = f_14 * smh_1023[k]
                    + f_3 * pc_z[k] * snh_1254[k];

        t_1675[k] = f_18 * smh_1046[k]
                    + f_4 * sng0_897[k]
                    - f_5 * sng1_897[k]
                    + f_3 * pc_y[k] * snh_1256[k];
    }

#pragma omp simd aligned(t_1676, t_1677, t_1678, pc_y, smh_1047, smh_1048, smh_1049, sng0_898, \
                         sng0_899, sng1_898, sng1_899, snh_1257, snh_1258, \
                         snh_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1676[k] = f_18 * smh_1047[k]
                    + f_6 * sng0_898[k]
                    - f_7 * sng1_898[k]
                    + f_3 * pc_y[k] * snh_1257[k];

        t_1677[k] = f_18 * smh_1048[k]
                    + f_8 * sng0_899[k]
                    - f_9 * sng1_899[k]
                    + f_3 * pc_y[k] * snh_1258[k];

        t_1678[k] = f_18 * smh_1049[k]
                    + f_3 * pc_y[k] * snh_1259[k];
    }

#pragma omp simd aligned(t_1679, t_1680, t_1681, t_1682, pc_x, pc_y, pc_z, smh_1028, smh_1029, \
                         smh_1050, sng0_899, sng0_900, sng1_899, sng1_900, snh_1259, \
                         snh_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1679[k] = f_14 * smh_1028[k]
                    + f_1 * sng0_899[k]
                    - f_2 * sng1_899[k]
                    + f_3 * pc_z[k] * snh_1259[k];

        t_1680[k] = f_1 * sng0_900[k]
                    - f_2 * sng1_900[k]
                    + f_3 * pc_x[k] * snh_1260[k];

        t_1681[k] = f_19 * smh_1050[k]
                    + f_3 * pc_y[k] * snh_1260[k];

        t_1682[k] = f_19 * smh_1029[k]
                    + f_3 * pc_z[k] * snh_1260[k];
    }

#pragma omp simd aligned(t_1683, t_1684, t_1685, pc_x, pc_y, smh_1052, sng0_903, sng0_905, \
                         sng1_903, sng1_905, snh_1262, snh_1263, \
                         snh_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1683[k] = f_4 * sng0_903[k]
                    - f_5 * sng1_903[k]
                    + f_3 * pc_x[k] * snh_1263[k];

        t_1684[k] = f_19 * smh_1052[k]
                    + f_3 * pc_y[k] * snh_1262[k];

        t_1685[k] = f_4 * sng0_905[k]
                    - f_5 * sng1_905[k]
                    + f_3 * pc_x[k] * snh_1265[k];
    }

#pragma omp simd aligned(t_1686, t_1687, t_1688, pc_x, pc_y, pc_z, smh_1032, smh_1055, \
                         sng0_906, sng1_906, snh_1263, snh_1265, \
                         snh_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1686[k] = f_6 * sng0_906[k]
                    - f_7 * sng1_906[k]
                    + f_3 * pc_x[k] * snh_1266[k];

        t_1687[k] = f_19 * smh_1032[k]
                    + f_3 * pc_z[k] * snh_1263[k];

        t_1688[k] = f_19 * smh_1055[k]
                    + f_3 * pc_y[k] * snh_1265[k];
    }

#pragma omp simd aligned(t_1689, t_1690, t_1691, pc_x, pc_z, smh_1035, sng0_909, sng0_910, \
                         sng1_909, sng1_910, snh_1266, snh_1269, \
                         snh_1270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1689[k] = f_6 * sng0_909[k]
                    - f_7 * sng1_909[k]
                    + f_3 * pc_x[k] * snh_1269[k];

        t_1690[k] = f_8 * sng0_910[k]
                    - f_9 * sng1_910[k]
                    + f_3 * pc_x[k] * snh_1270[k];

        t_1691[k] = f_19 * smh_1035[k]
                    + f_3 * pc_z[k] * snh_1266[k];
    }

#pragma omp simd aligned(t_1692, t_1693, t_1694, t_1695, pc_x, pc_y, smh_1059, sng0_912, \
                         sng0_914, sng1_912, sng1_914, snh_1269, snh_1272, snh_1274, \
                         snh_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1692[k] = f_8 * sng0_912[k]
                    - f_9 * sng1_912[k]
                    + f_3 * pc_x[k] * snh_1272[k];

        t_1693[k] = f_19 * smh_1059[k]
                    + f_3 * pc_y[k] * snh_1269[k];

        t_1694[k] = f_8 * sng0_914[k]
                    - f_9 * sng1_914[k]
                    + f_3 * pc_x[k] * snh_1274[k];

        t_1695[k] = f_3 * pc_x[k] * snh_1275[k];
    }

#pragma omp simd aligned(t_1696, t_1697, t_1698, t_1699, t_1700, pc_x, snh_1276, snh_1277, \
                         snh_1278, snh_1279, snh_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1696[k] = f_3 * pc_x[k] * snh_1276[k];

        t_1697[k] = f_3 * pc_x[k] * snh_1277[k];

        t_1698[k] = f_3 * pc_x[k] * snh_1278[k];

        t_1699[k] = f_3 * pc_x[k] * snh_1279[k];

        t_1700[k] = f_3 * pc_x[k] * snh_1280[k];
    }

#pragma omp simd aligned(t_1701, t_1702, t_1703, pc_y, pc_z, smh_1044, smh_1065, smh_1067, \
                         sng0_910, sng0_912, sng1_910, sng1_912, snh_1275, \
                         snh_1277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1701[k] = f_19 * smh_1065[k]
                    + f_1 * sng0_910[k]
                    - f_2 * sng1_910[k]
                    + f_3 * pc_y[k] * snh_1275[k];

        t_1702[k] = f_19 * smh_1044[k]
                    + f_3 * pc_z[k] * snh_1275[k];

        t_1703[k] = f_19 * smh_1067[k]
                    + f_4 * sng0_912[k]
                    - f_5 * sng1_912[k]
                    + f_3 * pc_y[k] * snh_1277[k];
    }

#pragma omp simd aligned(t_1704, t_1705, t_1706, pc_y, smh_1068, smh_1069, smh_1070, sng0_913, \
                         sng0_914, sng1_913, sng1_914, snh_1278, snh_1279, \
                         snh_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1704[k] = f_19 * smh_1068[k]
                    + f_6 * sng0_913[k]
                    - f_7 * sng1_913[k]
                    + f_3 * pc_y[k] * snh_1278[k];

        t_1705[k] = f_19 * smh_1069[k]
                    + f_8 * sng0_914[k]
                    - f_9 * sng1_914[k]
                    + f_3 * pc_y[k] * snh_1279[k];

        t_1706[k] = f_19 * smh_1070[k]
                    + f_3 * pc_y[k] * snh_1280[k];
    }

#pragma omp simd aligned(t_1707, t_1708, t_1709, t_1710, pc_x, pc_y, pc_z, smh_1049, smh_1050, \
                         smh_1071, sng0_914, sng0_915, sng1_914, sng1_915, snh_1280, \
                         snh_1281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1707[k] = f_19 * smh_1049[k]
                    + f_1 * sng0_914[k]
                    - f_2 * sng1_914[k]
                    + f_3 * pc_z[k] * snh_1280[k];

        t_1708[k] = f_1 * sng0_915[k]
                    - f_2 * sng1_915[k]
                    + f_3 * pc_x[k] * snh_1281[k];

        t_1709[k] = f_14 * smh_1071[k]
                    + f_3 * pc_y[k] * snh_1281[k];

        t_1710[k] = f_18 * smh_1050[k]
                    + f_3 * pc_z[k] * snh_1281[k];
    }

#pragma omp simd aligned(t_1711, t_1712, t_1713, pc_x, pc_y, smh_1073, sng0_918, sng0_920, \
                         sng1_918, sng1_920, snh_1283, snh_1284, \
                         snh_1286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1711[k] = f_4 * sng0_918[k]
                    - f_5 * sng1_918[k]
                    + f_3 * pc_x[k] * snh_1284[k];

        t_1712[k] = f_14 * smh_1073[k]
                    + f_3 * pc_y[k] * snh_1283[k];

        t_1713[k] = f_4 * sng0_920[k]
                    - f_5 * sng1_920[k]
                    + f_3 * pc_x[k] * snh_1286[k];
    }

#pragma omp simd aligned(t_1714, t_1715, t_1716, pc_x, pc_y, pc_z, smh_1053, smh_1076, \
                         sng0_921, sng1_921, snh_1284, snh_1286, \
                         snh_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1714[k] = f_6 * sng0_921[k]
                    - f_7 * sng1_921[k]
                    + f_3 * pc_x[k] * snh_1287[k];

        t_1715[k] = f_18 * smh_1053[k]
                    + f_3 * pc_z[k] * snh_1284[k];

        t_1716[k] = f_14 * smh_1076[k]
                    + f_3 * pc_y[k] * snh_1286[k];
    }

#pragma omp simd aligned(t_1717, t_1718, t_1719, pc_x, pc_z, smh_1056, sng0_924, sng0_925, \
                         sng1_924, sng1_925, snh_1287, snh_1290, \
                         snh_1291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1717[k] = f_6 * sng0_924[k]
                    - f_7 * sng1_924[k]
                    + f_3 * pc_x[k] * snh_1290[k];

        t_1718[k] = f_8 * sng0_925[k]
                    - f_9 * sng1_925[k]
                    + f_3 * pc_x[k] * snh_1291[k];

        t_1719[k] = f_18 * smh_1056[k]
                    + f_3 * pc_z[k] * snh_1287[k];
    }

#pragma omp simd aligned(t_1720, t_1721, t_1722, t_1723, pc_x, pc_y, smh_1080, sng0_927, \
                         sng0_929, sng1_927, sng1_929, snh_1290, snh_1293, snh_1295, \
                         snh_1296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1720[k] = f_8 * sng0_927[k]
                    - f_9 * sng1_927[k]
                    + f_3 * pc_x[k] * snh_1293[k];

        t_1721[k] = f_14 * smh_1080[k]
                    + f_3 * pc_y[k] * snh_1290[k];

        t_1722[k] = f_8 * sng0_929[k]
                    - f_9 * sng1_929[k]
                    + f_3 * pc_x[k] * snh_1295[k];

        t_1723[k] = f_3 * pc_x[k] * snh_1296[k];
    }

#pragma omp simd aligned(t_1724, t_1725, t_1726, t_1727, t_1728, pc_x, snh_1297, snh_1298, \
                         snh_1299, snh_1300, snh_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1724[k] = f_3 * pc_x[k] * snh_1297[k];

        t_1725[k] = f_3 * pc_x[k] * snh_1298[k];

        t_1726[k] = f_3 * pc_x[k] * snh_1299[k];

        t_1727[k] = f_3 * pc_x[k] * snh_1300[k];

        t_1728[k] = f_3 * pc_x[k] * snh_1301[k];
    }
}

static auto
compute_prim_sni_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smi0,
                                                           const size_t smh, const size_t smi1,
                                                           const size_t sng0, const size_t sng1,
                                                           const size_t snh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_15 = 4.5 / q;
    const auto f_16 = 4.0 / q;
    const auto f_17 = 3.5 / q;
    const auto f_18 = 3.0 / q;

    auto *t_1729 = buffer.data(target + 1729);
    auto *t_1730 = buffer.data(target + 1730);
    auto *t_1731 = buffer.data(target + 1731);
    auto *t_1732 = buffer.data(target + 1732);
    auto *t_1733 = buffer.data(target + 1733);
    auto *t_1734 = buffer.data(target + 1734);
    auto *t_1735 = buffer.data(target + 1735);
    auto *t_1736 = buffer.data(target + 1736);
    auto *t_1737 = buffer.data(target + 1737);
    auto *t_1738 = buffer.data(target + 1738);
    auto *t_1739 = buffer.data(target + 1739);
    auto *t_1740 = buffer.data(target + 1740);
    auto *t_1741 = buffer.data(target + 1741);
    auto *t_1742 = buffer.data(target + 1742);
    auto *t_1743 = buffer.data(target + 1743);
    auto *t_1744 = buffer.data(target + 1744);
    auto *t_1745 = buffer.data(target + 1745);
    auto *t_1746 = buffer.data(target + 1746);
    auto *t_1747 = buffer.data(target + 1747);
    auto *t_1748 = buffer.data(target + 1748);
    auto *t_1749 = buffer.data(target + 1749);
    auto *t_1750 = buffer.data(target + 1750);
    auto *t_1751 = buffer.data(target + 1751);
    auto *t_1752 = buffer.data(target + 1752);
    auto *t_1753 = buffer.data(target + 1753);
    auto *t_1754 = buffer.data(target + 1754);
    auto *t_1755 = buffer.data(target + 1755);
    auto *t_1756 = buffer.data(target + 1756);
    auto *t_1757 = buffer.data(target + 1757);
    auto *t_1758 = buffer.data(target + 1758);
    auto *t_1759 = buffer.data(target + 1759);
    auto *t_1760 = buffer.data(target + 1760);
    auto *t_1761 = buffer.data(target + 1761);
    auto *t_1762 = buffer.data(target + 1762);
    auto *t_1763 = buffer.data(target + 1763);
    auto *t_1764 = buffer.data(target + 1764);
    auto *t_1765 = buffer.data(target + 1765);
    auto *t_1766 = buffer.data(target + 1766);
    auto *t_1767 = buffer.data(target + 1767);
    auto *t_1768 = buffer.data(target + 1768);
    auto *t_1769 = buffer.data(target + 1769);
    auto *t_1770 = buffer.data(target + 1770);
    auto *t_1771 = buffer.data(target + 1771);
    auto *t_1772 = buffer.data(target + 1772);
    auto *t_1773 = buffer.data(target + 1773);
    auto *t_1774 = buffer.data(target + 1774);
    auto *t_1775 = buffer.data(target + 1775);
    auto *t_1776 = buffer.data(target + 1776);
    auto *t_1777 = buffer.data(target + 1777);
    auto *t_1778 = buffer.data(target + 1778);
    auto *t_1779 = buffer.data(target + 1779);
    auto *t_1780 = buffer.data(target + 1780);
    auto *t_1781 = buffer.data(target + 1781);
    auto *t_1782 = buffer.data(target + 1782);
    auto *t_1783 = buffer.data(target + 1783);
    auto *t_1784 = buffer.data(target + 1784);
    auto *t_1785 = buffer.data(target + 1785);
    auto *t_1786 = buffer.data(target + 1786);
    auto *t_1787 = buffer.data(target + 1787);
    auto *t_1788 = buffer.data(target + 1788);
    auto *t_1789 = buffer.data(target + 1789);
    auto *t_1790 = buffer.data(target + 1790);
    auto *t_1791 = buffer.data(target + 1791);
    auto *t_1792 = buffer.data(target + 1792);
    auto *t_1793 = buffer.data(target + 1793);
    auto *t_1794 = buffer.data(target + 1794);
    auto *t_1795 = buffer.data(target + 1795);
    auto *t_1796 = buffer.data(target + 1796);
    auto *t_1797 = buffer.data(target + 1797);
    auto *t_1798 = buffer.data(target + 1798);
    auto *t_1799 = buffer.data(target + 1799);
    auto *t_1800 = buffer.data(target + 1800);
    auto *t_1801 = buffer.data(target + 1801);
    auto *t_1802 = buffer.data(target + 1802);
    auto *t_1803 = buffer.data(target + 1803);
    auto *t_1804 = buffer.data(target + 1804);
    auto *t_1805 = buffer.data(target + 1805);
    auto *t_1806 = buffer.data(target + 1806);
    auto *t_1807 = buffer.data(target + 1807);
    auto *t_1808 = buffer.data(target + 1808);
    auto *t_1809 = buffer.data(target + 1809);
    auto *t_1810 = buffer.data(target + 1810);
    auto *t_1811 = buffer.data(target + 1811);
    auto *t_1812 = buffer.data(target + 1812);
    auto *t_1813 = buffer.data(target + 1813);
    auto *t_1814 = buffer.data(target + 1814);
    auto *t_1815 = buffer.data(target + 1815);
    auto *t_1816 = buffer.data(target + 1816);
    auto *t_1817 = buffer.data(target + 1817);
    auto *t_1818 = buffer.data(target + 1818);
    auto *t_1819 = buffer.data(target + 1819);
    auto *t_1820 = buffer.data(target + 1820);
    auto *t_1821 = buffer.data(target + 1821);
    auto *t_1822 = buffer.data(target + 1822);
    auto *t_1823 = buffer.data(target + 1823);
    auto *t_1824 = buffer.data(target + 1824);
    auto *t_1825 = buffer.data(target + 1825);
    auto *t_1826 = buffer.data(target + 1826);
    auto *t_1827 = buffer.data(target + 1827);
    auto *t_1828 = buffer.data(target + 1828);
    auto *t_1829 = buffer.data(target + 1829);
    auto *t_1830 = buffer.data(target + 1830);
    auto *t_1831 = buffer.data(target + 1831);
    auto *t_1832 = buffer.data(target + 1832);
    auto *t_1833 = buffer.data(target + 1833);
    auto *t_1834 = buffer.data(target + 1834);
    auto *t_1835 = buffer.data(target + 1835);
    auto *t_1836 = buffer.data(target + 1836);
    auto *t_1837 = buffer.data(target + 1837);
    auto *t_1838 = buffer.data(target + 1838);
    auto *t_1839 = buffer.data(target + 1839);
    auto *t_1840 = buffer.data(target + 1840);
    auto *t_1841 = buffer.data(target + 1841);
    auto *t_1842 = buffer.data(target + 1842);
    auto *t_1843 = buffer.data(target + 1843);
    auto *t_1844 = buffer.data(target + 1844);
    auto *t_1845 = buffer.data(target + 1845);
    auto *t_1846 = buffer.data(target + 1846);
    auto *t_1847 = buffer.data(target + 1847);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi0_1512 = buffer.data(smi0 + 1512);
    const auto *smi0_1517 = buffer.data(smi0 + 1517);
    const auto *smi0_1521 = buffer.data(smi0 + 1521);
    const auto *smi0_1526 = buffer.data(smi0 + 1526);
    const auto *smi0_1533 = buffer.data(smi0 + 1533);
    const auto *smi0_1535 = buffer.data(smi0 + 1535);
    const auto *smi0_1536 = buffer.data(smi0 + 1536);
    const auto *smi0_1537 = buffer.data(smi0 + 1537);
    const auto *smi0_1539 = buffer.data(smi0 + 1539);

    const auto *smh_1065 = buffer.data(smh + 1065);
    const auto *smh_1070 = buffer.data(smh + 1070);
    const auto *smh_1071 = buffer.data(smh + 1071);
    const auto *smh_1074 = buffer.data(smh + 1074);
    const auto *smh_1077 = buffer.data(smh + 1077);
    const auto *smh_1086 = buffer.data(smh + 1086);
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
    const auto *smh_1107 = buffer.data(smh + 1107);
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
    const auto *smh_1128 = buffer.data(smh + 1128);
    const auto *smh_1130 = buffer.data(smh + 1130);
    const auto *smh_1131 = buffer.data(smh + 1131);
    const auto *smh_1132 = buffer.data(smh + 1132);
    const auto *smh_1133 = buffer.data(smh + 1133);
    const auto *smh_1134 = buffer.data(smh + 1134);
    const auto *smh_1136 = buffer.data(smh + 1136);
    const auto *smh_1137 = buffer.data(smh + 1137);
    const auto *smh_1139 = buffer.data(smh + 1139);
    const auto *smh_1140 = buffer.data(smh + 1140);
    const auto *smh_1143 = buffer.data(smh + 1143);
    const auto *smh_1149 = buffer.data(smh + 1149);
    const auto *smh_1151 = buffer.data(smh + 1151);
    const auto *smh_1152 = buffer.data(smh + 1152);
    const auto *smh_1153 = buffer.data(smh + 1153);
    const auto *smh_1154 = buffer.data(smh + 1154);

    const auto *smi1_1512 = buffer.data(smi1 + 1512);
    const auto *smi1_1517 = buffer.data(smi1 + 1517);
    const auto *smi1_1521 = buffer.data(smi1 + 1521);
    const auto *smi1_1526 = buffer.data(smi1 + 1526);
    const auto *smi1_1533 = buffer.data(smi1 + 1533);
    const auto *smi1_1535 = buffer.data(smi1 + 1535);
    const auto *smi1_1536 = buffer.data(smi1 + 1536);
    const auto *smi1_1537 = buffer.data(smi1 + 1537);
    const auto *smi1_1539 = buffer.data(smi1 + 1539);

    const auto *sng0_925 = buffer.data(sng0 + 925);
    const auto *sng0_927 = buffer.data(sng0 + 927);
    const auto *sng0_928 = buffer.data(sng0 + 928);
    const auto *sng0_929 = buffer.data(sng0 + 929);
    const auto *sng0_930 = buffer.data(sng0 + 930);
    const auto *sng0_933 = buffer.data(sng0 + 933);
    const auto *sng0_935 = buffer.data(sng0 + 935);
    const auto *sng0_936 = buffer.data(sng0 + 936);
    const auto *sng0_939 = buffer.data(sng0 + 939);
    const auto *sng0_940 = buffer.data(sng0 + 940);
    const auto *sng0_942 = buffer.data(sng0 + 942);
    const auto *sng0_943 = buffer.data(sng0 + 943);
    const auto *sng0_944 = buffer.data(sng0 + 944);
    const auto *sng0_945 = buffer.data(sng0 + 945);
    const auto *sng0_948 = buffer.data(sng0 + 948);
    const auto *sng0_950 = buffer.data(sng0 + 950);
    const auto *sng0_951 = buffer.data(sng0 + 951);
    const auto *sng0_954 = buffer.data(sng0 + 954);
    const auto *sng0_955 = buffer.data(sng0 + 955);
    const auto *sng0_957 = buffer.data(sng0 + 957);
    const auto *sng0_958 = buffer.data(sng0 + 958);
    const auto *sng0_959 = buffer.data(sng0 + 959);
    const auto *sng0_963 = buffer.data(sng0 + 963);
    const auto *sng0_966 = buffer.data(sng0 + 966);
    const auto *sng0_970 = buffer.data(sng0 + 970);
    const auto *sng0_972 = buffer.data(sng0 + 972);
    const auto *sng0_975 = buffer.data(sng0 + 975);
    const auto *sng0_978 = buffer.data(sng0 + 978);
    const auto *sng0_980 = buffer.data(sng0 + 980);
    const auto *sng0_981 = buffer.data(sng0 + 981);
    const auto *sng0_984 = buffer.data(sng0 + 984);
    const auto *sng0_985 = buffer.data(sng0 + 985);
    const auto *sng0_987 = buffer.data(sng0 + 987);
    const auto *sng0_988 = buffer.data(sng0 + 988);
    const auto *sng0_989 = buffer.data(sng0 + 989);

    const auto *sng1_925 = buffer.data(sng1 + 925);
    const auto *sng1_927 = buffer.data(sng1 + 927);
    const auto *sng1_928 = buffer.data(sng1 + 928);
    const auto *sng1_929 = buffer.data(sng1 + 929);
    const auto *sng1_930 = buffer.data(sng1 + 930);
    const auto *sng1_933 = buffer.data(sng1 + 933);
    const auto *sng1_935 = buffer.data(sng1 + 935);
    const auto *sng1_936 = buffer.data(sng1 + 936);
    const auto *sng1_939 = buffer.data(sng1 + 939);
    const auto *sng1_940 = buffer.data(sng1 + 940);
    const auto *sng1_942 = buffer.data(sng1 + 942);
    const auto *sng1_943 = buffer.data(sng1 + 943);
    const auto *sng1_944 = buffer.data(sng1 + 944);
    const auto *sng1_945 = buffer.data(sng1 + 945);
    const auto *sng1_948 = buffer.data(sng1 + 948);
    const auto *sng1_950 = buffer.data(sng1 + 950);
    const auto *sng1_951 = buffer.data(sng1 + 951);
    const auto *sng1_954 = buffer.data(sng1 + 954);
    const auto *sng1_955 = buffer.data(sng1 + 955);
    const auto *sng1_957 = buffer.data(sng1 + 957);
    const auto *sng1_958 = buffer.data(sng1 + 958);
    const auto *sng1_959 = buffer.data(sng1 + 959);
    const auto *sng1_963 = buffer.data(sng1 + 963);
    const auto *sng1_966 = buffer.data(sng1 + 966);
    const auto *sng1_970 = buffer.data(sng1 + 970);
    const auto *sng1_972 = buffer.data(sng1 + 972);
    const auto *sng1_975 = buffer.data(sng1 + 975);
    const auto *sng1_978 = buffer.data(sng1 + 978);
    const auto *sng1_980 = buffer.data(sng1 + 980);
    const auto *sng1_981 = buffer.data(sng1 + 981);
    const auto *sng1_984 = buffer.data(sng1 + 984);
    const auto *sng1_985 = buffer.data(sng1 + 985);
    const auto *sng1_987 = buffer.data(sng1 + 987);
    const auto *sng1_988 = buffer.data(sng1 + 988);
    const auto *sng1_989 = buffer.data(sng1 + 989);

    const auto *snh_1296 = buffer.data(snh + 1296);
    const auto *snh_1298 = buffer.data(snh + 1298);
    const auto *snh_1299 = buffer.data(snh + 1299);
    const auto *snh_1300 = buffer.data(snh + 1300);
    const auto *snh_1301 = buffer.data(snh + 1301);
    const auto *snh_1302 = buffer.data(snh + 1302);
    const auto *snh_1304 = buffer.data(snh + 1304);
    const auto *snh_1305 = buffer.data(snh + 1305);
    const auto *snh_1307 = buffer.data(snh + 1307);
    const auto *snh_1308 = buffer.data(snh + 1308);
    const auto *snh_1311 = buffer.data(snh + 1311);
    const auto *snh_1312 = buffer.data(snh + 1312);
    const auto *snh_1314 = buffer.data(snh + 1314);
    const auto *snh_1316 = buffer.data(snh + 1316);
    const auto *snh_1317 = buffer.data(snh + 1317);
    const auto *snh_1318 = buffer.data(snh + 1318);
    const auto *snh_1319 = buffer.data(snh + 1319);
    const auto *snh_1320 = buffer.data(snh + 1320);
    const auto *snh_1321 = buffer.data(snh + 1321);
    const auto *snh_1322 = buffer.data(snh + 1322);
    const auto *snh_1323 = buffer.data(snh + 1323);
    const auto *snh_1325 = buffer.data(snh + 1325);
    const auto *snh_1326 = buffer.data(snh + 1326);
    const auto *snh_1328 = buffer.data(snh + 1328);
    const auto *snh_1329 = buffer.data(snh + 1329);
    const auto *snh_1332 = buffer.data(snh + 1332);
    const auto *snh_1333 = buffer.data(snh + 1333);
    const auto *snh_1335 = buffer.data(snh + 1335);
    const auto *snh_1337 = buffer.data(snh + 1337);
    const auto *snh_1338 = buffer.data(snh + 1338);
    const auto *snh_1339 = buffer.data(snh + 1339);
    const auto *snh_1340 = buffer.data(snh + 1340);
    const auto *snh_1341 = buffer.data(snh + 1341);
    const auto *snh_1342 = buffer.data(snh + 1342);
    const auto *snh_1343 = buffer.data(snh + 1343);
    const auto *snh_1344 = buffer.data(snh + 1344);
    const auto *snh_1346 = buffer.data(snh + 1346);
    const auto *snh_1347 = buffer.data(snh + 1347);
    const auto *snh_1349 = buffer.data(snh + 1349);
    const auto *snh_1350 = buffer.data(snh + 1350);
    const auto *snh_1353 = buffer.data(snh + 1353);
    const auto *snh_1354 = buffer.data(snh + 1354);
    const auto *snh_1356 = buffer.data(snh + 1356);
    const auto *snh_1359 = buffer.data(snh + 1359);
    const auto *snh_1360 = buffer.data(snh + 1360);
    const auto *snh_1361 = buffer.data(snh + 1361);
    const auto *snh_1362 = buffer.data(snh + 1362);
    const auto *snh_1363 = buffer.data(snh + 1363);
    const auto *snh_1364 = buffer.data(snh + 1364);
    const auto *snh_1365 = buffer.data(snh + 1365);
    const auto *snh_1367 = buffer.data(snh + 1367);
    const auto *snh_1368 = buffer.data(snh + 1368);
    const auto *snh_1370 = buffer.data(snh + 1370);
    const auto *snh_1371 = buffer.data(snh + 1371);
    const auto *snh_1374 = buffer.data(snh + 1374);
    const auto *snh_1375 = buffer.data(snh + 1375);
    const auto *snh_1377 = buffer.data(snh + 1377);
    const auto *snh_1379 = buffer.data(snh + 1379);
    const auto *snh_1380 = buffer.data(snh + 1380);
    const auto *snh_1381 = buffer.data(snh + 1381);
    const auto *snh_1382 = buffer.data(snh + 1382);
    const auto *snh_1383 = buffer.data(snh + 1383);
    const auto *snh_1384 = buffer.data(snh + 1384);
    const auto *snh_1385 = buffer.data(snh + 1385);

#pragma omp simd aligned(t_1729, t_1730, t_1731, pc_y, pc_z, smh_1065, smh_1086, smh_1088, \
                         sng0_925, sng0_927, sng1_925, sng1_927, snh_1296, \
                         snh_1298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1729[k] = f_14 * smh_1086[k]
                    + f_1 * sng0_925[k]
                    - f_2 * sng1_925[k]
                    + f_3 * pc_y[k] * snh_1296[k];

        t_1730[k] = f_18 * smh_1065[k]
                    + f_3 * pc_z[k] * snh_1296[k];

        t_1731[k] = f_14 * smh_1088[k]
                    + f_4 * sng0_927[k]
                    - f_5 * sng1_927[k]
                    + f_3 * pc_y[k] * snh_1298[k];
    }

#pragma omp simd aligned(t_1732, t_1733, t_1734, pc_y, smh_1089, smh_1090, smh_1091, sng0_928, \
                         sng0_929, sng1_928, sng1_929, snh_1299, snh_1300, \
                         snh_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1732[k] = f_14 * smh_1089[k]
                    + f_6 * sng0_928[k]
                    - f_7 * sng1_928[k]
                    + f_3 * pc_y[k] * snh_1299[k];

        t_1733[k] = f_14 * smh_1090[k]
                    + f_8 * sng0_929[k]
                    - f_9 * sng1_929[k]
                    + f_3 * pc_y[k] * snh_1300[k];

        t_1734[k] = f_14 * smh_1091[k]
                    + f_3 * pc_y[k] * snh_1301[k];
    }

#pragma omp simd aligned(t_1735, t_1736, t_1737, t_1738, pc_x, pc_y, pc_z, smh_1070, smh_1071, \
                         smh_1092, sng0_929, sng0_930, sng1_929, sng1_930, snh_1301, \
                         snh_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1735[k] = f_18 * smh_1070[k]
                    + f_1 * sng0_929[k]
                    - f_2 * sng1_929[k]
                    + f_3 * pc_z[k] * snh_1301[k];

        t_1736[k] = f_1 * sng0_930[k]
                    - f_2 * sng1_930[k]
                    + f_3 * pc_x[k] * snh_1302[k];

        t_1737[k] = f_13 * smh_1092[k]
                    + f_3 * pc_y[k] * snh_1302[k];

        t_1738[k] = f_17 * smh_1071[k]
                    + f_3 * pc_z[k] * snh_1302[k];
    }

#pragma omp simd aligned(t_1739, t_1740, t_1741, pc_x, pc_y, smh_1094, sng0_933, sng0_935, \
                         sng1_933, sng1_935, snh_1304, snh_1305, \
                         snh_1307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1739[k] = f_4 * sng0_933[k]
                    - f_5 * sng1_933[k]
                    + f_3 * pc_x[k] * snh_1305[k];

        t_1740[k] = f_13 * smh_1094[k]
                    + f_3 * pc_y[k] * snh_1304[k];

        t_1741[k] = f_4 * sng0_935[k]
                    - f_5 * sng1_935[k]
                    + f_3 * pc_x[k] * snh_1307[k];
    }

#pragma omp simd aligned(t_1742, t_1743, t_1744, pc_x, pc_y, pc_z, smh_1074, smh_1097, \
                         sng0_936, sng1_936, snh_1305, snh_1307, \
                         snh_1308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1742[k] = f_6 * sng0_936[k]
                    - f_7 * sng1_936[k]
                    + f_3 * pc_x[k] * snh_1308[k];

        t_1743[k] = f_17 * smh_1074[k]
                    + f_3 * pc_z[k] * snh_1305[k];

        t_1744[k] = f_13 * smh_1097[k]
                    + f_3 * pc_y[k] * snh_1307[k];
    }

#pragma omp simd aligned(t_1745, t_1746, t_1747, pc_x, pc_z, smh_1077, sng0_939, sng0_940, \
                         sng1_939, sng1_940, snh_1308, snh_1311, \
                         snh_1312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1745[k] = f_6 * sng0_939[k]
                    - f_7 * sng1_939[k]
                    + f_3 * pc_x[k] * snh_1311[k];

        t_1746[k] = f_8 * sng0_940[k]
                    - f_9 * sng1_940[k]
                    + f_3 * pc_x[k] * snh_1312[k];

        t_1747[k] = f_17 * smh_1077[k]
                    + f_3 * pc_z[k] * snh_1308[k];
    }

#pragma omp simd aligned(t_1748, t_1749, t_1750, t_1751, pc_x, pc_y, smh_1101, sng0_942, \
                         sng0_944, sng1_942, sng1_944, snh_1311, snh_1314, snh_1316, \
                         snh_1317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1748[k] = f_8 * sng0_942[k]
                    - f_9 * sng1_942[k]
                    + f_3 * pc_x[k] * snh_1314[k];

        t_1749[k] = f_13 * smh_1101[k]
                    + f_3 * pc_y[k] * snh_1311[k];

        t_1750[k] = f_8 * sng0_944[k]
                    - f_9 * sng1_944[k]
                    + f_3 * pc_x[k] * snh_1316[k];

        t_1751[k] = f_3 * pc_x[k] * snh_1317[k];
    }

#pragma omp simd aligned(t_1752, t_1753, t_1754, t_1755, t_1756, pc_x, snh_1318, snh_1319, \
                         snh_1320, snh_1321, snh_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1752[k] = f_3 * pc_x[k] * snh_1318[k];

        t_1753[k] = f_3 * pc_x[k] * snh_1319[k];

        t_1754[k] = f_3 * pc_x[k] * snh_1320[k];

        t_1755[k] = f_3 * pc_x[k] * snh_1321[k];

        t_1756[k] = f_3 * pc_x[k] * snh_1322[k];
    }

#pragma omp simd aligned(t_1757, t_1758, t_1759, pc_y, pc_z, smh_1086, smh_1107, smh_1109, \
                         sng0_940, sng0_942, sng1_940, sng1_942, snh_1317, \
                         snh_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1757[k] = f_13 * smh_1107[k]
                    + f_1 * sng0_940[k]
                    - f_2 * sng1_940[k]
                    + f_3 * pc_y[k] * snh_1317[k];

        t_1758[k] = f_17 * smh_1086[k]
                    + f_3 * pc_z[k] * snh_1317[k];

        t_1759[k] = f_13 * smh_1109[k]
                    + f_4 * sng0_942[k]
                    - f_5 * sng1_942[k]
                    + f_3 * pc_y[k] * snh_1319[k];
    }

#pragma omp simd aligned(t_1760, t_1761, t_1762, pc_y, smh_1110, smh_1111, smh_1112, sng0_943, \
                         sng0_944, sng1_943, sng1_944, snh_1320, snh_1321, \
                         snh_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1760[k] = f_13 * smh_1110[k]
                    + f_6 * sng0_943[k]
                    - f_7 * sng1_943[k]
                    + f_3 * pc_y[k] * snh_1320[k];

        t_1761[k] = f_13 * smh_1111[k]
                    + f_8 * sng0_944[k]
                    - f_9 * sng1_944[k]
                    + f_3 * pc_y[k] * snh_1321[k];

        t_1762[k] = f_13 * smh_1112[k]
                    + f_3 * pc_y[k] * snh_1322[k];
    }

#pragma omp simd aligned(t_1763, t_1764, t_1765, t_1766, pc_x, pc_y, pc_z, smh_1091, smh_1092, \
                         smh_1113, sng0_944, sng0_945, sng1_944, sng1_945, snh_1322, \
                         snh_1323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1763[k] = f_17 * smh_1091[k]
                    + f_1 * sng0_944[k]
                    - f_2 * sng1_944[k]
                    + f_3 * pc_z[k] * snh_1322[k];

        t_1764[k] = f_1 * sng0_945[k]
                    - f_2 * sng1_945[k]
                    + f_3 * pc_x[k] * snh_1323[k];

        t_1765[k] = f_12 * smh_1113[k]
                    + f_3 * pc_y[k] * snh_1323[k];

        t_1766[k] = f_16 * smh_1092[k]
                    + f_3 * pc_z[k] * snh_1323[k];
    }

#pragma omp simd aligned(t_1767, t_1768, t_1769, pc_x, pc_y, smh_1115, sng0_948, sng0_950, \
                         sng1_948, sng1_950, snh_1325, snh_1326, \
                         snh_1328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1767[k] = f_4 * sng0_948[k]
                    - f_5 * sng1_948[k]
                    + f_3 * pc_x[k] * snh_1326[k];

        t_1768[k] = f_12 * smh_1115[k]
                    + f_3 * pc_y[k] * snh_1325[k];

        t_1769[k] = f_4 * sng0_950[k]
                    - f_5 * sng1_950[k]
                    + f_3 * pc_x[k] * snh_1328[k];
    }

#pragma omp simd aligned(t_1770, t_1771, t_1772, pc_x, pc_y, pc_z, smh_1095, smh_1118, \
                         sng0_951, sng1_951, snh_1326, snh_1328, \
                         snh_1329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1770[k] = f_6 * sng0_951[k]
                    - f_7 * sng1_951[k]
                    + f_3 * pc_x[k] * snh_1329[k];

        t_1771[k] = f_16 * smh_1095[k]
                    + f_3 * pc_z[k] * snh_1326[k];

        t_1772[k] = f_12 * smh_1118[k]
                    + f_3 * pc_y[k] * snh_1328[k];
    }

#pragma omp simd aligned(t_1773, t_1774, t_1775, pc_x, pc_z, smh_1098, sng0_954, sng0_955, \
                         sng1_954, sng1_955, snh_1329, snh_1332, \
                         snh_1333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1773[k] = f_6 * sng0_954[k]
                    - f_7 * sng1_954[k]
                    + f_3 * pc_x[k] * snh_1332[k];

        t_1774[k] = f_8 * sng0_955[k]
                    - f_9 * sng1_955[k]
                    + f_3 * pc_x[k] * snh_1333[k];

        t_1775[k] = f_16 * smh_1098[k]
                    + f_3 * pc_z[k] * snh_1329[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, t_1779, pc_x, pc_y, smh_1122, sng0_957, \
                         sng0_959, sng1_957, sng1_959, snh_1332, snh_1335, snh_1337, \
                         snh_1338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = f_8 * sng0_957[k]
                    - f_9 * sng1_957[k]
                    + f_3 * pc_x[k] * snh_1335[k];

        t_1777[k] = f_12 * smh_1122[k]
                    + f_3 * pc_y[k] * snh_1332[k];

        t_1778[k] = f_8 * sng0_959[k]
                    - f_9 * sng1_959[k]
                    + f_3 * pc_x[k] * snh_1337[k];

        t_1779[k] = f_3 * pc_x[k] * snh_1338[k];
    }

#pragma omp simd aligned(t_1780, t_1781, t_1782, t_1783, t_1784, pc_x, snh_1339, snh_1340, \
                         snh_1341, snh_1342, snh_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1780[k] = f_3 * pc_x[k] * snh_1339[k];

        t_1781[k] = f_3 * pc_x[k] * snh_1340[k];

        t_1782[k] = f_3 * pc_x[k] * snh_1341[k];

        t_1783[k] = f_3 * pc_x[k] * snh_1342[k];

        t_1784[k] = f_3 * pc_x[k] * snh_1343[k];
    }

#pragma omp simd aligned(t_1785, t_1786, t_1787, pc_y, pc_z, smh_1107, smh_1128, smh_1130, \
                         sng0_955, sng0_957, sng1_955, sng1_957, snh_1338, \
                         snh_1340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1785[k] = f_12 * smh_1128[k]
                    + f_1 * sng0_955[k]
                    - f_2 * sng1_955[k]
                    + f_3 * pc_y[k] * snh_1338[k];

        t_1786[k] = f_16 * smh_1107[k]
                    + f_3 * pc_z[k] * snh_1338[k];

        t_1787[k] = f_12 * smh_1130[k]
                    + f_4 * sng0_957[k]
                    - f_5 * sng1_957[k]
                    + f_3 * pc_y[k] * snh_1340[k];
    }

#pragma omp simd aligned(t_1788, t_1789, t_1790, pc_y, smh_1131, smh_1132, smh_1133, sng0_958, \
                         sng0_959, sng1_958, sng1_959, snh_1341, snh_1342, \
                         snh_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1788[k] = f_12 * smh_1131[k]
                    + f_6 * sng0_958[k]
                    - f_7 * sng1_958[k]
                    + f_3 * pc_y[k] * snh_1341[k];

        t_1789[k] = f_12 * smh_1132[k]
                    + f_8 * sng0_959[k]
                    - f_9 * sng1_959[k]
                    + f_3 * pc_y[k] * snh_1342[k];

        t_1790[k] = f_12 * smh_1133[k]
                    + f_3 * pc_y[k] * snh_1343[k];
    }

#pragma omp simd aligned(t_1791, t_1792, t_1793, t_1794, pb_y, pc_y, pc_z, smi0_1512, \
                         smh_1112, smh_1113, smh_1134, smi1_1512, sng0_959, sng1_959, \
                         snh_1343, snh_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1791[k] = f_16 * smh_1112[k]
                    + f_1 * sng0_959[k]
                    - f_2 * sng1_959[k]
                    + f_3 * pc_z[k] * snh_1343[k];

        t_1792[k] = pb_y[k] * smi0_1512[k]
                    - f_10 * pc_y[k] * smi1_1512[k];

        t_1793[k] = f_11 * smh_1134[k]
                    + f_3 * pc_y[k] * snh_1344[k];

        t_1794[k] = f_15 * smh_1113[k]
                    + f_3 * pc_z[k] * snh_1344[k];
    }

#pragma omp simd aligned(t_1795, t_1796, t_1797, pb_y, pc_x, pc_y, smi0_1517, smh_1136, \
                         smi1_1517, sng0_963, sng1_963, snh_1346, \
                         snh_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1795[k] = f_4 * sng0_963[k]
                    - f_5 * sng1_963[k]
                    + f_3 * pc_x[k] * snh_1347[k];

        t_1796[k] = f_11 * smh_1136[k]
                    + f_3 * pc_y[k] * snh_1346[k];

        t_1797[k] = pb_y[k] * smi0_1517[k]
                    - f_10 * pc_y[k] * smi1_1517[k];
    }

#pragma omp simd aligned(t_1798, t_1799, t_1800, pc_x, pc_y, pc_z, smh_1116, smh_1139, \
                         sng0_966, sng1_966, snh_1347, snh_1349, \
                         snh_1350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1798[k] = f_6 * sng0_966[k]
                    - f_7 * sng1_966[k]
                    + f_3 * pc_x[k] * snh_1350[k];

        t_1799[k] = f_15 * smh_1116[k]
                    + f_3 * pc_z[k] * snh_1347[k];

        t_1800[k] = f_11 * smh_1139[k]
                    + f_3 * pc_y[k] * snh_1349[k];
    }

#pragma omp simd aligned(t_1801, t_1802, t_1803, pb_y, pc_x, pc_y, pc_z, smi0_1521, smh_1119, \
                         smi1_1521, sng0_970, sng1_970, snh_1350, \
                         snh_1354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1801[k] = pb_y[k] * smi0_1521[k]
                    - f_10 * pc_y[k] * smi1_1521[k];

        t_1802[k] = f_8 * sng0_970[k]
                    - f_9 * sng1_970[k]
                    + f_3 * pc_x[k] * snh_1354[k];

        t_1803[k] = f_15 * smh_1119[k]
                    + f_3 * pc_z[k] * snh_1350[k];
    }

#pragma omp simd aligned(t_1804, t_1805, t_1806, t_1807, pb_y, pc_x, pc_y, smi0_1526, \
                         smh_1143, smi1_1526, sng0_972, sng1_972, snh_1353, snh_1356, \
                         snh_1359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1804[k] = f_8 * sng0_972[k]
                    - f_9 * sng1_972[k]
                    + f_3 * pc_x[k] * snh_1356[k];

        t_1805[k] = f_11 * smh_1143[k]
                    + f_3 * pc_y[k] * snh_1353[k];

        t_1806[k] = pb_y[k] * smi0_1526[k]
                    - f_10 * pc_y[k] * smi1_1526[k];

        t_1807[k] = f_3 * pc_x[k] * snh_1359[k];
    }

#pragma omp simd aligned(t_1808, t_1809, t_1810, t_1811, t_1812, pc_x, snh_1360, snh_1361, \
                         snh_1362, snh_1363, snh_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1808[k] = f_3 * pc_x[k] * snh_1360[k];

        t_1809[k] = f_3 * pc_x[k] * snh_1361[k];

        t_1810[k] = f_3 * pc_x[k] * snh_1362[k];

        t_1811[k] = f_3 * pc_x[k] * snh_1363[k];

        t_1812[k] = f_3 * pc_x[k] * snh_1364[k];
    }

#pragma omp simd aligned(t_1813, t_1814, t_1815, pb_y, pc_y, pc_z, smi0_1533, smi0_1535, \
                         smh_1128, smh_1149, smh_1151, smi1_1533, smi1_1535, \
                         snh_1359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1813[k] = pb_y[k] * smi0_1533[k]
                    + f_18 * smh_1149[k]
                    - f_10 * pc_y[k] * smi1_1533[k];

        t_1814[k] = f_15 * smh_1128[k]
                    + f_3 * pc_z[k] * snh_1359[k];

        t_1815[k] = pb_y[k] * smi0_1535[k]
                    + f_14 * smh_1151[k]
                    - f_10 * pc_y[k] * smi1_1535[k];
    }

#pragma omp simd aligned(t_1816, t_1817, t_1818, t_1819, pb_y, pc_y, smi0_1536, smi0_1537, \
                         smi0_1539, smh_1152, smh_1153, smh_1154, smi1_1536, smi1_1537, \
                         smi1_1539, snh_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1816[k] = pb_y[k] * smi0_1536[k]
                    + f_13 * smh_1152[k]
                    - f_10 * pc_y[k] * smi1_1536[k];

        t_1817[k] = pb_y[k] * smi0_1537[k]
                    + f_12 * smh_1153[k]
                    - f_10 * pc_y[k] * smi1_1537[k];

        t_1818[k] = f_11 * smh_1154[k]
                    + f_3 * pc_y[k] * snh_1364[k];

        t_1819[k] = pb_y[k] * smi0_1539[k]
                    - f_10 * pc_y[k] * smi1_1539[k];
    }

#pragma omp simd aligned(t_1820, t_1821, t_1822, t_1823, t_1824, pc_x, pc_y, pc_z, smh_1134, \
                         sng0_975, sng0_978, sng1_975, sng1_978, snh_1365, snh_1367, \
                         snh_1368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1820[k] = f_1 * sng0_975[k]
                    - f_2 * sng1_975[k]
                    + f_3 * pc_x[k] * snh_1365[k];

        t_1821[k] = f_3 * pc_y[k] * snh_1365[k];

        t_1822[k] = f_0 * smh_1134[k]
                    + f_3 * pc_z[k] * snh_1365[k];

        t_1823[k] = f_4 * sng0_978[k]
                    - f_5 * sng1_978[k]
                    + f_3 * pc_x[k] * snh_1368[k];

        t_1824[k] = f_3 * pc_y[k] * snh_1367[k];
    }

#pragma omp simd aligned(t_1825, t_1826, t_1827, t_1828, pc_x, pc_y, pc_z, smh_1137, sng0_980, \
                         sng0_981, sng1_980, sng1_981, snh_1368, snh_1370, \
                         snh_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1825[k] = f_4 * sng0_980[k]
                    - f_5 * sng1_980[k]
                    + f_3 * pc_x[k] * snh_1370[k];

        t_1826[k] = f_6 * sng0_981[k]
                    - f_7 * sng1_981[k]
                    + f_3 * pc_x[k] * snh_1371[k];

        t_1827[k] = f_0 * smh_1137[k]
                    + f_3 * pc_z[k] * snh_1368[k];

        t_1828[k] = f_3 * pc_y[k] * snh_1370[k];
    }

#pragma omp simd aligned(t_1829, t_1830, t_1831, pc_x, pc_z, smh_1140, sng0_984, sng0_985, \
                         sng1_984, sng1_985, snh_1371, snh_1374, \
                         snh_1375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1829[k] = f_6 * sng0_984[k]
                    - f_7 * sng1_984[k]
                    + f_3 * pc_x[k] * snh_1374[k];

        t_1830[k] = f_8 * sng0_985[k]
                    - f_9 * sng1_985[k]
                    + f_3 * pc_x[k] * snh_1375[k];

        t_1831[k] = f_0 * smh_1140[k]
                    + f_3 * pc_z[k] * snh_1371[k];
    }

#pragma omp simd aligned(t_1832, t_1833, t_1834, t_1835, t_1836, pc_x, pc_y, sng0_987, \
                         sng0_989, sng1_987, sng1_989, snh_1374, snh_1377, snh_1379, snh_1380, \
                         snh_1381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1832[k] = f_8 * sng0_987[k]
                    - f_9 * sng1_987[k]
                    + f_3 * pc_x[k] * snh_1377[k];

        t_1833[k] = f_3 * pc_y[k] * snh_1374[k];

        t_1834[k] = f_8 * sng0_989[k]
                    - f_9 * sng1_989[k]
                    + f_3 * pc_x[k] * snh_1379[k];

        t_1835[k] = f_3 * pc_x[k] * snh_1380[k];

        t_1836[k] = f_3 * pc_x[k] * snh_1381[k];
    }

#pragma omp simd aligned(t_1837, t_1838, t_1839, t_1840, t_1841, pc_x, pc_y, sng0_985, \
                         sng1_985, snh_1380, snh_1382, snh_1383, snh_1384, \
                         snh_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1837[k] = f_3 * pc_x[k] * snh_1382[k];

        t_1838[k] = f_3 * pc_x[k] * snh_1383[k];

        t_1839[k] = f_3 * pc_x[k] * snh_1384[k];

        t_1840[k] = f_3 * pc_x[k] * snh_1385[k];

        t_1841[k] = f_1 * sng0_985[k]
                    - f_2 * sng1_985[k]
                    + f_3 * pc_y[k] * snh_1380[k];
    }

#pragma omp simd aligned(t_1842, t_1843, t_1844, pc_y, pc_z, smh_1149, sng0_987, sng0_988, \
                         sng1_987, sng1_988, snh_1380, snh_1382, \
                         snh_1383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1842[k] = f_0 * smh_1149[k]
                    + f_3 * pc_z[k] * snh_1380[k];

        t_1843[k] = f_4 * sng0_987[k]
                    - f_5 * sng1_987[k]
                    + f_3 * pc_y[k] * snh_1382[k];

        t_1844[k] = f_6 * sng0_988[k]
                    - f_7 * sng1_988[k]
                    + f_3 * pc_y[k] * snh_1383[k];
    }

#pragma omp simd aligned(t_1845, t_1846, t_1847, pc_y, pc_z, smh_1154, sng0_989, sng1_989, \
                         snh_1384, snh_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1845[k] = f_8 * sng0_989[k]
                    - f_9 * sng1_989[k]
                    + f_3 * pc_y[k] * snh_1384[k];

        t_1846[k] = f_3 * pc_y[k] * snh_1385[k];

        t_1847[k] = f_0 * smh_1154[k]
                    + f_1 * sng0_989[k]
                    - f_2 * sng1_989[k]
                    + f_3 * pc_z[k] * snh_1385[k];
    }
}

auto
compute_prim_sni_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t smi0, const size_t smh,
                                                   const size_t smi1, const size_t sng0,
                                                   const size_t sng1, const size_t snh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sni_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, smi0, smh,
                                                              smi1, sng0, sng1, snh, ncols,
                                                              gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, smi0,
                                                               smh, smi1, sng0, sng1, snh,
                                                               ncols, gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, smi0,
                                                               smh, smi1, sng0, sng1, snh,
                                                               ncols, gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece12(buffer, target, pb, pc, smi0,
                                                               smh, smi1, snh, ncols, gamma, p,
                                                               q);

    compute_prim_sni_three_center_electron_repulsion_0_piece13(buffer, target, pb, pc, smi0,
                                                               smh, smi1, sng0, sng1, snh,
                                                               ncols, gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece14(buffer, target, pc, smh, sng0,
                                                               sng1, snh, ncols, gamma, p, q);

    compute_prim_sni_three_center_electron_repulsion_0_piece15(buffer, target, pb, pc, smi0,
                                                               smh, smi1, sng0, sng1, snh,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
