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


#include "SimdThreeCenterElectronRepulsionVrrRecSMK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_0 = buffer.data(slk0 + 0);
    const auto *slk0_3 = buffer.data(slk0 + 3);
    const auto *slk0_5 = buffer.data(slk0 + 5);
    const auto *slk0_6 = buffer.data(slk0 + 6);
    const auto *slk0_9 = buffer.data(slk0 + 9);
    const auto *slk0_10 = buffer.data(slk0 + 10);
    const auto *slk0_12 = buffer.data(slk0 + 12);
    const auto *slk0_14 = buffer.data(slk0 + 14);
    const auto *slk0_15 = buffer.data(slk0 + 15);
    const auto *slk0_17 = buffer.data(slk0 + 17);
    const auto *slk0_18 = buffer.data(slk0 + 18);
    const auto *slk0_20 = buffer.data(slk0 + 20);
    const auto *slk0_28 = buffer.data(slk0 + 28);
    const auto *slk0_35 = buffer.data(slk0 + 35);

    const auto *sli_0 = buffer.data(sli + 0);
    const auto *sli_1 = buffer.data(sli + 1);
    const auto *sli_2 = buffer.data(sli + 2);
    const auto *sli_3 = buffer.data(sli + 3);
    const auto *sli_5 = buffer.data(sli + 5);
    const auto *sli_6 = buffer.data(sli + 6);
    const auto *sli_7 = buffer.data(sli + 7);
    const auto *sli_8 = buffer.data(sli + 8);
    const auto *sli_9 = buffer.data(sli + 9);
    const auto *sli_10 = buffer.data(sli + 10);
    const auto *sli_11 = buffer.data(sli + 11);
    const auto *sli_12 = buffer.data(sli + 12);
    const auto *sli_13 = buffer.data(sli + 13);
    const auto *sli_14 = buffer.data(sli + 14);
    const auto *sli_15 = buffer.data(sli + 15);
    const auto *sli_17 = buffer.data(sli + 17);
    const auto *sli_18 = buffer.data(sli + 18);
    const auto *sli_20 = buffer.data(sli + 20);
    const auto *sli_21 = buffer.data(sli + 21);
    const auto *sli_22 = buffer.data(sli + 22);
    const auto *sli_23 = buffer.data(sli + 23);
    const auto *sli_24 = buffer.data(sli + 24);
    const auto *sli_25 = buffer.data(sli + 25);
    const auto *sli_26 = buffer.data(sli + 26);
    const auto *sli_27 = buffer.data(sli + 27);
    const auto *sli_28 = buffer.data(sli + 28);
    const auto *sli_30 = buffer.data(sli + 30);
    const auto *sli_33 = buffer.data(sli + 33);
    const auto *sli_37 = buffer.data(sli + 37);
    const auto *sli_49 = buffer.data(sli + 49);
    const auto *sli_50 = buffer.data(sli + 50);
    const auto *sli_51 = buffer.data(sli + 51);
    const auto *sli_52 = buffer.data(sli + 52);
    const auto *sli_53 = buffer.data(sli + 53);
    const auto *sli_54 = buffer.data(sli + 54);
    const auto *sli_55 = buffer.data(sli + 55);
    const auto *sli_77 = buffer.data(sli + 77);
    const auto *sli_78 = buffer.data(sli + 78);
    const auto *sli_79 = buffer.data(sli + 79);
    const auto *sli_80 = buffer.data(sli + 80);
    const auto *sli_81 = buffer.data(sli + 81);
    const auto *sli_82 = buffer.data(sli + 82);
    const auto *sli_83 = buffer.data(sli + 83);
    const auto *sli_84 = buffer.data(sli + 84);
    const auto *sli_87 = buffer.data(sli + 87);
    const auto *sli_89 = buffer.data(sli + 89);
    const auto *sli_90 = buffer.data(sli + 90);
    const auto *sli_93 = buffer.data(sli + 93);
    const auto *sli_94 = buffer.data(sli + 94);
    const auto *sli_96 = buffer.data(sli + 96);

    const auto *slk1_0 = buffer.data(slk1 + 0);
    const auto *slk1_3 = buffer.data(slk1 + 3);
    const auto *slk1_5 = buffer.data(slk1 + 5);
    const auto *slk1_6 = buffer.data(slk1 + 6);
    const auto *slk1_9 = buffer.data(slk1 + 9);
    const auto *slk1_10 = buffer.data(slk1 + 10);
    const auto *slk1_12 = buffer.data(slk1 + 12);
    const auto *slk1_14 = buffer.data(slk1 + 14);
    const auto *slk1_15 = buffer.data(slk1 + 15);
    const auto *slk1_17 = buffer.data(slk1 + 17);
    const auto *slk1_18 = buffer.data(slk1 + 18);
    const auto *slk1_20 = buffer.data(slk1 + 20);
    const auto *slk1_28 = buffer.data(slk1 + 28);
    const auto *slk1_35 = buffer.data(slk1 + 35);

    const auto *smh0_0 = buffer.data(smh0 + 0);
    const auto *smh0_3 = buffer.data(smh0 + 3);
    const auto *smh0_5 = buffer.data(smh0 + 5);
    const auto *smh0_6 = buffer.data(smh0 + 6);
    const auto *smh0_9 = buffer.data(smh0 + 9);
    const auto *smh0_10 = buffer.data(smh0 + 10);
    const auto *smh0_12 = buffer.data(smh0 + 12);
    const auto *smh0_14 = buffer.data(smh0 + 14);
    const auto *smh0_15 = buffer.data(smh0 + 15);
    const auto *smh0_17 = buffer.data(smh0 + 17);
    const auto *smh0_18 = buffer.data(smh0 + 18);
    const auto *smh0_19 = buffer.data(smh0 + 19);
    const auto *smh0_20 = buffer.data(smh0 + 20);
    const auto *smh0_36 = buffer.data(smh0 + 36);
    const auto *smh0_38 = buffer.data(smh0 + 38);
    const auto *smh0_39 = buffer.data(smh0 + 39);
    const auto *smh0_40 = buffer.data(smh0 + 40);
    const auto *smh0_41 = buffer.data(smh0 + 41);
    const auto *smh0_59 = buffer.data(smh0 + 59);
    const auto *smh0_60 = buffer.data(smh0 + 60);
    const auto *smh0_61 = buffer.data(smh0 + 61);
    const auto *smh0_62 = buffer.data(smh0 + 62);
    const auto *smh0_63 = buffer.data(smh0 + 63);
    const auto *smh0_66 = buffer.data(smh0 + 66);
    const auto *smh0_68 = buffer.data(smh0 + 68);
    const auto *smh0_69 = buffer.data(smh0 + 69);
    const auto *smh0_72 = buffer.data(smh0 + 72);
    const auto *smh0_73 = buffer.data(smh0 + 73);
    const auto *smh0_75 = buffer.data(smh0 + 75);

    const auto *smh1_0 = buffer.data(smh1 + 0);
    const auto *smh1_3 = buffer.data(smh1 + 3);
    const auto *smh1_5 = buffer.data(smh1 + 5);
    const auto *smh1_6 = buffer.data(smh1 + 6);
    const auto *smh1_9 = buffer.data(smh1 + 9);
    const auto *smh1_10 = buffer.data(smh1 + 10);
    const auto *smh1_12 = buffer.data(smh1 + 12);
    const auto *smh1_14 = buffer.data(smh1 + 14);
    const auto *smh1_15 = buffer.data(smh1 + 15);
    const auto *smh1_17 = buffer.data(smh1 + 17);
    const auto *smh1_18 = buffer.data(smh1 + 18);
    const auto *smh1_19 = buffer.data(smh1 + 19);
    const auto *smh1_20 = buffer.data(smh1 + 20);
    const auto *smh1_36 = buffer.data(smh1 + 36);
    const auto *smh1_38 = buffer.data(smh1 + 38);
    const auto *smh1_39 = buffer.data(smh1 + 39);
    const auto *smh1_40 = buffer.data(smh1 + 40);
    const auto *smh1_41 = buffer.data(smh1 + 41);
    const auto *smh1_59 = buffer.data(smh1 + 59);
    const auto *smh1_60 = buffer.data(smh1 + 60);
    const auto *smh1_61 = buffer.data(smh1 + 61);
    const auto *smh1_62 = buffer.data(smh1 + 62);
    const auto *smh1_63 = buffer.data(smh1 + 63);
    const auto *smh1_66 = buffer.data(smh1 + 66);
    const auto *smh1_68 = buffer.data(smh1 + 68);
    const auto *smh1_69 = buffer.data(smh1 + 69);
    const auto *smh1_72 = buffer.data(smh1 + 72);
    const auto *smh1_73 = buffer.data(smh1 + 73);
    const auto *smh1_75 = buffer.data(smh1 + 75);

    const auto *smi_0 = buffer.data(smi + 0);
    const auto *smi_2 = buffer.data(smi + 2);
    const auto *smi_3 = buffer.data(smi + 3);
    const auto *smi_5 = buffer.data(smi + 5);
    const auto *smi_6 = buffer.data(smi + 6);
    const auto *smi_9 = buffer.data(smi + 9);
    const auto *smi_10 = buffer.data(smi + 10);
    const auto *smi_12 = buffer.data(smi + 12);
    const auto *smi_14 = buffer.data(smi + 14);
    const auto *smi_15 = buffer.data(smi + 15);
    const auto *smi_17 = buffer.data(smi + 17);
    const auto *smi_18 = buffer.data(smi + 18);
    const auto *smi_20 = buffer.data(smi + 20);
    const auto *smi_21 = buffer.data(smi + 21);
    const auto *smi_22 = buffer.data(smi + 22);
    const auto *smi_23 = buffer.data(smi + 23);
    const auto *smi_24 = buffer.data(smi + 24);
    const auto *smi_25 = buffer.data(smi + 25);
    const auto *smi_26 = buffer.data(smi + 26);
    const auto *smi_27 = buffer.data(smi + 27);
    const auto *smi_28 = buffer.data(smi + 28);
    const auto *smi_30 = buffer.data(smi + 30);
    const auto *smi_31 = buffer.data(smi + 31);
    const auto *smi_33 = buffer.data(smi + 33);
    const auto *smi_34 = buffer.data(smi + 34);
    const auto *smi_37 = buffer.data(smi + 37);
    const auto *smi_38 = buffer.data(smi + 38);
    const auto *smi_42 = buffer.data(smi + 42);
    const auto *smi_49 = buffer.data(smi + 49);
    const auto *smi_50 = buffer.data(smi + 50);
    const auto *smi_51 = buffer.data(smi + 51);
    const auto *smi_52 = buffer.data(smi + 52);
    const auto *smi_53 = buffer.data(smi + 53);
    const auto *smi_54 = buffer.data(smi + 54);
    const auto *smi_55 = buffer.data(smi + 55);
    const auto *smi_56 = buffer.data(smi + 56);
    const auto *smi_58 = buffer.data(smi + 58);
    const auto *smi_59 = buffer.data(smi + 59);
    const auto *smi_61 = buffer.data(smi + 61);
    const auto *smi_62 = buffer.data(smi + 62);
    const auto *smi_65 = buffer.data(smi + 65);
    const auto *smi_66 = buffer.data(smi + 66);
    const auto *smi_70 = buffer.data(smi + 70);
    const auto *smi_77 = buffer.data(smi + 77);
    const auto *smi_78 = buffer.data(smi + 78);
    const auto *smi_79 = buffer.data(smi + 79);
    const auto *smi_80 = buffer.data(smi + 80);
    const auto *smi_81 = buffer.data(smi + 81);
    const auto *smi_82 = buffer.data(smi + 82);
    const auto *smi_83 = buffer.data(smi + 83);
    const auto *smi_84 = buffer.data(smi + 84);
    const auto *smi_86 = buffer.data(smi + 86);
    const auto *smi_87 = buffer.data(smi + 87);
    const auto *smi_89 = buffer.data(smi + 89);
    const auto *smi_90 = buffer.data(smi + 90);
    const auto *smi_93 = buffer.data(smi + 93);
    const auto *smi_94 = buffer.data(smi + 94);
    const auto *smi_96 = buffer.data(smi + 96);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sli_0, sli_3, smh0_0, smh0_3, \
                         smh1_0, smh1_3, smi_0, smi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sli_0[k]
                 + f_1 * smh0_0[k]
                 - f_2 * smh1_0[k]
                 + f_3 * pc_x[k] * smi_0[k];

        t_1[k] = f_3 * pc_y[k] * smi_0[k];

        t_2[k] = f_3 * pc_z[k] * smi_0[k];

        t_3[k] = f_0 * sli_3[k]
                 + f_4 * smh0_3[k]
                 - f_5 * smh1_3[k]
                 + f_3 * pc_x[k] * smi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sli_5, sli_6, smh0_5, smh0_6, smh1_5, \
                         smh1_6, smi_2, smi_5, smi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * smi_2[k];

        t_5[k] = f_0 * sli_5[k]
                 + f_4 * smh0_5[k]
                 - f_5 * smh1_5[k]
                 + f_3 * pc_x[k] * smi_5[k];

        t_6[k] = f_0 * sli_6[k]
                 + f_6 * smh0_6[k]
                 - f_7 * smh1_6[k]
                 + f_3 * pc_x[k] * smi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sli_9, smh0_9, smh1_9, smi_3, smi_5, \
                         smi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * smi_3[k];

        t_8[k] = f_3 * pc_y[k] * smi_5[k];

        t_9[k] = f_0 * sli_9[k]
                 + f_6 * smh0_9[k]
                 - f_7 * smh1_9[k]
                 + f_3 * pc_x[k] * smi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sli_10, sli_12, smh0_10, smh0_12, \
                         smh1_10, smh1_12, smi_6, smi_10, smi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sli_10[k]
                  + f_8 * smh0_10[k]
                  - f_9 * smh1_10[k]
                  + f_3 * pc_x[k] * smi_10[k];

        t_11[k] = f_3 * pc_z[k] * smi_6[k];

        t_12[k] = f_0 * sli_12[k]
                  + f_8 * smh0_12[k]
                  - f_9 * smh1_12[k]
                  + f_3 * pc_x[k] * smi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sli_14, sli_15, smh0_14, smh0_15, \
                         smh1_14, smh1_15, smi_9, smi_14, smi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * smi_9[k];

        t_14[k] = f_0 * sli_14[k]
                  + f_8 * smh0_14[k]
                  - f_9 * smh1_14[k]
                  + f_3 * pc_x[k] * smi_14[k];

        t_15[k] = f_0 * sli_15[k]
                  + f_10 * smh0_15[k]
                  - f_11 * smh1_15[k]
                  + f_3 * pc_x[k] * smi_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sli_17, sli_18, smh0_17, smh0_18, \
                         smh1_17, smh1_18, smi_10, smi_17, smi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * smi_10[k];

        t_17[k] = f_0 * sli_17[k]
                  + f_10 * smh0_17[k]
                  - f_11 * smh1_17[k]
                  + f_3 * pc_x[k] * smi_17[k];

        t_18[k] = f_0 * sli_18[k]
                  + f_10 * smh0_18[k]
                  - f_11 * smh1_18[k]
                  + f_3 * pc_x[k] * smi_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, sli_20, sli_21, sli_22, smh0_20, \
                         smh1_20, smi_14, smi_20, smi_21, smi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * smi_14[k];

        t_20[k] = f_0 * sli_20[k]
                  + f_10 * smh0_20[k]
                  - f_11 * smh1_20[k]
                  + f_3 * pc_x[k] * smi_20[k];

        t_21[k] = f_0 * sli_21[k]
                  + f_3 * pc_x[k] * smi_21[k];

        t_22[k] = f_0 * sli_22[k]
                  + f_3 * pc_x[k] * smi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, sli_23, sli_24, sli_25, sli_26, \
                         sli_27, smi_23, smi_24, smi_25, smi_26, \
                         smi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * sli_23[k]
                  + f_3 * pc_x[k] * smi_23[k];

        t_24[k] = f_0 * sli_24[k]
                  + f_3 * pc_x[k] * smi_24[k];

        t_25[k] = f_0 * sli_25[k]
                  + f_3 * pc_x[k] * smi_25[k];

        t_26[k] = f_0 * sli_26[k]
                  + f_3 * pc_x[k] * smi_26[k];

        t_27[k] = f_0 * sli_27[k]
                  + f_3 * pc_x[k] * smi_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, smh0_15, smh0_17, smh0_18, \
                         smh1_15, smh1_17, smh1_18, smi_21, smi_23, \
                         smi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * smh0_15[k]
                  - f_2 * smh1_15[k]
                  + f_3 * pc_y[k] * smi_21[k];

        t_29[k] = f_3 * pc_z[k] * smi_21[k];

        t_30[k] = f_4 * smh0_17[k]
                  - f_5 * smh1_17[k]
                  + f_3 * pc_y[k] * smi_23[k];

        t_31[k] = f_6 * smh0_18[k]
                  - f_7 * smh1_18[k]
                  + f_3 * pc_y[k] * smi_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, smh0_19, smh0_20, smh1_19, \
                         smh1_20, smi_25, smi_26, smi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * smh0_19[k]
                  - f_9 * smh1_19[k]
                  + f_3 * pc_y[k] * smi_25[k];

        t_33[k] = f_10 * smh0_20[k]
                  - f_11 * smh1_20[k]
                  + f_3 * pc_y[k] * smi_26[k];

        t_34[k] = f_3 * pc_y[k] * smi_27[k];

        t_35[k] = f_1 * smh0_20[k]
                  - f_2 * smh1_20[k]
                  + f_3 * pc_z[k] * smi_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, slk0_0, slk0_3, sli_0, \
                         sli_1, slk1_0, slk1_3, smi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * slk0_0[k]
                  - f_12 * pc_y[k] * slk1_0[k];

        t_37[k] = f_13 * sli_0[k]
                  + f_3 * pc_y[k] * smi_28[k];

        t_38[k] = f_3 * pc_z[k] * smi_28[k];

        t_39[k] = pb_y[k] * slk0_3[k]
                  + f_14 * sli_1[k]
                  - f_12 * pc_y[k] * slk1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_y, pc_z, slk0_5, slk0_6, sli_2, \
                         sli_3, slk1_5, slk1_6, smi_30, smi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * sli_2[k]
                  + f_3 * pc_y[k] * smi_30[k];

        t_41[k] = pb_y[k] * slk0_5[k]
                  - f_12 * pc_y[k] * slk1_5[k];

        t_42[k] = pb_y[k] * slk0_6[k]
                  + f_15 * sli_3[k]
                  - f_12 * pc_y[k] * slk1_6[k];

        t_43[k] = f_3 * pc_z[k] * smi_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, slk0_9, slk0_10, sli_5, \
                         sli_6, slk1_9, slk1_10, smi_33, smi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * sli_5[k]
                  + f_3 * pc_y[k] * smi_33[k];

        t_45[k] = pb_y[k] * slk0_9[k]
                  - f_12 * pc_y[k] * slk1_9[k];

        t_46[k] = pb_y[k] * slk0_10[k]
                  + f_16 * sli_6[k]
                  - f_12 * pc_y[k] * slk1_10[k];

        t_47[k] = f_3 * pc_z[k] * smi_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, slk0_12, slk0_14, slk0_15, sli_8, \
                         sli_9, sli_10, slk1_12, slk1_14, slk1_15, \
                         smi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * slk0_12[k]
                  + f_14 * sli_8[k]
                  - f_12 * pc_y[k] * slk1_12[k];

        t_49[k] = f_13 * sli_9[k]
                  + f_3 * pc_y[k] * smi_37[k];

        t_50[k] = pb_y[k] * slk0_14[k]
                  - f_12 * pc_y[k] * slk1_14[k];

        t_51[k] = pb_y[k] * slk0_15[k]
                  + f_17 * sli_10[k]
                  - f_12 * pc_y[k] * slk1_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, pc_z, slk0_17, slk0_18, sli_12, \
                         sli_13, sli_14, slk1_17, slk1_18, smi_38, \
                         smi_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * smi_38[k];

        t_53[k] = pb_y[k] * slk0_17[k]
                  + f_15 * sli_12[k]
                  - f_12 * pc_y[k] * slk1_17[k];

        t_54[k] = pb_y[k] * slk0_18[k]
                  + f_14 * sli_13[k]
                  - f_12 * pc_y[k] * slk1_18[k];

        t_55[k] = f_13 * sli_14[k]
                  + f_3 * pc_y[k] * smi_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pc_x, pc_y, slk0_20, sli_49, sli_50, \
                         sli_51, slk1_20, smi_49, smi_50, smi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * slk0_20[k]
                  - f_12 * pc_y[k] * slk1_20[k];

        t_57[k] = f_18 * sli_49[k]
                  + f_3 * pc_x[k] * smi_49[k];

        t_58[k] = f_18 * sli_50[k]
                  + f_3 * pc_x[k] * smi_50[k];

        t_59[k] = f_18 * sli_51[k]
                  + f_3 * pc_x[k] * smi_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, sli_52, sli_53, sli_54, sli_55, smi_52, \
                         smi_53, smi_54, smi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_18 * sli_52[k]
                  + f_3 * pc_x[k] * smi_52[k];

        t_61[k] = f_18 * sli_53[k]
                  + f_3 * pc_x[k] * smi_53[k];

        t_62[k] = f_18 * sli_54[k]
                  + f_3 * pc_x[k] * smi_54[k];

        t_63[k] = f_18 * sli_55[k]
                  + f_3 * pc_x[k] * smi_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pc_y, pc_z, sli_21, sli_23, smh0_36, smh0_38, \
                         smh1_36, smh1_38, smi_49, smi_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_13 * sli_21[k]
                  + f_1 * smh0_36[k]
                  - f_2 * smh1_36[k]
                  + f_3 * pc_y[k] * smi_49[k];

        t_65[k] = f_3 * pc_z[k] * smi_49[k];

        t_66[k] = f_13 * sli_23[k]
                  + f_4 * smh0_38[k]
                  - f_5 * smh1_38[k]
                  + f_3 * pc_y[k] * smi_51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, sli_24, sli_25, sli_26, smh0_39, smh0_40, \
                         smh0_41, smh1_39, smh1_40, smh1_41, smi_52, smi_53, \
                         smi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * sli_24[k]
                  + f_6 * smh0_39[k]
                  - f_7 * smh1_39[k]
                  + f_3 * pc_y[k] * smi_52[k];

        t_68[k] = f_13 * sli_25[k]
                  + f_8 * smh0_40[k]
                  - f_9 * smh1_40[k]
                  + f_3 * pc_y[k] * smi_53[k];

        t_69[k] = f_13 * sli_26[k]
                  + f_10 * smh0_41[k]
                  - f_11 * smh1_41[k]
                  + f_3 * pc_y[k] * smi_54[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_y, pb_z, pc_y, pc_z, slk0_0, slk0_35, \
                         sli_27, slk1_0, slk1_35, smi_55, smi_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * sli_27[k]
                  + f_3 * pc_y[k] * smi_55[k];

        t_71[k] = pb_y[k] * slk0_35[k]
                  - f_12 * pc_y[k] * slk1_35[k];

        t_72[k] = pb_z[k] * slk0_0[k]
                  - f_12 * pc_z[k] * slk1_0[k];

        t_73[k] = f_3 * pc_y[k] * smi_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_y, pc_z, slk0_3, slk0_5, sli_0, \
                         sli_2, slk1_3, slk1_5, smi_56, smi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * sli_0[k]
                  + f_3 * pc_z[k] * smi_56[k];

        t_75[k] = pb_z[k] * slk0_3[k]
                  - f_12 * pc_z[k] * slk1_3[k];

        t_76[k] = f_3 * pc_y[k] * smi_58[k];

        t_77[k] = pb_z[k] * slk0_5[k]
                  + f_14 * sli_2[k]
                  - f_12 * pc_z[k] * slk1_5[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_z, pc_y, pc_z, slk0_6, slk0_9, sli_3, \
                         sli_5, slk1_6, slk1_9, smi_59, smi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_z[k] * slk0_6[k]
                  - f_12 * pc_z[k] * slk1_6[k];

        t_79[k] = f_13 * sli_3[k]
                  + f_3 * pc_z[k] * smi_59[k];

        t_80[k] = f_3 * pc_y[k] * smi_61[k];

        t_81[k] = pb_z[k] * slk0_9[k]
                  + f_15 * sli_5[k]
                  - f_12 * pc_z[k] * slk1_9[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_z, pc_y, pc_z, slk0_10, slk0_12, sli_6, \
                         sli_7, slk1_10, slk1_12, smi_62, smi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * slk0_10[k]
                  - f_12 * pc_z[k] * slk1_10[k];

        t_83[k] = f_13 * sli_6[k]
                  + f_3 * pc_z[k] * smi_62[k];

        t_84[k] = pb_z[k] * slk0_12[k]
                  + f_14 * sli_7[k]
                  - f_12 * pc_z[k] * slk1_12[k];

        t_85[k] = f_3 * pc_y[k] * smi_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_z, pc_z, slk0_14, slk0_15, slk0_17, sli_9, \
                         sli_10, sli_11, slk1_14, slk1_15, slk1_17, \
                         smi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * slk0_14[k]
                  + f_16 * sli_9[k]
                  - f_12 * pc_z[k] * slk1_14[k];

        t_87[k] = pb_z[k] * slk0_15[k]
                  - f_12 * pc_z[k] * slk1_15[k];

        t_88[k] = f_13 * sli_10[k]
                  + f_3 * pc_z[k] * smi_66[k];

        t_89[k] = pb_z[k] * slk0_17[k]
                  + f_14 * sli_11[k]
                  - f_12 * pc_z[k] * slk1_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_z, pc_y, pc_z, slk0_18, slk0_20, sli_12, sli_14, \
                         slk1_18, slk1_20, smi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * slk0_18[k]
                  + f_15 * sli_12[k]
                  - f_12 * pc_z[k] * slk1_18[k];

        t_91[k] = f_3 * pc_y[k] * smi_70[k];

        t_92[k] = pb_z[k] * slk0_20[k]
                  + f_17 * sli_14[k]
                  - f_12 * pc_z[k] * slk1_20[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, sli_77, sli_78, sli_79, sli_80, \
                         sli_81, smi_77, smi_78, smi_79, smi_80, \
                         smi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_18 * sli_77[k]
                  + f_3 * pc_x[k] * smi_77[k];

        t_94[k] = f_18 * sli_78[k]
                  + f_3 * pc_x[k] * smi_78[k];

        t_95[k] = f_18 * sli_79[k]
                  + f_3 * pc_x[k] * smi_79[k];

        t_96[k] = f_18 * sli_80[k]
                  + f_3 * pc_x[k] * smi_80[k];

        t_97[k] = f_18 * sli_81[k]
                  + f_3 * pc_x[k] * smi_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_z, slk0_28, sli_21, sli_82, \
                         sli_83, slk1_28, smi_77, smi_82, smi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_18 * sli_82[k]
                  + f_3 * pc_x[k] * smi_82[k];

        t_99[k] = f_18 * sli_83[k]
                  + f_3 * pc_x[k] * smi_83[k];

        t_100[k] = pb_z[k] * slk0_28[k]
                   - f_12 * pc_z[k] * slk1_28[k];

        t_101[k] = f_13 * sli_21[k]
                   + f_3 * pc_z[k] * smi_77[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, smh0_59, smh0_60, smh0_61, smh1_59, \
                         smh1_60, smh1_61, smi_79, smi_80, smi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * smh0_59[k]
                   - f_5 * smh1_59[k]
                   + f_3 * pc_y[k] * smi_79[k];

        t_103[k] = f_6 * smh0_60[k]
                   - f_7 * smh1_60[k]
                   + f_3 * pc_y[k] * smi_80[k];

        t_104[k] = f_8 * smh0_61[k]
                   - f_9 * smh1_61[k]
                   + f_3 * pc_y[k] * smi_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, sli_27, sli_84, \
                         smh0_62, smh0_63, smh1_62, smh1_63, smi_82, smi_83, \
                         smi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * smh0_62[k]
                   - f_11 * smh1_62[k]
                   + f_3 * pc_y[k] * smi_82[k];

        t_106[k] = f_3 * pc_y[k] * smi_83[k];

        t_107[k] = f_13 * sli_27[k]
                   + f_1 * smh0_62[k]
                   - f_2 * smh1_62[k]
                   + f_3 * pc_z[k] * smi_83[k];

        t_108[k] = f_19 * sli_84[k]
                   + f_1 * smh0_63[k]
                   - f_2 * smh1_63[k]
                   + f_3 * pc_x[k] * smi_84[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, pc_z, sli_28, sli_30, sli_87, \
                         smh0_66, smh1_66, smi_84, smi_86, smi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * sli_28[k]
                   + f_3 * pc_y[k] * smi_84[k];

        t_110[k] = f_3 * pc_z[k] * smi_84[k];

        t_111[k] = f_19 * sli_87[k]
                   + f_4 * smh0_66[k]
                   - f_5 * smh1_66[k]
                   + f_3 * pc_x[k] * smi_87[k];

        t_112[k] = f_14 * sli_30[k]
                   + f_3 * pc_y[k] * smi_86[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pc_x, pc_z, sli_89, sli_90, smh0_68, smh0_69, \
                         smh1_68, smh1_69, smi_87, smi_89, smi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_19 * sli_89[k]
                   + f_4 * smh0_68[k]
                   - f_5 * smh1_68[k]
                   + f_3 * pc_x[k] * smi_89[k];

        t_114[k] = f_19 * sli_90[k]
                   + f_6 * smh0_69[k]
                   - f_7 * smh1_69[k]
                   + f_3 * pc_x[k] * smi_90[k];

        t_115[k] = f_3 * pc_z[k] * smi_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pc_x, pc_y, sli_33, sli_93, sli_94, smh0_72, \
                         smh0_73, smh1_72, smh1_73, smi_89, smi_93, \
                         smi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_14 * sli_33[k]
                   + f_3 * pc_y[k] * smi_89[k];

        t_117[k] = f_19 * sli_93[k]
                   + f_6 * smh0_72[k]
                   - f_7 * smh1_72[k]
                   + f_3 * pc_x[k] * smi_93[k];

        t_118[k] = f_19 * sli_94[k]
                   + f_8 * smh0_73[k]
                   - f_9 * smh1_73[k]
                   + f_3 * pc_x[k] * smi_94[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, sli_37, sli_96, smh0_75, \
                         smh1_75, smi_90, smi_93, smi_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * pc_z[k] * smi_90[k];

        t_120[k] = f_19 * sli_96[k]
                   + f_8 * smh0_75[k]
                   - f_9 * smh1_75[k]
                   + f_3 * pc_x[k] * smi_96[k];

        t_121[k] = f_14 * sli_37[k]
                   + f_3 * pc_y[k] * smi_93[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_39 = buffer.data(slk0 + 39);
    const auto *slk0_42 = buffer.data(slk0 + 42);
    const auto *slk0_46 = buffer.data(slk0 + 46);
    const auto *slk0_51 = buffer.data(slk0 + 51);
    const auto *slk0_64 = buffer.data(slk0 + 64);
    const auto *slk0_72 = buffer.data(slk0 + 72);
    const auto *slk0_77 = buffer.data(slk0 + 77);
    const auto *slk0_81 = buffer.data(slk0 + 81);
    const auto *slk0_84 = buffer.data(slk0 + 84);
    const auto *slk0_86 = buffer.data(slk0 + 86);
    const auto *slk0_89 = buffer.data(slk0 + 89);
    const auto *slk0_90 = buffer.data(slk0 + 90);
    const auto *slk0_92 = buffer.data(slk0 + 92);
    const auto *slk0_107 = buffer.data(slk0 + 107);

    const auto *sli_28 = buffer.data(sli + 28);
    const auto *sli_31 = buffer.data(sli + 31);
    const auto *sli_34 = buffer.data(sli + 34);
    const auto *sli_38 = buffer.data(sli + 38);
    const auto *sli_42 = buffer.data(sli + 42);
    const auto *sli_49 = buffer.data(sli + 49);
    const auto *sli_51 = buffer.data(sli + 51);
    const auto *sli_52 = buffer.data(sli + 52);
    const auto *sli_53 = buffer.data(sli + 53);
    const auto *sli_54 = buffer.data(sli + 54);
    const auto *sli_55 = buffer.data(sli + 55);
    const auto *sli_56 = buffer.data(sli + 56);
    const auto *sli_58 = buffer.data(sli + 58);
    const auto *sli_59 = buffer.data(sli + 59);
    const auto *sli_61 = buffer.data(sli + 61);
    const auto *sli_62 = buffer.data(sli + 62);
    const auto *sli_64 = buffer.data(sli + 64);
    const auto *sli_65 = buffer.data(sli + 65);
    const auto *sli_66 = buffer.data(sli + 66);
    const auto *sli_68 = buffer.data(sli + 68);
    const auto *sli_69 = buffer.data(sli + 69);
    const auto *sli_70 = buffer.data(sli + 70);
    const auto *sli_77 = buffer.data(sli + 77);
    const auto *sli_79 = buffer.data(sli + 79);
    const auto *sli_80 = buffer.data(sli + 80);
    const auto *sli_81 = buffer.data(sli + 81);
    const auto *sli_82 = buffer.data(sli + 82);
    const auto *sli_83 = buffer.data(sli + 83);
    const auto *sli_84 = buffer.data(sli + 84);
    const auto *sli_86 = buffer.data(sli + 86);
    const auto *sli_89 = buffer.data(sli + 89);
    const auto *sli_93 = buffer.data(sli + 93);
    const auto *sli_98 = buffer.data(sli + 98);
    const auto *sli_99 = buffer.data(sli + 99);
    const auto *sli_101 = buffer.data(sli + 101);
    const auto *sli_102 = buffer.data(sli + 102);
    const auto *sli_104 = buffer.data(sli + 104);
    const auto *sli_105 = buffer.data(sli + 105);
    const auto *sli_106 = buffer.data(sli + 106);
    const auto *sli_107 = buffer.data(sli + 107);
    const auto *sli_108 = buffer.data(sli + 108);
    const auto *sli_109 = buffer.data(sli + 109);
    const auto *sli_110 = buffer.data(sli + 110);
    const auto *sli_111 = buffer.data(sli + 111);
    const auto *sli_133 = buffer.data(sli + 133);
    const auto *sli_134 = buffer.data(sli + 134);
    const auto *sli_135 = buffer.data(sli + 135);
    const auto *sli_136 = buffer.data(sli + 136);
    const auto *sli_137 = buffer.data(sli + 137);
    const auto *sli_138 = buffer.data(sli + 138);
    const auto *sli_139 = buffer.data(sli + 139);
    const auto *sli_140 = buffer.data(sli + 140);
    const auto *sli_143 = buffer.data(sli + 143);
    const auto *sli_145 = buffer.data(sli + 145);
    const auto *sli_146 = buffer.data(sli + 146);
    const auto *sli_149 = buffer.data(sli + 149);
    const auto *sli_150 = buffer.data(sli + 150);
    const auto *sli_152 = buffer.data(sli + 152);
    const auto *sli_154 = buffer.data(sli + 154);
    const auto *sli_155 = buffer.data(sli + 155);
    const auto *sli_157 = buffer.data(sli + 157);
    const auto *sli_158 = buffer.data(sli + 158);
    const auto *sli_160 = buffer.data(sli + 160);
    const auto *sli_161 = buffer.data(sli + 161);
    const auto *sli_162 = buffer.data(sli + 162);
    const auto *sli_163 = buffer.data(sli + 163);
    const auto *sli_164 = buffer.data(sli + 164);
    const auto *sli_165 = buffer.data(sli + 165);
    const auto *sli_166 = buffer.data(sli + 166);
    const auto *sli_167 = buffer.data(sli + 167);
    const auto *sli_168 = buffer.data(sli + 168);
    const auto *sli_171 = buffer.data(sli + 171);
    const auto *sli_173 = buffer.data(sli + 173);
    const auto *sli_174 = buffer.data(sli + 174);
    const auto *sli_177 = buffer.data(sli + 177);
    const auto *sli_178 = buffer.data(sli + 178);
    const auto *sli_180 = buffer.data(sli + 180);
    const auto *sli_182 = buffer.data(sli + 182);
    const auto *sli_183 = buffer.data(sli + 183);
    const auto *sli_185 = buffer.data(sli + 185);
    const auto *sli_186 = buffer.data(sli + 186);

    const auto *slk1_39 = buffer.data(slk1 + 39);
    const auto *slk1_42 = buffer.data(slk1 + 42);
    const auto *slk1_46 = buffer.data(slk1 + 46);
    const auto *slk1_51 = buffer.data(slk1 + 51);
    const auto *slk1_64 = buffer.data(slk1 + 64);
    const auto *slk1_72 = buffer.data(slk1 + 72);
    const auto *slk1_77 = buffer.data(slk1 + 77);
    const auto *slk1_81 = buffer.data(slk1 + 81);
    const auto *slk1_84 = buffer.data(slk1 + 84);
    const auto *slk1_86 = buffer.data(slk1 + 86);
    const auto *slk1_89 = buffer.data(slk1 + 89);
    const auto *slk1_90 = buffer.data(slk1 + 90);
    const auto *slk1_92 = buffer.data(slk1 + 92);
    const auto *slk1_107 = buffer.data(slk1 + 107);

    const auto *smh0_77 = buffer.data(smh0 + 77);
    const auto *smh0_78 = buffer.data(smh0 + 78);
    const auto *smh0_80 = buffer.data(smh0 + 80);
    const auto *smh0_81 = buffer.data(smh0 + 81);
    const auto *smh0_82 = buffer.data(smh0 + 82);
    const auto *smh0_83 = buffer.data(smh0 + 83);
    const auto *smh0_101 = buffer.data(smh0 + 101);
    const auto *smh0_102 = buffer.data(smh0 + 102);
    const auto *smh0_103 = buffer.data(smh0 + 103);
    const auto *smh0_104 = buffer.data(smh0 + 104);
    const auto *smh0_105 = buffer.data(smh0 + 105);
    const auto *smh0_108 = buffer.data(smh0 + 108);
    const auto *smh0_110 = buffer.data(smh0 + 110);
    const auto *smh0_111 = buffer.data(smh0 + 111);
    const auto *smh0_114 = buffer.data(smh0 + 114);
    const auto *smh0_115 = buffer.data(smh0 + 115);
    const auto *smh0_117 = buffer.data(smh0 + 117);
    const auto *smh0_119 = buffer.data(smh0 + 119);
    const auto *smh0_120 = buffer.data(smh0 + 120);
    const auto *smh0_122 = buffer.data(smh0 + 122);
    const auto *smh0_123 = buffer.data(smh0 + 123);
    const auto *smh0_124 = buffer.data(smh0 + 124);
    const auto *smh0_125 = buffer.data(smh0 + 125);
    const auto *smh0_126 = buffer.data(smh0 + 126);
    const auto *smh0_129 = buffer.data(smh0 + 129);
    const auto *smh0_131 = buffer.data(smh0 + 131);
    const auto *smh0_132 = buffer.data(smh0 + 132);
    const auto *smh0_135 = buffer.data(smh0 + 135);
    const auto *smh0_136 = buffer.data(smh0 + 136);
    const auto *smh0_138 = buffer.data(smh0 + 138);
    const auto *smh0_140 = buffer.data(smh0 + 140);
    const auto *smh0_141 = buffer.data(smh0 + 141);
    const auto *smh0_143 = buffer.data(smh0 + 143);
    const auto *smh0_144 = buffer.data(smh0 + 144);

    const auto *smh1_77 = buffer.data(smh1 + 77);
    const auto *smh1_78 = buffer.data(smh1 + 78);
    const auto *smh1_80 = buffer.data(smh1 + 80);
    const auto *smh1_81 = buffer.data(smh1 + 81);
    const auto *smh1_82 = buffer.data(smh1 + 82);
    const auto *smh1_83 = buffer.data(smh1 + 83);
    const auto *smh1_101 = buffer.data(smh1 + 101);
    const auto *smh1_102 = buffer.data(smh1 + 102);
    const auto *smh1_103 = buffer.data(smh1 + 103);
    const auto *smh1_104 = buffer.data(smh1 + 104);
    const auto *smh1_105 = buffer.data(smh1 + 105);
    const auto *smh1_108 = buffer.data(smh1 + 108);
    const auto *smh1_110 = buffer.data(smh1 + 110);
    const auto *smh1_111 = buffer.data(smh1 + 111);
    const auto *smh1_114 = buffer.data(smh1 + 114);
    const auto *smh1_115 = buffer.data(smh1 + 115);
    const auto *smh1_117 = buffer.data(smh1 + 117);
    const auto *smh1_119 = buffer.data(smh1 + 119);
    const auto *smh1_120 = buffer.data(smh1 + 120);
    const auto *smh1_122 = buffer.data(smh1 + 122);
    const auto *smh1_123 = buffer.data(smh1 + 123);
    const auto *smh1_124 = buffer.data(smh1 + 124);
    const auto *smh1_125 = buffer.data(smh1 + 125);
    const auto *smh1_126 = buffer.data(smh1 + 126);
    const auto *smh1_129 = buffer.data(smh1 + 129);
    const auto *smh1_131 = buffer.data(smh1 + 131);
    const auto *smh1_132 = buffer.data(smh1 + 132);
    const auto *smh1_135 = buffer.data(smh1 + 135);
    const auto *smh1_136 = buffer.data(smh1 + 136);
    const auto *smh1_138 = buffer.data(smh1 + 138);
    const auto *smh1_140 = buffer.data(smh1 + 140);
    const auto *smh1_141 = buffer.data(smh1 + 141);
    const auto *smh1_143 = buffer.data(smh1 + 143);
    const auto *smh1_144 = buffer.data(smh1 + 144);

    const auto *smi_94 = buffer.data(smi + 94);
    const auto *smi_98 = buffer.data(smi + 98);
    const auto *smi_99 = buffer.data(smi + 99);
    const auto *smi_101 = buffer.data(smi + 101);
    const auto *smi_102 = buffer.data(smi + 102);
    const auto *smi_104 = buffer.data(smi + 104);
    const auto *smi_105 = buffer.data(smi + 105);
    const auto *smi_106 = buffer.data(smi + 106);
    const auto *smi_107 = buffer.data(smi + 107);
    const auto *smi_108 = buffer.data(smi + 108);
    const auto *smi_109 = buffer.data(smi + 109);
    const auto *smi_110 = buffer.data(smi + 110);
    const auto *smi_111 = buffer.data(smi + 111);
    const auto *smi_112 = buffer.data(smi + 112);
    const auto *smi_114 = buffer.data(smi + 114);
    const auto *smi_115 = buffer.data(smi + 115);
    const auto *smi_117 = buffer.data(smi + 117);
    const auto *smi_118 = buffer.data(smi + 118);
    const auto *smi_121 = buffer.data(smi + 121);
    const auto *smi_122 = buffer.data(smi + 122);
    const auto *smi_126 = buffer.data(smi + 126);
    const auto *smi_133 = buffer.data(smi + 133);
    const auto *smi_134 = buffer.data(smi + 134);
    const auto *smi_135 = buffer.data(smi + 135);
    const auto *smi_136 = buffer.data(smi + 136);
    const auto *smi_137 = buffer.data(smi + 137);
    const auto *smi_138 = buffer.data(smi + 138);
    const auto *smi_139 = buffer.data(smi + 139);
    const auto *smi_140 = buffer.data(smi + 140);
    const auto *smi_142 = buffer.data(smi + 142);
    const auto *smi_143 = buffer.data(smi + 143);
    const auto *smi_145 = buffer.data(smi + 145);
    const auto *smi_146 = buffer.data(smi + 146);
    const auto *smi_149 = buffer.data(smi + 149);
    const auto *smi_150 = buffer.data(smi + 150);
    const auto *smi_152 = buffer.data(smi + 152);
    const auto *smi_154 = buffer.data(smi + 154);
    const auto *smi_155 = buffer.data(smi + 155);
    const auto *smi_157 = buffer.data(smi + 157);
    const auto *smi_158 = buffer.data(smi + 158);
    const auto *smi_160 = buffer.data(smi + 160);
    const auto *smi_161 = buffer.data(smi + 161);
    const auto *smi_162 = buffer.data(smi + 162);
    const auto *smi_163 = buffer.data(smi + 163);
    const auto *smi_164 = buffer.data(smi + 164);
    const auto *smi_165 = buffer.data(smi + 165);
    const auto *smi_166 = buffer.data(smi + 166);
    const auto *smi_167 = buffer.data(smi + 167);
    const auto *smi_168 = buffer.data(smi + 168);
    const auto *smi_170 = buffer.data(smi + 170);
    const auto *smi_171 = buffer.data(smi + 171);
    const auto *smi_173 = buffer.data(smi + 173);
    const auto *smi_174 = buffer.data(smi + 174);
    const auto *smi_177 = buffer.data(smi + 177);
    const auto *smi_178 = buffer.data(smi + 178);
    const auto *smi_180 = buffer.data(smi + 180);
    const auto *smi_182 = buffer.data(smi + 182);
    const auto *smi_183 = buffer.data(smi + 183);
    const auto *smi_185 = buffer.data(smi + 185);
    const auto *smi_186 = buffer.data(smi + 186);

#pragma omp simd aligned(t_122, t_123, t_124, pc_x, pc_z, sli_98, sli_99, smh0_77, smh0_78, \
                         smh1_77, smh1_78, smi_94, smi_98, smi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_19 * sli_98[k]
                   + f_8 * smh0_77[k]
                   - f_9 * smh1_77[k]
                   + f_3 * pc_x[k] * smi_98[k];

        t_123[k] = f_19 * sli_99[k]
                   + f_10 * smh0_78[k]
                   - f_11 * smh1_78[k]
                   + f_3 * pc_x[k] * smi_99[k];

        t_124[k] = f_3 * pc_z[k] * smi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_y, sli_42, sli_101, sli_102, smh0_80, \
                         smh0_81, smh1_80, smh1_81, smi_98, smi_101, \
                         smi_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_19 * sli_101[k]
                   + f_10 * smh0_80[k]
                   - f_11 * smh1_80[k]
                   + f_3 * pc_x[k] * smi_101[k];

        t_126[k] = f_19 * sli_102[k]
                   + f_10 * smh0_81[k]
                   - f_11 * smh1_81[k]
                   + f_3 * pc_x[k] * smi_102[k];

        t_127[k] = f_14 * sli_42[k]
                   + f_3 * pc_y[k] * smi_98[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, sli_104, sli_105, sli_106, sli_107, \
                         smh0_83, smh1_83, smi_104, smi_105, smi_106, \
                         smi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_19 * sli_104[k]
                   + f_10 * smh0_83[k]
                   - f_11 * smh1_83[k]
                   + f_3 * pc_x[k] * smi_104[k];

        t_129[k] = f_19 * sli_105[k]
                   + f_3 * pc_x[k] * smi_105[k];

        t_130[k] = f_19 * sli_106[k]
                   + f_3 * pc_x[k] * smi_106[k];

        t_131[k] = f_19 * sli_107[k]
                   + f_3 * pc_x[k] * smi_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, sli_108, sli_109, sli_110, sli_111, \
                         smi_108, smi_109, smi_110, smi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_19 * sli_108[k]
                   + f_3 * pc_x[k] * smi_108[k];

        t_133[k] = f_19 * sli_109[k]
                   + f_3 * pc_x[k] * smi_109[k];

        t_134[k] = f_19 * sli_110[k]
                   + f_3 * pc_x[k] * smi_110[k];

        t_135[k] = f_19 * sli_111[k]
                   + f_3 * pc_x[k] * smi_111[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, pc_z, sli_49, sli_51, smh0_78, smh0_80, \
                         smh1_78, smh1_80, smi_105, smi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * sli_49[k]
                   + f_1 * smh0_78[k]
                   - f_2 * smh1_78[k]
                   + f_3 * pc_y[k] * smi_105[k];

        t_137[k] = f_3 * pc_z[k] * smi_105[k];

        t_138[k] = f_14 * sli_51[k]
                   + f_4 * smh0_80[k]
                   - f_5 * smh1_80[k]
                   + f_3 * pc_y[k] * smi_107[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pc_y, sli_52, sli_53, sli_54, smh0_81, smh0_82, \
                         smh0_83, smh1_81, smh1_82, smh1_83, smi_108, smi_109, \
                         smi_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * sli_52[k]
                   + f_6 * smh0_81[k]
                   - f_7 * smh1_81[k]
                   + f_3 * pc_y[k] * smi_108[k];

        t_140[k] = f_14 * sli_53[k]
                   + f_8 * smh0_82[k]
                   - f_9 * smh1_82[k]
                   + f_3 * pc_y[k] * smi_109[k];

        t_141[k] = f_14 * sli_54[k]
                   + f_10 * smh0_83[k]
                   - f_11 * smh1_83[k]
                   + f_3 * pc_y[k] * smi_110[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_y, pc_y, pc_z, slk0_72, sli_55, \
                         sli_56, slk1_72, smh0_83, smh1_83, smi_111, \
                         smi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * sli_55[k]
                   + f_3 * pc_y[k] * smi_111[k];

        t_143[k] = f_1 * smh0_83[k]
                   - f_2 * smh1_83[k]
                   + f_3 * pc_z[k] * smi_111[k];

        t_144[k] = pb_y[k] * slk0_72[k]
                   - f_12 * pc_y[k] * slk1_72[k];

        t_145[k] = f_13 * sli_56[k]
                   + f_3 * pc_y[k] * smi_112[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_y, pb_z, pc_y, pc_z, slk0_39, slk0_77, \
                         sli_28, sli_58, slk1_39, slk1_77, smi_112, \
                         smi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * sli_28[k]
                   + f_3 * pc_z[k] * smi_112[k];

        t_147[k] = pb_z[k] * slk0_39[k]
                   - f_12 * pc_z[k] * slk1_39[k];

        t_148[k] = f_13 * sli_58[k]
                   + f_3 * pc_y[k] * smi_114[k];

        t_149[k] = pb_y[k] * slk0_77[k]
                   - f_12 * pc_y[k] * slk1_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, pb_z, pc_y, pc_z, slk0_42, slk0_81, \
                         sli_31, sli_61, slk1_42, slk1_81, smi_115, \
                         smi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * slk0_42[k]
                   - f_12 * pc_z[k] * slk1_42[k];

        t_151[k] = f_13 * sli_31[k]
                   + f_3 * pc_z[k] * smi_115[k];

        t_152[k] = f_13 * sli_61[k]
                   + f_3 * pc_y[k] * smi_117[k];

        t_153[k] = pb_y[k] * slk0_81[k]
                   - f_12 * pc_y[k] * slk1_81[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pb_z, pc_y, pc_z, slk0_46, slk0_84, \
                         sli_34, sli_64, slk1_46, slk1_84, smi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * slk0_46[k]
                   - f_12 * pc_z[k] * slk1_46[k];

        t_155[k] = f_13 * sli_34[k]
                   + f_3 * pc_z[k] * smi_118[k];

        t_156[k] = pb_y[k] * slk0_84[k]
                   + f_14 * sli_64[k]
                   - f_12 * pc_y[k] * slk1_84[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_y, pb_z, pc_y, pc_z, slk0_51, slk0_86, \
                         sli_38, sli_65, slk1_51, slk1_86, smi_121, \
                         smi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * sli_65[k]
                   + f_3 * pc_y[k] * smi_121[k];

        t_158[k] = pb_y[k] * slk0_86[k]
                   - f_12 * pc_y[k] * slk1_86[k];

        t_159[k] = pb_z[k] * slk0_51[k]
                   - f_12 * pc_z[k] * slk1_51[k];

        t_160[k] = f_13 * sli_38[k]
                   + f_3 * pc_z[k] * smi_122[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pc_y, slk0_89, slk0_90, slk0_92, \
                         sli_68, sli_69, sli_70, slk1_89, slk1_90, slk1_92, \
                         smi_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * slk0_89[k]
                   + f_15 * sli_68[k]
                   - f_12 * pc_y[k] * slk1_89[k];

        t_162[k] = pb_y[k] * slk0_90[k]
                   + f_14 * sli_69[k]
                   - f_12 * pc_y[k] * slk1_90[k];

        t_163[k] = f_13 * sli_70[k]
                   + f_3 * pc_y[k] * smi_126[k];

        t_164[k] = pb_y[k] * slk0_92[k]
                   - f_12 * pc_y[k] * slk1_92[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, sli_133, sli_134, sli_135, \
                         sli_136, sli_137, smi_133, smi_134, smi_135, smi_136, \
                         smi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_19 * sli_133[k]
                   + f_3 * pc_x[k] * smi_133[k];

        t_166[k] = f_19 * sli_134[k]
                   + f_3 * pc_x[k] * smi_134[k];

        t_167[k] = f_19 * sli_135[k]
                   + f_3 * pc_x[k] * smi_135[k];

        t_168[k] = f_19 * sli_136[k]
                   + f_3 * pc_x[k] * smi_136[k];

        t_169[k] = f_19 * sli_137[k]
                   + f_3 * pc_x[k] * smi_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_z, pc_x, pc_z, slk0_64, sli_49, \
                         sli_138, sli_139, slk1_64, smi_133, smi_138, \
                         smi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_19 * sli_138[k]
                   + f_3 * pc_x[k] * smi_138[k];

        t_171[k] = f_19 * sli_139[k]
                   + f_3 * pc_x[k] * smi_139[k];

        t_172[k] = pb_z[k] * slk0_64[k]
                   - f_12 * pc_z[k] * slk1_64[k];

        t_173[k] = f_13 * sli_49[k]
                   + f_3 * pc_z[k] * smi_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, sli_79, sli_80, sli_81, smh0_101, \
                         smh0_102, smh0_103, smh1_101, smh1_102, smh1_103, smi_135, smi_136, \
                         smi_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * sli_79[k]
                   + f_4 * smh0_101[k]
                   - f_5 * smh1_101[k]
                   + f_3 * pc_y[k] * smi_135[k];

        t_175[k] = f_13 * sli_80[k]
                   + f_6 * smh0_102[k]
                   - f_7 * smh1_102[k]
                   + f_3 * pc_y[k] * smi_136[k];

        t_176[k] = f_13 * sli_81[k]
                   + f_8 * smh0_103[k]
                   - f_9 * smh1_103[k]
                   + f_3 * pc_y[k] * smi_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_y, pc_y, slk0_107, sli_82, sli_83, slk1_107, \
                         smh0_104, smh1_104, smi_138, smi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * sli_82[k]
                   + f_10 * smh0_104[k]
                   - f_11 * smh1_104[k]
                   + f_3 * pc_y[k] * smi_138[k];

        t_178[k] = f_13 * sli_83[k]
                   + f_3 * pc_y[k] * smi_139[k];

        t_179[k] = pb_y[k] * slk0_107[k]
                   - f_12 * pc_y[k] * slk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, sli_56, sli_140, \
                         sli_143, smh0_105, smh0_108, smh1_105, smh1_108, smi_140, \
                         smi_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_19 * sli_140[k]
                   + f_1 * smh0_105[k]
                   - f_2 * smh1_105[k]
                   + f_3 * pc_x[k] * smi_140[k];

        t_181[k] = f_3 * pc_y[k] * smi_140[k];

        t_182[k] = f_14 * sli_56[k]
                   + f_3 * pc_z[k] * smi_140[k];

        t_183[k] = f_19 * sli_143[k]
                   + f_4 * smh0_108[k]
                   - f_5 * smh1_108[k]
                   + f_3 * pc_x[k] * smi_143[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, sli_145, sli_146, smh0_110, \
                         smh0_111, smh1_110, smh1_111, smi_142, smi_145, \
                         smi_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_y[k] * smi_142[k];

        t_185[k] = f_19 * sli_145[k]
                   + f_4 * smh0_110[k]
                   - f_5 * smh1_110[k]
                   + f_3 * pc_x[k] * smi_145[k];

        t_186[k] = f_19 * sli_146[k]
                   + f_6 * smh0_111[k]
                   - f_7 * smh1_111[k]
                   + f_3 * pc_x[k] * smi_146[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pc_x, pc_y, pc_z, sli_59, sli_149, smh0_114, \
                         smh1_114, smi_143, smi_145, smi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_14 * sli_59[k]
                   + f_3 * pc_z[k] * smi_143[k];

        t_188[k] = f_3 * pc_y[k] * smi_145[k];

        t_189[k] = f_19 * sli_149[k]
                   + f_6 * smh0_114[k]
                   - f_7 * smh1_114[k]
                   + f_3 * pc_x[k] * smi_149[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_z, sli_62, sli_150, sli_152, smh0_115, \
                         smh0_117, smh1_115, smh1_117, smi_146, smi_150, \
                         smi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_19 * sli_150[k]
                   + f_8 * smh0_115[k]
                   - f_9 * smh1_115[k]
                   + f_3 * pc_x[k] * smi_150[k];

        t_191[k] = f_14 * sli_62[k]
                   + f_3 * pc_z[k] * smi_146[k];

        t_192[k] = f_19 * sli_152[k]
                   + f_8 * smh0_117[k]
                   - f_9 * smh1_117[k]
                   + f_3 * pc_x[k] * smi_152[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, sli_154, sli_155, smh0_119, \
                         smh0_120, smh1_119, smh1_120, smi_149, smi_154, \
                         smi_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * smi_149[k];

        t_194[k] = f_19 * sli_154[k]
                   + f_8 * smh0_119[k]
                   - f_9 * smh1_119[k]
                   + f_3 * pc_x[k] * smi_154[k];

        t_195[k] = f_19 * sli_155[k]
                   + f_10 * smh0_120[k]
                   - f_11 * smh1_120[k]
                   + f_3 * pc_x[k] * smi_155[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pc_x, pc_z, sli_66, sli_157, sli_158, smh0_122, \
                         smh0_123, smh1_122, smh1_123, smi_150, smi_157, \
                         smi_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_14 * sli_66[k]
                   + f_3 * pc_z[k] * smi_150[k];

        t_197[k] = f_19 * sli_157[k]
                   + f_10 * smh0_122[k]
                   - f_11 * smh1_122[k]
                   + f_3 * pc_x[k] * smi_157[k];

        t_198[k] = f_19 * sli_158[k]
                   + f_10 * smh0_123[k]
                   - f_11 * smh1_123[k]
                   + f_3 * pc_x[k] * smi_158[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pc_x, pc_y, sli_160, sli_161, sli_162, \
                         smh0_125, smh1_125, smi_154, smi_160, smi_161, \
                         smi_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_3 * pc_y[k] * smi_154[k];

        t_200[k] = f_19 * sli_160[k]
                   + f_10 * smh0_125[k]
                   - f_11 * smh1_125[k]
                   + f_3 * pc_x[k] * smi_160[k];

        t_201[k] = f_19 * sli_161[k]
                   + f_3 * pc_x[k] * smi_161[k];

        t_202[k] = f_19 * sli_162[k]
                   + f_3 * pc_x[k] * smi_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pc_x, sli_163, sli_164, sli_165, \
                         sli_166, sli_167, smi_163, smi_164, smi_165, smi_166, \
                         smi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_19 * sli_163[k]
                   + f_3 * pc_x[k] * smi_163[k];

        t_204[k] = f_19 * sli_164[k]
                   + f_3 * pc_x[k] * smi_164[k];

        t_205[k] = f_19 * sli_165[k]
                   + f_3 * pc_x[k] * smi_165[k];

        t_206[k] = f_19 * sli_166[k]
                   + f_3 * pc_x[k] * smi_166[k];

        t_207[k] = f_19 * sli_167[k]
                   + f_3 * pc_x[k] * smi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_y, pc_z, sli_77, smh0_120, smh0_122, \
                         smh0_123, smh1_120, smh1_122, smh1_123, smi_161, smi_163, \
                         smi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * smh0_120[k]
                   - f_2 * smh1_120[k]
                   + f_3 * pc_y[k] * smi_161[k];

        t_209[k] = f_14 * sli_77[k]
                   + f_3 * pc_z[k] * smi_161[k];

        t_210[k] = f_4 * smh0_122[k]
                   - f_5 * smh1_122[k]
                   + f_3 * pc_y[k] * smi_163[k];

        t_211[k] = f_6 * smh0_123[k]
                   - f_7 * smh1_123[k]
                   + f_3 * pc_y[k] * smi_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_y, pc_z, sli_83, smh0_124, smh0_125, \
                         smh1_124, smh1_125, smi_165, smi_166, \
                         smi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * smh0_124[k]
                   - f_9 * smh1_124[k]
                   + f_3 * pc_y[k] * smi_165[k];

        t_213[k] = f_10 * smh0_125[k]
                   - f_11 * smh1_125[k]
                   + f_3 * pc_y[k] * smi_166[k];

        t_214[k] = f_3 * pc_y[k] * smi_167[k];

        t_215[k] = f_14 * sli_83[k]
                   + f_1 * smh0_125[k]
                   - f_2 * smh1_125[k]
                   + f_3 * pc_z[k] * smi_167[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, pc_z, sli_84, sli_168, \
                         sli_171, smh0_126, smh0_129, smh1_126, smh1_129, smi_168, \
                         smi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_20 * sli_168[k]
                   + f_1 * smh0_126[k]
                   - f_2 * smh1_126[k]
                   + f_3 * pc_x[k] * smi_168[k];

        t_217[k] = f_15 * sli_84[k]
                   + f_3 * pc_y[k] * smi_168[k];

        t_218[k] = f_3 * pc_z[k] * smi_168[k];

        t_219[k] = f_20 * sli_171[k]
                   + f_4 * smh0_129[k]
                   - f_5 * smh1_129[k]
                   + f_3 * pc_x[k] * smi_171[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_x, pc_y, sli_86, sli_173, sli_174, smh0_131, \
                         smh0_132, smh1_131, smh1_132, smi_170, smi_173, \
                         smi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * sli_86[k]
                   + f_3 * pc_y[k] * smi_170[k];

        t_221[k] = f_20 * sli_173[k]
                   + f_4 * smh0_131[k]
                   - f_5 * smh1_131[k]
                   + f_3 * pc_x[k] * smi_173[k];

        t_222[k] = f_20 * sli_174[k]
                   + f_6 * smh0_132[k]
                   - f_7 * smh1_132[k]
                   + f_3 * pc_x[k] * smi_174[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pc_x, pc_y, pc_z, sli_89, sli_177, smh0_135, \
                         smh1_135, smi_171, smi_173, smi_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_z[k] * smi_171[k];

        t_224[k] = f_15 * sli_89[k]
                   + f_3 * pc_y[k] * smi_173[k];

        t_225[k] = f_20 * sli_177[k]
                   + f_6 * smh0_135[k]
                   - f_7 * smh1_135[k]
                   + f_3 * pc_x[k] * smi_177[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_x, pc_z, sli_178, sli_180, smh0_136, \
                         smh0_138, smh1_136, smh1_138, smi_174, smi_178, \
                         smi_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_20 * sli_178[k]
                   + f_8 * smh0_136[k]
                   - f_9 * smh1_136[k]
                   + f_3 * pc_x[k] * smi_178[k];

        t_227[k] = f_3 * pc_z[k] * smi_174[k];

        t_228[k] = f_20 * sli_180[k]
                   + f_8 * smh0_138[k]
                   - f_9 * smh1_138[k]
                   + f_3 * pc_x[k] * smi_180[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, sli_93, sli_182, sli_183, smh0_140, \
                         smh0_141, smh1_140, smh1_141, smi_177, smi_182, \
                         smi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_15 * sli_93[k]
                   + f_3 * pc_y[k] * smi_177[k];

        t_230[k] = f_20 * sli_182[k]
                   + f_8 * smh0_140[k]
                   - f_9 * smh1_140[k]
                   + f_3 * pc_x[k] * smi_182[k];

        t_231[k] = f_20 * sli_183[k]
                   + f_10 * smh0_141[k]
                   - f_11 * smh1_141[k]
                   + f_3 * pc_x[k] * smi_183[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_z, sli_185, sli_186, smh0_143, \
                         smh0_144, smh1_143, smh1_144, smi_178, smi_185, \
                         smi_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * pc_z[k] * smi_178[k];

        t_233[k] = f_20 * sli_185[k]
                   + f_10 * smh0_143[k]
                   - f_11 * smh1_143[k]
                   + f_3 * pc_x[k] * smi_185[k];

        t_234[k] = f_20 * sli_186[k]
                   + f_10 * smh0_144[k]
                   - f_11 * smh1_144[k]
                   + f_3 * pc_x[k] * smi_186[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_108 = buffer.data(slk0 + 108);
    const auto *slk0_111 = buffer.data(slk0 + 111);
    const auto *slk0_114 = buffer.data(slk0 + 114);
    const auto *slk0_118 = buffer.data(slk0 + 118);
    const auto *slk0_120 = buffer.data(slk0 + 120);
    const auto *slk0_123 = buffer.data(slk0 + 123);
    const auto *slk0_125 = buffer.data(slk0 + 125);
    const auto *slk0_126 = buffer.data(slk0 + 126);
    const auto *slk0_136 = buffer.data(slk0 + 136);
    const auto *slk0_180 = buffer.data(slk0 + 180);
    const auto *slk0_183 = buffer.data(slk0 + 183);
    const auto *slk0_185 = buffer.data(slk0 + 185);
    const auto *slk0_186 = buffer.data(slk0 + 186);
    const auto *slk0_189 = buffer.data(slk0 + 189);
    const auto *slk0_190 = buffer.data(slk0 + 190);
    const auto *slk0_192 = buffer.data(slk0 + 192);
    const auto *slk0_194 = buffer.data(slk0 + 194);
    const auto *slk0_195 = buffer.data(slk0 + 195);
    const auto *slk0_197 = buffer.data(slk0 + 197);
    const auto *slk0_198 = buffer.data(slk0 + 198);
    const auto *slk0_200 = buffer.data(slk0 + 200);
    const auto *slk0_215 = buffer.data(slk0 + 215);

    const auto *sli_84 = buffer.data(sli + 84);
    const auto *sli_87 = buffer.data(sli + 87);
    const auto *sli_90 = buffer.data(sli + 90);
    const auto *sli_91 = buffer.data(sli + 91);
    const auto *sli_94 = buffer.data(sli + 94);
    const auto *sli_95 = buffer.data(sli + 95);
    const auto *sli_96 = buffer.data(sli + 96);
    const auto *sli_98 = buffer.data(sli + 98);
    const auto *sli_105 = buffer.data(sli + 105);
    const auto *sli_107 = buffer.data(sli + 107);
    const auto *sli_108 = buffer.data(sli + 108);
    const auto *sli_109 = buffer.data(sli + 109);
    const auto *sli_110 = buffer.data(sli + 110);
    const auto *sli_111 = buffer.data(sli + 111);
    const auto *sli_112 = buffer.data(sli + 112);
    const auto *sli_114 = buffer.data(sli + 114);
    const auto *sli_115 = buffer.data(sli + 115);
    const auto *sli_117 = buffer.data(sli + 117);
    const auto *sli_118 = buffer.data(sli + 118);
    const auto *sli_121 = buffer.data(sli + 121);
    const auto *sli_122 = buffer.data(sli + 122);
    const auto *sli_126 = buffer.data(sli + 126);
    const auto *sli_133 = buffer.data(sli + 133);
    const auto *sli_135 = buffer.data(sli + 135);
    const auto *sli_136 = buffer.data(sli + 136);
    const auto *sli_137 = buffer.data(sli + 137);
    const auto *sli_138 = buffer.data(sli + 138);
    const auto *sli_139 = buffer.data(sli + 139);
    const auto *sli_140 = buffer.data(sli + 140);
    const auto *sli_141 = buffer.data(sli + 141);
    const auto *sli_142 = buffer.data(sli + 142);
    const auto *sli_143 = buffer.data(sli + 143);
    const auto *sli_145 = buffer.data(sli + 145);
    const auto *sli_146 = buffer.data(sli + 146);
    const auto *sli_148 = buffer.data(sli + 148);
    const auto *sli_149 = buffer.data(sli + 149);
    const auto *sli_150 = buffer.data(sli + 150);
    const auto *sli_152 = buffer.data(sli + 152);
    const auto *sli_153 = buffer.data(sli + 153);
    const auto *sli_154 = buffer.data(sli + 154);
    const auto *sli_161 = buffer.data(sli + 161);
    const auto *sli_163 = buffer.data(sli + 163);
    const auto *sli_164 = buffer.data(sli + 164);
    const auto *sli_165 = buffer.data(sli + 165);
    const auto *sli_166 = buffer.data(sli + 166);
    const auto *sli_167 = buffer.data(sli + 167);
    const auto *sli_188 = buffer.data(sli + 188);
    const auto *sli_189 = buffer.data(sli + 189);
    const auto *sli_190 = buffer.data(sli + 190);
    const auto *sli_191 = buffer.data(sli + 191);
    const auto *sli_192 = buffer.data(sli + 192);
    const auto *sli_193 = buffer.data(sli + 193);
    const auto *sli_194 = buffer.data(sli + 194);
    const auto *sli_195 = buffer.data(sli + 195);
    const auto *sli_201 = buffer.data(sli + 201);
    const auto *sli_205 = buffer.data(sli + 205);
    const auto *sli_210 = buffer.data(sli + 210);
    const auto *sli_216 = buffer.data(sli + 216);
    const auto *sli_217 = buffer.data(sli + 217);
    const auto *sli_218 = buffer.data(sli + 218);
    const auto *sli_219 = buffer.data(sli + 219);
    const auto *sli_220 = buffer.data(sli + 220);
    const auto *sli_221 = buffer.data(sli + 221);
    const auto *sli_222 = buffer.data(sli + 222);
    const auto *sli_223 = buffer.data(sli + 223);
    const auto *sli_245 = buffer.data(sli + 245);
    const auto *sli_246 = buffer.data(sli + 246);
    const auto *sli_247 = buffer.data(sli + 247);
    const auto *sli_248 = buffer.data(sli + 248);
    const auto *sli_249 = buffer.data(sli + 249);
    const auto *sli_250 = buffer.data(sli + 250);
    const auto *sli_251 = buffer.data(sli + 251);
    const auto *sli_252 = buffer.data(sli + 252);
    const auto *sli_255 = buffer.data(sli + 255);
    const auto *sli_257 = buffer.data(sli + 257);
    const auto *sli_258 = buffer.data(sli + 258);
    const auto *sli_261 = buffer.data(sli + 261);
    const auto *sli_262 = buffer.data(sli + 262);
    const auto *sli_264 = buffer.data(sli + 264);
    const auto *sli_266 = buffer.data(sli + 266);
    const auto *sli_267 = buffer.data(sli + 267);
    const auto *sli_269 = buffer.data(sli + 269);
    const auto *sli_270 = buffer.data(sli + 270);
    const auto *sli_272 = buffer.data(sli + 272);
    const auto *sli_273 = buffer.data(sli + 273);
    const auto *sli_274 = buffer.data(sli + 274);

    const auto *slk1_108 = buffer.data(slk1 + 108);
    const auto *slk1_111 = buffer.data(slk1 + 111);
    const auto *slk1_114 = buffer.data(slk1 + 114);
    const auto *slk1_118 = buffer.data(slk1 + 118);
    const auto *slk1_120 = buffer.data(slk1 + 120);
    const auto *slk1_123 = buffer.data(slk1 + 123);
    const auto *slk1_125 = buffer.data(slk1 + 125);
    const auto *slk1_126 = buffer.data(slk1 + 126);
    const auto *slk1_136 = buffer.data(slk1 + 136);
    const auto *slk1_180 = buffer.data(slk1 + 180);
    const auto *slk1_183 = buffer.data(slk1 + 183);
    const auto *slk1_185 = buffer.data(slk1 + 185);
    const auto *slk1_186 = buffer.data(slk1 + 186);
    const auto *slk1_189 = buffer.data(slk1 + 189);
    const auto *slk1_190 = buffer.data(slk1 + 190);
    const auto *slk1_192 = buffer.data(slk1 + 192);
    const auto *slk1_194 = buffer.data(slk1 + 194);
    const auto *slk1_195 = buffer.data(slk1 + 195);
    const auto *slk1_197 = buffer.data(slk1 + 197);
    const auto *slk1_198 = buffer.data(slk1 + 198);
    const auto *slk1_200 = buffer.data(slk1 + 200);
    const auto *slk1_215 = buffer.data(slk1 + 215);

    const auto *smh0_141 = buffer.data(smh0 + 141);
    const auto *smh0_143 = buffer.data(smh0 + 143);
    const auto *smh0_144 = buffer.data(smh0 + 144);
    const auto *smh0_145 = buffer.data(smh0 + 145);
    const auto *smh0_146 = buffer.data(smh0 + 146);
    const auto *smh0_152 = buffer.data(smh0 + 152);
    const auto *smh0_156 = buffer.data(smh0 + 156);
    const auto *smh0_161 = buffer.data(smh0 + 161);
    const auto *smh0_164 = buffer.data(smh0 + 164);
    const auto *smh0_165 = buffer.data(smh0 + 165);
    const auto *smh0_166 = buffer.data(smh0 + 166);
    const auto *smh0_167 = buffer.data(smh0 + 167);
    const auto *smh0_183 = buffer.data(smh0 + 183);
    const auto *smh0_185 = buffer.data(smh0 + 185);
    const auto *smh0_186 = buffer.data(smh0 + 186);
    const auto *smh0_187 = buffer.data(smh0 + 187);
    const auto *smh0_188 = buffer.data(smh0 + 188);
    const auto *smh0_189 = buffer.data(smh0 + 189);
    const auto *smh0_192 = buffer.data(smh0 + 192);
    const auto *smh0_194 = buffer.data(smh0 + 194);
    const auto *smh0_195 = buffer.data(smh0 + 195);
    const auto *smh0_198 = buffer.data(smh0 + 198);
    const auto *smh0_199 = buffer.data(smh0 + 199);
    const auto *smh0_201 = buffer.data(smh0 + 201);
    const auto *smh0_203 = buffer.data(smh0 + 203);
    const auto *smh0_204 = buffer.data(smh0 + 204);
    const auto *smh0_206 = buffer.data(smh0 + 206);
    const auto *smh0_207 = buffer.data(smh0 + 207);
    const auto *smh0_209 = buffer.data(smh0 + 209);

    const auto *smh1_141 = buffer.data(smh1 + 141);
    const auto *smh1_143 = buffer.data(smh1 + 143);
    const auto *smh1_144 = buffer.data(smh1 + 144);
    const auto *smh1_145 = buffer.data(smh1 + 145);
    const auto *smh1_146 = buffer.data(smh1 + 146);
    const auto *smh1_152 = buffer.data(smh1 + 152);
    const auto *smh1_156 = buffer.data(smh1 + 156);
    const auto *smh1_161 = buffer.data(smh1 + 161);
    const auto *smh1_164 = buffer.data(smh1 + 164);
    const auto *smh1_165 = buffer.data(smh1 + 165);
    const auto *smh1_166 = buffer.data(smh1 + 166);
    const auto *smh1_167 = buffer.data(smh1 + 167);
    const auto *smh1_183 = buffer.data(smh1 + 183);
    const auto *smh1_185 = buffer.data(smh1 + 185);
    const auto *smh1_186 = buffer.data(smh1 + 186);
    const auto *smh1_187 = buffer.data(smh1 + 187);
    const auto *smh1_188 = buffer.data(smh1 + 188);
    const auto *smh1_189 = buffer.data(smh1 + 189);
    const auto *smh1_192 = buffer.data(smh1 + 192);
    const auto *smh1_194 = buffer.data(smh1 + 194);
    const auto *smh1_195 = buffer.data(smh1 + 195);
    const auto *smh1_198 = buffer.data(smh1 + 198);
    const auto *smh1_199 = buffer.data(smh1 + 199);
    const auto *smh1_201 = buffer.data(smh1 + 201);
    const auto *smh1_203 = buffer.data(smh1 + 203);
    const auto *smh1_204 = buffer.data(smh1 + 204);
    const auto *smh1_206 = buffer.data(smh1 + 206);
    const auto *smh1_207 = buffer.data(smh1 + 207);
    const auto *smh1_209 = buffer.data(smh1 + 209);

    const auto *smi_182 = buffer.data(smi + 182);
    const auto *smi_188 = buffer.data(smi + 188);
    const auto *smi_189 = buffer.data(smi + 189);
    const auto *smi_190 = buffer.data(smi + 190);
    const auto *smi_191 = buffer.data(smi + 191);
    const auto *smi_192 = buffer.data(smi + 192);
    const auto *smi_193 = buffer.data(smi + 193);
    const auto *smi_194 = buffer.data(smi + 194);
    const auto *smi_195 = buffer.data(smi + 195);
    const auto *smi_196 = buffer.data(smi + 196);
    const auto *smi_198 = buffer.data(smi + 198);
    const auto *smi_199 = buffer.data(smi + 199);
    const auto *smi_201 = buffer.data(smi + 201);
    const auto *smi_202 = buffer.data(smi + 202);
    const auto *smi_205 = buffer.data(smi + 205);
    const auto *smi_206 = buffer.data(smi + 206);
    const auto *smi_210 = buffer.data(smi + 210);
    const auto *smi_216 = buffer.data(smi + 216);
    const auto *smi_217 = buffer.data(smi + 217);
    const auto *smi_218 = buffer.data(smi + 218);
    const auto *smi_219 = buffer.data(smi + 219);
    const auto *smi_220 = buffer.data(smi + 220);
    const auto *smi_221 = buffer.data(smi + 221);
    const auto *smi_222 = buffer.data(smi + 222);
    const auto *smi_223 = buffer.data(smi + 223);
    const auto *smi_224 = buffer.data(smi + 224);
    const auto *smi_226 = buffer.data(smi + 226);
    const auto *smi_227 = buffer.data(smi + 227);
    const auto *smi_229 = buffer.data(smi + 229);
    const auto *smi_230 = buffer.data(smi + 230);
    const auto *smi_233 = buffer.data(smi + 233);
    const auto *smi_234 = buffer.data(smi + 234);
    const auto *smi_238 = buffer.data(smi + 238);
    const auto *smi_245 = buffer.data(smi + 245);
    const auto *smi_246 = buffer.data(smi + 246);
    const auto *smi_247 = buffer.data(smi + 247);
    const auto *smi_248 = buffer.data(smi + 248);
    const auto *smi_249 = buffer.data(smi + 249);
    const auto *smi_250 = buffer.data(smi + 250);
    const auto *smi_251 = buffer.data(smi + 251);
    const auto *smi_252 = buffer.data(smi + 252);
    const auto *smi_254 = buffer.data(smi + 254);
    const auto *smi_255 = buffer.data(smi + 255);
    const auto *smi_257 = buffer.data(smi + 257);
    const auto *smi_258 = buffer.data(smi + 258);
    const auto *smi_261 = buffer.data(smi + 261);
    const auto *smi_262 = buffer.data(smi + 262);
    const auto *smi_264 = buffer.data(smi + 264);
    const auto *smi_266 = buffer.data(smi + 266);
    const auto *smi_267 = buffer.data(smi + 267);
    const auto *smi_269 = buffer.data(smi + 269);
    const auto *smi_270 = buffer.data(smi + 270);
    const auto *smi_272 = buffer.data(smi + 272);
    const auto *smi_273 = buffer.data(smi + 273);
    const auto *smi_274 = buffer.data(smi + 274);

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, sli_98, sli_188, sli_189, \
                         sli_190, smh0_146, smh1_146, smi_182, smi_188, smi_189, \
                         smi_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_15 * sli_98[k]
                   + f_3 * pc_y[k] * smi_182[k];

        t_236[k] = f_20 * sli_188[k]
                   + f_10 * smh0_146[k]
                   - f_11 * smh1_146[k]
                   + f_3 * pc_x[k] * smi_188[k];

        t_237[k] = f_20 * sli_189[k]
                   + f_3 * pc_x[k] * smi_189[k];

        t_238[k] = f_20 * sli_190[k]
                   + f_3 * pc_x[k] * smi_190[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pc_x, sli_191, sli_192, sli_193, \
                         sli_194, sli_195, smi_191, smi_192, smi_193, smi_194, \
                         smi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_20 * sli_191[k]
                   + f_3 * pc_x[k] * smi_191[k];

        t_240[k] = f_20 * sli_192[k]
                   + f_3 * pc_x[k] * smi_192[k];

        t_241[k] = f_20 * sli_193[k]
                   + f_3 * pc_x[k] * smi_193[k];

        t_242[k] = f_20 * sli_194[k]
                   + f_3 * pc_x[k] * smi_194[k];

        t_243[k] = f_20 * sli_195[k]
                   + f_3 * pc_x[k] * smi_195[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, pc_y, pc_z, sli_105, sli_107, smh0_141, \
                         smh0_143, smh1_141, smh1_143, smi_189, \
                         smi_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_15 * sli_105[k]
                   + f_1 * smh0_141[k]
                   - f_2 * smh1_141[k]
                   + f_3 * pc_y[k] * smi_189[k];

        t_245[k] = f_3 * pc_z[k] * smi_189[k];

        t_246[k] = f_15 * sli_107[k]
                   + f_4 * smh0_143[k]
                   - f_5 * smh1_143[k]
                   + f_3 * pc_y[k] * smi_191[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, sli_108, sli_109, sli_110, smh0_144, \
                         smh0_145, smh0_146, smh1_144, smh1_145, smh1_146, smi_192, smi_193, \
                         smi_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_15 * sli_108[k]
                   + f_6 * smh0_144[k]
                   - f_7 * smh1_144[k]
                   + f_3 * pc_y[k] * smi_192[k];

        t_248[k] = f_15 * sli_109[k]
                   + f_8 * smh0_145[k]
                   - f_9 * smh1_145[k]
                   + f_3 * pc_y[k] * smi_193[k];

        t_249[k] = f_15 * sli_110[k]
                   + f_10 * smh0_146[k]
                   - f_11 * smh1_146[k]
                   + f_3 * pc_y[k] * smi_194[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pb_z, pc_y, pc_z, slk0_108, sli_111, \
                         sli_112, slk1_108, smh0_146, smh1_146, smi_195, \
                         smi_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * sli_111[k]
                   + f_3 * pc_y[k] * smi_195[k];

        t_251[k] = f_1 * smh0_146[k]
                   - f_2 * smh1_146[k]
                   + f_3 * pc_z[k] * smi_195[k];

        t_252[k] = pb_z[k] * slk0_108[k]
                   - f_12 * pc_z[k] * slk1_108[k];

        t_253[k] = f_14 * sli_112[k]
                   + f_3 * pc_y[k] * smi_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_z, pc_y, pc_z, slk0_111, sli_84, sli_114, \
                         slk1_111, smi_196, smi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * sli_84[k]
                   + f_3 * pc_z[k] * smi_196[k];

        t_255[k] = pb_z[k] * slk0_111[k]
                   - f_12 * pc_z[k] * slk1_111[k];

        t_256[k] = f_14 * sli_114[k]
                   + f_3 * pc_y[k] * smi_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_z, pc_x, pc_z, slk0_114, sli_87, sli_201, \
                         slk1_114, smh0_152, smh1_152, smi_199, \
                         smi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_20 * sli_201[k]
                   + f_4 * smh0_152[k]
                   - f_5 * smh1_152[k]
                   + f_3 * pc_x[k] * smi_201[k];

        t_258[k] = pb_z[k] * slk0_114[k]
                   - f_12 * pc_z[k] * slk1_114[k];

        t_259[k] = f_13 * sli_87[k]
                   + f_3 * pc_z[k] * smi_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_z, pc_x, pc_y, pc_z, slk0_118, sli_117, \
                         sli_205, slk1_118, smh0_156, smh1_156, smi_201, \
                         smi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * sli_117[k]
                   + f_3 * pc_y[k] * smi_201[k];

        t_261[k] = f_20 * sli_205[k]
                   + f_6 * smh0_156[k]
                   - f_7 * smh1_156[k]
                   + f_3 * pc_x[k] * smi_205[k];

        t_262[k] = pb_z[k] * slk0_118[k]
                   - f_12 * pc_z[k] * slk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_z, pc_y, pc_z, slk0_120, sli_90, sli_91, \
                         sli_121, slk1_120, smi_202, smi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * sli_90[k]
                   + f_3 * pc_z[k] * smi_202[k];

        t_264[k] = pb_z[k] * slk0_120[k]
                   + f_14 * sli_91[k]
                   - f_12 * pc_z[k] * slk1_120[k];

        t_265[k] = f_14 * sli_121[k]
                   + f_3 * pc_y[k] * smi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_z, pc_x, pc_z, slk0_123, sli_94, sli_210, \
                         slk1_123, smh0_161, smh1_161, smi_206, \
                         smi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_20 * sli_210[k]
                   + f_8 * smh0_161[k]
                   - f_9 * smh1_161[k]
                   + f_3 * pc_x[k] * smi_210[k];

        t_267[k] = pb_z[k] * slk0_123[k]
                   - f_12 * pc_z[k] * slk1_123[k];

        t_268[k] = f_13 * sli_94[k]
                   + f_3 * pc_z[k] * smi_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_z, pc_y, pc_z, slk0_125, slk0_126, sli_95, \
                         sli_96, sli_126, slk1_125, slk1_126, smi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pb_z[k] * slk0_125[k]
                   + f_14 * sli_95[k]
                   - f_12 * pc_z[k] * slk1_125[k];

        t_270[k] = pb_z[k] * slk0_126[k]
                   + f_15 * sli_96[k]
                   - f_12 * pc_z[k] * slk1_126[k];

        t_271[k] = f_14 * sli_126[k]
                   + f_3 * pc_y[k] * smi_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, sli_216, sli_217, sli_218, sli_219, \
                         smh0_167, smh1_167, smi_216, smi_217, smi_218, \
                         smi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_20 * sli_216[k]
                   + f_10 * smh0_167[k]
                   - f_11 * smh1_167[k]
                   + f_3 * pc_x[k] * smi_216[k];

        t_273[k] = f_20 * sli_217[k]
                   + f_3 * pc_x[k] * smi_217[k];

        t_274[k] = f_20 * sli_218[k]
                   + f_3 * pc_x[k] * smi_218[k];

        t_275[k] = f_20 * sli_219[k]
                   + f_3 * pc_x[k] * smi_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, sli_220, sli_221, sli_222, sli_223, \
                         smi_220, smi_221, smi_222, smi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_20 * sli_220[k]
                   + f_3 * pc_x[k] * smi_220[k];

        t_277[k] = f_20 * sli_221[k]
                   + f_3 * pc_x[k] * smi_221[k];

        t_278[k] = f_20 * sli_222[k]
                   + f_3 * pc_x[k] * smi_222[k];

        t_279[k] = f_20 * sli_223[k]
                   + f_3 * pc_x[k] * smi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_z, pc_y, pc_z, slk0_136, sli_105, sli_135, \
                         slk1_136, smh0_164, smh1_164, smi_217, \
                         smi_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_z[k] * slk0_136[k]
                   - f_12 * pc_z[k] * slk1_136[k];

        t_281[k] = f_13 * sli_105[k]
                   + f_3 * pc_z[k] * smi_217[k];

        t_282[k] = f_14 * sli_135[k]
                   + f_4 * smh0_164[k]
                   - f_5 * smh1_164[k]
                   + f_3 * pc_y[k] * smi_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, sli_136, sli_137, sli_138, smh0_165, \
                         smh0_166, smh0_167, smh1_165, smh1_166, smh1_167, smi_220, smi_221, \
                         smi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * sli_136[k]
                   + f_6 * smh0_165[k]
                   - f_7 * smh1_165[k]
                   + f_3 * pc_y[k] * smi_220[k];

        t_284[k] = f_14 * sli_137[k]
                   + f_8 * smh0_166[k]
                   - f_9 * smh1_166[k]
                   + f_3 * pc_y[k] * smi_221[k];

        t_285[k] = f_14 * sli_138[k]
                   + f_10 * smh0_167[k]
                   - f_11 * smh1_167[k]
                   + f_3 * pc_y[k] * smi_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, slk0_180, sli_111, \
                         sli_139, sli_140, slk1_180, smh0_167, smh1_167, smi_223, \
                         smi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * sli_139[k]
                   + f_3 * pc_y[k] * smi_223[k];

        t_287[k] = f_13 * sli_111[k]
                   + f_1 * smh0_167[k]
                   - f_2 * smh1_167[k]
                   + f_3 * pc_z[k] * smi_223[k];

        t_288[k] = pb_y[k] * slk0_180[k]
                   - f_12 * pc_y[k] * slk1_180[k];

        t_289[k] = f_13 * sli_140[k]
                   + f_3 * pc_y[k] * smi_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_y, pc_z, slk0_183, slk0_185, \
                         sli_112, sli_141, sli_142, slk1_183, slk1_185, smi_224, \
                         smi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * sli_112[k]
                   + f_3 * pc_z[k] * smi_224[k];

        t_291[k] = pb_y[k] * slk0_183[k]
                   + f_14 * sli_141[k]
                   - f_12 * pc_y[k] * slk1_183[k];

        t_292[k] = f_13 * sli_142[k]
                   + f_3 * pc_y[k] * smi_226[k];

        t_293[k] = pb_y[k] * slk0_185[k]
                   - f_12 * pc_y[k] * slk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pc_y, pc_z, slk0_186, slk0_189, \
                         sli_115, sli_143, sli_145, slk1_186, slk1_189, smi_227, \
                         smi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pb_y[k] * slk0_186[k]
                   + f_15 * sli_143[k]
                   - f_12 * pc_y[k] * slk1_186[k];

        t_295[k] = f_14 * sli_115[k]
                   + f_3 * pc_z[k] * smi_227[k];

        t_296[k] = f_13 * sli_145[k]
                   + f_3 * pc_y[k] * smi_229[k];

        t_297[k] = pb_y[k] * slk0_189[k]
                   - f_12 * pc_y[k] * slk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pb_y, pc_y, pc_z, slk0_190, slk0_192, sli_118, \
                         sli_146, sli_148, slk1_190, slk1_192, \
                         smi_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pb_y[k] * slk0_190[k]
                   + f_16 * sli_146[k]
                   - f_12 * pc_y[k] * slk1_190[k];

        t_299[k] = f_14 * sli_118[k]
                   + f_3 * pc_z[k] * smi_230[k];

        t_300[k] = pb_y[k] * slk0_192[k]
                   + f_14 * sli_148[k]
                   - f_12 * pc_y[k] * slk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_y, pc_y, pc_z, slk0_194, slk0_195, \
                         sli_122, sli_149, sli_150, slk1_194, slk1_195, smi_233, \
                         smi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * sli_149[k]
                   + f_3 * pc_y[k] * smi_233[k];

        t_302[k] = pb_y[k] * slk0_194[k]
                   - f_12 * pc_y[k] * slk1_194[k];

        t_303[k] = pb_y[k] * slk0_195[k]
                   + f_17 * sli_150[k]
                   - f_12 * pc_y[k] * slk1_195[k];

        t_304[k] = f_14 * sli_122[k]
                   + f_3 * pc_z[k] * smi_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_y, pc_y, slk0_197, slk0_198, slk0_200, \
                         sli_152, sli_153, sli_154, slk1_197, slk1_198, slk1_200, \
                         smi_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_y[k] * slk0_197[k]
                   + f_15 * sli_152[k]
                   - f_12 * pc_y[k] * slk1_197[k];

        t_306[k] = pb_y[k] * slk0_198[k]
                   + f_14 * sli_153[k]
                   - f_12 * pc_y[k] * slk1_198[k];

        t_307[k] = f_13 * sli_154[k]
                   + f_3 * pc_y[k] * smi_238[k];

        t_308[k] = pb_y[k] * slk0_200[k]
                   - f_12 * pc_y[k] * slk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, sli_245, sli_246, sli_247, \
                         sli_248, sli_249, smi_245, smi_246, smi_247, smi_248, \
                         smi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_20 * sli_245[k]
                   + f_3 * pc_x[k] * smi_245[k];

        t_310[k] = f_20 * sli_246[k]
                   + f_3 * pc_x[k] * smi_246[k];

        t_311[k] = f_20 * sli_247[k]
                   + f_3 * pc_x[k] * smi_247[k];

        t_312[k] = f_20 * sli_248[k]
                   + f_3 * pc_x[k] * smi_248[k];

        t_313[k] = f_20 * sli_249[k]
                   + f_3 * pc_x[k] * smi_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, sli_133, sli_161, \
                         sli_250, sli_251, smh0_183, smh1_183, smi_245, smi_250, \
                         smi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_20 * sli_250[k]
                   + f_3 * pc_x[k] * smi_250[k];

        t_315[k] = f_20 * sli_251[k]
                   + f_3 * pc_x[k] * smi_251[k];

        t_316[k] = f_13 * sli_161[k]
                   + f_1 * smh0_183[k]
                   - f_2 * smh1_183[k]
                   + f_3 * pc_y[k] * smi_245[k];

        t_317[k] = f_14 * sli_133[k]
                   + f_3 * pc_z[k] * smi_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, sli_163, sli_164, sli_165, smh0_185, \
                         smh0_186, smh0_187, smh1_185, smh1_186, smh1_187, smi_247, smi_248, \
                         smi_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * sli_163[k]
                   + f_4 * smh0_185[k]
                   - f_5 * smh1_185[k]
                   + f_3 * pc_y[k] * smi_247[k];

        t_319[k] = f_13 * sli_164[k]
                   + f_6 * smh0_186[k]
                   - f_7 * smh1_186[k]
                   + f_3 * pc_y[k] * smi_248[k];

        t_320[k] = f_13 * sli_165[k]
                   + f_8 * smh0_187[k]
                   - f_9 * smh1_187[k]
                   + f_3 * pc_y[k] * smi_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_y, pc_y, slk0_215, sli_166, sli_167, \
                         slk1_215, smh0_188, smh1_188, smi_250, \
                         smi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * sli_166[k]
                   + f_10 * smh0_188[k]
                   - f_11 * smh1_188[k]
                   + f_3 * pc_y[k] * smi_250[k];

        t_322[k] = f_13 * sli_167[k]
                   + f_3 * pc_y[k] * smi_251[k];

        t_323[k] = pb_y[k] * slk0_215[k]
                   - f_12 * pc_y[k] * slk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pc_x, pc_y, pc_z, sli_140, sli_252, \
                         sli_255, smh0_189, smh0_192, smh1_189, smh1_192, smi_252, \
                         smi_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_20 * sli_252[k]
                   + f_1 * smh0_189[k]
                   - f_2 * smh1_189[k]
                   + f_3 * pc_x[k] * smi_252[k];

        t_325[k] = f_3 * pc_y[k] * smi_252[k];

        t_326[k] = f_15 * sli_140[k]
                   + f_3 * pc_z[k] * smi_252[k];

        t_327[k] = f_20 * sli_255[k]
                   + f_4 * smh0_192[k]
                   - f_5 * smh1_192[k]
                   + f_3 * pc_x[k] * smi_255[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pc_x, pc_y, sli_257, sli_258, smh0_194, \
                         smh0_195, smh1_194, smh1_195, smi_254, smi_257, \
                         smi_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_3 * pc_y[k] * smi_254[k];

        t_329[k] = f_20 * sli_257[k]
                   + f_4 * smh0_194[k]
                   - f_5 * smh1_194[k]
                   + f_3 * pc_x[k] * smi_257[k];

        t_330[k] = f_20 * sli_258[k]
                   + f_6 * smh0_195[k]
                   - f_7 * smh1_195[k]
                   + f_3 * pc_x[k] * smi_258[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_x, pc_y, pc_z, sli_143, sli_261, smh0_198, \
                         smh1_198, smi_255, smi_257, smi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_15 * sli_143[k]
                   + f_3 * pc_z[k] * smi_255[k];

        t_332[k] = f_3 * pc_y[k] * smi_257[k];

        t_333[k] = f_20 * sli_261[k]
                   + f_6 * smh0_198[k]
                   - f_7 * smh1_198[k]
                   + f_3 * pc_x[k] * smi_261[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pc_x, pc_z, sli_146, sli_262, sli_264, smh0_199, \
                         smh0_201, smh1_199, smh1_201, smi_258, smi_262, \
                         smi_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_20 * sli_262[k]
                   + f_8 * smh0_199[k]
                   - f_9 * smh1_199[k]
                   + f_3 * pc_x[k] * smi_262[k];

        t_335[k] = f_15 * sli_146[k]
                   + f_3 * pc_z[k] * smi_258[k];

        t_336[k] = f_20 * sli_264[k]
                   + f_8 * smh0_201[k]
                   - f_9 * smh1_201[k]
                   + f_3 * pc_x[k] * smi_264[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pc_x, pc_y, sli_266, sli_267, smh0_203, \
                         smh0_204, smh1_203, smh1_204, smi_261, smi_266, \
                         smi_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_3 * pc_y[k] * smi_261[k];

        t_338[k] = f_20 * sli_266[k]
                   + f_8 * smh0_203[k]
                   - f_9 * smh1_203[k]
                   + f_3 * pc_x[k] * smi_266[k];

        t_339[k] = f_20 * sli_267[k]
                   + f_10 * smh0_204[k]
                   - f_11 * smh1_204[k]
                   + f_3 * pc_x[k] * smi_267[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pc_x, pc_z, sli_150, sli_269, sli_270, smh0_206, \
                         smh0_207, smh1_206, smh1_207, smi_262, smi_269, \
                         smi_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_15 * sli_150[k]
                   + f_3 * pc_z[k] * smi_262[k];

        t_341[k] = f_20 * sli_269[k]
                   + f_10 * smh0_206[k]
                   - f_11 * smh1_206[k]
                   + f_3 * pc_x[k] * smi_269[k];

        t_342[k] = f_20 * sli_270[k]
                   + f_10 * smh0_207[k]
                   - f_11 * smh1_207[k]
                   + f_3 * pc_x[k] * smi_270[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pc_x, pc_y, sli_272, sli_273, sli_274, \
                         smh0_209, smh1_209, smi_266, smi_272, smi_273, \
                         smi_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_3 * pc_y[k] * smi_266[k];

        t_344[k] = f_20 * sli_272[k]
                   + f_10 * smh0_209[k]
                   - f_11 * smh1_209[k]
                   + f_3 * pc_x[k] * smi_272[k];

        t_345[k] = f_20 * sli_273[k]
                   + f_3 * pc_x[k] * smi_273[k];

        t_346[k] = f_20 * sli_274[k]
                   + f_3 * pc_x[k] * smi_274[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_216 = buffer.data(slk0 + 216);
    const auto *slk0_219 = buffer.data(slk0 + 219);
    const auto *slk0_222 = buffer.data(slk0 + 222);
    const auto *slk0_226 = buffer.data(slk0 + 226);
    const auto *slk0_228 = buffer.data(slk0 + 228);
    const auto *slk0_231 = buffer.data(slk0 + 231);
    const auto *slk0_233 = buffer.data(slk0 + 233);
    const auto *slk0_234 = buffer.data(slk0 + 234);
    const auto *slk0_244 = buffer.data(slk0 + 244);

    const auto *sli_161 = buffer.data(sli + 161);
    const auto *sli_167 = buffer.data(sli + 167);
    const auto *sli_168 = buffer.data(sli + 168);
    const auto *sli_170 = buffer.data(sli + 170);
    const auto *sli_171 = buffer.data(sli + 171);
    const auto *sli_173 = buffer.data(sli + 173);
    const auto *sli_174 = buffer.data(sli + 174);
    const auto *sli_175 = buffer.data(sli + 175);
    const auto *sli_177 = buffer.data(sli + 177);
    const auto *sli_178 = buffer.data(sli + 178);
    const auto *sli_179 = buffer.data(sli + 179);
    const auto *sli_180 = buffer.data(sli + 180);
    const auto *sli_182 = buffer.data(sli + 182);
    const auto *sli_189 = buffer.data(sli + 189);
    const auto *sli_191 = buffer.data(sli + 191);
    const auto *sli_192 = buffer.data(sli + 192);
    const auto *sli_193 = buffer.data(sli + 193);
    const auto *sli_194 = buffer.data(sli + 194);
    const auto *sli_195 = buffer.data(sli + 195);
    const auto *sli_196 = buffer.data(sli + 196);
    const auto *sli_198 = buffer.data(sli + 198);
    const auto *sli_199 = buffer.data(sli + 199);
    const auto *sli_201 = buffer.data(sli + 201);
    const auto *sli_202 = buffer.data(sli + 202);
    const auto *sli_205 = buffer.data(sli + 205);
    const auto *sli_206 = buffer.data(sli + 206);
    const auto *sli_210 = buffer.data(sli + 210);
    const auto *sli_219 = buffer.data(sli + 219);
    const auto *sli_220 = buffer.data(sli + 220);
    const auto *sli_221 = buffer.data(sli + 221);
    const auto *sli_222 = buffer.data(sli + 222);
    const auto *sli_223 = buffer.data(sli + 223);
    const auto *sli_224 = buffer.data(sli + 224);
    const auto *sli_226 = buffer.data(sli + 226);
    const auto *sli_229 = buffer.data(sli + 229);
    const auto *sli_233 = buffer.data(sli + 233);
    const auto *sli_238 = buffer.data(sli + 238);
    const auto *sli_275 = buffer.data(sli + 275);
    const auto *sli_276 = buffer.data(sli + 276);
    const auto *sli_277 = buffer.data(sli + 277);
    const auto *sli_278 = buffer.data(sli + 278);
    const auto *sli_279 = buffer.data(sli + 279);
    const auto *sli_280 = buffer.data(sli + 280);
    const auto *sli_283 = buffer.data(sli + 283);
    const auto *sli_285 = buffer.data(sli + 285);
    const auto *sli_286 = buffer.data(sli + 286);
    const auto *sli_289 = buffer.data(sli + 289);
    const auto *sli_290 = buffer.data(sli + 290);
    const auto *sli_292 = buffer.data(sli + 292);
    const auto *sli_294 = buffer.data(sli + 294);
    const auto *sli_295 = buffer.data(sli + 295);
    const auto *sli_297 = buffer.data(sli + 297);
    const auto *sli_298 = buffer.data(sli + 298);
    const auto *sli_300 = buffer.data(sli + 300);
    const auto *sli_301 = buffer.data(sli + 301);
    const auto *sli_302 = buffer.data(sli + 302);
    const auto *sli_303 = buffer.data(sli + 303);
    const auto *sli_304 = buffer.data(sli + 304);
    const auto *sli_305 = buffer.data(sli + 305);
    const auto *sli_306 = buffer.data(sli + 306);
    const auto *sli_307 = buffer.data(sli + 307);
    const auto *sli_313 = buffer.data(sli + 313);
    const auto *sli_317 = buffer.data(sli + 317);
    const auto *sli_322 = buffer.data(sli + 322);
    const auto *sli_328 = buffer.data(sli + 328);
    const auto *sli_329 = buffer.data(sli + 329);
    const auto *sli_330 = buffer.data(sli + 330);
    const auto *sli_331 = buffer.data(sli + 331);
    const auto *sli_332 = buffer.data(sli + 332);
    const auto *sli_333 = buffer.data(sli + 333);
    const auto *sli_334 = buffer.data(sli + 334);
    const auto *sli_335 = buffer.data(sli + 335);
    const auto *sli_336 = buffer.data(sli + 336);
    const auto *sli_339 = buffer.data(sli + 339);
    const auto *sli_341 = buffer.data(sli + 341);
    const auto *sli_342 = buffer.data(sli + 342);
    const auto *sli_345 = buffer.data(sli + 345);
    const auto *sli_346 = buffer.data(sli + 346);
    const auto *sli_348 = buffer.data(sli + 348);
    const auto *sli_350 = buffer.data(sli + 350);
    const auto *sli_351 = buffer.data(sli + 351);
    const auto *sli_353 = buffer.data(sli + 353);
    const auto *sli_354 = buffer.data(sli + 354);
    const auto *sli_356 = buffer.data(sli + 356);
    const auto *sli_357 = buffer.data(sli + 357);
    const auto *sli_358 = buffer.data(sli + 358);
    const auto *sli_359 = buffer.data(sli + 359);

    const auto *slk1_216 = buffer.data(slk1 + 216);
    const auto *slk1_219 = buffer.data(slk1 + 219);
    const auto *slk1_222 = buffer.data(slk1 + 222);
    const auto *slk1_226 = buffer.data(slk1 + 226);
    const auto *slk1_228 = buffer.data(slk1 + 228);
    const auto *slk1_231 = buffer.data(slk1 + 231);
    const auto *slk1_233 = buffer.data(slk1 + 233);
    const auto *slk1_234 = buffer.data(slk1 + 234);
    const auto *slk1_244 = buffer.data(slk1 + 244);

    const auto *smh0_204 = buffer.data(smh0 + 204);
    const auto *smh0_206 = buffer.data(smh0 + 206);
    const auto *smh0_207 = buffer.data(smh0 + 207);
    const auto *smh0_208 = buffer.data(smh0 + 208);
    const auto *smh0_209 = buffer.data(smh0 + 209);
    const auto *smh0_210 = buffer.data(smh0 + 210);
    const auto *smh0_213 = buffer.data(smh0 + 213);
    const auto *smh0_215 = buffer.data(smh0 + 215);
    const auto *smh0_216 = buffer.data(smh0 + 216);
    const auto *smh0_219 = buffer.data(smh0 + 219);
    const auto *smh0_220 = buffer.data(smh0 + 220);
    const auto *smh0_222 = buffer.data(smh0 + 222);
    const auto *smh0_224 = buffer.data(smh0 + 224);
    const auto *smh0_225 = buffer.data(smh0 + 225);
    const auto *smh0_227 = buffer.data(smh0 + 227);
    const auto *smh0_228 = buffer.data(smh0 + 228);
    const auto *smh0_229 = buffer.data(smh0 + 229);
    const auto *smh0_230 = buffer.data(smh0 + 230);
    const auto *smh0_236 = buffer.data(smh0 + 236);
    const auto *smh0_240 = buffer.data(smh0 + 240);
    const auto *smh0_245 = buffer.data(smh0 + 245);
    const auto *smh0_248 = buffer.data(smh0 + 248);
    const auto *smh0_249 = buffer.data(smh0 + 249);
    const auto *smh0_250 = buffer.data(smh0 + 250);
    const auto *smh0_251 = buffer.data(smh0 + 251);
    const auto *smh0_252 = buffer.data(smh0 + 252);
    const auto *smh0_255 = buffer.data(smh0 + 255);
    const auto *smh0_257 = buffer.data(smh0 + 257);
    const auto *smh0_258 = buffer.data(smh0 + 258);
    const auto *smh0_261 = buffer.data(smh0 + 261);
    const auto *smh0_262 = buffer.data(smh0 + 262);
    const auto *smh0_264 = buffer.data(smh0 + 264);
    const auto *smh0_266 = buffer.data(smh0 + 266);
    const auto *smh0_267 = buffer.data(smh0 + 267);
    const auto *smh0_269 = buffer.data(smh0 + 269);
    const auto *smh0_270 = buffer.data(smh0 + 270);
    const auto *smh0_272 = buffer.data(smh0 + 272);

    const auto *smh1_204 = buffer.data(smh1 + 204);
    const auto *smh1_206 = buffer.data(smh1 + 206);
    const auto *smh1_207 = buffer.data(smh1 + 207);
    const auto *smh1_208 = buffer.data(smh1 + 208);
    const auto *smh1_209 = buffer.data(smh1 + 209);
    const auto *smh1_210 = buffer.data(smh1 + 210);
    const auto *smh1_213 = buffer.data(smh1 + 213);
    const auto *smh1_215 = buffer.data(smh1 + 215);
    const auto *smh1_216 = buffer.data(smh1 + 216);
    const auto *smh1_219 = buffer.data(smh1 + 219);
    const auto *smh1_220 = buffer.data(smh1 + 220);
    const auto *smh1_222 = buffer.data(smh1 + 222);
    const auto *smh1_224 = buffer.data(smh1 + 224);
    const auto *smh1_225 = buffer.data(smh1 + 225);
    const auto *smh1_227 = buffer.data(smh1 + 227);
    const auto *smh1_228 = buffer.data(smh1 + 228);
    const auto *smh1_229 = buffer.data(smh1 + 229);
    const auto *smh1_230 = buffer.data(smh1 + 230);
    const auto *smh1_236 = buffer.data(smh1 + 236);
    const auto *smh1_240 = buffer.data(smh1 + 240);
    const auto *smh1_245 = buffer.data(smh1 + 245);
    const auto *smh1_248 = buffer.data(smh1 + 248);
    const auto *smh1_249 = buffer.data(smh1 + 249);
    const auto *smh1_250 = buffer.data(smh1 + 250);
    const auto *smh1_251 = buffer.data(smh1 + 251);
    const auto *smh1_252 = buffer.data(smh1 + 252);
    const auto *smh1_255 = buffer.data(smh1 + 255);
    const auto *smh1_257 = buffer.data(smh1 + 257);
    const auto *smh1_258 = buffer.data(smh1 + 258);
    const auto *smh1_261 = buffer.data(smh1 + 261);
    const auto *smh1_262 = buffer.data(smh1 + 262);
    const auto *smh1_264 = buffer.data(smh1 + 264);
    const auto *smh1_266 = buffer.data(smh1 + 266);
    const auto *smh1_267 = buffer.data(smh1 + 267);
    const auto *smh1_269 = buffer.data(smh1 + 269);
    const auto *smh1_270 = buffer.data(smh1 + 270);
    const auto *smh1_272 = buffer.data(smh1 + 272);

    const auto *smi_273 = buffer.data(smi + 273);
    const auto *smi_275 = buffer.data(smi + 275);
    const auto *smi_276 = buffer.data(smi + 276);
    const auto *smi_277 = buffer.data(smi + 277);
    const auto *smi_278 = buffer.data(smi + 278);
    const auto *smi_279 = buffer.data(smi + 279);
    const auto *smi_280 = buffer.data(smi + 280);
    const auto *smi_282 = buffer.data(smi + 282);
    const auto *smi_283 = buffer.data(smi + 283);
    const auto *smi_285 = buffer.data(smi + 285);
    const auto *smi_286 = buffer.data(smi + 286);
    const auto *smi_289 = buffer.data(smi + 289);
    const auto *smi_290 = buffer.data(smi + 290);
    const auto *smi_292 = buffer.data(smi + 292);
    const auto *smi_294 = buffer.data(smi + 294);
    const auto *smi_295 = buffer.data(smi + 295);
    const auto *smi_297 = buffer.data(smi + 297);
    const auto *smi_298 = buffer.data(smi + 298);
    const auto *smi_300 = buffer.data(smi + 300);
    const auto *smi_301 = buffer.data(smi + 301);
    const auto *smi_302 = buffer.data(smi + 302);
    const auto *smi_303 = buffer.data(smi + 303);
    const auto *smi_304 = buffer.data(smi + 304);
    const auto *smi_305 = buffer.data(smi + 305);
    const auto *smi_306 = buffer.data(smi + 306);
    const auto *smi_307 = buffer.data(smi + 307);
    const auto *smi_308 = buffer.data(smi + 308);
    const auto *smi_310 = buffer.data(smi + 310);
    const auto *smi_311 = buffer.data(smi + 311);
    const auto *smi_313 = buffer.data(smi + 313);
    const auto *smi_314 = buffer.data(smi + 314);
    const auto *smi_317 = buffer.data(smi + 317);
    const auto *smi_318 = buffer.data(smi + 318);
    const auto *smi_322 = buffer.data(smi + 322);
    const auto *smi_328 = buffer.data(smi + 328);
    const auto *smi_329 = buffer.data(smi + 329);
    const auto *smi_330 = buffer.data(smi + 330);
    const auto *smi_331 = buffer.data(smi + 331);
    const auto *smi_332 = buffer.data(smi + 332);
    const auto *smi_333 = buffer.data(smi + 333);
    const auto *smi_334 = buffer.data(smi + 334);
    const auto *smi_335 = buffer.data(smi + 335);
    const auto *smi_336 = buffer.data(smi + 336);
    const auto *smi_338 = buffer.data(smi + 338);
    const auto *smi_339 = buffer.data(smi + 339);
    const auto *smi_341 = buffer.data(smi + 341);
    const auto *smi_342 = buffer.data(smi + 342);
    const auto *smi_345 = buffer.data(smi + 345);
    const auto *smi_346 = buffer.data(smi + 346);
    const auto *smi_348 = buffer.data(smi + 348);
    const auto *smi_350 = buffer.data(smi + 350);
    const auto *smi_351 = buffer.data(smi + 351);
    const auto *smi_353 = buffer.data(smi + 353);
    const auto *smi_354 = buffer.data(smi + 354);
    const auto *smi_356 = buffer.data(smi + 356);
    const auto *smi_357 = buffer.data(smi + 357);
    const auto *smi_358 = buffer.data(smi + 358);
    const auto *smi_359 = buffer.data(smi + 359);

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pc_x, sli_275, sli_276, sli_277, \
                         sli_278, sli_279, smi_275, smi_276, smi_277, smi_278, \
                         smi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_20 * sli_275[k]
                   + f_3 * pc_x[k] * smi_275[k];

        t_348[k] = f_20 * sli_276[k]
                   + f_3 * pc_x[k] * smi_276[k];

        t_349[k] = f_20 * sli_277[k]
                   + f_3 * pc_x[k] * smi_277[k];

        t_350[k] = f_20 * sli_278[k]
                   + f_3 * pc_x[k] * smi_278[k];

        t_351[k] = f_20 * sli_279[k]
                   + f_3 * pc_x[k] * smi_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_y, pc_z, sli_161, smh0_204, smh0_206, \
                         smh0_207, smh1_204, smh1_206, smh1_207, smi_273, smi_275, \
                         smi_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * smh0_204[k]
                   - f_2 * smh1_204[k]
                   + f_3 * pc_y[k] * smi_273[k];

        t_353[k] = f_15 * sli_161[k]
                   + f_3 * pc_z[k] * smi_273[k];

        t_354[k] = f_4 * smh0_206[k]
                   - f_5 * smh1_206[k]
                   + f_3 * pc_y[k] * smi_275[k];

        t_355[k] = f_6 * smh0_207[k]
                   - f_7 * smh1_207[k]
                   + f_3 * pc_y[k] * smi_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, sli_167, smh0_208, smh0_209, \
                         smh1_208, smh1_209, smi_277, smi_278, \
                         smi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_8 * smh0_208[k]
                   - f_9 * smh1_208[k]
                   + f_3 * pc_y[k] * smi_277[k];

        t_357[k] = f_10 * smh0_209[k]
                   - f_11 * smh1_209[k]
                   + f_3 * pc_y[k] * smi_278[k];

        t_358[k] = f_3 * pc_y[k] * smi_279[k];

        t_359[k] = f_15 * sli_167[k]
                   + f_1 * smh0_209[k]
                   - f_2 * smh1_209[k]
                   + f_3 * pc_z[k] * smi_279[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, sli_168, sli_280, \
                         sli_283, smh0_210, smh0_213, smh1_210, smh1_213, smi_280, \
                         smi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_17 * sli_280[k]
                   + f_1 * smh0_210[k]
                   - f_2 * smh1_210[k]
                   + f_3 * pc_x[k] * smi_280[k];

        t_361[k] = f_16 * sli_168[k]
                   + f_3 * pc_y[k] * smi_280[k];

        t_362[k] = f_3 * pc_z[k] * smi_280[k];

        t_363[k] = f_17 * sli_283[k]
                   + f_4 * smh0_213[k]
                   - f_5 * smh1_213[k]
                   + f_3 * pc_x[k] * smi_283[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pc_x, pc_y, sli_170, sli_285, sli_286, smh0_215, \
                         smh0_216, smh1_215, smh1_216, smi_282, smi_285, \
                         smi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_16 * sli_170[k]
                   + f_3 * pc_y[k] * smi_282[k];

        t_365[k] = f_17 * sli_285[k]
                   + f_4 * smh0_215[k]
                   - f_5 * smh1_215[k]
                   + f_3 * pc_x[k] * smi_285[k];

        t_366[k] = f_17 * sli_286[k]
                   + f_6 * smh0_216[k]
                   - f_7 * smh1_216[k]
                   + f_3 * pc_x[k] * smi_286[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pc_x, pc_y, pc_z, sli_173, sli_289, smh0_219, \
                         smh1_219, smi_283, smi_285, smi_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_3 * pc_z[k] * smi_283[k];

        t_368[k] = f_16 * sli_173[k]
                   + f_3 * pc_y[k] * smi_285[k];

        t_369[k] = f_17 * sli_289[k]
                   + f_6 * smh0_219[k]
                   - f_7 * smh1_219[k]
                   + f_3 * pc_x[k] * smi_289[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pc_x, pc_z, sli_290, sli_292, smh0_220, \
                         smh0_222, smh1_220, smh1_222, smi_286, smi_290, \
                         smi_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_17 * sli_290[k]
                   + f_8 * smh0_220[k]
                   - f_9 * smh1_220[k]
                   + f_3 * pc_x[k] * smi_290[k];

        t_371[k] = f_3 * pc_z[k] * smi_286[k];

        t_372[k] = f_17 * sli_292[k]
                   + f_8 * smh0_222[k]
                   - f_9 * smh1_222[k]
                   + f_3 * pc_x[k] * smi_292[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, sli_177, sli_294, sli_295, smh0_224, \
                         smh0_225, smh1_224, smh1_225, smi_289, smi_294, \
                         smi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * sli_177[k]
                   + f_3 * pc_y[k] * smi_289[k];

        t_374[k] = f_17 * sli_294[k]
                   + f_8 * smh0_224[k]
                   - f_9 * smh1_224[k]
                   + f_3 * pc_x[k] * smi_294[k];

        t_375[k] = f_17 * sli_295[k]
                   + f_10 * smh0_225[k]
                   - f_11 * smh1_225[k]
                   + f_3 * pc_x[k] * smi_295[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_z, sli_297, sli_298, smh0_227, \
                         smh0_228, smh1_227, smh1_228, smi_290, smi_297, \
                         smi_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_3 * pc_z[k] * smi_290[k];

        t_377[k] = f_17 * sli_297[k]
                   + f_10 * smh0_227[k]
                   - f_11 * smh1_227[k]
                   + f_3 * pc_x[k] * smi_297[k];

        t_378[k] = f_17 * sli_298[k]
                   + f_10 * smh0_228[k]
                   - f_11 * smh1_228[k]
                   + f_3 * pc_x[k] * smi_298[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, sli_182, sli_300, sli_301, \
                         sli_302, smh0_230, smh1_230, smi_294, smi_300, smi_301, \
                         smi_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_16 * sli_182[k]
                   + f_3 * pc_y[k] * smi_294[k];

        t_380[k] = f_17 * sli_300[k]
                   + f_10 * smh0_230[k]
                   - f_11 * smh1_230[k]
                   + f_3 * pc_x[k] * smi_300[k];

        t_381[k] = f_17 * sli_301[k]
                   + f_3 * pc_x[k] * smi_301[k];

        t_382[k] = f_17 * sli_302[k]
                   + f_3 * pc_x[k] * smi_302[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pc_x, sli_303, sli_304, sli_305, \
                         sli_306, sli_307, smi_303, smi_304, smi_305, smi_306, \
                         smi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_17 * sli_303[k]
                   + f_3 * pc_x[k] * smi_303[k];

        t_384[k] = f_17 * sli_304[k]
                   + f_3 * pc_x[k] * smi_304[k];

        t_385[k] = f_17 * sli_305[k]
                   + f_3 * pc_x[k] * smi_305[k];

        t_386[k] = f_17 * sli_306[k]
                   + f_3 * pc_x[k] * smi_306[k];

        t_387[k] = f_17 * sli_307[k]
                   + f_3 * pc_x[k] * smi_307[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, pc_z, sli_189, sli_191, smh0_225, \
                         smh0_227, smh1_225, smh1_227, smi_301, \
                         smi_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_16 * sli_189[k]
                   + f_1 * smh0_225[k]
                   - f_2 * smh1_225[k]
                   + f_3 * pc_y[k] * smi_301[k];

        t_389[k] = f_3 * pc_z[k] * smi_301[k];

        t_390[k] = f_16 * sli_191[k]
                   + f_4 * smh0_227[k]
                   - f_5 * smh1_227[k]
                   + f_3 * pc_y[k] * smi_303[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_y, sli_192, sli_193, sli_194, smh0_228, \
                         smh0_229, smh0_230, smh1_228, smh1_229, smh1_230, smi_304, smi_305, \
                         smi_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_16 * sli_192[k]
                   + f_6 * smh0_228[k]
                   - f_7 * smh1_228[k]
                   + f_3 * pc_y[k] * smi_304[k];

        t_392[k] = f_16 * sli_193[k]
                   + f_8 * smh0_229[k]
                   - f_9 * smh1_229[k]
                   + f_3 * pc_y[k] * smi_305[k];

        t_393[k] = f_16 * sli_194[k]
                   + f_10 * smh0_230[k]
                   - f_11 * smh1_230[k]
                   + f_3 * pc_y[k] * smi_306[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pb_z, pc_y, pc_z, slk0_216, sli_195, \
                         sli_196, slk1_216, smh0_230, smh1_230, smi_307, \
                         smi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * sli_195[k]
                   + f_3 * pc_y[k] * smi_307[k];

        t_395[k] = f_1 * smh0_230[k]
                   - f_2 * smh1_230[k]
                   + f_3 * pc_z[k] * smi_307[k];

        t_396[k] = pb_z[k] * slk0_216[k]
                   - f_12 * pc_z[k] * slk1_216[k];

        t_397[k] = f_15 * sli_196[k]
                   + f_3 * pc_y[k] * smi_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pb_z, pc_y, pc_z, slk0_219, sli_168, sli_198, \
                         slk1_219, smi_308, smi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * sli_168[k]
                   + f_3 * pc_z[k] * smi_308[k];

        t_399[k] = pb_z[k] * slk0_219[k]
                   - f_12 * pc_z[k] * slk1_219[k];

        t_400[k] = f_15 * sli_198[k]
                   + f_3 * pc_y[k] * smi_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pb_z, pc_x, pc_z, slk0_222, sli_171, sli_313, \
                         slk1_222, smh0_236, smh1_236, smi_311, \
                         smi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_17 * sli_313[k]
                   + f_4 * smh0_236[k]
                   - f_5 * smh1_236[k]
                   + f_3 * pc_x[k] * smi_313[k];

        t_402[k] = pb_z[k] * slk0_222[k]
                   - f_12 * pc_z[k] * slk1_222[k];

        t_403[k] = f_13 * sli_171[k]
                   + f_3 * pc_z[k] * smi_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pb_z, pc_x, pc_y, pc_z, slk0_226, sli_201, \
                         sli_317, slk1_226, smh0_240, smh1_240, smi_313, \
                         smi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * sli_201[k]
                   + f_3 * pc_y[k] * smi_313[k];

        t_405[k] = f_17 * sli_317[k]
                   + f_6 * smh0_240[k]
                   - f_7 * smh1_240[k]
                   + f_3 * pc_x[k] * smi_317[k];

        t_406[k] = pb_z[k] * slk0_226[k]
                   - f_12 * pc_z[k] * slk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_z, pc_y, pc_z, slk0_228, sli_174, sli_175, \
                         sli_205, slk1_228, smi_314, smi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * sli_174[k]
                   + f_3 * pc_z[k] * smi_314[k];

        t_408[k] = pb_z[k] * slk0_228[k]
                   + f_14 * sli_175[k]
                   - f_12 * pc_z[k] * slk1_228[k];

        t_409[k] = f_15 * sli_205[k]
                   + f_3 * pc_y[k] * smi_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pb_z, pc_x, pc_z, slk0_231, sli_178, sli_322, \
                         slk1_231, smh0_245, smh1_245, smi_318, \
                         smi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_17 * sli_322[k]
                   + f_8 * smh0_245[k]
                   - f_9 * smh1_245[k]
                   + f_3 * pc_x[k] * smi_322[k];

        t_411[k] = pb_z[k] * slk0_231[k]
                   - f_12 * pc_z[k] * slk1_231[k];

        t_412[k] = f_13 * sli_178[k]
                   + f_3 * pc_z[k] * smi_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pb_z, pc_y, pc_z, slk0_233, slk0_234, sli_179, \
                         sli_180, sli_210, slk1_233, slk1_234, \
                         smi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_z[k] * slk0_233[k]
                   + f_14 * sli_179[k]
                   - f_12 * pc_z[k] * slk1_233[k];

        t_414[k] = pb_z[k] * slk0_234[k]
                   + f_15 * sli_180[k]
                   - f_12 * pc_z[k] * slk1_234[k];

        t_415[k] = f_15 * sli_210[k]
                   + f_3 * pc_y[k] * smi_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, sli_328, sli_329, sli_330, sli_331, \
                         smh0_251, smh1_251, smi_328, smi_329, smi_330, \
                         smi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_17 * sli_328[k]
                   + f_10 * smh0_251[k]
                   - f_11 * smh1_251[k]
                   + f_3 * pc_x[k] * smi_328[k];

        t_417[k] = f_17 * sli_329[k]
                   + f_3 * pc_x[k] * smi_329[k];

        t_418[k] = f_17 * sli_330[k]
                   + f_3 * pc_x[k] * smi_330[k];

        t_419[k] = f_17 * sli_331[k]
                   + f_3 * pc_x[k] * smi_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, sli_332, sli_333, sli_334, sli_335, \
                         smi_332, smi_333, smi_334, smi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_17 * sli_332[k]
                   + f_3 * pc_x[k] * smi_332[k];

        t_421[k] = f_17 * sli_333[k]
                   + f_3 * pc_x[k] * smi_333[k];

        t_422[k] = f_17 * sli_334[k]
                   + f_3 * pc_x[k] * smi_334[k];

        t_423[k] = f_17 * sli_335[k]
                   + f_3 * pc_x[k] * smi_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_z, pc_y, pc_z, slk0_244, sli_189, sli_219, \
                         slk1_244, smh0_248, smh1_248, smi_329, \
                         smi_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_z[k] * slk0_244[k]
                   - f_12 * pc_z[k] * slk1_244[k];

        t_425[k] = f_13 * sli_189[k]
                   + f_3 * pc_z[k] * smi_329[k];

        t_426[k] = f_15 * sli_219[k]
                   + f_4 * smh0_248[k]
                   - f_5 * smh1_248[k]
                   + f_3 * pc_y[k] * smi_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, sli_220, sli_221, sli_222, smh0_249, \
                         smh0_250, smh0_251, smh1_249, smh1_250, smh1_251, smi_332, smi_333, \
                         smi_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * sli_220[k]
                   + f_6 * smh0_249[k]
                   - f_7 * smh1_249[k]
                   + f_3 * pc_y[k] * smi_332[k];

        t_428[k] = f_15 * sli_221[k]
                   + f_8 * smh0_250[k]
                   - f_9 * smh1_250[k]
                   + f_3 * pc_y[k] * smi_333[k];

        t_429[k] = f_15 * sli_222[k]
                   + f_10 * smh0_251[k]
                   - f_11 * smh1_251[k]
                   + f_3 * pc_y[k] * smi_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, sli_195, sli_223, sli_336, \
                         smh0_251, smh0_252, smh1_251, smh1_252, smi_335, \
                         smi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * sli_223[k]
                   + f_3 * pc_y[k] * smi_335[k];

        t_431[k] = f_13 * sli_195[k]
                   + f_1 * smh0_251[k]
                   - f_2 * smh1_251[k]
                   + f_3 * pc_z[k] * smi_335[k];

        t_432[k] = f_17 * sli_336[k]
                   + f_1 * smh0_252[k]
                   - f_2 * smh1_252[k]
                   + f_3 * pc_x[k] * smi_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, sli_196, sli_224, \
                         sli_226, sli_339, smh0_255, smh1_255, smi_336, smi_338, \
                         smi_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * sli_224[k]
                   + f_3 * pc_y[k] * smi_336[k];

        t_434[k] = f_14 * sli_196[k]
                   + f_3 * pc_z[k] * smi_336[k];

        t_435[k] = f_17 * sli_339[k]
                   + f_4 * smh0_255[k]
                   - f_5 * smh1_255[k]
                   + f_3 * pc_x[k] * smi_339[k];

        t_436[k] = f_14 * sli_226[k]
                   + f_3 * pc_y[k] * smi_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, sli_199, sli_341, sli_342, smh0_257, \
                         smh0_258, smh1_257, smh1_258, smi_339, smi_341, \
                         smi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_17 * sli_341[k]
                   + f_4 * smh0_257[k]
                   - f_5 * smh1_257[k]
                   + f_3 * pc_x[k] * smi_341[k];

        t_438[k] = f_17 * sli_342[k]
                   + f_6 * smh0_258[k]
                   - f_7 * smh1_258[k]
                   + f_3 * pc_x[k] * smi_342[k];

        t_439[k] = f_14 * sli_199[k]
                   + f_3 * pc_z[k] * smi_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, sli_229, sli_345, sli_346, smh0_261, \
                         smh0_262, smh1_261, smh1_262, smi_341, smi_345, \
                         smi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * sli_229[k]
                   + f_3 * pc_y[k] * smi_341[k];

        t_441[k] = f_17 * sli_345[k]
                   + f_6 * smh0_261[k]
                   - f_7 * smh1_261[k]
                   + f_3 * pc_x[k] * smi_345[k];

        t_442[k] = f_17 * sli_346[k]
                   + f_8 * smh0_262[k]
                   - f_9 * smh1_262[k]
                   + f_3 * pc_x[k] * smi_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, sli_202, sli_233, sli_348, \
                         smh0_264, smh1_264, smi_342, smi_345, \
                         smi_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * sli_202[k]
                   + f_3 * pc_z[k] * smi_342[k];

        t_444[k] = f_17 * sli_348[k]
                   + f_8 * smh0_264[k]
                   - f_9 * smh1_264[k]
                   + f_3 * pc_x[k] * smi_348[k];

        t_445[k] = f_14 * sli_233[k]
                   + f_3 * pc_y[k] * smi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, sli_206, sli_350, sli_351, smh0_266, \
                         smh0_267, smh1_266, smh1_267, smi_346, smi_350, \
                         smi_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_17 * sli_350[k]
                   + f_8 * smh0_266[k]
                   - f_9 * smh1_266[k]
                   + f_3 * pc_x[k] * smi_350[k];

        t_447[k] = f_17 * sli_351[k]
                   + f_10 * smh0_267[k]
                   - f_11 * smh1_267[k]
                   + f_3 * pc_x[k] * smi_351[k];

        t_448[k] = f_14 * sli_206[k]
                   + f_3 * pc_z[k] * smi_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, sli_238, sli_353, sli_354, smh0_269, \
                         smh0_270, smh1_269, smh1_270, smi_350, smi_353, \
                         smi_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * sli_353[k]
                   + f_10 * smh0_269[k]
                   - f_11 * smh1_269[k]
                   + f_3 * pc_x[k] * smi_353[k];

        t_450[k] = f_17 * sli_354[k]
                   + f_10 * smh0_270[k]
                   - f_11 * smh1_270[k]
                   + f_3 * pc_x[k] * smi_354[k];

        t_451[k] = f_14 * sli_238[k]
                   + f_3 * pc_y[k] * smi_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, sli_356, sli_357, sli_358, sli_359, \
                         smh0_272, smh1_272, smi_356, smi_357, smi_358, \
                         smi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_17 * sli_356[k]
                   + f_10 * smh0_272[k]
                   - f_11 * smh1_272[k]
                   + f_3 * pc_x[k] * smi_356[k];

        t_453[k] = f_17 * sli_357[k]
                   + f_3 * pc_x[k] * smi_357[k];

        t_454[k] = f_17 * sli_358[k]
                   + f_3 * pc_x[k] * smi_358[k];

        t_455[k] = f_17 * sli_359[k]
                   + f_3 * pc_x[k] * smi_359[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_324 = buffer.data(slk0 + 324);
    const auto *slk0_327 = buffer.data(slk0 + 327);
    const auto *slk0_329 = buffer.data(slk0 + 329);
    const auto *slk0_330 = buffer.data(slk0 + 330);
    const auto *slk0_333 = buffer.data(slk0 + 333);
    const auto *slk0_334 = buffer.data(slk0 + 334);
    const auto *slk0_336 = buffer.data(slk0 + 336);
    const auto *slk0_338 = buffer.data(slk0 + 338);
    const auto *slk0_339 = buffer.data(slk0 + 339);
    const auto *slk0_341 = buffer.data(slk0 + 341);
    const auto *slk0_342 = buffer.data(slk0 + 342);
    const auto *slk0_344 = buffer.data(slk0 + 344);
    const auto *slk0_359 = buffer.data(slk0 + 359);

    const auto *sli_217 = buffer.data(sli + 217);
    const auto *sli_223 = buffer.data(sli + 223);
    const auto *sli_224 = buffer.data(sli + 224);
    const auto *sli_227 = buffer.data(sli + 227);
    const auto *sli_230 = buffer.data(sli + 230);
    const auto *sli_234 = buffer.data(sli + 234);
    const auto *sli_245 = buffer.data(sli + 245);
    const auto *sli_247 = buffer.data(sli + 247);
    const auto *sli_248 = buffer.data(sli + 248);
    const auto *sli_249 = buffer.data(sli + 249);
    const auto *sli_250 = buffer.data(sli + 250);
    const auto *sli_251 = buffer.data(sli + 251);
    const auto *sli_252 = buffer.data(sli + 252);
    const auto *sli_253 = buffer.data(sli + 253);
    const auto *sli_254 = buffer.data(sli + 254);
    const auto *sli_255 = buffer.data(sli + 255);
    const auto *sli_257 = buffer.data(sli + 257);
    const auto *sli_258 = buffer.data(sli + 258);
    const auto *sli_260 = buffer.data(sli + 260);
    const auto *sli_261 = buffer.data(sli + 261);
    const auto *sli_262 = buffer.data(sli + 262);
    const auto *sli_264 = buffer.data(sli + 264);
    const auto *sli_265 = buffer.data(sli + 265);
    const auto *sli_266 = buffer.data(sli + 266);
    const auto *sli_273 = buffer.data(sli + 273);
    const auto *sli_275 = buffer.data(sli + 275);
    const auto *sli_276 = buffer.data(sli + 276);
    const auto *sli_277 = buffer.data(sli + 277);
    const auto *sli_278 = buffer.data(sli + 278);
    const auto *sli_279 = buffer.data(sli + 279);
    const auto *sli_280 = buffer.data(sli + 280);
    const auto *sli_282 = buffer.data(sli + 282);
    const auto *sli_285 = buffer.data(sli + 285);
    const auto *sli_289 = buffer.data(sli + 289);
    const auto *sli_294 = buffer.data(sli + 294);
    const auto *sli_301 = buffer.data(sli + 301);
    const auto *sli_303 = buffer.data(sli + 303);
    const auto *sli_360 = buffer.data(sli + 360);
    const auto *sli_361 = buffer.data(sli + 361);
    const auto *sli_362 = buffer.data(sli + 362);
    const auto *sli_363 = buffer.data(sli + 363);
    const auto *sli_385 = buffer.data(sli + 385);
    const auto *sli_386 = buffer.data(sli + 386);
    const auto *sli_387 = buffer.data(sli + 387);
    const auto *sli_388 = buffer.data(sli + 388);
    const auto *sli_389 = buffer.data(sli + 389);
    const auto *sli_390 = buffer.data(sli + 390);
    const auto *sli_391 = buffer.data(sli + 391);
    const auto *sli_392 = buffer.data(sli + 392);
    const auto *sli_395 = buffer.data(sli + 395);
    const auto *sli_397 = buffer.data(sli + 397);
    const auto *sli_398 = buffer.data(sli + 398);
    const auto *sli_401 = buffer.data(sli + 401);
    const auto *sli_402 = buffer.data(sli + 402);
    const auto *sli_404 = buffer.data(sli + 404);
    const auto *sli_406 = buffer.data(sli + 406);
    const auto *sli_407 = buffer.data(sli + 407);
    const auto *sli_409 = buffer.data(sli + 409);
    const auto *sli_410 = buffer.data(sli + 410);
    const auto *sli_412 = buffer.data(sli + 412);
    const auto *sli_413 = buffer.data(sli + 413);
    const auto *sli_414 = buffer.data(sli + 414);
    const auto *sli_415 = buffer.data(sli + 415);
    const auto *sli_416 = buffer.data(sli + 416);
    const auto *sli_417 = buffer.data(sli + 417);
    const auto *sli_418 = buffer.data(sli + 418);
    const auto *sli_419 = buffer.data(sli + 419);
    const auto *sli_420 = buffer.data(sli + 420);
    const auto *sli_423 = buffer.data(sli + 423);
    const auto *sli_425 = buffer.data(sli + 425);
    const auto *sli_426 = buffer.data(sli + 426);
    const auto *sli_429 = buffer.data(sli + 429);
    const auto *sli_430 = buffer.data(sli + 430);
    const auto *sli_432 = buffer.data(sli + 432);
    const auto *sli_434 = buffer.data(sli + 434);
    const auto *sli_435 = buffer.data(sli + 435);
    const auto *sli_437 = buffer.data(sli + 437);
    const auto *sli_438 = buffer.data(sli + 438);
    const auto *sli_440 = buffer.data(sli + 440);
    const auto *sli_441 = buffer.data(sli + 441);
    const auto *sli_442 = buffer.data(sli + 442);
    const auto *sli_443 = buffer.data(sli + 443);
    const auto *sli_444 = buffer.data(sli + 444);
    const auto *sli_445 = buffer.data(sli + 445);
    const auto *sli_446 = buffer.data(sli + 446);
    const auto *sli_447 = buffer.data(sli + 447);

    const auto *slk1_324 = buffer.data(slk1 + 324);
    const auto *slk1_327 = buffer.data(slk1 + 327);
    const auto *slk1_329 = buffer.data(slk1 + 329);
    const auto *slk1_330 = buffer.data(slk1 + 330);
    const auto *slk1_333 = buffer.data(slk1 + 333);
    const auto *slk1_334 = buffer.data(slk1 + 334);
    const auto *slk1_336 = buffer.data(slk1 + 336);
    const auto *slk1_338 = buffer.data(slk1 + 338);
    const auto *slk1_339 = buffer.data(slk1 + 339);
    const auto *slk1_341 = buffer.data(slk1 + 341);
    const auto *slk1_342 = buffer.data(slk1 + 342);
    const auto *slk1_344 = buffer.data(slk1 + 344);
    const auto *slk1_359 = buffer.data(slk1 + 359);

    const auto *smh0_267 = buffer.data(smh0 + 267);
    const auto *smh0_269 = buffer.data(smh0 + 269);
    const auto *smh0_270 = buffer.data(smh0 + 270);
    const auto *smh0_271 = buffer.data(smh0 + 271);
    const auto *smh0_272 = buffer.data(smh0 + 272);
    const auto *smh0_288 = buffer.data(smh0 + 288);
    const auto *smh0_290 = buffer.data(smh0 + 290);
    const auto *smh0_291 = buffer.data(smh0 + 291);
    const auto *smh0_292 = buffer.data(smh0 + 292);
    const auto *smh0_293 = buffer.data(smh0 + 293);
    const auto *smh0_294 = buffer.data(smh0 + 294);
    const auto *smh0_297 = buffer.data(smh0 + 297);
    const auto *smh0_299 = buffer.data(smh0 + 299);
    const auto *smh0_300 = buffer.data(smh0 + 300);
    const auto *smh0_303 = buffer.data(smh0 + 303);
    const auto *smh0_304 = buffer.data(smh0 + 304);
    const auto *smh0_306 = buffer.data(smh0 + 306);
    const auto *smh0_308 = buffer.data(smh0 + 308);
    const auto *smh0_309 = buffer.data(smh0 + 309);
    const auto *smh0_311 = buffer.data(smh0 + 311);
    const auto *smh0_312 = buffer.data(smh0 + 312);
    const auto *smh0_313 = buffer.data(smh0 + 313);
    const auto *smh0_314 = buffer.data(smh0 + 314);
    const auto *smh0_315 = buffer.data(smh0 + 315);
    const auto *smh0_318 = buffer.data(smh0 + 318);
    const auto *smh0_320 = buffer.data(smh0 + 320);
    const auto *smh0_321 = buffer.data(smh0 + 321);
    const auto *smh0_324 = buffer.data(smh0 + 324);
    const auto *smh0_325 = buffer.data(smh0 + 325);
    const auto *smh0_327 = buffer.data(smh0 + 327);
    const auto *smh0_329 = buffer.data(smh0 + 329);
    const auto *smh0_330 = buffer.data(smh0 + 330);
    const auto *smh0_332 = buffer.data(smh0 + 332);
    const auto *smh0_333 = buffer.data(smh0 + 333);
    const auto *smh0_335 = buffer.data(smh0 + 335);

    const auto *smh1_267 = buffer.data(smh1 + 267);
    const auto *smh1_269 = buffer.data(smh1 + 269);
    const auto *smh1_270 = buffer.data(smh1 + 270);
    const auto *smh1_271 = buffer.data(smh1 + 271);
    const auto *smh1_272 = buffer.data(smh1 + 272);
    const auto *smh1_288 = buffer.data(smh1 + 288);
    const auto *smh1_290 = buffer.data(smh1 + 290);
    const auto *smh1_291 = buffer.data(smh1 + 291);
    const auto *smh1_292 = buffer.data(smh1 + 292);
    const auto *smh1_293 = buffer.data(smh1 + 293);
    const auto *smh1_294 = buffer.data(smh1 + 294);
    const auto *smh1_297 = buffer.data(smh1 + 297);
    const auto *smh1_299 = buffer.data(smh1 + 299);
    const auto *smh1_300 = buffer.data(smh1 + 300);
    const auto *smh1_303 = buffer.data(smh1 + 303);
    const auto *smh1_304 = buffer.data(smh1 + 304);
    const auto *smh1_306 = buffer.data(smh1 + 306);
    const auto *smh1_308 = buffer.data(smh1 + 308);
    const auto *smh1_309 = buffer.data(smh1 + 309);
    const auto *smh1_311 = buffer.data(smh1 + 311);
    const auto *smh1_312 = buffer.data(smh1 + 312);
    const auto *smh1_313 = buffer.data(smh1 + 313);
    const auto *smh1_314 = buffer.data(smh1 + 314);
    const auto *smh1_315 = buffer.data(smh1 + 315);
    const auto *smh1_318 = buffer.data(smh1 + 318);
    const auto *smh1_320 = buffer.data(smh1 + 320);
    const auto *smh1_321 = buffer.data(smh1 + 321);
    const auto *smh1_324 = buffer.data(smh1 + 324);
    const auto *smh1_325 = buffer.data(smh1 + 325);
    const auto *smh1_327 = buffer.data(smh1 + 327);
    const auto *smh1_329 = buffer.data(smh1 + 329);
    const auto *smh1_330 = buffer.data(smh1 + 330);
    const auto *smh1_332 = buffer.data(smh1 + 332);
    const auto *smh1_333 = buffer.data(smh1 + 333);
    const auto *smh1_335 = buffer.data(smh1 + 335);

    const auto *smi_357 = buffer.data(smi + 357);
    const auto *smi_359 = buffer.data(smi + 359);
    const auto *smi_360 = buffer.data(smi + 360);
    const auto *smi_361 = buffer.data(smi + 361);
    const auto *smi_362 = buffer.data(smi + 362);
    const auto *smi_363 = buffer.data(smi + 363);
    const auto *smi_364 = buffer.data(smi + 364);
    const auto *smi_366 = buffer.data(smi + 366);
    const auto *smi_367 = buffer.data(smi + 367);
    const auto *smi_369 = buffer.data(smi + 369);
    const auto *smi_370 = buffer.data(smi + 370);
    const auto *smi_373 = buffer.data(smi + 373);
    const auto *smi_374 = buffer.data(smi + 374);
    const auto *smi_378 = buffer.data(smi + 378);
    const auto *smi_385 = buffer.data(smi + 385);
    const auto *smi_386 = buffer.data(smi + 386);
    const auto *smi_387 = buffer.data(smi + 387);
    const auto *smi_388 = buffer.data(smi + 388);
    const auto *smi_389 = buffer.data(smi + 389);
    const auto *smi_390 = buffer.data(smi + 390);
    const auto *smi_391 = buffer.data(smi + 391);
    const auto *smi_392 = buffer.data(smi + 392);
    const auto *smi_394 = buffer.data(smi + 394);
    const auto *smi_395 = buffer.data(smi + 395);
    const auto *smi_397 = buffer.data(smi + 397);
    const auto *smi_398 = buffer.data(smi + 398);
    const auto *smi_401 = buffer.data(smi + 401);
    const auto *smi_402 = buffer.data(smi + 402);
    const auto *smi_404 = buffer.data(smi + 404);
    const auto *smi_406 = buffer.data(smi + 406);
    const auto *smi_407 = buffer.data(smi + 407);
    const auto *smi_409 = buffer.data(smi + 409);
    const auto *smi_410 = buffer.data(smi + 410);
    const auto *smi_412 = buffer.data(smi + 412);
    const auto *smi_413 = buffer.data(smi + 413);
    const auto *smi_414 = buffer.data(smi + 414);
    const auto *smi_415 = buffer.data(smi + 415);
    const auto *smi_416 = buffer.data(smi + 416);
    const auto *smi_417 = buffer.data(smi + 417);
    const auto *smi_418 = buffer.data(smi + 418);
    const auto *smi_419 = buffer.data(smi + 419);
    const auto *smi_420 = buffer.data(smi + 420);
    const auto *smi_422 = buffer.data(smi + 422);
    const auto *smi_423 = buffer.data(smi + 423);
    const auto *smi_425 = buffer.data(smi + 425);
    const auto *smi_426 = buffer.data(smi + 426);
    const auto *smi_429 = buffer.data(smi + 429);
    const auto *smi_430 = buffer.data(smi + 430);
    const auto *smi_432 = buffer.data(smi + 432);
    const auto *smi_434 = buffer.data(smi + 434);
    const auto *smi_435 = buffer.data(smi + 435);
    const auto *smi_437 = buffer.data(smi + 437);
    const auto *smi_438 = buffer.data(smi + 438);
    const auto *smi_440 = buffer.data(smi + 440);
    const auto *smi_441 = buffer.data(smi + 441);
    const auto *smi_442 = buffer.data(smi + 442);
    const auto *smi_443 = buffer.data(smi + 443);
    const auto *smi_444 = buffer.data(smi + 444);
    const auto *smi_445 = buffer.data(smi + 445);
    const auto *smi_446 = buffer.data(smi + 446);
    const auto *smi_447 = buffer.data(smi + 447);

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, sli_360, sli_361, sli_362, sli_363, \
                         smi_360, smi_361, smi_362, smi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_17 * sli_360[k]
                   + f_3 * pc_x[k] * smi_360[k];

        t_457[k] = f_17 * sli_361[k]
                   + f_3 * pc_x[k] * smi_361[k];

        t_458[k] = f_17 * sli_362[k]
                   + f_3 * pc_x[k] * smi_362[k];

        t_459[k] = f_17 * sli_363[k]
                   + f_3 * pc_x[k] * smi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, sli_217, sli_245, sli_247, smh0_267, \
                         smh0_269, smh1_267, smh1_269, smi_357, \
                         smi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * sli_245[k]
                   + f_1 * smh0_267[k]
                   - f_2 * smh1_267[k]
                   + f_3 * pc_y[k] * smi_357[k];

        t_461[k] = f_14 * sli_217[k]
                   + f_3 * pc_z[k] * smi_357[k];

        t_462[k] = f_14 * sli_247[k]
                   + f_4 * smh0_269[k]
                   - f_5 * smh1_269[k]
                   + f_3 * pc_y[k] * smi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, sli_248, sli_249, sli_250, smh0_270, \
                         smh0_271, smh0_272, smh1_270, smh1_271, smh1_272, smi_360, smi_361, \
                         smi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * sli_248[k]
                   + f_6 * smh0_270[k]
                   - f_7 * smh1_270[k]
                   + f_3 * pc_y[k] * smi_360[k];

        t_464[k] = f_14 * sli_249[k]
                   + f_8 * smh0_271[k]
                   - f_9 * smh1_271[k]
                   + f_3 * pc_y[k] * smi_361[k];

        t_465[k] = f_14 * sli_250[k]
                   + f_10 * smh0_272[k]
                   - f_11 * smh1_272[k]
                   + f_3 * pc_y[k] * smi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_y, pc_y, pc_z, slk0_324, sli_223, \
                         sli_251, sli_252, slk1_324, smh0_272, smh1_272, smi_363, \
                         smi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * sli_251[k]
                   + f_3 * pc_y[k] * smi_363[k];

        t_467[k] = f_14 * sli_223[k]
                   + f_1 * smh0_272[k]
                   - f_2 * smh1_272[k]
                   + f_3 * pc_z[k] * smi_363[k];

        t_468[k] = pb_y[k] * slk0_324[k]
                   - f_12 * pc_y[k] * slk1_324[k];

        t_469[k] = f_13 * sli_252[k]
                   + f_3 * pc_y[k] * smi_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pb_y, pc_y, pc_z, slk0_327, slk0_329, \
                         sli_224, sli_253, sli_254, slk1_327, slk1_329, smi_364, \
                         smi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * sli_224[k]
                   + f_3 * pc_z[k] * smi_364[k];

        t_471[k] = pb_y[k] * slk0_327[k]
                   + f_14 * sli_253[k]
                   - f_12 * pc_y[k] * slk1_327[k];

        t_472[k] = f_13 * sli_254[k]
                   + f_3 * pc_y[k] * smi_366[k];

        t_473[k] = pb_y[k] * slk0_329[k]
                   - f_12 * pc_y[k] * slk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pc_y, pc_z, slk0_330, slk0_333, \
                         sli_227, sli_255, sli_257, slk1_330, slk1_333, smi_367, \
                         smi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pb_y[k] * slk0_330[k]
                   + f_15 * sli_255[k]
                   - f_12 * pc_y[k] * slk1_330[k];

        t_475[k] = f_15 * sli_227[k]
                   + f_3 * pc_z[k] * smi_367[k];

        t_476[k] = f_13 * sli_257[k]
                   + f_3 * pc_y[k] * smi_369[k];

        t_477[k] = pb_y[k] * slk0_333[k]
                   - f_12 * pc_y[k] * slk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pb_y, pc_y, pc_z, slk0_334, slk0_336, sli_230, \
                         sli_258, sli_260, slk1_334, slk1_336, \
                         smi_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pb_y[k] * slk0_334[k]
                   + f_16 * sli_258[k]
                   - f_12 * pc_y[k] * slk1_334[k];

        t_479[k] = f_15 * sli_230[k]
                   + f_3 * pc_z[k] * smi_370[k];

        t_480[k] = pb_y[k] * slk0_336[k]
                   + f_14 * sli_260[k]
                   - f_12 * pc_y[k] * slk1_336[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pb_y, pc_y, pc_z, slk0_338, slk0_339, \
                         sli_234, sli_261, sli_262, slk1_338, slk1_339, smi_373, \
                         smi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * sli_261[k]
                   + f_3 * pc_y[k] * smi_373[k];

        t_482[k] = pb_y[k] * slk0_338[k]
                   - f_12 * pc_y[k] * slk1_338[k];

        t_483[k] = pb_y[k] * slk0_339[k]
                   + f_17 * sli_262[k]
                   - f_12 * pc_y[k] * slk1_339[k];

        t_484[k] = f_15 * sli_234[k]
                   + f_3 * pc_z[k] * smi_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pb_y, pc_y, slk0_341, slk0_342, slk0_344, \
                         sli_264, sli_265, sli_266, slk1_341, slk1_342, slk1_344, \
                         smi_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pb_y[k] * slk0_341[k]
                   + f_15 * sli_264[k]
                   - f_12 * pc_y[k] * slk1_341[k];

        t_486[k] = pb_y[k] * slk0_342[k]
                   + f_14 * sli_265[k]
                   - f_12 * pc_y[k] * slk1_342[k];

        t_487[k] = f_13 * sli_266[k]
                   + f_3 * pc_y[k] * smi_378[k];

        t_488[k] = pb_y[k] * slk0_344[k]
                   - f_12 * pc_y[k] * slk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, sli_385, sli_386, sli_387, \
                         sli_388, sli_389, smi_385, smi_386, smi_387, smi_388, \
                         smi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_17 * sli_385[k]
                   + f_3 * pc_x[k] * smi_385[k];

        t_490[k] = f_17 * sli_386[k]
                   + f_3 * pc_x[k] * smi_386[k];

        t_491[k] = f_17 * sli_387[k]
                   + f_3 * pc_x[k] * smi_387[k];

        t_492[k] = f_17 * sli_388[k]
                   + f_3 * pc_x[k] * smi_388[k];

        t_493[k] = f_17 * sli_389[k]
                   + f_3 * pc_x[k] * smi_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, sli_245, sli_273, \
                         sli_390, sli_391, smh0_288, smh1_288, smi_385, smi_390, \
                         smi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_17 * sli_390[k]
                   + f_3 * pc_x[k] * smi_390[k];

        t_495[k] = f_17 * sli_391[k]
                   + f_3 * pc_x[k] * smi_391[k];

        t_496[k] = f_13 * sli_273[k]
                   + f_1 * smh0_288[k]
                   - f_2 * smh1_288[k]
                   + f_3 * pc_y[k] * smi_385[k];

        t_497[k] = f_15 * sli_245[k]
                   + f_3 * pc_z[k] * smi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, sli_275, sli_276, sli_277, smh0_290, \
                         smh0_291, smh0_292, smh1_290, smh1_291, smh1_292, smi_387, smi_388, \
                         smi_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * sli_275[k]
                   + f_4 * smh0_290[k]
                   - f_5 * smh1_290[k]
                   + f_3 * pc_y[k] * smi_387[k];

        t_499[k] = f_13 * sli_276[k]
                   + f_6 * smh0_291[k]
                   - f_7 * smh1_291[k]
                   + f_3 * pc_y[k] * smi_388[k];

        t_500[k] = f_13 * sli_277[k]
                   + f_8 * smh0_292[k]
                   - f_9 * smh1_292[k]
                   + f_3 * pc_y[k] * smi_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pb_y, pc_y, slk0_359, sli_278, sli_279, \
                         slk1_359, smh0_293, smh1_293, smi_390, \
                         smi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * sli_278[k]
                   + f_10 * smh0_293[k]
                   - f_11 * smh1_293[k]
                   + f_3 * pc_y[k] * smi_390[k];

        t_502[k] = f_13 * sli_279[k]
                   + f_3 * pc_y[k] * smi_391[k];

        t_503[k] = pb_y[k] * slk0_359[k]
                   - f_12 * pc_y[k] * slk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, sli_252, sli_392, \
                         sli_395, smh0_294, smh0_297, smh1_294, smh1_297, smi_392, \
                         smi_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_17 * sli_392[k]
                   + f_1 * smh0_294[k]
                   - f_2 * smh1_294[k]
                   + f_3 * pc_x[k] * smi_392[k];

        t_505[k] = f_3 * pc_y[k] * smi_392[k];

        t_506[k] = f_16 * sli_252[k]
                   + f_3 * pc_z[k] * smi_392[k];

        t_507[k] = f_17 * sli_395[k]
                   + f_4 * smh0_297[k]
                   - f_5 * smh1_297[k]
                   + f_3 * pc_x[k] * smi_395[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, pc_y, sli_397, sli_398, smh0_299, \
                         smh0_300, smh1_299, smh1_300, smi_394, smi_397, \
                         smi_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_3 * pc_y[k] * smi_394[k];

        t_509[k] = f_17 * sli_397[k]
                   + f_4 * smh0_299[k]
                   - f_5 * smh1_299[k]
                   + f_3 * pc_x[k] * smi_397[k];

        t_510[k] = f_17 * sli_398[k]
                   + f_6 * smh0_300[k]
                   - f_7 * smh1_300[k]
                   + f_3 * pc_x[k] * smi_398[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pc_x, pc_y, pc_z, sli_255, sli_401, smh0_303, \
                         smh1_303, smi_395, smi_397, smi_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_16 * sli_255[k]
                   + f_3 * pc_z[k] * smi_395[k];

        t_512[k] = f_3 * pc_y[k] * smi_397[k];

        t_513[k] = f_17 * sli_401[k]
                   + f_6 * smh0_303[k]
                   - f_7 * smh1_303[k]
                   + f_3 * pc_x[k] * smi_401[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pc_x, pc_z, sli_258, sli_402, sli_404, smh0_304, \
                         smh0_306, smh1_304, smh1_306, smi_398, smi_402, \
                         smi_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_17 * sli_402[k]
                   + f_8 * smh0_304[k]
                   - f_9 * smh1_304[k]
                   + f_3 * pc_x[k] * smi_402[k];

        t_515[k] = f_16 * sli_258[k]
                   + f_3 * pc_z[k] * smi_398[k];

        t_516[k] = f_17 * sli_404[k]
                   + f_8 * smh0_306[k]
                   - f_9 * smh1_306[k]
                   + f_3 * pc_x[k] * smi_404[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pc_x, pc_y, sli_406, sli_407, smh0_308, \
                         smh0_309, smh1_308, smh1_309, smi_401, smi_406, \
                         smi_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_3 * pc_y[k] * smi_401[k];

        t_518[k] = f_17 * sli_406[k]
                   + f_8 * smh0_308[k]
                   - f_9 * smh1_308[k]
                   + f_3 * pc_x[k] * smi_406[k];

        t_519[k] = f_17 * sli_407[k]
                   + f_10 * smh0_309[k]
                   - f_11 * smh1_309[k]
                   + f_3 * pc_x[k] * smi_407[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_x, pc_z, sli_262, sli_409, sli_410, smh0_311, \
                         smh0_312, smh1_311, smh1_312, smi_402, smi_409, \
                         smi_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_16 * sli_262[k]
                   + f_3 * pc_z[k] * smi_402[k];

        t_521[k] = f_17 * sli_409[k]
                   + f_10 * smh0_311[k]
                   - f_11 * smh1_311[k]
                   + f_3 * pc_x[k] * smi_409[k];

        t_522[k] = f_17 * sli_410[k]
                   + f_10 * smh0_312[k]
                   - f_11 * smh1_312[k]
                   + f_3 * pc_x[k] * smi_410[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pc_x, pc_y, sli_412, sli_413, sli_414, \
                         smh0_314, smh1_314, smi_406, smi_412, smi_413, \
                         smi_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_3 * pc_y[k] * smi_406[k];

        t_524[k] = f_17 * sli_412[k]
                   + f_10 * smh0_314[k]
                   - f_11 * smh1_314[k]
                   + f_3 * pc_x[k] * smi_412[k];

        t_525[k] = f_17 * sli_413[k]
                   + f_3 * pc_x[k] * smi_413[k];

        t_526[k] = f_17 * sli_414[k]
                   + f_3 * pc_x[k] * smi_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, pc_x, sli_415, sli_416, sli_417, \
                         sli_418, sli_419, smi_415, smi_416, smi_417, smi_418, \
                         smi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_17 * sli_415[k]
                   + f_3 * pc_x[k] * smi_415[k];

        t_528[k] = f_17 * sli_416[k]
                   + f_3 * pc_x[k] * smi_416[k];

        t_529[k] = f_17 * sli_417[k]
                   + f_3 * pc_x[k] * smi_417[k];

        t_530[k] = f_17 * sli_418[k]
                   + f_3 * pc_x[k] * smi_418[k];

        t_531[k] = f_17 * sli_419[k]
                   + f_3 * pc_x[k] * smi_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pc_y, pc_z, sli_273, smh0_309, smh0_311, \
                         smh0_312, smh1_309, smh1_311, smh1_312, smi_413, smi_415, \
                         smi_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_1 * smh0_309[k]
                   - f_2 * smh1_309[k]
                   + f_3 * pc_y[k] * smi_413[k];

        t_533[k] = f_16 * sli_273[k]
                   + f_3 * pc_z[k] * smi_413[k];

        t_534[k] = f_4 * smh0_311[k]
                   - f_5 * smh1_311[k]
                   + f_3 * pc_y[k] * smi_415[k];

        t_535[k] = f_6 * smh0_312[k]
                   - f_7 * smh1_312[k]
                   + f_3 * pc_y[k] * smi_416[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pc_y, pc_z, sli_279, smh0_313, smh0_314, \
                         smh1_313, smh1_314, smi_417, smi_418, \
                         smi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_8 * smh0_313[k]
                   - f_9 * smh1_313[k]
                   + f_3 * pc_y[k] * smi_417[k];

        t_537[k] = f_10 * smh0_314[k]
                   - f_11 * smh1_314[k]
                   + f_3 * pc_y[k] * smi_418[k];

        t_538[k] = f_3 * pc_y[k] * smi_419[k];

        t_539[k] = f_16 * sli_279[k]
                   + f_1 * smh0_314[k]
                   - f_2 * smh1_314[k]
                   + f_3 * pc_z[k] * smi_419[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pc_x, pc_y, pc_z, sli_280, sli_420, \
                         sli_423, smh0_315, smh0_318, smh1_315, smh1_318, smi_420, \
                         smi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_16 * sli_420[k]
                   + f_1 * smh0_315[k]
                   - f_2 * smh1_315[k]
                   + f_3 * pc_x[k] * smi_420[k];

        t_541[k] = f_17 * sli_280[k]
                   + f_3 * pc_y[k] * smi_420[k];

        t_542[k] = f_3 * pc_z[k] * smi_420[k];

        t_543[k] = f_16 * sli_423[k]
                   + f_4 * smh0_318[k]
                   - f_5 * smh1_318[k]
                   + f_3 * pc_x[k] * smi_423[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, pc_x, pc_y, sli_282, sli_425, sli_426, smh0_320, \
                         smh0_321, smh1_320, smh1_321, smi_422, smi_425, \
                         smi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_17 * sli_282[k]
                   + f_3 * pc_y[k] * smi_422[k];

        t_545[k] = f_16 * sli_425[k]
                   + f_4 * smh0_320[k]
                   - f_5 * smh1_320[k]
                   + f_3 * pc_x[k] * smi_425[k];

        t_546[k] = f_16 * sli_426[k]
                   + f_6 * smh0_321[k]
                   - f_7 * smh1_321[k]
                   + f_3 * pc_x[k] * smi_426[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pc_x, pc_y, pc_z, sli_285, sli_429, smh0_324, \
                         smh1_324, smi_423, smi_425, smi_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_3 * pc_z[k] * smi_423[k];

        t_548[k] = f_17 * sli_285[k]
                   + f_3 * pc_y[k] * smi_425[k];

        t_549[k] = f_16 * sli_429[k]
                   + f_6 * smh0_324[k]
                   - f_7 * smh1_324[k]
                   + f_3 * pc_x[k] * smi_429[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pc_x, pc_z, sli_430, sli_432, smh0_325, \
                         smh0_327, smh1_325, smh1_327, smi_426, smi_430, \
                         smi_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_16 * sli_430[k]
                   + f_8 * smh0_325[k]
                   - f_9 * smh1_325[k]
                   + f_3 * pc_x[k] * smi_430[k];

        t_551[k] = f_3 * pc_z[k] * smi_426[k];

        t_552[k] = f_16 * sli_432[k]
                   + f_8 * smh0_327[k]
                   - f_9 * smh1_327[k]
                   + f_3 * pc_x[k] * smi_432[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_x, pc_y, sli_289, sli_434, sli_435, smh0_329, \
                         smh0_330, smh1_329, smh1_330, smi_429, smi_434, \
                         smi_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_17 * sli_289[k]
                   + f_3 * pc_y[k] * smi_429[k];

        t_554[k] = f_16 * sli_434[k]
                   + f_8 * smh0_329[k]
                   - f_9 * smh1_329[k]
                   + f_3 * pc_x[k] * smi_434[k];

        t_555[k] = f_16 * sli_435[k]
                   + f_10 * smh0_330[k]
                   - f_11 * smh1_330[k]
                   + f_3 * pc_x[k] * smi_435[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_x, pc_z, sli_437, sli_438, smh0_332, \
                         smh0_333, smh1_332, smh1_333, smi_430, smi_437, \
                         smi_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_3 * pc_z[k] * smi_430[k];

        t_557[k] = f_16 * sli_437[k]
                   + f_10 * smh0_332[k]
                   - f_11 * smh1_332[k]
                   + f_3 * pc_x[k] * smi_437[k];

        t_558[k] = f_16 * sli_438[k]
                   + f_10 * smh0_333[k]
                   - f_11 * smh1_333[k]
                   + f_3 * pc_x[k] * smi_438[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pc_x, pc_y, sli_294, sli_440, sli_441, \
                         sli_442, smh0_335, smh1_335, smi_434, smi_440, smi_441, \
                         smi_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_17 * sli_294[k]
                   + f_3 * pc_y[k] * smi_434[k];

        t_560[k] = f_16 * sli_440[k]
                   + f_10 * smh0_335[k]
                   - f_11 * smh1_335[k]
                   + f_3 * pc_x[k] * smi_440[k];

        t_561[k] = f_16 * sli_441[k]
                   + f_3 * pc_x[k] * smi_441[k];

        t_562[k] = f_16 * sli_442[k]
                   + f_3 * pc_x[k] * smi_442[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, pc_x, sli_443, sli_444, sli_445, \
                         sli_446, sli_447, smi_443, smi_444, smi_445, smi_446, \
                         smi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_16 * sli_443[k]
                   + f_3 * pc_x[k] * smi_443[k];

        t_564[k] = f_16 * sli_444[k]
                   + f_3 * pc_x[k] * smi_444[k];

        t_565[k] = f_16 * sli_445[k]
                   + f_3 * pc_x[k] * smi_445[k];

        t_566[k] = f_16 * sli_446[k]
                   + f_3 * pc_x[k] * smi_446[k];

        t_567[k] = f_16 * sli_447[k]
                   + f_3 * pc_x[k] * smi_447[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pc_y, pc_z, sli_301, sli_303, smh0_330, \
                         smh0_332, smh1_330, smh1_332, smi_441, \
                         smi_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_17 * sli_301[k]
                   + f_1 * smh0_330[k]
                   - f_2 * smh1_330[k]
                   + f_3 * pc_y[k] * smi_441[k];

        t_569[k] = f_3 * pc_z[k] * smi_441[k];

        t_570[k] = f_17 * sli_303[k]
                   + f_4 * smh0_332[k]
                   - f_5 * smh1_332[k]
                   + f_3 * pc_y[k] * smi_443[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_360 = buffer.data(slk0 + 360);
    const auto *slk0_363 = buffer.data(slk0 + 363);
    const auto *slk0_366 = buffer.data(slk0 + 366);
    const auto *slk0_370 = buffer.data(slk0 + 370);
    const auto *slk0_372 = buffer.data(slk0 + 372);
    const auto *slk0_375 = buffer.data(slk0 + 375);
    const auto *slk0_377 = buffer.data(slk0 + 377);
    const auto *slk0_378 = buffer.data(slk0 + 378);
    const auto *slk0_388 = buffer.data(slk0 + 388);

    const auto *sli_280 = buffer.data(sli + 280);
    const auto *sli_283 = buffer.data(sli + 283);
    const auto *sli_286 = buffer.data(sli + 286);
    const auto *sli_287 = buffer.data(sli + 287);
    const auto *sli_290 = buffer.data(sli + 290);
    const auto *sli_291 = buffer.data(sli + 291);
    const auto *sli_292 = buffer.data(sli + 292);
    const auto *sli_301 = buffer.data(sli + 301);
    const auto *sli_304 = buffer.data(sli + 304);
    const auto *sli_305 = buffer.data(sli + 305);
    const auto *sli_306 = buffer.data(sli + 306);
    const auto *sli_307 = buffer.data(sli + 307);
    const auto *sli_308 = buffer.data(sli + 308);
    const auto *sli_310 = buffer.data(sli + 310);
    const auto *sli_311 = buffer.data(sli + 311);
    const auto *sli_313 = buffer.data(sli + 313);
    const auto *sli_314 = buffer.data(sli + 314);
    const auto *sli_317 = buffer.data(sli + 317);
    const auto *sli_318 = buffer.data(sli + 318);
    const auto *sli_322 = buffer.data(sli + 322);
    const auto *sli_329 = buffer.data(sli + 329);
    const auto *sli_331 = buffer.data(sli + 331);
    const auto *sli_332 = buffer.data(sli + 332);
    const auto *sli_333 = buffer.data(sli + 333);
    const auto *sli_334 = buffer.data(sli + 334);
    const auto *sli_335 = buffer.data(sli + 335);
    const auto *sli_336 = buffer.data(sli + 336);
    const auto *sli_338 = buffer.data(sli + 338);
    const auto *sli_339 = buffer.data(sli + 339);
    const auto *sli_341 = buffer.data(sli + 341);
    const auto *sli_342 = buffer.data(sli + 342);
    const auto *sli_345 = buffer.data(sli + 345);
    const auto *sli_346 = buffer.data(sli + 346);
    const auto *sli_350 = buffer.data(sli + 350);
    const auto *sli_357 = buffer.data(sli + 357);
    const auto *sli_359 = buffer.data(sli + 359);
    const auto *sli_360 = buffer.data(sli + 360);
    const auto *sli_361 = buffer.data(sli + 361);
    const auto *sli_362 = buffer.data(sli + 362);
    const auto *sli_363 = buffer.data(sli + 363);
    const auto *sli_364 = buffer.data(sli + 364);
    const auto *sli_366 = buffer.data(sli + 366);
    const auto *sli_369 = buffer.data(sli + 369);
    const auto *sli_373 = buffer.data(sli + 373);
    const auto *sli_378 = buffer.data(sli + 378);
    const auto *sli_385 = buffer.data(sli + 385);
    const auto *sli_387 = buffer.data(sli + 387);
    const auto *sli_453 = buffer.data(sli + 453);
    const auto *sli_457 = buffer.data(sli + 457);
    const auto *sli_462 = buffer.data(sli + 462);
    const auto *sli_468 = buffer.data(sli + 468);
    const auto *sli_469 = buffer.data(sli + 469);
    const auto *sli_470 = buffer.data(sli + 470);
    const auto *sli_471 = buffer.data(sli + 471);
    const auto *sli_472 = buffer.data(sli + 472);
    const auto *sli_473 = buffer.data(sli + 473);
    const auto *sli_474 = buffer.data(sli + 474);
    const auto *sli_475 = buffer.data(sli + 475);
    const auto *sli_476 = buffer.data(sli + 476);
    const auto *sli_479 = buffer.data(sli + 479);
    const auto *sli_481 = buffer.data(sli + 481);
    const auto *sli_482 = buffer.data(sli + 482);
    const auto *sli_485 = buffer.data(sli + 485);
    const auto *sli_486 = buffer.data(sli + 486);
    const auto *sli_488 = buffer.data(sli + 488);
    const auto *sli_490 = buffer.data(sli + 490);
    const auto *sli_491 = buffer.data(sli + 491);
    const auto *sli_493 = buffer.data(sli + 493);
    const auto *sli_494 = buffer.data(sli + 494);
    const auto *sli_496 = buffer.data(sli + 496);
    const auto *sli_497 = buffer.data(sli + 497);
    const auto *sli_498 = buffer.data(sli + 498);
    const auto *sli_499 = buffer.data(sli + 499);
    const auto *sli_500 = buffer.data(sli + 500);
    const auto *sli_501 = buffer.data(sli + 501);
    const auto *sli_502 = buffer.data(sli + 502);
    const auto *sli_503 = buffer.data(sli + 503);
    const auto *sli_504 = buffer.data(sli + 504);
    const auto *sli_507 = buffer.data(sli + 507);
    const auto *sli_509 = buffer.data(sli + 509);
    const auto *sli_510 = buffer.data(sli + 510);
    const auto *sli_513 = buffer.data(sli + 513);
    const auto *sli_514 = buffer.data(sli + 514);
    const auto *sli_516 = buffer.data(sli + 516);
    const auto *sli_518 = buffer.data(sli + 518);
    const auto *sli_519 = buffer.data(sli + 519);
    const auto *sli_521 = buffer.data(sli + 521);
    const auto *sli_522 = buffer.data(sli + 522);
    const auto *sli_524 = buffer.data(sli + 524);
    const auto *sli_525 = buffer.data(sli + 525);
    const auto *sli_526 = buffer.data(sli + 526);
    const auto *sli_527 = buffer.data(sli + 527);
    const auto *sli_528 = buffer.data(sli + 528);
    const auto *sli_529 = buffer.data(sli + 529);
    const auto *sli_530 = buffer.data(sli + 530);
    const auto *sli_531 = buffer.data(sli + 531);

    const auto *slk1_360 = buffer.data(slk1 + 360);
    const auto *slk1_363 = buffer.data(slk1 + 363);
    const auto *slk1_366 = buffer.data(slk1 + 366);
    const auto *slk1_370 = buffer.data(slk1 + 370);
    const auto *slk1_372 = buffer.data(slk1 + 372);
    const auto *slk1_375 = buffer.data(slk1 + 375);
    const auto *slk1_377 = buffer.data(slk1 + 377);
    const auto *slk1_378 = buffer.data(slk1 + 378);
    const auto *slk1_388 = buffer.data(slk1 + 388);

    const auto *smh0_333 = buffer.data(smh0 + 333);
    const auto *smh0_334 = buffer.data(smh0 + 334);
    const auto *smh0_335 = buffer.data(smh0 + 335);
    const auto *smh0_341 = buffer.data(smh0 + 341);
    const auto *smh0_345 = buffer.data(smh0 + 345);
    const auto *smh0_350 = buffer.data(smh0 + 350);
    const auto *smh0_353 = buffer.data(smh0 + 353);
    const auto *smh0_354 = buffer.data(smh0 + 354);
    const auto *smh0_355 = buffer.data(smh0 + 355);
    const auto *smh0_356 = buffer.data(smh0 + 356);
    const auto *smh0_357 = buffer.data(smh0 + 357);
    const auto *smh0_360 = buffer.data(smh0 + 360);
    const auto *smh0_362 = buffer.data(smh0 + 362);
    const auto *smh0_363 = buffer.data(smh0 + 363);
    const auto *smh0_366 = buffer.data(smh0 + 366);
    const auto *smh0_367 = buffer.data(smh0 + 367);
    const auto *smh0_369 = buffer.data(smh0 + 369);
    const auto *smh0_371 = buffer.data(smh0 + 371);
    const auto *smh0_372 = buffer.data(smh0 + 372);
    const auto *smh0_374 = buffer.data(smh0 + 374);
    const auto *smh0_375 = buffer.data(smh0 + 375);
    const auto *smh0_376 = buffer.data(smh0 + 376);
    const auto *smh0_377 = buffer.data(smh0 + 377);
    const auto *smh0_378 = buffer.data(smh0 + 378);
    const auto *smh0_381 = buffer.data(smh0 + 381);
    const auto *smh0_383 = buffer.data(smh0 + 383);
    const auto *smh0_384 = buffer.data(smh0 + 384);
    const auto *smh0_387 = buffer.data(smh0 + 387);
    const auto *smh0_388 = buffer.data(smh0 + 388);
    const auto *smh0_390 = buffer.data(smh0 + 390);
    const auto *smh0_392 = buffer.data(smh0 + 392);
    const auto *smh0_393 = buffer.data(smh0 + 393);
    const auto *smh0_395 = buffer.data(smh0 + 395);
    const auto *smh0_396 = buffer.data(smh0 + 396);
    const auto *smh0_398 = buffer.data(smh0 + 398);

    const auto *smh1_333 = buffer.data(smh1 + 333);
    const auto *smh1_334 = buffer.data(smh1 + 334);
    const auto *smh1_335 = buffer.data(smh1 + 335);
    const auto *smh1_341 = buffer.data(smh1 + 341);
    const auto *smh1_345 = buffer.data(smh1 + 345);
    const auto *smh1_350 = buffer.data(smh1 + 350);
    const auto *smh1_353 = buffer.data(smh1 + 353);
    const auto *smh1_354 = buffer.data(smh1 + 354);
    const auto *smh1_355 = buffer.data(smh1 + 355);
    const auto *smh1_356 = buffer.data(smh1 + 356);
    const auto *smh1_357 = buffer.data(smh1 + 357);
    const auto *smh1_360 = buffer.data(smh1 + 360);
    const auto *smh1_362 = buffer.data(smh1 + 362);
    const auto *smh1_363 = buffer.data(smh1 + 363);
    const auto *smh1_366 = buffer.data(smh1 + 366);
    const auto *smh1_367 = buffer.data(smh1 + 367);
    const auto *smh1_369 = buffer.data(smh1 + 369);
    const auto *smh1_371 = buffer.data(smh1 + 371);
    const auto *smh1_372 = buffer.data(smh1 + 372);
    const auto *smh1_374 = buffer.data(smh1 + 374);
    const auto *smh1_375 = buffer.data(smh1 + 375);
    const auto *smh1_376 = buffer.data(smh1 + 376);
    const auto *smh1_377 = buffer.data(smh1 + 377);
    const auto *smh1_378 = buffer.data(smh1 + 378);
    const auto *smh1_381 = buffer.data(smh1 + 381);
    const auto *smh1_383 = buffer.data(smh1 + 383);
    const auto *smh1_384 = buffer.data(smh1 + 384);
    const auto *smh1_387 = buffer.data(smh1 + 387);
    const auto *smh1_388 = buffer.data(smh1 + 388);
    const auto *smh1_390 = buffer.data(smh1 + 390);
    const auto *smh1_392 = buffer.data(smh1 + 392);
    const auto *smh1_393 = buffer.data(smh1 + 393);
    const auto *smh1_395 = buffer.data(smh1 + 395);
    const auto *smh1_396 = buffer.data(smh1 + 396);
    const auto *smh1_398 = buffer.data(smh1 + 398);

    const auto *smi_444 = buffer.data(smi + 444);
    const auto *smi_445 = buffer.data(smi + 445);
    const auto *smi_446 = buffer.data(smi + 446);
    const auto *smi_447 = buffer.data(smi + 447);
    const auto *smi_448 = buffer.data(smi + 448);
    const auto *smi_450 = buffer.data(smi + 450);
    const auto *smi_451 = buffer.data(smi + 451);
    const auto *smi_453 = buffer.data(smi + 453);
    const auto *smi_454 = buffer.data(smi + 454);
    const auto *smi_457 = buffer.data(smi + 457);
    const auto *smi_458 = buffer.data(smi + 458);
    const auto *smi_462 = buffer.data(smi + 462);
    const auto *smi_468 = buffer.data(smi + 468);
    const auto *smi_469 = buffer.data(smi + 469);
    const auto *smi_470 = buffer.data(smi + 470);
    const auto *smi_471 = buffer.data(smi + 471);
    const auto *smi_472 = buffer.data(smi + 472);
    const auto *smi_473 = buffer.data(smi + 473);
    const auto *smi_474 = buffer.data(smi + 474);
    const auto *smi_475 = buffer.data(smi + 475);
    const auto *smi_476 = buffer.data(smi + 476);
    const auto *smi_478 = buffer.data(smi + 478);
    const auto *smi_479 = buffer.data(smi + 479);
    const auto *smi_481 = buffer.data(smi + 481);
    const auto *smi_482 = buffer.data(smi + 482);
    const auto *smi_485 = buffer.data(smi + 485);
    const auto *smi_486 = buffer.data(smi + 486);
    const auto *smi_488 = buffer.data(smi + 488);
    const auto *smi_490 = buffer.data(smi + 490);
    const auto *smi_491 = buffer.data(smi + 491);
    const auto *smi_493 = buffer.data(smi + 493);
    const auto *smi_494 = buffer.data(smi + 494);
    const auto *smi_496 = buffer.data(smi + 496);
    const auto *smi_497 = buffer.data(smi + 497);
    const auto *smi_498 = buffer.data(smi + 498);
    const auto *smi_499 = buffer.data(smi + 499);
    const auto *smi_500 = buffer.data(smi + 500);
    const auto *smi_501 = buffer.data(smi + 501);
    const auto *smi_502 = buffer.data(smi + 502);
    const auto *smi_503 = buffer.data(smi + 503);
    const auto *smi_504 = buffer.data(smi + 504);
    const auto *smi_506 = buffer.data(smi + 506);
    const auto *smi_507 = buffer.data(smi + 507);
    const auto *smi_509 = buffer.data(smi + 509);
    const auto *smi_510 = buffer.data(smi + 510);
    const auto *smi_513 = buffer.data(smi + 513);
    const auto *smi_514 = buffer.data(smi + 514);
    const auto *smi_516 = buffer.data(smi + 516);
    const auto *smi_518 = buffer.data(smi + 518);
    const auto *smi_519 = buffer.data(smi + 519);
    const auto *smi_521 = buffer.data(smi + 521);
    const auto *smi_522 = buffer.data(smi + 522);
    const auto *smi_524 = buffer.data(smi + 524);
    const auto *smi_525 = buffer.data(smi + 525);
    const auto *smi_526 = buffer.data(smi + 526);
    const auto *smi_527 = buffer.data(smi + 527);
    const auto *smi_528 = buffer.data(smi + 528);
    const auto *smi_529 = buffer.data(smi + 529);
    const auto *smi_530 = buffer.data(smi + 530);
    const auto *smi_531 = buffer.data(smi + 531);

#pragma omp simd aligned(t_571, t_572, t_573, pc_y, sli_304, sli_305, sli_306, smh0_333, \
                         smh0_334, smh0_335, smh1_333, smh1_334, smh1_335, smi_444, smi_445, \
                         smi_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_17 * sli_304[k]
                   + f_6 * smh0_333[k]
                   - f_7 * smh1_333[k]
                   + f_3 * pc_y[k] * smi_444[k];

        t_572[k] = f_17 * sli_305[k]
                   + f_8 * smh0_334[k]
                   - f_9 * smh1_334[k]
                   + f_3 * pc_y[k] * smi_445[k];

        t_573[k] = f_17 * sli_306[k]
                   + f_10 * smh0_335[k]
                   - f_11 * smh1_335[k]
                   + f_3 * pc_y[k] * smi_446[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pb_z, pc_y, pc_z, slk0_360, sli_307, \
                         sli_308, slk1_360, smh0_335, smh1_335, smi_447, \
                         smi_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * sli_307[k]
                   + f_3 * pc_y[k] * smi_447[k];

        t_575[k] = f_1 * smh0_335[k]
                   - f_2 * smh1_335[k]
                   + f_3 * pc_z[k] * smi_447[k];

        t_576[k] = pb_z[k] * slk0_360[k]
                   - f_12 * pc_z[k] * slk1_360[k];

        t_577[k] = f_16 * sli_308[k]
                   + f_3 * pc_y[k] * smi_448[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pb_z, pc_y, pc_z, slk0_363, sli_280, sli_310, \
                         slk1_363, smi_448, smi_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * sli_280[k]
                   + f_3 * pc_z[k] * smi_448[k];

        t_579[k] = pb_z[k] * slk0_363[k]
                   - f_12 * pc_z[k] * slk1_363[k];

        t_580[k] = f_16 * sli_310[k]
                   + f_3 * pc_y[k] * smi_450[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pb_z, pc_x, pc_z, slk0_366, sli_283, sli_453, \
                         slk1_366, smh0_341, smh1_341, smi_451, \
                         smi_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * sli_453[k]
                   + f_4 * smh0_341[k]
                   - f_5 * smh1_341[k]
                   + f_3 * pc_x[k] * smi_453[k];

        t_582[k] = pb_z[k] * slk0_366[k]
                   - f_12 * pc_z[k] * slk1_366[k];

        t_583[k] = f_13 * sli_283[k]
                   + f_3 * pc_z[k] * smi_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pb_z, pc_x, pc_y, pc_z, slk0_370, sli_313, \
                         sli_457, slk1_370, smh0_345, smh1_345, smi_453, \
                         smi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * sli_313[k]
                   + f_3 * pc_y[k] * smi_453[k];

        t_585[k] = f_16 * sli_457[k]
                   + f_6 * smh0_345[k]
                   - f_7 * smh1_345[k]
                   + f_3 * pc_x[k] * smi_457[k];

        t_586[k] = pb_z[k] * slk0_370[k]
                   - f_12 * pc_z[k] * slk1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pb_z, pc_y, pc_z, slk0_372, sli_286, sli_287, \
                         sli_317, slk1_372, smi_454, smi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_13 * sli_286[k]
                   + f_3 * pc_z[k] * smi_454[k];

        t_588[k] = pb_z[k] * slk0_372[k]
                   + f_14 * sli_287[k]
                   - f_12 * pc_z[k] * slk1_372[k];

        t_589[k] = f_16 * sli_317[k]
                   + f_3 * pc_y[k] * smi_457[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pb_z, pc_x, pc_z, slk0_375, sli_290, sli_462, \
                         slk1_375, smh0_350, smh1_350, smi_458, \
                         smi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_16 * sli_462[k]
                   + f_8 * smh0_350[k]
                   - f_9 * smh1_350[k]
                   + f_3 * pc_x[k] * smi_462[k];

        t_591[k] = pb_z[k] * slk0_375[k]
                   - f_12 * pc_z[k] * slk1_375[k];

        t_592[k] = f_13 * sli_290[k]
                   + f_3 * pc_z[k] * smi_458[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pb_z, pc_y, pc_z, slk0_377, slk0_378, sli_291, \
                         sli_292, sli_322, slk1_377, slk1_378, \
                         smi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pb_z[k] * slk0_377[k]
                   + f_14 * sli_291[k]
                   - f_12 * pc_z[k] * slk1_377[k];

        t_594[k] = pb_z[k] * slk0_378[k]
                   + f_15 * sli_292[k]
                   - f_12 * pc_z[k] * slk1_378[k];

        t_595[k] = f_16 * sli_322[k]
                   + f_3 * pc_y[k] * smi_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, sli_468, sli_469, sli_470, sli_471, \
                         smh0_356, smh1_356, smi_468, smi_469, smi_470, \
                         smi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_16 * sli_468[k]
                   + f_10 * smh0_356[k]
                   - f_11 * smh1_356[k]
                   + f_3 * pc_x[k] * smi_468[k];

        t_597[k] = f_16 * sli_469[k]
                   + f_3 * pc_x[k] * smi_469[k];

        t_598[k] = f_16 * sli_470[k]
                   + f_3 * pc_x[k] * smi_470[k];

        t_599[k] = f_16 * sli_471[k]
                   + f_3 * pc_x[k] * smi_471[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, sli_472, sli_473, sli_474, sli_475, \
                         smi_472, smi_473, smi_474, smi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_16 * sli_472[k]
                   + f_3 * pc_x[k] * smi_472[k];

        t_601[k] = f_16 * sli_473[k]
                   + f_3 * pc_x[k] * smi_473[k];

        t_602[k] = f_16 * sli_474[k]
                   + f_3 * pc_x[k] * smi_474[k];

        t_603[k] = f_16 * sli_475[k]
                   + f_3 * pc_x[k] * smi_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pb_z, pc_y, pc_z, slk0_388, sli_301, sli_331, \
                         slk1_388, smh0_353, smh1_353, smi_469, \
                         smi_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pb_z[k] * slk0_388[k]
                   - f_12 * pc_z[k] * slk1_388[k];

        t_605[k] = f_13 * sli_301[k]
                   + f_3 * pc_z[k] * smi_469[k];

        t_606[k] = f_16 * sli_331[k]
                   + f_4 * smh0_353[k]
                   - f_5 * smh1_353[k]
                   + f_3 * pc_y[k] * smi_471[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, sli_332, sli_333, sli_334, smh0_354, \
                         smh0_355, smh0_356, smh1_354, smh1_355, smh1_356, smi_472, smi_473, \
                         smi_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_16 * sli_332[k]
                   + f_6 * smh0_354[k]
                   - f_7 * smh1_354[k]
                   + f_3 * pc_y[k] * smi_472[k];

        t_608[k] = f_16 * sli_333[k]
                   + f_8 * smh0_355[k]
                   - f_9 * smh1_355[k]
                   + f_3 * pc_y[k] * smi_473[k];

        t_609[k] = f_16 * sli_334[k]
                   + f_10 * smh0_356[k]
                   - f_11 * smh1_356[k]
                   + f_3 * pc_y[k] * smi_474[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_x, pc_y, pc_z, sli_307, sli_335, sli_476, \
                         smh0_356, smh0_357, smh1_356, smh1_357, smi_475, \
                         smi_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * sli_335[k]
                   + f_3 * pc_y[k] * smi_475[k];

        t_611[k] = f_13 * sli_307[k]
                   + f_1 * smh0_356[k]
                   - f_2 * smh1_356[k]
                   + f_3 * pc_z[k] * smi_475[k];

        t_612[k] = f_16 * sli_476[k]
                   + f_1 * smh0_357[k]
                   - f_2 * smh1_357[k]
                   + f_3 * pc_x[k] * smi_476[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pc_x, pc_y, pc_z, sli_308, sli_336, \
                         sli_338, sli_479, smh0_360, smh1_360, smi_476, smi_478, \
                         smi_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_15 * sli_336[k]
                   + f_3 * pc_y[k] * smi_476[k];

        t_614[k] = f_14 * sli_308[k]
                   + f_3 * pc_z[k] * smi_476[k];

        t_615[k] = f_16 * sli_479[k]
                   + f_4 * smh0_360[k]
                   - f_5 * smh1_360[k]
                   + f_3 * pc_x[k] * smi_479[k];

        t_616[k] = f_15 * sli_338[k]
                   + f_3 * pc_y[k] * smi_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pc_x, pc_z, sli_311, sli_481, sli_482, smh0_362, \
                         smh0_363, smh1_362, smh1_363, smi_479, smi_481, \
                         smi_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_16 * sli_481[k]
                   + f_4 * smh0_362[k]
                   - f_5 * smh1_362[k]
                   + f_3 * pc_x[k] * smi_481[k];

        t_618[k] = f_16 * sli_482[k]
                   + f_6 * smh0_363[k]
                   - f_7 * smh1_363[k]
                   + f_3 * pc_x[k] * smi_482[k];

        t_619[k] = f_14 * sli_311[k]
                   + f_3 * pc_z[k] * smi_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pc_x, pc_y, sli_341, sli_485, sli_486, smh0_366, \
                         smh0_367, smh1_366, smh1_367, smi_481, smi_485, \
                         smi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * sli_341[k]
                   + f_3 * pc_y[k] * smi_481[k];

        t_621[k] = f_16 * sli_485[k]
                   + f_6 * smh0_366[k]
                   - f_7 * smh1_366[k]
                   + f_3 * pc_x[k] * smi_485[k];

        t_622[k] = f_16 * sli_486[k]
                   + f_8 * smh0_367[k]
                   - f_9 * smh1_367[k]
                   + f_3 * pc_x[k] * smi_486[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, sli_314, sli_345, sli_488, \
                         smh0_369, smh1_369, smi_482, smi_485, \
                         smi_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * sli_314[k]
                   + f_3 * pc_z[k] * smi_482[k];

        t_624[k] = f_16 * sli_488[k]
                   + f_8 * smh0_369[k]
                   - f_9 * smh1_369[k]
                   + f_3 * pc_x[k] * smi_488[k];

        t_625[k] = f_15 * sli_345[k]
                   + f_3 * pc_y[k] * smi_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_x, pc_z, sli_318, sli_490, sli_491, smh0_371, \
                         smh0_372, smh1_371, smh1_372, smi_486, smi_490, \
                         smi_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_16 * sli_490[k]
                   + f_8 * smh0_371[k]
                   - f_9 * smh1_371[k]
                   + f_3 * pc_x[k] * smi_490[k];

        t_627[k] = f_16 * sli_491[k]
                   + f_10 * smh0_372[k]
                   - f_11 * smh1_372[k]
                   + f_3 * pc_x[k] * smi_491[k];

        t_628[k] = f_14 * sli_318[k]
                   + f_3 * pc_z[k] * smi_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, sli_350, sli_493, sli_494, smh0_374, \
                         smh0_375, smh1_374, smh1_375, smi_490, smi_493, \
                         smi_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_16 * sli_493[k]
                   + f_10 * smh0_374[k]
                   - f_11 * smh1_374[k]
                   + f_3 * pc_x[k] * smi_493[k];

        t_630[k] = f_16 * sli_494[k]
                   + f_10 * smh0_375[k]
                   - f_11 * smh1_375[k]
                   + f_3 * pc_x[k] * smi_494[k];

        t_631[k] = f_15 * sli_350[k]
                   + f_3 * pc_y[k] * smi_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, sli_496, sli_497, sli_498, sli_499, \
                         smh0_377, smh1_377, smi_496, smi_497, smi_498, \
                         smi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_16 * sli_496[k]
                   + f_10 * smh0_377[k]
                   - f_11 * smh1_377[k]
                   + f_3 * pc_x[k] * smi_496[k];

        t_633[k] = f_16 * sli_497[k]
                   + f_3 * pc_x[k] * smi_497[k];

        t_634[k] = f_16 * sli_498[k]
                   + f_3 * pc_x[k] * smi_498[k];

        t_635[k] = f_16 * sli_499[k]
                   + f_3 * pc_x[k] * smi_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, sli_500, sli_501, sli_502, sli_503, \
                         smi_500, smi_501, smi_502, smi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_16 * sli_500[k]
                   + f_3 * pc_x[k] * smi_500[k];

        t_637[k] = f_16 * sli_501[k]
                   + f_3 * pc_x[k] * smi_501[k];

        t_638[k] = f_16 * sli_502[k]
                   + f_3 * pc_x[k] * smi_502[k];

        t_639[k] = f_16 * sli_503[k]
                   + f_3 * pc_x[k] * smi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, sli_329, sli_357, sli_359, smh0_372, \
                         smh0_374, smh1_372, smh1_374, smi_497, \
                         smi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * sli_357[k]
                   + f_1 * smh0_372[k]
                   - f_2 * smh1_372[k]
                   + f_3 * pc_y[k] * smi_497[k];

        t_641[k] = f_14 * sli_329[k]
                   + f_3 * pc_z[k] * smi_497[k];

        t_642[k] = f_15 * sli_359[k]
                   + f_4 * smh0_374[k]
                   - f_5 * smh1_374[k]
                   + f_3 * pc_y[k] * smi_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, sli_360, sli_361, sli_362, smh0_375, \
                         smh0_376, smh0_377, smh1_375, smh1_376, smh1_377, smi_500, smi_501, \
                         smi_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * sli_360[k]
                   + f_6 * smh0_375[k]
                   - f_7 * smh1_375[k]
                   + f_3 * pc_y[k] * smi_500[k];

        t_644[k] = f_15 * sli_361[k]
                   + f_8 * smh0_376[k]
                   - f_9 * smh1_376[k]
                   + f_3 * pc_y[k] * smi_501[k];

        t_645[k] = f_15 * sli_362[k]
                   + f_10 * smh0_377[k]
                   - f_11 * smh1_377[k]
                   + f_3 * pc_y[k] * smi_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, sli_335, sli_363, sli_504, \
                         smh0_377, smh0_378, smh1_377, smh1_378, smi_503, \
                         smi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * sli_363[k]
                   + f_3 * pc_y[k] * smi_503[k];

        t_647[k] = f_14 * sli_335[k]
                   + f_1 * smh0_377[k]
                   - f_2 * smh1_377[k]
                   + f_3 * pc_z[k] * smi_503[k];

        t_648[k] = f_16 * sli_504[k]
                   + f_1 * smh0_378[k]
                   - f_2 * smh1_378[k]
                   + f_3 * pc_x[k] * smi_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, sli_336, sli_364, \
                         sli_366, sli_507, smh0_381, smh1_381, smi_504, smi_506, \
                         smi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_14 * sli_364[k]
                   + f_3 * pc_y[k] * smi_504[k];

        t_650[k] = f_15 * sli_336[k]
                   + f_3 * pc_z[k] * smi_504[k];

        t_651[k] = f_16 * sli_507[k]
                   + f_4 * smh0_381[k]
                   - f_5 * smh1_381[k]
                   + f_3 * pc_x[k] * smi_507[k];

        t_652[k] = f_14 * sli_366[k]
                   + f_3 * pc_y[k] * smi_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, sli_339, sli_509, sli_510, smh0_383, \
                         smh0_384, smh1_383, smh1_384, smi_507, smi_509, \
                         smi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_16 * sli_509[k]
                   + f_4 * smh0_383[k]
                   - f_5 * smh1_383[k]
                   + f_3 * pc_x[k] * smi_509[k];

        t_654[k] = f_16 * sli_510[k]
                   + f_6 * smh0_384[k]
                   - f_7 * smh1_384[k]
                   + f_3 * pc_x[k] * smi_510[k];

        t_655[k] = f_15 * sli_339[k]
                   + f_3 * pc_z[k] * smi_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, sli_369, sli_513, sli_514, smh0_387, \
                         smh0_388, smh1_387, smh1_388, smi_509, smi_513, \
                         smi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * sli_369[k]
                   + f_3 * pc_y[k] * smi_509[k];

        t_657[k] = f_16 * sli_513[k]
                   + f_6 * smh0_387[k]
                   - f_7 * smh1_387[k]
                   + f_3 * pc_x[k] * smi_513[k];

        t_658[k] = f_16 * sli_514[k]
                   + f_8 * smh0_388[k]
                   - f_9 * smh1_388[k]
                   + f_3 * pc_x[k] * smi_514[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pc_x, pc_y, pc_z, sli_342, sli_373, sli_516, \
                         smh0_390, smh1_390, smi_510, smi_513, \
                         smi_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * sli_342[k]
                   + f_3 * pc_z[k] * smi_510[k];

        t_660[k] = f_16 * sli_516[k]
                   + f_8 * smh0_390[k]
                   - f_9 * smh1_390[k]
                   + f_3 * pc_x[k] * smi_516[k];

        t_661[k] = f_14 * sli_373[k]
                   + f_3 * pc_y[k] * smi_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_z, sli_346, sli_518, sli_519, smh0_392, \
                         smh0_393, smh1_392, smh1_393, smi_514, smi_518, \
                         smi_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_16 * sli_518[k]
                   + f_8 * smh0_392[k]
                   - f_9 * smh1_392[k]
                   + f_3 * pc_x[k] * smi_518[k];

        t_663[k] = f_16 * sli_519[k]
                   + f_10 * smh0_393[k]
                   - f_11 * smh1_393[k]
                   + f_3 * pc_x[k] * smi_519[k];

        t_664[k] = f_15 * sli_346[k]
                   + f_3 * pc_z[k] * smi_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_x, pc_y, sli_378, sli_521, sli_522, smh0_395, \
                         smh0_396, smh1_395, smh1_396, smi_518, smi_521, \
                         smi_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_16 * sli_521[k]
                   + f_10 * smh0_395[k]
                   - f_11 * smh1_395[k]
                   + f_3 * pc_x[k] * smi_521[k];

        t_666[k] = f_16 * sli_522[k]
                   + f_10 * smh0_396[k]
                   - f_11 * smh1_396[k]
                   + f_3 * pc_x[k] * smi_522[k];

        t_667[k] = f_14 * sli_378[k]
                   + f_3 * pc_y[k] * smi_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pc_x, sli_524, sli_525, sli_526, sli_527, \
                         smh0_398, smh1_398, smi_524, smi_525, smi_526, \
                         smi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_16 * sli_524[k]
                   + f_10 * smh0_398[k]
                   - f_11 * smh1_398[k]
                   + f_3 * pc_x[k] * smi_524[k];

        t_669[k] = f_16 * sli_525[k]
                   + f_3 * pc_x[k] * smi_525[k];

        t_670[k] = f_16 * sli_526[k]
                   + f_3 * pc_x[k] * smi_526[k];

        t_671[k] = f_16 * sli_527[k]
                   + f_3 * pc_x[k] * smi_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, sli_528, sli_529, sli_530, sli_531, \
                         smi_528, smi_529, smi_530, smi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_16 * sli_528[k]
                   + f_3 * pc_x[k] * smi_528[k];

        t_673[k] = f_16 * sli_529[k]
                   + f_3 * pc_x[k] * smi_529[k];

        t_674[k] = f_16 * sli_530[k]
                   + f_3 * pc_x[k] * smi_530[k];

        t_675[k] = f_16 * sli_531[k]
                   + f_3 * pc_x[k] * smi_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, sli_357, sli_385, sli_387, smh0_393, \
                         smh0_395, smh1_393, smh1_395, smi_525, \
                         smi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * sli_385[k]
                   + f_1 * smh0_393[k]
                   - f_2 * smh1_393[k]
                   + f_3 * pc_y[k] * smi_525[k];

        t_677[k] = f_15 * sli_357[k]
                   + f_3 * pc_z[k] * smi_525[k];

        t_678[k] = f_14 * sli_387[k]
                   + f_4 * smh0_395[k]
                   - f_5 * smh1_395[k]
                   + f_3 * pc_y[k] * smi_527[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_20 = 3.0 / q;

    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_504 = buffer.data(slk0 + 504);
    const auto *slk0_507 = buffer.data(slk0 + 507);
    const auto *slk0_509 = buffer.data(slk0 + 509);
    const auto *slk0_510 = buffer.data(slk0 + 510);
    const auto *slk0_513 = buffer.data(slk0 + 513);
    const auto *slk0_514 = buffer.data(slk0 + 514);
    const auto *slk0_516 = buffer.data(slk0 + 516);
    const auto *slk0_518 = buffer.data(slk0 + 518);
    const auto *slk0_519 = buffer.data(slk0 + 519);
    const auto *slk0_521 = buffer.data(slk0 + 521);
    const auto *slk0_522 = buffer.data(slk0 + 522);
    const auto *slk0_524 = buffer.data(slk0 + 524);
    const auto *slk0_539 = buffer.data(slk0 + 539);

    const auto *sli_363 = buffer.data(sli + 363);
    const auto *sli_364 = buffer.data(sli + 364);
    const auto *sli_367 = buffer.data(sli + 367);
    const auto *sli_370 = buffer.data(sli + 370);
    const auto *sli_374 = buffer.data(sli + 374);
    const auto *sli_385 = buffer.data(sli + 385);
    const auto *sli_388 = buffer.data(sli + 388);
    const auto *sli_389 = buffer.data(sli + 389);
    const auto *sli_390 = buffer.data(sli + 390);
    const auto *sli_391 = buffer.data(sli + 391);
    const auto *sli_392 = buffer.data(sli + 392);
    const auto *sli_393 = buffer.data(sli + 393);
    const auto *sli_394 = buffer.data(sli + 394);
    const auto *sli_395 = buffer.data(sli + 395);
    const auto *sli_397 = buffer.data(sli + 397);
    const auto *sli_398 = buffer.data(sli + 398);
    const auto *sli_400 = buffer.data(sli + 400);
    const auto *sli_401 = buffer.data(sli + 401);
    const auto *sli_402 = buffer.data(sli + 402);
    const auto *sli_404 = buffer.data(sli + 404);
    const auto *sli_405 = buffer.data(sli + 405);
    const auto *sli_406 = buffer.data(sli + 406);
    const auto *sli_413 = buffer.data(sli + 413);
    const auto *sli_415 = buffer.data(sli + 415);
    const auto *sli_416 = buffer.data(sli + 416);
    const auto *sli_417 = buffer.data(sli + 417);
    const auto *sli_418 = buffer.data(sli + 418);
    const auto *sli_419 = buffer.data(sli + 419);
    const auto *sli_420 = buffer.data(sli + 420);
    const auto *sli_422 = buffer.data(sli + 422);
    const auto *sli_425 = buffer.data(sli + 425);
    const auto *sli_429 = buffer.data(sli + 429);
    const auto *sli_434 = buffer.data(sli + 434);
    const auto *sli_441 = buffer.data(sli + 441);
    const auto *sli_443 = buffer.data(sli + 443);
    const auto *sli_444 = buffer.data(sli + 444);
    const auto *sli_445 = buffer.data(sli + 445);
    const auto *sli_446 = buffer.data(sli + 446);
    const auto *sli_553 = buffer.data(sli + 553);
    const auto *sli_554 = buffer.data(sli + 554);
    const auto *sli_555 = buffer.data(sli + 555);
    const auto *sli_556 = buffer.data(sli + 556);
    const auto *sli_557 = buffer.data(sli + 557);
    const auto *sli_558 = buffer.data(sli + 558);
    const auto *sli_559 = buffer.data(sli + 559);
    const auto *sli_560 = buffer.data(sli + 560);
    const auto *sli_563 = buffer.data(sli + 563);
    const auto *sli_565 = buffer.data(sli + 565);
    const auto *sli_566 = buffer.data(sli + 566);
    const auto *sli_569 = buffer.data(sli + 569);
    const auto *sli_570 = buffer.data(sli + 570);
    const auto *sli_572 = buffer.data(sli + 572);
    const auto *sli_574 = buffer.data(sli + 574);
    const auto *sli_575 = buffer.data(sli + 575);
    const auto *sli_577 = buffer.data(sli + 577);
    const auto *sli_578 = buffer.data(sli + 578);
    const auto *sli_580 = buffer.data(sli + 580);
    const auto *sli_581 = buffer.data(sli + 581);
    const auto *sli_582 = buffer.data(sli + 582);
    const auto *sli_583 = buffer.data(sli + 583);
    const auto *sli_584 = buffer.data(sli + 584);
    const auto *sli_585 = buffer.data(sli + 585);
    const auto *sli_586 = buffer.data(sli + 586);
    const auto *sli_587 = buffer.data(sli + 587);
    const auto *sli_588 = buffer.data(sli + 588);
    const auto *sli_591 = buffer.data(sli + 591);
    const auto *sli_593 = buffer.data(sli + 593);
    const auto *sli_594 = buffer.data(sli + 594);
    const auto *sli_597 = buffer.data(sli + 597);
    const auto *sli_598 = buffer.data(sli + 598);
    const auto *sli_600 = buffer.data(sli + 600);
    const auto *sli_602 = buffer.data(sli + 602);
    const auto *sli_603 = buffer.data(sli + 603);
    const auto *sli_605 = buffer.data(sli + 605);
    const auto *sli_606 = buffer.data(sli + 606);
    const auto *sli_608 = buffer.data(sli + 608);
    const auto *sli_609 = buffer.data(sli + 609);
    const auto *sli_610 = buffer.data(sli + 610);
    const auto *sli_611 = buffer.data(sli + 611);
    const auto *sli_612 = buffer.data(sli + 612);
    const auto *sli_613 = buffer.data(sli + 613);
    const auto *sli_614 = buffer.data(sli + 614);
    const auto *sli_615 = buffer.data(sli + 615);

    const auto *slk1_504 = buffer.data(slk1 + 504);
    const auto *slk1_507 = buffer.data(slk1 + 507);
    const auto *slk1_509 = buffer.data(slk1 + 509);
    const auto *slk1_510 = buffer.data(slk1 + 510);
    const auto *slk1_513 = buffer.data(slk1 + 513);
    const auto *slk1_514 = buffer.data(slk1 + 514);
    const auto *slk1_516 = buffer.data(slk1 + 516);
    const auto *slk1_518 = buffer.data(slk1 + 518);
    const auto *slk1_519 = buffer.data(slk1 + 519);
    const auto *slk1_521 = buffer.data(slk1 + 521);
    const auto *slk1_522 = buffer.data(slk1 + 522);
    const auto *slk1_524 = buffer.data(slk1 + 524);
    const auto *slk1_539 = buffer.data(slk1 + 539);

    const auto *smh0_396 = buffer.data(smh0 + 396);
    const auto *smh0_397 = buffer.data(smh0 + 397);
    const auto *smh0_398 = buffer.data(smh0 + 398);
    const auto *smh0_414 = buffer.data(smh0 + 414);
    const auto *smh0_416 = buffer.data(smh0 + 416);
    const auto *smh0_417 = buffer.data(smh0 + 417);
    const auto *smh0_418 = buffer.data(smh0 + 418);
    const auto *smh0_419 = buffer.data(smh0 + 419);
    const auto *smh0_420 = buffer.data(smh0 + 420);
    const auto *smh0_423 = buffer.data(smh0 + 423);
    const auto *smh0_425 = buffer.data(smh0 + 425);
    const auto *smh0_426 = buffer.data(smh0 + 426);
    const auto *smh0_429 = buffer.data(smh0 + 429);
    const auto *smh0_430 = buffer.data(smh0 + 430);
    const auto *smh0_432 = buffer.data(smh0 + 432);
    const auto *smh0_434 = buffer.data(smh0 + 434);
    const auto *smh0_435 = buffer.data(smh0 + 435);
    const auto *smh0_437 = buffer.data(smh0 + 437);
    const auto *smh0_438 = buffer.data(smh0 + 438);
    const auto *smh0_439 = buffer.data(smh0 + 439);
    const auto *smh0_440 = buffer.data(smh0 + 440);
    const auto *smh0_441 = buffer.data(smh0 + 441);
    const auto *smh0_444 = buffer.data(smh0 + 444);
    const auto *smh0_446 = buffer.data(smh0 + 446);
    const auto *smh0_447 = buffer.data(smh0 + 447);
    const auto *smh0_450 = buffer.data(smh0 + 450);
    const auto *smh0_451 = buffer.data(smh0 + 451);
    const auto *smh0_453 = buffer.data(smh0 + 453);
    const auto *smh0_455 = buffer.data(smh0 + 455);
    const auto *smh0_456 = buffer.data(smh0 + 456);
    const auto *smh0_458 = buffer.data(smh0 + 458);
    const auto *smh0_459 = buffer.data(smh0 + 459);
    const auto *smh0_460 = buffer.data(smh0 + 460);
    const auto *smh0_461 = buffer.data(smh0 + 461);

    const auto *smh1_396 = buffer.data(smh1 + 396);
    const auto *smh1_397 = buffer.data(smh1 + 397);
    const auto *smh1_398 = buffer.data(smh1 + 398);
    const auto *smh1_414 = buffer.data(smh1 + 414);
    const auto *smh1_416 = buffer.data(smh1 + 416);
    const auto *smh1_417 = buffer.data(smh1 + 417);
    const auto *smh1_418 = buffer.data(smh1 + 418);
    const auto *smh1_419 = buffer.data(smh1 + 419);
    const auto *smh1_420 = buffer.data(smh1 + 420);
    const auto *smh1_423 = buffer.data(smh1 + 423);
    const auto *smh1_425 = buffer.data(smh1 + 425);
    const auto *smh1_426 = buffer.data(smh1 + 426);
    const auto *smh1_429 = buffer.data(smh1 + 429);
    const auto *smh1_430 = buffer.data(smh1 + 430);
    const auto *smh1_432 = buffer.data(smh1 + 432);
    const auto *smh1_434 = buffer.data(smh1 + 434);
    const auto *smh1_435 = buffer.data(smh1 + 435);
    const auto *smh1_437 = buffer.data(smh1 + 437);
    const auto *smh1_438 = buffer.data(smh1 + 438);
    const auto *smh1_439 = buffer.data(smh1 + 439);
    const auto *smh1_440 = buffer.data(smh1 + 440);
    const auto *smh1_441 = buffer.data(smh1 + 441);
    const auto *smh1_444 = buffer.data(smh1 + 444);
    const auto *smh1_446 = buffer.data(smh1 + 446);
    const auto *smh1_447 = buffer.data(smh1 + 447);
    const auto *smh1_450 = buffer.data(smh1 + 450);
    const auto *smh1_451 = buffer.data(smh1 + 451);
    const auto *smh1_453 = buffer.data(smh1 + 453);
    const auto *smh1_455 = buffer.data(smh1 + 455);
    const auto *smh1_456 = buffer.data(smh1 + 456);
    const auto *smh1_458 = buffer.data(smh1 + 458);
    const auto *smh1_459 = buffer.data(smh1 + 459);
    const auto *smh1_460 = buffer.data(smh1 + 460);
    const auto *smh1_461 = buffer.data(smh1 + 461);

    const auto *smi_528 = buffer.data(smi + 528);
    const auto *smi_529 = buffer.data(smi + 529);
    const auto *smi_530 = buffer.data(smi + 530);
    const auto *smi_531 = buffer.data(smi + 531);
    const auto *smi_532 = buffer.data(smi + 532);
    const auto *smi_534 = buffer.data(smi + 534);
    const auto *smi_535 = buffer.data(smi + 535);
    const auto *smi_537 = buffer.data(smi + 537);
    const auto *smi_538 = buffer.data(smi + 538);
    const auto *smi_541 = buffer.data(smi + 541);
    const auto *smi_542 = buffer.data(smi + 542);
    const auto *smi_546 = buffer.data(smi + 546);
    const auto *smi_553 = buffer.data(smi + 553);
    const auto *smi_554 = buffer.data(smi + 554);
    const auto *smi_555 = buffer.data(smi + 555);
    const auto *smi_556 = buffer.data(smi + 556);
    const auto *smi_557 = buffer.data(smi + 557);
    const auto *smi_558 = buffer.data(smi + 558);
    const auto *smi_559 = buffer.data(smi + 559);
    const auto *smi_560 = buffer.data(smi + 560);
    const auto *smi_562 = buffer.data(smi + 562);
    const auto *smi_563 = buffer.data(smi + 563);
    const auto *smi_565 = buffer.data(smi + 565);
    const auto *smi_566 = buffer.data(smi + 566);
    const auto *smi_569 = buffer.data(smi + 569);
    const auto *smi_570 = buffer.data(smi + 570);
    const auto *smi_572 = buffer.data(smi + 572);
    const auto *smi_574 = buffer.data(smi + 574);
    const auto *smi_575 = buffer.data(smi + 575);
    const auto *smi_577 = buffer.data(smi + 577);
    const auto *smi_578 = buffer.data(smi + 578);
    const auto *smi_580 = buffer.data(smi + 580);
    const auto *smi_581 = buffer.data(smi + 581);
    const auto *smi_582 = buffer.data(smi + 582);
    const auto *smi_583 = buffer.data(smi + 583);
    const auto *smi_584 = buffer.data(smi + 584);
    const auto *smi_585 = buffer.data(smi + 585);
    const auto *smi_586 = buffer.data(smi + 586);
    const auto *smi_587 = buffer.data(smi + 587);
    const auto *smi_588 = buffer.data(smi + 588);
    const auto *smi_590 = buffer.data(smi + 590);
    const auto *smi_591 = buffer.data(smi + 591);
    const auto *smi_593 = buffer.data(smi + 593);
    const auto *smi_594 = buffer.data(smi + 594);
    const auto *smi_597 = buffer.data(smi + 597);
    const auto *smi_598 = buffer.data(smi + 598);
    const auto *smi_600 = buffer.data(smi + 600);
    const auto *smi_602 = buffer.data(smi + 602);
    const auto *smi_603 = buffer.data(smi + 603);
    const auto *smi_605 = buffer.data(smi + 605);
    const auto *smi_606 = buffer.data(smi + 606);
    const auto *smi_608 = buffer.data(smi + 608);
    const auto *smi_609 = buffer.data(smi + 609);
    const auto *smi_610 = buffer.data(smi + 610);
    const auto *smi_611 = buffer.data(smi + 611);
    const auto *smi_612 = buffer.data(smi + 612);
    const auto *smi_613 = buffer.data(smi + 613);
    const auto *smi_614 = buffer.data(smi + 614);
    const auto *smi_615 = buffer.data(smi + 615);

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, sli_388, sli_389, sli_390, smh0_396, \
                         smh0_397, smh0_398, smh1_396, smh1_397, smh1_398, smi_528, smi_529, \
                         smi_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * sli_388[k]
                   + f_6 * smh0_396[k]
                   - f_7 * smh1_396[k]
                   + f_3 * pc_y[k] * smi_528[k];

        t_680[k] = f_14 * sli_389[k]
                   + f_8 * smh0_397[k]
                   - f_9 * smh1_397[k]
                   + f_3 * pc_y[k] * smi_529[k];

        t_681[k] = f_14 * sli_390[k]
                   + f_10 * smh0_398[k]
                   - f_11 * smh1_398[k]
                   + f_3 * pc_y[k] * smi_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pb_y, pc_y, pc_z, slk0_504, sli_363, \
                         sli_391, sli_392, slk1_504, smh0_398, smh1_398, smi_531, \
                         smi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * sli_391[k]
                   + f_3 * pc_y[k] * smi_531[k];

        t_683[k] = f_15 * sli_363[k]
                   + f_1 * smh0_398[k]
                   - f_2 * smh1_398[k]
                   + f_3 * pc_z[k] * smi_531[k];

        t_684[k] = pb_y[k] * slk0_504[k]
                   - f_12 * pc_y[k] * slk1_504[k];

        t_685[k] = f_13 * sli_392[k]
                   + f_3 * pc_y[k] * smi_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pb_y, pc_y, pc_z, slk0_507, slk0_509, \
                         sli_364, sli_393, sli_394, slk1_507, slk1_509, smi_532, \
                         smi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * sli_364[k]
                   + f_3 * pc_z[k] * smi_532[k];

        t_687[k] = pb_y[k] * slk0_507[k]
                   + f_14 * sli_393[k]
                   - f_12 * pc_y[k] * slk1_507[k];

        t_688[k] = f_13 * sli_394[k]
                   + f_3 * pc_y[k] * smi_534[k];

        t_689[k] = pb_y[k] * slk0_509[k]
                   - f_12 * pc_y[k] * slk1_509[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pb_y, pc_y, pc_z, slk0_510, slk0_513, \
                         sli_367, sli_395, sli_397, slk1_510, slk1_513, smi_535, \
                         smi_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pb_y[k] * slk0_510[k]
                   + f_15 * sli_395[k]
                   - f_12 * pc_y[k] * slk1_510[k];

        t_691[k] = f_16 * sli_367[k]
                   + f_3 * pc_z[k] * smi_535[k];

        t_692[k] = f_13 * sli_397[k]
                   + f_3 * pc_y[k] * smi_537[k];

        t_693[k] = pb_y[k] * slk0_513[k]
                   - f_12 * pc_y[k] * slk1_513[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pb_y, pc_y, pc_z, slk0_514, slk0_516, sli_370, \
                         sli_398, sli_400, slk1_514, slk1_516, \
                         smi_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pb_y[k] * slk0_514[k]
                   + f_16 * sli_398[k]
                   - f_12 * pc_y[k] * slk1_514[k];

        t_695[k] = f_16 * sli_370[k]
                   + f_3 * pc_z[k] * smi_538[k];

        t_696[k] = pb_y[k] * slk0_516[k]
                   + f_14 * sli_400[k]
                   - f_12 * pc_y[k] * slk1_516[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, pb_y, pc_y, pc_z, slk0_518, slk0_519, \
                         sli_374, sli_401, sli_402, slk1_518, slk1_519, smi_541, \
                         smi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * sli_401[k]
                   + f_3 * pc_y[k] * smi_541[k];

        t_698[k] = pb_y[k] * slk0_518[k]
                   - f_12 * pc_y[k] * slk1_518[k];

        t_699[k] = pb_y[k] * slk0_519[k]
                   + f_17 * sli_402[k]
                   - f_12 * pc_y[k] * slk1_519[k];

        t_700[k] = f_16 * sli_374[k]
                   + f_3 * pc_z[k] * smi_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pb_y, pc_y, slk0_521, slk0_522, slk0_524, \
                         sli_404, sli_405, sli_406, slk1_521, slk1_522, slk1_524, \
                         smi_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pb_y[k] * slk0_521[k]
                   + f_15 * sli_404[k]
                   - f_12 * pc_y[k] * slk1_521[k];

        t_702[k] = pb_y[k] * slk0_522[k]
                   + f_14 * sli_405[k]
                   - f_12 * pc_y[k] * slk1_522[k];

        t_703[k] = f_13 * sli_406[k]
                   + f_3 * pc_y[k] * smi_546[k];

        t_704[k] = pb_y[k] * slk0_524[k]
                   - f_12 * pc_y[k] * slk1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, sli_553, sli_554, sli_555, \
                         sli_556, sli_557, smi_553, smi_554, smi_555, smi_556, \
                         smi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_16 * sli_553[k]
                   + f_3 * pc_x[k] * smi_553[k];

        t_706[k] = f_16 * sli_554[k]
                   + f_3 * pc_x[k] * smi_554[k];

        t_707[k] = f_16 * sli_555[k]
                   + f_3 * pc_x[k] * smi_555[k];

        t_708[k] = f_16 * sli_556[k]
                   + f_3 * pc_x[k] * smi_556[k];

        t_709[k] = f_16 * sli_557[k]
                   + f_3 * pc_x[k] * smi_557[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, pc_y, pc_z, sli_385, sli_413, \
                         sli_558, sli_559, smh0_414, smh1_414, smi_553, smi_558, \
                         smi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_16 * sli_558[k]
                   + f_3 * pc_x[k] * smi_558[k];

        t_711[k] = f_16 * sli_559[k]
                   + f_3 * pc_x[k] * smi_559[k];

        t_712[k] = f_13 * sli_413[k]
                   + f_1 * smh0_414[k]
                   - f_2 * smh1_414[k]
                   + f_3 * pc_y[k] * smi_553[k];

        t_713[k] = f_16 * sli_385[k]
                   + f_3 * pc_z[k] * smi_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, sli_415, sli_416, sli_417, smh0_416, \
                         smh0_417, smh0_418, smh1_416, smh1_417, smh1_418, smi_555, smi_556, \
                         smi_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_13 * sli_415[k]
                   + f_4 * smh0_416[k]
                   - f_5 * smh1_416[k]
                   + f_3 * pc_y[k] * smi_555[k];

        t_715[k] = f_13 * sli_416[k]
                   + f_6 * smh0_417[k]
                   - f_7 * smh1_417[k]
                   + f_3 * pc_y[k] * smi_556[k];

        t_716[k] = f_13 * sli_417[k]
                   + f_8 * smh0_418[k]
                   - f_9 * smh1_418[k]
                   + f_3 * pc_y[k] * smi_557[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pb_y, pc_y, slk0_539, sli_418, sli_419, \
                         slk1_539, smh0_419, smh1_419, smi_558, \
                         smi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_13 * sli_418[k]
                   + f_10 * smh0_419[k]
                   - f_11 * smh1_419[k]
                   + f_3 * pc_y[k] * smi_558[k];

        t_718[k] = f_13 * sli_419[k]
                   + f_3 * pc_y[k] * smi_559[k];

        t_719[k] = pb_y[k] * slk0_539[k]
                   - f_12 * pc_y[k] * slk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pc_x, pc_y, pc_z, sli_392, sli_560, \
                         sli_563, smh0_420, smh0_423, smh1_420, smh1_423, smi_560, \
                         smi_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_16 * sli_560[k]
                   + f_1 * smh0_420[k]
                   - f_2 * smh1_420[k]
                   + f_3 * pc_x[k] * smi_560[k];

        t_721[k] = f_3 * pc_y[k] * smi_560[k];

        t_722[k] = f_17 * sli_392[k]
                   + f_3 * pc_z[k] * smi_560[k];

        t_723[k] = f_16 * sli_563[k]
                   + f_4 * smh0_423[k]
                   - f_5 * smh1_423[k]
                   + f_3 * pc_x[k] * smi_563[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pc_x, pc_y, sli_565, sli_566, smh0_425, \
                         smh0_426, smh1_425, smh1_426, smi_562, smi_565, \
                         smi_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_3 * pc_y[k] * smi_562[k];

        t_725[k] = f_16 * sli_565[k]
                   + f_4 * smh0_425[k]
                   - f_5 * smh1_425[k]
                   + f_3 * pc_x[k] * smi_565[k];

        t_726[k] = f_16 * sli_566[k]
                   + f_6 * smh0_426[k]
                   - f_7 * smh1_426[k]
                   + f_3 * pc_x[k] * smi_566[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_x, pc_y, pc_z, sli_395, sli_569, smh0_429, \
                         smh1_429, smi_563, smi_565, smi_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_17 * sli_395[k]
                   + f_3 * pc_z[k] * smi_563[k];

        t_728[k] = f_3 * pc_y[k] * smi_565[k];

        t_729[k] = f_16 * sli_569[k]
                   + f_6 * smh0_429[k]
                   - f_7 * smh1_429[k]
                   + f_3 * pc_x[k] * smi_569[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pc_x, pc_z, sli_398, sli_570, sli_572, smh0_430, \
                         smh0_432, smh1_430, smh1_432, smi_566, smi_570, \
                         smi_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_16 * sli_570[k]
                   + f_8 * smh0_430[k]
                   - f_9 * smh1_430[k]
                   + f_3 * pc_x[k] * smi_570[k];

        t_731[k] = f_17 * sli_398[k]
                   + f_3 * pc_z[k] * smi_566[k];

        t_732[k] = f_16 * sli_572[k]
                   + f_8 * smh0_432[k]
                   - f_9 * smh1_432[k]
                   + f_3 * pc_x[k] * smi_572[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pc_x, pc_y, sli_574, sli_575, smh0_434, \
                         smh0_435, smh1_434, smh1_435, smi_569, smi_574, \
                         smi_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_3 * pc_y[k] * smi_569[k];

        t_734[k] = f_16 * sli_574[k]
                   + f_8 * smh0_434[k]
                   - f_9 * smh1_434[k]
                   + f_3 * pc_x[k] * smi_574[k];

        t_735[k] = f_16 * sli_575[k]
                   + f_10 * smh0_435[k]
                   - f_11 * smh1_435[k]
                   + f_3 * pc_x[k] * smi_575[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pc_x, pc_z, sli_402, sli_577, sli_578, smh0_437, \
                         smh0_438, smh1_437, smh1_438, smi_570, smi_577, \
                         smi_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_17 * sli_402[k]
                   + f_3 * pc_z[k] * smi_570[k];

        t_737[k] = f_16 * sli_577[k]
                   + f_10 * smh0_437[k]
                   - f_11 * smh1_437[k]
                   + f_3 * pc_x[k] * smi_577[k];

        t_738[k] = f_16 * sli_578[k]
                   + f_10 * smh0_438[k]
                   - f_11 * smh1_438[k]
                   + f_3 * pc_x[k] * smi_578[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, pc_x, pc_y, sli_580, sli_581, sli_582, \
                         smh0_440, smh1_440, smi_574, smi_580, smi_581, \
                         smi_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * smi_574[k];

        t_740[k] = f_16 * sli_580[k]
                   + f_10 * smh0_440[k]
                   - f_11 * smh1_440[k]
                   + f_3 * pc_x[k] * smi_580[k];

        t_741[k] = f_16 * sli_581[k]
                   + f_3 * pc_x[k] * smi_581[k];

        t_742[k] = f_16 * sli_582[k]
                   + f_3 * pc_x[k] * smi_582[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, pc_x, sli_583, sli_584, sli_585, \
                         sli_586, sli_587, smi_583, smi_584, smi_585, smi_586, \
                         smi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_16 * sli_583[k]
                   + f_3 * pc_x[k] * smi_583[k];

        t_744[k] = f_16 * sli_584[k]
                   + f_3 * pc_x[k] * smi_584[k];

        t_745[k] = f_16 * sli_585[k]
                   + f_3 * pc_x[k] * smi_585[k];

        t_746[k] = f_16 * sli_586[k]
                   + f_3 * pc_x[k] * smi_586[k];

        t_747[k] = f_16 * sli_587[k]
                   + f_3 * pc_x[k] * smi_587[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_y, pc_z, sli_413, smh0_435, smh0_437, \
                         smh0_438, smh1_435, smh1_437, smh1_438, smi_581, smi_583, \
                         smi_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_1 * smh0_435[k]
                   - f_2 * smh1_435[k]
                   + f_3 * pc_y[k] * smi_581[k];

        t_749[k] = f_17 * sli_413[k]
                   + f_3 * pc_z[k] * smi_581[k];

        t_750[k] = f_4 * smh0_437[k]
                   - f_5 * smh1_437[k]
                   + f_3 * pc_y[k] * smi_583[k];

        t_751[k] = f_6 * smh0_438[k]
                   - f_7 * smh1_438[k]
                   + f_3 * pc_y[k] * smi_584[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pc_y, pc_z, sli_419, smh0_439, smh0_440, \
                         smh1_439, smh1_440, smi_585, smi_586, \
                         smi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_8 * smh0_439[k]
                   - f_9 * smh1_439[k]
                   + f_3 * pc_y[k] * smi_585[k];

        t_753[k] = f_10 * smh0_440[k]
                   - f_11 * smh1_440[k]
                   + f_3 * pc_y[k] * smi_586[k];

        t_754[k] = f_3 * pc_y[k] * smi_587[k];

        t_755[k] = f_17 * sli_419[k]
                   + f_1 * smh0_440[k]
                   - f_2 * smh1_440[k]
                   + f_3 * pc_z[k] * smi_587[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, pc_y, pc_z, sli_420, sli_588, \
                         sli_591, smh0_441, smh0_444, smh1_441, smh1_444, smi_588, \
                         smi_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_15 * sli_588[k]
                   + f_1 * smh0_441[k]
                   - f_2 * smh1_441[k]
                   + f_3 * pc_x[k] * smi_588[k];

        t_757[k] = f_20 * sli_420[k]
                   + f_3 * pc_y[k] * smi_588[k];

        t_758[k] = f_3 * pc_z[k] * smi_588[k];

        t_759[k] = f_15 * sli_591[k]
                   + f_4 * smh0_444[k]
                   - f_5 * smh1_444[k]
                   + f_3 * pc_x[k] * smi_591[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, pc_x, pc_y, sli_422, sli_593, sli_594, smh0_446, \
                         smh0_447, smh1_446, smh1_447, smi_590, smi_593, \
                         smi_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_20 * sli_422[k]
                   + f_3 * pc_y[k] * smi_590[k];

        t_761[k] = f_15 * sli_593[k]
                   + f_4 * smh0_446[k]
                   - f_5 * smh1_446[k]
                   + f_3 * pc_x[k] * smi_593[k];

        t_762[k] = f_15 * sli_594[k]
                   + f_6 * smh0_447[k]
                   - f_7 * smh1_447[k]
                   + f_3 * pc_x[k] * smi_594[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, pc_x, pc_y, pc_z, sli_425, sli_597, smh0_450, \
                         smh1_450, smi_591, smi_593, smi_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_3 * pc_z[k] * smi_591[k];

        t_764[k] = f_20 * sli_425[k]
                   + f_3 * pc_y[k] * smi_593[k];

        t_765[k] = f_15 * sli_597[k]
                   + f_6 * smh0_450[k]
                   - f_7 * smh1_450[k]
                   + f_3 * pc_x[k] * smi_597[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, pc_x, pc_z, sli_598, sli_600, smh0_451, \
                         smh0_453, smh1_451, smh1_453, smi_594, smi_598, \
                         smi_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_15 * sli_598[k]
                   + f_8 * smh0_451[k]
                   - f_9 * smh1_451[k]
                   + f_3 * pc_x[k] * smi_598[k];

        t_767[k] = f_3 * pc_z[k] * smi_594[k];

        t_768[k] = f_15 * sli_600[k]
                   + f_8 * smh0_453[k]
                   - f_9 * smh1_453[k]
                   + f_3 * pc_x[k] * smi_600[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, pc_x, pc_y, sli_429, sli_602, sli_603, smh0_455, \
                         smh0_456, smh1_455, smh1_456, smi_597, smi_602, \
                         smi_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_20 * sli_429[k]
                   + f_3 * pc_y[k] * smi_597[k];

        t_770[k] = f_15 * sli_602[k]
                   + f_8 * smh0_455[k]
                   - f_9 * smh1_455[k]
                   + f_3 * pc_x[k] * smi_602[k];

        t_771[k] = f_15 * sli_603[k]
                   + f_10 * smh0_456[k]
                   - f_11 * smh1_456[k]
                   + f_3 * pc_x[k] * smi_603[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pc_x, pc_z, sli_605, sli_606, smh0_458, \
                         smh0_459, smh1_458, smh1_459, smi_598, smi_605, \
                         smi_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * smi_598[k];

        t_773[k] = f_15 * sli_605[k]
                   + f_10 * smh0_458[k]
                   - f_11 * smh1_458[k]
                   + f_3 * pc_x[k] * smi_605[k];

        t_774[k] = f_15 * sli_606[k]
                   + f_10 * smh0_459[k]
                   - f_11 * smh1_459[k]
                   + f_3 * pc_x[k] * smi_606[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, pc_x, pc_y, sli_434, sli_608, sli_609, \
                         sli_610, smh0_461, smh1_461, smi_602, smi_608, smi_609, \
                         smi_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_20 * sli_434[k]
                   + f_3 * pc_y[k] * smi_602[k];

        t_776[k] = f_15 * sli_608[k]
                   + f_10 * smh0_461[k]
                   - f_11 * smh1_461[k]
                   + f_3 * pc_x[k] * smi_608[k];

        t_777[k] = f_15 * sli_609[k]
                   + f_3 * pc_x[k] * smi_609[k];

        t_778[k] = f_15 * sli_610[k]
                   + f_3 * pc_x[k] * smi_610[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, t_783, pc_x, sli_611, sli_612, sli_613, \
                         sli_614, sli_615, smi_611, smi_612, smi_613, smi_614, \
                         smi_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_15 * sli_611[k]
                   + f_3 * pc_x[k] * smi_611[k];

        t_780[k] = f_15 * sli_612[k]
                   + f_3 * pc_x[k] * smi_612[k];

        t_781[k] = f_15 * sli_613[k]
                   + f_3 * pc_x[k] * smi_613[k];

        t_782[k] = f_15 * sli_614[k]
                   + f_3 * pc_x[k] * smi_614[k];

        t_783[k] = f_15 * sli_615[k]
                   + f_3 * pc_x[k] * smi_615[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, pc_y, pc_z, sli_441, sli_443, smh0_456, \
                         smh0_458, smh1_456, smh1_458, smi_609, \
                         smi_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_20 * sli_441[k]
                   + f_1 * smh0_456[k]
                   - f_2 * smh1_456[k]
                   + f_3 * pc_y[k] * smi_609[k];

        t_785[k] = f_3 * pc_z[k] * smi_609[k];

        t_786[k] = f_20 * sli_443[k]
                   + f_4 * smh0_458[k]
                   - f_5 * smh1_458[k]
                   + f_3 * pc_y[k] * smi_611[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_y, sli_444, sli_445, sli_446, smh0_459, \
                         smh0_460, smh0_461, smh1_459, smh1_460, smh1_461, smi_612, smi_613, \
                         smi_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_20 * sli_444[k]
                   + f_6 * smh0_459[k]
                   - f_7 * smh1_459[k]
                   + f_3 * pc_y[k] * smi_612[k];

        t_788[k] = f_20 * sli_445[k]
                   + f_8 * smh0_460[k]
                   - f_9 * smh1_460[k]
                   + f_3 * pc_y[k] * smi_613[k];

        t_789[k] = f_20 * sli_446[k]
                   + f_10 * smh0_461[k]
                   - f_11 * smh1_461[k]
                   + f_3 * pc_y[k] * smi_614[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_20 = 3.0 / q;

    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_540 = buffer.data(slk0 + 540);
    const auto *slk0_543 = buffer.data(slk0 + 543);
    const auto *slk0_546 = buffer.data(slk0 + 546);
    const auto *slk0_550 = buffer.data(slk0 + 550);
    const auto *slk0_552 = buffer.data(slk0 + 552);
    const auto *slk0_555 = buffer.data(slk0 + 555);
    const auto *slk0_557 = buffer.data(slk0 + 557);
    const auto *slk0_558 = buffer.data(slk0 + 558);
    const auto *slk0_568 = buffer.data(slk0 + 568);

    const auto *sli_420 = buffer.data(sli + 420);
    const auto *sli_423 = buffer.data(sli + 423);
    const auto *sli_426 = buffer.data(sli + 426);
    const auto *sli_427 = buffer.data(sli + 427);
    const auto *sli_430 = buffer.data(sli + 430);
    const auto *sli_431 = buffer.data(sli + 431);
    const auto *sli_432 = buffer.data(sli + 432);
    const auto *sli_441 = buffer.data(sli + 441);
    const auto *sli_447 = buffer.data(sli + 447);
    const auto *sli_448 = buffer.data(sli + 448);
    const auto *sli_450 = buffer.data(sli + 450);
    const auto *sli_451 = buffer.data(sli + 451);
    const auto *sli_453 = buffer.data(sli + 453);
    const auto *sli_454 = buffer.data(sli + 454);
    const auto *sli_457 = buffer.data(sli + 457);
    const auto *sli_458 = buffer.data(sli + 458);
    const auto *sli_462 = buffer.data(sli + 462);
    const auto *sli_469 = buffer.data(sli + 469);
    const auto *sli_471 = buffer.data(sli + 471);
    const auto *sli_472 = buffer.data(sli + 472);
    const auto *sli_473 = buffer.data(sli + 473);
    const auto *sli_474 = buffer.data(sli + 474);
    const auto *sli_475 = buffer.data(sli + 475);
    const auto *sli_476 = buffer.data(sli + 476);
    const auto *sli_478 = buffer.data(sli + 478);
    const auto *sli_479 = buffer.data(sli + 479);
    const auto *sli_481 = buffer.data(sli + 481);
    const auto *sli_482 = buffer.data(sli + 482);
    const auto *sli_485 = buffer.data(sli + 485);
    const auto *sli_486 = buffer.data(sli + 486);
    const auto *sli_490 = buffer.data(sli + 490);
    const auto *sli_497 = buffer.data(sli + 497);
    const auto *sli_499 = buffer.data(sli + 499);
    const auto *sli_500 = buffer.data(sli + 500);
    const auto *sli_501 = buffer.data(sli + 501);
    const auto *sli_502 = buffer.data(sli + 502);
    const auto *sli_503 = buffer.data(sli + 503);
    const auto *sli_504 = buffer.data(sli + 504);
    const auto *sli_506 = buffer.data(sli + 506);
    const auto *sli_509 = buffer.data(sli + 509);
    const auto *sli_513 = buffer.data(sli + 513);
    const auto *sli_518 = buffer.data(sli + 518);
    const auto *sli_525 = buffer.data(sli + 525);
    const auto *sli_527 = buffer.data(sli + 527);
    const auto *sli_528 = buffer.data(sli + 528);
    const auto *sli_529 = buffer.data(sli + 529);
    const auto *sli_530 = buffer.data(sli + 530);
    const auto *sli_621 = buffer.data(sli + 621);
    const auto *sli_625 = buffer.data(sli + 625);
    const auto *sli_630 = buffer.data(sli + 630);
    const auto *sli_636 = buffer.data(sli + 636);
    const auto *sli_637 = buffer.data(sli + 637);
    const auto *sli_638 = buffer.data(sli + 638);
    const auto *sli_639 = buffer.data(sli + 639);
    const auto *sli_640 = buffer.data(sli + 640);
    const auto *sli_641 = buffer.data(sli + 641);
    const auto *sli_642 = buffer.data(sli + 642);
    const auto *sli_643 = buffer.data(sli + 643);
    const auto *sli_644 = buffer.data(sli + 644);
    const auto *sli_647 = buffer.data(sli + 647);
    const auto *sli_649 = buffer.data(sli + 649);
    const auto *sli_650 = buffer.data(sli + 650);
    const auto *sli_653 = buffer.data(sli + 653);
    const auto *sli_654 = buffer.data(sli + 654);
    const auto *sli_656 = buffer.data(sli + 656);
    const auto *sli_658 = buffer.data(sli + 658);
    const auto *sli_659 = buffer.data(sli + 659);
    const auto *sli_661 = buffer.data(sli + 661);
    const auto *sli_662 = buffer.data(sli + 662);
    const auto *sli_664 = buffer.data(sli + 664);
    const auto *sli_665 = buffer.data(sli + 665);
    const auto *sli_666 = buffer.data(sli + 666);
    const auto *sli_667 = buffer.data(sli + 667);
    const auto *sli_668 = buffer.data(sli + 668);
    const auto *sli_669 = buffer.data(sli + 669);
    const auto *sli_670 = buffer.data(sli + 670);
    const auto *sli_671 = buffer.data(sli + 671);
    const auto *sli_672 = buffer.data(sli + 672);
    const auto *sli_675 = buffer.data(sli + 675);
    const auto *sli_677 = buffer.data(sli + 677);
    const auto *sli_678 = buffer.data(sli + 678);
    const auto *sli_681 = buffer.data(sli + 681);
    const auto *sli_682 = buffer.data(sli + 682);
    const auto *sli_684 = buffer.data(sli + 684);
    const auto *sli_686 = buffer.data(sli + 686);
    const auto *sli_687 = buffer.data(sli + 687);
    const auto *sli_689 = buffer.data(sli + 689);
    const auto *sli_690 = buffer.data(sli + 690);
    const auto *sli_692 = buffer.data(sli + 692);
    const auto *sli_693 = buffer.data(sli + 693);
    const auto *sli_694 = buffer.data(sli + 694);
    const auto *sli_695 = buffer.data(sli + 695);
    const auto *sli_696 = buffer.data(sli + 696);
    const auto *sli_697 = buffer.data(sli + 697);
    const auto *sli_698 = buffer.data(sli + 698);
    const auto *sli_699 = buffer.data(sli + 699);

    const auto *slk1_540 = buffer.data(slk1 + 540);
    const auto *slk1_543 = buffer.data(slk1 + 543);
    const auto *slk1_546 = buffer.data(slk1 + 546);
    const auto *slk1_550 = buffer.data(slk1 + 550);
    const auto *slk1_552 = buffer.data(slk1 + 552);
    const auto *slk1_555 = buffer.data(slk1 + 555);
    const auto *slk1_557 = buffer.data(slk1 + 557);
    const auto *slk1_558 = buffer.data(slk1 + 558);
    const auto *slk1_568 = buffer.data(slk1 + 568);

    const auto *smh0_461 = buffer.data(smh0 + 461);
    const auto *smh0_467 = buffer.data(smh0 + 467);
    const auto *smh0_471 = buffer.data(smh0 + 471);
    const auto *smh0_476 = buffer.data(smh0 + 476);
    const auto *smh0_479 = buffer.data(smh0 + 479);
    const auto *smh0_480 = buffer.data(smh0 + 480);
    const auto *smh0_481 = buffer.data(smh0 + 481);
    const auto *smh0_482 = buffer.data(smh0 + 482);
    const auto *smh0_483 = buffer.data(smh0 + 483);
    const auto *smh0_486 = buffer.data(smh0 + 486);
    const auto *smh0_488 = buffer.data(smh0 + 488);
    const auto *smh0_489 = buffer.data(smh0 + 489);
    const auto *smh0_492 = buffer.data(smh0 + 492);
    const auto *smh0_493 = buffer.data(smh0 + 493);
    const auto *smh0_495 = buffer.data(smh0 + 495);
    const auto *smh0_497 = buffer.data(smh0 + 497);
    const auto *smh0_498 = buffer.data(smh0 + 498);
    const auto *smh0_500 = buffer.data(smh0 + 500);
    const auto *smh0_501 = buffer.data(smh0 + 501);
    const auto *smh0_502 = buffer.data(smh0 + 502);
    const auto *smh0_503 = buffer.data(smh0 + 503);
    const auto *smh0_504 = buffer.data(smh0 + 504);
    const auto *smh0_507 = buffer.data(smh0 + 507);
    const auto *smh0_509 = buffer.data(smh0 + 509);
    const auto *smh0_510 = buffer.data(smh0 + 510);
    const auto *smh0_513 = buffer.data(smh0 + 513);
    const auto *smh0_514 = buffer.data(smh0 + 514);
    const auto *smh0_516 = buffer.data(smh0 + 516);
    const auto *smh0_518 = buffer.data(smh0 + 518);
    const auto *smh0_519 = buffer.data(smh0 + 519);
    const auto *smh0_521 = buffer.data(smh0 + 521);
    const auto *smh0_522 = buffer.data(smh0 + 522);
    const auto *smh0_523 = buffer.data(smh0 + 523);
    const auto *smh0_524 = buffer.data(smh0 + 524);

    const auto *smh1_461 = buffer.data(smh1 + 461);
    const auto *smh1_467 = buffer.data(smh1 + 467);
    const auto *smh1_471 = buffer.data(smh1 + 471);
    const auto *smh1_476 = buffer.data(smh1 + 476);
    const auto *smh1_479 = buffer.data(smh1 + 479);
    const auto *smh1_480 = buffer.data(smh1 + 480);
    const auto *smh1_481 = buffer.data(smh1 + 481);
    const auto *smh1_482 = buffer.data(smh1 + 482);
    const auto *smh1_483 = buffer.data(smh1 + 483);
    const auto *smh1_486 = buffer.data(smh1 + 486);
    const auto *smh1_488 = buffer.data(smh1 + 488);
    const auto *smh1_489 = buffer.data(smh1 + 489);
    const auto *smh1_492 = buffer.data(smh1 + 492);
    const auto *smh1_493 = buffer.data(smh1 + 493);
    const auto *smh1_495 = buffer.data(smh1 + 495);
    const auto *smh1_497 = buffer.data(smh1 + 497);
    const auto *smh1_498 = buffer.data(smh1 + 498);
    const auto *smh1_500 = buffer.data(smh1 + 500);
    const auto *smh1_501 = buffer.data(smh1 + 501);
    const auto *smh1_502 = buffer.data(smh1 + 502);
    const auto *smh1_503 = buffer.data(smh1 + 503);
    const auto *smh1_504 = buffer.data(smh1 + 504);
    const auto *smh1_507 = buffer.data(smh1 + 507);
    const auto *smh1_509 = buffer.data(smh1 + 509);
    const auto *smh1_510 = buffer.data(smh1 + 510);
    const auto *smh1_513 = buffer.data(smh1 + 513);
    const auto *smh1_514 = buffer.data(smh1 + 514);
    const auto *smh1_516 = buffer.data(smh1 + 516);
    const auto *smh1_518 = buffer.data(smh1 + 518);
    const auto *smh1_519 = buffer.data(smh1 + 519);
    const auto *smh1_521 = buffer.data(smh1 + 521);
    const auto *smh1_522 = buffer.data(smh1 + 522);
    const auto *smh1_523 = buffer.data(smh1 + 523);
    const auto *smh1_524 = buffer.data(smh1 + 524);

    const auto *smi_615 = buffer.data(smi + 615);
    const auto *smi_616 = buffer.data(smi + 616);
    const auto *smi_618 = buffer.data(smi + 618);
    const auto *smi_619 = buffer.data(smi + 619);
    const auto *smi_621 = buffer.data(smi + 621);
    const auto *smi_622 = buffer.data(smi + 622);
    const auto *smi_625 = buffer.data(smi + 625);
    const auto *smi_626 = buffer.data(smi + 626);
    const auto *smi_630 = buffer.data(smi + 630);
    const auto *smi_636 = buffer.data(smi + 636);
    const auto *smi_637 = buffer.data(smi + 637);
    const auto *smi_638 = buffer.data(smi + 638);
    const auto *smi_639 = buffer.data(smi + 639);
    const auto *smi_640 = buffer.data(smi + 640);
    const auto *smi_641 = buffer.data(smi + 641);
    const auto *smi_642 = buffer.data(smi + 642);
    const auto *smi_643 = buffer.data(smi + 643);
    const auto *smi_644 = buffer.data(smi + 644);
    const auto *smi_646 = buffer.data(smi + 646);
    const auto *smi_647 = buffer.data(smi + 647);
    const auto *smi_649 = buffer.data(smi + 649);
    const auto *smi_650 = buffer.data(smi + 650);
    const auto *smi_653 = buffer.data(smi + 653);
    const auto *smi_654 = buffer.data(smi + 654);
    const auto *smi_656 = buffer.data(smi + 656);
    const auto *smi_658 = buffer.data(smi + 658);
    const auto *smi_659 = buffer.data(smi + 659);
    const auto *smi_661 = buffer.data(smi + 661);
    const auto *smi_662 = buffer.data(smi + 662);
    const auto *smi_664 = buffer.data(smi + 664);
    const auto *smi_665 = buffer.data(smi + 665);
    const auto *smi_666 = buffer.data(smi + 666);
    const auto *smi_667 = buffer.data(smi + 667);
    const auto *smi_668 = buffer.data(smi + 668);
    const auto *smi_669 = buffer.data(smi + 669);
    const auto *smi_670 = buffer.data(smi + 670);
    const auto *smi_671 = buffer.data(smi + 671);
    const auto *smi_672 = buffer.data(smi + 672);
    const auto *smi_674 = buffer.data(smi + 674);
    const auto *smi_675 = buffer.data(smi + 675);
    const auto *smi_677 = buffer.data(smi + 677);
    const auto *smi_678 = buffer.data(smi + 678);
    const auto *smi_681 = buffer.data(smi + 681);
    const auto *smi_682 = buffer.data(smi + 682);
    const auto *smi_684 = buffer.data(smi + 684);
    const auto *smi_686 = buffer.data(smi + 686);
    const auto *smi_687 = buffer.data(smi + 687);
    const auto *smi_689 = buffer.data(smi + 689);
    const auto *smi_690 = buffer.data(smi + 690);
    const auto *smi_692 = buffer.data(smi + 692);
    const auto *smi_693 = buffer.data(smi + 693);
    const auto *smi_694 = buffer.data(smi + 694);
    const auto *smi_695 = buffer.data(smi + 695);
    const auto *smi_696 = buffer.data(smi + 696);
    const auto *smi_697 = buffer.data(smi + 697);
    const auto *smi_698 = buffer.data(smi + 698);
    const auto *smi_699 = buffer.data(smi + 699);

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pb_z, pc_y, pc_z, slk0_540, sli_447, \
                         sli_448, slk1_540, smh0_461, smh1_461, smi_615, \
                         smi_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_20 * sli_447[k]
                   + f_3 * pc_y[k] * smi_615[k];

        t_791[k] = f_1 * smh0_461[k]
                   - f_2 * smh1_461[k]
                   + f_3 * pc_z[k] * smi_615[k];

        t_792[k] = pb_z[k] * slk0_540[k]
                   - f_12 * pc_z[k] * slk1_540[k];

        t_793[k] = f_17 * sli_448[k]
                   + f_3 * pc_y[k] * smi_616[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pb_z, pc_y, pc_z, slk0_543, sli_420, sli_450, \
                         slk1_543, smi_616, smi_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_13 * sli_420[k]
                   + f_3 * pc_z[k] * smi_616[k];

        t_795[k] = pb_z[k] * slk0_543[k]
                   - f_12 * pc_z[k] * slk1_543[k];

        t_796[k] = f_17 * sli_450[k]
                   + f_3 * pc_y[k] * smi_618[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pb_z, pc_x, pc_z, slk0_546, sli_423, sli_621, \
                         slk1_546, smh0_467, smh1_467, smi_619, \
                         smi_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_15 * sli_621[k]
                   + f_4 * smh0_467[k]
                   - f_5 * smh1_467[k]
                   + f_3 * pc_x[k] * smi_621[k];

        t_798[k] = pb_z[k] * slk0_546[k]
                   - f_12 * pc_z[k] * slk1_546[k];

        t_799[k] = f_13 * sli_423[k]
                   + f_3 * pc_z[k] * smi_619[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pb_z, pc_x, pc_y, pc_z, slk0_550, sli_453, \
                         sli_625, slk1_550, smh0_471, smh1_471, smi_621, \
                         smi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_17 * sli_453[k]
                   + f_3 * pc_y[k] * smi_621[k];

        t_801[k] = f_15 * sli_625[k]
                   + f_6 * smh0_471[k]
                   - f_7 * smh1_471[k]
                   + f_3 * pc_x[k] * smi_625[k];

        t_802[k] = pb_z[k] * slk0_550[k]
                   - f_12 * pc_z[k] * slk1_550[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pb_z, pc_y, pc_z, slk0_552, sli_426, sli_427, \
                         sli_457, slk1_552, smi_622, smi_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_13 * sli_426[k]
                   + f_3 * pc_z[k] * smi_622[k];

        t_804[k] = pb_z[k] * slk0_552[k]
                   + f_14 * sli_427[k]
                   - f_12 * pc_z[k] * slk1_552[k];

        t_805[k] = f_17 * sli_457[k]
                   + f_3 * pc_y[k] * smi_625[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pb_z, pc_x, pc_z, slk0_555, sli_430, sli_630, \
                         slk1_555, smh0_476, smh1_476, smi_626, \
                         smi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_15 * sli_630[k]
                   + f_8 * smh0_476[k]
                   - f_9 * smh1_476[k]
                   + f_3 * pc_x[k] * smi_630[k];

        t_807[k] = pb_z[k] * slk0_555[k]
                   - f_12 * pc_z[k] * slk1_555[k];

        t_808[k] = f_13 * sli_430[k]
                   + f_3 * pc_z[k] * smi_626[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pb_z, pc_y, pc_z, slk0_557, slk0_558, sli_431, \
                         sli_432, sli_462, slk1_557, slk1_558, \
                         smi_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pb_z[k] * slk0_557[k]
                   + f_14 * sli_431[k]
                   - f_12 * pc_z[k] * slk1_557[k];

        t_810[k] = pb_z[k] * slk0_558[k]
                   + f_15 * sli_432[k]
                   - f_12 * pc_z[k] * slk1_558[k];

        t_811[k] = f_17 * sli_462[k]
                   + f_3 * pc_y[k] * smi_630[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pc_x, sli_636, sli_637, sli_638, sli_639, \
                         smh0_482, smh1_482, smi_636, smi_637, smi_638, \
                         smi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_15 * sli_636[k]
                   + f_10 * smh0_482[k]
                   - f_11 * smh1_482[k]
                   + f_3 * pc_x[k] * smi_636[k];

        t_813[k] = f_15 * sli_637[k]
                   + f_3 * pc_x[k] * smi_637[k];

        t_814[k] = f_15 * sli_638[k]
                   + f_3 * pc_x[k] * smi_638[k];

        t_815[k] = f_15 * sli_639[k]
                   + f_3 * pc_x[k] * smi_639[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, pc_x, sli_640, sli_641, sli_642, sli_643, \
                         smi_640, smi_641, smi_642, smi_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_15 * sli_640[k]
                   + f_3 * pc_x[k] * smi_640[k];

        t_817[k] = f_15 * sli_641[k]
                   + f_3 * pc_x[k] * smi_641[k];

        t_818[k] = f_15 * sli_642[k]
                   + f_3 * pc_x[k] * smi_642[k];

        t_819[k] = f_15 * sli_643[k]
                   + f_3 * pc_x[k] * smi_643[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, pb_z, pc_y, pc_z, slk0_568, sli_441, sli_471, \
                         slk1_568, smh0_479, smh1_479, smi_637, \
                         smi_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = pb_z[k] * slk0_568[k]
                   - f_12 * pc_z[k] * slk1_568[k];

        t_821[k] = f_13 * sli_441[k]
                   + f_3 * pc_z[k] * smi_637[k];

        t_822[k] = f_17 * sli_471[k]
                   + f_4 * smh0_479[k]
                   - f_5 * smh1_479[k]
                   + f_3 * pc_y[k] * smi_639[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pc_y, sli_472, sli_473, sli_474, smh0_480, \
                         smh0_481, smh0_482, smh1_480, smh1_481, smh1_482, smi_640, smi_641, \
                         smi_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_17 * sli_472[k]
                   + f_6 * smh0_480[k]
                   - f_7 * smh1_480[k]
                   + f_3 * pc_y[k] * smi_640[k];

        t_824[k] = f_17 * sli_473[k]
                   + f_8 * smh0_481[k]
                   - f_9 * smh1_481[k]
                   + f_3 * pc_y[k] * smi_641[k];

        t_825[k] = f_17 * sli_474[k]
                   + f_10 * smh0_482[k]
                   - f_11 * smh1_482[k]
                   + f_3 * pc_y[k] * smi_642[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, pc_x, pc_y, pc_z, sli_447, sli_475, sli_644, \
                         smh0_482, smh0_483, smh1_482, smh1_483, smi_643, \
                         smi_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_17 * sli_475[k]
                   + f_3 * pc_y[k] * smi_643[k];

        t_827[k] = f_13 * sli_447[k]
                   + f_1 * smh0_482[k]
                   - f_2 * smh1_482[k]
                   + f_3 * pc_z[k] * smi_643[k];

        t_828[k] = f_15 * sli_644[k]
                   + f_1 * smh0_483[k]
                   - f_2 * smh1_483[k]
                   + f_3 * pc_x[k] * smi_644[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pc_x, pc_y, pc_z, sli_448, sli_476, \
                         sli_478, sli_647, smh0_486, smh1_486, smi_644, smi_646, \
                         smi_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_16 * sli_476[k]
                   + f_3 * pc_y[k] * smi_644[k];

        t_830[k] = f_14 * sli_448[k]
                   + f_3 * pc_z[k] * smi_644[k];

        t_831[k] = f_15 * sli_647[k]
                   + f_4 * smh0_486[k]
                   - f_5 * smh1_486[k]
                   + f_3 * pc_x[k] * smi_647[k];

        t_832[k] = f_16 * sli_478[k]
                   + f_3 * pc_y[k] * smi_646[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, pc_z, sli_451, sli_649, sli_650, smh0_488, \
                         smh0_489, smh1_488, smh1_489, smi_647, smi_649, \
                         smi_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_15 * sli_649[k]
                   + f_4 * smh0_488[k]
                   - f_5 * smh1_488[k]
                   + f_3 * pc_x[k] * smi_649[k];

        t_834[k] = f_15 * sli_650[k]
                   + f_6 * smh0_489[k]
                   - f_7 * smh1_489[k]
                   + f_3 * pc_x[k] * smi_650[k];

        t_835[k] = f_14 * sli_451[k]
                   + f_3 * pc_z[k] * smi_647[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pc_x, pc_y, sli_481, sli_653, sli_654, smh0_492, \
                         smh0_493, smh1_492, smh1_493, smi_649, smi_653, \
                         smi_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * sli_481[k]
                   + f_3 * pc_y[k] * smi_649[k];

        t_837[k] = f_15 * sli_653[k]
                   + f_6 * smh0_492[k]
                   - f_7 * smh1_492[k]
                   + f_3 * pc_x[k] * smi_653[k];

        t_838[k] = f_15 * sli_654[k]
                   + f_8 * smh0_493[k]
                   - f_9 * smh1_493[k]
                   + f_3 * pc_x[k] * smi_654[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pc_x, pc_y, pc_z, sli_454, sli_485, sli_656, \
                         smh0_495, smh1_495, smi_650, smi_653, \
                         smi_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_14 * sli_454[k]
                   + f_3 * pc_z[k] * smi_650[k];

        t_840[k] = f_15 * sli_656[k]
                   + f_8 * smh0_495[k]
                   - f_9 * smh1_495[k]
                   + f_3 * pc_x[k] * smi_656[k];

        t_841[k] = f_16 * sli_485[k]
                   + f_3 * pc_y[k] * smi_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, pc_z, sli_458, sli_658, sli_659, smh0_497, \
                         smh0_498, smh1_497, smh1_498, smi_654, smi_658, \
                         smi_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_15 * sli_658[k]
                   + f_8 * smh0_497[k]
                   - f_9 * smh1_497[k]
                   + f_3 * pc_x[k] * smi_658[k];

        t_843[k] = f_15 * sli_659[k]
                   + f_10 * smh0_498[k]
                   - f_11 * smh1_498[k]
                   + f_3 * pc_x[k] * smi_659[k];

        t_844[k] = f_14 * sli_458[k]
                   + f_3 * pc_z[k] * smi_654[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, sli_490, sli_661, sli_662, smh0_500, \
                         smh0_501, smh1_500, smh1_501, smi_658, smi_661, \
                         smi_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_15 * sli_661[k]
                   + f_10 * smh0_500[k]
                   - f_11 * smh1_500[k]
                   + f_3 * pc_x[k] * smi_661[k];

        t_846[k] = f_15 * sli_662[k]
                   + f_10 * smh0_501[k]
                   - f_11 * smh1_501[k]
                   + f_3 * pc_x[k] * smi_662[k];

        t_847[k] = f_16 * sli_490[k]
                   + f_3 * pc_y[k] * smi_658[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, sli_664, sli_665, sli_666, sli_667, \
                         smh0_503, smh1_503, smi_664, smi_665, smi_666, \
                         smi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_15 * sli_664[k]
                   + f_10 * smh0_503[k]
                   - f_11 * smh1_503[k]
                   + f_3 * pc_x[k] * smi_664[k];

        t_849[k] = f_15 * sli_665[k]
                   + f_3 * pc_x[k] * smi_665[k];

        t_850[k] = f_15 * sli_666[k]
                   + f_3 * pc_x[k] * smi_666[k];

        t_851[k] = f_15 * sli_667[k]
                   + f_3 * pc_x[k] * smi_667[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, sli_668, sli_669, sli_670, sli_671, \
                         smi_668, smi_669, smi_670, smi_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_15 * sli_668[k]
                   + f_3 * pc_x[k] * smi_668[k];

        t_853[k] = f_15 * sli_669[k]
                   + f_3 * pc_x[k] * smi_669[k];

        t_854[k] = f_15 * sli_670[k]
                   + f_3 * pc_x[k] * smi_670[k];

        t_855[k] = f_15 * sli_671[k]
                   + f_3 * pc_x[k] * smi_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, sli_469, sli_497, sli_499, smh0_498, \
                         smh0_500, smh1_498, smh1_500, smi_665, \
                         smi_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * sli_497[k]
                   + f_1 * smh0_498[k]
                   - f_2 * smh1_498[k]
                   + f_3 * pc_y[k] * smi_665[k];

        t_857[k] = f_14 * sli_469[k]
                   + f_3 * pc_z[k] * smi_665[k];

        t_858[k] = f_16 * sli_499[k]
                   + f_4 * smh0_500[k]
                   - f_5 * smh1_500[k]
                   + f_3 * pc_y[k] * smi_667[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_y, sli_500, sli_501, sli_502, smh0_501, \
                         smh0_502, smh0_503, smh1_501, smh1_502, smh1_503, smi_668, smi_669, \
                         smi_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * sli_500[k]
                   + f_6 * smh0_501[k]
                   - f_7 * smh1_501[k]
                   + f_3 * pc_y[k] * smi_668[k];

        t_860[k] = f_16 * sli_501[k]
                   + f_8 * smh0_502[k]
                   - f_9 * smh1_502[k]
                   + f_3 * pc_y[k] * smi_669[k];

        t_861[k] = f_16 * sli_502[k]
                   + f_10 * smh0_503[k]
                   - f_11 * smh1_503[k]
                   + f_3 * pc_y[k] * smi_670[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pc_x, pc_y, pc_z, sli_475, sli_503, sli_672, \
                         smh0_503, smh0_504, smh1_503, smh1_504, smi_671, \
                         smi_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * sli_503[k]
                   + f_3 * pc_y[k] * smi_671[k];

        t_863[k] = f_14 * sli_475[k]
                   + f_1 * smh0_503[k]
                   - f_2 * smh1_503[k]
                   + f_3 * pc_z[k] * smi_671[k];

        t_864[k] = f_15 * sli_672[k]
                   + f_1 * smh0_504[k]
                   - f_2 * smh1_504[k]
                   + f_3 * pc_x[k] * smi_672[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pc_x, pc_y, pc_z, sli_476, sli_504, \
                         sli_506, sli_675, smh0_507, smh1_507, smi_672, smi_674, \
                         smi_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_15 * sli_504[k]
                   + f_3 * pc_y[k] * smi_672[k];

        t_866[k] = f_15 * sli_476[k]
                   + f_3 * pc_z[k] * smi_672[k];

        t_867[k] = f_15 * sli_675[k]
                   + f_4 * smh0_507[k]
                   - f_5 * smh1_507[k]
                   + f_3 * pc_x[k] * smi_675[k];

        t_868[k] = f_15 * sli_506[k]
                   + f_3 * pc_y[k] * smi_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pc_x, pc_z, sli_479, sli_677, sli_678, smh0_509, \
                         smh0_510, smh1_509, smh1_510, smi_675, smi_677, \
                         smi_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_15 * sli_677[k]
                   + f_4 * smh0_509[k]
                   - f_5 * smh1_509[k]
                   + f_3 * pc_x[k] * smi_677[k];

        t_870[k] = f_15 * sli_678[k]
                   + f_6 * smh0_510[k]
                   - f_7 * smh1_510[k]
                   + f_3 * pc_x[k] * smi_678[k];

        t_871[k] = f_15 * sli_479[k]
                   + f_3 * pc_z[k] * smi_675[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pc_x, pc_y, sli_509, sli_681, sli_682, smh0_513, \
                         smh0_514, smh1_513, smh1_514, smi_677, smi_681, \
                         smi_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_15 * sli_509[k]
                   + f_3 * pc_y[k] * smi_677[k];

        t_873[k] = f_15 * sli_681[k]
                   + f_6 * smh0_513[k]
                   - f_7 * smh1_513[k]
                   + f_3 * pc_x[k] * smi_681[k];

        t_874[k] = f_15 * sli_682[k]
                   + f_8 * smh0_514[k]
                   - f_9 * smh1_514[k]
                   + f_3 * pc_x[k] * smi_682[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pc_x, pc_y, pc_z, sli_482, sli_513, sli_684, \
                         smh0_516, smh1_516, smi_678, smi_681, \
                         smi_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_15 * sli_482[k]
                   + f_3 * pc_z[k] * smi_678[k];

        t_876[k] = f_15 * sli_684[k]
                   + f_8 * smh0_516[k]
                   - f_9 * smh1_516[k]
                   + f_3 * pc_x[k] * smi_684[k];

        t_877[k] = f_15 * sli_513[k]
                   + f_3 * pc_y[k] * smi_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_x, pc_z, sli_486, sli_686, sli_687, smh0_518, \
                         smh0_519, smh1_518, smh1_519, smi_682, smi_686, \
                         smi_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_15 * sli_686[k]
                   + f_8 * smh0_518[k]
                   - f_9 * smh1_518[k]
                   + f_3 * pc_x[k] * smi_686[k];

        t_879[k] = f_15 * sli_687[k]
                   + f_10 * smh0_519[k]
                   - f_11 * smh1_519[k]
                   + f_3 * pc_x[k] * smi_687[k];

        t_880[k] = f_15 * sli_486[k]
                   + f_3 * pc_z[k] * smi_682[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pc_x, pc_y, sli_518, sli_689, sli_690, smh0_521, \
                         smh0_522, smh1_521, smh1_522, smi_686, smi_689, \
                         smi_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_15 * sli_689[k]
                   + f_10 * smh0_521[k]
                   - f_11 * smh1_521[k]
                   + f_3 * pc_x[k] * smi_689[k];

        t_882[k] = f_15 * sli_690[k]
                   + f_10 * smh0_522[k]
                   - f_11 * smh1_522[k]
                   + f_3 * pc_x[k] * smi_690[k];

        t_883[k] = f_15 * sli_518[k]
                   + f_3 * pc_y[k] * smi_686[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, pc_x, sli_692, sli_693, sli_694, sli_695, \
                         smh0_524, smh1_524, smi_692, smi_693, smi_694, \
                         smi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_15 * sli_692[k]
                   + f_10 * smh0_524[k]
                   - f_11 * smh1_524[k]
                   + f_3 * pc_x[k] * smi_692[k];

        t_885[k] = f_15 * sli_693[k]
                   + f_3 * pc_x[k] * smi_693[k];

        t_886[k] = f_15 * sli_694[k]
                   + f_3 * pc_x[k] * smi_694[k];

        t_887[k] = f_15 * sli_695[k]
                   + f_3 * pc_x[k] * smi_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pc_x, sli_696, sli_697, sli_698, sli_699, \
                         smi_696, smi_697, smi_698, smi_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_15 * sli_696[k]
                   + f_3 * pc_x[k] * smi_696[k];

        t_889[k] = f_15 * sli_697[k]
                   + f_3 * pc_x[k] * smi_697[k];

        t_890[k] = f_15 * sli_698[k]
                   + f_3 * pc_x[k] * smi_698[k];

        t_891[k] = f_15 * sli_699[k]
                   + f_3 * pc_x[k] * smi_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pc_y, pc_z, sli_497, sli_525, sli_527, smh0_519, \
                         smh0_521, smh1_519, smh1_521, smi_693, \
                         smi_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_15 * sli_525[k]
                   + f_1 * smh0_519[k]
                   - f_2 * smh1_519[k]
                   + f_3 * pc_y[k] * smi_693[k];

        t_893[k] = f_15 * sli_497[k]
                   + f_3 * pc_z[k] * smi_693[k];

        t_894[k] = f_15 * sli_527[k]
                   + f_4 * smh0_521[k]
                   - f_5 * smh1_521[k]
                   + f_3 * pc_y[k] * smi_695[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pc_y, sli_528, sli_529, sli_530, smh0_522, \
                         smh0_523, smh0_524, smh1_522, smh1_523, smh1_524, smi_696, smi_697, \
                         smi_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_15 * sli_528[k]
                   + f_6 * smh0_522[k]
                   - f_7 * smh1_522[k]
                   + f_3 * pc_y[k] * smi_696[k];

        t_896[k] = f_15 * sli_529[k]
                   + f_8 * smh0_523[k]
                   - f_9 * smh1_523[k]
                   + f_3 * pc_y[k] * smi_697[k];

        t_897[k] = f_15 * sli_530[k]
                   + f_10 * smh0_524[k]
                   - f_11 * smh1_524[k]
                   + f_3 * pc_y[k] * smi_698[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_720 = buffer.data(slk0 + 720);
    const auto *slk0_723 = buffer.data(slk0 + 723);
    const auto *slk0_725 = buffer.data(slk0 + 725);
    const auto *slk0_726 = buffer.data(slk0 + 726);
    const auto *slk0_729 = buffer.data(slk0 + 729);
    const auto *slk0_730 = buffer.data(slk0 + 730);
    const auto *slk0_732 = buffer.data(slk0 + 732);
    const auto *slk0_734 = buffer.data(slk0 + 734);
    const auto *slk0_735 = buffer.data(slk0 + 735);
    const auto *slk0_737 = buffer.data(slk0 + 737);
    const auto *slk0_738 = buffer.data(slk0 + 738);
    const auto *slk0_740 = buffer.data(slk0 + 740);
    const auto *slk0_755 = buffer.data(slk0 + 755);

    const auto *sli_503 = buffer.data(sli + 503);
    const auto *sli_504 = buffer.data(sli + 504);
    const auto *sli_507 = buffer.data(sli + 507);
    const auto *sli_510 = buffer.data(sli + 510);
    const auto *sli_514 = buffer.data(sli + 514);
    const auto *sli_525 = buffer.data(sli + 525);
    const auto *sli_531 = buffer.data(sli + 531);
    const auto *sli_532 = buffer.data(sli + 532);
    const auto *sli_534 = buffer.data(sli + 534);
    const auto *sli_535 = buffer.data(sli + 535);
    const auto *sli_537 = buffer.data(sli + 537);
    const auto *sli_538 = buffer.data(sli + 538);
    const auto *sli_541 = buffer.data(sli + 541);
    const auto *sli_542 = buffer.data(sli + 542);
    const auto *sli_546 = buffer.data(sli + 546);
    const auto *sli_553 = buffer.data(sli + 553);
    const auto *sli_555 = buffer.data(sli + 555);
    const auto *sli_556 = buffer.data(sli + 556);
    const auto *sli_557 = buffer.data(sli + 557);
    const auto *sli_558 = buffer.data(sli + 558);
    const auto *sli_559 = buffer.data(sli + 559);
    const auto *sli_560 = buffer.data(sli + 560);
    const auto *sli_561 = buffer.data(sli + 561);
    const auto *sli_562 = buffer.data(sli + 562);
    const auto *sli_563 = buffer.data(sli + 563);
    const auto *sli_565 = buffer.data(sli + 565);
    const auto *sli_566 = buffer.data(sli + 566);
    const auto *sli_568 = buffer.data(sli + 568);
    const auto *sli_569 = buffer.data(sli + 569);
    const auto *sli_570 = buffer.data(sli + 570);
    const auto *sli_572 = buffer.data(sli + 572);
    const auto *sli_573 = buffer.data(sli + 573);
    const auto *sli_574 = buffer.data(sli + 574);
    const auto *sli_581 = buffer.data(sli + 581);
    const auto *sli_583 = buffer.data(sli + 583);
    const auto *sli_584 = buffer.data(sli + 584);
    const auto *sli_585 = buffer.data(sli + 585);
    const auto *sli_586 = buffer.data(sli + 586);
    const auto *sli_587 = buffer.data(sli + 587);
    const auto *sli_700 = buffer.data(sli + 700);
    const auto *sli_703 = buffer.data(sli + 703);
    const auto *sli_705 = buffer.data(sli + 705);
    const auto *sli_706 = buffer.data(sli + 706);
    const auto *sli_709 = buffer.data(sli + 709);
    const auto *sli_710 = buffer.data(sli + 710);
    const auto *sli_712 = buffer.data(sli + 712);
    const auto *sli_714 = buffer.data(sli + 714);
    const auto *sli_715 = buffer.data(sli + 715);
    const auto *sli_717 = buffer.data(sli + 717);
    const auto *sli_718 = buffer.data(sli + 718);
    const auto *sli_720 = buffer.data(sli + 720);
    const auto *sli_721 = buffer.data(sli + 721);
    const auto *sli_722 = buffer.data(sli + 722);
    const auto *sli_723 = buffer.data(sli + 723);
    const auto *sli_724 = buffer.data(sli + 724);
    const auto *sli_725 = buffer.data(sli + 725);
    const auto *sli_726 = buffer.data(sli + 726);
    const auto *sli_727 = buffer.data(sli + 727);
    const auto *sli_749 = buffer.data(sli + 749);
    const auto *sli_750 = buffer.data(sli + 750);
    const auto *sli_751 = buffer.data(sli + 751);
    const auto *sli_752 = buffer.data(sli + 752);
    const auto *sli_753 = buffer.data(sli + 753);
    const auto *sli_754 = buffer.data(sli + 754);
    const auto *sli_755 = buffer.data(sli + 755);
    const auto *sli_756 = buffer.data(sli + 756);
    const auto *sli_759 = buffer.data(sli + 759);
    const auto *sli_761 = buffer.data(sli + 761);
    const auto *sli_762 = buffer.data(sli + 762);
    const auto *sli_765 = buffer.data(sli + 765);
    const auto *sli_766 = buffer.data(sli + 766);
    const auto *sli_768 = buffer.data(sli + 768);
    const auto *sli_770 = buffer.data(sli + 770);
    const auto *sli_771 = buffer.data(sli + 771);
    const auto *sli_773 = buffer.data(sli + 773);
    const auto *sli_774 = buffer.data(sli + 774);
    const auto *sli_776 = buffer.data(sli + 776);
    const auto *sli_777 = buffer.data(sli + 777);
    const auto *sli_778 = buffer.data(sli + 778);
    const auto *sli_779 = buffer.data(sli + 779);
    const auto *sli_780 = buffer.data(sli + 780);
    const auto *sli_781 = buffer.data(sli + 781);
    const auto *sli_782 = buffer.data(sli + 782);
    const auto *sli_783 = buffer.data(sli + 783);

    const auto *slk1_720 = buffer.data(slk1 + 720);
    const auto *slk1_723 = buffer.data(slk1 + 723);
    const auto *slk1_725 = buffer.data(slk1 + 725);
    const auto *slk1_726 = buffer.data(slk1 + 726);
    const auto *slk1_729 = buffer.data(slk1 + 729);
    const auto *slk1_730 = buffer.data(slk1 + 730);
    const auto *slk1_732 = buffer.data(slk1 + 732);
    const auto *slk1_734 = buffer.data(slk1 + 734);
    const auto *slk1_735 = buffer.data(slk1 + 735);
    const auto *slk1_737 = buffer.data(slk1 + 737);
    const auto *slk1_738 = buffer.data(slk1 + 738);
    const auto *slk1_740 = buffer.data(slk1 + 740);
    const auto *slk1_755 = buffer.data(slk1 + 755);

    const auto *smh0_524 = buffer.data(smh0 + 524);
    const auto *smh0_525 = buffer.data(smh0 + 525);
    const auto *smh0_528 = buffer.data(smh0 + 528);
    const auto *smh0_530 = buffer.data(smh0 + 530);
    const auto *smh0_531 = buffer.data(smh0 + 531);
    const auto *smh0_534 = buffer.data(smh0 + 534);
    const auto *smh0_535 = buffer.data(smh0 + 535);
    const auto *smh0_537 = buffer.data(smh0 + 537);
    const auto *smh0_539 = buffer.data(smh0 + 539);
    const auto *smh0_540 = buffer.data(smh0 + 540);
    const auto *smh0_542 = buffer.data(smh0 + 542);
    const auto *smh0_543 = buffer.data(smh0 + 543);
    const auto *smh0_544 = buffer.data(smh0 + 544);
    const auto *smh0_545 = buffer.data(smh0 + 545);
    const auto *smh0_561 = buffer.data(smh0 + 561);
    const auto *smh0_563 = buffer.data(smh0 + 563);
    const auto *smh0_564 = buffer.data(smh0 + 564);
    const auto *smh0_565 = buffer.data(smh0 + 565);
    const auto *smh0_566 = buffer.data(smh0 + 566);
    const auto *smh0_567 = buffer.data(smh0 + 567);
    const auto *smh0_570 = buffer.data(smh0 + 570);
    const auto *smh0_572 = buffer.data(smh0 + 572);
    const auto *smh0_573 = buffer.data(smh0 + 573);
    const auto *smh0_576 = buffer.data(smh0 + 576);
    const auto *smh0_577 = buffer.data(smh0 + 577);
    const auto *smh0_579 = buffer.data(smh0 + 579);
    const auto *smh0_581 = buffer.data(smh0 + 581);
    const auto *smh0_582 = buffer.data(smh0 + 582);
    const auto *smh0_584 = buffer.data(smh0 + 584);
    const auto *smh0_585 = buffer.data(smh0 + 585);
    const auto *smh0_586 = buffer.data(smh0 + 586);
    const auto *smh0_587 = buffer.data(smh0 + 587);

    const auto *smh1_524 = buffer.data(smh1 + 524);
    const auto *smh1_525 = buffer.data(smh1 + 525);
    const auto *smh1_528 = buffer.data(smh1 + 528);
    const auto *smh1_530 = buffer.data(smh1 + 530);
    const auto *smh1_531 = buffer.data(smh1 + 531);
    const auto *smh1_534 = buffer.data(smh1 + 534);
    const auto *smh1_535 = buffer.data(smh1 + 535);
    const auto *smh1_537 = buffer.data(smh1 + 537);
    const auto *smh1_539 = buffer.data(smh1 + 539);
    const auto *smh1_540 = buffer.data(smh1 + 540);
    const auto *smh1_542 = buffer.data(smh1 + 542);
    const auto *smh1_543 = buffer.data(smh1 + 543);
    const auto *smh1_544 = buffer.data(smh1 + 544);
    const auto *smh1_545 = buffer.data(smh1 + 545);
    const auto *smh1_561 = buffer.data(smh1 + 561);
    const auto *smh1_563 = buffer.data(smh1 + 563);
    const auto *smh1_564 = buffer.data(smh1 + 564);
    const auto *smh1_565 = buffer.data(smh1 + 565);
    const auto *smh1_566 = buffer.data(smh1 + 566);
    const auto *smh1_567 = buffer.data(smh1 + 567);
    const auto *smh1_570 = buffer.data(smh1 + 570);
    const auto *smh1_572 = buffer.data(smh1 + 572);
    const auto *smh1_573 = buffer.data(smh1 + 573);
    const auto *smh1_576 = buffer.data(smh1 + 576);
    const auto *smh1_577 = buffer.data(smh1 + 577);
    const auto *smh1_579 = buffer.data(smh1 + 579);
    const auto *smh1_581 = buffer.data(smh1 + 581);
    const auto *smh1_582 = buffer.data(smh1 + 582);
    const auto *smh1_584 = buffer.data(smh1 + 584);
    const auto *smh1_585 = buffer.data(smh1 + 585);
    const auto *smh1_586 = buffer.data(smh1 + 586);
    const auto *smh1_587 = buffer.data(smh1 + 587);

    const auto *smi_699 = buffer.data(smi + 699);
    const auto *smi_700 = buffer.data(smi + 700);
    const auto *smi_702 = buffer.data(smi + 702);
    const auto *smi_703 = buffer.data(smi + 703);
    const auto *smi_705 = buffer.data(smi + 705);
    const auto *smi_706 = buffer.data(smi + 706);
    const auto *smi_709 = buffer.data(smi + 709);
    const auto *smi_710 = buffer.data(smi + 710);
    const auto *smi_712 = buffer.data(smi + 712);
    const auto *smi_714 = buffer.data(smi + 714);
    const auto *smi_715 = buffer.data(smi + 715);
    const auto *smi_717 = buffer.data(smi + 717);
    const auto *smi_718 = buffer.data(smi + 718);
    const auto *smi_720 = buffer.data(smi + 720);
    const auto *smi_721 = buffer.data(smi + 721);
    const auto *smi_722 = buffer.data(smi + 722);
    const auto *smi_723 = buffer.data(smi + 723);
    const auto *smi_724 = buffer.data(smi + 724);
    const auto *smi_725 = buffer.data(smi + 725);
    const auto *smi_726 = buffer.data(smi + 726);
    const auto *smi_727 = buffer.data(smi + 727);
    const auto *smi_728 = buffer.data(smi + 728);
    const auto *smi_730 = buffer.data(smi + 730);
    const auto *smi_731 = buffer.data(smi + 731);
    const auto *smi_733 = buffer.data(smi + 733);
    const auto *smi_734 = buffer.data(smi + 734);
    const auto *smi_737 = buffer.data(smi + 737);
    const auto *smi_738 = buffer.data(smi + 738);
    const auto *smi_742 = buffer.data(smi + 742);
    const auto *smi_749 = buffer.data(smi + 749);
    const auto *smi_750 = buffer.data(smi + 750);
    const auto *smi_751 = buffer.data(smi + 751);
    const auto *smi_752 = buffer.data(smi + 752);
    const auto *smi_753 = buffer.data(smi + 753);
    const auto *smi_754 = buffer.data(smi + 754);
    const auto *smi_755 = buffer.data(smi + 755);
    const auto *smi_756 = buffer.data(smi + 756);
    const auto *smi_758 = buffer.data(smi + 758);
    const auto *smi_759 = buffer.data(smi + 759);
    const auto *smi_761 = buffer.data(smi + 761);
    const auto *smi_762 = buffer.data(smi + 762);
    const auto *smi_765 = buffer.data(smi + 765);
    const auto *smi_766 = buffer.data(smi + 766);
    const auto *smi_768 = buffer.data(smi + 768);
    const auto *smi_770 = buffer.data(smi + 770);
    const auto *smi_771 = buffer.data(smi + 771);
    const auto *smi_773 = buffer.data(smi + 773);
    const auto *smi_774 = buffer.data(smi + 774);
    const auto *smi_776 = buffer.data(smi + 776);
    const auto *smi_777 = buffer.data(smi + 777);
    const auto *smi_778 = buffer.data(smi + 778);
    const auto *smi_779 = buffer.data(smi + 779);
    const auto *smi_780 = buffer.data(smi + 780);
    const auto *smi_781 = buffer.data(smi + 781);
    const auto *smi_782 = buffer.data(smi + 782);
    const auto *smi_783 = buffer.data(smi + 783);

#pragma omp simd aligned(t_898, t_899, t_900, pc_x, pc_y, pc_z, sli_503, sli_531, sli_700, \
                         smh0_524, smh0_525, smh1_524, smh1_525, smi_699, \
                         smi_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * sli_531[k]
                   + f_3 * pc_y[k] * smi_699[k];

        t_899[k] = f_15 * sli_503[k]
                   + f_1 * smh0_524[k]
                   - f_2 * smh1_524[k]
                   + f_3 * pc_z[k] * smi_699[k];

        t_900[k] = f_15 * sli_700[k]
                   + f_1 * smh0_525[k]
                   - f_2 * smh1_525[k]
                   + f_3 * pc_x[k] * smi_700[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pc_x, pc_y, pc_z, sli_504, sli_532, \
                         sli_534, sli_703, smh0_528, smh1_528, smi_700, smi_702, \
                         smi_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_14 * sli_532[k]
                   + f_3 * pc_y[k] * smi_700[k];

        t_902[k] = f_16 * sli_504[k]
                   + f_3 * pc_z[k] * smi_700[k];

        t_903[k] = f_15 * sli_703[k]
                   + f_4 * smh0_528[k]
                   - f_5 * smh1_528[k]
                   + f_3 * pc_x[k] * smi_703[k];

        t_904[k] = f_14 * sli_534[k]
                   + f_3 * pc_y[k] * smi_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, sli_507, sli_705, sli_706, smh0_530, \
                         smh0_531, smh1_530, smh1_531, smi_703, smi_705, \
                         smi_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_15 * sli_705[k]
                   + f_4 * smh0_530[k]
                   - f_5 * smh1_530[k]
                   + f_3 * pc_x[k] * smi_705[k];

        t_906[k] = f_15 * sli_706[k]
                   + f_6 * smh0_531[k]
                   - f_7 * smh1_531[k]
                   + f_3 * pc_x[k] * smi_706[k];

        t_907[k] = f_16 * sli_507[k]
                   + f_3 * pc_z[k] * smi_703[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, sli_537, sli_709, sli_710, smh0_534, \
                         smh0_535, smh1_534, smh1_535, smi_705, smi_709, \
                         smi_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * sli_537[k]
                   + f_3 * pc_y[k] * smi_705[k];

        t_909[k] = f_15 * sli_709[k]
                   + f_6 * smh0_534[k]
                   - f_7 * smh1_534[k]
                   + f_3 * pc_x[k] * smi_709[k];

        t_910[k] = f_15 * sli_710[k]
                   + f_8 * smh0_535[k]
                   - f_9 * smh1_535[k]
                   + f_3 * pc_x[k] * smi_710[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pc_x, pc_y, pc_z, sli_510, sli_541, sli_712, \
                         smh0_537, smh1_537, smi_706, smi_709, \
                         smi_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_16 * sli_510[k]
                   + f_3 * pc_z[k] * smi_706[k];

        t_912[k] = f_15 * sli_712[k]
                   + f_8 * smh0_537[k]
                   - f_9 * smh1_537[k]
                   + f_3 * pc_x[k] * smi_712[k];

        t_913[k] = f_14 * sli_541[k]
                   + f_3 * pc_y[k] * smi_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pc_x, pc_z, sli_514, sli_714, sli_715, smh0_539, \
                         smh0_540, smh1_539, smh1_540, smi_710, smi_714, \
                         smi_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_15 * sli_714[k]
                   + f_8 * smh0_539[k]
                   - f_9 * smh1_539[k]
                   + f_3 * pc_x[k] * smi_714[k];

        t_915[k] = f_15 * sli_715[k]
                   + f_10 * smh0_540[k]
                   - f_11 * smh1_540[k]
                   + f_3 * pc_x[k] * smi_715[k];

        t_916[k] = f_16 * sli_514[k]
                   + f_3 * pc_z[k] * smi_710[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pc_x, pc_y, sli_546, sli_717, sli_718, smh0_542, \
                         smh0_543, smh1_542, smh1_543, smi_714, smi_717, \
                         smi_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_15 * sli_717[k]
                   + f_10 * smh0_542[k]
                   - f_11 * smh1_542[k]
                   + f_3 * pc_x[k] * smi_717[k];

        t_918[k] = f_15 * sli_718[k]
                   + f_10 * smh0_543[k]
                   - f_11 * smh1_543[k]
                   + f_3 * pc_x[k] * smi_718[k];

        t_919[k] = f_14 * sli_546[k]
                   + f_3 * pc_y[k] * smi_714[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pc_x, sli_720, sli_721, sli_722, sli_723, \
                         smh0_545, smh1_545, smi_720, smi_721, smi_722, \
                         smi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_15 * sli_720[k]
                   + f_10 * smh0_545[k]
                   - f_11 * smh1_545[k]
                   + f_3 * pc_x[k] * smi_720[k];

        t_921[k] = f_15 * sli_721[k]
                   + f_3 * pc_x[k] * smi_721[k];

        t_922[k] = f_15 * sli_722[k]
                   + f_3 * pc_x[k] * smi_722[k];

        t_923[k] = f_15 * sli_723[k]
                   + f_3 * pc_x[k] * smi_723[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, sli_724, sli_725, sli_726, sli_727, \
                         smi_724, smi_725, smi_726, smi_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_15 * sli_724[k]
                   + f_3 * pc_x[k] * smi_724[k];

        t_925[k] = f_15 * sli_725[k]
                   + f_3 * pc_x[k] * smi_725[k];

        t_926[k] = f_15 * sli_726[k]
                   + f_3 * pc_x[k] * smi_726[k];

        t_927[k] = f_15 * sli_727[k]
                   + f_3 * pc_x[k] * smi_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_y, pc_z, sli_525, sli_553, sli_555, smh0_540, \
                         smh0_542, smh1_540, smh1_542, smi_721, \
                         smi_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_14 * sli_553[k]
                   + f_1 * smh0_540[k]
                   - f_2 * smh1_540[k]
                   + f_3 * pc_y[k] * smi_721[k];

        t_929[k] = f_16 * sli_525[k]
                   + f_3 * pc_z[k] * smi_721[k];

        t_930[k] = f_14 * sli_555[k]
                   + f_4 * smh0_542[k]
                   - f_5 * smh1_542[k]
                   + f_3 * pc_y[k] * smi_723[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, pc_y, sli_556, sli_557, sli_558, smh0_543, \
                         smh0_544, smh0_545, smh1_543, smh1_544, smh1_545, smi_724, smi_725, \
                         smi_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_14 * sli_556[k]
                   + f_6 * smh0_543[k]
                   - f_7 * smh1_543[k]
                   + f_3 * pc_y[k] * smi_724[k];

        t_932[k] = f_14 * sli_557[k]
                   + f_8 * smh0_544[k]
                   - f_9 * smh1_544[k]
                   + f_3 * pc_y[k] * smi_725[k];

        t_933[k] = f_14 * sli_558[k]
                   + f_10 * smh0_545[k]
                   - f_11 * smh1_545[k]
                   + f_3 * pc_y[k] * smi_726[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pb_y, pc_y, pc_z, slk0_720, sli_531, \
                         sli_559, sli_560, slk1_720, smh0_545, smh1_545, smi_727, \
                         smi_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_14 * sli_559[k]
                   + f_3 * pc_y[k] * smi_727[k];

        t_935[k] = f_16 * sli_531[k]
                   + f_1 * smh0_545[k]
                   - f_2 * smh1_545[k]
                   + f_3 * pc_z[k] * smi_727[k];

        t_936[k] = pb_y[k] * slk0_720[k]
                   - f_12 * pc_y[k] * slk1_720[k];

        t_937[k] = f_13 * sli_560[k]
                   + f_3 * pc_y[k] * smi_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pb_y, pc_y, pc_z, slk0_723, slk0_725, \
                         sli_532, sli_561, sli_562, slk1_723, slk1_725, smi_728, \
                         smi_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_17 * sli_532[k]
                   + f_3 * pc_z[k] * smi_728[k];

        t_939[k] = pb_y[k] * slk0_723[k]
                   + f_14 * sli_561[k]
                   - f_12 * pc_y[k] * slk1_723[k];

        t_940[k] = f_13 * sli_562[k]
                   + f_3 * pc_y[k] * smi_730[k];

        t_941[k] = pb_y[k] * slk0_725[k]
                   - f_12 * pc_y[k] * slk1_725[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pb_y, pc_y, pc_z, slk0_726, slk0_729, \
                         sli_535, sli_563, sli_565, slk1_726, slk1_729, smi_731, \
                         smi_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = pb_y[k] * slk0_726[k]
                   + f_15 * sli_563[k]
                   - f_12 * pc_y[k] * slk1_726[k];

        t_943[k] = f_17 * sli_535[k]
                   + f_3 * pc_z[k] * smi_731[k];

        t_944[k] = f_13 * sli_565[k]
                   + f_3 * pc_y[k] * smi_733[k];

        t_945[k] = pb_y[k] * slk0_729[k]
                   - f_12 * pc_y[k] * slk1_729[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pb_y, pc_y, pc_z, slk0_730, slk0_732, sli_538, \
                         sli_566, sli_568, slk1_730, slk1_732, \
                         smi_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = pb_y[k] * slk0_730[k]
                   + f_16 * sli_566[k]
                   - f_12 * pc_y[k] * slk1_730[k];

        t_947[k] = f_17 * sli_538[k]
                   + f_3 * pc_z[k] * smi_734[k];

        t_948[k] = pb_y[k] * slk0_732[k]
                   + f_14 * sli_568[k]
                   - f_12 * pc_y[k] * slk1_732[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pb_y, pc_y, pc_z, slk0_734, slk0_735, \
                         sli_542, sli_569, sli_570, slk1_734, slk1_735, smi_737, \
                         smi_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_13 * sli_569[k]
                   + f_3 * pc_y[k] * smi_737[k];

        t_950[k] = pb_y[k] * slk0_734[k]
                   - f_12 * pc_y[k] * slk1_734[k];

        t_951[k] = pb_y[k] * slk0_735[k]
                   + f_17 * sli_570[k]
                   - f_12 * pc_y[k] * slk1_735[k];

        t_952[k] = f_17 * sli_542[k]
                   + f_3 * pc_z[k] * smi_738[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pb_y, pc_y, slk0_737, slk0_738, slk0_740, \
                         sli_572, sli_573, sli_574, slk1_737, slk1_738, slk1_740, \
                         smi_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = pb_y[k] * slk0_737[k]
                   + f_15 * sli_572[k]
                   - f_12 * pc_y[k] * slk1_737[k];

        t_954[k] = pb_y[k] * slk0_738[k]
                   + f_14 * sli_573[k]
                   - f_12 * pc_y[k] * slk1_738[k];

        t_955[k] = f_13 * sli_574[k]
                   + f_3 * pc_y[k] * smi_742[k];

        t_956[k] = pb_y[k] * slk0_740[k]
                   - f_12 * pc_y[k] * slk1_740[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, pc_x, sli_749, sli_750, sli_751, \
                         sli_752, sli_753, smi_749, smi_750, smi_751, smi_752, \
                         smi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_15 * sli_749[k]
                   + f_3 * pc_x[k] * smi_749[k];

        t_958[k] = f_15 * sli_750[k]
                   + f_3 * pc_x[k] * smi_750[k];

        t_959[k] = f_15 * sli_751[k]
                   + f_3 * pc_x[k] * smi_751[k];

        t_960[k] = f_15 * sli_752[k]
                   + f_3 * pc_x[k] * smi_752[k];

        t_961[k] = f_15 * sli_753[k]
                   + f_3 * pc_x[k] * smi_753[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pc_x, pc_y, pc_z, sli_553, sli_581, \
                         sli_754, sli_755, smh0_561, smh1_561, smi_749, smi_754, \
                         smi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_15 * sli_754[k]
                   + f_3 * pc_x[k] * smi_754[k];

        t_963[k] = f_15 * sli_755[k]
                   + f_3 * pc_x[k] * smi_755[k];

        t_964[k] = f_13 * sli_581[k]
                   + f_1 * smh0_561[k]
                   - f_2 * smh1_561[k]
                   + f_3 * pc_y[k] * smi_749[k];

        t_965[k] = f_17 * sli_553[k]
                   + f_3 * pc_z[k] * smi_749[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_y, sli_583, sli_584, sli_585, smh0_563, \
                         smh0_564, smh0_565, smh1_563, smh1_564, smh1_565, smi_751, smi_752, \
                         smi_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_13 * sli_583[k]
                   + f_4 * smh0_563[k]
                   - f_5 * smh1_563[k]
                   + f_3 * pc_y[k] * smi_751[k];

        t_967[k] = f_13 * sli_584[k]
                   + f_6 * smh0_564[k]
                   - f_7 * smh1_564[k]
                   + f_3 * pc_y[k] * smi_752[k];

        t_968[k] = f_13 * sli_585[k]
                   + f_8 * smh0_565[k]
                   - f_9 * smh1_565[k]
                   + f_3 * pc_y[k] * smi_753[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pb_y, pc_y, slk0_755, sli_586, sli_587, \
                         slk1_755, smh0_566, smh1_566, smi_754, \
                         smi_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_13 * sli_586[k]
                   + f_10 * smh0_566[k]
                   - f_11 * smh1_566[k]
                   + f_3 * pc_y[k] * smi_754[k];

        t_970[k] = f_13 * sli_587[k]
                   + f_3 * pc_y[k] * smi_755[k];

        t_971[k] = pb_y[k] * slk0_755[k]
                   - f_12 * pc_y[k] * slk1_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, pc_x, pc_y, pc_z, sli_560, sli_756, \
                         sli_759, smh0_567, smh0_570, smh1_567, smh1_570, smi_756, \
                         smi_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_15 * sli_756[k]
                   + f_1 * smh0_567[k]
                   - f_2 * smh1_567[k]
                   + f_3 * pc_x[k] * smi_756[k];

        t_973[k] = f_3 * pc_y[k] * smi_756[k];

        t_974[k] = f_20 * sli_560[k]
                   + f_3 * pc_z[k] * smi_756[k];

        t_975[k] = f_15 * sli_759[k]
                   + f_4 * smh0_570[k]
                   - f_5 * smh1_570[k]
                   + f_3 * pc_x[k] * smi_759[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_x, pc_y, sli_761, sli_762, smh0_572, \
                         smh0_573, smh1_572, smh1_573, smi_758, smi_761, \
                         smi_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_3 * pc_y[k] * smi_758[k];

        t_977[k] = f_15 * sli_761[k]
                   + f_4 * smh0_572[k]
                   - f_5 * smh1_572[k]
                   + f_3 * pc_x[k] * smi_761[k];

        t_978[k] = f_15 * sli_762[k]
                   + f_6 * smh0_573[k]
                   - f_7 * smh1_573[k]
                   + f_3 * pc_x[k] * smi_762[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, pc_x, pc_y, pc_z, sli_563, sli_765, smh0_576, \
                         smh1_576, smi_759, smi_761, smi_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_20 * sli_563[k]
                   + f_3 * pc_z[k] * smi_759[k];

        t_980[k] = f_3 * pc_y[k] * smi_761[k];

        t_981[k] = f_15 * sli_765[k]
                   + f_6 * smh0_576[k]
                   - f_7 * smh1_576[k]
                   + f_3 * pc_x[k] * smi_765[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, pc_x, pc_z, sli_566, sli_766, sli_768, smh0_577, \
                         smh0_579, smh1_577, smh1_579, smi_762, smi_766, \
                         smi_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_15 * sli_766[k]
                   + f_8 * smh0_577[k]
                   - f_9 * smh1_577[k]
                   + f_3 * pc_x[k] * smi_766[k];

        t_983[k] = f_20 * sli_566[k]
                   + f_3 * pc_z[k] * smi_762[k];

        t_984[k] = f_15 * sli_768[k]
                   + f_8 * smh0_579[k]
                   - f_9 * smh1_579[k]
                   + f_3 * pc_x[k] * smi_768[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pc_x, pc_y, sli_770, sli_771, smh0_581, \
                         smh0_582, smh1_581, smh1_582, smi_765, smi_770, \
                         smi_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_3 * pc_y[k] * smi_765[k];

        t_986[k] = f_15 * sli_770[k]
                   + f_8 * smh0_581[k]
                   - f_9 * smh1_581[k]
                   + f_3 * pc_x[k] * smi_770[k];

        t_987[k] = f_15 * sli_771[k]
                   + f_10 * smh0_582[k]
                   - f_11 * smh1_582[k]
                   + f_3 * pc_x[k] * smi_771[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, pc_x, pc_z, sli_570, sli_773, sli_774, smh0_584, \
                         smh0_585, smh1_584, smh1_585, smi_766, smi_773, \
                         smi_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_20 * sli_570[k]
                   + f_3 * pc_z[k] * smi_766[k];

        t_989[k] = f_15 * sli_773[k]
                   + f_10 * smh0_584[k]
                   - f_11 * smh1_584[k]
                   + f_3 * pc_x[k] * smi_773[k];

        t_990[k] = f_15 * sli_774[k]
                   + f_10 * smh0_585[k]
                   - f_11 * smh1_585[k]
                   + f_3 * pc_x[k] * smi_774[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pc_x, pc_y, sli_776, sli_777, sli_778, \
                         smh0_587, smh1_587, smi_770, smi_776, smi_777, \
                         smi_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_3 * pc_y[k] * smi_770[k];

        t_992[k] = f_15 * sli_776[k]
                   + f_10 * smh0_587[k]
                   - f_11 * smh1_587[k]
                   + f_3 * pc_x[k] * smi_776[k];

        t_993[k] = f_15 * sli_777[k]
                   + f_3 * pc_x[k] * smi_777[k];

        t_994[k] = f_15 * sli_778[k]
                   + f_3 * pc_x[k] * smi_778[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pc_x, sli_779, sli_780, sli_781, \
                         sli_782, sli_783, smi_779, smi_780, smi_781, smi_782, \
                         smi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_15 * sli_779[k]
                   + f_3 * pc_x[k] * smi_779[k];

        t_996[k] = f_15 * sli_780[k]
                   + f_3 * pc_x[k] * smi_780[k];

        t_997[k] = f_15 * sli_781[k]
                   + f_3 * pc_x[k] * smi_781[k];

        t_998[k] = f_15 * sli_782[k]
                   + f_3 * pc_x[k] * smi_782[k];

        t_999[k] = f_15 * sli_783[k]
                   + f_3 * pc_x[k] * smi_783[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pc_y, pc_z, sli_581, smh0_582, \
                         smh0_584, smh0_585, smh1_582, smh1_584, smh1_585, smi_777, smi_779, \
                         smi_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_1 * smh0_582[k]
                    - f_2 * smh1_582[k]
                    + f_3 * pc_y[k] * smi_777[k];

        t_1001[k] = f_20 * sli_581[k]
                    + f_3 * pc_z[k] * smi_777[k];

        t_1002[k] = f_4 * smh0_584[k]
                    - f_5 * smh1_584[k]
                    + f_3 * pc_y[k] * smi_779[k];

        t_1003[k] = f_6 * smh0_585[k]
                    - f_7 * smh1_585[k]
                    + f_3 * pc_y[k] * smi_780[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, sli_587, smh0_586, \
                         smh0_587, smh1_586, smh1_587, smi_781, smi_782, \
                         smi_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_8 * smh0_586[k]
                    - f_9 * smh1_586[k]
                    + f_3 * pc_y[k] * smi_781[k];

        t_1005[k] = f_10 * smh0_587[k]
                    - f_11 * smh1_587[k]
                    + f_3 * pc_y[k] * smi_782[k];

        t_1006[k] = f_3 * pc_y[k] * smi_783[k];

        t_1007[k] = f_20 * sli_587[k]
                    + f_1 * smh0_587[k]
                    - f_2 * smh1_587[k]
                    + f_3 * pc_z[k] * smi_783[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t slk0,
                                                          const size_t sli, const size_t slk1,
                                                          const size_t smh0, const size_t smh1,
                                                          const size_t smi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_756 = buffer.data(slk0 + 756);
    const auto *slk0_759 = buffer.data(slk0 + 759);
    const auto *slk0_762 = buffer.data(slk0 + 762);
    const auto *slk0_766 = buffer.data(slk0 + 766);
    const auto *slk0_768 = buffer.data(slk0 + 768);
    const auto *slk0_771 = buffer.data(slk0 + 771);
    const auto *slk0_773 = buffer.data(slk0 + 773);
    const auto *slk0_774 = buffer.data(slk0 + 774);
    const auto *slk0_784 = buffer.data(slk0 + 784);

    const auto *sli_588 = buffer.data(sli + 588);
    const auto *sli_590 = buffer.data(sli + 590);
    const auto *sli_591 = buffer.data(sli + 591);
    const auto *sli_593 = buffer.data(sli + 593);
    const auto *sli_594 = buffer.data(sli + 594);
    const auto *sli_595 = buffer.data(sli + 595);
    const auto *sli_597 = buffer.data(sli + 597);
    const auto *sli_598 = buffer.data(sli + 598);
    const auto *sli_599 = buffer.data(sli + 599);
    const auto *sli_600 = buffer.data(sli + 600);
    const auto *sli_602 = buffer.data(sli + 602);
    const auto *sli_609 = buffer.data(sli + 609);
    const auto *sli_611 = buffer.data(sli + 611);
    const auto *sli_612 = buffer.data(sli + 612);
    const auto *sli_613 = buffer.data(sli + 613);
    const auto *sli_614 = buffer.data(sli + 614);
    const auto *sli_615 = buffer.data(sli + 615);
    const auto *sli_616 = buffer.data(sli + 616);
    const auto *sli_618 = buffer.data(sli + 618);
    const auto *sli_619 = buffer.data(sli + 619);
    const auto *sli_621 = buffer.data(sli + 621);
    const auto *sli_622 = buffer.data(sli + 622);
    const auto *sli_625 = buffer.data(sli + 625);
    const auto *sli_626 = buffer.data(sli + 626);
    const auto *sli_630 = buffer.data(sli + 630);
    const auto *sli_637 = buffer.data(sli + 637);
    const auto *sli_639 = buffer.data(sli + 639);
    const auto *sli_640 = buffer.data(sli + 640);
    const auto *sli_641 = buffer.data(sli + 641);
    const auto *sli_642 = buffer.data(sli + 642);
    const auto *sli_643 = buffer.data(sli + 643);
    const auto *sli_644 = buffer.data(sli + 644);
    const auto *sli_646 = buffer.data(sli + 646);
    const auto *sli_649 = buffer.data(sli + 649);
    const auto *sli_653 = buffer.data(sli + 653);
    const auto *sli_658 = buffer.data(sli + 658);
    const auto *sli_665 = buffer.data(sli + 665);
    const auto *sli_667 = buffer.data(sli + 667);
    const auto *sli_668 = buffer.data(sli + 668);
    const auto *sli_669 = buffer.data(sli + 669);
    const auto *sli_670 = buffer.data(sli + 670);
    const auto *sli_671 = buffer.data(sli + 671);
    const auto *sli_784 = buffer.data(sli + 784);
    const auto *sli_787 = buffer.data(sli + 787);
    const auto *sli_789 = buffer.data(sli + 789);
    const auto *sli_790 = buffer.data(sli + 790);
    const auto *sli_793 = buffer.data(sli + 793);
    const auto *sli_794 = buffer.data(sli + 794);
    const auto *sli_796 = buffer.data(sli + 796);
    const auto *sli_798 = buffer.data(sli + 798);
    const auto *sli_799 = buffer.data(sli + 799);
    const auto *sli_801 = buffer.data(sli + 801);
    const auto *sli_802 = buffer.data(sli + 802);
    const auto *sli_804 = buffer.data(sli + 804);
    const auto *sli_805 = buffer.data(sli + 805);
    const auto *sli_806 = buffer.data(sli + 806);
    const auto *sli_807 = buffer.data(sli + 807);
    const auto *sli_808 = buffer.data(sli + 808);
    const auto *sli_809 = buffer.data(sli + 809);
    const auto *sli_810 = buffer.data(sli + 810);
    const auto *sli_811 = buffer.data(sli + 811);
    const auto *sli_817 = buffer.data(sli + 817);
    const auto *sli_821 = buffer.data(sli + 821);
    const auto *sli_826 = buffer.data(sli + 826);
    const auto *sli_832 = buffer.data(sli + 832);
    const auto *sli_833 = buffer.data(sli + 833);
    const auto *sli_834 = buffer.data(sli + 834);
    const auto *sli_835 = buffer.data(sli + 835);
    const auto *sli_836 = buffer.data(sli + 836);
    const auto *sli_837 = buffer.data(sli + 837);
    const auto *sli_838 = buffer.data(sli + 838);
    const auto *sli_839 = buffer.data(sli + 839);
    const auto *sli_840 = buffer.data(sli + 840);
    const auto *sli_843 = buffer.data(sli + 843);
    const auto *sli_845 = buffer.data(sli + 845);
    const auto *sli_846 = buffer.data(sli + 846);
    const auto *sli_849 = buffer.data(sli + 849);
    const auto *sli_850 = buffer.data(sli + 850);
    const auto *sli_852 = buffer.data(sli + 852);
    const auto *sli_854 = buffer.data(sli + 854);
    const auto *sli_855 = buffer.data(sli + 855);
    const auto *sli_857 = buffer.data(sli + 857);
    const auto *sli_858 = buffer.data(sli + 858);
    const auto *sli_860 = buffer.data(sli + 860);
    const auto *sli_861 = buffer.data(sli + 861);
    const auto *sli_862 = buffer.data(sli + 862);
    const auto *sli_863 = buffer.data(sli + 863);
    const auto *sli_864 = buffer.data(sli + 864);
    const auto *sli_865 = buffer.data(sli + 865);
    const auto *sli_866 = buffer.data(sli + 866);
    const auto *sli_867 = buffer.data(sli + 867);
    const auto *sli_868 = buffer.data(sli + 868);

    const auto *slk1_756 = buffer.data(slk1 + 756);
    const auto *slk1_759 = buffer.data(slk1 + 759);
    const auto *slk1_762 = buffer.data(slk1 + 762);
    const auto *slk1_766 = buffer.data(slk1 + 766);
    const auto *slk1_768 = buffer.data(slk1 + 768);
    const auto *slk1_771 = buffer.data(slk1 + 771);
    const auto *slk1_773 = buffer.data(slk1 + 773);
    const auto *slk1_774 = buffer.data(slk1 + 774);
    const auto *slk1_784 = buffer.data(slk1 + 784);

    const auto *smh0_588 = buffer.data(smh0 + 588);
    const auto *smh0_591 = buffer.data(smh0 + 591);
    const auto *smh0_593 = buffer.data(smh0 + 593);
    const auto *smh0_594 = buffer.data(smh0 + 594);
    const auto *smh0_597 = buffer.data(smh0 + 597);
    const auto *smh0_598 = buffer.data(smh0 + 598);
    const auto *smh0_600 = buffer.data(smh0 + 600);
    const auto *smh0_602 = buffer.data(smh0 + 602);
    const auto *smh0_603 = buffer.data(smh0 + 603);
    const auto *smh0_605 = buffer.data(smh0 + 605);
    const auto *smh0_606 = buffer.data(smh0 + 606);
    const auto *smh0_607 = buffer.data(smh0 + 607);
    const auto *smh0_608 = buffer.data(smh0 + 608);
    const auto *smh0_614 = buffer.data(smh0 + 614);
    const auto *smh0_618 = buffer.data(smh0 + 618);
    const auto *smh0_623 = buffer.data(smh0 + 623);
    const auto *smh0_626 = buffer.data(smh0 + 626);
    const auto *smh0_627 = buffer.data(smh0 + 627);
    const auto *smh0_628 = buffer.data(smh0 + 628);
    const auto *smh0_629 = buffer.data(smh0 + 629);
    const auto *smh0_630 = buffer.data(smh0 + 630);
    const auto *smh0_633 = buffer.data(smh0 + 633);
    const auto *smh0_635 = buffer.data(smh0 + 635);
    const auto *smh0_636 = buffer.data(smh0 + 636);
    const auto *smh0_639 = buffer.data(smh0 + 639);
    const auto *smh0_640 = buffer.data(smh0 + 640);
    const auto *smh0_642 = buffer.data(smh0 + 642);
    const auto *smh0_644 = buffer.data(smh0 + 644);
    const auto *smh0_645 = buffer.data(smh0 + 645);
    const auto *smh0_647 = buffer.data(smh0 + 647);
    const auto *smh0_648 = buffer.data(smh0 + 648);
    const auto *smh0_649 = buffer.data(smh0 + 649);
    const auto *smh0_650 = buffer.data(smh0 + 650);
    const auto *smh0_651 = buffer.data(smh0 + 651);

    const auto *smh1_588 = buffer.data(smh1 + 588);
    const auto *smh1_591 = buffer.data(smh1 + 591);
    const auto *smh1_593 = buffer.data(smh1 + 593);
    const auto *smh1_594 = buffer.data(smh1 + 594);
    const auto *smh1_597 = buffer.data(smh1 + 597);
    const auto *smh1_598 = buffer.data(smh1 + 598);
    const auto *smh1_600 = buffer.data(smh1 + 600);
    const auto *smh1_602 = buffer.data(smh1 + 602);
    const auto *smh1_603 = buffer.data(smh1 + 603);
    const auto *smh1_605 = buffer.data(smh1 + 605);
    const auto *smh1_606 = buffer.data(smh1 + 606);
    const auto *smh1_607 = buffer.data(smh1 + 607);
    const auto *smh1_608 = buffer.data(smh1 + 608);
    const auto *smh1_614 = buffer.data(smh1 + 614);
    const auto *smh1_618 = buffer.data(smh1 + 618);
    const auto *smh1_623 = buffer.data(smh1 + 623);
    const auto *smh1_626 = buffer.data(smh1 + 626);
    const auto *smh1_627 = buffer.data(smh1 + 627);
    const auto *smh1_628 = buffer.data(smh1 + 628);
    const auto *smh1_629 = buffer.data(smh1 + 629);
    const auto *smh1_630 = buffer.data(smh1 + 630);
    const auto *smh1_633 = buffer.data(smh1 + 633);
    const auto *smh1_635 = buffer.data(smh1 + 635);
    const auto *smh1_636 = buffer.data(smh1 + 636);
    const auto *smh1_639 = buffer.data(smh1 + 639);
    const auto *smh1_640 = buffer.data(smh1 + 640);
    const auto *smh1_642 = buffer.data(smh1 + 642);
    const auto *smh1_644 = buffer.data(smh1 + 644);
    const auto *smh1_645 = buffer.data(smh1 + 645);
    const auto *smh1_647 = buffer.data(smh1 + 647);
    const auto *smh1_648 = buffer.data(smh1 + 648);
    const auto *smh1_649 = buffer.data(smh1 + 649);
    const auto *smh1_650 = buffer.data(smh1 + 650);
    const auto *smh1_651 = buffer.data(smh1 + 651);

    const auto *smi_784 = buffer.data(smi + 784);
    const auto *smi_786 = buffer.data(smi + 786);
    const auto *smi_787 = buffer.data(smi + 787);
    const auto *smi_789 = buffer.data(smi + 789);
    const auto *smi_790 = buffer.data(smi + 790);
    const auto *smi_793 = buffer.data(smi + 793);
    const auto *smi_794 = buffer.data(smi + 794);
    const auto *smi_796 = buffer.data(smi + 796);
    const auto *smi_798 = buffer.data(smi + 798);
    const auto *smi_799 = buffer.data(smi + 799);
    const auto *smi_801 = buffer.data(smi + 801);
    const auto *smi_802 = buffer.data(smi + 802);
    const auto *smi_804 = buffer.data(smi + 804);
    const auto *smi_805 = buffer.data(smi + 805);
    const auto *smi_806 = buffer.data(smi + 806);
    const auto *smi_807 = buffer.data(smi + 807);
    const auto *smi_808 = buffer.data(smi + 808);
    const auto *smi_809 = buffer.data(smi + 809);
    const auto *smi_810 = buffer.data(smi + 810);
    const auto *smi_811 = buffer.data(smi + 811);
    const auto *smi_812 = buffer.data(smi + 812);
    const auto *smi_814 = buffer.data(smi + 814);
    const auto *smi_815 = buffer.data(smi + 815);
    const auto *smi_817 = buffer.data(smi + 817);
    const auto *smi_818 = buffer.data(smi + 818);
    const auto *smi_821 = buffer.data(smi + 821);
    const auto *smi_822 = buffer.data(smi + 822);
    const auto *smi_826 = buffer.data(smi + 826);
    const auto *smi_832 = buffer.data(smi + 832);
    const auto *smi_833 = buffer.data(smi + 833);
    const auto *smi_834 = buffer.data(smi + 834);
    const auto *smi_835 = buffer.data(smi + 835);
    const auto *smi_836 = buffer.data(smi + 836);
    const auto *smi_837 = buffer.data(smi + 837);
    const auto *smi_838 = buffer.data(smi + 838);
    const auto *smi_839 = buffer.data(smi + 839);
    const auto *smi_840 = buffer.data(smi + 840);
    const auto *smi_842 = buffer.data(smi + 842);
    const auto *smi_843 = buffer.data(smi + 843);
    const auto *smi_845 = buffer.data(smi + 845);
    const auto *smi_846 = buffer.data(smi + 846);
    const auto *smi_849 = buffer.data(smi + 849);
    const auto *smi_850 = buffer.data(smi + 850);
    const auto *smi_852 = buffer.data(smi + 852);
    const auto *smi_854 = buffer.data(smi + 854);
    const auto *smi_855 = buffer.data(smi + 855);
    const auto *smi_857 = buffer.data(smi + 857);
    const auto *smi_858 = buffer.data(smi + 858);
    const auto *smi_860 = buffer.data(smi + 860);
    const auto *smi_861 = buffer.data(smi + 861);
    const auto *smi_862 = buffer.data(smi + 862);
    const auto *smi_863 = buffer.data(smi + 863);
    const auto *smi_864 = buffer.data(smi + 864);
    const auto *smi_865 = buffer.data(smi + 865);
    const auto *smi_866 = buffer.data(smi + 866);
    const auto *smi_867 = buffer.data(smi + 867);
    const auto *smi_868 = buffer.data(smi + 868);

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, sli_588, sli_784, \
                         sli_787, smh0_588, smh0_591, smh1_588, smh1_591, smi_784, \
                         smi_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_14 * sli_784[k]
                    + f_1 * smh0_588[k]
                    - f_2 * smh1_588[k]
                    + f_3 * pc_x[k] * smi_784[k];

        t_1009[k] = f_19 * sli_588[k]
                    + f_3 * pc_y[k] * smi_784[k];

        t_1010[k] = f_3 * pc_z[k] * smi_784[k];

        t_1011[k] = f_14 * sli_787[k]
                    + f_4 * smh0_591[k]
                    - f_5 * smh1_591[k]
                    + f_3 * pc_x[k] * smi_787[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pc_x, pc_y, sli_590, sli_789, sli_790, \
                         smh0_593, smh0_594, smh1_593, smh1_594, smi_786, smi_789, \
                         smi_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_19 * sli_590[k]
                    + f_3 * pc_y[k] * smi_786[k];

        t_1013[k] = f_14 * sli_789[k]
                    + f_4 * smh0_593[k]
                    - f_5 * smh1_593[k]
                    + f_3 * pc_x[k] * smi_789[k];

        t_1014[k] = f_14 * sli_790[k]
                    + f_6 * smh0_594[k]
                    - f_7 * smh1_594[k]
                    + f_3 * pc_x[k] * smi_790[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, pc_x, pc_y, pc_z, sli_593, sli_793, smh0_597, \
                         smh1_597, smi_787, smi_789, smi_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = f_3 * pc_z[k] * smi_787[k];

        t_1016[k] = f_19 * sli_593[k]
                    + f_3 * pc_y[k] * smi_789[k];

        t_1017[k] = f_14 * sli_793[k]
                    + f_6 * smh0_597[k]
                    - f_7 * smh1_597[k]
                    + f_3 * pc_x[k] * smi_793[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, pc_x, pc_z, sli_794, sli_796, smh0_598, \
                         smh0_600, smh1_598, smh1_600, smi_790, smi_794, \
                         smi_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_14 * sli_794[k]
                    + f_8 * smh0_598[k]
                    - f_9 * smh1_598[k]
                    + f_3 * pc_x[k] * smi_794[k];

        t_1019[k] = f_3 * pc_z[k] * smi_790[k];

        t_1020[k] = f_14 * sli_796[k]
                    + f_8 * smh0_600[k]
                    - f_9 * smh1_600[k]
                    + f_3 * pc_x[k] * smi_796[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, pc_x, pc_y, sli_597, sli_798, sli_799, \
                         smh0_602, smh0_603, smh1_602, smh1_603, smi_793, smi_798, \
                         smi_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_19 * sli_597[k]
                    + f_3 * pc_y[k] * smi_793[k];

        t_1022[k] = f_14 * sli_798[k]
                    + f_8 * smh0_602[k]
                    - f_9 * smh1_602[k]
                    + f_3 * pc_x[k] * smi_798[k];

        t_1023[k] = f_14 * sli_799[k]
                    + f_10 * smh0_603[k]
                    - f_11 * smh1_603[k]
                    + f_3 * pc_x[k] * smi_799[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, pc_x, pc_z, sli_801, sli_802, smh0_605, \
                         smh0_606, smh1_605, smh1_606, smi_794, smi_801, \
                         smi_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_3 * pc_z[k] * smi_794[k];

        t_1025[k] = f_14 * sli_801[k]
                    + f_10 * smh0_605[k]
                    - f_11 * smh1_605[k]
                    + f_3 * pc_x[k] * smi_801[k];

        t_1026[k] = f_14 * sli_802[k]
                    + f_10 * smh0_606[k]
                    - f_11 * smh1_606[k]
                    + f_3 * pc_x[k] * smi_802[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, pc_x, pc_y, sli_602, sli_804, \
                         sli_805, sli_806, smh0_608, smh1_608, smi_798, smi_804, smi_805, \
                         smi_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_19 * sli_602[k]
                    + f_3 * pc_y[k] * smi_798[k];

        t_1028[k] = f_14 * sli_804[k]
                    + f_10 * smh0_608[k]
                    - f_11 * smh1_608[k]
                    + f_3 * pc_x[k] * smi_804[k];

        t_1029[k] = f_14 * sli_805[k]
                    + f_3 * pc_x[k] * smi_805[k];

        t_1030[k] = f_14 * sli_806[k]
                    + f_3 * pc_x[k] * smi_806[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, t_1034, t_1035, pc_x, sli_807, sli_808, \
                         sli_809, sli_810, sli_811, smi_807, smi_808, smi_809, smi_810, \
                         smi_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_14 * sli_807[k]
                    + f_3 * pc_x[k] * smi_807[k];

        t_1032[k] = f_14 * sli_808[k]
                    + f_3 * pc_x[k] * smi_808[k];

        t_1033[k] = f_14 * sli_809[k]
                    + f_3 * pc_x[k] * smi_809[k];

        t_1034[k] = f_14 * sli_810[k]
                    + f_3 * pc_x[k] * smi_810[k];

        t_1035[k] = f_14 * sli_811[k]
                    + f_3 * pc_x[k] * smi_811[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, pc_y, pc_z, sli_609, sli_611, smh0_603, \
                         smh0_605, smh1_603, smh1_605, smi_805, \
                         smi_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = f_19 * sli_609[k]
                    + f_1 * smh0_603[k]
                    - f_2 * smh1_603[k]
                    + f_3 * pc_y[k] * smi_805[k];

        t_1037[k] = f_3 * pc_z[k] * smi_805[k];

        t_1038[k] = f_19 * sli_611[k]
                    + f_4 * smh0_605[k]
                    - f_5 * smh1_605[k]
                    + f_3 * pc_y[k] * smi_807[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pc_y, sli_612, sli_613, sli_614, smh0_606, \
                         smh0_607, smh0_608, smh1_606, smh1_607, smh1_608, smi_808, smi_809, \
                         smi_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_19 * sli_612[k]
                    + f_6 * smh0_606[k]
                    - f_7 * smh1_606[k]
                    + f_3 * pc_y[k] * smi_808[k];

        t_1040[k] = f_19 * sli_613[k]
                    + f_8 * smh0_607[k]
                    - f_9 * smh1_607[k]
                    + f_3 * pc_y[k] * smi_809[k];

        t_1041[k] = f_19 * sli_614[k]
                    + f_10 * smh0_608[k]
                    - f_11 * smh1_608[k]
                    + f_3 * pc_y[k] * smi_810[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pb_z, pc_y, pc_z, slk0_756, sli_615, \
                         sli_616, slk1_756, smh0_608, smh1_608, smi_811, \
                         smi_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_19 * sli_615[k]
                    + f_3 * pc_y[k] * smi_811[k];

        t_1043[k] = f_1 * smh0_608[k]
                    - f_2 * smh1_608[k]
                    + f_3 * pc_z[k] * smi_811[k];

        t_1044[k] = pb_z[k] * slk0_756[k]
                    - f_12 * pc_z[k] * slk1_756[k];

        t_1045[k] = f_20 * sli_616[k]
                    + f_3 * pc_y[k] * smi_812[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pb_z, pc_y, pc_z, slk0_759, sli_588, sli_618, \
                         slk1_759, smi_812, smi_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_13 * sli_588[k]
                    + f_3 * pc_z[k] * smi_812[k];

        t_1047[k] = pb_z[k] * slk0_759[k]
                    - f_12 * pc_z[k] * slk1_759[k];

        t_1048[k] = f_20 * sli_618[k]
                    + f_3 * pc_y[k] * smi_814[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pb_z, pc_x, pc_z, slk0_762, sli_591, sli_817, \
                         slk1_762, smh0_614, smh1_614, smi_815, \
                         smi_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_14 * sli_817[k]
                    + f_4 * smh0_614[k]
                    - f_5 * smh1_614[k]
                    + f_3 * pc_x[k] * smi_817[k];

        t_1050[k] = pb_z[k] * slk0_762[k]
                    - f_12 * pc_z[k] * slk1_762[k];

        t_1051[k] = f_13 * sli_591[k]
                    + f_3 * pc_z[k] * smi_815[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pb_z, pc_x, pc_y, pc_z, slk0_766, sli_621, \
                         sli_821, slk1_766, smh0_618, smh1_618, smi_817, \
                         smi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_20 * sli_621[k]
                    + f_3 * pc_y[k] * smi_817[k];

        t_1053[k] = f_14 * sli_821[k]
                    + f_6 * smh0_618[k]
                    - f_7 * smh1_618[k]
                    + f_3 * pc_x[k] * smi_821[k];

        t_1054[k] = pb_z[k] * slk0_766[k]
                    - f_12 * pc_z[k] * slk1_766[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pb_z, pc_y, pc_z, slk0_768, sli_594, sli_595, \
                         sli_625, slk1_768, smi_818, smi_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_13 * sli_594[k]
                    + f_3 * pc_z[k] * smi_818[k];

        t_1056[k] = pb_z[k] * slk0_768[k]
                    + f_14 * sli_595[k]
                    - f_12 * pc_z[k] * slk1_768[k];

        t_1057[k] = f_20 * sli_625[k]
                    + f_3 * pc_y[k] * smi_821[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pb_z, pc_x, pc_z, slk0_771, sli_598, sli_826, \
                         slk1_771, smh0_623, smh1_623, smi_822, \
                         smi_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_14 * sli_826[k]
                    + f_8 * smh0_623[k]
                    - f_9 * smh1_623[k]
                    + f_3 * pc_x[k] * smi_826[k];

        t_1059[k] = pb_z[k] * slk0_771[k]
                    - f_12 * pc_z[k] * slk1_771[k];

        t_1060[k] = f_13 * sli_598[k]
                    + f_3 * pc_z[k] * smi_822[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pb_z, pc_y, pc_z, slk0_773, slk0_774, \
                         sli_599, sli_600, sli_630, slk1_773, slk1_774, \
                         smi_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = pb_z[k] * slk0_773[k]
                    + f_14 * sli_599[k]
                    - f_12 * pc_z[k] * slk1_773[k];

        t_1062[k] = pb_z[k] * slk0_774[k]
                    + f_15 * sli_600[k]
                    - f_12 * pc_z[k] * slk1_774[k];

        t_1063[k] = f_20 * sli_630[k]
                    + f_3 * pc_y[k] * smi_826[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, t_1067, pc_x, sli_832, sli_833, sli_834, \
                         sli_835, smh0_629, smh1_629, smi_832, smi_833, smi_834, \
                         smi_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_14 * sli_832[k]
                    + f_10 * smh0_629[k]
                    - f_11 * smh1_629[k]
                    + f_3 * pc_x[k] * smi_832[k];

        t_1065[k] = f_14 * sli_833[k]
                    + f_3 * pc_x[k] * smi_833[k];

        t_1066[k] = f_14 * sli_834[k]
                    + f_3 * pc_x[k] * smi_834[k];

        t_1067[k] = f_14 * sli_835[k]
                    + f_3 * pc_x[k] * smi_835[k];
    }

#pragma omp simd aligned(t_1068, t_1069, t_1070, t_1071, pc_x, sli_836, sli_837, sli_838, \
                         sli_839, smi_836, smi_837, smi_838, smi_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = f_14 * sli_836[k]
                    + f_3 * pc_x[k] * smi_836[k];

        t_1069[k] = f_14 * sli_837[k]
                    + f_3 * pc_x[k] * smi_837[k];

        t_1070[k] = f_14 * sli_838[k]
                    + f_3 * pc_x[k] * smi_838[k];

        t_1071[k] = f_14 * sli_839[k]
                    + f_3 * pc_x[k] * smi_839[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, pb_z, pc_y, pc_z, slk0_784, sli_609, sli_639, \
                         slk1_784, smh0_626, smh1_626, smi_833, \
                         smi_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = pb_z[k] * slk0_784[k]
                    - f_12 * pc_z[k] * slk1_784[k];

        t_1073[k] = f_13 * sli_609[k]
                    + f_3 * pc_z[k] * smi_833[k];

        t_1074[k] = f_20 * sli_639[k]
                    + f_4 * smh0_626[k]
                    - f_5 * smh1_626[k]
                    + f_3 * pc_y[k] * smi_835[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, pc_y, sli_640, sli_641, sli_642, smh0_627, \
                         smh0_628, smh0_629, smh1_627, smh1_628, smh1_629, smi_836, smi_837, \
                         smi_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_20 * sli_640[k]
                    + f_6 * smh0_627[k]
                    - f_7 * smh1_627[k]
                    + f_3 * pc_y[k] * smi_836[k];

        t_1076[k] = f_20 * sli_641[k]
                    + f_8 * smh0_628[k]
                    - f_9 * smh1_628[k]
                    + f_3 * pc_y[k] * smi_837[k];

        t_1077[k] = f_20 * sli_642[k]
                    + f_10 * smh0_629[k]
                    - f_11 * smh1_629[k]
                    + f_3 * pc_y[k] * smi_838[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, pc_x, pc_y, pc_z, sli_615, sli_643, sli_840, \
                         smh0_629, smh0_630, smh1_629, smh1_630, smi_839, \
                         smi_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = f_20 * sli_643[k]
                    + f_3 * pc_y[k] * smi_839[k];

        t_1079[k] = f_13 * sli_615[k]
                    + f_1 * smh0_629[k]
                    - f_2 * smh1_629[k]
                    + f_3 * pc_z[k] * smi_839[k];

        t_1080[k] = f_14 * sli_840[k]
                    + f_1 * smh0_630[k]
                    - f_2 * smh1_630[k]
                    + f_3 * pc_x[k] * smi_840[k];
    }

#pragma omp simd aligned(t_1081, t_1082, t_1083, t_1084, pc_x, pc_y, pc_z, sli_616, sli_644, \
                         sli_646, sli_843, smh0_633, smh1_633, smi_840, smi_842, \
                         smi_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1081[k] = f_17 * sli_644[k]
                    + f_3 * pc_y[k] * smi_840[k];

        t_1082[k] = f_14 * sli_616[k]
                    + f_3 * pc_z[k] * smi_840[k];

        t_1083[k] = f_14 * sli_843[k]
                    + f_4 * smh0_633[k]
                    - f_5 * smh1_633[k]
                    + f_3 * pc_x[k] * smi_843[k];

        t_1084[k] = f_17 * sli_646[k]
                    + f_3 * pc_y[k] * smi_842[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_x, pc_z, sli_619, sli_845, sli_846, \
                         smh0_635, smh0_636, smh1_635, smh1_636, smi_843, smi_845, \
                         smi_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_14 * sli_845[k]
                    + f_4 * smh0_635[k]
                    - f_5 * smh1_635[k]
                    + f_3 * pc_x[k] * smi_845[k];

        t_1086[k] = f_14 * sli_846[k]
                    + f_6 * smh0_636[k]
                    - f_7 * smh1_636[k]
                    + f_3 * pc_x[k] * smi_846[k];

        t_1087[k] = f_14 * sli_619[k]
                    + f_3 * pc_z[k] * smi_843[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_x, pc_y, sli_649, sli_849, sli_850, \
                         smh0_639, smh0_640, smh1_639, smh1_640, smi_845, smi_849, \
                         smi_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * sli_649[k]
                    + f_3 * pc_y[k] * smi_845[k];

        t_1089[k] = f_14 * sli_849[k]
                    + f_6 * smh0_639[k]
                    - f_7 * smh1_639[k]
                    + f_3 * pc_x[k] * smi_849[k];

        t_1090[k] = f_14 * sli_850[k]
                    + f_8 * smh0_640[k]
                    - f_9 * smh1_640[k]
                    + f_3 * pc_x[k] * smi_850[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_y, pc_z, sli_622, sli_653, sli_852, \
                         smh0_642, smh1_642, smi_846, smi_849, \
                         smi_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_14 * sli_622[k]
                    + f_3 * pc_z[k] * smi_846[k];

        t_1092[k] = f_14 * sli_852[k]
                    + f_8 * smh0_642[k]
                    - f_9 * smh1_642[k]
                    + f_3 * pc_x[k] * smi_852[k];

        t_1093[k] = f_17 * sli_653[k]
                    + f_3 * pc_y[k] * smi_849[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, pc_z, sli_626, sli_854, sli_855, \
                         smh0_644, smh0_645, smh1_644, smh1_645, smi_850, smi_854, \
                         smi_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_14 * sli_854[k]
                    + f_8 * smh0_644[k]
                    - f_9 * smh1_644[k]
                    + f_3 * pc_x[k] * smi_854[k];

        t_1095[k] = f_14 * sli_855[k]
                    + f_10 * smh0_645[k]
                    - f_11 * smh1_645[k]
                    + f_3 * pc_x[k] * smi_855[k];

        t_1096[k] = f_14 * sli_626[k]
                    + f_3 * pc_z[k] * smi_850[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_y, sli_658, sli_857, sli_858, \
                         smh0_647, smh0_648, smh1_647, smh1_648, smi_854, smi_857, \
                         smi_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_14 * sli_857[k]
                    + f_10 * smh0_647[k]
                    - f_11 * smh1_647[k]
                    + f_3 * pc_x[k] * smi_857[k];

        t_1098[k] = f_14 * sli_858[k]
                    + f_10 * smh0_648[k]
                    - f_11 * smh1_648[k]
                    + f_3 * pc_x[k] * smi_858[k];

        t_1099[k] = f_17 * sli_658[k]
                    + f_3 * pc_y[k] * smi_854[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pc_x, sli_860, sli_861, sli_862, \
                         sli_863, smh0_650, smh1_650, smi_860, smi_861, smi_862, \
                         smi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_14 * sli_860[k]
                    + f_10 * smh0_650[k]
                    - f_11 * smh1_650[k]
                    + f_3 * pc_x[k] * smi_860[k];

        t_1101[k] = f_14 * sli_861[k]
                    + f_3 * pc_x[k] * smi_861[k];

        t_1102[k] = f_14 * sli_862[k]
                    + f_3 * pc_x[k] * smi_862[k];

        t_1103[k] = f_14 * sli_863[k]
                    + f_3 * pc_x[k] * smi_863[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, sli_864, sli_865, sli_866, \
                         sli_867, smi_864, smi_865, smi_866, smi_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_14 * sli_864[k]
                    + f_3 * pc_x[k] * smi_864[k];

        t_1105[k] = f_14 * sli_865[k]
                    + f_3 * pc_x[k] * smi_865[k];

        t_1106[k] = f_14 * sli_866[k]
                    + f_3 * pc_x[k] * smi_866[k];

        t_1107[k] = f_14 * sli_867[k]
                    + f_3 * pc_x[k] * smi_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, sli_637, sli_665, sli_667, \
                         smh0_645, smh0_647, smh1_645, smh1_647, smi_861, \
                         smi_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_17 * sli_665[k]
                    + f_1 * smh0_645[k]
                    - f_2 * smh1_645[k]
                    + f_3 * pc_y[k] * smi_861[k];

        t_1109[k] = f_14 * sli_637[k]
                    + f_3 * pc_z[k] * smi_861[k];

        t_1110[k] = f_17 * sli_667[k]
                    + f_4 * smh0_647[k]
                    - f_5 * smh1_647[k]
                    + f_3 * pc_y[k] * smi_863[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, pc_y, sli_668, sli_669, sli_670, smh0_648, \
                         smh0_649, smh0_650, smh1_648, smh1_649, smh1_650, smi_864, smi_865, \
                         smi_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_17 * sli_668[k]
                    + f_6 * smh0_648[k]
                    - f_7 * smh1_648[k]
                    + f_3 * pc_y[k] * smi_864[k];

        t_1112[k] = f_17 * sli_669[k]
                    + f_8 * smh0_649[k]
                    - f_9 * smh1_649[k]
                    + f_3 * pc_y[k] * smi_865[k];

        t_1113[k] = f_17 * sli_670[k]
                    + f_10 * smh0_650[k]
                    - f_11 * smh1_650[k]
                    + f_3 * pc_y[k] * smi_866[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_x, pc_y, pc_z, sli_643, sli_671, sli_868, \
                         smh0_650, smh0_651, smh1_650, smh1_651, smi_867, \
                         smi_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_17 * sli_671[k]
                    + f_3 * pc_y[k] * smi_867[k];

        t_1115[k] = f_14 * sli_643[k]
                    + f_1 * smh0_650[k]
                    - f_2 * smh1_650[k]
                    + f_3 * pc_z[k] * smi_867[k];

        t_1116[k] = f_14 * sli_868[k]
                    + f_1 * smh0_651[k]
                    - f_2 * smh1_651[k]
                    + f_3 * pc_x[k] * smi_868[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t sli, const size_t smh0,
                                                           const size_t smh1, const size_t smi,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli_644 = buffer.data(sli + 644);
    const auto *sli_647 = buffer.data(sli + 647);
    const auto *sli_650 = buffer.data(sli + 650);
    const auto *sli_654 = buffer.data(sli + 654);
    const auto *sli_665 = buffer.data(sli + 665);
    const auto *sli_671 = buffer.data(sli + 671);
    const auto *sli_672 = buffer.data(sli + 672);
    const auto *sli_674 = buffer.data(sli + 674);
    const auto *sli_675 = buffer.data(sli + 675);
    const auto *sli_677 = buffer.data(sli + 677);
    const auto *sli_678 = buffer.data(sli + 678);
    const auto *sli_681 = buffer.data(sli + 681);
    const auto *sli_682 = buffer.data(sli + 682);
    const auto *sli_686 = buffer.data(sli + 686);
    const auto *sli_693 = buffer.data(sli + 693);
    const auto *sli_695 = buffer.data(sli + 695);
    const auto *sli_696 = buffer.data(sli + 696);
    const auto *sli_697 = buffer.data(sli + 697);
    const auto *sli_698 = buffer.data(sli + 698);
    const auto *sli_699 = buffer.data(sli + 699);
    const auto *sli_700 = buffer.data(sli + 700);
    const auto *sli_702 = buffer.data(sli + 702);
    const auto *sli_703 = buffer.data(sli + 703);
    const auto *sli_705 = buffer.data(sli + 705);
    const auto *sli_706 = buffer.data(sli + 706);
    const auto *sli_709 = buffer.data(sli + 709);
    const auto *sli_710 = buffer.data(sli + 710);
    const auto *sli_714 = buffer.data(sli + 714);
    const auto *sli_721 = buffer.data(sli + 721);
    const auto *sli_723 = buffer.data(sli + 723);
    const auto *sli_724 = buffer.data(sli + 724);
    const auto *sli_725 = buffer.data(sli + 725);
    const auto *sli_726 = buffer.data(sli + 726);
    const auto *sli_727 = buffer.data(sli + 727);
    const auto *sli_728 = buffer.data(sli + 728);
    const auto *sli_730 = buffer.data(sli + 730);
    const auto *sli_733 = buffer.data(sli + 733);
    const auto *sli_737 = buffer.data(sli + 737);
    const auto *sli_742 = buffer.data(sli + 742);
    const auto *sli_749 = buffer.data(sli + 749);
    const auto *sli_751 = buffer.data(sli + 751);
    const auto *sli_752 = buffer.data(sli + 752);
    const auto *sli_753 = buffer.data(sli + 753);
    const auto *sli_754 = buffer.data(sli + 754);
    const auto *sli_871 = buffer.data(sli + 871);
    const auto *sli_873 = buffer.data(sli + 873);
    const auto *sli_874 = buffer.data(sli + 874);
    const auto *sli_877 = buffer.data(sli + 877);
    const auto *sli_878 = buffer.data(sli + 878);
    const auto *sli_880 = buffer.data(sli + 880);
    const auto *sli_882 = buffer.data(sli + 882);
    const auto *sli_883 = buffer.data(sli + 883);
    const auto *sli_885 = buffer.data(sli + 885);
    const auto *sli_886 = buffer.data(sli + 886);
    const auto *sli_888 = buffer.data(sli + 888);
    const auto *sli_889 = buffer.data(sli + 889);
    const auto *sli_890 = buffer.data(sli + 890);
    const auto *sli_891 = buffer.data(sli + 891);
    const auto *sli_892 = buffer.data(sli + 892);
    const auto *sli_893 = buffer.data(sli + 893);
    const auto *sli_894 = buffer.data(sli + 894);
    const auto *sli_895 = buffer.data(sli + 895);
    const auto *sli_896 = buffer.data(sli + 896);
    const auto *sli_899 = buffer.data(sli + 899);
    const auto *sli_901 = buffer.data(sli + 901);
    const auto *sli_902 = buffer.data(sli + 902);
    const auto *sli_905 = buffer.data(sli + 905);
    const auto *sli_906 = buffer.data(sli + 906);
    const auto *sli_908 = buffer.data(sli + 908);
    const auto *sli_910 = buffer.data(sli + 910);
    const auto *sli_911 = buffer.data(sli + 911);
    const auto *sli_913 = buffer.data(sli + 913);
    const auto *sli_914 = buffer.data(sli + 914);
    const auto *sli_916 = buffer.data(sli + 916);
    const auto *sli_917 = buffer.data(sli + 917);
    const auto *sli_918 = buffer.data(sli + 918);
    const auto *sli_919 = buffer.data(sli + 919);
    const auto *sli_920 = buffer.data(sli + 920);
    const auto *sli_921 = buffer.data(sli + 921);
    const auto *sli_922 = buffer.data(sli + 922);
    const auto *sli_923 = buffer.data(sli + 923);
    const auto *sli_924 = buffer.data(sli + 924);
    const auto *sli_927 = buffer.data(sli + 927);
    const auto *sli_929 = buffer.data(sli + 929);
    const auto *sli_930 = buffer.data(sli + 930);
    const auto *sli_933 = buffer.data(sli + 933);
    const auto *sli_934 = buffer.data(sli + 934);
    const auto *sli_936 = buffer.data(sli + 936);
    const auto *sli_938 = buffer.data(sli + 938);
    const auto *sli_939 = buffer.data(sli + 939);
    const auto *sli_941 = buffer.data(sli + 941);
    const auto *sli_942 = buffer.data(sli + 942);
    const auto *sli_944 = buffer.data(sli + 944);
    const auto *sli_945 = buffer.data(sli + 945);
    const auto *sli_946 = buffer.data(sli + 946);
    const auto *sli_947 = buffer.data(sli + 947);
    const auto *sli_948 = buffer.data(sli + 948);
    const auto *sli_949 = buffer.data(sli + 949);
    const auto *sli_950 = buffer.data(sli + 950);
    const auto *sli_951 = buffer.data(sli + 951);

    const auto *smh0_654 = buffer.data(smh0 + 654);
    const auto *smh0_656 = buffer.data(smh0 + 656);
    const auto *smh0_657 = buffer.data(smh0 + 657);
    const auto *smh0_660 = buffer.data(smh0 + 660);
    const auto *smh0_661 = buffer.data(smh0 + 661);
    const auto *smh0_663 = buffer.data(smh0 + 663);
    const auto *smh0_665 = buffer.data(smh0 + 665);
    const auto *smh0_666 = buffer.data(smh0 + 666);
    const auto *smh0_668 = buffer.data(smh0 + 668);
    const auto *smh0_669 = buffer.data(smh0 + 669);
    const auto *smh0_670 = buffer.data(smh0 + 670);
    const auto *smh0_671 = buffer.data(smh0 + 671);
    const auto *smh0_672 = buffer.data(smh0 + 672);
    const auto *smh0_675 = buffer.data(smh0 + 675);
    const auto *smh0_677 = buffer.data(smh0 + 677);
    const auto *smh0_678 = buffer.data(smh0 + 678);
    const auto *smh0_681 = buffer.data(smh0 + 681);
    const auto *smh0_682 = buffer.data(smh0 + 682);
    const auto *smh0_684 = buffer.data(smh0 + 684);
    const auto *smh0_686 = buffer.data(smh0 + 686);
    const auto *smh0_687 = buffer.data(smh0 + 687);
    const auto *smh0_689 = buffer.data(smh0 + 689);
    const auto *smh0_690 = buffer.data(smh0 + 690);
    const auto *smh0_691 = buffer.data(smh0 + 691);
    const auto *smh0_692 = buffer.data(smh0 + 692);
    const auto *smh0_693 = buffer.data(smh0 + 693);
    const auto *smh0_696 = buffer.data(smh0 + 696);
    const auto *smh0_698 = buffer.data(smh0 + 698);
    const auto *smh0_699 = buffer.data(smh0 + 699);
    const auto *smh0_702 = buffer.data(smh0 + 702);
    const auto *smh0_703 = buffer.data(smh0 + 703);
    const auto *smh0_705 = buffer.data(smh0 + 705);
    const auto *smh0_707 = buffer.data(smh0 + 707);
    const auto *smh0_708 = buffer.data(smh0 + 708);
    const auto *smh0_710 = buffer.data(smh0 + 710);
    const auto *smh0_711 = buffer.data(smh0 + 711);
    const auto *smh0_712 = buffer.data(smh0 + 712);
    const auto *smh0_713 = buffer.data(smh0 + 713);

    const auto *smh1_654 = buffer.data(smh1 + 654);
    const auto *smh1_656 = buffer.data(smh1 + 656);
    const auto *smh1_657 = buffer.data(smh1 + 657);
    const auto *smh1_660 = buffer.data(smh1 + 660);
    const auto *smh1_661 = buffer.data(smh1 + 661);
    const auto *smh1_663 = buffer.data(smh1 + 663);
    const auto *smh1_665 = buffer.data(smh1 + 665);
    const auto *smh1_666 = buffer.data(smh1 + 666);
    const auto *smh1_668 = buffer.data(smh1 + 668);
    const auto *smh1_669 = buffer.data(smh1 + 669);
    const auto *smh1_670 = buffer.data(smh1 + 670);
    const auto *smh1_671 = buffer.data(smh1 + 671);
    const auto *smh1_672 = buffer.data(smh1 + 672);
    const auto *smh1_675 = buffer.data(smh1 + 675);
    const auto *smh1_677 = buffer.data(smh1 + 677);
    const auto *smh1_678 = buffer.data(smh1 + 678);
    const auto *smh1_681 = buffer.data(smh1 + 681);
    const auto *smh1_682 = buffer.data(smh1 + 682);
    const auto *smh1_684 = buffer.data(smh1 + 684);
    const auto *smh1_686 = buffer.data(smh1 + 686);
    const auto *smh1_687 = buffer.data(smh1 + 687);
    const auto *smh1_689 = buffer.data(smh1 + 689);
    const auto *smh1_690 = buffer.data(smh1 + 690);
    const auto *smh1_691 = buffer.data(smh1 + 691);
    const auto *smh1_692 = buffer.data(smh1 + 692);
    const auto *smh1_693 = buffer.data(smh1 + 693);
    const auto *smh1_696 = buffer.data(smh1 + 696);
    const auto *smh1_698 = buffer.data(smh1 + 698);
    const auto *smh1_699 = buffer.data(smh1 + 699);
    const auto *smh1_702 = buffer.data(smh1 + 702);
    const auto *smh1_703 = buffer.data(smh1 + 703);
    const auto *smh1_705 = buffer.data(smh1 + 705);
    const auto *smh1_707 = buffer.data(smh1 + 707);
    const auto *smh1_708 = buffer.data(smh1 + 708);
    const auto *smh1_710 = buffer.data(smh1 + 710);
    const auto *smh1_711 = buffer.data(smh1 + 711);
    const auto *smh1_712 = buffer.data(smh1 + 712);
    const auto *smh1_713 = buffer.data(smh1 + 713);

    const auto *smi_868 = buffer.data(smi + 868);
    const auto *smi_870 = buffer.data(smi + 870);
    const auto *smi_871 = buffer.data(smi + 871);
    const auto *smi_873 = buffer.data(smi + 873);
    const auto *smi_874 = buffer.data(smi + 874);
    const auto *smi_877 = buffer.data(smi + 877);
    const auto *smi_878 = buffer.data(smi + 878);
    const auto *smi_880 = buffer.data(smi + 880);
    const auto *smi_882 = buffer.data(smi + 882);
    const auto *smi_883 = buffer.data(smi + 883);
    const auto *smi_885 = buffer.data(smi + 885);
    const auto *smi_886 = buffer.data(smi + 886);
    const auto *smi_888 = buffer.data(smi + 888);
    const auto *smi_889 = buffer.data(smi + 889);
    const auto *smi_890 = buffer.data(smi + 890);
    const auto *smi_891 = buffer.data(smi + 891);
    const auto *smi_892 = buffer.data(smi + 892);
    const auto *smi_893 = buffer.data(smi + 893);
    const auto *smi_894 = buffer.data(smi + 894);
    const auto *smi_895 = buffer.data(smi + 895);
    const auto *smi_896 = buffer.data(smi + 896);
    const auto *smi_898 = buffer.data(smi + 898);
    const auto *smi_899 = buffer.data(smi + 899);
    const auto *smi_901 = buffer.data(smi + 901);
    const auto *smi_902 = buffer.data(smi + 902);
    const auto *smi_905 = buffer.data(smi + 905);
    const auto *smi_906 = buffer.data(smi + 906);
    const auto *smi_908 = buffer.data(smi + 908);
    const auto *smi_910 = buffer.data(smi + 910);
    const auto *smi_911 = buffer.data(smi + 911);
    const auto *smi_913 = buffer.data(smi + 913);
    const auto *smi_914 = buffer.data(smi + 914);
    const auto *smi_916 = buffer.data(smi + 916);
    const auto *smi_917 = buffer.data(smi + 917);
    const auto *smi_918 = buffer.data(smi + 918);
    const auto *smi_919 = buffer.data(smi + 919);
    const auto *smi_920 = buffer.data(smi + 920);
    const auto *smi_921 = buffer.data(smi + 921);
    const auto *smi_922 = buffer.data(smi + 922);
    const auto *smi_923 = buffer.data(smi + 923);
    const auto *smi_924 = buffer.data(smi + 924);
    const auto *smi_926 = buffer.data(smi + 926);
    const auto *smi_927 = buffer.data(smi + 927);
    const auto *smi_929 = buffer.data(smi + 929);
    const auto *smi_930 = buffer.data(smi + 930);
    const auto *smi_933 = buffer.data(smi + 933);
    const auto *smi_934 = buffer.data(smi + 934);
    const auto *smi_936 = buffer.data(smi + 936);
    const auto *smi_938 = buffer.data(smi + 938);
    const auto *smi_939 = buffer.data(smi + 939);
    const auto *smi_941 = buffer.data(smi + 941);
    const auto *smi_942 = buffer.data(smi + 942);
    const auto *smi_944 = buffer.data(smi + 944);
    const auto *smi_945 = buffer.data(smi + 945);
    const auto *smi_946 = buffer.data(smi + 946);
    const auto *smi_947 = buffer.data(smi + 947);
    const auto *smi_948 = buffer.data(smi + 948);
    const auto *smi_949 = buffer.data(smi + 949);
    const auto *smi_950 = buffer.data(smi + 950);
    const auto *smi_951 = buffer.data(smi + 951);

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, pc_x, pc_y, pc_z, sli_644, sli_672, \
                         sli_674, sli_871, smh0_654, smh1_654, smi_868, smi_870, \
                         smi_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_16 * sli_672[k]
                    + f_3 * pc_y[k] * smi_868[k];

        t_1118[k] = f_15 * sli_644[k]
                    + f_3 * pc_z[k] * smi_868[k];

        t_1119[k] = f_14 * sli_871[k]
                    + f_4 * smh0_654[k]
                    - f_5 * smh1_654[k]
                    + f_3 * pc_x[k] * smi_871[k];

        t_1120[k] = f_16 * sli_674[k]
                    + f_3 * pc_y[k] * smi_870[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_x, pc_z, sli_647, sli_873, sli_874, \
                         smh0_656, smh0_657, smh1_656, smh1_657, smi_871, smi_873, \
                         smi_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_14 * sli_873[k]
                    + f_4 * smh0_656[k]
                    - f_5 * smh1_656[k]
                    + f_3 * pc_x[k] * smi_873[k];

        t_1122[k] = f_14 * sli_874[k]
                    + f_6 * smh0_657[k]
                    - f_7 * smh1_657[k]
                    + f_3 * pc_x[k] * smi_874[k];

        t_1123[k] = f_15 * sli_647[k]
                    + f_3 * pc_z[k] * smi_871[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_y, sli_677, sli_877, sli_878, \
                         smh0_660, smh0_661, smh1_660, smh1_661, smi_873, smi_877, \
                         smi_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_16 * sli_677[k]
                    + f_3 * pc_y[k] * smi_873[k];

        t_1125[k] = f_14 * sli_877[k]
                    + f_6 * smh0_660[k]
                    - f_7 * smh1_660[k]
                    + f_3 * pc_x[k] * smi_877[k];

        t_1126[k] = f_14 * sli_878[k]
                    + f_8 * smh0_661[k]
                    - f_9 * smh1_661[k]
                    + f_3 * pc_x[k] * smi_878[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, pc_y, pc_z, sli_650, sli_681, sli_880, \
                         smh0_663, smh1_663, smi_874, smi_877, \
                         smi_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_15 * sli_650[k]
                    + f_3 * pc_z[k] * smi_874[k];

        t_1128[k] = f_14 * sli_880[k]
                    + f_8 * smh0_663[k]
                    - f_9 * smh1_663[k]
                    + f_3 * pc_x[k] * smi_880[k];

        t_1129[k] = f_16 * sli_681[k]
                    + f_3 * pc_y[k] * smi_877[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, pc_z, sli_654, sli_882, sli_883, \
                         smh0_665, smh0_666, smh1_665, smh1_666, smi_878, smi_882, \
                         smi_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_14 * sli_882[k]
                    + f_8 * smh0_665[k]
                    - f_9 * smh1_665[k]
                    + f_3 * pc_x[k] * smi_882[k];

        t_1131[k] = f_14 * sli_883[k]
                    + f_10 * smh0_666[k]
                    - f_11 * smh1_666[k]
                    + f_3 * pc_x[k] * smi_883[k];

        t_1132[k] = f_15 * sli_654[k]
                    + f_3 * pc_z[k] * smi_878[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, pc_y, sli_686, sli_885, sli_886, \
                         smh0_668, smh0_669, smh1_668, smh1_669, smi_882, smi_885, \
                         smi_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_14 * sli_885[k]
                    + f_10 * smh0_668[k]
                    - f_11 * smh1_668[k]
                    + f_3 * pc_x[k] * smi_885[k];

        t_1134[k] = f_14 * sli_886[k]
                    + f_10 * smh0_669[k]
                    - f_11 * smh1_669[k]
                    + f_3 * pc_x[k] * smi_886[k];

        t_1135[k] = f_16 * sli_686[k]
                    + f_3 * pc_y[k] * smi_882[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, pc_x, sli_888, sli_889, sli_890, \
                         sli_891, smh0_671, smh1_671, smi_888, smi_889, smi_890, \
                         smi_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_14 * sli_888[k]
                    + f_10 * smh0_671[k]
                    - f_11 * smh1_671[k]
                    + f_3 * pc_x[k] * smi_888[k];

        t_1137[k] = f_14 * sli_889[k]
                    + f_3 * pc_x[k] * smi_889[k];

        t_1138[k] = f_14 * sli_890[k]
                    + f_3 * pc_x[k] * smi_890[k];

        t_1139[k] = f_14 * sli_891[k]
                    + f_3 * pc_x[k] * smi_891[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, pc_x, sli_892, sli_893, sli_894, \
                         sli_895, smi_892, smi_893, smi_894, smi_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_14 * sli_892[k]
                    + f_3 * pc_x[k] * smi_892[k];

        t_1141[k] = f_14 * sli_893[k]
                    + f_3 * pc_x[k] * smi_893[k];

        t_1142[k] = f_14 * sli_894[k]
                    + f_3 * pc_x[k] * smi_894[k];

        t_1143[k] = f_14 * sli_895[k]
                    + f_3 * pc_x[k] * smi_895[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pc_y, pc_z, sli_665, sli_693, sli_695, \
                         smh0_666, smh0_668, smh1_666, smh1_668, smi_889, \
                         smi_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_16 * sli_693[k]
                    + f_1 * smh0_666[k]
                    - f_2 * smh1_666[k]
                    + f_3 * pc_y[k] * smi_889[k];

        t_1145[k] = f_15 * sli_665[k]
                    + f_3 * pc_z[k] * smi_889[k];

        t_1146[k] = f_16 * sli_695[k]
                    + f_4 * smh0_668[k]
                    - f_5 * smh1_668[k]
                    + f_3 * pc_y[k] * smi_891[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pc_y, sli_696, sli_697, sli_698, smh0_669, \
                         smh0_670, smh0_671, smh1_669, smh1_670, smh1_671, smi_892, smi_893, \
                         smi_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_16 * sli_696[k]
                    + f_6 * smh0_669[k]
                    - f_7 * smh1_669[k]
                    + f_3 * pc_y[k] * smi_892[k];

        t_1148[k] = f_16 * sli_697[k]
                    + f_8 * smh0_670[k]
                    - f_9 * smh1_670[k]
                    + f_3 * pc_y[k] * smi_893[k];

        t_1149[k] = f_16 * sli_698[k]
                    + f_10 * smh0_671[k]
                    - f_11 * smh1_671[k]
                    + f_3 * pc_y[k] * smi_894[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, sli_671, sli_699, sli_896, \
                         smh0_671, smh0_672, smh1_671, smh1_672, smi_895, \
                         smi_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_16 * sli_699[k]
                    + f_3 * pc_y[k] * smi_895[k];

        t_1151[k] = f_15 * sli_671[k]
                    + f_1 * smh0_671[k]
                    - f_2 * smh1_671[k]
                    + f_3 * pc_z[k] * smi_895[k];

        t_1152[k] = f_14 * sli_896[k]
                    + f_1 * smh0_672[k]
                    - f_2 * smh1_672[k]
                    + f_3 * pc_x[k] * smi_896[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pc_x, pc_y, pc_z, sli_672, sli_700, \
                         sli_702, sli_899, smh0_675, smh1_675, smi_896, smi_898, \
                         smi_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_15 * sli_700[k]
                    + f_3 * pc_y[k] * smi_896[k];

        t_1154[k] = f_16 * sli_672[k]
                    + f_3 * pc_z[k] * smi_896[k];

        t_1155[k] = f_14 * sli_899[k]
                    + f_4 * smh0_675[k]
                    - f_5 * smh1_675[k]
                    + f_3 * pc_x[k] * smi_899[k];

        t_1156[k] = f_15 * sli_702[k]
                    + f_3 * pc_y[k] * smi_898[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pc_x, pc_z, sli_675, sli_901, sli_902, \
                         smh0_677, smh0_678, smh1_677, smh1_678, smi_899, smi_901, \
                         smi_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_14 * sli_901[k]
                    + f_4 * smh0_677[k]
                    - f_5 * smh1_677[k]
                    + f_3 * pc_x[k] * smi_901[k];

        t_1158[k] = f_14 * sli_902[k]
                    + f_6 * smh0_678[k]
                    - f_7 * smh1_678[k]
                    + f_3 * pc_x[k] * smi_902[k];

        t_1159[k] = f_16 * sli_675[k]
                    + f_3 * pc_z[k] * smi_899[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, pc_y, sli_705, sli_905, sli_906, \
                         smh0_681, smh0_682, smh1_681, smh1_682, smi_901, smi_905, \
                         smi_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_15 * sli_705[k]
                    + f_3 * pc_y[k] * smi_901[k];

        t_1161[k] = f_14 * sli_905[k]
                    + f_6 * smh0_681[k]
                    - f_7 * smh1_681[k]
                    + f_3 * pc_x[k] * smi_905[k];

        t_1162[k] = f_14 * sli_906[k]
                    + f_8 * smh0_682[k]
                    - f_9 * smh1_682[k]
                    + f_3 * pc_x[k] * smi_906[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_x, pc_y, pc_z, sli_678, sli_709, sli_908, \
                         smh0_684, smh1_684, smi_902, smi_905, \
                         smi_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * sli_678[k]
                    + f_3 * pc_z[k] * smi_902[k];

        t_1164[k] = f_14 * sli_908[k]
                    + f_8 * smh0_684[k]
                    - f_9 * smh1_684[k]
                    + f_3 * pc_x[k] * smi_908[k];

        t_1165[k] = f_15 * sli_709[k]
                    + f_3 * pc_y[k] * smi_905[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_x, pc_z, sli_682, sli_910, sli_911, \
                         smh0_686, smh0_687, smh1_686, smh1_687, smi_906, smi_910, \
                         smi_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_14 * sli_910[k]
                    + f_8 * smh0_686[k]
                    - f_9 * smh1_686[k]
                    + f_3 * pc_x[k] * smi_910[k];

        t_1167[k] = f_14 * sli_911[k]
                    + f_10 * smh0_687[k]
                    - f_11 * smh1_687[k]
                    + f_3 * pc_x[k] * smi_911[k];

        t_1168[k] = f_16 * sli_682[k]
                    + f_3 * pc_z[k] * smi_906[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, pc_x, pc_y, sli_714, sli_913, sli_914, \
                         smh0_689, smh0_690, smh1_689, smh1_690, smi_910, smi_913, \
                         smi_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_14 * sli_913[k]
                    + f_10 * smh0_689[k]
                    - f_11 * smh1_689[k]
                    + f_3 * pc_x[k] * smi_913[k];

        t_1170[k] = f_14 * sli_914[k]
                    + f_10 * smh0_690[k]
                    - f_11 * smh1_690[k]
                    + f_3 * pc_x[k] * smi_914[k];

        t_1171[k] = f_15 * sli_714[k]
                    + f_3 * pc_y[k] * smi_910[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, pc_x, sli_916, sli_917, sli_918, \
                         sli_919, smh0_692, smh1_692, smi_916, smi_917, smi_918, \
                         smi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_14 * sli_916[k]
                    + f_10 * smh0_692[k]
                    - f_11 * smh1_692[k]
                    + f_3 * pc_x[k] * smi_916[k];

        t_1173[k] = f_14 * sli_917[k]
                    + f_3 * pc_x[k] * smi_917[k];

        t_1174[k] = f_14 * sli_918[k]
                    + f_3 * pc_x[k] * smi_918[k];

        t_1175[k] = f_14 * sli_919[k]
                    + f_3 * pc_x[k] * smi_919[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, t_1179, pc_x, sli_920, sli_921, sli_922, \
                         sli_923, smi_920, smi_921, smi_922, smi_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_14 * sli_920[k]
                    + f_3 * pc_x[k] * smi_920[k];

        t_1177[k] = f_14 * sli_921[k]
                    + f_3 * pc_x[k] * smi_921[k];

        t_1178[k] = f_14 * sli_922[k]
                    + f_3 * pc_x[k] * smi_922[k];

        t_1179[k] = f_14 * sli_923[k]
                    + f_3 * pc_x[k] * smi_923[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pc_y, pc_z, sli_693, sli_721, sli_723, \
                         smh0_687, smh0_689, smh1_687, smh1_689, smi_917, \
                         smi_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_15 * sli_721[k]
                    + f_1 * smh0_687[k]
                    - f_2 * smh1_687[k]
                    + f_3 * pc_y[k] * smi_917[k];

        t_1181[k] = f_16 * sli_693[k]
                    + f_3 * pc_z[k] * smi_917[k];

        t_1182[k] = f_15 * sli_723[k]
                    + f_4 * smh0_689[k]
                    - f_5 * smh1_689[k]
                    + f_3 * pc_y[k] * smi_919[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, pc_y, sli_724, sli_725, sli_726, smh0_690, \
                         smh0_691, smh0_692, smh1_690, smh1_691, smh1_692, smi_920, smi_921, \
                         smi_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_15 * sli_724[k]
                    + f_6 * smh0_690[k]
                    - f_7 * smh1_690[k]
                    + f_3 * pc_y[k] * smi_920[k];

        t_1184[k] = f_15 * sli_725[k]
                    + f_8 * smh0_691[k]
                    - f_9 * smh1_691[k]
                    + f_3 * pc_y[k] * smi_921[k];

        t_1185[k] = f_15 * sli_726[k]
                    + f_10 * smh0_692[k]
                    - f_11 * smh1_692[k]
                    + f_3 * pc_y[k] * smi_922[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, pc_x, pc_y, pc_z, sli_699, sli_727, sli_924, \
                         smh0_692, smh0_693, smh1_692, smh1_693, smi_923, \
                         smi_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_15 * sli_727[k]
                    + f_3 * pc_y[k] * smi_923[k];

        t_1187[k] = f_16 * sli_699[k]
                    + f_1 * smh0_692[k]
                    - f_2 * smh1_692[k]
                    + f_3 * pc_z[k] * smi_923[k];

        t_1188[k] = f_14 * sli_924[k]
                    + f_1 * smh0_693[k]
                    - f_2 * smh1_693[k]
                    + f_3 * pc_x[k] * smi_924[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pc_x, pc_y, pc_z, sli_700, sli_728, \
                         sli_730, sli_927, smh0_696, smh1_696, smi_924, smi_926, \
                         smi_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_14 * sli_728[k]
                    + f_3 * pc_y[k] * smi_924[k];

        t_1190[k] = f_17 * sli_700[k]
                    + f_3 * pc_z[k] * smi_924[k];

        t_1191[k] = f_14 * sli_927[k]
                    + f_4 * smh0_696[k]
                    - f_5 * smh1_696[k]
                    + f_3 * pc_x[k] * smi_927[k];

        t_1192[k] = f_14 * sli_730[k]
                    + f_3 * pc_y[k] * smi_926[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, pc_x, pc_z, sli_703, sli_929, sli_930, \
                         smh0_698, smh0_699, smh1_698, smh1_699, smi_927, smi_929, \
                         smi_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_14 * sli_929[k]
                    + f_4 * smh0_698[k]
                    - f_5 * smh1_698[k]
                    + f_3 * pc_x[k] * smi_929[k];

        t_1194[k] = f_14 * sli_930[k]
                    + f_6 * smh0_699[k]
                    - f_7 * smh1_699[k]
                    + f_3 * pc_x[k] * smi_930[k];

        t_1195[k] = f_17 * sli_703[k]
                    + f_3 * pc_z[k] * smi_927[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, sli_733, sli_933, sli_934, \
                         smh0_702, smh0_703, smh1_702, smh1_703, smi_929, smi_933, \
                         smi_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_14 * sli_733[k]
                    + f_3 * pc_y[k] * smi_929[k];

        t_1197[k] = f_14 * sli_933[k]
                    + f_6 * smh0_702[k]
                    - f_7 * smh1_702[k]
                    + f_3 * pc_x[k] * smi_933[k];

        t_1198[k] = f_14 * sli_934[k]
                    + f_8 * smh0_703[k]
                    - f_9 * smh1_703[k]
                    + f_3 * pc_x[k] * smi_934[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_x, pc_y, pc_z, sli_706, sli_737, sli_936, \
                         smh0_705, smh1_705, smi_930, smi_933, \
                         smi_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_17 * sli_706[k]
                    + f_3 * pc_z[k] * smi_930[k];

        t_1200[k] = f_14 * sli_936[k]
                    + f_8 * smh0_705[k]
                    - f_9 * smh1_705[k]
                    + f_3 * pc_x[k] * smi_936[k];

        t_1201[k] = f_14 * sli_737[k]
                    + f_3 * pc_y[k] * smi_933[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pc_x, pc_z, sli_710, sli_938, sli_939, \
                         smh0_707, smh0_708, smh1_707, smh1_708, smi_934, smi_938, \
                         smi_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_14 * sli_938[k]
                    + f_8 * smh0_707[k]
                    - f_9 * smh1_707[k]
                    + f_3 * pc_x[k] * smi_938[k];

        t_1203[k] = f_14 * sli_939[k]
                    + f_10 * smh0_708[k]
                    - f_11 * smh1_708[k]
                    + f_3 * pc_x[k] * smi_939[k];

        t_1204[k] = f_17 * sli_710[k]
                    + f_3 * pc_z[k] * smi_934[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, pc_x, pc_y, sli_742, sli_941, sli_942, \
                         smh0_710, smh0_711, smh1_710, smh1_711, smi_938, smi_941, \
                         smi_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_14 * sli_941[k]
                    + f_10 * smh0_710[k]
                    - f_11 * smh1_710[k]
                    + f_3 * pc_x[k] * smi_941[k];

        t_1206[k] = f_14 * sli_942[k]
                    + f_10 * smh0_711[k]
                    - f_11 * smh1_711[k]
                    + f_3 * pc_x[k] * smi_942[k];

        t_1207[k] = f_14 * sli_742[k]
                    + f_3 * pc_y[k] * smi_938[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, t_1211, pc_x, sli_944, sli_945, sli_946, \
                         sli_947, smh0_713, smh1_713, smi_944, smi_945, smi_946, \
                         smi_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_14 * sli_944[k]
                    + f_10 * smh0_713[k]
                    - f_11 * smh1_713[k]
                    + f_3 * pc_x[k] * smi_944[k];

        t_1209[k] = f_14 * sli_945[k]
                    + f_3 * pc_x[k] * smi_945[k];

        t_1210[k] = f_14 * sli_946[k]
                    + f_3 * pc_x[k] * smi_946[k];

        t_1211[k] = f_14 * sli_947[k]
                    + f_3 * pc_x[k] * smi_947[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, t_1215, pc_x, sli_948, sli_949, sli_950, \
                         sli_951, smi_948, smi_949, smi_950, smi_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_14 * sli_948[k]
                    + f_3 * pc_x[k] * smi_948[k];

        t_1213[k] = f_14 * sli_949[k]
                    + f_3 * pc_x[k] * smi_949[k];

        t_1214[k] = f_14 * sli_950[k]
                    + f_3 * pc_x[k] * smi_950[k];

        t_1215[k] = f_14 * sli_951[k]
                    + f_3 * pc_x[k] * smi_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pc_y, pc_z, sli_721, sli_749, sli_751, \
                         smh0_708, smh0_710, smh1_708, smh1_710, smi_945, \
                         smi_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_14 * sli_749[k]
                    + f_1 * smh0_708[k]
                    - f_2 * smh1_708[k]
                    + f_3 * pc_y[k] * smi_945[k];

        t_1217[k] = f_17 * sli_721[k]
                    + f_3 * pc_z[k] * smi_945[k];

        t_1218[k] = f_14 * sli_751[k]
                    + f_4 * smh0_710[k]
                    - f_5 * smh1_710[k]
                    + f_3 * pc_y[k] * smi_947[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, pc_y, sli_752, sli_753, sli_754, smh0_711, \
                         smh0_712, smh0_713, smh1_711, smh1_712, smh1_713, smi_948, smi_949, \
                         smi_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_14 * sli_752[k]
                    + f_6 * smh0_711[k]
                    - f_7 * smh1_711[k]
                    + f_3 * pc_y[k] * smi_948[k];

        t_1220[k] = f_14 * sli_753[k]
                    + f_8 * smh0_712[k]
                    - f_9 * smh1_712[k]
                    + f_3 * pc_y[k] * smi_949[k];

        t_1221[k] = f_14 * sli_754[k]
                    + f_10 * smh0_713[k]
                    - f_11 * smh1_713[k]
                    + f_3 * pc_y[k] * smi_950[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t slk0,
                                                           const size_t sli, const size_t slk1,
                                                           const size_t smh0, const size_t smh1,
                                                           const size_t smi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_972 = buffer.data(slk0 + 972);
    const auto *slk0_975 = buffer.data(slk0 + 975);
    const auto *slk0_977 = buffer.data(slk0 + 977);
    const auto *slk0_978 = buffer.data(slk0 + 978);
    const auto *slk0_981 = buffer.data(slk0 + 981);
    const auto *slk0_982 = buffer.data(slk0 + 982);
    const auto *slk0_984 = buffer.data(slk0 + 984);
    const auto *slk0_986 = buffer.data(slk0 + 986);
    const auto *slk0_987 = buffer.data(slk0 + 987);
    const auto *slk0_989 = buffer.data(slk0 + 989);
    const auto *slk0_990 = buffer.data(slk0 + 990);
    const auto *slk0_992 = buffer.data(slk0 + 992);
    const auto *slk0_1007 = buffer.data(slk0 + 1007);
    const auto *slk0_1008 = buffer.data(slk0 + 1008);
    const auto *slk0_1011 = buffer.data(slk0 + 1011);
    const auto *slk0_1014 = buffer.data(slk0 + 1014);
    const auto *slk0_1296 = buffer.data(slk0 + 1296);
    const auto *slk0_1299 = buffer.data(slk0 + 1299);
    const auto *slk0_1301 = buffer.data(slk0 + 1301);
    const auto *slk0_1302 = buffer.data(slk0 + 1302);
    const auto *slk0_1305 = buffer.data(slk0 + 1305);
    const auto *slk0_1306 = buffer.data(slk0 + 1306);
    const auto *slk0_1308 = buffer.data(slk0 + 1308);
    const auto *slk0_1310 = buffer.data(slk0 + 1310);
    const auto *slk0_1311 = buffer.data(slk0 + 1311);
    const auto *slk0_1313 = buffer.data(slk0 + 1313);
    const auto *slk0_1314 = buffer.data(slk0 + 1314);
    const auto *slk0_1316 = buffer.data(slk0 + 1316);
    const auto *slk0_1324 = buffer.data(slk0 + 1324);
    const auto *slk0_1326 = buffer.data(slk0 + 1326);
    const auto *slk0_1327 = buffer.data(slk0 + 1327);
    const auto *slk0_1328 = buffer.data(slk0 + 1328);
    const auto *slk0_1329 = buffer.data(slk0 + 1329);
    const auto *slk0_1331 = buffer.data(slk0 + 1331);
    const auto *slk0_1337 = buffer.data(slk0 + 1337);

    const auto *sli_727 = buffer.data(sli + 727);
    const auto *sli_728 = buffer.data(sli + 728);
    const auto *sli_731 = buffer.data(sli + 731);
    const auto *sli_734 = buffer.data(sli + 734);
    const auto *sli_738 = buffer.data(sli + 738);
    const auto *sli_749 = buffer.data(sli + 749);
    const auto *sli_755 = buffer.data(sli + 755);
    const auto *sli_756 = buffer.data(sli + 756);
    const auto *sli_757 = buffer.data(sli + 757);
    const auto *sli_758 = buffer.data(sli + 758);
    const auto *sli_759 = buffer.data(sli + 759);
    const auto *sli_761 = buffer.data(sli + 761);
    const auto *sli_762 = buffer.data(sli + 762);
    const auto *sli_764 = buffer.data(sli + 764);
    const auto *sli_765 = buffer.data(sli + 765);
    const auto *sli_766 = buffer.data(sli + 766);
    const auto *sli_768 = buffer.data(sli + 768);
    const auto *sli_769 = buffer.data(sli + 769);
    const auto *sli_770 = buffer.data(sli + 770);
    const auto *sli_777 = buffer.data(sli + 777);
    const auto *sli_779 = buffer.data(sli + 779);
    const auto *sli_780 = buffer.data(sli + 780);
    const auto *sli_781 = buffer.data(sli + 781);
    const auto *sli_782 = buffer.data(sli + 782);
    const auto *sli_783 = buffer.data(sli + 783);
    const auto *sli_784 = buffer.data(sli + 784);
    const auto *sli_786 = buffer.data(sli + 786);
    const auto *sli_789 = buffer.data(sli + 789);
    const auto *sli_793 = buffer.data(sli + 793);
    const auto *sli_798 = buffer.data(sli + 798);
    const auto *sli_811 = buffer.data(sli + 811);
    const auto *sli_812 = buffer.data(sli + 812);
    const auto *sli_814 = buffer.data(sli + 814);
    const auto *sli_973 = buffer.data(sli + 973);
    const auto *sli_974 = buffer.data(sli + 974);
    const auto *sli_975 = buffer.data(sli + 975);
    const auto *sli_976 = buffer.data(sli + 976);
    const auto *sli_977 = buffer.data(sli + 977);
    const auto *sli_978 = buffer.data(sli + 978);
    const auto *sli_979 = buffer.data(sli + 979);
    const auto *sli_980 = buffer.data(sli + 980);
    const auto *sli_983 = buffer.data(sli + 983);
    const auto *sli_985 = buffer.data(sli + 985);
    const auto *sli_986 = buffer.data(sli + 986);
    const auto *sli_989 = buffer.data(sli + 989);
    const auto *sli_990 = buffer.data(sli + 990);
    const auto *sli_992 = buffer.data(sli + 992);
    const auto *sli_994 = buffer.data(sli + 994);
    const auto *sli_995 = buffer.data(sli + 995);
    const auto *sli_997 = buffer.data(sli + 997);
    const auto *sli_998 = buffer.data(sli + 998);
    const auto *sli_1000 = buffer.data(sli + 1000);
    const auto *sli_1001 = buffer.data(sli + 1001);
    const auto *sli_1002 = buffer.data(sli + 1002);
    const auto *sli_1003 = buffer.data(sli + 1003);
    const auto *sli_1004 = buffer.data(sli + 1004);
    const auto *sli_1005 = buffer.data(sli + 1005);
    const auto *sli_1006 = buffer.data(sli + 1006);
    const auto *sli_1007 = buffer.data(sli + 1007);
    const auto *sli_1008 = buffer.data(sli + 1008);
    const auto *sli_1011 = buffer.data(sli + 1011);
    const auto *sli_1013 = buffer.data(sli + 1013);
    const auto *sli_1014 = buffer.data(sli + 1014);
    const auto *sli_1017 = buffer.data(sli + 1017);
    const auto *sli_1018 = buffer.data(sli + 1018);
    const auto *sli_1020 = buffer.data(sli + 1020);
    const auto *sli_1022 = buffer.data(sli + 1022);
    const auto *sli_1023 = buffer.data(sli + 1023);
    const auto *sli_1025 = buffer.data(sli + 1025);
    const auto *sli_1026 = buffer.data(sli + 1026);
    const auto *sli_1028 = buffer.data(sli + 1028);
    const auto *sli_1029 = buffer.data(sli + 1029);
    const auto *sli_1030 = buffer.data(sli + 1030);
    const auto *sli_1031 = buffer.data(sli + 1031);
    const auto *sli_1032 = buffer.data(sli + 1032);
    const auto *sli_1033 = buffer.data(sli + 1033);
    const auto *sli_1034 = buffer.data(sli + 1034);
    const auto *sli_1035 = buffer.data(sli + 1035);
    const auto *sli_1041 = buffer.data(sli + 1041);

    const auto *slk1_972 = buffer.data(slk1 + 972);
    const auto *slk1_975 = buffer.data(slk1 + 975);
    const auto *slk1_977 = buffer.data(slk1 + 977);
    const auto *slk1_978 = buffer.data(slk1 + 978);
    const auto *slk1_981 = buffer.data(slk1 + 981);
    const auto *slk1_982 = buffer.data(slk1 + 982);
    const auto *slk1_984 = buffer.data(slk1 + 984);
    const auto *slk1_986 = buffer.data(slk1 + 986);
    const auto *slk1_987 = buffer.data(slk1 + 987);
    const auto *slk1_989 = buffer.data(slk1 + 989);
    const auto *slk1_990 = buffer.data(slk1 + 990);
    const auto *slk1_992 = buffer.data(slk1 + 992);
    const auto *slk1_1007 = buffer.data(slk1 + 1007);
    const auto *slk1_1008 = buffer.data(slk1 + 1008);
    const auto *slk1_1011 = buffer.data(slk1 + 1011);
    const auto *slk1_1014 = buffer.data(slk1 + 1014);
    const auto *slk1_1296 = buffer.data(slk1 + 1296);
    const auto *slk1_1299 = buffer.data(slk1 + 1299);
    const auto *slk1_1301 = buffer.data(slk1 + 1301);
    const auto *slk1_1302 = buffer.data(slk1 + 1302);
    const auto *slk1_1305 = buffer.data(slk1 + 1305);
    const auto *slk1_1306 = buffer.data(slk1 + 1306);
    const auto *slk1_1308 = buffer.data(slk1 + 1308);
    const auto *slk1_1310 = buffer.data(slk1 + 1310);
    const auto *slk1_1311 = buffer.data(slk1 + 1311);
    const auto *slk1_1313 = buffer.data(slk1 + 1313);
    const auto *slk1_1314 = buffer.data(slk1 + 1314);
    const auto *slk1_1316 = buffer.data(slk1 + 1316);
    const auto *slk1_1324 = buffer.data(slk1 + 1324);
    const auto *slk1_1326 = buffer.data(slk1 + 1326);
    const auto *slk1_1327 = buffer.data(slk1 + 1327);
    const auto *slk1_1328 = buffer.data(slk1 + 1328);
    const auto *slk1_1329 = buffer.data(slk1 + 1329);
    const auto *slk1_1331 = buffer.data(slk1 + 1331);
    const auto *slk1_1337 = buffer.data(slk1 + 1337);

    const auto *smh0_713 = buffer.data(smh0 + 713);
    const auto *smh0_729 = buffer.data(smh0 + 729);
    const auto *smh0_731 = buffer.data(smh0 + 731);
    const auto *smh0_732 = buffer.data(smh0 + 732);
    const auto *smh0_733 = buffer.data(smh0 + 733);
    const auto *smh0_734 = buffer.data(smh0 + 734);
    const auto *smh0_735 = buffer.data(smh0 + 735);
    const auto *smh0_738 = buffer.data(smh0 + 738);
    const auto *smh0_740 = buffer.data(smh0 + 740);
    const auto *smh0_741 = buffer.data(smh0 + 741);
    const auto *smh0_744 = buffer.data(smh0 + 744);
    const auto *smh0_745 = buffer.data(smh0 + 745);
    const auto *smh0_747 = buffer.data(smh0 + 747);
    const auto *smh0_749 = buffer.data(smh0 + 749);
    const auto *smh0_750 = buffer.data(smh0 + 750);
    const auto *smh0_752 = buffer.data(smh0 + 752);
    const auto *smh0_753 = buffer.data(smh0 + 753);
    const auto *smh0_754 = buffer.data(smh0 + 754);
    const auto *smh0_755 = buffer.data(smh0 + 755);

    const auto *smh1_713 = buffer.data(smh1 + 713);
    const auto *smh1_729 = buffer.data(smh1 + 729);
    const auto *smh1_731 = buffer.data(smh1 + 731);
    const auto *smh1_732 = buffer.data(smh1 + 732);
    const auto *smh1_733 = buffer.data(smh1 + 733);
    const auto *smh1_734 = buffer.data(smh1 + 734);
    const auto *smh1_735 = buffer.data(smh1 + 735);
    const auto *smh1_738 = buffer.data(smh1 + 738);
    const auto *smh1_740 = buffer.data(smh1 + 740);
    const auto *smh1_741 = buffer.data(smh1 + 741);
    const auto *smh1_744 = buffer.data(smh1 + 744);
    const auto *smh1_745 = buffer.data(smh1 + 745);
    const auto *smh1_747 = buffer.data(smh1 + 747);
    const auto *smh1_749 = buffer.data(smh1 + 749);
    const auto *smh1_750 = buffer.data(smh1 + 750);
    const auto *smh1_752 = buffer.data(smh1 + 752);
    const auto *smh1_753 = buffer.data(smh1 + 753);
    const auto *smh1_754 = buffer.data(smh1 + 754);
    const auto *smh1_755 = buffer.data(smh1 + 755);

    const auto *smi_951 = buffer.data(smi + 951);
    const auto *smi_952 = buffer.data(smi + 952);
    const auto *smi_954 = buffer.data(smi + 954);
    const auto *smi_955 = buffer.data(smi + 955);
    const auto *smi_957 = buffer.data(smi + 957);
    const auto *smi_958 = buffer.data(smi + 958);
    const auto *smi_961 = buffer.data(smi + 961);
    const auto *smi_962 = buffer.data(smi + 962);
    const auto *smi_966 = buffer.data(smi + 966);
    const auto *smi_973 = buffer.data(smi + 973);
    const auto *smi_974 = buffer.data(smi + 974);
    const auto *smi_975 = buffer.data(smi + 975);
    const auto *smi_976 = buffer.data(smi + 976);
    const auto *smi_977 = buffer.data(smi + 977);
    const auto *smi_978 = buffer.data(smi + 978);
    const auto *smi_979 = buffer.data(smi + 979);
    const auto *smi_980 = buffer.data(smi + 980);
    const auto *smi_982 = buffer.data(smi + 982);
    const auto *smi_983 = buffer.data(smi + 983);
    const auto *smi_985 = buffer.data(smi + 985);
    const auto *smi_986 = buffer.data(smi + 986);
    const auto *smi_989 = buffer.data(smi + 989);
    const auto *smi_990 = buffer.data(smi + 990);
    const auto *smi_992 = buffer.data(smi + 992);
    const auto *smi_994 = buffer.data(smi + 994);
    const auto *smi_995 = buffer.data(smi + 995);
    const auto *smi_997 = buffer.data(smi + 997);
    const auto *smi_998 = buffer.data(smi + 998);
    const auto *smi_1000 = buffer.data(smi + 1000);
    const auto *smi_1001 = buffer.data(smi + 1001);
    const auto *smi_1002 = buffer.data(smi + 1002);
    const auto *smi_1003 = buffer.data(smi + 1003);
    const auto *smi_1004 = buffer.data(smi + 1004);
    const auto *smi_1005 = buffer.data(smi + 1005);
    const auto *smi_1006 = buffer.data(smi + 1006);
    const auto *smi_1007 = buffer.data(smi + 1007);
    const auto *smi_1008 = buffer.data(smi + 1008);
    const auto *smi_1010 = buffer.data(smi + 1010);
    const auto *smi_1011 = buffer.data(smi + 1011);
    const auto *smi_1013 = buffer.data(smi + 1013);
    const auto *smi_1014 = buffer.data(smi + 1014);
    const auto *smi_1017 = buffer.data(smi + 1017);
    const auto *smi_1018 = buffer.data(smi + 1018);
    const auto *smi_1022 = buffer.data(smi + 1022);
    const auto *smi_1029 = buffer.data(smi + 1029);
    const auto *smi_1030 = buffer.data(smi + 1030);
    const auto *smi_1031 = buffer.data(smi + 1031);
    const auto *smi_1032 = buffer.data(smi + 1032);
    const auto *smi_1033 = buffer.data(smi + 1033);
    const auto *smi_1034 = buffer.data(smi + 1034);
    const auto *smi_1035 = buffer.data(smi + 1035);
    const auto *smi_1036 = buffer.data(smi + 1036);
    const auto *smi_1038 = buffer.data(smi + 1038);

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pb_y, pc_y, pc_z, slk0_972, sli_727, \
                         sli_755, sli_756, slk1_972, smh0_713, smh1_713, smi_951, \
                         smi_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_14 * sli_755[k]
                    + f_3 * pc_y[k] * smi_951[k];

        t_1223[k] = f_17 * sli_727[k]
                    + f_1 * smh0_713[k]
                    - f_2 * smh1_713[k]
                    + f_3 * pc_z[k] * smi_951[k];

        t_1224[k] = pb_y[k] * slk0_972[k]
                    - f_12 * pc_y[k] * slk1_972[k];

        t_1225[k] = f_13 * sli_756[k]
                    + f_3 * pc_y[k] * smi_952[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, pb_y, pc_y, pc_z, slk0_975, slk0_977, \
                         sli_728, sli_757, sli_758, slk1_975, slk1_977, smi_952, \
                         smi_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_20 * sli_728[k]
                    + f_3 * pc_z[k] * smi_952[k];

        t_1227[k] = pb_y[k] * slk0_975[k]
                    + f_14 * sli_757[k]
                    - f_12 * pc_y[k] * slk1_975[k];

        t_1228[k] = f_13 * sli_758[k]
                    + f_3 * pc_y[k] * smi_954[k];

        t_1229[k] = pb_y[k] * slk0_977[k]
                    - f_12 * pc_y[k] * slk1_977[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pb_y, pc_y, pc_z, slk0_978, slk0_981, \
                         sli_731, sli_759, sli_761, slk1_978, slk1_981, smi_955, \
                         smi_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = pb_y[k] * slk0_978[k]
                    + f_15 * sli_759[k]
                    - f_12 * pc_y[k] * slk1_978[k];

        t_1231[k] = f_20 * sli_731[k]
                    + f_3 * pc_z[k] * smi_955[k];

        t_1232[k] = f_13 * sli_761[k]
                    + f_3 * pc_y[k] * smi_957[k];

        t_1233[k] = pb_y[k] * slk0_981[k]
                    - f_12 * pc_y[k] * slk1_981[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pb_y, pc_y, pc_z, slk0_982, slk0_984, \
                         sli_734, sli_762, sli_764, slk1_982, slk1_984, \
                         smi_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = pb_y[k] * slk0_982[k]
                    + f_16 * sli_762[k]
                    - f_12 * pc_y[k] * slk1_982[k];

        t_1235[k] = f_20 * sli_734[k]
                    + f_3 * pc_z[k] * smi_958[k];

        t_1236[k] = pb_y[k] * slk0_984[k]
                    + f_14 * sli_764[k]
                    - f_12 * pc_y[k] * slk1_984[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, pb_y, pc_y, pc_z, slk0_986, slk0_987, \
                         sli_738, sli_765, sli_766, slk1_986, slk1_987, smi_961, \
                         smi_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_13 * sli_765[k]
                    + f_3 * pc_y[k] * smi_961[k];

        t_1238[k] = pb_y[k] * slk0_986[k]
                    - f_12 * pc_y[k] * slk1_986[k];

        t_1239[k] = pb_y[k] * slk0_987[k]
                    + f_17 * sli_766[k]
                    - f_12 * pc_y[k] * slk1_987[k];

        t_1240[k] = f_20 * sli_738[k]
                    + f_3 * pc_z[k] * smi_962[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, t_1244, pb_y, pc_y, slk0_989, slk0_990, \
                         slk0_992, sli_768, sli_769, sli_770, slk1_989, slk1_990, slk1_992, \
                         smi_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = pb_y[k] * slk0_989[k]
                    + f_15 * sli_768[k]
                    - f_12 * pc_y[k] * slk1_989[k];

        t_1242[k] = pb_y[k] * slk0_990[k]
                    + f_14 * sli_769[k]
                    - f_12 * pc_y[k] * slk1_990[k];

        t_1243[k] = f_13 * sli_770[k]
                    + f_3 * pc_y[k] * smi_966[k];

        t_1244[k] = pb_y[k] * slk0_992[k]
                    - f_12 * pc_y[k] * slk1_992[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, t_1248, t_1249, pc_x, sli_973, sli_974, \
                         sli_975, sli_976, sli_977, smi_973, smi_974, smi_975, smi_976, \
                         smi_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = f_14 * sli_973[k]
                    + f_3 * pc_x[k] * smi_973[k];

        t_1246[k] = f_14 * sli_974[k]
                    + f_3 * pc_x[k] * smi_974[k];

        t_1247[k] = f_14 * sli_975[k]
                    + f_3 * pc_x[k] * smi_975[k];

        t_1248[k] = f_14 * sli_976[k]
                    + f_3 * pc_x[k] * smi_976[k];

        t_1249[k] = f_14 * sli_977[k]
                    + f_3 * pc_x[k] * smi_977[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pc_x, pc_y, pc_z, sli_749, sli_777, \
                         sli_978, sli_979, smh0_729, smh1_729, smi_973, smi_978, \
                         smi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_14 * sli_978[k]
                    + f_3 * pc_x[k] * smi_978[k];

        t_1251[k] = f_14 * sli_979[k]
                    + f_3 * pc_x[k] * smi_979[k];

        t_1252[k] = f_13 * sli_777[k]
                    + f_1 * smh0_729[k]
                    - f_2 * smh1_729[k]
                    + f_3 * pc_y[k] * smi_973[k];

        t_1253[k] = f_20 * sli_749[k]
                    + f_3 * pc_z[k] * smi_973[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, sli_779, sli_780, sli_781, smh0_731, \
                         smh0_732, smh0_733, smh1_731, smh1_732, smh1_733, smi_975, smi_976, \
                         smi_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_13 * sli_779[k]
                    + f_4 * smh0_731[k]
                    - f_5 * smh1_731[k]
                    + f_3 * pc_y[k] * smi_975[k];

        t_1255[k] = f_13 * sli_780[k]
                    + f_6 * smh0_732[k]
                    - f_7 * smh1_732[k]
                    + f_3 * pc_y[k] * smi_976[k];

        t_1256[k] = f_13 * sli_781[k]
                    + f_8 * smh0_733[k]
                    - f_9 * smh1_733[k]
                    + f_3 * pc_y[k] * smi_977[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pb_y, pc_y, slk0_1007, sli_782, sli_783, \
                         slk1_1007, smh0_734, smh1_734, smi_978, \
                         smi_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_13 * sli_782[k]
                    + f_10 * smh0_734[k]
                    - f_11 * smh1_734[k]
                    + f_3 * pc_y[k] * smi_978[k];

        t_1258[k] = f_13 * sli_783[k]
                    + f_3 * pc_y[k] * smi_979[k];

        t_1259[k] = pb_y[k] * slk0_1007[k]
                    - f_12 * pc_y[k] * slk1_1007[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pc_x, pc_y, pc_z, sli_756, sli_980, \
                         sli_983, smh0_735, smh0_738, smh1_735, smh1_738, smi_980, \
                         smi_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_14 * sli_980[k]
                    + f_1 * smh0_735[k]
                    - f_2 * smh1_735[k]
                    + f_3 * pc_x[k] * smi_980[k];

        t_1261[k] = f_3 * pc_y[k] * smi_980[k];

        t_1262[k] = f_19 * sli_756[k]
                    + f_3 * pc_z[k] * smi_980[k];

        t_1263[k] = f_14 * sli_983[k]
                    + f_4 * smh0_738[k]
                    - f_5 * smh1_738[k]
                    + f_3 * pc_x[k] * smi_983[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, pc_x, pc_y, sli_985, sli_986, smh0_740, \
                         smh0_741, smh1_740, smh1_741, smi_982, smi_985, \
                         smi_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_3 * pc_y[k] * smi_982[k];

        t_1265[k] = f_14 * sli_985[k]
                    + f_4 * smh0_740[k]
                    - f_5 * smh1_740[k]
                    + f_3 * pc_x[k] * smi_985[k];

        t_1266[k] = f_14 * sli_986[k]
                    + f_6 * smh0_741[k]
                    - f_7 * smh1_741[k]
                    + f_3 * pc_x[k] * smi_986[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, pc_x, pc_y, pc_z, sli_759, sli_989, smh0_744, \
                         smh1_744, smi_983, smi_985, smi_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_19 * sli_759[k]
                    + f_3 * pc_z[k] * smi_983[k];

        t_1268[k] = f_3 * pc_y[k] * smi_985[k];

        t_1269[k] = f_14 * sli_989[k]
                    + f_6 * smh0_744[k]
                    - f_7 * smh1_744[k]
                    + f_3 * pc_x[k] * smi_989[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, pc_x, pc_z, sli_762, sli_990, sli_992, \
                         smh0_745, smh0_747, smh1_745, smh1_747, smi_986, smi_990, \
                         smi_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = f_14 * sli_990[k]
                    + f_8 * smh0_745[k]
                    - f_9 * smh1_745[k]
                    + f_3 * pc_x[k] * smi_990[k];

        t_1271[k] = f_19 * sli_762[k]
                    + f_3 * pc_z[k] * smi_986[k];

        t_1272[k] = f_14 * sli_992[k]
                    + f_8 * smh0_747[k]
                    - f_9 * smh1_747[k]
                    + f_3 * pc_x[k] * smi_992[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, pc_x, pc_y, sli_994, sli_995, smh0_749, \
                         smh0_750, smh1_749, smh1_750, smi_989, smi_994, \
                         smi_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_3 * pc_y[k] * smi_989[k];

        t_1274[k] = f_14 * sli_994[k]
                    + f_8 * smh0_749[k]
                    - f_9 * smh1_749[k]
                    + f_3 * pc_x[k] * smi_994[k];

        t_1275[k] = f_14 * sli_995[k]
                    + f_10 * smh0_750[k]
                    - f_11 * smh1_750[k]
                    + f_3 * pc_x[k] * smi_995[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, pc_x, pc_z, sli_766, sli_997, sli_998, \
                         smh0_752, smh0_753, smh1_752, smh1_753, smi_990, smi_997, \
                         smi_998 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_19 * sli_766[k]
                    + f_3 * pc_z[k] * smi_990[k];

        t_1277[k] = f_14 * sli_997[k]
                    + f_10 * smh0_752[k]
                    - f_11 * smh1_752[k]
                    + f_3 * pc_x[k] * smi_997[k];

        t_1278[k] = f_14 * sli_998[k]
                    + f_10 * smh0_753[k]
                    - f_11 * smh1_753[k]
                    + f_3 * pc_x[k] * smi_998[k];
    }

#pragma omp simd aligned(t_1279, t_1280, t_1281, t_1282, pc_x, pc_y, sli_1000, sli_1001, \
                         sli_1002, smh0_755, smh1_755, smi_994, smi_1000, smi_1001, \
                         smi_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1279[k] = f_3 * pc_y[k] * smi_994[k];

        t_1280[k] = f_14 * sli_1000[k]
                    + f_10 * smh0_755[k]
                    - f_11 * smh1_755[k]
                    + f_3 * pc_x[k] * smi_1000[k];

        t_1281[k] = f_14 * sli_1001[k]
                    + f_3 * pc_x[k] * smi_1001[k];

        t_1282[k] = f_14 * sli_1002[k]
                    + f_3 * pc_x[k] * smi_1002[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, t_1286, t_1287, pc_x, sli_1003, sli_1004, \
                         sli_1005, sli_1006, sli_1007, smi_1003, smi_1004, smi_1005, smi_1006, \
                         smi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = f_14 * sli_1003[k]
                    + f_3 * pc_x[k] * smi_1003[k];

        t_1284[k] = f_14 * sli_1004[k]
                    + f_3 * pc_x[k] * smi_1004[k];

        t_1285[k] = f_14 * sli_1005[k]
                    + f_3 * pc_x[k] * smi_1005[k];

        t_1286[k] = f_14 * sli_1006[k]
                    + f_3 * pc_x[k] * smi_1006[k];

        t_1287[k] = f_14 * sli_1007[k]
                    + f_3 * pc_x[k] * smi_1007[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pc_y, pc_z, sli_777, smh0_750, \
                         smh0_752, smh0_753, smh1_750, smh1_752, smh1_753, smi_1001, smi_1003, \
                         smi_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = f_1 * smh0_750[k]
                    - f_2 * smh1_750[k]
                    + f_3 * pc_y[k] * smi_1001[k];

        t_1289[k] = f_19 * sli_777[k]
                    + f_3 * pc_z[k] * smi_1001[k];

        t_1290[k] = f_4 * smh0_752[k]
                    - f_5 * smh1_752[k]
                    + f_3 * pc_y[k] * smi_1003[k];

        t_1291[k] = f_6 * smh0_753[k]
                    - f_7 * smh1_753[k]
                    + f_3 * pc_y[k] * smi_1004[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, t_1295, pc_y, pc_z, sli_783, smh0_754, \
                         smh0_755, smh1_754, smh1_755, smi_1005, smi_1006, \
                         smi_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_8 * smh0_754[k]
                    - f_9 * smh1_754[k]
                    + f_3 * pc_y[k] * smi_1005[k];

        t_1293[k] = f_10 * smh0_755[k]
                    - f_11 * smh1_755[k]
                    + f_3 * pc_y[k] * smi_1006[k];

        t_1294[k] = f_3 * pc_y[k] * smi_1007[k];

        t_1295[k] = f_19 * sli_783[k]
                    + f_1 * smh0_755[k]
                    - f_2 * smh1_755[k]
                    + f_3 * pc_z[k] * smi_1007[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, t_1299, pb_x, pc_x, pc_y, pc_z, slk0_1296, \
                         slk0_1299, sli_784, sli_1008, sli_1011, slk1_1296, slk1_1299, \
                         smi_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = pb_x[k] * slk0_1296[k]
                    + f_19 * sli_1008[k]
                    - f_12 * pc_x[k] * slk1_1296[k];

        t_1297[k] = f_18 * sli_784[k]
                    + f_3 * pc_y[k] * smi_1008[k];

        t_1298[k] = f_3 * pc_z[k] * smi_1008[k];

        t_1299[k] = pb_x[k] * slk0_1299[k]
                    + f_17 * sli_1011[k]
                    - f_12 * pc_x[k] * slk1_1299[k];
    }

#pragma omp simd aligned(t_1300, t_1301, t_1302, pb_x, pc_x, pc_y, slk0_1301, slk0_1302, \
                         sli_786, sli_1013, sli_1014, slk1_1301, slk1_1302, \
                         smi_1010 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1300[k] = f_18 * sli_786[k]
                    + f_3 * pc_y[k] * smi_1010[k];

        t_1301[k] = pb_x[k] * slk0_1301[k]
                    + f_17 * sli_1013[k]
                    - f_12 * pc_x[k] * slk1_1301[k];

        t_1302[k] = pb_x[k] * slk0_1302[k]
                    + f_16 * sli_1014[k]
                    - f_12 * pc_x[k] * slk1_1302[k];
    }

#pragma omp simd aligned(t_1303, t_1304, t_1305, pb_x, pc_x, pc_y, pc_z, slk0_1305, sli_789, \
                         sli_1017, slk1_1305, smi_1011, smi_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1303[k] = f_3 * pc_z[k] * smi_1011[k];

        t_1304[k] = f_18 * sli_789[k]
                    + f_3 * pc_y[k] * smi_1013[k];

        t_1305[k] = pb_x[k] * slk0_1305[k]
                    + f_16 * sli_1017[k]
                    - f_12 * pc_x[k] * slk1_1305[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, pb_x, pc_x, pc_z, slk0_1306, slk0_1308, \
                         sli_1018, sli_1020, slk1_1306, slk1_1308, \
                         smi_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = pb_x[k] * slk0_1306[k]
                    + f_15 * sli_1018[k]
                    - f_12 * pc_x[k] * slk1_1306[k];

        t_1307[k] = f_3 * pc_z[k] * smi_1014[k];

        t_1308[k] = pb_x[k] * slk0_1308[k]
                    + f_15 * sli_1020[k]
                    - f_12 * pc_x[k] * slk1_1308[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, pb_x, pc_x, pc_y, slk0_1310, slk0_1311, \
                         sli_793, sli_1022, sli_1023, slk1_1310, slk1_1311, \
                         smi_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = f_18 * sli_793[k]
                    + f_3 * pc_y[k] * smi_1017[k];

        t_1310[k] = pb_x[k] * slk0_1310[k]
                    + f_15 * sli_1022[k]
                    - f_12 * pc_x[k] * slk1_1310[k];

        t_1311[k] = pb_x[k] * slk0_1311[k]
                    + f_14 * sli_1023[k]
                    - f_12 * pc_x[k] * slk1_1311[k];
    }

#pragma omp simd aligned(t_1312, t_1313, t_1314, pb_x, pc_x, pc_z, slk0_1313, slk0_1314, \
                         sli_1025, sli_1026, slk1_1313, slk1_1314, \
                         smi_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1312[k] = f_3 * pc_z[k] * smi_1018[k];

        t_1313[k] = pb_x[k] * slk0_1313[k]
                    + f_14 * sli_1025[k]
                    - f_12 * pc_x[k] * slk1_1313[k];

        t_1314[k] = pb_x[k] * slk0_1314[k]
                    + f_14 * sli_1026[k]
                    - f_12 * pc_x[k] * slk1_1314[k];
    }

#pragma omp simd aligned(t_1315, t_1316, t_1317, t_1318, pb_x, pc_x, pc_y, slk0_1316, sli_798, \
                         sli_1028, sli_1029, sli_1030, slk1_1316, smi_1022, smi_1029, \
                         smi_1030 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1315[k] = f_18 * sli_798[k]
                    + f_3 * pc_y[k] * smi_1022[k];

        t_1316[k] = pb_x[k] * slk0_1316[k]
                    + f_14 * sli_1028[k]
                    - f_12 * pc_x[k] * slk1_1316[k];

        t_1317[k] = f_13 * sli_1029[k]
                    + f_3 * pc_x[k] * smi_1029[k];

        t_1318[k] = f_13 * sli_1030[k]
                    + f_3 * pc_x[k] * smi_1030[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, t_1322, t_1323, pc_x, sli_1031, sli_1032, \
                         sli_1033, sli_1034, sli_1035, smi_1031, smi_1032, smi_1033, smi_1034, \
                         smi_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_13 * sli_1031[k]
                    + f_3 * pc_x[k] * smi_1031[k];

        t_1320[k] = f_13 * sli_1032[k]
                    + f_3 * pc_x[k] * smi_1032[k];

        t_1321[k] = f_13 * sli_1033[k]
                    + f_3 * pc_x[k] * smi_1033[k];

        t_1322[k] = f_13 * sli_1034[k]
                    + f_3 * pc_x[k] * smi_1034[k];

        t_1323[k] = f_13 * sli_1035[k]
                    + f_3 * pc_x[k] * smi_1035[k];
    }

#pragma omp simd aligned(t_1324, t_1325, t_1326, t_1327, pb_x, pc_x, pc_z, slk0_1324, \
                         slk0_1326, slk0_1327, slk1_1324, slk1_1326, slk1_1327, \
                         smi_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1324[k] = pb_x[k] * slk0_1324[k]
                    - f_12 * pc_x[k] * slk1_1324[k];

        t_1325[k] = f_3 * pc_z[k] * smi_1029[k];

        t_1326[k] = pb_x[k] * slk0_1326[k]
                    - f_12 * pc_x[k] * slk1_1326[k];

        t_1327[k] = pb_x[k] * slk0_1327[k]
                    - f_12 * pc_x[k] * slk1_1327[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, t_1331, pb_x, pc_x, pc_y, slk0_1328, \
                         slk0_1329, slk0_1331, sli_811, slk1_1328, slk1_1329, slk1_1331, \
                         smi_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = pb_x[k] * slk0_1328[k]
                    - f_12 * pc_x[k] * slk1_1328[k];

        t_1329[k] = pb_x[k] * slk0_1329[k]
                    - f_12 * pc_x[k] * slk1_1329[k];

        t_1330[k] = f_18 * sli_811[k]
                    + f_3 * pc_y[k] * smi_1035[k];

        t_1331[k] = pb_x[k] * slk0_1331[k]
                    - f_12 * pc_x[k] * slk1_1331[k];
    }

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, pb_z, pc_y, pc_z, slk0_1008, \
                         slk0_1011, sli_784, sli_812, slk1_1008, slk1_1011, \
                         smi_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = pb_z[k] * slk0_1008[k]
                    - f_12 * pc_z[k] * slk1_1008[k];

        t_1333[k] = f_19 * sli_812[k]
                    + f_3 * pc_y[k] * smi_1036[k];

        t_1334[k] = f_13 * sli_784[k]
                    + f_3 * pc_z[k] * smi_1036[k];

        t_1335[k] = pb_z[k] * slk0_1011[k]
                    - f_12 * pc_z[k] * slk1_1011[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, pb_x, pb_z, pc_x, pc_y, pc_z, slk0_1014, \
                         slk0_1337, sli_814, sli_1041, slk1_1014, slk1_1337, \
                         smi_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_19 * sli_814[k]
                    + f_3 * pc_y[k] * smi_1038[k];

        t_1337[k] = pb_x[k] * slk0_1337[k]
                    + f_17 * sli_1041[k]
                    - f_12 * pc_x[k] * slk1_1337[k];

        t_1338[k] = pb_z[k] * slk0_1014[k]
                    - f_12 * pc_z[k] * slk1_1014[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t slk0,
                                                           const size_t sli, const size_t slk1,
                                                           const size_t smi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_1018 = buffer.data(slk0 + 1018);
    const auto *slk0_1023 = buffer.data(slk0 + 1023);
    const auto *slk0_1341 = buffer.data(slk0 + 1341);
    const auto *slk0_1344 = buffer.data(slk0 + 1344);
    const auto *slk0_1346 = buffer.data(slk0 + 1346);
    const auto *slk0_1349 = buffer.data(slk0 + 1349);
    const auto *slk0_1350 = buffer.data(slk0 + 1350);
    const auto *slk0_1352 = buffer.data(slk0 + 1352);
    const auto *slk0_1360 = buffer.data(slk0 + 1360);
    const auto *slk0_1362 = buffer.data(slk0 + 1362);
    const auto *slk0_1363 = buffer.data(slk0 + 1363);
    const auto *slk0_1364 = buffer.data(slk0 + 1364);
    const auto *slk0_1365 = buffer.data(slk0 + 1365);
    const auto *slk0_1367 = buffer.data(slk0 + 1367);
    const auto *slk0_1368 = buffer.data(slk0 + 1368);
    const auto *slk0_1371 = buffer.data(slk0 + 1371);
    const auto *slk0_1373 = buffer.data(slk0 + 1373);
    const auto *slk0_1374 = buffer.data(slk0 + 1374);
    const auto *slk0_1377 = buffer.data(slk0 + 1377);
    const auto *slk0_1378 = buffer.data(slk0 + 1378);
    const auto *slk0_1380 = buffer.data(slk0 + 1380);
    const auto *slk0_1382 = buffer.data(slk0 + 1382);
    const auto *slk0_1383 = buffer.data(slk0 + 1383);
    const auto *slk0_1385 = buffer.data(slk0 + 1385);
    const auto *slk0_1386 = buffer.data(slk0 + 1386);
    const auto *slk0_1388 = buffer.data(slk0 + 1388);
    const auto *slk0_1396 = buffer.data(slk0 + 1396);
    const auto *slk0_1398 = buffer.data(slk0 + 1398);
    const auto *slk0_1399 = buffer.data(slk0 + 1399);
    const auto *slk0_1400 = buffer.data(slk0 + 1400);
    const auto *slk0_1401 = buffer.data(slk0 + 1401);
    const auto *slk0_1403 = buffer.data(slk0 + 1403);
    const auto *slk0_1404 = buffer.data(slk0 + 1404);
    const auto *slk0_1407 = buffer.data(slk0 + 1407);
    const auto *slk0_1409 = buffer.data(slk0 + 1409);
    const auto *slk0_1410 = buffer.data(slk0 + 1410);
    const auto *slk0_1413 = buffer.data(slk0 + 1413);
    const auto *slk0_1414 = buffer.data(slk0 + 1414);
    const auto *slk0_1416 = buffer.data(slk0 + 1416);
    const auto *slk0_1418 = buffer.data(slk0 + 1418);
    const auto *slk0_1419 = buffer.data(slk0 + 1419);
    const auto *slk0_1421 = buffer.data(slk0 + 1421);
    const auto *slk0_1422 = buffer.data(slk0 + 1422);
    const auto *slk0_1424 = buffer.data(slk0 + 1424);
    const auto *slk0_1432 = buffer.data(slk0 + 1432);
    const auto *slk0_1434 = buffer.data(slk0 + 1434);
    const auto *slk0_1435 = buffer.data(slk0 + 1435);
    const auto *slk0_1436 = buffer.data(slk0 + 1436);
    const auto *slk0_1437 = buffer.data(slk0 + 1437);
    const auto *slk0_1439 = buffer.data(slk0 + 1439);
    const auto *slk0_1440 = buffer.data(slk0 + 1440);
    const auto *slk0_1443 = buffer.data(slk0 + 1443);
    const auto *slk0_1445 = buffer.data(slk0 + 1445);
    const auto *slk0_1446 = buffer.data(slk0 + 1446);
    const auto *slk0_1449 = buffer.data(slk0 + 1449);
    const auto *slk0_1450 = buffer.data(slk0 + 1450);
    const auto *slk0_1452 = buffer.data(slk0 + 1452);
    const auto *slk0_1454 = buffer.data(slk0 + 1454);

    const auto *sli_787 = buffer.data(sli + 787);
    const auto *sli_790 = buffer.data(sli + 790);
    const auto *sli_794 = buffer.data(sli + 794);
    const auto *sli_805 = buffer.data(sli + 805);
    const auto *sli_812 = buffer.data(sli + 812);
    const auto *sli_815 = buffer.data(sli + 815);
    const auto *sli_817 = buffer.data(sli + 817);
    const auto *sli_818 = buffer.data(sli + 818);
    const auto *sli_821 = buffer.data(sli + 821);
    const auto *sli_822 = buffer.data(sli + 822);
    const auto *sli_826 = buffer.data(sli + 826);
    const auto *sli_833 = buffer.data(sli + 833);
    const auto *sli_839 = buffer.data(sli + 839);
    const auto *sli_840 = buffer.data(sli + 840);
    const auto *sli_842 = buffer.data(sli + 842);
    const auto *sli_843 = buffer.data(sli + 843);
    const auto *sli_845 = buffer.data(sli + 845);
    const auto *sli_846 = buffer.data(sli + 846);
    const auto *sli_849 = buffer.data(sli + 849);
    const auto *sli_850 = buffer.data(sli + 850);
    const auto *sli_854 = buffer.data(sli + 854);
    const auto *sli_861 = buffer.data(sli + 861);
    const auto *sli_867 = buffer.data(sli + 867);
    const auto *sli_868 = buffer.data(sli + 868);
    const auto *sli_870 = buffer.data(sli + 870);
    const auto *sli_871 = buffer.data(sli + 871);
    const auto *sli_873 = buffer.data(sli + 873);
    const auto *sli_874 = buffer.data(sli + 874);
    const auto *sli_877 = buffer.data(sli + 877);
    const auto *sli_882 = buffer.data(sli + 882);
    const auto *sli_895 = buffer.data(sli + 895);
    const auto *sli_896 = buffer.data(sli + 896);
    const auto *sli_898 = buffer.data(sli + 898);
    const auto *sli_901 = buffer.data(sli + 901);
    const auto *sli_905 = buffer.data(sli + 905);
    const auto *sli_1045 = buffer.data(sli + 1045);
    const auto *sli_1048 = buffer.data(sli + 1048);
    const auto *sli_1050 = buffer.data(sli + 1050);
    const auto *sli_1053 = buffer.data(sli + 1053);
    const auto *sli_1054 = buffer.data(sli + 1054);
    const auto *sli_1056 = buffer.data(sli + 1056);
    const auto *sli_1057 = buffer.data(sli + 1057);
    const auto *sli_1058 = buffer.data(sli + 1058);
    const auto *sli_1059 = buffer.data(sli + 1059);
    const auto *sli_1060 = buffer.data(sli + 1060);
    const auto *sli_1061 = buffer.data(sli + 1061);
    const auto *sli_1062 = buffer.data(sli + 1062);
    const auto *sli_1063 = buffer.data(sli + 1063);
    const auto *sli_1064 = buffer.data(sli + 1064);
    const auto *sli_1067 = buffer.data(sli + 1067);
    const auto *sli_1069 = buffer.data(sli + 1069);
    const auto *sli_1070 = buffer.data(sli + 1070);
    const auto *sli_1073 = buffer.data(sli + 1073);
    const auto *sli_1074 = buffer.data(sli + 1074);
    const auto *sli_1076 = buffer.data(sli + 1076);
    const auto *sli_1078 = buffer.data(sli + 1078);
    const auto *sli_1079 = buffer.data(sli + 1079);
    const auto *sli_1081 = buffer.data(sli + 1081);
    const auto *sli_1082 = buffer.data(sli + 1082);
    const auto *sli_1084 = buffer.data(sli + 1084);
    const auto *sli_1085 = buffer.data(sli + 1085);
    const auto *sli_1086 = buffer.data(sli + 1086);
    const auto *sli_1087 = buffer.data(sli + 1087);
    const auto *sli_1088 = buffer.data(sli + 1088);
    const auto *sli_1089 = buffer.data(sli + 1089);
    const auto *sli_1090 = buffer.data(sli + 1090);
    const auto *sli_1091 = buffer.data(sli + 1091);
    const auto *sli_1092 = buffer.data(sli + 1092);
    const auto *sli_1095 = buffer.data(sli + 1095);
    const auto *sli_1097 = buffer.data(sli + 1097);
    const auto *sli_1098 = buffer.data(sli + 1098);
    const auto *sli_1101 = buffer.data(sli + 1101);
    const auto *sli_1102 = buffer.data(sli + 1102);
    const auto *sli_1104 = buffer.data(sli + 1104);
    const auto *sli_1106 = buffer.data(sli + 1106);
    const auto *sli_1107 = buffer.data(sli + 1107);
    const auto *sli_1109 = buffer.data(sli + 1109);
    const auto *sli_1110 = buffer.data(sli + 1110);
    const auto *sli_1112 = buffer.data(sli + 1112);
    const auto *sli_1113 = buffer.data(sli + 1113);
    const auto *sli_1114 = buffer.data(sli + 1114);
    const auto *sli_1115 = buffer.data(sli + 1115);
    const auto *sli_1116 = buffer.data(sli + 1116);
    const auto *sli_1117 = buffer.data(sli + 1117);
    const auto *sli_1118 = buffer.data(sli + 1118);
    const auto *sli_1119 = buffer.data(sli + 1119);
    const auto *sli_1120 = buffer.data(sli + 1120);
    const auto *sli_1123 = buffer.data(sli + 1123);
    const auto *sli_1125 = buffer.data(sli + 1125);
    const auto *sli_1126 = buffer.data(sli + 1126);
    const auto *sli_1129 = buffer.data(sli + 1129);
    const auto *sli_1130 = buffer.data(sli + 1130);
    const auto *sli_1132 = buffer.data(sli + 1132);
    const auto *sli_1134 = buffer.data(sli + 1134);

    const auto *slk1_1018 = buffer.data(slk1 + 1018);
    const auto *slk1_1023 = buffer.data(slk1 + 1023);
    const auto *slk1_1341 = buffer.data(slk1 + 1341);
    const auto *slk1_1344 = buffer.data(slk1 + 1344);
    const auto *slk1_1346 = buffer.data(slk1 + 1346);
    const auto *slk1_1349 = buffer.data(slk1 + 1349);
    const auto *slk1_1350 = buffer.data(slk1 + 1350);
    const auto *slk1_1352 = buffer.data(slk1 + 1352);
    const auto *slk1_1360 = buffer.data(slk1 + 1360);
    const auto *slk1_1362 = buffer.data(slk1 + 1362);
    const auto *slk1_1363 = buffer.data(slk1 + 1363);
    const auto *slk1_1364 = buffer.data(slk1 + 1364);
    const auto *slk1_1365 = buffer.data(slk1 + 1365);
    const auto *slk1_1367 = buffer.data(slk1 + 1367);
    const auto *slk1_1368 = buffer.data(slk1 + 1368);
    const auto *slk1_1371 = buffer.data(slk1 + 1371);
    const auto *slk1_1373 = buffer.data(slk1 + 1373);
    const auto *slk1_1374 = buffer.data(slk1 + 1374);
    const auto *slk1_1377 = buffer.data(slk1 + 1377);
    const auto *slk1_1378 = buffer.data(slk1 + 1378);
    const auto *slk1_1380 = buffer.data(slk1 + 1380);
    const auto *slk1_1382 = buffer.data(slk1 + 1382);
    const auto *slk1_1383 = buffer.data(slk1 + 1383);
    const auto *slk1_1385 = buffer.data(slk1 + 1385);
    const auto *slk1_1386 = buffer.data(slk1 + 1386);
    const auto *slk1_1388 = buffer.data(slk1 + 1388);
    const auto *slk1_1396 = buffer.data(slk1 + 1396);
    const auto *slk1_1398 = buffer.data(slk1 + 1398);
    const auto *slk1_1399 = buffer.data(slk1 + 1399);
    const auto *slk1_1400 = buffer.data(slk1 + 1400);
    const auto *slk1_1401 = buffer.data(slk1 + 1401);
    const auto *slk1_1403 = buffer.data(slk1 + 1403);
    const auto *slk1_1404 = buffer.data(slk1 + 1404);
    const auto *slk1_1407 = buffer.data(slk1 + 1407);
    const auto *slk1_1409 = buffer.data(slk1 + 1409);
    const auto *slk1_1410 = buffer.data(slk1 + 1410);
    const auto *slk1_1413 = buffer.data(slk1 + 1413);
    const auto *slk1_1414 = buffer.data(slk1 + 1414);
    const auto *slk1_1416 = buffer.data(slk1 + 1416);
    const auto *slk1_1418 = buffer.data(slk1 + 1418);
    const auto *slk1_1419 = buffer.data(slk1 + 1419);
    const auto *slk1_1421 = buffer.data(slk1 + 1421);
    const auto *slk1_1422 = buffer.data(slk1 + 1422);
    const auto *slk1_1424 = buffer.data(slk1 + 1424);
    const auto *slk1_1432 = buffer.data(slk1 + 1432);
    const auto *slk1_1434 = buffer.data(slk1 + 1434);
    const auto *slk1_1435 = buffer.data(slk1 + 1435);
    const auto *slk1_1436 = buffer.data(slk1 + 1436);
    const auto *slk1_1437 = buffer.data(slk1 + 1437);
    const auto *slk1_1439 = buffer.data(slk1 + 1439);
    const auto *slk1_1440 = buffer.data(slk1 + 1440);
    const auto *slk1_1443 = buffer.data(slk1 + 1443);
    const auto *slk1_1445 = buffer.data(slk1 + 1445);
    const auto *slk1_1446 = buffer.data(slk1 + 1446);
    const auto *slk1_1449 = buffer.data(slk1 + 1449);
    const auto *slk1_1450 = buffer.data(slk1 + 1450);
    const auto *slk1_1452 = buffer.data(slk1 + 1452);
    const auto *slk1_1454 = buffer.data(slk1 + 1454);

    const auto *smi_1039 = buffer.data(smi + 1039);
    const auto *smi_1041 = buffer.data(smi + 1041);
    const auto *smi_1042 = buffer.data(smi + 1042);
    const auto *smi_1045 = buffer.data(smi + 1045);
    const auto *smi_1046 = buffer.data(smi + 1046);
    const auto *smi_1050 = buffer.data(smi + 1050);
    const auto *smi_1057 = buffer.data(smi + 1057);
    const auto *smi_1058 = buffer.data(smi + 1058);
    const auto *smi_1059 = buffer.data(smi + 1059);
    const auto *smi_1060 = buffer.data(smi + 1060);
    const auto *smi_1061 = buffer.data(smi + 1061);
    const auto *smi_1062 = buffer.data(smi + 1062);
    const auto *smi_1063 = buffer.data(smi + 1063);
    const auto *smi_1064 = buffer.data(smi + 1064);
    const auto *smi_1066 = buffer.data(smi + 1066);
    const auto *smi_1067 = buffer.data(smi + 1067);
    const auto *smi_1069 = buffer.data(smi + 1069);
    const auto *smi_1070 = buffer.data(smi + 1070);
    const auto *smi_1073 = buffer.data(smi + 1073);
    const auto *smi_1074 = buffer.data(smi + 1074);
    const auto *smi_1078 = buffer.data(smi + 1078);
    const auto *smi_1085 = buffer.data(smi + 1085);
    const auto *smi_1086 = buffer.data(smi + 1086);
    const auto *smi_1087 = buffer.data(smi + 1087);
    const auto *smi_1088 = buffer.data(smi + 1088);
    const auto *smi_1089 = buffer.data(smi + 1089);
    const auto *smi_1090 = buffer.data(smi + 1090);
    const auto *smi_1091 = buffer.data(smi + 1091);
    const auto *smi_1092 = buffer.data(smi + 1092);
    const auto *smi_1094 = buffer.data(smi + 1094);
    const auto *smi_1095 = buffer.data(smi + 1095);
    const auto *smi_1097 = buffer.data(smi + 1097);
    const auto *smi_1098 = buffer.data(smi + 1098);
    const auto *smi_1101 = buffer.data(smi + 1101);
    const auto *smi_1102 = buffer.data(smi + 1102);
    const auto *smi_1106 = buffer.data(smi + 1106);
    const auto *smi_1113 = buffer.data(smi + 1113);
    const auto *smi_1114 = buffer.data(smi + 1114);
    const auto *smi_1115 = buffer.data(smi + 1115);
    const auto *smi_1116 = buffer.data(smi + 1116);
    const auto *smi_1117 = buffer.data(smi + 1117);
    const auto *smi_1118 = buffer.data(smi + 1118);
    const auto *smi_1119 = buffer.data(smi + 1119);
    const auto *smi_1120 = buffer.data(smi + 1120);
    const auto *smi_1122 = buffer.data(smi + 1122);
    const auto *smi_1123 = buffer.data(smi + 1123);
    const auto *smi_1125 = buffer.data(smi + 1125);
    const auto *smi_1126 = buffer.data(smi + 1126);
    const auto *smi_1129 = buffer.data(smi + 1129);

#pragma omp simd aligned(t_1339, t_1340, t_1341, pb_x, pc_x, pc_y, pc_z, slk0_1341, sli_787, \
                         sli_817, sli_1045, slk1_1341, smi_1039, \
                         smi_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1339[k] = f_13 * sli_787[k]
                    + f_3 * pc_z[k] * smi_1039[k];

        t_1340[k] = f_19 * sli_817[k]
                    + f_3 * pc_y[k] * smi_1041[k];

        t_1341[k] = pb_x[k] * slk0_1341[k]
                    + f_16 * sli_1045[k]
                    - f_12 * pc_x[k] * slk1_1341[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, pb_x, pb_z, pc_x, pc_z, slk0_1018, slk0_1344, \
                         sli_790, sli_1048, slk1_1018, slk1_1344, \
                         smi_1042 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = pb_z[k] * slk0_1018[k]
                    - f_12 * pc_z[k] * slk1_1018[k];

        t_1343[k] = f_13 * sli_790[k]
                    + f_3 * pc_z[k] * smi_1042[k];

        t_1344[k] = pb_x[k] * slk0_1344[k]
                    + f_15 * sli_1048[k]
                    - f_12 * pc_x[k] * slk1_1344[k];
    }

#pragma omp simd aligned(t_1345, t_1346, t_1347, pb_x, pb_z, pc_x, pc_y, pc_z, slk0_1023, \
                         slk0_1346, sli_821, sli_1050, slk1_1023, slk1_1346, \
                         smi_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1345[k] = f_19 * sli_821[k]
                    + f_3 * pc_y[k] * smi_1045[k];

        t_1346[k] = pb_x[k] * slk0_1346[k]
                    + f_15 * sli_1050[k]
                    - f_12 * pc_x[k] * slk1_1346[k];

        t_1347[k] = pb_z[k] * slk0_1023[k]
                    - f_12 * pc_z[k] * slk1_1023[k];
    }

#pragma omp simd aligned(t_1348, t_1349, t_1350, pb_x, pc_x, pc_z, slk0_1349, slk0_1350, \
                         sli_794, sli_1053, sli_1054, slk1_1349, slk1_1350, \
                         smi_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1348[k] = f_13 * sli_794[k]
                    + f_3 * pc_z[k] * smi_1046[k];

        t_1349[k] = pb_x[k] * slk0_1349[k]
                    + f_14 * sli_1053[k]
                    - f_12 * pc_x[k] * slk1_1349[k];

        t_1350[k] = pb_x[k] * slk0_1350[k]
                    + f_14 * sli_1054[k]
                    - f_12 * pc_x[k] * slk1_1350[k];
    }

#pragma omp simd aligned(t_1351, t_1352, t_1353, t_1354, pb_x, pc_x, pc_y, slk0_1352, sli_826, \
                         sli_1056, sli_1057, sli_1058, slk1_1352, smi_1050, smi_1057, \
                         smi_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1351[k] = f_19 * sli_826[k]
                    + f_3 * pc_y[k] * smi_1050[k];

        t_1352[k] = pb_x[k] * slk0_1352[k]
                    + f_14 * sli_1056[k]
                    - f_12 * pc_x[k] * slk1_1352[k];

        t_1353[k] = f_13 * sli_1057[k]
                    + f_3 * pc_x[k] * smi_1057[k];

        t_1354[k] = f_13 * sli_1058[k]
                    + f_3 * pc_x[k] * smi_1058[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, t_1358, t_1359, pc_x, sli_1059, sli_1060, \
                         sli_1061, sli_1062, sli_1063, smi_1059, smi_1060, smi_1061, smi_1062, \
                         smi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_13 * sli_1059[k]
                    + f_3 * pc_x[k] * smi_1059[k];

        t_1356[k] = f_13 * sli_1060[k]
                    + f_3 * pc_x[k] * smi_1060[k];

        t_1357[k] = f_13 * sli_1061[k]
                    + f_3 * pc_x[k] * smi_1061[k];

        t_1358[k] = f_13 * sli_1062[k]
                    + f_3 * pc_x[k] * smi_1062[k];

        t_1359[k] = f_13 * sli_1063[k]
                    + f_3 * pc_x[k] * smi_1063[k];
    }

#pragma omp simd aligned(t_1360, t_1361, t_1362, t_1363, pb_x, pc_x, pc_z, slk0_1360, \
                         slk0_1362, slk0_1363, sli_805, slk1_1360, slk1_1362, slk1_1363, \
                         smi_1057 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1360[k] = pb_x[k] * slk0_1360[k]
                    - f_12 * pc_x[k] * slk1_1360[k];

        t_1361[k] = f_13 * sli_805[k]
                    + f_3 * pc_z[k] * smi_1057[k];

        t_1362[k] = pb_x[k] * slk0_1362[k]
                    - f_12 * pc_x[k] * slk1_1362[k];

        t_1363[k] = pb_x[k] * slk0_1363[k]
                    - f_12 * pc_x[k] * slk1_1363[k];
    }

#pragma omp simd aligned(t_1364, t_1365, t_1366, t_1367, pb_x, pc_x, pc_y, slk0_1364, \
                         slk0_1365, slk0_1367, sli_839, slk1_1364, slk1_1365, slk1_1367, \
                         smi_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1364[k] = pb_x[k] * slk0_1364[k]
                    - f_12 * pc_x[k] * slk1_1364[k];

        t_1365[k] = pb_x[k] * slk0_1365[k]
                    - f_12 * pc_x[k] * slk1_1365[k];

        t_1366[k] = f_19 * sli_839[k]
                    + f_3 * pc_y[k] * smi_1063[k];

        t_1367[k] = pb_x[k] * slk0_1367[k]
                    - f_12 * pc_x[k] * slk1_1367[k];
    }

#pragma omp simd aligned(t_1368, t_1369, t_1370, pb_x, pc_x, pc_y, pc_z, slk0_1368, sli_812, \
                         sli_840, sli_1064, slk1_1368, smi_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1368[k] = pb_x[k] * slk0_1368[k]
                    + f_19 * sli_1064[k]
                    - f_12 * pc_x[k] * slk1_1368[k];

        t_1369[k] = f_20 * sli_840[k]
                    + f_3 * pc_y[k] * smi_1064[k];

        t_1370[k] = f_14 * sli_812[k]
                    + f_3 * pc_z[k] * smi_1064[k];
    }

#pragma omp simd aligned(t_1371, t_1372, t_1373, pb_x, pc_x, pc_y, slk0_1371, slk0_1373, \
                         sli_842, sli_1067, sli_1069, slk1_1371, slk1_1373, \
                         smi_1066 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1371[k] = pb_x[k] * slk0_1371[k]
                    + f_17 * sli_1067[k]
                    - f_12 * pc_x[k] * slk1_1371[k];

        t_1372[k] = f_20 * sli_842[k]
                    + f_3 * pc_y[k] * smi_1066[k];

        t_1373[k] = pb_x[k] * slk0_1373[k]
                    + f_17 * sli_1069[k]
                    - f_12 * pc_x[k] * slk1_1373[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, pb_x, pc_x, pc_y, pc_z, slk0_1374, sli_815, \
                         sli_845, sli_1070, slk1_1374, smi_1067, \
                         smi_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = pb_x[k] * slk0_1374[k]
                    + f_16 * sli_1070[k]
                    - f_12 * pc_x[k] * slk1_1374[k];

        t_1375[k] = f_14 * sli_815[k]
                    + f_3 * pc_z[k] * smi_1067[k];

        t_1376[k] = f_20 * sli_845[k]
                    + f_3 * pc_y[k] * smi_1069[k];
    }

#pragma omp simd aligned(t_1377, t_1378, t_1379, pb_x, pc_x, pc_z, slk0_1377, slk0_1378, \
                         sli_818, sli_1073, sli_1074, slk1_1377, slk1_1378, \
                         smi_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1377[k] = pb_x[k] * slk0_1377[k]
                    + f_16 * sli_1073[k]
                    - f_12 * pc_x[k] * slk1_1377[k];

        t_1378[k] = pb_x[k] * slk0_1378[k]
                    + f_15 * sli_1074[k]
                    - f_12 * pc_x[k] * slk1_1378[k];

        t_1379[k] = f_14 * sli_818[k]
                    + f_3 * pc_z[k] * smi_1070[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, pb_x, pc_x, pc_y, slk0_1380, slk0_1382, \
                         sli_849, sli_1076, sli_1078, slk1_1380, slk1_1382, \
                         smi_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = pb_x[k] * slk0_1380[k]
                    + f_15 * sli_1076[k]
                    - f_12 * pc_x[k] * slk1_1380[k];

        t_1381[k] = f_20 * sli_849[k]
                    + f_3 * pc_y[k] * smi_1073[k];

        t_1382[k] = pb_x[k] * slk0_1382[k]
                    + f_15 * sli_1078[k]
                    - f_12 * pc_x[k] * slk1_1382[k];
    }

#pragma omp simd aligned(t_1383, t_1384, t_1385, pb_x, pc_x, pc_z, slk0_1383, slk0_1385, \
                         sli_822, sli_1079, sli_1081, slk1_1383, slk1_1385, \
                         smi_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1383[k] = pb_x[k] * slk0_1383[k]
                    + f_14 * sli_1079[k]
                    - f_12 * pc_x[k] * slk1_1383[k];

        t_1384[k] = f_14 * sli_822[k]
                    + f_3 * pc_z[k] * smi_1074[k];

        t_1385[k] = pb_x[k] * slk0_1385[k]
                    + f_14 * sli_1081[k]
                    - f_12 * pc_x[k] * slk1_1385[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, pb_x, pc_x, pc_y, slk0_1386, slk0_1388, \
                         sli_854, sli_1082, sli_1084, slk1_1386, slk1_1388, \
                         smi_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = pb_x[k] * slk0_1386[k]
                    + f_14 * sli_1082[k]
                    - f_12 * pc_x[k] * slk1_1386[k];

        t_1387[k] = f_20 * sli_854[k]
                    + f_3 * pc_y[k] * smi_1078[k];

        t_1388[k] = pb_x[k] * slk0_1388[k]
                    + f_14 * sli_1084[k]
                    - f_12 * pc_x[k] * slk1_1388[k];
    }

#pragma omp simd aligned(t_1389, t_1390, t_1391, t_1392, t_1393, pc_x, sli_1085, sli_1086, \
                         sli_1087, sli_1088, sli_1089, smi_1085, smi_1086, smi_1087, smi_1088, \
                         smi_1089 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1389[k] = f_13 * sli_1085[k]
                    + f_3 * pc_x[k] * smi_1085[k];

        t_1390[k] = f_13 * sli_1086[k]
                    + f_3 * pc_x[k] * smi_1086[k];

        t_1391[k] = f_13 * sli_1087[k]
                    + f_3 * pc_x[k] * smi_1087[k];

        t_1392[k] = f_13 * sli_1088[k]
                    + f_3 * pc_x[k] * smi_1088[k];

        t_1393[k] = f_13 * sli_1089[k]
                    + f_3 * pc_x[k] * smi_1089[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, t_1397, pb_x, pc_x, pc_z, slk0_1396, sli_833, \
                         sli_1090, sli_1091, slk1_1396, smi_1085, smi_1090, \
                         smi_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_13 * sli_1090[k]
                    + f_3 * pc_x[k] * smi_1090[k];

        t_1395[k] = f_13 * sli_1091[k]
                    + f_3 * pc_x[k] * smi_1091[k];

        t_1396[k] = pb_x[k] * slk0_1396[k]
                    - f_12 * pc_x[k] * slk1_1396[k];

        t_1397[k] = f_14 * sli_833[k]
                    + f_3 * pc_z[k] * smi_1085[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, t_1401, pb_x, pc_x, slk0_1398, slk0_1399, \
                         slk0_1400, slk0_1401, slk1_1398, slk1_1399, slk1_1400, \
                         slk1_1401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = pb_x[k] * slk0_1398[k]
                    - f_12 * pc_x[k] * slk1_1398[k];

        t_1399[k] = pb_x[k] * slk0_1399[k]
                    - f_12 * pc_x[k] * slk1_1399[k];

        t_1400[k] = pb_x[k] * slk0_1400[k]
                    - f_12 * pc_x[k] * slk1_1400[k];

        t_1401[k] = pb_x[k] * slk0_1401[k]
                    - f_12 * pc_x[k] * slk1_1401[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, t_1405, pb_x, pc_x, pc_y, slk0_1403, \
                         slk0_1404, sli_867, sli_868, sli_1092, slk1_1403, slk1_1404, \
                         smi_1091, smi_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_20 * sli_867[k]
                    + f_3 * pc_y[k] * smi_1091[k];

        t_1403[k] = pb_x[k] * slk0_1403[k]
                    - f_12 * pc_x[k] * slk1_1403[k];

        t_1404[k] = pb_x[k] * slk0_1404[k]
                    + f_19 * sli_1092[k]
                    - f_12 * pc_x[k] * slk1_1404[k];

        t_1405[k] = f_17 * sli_868[k]
                    + f_3 * pc_y[k] * smi_1092[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, pb_x, pc_x, pc_y, pc_z, slk0_1407, sli_840, \
                         sli_870, sli_1095, slk1_1407, smi_1092, \
                         smi_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_15 * sli_840[k]
                    + f_3 * pc_z[k] * smi_1092[k];

        t_1407[k] = pb_x[k] * slk0_1407[k]
                    + f_17 * sli_1095[k]
                    - f_12 * pc_x[k] * slk1_1407[k];

        t_1408[k] = f_17 * sli_870[k]
                    + f_3 * pc_y[k] * smi_1094[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pb_x, pc_x, pc_z, slk0_1409, slk0_1410, \
                         sli_843, sli_1097, sli_1098, slk1_1409, slk1_1410, \
                         smi_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = pb_x[k] * slk0_1409[k]
                    + f_17 * sli_1097[k]
                    - f_12 * pc_x[k] * slk1_1409[k];

        t_1410[k] = pb_x[k] * slk0_1410[k]
                    + f_16 * sli_1098[k]
                    - f_12 * pc_x[k] * slk1_1410[k];

        t_1411[k] = f_15 * sli_843[k]
                    + f_3 * pc_z[k] * smi_1095[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, pb_x, pc_x, pc_y, slk0_1413, slk0_1414, \
                         sli_873, sli_1101, sli_1102, slk1_1413, slk1_1414, \
                         smi_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_17 * sli_873[k]
                    + f_3 * pc_y[k] * smi_1097[k];

        t_1413[k] = pb_x[k] * slk0_1413[k]
                    + f_16 * sli_1101[k]
                    - f_12 * pc_x[k] * slk1_1413[k];

        t_1414[k] = pb_x[k] * slk0_1414[k]
                    + f_15 * sli_1102[k]
                    - f_12 * pc_x[k] * slk1_1414[k];
    }

#pragma omp simd aligned(t_1415, t_1416, t_1417, pb_x, pc_x, pc_y, pc_z, slk0_1416, sli_846, \
                         sli_877, sli_1104, slk1_1416, smi_1098, \
                         smi_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1415[k] = f_15 * sli_846[k]
                    + f_3 * pc_z[k] * smi_1098[k];

        t_1416[k] = pb_x[k] * slk0_1416[k]
                    + f_15 * sli_1104[k]
                    - f_12 * pc_x[k] * slk1_1416[k];

        t_1417[k] = f_17 * sli_877[k]
                    + f_3 * pc_y[k] * smi_1101[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, pb_x, pc_x, pc_z, slk0_1418, slk0_1419, \
                         sli_850, sli_1106, sli_1107, slk1_1418, slk1_1419, \
                         smi_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = pb_x[k] * slk0_1418[k]
                    + f_15 * sli_1106[k]
                    - f_12 * pc_x[k] * slk1_1418[k];

        t_1419[k] = pb_x[k] * slk0_1419[k]
                    + f_14 * sli_1107[k]
                    - f_12 * pc_x[k] * slk1_1419[k];

        t_1420[k] = f_15 * sli_850[k]
                    + f_3 * pc_z[k] * smi_1102[k];
    }

#pragma omp simd aligned(t_1421, t_1422, t_1423, pb_x, pc_x, pc_y, slk0_1421, slk0_1422, \
                         sli_882, sli_1109, sli_1110, slk1_1421, slk1_1422, \
                         smi_1106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1421[k] = pb_x[k] * slk0_1421[k]
                    + f_14 * sli_1109[k]
                    - f_12 * pc_x[k] * slk1_1421[k];

        t_1422[k] = pb_x[k] * slk0_1422[k]
                    + f_14 * sli_1110[k]
                    - f_12 * pc_x[k] * slk1_1422[k];

        t_1423[k] = f_17 * sli_882[k]
                    + f_3 * pc_y[k] * smi_1106[k];
    }

#pragma omp simd aligned(t_1424, t_1425, t_1426, t_1427, pb_x, pc_x, slk0_1424, sli_1112, \
                         sli_1113, sli_1114, sli_1115, slk1_1424, smi_1113, smi_1114, \
                         smi_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1424[k] = pb_x[k] * slk0_1424[k]
                    + f_14 * sli_1112[k]
                    - f_12 * pc_x[k] * slk1_1424[k];

        t_1425[k] = f_13 * sli_1113[k]
                    + f_3 * pc_x[k] * smi_1113[k];

        t_1426[k] = f_13 * sli_1114[k]
                    + f_3 * pc_x[k] * smi_1114[k];

        t_1427[k] = f_13 * sli_1115[k]
                    + f_3 * pc_x[k] * smi_1115[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, t_1431, pc_x, sli_1116, sli_1117, sli_1118, \
                         sli_1119, smi_1116, smi_1117, smi_1118, \
                         smi_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_13 * sli_1116[k]
                    + f_3 * pc_x[k] * smi_1116[k];

        t_1429[k] = f_13 * sli_1117[k]
                    + f_3 * pc_x[k] * smi_1117[k];

        t_1430[k] = f_13 * sli_1118[k]
                    + f_3 * pc_x[k] * smi_1118[k];

        t_1431[k] = f_13 * sli_1119[k]
                    + f_3 * pc_x[k] * smi_1119[k];
    }

#pragma omp simd aligned(t_1432, t_1433, t_1434, t_1435, pb_x, pc_x, pc_z, slk0_1432, \
                         slk0_1434, slk0_1435, sli_861, slk1_1432, slk1_1434, slk1_1435, \
                         smi_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1432[k] = pb_x[k] * slk0_1432[k]
                    - f_12 * pc_x[k] * slk1_1432[k];

        t_1433[k] = f_15 * sli_861[k]
                    + f_3 * pc_z[k] * smi_1113[k];

        t_1434[k] = pb_x[k] * slk0_1434[k]
                    - f_12 * pc_x[k] * slk1_1434[k];

        t_1435[k] = pb_x[k] * slk0_1435[k]
                    - f_12 * pc_x[k] * slk1_1435[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, t_1439, pb_x, pc_x, pc_y, slk0_1436, \
                         slk0_1437, slk0_1439, sli_895, slk1_1436, slk1_1437, slk1_1439, \
                         smi_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = pb_x[k] * slk0_1436[k]
                    - f_12 * pc_x[k] * slk1_1436[k];

        t_1437[k] = pb_x[k] * slk0_1437[k]
                    - f_12 * pc_x[k] * slk1_1437[k];

        t_1438[k] = f_17 * sli_895[k]
                    + f_3 * pc_y[k] * smi_1119[k];

        t_1439[k] = pb_x[k] * slk0_1439[k]
                    - f_12 * pc_x[k] * slk1_1439[k];
    }

#pragma omp simd aligned(t_1440, t_1441, t_1442, pb_x, pc_x, pc_y, pc_z, slk0_1440, sli_868, \
                         sli_896, sli_1120, slk1_1440, smi_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1440[k] = pb_x[k] * slk0_1440[k]
                    + f_19 * sli_1120[k]
                    - f_12 * pc_x[k] * slk1_1440[k];

        t_1441[k] = f_16 * sli_896[k]
                    + f_3 * pc_y[k] * smi_1120[k];

        t_1442[k] = f_16 * sli_868[k]
                    + f_3 * pc_z[k] * smi_1120[k];
    }

#pragma omp simd aligned(t_1443, t_1444, t_1445, pb_x, pc_x, pc_y, slk0_1443, slk0_1445, \
                         sli_898, sli_1123, sli_1125, slk1_1443, slk1_1445, \
                         smi_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1443[k] = pb_x[k] * slk0_1443[k]
                    + f_17 * sli_1123[k]
                    - f_12 * pc_x[k] * slk1_1443[k];

        t_1444[k] = f_16 * sli_898[k]
                    + f_3 * pc_y[k] * smi_1122[k];

        t_1445[k] = pb_x[k] * slk0_1445[k]
                    + f_17 * sli_1125[k]
                    - f_12 * pc_x[k] * slk1_1445[k];
    }

#pragma omp simd aligned(t_1446, t_1447, t_1448, pb_x, pc_x, pc_y, pc_z, slk0_1446, sli_871, \
                         sli_901, sli_1126, slk1_1446, smi_1123, \
                         smi_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1446[k] = pb_x[k] * slk0_1446[k]
                    + f_16 * sli_1126[k]
                    - f_12 * pc_x[k] * slk1_1446[k];

        t_1447[k] = f_16 * sli_871[k]
                    + f_3 * pc_z[k] * smi_1123[k];

        t_1448[k] = f_16 * sli_901[k]
                    + f_3 * pc_y[k] * smi_1125[k];
    }

#pragma omp simd aligned(t_1449, t_1450, t_1451, pb_x, pc_x, pc_z, slk0_1449, slk0_1450, \
                         sli_874, sli_1129, sli_1130, slk1_1449, slk1_1450, \
                         smi_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1449[k] = pb_x[k] * slk0_1449[k]
                    + f_16 * sli_1129[k]
                    - f_12 * pc_x[k] * slk1_1449[k];

        t_1450[k] = pb_x[k] * slk0_1450[k]
                    + f_15 * sli_1130[k]
                    - f_12 * pc_x[k] * slk1_1450[k];

        t_1451[k] = f_16 * sli_874[k]
                    + f_3 * pc_z[k] * smi_1126[k];
    }

#pragma omp simd aligned(t_1452, t_1453, t_1454, pb_x, pc_x, pc_y, slk0_1452, slk0_1454, \
                         sli_905, sli_1132, sli_1134, slk1_1452, slk1_1454, \
                         smi_1129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1452[k] = pb_x[k] * slk0_1452[k]
                    + f_15 * sli_1132[k]
                    - f_12 * pc_x[k] * slk1_1452[k];

        t_1453[k] = f_16 * sli_905[k]
                    + f_3 * pc_y[k] * smi_1129[k];

        t_1454[k] = pb_x[k] * slk0_1454[k]
                    + f_15 * sli_1134[k]
                    - f_12 * pc_x[k] * slk1_1454[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t slk0,
                                                           const size_t sli, const size_t slk1,
                                                           const size_t smi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_1260 = buffer.data(slk0 + 1260);
    const auto *slk0_1265 = buffer.data(slk0 + 1265);
    const auto *slk0_1269 = buffer.data(slk0 + 1269);
    const auto *slk0_1274 = buffer.data(slk0 + 1274);
    const auto *slk0_1280 = buffer.data(slk0 + 1280);
    const auto *slk0_1455 = buffer.data(slk0 + 1455);
    const auto *slk0_1457 = buffer.data(slk0 + 1457);
    const auto *slk0_1458 = buffer.data(slk0 + 1458);
    const auto *slk0_1460 = buffer.data(slk0 + 1460);
    const auto *slk0_1468 = buffer.data(slk0 + 1468);
    const auto *slk0_1470 = buffer.data(slk0 + 1470);
    const auto *slk0_1471 = buffer.data(slk0 + 1471);
    const auto *slk0_1472 = buffer.data(slk0 + 1472);
    const auto *slk0_1473 = buffer.data(slk0 + 1473);
    const auto *slk0_1475 = buffer.data(slk0 + 1475);
    const auto *slk0_1476 = buffer.data(slk0 + 1476);
    const auto *slk0_1479 = buffer.data(slk0 + 1479);
    const auto *slk0_1481 = buffer.data(slk0 + 1481);
    const auto *slk0_1482 = buffer.data(slk0 + 1482);
    const auto *slk0_1485 = buffer.data(slk0 + 1485);
    const auto *slk0_1486 = buffer.data(slk0 + 1486);
    const auto *slk0_1488 = buffer.data(slk0 + 1488);
    const auto *slk0_1490 = buffer.data(slk0 + 1490);
    const auto *slk0_1491 = buffer.data(slk0 + 1491);
    const auto *slk0_1493 = buffer.data(slk0 + 1493);
    const auto *slk0_1494 = buffer.data(slk0 + 1494);
    const auto *slk0_1496 = buffer.data(slk0 + 1496);
    const auto *slk0_1504 = buffer.data(slk0 + 1504);
    const auto *slk0_1506 = buffer.data(slk0 + 1506);
    const auto *slk0_1507 = buffer.data(slk0 + 1507);
    const auto *slk0_1508 = buffer.data(slk0 + 1508);
    const auto *slk0_1509 = buffer.data(slk0 + 1509);
    const auto *slk0_1511 = buffer.data(slk0 + 1511);
    const auto *slk0_1512 = buffer.data(slk0 + 1512);
    const auto *slk0_1515 = buffer.data(slk0 + 1515);
    const auto *slk0_1517 = buffer.data(slk0 + 1517);
    const auto *slk0_1518 = buffer.data(slk0 + 1518);
    const auto *slk0_1521 = buffer.data(slk0 + 1521);
    const auto *slk0_1522 = buffer.data(slk0 + 1522);
    const auto *slk0_1524 = buffer.data(slk0 + 1524);
    const auto *slk0_1526 = buffer.data(slk0 + 1526);
    const auto *slk0_1527 = buffer.data(slk0 + 1527);
    const auto *slk0_1529 = buffer.data(slk0 + 1529);
    const auto *slk0_1530 = buffer.data(slk0 + 1530);
    const auto *slk0_1532 = buffer.data(slk0 + 1532);
    const auto *slk0_1540 = buffer.data(slk0 + 1540);
    const auto *slk0_1542 = buffer.data(slk0 + 1542);
    const auto *slk0_1543 = buffer.data(slk0 + 1543);
    const auto *slk0_1544 = buffer.data(slk0 + 1544);
    const auto *slk0_1545 = buffer.data(slk0 + 1545);
    const auto *slk0_1547 = buffer.data(slk0 + 1547);
    const auto *slk0_1551 = buffer.data(slk0 + 1551);
    const auto *slk0_1554 = buffer.data(slk0 + 1554);
    const auto *slk0_1558 = buffer.data(slk0 + 1558);
    const auto *slk0_1560 = buffer.data(slk0 + 1560);
    const auto *slk0_1563 = buffer.data(slk0 + 1563);
    const auto *slk0_1565 = buffer.data(slk0 + 1565);
    const auto *slk0_1566 = buffer.data(slk0 + 1566);

    const auto *sli_878 = buffer.data(sli + 878);
    const auto *sli_889 = buffer.data(sli + 889);
    const auto *sli_896 = buffer.data(sli + 896);
    const auto *sli_899 = buffer.data(sli + 899);
    const auto *sli_902 = buffer.data(sli + 902);
    const auto *sli_906 = buffer.data(sli + 906);
    const auto *sli_910 = buffer.data(sli + 910);
    const auto *sli_917 = buffer.data(sli + 917);
    const auto *sli_923 = buffer.data(sli + 923);
    const auto *sli_924 = buffer.data(sli + 924);
    const auto *sli_926 = buffer.data(sli + 926);
    const auto *sli_927 = buffer.data(sli + 927);
    const auto *sli_929 = buffer.data(sli + 929);
    const auto *sli_930 = buffer.data(sli + 930);
    const auto *sli_933 = buffer.data(sli + 933);
    const auto *sli_934 = buffer.data(sli + 934);
    const auto *sli_938 = buffer.data(sli + 938);
    const auto *sli_945 = buffer.data(sli + 945);
    const auto *sli_951 = buffer.data(sli + 951);
    const auto *sli_952 = buffer.data(sli + 952);
    const auto *sli_954 = buffer.data(sli + 954);
    const auto *sli_955 = buffer.data(sli + 955);
    const auto *sli_957 = buffer.data(sli + 957);
    const auto *sli_958 = buffer.data(sli + 958);
    const auto *sli_961 = buffer.data(sli + 961);
    const auto *sli_962 = buffer.data(sli + 962);
    const auto *sli_966 = buffer.data(sli + 966);
    const auto *sli_979 = buffer.data(sli + 979);
    const auto *sli_980 = buffer.data(sli + 980);
    const auto *sli_982 = buffer.data(sli + 982);
    const auto *sli_985 = buffer.data(sli + 985);
    const auto *sli_989 = buffer.data(sli + 989);
    const auto *sli_994 = buffer.data(sli + 994);
    const auto *sli_1135 = buffer.data(sli + 1135);
    const auto *sli_1137 = buffer.data(sli + 1137);
    const auto *sli_1138 = buffer.data(sli + 1138);
    const auto *sli_1140 = buffer.data(sli + 1140);
    const auto *sli_1141 = buffer.data(sli + 1141);
    const auto *sli_1142 = buffer.data(sli + 1142);
    const auto *sli_1143 = buffer.data(sli + 1143);
    const auto *sli_1144 = buffer.data(sli + 1144);
    const auto *sli_1145 = buffer.data(sli + 1145);
    const auto *sli_1146 = buffer.data(sli + 1146);
    const auto *sli_1147 = buffer.data(sli + 1147);
    const auto *sli_1148 = buffer.data(sli + 1148);
    const auto *sli_1151 = buffer.data(sli + 1151);
    const auto *sli_1153 = buffer.data(sli + 1153);
    const auto *sli_1154 = buffer.data(sli + 1154);
    const auto *sli_1157 = buffer.data(sli + 1157);
    const auto *sli_1158 = buffer.data(sli + 1158);
    const auto *sli_1160 = buffer.data(sli + 1160);
    const auto *sli_1162 = buffer.data(sli + 1162);
    const auto *sli_1163 = buffer.data(sli + 1163);
    const auto *sli_1165 = buffer.data(sli + 1165);
    const auto *sli_1166 = buffer.data(sli + 1166);
    const auto *sli_1168 = buffer.data(sli + 1168);
    const auto *sli_1169 = buffer.data(sli + 1169);
    const auto *sli_1170 = buffer.data(sli + 1170);
    const auto *sli_1171 = buffer.data(sli + 1171);
    const auto *sli_1172 = buffer.data(sli + 1172);
    const auto *sli_1173 = buffer.data(sli + 1173);
    const auto *sli_1174 = buffer.data(sli + 1174);
    const auto *sli_1175 = buffer.data(sli + 1175);
    const auto *sli_1176 = buffer.data(sli + 1176);
    const auto *sli_1179 = buffer.data(sli + 1179);
    const auto *sli_1181 = buffer.data(sli + 1181);
    const auto *sli_1182 = buffer.data(sli + 1182);
    const auto *sli_1185 = buffer.data(sli + 1185);
    const auto *sli_1186 = buffer.data(sli + 1186);
    const auto *sli_1188 = buffer.data(sli + 1188);
    const auto *sli_1190 = buffer.data(sli + 1190);
    const auto *sli_1191 = buffer.data(sli + 1191);
    const auto *sli_1193 = buffer.data(sli + 1193);
    const auto *sli_1194 = buffer.data(sli + 1194);
    const auto *sli_1196 = buffer.data(sli + 1196);
    const auto *sli_1197 = buffer.data(sli + 1197);
    const auto *sli_1198 = buffer.data(sli + 1198);
    const auto *sli_1199 = buffer.data(sli + 1199);
    const auto *sli_1200 = buffer.data(sli + 1200);
    const auto *sli_1201 = buffer.data(sli + 1201);
    const auto *sli_1202 = buffer.data(sli + 1202);
    const auto *sli_1203 = buffer.data(sli + 1203);
    const auto *sli_1207 = buffer.data(sli + 1207);
    const auto *sli_1210 = buffer.data(sli + 1210);
    const auto *sli_1214 = buffer.data(sli + 1214);
    const auto *sli_1216 = buffer.data(sli + 1216);
    const auto *sli_1219 = buffer.data(sli + 1219);
    const auto *sli_1221 = buffer.data(sli + 1221);
    const auto *sli_1222 = buffer.data(sli + 1222);
    const auto *sli_1225 = buffer.data(sli + 1225);
    const auto *sli_1226 = buffer.data(sli + 1226);
    const auto *sli_1227 = buffer.data(sli + 1227);

    const auto *slk1_1260 = buffer.data(slk1 + 1260);
    const auto *slk1_1265 = buffer.data(slk1 + 1265);
    const auto *slk1_1269 = buffer.data(slk1 + 1269);
    const auto *slk1_1274 = buffer.data(slk1 + 1274);
    const auto *slk1_1280 = buffer.data(slk1 + 1280);
    const auto *slk1_1455 = buffer.data(slk1 + 1455);
    const auto *slk1_1457 = buffer.data(slk1 + 1457);
    const auto *slk1_1458 = buffer.data(slk1 + 1458);
    const auto *slk1_1460 = buffer.data(slk1 + 1460);
    const auto *slk1_1468 = buffer.data(slk1 + 1468);
    const auto *slk1_1470 = buffer.data(slk1 + 1470);
    const auto *slk1_1471 = buffer.data(slk1 + 1471);
    const auto *slk1_1472 = buffer.data(slk1 + 1472);
    const auto *slk1_1473 = buffer.data(slk1 + 1473);
    const auto *slk1_1475 = buffer.data(slk1 + 1475);
    const auto *slk1_1476 = buffer.data(slk1 + 1476);
    const auto *slk1_1479 = buffer.data(slk1 + 1479);
    const auto *slk1_1481 = buffer.data(slk1 + 1481);
    const auto *slk1_1482 = buffer.data(slk1 + 1482);
    const auto *slk1_1485 = buffer.data(slk1 + 1485);
    const auto *slk1_1486 = buffer.data(slk1 + 1486);
    const auto *slk1_1488 = buffer.data(slk1 + 1488);
    const auto *slk1_1490 = buffer.data(slk1 + 1490);
    const auto *slk1_1491 = buffer.data(slk1 + 1491);
    const auto *slk1_1493 = buffer.data(slk1 + 1493);
    const auto *slk1_1494 = buffer.data(slk1 + 1494);
    const auto *slk1_1496 = buffer.data(slk1 + 1496);
    const auto *slk1_1504 = buffer.data(slk1 + 1504);
    const auto *slk1_1506 = buffer.data(slk1 + 1506);
    const auto *slk1_1507 = buffer.data(slk1 + 1507);
    const auto *slk1_1508 = buffer.data(slk1 + 1508);
    const auto *slk1_1509 = buffer.data(slk1 + 1509);
    const auto *slk1_1511 = buffer.data(slk1 + 1511);
    const auto *slk1_1512 = buffer.data(slk1 + 1512);
    const auto *slk1_1515 = buffer.data(slk1 + 1515);
    const auto *slk1_1517 = buffer.data(slk1 + 1517);
    const auto *slk1_1518 = buffer.data(slk1 + 1518);
    const auto *slk1_1521 = buffer.data(slk1 + 1521);
    const auto *slk1_1522 = buffer.data(slk1 + 1522);
    const auto *slk1_1524 = buffer.data(slk1 + 1524);
    const auto *slk1_1526 = buffer.data(slk1 + 1526);
    const auto *slk1_1527 = buffer.data(slk1 + 1527);
    const auto *slk1_1529 = buffer.data(slk1 + 1529);
    const auto *slk1_1530 = buffer.data(slk1 + 1530);
    const auto *slk1_1532 = buffer.data(slk1 + 1532);
    const auto *slk1_1540 = buffer.data(slk1 + 1540);
    const auto *slk1_1542 = buffer.data(slk1 + 1542);
    const auto *slk1_1543 = buffer.data(slk1 + 1543);
    const auto *slk1_1544 = buffer.data(slk1 + 1544);
    const auto *slk1_1545 = buffer.data(slk1 + 1545);
    const auto *slk1_1547 = buffer.data(slk1 + 1547);
    const auto *slk1_1551 = buffer.data(slk1 + 1551);
    const auto *slk1_1554 = buffer.data(slk1 + 1554);
    const auto *slk1_1558 = buffer.data(slk1 + 1558);
    const auto *slk1_1560 = buffer.data(slk1 + 1560);
    const auto *slk1_1563 = buffer.data(slk1 + 1563);
    const auto *slk1_1565 = buffer.data(slk1 + 1565);
    const auto *slk1_1566 = buffer.data(slk1 + 1566);

    const auto *smi_1130 = buffer.data(smi + 1130);
    const auto *smi_1134 = buffer.data(smi + 1134);
    const auto *smi_1141 = buffer.data(smi + 1141);
    const auto *smi_1142 = buffer.data(smi + 1142);
    const auto *smi_1143 = buffer.data(smi + 1143);
    const auto *smi_1144 = buffer.data(smi + 1144);
    const auto *smi_1145 = buffer.data(smi + 1145);
    const auto *smi_1146 = buffer.data(smi + 1146);
    const auto *smi_1147 = buffer.data(smi + 1147);
    const auto *smi_1148 = buffer.data(smi + 1148);
    const auto *smi_1150 = buffer.data(smi + 1150);
    const auto *smi_1151 = buffer.data(smi + 1151);
    const auto *smi_1153 = buffer.data(smi + 1153);
    const auto *smi_1154 = buffer.data(smi + 1154);
    const auto *smi_1157 = buffer.data(smi + 1157);
    const auto *smi_1158 = buffer.data(smi + 1158);
    const auto *smi_1162 = buffer.data(smi + 1162);
    const auto *smi_1169 = buffer.data(smi + 1169);
    const auto *smi_1170 = buffer.data(smi + 1170);
    const auto *smi_1171 = buffer.data(smi + 1171);
    const auto *smi_1172 = buffer.data(smi + 1172);
    const auto *smi_1173 = buffer.data(smi + 1173);
    const auto *smi_1174 = buffer.data(smi + 1174);
    const auto *smi_1175 = buffer.data(smi + 1175);
    const auto *smi_1176 = buffer.data(smi + 1176);
    const auto *smi_1178 = buffer.data(smi + 1178);
    const auto *smi_1179 = buffer.data(smi + 1179);
    const auto *smi_1181 = buffer.data(smi + 1181);
    const auto *smi_1182 = buffer.data(smi + 1182);
    const auto *smi_1185 = buffer.data(smi + 1185);
    const auto *smi_1186 = buffer.data(smi + 1186);
    const auto *smi_1190 = buffer.data(smi + 1190);
    const auto *smi_1197 = buffer.data(smi + 1197);
    const auto *smi_1198 = buffer.data(smi + 1198);
    const auto *smi_1199 = buffer.data(smi + 1199);
    const auto *smi_1200 = buffer.data(smi + 1200);
    const auto *smi_1201 = buffer.data(smi + 1201);
    const auto *smi_1202 = buffer.data(smi + 1202);
    const auto *smi_1203 = buffer.data(smi + 1203);
    const auto *smi_1204 = buffer.data(smi + 1204);
    const auto *smi_1206 = buffer.data(smi + 1206);
    const auto *smi_1207 = buffer.data(smi + 1207);
    const auto *smi_1209 = buffer.data(smi + 1209);
    const auto *smi_1210 = buffer.data(smi + 1210);
    const auto *smi_1213 = buffer.data(smi + 1213);
    const auto *smi_1214 = buffer.data(smi + 1214);
    const auto *smi_1218 = buffer.data(smi + 1218);
    const auto *smi_1225 = buffer.data(smi + 1225);
    const auto *smi_1226 = buffer.data(smi + 1226);
    const auto *smi_1227 = buffer.data(smi + 1227);

#pragma omp simd aligned(t_1455, t_1456, t_1457, pb_x, pc_x, pc_z, slk0_1455, slk0_1457, \
                         sli_878, sli_1135, sli_1137, slk1_1455, slk1_1457, \
                         smi_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1455[k] = pb_x[k] * slk0_1455[k]
                    + f_14 * sli_1135[k]
                    - f_12 * pc_x[k] * slk1_1455[k];

        t_1456[k] = f_16 * sli_878[k]
                    + f_3 * pc_z[k] * smi_1130[k];

        t_1457[k] = pb_x[k] * slk0_1457[k]
                    + f_14 * sli_1137[k]
                    - f_12 * pc_x[k] * slk1_1457[k];
    }

#pragma omp simd aligned(t_1458, t_1459, t_1460, pb_x, pc_x, pc_y, slk0_1458, slk0_1460, \
                         sli_910, sli_1138, sli_1140, slk1_1458, slk1_1460, \
                         smi_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1458[k] = pb_x[k] * slk0_1458[k]
                    + f_14 * sli_1138[k]
                    - f_12 * pc_x[k] * slk1_1458[k];

        t_1459[k] = f_16 * sli_910[k]
                    + f_3 * pc_y[k] * smi_1134[k];

        t_1460[k] = pb_x[k] * slk0_1460[k]
                    + f_14 * sli_1140[k]
                    - f_12 * pc_x[k] * slk1_1460[k];
    }

#pragma omp simd aligned(t_1461, t_1462, t_1463, t_1464, t_1465, pc_x, sli_1141, sli_1142, \
                         sli_1143, sli_1144, sli_1145, smi_1141, smi_1142, smi_1143, smi_1144, \
                         smi_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1461[k] = f_13 * sli_1141[k]
                    + f_3 * pc_x[k] * smi_1141[k];

        t_1462[k] = f_13 * sli_1142[k]
                    + f_3 * pc_x[k] * smi_1142[k];

        t_1463[k] = f_13 * sli_1143[k]
                    + f_3 * pc_x[k] * smi_1143[k];

        t_1464[k] = f_13 * sli_1144[k]
                    + f_3 * pc_x[k] * smi_1144[k];

        t_1465[k] = f_13 * sli_1145[k]
                    + f_3 * pc_x[k] * smi_1145[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, t_1469, pb_x, pc_x, pc_z, slk0_1468, sli_889, \
                         sli_1146, sli_1147, slk1_1468, smi_1141, smi_1146, \
                         smi_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_13 * sli_1146[k]
                    + f_3 * pc_x[k] * smi_1146[k];

        t_1467[k] = f_13 * sli_1147[k]
                    + f_3 * pc_x[k] * smi_1147[k];

        t_1468[k] = pb_x[k] * slk0_1468[k]
                    - f_12 * pc_x[k] * slk1_1468[k];

        t_1469[k] = f_16 * sli_889[k]
                    + f_3 * pc_z[k] * smi_1141[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, pb_x, pc_x, slk0_1470, slk0_1471, \
                         slk0_1472, slk0_1473, slk1_1470, slk1_1471, slk1_1472, \
                         slk1_1473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = pb_x[k] * slk0_1470[k]
                    - f_12 * pc_x[k] * slk1_1470[k];

        t_1471[k] = pb_x[k] * slk0_1471[k]
                    - f_12 * pc_x[k] * slk1_1471[k];

        t_1472[k] = pb_x[k] * slk0_1472[k]
                    - f_12 * pc_x[k] * slk1_1472[k];

        t_1473[k] = pb_x[k] * slk0_1473[k]
                    - f_12 * pc_x[k] * slk1_1473[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, pb_x, pc_x, pc_y, slk0_1475, \
                         slk0_1476, sli_923, sli_924, sli_1148, slk1_1475, slk1_1476, \
                         smi_1147, smi_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_16 * sli_923[k]
                    + f_3 * pc_y[k] * smi_1147[k];

        t_1475[k] = pb_x[k] * slk0_1475[k]
                    - f_12 * pc_x[k] * slk1_1475[k];

        t_1476[k] = pb_x[k] * slk0_1476[k]
                    + f_19 * sli_1148[k]
                    - f_12 * pc_x[k] * slk1_1476[k];

        t_1477[k] = f_15 * sli_924[k]
                    + f_3 * pc_y[k] * smi_1148[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, pb_x, pc_x, pc_y, pc_z, slk0_1479, sli_896, \
                         sli_926, sli_1151, slk1_1479, smi_1148, \
                         smi_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_17 * sli_896[k]
                    + f_3 * pc_z[k] * smi_1148[k];

        t_1479[k] = pb_x[k] * slk0_1479[k]
                    + f_17 * sli_1151[k]
                    - f_12 * pc_x[k] * slk1_1479[k];

        t_1480[k] = f_15 * sli_926[k]
                    + f_3 * pc_y[k] * smi_1150[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, pb_x, pc_x, pc_z, slk0_1481, slk0_1482, \
                         sli_899, sli_1153, sli_1154, slk1_1481, slk1_1482, \
                         smi_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = pb_x[k] * slk0_1481[k]
                    + f_17 * sli_1153[k]
                    - f_12 * pc_x[k] * slk1_1481[k];

        t_1482[k] = pb_x[k] * slk0_1482[k]
                    + f_16 * sli_1154[k]
                    - f_12 * pc_x[k] * slk1_1482[k];

        t_1483[k] = f_17 * sli_899[k]
                    + f_3 * pc_z[k] * smi_1151[k];
    }

#pragma omp simd aligned(t_1484, t_1485, t_1486, pb_x, pc_x, pc_y, slk0_1485, slk0_1486, \
                         sli_929, sli_1157, sli_1158, slk1_1485, slk1_1486, \
                         smi_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1484[k] = f_15 * sli_929[k]
                    + f_3 * pc_y[k] * smi_1153[k];

        t_1485[k] = pb_x[k] * slk0_1485[k]
                    + f_16 * sli_1157[k]
                    - f_12 * pc_x[k] * slk1_1485[k];

        t_1486[k] = pb_x[k] * slk0_1486[k]
                    + f_15 * sli_1158[k]
                    - f_12 * pc_x[k] * slk1_1486[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pb_x, pc_x, pc_y, pc_z, slk0_1488, sli_902, \
                         sli_933, sli_1160, slk1_1488, smi_1154, \
                         smi_1157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_17 * sli_902[k]
                    + f_3 * pc_z[k] * smi_1154[k];

        t_1488[k] = pb_x[k] * slk0_1488[k]
                    + f_15 * sli_1160[k]
                    - f_12 * pc_x[k] * slk1_1488[k];

        t_1489[k] = f_15 * sli_933[k]
                    + f_3 * pc_y[k] * smi_1157[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pb_x, pc_x, pc_z, slk0_1490, slk0_1491, \
                         sli_906, sli_1162, sli_1163, slk1_1490, slk1_1491, \
                         smi_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = pb_x[k] * slk0_1490[k]
                    + f_15 * sli_1162[k]
                    - f_12 * pc_x[k] * slk1_1490[k];

        t_1491[k] = pb_x[k] * slk0_1491[k]
                    + f_14 * sli_1163[k]
                    - f_12 * pc_x[k] * slk1_1491[k];

        t_1492[k] = f_17 * sli_906[k]
                    + f_3 * pc_z[k] * smi_1158[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pb_x, pc_x, pc_y, slk0_1493, slk0_1494, \
                         sli_938, sli_1165, sli_1166, slk1_1493, slk1_1494, \
                         smi_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = pb_x[k] * slk0_1493[k]
                    + f_14 * sli_1165[k]
                    - f_12 * pc_x[k] * slk1_1493[k];

        t_1494[k] = pb_x[k] * slk0_1494[k]
                    + f_14 * sli_1166[k]
                    - f_12 * pc_x[k] * slk1_1494[k];

        t_1495[k] = f_15 * sli_938[k]
                    + f_3 * pc_y[k] * smi_1162[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, t_1499, pb_x, pc_x, slk0_1496, sli_1168, \
                         sli_1169, sli_1170, sli_1171, slk1_1496, smi_1169, smi_1170, \
                         smi_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = pb_x[k] * slk0_1496[k]
                    + f_14 * sli_1168[k]
                    - f_12 * pc_x[k] * slk1_1496[k];

        t_1497[k] = f_13 * sli_1169[k]
                    + f_3 * pc_x[k] * smi_1169[k];

        t_1498[k] = f_13 * sli_1170[k]
                    + f_3 * pc_x[k] * smi_1170[k];

        t_1499[k] = f_13 * sli_1171[k]
                    + f_3 * pc_x[k] * smi_1171[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, t_1503, pc_x, sli_1172, sli_1173, sli_1174, \
                         sli_1175, smi_1172, smi_1173, smi_1174, \
                         smi_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_13 * sli_1172[k]
                    + f_3 * pc_x[k] * smi_1172[k];

        t_1501[k] = f_13 * sli_1173[k]
                    + f_3 * pc_x[k] * smi_1173[k];

        t_1502[k] = f_13 * sli_1174[k]
                    + f_3 * pc_x[k] * smi_1174[k];

        t_1503[k] = f_13 * sli_1175[k]
                    + f_3 * pc_x[k] * smi_1175[k];
    }

#pragma omp simd aligned(t_1504, t_1505, t_1506, t_1507, pb_x, pc_x, pc_z, slk0_1504, \
                         slk0_1506, slk0_1507, sli_917, slk1_1504, slk1_1506, slk1_1507, \
                         smi_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1504[k] = pb_x[k] * slk0_1504[k]
                    - f_12 * pc_x[k] * slk1_1504[k];

        t_1505[k] = f_17 * sli_917[k]
                    + f_3 * pc_z[k] * smi_1169[k];

        t_1506[k] = pb_x[k] * slk0_1506[k]
                    - f_12 * pc_x[k] * slk1_1506[k];

        t_1507[k] = pb_x[k] * slk0_1507[k]
                    - f_12 * pc_x[k] * slk1_1507[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, t_1511, pb_x, pc_x, pc_y, slk0_1508, \
                         slk0_1509, slk0_1511, sli_951, slk1_1508, slk1_1509, slk1_1511, \
                         smi_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = pb_x[k] * slk0_1508[k]
                    - f_12 * pc_x[k] * slk1_1508[k];

        t_1509[k] = pb_x[k] * slk0_1509[k]
                    - f_12 * pc_x[k] * slk1_1509[k];

        t_1510[k] = f_15 * sli_951[k]
                    + f_3 * pc_y[k] * smi_1175[k];

        t_1511[k] = pb_x[k] * slk0_1511[k]
                    - f_12 * pc_x[k] * slk1_1511[k];
    }

#pragma omp simd aligned(t_1512, t_1513, t_1514, pb_x, pc_x, pc_y, pc_z, slk0_1512, sli_924, \
                         sli_952, sli_1176, slk1_1512, smi_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1512[k] = pb_x[k] * slk0_1512[k]
                    + f_19 * sli_1176[k]
                    - f_12 * pc_x[k] * slk1_1512[k];

        t_1513[k] = f_14 * sli_952[k]
                    + f_3 * pc_y[k] * smi_1176[k];

        t_1514[k] = f_20 * sli_924[k]
                    + f_3 * pc_z[k] * smi_1176[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, pb_x, pc_x, pc_y, slk0_1515, slk0_1517, \
                         sli_954, sli_1179, sli_1181, slk1_1515, slk1_1517, \
                         smi_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = pb_x[k] * slk0_1515[k]
                    + f_17 * sli_1179[k]
                    - f_12 * pc_x[k] * slk1_1515[k];

        t_1516[k] = f_14 * sli_954[k]
                    + f_3 * pc_y[k] * smi_1178[k];

        t_1517[k] = pb_x[k] * slk0_1517[k]
                    + f_17 * sli_1181[k]
                    - f_12 * pc_x[k] * slk1_1517[k];
    }

#pragma omp simd aligned(t_1518, t_1519, t_1520, pb_x, pc_x, pc_y, pc_z, slk0_1518, sli_927, \
                         sli_957, sli_1182, slk1_1518, smi_1179, \
                         smi_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = pb_x[k] * slk0_1518[k]
                    + f_16 * sli_1182[k]
                    - f_12 * pc_x[k] * slk1_1518[k];

        t_1519[k] = f_20 * sli_927[k]
                    + f_3 * pc_z[k] * smi_1179[k];

        t_1520[k] = f_14 * sli_957[k]
                    + f_3 * pc_y[k] * smi_1181[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, pb_x, pc_x, pc_z, slk0_1521, slk0_1522, \
                         sli_930, sli_1185, sli_1186, slk1_1521, slk1_1522, \
                         smi_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = pb_x[k] * slk0_1521[k]
                    + f_16 * sli_1185[k]
                    - f_12 * pc_x[k] * slk1_1521[k];

        t_1522[k] = pb_x[k] * slk0_1522[k]
                    + f_15 * sli_1186[k]
                    - f_12 * pc_x[k] * slk1_1522[k];

        t_1523[k] = f_20 * sli_930[k]
                    + f_3 * pc_z[k] * smi_1182[k];
    }

#pragma omp simd aligned(t_1524, t_1525, t_1526, pb_x, pc_x, pc_y, slk0_1524, slk0_1526, \
                         sli_961, sli_1188, sli_1190, slk1_1524, slk1_1526, \
                         smi_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1524[k] = pb_x[k] * slk0_1524[k]
                    + f_15 * sli_1188[k]
                    - f_12 * pc_x[k] * slk1_1524[k];

        t_1525[k] = f_14 * sli_961[k]
                    + f_3 * pc_y[k] * smi_1185[k];

        t_1526[k] = pb_x[k] * slk0_1526[k]
                    + f_15 * sli_1190[k]
                    - f_12 * pc_x[k] * slk1_1526[k];
    }

#pragma omp simd aligned(t_1527, t_1528, t_1529, pb_x, pc_x, pc_z, slk0_1527, slk0_1529, \
                         sli_934, sli_1191, sli_1193, slk1_1527, slk1_1529, \
                         smi_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1527[k] = pb_x[k] * slk0_1527[k]
                    + f_14 * sli_1191[k]
                    - f_12 * pc_x[k] * slk1_1527[k];

        t_1528[k] = f_20 * sli_934[k]
                    + f_3 * pc_z[k] * smi_1186[k];

        t_1529[k] = pb_x[k] * slk0_1529[k]
                    + f_14 * sli_1193[k]
                    - f_12 * pc_x[k] * slk1_1529[k];
    }

#pragma omp simd aligned(t_1530, t_1531, t_1532, pb_x, pc_x, pc_y, slk0_1530, slk0_1532, \
                         sli_966, sli_1194, sli_1196, slk1_1530, slk1_1532, \
                         smi_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1530[k] = pb_x[k] * slk0_1530[k]
                    + f_14 * sli_1194[k]
                    - f_12 * pc_x[k] * slk1_1530[k];

        t_1531[k] = f_14 * sli_966[k]
                    + f_3 * pc_y[k] * smi_1190[k];

        t_1532[k] = pb_x[k] * slk0_1532[k]
                    + f_14 * sli_1196[k]
                    - f_12 * pc_x[k] * slk1_1532[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, t_1536, t_1537, pc_x, sli_1197, sli_1198, \
                         sli_1199, sli_1200, sli_1201, smi_1197, smi_1198, smi_1199, smi_1200, \
                         smi_1201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = f_13 * sli_1197[k]
                    + f_3 * pc_x[k] * smi_1197[k];

        t_1534[k] = f_13 * sli_1198[k]
                    + f_3 * pc_x[k] * smi_1198[k];

        t_1535[k] = f_13 * sli_1199[k]
                    + f_3 * pc_x[k] * smi_1199[k];

        t_1536[k] = f_13 * sli_1200[k]
                    + f_3 * pc_x[k] * smi_1200[k];

        t_1537[k] = f_13 * sli_1201[k]
                    + f_3 * pc_x[k] * smi_1201[k];
    }

#pragma omp simd aligned(t_1538, t_1539, t_1540, t_1541, pb_x, pc_x, pc_z, slk0_1540, sli_945, \
                         sli_1202, sli_1203, slk1_1540, smi_1197, smi_1202, \
                         smi_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1538[k] = f_13 * sli_1202[k]
                    + f_3 * pc_x[k] * smi_1202[k];

        t_1539[k] = f_13 * sli_1203[k]
                    + f_3 * pc_x[k] * smi_1203[k];

        t_1540[k] = pb_x[k] * slk0_1540[k]
                    - f_12 * pc_x[k] * slk1_1540[k];

        t_1541[k] = f_20 * sli_945[k]
                    + f_3 * pc_z[k] * smi_1197[k];
    }

#pragma omp simd aligned(t_1542, t_1543, t_1544, t_1545, pb_x, pc_x, slk0_1542, slk0_1543, \
                         slk0_1544, slk0_1545, slk1_1542, slk1_1543, slk1_1544, \
                         slk1_1545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1542[k] = pb_x[k] * slk0_1542[k]
                    - f_12 * pc_x[k] * slk1_1542[k];

        t_1543[k] = pb_x[k] * slk0_1543[k]
                    - f_12 * pc_x[k] * slk1_1543[k];

        t_1544[k] = pb_x[k] * slk0_1544[k]
                    - f_12 * pc_x[k] * slk1_1544[k];

        t_1545[k] = pb_x[k] * slk0_1545[k]
                    - f_12 * pc_x[k] * slk1_1545[k];
    }

#pragma omp simd aligned(t_1546, t_1547, t_1548, t_1549, pb_x, pb_y, pc_x, pc_y, slk0_1260, \
                         slk0_1547, sli_979, sli_980, slk1_1260, slk1_1547, smi_1203, \
                         smi_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1546[k] = f_14 * sli_979[k]
                    + f_3 * pc_y[k] * smi_1203[k];

        t_1547[k] = pb_x[k] * slk0_1547[k]
                    - f_12 * pc_x[k] * slk1_1547[k];

        t_1548[k] = pb_y[k] * slk0_1260[k]
                    - f_12 * pc_y[k] * slk1_1260[k];

        t_1549[k] = f_13 * sli_980[k]
                    + f_3 * pc_y[k] * smi_1204[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, pb_x, pc_x, pc_y, pc_z, slk0_1551, sli_952, \
                         sli_982, sli_1207, slk1_1551, smi_1204, \
                         smi_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = f_19 * sli_952[k]
                    + f_3 * pc_z[k] * smi_1204[k];

        t_1551[k] = pb_x[k] * slk0_1551[k]
                    + f_17 * sli_1207[k]
                    - f_12 * pc_x[k] * slk1_1551[k];

        t_1552[k] = f_13 * sli_982[k]
                    + f_3 * pc_y[k] * smi_1206[k];
    }

#pragma omp simd aligned(t_1553, t_1554, t_1555, pb_x, pb_y, pc_x, pc_y, pc_z, slk0_1265, \
                         slk0_1554, sli_955, sli_1210, slk1_1265, slk1_1554, \
                         smi_1207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1553[k] = pb_y[k] * slk0_1265[k]
                    - f_12 * pc_y[k] * slk1_1265[k];

        t_1554[k] = pb_x[k] * slk0_1554[k]
                    + f_16 * sli_1210[k]
                    - f_12 * pc_x[k] * slk1_1554[k];

        t_1555[k] = f_19 * sli_955[k]
                    + f_3 * pc_z[k] * smi_1207[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, pb_x, pb_y, pc_x, pc_y, slk0_1269, slk0_1558, \
                         sli_985, sli_1214, slk1_1269, slk1_1558, \
                         smi_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_13 * sli_985[k]
                    + f_3 * pc_y[k] * smi_1209[k];

        t_1557[k] = pb_y[k] * slk0_1269[k]
                    - f_12 * pc_y[k] * slk1_1269[k];

        t_1558[k] = pb_x[k] * slk0_1558[k]
                    + f_15 * sli_1214[k]
                    - f_12 * pc_x[k] * slk1_1558[k];
    }

#pragma omp simd aligned(t_1559, t_1560, t_1561, pb_x, pc_x, pc_y, pc_z, slk0_1560, sli_958, \
                         sli_989, sli_1216, slk1_1560, smi_1210, \
                         smi_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1559[k] = f_19 * sli_958[k]
                    + f_3 * pc_z[k] * smi_1210[k];

        t_1560[k] = pb_x[k] * slk0_1560[k]
                    + f_15 * sli_1216[k]
                    - f_12 * pc_x[k] * slk1_1560[k];

        t_1561[k] = f_13 * sli_989[k]
                    + f_3 * pc_y[k] * smi_1213[k];
    }

#pragma omp simd aligned(t_1562, t_1563, t_1564, pb_x, pb_y, pc_x, pc_y, pc_z, slk0_1274, \
                         slk0_1563, sli_962, sli_1219, slk1_1274, slk1_1563, \
                         smi_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1562[k] = pb_y[k] * slk0_1274[k]
                    - f_12 * pc_y[k] * slk1_1274[k];

        t_1563[k] = pb_x[k] * slk0_1563[k]
                    + f_14 * sli_1219[k]
                    - f_12 * pc_x[k] * slk1_1563[k];

        t_1564[k] = f_19 * sli_962[k]
                    + f_3 * pc_z[k] * smi_1214[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, pb_x, pc_x, pc_y, slk0_1565, slk0_1566, \
                         sli_994, sli_1221, sli_1222, slk1_1565, slk1_1566, \
                         smi_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = pb_x[k] * slk0_1565[k]
                    + f_14 * sli_1221[k]
                    - f_12 * pc_x[k] * slk1_1565[k];

        t_1566[k] = pb_x[k] * slk0_1566[k]
                    + f_14 * sli_1222[k]
                    - f_12 * pc_x[k] * slk1_1566[k];

        t_1567[k] = f_13 * sli_994[k]
                    + f_3 * pc_y[k] * smi_1218[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, t_1571, pb_y, pc_x, pc_y, slk0_1280, \
                         sli_1225, sli_1226, sli_1227, slk1_1280, smi_1225, smi_1226, \
                         smi_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = pb_y[k] * slk0_1280[k]
                    - f_12 * pc_y[k] * slk1_1280[k];

        t_1569[k] = f_13 * sli_1225[k]
                    + f_3 * pc_x[k] * smi_1225[k];

        t_1570[k] = f_13 * sli_1226[k]
                    + f_3 * pc_x[k] * smi_1226[k];

        t_1571[k] = f_13 * sli_1227[k]
                    + f_3 * pc_x[k] * smi_1227[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t slk0,
                                                           const size_t sli, const size_t slk1,
                                                           const size_t smh0, const size_t smh1,
                                                           const size_t smi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_1296 = buffer.data(slk0 + 1296);
    const auto *slk0_1299 = buffer.data(slk0 + 1299);
    const auto *slk0_1302 = buffer.data(slk0 + 1302);
    const auto *slk0_1306 = buffer.data(slk0 + 1306);
    const auto *slk0_1311 = buffer.data(slk0 + 1311);
    const auto *slk0_1324 = buffer.data(slk0 + 1324);
    const auto *slk0_1326 = buffer.data(slk0 + 1326);
    const auto *slk0_1327 = buffer.data(slk0 + 1327);
    const auto *slk0_1328 = buffer.data(slk0 + 1328);
    const auto *slk0_1329 = buffer.data(slk0 + 1329);
    const auto *slk0_1576 = buffer.data(slk0 + 1576);
    const auto *slk0_1578 = buffer.data(slk0 + 1578);
    const auto *slk0_1579 = buffer.data(slk0 + 1579);
    const auto *slk0_1580 = buffer.data(slk0 + 1580);
    const auto *slk0_1581 = buffer.data(slk0 + 1581);
    const auto *slk0_1583 = buffer.data(slk0 + 1583);
    const auto *slk0_1584 = buffer.data(slk0 + 1584);
    const auto *slk0_1587 = buffer.data(slk0 + 1587);
    const auto *slk0_1589 = buffer.data(slk0 + 1589);
    const auto *slk0_1590 = buffer.data(slk0 + 1590);
    const auto *slk0_1593 = buffer.data(slk0 + 1593);
    const auto *slk0_1594 = buffer.data(slk0 + 1594);
    const auto *slk0_1596 = buffer.data(slk0 + 1596);
    const auto *slk0_1598 = buffer.data(slk0 + 1598);
    const auto *slk0_1599 = buffer.data(slk0 + 1599);
    const auto *slk0_1601 = buffer.data(slk0 + 1601);
    const auto *slk0_1602 = buffer.data(slk0 + 1602);
    const auto *slk0_1604 = buffer.data(slk0 + 1604);
    const auto *slk0_1612 = buffer.data(slk0 + 1612);
    const auto *slk0_1614 = buffer.data(slk0 + 1614);
    const auto *slk0_1615 = buffer.data(slk0 + 1615);
    const auto *slk0_1616 = buffer.data(slk0 + 1616);
    const auto *slk0_1617 = buffer.data(slk0 + 1617);
    const auto *slk0_1619 = buffer.data(slk0 + 1619);

    const auto *sli_973 = buffer.data(sli + 973);
    const auto *sli_980 = buffer.data(sli + 980);
    const auto *sli_983 = buffer.data(sli + 983);
    const auto *sli_986 = buffer.data(sli + 986);
    const auto *sli_990 = buffer.data(sli + 990);
    const auto *sli_1001 = buffer.data(sli + 1001);
    const auto *sli_1007 = buffer.data(sli + 1007);
    const auto *sli_1008 = buffer.data(sli + 1008);
    const auto *sli_1010 = buffer.data(sli + 1010);
    const auto *sli_1011 = buffer.data(sli + 1011);
    const auto *sli_1013 = buffer.data(sli + 1013);
    const auto *sli_1014 = buffer.data(sli + 1014);
    const auto *sli_1017 = buffer.data(sli + 1017);
    const auto *sli_1018 = buffer.data(sli + 1018);
    const auto *sli_1022 = buffer.data(sli + 1022);
    const auto *sli_1029 = buffer.data(sli + 1029);
    const auto *sli_1030 = buffer.data(sli + 1030);
    const auto *sli_1031 = buffer.data(sli + 1031);
    const auto *sli_1032 = buffer.data(sli + 1032);
    const auto *sli_1033 = buffer.data(sli + 1033);
    const auto *sli_1034 = buffer.data(sli + 1034);
    const auto *sli_1035 = buffer.data(sli + 1035);
    const auto *sli_1036 = buffer.data(sli + 1036);
    const auto *sli_1038 = buffer.data(sli + 1038);
    const auto *sli_1041 = buffer.data(sli + 1041);
    const auto *sli_1045 = buffer.data(sli + 1045);
    const auto *sli_1050 = buffer.data(sli + 1050);
    const auto *sli_1063 = buffer.data(sli + 1063);
    const auto *sli_1064 = buffer.data(sli + 1064);
    const auto *sli_1066 = buffer.data(sli + 1066);
    const auto *sli_1228 = buffer.data(sli + 1228);
    const auto *sli_1229 = buffer.data(sli + 1229);
    const auto *sli_1230 = buffer.data(sli + 1230);
    const auto *sli_1231 = buffer.data(sli + 1231);
    const auto *sli_1232 = buffer.data(sli + 1232);
    const auto *sli_1235 = buffer.data(sli + 1235);
    const auto *sli_1237 = buffer.data(sli + 1237);
    const auto *sli_1238 = buffer.data(sli + 1238);
    const auto *sli_1241 = buffer.data(sli + 1241);
    const auto *sli_1242 = buffer.data(sli + 1242);
    const auto *sli_1244 = buffer.data(sli + 1244);
    const auto *sli_1246 = buffer.data(sli + 1246);
    const auto *sli_1247 = buffer.data(sli + 1247);
    const auto *sli_1249 = buffer.data(sli + 1249);
    const auto *sli_1250 = buffer.data(sli + 1250);
    const auto *sli_1252 = buffer.data(sli + 1252);
    const auto *sli_1253 = buffer.data(sli + 1253);
    const auto *sli_1254 = buffer.data(sli + 1254);
    const auto *sli_1255 = buffer.data(sli + 1255);
    const auto *sli_1256 = buffer.data(sli + 1256);
    const auto *sli_1257 = buffer.data(sli + 1257);
    const auto *sli_1258 = buffer.data(sli + 1258);
    const auto *sli_1259 = buffer.data(sli + 1259);

    const auto *slk1_1296 = buffer.data(slk1 + 1296);
    const auto *slk1_1299 = buffer.data(slk1 + 1299);
    const auto *slk1_1302 = buffer.data(slk1 + 1302);
    const auto *slk1_1306 = buffer.data(slk1 + 1306);
    const auto *slk1_1311 = buffer.data(slk1 + 1311);
    const auto *slk1_1324 = buffer.data(slk1 + 1324);
    const auto *slk1_1326 = buffer.data(slk1 + 1326);
    const auto *slk1_1327 = buffer.data(slk1 + 1327);
    const auto *slk1_1328 = buffer.data(slk1 + 1328);
    const auto *slk1_1329 = buffer.data(slk1 + 1329);
    const auto *slk1_1576 = buffer.data(slk1 + 1576);
    const auto *slk1_1578 = buffer.data(slk1 + 1578);
    const auto *slk1_1579 = buffer.data(slk1 + 1579);
    const auto *slk1_1580 = buffer.data(slk1 + 1580);
    const auto *slk1_1581 = buffer.data(slk1 + 1581);
    const auto *slk1_1583 = buffer.data(slk1 + 1583);
    const auto *slk1_1584 = buffer.data(slk1 + 1584);
    const auto *slk1_1587 = buffer.data(slk1 + 1587);
    const auto *slk1_1589 = buffer.data(slk1 + 1589);
    const auto *slk1_1590 = buffer.data(slk1 + 1590);
    const auto *slk1_1593 = buffer.data(slk1 + 1593);
    const auto *slk1_1594 = buffer.data(slk1 + 1594);
    const auto *slk1_1596 = buffer.data(slk1 + 1596);
    const auto *slk1_1598 = buffer.data(slk1 + 1598);
    const auto *slk1_1599 = buffer.data(slk1 + 1599);
    const auto *slk1_1601 = buffer.data(slk1 + 1601);
    const auto *slk1_1602 = buffer.data(slk1 + 1602);
    const auto *slk1_1604 = buffer.data(slk1 + 1604);
    const auto *slk1_1612 = buffer.data(slk1 + 1612);
    const auto *slk1_1614 = buffer.data(slk1 + 1614);
    const auto *slk1_1615 = buffer.data(slk1 + 1615);
    const auto *slk1_1616 = buffer.data(slk1 + 1616);
    const auto *slk1_1617 = buffer.data(slk1 + 1617);
    const auto *slk1_1619 = buffer.data(slk1 + 1619);

    const auto *smh0_945 = buffer.data(smh0 + 945);
    const auto *smh0_948 = buffer.data(smh0 + 948);
    const auto *smh0_950 = buffer.data(smh0 + 950);
    const auto *smh0_951 = buffer.data(smh0 + 951);
    const auto *smh0_954 = buffer.data(smh0 + 954);
    const auto *smh0_955 = buffer.data(smh0 + 955);
    const auto *smh0_957 = buffer.data(smh0 + 957);
    const auto *smh0_959 = buffer.data(smh0 + 959);
    const auto *smh0_960 = buffer.data(smh0 + 960);
    const auto *smh0_962 = buffer.data(smh0 + 962);
    const auto *smh0_963 = buffer.data(smh0 + 963);
    const auto *smh0_964 = buffer.data(smh0 + 964);
    const auto *smh0_965 = buffer.data(smh0 + 965);
    const auto *smh0_971 = buffer.data(smh0 + 971);
    const auto *smh0_975 = buffer.data(smh0 + 975);
    const auto *smh0_978 = buffer.data(smh0 + 978);
    const auto *smh0_980 = buffer.data(smh0 + 980);
    const auto *smh0_983 = buffer.data(smh0 + 983);
    const auto *smh0_984 = buffer.data(smh0 + 984);
    const auto *smh0_986 = buffer.data(smh0 + 986);
    const auto *smh0_987 = buffer.data(smh0 + 987);
    const auto *smh0_990 = buffer.data(smh0 + 990);
    const auto *smh0_992 = buffer.data(smh0 + 992);

    const auto *smh1_945 = buffer.data(smh1 + 945);
    const auto *smh1_948 = buffer.data(smh1 + 948);
    const auto *smh1_950 = buffer.data(smh1 + 950);
    const auto *smh1_951 = buffer.data(smh1 + 951);
    const auto *smh1_954 = buffer.data(smh1 + 954);
    const auto *smh1_955 = buffer.data(smh1 + 955);
    const auto *smh1_957 = buffer.data(smh1 + 957);
    const auto *smh1_959 = buffer.data(smh1 + 959);
    const auto *smh1_960 = buffer.data(smh1 + 960);
    const auto *smh1_962 = buffer.data(smh1 + 962);
    const auto *smh1_963 = buffer.data(smh1 + 963);
    const auto *smh1_964 = buffer.data(smh1 + 964);
    const auto *smh1_965 = buffer.data(smh1 + 965);
    const auto *smh1_971 = buffer.data(smh1 + 971);
    const auto *smh1_975 = buffer.data(smh1 + 975);
    const auto *smh1_978 = buffer.data(smh1 + 978);
    const auto *smh1_980 = buffer.data(smh1 + 980);
    const auto *smh1_983 = buffer.data(smh1 + 983);
    const auto *smh1_984 = buffer.data(smh1 + 984);
    const auto *smh1_986 = buffer.data(smh1 + 986);
    const auto *smh1_987 = buffer.data(smh1 + 987);
    const auto *smh1_990 = buffer.data(smh1 + 990);
    const auto *smh1_992 = buffer.data(smh1 + 992);

    const auto *smi_1225 = buffer.data(smi + 1225);
    const auto *smi_1228 = buffer.data(smi + 1228);
    const auto *smi_1229 = buffer.data(smi + 1229);
    const auto *smi_1230 = buffer.data(smi + 1230);
    const auto *smi_1231 = buffer.data(smi + 1231);
    const auto *smi_1232 = buffer.data(smi + 1232);
    const auto *smi_1234 = buffer.data(smi + 1234);
    const auto *smi_1235 = buffer.data(smi + 1235);
    const auto *smi_1237 = buffer.data(smi + 1237);
    const auto *smi_1238 = buffer.data(smi + 1238);
    const auto *smi_1241 = buffer.data(smi + 1241);
    const auto *smi_1242 = buffer.data(smi + 1242);
    const auto *smi_1246 = buffer.data(smi + 1246);
    const auto *smi_1253 = buffer.data(smi + 1253);
    const auto *smi_1254 = buffer.data(smi + 1254);
    const auto *smi_1255 = buffer.data(smi + 1255);
    const auto *smi_1256 = buffer.data(smi + 1256);
    const auto *smi_1257 = buffer.data(smi + 1257);
    const auto *smi_1258 = buffer.data(smi + 1258);
    const auto *smi_1259 = buffer.data(smi + 1259);
    const auto *smi_1260 = buffer.data(smi + 1260);
    const auto *smi_1262 = buffer.data(smi + 1262);
    const auto *smi_1263 = buffer.data(smi + 1263);
    const auto *smi_1265 = buffer.data(smi + 1265);
    const auto *smi_1266 = buffer.data(smi + 1266);
    const auto *smi_1269 = buffer.data(smi + 1269);
    const auto *smi_1270 = buffer.data(smi + 1270);
    const auto *smi_1272 = buffer.data(smi + 1272);
    const auto *smi_1274 = buffer.data(smi + 1274);
    const auto *smi_1275 = buffer.data(smi + 1275);
    const auto *smi_1277 = buffer.data(smi + 1277);
    const auto *smi_1278 = buffer.data(smi + 1278);
    const auto *smi_1280 = buffer.data(smi + 1280);
    const auto *smi_1281 = buffer.data(smi + 1281);
    const auto *smi_1282 = buffer.data(smi + 1282);
    const auto *smi_1283 = buffer.data(smi + 1283);
    const auto *smi_1284 = buffer.data(smi + 1284);
    const auto *smi_1285 = buffer.data(smi + 1285);
    const auto *smi_1286 = buffer.data(smi + 1286);
    const auto *smi_1287 = buffer.data(smi + 1287);
    const auto *smi_1288 = buffer.data(smi + 1288);
    const auto *smi_1290 = buffer.data(smi + 1290);
    const auto *smi_1291 = buffer.data(smi + 1291);
    const auto *smi_1293 = buffer.data(smi + 1293);
    const auto *smi_1294 = buffer.data(smi + 1294);
    const auto *smi_1297 = buffer.data(smi + 1297);
    const auto *smi_1298 = buffer.data(smi + 1298);
    const auto *smi_1300 = buffer.data(smi + 1300);
    const auto *smi_1302 = buffer.data(smi + 1302);
    const auto *smi_1305 = buffer.data(smi + 1305);
    const auto *smi_1306 = buffer.data(smi + 1306);
    const auto *smi_1308 = buffer.data(smi + 1308);
    const auto *smi_1309 = buffer.data(smi + 1309);
    const auto *smi_1310 = buffer.data(smi + 1310);
    const auto *smi_1311 = buffer.data(smi + 1311);
    const auto *smi_1312 = buffer.data(smi + 1312);
    const auto *smi_1313 = buffer.data(smi + 1313);
    const auto *smi_1314 = buffer.data(smi + 1314);
    const auto *smi_1315 = buffer.data(smi + 1315);
    const auto *smi_1316 = buffer.data(smi + 1316);
    const auto *smi_1318 = buffer.data(smi + 1318);
    const auto *smi_1319 = buffer.data(smi + 1319);
    const auto *smi_1321 = buffer.data(smi + 1321);

#pragma omp simd aligned(t_1572, t_1573, t_1574, t_1575, pc_x, sli_1228, sli_1229, sli_1230, \
                         sli_1231, smi_1228, smi_1229, smi_1230, \
                         smi_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1572[k] = f_13 * sli_1228[k]
                    + f_3 * pc_x[k] * smi_1228[k];

        t_1573[k] = f_13 * sli_1229[k]
                    + f_3 * pc_x[k] * smi_1229[k];

        t_1574[k] = f_13 * sli_1230[k]
                    + f_3 * pc_x[k] * smi_1230[k];

        t_1575[k] = f_13 * sli_1231[k]
                    + f_3 * pc_x[k] * smi_1231[k];
    }

#pragma omp simd aligned(t_1576, t_1577, t_1578, t_1579, pb_x, pc_x, pc_z, slk0_1576, \
                         slk0_1578, slk0_1579, sli_973, slk1_1576, slk1_1578, slk1_1579, \
                         smi_1225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1576[k] = pb_x[k] * slk0_1576[k]
                    - f_12 * pc_x[k] * slk1_1576[k];

        t_1577[k] = f_19 * sli_973[k]
                    + f_3 * pc_z[k] * smi_1225[k];

        t_1578[k] = pb_x[k] * slk0_1578[k]
                    - f_12 * pc_x[k] * slk1_1578[k];

        t_1579[k] = pb_x[k] * slk0_1579[k]
                    - f_12 * pc_x[k] * slk1_1579[k];
    }

#pragma omp simd aligned(t_1580, t_1581, t_1582, t_1583, pb_x, pc_x, pc_y, slk0_1580, \
                         slk0_1581, slk0_1583, sli_1007, slk1_1580, slk1_1581, slk1_1583, \
                         smi_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1580[k] = pb_x[k] * slk0_1580[k]
                    - f_12 * pc_x[k] * slk1_1580[k];

        t_1581[k] = pb_x[k] * slk0_1581[k]
                    - f_12 * pc_x[k] * slk1_1581[k];

        t_1582[k] = f_13 * sli_1007[k]
                    + f_3 * pc_y[k] * smi_1231[k];

        t_1583[k] = pb_x[k] * slk0_1583[k]
                    - f_12 * pc_x[k] * slk1_1583[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, t_1587, pb_x, pc_x, pc_y, pc_z, slk0_1584, \
                         slk0_1587, sli_980, sli_1232, sli_1235, slk1_1584, slk1_1587, \
                         smi_1232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = pb_x[k] * slk0_1584[k]
                    + f_19 * sli_1232[k]
                    - f_12 * pc_x[k] * slk1_1584[k];

        t_1585[k] = f_3 * pc_y[k] * smi_1232[k];

        t_1586[k] = f_18 * sli_980[k]
                    + f_3 * pc_z[k] * smi_1232[k];

        t_1587[k] = pb_x[k] * slk0_1587[k]
                    + f_17 * sli_1235[k]
                    - f_12 * pc_x[k] * slk1_1587[k];
    }

#pragma omp simd aligned(t_1588, t_1589, t_1590, pb_x, pc_x, pc_y, slk0_1589, slk0_1590, \
                         sli_1237, sli_1238, slk1_1589, slk1_1590, \
                         smi_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1588[k] = f_3 * pc_y[k] * smi_1234[k];

        t_1589[k] = pb_x[k] * slk0_1589[k]
                    + f_17 * sli_1237[k]
                    - f_12 * pc_x[k] * slk1_1589[k];

        t_1590[k] = pb_x[k] * slk0_1590[k]
                    + f_16 * sli_1238[k]
                    - f_12 * pc_x[k] * slk1_1590[k];
    }

#pragma omp simd aligned(t_1591, t_1592, t_1593, pb_x, pc_x, pc_y, pc_z, slk0_1593, sli_983, \
                         sli_1241, slk1_1593, smi_1235, smi_1237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1591[k] = f_18 * sli_983[k]
                    + f_3 * pc_z[k] * smi_1235[k];

        t_1592[k] = f_3 * pc_y[k] * smi_1237[k];

        t_1593[k] = pb_x[k] * slk0_1593[k]
                    + f_16 * sli_1241[k]
                    - f_12 * pc_x[k] * slk1_1593[k];
    }

#pragma omp simd aligned(t_1594, t_1595, t_1596, pb_x, pc_x, pc_z, slk0_1594, slk0_1596, \
                         sli_986, sli_1242, sli_1244, slk1_1594, slk1_1596, \
                         smi_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1594[k] = pb_x[k] * slk0_1594[k]
                    + f_15 * sli_1242[k]
                    - f_12 * pc_x[k] * slk1_1594[k];

        t_1595[k] = f_18 * sli_986[k]
                    + f_3 * pc_z[k] * smi_1238[k];

        t_1596[k] = pb_x[k] * slk0_1596[k]
                    + f_15 * sli_1244[k]
                    - f_12 * pc_x[k] * slk1_1596[k];
    }

#pragma omp simd aligned(t_1597, t_1598, t_1599, pb_x, pc_x, pc_y, slk0_1598, slk0_1599, \
                         sli_1246, sli_1247, slk1_1598, slk1_1599, \
                         smi_1241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1597[k] = f_3 * pc_y[k] * smi_1241[k];

        t_1598[k] = pb_x[k] * slk0_1598[k]
                    + f_15 * sli_1246[k]
                    - f_12 * pc_x[k] * slk1_1598[k];

        t_1599[k] = pb_x[k] * slk0_1599[k]
                    + f_14 * sli_1247[k]
                    - f_12 * pc_x[k] * slk1_1599[k];
    }

#pragma omp simd aligned(t_1600, t_1601, t_1602, pb_x, pc_x, pc_z, slk0_1601, slk0_1602, \
                         sli_990, sli_1249, sli_1250, slk1_1601, slk1_1602, \
                         smi_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1600[k] = f_18 * sli_990[k]
                    + f_3 * pc_z[k] * smi_1242[k];

        t_1601[k] = pb_x[k] * slk0_1601[k]
                    + f_14 * sli_1249[k]
                    - f_12 * pc_x[k] * slk1_1601[k];

        t_1602[k] = pb_x[k] * slk0_1602[k]
                    + f_14 * sli_1250[k]
                    - f_12 * pc_x[k] * slk1_1602[k];
    }

#pragma omp simd aligned(t_1603, t_1604, t_1605, t_1606, pb_x, pc_x, pc_y, slk0_1604, \
                         sli_1252, sli_1253, sli_1254, slk1_1604, smi_1246, smi_1253, \
                         smi_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1603[k] = f_3 * pc_y[k] * smi_1246[k];

        t_1604[k] = pb_x[k] * slk0_1604[k]
                    + f_14 * sli_1252[k]
                    - f_12 * pc_x[k] * slk1_1604[k];

        t_1605[k] = f_13 * sli_1253[k]
                    + f_3 * pc_x[k] * smi_1253[k];

        t_1606[k] = f_13 * sli_1254[k]
                    + f_3 * pc_x[k] * smi_1254[k];
    }

#pragma omp simd aligned(t_1607, t_1608, t_1609, t_1610, t_1611, pc_x, sli_1255, sli_1256, \
                         sli_1257, sli_1258, sli_1259, smi_1255, smi_1256, smi_1257, smi_1258, \
                         smi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1607[k] = f_13 * sli_1255[k]
                    + f_3 * pc_x[k] * smi_1255[k];

        t_1608[k] = f_13 * sli_1256[k]
                    + f_3 * pc_x[k] * smi_1256[k];

        t_1609[k] = f_13 * sli_1257[k]
                    + f_3 * pc_x[k] * smi_1257[k];

        t_1610[k] = f_13 * sli_1258[k]
                    + f_3 * pc_x[k] * smi_1258[k];

        t_1611[k] = f_13 * sli_1259[k]
                    + f_3 * pc_x[k] * smi_1259[k];
    }

#pragma omp simd aligned(t_1612, t_1613, t_1614, t_1615, pb_x, pc_x, pc_z, slk0_1612, \
                         slk0_1614, slk0_1615, sli_1001, slk1_1612, slk1_1614, slk1_1615, \
                         smi_1253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1612[k] = pb_x[k] * slk0_1612[k]
                    - f_12 * pc_x[k] * slk1_1612[k];

        t_1613[k] = f_18 * sli_1001[k]
                    + f_3 * pc_z[k] * smi_1253[k];

        t_1614[k] = pb_x[k] * slk0_1614[k]
                    - f_12 * pc_x[k] * slk1_1614[k];

        t_1615[k] = pb_x[k] * slk0_1615[k]
                    - f_12 * pc_x[k] * slk1_1615[k];
    }

#pragma omp simd aligned(t_1616, t_1617, t_1618, t_1619, pb_x, pc_x, pc_y, slk0_1616, \
                         slk0_1617, slk0_1619, slk1_1616, slk1_1617, slk1_1619, \
                         smi_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1616[k] = pb_x[k] * slk0_1616[k]
                    - f_12 * pc_x[k] * slk1_1616[k];

        t_1617[k] = pb_x[k] * slk0_1617[k]
                    - f_12 * pc_x[k] * slk1_1617[k];

        t_1618[k] = f_3 * pc_y[k] * smi_1259[k];

        t_1619[k] = pb_x[k] * slk0_1619[k]
                    - f_12 * pc_x[k] * slk1_1619[k];
    }

#pragma omp simd aligned(t_1620, t_1621, t_1622, t_1623, pc_x, pc_y, pc_z, sli_1008, smh0_945, \
                         smh0_948, smh1_945, smh1_948, smi_1260, \
                         smi_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1620[k] = f_1 * smh0_945[k]
                    - f_2 * smh1_945[k]
                    + f_3 * pc_x[k] * smi_1260[k];

        t_1621[k] = f_0 * sli_1008[k]
                    + f_3 * pc_y[k] * smi_1260[k];

        t_1622[k] = f_3 * pc_z[k] * smi_1260[k];

        t_1623[k] = f_4 * smh0_948[k]
                    - f_5 * smh1_948[k]
                    + f_3 * pc_x[k] * smi_1263[k];
    }

#pragma omp simd aligned(t_1624, t_1625, t_1626, t_1627, pc_x, pc_y, pc_z, sli_1010, smh0_950, \
                         smh0_951, smh1_950, smh1_951, smi_1262, smi_1263, smi_1265, \
                         smi_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1624[k] = f_0 * sli_1010[k]
                    + f_3 * pc_y[k] * smi_1262[k];

        t_1625[k] = f_4 * smh0_950[k]
                    - f_5 * smh1_950[k]
                    + f_3 * pc_x[k] * smi_1265[k];

        t_1626[k] = f_6 * smh0_951[k]
                    - f_7 * smh1_951[k]
                    + f_3 * pc_x[k] * smi_1266[k];

        t_1627[k] = f_3 * pc_z[k] * smi_1263[k];
    }

#pragma omp simd aligned(t_1628, t_1629, t_1630, t_1631, pc_x, pc_y, pc_z, sli_1013, smh0_954, \
                         smh0_955, smh1_954, smh1_955, smi_1265, smi_1266, smi_1269, \
                         smi_1270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1628[k] = f_0 * sli_1013[k]
                    + f_3 * pc_y[k] * smi_1265[k];

        t_1629[k] = f_6 * smh0_954[k]
                    - f_7 * smh1_954[k]
                    + f_3 * pc_x[k] * smi_1269[k];

        t_1630[k] = f_8 * smh0_955[k]
                    - f_9 * smh1_955[k]
                    + f_3 * pc_x[k] * smi_1270[k];

        t_1631[k] = f_3 * pc_z[k] * smi_1266[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, pc_x, pc_y, sli_1017, smh0_957, smh0_959, \
                         smh1_957, smh1_959, smi_1269, smi_1272, \
                         smi_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_8 * smh0_957[k]
                    - f_9 * smh1_957[k]
                    + f_3 * pc_x[k] * smi_1272[k];

        t_1633[k] = f_0 * sli_1017[k]
                    + f_3 * pc_y[k] * smi_1269[k];

        t_1634[k] = f_8 * smh0_959[k]
                    - f_9 * smh1_959[k]
                    + f_3 * pc_x[k] * smi_1274[k];
    }

#pragma omp simd aligned(t_1635, t_1636, t_1637, t_1638, pc_x, pc_z, smh0_960, smh0_962, \
                         smh0_963, smh1_960, smh1_962, smh1_963, smi_1270, smi_1275, smi_1277, \
                         smi_1278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1635[k] = f_10 * smh0_960[k]
                    - f_11 * smh1_960[k]
                    + f_3 * pc_x[k] * smi_1275[k];

        t_1636[k] = f_3 * pc_z[k] * smi_1270[k];

        t_1637[k] = f_10 * smh0_962[k]
                    - f_11 * smh1_962[k]
                    + f_3 * pc_x[k] * smi_1277[k];

        t_1638[k] = f_10 * smh0_963[k]
                    - f_11 * smh1_963[k]
                    + f_3 * pc_x[k] * smi_1278[k];
    }

#pragma omp simd aligned(t_1639, t_1640, t_1641, t_1642, t_1643, pc_x, pc_y, sli_1022, \
                         smh0_965, smh1_965, smi_1274, smi_1280, smi_1281, smi_1282, \
                         smi_1283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1639[k] = f_0 * sli_1022[k]
                    + f_3 * pc_y[k] * smi_1274[k];

        t_1640[k] = f_10 * smh0_965[k]
                    - f_11 * smh1_965[k]
                    + f_3 * pc_x[k] * smi_1280[k];

        t_1641[k] = f_3 * pc_x[k] * smi_1281[k];

        t_1642[k] = f_3 * pc_x[k] * smi_1282[k];

        t_1643[k] = f_3 * pc_x[k] * smi_1283[k];
    }

#pragma omp simd aligned(t_1644, t_1645, t_1646, t_1647, t_1648, pc_x, pc_y, sli_1029, \
                         smh0_960, smh1_960, smi_1281, smi_1284, smi_1285, smi_1286, \
                         smi_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1644[k] = f_3 * pc_x[k] * smi_1284[k];

        t_1645[k] = f_3 * pc_x[k] * smi_1285[k];

        t_1646[k] = f_3 * pc_x[k] * smi_1286[k];

        t_1647[k] = f_3 * pc_x[k] * smi_1287[k];

        t_1648[k] = f_0 * sli_1029[k]
                    + f_1 * smh0_960[k]
                    - f_2 * smh1_960[k]
                    + f_3 * pc_y[k] * smi_1281[k];
    }

#pragma omp simd aligned(t_1649, t_1650, t_1651, pc_y, pc_z, sli_1031, sli_1032, smh0_962, \
                         smh0_963, smh1_962, smh1_963, smi_1281, smi_1283, \
                         smi_1284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1649[k] = f_3 * pc_z[k] * smi_1281[k];

        t_1650[k] = f_0 * sli_1031[k]
                    + f_4 * smh0_962[k]
                    - f_5 * smh1_962[k]
                    + f_3 * pc_y[k] * smi_1283[k];

        t_1651[k] = f_0 * sli_1032[k]
                    + f_6 * smh0_963[k]
                    - f_7 * smh1_963[k]
                    + f_3 * pc_y[k] * smi_1284[k];
    }

#pragma omp simd aligned(t_1652, t_1653, t_1654, t_1655, pc_y, pc_z, sli_1033, sli_1034, \
                         sli_1035, smh0_964, smh0_965, smh1_964, smh1_965, smi_1285, smi_1286, \
                         smi_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1652[k] = f_0 * sli_1033[k]
                    + f_8 * smh0_964[k]
                    - f_9 * smh1_964[k]
                    + f_3 * pc_y[k] * smi_1285[k];

        t_1653[k] = f_0 * sli_1034[k]
                    + f_10 * smh0_965[k]
                    - f_11 * smh1_965[k]
                    + f_3 * pc_y[k] * smi_1286[k];

        t_1654[k] = f_0 * sli_1035[k]
                    + f_3 * pc_y[k] * smi_1287[k];

        t_1655[k] = f_1 * smh0_965[k]
                    - f_2 * smh1_965[k]
                    + f_3 * pc_z[k] * smi_1287[k];
    }

#pragma omp simd aligned(t_1656, t_1657, t_1658, t_1659, pb_z, pc_y, pc_z, slk0_1296, \
                         slk0_1299, sli_1008, sli_1036, slk1_1296, slk1_1299, \
                         smi_1288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1656[k] = pb_z[k] * slk0_1296[k]
                    - f_12 * pc_z[k] * slk1_1296[k];

        t_1657[k] = f_18 * sli_1036[k]
                    + f_3 * pc_y[k] * smi_1288[k];

        t_1658[k] = f_13 * sli_1008[k]
                    + f_3 * pc_z[k] * smi_1288[k];

        t_1659[k] = pb_z[k] * slk0_1299[k]
                    - f_12 * pc_z[k] * slk1_1299[k];
    }

#pragma omp simd aligned(t_1660, t_1661, t_1662, pb_z, pc_x, pc_y, pc_z, slk0_1302, sli_1038, \
                         slk1_1302, smh0_971, smh1_971, smi_1290, \
                         smi_1293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1660[k] = f_18 * sli_1038[k]
                    + f_3 * pc_y[k] * smi_1290[k];

        t_1661[k] = f_4 * smh0_971[k]
                    - f_5 * smh1_971[k]
                    + f_3 * pc_x[k] * smi_1293[k];

        t_1662[k] = pb_z[k] * slk0_1302[k]
                    - f_12 * pc_z[k] * slk1_1302[k];
    }

#pragma omp simd aligned(t_1663, t_1664, t_1665, pc_x, pc_y, pc_z, sli_1011, sli_1041, \
                         smh0_975, smh1_975, smi_1291, smi_1293, \
                         smi_1297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1663[k] = f_13 * sli_1011[k]
                    + f_3 * pc_z[k] * smi_1291[k];

        t_1664[k] = f_18 * sli_1041[k]
                    + f_3 * pc_y[k] * smi_1293[k];

        t_1665[k] = f_6 * smh0_975[k]
                    - f_7 * smh1_975[k]
                    + f_3 * pc_x[k] * smi_1297[k];
    }

#pragma omp simd aligned(t_1666, t_1667, t_1668, pb_z, pc_x, pc_z, slk0_1306, sli_1014, \
                         slk1_1306, smh0_978, smh1_978, smi_1294, \
                         smi_1300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1666[k] = pb_z[k] * slk0_1306[k]
                    - f_12 * pc_z[k] * slk1_1306[k];

        t_1667[k] = f_13 * sli_1014[k]
                    + f_3 * pc_z[k] * smi_1294[k];

        t_1668[k] = f_8 * smh0_978[k]
                    - f_9 * smh1_978[k]
                    + f_3 * pc_x[k] * smi_1300[k];
    }

#pragma omp simd aligned(t_1669, t_1670, t_1671, pb_z, pc_x, pc_y, pc_z, slk0_1311, sli_1045, \
                         slk1_1311, smh0_980, smh1_980, smi_1297, \
                         smi_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1669[k] = f_18 * sli_1045[k]
                    + f_3 * pc_y[k] * smi_1297[k];

        t_1670[k] = f_8 * smh0_980[k]
                    - f_9 * smh1_980[k]
                    + f_3 * pc_x[k] * smi_1302[k];

        t_1671[k] = pb_z[k] * slk0_1311[k]
                    - f_12 * pc_z[k] * slk1_1311[k];
    }

#pragma omp simd aligned(t_1672, t_1673, t_1674, pc_x, pc_z, sli_1018, smh0_983, smh0_984, \
                         smh1_983, smh1_984, smi_1298, smi_1305, \
                         smi_1306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1672[k] = f_13 * sli_1018[k]
                    + f_3 * pc_z[k] * smi_1298[k];

        t_1673[k] = f_10 * smh0_983[k]
                    - f_11 * smh1_983[k]
                    + f_3 * pc_x[k] * smi_1305[k];

        t_1674[k] = f_10 * smh0_984[k]
                    - f_11 * smh1_984[k]
                    + f_3 * pc_x[k] * smi_1306[k];
    }

#pragma omp simd aligned(t_1675, t_1676, t_1677, t_1678, t_1679, pc_x, pc_y, sli_1050, \
                         smh0_986, smh1_986, smi_1302, smi_1308, smi_1309, smi_1310, \
                         smi_1311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1675[k] = f_18 * sli_1050[k]
                    + f_3 * pc_y[k] * smi_1302[k];

        t_1676[k] = f_10 * smh0_986[k]
                    - f_11 * smh1_986[k]
                    + f_3 * pc_x[k] * smi_1308[k];

        t_1677[k] = f_3 * pc_x[k] * smi_1309[k];

        t_1678[k] = f_3 * pc_x[k] * smi_1310[k];

        t_1679[k] = f_3 * pc_x[k] * smi_1311[k];
    }

#pragma omp simd aligned(t_1680, t_1681, t_1682, t_1683, t_1684, pb_z, pc_x, pc_z, slk0_1324, \
                         slk1_1324, smi_1312, smi_1313, smi_1314, \
                         smi_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1680[k] = f_3 * pc_x[k] * smi_1312[k];

        t_1681[k] = f_3 * pc_x[k] * smi_1313[k];

        t_1682[k] = f_3 * pc_x[k] * smi_1314[k];

        t_1683[k] = f_3 * pc_x[k] * smi_1315[k];

        t_1684[k] = pb_z[k] * slk0_1324[k]
                    - f_12 * pc_z[k] * slk1_1324[k];
    }

#pragma omp simd aligned(t_1685, t_1686, t_1687, pb_z, pc_z, slk0_1326, slk0_1327, sli_1029, \
                         sli_1030, sli_1031, slk1_1326, slk1_1327, \
                         smi_1309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1685[k] = f_13 * sli_1029[k]
                    + f_3 * pc_z[k] * smi_1309[k];

        t_1686[k] = pb_z[k] * slk0_1326[k]
                    + f_14 * sli_1030[k]
                    - f_12 * pc_z[k] * slk1_1326[k];

        t_1687[k] = pb_z[k] * slk0_1327[k]
                    + f_15 * sli_1031[k]
                    - f_12 * pc_z[k] * slk1_1327[k];
    }

#pragma omp simd aligned(t_1688, t_1689, t_1690, pb_z, pc_y, pc_z, slk0_1328, slk0_1329, \
                         sli_1032, sli_1033, sli_1063, slk1_1328, slk1_1329, \
                         smi_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1688[k] = pb_z[k] * slk0_1328[k]
                    + f_16 * sli_1032[k]
                    - f_12 * pc_z[k] * slk1_1328[k];

        t_1689[k] = pb_z[k] * slk0_1329[k]
                    + f_17 * sli_1033[k]
                    - f_12 * pc_z[k] * slk1_1329[k];

        t_1690[k] = f_18 * sli_1063[k]
                    + f_3 * pc_y[k] * smi_1315[k];
    }

#pragma omp simd aligned(t_1691, t_1692, t_1693, t_1694, pc_x, pc_y, pc_z, sli_1035, sli_1036, \
                         sli_1064, smh0_986, smh0_987, smh1_986, smh1_987, smi_1315, \
                         smi_1316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1691[k] = f_13 * sli_1035[k]
                    + f_1 * smh0_986[k]
                    - f_2 * smh1_986[k]
                    + f_3 * pc_z[k] * smi_1315[k];

        t_1692[k] = f_1 * smh0_987[k]
                    - f_2 * smh1_987[k]
                    + f_3 * pc_x[k] * smi_1316[k];

        t_1693[k] = f_19 * sli_1064[k]
                    + f_3 * pc_y[k] * smi_1316[k];

        t_1694[k] = f_14 * sli_1036[k]
                    + f_3 * pc_z[k] * smi_1316[k];
    }

#pragma omp simd aligned(t_1695, t_1696, t_1697, pc_x, pc_y, sli_1066, smh0_990, smh0_992, \
                         smh1_990, smh1_992, smi_1318, smi_1319, \
                         smi_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1695[k] = f_4 * smh0_990[k]
                    - f_5 * smh1_990[k]
                    + f_3 * pc_x[k] * smi_1319[k];

        t_1696[k] = f_19 * sli_1066[k]
                    + f_3 * pc_y[k] * smi_1318[k];

        t_1697[k] = f_4 * smh0_992[k]
                    - f_5 * smh1_992[k]
                    + f_3 * pc_x[k] * smi_1321[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t sli, const size_t smh0,
                                                           const size_t smh1, const size_t smi,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sli_1039 = buffer.data(sli + 1039);
    const auto *sli_1042 = buffer.data(sli + 1042);
    const auto *sli_1046 = buffer.data(sli + 1046);
    const auto *sli_1057 = buffer.data(sli + 1057);
    const auto *sli_1063 = buffer.data(sli + 1063);
    const auto *sli_1064 = buffer.data(sli + 1064);
    const auto *sli_1067 = buffer.data(sli + 1067);
    const auto *sli_1069 = buffer.data(sli + 1069);
    const auto *sli_1070 = buffer.data(sli + 1070);
    const auto *sli_1073 = buffer.data(sli + 1073);
    const auto *sli_1074 = buffer.data(sli + 1074);
    const auto *sli_1078 = buffer.data(sli + 1078);
    const auto *sli_1085 = buffer.data(sli + 1085);
    const auto *sli_1087 = buffer.data(sli + 1087);
    const auto *sli_1088 = buffer.data(sli + 1088);
    const auto *sli_1089 = buffer.data(sli + 1089);
    const auto *sli_1090 = buffer.data(sli + 1090);
    const auto *sli_1091 = buffer.data(sli + 1091);
    const auto *sli_1092 = buffer.data(sli + 1092);
    const auto *sli_1094 = buffer.data(sli + 1094);
    const auto *sli_1095 = buffer.data(sli + 1095);
    const auto *sli_1097 = buffer.data(sli + 1097);
    const auto *sli_1098 = buffer.data(sli + 1098);
    const auto *sli_1101 = buffer.data(sli + 1101);
    const auto *sli_1102 = buffer.data(sli + 1102);
    const auto *sli_1106 = buffer.data(sli + 1106);
    const auto *sli_1113 = buffer.data(sli + 1113);
    const auto *sli_1115 = buffer.data(sli + 1115);
    const auto *sli_1116 = buffer.data(sli + 1116);
    const auto *sli_1117 = buffer.data(sli + 1117);
    const auto *sli_1118 = buffer.data(sli + 1118);
    const auto *sli_1119 = buffer.data(sli + 1119);
    const auto *sli_1120 = buffer.data(sli + 1120);
    const auto *sli_1122 = buffer.data(sli + 1122);
    const auto *sli_1123 = buffer.data(sli + 1123);
    const auto *sli_1125 = buffer.data(sli + 1125);
    const auto *sli_1126 = buffer.data(sli + 1126);
    const auto *sli_1129 = buffer.data(sli + 1129);
    const auto *sli_1130 = buffer.data(sli + 1130);
    const auto *sli_1134 = buffer.data(sli + 1134);
    const auto *sli_1141 = buffer.data(sli + 1141);
    const auto *sli_1143 = buffer.data(sli + 1143);
    const auto *sli_1144 = buffer.data(sli + 1144);
    const auto *sli_1145 = buffer.data(sli + 1145);
    const auto *sli_1146 = buffer.data(sli + 1146);
    const auto *sli_1147 = buffer.data(sli + 1147);
    const auto *sli_1148 = buffer.data(sli + 1148);
    const auto *sli_1150 = buffer.data(sli + 1150);
    const auto *sli_1153 = buffer.data(sli + 1153);
    const auto *sli_1157 = buffer.data(sli + 1157);

    const auto *smh0_993 = buffer.data(smh0 + 993);
    const auto *smh0_996 = buffer.data(smh0 + 996);
    const auto *smh0_997 = buffer.data(smh0 + 997);
    const auto *smh0_999 = buffer.data(smh0 + 999);
    const auto *smh0_1001 = buffer.data(smh0 + 1001);
    const auto *smh0_1002 = buffer.data(smh0 + 1002);
    const auto *smh0_1004 = buffer.data(smh0 + 1004);
    const auto *smh0_1005 = buffer.data(smh0 + 1005);
    const auto *smh0_1006 = buffer.data(smh0 + 1006);
    const auto *smh0_1007 = buffer.data(smh0 + 1007);
    const auto *smh0_1008 = buffer.data(smh0 + 1008);
    const auto *smh0_1011 = buffer.data(smh0 + 1011);
    const auto *smh0_1013 = buffer.data(smh0 + 1013);
    const auto *smh0_1014 = buffer.data(smh0 + 1014);
    const auto *smh0_1017 = buffer.data(smh0 + 1017);
    const auto *smh0_1018 = buffer.data(smh0 + 1018);
    const auto *smh0_1020 = buffer.data(smh0 + 1020);
    const auto *smh0_1022 = buffer.data(smh0 + 1022);
    const auto *smh0_1023 = buffer.data(smh0 + 1023);
    const auto *smh0_1025 = buffer.data(smh0 + 1025);
    const auto *smh0_1026 = buffer.data(smh0 + 1026);
    const auto *smh0_1027 = buffer.data(smh0 + 1027);
    const auto *smh0_1028 = buffer.data(smh0 + 1028);
    const auto *smh0_1029 = buffer.data(smh0 + 1029);
    const auto *smh0_1032 = buffer.data(smh0 + 1032);
    const auto *smh0_1034 = buffer.data(smh0 + 1034);
    const auto *smh0_1035 = buffer.data(smh0 + 1035);
    const auto *smh0_1038 = buffer.data(smh0 + 1038);
    const auto *smh0_1039 = buffer.data(smh0 + 1039);
    const auto *smh0_1041 = buffer.data(smh0 + 1041);
    const auto *smh0_1043 = buffer.data(smh0 + 1043);
    const auto *smh0_1044 = buffer.data(smh0 + 1044);
    const auto *smh0_1046 = buffer.data(smh0 + 1046);
    const auto *smh0_1047 = buffer.data(smh0 + 1047);
    const auto *smh0_1048 = buffer.data(smh0 + 1048);
    const auto *smh0_1049 = buffer.data(smh0 + 1049);
    const auto *smh0_1050 = buffer.data(smh0 + 1050);
    const auto *smh0_1053 = buffer.data(smh0 + 1053);
    const auto *smh0_1055 = buffer.data(smh0 + 1055);
    const auto *smh0_1056 = buffer.data(smh0 + 1056);
    const auto *smh0_1059 = buffer.data(smh0 + 1059);
    const auto *smh0_1060 = buffer.data(smh0 + 1060);
    const auto *smh0_1062 = buffer.data(smh0 + 1062);
    const auto *smh0_1064 = buffer.data(smh0 + 1064);
    const auto *smh0_1065 = buffer.data(smh0 + 1065);
    const auto *smh0_1067 = buffer.data(smh0 + 1067);

    const auto *smh1_993 = buffer.data(smh1 + 993);
    const auto *smh1_996 = buffer.data(smh1 + 996);
    const auto *smh1_997 = buffer.data(smh1 + 997);
    const auto *smh1_999 = buffer.data(smh1 + 999);
    const auto *smh1_1001 = buffer.data(smh1 + 1001);
    const auto *smh1_1002 = buffer.data(smh1 + 1002);
    const auto *smh1_1004 = buffer.data(smh1 + 1004);
    const auto *smh1_1005 = buffer.data(smh1 + 1005);
    const auto *smh1_1006 = buffer.data(smh1 + 1006);
    const auto *smh1_1007 = buffer.data(smh1 + 1007);
    const auto *smh1_1008 = buffer.data(smh1 + 1008);
    const auto *smh1_1011 = buffer.data(smh1 + 1011);
    const auto *smh1_1013 = buffer.data(smh1 + 1013);
    const auto *smh1_1014 = buffer.data(smh1 + 1014);
    const auto *smh1_1017 = buffer.data(smh1 + 1017);
    const auto *smh1_1018 = buffer.data(smh1 + 1018);
    const auto *smh1_1020 = buffer.data(smh1 + 1020);
    const auto *smh1_1022 = buffer.data(smh1 + 1022);
    const auto *smh1_1023 = buffer.data(smh1 + 1023);
    const auto *smh1_1025 = buffer.data(smh1 + 1025);
    const auto *smh1_1026 = buffer.data(smh1 + 1026);
    const auto *smh1_1027 = buffer.data(smh1 + 1027);
    const auto *smh1_1028 = buffer.data(smh1 + 1028);
    const auto *smh1_1029 = buffer.data(smh1 + 1029);
    const auto *smh1_1032 = buffer.data(smh1 + 1032);
    const auto *smh1_1034 = buffer.data(smh1 + 1034);
    const auto *smh1_1035 = buffer.data(smh1 + 1035);
    const auto *smh1_1038 = buffer.data(smh1 + 1038);
    const auto *smh1_1039 = buffer.data(smh1 + 1039);
    const auto *smh1_1041 = buffer.data(smh1 + 1041);
    const auto *smh1_1043 = buffer.data(smh1 + 1043);
    const auto *smh1_1044 = buffer.data(smh1 + 1044);
    const auto *smh1_1046 = buffer.data(smh1 + 1046);
    const auto *smh1_1047 = buffer.data(smh1 + 1047);
    const auto *smh1_1048 = buffer.data(smh1 + 1048);
    const auto *smh1_1049 = buffer.data(smh1 + 1049);
    const auto *smh1_1050 = buffer.data(smh1 + 1050);
    const auto *smh1_1053 = buffer.data(smh1 + 1053);
    const auto *smh1_1055 = buffer.data(smh1 + 1055);
    const auto *smh1_1056 = buffer.data(smh1 + 1056);
    const auto *smh1_1059 = buffer.data(smh1 + 1059);
    const auto *smh1_1060 = buffer.data(smh1 + 1060);
    const auto *smh1_1062 = buffer.data(smh1 + 1062);
    const auto *smh1_1064 = buffer.data(smh1 + 1064);
    const auto *smh1_1065 = buffer.data(smh1 + 1065);
    const auto *smh1_1067 = buffer.data(smh1 + 1067);

    const auto *smi_1319 = buffer.data(smi + 1319);
    const auto *smi_1321 = buffer.data(smi + 1321);
    const auto *smi_1322 = buffer.data(smi + 1322);
    const auto *smi_1325 = buffer.data(smi + 1325);
    const auto *smi_1326 = buffer.data(smi + 1326);
    const auto *smi_1328 = buffer.data(smi + 1328);
    const auto *smi_1330 = buffer.data(smi + 1330);
    const auto *smi_1331 = buffer.data(smi + 1331);
    const auto *smi_1333 = buffer.data(smi + 1333);
    const auto *smi_1334 = buffer.data(smi + 1334);
    const auto *smi_1336 = buffer.data(smi + 1336);
    const auto *smi_1337 = buffer.data(smi + 1337);
    const auto *smi_1338 = buffer.data(smi + 1338);
    const auto *smi_1339 = buffer.data(smi + 1339);
    const auto *smi_1340 = buffer.data(smi + 1340);
    const auto *smi_1341 = buffer.data(smi + 1341);
    const auto *smi_1342 = buffer.data(smi + 1342);
    const auto *smi_1343 = buffer.data(smi + 1343);
    const auto *smi_1344 = buffer.data(smi + 1344);
    const auto *smi_1346 = buffer.data(smi + 1346);
    const auto *smi_1347 = buffer.data(smi + 1347);
    const auto *smi_1349 = buffer.data(smi + 1349);
    const auto *smi_1350 = buffer.data(smi + 1350);
    const auto *smi_1353 = buffer.data(smi + 1353);
    const auto *smi_1354 = buffer.data(smi + 1354);
    const auto *smi_1356 = buffer.data(smi + 1356);
    const auto *smi_1358 = buffer.data(smi + 1358);
    const auto *smi_1359 = buffer.data(smi + 1359);
    const auto *smi_1361 = buffer.data(smi + 1361);
    const auto *smi_1362 = buffer.data(smi + 1362);
    const auto *smi_1364 = buffer.data(smi + 1364);
    const auto *smi_1365 = buffer.data(smi + 1365);
    const auto *smi_1366 = buffer.data(smi + 1366);
    const auto *smi_1367 = buffer.data(smi + 1367);
    const auto *smi_1368 = buffer.data(smi + 1368);
    const auto *smi_1369 = buffer.data(smi + 1369);
    const auto *smi_1370 = buffer.data(smi + 1370);
    const auto *smi_1371 = buffer.data(smi + 1371);
    const auto *smi_1372 = buffer.data(smi + 1372);
    const auto *smi_1374 = buffer.data(smi + 1374);
    const auto *smi_1375 = buffer.data(smi + 1375);
    const auto *smi_1377 = buffer.data(smi + 1377);
    const auto *smi_1378 = buffer.data(smi + 1378);
    const auto *smi_1381 = buffer.data(smi + 1381);
    const auto *smi_1382 = buffer.data(smi + 1382);
    const auto *smi_1384 = buffer.data(smi + 1384);
    const auto *smi_1386 = buffer.data(smi + 1386);
    const auto *smi_1387 = buffer.data(smi + 1387);
    const auto *smi_1389 = buffer.data(smi + 1389);
    const auto *smi_1390 = buffer.data(smi + 1390);
    const auto *smi_1392 = buffer.data(smi + 1392);
    const auto *smi_1393 = buffer.data(smi + 1393);
    const auto *smi_1394 = buffer.data(smi + 1394);
    const auto *smi_1395 = buffer.data(smi + 1395);
    const auto *smi_1396 = buffer.data(smi + 1396);
    const auto *smi_1397 = buffer.data(smi + 1397);
    const auto *smi_1398 = buffer.data(smi + 1398);
    const auto *smi_1399 = buffer.data(smi + 1399);
    const auto *smi_1400 = buffer.data(smi + 1400);
    const auto *smi_1402 = buffer.data(smi + 1402);
    const auto *smi_1403 = buffer.data(smi + 1403);
    const auto *smi_1405 = buffer.data(smi + 1405);
    const auto *smi_1406 = buffer.data(smi + 1406);
    const auto *smi_1409 = buffer.data(smi + 1409);
    const auto *smi_1410 = buffer.data(smi + 1410);
    const auto *smi_1412 = buffer.data(smi + 1412);
    const auto *smi_1414 = buffer.data(smi + 1414);
    const auto *smi_1415 = buffer.data(smi + 1415);
    const auto *smi_1417 = buffer.data(smi + 1417);

#pragma omp simd aligned(t_1698, t_1699, t_1700, pc_x, pc_y, pc_z, sli_1039, sli_1069, \
                         smh0_993, smh1_993, smi_1319, smi_1321, \
                         smi_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_6 * smh0_993[k]
                    - f_7 * smh1_993[k]
                    + f_3 * pc_x[k] * smi_1322[k];

        t_1699[k] = f_14 * sli_1039[k]
                    + f_3 * pc_z[k] * smi_1319[k];

        t_1700[k] = f_19 * sli_1069[k]
                    + f_3 * pc_y[k] * smi_1321[k];
    }

#pragma omp simd aligned(t_1701, t_1702, t_1703, pc_x, pc_z, sli_1042, smh0_996, smh0_997, \
                         smh1_996, smh1_997, smi_1322, smi_1325, \
                         smi_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1701[k] = f_6 * smh0_996[k]
                    - f_7 * smh1_996[k]
                    + f_3 * pc_x[k] * smi_1325[k];

        t_1702[k] = f_8 * smh0_997[k]
                    - f_9 * smh1_997[k]
                    + f_3 * pc_x[k] * smi_1326[k];

        t_1703[k] = f_14 * sli_1042[k]
                    + f_3 * pc_z[k] * smi_1322[k];
    }

#pragma omp simd aligned(t_1704, t_1705, t_1706, pc_x, pc_y, sli_1073, smh0_999, smh0_1001, \
                         smh1_999, smh1_1001, smi_1325, smi_1328, \
                         smi_1330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1704[k] = f_8 * smh0_999[k]
                    - f_9 * smh1_999[k]
                    + f_3 * pc_x[k] * smi_1328[k];

        t_1705[k] = f_19 * sli_1073[k]
                    + f_3 * pc_y[k] * smi_1325[k];

        t_1706[k] = f_8 * smh0_1001[k]
                    - f_9 * smh1_1001[k]
                    + f_3 * pc_x[k] * smi_1330[k];
    }

#pragma omp simd aligned(t_1707, t_1708, t_1709, pc_x, pc_z, sli_1046, smh0_1002, smh0_1004, \
                         smh1_1002, smh1_1004, smi_1326, smi_1331, \
                         smi_1333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1707[k] = f_10 * smh0_1002[k]
                    - f_11 * smh1_1002[k]
                    + f_3 * pc_x[k] * smi_1331[k];

        t_1708[k] = f_14 * sli_1046[k]
                    + f_3 * pc_z[k] * smi_1326[k];

        t_1709[k] = f_10 * smh0_1004[k]
                    - f_11 * smh1_1004[k]
                    + f_3 * pc_x[k] * smi_1333[k];
    }

#pragma omp simd aligned(t_1710, t_1711, t_1712, t_1713, pc_x, pc_y, sli_1078, smh0_1005, \
                         smh0_1007, smh1_1005, smh1_1007, smi_1330, smi_1334, smi_1336, \
                         smi_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1710[k] = f_10 * smh0_1005[k]
                    - f_11 * smh1_1005[k]
                    + f_3 * pc_x[k] * smi_1334[k];

        t_1711[k] = f_19 * sli_1078[k]
                    + f_3 * pc_y[k] * smi_1330[k];

        t_1712[k] = f_10 * smh0_1007[k]
                    - f_11 * smh1_1007[k]
                    + f_3 * pc_x[k] * smi_1336[k];

        t_1713[k] = f_3 * pc_x[k] * smi_1337[k];
    }

#pragma omp simd aligned(t_1714, t_1715, t_1716, t_1717, t_1718, t_1719, pc_x, smi_1338, \
                         smi_1339, smi_1340, smi_1341, smi_1342, \
                         smi_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1714[k] = f_3 * pc_x[k] * smi_1338[k];

        t_1715[k] = f_3 * pc_x[k] * smi_1339[k];

        t_1716[k] = f_3 * pc_x[k] * smi_1340[k];

        t_1717[k] = f_3 * pc_x[k] * smi_1341[k];

        t_1718[k] = f_3 * pc_x[k] * smi_1342[k];

        t_1719[k] = f_3 * pc_x[k] * smi_1343[k];
    }

#pragma omp simd aligned(t_1720, t_1721, t_1722, pc_y, pc_z, sli_1057, sli_1085, sli_1087, \
                         smh0_1002, smh0_1004, smh1_1002, smh1_1004, smi_1337, \
                         smi_1339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1720[k] = f_19 * sli_1085[k]
                    + f_1 * smh0_1002[k]
                    - f_2 * smh1_1002[k]
                    + f_3 * pc_y[k] * smi_1337[k];

        t_1721[k] = f_14 * sli_1057[k]
                    + f_3 * pc_z[k] * smi_1337[k];

        t_1722[k] = f_19 * sli_1087[k]
                    + f_4 * smh0_1004[k]
                    - f_5 * smh1_1004[k]
                    + f_3 * pc_y[k] * smi_1339[k];
    }

#pragma omp simd aligned(t_1723, t_1724, t_1725, pc_y, sli_1088, sli_1089, sli_1090, \
                         smh0_1005, smh0_1006, smh0_1007, smh1_1005, smh1_1006, smh1_1007, \
                         smi_1340, smi_1341, smi_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1723[k] = f_19 * sli_1088[k]
                    + f_6 * smh0_1005[k]
                    - f_7 * smh1_1005[k]
                    + f_3 * pc_y[k] * smi_1340[k];

        t_1724[k] = f_19 * sli_1089[k]
                    + f_8 * smh0_1006[k]
                    - f_9 * smh1_1006[k]
                    + f_3 * pc_y[k] * smi_1341[k];

        t_1725[k] = f_19 * sli_1090[k]
                    + f_10 * smh0_1007[k]
                    - f_11 * smh1_1007[k]
                    + f_3 * pc_y[k] * smi_1342[k];
    }

#pragma omp simd aligned(t_1726, t_1727, t_1728, t_1729, pc_x, pc_y, pc_z, sli_1063, sli_1091, \
                         sli_1092, smh0_1007, smh0_1008, smh1_1007, smh1_1008, smi_1343, \
                         smi_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1726[k] = f_19 * sli_1091[k]
                    + f_3 * pc_y[k] * smi_1343[k];

        t_1727[k] = f_14 * sli_1063[k]
                    + f_1 * smh0_1007[k]
                    - f_2 * smh1_1007[k]
                    + f_3 * pc_z[k] * smi_1343[k];

        t_1728[k] = f_1 * smh0_1008[k]
                    - f_2 * smh1_1008[k]
                    + f_3 * pc_x[k] * smi_1344[k];

        t_1729[k] = f_20 * sli_1092[k]
                    + f_3 * pc_y[k] * smi_1344[k];
    }

#pragma omp simd aligned(t_1730, t_1731, t_1732, pc_x, pc_y, pc_z, sli_1064, sli_1094, \
                         smh0_1011, smh1_1011, smi_1344, smi_1346, \
                         smi_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1730[k] = f_15 * sli_1064[k]
                    + f_3 * pc_z[k] * smi_1344[k];

        t_1731[k] = f_4 * smh0_1011[k]
                    - f_5 * smh1_1011[k]
                    + f_3 * pc_x[k] * smi_1347[k];

        t_1732[k] = f_20 * sli_1094[k]
                    + f_3 * pc_y[k] * smi_1346[k];
    }

#pragma omp simd aligned(t_1733, t_1734, t_1735, t_1736, pc_x, pc_y, pc_z, sli_1067, sli_1097, \
                         smh0_1013, smh0_1014, smh1_1013, smh1_1014, smi_1347, smi_1349, \
                         smi_1350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1733[k] = f_4 * smh0_1013[k]
                    - f_5 * smh1_1013[k]
                    + f_3 * pc_x[k] * smi_1349[k];

        t_1734[k] = f_6 * smh0_1014[k]
                    - f_7 * smh1_1014[k]
                    + f_3 * pc_x[k] * smi_1350[k];

        t_1735[k] = f_15 * sli_1067[k]
                    + f_3 * pc_z[k] * smi_1347[k];

        t_1736[k] = f_20 * sli_1097[k]
                    + f_3 * pc_y[k] * smi_1349[k];
    }

#pragma omp simd aligned(t_1737, t_1738, t_1739, pc_x, pc_z, sli_1070, smh0_1017, smh0_1018, \
                         smh1_1017, smh1_1018, smi_1350, smi_1353, \
                         smi_1354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1737[k] = f_6 * smh0_1017[k]
                    - f_7 * smh1_1017[k]
                    + f_3 * pc_x[k] * smi_1353[k];

        t_1738[k] = f_8 * smh0_1018[k]
                    - f_9 * smh1_1018[k]
                    + f_3 * pc_x[k] * smi_1354[k];

        t_1739[k] = f_15 * sli_1070[k]
                    + f_3 * pc_z[k] * smi_1350[k];
    }

#pragma omp simd aligned(t_1740, t_1741, t_1742, pc_x, pc_y, sli_1101, smh0_1020, smh0_1022, \
                         smh1_1020, smh1_1022, smi_1353, smi_1356, \
                         smi_1358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1740[k] = f_8 * smh0_1020[k]
                    - f_9 * smh1_1020[k]
                    + f_3 * pc_x[k] * smi_1356[k];

        t_1741[k] = f_20 * sli_1101[k]
                    + f_3 * pc_y[k] * smi_1353[k];

        t_1742[k] = f_8 * smh0_1022[k]
                    - f_9 * smh1_1022[k]
                    + f_3 * pc_x[k] * smi_1358[k];
    }

#pragma omp simd aligned(t_1743, t_1744, t_1745, pc_x, pc_z, sli_1074, smh0_1023, smh0_1025, \
                         smh1_1023, smh1_1025, smi_1354, smi_1359, \
                         smi_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1743[k] = f_10 * smh0_1023[k]
                    - f_11 * smh1_1023[k]
                    + f_3 * pc_x[k] * smi_1359[k];

        t_1744[k] = f_15 * sli_1074[k]
                    + f_3 * pc_z[k] * smi_1354[k];

        t_1745[k] = f_10 * smh0_1025[k]
                    - f_11 * smh1_1025[k]
                    + f_3 * pc_x[k] * smi_1361[k];
    }

#pragma omp simd aligned(t_1746, t_1747, t_1748, t_1749, pc_x, pc_y, sli_1106, smh0_1026, \
                         smh0_1028, smh1_1026, smh1_1028, smi_1358, smi_1362, smi_1364, \
                         smi_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1746[k] = f_10 * smh0_1026[k]
                    - f_11 * smh1_1026[k]
                    + f_3 * pc_x[k] * smi_1362[k];

        t_1747[k] = f_20 * sli_1106[k]
                    + f_3 * pc_y[k] * smi_1358[k];

        t_1748[k] = f_10 * smh0_1028[k]
                    - f_11 * smh1_1028[k]
                    + f_3 * pc_x[k] * smi_1364[k];

        t_1749[k] = f_3 * pc_x[k] * smi_1365[k];
    }

#pragma omp simd aligned(t_1750, t_1751, t_1752, t_1753, t_1754, t_1755, pc_x, smi_1366, \
                         smi_1367, smi_1368, smi_1369, smi_1370, \
                         smi_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1750[k] = f_3 * pc_x[k] * smi_1366[k];

        t_1751[k] = f_3 * pc_x[k] * smi_1367[k];

        t_1752[k] = f_3 * pc_x[k] * smi_1368[k];

        t_1753[k] = f_3 * pc_x[k] * smi_1369[k];

        t_1754[k] = f_3 * pc_x[k] * smi_1370[k];

        t_1755[k] = f_3 * pc_x[k] * smi_1371[k];
    }

#pragma omp simd aligned(t_1756, t_1757, t_1758, pc_y, pc_z, sli_1085, sli_1113, sli_1115, \
                         smh0_1023, smh0_1025, smh1_1023, smh1_1025, smi_1365, \
                         smi_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1756[k] = f_20 * sli_1113[k]
                    + f_1 * smh0_1023[k]
                    - f_2 * smh1_1023[k]
                    + f_3 * pc_y[k] * smi_1365[k];

        t_1757[k] = f_15 * sli_1085[k]
                    + f_3 * pc_z[k] * smi_1365[k];

        t_1758[k] = f_20 * sli_1115[k]
                    + f_4 * smh0_1025[k]
                    - f_5 * smh1_1025[k]
                    + f_3 * pc_y[k] * smi_1367[k];
    }

#pragma omp simd aligned(t_1759, t_1760, t_1761, pc_y, sli_1116, sli_1117, sli_1118, \
                         smh0_1026, smh0_1027, smh0_1028, smh1_1026, smh1_1027, smh1_1028, \
                         smi_1368, smi_1369, smi_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1759[k] = f_20 * sli_1116[k]
                    + f_6 * smh0_1026[k]
                    - f_7 * smh1_1026[k]
                    + f_3 * pc_y[k] * smi_1368[k];

        t_1760[k] = f_20 * sli_1117[k]
                    + f_8 * smh0_1027[k]
                    - f_9 * smh1_1027[k]
                    + f_3 * pc_y[k] * smi_1369[k];

        t_1761[k] = f_20 * sli_1118[k]
                    + f_10 * smh0_1028[k]
                    - f_11 * smh1_1028[k]
                    + f_3 * pc_y[k] * smi_1370[k];
    }

#pragma omp simd aligned(t_1762, t_1763, t_1764, t_1765, pc_x, pc_y, pc_z, sli_1091, sli_1119, \
                         sli_1120, smh0_1028, smh0_1029, smh1_1028, smh1_1029, smi_1371, \
                         smi_1372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1762[k] = f_20 * sli_1119[k]
                    + f_3 * pc_y[k] * smi_1371[k];

        t_1763[k] = f_15 * sli_1091[k]
                    + f_1 * smh0_1028[k]
                    - f_2 * smh1_1028[k]
                    + f_3 * pc_z[k] * smi_1371[k];

        t_1764[k] = f_1 * smh0_1029[k]
                    - f_2 * smh1_1029[k]
                    + f_3 * pc_x[k] * smi_1372[k];

        t_1765[k] = f_17 * sli_1120[k]
                    + f_3 * pc_y[k] * smi_1372[k];
    }

#pragma omp simd aligned(t_1766, t_1767, t_1768, pc_x, pc_y, pc_z, sli_1092, sli_1122, \
                         smh0_1032, smh1_1032, smi_1372, smi_1374, \
                         smi_1375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1766[k] = f_16 * sli_1092[k]
                    + f_3 * pc_z[k] * smi_1372[k];

        t_1767[k] = f_4 * smh0_1032[k]
                    - f_5 * smh1_1032[k]
                    + f_3 * pc_x[k] * smi_1375[k];

        t_1768[k] = f_17 * sli_1122[k]
                    + f_3 * pc_y[k] * smi_1374[k];
    }

#pragma omp simd aligned(t_1769, t_1770, t_1771, t_1772, pc_x, pc_y, pc_z, sli_1095, sli_1125, \
                         smh0_1034, smh0_1035, smh1_1034, smh1_1035, smi_1375, smi_1377, \
                         smi_1378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1769[k] = f_4 * smh0_1034[k]
                    - f_5 * smh1_1034[k]
                    + f_3 * pc_x[k] * smi_1377[k];

        t_1770[k] = f_6 * smh0_1035[k]
                    - f_7 * smh1_1035[k]
                    + f_3 * pc_x[k] * smi_1378[k];

        t_1771[k] = f_16 * sli_1095[k]
                    + f_3 * pc_z[k] * smi_1375[k];

        t_1772[k] = f_17 * sli_1125[k]
                    + f_3 * pc_y[k] * smi_1377[k];
    }

#pragma omp simd aligned(t_1773, t_1774, t_1775, pc_x, pc_z, sli_1098, smh0_1038, smh0_1039, \
                         smh1_1038, smh1_1039, smi_1378, smi_1381, \
                         smi_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1773[k] = f_6 * smh0_1038[k]
                    - f_7 * smh1_1038[k]
                    + f_3 * pc_x[k] * smi_1381[k];

        t_1774[k] = f_8 * smh0_1039[k]
                    - f_9 * smh1_1039[k]
                    + f_3 * pc_x[k] * smi_1382[k];

        t_1775[k] = f_16 * sli_1098[k]
                    + f_3 * pc_z[k] * smi_1378[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, pc_x, pc_y, sli_1129, smh0_1041, smh0_1043, \
                         smh1_1041, smh1_1043, smi_1381, smi_1384, \
                         smi_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = f_8 * smh0_1041[k]
                    - f_9 * smh1_1041[k]
                    + f_3 * pc_x[k] * smi_1384[k];

        t_1777[k] = f_17 * sli_1129[k]
                    + f_3 * pc_y[k] * smi_1381[k];

        t_1778[k] = f_8 * smh0_1043[k]
                    - f_9 * smh1_1043[k]
                    + f_3 * pc_x[k] * smi_1386[k];
    }

#pragma omp simd aligned(t_1779, t_1780, t_1781, pc_x, pc_z, sli_1102, smh0_1044, smh0_1046, \
                         smh1_1044, smh1_1046, smi_1382, smi_1387, \
                         smi_1389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1779[k] = f_10 * smh0_1044[k]
                    - f_11 * smh1_1044[k]
                    + f_3 * pc_x[k] * smi_1387[k];

        t_1780[k] = f_16 * sli_1102[k]
                    + f_3 * pc_z[k] * smi_1382[k];

        t_1781[k] = f_10 * smh0_1046[k]
                    - f_11 * smh1_1046[k]
                    + f_3 * pc_x[k] * smi_1389[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, t_1785, pc_x, pc_y, sli_1134, smh0_1047, \
                         smh0_1049, smh1_1047, smh1_1049, smi_1386, smi_1390, smi_1392, \
                         smi_1393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = f_10 * smh0_1047[k]
                    - f_11 * smh1_1047[k]
                    + f_3 * pc_x[k] * smi_1390[k];

        t_1783[k] = f_17 * sli_1134[k]
                    + f_3 * pc_y[k] * smi_1386[k];

        t_1784[k] = f_10 * smh0_1049[k]
                    - f_11 * smh1_1049[k]
                    + f_3 * pc_x[k] * smi_1392[k];

        t_1785[k] = f_3 * pc_x[k] * smi_1393[k];
    }

#pragma omp simd aligned(t_1786, t_1787, t_1788, t_1789, t_1790, t_1791, pc_x, smi_1394, \
                         smi_1395, smi_1396, smi_1397, smi_1398, \
                         smi_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1786[k] = f_3 * pc_x[k] * smi_1394[k];

        t_1787[k] = f_3 * pc_x[k] * smi_1395[k];

        t_1788[k] = f_3 * pc_x[k] * smi_1396[k];

        t_1789[k] = f_3 * pc_x[k] * smi_1397[k];

        t_1790[k] = f_3 * pc_x[k] * smi_1398[k];

        t_1791[k] = f_3 * pc_x[k] * smi_1399[k];
    }

#pragma omp simd aligned(t_1792, t_1793, t_1794, pc_y, pc_z, sli_1113, sli_1141, sli_1143, \
                         smh0_1044, smh0_1046, smh1_1044, smh1_1046, smi_1393, \
                         smi_1395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1792[k] = f_17 * sli_1141[k]
                    + f_1 * smh0_1044[k]
                    - f_2 * smh1_1044[k]
                    + f_3 * pc_y[k] * smi_1393[k];

        t_1793[k] = f_16 * sli_1113[k]
                    + f_3 * pc_z[k] * smi_1393[k];

        t_1794[k] = f_17 * sli_1143[k]
                    + f_4 * smh0_1046[k]
                    - f_5 * smh1_1046[k]
                    + f_3 * pc_y[k] * smi_1395[k];
    }

#pragma omp simd aligned(t_1795, t_1796, t_1797, pc_y, sli_1144, sli_1145, sli_1146, \
                         smh0_1047, smh0_1048, smh0_1049, smh1_1047, smh1_1048, smh1_1049, \
                         smi_1396, smi_1397, smi_1398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1795[k] = f_17 * sli_1144[k]
                    + f_6 * smh0_1047[k]
                    - f_7 * smh1_1047[k]
                    + f_3 * pc_y[k] * smi_1396[k];

        t_1796[k] = f_17 * sli_1145[k]
                    + f_8 * smh0_1048[k]
                    - f_9 * smh1_1048[k]
                    + f_3 * pc_y[k] * smi_1397[k];

        t_1797[k] = f_17 * sli_1146[k]
                    + f_10 * smh0_1049[k]
                    - f_11 * smh1_1049[k]
                    + f_3 * pc_y[k] * smi_1398[k];
    }

#pragma omp simd aligned(t_1798, t_1799, t_1800, t_1801, pc_x, pc_y, pc_z, sli_1119, sli_1147, \
                         sli_1148, smh0_1049, smh0_1050, smh1_1049, smh1_1050, smi_1399, \
                         smi_1400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1798[k] = f_17 * sli_1147[k]
                    + f_3 * pc_y[k] * smi_1399[k];

        t_1799[k] = f_16 * sli_1119[k]
                    + f_1 * smh0_1049[k]
                    - f_2 * smh1_1049[k]
                    + f_3 * pc_z[k] * smi_1399[k];

        t_1800[k] = f_1 * smh0_1050[k]
                    - f_2 * smh1_1050[k]
                    + f_3 * pc_x[k] * smi_1400[k];

        t_1801[k] = f_16 * sli_1148[k]
                    + f_3 * pc_y[k] * smi_1400[k];
    }

#pragma omp simd aligned(t_1802, t_1803, t_1804, pc_x, pc_y, pc_z, sli_1120, sli_1150, \
                         smh0_1053, smh1_1053, smi_1400, smi_1402, \
                         smi_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1802[k] = f_17 * sli_1120[k]
                    + f_3 * pc_z[k] * smi_1400[k];

        t_1803[k] = f_4 * smh0_1053[k]
                    - f_5 * smh1_1053[k]
                    + f_3 * pc_x[k] * smi_1403[k];

        t_1804[k] = f_16 * sli_1150[k]
                    + f_3 * pc_y[k] * smi_1402[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, t_1808, pc_x, pc_y, pc_z, sli_1123, sli_1153, \
                         smh0_1055, smh0_1056, smh1_1055, smh1_1056, smi_1403, smi_1405, \
                         smi_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = f_4 * smh0_1055[k]
                    - f_5 * smh1_1055[k]
                    + f_3 * pc_x[k] * smi_1405[k];

        t_1806[k] = f_6 * smh0_1056[k]
                    - f_7 * smh1_1056[k]
                    + f_3 * pc_x[k] * smi_1406[k];

        t_1807[k] = f_17 * sli_1123[k]
                    + f_3 * pc_z[k] * smi_1403[k];

        t_1808[k] = f_16 * sli_1153[k]
                    + f_3 * pc_y[k] * smi_1405[k];
    }

#pragma omp simd aligned(t_1809, t_1810, t_1811, pc_x, pc_z, sli_1126, smh0_1059, smh0_1060, \
                         smh1_1059, smh1_1060, smi_1406, smi_1409, \
                         smi_1410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1809[k] = f_6 * smh0_1059[k]
                    - f_7 * smh1_1059[k]
                    + f_3 * pc_x[k] * smi_1409[k];

        t_1810[k] = f_8 * smh0_1060[k]
                    - f_9 * smh1_1060[k]
                    + f_3 * pc_x[k] * smi_1410[k];

        t_1811[k] = f_17 * sli_1126[k]
                    + f_3 * pc_z[k] * smi_1406[k];
    }

#pragma omp simd aligned(t_1812, t_1813, t_1814, pc_x, pc_y, sli_1157, smh0_1062, smh0_1064, \
                         smh1_1062, smh1_1064, smi_1409, smi_1412, \
                         smi_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1812[k] = f_8 * smh0_1062[k]
                    - f_9 * smh1_1062[k]
                    + f_3 * pc_x[k] * smi_1412[k];

        t_1813[k] = f_16 * sli_1157[k]
                    + f_3 * pc_y[k] * smi_1409[k];

        t_1814[k] = f_8 * smh0_1064[k]
                    - f_9 * smh1_1064[k]
                    + f_3 * pc_x[k] * smi_1414[k];
    }

#pragma omp simd aligned(t_1815, t_1816, t_1817, pc_x, pc_z, sli_1130, smh0_1065, smh0_1067, \
                         smh1_1065, smh1_1067, smi_1410, smi_1415, \
                         smi_1417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1815[k] = f_10 * smh0_1065[k]
                    - f_11 * smh1_1065[k]
                    + f_3 * pc_x[k] * smi_1415[k];

        t_1816[k] = f_17 * sli_1130[k]
                    + f_3 * pc_z[k] * smi_1410[k];

        t_1817[k] = f_10 * smh0_1067[k]
                    - f_11 * smh1_1067[k]
                    + f_3 * pc_x[k] * smi_1417[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t slk0,
                                                           const size_t sli, const size_t slk1,
                                                           const size_t smh0, const size_t smh1,
                                                           const size_t smi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

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
    auto *t_1848 = buffer.data(target + 1848);
    auto *t_1849 = buffer.data(target + 1849);
    auto *t_1850 = buffer.data(target + 1850);
    auto *t_1851 = buffer.data(target + 1851);
    auto *t_1852 = buffer.data(target + 1852);
    auto *t_1853 = buffer.data(target + 1853);
    auto *t_1854 = buffer.data(target + 1854);
    auto *t_1855 = buffer.data(target + 1855);
    auto *t_1856 = buffer.data(target + 1856);
    auto *t_1857 = buffer.data(target + 1857);
    auto *t_1858 = buffer.data(target + 1858);
    auto *t_1859 = buffer.data(target + 1859);
    auto *t_1860 = buffer.data(target + 1860);
    auto *t_1861 = buffer.data(target + 1861);
    auto *t_1862 = buffer.data(target + 1862);
    auto *t_1863 = buffer.data(target + 1863);
    auto *t_1864 = buffer.data(target + 1864);
    auto *t_1865 = buffer.data(target + 1865);
    auto *t_1866 = buffer.data(target + 1866);
    auto *t_1867 = buffer.data(target + 1867);
    auto *t_1868 = buffer.data(target + 1868);
    auto *t_1869 = buffer.data(target + 1869);
    auto *t_1870 = buffer.data(target + 1870);
    auto *t_1871 = buffer.data(target + 1871);
    auto *t_1872 = buffer.data(target + 1872);
    auto *t_1873 = buffer.data(target + 1873);
    auto *t_1874 = buffer.data(target + 1874);
    auto *t_1875 = buffer.data(target + 1875);
    auto *t_1876 = buffer.data(target + 1876);
    auto *t_1877 = buffer.data(target + 1877);
    auto *t_1878 = buffer.data(target + 1878);
    auto *t_1879 = buffer.data(target + 1879);
    auto *t_1880 = buffer.data(target + 1880);
    auto *t_1881 = buffer.data(target + 1881);
    auto *t_1882 = buffer.data(target + 1882);
    auto *t_1883 = buffer.data(target + 1883);
    auto *t_1884 = buffer.data(target + 1884);
    auto *t_1885 = buffer.data(target + 1885);
    auto *t_1886 = buffer.data(target + 1886);
    auto *t_1887 = buffer.data(target + 1887);
    auto *t_1888 = buffer.data(target + 1888);
    auto *t_1889 = buffer.data(target + 1889);
    auto *t_1890 = buffer.data(target + 1890);
    auto *t_1891 = buffer.data(target + 1891);
    auto *t_1892 = buffer.data(target + 1892);
    auto *t_1893 = buffer.data(target + 1893);
    auto *t_1894 = buffer.data(target + 1894);
    auto *t_1895 = buffer.data(target + 1895);
    auto *t_1896 = buffer.data(target + 1896);
    auto *t_1897 = buffer.data(target + 1897);
    auto *t_1898 = buffer.data(target + 1898);
    auto *t_1899 = buffer.data(target + 1899);
    auto *t_1900 = buffer.data(target + 1900);
    auto *t_1901 = buffer.data(target + 1901);
    auto *t_1902 = buffer.data(target + 1902);
    auto *t_1903 = buffer.data(target + 1903);
    auto *t_1904 = buffer.data(target + 1904);
    auto *t_1905 = buffer.data(target + 1905);
    auto *t_1906 = buffer.data(target + 1906);
    auto *t_1907 = buffer.data(target + 1907);
    auto *t_1908 = buffer.data(target + 1908);
    auto *t_1909 = buffer.data(target + 1909);
    auto *t_1910 = buffer.data(target + 1910);
    auto *t_1911 = buffer.data(target + 1911);
    auto *t_1912 = buffer.data(target + 1912);
    auto *t_1913 = buffer.data(target + 1913);
    auto *t_1914 = buffer.data(target + 1914);
    auto *t_1915 = buffer.data(target + 1915);
    auto *t_1916 = buffer.data(target + 1916);
    auto *t_1917 = buffer.data(target + 1917);
    auto *t_1918 = buffer.data(target + 1918);
    auto *t_1919 = buffer.data(target + 1919);
    auto *t_1920 = buffer.data(target + 1920);
    auto *t_1921 = buffer.data(target + 1921);
    auto *t_1922 = buffer.data(target + 1922);
    auto *t_1923 = buffer.data(target + 1923);
    auto *t_1924 = buffer.data(target + 1924);
    auto *t_1925 = buffer.data(target + 1925);
    auto *t_1926 = buffer.data(target + 1926);
    auto *t_1927 = buffer.data(target + 1927);
    auto *t_1928 = buffer.data(target + 1928);
    auto *t_1929 = buffer.data(target + 1929);
    auto *t_1930 = buffer.data(target + 1930);
    auto *t_1931 = buffer.data(target + 1931);
    auto *t_1932 = buffer.data(target + 1932);
    auto *t_1933 = buffer.data(target + 1933);
    auto *t_1934 = buffer.data(target + 1934);
    auto *t_1935 = buffer.data(target + 1935);
    auto *t_1936 = buffer.data(target + 1936);
    auto *t_1937 = buffer.data(target + 1937);
    auto *t_1938 = buffer.data(target + 1938);
    auto *t_1939 = buffer.data(target + 1939);
    auto *t_1940 = buffer.data(target + 1940);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_1584 = buffer.data(slk0 + 1584);
    const auto *slk0_1589 = buffer.data(slk0 + 1589);
    const auto *slk0_1593 = buffer.data(slk0 + 1593);
    const auto *slk0_1598 = buffer.data(slk0 + 1598);
    const auto *slk0_1604 = buffer.data(slk0 + 1604);
    const auto *slk0_1612 = buffer.data(slk0 + 1612);
    const auto *slk0_1614 = buffer.data(slk0 + 1614);
    const auto *slk0_1615 = buffer.data(slk0 + 1615);
    const auto *slk0_1616 = buffer.data(slk0 + 1616);

    const auto *sli_1141 = buffer.data(sli + 1141);
    const auto *sli_1147 = buffer.data(sli + 1147);
    const auto *sli_1148 = buffer.data(sli + 1148);
    const auto *sli_1151 = buffer.data(sli + 1151);
    const auto *sli_1154 = buffer.data(sli + 1154);
    const auto *sli_1158 = buffer.data(sli + 1158);
    const auto *sli_1162 = buffer.data(sli + 1162);
    const auto *sli_1169 = buffer.data(sli + 1169);
    const auto *sli_1171 = buffer.data(sli + 1171);
    const auto *sli_1172 = buffer.data(sli + 1172);
    const auto *sli_1173 = buffer.data(sli + 1173);
    const auto *sli_1174 = buffer.data(sli + 1174);
    const auto *sli_1175 = buffer.data(sli + 1175);
    const auto *sli_1176 = buffer.data(sli + 1176);
    const auto *sli_1178 = buffer.data(sli + 1178);
    const auto *sli_1179 = buffer.data(sli + 1179);
    const auto *sli_1181 = buffer.data(sli + 1181);
    const auto *sli_1182 = buffer.data(sli + 1182);
    const auto *sli_1185 = buffer.data(sli + 1185);
    const auto *sli_1186 = buffer.data(sli + 1186);
    const auto *sli_1190 = buffer.data(sli + 1190);
    const auto *sli_1197 = buffer.data(sli + 1197);
    const auto *sli_1199 = buffer.data(sli + 1199);
    const auto *sli_1200 = buffer.data(sli + 1200);
    const auto *sli_1201 = buffer.data(sli + 1201);
    const auto *sli_1202 = buffer.data(sli + 1202);
    const auto *sli_1203 = buffer.data(sli + 1203);
    const auto *sli_1204 = buffer.data(sli + 1204);
    const auto *sli_1206 = buffer.data(sli + 1206);
    const auto *sli_1207 = buffer.data(sli + 1207);
    const auto *sli_1209 = buffer.data(sli + 1209);
    const auto *sli_1210 = buffer.data(sli + 1210);
    const auto *sli_1213 = buffer.data(sli + 1213);
    const auto *sli_1214 = buffer.data(sli + 1214);
    const auto *sli_1218 = buffer.data(sli + 1218);
    const auto *sli_1225 = buffer.data(sli + 1225);
    const auto *sli_1227 = buffer.data(sli + 1227);
    const auto *sli_1228 = buffer.data(sli + 1228);
    const auto *sli_1229 = buffer.data(sli + 1229);
    const auto *sli_1230 = buffer.data(sli + 1230);
    const auto *sli_1231 = buffer.data(sli + 1231);
    const auto *sli_1232 = buffer.data(sli + 1232);
    const auto *sli_1234 = buffer.data(sli + 1234);
    const auto *sli_1237 = buffer.data(sli + 1237);
    const auto *sli_1241 = buffer.data(sli + 1241);
    const auto *sli_1246 = buffer.data(sli + 1246);
    const auto *sli_1253 = buffer.data(sli + 1253);
    const auto *sli_1255 = buffer.data(sli + 1255);
    const auto *sli_1256 = buffer.data(sli + 1256);
    const auto *sli_1257 = buffer.data(sli + 1257);

    const auto *slk1_1584 = buffer.data(slk1 + 1584);
    const auto *slk1_1589 = buffer.data(slk1 + 1589);
    const auto *slk1_1593 = buffer.data(slk1 + 1593);
    const auto *slk1_1598 = buffer.data(slk1 + 1598);
    const auto *slk1_1604 = buffer.data(slk1 + 1604);
    const auto *slk1_1612 = buffer.data(slk1 + 1612);
    const auto *slk1_1614 = buffer.data(slk1 + 1614);
    const auto *slk1_1615 = buffer.data(slk1 + 1615);
    const auto *slk1_1616 = buffer.data(slk1 + 1616);

    const auto *smh0_1065 = buffer.data(smh0 + 1065);
    const auto *smh0_1067 = buffer.data(smh0 + 1067);
    const auto *smh0_1068 = buffer.data(smh0 + 1068);
    const auto *smh0_1069 = buffer.data(smh0 + 1069);
    const auto *smh0_1070 = buffer.data(smh0 + 1070);
    const auto *smh0_1071 = buffer.data(smh0 + 1071);
    const auto *smh0_1074 = buffer.data(smh0 + 1074);
    const auto *smh0_1076 = buffer.data(smh0 + 1076);
    const auto *smh0_1077 = buffer.data(smh0 + 1077);
    const auto *smh0_1080 = buffer.data(smh0 + 1080);
    const auto *smh0_1081 = buffer.data(smh0 + 1081);
    const auto *smh0_1083 = buffer.data(smh0 + 1083);
    const auto *smh0_1085 = buffer.data(smh0 + 1085);
    const auto *smh0_1086 = buffer.data(smh0 + 1086);
    const auto *smh0_1088 = buffer.data(smh0 + 1088);
    const auto *smh0_1089 = buffer.data(smh0 + 1089);
    const auto *smh0_1090 = buffer.data(smh0 + 1090);
    const auto *smh0_1091 = buffer.data(smh0 + 1091);
    const auto *smh0_1092 = buffer.data(smh0 + 1092);
    const auto *smh0_1095 = buffer.data(smh0 + 1095);
    const auto *smh0_1097 = buffer.data(smh0 + 1097);
    const auto *smh0_1098 = buffer.data(smh0 + 1098);
    const auto *smh0_1101 = buffer.data(smh0 + 1101);
    const auto *smh0_1102 = buffer.data(smh0 + 1102);
    const auto *smh0_1104 = buffer.data(smh0 + 1104);
    const auto *smh0_1106 = buffer.data(smh0 + 1106);
    const auto *smh0_1107 = buffer.data(smh0 + 1107);
    const auto *smh0_1109 = buffer.data(smh0 + 1109);
    const auto *smh0_1110 = buffer.data(smh0 + 1110);
    const auto *smh0_1111 = buffer.data(smh0 + 1111);
    const auto *smh0_1112 = buffer.data(smh0 + 1112);
    const auto *smh0_1116 = buffer.data(smh0 + 1116);
    const auto *smh0_1119 = buffer.data(smh0 + 1119);
    const auto *smh0_1123 = buffer.data(smh0 + 1123);
    const auto *smh0_1125 = buffer.data(smh0 + 1125);
    const auto *smh0_1128 = buffer.data(smh0 + 1128);
    const auto *smh0_1130 = buffer.data(smh0 + 1130);
    const auto *smh0_1131 = buffer.data(smh0 + 1131);

    const auto *smh1_1065 = buffer.data(smh1 + 1065);
    const auto *smh1_1067 = buffer.data(smh1 + 1067);
    const auto *smh1_1068 = buffer.data(smh1 + 1068);
    const auto *smh1_1069 = buffer.data(smh1 + 1069);
    const auto *smh1_1070 = buffer.data(smh1 + 1070);
    const auto *smh1_1071 = buffer.data(smh1 + 1071);
    const auto *smh1_1074 = buffer.data(smh1 + 1074);
    const auto *smh1_1076 = buffer.data(smh1 + 1076);
    const auto *smh1_1077 = buffer.data(smh1 + 1077);
    const auto *smh1_1080 = buffer.data(smh1 + 1080);
    const auto *smh1_1081 = buffer.data(smh1 + 1081);
    const auto *smh1_1083 = buffer.data(smh1 + 1083);
    const auto *smh1_1085 = buffer.data(smh1 + 1085);
    const auto *smh1_1086 = buffer.data(smh1 + 1086);
    const auto *smh1_1088 = buffer.data(smh1 + 1088);
    const auto *smh1_1089 = buffer.data(smh1 + 1089);
    const auto *smh1_1090 = buffer.data(smh1 + 1090);
    const auto *smh1_1091 = buffer.data(smh1 + 1091);
    const auto *smh1_1092 = buffer.data(smh1 + 1092);
    const auto *smh1_1095 = buffer.data(smh1 + 1095);
    const auto *smh1_1097 = buffer.data(smh1 + 1097);
    const auto *smh1_1098 = buffer.data(smh1 + 1098);
    const auto *smh1_1101 = buffer.data(smh1 + 1101);
    const auto *smh1_1102 = buffer.data(smh1 + 1102);
    const auto *smh1_1104 = buffer.data(smh1 + 1104);
    const auto *smh1_1106 = buffer.data(smh1 + 1106);
    const auto *smh1_1107 = buffer.data(smh1 + 1107);
    const auto *smh1_1109 = buffer.data(smh1 + 1109);
    const auto *smh1_1110 = buffer.data(smh1 + 1110);
    const auto *smh1_1111 = buffer.data(smh1 + 1111);
    const auto *smh1_1112 = buffer.data(smh1 + 1112);
    const auto *smh1_1116 = buffer.data(smh1 + 1116);
    const auto *smh1_1119 = buffer.data(smh1 + 1119);
    const auto *smh1_1123 = buffer.data(smh1 + 1123);
    const auto *smh1_1125 = buffer.data(smh1 + 1125);
    const auto *smh1_1128 = buffer.data(smh1 + 1128);
    const auto *smh1_1130 = buffer.data(smh1 + 1130);
    const auto *smh1_1131 = buffer.data(smh1 + 1131);

    const auto *smi_1414 = buffer.data(smi + 1414);
    const auto *smi_1418 = buffer.data(smi + 1418);
    const auto *smi_1420 = buffer.data(smi + 1420);
    const auto *smi_1421 = buffer.data(smi + 1421);
    const auto *smi_1422 = buffer.data(smi + 1422);
    const auto *smi_1423 = buffer.data(smi + 1423);
    const auto *smi_1424 = buffer.data(smi + 1424);
    const auto *smi_1425 = buffer.data(smi + 1425);
    const auto *smi_1426 = buffer.data(smi + 1426);
    const auto *smi_1427 = buffer.data(smi + 1427);
    const auto *smi_1428 = buffer.data(smi + 1428);
    const auto *smi_1430 = buffer.data(smi + 1430);
    const auto *smi_1431 = buffer.data(smi + 1431);
    const auto *smi_1433 = buffer.data(smi + 1433);
    const auto *smi_1434 = buffer.data(smi + 1434);
    const auto *smi_1437 = buffer.data(smi + 1437);
    const auto *smi_1438 = buffer.data(smi + 1438);
    const auto *smi_1440 = buffer.data(smi + 1440);
    const auto *smi_1442 = buffer.data(smi + 1442);
    const auto *smi_1443 = buffer.data(smi + 1443);
    const auto *smi_1445 = buffer.data(smi + 1445);
    const auto *smi_1446 = buffer.data(smi + 1446);
    const auto *smi_1448 = buffer.data(smi + 1448);
    const auto *smi_1449 = buffer.data(smi + 1449);
    const auto *smi_1450 = buffer.data(smi + 1450);
    const auto *smi_1451 = buffer.data(smi + 1451);
    const auto *smi_1452 = buffer.data(smi + 1452);
    const auto *smi_1453 = buffer.data(smi + 1453);
    const auto *smi_1454 = buffer.data(smi + 1454);
    const auto *smi_1455 = buffer.data(smi + 1455);
    const auto *smi_1456 = buffer.data(smi + 1456);
    const auto *smi_1458 = buffer.data(smi + 1458);
    const auto *smi_1459 = buffer.data(smi + 1459);
    const auto *smi_1461 = buffer.data(smi + 1461);
    const auto *smi_1462 = buffer.data(smi + 1462);
    const auto *smi_1465 = buffer.data(smi + 1465);
    const auto *smi_1466 = buffer.data(smi + 1466);
    const auto *smi_1468 = buffer.data(smi + 1468);
    const auto *smi_1470 = buffer.data(smi + 1470);
    const auto *smi_1471 = buffer.data(smi + 1471);
    const auto *smi_1473 = buffer.data(smi + 1473);
    const auto *smi_1474 = buffer.data(smi + 1474);
    const auto *smi_1476 = buffer.data(smi + 1476);
    const auto *smi_1477 = buffer.data(smi + 1477);
    const auto *smi_1478 = buffer.data(smi + 1478);
    const auto *smi_1479 = buffer.data(smi + 1479);
    const auto *smi_1480 = buffer.data(smi + 1480);
    const auto *smi_1481 = buffer.data(smi + 1481);
    const auto *smi_1482 = buffer.data(smi + 1482);
    const auto *smi_1483 = buffer.data(smi + 1483);
    const auto *smi_1484 = buffer.data(smi + 1484);
    const auto *smi_1486 = buffer.data(smi + 1486);
    const auto *smi_1487 = buffer.data(smi + 1487);
    const auto *smi_1489 = buffer.data(smi + 1489);
    const auto *smi_1490 = buffer.data(smi + 1490);
    const auto *smi_1493 = buffer.data(smi + 1493);
    const auto *smi_1494 = buffer.data(smi + 1494);
    const auto *smi_1496 = buffer.data(smi + 1496);
    const auto *smi_1498 = buffer.data(smi + 1498);
    const auto *smi_1499 = buffer.data(smi + 1499);
    const auto *smi_1501 = buffer.data(smi + 1501);
    const auto *smi_1502 = buffer.data(smi + 1502);
    const auto *smi_1505 = buffer.data(smi + 1505);
    const auto *smi_1506 = buffer.data(smi + 1506);
    const auto *smi_1507 = buffer.data(smi + 1507);
    const auto *smi_1508 = buffer.data(smi + 1508);
    const auto *smi_1509 = buffer.data(smi + 1509);
    const auto *smi_1510 = buffer.data(smi + 1510);
    const auto *smi_1511 = buffer.data(smi + 1511);

#pragma omp simd aligned(t_1818, t_1819, t_1820, t_1821, pc_x, pc_y, sli_1162, smh0_1068, \
                         smh0_1070, smh1_1068, smh1_1070, smi_1414, smi_1418, smi_1420, \
                         smi_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1818[k] = f_10 * smh0_1068[k]
                    - f_11 * smh1_1068[k]
                    + f_3 * pc_x[k] * smi_1418[k];

        t_1819[k] = f_16 * sli_1162[k]
                    + f_3 * pc_y[k] * smi_1414[k];

        t_1820[k] = f_10 * smh0_1070[k]
                    - f_11 * smh1_1070[k]
                    + f_3 * pc_x[k] * smi_1420[k];

        t_1821[k] = f_3 * pc_x[k] * smi_1421[k];
    }

#pragma omp simd aligned(t_1822, t_1823, t_1824, t_1825, t_1826, t_1827, pc_x, smi_1422, \
                         smi_1423, smi_1424, smi_1425, smi_1426, \
                         smi_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1822[k] = f_3 * pc_x[k] * smi_1422[k];

        t_1823[k] = f_3 * pc_x[k] * smi_1423[k];

        t_1824[k] = f_3 * pc_x[k] * smi_1424[k];

        t_1825[k] = f_3 * pc_x[k] * smi_1425[k];

        t_1826[k] = f_3 * pc_x[k] * smi_1426[k];

        t_1827[k] = f_3 * pc_x[k] * smi_1427[k];
    }

#pragma omp simd aligned(t_1828, t_1829, t_1830, pc_y, pc_z, sli_1141, sli_1169, sli_1171, \
                         smh0_1065, smh0_1067, smh1_1065, smh1_1067, smi_1421, \
                         smi_1423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1828[k] = f_16 * sli_1169[k]
                    + f_1 * smh0_1065[k]
                    - f_2 * smh1_1065[k]
                    + f_3 * pc_y[k] * smi_1421[k];

        t_1829[k] = f_17 * sli_1141[k]
                    + f_3 * pc_z[k] * smi_1421[k];

        t_1830[k] = f_16 * sli_1171[k]
                    + f_4 * smh0_1067[k]
                    - f_5 * smh1_1067[k]
                    + f_3 * pc_y[k] * smi_1423[k];
    }

#pragma omp simd aligned(t_1831, t_1832, t_1833, pc_y, sli_1172, sli_1173, sli_1174, \
                         smh0_1068, smh0_1069, smh0_1070, smh1_1068, smh1_1069, smh1_1070, \
                         smi_1424, smi_1425, smi_1426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1831[k] = f_16 * sli_1172[k]
                    + f_6 * smh0_1068[k]
                    - f_7 * smh1_1068[k]
                    + f_3 * pc_y[k] * smi_1424[k];

        t_1832[k] = f_16 * sli_1173[k]
                    + f_8 * smh0_1069[k]
                    - f_9 * smh1_1069[k]
                    + f_3 * pc_y[k] * smi_1425[k];

        t_1833[k] = f_16 * sli_1174[k]
                    + f_10 * smh0_1070[k]
                    - f_11 * smh1_1070[k]
                    + f_3 * pc_y[k] * smi_1426[k];
    }

#pragma omp simd aligned(t_1834, t_1835, t_1836, t_1837, pc_x, pc_y, pc_z, sli_1147, sli_1175, \
                         sli_1176, smh0_1070, smh0_1071, smh1_1070, smh1_1071, smi_1427, \
                         smi_1428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1834[k] = f_16 * sli_1175[k]
                    + f_3 * pc_y[k] * smi_1427[k];

        t_1835[k] = f_17 * sli_1147[k]
                    + f_1 * smh0_1070[k]
                    - f_2 * smh1_1070[k]
                    + f_3 * pc_z[k] * smi_1427[k];

        t_1836[k] = f_1 * smh0_1071[k]
                    - f_2 * smh1_1071[k]
                    + f_3 * pc_x[k] * smi_1428[k];

        t_1837[k] = f_15 * sli_1176[k]
                    + f_3 * pc_y[k] * smi_1428[k];
    }

#pragma omp simd aligned(t_1838, t_1839, t_1840, pc_x, pc_y, pc_z, sli_1148, sli_1178, \
                         smh0_1074, smh1_1074, smi_1428, smi_1430, \
                         smi_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1838[k] = f_20 * sli_1148[k]
                    + f_3 * pc_z[k] * smi_1428[k];

        t_1839[k] = f_4 * smh0_1074[k]
                    - f_5 * smh1_1074[k]
                    + f_3 * pc_x[k] * smi_1431[k];

        t_1840[k] = f_15 * sli_1178[k]
                    + f_3 * pc_y[k] * smi_1430[k];
    }

#pragma omp simd aligned(t_1841, t_1842, t_1843, t_1844, pc_x, pc_y, pc_z, sli_1151, sli_1181, \
                         smh0_1076, smh0_1077, smh1_1076, smh1_1077, smi_1431, smi_1433, \
                         smi_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1841[k] = f_4 * smh0_1076[k]
                    - f_5 * smh1_1076[k]
                    + f_3 * pc_x[k] * smi_1433[k];

        t_1842[k] = f_6 * smh0_1077[k]
                    - f_7 * smh1_1077[k]
                    + f_3 * pc_x[k] * smi_1434[k];

        t_1843[k] = f_20 * sli_1151[k]
                    + f_3 * pc_z[k] * smi_1431[k];

        t_1844[k] = f_15 * sli_1181[k]
                    + f_3 * pc_y[k] * smi_1433[k];
    }

#pragma omp simd aligned(t_1845, t_1846, t_1847, pc_x, pc_z, sli_1154, smh0_1080, smh0_1081, \
                         smh1_1080, smh1_1081, smi_1434, smi_1437, \
                         smi_1438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1845[k] = f_6 * smh0_1080[k]
                    - f_7 * smh1_1080[k]
                    + f_3 * pc_x[k] * smi_1437[k];

        t_1846[k] = f_8 * smh0_1081[k]
                    - f_9 * smh1_1081[k]
                    + f_3 * pc_x[k] * smi_1438[k];

        t_1847[k] = f_20 * sli_1154[k]
                    + f_3 * pc_z[k] * smi_1434[k];
    }

#pragma omp simd aligned(t_1848, t_1849, t_1850, pc_x, pc_y, sli_1185, smh0_1083, smh0_1085, \
                         smh1_1083, smh1_1085, smi_1437, smi_1440, \
                         smi_1442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1848[k] = f_8 * smh0_1083[k]
                    - f_9 * smh1_1083[k]
                    + f_3 * pc_x[k] * smi_1440[k];

        t_1849[k] = f_15 * sli_1185[k]
                    + f_3 * pc_y[k] * smi_1437[k];

        t_1850[k] = f_8 * smh0_1085[k]
                    - f_9 * smh1_1085[k]
                    + f_3 * pc_x[k] * smi_1442[k];
    }

#pragma omp simd aligned(t_1851, t_1852, t_1853, pc_x, pc_z, sli_1158, smh0_1086, smh0_1088, \
                         smh1_1086, smh1_1088, smi_1438, smi_1443, \
                         smi_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1851[k] = f_10 * smh0_1086[k]
                    - f_11 * smh1_1086[k]
                    + f_3 * pc_x[k] * smi_1443[k];

        t_1852[k] = f_20 * sli_1158[k]
                    + f_3 * pc_z[k] * smi_1438[k];

        t_1853[k] = f_10 * smh0_1088[k]
                    - f_11 * smh1_1088[k]
                    + f_3 * pc_x[k] * smi_1445[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, t_1857, pc_x, pc_y, sli_1190, smh0_1089, \
                         smh0_1091, smh1_1089, smh1_1091, smi_1442, smi_1446, smi_1448, \
                         smi_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = f_10 * smh0_1089[k]
                    - f_11 * smh1_1089[k]
                    + f_3 * pc_x[k] * smi_1446[k];

        t_1855[k] = f_15 * sli_1190[k]
                    + f_3 * pc_y[k] * smi_1442[k];

        t_1856[k] = f_10 * smh0_1091[k]
                    - f_11 * smh1_1091[k]
                    + f_3 * pc_x[k] * smi_1448[k];

        t_1857[k] = f_3 * pc_x[k] * smi_1449[k];
    }

#pragma omp simd aligned(t_1858, t_1859, t_1860, t_1861, t_1862, t_1863, pc_x, smi_1450, \
                         smi_1451, smi_1452, smi_1453, smi_1454, \
                         smi_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1858[k] = f_3 * pc_x[k] * smi_1450[k];

        t_1859[k] = f_3 * pc_x[k] * smi_1451[k];

        t_1860[k] = f_3 * pc_x[k] * smi_1452[k];

        t_1861[k] = f_3 * pc_x[k] * smi_1453[k];

        t_1862[k] = f_3 * pc_x[k] * smi_1454[k];

        t_1863[k] = f_3 * pc_x[k] * smi_1455[k];
    }

#pragma omp simd aligned(t_1864, t_1865, t_1866, pc_y, pc_z, sli_1169, sli_1197, sli_1199, \
                         smh0_1086, smh0_1088, smh1_1086, smh1_1088, smi_1449, \
                         smi_1451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1864[k] = f_15 * sli_1197[k]
                    + f_1 * smh0_1086[k]
                    - f_2 * smh1_1086[k]
                    + f_3 * pc_y[k] * smi_1449[k];

        t_1865[k] = f_20 * sli_1169[k]
                    + f_3 * pc_z[k] * smi_1449[k];

        t_1866[k] = f_15 * sli_1199[k]
                    + f_4 * smh0_1088[k]
                    - f_5 * smh1_1088[k]
                    + f_3 * pc_y[k] * smi_1451[k];
    }

#pragma omp simd aligned(t_1867, t_1868, t_1869, pc_y, sli_1200, sli_1201, sli_1202, \
                         smh0_1089, smh0_1090, smh0_1091, smh1_1089, smh1_1090, smh1_1091, \
                         smi_1452, smi_1453, smi_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1867[k] = f_15 * sli_1200[k]
                    + f_6 * smh0_1089[k]
                    - f_7 * smh1_1089[k]
                    + f_3 * pc_y[k] * smi_1452[k];

        t_1868[k] = f_15 * sli_1201[k]
                    + f_8 * smh0_1090[k]
                    - f_9 * smh1_1090[k]
                    + f_3 * pc_y[k] * smi_1453[k];

        t_1869[k] = f_15 * sli_1202[k]
                    + f_10 * smh0_1091[k]
                    - f_11 * smh1_1091[k]
                    + f_3 * pc_y[k] * smi_1454[k];
    }

#pragma omp simd aligned(t_1870, t_1871, t_1872, t_1873, pc_x, pc_y, pc_z, sli_1175, sli_1203, \
                         sli_1204, smh0_1091, smh0_1092, smh1_1091, smh1_1092, smi_1455, \
                         smi_1456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1870[k] = f_15 * sli_1203[k]
                    + f_3 * pc_y[k] * smi_1455[k];

        t_1871[k] = f_20 * sli_1175[k]
                    + f_1 * smh0_1091[k]
                    - f_2 * smh1_1091[k]
                    + f_3 * pc_z[k] * smi_1455[k];

        t_1872[k] = f_1 * smh0_1092[k]
                    - f_2 * smh1_1092[k]
                    + f_3 * pc_x[k] * smi_1456[k];

        t_1873[k] = f_14 * sli_1204[k]
                    + f_3 * pc_y[k] * smi_1456[k];
    }

#pragma omp simd aligned(t_1874, t_1875, t_1876, pc_x, pc_y, pc_z, sli_1176, sli_1206, \
                         smh0_1095, smh1_1095, smi_1456, smi_1458, \
                         smi_1459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1874[k] = f_19 * sli_1176[k]
                    + f_3 * pc_z[k] * smi_1456[k];

        t_1875[k] = f_4 * smh0_1095[k]
                    - f_5 * smh1_1095[k]
                    + f_3 * pc_x[k] * smi_1459[k];

        t_1876[k] = f_14 * sli_1206[k]
                    + f_3 * pc_y[k] * smi_1458[k];
    }

#pragma omp simd aligned(t_1877, t_1878, t_1879, t_1880, pc_x, pc_y, pc_z, sli_1179, sli_1209, \
                         smh0_1097, smh0_1098, smh1_1097, smh1_1098, smi_1459, smi_1461, \
                         smi_1462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1877[k] = f_4 * smh0_1097[k]
                    - f_5 * smh1_1097[k]
                    + f_3 * pc_x[k] * smi_1461[k];

        t_1878[k] = f_6 * smh0_1098[k]
                    - f_7 * smh1_1098[k]
                    + f_3 * pc_x[k] * smi_1462[k];

        t_1879[k] = f_19 * sli_1179[k]
                    + f_3 * pc_z[k] * smi_1459[k];

        t_1880[k] = f_14 * sli_1209[k]
                    + f_3 * pc_y[k] * smi_1461[k];
    }

#pragma omp simd aligned(t_1881, t_1882, t_1883, pc_x, pc_z, sli_1182, smh0_1101, smh0_1102, \
                         smh1_1101, smh1_1102, smi_1462, smi_1465, \
                         smi_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1881[k] = f_6 * smh0_1101[k]
                    - f_7 * smh1_1101[k]
                    + f_3 * pc_x[k] * smi_1465[k];

        t_1882[k] = f_8 * smh0_1102[k]
                    - f_9 * smh1_1102[k]
                    + f_3 * pc_x[k] * smi_1466[k];

        t_1883[k] = f_19 * sli_1182[k]
                    + f_3 * pc_z[k] * smi_1462[k];
    }

#pragma omp simd aligned(t_1884, t_1885, t_1886, pc_x, pc_y, sli_1213, smh0_1104, smh0_1106, \
                         smh1_1104, smh1_1106, smi_1465, smi_1468, \
                         smi_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1884[k] = f_8 * smh0_1104[k]
                    - f_9 * smh1_1104[k]
                    + f_3 * pc_x[k] * smi_1468[k];

        t_1885[k] = f_14 * sli_1213[k]
                    + f_3 * pc_y[k] * smi_1465[k];

        t_1886[k] = f_8 * smh0_1106[k]
                    - f_9 * smh1_1106[k]
                    + f_3 * pc_x[k] * smi_1470[k];
    }

#pragma omp simd aligned(t_1887, t_1888, t_1889, pc_x, pc_z, sli_1186, smh0_1107, smh0_1109, \
                         smh1_1107, smh1_1109, smi_1466, smi_1471, \
                         smi_1473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1887[k] = f_10 * smh0_1107[k]
                    - f_11 * smh1_1107[k]
                    + f_3 * pc_x[k] * smi_1471[k];

        t_1888[k] = f_19 * sli_1186[k]
                    + f_3 * pc_z[k] * smi_1466[k];

        t_1889[k] = f_10 * smh0_1109[k]
                    - f_11 * smh1_1109[k]
                    + f_3 * pc_x[k] * smi_1473[k];
    }

#pragma omp simd aligned(t_1890, t_1891, t_1892, t_1893, pc_x, pc_y, sli_1218, smh0_1110, \
                         smh0_1112, smh1_1110, smh1_1112, smi_1470, smi_1474, smi_1476, \
                         smi_1477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1890[k] = f_10 * smh0_1110[k]
                    - f_11 * smh1_1110[k]
                    + f_3 * pc_x[k] * smi_1474[k];

        t_1891[k] = f_14 * sli_1218[k]
                    + f_3 * pc_y[k] * smi_1470[k];

        t_1892[k] = f_10 * smh0_1112[k]
                    - f_11 * smh1_1112[k]
                    + f_3 * pc_x[k] * smi_1476[k];

        t_1893[k] = f_3 * pc_x[k] * smi_1477[k];
    }

#pragma omp simd aligned(t_1894, t_1895, t_1896, t_1897, t_1898, t_1899, pc_x, smi_1478, \
                         smi_1479, smi_1480, smi_1481, smi_1482, \
                         smi_1483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1894[k] = f_3 * pc_x[k] * smi_1478[k];

        t_1895[k] = f_3 * pc_x[k] * smi_1479[k];

        t_1896[k] = f_3 * pc_x[k] * smi_1480[k];

        t_1897[k] = f_3 * pc_x[k] * smi_1481[k];

        t_1898[k] = f_3 * pc_x[k] * smi_1482[k];

        t_1899[k] = f_3 * pc_x[k] * smi_1483[k];
    }

#pragma omp simd aligned(t_1900, t_1901, t_1902, pc_y, pc_z, sli_1197, sli_1225, sli_1227, \
                         smh0_1107, smh0_1109, smh1_1107, smh1_1109, smi_1477, \
                         smi_1479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1900[k] = f_14 * sli_1225[k]
                    + f_1 * smh0_1107[k]
                    - f_2 * smh1_1107[k]
                    + f_3 * pc_y[k] * smi_1477[k];

        t_1901[k] = f_19 * sli_1197[k]
                    + f_3 * pc_z[k] * smi_1477[k];

        t_1902[k] = f_14 * sli_1227[k]
                    + f_4 * smh0_1109[k]
                    - f_5 * smh1_1109[k]
                    + f_3 * pc_y[k] * smi_1479[k];
    }

#pragma omp simd aligned(t_1903, t_1904, t_1905, pc_y, sli_1228, sli_1229, sli_1230, \
                         smh0_1110, smh0_1111, smh0_1112, smh1_1110, smh1_1111, smh1_1112, \
                         smi_1480, smi_1481, smi_1482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1903[k] = f_14 * sli_1228[k]
                    + f_6 * smh0_1110[k]
                    - f_7 * smh1_1110[k]
                    + f_3 * pc_y[k] * smi_1480[k];

        t_1904[k] = f_14 * sli_1229[k]
                    + f_8 * smh0_1111[k]
                    - f_9 * smh1_1111[k]
                    + f_3 * pc_y[k] * smi_1481[k];

        t_1905[k] = f_14 * sli_1230[k]
                    + f_10 * smh0_1112[k]
                    - f_11 * smh1_1112[k]
                    + f_3 * pc_y[k] * smi_1482[k];
    }

#pragma omp simd aligned(t_1906, t_1907, t_1908, t_1909, pb_y, pc_y, pc_z, slk0_1584, \
                         sli_1203, sli_1231, sli_1232, slk1_1584, smh0_1112, smh1_1112, \
                         smi_1483, smi_1484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1906[k] = f_14 * sli_1231[k]
                    + f_3 * pc_y[k] * smi_1483[k];

        t_1907[k] = f_19 * sli_1203[k]
                    + f_1 * smh0_1112[k]
                    - f_2 * smh1_1112[k]
                    + f_3 * pc_z[k] * smi_1483[k];

        t_1908[k] = pb_y[k] * slk0_1584[k]
                    - f_12 * pc_y[k] * slk1_1584[k];

        t_1909[k] = f_13 * sli_1232[k]
                    + f_3 * pc_y[k] * smi_1484[k];
    }

#pragma omp simd aligned(t_1910, t_1911, t_1912, pc_x, pc_y, pc_z, sli_1204, sli_1234, \
                         smh0_1116, smh1_1116, smi_1484, smi_1486, \
                         smi_1487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1910[k] = f_18 * sli_1204[k]
                    + f_3 * pc_z[k] * smi_1484[k];

        t_1911[k] = f_4 * smh0_1116[k]
                    - f_5 * smh1_1116[k]
                    + f_3 * pc_x[k] * smi_1487[k];

        t_1912[k] = f_13 * sli_1234[k]
                    + f_3 * pc_y[k] * smi_1486[k];
    }

#pragma omp simd aligned(t_1913, t_1914, t_1915, pb_y, pc_x, pc_y, pc_z, slk0_1589, sli_1207, \
                         slk1_1589, smh0_1119, smh1_1119, smi_1487, \
                         smi_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1913[k] = pb_y[k] * slk0_1589[k]
                    - f_12 * pc_y[k] * slk1_1589[k];

        t_1914[k] = f_6 * smh0_1119[k]
                    - f_7 * smh1_1119[k]
                    + f_3 * pc_x[k] * smi_1490[k];

        t_1915[k] = f_18 * sli_1207[k]
                    + f_3 * pc_z[k] * smi_1487[k];
    }

#pragma omp simd aligned(t_1916, t_1917, t_1918, pb_y, pc_x, pc_y, slk0_1593, sli_1237, \
                         slk1_1593, smh0_1123, smh1_1123, smi_1489, \
                         smi_1494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1916[k] = f_13 * sli_1237[k]
                    + f_3 * pc_y[k] * smi_1489[k];

        t_1917[k] = pb_y[k] * slk0_1593[k]
                    - f_12 * pc_y[k] * slk1_1593[k];

        t_1918[k] = f_8 * smh0_1123[k]
                    - f_9 * smh1_1123[k]
                    + f_3 * pc_x[k] * smi_1494[k];
    }

#pragma omp simd aligned(t_1919, t_1920, t_1921, pc_x, pc_y, pc_z, sli_1210, sli_1241, \
                         smh0_1125, smh1_1125, smi_1490, smi_1493, \
                         smi_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1919[k] = f_18 * sli_1210[k]
                    + f_3 * pc_z[k] * smi_1490[k];

        t_1920[k] = f_8 * smh0_1125[k]
                    - f_9 * smh1_1125[k]
                    + f_3 * pc_x[k] * smi_1496[k];

        t_1921[k] = f_13 * sli_1241[k]
                    + f_3 * pc_y[k] * smi_1493[k];
    }

#pragma omp simd aligned(t_1922, t_1923, t_1924, pb_y, pc_x, pc_y, pc_z, slk0_1598, sli_1214, \
                         slk1_1598, smh0_1128, smh1_1128, smi_1494, \
                         smi_1499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1922[k] = pb_y[k] * slk0_1598[k]
                    - f_12 * pc_y[k] * slk1_1598[k];

        t_1923[k] = f_10 * smh0_1128[k]
                    - f_11 * smh1_1128[k]
                    + f_3 * pc_x[k] * smi_1499[k];

        t_1924[k] = f_18 * sli_1214[k]
                    + f_3 * pc_z[k] * smi_1494[k];
    }

#pragma omp simd aligned(t_1925, t_1926, t_1927, pc_x, pc_y, sli_1246, smh0_1130, smh0_1131, \
                         smh1_1130, smh1_1131, smi_1498, smi_1501, \
                         smi_1502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1925[k] = f_10 * smh0_1130[k]
                    - f_11 * smh1_1130[k]
                    + f_3 * pc_x[k] * smi_1501[k];

        t_1926[k] = f_10 * smh0_1131[k]
                    - f_11 * smh1_1131[k]
                    + f_3 * pc_x[k] * smi_1502[k];

        t_1927[k] = f_13 * sli_1246[k]
                    + f_3 * pc_y[k] * smi_1498[k];
    }

#pragma omp simd aligned(t_1928, t_1929, t_1930, t_1931, t_1932, t_1933, pb_y, pc_x, pc_y, \
                         slk0_1604, slk1_1604, smi_1505, smi_1506, smi_1507, smi_1508, \
                         smi_1509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1928[k] = pb_y[k] * slk0_1604[k]
                    - f_12 * pc_y[k] * slk1_1604[k];

        t_1929[k] = f_3 * pc_x[k] * smi_1505[k];

        t_1930[k] = f_3 * pc_x[k] * smi_1506[k];

        t_1931[k] = f_3 * pc_x[k] * smi_1507[k];

        t_1932[k] = f_3 * pc_x[k] * smi_1508[k];

        t_1933[k] = f_3 * pc_x[k] * smi_1509[k];
    }

#pragma omp simd aligned(t_1934, t_1935, t_1936, t_1937, pb_y, pc_x, pc_y, pc_z, slk0_1612, \
                         sli_1225, sli_1253, slk1_1612, smi_1505, smi_1510, \
                         smi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1934[k] = f_3 * pc_x[k] * smi_1510[k];

        t_1935[k] = f_3 * pc_x[k] * smi_1511[k];

        t_1936[k] = pb_y[k] * slk0_1612[k]
                    + f_19 * sli_1253[k]
                    - f_12 * pc_y[k] * slk1_1612[k];

        t_1937[k] = f_18 * sli_1225[k]
                    + f_3 * pc_z[k] * smi_1505[k];
    }

#pragma omp simd aligned(t_1938, t_1939, t_1940, pb_y, pc_y, slk0_1614, slk0_1615, slk0_1616, \
                         sli_1255, sli_1256, sli_1257, slk1_1614, slk1_1615, \
                         slk1_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1938[k] = pb_y[k] * slk0_1614[k]
                    + f_17 * sli_1255[k]
                    - f_12 * pc_y[k] * slk1_1614[k];

        t_1939[k] = pb_y[k] * slk0_1615[k]
                    + f_16 * sli_1256[k]
                    - f_12 * pc_y[k] * slk1_1615[k];

        t_1940[k] = pb_y[k] * slk0_1616[k]
                    + f_15 * sli_1257[k]
                    - f_12 * pc_y[k] * slk1_1616[k];
    }
}

static auto
compute_prim_smk_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t slk0,
                                                           const size_t sli, const size_t slk1,
                                                           const size_t smh0, const size_t smh1,
                                                           const size_t smi, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 1.0 / q;

    auto *t_1941 = buffer.data(target + 1941);
    auto *t_1942 = buffer.data(target + 1942);
    auto *t_1943 = buffer.data(target + 1943);
    auto *t_1944 = buffer.data(target + 1944);
    auto *t_1945 = buffer.data(target + 1945);
    auto *t_1946 = buffer.data(target + 1946);
    auto *t_1947 = buffer.data(target + 1947);
    auto *t_1948 = buffer.data(target + 1948);
    auto *t_1949 = buffer.data(target + 1949);
    auto *t_1950 = buffer.data(target + 1950);
    auto *t_1951 = buffer.data(target + 1951);
    auto *t_1952 = buffer.data(target + 1952);
    auto *t_1953 = buffer.data(target + 1953);
    auto *t_1954 = buffer.data(target + 1954);
    auto *t_1955 = buffer.data(target + 1955);
    auto *t_1956 = buffer.data(target + 1956);
    auto *t_1957 = buffer.data(target + 1957);
    auto *t_1958 = buffer.data(target + 1958);
    auto *t_1959 = buffer.data(target + 1959);
    auto *t_1960 = buffer.data(target + 1960);
    auto *t_1961 = buffer.data(target + 1961);
    auto *t_1962 = buffer.data(target + 1962);
    auto *t_1963 = buffer.data(target + 1963);
    auto *t_1964 = buffer.data(target + 1964);
    auto *t_1965 = buffer.data(target + 1965);
    auto *t_1966 = buffer.data(target + 1966);
    auto *t_1967 = buffer.data(target + 1967);
    auto *t_1968 = buffer.data(target + 1968);
    auto *t_1969 = buffer.data(target + 1969);
    auto *t_1970 = buffer.data(target + 1970);
    auto *t_1971 = buffer.data(target + 1971);
    auto *t_1972 = buffer.data(target + 1972);
    auto *t_1973 = buffer.data(target + 1973);
    auto *t_1974 = buffer.data(target + 1974);
    auto *t_1975 = buffer.data(target + 1975);
    auto *t_1976 = buffer.data(target + 1976);
    auto *t_1977 = buffer.data(target + 1977);
    auto *t_1978 = buffer.data(target + 1978);
    auto *t_1979 = buffer.data(target + 1979);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk0_1617 = buffer.data(slk0 + 1617);
    const auto *slk0_1619 = buffer.data(slk0 + 1619);

    const auto *sli_1232 = buffer.data(sli + 1232);
    const auto *sli_1235 = buffer.data(sli + 1235);
    const auto *sli_1238 = buffer.data(sli + 1238);
    const auto *sli_1242 = buffer.data(sli + 1242);
    const auto *sli_1253 = buffer.data(sli + 1253);
    const auto *sli_1258 = buffer.data(sli + 1258);
    const auto *sli_1259 = buffer.data(sli + 1259);

    const auto *slk1_1617 = buffer.data(slk1 + 1617);
    const auto *slk1_1619 = buffer.data(slk1 + 1619);

    const auto *smh0_1134 = buffer.data(smh0 + 1134);
    const auto *smh0_1137 = buffer.data(smh0 + 1137);
    const auto *smh0_1139 = buffer.data(smh0 + 1139);
    const auto *smh0_1140 = buffer.data(smh0 + 1140);
    const auto *smh0_1143 = buffer.data(smh0 + 1143);
    const auto *smh0_1144 = buffer.data(smh0 + 1144);
    const auto *smh0_1146 = buffer.data(smh0 + 1146);
    const auto *smh0_1148 = buffer.data(smh0 + 1148);
    const auto *smh0_1149 = buffer.data(smh0 + 1149);
    const auto *smh0_1151 = buffer.data(smh0 + 1151);
    const auto *smh0_1152 = buffer.data(smh0 + 1152);
    const auto *smh0_1153 = buffer.data(smh0 + 1153);
    const auto *smh0_1154 = buffer.data(smh0 + 1154);

    const auto *smh1_1134 = buffer.data(smh1 + 1134);
    const auto *smh1_1137 = buffer.data(smh1 + 1137);
    const auto *smh1_1139 = buffer.data(smh1 + 1139);
    const auto *smh1_1140 = buffer.data(smh1 + 1140);
    const auto *smh1_1143 = buffer.data(smh1 + 1143);
    const auto *smh1_1144 = buffer.data(smh1 + 1144);
    const auto *smh1_1146 = buffer.data(smh1 + 1146);
    const auto *smh1_1148 = buffer.data(smh1 + 1148);
    const auto *smh1_1149 = buffer.data(smh1 + 1149);
    const auto *smh1_1151 = buffer.data(smh1 + 1151);
    const auto *smh1_1152 = buffer.data(smh1 + 1152);
    const auto *smh1_1153 = buffer.data(smh1 + 1153);
    const auto *smh1_1154 = buffer.data(smh1 + 1154);

    const auto *smi_1511 = buffer.data(smi + 1511);
    const auto *smi_1512 = buffer.data(smi + 1512);
    const auto *smi_1514 = buffer.data(smi + 1514);
    const auto *smi_1515 = buffer.data(smi + 1515);
    const auto *smi_1517 = buffer.data(smi + 1517);
    const auto *smi_1518 = buffer.data(smi + 1518);
    const auto *smi_1521 = buffer.data(smi + 1521);
    const auto *smi_1522 = buffer.data(smi + 1522);
    const auto *smi_1524 = buffer.data(smi + 1524);
    const auto *smi_1526 = buffer.data(smi + 1526);
    const auto *smi_1527 = buffer.data(smi + 1527);
    const auto *smi_1529 = buffer.data(smi + 1529);
    const auto *smi_1530 = buffer.data(smi + 1530);
    const auto *smi_1532 = buffer.data(smi + 1532);
    const auto *smi_1533 = buffer.data(smi + 1533);
    const auto *smi_1534 = buffer.data(smi + 1534);
    const auto *smi_1535 = buffer.data(smi + 1535);
    const auto *smi_1536 = buffer.data(smi + 1536);
    const auto *smi_1537 = buffer.data(smi + 1537);
    const auto *smi_1538 = buffer.data(smi + 1538);
    const auto *smi_1539 = buffer.data(smi + 1539);

#pragma omp simd aligned(t_1941, t_1942, t_1943, pb_y, pc_y, slk0_1617, slk0_1619, sli_1258, \
                         sli_1259, slk1_1617, slk1_1619, smi_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1941[k] = pb_y[k] * slk0_1617[k]
                    + f_14 * sli_1258[k]
                    - f_12 * pc_y[k] * slk1_1617[k];

        t_1942[k] = f_13 * sli_1259[k]
                    + f_3 * pc_y[k] * smi_1511[k];

        t_1943[k] = pb_y[k] * slk0_1619[k]
                    - f_12 * pc_y[k] * slk1_1619[k];
    }

#pragma omp simd aligned(t_1944, t_1945, t_1946, t_1947, t_1948, pc_x, pc_y, pc_z, sli_1232, \
                         smh0_1134, smh0_1137, smh1_1134, smh1_1137, smi_1512, smi_1514, \
                         smi_1515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1944[k] = f_1 * smh0_1134[k]
                    - f_2 * smh1_1134[k]
                    + f_3 * pc_x[k] * smi_1512[k];

        t_1945[k] = f_3 * pc_y[k] * smi_1512[k];

        t_1946[k] = f_0 * sli_1232[k]
                    + f_3 * pc_z[k] * smi_1512[k];

        t_1947[k] = f_4 * smh0_1137[k]
                    - f_5 * smh1_1137[k]
                    + f_3 * pc_x[k] * smi_1515[k];

        t_1948[k] = f_3 * pc_y[k] * smi_1514[k];
    }

#pragma omp simd aligned(t_1949, t_1950, t_1951, t_1952, pc_x, pc_y, pc_z, sli_1235, \
                         smh0_1139, smh0_1140, smh1_1139, smh1_1140, smi_1515, smi_1517, \
                         smi_1518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1949[k] = f_4 * smh0_1139[k]
                    - f_5 * smh1_1139[k]
                    + f_3 * pc_x[k] * smi_1517[k];

        t_1950[k] = f_6 * smh0_1140[k]
                    - f_7 * smh1_1140[k]
                    + f_3 * pc_x[k] * smi_1518[k];

        t_1951[k] = f_0 * sli_1235[k]
                    + f_3 * pc_z[k] * smi_1515[k];

        t_1952[k] = f_3 * pc_y[k] * smi_1517[k];
    }

#pragma omp simd aligned(t_1953, t_1954, t_1955, pc_x, pc_z, sli_1238, smh0_1143, smh0_1144, \
                         smh1_1143, smh1_1144, smi_1518, smi_1521, \
                         smi_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1953[k] = f_6 * smh0_1143[k]
                    - f_7 * smh1_1143[k]
                    + f_3 * pc_x[k] * smi_1521[k];

        t_1954[k] = f_8 * smh0_1144[k]
                    - f_9 * smh1_1144[k]
                    + f_3 * pc_x[k] * smi_1522[k];

        t_1955[k] = f_0 * sli_1238[k]
                    + f_3 * pc_z[k] * smi_1518[k];
    }

#pragma omp simd aligned(t_1956, t_1957, t_1958, t_1959, pc_x, pc_y, smh0_1146, smh0_1148, \
                         smh0_1149, smh1_1146, smh1_1148, smh1_1149, smi_1521, smi_1524, \
                         smi_1526, smi_1527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1956[k] = f_8 * smh0_1146[k]
                    - f_9 * smh1_1146[k]
                    + f_3 * pc_x[k] * smi_1524[k];

        t_1957[k] = f_3 * pc_y[k] * smi_1521[k];

        t_1958[k] = f_8 * smh0_1148[k]
                    - f_9 * smh1_1148[k]
                    + f_3 * pc_x[k] * smi_1526[k];

        t_1959[k] = f_10 * smh0_1149[k]
                    - f_11 * smh1_1149[k]
                    + f_3 * pc_x[k] * smi_1527[k];
    }

#pragma omp simd aligned(t_1960, t_1961, t_1962, t_1963, pc_x, pc_y, pc_z, sli_1242, \
                         smh0_1151, smh0_1152, smh1_1151, smh1_1152, smi_1522, smi_1526, \
                         smi_1529, smi_1530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1960[k] = f_0 * sli_1242[k]
                    + f_3 * pc_z[k] * smi_1522[k];

        t_1961[k] = f_10 * smh0_1151[k]
                    - f_11 * smh1_1151[k]
                    + f_3 * pc_x[k] * smi_1529[k];

        t_1962[k] = f_10 * smh0_1152[k]
                    - f_11 * smh1_1152[k]
                    + f_3 * pc_x[k] * smi_1530[k];

        t_1963[k] = f_3 * pc_y[k] * smi_1526[k];
    }

#pragma omp simd aligned(t_1964, t_1965, t_1966, t_1967, t_1968, t_1969, pc_x, smh0_1154, \
                         smh1_1154, smi_1532, smi_1533, smi_1534, smi_1535, smi_1536, \
                         smi_1537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1964[k] = f_10 * smh0_1154[k]
                    - f_11 * smh1_1154[k]
                    + f_3 * pc_x[k] * smi_1532[k];

        t_1965[k] = f_3 * pc_x[k] * smi_1533[k];

        t_1966[k] = f_3 * pc_x[k] * smi_1534[k];

        t_1967[k] = f_3 * pc_x[k] * smi_1535[k];

        t_1968[k] = f_3 * pc_x[k] * smi_1536[k];

        t_1969[k] = f_3 * pc_x[k] * smi_1537[k];
    }

#pragma omp simd aligned(t_1970, t_1971, t_1972, t_1973, pc_x, pc_y, pc_z, sli_1253, \
                         smh0_1149, smh1_1149, smi_1533, smi_1538, \
                         smi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1970[k] = f_3 * pc_x[k] * smi_1538[k];

        t_1971[k] = f_3 * pc_x[k] * smi_1539[k];

        t_1972[k] = f_1 * smh0_1149[k]
                    - f_2 * smh1_1149[k]
                    + f_3 * pc_y[k] * smi_1533[k];

        t_1973[k] = f_0 * sli_1253[k]
                    + f_3 * pc_z[k] * smi_1533[k];
    }

#pragma omp simd aligned(t_1974, t_1975, t_1976, pc_y, smh0_1151, smh0_1152, smh0_1153, \
                         smh1_1151, smh1_1152, smh1_1153, smi_1535, smi_1536, \
                         smi_1537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1974[k] = f_4 * smh0_1151[k]
                    - f_5 * smh1_1151[k]
                    + f_3 * pc_y[k] * smi_1535[k];

        t_1975[k] = f_6 * smh0_1152[k]
                    - f_7 * smh1_1152[k]
                    + f_3 * pc_y[k] * smi_1536[k];

        t_1976[k] = f_8 * smh0_1153[k]
                    - f_9 * smh1_1153[k]
                    + f_3 * pc_y[k] * smi_1537[k];
    }

#pragma omp simd aligned(t_1977, t_1978, t_1979, pc_y, pc_z, sli_1259, smh0_1154, smh1_1154, \
                         smi_1538, smi_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1977[k] = f_10 * smh0_1154[k]
                    - f_11 * smh1_1154[k]
                    + f_3 * pc_y[k] * smi_1538[k];

        t_1978[k] = f_3 * pc_y[k] * smi_1539[k];

        t_1979[k] = f_0 * sli_1259[k]
                    + f_1 * smh0_1154[k]
                    - f_2 * smh1_1154[k]
                    + f_3 * pc_z[k] * smi_1539[k];
    }
}

auto
compute_prim_smk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t slk0, const size_t sli,
                                                   const size_t slk1, const size_t smh0,
                                                   const size_t smh1, const size_t smi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_smk_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, slk0, sli,
                                                              slk1, smh0, smh1, smi, ncols,
                                                              gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece10(buffer, target, pc, sli, smh0,
                                                               smh1, smi, ncols, gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, slk0,
                                                               sli, slk1, smh0, smh1, smi,
                                                               ncols, gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece12(buffer, target, pb, pc, slk0,
                                                               sli, slk1, smi, ncols, gamma, p,
                                                               q);

    compute_prim_smk_three_center_electron_repulsion_0_piece13(buffer, target, pb, pc, slk0,
                                                               sli, slk1, smi, ncols, gamma, p,
                                                               q);

    compute_prim_smk_three_center_electron_repulsion_0_piece14(buffer, target, pb, pc, slk0,
                                                               sli, slk1, smh0, smh1, smi,
                                                               ncols, gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece15(buffer, target, pc, sli, smh0,
                                                               smh1, smi, ncols, gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece16(buffer, target, pb, pc, slk0,
                                                               sli, slk1, smh0, smh1, smi,
                                                               ncols, gamma, p, q);

    compute_prim_smk_three_center_electron_repulsion_0_piece17(buffer, target, pb, pc, slk0,
                                                               sli, slk1, smh0, smh1, smi,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
