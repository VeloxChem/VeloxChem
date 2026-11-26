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


#include "SimdThreeCenterElectronRepulsionVrrRecSNK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;

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

    const auto *smk0_0 = buffer.data(smk0 + 0);
    const auto *smk0_3 = buffer.data(smk0 + 3);
    const auto *smk0_5 = buffer.data(smk0 + 5);
    const auto *smk0_6 = buffer.data(smk0 + 6);
    const auto *smk0_9 = buffer.data(smk0 + 9);
    const auto *smk0_10 = buffer.data(smk0 + 10);
    const auto *smk0_12 = buffer.data(smk0 + 12);
    const auto *smk0_14 = buffer.data(smk0 + 14);
    const auto *smk0_15 = buffer.data(smk0 + 15);
    const auto *smk0_17 = buffer.data(smk0 + 17);
    const auto *smk0_18 = buffer.data(smk0 + 18);
    const auto *smk0_20 = buffer.data(smk0 + 20);
    const auto *smk0_28 = buffer.data(smk0 + 28);
    const auto *smk0_35 = buffer.data(smk0 + 35);

    const auto *smi_0 = buffer.data(smi + 0);
    const auto *smi_1 = buffer.data(smi + 1);
    const auto *smi_2 = buffer.data(smi + 2);
    const auto *smi_3 = buffer.data(smi + 3);
    const auto *smi_5 = buffer.data(smi + 5);
    const auto *smi_6 = buffer.data(smi + 6);
    const auto *smi_7 = buffer.data(smi + 7);
    const auto *smi_8 = buffer.data(smi + 8);
    const auto *smi_9 = buffer.data(smi + 9);
    const auto *smi_10 = buffer.data(smi + 10);
    const auto *smi_11 = buffer.data(smi + 11);
    const auto *smi_12 = buffer.data(smi + 12);
    const auto *smi_13 = buffer.data(smi + 13);
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
    const auto *smi_33 = buffer.data(smi + 33);
    const auto *smi_37 = buffer.data(smi + 37);
    const auto *smi_49 = buffer.data(smi + 49);
    const auto *smi_50 = buffer.data(smi + 50);
    const auto *smi_51 = buffer.data(smi + 51);
    const auto *smi_52 = buffer.data(smi + 52);
    const auto *smi_53 = buffer.data(smi + 53);
    const auto *smi_54 = buffer.data(smi + 54);
    const auto *smi_55 = buffer.data(smi + 55);
    const auto *smi_77 = buffer.data(smi + 77);
    const auto *smi_78 = buffer.data(smi + 78);
    const auto *smi_79 = buffer.data(smi + 79);
    const auto *smi_80 = buffer.data(smi + 80);
    const auto *smi_81 = buffer.data(smi + 81);
    const auto *smi_82 = buffer.data(smi + 82);
    const auto *smi_83 = buffer.data(smi + 83);
    const auto *smi_84 = buffer.data(smi + 84);
    const auto *smi_87 = buffer.data(smi + 87);
    const auto *smi_89 = buffer.data(smi + 89);
    const auto *smi_90 = buffer.data(smi + 90);
    const auto *smi_93 = buffer.data(smi + 93);
    const auto *smi_94 = buffer.data(smi + 94);
    const auto *smi_96 = buffer.data(smi + 96);

    const auto *smk1_0 = buffer.data(smk1 + 0);
    const auto *smk1_3 = buffer.data(smk1 + 3);
    const auto *smk1_5 = buffer.data(smk1 + 5);
    const auto *smk1_6 = buffer.data(smk1 + 6);
    const auto *smk1_9 = buffer.data(smk1 + 9);
    const auto *smk1_10 = buffer.data(smk1 + 10);
    const auto *smk1_12 = buffer.data(smk1 + 12);
    const auto *smk1_14 = buffer.data(smk1 + 14);
    const auto *smk1_15 = buffer.data(smk1 + 15);
    const auto *smk1_17 = buffer.data(smk1 + 17);
    const auto *smk1_18 = buffer.data(smk1 + 18);
    const auto *smk1_20 = buffer.data(smk1 + 20);
    const auto *smk1_28 = buffer.data(smk1 + 28);
    const auto *smk1_35 = buffer.data(smk1 + 35);

    const auto *snh0_0 = buffer.data(snh0 + 0);
    const auto *snh0_3 = buffer.data(snh0 + 3);
    const auto *snh0_5 = buffer.data(snh0 + 5);
    const auto *snh0_6 = buffer.data(snh0 + 6);
    const auto *snh0_9 = buffer.data(snh0 + 9);
    const auto *snh0_10 = buffer.data(snh0 + 10);
    const auto *snh0_12 = buffer.data(snh0 + 12);
    const auto *snh0_14 = buffer.data(snh0 + 14);
    const auto *snh0_15 = buffer.data(snh0 + 15);
    const auto *snh0_17 = buffer.data(snh0 + 17);
    const auto *snh0_18 = buffer.data(snh0 + 18);
    const auto *snh0_19 = buffer.data(snh0 + 19);
    const auto *snh0_20 = buffer.data(snh0 + 20);
    const auto *snh0_36 = buffer.data(snh0 + 36);
    const auto *snh0_38 = buffer.data(snh0 + 38);
    const auto *snh0_39 = buffer.data(snh0 + 39);
    const auto *snh0_40 = buffer.data(snh0 + 40);
    const auto *snh0_41 = buffer.data(snh0 + 41);
    const auto *snh0_59 = buffer.data(snh0 + 59);
    const auto *snh0_60 = buffer.data(snh0 + 60);
    const auto *snh0_61 = buffer.data(snh0 + 61);
    const auto *snh0_62 = buffer.data(snh0 + 62);
    const auto *snh0_63 = buffer.data(snh0 + 63);
    const auto *snh0_66 = buffer.data(snh0 + 66);
    const auto *snh0_68 = buffer.data(snh0 + 68);
    const auto *snh0_69 = buffer.data(snh0 + 69);
    const auto *snh0_72 = buffer.data(snh0 + 72);
    const auto *snh0_73 = buffer.data(snh0 + 73);
    const auto *snh0_75 = buffer.data(snh0 + 75);

    const auto *snh1_0 = buffer.data(snh1 + 0);
    const auto *snh1_3 = buffer.data(snh1 + 3);
    const auto *snh1_5 = buffer.data(snh1 + 5);
    const auto *snh1_6 = buffer.data(snh1 + 6);
    const auto *snh1_9 = buffer.data(snh1 + 9);
    const auto *snh1_10 = buffer.data(snh1 + 10);
    const auto *snh1_12 = buffer.data(snh1 + 12);
    const auto *snh1_14 = buffer.data(snh1 + 14);
    const auto *snh1_15 = buffer.data(snh1 + 15);
    const auto *snh1_17 = buffer.data(snh1 + 17);
    const auto *snh1_18 = buffer.data(snh1 + 18);
    const auto *snh1_19 = buffer.data(snh1 + 19);
    const auto *snh1_20 = buffer.data(snh1 + 20);
    const auto *snh1_36 = buffer.data(snh1 + 36);
    const auto *snh1_38 = buffer.data(snh1 + 38);
    const auto *snh1_39 = buffer.data(snh1 + 39);
    const auto *snh1_40 = buffer.data(snh1 + 40);
    const auto *snh1_41 = buffer.data(snh1 + 41);
    const auto *snh1_59 = buffer.data(snh1 + 59);
    const auto *snh1_60 = buffer.data(snh1 + 60);
    const auto *snh1_61 = buffer.data(snh1 + 61);
    const auto *snh1_62 = buffer.data(snh1 + 62);
    const auto *snh1_63 = buffer.data(snh1 + 63);
    const auto *snh1_66 = buffer.data(snh1 + 66);
    const auto *snh1_68 = buffer.data(snh1 + 68);
    const auto *snh1_69 = buffer.data(snh1 + 69);
    const auto *snh1_72 = buffer.data(snh1 + 72);
    const auto *snh1_73 = buffer.data(snh1 + 73);
    const auto *snh1_75 = buffer.data(snh1 + 75);

    const auto *sni_0 = buffer.data(sni + 0);
    const auto *sni_2 = buffer.data(sni + 2);
    const auto *sni_3 = buffer.data(sni + 3);
    const auto *sni_5 = buffer.data(sni + 5);
    const auto *sni_6 = buffer.data(sni + 6);
    const auto *sni_9 = buffer.data(sni + 9);
    const auto *sni_10 = buffer.data(sni + 10);
    const auto *sni_12 = buffer.data(sni + 12);
    const auto *sni_14 = buffer.data(sni + 14);
    const auto *sni_15 = buffer.data(sni + 15);
    const auto *sni_17 = buffer.data(sni + 17);
    const auto *sni_18 = buffer.data(sni + 18);
    const auto *sni_20 = buffer.data(sni + 20);
    const auto *sni_21 = buffer.data(sni + 21);
    const auto *sni_22 = buffer.data(sni + 22);
    const auto *sni_23 = buffer.data(sni + 23);
    const auto *sni_24 = buffer.data(sni + 24);
    const auto *sni_25 = buffer.data(sni + 25);
    const auto *sni_26 = buffer.data(sni + 26);
    const auto *sni_27 = buffer.data(sni + 27);
    const auto *sni_28 = buffer.data(sni + 28);
    const auto *sni_30 = buffer.data(sni + 30);
    const auto *sni_31 = buffer.data(sni + 31);
    const auto *sni_33 = buffer.data(sni + 33);
    const auto *sni_34 = buffer.data(sni + 34);
    const auto *sni_37 = buffer.data(sni + 37);
    const auto *sni_38 = buffer.data(sni + 38);
    const auto *sni_42 = buffer.data(sni + 42);
    const auto *sni_49 = buffer.data(sni + 49);
    const auto *sni_50 = buffer.data(sni + 50);
    const auto *sni_51 = buffer.data(sni + 51);
    const auto *sni_52 = buffer.data(sni + 52);
    const auto *sni_53 = buffer.data(sni + 53);
    const auto *sni_54 = buffer.data(sni + 54);
    const auto *sni_55 = buffer.data(sni + 55);
    const auto *sni_56 = buffer.data(sni + 56);
    const auto *sni_58 = buffer.data(sni + 58);
    const auto *sni_59 = buffer.data(sni + 59);
    const auto *sni_61 = buffer.data(sni + 61);
    const auto *sni_62 = buffer.data(sni + 62);
    const auto *sni_65 = buffer.data(sni + 65);
    const auto *sni_66 = buffer.data(sni + 66);
    const auto *sni_70 = buffer.data(sni + 70);
    const auto *sni_77 = buffer.data(sni + 77);
    const auto *sni_78 = buffer.data(sni + 78);
    const auto *sni_79 = buffer.data(sni + 79);
    const auto *sni_80 = buffer.data(sni + 80);
    const auto *sni_81 = buffer.data(sni + 81);
    const auto *sni_82 = buffer.data(sni + 82);
    const auto *sni_83 = buffer.data(sni + 83);
    const auto *sni_84 = buffer.data(sni + 84);
    const auto *sni_86 = buffer.data(sni + 86);
    const auto *sni_87 = buffer.data(sni + 87);
    const auto *sni_89 = buffer.data(sni + 89);
    const auto *sni_90 = buffer.data(sni + 90);
    const auto *sni_93 = buffer.data(sni + 93);
    const auto *sni_94 = buffer.data(sni + 94);
    const auto *sni_96 = buffer.data(sni + 96);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, smi_0, smi_3, snh0_0, snh0_3, \
                         snh1_0, snh1_3, sni_0, sni_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * smi_0[k]
                 + f_1 * snh0_0[k]
                 - f_2 * snh1_0[k]
                 + f_3 * pc_x[k] * sni_0[k];

        t_1[k] = f_3 * pc_y[k] * sni_0[k];

        t_2[k] = f_3 * pc_z[k] * sni_0[k];

        t_3[k] = f_0 * smi_3[k]
                 + f_4 * snh0_3[k]
                 - f_5 * snh1_3[k]
                 + f_3 * pc_x[k] * sni_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, smi_5, smi_6, snh0_5, snh0_6, snh1_5, \
                         snh1_6, sni_2, sni_5, sni_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sni_2[k];

        t_5[k] = f_0 * smi_5[k]
                 + f_4 * snh0_5[k]
                 - f_5 * snh1_5[k]
                 + f_3 * pc_x[k] * sni_5[k];

        t_6[k] = f_0 * smi_6[k]
                 + f_6 * snh0_6[k]
                 - f_7 * snh1_6[k]
                 + f_3 * pc_x[k] * sni_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, smi_9, snh0_9, snh1_9, sni_3, sni_5, \
                         sni_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sni_3[k];

        t_8[k] = f_3 * pc_y[k] * sni_5[k];

        t_9[k] = f_0 * smi_9[k]
                 + f_6 * snh0_9[k]
                 - f_7 * snh1_9[k]
                 + f_3 * pc_x[k] * sni_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, smi_10, smi_12, snh0_10, snh0_12, \
                         snh1_10, snh1_12, sni_6, sni_10, sni_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * smi_10[k]
                  + f_8 * snh0_10[k]
                  - f_9 * snh1_10[k]
                  + f_3 * pc_x[k] * sni_10[k];

        t_11[k] = f_3 * pc_z[k] * sni_6[k];

        t_12[k] = f_0 * smi_12[k]
                  + f_8 * snh0_12[k]
                  - f_9 * snh1_12[k]
                  + f_3 * pc_x[k] * sni_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, smi_14, smi_15, snh0_14, snh0_15, \
                         snh1_14, snh1_15, sni_9, sni_14, sni_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sni_9[k];

        t_14[k] = f_0 * smi_14[k]
                  + f_8 * snh0_14[k]
                  - f_9 * snh1_14[k]
                  + f_3 * pc_x[k] * sni_14[k];

        t_15[k] = f_0 * smi_15[k]
                  + f_10 * snh0_15[k]
                  - f_11 * snh1_15[k]
                  + f_3 * pc_x[k] * sni_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, smi_17, smi_18, snh0_17, snh0_18, \
                         snh1_17, snh1_18, sni_10, sni_17, sni_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sni_10[k];

        t_17[k] = f_0 * smi_17[k]
                  + f_10 * snh0_17[k]
                  - f_11 * snh1_17[k]
                  + f_3 * pc_x[k] * sni_17[k];

        t_18[k] = f_0 * smi_18[k]
                  + f_10 * snh0_18[k]
                  - f_11 * snh1_18[k]
                  + f_3 * pc_x[k] * sni_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, smi_20, smi_21, smi_22, snh0_20, \
                         snh1_20, sni_14, sni_20, sni_21, sni_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sni_14[k];

        t_20[k] = f_0 * smi_20[k]
                  + f_10 * snh0_20[k]
                  - f_11 * snh1_20[k]
                  + f_3 * pc_x[k] * sni_20[k];

        t_21[k] = f_0 * smi_21[k]
                  + f_3 * pc_x[k] * sni_21[k];

        t_22[k] = f_0 * smi_22[k]
                  + f_3 * pc_x[k] * sni_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, smi_23, smi_24, smi_25, smi_26, \
                         smi_27, sni_23, sni_24, sni_25, sni_26, \
                         sni_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * smi_23[k]
                  + f_3 * pc_x[k] * sni_23[k];

        t_24[k] = f_0 * smi_24[k]
                  + f_3 * pc_x[k] * sni_24[k];

        t_25[k] = f_0 * smi_25[k]
                  + f_3 * pc_x[k] * sni_25[k];

        t_26[k] = f_0 * smi_26[k]
                  + f_3 * pc_x[k] * sni_26[k];

        t_27[k] = f_0 * smi_27[k]
                  + f_3 * pc_x[k] * sni_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, snh0_15, snh0_17, snh0_18, \
                         snh1_15, snh1_17, snh1_18, sni_21, sni_23, \
                         sni_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * snh0_15[k]
                  - f_2 * snh1_15[k]
                  + f_3 * pc_y[k] * sni_21[k];

        t_29[k] = f_3 * pc_z[k] * sni_21[k];

        t_30[k] = f_4 * snh0_17[k]
                  - f_5 * snh1_17[k]
                  + f_3 * pc_y[k] * sni_23[k];

        t_31[k] = f_6 * snh0_18[k]
                  - f_7 * snh1_18[k]
                  + f_3 * pc_y[k] * sni_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, snh0_19, snh0_20, snh1_19, \
                         snh1_20, sni_25, sni_26, sni_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * snh0_19[k]
                  - f_9 * snh1_19[k]
                  + f_3 * pc_y[k] * sni_25[k];

        t_33[k] = f_10 * snh0_20[k]
                  - f_11 * snh1_20[k]
                  + f_3 * pc_y[k] * sni_26[k];

        t_34[k] = f_3 * pc_y[k] * sni_27[k];

        t_35[k] = f_1 * snh0_20[k]
                  - f_2 * snh1_20[k]
                  + f_3 * pc_z[k] * sni_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, smk0_0, smk0_3, smi_0, \
                         smi_1, smk1_0, smk1_3, sni_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * smk0_0[k]
                  - f_12 * pc_y[k] * smk1_0[k];

        t_37[k] = f_13 * smi_0[k]
                  + f_3 * pc_y[k] * sni_28[k];

        t_38[k] = f_3 * pc_z[k] * sni_28[k];

        t_39[k] = pb_y[k] * smk0_3[k]
                  + f_14 * smi_1[k]
                  - f_12 * pc_y[k] * smk1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_y, pc_z, smk0_5, smk0_6, smi_2, \
                         smi_3, smk1_5, smk1_6, sni_30, sni_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * smi_2[k]
                  + f_3 * pc_y[k] * sni_30[k];

        t_41[k] = pb_y[k] * smk0_5[k]
                  - f_12 * pc_y[k] * smk1_5[k];

        t_42[k] = pb_y[k] * smk0_6[k]
                  + f_15 * smi_3[k]
                  - f_12 * pc_y[k] * smk1_6[k];

        t_43[k] = f_3 * pc_z[k] * sni_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, smk0_9, smk0_10, smi_5, \
                         smi_6, smk1_9, smk1_10, sni_33, sni_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * smi_5[k]
                  + f_3 * pc_y[k] * sni_33[k];

        t_45[k] = pb_y[k] * smk0_9[k]
                  - f_12 * pc_y[k] * smk1_9[k];

        t_46[k] = pb_y[k] * smk0_10[k]
                  + f_16 * smi_6[k]
                  - f_12 * pc_y[k] * smk1_10[k];

        t_47[k] = f_3 * pc_z[k] * sni_34[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_y, pc_y, smk0_12, smk0_14, smk0_15, smi_8, \
                         smi_9, smi_10, smk1_12, smk1_14, smk1_15, \
                         sni_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pb_y[k] * smk0_12[k]
                  + f_14 * smi_8[k]
                  - f_12 * pc_y[k] * smk1_12[k];

        t_49[k] = f_13 * smi_9[k]
                  + f_3 * pc_y[k] * sni_37[k];

        t_50[k] = pb_y[k] * smk0_14[k]
                  - f_12 * pc_y[k] * smk1_14[k];

        t_51[k] = pb_y[k] * smk0_15[k]
                  + f_17 * smi_10[k]
                  - f_12 * pc_y[k] * smk1_15[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, pc_z, smk0_17, smk0_18, smi_12, \
                         smi_13, smi_14, smk1_17, smk1_18, sni_38, \
                         sni_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * sni_38[k];

        t_53[k] = pb_y[k] * smk0_17[k]
                  + f_15 * smi_12[k]
                  - f_12 * pc_y[k] * smk1_17[k];

        t_54[k] = pb_y[k] * smk0_18[k]
                  + f_14 * smi_13[k]
                  - f_12 * pc_y[k] * smk1_18[k];

        t_55[k] = f_13 * smi_14[k]
                  + f_3 * pc_y[k] * sni_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pc_x, pc_y, smk0_20, smi_49, smi_50, \
                         smi_51, smk1_20, sni_49, sni_50, sni_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * smk0_20[k]
                  - f_12 * pc_y[k] * smk1_20[k];

        t_57[k] = f_18 * smi_49[k]
                  + f_3 * pc_x[k] * sni_49[k];

        t_58[k] = f_18 * smi_50[k]
                  + f_3 * pc_x[k] * sni_50[k];

        t_59[k] = f_18 * smi_51[k]
                  + f_3 * pc_x[k] * sni_51[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, smi_52, smi_53, smi_54, smi_55, sni_52, \
                         sni_53, sni_54, sni_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_18 * smi_52[k]
                  + f_3 * pc_x[k] * sni_52[k];

        t_61[k] = f_18 * smi_53[k]
                  + f_3 * pc_x[k] * sni_53[k];

        t_62[k] = f_18 * smi_54[k]
                  + f_3 * pc_x[k] * sni_54[k];

        t_63[k] = f_18 * smi_55[k]
                  + f_3 * pc_x[k] * sni_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pc_y, pc_z, smi_21, smi_23, snh0_36, snh0_38, \
                         snh1_36, snh1_38, sni_49, sni_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_13 * smi_21[k]
                  + f_1 * snh0_36[k]
                  - f_2 * snh1_36[k]
                  + f_3 * pc_y[k] * sni_49[k];

        t_65[k] = f_3 * pc_z[k] * sni_49[k];

        t_66[k] = f_13 * smi_23[k]
                  + f_4 * snh0_38[k]
                  - f_5 * snh1_38[k]
                  + f_3 * pc_y[k] * sni_51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, smi_24, smi_25, smi_26, snh0_39, snh0_40, \
                         snh0_41, snh1_39, snh1_40, snh1_41, sni_52, sni_53, \
                         sni_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_13 * smi_24[k]
                  + f_6 * snh0_39[k]
                  - f_7 * snh1_39[k]
                  + f_3 * pc_y[k] * sni_52[k];

        t_68[k] = f_13 * smi_25[k]
                  + f_8 * snh0_40[k]
                  - f_9 * snh1_40[k]
                  + f_3 * pc_y[k] * sni_53[k];

        t_69[k] = f_13 * smi_26[k]
                  + f_10 * snh0_41[k]
                  - f_11 * snh1_41[k]
                  + f_3 * pc_y[k] * sni_54[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_y, pb_z, pc_y, pc_z, smk0_0, smk0_35, \
                         smi_27, smk1_0, smk1_35, sni_55, sni_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_13 * smi_27[k]
                  + f_3 * pc_y[k] * sni_55[k];

        t_71[k] = pb_y[k] * smk0_35[k]
                  - f_12 * pc_y[k] * smk1_35[k];

        t_72[k] = pb_z[k] * smk0_0[k]
                  - f_12 * pc_z[k] * smk1_0[k];

        t_73[k] = f_3 * pc_y[k] * sni_56[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_y, pc_z, smk0_3, smk0_5, smi_0, \
                         smi_2, smk1_3, smk1_5, sni_56, sni_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_13 * smi_0[k]
                  + f_3 * pc_z[k] * sni_56[k];

        t_75[k] = pb_z[k] * smk0_3[k]
                  - f_12 * pc_z[k] * smk1_3[k];

        t_76[k] = f_3 * pc_y[k] * sni_58[k];

        t_77[k] = pb_z[k] * smk0_5[k]
                  + f_14 * smi_2[k]
                  - f_12 * pc_z[k] * smk1_5[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_z, pc_y, pc_z, smk0_6, smk0_9, smi_3, \
                         smi_5, smk1_6, smk1_9, sni_59, sni_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_z[k] * smk0_6[k]
                  - f_12 * pc_z[k] * smk1_6[k];

        t_79[k] = f_13 * smi_3[k]
                  + f_3 * pc_z[k] * sni_59[k];

        t_80[k] = f_3 * pc_y[k] * sni_61[k];

        t_81[k] = pb_z[k] * smk0_9[k]
                  + f_15 * smi_5[k]
                  - f_12 * pc_z[k] * smk1_9[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_z, pc_y, pc_z, smk0_10, smk0_12, smi_6, \
                         smi_7, smk1_10, smk1_12, sni_62, sni_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * smk0_10[k]
                  - f_12 * pc_z[k] * smk1_10[k];

        t_83[k] = f_13 * smi_6[k]
                  + f_3 * pc_z[k] * sni_62[k];

        t_84[k] = pb_z[k] * smk0_12[k]
                  + f_14 * smi_7[k]
                  - f_12 * pc_z[k] * smk1_12[k];

        t_85[k] = f_3 * pc_y[k] * sni_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_z, pc_z, smk0_14, smk0_15, smk0_17, smi_9, \
                         smi_10, smi_11, smk1_14, smk1_15, smk1_17, \
                         sni_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * smk0_14[k]
                  + f_16 * smi_9[k]
                  - f_12 * pc_z[k] * smk1_14[k];

        t_87[k] = pb_z[k] * smk0_15[k]
                  - f_12 * pc_z[k] * smk1_15[k];

        t_88[k] = f_13 * smi_10[k]
                  + f_3 * pc_z[k] * sni_66[k];

        t_89[k] = pb_z[k] * smk0_17[k]
                  + f_14 * smi_11[k]
                  - f_12 * pc_z[k] * smk1_17[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_z, pc_y, pc_z, smk0_18, smk0_20, smi_12, smi_14, \
                         smk1_18, smk1_20, sni_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * smk0_18[k]
                  + f_15 * smi_12[k]
                  - f_12 * pc_z[k] * smk1_18[k];

        t_91[k] = f_3 * pc_y[k] * sni_70[k];

        t_92[k] = pb_z[k] * smk0_20[k]
                  + f_17 * smi_14[k]
                  - f_12 * pc_z[k] * smk1_20[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, smi_77, smi_78, smi_79, smi_80, \
                         smi_81, sni_77, sni_78, sni_79, sni_80, \
                         sni_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_18 * smi_77[k]
                  + f_3 * pc_x[k] * sni_77[k];

        t_94[k] = f_18 * smi_78[k]
                  + f_3 * pc_x[k] * sni_78[k];

        t_95[k] = f_18 * smi_79[k]
                  + f_3 * pc_x[k] * sni_79[k];

        t_96[k] = f_18 * smi_80[k]
                  + f_3 * pc_x[k] * sni_80[k];

        t_97[k] = f_18 * smi_81[k]
                  + f_3 * pc_x[k] * sni_81[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_z, smk0_28, smi_21, smi_82, \
                         smi_83, smk1_28, sni_77, sni_82, sni_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_18 * smi_82[k]
                  + f_3 * pc_x[k] * sni_82[k];

        t_99[k] = f_18 * smi_83[k]
                  + f_3 * pc_x[k] * sni_83[k];

        t_100[k] = pb_z[k] * smk0_28[k]
                   - f_12 * pc_z[k] * smk1_28[k];

        t_101[k] = f_13 * smi_21[k]
                   + f_3 * pc_z[k] * sni_77[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, snh0_59, snh0_60, snh0_61, snh1_59, \
                         snh1_60, snh1_61, sni_79, sni_80, sni_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * snh0_59[k]
                   - f_5 * snh1_59[k]
                   + f_3 * pc_y[k] * sni_79[k];

        t_103[k] = f_6 * snh0_60[k]
                   - f_7 * snh1_60[k]
                   + f_3 * pc_y[k] * sni_80[k];

        t_104[k] = f_8 * snh0_61[k]
                   - f_9 * snh1_61[k]
                   + f_3 * pc_y[k] * sni_81[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_x, pc_y, pc_z, smi_27, smi_84, \
                         snh0_62, snh0_63, snh1_62, snh1_63, sni_82, sni_83, \
                         sni_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_10 * snh0_62[k]
                   - f_11 * snh1_62[k]
                   + f_3 * pc_y[k] * sni_82[k];

        t_106[k] = f_3 * pc_y[k] * sni_83[k];

        t_107[k] = f_13 * smi_27[k]
                   + f_1 * snh0_62[k]
                   - f_2 * snh1_62[k]
                   + f_3 * pc_z[k] * sni_83[k];

        t_108[k] = f_19 * smi_84[k]
                   + f_1 * snh0_63[k]
                   - f_2 * snh1_63[k]
                   + f_3 * pc_x[k] * sni_84[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, pc_z, smi_28, smi_30, smi_87, \
                         snh0_66, snh1_66, sni_84, sni_86, sni_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_14 * smi_28[k]
                   + f_3 * pc_y[k] * sni_84[k];

        t_110[k] = f_3 * pc_z[k] * sni_84[k];

        t_111[k] = f_19 * smi_87[k]
                   + f_4 * snh0_66[k]
                   - f_5 * snh1_66[k]
                   + f_3 * pc_x[k] * sni_87[k];

        t_112[k] = f_14 * smi_30[k]
                   + f_3 * pc_y[k] * sni_86[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pc_x, pc_z, smi_89, smi_90, snh0_68, snh0_69, \
                         snh1_68, snh1_69, sni_87, sni_89, sni_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_19 * smi_89[k]
                   + f_4 * snh0_68[k]
                   - f_5 * snh1_68[k]
                   + f_3 * pc_x[k] * sni_89[k];

        t_114[k] = f_19 * smi_90[k]
                   + f_6 * snh0_69[k]
                   - f_7 * snh1_69[k]
                   + f_3 * pc_x[k] * sni_90[k];

        t_115[k] = f_3 * pc_z[k] * sni_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pc_x, pc_y, smi_33, smi_93, smi_94, snh0_72, \
                         snh0_73, snh1_72, snh1_73, sni_89, sni_93, \
                         sni_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_14 * smi_33[k]
                   + f_3 * pc_y[k] * sni_89[k];

        t_117[k] = f_19 * smi_93[k]
                   + f_6 * snh0_72[k]
                   - f_7 * snh1_72[k]
                   + f_3 * pc_x[k] * sni_93[k];

        t_118[k] = f_19 * smi_94[k]
                   + f_8 * snh0_73[k]
                   - f_9 * snh1_73[k]
                   + f_3 * pc_x[k] * sni_94[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, smi_37, smi_96, snh0_75, \
                         snh1_75, sni_90, sni_93, sni_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * pc_z[k] * sni_90[k];

        t_120[k] = f_19 * smi_96[k]
                   + f_8 * snh0_75[k]
                   - f_9 * snh1_75[k]
                   + f_3 * pc_x[k] * sni_96[k];

        t_121[k] = f_14 * smi_37[k]
                   + f_3 * pc_y[k] * sni_93[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;

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

    const auto *smk0_39 = buffer.data(smk0 + 39);
    const auto *smk0_42 = buffer.data(smk0 + 42);
    const auto *smk0_46 = buffer.data(smk0 + 46);
    const auto *smk0_51 = buffer.data(smk0 + 51);
    const auto *smk0_64 = buffer.data(smk0 + 64);
    const auto *smk0_72 = buffer.data(smk0 + 72);
    const auto *smk0_77 = buffer.data(smk0 + 77);
    const auto *smk0_81 = buffer.data(smk0 + 81);
    const auto *smk0_84 = buffer.data(smk0 + 84);
    const auto *smk0_86 = buffer.data(smk0 + 86);
    const auto *smk0_89 = buffer.data(smk0 + 89);
    const auto *smk0_90 = buffer.data(smk0 + 90);
    const auto *smk0_92 = buffer.data(smk0 + 92);
    const auto *smk0_107 = buffer.data(smk0 + 107);

    const auto *smi_28 = buffer.data(smi + 28);
    const auto *smi_31 = buffer.data(smi + 31);
    const auto *smi_34 = buffer.data(smi + 34);
    const auto *smi_38 = buffer.data(smi + 38);
    const auto *smi_42 = buffer.data(smi + 42);
    const auto *smi_49 = buffer.data(smi + 49);
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
    const auto *smi_64 = buffer.data(smi + 64);
    const auto *smi_65 = buffer.data(smi + 65);
    const auto *smi_66 = buffer.data(smi + 66);
    const auto *smi_68 = buffer.data(smi + 68);
    const auto *smi_69 = buffer.data(smi + 69);
    const auto *smi_70 = buffer.data(smi + 70);
    const auto *smi_77 = buffer.data(smi + 77);
    const auto *smi_79 = buffer.data(smi + 79);
    const auto *smi_80 = buffer.data(smi + 80);
    const auto *smi_81 = buffer.data(smi + 81);
    const auto *smi_82 = buffer.data(smi + 82);
    const auto *smi_83 = buffer.data(smi + 83);
    const auto *smi_84 = buffer.data(smi + 84);
    const auto *smi_86 = buffer.data(smi + 86);
    const auto *smi_89 = buffer.data(smi + 89);
    const auto *smi_93 = buffer.data(smi + 93);
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
    const auto *smi_133 = buffer.data(smi + 133);
    const auto *smi_134 = buffer.data(smi + 134);
    const auto *smi_135 = buffer.data(smi + 135);
    const auto *smi_136 = buffer.data(smi + 136);
    const auto *smi_137 = buffer.data(smi + 137);
    const auto *smi_138 = buffer.data(smi + 138);
    const auto *smi_139 = buffer.data(smi + 139);
    const auto *smi_140 = buffer.data(smi + 140);
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

    const auto *smk1_39 = buffer.data(smk1 + 39);
    const auto *smk1_42 = buffer.data(smk1 + 42);
    const auto *smk1_46 = buffer.data(smk1 + 46);
    const auto *smk1_51 = buffer.data(smk1 + 51);
    const auto *smk1_64 = buffer.data(smk1 + 64);
    const auto *smk1_72 = buffer.data(smk1 + 72);
    const auto *smk1_77 = buffer.data(smk1 + 77);
    const auto *smk1_81 = buffer.data(smk1 + 81);
    const auto *smk1_84 = buffer.data(smk1 + 84);
    const auto *smk1_86 = buffer.data(smk1 + 86);
    const auto *smk1_89 = buffer.data(smk1 + 89);
    const auto *smk1_90 = buffer.data(smk1 + 90);
    const auto *smk1_92 = buffer.data(smk1 + 92);
    const auto *smk1_107 = buffer.data(smk1 + 107);

    const auto *snh0_77 = buffer.data(snh0 + 77);
    const auto *snh0_78 = buffer.data(snh0 + 78);
    const auto *snh0_80 = buffer.data(snh0 + 80);
    const auto *snh0_81 = buffer.data(snh0 + 81);
    const auto *snh0_82 = buffer.data(snh0 + 82);
    const auto *snh0_83 = buffer.data(snh0 + 83);
    const auto *snh0_101 = buffer.data(snh0 + 101);
    const auto *snh0_102 = buffer.data(snh0 + 102);
    const auto *snh0_103 = buffer.data(snh0 + 103);
    const auto *snh0_104 = buffer.data(snh0 + 104);
    const auto *snh0_105 = buffer.data(snh0 + 105);
    const auto *snh0_108 = buffer.data(snh0 + 108);
    const auto *snh0_110 = buffer.data(snh0 + 110);
    const auto *snh0_111 = buffer.data(snh0 + 111);
    const auto *snh0_114 = buffer.data(snh0 + 114);
    const auto *snh0_115 = buffer.data(snh0 + 115);
    const auto *snh0_117 = buffer.data(snh0 + 117);
    const auto *snh0_119 = buffer.data(snh0 + 119);
    const auto *snh0_120 = buffer.data(snh0 + 120);
    const auto *snh0_122 = buffer.data(snh0 + 122);
    const auto *snh0_123 = buffer.data(snh0 + 123);
    const auto *snh0_124 = buffer.data(snh0 + 124);
    const auto *snh0_125 = buffer.data(snh0 + 125);
    const auto *snh0_126 = buffer.data(snh0 + 126);
    const auto *snh0_129 = buffer.data(snh0 + 129);
    const auto *snh0_131 = buffer.data(snh0 + 131);
    const auto *snh0_132 = buffer.data(snh0 + 132);
    const auto *snh0_135 = buffer.data(snh0 + 135);
    const auto *snh0_136 = buffer.data(snh0 + 136);
    const auto *snh0_138 = buffer.data(snh0 + 138);
    const auto *snh0_140 = buffer.data(snh0 + 140);
    const auto *snh0_141 = buffer.data(snh0 + 141);
    const auto *snh0_143 = buffer.data(snh0 + 143);
    const auto *snh0_144 = buffer.data(snh0 + 144);

    const auto *snh1_77 = buffer.data(snh1 + 77);
    const auto *snh1_78 = buffer.data(snh1 + 78);
    const auto *snh1_80 = buffer.data(snh1 + 80);
    const auto *snh1_81 = buffer.data(snh1 + 81);
    const auto *snh1_82 = buffer.data(snh1 + 82);
    const auto *snh1_83 = buffer.data(snh1 + 83);
    const auto *snh1_101 = buffer.data(snh1 + 101);
    const auto *snh1_102 = buffer.data(snh1 + 102);
    const auto *snh1_103 = buffer.data(snh1 + 103);
    const auto *snh1_104 = buffer.data(snh1 + 104);
    const auto *snh1_105 = buffer.data(snh1 + 105);
    const auto *snh1_108 = buffer.data(snh1 + 108);
    const auto *snh1_110 = buffer.data(snh1 + 110);
    const auto *snh1_111 = buffer.data(snh1 + 111);
    const auto *snh1_114 = buffer.data(snh1 + 114);
    const auto *snh1_115 = buffer.data(snh1 + 115);
    const auto *snh1_117 = buffer.data(snh1 + 117);
    const auto *snh1_119 = buffer.data(snh1 + 119);
    const auto *snh1_120 = buffer.data(snh1 + 120);
    const auto *snh1_122 = buffer.data(snh1 + 122);
    const auto *snh1_123 = buffer.data(snh1 + 123);
    const auto *snh1_124 = buffer.data(snh1 + 124);
    const auto *snh1_125 = buffer.data(snh1 + 125);
    const auto *snh1_126 = buffer.data(snh1 + 126);
    const auto *snh1_129 = buffer.data(snh1 + 129);
    const auto *snh1_131 = buffer.data(snh1 + 131);
    const auto *snh1_132 = buffer.data(snh1 + 132);
    const auto *snh1_135 = buffer.data(snh1 + 135);
    const auto *snh1_136 = buffer.data(snh1 + 136);
    const auto *snh1_138 = buffer.data(snh1 + 138);
    const auto *snh1_140 = buffer.data(snh1 + 140);
    const auto *snh1_141 = buffer.data(snh1 + 141);
    const auto *snh1_143 = buffer.data(snh1 + 143);
    const auto *snh1_144 = buffer.data(snh1 + 144);

    const auto *sni_94 = buffer.data(sni + 94);
    const auto *sni_98 = buffer.data(sni + 98);
    const auto *sni_99 = buffer.data(sni + 99);
    const auto *sni_101 = buffer.data(sni + 101);
    const auto *sni_102 = buffer.data(sni + 102);
    const auto *sni_104 = buffer.data(sni + 104);
    const auto *sni_105 = buffer.data(sni + 105);
    const auto *sni_106 = buffer.data(sni + 106);
    const auto *sni_107 = buffer.data(sni + 107);
    const auto *sni_108 = buffer.data(sni + 108);
    const auto *sni_109 = buffer.data(sni + 109);
    const auto *sni_110 = buffer.data(sni + 110);
    const auto *sni_111 = buffer.data(sni + 111);
    const auto *sni_112 = buffer.data(sni + 112);
    const auto *sni_114 = buffer.data(sni + 114);
    const auto *sni_115 = buffer.data(sni + 115);
    const auto *sni_117 = buffer.data(sni + 117);
    const auto *sni_118 = buffer.data(sni + 118);
    const auto *sni_121 = buffer.data(sni + 121);
    const auto *sni_122 = buffer.data(sni + 122);
    const auto *sni_126 = buffer.data(sni + 126);
    const auto *sni_133 = buffer.data(sni + 133);
    const auto *sni_134 = buffer.data(sni + 134);
    const auto *sni_135 = buffer.data(sni + 135);
    const auto *sni_136 = buffer.data(sni + 136);
    const auto *sni_137 = buffer.data(sni + 137);
    const auto *sni_138 = buffer.data(sni + 138);
    const auto *sni_139 = buffer.data(sni + 139);
    const auto *sni_140 = buffer.data(sni + 140);
    const auto *sni_142 = buffer.data(sni + 142);
    const auto *sni_143 = buffer.data(sni + 143);
    const auto *sni_145 = buffer.data(sni + 145);
    const auto *sni_146 = buffer.data(sni + 146);
    const auto *sni_149 = buffer.data(sni + 149);
    const auto *sni_150 = buffer.data(sni + 150);
    const auto *sni_152 = buffer.data(sni + 152);
    const auto *sni_154 = buffer.data(sni + 154);
    const auto *sni_155 = buffer.data(sni + 155);
    const auto *sni_157 = buffer.data(sni + 157);
    const auto *sni_158 = buffer.data(sni + 158);
    const auto *sni_160 = buffer.data(sni + 160);
    const auto *sni_161 = buffer.data(sni + 161);
    const auto *sni_162 = buffer.data(sni + 162);
    const auto *sni_163 = buffer.data(sni + 163);
    const auto *sni_164 = buffer.data(sni + 164);
    const auto *sni_165 = buffer.data(sni + 165);
    const auto *sni_166 = buffer.data(sni + 166);
    const auto *sni_167 = buffer.data(sni + 167);
    const auto *sni_168 = buffer.data(sni + 168);
    const auto *sni_170 = buffer.data(sni + 170);
    const auto *sni_171 = buffer.data(sni + 171);
    const auto *sni_173 = buffer.data(sni + 173);
    const auto *sni_174 = buffer.data(sni + 174);
    const auto *sni_177 = buffer.data(sni + 177);
    const auto *sni_178 = buffer.data(sni + 178);
    const auto *sni_180 = buffer.data(sni + 180);
    const auto *sni_182 = buffer.data(sni + 182);
    const auto *sni_183 = buffer.data(sni + 183);
    const auto *sni_185 = buffer.data(sni + 185);
    const auto *sni_186 = buffer.data(sni + 186);

#pragma omp simd aligned(t_122, t_123, t_124, pc_x, pc_z, smi_98, smi_99, snh0_77, snh0_78, \
                         snh1_77, snh1_78, sni_94, sni_98, sni_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_19 * smi_98[k]
                   + f_8 * snh0_77[k]
                   - f_9 * snh1_77[k]
                   + f_3 * pc_x[k] * sni_98[k];

        t_123[k] = f_19 * smi_99[k]
                   + f_10 * snh0_78[k]
                   - f_11 * snh1_78[k]
                   + f_3 * pc_x[k] * sni_99[k];

        t_124[k] = f_3 * pc_z[k] * sni_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_y, smi_42, smi_101, smi_102, snh0_80, \
                         snh0_81, snh1_80, snh1_81, sni_98, sni_101, \
                         sni_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_19 * smi_101[k]
                   + f_10 * snh0_80[k]
                   - f_11 * snh1_80[k]
                   + f_3 * pc_x[k] * sni_101[k];

        t_126[k] = f_19 * smi_102[k]
                   + f_10 * snh0_81[k]
                   - f_11 * snh1_81[k]
                   + f_3 * pc_x[k] * sni_102[k];

        t_127[k] = f_14 * smi_42[k]
                   + f_3 * pc_y[k] * sni_98[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, smi_104, smi_105, smi_106, smi_107, \
                         snh0_83, snh1_83, sni_104, sni_105, sni_106, \
                         sni_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_19 * smi_104[k]
                   + f_10 * snh0_83[k]
                   - f_11 * snh1_83[k]
                   + f_3 * pc_x[k] * sni_104[k];

        t_129[k] = f_19 * smi_105[k]
                   + f_3 * pc_x[k] * sni_105[k];

        t_130[k] = f_19 * smi_106[k]
                   + f_3 * pc_x[k] * sni_106[k];

        t_131[k] = f_19 * smi_107[k]
                   + f_3 * pc_x[k] * sni_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, smi_108, smi_109, smi_110, smi_111, \
                         sni_108, sni_109, sni_110, sni_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_19 * smi_108[k]
                   + f_3 * pc_x[k] * sni_108[k];

        t_133[k] = f_19 * smi_109[k]
                   + f_3 * pc_x[k] * sni_109[k];

        t_134[k] = f_19 * smi_110[k]
                   + f_3 * pc_x[k] * sni_110[k];

        t_135[k] = f_19 * smi_111[k]
                   + f_3 * pc_x[k] * sni_111[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, pc_z, smi_49, smi_51, snh0_78, snh0_80, \
                         snh1_78, snh1_80, sni_105, sni_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_14 * smi_49[k]
                   + f_1 * snh0_78[k]
                   - f_2 * snh1_78[k]
                   + f_3 * pc_y[k] * sni_105[k];

        t_137[k] = f_3 * pc_z[k] * sni_105[k];

        t_138[k] = f_14 * smi_51[k]
                   + f_4 * snh0_80[k]
                   - f_5 * snh1_80[k]
                   + f_3 * pc_y[k] * sni_107[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pc_y, smi_52, smi_53, smi_54, snh0_81, snh0_82, \
                         snh0_83, snh1_81, snh1_82, snh1_83, sni_108, sni_109, \
                         sni_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_14 * smi_52[k]
                   + f_6 * snh0_81[k]
                   - f_7 * snh1_81[k]
                   + f_3 * pc_y[k] * sni_108[k];

        t_140[k] = f_14 * smi_53[k]
                   + f_8 * snh0_82[k]
                   - f_9 * snh1_82[k]
                   + f_3 * pc_y[k] * sni_109[k];

        t_141[k] = f_14 * smi_54[k]
                   + f_10 * snh0_83[k]
                   - f_11 * snh1_83[k]
                   + f_3 * pc_y[k] * sni_110[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_y, pc_y, pc_z, smk0_72, smi_55, \
                         smi_56, smk1_72, snh0_83, snh1_83, sni_111, \
                         sni_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * smi_55[k]
                   + f_3 * pc_y[k] * sni_111[k];

        t_143[k] = f_1 * snh0_83[k]
                   - f_2 * snh1_83[k]
                   + f_3 * pc_z[k] * sni_111[k];

        t_144[k] = pb_y[k] * smk0_72[k]
                   - f_12 * pc_y[k] * smk1_72[k];

        t_145[k] = f_13 * smi_56[k]
                   + f_3 * pc_y[k] * sni_112[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_y, pb_z, pc_y, pc_z, smk0_39, smk0_77, \
                         smi_28, smi_58, smk1_39, smk1_77, sni_112, \
                         sni_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_13 * smi_28[k]
                   + f_3 * pc_z[k] * sni_112[k];

        t_147[k] = pb_z[k] * smk0_39[k]
                   - f_12 * pc_z[k] * smk1_39[k];

        t_148[k] = f_13 * smi_58[k]
                   + f_3 * pc_y[k] * sni_114[k];

        t_149[k] = pb_y[k] * smk0_77[k]
                   - f_12 * pc_y[k] * smk1_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, pb_z, pc_y, pc_z, smk0_42, smk0_81, \
                         smi_31, smi_61, smk1_42, smk1_81, sni_115, \
                         sni_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * smk0_42[k]
                   - f_12 * pc_z[k] * smk1_42[k];

        t_151[k] = f_13 * smi_31[k]
                   + f_3 * pc_z[k] * sni_115[k];

        t_152[k] = f_13 * smi_61[k]
                   + f_3 * pc_y[k] * sni_117[k];

        t_153[k] = pb_y[k] * smk0_81[k]
                   - f_12 * pc_y[k] * smk1_81[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pb_z, pc_y, pc_z, smk0_46, smk0_84, \
                         smi_34, smi_64, smk1_46, smk1_84, sni_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * smk0_46[k]
                   - f_12 * pc_z[k] * smk1_46[k];

        t_155[k] = f_13 * smi_34[k]
                   + f_3 * pc_z[k] * sni_118[k];

        t_156[k] = pb_y[k] * smk0_84[k]
                   + f_14 * smi_64[k]
                   - f_12 * pc_y[k] * smk1_84[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pb_y, pb_z, pc_y, pc_z, smk0_51, smk0_86, \
                         smi_38, smi_65, smk1_51, smk1_86, sni_121, \
                         sni_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_13 * smi_65[k]
                   + f_3 * pc_y[k] * sni_121[k];

        t_158[k] = pb_y[k] * smk0_86[k]
                   - f_12 * pc_y[k] * smk1_86[k];

        t_159[k] = pb_z[k] * smk0_51[k]
                   - f_12 * pc_z[k] * smk1_51[k];

        t_160[k] = f_13 * smi_38[k]
                   + f_3 * pc_z[k] * sni_122[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pc_y, smk0_89, smk0_90, smk0_92, \
                         smi_68, smi_69, smi_70, smk1_89, smk1_90, smk1_92, \
                         sni_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * smk0_89[k]
                   + f_15 * smi_68[k]
                   - f_12 * pc_y[k] * smk1_89[k];

        t_162[k] = pb_y[k] * smk0_90[k]
                   + f_14 * smi_69[k]
                   - f_12 * pc_y[k] * smk1_90[k];

        t_163[k] = f_13 * smi_70[k]
                   + f_3 * pc_y[k] * sni_126[k];

        t_164[k] = pb_y[k] * smk0_92[k]
                   - f_12 * pc_y[k] * smk1_92[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, pc_x, smi_133, smi_134, smi_135, \
                         smi_136, smi_137, sni_133, sni_134, sni_135, sni_136, \
                         sni_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_19 * smi_133[k]
                   + f_3 * pc_x[k] * sni_133[k];

        t_166[k] = f_19 * smi_134[k]
                   + f_3 * pc_x[k] * sni_134[k];

        t_167[k] = f_19 * smi_135[k]
                   + f_3 * pc_x[k] * sni_135[k];

        t_168[k] = f_19 * smi_136[k]
                   + f_3 * pc_x[k] * sni_136[k];

        t_169[k] = f_19 * smi_137[k]
                   + f_3 * pc_x[k] * sni_137[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_z, pc_x, pc_z, smk0_64, smi_49, \
                         smi_138, smi_139, smk1_64, sni_133, sni_138, \
                         sni_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_19 * smi_138[k]
                   + f_3 * pc_x[k] * sni_138[k];

        t_171[k] = f_19 * smi_139[k]
                   + f_3 * pc_x[k] * sni_139[k];

        t_172[k] = pb_z[k] * smk0_64[k]
                   - f_12 * pc_z[k] * smk1_64[k];

        t_173[k] = f_13 * smi_49[k]
                   + f_3 * pc_z[k] * sni_133[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, smi_79, smi_80, smi_81, snh0_101, \
                         snh0_102, snh0_103, snh1_101, snh1_102, snh1_103, sni_135, sni_136, \
                         sni_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * smi_79[k]
                   + f_4 * snh0_101[k]
                   - f_5 * snh1_101[k]
                   + f_3 * pc_y[k] * sni_135[k];

        t_175[k] = f_13 * smi_80[k]
                   + f_6 * snh0_102[k]
                   - f_7 * snh1_102[k]
                   + f_3 * pc_y[k] * sni_136[k];

        t_176[k] = f_13 * smi_81[k]
                   + f_8 * snh0_103[k]
                   - f_9 * snh1_103[k]
                   + f_3 * pc_y[k] * sni_137[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_y, pc_y, smk0_107, smi_82, smi_83, smk1_107, \
                         snh0_104, snh1_104, sni_138, sni_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_13 * smi_82[k]
                   + f_10 * snh0_104[k]
                   - f_11 * snh1_104[k]
                   + f_3 * pc_y[k] * sni_138[k];

        t_178[k] = f_13 * smi_83[k]
                   + f_3 * pc_y[k] * sni_139[k];

        t_179[k] = pb_y[k] * smk0_107[k]
                   - f_12 * pc_y[k] * smk1_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, smi_56, smi_140, \
                         smi_143, snh0_105, snh0_108, snh1_105, snh1_108, sni_140, \
                         sni_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_19 * smi_140[k]
                   + f_1 * snh0_105[k]
                   - f_2 * snh1_105[k]
                   + f_3 * pc_x[k] * sni_140[k];

        t_181[k] = f_3 * pc_y[k] * sni_140[k];

        t_182[k] = f_14 * smi_56[k]
                   + f_3 * pc_z[k] * sni_140[k];

        t_183[k] = f_19 * smi_143[k]
                   + f_4 * snh0_108[k]
                   - f_5 * snh1_108[k]
                   + f_3 * pc_x[k] * sni_143[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, smi_145, smi_146, snh0_110, \
                         snh0_111, snh1_110, snh1_111, sni_142, sni_145, \
                         sni_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_y[k] * sni_142[k];

        t_185[k] = f_19 * smi_145[k]
                   + f_4 * snh0_110[k]
                   - f_5 * snh1_110[k]
                   + f_3 * pc_x[k] * sni_145[k];

        t_186[k] = f_19 * smi_146[k]
                   + f_6 * snh0_111[k]
                   - f_7 * snh1_111[k]
                   + f_3 * pc_x[k] * sni_146[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pc_x, pc_y, pc_z, smi_59, smi_149, snh0_114, \
                         snh1_114, sni_143, sni_145, sni_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_14 * smi_59[k]
                   + f_3 * pc_z[k] * sni_143[k];

        t_188[k] = f_3 * pc_y[k] * sni_145[k];

        t_189[k] = f_19 * smi_149[k]
                   + f_6 * snh0_114[k]
                   - f_7 * snh1_114[k]
                   + f_3 * pc_x[k] * sni_149[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_z, smi_62, smi_150, smi_152, snh0_115, \
                         snh0_117, snh1_115, snh1_117, sni_146, sni_150, \
                         sni_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_19 * smi_150[k]
                   + f_8 * snh0_115[k]
                   - f_9 * snh1_115[k]
                   + f_3 * pc_x[k] * sni_150[k];

        t_191[k] = f_14 * smi_62[k]
                   + f_3 * pc_z[k] * sni_146[k];

        t_192[k] = f_19 * smi_152[k]
                   + f_8 * snh0_117[k]
                   - f_9 * snh1_117[k]
                   + f_3 * pc_x[k] * sni_152[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, smi_154, smi_155, snh0_119, \
                         snh0_120, snh1_119, snh1_120, sni_149, sni_154, \
                         sni_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * sni_149[k];

        t_194[k] = f_19 * smi_154[k]
                   + f_8 * snh0_119[k]
                   - f_9 * snh1_119[k]
                   + f_3 * pc_x[k] * sni_154[k];

        t_195[k] = f_19 * smi_155[k]
                   + f_10 * snh0_120[k]
                   - f_11 * snh1_120[k]
                   + f_3 * pc_x[k] * sni_155[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pc_x, pc_z, smi_66, smi_157, smi_158, snh0_122, \
                         snh0_123, snh1_122, snh1_123, sni_150, sni_157, \
                         sni_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_14 * smi_66[k]
                   + f_3 * pc_z[k] * sni_150[k];

        t_197[k] = f_19 * smi_157[k]
                   + f_10 * snh0_122[k]
                   - f_11 * snh1_122[k]
                   + f_3 * pc_x[k] * sni_157[k];

        t_198[k] = f_19 * smi_158[k]
                   + f_10 * snh0_123[k]
                   - f_11 * snh1_123[k]
                   + f_3 * pc_x[k] * sni_158[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pc_x, pc_y, smi_160, smi_161, smi_162, \
                         snh0_125, snh1_125, sni_154, sni_160, sni_161, \
                         sni_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_3 * pc_y[k] * sni_154[k];

        t_200[k] = f_19 * smi_160[k]
                   + f_10 * snh0_125[k]
                   - f_11 * snh1_125[k]
                   + f_3 * pc_x[k] * sni_160[k];

        t_201[k] = f_19 * smi_161[k]
                   + f_3 * pc_x[k] * sni_161[k];

        t_202[k] = f_19 * smi_162[k]
                   + f_3 * pc_x[k] * sni_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pc_x, smi_163, smi_164, smi_165, \
                         smi_166, smi_167, sni_163, sni_164, sni_165, sni_166, \
                         sni_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_19 * smi_163[k]
                   + f_3 * pc_x[k] * sni_163[k];

        t_204[k] = f_19 * smi_164[k]
                   + f_3 * pc_x[k] * sni_164[k];

        t_205[k] = f_19 * smi_165[k]
                   + f_3 * pc_x[k] * sni_165[k];

        t_206[k] = f_19 * smi_166[k]
                   + f_3 * pc_x[k] * sni_166[k];

        t_207[k] = f_19 * smi_167[k]
                   + f_3 * pc_x[k] * sni_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_y, pc_z, smi_77, snh0_120, snh0_122, \
                         snh0_123, snh1_120, snh1_122, snh1_123, sni_161, sni_163, \
                         sni_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * snh0_120[k]
                   - f_2 * snh1_120[k]
                   + f_3 * pc_y[k] * sni_161[k];

        t_209[k] = f_14 * smi_77[k]
                   + f_3 * pc_z[k] * sni_161[k];

        t_210[k] = f_4 * snh0_122[k]
                   - f_5 * snh1_122[k]
                   + f_3 * pc_y[k] * sni_163[k];

        t_211[k] = f_6 * snh0_123[k]
                   - f_7 * snh1_123[k]
                   + f_3 * pc_y[k] * sni_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_y, pc_z, smi_83, snh0_124, snh0_125, \
                         snh1_124, snh1_125, sni_165, sni_166, \
                         sni_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * snh0_124[k]
                   - f_9 * snh1_124[k]
                   + f_3 * pc_y[k] * sni_165[k];

        t_213[k] = f_10 * snh0_125[k]
                   - f_11 * snh1_125[k]
                   + f_3 * pc_y[k] * sni_166[k];

        t_214[k] = f_3 * pc_y[k] * sni_167[k];

        t_215[k] = f_14 * smi_83[k]
                   + f_1 * snh0_125[k]
                   - f_2 * snh1_125[k]
                   + f_3 * pc_z[k] * sni_167[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, pc_z, smi_84, smi_168, \
                         smi_171, snh0_126, snh0_129, snh1_126, snh1_129, sni_168, \
                         sni_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_20 * smi_168[k]
                   + f_1 * snh0_126[k]
                   - f_2 * snh1_126[k]
                   + f_3 * pc_x[k] * sni_168[k];

        t_217[k] = f_15 * smi_84[k]
                   + f_3 * pc_y[k] * sni_168[k];

        t_218[k] = f_3 * pc_z[k] * sni_168[k];

        t_219[k] = f_20 * smi_171[k]
                   + f_4 * snh0_129[k]
                   - f_5 * snh1_129[k]
                   + f_3 * pc_x[k] * sni_171[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_x, pc_y, smi_86, smi_173, smi_174, snh0_131, \
                         snh0_132, snh1_131, snh1_132, sni_170, sni_173, \
                         sni_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * smi_86[k]
                   + f_3 * pc_y[k] * sni_170[k];

        t_221[k] = f_20 * smi_173[k]
                   + f_4 * snh0_131[k]
                   - f_5 * snh1_131[k]
                   + f_3 * pc_x[k] * sni_173[k];

        t_222[k] = f_20 * smi_174[k]
                   + f_6 * snh0_132[k]
                   - f_7 * snh1_132[k]
                   + f_3 * pc_x[k] * sni_174[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pc_x, pc_y, pc_z, smi_89, smi_177, snh0_135, \
                         snh1_135, sni_171, sni_173, sni_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_z[k] * sni_171[k];

        t_224[k] = f_15 * smi_89[k]
                   + f_3 * pc_y[k] * sni_173[k];

        t_225[k] = f_20 * smi_177[k]
                   + f_6 * snh0_135[k]
                   - f_7 * snh1_135[k]
                   + f_3 * pc_x[k] * sni_177[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_x, pc_z, smi_178, smi_180, snh0_136, \
                         snh0_138, snh1_136, snh1_138, sni_174, sni_178, \
                         sni_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_20 * smi_178[k]
                   + f_8 * snh0_136[k]
                   - f_9 * snh1_136[k]
                   + f_3 * pc_x[k] * sni_178[k];

        t_227[k] = f_3 * pc_z[k] * sni_174[k];

        t_228[k] = f_20 * smi_180[k]
                   + f_8 * snh0_138[k]
                   - f_9 * snh1_138[k]
                   + f_3 * pc_x[k] * sni_180[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, smi_93, smi_182, smi_183, snh0_140, \
                         snh0_141, snh1_140, snh1_141, sni_177, sni_182, \
                         sni_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_15 * smi_93[k]
                   + f_3 * pc_y[k] * sni_177[k];

        t_230[k] = f_20 * smi_182[k]
                   + f_8 * snh0_140[k]
                   - f_9 * snh1_140[k]
                   + f_3 * pc_x[k] * sni_182[k];

        t_231[k] = f_20 * smi_183[k]
                   + f_10 * snh0_141[k]
                   - f_11 * snh1_141[k]
                   + f_3 * pc_x[k] * sni_183[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_z, smi_185, smi_186, snh0_143, \
                         snh0_144, snh1_143, snh1_144, sni_178, sni_185, \
                         sni_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * pc_z[k] * sni_178[k];

        t_233[k] = f_20 * smi_185[k]
                   + f_10 * snh0_143[k]
                   - f_11 * snh1_143[k]
                   + f_3 * pc_x[k] * sni_185[k];

        t_234[k] = f_20 * smi_186[k]
                   + f_10 * snh0_144[k]
                   - f_11 * snh1_144[k]
                   + f_3 * pc_x[k] * sni_186[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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
    const auto f_20 = 3.5 / q;

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

    const auto *smk0_108 = buffer.data(smk0 + 108);
    const auto *smk0_111 = buffer.data(smk0 + 111);
    const auto *smk0_114 = buffer.data(smk0 + 114);
    const auto *smk0_118 = buffer.data(smk0 + 118);
    const auto *smk0_120 = buffer.data(smk0 + 120);
    const auto *smk0_123 = buffer.data(smk0 + 123);
    const auto *smk0_125 = buffer.data(smk0 + 125);
    const auto *smk0_126 = buffer.data(smk0 + 126);
    const auto *smk0_136 = buffer.data(smk0 + 136);
    const auto *smk0_180 = buffer.data(smk0 + 180);
    const auto *smk0_183 = buffer.data(smk0 + 183);
    const auto *smk0_185 = buffer.data(smk0 + 185);
    const auto *smk0_186 = buffer.data(smk0 + 186);
    const auto *smk0_189 = buffer.data(smk0 + 189);
    const auto *smk0_190 = buffer.data(smk0 + 190);
    const auto *smk0_192 = buffer.data(smk0 + 192);
    const auto *smk0_194 = buffer.data(smk0 + 194);
    const auto *smk0_195 = buffer.data(smk0 + 195);
    const auto *smk0_197 = buffer.data(smk0 + 197);
    const auto *smk0_198 = buffer.data(smk0 + 198);
    const auto *smk0_200 = buffer.data(smk0 + 200);
    const auto *smk0_215 = buffer.data(smk0 + 215);

    const auto *smi_84 = buffer.data(smi + 84);
    const auto *smi_87 = buffer.data(smi + 87);
    const auto *smi_90 = buffer.data(smi + 90);
    const auto *smi_91 = buffer.data(smi + 91);
    const auto *smi_94 = buffer.data(smi + 94);
    const auto *smi_95 = buffer.data(smi + 95);
    const auto *smi_96 = buffer.data(smi + 96);
    const auto *smi_98 = buffer.data(smi + 98);
    const auto *smi_105 = buffer.data(smi + 105);
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
    const auto *smi_135 = buffer.data(smi + 135);
    const auto *smi_136 = buffer.data(smi + 136);
    const auto *smi_137 = buffer.data(smi + 137);
    const auto *smi_138 = buffer.data(smi + 138);
    const auto *smi_139 = buffer.data(smi + 139);
    const auto *smi_140 = buffer.data(smi + 140);
    const auto *smi_141 = buffer.data(smi + 141);
    const auto *smi_142 = buffer.data(smi + 142);
    const auto *smi_143 = buffer.data(smi + 143);
    const auto *smi_145 = buffer.data(smi + 145);
    const auto *smi_146 = buffer.data(smi + 146);
    const auto *smi_148 = buffer.data(smi + 148);
    const auto *smi_149 = buffer.data(smi + 149);
    const auto *smi_150 = buffer.data(smi + 150);
    const auto *smi_152 = buffer.data(smi + 152);
    const auto *smi_153 = buffer.data(smi + 153);
    const auto *smi_154 = buffer.data(smi + 154);
    const auto *smi_161 = buffer.data(smi + 161);
    const auto *smi_163 = buffer.data(smi + 163);
    const auto *smi_164 = buffer.data(smi + 164);
    const auto *smi_165 = buffer.data(smi + 165);
    const auto *smi_166 = buffer.data(smi + 166);
    const auto *smi_167 = buffer.data(smi + 167);
    const auto *smi_188 = buffer.data(smi + 188);
    const auto *smi_189 = buffer.data(smi + 189);
    const auto *smi_190 = buffer.data(smi + 190);
    const auto *smi_191 = buffer.data(smi + 191);
    const auto *smi_192 = buffer.data(smi + 192);
    const auto *smi_193 = buffer.data(smi + 193);
    const auto *smi_194 = buffer.data(smi + 194);
    const auto *smi_195 = buffer.data(smi + 195);
    const auto *smi_201 = buffer.data(smi + 201);
    const auto *smi_205 = buffer.data(smi + 205);
    const auto *smi_210 = buffer.data(smi + 210);
    const auto *smi_216 = buffer.data(smi + 216);
    const auto *smi_217 = buffer.data(smi + 217);
    const auto *smi_218 = buffer.data(smi + 218);
    const auto *smi_219 = buffer.data(smi + 219);
    const auto *smi_220 = buffer.data(smi + 220);
    const auto *smi_221 = buffer.data(smi + 221);
    const auto *smi_222 = buffer.data(smi + 222);
    const auto *smi_223 = buffer.data(smi + 223);
    const auto *smi_245 = buffer.data(smi + 245);
    const auto *smi_246 = buffer.data(smi + 246);
    const auto *smi_247 = buffer.data(smi + 247);
    const auto *smi_248 = buffer.data(smi + 248);
    const auto *smi_249 = buffer.data(smi + 249);
    const auto *smi_250 = buffer.data(smi + 250);
    const auto *smi_251 = buffer.data(smi + 251);
    const auto *smi_252 = buffer.data(smi + 252);
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

    const auto *smk1_108 = buffer.data(smk1 + 108);
    const auto *smk1_111 = buffer.data(smk1 + 111);
    const auto *smk1_114 = buffer.data(smk1 + 114);
    const auto *smk1_118 = buffer.data(smk1 + 118);
    const auto *smk1_120 = buffer.data(smk1 + 120);
    const auto *smk1_123 = buffer.data(smk1 + 123);
    const auto *smk1_125 = buffer.data(smk1 + 125);
    const auto *smk1_126 = buffer.data(smk1 + 126);
    const auto *smk1_136 = buffer.data(smk1 + 136);
    const auto *smk1_180 = buffer.data(smk1 + 180);
    const auto *smk1_183 = buffer.data(smk1 + 183);
    const auto *smk1_185 = buffer.data(smk1 + 185);
    const auto *smk1_186 = buffer.data(smk1 + 186);
    const auto *smk1_189 = buffer.data(smk1 + 189);
    const auto *smk1_190 = buffer.data(smk1 + 190);
    const auto *smk1_192 = buffer.data(smk1 + 192);
    const auto *smk1_194 = buffer.data(smk1 + 194);
    const auto *smk1_195 = buffer.data(smk1 + 195);
    const auto *smk1_197 = buffer.data(smk1 + 197);
    const auto *smk1_198 = buffer.data(smk1 + 198);
    const auto *smk1_200 = buffer.data(smk1 + 200);
    const auto *smk1_215 = buffer.data(smk1 + 215);

    const auto *snh0_141 = buffer.data(snh0 + 141);
    const auto *snh0_143 = buffer.data(snh0 + 143);
    const auto *snh0_144 = buffer.data(snh0 + 144);
    const auto *snh0_145 = buffer.data(snh0 + 145);
    const auto *snh0_146 = buffer.data(snh0 + 146);
    const auto *snh0_152 = buffer.data(snh0 + 152);
    const auto *snh0_156 = buffer.data(snh0 + 156);
    const auto *snh0_161 = buffer.data(snh0 + 161);
    const auto *snh0_164 = buffer.data(snh0 + 164);
    const auto *snh0_165 = buffer.data(snh0 + 165);
    const auto *snh0_166 = buffer.data(snh0 + 166);
    const auto *snh0_167 = buffer.data(snh0 + 167);
    const auto *snh0_183 = buffer.data(snh0 + 183);
    const auto *snh0_185 = buffer.data(snh0 + 185);
    const auto *snh0_186 = buffer.data(snh0 + 186);
    const auto *snh0_187 = buffer.data(snh0 + 187);
    const auto *snh0_188 = buffer.data(snh0 + 188);
    const auto *snh0_189 = buffer.data(snh0 + 189);
    const auto *snh0_192 = buffer.data(snh0 + 192);
    const auto *snh0_194 = buffer.data(snh0 + 194);
    const auto *snh0_195 = buffer.data(snh0 + 195);
    const auto *snh0_198 = buffer.data(snh0 + 198);
    const auto *snh0_199 = buffer.data(snh0 + 199);
    const auto *snh0_201 = buffer.data(snh0 + 201);
    const auto *snh0_203 = buffer.data(snh0 + 203);
    const auto *snh0_204 = buffer.data(snh0 + 204);
    const auto *snh0_206 = buffer.data(snh0 + 206);
    const auto *snh0_207 = buffer.data(snh0 + 207);
    const auto *snh0_209 = buffer.data(snh0 + 209);

    const auto *snh1_141 = buffer.data(snh1 + 141);
    const auto *snh1_143 = buffer.data(snh1 + 143);
    const auto *snh1_144 = buffer.data(snh1 + 144);
    const auto *snh1_145 = buffer.data(snh1 + 145);
    const auto *snh1_146 = buffer.data(snh1 + 146);
    const auto *snh1_152 = buffer.data(snh1 + 152);
    const auto *snh1_156 = buffer.data(snh1 + 156);
    const auto *snh1_161 = buffer.data(snh1 + 161);
    const auto *snh1_164 = buffer.data(snh1 + 164);
    const auto *snh1_165 = buffer.data(snh1 + 165);
    const auto *snh1_166 = buffer.data(snh1 + 166);
    const auto *snh1_167 = buffer.data(snh1 + 167);
    const auto *snh1_183 = buffer.data(snh1 + 183);
    const auto *snh1_185 = buffer.data(snh1 + 185);
    const auto *snh1_186 = buffer.data(snh1 + 186);
    const auto *snh1_187 = buffer.data(snh1 + 187);
    const auto *snh1_188 = buffer.data(snh1 + 188);
    const auto *snh1_189 = buffer.data(snh1 + 189);
    const auto *snh1_192 = buffer.data(snh1 + 192);
    const auto *snh1_194 = buffer.data(snh1 + 194);
    const auto *snh1_195 = buffer.data(snh1 + 195);
    const auto *snh1_198 = buffer.data(snh1 + 198);
    const auto *snh1_199 = buffer.data(snh1 + 199);
    const auto *snh1_201 = buffer.data(snh1 + 201);
    const auto *snh1_203 = buffer.data(snh1 + 203);
    const auto *snh1_204 = buffer.data(snh1 + 204);
    const auto *snh1_206 = buffer.data(snh1 + 206);
    const auto *snh1_207 = buffer.data(snh1 + 207);
    const auto *snh1_209 = buffer.data(snh1 + 209);

    const auto *sni_182 = buffer.data(sni + 182);
    const auto *sni_188 = buffer.data(sni + 188);
    const auto *sni_189 = buffer.data(sni + 189);
    const auto *sni_190 = buffer.data(sni + 190);
    const auto *sni_191 = buffer.data(sni + 191);
    const auto *sni_192 = buffer.data(sni + 192);
    const auto *sni_193 = buffer.data(sni + 193);
    const auto *sni_194 = buffer.data(sni + 194);
    const auto *sni_195 = buffer.data(sni + 195);
    const auto *sni_196 = buffer.data(sni + 196);
    const auto *sni_198 = buffer.data(sni + 198);
    const auto *sni_199 = buffer.data(sni + 199);
    const auto *sni_201 = buffer.data(sni + 201);
    const auto *sni_202 = buffer.data(sni + 202);
    const auto *sni_205 = buffer.data(sni + 205);
    const auto *sni_206 = buffer.data(sni + 206);
    const auto *sni_210 = buffer.data(sni + 210);
    const auto *sni_216 = buffer.data(sni + 216);
    const auto *sni_217 = buffer.data(sni + 217);
    const auto *sni_218 = buffer.data(sni + 218);
    const auto *sni_219 = buffer.data(sni + 219);
    const auto *sni_220 = buffer.data(sni + 220);
    const auto *sni_221 = buffer.data(sni + 221);
    const auto *sni_222 = buffer.data(sni + 222);
    const auto *sni_223 = buffer.data(sni + 223);
    const auto *sni_224 = buffer.data(sni + 224);
    const auto *sni_226 = buffer.data(sni + 226);
    const auto *sni_227 = buffer.data(sni + 227);
    const auto *sni_229 = buffer.data(sni + 229);
    const auto *sni_230 = buffer.data(sni + 230);
    const auto *sni_233 = buffer.data(sni + 233);
    const auto *sni_234 = buffer.data(sni + 234);
    const auto *sni_238 = buffer.data(sni + 238);
    const auto *sni_245 = buffer.data(sni + 245);
    const auto *sni_246 = buffer.data(sni + 246);
    const auto *sni_247 = buffer.data(sni + 247);
    const auto *sni_248 = buffer.data(sni + 248);
    const auto *sni_249 = buffer.data(sni + 249);
    const auto *sni_250 = buffer.data(sni + 250);
    const auto *sni_251 = buffer.data(sni + 251);
    const auto *sni_252 = buffer.data(sni + 252);
    const auto *sni_254 = buffer.data(sni + 254);
    const auto *sni_255 = buffer.data(sni + 255);
    const auto *sni_257 = buffer.data(sni + 257);
    const auto *sni_258 = buffer.data(sni + 258);
    const auto *sni_261 = buffer.data(sni + 261);
    const auto *sni_262 = buffer.data(sni + 262);
    const auto *sni_264 = buffer.data(sni + 264);
    const auto *sni_266 = buffer.data(sni + 266);
    const auto *sni_267 = buffer.data(sni + 267);
    const auto *sni_269 = buffer.data(sni + 269);
    const auto *sni_270 = buffer.data(sni + 270);
    const auto *sni_272 = buffer.data(sni + 272);
    const auto *sni_273 = buffer.data(sni + 273);
    const auto *sni_274 = buffer.data(sni + 274);

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, smi_98, smi_188, smi_189, \
                         smi_190, snh0_146, snh1_146, sni_182, sni_188, sni_189, \
                         sni_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_15 * smi_98[k]
                   + f_3 * pc_y[k] * sni_182[k];

        t_236[k] = f_20 * smi_188[k]
                   + f_10 * snh0_146[k]
                   - f_11 * snh1_146[k]
                   + f_3 * pc_x[k] * sni_188[k];

        t_237[k] = f_20 * smi_189[k]
                   + f_3 * pc_x[k] * sni_189[k];

        t_238[k] = f_20 * smi_190[k]
                   + f_3 * pc_x[k] * sni_190[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pc_x, smi_191, smi_192, smi_193, \
                         smi_194, smi_195, sni_191, sni_192, sni_193, sni_194, \
                         sni_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_20 * smi_191[k]
                   + f_3 * pc_x[k] * sni_191[k];

        t_240[k] = f_20 * smi_192[k]
                   + f_3 * pc_x[k] * sni_192[k];

        t_241[k] = f_20 * smi_193[k]
                   + f_3 * pc_x[k] * sni_193[k];

        t_242[k] = f_20 * smi_194[k]
                   + f_3 * pc_x[k] * sni_194[k];

        t_243[k] = f_20 * smi_195[k]
                   + f_3 * pc_x[k] * sni_195[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, pc_y, pc_z, smi_105, smi_107, snh0_141, \
                         snh0_143, snh1_141, snh1_143, sni_189, \
                         sni_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_15 * smi_105[k]
                   + f_1 * snh0_141[k]
                   - f_2 * snh1_141[k]
                   + f_3 * pc_y[k] * sni_189[k];

        t_245[k] = f_3 * pc_z[k] * sni_189[k];

        t_246[k] = f_15 * smi_107[k]
                   + f_4 * snh0_143[k]
                   - f_5 * snh1_143[k]
                   + f_3 * pc_y[k] * sni_191[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pc_y, smi_108, smi_109, smi_110, snh0_144, \
                         snh0_145, snh0_146, snh1_144, snh1_145, snh1_146, sni_192, sni_193, \
                         sni_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_15 * smi_108[k]
                   + f_6 * snh0_144[k]
                   - f_7 * snh1_144[k]
                   + f_3 * pc_y[k] * sni_192[k];

        t_248[k] = f_15 * smi_109[k]
                   + f_8 * snh0_145[k]
                   - f_9 * snh1_145[k]
                   + f_3 * pc_y[k] * sni_193[k];

        t_249[k] = f_15 * smi_110[k]
                   + f_10 * snh0_146[k]
                   - f_11 * snh1_146[k]
                   + f_3 * pc_y[k] * sni_194[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pb_z, pc_y, pc_z, smk0_108, smi_111, \
                         smi_112, smk1_108, snh0_146, snh1_146, sni_195, \
                         sni_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_15 * smi_111[k]
                   + f_3 * pc_y[k] * sni_195[k];

        t_251[k] = f_1 * snh0_146[k]
                   - f_2 * snh1_146[k]
                   + f_3 * pc_z[k] * sni_195[k];

        t_252[k] = pb_z[k] * smk0_108[k]
                   - f_12 * pc_z[k] * smk1_108[k];

        t_253[k] = f_14 * smi_112[k]
                   + f_3 * pc_y[k] * sni_196[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_z, pc_y, pc_z, smk0_111, smi_84, smi_114, \
                         smk1_111, sni_196, sni_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_13 * smi_84[k]
                   + f_3 * pc_z[k] * sni_196[k];

        t_255[k] = pb_z[k] * smk0_111[k]
                   - f_12 * pc_z[k] * smk1_111[k];

        t_256[k] = f_14 * smi_114[k]
                   + f_3 * pc_y[k] * sni_198[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_z, pc_x, pc_z, smk0_114, smi_87, smi_201, \
                         smk1_114, snh0_152, snh1_152, sni_199, \
                         sni_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_20 * smi_201[k]
                   + f_4 * snh0_152[k]
                   - f_5 * snh1_152[k]
                   + f_3 * pc_x[k] * sni_201[k];

        t_258[k] = pb_z[k] * smk0_114[k]
                   - f_12 * pc_z[k] * smk1_114[k];

        t_259[k] = f_13 * smi_87[k]
                   + f_3 * pc_z[k] * sni_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_z, pc_x, pc_y, pc_z, smk0_118, smi_117, \
                         smi_205, smk1_118, snh0_156, snh1_156, sni_201, \
                         sni_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * smi_117[k]
                   + f_3 * pc_y[k] * sni_201[k];

        t_261[k] = f_20 * smi_205[k]
                   + f_6 * snh0_156[k]
                   - f_7 * snh1_156[k]
                   + f_3 * pc_x[k] * sni_205[k];

        t_262[k] = pb_z[k] * smk0_118[k]
                   - f_12 * pc_z[k] * smk1_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_z, pc_y, pc_z, smk0_120, smi_90, smi_91, \
                         smi_121, smk1_120, sni_202, sni_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_13 * smi_90[k]
                   + f_3 * pc_z[k] * sni_202[k];

        t_264[k] = pb_z[k] * smk0_120[k]
                   + f_14 * smi_91[k]
                   - f_12 * pc_z[k] * smk1_120[k];

        t_265[k] = f_14 * smi_121[k]
                   + f_3 * pc_y[k] * sni_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_z, pc_x, pc_z, smk0_123, smi_94, smi_210, \
                         smk1_123, snh0_161, snh1_161, sni_206, \
                         sni_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_20 * smi_210[k]
                   + f_8 * snh0_161[k]
                   - f_9 * snh1_161[k]
                   + f_3 * pc_x[k] * sni_210[k];

        t_267[k] = pb_z[k] * smk0_123[k]
                   - f_12 * pc_z[k] * smk1_123[k];

        t_268[k] = f_13 * smi_94[k]
                   + f_3 * pc_z[k] * sni_206[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_z, pc_y, pc_z, smk0_125, smk0_126, smi_95, \
                         smi_96, smi_126, smk1_125, smk1_126, sni_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pb_z[k] * smk0_125[k]
                   + f_14 * smi_95[k]
                   - f_12 * pc_z[k] * smk1_125[k];

        t_270[k] = pb_z[k] * smk0_126[k]
                   + f_15 * smi_96[k]
                   - f_12 * pc_z[k] * smk1_126[k];

        t_271[k] = f_14 * smi_126[k]
                   + f_3 * pc_y[k] * sni_210[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, smi_216, smi_217, smi_218, smi_219, \
                         snh0_167, snh1_167, sni_216, sni_217, sni_218, \
                         sni_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_20 * smi_216[k]
                   + f_10 * snh0_167[k]
                   - f_11 * snh1_167[k]
                   + f_3 * pc_x[k] * sni_216[k];

        t_273[k] = f_20 * smi_217[k]
                   + f_3 * pc_x[k] * sni_217[k];

        t_274[k] = f_20 * smi_218[k]
                   + f_3 * pc_x[k] * sni_218[k];

        t_275[k] = f_20 * smi_219[k]
                   + f_3 * pc_x[k] * sni_219[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, smi_220, smi_221, smi_222, smi_223, \
                         sni_220, sni_221, sni_222, sni_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_20 * smi_220[k]
                   + f_3 * pc_x[k] * sni_220[k];

        t_277[k] = f_20 * smi_221[k]
                   + f_3 * pc_x[k] * sni_221[k];

        t_278[k] = f_20 * smi_222[k]
                   + f_3 * pc_x[k] * sni_222[k];

        t_279[k] = f_20 * smi_223[k]
                   + f_3 * pc_x[k] * sni_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_z, pc_y, pc_z, smk0_136, smi_105, smi_135, \
                         smk1_136, snh0_164, snh1_164, sni_217, \
                         sni_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pb_z[k] * smk0_136[k]
                   - f_12 * pc_z[k] * smk1_136[k];

        t_281[k] = f_13 * smi_105[k]
                   + f_3 * pc_z[k] * sni_217[k];

        t_282[k] = f_14 * smi_135[k]
                   + f_4 * snh0_164[k]
                   - f_5 * snh1_164[k]
                   + f_3 * pc_y[k] * sni_219[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pc_y, smi_136, smi_137, smi_138, snh0_165, \
                         snh0_166, snh0_167, snh1_165, snh1_166, snh1_167, sni_220, sni_221, \
                         sni_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_14 * smi_136[k]
                   + f_6 * snh0_165[k]
                   - f_7 * snh1_165[k]
                   + f_3 * pc_y[k] * sni_220[k];

        t_284[k] = f_14 * smi_137[k]
                   + f_8 * snh0_166[k]
                   - f_9 * snh1_166[k]
                   + f_3 * pc_y[k] * sni_221[k];

        t_285[k] = f_14 * smi_138[k]
                   + f_10 * snh0_167[k]
                   - f_11 * snh1_167[k]
                   + f_3 * pc_y[k] * sni_222[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, smk0_180, smi_111, \
                         smi_139, smi_140, smk1_180, snh0_167, snh1_167, sni_223, \
                         sni_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_14 * smi_139[k]
                   + f_3 * pc_y[k] * sni_223[k];

        t_287[k] = f_13 * smi_111[k]
                   + f_1 * snh0_167[k]
                   - f_2 * snh1_167[k]
                   + f_3 * pc_z[k] * sni_223[k];

        t_288[k] = pb_y[k] * smk0_180[k]
                   - f_12 * pc_y[k] * smk1_180[k];

        t_289[k] = f_13 * smi_140[k]
                   + f_3 * pc_y[k] * sni_224[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_y, pc_z, smk0_183, smk0_185, \
                         smi_112, smi_141, smi_142, smk1_183, smk1_185, sni_224, \
                         sni_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_14 * smi_112[k]
                   + f_3 * pc_z[k] * sni_224[k];

        t_291[k] = pb_y[k] * smk0_183[k]
                   + f_14 * smi_141[k]
                   - f_12 * pc_y[k] * smk1_183[k];

        t_292[k] = f_13 * smi_142[k]
                   + f_3 * pc_y[k] * sni_226[k];

        t_293[k] = pb_y[k] * smk0_185[k]
                   - f_12 * pc_y[k] * smk1_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pc_y, pc_z, smk0_186, smk0_189, \
                         smi_115, smi_143, smi_145, smk1_186, smk1_189, sni_227, \
                         sni_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pb_y[k] * smk0_186[k]
                   + f_15 * smi_143[k]
                   - f_12 * pc_y[k] * smk1_186[k];

        t_295[k] = f_14 * smi_115[k]
                   + f_3 * pc_z[k] * sni_227[k];

        t_296[k] = f_13 * smi_145[k]
                   + f_3 * pc_y[k] * sni_229[k];

        t_297[k] = pb_y[k] * smk0_189[k]
                   - f_12 * pc_y[k] * smk1_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pb_y, pc_y, pc_z, smk0_190, smk0_192, smi_118, \
                         smi_146, smi_148, smk1_190, smk1_192, \
                         sni_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = pb_y[k] * smk0_190[k]
                   + f_16 * smi_146[k]
                   - f_12 * pc_y[k] * smk1_190[k];

        t_299[k] = f_14 * smi_118[k]
                   + f_3 * pc_z[k] * sni_230[k];

        t_300[k] = pb_y[k] * smk0_192[k]
                   + f_14 * smi_148[k]
                   - f_12 * pc_y[k] * smk1_192[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_y, pc_y, pc_z, smk0_194, smk0_195, \
                         smi_122, smi_149, smi_150, smk1_194, smk1_195, sni_233, \
                         sni_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_13 * smi_149[k]
                   + f_3 * pc_y[k] * sni_233[k];

        t_302[k] = pb_y[k] * smk0_194[k]
                   - f_12 * pc_y[k] * smk1_194[k];

        t_303[k] = pb_y[k] * smk0_195[k]
                   + f_17 * smi_150[k]
                   - f_12 * pc_y[k] * smk1_195[k];

        t_304[k] = f_14 * smi_122[k]
                   + f_3 * pc_z[k] * sni_234[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_y, pc_y, smk0_197, smk0_198, smk0_200, \
                         smi_152, smi_153, smi_154, smk1_197, smk1_198, smk1_200, \
                         sni_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_y[k] * smk0_197[k]
                   + f_15 * smi_152[k]
                   - f_12 * pc_y[k] * smk1_197[k];

        t_306[k] = pb_y[k] * smk0_198[k]
                   + f_14 * smi_153[k]
                   - f_12 * pc_y[k] * smk1_198[k];

        t_307[k] = f_13 * smi_154[k]
                   + f_3 * pc_y[k] * sni_238[k];

        t_308[k] = pb_y[k] * smk0_200[k]
                   - f_12 * pc_y[k] * smk1_200[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pc_x, smi_245, smi_246, smi_247, \
                         smi_248, smi_249, sni_245, sni_246, sni_247, sni_248, \
                         sni_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_20 * smi_245[k]
                   + f_3 * pc_x[k] * sni_245[k];

        t_310[k] = f_20 * smi_246[k]
                   + f_3 * pc_x[k] * sni_246[k];

        t_311[k] = f_20 * smi_247[k]
                   + f_3 * pc_x[k] * sni_247[k];

        t_312[k] = f_20 * smi_248[k]
                   + f_3 * pc_x[k] * sni_248[k];

        t_313[k] = f_20 * smi_249[k]
                   + f_3 * pc_x[k] * sni_249[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, smi_133, smi_161, \
                         smi_250, smi_251, snh0_183, snh1_183, sni_245, sni_250, \
                         sni_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_20 * smi_250[k]
                   + f_3 * pc_x[k] * sni_250[k];

        t_315[k] = f_20 * smi_251[k]
                   + f_3 * pc_x[k] * sni_251[k];

        t_316[k] = f_13 * smi_161[k]
                   + f_1 * snh0_183[k]
                   - f_2 * snh1_183[k]
                   + f_3 * pc_y[k] * sni_245[k];

        t_317[k] = f_14 * smi_133[k]
                   + f_3 * pc_z[k] * sni_245[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_y, smi_163, smi_164, smi_165, snh0_185, \
                         snh0_186, snh0_187, snh1_185, snh1_186, snh1_187, sni_247, sni_248, \
                         sni_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_13 * smi_163[k]
                   + f_4 * snh0_185[k]
                   - f_5 * snh1_185[k]
                   + f_3 * pc_y[k] * sni_247[k];

        t_319[k] = f_13 * smi_164[k]
                   + f_6 * snh0_186[k]
                   - f_7 * snh1_186[k]
                   + f_3 * pc_y[k] * sni_248[k];

        t_320[k] = f_13 * smi_165[k]
                   + f_8 * snh0_187[k]
                   - f_9 * snh1_187[k]
                   + f_3 * pc_y[k] * sni_249[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_y, pc_y, smk0_215, smi_166, smi_167, \
                         smk1_215, snh0_188, snh1_188, sni_250, \
                         sni_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * smi_166[k]
                   + f_10 * snh0_188[k]
                   - f_11 * snh1_188[k]
                   + f_3 * pc_y[k] * sni_250[k];

        t_322[k] = f_13 * smi_167[k]
                   + f_3 * pc_y[k] * sni_251[k];

        t_323[k] = pb_y[k] * smk0_215[k]
                   - f_12 * pc_y[k] * smk1_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pc_x, pc_y, pc_z, smi_140, smi_252, \
                         smi_255, snh0_189, snh0_192, snh1_189, snh1_192, sni_252, \
                         sni_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_20 * smi_252[k]
                   + f_1 * snh0_189[k]
                   - f_2 * snh1_189[k]
                   + f_3 * pc_x[k] * sni_252[k];

        t_325[k] = f_3 * pc_y[k] * sni_252[k];

        t_326[k] = f_15 * smi_140[k]
                   + f_3 * pc_z[k] * sni_252[k];

        t_327[k] = f_20 * smi_255[k]
                   + f_4 * snh0_192[k]
                   - f_5 * snh1_192[k]
                   + f_3 * pc_x[k] * sni_255[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pc_x, pc_y, smi_257, smi_258, snh0_194, \
                         snh0_195, snh1_194, snh1_195, sni_254, sni_257, \
                         sni_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_3 * pc_y[k] * sni_254[k];

        t_329[k] = f_20 * smi_257[k]
                   + f_4 * snh0_194[k]
                   - f_5 * snh1_194[k]
                   + f_3 * pc_x[k] * sni_257[k];

        t_330[k] = f_20 * smi_258[k]
                   + f_6 * snh0_195[k]
                   - f_7 * snh1_195[k]
                   + f_3 * pc_x[k] * sni_258[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_x, pc_y, pc_z, smi_143, smi_261, snh0_198, \
                         snh1_198, sni_255, sni_257, sni_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_15 * smi_143[k]
                   + f_3 * pc_z[k] * sni_255[k];

        t_332[k] = f_3 * pc_y[k] * sni_257[k];

        t_333[k] = f_20 * smi_261[k]
                   + f_6 * snh0_198[k]
                   - f_7 * snh1_198[k]
                   + f_3 * pc_x[k] * sni_261[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pc_x, pc_z, smi_146, smi_262, smi_264, snh0_199, \
                         snh0_201, snh1_199, snh1_201, sni_258, sni_262, \
                         sni_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_20 * smi_262[k]
                   + f_8 * snh0_199[k]
                   - f_9 * snh1_199[k]
                   + f_3 * pc_x[k] * sni_262[k];

        t_335[k] = f_15 * smi_146[k]
                   + f_3 * pc_z[k] * sni_258[k];

        t_336[k] = f_20 * smi_264[k]
                   + f_8 * snh0_201[k]
                   - f_9 * snh1_201[k]
                   + f_3 * pc_x[k] * sni_264[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pc_x, pc_y, smi_266, smi_267, snh0_203, \
                         snh0_204, snh1_203, snh1_204, sni_261, sni_266, \
                         sni_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_3 * pc_y[k] * sni_261[k];

        t_338[k] = f_20 * smi_266[k]
                   + f_8 * snh0_203[k]
                   - f_9 * snh1_203[k]
                   + f_3 * pc_x[k] * sni_266[k];

        t_339[k] = f_20 * smi_267[k]
                   + f_10 * snh0_204[k]
                   - f_11 * snh1_204[k]
                   + f_3 * pc_x[k] * sni_267[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pc_x, pc_z, smi_150, smi_269, smi_270, snh0_206, \
                         snh0_207, snh1_206, snh1_207, sni_262, sni_269, \
                         sni_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_15 * smi_150[k]
                   + f_3 * pc_z[k] * sni_262[k];

        t_341[k] = f_20 * smi_269[k]
                   + f_10 * snh0_206[k]
                   - f_11 * snh1_206[k]
                   + f_3 * pc_x[k] * sni_269[k];

        t_342[k] = f_20 * smi_270[k]
                   + f_10 * snh0_207[k]
                   - f_11 * snh1_207[k]
                   + f_3 * pc_x[k] * sni_270[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pc_x, pc_y, smi_272, smi_273, smi_274, \
                         snh0_209, snh1_209, sni_266, sni_272, sni_273, \
                         sni_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_3 * pc_y[k] * sni_266[k];

        t_344[k] = f_20 * smi_272[k]
                   + f_10 * snh0_209[k]
                   - f_11 * snh1_209[k]
                   + f_3 * pc_x[k] * sni_272[k];

        t_345[k] = f_20 * smi_273[k]
                   + f_3 * pc_x[k] * sni_273[k];

        t_346[k] = f_20 * smi_274[k]
                   + f_3 * pc_x[k] * sni_274[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *smk0_216 = buffer.data(smk0 + 216);
    const auto *smk0_219 = buffer.data(smk0 + 219);
    const auto *smk0_222 = buffer.data(smk0 + 222);
    const auto *smk0_226 = buffer.data(smk0 + 226);
    const auto *smk0_228 = buffer.data(smk0 + 228);
    const auto *smk0_231 = buffer.data(smk0 + 231);
    const auto *smk0_233 = buffer.data(smk0 + 233);
    const auto *smk0_234 = buffer.data(smk0 + 234);
    const auto *smk0_244 = buffer.data(smk0 + 244);

    const auto *smi_161 = buffer.data(smi + 161);
    const auto *smi_167 = buffer.data(smi + 167);
    const auto *smi_168 = buffer.data(smi + 168);
    const auto *smi_170 = buffer.data(smi + 170);
    const auto *smi_171 = buffer.data(smi + 171);
    const auto *smi_173 = buffer.data(smi + 173);
    const auto *smi_174 = buffer.data(smi + 174);
    const auto *smi_175 = buffer.data(smi + 175);
    const auto *smi_177 = buffer.data(smi + 177);
    const auto *smi_178 = buffer.data(smi + 178);
    const auto *smi_179 = buffer.data(smi + 179);
    const auto *smi_180 = buffer.data(smi + 180);
    const auto *smi_182 = buffer.data(smi + 182);
    const auto *smi_189 = buffer.data(smi + 189);
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
    const auto *smi_219 = buffer.data(smi + 219);
    const auto *smi_220 = buffer.data(smi + 220);
    const auto *smi_221 = buffer.data(smi + 221);
    const auto *smi_222 = buffer.data(smi + 222);
    const auto *smi_223 = buffer.data(smi + 223);
    const auto *smi_224 = buffer.data(smi + 224);
    const auto *smi_226 = buffer.data(smi + 226);
    const auto *smi_229 = buffer.data(smi + 229);
    const auto *smi_233 = buffer.data(smi + 233);
    const auto *smi_238 = buffer.data(smi + 238);
    const auto *smi_275 = buffer.data(smi + 275);
    const auto *smi_276 = buffer.data(smi + 276);
    const auto *smi_277 = buffer.data(smi + 277);
    const auto *smi_278 = buffer.data(smi + 278);
    const auto *smi_279 = buffer.data(smi + 279);
    const auto *smi_280 = buffer.data(smi + 280);
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
    const auto *smi_313 = buffer.data(smi + 313);
    const auto *smi_317 = buffer.data(smi + 317);
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

    const auto *smk1_216 = buffer.data(smk1 + 216);
    const auto *smk1_219 = buffer.data(smk1 + 219);
    const auto *smk1_222 = buffer.data(smk1 + 222);
    const auto *smk1_226 = buffer.data(smk1 + 226);
    const auto *smk1_228 = buffer.data(smk1 + 228);
    const auto *smk1_231 = buffer.data(smk1 + 231);
    const auto *smk1_233 = buffer.data(smk1 + 233);
    const auto *smk1_234 = buffer.data(smk1 + 234);
    const auto *smk1_244 = buffer.data(smk1 + 244);

    const auto *snh0_204 = buffer.data(snh0 + 204);
    const auto *snh0_206 = buffer.data(snh0 + 206);
    const auto *snh0_207 = buffer.data(snh0 + 207);
    const auto *snh0_208 = buffer.data(snh0 + 208);
    const auto *snh0_209 = buffer.data(snh0 + 209);
    const auto *snh0_210 = buffer.data(snh0 + 210);
    const auto *snh0_213 = buffer.data(snh0 + 213);
    const auto *snh0_215 = buffer.data(snh0 + 215);
    const auto *snh0_216 = buffer.data(snh0 + 216);
    const auto *snh0_219 = buffer.data(snh0 + 219);
    const auto *snh0_220 = buffer.data(snh0 + 220);
    const auto *snh0_222 = buffer.data(snh0 + 222);
    const auto *snh0_224 = buffer.data(snh0 + 224);
    const auto *snh0_225 = buffer.data(snh0 + 225);
    const auto *snh0_227 = buffer.data(snh0 + 227);
    const auto *snh0_228 = buffer.data(snh0 + 228);
    const auto *snh0_229 = buffer.data(snh0 + 229);
    const auto *snh0_230 = buffer.data(snh0 + 230);
    const auto *snh0_236 = buffer.data(snh0 + 236);
    const auto *snh0_240 = buffer.data(snh0 + 240);
    const auto *snh0_245 = buffer.data(snh0 + 245);
    const auto *snh0_248 = buffer.data(snh0 + 248);
    const auto *snh0_249 = buffer.data(snh0 + 249);
    const auto *snh0_250 = buffer.data(snh0 + 250);
    const auto *snh0_251 = buffer.data(snh0 + 251);
    const auto *snh0_252 = buffer.data(snh0 + 252);
    const auto *snh0_255 = buffer.data(snh0 + 255);
    const auto *snh0_257 = buffer.data(snh0 + 257);
    const auto *snh0_258 = buffer.data(snh0 + 258);
    const auto *snh0_261 = buffer.data(snh0 + 261);
    const auto *snh0_262 = buffer.data(snh0 + 262);
    const auto *snh0_264 = buffer.data(snh0 + 264);
    const auto *snh0_266 = buffer.data(snh0 + 266);
    const auto *snh0_267 = buffer.data(snh0 + 267);
    const auto *snh0_269 = buffer.data(snh0 + 269);
    const auto *snh0_270 = buffer.data(snh0 + 270);
    const auto *snh0_272 = buffer.data(snh0 + 272);

    const auto *snh1_204 = buffer.data(snh1 + 204);
    const auto *snh1_206 = buffer.data(snh1 + 206);
    const auto *snh1_207 = buffer.data(snh1 + 207);
    const auto *snh1_208 = buffer.data(snh1 + 208);
    const auto *snh1_209 = buffer.data(snh1 + 209);
    const auto *snh1_210 = buffer.data(snh1 + 210);
    const auto *snh1_213 = buffer.data(snh1 + 213);
    const auto *snh1_215 = buffer.data(snh1 + 215);
    const auto *snh1_216 = buffer.data(snh1 + 216);
    const auto *snh1_219 = buffer.data(snh1 + 219);
    const auto *snh1_220 = buffer.data(snh1 + 220);
    const auto *snh1_222 = buffer.data(snh1 + 222);
    const auto *snh1_224 = buffer.data(snh1 + 224);
    const auto *snh1_225 = buffer.data(snh1 + 225);
    const auto *snh1_227 = buffer.data(snh1 + 227);
    const auto *snh1_228 = buffer.data(snh1 + 228);
    const auto *snh1_229 = buffer.data(snh1 + 229);
    const auto *snh1_230 = buffer.data(snh1 + 230);
    const auto *snh1_236 = buffer.data(snh1 + 236);
    const auto *snh1_240 = buffer.data(snh1 + 240);
    const auto *snh1_245 = buffer.data(snh1 + 245);
    const auto *snh1_248 = buffer.data(snh1 + 248);
    const auto *snh1_249 = buffer.data(snh1 + 249);
    const auto *snh1_250 = buffer.data(snh1 + 250);
    const auto *snh1_251 = buffer.data(snh1 + 251);
    const auto *snh1_252 = buffer.data(snh1 + 252);
    const auto *snh1_255 = buffer.data(snh1 + 255);
    const auto *snh1_257 = buffer.data(snh1 + 257);
    const auto *snh1_258 = buffer.data(snh1 + 258);
    const auto *snh1_261 = buffer.data(snh1 + 261);
    const auto *snh1_262 = buffer.data(snh1 + 262);
    const auto *snh1_264 = buffer.data(snh1 + 264);
    const auto *snh1_266 = buffer.data(snh1 + 266);
    const auto *snh1_267 = buffer.data(snh1 + 267);
    const auto *snh1_269 = buffer.data(snh1 + 269);
    const auto *snh1_270 = buffer.data(snh1 + 270);
    const auto *snh1_272 = buffer.data(snh1 + 272);

    const auto *sni_273 = buffer.data(sni + 273);
    const auto *sni_275 = buffer.data(sni + 275);
    const auto *sni_276 = buffer.data(sni + 276);
    const auto *sni_277 = buffer.data(sni + 277);
    const auto *sni_278 = buffer.data(sni + 278);
    const auto *sni_279 = buffer.data(sni + 279);
    const auto *sni_280 = buffer.data(sni + 280);
    const auto *sni_282 = buffer.data(sni + 282);
    const auto *sni_283 = buffer.data(sni + 283);
    const auto *sni_285 = buffer.data(sni + 285);
    const auto *sni_286 = buffer.data(sni + 286);
    const auto *sni_289 = buffer.data(sni + 289);
    const auto *sni_290 = buffer.data(sni + 290);
    const auto *sni_292 = buffer.data(sni + 292);
    const auto *sni_294 = buffer.data(sni + 294);
    const auto *sni_295 = buffer.data(sni + 295);
    const auto *sni_297 = buffer.data(sni + 297);
    const auto *sni_298 = buffer.data(sni + 298);
    const auto *sni_300 = buffer.data(sni + 300);
    const auto *sni_301 = buffer.data(sni + 301);
    const auto *sni_302 = buffer.data(sni + 302);
    const auto *sni_303 = buffer.data(sni + 303);
    const auto *sni_304 = buffer.data(sni + 304);
    const auto *sni_305 = buffer.data(sni + 305);
    const auto *sni_306 = buffer.data(sni + 306);
    const auto *sni_307 = buffer.data(sni + 307);
    const auto *sni_308 = buffer.data(sni + 308);
    const auto *sni_310 = buffer.data(sni + 310);
    const auto *sni_311 = buffer.data(sni + 311);
    const auto *sni_313 = buffer.data(sni + 313);
    const auto *sni_314 = buffer.data(sni + 314);
    const auto *sni_317 = buffer.data(sni + 317);
    const auto *sni_318 = buffer.data(sni + 318);
    const auto *sni_322 = buffer.data(sni + 322);
    const auto *sni_328 = buffer.data(sni + 328);
    const auto *sni_329 = buffer.data(sni + 329);
    const auto *sni_330 = buffer.data(sni + 330);
    const auto *sni_331 = buffer.data(sni + 331);
    const auto *sni_332 = buffer.data(sni + 332);
    const auto *sni_333 = buffer.data(sni + 333);
    const auto *sni_334 = buffer.data(sni + 334);
    const auto *sni_335 = buffer.data(sni + 335);
    const auto *sni_336 = buffer.data(sni + 336);
    const auto *sni_338 = buffer.data(sni + 338);
    const auto *sni_339 = buffer.data(sni + 339);
    const auto *sni_341 = buffer.data(sni + 341);
    const auto *sni_342 = buffer.data(sni + 342);
    const auto *sni_345 = buffer.data(sni + 345);
    const auto *sni_346 = buffer.data(sni + 346);
    const auto *sni_348 = buffer.data(sni + 348);
    const auto *sni_350 = buffer.data(sni + 350);
    const auto *sni_351 = buffer.data(sni + 351);
    const auto *sni_353 = buffer.data(sni + 353);
    const auto *sni_354 = buffer.data(sni + 354);
    const auto *sni_356 = buffer.data(sni + 356);
    const auto *sni_357 = buffer.data(sni + 357);
    const auto *sni_358 = buffer.data(sni + 358);
    const auto *sni_359 = buffer.data(sni + 359);

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pc_x, smi_275, smi_276, smi_277, \
                         smi_278, smi_279, sni_275, sni_276, sni_277, sni_278, \
                         sni_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_20 * smi_275[k]
                   + f_3 * pc_x[k] * sni_275[k];

        t_348[k] = f_20 * smi_276[k]
                   + f_3 * pc_x[k] * sni_276[k];

        t_349[k] = f_20 * smi_277[k]
                   + f_3 * pc_x[k] * sni_277[k];

        t_350[k] = f_20 * smi_278[k]
                   + f_3 * pc_x[k] * sni_278[k];

        t_351[k] = f_20 * smi_279[k]
                   + f_3 * pc_x[k] * sni_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_y, pc_z, smi_161, snh0_204, snh0_206, \
                         snh0_207, snh1_204, snh1_206, snh1_207, sni_273, sni_275, \
                         sni_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * snh0_204[k]
                   - f_2 * snh1_204[k]
                   + f_3 * pc_y[k] * sni_273[k];

        t_353[k] = f_15 * smi_161[k]
                   + f_3 * pc_z[k] * sni_273[k];

        t_354[k] = f_4 * snh0_206[k]
                   - f_5 * snh1_206[k]
                   + f_3 * pc_y[k] * sni_275[k];

        t_355[k] = f_6 * snh0_207[k]
                   - f_7 * snh1_207[k]
                   + f_3 * pc_y[k] * sni_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, smi_167, snh0_208, snh0_209, \
                         snh1_208, snh1_209, sni_277, sni_278, \
                         sni_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_8 * snh0_208[k]
                   - f_9 * snh1_208[k]
                   + f_3 * pc_y[k] * sni_277[k];

        t_357[k] = f_10 * snh0_209[k]
                   - f_11 * snh1_209[k]
                   + f_3 * pc_y[k] * sni_278[k];

        t_358[k] = f_3 * pc_y[k] * sni_279[k];

        t_359[k] = f_15 * smi_167[k]
                   + f_1 * snh0_209[k]
                   - f_2 * snh1_209[k]
                   + f_3 * pc_z[k] * sni_279[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, smi_168, smi_280, \
                         smi_283, snh0_210, snh0_213, snh1_210, snh1_213, sni_280, \
                         sni_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_21 * smi_280[k]
                   + f_1 * snh0_210[k]
                   - f_2 * snh1_210[k]
                   + f_3 * pc_x[k] * sni_280[k];

        t_361[k] = f_16 * smi_168[k]
                   + f_3 * pc_y[k] * sni_280[k];

        t_362[k] = f_3 * pc_z[k] * sni_280[k];

        t_363[k] = f_21 * smi_283[k]
                   + f_4 * snh0_213[k]
                   - f_5 * snh1_213[k]
                   + f_3 * pc_x[k] * sni_283[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pc_x, pc_y, smi_170, smi_285, smi_286, snh0_215, \
                         snh0_216, snh1_215, snh1_216, sni_282, sni_285, \
                         sni_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_16 * smi_170[k]
                   + f_3 * pc_y[k] * sni_282[k];

        t_365[k] = f_21 * smi_285[k]
                   + f_4 * snh0_215[k]
                   - f_5 * snh1_215[k]
                   + f_3 * pc_x[k] * sni_285[k];

        t_366[k] = f_21 * smi_286[k]
                   + f_6 * snh0_216[k]
                   - f_7 * snh1_216[k]
                   + f_3 * pc_x[k] * sni_286[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pc_x, pc_y, pc_z, smi_173, smi_289, snh0_219, \
                         snh1_219, sni_283, sni_285, sni_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_3 * pc_z[k] * sni_283[k];

        t_368[k] = f_16 * smi_173[k]
                   + f_3 * pc_y[k] * sni_285[k];

        t_369[k] = f_21 * smi_289[k]
                   + f_6 * snh0_219[k]
                   - f_7 * snh1_219[k]
                   + f_3 * pc_x[k] * sni_289[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pc_x, pc_z, smi_290, smi_292, snh0_220, \
                         snh0_222, snh1_220, snh1_222, sni_286, sni_290, \
                         sni_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_21 * smi_290[k]
                   + f_8 * snh0_220[k]
                   - f_9 * snh1_220[k]
                   + f_3 * pc_x[k] * sni_290[k];

        t_371[k] = f_3 * pc_z[k] * sni_286[k];

        t_372[k] = f_21 * smi_292[k]
                   + f_8 * snh0_222[k]
                   - f_9 * snh1_222[k]
                   + f_3 * pc_x[k] * sni_292[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, smi_177, smi_294, smi_295, snh0_224, \
                         snh0_225, snh1_224, snh1_225, sni_289, sni_294, \
                         sni_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_16 * smi_177[k]
                   + f_3 * pc_y[k] * sni_289[k];

        t_374[k] = f_21 * smi_294[k]
                   + f_8 * snh0_224[k]
                   - f_9 * snh1_224[k]
                   + f_3 * pc_x[k] * sni_294[k];

        t_375[k] = f_21 * smi_295[k]
                   + f_10 * snh0_225[k]
                   - f_11 * snh1_225[k]
                   + f_3 * pc_x[k] * sni_295[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_z, smi_297, smi_298, snh0_227, \
                         snh0_228, snh1_227, snh1_228, sni_290, sni_297, \
                         sni_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_3 * pc_z[k] * sni_290[k];

        t_377[k] = f_21 * smi_297[k]
                   + f_10 * snh0_227[k]
                   - f_11 * snh1_227[k]
                   + f_3 * pc_x[k] * sni_297[k];

        t_378[k] = f_21 * smi_298[k]
                   + f_10 * snh0_228[k]
                   - f_11 * snh1_228[k]
                   + f_3 * pc_x[k] * sni_298[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, smi_182, smi_300, smi_301, \
                         smi_302, snh0_230, snh1_230, sni_294, sni_300, sni_301, \
                         sni_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_16 * smi_182[k]
                   + f_3 * pc_y[k] * sni_294[k];

        t_380[k] = f_21 * smi_300[k]
                   + f_10 * snh0_230[k]
                   - f_11 * snh1_230[k]
                   + f_3 * pc_x[k] * sni_300[k];

        t_381[k] = f_21 * smi_301[k]
                   + f_3 * pc_x[k] * sni_301[k];

        t_382[k] = f_21 * smi_302[k]
                   + f_3 * pc_x[k] * sni_302[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pc_x, smi_303, smi_304, smi_305, \
                         smi_306, smi_307, sni_303, sni_304, sni_305, sni_306, \
                         sni_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_21 * smi_303[k]
                   + f_3 * pc_x[k] * sni_303[k];

        t_384[k] = f_21 * smi_304[k]
                   + f_3 * pc_x[k] * sni_304[k];

        t_385[k] = f_21 * smi_305[k]
                   + f_3 * pc_x[k] * sni_305[k];

        t_386[k] = f_21 * smi_306[k]
                   + f_3 * pc_x[k] * sni_306[k];

        t_387[k] = f_21 * smi_307[k]
                   + f_3 * pc_x[k] * sni_307[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, pc_z, smi_189, smi_191, snh0_225, \
                         snh0_227, snh1_225, snh1_227, sni_301, \
                         sni_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_16 * smi_189[k]
                   + f_1 * snh0_225[k]
                   - f_2 * snh1_225[k]
                   + f_3 * pc_y[k] * sni_301[k];

        t_389[k] = f_3 * pc_z[k] * sni_301[k];

        t_390[k] = f_16 * smi_191[k]
                   + f_4 * snh0_227[k]
                   - f_5 * snh1_227[k]
                   + f_3 * pc_y[k] * sni_303[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_y, smi_192, smi_193, smi_194, snh0_228, \
                         snh0_229, snh0_230, snh1_228, snh1_229, snh1_230, sni_304, sni_305, \
                         sni_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_16 * smi_192[k]
                   + f_6 * snh0_228[k]
                   - f_7 * snh1_228[k]
                   + f_3 * pc_y[k] * sni_304[k];

        t_392[k] = f_16 * smi_193[k]
                   + f_8 * snh0_229[k]
                   - f_9 * snh1_229[k]
                   + f_3 * pc_y[k] * sni_305[k];

        t_393[k] = f_16 * smi_194[k]
                   + f_10 * snh0_230[k]
                   - f_11 * snh1_230[k]
                   + f_3 * pc_y[k] * sni_306[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pb_z, pc_y, pc_z, smk0_216, smi_195, \
                         smi_196, smk1_216, snh0_230, snh1_230, sni_307, \
                         sni_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_16 * smi_195[k]
                   + f_3 * pc_y[k] * sni_307[k];

        t_395[k] = f_1 * snh0_230[k]
                   - f_2 * snh1_230[k]
                   + f_3 * pc_z[k] * sni_307[k];

        t_396[k] = pb_z[k] * smk0_216[k]
                   - f_12 * pc_z[k] * smk1_216[k];

        t_397[k] = f_15 * smi_196[k]
                   + f_3 * pc_y[k] * sni_308[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pb_z, pc_y, pc_z, smk0_219, smi_168, smi_198, \
                         smk1_219, sni_308, sni_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_13 * smi_168[k]
                   + f_3 * pc_z[k] * sni_308[k];

        t_399[k] = pb_z[k] * smk0_219[k]
                   - f_12 * pc_z[k] * smk1_219[k];

        t_400[k] = f_15 * smi_198[k]
                   + f_3 * pc_y[k] * sni_310[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pb_z, pc_x, pc_z, smk0_222, smi_171, smi_313, \
                         smk1_222, snh0_236, snh1_236, sni_311, \
                         sni_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_21 * smi_313[k]
                   + f_4 * snh0_236[k]
                   - f_5 * snh1_236[k]
                   + f_3 * pc_x[k] * sni_313[k];

        t_402[k] = pb_z[k] * smk0_222[k]
                   - f_12 * pc_z[k] * smk1_222[k];

        t_403[k] = f_13 * smi_171[k]
                   + f_3 * pc_z[k] * sni_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pb_z, pc_x, pc_y, pc_z, smk0_226, smi_201, \
                         smi_317, smk1_226, snh0_240, snh1_240, sni_313, \
                         sni_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_15 * smi_201[k]
                   + f_3 * pc_y[k] * sni_313[k];

        t_405[k] = f_21 * smi_317[k]
                   + f_6 * snh0_240[k]
                   - f_7 * snh1_240[k]
                   + f_3 * pc_x[k] * sni_317[k];

        t_406[k] = pb_z[k] * smk0_226[k]
                   - f_12 * pc_z[k] * smk1_226[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_z, pc_y, pc_z, smk0_228, smi_174, smi_175, \
                         smi_205, smk1_228, sni_314, sni_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * smi_174[k]
                   + f_3 * pc_z[k] * sni_314[k];

        t_408[k] = pb_z[k] * smk0_228[k]
                   + f_14 * smi_175[k]
                   - f_12 * pc_z[k] * smk1_228[k];

        t_409[k] = f_15 * smi_205[k]
                   + f_3 * pc_y[k] * sni_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pb_z, pc_x, pc_z, smk0_231, smi_178, smi_322, \
                         smk1_231, snh0_245, snh1_245, sni_318, \
                         sni_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_21 * smi_322[k]
                   + f_8 * snh0_245[k]
                   - f_9 * snh1_245[k]
                   + f_3 * pc_x[k] * sni_322[k];

        t_411[k] = pb_z[k] * smk0_231[k]
                   - f_12 * pc_z[k] * smk1_231[k];

        t_412[k] = f_13 * smi_178[k]
                   + f_3 * pc_z[k] * sni_318[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pb_z, pc_y, pc_z, smk0_233, smk0_234, smi_179, \
                         smi_180, smi_210, smk1_233, smk1_234, \
                         sni_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_z[k] * smk0_233[k]
                   + f_14 * smi_179[k]
                   - f_12 * pc_z[k] * smk1_233[k];

        t_414[k] = pb_z[k] * smk0_234[k]
                   + f_15 * smi_180[k]
                   - f_12 * pc_z[k] * smk1_234[k];

        t_415[k] = f_15 * smi_210[k]
                   + f_3 * pc_y[k] * sni_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_x, smi_328, smi_329, smi_330, smi_331, \
                         snh0_251, snh1_251, sni_328, sni_329, sni_330, \
                         sni_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_21 * smi_328[k]
                   + f_10 * snh0_251[k]
                   - f_11 * snh1_251[k]
                   + f_3 * pc_x[k] * sni_328[k];

        t_417[k] = f_21 * smi_329[k]
                   + f_3 * pc_x[k] * sni_329[k];

        t_418[k] = f_21 * smi_330[k]
                   + f_3 * pc_x[k] * sni_330[k];

        t_419[k] = f_21 * smi_331[k]
                   + f_3 * pc_x[k] * sni_331[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, smi_332, smi_333, smi_334, smi_335, \
                         sni_332, sni_333, sni_334, sni_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_21 * smi_332[k]
                   + f_3 * pc_x[k] * sni_332[k];

        t_421[k] = f_21 * smi_333[k]
                   + f_3 * pc_x[k] * sni_333[k];

        t_422[k] = f_21 * smi_334[k]
                   + f_3 * pc_x[k] * sni_334[k];

        t_423[k] = f_21 * smi_335[k]
                   + f_3 * pc_x[k] * sni_335[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_z, pc_y, pc_z, smk0_244, smi_189, smi_219, \
                         smk1_244, snh0_248, snh1_248, sni_329, \
                         sni_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_z[k] * smk0_244[k]
                   - f_12 * pc_z[k] * smk1_244[k];

        t_425[k] = f_13 * smi_189[k]
                   + f_3 * pc_z[k] * sni_329[k];

        t_426[k] = f_15 * smi_219[k]
                   + f_4 * snh0_248[k]
                   - f_5 * snh1_248[k]
                   + f_3 * pc_y[k] * sni_331[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_y, smi_220, smi_221, smi_222, snh0_249, \
                         snh0_250, snh0_251, snh1_249, snh1_250, snh1_251, sni_332, sni_333, \
                         sni_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_15 * smi_220[k]
                   + f_6 * snh0_249[k]
                   - f_7 * snh1_249[k]
                   + f_3 * pc_y[k] * sni_332[k];

        t_428[k] = f_15 * smi_221[k]
                   + f_8 * snh0_250[k]
                   - f_9 * snh1_250[k]
                   + f_3 * pc_y[k] * sni_333[k];

        t_429[k] = f_15 * smi_222[k]
                   + f_10 * snh0_251[k]
                   - f_11 * snh1_251[k]
                   + f_3 * pc_y[k] * sni_334[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_y, pc_z, smi_195, smi_223, smi_336, \
                         snh0_251, snh0_252, snh1_251, snh1_252, sni_335, \
                         sni_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_15 * smi_223[k]
                   + f_3 * pc_y[k] * sni_335[k];

        t_431[k] = f_13 * smi_195[k]
                   + f_1 * snh0_251[k]
                   - f_2 * snh1_251[k]
                   + f_3 * pc_z[k] * sni_335[k];

        t_432[k] = f_21 * smi_336[k]
                   + f_1 * snh0_252[k]
                   - f_2 * snh1_252[k]
                   + f_3 * pc_x[k] * sni_336[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, pc_z, smi_196, smi_224, \
                         smi_226, smi_339, snh0_255, snh1_255, sni_336, sni_338, \
                         sni_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_14 * smi_224[k]
                   + f_3 * pc_y[k] * sni_336[k];

        t_434[k] = f_14 * smi_196[k]
                   + f_3 * pc_z[k] * sni_336[k];

        t_435[k] = f_21 * smi_339[k]
                   + f_4 * snh0_255[k]
                   - f_5 * snh1_255[k]
                   + f_3 * pc_x[k] * sni_339[k];

        t_436[k] = f_14 * smi_226[k]
                   + f_3 * pc_y[k] * sni_338[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pc_x, pc_z, smi_199, smi_341, smi_342, snh0_257, \
                         snh0_258, snh1_257, snh1_258, sni_339, sni_341, \
                         sni_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_21 * smi_341[k]
                   + f_4 * snh0_257[k]
                   - f_5 * snh1_257[k]
                   + f_3 * pc_x[k] * sni_341[k];

        t_438[k] = f_21 * smi_342[k]
                   + f_6 * snh0_258[k]
                   - f_7 * snh1_258[k]
                   + f_3 * pc_x[k] * sni_342[k];

        t_439[k] = f_14 * smi_199[k]
                   + f_3 * pc_z[k] * sni_339[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pc_x, pc_y, smi_229, smi_345, smi_346, snh0_261, \
                         snh0_262, snh1_261, snh1_262, sni_341, sni_345, \
                         sni_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * smi_229[k]
                   + f_3 * pc_y[k] * sni_341[k];

        t_441[k] = f_21 * smi_345[k]
                   + f_6 * snh0_261[k]
                   - f_7 * snh1_261[k]
                   + f_3 * pc_x[k] * sni_345[k];

        t_442[k] = f_21 * smi_346[k]
                   + f_8 * snh0_262[k]
                   - f_9 * snh1_262[k]
                   + f_3 * pc_x[k] * sni_346[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, pc_x, pc_y, pc_z, smi_202, smi_233, smi_348, \
                         snh0_264, snh1_264, sni_342, sni_345, \
                         sni_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_14 * smi_202[k]
                   + f_3 * pc_z[k] * sni_342[k];

        t_444[k] = f_21 * smi_348[k]
                   + f_8 * snh0_264[k]
                   - f_9 * snh1_264[k]
                   + f_3 * pc_x[k] * sni_348[k];

        t_445[k] = f_14 * smi_233[k]
                   + f_3 * pc_y[k] * sni_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pc_x, pc_z, smi_206, smi_350, smi_351, snh0_266, \
                         snh0_267, snh1_266, snh1_267, sni_346, sni_350, \
                         sni_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_21 * smi_350[k]
                   + f_8 * snh0_266[k]
                   - f_9 * snh1_266[k]
                   + f_3 * pc_x[k] * sni_350[k];

        t_447[k] = f_21 * smi_351[k]
                   + f_10 * snh0_267[k]
                   - f_11 * snh1_267[k]
                   + f_3 * pc_x[k] * sni_351[k];

        t_448[k] = f_14 * smi_206[k]
                   + f_3 * pc_z[k] * sni_346[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_y, smi_238, smi_353, smi_354, snh0_269, \
                         snh0_270, snh1_269, snh1_270, sni_350, sni_353, \
                         sni_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_21 * smi_353[k]
                   + f_10 * snh0_269[k]
                   - f_11 * snh1_269[k]
                   + f_3 * pc_x[k] * sni_353[k];

        t_450[k] = f_21 * smi_354[k]
                   + f_10 * snh0_270[k]
                   - f_11 * snh1_270[k]
                   + f_3 * pc_x[k] * sni_354[k];

        t_451[k] = f_14 * smi_238[k]
                   + f_3 * pc_y[k] * sni_350[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, smi_356, smi_357, smi_358, smi_359, \
                         snh0_272, snh1_272, sni_356, sni_357, sni_358, \
                         sni_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_21 * smi_356[k]
                   + f_10 * snh0_272[k]
                   - f_11 * snh1_272[k]
                   + f_3 * pc_x[k] * sni_356[k];

        t_453[k] = f_21 * smi_357[k]
                   + f_3 * pc_x[k] * sni_357[k];

        t_454[k] = f_21 * smi_358[k]
                   + f_3 * pc_x[k] * sni_358[k];

        t_455[k] = f_21 * smi_359[k]
                   + f_3 * pc_x[k] * sni_359[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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
    const auto f_21 = 3.0 / q;

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

    const auto *smk0_324 = buffer.data(smk0 + 324);
    const auto *smk0_327 = buffer.data(smk0 + 327);
    const auto *smk0_329 = buffer.data(smk0 + 329);
    const auto *smk0_330 = buffer.data(smk0 + 330);
    const auto *smk0_333 = buffer.data(smk0 + 333);
    const auto *smk0_334 = buffer.data(smk0 + 334);
    const auto *smk0_336 = buffer.data(smk0 + 336);
    const auto *smk0_338 = buffer.data(smk0 + 338);
    const auto *smk0_339 = buffer.data(smk0 + 339);
    const auto *smk0_341 = buffer.data(smk0 + 341);
    const auto *smk0_342 = buffer.data(smk0 + 342);
    const auto *smk0_344 = buffer.data(smk0 + 344);
    const auto *smk0_359 = buffer.data(smk0 + 359);

    const auto *smi_217 = buffer.data(smi + 217);
    const auto *smi_223 = buffer.data(smi + 223);
    const auto *smi_224 = buffer.data(smi + 224);
    const auto *smi_227 = buffer.data(smi + 227);
    const auto *smi_230 = buffer.data(smi + 230);
    const auto *smi_234 = buffer.data(smi + 234);
    const auto *smi_245 = buffer.data(smi + 245);
    const auto *smi_247 = buffer.data(smi + 247);
    const auto *smi_248 = buffer.data(smi + 248);
    const auto *smi_249 = buffer.data(smi + 249);
    const auto *smi_250 = buffer.data(smi + 250);
    const auto *smi_251 = buffer.data(smi + 251);
    const auto *smi_252 = buffer.data(smi + 252);
    const auto *smi_253 = buffer.data(smi + 253);
    const auto *smi_254 = buffer.data(smi + 254);
    const auto *smi_255 = buffer.data(smi + 255);
    const auto *smi_257 = buffer.data(smi + 257);
    const auto *smi_258 = buffer.data(smi + 258);
    const auto *smi_260 = buffer.data(smi + 260);
    const auto *smi_261 = buffer.data(smi + 261);
    const auto *smi_262 = buffer.data(smi + 262);
    const auto *smi_264 = buffer.data(smi + 264);
    const auto *smi_265 = buffer.data(smi + 265);
    const auto *smi_266 = buffer.data(smi + 266);
    const auto *smi_273 = buffer.data(smi + 273);
    const auto *smi_275 = buffer.data(smi + 275);
    const auto *smi_276 = buffer.data(smi + 276);
    const auto *smi_277 = buffer.data(smi + 277);
    const auto *smi_278 = buffer.data(smi + 278);
    const auto *smi_279 = buffer.data(smi + 279);
    const auto *smi_280 = buffer.data(smi + 280);
    const auto *smi_282 = buffer.data(smi + 282);
    const auto *smi_285 = buffer.data(smi + 285);
    const auto *smi_289 = buffer.data(smi + 289);
    const auto *smi_294 = buffer.data(smi + 294);
    const auto *smi_301 = buffer.data(smi + 301);
    const auto *smi_303 = buffer.data(smi + 303);
    const auto *smi_360 = buffer.data(smi + 360);
    const auto *smi_361 = buffer.data(smi + 361);
    const auto *smi_362 = buffer.data(smi + 362);
    const auto *smi_363 = buffer.data(smi + 363);
    const auto *smi_385 = buffer.data(smi + 385);
    const auto *smi_386 = buffer.data(smi + 386);
    const auto *smi_387 = buffer.data(smi + 387);
    const auto *smi_388 = buffer.data(smi + 388);
    const auto *smi_389 = buffer.data(smi + 389);
    const auto *smi_390 = buffer.data(smi + 390);
    const auto *smi_391 = buffer.data(smi + 391);
    const auto *smi_392 = buffer.data(smi + 392);
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

    const auto *smk1_324 = buffer.data(smk1 + 324);
    const auto *smk1_327 = buffer.data(smk1 + 327);
    const auto *smk1_329 = buffer.data(smk1 + 329);
    const auto *smk1_330 = buffer.data(smk1 + 330);
    const auto *smk1_333 = buffer.data(smk1 + 333);
    const auto *smk1_334 = buffer.data(smk1 + 334);
    const auto *smk1_336 = buffer.data(smk1 + 336);
    const auto *smk1_338 = buffer.data(smk1 + 338);
    const auto *smk1_339 = buffer.data(smk1 + 339);
    const auto *smk1_341 = buffer.data(smk1 + 341);
    const auto *smk1_342 = buffer.data(smk1 + 342);
    const auto *smk1_344 = buffer.data(smk1 + 344);
    const auto *smk1_359 = buffer.data(smk1 + 359);

    const auto *snh0_267 = buffer.data(snh0 + 267);
    const auto *snh0_269 = buffer.data(snh0 + 269);
    const auto *snh0_270 = buffer.data(snh0 + 270);
    const auto *snh0_271 = buffer.data(snh0 + 271);
    const auto *snh0_272 = buffer.data(snh0 + 272);
    const auto *snh0_288 = buffer.data(snh0 + 288);
    const auto *snh0_290 = buffer.data(snh0 + 290);
    const auto *snh0_291 = buffer.data(snh0 + 291);
    const auto *snh0_292 = buffer.data(snh0 + 292);
    const auto *snh0_293 = buffer.data(snh0 + 293);
    const auto *snh0_294 = buffer.data(snh0 + 294);
    const auto *snh0_297 = buffer.data(snh0 + 297);
    const auto *snh0_299 = buffer.data(snh0 + 299);
    const auto *snh0_300 = buffer.data(snh0 + 300);
    const auto *snh0_303 = buffer.data(snh0 + 303);
    const auto *snh0_304 = buffer.data(snh0 + 304);
    const auto *snh0_306 = buffer.data(snh0 + 306);
    const auto *snh0_308 = buffer.data(snh0 + 308);
    const auto *snh0_309 = buffer.data(snh0 + 309);
    const auto *snh0_311 = buffer.data(snh0 + 311);
    const auto *snh0_312 = buffer.data(snh0 + 312);
    const auto *snh0_313 = buffer.data(snh0 + 313);
    const auto *snh0_314 = buffer.data(snh0 + 314);
    const auto *snh0_315 = buffer.data(snh0 + 315);
    const auto *snh0_318 = buffer.data(snh0 + 318);
    const auto *snh0_320 = buffer.data(snh0 + 320);
    const auto *snh0_321 = buffer.data(snh0 + 321);
    const auto *snh0_324 = buffer.data(snh0 + 324);
    const auto *snh0_325 = buffer.data(snh0 + 325);
    const auto *snh0_327 = buffer.data(snh0 + 327);
    const auto *snh0_329 = buffer.data(snh0 + 329);
    const auto *snh0_330 = buffer.data(snh0 + 330);
    const auto *snh0_332 = buffer.data(snh0 + 332);
    const auto *snh0_333 = buffer.data(snh0 + 333);
    const auto *snh0_335 = buffer.data(snh0 + 335);

    const auto *snh1_267 = buffer.data(snh1 + 267);
    const auto *snh1_269 = buffer.data(snh1 + 269);
    const auto *snh1_270 = buffer.data(snh1 + 270);
    const auto *snh1_271 = buffer.data(snh1 + 271);
    const auto *snh1_272 = buffer.data(snh1 + 272);
    const auto *snh1_288 = buffer.data(snh1 + 288);
    const auto *snh1_290 = buffer.data(snh1 + 290);
    const auto *snh1_291 = buffer.data(snh1 + 291);
    const auto *snh1_292 = buffer.data(snh1 + 292);
    const auto *snh1_293 = buffer.data(snh1 + 293);
    const auto *snh1_294 = buffer.data(snh1 + 294);
    const auto *snh1_297 = buffer.data(snh1 + 297);
    const auto *snh1_299 = buffer.data(snh1 + 299);
    const auto *snh1_300 = buffer.data(snh1 + 300);
    const auto *snh1_303 = buffer.data(snh1 + 303);
    const auto *snh1_304 = buffer.data(snh1 + 304);
    const auto *snh1_306 = buffer.data(snh1 + 306);
    const auto *snh1_308 = buffer.data(snh1 + 308);
    const auto *snh1_309 = buffer.data(snh1 + 309);
    const auto *snh1_311 = buffer.data(snh1 + 311);
    const auto *snh1_312 = buffer.data(snh1 + 312);
    const auto *snh1_313 = buffer.data(snh1 + 313);
    const auto *snh1_314 = buffer.data(snh1 + 314);
    const auto *snh1_315 = buffer.data(snh1 + 315);
    const auto *snh1_318 = buffer.data(snh1 + 318);
    const auto *snh1_320 = buffer.data(snh1 + 320);
    const auto *snh1_321 = buffer.data(snh1 + 321);
    const auto *snh1_324 = buffer.data(snh1 + 324);
    const auto *snh1_325 = buffer.data(snh1 + 325);
    const auto *snh1_327 = buffer.data(snh1 + 327);
    const auto *snh1_329 = buffer.data(snh1 + 329);
    const auto *snh1_330 = buffer.data(snh1 + 330);
    const auto *snh1_332 = buffer.data(snh1 + 332);
    const auto *snh1_333 = buffer.data(snh1 + 333);
    const auto *snh1_335 = buffer.data(snh1 + 335);

    const auto *sni_357 = buffer.data(sni + 357);
    const auto *sni_359 = buffer.data(sni + 359);
    const auto *sni_360 = buffer.data(sni + 360);
    const auto *sni_361 = buffer.data(sni + 361);
    const auto *sni_362 = buffer.data(sni + 362);
    const auto *sni_363 = buffer.data(sni + 363);
    const auto *sni_364 = buffer.data(sni + 364);
    const auto *sni_366 = buffer.data(sni + 366);
    const auto *sni_367 = buffer.data(sni + 367);
    const auto *sni_369 = buffer.data(sni + 369);
    const auto *sni_370 = buffer.data(sni + 370);
    const auto *sni_373 = buffer.data(sni + 373);
    const auto *sni_374 = buffer.data(sni + 374);
    const auto *sni_378 = buffer.data(sni + 378);
    const auto *sni_385 = buffer.data(sni + 385);
    const auto *sni_386 = buffer.data(sni + 386);
    const auto *sni_387 = buffer.data(sni + 387);
    const auto *sni_388 = buffer.data(sni + 388);
    const auto *sni_389 = buffer.data(sni + 389);
    const auto *sni_390 = buffer.data(sni + 390);
    const auto *sni_391 = buffer.data(sni + 391);
    const auto *sni_392 = buffer.data(sni + 392);
    const auto *sni_394 = buffer.data(sni + 394);
    const auto *sni_395 = buffer.data(sni + 395);
    const auto *sni_397 = buffer.data(sni + 397);
    const auto *sni_398 = buffer.data(sni + 398);
    const auto *sni_401 = buffer.data(sni + 401);
    const auto *sni_402 = buffer.data(sni + 402);
    const auto *sni_404 = buffer.data(sni + 404);
    const auto *sni_406 = buffer.data(sni + 406);
    const auto *sni_407 = buffer.data(sni + 407);
    const auto *sni_409 = buffer.data(sni + 409);
    const auto *sni_410 = buffer.data(sni + 410);
    const auto *sni_412 = buffer.data(sni + 412);
    const auto *sni_413 = buffer.data(sni + 413);
    const auto *sni_414 = buffer.data(sni + 414);
    const auto *sni_415 = buffer.data(sni + 415);
    const auto *sni_416 = buffer.data(sni + 416);
    const auto *sni_417 = buffer.data(sni + 417);
    const auto *sni_418 = buffer.data(sni + 418);
    const auto *sni_419 = buffer.data(sni + 419);
    const auto *sni_420 = buffer.data(sni + 420);
    const auto *sni_422 = buffer.data(sni + 422);
    const auto *sni_423 = buffer.data(sni + 423);
    const auto *sni_425 = buffer.data(sni + 425);
    const auto *sni_426 = buffer.data(sni + 426);
    const auto *sni_429 = buffer.data(sni + 429);
    const auto *sni_430 = buffer.data(sni + 430);
    const auto *sni_432 = buffer.data(sni + 432);
    const auto *sni_434 = buffer.data(sni + 434);
    const auto *sni_435 = buffer.data(sni + 435);
    const auto *sni_437 = buffer.data(sni + 437);
    const auto *sni_438 = buffer.data(sni + 438);
    const auto *sni_440 = buffer.data(sni + 440);
    const auto *sni_441 = buffer.data(sni + 441);
    const auto *sni_442 = buffer.data(sni + 442);
    const auto *sni_443 = buffer.data(sni + 443);
    const auto *sni_444 = buffer.data(sni + 444);
    const auto *sni_445 = buffer.data(sni + 445);
    const auto *sni_446 = buffer.data(sni + 446);
    const auto *sni_447 = buffer.data(sni + 447);

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, smi_360, smi_361, smi_362, smi_363, \
                         sni_360, sni_361, sni_362, sni_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_21 * smi_360[k]
                   + f_3 * pc_x[k] * sni_360[k];

        t_457[k] = f_21 * smi_361[k]
                   + f_3 * pc_x[k] * sni_361[k];

        t_458[k] = f_21 * smi_362[k]
                   + f_3 * pc_x[k] * sni_362[k];

        t_459[k] = f_21 * smi_363[k]
                   + f_3 * pc_x[k] * sni_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, smi_217, smi_245, smi_247, snh0_267, \
                         snh0_269, snh1_267, snh1_269, sni_357, \
                         sni_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * smi_245[k]
                   + f_1 * snh0_267[k]
                   - f_2 * snh1_267[k]
                   + f_3 * pc_y[k] * sni_357[k];

        t_461[k] = f_14 * smi_217[k]
                   + f_3 * pc_z[k] * sni_357[k];

        t_462[k] = f_14 * smi_247[k]
                   + f_4 * snh0_269[k]
                   - f_5 * snh1_269[k]
                   + f_3 * pc_y[k] * sni_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_y, smi_248, smi_249, smi_250, snh0_270, \
                         snh0_271, snh0_272, snh1_270, snh1_271, snh1_272, sni_360, sni_361, \
                         sni_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_14 * smi_248[k]
                   + f_6 * snh0_270[k]
                   - f_7 * snh1_270[k]
                   + f_3 * pc_y[k] * sni_360[k];

        t_464[k] = f_14 * smi_249[k]
                   + f_8 * snh0_271[k]
                   - f_9 * snh1_271[k]
                   + f_3 * pc_y[k] * sni_361[k];

        t_465[k] = f_14 * smi_250[k]
                   + f_10 * snh0_272[k]
                   - f_11 * snh1_272[k]
                   + f_3 * pc_y[k] * sni_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_y, pc_y, pc_z, smk0_324, smi_223, \
                         smi_251, smi_252, smk1_324, snh0_272, snh1_272, sni_363, \
                         sni_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * smi_251[k]
                   + f_3 * pc_y[k] * sni_363[k];

        t_467[k] = f_14 * smi_223[k]
                   + f_1 * snh0_272[k]
                   - f_2 * snh1_272[k]
                   + f_3 * pc_z[k] * sni_363[k];

        t_468[k] = pb_y[k] * smk0_324[k]
                   - f_12 * pc_y[k] * smk1_324[k];

        t_469[k] = f_13 * smi_252[k]
                   + f_3 * pc_y[k] * sni_364[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pb_y, pc_y, pc_z, smk0_327, smk0_329, \
                         smi_224, smi_253, smi_254, smk1_327, smk1_329, sni_364, \
                         sni_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_15 * smi_224[k]
                   + f_3 * pc_z[k] * sni_364[k];

        t_471[k] = pb_y[k] * smk0_327[k]
                   + f_14 * smi_253[k]
                   - f_12 * pc_y[k] * smk1_327[k];

        t_472[k] = f_13 * smi_254[k]
                   + f_3 * pc_y[k] * sni_366[k];

        t_473[k] = pb_y[k] * smk0_329[k]
                   - f_12 * pc_y[k] * smk1_329[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pc_y, pc_z, smk0_330, smk0_333, \
                         smi_227, smi_255, smi_257, smk1_330, smk1_333, sni_367, \
                         sni_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = pb_y[k] * smk0_330[k]
                   + f_15 * smi_255[k]
                   - f_12 * pc_y[k] * smk1_330[k];

        t_475[k] = f_15 * smi_227[k]
                   + f_3 * pc_z[k] * sni_367[k];

        t_476[k] = f_13 * smi_257[k]
                   + f_3 * pc_y[k] * sni_369[k];

        t_477[k] = pb_y[k] * smk0_333[k]
                   - f_12 * pc_y[k] * smk1_333[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pb_y, pc_y, pc_z, smk0_334, smk0_336, smi_230, \
                         smi_258, smi_260, smk1_334, smk1_336, \
                         sni_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = pb_y[k] * smk0_334[k]
                   + f_16 * smi_258[k]
                   - f_12 * pc_y[k] * smk1_334[k];

        t_479[k] = f_15 * smi_230[k]
                   + f_3 * pc_z[k] * sni_370[k];

        t_480[k] = pb_y[k] * smk0_336[k]
                   + f_14 * smi_260[k]
                   - f_12 * pc_y[k] * smk1_336[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pb_y, pc_y, pc_z, smk0_338, smk0_339, \
                         smi_234, smi_261, smi_262, smk1_338, smk1_339, sni_373, \
                         sni_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_13 * smi_261[k]
                   + f_3 * pc_y[k] * sni_373[k];

        t_482[k] = pb_y[k] * smk0_338[k]
                   - f_12 * pc_y[k] * smk1_338[k];

        t_483[k] = pb_y[k] * smk0_339[k]
                   + f_17 * smi_262[k]
                   - f_12 * pc_y[k] * smk1_339[k];

        t_484[k] = f_15 * smi_234[k]
                   + f_3 * pc_z[k] * sni_374[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pb_y, pc_y, smk0_341, smk0_342, smk0_344, \
                         smi_264, smi_265, smi_266, smk1_341, smk1_342, smk1_344, \
                         sni_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pb_y[k] * smk0_341[k]
                   + f_15 * smi_264[k]
                   - f_12 * pc_y[k] * smk1_341[k];

        t_486[k] = pb_y[k] * smk0_342[k]
                   + f_14 * smi_265[k]
                   - f_12 * pc_y[k] * smk1_342[k];

        t_487[k] = f_13 * smi_266[k]
                   + f_3 * pc_y[k] * sni_378[k];

        t_488[k] = pb_y[k] * smk0_344[k]
                   - f_12 * pc_y[k] * smk1_344[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, smi_385, smi_386, smi_387, \
                         smi_388, smi_389, sni_385, sni_386, sni_387, sni_388, \
                         sni_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_21 * smi_385[k]
                   + f_3 * pc_x[k] * sni_385[k];

        t_490[k] = f_21 * smi_386[k]
                   + f_3 * pc_x[k] * sni_386[k];

        t_491[k] = f_21 * smi_387[k]
                   + f_3 * pc_x[k] * sni_387[k];

        t_492[k] = f_21 * smi_388[k]
                   + f_3 * pc_x[k] * sni_388[k];

        t_493[k] = f_21 * smi_389[k]
                   + f_3 * pc_x[k] * sni_389[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pc_y, pc_z, smi_245, smi_273, \
                         smi_390, smi_391, snh0_288, snh1_288, sni_385, sni_390, \
                         sni_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_21 * smi_390[k]
                   + f_3 * pc_x[k] * sni_390[k];

        t_495[k] = f_21 * smi_391[k]
                   + f_3 * pc_x[k] * sni_391[k];

        t_496[k] = f_13 * smi_273[k]
                   + f_1 * snh0_288[k]
                   - f_2 * snh1_288[k]
                   + f_3 * pc_y[k] * sni_385[k];

        t_497[k] = f_15 * smi_245[k]
                   + f_3 * pc_z[k] * sni_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, smi_275, smi_276, smi_277, snh0_290, \
                         snh0_291, snh0_292, snh1_290, snh1_291, snh1_292, sni_387, sni_388, \
                         sni_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_13 * smi_275[k]
                   + f_4 * snh0_290[k]
                   - f_5 * snh1_290[k]
                   + f_3 * pc_y[k] * sni_387[k];

        t_499[k] = f_13 * smi_276[k]
                   + f_6 * snh0_291[k]
                   - f_7 * snh1_291[k]
                   + f_3 * pc_y[k] * sni_388[k];

        t_500[k] = f_13 * smi_277[k]
                   + f_8 * snh0_292[k]
                   - f_9 * snh1_292[k]
                   + f_3 * pc_y[k] * sni_389[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pb_y, pc_y, smk0_359, smi_278, smi_279, \
                         smk1_359, snh0_293, snh1_293, sni_390, \
                         sni_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * smi_278[k]
                   + f_10 * snh0_293[k]
                   - f_11 * snh1_293[k]
                   + f_3 * pc_y[k] * sni_390[k];

        t_502[k] = f_13 * smi_279[k]
                   + f_3 * pc_y[k] * sni_391[k];

        t_503[k] = pb_y[k] * smk0_359[k]
                   - f_12 * pc_y[k] * smk1_359[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, smi_252, smi_392, \
                         smi_395, snh0_294, snh0_297, snh1_294, snh1_297, sni_392, \
                         sni_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_21 * smi_392[k]
                   + f_1 * snh0_294[k]
                   - f_2 * snh1_294[k]
                   + f_3 * pc_x[k] * sni_392[k];

        t_505[k] = f_3 * pc_y[k] * sni_392[k];

        t_506[k] = f_16 * smi_252[k]
                   + f_3 * pc_z[k] * sni_392[k];

        t_507[k] = f_21 * smi_395[k]
                   + f_4 * snh0_297[k]
                   - f_5 * snh1_297[k]
                   + f_3 * pc_x[k] * sni_395[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, pc_y, smi_397, smi_398, snh0_299, \
                         snh0_300, snh1_299, snh1_300, sni_394, sni_397, \
                         sni_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_3 * pc_y[k] * sni_394[k];

        t_509[k] = f_21 * smi_397[k]
                   + f_4 * snh0_299[k]
                   - f_5 * snh1_299[k]
                   + f_3 * pc_x[k] * sni_397[k];

        t_510[k] = f_21 * smi_398[k]
                   + f_6 * snh0_300[k]
                   - f_7 * snh1_300[k]
                   + f_3 * pc_x[k] * sni_398[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pc_x, pc_y, pc_z, smi_255, smi_401, snh0_303, \
                         snh1_303, sni_395, sni_397, sni_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_16 * smi_255[k]
                   + f_3 * pc_z[k] * sni_395[k];

        t_512[k] = f_3 * pc_y[k] * sni_397[k];

        t_513[k] = f_21 * smi_401[k]
                   + f_6 * snh0_303[k]
                   - f_7 * snh1_303[k]
                   + f_3 * pc_x[k] * sni_401[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pc_x, pc_z, smi_258, smi_402, smi_404, snh0_304, \
                         snh0_306, snh1_304, snh1_306, sni_398, sni_402, \
                         sni_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_21 * smi_402[k]
                   + f_8 * snh0_304[k]
                   - f_9 * snh1_304[k]
                   + f_3 * pc_x[k] * sni_402[k];

        t_515[k] = f_16 * smi_258[k]
                   + f_3 * pc_z[k] * sni_398[k];

        t_516[k] = f_21 * smi_404[k]
                   + f_8 * snh0_306[k]
                   - f_9 * snh1_306[k]
                   + f_3 * pc_x[k] * sni_404[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pc_x, pc_y, smi_406, smi_407, snh0_308, \
                         snh0_309, snh1_308, snh1_309, sni_401, sni_406, \
                         sni_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_3 * pc_y[k] * sni_401[k];

        t_518[k] = f_21 * smi_406[k]
                   + f_8 * snh0_308[k]
                   - f_9 * snh1_308[k]
                   + f_3 * pc_x[k] * sni_406[k];

        t_519[k] = f_21 * smi_407[k]
                   + f_10 * snh0_309[k]
                   - f_11 * snh1_309[k]
                   + f_3 * pc_x[k] * sni_407[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_x, pc_z, smi_262, smi_409, smi_410, snh0_311, \
                         snh0_312, snh1_311, snh1_312, sni_402, sni_409, \
                         sni_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_16 * smi_262[k]
                   + f_3 * pc_z[k] * sni_402[k];

        t_521[k] = f_21 * smi_409[k]
                   + f_10 * snh0_311[k]
                   - f_11 * snh1_311[k]
                   + f_3 * pc_x[k] * sni_409[k];

        t_522[k] = f_21 * smi_410[k]
                   + f_10 * snh0_312[k]
                   - f_11 * snh1_312[k]
                   + f_3 * pc_x[k] * sni_410[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pc_x, pc_y, smi_412, smi_413, smi_414, \
                         snh0_314, snh1_314, sni_406, sni_412, sni_413, \
                         sni_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_3 * pc_y[k] * sni_406[k];

        t_524[k] = f_21 * smi_412[k]
                   + f_10 * snh0_314[k]
                   - f_11 * snh1_314[k]
                   + f_3 * pc_x[k] * sni_412[k];

        t_525[k] = f_21 * smi_413[k]
                   + f_3 * pc_x[k] * sni_413[k];

        t_526[k] = f_21 * smi_414[k]
                   + f_3 * pc_x[k] * sni_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, pc_x, smi_415, smi_416, smi_417, \
                         smi_418, smi_419, sni_415, sni_416, sni_417, sni_418, \
                         sni_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_21 * smi_415[k]
                   + f_3 * pc_x[k] * sni_415[k];

        t_528[k] = f_21 * smi_416[k]
                   + f_3 * pc_x[k] * sni_416[k];

        t_529[k] = f_21 * smi_417[k]
                   + f_3 * pc_x[k] * sni_417[k];

        t_530[k] = f_21 * smi_418[k]
                   + f_3 * pc_x[k] * sni_418[k];

        t_531[k] = f_21 * smi_419[k]
                   + f_3 * pc_x[k] * sni_419[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pc_y, pc_z, smi_273, snh0_309, snh0_311, \
                         snh0_312, snh1_309, snh1_311, snh1_312, sni_413, sni_415, \
                         sni_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_1 * snh0_309[k]
                   - f_2 * snh1_309[k]
                   + f_3 * pc_y[k] * sni_413[k];

        t_533[k] = f_16 * smi_273[k]
                   + f_3 * pc_z[k] * sni_413[k];

        t_534[k] = f_4 * snh0_311[k]
                   - f_5 * snh1_311[k]
                   + f_3 * pc_y[k] * sni_415[k];

        t_535[k] = f_6 * snh0_312[k]
                   - f_7 * snh1_312[k]
                   + f_3 * pc_y[k] * sni_416[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pc_y, pc_z, smi_279, snh0_313, snh0_314, \
                         snh1_313, snh1_314, sni_417, sni_418, \
                         sni_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_8 * snh0_313[k]
                   - f_9 * snh1_313[k]
                   + f_3 * pc_y[k] * sni_417[k];

        t_537[k] = f_10 * snh0_314[k]
                   - f_11 * snh1_314[k]
                   + f_3 * pc_y[k] * sni_418[k];

        t_538[k] = f_3 * pc_y[k] * sni_419[k];

        t_539[k] = f_16 * smi_279[k]
                   + f_1 * snh0_314[k]
                   - f_2 * snh1_314[k]
                   + f_3 * pc_z[k] * sni_419[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pc_x, pc_y, pc_z, smi_280, smi_420, \
                         smi_423, snh0_315, snh0_318, snh1_315, snh1_318, sni_420, \
                         sni_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_17 * smi_420[k]
                   + f_1 * snh0_315[k]
                   - f_2 * snh1_315[k]
                   + f_3 * pc_x[k] * sni_420[k];

        t_541[k] = f_17 * smi_280[k]
                   + f_3 * pc_y[k] * sni_420[k];

        t_542[k] = f_3 * pc_z[k] * sni_420[k];

        t_543[k] = f_17 * smi_423[k]
                   + f_4 * snh0_318[k]
                   - f_5 * snh1_318[k]
                   + f_3 * pc_x[k] * sni_423[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, pc_x, pc_y, smi_282, smi_425, smi_426, snh0_320, \
                         snh0_321, snh1_320, snh1_321, sni_422, sni_425, \
                         sni_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_17 * smi_282[k]
                   + f_3 * pc_y[k] * sni_422[k];

        t_545[k] = f_17 * smi_425[k]
                   + f_4 * snh0_320[k]
                   - f_5 * snh1_320[k]
                   + f_3 * pc_x[k] * sni_425[k];

        t_546[k] = f_17 * smi_426[k]
                   + f_6 * snh0_321[k]
                   - f_7 * snh1_321[k]
                   + f_3 * pc_x[k] * sni_426[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, pc_x, pc_y, pc_z, smi_285, smi_429, snh0_324, \
                         snh1_324, sni_423, sni_425, sni_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_3 * pc_z[k] * sni_423[k];

        t_548[k] = f_17 * smi_285[k]
                   + f_3 * pc_y[k] * sni_425[k];

        t_549[k] = f_17 * smi_429[k]
                   + f_6 * snh0_324[k]
                   - f_7 * snh1_324[k]
                   + f_3 * pc_x[k] * sni_429[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, pc_x, pc_z, smi_430, smi_432, snh0_325, \
                         snh0_327, snh1_325, snh1_327, sni_426, sni_430, \
                         sni_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_17 * smi_430[k]
                   + f_8 * snh0_325[k]
                   - f_9 * snh1_325[k]
                   + f_3 * pc_x[k] * sni_430[k];

        t_551[k] = f_3 * pc_z[k] * sni_426[k];

        t_552[k] = f_17 * smi_432[k]
                   + f_8 * snh0_327[k]
                   - f_9 * snh1_327[k]
                   + f_3 * pc_x[k] * sni_432[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_x, pc_y, smi_289, smi_434, smi_435, snh0_329, \
                         snh0_330, snh1_329, snh1_330, sni_429, sni_434, \
                         sni_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_17 * smi_289[k]
                   + f_3 * pc_y[k] * sni_429[k];

        t_554[k] = f_17 * smi_434[k]
                   + f_8 * snh0_329[k]
                   - f_9 * snh1_329[k]
                   + f_3 * pc_x[k] * sni_434[k];

        t_555[k] = f_17 * smi_435[k]
                   + f_10 * snh0_330[k]
                   - f_11 * snh1_330[k]
                   + f_3 * pc_x[k] * sni_435[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_x, pc_z, smi_437, smi_438, snh0_332, \
                         snh0_333, snh1_332, snh1_333, sni_430, sni_437, \
                         sni_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_3 * pc_z[k] * sni_430[k];

        t_557[k] = f_17 * smi_437[k]
                   + f_10 * snh0_332[k]
                   - f_11 * snh1_332[k]
                   + f_3 * pc_x[k] * sni_437[k];

        t_558[k] = f_17 * smi_438[k]
                   + f_10 * snh0_333[k]
                   - f_11 * snh1_333[k]
                   + f_3 * pc_x[k] * sni_438[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pc_x, pc_y, smi_294, smi_440, smi_441, \
                         smi_442, snh0_335, snh1_335, sni_434, sni_440, sni_441, \
                         sni_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_17 * smi_294[k]
                   + f_3 * pc_y[k] * sni_434[k];

        t_560[k] = f_17 * smi_440[k]
                   + f_10 * snh0_335[k]
                   - f_11 * snh1_335[k]
                   + f_3 * pc_x[k] * sni_440[k];

        t_561[k] = f_17 * smi_441[k]
                   + f_3 * pc_x[k] * sni_441[k];

        t_562[k] = f_17 * smi_442[k]
                   + f_3 * pc_x[k] * sni_442[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, pc_x, smi_443, smi_444, smi_445, \
                         smi_446, smi_447, sni_443, sni_444, sni_445, sni_446, \
                         sni_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_17 * smi_443[k]
                   + f_3 * pc_x[k] * sni_443[k];

        t_564[k] = f_17 * smi_444[k]
                   + f_3 * pc_x[k] * sni_444[k];

        t_565[k] = f_17 * smi_445[k]
                   + f_3 * pc_x[k] * sni_445[k];

        t_566[k] = f_17 * smi_446[k]
                   + f_3 * pc_x[k] * sni_446[k];

        t_567[k] = f_17 * smi_447[k]
                   + f_3 * pc_x[k] * sni_447[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pc_y, pc_z, smi_301, smi_303, snh0_330, \
                         snh0_332, snh1_330, snh1_332, sni_441, \
                         sni_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_17 * smi_301[k]
                   + f_1 * snh0_330[k]
                   - f_2 * snh1_330[k]
                   + f_3 * pc_y[k] * sni_441[k];

        t_569[k] = f_3 * pc_z[k] * sni_441[k];

        t_570[k] = f_17 * smi_303[k]
                   + f_4 * snh0_332[k]
                   - f_5 * snh1_332[k]
                   + f_3 * pc_y[k] * sni_443[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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

    const auto *smk0_360 = buffer.data(smk0 + 360);
    const auto *smk0_363 = buffer.data(smk0 + 363);
    const auto *smk0_366 = buffer.data(smk0 + 366);
    const auto *smk0_370 = buffer.data(smk0 + 370);
    const auto *smk0_372 = buffer.data(smk0 + 372);
    const auto *smk0_375 = buffer.data(smk0 + 375);
    const auto *smk0_377 = buffer.data(smk0 + 377);
    const auto *smk0_378 = buffer.data(smk0 + 378);
    const auto *smk0_388 = buffer.data(smk0 + 388);

    const auto *smi_280 = buffer.data(smi + 280);
    const auto *smi_283 = buffer.data(smi + 283);
    const auto *smi_286 = buffer.data(smi + 286);
    const auto *smi_287 = buffer.data(smi + 287);
    const auto *smi_290 = buffer.data(smi + 290);
    const auto *smi_291 = buffer.data(smi + 291);
    const auto *smi_292 = buffer.data(smi + 292);
    const auto *smi_301 = buffer.data(smi + 301);
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
    const auto *smi_329 = buffer.data(smi + 329);
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
    const auto *smi_350 = buffer.data(smi + 350);
    const auto *smi_357 = buffer.data(smi + 357);
    const auto *smi_359 = buffer.data(smi + 359);
    const auto *smi_360 = buffer.data(smi + 360);
    const auto *smi_361 = buffer.data(smi + 361);
    const auto *smi_362 = buffer.data(smi + 362);
    const auto *smi_363 = buffer.data(smi + 363);
    const auto *smi_364 = buffer.data(smi + 364);
    const auto *smi_366 = buffer.data(smi + 366);
    const auto *smi_369 = buffer.data(smi + 369);
    const auto *smi_373 = buffer.data(smi + 373);
    const auto *smi_378 = buffer.data(smi + 378);
    const auto *smi_385 = buffer.data(smi + 385);
    const auto *smi_387 = buffer.data(smi + 387);
    const auto *smi_453 = buffer.data(smi + 453);
    const auto *smi_457 = buffer.data(smi + 457);
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

    const auto *smk1_360 = buffer.data(smk1 + 360);
    const auto *smk1_363 = buffer.data(smk1 + 363);
    const auto *smk1_366 = buffer.data(smk1 + 366);
    const auto *smk1_370 = buffer.data(smk1 + 370);
    const auto *smk1_372 = buffer.data(smk1 + 372);
    const auto *smk1_375 = buffer.data(smk1 + 375);
    const auto *smk1_377 = buffer.data(smk1 + 377);
    const auto *smk1_378 = buffer.data(smk1 + 378);
    const auto *smk1_388 = buffer.data(smk1 + 388);

    const auto *snh0_333 = buffer.data(snh0 + 333);
    const auto *snh0_334 = buffer.data(snh0 + 334);
    const auto *snh0_335 = buffer.data(snh0 + 335);
    const auto *snh0_341 = buffer.data(snh0 + 341);
    const auto *snh0_345 = buffer.data(snh0 + 345);
    const auto *snh0_350 = buffer.data(snh0 + 350);
    const auto *snh0_353 = buffer.data(snh0 + 353);
    const auto *snh0_354 = buffer.data(snh0 + 354);
    const auto *snh0_355 = buffer.data(snh0 + 355);
    const auto *snh0_356 = buffer.data(snh0 + 356);
    const auto *snh0_357 = buffer.data(snh0 + 357);
    const auto *snh0_360 = buffer.data(snh0 + 360);
    const auto *snh0_362 = buffer.data(snh0 + 362);
    const auto *snh0_363 = buffer.data(snh0 + 363);
    const auto *snh0_366 = buffer.data(snh0 + 366);
    const auto *snh0_367 = buffer.data(snh0 + 367);
    const auto *snh0_369 = buffer.data(snh0 + 369);
    const auto *snh0_371 = buffer.data(snh0 + 371);
    const auto *snh0_372 = buffer.data(snh0 + 372);
    const auto *snh0_374 = buffer.data(snh0 + 374);
    const auto *snh0_375 = buffer.data(snh0 + 375);
    const auto *snh0_376 = buffer.data(snh0 + 376);
    const auto *snh0_377 = buffer.data(snh0 + 377);
    const auto *snh0_378 = buffer.data(snh0 + 378);
    const auto *snh0_381 = buffer.data(snh0 + 381);
    const auto *snh0_383 = buffer.data(snh0 + 383);
    const auto *snh0_384 = buffer.data(snh0 + 384);
    const auto *snh0_387 = buffer.data(snh0 + 387);
    const auto *snh0_388 = buffer.data(snh0 + 388);
    const auto *snh0_390 = buffer.data(snh0 + 390);
    const auto *snh0_392 = buffer.data(snh0 + 392);
    const auto *snh0_393 = buffer.data(snh0 + 393);
    const auto *snh0_395 = buffer.data(snh0 + 395);
    const auto *snh0_396 = buffer.data(snh0 + 396);
    const auto *snh0_398 = buffer.data(snh0 + 398);

    const auto *snh1_333 = buffer.data(snh1 + 333);
    const auto *snh1_334 = buffer.data(snh1 + 334);
    const auto *snh1_335 = buffer.data(snh1 + 335);
    const auto *snh1_341 = buffer.data(snh1 + 341);
    const auto *snh1_345 = buffer.data(snh1 + 345);
    const auto *snh1_350 = buffer.data(snh1 + 350);
    const auto *snh1_353 = buffer.data(snh1 + 353);
    const auto *snh1_354 = buffer.data(snh1 + 354);
    const auto *snh1_355 = buffer.data(snh1 + 355);
    const auto *snh1_356 = buffer.data(snh1 + 356);
    const auto *snh1_357 = buffer.data(snh1 + 357);
    const auto *snh1_360 = buffer.data(snh1 + 360);
    const auto *snh1_362 = buffer.data(snh1 + 362);
    const auto *snh1_363 = buffer.data(snh1 + 363);
    const auto *snh1_366 = buffer.data(snh1 + 366);
    const auto *snh1_367 = buffer.data(snh1 + 367);
    const auto *snh1_369 = buffer.data(snh1 + 369);
    const auto *snh1_371 = buffer.data(snh1 + 371);
    const auto *snh1_372 = buffer.data(snh1 + 372);
    const auto *snh1_374 = buffer.data(snh1 + 374);
    const auto *snh1_375 = buffer.data(snh1 + 375);
    const auto *snh1_376 = buffer.data(snh1 + 376);
    const auto *snh1_377 = buffer.data(snh1 + 377);
    const auto *snh1_378 = buffer.data(snh1 + 378);
    const auto *snh1_381 = buffer.data(snh1 + 381);
    const auto *snh1_383 = buffer.data(snh1 + 383);
    const auto *snh1_384 = buffer.data(snh1 + 384);
    const auto *snh1_387 = buffer.data(snh1 + 387);
    const auto *snh1_388 = buffer.data(snh1 + 388);
    const auto *snh1_390 = buffer.data(snh1 + 390);
    const auto *snh1_392 = buffer.data(snh1 + 392);
    const auto *snh1_393 = buffer.data(snh1 + 393);
    const auto *snh1_395 = buffer.data(snh1 + 395);
    const auto *snh1_396 = buffer.data(snh1 + 396);
    const auto *snh1_398 = buffer.data(snh1 + 398);

    const auto *sni_444 = buffer.data(sni + 444);
    const auto *sni_445 = buffer.data(sni + 445);
    const auto *sni_446 = buffer.data(sni + 446);
    const auto *sni_447 = buffer.data(sni + 447);
    const auto *sni_448 = buffer.data(sni + 448);
    const auto *sni_450 = buffer.data(sni + 450);
    const auto *sni_451 = buffer.data(sni + 451);
    const auto *sni_453 = buffer.data(sni + 453);
    const auto *sni_454 = buffer.data(sni + 454);
    const auto *sni_457 = buffer.data(sni + 457);
    const auto *sni_458 = buffer.data(sni + 458);
    const auto *sni_462 = buffer.data(sni + 462);
    const auto *sni_468 = buffer.data(sni + 468);
    const auto *sni_469 = buffer.data(sni + 469);
    const auto *sni_470 = buffer.data(sni + 470);
    const auto *sni_471 = buffer.data(sni + 471);
    const auto *sni_472 = buffer.data(sni + 472);
    const auto *sni_473 = buffer.data(sni + 473);
    const auto *sni_474 = buffer.data(sni + 474);
    const auto *sni_475 = buffer.data(sni + 475);
    const auto *sni_476 = buffer.data(sni + 476);
    const auto *sni_478 = buffer.data(sni + 478);
    const auto *sni_479 = buffer.data(sni + 479);
    const auto *sni_481 = buffer.data(sni + 481);
    const auto *sni_482 = buffer.data(sni + 482);
    const auto *sni_485 = buffer.data(sni + 485);
    const auto *sni_486 = buffer.data(sni + 486);
    const auto *sni_488 = buffer.data(sni + 488);
    const auto *sni_490 = buffer.data(sni + 490);
    const auto *sni_491 = buffer.data(sni + 491);
    const auto *sni_493 = buffer.data(sni + 493);
    const auto *sni_494 = buffer.data(sni + 494);
    const auto *sni_496 = buffer.data(sni + 496);
    const auto *sni_497 = buffer.data(sni + 497);
    const auto *sni_498 = buffer.data(sni + 498);
    const auto *sni_499 = buffer.data(sni + 499);
    const auto *sni_500 = buffer.data(sni + 500);
    const auto *sni_501 = buffer.data(sni + 501);
    const auto *sni_502 = buffer.data(sni + 502);
    const auto *sni_503 = buffer.data(sni + 503);
    const auto *sni_504 = buffer.data(sni + 504);
    const auto *sni_506 = buffer.data(sni + 506);
    const auto *sni_507 = buffer.data(sni + 507);
    const auto *sni_509 = buffer.data(sni + 509);
    const auto *sni_510 = buffer.data(sni + 510);
    const auto *sni_513 = buffer.data(sni + 513);
    const auto *sni_514 = buffer.data(sni + 514);
    const auto *sni_516 = buffer.data(sni + 516);
    const auto *sni_518 = buffer.data(sni + 518);
    const auto *sni_519 = buffer.data(sni + 519);
    const auto *sni_521 = buffer.data(sni + 521);
    const auto *sni_522 = buffer.data(sni + 522);
    const auto *sni_524 = buffer.data(sni + 524);
    const auto *sni_525 = buffer.data(sni + 525);
    const auto *sni_526 = buffer.data(sni + 526);
    const auto *sni_527 = buffer.data(sni + 527);
    const auto *sni_528 = buffer.data(sni + 528);
    const auto *sni_529 = buffer.data(sni + 529);
    const auto *sni_530 = buffer.data(sni + 530);
    const auto *sni_531 = buffer.data(sni + 531);

#pragma omp simd aligned(t_571, t_572, t_573, pc_y, smi_304, smi_305, smi_306, snh0_333, \
                         snh0_334, snh0_335, snh1_333, snh1_334, snh1_335, sni_444, sni_445, \
                         sni_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_17 * smi_304[k]
                   + f_6 * snh0_333[k]
                   - f_7 * snh1_333[k]
                   + f_3 * pc_y[k] * sni_444[k];

        t_572[k] = f_17 * smi_305[k]
                   + f_8 * snh0_334[k]
                   - f_9 * snh1_334[k]
                   + f_3 * pc_y[k] * sni_445[k];

        t_573[k] = f_17 * smi_306[k]
                   + f_10 * snh0_335[k]
                   - f_11 * snh1_335[k]
                   + f_3 * pc_y[k] * sni_446[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pb_z, pc_y, pc_z, smk0_360, smi_307, \
                         smi_308, smk1_360, snh0_335, snh1_335, sni_447, \
                         sni_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_17 * smi_307[k]
                   + f_3 * pc_y[k] * sni_447[k];

        t_575[k] = f_1 * snh0_335[k]
                   - f_2 * snh1_335[k]
                   + f_3 * pc_z[k] * sni_447[k];

        t_576[k] = pb_z[k] * smk0_360[k]
                   - f_12 * pc_z[k] * smk1_360[k];

        t_577[k] = f_16 * smi_308[k]
                   + f_3 * pc_y[k] * sni_448[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pb_z, pc_y, pc_z, smk0_363, smi_280, smi_310, \
                         smk1_363, sni_448, sni_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_13 * smi_280[k]
                   + f_3 * pc_z[k] * sni_448[k];

        t_579[k] = pb_z[k] * smk0_363[k]
                   - f_12 * pc_z[k] * smk1_363[k];

        t_580[k] = f_16 * smi_310[k]
                   + f_3 * pc_y[k] * sni_450[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pb_z, pc_x, pc_z, smk0_366, smi_283, smi_453, \
                         smk1_366, snh0_341, snh1_341, sni_451, \
                         sni_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_17 * smi_453[k]
                   + f_4 * snh0_341[k]
                   - f_5 * snh1_341[k]
                   + f_3 * pc_x[k] * sni_453[k];

        t_582[k] = pb_z[k] * smk0_366[k]
                   - f_12 * pc_z[k] * smk1_366[k];

        t_583[k] = f_13 * smi_283[k]
                   + f_3 * pc_z[k] * sni_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pb_z, pc_x, pc_y, pc_z, smk0_370, smi_313, \
                         smi_457, smk1_370, snh0_345, snh1_345, sni_453, \
                         sni_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * smi_313[k]
                   + f_3 * pc_y[k] * sni_453[k];

        t_585[k] = f_17 * smi_457[k]
                   + f_6 * snh0_345[k]
                   - f_7 * snh1_345[k]
                   + f_3 * pc_x[k] * sni_457[k];

        t_586[k] = pb_z[k] * smk0_370[k]
                   - f_12 * pc_z[k] * smk1_370[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, pb_z, pc_y, pc_z, smk0_372, smi_286, smi_287, \
                         smi_317, smk1_372, sni_454, sni_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_13 * smi_286[k]
                   + f_3 * pc_z[k] * sni_454[k];

        t_588[k] = pb_z[k] * smk0_372[k]
                   + f_14 * smi_287[k]
                   - f_12 * pc_z[k] * smk1_372[k];

        t_589[k] = f_16 * smi_317[k]
                   + f_3 * pc_y[k] * sni_457[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, pb_z, pc_x, pc_z, smk0_375, smi_290, smi_462, \
                         smk1_375, snh0_350, snh1_350, sni_458, \
                         sni_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_17 * smi_462[k]
                   + f_8 * snh0_350[k]
                   - f_9 * snh1_350[k]
                   + f_3 * pc_x[k] * sni_462[k];

        t_591[k] = pb_z[k] * smk0_375[k]
                   - f_12 * pc_z[k] * smk1_375[k];

        t_592[k] = f_13 * smi_290[k]
                   + f_3 * pc_z[k] * sni_458[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pb_z, pc_y, pc_z, smk0_377, smk0_378, smi_291, \
                         smi_292, smi_322, smk1_377, smk1_378, \
                         sni_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pb_z[k] * smk0_377[k]
                   + f_14 * smi_291[k]
                   - f_12 * pc_z[k] * smk1_377[k];

        t_594[k] = pb_z[k] * smk0_378[k]
                   + f_15 * smi_292[k]
                   - f_12 * pc_z[k] * smk1_378[k];

        t_595[k] = f_16 * smi_322[k]
                   + f_3 * pc_y[k] * sni_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, smi_468, smi_469, smi_470, smi_471, \
                         snh0_356, snh1_356, sni_468, sni_469, sni_470, \
                         sni_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * smi_468[k]
                   + f_10 * snh0_356[k]
                   - f_11 * snh1_356[k]
                   + f_3 * pc_x[k] * sni_468[k];

        t_597[k] = f_17 * smi_469[k]
                   + f_3 * pc_x[k] * sni_469[k];

        t_598[k] = f_17 * smi_470[k]
                   + f_3 * pc_x[k] * sni_470[k];

        t_599[k] = f_17 * smi_471[k]
                   + f_3 * pc_x[k] * sni_471[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, smi_472, smi_473, smi_474, smi_475, \
                         sni_472, sni_473, sni_474, sni_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_17 * smi_472[k]
                   + f_3 * pc_x[k] * sni_472[k];

        t_601[k] = f_17 * smi_473[k]
                   + f_3 * pc_x[k] * sni_473[k];

        t_602[k] = f_17 * smi_474[k]
                   + f_3 * pc_x[k] * sni_474[k];

        t_603[k] = f_17 * smi_475[k]
                   + f_3 * pc_x[k] * sni_475[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, pb_z, pc_y, pc_z, smk0_388, smi_301, smi_331, \
                         smk1_388, snh0_353, snh1_353, sni_469, \
                         sni_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pb_z[k] * smk0_388[k]
                   - f_12 * pc_z[k] * smk1_388[k];

        t_605[k] = f_13 * smi_301[k]
                   + f_3 * pc_z[k] * sni_469[k];

        t_606[k] = f_16 * smi_331[k]
                   + f_4 * snh0_353[k]
                   - f_5 * snh1_353[k]
                   + f_3 * pc_y[k] * sni_471[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, smi_332, smi_333, smi_334, snh0_354, \
                         snh0_355, snh0_356, snh1_354, snh1_355, snh1_356, sni_472, sni_473, \
                         sni_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_16 * smi_332[k]
                   + f_6 * snh0_354[k]
                   - f_7 * snh1_354[k]
                   + f_3 * pc_y[k] * sni_472[k];

        t_608[k] = f_16 * smi_333[k]
                   + f_8 * snh0_355[k]
                   - f_9 * snh1_355[k]
                   + f_3 * pc_y[k] * sni_473[k];

        t_609[k] = f_16 * smi_334[k]
                   + f_10 * snh0_356[k]
                   - f_11 * snh1_356[k]
                   + f_3 * pc_y[k] * sni_474[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_x, pc_y, pc_z, smi_307, smi_335, smi_476, \
                         snh0_356, snh0_357, snh1_356, snh1_357, sni_475, \
                         sni_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_16 * smi_335[k]
                   + f_3 * pc_y[k] * sni_475[k];

        t_611[k] = f_13 * smi_307[k]
                   + f_1 * snh0_356[k]
                   - f_2 * snh1_356[k]
                   + f_3 * pc_z[k] * sni_475[k];

        t_612[k] = f_17 * smi_476[k]
                   + f_1 * snh0_357[k]
                   - f_2 * snh1_357[k]
                   + f_3 * pc_x[k] * sni_476[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pc_x, pc_y, pc_z, smi_308, smi_336, \
                         smi_338, smi_479, snh0_360, snh1_360, sni_476, sni_478, \
                         sni_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_15 * smi_336[k]
                   + f_3 * pc_y[k] * sni_476[k];

        t_614[k] = f_14 * smi_308[k]
                   + f_3 * pc_z[k] * sni_476[k];

        t_615[k] = f_17 * smi_479[k]
                   + f_4 * snh0_360[k]
                   - f_5 * snh1_360[k]
                   + f_3 * pc_x[k] * sni_479[k];

        t_616[k] = f_15 * smi_338[k]
                   + f_3 * pc_y[k] * sni_478[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, pc_x, pc_z, smi_311, smi_481, smi_482, snh0_362, \
                         snh0_363, snh1_362, snh1_363, sni_479, sni_481, \
                         sni_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_17 * smi_481[k]
                   + f_4 * snh0_362[k]
                   - f_5 * snh1_362[k]
                   + f_3 * pc_x[k] * sni_481[k];

        t_618[k] = f_17 * smi_482[k]
                   + f_6 * snh0_363[k]
                   - f_7 * snh1_363[k]
                   + f_3 * pc_x[k] * sni_482[k];

        t_619[k] = f_14 * smi_311[k]
                   + f_3 * pc_z[k] * sni_479[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pc_x, pc_y, smi_341, smi_485, smi_486, snh0_366, \
                         snh0_367, snh1_366, snh1_367, sni_481, sni_485, \
                         sni_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_15 * smi_341[k]
                   + f_3 * pc_y[k] * sni_481[k];

        t_621[k] = f_17 * smi_485[k]
                   + f_6 * snh0_366[k]
                   - f_7 * snh1_366[k]
                   + f_3 * pc_x[k] * sni_485[k];

        t_622[k] = f_17 * smi_486[k]
                   + f_8 * snh0_367[k]
                   - f_9 * snh1_367[k]
                   + f_3 * pc_x[k] * sni_486[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, smi_314, smi_345, smi_488, \
                         snh0_369, snh1_369, sni_482, sni_485, \
                         sni_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_14 * smi_314[k]
                   + f_3 * pc_z[k] * sni_482[k];

        t_624[k] = f_17 * smi_488[k]
                   + f_8 * snh0_369[k]
                   - f_9 * snh1_369[k]
                   + f_3 * pc_x[k] * sni_488[k];

        t_625[k] = f_15 * smi_345[k]
                   + f_3 * pc_y[k] * sni_485[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_x, pc_z, smi_318, smi_490, smi_491, snh0_371, \
                         snh0_372, snh1_371, snh1_372, sni_486, sni_490, \
                         sni_491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_17 * smi_490[k]
                   + f_8 * snh0_371[k]
                   - f_9 * snh1_371[k]
                   + f_3 * pc_x[k] * sni_490[k];

        t_627[k] = f_17 * smi_491[k]
                   + f_10 * snh0_372[k]
                   - f_11 * snh1_372[k]
                   + f_3 * pc_x[k] * sni_491[k];

        t_628[k] = f_14 * smi_318[k]
                   + f_3 * pc_z[k] * sni_486[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_y, smi_350, smi_493, smi_494, snh0_374, \
                         snh0_375, snh1_374, snh1_375, sni_490, sni_493, \
                         sni_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_17 * smi_493[k]
                   + f_10 * snh0_374[k]
                   - f_11 * snh1_374[k]
                   + f_3 * pc_x[k] * sni_493[k];

        t_630[k] = f_17 * smi_494[k]
                   + f_10 * snh0_375[k]
                   - f_11 * snh1_375[k]
                   + f_3 * pc_x[k] * sni_494[k];

        t_631[k] = f_15 * smi_350[k]
                   + f_3 * pc_y[k] * sni_490[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, smi_496, smi_497, smi_498, smi_499, \
                         snh0_377, snh1_377, sni_496, sni_497, sni_498, \
                         sni_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_17 * smi_496[k]
                   + f_10 * snh0_377[k]
                   - f_11 * snh1_377[k]
                   + f_3 * pc_x[k] * sni_496[k];

        t_633[k] = f_17 * smi_497[k]
                   + f_3 * pc_x[k] * sni_497[k];

        t_634[k] = f_17 * smi_498[k]
                   + f_3 * pc_x[k] * sni_498[k];

        t_635[k] = f_17 * smi_499[k]
                   + f_3 * pc_x[k] * sni_499[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, smi_500, smi_501, smi_502, smi_503, \
                         sni_500, sni_501, sni_502, sni_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_17 * smi_500[k]
                   + f_3 * pc_x[k] * sni_500[k];

        t_637[k] = f_17 * smi_501[k]
                   + f_3 * pc_x[k] * sni_501[k];

        t_638[k] = f_17 * smi_502[k]
                   + f_3 * pc_x[k] * sni_502[k];

        t_639[k] = f_17 * smi_503[k]
                   + f_3 * pc_x[k] * sni_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, smi_329, smi_357, smi_359, snh0_372, \
                         snh0_374, snh1_372, snh1_374, sni_497, \
                         sni_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_15 * smi_357[k]
                   + f_1 * snh0_372[k]
                   - f_2 * snh1_372[k]
                   + f_3 * pc_y[k] * sni_497[k];

        t_641[k] = f_14 * smi_329[k]
                   + f_3 * pc_z[k] * sni_497[k];

        t_642[k] = f_15 * smi_359[k]
                   + f_4 * snh0_374[k]
                   - f_5 * snh1_374[k]
                   + f_3 * pc_y[k] * sni_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pc_y, smi_360, smi_361, smi_362, snh0_375, \
                         snh0_376, snh0_377, snh1_375, snh1_376, snh1_377, sni_500, sni_501, \
                         sni_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_15 * smi_360[k]
                   + f_6 * snh0_375[k]
                   - f_7 * snh1_375[k]
                   + f_3 * pc_y[k] * sni_500[k];

        t_644[k] = f_15 * smi_361[k]
                   + f_8 * snh0_376[k]
                   - f_9 * snh1_376[k]
                   + f_3 * pc_y[k] * sni_501[k];

        t_645[k] = f_15 * smi_362[k]
                   + f_10 * snh0_377[k]
                   - f_11 * snh1_377[k]
                   + f_3 * pc_y[k] * sni_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pc_x, pc_y, pc_z, smi_335, smi_363, smi_504, \
                         snh0_377, snh0_378, snh1_377, snh1_378, sni_503, \
                         sni_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_15 * smi_363[k]
                   + f_3 * pc_y[k] * sni_503[k];

        t_647[k] = f_14 * smi_335[k]
                   + f_1 * snh0_377[k]
                   - f_2 * snh1_377[k]
                   + f_3 * pc_z[k] * sni_503[k];

        t_648[k] = f_17 * smi_504[k]
                   + f_1 * snh0_378[k]
                   - f_2 * snh1_378[k]
                   + f_3 * pc_x[k] * sni_504[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pc_x, pc_y, pc_z, smi_336, smi_364, \
                         smi_366, smi_507, snh0_381, snh1_381, sni_504, sni_506, \
                         sni_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_14 * smi_364[k]
                   + f_3 * pc_y[k] * sni_504[k];

        t_650[k] = f_15 * smi_336[k]
                   + f_3 * pc_z[k] * sni_504[k];

        t_651[k] = f_17 * smi_507[k]
                   + f_4 * snh0_381[k]
                   - f_5 * snh1_381[k]
                   + f_3 * pc_x[k] * sni_507[k];

        t_652[k] = f_14 * smi_366[k]
                   + f_3 * pc_y[k] * sni_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, smi_339, smi_509, smi_510, snh0_383, \
                         snh0_384, snh1_383, snh1_384, sni_507, sni_509, \
                         sni_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_17 * smi_509[k]
                   + f_4 * snh0_383[k]
                   - f_5 * snh1_383[k]
                   + f_3 * pc_x[k] * sni_509[k];

        t_654[k] = f_17 * smi_510[k]
                   + f_6 * snh0_384[k]
                   - f_7 * snh1_384[k]
                   + f_3 * pc_x[k] * sni_510[k];

        t_655[k] = f_15 * smi_339[k]
                   + f_3 * pc_z[k] * sni_507[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, smi_369, smi_513, smi_514, snh0_387, \
                         snh0_388, snh1_387, snh1_388, sni_509, sni_513, \
                         sni_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * smi_369[k]
                   + f_3 * pc_y[k] * sni_509[k];

        t_657[k] = f_17 * smi_513[k]
                   + f_6 * snh0_387[k]
                   - f_7 * snh1_387[k]
                   + f_3 * pc_x[k] * sni_513[k];

        t_658[k] = f_17 * smi_514[k]
                   + f_8 * snh0_388[k]
                   - f_9 * snh1_388[k]
                   + f_3 * pc_x[k] * sni_514[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, pc_x, pc_y, pc_z, smi_342, smi_373, smi_516, \
                         snh0_390, snh1_390, sni_510, sni_513, \
                         sni_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_15 * smi_342[k]
                   + f_3 * pc_z[k] * sni_510[k];

        t_660[k] = f_17 * smi_516[k]
                   + f_8 * snh0_390[k]
                   - f_9 * snh1_390[k]
                   + f_3 * pc_x[k] * sni_516[k];

        t_661[k] = f_14 * smi_373[k]
                   + f_3 * pc_y[k] * sni_513[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_z, smi_346, smi_518, smi_519, snh0_392, \
                         snh0_393, snh1_392, snh1_393, sni_514, sni_518, \
                         sni_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_17 * smi_518[k]
                   + f_8 * snh0_392[k]
                   - f_9 * snh1_392[k]
                   + f_3 * pc_x[k] * sni_518[k];

        t_663[k] = f_17 * smi_519[k]
                   + f_10 * snh0_393[k]
                   - f_11 * snh1_393[k]
                   + f_3 * pc_x[k] * sni_519[k];

        t_664[k] = f_15 * smi_346[k]
                   + f_3 * pc_z[k] * sni_514[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, pc_x, pc_y, smi_378, smi_521, smi_522, snh0_395, \
                         snh0_396, snh1_395, snh1_396, sni_518, sni_521, \
                         sni_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_17 * smi_521[k]
                   + f_10 * snh0_395[k]
                   - f_11 * snh1_395[k]
                   + f_3 * pc_x[k] * sni_521[k];

        t_666[k] = f_17 * smi_522[k]
                   + f_10 * snh0_396[k]
                   - f_11 * snh1_396[k]
                   + f_3 * pc_x[k] * sni_522[k];

        t_667[k] = f_14 * smi_378[k]
                   + f_3 * pc_y[k] * sni_518[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pc_x, smi_524, smi_525, smi_526, smi_527, \
                         snh0_398, snh1_398, sni_524, sni_525, sni_526, \
                         sni_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_17 * smi_524[k]
                   + f_10 * snh0_398[k]
                   - f_11 * snh1_398[k]
                   + f_3 * pc_x[k] * sni_524[k];

        t_669[k] = f_17 * smi_525[k]
                   + f_3 * pc_x[k] * sni_525[k];

        t_670[k] = f_17 * smi_526[k]
                   + f_3 * pc_x[k] * sni_526[k];

        t_671[k] = f_17 * smi_527[k]
                   + f_3 * pc_x[k] * sni_527[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, smi_528, smi_529, smi_530, smi_531, \
                         sni_528, sni_529, sni_530, sni_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_17 * smi_528[k]
                   + f_3 * pc_x[k] * sni_528[k];

        t_673[k] = f_17 * smi_529[k]
                   + f_3 * pc_x[k] * sni_529[k];

        t_674[k] = f_17 * smi_530[k]
                   + f_3 * pc_x[k] * sni_530[k];

        t_675[k] = f_17 * smi_531[k]
                   + f_3 * pc_x[k] * sni_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_y, pc_z, smi_357, smi_385, smi_387, snh0_393, \
                         snh0_395, snh1_393, snh1_395, sni_525, \
                         sni_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_14 * smi_385[k]
                   + f_1 * snh0_393[k]
                   - f_2 * snh1_393[k]
                   + f_3 * pc_y[k] * sni_525[k];

        t_677[k] = f_15 * smi_357[k]
                   + f_3 * pc_z[k] * sni_525[k];

        t_678[k] = f_14 * smi_387[k]
                   + f_4 * snh0_395[k]
                   - f_5 * snh1_395[k]
                   + f_3 * pc_y[k] * sni_527[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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
    const auto f_21 = 3.0 / q;

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

    const auto *smk0_504 = buffer.data(smk0 + 504);
    const auto *smk0_507 = buffer.data(smk0 + 507);
    const auto *smk0_509 = buffer.data(smk0 + 509);
    const auto *smk0_510 = buffer.data(smk0 + 510);
    const auto *smk0_513 = buffer.data(smk0 + 513);
    const auto *smk0_514 = buffer.data(smk0 + 514);
    const auto *smk0_516 = buffer.data(smk0 + 516);
    const auto *smk0_518 = buffer.data(smk0 + 518);
    const auto *smk0_519 = buffer.data(smk0 + 519);
    const auto *smk0_521 = buffer.data(smk0 + 521);
    const auto *smk0_522 = buffer.data(smk0 + 522);
    const auto *smk0_524 = buffer.data(smk0 + 524);
    const auto *smk0_539 = buffer.data(smk0 + 539);

    const auto *smi_363 = buffer.data(smi + 363);
    const auto *smi_364 = buffer.data(smi + 364);
    const auto *smi_367 = buffer.data(smi + 367);
    const auto *smi_370 = buffer.data(smi + 370);
    const auto *smi_374 = buffer.data(smi + 374);
    const auto *smi_385 = buffer.data(smi + 385);
    const auto *smi_388 = buffer.data(smi + 388);
    const auto *smi_389 = buffer.data(smi + 389);
    const auto *smi_390 = buffer.data(smi + 390);
    const auto *smi_391 = buffer.data(smi + 391);
    const auto *smi_392 = buffer.data(smi + 392);
    const auto *smi_393 = buffer.data(smi + 393);
    const auto *smi_394 = buffer.data(smi + 394);
    const auto *smi_395 = buffer.data(smi + 395);
    const auto *smi_397 = buffer.data(smi + 397);
    const auto *smi_398 = buffer.data(smi + 398);
    const auto *smi_400 = buffer.data(smi + 400);
    const auto *smi_401 = buffer.data(smi + 401);
    const auto *smi_402 = buffer.data(smi + 402);
    const auto *smi_404 = buffer.data(smi + 404);
    const auto *smi_405 = buffer.data(smi + 405);
    const auto *smi_406 = buffer.data(smi + 406);
    const auto *smi_413 = buffer.data(smi + 413);
    const auto *smi_415 = buffer.data(smi + 415);
    const auto *smi_416 = buffer.data(smi + 416);
    const auto *smi_417 = buffer.data(smi + 417);
    const auto *smi_418 = buffer.data(smi + 418);
    const auto *smi_419 = buffer.data(smi + 419);
    const auto *smi_420 = buffer.data(smi + 420);
    const auto *smi_422 = buffer.data(smi + 422);
    const auto *smi_425 = buffer.data(smi + 425);
    const auto *smi_429 = buffer.data(smi + 429);
    const auto *smi_434 = buffer.data(smi + 434);
    const auto *smi_441 = buffer.data(smi + 441);
    const auto *smi_443 = buffer.data(smi + 443);
    const auto *smi_444 = buffer.data(smi + 444);
    const auto *smi_445 = buffer.data(smi + 445);
    const auto *smi_446 = buffer.data(smi + 446);
    const auto *smi_553 = buffer.data(smi + 553);
    const auto *smi_554 = buffer.data(smi + 554);
    const auto *smi_555 = buffer.data(smi + 555);
    const auto *smi_556 = buffer.data(smi + 556);
    const auto *smi_557 = buffer.data(smi + 557);
    const auto *smi_558 = buffer.data(smi + 558);
    const auto *smi_559 = buffer.data(smi + 559);
    const auto *smi_560 = buffer.data(smi + 560);
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

    const auto *smk1_504 = buffer.data(smk1 + 504);
    const auto *smk1_507 = buffer.data(smk1 + 507);
    const auto *smk1_509 = buffer.data(smk1 + 509);
    const auto *smk1_510 = buffer.data(smk1 + 510);
    const auto *smk1_513 = buffer.data(smk1 + 513);
    const auto *smk1_514 = buffer.data(smk1 + 514);
    const auto *smk1_516 = buffer.data(smk1 + 516);
    const auto *smk1_518 = buffer.data(smk1 + 518);
    const auto *smk1_519 = buffer.data(smk1 + 519);
    const auto *smk1_521 = buffer.data(smk1 + 521);
    const auto *smk1_522 = buffer.data(smk1 + 522);
    const auto *smk1_524 = buffer.data(smk1 + 524);
    const auto *smk1_539 = buffer.data(smk1 + 539);

    const auto *snh0_396 = buffer.data(snh0 + 396);
    const auto *snh0_397 = buffer.data(snh0 + 397);
    const auto *snh0_398 = buffer.data(snh0 + 398);
    const auto *snh0_414 = buffer.data(snh0 + 414);
    const auto *snh0_416 = buffer.data(snh0 + 416);
    const auto *snh0_417 = buffer.data(snh0 + 417);
    const auto *snh0_418 = buffer.data(snh0 + 418);
    const auto *snh0_419 = buffer.data(snh0 + 419);
    const auto *snh0_420 = buffer.data(snh0 + 420);
    const auto *snh0_423 = buffer.data(snh0 + 423);
    const auto *snh0_425 = buffer.data(snh0 + 425);
    const auto *snh0_426 = buffer.data(snh0 + 426);
    const auto *snh0_429 = buffer.data(snh0 + 429);
    const auto *snh0_430 = buffer.data(snh0 + 430);
    const auto *snh0_432 = buffer.data(snh0 + 432);
    const auto *snh0_434 = buffer.data(snh0 + 434);
    const auto *snh0_435 = buffer.data(snh0 + 435);
    const auto *snh0_437 = buffer.data(snh0 + 437);
    const auto *snh0_438 = buffer.data(snh0 + 438);
    const auto *snh0_439 = buffer.data(snh0 + 439);
    const auto *snh0_440 = buffer.data(snh0 + 440);
    const auto *snh0_441 = buffer.data(snh0 + 441);
    const auto *snh0_444 = buffer.data(snh0 + 444);
    const auto *snh0_446 = buffer.data(snh0 + 446);
    const auto *snh0_447 = buffer.data(snh0 + 447);
    const auto *snh0_450 = buffer.data(snh0 + 450);
    const auto *snh0_451 = buffer.data(snh0 + 451);
    const auto *snh0_453 = buffer.data(snh0 + 453);
    const auto *snh0_455 = buffer.data(snh0 + 455);
    const auto *snh0_456 = buffer.data(snh0 + 456);
    const auto *snh0_458 = buffer.data(snh0 + 458);
    const auto *snh0_459 = buffer.data(snh0 + 459);
    const auto *snh0_460 = buffer.data(snh0 + 460);
    const auto *snh0_461 = buffer.data(snh0 + 461);

    const auto *snh1_396 = buffer.data(snh1 + 396);
    const auto *snh1_397 = buffer.data(snh1 + 397);
    const auto *snh1_398 = buffer.data(snh1 + 398);
    const auto *snh1_414 = buffer.data(snh1 + 414);
    const auto *snh1_416 = buffer.data(snh1 + 416);
    const auto *snh1_417 = buffer.data(snh1 + 417);
    const auto *snh1_418 = buffer.data(snh1 + 418);
    const auto *snh1_419 = buffer.data(snh1 + 419);
    const auto *snh1_420 = buffer.data(snh1 + 420);
    const auto *snh1_423 = buffer.data(snh1 + 423);
    const auto *snh1_425 = buffer.data(snh1 + 425);
    const auto *snh1_426 = buffer.data(snh1 + 426);
    const auto *snh1_429 = buffer.data(snh1 + 429);
    const auto *snh1_430 = buffer.data(snh1 + 430);
    const auto *snh1_432 = buffer.data(snh1 + 432);
    const auto *snh1_434 = buffer.data(snh1 + 434);
    const auto *snh1_435 = buffer.data(snh1 + 435);
    const auto *snh1_437 = buffer.data(snh1 + 437);
    const auto *snh1_438 = buffer.data(snh1 + 438);
    const auto *snh1_439 = buffer.data(snh1 + 439);
    const auto *snh1_440 = buffer.data(snh1 + 440);
    const auto *snh1_441 = buffer.data(snh1 + 441);
    const auto *snh1_444 = buffer.data(snh1 + 444);
    const auto *snh1_446 = buffer.data(snh1 + 446);
    const auto *snh1_447 = buffer.data(snh1 + 447);
    const auto *snh1_450 = buffer.data(snh1 + 450);
    const auto *snh1_451 = buffer.data(snh1 + 451);
    const auto *snh1_453 = buffer.data(snh1 + 453);
    const auto *snh1_455 = buffer.data(snh1 + 455);
    const auto *snh1_456 = buffer.data(snh1 + 456);
    const auto *snh1_458 = buffer.data(snh1 + 458);
    const auto *snh1_459 = buffer.data(snh1 + 459);
    const auto *snh1_460 = buffer.data(snh1 + 460);
    const auto *snh1_461 = buffer.data(snh1 + 461);

    const auto *sni_528 = buffer.data(sni + 528);
    const auto *sni_529 = buffer.data(sni + 529);
    const auto *sni_530 = buffer.data(sni + 530);
    const auto *sni_531 = buffer.data(sni + 531);
    const auto *sni_532 = buffer.data(sni + 532);
    const auto *sni_534 = buffer.data(sni + 534);
    const auto *sni_535 = buffer.data(sni + 535);
    const auto *sni_537 = buffer.data(sni + 537);
    const auto *sni_538 = buffer.data(sni + 538);
    const auto *sni_541 = buffer.data(sni + 541);
    const auto *sni_542 = buffer.data(sni + 542);
    const auto *sni_546 = buffer.data(sni + 546);
    const auto *sni_553 = buffer.data(sni + 553);
    const auto *sni_554 = buffer.data(sni + 554);
    const auto *sni_555 = buffer.data(sni + 555);
    const auto *sni_556 = buffer.data(sni + 556);
    const auto *sni_557 = buffer.data(sni + 557);
    const auto *sni_558 = buffer.data(sni + 558);
    const auto *sni_559 = buffer.data(sni + 559);
    const auto *sni_560 = buffer.data(sni + 560);
    const auto *sni_562 = buffer.data(sni + 562);
    const auto *sni_563 = buffer.data(sni + 563);
    const auto *sni_565 = buffer.data(sni + 565);
    const auto *sni_566 = buffer.data(sni + 566);
    const auto *sni_569 = buffer.data(sni + 569);
    const auto *sni_570 = buffer.data(sni + 570);
    const auto *sni_572 = buffer.data(sni + 572);
    const auto *sni_574 = buffer.data(sni + 574);
    const auto *sni_575 = buffer.data(sni + 575);
    const auto *sni_577 = buffer.data(sni + 577);
    const auto *sni_578 = buffer.data(sni + 578);
    const auto *sni_580 = buffer.data(sni + 580);
    const auto *sni_581 = buffer.data(sni + 581);
    const auto *sni_582 = buffer.data(sni + 582);
    const auto *sni_583 = buffer.data(sni + 583);
    const auto *sni_584 = buffer.data(sni + 584);
    const auto *sni_585 = buffer.data(sni + 585);
    const auto *sni_586 = buffer.data(sni + 586);
    const auto *sni_587 = buffer.data(sni + 587);
    const auto *sni_588 = buffer.data(sni + 588);
    const auto *sni_590 = buffer.data(sni + 590);
    const auto *sni_591 = buffer.data(sni + 591);
    const auto *sni_593 = buffer.data(sni + 593);
    const auto *sni_594 = buffer.data(sni + 594);
    const auto *sni_597 = buffer.data(sni + 597);
    const auto *sni_598 = buffer.data(sni + 598);
    const auto *sni_600 = buffer.data(sni + 600);
    const auto *sni_602 = buffer.data(sni + 602);
    const auto *sni_603 = buffer.data(sni + 603);
    const auto *sni_605 = buffer.data(sni + 605);
    const auto *sni_606 = buffer.data(sni + 606);
    const auto *sni_608 = buffer.data(sni + 608);
    const auto *sni_609 = buffer.data(sni + 609);
    const auto *sni_610 = buffer.data(sni + 610);
    const auto *sni_611 = buffer.data(sni + 611);
    const auto *sni_612 = buffer.data(sni + 612);
    const auto *sni_613 = buffer.data(sni + 613);
    const auto *sni_614 = buffer.data(sni + 614);
    const auto *sni_615 = buffer.data(sni + 615);

#pragma omp simd aligned(t_679, t_680, t_681, pc_y, smi_388, smi_389, smi_390, snh0_396, \
                         snh0_397, snh0_398, snh1_396, snh1_397, snh1_398, sni_528, sni_529, \
                         sni_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_14 * smi_388[k]
                   + f_6 * snh0_396[k]
                   - f_7 * snh1_396[k]
                   + f_3 * pc_y[k] * sni_528[k];

        t_680[k] = f_14 * smi_389[k]
                   + f_8 * snh0_397[k]
                   - f_9 * snh1_397[k]
                   + f_3 * pc_y[k] * sni_529[k];

        t_681[k] = f_14 * smi_390[k]
                   + f_10 * snh0_398[k]
                   - f_11 * snh1_398[k]
                   + f_3 * pc_y[k] * sni_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, pb_y, pc_y, pc_z, smk0_504, smi_363, \
                         smi_391, smi_392, smk1_504, snh0_398, snh1_398, sni_531, \
                         sni_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_14 * smi_391[k]
                   + f_3 * pc_y[k] * sni_531[k];

        t_683[k] = f_15 * smi_363[k]
                   + f_1 * snh0_398[k]
                   - f_2 * snh1_398[k]
                   + f_3 * pc_z[k] * sni_531[k];

        t_684[k] = pb_y[k] * smk0_504[k]
                   - f_12 * pc_y[k] * smk1_504[k];

        t_685[k] = f_13 * smi_392[k]
                   + f_3 * pc_y[k] * sni_532[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pb_y, pc_y, pc_z, smk0_507, smk0_509, \
                         smi_364, smi_393, smi_394, smk1_507, smk1_509, sni_532, \
                         sni_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_16 * smi_364[k]
                   + f_3 * pc_z[k] * sni_532[k];

        t_687[k] = pb_y[k] * smk0_507[k]
                   + f_14 * smi_393[k]
                   - f_12 * pc_y[k] * smk1_507[k];

        t_688[k] = f_13 * smi_394[k]
                   + f_3 * pc_y[k] * sni_534[k];

        t_689[k] = pb_y[k] * smk0_509[k]
                   - f_12 * pc_y[k] * smk1_509[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pb_y, pc_y, pc_z, smk0_510, smk0_513, \
                         smi_367, smi_395, smi_397, smk1_510, smk1_513, sni_535, \
                         sni_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = pb_y[k] * smk0_510[k]
                   + f_15 * smi_395[k]
                   - f_12 * pc_y[k] * smk1_510[k];

        t_691[k] = f_16 * smi_367[k]
                   + f_3 * pc_z[k] * sni_535[k];

        t_692[k] = f_13 * smi_397[k]
                   + f_3 * pc_y[k] * sni_537[k];

        t_693[k] = pb_y[k] * smk0_513[k]
                   - f_12 * pc_y[k] * smk1_513[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pb_y, pc_y, pc_z, smk0_514, smk0_516, smi_370, \
                         smi_398, smi_400, smk1_514, smk1_516, \
                         sni_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pb_y[k] * smk0_514[k]
                   + f_16 * smi_398[k]
                   - f_12 * pc_y[k] * smk1_514[k];

        t_695[k] = f_16 * smi_370[k]
                   + f_3 * pc_z[k] * sni_538[k];

        t_696[k] = pb_y[k] * smk0_516[k]
                   + f_14 * smi_400[k]
                   - f_12 * pc_y[k] * smk1_516[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, t_700, pb_y, pc_y, pc_z, smk0_518, smk0_519, \
                         smi_374, smi_401, smi_402, smk1_518, smk1_519, sni_541, \
                         sni_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * smi_401[k]
                   + f_3 * pc_y[k] * sni_541[k];

        t_698[k] = pb_y[k] * smk0_518[k]
                   - f_12 * pc_y[k] * smk1_518[k];

        t_699[k] = pb_y[k] * smk0_519[k]
                   + f_17 * smi_402[k]
                   - f_12 * pc_y[k] * smk1_519[k];

        t_700[k] = f_16 * smi_374[k]
                   + f_3 * pc_z[k] * sni_542[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pb_y, pc_y, smk0_521, smk0_522, smk0_524, \
                         smi_404, smi_405, smi_406, smk1_521, smk1_522, smk1_524, \
                         sni_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = pb_y[k] * smk0_521[k]
                   + f_15 * smi_404[k]
                   - f_12 * pc_y[k] * smk1_521[k];

        t_702[k] = pb_y[k] * smk0_522[k]
                   + f_14 * smi_405[k]
                   - f_12 * pc_y[k] * smk1_522[k];

        t_703[k] = f_13 * smi_406[k]
                   + f_3 * pc_y[k] * sni_546[k];

        t_704[k] = pb_y[k] * smk0_524[k]
                   - f_12 * pc_y[k] * smk1_524[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, smi_553, smi_554, smi_555, \
                         smi_556, smi_557, sni_553, sni_554, sni_555, sni_556, \
                         sni_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_17 * smi_553[k]
                   + f_3 * pc_x[k] * sni_553[k];

        t_706[k] = f_17 * smi_554[k]
                   + f_3 * pc_x[k] * sni_554[k];

        t_707[k] = f_17 * smi_555[k]
                   + f_3 * pc_x[k] * sni_555[k];

        t_708[k] = f_17 * smi_556[k]
                   + f_3 * pc_x[k] * sni_556[k];

        t_709[k] = f_17 * smi_557[k]
                   + f_3 * pc_x[k] * sni_557[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, pc_y, pc_z, smi_385, smi_413, \
                         smi_558, smi_559, snh0_414, snh1_414, sni_553, sni_558, \
                         sni_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_17 * smi_558[k]
                   + f_3 * pc_x[k] * sni_558[k];

        t_711[k] = f_17 * smi_559[k]
                   + f_3 * pc_x[k] * sni_559[k];

        t_712[k] = f_13 * smi_413[k]
                   + f_1 * snh0_414[k]
                   - f_2 * snh1_414[k]
                   + f_3 * pc_y[k] * sni_553[k];

        t_713[k] = f_16 * smi_385[k]
                   + f_3 * pc_z[k] * sni_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, smi_415, smi_416, smi_417, snh0_416, \
                         snh0_417, snh0_418, snh1_416, snh1_417, snh1_418, sni_555, sni_556, \
                         sni_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_13 * smi_415[k]
                   + f_4 * snh0_416[k]
                   - f_5 * snh1_416[k]
                   + f_3 * pc_y[k] * sni_555[k];

        t_715[k] = f_13 * smi_416[k]
                   + f_6 * snh0_417[k]
                   - f_7 * snh1_417[k]
                   + f_3 * pc_y[k] * sni_556[k];

        t_716[k] = f_13 * smi_417[k]
                   + f_8 * snh0_418[k]
                   - f_9 * snh1_418[k]
                   + f_3 * pc_y[k] * sni_557[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pb_y, pc_y, smk0_539, smi_418, smi_419, \
                         smk1_539, snh0_419, snh1_419, sni_558, \
                         sni_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_13 * smi_418[k]
                   + f_10 * snh0_419[k]
                   - f_11 * snh1_419[k]
                   + f_3 * pc_y[k] * sni_558[k];

        t_718[k] = f_13 * smi_419[k]
                   + f_3 * pc_y[k] * sni_559[k];

        t_719[k] = pb_y[k] * smk0_539[k]
                   - f_12 * pc_y[k] * smk1_539[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pc_x, pc_y, pc_z, smi_392, smi_560, \
                         smi_563, snh0_420, snh0_423, snh1_420, snh1_423, sni_560, \
                         sni_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_17 * smi_560[k]
                   + f_1 * snh0_420[k]
                   - f_2 * snh1_420[k]
                   + f_3 * pc_x[k] * sni_560[k];

        t_721[k] = f_3 * pc_y[k] * sni_560[k];

        t_722[k] = f_17 * smi_392[k]
                   + f_3 * pc_z[k] * sni_560[k];

        t_723[k] = f_17 * smi_563[k]
                   + f_4 * snh0_423[k]
                   - f_5 * snh1_423[k]
                   + f_3 * pc_x[k] * sni_563[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, pc_x, pc_y, smi_565, smi_566, snh0_425, \
                         snh0_426, snh1_425, snh1_426, sni_562, sni_565, \
                         sni_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_3 * pc_y[k] * sni_562[k];

        t_725[k] = f_17 * smi_565[k]
                   + f_4 * snh0_425[k]
                   - f_5 * snh1_425[k]
                   + f_3 * pc_x[k] * sni_565[k];

        t_726[k] = f_17 * smi_566[k]
                   + f_6 * snh0_426[k]
                   - f_7 * snh1_426[k]
                   + f_3 * pc_x[k] * sni_566[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_x, pc_y, pc_z, smi_395, smi_569, snh0_429, \
                         snh1_429, sni_563, sni_565, sni_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_17 * smi_395[k]
                   + f_3 * pc_z[k] * sni_563[k];

        t_728[k] = f_3 * pc_y[k] * sni_565[k];

        t_729[k] = f_17 * smi_569[k]
                   + f_6 * snh0_429[k]
                   - f_7 * snh1_429[k]
                   + f_3 * pc_x[k] * sni_569[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pc_x, pc_z, smi_398, smi_570, smi_572, snh0_430, \
                         snh0_432, snh1_430, snh1_432, sni_566, sni_570, \
                         sni_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_17 * smi_570[k]
                   + f_8 * snh0_430[k]
                   - f_9 * snh1_430[k]
                   + f_3 * pc_x[k] * sni_570[k];

        t_731[k] = f_17 * smi_398[k]
                   + f_3 * pc_z[k] * sni_566[k];

        t_732[k] = f_17 * smi_572[k]
                   + f_8 * snh0_432[k]
                   - f_9 * snh1_432[k]
                   + f_3 * pc_x[k] * sni_572[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pc_x, pc_y, smi_574, smi_575, snh0_434, \
                         snh0_435, snh1_434, snh1_435, sni_569, sni_574, \
                         sni_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_3 * pc_y[k] * sni_569[k];

        t_734[k] = f_17 * smi_574[k]
                   + f_8 * snh0_434[k]
                   - f_9 * snh1_434[k]
                   + f_3 * pc_x[k] * sni_574[k];

        t_735[k] = f_17 * smi_575[k]
                   + f_10 * snh0_435[k]
                   - f_11 * snh1_435[k]
                   + f_3 * pc_x[k] * sni_575[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pc_x, pc_z, smi_402, smi_577, smi_578, snh0_437, \
                         snh0_438, snh1_437, snh1_438, sni_570, sni_577, \
                         sni_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_17 * smi_402[k]
                   + f_3 * pc_z[k] * sni_570[k];

        t_737[k] = f_17 * smi_577[k]
                   + f_10 * snh0_437[k]
                   - f_11 * snh1_437[k]
                   + f_3 * pc_x[k] * sni_577[k];

        t_738[k] = f_17 * smi_578[k]
                   + f_10 * snh0_438[k]
                   - f_11 * snh1_438[k]
                   + f_3 * pc_x[k] * sni_578[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, pc_x, pc_y, smi_580, smi_581, smi_582, \
                         snh0_440, snh1_440, sni_574, sni_580, sni_581, \
                         sni_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * sni_574[k];

        t_740[k] = f_17 * smi_580[k]
                   + f_10 * snh0_440[k]
                   - f_11 * snh1_440[k]
                   + f_3 * pc_x[k] * sni_580[k];

        t_741[k] = f_17 * smi_581[k]
                   + f_3 * pc_x[k] * sni_581[k];

        t_742[k] = f_17 * smi_582[k]
                   + f_3 * pc_x[k] * sni_582[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, pc_x, smi_583, smi_584, smi_585, \
                         smi_586, smi_587, sni_583, sni_584, sni_585, sni_586, \
                         sni_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_17 * smi_583[k]
                   + f_3 * pc_x[k] * sni_583[k];

        t_744[k] = f_17 * smi_584[k]
                   + f_3 * pc_x[k] * sni_584[k];

        t_745[k] = f_17 * smi_585[k]
                   + f_3 * pc_x[k] * sni_585[k];

        t_746[k] = f_17 * smi_586[k]
                   + f_3 * pc_x[k] * sni_586[k];

        t_747[k] = f_17 * smi_587[k]
                   + f_3 * pc_x[k] * sni_587[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_y, pc_z, smi_413, snh0_435, snh0_437, \
                         snh0_438, snh1_435, snh1_437, snh1_438, sni_581, sni_583, \
                         sni_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_1 * snh0_435[k]
                   - f_2 * snh1_435[k]
                   + f_3 * pc_y[k] * sni_581[k];

        t_749[k] = f_17 * smi_413[k]
                   + f_3 * pc_z[k] * sni_581[k];

        t_750[k] = f_4 * snh0_437[k]
                   - f_5 * snh1_437[k]
                   + f_3 * pc_y[k] * sni_583[k];

        t_751[k] = f_6 * snh0_438[k]
                   - f_7 * snh1_438[k]
                   + f_3 * pc_y[k] * sni_584[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pc_y, pc_z, smi_419, snh0_439, snh0_440, \
                         snh1_439, snh1_440, sni_585, sni_586, \
                         sni_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_8 * snh0_439[k]
                   - f_9 * snh1_439[k]
                   + f_3 * pc_y[k] * sni_585[k];

        t_753[k] = f_10 * snh0_440[k]
                   - f_11 * snh1_440[k]
                   + f_3 * pc_y[k] * sni_586[k];

        t_754[k] = f_3 * pc_y[k] * sni_587[k];

        t_755[k] = f_17 * smi_419[k]
                   + f_1 * snh0_440[k]
                   - f_2 * snh1_440[k]
                   + f_3 * pc_z[k] * sni_587[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, pc_y, pc_z, smi_420, smi_588, \
                         smi_591, snh0_441, snh0_444, snh1_441, snh1_444, sni_588, \
                         sni_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_16 * smi_588[k]
                   + f_1 * snh0_441[k]
                   - f_2 * snh1_441[k]
                   + f_3 * pc_x[k] * sni_588[k];

        t_757[k] = f_21 * smi_420[k]
                   + f_3 * pc_y[k] * sni_588[k];

        t_758[k] = f_3 * pc_z[k] * sni_588[k];

        t_759[k] = f_16 * smi_591[k]
                   + f_4 * snh0_444[k]
                   - f_5 * snh1_444[k]
                   + f_3 * pc_x[k] * sni_591[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, pc_x, pc_y, smi_422, smi_593, smi_594, snh0_446, \
                         snh0_447, snh1_446, snh1_447, sni_590, sni_593, \
                         sni_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_21 * smi_422[k]
                   + f_3 * pc_y[k] * sni_590[k];

        t_761[k] = f_16 * smi_593[k]
                   + f_4 * snh0_446[k]
                   - f_5 * snh1_446[k]
                   + f_3 * pc_x[k] * sni_593[k];

        t_762[k] = f_16 * smi_594[k]
                   + f_6 * snh0_447[k]
                   - f_7 * snh1_447[k]
                   + f_3 * pc_x[k] * sni_594[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, pc_x, pc_y, pc_z, smi_425, smi_597, snh0_450, \
                         snh1_450, sni_591, sni_593, sni_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_3 * pc_z[k] * sni_591[k];

        t_764[k] = f_21 * smi_425[k]
                   + f_3 * pc_y[k] * sni_593[k];

        t_765[k] = f_16 * smi_597[k]
                   + f_6 * snh0_450[k]
                   - f_7 * snh1_450[k]
                   + f_3 * pc_x[k] * sni_597[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, pc_x, pc_z, smi_598, smi_600, snh0_451, \
                         snh0_453, snh1_451, snh1_453, sni_594, sni_598, \
                         sni_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_16 * smi_598[k]
                   + f_8 * snh0_451[k]
                   - f_9 * snh1_451[k]
                   + f_3 * pc_x[k] * sni_598[k];

        t_767[k] = f_3 * pc_z[k] * sni_594[k];

        t_768[k] = f_16 * smi_600[k]
                   + f_8 * snh0_453[k]
                   - f_9 * snh1_453[k]
                   + f_3 * pc_x[k] * sni_600[k];
    }

#pragma omp simd aligned(t_769, t_770, t_771, pc_x, pc_y, smi_429, smi_602, smi_603, snh0_455, \
                         snh0_456, snh1_455, snh1_456, sni_597, sni_602, \
                         sni_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_769[k] = f_21 * smi_429[k]
                   + f_3 * pc_y[k] * sni_597[k];

        t_770[k] = f_16 * smi_602[k]
                   + f_8 * snh0_455[k]
                   - f_9 * snh1_455[k]
                   + f_3 * pc_x[k] * sni_602[k];

        t_771[k] = f_16 * smi_603[k]
                   + f_10 * snh0_456[k]
                   - f_11 * snh1_456[k]
                   + f_3 * pc_x[k] * sni_603[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, pc_x, pc_z, smi_605, smi_606, snh0_458, \
                         snh0_459, snh1_458, snh1_459, sni_598, sni_605, \
                         sni_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * sni_598[k];

        t_773[k] = f_16 * smi_605[k]
                   + f_10 * snh0_458[k]
                   - f_11 * snh1_458[k]
                   + f_3 * pc_x[k] * sni_605[k];

        t_774[k] = f_16 * smi_606[k]
                   + f_10 * snh0_459[k]
                   - f_11 * snh1_459[k]
                   + f_3 * pc_x[k] * sni_606[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, pc_x, pc_y, smi_434, smi_608, smi_609, \
                         smi_610, snh0_461, snh1_461, sni_602, sni_608, sni_609, \
                         sni_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_21 * smi_434[k]
                   + f_3 * pc_y[k] * sni_602[k];

        t_776[k] = f_16 * smi_608[k]
                   + f_10 * snh0_461[k]
                   - f_11 * snh1_461[k]
                   + f_3 * pc_x[k] * sni_608[k];

        t_777[k] = f_16 * smi_609[k]
                   + f_3 * pc_x[k] * sni_609[k];

        t_778[k] = f_16 * smi_610[k]
                   + f_3 * pc_x[k] * sni_610[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, t_782, t_783, pc_x, smi_611, smi_612, smi_613, \
                         smi_614, smi_615, sni_611, sni_612, sni_613, sni_614, \
                         sni_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_16 * smi_611[k]
                   + f_3 * pc_x[k] * sni_611[k];

        t_780[k] = f_16 * smi_612[k]
                   + f_3 * pc_x[k] * sni_612[k];

        t_781[k] = f_16 * smi_613[k]
                   + f_3 * pc_x[k] * sni_613[k];

        t_782[k] = f_16 * smi_614[k]
                   + f_3 * pc_x[k] * sni_614[k];

        t_783[k] = f_16 * smi_615[k]
                   + f_3 * pc_x[k] * sni_615[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, pc_y, pc_z, smi_441, smi_443, snh0_456, \
                         snh0_458, snh1_456, snh1_458, sni_609, \
                         sni_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_21 * smi_441[k]
                   + f_1 * snh0_456[k]
                   - f_2 * snh1_456[k]
                   + f_3 * pc_y[k] * sni_609[k];

        t_785[k] = f_3 * pc_z[k] * sni_609[k];

        t_786[k] = f_21 * smi_443[k]
                   + f_4 * snh0_458[k]
                   - f_5 * snh1_458[k]
                   + f_3 * pc_y[k] * sni_611[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pc_y, smi_444, smi_445, smi_446, snh0_459, \
                         snh0_460, snh0_461, snh1_459, snh1_460, snh1_461, sni_612, sni_613, \
                         sni_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_21 * smi_444[k]
                   + f_6 * snh0_459[k]
                   - f_7 * snh1_459[k]
                   + f_3 * pc_y[k] * sni_612[k];

        t_788[k] = f_21 * smi_445[k]
                   + f_8 * snh0_460[k]
                   - f_9 * snh1_460[k]
                   + f_3 * pc_y[k] * sni_613[k];

        t_789[k] = f_21 * smi_446[k]
                   + f_10 * snh0_461[k]
                   - f_11 * snh1_461[k]
                   + f_3 * pc_y[k] * sni_614[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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
    const auto f_21 = 3.0 / q;

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

    const auto *smk0_540 = buffer.data(smk0 + 540);
    const auto *smk0_543 = buffer.data(smk0 + 543);
    const auto *smk0_546 = buffer.data(smk0 + 546);
    const auto *smk0_550 = buffer.data(smk0 + 550);
    const auto *smk0_552 = buffer.data(smk0 + 552);
    const auto *smk0_555 = buffer.data(smk0 + 555);
    const auto *smk0_557 = buffer.data(smk0 + 557);
    const auto *smk0_558 = buffer.data(smk0 + 558);
    const auto *smk0_568 = buffer.data(smk0 + 568);

    const auto *smi_420 = buffer.data(smi + 420);
    const auto *smi_423 = buffer.data(smi + 423);
    const auto *smi_426 = buffer.data(smi + 426);
    const auto *smi_427 = buffer.data(smi + 427);
    const auto *smi_430 = buffer.data(smi + 430);
    const auto *smi_431 = buffer.data(smi + 431);
    const auto *smi_432 = buffer.data(smi + 432);
    const auto *smi_441 = buffer.data(smi + 441);
    const auto *smi_447 = buffer.data(smi + 447);
    const auto *smi_448 = buffer.data(smi + 448);
    const auto *smi_450 = buffer.data(smi + 450);
    const auto *smi_451 = buffer.data(smi + 451);
    const auto *smi_453 = buffer.data(smi + 453);
    const auto *smi_454 = buffer.data(smi + 454);
    const auto *smi_457 = buffer.data(smi + 457);
    const auto *smi_458 = buffer.data(smi + 458);
    const auto *smi_462 = buffer.data(smi + 462);
    const auto *smi_469 = buffer.data(smi + 469);
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
    const auto *smi_490 = buffer.data(smi + 490);
    const auto *smi_497 = buffer.data(smi + 497);
    const auto *smi_499 = buffer.data(smi + 499);
    const auto *smi_500 = buffer.data(smi + 500);
    const auto *smi_501 = buffer.data(smi + 501);
    const auto *smi_502 = buffer.data(smi + 502);
    const auto *smi_503 = buffer.data(smi + 503);
    const auto *smi_504 = buffer.data(smi + 504);
    const auto *smi_506 = buffer.data(smi + 506);
    const auto *smi_509 = buffer.data(smi + 509);
    const auto *smi_513 = buffer.data(smi + 513);
    const auto *smi_518 = buffer.data(smi + 518);
    const auto *smi_525 = buffer.data(smi + 525);
    const auto *smi_527 = buffer.data(smi + 527);
    const auto *smi_528 = buffer.data(smi + 528);
    const auto *smi_529 = buffer.data(smi + 529);
    const auto *smi_530 = buffer.data(smi + 530);
    const auto *smi_621 = buffer.data(smi + 621);
    const auto *smi_625 = buffer.data(smi + 625);
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

    const auto *smk1_540 = buffer.data(smk1 + 540);
    const auto *smk1_543 = buffer.data(smk1 + 543);
    const auto *smk1_546 = buffer.data(smk1 + 546);
    const auto *smk1_550 = buffer.data(smk1 + 550);
    const auto *smk1_552 = buffer.data(smk1 + 552);
    const auto *smk1_555 = buffer.data(smk1 + 555);
    const auto *smk1_557 = buffer.data(smk1 + 557);
    const auto *smk1_558 = buffer.data(smk1 + 558);
    const auto *smk1_568 = buffer.data(smk1 + 568);

    const auto *snh0_461 = buffer.data(snh0 + 461);
    const auto *snh0_467 = buffer.data(snh0 + 467);
    const auto *snh0_471 = buffer.data(snh0 + 471);
    const auto *snh0_476 = buffer.data(snh0 + 476);
    const auto *snh0_479 = buffer.data(snh0 + 479);
    const auto *snh0_480 = buffer.data(snh0 + 480);
    const auto *snh0_481 = buffer.data(snh0 + 481);
    const auto *snh0_482 = buffer.data(snh0 + 482);
    const auto *snh0_483 = buffer.data(snh0 + 483);
    const auto *snh0_486 = buffer.data(snh0 + 486);
    const auto *snh0_488 = buffer.data(snh0 + 488);
    const auto *snh0_489 = buffer.data(snh0 + 489);
    const auto *snh0_492 = buffer.data(snh0 + 492);
    const auto *snh0_493 = buffer.data(snh0 + 493);
    const auto *snh0_495 = buffer.data(snh0 + 495);
    const auto *snh0_497 = buffer.data(snh0 + 497);
    const auto *snh0_498 = buffer.data(snh0 + 498);
    const auto *snh0_500 = buffer.data(snh0 + 500);
    const auto *snh0_501 = buffer.data(snh0 + 501);
    const auto *snh0_502 = buffer.data(snh0 + 502);
    const auto *snh0_503 = buffer.data(snh0 + 503);
    const auto *snh0_504 = buffer.data(snh0 + 504);
    const auto *snh0_507 = buffer.data(snh0 + 507);
    const auto *snh0_509 = buffer.data(snh0 + 509);
    const auto *snh0_510 = buffer.data(snh0 + 510);
    const auto *snh0_513 = buffer.data(snh0 + 513);
    const auto *snh0_514 = buffer.data(snh0 + 514);
    const auto *snh0_516 = buffer.data(snh0 + 516);
    const auto *snh0_518 = buffer.data(snh0 + 518);
    const auto *snh0_519 = buffer.data(snh0 + 519);
    const auto *snh0_521 = buffer.data(snh0 + 521);
    const auto *snh0_522 = buffer.data(snh0 + 522);
    const auto *snh0_523 = buffer.data(snh0 + 523);
    const auto *snh0_524 = buffer.data(snh0 + 524);

    const auto *snh1_461 = buffer.data(snh1 + 461);
    const auto *snh1_467 = buffer.data(snh1 + 467);
    const auto *snh1_471 = buffer.data(snh1 + 471);
    const auto *snh1_476 = buffer.data(snh1 + 476);
    const auto *snh1_479 = buffer.data(snh1 + 479);
    const auto *snh1_480 = buffer.data(snh1 + 480);
    const auto *snh1_481 = buffer.data(snh1 + 481);
    const auto *snh1_482 = buffer.data(snh1 + 482);
    const auto *snh1_483 = buffer.data(snh1 + 483);
    const auto *snh1_486 = buffer.data(snh1 + 486);
    const auto *snh1_488 = buffer.data(snh1 + 488);
    const auto *snh1_489 = buffer.data(snh1 + 489);
    const auto *snh1_492 = buffer.data(snh1 + 492);
    const auto *snh1_493 = buffer.data(snh1 + 493);
    const auto *snh1_495 = buffer.data(snh1 + 495);
    const auto *snh1_497 = buffer.data(snh1 + 497);
    const auto *snh1_498 = buffer.data(snh1 + 498);
    const auto *snh1_500 = buffer.data(snh1 + 500);
    const auto *snh1_501 = buffer.data(snh1 + 501);
    const auto *snh1_502 = buffer.data(snh1 + 502);
    const auto *snh1_503 = buffer.data(snh1 + 503);
    const auto *snh1_504 = buffer.data(snh1 + 504);
    const auto *snh1_507 = buffer.data(snh1 + 507);
    const auto *snh1_509 = buffer.data(snh1 + 509);
    const auto *snh1_510 = buffer.data(snh1 + 510);
    const auto *snh1_513 = buffer.data(snh1 + 513);
    const auto *snh1_514 = buffer.data(snh1 + 514);
    const auto *snh1_516 = buffer.data(snh1 + 516);
    const auto *snh1_518 = buffer.data(snh1 + 518);
    const auto *snh1_519 = buffer.data(snh1 + 519);
    const auto *snh1_521 = buffer.data(snh1 + 521);
    const auto *snh1_522 = buffer.data(snh1 + 522);
    const auto *snh1_523 = buffer.data(snh1 + 523);
    const auto *snh1_524 = buffer.data(snh1 + 524);

    const auto *sni_615 = buffer.data(sni + 615);
    const auto *sni_616 = buffer.data(sni + 616);
    const auto *sni_618 = buffer.data(sni + 618);
    const auto *sni_619 = buffer.data(sni + 619);
    const auto *sni_621 = buffer.data(sni + 621);
    const auto *sni_622 = buffer.data(sni + 622);
    const auto *sni_625 = buffer.data(sni + 625);
    const auto *sni_626 = buffer.data(sni + 626);
    const auto *sni_630 = buffer.data(sni + 630);
    const auto *sni_636 = buffer.data(sni + 636);
    const auto *sni_637 = buffer.data(sni + 637);
    const auto *sni_638 = buffer.data(sni + 638);
    const auto *sni_639 = buffer.data(sni + 639);
    const auto *sni_640 = buffer.data(sni + 640);
    const auto *sni_641 = buffer.data(sni + 641);
    const auto *sni_642 = buffer.data(sni + 642);
    const auto *sni_643 = buffer.data(sni + 643);
    const auto *sni_644 = buffer.data(sni + 644);
    const auto *sni_646 = buffer.data(sni + 646);
    const auto *sni_647 = buffer.data(sni + 647);
    const auto *sni_649 = buffer.data(sni + 649);
    const auto *sni_650 = buffer.data(sni + 650);
    const auto *sni_653 = buffer.data(sni + 653);
    const auto *sni_654 = buffer.data(sni + 654);
    const auto *sni_656 = buffer.data(sni + 656);
    const auto *sni_658 = buffer.data(sni + 658);
    const auto *sni_659 = buffer.data(sni + 659);
    const auto *sni_661 = buffer.data(sni + 661);
    const auto *sni_662 = buffer.data(sni + 662);
    const auto *sni_664 = buffer.data(sni + 664);
    const auto *sni_665 = buffer.data(sni + 665);
    const auto *sni_666 = buffer.data(sni + 666);
    const auto *sni_667 = buffer.data(sni + 667);
    const auto *sni_668 = buffer.data(sni + 668);
    const auto *sni_669 = buffer.data(sni + 669);
    const auto *sni_670 = buffer.data(sni + 670);
    const auto *sni_671 = buffer.data(sni + 671);
    const auto *sni_672 = buffer.data(sni + 672);
    const auto *sni_674 = buffer.data(sni + 674);
    const auto *sni_675 = buffer.data(sni + 675);
    const auto *sni_677 = buffer.data(sni + 677);
    const auto *sni_678 = buffer.data(sni + 678);
    const auto *sni_681 = buffer.data(sni + 681);
    const auto *sni_682 = buffer.data(sni + 682);
    const auto *sni_684 = buffer.data(sni + 684);
    const auto *sni_686 = buffer.data(sni + 686);
    const auto *sni_687 = buffer.data(sni + 687);
    const auto *sni_689 = buffer.data(sni + 689);
    const auto *sni_690 = buffer.data(sni + 690);
    const auto *sni_692 = buffer.data(sni + 692);
    const auto *sni_693 = buffer.data(sni + 693);
    const auto *sni_694 = buffer.data(sni + 694);
    const auto *sni_695 = buffer.data(sni + 695);
    const auto *sni_696 = buffer.data(sni + 696);
    const auto *sni_697 = buffer.data(sni + 697);
    const auto *sni_698 = buffer.data(sni + 698);
    const auto *sni_699 = buffer.data(sni + 699);

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pb_z, pc_y, pc_z, smk0_540, smi_447, \
                         smi_448, smk1_540, snh0_461, snh1_461, sni_615, \
                         sni_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_21 * smi_447[k]
                   + f_3 * pc_y[k] * sni_615[k];

        t_791[k] = f_1 * snh0_461[k]
                   - f_2 * snh1_461[k]
                   + f_3 * pc_z[k] * sni_615[k];

        t_792[k] = pb_z[k] * smk0_540[k]
                   - f_12 * pc_z[k] * smk1_540[k];

        t_793[k] = f_17 * smi_448[k]
                   + f_3 * pc_y[k] * sni_616[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pb_z, pc_y, pc_z, smk0_543, smi_420, smi_450, \
                         smk1_543, sni_616, sni_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_13 * smi_420[k]
                   + f_3 * pc_z[k] * sni_616[k];

        t_795[k] = pb_z[k] * smk0_543[k]
                   - f_12 * pc_z[k] * smk1_543[k];

        t_796[k] = f_17 * smi_450[k]
                   + f_3 * pc_y[k] * sni_618[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pb_z, pc_x, pc_z, smk0_546, smi_423, smi_621, \
                         smk1_546, snh0_467, snh1_467, sni_619, \
                         sni_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_16 * smi_621[k]
                   + f_4 * snh0_467[k]
                   - f_5 * snh1_467[k]
                   + f_3 * pc_x[k] * sni_621[k];

        t_798[k] = pb_z[k] * smk0_546[k]
                   - f_12 * pc_z[k] * smk1_546[k];

        t_799[k] = f_13 * smi_423[k]
                   + f_3 * pc_z[k] * sni_619[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pb_z, pc_x, pc_y, pc_z, smk0_550, smi_453, \
                         smi_625, smk1_550, snh0_471, snh1_471, sni_621, \
                         sni_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_17 * smi_453[k]
                   + f_3 * pc_y[k] * sni_621[k];

        t_801[k] = f_16 * smi_625[k]
                   + f_6 * snh0_471[k]
                   - f_7 * snh1_471[k]
                   + f_3 * pc_x[k] * sni_625[k];

        t_802[k] = pb_z[k] * smk0_550[k]
                   - f_12 * pc_z[k] * smk1_550[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pb_z, pc_y, pc_z, smk0_552, smi_426, smi_427, \
                         smi_457, smk1_552, sni_622, sni_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_13 * smi_426[k]
                   + f_3 * pc_z[k] * sni_622[k];

        t_804[k] = pb_z[k] * smk0_552[k]
                   + f_14 * smi_427[k]
                   - f_12 * pc_z[k] * smk1_552[k];

        t_805[k] = f_17 * smi_457[k]
                   + f_3 * pc_y[k] * sni_625[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pb_z, pc_x, pc_z, smk0_555, smi_430, smi_630, \
                         smk1_555, snh0_476, snh1_476, sni_626, \
                         sni_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_16 * smi_630[k]
                   + f_8 * snh0_476[k]
                   - f_9 * snh1_476[k]
                   + f_3 * pc_x[k] * sni_630[k];

        t_807[k] = pb_z[k] * smk0_555[k]
                   - f_12 * pc_z[k] * smk1_555[k];

        t_808[k] = f_13 * smi_430[k]
                   + f_3 * pc_z[k] * sni_626[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pb_z, pc_y, pc_z, smk0_557, smk0_558, smi_431, \
                         smi_432, smi_462, smk1_557, smk1_558, \
                         sni_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pb_z[k] * smk0_557[k]
                   + f_14 * smi_431[k]
                   - f_12 * pc_z[k] * smk1_557[k];

        t_810[k] = pb_z[k] * smk0_558[k]
                   + f_15 * smi_432[k]
                   - f_12 * pc_z[k] * smk1_558[k];

        t_811[k] = f_17 * smi_462[k]
                   + f_3 * pc_y[k] * sni_630[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pc_x, smi_636, smi_637, smi_638, smi_639, \
                         snh0_482, snh1_482, sni_636, sni_637, sni_638, \
                         sni_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_16 * smi_636[k]
                   + f_10 * snh0_482[k]
                   - f_11 * snh1_482[k]
                   + f_3 * pc_x[k] * sni_636[k];

        t_813[k] = f_16 * smi_637[k]
                   + f_3 * pc_x[k] * sni_637[k];

        t_814[k] = f_16 * smi_638[k]
                   + f_3 * pc_x[k] * sni_638[k];

        t_815[k] = f_16 * smi_639[k]
                   + f_3 * pc_x[k] * sni_639[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, t_819, pc_x, smi_640, smi_641, smi_642, smi_643, \
                         sni_640, sni_641, sni_642, sni_643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_16 * smi_640[k]
                   + f_3 * pc_x[k] * sni_640[k];

        t_817[k] = f_16 * smi_641[k]
                   + f_3 * pc_x[k] * sni_641[k];

        t_818[k] = f_16 * smi_642[k]
                   + f_3 * pc_x[k] * sni_642[k];

        t_819[k] = f_16 * smi_643[k]
                   + f_3 * pc_x[k] * sni_643[k];
    }

#pragma omp simd aligned(t_820, t_821, t_822, pb_z, pc_y, pc_z, smk0_568, smi_441, smi_471, \
                         smk1_568, snh0_479, snh1_479, sni_637, \
                         sni_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = pb_z[k] * smk0_568[k]
                   - f_12 * pc_z[k] * smk1_568[k];

        t_821[k] = f_13 * smi_441[k]
                   + f_3 * pc_z[k] * sni_637[k];

        t_822[k] = f_17 * smi_471[k]
                   + f_4 * snh0_479[k]
                   - f_5 * snh1_479[k]
                   + f_3 * pc_y[k] * sni_639[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pc_y, smi_472, smi_473, smi_474, snh0_480, \
                         snh0_481, snh0_482, snh1_480, snh1_481, snh1_482, sni_640, sni_641, \
                         sni_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_17 * smi_472[k]
                   + f_6 * snh0_480[k]
                   - f_7 * snh1_480[k]
                   + f_3 * pc_y[k] * sni_640[k];

        t_824[k] = f_17 * smi_473[k]
                   + f_8 * snh0_481[k]
                   - f_9 * snh1_481[k]
                   + f_3 * pc_y[k] * sni_641[k];

        t_825[k] = f_17 * smi_474[k]
                   + f_10 * snh0_482[k]
                   - f_11 * snh1_482[k]
                   + f_3 * pc_y[k] * sni_642[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, pc_x, pc_y, pc_z, smi_447, smi_475, smi_644, \
                         snh0_482, snh0_483, snh1_482, snh1_483, sni_643, \
                         sni_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_17 * smi_475[k]
                   + f_3 * pc_y[k] * sni_643[k];

        t_827[k] = f_13 * smi_447[k]
                   + f_1 * snh0_482[k]
                   - f_2 * snh1_482[k]
                   + f_3 * pc_z[k] * sni_643[k];

        t_828[k] = f_16 * smi_644[k]
                   + f_1 * snh0_483[k]
                   - f_2 * snh1_483[k]
                   + f_3 * pc_x[k] * sni_644[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pc_x, pc_y, pc_z, smi_448, smi_476, \
                         smi_478, smi_647, snh0_486, snh1_486, sni_644, sni_646, \
                         sni_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_16 * smi_476[k]
                   + f_3 * pc_y[k] * sni_644[k];

        t_830[k] = f_14 * smi_448[k]
                   + f_3 * pc_z[k] * sni_644[k];

        t_831[k] = f_16 * smi_647[k]
                   + f_4 * snh0_486[k]
                   - f_5 * snh1_486[k]
                   + f_3 * pc_x[k] * sni_647[k];

        t_832[k] = f_16 * smi_478[k]
                   + f_3 * pc_y[k] * sni_646[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, pc_z, smi_451, smi_649, smi_650, snh0_488, \
                         snh0_489, snh1_488, snh1_489, sni_647, sni_649, \
                         sni_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_16 * smi_649[k]
                   + f_4 * snh0_488[k]
                   - f_5 * snh1_488[k]
                   + f_3 * pc_x[k] * sni_649[k];

        t_834[k] = f_16 * smi_650[k]
                   + f_6 * snh0_489[k]
                   - f_7 * snh1_489[k]
                   + f_3 * pc_x[k] * sni_650[k];

        t_835[k] = f_14 * smi_451[k]
                   + f_3 * pc_z[k] * sni_647[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pc_x, pc_y, smi_481, smi_653, smi_654, snh0_492, \
                         snh0_493, snh1_492, snh1_493, sni_649, sni_653, \
                         sni_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * smi_481[k]
                   + f_3 * pc_y[k] * sni_649[k];

        t_837[k] = f_16 * smi_653[k]
                   + f_6 * snh0_492[k]
                   - f_7 * snh1_492[k]
                   + f_3 * pc_x[k] * sni_653[k];

        t_838[k] = f_16 * smi_654[k]
                   + f_8 * snh0_493[k]
                   - f_9 * snh1_493[k]
                   + f_3 * pc_x[k] * sni_654[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pc_x, pc_y, pc_z, smi_454, smi_485, smi_656, \
                         snh0_495, snh1_495, sni_650, sni_653, \
                         sni_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_14 * smi_454[k]
                   + f_3 * pc_z[k] * sni_650[k];

        t_840[k] = f_16 * smi_656[k]
                   + f_8 * snh0_495[k]
                   - f_9 * snh1_495[k]
                   + f_3 * pc_x[k] * sni_656[k];

        t_841[k] = f_16 * smi_485[k]
                   + f_3 * pc_y[k] * sni_653[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, pc_z, smi_458, smi_658, smi_659, snh0_497, \
                         snh0_498, snh1_497, snh1_498, sni_654, sni_658, \
                         sni_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_16 * smi_658[k]
                   + f_8 * snh0_497[k]
                   - f_9 * snh1_497[k]
                   + f_3 * pc_x[k] * sni_658[k];

        t_843[k] = f_16 * smi_659[k]
                   + f_10 * snh0_498[k]
                   - f_11 * snh1_498[k]
                   + f_3 * pc_x[k] * sni_659[k];

        t_844[k] = f_14 * smi_458[k]
                   + f_3 * pc_z[k] * sni_654[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, smi_490, smi_661, smi_662, snh0_500, \
                         snh0_501, snh1_500, snh1_501, sni_658, sni_661, \
                         sni_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_16 * smi_661[k]
                   + f_10 * snh0_500[k]
                   - f_11 * snh1_500[k]
                   + f_3 * pc_x[k] * sni_661[k];

        t_846[k] = f_16 * smi_662[k]
                   + f_10 * snh0_501[k]
                   - f_11 * snh1_501[k]
                   + f_3 * pc_x[k] * sni_662[k];

        t_847[k] = f_16 * smi_490[k]
                   + f_3 * pc_y[k] * sni_658[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, t_851, pc_x, smi_664, smi_665, smi_666, smi_667, \
                         snh0_503, snh1_503, sni_664, sni_665, sni_666, \
                         sni_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * smi_664[k]
                   + f_10 * snh0_503[k]
                   - f_11 * snh1_503[k]
                   + f_3 * pc_x[k] * sni_664[k];

        t_849[k] = f_16 * smi_665[k]
                   + f_3 * pc_x[k] * sni_665[k];

        t_850[k] = f_16 * smi_666[k]
                   + f_3 * pc_x[k] * sni_666[k];

        t_851[k] = f_16 * smi_667[k]
                   + f_3 * pc_x[k] * sni_667[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, t_855, pc_x, smi_668, smi_669, smi_670, smi_671, \
                         sni_668, sni_669, sni_670, sni_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_16 * smi_668[k]
                   + f_3 * pc_x[k] * sni_668[k];

        t_853[k] = f_16 * smi_669[k]
                   + f_3 * pc_x[k] * sni_669[k];

        t_854[k] = f_16 * smi_670[k]
                   + f_3 * pc_x[k] * sni_670[k];

        t_855[k] = f_16 * smi_671[k]
                   + f_3 * pc_x[k] * sni_671[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, pc_y, pc_z, smi_469, smi_497, smi_499, snh0_498, \
                         snh0_500, snh1_498, snh1_500, sni_665, \
                         sni_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_16 * smi_497[k]
                   + f_1 * snh0_498[k]
                   - f_2 * snh1_498[k]
                   + f_3 * pc_y[k] * sni_665[k];

        t_857[k] = f_14 * smi_469[k]
                   + f_3 * pc_z[k] * sni_665[k];

        t_858[k] = f_16 * smi_499[k]
                   + f_4 * snh0_500[k]
                   - f_5 * snh1_500[k]
                   + f_3 * pc_y[k] * sni_667[k];
    }

#pragma omp simd aligned(t_859, t_860, t_861, pc_y, smi_500, smi_501, smi_502, snh0_501, \
                         snh0_502, snh0_503, snh1_501, snh1_502, snh1_503, sni_668, sni_669, \
                         sni_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_859[k] = f_16 * smi_500[k]
                   + f_6 * snh0_501[k]
                   - f_7 * snh1_501[k]
                   + f_3 * pc_y[k] * sni_668[k];

        t_860[k] = f_16 * smi_501[k]
                   + f_8 * snh0_502[k]
                   - f_9 * snh1_502[k]
                   + f_3 * pc_y[k] * sni_669[k];

        t_861[k] = f_16 * smi_502[k]
                   + f_10 * snh0_503[k]
                   - f_11 * snh1_503[k]
                   + f_3 * pc_y[k] * sni_670[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, pc_x, pc_y, pc_z, smi_475, smi_503, smi_672, \
                         snh0_503, snh0_504, snh1_503, snh1_504, sni_671, \
                         sni_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * smi_503[k]
                   + f_3 * pc_y[k] * sni_671[k];

        t_863[k] = f_14 * smi_475[k]
                   + f_1 * snh0_503[k]
                   - f_2 * snh1_503[k]
                   + f_3 * pc_z[k] * sni_671[k];

        t_864[k] = f_16 * smi_672[k]
                   + f_1 * snh0_504[k]
                   - f_2 * snh1_504[k]
                   + f_3 * pc_x[k] * sni_672[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pc_x, pc_y, pc_z, smi_476, smi_504, \
                         smi_506, smi_675, snh0_507, snh1_507, sni_672, sni_674, \
                         sni_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_15 * smi_504[k]
                   + f_3 * pc_y[k] * sni_672[k];

        t_866[k] = f_15 * smi_476[k]
                   + f_3 * pc_z[k] * sni_672[k];

        t_867[k] = f_16 * smi_675[k]
                   + f_4 * snh0_507[k]
                   - f_5 * snh1_507[k]
                   + f_3 * pc_x[k] * sni_675[k];

        t_868[k] = f_15 * smi_506[k]
                   + f_3 * pc_y[k] * sni_674[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pc_x, pc_z, smi_479, smi_677, smi_678, snh0_509, \
                         snh0_510, snh1_509, snh1_510, sni_675, sni_677, \
                         sni_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_16 * smi_677[k]
                   + f_4 * snh0_509[k]
                   - f_5 * snh1_509[k]
                   + f_3 * pc_x[k] * sni_677[k];

        t_870[k] = f_16 * smi_678[k]
                   + f_6 * snh0_510[k]
                   - f_7 * snh1_510[k]
                   + f_3 * pc_x[k] * sni_678[k];

        t_871[k] = f_15 * smi_479[k]
                   + f_3 * pc_z[k] * sni_675[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pc_x, pc_y, smi_509, smi_681, smi_682, snh0_513, \
                         snh0_514, snh1_513, snh1_514, sni_677, sni_681, \
                         sni_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_15 * smi_509[k]
                   + f_3 * pc_y[k] * sni_677[k];

        t_873[k] = f_16 * smi_681[k]
                   + f_6 * snh0_513[k]
                   - f_7 * snh1_513[k]
                   + f_3 * pc_x[k] * sni_681[k];

        t_874[k] = f_16 * smi_682[k]
                   + f_8 * snh0_514[k]
                   - f_9 * snh1_514[k]
                   + f_3 * pc_x[k] * sni_682[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pc_x, pc_y, pc_z, smi_482, smi_513, smi_684, \
                         snh0_516, snh1_516, sni_678, sni_681, \
                         sni_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_15 * smi_482[k]
                   + f_3 * pc_z[k] * sni_678[k];

        t_876[k] = f_16 * smi_684[k]
                   + f_8 * snh0_516[k]
                   - f_9 * snh1_516[k]
                   + f_3 * pc_x[k] * sni_684[k];

        t_877[k] = f_15 * smi_513[k]
                   + f_3 * pc_y[k] * sni_681[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_x, pc_z, smi_486, smi_686, smi_687, snh0_518, \
                         snh0_519, snh1_518, snh1_519, sni_682, sni_686, \
                         sni_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_16 * smi_686[k]
                   + f_8 * snh0_518[k]
                   - f_9 * snh1_518[k]
                   + f_3 * pc_x[k] * sni_686[k];

        t_879[k] = f_16 * smi_687[k]
                   + f_10 * snh0_519[k]
                   - f_11 * snh1_519[k]
                   + f_3 * pc_x[k] * sni_687[k];

        t_880[k] = f_15 * smi_486[k]
                   + f_3 * pc_z[k] * sni_682[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pc_x, pc_y, smi_518, smi_689, smi_690, snh0_521, \
                         snh0_522, snh1_521, snh1_522, sni_686, sni_689, \
                         sni_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_16 * smi_689[k]
                   + f_10 * snh0_521[k]
                   - f_11 * snh1_521[k]
                   + f_3 * pc_x[k] * sni_689[k];

        t_882[k] = f_16 * smi_690[k]
                   + f_10 * snh0_522[k]
                   - f_11 * snh1_522[k]
                   + f_3 * pc_x[k] * sni_690[k];

        t_883[k] = f_15 * smi_518[k]
                   + f_3 * pc_y[k] * sni_686[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, pc_x, smi_692, smi_693, smi_694, smi_695, \
                         snh0_524, snh1_524, sni_692, sni_693, sni_694, \
                         sni_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_16 * smi_692[k]
                   + f_10 * snh0_524[k]
                   - f_11 * snh1_524[k]
                   + f_3 * pc_x[k] * sni_692[k];

        t_885[k] = f_16 * smi_693[k]
                   + f_3 * pc_x[k] * sni_693[k];

        t_886[k] = f_16 * smi_694[k]
                   + f_3 * pc_x[k] * sni_694[k];

        t_887[k] = f_16 * smi_695[k]
                   + f_3 * pc_x[k] * sni_695[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pc_x, smi_696, smi_697, smi_698, smi_699, \
                         sni_696, sni_697, sni_698, sni_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_16 * smi_696[k]
                   + f_3 * pc_x[k] * sni_696[k];

        t_889[k] = f_16 * smi_697[k]
                   + f_3 * pc_x[k] * sni_697[k];

        t_890[k] = f_16 * smi_698[k]
                   + f_3 * pc_x[k] * sni_698[k];

        t_891[k] = f_16 * smi_699[k]
                   + f_3 * pc_x[k] * sni_699[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, pc_y, pc_z, smi_497, smi_525, smi_527, snh0_519, \
                         snh0_521, snh1_519, snh1_521, sni_693, \
                         sni_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_15 * smi_525[k]
                   + f_1 * snh0_519[k]
                   - f_2 * snh1_519[k]
                   + f_3 * pc_y[k] * sni_693[k];

        t_893[k] = f_15 * smi_497[k]
                   + f_3 * pc_z[k] * sni_693[k];

        t_894[k] = f_15 * smi_527[k]
                   + f_4 * snh0_521[k]
                   - f_5 * snh1_521[k]
                   + f_3 * pc_y[k] * sni_695[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, pc_y, smi_528, smi_529, smi_530, snh0_522, \
                         snh0_523, snh0_524, snh1_522, snh1_523, snh1_524, sni_696, sni_697, \
                         sni_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_15 * smi_528[k]
                   + f_6 * snh0_522[k]
                   - f_7 * snh1_522[k]
                   + f_3 * pc_y[k] * sni_696[k];

        t_896[k] = f_15 * smi_529[k]
                   + f_8 * snh0_523[k]
                   - f_9 * snh1_523[k]
                   + f_3 * pc_y[k] * sni_697[k];

        t_897[k] = f_15 * smi_530[k]
                   + f_10 * snh0_524[k]
                   - f_11 * snh1_524[k]
                   + f_3 * pc_y[k] * sni_698[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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
    const auto f_21 = 3.0 / q;

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

    const auto *smk0_720 = buffer.data(smk0 + 720);
    const auto *smk0_723 = buffer.data(smk0 + 723);
    const auto *smk0_725 = buffer.data(smk0 + 725);
    const auto *smk0_726 = buffer.data(smk0 + 726);
    const auto *smk0_729 = buffer.data(smk0 + 729);
    const auto *smk0_730 = buffer.data(smk0 + 730);
    const auto *smk0_732 = buffer.data(smk0 + 732);
    const auto *smk0_734 = buffer.data(smk0 + 734);
    const auto *smk0_735 = buffer.data(smk0 + 735);
    const auto *smk0_737 = buffer.data(smk0 + 737);
    const auto *smk0_738 = buffer.data(smk0 + 738);
    const auto *smk0_740 = buffer.data(smk0 + 740);
    const auto *smk0_755 = buffer.data(smk0 + 755);

    const auto *smi_503 = buffer.data(smi + 503);
    const auto *smi_504 = buffer.data(smi + 504);
    const auto *smi_507 = buffer.data(smi + 507);
    const auto *smi_510 = buffer.data(smi + 510);
    const auto *smi_514 = buffer.data(smi + 514);
    const auto *smi_525 = buffer.data(smi + 525);
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
    const auto *smi_555 = buffer.data(smi + 555);
    const auto *smi_556 = buffer.data(smi + 556);
    const auto *smi_557 = buffer.data(smi + 557);
    const auto *smi_558 = buffer.data(smi + 558);
    const auto *smi_559 = buffer.data(smi + 559);
    const auto *smi_560 = buffer.data(smi + 560);
    const auto *smi_561 = buffer.data(smi + 561);
    const auto *smi_562 = buffer.data(smi + 562);
    const auto *smi_563 = buffer.data(smi + 563);
    const auto *smi_565 = buffer.data(smi + 565);
    const auto *smi_566 = buffer.data(smi + 566);
    const auto *smi_568 = buffer.data(smi + 568);
    const auto *smi_569 = buffer.data(smi + 569);
    const auto *smi_570 = buffer.data(smi + 570);
    const auto *smi_572 = buffer.data(smi + 572);
    const auto *smi_573 = buffer.data(smi + 573);
    const auto *smi_574 = buffer.data(smi + 574);
    const auto *smi_581 = buffer.data(smi + 581);
    const auto *smi_583 = buffer.data(smi + 583);
    const auto *smi_584 = buffer.data(smi + 584);
    const auto *smi_585 = buffer.data(smi + 585);
    const auto *smi_586 = buffer.data(smi + 586);
    const auto *smi_587 = buffer.data(smi + 587);
    const auto *smi_700 = buffer.data(smi + 700);
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
    const auto *smi_749 = buffer.data(smi + 749);
    const auto *smi_750 = buffer.data(smi + 750);
    const auto *smi_751 = buffer.data(smi + 751);
    const auto *smi_752 = buffer.data(smi + 752);
    const auto *smi_753 = buffer.data(smi + 753);
    const auto *smi_754 = buffer.data(smi + 754);
    const auto *smi_755 = buffer.data(smi + 755);
    const auto *smi_756 = buffer.data(smi + 756);
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

    const auto *smk1_720 = buffer.data(smk1 + 720);
    const auto *smk1_723 = buffer.data(smk1 + 723);
    const auto *smk1_725 = buffer.data(smk1 + 725);
    const auto *smk1_726 = buffer.data(smk1 + 726);
    const auto *smk1_729 = buffer.data(smk1 + 729);
    const auto *smk1_730 = buffer.data(smk1 + 730);
    const auto *smk1_732 = buffer.data(smk1 + 732);
    const auto *smk1_734 = buffer.data(smk1 + 734);
    const auto *smk1_735 = buffer.data(smk1 + 735);
    const auto *smk1_737 = buffer.data(smk1 + 737);
    const auto *smk1_738 = buffer.data(smk1 + 738);
    const auto *smk1_740 = buffer.data(smk1 + 740);
    const auto *smk1_755 = buffer.data(smk1 + 755);

    const auto *snh0_524 = buffer.data(snh0 + 524);
    const auto *snh0_525 = buffer.data(snh0 + 525);
    const auto *snh0_528 = buffer.data(snh0 + 528);
    const auto *snh0_530 = buffer.data(snh0 + 530);
    const auto *snh0_531 = buffer.data(snh0 + 531);
    const auto *snh0_534 = buffer.data(snh0 + 534);
    const auto *snh0_535 = buffer.data(snh0 + 535);
    const auto *snh0_537 = buffer.data(snh0 + 537);
    const auto *snh0_539 = buffer.data(snh0 + 539);
    const auto *snh0_540 = buffer.data(snh0 + 540);
    const auto *snh0_542 = buffer.data(snh0 + 542);
    const auto *snh0_543 = buffer.data(snh0 + 543);
    const auto *snh0_544 = buffer.data(snh0 + 544);
    const auto *snh0_545 = buffer.data(snh0 + 545);
    const auto *snh0_561 = buffer.data(snh0 + 561);
    const auto *snh0_563 = buffer.data(snh0 + 563);
    const auto *snh0_564 = buffer.data(snh0 + 564);
    const auto *snh0_565 = buffer.data(snh0 + 565);
    const auto *snh0_566 = buffer.data(snh0 + 566);
    const auto *snh0_567 = buffer.data(snh0 + 567);
    const auto *snh0_570 = buffer.data(snh0 + 570);
    const auto *snh0_572 = buffer.data(snh0 + 572);
    const auto *snh0_573 = buffer.data(snh0 + 573);
    const auto *snh0_576 = buffer.data(snh0 + 576);
    const auto *snh0_577 = buffer.data(snh0 + 577);
    const auto *snh0_579 = buffer.data(snh0 + 579);
    const auto *snh0_581 = buffer.data(snh0 + 581);
    const auto *snh0_582 = buffer.data(snh0 + 582);
    const auto *snh0_584 = buffer.data(snh0 + 584);
    const auto *snh0_585 = buffer.data(snh0 + 585);
    const auto *snh0_586 = buffer.data(snh0 + 586);
    const auto *snh0_587 = buffer.data(snh0 + 587);

    const auto *snh1_524 = buffer.data(snh1 + 524);
    const auto *snh1_525 = buffer.data(snh1 + 525);
    const auto *snh1_528 = buffer.data(snh1 + 528);
    const auto *snh1_530 = buffer.data(snh1 + 530);
    const auto *snh1_531 = buffer.data(snh1 + 531);
    const auto *snh1_534 = buffer.data(snh1 + 534);
    const auto *snh1_535 = buffer.data(snh1 + 535);
    const auto *snh1_537 = buffer.data(snh1 + 537);
    const auto *snh1_539 = buffer.data(snh1 + 539);
    const auto *snh1_540 = buffer.data(snh1 + 540);
    const auto *snh1_542 = buffer.data(snh1 + 542);
    const auto *snh1_543 = buffer.data(snh1 + 543);
    const auto *snh1_544 = buffer.data(snh1 + 544);
    const auto *snh1_545 = buffer.data(snh1 + 545);
    const auto *snh1_561 = buffer.data(snh1 + 561);
    const auto *snh1_563 = buffer.data(snh1 + 563);
    const auto *snh1_564 = buffer.data(snh1 + 564);
    const auto *snh1_565 = buffer.data(snh1 + 565);
    const auto *snh1_566 = buffer.data(snh1 + 566);
    const auto *snh1_567 = buffer.data(snh1 + 567);
    const auto *snh1_570 = buffer.data(snh1 + 570);
    const auto *snh1_572 = buffer.data(snh1 + 572);
    const auto *snh1_573 = buffer.data(snh1 + 573);
    const auto *snh1_576 = buffer.data(snh1 + 576);
    const auto *snh1_577 = buffer.data(snh1 + 577);
    const auto *snh1_579 = buffer.data(snh1 + 579);
    const auto *snh1_581 = buffer.data(snh1 + 581);
    const auto *snh1_582 = buffer.data(snh1 + 582);
    const auto *snh1_584 = buffer.data(snh1 + 584);
    const auto *snh1_585 = buffer.data(snh1 + 585);
    const auto *snh1_586 = buffer.data(snh1 + 586);
    const auto *snh1_587 = buffer.data(snh1 + 587);

    const auto *sni_699 = buffer.data(sni + 699);
    const auto *sni_700 = buffer.data(sni + 700);
    const auto *sni_702 = buffer.data(sni + 702);
    const auto *sni_703 = buffer.data(sni + 703);
    const auto *sni_705 = buffer.data(sni + 705);
    const auto *sni_706 = buffer.data(sni + 706);
    const auto *sni_709 = buffer.data(sni + 709);
    const auto *sni_710 = buffer.data(sni + 710);
    const auto *sni_712 = buffer.data(sni + 712);
    const auto *sni_714 = buffer.data(sni + 714);
    const auto *sni_715 = buffer.data(sni + 715);
    const auto *sni_717 = buffer.data(sni + 717);
    const auto *sni_718 = buffer.data(sni + 718);
    const auto *sni_720 = buffer.data(sni + 720);
    const auto *sni_721 = buffer.data(sni + 721);
    const auto *sni_722 = buffer.data(sni + 722);
    const auto *sni_723 = buffer.data(sni + 723);
    const auto *sni_724 = buffer.data(sni + 724);
    const auto *sni_725 = buffer.data(sni + 725);
    const auto *sni_726 = buffer.data(sni + 726);
    const auto *sni_727 = buffer.data(sni + 727);
    const auto *sni_728 = buffer.data(sni + 728);
    const auto *sni_730 = buffer.data(sni + 730);
    const auto *sni_731 = buffer.data(sni + 731);
    const auto *sni_733 = buffer.data(sni + 733);
    const auto *sni_734 = buffer.data(sni + 734);
    const auto *sni_737 = buffer.data(sni + 737);
    const auto *sni_738 = buffer.data(sni + 738);
    const auto *sni_742 = buffer.data(sni + 742);
    const auto *sni_749 = buffer.data(sni + 749);
    const auto *sni_750 = buffer.data(sni + 750);
    const auto *sni_751 = buffer.data(sni + 751);
    const auto *sni_752 = buffer.data(sni + 752);
    const auto *sni_753 = buffer.data(sni + 753);
    const auto *sni_754 = buffer.data(sni + 754);
    const auto *sni_755 = buffer.data(sni + 755);
    const auto *sni_756 = buffer.data(sni + 756);
    const auto *sni_758 = buffer.data(sni + 758);
    const auto *sni_759 = buffer.data(sni + 759);
    const auto *sni_761 = buffer.data(sni + 761);
    const auto *sni_762 = buffer.data(sni + 762);
    const auto *sni_765 = buffer.data(sni + 765);
    const auto *sni_766 = buffer.data(sni + 766);
    const auto *sni_768 = buffer.data(sni + 768);
    const auto *sni_770 = buffer.data(sni + 770);
    const auto *sni_771 = buffer.data(sni + 771);
    const auto *sni_773 = buffer.data(sni + 773);
    const auto *sni_774 = buffer.data(sni + 774);
    const auto *sni_776 = buffer.data(sni + 776);
    const auto *sni_777 = buffer.data(sni + 777);
    const auto *sni_778 = buffer.data(sni + 778);
    const auto *sni_779 = buffer.data(sni + 779);
    const auto *sni_780 = buffer.data(sni + 780);
    const auto *sni_781 = buffer.data(sni + 781);
    const auto *sni_782 = buffer.data(sni + 782);
    const auto *sni_783 = buffer.data(sni + 783);

#pragma omp simd aligned(t_898, t_899, t_900, pc_x, pc_y, pc_z, smi_503, smi_531, smi_700, \
                         snh0_524, snh0_525, snh1_524, snh1_525, sni_699, \
                         sni_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * smi_531[k]
                   + f_3 * pc_y[k] * sni_699[k];

        t_899[k] = f_15 * smi_503[k]
                   + f_1 * snh0_524[k]
                   - f_2 * snh1_524[k]
                   + f_3 * pc_z[k] * sni_699[k];

        t_900[k] = f_16 * smi_700[k]
                   + f_1 * snh0_525[k]
                   - f_2 * snh1_525[k]
                   + f_3 * pc_x[k] * sni_700[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, t_904, pc_x, pc_y, pc_z, smi_504, smi_532, \
                         smi_534, smi_703, snh0_528, snh1_528, sni_700, sni_702, \
                         sni_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = f_14 * smi_532[k]
                   + f_3 * pc_y[k] * sni_700[k];

        t_902[k] = f_16 * smi_504[k]
                   + f_3 * pc_z[k] * sni_700[k];

        t_903[k] = f_16 * smi_703[k]
                   + f_4 * snh0_528[k]
                   - f_5 * snh1_528[k]
                   + f_3 * pc_x[k] * sni_703[k];

        t_904[k] = f_14 * smi_534[k]
                   + f_3 * pc_y[k] * sni_702[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, smi_507, smi_705, smi_706, snh0_530, \
                         snh0_531, snh1_530, snh1_531, sni_703, sni_705, \
                         sni_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_16 * smi_705[k]
                   + f_4 * snh0_530[k]
                   - f_5 * snh1_530[k]
                   + f_3 * pc_x[k] * sni_705[k];

        t_906[k] = f_16 * smi_706[k]
                   + f_6 * snh0_531[k]
                   - f_7 * snh1_531[k]
                   + f_3 * pc_x[k] * sni_706[k];

        t_907[k] = f_16 * smi_507[k]
                   + f_3 * pc_z[k] * sni_703[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, smi_537, smi_709, smi_710, snh0_534, \
                         snh0_535, snh1_534, snh1_535, sni_705, sni_709, \
                         sni_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * smi_537[k]
                   + f_3 * pc_y[k] * sni_705[k];

        t_909[k] = f_16 * smi_709[k]
                   + f_6 * snh0_534[k]
                   - f_7 * snh1_534[k]
                   + f_3 * pc_x[k] * sni_709[k];

        t_910[k] = f_16 * smi_710[k]
                   + f_8 * snh0_535[k]
                   - f_9 * snh1_535[k]
                   + f_3 * pc_x[k] * sni_710[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, pc_x, pc_y, pc_z, smi_510, smi_541, smi_712, \
                         snh0_537, snh1_537, sni_706, sni_709, \
                         sni_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_16 * smi_510[k]
                   + f_3 * pc_z[k] * sni_706[k];

        t_912[k] = f_16 * smi_712[k]
                   + f_8 * snh0_537[k]
                   - f_9 * snh1_537[k]
                   + f_3 * pc_x[k] * sni_712[k];

        t_913[k] = f_14 * smi_541[k]
                   + f_3 * pc_y[k] * sni_709[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pc_x, pc_z, smi_514, smi_714, smi_715, snh0_539, \
                         snh0_540, snh1_539, snh1_540, sni_710, sni_714, \
                         sni_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_16 * smi_714[k]
                   + f_8 * snh0_539[k]
                   - f_9 * snh1_539[k]
                   + f_3 * pc_x[k] * sni_714[k];

        t_915[k] = f_16 * smi_715[k]
                   + f_10 * snh0_540[k]
                   - f_11 * snh1_540[k]
                   + f_3 * pc_x[k] * sni_715[k];

        t_916[k] = f_16 * smi_514[k]
                   + f_3 * pc_z[k] * sni_710[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pc_x, pc_y, smi_546, smi_717, smi_718, snh0_542, \
                         snh0_543, snh1_542, snh1_543, sni_714, sni_717, \
                         sni_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_16 * smi_717[k]
                   + f_10 * snh0_542[k]
                   - f_11 * snh1_542[k]
                   + f_3 * pc_x[k] * sni_717[k];

        t_918[k] = f_16 * smi_718[k]
                   + f_10 * snh0_543[k]
                   - f_11 * snh1_543[k]
                   + f_3 * pc_x[k] * sni_718[k];

        t_919[k] = f_14 * smi_546[k]
                   + f_3 * pc_y[k] * sni_714[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pc_x, smi_720, smi_721, smi_722, smi_723, \
                         snh0_545, snh1_545, sni_720, sni_721, sni_722, \
                         sni_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = f_16 * smi_720[k]
                   + f_10 * snh0_545[k]
                   - f_11 * snh1_545[k]
                   + f_3 * pc_x[k] * sni_720[k];

        t_921[k] = f_16 * smi_721[k]
                   + f_3 * pc_x[k] * sni_721[k];

        t_922[k] = f_16 * smi_722[k]
                   + f_3 * pc_x[k] * sni_722[k];

        t_923[k] = f_16 * smi_723[k]
                   + f_3 * pc_x[k] * sni_723[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pc_x, smi_724, smi_725, smi_726, smi_727, \
                         sni_724, sni_725, sni_726, sni_727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_16 * smi_724[k]
                   + f_3 * pc_x[k] * sni_724[k];

        t_925[k] = f_16 * smi_725[k]
                   + f_3 * pc_x[k] * sni_725[k];

        t_926[k] = f_16 * smi_726[k]
                   + f_3 * pc_x[k] * sni_726[k];

        t_927[k] = f_16 * smi_727[k]
                   + f_3 * pc_x[k] * sni_727[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_y, pc_z, smi_525, smi_553, smi_555, snh0_540, \
                         snh0_542, snh1_540, snh1_542, sni_721, \
                         sni_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_14 * smi_553[k]
                   + f_1 * snh0_540[k]
                   - f_2 * snh1_540[k]
                   + f_3 * pc_y[k] * sni_721[k];

        t_929[k] = f_16 * smi_525[k]
                   + f_3 * pc_z[k] * sni_721[k];

        t_930[k] = f_14 * smi_555[k]
                   + f_4 * snh0_542[k]
                   - f_5 * snh1_542[k]
                   + f_3 * pc_y[k] * sni_723[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, pc_y, smi_556, smi_557, smi_558, snh0_543, \
                         snh0_544, snh0_545, snh1_543, snh1_544, snh1_545, sni_724, sni_725, \
                         sni_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_14 * smi_556[k]
                   + f_6 * snh0_543[k]
                   - f_7 * snh1_543[k]
                   + f_3 * pc_y[k] * sni_724[k];

        t_932[k] = f_14 * smi_557[k]
                   + f_8 * snh0_544[k]
                   - f_9 * snh1_544[k]
                   + f_3 * pc_y[k] * sni_725[k];

        t_933[k] = f_14 * smi_558[k]
                   + f_10 * snh0_545[k]
                   - f_11 * snh1_545[k]
                   + f_3 * pc_y[k] * sni_726[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, pb_y, pc_y, pc_z, smk0_720, smi_531, \
                         smi_559, smi_560, smk1_720, snh0_545, snh1_545, sni_727, \
                         sni_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_14 * smi_559[k]
                   + f_3 * pc_y[k] * sni_727[k];

        t_935[k] = f_16 * smi_531[k]
                   + f_1 * snh0_545[k]
                   - f_2 * snh1_545[k]
                   + f_3 * pc_z[k] * sni_727[k];

        t_936[k] = pb_y[k] * smk0_720[k]
                   - f_12 * pc_y[k] * smk1_720[k];

        t_937[k] = f_13 * smi_560[k]
                   + f_3 * pc_y[k] * sni_728[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pb_y, pc_y, pc_z, smk0_723, smk0_725, \
                         smi_532, smi_561, smi_562, smk1_723, smk1_725, sni_728, \
                         sni_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_17 * smi_532[k]
                   + f_3 * pc_z[k] * sni_728[k];

        t_939[k] = pb_y[k] * smk0_723[k]
                   + f_14 * smi_561[k]
                   - f_12 * pc_y[k] * smk1_723[k];

        t_940[k] = f_13 * smi_562[k]
                   + f_3 * pc_y[k] * sni_730[k];

        t_941[k] = pb_y[k] * smk0_725[k]
                   - f_12 * pc_y[k] * smk1_725[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pb_y, pc_y, pc_z, smk0_726, smk0_729, \
                         smi_535, smi_563, smi_565, smk1_726, smk1_729, sni_731, \
                         sni_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = pb_y[k] * smk0_726[k]
                   + f_15 * smi_563[k]
                   - f_12 * pc_y[k] * smk1_726[k];

        t_943[k] = f_17 * smi_535[k]
                   + f_3 * pc_z[k] * sni_731[k];

        t_944[k] = f_13 * smi_565[k]
                   + f_3 * pc_y[k] * sni_733[k];

        t_945[k] = pb_y[k] * smk0_729[k]
                   - f_12 * pc_y[k] * smk1_729[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pb_y, pc_y, pc_z, smk0_730, smk0_732, smi_538, \
                         smi_566, smi_568, smk1_730, smk1_732, \
                         sni_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = pb_y[k] * smk0_730[k]
                   + f_16 * smi_566[k]
                   - f_12 * pc_y[k] * smk1_730[k];

        t_947[k] = f_17 * smi_538[k]
                   + f_3 * pc_z[k] * sni_734[k];

        t_948[k] = pb_y[k] * smk0_732[k]
                   + f_14 * smi_568[k]
                   - f_12 * pc_y[k] * smk1_732[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pb_y, pc_y, pc_z, smk0_734, smk0_735, \
                         smi_542, smi_569, smi_570, smk1_734, smk1_735, sni_737, \
                         sni_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_13 * smi_569[k]
                   + f_3 * pc_y[k] * sni_737[k];

        t_950[k] = pb_y[k] * smk0_734[k]
                   - f_12 * pc_y[k] * smk1_734[k];

        t_951[k] = pb_y[k] * smk0_735[k]
                   + f_17 * smi_570[k]
                   - f_12 * pc_y[k] * smk1_735[k];

        t_952[k] = f_17 * smi_542[k]
                   + f_3 * pc_z[k] * sni_738[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pb_y, pc_y, smk0_737, smk0_738, smk0_740, \
                         smi_572, smi_573, smi_574, smk1_737, smk1_738, smk1_740, \
                         sni_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = pb_y[k] * smk0_737[k]
                   + f_15 * smi_572[k]
                   - f_12 * pc_y[k] * smk1_737[k];

        t_954[k] = pb_y[k] * smk0_738[k]
                   + f_14 * smi_573[k]
                   - f_12 * pc_y[k] * smk1_738[k];

        t_955[k] = f_13 * smi_574[k]
                   + f_3 * pc_y[k] * sni_742[k];

        t_956[k] = pb_y[k] * smk0_740[k]
                   - f_12 * pc_y[k] * smk1_740[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, t_961, pc_x, smi_749, smi_750, smi_751, \
                         smi_752, smi_753, sni_749, sni_750, sni_751, sni_752, \
                         sni_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_16 * smi_749[k]
                   + f_3 * pc_x[k] * sni_749[k];

        t_958[k] = f_16 * smi_750[k]
                   + f_3 * pc_x[k] * sni_750[k];

        t_959[k] = f_16 * smi_751[k]
                   + f_3 * pc_x[k] * sni_751[k];

        t_960[k] = f_16 * smi_752[k]
                   + f_3 * pc_x[k] * sni_752[k];

        t_961[k] = f_16 * smi_753[k]
                   + f_3 * pc_x[k] * sni_753[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pc_x, pc_y, pc_z, smi_553, smi_581, \
                         smi_754, smi_755, snh0_561, snh1_561, sni_749, sni_754, \
                         sni_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_16 * smi_754[k]
                   + f_3 * pc_x[k] * sni_754[k];

        t_963[k] = f_16 * smi_755[k]
                   + f_3 * pc_x[k] * sni_755[k];

        t_964[k] = f_13 * smi_581[k]
                   + f_1 * snh0_561[k]
                   - f_2 * snh1_561[k]
                   + f_3 * pc_y[k] * sni_749[k];

        t_965[k] = f_17 * smi_553[k]
                   + f_3 * pc_z[k] * sni_749[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_y, smi_583, smi_584, smi_585, snh0_563, \
                         snh0_564, snh0_565, snh1_563, snh1_564, snh1_565, sni_751, sni_752, \
                         sni_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_13 * smi_583[k]
                   + f_4 * snh0_563[k]
                   - f_5 * snh1_563[k]
                   + f_3 * pc_y[k] * sni_751[k];

        t_967[k] = f_13 * smi_584[k]
                   + f_6 * snh0_564[k]
                   - f_7 * snh1_564[k]
                   + f_3 * pc_y[k] * sni_752[k];

        t_968[k] = f_13 * smi_585[k]
                   + f_8 * snh0_565[k]
                   - f_9 * snh1_565[k]
                   + f_3 * pc_y[k] * sni_753[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pb_y, pc_y, smk0_755, smi_586, smi_587, \
                         smk1_755, snh0_566, snh1_566, sni_754, \
                         sni_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_13 * smi_586[k]
                   + f_10 * snh0_566[k]
                   - f_11 * snh1_566[k]
                   + f_3 * pc_y[k] * sni_754[k];

        t_970[k] = f_13 * smi_587[k]
                   + f_3 * pc_y[k] * sni_755[k];

        t_971[k] = pb_y[k] * smk0_755[k]
                   - f_12 * pc_y[k] * smk1_755[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, pc_x, pc_y, pc_z, smi_560, smi_756, \
                         smi_759, snh0_567, snh0_570, snh1_567, snh1_570, sni_756, \
                         sni_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_16 * smi_756[k]
                   + f_1 * snh0_567[k]
                   - f_2 * snh1_567[k]
                   + f_3 * pc_x[k] * sni_756[k];

        t_973[k] = f_3 * pc_y[k] * sni_756[k];

        t_974[k] = f_21 * smi_560[k]
                   + f_3 * pc_z[k] * sni_756[k];

        t_975[k] = f_16 * smi_759[k]
                   + f_4 * snh0_570[k]
                   - f_5 * snh1_570[k]
                   + f_3 * pc_x[k] * sni_759[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_x, pc_y, smi_761, smi_762, snh0_572, \
                         snh0_573, snh1_572, snh1_573, sni_758, sni_761, \
                         sni_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_3 * pc_y[k] * sni_758[k];

        t_977[k] = f_16 * smi_761[k]
                   + f_4 * snh0_572[k]
                   - f_5 * snh1_572[k]
                   + f_3 * pc_x[k] * sni_761[k];

        t_978[k] = f_16 * smi_762[k]
                   + f_6 * snh0_573[k]
                   - f_7 * snh1_573[k]
                   + f_3 * pc_x[k] * sni_762[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, pc_x, pc_y, pc_z, smi_563, smi_765, snh0_576, \
                         snh1_576, sni_759, sni_761, sni_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = f_21 * smi_563[k]
                   + f_3 * pc_z[k] * sni_759[k];

        t_980[k] = f_3 * pc_y[k] * sni_761[k];

        t_981[k] = f_16 * smi_765[k]
                   + f_6 * snh0_576[k]
                   - f_7 * snh1_576[k]
                   + f_3 * pc_x[k] * sni_765[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, pc_x, pc_z, smi_566, smi_766, smi_768, snh0_577, \
                         snh0_579, snh1_577, snh1_579, sni_762, sni_766, \
                         sni_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_16 * smi_766[k]
                   + f_8 * snh0_577[k]
                   - f_9 * snh1_577[k]
                   + f_3 * pc_x[k] * sni_766[k];

        t_983[k] = f_21 * smi_566[k]
                   + f_3 * pc_z[k] * sni_762[k];

        t_984[k] = f_16 * smi_768[k]
                   + f_8 * snh0_579[k]
                   - f_9 * snh1_579[k]
                   + f_3 * pc_x[k] * sni_768[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pc_x, pc_y, smi_770, smi_771, snh0_581, \
                         snh0_582, snh1_581, snh1_582, sni_765, sni_770, \
                         sni_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_3 * pc_y[k] * sni_765[k];

        t_986[k] = f_16 * smi_770[k]
                   + f_8 * snh0_581[k]
                   - f_9 * snh1_581[k]
                   + f_3 * pc_x[k] * sni_770[k];

        t_987[k] = f_16 * smi_771[k]
                   + f_10 * snh0_582[k]
                   - f_11 * snh1_582[k]
                   + f_3 * pc_x[k] * sni_771[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, pc_x, pc_z, smi_570, smi_773, smi_774, snh0_584, \
                         snh0_585, snh1_584, snh1_585, sni_766, sni_773, \
                         sni_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_21 * smi_570[k]
                   + f_3 * pc_z[k] * sni_766[k];

        t_989[k] = f_16 * smi_773[k]
                   + f_10 * snh0_584[k]
                   - f_11 * snh1_584[k]
                   + f_3 * pc_x[k] * sni_773[k];

        t_990[k] = f_16 * smi_774[k]
                   + f_10 * snh0_585[k]
                   - f_11 * snh1_585[k]
                   + f_3 * pc_x[k] * sni_774[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pc_x, pc_y, smi_776, smi_777, smi_778, \
                         snh0_587, snh1_587, sni_770, sni_776, sni_777, \
                         sni_778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_3 * pc_y[k] * sni_770[k];

        t_992[k] = f_16 * smi_776[k]
                   + f_10 * snh0_587[k]
                   - f_11 * snh1_587[k]
                   + f_3 * pc_x[k] * sni_776[k];

        t_993[k] = f_16 * smi_777[k]
                   + f_3 * pc_x[k] * sni_777[k];

        t_994[k] = f_16 * smi_778[k]
                   + f_3 * pc_x[k] * sni_778[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pc_x, smi_779, smi_780, smi_781, \
                         smi_782, smi_783, sni_779, sni_780, sni_781, sni_782, \
                         sni_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_16 * smi_779[k]
                   + f_3 * pc_x[k] * sni_779[k];

        t_996[k] = f_16 * smi_780[k]
                   + f_3 * pc_x[k] * sni_780[k];

        t_997[k] = f_16 * smi_781[k]
                   + f_3 * pc_x[k] * sni_781[k];

        t_998[k] = f_16 * smi_782[k]
                   + f_3 * pc_x[k] * sni_782[k];

        t_999[k] = f_16 * smi_783[k]
                   + f_3 * pc_x[k] * sni_783[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pc_y, pc_z, smi_581, snh0_582, \
                         snh0_584, snh0_585, snh1_582, snh1_584, snh1_585, sni_777, sni_779, \
                         sni_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_1 * snh0_582[k]
                    - f_2 * snh1_582[k]
                    + f_3 * pc_y[k] * sni_777[k];

        t_1001[k] = f_21 * smi_581[k]
                    + f_3 * pc_z[k] * sni_777[k];

        t_1002[k] = f_4 * snh0_584[k]
                    - f_5 * snh1_584[k]
                    + f_3 * pc_y[k] * sni_779[k];

        t_1003[k] = f_6 * snh0_585[k]
                    - f_7 * snh1_585[k]
                    + f_3 * pc_y[k] * sni_780[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, smi_587, snh0_586, \
                         snh0_587, snh1_586, snh1_587, sni_781, sni_782, \
                         sni_783 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_8 * snh0_586[k]
                    - f_9 * snh1_586[k]
                    + f_3 * pc_y[k] * sni_781[k];

        t_1005[k] = f_10 * snh0_587[k]
                    - f_11 * snh1_587[k]
                    + f_3 * pc_y[k] * sni_782[k];

        t_1006[k] = f_3 * pc_y[k] * sni_783[k];

        t_1007[k] = f_21 * smi_587[k]
                    + f_1 * snh0_587[k]
                    - f_2 * snh1_587[k]
                    + f_3 * pc_z[k] * sni_783[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smk0,
                                                          const size_t smi, const size_t smk1,
                                                          const size_t snh0, const size_t snh1,
                                                          const size_t sni, const size_t ncols,
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
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *smk0_756 = buffer.data(smk0 + 756);
    const auto *smk0_759 = buffer.data(smk0 + 759);
    const auto *smk0_762 = buffer.data(smk0 + 762);
    const auto *smk0_766 = buffer.data(smk0 + 766);
    const auto *smk0_768 = buffer.data(smk0 + 768);
    const auto *smk0_771 = buffer.data(smk0 + 771);
    const auto *smk0_773 = buffer.data(smk0 + 773);
    const auto *smk0_774 = buffer.data(smk0 + 774);
    const auto *smk0_784 = buffer.data(smk0 + 784);

    const auto *smi_588 = buffer.data(smi + 588);
    const auto *smi_590 = buffer.data(smi + 590);
    const auto *smi_591 = buffer.data(smi + 591);
    const auto *smi_593 = buffer.data(smi + 593);
    const auto *smi_594 = buffer.data(smi + 594);
    const auto *smi_595 = buffer.data(smi + 595);
    const auto *smi_597 = buffer.data(smi + 597);
    const auto *smi_598 = buffer.data(smi + 598);
    const auto *smi_599 = buffer.data(smi + 599);
    const auto *smi_600 = buffer.data(smi + 600);
    const auto *smi_602 = buffer.data(smi + 602);
    const auto *smi_609 = buffer.data(smi + 609);
    const auto *smi_611 = buffer.data(smi + 611);
    const auto *smi_612 = buffer.data(smi + 612);
    const auto *smi_613 = buffer.data(smi + 613);
    const auto *smi_614 = buffer.data(smi + 614);
    const auto *smi_615 = buffer.data(smi + 615);
    const auto *smi_616 = buffer.data(smi + 616);
    const auto *smi_618 = buffer.data(smi + 618);
    const auto *smi_619 = buffer.data(smi + 619);
    const auto *smi_621 = buffer.data(smi + 621);
    const auto *smi_622 = buffer.data(smi + 622);
    const auto *smi_625 = buffer.data(smi + 625);
    const auto *smi_626 = buffer.data(smi + 626);
    const auto *smi_630 = buffer.data(smi + 630);
    const auto *smi_637 = buffer.data(smi + 637);
    const auto *smi_639 = buffer.data(smi + 639);
    const auto *smi_640 = buffer.data(smi + 640);
    const auto *smi_641 = buffer.data(smi + 641);
    const auto *smi_642 = buffer.data(smi + 642);
    const auto *smi_643 = buffer.data(smi + 643);
    const auto *smi_644 = buffer.data(smi + 644);
    const auto *smi_646 = buffer.data(smi + 646);
    const auto *smi_649 = buffer.data(smi + 649);
    const auto *smi_653 = buffer.data(smi + 653);
    const auto *smi_658 = buffer.data(smi + 658);
    const auto *smi_665 = buffer.data(smi + 665);
    const auto *smi_667 = buffer.data(smi + 667);
    const auto *smi_668 = buffer.data(smi + 668);
    const auto *smi_669 = buffer.data(smi + 669);
    const auto *smi_670 = buffer.data(smi + 670);
    const auto *smi_671 = buffer.data(smi + 671);
    const auto *smi_784 = buffer.data(smi + 784);
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
    const auto *smi_817 = buffer.data(smi + 817);
    const auto *smi_821 = buffer.data(smi + 821);
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

    const auto *smk1_756 = buffer.data(smk1 + 756);
    const auto *smk1_759 = buffer.data(smk1 + 759);
    const auto *smk1_762 = buffer.data(smk1 + 762);
    const auto *smk1_766 = buffer.data(smk1 + 766);
    const auto *smk1_768 = buffer.data(smk1 + 768);
    const auto *smk1_771 = buffer.data(smk1 + 771);
    const auto *smk1_773 = buffer.data(smk1 + 773);
    const auto *smk1_774 = buffer.data(smk1 + 774);
    const auto *smk1_784 = buffer.data(smk1 + 784);

    const auto *snh0_588 = buffer.data(snh0 + 588);
    const auto *snh0_591 = buffer.data(snh0 + 591);
    const auto *snh0_593 = buffer.data(snh0 + 593);
    const auto *snh0_594 = buffer.data(snh0 + 594);
    const auto *snh0_597 = buffer.data(snh0 + 597);
    const auto *snh0_598 = buffer.data(snh0 + 598);
    const auto *snh0_600 = buffer.data(snh0 + 600);
    const auto *snh0_602 = buffer.data(snh0 + 602);
    const auto *snh0_603 = buffer.data(snh0 + 603);
    const auto *snh0_605 = buffer.data(snh0 + 605);
    const auto *snh0_606 = buffer.data(snh0 + 606);
    const auto *snh0_607 = buffer.data(snh0 + 607);
    const auto *snh0_608 = buffer.data(snh0 + 608);
    const auto *snh0_614 = buffer.data(snh0 + 614);
    const auto *snh0_618 = buffer.data(snh0 + 618);
    const auto *snh0_623 = buffer.data(snh0 + 623);
    const auto *snh0_626 = buffer.data(snh0 + 626);
    const auto *snh0_627 = buffer.data(snh0 + 627);
    const auto *snh0_628 = buffer.data(snh0 + 628);
    const auto *snh0_629 = buffer.data(snh0 + 629);
    const auto *snh0_630 = buffer.data(snh0 + 630);
    const auto *snh0_633 = buffer.data(snh0 + 633);
    const auto *snh0_635 = buffer.data(snh0 + 635);
    const auto *snh0_636 = buffer.data(snh0 + 636);
    const auto *snh0_639 = buffer.data(snh0 + 639);
    const auto *snh0_640 = buffer.data(snh0 + 640);
    const auto *snh0_642 = buffer.data(snh0 + 642);
    const auto *snh0_644 = buffer.data(snh0 + 644);
    const auto *snh0_645 = buffer.data(snh0 + 645);
    const auto *snh0_647 = buffer.data(snh0 + 647);
    const auto *snh0_648 = buffer.data(snh0 + 648);
    const auto *snh0_649 = buffer.data(snh0 + 649);
    const auto *snh0_650 = buffer.data(snh0 + 650);
    const auto *snh0_651 = buffer.data(snh0 + 651);

    const auto *snh1_588 = buffer.data(snh1 + 588);
    const auto *snh1_591 = buffer.data(snh1 + 591);
    const auto *snh1_593 = buffer.data(snh1 + 593);
    const auto *snh1_594 = buffer.data(snh1 + 594);
    const auto *snh1_597 = buffer.data(snh1 + 597);
    const auto *snh1_598 = buffer.data(snh1 + 598);
    const auto *snh1_600 = buffer.data(snh1 + 600);
    const auto *snh1_602 = buffer.data(snh1 + 602);
    const auto *snh1_603 = buffer.data(snh1 + 603);
    const auto *snh1_605 = buffer.data(snh1 + 605);
    const auto *snh1_606 = buffer.data(snh1 + 606);
    const auto *snh1_607 = buffer.data(snh1 + 607);
    const auto *snh1_608 = buffer.data(snh1 + 608);
    const auto *snh1_614 = buffer.data(snh1 + 614);
    const auto *snh1_618 = buffer.data(snh1 + 618);
    const auto *snh1_623 = buffer.data(snh1 + 623);
    const auto *snh1_626 = buffer.data(snh1 + 626);
    const auto *snh1_627 = buffer.data(snh1 + 627);
    const auto *snh1_628 = buffer.data(snh1 + 628);
    const auto *snh1_629 = buffer.data(snh1 + 629);
    const auto *snh1_630 = buffer.data(snh1 + 630);
    const auto *snh1_633 = buffer.data(snh1 + 633);
    const auto *snh1_635 = buffer.data(snh1 + 635);
    const auto *snh1_636 = buffer.data(snh1 + 636);
    const auto *snh1_639 = buffer.data(snh1 + 639);
    const auto *snh1_640 = buffer.data(snh1 + 640);
    const auto *snh1_642 = buffer.data(snh1 + 642);
    const auto *snh1_644 = buffer.data(snh1 + 644);
    const auto *snh1_645 = buffer.data(snh1 + 645);
    const auto *snh1_647 = buffer.data(snh1 + 647);
    const auto *snh1_648 = buffer.data(snh1 + 648);
    const auto *snh1_649 = buffer.data(snh1 + 649);
    const auto *snh1_650 = buffer.data(snh1 + 650);
    const auto *snh1_651 = buffer.data(snh1 + 651);

    const auto *sni_784 = buffer.data(sni + 784);
    const auto *sni_786 = buffer.data(sni + 786);
    const auto *sni_787 = buffer.data(sni + 787);
    const auto *sni_789 = buffer.data(sni + 789);
    const auto *sni_790 = buffer.data(sni + 790);
    const auto *sni_793 = buffer.data(sni + 793);
    const auto *sni_794 = buffer.data(sni + 794);
    const auto *sni_796 = buffer.data(sni + 796);
    const auto *sni_798 = buffer.data(sni + 798);
    const auto *sni_799 = buffer.data(sni + 799);
    const auto *sni_801 = buffer.data(sni + 801);
    const auto *sni_802 = buffer.data(sni + 802);
    const auto *sni_804 = buffer.data(sni + 804);
    const auto *sni_805 = buffer.data(sni + 805);
    const auto *sni_806 = buffer.data(sni + 806);
    const auto *sni_807 = buffer.data(sni + 807);
    const auto *sni_808 = buffer.data(sni + 808);
    const auto *sni_809 = buffer.data(sni + 809);
    const auto *sni_810 = buffer.data(sni + 810);
    const auto *sni_811 = buffer.data(sni + 811);
    const auto *sni_812 = buffer.data(sni + 812);
    const auto *sni_814 = buffer.data(sni + 814);
    const auto *sni_815 = buffer.data(sni + 815);
    const auto *sni_817 = buffer.data(sni + 817);
    const auto *sni_818 = buffer.data(sni + 818);
    const auto *sni_821 = buffer.data(sni + 821);
    const auto *sni_822 = buffer.data(sni + 822);
    const auto *sni_826 = buffer.data(sni + 826);
    const auto *sni_832 = buffer.data(sni + 832);
    const auto *sni_833 = buffer.data(sni + 833);
    const auto *sni_834 = buffer.data(sni + 834);
    const auto *sni_835 = buffer.data(sni + 835);
    const auto *sni_836 = buffer.data(sni + 836);
    const auto *sni_837 = buffer.data(sni + 837);
    const auto *sni_838 = buffer.data(sni + 838);
    const auto *sni_839 = buffer.data(sni + 839);
    const auto *sni_840 = buffer.data(sni + 840);
    const auto *sni_842 = buffer.data(sni + 842);
    const auto *sni_843 = buffer.data(sni + 843);
    const auto *sni_845 = buffer.data(sni + 845);
    const auto *sni_846 = buffer.data(sni + 846);
    const auto *sni_849 = buffer.data(sni + 849);
    const auto *sni_850 = buffer.data(sni + 850);
    const auto *sni_852 = buffer.data(sni + 852);
    const auto *sni_854 = buffer.data(sni + 854);
    const auto *sni_855 = buffer.data(sni + 855);
    const auto *sni_857 = buffer.data(sni + 857);
    const auto *sni_858 = buffer.data(sni + 858);
    const auto *sni_860 = buffer.data(sni + 860);
    const auto *sni_861 = buffer.data(sni + 861);
    const auto *sni_862 = buffer.data(sni + 862);
    const auto *sni_863 = buffer.data(sni + 863);
    const auto *sni_864 = buffer.data(sni + 864);
    const auto *sni_865 = buffer.data(sni + 865);
    const auto *sni_866 = buffer.data(sni + 866);
    const auto *sni_867 = buffer.data(sni + 867);
    const auto *sni_868 = buffer.data(sni + 868);

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, smi_588, smi_784, \
                         smi_787, snh0_588, snh0_591, snh1_588, snh1_591, sni_784, \
                         sni_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_15 * smi_784[k]
                    + f_1 * snh0_588[k]
                    - f_2 * snh1_588[k]
                    + f_3 * pc_x[k] * sni_784[k];

        t_1009[k] = f_20 * smi_588[k]
                    + f_3 * pc_y[k] * sni_784[k];

        t_1010[k] = f_3 * pc_z[k] * sni_784[k];

        t_1011[k] = f_15 * smi_787[k]
                    + f_4 * snh0_591[k]
                    - f_5 * snh1_591[k]
                    + f_3 * pc_x[k] * sni_787[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pc_x, pc_y, smi_590, smi_789, smi_790, \
                         snh0_593, snh0_594, snh1_593, snh1_594, sni_786, sni_789, \
                         sni_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_20 * smi_590[k]
                    + f_3 * pc_y[k] * sni_786[k];

        t_1013[k] = f_15 * smi_789[k]
                    + f_4 * snh0_593[k]
                    - f_5 * snh1_593[k]
                    + f_3 * pc_x[k] * sni_789[k];

        t_1014[k] = f_15 * smi_790[k]
                    + f_6 * snh0_594[k]
                    - f_7 * snh1_594[k]
                    + f_3 * pc_x[k] * sni_790[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, pc_x, pc_y, pc_z, smi_593, smi_793, snh0_597, \
                         snh1_597, sni_787, sni_789, sni_793 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = f_3 * pc_z[k] * sni_787[k];

        t_1016[k] = f_20 * smi_593[k]
                    + f_3 * pc_y[k] * sni_789[k];

        t_1017[k] = f_15 * smi_793[k]
                    + f_6 * snh0_597[k]
                    - f_7 * snh1_597[k]
                    + f_3 * pc_x[k] * sni_793[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, pc_x, pc_z, smi_794, smi_796, snh0_598, \
                         snh0_600, snh1_598, snh1_600, sni_790, sni_794, \
                         sni_796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_15 * smi_794[k]
                    + f_8 * snh0_598[k]
                    - f_9 * snh1_598[k]
                    + f_3 * pc_x[k] * sni_794[k];

        t_1019[k] = f_3 * pc_z[k] * sni_790[k];

        t_1020[k] = f_15 * smi_796[k]
                    + f_8 * snh0_600[k]
                    - f_9 * snh1_600[k]
                    + f_3 * pc_x[k] * sni_796[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, pc_x, pc_y, smi_597, smi_798, smi_799, \
                         snh0_602, snh0_603, snh1_602, snh1_603, sni_793, sni_798, \
                         sni_799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_20 * smi_597[k]
                    + f_3 * pc_y[k] * sni_793[k];

        t_1022[k] = f_15 * smi_798[k]
                    + f_8 * snh0_602[k]
                    - f_9 * snh1_602[k]
                    + f_3 * pc_x[k] * sni_798[k];

        t_1023[k] = f_15 * smi_799[k]
                    + f_10 * snh0_603[k]
                    - f_11 * snh1_603[k]
                    + f_3 * pc_x[k] * sni_799[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, pc_x, pc_z, smi_801, smi_802, snh0_605, \
                         snh0_606, snh1_605, snh1_606, sni_794, sni_801, \
                         sni_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_3 * pc_z[k] * sni_794[k];

        t_1025[k] = f_15 * smi_801[k]
                    + f_10 * snh0_605[k]
                    - f_11 * snh1_605[k]
                    + f_3 * pc_x[k] * sni_801[k];

        t_1026[k] = f_15 * smi_802[k]
                    + f_10 * snh0_606[k]
                    - f_11 * snh1_606[k]
                    + f_3 * pc_x[k] * sni_802[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, pc_x, pc_y, smi_602, smi_804, \
                         smi_805, smi_806, snh0_608, snh1_608, sni_798, sni_804, sni_805, \
                         sni_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_20 * smi_602[k]
                    + f_3 * pc_y[k] * sni_798[k];

        t_1028[k] = f_15 * smi_804[k]
                    + f_10 * snh0_608[k]
                    - f_11 * snh1_608[k]
                    + f_3 * pc_x[k] * sni_804[k];

        t_1029[k] = f_15 * smi_805[k]
                    + f_3 * pc_x[k] * sni_805[k];

        t_1030[k] = f_15 * smi_806[k]
                    + f_3 * pc_x[k] * sni_806[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, t_1034, t_1035, pc_x, smi_807, smi_808, \
                         smi_809, smi_810, smi_811, sni_807, sni_808, sni_809, sni_810, \
                         sni_811 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_15 * smi_807[k]
                    + f_3 * pc_x[k] * sni_807[k];

        t_1032[k] = f_15 * smi_808[k]
                    + f_3 * pc_x[k] * sni_808[k];

        t_1033[k] = f_15 * smi_809[k]
                    + f_3 * pc_x[k] * sni_809[k];

        t_1034[k] = f_15 * smi_810[k]
                    + f_3 * pc_x[k] * sni_810[k];

        t_1035[k] = f_15 * smi_811[k]
                    + f_3 * pc_x[k] * sni_811[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, pc_y, pc_z, smi_609, smi_611, snh0_603, \
                         snh0_605, snh1_603, snh1_605, sni_805, \
                         sni_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = f_20 * smi_609[k]
                    + f_1 * snh0_603[k]
                    - f_2 * snh1_603[k]
                    + f_3 * pc_y[k] * sni_805[k];

        t_1037[k] = f_3 * pc_z[k] * sni_805[k];

        t_1038[k] = f_20 * smi_611[k]
                    + f_4 * snh0_605[k]
                    - f_5 * snh1_605[k]
                    + f_3 * pc_y[k] * sni_807[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pc_y, smi_612, smi_613, smi_614, snh0_606, \
                         snh0_607, snh0_608, snh1_606, snh1_607, snh1_608, sni_808, sni_809, \
                         sni_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_20 * smi_612[k]
                    + f_6 * snh0_606[k]
                    - f_7 * snh1_606[k]
                    + f_3 * pc_y[k] * sni_808[k];

        t_1040[k] = f_20 * smi_613[k]
                    + f_8 * snh0_607[k]
                    - f_9 * snh1_607[k]
                    + f_3 * pc_y[k] * sni_809[k];

        t_1041[k] = f_20 * smi_614[k]
                    + f_10 * snh0_608[k]
                    - f_11 * snh1_608[k]
                    + f_3 * pc_y[k] * sni_810[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, t_1045, pb_z, pc_y, pc_z, smk0_756, smi_615, \
                         smi_616, smk1_756, snh0_608, snh1_608, sni_811, \
                         sni_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_20 * smi_615[k]
                    + f_3 * pc_y[k] * sni_811[k];

        t_1043[k] = f_1 * snh0_608[k]
                    - f_2 * snh1_608[k]
                    + f_3 * pc_z[k] * sni_811[k];

        t_1044[k] = pb_z[k] * smk0_756[k]
                    - f_12 * pc_z[k] * smk1_756[k];

        t_1045[k] = f_21 * smi_616[k]
                    + f_3 * pc_y[k] * sni_812[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pb_z, pc_y, pc_z, smk0_759, smi_588, smi_618, \
                         smk1_759, sni_812, sni_814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_13 * smi_588[k]
                    + f_3 * pc_z[k] * sni_812[k];

        t_1047[k] = pb_z[k] * smk0_759[k]
                    - f_12 * pc_z[k] * smk1_759[k];

        t_1048[k] = f_21 * smi_618[k]
                    + f_3 * pc_y[k] * sni_814[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pb_z, pc_x, pc_z, smk0_762, smi_591, smi_817, \
                         smk1_762, snh0_614, snh1_614, sni_815, \
                         sni_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_15 * smi_817[k]
                    + f_4 * snh0_614[k]
                    - f_5 * snh1_614[k]
                    + f_3 * pc_x[k] * sni_817[k];

        t_1050[k] = pb_z[k] * smk0_762[k]
                    - f_12 * pc_z[k] * smk1_762[k];

        t_1051[k] = f_13 * smi_591[k]
                    + f_3 * pc_z[k] * sni_815[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pb_z, pc_x, pc_y, pc_z, smk0_766, smi_621, \
                         smi_821, smk1_766, snh0_618, snh1_618, sni_817, \
                         sni_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_21 * smi_621[k]
                    + f_3 * pc_y[k] * sni_817[k];

        t_1053[k] = f_15 * smi_821[k]
                    + f_6 * snh0_618[k]
                    - f_7 * snh1_618[k]
                    + f_3 * pc_x[k] * sni_821[k];

        t_1054[k] = pb_z[k] * smk0_766[k]
                    - f_12 * pc_z[k] * smk1_766[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pb_z, pc_y, pc_z, smk0_768, smi_594, smi_595, \
                         smi_625, smk1_768, sni_818, sni_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_13 * smi_594[k]
                    + f_3 * pc_z[k] * sni_818[k];

        t_1056[k] = pb_z[k] * smk0_768[k]
                    + f_14 * smi_595[k]
                    - f_12 * pc_z[k] * smk1_768[k];

        t_1057[k] = f_21 * smi_625[k]
                    + f_3 * pc_y[k] * sni_821[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pb_z, pc_x, pc_z, smk0_771, smi_598, smi_826, \
                         smk1_771, snh0_623, snh1_623, sni_822, \
                         sni_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_15 * smi_826[k]
                    + f_8 * snh0_623[k]
                    - f_9 * snh1_623[k]
                    + f_3 * pc_x[k] * sni_826[k];

        t_1059[k] = pb_z[k] * smk0_771[k]
                    - f_12 * pc_z[k] * smk1_771[k];

        t_1060[k] = f_13 * smi_598[k]
                    + f_3 * pc_z[k] * sni_822[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pb_z, pc_y, pc_z, smk0_773, smk0_774, \
                         smi_599, smi_600, smi_630, smk1_773, smk1_774, \
                         sni_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = pb_z[k] * smk0_773[k]
                    + f_14 * smi_599[k]
                    - f_12 * pc_z[k] * smk1_773[k];

        t_1062[k] = pb_z[k] * smk0_774[k]
                    + f_15 * smi_600[k]
                    - f_12 * pc_z[k] * smk1_774[k];

        t_1063[k] = f_21 * smi_630[k]
                    + f_3 * pc_y[k] * sni_826[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, t_1067, pc_x, smi_832, smi_833, smi_834, \
                         smi_835, snh0_629, snh1_629, sni_832, sni_833, sni_834, \
                         sni_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_15 * smi_832[k]
                    + f_10 * snh0_629[k]
                    - f_11 * snh1_629[k]
                    + f_3 * pc_x[k] * sni_832[k];

        t_1065[k] = f_15 * smi_833[k]
                    + f_3 * pc_x[k] * sni_833[k];

        t_1066[k] = f_15 * smi_834[k]
                    + f_3 * pc_x[k] * sni_834[k];

        t_1067[k] = f_15 * smi_835[k]
                    + f_3 * pc_x[k] * sni_835[k];
    }

#pragma omp simd aligned(t_1068, t_1069, t_1070, t_1071, pc_x, smi_836, smi_837, smi_838, \
                         smi_839, sni_836, sni_837, sni_838, sni_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = f_15 * smi_836[k]
                    + f_3 * pc_x[k] * sni_836[k];

        t_1069[k] = f_15 * smi_837[k]
                    + f_3 * pc_x[k] * sni_837[k];

        t_1070[k] = f_15 * smi_838[k]
                    + f_3 * pc_x[k] * sni_838[k];

        t_1071[k] = f_15 * smi_839[k]
                    + f_3 * pc_x[k] * sni_839[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, pb_z, pc_y, pc_z, smk0_784, smi_609, smi_639, \
                         smk1_784, snh0_626, snh1_626, sni_833, \
                         sni_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = pb_z[k] * smk0_784[k]
                    - f_12 * pc_z[k] * smk1_784[k];

        t_1073[k] = f_13 * smi_609[k]
                    + f_3 * pc_z[k] * sni_833[k];

        t_1074[k] = f_21 * smi_639[k]
                    + f_4 * snh0_626[k]
                    - f_5 * snh1_626[k]
                    + f_3 * pc_y[k] * sni_835[k];
    }

#pragma omp simd aligned(t_1075, t_1076, t_1077, pc_y, smi_640, smi_641, smi_642, snh0_627, \
                         snh0_628, snh0_629, snh1_627, snh1_628, snh1_629, sni_836, sni_837, \
                         sni_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1075[k] = f_21 * smi_640[k]
                    + f_6 * snh0_627[k]
                    - f_7 * snh1_627[k]
                    + f_3 * pc_y[k] * sni_836[k];

        t_1076[k] = f_21 * smi_641[k]
                    + f_8 * snh0_628[k]
                    - f_9 * snh1_628[k]
                    + f_3 * pc_y[k] * sni_837[k];

        t_1077[k] = f_21 * smi_642[k]
                    + f_10 * snh0_629[k]
                    - f_11 * snh1_629[k]
                    + f_3 * pc_y[k] * sni_838[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, pc_x, pc_y, pc_z, smi_615, smi_643, smi_840, \
                         snh0_629, snh0_630, snh1_629, snh1_630, sni_839, \
                         sni_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = f_21 * smi_643[k]
                    + f_3 * pc_y[k] * sni_839[k];

        t_1079[k] = f_13 * smi_615[k]
                    + f_1 * snh0_629[k]
                    - f_2 * snh1_629[k]
                    + f_3 * pc_z[k] * sni_839[k];

        t_1080[k] = f_15 * smi_840[k]
                    + f_1 * snh0_630[k]
                    - f_2 * snh1_630[k]
                    + f_3 * pc_x[k] * sni_840[k];
    }

#pragma omp simd aligned(t_1081, t_1082, t_1083, t_1084, pc_x, pc_y, pc_z, smi_616, smi_644, \
                         smi_646, smi_843, snh0_633, snh1_633, sni_840, sni_842, \
                         sni_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1081[k] = f_17 * smi_644[k]
                    + f_3 * pc_y[k] * sni_840[k];

        t_1082[k] = f_14 * smi_616[k]
                    + f_3 * pc_z[k] * sni_840[k];

        t_1083[k] = f_15 * smi_843[k]
                    + f_4 * snh0_633[k]
                    - f_5 * snh1_633[k]
                    + f_3 * pc_x[k] * sni_843[k];

        t_1084[k] = f_17 * smi_646[k]
                    + f_3 * pc_y[k] * sni_842[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_x, pc_z, smi_619, smi_845, smi_846, \
                         snh0_635, snh0_636, snh1_635, snh1_636, sni_843, sni_845, \
                         sni_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_15 * smi_845[k]
                    + f_4 * snh0_635[k]
                    - f_5 * snh1_635[k]
                    + f_3 * pc_x[k] * sni_845[k];

        t_1086[k] = f_15 * smi_846[k]
                    + f_6 * snh0_636[k]
                    - f_7 * snh1_636[k]
                    + f_3 * pc_x[k] * sni_846[k];

        t_1087[k] = f_14 * smi_619[k]
                    + f_3 * pc_z[k] * sni_843[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_x, pc_y, smi_649, smi_849, smi_850, \
                         snh0_639, snh0_640, snh1_639, snh1_640, sni_845, sni_849, \
                         sni_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * smi_649[k]
                    + f_3 * pc_y[k] * sni_845[k];

        t_1089[k] = f_15 * smi_849[k]
                    + f_6 * snh0_639[k]
                    - f_7 * snh1_639[k]
                    + f_3 * pc_x[k] * sni_849[k];

        t_1090[k] = f_15 * smi_850[k]
                    + f_8 * snh0_640[k]
                    - f_9 * snh1_640[k]
                    + f_3 * pc_x[k] * sni_850[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_y, pc_z, smi_622, smi_653, smi_852, \
                         snh0_642, snh1_642, sni_846, sni_849, \
                         sni_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_14 * smi_622[k]
                    + f_3 * pc_z[k] * sni_846[k];

        t_1092[k] = f_15 * smi_852[k]
                    + f_8 * snh0_642[k]
                    - f_9 * snh1_642[k]
                    + f_3 * pc_x[k] * sni_852[k];

        t_1093[k] = f_17 * smi_653[k]
                    + f_3 * pc_y[k] * sni_849[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, pc_z, smi_626, smi_854, smi_855, \
                         snh0_644, snh0_645, snh1_644, snh1_645, sni_850, sni_854, \
                         sni_855 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_15 * smi_854[k]
                    + f_8 * snh0_644[k]
                    - f_9 * snh1_644[k]
                    + f_3 * pc_x[k] * sni_854[k];

        t_1095[k] = f_15 * smi_855[k]
                    + f_10 * snh0_645[k]
                    - f_11 * snh1_645[k]
                    + f_3 * pc_x[k] * sni_855[k];

        t_1096[k] = f_14 * smi_626[k]
                    + f_3 * pc_z[k] * sni_850[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_y, smi_658, smi_857, smi_858, \
                         snh0_647, snh0_648, snh1_647, snh1_648, sni_854, sni_857, \
                         sni_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_15 * smi_857[k]
                    + f_10 * snh0_647[k]
                    - f_11 * snh1_647[k]
                    + f_3 * pc_x[k] * sni_857[k];

        t_1098[k] = f_15 * smi_858[k]
                    + f_10 * snh0_648[k]
                    - f_11 * snh1_648[k]
                    + f_3 * pc_x[k] * sni_858[k];

        t_1099[k] = f_17 * smi_658[k]
                    + f_3 * pc_y[k] * sni_854[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, t_1103, pc_x, smi_860, smi_861, smi_862, \
                         smi_863, snh0_650, snh1_650, sni_860, sni_861, sni_862, \
                         sni_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_15 * smi_860[k]
                    + f_10 * snh0_650[k]
                    - f_11 * snh1_650[k]
                    + f_3 * pc_x[k] * sni_860[k];

        t_1101[k] = f_15 * smi_861[k]
                    + f_3 * pc_x[k] * sni_861[k];

        t_1102[k] = f_15 * smi_862[k]
                    + f_3 * pc_x[k] * sni_862[k];

        t_1103[k] = f_15 * smi_863[k]
                    + f_3 * pc_x[k] * sni_863[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, t_1107, pc_x, smi_864, smi_865, smi_866, \
                         smi_867, sni_864, sni_865, sni_866, sni_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_15 * smi_864[k]
                    + f_3 * pc_x[k] * sni_864[k];

        t_1105[k] = f_15 * smi_865[k]
                    + f_3 * pc_x[k] * sni_865[k];

        t_1106[k] = f_15 * smi_866[k]
                    + f_3 * pc_x[k] * sni_866[k];

        t_1107[k] = f_15 * smi_867[k]
                    + f_3 * pc_x[k] * sni_867[k];
    }

#pragma omp simd aligned(t_1108, t_1109, t_1110, pc_y, pc_z, smi_637, smi_665, smi_667, \
                         snh0_645, snh0_647, snh1_645, snh1_647, sni_861, \
                         sni_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1108[k] = f_17 * smi_665[k]
                    + f_1 * snh0_645[k]
                    - f_2 * snh1_645[k]
                    + f_3 * pc_y[k] * sni_861[k];

        t_1109[k] = f_14 * smi_637[k]
                    + f_3 * pc_z[k] * sni_861[k];

        t_1110[k] = f_17 * smi_667[k]
                    + f_4 * snh0_647[k]
                    - f_5 * snh1_647[k]
                    + f_3 * pc_y[k] * sni_863[k];
    }

#pragma omp simd aligned(t_1111, t_1112, t_1113, pc_y, smi_668, smi_669, smi_670, snh0_648, \
                         snh0_649, snh0_650, snh1_648, snh1_649, snh1_650, sni_864, sni_865, \
                         sni_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1111[k] = f_17 * smi_668[k]
                    + f_6 * snh0_648[k]
                    - f_7 * snh1_648[k]
                    + f_3 * pc_y[k] * sni_864[k];

        t_1112[k] = f_17 * smi_669[k]
                    + f_8 * snh0_649[k]
                    - f_9 * snh1_649[k]
                    + f_3 * pc_y[k] * sni_865[k];

        t_1113[k] = f_17 * smi_670[k]
                    + f_10 * snh0_650[k]
                    - f_11 * snh1_650[k]
                    + f_3 * pc_y[k] * sni_866[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_x, pc_y, pc_z, smi_643, smi_671, smi_868, \
                         snh0_650, snh0_651, snh1_650, snh1_651, sni_867, \
                         sni_868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_17 * smi_671[k]
                    + f_3 * pc_y[k] * sni_867[k];

        t_1115[k] = f_14 * smi_643[k]
                    + f_1 * snh0_650[k]
                    - f_2 * snh1_650[k]
                    + f_3 * pc_z[k] * sni_867[k];

        t_1116[k] = f_15 * smi_868[k]
                    + f_1 * snh0_651[k]
                    - f_2 * snh1_651[k]
                    + f_3 * pc_x[k] * sni_868[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t smi, const size_t snh0,
                                                           const size_t snh1, const size_t sni,
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

    const auto *smi_644 = buffer.data(smi + 644);
    const auto *smi_647 = buffer.data(smi + 647);
    const auto *smi_650 = buffer.data(smi + 650);
    const auto *smi_654 = buffer.data(smi + 654);
    const auto *smi_665 = buffer.data(smi + 665);
    const auto *smi_671 = buffer.data(smi + 671);
    const auto *smi_672 = buffer.data(smi + 672);
    const auto *smi_674 = buffer.data(smi + 674);
    const auto *smi_675 = buffer.data(smi + 675);
    const auto *smi_677 = buffer.data(smi + 677);
    const auto *smi_678 = buffer.data(smi + 678);
    const auto *smi_681 = buffer.data(smi + 681);
    const auto *smi_682 = buffer.data(smi + 682);
    const auto *smi_686 = buffer.data(smi + 686);
    const auto *smi_693 = buffer.data(smi + 693);
    const auto *smi_695 = buffer.data(smi + 695);
    const auto *smi_696 = buffer.data(smi + 696);
    const auto *smi_697 = buffer.data(smi + 697);
    const auto *smi_698 = buffer.data(smi + 698);
    const auto *smi_699 = buffer.data(smi + 699);
    const auto *smi_700 = buffer.data(smi + 700);
    const auto *smi_702 = buffer.data(smi + 702);
    const auto *smi_703 = buffer.data(smi + 703);
    const auto *smi_705 = buffer.data(smi + 705);
    const auto *smi_706 = buffer.data(smi + 706);
    const auto *smi_709 = buffer.data(smi + 709);
    const auto *smi_710 = buffer.data(smi + 710);
    const auto *smi_714 = buffer.data(smi + 714);
    const auto *smi_721 = buffer.data(smi + 721);
    const auto *smi_723 = buffer.data(smi + 723);
    const auto *smi_724 = buffer.data(smi + 724);
    const auto *smi_725 = buffer.data(smi + 725);
    const auto *smi_726 = buffer.data(smi + 726);
    const auto *smi_727 = buffer.data(smi + 727);
    const auto *smi_728 = buffer.data(smi + 728);
    const auto *smi_730 = buffer.data(smi + 730);
    const auto *smi_733 = buffer.data(smi + 733);
    const auto *smi_737 = buffer.data(smi + 737);
    const auto *smi_742 = buffer.data(smi + 742);
    const auto *smi_749 = buffer.data(smi + 749);
    const auto *smi_751 = buffer.data(smi + 751);
    const auto *smi_752 = buffer.data(smi + 752);
    const auto *smi_753 = buffer.data(smi + 753);
    const auto *smi_754 = buffer.data(smi + 754);
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

    const auto *snh0_654 = buffer.data(snh0 + 654);
    const auto *snh0_656 = buffer.data(snh0 + 656);
    const auto *snh0_657 = buffer.data(snh0 + 657);
    const auto *snh0_660 = buffer.data(snh0 + 660);
    const auto *snh0_661 = buffer.data(snh0 + 661);
    const auto *snh0_663 = buffer.data(snh0 + 663);
    const auto *snh0_665 = buffer.data(snh0 + 665);
    const auto *snh0_666 = buffer.data(snh0 + 666);
    const auto *snh0_668 = buffer.data(snh0 + 668);
    const auto *snh0_669 = buffer.data(snh0 + 669);
    const auto *snh0_670 = buffer.data(snh0 + 670);
    const auto *snh0_671 = buffer.data(snh0 + 671);
    const auto *snh0_672 = buffer.data(snh0 + 672);
    const auto *snh0_675 = buffer.data(snh0 + 675);
    const auto *snh0_677 = buffer.data(snh0 + 677);
    const auto *snh0_678 = buffer.data(snh0 + 678);
    const auto *snh0_681 = buffer.data(snh0 + 681);
    const auto *snh0_682 = buffer.data(snh0 + 682);
    const auto *snh0_684 = buffer.data(snh0 + 684);
    const auto *snh0_686 = buffer.data(snh0 + 686);
    const auto *snh0_687 = buffer.data(snh0 + 687);
    const auto *snh0_689 = buffer.data(snh0 + 689);
    const auto *snh0_690 = buffer.data(snh0 + 690);
    const auto *snh0_691 = buffer.data(snh0 + 691);
    const auto *snh0_692 = buffer.data(snh0 + 692);
    const auto *snh0_693 = buffer.data(snh0 + 693);
    const auto *snh0_696 = buffer.data(snh0 + 696);
    const auto *snh0_698 = buffer.data(snh0 + 698);
    const auto *snh0_699 = buffer.data(snh0 + 699);
    const auto *snh0_702 = buffer.data(snh0 + 702);
    const auto *snh0_703 = buffer.data(snh0 + 703);
    const auto *snh0_705 = buffer.data(snh0 + 705);
    const auto *snh0_707 = buffer.data(snh0 + 707);
    const auto *snh0_708 = buffer.data(snh0 + 708);
    const auto *snh0_710 = buffer.data(snh0 + 710);
    const auto *snh0_711 = buffer.data(snh0 + 711);
    const auto *snh0_712 = buffer.data(snh0 + 712);
    const auto *snh0_713 = buffer.data(snh0 + 713);

    const auto *snh1_654 = buffer.data(snh1 + 654);
    const auto *snh1_656 = buffer.data(snh1 + 656);
    const auto *snh1_657 = buffer.data(snh1 + 657);
    const auto *snh1_660 = buffer.data(snh1 + 660);
    const auto *snh1_661 = buffer.data(snh1 + 661);
    const auto *snh1_663 = buffer.data(snh1 + 663);
    const auto *snh1_665 = buffer.data(snh1 + 665);
    const auto *snh1_666 = buffer.data(snh1 + 666);
    const auto *snh1_668 = buffer.data(snh1 + 668);
    const auto *snh1_669 = buffer.data(snh1 + 669);
    const auto *snh1_670 = buffer.data(snh1 + 670);
    const auto *snh1_671 = buffer.data(snh1 + 671);
    const auto *snh1_672 = buffer.data(snh1 + 672);
    const auto *snh1_675 = buffer.data(snh1 + 675);
    const auto *snh1_677 = buffer.data(snh1 + 677);
    const auto *snh1_678 = buffer.data(snh1 + 678);
    const auto *snh1_681 = buffer.data(snh1 + 681);
    const auto *snh1_682 = buffer.data(snh1 + 682);
    const auto *snh1_684 = buffer.data(snh1 + 684);
    const auto *snh1_686 = buffer.data(snh1 + 686);
    const auto *snh1_687 = buffer.data(snh1 + 687);
    const auto *snh1_689 = buffer.data(snh1 + 689);
    const auto *snh1_690 = buffer.data(snh1 + 690);
    const auto *snh1_691 = buffer.data(snh1 + 691);
    const auto *snh1_692 = buffer.data(snh1 + 692);
    const auto *snh1_693 = buffer.data(snh1 + 693);
    const auto *snh1_696 = buffer.data(snh1 + 696);
    const auto *snh1_698 = buffer.data(snh1 + 698);
    const auto *snh1_699 = buffer.data(snh1 + 699);
    const auto *snh1_702 = buffer.data(snh1 + 702);
    const auto *snh1_703 = buffer.data(snh1 + 703);
    const auto *snh1_705 = buffer.data(snh1 + 705);
    const auto *snh1_707 = buffer.data(snh1 + 707);
    const auto *snh1_708 = buffer.data(snh1 + 708);
    const auto *snh1_710 = buffer.data(snh1 + 710);
    const auto *snh1_711 = buffer.data(snh1 + 711);
    const auto *snh1_712 = buffer.data(snh1 + 712);
    const auto *snh1_713 = buffer.data(snh1 + 713);

    const auto *sni_868 = buffer.data(sni + 868);
    const auto *sni_870 = buffer.data(sni + 870);
    const auto *sni_871 = buffer.data(sni + 871);
    const auto *sni_873 = buffer.data(sni + 873);
    const auto *sni_874 = buffer.data(sni + 874);
    const auto *sni_877 = buffer.data(sni + 877);
    const auto *sni_878 = buffer.data(sni + 878);
    const auto *sni_880 = buffer.data(sni + 880);
    const auto *sni_882 = buffer.data(sni + 882);
    const auto *sni_883 = buffer.data(sni + 883);
    const auto *sni_885 = buffer.data(sni + 885);
    const auto *sni_886 = buffer.data(sni + 886);
    const auto *sni_888 = buffer.data(sni + 888);
    const auto *sni_889 = buffer.data(sni + 889);
    const auto *sni_890 = buffer.data(sni + 890);
    const auto *sni_891 = buffer.data(sni + 891);
    const auto *sni_892 = buffer.data(sni + 892);
    const auto *sni_893 = buffer.data(sni + 893);
    const auto *sni_894 = buffer.data(sni + 894);
    const auto *sni_895 = buffer.data(sni + 895);
    const auto *sni_896 = buffer.data(sni + 896);
    const auto *sni_898 = buffer.data(sni + 898);
    const auto *sni_899 = buffer.data(sni + 899);
    const auto *sni_901 = buffer.data(sni + 901);
    const auto *sni_902 = buffer.data(sni + 902);
    const auto *sni_905 = buffer.data(sni + 905);
    const auto *sni_906 = buffer.data(sni + 906);
    const auto *sni_908 = buffer.data(sni + 908);
    const auto *sni_910 = buffer.data(sni + 910);
    const auto *sni_911 = buffer.data(sni + 911);
    const auto *sni_913 = buffer.data(sni + 913);
    const auto *sni_914 = buffer.data(sni + 914);
    const auto *sni_916 = buffer.data(sni + 916);
    const auto *sni_917 = buffer.data(sni + 917);
    const auto *sni_918 = buffer.data(sni + 918);
    const auto *sni_919 = buffer.data(sni + 919);
    const auto *sni_920 = buffer.data(sni + 920);
    const auto *sni_921 = buffer.data(sni + 921);
    const auto *sni_922 = buffer.data(sni + 922);
    const auto *sni_923 = buffer.data(sni + 923);
    const auto *sni_924 = buffer.data(sni + 924);
    const auto *sni_926 = buffer.data(sni + 926);
    const auto *sni_927 = buffer.data(sni + 927);
    const auto *sni_929 = buffer.data(sni + 929);
    const auto *sni_930 = buffer.data(sni + 930);
    const auto *sni_933 = buffer.data(sni + 933);
    const auto *sni_934 = buffer.data(sni + 934);
    const auto *sni_936 = buffer.data(sni + 936);
    const auto *sni_938 = buffer.data(sni + 938);
    const auto *sni_939 = buffer.data(sni + 939);
    const auto *sni_941 = buffer.data(sni + 941);
    const auto *sni_942 = buffer.data(sni + 942);
    const auto *sni_944 = buffer.data(sni + 944);
    const auto *sni_945 = buffer.data(sni + 945);
    const auto *sni_946 = buffer.data(sni + 946);
    const auto *sni_947 = buffer.data(sni + 947);
    const auto *sni_948 = buffer.data(sni + 948);
    const auto *sni_949 = buffer.data(sni + 949);
    const auto *sni_950 = buffer.data(sni + 950);
    const auto *sni_951 = buffer.data(sni + 951);

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, pc_x, pc_y, pc_z, smi_644, smi_672, \
                         smi_674, smi_871, snh0_654, snh1_654, sni_868, sni_870, \
                         sni_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_16 * smi_672[k]
                    + f_3 * pc_y[k] * sni_868[k];

        t_1118[k] = f_15 * smi_644[k]
                    + f_3 * pc_z[k] * sni_868[k];

        t_1119[k] = f_15 * smi_871[k]
                    + f_4 * snh0_654[k]
                    - f_5 * snh1_654[k]
                    + f_3 * pc_x[k] * sni_871[k];

        t_1120[k] = f_16 * smi_674[k]
                    + f_3 * pc_y[k] * sni_870[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_x, pc_z, smi_647, smi_873, smi_874, \
                         snh0_656, snh0_657, snh1_656, snh1_657, sni_871, sni_873, \
                         sni_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_15 * smi_873[k]
                    + f_4 * snh0_656[k]
                    - f_5 * snh1_656[k]
                    + f_3 * pc_x[k] * sni_873[k];

        t_1122[k] = f_15 * smi_874[k]
                    + f_6 * snh0_657[k]
                    - f_7 * snh1_657[k]
                    + f_3 * pc_x[k] * sni_874[k];

        t_1123[k] = f_15 * smi_647[k]
                    + f_3 * pc_z[k] * sni_871[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_y, smi_677, smi_877, smi_878, \
                         snh0_660, snh0_661, snh1_660, snh1_661, sni_873, sni_877, \
                         sni_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_16 * smi_677[k]
                    + f_3 * pc_y[k] * sni_873[k];

        t_1125[k] = f_15 * smi_877[k]
                    + f_6 * snh0_660[k]
                    - f_7 * snh1_660[k]
                    + f_3 * pc_x[k] * sni_877[k];

        t_1126[k] = f_15 * smi_878[k]
                    + f_8 * snh0_661[k]
                    - f_9 * snh1_661[k]
                    + f_3 * pc_x[k] * sni_878[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, pc_y, pc_z, smi_650, smi_681, smi_880, \
                         snh0_663, snh1_663, sni_874, sni_877, \
                         sni_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_15 * smi_650[k]
                    + f_3 * pc_z[k] * sni_874[k];

        t_1128[k] = f_15 * smi_880[k]
                    + f_8 * snh0_663[k]
                    - f_9 * snh1_663[k]
                    + f_3 * pc_x[k] * sni_880[k];

        t_1129[k] = f_16 * smi_681[k]
                    + f_3 * pc_y[k] * sni_877[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, pc_z, smi_654, smi_882, smi_883, \
                         snh0_665, snh0_666, snh1_665, snh1_666, sni_878, sni_882, \
                         sni_883 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_15 * smi_882[k]
                    + f_8 * snh0_665[k]
                    - f_9 * snh1_665[k]
                    + f_3 * pc_x[k] * sni_882[k];

        t_1131[k] = f_15 * smi_883[k]
                    + f_10 * snh0_666[k]
                    - f_11 * snh1_666[k]
                    + f_3 * pc_x[k] * sni_883[k];

        t_1132[k] = f_15 * smi_654[k]
                    + f_3 * pc_z[k] * sni_878[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, pc_y, smi_686, smi_885, smi_886, \
                         snh0_668, snh0_669, snh1_668, snh1_669, sni_882, sni_885, \
                         sni_886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_15 * smi_885[k]
                    + f_10 * snh0_668[k]
                    - f_11 * snh1_668[k]
                    + f_3 * pc_x[k] * sni_885[k];

        t_1134[k] = f_15 * smi_886[k]
                    + f_10 * snh0_669[k]
                    - f_11 * snh1_669[k]
                    + f_3 * pc_x[k] * sni_886[k];

        t_1135[k] = f_16 * smi_686[k]
                    + f_3 * pc_y[k] * sni_882[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, t_1139, pc_x, smi_888, smi_889, smi_890, \
                         smi_891, snh0_671, snh1_671, sni_888, sni_889, sni_890, \
                         sni_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_15 * smi_888[k]
                    + f_10 * snh0_671[k]
                    - f_11 * snh1_671[k]
                    + f_3 * pc_x[k] * sni_888[k];

        t_1137[k] = f_15 * smi_889[k]
                    + f_3 * pc_x[k] * sni_889[k];

        t_1138[k] = f_15 * smi_890[k]
                    + f_3 * pc_x[k] * sni_890[k];

        t_1139[k] = f_15 * smi_891[k]
                    + f_3 * pc_x[k] * sni_891[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, t_1143, pc_x, smi_892, smi_893, smi_894, \
                         smi_895, sni_892, sni_893, sni_894, sni_895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_15 * smi_892[k]
                    + f_3 * pc_x[k] * sni_892[k];

        t_1141[k] = f_15 * smi_893[k]
                    + f_3 * pc_x[k] * sni_893[k];

        t_1142[k] = f_15 * smi_894[k]
                    + f_3 * pc_x[k] * sni_894[k];

        t_1143[k] = f_15 * smi_895[k]
                    + f_3 * pc_x[k] * sni_895[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pc_y, pc_z, smi_665, smi_693, smi_695, \
                         snh0_666, snh0_668, snh1_666, snh1_668, sni_889, \
                         sni_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_16 * smi_693[k]
                    + f_1 * snh0_666[k]
                    - f_2 * snh1_666[k]
                    + f_3 * pc_y[k] * sni_889[k];

        t_1145[k] = f_15 * smi_665[k]
                    + f_3 * pc_z[k] * sni_889[k];

        t_1146[k] = f_16 * smi_695[k]
                    + f_4 * snh0_668[k]
                    - f_5 * snh1_668[k]
                    + f_3 * pc_y[k] * sni_891[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pc_y, smi_696, smi_697, smi_698, snh0_669, \
                         snh0_670, snh0_671, snh1_669, snh1_670, snh1_671, sni_892, sni_893, \
                         sni_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_16 * smi_696[k]
                    + f_6 * snh0_669[k]
                    - f_7 * snh1_669[k]
                    + f_3 * pc_y[k] * sni_892[k];

        t_1148[k] = f_16 * smi_697[k]
                    + f_8 * snh0_670[k]
                    - f_9 * snh1_670[k]
                    + f_3 * pc_y[k] * sni_893[k];

        t_1149[k] = f_16 * smi_698[k]
                    + f_10 * snh0_671[k]
                    - f_11 * snh1_671[k]
                    + f_3 * pc_y[k] * sni_894[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, smi_671, smi_699, smi_896, \
                         snh0_671, snh0_672, snh1_671, snh1_672, sni_895, \
                         sni_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_16 * smi_699[k]
                    + f_3 * pc_y[k] * sni_895[k];

        t_1151[k] = f_15 * smi_671[k]
                    + f_1 * snh0_671[k]
                    - f_2 * snh1_671[k]
                    + f_3 * pc_z[k] * sni_895[k];

        t_1152[k] = f_15 * smi_896[k]
                    + f_1 * snh0_672[k]
                    - f_2 * snh1_672[k]
                    + f_3 * pc_x[k] * sni_896[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pc_x, pc_y, pc_z, smi_672, smi_700, \
                         smi_702, smi_899, snh0_675, snh1_675, sni_896, sni_898, \
                         sni_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_15 * smi_700[k]
                    + f_3 * pc_y[k] * sni_896[k];

        t_1154[k] = f_16 * smi_672[k]
                    + f_3 * pc_z[k] * sni_896[k];

        t_1155[k] = f_15 * smi_899[k]
                    + f_4 * snh0_675[k]
                    - f_5 * snh1_675[k]
                    + f_3 * pc_x[k] * sni_899[k];

        t_1156[k] = f_15 * smi_702[k]
                    + f_3 * pc_y[k] * sni_898[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, pc_x, pc_z, smi_675, smi_901, smi_902, \
                         snh0_677, snh0_678, snh1_677, snh1_678, sni_899, sni_901, \
                         sni_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_15 * smi_901[k]
                    + f_4 * snh0_677[k]
                    - f_5 * snh1_677[k]
                    + f_3 * pc_x[k] * sni_901[k];

        t_1158[k] = f_15 * smi_902[k]
                    + f_6 * snh0_678[k]
                    - f_7 * snh1_678[k]
                    + f_3 * pc_x[k] * sni_902[k];

        t_1159[k] = f_16 * smi_675[k]
                    + f_3 * pc_z[k] * sni_899[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, pc_y, smi_705, smi_905, smi_906, \
                         snh0_681, snh0_682, snh1_681, snh1_682, sni_901, sni_905, \
                         sni_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_15 * smi_705[k]
                    + f_3 * pc_y[k] * sni_901[k];

        t_1161[k] = f_15 * smi_905[k]
                    + f_6 * snh0_681[k]
                    - f_7 * snh1_681[k]
                    + f_3 * pc_x[k] * sni_905[k];

        t_1162[k] = f_15 * smi_906[k]
                    + f_8 * snh0_682[k]
                    - f_9 * snh1_682[k]
                    + f_3 * pc_x[k] * sni_906[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_x, pc_y, pc_z, smi_678, smi_709, smi_908, \
                         snh0_684, snh1_684, sni_902, sni_905, \
                         sni_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * smi_678[k]
                    + f_3 * pc_z[k] * sni_902[k];

        t_1164[k] = f_15 * smi_908[k]
                    + f_8 * snh0_684[k]
                    - f_9 * snh1_684[k]
                    + f_3 * pc_x[k] * sni_908[k];

        t_1165[k] = f_15 * smi_709[k]
                    + f_3 * pc_y[k] * sni_905[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_x, pc_z, smi_682, smi_910, smi_911, \
                         snh0_686, snh0_687, snh1_686, snh1_687, sni_906, sni_910, \
                         sni_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_15 * smi_910[k]
                    + f_8 * snh0_686[k]
                    - f_9 * snh1_686[k]
                    + f_3 * pc_x[k] * sni_910[k];

        t_1167[k] = f_15 * smi_911[k]
                    + f_10 * snh0_687[k]
                    - f_11 * snh1_687[k]
                    + f_3 * pc_x[k] * sni_911[k];

        t_1168[k] = f_16 * smi_682[k]
                    + f_3 * pc_z[k] * sni_906[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, pc_x, pc_y, smi_714, smi_913, smi_914, \
                         snh0_689, snh0_690, snh1_689, snh1_690, sni_910, sni_913, \
                         sni_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_15 * smi_913[k]
                    + f_10 * snh0_689[k]
                    - f_11 * snh1_689[k]
                    + f_3 * pc_x[k] * sni_913[k];

        t_1170[k] = f_15 * smi_914[k]
                    + f_10 * snh0_690[k]
                    - f_11 * snh1_690[k]
                    + f_3 * pc_x[k] * sni_914[k];

        t_1171[k] = f_15 * smi_714[k]
                    + f_3 * pc_y[k] * sni_910[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, t_1175, pc_x, smi_916, smi_917, smi_918, \
                         smi_919, snh0_692, snh1_692, sni_916, sni_917, sni_918, \
                         sni_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_15 * smi_916[k]
                    + f_10 * snh0_692[k]
                    - f_11 * snh1_692[k]
                    + f_3 * pc_x[k] * sni_916[k];

        t_1173[k] = f_15 * smi_917[k]
                    + f_3 * pc_x[k] * sni_917[k];

        t_1174[k] = f_15 * smi_918[k]
                    + f_3 * pc_x[k] * sni_918[k];

        t_1175[k] = f_15 * smi_919[k]
                    + f_3 * pc_x[k] * sni_919[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, t_1179, pc_x, smi_920, smi_921, smi_922, \
                         smi_923, sni_920, sni_921, sni_922, sni_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_15 * smi_920[k]
                    + f_3 * pc_x[k] * sni_920[k];

        t_1177[k] = f_15 * smi_921[k]
                    + f_3 * pc_x[k] * sni_921[k];

        t_1178[k] = f_15 * smi_922[k]
                    + f_3 * pc_x[k] * sni_922[k];

        t_1179[k] = f_15 * smi_923[k]
                    + f_3 * pc_x[k] * sni_923[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pc_y, pc_z, smi_693, smi_721, smi_723, \
                         snh0_687, snh0_689, snh1_687, snh1_689, sni_917, \
                         sni_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_15 * smi_721[k]
                    + f_1 * snh0_687[k]
                    - f_2 * snh1_687[k]
                    + f_3 * pc_y[k] * sni_917[k];

        t_1181[k] = f_16 * smi_693[k]
                    + f_3 * pc_z[k] * sni_917[k];

        t_1182[k] = f_15 * smi_723[k]
                    + f_4 * snh0_689[k]
                    - f_5 * snh1_689[k]
                    + f_3 * pc_y[k] * sni_919[k];
    }

#pragma omp simd aligned(t_1183, t_1184, t_1185, pc_y, smi_724, smi_725, smi_726, snh0_690, \
                         snh0_691, snh0_692, snh1_690, snh1_691, snh1_692, sni_920, sni_921, \
                         sni_922 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_15 * smi_724[k]
                    + f_6 * snh0_690[k]
                    - f_7 * snh1_690[k]
                    + f_3 * pc_y[k] * sni_920[k];

        t_1184[k] = f_15 * smi_725[k]
                    + f_8 * snh0_691[k]
                    - f_9 * snh1_691[k]
                    + f_3 * pc_y[k] * sni_921[k];

        t_1185[k] = f_15 * smi_726[k]
                    + f_10 * snh0_692[k]
                    - f_11 * snh1_692[k]
                    + f_3 * pc_y[k] * sni_922[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, pc_x, pc_y, pc_z, smi_699, smi_727, smi_924, \
                         snh0_692, snh0_693, snh1_692, snh1_693, sni_923, \
                         sni_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_15 * smi_727[k]
                    + f_3 * pc_y[k] * sni_923[k];

        t_1187[k] = f_16 * smi_699[k]
                    + f_1 * snh0_692[k]
                    - f_2 * snh1_692[k]
                    + f_3 * pc_z[k] * sni_923[k];

        t_1188[k] = f_15 * smi_924[k]
                    + f_1 * snh0_693[k]
                    - f_2 * snh1_693[k]
                    + f_3 * pc_x[k] * sni_924[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pc_x, pc_y, pc_z, smi_700, smi_728, \
                         smi_730, smi_927, snh0_696, snh1_696, sni_924, sni_926, \
                         sni_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_14 * smi_728[k]
                    + f_3 * pc_y[k] * sni_924[k];

        t_1190[k] = f_17 * smi_700[k]
                    + f_3 * pc_z[k] * sni_924[k];

        t_1191[k] = f_15 * smi_927[k]
                    + f_4 * snh0_696[k]
                    - f_5 * snh1_696[k]
                    + f_3 * pc_x[k] * sni_927[k];

        t_1192[k] = f_14 * smi_730[k]
                    + f_3 * pc_y[k] * sni_926[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, pc_x, pc_z, smi_703, smi_929, smi_930, \
                         snh0_698, snh0_699, snh1_698, snh1_699, sni_927, sni_929, \
                         sni_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_15 * smi_929[k]
                    + f_4 * snh0_698[k]
                    - f_5 * snh1_698[k]
                    + f_3 * pc_x[k] * sni_929[k];

        t_1194[k] = f_15 * smi_930[k]
                    + f_6 * snh0_699[k]
                    - f_7 * snh1_699[k]
                    + f_3 * pc_x[k] * sni_930[k];

        t_1195[k] = f_17 * smi_703[k]
                    + f_3 * pc_z[k] * sni_927[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, smi_733, smi_933, smi_934, \
                         snh0_702, snh0_703, snh1_702, snh1_703, sni_929, sni_933, \
                         sni_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_14 * smi_733[k]
                    + f_3 * pc_y[k] * sni_929[k];

        t_1197[k] = f_15 * smi_933[k]
                    + f_6 * snh0_702[k]
                    - f_7 * snh1_702[k]
                    + f_3 * pc_x[k] * sni_933[k];

        t_1198[k] = f_15 * smi_934[k]
                    + f_8 * snh0_703[k]
                    - f_9 * snh1_703[k]
                    + f_3 * pc_x[k] * sni_934[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_x, pc_y, pc_z, smi_706, smi_737, smi_936, \
                         snh0_705, snh1_705, sni_930, sni_933, \
                         sni_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_17 * smi_706[k]
                    + f_3 * pc_z[k] * sni_930[k];

        t_1200[k] = f_15 * smi_936[k]
                    + f_8 * snh0_705[k]
                    - f_9 * snh1_705[k]
                    + f_3 * pc_x[k] * sni_936[k];

        t_1201[k] = f_14 * smi_737[k]
                    + f_3 * pc_y[k] * sni_933[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pc_x, pc_z, smi_710, smi_938, smi_939, \
                         snh0_707, snh0_708, snh1_707, snh1_708, sni_934, sni_938, \
                         sni_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_15 * smi_938[k]
                    + f_8 * snh0_707[k]
                    - f_9 * snh1_707[k]
                    + f_3 * pc_x[k] * sni_938[k];

        t_1203[k] = f_15 * smi_939[k]
                    + f_10 * snh0_708[k]
                    - f_11 * snh1_708[k]
                    + f_3 * pc_x[k] * sni_939[k];

        t_1204[k] = f_17 * smi_710[k]
                    + f_3 * pc_z[k] * sni_934[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, pc_x, pc_y, smi_742, smi_941, smi_942, \
                         snh0_710, snh0_711, snh1_710, snh1_711, sni_938, sni_941, \
                         sni_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_15 * smi_941[k]
                    + f_10 * snh0_710[k]
                    - f_11 * snh1_710[k]
                    + f_3 * pc_x[k] * sni_941[k];

        t_1206[k] = f_15 * smi_942[k]
                    + f_10 * snh0_711[k]
                    - f_11 * snh1_711[k]
                    + f_3 * pc_x[k] * sni_942[k];

        t_1207[k] = f_14 * smi_742[k]
                    + f_3 * pc_y[k] * sni_938[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, t_1211, pc_x, smi_944, smi_945, smi_946, \
                         smi_947, snh0_713, snh1_713, sni_944, sni_945, sni_946, \
                         sni_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_15 * smi_944[k]
                    + f_10 * snh0_713[k]
                    - f_11 * snh1_713[k]
                    + f_3 * pc_x[k] * sni_944[k];

        t_1209[k] = f_15 * smi_945[k]
                    + f_3 * pc_x[k] * sni_945[k];

        t_1210[k] = f_15 * smi_946[k]
                    + f_3 * pc_x[k] * sni_946[k];

        t_1211[k] = f_15 * smi_947[k]
                    + f_3 * pc_x[k] * sni_947[k];
    }

#pragma omp simd aligned(t_1212, t_1213, t_1214, t_1215, pc_x, smi_948, smi_949, smi_950, \
                         smi_951, sni_948, sni_949, sni_950, sni_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1212[k] = f_15 * smi_948[k]
                    + f_3 * pc_x[k] * sni_948[k];

        t_1213[k] = f_15 * smi_949[k]
                    + f_3 * pc_x[k] * sni_949[k];

        t_1214[k] = f_15 * smi_950[k]
                    + f_3 * pc_x[k] * sni_950[k];

        t_1215[k] = f_15 * smi_951[k]
                    + f_3 * pc_x[k] * sni_951[k];
    }

#pragma omp simd aligned(t_1216, t_1217, t_1218, pc_y, pc_z, smi_721, smi_749, smi_751, \
                         snh0_708, snh0_710, snh1_708, snh1_710, sni_945, \
                         sni_947 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1216[k] = f_14 * smi_749[k]
                    + f_1 * snh0_708[k]
                    - f_2 * snh1_708[k]
                    + f_3 * pc_y[k] * sni_945[k];

        t_1217[k] = f_17 * smi_721[k]
                    + f_3 * pc_z[k] * sni_945[k];

        t_1218[k] = f_14 * smi_751[k]
                    + f_4 * snh0_710[k]
                    - f_5 * snh1_710[k]
                    + f_3 * pc_y[k] * sni_947[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, pc_y, smi_752, smi_753, smi_754, snh0_711, \
                         snh0_712, snh0_713, snh1_711, snh1_712, snh1_713, sni_948, sni_949, \
                         sni_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = f_14 * smi_752[k]
                    + f_6 * snh0_711[k]
                    - f_7 * snh1_711[k]
                    + f_3 * pc_y[k] * sni_948[k];

        t_1220[k] = f_14 * smi_753[k]
                    + f_8 * snh0_712[k]
                    - f_9 * snh1_712[k]
                    + f_3 * pc_y[k] * sni_949[k];

        t_1221[k] = f_14 * smi_754[k]
                    + f_10 * snh0_713[k]
                    - f_11 * snh1_713[k]
                    + f_3 * pc_y[k] * sni_950[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smk0,
                                                           const size_t smi, const size_t smk1,
                                                           const size_t snh0, const size_t snh1,
                                                           const size_t sni, const size_t ncols,
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
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smk0_972 = buffer.data(smk0 + 972);
    const auto *smk0_975 = buffer.data(smk0 + 975);
    const auto *smk0_977 = buffer.data(smk0 + 977);
    const auto *smk0_978 = buffer.data(smk0 + 978);
    const auto *smk0_981 = buffer.data(smk0 + 981);
    const auto *smk0_982 = buffer.data(smk0 + 982);
    const auto *smk0_984 = buffer.data(smk0 + 984);
    const auto *smk0_986 = buffer.data(smk0 + 986);
    const auto *smk0_987 = buffer.data(smk0 + 987);
    const auto *smk0_989 = buffer.data(smk0 + 989);
    const auto *smk0_990 = buffer.data(smk0 + 990);
    const auto *smk0_992 = buffer.data(smk0 + 992);
    const auto *smk0_1007 = buffer.data(smk0 + 1007);
    const auto *smk0_1008 = buffer.data(smk0 + 1008);
    const auto *smk0_1011 = buffer.data(smk0 + 1011);

    const auto *smi_727 = buffer.data(smi + 727);
    const auto *smi_728 = buffer.data(smi + 728);
    const auto *smi_731 = buffer.data(smi + 731);
    const auto *smi_734 = buffer.data(smi + 734);
    const auto *smi_738 = buffer.data(smi + 738);
    const auto *smi_749 = buffer.data(smi + 749);
    const auto *smi_755 = buffer.data(smi + 755);
    const auto *smi_756 = buffer.data(smi + 756);
    const auto *smi_757 = buffer.data(smi + 757);
    const auto *smi_758 = buffer.data(smi + 758);
    const auto *smi_759 = buffer.data(smi + 759);
    const auto *smi_761 = buffer.data(smi + 761);
    const auto *smi_762 = buffer.data(smi + 762);
    const auto *smi_764 = buffer.data(smi + 764);
    const auto *smi_765 = buffer.data(smi + 765);
    const auto *smi_766 = buffer.data(smi + 766);
    const auto *smi_768 = buffer.data(smi + 768);
    const auto *smi_769 = buffer.data(smi + 769);
    const auto *smi_770 = buffer.data(smi + 770);
    const auto *smi_777 = buffer.data(smi + 777);
    const auto *smi_779 = buffer.data(smi + 779);
    const auto *smi_780 = buffer.data(smi + 780);
    const auto *smi_781 = buffer.data(smi + 781);
    const auto *smi_782 = buffer.data(smi + 782);
    const auto *smi_783 = buffer.data(smi + 783);
    const auto *smi_784 = buffer.data(smi + 784);
    const auto *smi_786 = buffer.data(smi + 786);
    const auto *smi_789 = buffer.data(smi + 789);
    const auto *smi_793 = buffer.data(smi + 793);
    const auto *smi_798 = buffer.data(smi + 798);
    const auto *smi_805 = buffer.data(smi + 805);
    const auto *smi_807 = buffer.data(smi + 807);
    const auto *smi_808 = buffer.data(smi + 808);
    const auto *smi_809 = buffer.data(smi + 809);
    const auto *smi_810 = buffer.data(smi + 810);
    const auto *smi_811 = buffer.data(smi + 811);
    const auto *smi_812 = buffer.data(smi + 812);
    const auto *smi_814 = buffer.data(smi + 814);
    const auto *smi_973 = buffer.data(smi + 973);
    const auto *smi_974 = buffer.data(smi + 974);
    const auto *smi_975 = buffer.data(smi + 975);
    const auto *smi_976 = buffer.data(smi + 976);
    const auto *smi_977 = buffer.data(smi + 977);
    const auto *smi_978 = buffer.data(smi + 978);
    const auto *smi_979 = buffer.data(smi + 979);
    const auto *smi_980 = buffer.data(smi + 980);
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
    const auto *smi_1011 = buffer.data(smi + 1011);
    const auto *smi_1013 = buffer.data(smi + 1013);
    const auto *smi_1014 = buffer.data(smi + 1014);
    const auto *smi_1017 = buffer.data(smi + 1017);
    const auto *smi_1018 = buffer.data(smi + 1018);
    const auto *smi_1020 = buffer.data(smi + 1020);
    const auto *smi_1022 = buffer.data(smi + 1022);
    const auto *smi_1023 = buffer.data(smi + 1023);
    const auto *smi_1025 = buffer.data(smi + 1025);
    const auto *smi_1026 = buffer.data(smi + 1026);
    const auto *smi_1028 = buffer.data(smi + 1028);
    const auto *smi_1029 = buffer.data(smi + 1029);
    const auto *smi_1030 = buffer.data(smi + 1030);
    const auto *smi_1031 = buffer.data(smi + 1031);
    const auto *smi_1032 = buffer.data(smi + 1032);
    const auto *smi_1033 = buffer.data(smi + 1033);
    const auto *smi_1034 = buffer.data(smi + 1034);
    const auto *smi_1035 = buffer.data(smi + 1035);

    const auto *smk1_972 = buffer.data(smk1 + 972);
    const auto *smk1_975 = buffer.data(smk1 + 975);
    const auto *smk1_977 = buffer.data(smk1 + 977);
    const auto *smk1_978 = buffer.data(smk1 + 978);
    const auto *smk1_981 = buffer.data(smk1 + 981);
    const auto *smk1_982 = buffer.data(smk1 + 982);
    const auto *smk1_984 = buffer.data(smk1 + 984);
    const auto *smk1_986 = buffer.data(smk1 + 986);
    const auto *smk1_987 = buffer.data(smk1 + 987);
    const auto *smk1_989 = buffer.data(smk1 + 989);
    const auto *smk1_990 = buffer.data(smk1 + 990);
    const auto *smk1_992 = buffer.data(smk1 + 992);
    const auto *smk1_1007 = buffer.data(smk1 + 1007);
    const auto *smk1_1008 = buffer.data(smk1 + 1008);
    const auto *smk1_1011 = buffer.data(smk1 + 1011);

    const auto *snh0_713 = buffer.data(snh0 + 713);
    const auto *snh0_729 = buffer.data(snh0 + 729);
    const auto *snh0_731 = buffer.data(snh0 + 731);
    const auto *snh0_732 = buffer.data(snh0 + 732);
    const auto *snh0_733 = buffer.data(snh0 + 733);
    const auto *snh0_734 = buffer.data(snh0 + 734);
    const auto *snh0_735 = buffer.data(snh0 + 735);
    const auto *snh0_738 = buffer.data(snh0 + 738);
    const auto *snh0_740 = buffer.data(snh0 + 740);
    const auto *snh0_741 = buffer.data(snh0 + 741);
    const auto *snh0_744 = buffer.data(snh0 + 744);
    const auto *snh0_745 = buffer.data(snh0 + 745);
    const auto *snh0_747 = buffer.data(snh0 + 747);
    const auto *snh0_749 = buffer.data(snh0 + 749);
    const auto *snh0_750 = buffer.data(snh0 + 750);
    const auto *snh0_752 = buffer.data(snh0 + 752);
    const auto *snh0_753 = buffer.data(snh0 + 753);
    const auto *snh0_754 = buffer.data(snh0 + 754);
    const auto *snh0_755 = buffer.data(snh0 + 755);
    const auto *snh0_756 = buffer.data(snh0 + 756);
    const auto *snh0_759 = buffer.data(snh0 + 759);
    const auto *snh0_761 = buffer.data(snh0 + 761);
    const auto *snh0_762 = buffer.data(snh0 + 762);
    const auto *snh0_765 = buffer.data(snh0 + 765);
    const auto *snh0_766 = buffer.data(snh0 + 766);
    const auto *snh0_768 = buffer.data(snh0 + 768);
    const auto *snh0_770 = buffer.data(snh0 + 770);
    const auto *snh0_771 = buffer.data(snh0 + 771);
    const auto *snh0_773 = buffer.data(snh0 + 773);
    const auto *snh0_774 = buffer.data(snh0 + 774);
    const auto *snh0_775 = buffer.data(snh0 + 775);
    const auto *snh0_776 = buffer.data(snh0 + 776);

    const auto *snh1_713 = buffer.data(snh1 + 713);
    const auto *snh1_729 = buffer.data(snh1 + 729);
    const auto *snh1_731 = buffer.data(snh1 + 731);
    const auto *snh1_732 = buffer.data(snh1 + 732);
    const auto *snh1_733 = buffer.data(snh1 + 733);
    const auto *snh1_734 = buffer.data(snh1 + 734);
    const auto *snh1_735 = buffer.data(snh1 + 735);
    const auto *snh1_738 = buffer.data(snh1 + 738);
    const auto *snh1_740 = buffer.data(snh1 + 740);
    const auto *snh1_741 = buffer.data(snh1 + 741);
    const auto *snh1_744 = buffer.data(snh1 + 744);
    const auto *snh1_745 = buffer.data(snh1 + 745);
    const auto *snh1_747 = buffer.data(snh1 + 747);
    const auto *snh1_749 = buffer.data(snh1 + 749);
    const auto *snh1_750 = buffer.data(snh1 + 750);
    const auto *snh1_752 = buffer.data(snh1 + 752);
    const auto *snh1_753 = buffer.data(snh1 + 753);
    const auto *snh1_754 = buffer.data(snh1 + 754);
    const auto *snh1_755 = buffer.data(snh1 + 755);
    const auto *snh1_756 = buffer.data(snh1 + 756);
    const auto *snh1_759 = buffer.data(snh1 + 759);
    const auto *snh1_761 = buffer.data(snh1 + 761);
    const auto *snh1_762 = buffer.data(snh1 + 762);
    const auto *snh1_765 = buffer.data(snh1 + 765);
    const auto *snh1_766 = buffer.data(snh1 + 766);
    const auto *snh1_768 = buffer.data(snh1 + 768);
    const auto *snh1_770 = buffer.data(snh1 + 770);
    const auto *snh1_771 = buffer.data(snh1 + 771);
    const auto *snh1_773 = buffer.data(snh1 + 773);
    const auto *snh1_774 = buffer.data(snh1 + 774);
    const auto *snh1_775 = buffer.data(snh1 + 775);
    const auto *snh1_776 = buffer.data(snh1 + 776);

    const auto *sni_951 = buffer.data(sni + 951);
    const auto *sni_952 = buffer.data(sni + 952);
    const auto *sni_954 = buffer.data(sni + 954);
    const auto *sni_955 = buffer.data(sni + 955);
    const auto *sni_957 = buffer.data(sni + 957);
    const auto *sni_958 = buffer.data(sni + 958);
    const auto *sni_961 = buffer.data(sni + 961);
    const auto *sni_962 = buffer.data(sni + 962);
    const auto *sni_966 = buffer.data(sni + 966);
    const auto *sni_973 = buffer.data(sni + 973);
    const auto *sni_974 = buffer.data(sni + 974);
    const auto *sni_975 = buffer.data(sni + 975);
    const auto *sni_976 = buffer.data(sni + 976);
    const auto *sni_977 = buffer.data(sni + 977);
    const auto *sni_978 = buffer.data(sni + 978);
    const auto *sni_979 = buffer.data(sni + 979);
    const auto *sni_980 = buffer.data(sni + 980);
    const auto *sni_982 = buffer.data(sni + 982);
    const auto *sni_983 = buffer.data(sni + 983);
    const auto *sni_985 = buffer.data(sni + 985);
    const auto *sni_986 = buffer.data(sni + 986);
    const auto *sni_989 = buffer.data(sni + 989);
    const auto *sni_990 = buffer.data(sni + 990);
    const auto *sni_992 = buffer.data(sni + 992);
    const auto *sni_994 = buffer.data(sni + 994);
    const auto *sni_995 = buffer.data(sni + 995);
    const auto *sni_997 = buffer.data(sni + 997);
    const auto *sni_998 = buffer.data(sni + 998);
    const auto *sni_1000 = buffer.data(sni + 1000);
    const auto *sni_1001 = buffer.data(sni + 1001);
    const auto *sni_1002 = buffer.data(sni + 1002);
    const auto *sni_1003 = buffer.data(sni + 1003);
    const auto *sni_1004 = buffer.data(sni + 1004);
    const auto *sni_1005 = buffer.data(sni + 1005);
    const auto *sni_1006 = buffer.data(sni + 1006);
    const auto *sni_1007 = buffer.data(sni + 1007);
    const auto *sni_1008 = buffer.data(sni + 1008);
    const auto *sni_1010 = buffer.data(sni + 1010);
    const auto *sni_1011 = buffer.data(sni + 1011);
    const auto *sni_1013 = buffer.data(sni + 1013);
    const auto *sni_1014 = buffer.data(sni + 1014);
    const auto *sni_1017 = buffer.data(sni + 1017);
    const auto *sni_1018 = buffer.data(sni + 1018);
    const auto *sni_1020 = buffer.data(sni + 1020);
    const auto *sni_1022 = buffer.data(sni + 1022);
    const auto *sni_1023 = buffer.data(sni + 1023);
    const auto *sni_1025 = buffer.data(sni + 1025);
    const auto *sni_1026 = buffer.data(sni + 1026);
    const auto *sni_1028 = buffer.data(sni + 1028);
    const auto *sni_1029 = buffer.data(sni + 1029);
    const auto *sni_1030 = buffer.data(sni + 1030);
    const auto *sni_1031 = buffer.data(sni + 1031);
    const auto *sni_1032 = buffer.data(sni + 1032);
    const auto *sni_1033 = buffer.data(sni + 1033);
    const auto *sni_1034 = buffer.data(sni + 1034);
    const auto *sni_1035 = buffer.data(sni + 1035);
    const auto *sni_1036 = buffer.data(sni + 1036);
    const auto *sni_1038 = buffer.data(sni + 1038);

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pb_y, pc_y, pc_z, smk0_972, smi_727, \
                         smi_755, smi_756, smk1_972, snh0_713, snh1_713, sni_951, \
                         sni_952 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_14 * smi_755[k]
                    + f_3 * pc_y[k] * sni_951[k];

        t_1223[k] = f_17 * smi_727[k]
                    + f_1 * snh0_713[k]
                    - f_2 * snh1_713[k]
                    + f_3 * pc_z[k] * sni_951[k];

        t_1224[k] = pb_y[k] * smk0_972[k]
                    - f_12 * pc_y[k] * smk1_972[k];

        t_1225[k] = f_13 * smi_756[k]
                    + f_3 * pc_y[k] * sni_952[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, pb_y, pc_y, pc_z, smk0_975, smk0_977, \
                         smi_728, smi_757, smi_758, smk1_975, smk1_977, sni_952, \
                         sni_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_21 * smi_728[k]
                    + f_3 * pc_z[k] * sni_952[k];

        t_1227[k] = pb_y[k] * smk0_975[k]
                    + f_14 * smi_757[k]
                    - f_12 * pc_y[k] * smk1_975[k];

        t_1228[k] = f_13 * smi_758[k]
                    + f_3 * pc_y[k] * sni_954[k];

        t_1229[k] = pb_y[k] * smk0_977[k]
                    - f_12 * pc_y[k] * smk1_977[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pb_y, pc_y, pc_z, smk0_978, smk0_981, \
                         smi_731, smi_759, smi_761, smk1_978, smk1_981, sni_955, \
                         sni_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = pb_y[k] * smk0_978[k]
                    + f_15 * smi_759[k]
                    - f_12 * pc_y[k] * smk1_978[k];

        t_1231[k] = f_21 * smi_731[k]
                    + f_3 * pc_z[k] * sni_955[k];

        t_1232[k] = f_13 * smi_761[k]
                    + f_3 * pc_y[k] * sni_957[k];

        t_1233[k] = pb_y[k] * smk0_981[k]
                    - f_12 * pc_y[k] * smk1_981[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, pb_y, pc_y, pc_z, smk0_982, smk0_984, \
                         smi_734, smi_762, smi_764, smk1_982, smk1_984, \
                         sni_958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = pb_y[k] * smk0_982[k]
                    + f_16 * smi_762[k]
                    - f_12 * pc_y[k] * smk1_982[k];

        t_1235[k] = f_21 * smi_734[k]
                    + f_3 * pc_z[k] * sni_958[k];

        t_1236[k] = pb_y[k] * smk0_984[k]
                    + f_14 * smi_764[k]
                    - f_12 * pc_y[k] * smk1_984[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, pb_y, pc_y, pc_z, smk0_986, smk0_987, \
                         smi_738, smi_765, smi_766, smk1_986, smk1_987, sni_961, \
                         sni_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_13 * smi_765[k]
                    + f_3 * pc_y[k] * sni_961[k];

        t_1238[k] = pb_y[k] * smk0_986[k]
                    - f_12 * pc_y[k] * smk1_986[k];

        t_1239[k] = pb_y[k] * smk0_987[k]
                    + f_17 * smi_766[k]
                    - f_12 * pc_y[k] * smk1_987[k];

        t_1240[k] = f_21 * smi_738[k]
                    + f_3 * pc_z[k] * sni_962[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, t_1244, pb_y, pc_y, smk0_989, smk0_990, \
                         smk0_992, smi_768, smi_769, smi_770, smk1_989, smk1_990, smk1_992, \
                         sni_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = pb_y[k] * smk0_989[k]
                    + f_15 * smi_768[k]
                    - f_12 * pc_y[k] * smk1_989[k];

        t_1242[k] = pb_y[k] * smk0_990[k]
                    + f_14 * smi_769[k]
                    - f_12 * pc_y[k] * smk1_990[k];

        t_1243[k] = f_13 * smi_770[k]
                    + f_3 * pc_y[k] * sni_966[k];

        t_1244[k] = pb_y[k] * smk0_992[k]
                    - f_12 * pc_y[k] * smk1_992[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, t_1248, t_1249, pc_x, smi_973, smi_974, \
                         smi_975, smi_976, smi_977, sni_973, sni_974, sni_975, sni_976, \
                         sni_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = f_15 * smi_973[k]
                    + f_3 * pc_x[k] * sni_973[k];

        t_1246[k] = f_15 * smi_974[k]
                    + f_3 * pc_x[k] * sni_974[k];

        t_1247[k] = f_15 * smi_975[k]
                    + f_3 * pc_x[k] * sni_975[k];

        t_1248[k] = f_15 * smi_976[k]
                    + f_3 * pc_x[k] * sni_976[k];

        t_1249[k] = f_15 * smi_977[k]
                    + f_3 * pc_x[k] * sni_977[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pc_x, pc_y, pc_z, smi_749, smi_777, \
                         smi_978, smi_979, snh0_729, snh1_729, sni_973, sni_978, \
                         sni_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_15 * smi_978[k]
                    + f_3 * pc_x[k] * sni_978[k];

        t_1251[k] = f_15 * smi_979[k]
                    + f_3 * pc_x[k] * sni_979[k];

        t_1252[k] = f_13 * smi_777[k]
                    + f_1 * snh0_729[k]
                    - f_2 * snh1_729[k]
                    + f_3 * pc_y[k] * sni_973[k];

        t_1253[k] = f_21 * smi_749[k]
                    + f_3 * pc_z[k] * sni_973[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, smi_779, smi_780, smi_781, snh0_731, \
                         snh0_732, snh0_733, snh1_731, snh1_732, snh1_733, sni_975, sni_976, \
                         sni_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_13 * smi_779[k]
                    + f_4 * snh0_731[k]
                    - f_5 * snh1_731[k]
                    + f_3 * pc_y[k] * sni_975[k];

        t_1255[k] = f_13 * smi_780[k]
                    + f_6 * snh0_732[k]
                    - f_7 * snh1_732[k]
                    + f_3 * pc_y[k] * sni_976[k];

        t_1256[k] = f_13 * smi_781[k]
                    + f_8 * snh0_733[k]
                    - f_9 * snh1_733[k]
                    + f_3 * pc_y[k] * sni_977[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pb_y, pc_y, smk0_1007, smi_782, smi_783, \
                         smk1_1007, snh0_734, snh1_734, sni_978, \
                         sni_979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_13 * smi_782[k]
                    + f_10 * snh0_734[k]
                    - f_11 * snh1_734[k]
                    + f_3 * pc_y[k] * sni_978[k];

        t_1258[k] = f_13 * smi_783[k]
                    + f_3 * pc_y[k] * sni_979[k];

        t_1259[k] = pb_y[k] * smk0_1007[k]
                    - f_12 * pc_y[k] * smk1_1007[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pc_x, pc_y, pc_z, smi_756, smi_980, \
                         smi_983, snh0_735, snh0_738, snh1_735, snh1_738, sni_980, \
                         sni_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_15 * smi_980[k]
                    + f_1 * snh0_735[k]
                    - f_2 * snh1_735[k]
                    + f_3 * pc_x[k] * sni_980[k];

        t_1261[k] = f_3 * pc_y[k] * sni_980[k];

        t_1262[k] = f_20 * smi_756[k]
                    + f_3 * pc_z[k] * sni_980[k];

        t_1263[k] = f_15 * smi_983[k]
                    + f_4 * snh0_738[k]
                    - f_5 * snh1_738[k]
                    + f_3 * pc_x[k] * sni_983[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, pc_x, pc_y, smi_985, smi_986, snh0_740, \
                         snh0_741, snh1_740, snh1_741, sni_982, sni_985, \
                         sni_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_3 * pc_y[k] * sni_982[k];

        t_1265[k] = f_15 * smi_985[k]
                    + f_4 * snh0_740[k]
                    - f_5 * snh1_740[k]
                    + f_3 * pc_x[k] * sni_985[k];

        t_1266[k] = f_15 * smi_986[k]
                    + f_6 * snh0_741[k]
                    - f_7 * snh1_741[k]
                    + f_3 * pc_x[k] * sni_986[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, pc_x, pc_y, pc_z, smi_759, smi_989, snh0_744, \
                         snh1_744, sni_983, sni_985, sni_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_20 * smi_759[k]
                    + f_3 * pc_z[k] * sni_983[k];

        t_1268[k] = f_3 * pc_y[k] * sni_985[k];

        t_1269[k] = f_15 * smi_989[k]
                    + f_6 * snh0_744[k]
                    - f_7 * snh1_744[k]
                    + f_3 * pc_x[k] * sni_989[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, pc_x, pc_z, smi_762, smi_990, smi_992, \
                         snh0_745, snh0_747, snh1_745, snh1_747, sni_986, sni_990, \
                         sni_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = f_15 * smi_990[k]
                    + f_8 * snh0_745[k]
                    - f_9 * snh1_745[k]
                    + f_3 * pc_x[k] * sni_990[k];

        t_1271[k] = f_20 * smi_762[k]
                    + f_3 * pc_z[k] * sni_986[k];

        t_1272[k] = f_15 * smi_992[k]
                    + f_8 * snh0_747[k]
                    - f_9 * snh1_747[k]
                    + f_3 * pc_x[k] * sni_992[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, pc_x, pc_y, smi_994, smi_995, snh0_749, \
                         snh0_750, snh1_749, snh1_750, sni_989, sni_994, \
                         sni_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_3 * pc_y[k] * sni_989[k];

        t_1274[k] = f_15 * smi_994[k]
                    + f_8 * snh0_749[k]
                    - f_9 * snh1_749[k]
                    + f_3 * pc_x[k] * sni_994[k];

        t_1275[k] = f_15 * smi_995[k]
                    + f_10 * snh0_750[k]
                    - f_11 * snh1_750[k]
                    + f_3 * pc_x[k] * sni_995[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, pc_x, pc_z, smi_766, smi_997, smi_998, \
                         snh0_752, snh0_753, snh1_752, snh1_753, sni_990, sni_997, \
                         sni_998 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_20 * smi_766[k]
                    + f_3 * pc_z[k] * sni_990[k];

        t_1277[k] = f_15 * smi_997[k]
                    + f_10 * snh0_752[k]
                    - f_11 * snh1_752[k]
                    + f_3 * pc_x[k] * sni_997[k];

        t_1278[k] = f_15 * smi_998[k]
                    + f_10 * snh0_753[k]
                    - f_11 * snh1_753[k]
                    + f_3 * pc_x[k] * sni_998[k];
    }

#pragma omp simd aligned(t_1279, t_1280, t_1281, t_1282, pc_x, pc_y, smi_1000, smi_1001, \
                         smi_1002, snh0_755, snh1_755, sni_994, sni_1000, sni_1001, \
                         sni_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1279[k] = f_3 * pc_y[k] * sni_994[k];

        t_1280[k] = f_15 * smi_1000[k]
                    + f_10 * snh0_755[k]
                    - f_11 * snh1_755[k]
                    + f_3 * pc_x[k] * sni_1000[k];

        t_1281[k] = f_15 * smi_1001[k]
                    + f_3 * pc_x[k] * sni_1001[k];

        t_1282[k] = f_15 * smi_1002[k]
                    + f_3 * pc_x[k] * sni_1002[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, t_1286, t_1287, pc_x, smi_1003, smi_1004, \
                         smi_1005, smi_1006, smi_1007, sni_1003, sni_1004, sni_1005, sni_1006, \
                         sni_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = f_15 * smi_1003[k]
                    + f_3 * pc_x[k] * sni_1003[k];

        t_1284[k] = f_15 * smi_1004[k]
                    + f_3 * pc_x[k] * sni_1004[k];

        t_1285[k] = f_15 * smi_1005[k]
                    + f_3 * pc_x[k] * sni_1005[k];

        t_1286[k] = f_15 * smi_1006[k]
                    + f_3 * pc_x[k] * sni_1006[k];

        t_1287[k] = f_15 * smi_1007[k]
                    + f_3 * pc_x[k] * sni_1007[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pc_y, pc_z, smi_777, snh0_750, \
                         snh0_752, snh0_753, snh1_750, snh1_752, snh1_753, sni_1001, sni_1003, \
                         sni_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = f_1 * snh0_750[k]
                    - f_2 * snh1_750[k]
                    + f_3 * pc_y[k] * sni_1001[k];

        t_1289[k] = f_20 * smi_777[k]
                    + f_3 * pc_z[k] * sni_1001[k];

        t_1290[k] = f_4 * snh0_752[k]
                    - f_5 * snh1_752[k]
                    + f_3 * pc_y[k] * sni_1003[k];

        t_1291[k] = f_6 * snh0_753[k]
                    - f_7 * snh1_753[k]
                    + f_3 * pc_y[k] * sni_1004[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, t_1295, pc_y, pc_z, smi_783, snh0_754, \
                         snh0_755, snh1_754, snh1_755, sni_1005, sni_1006, \
                         sni_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_8 * snh0_754[k]
                    - f_9 * snh1_754[k]
                    + f_3 * pc_y[k] * sni_1005[k];

        t_1293[k] = f_10 * snh0_755[k]
                    - f_11 * snh1_755[k]
                    + f_3 * pc_y[k] * sni_1006[k];

        t_1294[k] = f_3 * pc_y[k] * sni_1007[k];

        t_1295[k] = f_20 * smi_783[k]
                    + f_1 * snh0_755[k]
                    - f_2 * snh1_755[k]
                    + f_3 * pc_z[k] * sni_1007[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, t_1299, pc_x, pc_y, pc_z, smi_784, smi_1008, \
                         smi_1011, snh0_756, snh0_759, snh1_756, snh1_759, sni_1008, \
                         sni_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_14 * smi_1008[k]
                    + f_1 * snh0_756[k]
                    - f_2 * snh1_756[k]
                    + f_3 * pc_x[k] * sni_1008[k];

        t_1297[k] = f_19 * smi_784[k]
                    + f_3 * pc_y[k] * sni_1008[k];

        t_1298[k] = f_3 * pc_z[k] * sni_1008[k];

        t_1299[k] = f_14 * smi_1011[k]
                    + f_4 * snh0_759[k]
                    - f_5 * snh1_759[k]
                    + f_3 * pc_x[k] * sni_1011[k];
    }

#pragma omp simd aligned(t_1300, t_1301, t_1302, pc_x, pc_y, smi_786, smi_1013, smi_1014, \
                         snh0_761, snh0_762, snh1_761, snh1_762, sni_1010, sni_1013, \
                         sni_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1300[k] = f_19 * smi_786[k]
                    + f_3 * pc_y[k] * sni_1010[k];

        t_1301[k] = f_14 * smi_1013[k]
                    + f_4 * snh0_761[k]
                    - f_5 * snh1_761[k]
                    + f_3 * pc_x[k] * sni_1013[k];

        t_1302[k] = f_14 * smi_1014[k]
                    + f_6 * snh0_762[k]
                    - f_7 * snh1_762[k]
                    + f_3 * pc_x[k] * sni_1014[k];
    }

#pragma omp simd aligned(t_1303, t_1304, t_1305, pc_x, pc_y, pc_z, smi_789, smi_1017, \
                         snh0_765, snh1_765, sni_1011, sni_1013, \
                         sni_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1303[k] = f_3 * pc_z[k] * sni_1011[k];

        t_1304[k] = f_19 * smi_789[k]
                    + f_3 * pc_y[k] * sni_1013[k];

        t_1305[k] = f_14 * smi_1017[k]
                    + f_6 * snh0_765[k]
                    - f_7 * snh1_765[k]
                    + f_3 * pc_x[k] * sni_1017[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, pc_x, pc_z, smi_1018, smi_1020, snh0_766, \
                         snh0_768, snh1_766, snh1_768, sni_1014, sni_1018, \
                         sni_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_14 * smi_1018[k]
                    + f_8 * snh0_766[k]
                    - f_9 * snh1_766[k]
                    + f_3 * pc_x[k] * sni_1018[k];

        t_1307[k] = f_3 * pc_z[k] * sni_1014[k];

        t_1308[k] = f_14 * smi_1020[k]
                    + f_8 * snh0_768[k]
                    - f_9 * snh1_768[k]
                    + f_3 * pc_x[k] * sni_1020[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, pc_x, pc_y, smi_793, smi_1022, smi_1023, \
                         snh0_770, snh0_771, snh1_770, snh1_771, sni_1017, sni_1022, \
                         sni_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = f_19 * smi_793[k]
                    + f_3 * pc_y[k] * sni_1017[k];

        t_1310[k] = f_14 * smi_1022[k]
                    + f_8 * snh0_770[k]
                    - f_9 * snh1_770[k]
                    + f_3 * pc_x[k] * sni_1022[k];

        t_1311[k] = f_14 * smi_1023[k]
                    + f_10 * snh0_771[k]
                    - f_11 * snh1_771[k]
                    + f_3 * pc_x[k] * sni_1023[k];
    }

#pragma omp simd aligned(t_1312, t_1313, t_1314, pc_x, pc_z, smi_1025, smi_1026, snh0_773, \
                         snh0_774, snh1_773, snh1_774, sni_1018, sni_1025, \
                         sni_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1312[k] = f_3 * pc_z[k] * sni_1018[k];

        t_1313[k] = f_14 * smi_1025[k]
                    + f_10 * snh0_773[k]
                    - f_11 * snh1_773[k]
                    + f_3 * pc_x[k] * sni_1025[k];

        t_1314[k] = f_14 * smi_1026[k]
                    + f_10 * snh0_774[k]
                    - f_11 * snh1_774[k]
                    + f_3 * pc_x[k] * sni_1026[k];
    }

#pragma omp simd aligned(t_1315, t_1316, t_1317, t_1318, pc_x, pc_y, smi_798, smi_1028, \
                         smi_1029, smi_1030, snh0_776, snh1_776, sni_1022, sni_1028, sni_1029, \
                         sni_1030 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1315[k] = f_19 * smi_798[k]
                    + f_3 * pc_y[k] * sni_1022[k];

        t_1316[k] = f_14 * smi_1028[k]
                    + f_10 * snh0_776[k]
                    - f_11 * snh1_776[k]
                    + f_3 * pc_x[k] * sni_1028[k];

        t_1317[k] = f_14 * smi_1029[k]
                    + f_3 * pc_x[k] * sni_1029[k];

        t_1318[k] = f_14 * smi_1030[k]
                    + f_3 * pc_x[k] * sni_1030[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, t_1322, t_1323, pc_x, smi_1031, smi_1032, \
                         smi_1033, smi_1034, smi_1035, sni_1031, sni_1032, sni_1033, sni_1034, \
                         sni_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_14 * smi_1031[k]
                    + f_3 * pc_x[k] * sni_1031[k];

        t_1320[k] = f_14 * smi_1032[k]
                    + f_3 * pc_x[k] * sni_1032[k];

        t_1321[k] = f_14 * smi_1033[k]
                    + f_3 * pc_x[k] * sni_1033[k];

        t_1322[k] = f_14 * smi_1034[k]
                    + f_3 * pc_x[k] * sni_1034[k];

        t_1323[k] = f_14 * smi_1035[k]
                    + f_3 * pc_x[k] * sni_1035[k];
    }

#pragma omp simd aligned(t_1324, t_1325, t_1326, pc_y, pc_z, smi_805, smi_807, snh0_771, \
                         snh0_773, snh1_771, snh1_773, sni_1029, \
                         sni_1031 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1324[k] = f_19 * smi_805[k]
                    + f_1 * snh0_771[k]
                    - f_2 * snh1_771[k]
                    + f_3 * pc_y[k] * sni_1029[k];

        t_1325[k] = f_3 * pc_z[k] * sni_1029[k];

        t_1326[k] = f_19 * smi_807[k]
                    + f_4 * snh0_773[k]
                    - f_5 * snh1_773[k]
                    + f_3 * pc_y[k] * sni_1031[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, pc_y, smi_808, smi_809, smi_810, snh0_774, \
                         snh0_775, snh0_776, snh1_774, snh1_775, snh1_776, sni_1032, sni_1033, \
                         sni_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = f_19 * smi_808[k]
                    + f_6 * snh0_774[k]
                    - f_7 * snh1_774[k]
                    + f_3 * pc_y[k] * sni_1032[k];

        t_1328[k] = f_19 * smi_809[k]
                    + f_8 * snh0_775[k]
                    - f_9 * snh1_775[k]
                    + f_3 * pc_y[k] * sni_1033[k];

        t_1329[k] = f_19 * smi_810[k]
                    + f_10 * snh0_776[k]
                    - f_11 * snh1_776[k]
                    + f_3 * pc_y[k] * sni_1034[k];
    }

#pragma omp simd aligned(t_1330, t_1331, t_1332, t_1333, pb_z, pc_y, pc_z, smk0_1008, smi_811, \
                         smi_812, smk1_1008, snh0_776, snh1_776, sni_1035, \
                         sni_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1330[k] = f_19 * smi_811[k]
                    + f_3 * pc_y[k] * sni_1035[k];

        t_1331[k] = f_1 * snh0_776[k]
                    - f_2 * snh1_776[k]
                    + f_3 * pc_z[k] * sni_1035[k];

        t_1332[k] = pb_z[k] * smk0_1008[k]
                    - f_12 * pc_z[k] * smk1_1008[k];

        t_1333[k] = f_20 * smi_812[k]
                    + f_3 * pc_y[k] * sni_1036[k];
    }

#pragma omp simd aligned(t_1334, t_1335, t_1336, pb_z, pc_y, pc_z, smk0_1011, smi_784, \
                         smi_814, smk1_1011, sni_1036, sni_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1334[k] = f_13 * smi_784[k]
                    + f_3 * pc_z[k] * sni_1036[k];

        t_1335[k] = pb_z[k] * smk0_1011[k]
                    - f_12 * pc_z[k] * smk1_1011[k];

        t_1336[k] = f_20 * smi_814[k]
                    + f_3 * pc_y[k] * sni_1038[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smk0,
                                                           const size_t smi, const size_t smk1,
                                                           const size_t snh0, const size_t snh1,
                                                           const size_t sni, const size_t ncols,
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
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smk0_1014 = buffer.data(smk0 + 1014);
    const auto *smk0_1018 = buffer.data(smk0 + 1018);
    const auto *smk0_1020 = buffer.data(smk0 + 1020);
    const auto *smk0_1023 = buffer.data(smk0 + 1023);
    const auto *smk0_1025 = buffer.data(smk0 + 1025);
    const auto *smk0_1026 = buffer.data(smk0 + 1026);
    const auto *smk0_1036 = buffer.data(smk0 + 1036);

    const auto *smi_787 = buffer.data(smi + 787);
    const auto *smi_790 = buffer.data(smi + 790);
    const auto *smi_791 = buffer.data(smi + 791);
    const auto *smi_794 = buffer.data(smi + 794);
    const auto *smi_795 = buffer.data(smi + 795);
    const auto *smi_796 = buffer.data(smi + 796);
    const auto *smi_805 = buffer.data(smi + 805);
    const auto *smi_811 = buffer.data(smi + 811);
    const auto *smi_812 = buffer.data(smi + 812);
    const auto *smi_815 = buffer.data(smi + 815);
    const auto *smi_817 = buffer.data(smi + 817);
    const auto *smi_818 = buffer.data(smi + 818);
    const auto *smi_821 = buffer.data(smi + 821);
    const auto *smi_822 = buffer.data(smi + 822);
    const auto *smi_826 = buffer.data(smi + 826);
    const auto *smi_833 = buffer.data(smi + 833);
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
    const auto *smi_854 = buffer.data(smi + 854);
    const auto *smi_861 = buffer.data(smi + 861);
    const auto *smi_863 = buffer.data(smi + 863);
    const auto *smi_864 = buffer.data(smi + 864);
    const auto *smi_865 = buffer.data(smi + 865);
    const auto *smi_866 = buffer.data(smi + 866);
    const auto *smi_867 = buffer.data(smi + 867);
    const auto *smi_868 = buffer.data(smi + 868);
    const auto *smi_870 = buffer.data(smi + 870);
    const auto *smi_873 = buffer.data(smi + 873);
    const auto *smi_877 = buffer.data(smi + 877);
    const auto *smi_882 = buffer.data(smi + 882);
    const auto *smi_889 = buffer.data(smi + 889);
    const auto *smi_891 = buffer.data(smi + 891);
    const auto *smi_892 = buffer.data(smi + 892);
    const auto *smi_893 = buffer.data(smi + 893);
    const auto *smi_894 = buffer.data(smi + 894);
    const auto *smi_895 = buffer.data(smi + 895);
    const auto *smi_1041 = buffer.data(smi + 1041);
    const auto *smi_1045 = buffer.data(smi + 1045);
    const auto *smi_1050 = buffer.data(smi + 1050);
    const auto *smi_1056 = buffer.data(smi + 1056);
    const auto *smi_1057 = buffer.data(smi + 1057);
    const auto *smi_1058 = buffer.data(smi + 1058);
    const auto *smi_1059 = buffer.data(smi + 1059);
    const auto *smi_1060 = buffer.data(smi + 1060);
    const auto *smi_1061 = buffer.data(smi + 1061);
    const auto *smi_1062 = buffer.data(smi + 1062);
    const auto *smi_1063 = buffer.data(smi + 1063);
    const auto *smi_1064 = buffer.data(smi + 1064);
    const auto *smi_1067 = buffer.data(smi + 1067);
    const auto *smi_1069 = buffer.data(smi + 1069);
    const auto *smi_1070 = buffer.data(smi + 1070);
    const auto *smi_1073 = buffer.data(smi + 1073);
    const auto *smi_1074 = buffer.data(smi + 1074);
    const auto *smi_1076 = buffer.data(smi + 1076);
    const auto *smi_1078 = buffer.data(smi + 1078);
    const auto *smi_1079 = buffer.data(smi + 1079);
    const auto *smi_1081 = buffer.data(smi + 1081);
    const auto *smi_1082 = buffer.data(smi + 1082);
    const auto *smi_1084 = buffer.data(smi + 1084);
    const auto *smi_1085 = buffer.data(smi + 1085);
    const auto *smi_1086 = buffer.data(smi + 1086);
    const auto *smi_1087 = buffer.data(smi + 1087);
    const auto *smi_1088 = buffer.data(smi + 1088);
    const auto *smi_1089 = buffer.data(smi + 1089);
    const auto *smi_1090 = buffer.data(smi + 1090);
    const auto *smi_1091 = buffer.data(smi + 1091);
    const auto *smi_1092 = buffer.data(smi + 1092);
    const auto *smi_1095 = buffer.data(smi + 1095);
    const auto *smi_1097 = buffer.data(smi + 1097);
    const auto *smi_1098 = buffer.data(smi + 1098);
    const auto *smi_1101 = buffer.data(smi + 1101);
    const auto *smi_1102 = buffer.data(smi + 1102);
    const auto *smi_1104 = buffer.data(smi + 1104);
    const auto *smi_1106 = buffer.data(smi + 1106);
    const auto *smi_1107 = buffer.data(smi + 1107);
    const auto *smi_1109 = buffer.data(smi + 1109);
    const auto *smi_1110 = buffer.data(smi + 1110);
    const auto *smi_1112 = buffer.data(smi + 1112);
    const auto *smi_1113 = buffer.data(smi + 1113);
    const auto *smi_1114 = buffer.data(smi + 1114);
    const auto *smi_1115 = buffer.data(smi + 1115);
    const auto *smi_1116 = buffer.data(smi + 1116);
    const auto *smi_1117 = buffer.data(smi + 1117);
    const auto *smi_1118 = buffer.data(smi + 1118);
    const auto *smi_1119 = buffer.data(smi + 1119);
    const auto *smi_1120 = buffer.data(smi + 1120);

    const auto *smk1_1014 = buffer.data(smk1 + 1014);
    const auto *smk1_1018 = buffer.data(smk1 + 1018);
    const auto *smk1_1020 = buffer.data(smk1 + 1020);
    const auto *smk1_1023 = buffer.data(smk1 + 1023);
    const auto *smk1_1025 = buffer.data(smk1 + 1025);
    const auto *smk1_1026 = buffer.data(smk1 + 1026);
    const auto *smk1_1036 = buffer.data(smk1 + 1036);

    const auto *snh0_782 = buffer.data(snh0 + 782);
    const auto *snh0_786 = buffer.data(snh0 + 786);
    const auto *snh0_791 = buffer.data(snh0 + 791);
    const auto *snh0_794 = buffer.data(snh0 + 794);
    const auto *snh0_795 = buffer.data(snh0 + 795);
    const auto *snh0_796 = buffer.data(snh0 + 796);
    const auto *snh0_797 = buffer.data(snh0 + 797);
    const auto *snh0_798 = buffer.data(snh0 + 798);
    const auto *snh0_801 = buffer.data(snh0 + 801);
    const auto *snh0_803 = buffer.data(snh0 + 803);
    const auto *snh0_804 = buffer.data(snh0 + 804);
    const auto *snh0_807 = buffer.data(snh0 + 807);
    const auto *snh0_808 = buffer.data(snh0 + 808);
    const auto *snh0_810 = buffer.data(snh0 + 810);
    const auto *snh0_812 = buffer.data(snh0 + 812);
    const auto *snh0_813 = buffer.data(snh0 + 813);
    const auto *snh0_815 = buffer.data(snh0 + 815);
    const auto *snh0_816 = buffer.data(snh0 + 816);
    const auto *snh0_817 = buffer.data(snh0 + 817);
    const auto *snh0_818 = buffer.data(snh0 + 818);
    const auto *snh0_819 = buffer.data(snh0 + 819);
    const auto *snh0_822 = buffer.data(snh0 + 822);
    const auto *snh0_824 = buffer.data(snh0 + 824);
    const auto *snh0_825 = buffer.data(snh0 + 825);
    const auto *snh0_828 = buffer.data(snh0 + 828);
    const auto *snh0_829 = buffer.data(snh0 + 829);
    const auto *snh0_831 = buffer.data(snh0 + 831);
    const auto *snh0_833 = buffer.data(snh0 + 833);
    const auto *snh0_834 = buffer.data(snh0 + 834);
    const auto *snh0_836 = buffer.data(snh0 + 836);
    const auto *snh0_837 = buffer.data(snh0 + 837);
    const auto *snh0_838 = buffer.data(snh0 + 838);
    const auto *snh0_839 = buffer.data(snh0 + 839);
    const auto *snh0_840 = buffer.data(snh0 + 840);

    const auto *snh1_782 = buffer.data(snh1 + 782);
    const auto *snh1_786 = buffer.data(snh1 + 786);
    const auto *snh1_791 = buffer.data(snh1 + 791);
    const auto *snh1_794 = buffer.data(snh1 + 794);
    const auto *snh1_795 = buffer.data(snh1 + 795);
    const auto *snh1_796 = buffer.data(snh1 + 796);
    const auto *snh1_797 = buffer.data(snh1 + 797);
    const auto *snh1_798 = buffer.data(snh1 + 798);
    const auto *snh1_801 = buffer.data(snh1 + 801);
    const auto *snh1_803 = buffer.data(snh1 + 803);
    const auto *snh1_804 = buffer.data(snh1 + 804);
    const auto *snh1_807 = buffer.data(snh1 + 807);
    const auto *snh1_808 = buffer.data(snh1 + 808);
    const auto *snh1_810 = buffer.data(snh1 + 810);
    const auto *snh1_812 = buffer.data(snh1 + 812);
    const auto *snh1_813 = buffer.data(snh1 + 813);
    const auto *snh1_815 = buffer.data(snh1 + 815);
    const auto *snh1_816 = buffer.data(snh1 + 816);
    const auto *snh1_817 = buffer.data(snh1 + 817);
    const auto *snh1_818 = buffer.data(snh1 + 818);
    const auto *snh1_819 = buffer.data(snh1 + 819);
    const auto *snh1_822 = buffer.data(snh1 + 822);
    const auto *snh1_824 = buffer.data(snh1 + 824);
    const auto *snh1_825 = buffer.data(snh1 + 825);
    const auto *snh1_828 = buffer.data(snh1 + 828);
    const auto *snh1_829 = buffer.data(snh1 + 829);
    const auto *snh1_831 = buffer.data(snh1 + 831);
    const auto *snh1_833 = buffer.data(snh1 + 833);
    const auto *snh1_834 = buffer.data(snh1 + 834);
    const auto *snh1_836 = buffer.data(snh1 + 836);
    const auto *snh1_837 = buffer.data(snh1 + 837);
    const auto *snh1_838 = buffer.data(snh1 + 838);
    const auto *snh1_839 = buffer.data(snh1 + 839);
    const auto *snh1_840 = buffer.data(snh1 + 840);

    const auto *sni_1039 = buffer.data(sni + 1039);
    const auto *sni_1041 = buffer.data(sni + 1041);
    const auto *sni_1042 = buffer.data(sni + 1042);
    const auto *sni_1045 = buffer.data(sni + 1045);
    const auto *sni_1046 = buffer.data(sni + 1046);
    const auto *sni_1050 = buffer.data(sni + 1050);
    const auto *sni_1056 = buffer.data(sni + 1056);
    const auto *sni_1057 = buffer.data(sni + 1057);
    const auto *sni_1058 = buffer.data(sni + 1058);
    const auto *sni_1059 = buffer.data(sni + 1059);
    const auto *sni_1060 = buffer.data(sni + 1060);
    const auto *sni_1061 = buffer.data(sni + 1061);
    const auto *sni_1062 = buffer.data(sni + 1062);
    const auto *sni_1063 = buffer.data(sni + 1063);
    const auto *sni_1064 = buffer.data(sni + 1064);
    const auto *sni_1066 = buffer.data(sni + 1066);
    const auto *sni_1067 = buffer.data(sni + 1067);
    const auto *sni_1069 = buffer.data(sni + 1069);
    const auto *sni_1070 = buffer.data(sni + 1070);
    const auto *sni_1073 = buffer.data(sni + 1073);
    const auto *sni_1074 = buffer.data(sni + 1074);
    const auto *sni_1076 = buffer.data(sni + 1076);
    const auto *sni_1078 = buffer.data(sni + 1078);
    const auto *sni_1079 = buffer.data(sni + 1079);
    const auto *sni_1081 = buffer.data(sni + 1081);
    const auto *sni_1082 = buffer.data(sni + 1082);
    const auto *sni_1084 = buffer.data(sni + 1084);
    const auto *sni_1085 = buffer.data(sni + 1085);
    const auto *sni_1086 = buffer.data(sni + 1086);
    const auto *sni_1087 = buffer.data(sni + 1087);
    const auto *sni_1088 = buffer.data(sni + 1088);
    const auto *sni_1089 = buffer.data(sni + 1089);
    const auto *sni_1090 = buffer.data(sni + 1090);
    const auto *sni_1091 = buffer.data(sni + 1091);
    const auto *sni_1092 = buffer.data(sni + 1092);
    const auto *sni_1094 = buffer.data(sni + 1094);
    const auto *sni_1095 = buffer.data(sni + 1095);
    const auto *sni_1097 = buffer.data(sni + 1097);
    const auto *sni_1098 = buffer.data(sni + 1098);
    const auto *sni_1101 = buffer.data(sni + 1101);
    const auto *sni_1102 = buffer.data(sni + 1102);
    const auto *sni_1104 = buffer.data(sni + 1104);
    const auto *sni_1106 = buffer.data(sni + 1106);
    const auto *sni_1107 = buffer.data(sni + 1107);
    const auto *sni_1109 = buffer.data(sni + 1109);
    const auto *sni_1110 = buffer.data(sni + 1110);
    const auto *sni_1112 = buffer.data(sni + 1112);
    const auto *sni_1113 = buffer.data(sni + 1113);
    const auto *sni_1114 = buffer.data(sni + 1114);
    const auto *sni_1115 = buffer.data(sni + 1115);
    const auto *sni_1116 = buffer.data(sni + 1116);
    const auto *sni_1117 = buffer.data(sni + 1117);
    const auto *sni_1118 = buffer.data(sni + 1118);
    const auto *sni_1119 = buffer.data(sni + 1119);
    const auto *sni_1120 = buffer.data(sni + 1120);

#pragma omp simd aligned(t_1337, t_1338, t_1339, pb_z, pc_x, pc_z, smk0_1014, smi_787, \
                         smi_1041, smk1_1014, snh0_782, snh1_782, sni_1039, \
                         sni_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1337[k] = f_14 * smi_1041[k]
                    + f_4 * snh0_782[k]
                    - f_5 * snh1_782[k]
                    + f_3 * pc_x[k] * sni_1041[k];

        t_1338[k] = pb_z[k] * smk0_1014[k]
                    - f_12 * pc_z[k] * smk1_1014[k];

        t_1339[k] = f_13 * smi_787[k]
                    + f_3 * pc_z[k] * sni_1039[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, pb_z, pc_x, pc_y, pc_z, smk0_1018, smi_817, \
                         smi_1045, smk1_1018, snh0_786, snh1_786, sni_1041, \
                         sni_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_20 * smi_817[k]
                    + f_3 * pc_y[k] * sni_1041[k];

        t_1341[k] = f_14 * smi_1045[k]
                    + f_6 * snh0_786[k]
                    - f_7 * snh1_786[k]
                    + f_3 * pc_x[k] * sni_1045[k];

        t_1342[k] = pb_z[k] * smk0_1018[k]
                    - f_12 * pc_z[k] * smk1_1018[k];
    }

#pragma omp simd aligned(t_1343, t_1344, t_1345, pb_z, pc_y, pc_z, smk0_1020, smi_790, \
                         smi_791, smi_821, smk1_1020, sni_1042, \
                         sni_1045 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1343[k] = f_13 * smi_790[k]
                    + f_3 * pc_z[k] * sni_1042[k];

        t_1344[k] = pb_z[k] * smk0_1020[k]
                    + f_14 * smi_791[k]
                    - f_12 * pc_z[k] * smk1_1020[k];

        t_1345[k] = f_20 * smi_821[k]
                    + f_3 * pc_y[k] * sni_1045[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pb_z, pc_x, pc_z, smk0_1023, smi_794, \
                         smi_1050, smk1_1023, snh0_791, snh1_791, sni_1046, \
                         sni_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_14 * smi_1050[k]
                    + f_8 * snh0_791[k]
                    - f_9 * snh1_791[k]
                    + f_3 * pc_x[k] * sni_1050[k];

        t_1347[k] = pb_z[k] * smk0_1023[k]
                    - f_12 * pc_z[k] * smk1_1023[k];

        t_1348[k] = f_13 * smi_794[k]
                    + f_3 * pc_z[k] * sni_1046[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pb_z, pc_y, pc_z, smk0_1025, smk0_1026, \
                         smi_795, smi_796, smi_826, smk1_1025, smk1_1026, \
                         sni_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = pb_z[k] * smk0_1025[k]
                    + f_14 * smi_795[k]
                    - f_12 * pc_z[k] * smk1_1025[k];

        t_1350[k] = pb_z[k] * smk0_1026[k]
                    + f_15 * smi_796[k]
                    - f_12 * pc_z[k] * smk1_1026[k];

        t_1351[k] = f_20 * smi_826[k]
                    + f_3 * pc_y[k] * sni_1050[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, t_1355, pc_x, smi_1056, smi_1057, smi_1058, \
                         smi_1059, snh0_797, snh1_797, sni_1056, sni_1057, sni_1058, \
                         sni_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_14 * smi_1056[k]
                    + f_10 * snh0_797[k]
                    - f_11 * snh1_797[k]
                    + f_3 * pc_x[k] * sni_1056[k];

        t_1353[k] = f_14 * smi_1057[k]
                    + f_3 * pc_x[k] * sni_1057[k];

        t_1354[k] = f_14 * smi_1058[k]
                    + f_3 * pc_x[k] * sni_1058[k];

        t_1355[k] = f_14 * smi_1059[k]
                    + f_3 * pc_x[k] * sni_1059[k];
    }

#pragma omp simd aligned(t_1356, t_1357, t_1358, t_1359, pc_x, smi_1060, smi_1061, smi_1062, \
                         smi_1063, sni_1060, sni_1061, sni_1062, \
                         sni_1063 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1356[k] = f_14 * smi_1060[k]
                    + f_3 * pc_x[k] * sni_1060[k];

        t_1357[k] = f_14 * smi_1061[k]
                    + f_3 * pc_x[k] * sni_1061[k];

        t_1358[k] = f_14 * smi_1062[k]
                    + f_3 * pc_x[k] * sni_1062[k];

        t_1359[k] = f_14 * smi_1063[k]
                    + f_3 * pc_x[k] * sni_1063[k];
    }

#pragma omp simd aligned(t_1360, t_1361, t_1362, pb_z, pc_y, pc_z, smk0_1036, smi_805, \
                         smi_835, smk1_1036, snh0_794, snh1_794, sni_1057, \
                         sni_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1360[k] = pb_z[k] * smk0_1036[k]
                    - f_12 * pc_z[k] * smk1_1036[k];

        t_1361[k] = f_13 * smi_805[k]
                    + f_3 * pc_z[k] * sni_1057[k];

        t_1362[k] = f_20 * smi_835[k]
                    + f_4 * snh0_794[k]
                    - f_5 * snh1_794[k]
                    + f_3 * pc_y[k] * sni_1059[k];
    }

#pragma omp simd aligned(t_1363, t_1364, t_1365, pc_y, smi_836, smi_837, smi_838, snh0_795, \
                         snh0_796, snh0_797, snh1_795, snh1_796, snh1_797, sni_1060, sni_1061, \
                         sni_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1363[k] = f_20 * smi_836[k]
                    + f_6 * snh0_795[k]
                    - f_7 * snh1_795[k]
                    + f_3 * pc_y[k] * sni_1060[k];

        t_1364[k] = f_20 * smi_837[k]
                    + f_8 * snh0_796[k]
                    - f_9 * snh1_796[k]
                    + f_3 * pc_y[k] * sni_1061[k];

        t_1365[k] = f_20 * smi_838[k]
                    + f_10 * snh0_797[k]
                    - f_11 * snh1_797[k]
                    + f_3 * pc_y[k] * sni_1062[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, pc_x, pc_y, pc_z, smi_811, smi_839, smi_1064, \
                         snh0_797, snh0_798, snh1_797, snh1_798, sni_1063, \
                         sni_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_20 * smi_839[k]
                    + f_3 * pc_y[k] * sni_1063[k];

        t_1367[k] = f_13 * smi_811[k]
                    + f_1 * snh0_797[k]
                    - f_2 * snh1_797[k]
                    + f_3 * pc_z[k] * sni_1063[k];

        t_1368[k] = f_14 * smi_1064[k]
                    + f_1 * snh0_798[k]
                    - f_2 * snh1_798[k]
                    + f_3 * pc_x[k] * sni_1064[k];
    }

#pragma omp simd aligned(t_1369, t_1370, t_1371, t_1372, pc_x, pc_y, pc_z, smi_812, smi_840, \
                         smi_842, smi_1067, snh0_801, snh1_801, sni_1064, sni_1066, \
                         sni_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1369[k] = f_21 * smi_840[k]
                    + f_3 * pc_y[k] * sni_1064[k];

        t_1370[k] = f_14 * smi_812[k]
                    + f_3 * pc_z[k] * sni_1064[k];

        t_1371[k] = f_14 * smi_1067[k]
                    + f_4 * snh0_801[k]
                    - f_5 * snh1_801[k]
                    + f_3 * pc_x[k] * sni_1067[k];

        t_1372[k] = f_21 * smi_842[k]
                    + f_3 * pc_y[k] * sni_1066[k];
    }

#pragma omp simd aligned(t_1373, t_1374, t_1375, pc_x, pc_z, smi_815, smi_1069, smi_1070, \
                         snh0_803, snh0_804, snh1_803, snh1_804, sni_1067, sni_1069, \
                         sni_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1373[k] = f_14 * smi_1069[k]
                    + f_4 * snh0_803[k]
                    - f_5 * snh1_803[k]
                    + f_3 * pc_x[k] * sni_1069[k];

        t_1374[k] = f_14 * smi_1070[k]
                    + f_6 * snh0_804[k]
                    - f_7 * snh1_804[k]
                    + f_3 * pc_x[k] * sni_1070[k];

        t_1375[k] = f_14 * smi_815[k]
                    + f_3 * pc_z[k] * sni_1067[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, pc_x, pc_y, smi_845, smi_1073, smi_1074, \
                         snh0_807, snh0_808, snh1_807, snh1_808, sni_1069, sni_1073, \
                         sni_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_21 * smi_845[k]
                    + f_3 * pc_y[k] * sni_1069[k];

        t_1377[k] = f_14 * smi_1073[k]
                    + f_6 * snh0_807[k]
                    - f_7 * snh1_807[k]
                    + f_3 * pc_x[k] * sni_1073[k];

        t_1378[k] = f_14 * smi_1074[k]
                    + f_8 * snh0_808[k]
                    - f_9 * snh1_808[k]
                    + f_3 * pc_x[k] * sni_1074[k];
    }

#pragma omp simd aligned(t_1379, t_1380, t_1381, pc_x, pc_y, pc_z, smi_818, smi_849, smi_1076, \
                         snh0_810, snh1_810, sni_1070, sni_1073, \
                         sni_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1379[k] = f_14 * smi_818[k]
                    + f_3 * pc_z[k] * sni_1070[k];

        t_1380[k] = f_14 * smi_1076[k]
                    + f_8 * snh0_810[k]
                    - f_9 * snh1_810[k]
                    + f_3 * pc_x[k] * sni_1076[k];

        t_1381[k] = f_21 * smi_849[k]
                    + f_3 * pc_y[k] * sni_1073[k];
    }

#pragma omp simd aligned(t_1382, t_1383, t_1384, pc_x, pc_z, smi_822, smi_1078, smi_1079, \
                         snh0_812, snh0_813, snh1_812, snh1_813, sni_1074, sni_1078, \
                         sni_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1382[k] = f_14 * smi_1078[k]
                    + f_8 * snh0_812[k]
                    - f_9 * snh1_812[k]
                    + f_3 * pc_x[k] * sni_1078[k];

        t_1383[k] = f_14 * smi_1079[k]
                    + f_10 * snh0_813[k]
                    - f_11 * snh1_813[k]
                    + f_3 * pc_x[k] * sni_1079[k];

        t_1384[k] = f_14 * smi_822[k]
                    + f_3 * pc_z[k] * sni_1074[k];
    }

#pragma omp simd aligned(t_1385, t_1386, t_1387, pc_x, pc_y, smi_854, smi_1081, smi_1082, \
                         snh0_815, snh0_816, snh1_815, snh1_816, sni_1078, sni_1081, \
                         sni_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1385[k] = f_14 * smi_1081[k]
                    + f_10 * snh0_815[k]
                    - f_11 * snh1_815[k]
                    + f_3 * pc_x[k] * sni_1081[k];

        t_1386[k] = f_14 * smi_1082[k]
                    + f_10 * snh0_816[k]
                    - f_11 * snh1_816[k]
                    + f_3 * pc_x[k] * sni_1082[k];

        t_1387[k] = f_21 * smi_854[k]
                    + f_3 * pc_y[k] * sni_1078[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, t_1391, pc_x, smi_1084, smi_1085, smi_1086, \
                         smi_1087, snh0_818, snh1_818, sni_1084, sni_1085, sni_1086, \
                         sni_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = f_14 * smi_1084[k]
                    + f_10 * snh0_818[k]
                    - f_11 * snh1_818[k]
                    + f_3 * pc_x[k] * sni_1084[k];

        t_1389[k] = f_14 * smi_1085[k]
                    + f_3 * pc_x[k] * sni_1085[k];

        t_1390[k] = f_14 * smi_1086[k]
                    + f_3 * pc_x[k] * sni_1086[k];

        t_1391[k] = f_14 * smi_1087[k]
                    + f_3 * pc_x[k] * sni_1087[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, t_1395, pc_x, smi_1088, smi_1089, smi_1090, \
                         smi_1091, sni_1088, sni_1089, sni_1090, \
                         sni_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_14 * smi_1088[k]
                    + f_3 * pc_x[k] * sni_1088[k];

        t_1393[k] = f_14 * smi_1089[k]
                    + f_3 * pc_x[k] * sni_1089[k];

        t_1394[k] = f_14 * smi_1090[k]
                    + f_3 * pc_x[k] * sni_1090[k];

        t_1395[k] = f_14 * smi_1091[k]
                    + f_3 * pc_x[k] * sni_1091[k];
    }

#pragma omp simd aligned(t_1396, t_1397, t_1398, pc_y, pc_z, smi_833, smi_861, smi_863, \
                         snh0_813, snh0_815, snh1_813, snh1_815, sni_1085, \
                         sni_1087 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1396[k] = f_21 * smi_861[k]
                    + f_1 * snh0_813[k]
                    - f_2 * snh1_813[k]
                    + f_3 * pc_y[k] * sni_1085[k];

        t_1397[k] = f_14 * smi_833[k]
                    + f_3 * pc_z[k] * sni_1085[k];

        t_1398[k] = f_21 * smi_863[k]
                    + f_4 * snh0_815[k]
                    - f_5 * snh1_815[k]
                    + f_3 * pc_y[k] * sni_1087[k];
    }

#pragma omp simd aligned(t_1399, t_1400, t_1401, pc_y, smi_864, smi_865, smi_866, snh0_816, \
                         snh0_817, snh0_818, snh1_816, snh1_817, snh1_818, sni_1088, sni_1089, \
                         sni_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1399[k] = f_21 * smi_864[k]
                    + f_6 * snh0_816[k]
                    - f_7 * snh1_816[k]
                    + f_3 * pc_y[k] * sni_1088[k];

        t_1400[k] = f_21 * smi_865[k]
                    + f_8 * snh0_817[k]
                    - f_9 * snh1_817[k]
                    + f_3 * pc_y[k] * sni_1089[k];

        t_1401[k] = f_21 * smi_866[k]
                    + f_10 * snh0_818[k]
                    - f_11 * snh1_818[k]
                    + f_3 * pc_y[k] * sni_1090[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, pc_x, pc_y, pc_z, smi_839, smi_867, smi_1092, \
                         snh0_818, snh0_819, snh1_818, snh1_819, sni_1091, \
                         sni_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_21 * smi_867[k]
                    + f_3 * pc_y[k] * sni_1091[k];

        t_1403[k] = f_14 * smi_839[k]
                    + f_1 * snh0_818[k]
                    - f_2 * snh1_818[k]
                    + f_3 * pc_z[k] * sni_1091[k];

        t_1404[k] = f_14 * smi_1092[k]
                    + f_1 * snh0_819[k]
                    - f_2 * snh1_819[k]
                    + f_3 * pc_x[k] * sni_1092[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, t_1408, pc_x, pc_y, pc_z, smi_840, smi_868, \
                         smi_870, smi_1095, snh0_822, snh1_822, sni_1092, sni_1094, \
                         sni_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = f_17 * smi_868[k]
                    + f_3 * pc_y[k] * sni_1092[k];

        t_1406[k] = f_15 * smi_840[k]
                    + f_3 * pc_z[k] * sni_1092[k];

        t_1407[k] = f_14 * smi_1095[k]
                    + f_4 * snh0_822[k]
                    - f_5 * snh1_822[k]
                    + f_3 * pc_x[k] * sni_1095[k];

        t_1408[k] = f_17 * smi_870[k]
                    + f_3 * pc_y[k] * sni_1094[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pc_x, pc_z, smi_843, smi_1097, smi_1098, \
                         snh0_824, snh0_825, snh1_824, snh1_825, sni_1095, sni_1097, \
                         sni_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_14 * smi_1097[k]
                    + f_4 * snh0_824[k]
                    - f_5 * snh1_824[k]
                    + f_3 * pc_x[k] * sni_1097[k];

        t_1410[k] = f_14 * smi_1098[k]
                    + f_6 * snh0_825[k]
                    - f_7 * snh1_825[k]
                    + f_3 * pc_x[k] * sni_1098[k];

        t_1411[k] = f_15 * smi_843[k]
                    + f_3 * pc_z[k] * sni_1095[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, pc_x, pc_y, smi_873, smi_1101, smi_1102, \
                         snh0_828, snh0_829, snh1_828, snh1_829, sni_1097, sni_1101, \
                         sni_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_17 * smi_873[k]
                    + f_3 * pc_y[k] * sni_1097[k];

        t_1413[k] = f_14 * smi_1101[k]
                    + f_6 * snh0_828[k]
                    - f_7 * snh1_828[k]
                    + f_3 * pc_x[k] * sni_1101[k];

        t_1414[k] = f_14 * smi_1102[k]
                    + f_8 * snh0_829[k]
                    - f_9 * snh1_829[k]
                    + f_3 * pc_x[k] * sni_1102[k];
    }

#pragma omp simd aligned(t_1415, t_1416, t_1417, pc_x, pc_y, pc_z, smi_846, smi_877, smi_1104, \
                         snh0_831, snh1_831, sni_1098, sni_1101, \
                         sni_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1415[k] = f_15 * smi_846[k]
                    + f_3 * pc_z[k] * sni_1098[k];

        t_1416[k] = f_14 * smi_1104[k]
                    + f_8 * snh0_831[k]
                    - f_9 * snh1_831[k]
                    + f_3 * pc_x[k] * sni_1104[k];

        t_1417[k] = f_17 * smi_877[k]
                    + f_3 * pc_y[k] * sni_1101[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, pc_x, pc_z, smi_850, smi_1106, smi_1107, \
                         snh0_833, snh0_834, snh1_833, snh1_834, sni_1102, sni_1106, \
                         sni_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_14 * smi_1106[k]
                    + f_8 * snh0_833[k]
                    - f_9 * snh1_833[k]
                    + f_3 * pc_x[k] * sni_1106[k];

        t_1419[k] = f_14 * smi_1107[k]
                    + f_10 * snh0_834[k]
                    - f_11 * snh1_834[k]
                    + f_3 * pc_x[k] * sni_1107[k];

        t_1420[k] = f_15 * smi_850[k]
                    + f_3 * pc_z[k] * sni_1102[k];
    }

#pragma omp simd aligned(t_1421, t_1422, t_1423, pc_x, pc_y, smi_882, smi_1109, smi_1110, \
                         snh0_836, snh0_837, snh1_836, snh1_837, sni_1106, sni_1109, \
                         sni_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1421[k] = f_14 * smi_1109[k]
                    + f_10 * snh0_836[k]
                    - f_11 * snh1_836[k]
                    + f_3 * pc_x[k] * sni_1109[k];

        t_1422[k] = f_14 * smi_1110[k]
                    + f_10 * snh0_837[k]
                    - f_11 * snh1_837[k]
                    + f_3 * pc_x[k] * sni_1110[k];

        t_1423[k] = f_17 * smi_882[k]
                    + f_3 * pc_y[k] * sni_1106[k];
    }

#pragma omp simd aligned(t_1424, t_1425, t_1426, t_1427, pc_x, smi_1112, smi_1113, smi_1114, \
                         smi_1115, snh0_839, snh1_839, sni_1112, sni_1113, sni_1114, \
                         sni_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1424[k] = f_14 * smi_1112[k]
                    + f_10 * snh0_839[k]
                    - f_11 * snh1_839[k]
                    + f_3 * pc_x[k] * sni_1112[k];

        t_1425[k] = f_14 * smi_1113[k]
                    + f_3 * pc_x[k] * sni_1113[k];

        t_1426[k] = f_14 * smi_1114[k]
                    + f_3 * pc_x[k] * sni_1114[k];

        t_1427[k] = f_14 * smi_1115[k]
                    + f_3 * pc_x[k] * sni_1115[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, t_1431, pc_x, smi_1116, smi_1117, smi_1118, \
                         smi_1119, sni_1116, sni_1117, sni_1118, \
                         sni_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_14 * smi_1116[k]
                    + f_3 * pc_x[k] * sni_1116[k];

        t_1429[k] = f_14 * smi_1117[k]
                    + f_3 * pc_x[k] * sni_1117[k];

        t_1430[k] = f_14 * smi_1118[k]
                    + f_3 * pc_x[k] * sni_1118[k];

        t_1431[k] = f_14 * smi_1119[k]
                    + f_3 * pc_x[k] * sni_1119[k];
    }

#pragma omp simd aligned(t_1432, t_1433, t_1434, pc_y, pc_z, smi_861, smi_889, smi_891, \
                         snh0_834, snh0_836, snh1_834, snh1_836, sni_1113, \
                         sni_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1432[k] = f_17 * smi_889[k]
                    + f_1 * snh0_834[k]
                    - f_2 * snh1_834[k]
                    + f_3 * pc_y[k] * sni_1113[k];

        t_1433[k] = f_15 * smi_861[k]
                    + f_3 * pc_z[k] * sni_1113[k];

        t_1434[k] = f_17 * smi_891[k]
                    + f_4 * snh0_836[k]
                    - f_5 * snh1_836[k]
                    + f_3 * pc_y[k] * sni_1115[k];
    }

#pragma omp simd aligned(t_1435, t_1436, t_1437, pc_y, smi_892, smi_893, smi_894, snh0_837, \
                         snh0_838, snh0_839, snh1_837, snh1_838, snh1_839, sni_1116, sni_1117, \
                         sni_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1435[k] = f_17 * smi_892[k]
                    + f_6 * snh0_837[k]
                    - f_7 * snh1_837[k]
                    + f_3 * pc_y[k] * sni_1116[k];

        t_1436[k] = f_17 * smi_893[k]
                    + f_8 * snh0_838[k]
                    - f_9 * snh1_838[k]
                    + f_3 * pc_y[k] * sni_1117[k];

        t_1437[k] = f_17 * smi_894[k]
                    + f_10 * snh0_839[k]
                    - f_11 * snh1_839[k]
                    + f_3 * pc_y[k] * sni_1118[k];
    }

#pragma omp simd aligned(t_1438, t_1439, t_1440, pc_x, pc_y, pc_z, smi_867, smi_895, smi_1120, \
                         snh0_839, snh0_840, snh1_839, snh1_840, sni_1119, \
                         sni_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1438[k] = f_17 * smi_895[k]
                    + f_3 * pc_y[k] * sni_1119[k];

        t_1439[k] = f_15 * smi_867[k]
                    + f_1 * snh0_839[k]
                    - f_2 * snh1_839[k]
                    + f_3 * pc_z[k] * sni_1119[k];

        t_1440[k] = f_14 * smi_1120[k]
                    + f_1 * snh0_840[k]
                    - f_2 * snh1_840[k]
                    + f_3 * pc_x[k] * sni_1120[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t smi, const size_t snh0,
                                                           const size_t snh1, const size_t sni,
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
    const auto f_21 = 3.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi_868 = buffer.data(smi + 868);
    const auto *smi_871 = buffer.data(smi + 871);
    const auto *smi_874 = buffer.data(smi + 874);
    const auto *smi_878 = buffer.data(smi + 878);
    const auto *smi_889 = buffer.data(smi + 889);
    const auto *smi_895 = buffer.data(smi + 895);
    const auto *smi_896 = buffer.data(smi + 896);
    const auto *smi_898 = buffer.data(smi + 898);
    const auto *smi_899 = buffer.data(smi + 899);
    const auto *smi_901 = buffer.data(smi + 901);
    const auto *smi_902 = buffer.data(smi + 902);
    const auto *smi_905 = buffer.data(smi + 905);
    const auto *smi_906 = buffer.data(smi + 906);
    const auto *smi_910 = buffer.data(smi + 910);
    const auto *smi_917 = buffer.data(smi + 917);
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
    const auto *smi_938 = buffer.data(smi + 938);
    const auto *smi_945 = buffer.data(smi + 945);
    const auto *smi_947 = buffer.data(smi + 947);
    const auto *smi_948 = buffer.data(smi + 948);
    const auto *smi_949 = buffer.data(smi + 949);
    const auto *smi_950 = buffer.data(smi + 950);
    const auto *smi_951 = buffer.data(smi + 951);
    const auto *smi_952 = buffer.data(smi + 952);
    const auto *smi_954 = buffer.data(smi + 954);
    const auto *smi_957 = buffer.data(smi + 957);
    const auto *smi_961 = buffer.data(smi + 961);
    const auto *smi_966 = buffer.data(smi + 966);
    const auto *smi_973 = buffer.data(smi + 973);
    const auto *smi_975 = buffer.data(smi + 975);
    const auto *smi_976 = buffer.data(smi + 976);
    const auto *smi_977 = buffer.data(smi + 977);
    const auto *smi_978 = buffer.data(smi + 978);
    const auto *smi_1123 = buffer.data(smi + 1123);
    const auto *smi_1125 = buffer.data(smi + 1125);
    const auto *smi_1126 = buffer.data(smi + 1126);
    const auto *smi_1129 = buffer.data(smi + 1129);
    const auto *smi_1130 = buffer.data(smi + 1130);
    const auto *smi_1132 = buffer.data(smi + 1132);
    const auto *smi_1134 = buffer.data(smi + 1134);
    const auto *smi_1135 = buffer.data(smi + 1135);
    const auto *smi_1137 = buffer.data(smi + 1137);
    const auto *smi_1138 = buffer.data(smi + 1138);
    const auto *smi_1140 = buffer.data(smi + 1140);
    const auto *smi_1141 = buffer.data(smi + 1141);
    const auto *smi_1142 = buffer.data(smi + 1142);
    const auto *smi_1143 = buffer.data(smi + 1143);
    const auto *smi_1144 = buffer.data(smi + 1144);
    const auto *smi_1145 = buffer.data(smi + 1145);
    const auto *smi_1146 = buffer.data(smi + 1146);
    const auto *smi_1147 = buffer.data(smi + 1147);
    const auto *smi_1148 = buffer.data(smi + 1148);
    const auto *smi_1151 = buffer.data(smi + 1151);
    const auto *smi_1153 = buffer.data(smi + 1153);
    const auto *smi_1154 = buffer.data(smi + 1154);
    const auto *smi_1157 = buffer.data(smi + 1157);
    const auto *smi_1158 = buffer.data(smi + 1158);
    const auto *smi_1160 = buffer.data(smi + 1160);
    const auto *smi_1162 = buffer.data(smi + 1162);
    const auto *smi_1163 = buffer.data(smi + 1163);
    const auto *smi_1165 = buffer.data(smi + 1165);
    const auto *smi_1166 = buffer.data(smi + 1166);
    const auto *smi_1168 = buffer.data(smi + 1168);
    const auto *smi_1169 = buffer.data(smi + 1169);
    const auto *smi_1170 = buffer.data(smi + 1170);
    const auto *smi_1171 = buffer.data(smi + 1171);
    const auto *smi_1172 = buffer.data(smi + 1172);
    const auto *smi_1173 = buffer.data(smi + 1173);
    const auto *smi_1174 = buffer.data(smi + 1174);
    const auto *smi_1175 = buffer.data(smi + 1175);
    const auto *smi_1176 = buffer.data(smi + 1176);
    const auto *smi_1179 = buffer.data(smi + 1179);
    const auto *smi_1181 = buffer.data(smi + 1181);
    const auto *smi_1182 = buffer.data(smi + 1182);
    const auto *smi_1185 = buffer.data(smi + 1185);
    const auto *smi_1186 = buffer.data(smi + 1186);
    const auto *smi_1188 = buffer.data(smi + 1188);
    const auto *smi_1190 = buffer.data(smi + 1190);
    const auto *smi_1191 = buffer.data(smi + 1191);
    const auto *smi_1193 = buffer.data(smi + 1193);
    const auto *smi_1194 = buffer.data(smi + 1194);
    const auto *smi_1196 = buffer.data(smi + 1196);
    const auto *smi_1197 = buffer.data(smi + 1197);
    const auto *smi_1198 = buffer.data(smi + 1198);
    const auto *smi_1199 = buffer.data(smi + 1199);
    const auto *smi_1200 = buffer.data(smi + 1200);
    const auto *smi_1201 = buffer.data(smi + 1201);
    const auto *smi_1202 = buffer.data(smi + 1202);
    const auto *smi_1203 = buffer.data(smi + 1203);

    const auto *snh0_843 = buffer.data(snh0 + 843);
    const auto *snh0_845 = buffer.data(snh0 + 845);
    const auto *snh0_846 = buffer.data(snh0 + 846);
    const auto *snh0_849 = buffer.data(snh0 + 849);
    const auto *snh0_850 = buffer.data(snh0 + 850);
    const auto *snh0_852 = buffer.data(snh0 + 852);
    const auto *snh0_854 = buffer.data(snh0 + 854);
    const auto *snh0_855 = buffer.data(snh0 + 855);
    const auto *snh0_857 = buffer.data(snh0 + 857);
    const auto *snh0_858 = buffer.data(snh0 + 858);
    const auto *snh0_859 = buffer.data(snh0 + 859);
    const auto *snh0_860 = buffer.data(snh0 + 860);
    const auto *snh0_861 = buffer.data(snh0 + 861);
    const auto *snh0_864 = buffer.data(snh0 + 864);
    const auto *snh0_866 = buffer.data(snh0 + 866);
    const auto *snh0_867 = buffer.data(snh0 + 867);
    const auto *snh0_870 = buffer.data(snh0 + 870);
    const auto *snh0_871 = buffer.data(snh0 + 871);
    const auto *snh0_873 = buffer.data(snh0 + 873);
    const auto *snh0_875 = buffer.data(snh0 + 875);
    const auto *snh0_876 = buffer.data(snh0 + 876);
    const auto *snh0_878 = buffer.data(snh0 + 878);
    const auto *snh0_879 = buffer.data(snh0 + 879);
    const auto *snh0_880 = buffer.data(snh0 + 880);
    const auto *snh0_881 = buffer.data(snh0 + 881);
    const auto *snh0_882 = buffer.data(snh0 + 882);
    const auto *snh0_885 = buffer.data(snh0 + 885);
    const auto *snh0_887 = buffer.data(snh0 + 887);
    const auto *snh0_888 = buffer.data(snh0 + 888);
    const auto *snh0_891 = buffer.data(snh0 + 891);
    const auto *snh0_892 = buffer.data(snh0 + 892);
    const auto *snh0_894 = buffer.data(snh0 + 894);
    const auto *snh0_896 = buffer.data(snh0 + 896);
    const auto *snh0_897 = buffer.data(snh0 + 897);
    const auto *snh0_899 = buffer.data(snh0 + 899);
    const auto *snh0_900 = buffer.data(snh0 + 900);
    const auto *snh0_901 = buffer.data(snh0 + 901);
    const auto *snh0_902 = buffer.data(snh0 + 902);

    const auto *snh1_843 = buffer.data(snh1 + 843);
    const auto *snh1_845 = buffer.data(snh1 + 845);
    const auto *snh1_846 = buffer.data(snh1 + 846);
    const auto *snh1_849 = buffer.data(snh1 + 849);
    const auto *snh1_850 = buffer.data(snh1 + 850);
    const auto *snh1_852 = buffer.data(snh1 + 852);
    const auto *snh1_854 = buffer.data(snh1 + 854);
    const auto *snh1_855 = buffer.data(snh1 + 855);
    const auto *snh1_857 = buffer.data(snh1 + 857);
    const auto *snh1_858 = buffer.data(snh1 + 858);
    const auto *snh1_859 = buffer.data(snh1 + 859);
    const auto *snh1_860 = buffer.data(snh1 + 860);
    const auto *snh1_861 = buffer.data(snh1 + 861);
    const auto *snh1_864 = buffer.data(snh1 + 864);
    const auto *snh1_866 = buffer.data(snh1 + 866);
    const auto *snh1_867 = buffer.data(snh1 + 867);
    const auto *snh1_870 = buffer.data(snh1 + 870);
    const auto *snh1_871 = buffer.data(snh1 + 871);
    const auto *snh1_873 = buffer.data(snh1 + 873);
    const auto *snh1_875 = buffer.data(snh1 + 875);
    const auto *snh1_876 = buffer.data(snh1 + 876);
    const auto *snh1_878 = buffer.data(snh1 + 878);
    const auto *snh1_879 = buffer.data(snh1 + 879);
    const auto *snh1_880 = buffer.data(snh1 + 880);
    const auto *snh1_881 = buffer.data(snh1 + 881);
    const auto *snh1_882 = buffer.data(snh1 + 882);
    const auto *snh1_885 = buffer.data(snh1 + 885);
    const auto *snh1_887 = buffer.data(snh1 + 887);
    const auto *snh1_888 = buffer.data(snh1 + 888);
    const auto *snh1_891 = buffer.data(snh1 + 891);
    const auto *snh1_892 = buffer.data(snh1 + 892);
    const auto *snh1_894 = buffer.data(snh1 + 894);
    const auto *snh1_896 = buffer.data(snh1 + 896);
    const auto *snh1_897 = buffer.data(snh1 + 897);
    const auto *snh1_899 = buffer.data(snh1 + 899);
    const auto *snh1_900 = buffer.data(snh1 + 900);
    const auto *snh1_901 = buffer.data(snh1 + 901);
    const auto *snh1_902 = buffer.data(snh1 + 902);

    const auto *sni_1120 = buffer.data(sni + 1120);
    const auto *sni_1122 = buffer.data(sni + 1122);
    const auto *sni_1123 = buffer.data(sni + 1123);
    const auto *sni_1125 = buffer.data(sni + 1125);
    const auto *sni_1126 = buffer.data(sni + 1126);
    const auto *sni_1129 = buffer.data(sni + 1129);
    const auto *sni_1130 = buffer.data(sni + 1130);
    const auto *sni_1132 = buffer.data(sni + 1132);
    const auto *sni_1134 = buffer.data(sni + 1134);
    const auto *sni_1135 = buffer.data(sni + 1135);
    const auto *sni_1137 = buffer.data(sni + 1137);
    const auto *sni_1138 = buffer.data(sni + 1138);
    const auto *sni_1140 = buffer.data(sni + 1140);
    const auto *sni_1141 = buffer.data(sni + 1141);
    const auto *sni_1142 = buffer.data(sni + 1142);
    const auto *sni_1143 = buffer.data(sni + 1143);
    const auto *sni_1144 = buffer.data(sni + 1144);
    const auto *sni_1145 = buffer.data(sni + 1145);
    const auto *sni_1146 = buffer.data(sni + 1146);
    const auto *sni_1147 = buffer.data(sni + 1147);
    const auto *sni_1148 = buffer.data(sni + 1148);
    const auto *sni_1150 = buffer.data(sni + 1150);
    const auto *sni_1151 = buffer.data(sni + 1151);
    const auto *sni_1153 = buffer.data(sni + 1153);
    const auto *sni_1154 = buffer.data(sni + 1154);
    const auto *sni_1157 = buffer.data(sni + 1157);
    const auto *sni_1158 = buffer.data(sni + 1158);
    const auto *sni_1160 = buffer.data(sni + 1160);
    const auto *sni_1162 = buffer.data(sni + 1162);
    const auto *sni_1163 = buffer.data(sni + 1163);
    const auto *sni_1165 = buffer.data(sni + 1165);
    const auto *sni_1166 = buffer.data(sni + 1166);
    const auto *sni_1168 = buffer.data(sni + 1168);
    const auto *sni_1169 = buffer.data(sni + 1169);
    const auto *sni_1170 = buffer.data(sni + 1170);
    const auto *sni_1171 = buffer.data(sni + 1171);
    const auto *sni_1172 = buffer.data(sni + 1172);
    const auto *sni_1173 = buffer.data(sni + 1173);
    const auto *sni_1174 = buffer.data(sni + 1174);
    const auto *sni_1175 = buffer.data(sni + 1175);
    const auto *sni_1176 = buffer.data(sni + 1176);
    const auto *sni_1178 = buffer.data(sni + 1178);
    const auto *sni_1179 = buffer.data(sni + 1179);
    const auto *sni_1181 = buffer.data(sni + 1181);
    const auto *sni_1182 = buffer.data(sni + 1182);
    const auto *sni_1185 = buffer.data(sni + 1185);
    const auto *sni_1186 = buffer.data(sni + 1186);
    const auto *sni_1188 = buffer.data(sni + 1188);
    const auto *sni_1190 = buffer.data(sni + 1190);
    const auto *sni_1191 = buffer.data(sni + 1191);
    const auto *sni_1193 = buffer.data(sni + 1193);
    const auto *sni_1194 = buffer.data(sni + 1194);
    const auto *sni_1196 = buffer.data(sni + 1196);
    const auto *sni_1197 = buffer.data(sni + 1197);
    const auto *sni_1198 = buffer.data(sni + 1198);
    const auto *sni_1199 = buffer.data(sni + 1199);
    const auto *sni_1200 = buffer.data(sni + 1200);
    const auto *sni_1201 = buffer.data(sni + 1201);
    const auto *sni_1202 = buffer.data(sni + 1202);
    const auto *sni_1203 = buffer.data(sni + 1203);

#pragma omp simd aligned(t_1441, t_1442, t_1443, t_1444, pc_x, pc_y, pc_z, smi_868, smi_896, \
                         smi_898, smi_1123, snh0_843, snh1_843, sni_1120, sni_1122, \
                         sni_1123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1441[k] = f_16 * smi_896[k]
                    + f_3 * pc_y[k] * sni_1120[k];

        t_1442[k] = f_16 * smi_868[k]
                    + f_3 * pc_z[k] * sni_1120[k];

        t_1443[k] = f_14 * smi_1123[k]
                    + f_4 * snh0_843[k]
                    - f_5 * snh1_843[k]
                    + f_3 * pc_x[k] * sni_1123[k];

        t_1444[k] = f_16 * smi_898[k]
                    + f_3 * pc_y[k] * sni_1122[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pc_x, pc_z, smi_871, smi_1125, smi_1126, \
                         snh0_845, snh0_846, snh1_845, snh1_846, sni_1123, sni_1125, \
                         sni_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_14 * smi_1125[k]
                    + f_4 * snh0_845[k]
                    - f_5 * snh1_845[k]
                    + f_3 * pc_x[k] * sni_1125[k];

        t_1446[k] = f_14 * smi_1126[k]
                    + f_6 * snh0_846[k]
                    - f_7 * snh1_846[k]
                    + f_3 * pc_x[k] * sni_1126[k];

        t_1447[k] = f_16 * smi_871[k]
                    + f_3 * pc_z[k] * sni_1123[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, pc_y, smi_901, smi_1129, smi_1130, \
                         snh0_849, snh0_850, snh1_849, snh1_850, sni_1125, sni_1129, \
                         sni_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_16 * smi_901[k]
                    + f_3 * pc_y[k] * sni_1125[k];

        t_1449[k] = f_14 * smi_1129[k]
                    + f_6 * snh0_849[k]
                    - f_7 * snh1_849[k]
                    + f_3 * pc_x[k] * sni_1129[k];

        t_1450[k] = f_14 * smi_1130[k]
                    + f_8 * snh0_850[k]
                    - f_9 * snh1_850[k]
                    + f_3 * pc_x[k] * sni_1130[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_x, pc_y, pc_z, smi_874, smi_905, smi_1132, \
                         snh0_852, snh1_852, sni_1126, sni_1129, \
                         sni_1132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_16 * smi_874[k]
                    + f_3 * pc_z[k] * sni_1126[k];

        t_1452[k] = f_14 * smi_1132[k]
                    + f_8 * snh0_852[k]
                    - f_9 * snh1_852[k]
                    + f_3 * pc_x[k] * sni_1132[k];

        t_1453[k] = f_16 * smi_905[k]
                    + f_3 * pc_y[k] * sni_1129[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, pc_z, smi_878, smi_1134, smi_1135, \
                         snh0_854, snh0_855, snh1_854, snh1_855, sni_1130, sni_1134, \
                         sni_1135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_14 * smi_1134[k]
                    + f_8 * snh0_854[k]
                    - f_9 * snh1_854[k]
                    + f_3 * pc_x[k] * sni_1134[k];

        t_1455[k] = f_14 * smi_1135[k]
                    + f_10 * snh0_855[k]
                    - f_11 * snh1_855[k]
                    + f_3 * pc_x[k] * sni_1135[k];

        t_1456[k] = f_16 * smi_878[k]
                    + f_3 * pc_z[k] * sni_1130[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pc_x, pc_y, smi_910, smi_1137, smi_1138, \
                         snh0_857, snh0_858, snh1_857, snh1_858, sni_1134, sni_1137, \
                         sni_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_14 * smi_1137[k]
                    + f_10 * snh0_857[k]
                    - f_11 * snh1_857[k]
                    + f_3 * pc_x[k] * sni_1137[k];

        t_1458[k] = f_14 * smi_1138[k]
                    + f_10 * snh0_858[k]
                    - f_11 * snh1_858[k]
                    + f_3 * pc_x[k] * sni_1138[k];

        t_1459[k] = f_16 * smi_910[k]
                    + f_3 * pc_y[k] * sni_1134[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, t_1463, pc_x, smi_1140, smi_1141, smi_1142, \
                         smi_1143, snh0_860, snh1_860, sni_1140, sni_1141, sni_1142, \
                         sni_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_14 * smi_1140[k]
                    + f_10 * snh0_860[k]
                    - f_11 * snh1_860[k]
                    + f_3 * pc_x[k] * sni_1140[k];

        t_1461[k] = f_14 * smi_1141[k]
                    + f_3 * pc_x[k] * sni_1141[k];

        t_1462[k] = f_14 * smi_1142[k]
                    + f_3 * pc_x[k] * sni_1142[k];

        t_1463[k] = f_14 * smi_1143[k]
                    + f_3 * pc_x[k] * sni_1143[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, t_1467, pc_x, smi_1144, smi_1145, smi_1146, \
                         smi_1147, sni_1144, sni_1145, sni_1146, \
                         sni_1147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_14 * smi_1144[k]
                    + f_3 * pc_x[k] * sni_1144[k];

        t_1465[k] = f_14 * smi_1145[k]
                    + f_3 * pc_x[k] * sni_1145[k];

        t_1466[k] = f_14 * smi_1146[k]
                    + f_3 * pc_x[k] * sni_1146[k];

        t_1467[k] = f_14 * smi_1147[k]
                    + f_3 * pc_x[k] * sni_1147[k];
    }

#pragma omp simd aligned(t_1468, t_1469, t_1470, pc_y, pc_z, smi_889, smi_917, smi_919, \
                         snh0_855, snh0_857, snh1_855, snh1_857, sni_1141, \
                         sni_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1468[k] = f_16 * smi_917[k]
                    + f_1 * snh0_855[k]
                    - f_2 * snh1_855[k]
                    + f_3 * pc_y[k] * sni_1141[k];

        t_1469[k] = f_16 * smi_889[k]
                    + f_3 * pc_z[k] * sni_1141[k];

        t_1470[k] = f_16 * smi_919[k]
                    + f_4 * snh0_857[k]
                    - f_5 * snh1_857[k]
                    + f_3 * pc_y[k] * sni_1143[k];
    }

#pragma omp simd aligned(t_1471, t_1472, t_1473, pc_y, smi_920, smi_921, smi_922, snh0_858, \
                         snh0_859, snh0_860, snh1_858, snh1_859, snh1_860, sni_1144, sni_1145, \
                         sni_1146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1471[k] = f_16 * smi_920[k]
                    + f_6 * snh0_858[k]
                    - f_7 * snh1_858[k]
                    + f_3 * pc_y[k] * sni_1144[k];

        t_1472[k] = f_16 * smi_921[k]
                    + f_8 * snh0_859[k]
                    - f_9 * snh1_859[k]
                    + f_3 * pc_y[k] * sni_1145[k];

        t_1473[k] = f_16 * smi_922[k]
                    + f_10 * snh0_860[k]
                    - f_11 * snh1_860[k]
                    + f_3 * pc_y[k] * sni_1146[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, pc_x, pc_y, pc_z, smi_895, smi_923, smi_1148, \
                         snh0_860, snh0_861, snh1_860, snh1_861, sni_1147, \
                         sni_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_16 * smi_923[k]
                    + f_3 * pc_y[k] * sni_1147[k];

        t_1475[k] = f_16 * smi_895[k]
                    + f_1 * snh0_860[k]
                    - f_2 * snh1_860[k]
                    + f_3 * pc_z[k] * sni_1147[k];

        t_1476[k] = f_14 * smi_1148[k]
                    + f_1 * snh0_861[k]
                    - f_2 * snh1_861[k]
                    + f_3 * pc_x[k] * sni_1148[k];
    }

#pragma omp simd aligned(t_1477, t_1478, t_1479, t_1480, pc_x, pc_y, pc_z, smi_896, smi_924, \
                         smi_926, smi_1151, snh0_864, snh1_864, sni_1148, sni_1150, \
                         sni_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1477[k] = f_15 * smi_924[k]
                    + f_3 * pc_y[k] * sni_1148[k];

        t_1478[k] = f_17 * smi_896[k]
                    + f_3 * pc_z[k] * sni_1148[k];

        t_1479[k] = f_14 * smi_1151[k]
                    + f_4 * snh0_864[k]
                    - f_5 * snh1_864[k]
                    + f_3 * pc_x[k] * sni_1151[k];

        t_1480[k] = f_15 * smi_926[k]
                    + f_3 * pc_y[k] * sni_1150[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, pc_x, pc_z, smi_899, smi_1153, smi_1154, \
                         snh0_866, snh0_867, snh1_866, snh1_867, sni_1151, sni_1153, \
                         sni_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_14 * smi_1153[k]
                    + f_4 * snh0_866[k]
                    - f_5 * snh1_866[k]
                    + f_3 * pc_x[k] * sni_1153[k];

        t_1482[k] = f_14 * smi_1154[k]
                    + f_6 * snh0_867[k]
                    - f_7 * snh1_867[k]
                    + f_3 * pc_x[k] * sni_1154[k];

        t_1483[k] = f_17 * smi_899[k]
                    + f_3 * pc_z[k] * sni_1151[k];
    }

#pragma omp simd aligned(t_1484, t_1485, t_1486, pc_x, pc_y, smi_929, smi_1157, smi_1158, \
                         snh0_870, snh0_871, snh1_870, snh1_871, sni_1153, sni_1157, \
                         sni_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1484[k] = f_15 * smi_929[k]
                    + f_3 * pc_y[k] * sni_1153[k];

        t_1485[k] = f_14 * smi_1157[k]
                    + f_6 * snh0_870[k]
                    - f_7 * snh1_870[k]
                    + f_3 * pc_x[k] * sni_1157[k];

        t_1486[k] = f_14 * smi_1158[k]
                    + f_8 * snh0_871[k]
                    - f_9 * snh1_871[k]
                    + f_3 * pc_x[k] * sni_1158[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pc_x, pc_y, pc_z, smi_902, smi_933, smi_1160, \
                         snh0_873, snh1_873, sni_1154, sni_1157, \
                         sni_1160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_17 * smi_902[k]
                    + f_3 * pc_z[k] * sni_1154[k];

        t_1488[k] = f_14 * smi_1160[k]
                    + f_8 * snh0_873[k]
                    - f_9 * snh1_873[k]
                    + f_3 * pc_x[k] * sni_1160[k];

        t_1489[k] = f_15 * smi_933[k]
                    + f_3 * pc_y[k] * sni_1157[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pc_x, pc_z, smi_906, smi_1162, smi_1163, \
                         snh0_875, snh0_876, snh1_875, snh1_876, sni_1158, sni_1162, \
                         sni_1163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = f_14 * smi_1162[k]
                    + f_8 * snh0_875[k]
                    - f_9 * snh1_875[k]
                    + f_3 * pc_x[k] * sni_1162[k];

        t_1491[k] = f_14 * smi_1163[k]
                    + f_10 * snh0_876[k]
                    - f_11 * snh1_876[k]
                    + f_3 * pc_x[k] * sni_1163[k];

        t_1492[k] = f_17 * smi_906[k]
                    + f_3 * pc_z[k] * sni_1158[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pc_x, pc_y, smi_938, smi_1165, smi_1166, \
                         snh0_878, snh0_879, snh1_878, snh1_879, sni_1162, sni_1165, \
                         sni_1166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = f_14 * smi_1165[k]
                    + f_10 * snh0_878[k]
                    - f_11 * snh1_878[k]
                    + f_3 * pc_x[k] * sni_1165[k];

        t_1494[k] = f_14 * smi_1166[k]
                    + f_10 * snh0_879[k]
                    - f_11 * snh1_879[k]
                    + f_3 * pc_x[k] * sni_1166[k];

        t_1495[k] = f_15 * smi_938[k]
                    + f_3 * pc_y[k] * sni_1162[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, t_1499, pc_x, smi_1168, smi_1169, smi_1170, \
                         smi_1171, snh0_881, snh1_881, sni_1168, sni_1169, sni_1170, \
                         sni_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = f_14 * smi_1168[k]
                    + f_10 * snh0_881[k]
                    - f_11 * snh1_881[k]
                    + f_3 * pc_x[k] * sni_1168[k];

        t_1497[k] = f_14 * smi_1169[k]
                    + f_3 * pc_x[k] * sni_1169[k];

        t_1498[k] = f_14 * smi_1170[k]
                    + f_3 * pc_x[k] * sni_1170[k];

        t_1499[k] = f_14 * smi_1171[k]
                    + f_3 * pc_x[k] * sni_1171[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, t_1503, pc_x, smi_1172, smi_1173, smi_1174, \
                         smi_1175, sni_1172, sni_1173, sni_1174, \
                         sni_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_14 * smi_1172[k]
                    + f_3 * pc_x[k] * sni_1172[k];

        t_1501[k] = f_14 * smi_1173[k]
                    + f_3 * pc_x[k] * sni_1173[k];

        t_1502[k] = f_14 * smi_1174[k]
                    + f_3 * pc_x[k] * sni_1174[k];

        t_1503[k] = f_14 * smi_1175[k]
                    + f_3 * pc_x[k] * sni_1175[k];
    }

#pragma omp simd aligned(t_1504, t_1505, t_1506, pc_y, pc_z, smi_917, smi_945, smi_947, \
                         snh0_876, snh0_878, snh1_876, snh1_878, sni_1169, \
                         sni_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1504[k] = f_15 * smi_945[k]
                    + f_1 * snh0_876[k]
                    - f_2 * snh1_876[k]
                    + f_3 * pc_y[k] * sni_1169[k];

        t_1505[k] = f_17 * smi_917[k]
                    + f_3 * pc_z[k] * sni_1169[k];

        t_1506[k] = f_15 * smi_947[k]
                    + f_4 * snh0_878[k]
                    - f_5 * snh1_878[k]
                    + f_3 * pc_y[k] * sni_1171[k];
    }

#pragma omp simd aligned(t_1507, t_1508, t_1509, pc_y, smi_948, smi_949, smi_950, snh0_879, \
                         snh0_880, snh0_881, snh1_879, snh1_880, snh1_881, sni_1172, sni_1173, \
                         sni_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1507[k] = f_15 * smi_948[k]
                    + f_6 * snh0_879[k]
                    - f_7 * snh1_879[k]
                    + f_3 * pc_y[k] * sni_1172[k];

        t_1508[k] = f_15 * smi_949[k]
                    + f_8 * snh0_880[k]
                    - f_9 * snh1_880[k]
                    + f_3 * pc_y[k] * sni_1173[k];

        t_1509[k] = f_15 * smi_950[k]
                    + f_10 * snh0_881[k]
                    - f_11 * snh1_881[k]
                    + f_3 * pc_y[k] * sni_1174[k];
    }

#pragma omp simd aligned(t_1510, t_1511, t_1512, pc_x, pc_y, pc_z, smi_923, smi_951, smi_1176, \
                         snh0_881, snh0_882, snh1_881, snh1_882, sni_1175, \
                         sni_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1510[k] = f_15 * smi_951[k]
                    + f_3 * pc_y[k] * sni_1175[k];

        t_1511[k] = f_17 * smi_923[k]
                    + f_1 * snh0_881[k]
                    - f_2 * snh1_881[k]
                    + f_3 * pc_z[k] * sni_1175[k];

        t_1512[k] = f_14 * smi_1176[k]
                    + f_1 * snh0_882[k]
                    - f_2 * snh1_882[k]
                    + f_3 * pc_x[k] * sni_1176[k];
    }

#pragma omp simd aligned(t_1513, t_1514, t_1515, t_1516, pc_x, pc_y, pc_z, smi_924, smi_952, \
                         smi_954, smi_1179, snh0_885, snh1_885, sni_1176, sni_1178, \
                         sni_1179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1513[k] = f_14 * smi_952[k]
                    + f_3 * pc_y[k] * sni_1176[k];

        t_1514[k] = f_21 * smi_924[k]
                    + f_3 * pc_z[k] * sni_1176[k];

        t_1515[k] = f_14 * smi_1179[k]
                    + f_4 * snh0_885[k]
                    - f_5 * snh1_885[k]
                    + f_3 * pc_x[k] * sni_1179[k];

        t_1516[k] = f_14 * smi_954[k]
                    + f_3 * pc_y[k] * sni_1178[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, pc_x, pc_z, smi_927, smi_1181, smi_1182, \
                         snh0_887, snh0_888, snh1_887, snh1_888, sni_1179, sni_1181, \
                         sni_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = f_14 * smi_1181[k]
                    + f_4 * snh0_887[k]
                    - f_5 * snh1_887[k]
                    + f_3 * pc_x[k] * sni_1181[k];

        t_1518[k] = f_14 * smi_1182[k]
                    + f_6 * snh0_888[k]
                    - f_7 * snh1_888[k]
                    + f_3 * pc_x[k] * sni_1182[k];

        t_1519[k] = f_21 * smi_927[k]
                    + f_3 * pc_z[k] * sni_1179[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, pc_x, pc_y, smi_957, smi_1185, smi_1186, \
                         snh0_891, snh0_892, snh1_891, snh1_892, sni_1181, sni_1185, \
                         sni_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_14 * smi_957[k]
                    + f_3 * pc_y[k] * sni_1181[k];

        t_1521[k] = f_14 * smi_1185[k]
                    + f_6 * snh0_891[k]
                    - f_7 * snh1_891[k]
                    + f_3 * pc_x[k] * sni_1185[k];

        t_1522[k] = f_14 * smi_1186[k]
                    + f_8 * snh0_892[k]
                    - f_9 * snh1_892[k]
                    + f_3 * pc_x[k] * sni_1186[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pc_x, pc_y, pc_z, smi_930, smi_961, smi_1188, \
                         snh0_894, snh1_894, sni_1182, sni_1185, \
                         sni_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_21 * smi_930[k]
                    + f_3 * pc_z[k] * sni_1182[k];

        t_1524[k] = f_14 * smi_1188[k]
                    + f_8 * snh0_894[k]
                    - f_9 * snh1_894[k]
                    + f_3 * pc_x[k] * sni_1188[k];

        t_1525[k] = f_14 * smi_961[k]
                    + f_3 * pc_y[k] * sni_1185[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, pc_x, pc_z, smi_934, smi_1190, smi_1191, \
                         snh0_896, snh0_897, snh1_896, snh1_897, sni_1186, sni_1190, \
                         sni_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_14 * smi_1190[k]
                    + f_8 * snh0_896[k]
                    - f_9 * snh1_896[k]
                    + f_3 * pc_x[k] * sni_1190[k];

        t_1527[k] = f_14 * smi_1191[k]
                    + f_10 * snh0_897[k]
                    - f_11 * snh1_897[k]
                    + f_3 * pc_x[k] * sni_1191[k];

        t_1528[k] = f_21 * smi_934[k]
                    + f_3 * pc_z[k] * sni_1186[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, pc_x, pc_y, smi_966, smi_1193, smi_1194, \
                         snh0_899, snh0_900, snh1_899, snh1_900, sni_1190, sni_1193, \
                         sni_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_14 * smi_1193[k]
                    + f_10 * snh0_899[k]
                    - f_11 * snh1_899[k]
                    + f_3 * pc_x[k] * sni_1193[k];

        t_1530[k] = f_14 * smi_1194[k]
                    + f_10 * snh0_900[k]
                    - f_11 * snh1_900[k]
                    + f_3 * pc_x[k] * sni_1194[k];

        t_1531[k] = f_14 * smi_966[k]
                    + f_3 * pc_y[k] * sni_1190[k];
    }

#pragma omp simd aligned(t_1532, t_1533, t_1534, t_1535, pc_x, smi_1196, smi_1197, smi_1198, \
                         smi_1199, snh0_902, snh1_902, sni_1196, sni_1197, sni_1198, \
                         sni_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1532[k] = f_14 * smi_1196[k]
                    + f_10 * snh0_902[k]
                    - f_11 * snh1_902[k]
                    + f_3 * pc_x[k] * sni_1196[k];

        t_1533[k] = f_14 * smi_1197[k]
                    + f_3 * pc_x[k] * sni_1197[k];

        t_1534[k] = f_14 * smi_1198[k]
                    + f_3 * pc_x[k] * sni_1198[k];

        t_1535[k] = f_14 * smi_1199[k]
                    + f_3 * pc_x[k] * sni_1199[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, t_1539, pc_x, smi_1200, smi_1201, smi_1202, \
                         smi_1203, sni_1200, sni_1201, sni_1202, \
                         sni_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_14 * smi_1200[k]
                    + f_3 * pc_x[k] * sni_1200[k];

        t_1537[k] = f_14 * smi_1201[k]
                    + f_3 * pc_x[k] * sni_1201[k];

        t_1538[k] = f_14 * smi_1202[k]
                    + f_3 * pc_x[k] * sni_1202[k];

        t_1539[k] = f_14 * smi_1203[k]
                    + f_3 * pc_x[k] * sni_1203[k];
    }

#pragma omp simd aligned(t_1540, t_1541, t_1542, pc_y, pc_z, smi_945, smi_973, smi_975, \
                         snh0_897, snh0_899, snh1_897, snh1_899, sni_1197, \
                         sni_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1540[k] = f_14 * smi_973[k]
                    + f_1 * snh0_897[k]
                    - f_2 * snh1_897[k]
                    + f_3 * pc_y[k] * sni_1197[k];

        t_1541[k] = f_21 * smi_945[k]
                    + f_3 * pc_z[k] * sni_1197[k];

        t_1542[k] = f_14 * smi_975[k]
                    + f_4 * snh0_899[k]
                    - f_5 * snh1_899[k]
                    + f_3 * pc_y[k] * sni_1199[k];
    }

#pragma omp simd aligned(t_1543, t_1544, t_1545, pc_y, smi_976, smi_977, smi_978, snh0_900, \
                         snh0_901, snh0_902, snh1_900, snh1_901, snh1_902, sni_1200, sni_1201, \
                         sni_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1543[k] = f_14 * smi_976[k]
                    + f_6 * snh0_900[k]
                    - f_7 * snh1_900[k]
                    + f_3 * pc_y[k] * sni_1200[k];

        t_1544[k] = f_14 * smi_977[k]
                    + f_8 * snh0_901[k]
                    - f_9 * snh1_901[k]
                    + f_3 * pc_y[k] * sni_1201[k];

        t_1545[k] = f_14 * smi_978[k]
                    + f_10 * snh0_902[k]
                    - f_11 * snh1_902[k]
                    + f_3 * pc_y[k] * sni_1202[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smk0,
                                                           const size_t smi, const size_t smk1,
                                                           const size_t snh0, const size_t snh1,
                                                           const size_t sni, const size_t ncols,
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smk0_1260 = buffer.data(smk0 + 1260);
    const auto *smk0_1263 = buffer.data(smk0 + 1263);
    const auto *smk0_1265 = buffer.data(smk0 + 1265);
    const auto *smk0_1266 = buffer.data(smk0 + 1266);
    const auto *smk0_1269 = buffer.data(smk0 + 1269);
    const auto *smk0_1270 = buffer.data(smk0 + 1270);
    const auto *smk0_1272 = buffer.data(smk0 + 1272);
    const auto *smk0_1274 = buffer.data(smk0 + 1274);
    const auto *smk0_1275 = buffer.data(smk0 + 1275);
    const auto *smk0_1277 = buffer.data(smk0 + 1277);
    const auto *smk0_1278 = buffer.data(smk0 + 1278);
    const auto *smk0_1280 = buffer.data(smk0 + 1280);
    const auto *smk0_1295 = buffer.data(smk0 + 1295);
    const auto *smk0_1296 = buffer.data(smk0 + 1296);
    const auto *smk0_1299 = buffer.data(smk0 + 1299);
    const auto *smk0_1302 = buffer.data(smk0 + 1302);
    const auto *smk0_1620 = buffer.data(smk0 + 1620);
    const auto *smk0_1623 = buffer.data(smk0 + 1623);
    const auto *smk0_1625 = buffer.data(smk0 + 1625);
    const auto *smk0_1626 = buffer.data(smk0 + 1626);
    const auto *smk0_1629 = buffer.data(smk0 + 1629);
    const auto *smk0_1630 = buffer.data(smk0 + 1630);
    const auto *smk0_1632 = buffer.data(smk0 + 1632);
    const auto *smk0_1634 = buffer.data(smk0 + 1634);
    const auto *smk0_1635 = buffer.data(smk0 + 1635);
    const auto *smk0_1637 = buffer.data(smk0 + 1637);
    const auto *smk0_1638 = buffer.data(smk0 + 1638);
    const auto *smk0_1640 = buffer.data(smk0 + 1640);
    const auto *smk0_1648 = buffer.data(smk0 + 1648);
    const auto *smk0_1650 = buffer.data(smk0 + 1650);
    const auto *smk0_1651 = buffer.data(smk0 + 1651);
    const auto *smk0_1652 = buffer.data(smk0 + 1652);
    const auto *smk0_1653 = buffer.data(smk0 + 1653);
    const auto *smk0_1655 = buffer.data(smk0 + 1655);
    const auto *smk0_1661 = buffer.data(smk0 + 1661);

    const auto *smi_951 = buffer.data(smi + 951);
    const auto *smi_952 = buffer.data(smi + 952);
    const auto *smi_955 = buffer.data(smi + 955);
    const auto *smi_958 = buffer.data(smi + 958);
    const auto *smi_962 = buffer.data(smi + 962);
    const auto *smi_973 = buffer.data(smi + 973);
    const auto *smi_979 = buffer.data(smi + 979);
    const auto *smi_980 = buffer.data(smi + 980);
    const auto *smi_981 = buffer.data(smi + 981);
    const auto *smi_982 = buffer.data(smi + 982);
    const auto *smi_983 = buffer.data(smi + 983);
    const auto *smi_985 = buffer.data(smi + 985);
    const auto *smi_986 = buffer.data(smi + 986);
    const auto *smi_988 = buffer.data(smi + 988);
    const auto *smi_989 = buffer.data(smi + 989);
    const auto *smi_990 = buffer.data(smi + 990);
    const auto *smi_992 = buffer.data(smi + 992);
    const auto *smi_993 = buffer.data(smi + 993);
    const auto *smi_994 = buffer.data(smi + 994);
    const auto *smi_1001 = buffer.data(smi + 1001);
    const auto *smi_1003 = buffer.data(smi + 1003);
    const auto *smi_1004 = buffer.data(smi + 1004);
    const auto *smi_1005 = buffer.data(smi + 1005);
    const auto *smi_1006 = buffer.data(smi + 1006);
    const auto *smi_1007 = buffer.data(smi + 1007);
    const auto *smi_1008 = buffer.data(smi + 1008);
    const auto *smi_1010 = buffer.data(smi + 1010);
    const auto *smi_1013 = buffer.data(smi + 1013);
    const auto *smi_1017 = buffer.data(smi + 1017);
    const auto *smi_1022 = buffer.data(smi + 1022);
    const auto *smi_1035 = buffer.data(smi + 1035);
    const auto *smi_1036 = buffer.data(smi + 1036);
    const auto *smi_1038 = buffer.data(smi + 1038);
    const auto *smi_1225 = buffer.data(smi + 1225);
    const auto *smi_1226 = buffer.data(smi + 1226);
    const auto *smi_1227 = buffer.data(smi + 1227);
    const auto *smi_1228 = buffer.data(smi + 1228);
    const auto *smi_1229 = buffer.data(smi + 1229);
    const auto *smi_1230 = buffer.data(smi + 1230);
    const auto *smi_1231 = buffer.data(smi + 1231);
    const auto *smi_1232 = buffer.data(smi + 1232);
    const auto *smi_1235 = buffer.data(smi + 1235);
    const auto *smi_1237 = buffer.data(smi + 1237);
    const auto *smi_1238 = buffer.data(smi + 1238);
    const auto *smi_1241 = buffer.data(smi + 1241);
    const auto *smi_1242 = buffer.data(smi + 1242);
    const auto *smi_1244 = buffer.data(smi + 1244);
    const auto *smi_1246 = buffer.data(smi + 1246);
    const auto *smi_1247 = buffer.data(smi + 1247);
    const auto *smi_1249 = buffer.data(smi + 1249);
    const auto *smi_1250 = buffer.data(smi + 1250);
    const auto *smi_1252 = buffer.data(smi + 1252);
    const auto *smi_1253 = buffer.data(smi + 1253);
    const auto *smi_1254 = buffer.data(smi + 1254);
    const auto *smi_1255 = buffer.data(smi + 1255);
    const auto *smi_1256 = buffer.data(smi + 1256);
    const auto *smi_1257 = buffer.data(smi + 1257);
    const auto *smi_1258 = buffer.data(smi + 1258);
    const auto *smi_1259 = buffer.data(smi + 1259);
    const auto *smi_1260 = buffer.data(smi + 1260);
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
    const auto *smi_1293 = buffer.data(smi + 1293);

    const auto *smk1_1260 = buffer.data(smk1 + 1260);
    const auto *smk1_1263 = buffer.data(smk1 + 1263);
    const auto *smk1_1265 = buffer.data(smk1 + 1265);
    const auto *smk1_1266 = buffer.data(smk1 + 1266);
    const auto *smk1_1269 = buffer.data(smk1 + 1269);
    const auto *smk1_1270 = buffer.data(smk1 + 1270);
    const auto *smk1_1272 = buffer.data(smk1 + 1272);
    const auto *smk1_1274 = buffer.data(smk1 + 1274);
    const auto *smk1_1275 = buffer.data(smk1 + 1275);
    const auto *smk1_1277 = buffer.data(smk1 + 1277);
    const auto *smk1_1278 = buffer.data(smk1 + 1278);
    const auto *smk1_1280 = buffer.data(smk1 + 1280);
    const auto *smk1_1295 = buffer.data(smk1 + 1295);
    const auto *smk1_1296 = buffer.data(smk1 + 1296);
    const auto *smk1_1299 = buffer.data(smk1 + 1299);
    const auto *smk1_1302 = buffer.data(smk1 + 1302);
    const auto *smk1_1620 = buffer.data(smk1 + 1620);
    const auto *smk1_1623 = buffer.data(smk1 + 1623);
    const auto *smk1_1625 = buffer.data(smk1 + 1625);
    const auto *smk1_1626 = buffer.data(smk1 + 1626);
    const auto *smk1_1629 = buffer.data(smk1 + 1629);
    const auto *smk1_1630 = buffer.data(smk1 + 1630);
    const auto *smk1_1632 = buffer.data(smk1 + 1632);
    const auto *smk1_1634 = buffer.data(smk1 + 1634);
    const auto *smk1_1635 = buffer.data(smk1 + 1635);
    const auto *smk1_1637 = buffer.data(smk1 + 1637);
    const auto *smk1_1638 = buffer.data(smk1 + 1638);
    const auto *smk1_1640 = buffer.data(smk1 + 1640);
    const auto *smk1_1648 = buffer.data(smk1 + 1648);
    const auto *smk1_1650 = buffer.data(smk1 + 1650);
    const auto *smk1_1651 = buffer.data(smk1 + 1651);
    const auto *smk1_1652 = buffer.data(smk1 + 1652);
    const auto *smk1_1653 = buffer.data(smk1 + 1653);
    const auto *smk1_1655 = buffer.data(smk1 + 1655);
    const auto *smk1_1661 = buffer.data(smk1 + 1661);

    const auto *snh0_902 = buffer.data(snh0 + 902);
    const auto *snh0_918 = buffer.data(snh0 + 918);
    const auto *snh0_920 = buffer.data(snh0 + 920);
    const auto *snh0_921 = buffer.data(snh0 + 921);
    const auto *snh0_922 = buffer.data(snh0 + 922);
    const auto *snh0_923 = buffer.data(snh0 + 923);
    const auto *snh0_924 = buffer.data(snh0 + 924);
    const auto *snh0_927 = buffer.data(snh0 + 927);
    const auto *snh0_929 = buffer.data(snh0 + 929);
    const auto *snh0_930 = buffer.data(snh0 + 930);
    const auto *snh0_933 = buffer.data(snh0 + 933);
    const auto *snh0_934 = buffer.data(snh0 + 934);
    const auto *snh0_936 = buffer.data(snh0 + 936);
    const auto *snh0_938 = buffer.data(snh0 + 938);
    const auto *snh0_939 = buffer.data(snh0 + 939);
    const auto *snh0_941 = buffer.data(snh0 + 941);
    const auto *snh0_942 = buffer.data(snh0 + 942);
    const auto *snh0_943 = buffer.data(snh0 + 943);
    const auto *snh0_944 = buffer.data(snh0 + 944);

    const auto *snh1_902 = buffer.data(snh1 + 902);
    const auto *snh1_918 = buffer.data(snh1 + 918);
    const auto *snh1_920 = buffer.data(snh1 + 920);
    const auto *snh1_921 = buffer.data(snh1 + 921);
    const auto *snh1_922 = buffer.data(snh1 + 922);
    const auto *snh1_923 = buffer.data(snh1 + 923);
    const auto *snh1_924 = buffer.data(snh1 + 924);
    const auto *snh1_927 = buffer.data(snh1 + 927);
    const auto *snh1_929 = buffer.data(snh1 + 929);
    const auto *snh1_930 = buffer.data(snh1 + 930);
    const auto *snh1_933 = buffer.data(snh1 + 933);
    const auto *snh1_934 = buffer.data(snh1 + 934);
    const auto *snh1_936 = buffer.data(snh1 + 936);
    const auto *snh1_938 = buffer.data(snh1 + 938);
    const auto *snh1_939 = buffer.data(snh1 + 939);
    const auto *snh1_941 = buffer.data(snh1 + 941);
    const auto *snh1_942 = buffer.data(snh1 + 942);
    const auto *snh1_943 = buffer.data(snh1 + 943);
    const auto *snh1_944 = buffer.data(snh1 + 944);

    const auto *sni_1203 = buffer.data(sni + 1203);
    const auto *sni_1204 = buffer.data(sni + 1204);
    const auto *sni_1206 = buffer.data(sni + 1206);
    const auto *sni_1207 = buffer.data(sni + 1207);
    const auto *sni_1209 = buffer.data(sni + 1209);
    const auto *sni_1210 = buffer.data(sni + 1210);
    const auto *sni_1213 = buffer.data(sni + 1213);
    const auto *sni_1214 = buffer.data(sni + 1214);
    const auto *sni_1218 = buffer.data(sni + 1218);
    const auto *sni_1225 = buffer.data(sni + 1225);
    const auto *sni_1226 = buffer.data(sni + 1226);
    const auto *sni_1227 = buffer.data(sni + 1227);
    const auto *sni_1228 = buffer.data(sni + 1228);
    const auto *sni_1229 = buffer.data(sni + 1229);
    const auto *sni_1230 = buffer.data(sni + 1230);
    const auto *sni_1231 = buffer.data(sni + 1231);
    const auto *sni_1232 = buffer.data(sni + 1232);
    const auto *sni_1234 = buffer.data(sni + 1234);
    const auto *sni_1235 = buffer.data(sni + 1235);
    const auto *sni_1237 = buffer.data(sni + 1237);
    const auto *sni_1238 = buffer.data(sni + 1238);
    const auto *sni_1241 = buffer.data(sni + 1241);
    const auto *sni_1242 = buffer.data(sni + 1242);
    const auto *sni_1244 = buffer.data(sni + 1244);
    const auto *sni_1246 = buffer.data(sni + 1246);
    const auto *sni_1247 = buffer.data(sni + 1247);
    const auto *sni_1249 = buffer.data(sni + 1249);
    const auto *sni_1250 = buffer.data(sni + 1250);
    const auto *sni_1252 = buffer.data(sni + 1252);
    const auto *sni_1253 = buffer.data(sni + 1253);
    const auto *sni_1254 = buffer.data(sni + 1254);
    const auto *sni_1255 = buffer.data(sni + 1255);
    const auto *sni_1256 = buffer.data(sni + 1256);
    const auto *sni_1257 = buffer.data(sni + 1257);
    const auto *sni_1258 = buffer.data(sni + 1258);
    const auto *sni_1259 = buffer.data(sni + 1259);
    const auto *sni_1260 = buffer.data(sni + 1260);
    const auto *sni_1262 = buffer.data(sni + 1262);
    const auto *sni_1263 = buffer.data(sni + 1263);
    const auto *sni_1265 = buffer.data(sni + 1265);
    const auto *sni_1266 = buffer.data(sni + 1266);
    const auto *sni_1269 = buffer.data(sni + 1269);
    const auto *sni_1270 = buffer.data(sni + 1270);
    const auto *sni_1274 = buffer.data(sni + 1274);
    const auto *sni_1281 = buffer.data(sni + 1281);
    const auto *sni_1282 = buffer.data(sni + 1282);
    const auto *sni_1283 = buffer.data(sni + 1283);
    const auto *sni_1284 = buffer.data(sni + 1284);
    const auto *sni_1285 = buffer.data(sni + 1285);
    const auto *sni_1286 = buffer.data(sni + 1286);
    const auto *sni_1287 = buffer.data(sni + 1287);
    const auto *sni_1288 = buffer.data(sni + 1288);
    const auto *sni_1290 = buffer.data(sni + 1290);

#pragma omp simd aligned(t_1546, t_1547, t_1548, t_1549, pb_y, pc_y, pc_z, smk0_1260, smi_951, \
                         smi_979, smi_980, smk1_1260, snh0_902, snh1_902, sni_1203, \
                         sni_1204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1546[k] = f_14 * smi_979[k]
                    + f_3 * pc_y[k] * sni_1203[k];

        t_1547[k] = f_21 * smi_951[k]
                    + f_1 * snh0_902[k]
                    - f_2 * snh1_902[k]
                    + f_3 * pc_z[k] * sni_1203[k];

        t_1548[k] = pb_y[k] * smk0_1260[k]
                    - f_12 * pc_y[k] * smk1_1260[k];

        t_1549[k] = f_13 * smi_980[k]
                    + f_3 * pc_y[k] * sni_1204[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, t_1553, pb_y, pc_y, pc_z, smk0_1263, \
                         smk0_1265, smi_952, smi_981, smi_982, smk1_1263, smk1_1265, sni_1204, \
                         sni_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = f_20 * smi_952[k]
                    + f_3 * pc_z[k] * sni_1204[k];

        t_1551[k] = pb_y[k] * smk0_1263[k]
                    + f_14 * smi_981[k]
                    - f_12 * pc_y[k] * smk1_1263[k];

        t_1552[k] = f_13 * smi_982[k]
                    + f_3 * pc_y[k] * sni_1206[k];

        t_1553[k] = pb_y[k] * smk0_1265[k]
                    - f_12 * pc_y[k] * smk1_1265[k];
    }

#pragma omp simd aligned(t_1554, t_1555, t_1556, t_1557, pb_y, pc_y, pc_z, smk0_1266, \
                         smk0_1269, smi_955, smi_983, smi_985, smk1_1266, smk1_1269, sni_1207, \
                         sni_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1554[k] = pb_y[k] * smk0_1266[k]
                    + f_15 * smi_983[k]
                    - f_12 * pc_y[k] * smk1_1266[k];

        t_1555[k] = f_20 * smi_955[k]
                    + f_3 * pc_z[k] * sni_1207[k];

        t_1556[k] = f_13 * smi_985[k]
                    + f_3 * pc_y[k] * sni_1209[k];

        t_1557[k] = pb_y[k] * smk0_1269[k]
                    - f_12 * pc_y[k] * smk1_1269[k];
    }

#pragma omp simd aligned(t_1558, t_1559, t_1560, pb_y, pc_y, pc_z, smk0_1270, smk0_1272, \
                         smi_958, smi_986, smi_988, smk1_1270, smk1_1272, \
                         sni_1210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1558[k] = pb_y[k] * smk0_1270[k]
                    + f_16 * smi_986[k]
                    - f_12 * pc_y[k] * smk1_1270[k];

        t_1559[k] = f_20 * smi_958[k]
                    + f_3 * pc_z[k] * sni_1210[k];

        t_1560[k] = pb_y[k] * smk0_1272[k]
                    + f_14 * smi_988[k]
                    - f_12 * pc_y[k] * smk1_1272[k];
    }

#pragma omp simd aligned(t_1561, t_1562, t_1563, t_1564, pb_y, pc_y, pc_z, smk0_1274, \
                         smk0_1275, smi_962, smi_989, smi_990, smk1_1274, smk1_1275, sni_1213, \
                         sni_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1561[k] = f_13 * smi_989[k]
                    + f_3 * pc_y[k] * sni_1213[k];

        t_1562[k] = pb_y[k] * smk0_1274[k]
                    - f_12 * pc_y[k] * smk1_1274[k];

        t_1563[k] = pb_y[k] * smk0_1275[k]
                    + f_17 * smi_990[k]
                    - f_12 * pc_y[k] * smk1_1275[k];

        t_1564[k] = f_20 * smi_962[k]
                    + f_3 * pc_z[k] * sni_1214[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, t_1568, pb_y, pc_y, smk0_1277, smk0_1278, \
                         smk0_1280, smi_992, smi_993, smi_994, smk1_1277, smk1_1278, \
                         smk1_1280, sni_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = pb_y[k] * smk0_1277[k]
                    + f_15 * smi_992[k]
                    - f_12 * pc_y[k] * smk1_1277[k];

        t_1566[k] = pb_y[k] * smk0_1278[k]
                    + f_14 * smi_993[k]
                    - f_12 * pc_y[k] * smk1_1278[k];

        t_1567[k] = f_13 * smi_994[k]
                    + f_3 * pc_y[k] * sni_1218[k];

        t_1568[k] = pb_y[k] * smk0_1280[k]
                    - f_12 * pc_y[k] * smk1_1280[k];
    }

#pragma omp simd aligned(t_1569, t_1570, t_1571, t_1572, t_1573, pc_x, smi_1225, smi_1226, \
                         smi_1227, smi_1228, smi_1229, sni_1225, sni_1226, sni_1227, sni_1228, \
                         sni_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1569[k] = f_14 * smi_1225[k]
                    + f_3 * pc_x[k] * sni_1225[k];

        t_1570[k] = f_14 * smi_1226[k]
                    + f_3 * pc_x[k] * sni_1226[k];

        t_1571[k] = f_14 * smi_1227[k]
                    + f_3 * pc_x[k] * sni_1227[k];

        t_1572[k] = f_14 * smi_1228[k]
                    + f_3 * pc_x[k] * sni_1228[k];

        t_1573[k] = f_14 * smi_1229[k]
                    + f_3 * pc_x[k] * sni_1229[k];
    }

#pragma omp simd aligned(t_1574, t_1575, t_1576, t_1577, pc_x, pc_y, pc_z, smi_973, smi_1001, \
                         smi_1230, smi_1231, snh0_918, snh1_918, sni_1225, sni_1230, \
                         sni_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = f_14 * smi_1230[k]
                    + f_3 * pc_x[k] * sni_1230[k];

        t_1575[k] = f_14 * smi_1231[k]
                    + f_3 * pc_x[k] * sni_1231[k];

        t_1576[k] = f_13 * smi_1001[k]
                    + f_1 * snh0_918[k]
                    - f_2 * snh1_918[k]
                    + f_3 * pc_y[k] * sni_1225[k];

        t_1577[k] = f_20 * smi_973[k]
                    + f_3 * pc_z[k] * sni_1225[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pc_y, smi_1003, smi_1004, smi_1005, snh0_920, \
                         snh0_921, snh0_922, snh1_920, snh1_921, snh1_922, sni_1227, sni_1228, \
                         sni_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = f_13 * smi_1003[k]
                    + f_4 * snh0_920[k]
                    - f_5 * snh1_920[k]
                    + f_3 * pc_y[k] * sni_1227[k];

        t_1579[k] = f_13 * smi_1004[k]
                    + f_6 * snh0_921[k]
                    - f_7 * snh1_921[k]
                    + f_3 * pc_y[k] * sni_1228[k];

        t_1580[k] = f_13 * smi_1005[k]
                    + f_8 * snh0_922[k]
                    - f_9 * snh1_922[k]
                    + f_3 * pc_y[k] * sni_1229[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, pb_y, pc_y, smk0_1295, smi_1006, smi_1007, \
                         smk1_1295, snh0_923, snh1_923, sni_1230, \
                         sni_1231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_13 * smi_1006[k]
                    + f_10 * snh0_923[k]
                    - f_11 * snh1_923[k]
                    + f_3 * pc_y[k] * sni_1230[k];

        t_1582[k] = f_13 * smi_1007[k]
                    + f_3 * pc_y[k] * sni_1231[k];

        t_1583[k] = pb_y[k] * smk0_1295[k]
                    - f_12 * pc_y[k] * smk1_1295[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, t_1587, pc_x, pc_y, pc_z, smi_980, smi_1232, \
                         smi_1235, snh0_924, snh0_927, snh1_924, snh1_927, sni_1232, \
                         sni_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = f_14 * smi_1232[k]
                    + f_1 * snh0_924[k]
                    - f_2 * snh1_924[k]
                    + f_3 * pc_x[k] * sni_1232[k];

        t_1585[k] = f_3 * pc_y[k] * sni_1232[k];

        t_1586[k] = f_19 * smi_980[k]
                    + f_3 * pc_z[k] * sni_1232[k];

        t_1587[k] = f_14 * smi_1235[k]
                    + f_4 * snh0_927[k]
                    - f_5 * snh1_927[k]
                    + f_3 * pc_x[k] * sni_1235[k];
    }

#pragma omp simd aligned(t_1588, t_1589, t_1590, pc_x, pc_y, smi_1237, smi_1238, snh0_929, \
                         snh0_930, snh1_929, snh1_930, sni_1234, sni_1237, \
                         sni_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1588[k] = f_3 * pc_y[k] * sni_1234[k];

        t_1589[k] = f_14 * smi_1237[k]
                    + f_4 * snh0_929[k]
                    - f_5 * snh1_929[k]
                    + f_3 * pc_x[k] * sni_1237[k];

        t_1590[k] = f_14 * smi_1238[k]
                    + f_6 * snh0_930[k]
                    - f_7 * snh1_930[k]
                    + f_3 * pc_x[k] * sni_1238[k];
    }

#pragma omp simd aligned(t_1591, t_1592, t_1593, pc_x, pc_y, pc_z, smi_983, smi_1241, \
                         snh0_933, snh1_933, sni_1235, sni_1237, \
                         sni_1241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1591[k] = f_19 * smi_983[k]
                    + f_3 * pc_z[k] * sni_1235[k];

        t_1592[k] = f_3 * pc_y[k] * sni_1237[k];

        t_1593[k] = f_14 * smi_1241[k]
                    + f_6 * snh0_933[k]
                    - f_7 * snh1_933[k]
                    + f_3 * pc_x[k] * sni_1241[k];
    }

#pragma omp simd aligned(t_1594, t_1595, t_1596, pc_x, pc_z, smi_986, smi_1242, smi_1244, \
                         snh0_934, snh0_936, snh1_934, snh1_936, sni_1238, sni_1242, \
                         sni_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1594[k] = f_14 * smi_1242[k]
                    + f_8 * snh0_934[k]
                    - f_9 * snh1_934[k]
                    + f_3 * pc_x[k] * sni_1242[k];

        t_1595[k] = f_19 * smi_986[k]
                    + f_3 * pc_z[k] * sni_1238[k];

        t_1596[k] = f_14 * smi_1244[k]
                    + f_8 * snh0_936[k]
                    - f_9 * snh1_936[k]
                    + f_3 * pc_x[k] * sni_1244[k];
    }

#pragma omp simd aligned(t_1597, t_1598, t_1599, pc_x, pc_y, smi_1246, smi_1247, snh0_938, \
                         snh0_939, snh1_938, snh1_939, sni_1241, sni_1246, \
                         sni_1247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1597[k] = f_3 * pc_y[k] * sni_1241[k];

        t_1598[k] = f_14 * smi_1246[k]
                    + f_8 * snh0_938[k]
                    - f_9 * snh1_938[k]
                    + f_3 * pc_x[k] * sni_1246[k];

        t_1599[k] = f_14 * smi_1247[k]
                    + f_10 * snh0_939[k]
                    - f_11 * snh1_939[k]
                    + f_3 * pc_x[k] * sni_1247[k];
    }

#pragma omp simd aligned(t_1600, t_1601, t_1602, pc_x, pc_z, smi_990, smi_1249, smi_1250, \
                         snh0_941, snh0_942, snh1_941, snh1_942, sni_1242, sni_1249, \
                         sni_1250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1600[k] = f_19 * smi_990[k]
                    + f_3 * pc_z[k] * sni_1242[k];

        t_1601[k] = f_14 * smi_1249[k]
                    + f_10 * snh0_941[k]
                    - f_11 * snh1_941[k]
                    + f_3 * pc_x[k] * sni_1249[k];

        t_1602[k] = f_14 * smi_1250[k]
                    + f_10 * snh0_942[k]
                    - f_11 * snh1_942[k]
                    + f_3 * pc_x[k] * sni_1250[k];
    }

#pragma omp simd aligned(t_1603, t_1604, t_1605, t_1606, pc_x, pc_y, smi_1252, smi_1253, \
                         smi_1254, snh0_944, snh1_944, sni_1246, sni_1252, sni_1253, \
                         sni_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1603[k] = f_3 * pc_y[k] * sni_1246[k];

        t_1604[k] = f_14 * smi_1252[k]
                    + f_10 * snh0_944[k]
                    - f_11 * snh1_944[k]
                    + f_3 * pc_x[k] * sni_1252[k];

        t_1605[k] = f_14 * smi_1253[k]
                    + f_3 * pc_x[k] * sni_1253[k];

        t_1606[k] = f_14 * smi_1254[k]
                    + f_3 * pc_x[k] * sni_1254[k];
    }

#pragma omp simd aligned(t_1607, t_1608, t_1609, t_1610, t_1611, pc_x, smi_1255, smi_1256, \
                         smi_1257, smi_1258, smi_1259, sni_1255, sni_1256, sni_1257, sni_1258, \
                         sni_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1607[k] = f_14 * smi_1255[k]
                    + f_3 * pc_x[k] * sni_1255[k];

        t_1608[k] = f_14 * smi_1256[k]
                    + f_3 * pc_x[k] * sni_1256[k];

        t_1609[k] = f_14 * smi_1257[k]
                    + f_3 * pc_x[k] * sni_1257[k];

        t_1610[k] = f_14 * smi_1258[k]
                    + f_3 * pc_x[k] * sni_1258[k];

        t_1611[k] = f_14 * smi_1259[k]
                    + f_3 * pc_x[k] * sni_1259[k];
    }

#pragma omp simd aligned(t_1612, t_1613, t_1614, t_1615, pc_y, pc_z, smi_1001, snh0_939, \
                         snh0_941, snh0_942, snh1_939, snh1_941, snh1_942, sni_1253, sni_1255, \
                         sni_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1612[k] = f_1 * snh0_939[k]
                    - f_2 * snh1_939[k]
                    + f_3 * pc_y[k] * sni_1253[k];

        t_1613[k] = f_19 * smi_1001[k]
                    + f_3 * pc_z[k] * sni_1253[k];

        t_1614[k] = f_4 * snh0_941[k]
                    - f_5 * snh1_941[k]
                    + f_3 * pc_y[k] * sni_1255[k];

        t_1615[k] = f_6 * snh0_942[k]
                    - f_7 * snh1_942[k]
                    + f_3 * pc_y[k] * sni_1256[k];
    }

#pragma omp simd aligned(t_1616, t_1617, t_1618, t_1619, pc_y, pc_z, smi_1007, snh0_943, \
                         snh0_944, snh1_943, snh1_944, sni_1257, sni_1258, \
                         sni_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1616[k] = f_8 * snh0_943[k]
                    - f_9 * snh1_943[k]
                    + f_3 * pc_y[k] * sni_1257[k];

        t_1617[k] = f_10 * snh0_944[k]
                    - f_11 * snh1_944[k]
                    + f_3 * pc_y[k] * sni_1258[k];

        t_1618[k] = f_3 * pc_y[k] * sni_1259[k];

        t_1619[k] = f_19 * smi_1007[k]
                    + f_1 * snh0_944[k]
                    - f_2 * snh1_944[k]
                    + f_3 * pc_z[k] * sni_1259[k];
    }

#pragma omp simd aligned(t_1620, t_1621, t_1622, t_1623, pb_x, pc_x, pc_y, pc_z, smk0_1620, \
                         smk0_1623, smi_1008, smi_1260, smi_1263, smk1_1620, smk1_1623, \
                         sni_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1620[k] = pb_x[k] * smk0_1620[k]
                    + f_20 * smi_1260[k]
                    - f_12 * pc_x[k] * smk1_1620[k];

        t_1621[k] = f_18 * smi_1008[k]
                    + f_3 * pc_y[k] * sni_1260[k];

        t_1622[k] = f_3 * pc_z[k] * sni_1260[k];

        t_1623[k] = pb_x[k] * smk0_1623[k]
                    + f_17 * smi_1263[k]
                    - f_12 * pc_x[k] * smk1_1623[k];
    }

#pragma omp simd aligned(t_1624, t_1625, t_1626, pb_x, pc_x, pc_y, smk0_1625, smk0_1626, \
                         smi_1010, smi_1265, smi_1266, smk1_1625, smk1_1626, \
                         sni_1262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1624[k] = f_18 * smi_1010[k]
                    + f_3 * pc_y[k] * sni_1262[k];

        t_1625[k] = pb_x[k] * smk0_1625[k]
                    + f_17 * smi_1265[k]
                    - f_12 * pc_x[k] * smk1_1625[k];

        t_1626[k] = pb_x[k] * smk0_1626[k]
                    + f_16 * smi_1266[k]
                    - f_12 * pc_x[k] * smk1_1626[k];
    }

#pragma omp simd aligned(t_1627, t_1628, t_1629, pb_x, pc_x, pc_y, pc_z, smk0_1629, smi_1013, \
                         smi_1269, smk1_1629, sni_1263, sni_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1627[k] = f_3 * pc_z[k] * sni_1263[k];

        t_1628[k] = f_18 * smi_1013[k]
                    + f_3 * pc_y[k] * sni_1265[k];

        t_1629[k] = pb_x[k] * smk0_1629[k]
                    + f_16 * smi_1269[k]
                    - f_12 * pc_x[k] * smk1_1629[k];
    }

#pragma omp simd aligned(t_1630, t_1631, t_1632, pb_x, pc_x, pc_z, smk0_1630, smk0_1632, \
                         smi_1270, smi_1272, smk1_1630, smk1_1632, \
                         sni_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1630[k] = pb_x[k] * smk0_1630[k]
                    + f_15 * smi_1270[k]
                    - f_12 * pc_x[k] * smk1_1630[k];

        t_1631[k] = f_3 * pc_z[k] * sni_1266[k];

        t_1632[k] = pb_x[k] * smk0_1632[k]
                    + f_15 * smi_1272[k]
                    - f_12 * pc_x[k] * smk1_1632[k];
    }

#pragma omp simd aligned(t_1633, t_1634, t_1635, pb_x, pc_x, pc_y, smk0_1634, smk0_1635, \
                         smi_1017, smi_1274, smi_1275, smk1_1634, smk1_1635, \
                         sni_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1633[k] = f_18 * smi_1017[k]
                    + f_3 * pc_y[k] * sni_1269[k];

        t_1634[k] = pb_x[k] * smk0_1634[k]
                    + f_15 * smi_1274[k]
                    - f_12 * pc_x[k] * smk1_1634[k];

        t_1635[k] = pb_x[k] * smk0_1635[k]
                    + f_14 * smi_1275[k]
                    - f_12 * pc_x[k] * smk1_1635[k];
    }

#pragma omp simd aligned(t_1636, t_1637, t_1638, pb_x, pc_x, pc_z, smk0_1637, smk0_1638, \
                         smi_1277, smi_1278, smk1_1637, smk1_1638, \
                         sni_1270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1636[k] = f_3 * pc_z[k] * sni_1270[k];

        t_1637[k] = pb_x[k] * smk0_1637[k]
                    + f_14 * smi_1277[k]
                    - f_12 * pc_x[k] * smk1_1637[k];

        t_1638[k] = pb_x[k] * smk0_1638[k]
                    + f_14 * smi_1278[k]
                    - f_12 * pc_x[k] * smk1_1638[k];
    }

#pragma omp simd aligned(t_1639, t_1640, t_1641, t_1642, pb_x, pc_x, pc_y, smk0_1640, \
                         smi_1022, smi_1280, smi_1281, smi_1282, smk1_1640, sni_1274, \
                         sni_1281, sni_1282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1639[k] = f_18 * smi_1022[k]
                    + f_3 * pc_y[k] * sni_1274[k];

        t_1640[k] = pb_x[k] * smk0_1640[k]
                    + f_14 * smi_1280[k]
                    - f_12 * pc_x[k] * smk1_1640[k];

        t_1641[k] = f_13 * smi_1281[k]
                    + f_3 * pc_x[k] * sni_1281[k];

        t_1642[k] = f_13 * smi_1282[k]
                    + f_3 * pc_x[k] * sni_1282[k];
    }

#pragma omp simd aligned(t_1643, t_1644, t_1645, t_1646, t_1647, pc_x, smi_1283, smi_1284, \
                         smi_1285, smi_1286, smi_1287, sni_1283, sni_1284, sni_1285, sni_1286, \
                         sni_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1643[k] = f_13 * smi_1283[k]
                    + f_3 * pc_x[k] * sni_1283[k];

        t_1644[k] = f_13 * smi_1284[k]
                    + f_3 * pc_x[k] * sni_1284[k];

        t_1645[k] = f_13 * smi_1285[k]
                    + f_3 * pc_x[k] * sni_1285[k];

        t_1646[k] = f_13 * smi_1286[k]
                    + f_3 * pc_x[k] * sni_1286[k];

        t_1647[k] = f_13 * smi_1287[k]
                    + f_3 * pc_x[k] * sni_1287[k];
    }

#pragma omp simd aligned(t_1648, t_1649, t_1650, t_1651, pb_x, pc_x, pc_z, smk0_1648, \
                         smk0_1650, smk0_1651, smk1_1648, smk1_1650, smk1_1651, \
                         sni_1281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1648[k] = pb_x[k] * smk0_1648[k]
                    - f_12 * pc_x[k] * smk1_1648[k];

        t_1649[k] = f_3 * pc_z[k] * sni_1281[k];

        t_1650[k] = pb_x[k] * smk0_1650[k]
                    - f_12 * pc_x[k] * smk1_1650[k];

        t_1651[k] = pb_x[k] * smk0_1651[k]
                    - f_12 * pc_x[k] * smk1_1651[k];
    }

#pragma omp simd aligned(t_1652, t_1653, t_1654, t_1655, pb_x, pc_x, pc_y, smk0_1652, \
                         smk0_1653, smk0_1655, smi_1035, smk1_1652, smk1_1653, smk1_1655, \
                         sni_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1652[k] = pb_x[k] * smk0_1652[k]
                    - f_12 * pc_x[k] * smk1_1652[k];

        t_1653[k] = pb_x[k] * smk0_1653[k]
                    - f_12 * pc_x[k] * smk1_1653[k];

        t_1654[k] = f_18 * smi_1035[k]
                    + f_3 * pc_y[k] * sni_1287[k];

        t_1655[k] = pb_x[k] * smk0_1655[k]
                    - f_12 * pc_x[k] * smk1_1655[k];
    }

#pragma omp simd aligned(t_1656, t_1657, t_1658, t_1659, pb_z, pc_y, pc_z, smk0_1296, \
                         smk0_1299, smi_1008, smi_1036, smk1_1296, smk1_1299, \
                         sni_1288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1656[k] = pb_z[k] * smk0_1296[k]
                    - f_12 * pc_z[k] * smk1_1296[k];

        t_1657[k] = f_19 * smi_1036[k]
                    + f_3 * pc_y[k] * sni_1288[k];

        t_1658[k] = f_13 * smi_1008[k]
                    + f_3 * pc_z[k] * sni_1288[k];

        t_1659[k] = pb_z[k] * smk0_1299[k]
                    - f_12 * pc_z[k] * smk1_1299[k];
    }

#pragma omp simd aligned(t_1660, t_1661, t_1662, pb_x, pb_z, pc_x, pc_y, pc_z, smk0_1302, \
                         smk0_1661, smi_1038, smi_1293, smk1_1302, smk1_1661, \
                         sni_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1660[k] = f_19 * smi_1038[k]
                    + f_3 * pc_y[k] * sni_1290[k];

        t_1661[k] = pb_x[k] * smk0_1661[k]
                    + f_17 * smi_1293[k]
                    - f_12 * pc_x[k] * smk1_1661[k];

        t_1662[k] = pb_z[k] * smk0_1302[k]
                    - f_12 * pc_z[k] * smk1_1302[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smk0,
                                                           const size_t smi, const size_t smk1,
                                                           const size_t sni, const size_t ncols,
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
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smk0_1306 = buffer.data(smk0 + 1306);
    const auto *smk0_1311 = buffer.data(smk0 + 1311);
    const auto *smk0_1665 = buffer.data(smk0 + 1665);
    const auto *smk0_1668 = buffer.data(smk0 + 1668);
    const auto *smk0_1670 = buffer.data(smk0 + 1670);
    const auto *smk0_1673 = buffer.data(smk0 + 1673);
    const auto *smk0_1674 = buffer.data(smk0 + 1674);
    const auto *smk0_1676 = buffer.data(smk0 + 1676);
    const auto *smk0_1684 = buffer.data(smk0 + 1684);
    const auto *smk0_1686 = buffer.data(smk0 + 1686);
    const auto *smk0_1687 = buffer.data(smk0 + 1687);
    const auto *smk0_1688 = buffer.data(smk0 + 1688);
    const auto *smk0_1689 = buffer.data(smk0 + 1689);
    const auto *smk0_1691 = buffer.data(smk0 + 1691);
    const auto *smk0_1692 = buffer.data(smk0 + 1692);
    const auto *smk0_1695 = buffer.data(smk0 + 1695);
    const auto *smk0_1697 = buffer.data(smk0 + 1697);
    const auto *smk0_1698 = buffer.data(smk0 + 1698);
    const auto *smk0_1701 = buffer.data(smk0 + 1701);
    const auto *smk0_1702 = buffer.data(smk0 + 1702);
    const auto *smk0_1704 = buffer.data(smk0 + 1704);
    const auto *smk0_1706 = buffer.data(smk0 + 1706);
    const auto *smk0_1707 = buffer.data(smk0 + 1707);
    const auto *smk0_1709 = buffer.data(smk0 + 1709);
    const auto *smk0_1710 = buffer.data(smk0 + 1710);
    const auto *smk0_1712 = buffer.data(smk0 + 1712);
    const auto *smk0_1720 = buffer.data(smk0 + 1720);
    const auto *smk0_1722 = buffer.data(smk0 + 1722);
    const auto *smk0_1723 = buffer.data(smk0 + 1723);
    const auto *smk0_1724 = buffer.data(smk0 + 1724);
    const auto *smk0_1725 = buffer.data(smk0 + 1725);
    const auto *smk0_1727 = buffer.data(smk0 + 1727);
    const auto *smk0_1728 = buffer.data(smk0 + 1728);
    const auto *smk0_1731 = buffer.data(smk0 + 1731);
    const auto *smk0_1733 = buffer.data(smk0 + 1733);
    const auto *smk0_1734 = buffer.data(smk0 + 1734);
    const auto *smk0_1737 = buffer.data(smk0 + 1737);
    const auto *smk0_1738 = buffer.data(smk0 + 1738);
    const auto *smk0_1740 = buffer.data(smk0 + 1740);
    const auto *smk0_1742 = buffer.data(smk0 + 1742);
    const auto *smk0_1743 = buffer.data(smk0 + 1743);
    const auto *smk0_1745 = buffer.data(smk0 + 1745);
    const auto *smk0_1746 = buffer.data(smk0 + 1746);
    const auto *smk0_1748 = buffer.data(smk0 + 1748);
    const auto *smk0_1756 = buffer.data(smk0 + 1756);
    const auto *smk0_1758 = buffer.data(smk0 + 1758);
    const auto *smk0_1759 = buffer.data(smk0 + 1759);
    const auto *smk0_1760 = buffer.data(smk0 + 1760);
    const auto *smk0_1761 = buffer.data(smk0 + 1761);
    const auto *smk0_1763 = buffer.data(smk0 + 1763);
    const auto *smk0_1764 = buffer.data(smk0 + 1764);
    const auto *smk0_1767 = buffer.data(smk0 + 1767);
    const auto *smk0_1769 = buffer.data(smk0 + 1769);
    const auto *smk0_1770 = buffer.data(smk0 + 1770);
    const auto *smk0_1773 = buffer.data(smk0 + 1773);
    const auto *smk0_1774 = buffer.data(smk0 + 1774);
    const auto *smk0_1776 = buffer.data(smk0 + 1776);
    const auto *smk0_1778 = buffer.data(smk0 + 1778);

    const auto *smi_1011 = buffer.data(smi + 1011);
    const auto *smi_1014 = buffer.data(smi + 1014);
    const auto *smi_1018 = buffer.data(smi + 1018);
    const auto *smi_1029 = buffer.data(smi + 1029);
    const auto *smi_1036 = buffer.data(smi + 1036);
    const auto *smi_1039 = buffer.data(smi + 1039);
    const auto *smi_1041 = buffer.data(smi + 1041);
    const auto *smi_1042 = buffer.data(smi + 1042);
    const auto *smi_1045 = buffer.data(smi + 1045);
    const auto *smi_1046 = buffer.data(smi + 1046);
    const auto *smi_1050 = buffer.data(smi + 1050);
    const auto *smi_1057 = buffer.data(smi + 1057);
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
    const auto *smi_1091 = buffer.data(smi + 1091);
    const auto *smi_1092 = buffer.data(smi + 1092);
    const auto *smi_1094 = buffer.data(smi + 1094);
    const auto *smi_1095 = buffer.data(smi + 1095);
    const auto *smi_1097 = buffer.data(smi + 1097);
    const auto *smi_1098 = buffer.data(smi + 1098);
    const auto *smi_1101 = buffer.data(smi + 1101);
    const auto *smi_1106 = buffer.data(smi + 1106);
    const auto *smi_1119 = buffer.data(smi + 1119);
    const auto *smi_1120 = buffer.data(smi + 1120);
    const auto *smi_1122 = buffer.data(smi + 1122);
    const auto *smi_1125 = buffer.data(smi + 1125);
    const auto *smi_1129 = buffer.data(smi + 1129);
    const auto *smi_1297 = buffer.data(smi + 1297);
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
    const auto *smi_1375 = buffer.data(smi + 1375);
    const auto *smi_1377 = buffer.data(smi + 1377);
    const auto *smi_1378 = buffer.data(smi + 1378);
    const auto *smi_1381 = buffer.data(smi + 1381);
    const auto *smi_1382 = buffer.data(smi + 1382);
    const auto *smi_1384 = buffer.data(smi + 1384);
    const auto *smi_1386 = buffer.data(smi + 1386);

    const auto *smk1_1306 = buffer.data(smk1 + 1306);
    const auto *smk1_1311 = buffer.data(smk1 + 1311);
    const auto *smk1_1665 = buffer.data(smk1 + 1665);
    const auto *smk1_1668 = buffer.data(smk1 + 1668);
    const auto *smk1_1670 = buffer.data(smk1 + 1670);
    const auto *smk1_1673 = buffer.data(smk1 + 1673);
    const auto *smk1_1674 = buffer.data(smk1 + 1674);
    const auto *smk1_1676 = buffer.data(smk1 + 1676);
    const auto *smk1_1684 = buffer.data(smk1 + 1684);
    const auto *smk1_1686 = buffer.data(smk1 + 1686);
    const auto *smk1_1687 = buffer.data(smk1 + 1687);
    const auto *smk1_1688 = buffer.data(smk1 + 1688);
    const auto *smk1_1689 = buffer.data(smk1 + 1689);
    const auto *smk1_1691 = buffer.data(smk1 + 1691);
    const auto *smk1_1692 = buffer.data(smk1 + 1692);
    const auto *smk1_1695 = buffer.data(smk1 + 1695);
    const auto *smk1_1697 = buffer.data(smk1 + 1697);
    const auto *smk1_1698 = buffer.data(smk1 + 1698);
    const auto *smk1_1701 = buffer.data(smk1 + 1701);
    const auto *smk1_1702 = buffer.data(smk1 + 1702);
    const auto *smk1_1704 = buffer.data(smk1 + 1704);
    const auto *smk1_1706 = buffer.data(smk1 + 1706);
    const auto *smk1_1707 = buffer.data(smk1 + 1707);
    const auto *smk1_1709 = buffer.data(smk1 + 1709);
    const auto *smk1_1710 = buffer.data(smk1 + 1710);
    const auto *smk1_1712 = buffer.data(smk1 + 1712);
    const auto *smk1_1720 = buffer.data(smk1 + 1720);
    const auto *smk1_1722 = buffer.data(smk1 + 1722);
    const auto *smk1_1723 = buffer.data(smk1 + 1723);
    const auto *smk1_1724 = buffer.data(smk1 + 1724);
    const auto *smk1_1725 = buffer.data(smk1 + 1725);
    const auto *smk1_1727 = buffer.data(smk1 + 1727);
    const auto *smk1_1728 = buffer.data(smk1 + 1728);
    const auto *smk1_1731 = buffer.data(smk1 + 1731);
    const auto *smk1_1733 = buffer.data(smk1 + 1733);
    const auto *smk1_1734 = buffer.data(smk1 + 1734);
    const auto *smk1_1737 = buffer.data(smk1 + 1737);
    const auto *smk1_1738 = buffer.data(smk1 + 1738);
    const auto *smk1_1740 = buffer.data(smk1 + 1740);
    const auto *smk1_1742 = buffer.data(smk1 + 1742);
    const auto *smk1_1743 = buffer.data(smk1 + 1743);
    const auto *smk1_1745 = buffer.data(smk1 + 1745);
    const auto *smk1_1746 = buffer.data(smk1 + 1746);
    const auto *smk1_1748 = buffer.data(smk1 + 1748);
    const auto *smk1_1756 = buffer.data(smk1 + 1756);
    const auto *smk1_1758 = buffer.data(smk1 + 1758);
    const auto *smk1_1759 = buffer.data(smk1 + 1759);
    const auto *smk1_1760 = buffer.data(smk1 + 1760);
    const auto *smk1_1761 = buffer.data(smk1 + 1761);
    const auto *smk1_1763 = buffer.data(smk1 + 1763);
    const auto *smk1_1764 = buffer.data(smk1 + 1764);
    const auto *smk1_1767 = buffer.data(smk1 + 1767);
    const auto *smk1_1769 = buffer.data(smk1 + 1769);
    const auto *smk1_1770 = buffer.data(smk1 + 1770);
    const auto *smk1_1773 = buffer.data(smk1 + 1773);
    const auto *smk1_1774 = buffer.data(smk1 + 1774);
    const auto *smk1_1776 = buffer.data(smk1 + 1776);
    const auto *smk1_1778 = buffer.data(smk1 + 1778);

    const auto *sni_1291 = buffer.data(sni + 1291);
    const auto *sni_1293 = buffer.data(sni + 1293);
    const auto *sni_1294 = buffer.data(sni + 1294);
    const auto *sni_1297 = buffer.data(sni + 1297);
    const auto *sni_1298 = buffer.data(sni + 1298);
    const auto *sni_1302 = buffer.data(sni + 1302);
    const auto *sni_1309 = buffer.data(sni + 1309);
    const auto *sni_1310 = buffer.data(sni + 1310);
    const auto *sni_1311 = buffer.data(sni + 1311);
    const auto *sni_1312 = buffer.data(sni + 1312);
    const auto *sni_1313 = buffer.data(sni + 1313);
    const auto *sni_1314 = buffer.data(sni + 1314);
    const auto *sni_1315 = buffer.data(sni + 1315);
    const auto *sni_1316 = buffer.data(sni + 1316);
    const auto *sni_1318 = buffer.data(sni + 1318);
    const auto *sni_1319 = buffer.data(sni + 1319);
    const auto *sni_1321 = buffer.data(sni + 1321);
    const auto *sni_1322 = buffer.data(sni + 1322);
    const auto *sni_1325 = buffer.data(sni + 1325);
    const auto *sni_1326 = buffer.data(sni + 1326);
    const auto *sni_1330 = buffer.data(sni + 1330);
    const auto *sni_1337 = buffer.data(sni + 1337);
    const auto *sni_1338 = buffer.data(sni + 1338);
    const auto *sni_1339 = buffer.data(sni + 1339);
    const auto *sni_1340 = buffer.data(sni + 1340);
    const auto *sni_1341 = buffer.data(sni + 1341);
    const auto *sni_1342 = buffer.data(sni + 1342);
    const auto *sni_1343 = buffer.data(sni + 1343);
    const auto *sni_1344 = buffer.data(sni + 1344);
    const auto *sni_1346 = buffer.data(sni + 1346);
    const auto *sni_1347 = buffer.data(sni + 1347);
    const auto *sni_1349 = buffer.data(sni + 1349);
    const auto *sni_1350 = buffer.data(sni + 1350);
    const auto *sni_1353 = buffer.data(sni + 1353);
    const auto *sni_1354 = buffer.data(sni + 1354);
    const auto *sni_1358 = buffer.data(sni + 1358);
    const auto *sni_1365 = buffer.data(sni + 1365);
    const auto *sni_1366 = buffer.data(sni + 1366);
    const auto *sni_1367 = buffer.data(sni + 1367);
    const auto *sni_1368 = buffer.data(sni + 1368);
    const auto *sni_1369 = buffer.data(sni + 1369);
    const auto *sni_1370 = buffer.data(sni + 1370);
    const auto *sni_1371 = buffer.data(sni + 1371);
    const auto *sni_1372 = buffer.data(sni + 1372);
    const auto *sni_1374 = buffer.data(sni + 1374);
    const auto *sni_1375 = buffer.data(sni + 1375);
    const auto *sni_1377 = buffer.data(sni + 1377);
    const auto *sni_1378 = buffer.data(sni + 1378);
    const auto *sni_1381 = buffer.data(sni + 1381);

#pragma omp simd aligned(t_1663, t_1664, t_1665, pb_x, pc_x, pc_y, pc_z, smk0_1665, smi_1011, \
                         smi_1041, smi_1297, smk1_1665, sni_1291, \
                         sni_1293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1663[k] = f_13 * smi_1011[k]
                    + f_3 * pc_z[k] * sni_1291[k];

        t_1664[k] = f_19 * smi_1041[k]
                    + f_3 * pc_y[k] * sni_1293[k];

        t_1665[k] = pb_x[k] * smk0_1665[k]
                    + f_16 * smi_1297[k]
                    - f_12 * pc_x[k] * smk1_1665[k];
    }

#pragma omp simd aligned(t_1666, t_1667, t_1668, pb_x, pb_z, pc_x, pc_z, smk0_1306, smk0_1668, \
                         smi_1014, smi_1300, smk1_1306, smk1_1668, \
                         sni_1294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1666[k] = pb_z[k] * smk0_1306[k]
                    - f_12 * pc_z[k] * smk1_1306[k];

        t_1667[k] = f_13 * smi_1014[k]
                    + f_3 * pc_z[k] * sni_1294[k];

        t_1668[k] = pb_x[k] * smk0_1668[k]
                    + f_15 * smi_1300[k]
                    - f_12 * pc_x[k] * smk1_1668[k];
    }

#pragma omp simd aligned(t_1669, t_1670, t_1671, pb_x, pb_z, pc_x, pc_y, pc_z, smk0_1311, \
                         smk0_1670, smi_1045, smi_1302, smk1_1311, smk1_1670, \
                         sni_1297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1669[k] = f_19 * smi_1045[k]
                    + f_3 * pc_y[k] * sni_1297[k];

        t_1670[k] = pb_x[k] * smk0_1670[k]
                    + f_15 * smi_1302[k]
                    - f_12 * pc_x[k] * smk1_1670[k];

        t_1671[k] = pb_z[k] * smk0_1311[k]
                    - f_12 * pc_z[k] * smk1_1311[k];
    }

#pragma omp simd aligned(t_1672, t_1673, t_1674, pb_x, pc_x, pc_z, smk0_1673, smk0_1674, \
                         smi_1018, smi_1305, smi_1306, smk1_1673, smk1_1674, \
                         sni_1298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1672[k] = f_13 * smi_1018[k]
                    + f_3 * pc_z[k] * sni_1298[k];

        t_1673[k] = pb_x[k] * smk0_1673[k]
                    + f_14 * smi_1305[k]
                    - f_12 * pc_x[k] * smk1_1673[k];

        t_1674[k] = pb_x[k] * smk0_1674[k]
                    + f_14 * smi_1306[k]
                    - f_12 * pc_x[k] * smk1_1674[k];
    }

#pragma omp simd aligned(t_1675, t_1676, t_1677, t_1678, pb_x, pc_x, pc_y, smk0_1676, \
                         smi_1050, smi_1308, smi_1309, smi_1310, smk1_1676, sni_1302, \
                         sni_1309, sni_1310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1675[k] = f_19 * smi_1050[k]
                    + f_3 * pc_y[k] * sni_1302[k];

        t_1676[k] = pb_x[k] * smk0_1676[k]
                    + f_14 * smi_1308[k]
                    - f_12 * pc_x[k] * smk1_1676[k];

        t_1677[k] = f_13 * smi_1309[k]
                    + f_3 * pc_x[k] * sni_1309[k];

        t_1678[k] = f_13 * smi_1310[k]
                    + f_3 * pc_x[k] * sni_1310[k];
    }

#pragma omp simd aligned(t_1679, t_1680, t_1681, t_1682, t_1683, pc_x, smi_1311, smi_1312, \
                         smi_1313, smi_1314, smi_1315, sni_1311, sni_1312, sni_1313, sni_1314, \
                         sni_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1679[k] = f_13 * smi_1311[k]
                    + f_3 * pc_x[k] * sni_1311[k];

        t_1680[k] = f_13 * smi_1312[k]
                    + f_3 * pc_x[k] * sni_1312[k];

        t_1681[k] = f_13 * smi_1313[k]
                    + f_3 * pc_x[k] * sni_1313[k];

        t_1682[k] = f_13 * smi_1314[k]
                    + f_3 * pc_x[k] * sni_1314[k];

        t_1683[k] = f_13 * smi_1315[k]
                    + f_3 * pc_x[k] * sni_1315[k];
    }

#pragma omp simd aligned(t_1684, t_1685, t_1686, t_1687, pb_x, pc_x, pc_z, smk0_1684, \
                         smk0_1686, smk0_1687, smi_1029, smk1_1684, smk1_1686, smk1_1687, \
                         sni_1309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1684[k] = pb_x[k] * smk0_1684[k]
                    - f_12 * pc_x[k] * smk1_1684[k];

        t_1685[k] = f_13 * smi_1029[k]
                    + f_3 * pc_z[k] * sni_1309[k];

        t_1686[k] = pb_x[k] * smk0_1686[k]
                    - f_12 * pc_x[k] * smk1_1686[k];

        t_1687[k] = pb_x[k] * smk0_1687[k]
                    - f_12 * pc_x[k] * smk1_1687[k];
    }

#pragma omp simd aligned(t_1688, t_1689, t_1690, t_1691, pb_x, pc_x, pc_y, smk0_1688, \
                         smk0_1689, smk0_1691, smi_1063, smk1_1688, smk1_1689, smk1_1691, \
                         sni_1315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1688[k] = pb_x[k] * smk0_1688[k]
                    - f_12 * pc_x[k] * smk1_1688[k];

        t_1689[k] = pb_x[k] * smk0_1689[k]
                    - f_12 * pc_x[k] * smk1_1689[k];

        t_1690[k] = f_19 * smi_1063[k]
                    + f_3 * pc_y[k] * sni_1315[k];

        t_1691[k] = pb_x[k] * smk0_1691[k]
                    - f_12 * pc_x[k] * smk1_1691[k];
    }

#pragma omp simd aligned(t_1692, t_1693, t_1694, pb_x, pc_x, pc_y, pc_z, smk0_1692, smi_1036, \
                         smi_1064, smi_1316, smk1_1692, sni_1316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1692[k] = pb_x[k] * smk0_1692[k]
                    + f_20 * smi_1316[k]
                    - f_12 * pc_x[k] * smk1_1692[k];

        t_1693[k] = f_20 * smi_1064[k]
                    + f_3 * pc_y[k] * sni_1316[k];

        t_1694[k] = f_14 * smi_1036[k]
                    + f_3 * pc_z[k] * sni_1316[k];
    }

#pragma omp simd aligned(t_1695, t_1696, t_1697, pb_x, pc_x, pc_y, smk0_1695, smk0_1697, \
                         smi_1066, smi_1319, smi_1321, smk1_1695, smk1_1697, \
                         sni_1318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1695[k] = pb_x[k] * smk0_1695[k]
                    + f_17 * smi_1319[k]
                    - f_12 * pc_x[k] * smk1_1695[k];

        t_1696[k] = f_20 * smi_1066[k]
                    + f_3 * pc_y[k] * sni_1318[k];

        t_1697[k] = pb_x[k] * smk0_1697[k]
                    + f_17 * smi_1321[k]
                    - f_12 * pc_x[k] * smk1_1697[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, pb_x, pc_x, pc_y, pc_z, smk0_1698, smi_1039, \
                         smi_1069, smi_1322, smk1_1698, sni_1319, \
                         sni_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = pb_x[k] * smk0_1698[k]
                    + f_16 * smi_1322[k]
                    - f_12 * pc_x[k] * smk1_1698[k];

        t_1699[k] = f_14 * smi_1039[k]
                    + f_3 * pc_z[k] * sni_1319[k];

        t_1700[k] = f_20 * smi_1069[k]
                    + f_3 * pc_y[k] * sni_1321[k];
    }

#pragma omp simd aligned(t_1701, t_1702, t_1703, pb_x, pc_x, pc_z, smk0_1701, smk0_1702, \
                         smi_1042, smi_1325, smi_1326, smk1_1701, smk1_1702, \
                         sni_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1701[k] = pb_x[k] * smk0_1701[k]
                    + f_16 * smi_1325[k]
                    - f_12 * pc_x[k] * smk1_1701[k];

        t_1702[k] = pb_x[k] * smk0_1702[k]
                    + f_15 * smi_1326[k]
                    - f_12 * pc_x[k] * smk1_1702[k];

        t_1703[k] = f_14 * smi_1042[k]
                    + f_3 * pc_z[k] * sni_1322[k];
    }

#pragma omp simd aligned(t_1704, t_1705, t_1706, pb_x, pc_x, pc_y, smk0_1704, smk0_1706, \
                         smi_1073, smi_1328, smi_1330, smk1_1704, smk1_1706, \
                         sni_1325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1704[k] = pb_x[k] * smk0_1704[k]
                    + f_15 * smi_1328[k]
                    - f_12 * pc_x[k] * smk1_1704[k];

        t_1705[k] = f_20 * smi_1073[k]
                    + f_3 * pc_y[k] * sni_1325[k];

        t_1706[k] = pb_x[k] * smk0_1706[k]
                    + f_15 * smi_1330[k]
                    - f_12 * pc_x[k] * smk1_1706[k];
    }

#pragma omp simd aligned(t_1707, t_1708, t_1709, pb_x, pc_x, pc_z, smk0_1707, smk0_1709, \
                         smi_1046, smi_1331, smi_1333, smk1_1707, smk1_1709, \
                         sni_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1707[k] = pb_x[k] * smk0_1707[k]
                    + f_14 * smi_1331[k]
                    - f_12 * pc_x[k] * smk1_1707[k];

        t_1708[k] = f_14 * smi_1046[k]
                    + f_3 * pc_z[k] * sni_1326[k];

        t_1709[k] = pb_x[k] * smk0_1709[k]
                    + f_14 * smi_1333[k]
                    - f_12 * pc_x[k] * smk1_1709[k];
    }

#pragma omp simd aligned(t_1710, t_1711, t_1712, pb_x, pc_x, pc_y, smk0_1710, smk0_1712, \
                         smi_1078, smi_1334, smi_1336, smk1_1710, smk1_1712, \
                         sni_1330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1710[k] = pb_x[k] * smk0_1710[k]
                    + f_14 * smi_1334[k]
                    - f_12 * pc_x[k] * smk1_1710[k];

        t_1711[k] = f_20 * smi_1078[k]
                    + f_3 * pc_y[k] * sni_1330[k];

        t_1712[k] = pb_x[k] * smk0_1712[k]
                    + f_14 * smi_1336[k]
                    - f_12 * pc_x[k] * smk1_1712[k];
    }

#pragma omp simd aligned(t_1713, t_1714, t_1715, t_1716, t_1717, pc_x, smi_1337, smi_1338, \
                         smi_1339, smi_1340, smi_1341, sni_1337, sni_1338, sni_1339, sni_1340, \
                         sni_1341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1713[k] = f_13 * smi_1337[k]
                    + f_3 * pc_x[k] * sni_1337[k];

        t_1714[k] = f_13 * smi_1338[k]
                    + f_3 * pc_x[k] * sni_1338[k];

        t_1715[k] = f_13 * smi_1339[k]
                    + f_3 * pc_x[k] * sni_1339[k];

        t_1716[k] = f_13 * smi_1340[k]
                    + f_3 * pc_x[k] * sni_1340[k];

        t_1717[k] = f_13 * smi_1341[k]
                    + f_3 * pc_x[k] * sni_1341[k];
    }

#pragma omp simd aligned(t_1718, t_1719, t_1720, t_1721, pb_x, pc_x, pc_z, smk0_1720, \
                         smi_1057, smi_1342, smi_1343, smk1_1720, sni_1337, sni_1342, \
                         sni_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1718[k] = f_13 * smi_1342[k]
                    + f_3 * pc_x[k] * sni_1342[k];

        t_1719[k] = f_13 * smi_1343[k]
                    + f_3 * pc_x[k] * sni_1343[k];

        t_1720[k] = pb_x[k] * smk0_1720[k]
                    - f_12 * pc_x[k] * smk1_1720[k];

        t_1721[k] = f_14 * smi_1057[k]
                    + f_3 * pc_z[k] * sni_1337[k];
    }

#pragma omp simd aligned(t_1722, t_1723, t_1724, t_1725, pb_x, pc_x, smk0_1722, smk0_1723, \
                         smk0_1724, smk0_1725, smk1_1722, smk1_1723, smk1_1724, \
                         smk1_1725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1722[k] = pb_x[k] * smk0_1722[k]
                    - f_12 * pc_x[k] * smk1_1722[k];

        t_1723[k] = pb_x[k] * smk0_1723[k]
                    - f_12 * pc_x[k] * smk1_1723[k];

        t_1724[k] = pb_x[k] * smk0_1724[k]
                    - f_12 * pc_x[k] * smk1_1724[k];

        t_1725[k] = pb_x[k] * smk0_1725[k]
                    - f_12 * pc_x[k] * smk1_1725[k];
    }

#pragma omp simd aligned(t_1726, t_1727, t_1728, t_1729, pb_x, pc_x, pc_y, smk0_1727, \
                         smk0_1728, smi_1091, smi_1092, smi_1344, smk1_1727, smk1_1728, \
                         sni_1343, sni_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1726[k] = f_20 * smi_1091[k]
                    + f_3 * pc_y[k] * sni_1343[k];

        t_1727[k] = pb_x[k] * smk0_1727[k]
                    - f_12 * pc_x[k] * smk1_1727[k];

        t_1728[k] = pb_x[k] * smk0_1728[k]
                    + f_20 * smi_1344[k]
                    - f_12 * pc_x[k] * smk1_1728[k];

        t_1729[k] = f_21 * smi_1092[k]
                    + f_3 * pc_y[k] * sni_1344[k];
    }

#pragma omp simd aligned(t_1730, t_1731, t_1732, pb_x, pc_x, pc_y, pc_z, smk0_1731, smi_1064, \
                         smi_1094, smi_1347, smk1_1731, sni_1344, \
                         sni_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1730[k] = f_15 * smi_1064[k]
                    + f_3 * pc_z[k] * sni_1344[k];

        t_1731[k] = pb_x[k] * smk0_1731[k]
                    + f_17 * smi_1347[k]
                    - f_12 * pc_x[k] * smk1_1731[k];

        t_1732[k] = f_21 * smi_1094[k]
                    + f_3 * pc_y[k] * sni_1346[k];
    }

#pragma omp simd aligned(t_1733, t_1734, t_1735, pb_x, pc_x, pc_z, smk0_1733, smk0_1734, \
                         smi_1067, smi_1349, smi_1350, smk1_1733, smk1_1734, \
                         sni_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1733[k] = pb_x[k] * smk0_1733[k]
                    + f_17 * smi_1349[k]
                    - f_12 * pc_x[k] * smk1_1733[k];

        t_1734[k] = pb_x[k] * smk0_1734[k]
                    + f_16 * smi_1350[k]
                    - f_12 * pc_x[k] * smk1_1734[k];

        t_1735[k] = f_15 * smi_1067[k]
                    + f_3 * pc_z[k] * sni_1347[k];
    }

#pragma omp simd aligned(t_1736, t_1737, t_1738, pb_x, pc_x, pc_y, smk0_1737, smk0_1738, \
                         smi_1097, smi_1353, smi_1354, smk1_1737, smk1_1738, \
                         sni_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1736[k] = f_21 * smi_1097[k]
                    + f_3 * pc_y[k] * sni_1349[k];

        t_1737[k] = pb_x[k] * smk0_1737[k]
                    + f_16 * smi_1353[k]
                    - f_12 * pc_x[k] * smk1_1737[k];

        t_1738[k] = pb_x[k] * smk0_1738[k]
                    + f_15 * smi_1354[k]
                    - f_12 * pc_x[k] * smk1_1738[k];
    }

#pragma omp simd aligned(t_1739, t_1740, t_1741, pb_x, pc_x, pc_y, pc_z, smk0_1740, smi_1070, \
                         smi_1101, smi_1356, smk1_1740, sni_1350, \
                         sni_1353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1739[k] = f_15 * smi_1070[k]
                    + f_3 * pc_z[k] * sni_1350[k];

        t_1740[k] = pb_x[k] * smk0_1740[k]
                    + f_15 * smi_1356[k]
                    - f_12 * pc_x[k] * smk1_1740[k];

        t_1741[k] = f_21 * smi_1101[k]
                    + f_3 * pc_y[k] * sni_1353[k];
    }

#pragma omp simd aligned(t_1742, t_1743, t_1744, pb_x, pc_x, pc_z, smk0_1742, smk0_1743, \
                         smi_1074, smi_1358, smi_1359, smk1_1742, smk1_1743, \
                         sni_1354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1742[k] = pb_x[k] * smk0_1742[k]
                    + f_15 * smi_1358[k]
                    - f_12 * pc_x[k] * smk1_1742[k];

        t_1743[k] = pb_x[k] * smk0_1743[k]
                    + f_14 * smi_1359[k]
                    - f_12 * pc_x[k] * smk1_1743[k];

        t_1744[k] = f_15 * smi_1074[k]
                    + f_3 * pc_z[k] * sni_1354[k];
    }

#pragma omp simd aligned(t_1745, t_1746, t_1747, pb_x, pc_x, pc_y, smk0_1745, smk0_1746, \
                         smi_1106, smi_1361, smi_1362, smk1_1745, smk1_1746, \
                         sni_1358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1745[k] = pb_x[k] * smk0_1745[k]
                    + f_14 * smi_1361[k]
                    - f_12 * pc_x[k] * smk1_1745[k];

        t_1746[k] = pb_x[k] * smk0_1746[k]
                    + f_14 * smi_1362[k]
                    - f_12 * pc_x[k] * smk1_1746[k];

        t_1747[k] = f_21 * smi_1106[k]
                    + f_3 * pc_y[k] * sni_1358[k];
    }

#pragma omp simd aligned(t_1748, t_1749, t_1750, t_1751, pb_x, pc_x, smk0_1748, smi_1364, \
                         smi_1365, smi_1366, smi_1367, smk1_1748, sni_1365, sni_1366, \
                         sni_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1748[k] = pb_x[k] * smk0_1748[k]
                    + f_14 * smi_1364[k]
                    - f_12 * pc_x[k] * smk1_1748[k];

        t_1749[k] = f_13 * smi_1365[k]
                    + f_3 * pc_x[k] * sni_1365[k];

        t_1750[k] = f_13 * smi_1366[k]
                    + f_3 * pc_x[k] * sni_1366[k];

        t_1751[k] = f_13 * smi_1367[k]
                    + f_3 * pc_x[k] * sni_1367[k];
    }

#pragma omp simd aligned(t_1752, t_1753, t_1754, t_1755, pc_x, smi_1368, smi_1369, smi_1370, \
                         smi_1371, sni_1368, sni_1369, sni_1370, \
                         sni_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1752[k] = f_13 * smi_1368[k]
                    + f_3 * pc_x[k] * sni_1368[k];

        t_1753[k] = f_13 * smi_1369[k]
                    + f_3 * pc_x[k] * sni_1369[k];

        t_1754[k] = f_13 * smi_1370[k]
                    + f_3 * pc_x[k] * sni_1370[k];

        t_1755[k] = f_13 * smi_1371[k]
                    + f_3 * pc_x[k] * sni_1371[k];
    }

#pragma omp simd aligned(t_1756, t_1757, t_1758, t_1759, pb_x, pc_x, pc_z, smk0_1756, \
                         smk0_1758, smk0_1759, smi_1085, smk1_1756, smk1_1758, smk1_1759, \
                         sni_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1756[k] = pb_x[k] * smk0_1756[k]
                    - f_12 * pc_x[k] * smk1_1756[k];

        t_1757[k] = f_15 * smi_1085[k]
                    + f_3 * pc_z[k] * sni_1365[k];

        t_1758[k] = pb_x[k] * smk0_1758[k]
                    - f_12 * pc_x[k] * smk1_1758[k];

        t_1759[k] = pb_x[k] * smk0_1759[k]
                    - f_12 * pc_x[k] * smk1_1759[k];
    }

#pragma omp simd aligned(t_1760, t_1761, t_1762, t_1763, pb_x, pc_x, pc_y, smk0_1760, \
                         smk0_1761, smk0_1763, smi_1119, smk1_1760, smk1_1761, smk1_1763, \
                         sni_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1760[k] = pb_x[k] * smk0_1760[k]
                    - f_12 * pc_x[k] * smk1_1760[k];

        t_1761[k] = pb_x[k] * smk0_1761[k]
                    - f_12 * pc_x[k] * smk1_1761[k];

        t_1762[k] = f_21 * smi_1119[k]
                    + f_3 * pc_y[k] * sni_1371[k];

        t_1763[k] = pb_x[k] * smk0_1763[k]
                    - f_12 * pc_x[k] * smk1_1763[k];
    }

#pragma omp simd aligned(t_1764, t_1765, t_1766, pb_x, pc_x, pc_y, pc_z, smk0_1764, smi_1092, \
                         smi_1120, smi_1372, smk1_1764, sni_1372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1764[k] = pb_x[k] * smk0_1764[k]
                    + f_20 * smi_1372[k]
                    - f_12 * pc_x[k] * smk1_1764[k];

        t_1765[k] = f_17 * smi_1120[k]
                    + f_3 * pc_y[k] * sni_1372[k];

        t_1766[k] = f_16 * smi_1092[k]
                    + f_3 * pc_z[k] * sni_1372[k];
    }

#pragma omp simd aligned(t_1767, t_1768, t_1769, pb_x, pc_x, pc_y, smk0_1767, smk0_1769, \
                         smi_1122, smi_1375, smi_1377, smk1_1767, smk1_1769, \
                         sni_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1767[k] = pb_x[k] * smk0_1767[k]
                    + f_17 * smi_1375[k]
                    - f_12 * pc_x[k] * smk1_1767[k];

        t_1768[k] = f_17 * smi_1122[k]
                    + f_3 * pc_y[k] * sni_1374[k];

        t_1769[k] = pb_x[k] * smk0_1769[k]
                    + f_17 * smi_1377[k]
                    - f_12 * pc_x[k] * smk1_1769[k];
    }

#pragma omp simd aligned(t_1770, t_1771, t_1772, pb_x, pc_x, pc_y, pc_z, smk0_1770, smi_1095, \
                         smi_1125, smi_1378, smk1_1770, sni_1375, \
                         sni_1377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1770[k] = pb_x[k] * smk0_1770[k]
                    + f_16 * smi_1378[k]
                    - f_12 * pc_x[k] * smk1_1770[k];

        t_1771[k] = f_16 * smi_1095[k]
                    + f_3 * pc_z[k] * sni_1375[k];

        t_1772[k] = f_17 * smi_1125[k]
                    + f_3 * pc_y[k] * sni_1377[k];
    }

#pragma omp simd aligned(t_1773, t_1774, t_1775, pb_x, pc_x, pc_z, smk0_1773, smk0_1774, \
                         smi_1098, smi_1381, smi_1382, smk1_1773, smk1_1774, \
                         sni_1378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1773[k] = pb_x[k] * smk0_1773[k]
                    + f_16 * smi_1381[k]
                    - f_12 * pc_x[k] * smk1_1773[k];

        t_1774[k] = pb_x[k] * smk0_1774[k]
                    + f_15 * smi_1382[k]
                    - f_12 * pc_x[k] * smk1_1774[k];

        t_1775[k] = f_16 * smi_1098[k]
                    + f_3 * pc_z[k] * sni_1378[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, pb_x, pc_x, pc_y, smk0_1776, smk0_1778, \
                         smi_1129, smi_1384, smi_1386, smk1_1776, smk1_1778, \
                         sni_1381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = pb_x[k] * smk0_1776[k]
                    + f_15 * smi_1384[k]
                    - f_12 * pc_x[k] * smk1_1776[k];

        t_1777[k] = f_17 * smi_1129[k]
                    + f_3 * pc_y[k] * sni_1381[k];

        t_1778[k] = pb_x[k] * smk0_1778[k]
                    + f_15 * smi_1386[k]
                    - f_12 * pc_x[k] * smk1_1778[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smk0,
                                                           const size_t smi, const size_t smk1,
                                                           const size_t sni, const size_t ncols,
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
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smk0_1779 = buffer.data(smk0 + 1779);
    const auto *smk0_1781 = buffer.data(smk0 + 1781);
    const auto *smk0_1782 = buffer.data(smk0 + 1782);
    const auto *smk0_1784 = buffer.data(smk0 + 1784);
    const auto *smk0_1792 = buffer.data(smk0 + 1792);
    const auto *smk0_1794 = buffer.data(smk0 + 1794);
    const auto *smk0_1795 = buffer.data(smk0 + 1795);
    const auto *smk0_1796 = buffer.data(smk0 + 1796);
    const auto *smk0_1797 = buffer.data(smk0 + 1797);
    const auto *smk0_1799 = buffer.data(smk0 + 1799);
    const auto *smk0_1800 = buffer.data(smk0 + 1800);
    const auto *smk0_1803 = buffer.data(smk0 + 1803);
    const auto *smk0_1805 = buffer.data(smk0 + 1805);
    const auto *smk0_1806 = buffer.data(smk0 + 1806);
    const auto *smk0_1809 = buffer.data(smk0 + 1809);
    const auto *smk0_1810 = buffer.data(smk0 + 1810);
    const auto *smk0_1812 = buffer.data(smk0 + 1812);
    const auto *smk0_1814 = buffer.data(smk0 + 1814);
    const auto *smk0_1815 = buffer.data(smk0 + 1815);
    const auto *smk0_1817 = buffer.data(smk0 + 1817);
    const auto *smk0_1818 = buffer.data(smk0 + 1818);
    const auto *smk0_1820 = buffer.data(smk0 + 1820);
    const auto *smk0_1828 = buffer.data(smk0 + 1828);
    const auto *smk0_1830 = buffer.data(smk0 + 1830);
    const auto *smk0_1831 = buffer.data(smk0 + 1831);
    const auto *smk0_1832 = buffer.data(smk0 + 1832);
    const auto *smk0_1833 = buffer.data(smk0 + 1833);
    const auto *smk0_1835 = buffer.data(smk0 + 1835);
    const auto *smk0_1836 = buffer.data(smk0 + 1836);
    const auto *smk0_1839 = buffer.data(smk0 + 1839);
    const auto *smk0_1841 = buffer.data(smk0 + 1841);
    const auto *smk0_1842 = buffer.data(smk0 + 1842);
    const auto *smk0_1845 = buffer.data(smk0 + 1845);
    const auto *smk0_1846 = buffer.data(smk0 + 1846);
    const auto *smk0_1848 = buffer.data(smk0 + 1848);
    const auto *smk0_1850 = buffer.data(smk0 + 1850);
    const auto *smk0_1851 = buffer.data(smk0 + 1851);
    const auto *smk0_1853 = buffer.data(smk0 + 1853);
    const auto *smk0_1854 = buffer.data(smk0 + 1854);
    const auto *smk0_1856 = buffer.data(smk0 + 1856);
    const auto *smk0_1864 = buffer.data(smk0 + 1864);
    const auto *smk0_1866 = buffer.data(smk0 + 1866);
    const auto *smk0_1867 = buffer.data(smk0 + 1867);
    const auto *smk0_1868 = buffer.data(smk0 + 1868);
    const auto *smk0_1869 = buffer.data(smk0 + 1869);
    const auto *smk0_1871 = buffer.data(smk0 + 1871);
    const auto *smk0_1872 = buffer.data(smk0 + 1872);
    const auto *smk0_1875 = buffer.data(smk0 + 1875);
    const auto *smk0_1877 = buffer.data(smk0 + 1877);
    const auto *smk0_1878 = buffer.data(smk0 + 1878);
    const auto *smk0_1881 = buffer.data(smk0 + 1881);
    const auto *smk0_1882 = buffer.data(smk0 + 1882);
    const auto *smk0_1884 = buffer.data(smk0 + 1884);
    const auto *smk0_1886 = buffer.data(smk0 + 1886);
    const auto *smk0_1887 = buffer.data(smk0 + 1887);
    const auto *smk0_1889 = buffer.data(smk0 + 1889);
    const auto *smk0_1890 = buffer.data(smk0 + 1890);
    const auto *smk0_1892 = buffer.data(smk0 + 1892);

    const auto *smi_1102 = buffer.data(smi + 1102);
    const auto *smi_1113 = buffer.data(smi + 1113);
    const auto *smi_1120 = buffer.data(smi + 1120);
    const auto *smi_1123 = buffer.data(smi + 1123);
    const auto *smi_1126 = buffer.data(smi + 1126);
    const auto *smi_1130 = buffer.data(smi + 1130);
    const auto *smi_1134 = buffer.data(smi + 1134);
    const auto *smi_1141 = buffer.data(smi + 1141);
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
    const auto *smi_1175 = buffer.data(smi + 1175);
    const auto *smi_1176 = buffer.data(smi + 1176);
    const auto *smi_1178 = buffer.data(smi + 1178);
    const auto *smi_1179 = buffer.data(smi + 1179);
    const auto *smi_1181 = buffer.data(smi + 1181);
    const auto *smi_1182 = buffer.data(smi + 1182);
    const auto *smi_1185 = buffer.data(smi + 1185);
    const auto *smi_1186 = buffer.data(smi + 1186);
    const auto *smi_1190 = buffer.data(smi + 1190);
    const auto *smi_1203 = buffer.data(smi + 1203);
    const auto *smi_1204 = buffer.data(smi + 1204);
    const auto *smi_1206 = buffer.data(smi + 1206);
    const auto *smi_1209 = buffer.data(smi + 1209);
    const auto *smi_1213 = buffer.data(smi + 1213);
    const auto *smi_1218 = buffer.data(smi + 1218);
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
    const auto *smi_1403 = buffer.data(smi + 1403);
    const auto *smi_1405 = buffer.data(smi + 1405);
    const auto *smi_1406 = buffer.data(smi + 1406);
    const auto *smi_1409 = buffer.data(smi + 1409);
    const auto *smi_1410 = buffer.data(smi + 1410);
    const auto *smi_1412 = buffer.data(smi + 1412);
    const auto *smi_1414 = buffer.data(smi + 1414);
    const auto *smi_1415 = buffer.data(smi + 1415);
    const auto *smi_1417 = buffer.data(smi + 1417);
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

    const auto *smk1_1779 = buffer.data(smk1 + 1779);
    const auto *smk1_1781 = buffer.data(smk1 + 1781);
    const auto *smk1_1782 = buffer.data(smk1 + 1782);
    const auto *smk1_1784 = buffer.data(smk1 + 1784);
    const auto *smk1_1792 = buffer.data(smk1 + 1792);
    const auto *smk1_1794 = buffer.data(smk1 + 1794);
    const auto *smk1_1795 = buffer.data(smk1 + 1795);
    const auto *smk1_1796 = buffer.data(smk1 + 1796);
    const auto *smk1_1797 = buffer.data(smk1 + 1797);
    const auto *smk1_1799 = buffer.data(smk1 + 1799);
    const auto *smk1_1800 = buffer.data(smk1 + 1800);
    const auto *smk1_1803 = buffer.data(smk1 + 1803);
    const auto *smk1_1805 = buffer.data(smk1 + 1805);
    const auto *smk1_1806 = buffer.data(smk1 + 1806);
    const auto *smk1_1809 = buffer.data(smk1 + 1809);
    const auto *smk1_1810 = buffer.data(smk1 + 1810);
    const auto *smk1_1812 = buffer.data(smk1 + 1812);
    const auto *smk1_1814 = buffer.data(smk1 + 1814);
    const auto *smk1_1815 = buffer.data(smk1 + 1815);
    const auto *smk1_1817 = buffer.data(smk1 + 1817);
    const auto *smk1_1818 = buffer.data(smk1 + 1818);
    const auto *smk1_1820 = buffer.data(smk1 + 1820);
    const auto *smk1_1828 = buffer.data(smk1 + 1828);
    const auto *smk1_1830 = buffer.data(smk1 + 1830);
    const auto *smk1_1831 = buffer.data(smk1 + 1831);
    const auto *smk1_1832 = buffer.data(smk1 + 1832);
    const auto *smk1_1833 = buffer.data(smk1 + 1833);
    const auto *smk1_1835 = buffer.data(smk1 + 1835);
    const auto *smk1_1836 = buffer.data(smk1 + 1836);
    const auto *smk1_1839 = buffer.data(smk1 + 1839);
    const auto *smk1_1841 = buffer.data(smk1 + 1841);
    const auto *smk1_1842 = buffer.data(smk1 + 1842);
    const auto *smk1_1845 = buffer.data(smk1 + 1845);
    const auto *smk1_1846 = buffer.data(smk1 + 1846);
    const auto *smk1_1848 = buffer.data(smk1 + 1848);
    const auto *smk1_1850 = buffer.data(smk1 + 1850);
    const auto *smk1_1851 = buffer.data(smk1 + 1851);
    const auto *smk1_1853 = buffer.data(smk1 + 1853);
    const auto *smk1_1854 = buffer.data(smk1 + 1854);
    const auto *smk1_1856 = buffer.data(smk1 + 1856);
    const auto *smk1_1864 = buffer.data(smk1 + 1864);
    const auto *smk1_1866 = buffer.data(smk1 + 1866);
    const auto *smk1_1867 = buffer.data(smk1 + 1867);
    const auto *smk1_1868 = buffer.data(smk1 + 1868);
    const auto *smk1_1869 = buffer.data(smk1 + 1869);
    const auto *smk1_1871 = buffer.data(smk1 + 1871);
    const auto *smk1_1872 = buffer.data(smk1 + 1872);
    const auto *smk1_1875 = buffer.data(smk1 + 1875);
    const auto *smk1_1877 = buffer.data(smk1 + 1877);
    const auto *smk1_1878 = buffer.data(smk1 + 1878);
    const auto *smk1_1881 = buffer.data(smk1 + 1881);
    const auto *smk1_1882 = buffer.data(smk1 + 1882);
    const auto *smk1_1884 = buffer.data(smk1 + 1884);
    const auto *smk1_1886 = buffer.data(smk1 + 1886);
    const auto *smk1_1887 = buffer.data(smk1 + 1887);
    const auto *smk1_1889 = buffer.data(smk1 + 1889);
    const auto *smk1_1890 = buffer.data(smk1 + 1890);
    const auto *smk1_1892 = buffer.data(smk1 + 1892);

    const auto *sni_1382 = buffer.data(sni + 1382);
    const auto *sni_1386 = buffer.data(sni + 1386);
    const auto *sni_1393 = buffer.data(sni + 1393);
    const auto *sni_1394 = buffer.data(sni + 1394);
    const auto *sni_1395 = buffer.data(sni + 1395);
    const auto *sni_1396 = buffer.data(sni + 1396);
    const auto *sni_1397 = buffer.data(sni + 1397);
    const auto *sni_1398 = buffer.data(sni + 1398);
    const auto *sni_1399 = buffer.data(sni + 1399);
    const auto *sni_1400 = buffer.data(sni + 1400);
    const auto *sni_1402 = buffer.data(sni + 1402);
    const auto *sni_1403 = buffer.data(sni + 1403);
    const auto *sni_1405 = buffer.data(sni + 1405);
    const auto *sni_1406 = buffer.data(sni + 1406);
    const auto *sni_1409 = buffer.data(sni + 1409);
    const auto *sni_1410 = buffer.data(sni + 1410);
    const auto *sni_1414 = buffer.data(sni + 1414);
    const auto *sni_1421 = buffer.data(sni + 1421);
    const auto *sni_1422 = buffer.data(sni + 1422);
    const auto *sni_1423 = buffer.data(sni + 1423);
    const auto *sni_1424 = buffer.data(sni + 1424);
    const auto *sni_1425 = buffer.data(sni + 1425);
    const auto *sni_1426 = buffer.data(sni + 1426);
    const auto *sni_1427 = buffer.data(sni + 1427);
    const auto *sni_1428 = buffer.data(sni + 1428);
    const auto *sni_1430 = buffer.data(sni + 1430);
    const auto *sni_1431 = buffer.data(sni + 1431);
    const auto *sni_1433 = buffer.data(sni + 1433);
    const auto *sni_1434 = buffer.data(sni + 1434);
    const auto *sni_1437 = buffer.data(sni + 1437);
    const auto *sni_1438 = buffer.data(sni + 1438);
    const auto *sni_1442 = buffer.data(sni + 1442);
    const auto *sni_1449 = buffer.data(sni + 1449);
    const auto *sni_1450 = buffer.data(sni + 1450);
    const auto *sni_1451 = buffer.data(sni + 1451);
    const auto *sni_1452 = buffer.data(sni + 1452);
    const auto *sni_1453 = buffer.data(sni + 1453);
    const auto *sni_1454 = buffer.data(sni + 1454);
    const auto *sni_1455 = buffer.data(sni + 1455);
    const auto *sni_1456 = buffer.data(sni + 1456);
    const auto *sni_1458 = buffer.data(sni + 1458);
    const auto *sni_1459 = buffer.data(sni + 1459);
    const auto *sni_1461 = buffer.data(sni + 1461);
    const auto *sni_1462 = buffer.data(sni + 1462);
    const auto *sni_1465 = buffer.data(sni + 1465);
    const auto *sni_1466 = buffer.data(sni + 1466);
    const auto *sni_1470 = buffer.data(sni + 1470);
    const auto *sni_1477 = buffer.data(sni + 1477);
    const auto *sni_1478 = buffer.data(sni + 1478);
    const auto *sni_1479 = buffer.data(sni + 1479);

#pragma omp simd aligned(t_1779, t_1780, t_1781, pb_x, pc_x, pc_z, smk0_1779, smk0_1781, \
                         smi_1102, smi_1387, smi_1389, smk1_1779, smk1_1781, \
                         sni_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1779[k] = pb_x[k] * smk0_1779[k]
                    + f_14 * smi_1387[k]
                    - f_12 * pc_x[k] * smk1_1779[k];

        t_1780[k] = f_16 * smi_1102[k]
                    + f_3 * pc_z[k] * sni_1382[k];

        t_1781[k] = pb_x[k] * smk0_1781[k]
                    + f_14 * smi_1389[k]
                    - f_12 * pc_x[k] * smk1_1781[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, pb_x, pc_x, pc_y, smk0_1782, smk0_1784, \
                         smi_1134, smi_1390, smi_1392, smk1_1782, smk1_1784, \
                         sni_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = pb_x[k] * smk0_1782[k]
                    + f_14 * smi_1390[k]
                    - f_12 * pc_x[k] * smk1_1782[k];

        t_1783[k] = f_17 * smi_1134[k]
                    + f_3 * pc_y[k] * sni_1386[k];

        t_1784[k] = pb_x[k] * smk0_1784[k]
                    + f_14 * smi_1392[k]
                    - f_12 * pc_x[k] * smk1_1784[k];
    }

#pragma omp simd aligned(t_1785, t_1786, t_1787, t_1788, t_1789, pc_x, smi_1393, smi_1394, \
                         smi_1395, smi_1396, smi_1397, sni_1393, sni_1394, sni_1395, sni_1396, \
                         sni_1397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1785[k] = f_13 * smi_1393[k]
                    + f_3 * pc_x[k] * sni_1393[k];

        t_1786[k] = f_13 * smi_1394[k]
                    + f_3 * pc_x[k] * sni_1394[k];

        t_1787[k] = f_13 * smi_1395[k]
                    + f_3 * pc_x[k] * sni_1395[k];

        t_1788[k] = f_13 * smi_1396[k]
                    + f_3 * pc_x[k] * sni_1396[k];

        t_1789[k] = f_13 * smi_1397[k]
                    + f_3 * pc_x[k] * sni_1397[k];
    }

#pragma omp simd aligned(t_1790, t_1791, t_1792, t_1793, pb_x, pc_x, pc_z, smk0_1792, \
                         smi_1113, smi_1398, smi_1399, smk1_1792, sni_1393, sni_1398, \
                         sni_1399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1790[k] = f_13 * smi_1398[k]
                    + f_3 * pc_x[k] * sni_1398[k];

        t_1791[k] = f_13 * smi_1399[k]
                    + f_3 * pc_x[k] * sni_1399[k];

        t_1792[k] = pb_x[k] * smk0_1792[k]
                    - f_12 * pc_x[k] * smk1_1792[k];

        t_1793[k] = f_16 * smi_1113[k]
                    + f_3 * pc_z[k] * sni_1393[k];
    }

#pragma omp simd aligned(t_1794, t_1795, t_1796, t_1797, pb_x, pc_x, smk0_1794, smk0_1795, \
                         smk0_1796, smk0_1797, smk1_1794, smk1_1795, smk1_1796, \
                         smk1_1797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1794[k] = pb_x[k] * smk0_1794[k]
                    - f_12 * pc_x[k] * smk1_1794[k];

        t_1795[k] = pb_x[k] * smk0_1795[k]
                    - f_12 * pc_x[k] * smk1_1795[k];

        t_1796[k] = pb_x[k] * smk0_1796[k]
                    - f_12 * pc_x[k] * smk1_1796[k];

        t_1797[k] = pb_x[k] * smk0_1797[k]
                    - f_12 * pc_x[k] * smk1_1797[k];
    }

#pragma omp simd aligned(t_1798, t_1799, t_1800, t_1801, pb_x, pc_x, pc_y, smk0_1799, \
                         smk0_1800, smi_1147, smi_1148, smi_1400, smk1_1799, smk1_1800, \
                         sni_1399, sni_1400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1798[k] = f_17 * smi_1147[k]
                    + f_3 * pc_y[k] * sni_1399[k];

        t_1799[k] = pb_x[k] * smk0_1799[k]
                    - f_12 * pc_x[k] * smk1_1799[k];

        t_1800[k] = pb_x[k] * smk0_1800[k]
                    + f_20 * smi_1400[k]
                    - f_12 * pc_x[k] * smk1_1800[k];

        t_1801[k] = f_16 * smi_1148[k]
                    + f_3 * pc_y[k] * sni_1400[k];
    }

#pragma omp simd aligned(t_1802, t_1803, t_1804, pb_x, pc_x, pc_y, pc_z, smk0_1803, smi_1120, \
                         smi_1150, smi_1403, smk1_1803, sni_1400, \
                         sni_1402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1802[k] = f_17 * smi_1120[k]
                    + f_3 * pc_z[k] * sni_1400[k];

        t_1803[k] = pb_x[k] * smk0_1803[k]
                    + f_17 * smi_1403[k]
                    - f_12 * pc_x[k] * smk1_1803[k];

        t_1804[k] = f_16 * smi_1150[k]
                    + f_3 * pc_y[k] * sni_1402[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, pb_x, pc_x, pc_z, smk0_1805, smk0_1806, \
                         smi_1123, smi_1405, smi_1406, smk1_1805, smk1_1806, \
                         sni_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = pb_x[k] * smk0_1805[k]
                    + f_17 * smi_1405[k]
                    - f_12 * pc_x[k] * smk1_1805[k];

        t_1806[k] = pb_x[k] * smk0_1806[k]
                    + f_16 * smi_1406[k]
                    - f_12 * pc_x[k] * smk1_1806[k];

        t_1807[k] = f_17 * smi_1123[k]
                    + f_3 * pc_z[k] * sni_1403[k];
    }

#pragma omp simd aligned(t_1808, t_1809, t_1810, pb_x, pc_x, pc_y, smk0_1809, smk0_1810, \
                         smi_1153, smi_1409, smi_1410, smk1_1809, smk1_1810, \
                         sni_1405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1808[k] = f_16 * smi_1153[k]
                    + f_3 * pc_y[k] * sni_1405[k];

        t_1809[k] = pb_x[k] * smk0_1809[k]
                    + f_16 * smi_1409[k]
                    - f_12 * pc_x[k] * smk1_1809[k];

        t_1810[k] = pb_x[k] * smk0_1810[k]
                    + f_15 * smi_1410[k]
                    - f_12 * pc_x[k] * smk1_1810[k];
    }

#pragma omp simd aligned(t_1811, t_1812, t_1813, pb_x, pc_x, pc_y, pc_z, smk0_1812, smi_1126, \
                         smi_1157, smi_1412, smk1_1812, sni_1406, \
                         sni_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1811[k] = f_17 * smi_1126[k]
                    + f_3 * pc_z[k] * sni_1406[k];

        t_1812[k] = pb_x[k] * smk0_1812[k]
                    + f_15 * smi_1412[k]
                    - f_12 * pc_x[k] * smk1_1812[k];

        t_1813[k] = f_16 * smi_1157[k]
                    + f_3 * pc_y[k] * sni_1409[k];
    }

#pragma omp simd aligned(t_1814, t_1815, t_1816, pb_x, pc_x, pc_z, smk0_1814, smk0_1815, \
                         smi_1130, smi_1414, smi_1415, smk1_1814, smk1_1815, \
                         sni_1410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1814[k] = pb_x[k] * smk0_1814[k]
                    + f_15 * smi_1414[k]
                    - f_12 * pc_x[k] * smk1_1814[k];

        t_1815[k] = pb_x[k] * smk0_1815[k]
                    + f_14 * smi_1415[k]
                    - f_12 * pc_x[k] * smk1_1815[k];

        t_1816[k] = f_17 * smi_1130[k]
                    + f_3 * pc_z[k] * sni_1410[k];
    }

#pragma omp simd aligned(t_1817, t_1818, t_1819, pb_x, pc_x, pc_y, smk0_1817, smk0_1818, \
                         smi_1162, smi_1417, smi_1418, smk1_1817, smk1_1818, \
                         sni_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1817[k] = pb_x[k] * smk0_1817[k]
                    + f_14 * smi_1417[k]
                    - f_12 * pc_x[k] * smk1_1817[k];

        t_1818[k] = pb_x[k] * smk0_1818[k]
                    + f_14 * smi_1418[k]
                    - f_12 * pc_x[k] * smk1_1818[k];

        t_1819[k] = f_16 * smi_1162[k]
                    + f_3 * pc_y[k] * sni_1414[k];
    }

#pragma omp simd aligned(t_1820, t_1821, t_1822, t_1823, pb_x, pc_x, smk0_1820, smi_1420, \
                         smi_1421, smi_1422, smi_1423, smk1_1820, sni_1421, sni_1422, \
                         sni_1423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1820[k] = pb_x[k] * smk0_1820[k]
                    + f_14 * smi_1420[k]
                    - f_12 * pc_x[k] * smk1_1820[k];

        t_1821[k] = f_13 * smi_1421[k]
                    + f_3 * pc_x[k] * sni_1421[k];

        t_1822[k] = f_13 * smi_1422[k]
                    + f_3 * pc_x[k] * sni_1422[k];

        t_1823[k] = f_13 * smi_1423[k]
                    + f_3 * pc_x[k] * sni_1423[k];
    }

#pragma omp simd aligned(t_1824, t_1825, t_1826, t_1827, pc_x, smi_1424, smi_1425, smi_1426, \
                         smi_1427, sni_1424, sni_1425, sni_1426, \
                         sni_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1824[k] = f_13 * smi_1424[k]
                    + f_3 * pc_x[k] * sni_1424[k];

        t_1825[k] = f_13 * smi_1425[k]
                    + f_3 * pc_x[k] * sni_1425[k];

        t_1826[k] = f_13 * smi_1426[k]
                    + f_3 * pc_x[k] * sni_1426[k];

        t_1827[k] = f_13 * smi_1427[k]
                    + f_3 * pc_x[k] * sni_1427[k];
    }

#pragma omp simd aligned(t_1828, t_1829, t_1830, t_1831, pb_x, pc_x, pc_z, smk0_1828, \
                         smk0_1830, smk0_1831, smi_1141, smk1_1828, smk1_1830, smk1_1831, \
                         sni_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1828[k] = pb_x[k] * smk0_1828[k]
                    - f_12 * pc_x[k] * smk1_1828[k];

        t_1829[k] = f_17 * smi_1141[k]
                    + f_3 * pc_z[k] * sni_1421[k];

        t_1830[k] = pb_x[k] * smk0_1830[k]
                    - f_12 * pc_x[k] * smk1_1830[k];

        t_1831[k] = pb_x[k] * smk0_1831[k]
                    - f_12 * pc_x[k] * smk1_1831[k];
    }

#pragma omp simd aligned(t_1832, t_1833, t_1834, t_1835, pb_x, pc_x, pc_y, smk0_1832, \
                         smk0_1833, smk0_1835, smi_1175, smk1_1832, smk1_1833, smk1_1835, \
                         sni_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1832[k] = pb_x[k] * smk0_1832[k]
                    - f_12 * pc_x[k] * smk1_1832[k];

        t_1833[k] = pb_x[k] * smk0_1833[k]
                    - f_12 * pc_x[k] * smk1_1833[k];

        t_1834[k] = f_16 * smi_1175[k]
                    + f_3 * pc_y[k] * sni_1427[k];

        t_1835[k] = pb_x[k] * smk0_1835[k]
                    - f_12 * pc_x[k] * smk1_1835[k];
    }

#pragma omp simd aligned(t_1836, t_1837, t_1838, pb_x, pc_x, pc_y, pc_z, smk0_1836, smi_1148, \
                         smi_1176, smi_1428, smk1_1836, sni_1428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1836[k] = pb_x[k] * smk0_1836[k]
                    + f_20 * smi_1428[k]
                    - f_12 * pc_x[k] * smk1_1836[k];

        t_1837[k] = f_15 * smi_1176[k]
                    + f_3 * pc_y[k] * sni_1428[k];

        t_1838[k] = f_21 * smi_1148[k]
                    + f_3 * pc_z[k] * sni_1428[k];
    }

#pragma omp simd aligned(t_1839, t_1840, t_1841, pb_x, pc_x, pc_y, smk0_1839, smk0_1841, \
                         smi_1178, smi_1431, smi_1433, smk1_1839, smk1_1841, \
                         sni_1430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1839[k] = pb_x[k] * smk0_1839[k]
                    + f_17 * smi_1431[k]
                    - f_12 * pc_x[k] * smk1_1839[k];

        t_1840[k] = f_15 * smi_1178[k]
                    + f_3 * pc_y[k] * sni_1430[k];

        t_1841[k] = pb_x[k] * smk0_1841[k]
                    + f_17 * smi_1433[k]
                    - f_12 * pc_x[k] * smk1_1841[k];
    }

#pragma omp simd aligned(t_1842, t_1843, t_1844, pb_x, pc_x, pc_y, pc_z, smk0_1842, smi_1151, \
                         smi_1181, smi_1434, smk1_1842, sni_1431, \
                         sni_1433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1842[k] = pb_x[k] * smk0_1842[k]
                    + f_16 * smi_1434[k]
                    - f_12 * pc_x[k] * smk1_1842[k];

        t_1843[k] = f_21 * smi_1151[k]
                    + f_3 * pc_z[k] * sni_1431[k];

        t_1844[k] = f_15 * smi_1181[k]
                    + f_3 * pc_y[k] * sni_1433[k];
    }

#pragma omp simd aligned(t_1845, t_1846, t_1847, pb_x, pc_x, pc_z, smk0_1845, smk0_1846, \
                         smi_1154, smi_1437, smi_1438, smk1_1845, smk1_1846, \
                         sni_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1845[k] = pb_x[k] * smk0_1845[k]
                    + f_16 * smi_1437[k]
                    - f_12 * pc_x[k] * smk1_1845[k];

        t_1846[k] = pb_x[k] * smk0_1846[k]
                    + f_15 * smi_1438[k]
                    - f_12 * pc_x[k] * smk1_1846[k];

        t_1847[k] = f_21 * smi_1154[k]
                    + f_3 * pc_z[k] * sni_1434[k];
    }

#pragma omp simd aligned(t_1848, t_1849, t_1850, pb_x, pc_x, pc_y, smk0_1848, smk0_1850, \
                         smi_1185, smi_1440, smi_1442, smk1_1848, smk1_1850, \
                         sni_1437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1848[k] = pb_x[k] * smk0_1848[k]
                    + f_15 * smi_1440[k]
                    - f_12 * pc_x[k] * smk1_1848[k];

        t_1849[k] = f_15 * smi_1185[k]
                    + f_3 * pc_y[k] * sni_1437[k];

        t_1850[k] = pb_x[k] * smk0_1850[k]
                    + f_15 * smi_1442[k]
                    - f_12 * pc_x[k] * smk1_1850[k];
    }

#pragma omp simd aligned(t_1851, t_1852, t_1853, pb_x, pc_x, pc_z, smk0_1851, smk0_1853, \
                         smi_1158, smi_1443, smi_1445, smk1_1851, smk1_1853, \
                         sni_1438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1851[k] = pb_x[k] * smk0_1851[k]
                    + f_14 * smi_1443[k]
                    - f_12 * pc_x[k] * smk1_1851[k];

        t_1852[k] = f_21 * smi_1158[k]
                    + f_3 * pc_z[k] * sni_1438[k];

        t_1853[k] = pb_x[k] * smk0_1853[k]
                    + f_14 * smi_1445[k]
                    - f_12 * pc_x[k] * smk1_1853[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, pb_x, pc_x, pc_y, smk0_1854, smk0_1856, \
                         smi_1190, smi_1446, smi_1448, smk1_1854, smk1_1856, \
                         sni_1442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = pb_x[k] * smk0_1854[k]
                    + f_14 * smi_1446[k]
                    - f_12 * pc_x[k] * smk1_1854[k];

        t_1855[k] = f_15 * smi_1190[k]
                    + f_3 * pc_y[k] * sni_1442[k];

        t_1856[k] = pb_x[k] * smk0_1856[k]
                    + f_14 * smi_1448[k]
                    - f_12 * pc_x[k] * smk1_1856[k];
    }

#pragma omp simd aligned(t_1857, t_1858, t_1859, t_1860, t_1861, pc_x, smi_1449, smi_1450, \
                         smi_1451, smi_1452, smi_1453, sni_1449, sni_1450, sni_1451, sni_1452, \
                         sni_1453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1857[k] = f_13 * smi_1449[k]
                    + f_3 * pc_x[k] * sni_1449[k];

        t_1858[k] = f_13 * smi_1450[k]
                    + f_3 * pc_x[k] * sni_1450[k];

        t_1859[k] = f_13 * smi_1451[k]
                    + f_3 * pc_x[k] * sni_1451[k];

        t_1860[k] = f_13 * smi_1452[k]
                    + f_3 * pc_x[k] * sni_1452[k];

        t_1861[k] = f_13 * smi_1453[k]
                    + f_3 * pc_x[k] * sni_1453[k];
    }

#pragma omp simd aligned(t_1862, t_1863, t_1864, t_1865, pb_x, pc_x, pc_z, smk0_1864, \
                         smi_1169, smi_1454, smi_1455, smk1_1864, sni_1449, sni_1454, \
                         sni_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1862[k] = f_13 * smi_1454[k]
                    + f_3 * pc_x[k] * sni_1454[k];

        t_1863[k] = f_13 * smi_1455[k]
                    + f_3 * pc_x[k] * sni_1455[k];

        t_1864[k] = pb_x[k] * smk0_1864[k]
                    - f_12 * pc_x[k] * smk1_1864[k];

        t_1865[k] = f_21 * smi_1169[k]
                    + f_3 * pc_z[k] * sni_1449[k];
    }

#pragma omp simd aligned(t_1866, t_1867, t_1868, t_1869, pb_x, pc_x, smk0_1866, smk0_1867, \
                         smk0_1868, smk0_1869, smk1_1866, smk1_1867, smk1_1868, \
                         smk1_1869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1866[k] = pb_x[k] * smk0_1866[k]
                    - f_12 * pc_x[k] * smk1_1866[k];

        t_1867[k] = pb_x[k] * smk0_1867[k]
                    - f_12 * pc_x[k] * smk1_1867[k];

        t_1868[k] = pb_x[k] * smk0_1868[k]
                    - f_12 * pc_x[k] * smk1_1868[k];

        t_1869[k] = pb_x[k] * smk0_1869[k]
                    - f_12 * pc_x[k] * smk1_1869[k];
    }

#pragma omp simd aligned(t_1870, t_1871, t_1872, t_1873, pb_x, pc_x, pc_y, smk0_1871, \
                         smk0_1872, smi_1203, smi_1204, smi_1456, smk1_1871, smk1_1872, \
                         sni_1455, sni_1456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1870[k] = f_15 * smi_1203[k]
                    + f_3 * pc_y[k] * sni_1455[k];

        t_1871[k] = pb_x[k] * smk0_1871[k]
                    - f_12 * pc_x[k] * smk1_1871[k];

        t_1872[k] = pb_x[k] * smk0_1872[k]
                    + f_20 * smi_1456[k]
                    - f_12 * pc_x[k] * smk1_1872[k];

        t_1873[k] = f_14 * smi_1204[k]
                    + f_3 * pc_y[k] * sni_1456[k];
    }

#pragma omp simd aligned(t_1874, t_1875, t_1876, pb_x, pc_x, pc_y, pc_z, smk0_1875, smi_1176, \
                         smi_1206, smi_1459, smk1_1875, sni_1456, \
                         sni_1458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1874[k] = f_20 * smi_1176[k]
                    + f_3 * pc_z[k] * sni_1456[k];

        t_1875[k] = pb_x[k] * smk0_1875[k]
                    + f_17 * smi_1459[k]
                    - f_12 * pc_x[k] * smk1_1875[k];

        t_1876[k] = f_14 * smi_1206[k]
                    + f_3 * pc_y[k] * sni_1458[k];
    }

#pragma omp simd aligned(t_1877, t_1878, t_1879, pb_x, pc_x, pc_z, smk0_1877, smk0_1878, \
                         smi_1179, smi_1461, smi_1462, smk1_1877, smk1_1878, \
                         sni_1459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1877[k] = pb_x[k] * smk0_1877[k]
                    + f_17 * smi_1461[k]
                    - f_12 * pc_x[k] * smk1_1877[k];

        t_1878[k] = pb_x[k] * smk0_1878[k]
                    + f_16 * smi_1462[k]
                    - f_12 * pc_x[k] * smk1_1878[k];

        t_1879[k] = f_20 * smi_1179[k]
                    + f_3 * pc_z[k] * sni_1459[k];
    }

#pragma omp simd aligned(t_1880, t_1881, t_1882, pb_x, pc_x, pc_y, smk0_1881, smk0_1882, \
                         smi_1209, smi_1465, smi_1466, smk1_1881, smk1_1882, \
                         sni_1461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1880[k] = f_14 * smi_1209[k]
                    + f_3 * pc_y[k] * sni_1461[k];

        t_1881[k] = pb_x[k] * smk0_1881[k]
                    + f_16 * smi_1465[k]
                    - f_12 * pc_x[k] * smk1_1881[k];

        t_1882[k] = pb_x[k] * smk0_1882[k]
                    + f_15 * smi_1466[k]
                    - f_12 * pc_x[k] * smk1_1882[k];
    }

#pragma omp simd aligned(t_1883, t_1884, t_1885, pb_x, pc_x, pc_y, pc_z, smk0_1884, smi_1182, \
                         smi_1213, smi_1468, smk1_1884, sni_1462, \
                         sni_1465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1883[k] = f_20 * smi_1182[k]
                    + f_3 * pc_z[k] * sni_1462[k];

        t_1884[k] = pb_x[k] * smk0_1884[k]
                    + f_15 * smi_1468[k]
                    - f_12 * pc_x[k] * smk1_1884[k];

        t_1885[k] = f_14 * smi_1213[k]
                    + f_3 * pc_y[k] * sni_1465[k];
    }

#pragma omp simd aligned(t_1886, t_1887, t_1888, pb_x, pc_x, pc_z, smk0_1886, smk0_1887, \
                         smi_1186, smi_1470, smi_1471, smk1_1886, smk1_1887, \
                         sni_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1886[k] = pb_x[k] * smk0_1886[k]
                    + f_15 * smi_1470[k]
                    - f_12 * pc_x[k] * smk1_1886[k];

        t_1887[k] = pb_x[k] * smk0_1887[k]
                    + f_14 * smi_1471[k]
                    - f_12 * pc_x[k] * smk1_1887[k];

        t_1888[k] = f_20 * smi_1186[k]
                    + f_3 * pc_z[k] * sni_1466[k];
    }

#pragma omp simd aligned(t_1889, t_1890, t_1891, pb_x, pc_x, pc_y, smk0_1889, smk0_1890, \
                         smi_1218, smi_1473, smi_1474, smk1_1889, smk1_1890, \
                         sni_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1889[k] = pb_x[k] * smk0_1889[k]
                    + f_14 * smi_1473[k]
                    - f_12 * pc_x[k] * smk1_1889[k];

        t_1890[k] = pb_x[k] * smk0_1890[k]
                    + f_14 * smi_1474[k]
                    - f_12 * pc_x[k] * smk1_1890[k];

        t_1891[k] = f_14 * smi_1218[k]
                    + f_3 * pc_y[k] * sni_1470[k];
    }

#pragma omp simd aligned(t_1892, t_1893, t_1894, t_1895, pb_x, pc_x, smk0_1892, smi_1476, \
                         smi_1477, smi_1478, smi_1479, smk1_1892, sni_1477, sni_1478, \
                         sni_1479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1892[k] = pb_x[k] * smk0_1892[k]
                    + f_14 * smi_1476[k]
                    - f_12 * pc_x[k] * smk1_1892[k];

        t_1893[k] = f_13 * smi_1477[k]
                    + f_3 * pc_x[k] * sni_1477[k];

        t_1894[k] = f_13 * smi_1478[k]
                    + f_3 * pc_x[k] * sni_1478[k];

        t_1895[k] = f_13 * smi_1479[k]
                    + f_3 * pc_x[k] * sni_1479[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smk0,
                                                           const size_t smi, const size_t smk1,
                                                           const size_t snh0, const size_t snh1,
                                                           const size_t sni, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;

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
    auto *t_1980 = buffer.data(target + 1980);
    auto *t_1981 = buffer.data(target + 1981);
    auto *t_1982 = buffer.data(target + 1982);
    auto *t_1983 = buffer.data(target + 1983);
    auto *t_1984 = buffer.data(target + 1984);
    auto *t_1985 = buffer.data(target + 1985);
    auto *t_1986 = buffer.data(target + 1986);
    auto *t_1987 = buffer.data(target + 1987);
    auto *t_1988 = buffer.data(target + 1988);
    auto *t_1989 = buffer.data(target + 1989);
    auto *t_1990 = buffer.data(target + 1990);
    auto *t_1991 = buffer.data(target + 1991);
    auto *t_1992 = buffer.data(target + 1992);
    auto *t_1993 = buffer.data(target + 1993);
    auto *t_1994 = buffer.data(target + 1994);
    auto *t_1995 = buffer.data(target + 1995);
    auto *t_1996 = buffer.data(target + 1996);
    auto *t_1997 = buffer.data(target + 1997);
    auto *t_1998 = buffer.data(target + 1998);
    auto *t_1999 = buffer.data(target + 1999);
    auto *t_2000 = buffer.data(target + 2000);
    auto *t_2001 = buffer.data(target + 2001);
    auto *t_2002 = buffer.data(target + 2002);
    auto *t_2003 = buffer.data(target + 2003);
    auto *t_2004 = buffer.data(target + 2004);
    auto *t_2005 = buffer.data(target + 2005);
    auto *t_2006 = buffer.data(target + 2006);
    auto *t_2007 = buffer.data(target + 2007);
    auto *t_2008 = buffer.data(target + 2008);
    auto *t_2009 = buffer.data(target + 2009);
    auto *t_2010 = buffer.data(target + 2010);
    auto *t_2011 = buffer.data(target + 2011);
    auto *t_2012 = buffer.data(target + 2012);
    auto *t_2013 = buffer.data(target + 2013);
    auto *t_2014 = buffer.data(target + 2014);
    auto *t_2015 = buffer.data(target + 2015);
    auto *t_2016 = buffer.data(target + 2016);
    auto *t_2017 = buffer.data(target + 2017);
    auto *t_2018 = buffer.data(target + 2018);
    auto *t_2019 = buffer.data(target + 2019);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smk0_1584 = buffer.data(smk0 + 1584);
    const auto *smk0_1589 = buffer.data(smk0 + 1589);
    const auto *smk0_1593 = buffer.data(smk0 + 1593);
    const auto *smk0_1598 = buffer.data(smk0 + 1598);
    const auto *smk0_1604 = buffer.data(smk0 + 1604);
    const auto *smk0_1620 = buffer.data(smk0 + 1620);
    const auto *smk0_1623 = buffer.data(smk0 + 1623);
    const auto *smk0_1900 = buffer.data(smk0 + 1900);
    const auto *smk0_1902 = buffer.data(smk0 + 1902);
    const auto *smk0_1903 = buffer.data(smk0 + 1903);
    const auto *smk0_1904 = buffer.data(smk0 + 1904);
    const auto *smk0_1905 = buffer.data(smk0 + 1905);
    const auto *smk0_1907 = buffer.data(smk0 + 1907);
    const auto *smk0_1911 = buffer.data(smk0 + 1911);
    const auto *smk0_1914 = buffer.data(smk0 + 1914);
    const auto *smk0_1918 = buffer.data(smk0 + 1918);
    const auto *smk0_1920 = buffer.data(smk0 + 1920);
    const auto *smk0_1923 = buffer.data(smk0 + 1923);
    const auto *smk0_1925 = buffer.data(smk0 + 1925);
    const auto *smk0_1926 = buffer.data(smk0 + 1926);
    const auto *smk0_1936 = buffer.data(smk0 + 1936);
    const auto *smk0_1938 = buffer.data(smk0 + 1938);
    const auto *smk0_1939 = buffer.data(smk0 + 1939);
    const auto *smk0_1940 = buffer.data(smk0 + 1940);
    const auto *smk0_1941 = buffer.data(smk0 + 1941);
    const auto *smk0_1943 = buffer.data(smk0 + 1943);
    const auto *smk0_1944 = buffer.data(smk0 + 1944);
    const auto *smk0_1947 = buffer.data(smk0 + 1947);
    const auto *smk0_1949 = buffer.data(smk0 + 1949);
    const auto *smk0_1950 = buffer.data(smk0 + 1950);
    const auto *smk0_1953 = buffer.data(smk0 + 1953);
    const auto *smk0_1954 = buffer.data(smk0 + 1954);
    const auto *smk0_1956 = buffer.data(smk0 + 1956);
    const auto *smk0_1958 = buffer.data(smk0 + 1958);
    const auto *smk0_1959 = buffer.data(smk0 + 1959);
    const auto *smk0_1961 = buffer.data(smk0 + 1961);
    const auto *smk0_1962 = buffer.data(smk0 + 1962);
    const auto *smk0_1964 = buffer.data(smk0 + 1964);
    const auto *smk0_1972 = buffer.data(smk0 + 1972);
    const auto *smk0_1974 = buffer.data(smk0 + 1974);
    const auto *smk0_1975 = buffer.data(smk0 + 1975);
    const auto *smk0_1976 = buffer.data(smk0 + 1976);
    const auto *smk0_1977 = buffer.data(smk0 + 1977);
    const auto *smk0_1979 = buffer.data(smk0 + 1979);

    const auto *smi_1197 = buffer.data(smi + 1197);
    const auto *smi_1204 = buffer.data(smi + 1204);
    const auto *smi_1207 = buffer.data(smi + 1207);
    const auto *smi_1210 = buffer.data(smi + 1210);
    const auto *smi_1214 = buffer.data(smi + 1214);
    const auto *smi_1225 = buffer.data(smi + 1225);
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
    const auto *smi_1259 = buffer.data(smi + 1259);
    const auto *smi_1260 = buffer.data(smi + 1260);
    const auto *smi_1262 = buffer.data(smi + 1262);
    const auto *smi_1265 = buffer.data(smi + 1265);
    const auto *smi_1269 = buffer.data(smi + 1269);
    const auto *smi_1274 = buffer.data(smi + 1274);
    const auto *smi_1281 = buffer.data(smi + 1281);
    const auto *smi_1283 = buffer.data(smi + 1283);
    const auto *smi_1284 = buffer.data(smi + 1284);
    const auto *smi_1285 = buffer.data(smi + 1285);
    const auto *smi_1286 = buffer.data(smi + 1286);
    const auto *smi_1287 = buffer.data(smi + 1287);
    const auto *smi_1288 = buffer.data(smi + 1288);
    const auto *smi_1480 = buffer.data(smi + 1480);
    const auto *smi_1481 = buffer.data(smi + 1481);
    const auto *smi_1482 = buffer.data(smi + 1482);
    const auto *smi_1483 = buffer.data(smi + 1483);
    const auto *smi_1487 = buffer.data(smi + 1487);
    const auto *smi_1490 = buffer.data(smi + 1490);
    const auto *smi_1494 = buffer.data(smi + 1494);
    const auto *smi_1496 = buffer.data(smi + 1496);
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
    const auto *smi_1512 = buffer.data(smi + 1512);
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

    const auto *smk1_1584 = buffer.data(smk1 + 1584);
    const auto *smk1_1589 = buffer.data(smk1 + 1589);
    const auto *smk1_1593 = buffer.data(smk1 + 1593);
    const auto *smk1_1598 = buffer.data(smk1 + 1598);
    const auto *smk1_1604 = buffer.data(smk1 + 1604);
    const auto *smk1_1620 = buffer.data(smk1 + 1620);
    const auto *smk1_1623 = buffer.data(smk1 + 1623);
    const auto *smk1_1900 = buffer.data(smk1 + 1900);
    const auto *smk1_1902 = buffer.data(smk1 + 1902);
    const auto *smk1_1903 = buffer.data(smk1 + 1903);
    const auto *smk1_1904 = buffer.data(smk1 + 1904);
    const auto *smk1_1905 = buffer.data(smk1 + 1905);
    const auto *smk1_1907 = buffer.data(smk1 + 1907);
    const auto *smk1_1911 = buffer.data(smk1 + 1911);
    const auto *smk1_1914 = buffer.data(smk1 + 1914);
    const auto *smk1_1918 = buffer.data(smk1 + 1918);
    const auto *smk1_1920 = buffer.data(smk1 + 1920);
    const auto *smk1_1923 = buffer.data(smk1 + 1923);
    const auto *smk1_1925 = buffer.data(smk1 + 1925);
    const auto *smk1_1926 = buffer.data(smk1 + 1926);
    const auto *smk1_1936 = buffer.data(smk1 + 1936);
    const auto *smk1_1938 = buffer.data(smk1 + 1938);
    const auto *smk1_1939 = buffer.data(smk1 + 1939);
    const auto *smk1_1940 = buffer.data(smk1 + 1940);
    const auto *smk1_1941 = buffer.data(smk1 + 1941);
    const auto *smk1_1943 = buffer.data(smk1 + 1943);
    const auto *smk1_1944 = buffer.data(smk1 + 1944);
    const auto *smk1_1947 = buffer.data(smk1 + 1947);
    const auto *smk1_1949 = buffer.data(smk1 + 1949);
    const auto *smk1_1950 = buffer.data(smk1 + 1950);
    const auto *smk1_1953 = buffer.data(smk1 + 1953);
    const auto *smk1_1954 = buffer.data(smk1 + 1954);
    const auto *smk1_1956 = buffer.data(smk1 + 1956);
    const auto *smk1_1958 = buffer.data(smk1 + 1958);
    const auto *smk1_1959 = buffer.data(smk1 + 1959);
    const auto *smk1_1961 = buffer.data(smk1 + 1961);
    const auto *smk1_1962 = buffer.data(smk1 + 1962);
    const auto *smk1_1964 = buffer.data(smk1 + 1964);
    const auto *smk1_1972 = buffer.data(smk1 + 1972);
    const auto *smk1_1974 = buffer.data(smk1 + 1974);
    const auto *smk1_1975 = buffer.data(smk1 + 1975);
    const auto *smk1_1976 = buffer.data(smk1 + 1976);
    const auto *smk1_1977 = buffer.data(smk1 + 1977);
    const auto *smk1_1979 = buffer.data(smk1 + 1979);

    const auto *snh0_1155 = buffer.data(snh0 + 1155);
    const auto *snh0_1158 = buffer.data(snh0 + 1158);
    const auto *snh0_1160 = buffer.data(snh0 + 1160);
    const auto *snh0_1161 = buffer.data(snh0 + 1161);
    const auto *snh0_1164 = buffer.data(snh0 + 1164);
    const auto *snh0_1165 = buffer.data(snh0 + 1165);
    const auto *snh0_1167 = buffer.data(snh0 + 1167);
    const auto *snh0_1169 = buffer.data(snh0 + 1169);
    const auto *snh0_1170 = buffer.data(snh0 + 1170);
    const auto *snh0_1172 = buffer.data(snh0 + 1172);
    const auto *snh0_1173 = buffer.data(snh0 + 1173);
    const auto *snh0_1174 = buffer.data(snh0 + 1174);
    const auto *snh0_1175 = buffer.data(snh0 + 1175);

    const auto *snh1_1155 = buffer.data(snh1 + 1155);
    const auto *snh1_1158 = buffer.data(snh1 + 1158);
    const auto *snh1_1160 = buffer.data(snh1 + 1160);
    const auto *snh1_1161 = buffer.data(snh1 + 1161);
    const auto *snh1_1164 = buffer.data(snh1 + 1164);
    const auto *snh1_1165 = buffer.data(snh1 + 1165);
    const auto *snh1_1167 = buffer.data(snh1 + 1167);
    const auto *snh1_1169 = buffer.data(snh1 + 1169);
    const auto *snh1_1170 = buffer.data(snh1 + 1170);
    const auto *snh1_1172 = buffer.data(snh1 + 1172);
    const auto *snh1_1173 = buffer.data(snh1 + 1173);
    const auto *snh1_1174 = buffer.data(snh1 + 1174);
    const auto *snh1_1175 = buffer.data(snh1 + 1175);

    const auto *sni_1477 = buffer.data(sni + 1477);
    const auto *sni_1480 = buffer.data(sni + 1480);
    const auto *sni_1481 = buffer.data(sni + 1481);
    const auto *sni_1482 = buffer.data(sni + 1482);
    const auto *sni_1483 = buffer.data(sni + 1483);
    const auto *sni_1484 = buffer.data(sni + 1484);
    const auto *sni_1486 = buffer.data(sni + 1486);
    const auto *sni_1487 = buffer.data(sni + 1487);
    const auto *sni_1489 = buffer.data(sni + 1489);
    const auto *sni_1490 = buffer.data(sni + 1490);
    const auto *sni_1493 = buffer.data(sni + 1493);
    const auto *sni_1494 = buffer.data(sni + 1494);
    const auto *sni_1498 = buffer.data(sni + 1498);
    const auto *sni_1505 = buffer.data(sni + 1505);
    const auto *sni_1506 = buffer.data(sni + 1506);
    const auto *sni_1507 = buffer.data(sni + 1507);
    const auto *sni_1508 = buffer.data(sni + 1508);
    const auto *sni_1509 = buffer.data(sni + 1509);
    const auto *sni_1510 = buffer.data(sni + 1510);
    const auto *sni_1511 = buffer.data(sni + 1511);
    const auto *sni_1512 = buffer.data(sni + 1512);
    const auto *sni_1514 = buffer.data(sni + 1514);
    const auto *sni_1515 = buffer.data(sni + 1515);
    const auto *sni_1517 = buffer.data(sni + 1517);
    const auto *sni_1518 = buffer.data(sni + 1518);
    const auto *sni_1521 = buffer.data(sni + 1521);
    const auto *sni_1522 = buffer.data(sni + 1522);
    const auto *sni_1526 = buffer.data(sni + 1526);
    const auto *sni_1533 = buffer.data(sni + 1533);
    const auto *sni_1534 = buffer.data(sni + 1534);
    const auto *sni_1535 = buffer.data(sni + 1535);
    const auto *sni_1536 = buffer.data(sni + 1536);
    const auto *sni_1537 = buffer.data(sni + 1537);
    const auto *sni_1538 = buffer.data(sni + 1538);
    const auto *sni_1539 = buffer.data(sni + 1539);
    const auto *sni_1540 = buffer.data(sni + 1540);
    const auto *sni_1542 = buffer.data(sni + 1542);
    const auto *sni_1543 = buffer.data(sni + 1543);
    const auto *sni_1545 = buffer.data(sni + 1545);
    const auto *sni_1546 = buffer.data(sni + 1546);
    const auto *sni_1549 = buffer.data(sni + 1549);
    const auto *sni_1550 = buffer.data(sni + 1550);
    const auto *sni_1552 = buffer.data(sni + 1552);
    const auto *sni_1554 = buffer.data(sni + 1554);
    const auto *sni_1555 = buffer.data(sni + 1555);
    const auto *sni_1557 = buffer.data(sni + 1557);
    const auto *sni_1558 = buffer.data(sni + 1558);
    const auto *sni_1560 = buffer.data(sni + 1560);
    const auto *sni_1561 = buffer.data(sni + 1561);
    const auto *sni_1562 = buffer.data(sni + 1562);
    const auto *sni_1563 = buffer.data(sni + 1563);
    const auto *sni_1564 = buffer.data(sni + 1564);
    const auto *sni_1565 = buffer.data(sni + 1565);
    const auto *sni_1566 = buffer.data(sni + 1566);
    const auto *sni_1567 = buffer.data(sni + 1567);
    const auto *sni_1568 = buffer.data(sni + 1568);

#pragma omp simd aligned(t_1896, t_1897, t_1898, t_1899, pc_x, smi_1480, smi_1481, smi_1482, \
                         smi_1483, sni_1480, sni_1481, sni_1482, \
                         sni_1483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1896[k] = f_13 * smi_1480[k]
                    + f_3 * pc_x[k] * sni_1480[k];

        t_1897[k] = f_13 * smi_1481[k]
                    + f_3 * pc_x[k] * sni_1481[k];

        t_1898[k] = f_13 * smi_1482[k]
                    + f_3 * pc_x[k] * sni_1482[k];

        t_1899[k] = f_13 * smi_1483[k]
                    + f_3 * pc_x[k] * sni_1483[k];
    }

#pragma omp simd aligned(t_1900, t_1901, t_1902, t_1903, pb_x, pc_x, pc_z, smk0_1900, \
                         smk0_1902, smk0_1903, smi_1197, smk1_1900, smk1_1902, smk1_1903, \
                         sni_1477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1900[k] = pb_x[k] * smk0_1900[k]
                    - f_12 * pc_x[k] * smk1_1900[k];

        t_1901[k] = f_20 * smi_1197[k]
                    + f_3 * pc_z[k] * sni_1477[k];

        t_1902[k] = pb_x[k] * smk0_1902[k]
                    - f_12 * pc_x[k] * smk1_1902[k];

        t_1903[k] = pb_x[k] * smk0_1903[k]
                    - f_12 * pc_x[k] * smk1_1903[k];
    }

#pragma omp simd aligned(t_1904, t_1905, t_1906, t_1907, pb_x, pc_x, pc_y, smk0_1904, \
                         smk0_1905, smk0_1907, smi_1231, smk1_1904, smk1_1905, smk1_1907, \
                         sni_1483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1904[k] = pb_x[k] * smk0_1904[k]
                    - f_12 * pc_x[k] * smk1_1904[k];

        t_1905[k] = pb_x[k] * smk0_1905[k]
                    - f_12 * pc_x[k] * smk1_1905[k];

        t_1906[k] = f_14 * smi_1231[k]
                    + f_3 * pc_y[k] * sni_1483[k];

        t_1907[k] = pb_x[k] * smk0_1907[k]
                    - f_12 * pc_x[k] * smk1_1907[k];
    }

#pragma omp simd aligned(t_1908, t_1909, t_1910, pb_y, pc_y, pc_z, smk0_1584, smi_1204, \
                         smi_1232, smk1_1584, sni_1484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1908[k] = pb_y[k] * smk0_1584[k]
                    - f_12 * pc_y[k] * smk1_1584[k];

        t_1909[k] = f_13 * smi_1232[k]
                    + f_3 * pc_y[k] * sni_1484[k];

        t_1910[k] = f_19 * smi_1204[k]
                    + f_3 * pc_z[k] * sni_1484[k];
    }

#pragma omp simd aligned(t_1911, t_1912, t_1913, pb_x, pb_y, pc_x, pc_y, smk0_1589, smk0_1911, \
                         smi_1234, smi_1487, smk1_1589, smk1_1911, \
                         sni_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1911[k] = pb_x[k] * smk0_1911[k]
                    + f_17 * smi_1487[k]
                    - f_12 * pc_x[k] * smk1_1911[k];

        t_1912[k] = f_13 * smi_1234[k]
                    + f_3 * pc_y[k] * sni_1486[k];

        t_1913[k] = pb_y[k] * smk0_1589[k]
                    - f_12 * pc_y[k] * smk1_1589[k];
    }

#pragma omp simd aligned(t_1914, t_1915, t_1916, pb_x, pc_x, pc_y, pc_z, smk0_1914, smi_1207, \
                         smi_1237, smi_1490, smk1_1914, sni_1487, \
                         sni_1489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1914[k] = pb_x[k] * smk0_1914[k]
                    + f_16 * smi_1490[k]
                    - f_12 * pc_x[k] * smk1_1914[k];

        t_1915[k] = f_19 * smi_1207[k]
                    + f_3 * pc_z[k] * sni_1487[k];

        t_1916[k] = f_13 * smi_1237[k]
                    + f_3 * pc_y[k] * sni_1489[k];
    }

#pragma omp simd aligned(t_1917, t_1918, t_1919, pb_x, pb_y, pc_x, pc_y, pc_z, smk0_1593, \
                         smk0_1918, smi_1210, smi_1494, smk1_1593, smk1_1918, \
                         sni_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1917[k] = pb_y[k] * smk0_1593[k]
                    - f_12 * pc_y[k] * smk1_1593[k];

        t_1918[k] = pb_x[k] * smk0_1918[k]
                    + f_15 * smi_1494[k]
                    - f_12 * pc_x[k] * smk1_1918[k];

        t_1919[k] = f_19 * smi_1210[k]
                    + f_3 * pc_z[k] * sni_1490[k];
    }

#pragma omp simd aligned(t_1920, t_1921, t_1922, pb_x, pb_y, pc_x, pc_y, smk0_1598, smk0_1920, \
                         smi_1241, smi_1496, smk1_1598, smk1_1920, \
                         sni_1493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1920[k] = pb_x[k] * smk0_1920[k]
                    + f_15 * smi_1496[k]
                    - f_12 * pc_x[k] * smk1_1920[k];

        t_1921[k] = f_13 * smi_1241[k]
                    + f_3 * pc_y[k] * sni_1493[k];

        t_1922[k] = pb_y[k] * smk0_1598[k]
                    - f_12 * pc_y[k] * smk1_1598[k];
    }

#pragma omp simd aligned(t_1923, t_1924, t_1925, pb_x, pc_x, pc_z, smk0_1923, smk0_1925, \
                         smi_1214, smi_1499, smi_1501, smk1_1923, smk1_1925, \
                         sni_1494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1923[k] = pb_x[k] * smk0_1923[k]
                    + f_14 * smi_1499[k]
                    - f_12 * pc_x[k] * smk1_1923[k];

        t_1924[k] = f_19 * smi_1214[k]
                    + f_3 * pc_z[k] * sni_1494[k];

        t_1925[k] = pb_x[k] * smk0_1925[k]
                    + f_14 * smi_1501[k]
                    - f_12 * pc_x[k] * smk1_1925[k];
    }

#pragma omp simd aligned(t_1926, t_1927, t_1928, pb_x, pb_y, pc_x, pc_y, smk0_1604, smk0_1926, \
                         smi_1246, smi_1502, smk1_1604, smk1_1926, \
                         sni_1498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1926[k] = pb_x[k] * smk0_1926[k]
                    + f_14 * smi_1502[k]
                    - f_12 * pc_x[k] * smk1_1926[k];

        t_1927[k] = f_13 * smi_1246[k]
                    + f_3 * pc_y[k] * sni_1498[k];

        t_1928[k] = pb_y[k] * smk0_1604[k]
                    - f_12 * pc_y[k] * smk1_1604[k];
    }

#pragma omp simd aligned(t_1929, t_1930, t_1931, t_1932, t_1933, pc_x, smi_1505, smi_1506, \
                         smi_1507, smi_1508, smi_1509, sni_1505, sni_1506, sni_1507, sni_1508, \
                         sni_1509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1929[k] = f_13 * smi_1505[k]
                    + f_3 * pc_x[k] * sni_1505[k];

        t_1930[k] = f_13 * smi_1506[k]
                    + f_3 * pc_x[k] * sni_1506[k];

        t_1931[k] = f_13 * smi_1507[k]
                    + f_3 * pc_x[k] * sni_1507[k];

        t_1932[k] = f_13 * smi_1508[k]
                    + f_3 * pc_x[k] * sni_1508[k];

        t_1933[k] = f_13 * smi_1509[k]
                    + f_3 * pc_x[k] * sni_1509[k];
    }

#pragma omp simd aligned(t_1934, t_1935, t_1936, t_1937, pb_x, pc_x, pc_z, smk0_1936, \
                         smi_1225, smi_1510, smi_1511, smk1_1936, sni_1505, sni_1510, \
                         sni_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1934[k] = f_13 * smi_1510[k]
                    + f_3 * pc_x[k] * sni_1510[k];

        t_1935[k] = f_13 * smi_1511[k]
                    + f_3 * pc_x[k] * sni_1511[k];

        t_1936[k] = pb_x[k] * smk0_1936[k]
                    - f_12 * pc_x[k] * smk1_1936[k];

        t_1937[k] = f_19 * smi_1225[k]
                    + f_3 * pc_z[k] * sni_1505[k];
    }

#pragma omp simd aligned(t_1938, t_1939, t_1940, t_1941, pb_x, pc_x, smk0_1938, smk0_1939, \
                         smk0_1940, smk0_1941, smk1_1938, smk1_1939, smk1_1940, \
                         smk1_1941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1938[k] = pb_x[k] * smk0_1938[k]
                    - f_12 * pc_x[k] * smk1_1938[k];

        t_1939[k] = pb_x[k] * smk0_1939[k]
                    - f_12 * pc_x[k] * smk1_1939[k];

        t_1940[k] = pb_x[k] * smk0_1940[k]
                    - f_12 * pc_x[k] * smk1_1940[k];

        t_1941[k] = pb_x[k] * smk0_1941[k]
                    - f_12 * pc_x[k] * smk1_1941[k];
    }

#pragma omp simd aligned(t_1942, t_1943, t_1944, t_1945, pb_x, pc_x, pc_y, smk0_1943, \
                         smk0_1944, smi_1259, smi_1512, smk1_1943, smk1_1944, sni_1511, \
                         sni_1512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1942[k] = f_13 * smi_1259[k]
                    + f_3 * pc_y[k] * sni_1511[k];

        t_1943[k] = pb_x[k] * smk0_1943[k]
                    - f_12 * pc_x[k] * smk1_1943[k];

        t_1944[k] = pb_x[k] * smk0_1944[k]
                    + f_20 * smi_1512[k]
                    - f_12 * pc_x[k] * smk1_1944[k];

        t_1945[k] = f_3 * pc_y[k] * sni_1512[k];
    }

#pragma omp simd aligned(t_1946, t_1947, t_1948, pb_x, pc_x, pc_y, pc_z, smk0_1947, smi_1232, \
                         smi_1515, smk1_1947, sni_1512, sni_1514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1946[k] = f_18 * smi_1232[k]
                    + f_3 * pc_z[k] * sni_1512[k];

        t_1947[k] = pb_x[k] * smk0_1947[k]
                    + f_17 * smi_1515[k]
                    - f_12 * pc_x[k] * smk1_1947[k];

        t_1948[k] = f_3 * pc_y[k] * sni_1514[k];
    }

#pragma omp simd aligned(t_1949, t_1950, t_1951, pb_x, pc_x, pc_z, smk0_1949, smk0_1950, \
                         smi_1235, smi_1517, smi_1518, smk1_1949, smk1_1950, \
                         sni_1515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1949[k] = pb_x[k] * smk0_1949[k]
                    + f_17 * smi_1517[k]
                    - f_12 * pc_x[k] * smk1_1949[k];

        t_1950[k] = pb_x[k] * smk0_1950[k]
                    + f_16 * smi_1518[k]
                    - f_12 * pc_x[k] * smk1_1950[k];

        t_1951[k] = f_18 * smi_1235[k]
                    + f_3 * pc_z[k] * sni_1515[k];
    }

#pragma omp simd aligned(t_1952, t_1953, t_1954, pb_x, pc_x, pc_y, smk0_1953, smk0_1954, \
                         smi_1521, smi_1522, smk1_1953, smk1_1954, \
                         sni_1517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1952[k] = f_3 * pc_y[k] * sni_1517[k];

        t_1953[k] = pb_x[k] * smk0_1953[k]
                    + f_16 * smi_1521[k]
                    - f_12 * pc_x[k] * smk1_1953[k];

        t_1954[k] = pb_x[k] * smk0_1954[k]
                    + f_15 * smi_1522[k]
                    - f_12 * pc_x[k] * smk1_1954[k];
    }

#pragma omp simd aligned(t_1955, t_1956, t_1957, pb_x, pc_x, pc_y, pc_z, smk0_1956, smi_1238, \
                         smi_1524, smk1_1956, sni_1518, sni_1521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1955[k] = f_18 * smi_1238[k]
                    + f_3 * pc_z[k] * sni_1518[k];

        t_1956[k] = pb_x[k] * smk0_1956[k]
                    + f_15 * smi_1524[k]
                    - f_12 * pc_x[k] * smk1_1956[k];

        t_1957[k] = f_3 * pc_y[k] * sni_1521[k];
    }

#pragma omp simd aligned(t_1958, t_1959, t_1960, pb_x, pc_x, pc_z, smk0_1958, smk0_1959, \
                         smi_1242, smi_1526, smi_1527, smk1_1958, smk1_1959, \
                         sni_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1958[k] = pb_x[k] * smk0_1958[k]
                    + f_15 * smi_1526[k]
                    - f_12 * pc_x[k] * smk1_1958[k];

        t_1959[k] = pb_x[k] * smk0_1959[k]
                    + f_14 * smi_1527[k]
                    - f_12 * pc_x[k] * smk1_1959[k];

        t_1960[k] = f_18 * smi_1242[k]
                    + f_3 * pc_z[k] * sni_1522[k];
    }

#pragma omp simd aligned(t_1961, t_1962, t_1963, pb_x, pc_x, pc_y, smk0_1961, smk0_1962, \
                         smi_1529, smi_1530, smk1_1961, smk1_1962, \
                         sni_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1961[k] = pb_x[k] * smk0_1961[k]
                    + f_14 * smi_1529[k]
                    - f_12 * pc_x[k] * smk1_1961[k];

        t_1962[k] = pb_x[k] * smk0_1962[k]
                    + f_14 * smi_1530[k]
                    - f_12 * pc_x[k] * smk1_1962[k];

        t_1963[k] = f_3 * pc_y[k] * sni_1526[k];
    }

#pragma omp simd aligned(t_1964, t_1965, t_1966, t_1967, pb_x, pc_x, smk0_1964, smi_1532, \
                         smi_1533, smi_1534, smi_1535, smk1_1964, sni_1533, sni_1534, \
                         sni_1535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1964[k] = pb_x[k] * smk0_1964[k]
                    + f_14 * smi_1532[k]
                    - f_12 * pc_x[k] * smk1_1964[k];

        t_1965[k] = f_13 * smi_1533[k]
                    + f_3 * pc_x[k] * sni_1533[k];

        t_1966[k] = f_13 * smi_1534[k]
                    + f_3 * pc_x[k] * sni_1534[k];

        t_1967[k] = f_13 * smi_1535[k]
                    + f_3 * pc_x[k] * sni_1535[k];
    }

#pragma omp simd aligned(t_1968, t_1969, t_1970, t_1971, pc_x, smi_1536, smi_1537, smi_1538, \
                         smi_1539, sni_1536, sni_1537, sni_1538, \
                         sni_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1968[k] = f_13 * smi_1536[k]
                    + f_3 * pc_x[k] * sni_1536[k];

        t_1969[k] = f_13 * smi_1537[k]
                    + f_3 * pc_x[k] * sni_1537[k];

        t_1970[k] = f_13 * smi_1538[k]
                    + f_3 * pc_x[k] * sni_1538[k];

        t_1971[k] = f_13 * smi_1539[k]
                    + f_3 * pc_x[k] * sni_1539[k];
    }

#pragma omp simd aligned(t_1972, t_1973, t_1974, t_1975, pb_x, pc_x, pc_z, smk0_1972, \
                         smk0_1974, smk0_1975, smi_1253, smk1_1972, smk1_1974, smk1_1975, \
                         sni_1533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1972[k] = pb_x[k] * smk0_1972[k]
                    - f_12 * pc_x[k] * smk1_1972[k];

        t_1973[k] = f_18 * smi_1253[k]
                    + f_3 * pc_z[k] * sni_1533[k];

        t_1974[k] = pb_x[k] * smk0_1974[k]
                    - f_12 * pc_x[k] * smk1_1974[k];

        t_1975[k] = pb_x[k] * smk0_1975[k]
                    - f_12 * pc_x[k] * smk1_1975[k];
    }

#pragma omp simd aligned(t_1976, t_1977, t_1978, t_1979, pb_x, pc_x, pc_y, smk0_1976, \
                         smk0_1977, smk0_1979, smk1_1976, smk1_1977, smk1_1979, \
                         sni_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1976[k] = pb_x[k] * smk0_1976[k]
                    - f_12 * pc_x[k] * smk1_1976[k];

        t_1977[k] = pb_x[k] * smk0_1977[k]
                    - f_12 * pc_x[k] * smk1_1977[k];

        t_1978[k] = f_3 * pc_y[k] * sni_1539[k];

        t_1979[k] = pb_x[k] * smk0_1979[k]
                    - f_12 * pc_x[k] * smk1_1979[k];
    }

#pragma omp simd aligned(t_1980, t_1981, t_1982, t_1983, pc_x, pc_y, pc_z, smi_1260, \
                         snh0_1155, snh0_1158, snh1_1155, snh1_1158, sni_1540, \
                         sni_1543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1980[k] = f_1 * snh0_1155[k]
                    - f_2 * snh1_1155[k]
                    + f_3 * pc_x[k] * sni_1540[k];

        t_1981[k] = f_0 * smi_1260[k]
                    + f_3 * pc_y[k] * sni_1540[k];

        t_1982[k] = f_3 * pc_z[k] * sni_1540[k];

        t_1983[k] = f_4 * snh0_1158[k]
                    - f_5 * snh1_1158[k]
                    + f_3 * pc_x[k] * sni_1543[k];
    }

#pragma omp simd aligned(t_1984, t_1985, t_1986, t_1987, pc_x, pc_y, pc_z, smi_1262, \
                         snh0_1160, snh0_1161, snh1_1160, snh1_1161, sni_1542, sni_1543, \
                         sni_1545, sni_1546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1984[k] = f_0 * smi_1262[k]
                    + f_3 * pc_y[k] * sni_1542[k];

        t_1985[k] = f_4 * snh0_1160[k]
                    - f_5 * snh1_1160[k]
                    + f_3 * pc_x[k] * sni_1545[k];

        t_1986[k] = f_6 * snh0_1161[k]
                    - f_7 * snh1_1161[k]
                    + f_3 * pc_x[k] * sni_1546[k];

        t_1987[k] = f_3 * pc_z[k] * sni_1543[k];
    }

#pragma omp simd aligned(t_1988, t_1989, t_1990, t_1991, pc_x, pc_y, pc_z, smi_1265, \
                         snh0_1164, snh0_1165, snh1_1164, snh1_1165, sni_1545, sni_1546, \
                         sni_1549, sni_1550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1988[k] = f_0 * smi_1265[k]
                    + f_3 * pc_y[k] * sni_1545[k];

        t_1989[k] = f_6 * snh0_1164[k]
                    - f_7 * snh1_1164[k]
                    + f_3 * pc_x[k] * sni_1549[k];

        t_1990[k] = f_8 * snh0_1165[k]
                    - f_9 * snh1_1165[k]
                    + f_3 * pc_x[k] * sni_1550[k];

        t_1991[k] = f_3 * pc_z[k] * sni_1546[k];
    }

#pragma omp simd aligned(t_1992, t_1993, t_1994, pc_x, pc_y, smi_1269, snh0_1167, snh0_1169, \
                         snh1_1167, snh1_1169, sni_1549, sni_1552, \
                         sni_1554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1992[k] = f_8 * snh0_1167[k]
                    - f_9 * snh1_1167[k]
                    + f_3 * pc_x[k] * sni_1552[k];

        t_1993[k] = f_0 * smi_1269[k]
                    + f_3 * pc_y[k] * sni_1549[k];

        t_1994[k] = f_8 * snh0_1169[k]
                    - f_9 * snh1_1169[k]
                    + f_3 * pc_x[k] * sni_1554[k];
    }

#pragma omp simd aligned(t_1995, t_1996, t_1997, t_1998, pc_x, pc_z, snh0_1170, snh0_1172, \
                         snh0_1173, snh1_1170, snh1_1172, snh1_1173, sni_1550, sni_1555, \
                         sni_1557, sni_1558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1995[k] = f_10 * snh0_1170[k]
                    - f_11 * snh1_1170[k]
                    + f_3 * pc_x[k] * sni_1555[k];

        t_1996[k] = f_3 * pc_z[k] * sni_1550[k];

        t_1997[k] = f_10 * snh0_1172[k]
                    - f_11 * snh1_1172[k]
                    + f_3 * pc_x[k] * sni_1557[k];

        t_1998[k] = f_10 * snh0_1173[k]
                    - f_11 * snh1_1173[k]
                    + f_3 * pc_x[k] * sni_1558[k];
    }

#pragma omp simd aligned(t_1999, t_2000, t_2001, t_2002, t_2003, pc_x, pc_y, smi_1274, \
                         snh0_1175, snh1_1175, sni_1554, sni_1560, sni_1561, sni_1562, \
                         sni_1563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1999[k] = f_0 * smi_1274[k]
                    + f_3 * pc_y[k] * sni_1554[k];

        t_2000[k] = f_10 * snh0_1175[k]
                    - f_11 * snh1_1175[k]
                    + f_3 * pc_x[k] * sni_1560[k];

        t_2001[k] = f_3 * pc_x[k] * sni_1561[k];

        t_2002[k] = f_3 * pc_x[k] * sni_1562[k];

        t_2003[k] = f_3 * pc_x[k] * sni_1563[k];
    }

#pragma omp simd aligned(t_2004, t_2005, t_2006, t_2007, t_2008, pc_x, pc_y, smi_1281, \
                         snh0_1170, snh1_1170, sni_1561, sni_1564, sni_1565, sni_1566, \
                         sni_1567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2004[k] = f_3 * pc_x[k] * sni_1564[k];

        t_2005[k] = f_3 * pc_x[k] * sni_1565[k];

        t_2006[k] = f_3 * pc_x[k] * sni_1566[k];

        t_2007[k] = f_3 * pc_x[k] * sni_1567[k];

        t_2008[k] = f_0 * smi_1281[k]
                    + f_1 * snh0_1170[k]
                    - f_2 * snh1_1170[k]
                    + f_3 * pc_y[k] * sni_1561[k];
    }

#pragma omp simd aligned(t_2009, t_2010, t_2011, pc_y, pc_z, smi_1283, smi_1284, snh0_1172, \
                         snh0_1173, snh1_1172, snh1_1173, sni_1561, sni_1563, \
                         sni_1564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2009[k] = f_3 * pc_z[k] * sni_1561[k];

        t_2010[k] = f_0 * smi_1283[k]
                    + f_4 * snh0_1172[k]
                    - f_5 * snh1_1172[k]
                    + f_3 * pc_y[k] * sni_1563[k];

        t_2011[k] = f_0 * smi_1284[k]
                    + f_6 * snh0_1173[k]
                    - f_7 * snh1_1173[k]
                    + f_3 * pc_y[k] * sni_1564[k];
    }

#pragma omp simd aligned(t_2012, t_2013, t_2014, t_2015, pc_y, pc_z, smi_1285, smi_1286, \
                         smi_1287, snh0_1174, snh0_1175, snh1_1174, snh1_1175, sni_1565, \
                         sni_1566, sni_1567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2012[k] = f_0 * smi_1285[k]
                    + f_8 * snh0_1174[k]
                    - f_9 * snh1_1174[k]
                    + f_3 * pc_y[k] * sni_1565[k];

        t_2013[k] = f_0 * smi_1286[k]
                    + f_10 * snh0_1175[k]
                    - f_11 * snh1_1175[k]
                    + f_3 * pc_y[k] * sni_1566[k];

        t_2014[k] = f_0 * smi_1287[k]
                    + f_3 * pc_y[k] * sni_1567[k];

        t_2015[k] = f_1 * snh0_1175[k]
                    - f_2 * snh1_1175[k]
                    + f_3 * pc_z[k] * sni_1567[k];
    }

#pragma omp simd aligned(t_2016, t_2017, t_2018, t_2019, pb_z, pc_y, pc_z, smk0_1620, \
                         smk0_1623, smi_1260, smi_1288, smk1_1620, smk1_1623, \
                         sni_1568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2016[k] = pb_z[k] * smk0_1620[k]
                    - f_12 * pc_z[k] * smk1_1620[k];

        t_2017[k] = f_18 * smi_1288[k]
                    + f_3 * pc_y[k] * sni_1568[k];

        t_2018[k] = f_13 * smi_1260[k]
                    + f_3 * pc_z[k] * sni_1568[k];

        t_2019[k] = pb_z[k] * smk0_1623[k]
                    - f_12 * pc_z[k] * smk1_1623[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece18(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smk0,
                                                           const size_t smi, const size_t smk1,
                                                           const size_t snh0, const size_t snh1,
                                                           const size_t sni, const size_t ncols,
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

    auto *t_2020 = buffer.data(target + 2020);
    auto *t_2021 = buffer.data(target + 2021);
    auto *t_2022 = buffer.data(target + 2022);
    auto *t_2023 = buffer.data(target + 2023);
    auto *t_2024 = buffer.data(target + 2024);
    auto *t_2025 = buffer.data(target + 2025);
    auto *t_2026 = buffer.data(target + 2026);
    auto *t_2027 = buffer.data(target + 2027);
    auto *t_2028 = buffer.data(target + 2028);
    auto *t_2029 = buffer.data(target + 2029);
    auto *t_2030 = buffer.data(target + 2030);
    auto *t_2031 = buffer.data(target + 2031);
    auto *t_2032 = buffer.data(target + 2032);
    auto *t_2033 = buffer.data(target + 2033);
    auto *t_2034 = buffer.data(target + 2034);
    auto *t_2035 = buffer.data(target + 2035);
    auto *t_2036 = buffer.data(target + 2036);
    auto *t_2037 = buffer.data(target + 2037);
    auto *t_2038 = buffer.data(target + 2038);
    auto *t_2039 = buffer.data(target + 2039);
    auto *t_2040 = buffer.data(target + 2040);
    auto *t_2041 = buffer.data(target + 2041);
    auto *t_2042 = buffer.data(target + 2042);
    auto *t_2043 = buffer.data(target + 2043);
    auto *t_2044 = buffer.data(target + 2044);
    auto *t_2045 = buffer.data(target + 2045);
    auto *t_2046 = buffer.data(target + 2046);
    auto *t_2047 = buffer.data(target + 2047);
    auto *t_2048 = buffer.data(target + 2048);
    auto *t_2049 = buffer.data(target + 2049);
    auto *t_2050 = buffer.data(target + 2050);
    auto *t_2051 = buffer.data(target + 2051);
    auto *t_2052 = buffer.data(target + 2052);
    auto *t_2053 = buffer.data(target + 2053);
    auto *t_2054 = buffer.data(target + 2054);
    auto *t_2055 = buffer.data(target + 2055);
    auto *t_2056 = buffer.data(target + 2056);
    auto *t_2057 = buffer.data(target + 2057);
    auto *t_2058 = buffer.data(target + 2058);
    auto *t_2059 = buffer.data(target + 2059);
    auto *t_2060 = buffer.data(target + 2060);
    auto *t_2061 = buffer.data(target + 2061);
    auto *t_2062 = buffer.data(target + 2062);
    auto *t_2063 = buffer.data(target + 2063);
    auto *t_2064 = buffer.data(target + 2064);
    auto *t_2065 = buffer.data(target + 2065);
    auto *t_2066 = buffer.data(target + 2066);
    auto *t_2067 = buffer.data(target + 2067);
    auto *t_2068 = buffer.data(target + 2068);
    auto *t_2069 = buffer.data(target + 2069);
    auto *t_2070 = buffer.data(target + 2070);
    auto *t_2071 = buffer.data(target + 2071);
    auto *t_2072 = buffer.data(target + 2072);
    auto *t_2073 = buffer.data(target + 2073);
    auto *t_2074 = buffer.data(target + 2074);
    auto *t_2075 = buffer.data(target + 2075);
    auto *t_2076 = buffer.data(target + 2076);
    auto *t_2077 = buffer.data(target + 2077);
    auto *t_2078 = buffer.data(target + 2078);
    auto *t_2079 = buffer.data(target + 2079);
    auto *t_2080 = buffer.data(target + 2080);
    auto *t_2081 = buffer.data(target + 2081);
    auto *t_2082 = buffer.data(target + 2082);
    auto *t_2083 = buffer.data(target + 2083);
    auto *t_2084 = buffer.data(target + 2084);
    auto *t_2085 = buffer.data(target + 2085);
    auto *t_2086 = buffer.data(target + 2086);
    auto *t_2087 = buffer.data(target + 2087);
    auto *t_2088 = buffer.data(target + 2088);
    auto *t_2089 = buffer.data(target + 2089);
    auto *t_2090 = buffer.data(target + 2090);
    auto *t_2091 = buffer.data(target + 2091);
    auto *t_2092 = buffer.data(target + 2092);
    auto *t_2093 = buffer.data(target + 2093);
    auto *t_2094 = buffer.data(target + 2094);
    auto *t_2095 = buffer.data(target + 2095);
    auto *t_2096 = buffer.data(target + 2096);
    auto *t_2097 = buffer.data(target + 2097);
    auto *t_2098 = buffer.data(target + 2098);
    auto *t_2099 = buffer.data(target + 2099);
    auto *t_2100 = buffer.data(target + 2100);
    auto *t_2101 = buffer.data(target + 2101);
    auto *t_2102 = buffer.data(target + 2102);
    auto *t_2103 = buffer.data(target + 2103);
    auto *t_2104 = buffer.data(target + 2104);
    auto *t_2105 = buffer.data(target + 2105);
    auto *t_2106 = buffer.data(target + 2106);
    auto *t_2107 = buffer.data(target + 2107);
    auto *t_2108 = buffer.data(target + 2108);
    auto *t_2109 = buffer.data(target + 2109);
    auto *t_2110 = buffer.data(target + 2110);
    auto *t_2111 = buffer.data(target + 2111);
    auto *t_2112 = buffer.data(target + 2112);
    auto *t_2113 = buffer.data(target + 2113);
    auto *t_2114 = buffer.data(target + 2114);
    auto *t_2115 = buffer.data(target + 2115);
    auto *t_2116 = buffer.data(target + 2116);
    auto *t_2117 = buffer.data(target + 2117);
    auto *t_2118 = buffer.data(target + 2118);
    auto *t_2119 = buffer.data(target + 2119);
    auto *t_2120 = buffer.data(target + 2120);
    auto *t_2121 = buffer.data(target + 2121);
    auto *t_2122 = buffer.data(target + 2122);
    auto *t_2123 = buffer.data(target + 2123);
    auto *t_2124 = buffer.data(target + 2124);
    auto *t_2125 = buffer.data(target + 2125);
    auto *t_2126 = buffer.data(target + 2126);
    auto *t_2127 = buffer.data(target + 2127);
    auto *t_2128 = buffer.data(target + 2128);
    auto *t_2129 = buffer.data(target + 2129);
    auto *t_2130 = buffer.data(target + 2130);
    auto *t_2131 = buffer.data(target + 2131);
    auto *t_2132 = buffer.data(target + 2132);
    auto *t_2133 = buffer.data(target + 2133);
    auto *t_2134 = buffer.data(target + 2134);
    auto *t_2135 = buffer.data(target + 2135);
    auto *t_2136 = buffer.data(target + 2136);
    auto *t_2137 = buffer.data(target + 2137);
    auto *t_2138 = buffer.data(target + 2138);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smk0_1626 = buffer.data(smk0 + 1626);
    const auto *smk0_1630 = buffer.data(smk0 + 1630);
    const auto *smk0_1635 = buffer.data(smk0 + 1635);
    const auto *smk0_1648 = buffer.data(smk0 + 1648);
    const auto *smk0_1650 = buffer.data(smk0 + 1650);
    const auto *smk0_1651 = buffer.data(smk0 + 1651);
    const auto *smk0_1652 = buffer.data(smk0 + 1652);
    const auto *smk0_1653 = buffer.data(smk0 + 1653);

    const auto *smi_1263 = buffer.data(smi + 1263);
    const auto *smi_1266 = buffer.data(smi + 1266);
    const auto *smi_1270 = buffer.data(smi + 1270);
    const auto *smi_1281 = buffer.data(smi + 1281);
    const auto *smi_1282 = buffer.data(smi + 1282);
    const auto *smi_1283 = buffer.data(smi + 1283);
    const auto *smi_1284 = buffer.data(smi + 1284);
    const auto *smi_1285 = buffer.data(smi + 1285);
    const auto *smi_1287 = buffer.data(smi + 1287);
    const auto *smi_1288 = buffer.data(smi + 1288);
    const auto *smi_1290 = buffer.data(smi + 1290);
    const auto *smi_1291 = buffer.data(smi + 1291);
    const auto *smi_1293 = buffer.data(smi + 1293);
    const auto *smi_1294 = buffer.data(smi + 1294);
    const auto *smi_1297 = buffer.data(smi + 1297);
    const auto *smi_1298 = buffer.data(smi + 1298);
    const auto *smi_1302 = buffer.data(smi + 1302);
    const auto *smi_1309 = buffer.data(smi + 1309);
    const auto *smi_1315 = buffer.data(smi + 1315);
    const auto *smi_1316 = buffer.data(smi + 1316);
    const auto *smi_1318 = buffer.data(smi + 1318);
    const auto *smi_1319 = buffer.data(smi + 1319);
    const auto *smi_1321 = buffer.data(smi + 1321);
    const auto *smi_1322 = buffer.data(smi + 1322);
    const auto *smi_1325 = buffer.data(smi + 1325);
    const auto *smi_1326 = buffer.data(smi + 1326);
    const auto *smi_1330 = buffer.data(smi + 1330);
    const auto *smi_1337 = buffer.data(smi + 1337);
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
    const auto *smi_1358 = buffer.data(smi + 1358);
    const auto *smi_1365 = buffer.data(smi + 1365);
    const auto *smi_1367 = buffer.data(smi + 1367);
    const auto *smi_1368 = buffer.data(smi + 1368);
    const auto *smi_1369 = buffer.data(smi + 1369);
    const auto *smi_1370 = buffer.data(smi + 1370);
    const auto *smi_1371 = buffer.data(smi + 1371);
    const auto *smi_1372 = buffer.data(smi + 1372);
    const auto *smi_1374 = buffer.data(smi + 1374);
    const auto *smi_1377 = buffer.data(smi + 1377);
    const auto *smi_1381 = buffer.data(smi + 1381);

    const auto *smk1_1626 = buffer.data(smk1 + 1626);
    const auto *smk1_1630 = buffer.data(smk1 + 1630);
    const auto *smk1_1635 = buffer.data(smk1 + 1635);
    const auto *smk1_1648 = buffer.data(smk1 + 1648);
    const auto *smk1_1650 = buffer.data(smk1 + 1650);
    const auto *smk1_1651 = buffer.data(smk1 + 1651);
    const auto *smk1_1652 = buffer.data(smk1 + 1652);
    const auto *smk1_1653 = buffer.data(smk1 + 1653);

    const auto *snh0_1181 = buffer.data(snh0 + 1181);
    const auto *snh0_1185 = buffer.data(snh0 + 1185);
    const auto *snh0_1188 = buffer.data(snh0 + 1188);
    const auto *snh0_1190 = buffer.data(snh0 + 1190);
    const auto *snh0_1193 = buffer.data(snh0 + 1193);
    const auto *snh0_1194 = buffer.data(snh0 + 1194);
    const auto *snh0_1196 = buffer.data(snh0 + 1196);
    const auto *snh0_1197 = buffer.data(snh0 + 1197);
    const auto *snh0_1200 = buffer.data(snh0 + 1200);
    const auto *snh0_1202 = buffer.data(snh0 + 1202);
    const auto *snh0_1203 = buffer.data(snh0 + 1203);
    const auto *snh0_1206 = buffer.data(snh0 + 1206);
    const auto *snh0_1207 = buffer.data(snh0 + 1207);
    const auto *snh0_1209 = buffer.data(snh0 + 1209);
    const auto *snh0_1211 = buffer.data(snh0 + 1211);
    const auto *snh0_1212 = buffer.data(snh0 + 1212);
    const auto *snh0_1214 = buffer.data(snh0 + 1214);
    const auto *snh0_1215 = buffer.data(snh0 + 1215);
    const auto *snh0_1216 = buffer.data(snh0 + 1216);
    const auto *snh0_1217 = buffer.data(snh0 + 1217);
    const auto *snh0_1218 = buffer.data(snh0 + 1218);
    const auto *snh0_1221 = buffer.data(snh0 + 1221);
    const auto *snh0_1223 = buffer.data(snh0 + 1223);
    const auto *snh0_1224 = buffer.data(snh0 + 1224);
    const auto *snh0_1227 = buffer.data(snh0 + 1227);
    const auto *snh0_1228 = buffer.data(snh0 + 1228);
    const auto *snh0_1230 = buffer.data(snh0 + 1230);
    const auto *snh0_1232 = buffer.data(snh0 + 1232);
    const auto *snh0_1233 = buffer.data(snh0 + 1233);
    const auto *snh0_1235 = buffer.data(snh0 + 1235);
    const auto *snh0_1236 = buffer.data(snh0 + 1236);
    const auto *snh0_1237 = buffer.data(snh0 + 1237);
    const auto *snh0_1238 = buffer.data(snh0 + 1238);
    const auto *snh0_1239 = buffer.data(snh0 + 1239);
    const auto *snh0_1242 = buffer.data(snh0 + 1242);
    const auto *snh0_1244 = buffer.data(snh0 + 1244);
    const auto *snh0_1245 = buffer.data(snh0 + 1245);
    const auto *snh0_1248 = buffer.data(snh0 + 1248);
    const auto *snh0_1249 = buffer.data(snh0 + 1249);
    const auto *snh0_1251 = buffer.data(snh0 + 1251);
    const auto *snh0_1253 = buffer.data(snh0 + 1253);

    const auto *snh1_1181 = buffer.data(snh1 + 1181);
    const auto *snh1_1185 = buffer.data(snh1 + 1185);
    const auto *snh1_1188 = buffer.data(snh1 + 1188);
    const auto *snh1_1190 = buffer.data(snh1 + 1190);
    const auto *snh1_1193 = buffer.data(snh1 + 1193);
    const auto *snh1_1194 = buffer.data(snh1 + 1194);
    const auto *snh1_1196 = buffer.data(snh1 + 1196);
    const auto *snh1_1197 = buffer.data(snh1 + 1197);
    const auto *snh1_1200 = buffer.data(snh1 + 1200);
    const auto *snh1_1202 = buffer.data(snh1 + 1202);
    const auto *snh1_1203 = buffer.data(snh1 + 1203);
    const auto *snh1_1206 = buffer.data(snh1 + 1206);
    const auto *snh1_1207 = buffer.data(snh1 + 1207);
    const auto *snh1_1209 = buffer.data(snh1 + 1209);
    const auto *snh1_1211 = buffer.data(snh1 + 1211);
    const auto *snh1_1212 = buffer.data(snh1 + 1212);
    const auto *snh1_1214 = buffer.data(snh1 + 1214);
    const auto *snh1_1215 = buffer.data(snh1 + 1215);
    const auto *snh1_1216 = buffer.data(snh1 + 1216);
    const auto *snh1_1217 = buffer.data(snh1 + 1217);
    const auto *snh1_1218 = buffer.data(snh1 + 1218);
    const auto *snh1_1221 = buffer.data(snh1 + 1221);
    const auto *snh1_1223 = buffer.data(snh1 + 1223);
    const auto *snh1_1224 = buffer.data(snh1 + 1224);
    const auto *snh1_1227 = buffer.data(snh1 + 1227);
    const auto *snh1_1228 = buffer.data(snh1 + 1228);
    const auto *snh1_1230 = buffer.data(snh1 + 1230);
    const auto *snh1_1232 = buffer.data(snh1 + 1232);
    const auto *snh1_1233 = buffer.data(snh1 + 1233);
    const auto *snh1_1235 = buffer.data(snh1 + 1235);
    const auto *snh1_1236 = buffer.data(snh1 + 1236);
    const auto *snh1_1237 = buffer.data(snh1 + 1237);
    const auto *snh1_1238 = buffer.data(snh1 + 1238);
    const auto *snh1_1239 = buffer.data(snh1 + 1239);
    const auto *snh1_1242 = buffer.data(snh1 + 1242);
    const auto *snh1_1244 = buffer.data(snh1 + 1244);
    const auto *snh1_1245 = buffer.data(snh1 + 1245);
    const auto *snh1_1248 = buffer.data(snh1 + 1248);
    const auto *snh1_1249 = buffer.data(snh1 + 1249);
    const auto *snh1_1251 = buffer.data(snh1 + 1251);
    const auto *snh1_1253 = buffer.data(snh1 + 1253);

    const auto *sni_1570 = buffer.data(sni + 1570);
    const auto *sni_1571 = buffer.data(sni + 1571);
    const auto *sni_1573 = buffer.data(sni + 1573);
    const auto *sni_1574 = buffer.data(sni + 1574);
    const auto *sni_1577 = buffer.data(sni + 1577);
    const auto *sni_1578 = buffer.data(sni + 1578);
    const auto *sni_1580 = buffer.data(sni + 1580);
    const auto *sni_1582 = buffer.data(sni + 1582);
    const auto *sni_1585 = buffer.data(sni + 1585);
    const auto *sni_1586 = buffer.data(sni + 1586);
    const auto *sni_1588 = buffer.data(sni + 1588);
    const auto *sni_1589 = buffer.data(sni + 1589);
    const auto *sni_1590 = buffer.data(sni + 1590);
    const auto *sni_1591 = buffer.data(sni + 1591);
    const auto *sni_1592 = buffer.data(sni + 1592);
    const auto *sni_1593 = buffer.data(sni + 1593);
    const auto *sni_1594 = buffer.data(sni + 1594);
    const auto *sni_1595 = buffer.data(sni + 1595);
    const auto *sni_1596 = buffer.data(sni + 1596);
    const auto *sni_1598 = buffer.data(sni + 1598);
    const auto *sni_1599 = buffer.data(sni + 1599);
    const auto *sni_1601 = buffer.data(sni + 1601);
    const auto *sni_1602 = buffer.data(sni + 1602);
    const auto *sni_1605 = buffer.data(sni + 1605);
    const auto *sni_1606 = buffer.data(sni + 1606);
    const auto *sni_1608 = buffer.data(sni + 1608);
    const auto *sni_1610 = buffer.data(sni + 1610);
    const auto *sni_1611 = buffer.data(sni + 1611);
    const auto *sni_1613 = buffer.data(sni + 1613);
    const auto *sni_1614 = buffer.data(sni + 1614);
    const auto *sni_1616 = buffer.data(sni + 1616);
    const auto *sni_1617 = buffer.data(sni + 1617);
    const auto *sni_1618 = buffer.data(sni + 1618);
    const auto *sni_1619 = buffer.data(sni + 1619);
    const auto *sni_1620 = buffer.data(sni + 1620);
    const auto *sni_1621 = buffer.data(sni + 1621);
    const auto *sni_1622 = buffer.data(sni + 1622);
    const auto *sni_1623 = buffer.data(sni + 1623);
    const auto *sni_1624 = buffer.data(sni + 1624);
    const auto *sni_1626 = buffer.data(sni + 1626);
    const auto *sni_1627 = buffer.data(sni + 1627);
    const auto *sni_1629 = buffer.data(sni + 1629);
    const auto *sni_1630 = buffer.data(sni + 1630);
    const auto *sni_1633 = buffer.data(sni + 1633);
    const auto *sni_1634 = buffer.data(sni + 1634);
    const auto *sni_1636 = buffer.data(sni + 1636);
    const auto *sni_1638 = buffer.data(sni + 1638);
    const auto *sni_1639 = buffer.data(sni + 1639);
    const auto *sni_1641 = buffer.data(sni + 1641);
    const auto *sni_1642 = buffer.data(sni + 1642);
    const auto *sni_1644 = buffer.data(sni + 1644);
    const auto *sni_1645 = buffer.data(sni + 1645);
    const auto *sni_1646 = buffer.data(sni + 1646);
    const auto *sni_1647 = buffer.data(sni + 1647);
    const auto *sni_1648 = buffer.data(sni + 1648);
    const auto *sni_1649 = buffer.data(sni + 1649);
    const auto *sni_1650 = buffer.data(sni + 1650);
    const auto *sni_1651 = buffer.data(sni + 1651);
    const auto *sni_1652 = buffer.data(sni + 1652);
    const auto *sni_1654 = buffer.data(sni + 1654);
    const auto *sni_1655 = buffer.data(sni + 1655);
    const auto *sni_1657 = buffer.data(sni + 1657);
    const auto *sni_1658 = buffer.data(sni + 1658);
    const auto *sni_1661 = buffer.data(sni + 1661);
    const auto *sni_1662 = buffer.data(sni + 1662);
    const auto *sni_1664 = buffer.data(sni + 1664);
    const auto *sni_1666 = buffer.data(sni + 1666);

#pragma omp simd aligned(t_2020, t_2021, t_2022, pb_z, pc_x, pc_y, pc_z, smk0_1626, smi_1290, \
                         smk1_1626, snh0_1181, snh1_1181, sni_1570, \
                         sni_1573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2020[k] = f_18 * smi_1290[k]
                    + f_3 * pc_y[k] * sni_1570[k];

        t_2021[k] = f_4 * snh0_1181[k]
                    - f_5 * snh1_1181[k]
                    + f_3 * pc_x[k] * sni_1573[k];

        t_2022[k] = pb_z[k] * smk0_1626[k]
                    - f_12 * pc_z[k] * smk1_1626[k];
    }

#pragma omp simd aligned(t_2023, t_2024, t_2025, pc_x, pc_y, pc_z, smi_1263, smi_1293, \
                         snh0_1185, snh1_1185, sni_1571, sni_1573, \
                         sni_1577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2023[k] = f_13 * smi_1263[k]
                    + f_3 * pc_z[k] * sni_1571[k];

        t_2024[k] = f_18 * smi_1293[k]
                    + f_3 * pc_y[k] * sni_1573[k];

        t_2025[k] = f_6 * snh0_1185[k]
                    - f_7 * snh1_1185[k]
                    + f_3 * pc_x[k] * sni_1577[k];
    }

#pragma omp simd aligned(t_2026, t_2027, t_2028, pb_z, pc_x, pc_z, smk0_1630, smi_1266, \
                         smk1_1630, snh0_1188, snh1_1188, sni_1574, \
                         sni_1580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2026[k] = pb_z[k] * smk0_1630[k]
                    - f_12 * pc_z[k] * smk1_1630[k];

        t_2027[k] = f_13 * smi_1266[k]
                    + f_3 * pc_z[k] * sni_1574[k];

        t_2028[k] = f_8 * snh0_1188[k]
                    - f_9 * snh1_1188[k]
                    + f_3 * pc_x[k] * sni_1580[k];
    }

#pragma omp simd aligned(t_2029, t_2030, t_2031, pb_z, pc_x, pc_y, pc_z, smk0_1635, smi_1297, \
                         smk1_1635, snh0_1190, snh1_1190, sni_1577, \
                         sni_1582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2029[k] = f_18 * smi_1297[k]
                    + f_3 * pc_y[k] * sni_1577[k];

        t_2030[k] = f_8 * snh0_1190[k]
                    - f_9 * snh1_1190[k]
                    + f_3 * pc_x[k] * sni_1582[k];

        t_2031[k] = pb_z[k] * smk0_1635[k]
                    - f_12 * pc_z[k] * smk1_1635[k];
    }

#pragma omp simd aligned(t_2032, t_2033, t_2034, pc_x, pc_z, smi_1270, snh0_1193, snh0_1194, \
                         snh1_1193, snh1_1194, sni_1578, sni_1585, \
                         sni_1586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2032[k] = f_13 * smi_1270[k]
                    + f_3 * pc_z[k] * sni_1578[k];

        t_2033[k] = f_10 * snh0_1193[k]
                    - f_11 * snh1_1193[k]
                    + f_3 * pc_x[k] * sni_1585[k];

        t_2034[k] = f_10 * snh0_1194[k]
                    - f_11 * snh1_1194[k]
                    + f_3 * pc_x[k] * sni_1586[k];
    }

#pragma omp simd aligned(t_2035, t_2036, t_2037, t_2038, t_2039, pc_x, pc_y, smi_1302, \
                         snh0_1196, snh1_1196, sni_1582, sni_1588, sni_1589, sni_1590, \
                         sni_1591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2035[k] = f_18 * smi_1302[k]
                    + f_3 * pc_y[k] * sni_1582[k];

        t_2036[k] = f_10 * snh0_1196[k]
                    - f_11 * snh1_1196[k]
                    + f_3 * pc_x[k] * sni_1588[k];

        t_2037[k] = f_3 * pc_x[k] * sni_1589[k];

        t_2038[k] = f_3 * pc_x[k] * sni_1590[k];

        t_2039[k] = f_3 * pc_x[k] * sni_1591[k];
    }

#pragma omp simd aligned(t_2040, t_2041, t_2042, t_2043, t_2044, pb_z, pc_x, pc_z, smk0_1648, \
                         smk1_1648, sni_1592, sni_1593, sni_1594, \
                         sni_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2040[k] = f_3 * pc_x[k] * sni_1592[k];

        t_2041[k] = f_3 * pc_x[k] * sni_1593[k];

        t_2042[k] = f_3 * pc_x[k] * sni_1594[k];

        t_2043[k] = f_3 * pc_x[k] * sni_1595[k];

        t_2044[k] = pb_z[k] * smk0_1648[k]
                    - f_12 * pc_z[k] * smk1_1648[k];
    }

#pragma omp simd aligned(t_2045, t_2046, t_2047, pb_z, pc_z, smk0_1650, smk0_1651, smi_1281, \
                         smi_1282, smi_1283, smk1_1650, smk1_1651, \
                         sni_1589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2045[k] = f_13 * smi_1281[k]
                    + f_3 * pc_z[k] * sni_1589[k];

        t_2046[k] = pb_z[k] * smk0_1650[k]
                    + f_14 * smi_1282[k]
                    - f_12 * pc_z[k] * smk1_1650[k];

        t_2047[k] = pb_z[k] * smk0_1651[k]
                    + f_15 * smi_1283[k]
                    - f_12 * pc_z[k] * smk1_1651[k];
    }

#pragma omp simd aligned(t_2048, t_2049, t_2050, pb_z, pc_y, pc_z, smk0_1652, smk0_1653, \
                         smi_1284, smi_1285, smi_1315, smk1_1652, smk1_1653, \
                         sni_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2048[k] = pb_z[k] * smk0_1652[k]
                    + f_16 * smi_1284[k]
                    - f_12 * pc_z[k] * smk1_1652[k];

        t_2049[k] = pb_z[k] * smk0_1653[k]
                    + f_17 * smi_1285[k]
                    - f_12 * pc_z[k] * smk1_1653[k];

        t_2050[k] = f_18 * smi_1315[k]
                    + f_3 * pc_y[k] * sni_1595[k];
    }

#pragma omp simd aligned(t_2051, t_2052, t_2053, t_2054, pc_x, pc_y, pc_z, smi_1287, smi_1288, \
                         smi_1316, snh0_1196, snh0_1197, snh1_1196, snh1_1197, sni_1595, \
                         sni_1596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2051[k] = f_13 * smi_1287[k]
                    + f_1 * snh0_1196[k]
                    - f_2 * snh1_1196[k]
                    + f_3 * pc_z[k] * sni_1595[k];

        t_2052[k] = f_1 * snh0_1197[k]
                    - f_2 * snh1_1197[k]
                    + f_3 * pc_x[k] * sni_1596[k];

        t_2053[k] = f_19 * smi_1316[k]
                    + f_3 * pc_y[k] * sni_1596[k];

        t_2054[k] = f_14 * smi_1288[k]
                    + f_3 * pc_z[k] * sni_1596[k];
    }

#pragma omp simd aligned(t_2055, t_2056, t_2057, pc_x, pc_y, smi_1318, snh0_1200, snh0_1202, \
                         snh1_1200, snh1_1202, sni_1598, sni_1599, \
                         sni_1601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2055[k] = f_4 * snh0_1200[k]
                    - f_5 * snh1_1200[k]
                    + f_3 * pc_x[k] * sni_1599[k];

        t_2056[k] = f_19 * smi_1318[k]
                    + f_3 * pc_y[k] * sni_1598[k];

        t_2057[k] = f_4 * snh0_1202[k]
                    - f_5 * snh1_1202[k]
                    + f_3 * pc_x[k] * sni_1601[k];
    }

#pragma omp simd aligned(t_2058, t_2059, t_2060, pc_x, pc_y, pc_z, smi_1291, smi_1321, \
                         snh0_1203, snh1_1203, sni_1599, sni_1601, \
                         sni_1602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2058[k] = f_6 * snh0_1203[k]
                    - f_7 * snh1_1203[k]
                    + f_3 * pc_x[k] * sni_1602[k];

        t_2059[k] = f_14 * smi_1291[k]
                    + f_3 * pc_z[k] * sni_1599[k];

        t_2060[k] = f_19 * smi_1321[k]
                    + f_3 * pc_y[k] * sni_1601[k];
    }

#pragma omp simd aligned(t_2061, t_2062, t_2063, pc_x, pc_z, smi_1294, snh0_1206, snh0_1207, \
                         snh1_1206, snh1_1207, sni_1602, sni_1605, \
                         sni_1606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2061[k] = f_6 * snh0_1206[k]
                    - f_7 * snh1_1206[k]
                    + f_3 * pc_x[k] * sni_1605[k];

        t_2062[k] = f_8 * snh0_1207[k]
                    - f_9 * snh1_1207[k]
                    + f_3 * pc_x[k] * sni_1606[k];

        t_2063[k] = f_14 * smi_1294[k]
                    + f_3 * pc_z[k] * sni_1602[k];
    }

#pragma omp simd aligned(t_2064, t_2065, t_2066, pc_x, pc_y, smi_1325, snh0_1209, snh0_1211, \
                         snh1_1209, snh1_1211, sni_1605, sni_1608, \
                         sni_1610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2064[k] = f_8 * snh0_1209[k]
                    - f_9 * snh1_1209[k]
                    + f_3 * pc_x[k] * sni_1608[k];

        t_2065[k] = f_19 * smi_1325[k]
                    + f_3 * pc_y[k] * sni_1605[k];

        t_2066[k] = f_8 * snh0_1211[k]
                    - f_9 * snh1_1211[k]
                    + f_3 * pc_x[k] * sni_1610[k];
    }

#pragma omp simd aligned(t_2067, t_2068, t_2069, pc_x, pc_z, smi_1298, snh0_1212, snh0_1214, \
                         snh1_1212, snh1_1214, sni_1606, sni_1611, \
                         sni_1613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2067[k] = f_10 * snh0_1212[k]
                    - f_11 * snh1_1212[k]
                    + f_3 * pc_x[k] * sni_1611[k];

        t_2068[k] = f_14 * smi_1298[k]
                    + f_3 * pc_z[k] * sni_1606[k];

        t_2069[k] = f_10 * snh0_1214[k]
                    - f_11 * snh1_1214[k]
                    + f_3 * pc_x[k] * sni_1613[k];
    }

#pragma omp simd aligned(t_2070, t_2071, t_2072, t_2073, pc_x, pc_y, smi_1330, snh0_1215, \
                         snh0_1217, snh1_1215, snh1_1217, sni_1610, sni_1614, sni_1616, \
                         sni_1617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2070[k] = f_10 * snh0_1215[k]
                    - f_11 * snh1_1215[k]
                    + f_3 * pc_x[k] * sni_1614[k];

        t_2071[k] = f_19 * smi_1330[k]
                    + f_3 * pc_y[k] * sni_1610[k];

        t_2072[k] = f_10 * snh0_1217[k]
                    - f_11 * snh1_1217[k]
                    + f_3 * pc_x[k] * sni_1616[k];

        t_2073[k] = f_3 * pc_x[k] * sni_1617[k];
    }

#pragma omp simd aligned(t_2074, t_2075, t_2076, t_2077, t_2078, t_2079, pc_x, sni_1618, \
                         sni_1619, sni_1620, sni_1621, sni_1622, \
                         sni_1623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2074[k] = f_3 * pc_x[k] * sni_1618[k];

        t_2075[k] = f_3 * pc_x[k] * sni_1619[k];

        t_2076[k] = f_3 * pc_x[k] * sni_1620[k];

        t_2077[k] = f_3 * pc_x[k] * sni_1621[k];

        t_2078[k] = f_3 * pc_x[k] * sni_1622[k];

        t_2079[k] = f_3 * pc_x[k] * sni_1623[k];
    }

#pragma omp simd aligned(t_2080, t_2081, t_2082, pc_y, pc_z, smi_1309, smi_1337, smi_1339, \
                         snh0_1212, snh0_1214, snh1_1212, snh1_1214, sni_1617, \
                         sni_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2080[k] = f_19 * smi_1337[k]
                    + f_1 * snh0_1212[k]
                    - f_2 * snh1_1212[k]
                    + f_3 * pc_y[k] * sni_1617[k];

        t_2081[k] = f_14 * smi_1309[k]
                    + f_3 * pc_z[k] * sni_1617[k];

        t_2082[k] = f_19 * smi_1339[k]
                    + f_4 * snh0_1214[k]
                    - f_5 * snh1_1214[k]
                    + f_3 * pc_y[k] * sni_1619[k];
    }

#pragma omp simd aligned(t_2083, t_2084, t_2085, pc_y, smi_1340, smi_1341, smi_1342, \
                         snh0_1215, snh0_1216, snh0_1217, snh1_1215, snh1_1216, snh1_1217, \
                         sni_1620, sni_1621, sni_1622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2083[k] = f_19 * smi_1340[k]
                    + f_6 * snh0_1215[k]
                    - f_7 * snh1_1215[k]
                    + f_3 * pc_y[k] * sni_1620[k];

        t_2084[k] = f_19 * smi_1341[k]
                    + f_8 * snh0_1216[k]
                    - f_9 * snh1_1216[k]
                    + f_3 * pc_y[k] * sni_1621[k];

        t_2085[k] = f_19 * smi_1342[k]
                    + f_10 * snh0_1217[k]
                    - f_11 * snh1_1217[k]
                    + f_3 * pc_y[k] * sni_1622[k];
    }

#pragma omp simd aligned(t_2086, t_2087, t_2088, t_2089, pc_x, pc_y, pc_z, smi_1315, smi_1343, \
                         smi_1344, snh0_1217, snh0_1218, snh1_1217, snh1_1218, sni_1623, \
                         sni_1624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2086[k] = f_19 * smi_1343[k]
                    + f_3 * pc_y[k] * sni_1623[k];

        t_2087[k] = f_14 * smi_1315[k]
                    + f_1 * snh0_1217[k]
                    - f_2 * snh1_1217[k]
                    + f_3 * pc_z[k] * sni_1623[k];

        t_2088[k] = f_1 * snh0_1218[k]
                    - f_2 * snh1_1218[k]
                    + f_3 * pc_x[k] * sni_1624[k];

        t_2089[k] = f_20 * smi_1344[k]
                    + f_3 * pc_y[k] * sni_1624[k];
    }

#pragma omp simd aligned(t_2090, t_2091, t_2092, pc_x, pc_y, pc_z, smi_1316, smi_1346, \
                         snh0_1221, snh1_1221, sni_1624, sni_1626, \
                         sni_1627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2090[k] = f_15 * smi_1316[k]
                    + f_3 * pc_z[k] * sni_1624[k];

        t_2091[k] = f_4 * snh0_1221[k]
                    - f_5 * snh1_1221[k]
                    + f_3 * pc_x[k] * sni_1627[k];

        t_2092[k] = f_20 * smi_1346[k]
                    + f_3 * pc_y[k] * sni_1626[k];
    }

#pragma omp simd aligned(t_2093, t_2094, t_2095, t_2096, pc_x, pc_y, pc_z, smi_1319, smi_1349, \
                         snh0_1223, snh0_1224, snh1_1223, snh1_1224, sni_1627, sni_1629, \
                         sni_1630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2093[k] = f_4 * snh0_1223[k]
                    - f_5 * snh1_1223[k]
                    + f_3 * pc_x[k] * sni_1629[k];

        t_2094[k] = f_6 * snh0_1224[k]
                    - f_7 * snh1_1224[k]
                    + f_3 * pc_x[k] * sni_1630[k];

        t_2095[k] = f_15 * smi_1319[k]
                    + f_3 * pc_z[k] * sni_1627[k];

        t_2096[k] = f_20 * smi_1349[k]
                    + f_3 * pc_y[k] * sni_1629[k];
    }

#pragma omp simd aligned(t_2097, t_2098, t_2099, pc_x, pc_z, smi_1322, snh0_1227, snh0_1228, \
                         snh1_1227, snh1_1228, sni_1630, sni_1633, \
                         sni_1634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2097[k] = f_6 * snh0_1227[k]
                    - f_7 * snh1_1227[k]
                    + f_3 * pc_x[k] * sni_1633[k];

        t_2098[k] = f_8 * snh0_1228[k]
                    - f_9 * snh1_1228[k]
                    + f_3 * pc_x[k] * sni_1634[k];

        t_2099[k] = f_15 * smi_1322[k]
                    + f_3 * pc_z[k] * sni_1630[k];
    }

#pragma omp simd aligned(t_2100, t_2101, t_2102, pc_x, pc_y, smi_1353, snh0_1230, snh0_1232, \
                         snh1_1230, snh1_1232, sni_1633, sni_1636, \
                         sni_1638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2100[k] = f_8 * snh0_1230[k]
                    - f_9 * snh1_1230[k]
                    + f_3 * pc_x[k] * sni_1636[k];

        t_2101[k] = f_20 * smi_1353[k]
                    + f_3 * pc_y[k] * sni_1633[k];

        t_2102[k] = f_8 * snh0_1232[k]
                    - f_9 * snh1_1232[k]
                    + f_3 * pc_x[k] * sni_1638[k];
    }

#pragma omp simd aligned(t_2103, t_2104, t_2105, pc_x, pc_z, smi_1326, snh0_1233, snh0_1235, \
                         snh1_1233, snh1_1235, sni_1634, sni_1639, \
                         sni_1641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2103[k] = f_10 * snh0_1233[k]
                    - f_11 * snh1_1233[k]
                    + f_3 * pc_x[k] * sni_1639[k];

        t_2104[k] = f_15 * smi_1326[k]
                    + f_3 * pc_z[k] * sni_1634[k];

        t_2105[k] = f_10 * snh0_1235[k]
                    - f_11 * snh1_1235[k]
                    + f_3 * pc_x[k] * sni_1641[k];
    }

#pragma omp simd aligned(t_2106, t_2107, t_2108, t_2109, pc_x, pc_y, smi_1358, snh0_1236, \
                         snh0_1238, snh1_1236, snh1_1238, sni_1638, sni_1642, sni_1644, \
                         sni_1645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2106[k] = f_10 * snh0_1236[k]
                    - f_11 * snh1_1236[k]
                    + f_3 * pc_x[k] * sni_1642[k];

        t_2107[k] = f_20 * smi_1358[k]
                    + f_3 * pc_y[k] * sni_1638[k];

        t_2108[k] = f_10 * snh0_1238[k]
                    - f_11 * snh1_1238[k]
                    + f_3 * pc_x[k] * sni_1644[k];

        t_2109[k] = f_3 * pc_x[k] * sni_1645[k];
    }

#pragma omp simd aligned(t_2110, t_2111, t_2112, t_2113, t_2114, t_2115, pc_x, sni_1646, \
                         sni_1647, sni_1648, sni_1649, sni_1650, \
                         sni_1651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2110[k] = f_3 * pc_x[k] * sni_1646[k];

        t_2111[k] = f_3 * pc_x[k] * sni_1647[k];

        t_2112[k] = f_3 * pc_x[k] * sni_1648[k];

        t_2113[k] = f_3 * pc_x[k] * sni_1649[k];

        t_2114[k] = f_3 * pc_x[k] * sni_1650[k];

        t_2115[k] = f_3 * pc_x[k] * sni_1651[k];
    }

#pragma omp simd aligned(t_2116, t_2117, t_2118, pc_y, pc_z, smi_1337, smi_1365, smi_1367, \
                         snh0_1233, snh0_1235, snh1_1233, snh1_1235, sni_1645, \
                         sni_1647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2116[k] = f_20 * smi_1365[k]
                    + f_1 * snh0_1233[k]
                    - f_2 * snh1_1233[k]
                    + f_3 * pc_y[k] * sni_1645[k];

        t_2117[k] = f_15 * smi_1337[k]
                    + f_3 * pc_z[k] * sni_1645[k];

        t_2118[k] = f_20 * smi_1367[k]
                    + f_4 * snh0_1235[k]
                    - f_5 * snh1_1235[k]
                    + f_3 * pc_y[k] * sni_1647[k];
    }

#pragma omp simd aligned(t_2119, t_2120, t_2121, pc_y, smi_1368, smi_1369, smi_1370, \
                         snh0_1236, snh0_1237, snh0_1238, snh1_1236, snh1_1237, snh1_1238, \
                         sni_1648, sni_1649, sni_1650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2119[k] = f_20 * smi_1368[k]
                    + f_6 * snh0_1236[k]
                    - f_7 * snh1_1236[k]
                    + f_3 * pc_y[k] * sni_1648[k];

        t_2120[k] = f_20 * smi_1369[k]
                    + f_8 * snh0_1237[k]
                    - f_9 * snh1_1237[k]
                    + f_3 * pc_y[k] * sni_1649[k];

        t_2121[k] = f_20 * smi_1370[k]
                    + f_10 * snh0_1238[k]
                    - f_11 * snh1_1238[k]
                    + f_3 * pc_y[k] * sni_1650[k];
    }

#pragma omp simd aligned(t_2122, t_2123, t_2124, t_2125, pc_x, pc_y, pc_z, smi_1343, smi_1371, \
                         smi_1372, snh0_1238, snh0_1239, snh1_1238, snh1_1239, sni_1651, \
                         sni_1652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2122[k] = f_20 * smi_1371[k]
                    + f_3 * pc_y[k] * sni_1651[k];

        t_2123[k] = f_15 * smi_1343[k]
                    + f_1 * snh0_1238[k]
                    - f_2 * snh1_1238[k]
                    + f_3 * pc_z[k] * sni_1651[k];

        t_2124[k] = f_1 * snh0_1239[k]
                    - f_2 * snh1_1239[k]
                    + f_3 * pc_x[k] * sni_1652[k];

        t_2125[k] = f_21 * smi_1372[k]
                    + f_3 * pc_y[k] * sni_1652[k];
    }

#pragma omp simd aligned(t_2126, t_2127, t_2128, pc_x, pc_y, pc_z, smi_1344, smi_1374, \
                         snh0_1242, snh1_1242, sni_1652, sni_1654, \
                         sni_1655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2126[k] = f_16 * smi_1344[k]
                    + f_3 * pc_z[k] * sni_1652[k];

        t_2127[k] = f_4 * snh0_1242[k]
                    - f_5 * snh1_1242[k]
                    + f_3 * pc_x[k] * sni_1655[k];

        t_2128[k] = f_21 * smi_1374[k]
                    + f_3 * pc_y[k] * sni_1654[k];
    }

#pragma omp simd aligned(t_2129, t_2130, t_2131, t_2132, pc_x, pc_y, pc_z, smi_1347, smi_1377, \
                         snh0_1244, snh0_1245, snh1_1244, snh1_1245, sni_1655, sni_1657, \
                         sni_1658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2129[k] = f_4 * snh0_1244[k]
                    - f_5 * snh1_1244[k]
                    + f_3 * pc_x[k] * sni_1657[k];

        t_2130[k] = f_6 * snh0_1245[k]
                    - f_7 * snh1_1245[k]
                    + f_3 * pc_x[k] * sni_1658[k];

        t_2131[k] = f_16 * smi_1347[k]
                    + f_3 * pc_z[k] * sni_1655[k];

        t_2132[k] = f_21 * smi_1377[k]
                    + f_3 * pc_y[k] * sni_1657[k];
    }

#pragma omp simd aligned(t_2133, t_2134, t_2135, pc_x, pc_z, smi_1350, snh0_1248, snh0_1249, \
                         snh1_1248, snh1_1249, sni_1658, sni_1661, \
                         sni_1662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2133[k] = f_6 * snh0_1248[k]
                    - f_7 * snh1_1248[k]
                    + f_3 * pc_x[k] * sni_1661[k];

        t_2134[k] = f_8 * snh0_1249[k]
                    - f_9 * snh1_1249[k]
                    + f_3 * pc_x[k] * sni_1662[k];

        t_2135[k] = f_16 * smi_1350[k]
                    + f_3 * pc_z[k] * sni_1658[k];
    }

#pragma omp simd aligned(t_2136, t_2137, t_2138, pc_x, pc_y, smi_1381, snh0_1251, snh0_1253, \
                         snh1_1251, snh1_1253, sni_1661, sni_1664, \
                         sni_1666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2136[k] = f_8 * snh0_1251[k]
                    - f_9 * snh1_1251[k]
                    + f_3 * pc_x[k] * sni_1664[k];

        t_2137[k] = f_21 * smi_1381[k]
                    + f_3 * pc_y[k] * sni_1661[k];

        t_2138[k] = f_8 * snh0_1253[k]
                    - f_9 * snh1_1253[k]
                    + f_3 * pc_x[k] * sni_1666[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece19(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t smi, const size_t snh0,
                                                           const size_t snh1, const size_t sni,
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
    const auto f_15 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_20 = 3.5 / q;
    const auto f_21 = 3.0 / q;

    auto *t_2139 = buffer.data(target + 2139);
    auto *t_2140 = buffer.data(target + 2140);
    auto *t_2141 = buffer.data(target + 2141);
    auto *t_2142 = buffer.data(target + 2142);
    auto *t_2143 = buffer.data(target + 2143);
    auto *t_2144 = buffer.data(target + 2144);
    auto *t_2145 = buffer.data(target + 2145);
    auto *t_2146 = buffer.data(target + 2146);
    auto *t_2147 = buffer.data(target + 2147);
    auto *t_2148 = buffer.data(target + 2148);
    auto *t_2149 = buffer.data(target + 2149);
    auto *t_2150 = buffer.data(target + 2150);
    auto *t_2151 = buffer.data(target + 2151);
    auto *t_2152 = buffer.data(target + 2152);
    auto *t_2153 = buffer.data(target + 2153);
    auto *t_2154 = buffer.data(target + 2154);
    auto *t_2155 = buffer.data(target + 2155);
    auto *t_2156 = buffer.data(target + 2156);
    auto *t_2157 = buffer.data(target + 2157);
    auto *t_2158 = buffer.data(target + 2158);
    auto *t_2159 = buffer.data(target + 2159);
    auto *t_2160 = buffer.data(target + 2160);
    auto *t_2161 = buffer.data(target + 2161);
    auto *t_2162 = buffer.data(target + 2162);
    auto *t_2163 = buffer.data(target + 2163);
    auto *t_2164 = buffer.data(target + 2164);
    auto *t_2165 = buffer.data(target + 2165);
    auto *t_2166 = buffer.data(target + 2166);
    auto *t_2167 = buffer.data(target + 2167);
    auto *t_2168 = buffer.data(target + 2168);
    auto *t_2169 = buffer.data(target + 2169);
    auto *t_2170 = buffer.data(target + 2170);
    auto *t_2171 = buffer.data(target + 2171);
    auto *t_2172 = buffer.data(target + 2172);
    auto *t_2173 = buffer.data(target + 2173);
    auto *t_2174 = buffer.data(target + 2174);
    auto *t_2175 = buffer.data(target + 2175);
    auto *t_2176 = buffer.data(target + 2176);
    auto *t_2177 = buffer.data(target + 2177);
    auto *t_2178 = buffer.data(target + 2178);
    auto *t_2179 = buffer.data(target + 2179);
    auto *t_2180 = buffer.data(target + 2180);
    auto *t_2181 = buffer.data(target + 2181);
    auto *t_2182 = buffer.data(target + 2182);
    auto *t_2183 = buffer.data(target + 2183);
    auto *t_2184 = buffer.data(target + 2184);
    auto *t_2185 = buffer.data(target + 2185);
    auto *t_2186 = buffer.data(target + 2186);
    auto *t_2187 = buffer.data(target + 2187);
    auto *t_2188 = buffer.data(target + 2188);
    auto *t_2189 = buffer.data(target + 2189);
    auto *t_2190 = buffer.data(target + 2190);
    auto *t_2191 = buffer.data(target + 2191);
    auto *t_2192 = buffer.data(target + 2192);
    auto *t_2193 = buffer.data(target + 2193);
    auto *t_2194 = buffer.data(target + 2194);
    auto *t_2195 = buffer.data(target + 2195);
    auto *t_2196 = buffer.data(target + 2196);
    auto *t_2197 = buffer.data(target + 2197);
    auto *t_2198 = buffer.data(target + 2198);
    auto *t_2199 = buffer.data(target + 2199);
    auto *t_2200 = buffer.data(target + 2200);
    auto *t_2201 = buffer.data(target + 2201);
    auto *t_2202 = buffer.data(target + 2202);
    auto *t_2203 = buffer.data(target + 2203);
    auto *t_2204 = buffer.data(target + 2204);
    auto *t_2205 = buffer.data(target + 2205);
    auto *t_2206 = buffer.data(target + 2206);
    auto *t_2207 = buffer.data(target + 2207);
    auto *t_2208 = buffer.data(target + 2208);
    auto *t_2209 = buffer.data(target + 2209);
    auto *t_2210 = buffer.data(target + 2210);
    auto *t_2211 = buffer.data(target + 2211);
    auto *t_2212 = buffer.data(target + 2212);
    auto *t_2213 = buffer.data(target + 2213);
    auto *t_2214 = buffer.data(target + 2214);
    auto *t_2215 = buffer.data(target + 2215);
    auto *t_2216 = buffer.data(target + 2216);
    auto *t_2217 = buffer.data(target + 2217);
    auto *t_2218 = buffer.data(target + 2218);
    auto *t_2219 = buffer.data(target + 2219);
    auto *t_2220 = buffer.data(target + 2220);
    auto *t_2221 = buffer.data(target + 2221);
    auto *t_2222 = buffer.data(target + 2222);
    auto *t_2223 = buffer.data(target + 2223);
    auto *t_2224 = buffer.data(target + 2224);
    auto *t_2225 = buffer.data(target + 2225);
    auto *t_2226 = buffer.data(target + 2226);
    auto *t_2227 = buffer.data(target + 2227);
    auto *t_2228 = buffer.data(target + 2228);
    auto *t_2229 = buffer.data(target + 2229);
    auto *t_2230 = buffer.data(target + 2230);
    auto *t_2231 = buffer.data(target + 2231);
    auto *t_2232 = buffer.data(target + 2232);
    auto *t_2233 = buffer.data(target + 2233);
    auto *t_2234 = buffer.data(target + 2234);
    auto *t_2235 = buffer.data(target + 2235);
    auto *t_2236 = buffer.data(target + 2236);
    auto *t_2237 = buffer.data(target + 2237);
    auto *t_2238 = buffer.data(target + 2238);
    auto *t_2239 = buffer.data(target + 2239);
    auto *t_2240 = buffer.data(target + 2240);
    auto *t_2241 = buffer.data(target + 2241);
    auto *t_2242 = buffer.data(target + 2242);
    auto *t_2243 = buffer.data(target + 2243);
    auto *t_2244 = buffer.data(target + 2244);
    auto *t_2245 = buffer.data(target + 2245);
    auto *t_2246 = buffer.data(target + 2246);
    auto *t_2247 = buffer.data(target + 2247);
    auto *t_2248 = buffer.data(target + 2248);
    auto *t_2249 = buffer.data(target + 2249);
    auto *t_2250 = buffer.data(target + 2250);
    auto *t_2251 = buffer.data(target + 2251);
    auto *t_2252 = buffer.data(target + 2252);
    auto *t_2253 = buffer.data(target + 2253);
    auto *t_2254 = buffer.data(target + 2254);
    auto *t_2255 = buffer.data(target + 2255);
    auto *t_2256 = buffer.data(target + 2256);
    auto *t_2257 = buffer.data(target + 2257);
    auto *t_2258 = buffer.data(target + 2258);
    auto *t_2259 = buffer.data(target + 2259);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smi_1354 = buffer.data(smi + 1354);
    const auto *smi_1365 = buffer.data(smi + 1365);
    const auto *smi_1371 = buffer.data(smi + 1371);
    const auto *smi_1372 = buffer.data(smi + 1372);
    const auto *smi_1375 = buffer.data(smi + 1375);
    const auto *smi_1378 = buffer.data(smi + 1378);
    const auto *smi_1382 = buffer.data(smi + 1382);
    const auto *smi_1386 = buffer.data(smi + 1386);
    const auto *smi_1393 = buffer.data(smi + 1393);
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
    const auto *smi_1414 = buffer.data(smi + 1414);
    const auto *smi_1421 = buffer.data(smi + 1421);
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
    const auto *smi_1442 = buffer.data(smi + 1442);
    const auto *smi_1449 = buffer.data(smi + 1449);
    const auto *smi_1451 = buffer.data(smi + 1451);
    const auto *smi_1452 = buffer.data(smi + 1452);
    const auto *smi_1453 = buffer.data(smi + 1453);
    const auto *smi_1454 = buffer.data(smi + 1454);
    const auto *smi_1455 = buffer.data(smi + 1455);
    const auto *smi_1456 = buffer.data(smi + 1456);
    const auto *smi_1458 = buffer.data(smi + 1458);
    const auto *smi_1461 = buffer.data(smi + 1461);
    const auto *smi_1465 = buffer.data(smi + 1465);
    const auto *smi_1470 = buffer.data(smi + 1470);

    const auto *snh0_1254 = buffer.data(snh0 + 1254);
    const auto *snh0_1256 = buffer.data(snh0 + 1256);
    const auto *snh0_1257 = buffer.data(snh0 + 1257);
    const auto *snh0_1258 = buffer.data(snh0 + 1258);
    const auto *snh0_1259 = buffer.data(snh0 + 1259);
    const auto *snh0_1260 = buffer.data(snh0 + 1260);
    const auto *snh0_1263 = buffer.data(snh0 + 1263);
    const auto *snh0_1265 = buffer.data(snh0 + 1265);
    const auto *snh0_1266 = buffer.data(snh0 + 1266);
    const auto *snh0_1269 = buffer.data(snh0 + 1269);
    const auto *snh0_1270 = buffer.data(snh0 + 1270);
    const auto *snh0_1272 = buffer.data(snh0 + 1272);
    const auto *snh0_1274 = buffer.data(snh0 + 1274);
    const auto *snh0_1275 = buffer.data(snh0 + 1275);
    const auto *snh0_1277 = buffer.data(snh0 + 1277);
    const auto *snh0_1278 = buffer.data(snh0 + 1278);
    const auto *snh0_1279 = buffer.data(snh0 + 1279);
    const auto *snh0_1280 = buffer.data(snh0 + 1280);
    const auto *snh0_1281 = buffer.data(snh0 + 1281);
    const auto *snh0_1284 = buffer.data(snh0 + 1284);
    const auto *snh0_1286 = buffer.data(snh0 + 1286);
    const auto *snh0_1287 = buffer.data(snh0 + 1287);
    const auto *snh0_1290 = buffer.data(snh0 + 1290);
    const auto *snh0_1291 = buffer.data(snh0 + 1291);
    const auto *snh0_1293 = buffer.data(snh0 + 1293);
    const auto *snh0_1295 = buffer.data(snh0 + 1295);
    const auto *snh0_1296 = buffer.data(snh0 + 1296);
    const auto *snh0_1298 = buffer.data(snh0 + 1298);
    const auto *snh0_1299 = buffer.data(snh0 + 1299);
    const auto *snh0_1300 = buffer.data(snh0 + 1300);
    const auto *snh0_1301 = buffer.data(snh0 + 1301);
    const auto *snh0_1302 = buffer.data(snh0 + 1302);
    const auto *snh0_1305 = buffer.data(snh0 + 1305);
    const auto *snh0_1307 = buffer.data(snh0 + 1307);
    const auto *snh0_1308 = buffer.data(snh0 + 1308);
    const auto *snh0_1311 = buffer.data(snh0 + 1311);
    const auto *snh0_1312 = buffer.data(snh0 + 1312);
    const auto *snh0_1314 = buffer.data(snh0 + 1314);
    const auto *snh0_1316 = buffer.data(snh0 + 1316);
    const auto *snh0_1317 = buffer.data(snh0 + 1317);
    const auto *snh0_1319 = buffer.data(snh0 + 1319);
    const auto *snh0_1320 = buffer.data(snh0 + 1320);
    const auto *snh0_1322 = buffer.data(snh0 + 1322);

    const auto *snh1_1254 = buffer.data(snh1 + 1254);
    const auto *snh1_1256 = buffer.data(snh1 + 1256);
    const auto *snh1_1257 = buffer.data(snh1 + 1257);
    const auto *snh1_1258 = buffer.data(snh1 + 1258);
    const auto *snh1_1259 = buffer.data(snh1 + 1259);
    const auto *snh1_1260 = buffer.data(snh1 + 1260);
    const auto *snh1_1263 = buffer.data(snh1 + 1263);
    const auto *snh1_1265 = buffer.data(snh1 + 1265);
    const auto *snh1_1266 = buffer.data(snh1 + 1266);
    const auto *snh1_1269 = buffer.data(snh1 + 1269);
    const auto *snh1_1270 = buffer.data(snh1 + 1270);
    const auto *snh1_1272 = buffer.data(snh1 + 1272);
    const auto *snh1_1274 = buffer.data(snh1 + 1274);
    const auto *snh1_1275 = buffer.data(snh1 + 1275);
    const auto *snh1_1277 = buffer.data(snh1 + 1277);
    const auto *snh1_1278 = buffer.data(snh1 + 1278);
    const auto *snh1_1279 = buffer.data(snh1 + 1279);
    const auto *snh1_1280 = buffer.data(snh1 + 1280);
    const auto *snh1_1281 = buffer.data(snh1 + 1281);
    const auto *snh1_1284 = buffer.data(snh1 + 1284);
    const auto *snh1_1286 = buffer.data(snh1 + 1286);
    const auto *snh1_1287 = buffer.data(snh1 + 1287);
    const auto *snh1_1290 = buffer.data(snh1 + 1290);
    const auto *snh1_1291 = buffer.data(snh1 + 1291);
    const auto *snh1_1293 = buffer.data(snh1 + 1293);
    const auto *snh1_1295 = buffer.data(snh1 + 1295);
    const auto *snh1_1296 = buffer.data(snh1 + 1296);
    const auto *snh1_1298 = buffer.data(snh1 + 1298);
    const auto *snh1_1299 = buffer.data(snh1 + 1299);
    const auto *snh1_1300 = buffer.data(snh1 + 1300);
    const auto *snh1_1301 = buffer.data(snh1 + 1301);
    const auto *snh1_1302 = buffer.data(snh1 + 1302);
    const auto *snh1_1305 = buffer.data(snh1 + 1305);
    const auto *snh1_1307 = buffer.data(snh1 + 1307);
    const auto *snh1_1308 = buffer.data(snh1 + 1308);
    const auto *snh1_1311 = buffer.data(snh1 + 1311);
    const auto *snh1_1312 = buffer.data(snh1 + 1312);
    const auto *snh1_1314 = buffer.data(snh1 + 1314);
    const auto *snh1_1316 = buffer.data(snh1 + 1316);
    const auto *snh1_1317 = buffer.data(snh1 + 1317);
    const auto *snh1_1319 = buffer.data(snh1 + 1319);
    const auto *snh1_1320 = buffer.data(snh1 + 1320);
    const auto *snh1_1322 = buffer.data(snh1 + 1322);

    const auto *sni_1662 = buffer.data(sni + 1662);
    const auto *sni_1666 = buffer.data(sni + 1666);
    const auto *sni_1667 = buffer.data(sni + 1667);
    const auto *sni_1669 = buffer.data(sni + 1669);
    const auto *sni_1670 = buffer.data(sni + 1670);
    const auto *sni_1672 = buffer.data(sni + 1672);
    const auto *sni_1673 = buffer.data(sni + 1673);
    const auto *sni_1674 = buffer.data(sni + 1674);
    const auto *sni_1675 = buffer.data(sni + 1675);
    const auto *sni_1676 = buffer.data(sni + 1676);
    const auto *sni_1677 = buffer.data(sni + 1677);
    const auto *sni_1678 = buffer.data(sni + 1678);
    const auto *sni_1679 = buffer.data(sni + 1679);
    const auto *sni_1680 = buffer.data(sni + 1680);
    const auto *sni_1682 = buffer.data(sni + 1682);
    const auto *sni_1683 = buffer.data(sni + 1683);
    const auto *sni_1685 = buffer.data(sni + 1685);
    const auto *sni_1686 = buffer.data(sni + 1686);
    const auto *sni_1689 = buffer.data(sni + 1689);
    const auto *sni_1690 = buffer.data(sni + 1690);
    const auto *sni_1692 = buffer.data(sni + 1692);
    const auto *sni_1694 = buffer.data(sni + 1694);
    const auto *sni_1695 = buffer.data(sni + 1695);
    const auto *sni_1697 = buffer.data(sni + 1697);
    const auto *sni_1698 = buffer.data(sni + 1698);
    const auto *sni_1700 = buffer.data(sni + 1700);
    const auto *sni_1701 = buffer.data(sni + 1701);
    const auto *sni_1702 = buffer.data(sni + 1702);
    const auto *sni_1703 = buffer.data(sni + 1703);
    const auto *sni_1704 = buffer.data(sni + 1704);
    const auto *sni_1705 = buffer.data(sni + 1705);
    const auto *sni_1706 = buffer.data(sni + 1706);
    const auto *sni_1707 = buffer.data(sni + 1707);
    const auto *sni_1708 = buffer.data(sni + 1708);
    const auto *sni_1710 = buffer.data(sni + 1710);
    const auto *sni_1711 = buffer.data(sni + 1711);
    const auto *sni_1713 = buffer.data(sni + 1713);
    const auto *sni_1714 = buffer.data(sni + 1714);
    const auto *sni_1717 = buffer.data(sni + 1717);
    const auto *sni_1718 = buffer.data(sni + 1718);
    const auto *sni_1720 = buffer.data(sni + 1720);
    const auto *sni_1722 = buffer.data(sni + 1722);
    const auto *sni_1723 = buffer.data(sni + 1723);
    const auto *sni_1725 = buffer.data(sni + 1725);
    const auto *sni_1726 = buffer.data(sni + 1726);
    const auto *sni_1728 = buffer.data(sni + 1728);
    const auto *sni_1729 = buffer.data(sni + 1729);
    const auto *sni_1730 = buffer.data(sni + 1730);
    const auto *sni_1731 = buffer.data(sni + 1731);
    const auto *sni_1732 = buffer.data(sni + 1732);
    const auto *sni_1733 = buffer.data(sni + 1733);
    const auto *sni_1734 = buffer.data(sni + 1734);
    const auto *sni_1735 = buffer.data(sni + 1735);
    const auto *sni_1736 = buffer.data(sni + 1736);
    const auto *sni_1738 = buffer.data(sni + 1738);
    const auto *sni_1739 = buffer.data(sni + 1739);
    const auto *sni_1741 = buffer.data(sni + 1741);
    const auto *sni_1742 = buffer.data(sni + 1742);
    const auto *sni_1745 = buffer.data(sni + 1745);
    const auto *sni_1746 = buffer.data(sni + 1746);
    const auto *sni_1748 = buffer.data(sni + 1748);
    const auto *sni_1750 = buffer.data(sni + 1750);
    const auto *sni_1751 = buffer.data(sni + 1751);
    const auto *sni_1753 = buffer.data(sni + 1753);
    const auto *sni_1754 = buffer.data(sni + 1754);
    const auto *sni_1756 = buffer.data(sni + 1756);
    const auto *sni_1757 = buffer.data(sni + 1757);
    const auto *sni_1758 = buffer.data(sni + 1758);
    const auto *sni_1759 = buffer.data(sni + 1759);
    const auto *sni_1760 = buffer.data(sni + 1760);
    const auto *sni_1761 = buffer.data(sni + 1761);
    const auto *sni_1762 = buffer.data(sni + 1762);
    const auto *sni_1763 = buffer.data(sni + 1763);

#pragma omp simd aligned(t_2139, t_2140, t_2141, pc_x, pc_z, smi_1354, snh0_1254, snh0_1256, \
                         snh1_1254, snh1_1256, sni_1662, sni_1667, \
                         sni_1669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2139[k] = f_10 * snh0_1254[k]
                    - f_11 * snh1_1254[k]
                    + f_3 * pc_x[k] * sni_1667[k];

        t_2140[k] = f_16 * smi_1354[k]
                    + f_3 * pc_z[k] * sni_1662[k];

        t_2141[k] = f_10 * snh0_1256[k]
                    - f_11 * snh1_1256[k]
                    + f_3 * pc_x[k] * sni_1669[k];
    }

#pragma omp simd aligned(t_2142, t_2143, t_2144, t_2145, pc_x, pc_y, smi_1386, snh0_1257, \
                         snh0_1259, snh1_1257, snh1_1259, sni_1666, sni_1670, sni_1672, \
                         sni_1673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2142[k] = f_10 * snh0_1257[k]
                    - f_11 * snh1_1257[k]
                    + f_3 * pc_x[k] * sni_1670[k];

        t_2143[k] = f_21 * smi_1386[k]
                    + f_3 * pc_y[k] * sni_1666[k];

        t_2144[k] = f_10 * snh0_1259[k]
                    - f_11 * snh1_1259[k]
                    + f_3 * pc_x[k] * sni_1672[k];

        t_2145[k] = f_3 * pc_x[k] * sni_1673[k];
    }

#pragma omp simd aligned(t_2146, t_2147, t_2148, t_2149, t_2150, t_2151, pc_x, sni_1674, \
                         sni_1675, sni_1676, sni_1677, sni_1678, \
                         sni_1679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2146[k] = f_3 * pc_x[k] * sni_1674[k];

        t_2147[k] = f_3 * pc_x[k] * sni_1675[k];

        t_2148[k] = f_3 * pc_x[k] * sni_1676[k];

        t_2149[k] = f_3 * pc_x[k] * sni_1677[k];

        t_2150[k] = f_3 * pc_x[k] * sni_1678[k];

        t_2151[k] = f_3 * pc_x[k] * sni_1679[k];
    }

#pragma omp simd aligned(t_2152, t_2153, t_2154, pc_y, pc_z, smi_1365, smi_1393, smi_1395, \
                         snh0_1254, snh0_1256, snh1_1254, snh1_1256, sni_1673, \
                         sni_1675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2152[k] = f_21 * smi_1393[k]
                    + f_1 * snh0_1254[k]
                    - f_2 * snh1_1254[k]
                    + f_3 * pc_y[k] * sni_1673[k];

        t_2153[k] = f_16 * smi_1365[k]
                    + f_3 * pc_z[k] * sni_1673[k];

        t_2154[k] = f_21 * smi_1395[k]
                    + f_4 * snh0_1256[k]
                    - f_5 * snh1_1256[k]
                    + f_3 * pc_y[k] * sni_1675[k];
    }

#pragma omp simd aligned(t_2155, t_2156, t_2157, pc_y, smi_1396, smi_1397, smi_1398, \
                         snh0_1257, snh0_1258, snh0_1259, snh1_1257, snh1_1258, snh1_1259, \
                         sni_1676, sni_1677, sni_1678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2155[k] = f_21 * smi_1396[k]
                    + f_6 * snh0_1257[k]
                    - f_7 * snh1_1257[k]
                    + f_3 * pc_y[k] * sni_1676[k];

        t_2156[k] = f_21 * smi_1397[k]
                    + f_8 * snh0_1258[k]
                    - f_9 * snh1_1258[k]
                    + f_3 * pc_y[k] * sni_1677[k];

        t_2157[k] = f_21 * smi_1398[k]
                    + f_10 * snh0_1259[k]
                    - f_11 * snh1_1259[k]
                    + f_3 * pc_y[k] * sni_1678[k];
    }

#pragma omp simd aligned(t_2158, t_2159, t_2160, t_2161, pc_x, pc_y, pc_z, smi_1371, smi_1399, \
                         smi_1400, snh0_1259, snh0_1260, snh1_1259, snh1_1260, sni_1679, \
                         sni_1680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2158[k] = f_21 * smi_1399[k]
                    + f_3 * pc_y[k] * sni_1679[k];

        t_2159[k] = f_16 * smi_1371[k]
                    + f_1 * snh0_1259[k]
                    - f_2 * snh1_1259[k]
                    + f_3 * pc_z[k] * sni_1679[k];

        t_2160[k] = f_1 * snh0_1260[k]
                    - f_2 * snh1_1260[k]
                    + f_3 * pc_x[k] * sni_1680[k];

        t_2161[k] = f_17 * smi_1400[k]
                    + f_3 * pc_y[k] * sni_1680[k];
    }

#pragma omp simd aligned(t_2162, t_2163, t_2164, pc_x, pc_y, pc_z, smi_1372, smi_1402, \
                         snh0_1263, snh1_1263, sni_1680, sni_1682, \
                         sni_1683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2162[k] = f_17 * smi_1372[k]
                    + f_3 * pc_z[k] * sni_1680[k];

        t_2163[k] = f_4 * snh0_1263[k]
                    - f_5 * snh1_1263[k]
                    + f_3 * pc_x[k] * sni_1683[k];

        t_2164[k] = f_17 * smi_1402[k]
                    + f_3 * pc_y[k] * sni_1682[k];
    }

#pragma omp simd aligned(t_2165, t_2166, t_2167, t_2168, pc_x, pc_y, pc_z, smi_1375, smi_1405, \
                         snh0_1265, snh0_1266, snh1_1265, snh1_1266, sni_1683, sni_1685, \
                         sni_1686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2165[k] = f_4 * snh0_1265[k]
                    - f_5 * snh1_1265[k]
                    + f_3 * pc_x[k] * sni_1685[k];

        t_2166[k] = f_6 * snh0_1266[k]
                    - f_7 * snh1_1266[k]
                    + f_3 * pc_x[k] * sni_1686[k];

        t_2167[k] = f_17 * smi_1375[k]
                    + f_3 * pc_z[k] * sni_1683[k];

        t_2168[k] = f_17 * smi_1405[k]
                    + f_3 * pc_y[k] * sni_1685[k];
    }

#pragma omp simd aligned(t_2169, t_2170, t_2171, pc_x, pc_z, smi_1378, snh0_1269, snh0_1270, \
                         snh1_1269, snh1_1270, sni_1686, sni_1689, \
                         sni_1690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2169[k] = f_6 * snh0_1269[k]
                    - f_7 * snh1_1269[k]
                    + f_3 * pc_x[k] * sni_1689[k];

        t_2170[k] = f_8 * snh0_1270[k]
                    - f_9 * snh1_1270[k]
                    + f_3 * pc_x[k] * sni_1690[k];

        t_2171[k] = f_17 * smi_1378[k]
                    + f_3 * pc_z[k] * sni_1686[k];
    }

#pragma omp simd aligned(t_2172, t_2173, t_2174, pc_x, pc_y, smi_1409, snh0_1272, snh0_1274, \
                         snh1_1272, snh1_1274, sni_1689, sni_1692, \
                         sni_1694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2172[k] = f_8 * snh0_1272[k]
                    - f_9 * snh1_1272[k]
                    + f_3 * pc_x[k] * sni_1692[k];

        t_2173[k] = f_17 * smi_1409[k]
                    + f_3 * pc_y[k] * sni_1689[k];

        t_2174[k] = f_8 * snh0_1274[k]
                    - f_9 * snh1_1274[k]
                    + f_3 * pc_x[k] * sni_1694[k];
    }

#pragma omp simd aligned(t_2175, t_2176, t_2177, pc_x, pc_z, smi_1382, snh0_1275, snh0_1277, \
                         snh1_1275, snh1_1277, sni_1690, sni_1695, \
                         sni_1697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2175[k] = f_10 * snh0_1275[k]
                    - f_11 * snh1_1275[k]
                    + f_3 * pc_x[k] * sni_1695[k];

        t_2176[k] = f_17 * smi_1382[k]
                    + f_3 * pc_z[k] * sni_1690[k];

        t_2177[k] = f_10 * snh0_1277[k]
                    - f_11 * snh1_1277[k]
                    + f_3 * pc_x[k] * sni_1697[k];
    }

#pragma omp simd aligned(t_2178, t_2179, t_2180, t_2181, pc_x, pc_y, smi_1414, snh0_1278, \
                         snh0_1280, snh1_1278, snh1_1280, sni_1694, sni_1698, sni_1700, \
                         sni_1701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2178[k] = f_10 * snh0_1278[k]
                    - f_11 * snh1_1278[k]
                    + f_3 * pc_x[k] * sni_1698[k];

        t_2179[k] = f_17 * smi_1414[k]
                    + f_3 * pc_y[k] * sni_1694[k];

        t_2180[k] = f_10 * snh0_1280[k]
                    - f_11 * snh1_1280[k]
                    + f_3 * pc_x[k] * sni_1700[k];

        t_2181[k] = f_3 * pc_x[k] * sni_1701[k];
    }

#pragma omp simd aligned(t_2182, t_2183, t_2184, t_2185, t_2186, t_2187, pc_x, sni_1702, \
                         sni_1703, sni_1704, sni_1705, sni_1706, \
                         sni_1707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2182[k] = f_3 * pc_x[k] * sni_1702[k];

        t_2183[k] = f_3 * pc_x[k] * sni_1703[k];

        t_2184[k] = f_3 * pc_x[k] * sni_1704[k];

        t_2185[k] = f_3 * pc_x[k] * sni_1705[k];

        t_2186[k] = f_3 * pc_x[k] * sni_1706[k];

        t_2187[k] = f_3 * pc_x[k] * sni_1707[k];
    }

#pragma omp simd aligned(t_2188, t_2189, t_2190, pc_y, pc_z, smi_1393, smi_1421, smi_1423, \
                         snh0_1275, snh0_1277, snh1_1275, snh1_1277, sni_1701, \
                         sni_1703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2188[k] = f_17 * smi_1421[k]
                    + f_1 * snh0_1275[k]
                    - f_2 * snh1_1275[k]
                    + f_3 * pc_y[k] * sni_1701[k];

        t_2189[k] = f_17 * smi_1393[k]
                    + f_3 * pc_z[k] * sni_1701[k];

        t_2190[k] = f_17 * smi_1423[k]
                    + f_4 * snh0_1277[k]
                    - f_5 * snh1_1277[k]
                    + f_3 * pc_y[k] * sni_1703[k];
    }

#pragma omp simd aligned(t_2191, t_2192, t_2193, pc_y, smi_1424, smi_1425, smi_1426, \
                         snh0_1278, snh0_1279, snh0_1280, snh1_1278, snh1_1279, snh1_1280, \
                         sni_1704, sni_1705, sni_1706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2191[k] = f_17 * smi_1424[k]
                    + f_6 * snh0_1278[k]
                    - f_7 * snh1_1278[k]
                    + f_3 * pc_y[k] * sni_1704[k];

        t_2192[k] = f_17 * smi_1425[k]
                    + f_8 * snh0_1279[k]
                    - f_9 * snh1_1279[k]
                    + f_3 * pc_y[k] * sni_1705[k];

        t_2193[k] = f_17 * smi_1426[k]
                    + f_10 * snh0_1280[k]
                    - f_11 * snh1_1280[k]
                    + f_3 * pc_y[k] * sni_1706[k];
    }

#pragma omp simd aligned(t_2194, t_2195, t_2196, t_2197, pc_x, pc_y, pc_z, smi_1399, smi_1427, \
                         smi_1428, snh0_1280, snh0_1281, snh1_1280, snh1_1281, sni_1707, \
                         sni_1708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2194[k] = f_17 * smi_1427[k]
                    + f_3 * pc_y[k] * sni_1707[k];

        t_2195[k] = f_17 * smi_1399[k]
                    + f_1 * snh0_1280[k]
                    - f_2 * snh1_1280[k]
                    + f_3 * pc_z[k] * sni_1707[k];

        t_2196[k] = f_1 * snh0_1281[k]
                    - f_2 * snh1_1281[k]
                    + f_3 * pc_x[k] * sni_1708[k];

        t_2197[k] = f_16 * smi_1428[k]
                    + f_3 * pc_y[k] * sni_1708[k];
    }

#pragma omp simd aligned(t_2198, t_2199, t_2200, pc_x, pc_y, pc_z, smi_1400, smi_1430, \
                         snh0_1284, snh1_1284, sni_1708, sni_1710, \
                         sni_1711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2198[k] = f_21 * smi_1400[k]
                    + f_3 * pc_z[k] * sni_1708[k];

        t_2199[k] = f_4 * snh0_1284[k]
                    - f_5 * snh1_1284[k]
                    + f_3 * pc_x[k] * sni_1711[k];

        t_2200[k] = f_16 * smi_1430[k]
                    + f_3 * pc_y[k] * sni_1710[k];
    }

#pragma omp simd aligned(t_2201, t_2202, t_2203, t_2204, pc_x, pc_y, pc_z, smi_1403, smi_1433, \
                         snh0_1286, snh0_1287, snh1_1286, snh1_1287, sni_1711, sni_1713, \
                         sni_1714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2201[k] = f_4 * snh0_1286[k]
                    - f_5 * snh1_1286[k]
                    + f_3 * pc_x[k] * sni_1713[k];

        t_2202[k] = f_6 * snh0_1287[k]
                    - f_7 * snh1_1287[k]
                    + f_3 * pc_x[k] * sni_1714[k];

        t_2203[k] = f_21 * smi_1403[k]
                    + f_3 * pc_z[k] * sni_1711[k];

        t_2204[k] = f_16 * smi_1433[k]
                    + f_3 * pc_y[k] * sni_1713[k];
    }

#pragma omp simd aligned(t_2205, t_2206, t_2207, pc_x, pc_z, smi_1406, snh0_1290, snh0_1291, \
                         snh1_1290, snh1_1291, sni_1714, sni_1717, \
                         sni_1718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2205[k] = f_6 * snh0_1290[k]
                    - f_7 * snh1_1290[k]
                    + f_3 * pc_x[k] * sni_1717[k];

        t_2206[k] = f_8 * snh0_1291[k]
                    - f_9 * snh1_1291[k]
                    + f_3 * pc_x[k] * sni_1718[k];

        t_2207[k] = f_21 * smi_1406[k]
                    + f_3 * pc_z[k] * sni_1714[k];
    }

#pragma omp simd aligned(t_2208, t_2209, t_2210, pc_x, pc_y, smi_1437, snh0_1293, snh0_1295, \
                         snh1_1293, snh1_1295, sni_1717, sni_1720, \
                         sni_1722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2208[k] = f_8 * snh0_1293[k]
                    - f_9 * snh1_1293[k]
                    + f_3 * pc_x[k] * sni_1720[k];

        t_2209[k] = f_16 * smi_1437[k]
                    + f_3 * pc_y[k] * sni_1717[k];

        t_2210[k] = f_8 * snh0_1295[k]
                    - f_9 * snh1_1295[k]
                    + f_3 * pc_x[k] * sni_1722[k];
    }

#pragma omp simd aligned(t_2211, t_2212, t_2213, pc_x, pc_z, smi_1410, snh0_1296, snh0_1298, \
                         snh1_1296, snh1_1298, sni_1718, sni_1723, \
                         sni_1725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2211[k] = f_10 * snh0_1296[k]
                    - f_11 * snh1_1296[k]
                    + f_3 * pc_x[k] * sni_1723[k];

        t_2212[k] = f_21 * smi_1410[k]
                    + f_3 * pc_z[k] * sni_1718[k];

        t_2213[k] = f_10 * snh0_1298[k]
                    - f_11 * snh1_1298[k]
                    + f_3 * pc_x[k] * sni_1725[k];
    }

#pragma omp simd aligned(t_2214, t_2215, t_2216, t_2217, pc_x, pc_y, smi_1442, snh0_1299, \
                         snh0_1301, snh1_1299, snh1_1301, sni_1722, sni_1726, sni_1728, \
                         sni_1729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2214[k] = f_10 * snh0_1299[k]
                    - f_11 * snh1_1299[k]
                    + f_3 * pc_x[k] * sni_1726[k];

        t_2215[k] = f_16 * smi_1442[k]
                    + f_3 * pc_y[k] * sni_1722[k];

        t_2216[k] = f_10 * snh0_1301[k]
                    - f_11 * snh1_1301[k]
                    + f_3 * pc_x[k] * sni_1728[k];

        t_2217[k] = f_3 * pc_x[k] * sni_1729[k];
    }

#pragma omp simd aligned(t_2218, t_2219, t_2220, t_2221, t_2222, t_2223, pc_x, sni_1730, \
                         sni_1731, sni_1732, sni_1733, sni_1734, \
                         sni_1735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2218[k] = f_3 * pc_x[k] * sni_1730[k];

        t_2219[k] = f_3 * pc_x[k] * sni_1731[k];

        t_2220[k] = f_3 * pc_x[k] * sni_1732[k];

        t_2221[k] = f_3 * pc_x[k] * sni_1733[k];

        t_2222[k] = f_3 * pc_x[k] * sni_1734[k];

        t_2223[k] = f_3 * pc_x[k] * sni_1735[k];
    }

#pragma omp simd aligned(t_2224, t_2225, t_2226, pc_y, pc_z, smi_1421, smi_1449, smi_1451, \
                         snh0_1296, snh0_1298, snh1_1296, snh1_1298, sni_1729, \
                         sni_1731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2224[k] = f_16 * smi_1449[k]
                    + f_1 * snh0_1296[k]
                    - f_2 * snh1_1296[k]
                    + f_3 * pc_y[k] * sni_1729[k];

        t_2225[k] = f_21 * smi_1421[k]
                    + f_3 * pc_z[k] * sni_1729[k];

        t_2226[k] = f_16 * smi_1451[k]
                    + f_4 * snh0_1298[k]
                    - f_5 * snh1_1298[k]
                    + f_3 * pc_y[k] * sni_1731[k];
    }

#pragma omp simd aligned(t_2227, t_2228, t_2229, pc_y, smi_1452, smi_1453, smi_1454, \
                         snh0_1299, snh0_1300, snh0_1301, snh1_1299, snh1_1300, snh1_1301, \
                         sni_1732, sni_1733, sni_1734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2227[k] = f_16 * smi_1452[k]
                    + f_6 * snh0_1299[k]
                    - f_7 * snh1_1299[k]
                    + f_3 * pc_y[k] * sni_1732[k];

        t_2228[k] = f_16 * smi_1453[k]
                    + f_8 * snh0_1300[k]
                    - f_9 * snh1_1300[k]
                    + f_3 * pc_y[k] * sni_1733[k];

        t_2229[k] = f_16 * smi_1454[k]
                    + f_10 * snh0_1301[k]
                    - f_11 * snh1_1301[k]
                    + f_3 * pc_y[k] * sni_1734[k];
    }

#pragma omp simd aligned(t_2230, t_2231, t_2232, t_2233, pc_x, pc_y, pc_z, smi_1427, smi_1455, \
                         smi_1456, snh0_1301, snh0_1302, snh1_1301, snh1_1302, sni_1735, \
                         sni_1736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2230[k] = f_16 * smi_1455[k]
                    + f_3 * pc_y[k] * sni_1735[k];

        t_2231[k] = f_21 * smi_1427[k]
                    + f_1 * snh0_1301[k]
                    - f_2 * snh1_1301[k]
                    + f_3 * pc_z[k] * sni_1735[k];

        t_2232[k] = f_1 * snh0_1302[k]
                    - f_2 * snh1_1302[k]
                    + f_3 * pc_x[k] * sni_1736[k];

        t_2233[k] = f_15 * smi_1456[k]
                    + f_3 * pc_y[k] * sni_1736[k];
    }

#pragma omp simd aligned(t_2234, t_2235, t_2236, pc_x, pc_y, pc_z, smi_1428, smi_1458, \
                         snh0_1305, snh1_1305, sni_1736, sni_1738, \
                         sni_1739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2234[k] = f_20 * smi_1428[k]
                    + f_3 * pc_z[k] * sni_1736[k];

        t_2235[k] = f_4 * snh0_1305[k]
                    - f_5 * snh1_1305[k]
                    + f_3 * pc_x[k] * sni_1739[k];

        t_2236[k] = f_15 * smi_1458[k]
                    + f_3 * pc_y[k] * sni_1738[k];
    }

#pragma omp simd aligned(t_2237, t_2238, t_2239, t_2240, pc_x, pc_y, pc_z, smi_1431, smi_1461, \
                         snh0_1307, snh0_1308, snh1_1307, snh1_1308, sni_1739, sni_1741, \
                         sni_1742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2237[k] = f_4 * snh0_1307[k]
                    - f_5 * snh1_1307[k]
                    + f_3 * pc_x[k] * sni_1741[k];

        t_2238[k] = f_6 * snh0_1308[k]
                    - f_7 * snh1_1308[k]
                    + f_3 * pc_x[k] * sni_1742[k];

        t_2239[k] = f_20 * smi_1431[k]
                    + f_3 * pc_z[k] * sni_1739[k];

        t_2240[k] = f_15 * smi_1461[k]
                    + f_3 * pc_y[k] * sni_1741[k];
    }

#pragma omp simd aligned(t_2241, t_2242, t_2243, pc_x, pc_z, smi_1434, snh0_1311, snh0_1312, \
                         snh1_1311, snh1_1312, sni_1742, sni_1745, \
                         sni_1746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2241[k] = f_6 * snh0_1311[k]
                    - f_7 * snh1_1311[k]
                    + f_3 * pc_x[k] * sni_1745[k];

        t_2242[k] = f_8 * snh0_1312[k]
                    - f_9 * snh1_1312[k]
                    + f_3 * pc_x[k] * sni_1746[k];

        t_2243[k] = f_20 * smi_1434[k]
                    + f_3 * pc_z[k] * sni_1742[k];
    }

#pragma omp simd aligned(t_2244, t_2245, t_2246, pc_x, pc_y, smi_1465, snh0_1314, snh0_1316, \
                         snh1_1314, snh1_1316, sni_1745, sni_1748, \
                         sni_1750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2244[k] = f_8 * snh0_1314[k]
                    - f_9 * snh1_1314[k]
                    + f_3 * pc_x[k] * sni_1748[k];

        t_2245[k] = f_15 * smi_1465[k]
                    + f_3 * pc_y[k] * sni_1745[k];

        t_2246[k] = f_8 * snh0_1316[k]
                    - f_9 * snh1_1316[k]
                    + f_3 * pc_x[k] * sni_1750[k];
    }

#pragma omp simd aligned(t_2247, t_2248, t_2249, pc_x, pc_z, smi_1438, snh0_1317, snh0_1319, \
                         snh1_1317, snh1_1319, sni_1746, sni_1751, \
                         sni_1753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2247[k] = f_10 * snh0_1317[k]
                    - f_11 * snh1_1317[k]
                    + f_3 * pc_x[k] * sni_1751[k];

        t_2248[k] = f_20 * smi_1438[k]
                    + f_3 * pc_z[k] * sni_1746[k];

        t_2249[k] = f_10 * snh0_1319[k]
                    - f_11 * snh1_1319[k]
                    + f_3 * pc_x[k] * sni_1753[k];
    }

#pragma omp simd aligned(t_2250, t_2251, t_2252, t_2253, pc_x, pc_y, smi_1470, snh0_1320, \
                         snh0_1322, snh1_1320, snh1_1322, sni_1750, sni_1754, sni_1756, \
                         sni_1757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2250[k] = f_10 * snh0_1320[k]
                    - f_11 * snh1_1320[k]
                    + f_3 * pc_x[k] * sni_1754[k];

        t_2251[k] = f_15 * smi_1470[k]
                    + f_3 * pc_y[k] * sni_1750[k];

        t_2252[k] = f_10 * snh0_1322[k]
                    - f_11 * snh1_1322[k]
                    + f_3 * pc_x[k] * sni_1756[k];

        t_2253[k] = f_3 * pc_x[k] * sni_1757[k];
    }

#pragma omp simd aligned(t_2254, t_2255, t_2256, t_2257, t_2258, t_2259, pc_x, sni_1758, \
                         sni_1759, sni_1760, sni_1761, sni_1762, \
                         sni_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2254[k] = f_3 * pc_x[k] * sni_1758[k];

        t_2255[k] = f_3 * pc_x[k] * sni_1759[k];

        t_2256[k] = f_3 * pc_x[k] * sni_1760[k];

        t_2257[k] = f_3 * pc_x[k] * sni_1761[k];

        t_2258[k] = f_3 * pc_x[k] * sni_1762[k];

        t_2259[k] = f_3 * pc_x[k] * sni_1763[k];
    }
}

static auto
compute_prim_snk_three_center_electron_repulsion_0_piece20(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t smk0,
                                                           const size_t smi, const size_t smk1,
                                                           const size_t snh0, const size_t snh1,
                                                           const size_t sni, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
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
    const auto f_18 = 4.5 / q;
    const auto f_19 = 4.0 / q;
    const auto f_20 = 3.5 / q;

    auto *t_2260 = buffer.data(target + 2260);
    auto *t_2261 = buffer.data(target + 2261);
    auto *t_2262 = buffer.data(target + 2262);
    auto *t_2263 = buffer.data(target + 2263);
    auto *t_2264 = buffer.data(target + 2264);
    auto *t_2265 = buffer.data(target + 2265);
    auto *t_2266 = buffer.data(target + 2266);
    auto *t_2267 = buffer.data(target + 2267);
    auto *t_2268 = buffer.data(target + 2268);
    auto *t_2269 = buffer.data(target + 2269);
    auto *t_2270 = buffer.data(target + 2270);
    auto *t_2271 = buffer.data(target + 2271);
    auto *t_2272 = buffer.data(target + 2272);
    auto *t_2273 = buffer.data(target + 2273);
    auto *t_2274 = buffer.data(target + 2274);
    auto *t_2275 = buffer.data(target + 2275);
    auto *t_2276 = buffer.data(target + 2276);
    auto *t_2277 = buffer.data(target + 2277);
    auto *t_2278 = buffer.data(target + 2278);
    auto *t_2279 = buffer.data(target + 2279);
    auto *t_2280 = buffer.data(target + 2280);
    auto *t_2281 = buffer.data(target + 2281);
    auto *t_2282 = buffer.data(target + 2282);
    auto *t_2283 = buffer.data(target + 2283);
    auto *t_2284 = buffer.data(target + 2284);
    auto *t_2285 = buffer.data(target + 2285);
    auto *t_2286 = buffer.data(target + 2286);
    auto *t_2287 = buffer.data(target + 2287);
    auto *t_2288 = buffer.data(target + 2288);
    auto *t_2289 = buffer.data(target + 2289);
    auto *t_2290 = buffer.data(target + 2290);
    auto *t_2291 = buffer.data(target + 2291);
    auto *t_2292 = buffer.data(target + 2292);
    auto *t_2293 = buffer.data(target + 2293);
    auto *t_2294 = buffer.data(target + 2294);
    auto *t_2295 = buffer.data(target + 2295);
    auto *t_2296 = buffer.data(target + 2296);
    auto *t_2297 = buffer.data(target + 2297);
    auto *t_2298 = buffer.data(target + 2298);
    auto *t_2299 = buffer.data(target + 2299);
    auto *t_2300 = buffer.data(target + 2300);
    auto *t_2301 = buffer.data(target + 2301);
    auto *t_2302 = buffer.data(target + 2302);
    auto *t_2303 = buffer.data(target + 2303);
    auto *t_2304 = buffer.data(target + 2304);
    auto *t_2305 = buffer.data(target + 2305);
    auto *t_2306 = buffer.data(target + 2306);
    auto *t_2307 = buffer.data(target + 2307);
    auto *t_2308 = buffer.data(target + 2308);
    auto *t_2309 = buffer.data(target + 2309);
    auto *t_2310 = buffer.data(target + 2310);
    auto *t_2311 = buffer.data(target + 2311);
    auto *t_2312 = buffer.data(target + 2312);
    auto *t_2313 = buffer.data(target + 2313);
    auto *t_2314 = buffer.data(target + 2314);
    auto *t_2315 = buffer.data(target + 2315);
    auto *t_2316 = buffer.data(target + 2316);
    auto *t_2317 = buffer.data(target + 2317);
    auto *t_2318 = buffer.data(target + 2318);
    auto *t_2319 = buffer.data(target + 2319);
    auto *t_2320 = buffer.data(target + 2320);
    auto *t_2321 = buffer.data(target + 2321);
    auto *t_2322 = buffer.data(target + 2322);
    auto *t_2323 = buffer.data(target + 2323);
    auto *t_2324 = buffer.data(target + 2324);
    auto *t_2325 = buffer.data(target + 2325);
    auto *t_2326 = buffer.data(target + 2326);
    auto *t_2327 = buffer.data(target + 2327);
    auto *t_2328 = buffer.data(target + 2328);
    auto *t_2329 = buffer.data(target + 2329);
    auto *t_2330 = buffer.data(target + 2330);
    auto *t_2331 = buffer.data(target + 2331);
    auto *t_2332 = buffer.data(target + 2332);
    auto *t_2333 = buffer.data(target + 2333);
    auto *t_2334 = buffer.data(target + 2334);
    auto *t_2335 = buffer.data(target + 2335);
    auto *t_2336 = buffer.data(target + 2336);
    auto *t_2337 = buffer.data(target + 2337);
    auto *t_2338 = buffer.data(target + 2338);
    auto *t_2339 = buffer.data(target + 2339);
    auto *t_2340 = buffer.data(target + 2340);
    auto *t_2341 = buffer.data(target + 2341);
    auto *t_2342 = buffer.data(target + 2342);
    auto *t_2343 = buffer.data(target + 2343);
    auto *t_2344 = buffer.data(target + 2344);
    auto *t_2345 = buffer.data(target + 2345);
    auto *t_2346 = buffer.data(target + 2346);
    auto *t_2347 = buffer.data(target + 2347);
    auto *t_2348 = buffer.data(target + 2348);
    auto *t_2349 = buffer.data(target + 2349);
    auto *t_2350 = buffer.data(target + 2350);
    auto *t_2351 = buffer.data(target + 2351);
    auto *t_2352 = buffer.data(target + 2352);
    auto *t_2353 = buffer.data(target + 2353);
    auto *t_2354 = buffer.data(target + 2354);
    auto *t_2355 = buffer.data(target + 2355);
    auto *t_2356 = buffer.data(target + 2356);
    auto *t_2357 = buffer.data(target + 2357);
    auto *t_2358 = buffer.data(target + 2358);
    auto *t_2359 = buffer.data(target + 2359);
    auto *t_2360 = buffer.data(target + 2360);
    auto *t_2361 = buffer.data(target + 2361);
    auto *t_2362 = buffer.data(target + 2362);
    auto *t_2363 = buffer.data(target + 2363);
    auto *t_2364 = buffer.data(target + 2364);
    auto *t_2365 = buffer.data(target + 2365);
    auto *t_2366 = buffer.data(target + 2366);
    auto *t_2367 = buffer.data(target + 2367);
    auto *t_2368 = buffer.data(target + 2368);
    auto *t_2369 = buffer.data(target + 2369);
    auto *t_2370 = buffer.data(target + 2370);
    auto *t_2371 = buffer.data(target + 2371);
    auto *t_2372 = buffer.data(target + 2372);
    auto *t_2373 = buffer.data(target + 2373);
    auto *t_2374 = buffer.data(target + 2374);
    auto *t_2375 = buffer.data(target + 2375);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smk0_1944 = buffer.data(smk0 + 1944);
    const auto *smk0_1949 = buffer.data(smk0 + 1949);
    const auto *smk0_1953 = buffer.data(smk0 + 1953);
    const auto *smk0_1958 = buffer.data(smk0 + 1958);
    const auto *smk0_1964 = buffer.data(smk0 + 1964);
    const auto *smk0_1972 = buffer.data(smk0 + 1972);
    const auto *smk0_1974 = buffer.data(smk0 + 1974);
    const auto *smk0_1975 = buffer.data(smk0 + 1975);
    const auto *smk0_1976 = buffer.data(smk0 + 1976);
    const auto *smk0_1977 = buffer.data(smk0 + 1977);
    const auto *smk0_1979 = buffer.data(smk0 + 1979);

    const auto *smi_1449 = buffer.data(smi + 1449);
    const auto *smi_1455 = buffer.data(smi + 1455);
    const auto *smi_1456 = buffer.data(smi + 1456);
    const auto *smi_1459 = buffer.data(smi + 1459);
    const auto *smi_1462 = buffer.data(smi + 1462);
    const auto *smi_1466 = buffer.data(smi + 1466);
    const auto *smi_1477 = buffer.data(smi + 1477);
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
    const auto *smi_1498 = buffer.data(smi + 1498);
    const auto *smi_1505 = buffer.data(smi + 1505);
    const auto *smi_1507 = buffer.data(smi + 1507);
    const auto *smi_1508 = buffer.data(smi + 1508);
    const auto *smi_1509 = buffer.data(smi + 1509);
    const auto *smi_1510 = buffer.data(smi + 1510);
    const auto *smi_1511 = buffer.data(smi + 1511);
    const auto *smi_1512 = buffer.data(smi + 1512);
    const auto *smi_1514 = buffer.data(smi + 1514);
    const auto *smi_1515 = buffer.data(smi + 1515);
    const auto *smi_1517 = buffer.data(smi + 1517);
    const auto *smi_1518 = buffer.data(smi + 1518);
    const auto *smi_1521 = buffer.data(smi + 1521);
    const auto *smi_1522 = buffer.data(smi + 1522);
    const auto *smi_1526 = buffer.data(smi + 1526);
    const auto *smi_1533 = buffer.data(smi + 1533);
    const auto *smi_1535 = buffer.data(smi + 1535);
    const auto *smi_1536 = buffer.data(smi + 1536);
    const auto *smi_1537 = buffer.data(smi + 1537);
    const auto *smi_1538 = buffer.data(smi + 1538);
    const auto *smi_1539 = buffer.data(smi + 1539);

    const auto *smk1_1944 = buffer.data(smk1 + 1944);
    const auto *smk1_1949 = buffer.data(smk1 + 1949);
    const auto *smk1_1953 = buffer.data(smk1 + 1953);
    const auto *smk1_1958 = buffer.data(smk1 + 1958);
    const auto *smk1_1964 = buffer.data(smk1 + 1964);
    const auto *smk1_1972 = buffer.data(smk1 + 1972);
    const auto *smk1_1974 = buffer.data(smk1 + 1974);
    const auto *smk1_1975 = buffer.data(smk1 + 1975);
    const auto *smk1_1976 = buffer.data(smk1 + 1976);
    const auto *smk1_1977 = buffer.data(smk1 + 1977);
    const auto *smk1_1979 = buffer.data(smk1 + 1979);

    const auto *snh0_1317 = buffer.data(snh0 + 1317);
    const auto *snh0_1319 = buffer.data(snh0 + 1319);
    const auto *snh0_1320 = buffer.data(snh0 + 1320);
    const auto *snh0_1321 = buffer.data(snh0 + 1321);
    const auto *snh0_1322 = buffer.data(snh0 + 1322);
    const auto *snh0_1323 = buffer.data(snh0 + 1323);
    const auto *snh0_1326 = buffer.data(snh0 + 1326);
    const auto *snh0_1328 = buffer.data(snh0 + 1328);
    const auto *snh0_1329 = buffer.data(snh0 + 1329);
    const auto *snh0_1332 = buffer.data(snh0 + 1332);
    const auto *snh0_1333 = buffer.data(snh0 + 1333);
    const auto *snh0_1335 = buffer.data(snh0 + 1335);
    const auto *snh0_1337 = buffer.data(snh0 + 1337);
    const auto *snh0_1338 = buffer.data(snh0 + 1338);
    const auto *snh0_1340 = buffer.data(snh0 + 1340);
    const auto *snh0_1341 = buffer.data(snh0 + 1341);
    const auto *snh0_1342 = buffer.data(snh0 + 1342);
    const auto *snh0_1343 = buffer.data(snh0 + 1343);
    const auto *snh0_1347 = buffer.data(snh0 + 1347);
    const auto *snh0_1350 = buffer.data(snh0 + 1350);
    const auto *snh0_1354 = buffer.data(snh0 + 1354);
    const auto *snh0_1356 = buffer.data(snh0 + 1356);
    const auto *snh0_1359 = buffer.data(snh0 + 1359);
    const auto *snh0_1361 = buffer.data(snh0 + 1361);
    const auto *snh0_1362 = buffer.data(snh0 + 1362);
    const auto *snh0_1365 = buffer.data(snh0 + 1365);
    const auto *snh0_1368 = buffer.data(snh0 + 1368);
    const auto *snh0_1370 = buffer.data(snh0 + 1370);
    const auto *snh0_1371 = buffer.data(snh0 + 1371);
    const auto *snh0_1374 = buffer.data(snh0 + 1374);
    const auto *snh0_1375 = buffer.data(snh0 + 1375);
    const auto *snh0_1377 = buffer.data(snh0 + 1377);
    const auto *snh0_1379 = buffer.data(snh0 + 1379);
    const auto *snh0_1380 = buffer.data(snh0 + 1380);
    const auto *snh0_1382 = buffer.data(snh0 + 1382);
    const auto *snh0_1383 = buffer.data(snh0 + 1383);
    const auto *snh0_1384 = buffer.data(snh0 + 1384);
    const auto *snh0_1385 = buffer.data(snh0 + 1385);

    const auto *snh1_1317 = buffer.data(snh1 + 1317);
    const auto *snh1_1319 = buffer.data(snh1 + 1319);
    const auto *snh1_1320 = buffer.data(snh1 + 1320);
    const auto *snh1_1321 = buffer.data(snh1 + 1321);
    const auto *snh1_1322 = buffer.data(snh1 + 1322);
    const auto *snh1_1323 = buffer.data(snh1 + 1323);
    const auto *snh1_1326 = buffer.data(snh1 + 1326);
    const auto *snh1_1328 = buffer.data(snh1 + 1328);
    const auto *snh1_1329 = buffer.data(snh1 + 1329);
    const auto *snh1_1332 = buffer.data(snh1 + 1332);
    const auto *snh1_1333 = buffer.data(snh1 + 1333);
    const auto *snh1_1335 = buffer.data(snh1 + 1335);
    const auto *snh1_1337 = buffer.data(snh1 + 1337);
    const auto *snh1_1338 = buffer.data(snh1 + 1338);
    const auto *snh1_1340 = buffer.data(snh1 + 1340);
    const auto *snh1_1341 = buffer.data(snh1 + 1341);
    const auto *snh1_1342 = buffer.data(snh1 + 1342);
    const auto *snh1_1343 = buffer.data(snh1 + 1343);
    const auto *snh1_1347 = buffer.data(snh1 + 1347);
    const auto *snh1_1350 = buffer.data(snh1 + 1350);
    const auto *snh1_1354 = buffer.data(snh1 + 1354);
    const auto *snh1_1356 = buffer.data(snh1 + 1356);
    const auto *snh1_1359 = buffer.data(snh1 + 1359);
    const auto *snh1_1361 = buffer.data(snh1 + 1361);
    const auto *snh1_1362 = buffer.data(snh1 + 1362);
    const auto *snh1_1365 = buffer.data(snh1 + 1365);
    const auto *snh1_1368 = buffer.data(snh1 + 1368);
    const auto *snh1_1370 = buffer.data(snh1 + 1370);
    const auto *snh1_1371 = buffer.data(snh1 + 1371);
    const auto *snh1_1374 = buffer.data(snh1 + 1374);
    const auto *snh1_1375 = buffer.data(snh1 + 1375);
    const auto *snh1_1377 = buffer.data(snh1 + 1377);
    const auto *snh1_1379 = buffer.data(snh1 + 1379);
    const auto *snh1_1380 = buffer.data(snh1 + 1380);
    const auto *snh1_1382 = buffer.data(snh1 + 1382);
    const auto *snh1_1383 = buffer.data(snh1 + 1383);
    const auto *snh1_1384 = buffer.data(snh1 + 1384);
    const auto *snh1_1385 = buffer.data(snh1 + 1385);

    const auto *sni_1757 = buffer.data(sni + 1757);
    const auto *sni_1759 = buffer.data(sni + 1759);
    const auto *sni_1760 = buffer.data(sni + 1760);
    const auto *sni_1761 = buffer.data(sni + 1761);
    const auto *sni_1762 = buffer.data(sni + 1762);
    const auto *sni_1763 = buffer.data(sni + 1763);
    const auto *sni_1764 = buffer.data(sni + 1764);
    const auto *sni_1766 = buffer.data(sni + 1766);
    const auto *sni_1767 = buffer.data(sni + 1767);
    const auto *sni_1769 = buffer.data(sni + 1769);
    const auto *sni_1770 = buffer.data(sni + 1770);
    const auto *sni_1773 = buffer.data(sni + 1773);
    const auto *sni_1774 = buffer.data(sni + 1774);
    const auto *sni_1776 = buffer.data(sni + 1776);
    const auto *sni_1778 = buffer.data(sni + 1778);
    const auto *sni_1779 = buffer.data(sni + 1779);
    const auto *sni_1781 = buffer.data(sni + 1781);
    const auto *sni_1782 = buffer.data(sni + 1782);
    const auto *sni_1784 = buffer.data(sni + 1784);
    const auto *sni_1785 = buffer.data(sni + 1785);
    const auto *sni_1786 = buffer.data(sni + 1786);
    const auto *sni_1787 = buffer.data(sni + 1787);
    const auto *sni_1788 = buffer.data(sni + 1788);
    const auto *sni_1789 = buffer.data(sni + 1789);
    const auto *sni_1790 = buffer.data(sni + 1790);
    const auto *sni_1791 = buffer.data(sni + 1791);
    const auto *sni_1792 = buffer.data(sni + 1792);
    const auto *sni_1794 = buffer.data(sni + 1794);
    const auto *sni_1795 = buffer.data(sni + 1795);
    const auto *sni_1797 = buffer.data(sni + 1797);
    const auto *sni_1798 = buffer.data(sni + 1798);
    const auto *sni_1801 = buffer.data(sni + 1801);
    const auto *sni_1802 = buffer.data(sni + 1802);
    const auto *sni_1804 = buffer.data(sni + 1804);
    const auto *sni_1806 = buffer.data(sni + 1806);
    const auto *sni_1807 = buffer.data(sni + 1807);
    const auto *sni_1809 = buffer.data(sni + 1809);
    const auto *sni_1810 = buffer.data(sni + 1810);
    const auto *sni_1813 = buffer.data(sni + 1813);
    const auto *sni_1814 = buffer.data(sni + 1814);
    const auto *sni_1815 = buffer.data(sni + 1815);
    const auto *sni_1816 = buffer.data(sni + 1816);
    const auto *sni_1817 = buffer.data(sni + 1817);
    const auto *sni_1818 = buffer.data(sni + 1818);
    const auto *sni_1819 = buffer.data(sni + 1819);
    const auto *sni_1820 = buffer.data(sni + 1820);
    const auto *sni_1822 = buffer.data(sni + 1822);
    const auto *sni_1823 = buffer.data(sni + 1823);
    const auto *sni_1825 = buffer.data(sni + 1825);
    const auto *sni_1826 = buffer.data(sni + 1826);
    const auto *sni_1829 = buffer.data(sni + 1829);
    const auto *sni_1830 = buffer.data(sni + 1830);
    const auto *sni_1832 = buffer.data(sni + 1832);
    const auto *sni_1834 = buffer.data(sni + 1834);
    const auto *sni_1835 = buffer.data(sni + 1835);
    const auto *sni_1837 = buffer.data(sni + 1837);
    const auto *sni_1838 = buffer.data(sni + 1838);
    const auto *sni_1840 = buffer.data(sni + 1840);
    const auto *sni_1841 = buffer.data(sni + 1841);
    const auto *sni_1842 = buffer.data(sni + 1842);
    const auto *sni_1843 = buffer.data(sni + 1843);
    const auto *sni_1844 = buffer.data(sni + 1844);
    const auto *sni_1845 = buffer.data(sni + 1845);
    const auto *sni_1846 = buffer.data(sni + 1846);
    const auto *sni_1847 = buffer.data(sni + 1847);

#pragma omp simd aligned(t_2260, t_2261, t_2262, pc_y, pc_z, smi_1449, smi_1477, smi_1479, \
                         snh0_1317, snh0_1319, snh1_1317, snh1_1319, sni_1757, \
                         sni_1759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2260[k] = f_15 * smi_1477[k]
                    + f_1 * snh0_1317[k]
                    - f_2 * snh1_1317[k]
                    + f_3 * pc_y[k] * sni_1757[k];

        t_2261[k] = f_20 * smi_1449[k]
                    + f_3 * pc_z[k] * sni_1757[k];

        t_2262[k] = f_15 * smi_1479[k]
                    + f_4 * snh0_1319[k]
                    - f_5 * snh1_1319[k]
                    + f_3 * pc_y[k] * sni_1759[k];
    }

#pragma omp simd aligned(t_2263, t_2264, t_2265, pc_y, smi_1480, smi_1481, smi_1482, \
                         snh0_1320, snh0_1321, snh0_1322, snh1_1320, snh1_1321, snh1_1322, \
                         sni_1760, sni_1761, sni_1762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2263[k] = f_15 * smi_1480[k]
                    + f_6 * snh0_1320[k]
                    - f_7 * snh1_1320[k]
                    + f_3 * pc_y[k] * sni_1760[k];

        t_2264[k] = f_15 * smi_1481[k]
                    + f_8 * snh0_1321[k]
                    - f_9 * snh1_1321[k]
                    + f_3 * pc_y[k] * sni_1761[k];

        t_2265[k] = f_15 * smi_1482[k]
                    + f_10 * snh0_1322[k]
                    - f_11 * snh1_1322[k]
                    + f_3 * pc_y[k] * sni_1762[k];
    }

#pragma omp simd aligned(t_2266, t_2267, t_2268, t_2269, pc_x, pc_y, pc_z, smi_1455, smi_1483, \
                         smi_1484, snh0_1322, snh0_1323, snh1_1322, snh1_1323, sni_1763, \
                         sni_1764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2266[k] = f_15 * smi_1483[k]
                    + f_3 * pc_y[k] * sni_1763[k];

        t_2267[k] = f_20 * smi_1455[k]
                    + f_1 * snh0_1322[k]
                    - f_2 * snh1_1322[k]
                    + f_3 * pc_z[k] * sni_1763[k];

        t_2268[k] = f_1 * snh0_1323[k]
                    - f_2 * snh1_1323[k]
                    + f_3 * pc_x[k] * sni_1764[k];

        t_2269[k] = f_14 * smi_1484[k]
                    + f_3 * pc_y[k] * sni_1764[k];
    }

#pragma omp simd aligned(t_2270, t_2271, t_2272, pc_x, pc_y, pc_z, smi_1456, smi_1486, \
                         snh0_1326, snh1_1326, sni_1764, sni_1766, \
                         sni_1767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2270[k] = f_19 * smi_1456[k]
                    + f_3 * pc_z[k] * sni_1764[k];

        t_2271[k] = f_4 * snh0_1326[k]
                    - f_5 * snh1_1326[k]
                    + f_3 * pc_x[k] * sni_1767[k];

        t_2272[k] = f_14 * smi_1486[k]
                    + f_3 * pc_y[k] * sni_1766[k];
    }

#pragma omp simd aligned(t_2273, t_2274, t_2275, t_2276, pc_x, pc_y, pc_z, smi_1459, smi_1489, \
                         snh0_1328, snh0_1329, snh1_1328, snh1_1329, sni_1767, sni_1769, \
                         sni_1770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2273[k] = f_4 * snh0_1328[k]
                    - f_5 * snh1_1328[k]
                    + f_3 * pc_x[k] * sni_1769[k];

        t_2274[k] = f_6 * snh0_1329[k]
                    - f_7 * snh1_1329[k]
                    + f_3 * pc_x[k] * sni_1770[k];

        t_2275[k] = f_19 * smi_1459[k]
                    + f_3 * pc_z[k] * sni_1767[k];

        t_2276[k] = f_14 * smi_1489[k]
                    + f_3 * pc_y[k] * sni_1769[k];
    }

#pragma omp simd aligned(t_2277, t_2278, t_2279, pc_x, pc_z, smi_1462, snh0_1332, snh0_1333, \
                         snh1_1332, snh1_1333, sni_1770, sni_1773, \
                         sni_1774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2277[k] = f_6 * snh0_1332[k]
                    - f_7 * snh1_1332[k]
                    + f_3 * pc_x[k] * sni_1773[k];

        t_2278[k] = f_8 * snh0_1333[k]
                    - f_9 * snh1_1333[k]
                    + f_3 * pc_x[k] * sni_1774[k];

        t_2279[k] = f_19 * smi_1462[k]
                    + f_3 * pc_z[k] * sni_1770[k];
    }

#pragma omp simd aligned(t_2280, t_2281, t_2282, pc_x, pc_y, smi_1493, snh0_1335, snh0_1337, \
                         snh1_1335, snh1_1337, sni_1773, sni_1776, \
                         sni_1778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2280[k] = f_8 * snh0_1335[k]
                    - f_9 * snh1_1335[k]
                    + f_3 * pc_x[k] * sni_1776[k];

        t_2281[k] = f_14 * smi_1493[k]
                    + f_3 * pc_y[k] * sni_1773[k];

        t_2282[k] = f_8 * snh0_1337[k]
                    - f_9 * snh1_1337[k]
                    + f_3 * pc_x[k] * sni_1778[k];
    }

#pragma omp simd aligned(t_2283, t_2284, t_2285, pc_x, pc_z, smi_1466, snh0_1338, snh0_1340, \
                         snh1_1338, snh1_1340, sni_1774, sni_1779, \
                         sni_1781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2283[k] = f_10 * snh0_1338[k]
                    - f_11 * snh1_1338[k]
                    + f_3 * pc_x[k] * sni_1779[k];

        t_2284[k] = f_19 * smi_1466[k]
                    + f_3 * pc_z[k] * sni_1774[k];

        t_2285[k] = f_10 * snh0_1340[k]
                    - f_11 * snh1_1340[k]
                    + f_3 * pc_x[k] * sni_1781[k];
    }

#pragma omp simd aligned(t_2286, t_2287, t_2288, t_2289, pc_x, pc_y, smi_1498, snh0_1341, \
                         snh0_1343, snh1_1341, snh1_1343, sni_1778, sni_1782, sni_1784, \
                         sni_1785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2286[k] = f_10 * snh0_1341[k]
                    - f_11 * snh1_1341[k]
                    + f_3 * pc_x[k] * sni_1782[k];

        t_2287[k] = f_14 * smi_1498[k]
                    + f_3 * pc_y[k] * sni_1778[k];

        t_2288[k] = f_10 * snh0_1343[k]
                    - f_11 * snh1_1343[k]
                    + f_3 * pc_x[k] * sni_1784[k];

        t_2289[k] = f_3 * pc_x[k] * sni_1785[k];
    }

#pragma omp simd aligned(t_2290, t_2291, t_2292, t_2293, t_2294, t_2295, pc_x, sni_1786, \
                         sni_1787, sni_1788, sni_1789, sni_1790, \
                         sni_1791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2290[k] = f_3 * pc_x[k] * sni_1786[k];

        t_2291[k] = f_3 * pc_x[k] * sni_1787[k];

        t_2292[k] = f_3 * pc_x[k] * sni_1788[k];

        t_2293[k] = f_3 * pc_x[k] * sni_1789[k];

        t_2294[k] = f_3 * pc_x[k] * sni_1790[k];

        t_2295[k] = f_3 * pc_x[k] * sni_1791[k];
    }

#pragma omp simd aligned(t_2296, t_2297, t_2298, pc_y, pc_z, smi_1477, smi_1505, smi_1507, \
                         snh0_1338, snh0_1340, snh1_1338, snh1_1340, sni_1785, \
                         sni_1787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2296[k] = f_14 * smi_1505[k]
                    + f_1 * snh0_1338[k]
                    - f_2 * snh1_1338[k]
                    + f_3 * pc_y[k] * sni_1785[k];

        t_2297[k] = f_19 * smi_1477[k]
                    + f_3 * pc_z[k] * sni_1785[k];

        t_2298[k] = f_14 * smi_1507[k]
                    + f_4 * snh0_1340[k]
                    - f_5 * snh1_1340[k]
                    + f_3 * pc_y[k] * sni_1787[k];
    }

#pragma omp simd aligned(t_2299, t_2300, t_2301, pc_y, smi_1508, smi_1509, smi_1510, \
                         snh0_1341, snh0_1342, snh0_1343, snh1_1341, snh1_1342, snh1_1343, \
                         sni_1788, sni_1789, sni_1790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2299[k] = f_14 * smi_1508[k]
                    + f_6 * snh0_1341[k]
                    - f_7 * snh1_1341[k]
                    + f_3 * pc_y[k] * sni_1788[k];

        t_2300[k] = f_14 * smi_1509[k]
                    + f_8 * snh0_1342[k]
                    - f_9 * snh1_1342[k]
                    + f_3 * pc_y[k] * sni_1789[k];

        t_2301[k] = f_14 * smi_1510[k]
                    + f_10 * snh0_1343[k]
                    - f_11 * snh1_1343[k]
                    + f_3 * pc_y[k] * sni_1790[k];
    }

#pragma omp simd aligned(t_2302, t_2303, t_2304, t_2305, pb_y, pc_y, pc_z, smk0_1944, \
                         smi_1483, smi_1511, smi_1512, smk1_1944, snh0_1343, snh1_1343, \
                         sni_1791, sni_1792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2302[k] = f_14 * smi_1511[k]
                    + f_3 * pc_y[k] * sni_1791[k];

        t_2303[k] = f_19 * smi_1483[k]
                    + f_1 * snh0_1343[k]
                    - f_2 * snh1_1343[k]
                    + f_3 * pc_z[k] * sni_1791[k];

        t_2304[k] = pb_y[k] * smk0_1944[k]
                    - f_12 * pc_y[k] * smk1_1944[k];

        t_2305[k] = f_13 * smi_1512[k]
                    + f_3 * pc_y[k] * sni_1792[k];
    }

#pragma omp simd aligned(t_2306, t_2307, t_2308, pc_x, pc_y, pc_z, smi_1484, smi_1514, \
                         snh0_1347, snh1_1347, sni_1792, sni_1794, \
                         sni_1795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2306[k] = f_18 * smi_1484[k]
                    + f_3 * pc_z[k] * sni_1792[k];

        t_2307[k] = f_4 * snh0_1347[k]
                    - f_5 * snh1_1347[k]
                    + f_3 * pc_x[k] * sni_1795[k];

        t_2308[k] = f_13 * smi_1514[k]
                    + f_3 * pc_y[k] * sni_1794[k];
    }

#pragma omp simd aligned(t_2309, t_2310, t_2311, pb_y, pc_x, pc_y, pc_z, smk0_1949, smi_1487, \
                         smk1_1949, snh0_1350, snh1_1350, sni_1795, \
                         sni_1798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2309[k] = pb_y[k] * smk0_1949[k]
                    - f_12 * pc_y[k] * smk1_1949[k];

        t_2310[k] = f_6 * snh0_1350[k]
                    - f_7 * snh1_1350[k]
                    + f_3 * pc_x[k] * sni_1798[k];

        t_2311[k] = f_18 * smi_1487[k]
                    + f_3 * pc_z[k] * sni_1795[k];
    }

#pragma omp simd aligned(t_2312, t_2313, t_2314, pb_y, pc_x, pc_y, smk0_1953, smi_1517, \
                         smk1_1953, snh0_1354, snh1_1354, sni_1797, \
                         sni_1802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2312[k] = f_13 * smi_1517[k]
                    + f_3 * pc_y[k] * sni_1797[k];

        t_2313[k] = pb_y[k] * smk0_1953[k]
                    - f_12 * pc_y[k] * smk1_1953[k];

        t_2314[k] = f_8 * snh0_1354[k]
                    - f_9 * snh1_1354[k]
                    + f_3 * pc_x[k] * sni_1802[k];
    }

#pragma omp simd aligned(t_2315, t_2316, t_2317, pc_x, pc_y, pc_z, smi_1490, smi_1521, \
                         snh0_1356, snh1_1356, sni_1798, sni_1801, \
                         sni_1804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2315[k] = f_18 * smi_1490[k]
                    + f_3 * pc_z[k] * sni_1798[k];

        t_2316[k] = f_8 * snh0_1356[k]
                    - f_9 * snh1_1356[k]
                    + f_3 * pc_x[k] * sni_1804[k];

        t_2317[k] = f_13 * smi_1521[k]
                    + f_3 * pc_y[k] * sni_1801[k];
    }

#pragma omp simd aligned(t_2318, t_2319, t_2320, pb_y, pc_x, pc_y, pc_z, smk0_1958, smi_1494, \
                         smk1_1958, snh0_1359, snh1_1359, sni_1802, \
                         sni_1807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2318[k] = pb_y[k] * smk0_1958[k]
                    - f_12 * pc_y[k] * smk1_1958[k];

        t_2319[k] = f_10 * snh0_1359[k]
                    - f_11 * snh1_1359[k]
                    + f_3 * pc_x[k] * sni_1807[k];

        t_2320[k] = f_18 * smi_1494[k]
                    + f_3 * pc_z[k] * sni_1802[k];
    }

#pragma omp simd aligned(t_2321, t_2322, t_2323, pc_x, pc_y, smi_1526, snh0_1361, snh0_1362, \
                         snh1_1361, snh1_1362, sni_1806, sni_1809, \
                         sni_1810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2321[k] = f_10 * snh0_1361[k]
                    - f_11 * snh1_1361[k]
                    + f_3 * pc_x[k] * sni_1809[k];

        t_2322[k] = f_10 * snh0_1362[k]
                    - f_11 * snh1_1362[k]
                    + f_3 * pc_x[k] * sni_1810[k];

        t_2323[k] = f_13 * smi_1526[k]
                    + f_3 * pc_y[k] * sni_1806[k];
    }

#pragma omp simd aligned(t_2324, t_2325, t_2326, t_2327, t_2328, t_2329, pb_y, pc_x, pc_y, \
                         smk0_1964, smk1_1964, sni_1813, sni_1814, sni_1815, sni_1816, \
                         sni_1817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2324[k] = pb_y[k] * smk0_1964[k]
                    - f_12 * pc_y[k] * smk1_1964[k];

        t_2325[k] = f_3 * pc_x[k] * sni_1813[k];

        t_2326[k] = f_3 * pc_x[k] * sni_1814[k];

        t_2327[k] = f_3 * pc_x[k] * sni_1815[k];

        t_2328[k] = f_3 * pc_x[k] * sni_1816[k];

        t_2329[k] = f_3 * pc_x[k] * sni_1817[k];
    }

#pragma omp simd aligned(t_2330, t_2331, t_2332, t_2333, pb_y, pc_x, pc_y, pc_z, smk0_1972, \
                         smi_1505, smi_1533, smk1_1972, sni_1813, sni_1818, \
                         sni_1819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2330[k] = f_3 * pc_x[k] * sni_1818[k];

        t_2331[k] = f_3 * pc_x[k] * sni_1819[k];

        t_2332[k] = pb_y[k] * smk0_1972[k]
                    + f_20 * smi_1533[k]
                    - f_12 * pc_y[k] * smk1_1972[k];

        t_2333[k] = f_18 * smi_1505[k]
                    + f_3 * pc_z[k] * sni_1813[k];
    }

#pragma omp simd aligned(t_2334, t_2335, t_2336, pb_y, pc_y, smk0_1974, smk0_1975, smk0_1976, \
                         smi_1535, smi_1536, smi_1537, smk1_1974, smk1_1975, \
                         smk1_1976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2334[k] = pb_y[k] * smk0_1974[k]
                    + f_17 * smi_1535[k]
                    - f_12 * pc_y[k] * smk1_1974[k];

        t_2335[k] = pb_y[k] * smk0_1975[k]
                    + f_16 * smi_1536[k]
                    - f_12 * pc_y[k] * smk1_1975[k];

        t_2336[k] = pb_y[k] * smk0_1976[k]
                    + f_15 * smi_1537[k]
                    - f_12 * pc_y[k] * smk1_1976[k];
    }

#pragma omp simd aligned(t_2337, t_2338, t_2339, pb_y, pc_y, smk0_1977, smk0_1979, smi_1538, \
                         smi_1539, smk1_1977, smk1_1979, sni_1819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2337[k] = pb_y[k] * smk0_1977[k]
                    + f_14 * smi_1538[k]
                    - f_12 * pc_y[k] * smk1_1977[k];

        t_2338[k] = f_13 * smi_1539[k]
                    + f_3 * pc_y[k] * sni_1819[k];

        t_2339[k] = pb_y[k] * smk0_1979[k]
                    - f_12 * pc_y[k] * smk1_1979[k];
    }

#pragma omp simd aligned(t_2340, t_2341, t_2342, t_2343, t_2344, pc_x, pc_y, pc_z, smi_1512, \
                         snh0_1365, snh0_1368, snh1_1365, snh1_1368, sni_1820, sni_1822, \
                         sni_1823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2340[k] = f_1 * snh0_1365[k]
                    - f_2 * snh1_1365[k]
                    + f_3 * pc_x[k] * sni_1820[k];

        t_2341[k] = f_3 * pc_y[k] * sni_1820[k];

        t_2342[k] = f_0 * smi_1512[k]
                    + f_3 * pc_z[k] * sni_1820[k];

        t_2343[k] = f_4 * snh0_1368[k]
                    - f_5 * snh1_1368[k]
                    + f_3 * pc_x[k] * sni_1823[k];

        t_2344[k] = f_3 * pc_y[k] * sni_1822[k];
    }

#pragma omp simd aligned(t_2345, t_2346, t_2347, t_2348, pc_x, pc_y, pc_z, smi_1515, \
                         snh0_1370, snh0_1371, snh1_1370, snh1_1371, sni_1823, sni_1825, \
                         sni_1826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2345[k] = f_4 * snh0_1370[k]
                    - f_5 * snh1_1370[k]
                    + f_3 * pc_x[k] * sni_1825[k];

        t_2346[k] = f_6 * snh0_1371[k]
                    - f_7 * snh1_1371[k]
                    + f_3 * pc_x[k] * sni_1826[k];

        t_2347[k] = f_0 * smi_1515[k]
                    + f_3 * pc_z[k] * sni_1823[k];

        t_2348[k] = f_3 * pc_y[k] * sni_1825[k];
    }

#pragma omp simd aligned(t_2349, t_2350, t_2351, pc_x, pc_z, smi_1518, snh0_1374, snh0_1375, \
                         snh1_1374, snh1_1375, sni_1826, sni_1829, \
                         sni_1830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2349[k] = f_6 * snh0_1374[k]
                    - f_7 * snh1_1374[k]
                    + f_3 * pc_x[k] * sni_1829[k];

        t_2350[k] = f_8 * snh0_1375[k]
                    - f_9 * snh1_1375[k]
                    + f_3 * pc_x[k] * sni_1830[k];

        t_2351[k] = f_0 * smi_1518[k]
                    + f_3 * pc_z[k] * sni_1826[k];
    }

#pragma omp simd aligned(t_2352, t_2353, t_2354, t_2355, pc_x, pc_y, snh0_1377, snh0_1379, \
                         snh0_1380, snh1_1377, snh1_1379, snh1_1380, sni_1829, sni_1832, \
                         sni_1834, sni_1835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2352[k] = f_8 * snh0_1377[k]
                    - f_9 * snh1_1377[k]
                    + f_3 * pc_x[k] * sni_1832[k];

        t_2353[k] = f_3 * pc_y[k] * sni_1829[k];

        t_2354[k] = f_8 * snh0_1379[k]
                    - f_9 * snh1_1379[k]
                    + f_3 * pc_x[k] * sni_1834[k];

        t_2355[k] = f_10 * snh0_1380[k]
                    - f_11 * snh1_1380[k]
                    + f_3 * pc_x[k] * sni_1835[k];
    }

#pragma omp simd aligned(t_2356, t_2357, t_2358, t_2359, pc_x, pc_y, pc_z, smi_1522, \
                         snh0_1382, snh0_1383, snh1_1382, snh1_1383, sni_1830, sni_1834, \
                         sni_1837, sni_1838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2356[k] = f_0 * smi_1522[k]
                    + f_3 * pc_z[k] * sni_1830[k];

        t_2357[k] = f_10 * snh0_1382[k]
                    - f_11 * snh1_1382[k]
                    + f_3 * pc_x[k] * sni_1837[k];

        t_2358[k] = f_10 * snh0_1383[k]
                    - f_11 * snh1_1383[k]
                    + f_3 * pc_x[k] * sni_1838[k];

        t_2359[k] = f_3 * pc_y[k] * sni_1834[k];
    }

#pragma omp simd aligned(t_2360, t_2361, t_2362, t_2363, t_2364, t_2365, pc_x, snh0_1385, \
                         snh1_1385, sni_1840, sni_1841, sni_1842, sni_1843, sni_1844, \
                         sni_1845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2360[k] = f_10 * snh0_1385[k]
                    - f_11 * snh1_1385[k]
                    + f_3 * pc_x[k] * sni_1840[k];

        t_2361[k] = f_3 * pc_x[k] * sni_1841[k];

        t_2362[k] = f_3 * pc_x[k] * sni_1842[k];

        t_2363[k] = f_3 * pc_x[k] * sni_1843[k];

        t_2364[k] = f_3 * pc_x[k] * sni_1844[k];

        t_2365[k] = f_3 * pc_x[k] * sni_1845[k];
    }

#pragma omp simd aligned(t_2366, t_2367, t_2368, t_2369, pc_x, pc_y, pc_z, smi_1533, \
                         snh0_1380, snh1_1380, sni_1841, sni_1846, \
                         sni_1847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2366[k] = f_3 * pc_x[k] * sni_1846[k];

        t_2367[k] = f_3 * pc_x[k] * sni_1847[k];

        t_2368[k] = f_1 * snh0_1380[k]
                    - f_2 * snh1_1380[k]
                    + f_3 * pc_y[k] * sni_1841[k];

        t_2369[k] = f_0 * smi_1533[k]
                    + f_3 * pc_z[k] * sni_1841[k];
    }

#pragma omp simd aligned(t_2370, t_2371, t_2372, pc_y, snh0_1382, snh0_1383, snh0_1384, \
                         snh1_1382, snh1_1383, snh1_1384, sni_1843, sni_1844, \
                         sni_1845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2370[k] = f_4 * snh0_1382[k]
                    - f_5 * snh1_1382[k]
                    + f_3 * pc_y[k] * sni_1843[k];

        t_2371[k] = f_6 * snh0_1383[k]
                    - f_7 * snh1_1383[k]
                    + f_3 * pc_y[k] * sni_1844[k];

        t_2372[k] = f_8 * snh0_1384[k]
                    - f_9 * snh1_1384[k]
                    + f_3 * pc_y[k] * sni_1845[k];
    }

#pragma omp simd aligned(t_2373, t_2374, t_2375, pc_y, pc_z, smi_1539, snh0_1385, snh1_1385, \
                         sni_1846, sni_1847 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2373[k] = f_10 * snh0_1385[k]
                    - f_11 * snh1_1385[k]
                    + f_3 * pc_y[k] * sni_1846[k];

        t_2374[k] = f_3 * pc_y[k] * sni_1847[k];

        t_2375[k] = f_0 * smi_1539[k]
                    + f_1 * snh0_1385[k]
                    - f_2 * snh1_1385[k]
                    + f_3 * pc_z[k] * sni_1847[k];
    }
}

auto
compute_prim_snk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t smk0, const size_t smi,
                                                   const size_t smk1, const size_t snh0,
                                                   const size_t snh1, const size_t sni,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_snk_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, smk0, smi,
                                                              smk1, snh0, snh1, sni, ncols,
                                                              gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece10(buffer, target, pc, smi, snh0,
                                                               snh1, sni, ncols, gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, smk0,
                                                               smi, smk1, snh0, snh1, sni,
                                                               ncols, gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece12(buffer, target, pb, pc, smk0,
                                                               smi, smk1, snh0, snh1, sni,
                                                               ncols, gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece13(buffer, target, pc, smi, snh0,
                                                               snh1, sni, ncols, gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece14(buffer, target, pb, pc, smk0,
                                                               smi, smk1, snh0, snh1, sni,
                                                               ncols, gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece15(buffer, target, pb, pc, smk0,
                                                               smi, smk1, sni, ncols, gamma, p,
                                                               q);

    compute_prim_snk_three_center_electron_repulsion_0_piece16(buffer, target, pb, pc, smk0,
                                                               smi, smk1, sni, ncols, gamma, p,
                                                               q);

    compute_prim_snk_three_center_electron_repulsion_0_piece17(buffer, target, pb, pc, smk0,
                                                               smi, smk1, snh0, snh1, sni,
                                                               ncols, gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece18(buffer, target, pb, pc, smk0,
                                                               smi, smk1, snh0, snh1, sni,
                                                               ncols, gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece19(buffer, target, pc, smi, snh0,
                                                               snh1, sni, ncols, gamma, p, q);

    compute_prim_snk_three_center_electron_repulsion_0_piece20(buffer, target, pb, pc, smk0,
                                                               smi, smk1, snh0, snh1, sni,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
