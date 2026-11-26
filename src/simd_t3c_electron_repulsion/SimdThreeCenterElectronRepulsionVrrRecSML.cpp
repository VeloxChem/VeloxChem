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


#include "SimdThreeCenterElectronRepulsionVrrRecSML.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_0 = buffer.data(sll0 + 0);
    const auto *sll0_3 = buffer.data(sll0 + 3);
    const auto *sll0_5 = buffer.data(sll0 + 5);
    const auto *sll0_6 = buffer.data(sll0 + 6);
    const auto *sll0_9 = buffer.data(sll0 + 9);
    const auto *sll0_10 = buffer.data(sll0 + 10);
    const auto *sll0_12 = buffer.data(sll0 + 12);
    const auto *sll0_14 = buffer.data(sll0 + 14);
    const auto *sll0_15 = buffer.data(sll0 + 15);
    const auto *sll0_17 = buffer.data(sll0 + 17);
    const auto *sll0_18 = buffer.data(sll0 + 18);
    const auto *sll0_20 = buffer.data(sll0 + 20);
    const auto *sll0_21 = buffer.data(sll0 + 21);
    const auto *sll0_23 = buffer.data(sll0 + 23);
    const auto *sll0_24 = buffer.data(sll0 + 24);
    const auto *sll0_25 = buffer.data(sll0 + 25);
    const auto *sll0_27 = buffer.data(sll0 + 27);
    const auto *sll0_44 = buffer.data(sll0 + 44);

    const auto *slk_0 = buffer.data(slk + 0);
    const auto *slk_1 = buffer.data(slk + 1);
    const auto *slk_2 = buffer.data(slk + 2);
    const auto *slk_3 = buffer.data(slk + 3);
    const auto *slk_5 = buffer.data(slk + 5);
    const auto *slk_6 = buffer.data(slk + 6);
    const auto *slk_7 = buffer.data(slk + 7);
    const auto *slk_8 = buffer.data(slk + 8);
    const auto *slk_9 = buffer.data(slk + 9);
    const auto *slk_10 = buffer.data(slk + 10);
    const auto *slk_11 = buffer.data(slk + 11);
    const auto *slk_12 = buffer.data(slk + 12);
    const auto *slk_13 = buffer.data(slk + 13);
    const auto *slk_14 = buffer.data(slk + 14);
    const auto *slk_15 = buffer.data(slk + 15);
    const auto *slk_16 = buffer.data(slk + 16);
    const auto *slk_17 = buffer.data(slk + 17);
    const auto *slk_18 = buffer.data(slk + 18);
    const auto *slk_19 = buffer.data(slk + 19);
    const auto *slk_20 = buffer.data(slk + 20);
    const auto *slk_21 = buffer.data(slk + 21);
    const auto *slk_23 = buffer.data(slk + 23);
    const auto *slk_24 = buffer.data(slk + 24);
    const auto *slk_25 = buffer.data(slk + 25);
    const auto *slk_27 = buffer.data(slk + 27);
    const auto *slk_28 = buffer.data(slk + 28);
    const auto *slk_29 = buffer.data(slk + 29);
    const auto *slk_30 = buffer.data(slk + 30);
    const auto *slk_31 = buffer.data(slk + 31);
    const auto *slk_32 = buffer.data(slk + 32);
    const auto *slk_33 = buffer.data(slk + 33);
    const auto *slk_34 = buffer.data(slk + 34);
    const auto *slk_35 = buffer.data(slk + 35);
    const auto *slk_64 = buffer.data(slk + 64);
    const auto *slk_65 = buffer.data(slk + 65);
    const auto *slk_66 = buffer.data(slk + 66);
    const auto *slk_67 = buffer.data(slk + 67);
    const auto *slk_68 = buffer.data(slk + 68);
    const auto *slk_69 = buffer.data(slk + 69);
    const auto *slk_70 = buffer.data(slk + 70);
    const auto *slk_71 = buffer.data(slk + 71);

    const auto *sll1_0 = buffer.data(sll1 + 0);
    const auto *sll1_3 = buffer.data(sll1 + 3);
    const auto *sll1_5 = buffer.data(sll1 + 5);
    const auto *sll1_6 = buffer.data(sll1 + 6);
    const auto *sll1_9 = buffer.data(sll1 + 9);
    const auto *sll1_10 = buffer.data(sll1 + 10);
    const auto *sll1_12 = buffer.data(sll1 + 12);
    const auto *sll1_14 = buffer.data(sll1 + 14);
    const auto *sll1_15 = buffer.data(sll1 + 15);
    const auto *sll1_17 = buffer.data(sll1 + 17);
    const auto *sll1_18 = buffer.data(sll1 + 18);
    const auto *sll1_20 = buffer.data(sll1 + 20);
    const auto *sll1_21 = buffer.data(sll1 + 21);
    const auto *sll1_23 = buffer.data(sll1 + 23);
    const auto *sll1_24 = buffer.data(sll1 + 24);
    const auto *sll1_25 = buffer.data(sll1 + 25);
    const auto *sll1_27 = buffer.data(sll1 + 27);
    const auto *sll1_44 = buffer.data(sll1 + 44);

    const auto *smi0_0 = buffer.data(smi0 + 0);
    const auto *smi0_3 = buffer.data(smi0 + 3);
    const auto *smi0_5 = buffer.data(smi0 + 5);
    const auto *smi0_6 = buffer.data(smi0 + 6);
    const auto *smi0_9 = buffer.data(smi0 + 9);
    const auto *smi0_10 = buffer.data(smi0 + 10);
    const auto *smi0_12 = buffer.data(smi0 + 12);
    const auto *smi0_14 = buffer.data(smi0 + 14);
    const auto *smi0_15 = buffer.data(smi0 + 15);
    const auto *smi0_17 = buffer.data(smi0 + 17);
    const auto *smi0_18 = buffer.data(smi0 + 18);
    const auto *smi0_20 = buffer.data(smi0 + 20);
    const auto *smi0_21 = buffer.data(smi0 + 21);
    const auto *smi0_23 = buffer.data(smi0 + 23);
    const auto *smi0_24 = buffer.data(smi0 + 24);
    const auto *smi0_25 = buffer.data(smi0 + 25);
    const auto *smi0_26 = buffer.data(smi0 + 26);
    const auto *smi0_27 = buffer.data(smi0 + 27);
    const auto *smi0_49 = buffer.data(smi0 + 49);
    const auto *smi0_51 = buffer.data(smi0 + 51);
    const auto *smi0_52 = buffer.data(smi0 + 52);
    const auto *smi0_53 = buffer.data(smi0 + 53);
    const auto *smi0_54 = buffer.data(smi0 + 54);
    const auto *smi0_55 = buffer.data(smi0 + 55);

    const auto *smi1_0 = buffer.data(smi1 + 0);
    const auto *smi1_3 = buffer.data(smi1 + 3);
    const auto *smi1_5 = buffer.data(smi1 + 5);
    const auto *smi1_6 = buffer.data(smi1 + 6);
    const auto *smi1_9 = buffer.data(smi1 + 9);
    const auto *smi1_10 = buffer.data(smi1 + 10);
    const auto *smi1_12 = buffer.data(smi1 + 12);
    const auto *smi1_14 = buffer.data(smi1 + 14);
    const auto *smi1_15 = buffer.data(smi1 + 15);
    const auto *smi1_17 = buffer.data(smi1 + 17);
    const auto *smi1_18 = buffer.data(smi1 + 18);
    const auto *smi1_20 = buffer.data(smi1 + 20);
    const auto *smi1_21 = buffer.data(smi1 + 21);
    const auto *smi1_23 = buffer.data(smi1 + 23);
    const auto *smi1_24 = buffer.data(smi1 + 24);
    const auto *smi1_25 = buffer.data(smi1 + 25);
    const auto *smi1_26 = buffer.data(smi1 + 26);
    const auto *smi1_27 = buffer.data(smi1 + 27);
    const auto *smi1_49 = buffer.data(smi1 + 49);
    const auto *smi1_51 = buffer.data(smi1 + 51);
    const auto *smi1_52 = buffer.data(smi1 + 52);
    const auto *smi1_53 = buffer.data(smi1 + 53);
    const auto *smi1_54 = buffer.data(smi1 + 54);
    const auto *smi1_55 = buffer.data(smi1 + 55);

    const auto *smk_0 = buffer.data(smk + 0);
    const auto *smk_2 = buffer.data(smk + 2);
    const auto *smk_3 = buffer.data(smk + 3);
    const auto *smk_5 = buffer.data(smk + 5);
    const auto *smk_6 = buffer.data(smk + 6);
    const auto *smk_9 = buffer.data(smk + 9);
    const auto *smk_10 = buffer.data(smk + 10);
    const auto *smk_12 = buffer.data(smk + 12);
    const auto *smk_14 = buffer.data(smk + 14);
    const auto *smk_15 = buffer.data(smk + 15);
    const auto *smk_17 = buffer.data(smk + 17);
    const auto *smk_18 = buffer.data(smk + 18);
    const auto *smk_20 = buffer.data(smk + 20);
    const auto *smk_21 = buffer.data(smk + 21);
    const auto *smk_23 = buffer.data(smk + 23);
    const auto *smk_24 = buffer.data(smk + 24);
    const auto *smk_25 = buffer.data(smk + 25);
    const auto *smk_27 = buffer.data(smk + 27);
    const auto *smk_28 = buffer.data(smk + 28);
    const auto *smk_29 = buffer.data(smk + 29);
    const auto *smk_30 = buffer.data(smk + 30);
    const auto *smk_31 = buffer.data(smk + 31);
    const auto *smk_32 = buffer.data(smk + 32);
    const auto *smk_33 = buffer.data(smk + 33);
    const auto *smk_34 = buffer.data(smk + 34);
    const auto *smk_35 = buffer.data(smk + 35);
    const auto *smk_36 = buffer.data(smk + 36);
    const auto *smk_38 = buffer.data(smk + 38);
    const auto *smk_39 = buffer.data(smk + 39);
    const auto *smk_41 = buffer.data(smk + 41);
    const auto *smk_42 = buffer.data(smk + 42);
    const auto *smk_45 = buffer.data(smk + 45);
    const auto *smk_46 = buffer.data(smk + 46);
    const auto *smk_50 = buffer.data(smk + 50);
    const auto *smk_51 = buffer.data(smk + 51);
    const auto *smk_56 = buffer.data(smk + 56);
    const auto *smk_64 = buffer.data(smk + 64);
    const auto *smk_65 = buffer.data(smk + 65);
    const auto *smk_66 = buffer.data(smk + 66);
    const auto *smk_67 = buffer.data(smk + 67);
    const auto *smk_68 = buffer.data(smk + 68);
    const auto *smk_69 = buffer.data(smk + 69);
    const auto *smk_70 = buffer.data(smk + 70);
    const auto *smk_71 = buffer.data(smk + 71);
    const auto *smk_72 = buffer.data(smk + 72);
    const auto *smk_74 = buffer.data(smk + 74);
    const auto *smk_75 = buffer.data(smk + 75);
    const auto *smk_77 = buffer.data(smk + 77);
    const auto *smk_78 = buffer.data(smk + 78);
    const auto *smk_81 = buffer.data(smk + 81);
    const auto *smk_82 = buffer.data(smk + 82);
    const auto *smk_86 = buffer.data(smk + 86);
    const auto *smk_87 = buffer.data(smk + 87);
    const auto *smk_92 = buffer.data(smk + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, slk_0, slk_3, smi0_0, smi0_3, \
                         smi1_0, smi1_3, smk_0, smk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * slk_0[k]
                 + f_1 * smi0_0[k]
                 - f_2 * smi1_0[k]
                 + f_3 * pc_x[k] * smk_0[k];

        t_1[k] = f_3 * pc_y[k] * smk_0[k];

        t_2[k] = f_3 * pc_z[k] * smk_0[k];

        t_3[k] = f_0 * slk_3[k]
                 + f_4 * smi0_3[k]
                 - f_5 * smi1_3[k]
                 + f_3 * pc_x[k] * smk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, slk_5, slk_6, smi0_5, smi0_6, smi1_5, \
                         smi1_6, smk_2, smk_5, smk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * smk_2[k];

        t_5[k] = f_0 * slk_5[k]
                 + f_4 * smi0_5[k]
                 - f_5 * smi1_5[k]
                 + f_3 * pc_x[k] * smk_5[k];

        t_6[k] = f_0 * slk_6[k]
                 + f_6 * smi0_6[k]
                 - f_7 * smi1_6[k]
                 + f_3 * pc_x[k] * smk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, slk_9, smi0_9, smi1_9, smk_3, smk_5, \
                         smk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * smk_3[k];

        t_8[k] = f_3 * pc_y[k] * smk_5[k];

        t_9[k] = f_0 * slk_9[k]
                 + f_6 * smi0_9[k]
                 - f_7 * smi1_9[k]
                 + f_3 * pc_x[k] * smk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, slk_10, slk_12, smi0_10, smi0_12, \
                         smi1_10, smi1_12, smk_6, smk_10, smk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * slk_10[k]
                  + f_8 * smi0_10[k]
                  - f_9 * smi1_10[k]
                  + f_3 * pc_x[k] * smk_10[k];

        t_11[k] = f_3 * pc_z[k] * smk_6[k];

        t_12[k] = f_0 * slk_12[k]
                  + f_8 * smi0_12[k]
                  - f_9 * smi1_12[k]
                  + f_3 * pc_x[k] * smk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, slk_14, slk_15, smi0_14, smi0_15, \
                         smi1_14, smi1_15, smk_9, smk_14, smk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * smk_9[k];

        t_14[k] = f_0 * slk_14[k]
                  + f_8 * smi0_14[k]
                  - f_9 * smi1_14[k]
                  + f_3 * pc_x[k] * smk_14[k];

        t_15[k] = f_0 * slk_15[k]
                  + f_10 * smi0_15[k]
                  - f_11 * smi1_15[k]
                  + f_3 * pc_x[k] * smk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, slk_17, slk_18, smi0_17, smi0_18, \
                         smi1_17, smi1_18, smk_10, smk_17, smk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * smk_10[k];

        t_17[k] = f_0 * slk_17[k]
                  + f_10 * smi0_17[k]
                  - f_11 * smi1_17[k]
                  + f_3 * pc_x[k] * smk_17[k];

        t_18[k] = f_0 * slk_18[k]
                  + f_10 * smi0_18[k]
                  - f_11 * smi1_18[k]
                  + f_3 * pc_x[k] * smk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pc_x, pc_y, slk_20, slk_21, smi0_20, smi0_21, \
                         smi1_20, smi1_21, smk_14, smk_20, smk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * smk_14[k];

        t_20[k] = f_0 * slk_20[k]
                  + f_10 * smi0_20[k]
                  - f_11 * smi1_20[k]
                  + f_3 * pc_x[k] * smk_20[k];

        t_21[k] = f_0 * slk_21[k]
                  + f_12 * smi0_21[k]
                  - f_13 * smi1_21[k]
                  + f_3 * pc_x[k] * smk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pc_x, pc_z, slk_23, slk_24, smi0_23, smi0_24, \
                         smi1_23, smi1_24, smk_15, smk_23, smk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * smk_15[k];

        t_23[k] = f_0 * slk_23[k]
                  + f_12 * smi0_23[k]
                  - f_13 * smi1_23[k]
                  + f_3 * pc_x[k] * smk_23[k];

        t_24[k] = f_0 * slk_24[k]
                  + f_12 * smi0_24[k]
                  - f_13 * smi1_24[k]
                  + f_3 * pc_x[k] * smk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pc_x, pc_y, slk_25, slk_27, smi0_25, smi0_27, \
                         smi1_25, smi1_27, smk_20, smk_25, smk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * slk_25[k]
                  + f_12 * smi0_25[k]
                  - f_13 * smi1_25[k]
                  + f_3 * pc_x[k] * smk_25[k];

        t_26[k] = f_3 * pc_y[k] * smk_20[k];

        t_27[k] = f_0 * slk_27[k]
                  + f_12 * smi0_27[k]
                  - f_13 * smi1_27[k]
                  + f_3 * pc_x[k] * smk_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pc_x, slk_28, slk_29, slk_30, slk_31, \
                         slk_32, smk_28, smk_29, smk_30, smk_31, \
                         smk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * slk_28[k]
                  + f_3 * pc_x[k] * smk_28[k];

        t_29[k] = f_0 * slk_29[k]
                  + f_3 * pc_x[k] * smk_29[k];

        t_30[k] = f_0 * slk_30[k]
                  + f_3 * pc_x[k] * smk_30[k];

        t_31[k] = f_0 * slk_31[k]
                  + f_3 * pc_x[k] * smk_31[k];

        t_32[k] = f_0 * slk_32[k]
                  + f_3 * pc_x[k] * smk_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pc_x, pc_y, slk_33, slk_34, slk_35, smi0_21, \
                         smi1_21, smk_28, smk_33, smk_34, smk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * slk_33[k]
                  + f_3 * pc_x[k] * smk_33[k];

        t_34[k] = f_0 * slk_34[k]
                  + f_3 * pc_x[k] * smk_34[k];

        t_35[k] = f_0 * slk_35[k]
                  + f_3 * pc_x[k] * smk_35[k];

        t_36[k] = f_1 * smi0_21[k]
                  - f_2 * smi1_21[k]
                  + f_3 * pc_y[k] * smk_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pc_y, pc_z, smi0_23, smi0_24, smi0_25, \
                         smi1_23, smi1_24, smi1_25, smk_28, smk_30, smk_31, \
                         smk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * smk_28[k];

        t_38[k] = f_4 * smi0_23[k]
                  - f_5 * smi1_23[k]
                  + f_3 * pc_y[k] * smk_30[k];

        t_39[k] = f_6 * smi0_24[k]
                  - f_7 * smi1_24[k]
                  + f_3 * pc_y[k] * smk_31[k];

        t_40[k] = f_8 * smi0_25[k]
                  - f_9 * smi1_25[k]
                  + f_3 * pc_y[k] * smk_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, smi0_26, smi0_27, smi1_26, \
                         smi1_27, smk_33, smk_34, smk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * smi0_26[k]
                  - f_11 * smi1_26[k]
                  + f_3 * pc_y[k] * smk_33[k];

        t_42[k] = f_12 * smi0_27[k]
                  - f_13 * smi1_27[k]
                  + f_3 * pc_y[k] * smk_34[k];

        t_43[k] = f_3 * pc_y[k] * smk_35[k];

        t_44[k] = f_1 * smi0_27[k]
                  - f_2 * smi1_27[k]
                  + f_3 * pc_z[k] * smk_35[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_y, pc_y, pc_z, sll0_0, sll0_3, slk_0, \
                         slk_1, sll1_0, sll1_3, smk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * sll0_0[k]
                  - f_14 * pc_y[k] * sll1_0[k];

        t_46[k] = f_15 * slk_0[k]
                  + f_3 * pc_y[k] * smk_36[k];

        t_47[k] = f_3 * pc_z[k] * smk_36[k];

        t_48[k] = pb_y[k] * sll0_3[k]
                  + f_16 * slk_1[k]
                  - f_14 * pc_y[k] * sll1_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_y, pc_y, pc_z, sll0_5, sll0_6, slk_2, \
                         slk_3, sll1_5, sll1_6, smk_38, smk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_15 * slk_2[k]
                  + f_3 * pc_y[k] * smk_38[k];

        t_50[k] = pb_y[k] * sll0_5[k]
                  - f_14 * pc_y[k] * sll1_5[k];

        t_51[k] = pb_y[k] * sll0_6[k]
                  + f_17 * slk_3[k]
                  - f_14 * pc_y[k] * sll1_6[k];

        t_52[k] = f_3 * pc_z[k] * smk_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_y, pc_y, pc_z, sll0_9, sll0_10, slk_5, \
                         slk_6, sll1_9, sll1_10, smk_41, smk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_15 * slk_5[k]
                  + f_3 * pc_y[k] * smk_41[k];

        t_54[k] = pb_y[k] * sll0_9[k]
                  - f_14 * pc_y[k] * sll1_9[k];

        t_55[k] = pb_y[k] * sll0_10[k]
                  + f_18 * slk_6[k]
                  - f_14 * pc_y[k] * sll1_10[k];

        t_56[k] = f_3 * pc_z[k] * smk_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, sll0_12, sll0_14, sll0_15, slk_8, \
                         slk_9, slk_10, sll1_12, sll1_14, sll1_15, \
                         smk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * sll0_12[k]
                  + f_16 * slk_8[k]
                  - f_14 * pc_y[k] * sll1_12[k];

        t_58[k] = f_15 * slk_9[k]
                  + f_3 * pc_y[k] * smk_45[k];

        t_59[k] = pb_y[k] * sll0_14[k]
                  - f_14 * pc_y[k] * sll1_14[k];

        t_60[k] = pb_y[k] * sll0_15[k]
                  + f_19 * slk_10[k]
                  - f_14 * pc_y[k] * sll1_15[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_y, pc_y, pc_z, sll0_17, sll0_18, slk_12, \
                         slk_13, slk_14, sll1_17, sll1_18, smk_46, \
                         smk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pc_z[k] * smk_46[k];

        t_62[k] = pb_y[k] * sll0_17[k]
                  + f_17 * slk_12[k]
                  - f_14 * pc_y[k] * sll1_17[k];

        t_63[k] = pb_y[k] * sll0_18[k]
                  + f_16 * slk_13[k]
                  - f_14 * pc_y[k] * sll1_18[k];

        t_64[k] = f_15 * slk_14[k]
                  + f_3 * pc_y[k] * smk_50[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_y, pc_z, sll0_20, sll0_21, sll0_23, \
                         slk_15, slk_17, sll1_20, sll1_21, sll1_23, \
                         smk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sll0_20[k]
                  - f_14 * pc_y[k] * sll1_20[k];

        t_66[k] = pb_y[k] * sll0_21[k]
                  + f_20 * slk_15[k]
                  - f_14 * pc_y[k] * sll1_21[k];

        t_67[k] = f_3 * pc_z[k] * smk_51[k];

        t_68[k] = pb_y[k] * sll0_23[k]
                  + f_18 * slk_17[k]
                  - f_14 * pc_y[k] * sll1_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_y, pc_y, sll0_24, sll0_25, sll0_27, \
                         slk_18, slk_19, slk_20, sll1_24, sll1_25, sll1_27, \
                         smk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * sll0_24[k]
                  + f_17 * slk_18[k]
                  - f_14 * pc_y[k] * sll1_24[k];

        t_70[k] = pb_y[k] * sll0_25[k]
                  + f_16 * slk_19[k]
                  - f_14 * pc_y[k] * sll1_25[k];

        t_71[k] = f_15 * slk_20[k]
                  + f_3 * pc_y[k] * smk_56[k];

        t_72[k] = pb_y[k] * sll0_27[k]
                  - f_14 * pc_y[k] * sll1_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pc_x, slk_64, slk_65, slk_66, slk_67, \
                         slk_68, smk_64, smk_65, smk_66, smk_67, \
                         smk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_21 * slk_64[k]
                  + f_3 * pc_x[k] * smk_64[k];

        t_74[k] = f_21 * slk_65[k]
                  + f_3 * pc_x[k] * smk_65[k];

        t_75[k] = f_21 * slk_66[k]
                  + f_3 * pc_x[k] * smk_66[k];

        t_76[k] = f_21 * slk_67[k]
                  + f_3 * pc_x[k] * smk_67[k];

        t_77[k] = f_21 * slk_68[k]
                  + f_3 * pc_x[k] * smk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_x, pc_y, slk_28, slk_69, slk_70, slk_71, \
                         smi0_49, smi1_49, smk_64, smk_69, smk_70, \
                         smk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_21 * slk_69[k]
                  + f_3 * pc_x[k] * smk_69[k];

        t_79[k] = f_21 * slk_70[k]
                  + f_3 * pc_x[k] * smk_70[k];

        t_80[k] = f_21 * slk_71[k]
                  + f_3 * pc_x[k] * smk_71[k];

        t_81[k] = f_15 * slk_28[k]
                  + f_1 * smi0_49[k]
                  - f_2 * smi1_49[k]
                  + f_3 * pc_y[k] * smk_64[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pc_y, pc_z, slk_30, slk_31, smi0_51, smi0_52, \
                         smi1_51, smi1_52, smk_64, smk_66, smk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_z[k] * smk_64[k];

        t_83[k] = f_15 * slk_30[k]
                  + f_4 * smi0_51[k]
                  - f_5 * smi1_51[k]
                  + f_3 * pc_y[k] * smk_66[k];

        t_84[k] = f_15 * slk_31[k]
                  + f_6 * smi0_52[k]
                  - f_7 * smi1_52[k]
                  + f_3 * pc_y[k] * smk_67[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_y, slk_32, slk_33, slk_34, smi0_53, smi0_54, \
                         smi0_55, smi1_53, smi1_54, smi1_55, smk_68, smk_69, \
                         smk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_15 * slk_32[k]
                  + f_8 * smi0_53[k]
                  - f_9 * smi1_53[k]
                  + f_3 * pc_y[k] * smk_68[k];

        t_86[k] = f_15 * slk_33[k]
                  + f_10 * smi0_54[k]
                  - f_11 * smi1_54[k]
                  + f_3 * pc_y[k] * smk_69[k];

        t_87[k] = f_15 * slk_34[k]
                  + f_12 * smi0_55[k]
                  - f_13 * smi1_55[k]
                  + f_3 * pc_y[k] * smk_70[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, pc_y, pc_z, sll0_0, sll0_44, \
                         slk_35, sll1_0, sll1_44, smk_71, smk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * slk_35[k]
                  + f_3 * pc_y[k] * smk_71[k];

        t_89[k] = pb_y[k] * sll0_44[k]
                  - f_14 * pc_y[k] * sll1_44[k];

        t_90[k] = pb_z[k] * sll0_0[k]
                  - f_14 * pc_z[k] * sll1_0[k];

        t_91[k] = f_3 * pc_y[k] * smk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_z, pc_y, pc_z, sll0_3, sll0_5, slk_0, \
                         slk_2, sll1_3, sll1_5, smk_72, smk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * slk_0[k]
                  + f_3 * pc_z[k] * smk_72[k];

        t_93[k] = pb_z[k] * sll0_3[k]
                  - f_14 * pc_z[k] * sll1_3[k];

        t_94[k] = f_3 * pc_y[k] * smk_74[k];

        t_95[k] = pb_z[k] * sll0_5[k]
                  + f_16 * slk_2[k]
                  - f_14 * pc_z[k] * sll1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_z, pc_y, pc_z, sll0_6, sll0_9, slk_3, \
                         slk_5, sll1_6, sll1_9, smk_75, smk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_z[k] * sll0_6[k]
                  - f_14 * pc_z[k] * sll1_6[k];

        t_97[k] = f_15 * slk_3[k]
                  + f_3 * pc_z[k] * smk_75[k];

        t_98[k] = f_3 * pc_y[k] * smk_77[k];

        t_99[k] = pb_z[k] * sll0_9[k]
                  + f_17 * slk_5[k]
                  - f_14 * pc_z[k] * sll1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_z, pc_y, pc_z, sll0_10, sll0_12, \
                         slk_6, slk_7, sll1_10, sll1_12, smk_78, \
                         smk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_z[k] * sll0_10[k]
                   - f_14 * pc_z[k] * sll1_10[k];

        t_101[k] = f_15 * slk_6[k]
                   + f_3 * pc_z[k] * smk_78[k];

        t_102[k] = pb_z[k] * sll0_12[k]
                   + f_16 * slk_7[k]
                   - f_14 * pc_z[k] * sll1_12[k];

        t_103[k] = f_3 * pc_y[k] * smk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_z, pc_z, sll0_14, sll0_15, sll0_17, \
                         slk_9, slk_10, slk_11, sll1_14, sll1_15, sll1_17, \
                         smk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * sll0_14[k]
                   + f_18 * slk_9[k]
                   - f_14 * pc_z[k] * sll1_14[k];

        t_105[k] = pb_z[k] * sll0_15[k]
                   - f_14 * pc_z[k] * sll1_15[k];

        t_106[k] = f_15 * slk_10[k]
                   + f_3 * pc_z[k] * smk_82[k];

        t_107[k] = pb_z[k] * sll0_17[k]
                   + f_16 * slk_11[k]
                   - f_14 * pc_z[k] * sll1_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_z, pc_y, pc_z, sll0_18, sll0_20, \
                         sll0_21, slk_12, slk_14, sll1_18, sll1_20, sll1_21, \
                         smk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * sll0_18[k]
                   + f_17 * slk_12[k]
                   - f_14 * pc_z[k] * sll1_18[k];

        t_109[k] = f_3 * pc_y[k] * smk_86[k];

        t_110[k] = pb_z[k] * sll0_20[k]
                   + f_19 * slk_14[k]
                   - f_14 * pc_z[k] * sll1_20[k];

        t_111[k] = pb_z[k] * sll0_21[k]
                   - f_14 * pc_z[k] * sll1_21[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_z, pc_z, sll0_23, sll0_24, slk_15, slk_16, \
                         slk_17, sll1_23, sll1_24, smk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_15 * slk_15[k]
                   + f_3 * pc_z[k] * smk_87[k];

        t_113[k] = pb_z[k] * sll0_23[k]
                   + f_16 * slk_16[k]
                   - f_14 * pc_z[k] * sll1_23[k];

        t_114[k] = pb_z[k] * sll0_24[k]
                   + f_17 * slk_17[k]
                   - f_14 * pc_z[k] * sll1_24[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_z, pc_y, pc_z, sll0_25, sll0_27, slk_18, \
                         slk_20, sll1_25, sll1_27, smk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * sll0_25[k]
                   + f_18 * slk_18[k]
                   - f_14 * pc_z[k] * sll1_25[k];

        t_116[k] = f_3 * pc_y[k] * smk_92[k];

        t_117[k] = pb_z[k] * sll0_27[k]
                   + f_20 * slk_20[k]
                   - f_14 * pc_z[k] * sll1_27[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_36 = buffer.data(sll0 + 36);
    const auto *sll0_48 = buffer.data(sll0 + 48);
    const auto *sll0_51 = buffer.data(sll0 + 51);
    const auto *sll0_55 = buffer.data(sll0 + 55);
    const auto *sll0_60 = buffer.data(sll0 + 60);
    const auto *sll0_66 = buffer.data(sll0 + 66);
    const auto *sll0_81 = buffer.data(sll0 + 81);
    const auto *sll0_90 = buffer.data(sll0 + 90);
    const auto *sll0_95 = buffer.data(sll0 + 95);
    const auto *sll0_99 = buffer.data(sll0 + 99);
    const auto *sll0_102 = buffer.data(sll0 + 102);
    const auto *sll0_104 = buffer.data(sll0 + 104);
    const auto *sll0_107 = buffer.data(sll0 + 107);
    const auto *sll0_108 = buffer.data(sll0 + 108);
    const auto *sll0_110 = buffer.data(sll0 + 110);
    const auto *sll0_113 = buffer.data(sll0 + 113);
    const auto *sll0_114 = buffer.data(sll0 + 114);
    const auto *sll0_115 = buffer.data(sll0 + 115);
    const auto *sll0_117 = buffer.data(sll0 + 117);
    const auto *sll0_134 = buffer.data(sll0 + 134);

    const auto *slk_28 = buffer.data(slk + 28);
    const auto *slk_35 = buffer.data(slk + 35);
    const auto *slk_36 = buffer.data(slk + 36);
    const auto *slk_38 = buffer.data(slk + 38);
    const auto *slk_39 = buffer.data(slk + 39);
    const auto *slk_41 = buffer.data(slk + 41);
    const auto *slk_42 = buffer.data(slk + 42);
    const auto *slk_45 = buffer.data(slk + 45);
    const auto *slk_46 = buffer.data(slk + 46);
    const auto *slk_50 = buffer.data(slk + 50);
    const auto *slk_51 = buffer.data(slk + 51);
    const auto *slk_56 = buffer.data(slk + 56);
    const auto *slk_64 = buffer.data(slk + 64);
    const auto *slk_66 = buffer.data(slk + 66);
    const auto *slk_67 = buffer.data(slk + 67);
    const auto *slk_68 = buffer.data(slk + 68);
    const auto *slk_69 = buffer.data(slk + 69);
    const auto *slk_70 = buffer.data(slk + 70);
    const auto *slk_71 = buffer.data(slk + 71);
    const auto *slk_72 = buffer.data(slk + 72);
    const auto *slk_74 = buffer.data(slk + 74);
    const auto *slk_75 = buffer.data(slk + 75);
    const auto *slk_77 = buffer.data(slk + 77);
    const auto *slk_80 = buffer.data(slk + 80);
    const auto *slk_81 = buffer.data(slk + 81);
    const auto *slk_84 = buffer.data(slk + 84);
    const auto *slk_85 = buffer.data(slk + 85);
    const auto *slk_86 = buffer.data(slk + 86);
    const auto *slk_89 = buffer.data(slk + 89);
    const auto *slk_90 = buffer.data(slk + 90);
    const auto *slk_91 = buffer.data(slk + 91);
    const auto *slk_92 = buffer.data(slk + 92);
    const auto *slk_100 = buffer.data(slk + 100);
    const auto *slk_101 = buffer.data(slk + 101);
    const auto *slk_102 = buffer.data(slk + 102);
    const auto *slk_103 = buffer.data(slk + 103);
    const auto *slk_104 = buffer.data(slk + 104);
    const auto *slk_105 = buffer.data(slk + 105);
    const auto *slk_106 = buffer.data(slk + 106);
    const auto *slk_107 = buffer.data(slk + 107);
    const auto *slk_108 = buffer.data(slk + 108);
    const auto *slk_111 = buffer.data(slk + 111);
    const auto *slk_113 = buffer.data(slk + 113);
    const auto *slk_114 = buffer.data(slk + 114);
    const auto *slk_117 = buffer.data(slk + 117);
    const auto *slk_118 = buffer.data(slk + 118);
    const auto *slk_120 = buffer.data(slk + 120);
    const auto *slk_122 = buffer.data(slk + 122);
    const auto *slk_123 = buffer.data(slk + 123);
    const auto *slk_125 = buffer.data(slk + 125);
    const auto *slk_126 = buffer.data(slk + 126);
    const auto *slk_128 = buffer.data(slk + 128);
    const auto *slk_129 = buffer.data(slk + 129);
    const auto *slk_131 = buffer.data(slk + 131);
    const auto *slk_132 = buffer.data(slk + 132);
    const auto *slk_133 = buffer.data(slk + 133);
    const auto *slk_135 = buffer.data(slk + 135);
    const auto *slk_136 = buffer.data(slk + 136);
    const auto *slk_137 = buffer.data(slk + 137);
    const auto *slk_138 = buffer.data(slk + 138);
    const auto *slk_139 = buffer.data(slk + 139);
    const auto *slk_140 = buffer.data(slk + 140);
    const auto *slk_141 = buffer.data(slk + 141);
    const auto *slk_142 = buffer.data(slk + 142);
    const auto *slk_143 = buffer.data(slk + 143);
    const auto *slk_172 = buffer.data(slk + 172);
    const auto *slk_173 = buffer.data(slk + 173);
    const auto *slk_174 = buffer.data(slk + 174);
    const auto *slk_175 = buffer.data(slk + 175);
    const auto *slk_176 = buffer.data(slk + 176);
    const auto *slk_177 = buffer.data(slk + 177);
    const auto *slk_178 = buffer.data(slk + 178);
    const auto *slk_179 = buffer.data(slk + 179);
    const auto *slk_180 = buffer.data(slk + 180);
    const auto *slk_183 = buffer.data(slk + 183);
    const auto *slk_185 = buffer.data(slk + 185);
    const auto *slk_186 = buffer.data(slk + 186);

    const auto *sll1_36 = buffer.data(sll1 + 36);
    const auto *sll1_48 = buffer.data(sll1 + 48);
    const auto *sll1_51 = buffer.data(sll1 + 51);
    const auto *sll1_55 = buffer.data(sll1 + 55);
    const auto *sll1_60 = buffer.data(sll1 + 60);
    const auto *sll1_66 = buffer.data(sll1 + 66);
    const auto *sll1_81 = buffer.data(sll1 + 81);
    const auto *sll1_90 = buffer.data(sll1 + 90);
    const auto *sll1_95 = buffer.data(sll1 + 95);
    const auto *sll1_99 = buffer.data(sll1 + 99);
    const auto *sll1_102 = buffer.data(sll1 + 102);
    const auto *sll1_104 = buffer.data(sll1 + 104);
    const auto *sll1_107 = buffer.data(sll1 + 107);
    const auto *sll1_108 = buffer.data(sll1 + 108);
    const auto *sll1_110 = buffer.data(sll1 + 110);
    const auto *sll1_113 = buffer.data(sll1 + 113);
    const auto *sll1_114 = buffer.data(sll1 + 114);
    const auto *sll1_115 = buffer.data(sll1 + 115);
    const auto *sll1_117 = buffer.data(sll1 + 117);
    const auto *sll1_134 = buffer.data(sll1 + 134);

    const auto *smi0_79 = buffer.data(smi0 + 79);
    const auto *smi0_80 = buffer.data(smi0 + 80);
    const auto *smi0_81 = buffer.data(smi0 + 81);
    const auto *smi0_82 = buffer.data(smi0 + 82);
    const auto *smi0_83 = buffer.data(smi0 + 83);
    const auto *smi0_84 = buffer.data(smi0 + 84);
    const auto *smi0_87 = buffer.data(smi0 + 87);
    const auto *smi0_89 = buffer.data(smi0 + 89);
    const auto *smi0_90 = buffer.data(smi0 + 90);
    const auto *smi0_93 = buffer.data(smi0 + 93);
    const auto *smi0_94 = buffer.data(smi0 + 94);
    const auto *smi0_96 = buffer.data(smi0 + 96);
    const auto *smi0_98 = buffer.data(smi0 + 98);
    const auto *smi0_99 = buffer.data(smi0 + 99);
    const auto *smi0_101 = buffer.data(smi0 + 101);
    const auto *smi0_102 = buffer.data(smi0 + 102);
    const auto *smi0_104 = buffer.data(smi0 + 104);
    const auto *smi0_105 = buffer.data(smi0 + 105);
    const auto *smi0_107 = buffer.data(smi0 + 107);
    const auto *smi0_108 = buffer.data(smi0 + 108);
    const auto *smi0_109 = buffer.data(smi0 + 109);
    const auto *smi0_110 = buffer.data(smi0 + 110);
    const auto *smi0_111 = buffer.data(smi0 + 111);
    const auto *smi0_135 = buffer.data(smi0 + 135);
    const auto *smi0_136 = buffer.data(smi0 + 136);
    const auto *smi0_137 = buffer.data(smi0 + 137);
    const auto *smi0_138 = buffer.data(smi0 + 138);
    const auto *smi0_139 = buffer.data(smi0 + 139);
    const auto *smi0_140 = buffer.data(smi0 + 140);
    const auto *smi0_143 = buffer.data(smi0 + 143);
    const auto *smi0_145 = buffer.data(smi0 + 145);
    const auto *smi0_146 = buffer.data(smi0 + 146);

    const auto *smi1_79 = buffer.data(smi1 + 79);
    const auto *smi1_80 = buffer.data(smi1 + 80);
    const auto *smi1_81 = buffer.data(smi1 + 81);
    const auto *smi1_82 = buffer.data(smi1 + 82);
    const auto *smi1_83 = buffer.data(smi1 + 83);
    const auto *smi1_84 = buffer.data(smi1 + 84);
    const auto *smi1_87 = buffer.data(smi1 + 87);
    const auto *smi1_89 = buffer.data(smi1 + 89);
    const auto *smi1_90 = buffer.data(smi1 + 90);
    const auto *smi1_93 = buffer.data(smi1 + 93);
    const auto *smi1_94 = buffer.data(smi1 + 94);
    const auto *smi1_96 = buffer.data(smi1 + 96);
    const auto *smi1_98 = buffer.data(smi1 + 98);
    const auto *smi1_99 = buffer.data(smi1 + 99);
    const auto *smi1_101 = buffer.data(smi1 + 101);
    const auto *smi1_102 = buffer.data(smi1 + 102);
    const auto *smi1_104 = buffer.data(smi1 + 104);
    const auto *smi1_105 = buffer.data(smi1 + 105);
    const auto *smi1_107 = buffer.data(smi1 + 107);
    const auto *smi1_108 = buffer.data(smi1 + 108);
    const auto *smi1_109 = buffer.data(smi1 + 109);
    const auto *smi1_110 = buffer.data(smi1 + 110);
    const auto *smi1_111 = buffer.data(smi1 + 111);
    const auto *smi1_135 = buffer.data(smi1 + 135);
    const auto *smi1_136 = buffer.data(smi1 + 136);
    const auto *smi1_137 = buffer.data(smi1 + 137);
    const auto *smi1_138 = buffer.data(smi1 + 138);
    const auto *smi1_139 = buffer.data(smi1 + 139);
    const auto *smi1_140 = buffer.data(smi1 + 140);
    const auto *smi1_143 = buffer.data(smi1 + 143);
    const auto *smi1_145 = buffer.data(smi1 + 145);
    const auto *smi1_146 = buffer.data(smi1 + 146);

    const auto *smk_100 = buffer.data(smk + 100);
    const auto *smk_101 = buffer.data(smk + 101);
    const auto *smk_102 = buffer.data(smk + 102);
    const auto *smk_103 = buffer.data(smk + 103);
    const auto *smk_104 = buffer.data(smk + 104);
    const auto *smk_105 = buffer.data(smk + 105);
    const auto *smk_106 = buffer.data(smk + 106);
    const auto *smk_107 = buffer.data(smk + 107);
    const auto *smk_108 = buffer.data(smk + 108);
    const auto *smk_110 = buffer.data(smk + 110);
    const auto *smk_111 = buffer.data(smk + 111);
    const auto *smk_113 = buffer.data(smk + 113);
    const auto *smk_114 = buffer.data(smk + 114);
    const auto *smk_117 = buffer.data(smk + 117);
    const auto *smk_118 = buffer.data(smk + 118);
    const auto *smk_120 = buffer.data(smk + 120);
    const auto *smk_122 = buffer.data(smk + 122);
    const auto *smk_123 = buffer.data(smk + 123);
    const auto *smk_125 = buffer.data(smk + 125);
    const auto *smk_126 = buffer.data(smk + 126);
    const auto *smk_128 = buffer.data(smk + 128);
    const auto *smk_129 = buffer.data(smk + 129);
    const auto *smk_131 = buffer.data(smk + 131);
    const auto *smk_132 = buffer.data(smk + 132);
    const auto *smk_133 = buffer.data(smk + 133);
    const auto *smk_135 = buffer.data(smk + 135);
    const auto *smk_136 = buffer.data(smk + 136);
    const auto *smk_137 = buffer.data(smk + 137);
    const auto *smk_138 = buffer.data(smk + 138);
    const auto *smk_139 = buffer.data(smk + 139);
    const auto *smk_140 = buffer.data(smk + 140);
    const auto *smk_141 = buffer.data(smk + 141);
    const auto *smk_142 = buffer.data(smk + 142);
    const auto *smk_143 = buffer.data(smk + 143);
    const auto *smk_144 = buffer.data(smk + 144);
    const auto *smk_146 = buffer.data(smk + 146);
    const auto *smk_147 = buffer.data(smk + 147);
    const auto *smk_149 = buffer.data(smk + 149);
    const auto *smk_150 = buffer.data(smk + 150);
    const auto *smk_153 = buffer.data(smk + 153);
    const auto *smk_154 = buffer.data(smk + 154);
    const auto *smk_158 = buffer.data(smk + 158);
    const auto *smk_159 = buffer.data(smk + 159);
    const auto *smk_164 = buffer.data(smk + 164);
    const auto *smk_172 = buffer.data(smk + 172);
    const auto *smk_173 = buffer.data(smk + 173);
    const auto *smk_174 = buffer.data(smk + 174);
    const auto *smk_175 = buffer.data(smk + 175);
    const auto *smk_176 = buffer.data(smk + 176);
    const auto *smk_177 = buffer.data(smk + 177);
    const auto *smk_178 = buffer.data(smk + 178);
    const auto *smk_179 = buffer.data(smk + 179);
    const auto *smk_180 = buffer.data(smk + 180);
    const auto *smk_182 = buffer.data(smk + 182);
    const auto *smk_183 = buffer.data(smk + 183);
    const auto *smk_185 = buffer.data(smk + 185);
    const auto *smk_186 = buffer.data(smk + 186);

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pc_x, slk_100, slk_101, slk_102, \
                         slk_103, slk_104, smk_100, smk_101, smk_102, smk_103, \
                         smk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_21 * slk_100[k]
                   + f_3 * pc_x[k] * smk_100[k];

        t_119[k] = f_21 * slk_101[k]
                   + f_3 * pc_x[k] * smk_101[k];

        t_120[k] = f_21 * slk_102[k]
                   + f_3 * pc_x[k] * smk_102[k];

        t_121[k] = f_21 * slk_103[k]
                   + f_3 * pc_x[k] * smk_103[k];

        t_122[k] = f_21 * slk_104[k]
                   + f_3 * pc_x[k] * smk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_z, pc_x, pc_z, sll0_36, slk_105, \
                         slk_106, slk_107, sll1_36, smk_105, smk_106, \
                         smk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_21 * slk_105[k]
                   + f_3 * pc_x[k] * smk_105[k];

        t_124[k] = f_21 * slk_106[k]
                   + f_3 * pc_x[k] * smk_106[k];

        t_125[k] = f_21 * slk_107[k]
                   + f_3 * pc_x[k] * smk_107[k];

        t_126[k] = pb_z[k] * sll0_36[k]
                   - f_14 * pc_z[k] * sll1_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, pc_z, slk_28, smi0_79, smi0_80, smi1_79, \
                         smi1_80, smk_100, smk_102, smk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_15 * slk_28[k]
                   + f_3 * pc_z[k] * smk_100[k];

        t_128[k] = f_4 * smi0_79[k]
                   - f_5 * smi1_79[k]
                   + f_3 * pc_y[k] * smk_102[k];

        t_129[k] = f_6 * smi0_80[k]
                   - f_7 * smi1_80[k]
                   + f_3 * pc_y[k] * smk_103[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_y, smi0_81, smi0_82, smi0_83, smi1_81, \
                         smi1_82, smi1_83, smk_104, smk_105, smk_106, \
                         smk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * smi0_81[k]
                   - f_9 * smi1_81[k]
                   + f_3 * pc_y[k] * smk_104[k];

        t_131[k] = f_10 * smi0_82[k]
                   - f_11 * smi1_82[k]
                   + f_3 * pc_y[k] * smk_105[k];

        t_132[k] = f_12 * smi0_83[k]
                   - f_13 * smi1_83[k]
                   + f_3 * pc_y[k] * smk_106[k];

        t_133[k] = f_3 * pc_y[k] * smk_107[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, slk_35, slk_36, \
                         slk_108, smi0_83, smi0_84, smi1_83, smi1_84, smk_107, \
                         smk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * slk_35[k]
                   + f_1 * smi0_83[k]
                   - f_2 * smi1_83[k]
                   + f_3 * pc_z[k] * smk_107[k];

        t_135[k] = f_22 * slk_108[k]
                   + f_1 * smi0_84[k]
                   - f_2 * smi1_84[k]
                   + f_3 * pc_x[k] * smk_108[k];

        t_136[k] = f_16 * slk_36[k]
                   + f_3 * pc_y[k] * smk_108[k];

        t_137[k] = f_3 * pc_z[k] * smk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_y, slk_38, slk_111, slk_113, smi0_87, \
                         smi0_89, smi1_87, smi1_89, smk_110, smk_111, \
                         smk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_22 * slk_111[k]
                   + f_4 * smi0_87[k]
                   - f_5 * smi1_87[k]
                   + f_3 * pc_x[k] * smk_111[k];

        t_139[k] = f_16 * slk_38[k]
                   + f_3 * pc_y[k] * smk_110[k];

        t_140[k] = f_22 * slk_113[k]
                   + f_4 * smi0_89[k]
                   - f_5 * smi1_89[k]
                   + f_3 * pc_x[k] * smk_113[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_x, pc_y, pc_z, slk_41, slk_114, smi0_90, \
                         smi1_90, smk_111, smk_113, smk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_22 * slk_114[k]
                   + f_6 * smi0_90[k]
                   - f_7 * smi1_90[k]
                   + f_3 * pc_x[k] * smk_114[k];

        t_142[k] = f_3 * pc_z[k] * smk_111[k];

        t_143[k] = f_16 * slk_41[k]
                   + f_3 * pc_y[k] * smk_113[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_z, slk_117, slk_118, smi0_93, smi0_94, \
                         smi1_93, smi1_94, smk_114, smk_117, smk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_22 * slk_117[k]
                   + f_6 * smi0_93[k]
                   - f_7 * smi1_93[k]
                   + f_3 * pc_x[k] * smk_117[k];

        t_145[k] = f_22 * slk_118[k]
                   + f_8 * smi0_94[k]
                   - f_9 * smi1_94[k]
                   + f_3 * pc_x[k] * smk_118[k];

        t_146[k] = f_3 * pc_z[k] * smk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_x, pc_y, slk_45, slk_120, slk_122, smi0_96, \
                         smi0_98, smi1_96, smi1_98, smk_117, smk_120, \
                         smk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_22 * slk_120[k]
                   + f_8 * smi0_96[k]
                   - f_9 * smi1_96[k]
                   + f_3 * pc_x[k] * smk_120[k];

        t_148[k] = f_16 * slk_45[k]
                   + f_3 * pc_y[k] * smk_117[k];

        t_149[k] = f_22 * slk_122[k]
                   + f_8 * smi0_98[k]
                   - f_9 * smi1_98[k]
                   + f_3 * pc_x[k] * smk_122[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, slk_123, slk_125, smi0_99, smi0_101, \
                         smi1_99, smi1_101, smk_118, smk_123, smk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_22 * slk_123[k]
                   + f_10 * smi0_99[k]
                   - f_11 * smi1_99[k]
                   + f_3 * pc_x[k] * smk_123[k];

        t_151[k] = f_3 * pc_z[k] * smk_118[k];

        t_152[k] = f_22 * slk_125[k]
                   + f_10 * smi0_101[k]
                   - f_11 * smi1_101[k]
                   + f_3 * pc_x[k] * smk_125[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_x, pc_y, slk_50, slk_126, slk_128, smi0_102, \
                         smi0_104, smi1_102, smi1_104, smk_122, smk_126, \
                         smk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_22 * slk_126[k]
                   + f_10 * smi0_102[k]
                   - f_11 * smi1_102[k]
                   + f_3 * pc_x[k] * smk_126[k];

        t_154[k] = f_16 * slk_50[k]
                   + f_3 * pc_y[k] * smk_122[k];

        t_155[k] = f_22 * slk_128[k]
                   + f_10 * smi0_104[k]
                   - f_11 * smi1_104[k]
                   + f_3 * pc_x[k] * smk_128[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, slk_129, slk_131, smi0_105, \
                         smi0_107, smi1_105, smi1_107, smk_123, smk_129, \
                         smk_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_22 * slk_129[k]
                   + f_12 * smi0_105[k]
                   - f_13 * smi1_105[k]
                   + f_3 * pc_x[k] * smk_129[k];

        t_157[k] = f_3 * pc_z[k] * smk_123[k];

        t_158[k] = f_22 * slk_131[k]
                   + f_12 * smi0_107[k]
                   - f_13 * smi1_107[k]
                   + f_3 * pc_x[k] * smk_131[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, slk_56, slk_132, slk_133, smi0_108, \
                         smi0_109, smi1_108, smi1_109, smk_128, smk_132, \
                         smk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_22 * slk_132[k]
                   + f_12 * smi0_108[k]
                   - f_13 * smi1_108[k]
                   + f_3 * pc_x[k] * smk_132[k];

        t_160[k] = f_22 * slk_133[k]
                   + f_12 * smi0_109[k]
                   - f_13 * smi1_109[k]
                   + f_3 * pc_x[k] * smk_133[k];

        t_161[k] = f_16 * slk_56[k]
                   + f_3 * pc_y[k] * smk_128[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, slk_135, slk_136, slk_137, slk_138, \
                         smi0_111, smi1_111, smk_135, smk_136, smk_137, \
                         smk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_22 * slk_135[k]
                   + f_12 * smi0_111[k]
                   - f_13 * smi1_111[k]
                   + f_3 * pc_x[k] * smk_135[k];

        t_163[k] = f_22 * slk_136[k]
                   + f_3 * pc_x[k] * smk_136[k];

        t_164[k] = f_22 * slk_137[k]
                   + f_3 * pc_x[k] * smk_137[k];

        t_165[k] = f_22 * slk_138[k]
                   + f_3 * pc_x[k] * smk_138[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, slk_139, slk_140, slk_141, \
                         slk_142, slk_143, smk_139, smk_140, smk_141, smk_142, \
                         smk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_22 * slk_139[k]
                   + f_3 * pc_x[k] * smk_139[k];

        t_167[k] = f_22 * slk_140[k]
                   + f_3 * pc_x[k] * smk_140[k];

        t_168[k] = f_22 * slk_141[k]
                   + f_3 * pc_x[k] * smk_141[k];

        t_169[k] = f_22 * slk_142[k]
                   + f_3 * pc_x[k] * smk_142[k];

        t_170[k] = f_22 * slk_143[k]
                   + f_3 * pc_x[k] * smk_143[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pc_y, pc_z, slk_64, slk_66, smi0_105, smi0_107, \
                         smi1_105, smi1_107, smk_136, smk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_16 * slk_64[k]
                   + f_1 * smi0_105[k]
                   - f_2 * smi1_105[k]
                   + f_3 * pc_y[k] * smk_136[k];

        t_172[k] = f_3 * pc_z[k] * smk_136[k];

        t_173[k] = f_16 * slk_66[k]
                   + f_4 * smi0_107[k]
                   - f_5 * smi1_107[k]
                   + f_3 * pc_y[k] * smk_138[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, slk_67, slk_68, slk_69, smi0_108, \
                         smi0_109, smi0_110, smi1_108, smi1_109, smi1_110, smk_139, smk_140, \
                         smk_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_16 * slk_67[k]
                   + f_6 * smi0_108[k]
                   - f_7 * smi1_108[k]
                   + f_3 * pc_y[k] * smk_139[k];

        t_175[k] = f_16 * slk_68[k]
                   + f_8 * smi0_109[k]
                   - f_9 * smi1_109[k]
                   + f_3 * pc_y[k] * smk_140[k];

        t_176[k] = f_16 * slk_69[k]
                   + f_10 * smi0_110[k]
                   - f_11 * smi1_110[k]
                   + f_3 * pc_y[k] * smk_141[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pb_y, pc_y, pc_z, sll0_90, slk_70, \
                         slk_71, sll1_90, smi0_111, smi1_111, smk_142, \
                         smk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_16 * slk_70[k]
                   + f_12 * smi0_111[k]
                   - f_13 * smi1_111[k]
                   + f_3 * pc_y[k] * smk_142[k];

        t_178[k] = f_16 * slk_71[k]
                   + f_3 * pc_y[k] * smk_143[k];

        t_179[k] = f_1 * smi0_111[k]
                   - f_2 * smi1_111[k]
                   + f_3 * pc_z[k] * smk_143[k];

        t_180[k] = pb_y[k] * sll0_90[k]
                   - f_14 * pc_y[k] * sll1_90[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_z, pc_y, pc_z, sll0_48, slk_36, \
                         slk_72, slk_74, sll1_48, smk_144, smk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_15 * slk_72[k]
                   + f_3 * pc_y[k] * smk_144[k];

        t_182[k] = f_15 * slk_36[k]
                   + f_3 * pc_z[k] * smk_144[k];

        t_183[k] = pb_z[k] * sll0_48[k]
                   - f_14 * pc_z[k] * sll1_48[k];

        t_184[k] = f_15 * slk_74[k]
                   + f_3 * pc_y[k] * smk_146[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_y, pb_z, pc_y, pc_z, sll0_51, sll0_95, \
                         slk_39, slk_77, sll1_51, sll1_95, smk_147, \
                         smk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_y[k] * sll0_95[k]
                   - f_14 * pc_y[k] * sll1_95[k];

        t_186[k] = pb_z[k] * sll0_51[k]
                   - f_14 * pc_z[k] * sll1_51[k];

        t_187[k] = f_15 * slk_39[k]
                   + f_3 * pc_z[k] * smk_147[k];

        t_188[k] = f_15 * slk_77[k]
                   + f_3 * pc_y[k] * smk_149[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, pc_y, pc_z, sll0_55, sll0_99, \
                         slk_42, sll1_55, sll1_99, smk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_y[k] * sll0_99[k]
                   - f_14 * pc_y[k] * sll1_99[k];

        t_190[k] = pb_z[k] * sll0_55[k]
                   - f_14 * pc_z[k] * sll1_55[k];

        t_191[k] = f_15 * slk_42[k]
                   + f_3 * pc_z[k] * smk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pc_y, sll0_102, sll0_104, slk_80, slk_81, \
                         sll1_102, sll1_104, smk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_y[k] * sll0_102[k]
                   + f_16 * slk_80[k]
                   - f_14 * pc_y[k] * sll1_102[k];

        t_193[k] = f_15 * slk_81[k]
                   + f_3 * pc_y[k] * smk_153[k];

        t_194[k] = pb_y[k] * sll0_104[k]
                   - f_14 * pc_y[k] * sll1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_y, pb_z, pc_y, pc_z, sll0_60, sll0_107, \
                         slk_46, slk_84, sll1_60, sll1_107, smk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_z[k] * sll0_60[k]
                   - f_14 * pc_z[k] * sll1_60[k];

        t_196[k] = f_15 * slk_46[k]
                   + f_3 * pc_z[k] * smk_154[k];

        t_197[k] = pb_y[k] * sll0_107[k]
                   + f_17 * slk_84[k]
                   - f_14 * pc_y[k] * sll1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_y, pc_y, sll0_108, sll0_110, slk_85, slk_86, \
                         sll1_108, sll1_110, smk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_y[k] * sll0_108[k]
                   + f_16 * slk_85[k]
                   - f_14 * pc_y[k] * sll1_108[k];

        t_199[k] = f_15 * slk_86[k]
                   + f_3 * pc_y[k] * smk_158[k];

        t_200[k] = pb_y[k] * sll0_110[k]
                   - f_14 * pc_y[k] * sll1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_y, pb_z, pc_y, pc_z, sll0_66, sll0_113, \
                         slk_51, slk_89, sll1_66, sll1_113, smk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_z[k] * sll0_66[k]
                   - f_14 * pc_z[k] * sll1_66[k];

        t_202[k] = f_15 * slk_51[k]
                   + f_3 * pc_z[k] * smk_159[k];

        t_203[k] = pb_y[k] * sll0_113[k]
                   + f_18 * slk_89[k]
                   - f_14 * pc_y[k] * sll1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pc_y, sll0_114, sll0_115, sll0_117, \
                         slk_90, slk_91, slk_92, sll1_114, sll1_115, sll1_117, \
                         smk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_y[k] * sll0_114[k]
                   + f_17 * slk_90[k]
                   - f_14 * pc_y[k] * sll1_114[k];

        t_205[k] = pb_y[k] * sll0_115[k]
                   + f_16 * slk_91[k]
                   - f_14 * pc_y[k] * sll1_115[k];

        t_206[k] = f_15 * slk_92[k]
                   + f_3 * pc_y[k] * smk_164[k];

        t_207[k] = pb_y[k] * sll0_117[k]
                   - f_14 * pc_y[k] * sll1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, slk_172, slk_173, slk_174, \
                         slk_175, slk_176, smk_172, smk_173, smk_174, smk_175, \
                         smk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_22 * slk_172[k]
                   + f_3 * pc_x[k] * smk_172[k];

        t_209[k] = f_22 * slk_173[k]
                   + f_3 * pc_x[k] * smk_173[k];

        t_210[k] = f_22 * slk_174[k]
                   + f_3 * pc_x[k] * smk_174[k];

        t_211[k] = f_22 * slk_175[k]
                   + f_3 * pc_x[k] * smk_175[k];

        t_212[k] = f_22 * slk_176[k]
                   + f_3 * pc_x[k] * smk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_z, pc_x, pc_z, sll0_81, slk_177, \
                         slk_178, slk_179, sll1_81, smk_177, smk_178, \
                         smk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_22 * slk_177[k]
                   + f_3 * pc_x[k] * smk_177[k];

        t_214[k] = f_22 * slk_178[k]
                   + f_3 * pc_x[k] * smk_178[k];

        t_215[k] = f_22 * slk_179[k]
                   + f_3 * pc_x[k] * smk_179[k];

        t_216[k] = pb_z[k] * sll0_81[k]
                   - f_14 * pc_z[k] * sll1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, slk_64, slk_102, slk_103, smi0_135, \
                         smi0_136, smi1_135, smi1_136, smk_172, smk_174, \
                         smk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * slk_64[k]
                   + f_3 * pc_z[k] * smk_172[k];

        t_218[k] = f_15 * slk_102[k]
                   + f_4 * smi0_135[k]
                   - f_5 * smi1_135[k]
                   + f_3 * pc_y[k] * smk_174[k];

        t_219[k] = f_15 * slk_103[k]
                   + f_6 * smi0_136[k]
                   - f_7 * smi1_136[k]
                   + f_3 * pc_y[k] * smk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, slk_104, slk_105, slk_106, smi0_137, \
                         smi0_138, smi0_139, smi1_137, smi1_138, smi1_139, smk_176, smk_177, \
                         smk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * slk_104[k]
                   + f_8 * smi0_137[k]
                   - f_9 * smi1_137[k]
                   + f_3 * pc_y[k] * smk_176[k];

        t_221[k] = f_15 * slk_105[k]
                   + f_10 * smi0_138[k]
                   - f_11 * smi1_138[k]
                   + f_3 * pc_y[k] * smk_177[k];

        t_222[k] = f_15 * slk_106[k]
                   + f_12 * smi0_139[k]
                   - f_13 * smi1_139[k]
                   + f_3 * pc_y[k] * smk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pb_y, pc_x, pc_y, sll0_134, slk_107, \
                         slk_180, sll1_134, smi0_140, smi1_140, smk_179, \
                         smk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * slk_107[k]
                   + f_3 * pc_y[k] * smk_179[k];

        t_224[k] = pb_y[k] * sll0_134[k]
                   - f_14 * pc_y[k] * sll1_134[k];

        t_225[k] = f_22 * slk_180[k]
                   + f_1 * smi0_140[k]
                   - f_2 * smi1_140[k]
                   + f_3 * pc_x[k] * smk_180[k];

        t_226[k] = f_3 * pc_y[k] * smk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, slk_72, slk_183, smi0_143, \
                         smi1_143, smk_180, smk_182, smk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * slk_72[k]
                   + f_3 * pc_z[k] * smk_180[k];

        t_228[k] = f_22 * slk_183[k]
                   + f_4 * smi0_143[k]
                   - f_5 * smi1_143[k]
                   + f_3 * pc_x[k] * smk_183[k];

        t_229[k] = f_3 * pc_y[k] * smk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pc_x, pc_z, slk_75, slk_185, slk_186, smi0_145, \
                         smi0_146, smi1_145, smi1_146, smk_183, smk_185, \
                         smk_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_22 * slk_185[k]
                   + f_4 * smi0_145[k]
                   - f_5 * smi1_145[k]
                   + f_3 * pc_x[k] * smk_185[k];

        t_231[k] = f_22 * slk_186[k]
                   + f_6 * smi0_146[k]
                   - f_7 * smi1_146[k]
                   + f_3 * pc_x[k] * smk_186[k];

        t_232[k] = f_16 * slk_75[k]
                   + f_3 * pc_z[k] * smk_183[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_135 = buffer.data(sll0 + 135);
    const auto *sll0_138 = buffer.data(sll0 + 138);
    const auto *sll0_141 = buffer.data(sll0 + 141);
    const auto *sll0_145 = buffer.data(sll0 + 145);
    const auto *sll0_147 = buffer.data(sll0 + 147);
    const auto *sll0_150 = buffer.data(sll0 + 150);
    const auto *sll0_152 = buffer.data(sll0 + 152);
    const auto *sll0_153 = buffer.data(sll0 + 153);
    const auto *sll0_156 = buffer.data(sll0 + 156);
    const auto *sll0_158 = buffer.data(sll0 + 158);
    const auto *sll0_159 = buffer.data(sll0 + 159);
    const auto *sll0_160 = buffer.data(sll0 + 160);

    const auto *slk_78 = buffer.data(slk + 78);
    const auto *slk_82 = buffer.data(slk + 82);
    const auto *slk_87 = buffer.data(slk + 87);
    const auto *slk_100 = buffer.data(slk + 100);
    const auto *slk_107 = buffer.data(slk + 107);
    const auto *slk_108 = buffer.data(slk + 108);
    const auto *slk_110 = buffer.data(slk + 110);
    const auto *slk_111 = buffer.data(slk + 111);
    const auto *slk_113 = buffer.data(slk + 113);
    const auto *slk_114 = buffer.data(slk + 114);
    const auto *slk_115 = buffer.data(slk + 115);
    const auto *slk_117 = buffer.data(slk + 117);
    const auto *slk_118 = buffer.data(slk + 118);
    const auto *slk_119 = buffer.data(slk + 119);
    const auto *slk_120 = buffer.data(slk + 120);
    const auto *slk_122 = buffer.data(slk + 122);
    const auto *slk_123 = buffer.data(slk + 123);
    const auto *slk_124 = buffer.data(slk + 124);
    const auto *slk_125 = buffer.data(slk + 125);
    const auto *slk_126 = buffer.data(slk + 126);
    const auto *slk_128 = buffer.data(slk + 128);
    const auto *slk_136 = buffer.data(slk + 136);
    const auto *slk_138 = buffer.data(slk + 138);
    const auto *slk_139 = buffer.data(slk + 139);
    const auto *slk_140 = buffer.data(slk + 140);
    const auto *slk_141 = buffer.data(slk + 141);
    const auto *slk_142 = buffer.data(slk + 142);
    const auto *slk_143 = buffer.data(slk + 143);
    const auto *slk_144 = buffer.data(slk + 144);
    const auto *slk_146 = buffer.data(slk + 146);
    const auto *slk_149 = buffer.data(slk + 149);
    const auto *slk_153 = buffer.data(slk + 153);
    const auto *slk_158 = buffer.data(slk + 158);
    const auto *slk_189 = buffer.data(slk + 189);
    const auto *slk_190 = buffer.data(slk + 190);
    const auto *slk_192 = buffer.data(slk + 192);
    const auto *slk_194 = buffer.data(slk + 194);
    const auto *slk_195 = buffer.data(slk + 195);
    const auto *slk_197 = buffer.data(slk + 197);
    const auto *slk_198 = buffer.data(slk + 198);
    const auto *slk_200 = buffer.data(slk + 200);
    const auto *slk_201 = buffer.data(slk + 201);
    const auto *slk_203 = buffer.data(slk + 203);
    const auto *slk_204 = buffer.data(slk + 204);
    const auto *slk_205 = buffer.data(slk + 205);
    const auto *slk_207 = buffer.data(slk + 207);
    const auto *slk_208 = buffer.data(slk + 208);
    const auto *slk_209 = buffer.data(slk + 209);
    const auto *slk_210 = buffer.data(slk + 210);
    const auto *slk_211 = buffer.data(slk + 211);
    const auto *slk_212 = buffer.data(slk + 212);
    const auto *slk_213 = buffer.data(slk + 213);
    const auto *slk_214 = buffer.data(slk + 214);
    const auto *slk_215 = buffer.data(slk + 215);
    const auto *slk_216 = buffer.data(slk + 216);
    const auto *slk_219 = buffer.data(slk + 219);
    const auto *slk_221 = buffer.data(slk + 221);
    const auto *slk_222 = buffer.data(slk + 222);
    const auto *slk_225 = buffer.data(slk + 225);
    const auto *slk_226 = buffer.data(slk + 226);
    const auto *slk_228 = buffer.data(slk + 228);
    const auto *slk_230 = buffer.data(slk + 230);
    const auto *slk_231 = buffer.data(slk + 231);
    const auto *slk_233 = buffer.data(slk + 233);
    const auto *slk_234 = buffer.data(slk + 234);
    const auto *slk_236 = buffer.data(slk + 236);
    const auto *slk_237 = buffer.data(slk + 237);
    const auto *slk_239 = buffer.data(slk + 239);
    const auto *slk_240 = buffer.data(slk + 240);
    const auto *slk_241 = buffer.data(slk + 241);
    const auto *slk_243 = buffer.data(slk + 243);
    const auto *slk_244 = buffer.data(slk + 244);
    const auto *slk_245 = buffer.data(slk + 245);
    const auto *slk_246 = buffer.data(slk + 246);
    const auto *slk_247 = buffer.data(slk + 247);
    const auto *slk_248 = buffer.data(slk + 248);
    const auto *slk_249 = buffer.data(slk + 249);
    const auto *slk_250 = buffer.data(slk + 250);
    const auto *slk_251 = buffer.data(slk + 251);
    const auto *slk_257 = buffer.data(slk + 257);
    const auto *slk_261 = buffer.data(slk + 261);
    const auto *slk_266 = buffer.data(slk + 266);
    const auto *slk_272 = buffer.data(slk + 272);

    const auto *sll1_135 = buffer.data(sll1 + 135);
    const auto *sll1_138 = buffer.data(sll1 + 138);
    const auto *sll1_141 = buffer.data(sll1 + 141);
    const auto *sll1_145 = buffer.data(sll1 + 145);
    const auto *sll1_147 = buffer.data(sll1 + 147);
    const auto *sll1_150 = buffer.data(sll1 + 150);
    const auto *sll1_152 = buffer.data(sll1 + 152);
    const auto *sll1_153 = buffer.data(sll1 + 153);
    const auto *sll1_156 = buffer.data(sll1 + 156);
    const auto *sll1_158 = buffer.data(sll1 + 158);
    const auto *sll1_159 = buffer.data(sll1 + 159);
    const auto *sll1_160 = buffer.data(sll1 + 160);

    const auto *smi0_149 = buffer.data(smi0 + 149);
    const auto *smi0_150 = buffer.data(smi0 + 150);
    const auto *smi0_152 = buffer.data(smi0 + 152);
    const auto *smi0_154 = buffer.data(smi0 + 154);
    const auto *smi0_155 = buffer.data(smi0 + 155);
    const auto *smi0_157 = buffer.data(smi0 + 157);
    const auto *smi0_158 = buffer.data(smi0 + 158);
    const auto *smi0_160 = buffer.data(smi0 + 160);
    const auto *smi0_161 = buffer.data(smi0 + 161);
    const auto *smi0_163 = buffer.data(smi0 + 163);
    const auto *smi0_164 = buffer.data(smi0 + 164);
    const auto *smi0_165 = buffer.data(smi0 + 165);
    const auto *smi0_166 = buffer.data(smi0 + 166);
    const auto *smi0_167 = buffer.data(smi0 + 167);
    const auto *smi0_168 = buffer.data(smi0 + 168);
    const auto *smi0_171 = buffer.data(smi0 + 171);
    const auto *smi0_173 = buffer.data(smi0 + 173);
    const auto *smi0_174 = buffer.data(smi0 + 174);
    const auto *smi0_177 = buffer.data(smi0 + 177);
    const auto *smi0_178 = buffer.data(smi0 + 178);
    const auto *smi0_180 = buffer.data(smi0 + 180);
    const auto *smi0_182 = buffer.data(smi0 + 182);
    const auto *smi0_183 = buffer.data(smi0 + 183);
    const auto *smi0_185 = buffer.data(smi0 + 185);
    const auto *smi0_186 = buffer.data(smi0 + 186);
    const auto *smi0_188 = buffer.data(smi0 + 188);
    const auto *smi0_189 = buffer.data(smi0 + 189);
    const auto *smi0_191 = buffer.data(smi0 + 191);
    const auto *smi0_192 = buffer.data(smi0 + 192);
    const auto *smi0_193 = buffer.data(smi0 + 193);
    const auto *smi0_194 = buffer.data(smi0 + 194);
    const auto *smi0_195 = buffer.data(smi0 + 195);
    const auto *smi0_201 = buffer.data(smi0 + 201);
    const auto *smi0_205 = buffer.data(smi0 + 205);
    const auto *smi0_210 = buffer.data(smi0 + 210);
    const auto *smi0_216 = buffer.data(smi0 + 216);

    const auto *smi1_149 = buffer.data(smi1 + 149);
    const auto *smi1_150 = buffer.data(smi1 + 150);
    const auto *smi1_152 = buffer.data(smi1 + 152);
    const auto *smi1_154 = buffer.data(smi1 + 154);
    const auto *smi1_155 = buffer.data(smi1 + 155);
    const auto *smi1_157 = buffer.data(smi1 + 157);
    const auto *smi1_158 = buffer.data(smi1 + 158);
    const auto *smi1_160 = buffer.data(smi1 + 160);
    const auto *smi1_161 = buffer.data(smi1 + 161);
    const auto *smi1_163 = buffer.data(smi1 + 163);
    const auto *smi1_164 = buffer.data(smi1 + 164);
    const auto *smi1_165 = buffer.data(smi1 + 165);
    const auto *smi1_166 = buffer.data(smi1 + 166);
    const auto *smi1_167 = buffer.data(smi1 + 167);
    const auto *smi1_168 = buffer.data(smi1 + 168);
    const auto *smi1_171 = buffer.data(smi1 + 171);
    const auto *smi1_173 = buffer.data(smi1 + 173);
    const auto *smi1_174 = buffer.data(smi1 + 174);
    const auto *smi1_177 = buffer.data(smi1 + 177);
    const auto *smi1_178 = buffer.data(smi1 + 178);
    const auto *smi1_180 = buffer.data(smi1 + 180);
    const auto *smi1_182 = buffer.data(smi1 + 182);
    const auto *smi1_183 = buffer.data(smi1 + 183);
    const auto *smi1_185 = buffer.data(smi1 + 185);
    const auto *smi1_186 = buffer.data(smi1 + 186);
    const auto *smi1_188 = buffer.data(smi1 + 188);
    const auto *smi1_189 = buffer.data(smi1 + 189);
    const auto *smi1_191 = buffer.data(smi1 + 191);
    const auto *smi1_192 = buffer.data(smi1 + 192);
    const auto *smi1_193 = buffer.data(smi1 + 193);
    const auto *smi1_194 = buffer.data(smi1 + 194);
    const auto *smi1_195 = buffer.data(smi1 + 195);
    const auto *smi1_201 = buffer.data(smi1 + 201);
    const auto *smi1_205 = buffer.data(smi1 + 205);
    const auto *smi1_210 = buffer.data(smi1 + 210);
    const auto *smi1_216 = buffer.data(smi1 + 216);

    const auto *smk_185 = buffer.data(smk + 185);
    const auto *smk_186 = buffer.data(smk + 186);
    const auto *smk_189 = buffer.data(smk + 189);
    const auto *smk_190 = buffer.data(smk + 190);
    const auto *smk_192 = buffer.data(smk + 192);
    const auto *smk_194 = buffer.data(smk + 194);
    const auto *smk_195 = buffer.data(smk + 195);
    const auto *smk_197 = buffer.data(smk + 197);
    const auto *smk_198 = buffer.data(smk + 198);
    const auto *smk_200 = buffer.data(smk + 200);
    const auto *smk_201 = buffer.data(smk + 201);
    const auto *smk_203 = buffer.data(smk + 203);
    const auto *smk_204 = buffer.data(smk + 204);
    const auto *smk_205 = buffer.data(smk + 205);
    const auto *smk_207 = buffer.data(smk + 207);
    const auto *smk_208 = buffer.data(smk + 208);
    const auto *smk_209 = buffer.data(smk + 209);
    const auto *smk_210 = buffer.data(smk + 210);
    const auto *smk_211 = buffer.data(smk + 211);
    const auto *smk_212 = buffer.data(smk + 212);
    const auto *smk_213 = buffer.data(smk + 213);
    const auto *smk_214 = buffer.data(smk + 214);
    const auto *smk_215 = buffer.data(smk + 215);
    const auto *smk_216 = buffer.data(smk + 216);
    const auto *smk_218 = buffer.data(smk + 218);
    const auto *smk_219 = buffer.data(smk + 219);
    const auto *smk_221 = buffer.data(smk + 221);
    const auto *smk_222 = buffer.data(smk + 222);
    const auto *smk_225 = buffer.data(smk + 225);
    const auto *smk_226 = buffer.data(smk + 226);
    const auto *smk_228 = buffer.data(smk + 228);
    const auto *smk_230 = buffer.data(smk + 230);
    const auto *smk_231 = buffer.data(smk + 231);
    const auto *smk_233 = buffer.data(smk + 233);
    const auto *smk_234 = buffer.data(smk + 234);
    const auto *smk_236 = buffer.data(smk + 236);
    const auto *smk_237 = buffer.data(smk + 237);
    const auto *smk_239 = buffer.data(smk + 239);
    const auto *smk_240 = buffer.data(smk + 240);
    const auto *smk_241 = buffer.data(smk + 241);
    const auto *smk_243 = buffer.data(smk + 243);
    const auto *smk_244 = buffer.data(smk + 244);
    const auto *smk_245 = buffer.data(smk + 245);
    const auto *smk_246 = buffer.data(smk + 246);
    const auto *smk_247 = buffer.data(smk + 247);
    const auto *smk_248 = buffer.data(smk + 248);
    const auto *smk_249 = buffer.data(smk + 249);
    const auto *smk_250 = buffer.data(smk + 250);
    const auto *smk_251 = buffer.data(smk + 251);
    const auto *smk_252 = buffer.data(smk + 252);
    const auto *smk_254 = buffer.data(smk + 254);
    const auto *smk_255 = buffer.data(smk + 255);
    const auto *smk_257 = buffer.data(smk + 257);
    const auto *smk_258 = buffer.data(smk + 258);
    const auto *smk_261 = buffer.data(smk + 261);
    const auto *smk_262 = buffer.data(smk + 262);
    const auto *smk_266 = buffer.data(smk + 266);
    const auto *smk_267 = buffer.data(smk + 267);
    const auto *smk_272 = buffer.data(smk + 272);

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, pc_y, slk_189, slk_190, smi0_149, \
                         smi0_150, smi1_149, smi1_150, smk_185, smk_189, \
                         smk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * pc_y[k] * smk_185[k];

        t_234[k] = f_22 * slk_189[k]
                   + f_6 * smi0_149[k]
                   - f_7 * smi1_149[k]
                   + f_3 * pc_x[k] * smk_189[k];

        t_235[k] = f_22 * slk_190[k]
                   + f_8 * smi0_150[k]
                   - f_9 * smi1_150[k]
                   + f_3 * pc_x[k] * smk_190[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pc_x, pc_y, pc_z, slk_78, slk_192, smi0_152, \
                         smi1_152, smk_186, smk_189, smk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_16 * slk_78[k]
                   + f_3 * pc_z[k] * smk_186[k];

        t_237[k] = f_22 * slk_192[k]
                   + f_8 * smi0_152[k]
                   - f_9 * smi1_152[k]
                   + f_3 * pc_x[k] * smk_192[k];

        t_238[k] = f_3 * pc_y[k] * smk_189[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, slk_82, slk_194, slk_195, smi0_154, \
                         smi0_155, smi1_154, smi1_155, smk_190, smk_194, \
                         smk_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_22 * slk_194[k]
                   + f_8 * smi0_154[k]
                   - f_9 * smi1_154[k]
                   + f_3 * pc_x[k] * smk_194[k];

        t_240[k] = f_22 * slk_195[k]
                   + f_10 * smi0_155[k]
                   - f_11 * smi1_155[k]
                   + f_3 * pc_x[k] * smk_195[k];

        t_241[k] = f_16 * slk_82[k]
                   + f_3 * pc_z[k] * smk_190[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pc_x, pc_y, slk_197, slk_198, smi0_157, \
                         smi0_158, smi1_157, smi1_158, smk_194, smk_197, \
                         smk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_22 * slk_197[k]
                   + f_10 * smi0_157[k]
                   - f_11 * smi1_157[k]
                   + f_3 * pc_x[k] * smk_197[k];

        t_243[k] = f_22 * slk_198[k]
                   + f_10 * smi0_158[k]
                   - f_11 * smi1_158[k]
                   + f_3 * pc_x[k] * smk_198[k];

        t_244[k] = f_3 * pc_y[k] * smk_194[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_x, pc_z, slk_87, slk_200, slk_201, smi0_160, \
                         smi0_161, smi1_160, smi1_161, smk_195, smk_200, \
                         smk_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_22 * slk_200[k]
                   + f_10 * smi0_160[k]
                   - f_11 * smi1_160[k]
                   + f_3 * pc_x[k] * smk_200[k];

        t_246[k] = f_22 * slk_201[k]
                   + f_12 * smi0_161[k]
                   - f_13 * smi1_161[k]
                   + f_3 * pc_x[k] * smk_201[k];

        t_247[k] = f_16 * slk_87[k]
                   + f_3 * pc_z[k] * smk_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, slk_203, slk_204, slk_205, smi0_163, \
                         smi0_164, smi0_165, smi1_163, smi1_164, smi1_165, smk_203, smk_204, \
                         smk_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_22 * slk_203[k]
                   + f_12 * smi0_163[k]
                   - f_13 * smi1_163[k]
                   + f_3 * pc_x[k] * smk_203[k];

        t_249[k] = f_22 * slk_204[k]
                   + f_12 * smi0_164[k]
                   - f_13 * smi1_164[k]
                   + f_3 * pc_x[k] * smk_204[k];

        t_250[k] = f_22 * slk_205[k]
                   + f_12 * smi0_165[k]
                   - f_13 * smi1_165[k]
                   + f_3 * pc_x[k] * smk_205[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, slk_207, slk_208, slk_209, \
                         smi0_167, smi1_167, smk_200, smk_207, smk_208, \
                         smk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_3 * pc_y[k] * smk_200[k];

        t_252[k] = f_22 * slk_207[k]
                   + f_12 * smi0_167[k]
                   - f_13 * smi1_167[k]
                   + f_3 * pc_x[k] * smk_207[k];

        t_253[k] = f_22 * slk_208[k]
                   + f_3 * pc_x[k] * smk_208[k];

        t_254[k] = f_22 * slk_209[k]
                   + f_3 * pc_x[k] * smk_209[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, slk_210, slk_211, slk_212, \
                         slk_213, slk_214, smk_210, smk_211, smk_212, smk_213, \
                         smk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_22 * slk_210[k]
                   + f_3 * pc_x[k] * smk_210[k];

        t_256[k] = f_22 * slk_211[k]
                   + f_3 * pc_x[k] * smk_211[k];

        t_257[k] = f_22 * slk_212[k]
                   + f_3 * pc_x[k] * smk_212[k];

        t_258[k] = f_22 * slk_213[k]
                   + f_3 * pc_x[k] * smk_213[k];

        t_259[k] = f_22 * slk_214[k]
                   + f_3 * pc_x[k] * smk_214[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, pc_z, slk_100, slk_215, \
                         smi0_161, smi0_163, smi1_161, smi1_163, smk_208, smk_210, \
                         smk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_22 * slk_215[k]
                   + f_3 * pc_x[k] * smk_215[k];

        t_261[k] = f_1 * smi0_161[k]
                   - f_2 * smi1_161[k]
                   + f_3 * pc_y[k] * smk_208[k];

        t_262[k] = f_16 * slk_100[k]
                   + f_3 * pc_z[k] * smk_208[k];

        t_263[k] = f_4 * smi0_163[k]
                   - f_5 * smi1_163[k]
                   + f_3 * pc_y[k] * smk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, smi0_164, smi0_165, smi0_166, smi1_164, \
                         smi1_165, smi1_166, smk_211, smk_212, \
                         smk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_6 * smi0_164[k]
                   - f_7 * smi1_164[k]
                   + f_3 * pc_y[k] * smk_211[k];

        t_265[k] = f_8 * smi0_165[k]
                   - f_9 * smi1_165[k]
                   + f_3 * pc_y[k] * smk_212[k];

        t_266[k] = f_10 * smi0_166[k]
                   - f_11 * smi1_166[k]
                   + f_3 * pc_y[k] * smk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, slk_107, slk_216, \
                         smi0_167, smi0_168, smi1_167, smi1_168, smk_214, smk_215, \
                         smk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * smi0_167[k]
                   - f_13 * smi1_167[k]
                   + f_3 * pc_y[k] * smk_214[k];

        t_268[k] = f_3 * pc_y[k] * smk_215[k];

        t_269[k] = f_16 * slk_107[k]
                   + f_1 * smi0_167[k]
                   - f_2 * smi1_167[k]
                   + f_3 * pc_z[k] * smk_215[k];

        t_270[k] = f_20 * slk_216[k]
                   + f_1 * smi0_168[k]
                   - f_2 * smi1_168[k]
                   + f_3 * pc_x[k] * smk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, slk_108, slk_110, \
                         slk_219, smi0_171, smi1_171, smk_216, smk_218, \
                         smk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * slk_108[k]
                   + f_3 * pc_y[k] * smk_216[k];

        t_272[k] = f_3 * pc_z[k] * smk_216[k];

        t_273[k] = f_20 * slk_219[k]
                   + f_4 * smi0_171[k]
                   - f_5 * smi1_171[k]
                   + f_3 * pc_x[k] * smk_219[k];

        t_274[k] = f_17 * slk_110[k]
                   + f_3 * pc_y[k] * smk_218[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, slk_221, slk_222, smi0_173, \
                         smi0_174, smi1_173, smi1_174, smk_219, smk_221, \
                         smk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_20 * slk_221[k]
                   + f_4 * smi0_173[k]
                   - f_5 * smi1_173[k]
                   + f_3 * pc_x[k] * smk_221[k];

        t_276[k] = f_20 * slk_222[k]
                   + f_6 * smi0_174[k]
                   - f_7 * smi1_174[k]
                   + f_3 * pc_x[k] * smk_222[k];

        t_277[k] = f_3 * pc_z[k] * smk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pc_x, pc_y, slk_113, slk_225, slk_226, smi0_177, \
                         smi0_178, smi1_177, smi1_178, smk_221, smk_225, \
                         smk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * slk_113[k]
                   + f_3 * pc_y[k] * smk_221[k];

        t_279[k] = f_20 * slk_225[k]
                   + f_6 * smi0_177[k]
                   - f_7 * smi1_177[k]
                   + f_3 * pc_x[k] * smk_225[k];

        t_280[k] = f_20 * slk_226[k]
                   + f_8 * smi0_178[k]
                   - f_9 * smi1_178[k]
                   + f_3 * pc_x[k] * smk_226[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pc_x, pc_y, pc_z, slk_117, slk_228, smi0_180, \
                         smi1_180, smk_222, smk_225, smk_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * pc_z[k] * smk_222[k];

        t_282[k] = f_20 * slk_228[k]
                   + f_8 * smi0_180[k]
                   - f_9 * smi1_180[k]
                   + f_3 * pc_x[k] * smk_228[k];

        t_283[k] = f_17 * slk_117[k]
                   + f_3 * pc_y[k] * smk_225[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_z, slk_230, slk_231, smi0_182, \
                         smi0_183, smi1_182, smi1_183, smk_226, smk_230, \
                         smk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_20 * slk_230[k]
                   + f_8 * smi0_182[k]
                   - f_9 * smi1_182[k]
                   + f_3 * pc_x[k] * smk_230[k];

        t_285[k] = f_20 * slk_231[k]
                   + f_10 * smi0_183[k]
                   - f_11 * smi1_183[k]
                   + f_3 * pc_x[k] * smk_231[k];

        t_286[k] = f_3 * pc_z[k] * smk_226[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, slk_122, slk_233, slk_234, smi0_185, \
                         smi0_186, smi1_185, smi1_186, smk_230, smk_233, \
                         smk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_20 * slk_233[k]
                   + f_10 * smi0_185[k]
                   - f_11 * smi1_185[k]
                   + f_3 * pc_x[k] * smk_233[k];

        t_288[k] = f_20 * slk_234[k]
                   + f_10 * smi0_186[k]
                   - f_11 * smi1_186[k]
                   + f_3 * pc_x[k] * smk_234[k];

        t_289[k] = f_17 * slk_122[k]
                   + f_3 * pc_y[k] * smk_230[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, slk_236, slk_237, smi0_188, \
                         smi0_189, smi1_188, smi1_189, smk_231, smk_236, \
                         smk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_20 * slk_236[k]
                   + f_10 * smi0_188[k]
                   - f_11 * smi1_188[k]
                   + f_3 * pc_x[k] * smk_236[k];

        t_291[k] = f_20 * slk_237[k]
                   + f_12 * smi0_189[k]
                   - f_13 * smi1_189[k]
                   + f_3 * pc_x[k] * smk_237[k];

        t_292[k] = f_3 * pc_z[k] * smk_231[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pc_x, slk_239, slk_240, slk_241, smi0_191, \
                         smi0_192, smi0_193, smi1_191, smi1_192, smi1_193, smk_239, smk_240, \
                         smk_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_20 * slk_239[k]
                   + f_12 * smi0_191[k]
                   - f_13 * smi1_191[k]
                   + f_3 * pc_x[k] * smk_239[k];

        t_294[k] = f_20 * slk_240[k]
                   + f_12 * smi0_192[k]
                   - f_13 * smi1_192[k]
                   + f_3 * pc_x[k] * smk_240[k];

        t_295[k] = f_20 * slk_241[k]
                   + f_12 * smi0_193[k]
                   - f_13 * smi1_193[k]
                   + f_3 * pc_x[k] * smk_241[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pc_x, pc_y, slk_128, slk_243, slk_244, \
                         slk_245, smi0_195, smi1_195, smk_236, smk_243, smk_244, \
                         smk_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_17 * slk_128[k]
                   + f_3 * pc_y[k] * smk_236[k];

        t_297[k] = f_20 * slk_243[k]
                   + f_12 * smi0_195[k]
                   - f_13 * smi1_195[k]
                   + f_3 * pc_x[k] * smk_243[k];

        t_298[k] = f_20 * slk_244[k]
                   + f_3 * pc_x[k] * smk_244[k];

        t_299[k] = f_20 * slk_245[k]
                   + f_3 * pc_x[k] * smk_245[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, slk_246, slk_247, slk_248, \
                         slk_249, slk_250, smk_246, smk_247, smk_248, smk_249, \
                         smk_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_20 * slk_246[k]
                   + f_3 * pc_x[k] * smk_246[k];

        t_301[k] = f_20 * slk_247[k]
                   + f_3 * pc_x[k] * smk_247[k];

        t_302[k] = f_20 * slk_248[k]
                   + f_3 * pc_x[k] * smk_248[k];

        t_303[k] = f_20 * slk_249[k]
                   + f_3 * pc_x[k] * smk_249[k];

        t_304[k] = f_20 * slk_250[k]
                   + f_3 * pc_x[k] * smk_250[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pc_x, pc_y, pc_z, slk_136, slk_251, smi0_189, \
                         smi1_189, smk_244, smk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_20 * slk_251[k]
                   + f_3 * pc_x[k] * smk_251[k];

        t_306[k] = f_17 * slk_136[k]
                   + f_1 * smi0_189[k]
                   - f_2 * smi1_189[k]
                   + f_3 * pc_y[k] * smk_244[k];

        t_307[k] = f_3 * pc_z[k] * smk_244[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_y, slk_138, slk_139, slk_140, smi0_191, \
                         smi0_192, smi0_193, smi1_191, smi1_192, smi1_193, smk_246, smk_247, \
                         smk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_17 * slk_138[k]
                   + f_4 * smi0_191[k]
                   - f_5 * smi1_191[k]
                   + f_3 * pc_y[k] * smk_246[k];

        t_309[k] = f_17 * slk_139[k]
                   + f_6 * smi0_192[k]
                   - f_7 * smi1_192[k]
                   + f_3 * pc_y[k] * smk_247[k];

        t_310[k] = f_17 * slk_140[k]
                   + f_8 * smi0_193[k]
                   - f_9 * smi1_193[k]
                   + f_3 * pc_y[k] * smk_248[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, slk_141, slk_142, slk_143, \
                         smi0_194, smi0_195, smi1_194, smi1_195, smk_249, smk_250, \
                         smk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_17 * slk_141[k]
                   + f_10 * smi0_194[k]
                   - f_11 * smi1_194[k]
                   + f_3 * pc_y[k] * smk_249[k];

        t_312[k] = f_17 * slk_142[k]
                   + f_12 * smi0_195[k]
                   - f_13 * smi1_195[k]
                   + f_3 * pc_y[k] * smk_250[k];

        t_313[k] = f_17 * slk_143[k]
                   + f_3 * pc_y[k] * smk_251[k];

        t_314[k] = f_1 * smi0_195[k]
                   - f_2 * smi1_195[k]
                   + f_3 * pc_z[k] * smk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pb_z, pc_y, pc_z, sll0_135, sll0_138, \
                         slk_108, slk_144, sll1_135, sll1_138, \
                         smk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pb_z[k] * sll0_135[k]
                   - f_14 * pc_z[k] * sll1_135[k];

        t_316[k] = f_16 * slk_144[k]
                   + f_3 * pc_y[k] * smk_252[k];

        t_317[k] = f_15 * slk_108[k]
                   + f_3 * pc_z[k] * smk_252[k];

        t_318[k] = pb_z[k] * sll0_138[k]
                   - f_14 * pc_z[k] * sll1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pb_z, pc_x, pc_y, pc_z, sll0_141, slk_146, \
                         slk_257, sll1_141, smi0_201, smi1_201, smk_254, \
                         smk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * slk_146[k]
                   + f_3 * pc_y[k] * smk_254[k];

        t_320[k] = f_20 * slk_257[k]
                   + f_4 * smi0_201[k]
                   - f_5 * smi1_201[k]
                   + f_3 * pc_x[k] * smk_257[k];

        t_321[k] = pb_z[k] * sll0_141[k]
                   - f_14 * pc_z[k] * sll1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, slk_111, slk_149, slk_261, \
                         smi0_205, smi1_205, smk_255, smk_257, \
                         smk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * slk_111[k]
                   + f_3 * pc_z[k] * smk_255[k];

        t_323[k] = f_16 * slk_149[k]
                   + f_3 * pc_y[k] * smk_257[k];

        t_324[k] = f_20 * slk_261[k]
                   + f_6 * smi0_205[k]
                   - f_7 * smi1_205[k]
                   + f_3 * pc_x[k] * smk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pb_z, pc_y, pc_z, sll0_145, sll0_147, \
                         slk_114, slk_115, slk_153, sll1_145, sll1_147, smk_258, \
                         smk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pb_z[k] * sll0_145[k]
                   - f_14 * pc_z[k] * sll1_145[k];

        t_326[k] = f_15 * slk_114[k]
                   + f_3 * pc_z[k] * smk_258[k];

        t_327[k] = pb_z[k] * sll0_147[k]
                   + f_16 * slk_115[k]
                   - f_14 * pc_z[k] * sll1_147[k];

        t_328[k] = f_16 * slk_153[k]
                   + f_3 * pc_y[k] * smk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pb_z, pc_x, pc_z, sll0_150, slk_118, slk_266, \
                         sll1_150, smi0_210, smi1_210, smk_262, \
                         smk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_20 * slk_266[k]
                   + f_8 * smi0_210[k]
                   - f_9 * smi1_210[k]
                   + f_3 * pc_x[k] * smk_266[k];

        t_330[k] = pb_z[k] * sll0_150[k]
                   - f_14 * pc_z[k] * sll1_150[k];

        t_331[k] = f_15 * slk_118[k]
                   + f_3 * pc_z[k] * smk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_z, pc_y, pc_z, sll0_152, sll0_153, slk_119, \
                         slk_120, slk_158, sll1_152, sll1_153, \
                         smk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_z[k] * sll0_152[k]
                   + f_16 * slk_119[k]
                   - f_14 * pc_z[k] * sll1_152[k];

        t_333[k] = pb_z[k] * sll0_153[k]
                   + f_17 * slk_120[k]
                   - f_14 * pc_z[k] * sll1_153[k];

        t_334[k] = f_16 * slk_158[k]
                   + f_3 * pc_y[k] * smk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pb_z, pc_x, pc_z, sll0_156, slk_123, slk_272, \
                         sll1_156, smi0_216, smi1_216, smk_267, \
                         smk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_20 * slk_272[k]
                   + f_10 * smi0_216[k]
                   - f_11 * smi1_216[k]
                   + f_3 * pc_x[k] * smk_272[k];

        t_336[k] = pb_z[k] * sll0_156[k]
                   - f_14 * pc_z[k] * sll1_156[k];

        t_337[k] = f_15 * slk_123[k]
                   + f_3 * pc_z[k] * smk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_z, sll0_158, sll0_159, sll0_160, \
                         slk_124, slk_125, slk_126, sll1_158, sll1_159, \
                         sll1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pb_z[k] * sll0_158[k]
                   + f_16 * slk_124[k]
                   - f_14 * pc_z[k] * sll1_158[k];

        t_339[k] = pb_z[k] * sll0_159[k]
                   + f_17 * slk_125[k]
                   - f_14 * pc_z[k] * sll1_159[k];

        t_340[k] = pb_z[k] * sll0_160[k]
                   + f_18 * slk_126[k]
                   - f_14 * pc_z[k] * sll1_160[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_171 = buffer.data(sll0 + 171);
    const auto *sll0_225 = buffer.data(sll0 + 225);
    const auto *sll0_228 = buffer.data(sll0 + 228);
    const auto *sll0_230 = buffer.data(sll0 + 230);
    const auto *sll0_231 = buffer.data(sll0 + 231);
    const auto *sll0_234 = buffer.data(sll0 + 234);
    const auto *sll0_235 = buffer.data(sll0 + 235);
    const auto *sll0_237 = buffer.data(sll0 + 237);
    const auto *sll0_239 = buffer.data(sll0 + 239);
    const auto *sll0_240 = buffer.data(sll0 + 240);
    const auto *sll0_242 = buffer.data(sll0 + 242);
    const auto *sll0_243 = buffer.data(sll0 + 243);
    const auto *sll0_245 = buffer.data(sll0 + 245);
    const auto *sll0_246 = buffer.data(sll0 + 246);
    const auto *sll0_248 = buffer.data(sll0 + 248);
    const auto *sll0_249 = buffer.data(sll0 + 249);
    const auto *sll0_250 = buffer.data(sll0 + 250);
    const auto *sll0_252 = buffer.data(sll0 + 252);
    const auto *sll0_269 = buffer.data(sll0 + 269);

    const auto *slk_136 = buffer.data(slk + 136);
    const auto *slk_143 = buffer.data(slk + 143);
    const auto *slk_144 = buffer.data(slk + 144);
    const auto *slk_147 = buffer.data(slk + 147);
    const auto *slk_150 = buffer.data(slk + 150);
    const auto *slk_154 = buffer.data(slk + 154);
    const auto *slk_159 = buffer.data(slk + 159);
    const auto *slk_164 = buffer.data(slk + 164);
    const auto *slk_172 = buffer.data(slk + 172);
    const auto *slk_174 = buffer.data(slk + 174);
    const auto *slk_175 = buffer.data(slk + 175);
    const auto *slk_176 = buffer.data(slk + 176);
    const auto *slk_177 = buffer.data(slk + 177);
    const auto *slk_178 = buffer.data(slk + 178);
    const auto *slk_179 = buffer.data(slk + 179);
    const auto *slk_180 = buffer.data(slk + 180);
    const auto *slk_181 = buffer.data(slk + 181);
    const auto *slk_182 = buffer.data(slk + 182);
    const auto *slk_183 = buffer.data(slk + 183);
    const auto *slk_185 = buffer.data(slk + 185);
    const auto *slk_186 = buffer.data(slk + 186);
    const auto *slk_188 = buffer.data(slk + 188);
    const auto *slk_189 = buffer.data(slk + 189);
    const auto *slk_190 = buffer.data(slk + 190);
    const auto *slk_192 = buffer.data(slk + 192);
    const auto *slk_193 = buffer.data(slk + 193);
    const auto *slk_194 = buffer.data(slk + 194);
    const auto *slk_195 = buffer.data(slk + 195);
    const auto *slk_197 = buffer.data(slk + 197);
    const auto *slk_198 = buffer.data(slk + 198);
    const auto *slk_199 = buffer.data(slk + 199);
    const auto *slk_200 = buffer.data(slk + 200);
    const auto *slk_208 = buffer.data(slk + 208);
    const auto *slk_210 = buffer.data(slk + 210);
    const auto *slk_211 = buffer.data(slk + 211);
    const auto *slk_212 = buffer.data(slk + 212);
    const auto *slk_213 = buffer.data(slk + 213);
    const auto *slk_214 = buffer.data(slk + 214);
    const auto *slk_215 = buffer.data(slk + 215);
    const auto *slk_216 = buffer.data(slk + 216);
    const auto *slk_218 = buffer.data(slk + 218);
    const auto *slk_279 = buffer.data(slk + 279);
    const auto *slk_280 = buffer.data(slk + 280);
    const auto *slk_281 = buffer.data(slk + 281);
    const auto *slk_282 = buffer.data(slk + 282);
    const auto *slk_283 = buffer.data(slk + 283);
    const auto *slk_284 = buffer.data(slk + 284);
    const auto *slk_285 = buffer.data(slk + 285);
    const auto *slk_286 = buffer.data(slk + 286);
    const auto *slk_287 = buffer.data(slk + 287);
    const auto *slk_316 = buffer.data(slk + 316);
    const auto *slk_317 = buffer.data(slk + 317);
    const auto *slk_318 = buffer.data(slk + 318);
    const auto *slk_319 = buffer.data(slk + 319);
    const auto *slk_320 = buffer.data(slk + 320);
    const auto *slk_321 = buffer.data(slk + 321);
    const auto *slk_322 = buffer.data(slk + 322);
    const auto *slk_323 = buffer.data(slk + 323);
    const auto *slk_324 = buffer.data(slk + 324);
    const auto *slk_327 = buffer.data(slk + 327);
    const auto *slk_329 = buffer.data(slk + 329);
    const auto *slk_330 = buffer.data(slk + 330);
    const auto *slk_333 = buffer.data(slk + 333);
    const auto *slk_334 = buffer.data(slk + 334);
    const auto *slk_336 = buffer.data(slk + 336);
    const auto *slk_338 = buffer.data(slk + 338);
    const auto *slk_339 = buffer.data(slk + 339);
    const auto *slk_341 = buffer.data(slk + 341);
    const auto *slk_342 = buffer.data(slk + 342);
    const auto *slk_344 = buffer.data(slk + 344);
    const auto *slk_345 = buffer.data(slk + 345);
    const auto *slk_347 = buffer.data(slk + 347);
    const auto *slk_348 = buffer.data(slk + 348);
    const auto *slk_349 = buffer.data(slk + 349);
    const auto *slk_351 = buffer.data(slk + 351);
    const auto *slk_352 = buffer.data(slk + 352);
    const auto *slk_353 = buffer.data(slk + 353);
    const auto *slk_354 = buffer.data(slk + 354);
    const auto *slk_355 = buffer.data(slk + 355);
    const auto *slk_356 = buffer.data(slk + 356);
    const auto *slk_357 = buffer.data(slk + 357);
    const auto *slk_358 = buffer.data(slk + 358);
    const auto *slk_359 = buffer.data(slk + 359);
    const auto *slk_360 = buffer.data(slk + 360);
    const auto *slk_363 = buffer.data(slk + 363);
    const auto *slk_365 = buffer.data(slk + 365);

    const auto *sll1_171 = buffer.data(sll1 + 171);
    const auto *sll1_225 = buffer.data(sll1 + 225);
    const auto *sll1_228 = buffer.data(sll1 + 228);
    const auto *sll1_230 = buffer.data(sll1 + 230);
    const auto *sll1_231 = buffer.data(sll1 + 231);
    const auto *sll1_234 = buffer.data(sll1 + 234);
    const auto *sll1_235 = buffer.data(sll1 + 235);
    const auto *sll1_237 = buffer.data(sll1 + 237);
    const auto *sll1_239 = buffer.data(sll1 + 239);
    const auto *sll1_240 = buffer.data(sll1 + 240);
    const auto *sll1_242 = buffer.data(sll1 + 242);
    const auto *sll1_243 = buffer.data(sll1 + 243);
    const auto *sll1_245 = buffer.data(sll1 + 245);
    const auto *sll1_246 = buffer.data(sll1 + 246);
    const auto *sll1_248 = buffer.data(sll1 + 248);
    const auto *sll1_249 = buffer.data(sll1 + 249);
    const auto *sll1_250 = buffer.data(sll1 + 250);
    const auto *sll1_252 = buffer.data(sll1 + 252);
    const auto *sll1_269 = buffer.data(sll1 + 269);

    const auto *smi0_219 = buffer.data(smi0 + 219);
    const auto *smi0_220 = buffer.data(smi0 + 220);
    const auto *smi0_221 = buffer.data(smi0 + 221);
    const auto *smi0_222 = buffer.data(smi0 + 222);
    const auto *smi0_223 = buffer.data(smi0 + 223);
    const auto *smi0_245 = buffer.data(smi0 + 245);
    const auto *smi0_247 = buffer.data(smi0 + 247);
    const auto *smi0_248 = buffer.data(smi0 + 248);
    const auto *smi0_249 = buffer.data(smi0 + 249);
    const auto *smi0_250 = buffer.data(smi0 + 250);
    const auto *smi0_251 = buffer.data(smi0 + 251);
    const auto *smi0_252 = buffer.data(smi0 + 252);
    const auto *smi0_255 = buffer.data(smi0 + 255);
    const auto *smi0_257 = buffer.data(smi0 + 257);
    const auto *smi0_258 = buffer.data(smi0 + 258);
    const auto *smi0_261 = buffer.data(smi0 + 261);
    const auto *smi0_262 = buffer.data(smi0 + 262);
    const auto *smi0_264 = buffer.data(smi0 + 264);
    const auto *smi0_266 = buffer.data(smi0 + 266);
    const auto *smi0_267 = buffer.data(smi0 + 267);
    const auto *smi0_269 = buffer.data(smi0 + 269);
    const auto *smi0_270 = buffer.data(smi0 + 270);
    const auto *smi0_272 = buffer.data(smi0 + 272);
    const auto *smi0_273 = buffer.data(smi0 + 273);
    const auto *smi0_275 = buffer.data(smi0 + 275);
    const auto *smi0_276 = buffer.data(smi0 + 276);
    const auto *smi0_277 = buffer.data(smi0 + 277);
    const auto *smi0_278 = buffer.data(smi0 + 278);
    const auto *smi0_279 = buffer.data(smi0 + 279);
    const auto *smi0_280 = buffer.data(smi0 + 280);
    const auto *smi0_283 = buffer.data(smi0 + 283);
    const auto *smi0_285 = buffer.data(smi0 + 285);

    const auto *smi1_219 = buffer.data(smi1 + 219);
    const auto *smi1_220 = buffer.data(smi1 + 220);
    const auto *smi1_221 = buffer.data(smi1 + 221);
    const auto *smi1_222 = buffer.data(smi1 + 222);
    const auto *smi1_223 = buffer.data(smi1 + 223);
    const auto *smi1_245 = buffer.data(smi1 + 245);
    const auto *smi1_247 = buffer.data(smi1 + 247);
    const auto *smi1_248 = buffer.data(smi1 + 248);
    const auto *smi1_249 = buffer.data(smi1 + 249);
    const auto *smi1_250 = buffer.data(smi1 + 250);
    const auto *smi1_251 = buffer.data(smi1 + 251);
    const auto *smi1_252 = buffer.data(smi1 + 252);
    const auto *smi1_255 = buffer.data(smi1 + 255);
    const auto *smi1_257 = buffer.data(smi1 + 257);
    const auto *smi1_258 = buffer.data(smi1 + 258);
    const auto *smi1_261 = buffer.data(smi1 + 261);
    const auto *smi1_262 = buffer.data(smi1 + 262);
    const auto *smi1_264 = buffer.data(smi1 + 264);
    const auto *smi1_266 = buffer.data(smi1 + 266);
    const auto *smi1_267 = buffer.data(smi1 + 267);
    const auto *smi1_269 = buffer.data(smi1 + 269);
    const auto *smi1_270 = buffer.data(smi1 + 270);
    const auto *smi1_272 = buffer.data(smi1 + 272);
    const auto *smi1_273 = buffer.data(smi1 + 273);
    const auto *smi1_275 = buffer.data(smi1 + 275);
    const auto *smi1_276 = buffer.data(smi1 + 276);
    const auto *smi1_277 = buffer.data(smi1 + 277);
    const auto *smi1_278 = buffer.data(smi1 + 278);
    const auto *smi1_279 = buffer.data(smi1 + 279);
    const auto *smi1_280 = buffer.data(smi1 + 280);
    const auto *smi1_283 = buffer.data(smi1 + 283);
    const auto *smi1_285 = buffer.data(smi1 + 285);

    const auto *smk_272 = buffer.data(smk + 272);
    const auto *smk_279 = buffer.data(smk + 279);
    const auto *smk_280 = buffer.data(smk + 280);
    const auto *smk_281 = buffer.data(smk + 281);
    const auto *smk_282 = buffer.data(smk + 282);
    const auto *smk_283 = buffer.data(smk + 283);
    const auto *smk_284 = buffer.data(smk + 284);
    const auto *smk_285 = buffer.data(smk + 285);
    const auto *smk_286 = buffer.data(smk + 286);
    const auto *smk_287 = buffer.data(smk + 287);
    const auto *smk_288 = buffer.data(smk + 288);
    const auto *smk_290 = buffer.data(smk + 290);
    const auto *smk_291 = buffer.data(smk + 291);
    const auto *smk_293 = buffer.data(smk + 293);
    const auto *smk_294 = buffer.data(smk + 294);
    const auto *smk_297 = buffer.data(smk + 297);
    const auto *smk_298 = buffer.data(smk + 298);
    const auto *smk_302 = buffer.data(smk + 302);
    const auto *smk_303 = buffer.data(smk + 303);
    const auto *smk_308 = buffer.data(smk + 308);
    const auto *smk_316 = buffer.data(smk + 316);
    const auto *smk_317 = buffer.data(smk + 317);
    const auto *smk_318 = buffer.data(smk + 318);
    const auto *smk_319 = buffer.data(smk + 319);
    const auto *smk_320 = buffer.data(smk + 320);
    const auto *smk_321 = buffer.data(smk + 321);
    const auto *smk_322 = buffer.data(smk + 322);
    const auto *smk_323 = buffer.data(smk + 323);
    const auto *smk_324 = buffer.data(smk + 324);
    const auto *smk_326 = buffer.data(smk + 326);
    const auto *smk_327 = buffer.data(smk + 327);
    const auto *smk_329 = buffer.data(smk + 329);
    const auto *smk_330 = buffer.data(smk + 330);
    const auto *smk_333 = buffer.data(smk + 333);
    const auto *smk_334 = buffer.data(smk + 334);
    const auto *smk_336 = buffer.data(smk + 336);
    const auto *smk_338 = buffer.data(smk + 338);
    const auto *smk_339 = buffer.data(smk + 339);
    const auto *smk_341 = buffer.data(smk + 341);
    const auto *smk_342 = buffer.data(smk + 342);
    const auto *smk_344 = buffer.data(smk + 344);
    const auto *smk_345 = buffer.data(smk + 345);
    const auto *smk_347 = buffer.data(smk + 347);
    const auto *smk_348 = buffer.data(smk + 348);
    const auto *smk_349 = buffer.data(smk + 349);
    const auto *smk_351 = buffer.data(smk + 351);
    const auto *smk_352 = buffer.data(smk + 352);
    const auto *smk_353 = buffer.data(smk + 353);
    const auto *smk_354 = buffer.data(smk + 354);
    const auto *smk_355 = buffer.data(smk + 355);
    const auto *smk_356 = buffer.data(smk + 356);
    const auto *smk_357 = buffer.data(smk + 357);
    const auto *smk_358 = buffer.data(smk + 358);
    const auto *smk_359 = buffer.data(smk + 359);
    const auto *smk_360 = buffer.data(smk + 360);
    const auto *smk_362 = buffer.data(smk + 362);
    const auto *smk_363 = buffer.data(smk + 363);
    const auto *smk_365 = buffer.data(smk + 365);

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, slk_164, slk_279, slk_280, \
                         slk_281, smi0_223, smi1_223, smk_272, smk_279, smk_280, \
                         smk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * slk_164[k]
                   + f_3 * pc_y[k] * smk_272[k];

        t_342[k] = f_20 * slk_279[k]
                   + f_12 * smi0_223[k]
                   - f_13 * smi1_223[k]
                   + f_3 * pc_x[k] * smk_279[k];

        t_343[k] = f_20 * slk_280[k]
                   + f_3 * pc_x[k] * smk_280[k];

        t_344[k] = f_20 * slk_281[k]
                   + f_3 * pc_x[k] * smk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, slk_282, slk_283, slk_284, \
                         slk_285, slk_286, smk_282, smk_283, smk_284, smk_285, \
                         smk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_20 * slk_282[k]
                   + f_3 * pc_x[k] * smk_282[k];

        t_346[k] = f_20 * slk_283[k]
                   + f_3 * pc_x[k] * smk_283[k];

        t_347[k] = f_20 * slk_284[k]
                   + f_3 * pc_x[k] * smk_284[k];

        t_348[k] = f_20 * slk_285[k]
                   + f_3 * pc_x[k] * smk_285[k];

        t_349[k] = f_20 * slk_286[k]
                   + f_3 * pc_x[k] * smk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pb_z, pc_x, pc_z, sll0_171, slk_136, slk_287, \
                         sll1_171, smk_280, smk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_20 * slk_287[k]
                   + f_3 * pc_x[k] * smk_287[k];

        t_351[k] = pb_z[k] * sll0_171[k]
                   - f_14 * pc_z[k] * sll1_171[k];

        t_352[k] = f_15 * slk_136[k]
                   + f_3 * pc_z[k] * smk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, slk_174, slk_175, slk_176, smi0_219, \
                         smi0_220, smi0_221, smi1_219, smi1_220, smi1_221, smk_282, smk_283, \
                         smk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * slk_174[k]
                   + f_4 * smi0_219[k]
                   - f_5 * smi1_219[k]
                   + f_3 * pc_y[k] * smk_282[k];

        t_354[k] = f_16 * slk_175[k]
                   + f_6 * smi0_220[k]
                   - f_7 * smi1_220[k]
                   + f_3 * pc_y[k] * smk_283[k];

        t_355[k] = f_16 * slk_176[k]
                   + f_8 * smi0_221[k]
                   - f_9 * smi1_221[k]
                   + f_3 * pc_y[k] * smk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, slk_177, slk_178, slk_179, smi0_222, \
                         smi0_223, smi1_222, smi1_223, smk_285, smk_286, \
                         smk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * slk_177[k]
                   + f_10 * smi0_222[k]
                   - f_11 * smi1_222[k]
                   + f_3 * pc_y[k] * smk_285[k];

        t_357[k] = f_16 * slk_178[k]
                   + f_12 * smi0_223[k]
                   - f_13 * smi1_223[k]
                   + f_3 * pc_y[k] * smk_286[k];

        t_358[k] = f_16 * slk_179[k]
                   + f_3 * pc_y[k] * smk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pb_y, pc_y, pc_z, sll0_225, slk_143, \
                         slk_144, slk_180, sll1_225, smi0_223, smi1_223, smk_287, \
                         smk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * slk_143[k]
                   + f_1 * smi0_223[k]
                   - f_2 * smi1_223[k]
                   + f_3 * pc_z[k] * smk_287[k];

        t_360[k] = pb_y[k] * sll0_225[k]
                   - f_14 * pc_y[k] * sll1_225[k];

        t_361[k] = f_15 * slk_180[k]
                   + f_3 * pc_y[k] * smk_288[k];

        t_362[k] = f_16 * slk_144[k]
                   + f_3 * pc_z[k] * smk_288[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pb_y, pc_y, sll0_228, sll0_230, sll0_231, \
                         slk_181, slk_182, slk_183, sll1_228, sll1_230, sll1_231, \
                         smk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pb_y[k] * sll0_228[k]
                   + f_16 * slk_181[k]
                   - f_14 * pc_y[k] * sll1_228[k];

        t_364[k] = f_15 * slk_182[k]
                   + f_3 * pc_y[k] * smk_290[k];

        t_365[k] = pb_y[k] * sll0_230[k]
                   - f_14 * pc_y[k] * sll1_230[k];

        t_366[k] = pb_y[k] * sll0_231[k]
                   + f_17 * slk_183[k]
                   - f_14 * pc_y[k] * sll1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pb_y, pc_y, pc_z, sll0_234, sll0_235, \
                         slk_147, slk_185, slk_186, sll1_234, sll1_235, smk_291, \
                         smk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * slk_147[k]
                   + f_3 * pc_z[k] * smk_291[k];

        t_368[k] = f_15 * slk_185[k]
                   + f_3 * pc_y[k] * smk_293[k];

        t_369[k] = pb_y[k] * sll0_234[k]
                   - f_14 * pc_y[k] * sll1_234[k];

        t_370[k] = pb_y[k] * sll0_235[k]
                   + f_18 * slk_186[k]
                   - f_14 * pc_y[k] * sll1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pb_y, pc_y, pc_z, sll0_237, sll0_239, \
                         slk_150, slk_188, slk_189, sll1_237, sll1_239, smk_294, \
                         smk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * slk_150[k]
                   + f_3 * pc_z[k] * smk_294[k];

        t_372[k] = pb_y[k] * sll0_237[k]
                   + f_16 * slk_188[k]
                   - f_14 * pc_y[k] * sll1_237[k];

        t_373[k] = f_15 * slk_189[k]
                   + f_3 * pc_y[k] * smk_297[k];

        t_374[k] = pb_y[k] * sll0_239[k]
                   - f_14 * pc_y[k] * sll1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pb_y, pc_y, pc_z, sll0_240, sll0_242, slk_154, \
                         slk_190, slk_192, sll1_240, sll1_242, \
                         smk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pb_y[k] * sll0_240[k]
                   + f_19 * slk_190[k]
                   - f_14 * pc_y[k] * sll1_240[k];

        t_376[k] = f_16 * slk_154[k]
                   + f_3 * pc_z[k] * smk_298[k];

        t_377[k] = pb_y[k] * sll0_242[k]
                   + f_17 * slk_192[k]
                   - f_14 * pc_y[k] * sll1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pb_y, pc_y, sll0_243, sll0_245, sll0_246, \
                         slk_193, slk_194, slk_195, sll1_243, sll1_245, sll1_246, \
                         smk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_y[k] * sll0_243[k]
                   + f_16 * slk_193[k]
                   - f_14 * pc_y[k] * sll1_243[k];

        t_379[k] = f_15 * slk_194[k]
                   + f_3 * pc_y[k] * smk_302[k];

        t_380[k] = pb_y[k] * sll0_245[k]
                   - f_14 * pc_y[k] * sll1_245[k];

        t_381[k] = pb_y[k] * sll0_246[k]
                   + f_20 * slk_195[k]
                   - f_14 * pc_y[k] * sll1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pb_y, pc_y, pc_z, sll0_248, sll0_249, slk_159, \
                         slk_197, slk_198, sll1_248, sll1_249, \
                         smk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * slk_159[k]
                   + f_3 * pc_z[k] * smk_303[k];

        t_383[k] = pb_y[k] * sll0_248[k]
                   + f_18 * slk_197[k]
                   - f_14 * pc_y[k] * sll1_248[k];

        t_384[k] = pb_y[k] * sll0_249[k]
                   + f_17 * slk_198[k]
                   - f_14 * pc_y[k] * sll1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pb_y, pc_x, pc_y, sll0_250, sll0_252, \
                         slk_199, slk_200, slk_316, sll1_250, sll1_252, smk_308, \
                         smk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pb_y[k] * sll0_250[k]
                   + f_16 * slk_199[k]
                   - f_14 * pc_y[k] * sll1_250[k];

        t_386[k] = f_15 * slk_200[k]
                   + f_3 * pc_y[k] * smk_308[k];

        t_387[k] = pb_y[k] * sll0_252[k]
                   - f_14 * pc_y[k] * sll1_252[k];

        t_388[k] = f_20 * slk_316[k]
                   + f_3 * pc_x[k] * smk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, slk_317, slk_318, slk_319, \
                         slk_320, slk_321, smk_317, smk_318, smk_319, smk_320, \
                         smk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_20 * slk_317[k]
                   + f_3 * pc_x[k] * smk_317[k];

        t_390[k] = f_20 * slk_318[k]
                   + f_3 * pc_x[k] * smk_318[k];

        t_391[k] = f_20 * slk_319[k]
                   + f_3 * pc_x[k] * smk_319[k];

        t_392[k] = f_20 * slk_320[k]
                   + f_3 * pc_x[k] * smk_320[k];

        t_393[k] = f_20 * slk_321[k]
                   + f_3 * pc_x[k] * smk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, slk_172, slk_208, \
                         slk_322, slk_323, smi0_245, smi1_245, smk_316, smk_322, \
                         smk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_20 * slk_322[k]
                   + f_3 * pc_x[k] * smk_322[k];

        t_395[k] = f_20 * slk_323[k]
                   + f_3 * pc_x[k] * smk_323[k];

        t_396[k] = f_15 * slk_208[k]
                   + f_1 * smi0_245[k]
                   - f_2 * smi1_245[k]
                   + f_3 * pc_y[k] * smk_316[k];

        t_397[k] = f_16 * slk_172[k]
                   + f_3 * pc_z[k] * smk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, slk_210, slk_211, slk_212, smi0_247, \
                         smi0_248, smi0_249, smi1_247, smi1_248, smi1_249, smk_318, smk_319, \
                         smk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * slk_210[k]
                   + f_4 * smi0_247[k]
                   - f_5 * smi1_247[k]
                   + f_3 * pc_y[k] * smk_318[k];

        t_399[k] = f_15 * slk_211[k]
                   + f_6 * smi0_248[k]
                   - f_7 * smi1_248[k]
                   + f_3 * pc_y[k] * smk_319[k];

        t_400[k] = f_15 * slk_212[k]
                   + f_8 * smi0_249[k]
                   - f_9 * smi1_249[k]
                   + f_3 * pc_y[k] * smk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, slk_213, slk_214, slk_215, smi0_250, \
                         smi0_251, smi1_250, smi1_251, smk_321, smk_322, \
                         smk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * slk_213[k]
                   + f_10 * smi0_250[k]
                   - f_11 * smi1_250[k]
                   + f_3 * pc_y[k] * smk_321[k];

        t_402[k] = f_15 * slk_214[k]
                   + f_12 * smi0_251[k]
                   - f_13 * smi1_251[k]
                   + f_3 * pc_y[k] * smk_322[k];

        t_403[k] = f_15 * slk_215[k]
                   + f_3 * pc_y[k] * smk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pb_y, pc_x, pc_y, pc_z, sll0_269, \
                         slk_180, slk_324, sll1_269, smi0_252, smi1_252, \
                         smk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pb_y[k] * sll0_269[k]
                   - f_14 * pc_y[k] * sll1_269[k];

        t_405[k] = f_20 * slk_324[k]
                   + f_1 * smi0_252[k]
                   - f_2 * smi1_252[k]
                   + f_3 * pc_x[k] * smk_324[k];

        t_406[k] = f_3 * pc_y[k] * smk_324[k];

        t_407[k] = f_17 * slk_180[k]
                   + f_3 * pc_z[k] * smk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, slk_327, slk_329, smi0_255, \
                         smi0_257, smi1_255, smi1_257, smk_326, smk_327, \
                         smk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_20 * slk_327[k]
                   + f_4 * smi0_255[k]
                   - f_5 * smi1_255[k]
                   + f_3 * pc_x[k] * smk_327[k];

        t_409[k] = f_3 * pc_y[k] * smk_326[k];

        t_410[k] = f_20 * slk_329[k]
                   + f_4 * smi0_257[k]
                   - f_5 * smi1_257[k]
                   + f_3 * pc_x[k] * smk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_x, pc_y, pc_z, slk_183, slk_330, smi0_258, \
                         smi1_258, smk_327, smk_329, smk_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_20 * slk_330[k]
                   + f_6 * smi0_258[k]
                   - f_7 * smi1_258[k]
                   + f_3 * pc_x[k] * smk_330[k];

        t_412[k] = f_17 * slk_183[k]
                   + f_3 * pc_z[k] * smk_327[k];

        t_413[k] = f_3 * pc_y[k] * smk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_z, slk_186, slk_333, slk_334, smi0_261, \
                         smi0_262, smi1_261, smi1_262, smk_330, smk_333, \
                         smk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_20 * slk_333[k]
                   + f_6 * smi0_261[k]
                   - f_7 * smi1_261[k]
                   + f_3 * pc_x[k] * smk_333[k];

        t_415[k] = f_20 * slk_334[k]
                   + f_8 * smi0_262[k]
                   - f_9 * smi1_262[k]
                   + f_3 * pc_x[k] * smk_334[k];

        t_416[k] = f_17 * slk_186[k]
                   + f_3 * pc_z[k] * smk_330[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, slk_336, slk_338, smi0_264, \
                         smi0_266, smi1_264, smi1_266, smk_333, smk_336, \
                         smk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_20 * slk_336[k]
                   + f_8 * smi0_264[k]
                   - f_9 * smi1_264[k]
                   + f_3 * pc_x[k] * smk_336[k];

        t_418[k] = f_3 * pc_y[k] * smk_333[k];

        t_419[k] = f_20 * slk_338[k]
                   + f_8 * smi0_266[k]
                   - f_9 * smi1_266[k]
                   + f_3 * pc_x[k] * smk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_x, pc_z, slk_190, slk_339, slk_341, smi0_267, \
                         smi0_269, smi1_267, smi1_269, smk_334, smk_339, \
                         smk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_20 * slk_339[k]
                   + f_10 * smi0_267[k]
                   - f_11 * smi1_267[k]
                   + f_3 * pc_x[k] * smk_339[k];

        t_421[k] = f_17 * slk_190[k]
                   + f_3 * pc_z[k] * smk_334[k];

        t_422[k] = f_20 * slk_341[k]
                   + f_10 * smi0_269[k]
                   - f_11 * smi1_269[k]
                   + f_3 * pc_x[k] * smk_341[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, slk_342, slk_344, smi0_270, \
                         smi0_272, smi1_270, smi1_272, smk_338, smk_342, \
                         smk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_20 * slk_342[k]
                   + f_10 * smi0_270[k]
                   - f_11 * smi1_270[k]
                   + f_3 * pc_x[k] * smk_342[k];

        t_424[k] = f_3 * pc_y[k] * smk_338[k];

        t_425[k] = f_20 * slk_344[k]
                   + f_10 * smi0_272[k]
                   - f_11 * smi1_272[k]
                   + f_3 * pc_x[k] * smk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_x, pc_z, slk_195, slk_345, slk_347, smi0_273, \
                         smi0_275, smi1_273, smi1_275, smk_339, smk_345, \
                         smk_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_20 * slk_345[k]
                   + f_12 * smi0_273[k]
                   - f_13 * smi1_273[k]
                   + f_3 * pc_x[k] * smk_345[k];

        t_427[k] = f_17 * slk_195[k]
                   + f_3 * pc_z[k] * smk_339[k];

        t_428[k] = f_20 * slk_347[k]
                   + f_12 * smi0_275[k]
                   - f_13 * smi1_275[k]
                   + f_3 * pc_x[k] * smk_347[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_x, pc_y, slk_348, slk_349, smi0_276, \
                         smi0_277, smi1_276, smi1_277, smk_344, smk_348, \
                         smk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_20 * slk_348[k]
                   + f_12 * smi0_276[k]
                   - f_13 * smi1_276[k]
                   + f_3 * pc_x[k] * smk_348[k];

        t_430[k] = f_20 * slk_349[k]
                   + f_12 * smi0_277[k]
                   - f_13 * smi1_277[k]
                   + f_3 * pc_x[k] * smk_349[k];

        t_431[k] = f_3 * pc_y[k] * smk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, slk_351, slk_352, slk_353, slk_354, \
                         smi0_279, smi1_279, smk_351, smk_352, smk_353, \
                         smk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_20 * slk_351[k]
                   + f_12 * smi0_279[k]
                   - f_13 * smi1_279[k]
                   + f_3 * pc_x[k] * smk_351[k];

        t_433[k] = f_20 * slk_352[k]
                   + f_3 * pc_x[k] * smk_352[k];

        t_434[k] = f_20 * slk_353[k]
                   + f_3 * pc_x[k] * smk_353[k];

        t_435[k] = f_20 * slk_354[k]
                   + f_3 * pc_x[k] * smk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, slk_355, slk_356, slk_357, \
                         slk_358, slk_359, smk_355, smk_356, smk_357, smk_358, \
                         smk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_20 * slk_355[k]
                   + f_3 * pc_x[k] * smk_355[k];

        t_437[k] = f_20 * slk_356[k]
                   + f_3 * pc_x[k] * smk_356[k];

        t_438[k] = f_20 * slk_357[k]
                   + f_3 * pc_x[k] * smk_357[k];

        t_439[k] = f_20 * slk_358[k]
                   + f_3 * pc_x[k] * smk_358[k];

        t_440[k] = f_20 * slk_359[k]
                   + f_3 * pc_x[k] * smk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, slk_208, smi0_273, smi0_275, \
                         smi0_276, smi1_273, smi1_275, smi1_276, smk_352, smk_354, \
                         smk_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * smi0_273[k]
                   - f_2 * smi1_273[k]
                   + f_3 * pc_y[k] * smk_352[k];

        t_442[k] = f_17 * slk_208[k]
                   + f_3 * pc_z[k] * smk_352[k];

        t_443[k] = f_4 * smi0_275[k]
                   - f_5 * smi1_275[k]
                   + f_3 * pc_y[k] * smk_354[k];

        t_444[k] = f_6 * smi0_276[k]
                   - f_7 * smi1_276[k]
                   + f_3 * pc_y[k] * smk_355[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_y, smi0_277, smi0_278, smi0_279, \
                         smi1_277, smi1_278, smi1_279, smk_356, smk_357, smk_358, \
                         smk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * smi0_277[k]
                   - f_9 * smi1_277[k]
                   + f_3 * pc_y[k] * smk_356[k];

        t_446[k] = f_10 * smi0_278[k]
                   - f_11 * smi1_278[k]
                   + f_3 * pc_y[k] * smk_357[k];

        t_447[k] = f_12 * smi0_279[k]
                   - f_13 * smi1_279[k]
                   + f_3 * pc_y[k] * smk_358[k];

        t_448[k] = f_3 * pc_y[k] * smk_359[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, slk_215, slk_216, \
                         slk_360, smi0_279, smi0_280, smi1_279, smi1_280, smk_359, \
                         smk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * slk_215[k]
                   + f_1 * smi0_279[k]
                   - f_2 * smi1_279[k]
                   + f_3 * pc_z[k] * smk_359[k];

        t_450[k] = f_19 * slk_360[k]
                   + f_1 * smi0_280[k]
                   - f_2 * smi1_280[k]
                   + f_3 * pc_x[k] * smk_360[k];

        t_451[k] = f_18 * slk_216[k]
                   + f_3 * pc_y[k] * smk_360[k];

        t_452[k] = f_3 * pc_z[k] * smk_360[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pc_x, pc_y, slk_218, slk_363, slk_365, smi0_283, \
                         smi0_285, smi1_283, smi1_285, smk_362, smk_363, \
                         smk_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_19 * slk_363[k]
                   + f_4 * smi0_283[k]
                   - f_5 * smi1_283[k]
                   + f_3 * pc_x[k] * smk_363[k];

        t_454[k] = f_18 * slk_218[k]
                   + f_3 * pc_y[k] * smk_362[k];

        t_455[k] = f_19 * slk_365[k]
                   + f_4 * smi0_285[k]
                   - f_5 * smi1_285[k]
                   + f_3 * pc_x[k] * smk_365[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_270 = buffer.data(sll0 + 270);
    const auto *sll0_273 = buffer.data(sll0 + 273);
    const auto *sll0_276 = buffer.data(sll0 + 276);
    const auto *sll0_280 = buffer.data(sll0 + 280);
    const auto *sll0_282 = buffer.data(sll0 + 282);
    const auto *sll0_285 = buffer.data(sll0 + 285);
    const auto *sll0_287 = buffer.data(sll0 + 287);
    const auto *sll0_288 = buffer.data(sll0 + 288);
    const auto *sll0_291 = buffer.data(sll0 + 291);
    const auto *sll0_293 = buffer.data(sll0 + 293);
    const auto *sll0_294 = buffer.data(sll0 + 294);
    const auto *sll0_295 = buffer.data(sll0 + 295);
    const auto *sll0_306 = buffer.data(sll0 + 306);

    const auto *slk_216 = buffer.data(slk + 216);
    const auto *slk_219 = buffer.data(slk + 219);
    const auto *slk_221 = buffer.data(slk + 221);
    const auto *slk_222 = buffer.data(slk + 222);
    const auto *slk_223 = buffer.data(slk + 223);
    const auto *slk_225 = buffer.data(slk + 225);
    const auto *slk_226 = buffer.data(slk + 226);
    const auto *slk_227 = buffer.data(slk + 227);
    const auto *slk_228 = buffer.data(slk + 228);
    const auto *slk_230 = buffer.data(slk + 230);
    const auto *slk_231 = buffer.data(slk + 231);
    const auto *slk_232 = buffer.data(slk + 232);
    const auto *slk_233 = buffer.data(slk + 233);
    const auto *slk_234 = buffer.data(slk + 234);
    const auto *slk_236 = buffer.data(slk + 236);
    const auto *slk_244 = buffer.data(slk + 244);
    const auto *slk_246 = buffer.data(slk + 246);
    const auto *slk_247 = buffer.data(slk + 247);
    const auto *slk_248 = buffer.data(slk + 248);
    const auto *slk_249 = buffer.data(slk + 249);
    const auto *slk_250 = buffer.data(slk + 250);
    const auto *slk_251 = buffer.data(slk + 251);
    const auto *slk_252 = buffer.data(slk + 252);
    const auto *slk_254 = buffer.data(slk + 254);
    const auto *slk_255 = buffer.data(slk + 255);
    const auto *slk_257 = buffer.data(slk + 257);
    const auto *slk_258 = buffer.data(slk + 258);
    const auto *slk_261 = buffer.data(slk + 261);
    const auto *slk_262 = buffer.data(slk + 262);
    const auto *slk_266 = buffer.data(slk + 266);
    const auto *slk_267 = buffer.data(slk + 267);
    const auto *slk_272 = buffer.data(slk + 272);
    const auto *slk_282 = buffer.data(slk + 282);
    const auto *slk_283 = buffer.data(slk + 283);
    const auto *slk_284 = buffer.data(slk + 284);
    const auto *slk_285 = buffer.data(slk + 285);
    const auto *slk_286 = buffer.data(slk + 286);
    const auto *slk_287 = buffer.data(slk + 287);
    const auto *slk_288 = buffer.data(slk + 288);
    const auto *slk_290 = buffer.data(slk + 290);
    const auto *slk_293 = buffer.data(slk + 293);
    const auto *slk_297 = buffer.data(slk + 297);
    const auto *slk_302 = buffer.data(slk + 302);
    const auto *slk_366 = buffer.data(slk + 366);
    const auto *slk_369 = buffer.data(slk + 369);
    const auto *slk_370 = buffer.data(slk + 370);
    const auto *slk_372 = buffer.data(slk + 372);
    const auto *slk_374 = buffer.data(slk + 374);
    const auto *slk_375 = buffer.data(slk + 375);
    const auto *slk_377 = buffer.data(slk + 377);
    const auto *slk_378 = buffer.data(slk + 378);
    const auto *slk_380 = buffer.data(slk + 380);
    const auto *slk_381 = buffer.data(slk + 381);
    const auto *slk_383 = buffer.data(slk + 383);
    const auto *slk_384 = buffer.data(slk + 384);
    const auto *slk_385 = buffer.data(slk + 385);
    const auto *slk_387 = buffer.data(slk + 387);
    const auto *slk_388 = buffer.data(slk + 388);
    const auto *slk_389 = buffer.data(slk + 389);
    const auto *slk_390 = buffer.data(slk + 390);
    const auto *slk_391 = buffer.data(slk + 391);
    const auto *slk_392 = buffer.data(slk + 392);
    const auto *slk_393 = buffer.data(slk + 393);
    const auto *slk_394 = buffer.data(slk + 394);
    const auto *slk_395 = buffer.data(slk + 395);
    const auto *slk_401 = buffer.data(slk + 401);
    const auto *slk_405 = buffer.data(slk + 405);
    const auto *slk_410 = buffer.data(slk + 410);
    const auto *slk_416 = buffer.data(slk + 416);
    const auto *slk_423 = buffer.data(slk + 423);
    const auto *slk_424 = buffer.data(slk + 424);
    const auto *slk_425 = buffer.data(slk + 425);
    const auto *slk_426 = buffer.data(slk + 426);
    const auto *slk_427 = buffer.data(slk + 427);
    const auto *slk_428 = buffer.data(slk + 428);
    const auto *slk_429 = buffer.data(slk + 429);
    const auto *slk_430 = buffer.data(slk + 430);
    const auto *slk_431 = buffer.data(slk + 431);
    const auto *slk_432 = buffer.data(slk + 432);
    const auto *slk_435 = buffer.data(slk + 435);
    const auto *slk_437 = buffer.data(slk + 437);
    const auto *slk_438 = buffer.data(slk + 438);
    const auto *slk_441 = buffer.data(slk + 441);
    const auto *slk_442 = buffer.data(slk + 442);
    const auto *slk_444 = buffer.data(slk + 444);
    const auto *slk_446 = buffer.data(slk + 446);
    const auto *slk_447 = buffer.data(slk + 447);
    const auto *slk_449 = buffer.data(slk + 449);
    const auto *slk_450 = buffer.data(slk + 450);
    const auto *slk_452 = buffer.data(slk + 452);
    const auto *slk_453 = buffer.data(slk + 453);

    const auto *sll1_270 = buffer.data(sll1 + 270);
    const auto *sll1_273 = buffer.data(sll1 + 273);
    const auto *sll1_276 = buffer.data(sll1 + 276);
    const auto *sll1_280 = buffer.data(sll1 + 280);
    const auto *sll1_282 = buffer.data(sll1 + 282);
    const auto *sll1_285 = buffer.data(sll1 + 285);
    const auto *sll1_287 = buffer.data(sll1 + 287);
    const auto *sll1_288 = buffer.data(sll1 + 288);
    const auto *sll1_291 = buffer.data(sll1 + 291);
    const auto *sll1_293 = buffer.data(sll1 + 293);
    const auto *sll1_294 = buffer.data(sll1 + 294);
    const auto *sll1_295 = buffer.data(sll1 + 295);
    const auto *sll1_306 = buffer.data(sll1 + 306);

    const auto *smi0_286 = buffer.data(smi0 + 286);
    const auto *smi0_289 = buffer.data(smi0 + 289);
    const auto *smi0_290 = buffer.data(smi0 + 290);
    const auto *smi0_292 = buffer.data(smi0 + 292);
    const auto *smi0_294 = buffer.data(smi0 + 294);
    const auto *smi0_295 = buffer.data(smi0 + 295);
    const auto *smi0_297 = buffer.data(smi0 + 297);
    const auto *smi0_298 = buffer.data(smi0 + 298);
    const auto *smi0_300 = buffer.data(smi0 + 300);
    const auto *smi0_301 = buffer.data(smi0 + 301);
    const auto *smi0_303 = buffer.data(smi0 + 303);
    const auto *smi0_304 = buffer.data(smi0 + 304);
    const auto *smi0_305 = buffer.data(smi0 + 305);
    const auto *smi0_306 = buffer.data(smi0 + 306);
    const auto *smi0_307 = buffer.data(smi0 + 307);
    const auto *smi0_313 = buffer.data(smi0 + 313);
    const auto *smi0_317 = buffer.data(smi0 + 317);
    const auto *smi0_322 = buffer.data(smi0 + 322);
    const auto *smi0_328 = buffer.data(smi0 + 328);
    const auto *smi0_331 = buffer.data(smi0 + 331);
    const auto *smi0_332 = buffer.data(smi0 + 332);
    const auto *smi0_333 = buffer.data(smi0 + 333);
    const auto *smi0_334 = buffer.data(smi0 + 334);
    const auto *smi0_335 = buffer.data(smi0 + 335);
    const auto *smi0_336 = buffer.data(smi0 + 336);
    const auto *smi0_339 = buffer.data(smi0 + 339);
    const auto *smi0_341 = buffer.data(smi0 + 341);
    const auto *smi0_342 = buffer.data(smi0 + 342);
    const auto *smi0_345 = buffer.data(smi0 + 345);
    const auto *smi0_346 = buffer.data(smi0 + 346);
    const auto *smi0_348 = buffer.data(smi0 + 348);
    const auto *smi0_350 = buffer.data(smi0 + 350);
    const auto *smi0_351 = buffer.data(smi0 + 351);
    const auto *smi0_353 = buffer.data(smi0 + 353);
    const auto *smi0_354 = buffer.data(smi0 + 354);
    const auto *smi0_356 = buffer.data(smi0 + 356);
    const auto *smi0_357 = buffer.data(smi0 + 357);

    const auto *smi1_286 = buffer.data(smi1 + 286);
    const auto *smi1_289 = buffer.data(smi1 + 289);
    const auto *smi1_290 = buffer.data(smi1 + 290);
    const auto *smi1_292 = buffer.data(smi1 + 292);
    const auto *smi1_294 = buffer.data(smi1 + 294);
    const auto *smi1_295 = buffer.data(smi1 + 295);
    const auto *smi1_297 = buffer.data(smi1 + 297);
    const auto *smi1_298 = buffer.data(smi1 + 298);
    const auto *smi1_300 = buffer.data(smi1 + 300);
    const auto *smi1_301 = buffer.data(smi1 + 301);
    const auto *smi1_303 = buffer.data(smi1 + 303);
    const auto *smi1_304 = buffer.data(smi1 + 304);
    const auto *smi1_305 = buffer.data(smi1 + 305);
    const auto *smi1_306 = buffer.data(smi1 + 306);
    const auto *smi1_307 = buffer.data(smi1 + 307);
    const auto *smi1_313 = buffer.data(smi1 + 313);
    const auto *smi1_317 = buffer.data(smi1 + 317);
    const auto *smi1_322 = buffer.data(smi1 + 322);
    const auto *smi1_328 = buffer.data(smi1 + 328);
    const auto *smi1_331 = buffer.data(smi1 + 331);
    const auto *smi1_332 = buffer.data(smi1 + 332);
    const auto *smi1_333 = buffer.data(smi1 + 333);
    const auto *smi1_334 = buffer.data(smi1 + 334);
    const auto *smi1_335 = buffer.data(smi1 + 335);
    const auto *smi1_336 = buffer.data(smi1 + 336);
    const auto *smi1_339 = buffer.data(smi1 + 339);
    const auto *smi1_341 = buffer.data(smi1 + 341);
    const auto *smi1_342 = buffer.data(smi1 + 342);
    const auto *smi1_345 = buffer.data(smi1 + 345);
    const auto *smi1_346 = buffer.data(smi1 + 346);
    const auto *smi1_348 = buffer.data(smi1 + 348);
    const auto *smi1_350 = buffer.data(smi1 + 350);
    const auto *smi1_351 = buffer.data(smi1 + 351);
    const auto *smi1_353 = buffer.data(smi1 + 353);
    const auto *smi1_354 = buffer.data(smi1 + 354);
    const auto *smi1_356 = buffer.data(smi1 + 356);
    const auto *smi1_357 = buffer.data(smi1 + 357);

    const auto *smk_363 = buffer.data(smk + 363);
    const auto *smk_365 = buffer.data(smk + 365);
    const auto *smk_366 = buffer.data(smk + 366);
    const auto *smk_369 = buffer.data(smk + 369);
    const auto *smk_370 = buffer.data(smk + 370);
    const auto *smk_372 = buffer.data(smk + 372);
    const auto *smk_374 = buffer.data(smk + 374);
    const auto *smk_375 = buffer.data(smk + 375);
    const auto *smk_377 = buffer.data(smk + 377);
    const auto *smk_378 = buffer.data(smk + 378);
    const auto *smk_380 = buffer.data(smk + 380);
    const auto *smk_381 = buffer.data(smk + 381);
    const auto *smk_383 = buffer.data(smk + 383);
    const auto *smk_384 = buffer.data(smk + 384);
    const auto *smk_385 = buffer.data(smk + 385);
    const auto *smk_387 = buffer.data(smk + 387);
    const auto *smk_388 = buffer.data(smk + 388);
    const auto *smk_389 = buffer.data(smk + 389);
    const auto *smk_390 = buffer.data(smk + 390);
    const auto *smk_391 = buffer.data(smk + 391);
    const auto *smk_392 = buffer.data(smk + 392);
    const auto *smk_393 = buffer.data(smk + 393);
    const auto *smk_394 = buffer.data(smk + 394);
    const auto *smk_395 = buffer.data(smk + 395);
    const auto *smk_396 = buffer.data(smk + 396);
    const auto *smk_398 = buffer.data(smk + 398);
    const auto *smk_399 = buffer.data(smk + 399);
    const auto *smk_401 = buffer.data(smk + 401);
    const auto *smk_402 = buffer.data(smk + 402);
    const auto *smk_405 = buffer.data(smk + 405);
    const auto *smk_406 = buffer.data(smk + 406);
    const auto *smk_410 = buffer.data(smk + 410);
    const auto *smk_411 = buffer.data(smk + 411);
    const auto *smk_416 = buffer.data(smk + 416);
    const auto *smk_423 = buffer.data(smk + 423);
    const auto *smk_424 = buffer.data(smk + 424);
    const auto *smk_425 = buffer.data(smk + 425);
    const auto *smk_426 = buffer.data(smk + 426);
    const auto *smk_427 = buffer.data(smk + 427);
    const auto *smk_428 = buffer.data(smk + 428);
    const auto *smk_429 = buffer.data(smk + 429);
    const auto *smk_430 = buffer.data(smk + 430);
    const auto *smk_431 = buffer.data(smk + 431);
    const auto *smk_432 = buffer.data(smk + 432);
    const auto *smk_434 = buffer.data(smk + 434);
    const auto *smk_435 = buffer.data(smk + 435);
    const auto *smk_437 = buffer.data(smk + 437);
    const auto *smk_438 = buffer.data(smk + 438);
    const auto *smk_441 = buffer.data(smk + 441);
    const auto *smk_442 = buffer.data(smk + 442);
    const auto *smk_444 = buffer.data(smk + 444);
    const auto *smk_446 = buffer.data(smk + 446);
    const auto *smk_447 = buffer.data(smk + 447);
    const auto *smk_449 = buffer.data(smk + 449);
    const auto *smk_450 = buffer.data(smk + 450);
    const auto *smk_452 = buffer.data(smk + 452);
    const auto *smk_453 = buffer.data(smk + 453);

#pragma omp simd aligned(t_456, t_457, t_458, pc_x, pc_y, pc_z, slk_221, slk_366, smi0_286, \
                         smi1_286, smk_363, smk_365, smk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_19 * slk_366[k]
                   + f_6 * smi0_286[k]
                   - f_7 * smi1_286[k]
                   + f_3 * pc_x[k] * smk_366[k];

        t_457[k] = f_3 * pc_z[k] * smk_363[k];

        t_458[k] = f_18 * slk_221[k]
                   + f_3 * pc_y[k] * smk_365[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_z, slk_369, slk_370, smi0_289, \
                         smi0_290, smi1_289, smi1_290, smk_366, smk_369, \
                         smk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_19 * slk_369[k]
                   + f_6 * smi0_289[k]
                   - f_7 * smi1_289[k]
                   + f_3 * pc_x[k] * smk_369[k];

        t_460[k] = f_19 * slk_370[k]
                   + f_8 * smi0_290[k]
                   - f_9 * smi1_290[k]
                   + f_3 * pc_x[k] * smk_370[k];

        t_461[k] = f_3 * pc_z[k] * smk_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_x, pc_y, slk_225, slk_372, slk_374, smi0_292, \
                         smi0_294, smi1_292, smi1_294, smk_369, smk_372, \
                         smk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_19 * slk_372[k]
                   + f_8 * smi0_292[k]
                   - f_9 * smi1_292[k]
                   + f_3 * pc_x[k] * smk_372[k];

        t_463[k] = f_18 * slk_225[k]
                   + f_3 * pc_y[k] * smk_369[k];

        t_464[k] = f_19 * slk_374[k]
                   + f_8 * smi0_294[k]
                   - f_9 * smi1_294[k]
                   + f_3 * pc_x[k] * smk_374[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_z, slk_375, slk_377, smi0_295, \
                         smi0_297, smi1_295, smi1_297, smk_370, smk_375, \
                         smk_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_19 * slk_375[k]
                   + f_10 * smi0_295[k]
                   - f_11 * smi1_295[k]
                   + f_3 * pc_x[k] * smk_375[k];

        t_466[k] = f_3 * pc_z[k] * smk_370[k];

        t_467[k] = f_19 * slk_377[k]
                   + f_10 * smi0_297[k]
                   - f_11 * smi1_297[k]
                   + f_3 * pc_x[k] * smk_377[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, slk_230, slk_378, slk_380, smi0_298, \
                         smi0_300, smi1_298, smi1_300, smk_374, smk_378, \
                         smk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_19 * slk_378[k]
                   + f_10 * smi0_298[k]
                   - f_11 * smi1_298[k]
                   + f_3 * pc_x[k] * smk_378[k];

        t_469[k] = f_18 * slk_230[k]
                   + f_3 * pc_y[k] * smk_374[k];

        t_470[k] = f_19 * slk_380[k]
                   + f_10 * smi0_300[k]
                   - f_11 * smi1_300[k]
                   + f_3 * pc_x[k] * smk_380[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_x, pc_z, slk_381, slk_383, smi0_301, \
                         smi0_303, smi1_301, smi1_303, smk_375, smk_381, \
                         smk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_19 * slk_381[k]
                   + f_12 * smi0_301[k]
                   - f_13 * smi1_301[k]
                   + f_3 * pc_x[k] * smk_381[k];

        t_472[k] = f_3 * pc_z[k] * smk_375[k];

        t_473[k] = f_19 * slk_383[k]
                   + f_12 * smi0_303[k]
                   - f_13 * smi1_303[k]
                   + f_3 * pc_x[k] * smk_383[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, slk_236, slk_384, slk_385, smi0_304, \
                         smi0_305, smi1_304, smi1_305, smk_380, smk_384, \
                         smk_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_19 * slk_384[k]
                   + f_12 * smi0_304[k]
                   - f_13 * smi1_304[k]
                   + f_3 * pc_x[k] * smk_384[k];

        t_475[k] = f_19 * slk_385[k]
                   + f_12 * smi0_305[k]
                   - f_13 * smi1_305[k]
                   + f_3 * pc_x[k] * smk_385[k];

        t_476[k] = f_18 * slk_236[k]
                   + f_3 * pc_y[k] * smk_380[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pc_x, slk_387, slk_388, slk_389, slk_390, \
                         smi0_307, smi1_307, smk_387, smk_388, smk_389, \
                         smk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_19 * slk_387[k]
                   + f_12 * smi0_307[k]
                   - f_13 * smi1_307[k]
                   + f_3 * pc_x[k] * smk_387[k];

        t_478[k] = f_19 * slk_388[k]
                   + f_3 * pc_x[k] * smk_388[k];

        t_479[k] = f_19 * slk_389[k]
                   + f_3 * pc_x[k] * smk_389[k];

        t_480[k] = f_19 * slk_390[k]
                   + f_3 * pc_x[k] * smk_390[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, pc_x, slk_391, slk_392, slk_393, \
                         slk_394, slk_395, smk_391, smk_392, smk_393, smk_394, \
                         smk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_19 * slk_391[k]
                   + f_3 * pc_x[k] * smk_391[k];

        t_482[k] = f_19 * slk_392[k]
                   + f_3 * pc_x[k] * smk_392[k];

        t_483[k] = f_19 * slk_393[k]
                   + f_3 * pc_x[k] * smk_393[k];

        t_484[k] = f_19 * slk_394[k]
                   + f_3 * pc_x[k] * smk_394[k];

        t_485[k] = f_19 * slk_395[k]
                   + f_3 * pc_x[k] * smk_395[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pc_y, pc_z, slk_244, slk_246, smi0_301, \
                         smi0_303, smi1_301, smi1_303, smk_388, \
                         smk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_18 * slk_244[k]
                   + f_1 * smi0_301[k]
                   - f_2 * smi1_301[k]
                   + f_3 * pc_y[k] * smk_388[k];

        t_487[k] = f_3 * pc_z[k] * smk_388[k];

        t_488[k] = f_18 * slk_246[k]
                   + f_4 * smi0_303[k]
                   - f_5 * smi1_303[k]
                   + f_3 * pc_y[k] * smk_390[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_y, slk_247, slk_248, slk_249, smi0_304, \
                         smi0_305, smi0_306, smi1_304, smi1_305, smi1_306, smk_391, smk_392, \
                         smk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_18 * slk_247[k]
                   + f_6 * smi0_304[k]
                   - f_7 * smi1_304[k]
                   + f_3 * pc_y[k] * smk_391[k];

        t_490[k] = f_18 * slk_248[k]
                   + f_8 * smi0_305[k]
                   - f_9 * smi1_305[k]
                   + f_3 * pc_y[k] * smk_392[k];

        t_491[k] = f_18 * slk_249[k]
                   + f_10 * smi0_306[k]
                   - f_11 * smi1_306[k]
                   + f_3 * pc_y[k] * smk_393[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pb_z, pc_y, pc_z, sll0_270, slk_250, \
                         slk_251, sll1_270, smi0_307, smi1_307, smk_394, \
                         smk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_18 * slk_250[k]
                   + f_12 * smi0_307[k]
                   - f_13 * smi1_307[k]
                   + f_3 * pc_y[k] * smk_394[k];

        t_493[k] = f_18 * slk_251[k]
                   + f_3 * pc_y[k] * smk_395[k];

        t_494[k] = f_1 * smi0_307[k]
                   - f_2 * smi1_307[k]
                   + f_3 * pc_z[k] * smk_395[k];

        t_495[k] = pb_z[k] * sll0_270[k]
                   - f_14 * pc_z[k] * sll1_270[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_z, pc_y, pc_z, sll0_273, slk_216, \
                         slk_252, slk_254, sll1_273, smk_396, smk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_17 * slk_252[k]
                   + f_3 * pc_y[k] * smk_396[k];

        t_497[k] = f_15 * slk_216[k]
                   + f_3 * pc_z[k] * smk_396[k];

        t_498[k] = pb_z[k] * sll0_273[k]
                   - f_14 * pc_z[k] * sll1_273[k];

        t_499[k] = f_17 * slk_254[k]
                   + f_3 * pc_y[k] * smk_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pb_z, pc_x, pc_z, sll0_276, slk_219, slk_401, \
                         sll1_276, smi0_313, smi1_313, smk_399, \
                         smk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_19 * slk_401[k]
                   + f_4 * smi0_313[k]
                   - f_5 * smi1_313[k]
                   + f_3 * pc_x[k] * smk_401[k];

        t_501[k] = pb_z[k] * sll0_276[k]
                   - f_14 * pc_z[k] * sll1_276[k];

        t_502[k] = f_15 * slk_219[k]
                   + f_3 * pc_z[k] * smk_399[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pb_z, pc_x, pc_y, pc_z, sll0_280, slk_257, \
                         slk_405, sll1_280, smi0_317, smi1_317, smk_401, \
                         smk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_17 * slk_257[k]
                   + f_3 * pc_y[k] * smk_401[k];

        t_504[k] = f_19 * slk_405[k]
                   + f_6 * smi0_317[k]
                   - f_7 * smi1_317[k]
                   + f_3 * pc_x[k] * smk_405[k];

        t_505[k] = pb_z[k] * sll0_280[k]
                   - f_14 * pc_z[k] * sll1_280[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pb_z, pc_y, pc_z, sll0_282, slk_222, slk_223, \
                         slk_261, sll1_282, smk_402, smk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_15 * slk_222[k]
                   + f_3 * pc_z[k] * smk_402[k];

        t_507[k] = pb_z[k] * sll0_282[k]
                   + f_16 * slk_223[k]
                   - f_14 * pc_z[k] * sll1_282[k];

        t_508[k] = f_17 * slk_261[k]
                   + f_3 * pc_y[k] * smk_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pb_z, pc_x, pc_z, sll0_285, slk_226, slk_410, \
                         sll1_285, smi0_322, smi1_322, smk_406, \
                         smk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_19 * slk_410[k]
                   + f_8 * smi0_322[k]
                   - f_9 * smi1_322[k]
                   + f_3 * pc_x[k] * smk_410[k];

        t_510[k] = pb_z[k] * sll0_285[k]
                   - f_14 * pc_z[k] * sll1_285[k];

        t_511[k] = f_15 * slk_226[k]
                   + f_3 * pc_z[k] * smk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_z, pc_y, pc_z, sll0_287, sll0_288, slk_227, \
                         slk_228, slk_266, sll1_287, sll1_288, \
                         smk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_z[k] * sll0_287[k]
                   + f_16 * slk_227[k]
                   - f_14 * pc_z[k] * sll1_287[k];

        t_513[k] = pb_z[k] * sll0_288[k]
                   + f_17 * slk_228[k]
                   - f_14 * pc_z[k] * sll1_288[k];

        t_514[k] = f_17 * slk_266[k]
                   + f_3 * pc_y[k] * smk_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pb_z, pc_x, pc_z, sll0_291, slk_231, slk_416, \
                         sll1_291, smi0_328, smi1_328, smk_411, \
                         smk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_19 * slk_416[k]
                   + f_10 * smi0_328[k]
                   - f_11 * smi1_328[k]
                   + f_3 * pc_x[k] * smk_416[k];

        t_516[k] = pb_z[k] * sll0_291[k]
                   - f_14 * pc_z[k] * sll1_291[k];

        t_517[k] = f_15 * slk_231[k]
                   + f_3 * pc_z[k] * smk_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pb_z, pc_z, sll0_293, sll0_294, sll0_295, \
                         slk_232, slk_233, slk_234, sll1_293, sll1_294, \
                         sll1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pb_z[k] * sll0_293[k]
                   + f_16 * slk_232[k]
                   - f_14 * pc_z[k] * sll1_293[k];

        t_519[k] = pb_z[k] * sll0_294[k]
                   + f_17 * slk_233[k]
                   - f_14 * pc_z[k] * sll1_294[k];

        t_520[k] = pb_z[k] * sll0_295[k]
                   + f_18 * slk_234[k]
                   - f_14 * pc_z[k] * sll1_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, slk_272, slk_423, slk_424, \
                         slk_425, smi0_335, smi1_335, smk_416, smk_423, smk_424, \
                         smk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * slk_272[k]
                   + f_3 * pc_y[k] * smk_416[k];

        t_522[k] = f_19 * slk_423[k]
                   + f_12 * smi0_335[k]
                   - f_13 * smi1_335[k]
                   + f_3 * pc_x[k] * smk_423[k];

        t_523[k] = f_19 * slk_424[k]
                   + f_3 * pc_x[k] * smk_424[k];

        t_524[k] = f_19 * slk_425[k]
                   + f_3 * pc_x[k] * smk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, slk_426, slk_427, slk_428, \
                         slk_429, slk_430, smk_426, smk_427, smk_428, smk_429, \
                         smk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_19 * slk_426[k]
                   + f_3 * pc_x[k] * smk_426[k];

        t_526[k] = f_19 * slk_427[k]
                   + f_3 * pc_x[k] * smk_427[k];

        t_527[k] = f_19 * slk_428[k]
                   + f_3 * pc_x[k] * smk_428[k];

        t_528[k] = f_19 * slk_429[k]
                   + f_3 * pc_x[k] * smk_429[k];

        t_529[k] = f_19 * slk_430[k]
                   + f_3 * pc_x[k] * smk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pb_z, pc_x, pc_z, sll0_306, slk_244, slk_431, \
                         sll1_306, smk_424, smk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_19 * slk_431[k]
                   + f_3 * pc_x[k] * smk_431[k];

        t_531[k] = pb_z[k] * sll0_306[k]
                   - f_14 * pc_z[k] * sll1_306[k];

        t_532[k] = f_15 * slk_244[k]
                   + f_3 * pc_z[k] * smk_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pc_y, slk_282, slk_283, slk_284, smi0_331, \
                         smi0_332, smi0_333, smi1_331, smi1_332, smi1_333, smk_426, smk_427, \
                         smk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * slk_282[k]
                   + f_4 * smi0_331[k]
                   - f_5 * smi1_331[k]
                   + f_3 * pc_y[k] * smk_426[k];

        t_534[k] = f_17 * slk_283[k]
                   + f_6 * smi0_332[k]
                   - f_7 * smi1_332[k]
                   + f_3 * pc_y[k] * smk_427[k];

        t_535[k] = f_17 * slk_284[k]
                   + f_8 * smi0_333[k]
                   - f_9 * smi1_333[k]
                   + f_3 * pc_y[k] * smk_428[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pc_y, slk_285, slk_286, slk_287, smi0_334, \
                         smi0_335, smi1_334, smi1_335, smk_429, smk_430, \
                         smk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_17 * slk_285[k]
                   + f_10 * smi0_334[k]
                   - f_11 * smi1_334[k]
                   + f_3 * pc_y[k] * smk_429[k];

        t_537[k] = f_17 * slk_286[k]
                   + f_12 * smi0_335[k]
                   - f_13 * smi1_335[k]
                   + f_3 * pc_y[k] * smk_430[k];

        t_538[k] = f_17 * slk_287[k]
                   + f_3 * pc_y[k] * smk_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pc_x, pc_y, pc_z, slk_251, slk_288, slk_432, \
                         smi0_335, smi0_336, smi1_335, smi1_336, smk_431, \
                         smk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * slk_251[k]
                   + f_1 * smi0_335[k]
                   - f_2 * smi1_335[k]
                   + f_3 * pc_z[k] * smk_431[k];

        t_540[k] = f_19 * slk_432[k]
                   + f_1 * smi0_336[k]
                   - f_2 * smi1_336[k]
                   + f_3 * pc_x[k] * smk_432[k];

        t_541[k] = f_16 * slk_288[k]
                   + f_3 * pc_y[k] * smk_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, pc_y, pc_z, slk_252, slk_290, slk_435, \
                         smi0_339, smi1_339, smk_432, smk_434, \
                         smk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * slk_252[k]
                   + f_3 * pc_z[k] * smk_432[k];

        t_543[k] = f_19 * slk_435[k]
                   + f_4 * smi0_339[k]
                   - f_5 * smi1_339[k]
                   + f_3 * pc_x[k] * smk_435[k];

        t_544[k] = f_16 * slk_290[k]
                   + f_3 * pc_y[k] * smk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, slk_255, slk_437, slk_438, smi0_341, \
                         smi0_342, smi1_341, smi1_342, smk_435, smk_437, \
                         smk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_19 * slk_437[k]
                   + f_4 * smi0_341[k]
                   - f_5 * smi1_341[k]
                   + f_3 * pc_x[k] * smk_437[k];

        t_546[k] = f_19 * slk_438[k]
                   + f_6 * smi0_342[k]
                   - f_7 * smi1_342[k]
                   + f_3 * pc_x[k] * smk_438[k];

        t_547[k] = f_16 * slk_255[k]
                   + f_3 * pc_z[k] * smk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pc_x, pc_y, slk_293, slk_441, slk_442, smi0_345, \
                         smi0_346, smi1_345, smi1_346, smk_437, smk_441, \
                         smk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * slk_293[k]
                   + f_3 * pc_y[k] * smk_437[k];

        t_549[k] = f_19 * slk_441[k]
                   + f_6 * smi0_345[k]
                   - f_7 * smi1_345[k]
                   + f_3 * pc_x[k] * smk_441[k];

        t_550[k] = f_19 * slk_442[k]
                   + f_8 * smi0_346[k]
                   - f_9 * smi1_346[k]
                   + f_3 * pc_x[k] * smk_442[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pc_x, pc_y, pc_z, slk_258, slk_297, slk_444, \
                         smi0_348, smi1_348, smk_438, smk_441, \
                         smk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * slk_258[k]
                   + f_3 * pc_z[k] * smk_438[k];

        t_552[k] = f_19 * slk_444[k]
                   + f_8 * smi0_348[k]
                   - f_9 * smi1_348[k]
                   + f_3 * pc_x[k] * smk_444[k];

        t_553[k] = f_16 * slk_297[k]
                   + f_3 * pc_y[k] * smk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_x, pc_z, slk_262, slk_446, slk_447, smi0_350, \
                         smi0_351, smi1_350, smi1_351, smk_442, smk_446, \
                         smk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_19 * slk_446[k]
                   + f_8 * smi0_350[k]
                   - f_9 * smi1_350[k]
                   + f_3 * pc_x[k] * smk_446[k];

        t_555[k] = f_19 * slk_447[k]
                   + f_10 * smi0_351[k]
                   - f_11 * smi1_351[k]
                   + f_3 * pc_x[k] * smk_447[k];

        t_556[k] = f_16 * slk_262[k]
                   + f_3 * pc_z[k] * smk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, pc_y, slk_302, slk_449, slk_450, smi0_353, \
                         smi0_354, smi1_353, smi1_354, smk_446, smk_449, \
                         smk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_19 * slk_449[k]
                   + f_10 * smi0_353[k]
                   - f_11 * smi1_353[k]
                   + f_3 * pc_x[k] * smk_449[k];

        t_558[k] = f_19 * slk_450[k]
                   + f_10 * smi0_354[k]
                   - f_11 * smi1_354[k]
                   + f_3 * pc_x[k] * smk_450[k];

        t_559[k] = f_16 * slk_302[k]
                   + f_3 * pc_y[k] * smk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pc_x, pc_z, slk_267, slk_452, slk_453, smi0_356, \
                         smi0_357, smi1_356, smi1_357, smk_447, smk_452, \
                         smk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_19 * slk_452[k]
                   + f_10 * smi0_356[k]
                   - f_11 * smi1_356[k]
                   + f_3 * pc_x[k] * smk_452[k];

        t_561[k] = f_19 * slk_453[k]
                   + f_12 * smi0_357[k]
                   - f_13 * smi1_357[k]
                   + f_3 * pc_x[k] * smk_453[k];

        t_562[k] = f_16 * slk_267[k]
                   + f_3 * pc_z[k] * smk_447[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_405 = buffer.data(sll0 + 405);
    const auto *sll0_408 = buffer.data(sll0 + 408);
    const auto *sll0_410 = buffer.data(sll0 + 410);
    const auto *sll0_411 = buffer.data(sll0 + 411);
    const auto *sll0_414 = buffer.data(sll0 + 414);
    const auto *sll0_415 = buffer.data(sll0 + 415);
    const auto *sll0_417 = buffer.data(sll0 + 417);
    const auto *sll0_419 = buffer.data(sll0 + 419);
    const auto *sll0_420 = buffer.data(sll0 + 420);
    const auto *sll0_422 = buffer.data(sll0 + 422);
    const auto *sll0_423 = buffer.data(sll0 + 423);
    const auto *sll0_425 = buffer.data(sll0 + 425);
    const auto *sll0_426 = buffer.data(sll0 + 426);
    const auto *sll0_428 = buffer.data(sll0 + 428);
    const auto *sll0_429 = buffer.data(sll0 + 429);
    const auto *sll0_430 = buffer.data(sll0 + 430);
    const auto *sll0_432 = buffer.data(sll0 + 432);
    const auto *sll0_449 = buffer.data(sll0 + 449);

    const auto *slk_280 = buffer.data(slk + 280);
    const auto *slk_287 = buffer.data(slk + 287);
    const auto *slk_288 = buffer.data(slk + 288);
    const auto *slk_291 = buffer.data(slk + 291);
    const auto *slk_294 = buffer.data(slk + 294);
    const auto *slk_298 = buffer.data(slk + 298);
    const auto *slk_303 = buffer.data(slk + 303);
    const auto *slk_308 = buffer.data(slk + 308);
    const auto *slk_316 = buffer.data(slk + 316);
    const auto *slk_318 = buffer.data(slk + 318);
    const auto *slk_319 = buffer.data(slk + 319);
    const auto *slk_320 = buffer.data(slk + 320);
    const auto *slk_321 = buffer.data(slk + 321);
    const auto *slk_322 = buffer.data(slk + 322);
    const auto *slk_323 = buffer.data(slk + 323);
    const auto *slk_324 = buffer.data(slk + 324);
    const auto *slk_325 = buffer.data(slk + 325);
    const auto *slk_326 = buffer.data(slk + 326);
    const auto *slk_327 = buffer.data(slk + 327);
    const auto *slk_329 = buffer.data(slk + 329);
    const auto *slk_330 = buffer.data(slk + 330);
    const auto *slk_332 = buffer.data(slk + 332);
    const auto *slk_333 = buffer.data(slk + 333);
    const auto *slk_334 = buffer.data(slk + 334);
    const auto *slk_336 = buffer.data(slk + 336);
    const auto *slk_337 = buffer.data(slk + 337);
    const auto *slk_338 = buffer.data(slk + 338);
    const auto *slk_339 = buffer.data(slk + 339);
    const auto *slk_341 = buffer.data(slk + 341);
    const auto *slk_342 = buffer.data(slk + 342);
    const auto *slk_343 = buffer.data(slk + 343);
    const auto *slk_344 = buffer.data(slk + 344);
    const auto *slk_352 = buffer.data(slk + 352);
    const auto *slk_354 = buffer.data(slk + 354);
    const auto *slk_355 = buffer.data(slk + 355);
    const auto *slk_356 = buffer.data(slk + 356);
    const auto *slk_357 = buffer.data(slk + 357);
    const auto *slk_358 = buffer.data(slk + 358);
    const auto *slk_359 = buffer.data(slk + 359);
    const auto *slk_455 = buffer.data(slk + 455);
    const auto *slk_456 = buffer.data(slk + 456);
    const auto *slk_457 = buffer.data(slk + 457);
    const auto *slk_459 = buffer.data(slk + 459);
    const auto *slk_460 = buffer.data(slk + 460);
    const auto *slk_461 = buffer.data(slk + 461);
    const auto *slk_462 = buffer.data(slk + 462);
    const auto *slk_463 = buffer.data(slk + 463);
    const auto *slk_464 = buffer.data(slk + 464);
    const auto *slk_465 = buffer.data(slk + 465);
    const auto *slk_466 = buffer.data(slk + 466);
    const auto *slk_467 = buffer.data(slk + 467);
    const auto *slk_496 = buffer.data(slk + 496);
    const auto *slk_497 = buffer.data(slk + 497);
    const auto *slk_498 = buffer.data(slk + 498);
    const auto *slk_499 = buffer.data(slk + 499);
    const auto *slk_500 = buffer.data(slk + 500);
    const auto *slk_501 = buffer.data(slk + 501);
    const auto *slk_502 = buffer.data(slk + 502);
    const auto *slk_503 = buffer.data(slk + 503);
    const auto *slk_504 = buffer.data(slk + 504);
    const auto *slk_507 = buffer.data(slk + 507);
    const auto *slk_509 = buffer.data(slk + 509);
    const auto *slk_510 = buffer.data(slk + 510);
    const auto *slk_513 = buffer.data(slk + 513);
    const auto *slk_514 = buffer.data(slk + 514);
    const auto *slk_516 = buffer.data(slk + 516);
    const auto *slk_518 = buffer.data(slk + 518);
    const auto *slk_519 = buffer.data(slk + 519);
    const auto *slk_521 = buffer.data(slk + 521);
    const auto *slk_522 = buffer.data(slk + 522);
    const auto *slk_524 = buffer.data(slk + 524);
    const auto *slk_525 = buffer.data(slk + 525);
    const auto *slk_527 = buffer.data(slk + 527);
    const auto *slk_528 = buffer.data(slk + 528);
    const auto *slk_529 = buffer.data(slk + 529);
    const auto *slk_531 = buffer.data(slk + 531);
    const auto *slk_532 = buffer.data(slk + 532);
    const auto *slk_533 = buffer.data(slk + 533);
    const auto *slk_534 = buffer.data(slk + 534);
    const auto *slk_535 = buffer.data(slk + 535);
    const auto *slk_536 = buffer.data(slk + 536);
    const auto *slk_537 = buffer.data(slk + 537);
    const auto *slk_538 = buffer.data(slk + 538);
    const auto *slk_539 = buffer.data(slk + 539);

    const auto *sll1_405 = buffer.data(sll1 + 405);
    const auto *sll1_408 = buffer.data(sll1 + 408);
    const auto *sll1_410 = buffer.data(sll1 + 410);
    const auto *sll1_411 = buffer.data(sll1 + 411);
    const auto *sll1_414 = buffer.data(sll1 + 414);
    const auto *sll1_415 = buffer.data(sll1 + 415);
    const auto *sll1_417 = buffer.data(sll1 + 417);
    const auto *sll1_419 = buffer.data(sll1 + 419);
    const auto *sll1_420 = buffer.data(sll1 + 420);
    const auto *sll1_422 = buffer.data(sll1 + 422);
    const auto *sll1_423 = buffer.data(sll1 + 423);
    const auto *sll1_425 = buffer.data(sll1 + 425);
    const auto *sll1_426 = buffer.data(sll1 + 426);
    const auto *sll1_428 = buffer.data(sll1 + 428);
    const auto *sll1_429 = buffer.data(sll1 + 429);
    const auto *sll1_430 = buffer.data(sll1 + 430);
    const auto *sll1_432 = buffer.data(sll1 + 432);
    const auto *sll1_449 = buffer.data(sll1 + 449);

    const auto *smi0_357 = buffer.data(smi0 + 357);
    const auto *smi0_359 = buffer.data(smi0 + 359);
    const auto *smi0_360 = buffer.data(smi0 + 360);
    const auto *smi0_361 = buffer.data(smi0 + 361);
    const auto *smi0_362 = buffer.data(smi0 + 362);
    const auto *smi0_363 = buffer.data(smi0 + 363);
    const auto *smi0_385 = buffer.data(smi0 + 385);
    const auto *smi0_387 = buffer.data(smi0 + 387);
    const auto *smi0_388 = buffer.data(smi0 + 388);
    const auto *smi0_389 = buffer.data(smi0 + 389);
    const auto *smi0_390 = buffer.data(smi0 + 390);
    const auto *smi0_391 = buffer.data(smi0 + 391);
    const auto *smi0_392 = buffer.data(smi0 + 392);
    const auto *smi0_395 = buffer.data(smi0 + 395);
    const auto *smi0_397 = buffer.data(smi0 + 397);
    const auto *smi0_398 = buffer.data(smi0 + 398);
    const auto *smi0_401 = buffer.data(smi0 + 401);
    const auto *smi0_402 = buffer.data(smi0 + 402);
    const auto *smi0_404 = buffer.data(smi0 + 404);
    const auto *smi0_406 = buffer.data(smi0 + 406);
    const auto *smi0_407 = buffer.data(smi0 + 407);
    const auto *smi0_409 = buffer.data(smi0 + 409);
    const auto *smi0_410 = buffer.data(smi0 + 410);
    const auto *smi0_412 = buffer.data(smi0 + 412);
    const auto *smi0_413 = buffer.data(smi0 + 413);
    const auto *smi0_415 = buffer.data(smi0 + 415);
    const auto *smi0_416 = buffer.data(smi0 + 416);
    const auto *smi0_417 = buffer.data(smi0 + 417);
    const auto *smi0_418 = buffer.data(smi0 + 418);
    const auto *smi0_419 = buffer.data(smi0 + 419);

    const auto *smi1_357 = buffer.data(smi1 + 357);
    const auto *smi1_359 = buffer.data(smi1 + 359);
    const auto *smi1_360 = buffer.data(smi1 + 360);
    const auto *smi1_361 = buffer.data(smi1 + 361);
    const auto *smi1_362 = buffer.data(smi1 + 362);
    const auto *smi1_363 = buffer.data(smi1 + 363);
    const auto *smi1_385 = buffer.data(smi1 + 385);
    const auto *smi1_387 = buffer.data(smi1 + 387);
    const auto *smi1_388 = buffer.data(smi1 + 388);
    const auto *smi1_389 = buffer.data(smi1 + 389);
    const auto *smi1_390 = buffer.data(smi1 + 390);
    const auto *smi1_391 = buffer.data(smi1 + 391);
    const auto *smi1_392 = buffer.data(smi1 + 392);
    const auto *smi1_395 = buffer.data(smi1 + 395);
    const auto *smi1_397 = buffer.data(smi1 + 397);
    const auto *smi1_398 = buffer.data(smi1 + 398);
    const auto *smi1_401 = buffer.data(smi1 + 401);
    const auto *smi1_402 = buffer.data(smi1 + 402);
    const auto *smi1_404 = buffer.data(smi1 + 404);
    const auto *smi1_406 = buffer.data(smi1 + 406);
    const auto *smi1_407 = buffer.data(smi1 + 407);
    const auto *smi1_409 = buffer.data(smi1 + 409);
    const auto *smi1_410 = buffer.data(smi1 + 410);
    const auto *smi1_412 = buffer.data(smi1 + 412);
    const auto *smi1_413 = buffer.data(smi1 + 413);
    const auto *smi1_415 = buffer.data(smi1 + 415);
    const auto *smi1_416 = buffer.data(smi1 + 416);
    const auto *smi1_417 = buffer.data(smi1 + 417);
    const auto *smi1_418 = buffer.data(smi1 + 418);
    const auto *smi1_419 = buffer.data(smi1 + 419);

    const auto *smk_452 = buffer.data(smk + 452);
    const auto *smk_455 = buffer.data(smk + 455);
    const auto *smk_456 = buffer.data(smk + 456);
    const auto *smk_457 = buffer.data(smk + 457);
    const auto *smk_459 = buffer.data(smk + 459);
    const auto *smk_460 = buffer.data(smk + 460);
    const auto *smk_461 = buffer.data(smk + 461);
    const auto *smk_462 = buffer.data(smk + 462);
    const auto *smk_463 = buffer.data(smk + 463);
    const auto *smk_464 = buffer.data(smk + 464);
    const auto *smk_465 = buffer.data(smk + 465);
    const auto *smk_466 = buffer.data(smk + 466);
    const auto *smk_467 = buffer.data(smk + 467);
    const auto *smk_468 = buffer.data(smk + 468);
    const auto *smk_470 = buffer.data(smk + 470);
    const auto *smk_471 = buffer.data(smk + 471);
    const auto *smk_473 = buffer.data(smk + 473);
    const auto *smk_474 = buffer.data(smk + 474);
    const auto *smk_477 = buffer.data(smk + 477);
    const auto *smk_478 = buffer.data(smk + 478);
    const auto *smk_482 = buffer.data(smk + 482);
    const auto *smk_483 = buffer.data(smk + 483);
    const auto *smk_488 = buffer.data(smk + 488);
    const auto *smk_496 = buffer.data(smk + 496);
    const auto *smk_497 = buffer.data(smk + 497);
    const auto *smk_498 = buffer.data(smk + 498);
    const auto *smk_499 = buffer.data(smk + 499);
    const auto *smk_500 = buffer.data(smk + 500);
    const auto *smk_501 = buffer.data(smk + 501);
    const auto *smk_502 = buffer.data(smk + 502);
    const auto *smk_503 = buffer.data(smk + 503);
    const auto *smk_504 = buffer.data(smk + 504);
    const auto *smk_506 = buffer.data(smk + 506);
    const auto *smk_507 = buffer.data(smk + 507);
    const auto *smk_509 = buffer.data(smk + 509);
    const auto *smk_510 = buffer.data(smk + 510);
    const auto *smk_513 = buffer.data(smk + 513);
    const auto *smk_514 = buffer.data(smk + 514);
    const auto *smk_516 = buffer.data(smk + 516);
    const auto *smk_518 = buffer.data(smk + 518);
    const auto *smk_519 = buffer.data(smk + 519);
    const auto *smk_521 = buffer.data(smk + 521);
    const auto *smk_522 = buffer.data(smk + 522);
    const auto *smk_524 = buffer.data(smk + 524);
    const auto *smk_525 = buffer.data(smk + 525);
    const auto *smk_527 = buffer.data(smk + 527);
    const auto *smk_528 = buffer.data(smk + 528);
    const auto *smk_529 = buffer.data(smk + 529);
    const auto *smk_531 = buffer.data(smk + 531);
    const auto *smk_532 = buffer.data(smk + 532);
    const auto *smk_533 = buffer.data(smk + 533);
    const auto *smk_534 = buffer.data(smk + 534);
    const auto *smk_535 = buffer.data(smk + 535);
    const auto *smk_536 = buffer.data(smk + 536);
    const auto *smk_537 = buffer.data(smk + 537);
    const auto *smk_538 = buffer.data(smk + 538);
    const auto *smk_539 = buffer.data(smk + 539);

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, slk_455, slk_456, slk_457, smi0_359, \
                         smi0_360, smi0_361, smi1_359, smi1_360, smi1_361, smk_455, smk_456, \
                         smk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_19 * slk_455[k]
                   + f_12 * smi0_359[k]
                   - f_13 * smi1_359[k]
                   + f_3 * pc_x[k] * smk_455[k];

        t_564[k] = f_19 * slk_456[k]
                   + f_12 * smi0_360[k]
                   - f_13 * smi1_360[k]
                   + f_3 * pc_x[k] * smk_456[k];

        t_565[k] = f_19 * slk_457[k]
                   + f_12 * smi0_361[k]
                   - f_13 * smi1_361[k]
                   + f_3 * pc_x[k] * smk_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pc_x, pc_y, slk_308, slk_459, slk_460, \
                         slk_461, smi0_363, smi1_363, smk_452, smk_459, smk_460, \
                         smk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * slk_308[k]
                   + f_3 * pc_y[k] * smk_452[k];

        t_567[k] = f_19 * slk_459[k]
                   + f_12 * smi0_363[k]
                   - f_13 * smi1_363[k]
                   + f_3 * pc_x[k] * smk_459[k];

        t_568[k] = f_19 * slk_460[k]
                   + f_3 * pc_x[k] * smk_460[k];

        t_569[k] = f_19 * slk_461[k]
                   + f_3 * pc_x[k] * smk_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, slk_462, slk_463, slk_464, \
                         slk_465, slk_466, smk_462, smk_463, smk_464, smk_465, \
                         smk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_19 * slk_462[k]
                   + f_3 * pc_x[k] * smk_462[k];

        t_571[k] = f_19 * slk_463[k]
                   + f_3 * pc_x[k] * smk_463[k];

        t_572[k] = f_19 * slk_464[k]
                   + f_3 * pc_x[k] * smk_464[k];

        t_573[k] = f_19 * slk_465[k]
                   + f_3 * pc_x[k] * smk_465[k];

        t_574[k] = f_19 * slk_466[k]
                   + f_3 * pc_x[k] * smk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, pc_z, slk_280, slk_316, slk_467, \
                         smi0_357, smi1_357, smk_460, smk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_19 * slk_467[k]
                   + f_3 * pc_x[k] * smk_467[k];

        t_576[k] = f_16 * slk_316[k]
                   + f_1 * smi0_357[k]
                   - f_2 * smi1_357[k]
                   + f_3 * pc_y[k] * smk_460[k];

        t_577[k] = f_16 * slk_280[k]
                   + f_3 * pc_z[k] * smk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, slk_318, slk_319, slk_320, smi0_359, \
                         smi0_360, smi0_361, smi1_359, smi1_360, smi1_361, smk_462, smk_463, \
                         smk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * slk_318[k]
                   + f_4 * smi0_359[k]
                   - f_5 * smi1_359[k]
                   + f_3 * pc_y[k] * smk_462[k];

        t_579[k] = f_16 * slk_319[k]
                   + f_6 * smi0_360[k]
                   - f_7 * smi1_360[k]
                   + f_3 * pc_y[k] * smk_463[k];

        t_580[k] = f_16 * slk_320[k]
                   + f_8 * smi0_361[k]
                   - f_9 * smi1_361[k]
                   + f_3 * pc_y[k] * smk_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, slk_321, slk_322, slk_323, smi0_362, \
                         smi0_363, smi1_362, smi1_363, smk_465, smk_466, \
                         smk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * slk_321[k]
                   + f_10 * smi0_362[k]
                   - f_11 * smi1_362[k]
                   + f_3 * pc_y[k] * smk_465[k];

        t_582[k] = f_16 * slk_322[k]
                   + f_12 * smi0_363[k]
                   - f_13 * smi1_363[k]
                   + f_3 * pc_y[k] * smk_466[k];

        t_583[k] = f_16 * slk_323[k]
                   + f_3 * pc_y[k] * smk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_y, pc_y, pc_z, sll0_405, slk_287, \
                         slk_288, slk_324, sll1_405, smi0_363, smi1_363, smk_467, \
                         smk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * slk_287[k]
                   + f_1 * smi0_363[k]
                   - f_2 * smi1_363[k]
                   + f_3 * pc_z[k] * smk_467[k];

        t_585[k] = pb_y[k] * sll0_405[k]
                   - f_14 * pc_y[k] * sll1_405[k];

        t_586[k] = f_15 * slk_324[k]
                   + f_3 * pc_y[k] * smk_468[k];

        t_587[k] = f_17 * slk_288[k]
                   + f_3 * pc_z[k] * smk_468[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pb_y, pc_y, sll0_408, sll0_410, sll0_411, \
                         slk_325, slk_326, slk_327, sll1_408, sll1_410, sll1_411, \
                         smk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pb_y[k] * sll0_408[k]
                   + f_16 * slk_325[k]
                   - f_14 * pc_y[k] * sll1_408[k];

        t_589[k] = f_15 * slk_326[k]
                   + f_3 * pc_y[k] * smk_470[k];

        t_590[k] = pb_y[k] * sll0_410[k]
                   - f_14 * pc_y[k] * sll1_410[k];

        t_591[k] = pb_y[k] * sll0_411[k]
                   + f_17 * slk_327[k]
                   - f_14 * pc_y[k] * sll1_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pb_y, pc_y, pc_z, sll0_414, sll0_415, \
                         slk_291, slk_329, slk_330, sll1_414, sll1_415, smk_471, \
                         smk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * slk_291[k]
                   + f_3 * pc_z[k] * smk_471[k];

        t_593[k] = f_15 * slk_329[k]
                   + f_3 * pc_y[k] * smk_473[k];

        t_594[k] = pb_y[k] * sll0_414[k]
                   - f_14 * pc_y[k] * sll1_414[k];

        t_595[k] = pb_y[k] * sll0_415[k]
                   + f_18 * slk_330[k]
                   - f_14 * pc_y[k] * sll1_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pb_y, pc_y, pc_z, sll0_417, sll0_419, \
                         slk_294, slk_332, slk_333, sll1_417, sll1_419, smk_474, \
                         smk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * slk_294[k]
                   + f_3 * pc_z[k] * smk_474[k];

        t_597[k] = pb_y[k] * sll0_417[k]
                   + f_16 * slk_332[k]
                   - f_14 * pc_y[k] * sll1_417[k];

        t_598[k] = f_15 * slk_333[k]
                   + f_3 * pc_y[k] * smk_477[k];

        t_599[k] = pb_y[k] * sll0_419[k]
                   - f_14 * pc_y[k] * sll1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pb_y, pc_y, pc_z, sll0_420, sll0_422, slk_298, \
                         slk_334, slk_336, sll1_420, sll1_422, \
                         smk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pb_y[k] * sll0_420[k]
                   + f_19 * slk_334[k]
                   - f_14 * pc_y[k] * sll1_420[k];

        t_601[k] = f_17 * slk_298[k]
                   + f_3 * pc_z[k] * smk_478[k];

        t_602[k] = pb_y[k] * sll0_422[k]
                   + f_17 * slk_336[k]
                   - f_14 * pc_y[k] * sll1_422[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pb_y, pc_y, sll0_423, sll0_425, sll0_426, \
                         slk_337, slk_338, slk_339, sll1_423, sll1_425, sll1_426, \
                         smk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pb_y[k] * sll0_423[k]
                   + f_16 * slk_337[k]
                   - f_14 * pc_y[k] * sll1_423[k];

        t_604[k] = f_15 * slk_338[k]
                   + f_3 * pc_y[k] * smk_482[k];

        t_605[k] = pb_y[k] * sll0_425[k]
                   - f_14 * pc_y[k] * sll1_425[k];

        t_606[k] = pb_y[k] * sll0_426[k]
                   + f_20 * slk_339[k]
                   - f_14 * pc_y[k] * sll1_426[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pb_y, pc_y, pc_z, sll0_428, sll0_429, slk_303, \
                         slk_341, slk_342, sll1_428, sll1_429, \
                         smk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_17 * slk_303[k]
                   + f_3 * pc_z[k] * smk_483[k];

        t_608[k] = pb_y[k] * sll0_428[k]
                   + f_18 * slk_341[k]
                   - f_14 * pc_y[k] * sll1_428[k];

        t_609[k] = pb_y[k] * sll0_429[k]
                   + f_17 * slk_342[k]
                   - f_14 * pc_y[k] * sll1_429[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pb_y, pc_x, pc_y, sll0_430, sll0_432, \
                         slk_343, slk_344, slk_496, sll1_430, sll1_432, smk_488, \
                         smk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pb_y[k] * sll0_430[k]
                   + f_16 * slk_343[k]
                   - f_14 * pc_y[k] * sll1_430[k];

        t_611[k] = f_15 * slk_344[k]
                   + f_3 * pc_y[k] * smk_488[k];

        t_612[k] = pb_y[k] * sll0_432[k]
                   - f_14 * pc_y[k] * sll1_432[k];

        t_613[k] = f_19 * slk_496[k]
                   + f_3 * pc_x[k] * smk_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, slk_497, slk_498, slk_499, \
                         slk_500, slk_501, smk_497, smk_498, smk_499, smk_500, \
                         smk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_19 * slk_497[k]
                   + f_3 * pc_x[k] * smk_497[k];

        t_615[k] = f_19 * slk_498[k]
                   + f_3 * pc_x[k] * smk_498[k];

        t_616[k] = f_19 * slk_499[k]
                   + f_3 * pc_x[k] * smk_499[k];

        t_617[k] = f_19 * slk_500[k]
                   + f_3 * pc_x[k] * smk_500[k];

        t_618[k] = f_19 * slk_501[k]
                   + f_3 * pc_x[k] * smk_501[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, slk_316, slk_352, \
                         slk_502, slk_503, smi0_385, smi1_385, smk_496, smk_502, \
                         smk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_19 * slk_502[k]
                   + f_3 * pc_x[k] * smk_502[k];

        t_620[k] = f_19 * slk_503[k]
                   + f_3 * pc_x[k] * smk_503[k];

        t_621[k] = f_15 * slk_352[k]
                   + f_1 * smi0_385[k]
                   - f_2 * smi1_385[k]
                   + f_3 * pc_y[k] * smk_496[k];

        t_622[k] = f_17 * slk_316[k]
                   + f_3 * pc_z[k] * smk_496[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_y, slk_354, slk_355, slk_356, smi0_387, \
                         smi0_388, smi0_389, smi1_387, smi1_388, smi1_389, smk_498, smk_499, \
                         smk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_15 * slk_354[k]
                   + f_4 * smi0_387[k]
                   - f_5 * smi1_387[k]
                   + f_3 * pc_y[k] * smk_498[k];

        t_624[k] = f_15 * slk_355[k]
                   + f_6 * smi0_388[k]
                   - f_7 * smi1_388[k]
                   + f_3 * pc_y[k] * smk_499[k];

        t_625[k] = f_15 * slk_356[k]
                   + f_8 * smi0_389[k]
                   - f_9 * smi1_389[k]
                   + f_3 * pc_y[k] * smk_500[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, slk_357, slk_358, slk_359, smi0_390, \
                         smi0_391, smi1_390, smi1_391, smk_501, smk_502, \
                         smk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * slk_357[k]
                   + f_10 * smi0_390[k]
                   - f_11 * smi1_390[k]
                   + f_3 * pc_y[k] * smk_501[k];

        t_627[k] = f_15 * slk_358[k]
                   + f_12 * smi0_391[k]
                   - f_13 * smi1_391[k]
                   + f_3 * pc_y[k] * smk_502[k];

        t_628[k] = f_15 * slk_359[k]
                   + f_3 * pc_y[k] * smk_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pb_y, pc_x, pc_y, pc_z, sll0_449, \
                         slk_324, slk_504, sll1_449, smi0_392, smi1_392, \
                         smk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pb_y[k] * sll0_449[k]
                   - f_14 * pc_y[k] * sll1_449[k];

        t_630[k] = f_19 * slk_504[k]
                   + f_1 * smi0_392[k]
                   - f_2 * smi1_392[k]
                   + f_3 * pc_x[k] * smk_504[k];

        t_631[k] = f_3 * pc_y[k] * smk_504[k];

        t_632[k] = f_18 * slk_324[k]
                   + f_3 * pc_z[k] * smk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, slk_507, slk_509, smi0_395, \
                         smi0_397, smi1_395, smi1_397, smk_506, smk_507, \
                         smk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_19 * slk_507[k]
                   + f_4 * smi0_395[k]
                   - f_5 * smi1_395[k]
                   + f_3 * pc_x[k] * smk_507[k];

        t_634[k] = f_3 * pc_y[k] * smk_506[k];

        t_635[k] = f_19 * slk_509[k]
                   + f_4 * smi0_397[k]
                   - f_5 * smi1_397[k]
                   + f_3 * pc_x[k] * smk_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_x, pc_y, pc_z, slk_327, slk_510, smi0_398, \
                         smi1_398, smk_507, smk_509, smk_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_19 * slk_510[k]
                   + f_6 * smi0_398[k]
                   - f_7 * smi1_398[k]
                   + f_3 * pc_x[k] * smk_510[k];

        t_637[k] = f_18 * slk_327[k]
                   + f_3 * pc_z[k] * smk_507[k];

        t_638[k] = f_3 * pc_y[k] * smk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_z, slk_330, slk_513, slk_514, smi0_401, \
                         smi0_402, smi1_401, smi1_402, smk_510, smk_513, \
                         smk_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_19 * slk_513[k]
                   + f_6 * smi0_401[k]
                   - f_7 * smi1_401[k]
                   + f_3 * pc_x[k] * smk_513[k];

        t_640[k] = f_19 * slk_514[k]
                   + f_8 * smi0_402[k]
                   - f_9 * smi1_402[k]
                   + f_3 * pc_x[k] * smk_514[k];

        t_641[k] = f_18 * slk_330[k]
                   + f_3 * pc_z[k] * smk_510[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, slk_516, slk_518, smi0_404, \
                         smi0_406, smi1_404, smi1_406, smk_513, smk_516, \
                         smk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_19 * slk_516[k]
                   + f_8 * smi0_404[k]
                   - f_9 * smi1_404[k]
                   + f_3 * pc_x[k] * smk_516[k];

        t_643[k] = f_3 * pc_y[k] * smk_513[k];

        t_644[k] = f_19 * slk_518[k]
                   + f_8 * smi0_406[k]
                   - f_9 * smi1_406[k]
                   + f_3 * pc_x[k] * smk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_x, pc_z, slk_334, slk_519, slk_521, smi0_407, \
                         smi0_409, smi1_407, smi1_409, smk_514, smk_519, \
                         smk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_19 * slk_519[k]
                   + f_10 * smi0_407[k]
                   - f_11 * smi1_407[k]
                   + f_3 * pc_x[k] * smk_519[k];

        t_646[k] = f_18 * slk_334[k]
                   + f_3 * pc_z[k] * smk_514[k];

        t_647[k] = f_19 * slk_521[k]
                   + f_10 * smi0_409[k]
                   - f_11 * smi1_409[k]
                   + f_3 * pc_x[k] * smk_521[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, pc_y, slk_522, slk_524, smi0_410, \
                         smi0_412, smi1_410, smi1_412, smk_518, smk_522, \
                         smk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_19 * slk_522[k]
                   + f_10 * smi0_410[k]
                   - f_11 * smi1_410[k]
                   + f_3 * pc_x[k] * smk_522[k];

        t_649[k] = f_3 * pc_y[k] * smk_518[k];

        t_650[k] = f_19 * slk_524[k]
                   + f_10 * smi0_412[k]
                   - f_11 * smi1_412[k]
                   + f_3 * pc_x[k] * smk_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_x, pc_z, slk_339, slk_525, slk_527, smi0_413, \
                         smi0_415, smi1_413, smi1_415, smk_519, smk_525, \
                         smk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_19 * slk_525[k]
                   + f_12 * smi0_413[k]
                   - f_13 * smi1_413[k]
                   + f_3 * pc_x[k] * smk_525[k];

        t_652[k] = f_18 * slk_339[k]
                   + f_3 * pc_z[k] * smk_519[k];

        t_653[k] = f_19 * slk_527[k]
                   + f_12 * smi0_415[k]
                   - f_13 * smi1_415[k]
                   + f_3 * pc_x[k] * smk_527[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, pc_y, slk_528, slk_529, smi0_416, \
                         smi0_417, smi1_416, smi1_417, smk_524, smk_528, \
                         smk_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_19 * slk_528[k]
                   + f_12 * smi0_416[k]
                   - f_13 * smi1_416[k]
                   + f_3 * pc_x[k] * smk_528[k];

        t_655[k] = f_19 * slk_529[k]
                   + f_12 * smi0_417[k]
                   - f_13 * smi1_417[k]
                   + f_3 * pc_x[k] * smk_529[k];

        t_656[k] = f_3 * pc_y[k] * smk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pc_x, slk_531, slk_532, slk_533, slk_534, \
                         smi0_419, smi1_419, smk_531, smk_532, smk_533, \
                         smk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_19 * slk_531[k]
                   + f_12 * smi0_419[k]
                   - f_13 * smi1_419[k]
                   + f_3 * pc_x[k] * smk_531[k];

        t_658[k] = f_19 * slk_532[k]
                   + f_3 * pc_x[k] * smk_532[k];

        t_659[k] = f_19 * slk_533[k]
                   + f_3 * pc_x[k] * smk_533[k];

        t_660[k] = f_19 * slk_534[k]
                   + f_3 * pc_x[k] * smk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, slk_535, slk_536, slk_537, \
                         slk_538, slk_539, smk_535, smk_536, smk_537, smk_538, \
                         smk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_19 * slk_535[k]
                   + f_3 * pc_x[k] * smk_535[k];

        t_662[k] = f_19 * slk_536[k]
                   + f_3 * pc_x[k] * smk_536[k];

        t_663[k] = f_19 * slk_537[k]
                   + f_3 * pc_x[k] * smk_537[k];

        t_664[k] = f_19 * slk_538[k]
                   + f_3 * pc_x[k] * smk_538[k];

        t_665[k] = f_19 * slk_539[k]
                   + f_3 * pc_x[k] * smk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pc_y, pc_z, slk_352, smi0_413, smi0_415, \
                         smi0_416, smi1_413, smi1_415, smi1_416, smk_532, smk_534, \
                         smk_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * smi0_413[k]
                   - f_2 * smi1_413[k]
                   + f_3 * pc_y[k] * smk_532[k];

        t_667[k] = f_18 * slk_352[k]
                   + f_3 * pc_z[k] * smk_532[k];

        t_668[k] = f_4 * smi0_415[k]
                   - f_5 * smi1_415[k]
                   + f_3 * pc_y[k] * smk_534[k];

        t_669[k] = f_6 * smi0_416[k]
                   - f_7 * smi1_416[k]
                   + f_3 * pc_y[k] * smk_535[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pc_y, smi0_417, smi0_418, smi0_419, \
                         smi1_417, smi1_418, smi1_419, smk_536, smk_537, smk_538, \
                         smk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_8 * smi0_417[k]
                   - f_9 * smi1_417[k]
                   + f_3 * pc_y[k] * smk_536[k];

        t_671[k] = f_10 * smi0_418[k]
                   - f_11 * smi1_418[k]
                   + f_3 * pc_y[k] * smk_537[k];

        t_672[k] = f_12 * smi0_419[k]
                   - f_13 * smi1_419[k]
                   + f_3 * pc_y[k] * smk_538[k];

        t_673[k] = f_3 * pc_y[k] * smk_539[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_450 = buffer.data(sll0 + 450);
    const auto *sll0_453 = buffer.data(sll0 + 453);
    const auto *sll0_456 = buffer.data(sll0 + 456);
    const auto *sll0_460 = buffer.data(sll0 + 460);
    const auto *sll0_462 = buffer.data(sll0 + 462);
    const auto *sll0_465 = buffer.data(sll0 + 465);
    const auto *sll0_467 = buffer.data(sll0 + 467);
    const auto *sll0_468 = buffer.data(sll0 + 468);
    const auto *sll0_471 = buffer.data(sll0 + 471);
    const auto *sll0_473 = buffer.data(sll0 + 473);
    const auto *sll0_474 = buffer.data(sll0 + 474);
    const auto *sll0_475 = buffer.data(sll0 + 475);
    const auto *sll0_486 = buffer.data(sll0 + 486);

    const auto *slk_359 = buffer.data(slk + 359);
    const auto *slk_360 = buffer.data(slk + 360);
    const auto *slk_362 = buffer.data(slk + 362);
    const auto *slk_363 = buffer.data(slk + 363);
    const auto *slk_365 = buffer.data(slk + 365);
    const auto *slk_366 = buffer.data(slk + 366);
    const auto *slk_367 = buffer.data(slk + 367);
    const auto *slk_369 = buffer.data(slk + 369);
    const auto *slk_370 = buffer.data(slk + 370);
    const auto *slk_371 = buffer.data(slk + 371);
    const auto *slk_372 = buffer.data(slk + 372);
    const auto *slk_374 = buffer.data(slk + 374);
    const auto *slk_375 = buffer.data(slk + 375);
    const auto *slk_376 = buffer.data(slk + 376);
    const auto *slk_377 = buffer.data(slk + 377);
    const auto *slk_378 = buffer.data(slk + 378);
    const auto *slk_380 = buffer.data(slk + 380);
    const auto *slk_388 = buffer.data(slk + 388);
    const auto *slk_390 = buffer.data(slk + 390);
    const auto *slk_391 = buffer.data(slk + 391);
    const auto *slk_392 = buffer.data(slk + 392);
    const auto *slk_393 = buffer.data(slk + 393);
    const auto *slk_394 = buffer.data(slk + 394);
    const auto *slk_395 = buffer.data(slk + 395);
    const auto *slk_396 = buffer.data(slk + 396);
    const auto *slk_398 = buffer.data(slk + 398);
    const auto *slk_399 = buffer.data(slk + 399);
    const auto *slk_401 = buffer.data(slk + 401);
    const auto *slk_402 = buffer.data(slk + 402);
    const auto *slk_405 = buffer.data(slk + 405);
    const auto *slk_406 = buffer.data(slk + 406);
    const auto *slk_410 = buffer.data(slk + 410);
    const auto *slk_416 = buffer.data(slk + 416);
    const auto *slk_426 = buffer.data(slk + 426);
    const auto *slk_427 = buffer.data(slk + 427);
    const auto *slk_428 = buffer.data(slk + 428);
    const auto *slk_429 = buffer.data(slk + 429);
    const auto *slk_430 = buffer.data(slk + 430);
    const auto *slk_431 = buffer.data(slk + 431);
    const auto *slk_432 = buffer.data(slk + 432);
    const auto *slk_434 = buffer.data(slk + 434);
    const auto *slk_437 = buffer.data(slk + 437);
    const auto *slk_441 = buffer.data(slk + 441);
    const auto *slk_540 = buffer.data(slk + 540);
    const auto *slk_543 = buffer.data(slk + 543);
    const auto *slk_545 = buffer.data(slk + 545);
    const auto *slk_546 = buffer.data(slk + 546);
    const auto *slk_549 = buffer.data(slk + 549);
    const auto *slk_550 = buffer.data(slk + 550);
    const auto *slk_552 = buffer.data(slk + 552);
    const auto *slk_554 = buffer.data(slk + 554);
    const auto *slk_555 = buffer.data(slk + 555);
    const auto *slk_557 = buffer.data(slk + 557);
    const auto *slk_558 = buffer.data(slk + 558);
    const auto *slk_560 = buffer.data(slk + 560);
    const auto *slk_561 = buffer.data(slk + 561);
    const auto *slk_563 = buffer.data(slk + 563);
    const auto *slk_564 = buffer.data(slk + 564);
    const auto *slk_565 = buffer.data(slk + 565);
    const auto *slk_567 = buffer.data(slk + 567);
    const auto *slk_568 = buffer.data(slk + 568);
    const auto *slk_569 = buffer.data(slk + 569);
    const auto *slk_570 = buffer.data(slk + 570);
    const auto *slk_571 = buffer.data(slk + 571);
    const auto *slk_572 = buffer.data(slk + 572);
    const auto *slk_573 = buffer.data(slk + 573);
    const auto *slk_574 = buffer.data(slk + 574);
    const auto *slk_575 = buffer.data(slk + 575);
    const auto *slk_581 = buffer.data(slk + 581);
    const auto *slk_585 = buffer.data(slk + 585);
    const auto *slk_590 = buffer.data(slk + 590);
    const auto *slk_596 = buffer.data(slk + 596);
    const auto *slk_603 = buffer.data(slk + 603);
    const auto *slk_604 = buffer.data(slk + 604);
    const auto *slk_605 = buffer.data(slk + 605);
    const auto *slk_606 = buffer.data(slk + 606);
    const auto *slk_607 = buffer.data(slk + 607);
    const auto *slk_608 = buffer.data(slk + 608);
    const auto *slk_609 = buffer.data(slk + 609);
    const auto *slk_610 = buffer.data(slk + 610);
    const auto *slk_611 = buffer.data(slk + 611);
    const auto *slk_612 = buffer.data(slk + 612);
    const auto *slk_615 = buffer.data(slk + 615);
    const auto *slk_617 = buffer.data(slk + 617);
    const auto *slk_618 = buffer.data(slk + 618);
    const auto *slk_621 = buffer.data(slk + 621);
    const auto *slk_622 = buffer.data(slk + 622);
    const auto *slk_624 = buffer.data(slk + 624);
    const auto *slk_626 = buffer.data(slk + 626);
    const auto *slk_627 = buffer.data(slk + 627);

    const auto *sll1_450 = buffer.data(sll1 + 450);
    const auto *sll1_453 = buffer.data(sll1 + 453);
    const auto *sll1_456 = buffer.data(sll1 + 456);
    const auto *sll1_460 = buffer.data(sll1 + 460);
    const auto *sll1_462 = buffer.data(sll1 + 462);
    const auto *sll1_465 = buffer.data(sll1 + 465);
    const auto *sll1_467 = buffer.data(sll1 + 467);
    const auto *sll1_468 = buffer.data(sll1 + 468);
    const auto *sll1_471 = buffer.data(sll1 + 471);
    const auto *sll1_473 = buffer.data(sll1 + 473);
    const auto *sll1_474 = buffer.data(sll1 + 474);
    const auto *sll1_475 = buffer.data(sll1 + 475);
    const auto *sll1_486 = buffer.data(sll1 + 486);

    const auto *smi0_419 = buffer.data(smi0 + 419);
    const auto *smi0_420 = buffer.data(smi0 + 420);
    const auto *smi0_423 = buffer.data(smi0 + 423);
    const auto *smi0_425 = buffer.data(smi0 + 425);
    const auto *smi0_426 = buffer.data(smi0 + 426);
    const auto *smi0_429 = buffer.data(smi0 + 429);
    const auto *smi0_430 = buffer.data(smi0 + 430);
    const auto *smi0_432 = buffer.data(smi0 + 432);
    const auto *smi0_434 = buffer.data(smi0 + 434);
    const auto *smi0_435 = buffer.data(smi0 + 435);
    const auto *smi0_437 = buffer.data(smi0 + 437);
    const auto *smi0_438 = buffer.data(smi0 + 438);
    const auto *smi0_440 = buffer.data(smi0 + 440);
    const auto *smi0_441 = buffer.data(smi0 + 441);
    const auto *smi0_443 = buffer.data(smi0 + 443);
    const auto *smi0_444 = buffer.data(smi0 + 444);
    const auto *smi0_445 = buffer.data(smi0 + 445);
    const auto *smi0_446 = buffer.data(smi0 + 446);
    const auto *smi0_447 = buffer.data(smi0 + 447);
    const auto *smi0_453 = buffer.data(smi0 + 453);
    const auto *smi0_457 = buffer.data(smi0 + 457);
    const auto *smi0_462 = buffer.data(smi0 + 462);
    const auto *smi0_468 = buffer.data(smi0 + 468);
    const auto *smi0_471 = buffer.data(smi0 + 471);
    const auto *smi0_472 = buffer.data(smi0 + 472);
    const auto *smi0_473 = buffer.data(smi0 + 473);
    const auto *smi0_474 = buffer.data(smi0 + 474);
    const auto *smi0_475 = buffer.data(smi0 + 475);
    const auto *smi0_476 = buffer.data(smi0 + 476);
    const auto *smi0_479 = buffer.data(smi0 + 479);
    const auto *smi0_481 = buffer.data(smi0 + 481);
    const auto *smi0_482 = buffer.data(smi0 + 482);
    const auto *smi0_485 = buffer.data(smi0 + 485);
    const auto *smi0_486 = buffer.data(smi0 + 486);
    const auto *smi0_488 = buffer.data(smi0 + 488);
    const auto *smi0_490 = buffer.data(smi0 + 490);
    const auto *smi0_491 = buffer.data(smi0 + 491);

    const auto *smi1_419 = buffer.data(smi1 + 419);
    const auto *smi1_420 = buffer.data(smi1 + 420);
    const auto *smi1_423 = buffer.data(smi1 + 423);
    const auto *smi1_425 = buffer.data(smi1 + 425);
    const auto *smi1_426 = buffer.data(smi1 + 426);
    const auto *smi1_429 = buffer.data(smi1 + 429);
    const auto *smi1_430 = buffer.data(smi1 + 430);
    const auto *smi1_432 = buffer.data(smi1 + 432);
    const auto *smi1_434 = buffer.data(smi1 + 434);
    const auto *smi1_435 = buffer.data(smi1 + 435);
    const auto *smi1_437 = buffer.data(smi1 + 437);
    const auto *smi1_438 = buffer.data(smi1 + 438);
    const auto *smi1_440 = buffer.data(smi1 + 440);
    const auto *smi1_441 = buffer.data(smi1 + 441);
    const auto *smi1_443 = buffer.data(smi1 + 443);
    const auto *smi1_444 = buffer.data(smi1 + 444);
    const auto *smi1_445 = buffer.data(smi1 + 445);
    const auto *smi1_446 = buffer.data(smi1 + 446);
    const auto *smi1_447 = buffer.data(smi1 + 447);
    const auto *smi1_453 = buffer.data(smi1 + 453);
    const auto *smi1_457 = buffer.data(smi1 + 457);
    const auto *smi1_462 = buffer.data(smi1 + 462);
    const auto *smi1_468 = buffer.data(smi1 + 468);
    const auto *smi1_471 = buffer.data(smi1 + 471);
    const auto *smi1_472 = buffer.data(smi1 + 472);
    const auto *smi1_473 = buffer.data(smi1 + 473);
    const auto *smi1_474 = buffer.data(smi1 + 474);
    const auto *smi1_475 = buffer.data(smi1 + 475);
    const auto *smi1_476 = buffer.data(smi1 + 476);
    const auto *smi1_479 = buffer.data(smi1 + 479);
    const auto *smi1_481 = buffer.data(smi1 + 481);
    const auto *smi1_482 = buffer.data(smi1 + 482);
    const auto *smi1_485 = buffer.data(smi1 + 485);
    const auto *smi1_486 = buffer.data(smi1 + 486);
    const auto *smi1_488 = buffer.data(smi1 + 488);
    const auto *smi1_490 = buffer.data(smi1 + 490);
    const auto *smi1_491 = buffer.data(smi1 + 491);

    const auto *smk_539 = buffer.data(smk + 539);
    const auto *smk_540 = buffer.data(smk + 540);
    const auto *smk_542 = buffer.data(smk + 542);
    const auto *smk_543 = buffer.data(smk + 543);
    const auto *smk_545 = buffer.data(smk + 545);
    const auto *smk_546 = buffer.data(smk + 546);
    const auto *smk_549 = buffer.data(smk + 549);
    const auto *smk_550 = buffer.data(smk + 550);
    const auto *smk_552 = buffer.data(smk + 552);
    const auto *smk_554 = buffer.data(smk + 554);
    const auto *smk_555 = buffer.data(smk + 555);
    const auto *smk_557 = buffer.data(smk + 557);
    const auto *smk_558 = buffer.data(smk + 558);
    const auto *smk_560 = buffer.data(smk + 560);
    const auto *smk_561 = buffer.data(smk + 561);
    const auto *smk_563 = buffer.data(smk + 563);
    const auto *smk_564 = buffer.data(smk + 564);
    const auto *smk_565 = buffer.data(smk + 565);
    const auto *smk_567 = buffer.data(smk + 567);
    const auto *smk_568 = buffer.data(smk + 568);
    const auto *smk_569 = buffer.data(smk + 569);
    const auto *smk_570 = buffer.data(smk + 570);
    const auto *smk_571 = buffer.data(smk + 571);
    const auto *smk_572 = buffer.data(smk + 572);
    const auto *smk_573 = buffer.data(smk + 573);
    const auto *smk_574 = buffer.data(smk + 574);
    const auto *smk_575 = buffer.data(smk + 575);
    const auto *smk_576 = buffer.data(smk + 576);
    const auto *smk_578 = buffer.data(smk + 578);
    const auto *smk_579 = buffer.data(smk + 579);
    const auto *smk_581 = buffer.data(smk + 581);
    const auto *smk_582 = buffer.data(smk + 582);
    const auto *smk_585 = buffer.data(smk + 585);
    const auto *smk_586 = buffer.data(smk + 586);
    const auto *smk_590 = buffer.data(smk + 590);
    const auto *smk_591 = buffer.data(smk + 591);
    const auto *smk_596 = buffer.data(smk + 596);
    const auto *smk_603 = buffer.data(smk + 603);
    const auto *smk_604 = buffer.data(smk + 604);
    const auto *smk_605 = buffer.data(smk + 605);
    const auto *smk_606 = buffer.data(smk + 606);
    const auto *smk_607 = buffer.data(smk + 607);
    const auto *smk_608 = buffer.data(smk + 608);
    const auto *smk_609 = buffer.data(smk + 609);
    const auto *smk_610 = buffer.data(smk + 610);
    const auto *smk_611 = buffer.data(smk + 611);
    const auto *smk_612 = buffer.data(smk + 612);
    const auto *smk_614 = buffer.data(smk + 614);
    const auto *smk_615 = buffer.data(smk + 615);
    const auto *smk_617 = buffer.data(smk + 617);
    const auto *smk_618 = buffer.data(smk + 618);
    const auto *smk_621 = buffer.data(smk + 621);
    const auto *smk_622 = buffer.data(smk + 622);
    const auto *smk_624 = buffer.data(smk + 624);
    const auto *smk_626 = buffer.data(smk + 626);
    const auto *smk_627 = buffer.data(smk + 627);

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pc_x, pc_y, pc_z, slk_359, slk_360, \
                         slk_540, smi0_419, smi0_420, smi1_419, smi1_420, smk_539, \
                         smk_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_18 * slk_359[k]
                   + f_1 * smi0_419[k]
                   - f_2 * smi1_419[k]
                   + f_3 * pc_z[k] * smk_539[k];

        t_675[k] = f_18 * slk_540[k]
                   + f_1 * smi0_420[k]
                   - f_2 * smi1_420[k]
                   + f_3 * pc_x[k] * smk_540[k];

        t_676[k] = f_19 * slk_360[k]
                   + f_3 * pc_y[k] * smk_540[k];

        t_677[k] = f_3 * pc_z[k] * smk_540[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pc_x, pc_y, slk_362, slk_543, slk_545, smi0_423, \
                         smi0_425, smi1_423, smi1_425, smk_542, smk_543, \
                         smk_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_18 * slk_543[k]
                   + f_4 * smi0_423[k]
                   - f_5 * smi1_423[k]
                   + f_3 * pc_x[k] * smk_543[k];

        t_679[k] = f_19 * slk_362[k]
                   + f_3 * pc_y[k] * smk_542[k];

        t_680[k] = f_18 * slk_545[k]
                   + f_4 * smi0_425[k]
                   - f_5 * smi1_425[k]
                   + f_3 * pc_x[k] * smk_545[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, pc_x, pc_y, pc_z, slk_365, slk_546, smi0_426, \
                         smi1_426, smk_543, smk_545, smk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_18 * slk_546[k]
                   + f_6 * smi0_426[k]
                   - f_7 * smi1_426[k]
                   + f_3 * pc_x[k] * smk_546[k];

        t_682[k] = f_3 * pc_z[k] * smk_543[k];

        t_683[k] = f_19 * slk_365[k]
                   + f_3 * pc_y[k] * smk_545[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, pc_x, pc_z, slk_549, slk_550, smi0_429, \
                         smi0_430, smi1_429, smi1_430, smk_546, smk_549, \
                         smk_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_18 * slk_549[k]
                   + f_6 * smi0_429[k]
                   - f_7 * smi1_429[k]
                   + f_3 * pc_x[k] * smk_549[k];

        t_685[k] = f_18 * slk_550[k]
                   + f_8 * smi0_430[k]
                   - f_9 * smi1_430[k]
                   + f_3 * pc_x[k] * smk_550[k];

        t_686[k] = f_3 * pc_z[k] * smk_546[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pc_x, pc_y, slk_369, slk_552, slk_554, smi0_432, \
                         smi0_434, smi1_432, smi1_434, smk_549, smk_552, \
                         smk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_18 * slk_552[k]
                   + f_8 * smi0_432[k]
                   - f_9 * smi1_432[k]
                   + f_3 * pc_x[k] * smk_552[k];

        t_688[k] = f_19 * slk_369[k]
                   + f_3 * pc_y[k] * smk_549[k];

        t_689[k] = f_18 * slk_554[k]
                   + f_8 * smi0_434[k]
                   - f_9 * smi1_434[k]
                   + f_3 * pc_x[k] * smk_554[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, pc_x, pc_z, slk_555, slk_557, smi0_435, \
                         smi0_437, smi1_435, smi1_437, smk_550, smk_555, \
                         smk_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_18 * slk_555[k]
                   + f_10 * smi0_435[k]
                   - f_11 * smi1_435[k]
                   + f_3 * pc_x[k] * smk_555[k];

        t_691[k] = f_3 * pc_z[k] * smk_550[k];

        t_692[k] = f_18 * slk_557[k]
                   + f_10 * smi0_437[k]
                   - f_11 * smi1_437[k]
                   + f_3 * pc_x[k] * smk_557[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pc_x, pc_y, slk_374, slk_558, slk_560, smi0_438, \
                         smi0_440, smi1_438, smi1_440, smk_554, smk_558, \
                         smk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_18 * slk_558[k]
                   + f_10 * smi0_438[k]
                   - f_11 * smi1_438[k]
                   + f_3 * pc_x[k] * smk_558[k];

        t_694[k] = f_19 * slk_374[k]
                   + f_3 * pc_y[k] * smk_554[k];

        t_695[k] = f_18 * slk_560[k]
                   + f_10 * smi0_440[k]
                   - f_11 * smi1_440[k]
                   + f_3 * pc_x[k] * smk_560[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pc_x, pc_z, slk_561, slk_563, smi0_441, \
                         smi0_443, smi1_441, smi1_443, smk_555, smk_561, \
                         smk_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_18 * slk_561[k]
                   + f_12 * smi0_441[k]
                   - f_13 * smi1_441[k]
                   + f_3 * pc_x[k] * smk_561[k];

        t_697[k] = f_3 * pc_z[k] * smk_555[k];

        t_698[k] = f_18 * slk_563[k]
                   + f_12 * smi0_443[k]
                   - f_13 * smi1_443[k]
                   + f_3 * pc_x[k] * smk_563[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, pc_x, pc_y, slk_380, slk_564, slk_565, smi0_444, \
                         smi0_445, smi1_444, smi1_445, smk_560, smk_564, \
                         smk_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_18 * slk_564[k]
                   + f_12 * smi0_444[k]
                   - f_13 * smi1_444[k]
                   + f_3 * pc_x[k] * smk_564[k];

        t_700[k] = f_18 * slk_565[k]
                   + f_12 * smi0_445[k]
                   - f_13 * smi1_445[k]
                   + f_3 * pc_x[k] * smk_565[k];

        t_701[k] = f_19 * slk_380[k]
                   + f_3 * pc_y[k] * smk_560[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, pc_x, slk_567, slk_568, slk_569, slk_570, \
                         smi0_447, smi1_447, smk_567, smk_568, smk_569, \
                         smk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = f_18 * slk_567[k]
                   + f_12 * smi0_447[k]
                   - f_13 * smi1_447[k]
                   + f_3 * pc_x[k] * smk_567[k];

        t_703[k] = f_18 * slk_568[k]
                   + f_3 * pc_x[k] * smk_568[k];

        t_704[k] = f_18 * slk_569[k]
                   + f_3 * pc_x[k] * smk_569[k];

        t_705[k] = f_18 * slk_570[k]
                   + f_3 * pc_x[k] * smk_570[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, t_710, pc_x, slk_571, slk_572, slk_573, \
                         slk_574, slk_575, smk_571, smk_572, smk_573, smk_574, \
                         smk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_18 * slk_571[k]
                   + f_3 * pc_x[k] * smk_571[k];

        t_707[k] = f_18 * slk_572[k]
                   + f_3 * pc_x[k] * smk_572[k];

        t_708[k] = f_18 * slk_573[k]
                   + f_3 * pc_x[k] * smk_573[k];

        t_709[k] = f_18 * slk_574[k]
                   + f_3 * pc_x[k] * smk_574[k];

        t_710[k] = f_18 * slk_575[k]
                   + f_3 * pc_x[k] * smk_575[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, pc_y, pc_z, slk_388, slk_390, smi0_441, \
                         smi0_443, smi1_441, smi1_443, smk_568, \
                         smk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_19 * slk_388[k]
                   + f_1 * smi0_441[k]
                   - f_2 * smi1_441[k]
                   + f_3 * pc_y[k] * smk_568[k];

        t_712[k] = f_3 * pc_z[k] * smk_568[k];

        t_713[k] = f_19 * slk_390[k]
                   + f_4 * smi0_443[k]
                   - f_5 * smi1_443[k]
                   + f_3 * pc_y[k] * smk_570[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, slk_391, slk_392, slk_393, smi0_444, \
                         smi0_445, smi0_446, smi1_444, smi1_445, smi1_446, smk_571, smk_572, \
                         smk_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_19 * slk_391[k]
                   + f_6 * smi0_444[k]
                   - f_7 * smi1_444[k]
                   + f_3 * pc_y[k] * smk_571[k];

        t_715[k] = f_19 * slk_392[k]
                   + f_8 * smi0_445[k]
                   - f_9 * smi1_445[k]
                   + f_3 * pc_y[k] * smk_572[k];

        t_716[k] = f_19 * slk_393[k]
                   + f_10 * smi0_446[k]
                   - f_11 * smi1_446[k]
                   + f_3 * pc_y[k] * smk_573[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, pb_z, pc_y, pc_z, sll0_450, slk_394, \
                         slk_395, sll1_450, smi0_447, smi1_447, smk_574, \
                         smk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_19 * slk_394[k]
                   + f_12 * smi0_447[k]
                   - f_13 * smi1_447[k]
                   + f_3 * pc_y[k] * smk_574[k];

        t_718[k] = f_19 * slk_395[k]
                   + f_3 * pc_y[k] * smk_575[k];

        t_719[k] = f_1 * smi0_447[k]
                   - f_2 * smi1_447[k]
                   + f_3 * pc_z[k] * smk_575[k];

        t_720[k] = pb_z[k] * sll0_450[k]
                   - f_14 * pc_z[k] * sll1_450[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pb_z, pc_y, pc_z, sll0_453, slk_360, \
                         slk_396, slk_398, sll1_453, smk_576, smk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_18 * slk_396[k]
                   + f_3 * pc_y[k] * smk_576[k];

        t_722[k] = f_15 * slk_360[k]
                   + f_3 * pc_z[k] * smk_576[k];

        t_723[k] = pb_z[k] * sll0_453[k]
                   - f_14 * pc_z[k] * sll1_453[k];

        t_724[k] = f_18 * slk_398[k]
                   + f_3 * pc_y[k] * smk_578[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, pb_z, pc_x, pc_z, sll0_456, slk_363, slk_581, \
                         sll1_456, smi0_453, smi1_453, smk_579, \
                         smk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_18 * slk_581[k]
                   + f_4 * smi0_453[k]
                   - f_5 * smi1_453[k]
                   + f_3 * pc_x[k] * smk_581[k];

        t_726[k] = pb_z[k] * sll0_456[k]
                   - f_14 * pc_z[k] * sll1_456[k];

        t_727[k] = f_15 * slk_363[k]
                   + f_3 * pc_z[k] * smk_579[k];
    }

#pragma omp simd aligned(t_728, t_729, t_730, pb_z, pc_x, pc_y, pc_z, sll0_460, slk_401, \
                         slk_585, sll1_460, smi0_457, smi1_457, smk_581, \
                         smk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_728[k] = f_18 * slk_401[k]
                   + f_3 * pc_y[k] * smk_581[k];

        t_729[k] = f_18 * slk_585[k]
                   + f_6 * smi0_457[k]
                   - f_7 * smi1_457[k]
                   + f_3 * pc_x[k] * smk_585[k];

        t_730[k] = pb_z[k] * sll0_460[k]
                   - f_14 * pc_z[k] * sll1_460[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, pb_z, pc_y, pc_z, sll0_462, slk_366, slk_367, \
                         slk_405, sll1_462, smk_582, smk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_15 * slk_366[k]
                   + f_3 * pc_z[k] * smk_582[k];

        t_732[k] = pb_z[k] * sll0_462[k]
                   + f_16 * slk_367[k]
                   - f_14 * pc_z[k] * sll1_462[k];

        t_733[k] = f_18 * slk_405[k]
                   + f_3 * pc_y[k] * smk_585[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pb_z, pc_x, pc_z, sll0_465, slk_370, slk_590, \
                         sll1_465, smi0_462, smi1_462, smk_586, \
                         smk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_18 * slk_590[k]
                   + f_8 * smi0_462[k]
                   - f_9 * smi1_462[k]
                   + f_3 * pc_x[k] * smk_590[k];

        t_735[k] = pb_z[k] * sll0_465[k]
                   - f_14 * pc_z[k] * sll1_465[k];

        t_736[k] = f_15 * slk_370[k]
                   + f_3 * pc_z[k] * smk_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pb_z, pc_y, pc_z, sll0_467, sll0_468, slk_371, \
                         slk_372, slk_410, sll1_467, sll1_468, \
                         smk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pb_z[k] * sll0_467[k]
                   + f_16 * slk_371[k]
                   - f_14 * pc_z[k] * sll1_467[k];

        t_738[k] = pb_z[k] * sll0_468[k]
                   + f_17 * slk_372[k]
                   - f_14 * pc_z[k] * sll1_468[k];

        t_739[k] = f_18 * slk_410[k]
                   + f_3 * pc_y[k] * smk_590[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pb_z, pc_x, pc_z, sll0_471, slk_375, slk_596, \
                         sll1_471, smi0_468, smi1_468, smk_591, \
                         smk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_18 * slk_596[k]
                   + f_10 * smi0_468[k]
                   - f_11 * smi1_468[k]
                   + f_3 * pc_x[k] * smk_596[k];

        t_741[k] = pb_z[k] * sll0_471[k]
                   - f_14 * pc_z[k] * sll1_471[k];

        t_742[k] = f_15 * slk_375[k]
                   + f_3 * pc_z[k] * smk_591[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, pb_z, pc_z, sll0_473, sll0_474, sll0_475, \
                         slk_376, slk_377, slk_378, sll1_473, sll1_474, \
                         sll1_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = pb_z[k] * sll0_473[k]
                   + f_16 * slk_376[k]
                   - f_14 * pc_z[k] * sll1_473[k];

        t_744[k] = pb_z[k] * sll0_474[k]
                   + f_17 * slk_377[k]
                   - f_14 * pc_z[k] * sll1_474[k];

        t_745[k] = pb_z[k] * sll0_475[k]
                   + f_18 * slk_378[k]
                   - f_14 * pc_z[k] * sll1_475[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, slk_416, slk_603, slk_604, \
                         slk_605, smi0_475, smi1_475, smk_596, smk_603, smk_604, \
                         smk_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_18 * slk_416[k]
                   + f_3 * pc_y[k] * smk_596[k];

        t_747[k] = f_18 * slk_603[k]
                   + f_12 * smi0_475[k]
                   - f_13 * smi1_475[k]
                   + f_3 * pc_x[k] * smk_603[k];

        t_748[k] = f_18 * slk_604[k]
                   + f_3 * pc_x[k] * smk_604[k];

        t_749[k] = f_18 * slk_605[k]
                   + f_3 * pc_x[k] * smk_605[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, pc_x, slk_606, slk_607, slk_608, \
                         slk_609, slk_610, smk_606, smk_607, smk_608, smk_609, \
                         smk_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_18 * slk_606[k]
                   + f_3 * pc_x[k] * smk_606[k];

        t_751[k] = f_18 * slk_607[k]
                   + f_3 * pc_x[k] * smk_607[k];

        t_752[k] = f_18 * slk_608[k]
                   + f_3 * pc_x[k] * smk_608[k];

        t_753[k] = f_18 * slk_609[k]
                   + f_3 * pc_x[k] * smk_609[k];

        t_754[k] = f_18 * slk_610[k]
                   + f_3 * pc_x[k] * smk_610[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, pb_z, pc_x, pc_z, sll0_486, slk_388, slk_611, \
                         sll1_486, smk_604, smk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_18 * slk_611[k]
                   + f_3 * pc_x[k] * smk_611[k];

        t_756[k] = pb_z[k] * sll0_486[k]
                   - f_14 * pc_z[k] * sll1_486[k];

        t_757[k] = f_15 * slk_388[k]
                   + f_3 * pc_z[k] * smk_604[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_y, slk_426, slk_427, slk_428, smi0_471, \
                         smi0_472, smi0_473, smi1_471, smi1_472, smi1_473, smk_606, smk_607, \
                         smk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_18 * slk_426[k]
                   + f_4 * smi0_471[k]
                   - f_5 * smi1_471[k]
                   + f_3 * pc_y[k] * smk_606[k];

        t_759[k] = f_18 * slk_427[k]
                   + f_6 * smi0_472[k]
                   - f_7 * smi1_472[k]
                   + f_3 * pc_y[k] * smk_607[k];

        t_760[k] = f_18 * slk_428[k]
                   + f_8 * smi0_473[k]
                   - f_9 * smi1_473[k]
                   + f_3 * pc_y[k] * smk_608[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_y, slk_429, slk_430, slk_431, smi0_474, \
                         smi0_475, smi1_474, smi1_475, smk_609, smk_610, \
                         smk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_18 * slk_429[k]
                   + f_10 * smi0_474[k]
                   - f_11 * smi1_474[k]
                   + f_3 * pc_y[k] * smk_609[k];

        t_762[k] = f_18 * slk_430[k]
                   + f_12 * smi0_475[k]
                   - f_13 * smi1_475[k]
                   + f_3 * pc_y[k] * smk_610[k];

        t_763[k] = f_18 * slk_431[k]
                   + f_3 * pc_y[k] * smk_611[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pc_x, pc_y, pc_z, slk_395, slk_432, slk_612, \
                         smi0_475, smi0_476, smi1_475, smi1_476, smk_611, \
                         smk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_15 * slk_395[k]
                   + f_1 * smi0_475[k]
                   - f_2 * smi1_475[k]
                   + f_3 * pc_z[k] * smk_611[k];

        t_765[k] = f_18 * slk_612[k]
                   + f_1 * smi0_476[k]
                   - f_2 * smi1_476[k]
                   + f_3 * pc_x[k] * smk_612[k];

        t_766[k] = f_17 * slk_432[k]
                   + f_3 * pc_y[k] * smk_612[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_x, pc_y, pc_z, slk_396, slk_434, slk_615, \
                         smi0_479, smi1_479, smk_612, smk_614, \
                         smk_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_16 * slk_396[k]
                   + f_3 * pc_z[k] * smk_612[k];

        t_768[k] = f_18 * slk_615[k]
                   + f_4 * smi0_479[k]
                   - f_5 * smi1_479[k]
                   + f_3 * pc_x[k] * smk_615[k];

        t_769[k] = f_17 * slk_434[k]
                   + f_3 * pc_y[k] * smk_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pc_x, pc_z, slk_399, slk_617, slk_618, smi0_481, \
                         smi0_482, smi1_481, smi1_482, smk_615, smk_617, \
                         smk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_18 * slk_617[k]
                   + f_4 * smi0_481[k]
                   - f_5 * smi1_481[k]
                   + f_3 * pc_x[k] * smk_617[k];

        t_771[k] = f_18 * slk_618[k]
                   + f_6 * smi0_482[k]
                   - f_7 * smi1_482[k]
                   + f_3 * pc_x[k] * smk_618[k];

        t_772[k] = f_16 * slk_399[k]
                   + f_3 * pc_z[k] * smk_615[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pc_x, pc_y, slk_437, slk_621, slk_622, smi0_485, \
                         smi0_486, smi1_485, smi1_486, smk_617, smk_621, \
                         smk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_17 * slk_437[k]
                   + f_3 * pc_y[k] * smk_617[k];

        t_774[k] = f_18 * slk_621[k]
                   + f_6 * smi0_485[k]
                   - f_7 * smi1_485[k]
                   + f_3 * pc_x[k] * smk_621[k];

        t_775[k] = f_18 * slk_622[k]
                   + f_8 * smi0_486[k]
                   - f_9 * smi1_486[k]
                   + f_3 * pc_x[k] * smk_622[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pc_x, pc_y, pc_z, slk_402, slk_441, slk_624, \
                         smi0_488, smi1_488, smk_618, smk_621, \
                         smk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_16 * slk_402[k]
                   + f_3 * pc_z[k] * smk_618[k];

        t_777[k] = f_18 * slk_624[k]
                   + f_8 * smi0_488[k]
                   - f_9 * smi1_488[k]
                   + f_3 * pc_x[k] * smk_624[k];

        t_778[k] = f_17 * slk_441[k]
                   + f_3 * pc_y[k] * smk_621[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pc_x, pc_z, slk_406, slk_626, slk_627, smi0_490, \
                         smi0_491, smi1_490, smi1_491, smk_622, smk_626, \
                         smk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_18 * slk_626[k]
                   + f_8 * smi0_490[k]
                   - f_9 * smi1_490[k]
                   + f_3 * pc_x[k] * smk_626[k];

        t_780[k] = f_18 * slk_627[k]
                   + f_10 * smi0_491[k]
                   - f_11 * smi1_491[k]
                   + f_3 * pc_x[k] * smk_627[k];

        t_781[k] = f_16 * slk_406[k]
                   + f_3 * pc_z[k] * smk_622[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_630 = buffer.data(sll0 + 630);
    const auto *sll0_633 = buffer.data(sll0 + 633);
    const auto *sll0_635 = buffer.data(sll0 + 635);
    const auto *sll0_636 = buffer.data(sll0 + 636);
    const auto *sll0_639 = buffer.data(sll0 + 639);
    const auto *sll0_640 = buffer.data(sll0 + 640);
    const auto *sll0_642 = buffer.data(sll0 + 642);
    const auto *sll0_644 = buffer.data(sll0 + 644);
    const auto *sll0_645 = buffer.data(sll0 + 645);
    const auto *sll0_647 = buffer.data(sll0 + 647);
    const auto *sll0_648 = buffer.data(sll0 + 648);
    const auto *sll0_650 = buffer.data(sll0 + 650);
    const auto *sll0_651 = buffer.data(sll0 + 651);
    const auto *sll0_653 = buffer.data(sll0 + 653);
    const auto *sll0_654 = buffer.data(sll0 + 654);
    const auto *sll0_655 = buffer.data(sll0 + 655);
    const auto *sll0_657 = buffer.data(sll0 + 657);

    const auto *slk_411 = buffer.data(slk + 411);
    const auto *slk_424 = buffer.data(slk + 424);
    const auto *slk_431 = buffer.data(slk + 431);
    const auto *slk_432 = buffer.data(slk + 432);
    const auto *slk_435 = buffer.data(slk + 435);
    const auto *slk_438 = buffer.data(slk + 438);
    const auto *slk_442 = buffer.data(slk + 442);
    const auto *slk_446 = buffer.data(slk + 446);
    const auto *slk_447 = buffer.data(slk + 447);
    const auto *slk_452 = buffer.data(slk + 452);
    const auto *slk_460 = buffer.data(slk + 460);
    const auto *slk_462 = buffer.data(slk + 462);
    const auto *slk_463 = buffer.data(slk + 463);
    const auto *slk_464 = buffer.data(slk + 464);
    const auto *slk_465 = buffer.data(slk + 465);
    const auto *slk_466 = buffer.data(slk + 466);
    const auto *slk_467 = buffer.data(slk + 467);
    const auto *slk_468 = buffer.data(slk + 468);
    const auto *slk_470 = buffer.data(slk + 470);
    const auto *slk_471 = buffer.data(slk + 471);
    const auto *slk_473 = buffer.data(slk + 473);
    const auto *slk_474 = buffer.data(slk + 474);
    const auto *slk_477 = buffer.data(slk + 477);
    const auto *slk_478 = buffer.data(slk + 478);
    const auto *slk_482 = buffer.data(slk + 482);
    const auto *slk_483 = buffer.data(slk + 483);
    const auto *slk_488 = buffer.data(slk + 488);
    const auto *slk_496 = buffer.data(slk + 496);
    const auto *slk_498 = buffer.data(slk + 498);
    const auto *slk_499 = buffer.data(slk + 499);
    const auto *slk_500 = buffer.data(slk + 500);
    const auto *slk_501 = buffer.data(slk + 501);
    const auto *slk_502 = buffer.data(slk + 502);
    const auto *slk_503 = buffer.data(slk + 503);
    const auto *slk_504 = buffer.data(slk + 504);
    const auto *slk_505 = buffer.data(slk + 505);
    const auto *slk_506 = buffer.data(slk + 506);
    const auto *slk_507 = buffer.data(slk + 507);
    const auto *slk_509 = buffer.data(slk + 509);
    const auto *slk_510 = buffer.data(slk + 510);
    const auto *slk_512 = buffer.data(slk + 512);
    const auto *slk_513 = buffer.data(slk + 513);
    const auto *slk_514 = buffer.data(slk + 514);
    const auto *slk_516 = buffer.data(slk + 516);
    const auto *slk_517 = buffer.data(slk + 517);
    const auto *slk_518 = buffer.data(slk + 518);
    const auto *slk_519 = buffer.data(slk + 519);
    const auto *slk_521 = buffer.data(slk + 521);
    const auto *slk_522 = buffer.data(slk + 522);
    const auto *slk_523 = buffer.data(slk + 523);
    const auto *slk_524 = buffer.data(slk + 524);
    const auto *slk_629 = buffer.data(slk + 629);
    const auto *slk_630 = buffer.data(slk + 630);
    const auto *slk_632 = buffer.data(slk + 632);
    const auto *slk_633 = buffer.data(slk + 633);
    const auto *slk_635 = buffer.data(slk + 635);
    const auto *slk_636 = buffer.data(slk + 636);
    const auto *slk_637 = buffer.data(slk + 637);
    const auto *slk_639 = buffer.data(slk + 639);
    const auto *slk_640 = buffer.data(slk + 640);
    const auto *slk_641 = buffer.data(slk + 641);
    const auto *slk_642 = buffer.data(slk + 642);
    const auto *slk_643 = buffer.data(slk + 643);
    const auto *slk_644 = buffer.data(slk + 644);
    const auto *slk_645 = buffer.data(slk + 645);
    const auto *slk_646 = buffer.data(slk + 646);
    const auto *slk_647 = buffer.data(slk + 647);
    const auto *slk_648 = buffer.data(slk + 648);
    const auto *slk_651 = buffer.data(slk + 651);
    const auto *slk_653 = buffer.data(slk + 653);
    const auto *slk_654 = buffer.data(slk + 654);
    const auto *slk_657 = buffer.data(slk + 657);
    const auto *slk_658 = buffer.data(slk + 658);
    const auto *slk_660 = buffer.data(slk + 660);
    const auto *slk_662 = buffer.data(slk + 662);
    const auto *slk_663 = buffer.data(slk + 663);
    const auto *slk_665 = buffer.data(slk + 665);
    const auto *slk_666 = buffer.data(slk + 666);
    const auto *slk_668 = buffer.data(slk + 668);
    const auto *slk_669 = buffer.data(slk + 669);
    const auto *slk_671 = buffer.data(slk + 671);
    const auto *slk_672 = buffer.data(slk + 672);
    const auto *slk_673 = buffer.data(slk + 673);
    const auto *slk_675 = buffer.data(slk + 675);
    const auto *slk_676 = buffer.data(slk + 676);
    const auto *slk_677 = buffer.data(slk + 677);
    const auto *slk_678 = buffer.data(slk + 678);
    const auto *slk_679 = buffer.data(slk + 679);
    const auto *slk_680 = buffer.data(slk + 680);
    const auto *slk_681 = buffer.data(slk + 681);
    const auto *slk_682 = buffer.data(slk + 682);
    const auto *slk_683 = buffer.data(slk + 683);
    const auto *slk_712 = buffer.data(slk + 712);
    const auto *slk_713 = buffer.data(slk + 713);
    const auto *slk_714 = buffer.data(slk + 714);
    const auto *slk_715 = buffer.data(slk + 715);
    const auto *slk_716 = buffer.data(slk + 716);
    const auto *slk_717 = buffer.data(slk + 717);

    const auto *sll1_630 = buffer.data(sll1 + 630);
    const auto *sll1_633 = buffer.data(sll1 + 633);
    const auto *sll1_635 = buffer.data(sll1 + 635);
    const auto *sll1_636 = buffer.data(sll1 + 636);
    const auto *sll1_639 = buffer.data(sll1 + 639);
    const auto *sll1_640 = buffer.data(sll1 + 640);
    const auto *sll1_642 = buffer.data(sll1 + 642);
    const auto *sll1_644 = buffer.data(sll1 + 644);
    const auto *sll1_645 = buffer.data(sll1 + 645);
    const auto *sll1_647 = buffer.data(sll1 + 647);
    const auto *sll1_648 = buffer.data(sll1 + 648);
    const auto *sll1_650 = buffer.data(sll1 + 650);
    const auto *sll1_651 = buffer.data(sll1 + 651);
    const auto *sll1_653 = buffer.data(sll1 + 653);
    const auto *sll1_654 = buffer.data(sll1 + 654);
    const auto *sll1_655 = buffer.data(sll1 + 655);
    const auto *sll1_657 = buffer.data(sll1 + 657);

    const auto *smi0_493 = buffer.data(smi0 + 493);
    const auto *smi0_494 = buffer.data(smi0 + 494);
    const auto *smi0_496 = buffer.data(smi0 + 496);
    const auto *smi0_497 = buffer.data(smi0 + 497);
    const auto *smi0_499 = buffer.data(smi0 + 499);
    const auto *smi0_500 = buffer.data(smi0 + 500);
    const auto *smi0_501 = buffer.data(smi0 + 501);
    const auto *smi0_502 = buffer.data(smi0 + 502);
    const auto *smi0_503 = buffer.data(smi0 + 503);
    const auto *smi0_504 = buffer.data(smi0 + 504);
    const auto *smi0_507 = buffer.data(smi0 + 507);
    const auto *smi0_509 = buffer.data(smi0 + 509);
    const auto *smi0_510 = buffer.data(smi0 + 510);
    const auto *smi0_513 = buffer.data(smi0 + 513);
    const auto *smi0_514 = buffer.data(smi0 + 514);
    const auto *smi0_516 = buffer.data(smi0 + 516);
    const auto *smi0_518 = buffer.data(smi0 + 518);
    const auto *smi0_519 = buffer.data(smi0 + 519);
    const auto *smi0_521 = buffer.data(smi0 + 521);
    const auto *smi0_522 = buffer.data(smi0 + 522);
    const auto *smi0_524 = buffer.data(smi0 + 524);
    const auto *smi0_525 = buffer.data(smi0 + 525);
    const auto *smi0_527 = buffer.data(smi0 + 527);
    const auto *smi0_528 = buffer.data(smi0 + 528);
    const auto *smi0_529 = buffer.data(smi0 + 529);
    const auto *smi0_530 = buffer.data(smi0 + 530);
    const auto *smi0_531 = buffer.data(smi0 + 531);

    const auto *smi1_493 = buffer.data(smi1 + 493);
    const auto *smi1_494 = buffer.data(smi1 + 494);
    const auto *smi1_496 = buffer.data(smi1 + 496);
    const auto *smi1_497 = buffer.data(smi1 + 497);
    const auto *smi1_499 = buffer.data(smi1 + 499);
    const auto *smi1_500 = buffer.data(smi1 + 500);
    const auto *smi1_501 = buffer.data(smi1 + 501);
    const auto *smi1_502 = buffer.data(smi1 + 502);
    const auto *smi1_503 = buffer.data(smi1 + 503);
    const auto *smi1_504 = buffer.data(smi1 + 504);
    const auto *smi1_507 = buffer.data(smi1 + 507);
    const auto *smi1_509 = buffer.data(smi1 + 509);
    const auto *smi1_510 = buffer.data(smi1 + 510);
    const auto *smi1_513 = buffer.data(smi1 + 513);
    const auto *smi1_514 = buffer.data(smi1 + 514);
    const auto *smi1_516 = buffer.data(smi1 + 516);
    const auto *smi1_518 = buffer.data(smi1 + 518);
    const auto *smi1_519 = buffer.data(smi1 + 519);
    const auto *smi1_521 = buffer.data(smi1 + 521);
    const auto *smi1_522 = buffer.data(smi1 + 522);
    const auto *smi1_524 = buffer.data(smi1 + 524);
    const auto *smi1_525 = buffer.data(smi1 + 525);
    const auto *smi1_527 = buffer.data(smi1 + 527);
    const auto *smi1_528 = buffer.data(smi1 + 528);
    const auto *smi1_529 = buffer.data(smi1 + 529);
    const auto *smi1_530 = buffer.data(smi1 + 530);
    const auto *smi1_531 = buffer.data(smi1 + 531);

    const auto *smk_626 = buffer.data(smk + 626);
    const auto *smk_627 = buffer.data(smk + 627);
    const auto *smk_629 = buffer.data(smk + 629);
    const auto *smk_630 = buffer.data(smk + 630);
    const auto *smk_632 = buffer.data(smk + 632);
    const auto *smk_633 = buffer.data(smk + 633);
    const auto *smk_635 = buffer.data(smk + 635);
    const auto *smk_636 = buffer.data(smk + 636);
    const auto *smk_637 = buffer.data(smk + 637);
    const auto *smk_639 = buffer.data(smk + 639);
    const auto *smk_640 = buffer.data(smk + 640);
    const auto *smk_641 = buffer.data(smk + 641);
    const auto *smk_642 = buffer.data(smk + 642);
    const auto *smk_643 = buffer.data(smk + 643);
    const auto *smk_644 = buffer.data(smk + 644);
    const auto *smk_645 = buffer.data(smk + 645);
    const auto *smk_646 = buffer.data(smk + 646);
    const auto *smk_647 = buffer.data(smk + 647);
    const auto *smk_648 = buffer.data(smk + 648);
    const auto *smk_650 = buffer.data(smk + 650);
    const auto *smk_651 = buffer.data(smk + 651);
    const auto *smk_653 = buffer.data(smk + 653);
    const auto *smk_654 = buffer.data(smk + 654);
    const auto *smk_657 = buffer.data(smk + 657);
    const auto *smk_658 = buffer.data(smk + 658);
    const auto *smk_660 = buffer.data(smk + 660);
    const auto *smk_662 = buffer.data(smk + 662);
    const auto *smk_663 = buffer.data(smk + 663);
    const auto *smk_665 = buffer.data(smk + 665);
    const auto *smk_666 = buffer.data(smk + 666);
    const auto *smk_668 = buffer.data(smk + 668);
    const auto *smk_669 = buffer.data(smk + 669);
    const auto *smk_671 = buffer.data(smk + 671);
    const auto *smk_672 = buffer.data(smk + 672);
    const auto *smk_673 = buffer.data(smk + 673);
    const auto *smk_675 = buffer.data(smk + 675);
    const auto *smk_676 = buffer.data(smk + 676);
    const auto *smk_677 = buffer.data(smk + 677);
    const auto *smk_678 = buffer.data(smk + 678);
    const auto *smk_679 = buffer.data(smk + 679);
    const auto *smk_680 = buffer.data(smk + 680);
    const auto *smk_681 = buffer.data(smk + 681);
    const auto *smk_682 = buffer.data(smk + 682);
    const auto *smk_683 = buffer.data(smk + 683);
    const auto *smk_684 = buffer.data(smk + 684);
    const auto *smk_686 = buffer.data(smk + 686);
    const auto *smk_687 = buffer.data(smk + 687);
    const auto *smk_689 = buffer.data(smk + 689);
    const auto *smk_690 = buffer.data(smk + 690);
    const auto *smk_693 = buffer.data(smk + 693);
    const auto *smk_694 = buffer.data(smk + 694);
    const auto *smk_698 = buffer.data(smk + 698);
    const auto *smk_699 = buffer.data(smk + 699);
    const auto *smk_704 = buffer.data(smk + 704);
    const auto *smk_712 = buffer.data(smk + 712);
    const auto *smk_713 = buffer.data(smk + 713);
    const auto *smk_714 = buffer.data(smk + 714);
    const auto *smk_715 = buffer.data(smk + 715);
    const auto *smk_716 = buffer.data(smk + 716);
    const auto *smk_717 = buffer.data(smk + 717);

#pragma omp simd aligned(t_782, t_783, t_784, pc_x, pc_y, slk_446, slk_629, slk_630, smi0_493, \
                         smi0_494, smi1_493, smi1_494, smk_626, smk_629, \
                         smk_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_18 * slk_629[k]
                   + f_10 * smi0_493[k]
                   - f_11 * smi1_493[k]
                   + f_3 * pc_x[k] * smk_629[k];

        t_783[k] = f_18 * slk_630[k]
                   + f_10 * smi0_494[k]
                   - f_11 * smi1_494[k]
                   + f_3 * pc_x[k] * smk_630[k];

        t_784[k] = f_17 * slk_446[k]
                   + f_3 * pc_y[k] * smk_626[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pc_x, pc_z, slk_411, slk_632, slk_633, smi0_496, \
                         smi0_497, smi1_496, smi1_497, smk_627, smk_632, \
                         smk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_18 * slk_632[k]
                   + f_10 * smi0_496[k]
                   - f_11 * smi1_496[k]
                   + f_3 * pc_x[k] * smk_632[k];

        t_786[k] = f_18 * slk_633[k]
                   + f_12 * smi0_497[k]
                   - f_13 * smi1_497[k]
                   + f_3 * pc_x[k] * smk_633[k];

        t_787[k] = f_16 * slk_411[k]
                   + f_3 * pc_z[k] * smk_627[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, slk_635, slk_636, slk_637, smi0_499, \
                         smi0_500, smi0_501, smi1_499, smi1_500, smi1_501, smk_635, smk_636, \
                         smk_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_18 * slk_635[k]
                   + f_12 * smi0_499[k]
                   - f_13 * smi1_499[k]
                   + f_3 * pc_x[k] * smk_635[k];

        t_789[k] = f_18 * slk_636[k]
                   + f_12 * smi0_500[k]
                   - f_13 * smi1_500[k]
                   + f_3 * pc_x[k] * smk_636[k];

        t_790[k] = f_18 * slk_637[k]
                   + f_12 * smi0_501[k]
                   - f_13 * smi1_501[k]
                   + f_3 * pc_x[k] * smk_637[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pc_x, pc_y, slk_452, slk_639, slk_640, \
                         slk_641, smi0_503, smi1_503, smk_632, smk_639, smk_640, \
                         smk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_17 * slk_452[k]
                   + f_3 * pc_y[k] * smk_632[k];

        t_792[k] = f_18 * slk_639[k]
                   + f_12 * smi0_503[k]
                   - f_13 * smi1_503[k]
                   + f_3 * pc_x[k] * smk_639[k];

        t_793[k] = f_18 * slk_640[k]
                   + f_3 * pc_x[k] * smk_640[k];

        t_794[k] = f_18 * slk_641[k]
                   + f_3 * pc_x[k] * smk_641[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, pc_x, slk_642, slk_643, slk_644, \
                         slk_645, slk_646, smk_642, smk_643, smk_644, smk_645, \
                         smk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_18 * slk_642[k]
                   + f_3 * pc_x[k] * smk_642[k];

        t_796[k] = f_18 * slk_643[k]
                   + f_3 * pc_x[k] * smk_643[k];

        t_797[k] = f_18 * slk_644[k]
                   + f_3 * pc_x[k] * smk_644[k];

        t_798[k] = f_18 * slk_645[k]
                   + f_3 * pc_x[k] * smk_645[k];

        t_799[k] = f_18 * slk_646[k]
                   + f_3 * pc_x[k] * smk_646[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, slk_424, slk_460, slk_647, \
                         smi0_497, smi1_497, smk_640, smk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_18 * slk_647[k]
                   + f_3 * pc_x[k] * smk_647[k];

        t_801[k] = f_17 * slk_460[k]
                   + f_1 * smi0_497[k]
                   - f_2 * smi1_497[k]
                   + f_3 * pc_y[k] * smk_640[k];

        t_802[k] = f_16 * slk_424[k]
                   + f_3 * pc_z[k] * smk_640[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_y, slk_462, slk_463, slk_464, smi0_499, \
                         smi0_500, smi0_501, smi1_499, smi1_500, smi1_501, smk_642, smk_643, \
                         smk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_17 * slk_462[k]
                   + f_4 * smi0_499[k]
                   - f_5 * smi1_499[k]
                   + f_3 * pc_y[k] * smk_642[k];

        t_804[k] = f_17 * slk_463[k]
                   + f_6 * smi0_500[k]
                   - f_7 * smi1_500[k]
                   + f_3 * pc_y[k] * smk_643[k];

        t_805[k] = f_17 * slk_464[k]
                   + f_8 * smi0_501[k]
                   - f_9 * smi1_501[k]
                   + f_3 * pc_y[k] * smk_644[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pc_y, slk_465, slk_466, slk_467, smi0_502, \
                         smi0_503, smi1_502, smi1_503, smk_645, smk_646, \
                         smk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * slk_465[k]
                   + f_10 * smi0_502[k]
                   - f_11 * smi1_502[k]
                   + f_3 * pc_y[k] * smk_645[k];

        t_807[k] = f_17 * slk_466[k]
                   + f_12 * smi0_503[k]
                   - f_13 * smi1_503[k]
                   + f_3 * pc_y[k] * smk_646[k];

        t_808[k] = f_17 * slk_467[k]
                   + f_3 * pc_y[k] * smk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pc_x, pc_y, pc_z, slk_431, slk_468, slk_648, \
                         smi0_503, smi0_504, smi1_503, smi1_504, smk_647, \
                         smk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_16 * slk_431[k]
                   + f_1 * smi0_503[k]
                   - f_2 * smi1_503[k]
                   + f_3 * pc_z[k] * smk_647[k];

        t_810[k] = f_18 * slk_648[k]
                   + f_1 * smi0_504[k]
                   - f_2 * smi1_504[k]
                   + f_3 * pc_x[k] * smk_648[k];

        t_811[k] = f_16 * slk_468[k]
                   + f_3 * pc_y[k] * smk_648[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pc_x, pc_y, pc_z, slk_432, slk_470, slk_651, \
                         smi0_507, smi1_507, smk_648, smk_650, \
                         smk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_17 * slk_432[k]
                   + f_3 * pc_z[k] * smk_648[k];

        t_813[k] = f_18 * slk_651[k]
                   + f_4 * smi0_507[k]
                   - f_5 * smi1_507[k]
                   + f_3 * pc_x[k] * smk_651[k];

        t_814[k] = f_16 * slk_470[k]
                   + f_3 * pc_y[k] * smk_650[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pc_x, pc_z, slk_435, slk_653, slk_654, smi0_509, \
                         smi0_510, smi1_509, smi1_510, smk_651, smk_653, \
                         smk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_18 * slk_653[k]
                   + f_4 * smi0_509[k]
                   - f_5 * smi1_509[k]
                   + f_3 * pc_x[k] * smk_653[k];

        t_816[k] = f_18 * slk_654[k]
                   + f_6 * smi0_510[k]
                   - f_7 * smi1_510[k]
                   + f_3 * pc_x[k] * smk_654[k];

        t_817[k] = f_17 * slk_435[k]
                   + f_3 * pc_z[k] * smk_651[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pc_x, pc_y, slk_473, slk_657, slk_658, smi0_513, \
                         smi0_514, smi1_513, smi1_514, smk_653, smk_657, \
                         smk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_16 * slk_473[k]
                   + f_3 * pc_y[k] * smk_653[k];

        t_819[k] = f_18 * slk_657[k]
                   + f_6 * smi0_513[k]
                   - f_7 * smi1_513[k]
                   + f_3 * pc_x[k] * smk_657[k];

        t_820[k] = f_18 * slk_658[k]
                   + f_8 * smi0_514[k]
                   - f_9 * smi1_514[k]
                   + f_3 * pc_x[k] * smk_658[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pc_x, pc_y, pc_z, slk_438, slk_477, slk_660, \
                         smi0_516, smi1_516, smk_654, smk_657, \
                         smk_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_17 * slk_438[k]
                   + f_3 * pc_z[k] * smk_654[k];

        t_822[k] = f_18 * slk_660[k]
                   + f_8 * smi0_516[k]
                   - f_9 * smi1_516[k]
                   + f_3 * pc_x[k] * smk_660[k];

        t_823[k] = f_16 * slk_477[k]
                   + f_3 * pc_y[k] * smk_657[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, slk_442, slk_662, slk_663, smi0_518, \
                         smi0_519, smi1_518, smi1_519, smk_658, smk_662, \
                         smk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_18 * slk_662[k]
                   + f_8 * smi0_518[k]
                   - f_9 * smi1_518[k]
                   + f_3 * pc_x[k] * smk_662[k];

        t_825[k] = f_18 * slk_663[k]
                   + f_10 * smi0_519[k]
                   - f_11 * smi1_519[k]
                   + f_3 * pc_x[k] * smk_663[k];

        t_826[k] = f_17 * slk_442[k]
                   + f_3 * pc_z[k] * smk_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, pc_x, pc_y, slk_482, slk_665, slk_666, smi0_521, \
                         smi0_522, smi1_521, smi1_522, smk_662, smk_665, \
                         smk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_18 * slk_665[k]
                   + f_10 * smi0_521[k]
                   - f_11 * smi1_521[k]
                   + f_3 * pc_x[k] * smk_665[k];

        t_828[k] = f_18 * slk_666[k]
                   + f_10 * smi0_522[k]
                   - f_11 * smi1_522[k]
                   + f_3 * pc_x[k] * smk_666[k];

        t_829[k] = f_16 * slk_482[k]
                   + f_3 * pc_y[k] * smk_662[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, pc_x, pc_z, slk_447, slk_668, slk_669, smi0_524, \
                         smi0_525, smi1_524, smi1_525, smk_663, smk_668, \
                         smk_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_18 * slk_668[k]
                   + f_10 * smi0_524[k]
                   - f_11 * smi1_524[k]
                   + f_3 * pc_x[k] * smk_668[k];

        t_831[k] = f_18 * slk_669[k]
                   + f_12 * smi0_525[k]
                   - f_13 * smi1_525[k]
                   + f_3 * pc_x[k] * smk_669[k];

        t_832[k] = f_17 * slk_447[k]
                   + f_3 * pc_z[k] * smk_663[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, slk_671, slk_672, slk_673, smi0_527, \
                         smi0_528, smi0_529, smi1_527, smi1_528, smi1_529, smk_671, smk_672, \
                         smk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_18 * slk_671[k]
                   + f_12 * smi0_527[k]
                   - f_13 * smi1_527[k]
                   + f_3 * pc_x[k] * smk_671[k];

        t_834[k] = f_18 * slk_672[k]
                   + f_12 * smi0_528[k]
                   - f_13 * smi1_528[k]
                   + f_3 * pc_x[k] * smk_672[k];

        t_835[k] = f_18 * slk_673[k]
                   + f_12 * smi0_529[k]
                   - f_13 * smi1_529[k]
                   + f_3 * pc_x[k] * smk_673[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pc_x, pc_y, slk_488, slk_675, slk_676, \
                         slk_677, smi0_531, smi1_531, smk_668, smk_675, smk_676, \
                         smk_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * slk_488[k]
                   + f_3 * pc_y[k] * smk_668[k];

        t_837[k] = f_18 * slk_675[k]
                   + f_12 * smi0_531[k]
                   - f_13 * smi1_531[k]
                   + f_3 * pc_x[k] * smk_675[k];

        t_838[k] = f_18 * slk_676[k]
                   + f_3 * pc_x[k] * smk_676[k];

        t_839[k] = f_18 * slk_677[k]
                   + f_3 * pc_x[k] * smk_677[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, pc_x, slk_678, slk_679, slk_680, \
                         slk_681, slk_682, smk_678, smk_679, smk_680, smk_681, \
                         smk_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_18 * slk_678[k]
                   + f_3 * pc_x[k] * smk_678[k];

        t_841[k] = f_18 * slk_679[k]
                   + f_3 * pc_x[k] * smk_679[k];

        t_842[k] = f_18 * slk_680[k]
                   + f_3 * pc_x[k] * smk_680[k];

        t_843[k] = f_18 * slk_681[k]
                   + f_3 * pc_x[k] * smk_681[k];

        t_844[k] = f_18 * slk_682[k]
                   + f_3 * pc_x[k] * smk_682[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, pc_z, slk_460, slk_496, slk_683, \
                         smi0_525, smi1_525, smk_676, smk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_18 * slk_683[k]
                   + f_3 * pc_x[k] * smk_683[k];

        t_846[k] = f_16 * slk_496[k]
                   + f_1 * smi0_525[k]
                   - f_2 * smi1_525[k]
                   + f_3 * pc_y[k] * smk_676[k];

        t_847[k] = f_17 * slk_460[k]
                   + f_3 * pc_z[k] * smk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pc_y, slk_498, slk_499, slk_500, smi0_527, \
                         smi0_528, smi0_529, smi1_527, smi1_528, smi1_529, smk_678, smk_679, \
                         smk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * slk_498[k]
                   + f_4 * smi0_527[k]
                   - f_5 * smi1_527[k]
                   + f_3 * pc_y[k] * smk_678[k];

        t_849[k] = f_16 * slk_499[k]
                   + f_6 * smi0_528[k]
                   - f_7 * smi1_528[k]
                   + f_3 * pc_y[k] * smk_679[k];

        t_850[k] = f_16 * slk_500[k]
                   + f_8 * smi0_529[k]
                   - f_9 * smi1_529[k]
                   + f_3 * pc_y[k] * smk_680[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pc_y, slk_501, slk_502, slk_503, smi0_530, \
                         smi0_531, smi1_530, smi1_531, smk_681, smk_682, \
                         smk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_16 * slk_501[k]
                   + f_10 * smi0_530[k]
                   - f_11 * smi1_530[k]
                   + f_3 * pc_y[k] * smk_681[k];

        t_852[k] = f_16 * slk_502[k]
                   + f_12 * smi0_531[k]
                   - f_13 * smi1_531[k]
                   + f_3 * pc_y[k] * smk_682[k];

        t_853[k] = f_16 * slk_503[k]
                   + f_3 * pc_y[k] * smk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pb_y, pc_y, pc_z, sll0_630, slk_467, \
                         slk_468, slk_504, sll1_630, smi0_531, smi1_531, smk_683, \
                         smk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_17 * slk_467[k]
                   + f_1 * smi0_531[k]
                   - f_2 * smi1_531[k]
                   + f_3 * pc_z[k] * smk_683[k];

        t_855[k] = pb_y[k] * sll0_630[k]
                   - f_14 * pc_y[k] * sll1_630[k];

        t_856[k] = f_15 * slk_504[k]
                   + f_3 * pc_y[k] * smk_684[k];

        t_857[k] = f_18 * slk_468[k]
                   + f_3 * pc_z[k] * smk_684[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pb_y, pc_y, sll0_633, sll0_635, sll0_636, \
                         slk_505, slk_506, slk_507, sll1_633, sll1_635, sll1_636, \
                         smk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pb_y[k] * sll0_633[k]
                   + f_16 * slk_505[k]
                   - f_14 * pc_y[k] * sll1_633[k];

        t_859[k] = f_15 * slk_506[k]
                   + f_3 * pc_y[k] * smk_686[k];

        t_860[k] = pb_y[k] * sll0_635[k]
                   - f_14 * pc_y[k] * sll1_635[k];

        t_861[k] = pb_y[k] * sll0_636[k]
                   + f_17 * slk_507[k]
                   - f_14 * pc_y[k] * sll1_636[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pb_y, pc_y, pc_z, sll0_639, sll0_640, \
                         slk_471, slk_509, slk_510, sll1_639, sll1_640, smk_687, \
                         smk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_18 * slk_471[k]
                   + f_3 * pc_z[k] * smk_687[k];

        t_863[k] = f_15 * slk_509[k]
                   + f_3 * pc_y[k] * smk_689[k];

        t_864[k] = pb_y[k] * sll0_639[k]
                   - f_14 * pc_y[k] * sll1_639[k];

        t_865[k] = pb_y[k] * sll0_640[k]
                   + f_18 * slk_510[k]
                   - f_14 * pc_y[k] * sll1_640[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pb_y, pc_y, pc_z, sll0_642, sll0_644, \
                         slk_474, slk_512, slk_513, sll1_642, sll1_644, smk_690, \
                         smk_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * slk_474[k]
                   + f_3 * pc_z[k] * smk_690[k];

        t_867[k] = pb_y[k] * sll0_642[k]
                   + f_16 * slk_512[k]
                   - f_14 * pc_y[k] * sll1_642[k];

        t_868[k] = f_15 * slk_513[k]
                   + f_3 * pc_y[k] * smk_693[k];

        t_869[k] = pb_y[k] * sll0_644[k]
                   - f_14 * pc_y[k] * sll1_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pb_y, pc_y, pc_z, sll0_645, sll0_647, slk_478, \
                         slk_514, slk_516, sll1_645, sll1_647, \
                         smk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pb_y[k] * sll0_645[k]
                   + f_19 * slk_514[k]
                   - f_14 * pc_y[k] * sll1_645[k];

        t_871[k] = f_18 * slk_478[k]
                   + f_3 * pc_z[k] * smk_694[k];

        t_872[k] = pb_y[k] * sll0_647[k]
                   + f_17 * slk_516[k]
                   - f_14 * pc_y[k] * sll1_647[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pb_y, pc_y, sll0_648, sll0_650, sll0_651, \
                         slk_517, slk_518, slk_519, sll1_648, sll1_650, sll1_651, \
                         smk_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pb_y[k] * sll0_648[k]
                   + f_16 * slk_517[k]
                   - f_14 * pc_y[k] * sll1_648[k];

        t_874[k] = f_15 * slk_518[k]
                   + f_3 * pc_y[k] * smk_698[k];

        t_875[k] = pb_y[k] * sll0_650[k]
                   - f_14 * pc_y[k] * sll1_650[k];

        t_876[k] = pb_y[k] * sll0_651[k]
                   + f_20 * slk_519[k]
                   - f_14 * pc_y[k] * sll1_651[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pb_y, pc_y, pc_z, sll0_653, sll0_654, slk_483, \
                         slk_521, slk_522, sll1_653, sll1_654, \
                         smk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_18 * slk_483[k]
                   + f_3 * pc_z[k] * smk_699[k];

        t_878[k] = pb_y[k] * sll0_653[k]
                   + f_18 * slk_521[k]
                   - f_14 * pc_y[k] * sll1_653[k];

        t_879[k] = pb_y[k] * sll0_654[k]
                   + f_17 * slk_522[k]
                   - f_14 * pc_y[k] * sll1_654[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pb_y, pc_x, pc_y, sll0_655, sll0_657, \
                         slk_523, slk_524, slk_712, sll1_655, sll1_657, smk_704, \
                         smk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = pb_y[k] * sll0_655[k]
                   + f_16 * slk_523[k]
                   - f_14 * pc_y[k] * sll1_655[k];

        t_881[k] = f_15 * slk_524[k]
                   + f_3 * pc_y[k] * smk_704[k];

        t_882[k] = pb_y[k] * sll0_657[k]
                   - f_14 * pc_y[k] * sll1_657[k];

        t_883[k] = f_18 * slk_712[k]
                   + f_3 * pc_x[k] * smk_712[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, pc_x, slk_713, slk_714, slk_715, \
                         slk_716, slk_717, smk_713, smk_714, smk_715, smk_716, \
                         smk_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_18 * slk_713[k]
                   + f_3 * pc_x[k] * smk_713[k];

        t_885[k] = f_18 * slk_714[k]
                   + f_3 * pc_x[k] * smk_714[k];

        t_886[k] = f_18 * slk_715[k]
                   + f_3 * pc_x[k] * smk_715[k];

        t_887[k] = f_18 * slk_716[k]
                   + f_3 * pc_x[k] * smk_716[k];

        t_888[k] = f_18 * slk_717[k]
                   + f_3 * pc_x[k] * smk_717[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_674 = buffer.data(sll0 + 674);
    const auto *sll0_675 = buffer.data(sll0 + 675);
    const auto *sll0_678 = buffer.data(sll0 + 678);
    const auto *sll0_681 = buffer.data(sll0 + 681);

    const auto *slk_496 = buffer.data(slk + 496);
    const auto *slk_504 = buffer.data(slk + 504);
    const auto *slk_507 = buffer.data(slk + 507);
    const auto *slk_510 = buffer.data(slk + 510);
    const auto *slk_514 = buffer.data(slk + 514);
    const auto *slk_519 = buffer.data(slk + 519);
    const auto *slk_532 = buffer.data(slk + 532);
    const auto *slk_534 = buffer.data(slk + 534);
    const auto *slk_535 = buffer.data(slk + 535);
    const auto *slk_536 = buffer.data(slk + 536);
    const auto *slk_537 = buffer.data(slk + 537);
    const auto *slk_538 = buffer.data(slk + 538);
    const auto *slk_539 = buffer.data(slk + 539);
    const auto *slk_540 = buffer.data(slk + 540);
    const auto *slk_542 = buffer.data(slk + 542);
    const auto *slk_543 = buffer.data(slk + 543);
    const auto *slk_545 = buffer.data(slk + 545);
    const auto *slk_549 = buffer.data(slk + 549);
    const auto *slk_554 = buffer.data(slk + 554);
    const auto *slk_560 = buffer.data(slk + 560);
    const auto *slk_568 = buffer.data(slk + 568);
    const auto *slk_570 = buffer.data(slk + 570);
    const auto *slk_571 = buffer.data(slk + 571);
    const auto *slk_572 = buffer.data(slk + 572);
    const auto *slk_573 = buffer.data(slk + 573);
    const auto *slk_574 = buffer.data(slk + 574);
    const auto *slk_575 = buffer.data(slk + 575);
    const auto *slk_576 = buffer.data(slk + 576);
    const auto *slk_578 = buffer.data(slk + 578);
    const auto *slk_718 = buffer.data(slk + 718);
    const auto *slk_719 = buffer.data(slk + 719);
    const auto *slk_720 = buffer.data(slk + 720);
    const auto *slk_723 = buffer.data(slk + 723);
    const auto *slk_725 = buffer.data(slk + 725);
    const auto *slk_726 = buffer.data(slk + 726);
    const auto *slk_729 = buffer.data(slk + 729);
    const auto *slk_730 = buffer.data(slk + 730);
    const auto *slk_732 = buffer.data(slk + 732);
    const auto *slk_734 = buffer.data(slk + 734);
    const auto *slk_735 = buffer.data(slk + 735);
    const auto *slk_737 = buffer.data(slk + 737);
    const auto *slk_738 = buffer.data(slk + 738);
    const auto *slk_740 = buffer.data(slk + 740);
    const auto *slk_741 = buffer.data(slk + 741);
    const auto *slk_743 = buffer.data(slk + 743);
    const auto *slk_744 = buffer.data(slk + 744);
    const auto *slk_745 = buffer.data(slk + 745);
    const auto *slk_747 = buffer.data(slk + 747);
    const auto *slk_748 = buffer.data(slk + 748);
    const auto *slk_749 = buffer.data(slk + 749);
    const auto *slk_750 = buffer.data(slk + 750);
    const auto *slk_751 = buffer.data(slk + 751);
    const auto *slk_752 = buffer.data(slk + 752);
    const auto *slk_753 = buffer.data(slk + 753);
    const auto *slk_754 = buffer.data(slk + 754);
    const auto *slk_755 = buffer.data(slk + 755);
    const auto *slk_756 = buffer.data(slk + 756);
    const auto *slk_759 = buffer.data(slk + 759);
    const auto *slk_761 = buffer.data(slk + 761);
    const auto *slk_762 = buffer.data(slk + 762);
    const auto *slk_765 = buffer.data(slk + 765);
    const auto *slk_766 = buffer.data(slk + 766);
    const auto *slk_768 = buffer.data(slk + 768);
    const auto *slk_770 = buffer.data(slk + 770);
    const auto *slk_771 = buffer.data(slk + 771);
    const auto *slk_773 = buffer.data(slk + 773);
    const auto *slk_774 = buffer.data(slk + 774);
    const auto *slk_776 = buffer.data(slk + 776);
    const auto *slk_777 = buffer.data(slk + 777);
    const auto *slk_779 = buffer.data(slk + 779);
    const auto *slk_780 = buffer.data(slk + 780);
    const auto *slk_781 = buffer.data(slk + 781);
    const auto *slk_783 = buffer.data(slk + 783);
    const auto *slk_784 = buffer.data(slk + 784);
    const auto *slk_785 = buffer.data(slk + 785);
    const auto *slk_786 = buffer.data(slk + 786);
    const auto *slk_787 = buffer.data(slk + 787);
    const auto *slk_788 = buffer.data(slk + 788);
    const auto *slk_789 = buffer.data(slk + 789);
    const auto *slk_790 = buffer.data(slk + 790);
    const auto *slk_791 = buffer.data(slk + 791);
    const auto *slk_797 = buffer.data(slk + 797);

    const auto *sll1_674 = buffer.data(sll1 + 674);
    const auto *sll1_675 = buffer.data(sll1 + 675);
    const auto *sll1_678 = buffer.data(sll1 + 678);
    const auto *sll1_681 = buffer.data(sll1 + 681);

    const auto *smi0_553 = buffer.data(smi0 + 553);
    const auto *smi0_555 = buffer.data(smi0 + 555);
    const auto *smi0_556 = buffer.data(smi0 + 556);
    const auto *smi0_557 = buffer.data(smi0 + 557);
    const auto *smi0_558 = buffer.data(smi0 + 558);
    const auto *smi0_559 = buffer.data(smi0 + 559);
    const auto *smi0_560 = buffer.data(smi0 + 560);
    const auto *smi0_563 = buffer.data(smi0 + 563);
    const auto *smi0_565 = buffer.data(smi0 + 565);
    const auto *smi0_566 = buffer.data(smi0 + 566);
    const auto *smi0_569 = buffer.data(smi0 + 569);
    const auto *smi0_570 = buffer.data(smi0 + 570);
    const auto *smi0_572 = buffer.data(smi0 + 572);
    const auto *smi0_574 = buffer.data(smi0 + 574);
    const auto *smi0_575 = buffer.data(smi0 + 575);
    const auto *smi0_577 = buffer.data(smi0 + 577);
    const auto *smi0_578 = buffer.data(smi0 + 578);
    const auto *smi0_580 = buffer.data(smi0 + 580);
    const auto *smi0_581 = buffer.data(smi0 + 581);
    const auto *smi0_583 = buffer.data(smi0 + 583);
    const auto *smi0_584 = buffer.data(smi0 + 584);
    const auto *smi0_585 = buffer.data(smi0 + 585);
    const auto *smi0_586 = buffer.data(smi0 + 586);
    const auto *smi0_587 = buffer.data(smi0 + 587);
    const auto *smi0_588 = buffer.data(smi0 + 588);
    const auto *smi0_591 = buffer.data(smi0 + 591);
    const auto *smi0_593 = buffer.data(smi0 + 593);
    const auto *smi0_594 = buffer.data(smi0 + 594);
    const auto *smi0_597 = buffer.data(smi0 + 597);
    const auto *smi0_598 = buffer.data(smi0 + 598);
    const auto *smi0_600 = buffer.data(smi0 + 600);
    const auto *smi0_602 = buffer.data(smi0 + 602);
    const auto *smi0_603 = buffer.data(smi0 + 603);
    const auto *smi0_605 = buffer.data(smi0 + 605);
    const auto *smi0_606 = buffer.data(smi0 + 606);
    const auto *smi0_608 = buffer.data(smi0 + 608);
    const auto *smi0_609 = buffer.data(smi0 + 609);
    const auto *smi0_611 = buffer.data(smi0 + 611);
    const auto *smi0_612 = buffer.data(smi0 + 612);
    const auto *smi0_613 = buffer.data(smi0 + 613);
    const auto *smi0_614 = buffer.data(smi0 + 614);
    const auto *smi0_615 = buffer.data(smi0 + 615);
    const auto *smi0_621 = buffer.data(smi0 + 621);

    const auto *smi1_553 = buffer.data(smi1 + 553);
    const auto *smi1_555 = buffer.data(smi1 + 555);
    const auto *smi1_556 = buffer.data(smi1 + 556);
    const auto *smi1_557 = buffer.data(smi1 + 557);
    const auto *smi1_558 = buffer.data(smi1 + 558);
    const auto *smi1_559 = buffer.data(smi1 + 559);
    const auto *smi1_560 = buffer.data(smi1 + 560);
    const auto *smi1_563 = buffer.data(smi1 + 563);
    const auto *smi1_565 = buffer.data(smi1 + 565);
    const auto *smi1_566 = buffer.data(smi1 + 566);
    const auto *smi1_569 = buffer.data(smi1 + 569);
    const auto *smi1_570 = buffer.data(smi1 + 570);
    const auto *smi1_572 = buffer.data(smi1 + 572);
    const auto *smi1_574 = buffer.data(smi1 + 574);
    const auto *smi1_575 = buffer.data(smi1 + 575);
    const auto *smi1_577 = buffer.data(smi1 + 577);
    const auto *smi1_578 = buffer.data(smi1 + 578);
    const auto *smi1_580 = buffer.data(smi1 + 580);
    const auto *smi1_581 = buffer.data(smi1 + 581);
    const auto *smi1_583 = buffer.data(smi1 + 583);
    const auto *smi1_584 = buffer.data(smi1 + 584);
    const auto *smi1_585 = buffer.data(smi1 + 585);
    const auto *smi1_586 = buffer.data(smi1 + 586);
    const auto *smi1_587 = buffer.data(smi1 + 587);
    const auto *smi1_588 = buffer.data(smi1 + 588);
    const auto *smi1_591 = buffer.data(smi1 + 591);
    const auto *smi1_593 = buffer.data(smi1 + 593);
    const auto *smi1_594 = buffer.data(smi1 + 594);
    const auto *smi1_597 = buffer.data(smi1 + 597);
    const auto *smi1_598 = buffer.data(smi1 + 598);
    const auto *smi1_600 = buffer.data(smi1 + 600);
    const auto *smi1_602 = buffer.data(smi1 + 602);
    const auto *smi1_603 = buffer.data(smi1 + 603);
    const auto *smi1_605 = buffer.data(smi1 + 605);
    const auto *smi1_606 = buffer.data(smi1 + 606);
    const auto *smi1_608 = buffer.data(smi1 + 608);
    const auto *smi1_609 = buffer.data(smi1 + 609);
    const auto *smi1_611 = buffer.data(smi1 + 611);
    const auto *smi1_612 = buffer.data(smi1 + 612);
    const auto *smi1_613 = buffer.data(smi1 + 613);
    const auto *smi1_614 = buffer.data(smi1 + 614);
    const auto *smi1_615 = buffer.data(smi1 + 615);
    const auto *smi1_621 = buffer.data(smi1 + 621);

    const auto *smk_712 = buffer.data(smk + 712);
    const auto *smk_714 = buffer.data(smk + 714);
    const auto *smk_715 = buffer.data(smk + 715);
    const auto *smk_716 = buffer.data(smk + 716);
    const auto *smk_717 = buffer.data(smk + 717);
    const auto *smk_718 = buffer.data(smk + 718);
    const auto *smk_719 = buffer.data(smk + 719);
    const auto *smk_720 = buffer.data(smk + 720);
    const auto *smk_722 = buffer.data(smk + 722);
    const auto *smk_723 = buffer.data(smk + 723);
    const auto *smk_725 = buffer.data(smk + 725);
    const auto *smk_726 = buffer.data(smk + 726);
    const auto *smk_729 = buffer.data(smk + 729);
    const auto *smk_730 = buffer.data(smk + 730);
    const auto *smk_732 = buffer.data(smk + 732);
    const auto *smk_734 = buffer.data(smk + 734);
    const auto *smk_735 = buffer.data(smk + 735);
    const auto *smk_737 = buffer.data(smk + 737);
    const auto *smk_738 = buffer.data(smk + 738);
    const auto *smk_740 = buffer.data(smk + 740);
    const auto *smk_741 = buffer.data(smk + 741);
    const auto *smk_743 = buffer.data(smk + 743);
    const auto *smk_744 = buffer.data(smk + 744);
    const auto *smk_745 = buffer.data(smk + 745);
    const auto *smk_747 = buffer.data(smk + 747);
    const auto *smk_748 = buffer.data(smk + 748);
    const auto *smk_749 = buffer.data(smk + 749);
    const auto *smk_750 = buffer.data(smk + 750);
    const auto *smk_751 = buffer.data(smk + 751);
    const auto *smk_752 = buffer.data(smk + 752);
    const auto *smk_753 = buffer.data(smk + 753);
    const auto *smk_754 = buffer.data(smk + 754);
    const auto *smk_755 = buffer.data(smk + 755);
    const auto *smk_756 = buffer.data(smk + 756);
    const auto *smk_758 = buffer.data(smk + 758);
    const auto *smk_759 = buffer.data(smk + 759);
    const auto *smk_761 = buffer.data(smk + 761);
    const auto *smk_762 = buffer.data(smk + 762);
    const auto *smk_765 = buffer.data(smk + 765);
    const auto *smk_766 = buffer.data(smk + 766);
    const auto *smk_768 = buffer.data(smk + 768);
    const auto *smk_770 = buffer.data(smk + 770);
    const auto *smk_771 = buffer.data(smk + 771);
    const auto *smk_773 = buffer.data(smk + 773);
    const auto *smk_774 = buffer.data(smk + 774);
    const auto *smk_776 = buffer.data(smk + 776);
    const auto *smk_777 = buffer.data(smk + 777);
    const auto *smk_779 = buffer.data(smk + 779);
    const auto *smk_780 = buffer.data(smk + 780);
    const auto *smk_781 = buffer.data(smk + 781);
    const auto *smk_783 = buffer.data(smk + 783);
    const auto *smk_784 = buffer.data(smk + 784);
    const auto *smk_785 = buffer.data(smk + 785);
    const auto *smk_786 = buffer.data(smk + 786);
    const auto *smk_787 = buffer.data(smk + 787);
    const auto *smk_788 = buffer.data(smk + 788);
    const auto *smk_789 = buffer.data(smk + 789);
    const auto *smk_790 = buffer.data(smk + 790);
    const auto *smk_791 = buffer.data(smk + 791);
    const auto *smk_792 = buffer.data(smk + 792);
    const auto *smk_794 = buffer.data(smk + 794);
    const auto *smk_795 = buffer.data(smk + 795);
    const auto *smk_797 = buffer.data(smk + 797);

#pragma omp simd aligned(t_889, t_890, t_891, t_892, pc_x, pc_y, pc_z, slk_496, slk_532, \
                         slk_718, slk_719, smi0_553, smi1_553, smk_712, smk_718, \
                         smk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_18 * slk_718[k]
                   + f_3 * pc_x[k] * smk_718[k];

        t_890[k] = f_18 * slk_719[k]
                   + f_3 * pc_x[k] * smk_719[k];

        t_891[k] = f_15 * slk_532[k]
                   + f_1 * smi0_553[k]
                   - f_2 * smi1_553[k]
                   + f_3 * pc_y[k] * smk_712[k];

        t_892[k] = f_18 * slk_496[k]
                   + f_3 * pc_z[k] * smk_712[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, slk_534, slk_535, slk_536, smi0_555, \
                         smi0_556, smi0_557, smi1_555, smi1_556, smi1_557, smk_714, smk_715, \
                         smk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_15 * slk_534[k]
                   + f_4 * smi0_555[k]
                   - f_5 * smi1_555[k]
                   + f_3 * pc_y[k] * smk_714[k];

        t_894[k] = f_15 * slk_535[k]
                   + f_6 * smi0_556[k]
                   - f_7 * smi1_556[k]
                   + f_3 * pc_y[k] * smk_715[k];

        t_895[k] = f_15 * slk_536[k]
                   + f_8 * smi0_557[k]
                   - f_9 * smi1_557[k]
                   + f_3 * pc_y[k] * smk_716[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_y, slk_537, slk_538, slk_539, smi0_558, \
                         smi0_559, smi1_558, smi1_559, smk_717, smk_718, \
                         smk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_15 * slk_537[k]
                   + f_10 * smi0_558[k]
                   - f_11 * smi1_558[k]
                   + f_3 * pc_y[k] * smk_717[k];

        t_897[k] = f_15 * slk_538[k]
                   + f_12 * smi0_559[k]
                   - f_13 * smi1_559[k]
                   + f_3 * pc_y[k] * smk_718[k];

        t_898[k] = f_15 * slk_539[k]
                   + f_3 * pc_y[k] * smk_719[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pb_y, pc_x, pc_y, pc_z, sll0_674, \
                         slk_504, slk_720, sll1_674, smi0_560, smi1_560, \
                         smk_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pb_y[k] * sll0_674[k]
                   - f_14 * pc_y[k] * sll1_674[k];

        t_900[k] = f_18 * slk_720[k]
                   + f_1 * smi0_560[k]
                   - f_2 * smi1_560[k]
                   + f_3 * pc_x[k] * smk_720[k];

        t_901[k] = f_3 * pc_y[k] * smk_720[k];

        t_902[k] = f_19 * slk_504[k]
                   + f_3 * pc_z[k] * smk_720[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pc_x, pc_y, slk_723, slk_725, smi0_563, \
                         smi0_565, smi1_563, smi1_565, smk_722, smk_723, \
                         smk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_18 * slk_723[k]
                   + f_4 * smi0_563[k]
                   - f_5 * smi1_563[k]
                   + f_3 * pc_x[k] * smk_723[k];

        t_904[k] = f_3 * pc_y[k] * smk_722[k];

        t_905[k] = f_18 * slk_725[k]
                   + f_4 * smi0_565[k]
                   - f_5 * smi1_565[k]
                   + f_3 * pc_x[k] * smk_725[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pc_x, pc_y, pc_z, slk_507, slk_726, smi0_566, \
                         smi1_566, smk_723, smk_725, smk_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_18 * slk_726[k]
                   + f_6 * smi0_566[k]
                   - f_7 * smi1_566[k]
                   + f_3 * pc_x[k] * smk_726[k];

        t_907[k] = f_19 * slk_507[k]
                   + f_3 * pc_z[k] * smk_723[k];

        t_908[k] = f_3 * pc_y[k] * smk_725[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, pc_z, slk_510, slk_729, slk_730, smi0_569, \
                         smi0_570, smi1_569, smi1_570, smk_726, smk_729, \
                         smk_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_18 * slk_729[k]
                   + f_6 * smi0_569[k]
                   - f_7 * smi1_569[k]
                   + f_3 * pc_x[k] * smk_729[k];

        t_910[k] = f_18 * slk_730[k]
                   + f_8 * smi0_570[k]
                   - f_9 * smi1_570[k]
                   + f_3 * pc_x[k] * smk_730[k];

        t_911[k] = f_19 * slk_510[k]
                   + f_3 * pc_z[k] * smk_726[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_x, pc_y, slk_732, slk_734, smi0_572, \
                         smi0_574, smi1_572, smi1_574, smk_729, smk_732, \
                         smk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_18 * slk_732[k]
                   + f_8 * smi0_572[k]
                   - f_9 * smi1_572[k]
                   + f_3 * pc_x[k] * smk_732[k];

        t_913[k] = f_3 * pc_y[k] * smk_729[k];

        t_914[k] = f_18 * slk_734[k]
                   + f_8 * smi0_574[k]
                   - f_9 * smi1_574[k]
                   + f_3 * pc_x[k] * smk_734[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_x, pc_z, slk_514, slk_735, slk_737, smi0_575, \
                         smi0_577, smi1_575, smi1_577, smk_730, smk_735, \
                         smk_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_18 * slk_735[k]
                   + f_10 * smi0_575[k]
                   - f_11 * smi1_575[k]
                   + f_3 * pc_x[k] * smk_735[k];

        t_916[k] = f_19 * slk_514[k]
                   + f_3 * pc_z[k] * smk_730[k];

        t_917[k] = f_18 * slk_737[k]
                   + f_10 * smi0_577[k]
                   - f_11 * smi1_577[k]
                   + f_3 * pc_x[k] * smk_737[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_x, pc_y, slk_738, slk_740, smi0_578, \
                         smi0_580, smi1_578, smi1_580, smk_734, smk_738, \
                         smk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_18 * slk_738[k]
                   + f_10 * smi0_578[k]
                   - f_11 * smi1_578[k]
                   + f_3 * pc_x[k] * smk_738[k];

        t_919[k] = f_3 * pc_y[k] * smk_734[k];

        t_920[k] = f_18 * slk_740[k]
                   + f_10 * smi0_580[k]
                   - f_11 * smi1_580[k]
                   + f_3 * pc_x[k] * smk_740[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pc_x, pc_z, slk_519, slk_741, slk_743, smi0_581, \
                         smi0_583, smi1_581, smi1_583, smk_735, smk_741, \
                         smk_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_18 * slk_741[k]
                   + f_12 * smi0_581[k]
                   - f_13 * smi1_581[k]
                   + f_3 * pc_x[k] * smk_741[k];

        t_922[k] = f_19 * slk_519[k]
                   + f_3 * pc_z[k] * smk_735[k];

        t_923[k] = f_18 * slk_743[k]
                   + f_12 * smi0_583[k]
                   - f_13 * smi1_583[k]
                   + f_3 * pc_x[k] * smk_743[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_x, pc_y, slk_744, slk_745, smi0_584, \
                         smi0_585, smi1_584, smi1_585, smk_740, smk_744, \
                         smk_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_18 * slk_744[k]
                   + f_12 * smi0_584[k]
                   - f_13 * smi1_584[k]
                   + f_3 * pc_x[k] * smk_744[k];

        t_925[k] = f_18 * slk_745[k]
                   + f_12 * smi0_585[k]
                   - f_13 * smi1_585[k]
                   + f_3 * pc_x[k] * smk_745[k];

        t_926[k] = f_3 * pc_y[k] * smk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pc_x, slk_747, slk_748, slk_749, slk_750, \
                         smi0_587, smi1_587, smk_747, smk_748, smk_749, \
                         smk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_18 * slk_747[k]
                   + f_12 * smi0_587[k]
                   - f_13 * smi1_587[k]
                   + f_3 * pc_x[k] * smk_747[k];

        t_928[k] = f_18 * slk_748[k]
                   + f_3 * pc_x[k] * smk_748[k];

        t_929[k] = f_18 * slk_749[k]
                   + f_3 * pc_x[k] * smk_749[k];

        t_930[k] = f_18 * slk_750[k]
                   + f_3 * pc_x[k] * smk_750[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pc_x, slk_751, slk_752, slk_753, \
                         slk_754, slk_755, smk_751, smk_752, smk_753, smk_754, \
                         smk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_18 * slk_751[k]
                   + f_3 * pc_x[k] * smk_751[k];

        t_932[k] = f_18 * slk_752[k]
                   + f_3 * pc_x[k] * smk_752[k];

        t_933[k] = f_18 * slk_753[k]
                   + f_3 * pc_x[k] * smk_753[k];

        t_934[k] = f_18 * slk_754[k]
                   + f_3 * pc_x[k] * smk_754[k];

        t_935[k] = f_18 * slk_755[k]
                   + f_3 * pc_x[k] * smk_755[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, pc_y, pc_z, slk_532, smi0_581, smi0_583, \
                         smi0_584, smi1_581, smi1_583, smi1_584, smk_748, smk_750, \
                         smk_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_1 * smi0_581[k]
                   - f_2 * smi1_581[k]
                   + f_3 * pc_y[k] * smk_748[k];

        t_937[k] = f_19 * slk_532[k]
                   + f_3 * pc_z[k] * smk_748[k];

        t_938[k] = f_4 * smi0_583[k]
                   - f_5 * smi1_583[k]
                   + f_3 * pc_y[k] * smk_750[k];

        t_939[k] = f_6 * smi0_584[k]
                   - f_7 * smi1_584[k]
                   + f_3 * pc_y[k] * smk_751[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, pc_y, smi0_585, smi0_586, smi0_587, \
                         smi1_585, smi1_586, smi1_587, smk_752, smk_753, smk_754, \
                         smk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_8 * smi0_585[k]
                   - f_9 * smi1_585[k]
                   + f_3 * pc_y[k] * smk_752[k];

        t_941[k] = f_10 * smi0_586[k]
                   - f_11 * smi1_586[k]
                   + f_3 * pc_y[k] * smk_753[k];

        t_942[k] = f_12 * smi0_587[k]
                   - f_13 * smi1_587[k]
                   + f_3 * pc_y[k] * smk_754[k];

        t_943[k] = f_3 * pc_y[k] * smk_755[k];
    }

#pragma omp simd aligned(t_944, t_945, t_946, t_947, pc_x, pc_y, pc_z, slk_539, slk_540, \
                         slk_756, smi0_587, smi0_588, smi1_587, smi1_588, smk_755, \
                         smk_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_19 * slk_539[k]
                   + f_1 * smi0_587[k]
                   - f_2 * smi1_587[k]
                   + f_3 * pc_z[k] * smk_755[k];

        t_945[k] = f_17 * slk_756[k]
                   + f_1 * smi0_588[k]
                   - f_2 * smi1_588[k]
                   + f_3 * pc_x[k] * smk_756[k];

        t_946[k] = f_20 * slk_540[k]
                   + f_3 * pc_y[k] * smk_756[k];

        t_947[k] = f_3 * pc_z[k] * smk_756[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, pc_x, pc_y, slk_542, slk_759, slk_761, smi0_591, \
                         smi0_593, smi1_591, smi1_593, smk_758, smk_759, \
                         smk_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = f_17 * slk_759[k]
                   + f_4 * smi0_591[k]
                   - f_5 * smi1_591[k]
                   + f_3 * pc_x[k] * smk_759[k];

        t_949[k] = f_20 * slk_542[k]
                   + f_3 * pc_y[k] * smk_758[k];

        t_950[k] = f_17 * slk_761[k]
                   + f_4 * smi0_593[k]
                   - f_5 * smi1_593[k]
                   + f_3 * pc_x[k] * smk_761[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, pc_x, pc_y, pc_z, slk_545, slk_762, smi0_594, \
                         smi1_594, smk_759, smk_761, smk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_17 * slk_762[k]
                   + f_6 * smi0_594[k]
                   - f_7 * smi1_594[k]
                   + f_3 * pc_x[k] * smk_762[k];

        t_952[k] = f_3 * pc_z[k] * smk_759[k];

        t_953[k] = f_20 * slk_545[k]
                   + f_3 * pc_y[k] * smk_761[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, pc_x, pc_z, slk_765, slk_766, smi0_597, \
                         smi0_598, smi1_597, smi1_598, smk_762, smk_765, \
                         smk_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_17 * slk_765[k]
                   + f_6 * smi0_597[k]
                   - f_7 * smi1_597[k]
                   + f_3 * pc_x[k] * smk_765[k];

        t_955[k] = f_17 * slk_766[k]
                   + f_8 * smi0_598[k]
                   - f_9 * smi1_598[k]
                   + f_3 * pc_x[k] * smk_766[k];

        t_956[k] = f_3 * pc_z[k] * smk_762[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, pc_x, pc_y, slk_549, slk_768, slk_770, smi0_600, \
                         smi0_602, smi1_600, smi1_602, smk_765, smk_768, \
                         smk_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_17 * slk_768[k]
                   + f_8 * smi0_600[k]
                   - f_9 * smi1_600[k]
                   + f_3 * pc_x[k] * smk_768[k];

        t_958[k] = f_20 * slk_549[k]
                   + f_3 * pc_y[k] * smk_765[k];

        t_959[k] = f_17 * slk_770[k]
                   + f_8 * smi0_602[k]
                   - f_9 * smi1_602[k]
                   + f_3 * pc_x[k] * smk_770[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, pc_x, pc_z, slk_771, slk_773, smi0_603, \
                         smi0_605, smi1_603, smi1_605, smk_766, smk_771, \
                         smk_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_17 * slk_771[k]
                   + f_10 * smi0_603[k]
                   - f_11 * smi1_603[k]
                   + f_3 * pc_x[k] * smk_771[k];

        t_961[k] = f_3 * pc_z[k] * smk_766[k];

        t_962[k] = f_17 * slk_773[k]
                   + f_10 * smi0_605[k]
                   - f_11 * smi1_605[k]
                   + f_3 * pc_x[k] * smk_773[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, pc_x, pc_y, slk_554, slk_774, slk_776, smi0_606, \
                         smi0_608, smi1_606, smi1_608, smk_770, smk_774, \
                         smk_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = f_17 * slk_774[k]
                   + f_10 * smi0_606[k]
                   - f_11 * smi1_606[k]
                   + f_3 * pc_x[k] * smk_774[k];

        t_964[k] = f_20 * slk_554[k]
                   + f_3 * pc_y[k] * smk_770[k];

        t_965[k] = f_17 * slk_776[k]
                   + f_10 * smi0_608[k]
                   - f_11 * smi1_608[k]
                   + f_3 * pc_x[k] * smk_776[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_x, pc_z, slk_777, slk_779, smi0_609, \
                         smi0_611, smi1_609, smi1_611, smk_771, smk_777, \
                         smk_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_17 * slk_777[k]
                   + f_12 * smi0_609[k]
                   - f_13 * smi1_609[k]
                   + f_3 * pc_x[k] * smk_777[k];

        t_967[k] = f_3 * pc_z[k] * smk_771[k];

        t_968[k] = f_17 * slk_779[k]
                   + f_12 * smi0_611[k]
                   - f_13 * smi1_611[k]
                   + f_3 * pc_x[k] * smk_779[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pc_x, pc_y, slk_560, slk_780, slk_781, smi0_612, \
                         smi0_613, smi1_612, smi1_613, smk_776, smk_780, \
                         smk_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_17 * slk_780[k]
                   + f_12 * smi0_612[k]
                   - f_13 * smi1_612[k]
                   + f_3 * pc_x[k] * smk_780[k];

        t_970[k] = f_17 * slk_781[k]
                   + f_12 * smi0_613[k]
                   - f_13 * smi1_613[k]
                   + f_3 * pc_x[k] * smk_781[k];

        t_971[k] = f_20 * slk_560[k]
                   + f_3 * pc_y[k] * smk_776[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, pc_x, slk_783, slk_784, slk_785, slk_786, \
                         smi0_615, smi1_615, smk_783, smk_784, smk_785, \
                         smk_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_17 * slk_783[k]
                   + f_12 * smi0_615[k]
                   - f_13 * smi1_615[k]
                   + f_3 * pc_x[k] * smk_783[k];

        t_973[k] = f_17 * slk_784[k]
                   + f_3 * pc_x[k] * smk_784[k];

        t_974[k] = f_17 * slk_785[k]
                   + f_3 * pc_x[k] * smk_785[k];

        t_975[k] = f_17 * slk_786[k]
                   + f_3 * pc_x[k] * smk_786[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, t_979, t_980, pc_x, slk_787, slk_788, slk_789, \
                         slk_790, slk_791, smk_787, smk_788, smk_789, smk_790, \
                         smk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_17 * slk_787[k]
                   + f_3 * pc_x[k] * smk_787[k];

        t_977[k] = f_17 * slk_788[k]
                   + f_3 * pc_x[k] * smk_788[k];

        t_978[k] = f_17 * slk_789[k]
                   + f_3 * pc_x[k] * smk_789[k];

        t_979[k] = f_17 * slk_790[k]
                   + f_3 * pc_x[k] * smk_790[k];

        t_980[k] = f_17 * slk_791[k]
                   + f_3 * pc_x[k] * smk_791[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_y, pc_z, slk_568, slk_570, smi0_609, \
                         smi0_611, smi1_609, smi1_611, smk_784, \
                         smk_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_20 * slk_568[k]
                   + f_1 * smi0_609[k]
                   - f_2 * smi1_609[k]
                   + f_3 * pc_y[k] * smk_784[k];

        t_982[k] = f_3 * pc_z[k] * smk_784[k];

        t_983[k] = f_20 * slk_570[k]
                   + f_4 * smi0_611[k]
                   - f_5 * smi1_611[k]
                   + f_3 * pc_y[k] * smk_786[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pc_y, slk_571, slk_572, slk_573, smi0_612, \
                         smi0_613, smi0_614, smi1_612, smi1_613, smi1_614, smk_787, smk_788, \
                         smk_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_20 * slk_571[k]
                   + f_6 * smi0_612[k]
                   - f_7 * smi1_612[k]
                   + f_3 * pc_y[k] * smk_787[k];

        t_985[k] = f_20 * slk_572[k]
                   + f_8 * smi0_613[k]
                   - f_9 * smi1_613[k]
                   + f_3 * pc_y[k] * smk_788[k];

        t_986[k] = f_20 * slk_573[k]
                   + f_10 * smi0_614[k]
                   - f_11 * smi1_614[k]
                   + f_3 * pc_y[k] * smk_789[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, t_990, pb_z, pc_y, pc_z, sll0_675, slk_574, \
                         slk_575, sll1_675, smi0_615, smi1_615, smk_790, \
                         smk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_20 * slk_574[k]
                   + f_12 * smi0_615[k]
                   - f_13 * smi1_615[k]
                   + f_3 * pc_y[k] * smk_790[k];

        t_988[k] = f_20 * slk_575[k]
                   + f_3 * pc_y[k] * smk_791[k];

        t_989[k] = f_1 * smi0_615[k]
                   - f_2 * smi1_615[k]
                   + f_3 * pc_z[k] * smk_791[k];

        t_990[k] = pb_z[k] * sll0_675[k]
                   - f_14 * pc_z[k] * sll1_675[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pb_z, pc_y, pc_z, sll0_678, slk_540, \
                         slk_576, slk_578, sll1_678, smk_792, smk_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_19 * slk_576[k]
                   + f_3 * pc_y[k] * smk_792[k];

        t_992[k] = f_15 * slk_540[k]
                   + f_3 * pc_z[k] * smk_792[k];

        t_993[k] = pb_z[k] * sll0_678[k]
                   - f_14 * pc_z[k] * sll1_678[k];

        t_994[k] = f_19 * slk_578[k]
                   + f_3 * pc_y[k] * smk_794[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, pb_z, pc_x, pc_z, sll0_681, slk_543, slk_797, \
                         sll1_681, smi0_621, smi1_621, smk_795, \
                         smk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_17 * slk_797[k]
                   + f_4 * smi0_621[k]
                   - f_5 * smi1_621[k]
                   + f_3 * pc_x[k] * smk_797[k];

        t_996[k] = pb_z[k] * sll0_681[k]
                   - f_14 * pc_z[k] * sll1_681[k];

        t_997[k] = f_15 * slk_543[k]
                   + f_3 * pc_z[k] * smk_795[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sll0,
                                                          const size_t slk, const size_t sll1,
                                                          const size_t smi0, const size_t smi1,
                                                          const size_t smk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_685 = buffer.data(sll0 + 685);
    const auto *sll0_687 = buffer.data(sll0 + 687);
    const auto *sll0_690 = buffer.data(sll0 + 690);
    const auto *sll0_692 = buffer.data(sll0 + 692);
    const auto *sll0_693 = buffer.data(sll0 + 693);
    const auto *sll0_696 = buffer.data(sll0 + 696);
    const auto *sll0_698 = buffer.data(sll0 + 698);
    const auto *sll0_699 = buffer.data(sll0 + 699);
    const auto *sll0_700 = buffer.data(sll0 + 700);
    const auto *sll0_711 = buffer.data(sll0 + 711);

    const auto *slk_546 = buffer.data(slk + 546);
    const auto *slk_547 = buffer.data(slk + 547);
    const auto *slk_550 = buffer.data(slk + 550);
    const auto *slk_551 = buffer.data(slk + 551);
    const auto *slk_552 = buffer.data(slk + 552);
    const auto *slk_555 = buffer.data(slk + 555);
    const auto *slk_556 = buffer.data(slk + 556);
    const auto *slk_557 = buffer.data(slk + 557);
    const auto *slk_558 = buffer.data(slk + 558);
    const auto *slk_568 = buffer.data(slk + 568);
    const auto *slk_575 = buffer.data(slk + 575);
    const auto *slk_576 = buffer.data(slk + 576);
    const auto *slk_579 = buffer.data(slk + 579);
    const auto *slk_581 = buffer.data(slk + 581);
    const auto *slk_582 = buffer.data(slk + 582);
    const auto *slk_585 = buffer.data(slk + 585);
    const auto *slk_586 = buffer.data(slk + 586);
    const auto *slk_590 = buffer.data(slk + 590);
    const auto *slk_591 = buffer.data(slk + 591);
    const auto *slk_596 = buffer.data(slk + 596);
    const auto *slk_604 = buffer.data(slk + 604);
    const auto *slk_606 = buffer.data(slk + 606);
    const auto *slk_607 = buffer.data(slk + 607);
    const auto *slk_608 = buffer.data(slk + 608);
    const auto *slk_609 = buffer.data(slk + 609);
    const auto *slk_610 = buffer.data(slk + 610);
    const auto *slk_611 = buffer.data(slk + 611);
    const auto *slk_612 = buffer.data(slk + 612);
    const auto *slk_614 = buffer.data(slk + 614);
    const auto *slk_615 = buffer.data(slk + 615);
    const auto *slk_617 = buffer.data(slk + 617);
    const auto *slk_618 = buffer.data(slk + 618);
    const auto *slk_621 = buffer.data(slk + 621);
    const auto *slk_622 = buffer.data(slk + 622);
    const auto *slk_626 = buffer.data(slk + 626);
    const auto *slk_627 = buffer.data(slk + 627);
    const auto *slk_632 = buffer.data(slk + 632);
    const auto *slk_640 = buffer.data(slk + 640);
    const auto *slk_642 = buffer.data(slk + 642);
    const auto *slk_643 = buffer.data(slk + 643);
    const auto *slk_644 = buffer.data(slk + 644);
    const auto *slk_645 = buffer.data(slk + 645);
    const auto *slk_646 = buffer.data(slk + 646);
    const auto *slk_647 = buffer.data(slk + 647);
    const auto *slk_648 = buffer.data(slk + 648);
    const auto *slk_650 = buffer.data(slk + 650);
    const auto *slk_653 = buffer.data(slk + 653);
    const auto *slk_657 = buffer.data(slk + 657);
    const auto *slk_662 = buffer.data(slk + 662);
    const auto *slk_801 = buffer.data(slk + 801);
    const auto *slk_806 = buffer.data(slk + 806);
    const auto *slk_812 = buffer.data(slk + 812);
    const auto *slk_819 = buffer.data(slk + 819);
    const auto *slk_820 = buffer.data(slk + 820);
    const auto *slk_821 = buffer.data(slk + 821);
    const auto *slk_822 = buffer.data(slk + 822);
    const auto *slk_823 = buffer.data(slk + 823);
    const auto *slk_824 = buffer.data(slk + 824);
    const auto *slk_825 = buffer.data(slk + 825);
    const auto *slk_826 = buffer.data(slk + 826);
    const auto *slk_827 = buffer.data(slk + 827);
    const auto *slk_828 = buffer.data(slk + 828);
    const auto *slk_831 = buffer.data(slk + 831);
    const auto *slk_833 = buffer.data(slk + 833);
    const auto *slk_834 = buffer.data(slk + 834);
    const auto *slk_837 = buffer.data(slk + 837);
    const auto *slk_838 = buffer.data(slk + 838);
    const auto *slk_840 = buffer.data(slk + 840);
    const auto *slk_842 = buffer.data(slk + 842);
    const auto *slk_843 = buffer.data(slk + 843);
    const auto *slk_845 = buffer.data(slk + 845);
    const auto *slk_846 = buffer.data(slk + 846);
    const auto *slk_848 = buffer.data(slk + 848);
    const auto *slk_849 = buffer.data(slk + 849);
    const auto *slk_851 = buffer.data(slk + 851);
    const auto *slk_852 = buffer.data(slk + 852);
    const auto *slk_853 = buffer.data(slk + 853);
    const auto *slk_855 = buffer.data(slk + 855);
    const auto *slk_856 = buffer.data(slk + 856);
    const auto *slk_857 = buffer.data(slk + 857);
    const auto *slk_858 = buffer.data(slk + 858);
    const auto *slk_859 = buffer.data(slk + 859);
    const auto *slk_860 = buffer.data(slk + 860);
    const auto *slk_861 = buffer.data(slk + 861);
    const auto *slk_862 = buffer.data(slk + 862);
    const auto *slk_863 = buffer.data(slk + 863);
    const auto *slk_864 = buffer.data(slk + 864);
    const auto *slk_867 = buffer.data(slk + 867);
    const auto *slk_869 = buffer.data(slk + 869);
    const auto *slk_870 = buffer.data(slk + 870);
    const auto *slk_873 = buffer.data(slk + 873);
    const auto *slk_874 = buffer.data(slk + 874);
    const auto *slk_876 = buffer.data(slk + 876);
    const auto *slk_878 = buffer.data(slk + 878);
    const auto *slk_879 = buffer.data(slk + 879);
    const auto *slk_881 = buffer.data(slk + 881);
    const auto *slk_882 = buffer.data(slk + 882);
    const auto *slk_884 = buffer.data(slk + 884);
    const auto *slk_885 = buffer.data(slk + 885);

    const auto *sll1_685 = buffer.data(sll1 + 685);
    const auto *sll1_687 = buffer.data(sll1 + 687);
    const auto *sll1_690 = buffer.data(sll1 + 690);
    const auto *sll1_692 = buffer.data(sll1 + 692);
    const auto *sll1_693 = buffer.data(sll1 + 693);
    const auto *sll1_696 = buffer.data(sll1 + 696);
    const auto *sll1_698 = buffer.data(sll1 + 698);
    const auto *sll1_699 = buffer.data(sll1 + 699);
    const auto *sll1_700 = buffer.data(sll1 + 700);
    const auto *sll1_711 = buffer.data(sll1 + 711);

    const auto *smi0_625 = buffer.data(smi0 + 625);
    const auto *smi0_630 = buffer.data(smi0 + 630);
    const auto *smi0_636 = buffer.data(smi0 + 636);
    const auto *smi0_639 = buffer.data(smi0 + 639);
    const auto *smi0_640 = buffer.data(smi0 + 640);
    const auto *smi0_641 = buffer.data(smi0 + 641);
    const auto *smi0_642 = buffer.data(smi0 + 642);
    const auto *smi0_643 = buffer.data(smi0 + 643);
    const auto *smi0_644 = buffer.data(smi0 + 644);
    const auto *smi0_647 = buffer.data(smi0 + 647);
    const auto *smi0_649 = buffer.data(smi0 + 649);
    const auto *smi0_650 = buffer.data(smi0 + 650);
    const auto *smi0_653 = buffer.data(smi0 + 653);
    const auto *smi0_654 = buffer.data(smi0 + 654);
    const auto *smi0_656 = buffer.data(smi0 + 656);
    const auto *smi0_658 = buffer.data(smi0 + 658);
    const auto *smi0_659 = buffer.data(smi0 + 659);
    const auto *smi0_661 = buffer.data(smi0 + 661);
    const auto *smi0_662 = buffer.data(smi0 + 662);
    const auto *smi0_664 = buffer.data(smi0 + 664);
    const auto *smi0_665 = buffer.data(smi0 + 665);
    const auto *smi0_667 = buffer.data(smi0 + 667);
    const auto *smi0_668 = buffer.data(smi0 + 668);
    const auto *smi0_669 = buffer.data(smi0 + 669);
    const auto *smi0_670 = buffer.data(smi0 + 670);
    const auto *smi0_671 = buffer.data(smi0 + 671);
    const auto *smi0_672 = buffer.data(smi0 + 672);
    const auto *smi0_675 = buffer.data(smi0 + 675);
    const auto *smi0_677 = buffer.data(smi0 + 677);
    const auto *smi0_678 = buffer.data(smi0 + 678);
    const auto *smi0_681 = buffer.data(smi0 + 681);
    const auto *smi0_682 = buffer.data(smi0 + 682);
    const auto *smi0_684 = buffer.data(smi0 + 684);
    const auto *smi0_686 = buffer.data(smi0 + 686);
    const auto *smi0_687 = buffer.data(smi0 + 687);
    const auto *smi0_689 = buffer.data(smi0 + 689);
    const auto *smi0_690 = buffer.data(smi0 + 690);
    const auto *smi0_692 = buffer.data(smi0 + 692);
    const auto *smi0_693 = buffer.data(smi0 + 693);

    const auto *smi1_625 = buffer.data(smi1 + 625);
    const auto *smi1_630 = buffer.data(smi1 + 630);
    const auto *smi1_636 = buffer.data(smi1 + 636);
    const auto *smi1_639 = buffer.data(smi1 + 639);
    const auto *smi1_640 = buffer.data(smi1 + 640);
    const auto *smi1_641 = buffer.data(smi1 + 641);
    const auto *smi1_642 = buffer.data(smi1 + 642);
    const auto *smi1_643 = buffer.data(smi1 + 643);
    const auto *smi1_644 = buffer.data(smi1 + 644);
    const auto *smi1_647 = buffer.data(smi1 + 647);
    const auto *smi1_649 = buffer.data(smi1 + 649);
    const auto *smi1_650 = buffer.data(smi1 + 650);
    const auto *smi1_653 = buffer.data(smi1 + 653);
    const auto *smi1_654 = buffer.data(smi1 + 654);
    const auto *smi1_656 = buffer.data(smi1 + 656);
    const auto *smi1_658 = buffer.data(smi1 + 658);
    const auto *smi1_659 = buffer.data(smi1 + 659);
    const auto *smi1_661 = buffer.data(smi1 + 661);
    const auto *smi1_662 = buffer.data(smi1 + 662);
    const auto *smi1_664 = buffer.data(smi1 + 664);
    const auto *smi1_665 = buffer.data(smi1 + 665);
    const auto *smi1_667 = buffer.data(smi1 + 667);
    const auto *smi1_668 = buffer.data(smi1 + 668);
    const auto *smi1_669 = buffer.data(smi1 + 669);
    const auto *smi1_670 = buffer.data(smi1 + 670);
    const auto *smi1_671 = buffer.data(smi1 + 671);
    const auto *smi1_672 = buffer.data(smi1 + 672);
    const auto *smi1_675 = buffer.data(smi1 + 675);
    const auto *smi1_677 = buffer.data(smi1 + 677);
    const auto *smi1_678 = buffer.data(smi1 + 678);
    const auto *smi1_681 = buffer.data(smi1 + 681);
    const auto *smi1_682 = buffer.data(smi1 + 682);
    const auto *smi1_684 = buffer.data(smi1 + 684);
    const auto *smi1_686 = buffer.data(smi1 + 686);
    const auto *smi1_687 = buffer.data(smi1 + 687);
    const auto *smi1_689 = buffer.data(smi1 + 689);
    const auto *smi1_690 = buffer.data(smi1 + 690);
    const auto *smi1_692 = buffer.data(smi1 + 692);
    const auto *smi1_693 = buffer.data(smi1 + 693);

    const auto *smk_797 = buffer.data(smk + 797);
    const auto *smk_798 = buffer.data(smk + 798);
    const auto *smk_801 = buffer.data(smk + 801);
    const auto *smk_802 = buffer.data(smk + 802);
    const auto *smk_806 = buffer.data(smk + 806);
    const auto *smk_807 = buffer.data(smk + 807);
    const auto *smk_812 = buffer.data(smk + 812);
    const auto *smk_819 = buffer.data(smk + 819);
    const auto *smk_820 = buffer.data(smk + 820);
    const auto *smk_821 = buffer.data(smk + 821);
    const auto *smk_822 = buffer.data(smk + 822);
    const auto *smk_823 = buffer.data(smk + 823);
    const auto *smk_824 = buffer.data(smk + 824);
    const auto *smk_825 = buffer.data(smk + 825);
    const auto *smk_826 = buffer.data(smk + 826);
    const auto *smk_827 = buffer.data(smk + 827);
    const auto *smk_828 = buffer.data(smk + 828);
    const auto *smk_830 = buffer.data(smk + 830);
    const auto *smk_831 = buffer.data(smk + 831);
    const auto *smk_833 = buffer.data(smk + 833);
    const auto *smk_834 = buffer.data(smk + 834);
    const auto *smk_837 = buffer.data(smk + 837);
    const auto *smk_838 = buffer.data(smk + 838);
    const auto *smk_840 = buffer.data(smk + 840);
    const auto *smk_842 = buffer.data(smk + 842);
    const auto *smk_843 = buffer.data(smk + 843);
    const auto *smk_845 = buffer.data(smk + 845);
    const auto *smk_846 = buffer.data(smk + 846);
    const auto *smk_848 = buffer.data(smk + 848);
    const auto *smk_849 = buffer.data(smk + 849);
    const auto *smk_851 = buffer.data(smk + 851);
    const auto *smk_852 = buffer.data(smk + 852);
    const auto *smk_853 = buffer.data(smk + 853);
    const auto *smk_855 = buffer.data(smk + 855);
    const auto *smk_856 = buffer.data(smk + 856);
    const auto *smk_857 = buffer.data(smk + 857);
    const auto *smk_858 = buffer.data(smk + 858);
    const auto *smk_859 = buffer.data(smk + 859);
    const auto *smk_860 = buffer.data(smk + 860);
    const auto *smk_861 = buffer.data(smk + 861);
    const auto *smk_862 = buffer.data(smk + 862);
    const auto *smk_863 = buffer.data(smk + 863);
    const auto *smk_864 = buffer.data(smk + 864);
    const auto *smk_866 = buffer.data(smk + 866);
    const auto *smk_867 = buffer.data(smk + 867);
    const auto *smk_869 = buffer.data(smk + 869);
    const auto *smk_870 = buffer.data(smk + 870);
    const auto *smk_873 = buffer.data(smk + 873);
    const auto *smk_874 = buffer.data(smk + 874);
    const auto *smk_876 = buffer.data(smk + 876);
    const auto *smk_878 = buffer.data(smk + 878);
    const auto *smk_879 = buffer.data(smk + 879);
    const auto *smk_881 = buffer.data(smk + 881);
    const auto *smk_882 = buffer.data(smk + 882);
    const auto *smk_884 = buffer.data(smk + 884);
    const auto *smk_885 = buffer.data(smk + 885);

#pragma omp simd aligned(t_998, t_999, t_1000, pb_z, pc_x, pc_y, pc_z, sll0_685, slk_581, \
                         slk_801, sll1_685, smi0_625, smi1_625, smk_797, \
                         smk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_19 * slk_581[k]
                   + f_3 * pc_y[k] * smk_797[k];

        t_999[k] = f_17 * slk_801[k]
                   + f_6 * smi0_625[k]
                   - f_7 * smi1_625[k]
                   + f_3 * pc_x[k] * smk_801[k];

        t_1000[k] = pb_z[k] * sll0_685[k]
                    - f_14 * pc_z[k] * sll1_685[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, pb_z, pc_y, pc_z, sll0_687, slk_546, slk_547, \
                         slk_585, sll1_687, smk_798, smk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_15 * slk_546[k]
                    + f_3 * pc_z[k] * smk_798[k];

        t_1002[k] = pb_z[k] * sll0_687[k]
                    + f_16 * slk_547[k]
                    - f_14 * pc_z[k] * sll1_687[k];

        t_1003[k] = f_19 * slk_585[k]
                    + f_3 * pc_y[k] * smk_801[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, pb_z, pc_x, pc_z, sll0_690, slk_550, slk_806, \
                         sll1_690, smi0_630, smi1_630, smk_802, \
                         smk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_17 * slk_806[k]
                    + f_8 * smi0_630[k]
                    - f_9 * smi1_630[k]
                    + f_3 * pc_x[k] * smk_806[k];

        t_1005[k] = pb_z[k] * sll0_690[k]
                    - f_14 * pc_z[k] * sll1_690[k];

        t_1006[k] = f_15 * slk_550[k]
                    + f_3 * pc_z[k] * smk_802[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, pb_z, pc_y, pc_z, sll0_692, sll0_693, \
                         slk_551, slk_552, slk_590, sll1_692, sll1_693, \
                         smk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = pb_z[k] * sll0_692[k]
                    + f_16 * slk_551[k]
                    - f_14 * pc_z[k] * sll1_692[k];

        t_1008[k] = pb_z[k] * sll0_693[k]
                    + f_17 * slk_552[k]
                    - f_14 * pc_z[k] * sll1_693[k];

        t_1009[k] = f_19 * slk_590[k]
                    + f_3 * pc_y[k] * smk_806[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, pb_z, pc_x, pc_z, sll0_696, slk_555, slk_812, \
                         sll1_696, smi0_636, smi1_636, smk_807, \
                         smk_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_17 * slk_812[k]
                    + f_10 * smi0_636[k]
                    - f_11 * smi1_636[k]
                    + f_3 * pc_x[k] * smk_812[k];

        t_1011[k] = pb_z[k] * sll0_696[k]
                    - f_14 * pc_z[k] * sll1_696[k];

        t_1012[k] = f_15 * slk_555[k]
                    + f_3 * pc_z[k] * smk_807[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pb_z, pc_z, sll0_698, sll0_699, sll0_700, \
                         slk_556, slk_557, slk_558, sll1_698, sll1_699, \
                         sll1_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = pb_z[k] * sll0_698[k]
                    + f_16 * slk_556[k]
                    - f_14 * pc_z[k] * sll1_698[k];

        t_1014[k] = pb_z[k] * sll0_699[k]
                    + f_17 * slk_557[k]
                    - f_14 * pc_z[k] * sll1_699[k];

        t_1015[k] = pb_z[k] * sll0_700[k]
                    + f_18 * slk_558[k]
                    - f_14 * pc_z[k] * sll1_700[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, slk_596, slk_819, \
                         slk_820, slk_821, smi0_643, smi1_643, smk_812, smk_819, smk_820, \
                         smk_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_19 * slk_596[k]
                    + f_3 * pc_y[k] * smk_812[k];

        t_1017[k] = f_17 * slk_819[k]
                    + f_12 * smi0_643[k]
                    - f_13 * smi1_643[k]
                    + f_3 * pc_x[k] * smk_819[k];

        t_1018[k] = f_17 * slk_820[k]
                    + f_3 * pc_x[k] * smk_820[k];

        t_1019[k] = f_17 * slk_821[k]
                    + f_3 * pc_x[k] * smk_821[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, t_1024, pc_x, slk_822, slk_823, \
                         slk_824, slk_825, slk_826, smk_822, smk_823, smk_824, smk_825, \
                         smk_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_17 * slk_822[k]
                    + f_3 * pc_x[k] * smk_822[k];

        t_1021[k] = f_17 * slk_823[k]
                    + f_3 * pc_x[k] * smk_823[k];

        t_1022[k] = f_17 * slk_824[k]
                    + f_3 * pc_x[k] * smk_824[k];

        t_1023[k] = f_17 * slk_825[k]
                    + f_3 * pc_x[k] * smk_825[k];

        t_1024[k] = f_17 * slk_826[k]
                    + f_3 * pc_x[k] * smk_826[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, pb_z, pc_x, pc_z, sll0_711, slk_568, slk_827, \
                         sll1_711, smk_820, smk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_17 * slk_827[k]
                    + f_3 * pc_x[k] * smk_827[k];

        t_1026[k] = pb_z[k] * sll0_711[k]
                    - f_14 * pc_z[k] * sll1_711[k];

        t_1027[k] = f_15 * slk_568[k]
                    + f_3 * pc_z[k] * smk_820[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, pc_y, slk_606, slk_607, slk_608, smi0_639, \
                         smi0_640, smi0_641, smi1_639, smi1_640, smi1_641, smk_822, smk_823, \
                         smk_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_19 * slk_606[k]
                    + f_4 * smi0_639[k]
                    - f_5 * smi1_639[k]
                    + f_3 * pc_y[k] * smk_822[k];

        t_1029[k] = f_19 * slk_607[k]
                    + f_6 * smi0_640[k]
                    - f_7 * smi1_640[k]
                    + f_3 * pc_y[k] * smk_823[k];

        t_1030[k] = f_19 * slk_608[k]
                    + f_8 * smi0_641[k]
                    - f_9 * smi1_641[k]
                    + f_3 * pc_y[k] * smk_824[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, pc_y, slk_609, slk_610, slk_611, smi0_642, \
                         smi0_643, smi1_642, smi1_643, smk_825, smk_826, \
                         smk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_19 * slk_609[k]
                    + f_10 * smi0_642[k]
                    - f_11 * smi1_642[k]
                    + f_3 * pc_y[k] * smk_825[k];

        t_1032[k] = f_19 * slk_610[k]
                    + f_12 * smi0_643[k]
                    - f_13 * smi1_643[k]
                    + f_3 * pc_y[k] * smk_826[k];

        t_1033[k] = f_19 * slk_611[k]
                    + f_3 * pc_y[k] * smk_827[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, pc_x, pc_y, pc_z, slk_575, slk_612, slk_828, \
                         smi0_643, smi0_644, smi1_643, smi1_644, smk_827, \
                         smk_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_15 * slk_575[k]
                    + f_1 * smi0_643[k]
                    - f_2 * smi1_643[k]
                    + f_3 * pc_z[k] * smk_827[k];

        t_1035[k] = f_17 * slk_828[k]
                    + f_1 * smi0_644[k]
                    - f_2 * smi1_644[k]
                    + f_3 * pc_x[k] * smk_828[k];

        t_1036[k] = f_18 * slk_612[k]
                    + f_3 * pc_y[k] * smk_828[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, pc_x, pc_y, pc_z, slk_576, slk_614, slk_831, \
                         smi0_647, smi1_647, smk_828, smk_830, \
                         smk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_16 * slk_576[k]
                    + f_3 * pc_z[k] * smk_828[k];

        t_1038[k] = f_17 * slk_831[k]
                    + f_4 * smi0_647[k]
                    - f_5 * smi1_647[k]
                    + f_3 * pc_x[k] * smk_831[k];

        t_1039[k] = f_18 * slk_614[k]
                    + f_3 * pc_y[k] * smk_830[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pc_x, pc_z, slk_579, slk_833, slk_834, \
                         smi0_649, smi0_650, smi1_649, smi1_650, smk_831, smk_833, \
                         smk_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_17 * slk_833[k]
                    + f_4 * smi0_649[k]
                    - f_5 * smi1_649[k]
                    + f_3 * pc_x[k] * smk_833[k];

        t_1041[k] = f_17 * slk_834[k]
                    + f_6 * smi0_650[k]
                    - f_7 * smi1_650[k]
                    + f_3 * pc_x[k] * smk_834[k];

        t_1042[k] = f_16 * slk_579[k]
                    + f_3 * pc_z[k] * smk_831[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pc_x, pc_y, slk_617, slk_837, slk_838, \
                         smi0_653, smi0_654, smi1_653, smi1_654, smk_833, smk_837, \
                         smk_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_18 * slk_617[k]
                    + f_3 * pc_y[k] * smk_833[k];

        t_1044[k] = f_17 * slk_837[k]
                    + f_6 * smi0_653[k]
                    - f_7 * smi1_653[k]
                    + f_3 * pc_x[k] * smk_837[k];

        t_1045[k] = f_17 * slk_838[k]
                    + f_8 * smi0_654[k]
                    - f_9 * smi1_654[k]
                    + f_3 * pc_x[k] * smk_838[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pc_x, pc_y, pc_z, slk_582, slk_621, slk_840, \
                         smi0_656, smi1_656, smk_834, smk_837, \
                         smk_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_16 * slk_582[k]
                    + f_3 * pc_z[k] * smk_834[k];

        t_1047[k] = f_17 * slk_840[k]
                    + f_8 * smi0_656[k]
                    - f_9 * smi1_656[k]
                    + f_3 * pc_x[k] * smk_840[k];

        t_1048[k] = f_18 * slk_621[k]
                    + f_3 * pc_y[k] * smk_837[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pc_x, pc_z, slk_586, slk_842, slk_843, \
                         smi0_658, smi0_659, smi1_658, smi1_659, smk_838, smk_842, \
                         smk_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_17 * slk_842[k]
                    + f_8 * smi0_658[k]
                    - f_9 * smi1_658[k]
                    + f_3 * pc_x[k] * smk_842[k];

        t_1050[k] = f_17 * slk_843[k]
                    + f_10 * smi0_659[k]
                    - f_11 * smi1_659[k]
                    + f_3 * pc_x[k] * smk_843[k];

        t_1051[k] = f_16 * slk_586[k]
                    + f_3 * pc_z[k] * smk_838[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pc_x, pc_y, slk_626, slk_845, slk_846, \
                         smi0_661, smi0_662, smi1_661, smi1_662, smk_842, smk_845, \
                         smk_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_17 * slk_845[k]
                    + f_10 * smi0_661[k]
                    - f_11 * smi1_661[k]
                    + f_3 * pc_x[k] * smk_845[k];

        t_1053[k] = f_17 * slk_846[k]
                    + f_10 * smi0_662[k]
                    - f_11 * smi1_662[k]
                    + f_3 * pc_x[k] * smk_846[k];

        t_1054[k] = f_18 * slk_626[k]
                    + f_3 * pc_y[k] * smk_842[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, pc_z, slk_591, slk_848, slk_849, \
                         smi0_664, smi0_665, smi1_664, smi1_665, smk_843, smk_848, \
                         smk_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_17 * slk_848[k]
                    + f_10 * smi0_664[k]
                    - f_11 * smi1_664[k]
                    + f_3 * pc_x[k] * smk_848[k];

        t_1056[k] = f_17 * slk_849[k]
                    + f_12 * smi0_665[k]
                    - f_13 * smi1_665[k]
                    + f_3 * pc_x[k] * smk_849[k];

        t_1057[k] = f_16 * slk_591[k]
                    + f_3 * pc_z[k] * smk_843[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_x, slk_851, slk_852, slk_853, smi0_667, \
                         smi0_668, smi0_669, smi1_667, smi1_668, smi1_669, smk_851, smk_852, \
                         smk_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_17 * slk_851[k]
                    + f_12 * smi0_667[k]
                    - f_13 * smi1_667[k]
                    + f_3 * pc_x[k] * smk_851[k];

        t_1059[k] = f_17 * slk_852[k]
                    + f_12 * smi0_668[k]
                    - f_13 * smi1_668[k]
                    + f_3 * pc_x[k] * smk_852[k];

        t_1060[k] = f_17 * slk_853[k]
                    + f_12 * smi0_669[k]
                    - f_13 * smi1_669[k]
                    + f_3 * pc_x[k] * smk_853[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, pc_x, pc_y, slk_632, slk_855, \
                         slk_856, slk_857, smi0_671, smi1_671, smk_848, smk_855, smk_856, \
                         smk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_18 * slk_632[k]
                    + f_3 * pc_y[k] * smk_848[k];

        t_1062[k] = f_17 * slk_855[k]
                    + f_12 * smi0_671[k]
                    - f_13 * smi1_671[k]
                    + f_3 * pc_x[k] * smk_855[k];

        t_1063[k] = f_17 * slk_856[k]
                    + f_3 * pc_x[k] * smk_856[k];

        t_1064[k] = f_17 * slk_857[k]
                    + f_3 * pc_x[k] * smk_857[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, pc_x, slk_858, slk_859, \
                         slk_860, slk_861, slk_862, smk_858, smk_859, smk_860, smk_861, \
                         smk_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_17 * slk_858[k]
                    + f_3 * pc_x[k] * smk_858[k];

        t_1066[k] = f_17 * slk_859[k]
                    + f_3 * pc_x[k] * smk_859[k];

        t_1067[k] = f_17 * slk_860[k]
                    + f_3 * pc_x[k] * smk_860[k];

        t_1068[k] = f_17 * slk_861[k]
                    + f_3 * pc_x[k] * smk_861[k];

        t_1069[k] = f_17 * slk_862[k]
                    + f_3 * pc_x[k] * smk_862[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, slk_604, slk_640, slk_863, \
                         smi0_665, smi1_665, smk_856, smk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_17 * slk_863[k]
                    + f_3 * pc_x[k] * smk_863[k];

        t_1071[k] = f_18 * slk_640[k]
                    + f_1 * smi0_665[k]
                    - f_2 * smi1_665[k]
                    + f_3 * pc_y[k] * smk_856[k];

        t_1072[k] = f_16 * slk_604[k]
                    + f_3 * pc_z[k] * smk_856[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_y, slk_642, slk_643, slk_644, smi0_667, \
                         smi0_668, smi0_669, smi1_667, smi1_668, smi1_669, smk_858, smk_859, \
                         smk_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_18 * slk_642[k]
                    + f_4 * smi0_667[k]
                    - f_5 * smi1_667[k]
                    + f_3 * pc_y[k] * smk_858[k];

        t_1074[k] = f_18 * slk_643[k]
                    + f_6 * smi0_668[k]
                    - f_7 * smi1_668[k]
                    + f_3 * pc_y[k] * smk_859[k];

        t_1075[k] = f_18 * slk_644[k]
                    + f_8 * smi0_669[k]
                    - f_9 * smi1_669[k]
                    + f_3 * pc_y[k] * smk_860[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_y, slk_645, slk_646, slk_647, smi0_670, \
                         smi0_671, smi1_670, smi1_671, smk_861, smk_862, \
                         smk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_18 * slk_645[k]
                    + f_10 * smi0_670[k]
                    - f_11 * smi1_670[k]
                    + f_3 * pc_y[k] * smk_861[k];

        t_1077[k] = f_18 * slk_646[k]
                    + f_12 * smi0_671[k]
                    - f_13 * smi1_671[k]
                    + f_3 * pc_y[k] * smk_862[k];

        t_1078[k] = f_18 * slk_647[k]
                    + f_3 * pc_y[k] * smk_863[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, pc_x, pc_y, pc_z, slk_611, slk_648, slk_864, \
                         smi0_671, smi0_672, smi1_671, smi1_672, smk_863, \
                         smk_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_16 * slk_611[k]
                    + f_1 * smi0_671[k]
                    - f_2 * smi1_671[k]
                    + f_3 * pc_z[k] * smk_863[k];

        t_1080[k] = f_17 * slk_864[k]
                    + f_1 * smi0_672[k]
                    - f_2 * smi1_672[k]
                    + f_3 * pc_x[k] * smk_864[k];

        t_1081[k] = f_17 * slk_648[k]
                    + f_3 * pc_y[k] * smk_864[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, pc_x, pc_y, pc_z, slk_612, slk_650, slk_867, \
                         smi0_675, smi1_675, smk_864, smk_866, \
                         smk_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_17 * slk_612[k]
                    + f_3 * pc_z[k] * smk_864[k];

        t_1083[k] = f_17 * slk_867[k]
                    + f_4 * smi0_675[k]
                    - f_5 * smi1_675[k]
                    + f_3 * pc_x[k] * smk_867[k];

        t_1084[k] = f_17 * slk_650[k]
                    + f_3 * pc_y[k] * smk_866[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_x, pc_z, slk_615, slk_869, slk_870, \
                         smi0_677, smi0_678, smi1_677, smi1_678, smk_867, smk_869, \
                         smk_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_17 * slk_869[k]
                    + f_4 * smi0_677[k]
                    - f_5 * smi1_677[k]
                    + f_3 * pc_x[k] * smk_869[k];

        t_1086[k] = f_17 * slk_870[k]
                    + f_6 * smi0_678[k]
                    - f_7 * smi1_678[k]
                    + f_3 * pc_x[k] * smk_870[k];

        t_1087[k] = f_17 * slk_615[k]
                    + f_3 * pc_z[k] * smk_867[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_x, pc_y, slk_653, slk_873, slk_874, \
                         smi0_681, smi0_682, smi1_681, smi1_682, smk_869, smk_873, \
                         smk_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * slk_653[k]
                    + f_3 * pc_y[k] * smk_869[k];

        t_1089[k] = f_17 * slk_873[k]
                    + f_6 * smi0_681[k]
                    - f_7 * smi1_681[k]
                    + f_3 * pc_x[k] * smk_873[k];

        t_1090[k] = f_17 * slk_874[k]
                    + f_8 * smi0_682[k]
                    - f_9 * smi1_682[k]
                    + f_3 * pc_x[k] * smk_874[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_y, pc_z, slk_618, slk_657, slk_876, \
                         smi0_684, smi1_684, smk_870, smk_873, \
                         smk_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_17 * slk_618[k]
                    + f_3 * pc_z[k] * smk_870[k];

        t_1092[k] = f_17 * slk_876[k]
                    + f_8 * smi0_684[k]
                    - f_9 * smi1_684[k]
                    + f_3 * pc_x[k] * smk_876[k];

        t_1093[k] = f_17 * slk_657[k]
                    + f_3 * pc_y[k] * smk_873[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, pc_z, slk_622, slk_878, slk_879, \
                         smi0_686, smi0_687, smi1_686, smi1_687, smk_874, smk_878, \
                         smk_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_17 * slk_878[k]
                    + f_8 * smi0_686[k]
                    - f_9 * smi1_686[k]
                    + f_3 * pc_x[k] * smk_878[k];

        t_1095[k] = f_17 * slk_879[k]
                    + f_10 * smi0_687[k]
                    - f_11 * smi1_687[k]
                    + f_3 * pc_x[k] * smk_879[k];

        t_1096[k] = f_17 * slk_622[k]
                    + f_3 * pc_z[k] * smk_874[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_y, slk_662, slk_881, slk_882, \
                         smi0_689, smi0_690, smi1_689, smi1_690, smk_878, smk_881, \
                         smk_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_17 * slk_881[k]
                    + f_10 * smi0_689[k]
                    - f_11 * smi1_689[k]
                    + f_3 * pc_x[k] * smk_881[k];

        t_1098[k] = f_17 * slk_882[k]
                    + f_10 * smi0_690[k]
                    - f_11 * smi1_690[k]
                    + f_3 * pc_x[k] * smk_882[k];

        t_1099[k] = f_17 * slk_662[k]
                    + f_3 * pc_y[k] * smk_878[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_z, slk_627, slk_884, slk_885, \
                         smi0_692, smi0_693, smi1_692, smi1_693, smk_879, smk_884, \
                         smk_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_17 * slk_884[k]
                    + f_10 * smi0_692[k]
                    - f_11 * smi1_692[k]
                    + f_3 * pc_x[k] * smk_884[k];

        t_1101[k] = f_17 * slk_885[k]
                    + f_12 * smi0_693[k]
                    - f_13 * smi1_693[k]
                    + f_3 * pc_x[k] * smk_885[k];

        t_1102[k] = f_17 * slk_627[k]
                    + f_3 * pc_z[k] * smk_879[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smi0, const size_t smi1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_900 = buffer.data(sll0 + 900);
    const auto *sll0_903 = buffer.data(sll0 + 903);
    const auto *sll0_905 = buffer.data(sll0 + 905);
    const auto *sll0_906 = buffer.data(sll0 + 906);
    const auto *sll0_909 = buffer.data(sll0 + 909);
    const auto *sll0_910 = buffer.data(sll0 + 910);
    const auto *sll0_912 = buffer.data(sll0 + 912);
    const auto *sll0_914 = buffer.data(sll0 + 914);
    const auto *sll0_915 = buffer.data(sll0 + 915);
    const auto *sll0_917 = buffer.data(sll0 + 917);
    const auto *sll0_918 = buffer.data(sll0 + 918);
    const auto *sll0_920 = buffer.data(sll0 + 920);
    const auto *sll0_921 = buffer.data(sll0 + 921);
    const auto *sll0_923 = buffer.data(sll0 + 923);
    const auto *sll0_924 = buffer.data(sll0 + 924);
    const auto *sll0_925 = buffer.data(sll0 + 925);
    const auto *sll0_927 = buffer.data(sll0 + 927);

    const auto *slk_640 = buffer.data(slk + 640);
    const auto *slk_647 = buffer.data(slk + 647);
    const auto *slk_648 = buffer.data(slk + 648);
    const auto *slk_651 = buffer.data(slk + 651);
    const auto *slk_654 = buffer.data(slk + 654);
    const auto *slk_658 = buffer.data(slk + 658);
    const auto *slk_663 = buffer.data(slk + 663);
    const auto *slk_668 = buffer.data(slk + 668);
    const auto *slk_676 = buffer.data(slk + 676);
    const auto *slk_678 = buffer.data(slk + 678);
    const auto *slk_679 = buffer.data(slk + 679);
    const auto *slk_680 = buffer.data(slk + 680);
    const auto *slk_681 = buffer.data(slk + 681);
    const auto *slk_682 = buffer.data(slk + 682);
    const auto *slk_683 = buffer.data(slk + 683);
    const auto *slk_684 = buffer.data(slk + 684);
    const auto *slk_686 = buffer.data(slk + 686);
    const auto *slk_687 = buffer.data(slk + 687);
    const auto *slk_689 = buffer.data(slk + 689);
    const auto *slk_690 = buffer.data(slk + 690);
    const auto *slk_693 = buffer.data(slk + 693);
    const auto *slk_694 = buffer.data(slk + 694);
    const auto *slk_698 = buffer.data(slk + 698);
    const auto *slk_699 = buffer.data(slk + 699);
    const auto *slk_704 = buffer.data(slk + 704);
    const auto *slk_712 = buffer.data(slk + 712);
    const auto *slk_714 = buffer.data(slk + 714);
    const auto *slk_715 = buffer.data(slk + 715);
    const auto *slk_716 = buffer.data(slk + 716);
    const auto *slk_717 = buffer.data(slk + 717);
    const auto *slk_718 = buffer.data(slk + 718);
    const auto *slk_719 = buffer.data(slk + 719);
    const auto *slk_720 = buffer.data(slk + 720);
    const auto *slk_721 = buffer.data(slk + 721);
    const auto *slk_722 = buffer.data(slk + 722);
    const auto *slk_723 = buffer.data(slk + 723);
    const auto *slk_725 = buffer.data(slk + 725);
    const auto *slk_726 = buffer.data(slk + 726);
    const auto *slk_728 = buffer.data(slk + 728);
    const auto *slk_729 = buffer.data(slk + 729);
    const auto *slk_730 = buffer.data(slk + 730);
    const auto *slk_732 = buffer.data(slk + 732);
    const auto *slk_733 = buffer.data(slk + 733);
    const auto *slk_734 = buffer.data(slk + 734);
    const auto *slk_735 = buffer.data(slk + 735);
    const auto *slk_737 = buffer.data(slk + 737);
    const auto *slk_738 = buffer.data(slk + 738);
    const auto *slk_739 = buffer.data(slk + 739);
    const auto *slk_740 = buffer.data(slk + 740);
    const auto *slk_748 = buffer.data(slk + 748);
    const auto *slk_750 = buffer.data(slk + 750);
    const auto *slk_751 = buffer.data(slk + 751);
    const auto *slk_752 = buffer.data(slk + 752);
    const auto *slk_887 = buffer.data(slk + 887);
    const auto *slk_888 = buffer.data(slk + 888);
    const auto *slk_889 = buffer.data(slk + 889);
    const auto *slk_891 = buffer.data(slk + 891);
    const auto *slk_892 = buffer.data(slk + 892);
    const auto *slk_893 = buffer.data(slk + 893);
    const auto *slk_894 = buffer.data(slk + 894);
    const auto *slk_895 = buffer.data(slk + 895);
    const auto *slk_896 = buffer.data(slk + 896);
    const auto *slk_897 = buffer.data(slk + 897);
    const auto *slk_898 = buffer.data(slk + 898);
    const auto *slk_899 = buffer.data(slk + 899);
    const auto *slk_900 = buffer.data(slk + 900);
    const auto *slk_903 = buffer.data(slk + 903);
    const auto *slk_905 = buffer.data(slk + 905);
    const auto *slk_906 = buffer.data(slk + 906);
    const auto *slk_909 = buffer.data(slk + 909);
    const auto *slk_910 = buffer.data(slk + 910);
    const auto *slk_912 = buffer.data(slk + 912);
    const auto *slk_914 = buffer.data(slk + 914);
    const auto *slk_915 = buffer.data(slk + 915);
    const auto *slk_917 = buffer.data(slk + 917);
    const auto *slk_918 = buffer.data(slk + 918);
    const auto *slk_920 = buffer.data(slk + 920);
    const auto *slk_921 = buffer.data(slk + 921);
    const auto *slk_923 = buffer.data(slk + 923);
    const auto *slk_924 = buffer.data(slk + 924);
    const auto *slk_925 = buffer.data(slk + 925);
    const auto *slk_927 = buffer.data(slk + 927);
    const auto *slk_928 = buffer.data(slk + 928);
    const auto *slk_929 = buffer.data(slk + 929);
    const auto *slk_930 = buffer.data(slk + 930);
    const auto *slk_931 = buffer.data(slk + 931);
    const auto *slk_932 = buffer.data(slk + 932);
    const auto *slk_933 = buffer.data(slk + 933);
    const auto *slk_934 = buffer.data(slk + 934);
    const auto *slk_935 = buffer.data(slk + 935);
    const auto *slk_964 = buffer.data(slk + 964);
    const auto *slk_965 = buffer.data(slk + 965);
    const auto *slk_966 = buffer.data(slk + 966);
    const auto *slk_967 = buffer.data(slk + 967);
    const auto *slk_968 = buffer.data(slk + 968);
    const auto *slk_969 = buffer.data(slk + 969);
    const auto *slk_970 = buffer.data(slk + 970);
    const auto *slk_971 = buffer.data(slk + 971);

    const auto *sll1_900 = buffer.data(sll1 + 900);
    const auto *sll1_903 = buffer.data(sll1 + 903);
    const auto *sll1_905 = buffer.data(sll1 + 905);
    const auto *sll1_906 = buffer.data(sll1 + 906);
    const auto *sll1_909 = buffer.data(sll1 + 909);
    const auto *sll1_910 = buffer.data(sll1 + 910);
    const auto *sll1_912 = buffer.data(sll1 + 912);
    const auto *sll1_914 = buffer.data(sll1 + 914);
    const auto *sll1_915 = buffer.data(sll1 + 915);
    const auto *sll1_917 = buffer.data(sll1 + 917);
    const auto *sll1_918 = buffer.data(sll1 + 918);
    const auto *sll1_920 = buffer.data(sll1 + 920);
    const auto *sll1_921 = buffer.data(sll1 + 921);
    const auto *sll1_923 = buffer.data(sll1 + 923);
    const auto *sll1_924 = buffer.data(sll1 + 924);
    const auto *sll1_925 = buffer.data(sll1 + 925);
    const auto *sll1_927 = buffer.data(sll1 + 927);

    const auto *smi0_693 = buffer.data(smi0 + 693);
    const auto *smi0_695 = buffer.data(smi0 + 695);
    const auto *smi0_696 = buffer.data(smi0 + 696);
    const auto *smi0_697 = buffer.data(smi0 + 697);
    const auto *smi0_698 = buffer.data(smi0 + 698);
    const auto *smi0_699 = buffer.data(smi0 + 699);
    const auto *smi0_700 = buffer.data(smi0 + 700);
    const auto *smi0_703 = buffer.data(smi0 + 703);
    const auto *smi0_705 = buffer.data(smi0 + 705);
    const auto *smi0_706 = buffer.data(smi0 + 706);
    const auto *smi0_709 = buffer.data(smi0 + 709);
    const auto *smi0_710 = buffer.data(smi0 + 710);
    const auto *smi0_712 = buffer.data(smi0 + 712);
    const auto *smi0_714 = buffer.data(smi0 + 714);
    const auto *smi0_715 = buffer.data(smi0 + 715);
    const auto *smi0_717 = buffer.data(smi0 + 717);
    const auto *smi0_718 = buffer.data(smi0 + 718);
    const auto *smi0_720 = buffer.data(smi0 + 720);
    const auto *smi0_721 = buffer.data(smi0 + 721);
    const auto *smi0_723 = buffer.data(smi0 + 723);
    const auto *smi0_724 = buffer.data(smi0 + 724);
    const auto *smi0_725 = buffer.data(smi0 + 725);
    const auto *smi0_726 = buffer.data(smi0 + 726);
    const auto *smi0_727 = buffer.data(smi0 + 727);
    const auto *smi0_749 = buffer.data(smi0 + 749);
    const auto *smi0_751 = buffer.data(smi0 + 751);
    const auto *smi0_752 = buffer.data(smi0 + 752);
    const auto *smi0_753 = buffer.data(smi0 + 753);

    const auto *smi1_693 = buffer.data(smi1 + 693);
    const auto *smi1_695 = buffer.data(smi1 + 695);
    const auto *smi1_696 = buffer.data(smi1 + 696);
    const auto *smi1_697 = buffer.data(smi1 + 697);
    const auto *smi1_698 = buffer.data(smi1 + 698);
    const auto *smi1_699 = buffer.data(smi1 + 699);
    const auto *smi1_700 = buffer.data(smi1 + 700);
    const auto *smi1_703 = buffer.data(smi1 + 703);
    const auto *smi1_705 = buffer.data(smi1 + 705);
    const auto *smi1_706 = buffer.data(smi1 + 706);
    const auto *smi1_709 = buffer.data(smi1 + 709);
    const auto *smi1_710 = buffer.data(smi1 + 710);
    const auto *smi1_712 = buffer.data(smi1 + 712);
    const auto *smi1_714 = buffer.data(smi1 + 714);
    const auto *smi1_715 = buffer.data(smi1 + 715);
    const auto *smi1_717 = buffer.data(smi1 + 717);
    const auto *smi1_718 = buffer.data(smi1 + 718);
    const auto *smi1_720 = buffer.data(smi1 + 720);
    const auto *smi1_721 = buffer.data(smi1 + 721);
    const auto *smi1_723 = buffer.data(smi1 + 723);
    const auto *smi1_724 = buffer.data(smi1 + 724);
    const auto *smi1_725 = buffer.data(smi1 + 725);
    const auto *smi1_726 = buffer.data(smi1 + 726);
    const auto *smi1_727 = buffer.data(smi1 + 727);
    const auto *smi1_749 = buffer.data(smi1 + 749);
    const auto *smi1_751 = buffer.data(smi1 + 751);
    const auto *smi1_752 = buffer.data(smi1 + 752);
    const auto *smi1_753 = buffer.data(smi1 + 753);

    const auto *smk_884 = buffer.data(smk + 884);
    const auto *smk_887 = buffer.data(smk + 887);
    const auto *smk_888 = buffer.data(smk + 888);
    const auto *smk_889 = buffer.data(smk + 889);
    const auto *smk_891 = buffer.data(smk + 891);
    const auto *smk_892 = buffer.data(smk + 892);
    const auto *smk_893 = buffer.data(smk + 893);
    const auto *smk_894 = buffer.data(smk + 894);
    const auto *smk_895 = buffer.data(smk + 895);
    const auto *smk_896 = buffer.data(smk + 896);
    const auto *smk_897 = buffer.data(smk + 897);
    const auto *smk_898 = buffer.data(smk + 898);
    const auto *smk_899 = buffer.data(smk + 899);
    const auto *smk_900 = buffer.data(smk + 900);
    const auto *smk_902 = buffer.data(smk + 902);
    const auto *smk_903 = buffer.data(smk + 903);
    const auto *smk_905 = buffer.data(smk + 905);
    const auto *smk_906 = buffer.data(smk + 906);
    const auto *smk_909 = buffer.data(smk + 909);
    const auto *smk_910 = buffer.data(smk + 910);
    const auto *smk_912 = buffer.data(smk + 912);
    const auto *smk_914 = buffer.data(smk + 914);
    const auto *smk_915 = buffer.data(smk + 915);
    const auto *smk_917 = buffer.data(smk + 917);
    const auto *smk_918 = buffer.data(smk + 918);
    const auto *smk_920 = buffer.data(smk + 920);
    const auto *smk_921 = buffer.data(smk + 921);
    const auto *smk_923 = buffer.data(smk + 923);
    const auto *smk_924 = buffer.data(smk + 924);
    const auto *smk_925 = buffer.data(smk + 925);
    const auto *smk_927 = buffer.data(smk + 927);
    const auto *smk_928 = buffer.data(smk + 928);
    const auto *smk_929 = buffer.data(smk + 929);
    const auto *smk_930 = buffer.data(smk + 930);
    const auto *smk_931 = buffer.data(smk + 931);
    const auto *smk_932 = buffer.data(smk + 932);
    const auto *smk_933 = buffer.data(smk + 933);
    const auto *smk_934 = buffer.data(smk + 934);
    const auto *smk_935 = buffer.data(smk + 935);
    const auto *smk_936 = buffer.data(smk + 936);
    const auto *smk_938 = buffer.data(smk + 938);
    const auto *smk_939 = buffer.data(smk + 939);
    const auto *smk_941 = buffer.data(smk + 941);
    const auto *smk_942 = buffer.data(smk + 942);
    const auto *smk_945 = buffer.data(smk + 945);
    const auto *smk_946 = buffer.data(smk + 946);
    const auto *smk_950 = buffer.data(smk + 950);
    const auto *smk_951 = buffer.data(smk + 951);
    const auto *smk_956 = buffer.data(smk + 956);
    const auto *smk_964 = buffer.data(smk + 964);
    const auto *smk_965 = buffer.data(smk + 965);
    const auto *smk_966 = buffer.data(smk + 966);
    const auto *smk_967 = buffer.data(smk + 967);
    const auto *smk_968 = buffer.data(smk + 968);
    const auto *smk_969 = buffer.data(smk + 969);
    const auto *smk_970 = buffer.data(smk + 970);
    const auto *smk_971 = buffer.data(smk + 971);

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, slk_887, slk_888, slk_889, smi0_695, \
                         smi0_696, smi0_697, smi1_695, smi1_696, smi1_697, smk_887, smk_888, \
                         smk_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_17 * slk_887[k]
                    + f_12 * smi0_695[k]
                    - f_13 * smi1_695[k]
                    + f_3 * pc_x[k] * smk_887[k];

        t_1104[k] = f_17 * slk_888[k]
                    + f_12 * smi0_696[k]
                    - f_13 * smi1_696[k]
                    + f_3 * pc_x[k] * smk_888[k];

        t_1105[k] = f_17 * slk_889[k]
                    + f_12 * smi0_697[k]
                    - f_13 * smi1_697[k]
                    + f_3 * pc_x[k] * smk_889[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, pc_y, slk_668, slk_891, \
                         slk_892, slk_893, smi0_699, smi1_699, smk_884, smk_891, smk_892, \
                         smk_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_17 * slk_668[k]
                    + f_3 * pc_y[k] * smk_884[k];

        t_1107[k] = f_17 * slk_891[k]
                    + f_12 * smi0_699[k]
                    - f_13 * smi1_699[k]
                    + f_3 * pc_x[k] * smk_891[k];

        t_1108[k] = f_17 * slk_892[k]
                    + f_3 * pc_x[k] * smk_892[k];

        t_1109[k] = f_17 * slk_893[k]
                    + f_3 * pc_x[k] * smk_893[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, pc_x, slk_894, slk_895, \
                         slk_896, slk_897, slk_898, smk_894, smk_895, smk_896, smk_897, \
                         smk_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_17 * slk_894[k]
                    + f_3 * pc_x[k] * smk_894[k];

        t_1111[k] = f_17 * slk_895[k]
                    + f_3 * pc_x[k] * smk_895[k];

        t_1112[k] = f_17 * slk_896[k]
                    + f_3 * pc_x[k] * smk_896[k];

        t_1113[k] = f_17 * slk_897[k]
                    + f_3 * pc_x[k] * smk_897[k];

        t_1114[k] = f_17 * slk_898[k]
                    + f_3 * pc_x[k] * smk_898[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, pc_x, pc_y, pc_z, slk_640, slk_676, slk_899, \
                         smi0_693, smi1_693, smk_892, smk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_17 * slk_899[k]
                    + f_3 * pc_x[k] * smk_899[k];

        t_1116[k] = f_17 * slk_676[k]
                    + f_1 * smi0_693[k]
                    - f_2 * smi1_693[k]
                    + f_3 * pc_y[k] * smk_892[k];

        t_1117[k] = f_17 * slk_640[k]
                    + f_3 * pc_z[k] * smk_892[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pc_y, slk_678, slk_679, slk_680, smi0_695, \
                         smi0_696, smi0_697, smi1_695, smi1_696, smi1_697, smk_894, smk_895, \
                         smk_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_17 * slk_678[k]
                    + f_4 * smi0_695[k]
                    - f_5 * smi1_695[k]
                    + f_3 * pc_y[k] * smk_894[k];

        t_1119[k] = f_17 * slk_679[k]
                    + f_6 * smi0_696[k]
                    - f_7 * smi1_696[k]
                    + f_3 * pc_y[k] * smk_895[k];

        t_1120[k] = f_17 * slk_680[k]
                    + f_8 * smi0_697[k]
                    - f_9 * smi1_697[k]
                    + f_3 * pc_y[k] * smk_896[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_y, slk_681, slk_682, slk_683, smi0_698, \
                         smi0_699, smi1_698, smi1_699, smk_897, smk_898, \
                         smk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_17 * slk_681[k]
                    + f_10 * smi0_698[k]
                    - f_11 * smi1_698[k]
                    + f_3 * pc_y[k] * smk_897[k];

        t_1122[k] = f_17 * slk_682[k]
                    + f_12 * smi0_699[k]
                    - f_13 * smi1_699[k]
                    + f_3 * pc_y[k] * smk_898[k];

        t_1123[k] = f_17 * slk_683[k]
                    + f_3 * pc_y[k] * smk_899[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_y, pc_z, slk_647, slk_684, slk_900, \
                         smi0_699, smi0_700, smi1_699, smi1_700, smk_899, \
                         smk_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_17 * slk_647[k]
                    + f_1 * smi0_699[k]
                    - f_2 * smi1_699[k]
                    + f_3 * pc_z[k] * smk_899[k];

        t_1125[k] = f_17 * slk_900[k]
                    + f_1 * smi0_700[k]
                    - f_2 * smi1_700[k]
                    + f_3 * pc_x[k] * smk_900[k];

        t_1126[k] = f_16 * slk_684[k]
                    + f_3 * pc_y[k] * smk_900[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, pc_y, pc_z, slk_648, slk_686, slk_903, \
                         smi0_703, smi1_703, smk_900, smk_902, \
                         smk_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_18 * slk_648[k]
                    + f_3 * pc_z[k] * smk_900[k];

        t_1128[k] = f_17 * slk_903[k]
                    + f_4 * smi0_703[k]
                    - f_5 * smi1_703[k]
                    + f_3 * pc_x[k] * smk_903[k];

        t_1129[k] = f_16 * slk_686[k]
                    + f_3 * pc_y[k] * smk_902[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, pc_z, slk_651, slk_905, slk_906, \
                         smi0_705, smi0_706, smi1_705, smi1_706, smk_903, smk_905, \
                         smk_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_17 * slk_905[k]
                    + f_4 * smi0_705[k]
                    - f_5 * smi1_705[k]
                    + f_3 * pc_x[k] * smk_905[k];

        t_1131[k] = f_17 * slk_906[k]
                    + f_6 * smi0_706[k]
                    - f_7 * smi1_706[k]
                    + f_3 * pc_x[k] * smk_906[k];

        t_1132[k] = f_18 * slk_651[k]
                    + f_3 * pc_z[k] * smk_903[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, pc_y, slk_689, slk_909, slk_910, \
                         smi0_709, smi0_710, smi1_709, smi1_710, smk_905, smk_909, \
                         smk_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_16 * slk_689[k]
                    + f_3 * pc_y[k] * smk_905[k];

        t_1134[k] = f_17 * slk_909[k]
                    + f_6 * smi0_709[k]
                    - f_7 * smi1_709[k]
                    + f_3 * pc_x[k] * smk_909[k];

        t_1135[k] = f_17 * slk_910[k]
                    + f_8 * smi0_710[k]
                    - f_9 * smi1_710[k]
                    + f_3 * pc_x[k] * smk_910[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, pc_x, pc_y, pc_z, slk_654, slk_693, slk_912, \
                         smi0_712, smi1_712, smk_906, smk_909, \
                         smk_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_18 * slk_654[k]
                    + f_3 * pc_z[k] * smk_906[k];

        t_1137[k] = f_17 * slk_912[k]
                    + f_8 * smi0_712[k]
                    - f_9 * smi1_712[k]
                    + f_3 * pc_x[k] * smk_912[k];

        t_1138[k] = f_16 * slk_693[k]
                    + f_3 * pc_y[k] * smk_909[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, pc_x, pc_z, slk_658, slk_914, slk_915, \
                         smi0_714, smi0_715, smi1_714, smi1_715, smk_910, smk_914, \
                         smk_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_17 * slk_914[k]
                    + f_8 * smi0_714[k]
                    - f_9 * smi1_714[k]
                    + f_3 * pc_x[k] * smk_914[k];

        t_1140[k] = f_17 * slk_915[k]
                    + f_10 * smi0_715[k]
                    - f_11 * smi1_715[k]
                    + f_3 * pc_x[k] * smk_915[k];

        t_1141[k] = f_18 * slk_658[k]
                    + f_3 * pc_z[k] * smk_910[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, pc_x, pc_y, slk_698, slk_917, slk_918, \
                         smi0_717, smi0_718, smi1_717, smi1_718, smk_914, smk_917, \
                         smk_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_17 * slk_917[k]
                    + f_10 * smi0_717[k]
                    - f_11 * smi1_717[k]
                    + f_3 * pc_x[k] * smk_917[k];

        t_1143[k] = f_17 * slk_918[k]
                    + f_10 * smi0_718[k]
                    - f_11 * smi1_718[k]
                    + f_3 * pc_x[k] * smk_918[k];

        t_1144[k] = f_16 * slk_698[k]
                    + f_3 * pc_y[k] * smk_914[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, pc_x, pc_z, slk_663, slk_920, slk_921, \
                         smi0_720, smi0_721, smi1_720, smi1_721, smk_915, smk_920, \
                         smk_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_17 * slk_920[k]
                    + f_10 * smi0_720[k]
                    - f_11 * smi1_720[k]
                    + f_3 * pc_x[k] * smk_920[k];

        t_1146[k] = f_17 * slk_921[k]
                    + f_12 * smi0_721[k]
                    - f_13 * smi1_721[k]
                    + f_3 * pc_x[k] * smk_921[k];

        t_1147[k] = f_18 * slk_663[k]
                    + f_3 * pc_z[k] * smk_915[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, pc_x, slk_923, slk_924, slk_925, smi0_723, \
                         smi0_724, smi0_725, smi1_723, smi1_724, smi1_725, smk_923, smk_924, \
                         smk_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = f_17 * slk_923[k]
                    + f_12 * smi0_723[k]
                    - f_13 * smi1_723[k]
                    + f_3 * pc_x[k] * smk_923[k];

        t_1149[k] = f_17 * slk_924[k]
                    + f_12 * smi0_724[k]
                    - f_13 * smi1_724[k]
                    + f_3 * pc_x[k] * smk_924[k];

        t_1150[k] = f_17 * slk_925[k]
                    + f_12 * smi0_725[k]
                    - f_13 * smi1_725[k]
                    + f_3 * pc_x[k] * smk_925[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, pc_x, pc_y, slk_704, slk_927, \
                         slk_928, slk_929, smi0_727, smi1_727, smk_920, smk_927, smk_928, \
                         smk_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_16 * slk_704[k]
                    + f_3 * pc_y[k] * smk_920[k];

        t_1152[k] = f_17 * slk_927[k]
                    + f_12 * smi0_727[k]
                    - f_13 * smi1_727[k]
                    + f_3 * pc_x[k] * smk_927[k];

        t_1153[k] = f_17 * slk_928[k]
                    + f_3 * pc_x[k] * smk_928[k];

        t_1154[k] = f_17 * slk_929[k]
                    + f_3 * pc_x[k] * smk_929[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, pc_x, slk_930, slk_931, \
                         slk_932, slk_933, slk_934, smk_930, smk_931, smk_932, smk_933, \
                         smk_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_17 * slk_930[k]
                    + f_3 * pc_x[k] * smk_930[k];

        t_1156[k] = f_17 * slk_931[k]
                    + f_3 * pc_x[k] * smk_931[k];

        t_1157[k] = f_17 * slk_932[k]
                    + f_3 * pc_x[k] * smk_932[k];

        t_1158[k] = f_17 * slk_933[k]
                    + f_3 * pc_x[k] * smk_933[k];

        t_1159[k] = f_17 * slk_934[k]
                    + f_3 * pc_x[k] * smk_934[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, pc_y, pc_z, slk_676, slk_712, slk_935, \
                         smi0_721, smi1_721, smk_928, smk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_17 * slk_935[k]
                    + f_3 * pc_x[k] * smk_935[k];

        t_1161[k] = f_16 * slk_712[k]
                    + f_1 * smi0_721[k]
                    - f_2 * smi1_721[k]
                    + f_3 * pc_y[k] * smk_928[k];

        t_1162[k] = f_18 * slk_676[k]
                    + f_3 * pc_z[k] * smk_928[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_y, slk_714, slk_715, slk_716, smi0_723, \
                         smi0_724, smi0_725, smi1_723, smi1_724, smi1_725, smk_930, smk_931, \
                         smk_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * slk_714[k]
                    + f_4 * smi0_723[k]
                    - f_5 * smi1_723[k]
                    + f_3 * pc_y[k] * smk_930[k];

        t_1164[k] = f_16 * slk_715[k]
                    + f_6 * smi0_724[k]
                    - f_7 * smi1_724[k]
                    + f_3 * pc_y[k] * smk_931[k];

        t_1165[k] = f_16 * slk_716[k]
                    + f_8 * smi0_725[k]
                    - f_9 * smi1_725[k]
                    + f_3 * pc_y[k] * smk_932[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_y, slk_717, slk_718, slk_719, smi0_726, \
                         smi0_727, smi1_726, smi1_727, smk_933, smk_934, \
                         smk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_16 * slk_717[k]
                    + f_10 * smi0_726[k]
                    - f_11 * smi1_726[k]
                    + f_3 * pc_y[k] * smk_933[k];

        t_1167[k] = f_16 * slk_718[k]
                    + f_12 * smi0_727[k]
                    - f_13 * smi1_727[k]
                    + f_3 * pc_y[k] * smk_934[k];

        t_1168[k] = f_16 * slk_719[k]
                    + f_3 * pc_y[k] * smk_935[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, pb_y, pc_y, pc_z, sll0_900, slk_683, \
                         slk_684, slk_720, sll1_900, smi0_727, smi1_727, smk_935, \
                         smk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_18 * slk_683[k]
                    + f_1 * smi0_727[k]
                    - f_2 * smi1_727[k]
                    + f_3 * pc_z[k] * smk_935[k];

        t_1170[k] = pb_y[k] * sll0_900[k]
                    - f_14 * pc_y[k] * sll1_900[k];

        t_1171[k] = f_15 * slk_720[k]
                    + f_3 * pc_y[k] * smk_936[k];

        t_1172[k] = f_19 * slk_684[k]
                    + f_3 * pc_z[k] * smk_936[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, t_1176, pb_y, pc_y, sll0_903, sll0_905, \
                         sll0_906, slk_721, slk_722, slk_723, sll1_903, sll1_905, sll1_906, \
                         smk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = pb_y[k] * sll0_903[k]
                    + f_16 * slk_721[k]
                    - f_14 * pc_y[k] * sll1_903[k];

        t_1174[k] = f_15 * slk_722[k]
                    + f_3 * pc_y[k] * smk_938[k];

        t_1175[k] = pb_y[k] * sll0_905[k]
                    - f_14 * pc_y[k] * sll1_905[k];

        t_1176[k] = pb_y[k] * sll0_906[k]
                    + f_17 * slk_723[k]
                    - f_14 * pc_y[k] * sll1_906[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, pb_y, pc_y, pc_z, sll0_909, sll0_910, \
                         slk_687, slk_725, slk_726, sll1_909, sll1_910, smk_939, \
                         smk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_19 * slk_687[k]
                    + f_3 * pc_z[k] * smk_939[k];

        t_1178[k] = f_15 * slk_725[k]
                    + f_3 * pc_y[k] * smk_941[k];

        t_1179[k] = pb_y[k] * sll0_909[k]
                    - f_14 * pc_y[k] * sll1_909[k];

        t_1180[k] = pb_y[k] * sll0_910[k]
                    + f_18 * slk_726[k]
                    - f_14 * pc_y[k] * sll1_910[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, t_1184, pb_y, pc_y, pc_z, sll0_912, sll0_914, \
                         slk_690, slk_728, slk_729, sll1_912, sll1_914, smk_942, \
                         smk_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_19 * slk_690[k]
                    + f_3 * pc_z[k] * smk_942[k];

        t_1182[k] = pb_y[k] * sll0_912[k]
                    + f_16 * slk_728[k]
                    - f_14 * pc_y[k] * sll1_912[k];

        t_1183[k] = f_15 * slk_729[k]
                    + f_3 * pc_y[k] * smk_945[k];

        t_1184[k] = pb_y[k] * sll0_914[k]
                    - f_14 * pc_y[k] * sll1_914[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pb_y, pc_y, pc_z, sll0_915, sll0_917, \
                         slk_694, slk_730, slk_732, sll1_915, sll1_917, \
                         smk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = pb_y[k] * sll0_915[k]
                    + f_19 * slk_730[k]
                    - f_14 * pc_y[k] * sll1_915[k];

        t_1186[k] = f_19 * slk_694[k]
                    + f_3 * pc_z[k] * smk_946[k];

        t_1187[k] = pb_y[k] * sll0_917[k]
                    + f_17 * slk_732[k]
                    - f_14 * pc_y[k] * sll1_917[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, t_1191, pb_y, pc_y, sll0_918, sll0_920, \
                         sll0_921, slk_733, slk_734, slk_735, sll1_918, sll1_920, sll1_921, \
                         smk_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = pb_y[k] * sll0_918[k]
                    + f_16 * slk_733[k]
                    - f_14 * pc_y[k] * sll1_918[k];

        t_1189[k] = f_15 * slk_734[k]
                    + f_3 * pc_y[k] * smk_950[k];

        t_1190[k] = pb_y[k] * sll0_920[k]
                    - f_14 * pc_y[k] * sll1_920[k];

        t_1191[k] = pb_y[k] * sll0_921[k]
                    + f_20 * slk_735[k]
                    - f_14 * pc_y[k] * sll1_921[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, pb_y, pc_y, pc_z, sll0_923, sll0_924, \
                         slk_699, slk_737, slk_738, sll1_923, sll1_924, \
                         smk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = f_19 * slk_699[k]
                    + f_3 * pc_z[k] * smk_951[k];

        t_1193[k] = pb_y[k] * sll0_923[k]
                    + f_18 * slk_737[k]
                    - f_14 * pc_y[k] * sll1_923[k];

        t_1194[k] = pb_y[k] * sll0_924[k]
                    + f_17 * slk_738[k]
                    - f_14 * pc_y[k] * sll1_924[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, pb_y, pc_x, pc_y, sll0_925, sll0_927, \
                         slk_739, slk_740, slk_964, sll1_925, sll1_927, smk_956, \
                         smk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = pb_y[k] * sll0_925[k]
                    + f_16 * slk_739[k]
                    - f_14 * pc_y[k] * sll1_925[k];

        t_1196[k] = f_15 * slk_740[k]
                    + f_3 * pc_y[k] * smk_956[k];

        t_1197[k] = pb_y[k] * sll0_927[k]
                    - f_14 * pc_y[k] * sll1_927[k];

        t_1198[k] = f_17 * slk_964[k]
                    + f_3 * pc_x[k] * smk_964[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, t_1202, t_1203, pc_x, slk_965, slk_966, \
                         slk_967, slk_968, slk_969, smk_965, smk_966, smk_967, smk_968, \
                         smk_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_17 * slk_965[k]
                    + f_3 * pc_x[k] * smk_965[k];

        t_1200[k] = f_17 * slk_966[k]
                    + f_3 * pc_x[k] * smk_966[k];

        t_1201[k] = f_17 * slk_967[k]
                    + f_3 * pc_x[k] * smk_967[k];

        t_1202[k] = f_17 * slk_968[k]
                    + f_3 * pc_x[k] * smk_968[k];

        t_1203[k] = f_17 * slk_969[k]
                    + f_3 * pc_x[k] * smk_969[k];
    }

#pragma omp simd aligned(t_1204, t_1205, t_1206, t_1207, pc_x, pc_y, pc_z, slk_712, slk_748, \
                         slk_970, slk_971, smi0_749, smi1_749, smk_964, smk_970, \
                         smk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1204[k] = f_17 * slk_970[k]
                    + f_3 * pc_x[k] * smk_970[k];

        t_1205[k] = f_17 * slk_971[k]
                    + f_3 * pc_x[k] * smk_971[k];

        t_1206[k] = f_15 * slk_748[k]
                    + f_1 * smi0_749[k]
                    - f_2 * smi1_749[k]
                    + f_3 * pc_y[k] * smk_964[k];

        t_1207[k] = f_19 * slk_712[k]
                    + f_3 * pc_z[k] * smk_964[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, pc_y, slk_750, slk_751, slk_752, smi0_751, \
                         smi0_752, smi0_753, smi1_751, smi1_752, smi1_753, smk_966, smk_967, \
                         smk_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_15 * slk_750[k]
                    + f_4 * smi0_751[k]
                    - f_5 * smi1_751[k]
                    + f_3 * pc_y[k] * smk_966[k];

        t_1209[k] = f_15 * slk_751[k]
                    + f_6 * smi0_752[k]
                    - f_7 * smi1_752[k]
                    + f_3 * pc_y[k] * smk_967[k];

        t_1210[k] = f_15 * slk_752[k]
                    + f_8 * smi0_753[k]
                    - f_9 * smi1_753[k]
                    + f_3 * pc_y[k] * smk_968[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smi0, const size_t smi1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_944 = buffer.data(sll0 + 944);
    const auto *sll0_945 = buffer.data(sll0 + 945);
    const auto *sll0_948 = buffer.data(sll0 + 948);
    const auto *sll0_951 = buffer.data(sll0 + 951);
    const auto *sll0_955 = buffer.data(sll0 + 955);
    const auto *sll0_957 = buffer.data(sll0 + 957);
    const auto *sll0_960 = buffer.data(sll0 + 960);

    const auto *slk_720 = buffer.data(slk + 720);
    const auto *slk_723 = buffer.data(slk + 723);
    const auto *slk_726 = buffer.data(slk + 726);
    const auto *slk_730 = buffer.data(slk + 730);
    const auto *slk_735 = buffer.data(slk + 735);
    const auto *slk_748 = buffer.data(slk + 748);
    const auto *slk_753 = buffer.data(slk + 753);
    const auto *slk_754 = buffer.data(slk + 754);
    const auto *slk_755 = buffer.data(slk + 755);
    const auto *slk_756 = buffer.data(slk + 756);
    const auto *slk_758 = buffer.data(slk + 758);
    const auto *slk_759 = buffer.data(slk + 759);
    const auto *slk_761 = buffer.data(slk + 761);
    const auto *slk_762 = buffer.data(slk + 762);
    const auto *slk_763 = buffer.data(slk + 763);
    const auto *slk_765 = buffer.data(slk + 765);
    const auto *slk_766 = buffer.data(slk + 766);
    const auto *slk_770 = buffer.data(slk + 770);
    const auto *slk_776 = buffer.data(slk + 776);
    const auto *slk_784 = buffer.data(slk + 784);
    const auto *slk_786 = buffer.data(slk + 786);
    const auto *slk_787 = buffer.data(slk + 787);
    const auto *slk_788 = buffer.data(slk + 788);
    const auto *slk_789 = buffer.data(slk + 789);
    const auto *slk_790 = buffer.data(slk + 790);
    const auto *slk_791 = buffer.data(slk + 791);
    const auto *slk_792 = buffer.data(slk + 792);
    const auto *slk_794 = buffer.data(slk + 794);
    const auto *slk_797 = buffer.data(slk + 797);
    const auto *slk_801 = buffer.data(slk + 801);
    const auto *slk_972 = buffer.data(slk + 972);
    const auto *slk_975 = buffer.data(slk + 975);
    const auto *slk_977 = buffer.data(slk + 977);
    const auto *slk_978 = buffer.data(slk + 978);
    const auto *slk_981 = buffer.data(slk + 981);
    const auto *slk_982 = buffer.data(slk + 982);
    const auto *slk_984 = buffer.data(slk + 984);
    const auto *slk_986 = buffer.data(slk + 986);
    const auto *slk_987 = buffer.data(slk + 987);
    const auto *slk_989 = buffer.data(slk + 989);
    const auto *slk_990 = buffer.data(slk + 990);
    const auto *slk_992 = buffer.data(slk + 992);
    const auto *slk_993 = buffer.data(slk + 993);
    const auto *slk_995 = buffer.data(slk + 995);
    const auto *slk_996 = buffer.data(slk + 996);
    const auto *slk_997 = buffer.data(slk + 997);
    const auto *slk_999 = buffer.data(slk + 999);
    const auto *slk_1000 = buffer.data(slk + 1000);
    const auto *slk_1001 = buffer.data(slk + 1001);
    const auto *slk_1002 = buffer.data(slk + 1002);
    const auto *slk_1003 = buffer.data(slk + 1003);
    const auto *slk_1004 = buffer.data(slk + 1004);
    const auto *slk_1005 = buffer.data(slk + 1005);
    const auto *slk_1006 = buffer.data(slk + 1006);
    const auto *slk_1007 = buffer.data(slk + 1007);
    const auto *slk_1008 = buffer.data(slk + 1008);
    const auto *slk_1011 = buffer.data(slk + 1011);
    const auto *slk_1013 = buffer.data(slk + 1013);
    const auto *slk_1014 = buffer.data(slk + 1014);
    const auto *slk_1017 = buffer.data(slk + 1017);
    const auto *slk_1018 = buffer.data(slk + 1018);
    const auto *slk_1020 = buffer.data(slk + 1020);
    const auto *slk_1022 = buffer.data(slk + 1022);
    const auto *slk_1023 = buffer.data(slk + 1023);
    const auto *slk_1025 = buffer.data(slk + 1025);
    const auto *slk_1026 = buffer.data(slk + 1026);
    const auto *slk_1028 = buffer.data(slk + 1028);
    const auto *slk_1029 = buffer.data(slk + 1029);
    const auto *slk_1031 = buffer.data(slk + 1031);
    const auto *slk_1032 = buffer.data(slk + 1032);
    const auto *slk_1033 = buffer.data(slk + 1033);
    const auto *slk_1035 = buffer.data(slk + 1035);
    const auto *slk_1036 = buffer.data(slk + 1036);
    const auto *slk_1037 = buffer.data(slk + 1037);
    const auto *slk_1038 = buffer.data(slk + 1038);
    const auto *slk_1039 = buffer.data(slk + 1039);
    const auto *slk_1040 = buffer.data(slk + 1040);
    const auto *slk_1041 = buffer.data(slk + 1041);
    const auto *slk_1042 = buffer.data(slk + 1042);
    const auto *slk_1043 = buffer.data(slk + 1043);
    const auto *slk_1049 = buffer.data(slk + 1049);
    const auto *slk_1053 = buffer.data(slk + 1053);
    const auto *slk_1058 = buffer.data(slk + 1058);

    const auto *sll1_944 = buffer.data(sll1 + 944);
    const auto *sll1_945 = buffer.data(sll1 + 945);
    const auto *sll1_948 = buffer.data(sll1 + 948);
    const auto *sll1_951 = buffer.data(sll1 + 951);
    const auto *sll1_955 = buffer.data(sll1 + 955);
    const auto *sll1_957 = buffer.data(sll1 + 957);
    const auto *sll1_960 = buffer.data(sll1 + 960);

    const auto *smi0_754 = buffer.data(smi0 + 754);
    const auto *smi0_755 = buffer.data(smi0 + 755);
    const auto *smi0_756 = buffer.data(smi0 + 756);
    const auto *smi0_759 = buffer.data(smi0 + 759);
    const auto *smi0_761 = buffer.data(smi0 + 761);
    const auto *smi0_762 = buffer.data(smi0 + 762);
    const auto *smi0_765 = buffer.data(smi0 + 765);
    const auto *smi0_766 = buffer.data(smi0 + 766);
    const auto *smi0_768 = buffer.data(smi0 + 768);
    const auto *smi0_770 = buffer.data(smi0 + 770);
    const auto *smi0_771 = buffer.data(smi0 + 771);
    const auto *smi0_773 = buffer.data(smi0 + 773);
    const auto *smi0_774 = buffer.data(smi0 + 774);
    const auto *smi0_776 = buffer.data(smi0 + 776);
    const auto *smi0_777 = buffer.data(smi0 + 777);
    const auto *smi0_779 = buffer.data(smi0 + 779);
    const auto *smi0_780 = buffer.data(smi0 + 780);
    const auto *smi0_781 = buffer.data(smi0 + 781);
    const auto *smi0_782 = buffer.data(smi0 + 782);
    const auto *smi0_783 = buffer.data(smi0 + 783);
    const auto *smi0_784 = buffer.data(smi0 + 784);
    const auto *smi0_787 = buffer.data(smi0 + 787);
    const auto *smi0_789 = buffer.data(smi0 + 789);
    const auto *smi0_790 = buffer.data(smi0 + 790);
    const auto *smi0_793 = buffer.data(smi0 + 793);
    const auto *smi0_794 = buffer.data(smi0 + 794);
    const auto *smi0_796 = buffer.data(smi0 + 796);
    const auto *smi0_798 = buffer.data(smi0 + 798);
    const auto *smi0_799 = buffer.data(smi0 + 799);
    const auto *smi0_801 = buffer.data(smi0 + 801);
    const auto *smi0_802 = buffer.data(smi0 + 802);
    const auto *smi0_804 = buffer.data(smi0 + 804);
    const auto *smi0_805 = buffer.data(smi0 + 805);
    const auto *smi0_807 = buffer.data(smi0 + 807);
    const auto *smi0_808 = buffer.data(smi0 + 808);
    const auto *smi0_809 = buffer.data(smi0 + 809);
    const auto *smi0_810 = buffer.data(smi0 + 810);
    const auto *smi0_811 = buffer.data(smi0 + 811);
    const auto *smi0_817 = buffer.data(smi0 + 817);
    const auto *smi0_821 = buffer.data(smi0 + 821);
    const auto *smi0_826 = buffer.data(smi0 + 826);

    const auto *smi1_754 = buffer.data(smi1 + 754);
    const auto *smi1_755 = buffer.data(smi1 + 755);
    const auto *smi1_756 = buffer.data(smi1 + 756);
    const auto *smi1_759 = buffer.data(smi1 + 759);
    const auto *smi1_761 = buffer.data(smi1 + 761);
    const auto *smi1_762 = buffer.data(smi1 + 762);
    const auto *smi1_765 = buffer.data(smi1 + 765);
    const auto *smi1_766 = buffer.data(smi1 + 766);
    const auto *smi1_768 = buffer.data(smi1 + 768);
    const auto *smi1_770 = buffer.data(smi1 + 770);
    const auto *smi1_771 = buffer.data(smi1 + 771);
    const auto *smi1_773 = buffer.data(smi1 + 773);
    const auto *smi1_774 = buffer.data(smi1 + 774);
    const auto *smi1_776 = buffer.data(smi1 + 776);
    const auto *smi1_777 = buffer.data(smi1 + 777);
    const auto *smi1_779 = buffer.data(smi1 + 779);
    const auto *smi1_780 = buffer.data(smi1 + 780);
    const auto *smi1_781 = buffer.data(smi1 + 781);
    const auto *smi1_782 = buffer.data(smi1 + 782);
    const auto *smi1_783 = buffer.data(smi1 + 783);
    const auto *smi1_784 = buffer.data(smi1 + 784);
    const auto *smi1_787 = buffer.data(smi1 + 787);
    const auto *smi1_789 = buffer.data(smi1 + 789);
    const auto *smi1_790 = buffer.data(smi1 + 790);
    const auto *smi1_793 = buffer.data(smi1 + 793);
    const auto *smi1_794 = buffer.data(smi1 + 794);
    const auto *smi1_796 = buffer.data(smi1 + 796);
    const auto *smi1_798 = buffer.data(smi1 + 798);
    const auto *smi1_799 = buffer.data(smi1 + 799);
    const auto *smi1_801 = buffer.data(smi1 + 801);
    const auto *smi1_802 = buffer.data(smi1 + 802);
    const auto *smi1_804 = buffer.data(smi1 + 804);
    const auto *smi1_805 = buffer.data(smi1 + 805);
    const auto *smi1_807 = buffer.data(smi1 + 807);
    const auto *smi1_808 = buffer.data(smi1 + 808);
    const auto *smi1_809 = buffer.data(smi1 + 809);
    const auto *smi1_810 = buffer.data(smi1 + 810);
    const auto *smi1_811 = buffer.data(smi1 + 811);
    const auto *smi1_817 = buffer.data(smi1 + 817);
    const auto *smi1_821 = buffer.data(smi1 + 821);
    const auto *smi1_826 = buffer.data(smi1 + 826);

    const auto *smk_969 = buffer.data(smk + 969);
    const auto *smk_970 = buffer.data(smk + 970);
    const auto *smk_971 = buffer.data(smk + 971);
    const auto *smk_972 = buffer.data(smk + 972);
    const auto *smk_974 = buffer.data(smk + 974);
    const auto *smk_975 = buffer.data(smk + 975);
    const auto *smk_977 = buffer.data(smk + 977);
    const auto *smk_978 = buffer.data(smk + 978);
    const auto *smk_981 = buffer.data(smk + 981);
    const auto *smk_982 = buffer.data(smk + 982);
    const auto *smk_984 = buffer.data(smk + 984);
    const auto *smk_986 = buffer.data(smk + 986);
    const auto *smk_987 = buffer.data(smk + 987);
    const auto *smk_989 = buffer.data(smk + 989);
    const auto *smk_990 = buffer.data(smk + 990);
    const auto *smk_992 = buffer.data(smk + 992);
    const auto *smk_993 = buffer.data(smk + 993);
    const auto *smk_995 = buffer.data(smk + 995);
    const auto *smk_996 = buffer.data(smk + 996);
    const auto *smk_997 = buffer.data(smk + 997);
    const auto *smk_999 = buffer.data(smk + 999);
    const auto *smk_1000 = buffer.data(smk + 1000);
    const auto *smk_1001 = buffer.data(smk + 1001);
    const auto *smk_1002 = buffer.data(smk + 1002);
    const auto *smk_1003 = buffer.data(smk + 1003);
    const auto *smk_1004 = buffer.data(smk + 1004);
    const auto *smk_1005 = buffer.data(smk + 1005);
    const auto *smk_1006 = buffer.data(smk + 1006);
    const auto *smk_1007 = buffer.data(smk + 1007);
    const auto *smk_1008 = buffer.data(smk + 1008);
    const auto *smk_1010 = buffer.data(smk + 1010);
    const auto *smk_1011 = buffer.data(smk + 1011);
    const auto *smk_1013 = buffer.data(smk + 1013);
    const auto *smk_1014 = buffer.data(smk + 1014);
    const auto *smk_1017 = buffer.data(smk + 1017);
    const auto *smk_1018 = buffer.data(smk + 1018);
    const auto *smk_1020 = buffer.data(smk + 1020);
    const auto *smk_1022 = buffer.data(smk + 1022);
    const auto *smk_1023 = buffer.data(smk + 1023);
    const auto *smk_1025 = buffer.data(smk + 1025);
    const auto *smk_1026 = buffer.data(smk + 1026);
    const auto *smk_1028 = buffer.data(smk + 1028);
    const auto *smk_1029 = buffer.data(smk + 1029);
    const auto *smk_1031 = buffer.data(smk + 1031);
    const auto *smk_1032 = buffer.data(smk + 1032);
    const auto *smk_1033 = buffer.data(smk + 1033);
    const auto *smk_1035 = buffer.data(smk + 1035);
    const auto *smk_1036 = buffer.data(smk + 1036);
    const auto *smk_1037 = buffer.data(smk + 1037);
    const auto *smk_1038 = buffer.data(smk + 1038);
    const auto *smk_1039 = buffer.data(smk + 1039);
    const auto *smk_1040 = buffer.data(smk + 1040);
    const auto *smk_1041 = buffer.data(smk + 1041);
    const auto *smk_1042 = buffer.data(smk + 1042);
    const auto *smk_1043 = buffer.data(smk + 1043);
    const auto *smk_1044 = buffer.data(smk + 1044);
    const auto *smk_1046 = buffer.data(smk + 1046);
    const auto *smk_1047 = buffer.data(smk + 1047);
    const auto *smk_1049 = buffer.data(smk + 1049);
    const auto *smk_1050 = buffer.data(smk + 1050);
    const auto *smk_1053 = buffer.data(smk + 1053);
    const auto *smk_1054 = buffer.data(smk + 1054);
    const auto *smk_1058 = buffer.data(smk + 1058);

#pragma omp simd aligned(t_1211, t_1212, t_1213, pc_y, slk_753, slk_754, slk_755, smi0_754, \
                         smi0_755, smi1_754, smi1_755, smk_969, smk_970, \
                         smk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = f_15 * slk_753[k]
                    + f_10 * smi0_754[k]
                    - f_11 * smi1_754[k]
                    + f_3 * pc_y[k] * smk_969[k];

        t_1212[k] = f_15 * slk_754[k]
                    + f_12 * smi0_755[k]
                    - f_13 * smi1_755[k]
                    + f_3 * pc_y[k] * smk_970[k];

        t_1213[k] = f_15 * slk_755[k]
                    + f_3 * pc_y[k] * smk_971[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pb_y, pc_x, pc_y, pc_z, sll0_944, \
                         slk_720, slk_972, sll1_944, smi0_756, smi1_756, \
                         smk_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pb_y[k] * sll0_944[k]
                    - f_14 * pc_y[k] * sll1_944[k];

        t_1215[k] = f_17 * slk_972[k]
                    + f_1 * smi0_756[k]
                    - f_2 * smi1_756[k]
                    + f_3 * pc_x[k] * smk_972[k];

        t_1216[k] = f_3 * pc_y[k] * smk_972[k];

        t_1217[k] = f_20 * slk_720[k]
                    + f_3 * pc_z[k] * smk_972[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, pc_x, pc_y, slk_975, slk_977, smi0_759, \
                         smi0_761, smi1_759, smi1_761, smk_974, smk_975, \
                         smk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = f_17 * slk_975[k]
                    + f_4 * smi0_759[k]
                    - f_5 * smi1_759[k]
                    + f_3 * pc_x[k] * smk_975[k];

        t_1219[k] = f_3 * pc_y[k] * smk_974[k];

        t_1220[k] = f_17 * slk_977[k]
                    + f_4 * smi0_761[k]
                    - f_5 * smi1_761[k]
                    + f_3 * pc_x[k] * smk_977[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, pc_x, pc_y, pc_z, slk_723, slk_978, smi0_762, \
                         smi1_762, smk_975, smk_977, smk_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_17 * slk_978[k]
                    + f_6 * smi0_762[k]
                    - f_7 * smi1_762[k]
                    + f_3 * pc_x[k] * smk_978[k];

        t_1222[k] = f_20 * slk_723[k]
                    + f_3 * pc_z[k] * smk_975[k];

        t_1223[k] = f_3 * pc_y[k] * smk_977[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pc_x, pc_z, slk_726, slk_981, slk_982, \
                         smi0_765, smi0_766, smi1_765, smi1_766, smk_978, smk_981, \
                         smk_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = f_17 * slk_981[k]
                    + f_6 * smi0_765[k]
                    - f_7 * smi1_765[k]
                    + f_3 * pc_x[k] * smk_981[k];

        t_1225[k] = f_17 * slk_982[k]
                    + f_8 * smi0_766[k]
                    - f_9 * smi1_766[k]
                    + f_3 * pc_x[k] * smk_982[k];

        t_1226[k] = f_20 * slk_726[k]
                    + f_3 * pc_z[k] * smk_978[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pc_x, pc_y, slk_984, slk_986, smi0_768, \
                         smi0_770, smi1_768, smi1_770, smk_981, smk_984, \
                         smk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_17 * slk_984[k]
                    + f_8 * smi0_768[k]
                    - f_9 * smi1_768[k]
                    + f_3 * pc_x[k] * smk_984[k];

        t_1228[k] = f_3 * pc_y[k] * smk_981[k];

        t_1229[k] = f_17 * slk_986[k]
                    + f_8 * smi0_770[k]
                    - f_9 * smi1_770[k]
                    + f_3 * pc_x[k] * smk_986[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, pc_x, pc_z, slk_730, slk_987, slk_989, \
                         smi0_771, smi0_773, smi1_771, smi1_773, smk_982, smk_987, \
                         smk_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_17 * slk_987[k]
                    + f_10 * smi0_771[k]
                    - f_11 * smi1_771[k]
                    + f_3 * pc_x[k] * smk_987[k];

        t_1231[k] = f_20 * slk_730[k]
                    + f_3 * pc_z[k] * smk_982[k];

        t_1232[k] = f_17 * slk_989[k]
                    + f_10 * smi0_773[k]
                    - f_11 * smi1_773[k]
                    + f_3 * pc_x[k] * smk_989[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, pc_x, pc_y, slk_990, slk_992, smi0_774, \
                         smi0_776, smi1_774, smi1_776, smk_986, smk_990, \
                         smk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = f_17 * slk_990[k]
                    + f_10 * smi0_774[k]
                    - f_11 * smi1_774[k]
                    + f_3 * pc_x[k] * smk_990[k];

        t_1234[k] = f_3 * pc_y[k] * smk_986[k];

        t_1235[k] = f_17 * slk_992[k]
                    + f_10 * smi0_776[k]
                    - f_11 * smi1_776[k]
                    + f_3 * pc_x[k] * smk_992[k];
    }

#pragma omp simd aligned(t_1236, t_1237, t_1238, pc_x, pc_z, slk_735, slk_993, slk_995, \
                         smi0_777, smi0_779, smi1_777, smi1_779, smk_987, smk_993, \
                         smk_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1236[k] = f_17 * slk_993[k]
                    + f_12 * smi0_777[k]
                    - f_13 * smi1_777[k]
                    + f_3 * pc_x[k] * smk_993[k];

        t_1237[k] = f_20 * slk_735[k]
                    + f_3 * pc_z[k] * smk_987[k];

        t_1238[k] = f_17 * slk_995[k]
                    + f_12 * smi0_779[k]
                    - f_13 * smi1_779[k]
                    + f_3 * pc_x[k] * smk_995[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, pc_x, pc_y, slk_996, slk_997, smi0_780, \
                         smi0_781, smi1_780, smi1_781, smk_992, smk_996, \
                         smk_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_17 * slk_996[k]
                    + f_12 * smi0_780[k]
                    - f_13 * smi1_780[k]
                    + f_3 * pc_x[k] * smk_996[k];

        t_1240[k] = f_17 * slk_997[k]
                    + f_12 * smi0_781[k]
                    - f_13 * smi1_781[k]
                    + f_3 * pc_x[k] * smk_997[k];

        t_1241[k] = f_3 * pc_y[k] * smk_992[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, pc_x, slk_999, slk_1000, slk_1001, \
                         slk_1002, smi0_783, smi1_783, smk_999, smk_1000, smk_1001, \
                         smk_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = f_17 * slk_999[k]
                    + f_12 * smi0_783[k]
                    - f_13 * smi1_783[k]
                    + f_3 * pc_x[k] * smk_999[k];

        t_1243[k] = f_17 * slk_1000[k]
                    + f_3 * pc_x[k] * smk_1000[k];

        t_1244[k] = f_17 * slk_1001[k]
                    + f_3 * pc_x[k] * smk_1001[k];

        t_1245[k] = f_17 * slk_1002[k]
                    + f_3 * pc_x[k] * smk_1002[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, t_1250, pc_x, slk_1003, slk_1004, \
                         slk_1005, slk_1006, slk_1007, smk_1003, smk_1004, smk_1005, smk_1006, \
                         smk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = f_17 * slk_1003[k]
                    + f_3 * pc_x[k] * smk_1003[k];

        t_1247[k] = f_17 * slk_1004[k]
                    + f_3 * pc_x[k] * smk_1004[k];

        t_1248[k] = f_17 * slk_1005[k]
                    + f_3 * pc_x[k] * smk_1005[k];

        t_1249[k] = f_17 * slk_1006[k]
                    + f_3 * pc_x[k] * smk_1006[k];

        t_1250[k] = f_17 * slk_1007[k]
                    + f_3 * pc_x[k] * smk_1007[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, pc_y, pc_z, slk_748, smi0_777, \
                         smi0_779, smi0_780, smi1_777, smi1_779, smi1_780, smk_1000, smk_1002, \
                         smk_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_1 * smi0_777[k]
                    - f_2 * smi1_777[k]
                    + f_3 * pc_y[k] * smk_1000[k];

        t_1252[k] = f_20 * slk_748[k]
                    + f_3 * pc_z[k] * smk_1000[k];

        t_1253[k] = f_4 * smi0_779[k]
                    - f_5 * smi1_779[k]
                    + f_3 * pc_y[k] * smk_1002[k];

        t_1254[k] = f_6 * smi0_780[k]
                    - f_7 * smi1_780[k]
                    + f_3 * pc_y[k] * smk_1003[k];
    }

#pragma omp simd aligned(t_1255, t_1256, t_1257, t_1258, pc_y, smi0_781, smi0_782, smi0_783, \
                         smi1_781, smi1_782, smi1_783, smk_1004, smk_1005, smk_1006, \
                         smk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1255[k] = f_8 * smi0_781[k]
                    - f_9 * smi1_781[k]
                    + f_3 * pc_y[k] * smk_1004[k];

        t_1256[k] = f_10 * smi0_782[k]
                    - f_11 * smi1_782[k]
                    + f_3 * pc_y[k] * smk_1005[k];

        t_1257[k] = f_12 * smi0_783[k]
                    - f_13 * smi1_783[k]
                    + f_3 * pc_y[k] * smk_1006[k];

        t_1258[k] = f_3 * pc_y[k] * smk_1007[k];
    }

#pragma omp simd aligned(t_1259, t_1260, t_1261, t_1262, pc_x, pc_y, pc_z, slk_755, slk_756, \
                         slk_1008, smi0_783, smi0_784, smi1_783, smi1_784, smk_1007, \
                         smk_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1259[k] = f_20 * slk_755[k]
                    + f_1 * smi0_783[k]
                    - f_2 * smi1_783[k]
                    + f_3 * pc_z[k] * smk_1007[k];

        t_1260[k] = f_16 * slk_1008[k]
                    + f_1 * smi0_784[k]
                    - f_2 * smi1_784[k]
                    + f_3 * pc_x[k] * smk_1008[k];

        t_1261[k] = f_22 * slk_756[k]
                    + f_3 * pc_y[k] * smk_1008[k];

        t_1262[k] = f_3 * pc_z[k] * smk_1008[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, pc_x, pc_y, slk_758, slk_1011, slk_1013, \
                         smi0_787, smi0_789, smi1_787, smi1_789, smk_1010, smk_1011, \
                         smk_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = f_16 * slk_1011[k]
                    + f_4 * smi0_787[k]
                    - f_5 * smi1_787[k]
                    + f_3 * pc_x[k] * smk_1011[k];

        t_1264[k] = f_22 * slk_758[k]
                    + f_3 * pc_y[k] * smk_1010[k];

        t_1265[k] = f_16 * slk_1013[k]
                    + f_4 * smi0_789[k]
                    - f_5 * smi1_789[k]
                    + f_3 * pc_x[k] * smk_1013[k];
    }

#pragma omp simd aligned(t_1266, t_1267, t_1268, pc_x, pc_y, pc_z, slk_761, slk_1014, \
                         smi0_790, smi1_790, smk_1011, smk_1013, \
                         smk_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1266[k] = f_16 * slk_1014[k]
                    + f_6 * smi0_790[k]
                    - f_7 * smi1_790[k]
                    + f_3 * pc_x[k] * smk_1014[k];

        t_1267[k] = f_3 * pc_z[k] * smk_1011[k];

        t_1268[k] = f_22 * slk_761[k]
                    + f_3 * pc_y[k] * smk_1013[k];
    }

#pragma omp simd aligned(t_1269, t_1270, t_1271, pc_x, pc_z, slk_1017, slk_1018, smi0_793, \
                         smi0_794, smi1_793, smi1_794, smk_1014, smk_1017, \
                         smk_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1269[k] = f_16 * slk_1017[k]
                    + f_6 * smi0_793[k]
                    - f_7 * smi1_793[k]
                    + f_3 * pc_x[k] * smk_1017[k];

        t_1270[k] = f_16 * slk_1018[k]
                    + f_8 * smi0_794[k]
                    - f_9 * smi1_794[k]
                    + f_3 * pc_x[k] * smk_1018[k];

        t_1271[k] = f_3 * pc_z[k] * smk_1014[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, pc_x, pc_y, slk_765, slk_1020, slk_1022, \
                         smi0_796, smi0_798, smi1_796, smi1_798, smk_1017, smk_1020, \
                         smk_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_16 * slk_1020[k]
                    + f_8 * smi0_796[k]
                    - f_9 * smi1_796[k]
                    + f_3 * pc_x[k] * smk_1020[k];

        t_1273[k] = f_22 * slk_765[k]
                    + f_3 * pc_y[k] * smk_1017[k];

        t_1274[k] = f_16 * slk_1022[k]
                    + f_8 * smi0_798[k]
                    - f_9 * smi1_798[k]
                    + f_3 * pc_x[k] * smk_1022[k];
    }

#pragma omp simd aligned(t_1275, t_1276, t_1277, pc_x, pc_z, slk_1023, slk_1025, smi0_799, \
                         smi0_801, smi1_799, smi1_801, smk_1018, smk_1023, \
                         smk_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1275[k] = f_16 * slk_1023[k]
                    + f_10 * smi0_799[k]
                    - f_11 * smi1_799[k]
                    + f_3 * pc_x[k] * smk_1023[k];

        t_1276[k] = f_3 * pc_z[k] * smk_1018[k];

        t_1277[k] = f_16 * slk_1025[k]
                    + f_10 * smi0_801[k]
                    - f_11 * smi1_801[k]
                    + f_3 * pc_x[k] * smk_1025[k];
    }

#pragma omp simd aligned(t_1278, t_1279, t_1280, pc_x, pc_y, slk_770, slk_1026, slk_1028, \
                         smi0_802, smi0_804, smi1_802, smi1_804, smk_1022, smk_1026, \
                         smk_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1278[k] = f_16 * slk_1026[k]
                    + f_10 * smi0_802[k]
                    - f_11 * smi1_802[k]
                    + f_3 * pc_x[k] * smk_1026[k];

        t_1279[k] = f_22 * slk_770[k]
                    + f_3 * pc_y[k] * smk_1022[k];

        t_1280[k] = f_16 * slk_1028[k]
                    + f_10 * smi0_804[k]
                    - f_11 * smi1_804[k]
                    + f_3 * pc_x[k] * smk_1028[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, pc_x, pc_z, slk_1029, slk_1031, smi0_805, \
                         smi0_807, smi1_805, smi1_807, smk_1023, smk_1029, \
                         smk_1031 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_16 * slk_1029[k]
                    + f_12 * smi0_805[k]
                    - f_13 * smi1_805[k]
                    + f_3 * pc_x[k] * smk_1029[k];

        t_1282[k] = f_3 * pc_z[k] * smk_1023[k];

        t_1283[k] = f_16 * slk_1031[k]
                    + f_12 * smi0_807[k]
                    - f_13 * smi1_807[k]
                    + f_3 * pc_x[k] * smk_1031[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, pc_x, pc_y, slk_776, slk_1032, slk_1033, \
                         smi0_808, smi0_809, smi1_808, smi1_809, smk_1028, smk_1032, \
                         smk_1033 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_16 * slk_1032[k]
                    + f_12 * smi0_808[k]
                    - f_13 * smi1_808[k]
                    + f_3 * pc_x[k] * smk_1032[k];

        t_1285[k] = f_16 * slk_1033[k]
                    + f_12 * smi0_809[k]
                    - f_13 * smi1_809[k]
                    + f_3 * pc_x[k] * smk_1033[k];

        t_1286[k] = f_22 * slk_776[k]
                    + f_3 * pc_y[k] * smk_1028[k];
    }

#pragma omp simd aligned(t_1287, t_1288, t_1289, t_1290, pc_x, slk_1035, slk_1036, slk_1037, \
                         slk_1038, smi0_811, smi1_811, smk_1035, smk_1036, smk_1037, \
                         smk_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1287[k] = f_16 * slk_1035[k]
                    + f_12 * smi0_811[k]
                    - f_13 * smi1_811[k]
                    + f_3 * pc_x[k] * smk_1035[k];

        t_1288[k] = f_16 * slk_1036[k]
                    + f_3 * pc_x[k] * smk_1036[k];

        t_1289[k] = f_16 * slk_1037[k]
                    + f_3 * pc_x[k] * smk_1037[k];

        t_1290[k] = f_16 * slk_1038[k]
                    + f_3 * pc_x[k] * smk_1038[k];
    }

#pragma omp simd aligned(t_1291, t_1292, t_1293, t_1294, t_1295, pc_x, slk_1039, slk_1040, \
                         slk_1041, slk_1042, slk_1043, smk_1039, smk_1040, smk_1041, smk_1042, \
                         smk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1291[k] = f_16 * slk_1039[k]
                    + f_3 * pc_x[k] * smk_1039[k];

        t_1292[k] = f_16 * slk_1040[k]
                    + f_3 * pc_x[k] * smk_1040[k];

        t_1293[k] = f_16 * slk_1041[k]
                    + f_3 * pc_x[k] * smk_1041[k];

        t_1294[k] = f_16 * slk_1042[k]
                    + f_3 * pc_x[k] * smk_1042[k];

        t_1295[k] = f_16 * slk_1043[k]
                    + f_3 * pc_x[k] * smk_1043[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, pc_y, pc_z, slk_784, slk_786, smi0_805, \
                         smi0_807, smi1_805, smi1_807, smk_1036, \
                         smk_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_22 * slk_784[k]
                    + f_1 * smi0_805[k]
                    - f_2 * smi1_805[k]
                    + f_3 * pc_y[k] * smk_1036[k];

        t_1297[k] = f_3 * pc_z[k] * smk_1036[k];

        t_1298[k] = f_22 * slk_786[k]
                    + f_4 * smi0_807[k]
                    - f_5 * smi1_807[k]
                    + f_3 * pc_y[k] * smk_1038[k];
    }

#pragma omp simd aligned(t_1299, t_1300, t_1301, pc_y, slk_787, slk_788, slk_789, smi0_808, \
                         smi0_809, smi0_810, smi1_808, smi1_809, smi1_810, smk_1039, smk_1040, \
                         smk_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1299[k] = f_22 * slk_787[k]
                    + f_6 * smi0_808[k]
                    - f_7 * smi1_808[k]
                    + f_3 * pc_y[k] * smk_1039[k];

        t_1300[k] = f_22 * slk_788[k]
                    + f_8 * smi0_809[k]
                    - f_9 * smi1_809[k]
                    + f_3 * pc_y[k] * smk_1040[k];

        t_1301[k] = f_22 * slk_789[k]
                    + f_10 * smi0_810[k]
                    - f_11 * smi1_810[k]
                    + f_3 * pc_y[k] * smk_1041[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, t_1305, pb_z, pc_y, pc_z, sll0_945, slk_790, \
                         slk_791, sll1_945, smi0_811, smi1_811, smk_1042, \
                         smk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = f_22 * slk_790[k]
                    + f_12 * smi0_811[k]
                    - f_13 * smi1_811[k]
                    + f_3 * pc_y[k] * smk_1042[k];

        t_1303[k] = f_22 * slk_791[k]
                    + f_3 * pc_y[k] * smk_1043[k];

        t_1304[k] = f_1 * smi0_811[k]
                    - f_2 * smi1_811[k]
                    + f_3 * pc_z[k] * smk_1043[k];

        t_1305[k] = pb_z[k] * sll0_945[k]
                    - f_14 * pc_z[k] * sll1_945[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, t_1309, pb_z, pc_y, pc_z, sll0_948, slk_756, \
                         slk_792, slk_794, sll1_948, smk_1044, \
                         smk_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_20 * slk_792[k]
                    + f_3 * pc_y[k] * smk_1044[k];

        t_1307[k] = f_15 * slk_756[k]
                    + f_3 * pc_z[k] * smk_1044[k];

        t_1308[k] = pb_z[k] * sll0_948[k]
                    - f_14 * pc_z[k] * sll1_948[k];

        t_1309[k] = f_20 * slk_794[k]
                    + f_3 * pc_y[k] * smk_1046[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, pb_z, pc_x, pc_z, sll0_951, slk_759, \
                         slk_1049, sll1_951, smi0_817, smi1_817, smk_1047, \
                         smk_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_16 * slk_1049[k]
                    + f_4 * smi0_817[k]
                    - f_5 * smi1_817[k]
                    + f_3 * pc_x[k] * smk_1049[k];

        t_1311[k] = pb_z[k] * sll0_951[k]
                    - f_14 * pc_z[k] * sll1_951[k];

        t_1312[k] = f_15 * slk_759[k]
                    + f_3 * pc_z[k] * smk_1047[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, pb_z, pc_x, pc_y, pc_z, sll0_955, slk_797, \
                         slk_1053, sll1_955, smi0_821, smi1_821, smk_1049, \
                         smk_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_20 * slk_797[k]
                    + f_3 * pc_y[k] * smk_1049[k];

        t_1314[k] = f_16 * slk_1053[k]
                    + f_6 * smi0_821[k]
                    - f_7 * smi1_821[k]
                    + f_3 * pc_x[k] * smk_1053[k];

        t_1315[k] = pb_z[k] * sll0_955[k]
                    - f_14 * pc_z[k] * sll1_955[k];
    }

#pragma omp simd aligned(t_1316, t_1317, t_1318, pb_z, pc_y, pc_z, sll0_957, slk_762, slk_763, \
                         slk_801, sll1_957, smk_1050, smk_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1316[k] = f_15 * slk_762[k]
                    + f_3 * pc_z[k] * smk_1050[k];

        t_1317[k] = pb_z[k] * sll0_957[k]
                    + f_16 * slk_763[k]
                    - f_14 * pc_z[k] * sll1_957[k];

        t_1318[k] = f_20 * slk_801[k]
                    + f_3 * pc_y[k] * smk_1053[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, pb_z, pc_x, pc_z, sll0_960, slk_766, \
                         slk_1058, sll1_960, smi0_826, smi1_826, smk_1054, \
                         smk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_16 * slk_1058[k]
                    + f_8 * smi0_826[k]
                    - f_9 * smi1_826[k]
                    + f_3 * pc_x[k] * smk_1058[k];

        t_1320[k] = pb_z[k] * sll0_960[k]
                    - f_14 * pc_z[k] * sll1_960[k];

        t_1321[k] = f_15 * slk_766[k]
                    + f_3 * pc_z[k] * smk_1054[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smi0, const size_t smi1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_962 = buffer.data(sll0 + 962);
    const auto *sll0_963 = buffer.data(sll0 + 963);
    const auto *sll0_966 = buffer.data(sll0 + 966);
    const auto *sll0_968 = buffer.data(sll0 + 968);
    const auto *sll0_969 = buffer.data(sll0 + 969);
    const auto *sll0_970 = buffer.data(sll0 + 970);
    const auto *sll0_981 = buffer.data(sll0 + 981);

    const auto *slk_767 = buffer.data(slk + 767);
    const auto *slk_768 = buffer.data(slk + 768);
    const auto *slk_771 = buffer.data(slk + 771);
    const auto *slk_772 = buffer.data(slk + 772);
    const auto *slk_773 = buffer.data(slk + 773);
    const auto *slk_774 = buffer.data(slk + 774);
    const auto *slk_784 = buffer.data(slk + 784);
    const auto *slk_791 = buffer.data(slk + 791);
    const auto *slk_792 = buffer.data(slk + 792);
    const auto *slk_795 = buffer.data(slk + 795);
    const auto *slk_798 = buffer.data(slk + 798);
    const auto *slk_802 = buffer.data(slk + 802);
    const auto *slk_806 = buffer.data(slk + 806);
    const auto *slk_807 = buffer.data(slk + 807);
    const auto *slk_812 = buffer.data(slk + 812);
    const auto *slk_820 = buffer.data(slk + 820);
    const auto *slk_822 = buffer.data(slk + 822);
    const auto *slk_823 = buffer.data(slk + 823);
    const auto *slk_824 = buffer.data(slk + 824);
    const auto *slk_825 = buffer.data(slk + 825);
    const auto *slk_826 = buffer.data(slk + 826);
    const auto *slk_827 = buffer.data(slk + 827);
    const auto *slk_828 = buffer.data(slk + 828);
    const auto *slk_830 = buffer.data(slk + 830);
    const auto *slk_831 = buffer.data(slk + 831);
    const auto *slk_833 = buffer.data(slk + 833);
    const auto *slk_834 = buffer.data(slk + 834);
    const auto *slk_837 = buffer.data(slk + 837);
    const auto *slk_838 = buffer.data(slk + 838);
    const auto *slk_842 = buffer.data(slk + 842);
    const auto *slk_843 = buffer.data(slk + 843);
    const auto *slk_848 = buffer.data(slk + 848);
    const auto *slk_856 = buffer.data(slk + 856);
    const auto *slk_858 = buffer.data(slk + 858);
    const auto *slk_859 = buffer.data(slk + 859);
    const auto *slk_860 = buffer.data(slk + 860);
    const auto *slk_861 = buffer.data(slk + 861);
    const auto *slk_862 = buffer.data(slk + 862);
    const auto *slk_863 = buffer.data(slk + 863);
    const auto *slk_864 = buffer.data(slk + 864);
    const auto *slk_866 = buffer.data(slk + 866);
    const auto *slk_869 = buffer.data(slk + 869);
    const auto *slk_873 = buffer.data(slk + 873);
    const auto *slk_878 = buffer.data(slk + 878);
    const auto *slk_884 = buffer.data(slk + 884);
    const auto *slk_1064 = buffer.data(slk + 1064);
    const auto *slk_1071 = buffer.data(slk + 1071);
    const auto *slk_1072 = buffer.data(slk + 1072);
    const auto *slk_1073 = buffer.data(slk + 1073);
    const auto *slk_1074 = buffer.data(slk + 1074);
    const auto *slk_1075 = buffer.data(slk + 1075);
    const auto *slk_1076 = buffer.data(slk + 1076);
    const auto *slk_1077 = buffer.data(slk + 1077);
    const auto *slk_1078 = buffer.data(slk + 1078);
    const auto *slk_1079 = buffer.data(slk + 1079);
    const auto *slk_1080 = buffer.data(slk + 1080);
    const auto *slk_1083 = buffer.data(slk + 1083);
    const auto *slk_1085 = buffer.data(slk + 1085);
    const auto *slk_1086 = buffer.data(slk + 1086);
    const auto *slk_1089 = buffer.data(slk + 1089);
    const auto *slk_1090 = buffer.data(slk + 1090);
    const auto *slk_1092 = buffer.data(slk + 1092);
    const auto *slk_1094 = buffer.data(slk + 1094);
    const auto *slk_1095 = buffer.data(slk + 1095);
    const auto *slk_1097 = buffer.data(slk + 1097);
    const auto *slk_1098 = buffer.data(slk + 1098);
    const auto *slk_1100 = buffer.data(slk + 1100);
    const auto *slk_1101 = buffer.data(slk + 1101);
    const auto *slk_1103 = buffer.data(slk + 1103);
    const auto *slk_1104 = buffer.data(slk + 1104);
    const auto *slk_1105 = buffer.data(slk + 1105);
    const auto *slk_1107 = buffer.data(slk + 1107);
    const auto *slk_1108 = buffer.data(slk + 1108);
    const auto *slk_1109 = buffer.data(slk + 1109);
    const auto *slk_1110 = buffer.data(slk + 1110);
    const auto *slk_1111 = buffer.data(slk + 1111);
    const auto *slk_1112 = buffer.data(slk + 1112);
    const auto *slk_1113 = buffer.data(slk + 1113);
    const auto *slk_1114 = buffer.data(slk + 1114);
    const auto *slk_1115 = buffer.data(slk + 1115);
    const auto *slk_1116 = buffer.data(slk + 1116);
    const auto *slk_1119 = buffer.data(slk + 1119);
    const auto *slk_1121 = buffer.data(slk + 1121);
    const auto *slk_1122 = buffer.data(slk + 1122);
    const auto *slk_1125 = buffer.data(slk + 1125);
    const auto *slk_1126 = buffer.data(slk + 1126);
    const auto *slk_1128 = buffer.data(slk + 1128);
    const auto *slk_1130 = buffer.data(slk + 1130);
    const auto *slk_1131 = buffer.data(slk + 1131);
    const auto *slk_1133 = buffer.data(slk + 1133);
    const auto *slk_1134 = buffer.data(slk + 1134);
    const auto *slk_1136 = buffer.data(slk + 1136);
    const auto *slk_1137 = buffer.data(slk + 1137);
    const auto *slk_1139 = buffer.data(slk + 1139);
    const auto *slk_1140 = buffer.data(slk + 1140);
    const auto *slk_1141 = buffer.data(slk + 1141);
    const auto *slk_1143 = buffer.data(slk + 1143);
    const auto *slk_1144 = buffer.data(slk + 1144);
    const auto *slk_1145 = buffer.data(slk + 1145);

    const auto *sll1_962 = buffer.data(sll1 + 962);
    const auto *sll1_963 = buffer.data(sll1 + 963);
    const auto *sll1_966 = buffer.data(sll1 + 966);
    const auto *sll1_968 = buffer.data(sll1 + 968);
    const auto *sll1_969 = buffer.data(sll1 + 969);
    const auto *sll1_970 = buffer.data(sll1 + 970);
    const auto *sll1_981 = buffer.data(sll1 + 981);

    const auto *smi0_832 = buffer.data(smi0 + 832);
    const auto *smi0_835 = buffer.data(smi0 + 835);
    const auto *smi0_836 = buffer.data(smi0 + 836);
    const auto *smi0_837 = buffer.data(smi0 + 837);
    const auto *smi0_838 = buffer.data(smi0 + 838);
    const auto *smi0_839 = buffer.data(smi0 + 839);
    const auto *smi0_840 = buffer.data(smi0 + 840);
    const auto *smi0_843 = buffer.data(smi0 + 843);
    const auto *smi0_845 = buffer.data(smi0 + 845);
    const auto *smi0_846 = buffer.data(smi0 + 846);
    const auto *smi0_849 = buffer.data(smi0 + 849);
    const auto *smi0_850 = buffer.data(smi0 + 850);
    const auto *smi0_852 = buffer.data(smi0 + 852);
    const auto *smi0_854 = buffer.data(smi0 + 854);
    const auto *smi0_855 = buffer.data(smi0 + 855);
    const auto *smi0_857 = buffer.data(smi0 + 857);
    const auto *smi0_858 = buffer.data(smi0 + 858);
    const auto *smi0_860 = buffer.data(smi0 + 860);
    const auto *smi0_861 = buffer.data(smi0 + 861);
    const auto *smi0_863 = buffer.data(smi0 + 863);
    const auto *smi0_864 = buffer.data(smi0 + 864);
    const auto *smi0_865 = buffer.data(smi0 + 865);
    const auto *smi0_866 = buffer.data(smi0 + 866);
    const auto *smi0_867 = buffer.data(smi0 + 867);
    const auto *smi0_868 = buffer.data(smi0 + 868);
    const auto *smi0_871 = buffer.data(smi0 + 871);
    const auto *smi0_873 = buffer.data(smi0 + 873);
    const auto *smi0_874 = buffer.data(smi0 + 874);
    const auto *smi0_877 = buffer.data(smi0 + 877);
    const auto *smi0_878 = buffer.data(smi0 + 878);
    const auto *smi0_880 = buffer.data(smi0 + 880);
    const auto *smi0_882 = buffer.data(smi0 + 882);
    const auto *smi0_883 = buffer.data(smi0 + 883);
    const auto *smi0_885 = buffer.data(smi0 + 885);
    const auto *smi0_886 = buffer.data(smi0 + 886);
    const auto *smi0_888 = buffer.data(smi0 + 888);
    const auto *smi0_889 = buffer.data(smi0 + 889);
    const auto *smi0_891 = buffer.data(smi0 + 891);
    const auto *smi0_892 = buffer.data(smi0 + 892);
    const auto *smi0_893 = buffer.data(smi0 + 893);
    const auto *smi0_895 = buffer.data(smi0 + 895);

    const auto *smi1_832 = buffer.data(smi1 + 832);
    const auto *smi1_835 = buffer.data(smi1 + 835);
    const auto *smi1_836 = buffer.data(smi1 + 836);
    const auto *smi1_837 = buffer.data(smi1 + 837);
    const auto *smi1_838 = buffer.data(smi1 + 838);
    const auto *smi1_839 = buffer.data(smi1 + 839);
    const auto *smi1_840 = buffer.data(smi1 + 840);
    const auto *smi1_843 = buffer.data(smi1 + 843);
    const auto *smi1_845 = buffer.data(smi1 + 845);
    const auto *smi1_846 = buffer.data(smi1 + 846);
    const auto *smi1_849 = buffer.data(smi1 + 849);
    const auto *smi1_850 = buffer.data(smi1 + 850);
    const auto *smi1_852 = buffer.data(smi1 + 852);
    const auto *smi1_854 = buffer.data(smi1 + 854);
    const auto *smi1_855 = buffer.data(smi1 + 855);
    const auto *smi1_857 = buffer.data(smi1 + 857);
    const auto *smi1_858 = buffer.data(smi1 + 858);
    const auto *smi1_860 = buffer.data(smi1 + 860);
    const auto *smi1_861 = buffer.data(smi1 + 861);
    const auto *smi1_863 = buffer.data(smi1 + 863);
    const auto *smi1_864 = buffer.data(smi1 + 864);
    const auto *smi1_865 = buffer.data(smi1 + 865);
    const auto *smi1_866 = buffer.data(smi1 + 866);
    const auto *smi1_867 = buffer.data(smi1 + 867);
    const auto *smi1_868 = buffer.data(smi1 + 868);
    const auto *smi1_871 = buffer.data(smi1 + 871);
    const auto *smi1_873 = buffer.data(smi1 + 873);
    const auto *smi1_874 = buffer.data(smi1 + 874);
    const auto *smi1_877 = buffer.data(smi1 + 877);
    const auto *smi1_878 = buffer.data(smi1 + 878);
    const auto *smi1_880 = buffer.data(smi1 + 880);
    const auto *smi1_882 = buffer.data(smi1 + 882);
    const auto *smi1_883 = buffer.data(smi1 + 883);
    const auto *smi1_885 = buffer.data(smi1 + 885);
    const auto *smi1_886 = buffer.data(smi1 + 886);
    const auto *smi1_888 = buffer.data(smi1 + 888);
    const auto *smi1_889 = buffer.data(smi1 + 889);
    const auto *smi1_891 = buffer.data(smi1 + 891);
    const auto *smi1_892 = buffer.data(smi1 + 892);
    const auto *smi1_893 = buffer.data(smi1 + 893);
    const auto *smi1_895 = buffer.data(smi1 + 895);

    const auto *smk_1058 = buffer.data(smk + 1058);
    const auto *smk_1059 = buffer.data(smk + 1059);
    const auto *smk_1064 = buffer.data(smk + 1064);
    const auto *smk_1071 = buffer.data(smk + 1071);
    const auto *smk_1072 = buffer.data(smk + 1072);
    const auto *smk_1073 = buffer.data(smk + 1073);
    const auto *smk_1074 = buffer.data(smk + 1074);
    const auto *smk_1075 = buffer.data(smk + 1075);
    const auto *smk_1076 = buffer.data(smk + 1076);
    const auto *smk_1077 = buffer.data(smk + 1077);
    const auto *smk_1078 = buffer.data(smk + 1078);
    const auto *smk_1079 = buffer.data(smk + 1079);
    const auto *smk_1080 = buffer.data(smk + 1080);
    const auto *smk_1082 = buffer.data(smk + 1082);
    const auto *smk_1083 = buffer.data(smk + 1083);
    const auto *smk_1085 = buffer.data(smk + 1085);
    const auto *smk_1086 = buffer.data(smk + 1086);
    const auto *smk_1089 = buffer.data(smk + 1089);
    const auto *smk_1090 = buffer.data(smk + 1090);
    const auto *smk_1092 = buffer.data(smk + 1092);
    const auto *smk_1094 = buffer.data(smk + 1094);
    const auto *smk_1095 = buffer.data(smk + 1095);
    const auto *smk_1097 = buffer.data(smk + 1097);
    const auto *smk_1098 = buffer.data(smk + 1098);
    const auto *smk_1100 = buffer.data(smk + 1100);
    const auto *smk_1101 = buffer.data(smk + 1101);
    const auto *smk_1103 = buffer.data(smk + 1103);
    const auto *smk_1104 = buffer.data(smk + 1104);
    const auto *smk_1105 = buffer.data(smk + 1105);
    const auto *smk_1107 = buffer.data(smk + 1107);
    const auto *smk_1108 = buffer.data(smk + 1108);
    const auto *smk_1109 = buffer.data(smk + 1109);
    const auto *smk_1110 = buffer.data(smk + 1110);
    const auto *smk_1111 = buffer.data(smk + 1111);
    const auto *smk_1112 = buffer.data(smk + 1112);
    const auto *smk_1113 = buffer.data(smk + 1113);
    const auto *smk_1114 = buffer.data(smk + 1114);
    const auto *smk_1115 = buffer.data(smk + 1115);
    const auto *smk_1116 = buffer.data(smk + 1116);
    const auto *smk_1118 = buffer.data(smk + 1118);
    const auto *smk_1119 = buffer.data(smk + 1119);
    const auto *smk_1121 = buffer.data(smk + 1121);
    const auto *smk_1122 = buffer.data(smk + 1122);
    const auto *smk_1125 = buffer.data(smk + 1125);
    const auto *smk_1126 = buffer.data(smk + 1126);
    const auto *smk_1128 = buffer.data(smk + 1128);
    const auto *smk_1130 = buffer.data(smk + 1130);
    const auto *smk_1131 = buffer.data(smk + 1131);
    const auto *smk_1133 = buffer.data(smk + 1133);
    const auto *smk_1134 = buffer.data(smk + 1134);
    const auto *smk_1136 = buffer.data(smk + 1136);
    const auto *smk_1137 = buffer.data(smk + 1137);
    const auto *smk_1139 = buffer.data(smk + 1139);
    const auto *smk_1140 = buffer.data(smk + 1140);
    const auto *smk_1141 = buffer.data(smk + 1141);
    const auto *smk_1143 = buffer.data(smk + 1143);
    const auto *smk_1144 = buffer.data(smk + 1144);
    const auto *smk_1145 = buffer.data(smk + 1145);

#pragma omp simd aligned(t_1322, t_1323, t_1324, pb_z, pc_y, pc_z, sll0_962, sll0_963, \
                         slk_767, slk_768, slk_806, sll1_962, sll1_963, \
                         smk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = pb_z[k] * sll0_962[k]
                    + f_16 * slk_767[k]
                    - f_14 * pc_z[k] * sll1_962[k];

        t_1323[k] = pb_z[k] * sll0_963[k]
                    + f_17 * slk_768[k]
                    - f_14 * pc_z[k] * sll1_963[k];

        t_1324[k] = f_20 * slk_806[k]
                    + f_3 * pc_y[k] * smk_1058[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pb_z, pc_x, pc_z, sll0_966, slk_771, \
                         slk_1064, sll1_966, smi0_832, smi1_832, smk_1059, \
                         smk_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_16 * slk_1064[k]
                    + f_10 * smi0_832[k]
                    - f_11 * smi1_832[k]
                    + f_3 * pc_x[k] * smk_1064[k];

        t_1326[k] = pb_z[k] * sll0_966[k]
                    - f_14 * pc_z[k] * sll1_966[k];

        t_1327[k] = f_15 * slk_771[k]
                    + f_3 * pc_z[k] * smk_1059[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pb_z, pc_z, sll0_968, sll0_969, sll0_970, \
                         slk_772, slk_773, slk_774, sll1_968, sll1_969, \
                         sll1_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = pb_z[k] * sll0_968[k]
                    + f_16 * slk_772[k]
                    - f_14 * pc_z[k] * sll1_968[k];

        t_1329[k] = pb_z[k] * sll0_969[k]
                    + f_17 * slk_773[k]
                    - f_14 * pc_z[k] * sll1_969[k];

        t_1330[k] = pb_z[k] * sll0_970[k]
                    + f_18 * slk_774[k]
                    - f_14 * pc_z[k] * sll1_970[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, pc_x, pc_y, slk_812, slk_1071, \
                         slk_1072, slk_1073, smi0_839, smi1_839, smk_1064, smk_1071, smk_1072, \
                         smk_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_20 * slk_812[k]
                    + f_3 * pc_y[k] * smk_1064[k];

        t_1332[k] = f_16 * slk_1071[k]
                    + f_12 * smi0_839[k]
                    - f_13 * smi1_839[k]
                    + f_3 * pc_x[k] * smk_1071[k];

        t_1333[k] = f_16 * slk_1072[k]
                    + f_3 * pc_x[k] * smk_1072[k];

        t_1334[k] = f_16 * slk_1073[k]
                    + f_3 * pc_x[k] * smk_1073[k];
    }

#pragma omp simd aligned(t_1335, t_1336, t_1337, t_1338, t_1339, pc_x, slk_1074, slk_1075, \
                         slk_1076, slk_1077, slk_1078, smk_1074, smk_1075, smk_1076, smk_1077, \
                         smk_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1335[k] = f_16 * slk_1074[k]
                    + f_3 * pc_x[k] * smk_1074[k];

        t_1336[k] = f_16 * slk_1075[k]
                    + f_3 * pc_x[k] * smk_1075[k];

        t_1337[k] = f_16 * slk_1076[k]
                    + f_3 * pc_x[k] * smk_1076[k];

        t_1338[k] = f_16 * slk_1077[k]
                    + f_3 * pc_x[k] * smk_1077[k];

        t_1339[k] = f_16 * slk_1078[k]
                    + f_3 * pc_x[k] * smk_1078[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, pb_z, pc_x, pc_z, sll0_981, slk_784, \
                         slk_1079, sll1_981, smk_1072, smk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_16 * slk_1079[k]
                    + f_3 * pc_x[k] * smk_1079[k];

        t_1341[k] = pb_z[k] * sll0_981[k]
                    - f_14 * pc_z[k] * sll1_981[k];

        t_1342[k] = f_15 * slk_784[k]
                    + f_3 * pc_z[k] * smk_1072[k];
    }

#pragma omp simd aligned(t_1343, t_1344, t_1345, pc_y, slk_822, slk_823, slk_824, smi0_835, \
                         smi0_836, smi0_837, smi1_835, smi1_836, smi1_837, smk_1074, smk_1075, \
                         smk_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1343[k] = f_20 * slk_822[k]
                    + f_4 * smi0_835[k]
                    - f_5 * smi1_835[k]
                    + f_3 * pc_y[k] * smk_1074[k];

        t_1344[k] = f_20 * slk_823[k]
                    + f_6 * smi0_836[k]
                    - f_7 * smi1_836[k]
                    + f_3 * pc_y[k] * smk_1075[k];

        t_1345[k] = f_20 * slk_824[k]
                    + f_8 * smi0_837[k]
                    - f_9 * smi1_837[k]
                    + f_3 * pc_y[k] * smk_1076[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pc_y, slk_825, slk_826, slk_827, smi0_838, \
                         smi0_839, smi1_838, smi1_839, smk_1077, smk_1078, \
                         smk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_20 * slk_825[k]
                    + f_10 * smi0_838[k]
                    - f_11 * smi1_838[k]
                    + f_3 * pc_y[k] * smk_1077[k];

        t_1347[k] = f_20 * slk_826[k]
                    + f_12 * smi0_839[k]
                    - f_13 * smi1_839[k]
                    + f_3 * pc_y[k] * smk_1078[k];

        t_1348[k] = f_20 * slk_827[k]
                    + f_3 * pc_y[k] * smk_1079[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pc_x, pc_y, pc_z, slk_791, slk_828, slk_1080, \
                         smi0_839, smi0_840, smi1_839, smi1_840, smk_1079, \
                         smk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = f_15 * slk_791[k]
                    + f_1 * smi0_839[k]
                    - f_2 * smi1_839[k]
                    + f_3 * pc_z[k] * smk_1079[k];

        t_1350[k] = f_16 * slk_1080[k]
                    + f_1 * smi0_840[k]
                    - f_2 * smi1_840[k]
                    + f_3 * pc_x[k] * smk_1080[k];

        t_1351[k] = f_19 * slk_828[k]
                    + f_3 * pc_y[k] * smk_1080[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, pc_x, pc_y, pc_z, slk_792, slk_830, slk_1083, \
                         smi0_843, smi1_843, smk_1080, smk_1082, \
                         smk_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_16 * slk_792[k]
                    + f_3 * pc_z[k] * smk_1080[k];

        t_1353[k] = f_16 * slk_1083[k]
                    + f_4 * smi0_843[k]
                    - f_5 * smi1_843[k]
                    + f_3 * pc_x[k] * smk_1083[k];

        t_1354[k] = f_19 * slk_830[k]
                    + f_3 * pc_y[k] * smk_1082[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pc_x, pc_z, slk_795, slk_1085, slk_1086, \
                         smi0_845, smi0_846, smi1_845, smi1_846, smk_1083, smk_1085, \
                         smk_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_16 * slk_1085[k]
                    + f_4 * smi0_845[k]
                    - f_5 * smi1_845[k]
                    + f_3 * pc_x[k] * smk_1085[k];

        t_1356[k] = f_16 * slk_1086[k]
                    + f_6 * smi0_846[k]
                    - f_7 * smi1_846[k]
                    + f_3 * pc_x[k] * smk_1086[k];

        t_1357[k] = f_16 * slk_795[k]
                    + f_3 * pc_z[k] * smk_1083[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, pc_x, pc_y, slk_833, slk_1089, slk_1090, \
                         smi0_849, smi0_850, smi1_849, smi1_850, smk_1085, smk_1089, \
                         smk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_19 * slk_833[k]
                    + f_3 * pc_y[k] * smk_1085[k];

        t_1359[k] = f_16 * slk_1089[k]
                    + f_6 * smi0_849[k]
                    - f_7 * smi1_849[k]
                    + f_3 * pc_x[k] * smk_1089[k];

        t_1360[k] = f_16 * slk_1090[k]
                    + f_8 * smi0_850[k]
                    - f_9 * smi1_850[k]
                    + f_3 * pc_x[k] * smk_1090[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, pc_x, pc_y, pc_z, slk_798, slk_837, slk_1092, \
                         smi0_852, smi1_852, smk_1086, smk_1089, \
                         smk_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = f_16 * slk_798[k]
                    + f_3 * pc_z[k] * smk_1086[k];

        t_1362[k] = f_16 * slk_1092[k]
                    + f_8 * smi0_852[k]
                    - f_9 * smi1_852[k]
                    + f_3 * pc_x[k] * smk_1092[k];

        t_1363[k] = f_19 * slk_837[k]
                    + f_3 * pc_y[k] * smk_1089[k];
    }

#pragma omp simd aligned(t_1364, t_1365, t_1366, pc_x, pc_z, slk_802, slk_1094, slk_1095, \
                         smi0_854, smi0_855, smi1_854, smi1_855, smk_1090, smk_1094, \
                         smk_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1364[k] = f_16 * slk_1094[k]
                    + f_8 * smi0_854[k]
                    - f_9 * smi1_854[k]
                    + f_3 * pc_x[k] * smk_1094[k];

        t_1365[k] = f_16 * slk_1095[k]
                    + f_10 * smi0_855[k]
                    - f_11 * smi1_855[k]
                    + f_3 * pc_x[k] * smk_1095[k];

        t_1366[k] = f_16 * slk_802[k]
                    + f_3 * pc_z[k] * smk_1090[k];
    }

#pragma omp simd aligned(t_1367, t_1368, t_1369, pc_x, pc_y, slk_842, slk_1097, slk_1098, \
                         smi0_857, smi0_858, smi1_857, smi1_858, smk_1094, smk_1097, \
                         smk_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1367[k] = f_16 * slk_1097[k]
                    + f_10 * smi0_857[k]
                    - f_11 * smi1_857[k]
                    + f_3 * pc_x[k] * smk_1097[k];

        t_1368[k] = f_16 * slk_1098[k]
                    + f_10 * smi0_858[k]
                    - f_11 * smi1_858[k]
                    + f_3 * pc_x[k] * smk_1098[k];

        t_1369[k] = f_19 * slk_842[k]
                    + f_3 * pc_y[k] * smk_1094[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, pc_x, pc_z, slk_807, slk_1100, slk_1101, \
                         smi0_860, smi0_861, smi1_860, smi1_861, smk_1095, smk_1100, \
                         smk_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_16 * slk_1100[k]
                    + f_10 * smi0_860[k]
                    - f_11 * smi1_860[k]
                    + f_3 * pc_x[k] * smk_1100[k];

        t_1371[k] = f_16 * slk_1101[k]
                    + f_12 * smi0_861[k]
                    - f_13 * smi1_861[k]
                    + f_3 * pc_x[k] * smk_1101[k];

        t_1372[k] = f_16 * slk_807[k]
                    + f_3 * pc_z[k] * smk_1095[k];
    }

#pragma omp simd aligned(t_1373, t_1374, t_1375, pc_x, slk_1103, slk_1104, slk_1105, smi0_863, \
                         smi0_864, smi0_865, smi1_863, smi1_864, smi1_865, smk_1103, smk_1104, \
                         smk_1105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1373[k] = f_16 * slk_1103[k]
                    + f_12 * smi0_863[k]
                    - f_13 * smi1_863[k]
                    + f_3 * pc_x[k] * smk_1103[k];

        t_1374[k] = f_16 * slk_1104[k]
                    + f_12 * smi0_864[k]
                    - f_13 * smi1_864[k]
                    + f_3 * pc_x[k] * smk_1104[k];

        t_1375[k] = f_16 * slk_1105[k]
                    + f_12 * smi0_865[k]
                    - f_13 * smi1_865[k]
                    + f_3 * pc_x[k] * smk_1105[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, t_1379, pc_x, pc_y, slk_848, slk_1107, \
                         slk_1108, slk_1109, smi0_867, smi1_867, smk_1100, smk_1107, smk_1108, \
                         smk_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_19 * slk_848[k]
                    + f_3 * pc_y[k] * smk_1100[k];

        t_1377[k] = f_16 * slk_1107[k]
                    + f_12 * smi0_867[k]
                    - f_13 * smi1_867[k]
                    + f_3 * pc_x[k] * smk_1107[k];

        t_1378[k] = f_16 * slk_1108[k]
                    + f_3 * pc_x[k] * smk_1108[k];

        t_1379[k] = f_16 * slk_1109[k]
                    + f_3 * pc_x[k] * smk_1109[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, t_1383, t_1384, pc_x, slk_1110, slk_1111, \
                         slk_1112, slk_1113, slk_1114, smk_1110, smk_1111, smk_1112, smk_1113, \
                         smk_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_16 * slk_1110[k]
                    + f_3 * pc_x[k] * smk_1110[k];

        t_1381[k] = f_16 * slk_1111[k]
                    + f_3 * pc_x[k] * smk_1111[k];

        t_1382[k] = f_16 * slk_1112[k]
                    + f_3 * pc_x[k] * smk_1112[k];

        t_1383[k] = f_16 * slk_1113[k]
                    + f_3 * pc_x[k] * smk_1113[k];

        t_1384[k] = f_16 * slk_1114[k]
                    + f_3 * pc_x[k] * smk_1114[k];
    }

#pragma omp simd aligned(t_1385, t_1386, t_1387, pc_x, pc_y, pc_z, slk_820, slk_856, slk_1115, \
                         smi0_861, smi1_861, smk_1108, smk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1385[k] = f_16 * slk_1115[k]
                    + f_3 * pc_x[k] * smk_1115[k];

        t_1386[k] = f_19 * slk_856[k]
                    + f_1 * smi0_861[k]
                    - f_2 * smi1_861[k]
                    + f_3 * pc_y[k] * smk_1108[k];

        t_1387[k] = f_16 * slk_820[k]
                    + f_3 * pc_z[k] * smk_1108[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, pc_y, slk_858, slk_859, slk_860, smi0_863, \
                         smi0_864, smi0_865, smi1_863, smi1_864, smi1_865, smk_1110, smk_1111, \
                         smk_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = f_19 * slk_858[k]
                    + f_4 * smi0_863[k]
                    - f_5 * smi1_863[k]
                    + f_3 * pc_y[k] * smk_1110[k];

        t_1389[k] = f_19 * slk_859[k]
                    + f_6 * smi0_864[k]
                    - f_7 * smi1_864[k]
                    + f_3 * pc_y[k] * smk_1111[k];

        t_1390[k] = f_19 * slk_860[k]
                    + f_8 * smi0_865[k]
                    - f_9 * smi1_865[k]
                    + f_3 * pc_y[k] * smk_1112[k];
    }

#pragma omp simd aligned(t_1391, t_1392, t_1393, pc_y, slk_861, slk_862, slk_863, smi0_866, \
                         smi0_867, smi1_866, smi1_867, smk_1113, smk_1114, \
                         smk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1391[k] = f_19 * slk_861[k]
                    + f_10 * smi0_866[k]
                    - f_11 * smi1_866[k]
                    + f_3 * pc_y[k] * smk_1113[k];

        t_1392[k] = f_19 * slk_862[k]
                    + f_12 * smi0_867[k]
                    - f_13 * smi1_867[k]
                    + f_3 * pc_y[k] * smk_1114[k];

        t_1393[k] = f_19 * slk_863[k]
                    + f_3 * pc_y[k] * smk_1115[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, pc_x, pc_y, pc_z, slk_827, slk_864, slk_1116, \
                         smi0_867, smi0_868, smi1_867, smi1_868, smk_1115, \
                         smk_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_16 * slk_827[k]
                    + f_1 * smi0_867[k]
                    - f_2 * smi1_867[k]
                    + f_3 * pc_z[k] * smk_1115[k];

        t_1395[k] = f_16 * slk_1116[k]
                    + f_1 * smi0_868[k]
                    - f_2 * smi1_868[k]
                    + f_3 * pc_x[k] * smk_1116[k];

        t_1396[k] = f_18 * slk_864[k]
                    + f_3 * pc_y[k] * smk_1116[k];
    }

#pragma omp simd aligned(t_1397, t_1398, t_1399, pc_x, pc_y, pc_z, slk_828, slk_866, slk_1119, \
                         smi0_871, smi1_871, smk_1116, smk_1118, \
                         smk_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1397[k] = f_17 * slk_828[k]
                    + f_3 * pc_z[k] * smk_1116[k];

        t_1398[k] = f_16 * slk_1119[k]
                    + f_4 * smi0_871[k]
                    - f_5 * smi1_871[k]
                    + f_3 * pc_x[k] * smk_1119[k];

        t_1399[k] = f_18 * slk_866[k]
                    + f_3 * pc_y[k] * smk_1118[k];
    }

#pragma omp simd aligned(t_1400, t_1401, t_1402, pc_x, pc_z, slk_831, slk_1121, slk_1122, \
                         smi0_873, smi0_874, smi1_873, smi1_874, smk_1119, smk_1121, \
                         smk_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1400[k] = f_16 * slk_1121[k]
                    + f_4 * smi0_873[k]
                    - f_5 * smi1_873[k]
                    + f_3 * pc_x[k] * smk_1121[k];

        t_1401[k] = f_16 * slk_1122[k]
                    + f_6 * smi0_874[k]
                    - f_7 * smi1_874[k]
                    + f_3 * pc_x[k] * smk_1122[k];

        t_1402[k] = f_17 * slk_831[k]
                    + f_3 * pc_z[k] * smk_1119[k];
    }

#pragma omp simd aligned(t_1403, t_1404, t_1405, pc_x, pc_y, slk_869, slk_1125, slk_1126, \
                         smi0_877, smi0_878, smi1_877, smi1_878, smk_1121, smk_1125, \
                         smk_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1403[k] = f_18 * slk_869[k]
                    + f_3 * pc_y[k] * smk_1121[k];

        t_1404[k] = f_16 * slk_1125[k]
                    + f_6 * smi0_877[k]
                    - f_7 * smi1_877[k]
                    + f_3 * pc_x[k] * smk_1125[k];

        t_1405[k] = f_16 * slk_1126[k]
                    + f_8 * smi0_878[k]
                    - f_9 * smi1_878[k]
                    + f_3 * pc_x[k] * smk_1126[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, pc_x, pc_y, pc_z, slk_834, slk_873, slk_1128, \
                         smi0_880, smi1_880, smk_1122, smk_1125, \
                         smk_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_17 * slk_834[k]
                    + f_3 * pc_z[k] * smk_1122[k];

        t_1407[k] = f_16 * slk_1128[k]
                    + f_8 * smi0_880[k]
                    - f_9 * smi1_880[k]
                    + f_3 * pc_x[k] * smk_1128[k];

        t_1408[k] = f_18 * slk_873[k]
                    + f_3 * pc_y[k] * smk_1125[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pc_x, pc_z, slk_838, slk_1130, slk_1131, \
                         smi0_882, smi0_883, smi1_882, smi1_883, smk_1126, smk_1130, \
                         smk_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_16 * slk_1130[k]
                    + f_8 * smi0_882[k]
                    - f_9 * smi1_882[k]
                    + f_3 * pc_x[k] * smk_1130[k];

        t_1410[k] = f_16 * slk_1131[k]
                    + f_10 * smi0_883[k]
                    - f_11 * smi1_883[k]
                    + f_3 * pc_x[k] * smk_1131[k];

        t_1411[k] = f_17 * slk_838[k]
                    + f_3 * pc_z[k] * smk_1126[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, pc_x, pc_y, slk_878, slk_1133, slk_1134, \
                         smi0_885, smi0_886, smi1_885, smi1_886, smk_1130, smk_1133, \
                         smk_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_16 * slk_1133[k]
                    + f_10 * smi0_885[k]
                    - f_11 * smi1_885[k]
                    + f_3 * pc_x[k] * smk_1133[k];

        t_1413[k] = f_16 * slk_1134[k]
                    + f_10 * smi0_886[k]
                    - f_11 * smi1_886[k]
                    + f_3 * pc_x[k] * smk_1134[k];

        t_1414[k] = f_18 * slk_878[k]
                    + f_3 * pc_y[k] * smk_1130[k];
    }

#pragma omp simd aligned(t_1415, t_1416, t_1417, pc_x, pc_z, slk_843, slk_1136, slk_1137, \
                         smi0_888, smi0_889, smi1_888, smi1_889, smk_1131, smk_1136, \
                         smk_1137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1415[k] = f_16 * slk_1136[k]
                    + f_10 * smi0_888[k]
                    - f_11 * smi1_888[k]
                    + f_3 * pc_x[k] * smk_1136[k];

        t_1416[k] = f_16 * slk_1137[k]
                    + f_12 * smi0_889[k]
                    - f_13 * smi1_889[k]
                    + f_3 * pc_x[k] * smk_1137[k];

        t_1417[k] = f_17 * slk_843[k]
                    + f_3 * pc_z[k] * smk_1131[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, pc_x, slk_1139, slk_1140, slk_1141, smi0_891, \
                         smi0_892, smi0_893, smi1_891, smi1_892, smi1_893, smk_1139, smk_1140, \
                         smk_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_16 * slk_1139[k]
                    + f_12 * smi0_891[k]
                    - f_13 * smi1_891[k]
                    + f_3 * pc_x[k] * smk_1139[k];

        t_1419[k] = f_16 * slk_1140[k]
                    + f_12 * smi0_892[k]
                    - f_13 * smi1_892[k]
                    + f_3 * pc_x[k] * smk_1140[k];

        t_1420[k] = f_16 * slk_1141[k]
                    + f_12 * smi0_893[k]
                    - f_13 * smi1_893[k]
                    + f_3 * pc_x[k] * smk_1141[k];
    }

#pragma omp simd aligned(t_1421, t_1422, t_1423, t_1424, pc_x, pc_y, slk_884, slk_1143, \
                         slk_1144, slk_1145, smi0_895, smi1_895, smk_1136, smk_1143, smk_1144, \
                         smk_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1421[k] = f_18 * slk_884[k]
                    + f_3 * pc_y[k] * smk_1136[k];

        t_1422[k] = f_16 * slk_1143[k]
                    + f_12 * smi0_895[k]
                    - f_13 * smi1_895[k]
                    + f_3 * pc_x[k] * smk_1143[k];

        t_1423[k] = f_16 * slk_1144[k]
                    + f_3 * pc_x[k] * smk_1144[k];

        t_1424[k] = f_16 * slk_1145[k]
                    + f_3 * pc_x[k] * smk_1145[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t slk, const size_t smi0,
                                                           const size_t smi1, const size_t smk,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk_856 = buffer.data(slk + 856);
    const auto *slk_863 = buffer.data(slk + 863);
    const auto *slk_864 = buffer.data(slk + 864);
    const auto *slk_867 = buffer.data(slk + 867);
    const auto *slk_870 = buffer.data(slk + 870);
    const auto *slk_874 = buffer.data(slk + 874);
    const auto *slk_879 = buffer.data(slk + 879);
    const auto *slk_892 = buffer.data(slk + 892);
    const auto *slk_894 = buffer.data(slk + 894);
    const auto *slk_895 = buffer.data(slk + 895);
    const auto *slk_896 = buffer.data(slk + 896);
    const auto *slk_897 = buffer.data(slk + 897);
    const auto *slk_898 = buffer.data(slk + 898);
    const auto *slk_899 = buffer.data(slk + 899);
    const auto *slk_900 = buffer.data(slk + 900);
    const auto *slk_902 = buffer.data(slk + 902);
    const auto *slk_903 = buffer.data(slk + 903);
    const auto *slk_905 = buffer.data(slk + 905);
    const auto *slk_906 = buffer.data(slk + 906);
    const auto *slk_909 = buffer.data(slk + 909);
    const auto *slk_910 = buffer.data(slk + 910);
    const auto *slk_914 = buffer.data(slk + 914);
    const auto *slk_915 = buffer.data(slk + 915);
    const auto *slk_920 = buffer.data(slk + 920);
    const auto *slk_928 = buffer.data(slk + 928);
    const auto *slk_930 = buffer.data(slk + 930);
    const auto *slk_931 = buffer.data(slk + 931);
    const auto *slk_932 = buffer.data(slk + 932);
    const auto *slk_933 = buffer.data(slk + 933);
    const auto *slk_934 = buffer.data(slk + 934);
    const auto *slk_935 = buffer.data(slk + 935);
    const auto *slk_936 = buffer.data(slk + 936);
    const auto *slk_938 = buffer.data(slk + 938);
    const auto *slk_941 = buffer.data(slk + 941);
    const auto *slk_945 = buffer.data(slk + 945);
    const auto *slk_950 = buffer.data(slk + 950);
    const auto *slk_956 = buffer.data(slk + 956);
    const auto *slk_964 = buffer.data(slk + 964);
    const auto *slk_966 = buffer.data(slk + 966);
    const auto *slk_967 = buffer.data(slk + 967);
    const auto *slk_968 = buffer.data(slk + 968);
    const auto *slk_969 = buffer.data(slk + 969);
    const auto *slk_970 = buffer.data(slk + 970);
    const auto *slk_971 = buffer.data(slk + 971);
    const auto *slk_1146 = buffer.data(slk + 1146);
    const auto *slk_1147 = buffer.data(slk + 1147);
    const auto *slk_1148 = buffer.data(slk + 1148);
    const auto *slk_1149 = buffer.data(slk + 1149);
    const auto *slk_1150 = buffer.data(slk + 1150);
    const auto *slk_1151 = buffer.data(slk + 1151);
    const auto *slk_1152 = buffer.data(slk + 1152);
    const auto *slk_1155 = buffer.data(slk + 1155);
    const auto *slk_1157 = buffer.data(slk + 1157);
    const auto *slk_1158 = buffer.data(slk + 1158);
    const auto *slk_1161 = buffer.data(slk + 1161);
    const auto *slk_1162 = buffer.data(slk + 1162);
    const auto *slk_1164 = buffer.data(slk + 1164);
    const auto *slk_1166 = buffer.data(slk + 1166);
    const auto *slk_1167 = buffer.data(slk + 1167);
    const auto *slk_1169 = buffer.data(slk + 1169);
    const auto *slk_1170 = buffer.data(slk + 1170);
    const auto *slk_1172 = buffer.data(slk + 1172);
    const auto *slk_1173 = buffer.data(slk + 1173);
    const auto *slk_1175 = buffer.data(slk + 1175);
    const auto *slk_1176 = buffer.data(slk + 1176);
    const auto *slk_1177 = buffer.data(slk + 1177);
    const auto *slk_1179 = buffer.data(slk + 1179);
    const auto *slk_1180 = buffer.data(slk + 1180);
    const auto *slk_1181 = buffer.data(slk + 1181);
    const auto *slk_1182 = buffer.data(slk + 1182);
    const auto *slk_1183 = buffer.data(slk + 1183);
    const auto *slk_1184 = buffer.data(slk + 1184);
    const auto *slk_1185 = buffer.data(slk + 1185);
    const auto *slk_1186 = buffer.data(slk + 1186);
    const auto *slk_1187 = buffer.data(slk + 1187);
    const auto *slk_1188 = buffer.data(slk + 1188);
    const auto *slk_1191 = buffer.data(slk + 1191);
    const auto *slk_1193 = buffer.data(slk + 1193);
    const auto *slk_1194 = buffer.data(slk + 1194);
    const auto *slk_1197 = buffer.data(slk + 1197);
    const auto *slk_1198 = buffer.data(slk + 1198);
    const auto *slk_1200 = buffer.data(slk + 1200);
    const auto *slk_1202 = buffer.data(slk + 1202);
    const auto *slk_1203 = buffer.data(slk + 1203);
    const auto *slk_1205 = buffer.data(slk + 1205);
    const auto *slk_1206 = buffer.data(slk + 1206);
    const auto *slk_1208 = buffer.data(slk + 1208);
    const auto *slk_1209 = buffer.data(slk + 1209);
    const auto *slk_1211 = buffer.data(slk + 1211);
    const auto *slk_1212 = buffer.data(slk + 1212);
    const auto *slk_1213 = buffer.data(slk + 1213);
    const auto *slk_1215 = buffer.data(slk + 1215);
    const auto *slk_1216 = buffer.data(slk + 1216);
    const auto *slk_1217 = buffer.data(slk + 1217);
    const auto *slk_1218 = buffer.data(slk + 1218);
    const auto *slk_1219 = buffer.data(slk + 1219);
    const auto *slk_1220 = buffer.data(slk + 1220);
    const auto *slk_1221 = buffer.data(slk + 1221);
    const auto *slk_1222 = buffer.data(slk + 1222);
    const auto *slk_1223 = buffer.data(slk + 1223);

    const auto *smi0_889 = buffer.data(smi0 + 889);
    const auto *smi0_891 = buffer.data(smi0 + 891);
    const auto *smi0_892 = buffer.data(smi0 + 892);
    const auto *smi0_893 = buffer.data(smi0 + 893);
    const auto *smi0_894 = buffer.data(smi0 + 894);
    const auto *smi0_895 = buffer.data(smi0 + 895);
    const auto *smi0_896 = buffer.data(smi0 + 896);
    const auto *smi0_899 = buffer.data(smi0 + 899);
    const auto *smi0_901 = buffer.data(smi0 + 901);
    const auto *smi0_902 = buffer.data(smi0 + 902);
    const auto *smi0_905 = buffer.data(smi0 + 905);
    const auto *smi0_906 = buffer.data(smi0 + 906);
    const auto *smi0_908 = buffer.data(smi0 + 908);
    const auto *smi0_910 = buffer.data(smi0 + 910);
    const auto *smi0_911 = buffer.data(smi0 + 911);
    const auto *smi0_913 = buffer.data(smi0 + 913);
    const auto *smi0_914 = buffer.data(smi0 + 914);
    const auto *smi0_916 = buffer.data(smi0 + 916);
    const auto *smi0_917 = buffer.data(smi0 + 917);
    const auto *smi0_919 = buffer.data(smi0 + 919);
    const auto *smi0_920 = buffer.data(smi0 + 920);
    const auto *smi0_921 = buffer.data(smi0 + 921);
    const auto *smi0_922 = buffer.data(smi0 + 922);
    const auto *smi0_923 = buffer.data(smi0 + 923);
    const auto *smi0_924 = buffer.data(smi0 + 924);
    const auto *smi0_927 = buffer.data(smi0 + 927);
    const auto *smi0_929 = buffer.data(smi0 + 929);
    const auto *smi0_930 = buffer.data(smi0 + 930);
    const auto *smi0_933 = buffer.data(smi0 + 933);
    const auto *smi0_934 = buffer.data(smi0 + 934);
    const auto *smi0_936 = buffer.data(smi0 + 936);
    const auto *smi0_938 = buffer.data(smi0 + 938);
    const auto *smi0_939 = buffer.data(smi0 + 939);
    const auto *smi0_941 = buffer.data(smi0 + 941);
    const auto *smi0_942 = buffer.data(smi0 + 942);
    const auto *smi0_944 = buffer.data(smi0 + 944);
    const auto *smi0_945 = buffer.data(smi0 + 945);
    const auto *smi0_947 = buffer.data(smi0 + 947);
    const auto *smi0_948 = buffer.data(smi0 + 948);
    const auto *smi0_949 = buffer.data(smi0 + 949);
    const auto *smi0_950 = buffer.data(smi0 + 950);
    const auto *smi0_951 = buffer.data(smi0 + 951);

    const auto *smi1_889 = buffer.data(smi1 + 889);
    const auto *smi1_891 = buffer.data(smi1 + 891);
    const auto *smi1_892 = buffer.data(smi1 + 892);
    const auto *smi1_893 = buffer.data(smi1 + 893);
    const auto *smi1_894 = buffer.data(smi1 + 894);
    const auto *smi1_895 = buffer.data(smi1 + 895);
    const auto *smi1_896 = buffer.data(smi1 + 896);
    const auto *smi1_899 = buffer.data(smi1 + 899);
    const auto *smi1_901 = buffer.data(smi1 + 901);
    const auto *smi1_902 = buffer.data(smi1 + 902);
    const auto *smi1_905 = buffer.data(smi1 + 905);
    const auto *smi1_906 = buffer.data(smi1 + 906);
    const auto *smi1_908 = buffer.data(smi1 + 908);
    const auto *smi1_910 = buffer.data(smi1 + 910);
    const auto *smi1_911 = buffer.data(smi1 + 911);
    const auto *smi1_913 = buffer.data(smi1 + 913);
    const auto *smi1_914 = buffer.data(smi1 + 914);
    const auto *smi1_916 = buffer.data(smi1 + 916);
    const auto *smi1_917 = buffer.data(smi1 + 917);
    const auto *smi1_919 = buffer.data(smi1 + 919);
    const auto *smi1_920 = buffer.data(smi1 + 920);
    const auto *smi1_921 = buffer.data(smi1 + 921);
    const auto *smi1_922 = buffer.data(smi1 + 922);
    const auto *smi1_923 = buffer.data(smi1 + 923);
    const auto *smi1_924 = buffer.data(smi1 + 924);
    const auto *smi1_927 = buffer.data(smi1 + 927);
    const auto *smi1_929 = buffer.data(smi1 + 929);
    const auto *smi1_930 = buffer.data(smi1 + 930);
    const auto *smi1_933 = buffer.data(smi1 + 933);
    const auto *smi1_934 = buffer.data(smi1 + 934);
    const auto *smi1_936 = buffer.data(smi1 + 936);
    const auto *smi1_938 = buffer.data(smi1 + 938);
    const auto *smi1_939 = buffer.data(smi1 + 939);
    const auto *smi1_941 = buffer.data(smi1 + 941);
    const auto *smi1_942 = buffer.data(smi1 + 942);
    const auto *smi1_944 = buffer.data(smi1 + 944);
    const auto *smi1_945 = buffer.data(smi1 + 945);
    const auto *smi1_947 = buffer.data(smi1 + 947);
    const auto *smi1_948 = buffer.data(smi1 + 948);
    const auto *smi1_949 = buffer.data(smi1 + 949);
    const auto *smi1_950 = buffer.data(smi1 + 950);
    const auto *smi1_951 = buffer.data(smi1 + 951);

    const auto *smk_1144 = buffer.data(smk + 1144);
    const auto *smk_1146 = buffer.data(smk + 1146);
    const auto *smk_1147 = buffer.data(smk + 1147);
    const auto *smk_1148 = buffer.data(smk + 1148);
    const auto *smk_1149 = buffer.data(smk + 1149);
    const auto *smk_1150 = buffer.data(smk + 1150);
    const auto *smk_1151 = buffer.data(smk + 1151);
    const auto *smk_1152 = buffer.data(smk + 1152);
    const auto *smk_1154 = buffer.data(smk + 1154);
    const auto *smk_1155 = buffer.data(smk + 1155);
    const auto *smk_1157 = buffer.data(smk + 1157);
    const auto *smk_1158 = buffer.data(smk + 1158);
    const auto *smk_1161 = buffer.data(smk + 1161);
    const auto *smk_1162 = buffer.data(smk + 1162);
    const auto *smk_1164 = buffer.data(smk + 1164);
    const auto *smk_1166 = buffer.data(smk + 1166);
    const auto *smk_1167 = buffer.data(smk + 1167);
    const auto *smk_1169 = buffer.data(smk + 1169);
    const auto *smk_1170 = buffer.data(smk + 1170);
    const auto *smk_1172 = buffer.data(smk + 1172);
    const auto *smk_1173 = buffer.data(smk + 1173);
    const auto *smk_1175 = buffer.data(smk + 1175);
    const auto *smk_1176 = buffer.data(smk + 1176);
    const auto *smk_1177 = buffer.data(smk + 1177);
    const auto *smk_1179 = buffer.data(smk + 1179);
    const auto *smk_1180 = buffer.data(smk + 1180);
    const auto *smk_1181 = buffer.data(smk + 1181);
    const auto *smk_1182 = buffer.data(smk + 1182);
    const auto *smk_1183 = buffer.data(smk + 1183);
    const auto *smk_1184 = buffer.data(smk + 1184);
    const auto *smk_1185 = buffer.data(smk + 1185);
    const auto *smk_1186 = buffer.data(smk + 1186);
    const auto *smk_1187 = buffer.data(smk + 1187);
    const auto *smk_1188 = buffer.data(smk + 1188);
    const auto *smk_1190 = buffer.data(smk + 1190);
    const auto *smk_1191 = buffer.data(smk + 1191);
    const auto *smk_1193 = buffer.data(smk + 1193);
    const auto *smk_1194 = buffer.data(smk + 1194);
    const auto *smk_1197 = buffer.data(smk + 1197);
    const auto *smk_1198 = buffer.data(smk + 1198);
    const auto *smk_1200 = buffer.data(smk + 1200);
    const auto *smk_1202 = buffer.data(smk + 1202);
    const auto *smk_1203 = buffer.data(smk + 1203);
    const auto *smk_1205 = buffer.data(smk + 1205);
    const auto *smk_1206 = buffer.data(smk + 1206);
    const auto *smk_1208 = buffer.data(smk + 1208);
    const auto *smk_1209 = buffer.data(smk + 1209);
    const auto *smk_1211 = buffer.data(smk + 1211);
    const auto *smk_1212 = buffer.data(smk + 1212);
    const auto *smk_1213 = buffer.data(smk + 1213);
    const auto *smk_1215 = buffer.data(smk + 1215);
    const auto *smk_1216 = buffer.data(smk + 1216);
    const auto *smk_1217 = buffer.data(smk + 1217);
    const auto *smk_1218 = buffer.data(smk + 1218);
    const auto *smk_1219 = buffer.data(smk + 1219);
    const auto *smk_1220 = buffer.data(smk + 1220);
    const auto *smk_1221 = buffer.data(smk + 1221);
    const auto *smk_1222 = buffer.data(smk + 1222);
    const auto *smk_1223 = buffer.data(smk + 1223);

#pragma omp simd aligned(t_1425, t_1426, t_1427, t_1428, t_1429, pc_x, slk_1146, slk_1147, \
                         slk_1148, slk_1149, slk_1150, smk_1146, smk_1147, smk_1148, smk_1149, \
                         smk_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1425[k] = f_16 * slk_1146[k]
                    + f_3 * pc_x[k] * smk_1146[k];

        t_1426[k] = f_16 * slk_1147[k]
                    + f_3 * pc_x[k] * smk_1147[k];

        t_1427[k] = f_16 * slk_1148[k]
                    + f_3 * pc_x[k] * smk_1148[k];

        t_1428[k] = f_16 * slk_1149[k]
                    + f_3 * pc_x[k] * smk_1149[k];

        t_1429[k] = f_16 * slk_1150[k]
                    + f_3 * pc_x[k] * smk_1150[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, pc_x, pc_y, pc_z, slk_856, slk_892, slk_1151, \
                         smi0_889, smi1_889, smk_1144, smk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_16 * slk_1151[k]
                    + f_3 * pc_x[k] * smk_1151[k];

        t_1431[k] = f_18 * slk_892[k]
                    + f_1 * smi0_889[k]
                    - f_2 * smi1_889[k]
                    + f_3 * pc_y[k] * smk_1144[k];

        t_1432[k] = f_17 * slk_856[k]
                    + f_3 * pc_z[k] * smk_1144[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pc_y, slk_894, slk_895, slk_896, smi0_891, \
                         smi0_892, smi0_893, smi1_891, smi1_892, smi1_893, smk_1146, smk_1147, \
                         smk_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_18 * slk_894[k]
                    + f_4 * smi0_891[k]
                    - f_5 * smi1_891[k]
                    + f_3 * pc_y[k] * smk_1146[k];

        t_1434[k] = f_18 * slk_895[k]
                    + f_6 * smi0_892[k]
                    - f_7 * smi1_892[k]
                    + f_3 * pc_y[k] * smk_1147[k];

        t_1435[k] = f_18 * slk_896[k]
                    + f_8 * smi0_893[k]
                    - f_9 * smi1_893[k]
                    + f_3 * pc_y[k] * smk_1148[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, pc_y, slk_897, slk_898, slk_899, smi0_894, \
                         smi0_895, smi1_894, smi1_895, smk_1149, smk_1150, \
                         smk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_18 * slk_897[k]
                    + f_10 * smi0_894[k]
                    - f_11 * smi1_894[k]
                    + f_3 * pc_y[k] * smk_1149[k];

        t_1437[k] = f_18 * slk_898[k]
                    + f_12 * smi0_895[k]
                    - f_13 * smi1_895[k]
                    + f_3 * pc_y[k] * smk_1150[k];

        t_1438[k] = f_18 * slk_899[k]
                    + f_3 * pc_y[k] * smk_1151[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, pc_x, pc_y, pc_z, slk_863, slk_900, slk_1152, \
                         smi0_895, smi0_896, smi1_895, smi1_896, smk_1151, \
                         smk_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = f_17 * slk_863[k]
                    + f_1 * smi0_895[k]
                    - f_2 * smi1_895[k]
                    + f_3 * pc_z[k] * smk_1151[k];

        t_1440[k] = f_16 * slk_1152[k]
                    + f_1 * smi0_896[k]
                    - f_2 * smi1_896[k]
                    + f_3 * pc_x[k] * smk_1152[k];

        t_1441[k] = f_17 * slk_900[k]
                    + f_3 * pc_y[k] * smk_1152[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, pc_x, pc_y, pc_z, slk_864, slk_902, slk_1155, \
                         smi0_899, smi1_899, smk_1152, smk_1154, \
                         smk_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_18 * slk_864[k]
                    + f_3 * pc_z[k] * smk_1152[k];

        t_1443[k] = f_16 * slk_1155[k]
                    + f_4 * smi0_899[k]
                    - f_5 * smi1_899[k]
                    + f_3 * pc_x[k] * smk_1155[k];

        t_1444[k] = f_17 * slk_902[k]
                    + f_3 * pc_y[k] * smk_1154[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pc_x, pc_z, slk_867, slk_1157, slk_1158, \
                         smi0_901, smi0_902, smi1_901, smi1_902, smk_1155, smk_1157, \
                         smk_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_16 * slk_1157[k]
                    + f_4 * smi0_901[k]
                    - f_5 * smi1_901[k]
                    + f_3 * pc_x[k] * smk_1157[k];

        t_1446[k] = f_16 * slk_1158[k]
                    + f_6 * smi0_902[k]
                    - f_7 * smi1_902[k]
                    + f_3 * pc_x[k] * smk_1158[k];

        t_1447[k] = f_18 * slk_867[k]
                    + f_3 * pc_z[k] * smk_1155[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, pc_y, slk_905, slk_1161, slk_1162, \
                         smi0_905, smi0_906, smi1_905, smi1_906, smk_1157, smk_1161, \
                         smk_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_17 * slk_905[k]
                    + f_3 * pc_y[k] * smk_1157[k];

        t_1449[k] = f_16 * slk_1161[k]
                    + f_6 * smi0_905[k]
                    - f_7 * smi1_905[k]
                    + f_3 * pc_x[k] * smk_1161[k];

        t_1450[k] = f_16 * slk_1162[k]
                    + f_8 * smi0_906[k]
                    - f_9 * smi1_906[k]
                    + f_3 * pc_x[k] * smk_1162[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_x, pc_y, pc_z, slk_870, slk_909, slk_1164, \
                         smi0_908, smi1_908, smk_1158, smk_1161, \
                         smk_1164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_18 * slk_870[k]
                    + f_3 * pc_z[k] * smk_1158[k];

        t_1452[k] = f_16 * slk_1164[k]
                    + f_8 * smi0_908[k]
                    - f_9 * smi1_908[k]
                    + f_3 * pc_x[k] * smk_1164[k];

        t_1453[k] = f_17 * slk_909[k]
                    + f_3 * pc_y[k] * smk_1161[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, pc_z, slk_874, slk_1166, slk_1167, \
                         smi0_910, smi0_911, smi1_910, smi1_911, smk_1162, smk_1166, \
                         smk_1167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_16 * slk_1166[k]
                    + f_8 * smi0_910[k]
                    - f_9 * smi1_910[k]
                    + f_3 * pc_x[k] * smk_1166[k];

        t_1455[k] = f_16 * slk_1167[k]
                    + f_10 * smi0_911[k]
                    - f_11 * smi1_911[k]
                    + f_3 * pc_x[k] * smk_1167[k];

        t_1456[k] = f_18 * slk_874[k]
                    + f_3 * pc_z[k] * smk_1162[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pc_x, pc_y, slk_914, slk_1169, slk_1170, \
                         smi0_913, smi0_914, smi1_913, smi1_914, smk_1166, smk_1169, \
                         smk_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_16 * slk_1169[k]
                    + f_10 * smi0_913[k]
                    - f_11 * smi1_913[k]
                    + f_3 * pc_x[k] * smk_1169[k];

        t_1458[k] = f_16 * slk_1170[k]
                    + f_10 * smi0_914[k]
                    - f_11 * smi1_914[k]
                    + f_3 * pc_x[k] * smk_1170[k];

        t_1459[k] = f_17 * slk_914[k]
                    + f_3 * pc_y[k] * smk_1166[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, pc_x, pc_z, slk_879, slk_1172, slk_1173, \
                         smi0_916, smi0_917, smi1_916, smi1_917, smk_1167, smk_1172, \
                         smk_1173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_16 * slk_1172[k]
                    + f_10 * smi0_916[k]
                    - f_11 * smi1_916[k]
                    + f_3 * pc_x[k] * smk_1172[k];

        t_1461[k] = f_16 * slk_1173[k]
                    + f_12 * smi0_917[k]
                    - f_13 * smi1_917[k]
                    + f_3 * pc_x[k] * smk_1173[k];

        t_1462[k] = f_18 * slk_879[k]
                    + f_3 * pc_z[k] * smk_1167[k];
    }

#pragma omp simd aligned(t_1463, t_1464, t_1465, pc_x, slk_1175, slk_1176, slk_1177, smi0_919, \
                         smi0_920, smi0_921, smi1_919, smi1_920, smi1_921, smk_1175, smk_1176, \
                         smk_1177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1463[k] = f_16 * slk_1175[k]
                    + f_12 * smi0_919[k]
                    - f_13 * smi1_919[k]
                    + f_3 * pc_x[k] * smk_1175[k];

        t_1464[k] = f_16 * slk_1176[k]
                    + f_12 * smi0_920[k]
                    - f_13 * smi1_920[k]
                    + f_3 * pc_x[k] * smk_1176[k];

        t_1465[k] = f_16 * slk_1177[k]
                    + f_12 * smi0_921[k]
                    - f_13 * smi1_921[k]
                    + f_3 * pc_x[k] * smk_1177[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, t_1469, pc_x, pc_y, slk_920, slk_1179, \
                         slk_1180, slk_1181, smi0_923, smi1_923, smk_1172, smk_1179, smk_1180, \
                         smk_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_17 * slk_920[k]
                    + f_3 * pc_y[k] * smk_1172[k];

        t_1467[k] = f_16 * slk_1179[k]
                    + f_12 * smi0_923[k]
                    - f_13 * smi1_923[k]
                    + f_3 * pc_x[k] * smk_1179[k];

        t_1468[k] = f_16 * slk_1180[k]
                    + f_3 * pc_x[k] * smk_1180[k];

        t_1469[k] = f_16 * slk_1181[k]
                    + f_3 * pc_x[k] * smk_1181[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, t_1474, pc_x, slk_1182, slk_1183, \
                         slk_1184, slk_1185, slk_1186, smk_1182, smk_1183, smk_1184, smk_1185, \
                         smk_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = f_16 * slk_1182[k]
                    + f_3 * pc_x[k] * smk_1182[k];

        t_1471[k] = f_16 * slk_1183[k]
                    + f_3 * pc_x[k] * smk_1183[k];

        t_1472[k] = f_16 * slk_1184[k]
                    + f_3 * pc_x[k] * smk_1184[k];

        t_1473[k] = f_16 * slk_1185[k]
                    + f_3 * pc_x[k] * smk_1185[k];

        t_1474[k] = f_16 * slk_1186[k]
                    + f_3 * pc_x[k] * smk_1186[k];
    }

#pragma omp simd aligned(t_1475, t_1476, t_1477, pc_x, pc_y, pc_z, slk_892, slk_928, slk_1187, \
                         smi0_917, smi1_917, smk_1180, smk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1475[k] = f_16 * slk_1187[k]
                    + f_3 * pc_x[k] * smk_1187[k];

        t_1476[k] = f_17 * slk_928[k]
                    + f_1 * smi0_917[k]
                    - f_2 * smi1_917[k]
                    + f_3 * pc_y[k] * smk_1180[k];

        t_1477[k] = f_18 * slk_892[k]
                    + f_3 * pc_z[k] * smk_1180[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, pc_y, slk_930, slk_931, slk_932, smi0_919, \
                         smi0_920, smi0_921, smi1_919, smi1_920, smi1_921, smk_1182, smk_1183, \
                         smk_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_17 * slk_930[k]
                    + f_4 * smi0_919[k]
                    - f_5 * smi1_919[k]
                    + f_3 * pc_y[k] * smk_1182[k];

        t_1479[k] = f_17 * slk_931[k]
                    + f_6 * smi0_920[k]
                    - f_7 * smi1_920[k]
                    + f_3 * pc_y[k] * smk_1183[k];

        t_1480[k] = f_17 * slk_932[k]
                    + f_8 * smi0_921[k]
                    - f_9 * smi1_921[k]
                    + f_3 * pc_y[k] * smk_1184[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, pc_y, slk_933, slk_934, slk_935, smi0_922, \
                         smi0_923, smi1_922, smi1_923, smk_1185, smk_1186, \
                         smk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_17 * slk_933[k]
                    + f_10 * smi0_922[k]
                    - f_11 * smi1_922[k]
                    + f_3 * pc_y[k] * smk_1185[k];

        t_1482[k] = f_17 * slk_934[k]
                    + f_12 * smi0_923[k]
                    - f_13 * smi1_923[k]
                    + f_3 * pc_y[k] * smk_1186[k];

        t_1483[k] = f_17 * slk_935[k]
                    + f_3 * pc_y[k] * smk_1187[k];
    }

#pragma omp simd aligned(t_1484, t_1485, t_1486, pc_x, pc_y, pc_z, slk_899, slk_936, slk_1188, \
                         smi0_923, smi0_924, smi1_923, smi1_924, smk_1187, \
                         smk_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1484[k] = f_18 * slk_899[k]
                    + f_1 * smi0_923[k]
                    - f_2 * smi1_923[k]
                    + f_3 * pc_z[k] * smk_1187[k];

        t_1485[k] = f_16 * slk_1188[k]
                    + f_1 * smi0_924[k]
                    - f_2 * smi1_924[k]
                    + f_3 * pc_x[k] * smk_1188[k];

        t_1486[k] = f_16 * slk_936[k]
                    + f_3 * pc_y[k] * smk_1188[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pc_x, pc_y, pc_z, slk_900, slk_938, slk_1191, \
                         smi0_927, smi1_927, smk_1188, smk_1190, \
                         smk_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_19 * slk_900[k]
                    + f_3 * pc_z[k] * smk_1188[k];

        t_1488[k] = f_16 * slk_1191[k]
                    + f_4 * smi0_927[k]
                    - f_5 * smi1_927[k]
                    + f_3 * pc_x[k] * smk_1191[k];

        t_1489[k] = f_16 * slk_938[k]
                    + f_3 * pc_y[k] * smk_1190[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pc_x, pc_z, slk_903, slk_1193, slk_1194, \
                         smi0_929, smi0_930, smi1_929, smi1_930, smk_1191, smk_1193, \
                         smk_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = f_16 * slk_1193[k]
                    + f_4 * smi0_929[k]
                    - f_5 * smi1_929[k]
                    + f_3 * pc_x[k] * smk_1193[k];

        t_1491[k] = f_16 * slk_1194[k]
                    + f_6 * smi0_930[k]
                    - f_7 * smi1_930[k]
                    + f_3 * pc_x[k] * smk_1194[k];

        t_1492[k] = f_19 * slk_903[k]
                    + f_3 * pc_z[k] * smk_1191[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pc_x, pc_y, slk_941, slk_1197, slk_1198, \
                         smi0_933, smi0_934, smi1_933, smi1_934, smk_1193, smk_1197, \
                         smk_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = f_16 * slk_941[k]
                    + f_3 * pc_y[k] * smk_1193[k];

        t_1494[k] = f_16 * slk_1197[k]
                    + f_6 * smi0_933[k]
                    - f_7 * smi1_933[k]
                    + f_3 * pc_x[k] * smk_1197[k];

        t_1495[k] = f_16 * slk_1198[k]
                    + f_8 * smi0_934[k]
                    - f_9 * smi1_934[k]
                    + f_3 * pc_x[k] * smk_1198[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, pc_x, pc_y, pc_z, slk_906, slk_945, slk_1200, \
                         smi0_936, smi1_936, smk_1194, smk_1197, \
                         smk_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = f_19 * slk_906[k]
                    + f_3 * pc_z[k] * smk_1194[k];

        t_1497[k] = f_16 * slk_1200[k]
                    + f_8 * smi0_936[k]
                    - f_9 * smi1_936[k]
                    + f_3 * pc_x[k] * smk_1200[k];

        t_1498[k] = f_16 * slk_945[k]
                    + f_3 * pc_y[k] * smk_1197[k];
    }

#pragma omp simd aligned(t_1499, t_1500, t_1501, pc_x, pc_z, slk_910, slk_1202, slk_1203, \
                         smi0_938, smi0_939, smi1_938, smi1_939, smk_1198, smk_1202, \
                         smk_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1499[k] = f_16 * slk_1202[k]
                    + f_8 * smi0_938[k]
                    - f_9 * smi1_938[k]
                    + f_3 * pc_x[k] * smk_1202[k];

        t_1500[k] = f_16 * slk_1203[k]
                    + f_10 * smi0_939[k]
                    - f_11 * smi1_939[k]
                    + f_3 * pc_x[k] * smk_1203[k];

        t_1501[k] = f_19 * slk_910[k]
                    + f_3 * pc_z[k] * smk_1198[k];
    }

#pragma omp simd aligned(t_1502, t_1503, t_1504, pc_x, pc_y, slk_950, slk_1205, slk_1206, \
                         smi0_941, smi0_942, smi1_941, smi1_942, smk_1202, smk_1205, \
                         smk_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = f_16 * slk_1205[k]
                    + f_10 * smi0_941[k]
                    - f_11 * smi1_941[k]
                    + f_3 * pc_x[k] * smk_1205[k];

        t_1503[k] = f_16 * slk_1206[k]
                    + f_10 * smi0_942[k]
                    - f_11 * smi1_942[k]
                    + f_3 * pc_x[k] * smk_1206[k];

        t_1504[k] = f_16 * slk_950[k]
                    + f_3 * pc_y[k] * smk_1202[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, pc_x, pc_z, slk_915, slk_1208, slk_1209, \
                         smi0_944, smi0_945, smi1_944, smi1_945, smk_1203, smk_1208, \
                         smk_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = f_16 * slk_1208[k]
                    + f_10 * smi0_944[k]
                    - f_11 * smi1_944[k]
                    + f_3 * pc_x[k] * smk_1208[k];

        t_1506[k] = f_16 * slk_1209[k]
                    + f_12 * smi0_945[k]
                    - f_13 * smi1_945[k]
                    + f_3 * pc_x[k] * smk_1209[k];

        t_1507[k] = f_19 * slk_915[k]
                    + f_3 * pc_z[k] * smk_1203[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, pc_x, slk_1211, slk_1212, slk_1213, smi0_947, \
                         smi0_948, smi0_949, smi1_947, smi1_948, smi1_949, smk_1211, smk_1212, \
                         smk_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = f_16 * slk_1211[k]
                    + f_12 * smi0_947[k]
                    - f_13 * smi1_947[k]
                    + f_3 * pc_x[k] * smk_1211[k];

        t_1509[k] = f_16 * slk_1212[k]
                    + f_12 * smi0_948[k]
                    - f_13 * smi1_948[k]
                    + f_3 * pc_x[k] * smk_1212[k];

        t_1510[k] = f_16 * slk_1213[k]
                    + f_12 * smi0_949[k]
                    - f_13 * smi1_949[k]
                    + f_3 * pc_x[k] * smk_1213[k];
    }

#pragma omp simd aligned(t_1511, t_1512, t_1513, t_1514, pc_x, pc_y, slk_956, slk_1215, \
                         slk_1216, slk_1217, smi0_951, smi1_951, smk_1208, smk_1215, smk_1216, \
                         smk_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1511[k] = f_16 * slk_956[k]
                    + f_3 * pc_y[k] * smk_1208[k];

        t_1512[k] = f_16 * slk_1215[k]
                    + f_12 * smi0_951[k]
                    - f_13 * smi1_951[k]
                    + f_3 * pc_x[k] * smk_1215[k];

        t_1513[k] = f_16 * slk_1216[k]
                    + f_3 * pc_x[k] * smk_1216[k];

        t_1514[k] = f_16 * slk_1217[k]
                    + f_3 * pc_x[k] * smk_1217[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, t_1518, t_1519, pc_x, slk_1218, slk_1219, \
                         slk_1220, slk_1221, slk_1222, smk_1218, smk_1219, smk_1220, smk_1221, \
                         smk_1222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = f_16 * slk_1218[k]
                    + f_3 * pc_x[k] * smk_1218[k];

        t_1516[k] = f_16 * slk_1219[k]
                    + f_3 * pc_x[k] * smk_1219[k];

        t_1517[k] = f_16 * slk_1220[k]
                    + f_3 * pc_x[k] * smk_1220[k];

        t_1518[k] = f_16 * slk_1221[k]
                    + f_3 * pc_x[k] * smk_1221[k];

        t_1519[k] = f_16 * slk_1222[k]
                    + f_3 * pc_x[k] * smk_1222[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, pc_x, pc_y, pc_z, slk_928, slk_964, slk_1223, \
                         smi0_945, smi1_945, smk_1216, smk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_16 * slk_1223[k]
                    + f_3 * pc_x[k] * smk_1223[k];

        t_1521[k] = f_16 * slk_964[k]
                    + f_1 * smi0_945[k]
                    - f_2 * smi1_945[k]
                    + f_3 * pc_y[k] * smk_1216[k];

        t_1522[k] = f_19 * slk_928[k]
                    + f_3 * pc_z[k] * smk_1216[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pc_y, slk_966, slk_967, slk_968, smi0_947, \
                         smi0_948, smi0_949, smi1_947, smi1_948, smi1_949, smk_1218, smk_1219, \
                         smk_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_16 * slk_966[k]
                    + f_4 * smi0_947[k]
                    - f_5 * smi1_947[k]
                    + f_3 * pc_y[k] * smk_1218[k];

        t_1524[k] = f_16 * slk_967[k]
                    + f_6 * smi0_948[k]
                    - f_7 * smi1_948[k]
                    + f_3 * pc_y[k] * smk_1219[k];

        t_1525[k] = f_16 * slk_968[k]
                    + f_8 * smi0_949[k]
                    - f_9 * smi1_949[k]
                    + f_3 * pc_y[k] * smk_1220[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, pc_y, slk_969, slk_970, slk_971, smi0_950, \
                         smi0_951, smi1_950, smi1_951, smk_1221, smk_1222, \
                         smk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_16 * slk_969[k]
                    + f_10 * smi0_950[k]
                    - f_11 * smi1_950[k]
                    + f_3 * pc_y[k] * smk_1221[k];

        t_1527[k] = f_16 * slk_970[k]
                    + f_12 * smi0_951[k]
                    - f_13 * smi1_951[k]
                    + f_3 * pc_y[k] * smk_1222[k];

        t_1528[k] = f_16 * slk_971[k]
                    + f_3 * pc_y[k] * smk_1223[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smi0, const size_t smi1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_1215 = buffer.data(sll0 + 1215);
    const auto *sll0_1218 = buffer.data(sll0 + 1218);
    const auto *sll0_1220 = buffer.data(sll0 + 1220);
    const auto *sll0_1221 = buffer.data(sll0 + 1221);
    const auto *sll0_1224 = buffer.data(sll0 + 1224);
    const auto *sll0_1225 = buffer.data(sll0 + 1225);
    const auto *sll0_1227 = buffer.data(sll0 + 1227);
    const auto *sll0_1229 = buffer.data(sll0 + 1229);
    const auto *sll0_1230 = buffer.data(sll0 + 1230);
    const auto *sll0_1232 = buffer.data(sll0 + 1232);
    const auto *sll0_1233 = buffer.data(sll0 + 1233);
    const auto *sll0_1235 = buffer.data(sll0 + 1235);
    const auto *sll0_1236 = buffer.data(sll0 + 1236);
    const auto *sll0_1238 = buffer.data(sll0 + 1238);
    const auto *sll0_1239 = buffer.data(sll0 + 1239);
    const auto *sll0_1240 = buffer.data(sll0 + 1240);
    const auto *sll0_1242 = buffer.data(sll0 + 1242);
    const auto *sll0_1259 = buffer.data(sll0 + 1259);
    const auto *sll0_1620 = buffer.data(sll0 + 1620);
    const auto *sll0_1623 = buffer.data(sll0 + 1623);
    const auto *sll0_1625 = buffer.data(sll0 + 1625);
    const auto *sll0_1626 = buffer.data(sll0 + 1626);
    const auto *sll0_1629 = buffer.data(sll0 + 1629);
    const auto *sll0_1630 = buffer.data(sll0 + 1630);
    const auto *sll0_1632 = buffer.data(sll0 + 1632);
    const auto *sll0_1634 = buffer.data(sll0 + 1634);
    const auto *sll0_1635 = buffer.data(sll0 + 1635);
    const auto *sll0_1637 = buffer.data(sll0 + 1637);
    const auto *sll0_1638 = buffer.data(sll0 + 1638);
    const auto *sll0_1640 = buffer.data(sll0 + 1640);
    const auto *sll0_1641 = buffer.data(sll0 + 1641);

    const auto *slk_935 = buffer.data(slk + 935);
    const auto *slk_936 = buffer.data(slk + 936);
    const auto *slk_939 = buffer.data(slk + 939);
    const auto *slk_942 = buffer.data(slk + 942);
    const auto *slk_946 = buffer.data(slk + 946);
    const auto *slk_951 = buffer.data(slk + 951);
    const auto *slk_964 = buffer.data(slk + 964);
    const auto *slk_972 = buffer.data(slk + 972);
    const auto *slk_973 = buffer.data(slk + 973);
    const auto *slk_974 = buffer.data(slk + 974);
    const auto *slk_975 = buffer.data(slk + 975);
    const auto *slk_977 = buffer.data(slk + 977);
    const auto *slk_978 = buffer.data(slk + 978);
    const auto *slk_980 = buffer.data(slk + 980);
    const auto *slk_981 = buffer.data(slk + 981);
    const auto *slk_982 = buffer.data(slk + 982);
    const auto *slk_984 = buffer.data(slk + 984);
    const auto *slk_985 = buffer.data(slk + 985);
    const auto *slk_986 = buffer.data(slk + 986);
    const auto *slk_987 = buffer.data(slk + 987);
    const auto *slk_989 = buffer.data(slk + 989);
    const auto *slk_990 = buffer.data(slk + 990);
    const auto *slk_991 = buffer.data(slk + 991);
    const auto *slk_992 = buffer.data(slk + 992);
    const auto *slk_1000 = buffer.data(slk + 1000);
    const auto *slk_1002 = buffer.data(slk + 1002);
    const auto *slk_1003 = buffer.data(slk + 1003);
    const auto *slk_1004 = buffer.data(slk + 1004);
    const auto *slk_1005 = buffer.data(slk + 1005);
    const auto *slk_1006 = buffer.data(slk + 1006);
    const auto *slk_1007 = buffer.data(slk + 1007);
    const auto *slk_1008 = buffer.data(slk + 1008);
    const auto *slk_1010 = buffer.data(slk + 1010);
    const auto *slk_1013 = buffer.data(slk + 1013);
    const auto *slk_1017 = buffer.data(slk + 1017);
    const auto *slk_1022 = buffer.data(slk + 1022);
    const auto *slk_1252 = buffer.data(slk + 1252);
    const auto *slk_1253 = buffer.data(slk + 1253);
    const auto *slk_1254 = buffer.data(slk + 1254);
    const auto *slk_1255 = buffer.data(slk + 1255);
    const auto *slk_1256 = buffer.data(slk + 1256);
    const auto *slk_1257 = buffer.data(slk + 1257);
    const auto *slk_1258 = buffer.data(slk + 1258);
    const auto *slk_1259 = buffer.data(slk + 1259);
    const auto *slk_1260 = buffer.data(slk + 1260);
    const auto *slk_1263 = buffer.data(slk + 1263);
    const auto *slk_1265 = buffer.data(slk + 1265);
    const auto *slk_1266 = buffer.data(slk + 1266);
    const auto *slk_1269 = buffer.data(slk + 1269);
    const auto *slk_1270 = buffer.data(slk + 1270);
    const auto *slk_1272 = buffer.data(slk + 1272);
    const auto *slk_1274 = buffer.data(slk + 1274);
    const auto *slk_1275 = buffer.data(slk + 1275);
    const auto *slk_1277 = buffer.data(slk + 1277);
    const auto *slk_1278 = buffer.data(slk + 1278);
    const auto *slk_1280 = buffer.data(slk + 1280);
    const auto *slk_1281 = buffer.data(slk + 1281);
    const auto *slk_1283 = buffer.data(slk + 1283);
    const auto *slk_1284 = buffer.data(slk + 1284);
    const auto *slk_1285 = buffer.data(slk + 1285);
    const auto *slk_1287 = buffer.data(slk + 1287);
    const auto *slk_1288 = buffer.data(slk + 1288);
    const auto *slk_1289 = buffer.data(slk + 1289);
    const auto *slk_1290 = buffer.data(slk + 1290);
    const auto *slk_1291 = buffer.data(slk + 1291);
    const auto *slk_1292 = buffer.data(slk + 1292);
    const auto *slk_1293 = buffer.data(slk + 1293);
    const auto *slk_1294 = buffer.data(slk + 1294);
    const auto *slk_1295 = buffer.data(slk + 1295);
    const auto *slk_1296 = buffer.data(slk + 1296);
    const auto *slk_1299 = buffer.data(slk + 1299);
    const auto *slk_1301 = buffer.data(slk + 1301);
    const auto *slk_1302 = buffer.data(slk + 1302);
    const auto *slk_1305 = buffer.data(slk + 1305);
    const auto *slk_1306 = buffer.data(slk + 1306);
    const auto *slk_1308 = buffer.data(slk + 1308);
    const auto *slk_1310 = buffer.data(slk + 1310);
    const auto *slk_1311 = buffer.data(slk + 1311);
    const auto *slk_1313 = buffer.data(slk + 1313);
    const auto *slk_1314 = buffer.data(slk + 1314);
    const auto *slk_1316 = buffer.data(slk + 1316);
    const auto *slk_1317 = buffer.data(slk + 1317);

    const auto *sll1_1215 = buffer.data(sll1 + 1215);
    const auto *sll1_1218 = buffer.data(sll1 + 1218);
    const auto *sll1_1220 = buffer.data(sll1 + 1220);
    const auto *sll1_1221 = buffer.data(sll1 + 1221);
    const auto *sll1_1224 = buffer.data(sll1 + 1224);
    const auto *sll1_1225 = buffer.data(sll1 + 1225);
    const auto *sll1_1227 = buffer.data(sll1 + 1227);
    const auto *sll1_1229 = buffer.data(sll1 + 1229);
    const auto *sll1_1230 = buffer.data(sll1 + 1230);
    const auto *sll1_1232 = buffer.data(sll1 + 1232);
    const auto *sll1_1233 = buffer.data(sll1 + 1233);
    const auto *sll1_1235 = buffer.data(sll1 + 1235);
    const auto *sll1_1236 = buffer.data(sll1 + 1236);
    const auto *sll1_1238 = buffer.data(sll1 + 1238);
    const auto *sll1_1239 = buffer.data(sll1 + 1239);
    const auto *sll1_1240 = buffer.data(sll1 + 1240);
    const auto *sll1_1242 = buffer.data(sll1 + 1242);
    const auto *sll1_1259 = buffer.data(sll1 + 1259);
    const auto *sll1_1620 = buffer.data(sll1 + 1620);
    const auto *sll1_1623 = buffer.data(sll1 + 1623);
    const auto *sll1_1625 = buffer.data(sll1 + 1625);
    const auto *sll1_1626 = buffer.data(sll1 + 1626);
    const auto *sll1_1629 = buffer.data(sll1 + 1629);
    const auto *sll1_1630 = buffer.data(sll1 + 1630);
    const auto *sll1_1632 = buffer.data(sll1 + 1632);
    const auto *sll1_1634 = buffer.data(sll1 + 1634);
    const auto *sll1_1635 = buffer.data(sll1 + 1635);
    const auto *sll1_1637 = buffer.data(sll1 + 1637);
    const auto *sll1_1638 = buffer.data(sll1 + 1638);
    const auto *sll1_1640 = buffer.data(sll1 + 1640);
    const auto *sll1_1641 = buffer.data(sll1 + 1641);

    const auto *smi0_951 = buffer.data(smi0 + 951);
    const auto *smi0_973 = buffer.data(smi0 + 973);
    const auto *smi0_975 = buffer.data(smi0 + 975);
    const auto *smi0_976 = buffer.data(smi0 + 976);
    const auto *smi0_977 = buffer.data(smi0 + 977);
    const auto *smi0_978 = buffer.data(smi0 + 978);
    const auto *smi0_979 = buffer.data(smi0 + 979);
    const auto *smi0_980 = buffer.data(smi0 + 980);
    const auto *smi0_983 = buffer.data(smi0 + 983);
    const auto *smi0_985 = buffer.data(smi0 + 985);
    const auto *smi0_986 = buffer.data(smi0 + 986);
    const auto *smi0_989 = buffer.data(smi0 + 989);
    const auto *smi0_990 = buffer.data(smi0 + 990);
    const auto *smi0_992 = buffer.data(smi0 + 992);
    const auto *smi0_994 = buffer.data(smi0 + 994);
    const auto *smi0_995 = buffer.data(smi0 + 995);
    const auto *smi0_997 = buffer.data(smi0 + 997);
    const auto *smi0_998 = buffer.data(smi0 + 998);
    const auto *smi0_1000 = buffer.data(smi0 + 1000);
    const auto *smi0_1001 = buffer.data(smi0 + 1001);
    const auto *smi0_1003 = buffer.data(smi0 + 1003);
    const auto *smi0_1004 = buffer.data(smi0 + 1004);
    const auto *smi0_1005 = buffer.data(smi0 + 1005);
    const auto *smi0_1006 = buffer.data(smi0 + 1006);
    const auto *smi0_1007 = buffer.data(smi0 + 1007);

    const auto *smi1_951 = buffer.data(smi1 + 951);
    const auto *smi1_973 = buffer.data(smi1 + 973);
    const auto *smi1_975 = buffer.data(smi1 + 975);
    const auto *smi1_976 = buffer.data(smi1 + 976);
    const auto *smi1_977 = buffer.data(smi1 + 977);
    const auto *smi1_978 = buffer.data(smi1 + 978);
    const auto *smi1_979 = buffer.data(smi1 + 979);
    const auto *smi1_980 = buffer.data(smi1 + 980);
    const auto *smi1_983 = buffer.data(smi1 + 983);
    const auto *smi1_985 = buffer.data(smi1 + 985);
    const auto *smi1_986 = buffer.data(smi1 + 986);
    const auto *smi1_989 = buffer.data(smi1 + 989);
    const auto *smi1_990 = buffer.data(smi1 + 990);
    const auto *smi1_992 = buffer.data(smi1 + 992);
    const auto *smi1_994 = buffer.data(smi1 + 994);
    const auto *smi1_995 = buffer.data(smi1 + 995);
    const auto *smi1_997 = buffer.data(smi1 + 997);
    const auto *smi1_998 = buffer.data(smi1 + 998);
    const auto *smi1_1000 = buffer.data(smi1 + 1000);
    const auto *smi1_1001 = buffer.data(smi1 + 1001);
    const auto *smi1_1003 = buffer.data(smi1 + 1003);
    const auto *smi1_1004 = buffer.data(smi1 + 1004);
    const auto *smi1_1005 = buffer.data(smi1 + 1005);
    const auto *smi1_1006 = buffer.data(smi1 + 1006);
    const auto *smi1_1007 = buffer.data(smi1 + 1007);

    const auto *smk_1223 = buffer.data(smk + 1223);
    const auto *smk_1224 = buffer.data(smk + 1224);
    const auto *smk_1226 = buffer.data(smk + 1226);
    const auto *smk_1227 = buffer.data(smk + 1227);
    const auto *smk_1229 = buffer.data(smk + 1229);
    const auto *smk_1230 = buffer.data(smk + 1230);
    const auto *smk_1233 = buffer.data(smk + 1233);
    const auto *smk_1234 = buffer.data(smk + 1234);
    const auto *smk_1238 = buffer.data(smk + 1238);
    const auto *smk_1239 = buffer.data(smk + 1239);
    const auto *smk_1244 = buffer.data(smk + 1244);
    const auto *smk_1252 = buffer.data(smk + 1252);
    const auto *smk_1253 = buffer.data(smk + 1253);
    const auto *smk_1254 = buffer.data(smk + 1254);
    const auto *smk_1255 = buffer.data(smk + 1255);
    const auto *smk_1256 = buffer.data(smk + 1256);
    const auto *smk_1257 = buffer.data(smk + 1257);
    const auto *smk_1258 = buffer.data(smk + 1258);
    const auto *smk_1259 = buffer.data(smk + 1259);
    const auto *smk_1260 = buffer.data(smk + 1260);
    const auto *smk_1262 = buffer.data(smk + 1262);
    const auto *smk_1263 = buffer.data(smk + 1263);
    const auto *smk_1265 = buffer.data(smk + 1265);
    const auto *smk_1266 = buffer.data(smk + 1266);
    const auto *smk_1269 = buffer.data(smk + 1269);
    const auto *smk_1270 = buffer.data(smk + 1270);
    const auto *smk_1272 = buffer.data(smk + 1272);
    const auto *smk_1274 = buffer.data(smk + 1274);
    const auto *smk_1275 = buffer.data(smk + 1275);
    const auto *smk_1277 = buffer.data(smk + 1277);
    const auto *smk_1278 = buffer.data(smk + 1278);
    const auto *smk_1280 = buffer.data(smk + 1280);
    const auto *smk_1281 = buffer.data(smk + 1281);
    const auto *smk_1283 = buffer.data(smk + 1283);
    const auto *smk_1284 = buffer.data(smk + 1284);
    const auto *smk_1285 = buffer.data(smk + 1285);
    const auto *smk_1287 = buffer.data(smk + 1287);
    const auto *smk_1288 = buffer.data(smk + 1288);
    const auto *smk_1289 = buffer.data(smk + 1289);
    const auto *smk_1290 = buffer.data(smk + 1290);
    const auto *smk_1291 = buffer.data(smk + 1291);
    const auto *smk_1292 = buffer.data(smk + 1292);
    const auto *smk_1293 = buffer.data(smk + 1293);
    const auto *smk_1294 = buffer.data(smk + 1294);
    const auto *smk_1295 = buffer.data(smk + 1295);
    const auto *smk_1296 = buffer.data(smk + 1296);
    const auto *smk_1298 = buffer.data(smk + 1298);
    const auto *smk_1299 = buffer.data(smk + 1299);
    const auto *smk_1301 = buffer.data(smk + 1301);
    const auto *smk_1302 = buffer.data(smk + 1302);
    const auto *smk_1305 = buffer.data(smk + 1305);
    const auto *smk_1306 = buffer.data(smk + 1306);
    const auto *smk_1310 = buffer.data(smk + 1310);
    const auto *smk_1311 = buffer.data(smk + 1311);

#pragma omp simd aligned(t_1529, t_1530, t_1531, t_1532, pb_y, pc_y, pc_z, sll0_1215, slk_935, \
                         slk_936, slk_972, sll1_1215, smi0_951, smi1_951, smk_1223, \
                         smk_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_19 * slk_935[k]
                    + f_1 * smi0_951[k]
                    - f_2 * smi1_951[k]
                    + f_3 * pc_z[k] * smk_1223[k];

        t_1530[k] = pb_y[k] * sll0_1215[k]
                    - f_14 * pc_y[k] * sll1_1215[k];

        t_1531[k] = f_15 * slk_972[k]
                    + f_3 * pc_y[k] * smk_1224[k];

        t_1532[k] = f_20 * slk_936[k]
                    + f_3 * pc_z[k] * smk_1224[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, t_1536, pb_y, pc_y, sll0_1218, sll0_1220, \
                         sll0_1221, slk_973, slk_974, slk_975, sll1_1218, sll1_1220, \
                         sll1_1221, smk_1226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = pb_y[k] * sll0_1218[k]
                    + f_16 * slk_973[k]
                    - f_14 * pc_y[k] * sll1_1218[k];

        t_1534[k] = f_15 * slk_974[k]
                    + f_3 * pc_y[k] * smk_1226[k];

        t_1535[k] = pb_y[k] * sll0_1220[k]
                    - f_14 * pc_y[k] * sll1_1220[k];

        t_1536[k] = pb_y[k] * sll0_1221[k]
                    + f_17 * slk_975[k]
                    - f_14 * pc_y[k] * sll1_1221[k];
    }

#pragma omp simd aligned(t_1537, t_1538, t_1539, t_1540, pb_y, pc_y, pc_z, sll0_1224, \
                         sll0_1225, slk_939, slk_977, slk_978, sll1_1224, sll1_1225, smk_1227, \
                         smk_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1537[k] = f_20 * slk_939[k]
                    + f_3 * pc_z[k] * smk_1227[k];

        t_1538[k] = f_15 * slk_977[k]
                    + f_3 * pc_y[k] * smk_1229[k];

        t_1539[k] = pb_y[k] * sll0_1224[k]
                    - f_14 * pc_y[k] * sll1_1224[k];

        t_1540[k] = pb_y[k] * sll0_1225[k]
                    + f_18 * slk_978[k]
                    - f_14 * pc_y[k] * sll1_1225[k];
    }

#pragma omp simd aligned(t_1541, t_1542, t_1543, t_1544, pb_y, pc_y, pc_z, sll0_1227, \
                         sll0_1229, slk_942, slk_980, slk_981, sll1_1227, sll1_1229, smk_1230, \
                         smk_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1541[k] = f_20 * slk_942[k]
                    + f_3 * pc_z[k] * smk_1230[k];

        t_1542[k] = pb_y[k] * sll0_1227[k]
                    + f_16 * slk_980[k]
                    - f_14 * pc_y[k] * sll1_1227[k];

        t_1543[k] = f_15 * slk_981[k]
                    + f_3 * pc_y[k] * smk_1233[k];

        t_1544[k] = pb_y[k] * sll0_1229[k]
                    - f_14 * pc_y[k] * sll1_1229[k];
    }

#pragma omp simd aligned(t_1545, t_1546, t_1547, pb_y, pc_y, pc_z, sll0_1230, sll0_1232, \
                         slk_946, slk_982, slk_984, sll1_1230, sll1_1232, \
                         smk_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1545[k] = pb_y[k] * sll0_1230[k]
                    + f_19 * slk_982[k]
                    - f_14 * pc_y[k] * sll1_1230[k];

        t_1546[k] = f_20 * slk_946[k]
                    + f_3 * pc_z[k] * smk_1234[k];

        t_1547[k] = pb_y[k] * sll0_1232[k]
                    + f_17 * slk_984[k]
                    - f_14 * pc_y[k] * sll1_1232[k];
    }

#pragma omp simd aligned(t_1548, t_1549, t_1550, t_1551, pb_y, pc_y, sll0_1233, sll0_1235, \
                         sll0_1236, slk_985, slk_986, slk_987, sll1_1233, sll1_1235, \
                         sll1_1236, smk_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1548[k] = pb_y[k] * sll0_1233[k]
                    + f_16 * slk_985[k]
                    - f_14 * pc_y[k] * sll1_1233[k];

        t_1549[k] = f_15 * slk_986[k]
                    + f_3 * pc_y[k] * smk_1238[k];

        t_1550[k] = pb_y[k] * sll0_1235[k]
                    - f_14 * pc_y[k] * sll1_1235[k];

        t_1551[k] = pb_y[k] * sll0_1236[k]
                    + f_20 * slk_987[k]
                    - f_14 * pc_y[k] * sll1_1236[k];
    }

#pragma omp simd aligned(t_1552, t_1553, t_1554, pb_y, pc_y, pc_z, sll0_1238, sll0_1239, \
                         slk_951, slk_989, slk_990, sll1_1238, sll1_1239, \
                         smk_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1552[k] = f_20 * slk_951[k]
                    + f_3 * pc_z[k] * smk_1239[k];

        t_1553[k] = pb_y[k] * sll0_1238[k]
                    + f_18 * slk_989[k]
                    - f_14 * pc_y[k] * sll1_1238[k];

        t_1554[k] = pb_y[k] * sll0_1239[k]
                    + f_17 * slk_990[k]
                    - f_14 * pc_y[k] * sll1_1239[k];
    }

#pragma omp simd aligned(t_1555, t_1556, t_1557, t_1558, pb_y, pc_x, pc_y, sll0_1240, \
                         sll0_1242, slk_991, slk_992, slk_1252, sll1_1240, sll1_1242, \
                         smk_1244, smk_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1555[k] = pb_y[k] * sll0_1240[k]
                    + f_16 * slk_991[k]
                    - f_14 * pc_y[k] * sll1_1240[k];

        t_1556[k] = f_15 * slk_992[k]
                    + f_3 * pc_y[k] * smk_1244[k];

        t_1557[k] = pb_y[k] * sll0_1242[k]
                    - f_14 * pc_y[k] * sll1_1242[k];

        t_1558[k] = f_16 * slk_1252[k]
                    + f_3 * pc_x[k] * smk_1252[k];
    }

#pragma omp simd aligned(t_1559, t_1560, t_1561, t_1562, t_1563, pc_x, slk_1253, slk_1254, \
                         slk_1255, slk_1256, slk_1257, smk_1253, smk_1254, smk_1255, smk_1256, \
                         smk_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1559[k] = f_16 * slk_1253[k]
                    + f_3 * pc_x[k] * smk_1253[k];

        t_1560[k] = f_16 * slk_1254[k]
                    + f_3 * pc_x[k] * smk_1254[k];

        t_1561[k] = f_16 * slk_1255[k]
                    + f_3 * pc_x[k] * smk_1255[k];

        t_1562[k] = f_16 * slk_1256[k]
                    + f_3 * pc_x[k] * smk_1256[k];

        t_1563[k] = f_16 * slk_1257[k]
                    + f_3 * pc_x[k] * smk_1257[k];
    }

#pragma omp simd aligned(t_1564, t_1565, t_1566, t_1567, pc_x, pc_y, pc_z, slk_964, slk_1000, \
                         slk_1258, slk_1259, smi0_973, smi1_973, smk_1252, smk_1258, \
                         smk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1564[k] = f_16 * slk_1258[k]
                    + f_3 * pc_x[k] * smk_1258[k];

        t_1565[k] = f_16 * slk_1259[k]
                    + f_3 * pc_x[k] * smk_1259[k];

        t_1566[k] = f_15 * slk_1000[k]
                    + f_1 * smi0_973[k]
                    - f_2 * smi1_973[k]
                    + f_3 * pc_y[k] * smk_1252[k];

        t_1567[k] = f_20 * slk_964[k]
                    + f_3 * pc_z[k] * smk_1252[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, pc_y, slk_1002, slk_1003, slk_1004, smi0_975, \
                         smi0_976, smi0_977, smi1_975, smi1_976, smi1_977, smk_1254, smk_1255, \
                         smk_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = f_15 * slk_1002[k]
                    + f_4 * smi0_975[k]
                    - f_5 * smi1_975[k]
                    + f_3 * pc_y[k] * smk_1254[k];

        t_1569[k] = f_15 * slk_1003[k]
                    + f_6 * smi0_976[k]
                    - f_7 * smi1_976[k]
                    + f_3 * pc_y[k] * smk_1255[k];

        t_1570[k] = f_15 * slk_1004[k]
                    + f_8 * smi0_977[k]
                    - f_9 * smi1_977[k]
                    + f_3 * pc_y[k] * smk_1256[k];
    }

#pragma omp simd aligned(t_1571, t_1572, t_1573, pc_y, slk_1005, slk_1006, slk_1007, smi0_978, \
                         smi0_979, smi1_978, smi1_979, smk_1257, smk_1258, \
                         smk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1571[k] = f_15 * slk_1005[k]
                    + f_10 * smi0_978[k]
                    - f_11 * smi1_978[k]
                    + f_3 * pc_y[k] * smk_1257[k];

        t_1572[k] = f_15 * slk_1006[k]
                    + f_12 * smi0_979[k]
                    - f_13 * smi1_979[k]
                    + f_3 * pc_y[k] * smk_1258[k];

        t_1573[k] = f_15 * slk_1007[k]
                    + f_3 * pc_y[k] * smk_1259[k];
    }

#pragma omp simd aligned(t_1574, t_1575, t_1576, t_1577, pb_y, pc_x, pc_y, pc_z, sll0_1259, \
                         slk_972, slk_1260, sll1_1259, smi0_980, smi1_980, \
                         smk_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1574[k] = pb_y[k] * sll0_1259[k]
                    - f_14 * pc_y[k] * sll1_1259[k];

        t_1575[k] = f_16 * slk_1260[k]
                    + f_1 * smi0_980[k]
                    - f_2 * smi1_980[k]
                    + f_3 * pc_x[k] * smk_1260[k];

        t_1576[k] = f_3 * pc_y[k] * smk_1260[k];

        t_1577[k] = f_22 * slk_972[k]
                    + f_3 * pc_z[k] * smk_1260[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pc_x, pc_y, slk_1263, slk_1265, smi0_983, \
                         smi0_985, smi1_983, smi1_985, smk_1262, smk_1263, \
                         smk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = f_16 * slk_1263[k]
                    + f_4 * smi0_983[k]
                    - f_5 * smi1_983[k]
                    + f_3 * pc_x[k] * smk_1263[k];

        t_1579[k] = f_3 * pc_y[k] * smk_1262[k];

        t_1580[k] = f_16 * slk_1265[k]
                    + f_4 * smi0_985[k]
                    - f_5 * smi1_985[k]
                    + f_3 * pc_x[k] * smk_1265[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, pc_x, pc_y, pc_z, slk_975, slk_1266, \
                         smi0_986, smi1_986, smk_1263, smk_1265, \
                         smk_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_16 * slk_1266[k]
                    + f_6 * smi0_986[k]
                    - f_7 * smi1_986[k]
                    + f_3 * pc_x[k] * smk_1266[k];

        t_1582[k] = f_22 * slk_975[k]
                    + f_3 * pc_z[k] * smk_1263[k];

        t_1583[k] = f_3 * pc_y[k] * smk_1265[k];
    }

#pragma omp simd aligned(t_1584, t_1585, t_1586, pc_x, pc_z, slk_978, slk_1269, slk_1270, \
                         smi0_989, smi0_990, smi1_989, smi1_990, smk_1266, smk_1269, \
                         smk_1270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1584[k] = f_16 * slk_1269[k]
                    + f_6 * smi0_989[k]
                    - f_7 * smi1_989[k]
                    + f_3 * pc_x[k] * smk_1269[k];

        t_1585[k] = f_16 * slk_1270[k]
                    + f_8 * smi0_990[k]
                    - f_9 * smi1_990[k]
                    + f_3 * pc_x[k] * smk_1270[k];

        t_1586[k] = f_22 * slk_978[k]
                    + f_3 * pc_z[k] * smk_1266[k];
    }

#pragma omp simd aligned(t_1587, t_1588, t_1589, pc_x, pc_y, slk_1272, slk_1274, smi0_992, \
                         smi0_994, smi1_992, smi1_994, smk_1269, smk_1272, \
                         smk_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1587[k] = f_16 * slk_1272[k]
                    + f_8 * smi0_992[k]
                    - f_9 * smi1_992[k]
                    + f_3 * pc_x[k] * smk_1272[k];

        t_1588[k] = f_3 * pc_y[k] * smk_1269[k];

        t_1589[k] = f_16 * slk_1274[k]
                    + f_8 * smi0_994[k]
                    - f_9 * smi1_994[k]
                    + f_3 * pc_x[k] * smk_1274[k];
    }

#pragma omp simd aligned(t_1590, t_1591, t_1592, pc_x, pc_z, slk_982, slk_1275, slk_1277, \
                         smi0_995, smi0_997, smi1_995, smi1_997, smk_1270, smk_1275, \
                         smk_1277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1590[k] = f_16 * slk_1275[k]
                    + f_10 * smi0_995[k]
                    - f_11 * smi1_995[k]
                    + f_3 * pc_x[k] * smk_1275[k];

        t_1591[k] = f_22 * slk_982[k]
                    + f_3 * pc_z[k] * smk_1270[k];

        t_1592[k] = f_16 * slk_1277[k]
                    + f_10 * smi0_997[k]
                    - f_11 * smi1_997[k]
                    + f_3 * pc_x[k] * smk_1277[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, pc_x, pc_y, slk_1278, slk_1280, smi0_998, \
                         smi0_1000, smi1_998, smi1_1000, smk_1274, smk_1278, \
                         smk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = f_16 * slk_1278[k]
                    + f_10 * smi0_998[k]
                    - f_11 * smi1_998[k]
                    + f_3 * pc_x[k] * smk_1278[k];

        t_1594[k] = f_3 * pc_y[k] * smk_1274[k];

        t_1595[k] = f_16 * slk_1280[k]
                    + f_10 * smi0_1000[k]
                    - f_11 * smi1_1000[k]
                    + f_3 * pc_x[k] * smk_1280[k];
    }

#pragma omp simd aligned(t_1596, t_1597, t_1598, pc_x, pc_z, slk_987, slk_1281, slk_1283, \
                         smi0_1001, smi0_1003, smi1_1001, smi1_1003, smk_1275, smk_1281, \
                         smk_1283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = f_16 * slk_1281[k]
                    + f_12 * smi0_1001[k]
                    - f_13 * smi1_1001[k]
                    + f_3 * pc_x[k] * smk_1281[k];

        t_1597[k] = f_22 * slk_987[k]
                    + f_3 * pc_z[k] * smk_1275[k];

        t_1598[k] = f_16 * slk_1283[k]
                    + f_12 * smi0_1003[k]
                    - f_13 * smi1_1003[k]
                    + f_3 * pc_x[k] * smk_1283[k];
    }

#pragma omp simd aligned(t_1599, t_1600, t_1601, pc_x, pc_y, slk_1284, slk_1285, smi0_1004, \
                         smi0_1005, smi1_1004, smi1_1005, smk_1280, smk_1284, \
                         smk_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1599[k] = f_16 * slk_1284[k]
                    + f_12 * smi0_1004[k]
                    - f_13 * smi1_1004[k]
                    + f_3 * pc_x[k] * smk_1284[k];

        t_1600[k] = f_16 * slk_1285[k]
                    + f_12 * smi0_1005[k]
                    - f_13 * smi1_1005[k]
                    + f_3 * pc_x[k] * smk_1285[k];

        t_1601[k] = f_3 * pc_y[k] * smk_1280[k];
    }

#pragma omp simd aligned(t_1602, t_1603, t_1604, t_1605, pc_x, slk_1287, slk_1288, slk_1289, \
                         slk_1290, smi0_1007, smi1_1007, smk_1287, smk_1288, smk_1289, \
                         smk_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1602[k] = f_16 * slk_1287[k]
                    + f_12 * smi0_1007[k]
                    - f_13 * smi1_1007[k]
                    + f_3 * pc_x[k] * smk_1287[k];

        t_1603[k] = f_16 * slk_1288[k]
                    + f_3 * pc_x[k] * smk_1288[k];

        t_1604[k] = f_16 * slk_1289[k]
                    + f_3 * pc_x[k] * smk_1289[k];

        t_1605[k] = f_16 * slk_1290[k]
                    + f_3 * pc_x[k] * smk_1290[k];
    }

#pragma omp simd aligned(t_1606, t_1607, t_1608, t_1609, t_1610, pc_x, slk_1291, slk_1292, \
                         slk_1293, slk_1294, slk_1295, smk_1291, smk_1292, smk_1293, smk_1294, \
                         smk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1606[k] = f_16 * slk_1291[k]
                    + f_3 * pc_x[k] * smk_1291[k];

        t_1607[k] = f_16 * slk_1292[k]
                    + f_3 * pc_x[k] * smk_1292[k];

        t_1608[k] = f_16 * slk_1293[k]
                    + f_3 * pc_x[k] * smk_1293[k];

        t_1609[k] = f_16 * slk_1294[k]
                    + f_3 * pc_x[k] * smk_1294[k];

        t_1610[k] = f_16 * slk_1295[k]
                    + f_3 * pc_x[k] * smk_1295[k];
    }

#pragma omp simd aligned(t_1611, t_1612, t_1613, t_1614, pc_y, pc_z, slk_1000, smi0_1001, \
                         smi0_1003, smi0_1004, smi1_1001, smi1_1003, smi1_1004, smk_1288, \
                         smk_1290, smk_1291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1611[k] = f_1 * smi0_1001[k]
                    - f_2 * smi1_1001[k]
                    + f_3 * pc_y[k] * smk_1288[k];

        t_1612[k] = f_22 * slk_1000[k]
                    + f_3 * pc_z[k] * smk_1288[k];

        t_1613[k] = f_4 * smi0_1003[k]
                    - f_5 * smi1_1003[k]
                    + f_3 * pc_y[k] * smk_1290[k];

        t_1614[k] = f_6 * smi0_1004[k]
                    - f_7 * smi1_1004[k]
                    + f_3 * pc_y[k] * smk_1291[k];
    }

#pragma omp simd aligned(t_1615, t_1616, t_1617, t_1618, pc_y, smi0_1005, smi0_1006, \
                         smi0_1007, smi1_1005, smi1_1006, smi1_1007, smk_1292, smk_1293, \
                         smk_1294, smk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1615[k] = f_8 * smi0_1005[k]
                    - f_9 * smi1_1005[k]
                    + f_3 * pc_y[k] * smk_1292[k];

        t_1616[k] = f_10 * smi0_1006[k]
                    - f_11 * smi1_1006[k]
                    + f_3 * pc_y[k] * smk_1293[k];

        t_1617[k] = f_12 * smi0_1007[k]
                    - f_13 * smi1_1007[k]
                    + f_3 * pc_y[k] * smk_1294[k];

        t_1618[k] = f_3 * pc_y[k] * smk_1295[k];
    }

#pragma omp simd aligned(t_1619, t_1620, t_1621, pb_x, pc_x, pc_y, pc_z, sll0_1620, slk_1007, \
                         slk_1008, slk_1296, sll1_1620, smi0_1007, smi1_1007, smk_1295, \
                         smk_1296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1619[k] = f_22 * slk_1007[k]
                    + f_1 * smi0_1007[k]
                    - f_2 * smi1_1007[k]
                    + f_3 * pc_z[k] * smk_1295[k];

        t_1620[k] = pb_x[k] * sll0_1620[k]
                    + f_21 * slk_1296[k]
                    - f_14 * pc_x[k] * sll1_1620[k];

        t_1621[k] = f_21 * slk_1008[k]
                    + f_3 * pc_y[k] * smk_1296[k];
    }

#pragma omp simd aligned(t_1622, t_1623, t_1624, pb_x, pc_x, pc_y, pc_z, sll0_1623, slk_1010, \
                         slk_1299, sll1_1623, smk_1296, smk_1298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1622[k] = f_3 * pc_z[k] * smk_1296[k];

        t_1623[k] = pb_x[k] * sll0_1623[k]
                    + f_20 * slk_1299[k]
                    - f_14 * pc_x[k] * sll1_1623[k];

        t_1624[k] = f_21 * slk_1010[k]
                    + f_3 * pc_y[k] * smk_1298[k];
    }

#pragma omp simd aligned(t_1625, t_1626, t_1627, pb_x, pc_x, pc_z, sll0_1625, sll0_1626, \
                         slk_1301, slk_1302, sll1_1625, sll1_1626, \
                         smk_1299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1625[k] = pb_x[k] * sll0_1625[k]
                    + f_20 * slk_1301[k]
                    - f_14 * pc_x[k] * sll1_1625[k];

        t_1626[k] = pb_x[k] * sll0_1626[k]
                    + f_19 * slk_1302[k]
                    - f_14 * pc_x[k] * sll1_1626[k];

        t_1627[k] = f_3 * pc_z[k] * smk_1299[k];
    }

#pragma omp simd aligned(t_1628, t_1629, t_1630, pb_x, pc_x, pc_y, sll0_1629, sll0_1630, \
                         slk_1013, slk_1305, slk_1306, sll1_1629, sll1_1630, \
                         smk_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1628[k] = f_21 * slk_1013[k]
                    + f_3 * pc_y[k] * smk_1301[k];

        t_1629[k] = pb_x[k] * sll0_1629[k]
                    + f_19 * slk_1305[k]
                    - f_14 * pc_x[k] * sll1_1629[k];

        t_1630[k] = pb_x[k] * sll0_1630[k]
                    + f_18 * slk_1306[k]
                    - f_14 * pc_x[k] * sll1_1630[k];
    }

#pragma omp simd aligned(t_1631, t_1632, t_1633, pb_x, pc_x, pc_y, pc_z, sll0_1632, slk_1017, \
                         slk_1308, sll1_1632, smk_1302, smk_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1631[k] = f_3 * pc_z[k] * smk_1302[k];

        t_1632[k] = pb_x[k] * sll0_1632[k]
                    + f_18 * slk_1308[k]
                    - f_14 * pc_x[k] * sll1_1632[k];

        t_1633[k] = f_21 * slk_1017[k]
                    + f_3 * pc_y[k] * smk_1305[k];
    }

#pragma omp simd aligned(t_1634, t_1635, t_1636, pb_x, pc_x, pc_z, sll0_1634, sll0_1635, \
                         slk_1310, slk_1311, sll1_1634, sll1_1635, \
                         smk_1306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1634[k] = pb_x[k] * sll0_1634[k]
                    + f_18 * slk_1310[k]
                    - f_14 * pc_x[k] * sll1_1634[k];

        t_1635[k] = pb_x[k] * sll0_1635[k]
                    + f_17 * slk_1311[k]
                    - f_14 * pc_x[k] * sll1_1635[k];

        t_1636[k] = f_3 * pc_z[k] * smk_1306[k];
    }

#pragma omp simd aligned(t_1637, t_1638, t_1639, pb_x, pc_x, pc_y, sll0_1637, sll0_1638, \
                         slk_1022, slk_1313, slk_1314, sll1_1637, sll1_1638, \
                         smk_1310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1637[k] = pb_x[k] * sll0_1637[k]
                    + f_17 * slk_1313[k]
                    - f_14 * pc_x[k] * sll1_1637[k];

        t_1638[k] = pb_x[k] * sll0_1638[k]
                    + f_17 * slk_1314[k]
                    - f_14 * pc_x[k] * sll1_1638[k];

        t_1639[k] = f_21 * slk_1022[k]
                    + f_3 * pc_y[k] * smk_1310[k];
    }

#pragma omp simd aligned(t_1640, t_1641, t_1642, pb_x, pc_x, pc_z, sll0_1640, sll0_1641, \
                         slk_1316, slk_1317, sll1_1640, sll1_1641, \
                         smk_1311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1640[k] = pb_x[k] * sll0_1640[k]
                    + f_17 * slk_1316[k]
                    - f_14 * pc_x[k] * sll1_1640[k];

        t_1641[k] = pb_x[k] * sll0_1641[k]
                    + f_16 * slk_1317[k]
                    - f_14 * pc_x[k] * sll1_1641[k];

        t_1642[k] = f_3 * pc_z[k] * smk_1311[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_1260 = buffer.data(sll0 + 1260);
    const auto *sll0_1263 = buffer.data(sll0 + 1263);
    const auto *sll0_1266 = buffer.data(sll0 + 1266);
    const auto *sll0_1270 = buffer.data(sll0 + 1270);
    const auto *sll0_1275 = buffer.data(sll0 + 1275);
    const auto *sll0_1281 = buffer.data(sll0 + 1281);
    const auto *sll0_1643 = buffer.data(sll0 + 1643);
    const auto *sll0_1644 = buffer.data(sll0 + 1644);
    const auto *sll0_1645 = buffer.data(sll0 + 1645);
    const auto *sll0_1647 = buffer.data(sll0 + 1647);
    const auto *sll0_1656 = buffer.data(sll0 + 1656);
    const auto *sll0_1658 = buffer.data(sll0 + 1658);
    const auto *sll0_1659 = buffer.data(sll0 + 1659);
    const auto *sll0_1660 = buffer.data(sll0 + 1660);
    const auto *sll0_1661 = buffer.data(sll0 + 1661);
    const auto *sll0_1662 = buffer.data(sll0 + 1662);
    const auto *sll0_1664 = buffer.data(sll0 + 1664);
    const auto *sll0_1670 = buffer.data(sll0 + 1670);
    const auto *sll0_1674 = buffer.data(sll0 + 1674);
    const auto *sll0_1677 = buffer.data(sll0 + 1677);
    const auto *sll0_1679 = buffer.data(sll0 + 1679);
    const auto *sll0_1682 = buffer.data(sll0 + 1682);
    const auto *sll0_1683 = buffer.data(sll0 + 1683);
    const auto *sll0_1685 = buffer.data(sll0 + 1685);
    const auto *sll0_1688 = buffer.data(sll0 + 1688);
    const auto *sll0_1689 = buffer.data(sll0 + 1689);
    const auto *sll0_1690 = buffer.data(sll0 + 1690);
    const auto *sll0_1692 = buffer.data(sll0 + 1692);
    const auto *sll0_1701 = buffer.data(sll0 + 1701);
    const auto *sll0_1703 = buffer.data(sll0 + 1703);
    const auto *sll0_1704 = buffer.data(sll0 + 1704);
    const auto *sll0_1705 = buffer.data(sll0 + 1705);
    const auto *sll0_1706 = buffer.data(sll0 + 1706);
    const auto *sll0_1707 = buffer.data(sll0 + 1707);
    const auto *sll0_1709 = buffer.data(sll0 + 1709);
    const auto *sll0_1710 = buffer.data(sll0 + 1710);
    const auto *sll0_1713 = buffer.data(sll0 + 1713);
    const auto *sll0_1715 = buffer.data(sll0 + 1715);
    const auto *sll0_1716 = buffer.data(sll0 + 1716);
    const auto *sll0_1719 = buffer.data(sll0 + 1719);
    const auto *sll0_1720 = buffer.data(sll0 + 1720);
    const auto *sll0_1722 = buffer.data(sll0 + 1722);
    const auto *sll0_1724 = buffer.data(sll0 + 1724);
    const auto *sll0_1725 = buffer.data(sll0 + 1725);
    const auto *sll0_1727 = buffer.data(sll0 + 1727);
    const auto *sll0_1728 = buffer.data(sll0 + 1728);
    const auto *sll0_1730 = buffer.data(sll0 + 1730);
    const auto *sll0_1731 = buffer.data(sll0 + 1731);
    const auto *sll0_1733 = buffer.data(sll0 + 1733);
    const auto *sll0_1734 = buffer.data(sll0 + 1734);
    const auto *sll0_1735 = buffer.data(sll0 + 1735);
    const auto *sll0_1737 = buffer.data(sll0 + 1737);
    const auto *sll0_1746 = buffer.data(sll0 + 1746);
    const auto *sll0_1748 = buffer.data(sll0 + 1748);
    const auto *sll0_1749 = buffer.data(sll0 + 1749);
    const auto *sll0_1750 = buffer.data(sll0 + 1750);
    const auto *sll0_1751 = buffer.data(sll0 + 1751);
    const auto *sll0_1752 = buffer.data(sll0 + 1752);
    const auto *sll0_1754 = buffer.data(sll0 + 1754);
    const auto *sll0_1755 = buffer.data(sll0 + 1755);
    const auto *sll0_1758 = buffer.data(sll0 + 1758);
    const auto *sll0_1760 = buffer.data(sll0 + 1760);
    const auto *sll0_1761 = buffer.data(sll0 + 1761);

    const auto *slk_1008 = buffer.data(slk + 1008);
    const auto *slk_1011 = buffer.data(slk + 1011);
    const auto *slk_1014 = buffer.data(slk + 1014);
    const auto *slk_1018 = buffer.data(slk + 1018);
    const auto *slk_1023 = buffer.data(slk + 1023);
    const auto *slk_1028 = buffer.data(slk + 1028);
    const auto *slk_1036 = buffer.data(slk + 1036);
    const auto *slk_1043 = buffer.data(slk + 1043);
    const auto *slk_1044 = buffer.data(slk + 1044);
    const auto *slk_1046 = buffer.data(slk + 1046);
    const auto *slk_1047 = buffer.data(slk + 1047);
    const auto *slk_1049 = buffer.data(slk + 1049);
    const auto *slk_1050 = buffer.data(slk + 1050);
    const auto *slk_1053 = buffer.data(slk + 1053);
    const auto *slk_1054 = buffer.data(slk + 1054);
    const auto *slk_1058 = buffer.data(slk + 1058);
    const auto *slk_1059 = buffer.data(slk + 1059);
    const auto *slk_1064 = buffer.data(slk + 1064);
    const auto *slk_1072 = buffer.data(slk + 1072);
    const auto *slk_1079 = buffer.data(slk + 1079);
    const auto *slk_1080 = buffer.data(slk + 1080);
    const auto *slk_1082 = buffer.data(slk + 1082);
    const auto *slk_1083 = buffer.data(slk + 1083);
    const auto *slk_1085 = buffer.data(slk + 1085);
    const auto *slk_1089 = buffer.data(slk + 1089);
    const auto *slk_1094 = buffer.data(slk + 1094);
    const auto *slk_1100 = buffer.data(slk + 1100);
    const auto *slk_1115 = buffer.data(slk + 1115);
    const auto *slk_1116 = buffer.data(slk + 1116);
    const auto *slk_1118 = buffer.data(slk + 1118);
    const auto *slk_1319 = buffer.data(slk + 1319);
    const auto *slk_1320 = buffer.data(slk + 1320);
    const auto *slk_1321 = buffer.data(slk + 1321);
    const auto *slk_1323 = buffer.data(slk + 1323);
    const auto *slk_1324 = buffer.data(slk + 1324);
    const auto *slk_1325 = buffer.data(slk + 1325);
    const auto *slk_1326 = buffer.data(slk + 1326);
    const auto *slk_1327 = buffer.data(slk + 1327);
    const auto *slk_1328 = buffer.data(slk + 1328);
    const auto *slk_1329 = buffer.data(slk + 1329);
    const auto *slk_1330 = buffer.data(slk + 1330);
    const auto *slk_1331 = buffer.data(slk + 1331);
    const auto *slk_1337 = buffer.data(slk + 1337);
    const auto *slk_1341 = buffer.data(slk + 1341);
    const auto *slk_1344 = buffer.data(slk + 1344);
    const auto *slk_1346 = buffer.data(slk + 1346);
    const auto *slk_1349 = buffer.data(slk + 1349);
    const auto *slk_1350 = buffer.data(slk + 1350);
    const auto *slk_1352 = buffer.data(slk + 1352);
    const auto *slk_1355 = buffer.data(slk + 1355);
    const auto *slk_1356 = buffer.data(slk + 1356);
    const auto *slk_1357 = buffer.data(slk + 1357);
    const auto *slk_1359 = buffer.data(slk + 1359);
    const auto *slk_1360 = buffer.data(slk + 1360);
    const auto *slk_1361 = buffer.data(slk + 1361);
    const auto *slk_1362 = buffer.data(slk + 1362);
    const auto *slk_1363 = buffer.data(slk + 1363);
    const auto *slk_1364 = buffer.data(slk + 1364);
    const auto *slk_1365 = buffer.data(slk + 1365);
    const auto *slk_1366 = buffer.data(slk + 1366);
    const auto *slk_1367 = buffer.data(slk + 1367);
    const auto *slk_1368 = buffer.data(slk + 1368);
    const auto *slk_1371 = buffer.data(slk + 1371);
    const auto *slk_1373 = buffer.data(slk + 1373);
    const auto *slk_1374 = buffer.data(slk + 1374);
    const auto *slk_1377 = buffer.data(slk + 1377);
    const auto *slk_1378 = buffer.data(slk + 1378);
    const auto *slk_1380 = buffer.data(slk + 1380);
    const auto *slk_1382 = buffer.data(slk + 1382);
    const auto *slk_1383 = buffer.data(slk + 1383);
    const auto *slk_1385 = buffer.data(slk + 1385);
    const auto *slk_1386 = buffer.data(slk + 1386);
    const auto *slk_1388 = buffer.data(slk + 1388);
    const auto *slk_1389 = buffer.data(slk + 1389);
    const auto *slk_1391 = buffer.data(slk + 1391);
    const auto *slk_1392 = buffer.data(slk + 1392);
    const auto *slk_1393 = buffer.data(slk + 1393);
    const auto *slk_1395 = buffer.data(slk + 1395);
    const auto *slk_1396 = buffer.data(slk + 1396);
    const auto *slk_1397 = buffer.data(slk + 1397);
    const auto *slk_1398 = buffer.data(slk + 1398);
    const auto *slk_1399 = buffer.data(slk + 1399);
    const auto *slk_1400 = buffer.data(slk + 1400);
    const auto *slk_1401 = buffer.data(slk + 1401);
    const auto *slk_1402 = buffer.data(slk + 1402);
    const auto *slk_1403 = buffer.data(slk + 1403);
    const auto *slk_1404 = buffer.data(slk + 1404);
    const auto *slk_1407 = buffer.data(slk + 1407);
    const auto *slk_1409 = buffer.data(slk + 1409);
    const auto *slk_1410 = buffer.data(slk + 1410);

    const auto *sll1_1260 = buffer.data(sll1 + 1260);
    const auto *sll1_1263 = buffer.data(sll1 + 1263);
    const auto *sll1_1266 = buffer.data(sll1 + 1266);
    const auto *sll1_1270 = buffer.data(sll1 + 1270);
    const auto *sll1_1275 = buffer.data(sll1 + 1275);
    const auto *sll1_1281 = buffer.data(sll1 + 1281);
    const auto *sll1_1643 = buffer.data(sll1 + 1643);
    const auto *sll1_1644 = buffer.data(sll1 + 1644);
    const auto *sll1_1645 = buffer.data(sll1 + 1645);
    const auto *sll1_1647 = buffer.data(sll1 + 1647);
    const auto *sll1_1656 = buffer.data(sll1 + 1656);
    const auto *sll1_1658 = buffer.data(sll1 + 1658);
    const auto *sll1_1659 = buffer.data(sll1 + 1659);
    const auto *sll1_1660 = buffer.data(sll1 + 1660);
    const auto *sll1_1661 = buffer.data(sll1 + 1661);
    const auto *sll1_1662 = buffer.data(sll1 + 1662);
    const auto *sll1_1664 = buffer.data(sll1 + 1664);
    const auto *sll1_1670 = buffer.data(sll1 + 1670);
    const auto *sll1_1674 = buffer.data(sll1 + 1674);
    const auto *sll1_1677 = buffer.data(sll1 + 1677);
    const auto *sll1_1679 = buffer.data(sll1 + 1679);
    const auto *sll1_1682 = buffer.data(sll1 + 1682);
    const auto *sll1_1683 = buffer.data(sll1 + 1683);
    const auto *sll1_1685 = buffer.data(sll1 + 1685);
    const auto *sll1_1688 = buffer.data(sll1 + 1688);
    const auto *sll1_1689 = buffer.data(sll1 + 1689);
    const auto *sll1_1690 = buffer.data(sll1 + 1690);
    const auto *sll1_1692 = buffer.data(sll1 + 1692);
    const auto *sll1_1701 = buffer.data(sll1 + 1701);
    const auto *sll1_1703 = buffer.data(sll1 + 1703);
    const auto *sll1_1704 = buffer.data(sll1 + 1704);
    const auto *sll1_1705 = buffer.data(sll1 + 1705);
    const auto *sll1_1706 = buffer.data(sll1 + 1706);
    const auto *sll1_1707 = buffer.data(sll1 + 1707);
    const auto *sll1_1709 = buffer.data(sll1 + 1709);
    const auto *sll1_1710 = buffer.data(sll1 + 1710);
    const auto *sll1_1713 = buffer.data(sll1 + 1713);
    const auto *sll1_1715 = buffer.data(sll1 + 1715);
    const auto *sll1_1716 = buffer.data(sll1 + 1716);
    const auto *sll1_1719 = buffer.data(sll1 + 1719);
    const auto *sll1_1720 = buffer.data(sll1 + 1720);
    const auto *sll1_1722 = buffer.data(sll1 + 1722);
    const auto *sll1_1724 = buffer.data(sll1 + 1724);
    const auto *sll1_1725 = buffer.data(sll1 + 1725);
    const auto *sll1_1727 = buffer.data(sll1 + 1727);
    const auto *sll1_1728 = buffer.data(sll1 + 1728);
    const auto *sll1_1730 = buffer.data(sll1 + 1730);
    const auto *sll1_1731 = buffer.data(sll1 + 1731);
    const auto *sll1_1733 = buffer.data(sll1 + 1733);
    const auto *sll1_1734 = buffer.data(sll1 + 1734);
    const auto *sll1_1735 = buffer.data(sll1 + 1735);
    const auto *sll1_1737 = buffer.data(sll1 + 1737);
    const auto *sll1_1746 = buffer.data(sll1 + 1746);
    const auto *sll1_1748 = buffer.data(sll1 + 1748);
    const auto *sll1_1749 = buffer.data(sll1 + 1749);
    const auto *sll1_1750 = buffer.data(sll1 + 1750);
    const auto *sll1_1751 = buffer.data(sll1 + 1751);
    const auto *sll1_1752 = buffer.data(sll1 + 1752);
    const auto *sll1_1754 = buffer.data(sll1 + 1754);
    const auto *sll1_1755 = buffer.data(sll1 + 1755);
    const auto *sll1_1758 = buffer.data(sll1 + 1758);
    const auto *sll1_1760 = buffer.data(sll1 + 1760);
    const auto *sll1_1761 = buffer.data(sll1 + 1761);

    const auto *smk_1316 = buffer.data(smk + 1316);
    const auto *smk_1324 = buffer.data(smk + 1324);
    const auto *smk_1325 = buffer.data(smk + 1325);
    const auto *smk_1326 = buffer.data(smk + 1326);
    const auto *smk_1327 = buffer.data(smk + 1327);
    const auto *smk_1328 = buffer.data(smk + 1328);
    const auto *smk_1329 = buffer.data(smk + 1329);
    const auto *smk_1330 = buffer.data(smk + 1330);
    const auto *smk_1331 = buffer.data(smk + 1331);
    const auto *smk_1332 = buffer.data(smk + 1332);
    const auto *smk_1334 = buffer.data(smk + 1334);
    const auto *smk_1335 = buffer.data(smk + 1335);
    const auto *smk_1337 = buffer.data(smk + 1337);
    const auto *smk_1338 = buffer.data(smk + 1338);
    const auto *smk_1341 = buffer.data(smk + 1341);
    const auto *smk_1342 = buffer.data(smk + 1342);
    const auto *smk_1346 = buffer.data(smk + 1346);
    const auto *smk_1347 = buffer.data(smk + 1347);
    const auto *smk_1352 = buffer.data(smk + 1352);
    const auto *smk_1360 = buffer.data(smk + 1360);
    const auto *smk_1361 = buffer.data(smk + 1361);
    const auto *smk_1362 = buffer.data(smk + 1362);
    const auto *smk_1363 = buffer.data(smk + 1363);
    const auto *smk_1364 = buffer.data(smk + 1364);
    const auto *smk_1365 = buffer.data(smk + 1365);
    const auto *smk_1366 = buffer.data(smk + 1366);
    const auto *smk_1367 = buffer.data(smk + 1367);
    const auto *smk_1368 = buffer.data(smk + 1368);
    const auto *smk_1370 = buffer.data(smk + 1370);
    const auto *smk_1371 = buffer.data(smk + 1371);
    const auto *smk_1373 = buffer.data(smk + 1373);
    const auto *smk_1374 = buffer.data(smk + 1374);
    const auto *smk_1377 = buffer.data(smk + 1377);
    const auto *smk_1378 = buffer.data(smk + 1378);
    const auto *smk_1382 = buffer.data(smk + 1382);
    const auto *smk_1383 = buffer.data(smk + 1383);
    const auto *smk_1388 = buffer.data(smk + 1388);
    const auto *smk_1396 = buffer.data(smk + 1396);
    const auto *smk_1397 = buffer.data(smk + 1397);
    const auto *smk_1398 = buffer.data(smk + 1398);
    const auto *smk_1399 = buffer.data(smk + 1399);
    const auto *smk_1400 = buffer.data(smk + 1400);
    const auto *smk_1401 = buffer.data(smk + 1401);
    const auto *smk_1402 = buffer.data(smk + 1402);
    const auto *smk_1403 = buffer.data(smk + 1403);
    const auto *smk_1404 = buffer.data(smk + 1404);
    const auto *smk_1406 = buffer.data(smk + 1406);
    const auto *smk_1407 = buffer.data(smk + 1407);

#pragma omp simd aligned(t_1643, t_1644, t_1645, pb_x, pc_x, sll0_1643, sll0_1644, sll0_1645, \
                         slk_1319, slk_1320, slk_1321, sll1_1643, sll1_1644, \
                         sll1_1645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1643[k] = pb_x[k] * sll0_1643[k]
                    + f_16 * slk_1319[k]
                    - f_14 * pc_x[k] * sll1_1643[k];

        t_1644[k] = pb_x[k] * sll0_1644[k]
                    + f_16 * slk_1320[k]
                    - f_14 * pc_x[k] * sll1_1644[k];

        t_1645[k] = pb_x[k] * sll0_1645[k]
                    + f_16 * slk_1321[k]
                    - f_14 * pc_x[k] * sll1_1645[k];
    }

#pragma omp simd aligned(t_1646, t_1647, t_1648, t_1649, pb_x, pc_x, pc_y, sll0_1647, \
                         slk_1028, slk_1323, slk_1324, slk_1325, sll1_1647, smk_1316, \
                         smk_1324, smk_1325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1646[k] = f_21 * slk_1028[k]
                    + f_3 * pc_y[k] * smk_1316[k];

        t_1647[k] = pb_x[k] * sll0_1647[k]
                    + f_16 * slk_1323[k]
                    - f_14 * pc_x[k] * sll1_1647[k];

        t_1648[k] = f_15 * slk_1324[k]
                    + f_3 * pc_x[k] * smk_1324[k];

        t_1649[k] = f_15 * slk_1325[k]
                    + f_3 * pc_x[k] * smk_1325[k];
    }

#pragma omp simd aligned(t_1650, t_1651, t_1652, t_1653, t_1654, pc_x, slk_1326, slk_1327, \
                         slk_1328, slk_1329, slk_1330, smk_1326, smk_1327, smk_1328, smk_1329, \
                         smk_1330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1650[k] = f_15 * slk_1326[k]
                    + f_3 * pc_x[k] * smk_1326[k];

        t_1651[k] = f_15 * slk_1327[k]
                    + f_3 * pc_x[k] * smk_1327[k];

        t_1652[k] = f_15 * slk_1328[k]
                    + f_3 * pc_x[k] * smk_1328[k];

        t_1653[k] = f_15 * slk_1329[k]
                    + f_3 * pc_x[k] * smk_1329[k];

        t_1654[k] = f_15 * slk_1330[k]
                    + f_3 * pc_x[k] * smk_1330[k];
    }

#pragma omp simd aligned(t_1655, t_1656, t_1657, t_1658, pb_x, pc_x, pc_z, sll0_1656, \
                         sll0_1658, slk_1331, sll1_1656, sll1_1658, smk_1324, \
                         smk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1655[k] = f_15 * slk_1331[k]
                    + f_3 * pc_x[k] * smk_1331[k];

        t_1656[k] = pb_x[k] * sll0_1656[k]
                    - f_14 * pc_x[k] * sll1_1656[k];

        t_1657[k] = f_3 * pc_z[k] * smk_1324[k];

        t_1658[k] = pb_x[k] * sll0_1658[k]
                    - f_14 * pc_x[k] * sll1_1658[k];
    }

#pragma omp simd aligned(t_1659, t_1660, t_1661, t_1662, pb_x, pc_x, sll0_1659, sll0_1660, \
                         sll0_1661, sll0_1662, sll1_1659, sll1_1660, sll1_1661, \
                         sll1_1662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1659[k] = pb_x[k] * sll0_1659[k]
                    - f_14 * pc_x[k] * sll1_1659[k];

        t_1660[k] = pb_x[k] * sll0_1660[k]
                    - f_14 * pc_x[k] * sll1_1660[k];

        t_1661[k] = pb_x[k] * sll0_1661[k]
                    - f_14 * pc_x[k] * sll1_1661[k];

        t_1662[k] = pb_x[k] * sll0_1662[k]
                    - f_14 * pc_x[k] * sll1_1662[k];
    }

#pragma omp simd aligned(t_1663, t_1664, t_1665, pb_x, pb_z, pc_x, pc_y, pc_z, sll0_1260, \
                         sll0_1664, slk_1043, sll1_1260, sll1_1664, \
                         smk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1663[k] = f_21 * slk_1043[k]
                    + f_3 * pc_y[k] * smk_1331[k];

        t_1664[k] = pb_x[k] * sll0_1664[k]
                    - f_14 * pc_x[k] * sll1_1664[k];

        t_1665[k] = pb_z[k] * sll0_1260[k]
                    - f_14 * pc_z[k] * sll1_1260[k];
    }

#pragma omp simd aligned(t_1666, t_1667, t_1668, t_1669, pb_z, pc_y, pc_z, sll0_1263, \
                         slk_1008, slk_1044, slk_1046, sll1_1263, smk_1332, \
                         smk_1334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1666[k] = f_22 * slk_1044[k]
                    + f_3 * pc_y[k] * smk_1332[k];

        t_1667[k] = f_15 * slk_1008[k]
                    + f_3 * pc_z[k] * smk_1332[k];

        t_1668[k] = pb_z[k] * sll0_1263[k]
                    - f_14 * pc_z[k] * sll1_1263[k];

        t_1669[k] = f_22 * slk_1046[k]
                    + f_3 * pc_y[k] * smk_1334[k];
    }

#pragma omp simd aligned(t_1670, t_1671, t_1672, pb_x, pb_z, pc_x, pc_z, sll0_1266, sll0_1670, \
                         slk_1011, slk_1337, sll1_1266, sll1_1670, \
                         smk_1335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1670[k] = pb_x[k] * sll0_1670[k]
                    + f_20 * slk_1337[k]
                    - f_14 * pc_x[k] * sll1_1670[k];

        t_1671[k] = pb_z[k] * sll0_1266[k]
                    - f_14 * pc_z[k] * sll1_1266[k];

        t_1672[k] = f_15 * slk_1011[k]
                    + f_3 * pc_z[k] * smk_1335[k];
    }

#pragma omp simd aligned(t_1673, t_1674, t_1675, pb_x, pb_z, pc_x, pc_y, pc_z, sll0_1270, \
                         sll0_1674, slk_1049, slk_1341, sll1_1270, sll1_1674, \
                         smk_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1673[k] = f_22 * slk_1049[k]
                    + f_3 * pc_y[k] * smk_1337[k];

        t_1674[k] = pb_x[k] * sll0_1674[k]
                    + f_19 * slk_1341[k]
                    - f_14 * pc_x[k] * sll1_1674[k];

        t_1675[k] = pb_z[k] * sll0_1270[k]
                    - f_14 * pc_z[k] * sll1_1270[k];
    }

#pragma omp simd aligned(t_1676, t_1677, t_1678, pb_x, pc_x, pc_y, pc_z, sll0_1677, slk_1014, \
                         slk_1053, slk_1344, sll1_1677, smk_1338, \
                         smk_1341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1676[k] = f_15 * slk_1014[k]
                    + f_3 * pc_z[k] * smk_1338[k];

        t_1677[k] = pb_x[k] * sll0_1677[k]
                    + f_18 * slk_1344[k]
                    - f_14 * pc_x[k] * sll1_1677[k];

        t_1678[k] = f_22 * slk_1053[k]
                    + f_3 * pc_y[k] * smk_1341[k];
    }

#pragma omp simd aligned(t_1679, t_1680, t_1681, pb_x, pb_z, pc_x, pc_z, sll0_1275, sll0_1679, \
                         slk_1018, slk_1346, sll1_1275, sll1_1679, \
                         smk_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1679[k] = pb_x[k] * sll0_1679[k]
                    + f_18 * slk_1346[k]
                    - f_14 * pc_x[k] * sll1_1679[k];

        t_1680[k] = pb_z[k] * sll0_1275[k]
                    - f_14 * pc_z[k] * sll1_1275[k];

        t_1681[k] = f_15 * slk_1018[k]
                    + f_3 * pc_z[k] * smk_1342[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, pb_x, pc_x, pc_y, sll0_1682, sll0_1683, \
                         slk_1058, slk_1349, slk_1350, sll1_1682, sll1_1683, \
                         smk_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = pb_x[k] * sll0_1682[k]
                    + f_17 * slk_1349[k]
                    - f_14 * pc_x[k] * sll1_1682[k];

        t_1683[k] = pb_x[k] * sll0_1683[k]
                    + f_17 * slk_1350[k]
                    - f_14 * pc_x[k] * sll1_1683[k];

        t_1684[k] = f_22 * slk_1058[k]
                    + f_3 * pc_y[k] * smk_1346[k];
    }

#pragma omp simd aligned(t_1685, t_1686, t_1687, pb_x, pb_z, pc_x, pc_z, sll0_1281, sll0_1685, \
                         slk_1023, slk_1352, sll1_1281, sll1_1685, \
                         smk_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1685[k] = pb_x[k] * sll0_1685[k]
                    + f_17 * slk_1352[k]
                    - f_14 * pc_x[k] * sll1_1685[k];

        t_1686[k] = pb_z[k] * sll0_1281[k]
                    - f_14 * pc_z[k] * sll1_1281[k];

        t_1687[k] = f_15 * slk_1023[k]
                    + f_3 * pc_z[k] * smk_1347[k];
    }

#pragma omp simd aligned(t_1688, t_1689, t_1690, pb_x, pc_x, sll0_1688, sll0_1689, sll0_1690, \
                         slk_1355, slk_1356, slk_1357, sll1_1688, sll1_1689, \
                         sll1_1690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1688[k] = pb_x[k] * sll0_1688[k]
                    + f_16 * slk_1355[k]
                    - f_14 * pc_x[k] * sll1_1688[k];

        t_1689[k] = pb_x[k] * sll0_1689[k]
                    + f_16 * slk_1356[k]
                    - f_14 * pc_x[k] * sll1_1689[k];

        t_1690[k] = pb_x[k] * sll0_1690[k]
                    + f_16 * slk_1357[k]
                    - f_14 * pc_x[k] * sll1_1690[k];
    }

#pragma omp simd aligned(t_1691, t_1692, t_1693, t_1694, pb_x, pc_x, pc_y, sll0_1692, \
                         slk_1064, slk_1359, slk_1360, slk_1361, sll1_1692, smk_1352, \
                         smk_1360, smk_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1691[k] = f_22 * slk_1064[k]
                    + f_3 * pc_y[k] * smk_1352[k];

        t_1692[k] = pb_x[k] * sll0_1692[k]
                    + f_16 * slk_1359[k]
                    - f_14 * pc_x[k] * sll1_1692[k];

        t_1693[k] = f_15 * slk_1360[k]
                    + f_3 * pc_x[k] * smk_1360[k];

        t_1694[k] = f_15 * slk_1361[k]
                    + f_3 * pc_x[k] * smk_1361[k];
    }

#pragma omp simd aligned(t_1695, t_1696, t_1697, t_1698, t_1699, pc_x, slk_1362, slk_1363, \
                         slk_1364, slk_1365, slk_1366, smk_1362, smk_1363, smk_1364, smk_1365, \
                         smk_1366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1695[k] = f_15 * slk_1362[k]
                    + f_3 * pc_x[k] * smk_1362[k];

        t_1696[k] = f_15 * slk_1363[k]
                    + f_3 * pc_x[k] * smk_1363[k];

        t_1697[k] = f_15 * slk_1364[k]
                    + f_3 * pc_x[k] * smk_1364[k];

        t_1698[k] = f_15 * slk_1365[k]
                    + f_3 * pc_x[k] * smk_1365[k];

        t_1699[k] = f_15 * slk_1366[k]
                    + f_3 * pc_x[k] * smk_1366[k];
    }

#pragma omp simd aligned(t_1700, t_1701, t_1702, t_1703, pb_x, pc_x, pc_z, sll0_1701, \
                         sll0_1703, slk_1036, slk_1367, sll1_1701, sll1_1703, smk_1360, \
                         smk_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1700[k] = f_15 * slk_1367[k]
                    + f_3 * pc_x[k] * smk_1367[k];

        t_1701[k] = pb_x[k] * sll0_1701[k]
                    - f_14 * pc_x[k] * sll1_1701[k];

        t_1702[k] = f_15 * slk_1036[k]
                    + f_3 * pc_z[k] * smk_1360[k];

        t_1703[k] = pb_x[k] * sll0_1703[k]
                    - f_14 * pc_x[k] * sll1_1703[k];
    }

#pragma omp simd aligned(t_1704, t_1705, t_1706, t_1707, pb_x, pc_x, sll0_1704, sll0_1705, \
                         sll0_1706, sll0_1707, sll1_1704, sll1_1705, sll1_1706, \
                         sll1_1707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1704[k] = pb_x[k] * sll0_1704[k]
                    - f_14 * pc_x[k] * sll1_1704[k];

        t_1705[k] = pb_x[k] * sll0_1705[k]
                    - f_14 * pc_x[k] * sll1_1705[k];

        t_1706[k] = pb_x[k] * sll0_1706[k]
                    - f_14 * pc_x[k] * sll1_1706[k];

        t_1707[k] = pb_x[k] * sll0_1707[k]
                    - f_14 * pc_x[k] * sll1_1707[k];
    }

#pragma omp simd aligned(t_1708, t_1709, t_1710, t_1711, pb_x, pc_x, pc_y, sll0_1709, \
                         sll0_1710, slk_1079, slk_1080, slk_1368, sll1_1709, sll1_1710, \
                         smk_1367, smk_1368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1708[k] = f_22 * slk_1079[k]
                    + f_3 * pc_y[k] * smk_1367[k];

        t_1709[k] = pb_x[k] * sll0_1709[k]
                    - f_14 * pc_x[k] * sll1_1709[k];

        t_1710[k] = pb_x[k] * sll0_1710[k]
                    + f_21 * slk_1368[k]
                    - f_14 * pc_x[k] * sll1_1710[k];

        t_1711[k] = f_20 * slk_1080[k]
                    + f_3 * pc_y[k] * smk_1368[k];
    }

#pragma omp simd aligned(t_1712, t_1713, t_1714, pb_x, pc_x, pc_y, pc_z, sll0_1713, slk_1044, \
                         slk_1082, slk_1371, sll1_1713, smk_1368, \
                         smk_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1712[k] = f_16 * slk_1044[k]
                    + f_3 * pc_z[k] * smk_1368[k];

        t_1713[k] = pb_x[k] * sll0_1713[k]
                    + f_20 * slk_1371[k]
                    - f_14 * pc_x[k] * sll1_1713[k];

        t_1714[k] = f_20 * slk_1082[k]
                    + f_3 * pc_y[k] * smk_1370[k];
    }

#pragma omp simd aligned(t_1715, t_1716, t_1717, pb_x, pc_x, pc_z, sll0_1715, sll0_1716, \
                         slk_1047, slk_1373, slk_1374, sll1_1715, sll1_1716, \
                         smk_1371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1715[k] = pb_x[k] * sll0_1715[k]
                    + f_20 * slk_1373[k]
                    - f_14 * pc_x[k] * sll1_1715[k];

        t_1716[k] = pb_x[k] * sll0_1716[k]
                    + f_19 * slk_1374[k]
                    - f_14 * pc_x[k] * sll1_1716[k];

        t_1717[k] = f_16 * slk_1047[k]
                    + f_3 * pc_z[k] * smk_1371[k];
    }

#pragma omp simd aligned(t_1718, t_1719, t_1720, pb_x, pc_x, pc_y, sll0_1719, sll0_1720, \
                         slk_1085, slk_1377, slk_1378, sll1_1719, sll1_1720, \
                         smk_1373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1718[k] = f_20 * slk_1085[k]
                    + f_3 * pc_y[k] * smk_1373[k];

        t_1719[k] = pb_x[k] * sll0_1719[k]
                    + f_19 * slk_1377[k]
                    - f_14 * pc_x[k] * sll1_1719[k];

        t_1720[k] = pb_x[k] * sll0_1720[k]
                    + f_18 * slk_1378[k]
                    - f_14 * pc_x[k] * sll1_1720[k];
    }

#pragma omp simd aligned(t_1721, t_1722, t_1723, pb_x, pc_x, pc_y, pc_z, sll0_1722, slk_1050, \
                         slk_1089, slk_1380, sll1_1722, smk_1374, \
                         smk_1377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1721[k] = f_16 * slk_1050[k]
                    + f_3 * pc_z[k] * smk_1374[k];

        t_1722[k] = pb_x[k] * sll0_1722[k]
                    + f_18 * slk_1380[k]
                    - f_14 * pc_x[k] * sll1_1722[k];

        t_1723[k] = f_20 * slk_1089[k]
                    + f_3 * pc_y[k] * smk_1377[k];
    }

#pragma omp simd aligned(t_1724, t_1725, t_1726, pb_x, pc_x, pc_z, sll0_1724, sll0_1725, \
                         slk_1054, slk_1382, slk_1383, sll1_1724, sll1_1725, \
                         smk_1378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1724[k] = pb_x[k] * sll0_1724[k]
                    + f_18 * slk_1382[k]
                    - f_14 * pc_x[k] * sll1_1724[k];

        t_1725[k] = pb_x[k] * sll0_1725[k]
                    + f_17 * slk_1383[k]
                    - f_14 * pc_x[k] * sll1_1725[k];

        t_1726[k] = f_16 * slk_1054[k]
                    + f_3 * pc_z[k] * smk_1378[k];
    }

#pragma omp simd aligned(t_1727, t_1728, t_1729, pb_x, pc_x, pc_y, sll0_1727, sll0_1728, \
                         slk_1094, slk_1385, slk_1386, sll1_1727, sll1_1728, \
                         smk_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1727[k] = pb_x[k] * sll0_1727[k]
                    + f_17 * slk_1385[k]
                    - f_14 * pc_x[k] * sll1_1727[k];

        t_1728[k] = pb_x[k] * sll0_1728[k]
                    + f_17 * slk_1386[k]
                    - f_14 * pc_x[k] * sll1_1728[k];

        t_1729[k] = f_20 * slk_1094[k]
                    + f_3 * pc_y[k] * smk_1382[k];
    }

#pragma omp simd aligned(t_1730, t_1731, t_1732, pb_x, pc_x, pc_z, sll0_1730, sll0_1731, \
                         slk_1059, slk_1388, slk_1389, sll1_1730, sll1_1731, \
                         smk_1383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1730[k] = pb_x[k] * sll0_1730[k]
                    + f_17 * slk_1388[k]
                    - f_14 * pc_x[k] * sll1_1730[k];

        t_1731[k] = pb_x[k] * sll0_1731[k]
                    + f_16 * slk_1389[k]
                    - f_14 * pc_x[k] * sll1_1731[k];

        t_1732[k] = f_16 * slk_1059[k]
                    + f_3 * pc_z[k] * smk_1383[k];
    }

#pragma omp simd aligned(t_1733, t_1734, t_1735, pb_x, pc_x, sll0_1733, sll0_1734, sll0_1735, \
                         slk_1391, slk_1392, slk_1393, sll1_1733, sll1_1734, \
                         sll1_1735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1733[k] = pb_x[k] * sll0_1733[k]
                    + f_16 * slk_1391[k]
                    - f_14 * pc_x[k] * sll1_1733[k];

        t_1734[k] = pb_x[k] * sll0_1734[k]
                    + f_16 * slk_1392[k]
                    - f_14 * pc_x[k] * sll1_1734[k];

        t_1735[k] = pb_x[k] * sll0_1735[k]
                    + f_16 * slk_1393[k]
                    - f_14 * pc_x[k] * sll1_1735[k];
    }

#pragma omp simd aligned(t_1736, t_1737, t_1738, t_1739, pb_x, pc_x, pc_y, sll0_1737, \
                         slk_1100, slk_1395, slk_1396, slk_1397, sll1_1737, smk_1388, \
                         smk_1396, smk_1397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1736[k] = f_20 * slk_1100[k]
                    + f_3 * pc_y[k] * smk_1388[k];

        t_1737[k] = pb_x[k] * sll0_1737[k]
                    + f_16 * slk_1395[k]
                    - f_14 * pc_x[k] * sll1_1737[k];

        t_1738[k] = f_15 * slk_1396[k]
                    + f_3 * pc_x[k] * smk_1396[k];

        t_1739[k] = f_15 * slk_1397[k]
                    + f_3 * pc_x[k] * smk_1397[k];
    }

#pragma omp simd aligned(t_1740, t_1741, t_1742, t_1743, t_1744, pc_x, slk_1398, slk_1399, \
                         slk_1400, slk_1401, slk_1402, smk_1398, smk_1399, smk_1400, smk_1401, \
                         smk_1402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1740[k] = f_15 * slk_1398[k]
                    + f_3 * pc_x[k] * smk_1398[k];

        t_1741[k] = f_15 * slk_1399[k]
                    + f_3 * pc_x[k] * smk_1399[k];

        t_1742[k] = f_15 * slk_1400[k]
                    + f_3 * pc_x[k] * smk_1400[k];

        t_1743[k] = f_15 * slk_1401[k]
                    + f_3 * pc_x[k] * smk_1401[k];

        t_1744[k] = f_15 * slk_1402[k]
                    + f_3 * pc_x[k] * smk_1402[k];
    }

#pragma omp simd aligned(t_1745, t_1746, t_1747, t_1748, pb_x, pc_x, pc_z, sll0_1746, \
                         sll0_1748, slk_1072, slk_1403, sll1_1746, sll1_1748, smk_1396, \
                         smk_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1745[k] = f_15 * slk_1403[k]
                    + f_3 * pc_x[k] * smk_1403[k];

        t_1746[k] = pb_x[k] * sll0_1746[k]
                    - f_14 * pc_x[k] * sll1_1746[k];

        t_1747[k] = f_16 * slk_1072[k]
                    + f_3 * pc_z[k] * smk_1396[k];

        t_1748[k] = pb_x[k] * sll0_1748[k]
                    - f_14 * pc_x[k] * sll1_1748[k];
    }

#pragma omp simd aligned(t_1749, t_1750, t_1751, t_1752, pb_x, pc_x, sll0_1749, sll0_1750, \
                         sll0_1751, sll0_1752, sll1_1749, sll1_1750, sll1_1751, \
                         sll1_1752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1749[k] = pb_x[k] * sll0_1749[k]
                    - f_14 * pc_x[k] * sll1_1749[k];

        t_1750[k] = pb_x[k] * sll0_1750[k]
                    - f_14 * pc_x[k] * sll1_1750[k];

        t_1751[k] = pb_x[k] * sll0_1751[k]
                    - f_14 * pc_x[k] * sll1_1751[k];

        t_1752[k] = pb_x[k] * sll0_1752[k]
                    - f_14 * pc_x[k] * sll1_1752[k];
    }

#pragma omp simd aligned(t_1753, t_1754, t_1755, t_1756, pb_x, pc_x, pc_y, sll0_1754, \
                         sll0_1755, slk_1115, slk_1116, slk_1404, sll1_1754, sll1_1755, \
                         smk_1403, smk_1404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1753[k] = f_20 * slk_1115[k]
                    + f_3 * pc_y[k] * smk_1403[k];

        t_1754[k] = pb_x[k] * sll0_1754[k]
                    - f_14 * pc_x[k] * sll1_1754[k];

        t_1755[k] = pb_x[k] * sll0_1755[k]
                    + f_21 * slk_1404[k]
                    - f_14 * pc_x[k] * sll1_1755[k];

        t_1756[k] = f_19 * slk_1116[k]
                    + f_3 * pc_y[k] * smk_1404[k];
    }

#pragma omp simd aligned(t_1757, t_1758, t_1759, pb_x, pc_x, pc_y, pc_z, sll0_1758, slk_1080, \
                         slk_1118, slk_1407, sll1_1758, smk_1404, \
                         smk_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1757[k] = f_17 * slk_1080[k]
                    + f_3 * pc_z[k] * smk_1404[k];

        t_1758[k] = pb_x[k] * sll0_1758[k]
                    + f_20 * slk_1407[k]
                    - f_14 * pc_x[k] * sll1_1758[k];

        t_1759[k] = f_19 * slk_1118[k]
                    + f_3 * pc_y[k] * smk_1406[k];
    }

#pragma omp simd aligned(t_1760, t_1761, t_1762, pb_x, pc_x, pc_z, sll0_1760, sll0_1761, \
                         slk_1083, slk_1409, slk_1410, sll1_1760, sll1_1761, \
                         smk_1407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1760[k] = pb_x[k] * sll0_1760[k]
                    + f_20 * slk_1409[k]
                    - f_14 * pc_x[k] * sll1_1760[k];

        t_1761[k] = pb_x[k] * sll0_1761[k]
                    + f_19 * slk_1410[k]
                    - f_14 * pc_x[k] * sll1_1761[k];

        t_1762[k] = f_17 * slk_1083[k]
                    + f_3 * pc_z[k] * smk_1407[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_1764 = buffer.data(sll0 + 1764);
    const auto *sll0_1765 = buffer.data(sll0 + 1765);
    const auto *sll0_1767 = buffer.data(sll0 + 1767);
    const auto *sll0_1769 = buffer.data(sll0 + 1769);
    const auto *sll0_1770 = buffer.data(sll0 + 1770);
    const auto *sll0_1772 = buffer.data(sll0 + 1772);
    const auto *sll0_1773 = buffer.data(sll0 + 1773);
    const auto *sll0_1775 = buffer.data(sll0 + 1775);
    const auto *sll0_1776 = buffer.data(sll0 + 1776);
    const auto *sll0_1778 = buffer.data(sll0 + 1778);
    const auto *sll0_1779 = buffer.data(sll0 + 1779);
    const auto *sll0_1780 = buffer.data(sll0 + 1780);
    const auto *sll0_1782 = buffer.data(sll0 + 1782);
    const auto *sll0_1791 = buffer.data(sll0 + 1791);
    const auto *sll0_1793 = buffer.data(sll0 + 1793);
    const auto *sll0_1794 = buffer.data(sll0 + 1794);
    const auto *sll0_1795 = buffer.data(sll0 + 1795);
    const auto *sll0_1796 = buffer.data(sll0 + 1796);
    const auto *sll0_1797 = buffer.data(sll0 + 1797);
    const auto *sll0_1799 = buffer.data(sll0 + 1799);
    const auto *sll0_1800 = buffer.data(sll0 + 1800);
    const auto *sll0_1803 = buffer.data(sll0 + 1803);
    const auto *sll0_1805 = buffer.data(sll0 + 1805);
    const auto *sll0_1806 = buffer.data(sll0 + 1806);
    const auto *sll0_1809 = buffer.data(sll0 + 1809);
    const auto *sll0_1810 = buffer.data(sll0 + 1810);
    const auto *sll0_1812 = buffer.data(sll0 + 1812);
    const auto *sll0_1814 = buffer.data(sll0 + 1814);
    const auto *sll0_1815 = buffer.data(sll0 + 1815);
    const auto *sll0_1817 = buffer.data(sll0 + 1817);
    const auto *sll0_1818 = buffer.data(sll0 + 1818);
    const auto *sll0_1820 = buffer.data(sll0 + 1820);
    const auto *sll0_1821 = buffer.data(sll0 + 1821);
    const auto *sll0_1823 = buffer.data(sll0 + 1823);
    const auto *sll0_1824 = buffer.data(sll0 + 1824);
    const auto *sll0_1825 = buffer.data(sll0 + 1825);
    const auto *sll0_1827 = buffer.data(sll0 + 1827);
    const auto *sll0_1836 = buffer.data(sll0 + 1836);
    const auto *sll0_1838 = buffer.data(sll0 + 1838);
    const auto *sll0_1839 = buffer.data(sll0 + 1839);
    const auto *sll0_1840 = buffer.data(sll0 + 1840);
    const auto *sll0_1841 = buffer.data(sll0 + 1841);
    const auto *sll0_1842 = buffer.data(sll0 + 1842);
    const auto *sll0_1844 = buffer.data(sll0 + 1844);
    const auto *sll0_1845 = buffer.data(sll0 + 1845);
    const auto *sll0_1848 = buffer.data(sll0 + 1848);
    const auto *sll0_1850 = buffer.data(sll0 + 1850);
    const auto *sll0_1851 = buffer.data(sll0 + 1851);
    const auto *sll0_1854 = buffer.data(sll0 + 1854);
    const auto *sll0_1855 = buffer.data(sll0 + 1855);
    const auto *sll0_1857 = buffer.data(sll0 + 1857);
    const auto *sll0_1859 = buffer.data(sll0 + 1859);
    const auto *sll0_1860 = buffer.data(sll0 + 1860);
    const auto *sll0_1862 = buffer.data(sll0 + 1862);
    const auto *sll0_1863 = buffer.data(sll0 + 1863);
    const auto *sll0_1865 = buffer.data(sll0 + 1865);
    const auto *sll0_1866 = buffer.data(sll0 + 1866);
    const auto *sll0_1868 = buffer.data(sll0 + 1868);
    const auto *sll0_1869 = buffer.data(sll0 + 1869);
    const auto *sll0_1870 = buffer.data(sll0 + 1870);
    const auto *sll0_1872 = buffer.data(sll0 + 1872);

    const auto *slk_1086 = buffer.data(slk + 1086);
    const auto *slk_1090 = buffer.data(slk + 1090);
    const auto *slk_1095 = buffer.data(slk + 1095);
    const auto *slk_1108 = buffer.data(slk + 1108);
    const auto *slk_1116 = buffer.data(slk + 1116);
    const auto *slk_1119 = buffer.data(slk + 1119);
    const auto *slk_1121 = buffer.data(slk + 1121);
    const auto *slk_1122 = buffer.data(slk + 1122);
    const auto *slk_1125 = buffer.data(slk + 1125);
    const auto *slk_1126 = buffer.data(slk + 1126);
    const auto *slk_1130 = buffer.data(slk + 1130);
    const auto *slk_1131 = buffer.data(slk + 1131);
    const auto *slk_1136 = buffer.data(slk + 1136);
    const auto *slk_1144 = buffer.data(slk + 1144);
    const auto *slk_1151 = buffer.data(slk + 1151);
    const auto *slk_1152 = buffer.data(slk + 1152);
    const auto *slk_1154 = buffer.data(slk + 1154);
    const auto *slk_1155 = buffer.data(slk + 1155);
    const auto *slk_1157 = buffer.data(slk + 1157);
    const auto *slk_1158 = buffer.data(slk + 1158);
    const auto *slk_1161 = buffer.data(slk + 1161);
    const auto *slk_1162 = buffer.data(slk + 1162);
    const auto *slk_1166 = buffer.data(slk + 1166);
    const auto *slk_1167 = buffer.data(slk + 1167);
    const auto *slk_1172 = buffer.data(slk + 1172);
    const auto *slk_1187 = buffer.data(slk + 1187);
    const auto *slk_1188 = buffer.data(slk + 1188);
    const auto *slk_1190 = buffer.data(slk + 1190);
    const auto *slk_1193 = buffer.data(slk + 1193);
    const auto *slk_1197 = buffer.data(slk + 1197);
    const auto *slk_1202 = buffer.data(slk + 1202);
    const auto *slk_1208 = buffer.data(slk + 1208);
    const auto *slk_1413 = buffer.data(slk + 1413);
    const auto *slk_1414 = buffer.data(slk + 1414);
    const auto *slk_1416 = buffer.data(slk + 1416);
    const auto *slk_1418 = buffer.data(slk + 1418);
    const auto *slk_1419 = buffer.data(slk + 1419);
    const auto *slk_1421 = buffer.data(slk + 1421);
    const auto *slk_1422 = buffer.data(slk + 1422);
    const auto *slk_1424 = buffer.data(slk + 1424);
    const auto *slk_1425 = buffer.data(slk + 1425);
    const auto *slk_1427 = buffer.data(slk + 1427);
    const auto *slk_1428 = buffer.data(slk + 1428);
    const auto *slk_1429 = buffer.data(slk + 1429);
    const auto *slk_1431 = buffer.data(slk + 1431);
    const auto *slk_1432 = buffer.data(slk + 1432);
    const auto *slk_1433 = buffer.data(slk + 1433);
    const auto *slk_1434 = buffer.data(slk + 1434);
    const auto *slk_1435 = buffer.data(slk + 1435);
    const auto *slk_1436 = buffer.data(slk + 1436);
    const auto *slk_1437 = buffer.data(slk + 1437);
    const auto *slk_1438 = buffer.data(slk + 1438);
    const auto *slk_1439 = buffer.data(slk + 1439);
    const auto *slk_1440 = buffer.data(slk + 1440);
    const auto *slk_1443 = buffer.data(slk + 1443);
    const auto *slk_1445 = buffer.data(slk + 1445);
    const auto *slk_1446 = buffer.data(slk + 1446);
    const auto *slk_1449 = buffer.data(slk + 1449);
    const auto *slk_1450 = buffer.data(slk + 1450);
    const auto *slk_1452 = buffer.data(slk + 1452);
    const auto *slk_1454 = buffer.data(slk + 1454);
    const auto *slk_1455 = buffer.data(slk + 1455);
    const auto *slk_1457 = buffer.data(slk + 1457);
    const auto *slk_1458 = buffer.data(slk + 1458);
    const auto *slk_1460 = buffer.data(slk + 1460);
    const auto *slk_1461 = buffer.data(slk + 1461);
    const auto *slk_1463 = buffer.data(slk + 1463);
    const auto *slk_1464 = buffer.data(slk + 1464);
    const auto *slk_1465 = buffer.data(slk + 1465);
    const auto *slk_1467 = buffer.data(slk + 1467);
    const auto *slk_1468 = buffer.data(slk + 1468);
    const auto *slk_1469 = buffer.data(slk + 1469);
    const auto *slk_1470 = buffer.data(slk + 1470);
    const auto *slk_1471 = buffer.data(slk + 1471);
    const auto *slk_1472 = buffer.data(slk + 1472);
    const auto *slk_1473 = buffer.data(slk + 1473);
    const auto *slk_1474 = buffer.data(slk + 1474);
    const auto *slk_1475 = buffer.data(slk + 1475);
    const auto *slk_1476 = buffer.data(slk + 1476);
    const auto *slk_1479 = buffer.data(slk + 1479);
    const auto *slk_1481 = buffer.data(slk + 1481);
    const auto *slk_1482 = buffer.data(slk + 1482);
    const auto *slk_1485 = buffer.data(slk + 1485);
    const auto *slk_1486 = buffer.data(slk + 1486);
    const auto *slk_1488 = buffer.data(slk + 1488);
    const auto *slk_1490 = buffer.data(slk + 1490);
    const auto *slk_1491 = buffer.data(slk + 1491);
    const auto *slk_1493 = buffer.data(slk + 1493);
    const auto *slk_1494 = buffer.data(slk + 1494);
    const auto *slk_1496 = buffer.data(slk + 1496);
    const auto *slk_1497 = buffer.data(slk + 1497);
    const auto *slk_1499 = buffer.data(slk + 1499);
    const auto *slk_1500 = buffer.data(slk + 1500);
    const auto *slk_1501 = buffer.data(slk + 1501);
    const auto *slk_1503 = buffer.data(slk + 1503);
    const auto *slk_1504 = buffer.data(slk + 1504);
    const auto *slk_1505 = buffer.data(slk + 1505);
    const auto *slk_1506 = buffer.data(slk + 1506);
    const auto *slk_1507 = buffer.data(slk + 1507);
    const auto *slk_1508 = buffer.data(slk + 1508);
    const auto *slk_1509 = buffer.data(slk + 1509);
    const auto *slk_1510 = buffer.data(slk + 1510);

    const auto *sll1_1764 = buffer.data(sll1 + 1764);
    const auto *sll1_1765 = buffer.data(sll1 + 1765);
    const auto *sll1_1767 = buffer.data(sll1 + 1767);
    const auto *sll1_1769 = buffer.data(sll1 + 1769);
    const auto *sll1_1770 = buffer.data(sll1 + 1770);
    const auto *sll1_1772 = buffer.data(sll1 + 1772);
    const auto *sll1_1773 = buffer.data(sll1 + 1773);
    const auto *sll1_1775 = buffer.data(sll1 + 1775);
    const auto *sll1_1776 = buffer.data(sll1 + 1776);
    const auto *sll1_1778 = buffer.data(sll1 + 1778);
    const auto *sll1_1779 = buffer.data(sll1 + 1779);
    const auto *sll1_1780 = buffer.data(sll1 + 1780);
    const auto *sll1_1782 = buffer.data(sll1 + 1782);
    const auto *sll1_1791 = buffer.data(sll1 + 1791);
    const auto *sll1_1793 = buffer.data(sll1 + 1793);
    const auto *sll1_1794 = buffer.data(sll1 + 1794);
    const auto *sll1_1795 = buffer.data(sll1 + 1795);
    const auto *sll1_1796 = buffer.data(sll1 + 1796);
    const auto *sll1_1797 = buffer.data(sll1 + 1797);
    const auto *sll1_1799 = buffer.data(sll1 + 1799);
    const auto *sll1_1800 = buffer.data(sll1 + 1800);
    const auto *sll1_1803 = buffer.data(sll1 + 1803);
    const auto *sll1_1805 = buffer.data(sll1 + 1805);
    const auto *sll1_1806 = buffer.data(sll1 + 1806);
    const auto *sll1_1809 = buffer.data(sll1 + 1809);
    const auto *sll1_1810 = buffer.data(sll1 + 1810);
    const auto *sll1_1812 = buffer.data(sll1 + 1812);
    const auto *sll1_1814 = buffer.data(sll1 + 1814);
    const auto *sll1_1815 = buffer.data(sll1 + 1815);
    const auto *sll1_1817 = buffer.data(sll1 + 1817);
    const auto *sll1_1818 = buffer.data(sll1 + 1818);
    const auto *sll1_1820 = buffer.data(sll1 + 1820);
    const auto *sll1_1821 = buffer.data(sll1 + 1821);
    const auto *sll1_1823 = buffer.data(sll1 + 1823);
    const auto *sll1_1824 = buffer.data(sll1 + 1824);
    const auto *sll1_1825 = buffer.data(sll1 + 1825);
    const auto *sll1_1827 = buffer.data(sll1 + 1827);
    const auto *sll1_1836 = buffer.data(sll1 + 1836);
    const auto *sll1_1838 = buffer.data(sll1 + 1838);
    const auto *sll1_1839 = buffer.data(sll1 + 1839);
    const auto *sll1_1840 = buffer.data(sll1 + 1840);
    const auto *sll1_1841 = buffer.data(sll1 + 1841);
    const auto *sll1_1842 = buffer.data(sll1 + 1842);
    const auto *sll1_1844 = buffer.data(sll1 + 1844);
    const auto *sll1_1845 = buffer.data(sll1 + 1845);
    const auto *sll1_1848 = buffer.data(sll1 + 1848);
    const auto *sll1_1850 = buffer.data(sll1 + 1850);
    const auto *sll1_1851 = buffer.data(sll1 + 1851);
    const auto *sll1_1854 = buffer.data(sll1 + 1854);
    const auto *sll1_1855 = buffer.data(sll1 + 1855);
    const auto *sll1_1857 = buffer.data(sll1 + 1857);
    const auto *sll1_1859 = buffer.data(sll1 + 1859);
    const auto *sll1_1860 = buffer.data(sll1 + 1860);
    const auto *sll1_1862 = buffer.data(sll1 + 1862);
    const auto *sll1_1863 = buffer.data(sll1 + 1863);
    const auto *sll1_1865 = buffer.data(sll1 + 1865);
    const auto *sll1_1866 = buffer.data(sll1 + 1866);
    const auto *sll1_1868 = buffer.data(sll1 + 1868);
    const auto *sll1_1869 = buffer.data(sll1 + 1869);
    const auto *sll1_1870 = buffer.data(sll1 + 1870);
    const auto *sll1_1872 = buffer.data(sll1 + 1872);

    const auto *smk_1409 = buffer.data(smk + 1409);
    const auto *smk_1410 = buffer.data(smk + 1410);
    const auto *smk_1413 = buffer.data(smk + 1413);
    const auto *smk_1414 = buffer.data(smk + 1414);
    const auto *smk_1418 = buffer.data(smk + 1418);
    const auto *smk_1419 = buffer.data(smk + 1419);
    const auto *smk_1424 = buffer.data(smk + 1424);
    const auto *smk_1432 = buffer.data(smk + 1432);
    const auto *smk_1433 = buffer.data(smk + 1433);
    const auto *smk_1434 = buffer.data(smk + 1434);
    const auto *smk_1435 = buffer.data(smk + 1435);
    const auto *smk_1436 = buffer.data(smk + 1436);
    const auto *smk_1437 = buffer.data(smk + 1437);
    const auto *smk_1438 = buffer.data(smk + 1438);
    const auto *smk_1439 = buffer.data(smk + 1439);
    const auto *smk_1440 = buffer.data(smk + 1440);
    const auto *smk_1442 = buffer.data(smk + 1442);
    const auto *smk_1443 = buffer.data(smk + 1443);
    const auto *smk_1445 = buffer.data(smk + 1445);
    const auto *smk_1446 = buffer.data(smk + 1446);
    const auto *smk_1449 = buffer.data(smk + 1449);
    const auto *smk_1450 = buffer.data(smk + 1450);
    const auto *smk_1454 = buffer.data(smk + 1454);
    const auto *smk_1455 = buffer.data(smk + 1455);
    const auto *smk_1460 = buffer.data(smk + 1460);
    const auto *smk_1468 = buffer.data(smk + 1468);
    const auto *smk_1469 = buffer.data(smk + 1469);
    const auto *smk_1470 = buffer.data(smk + 1470);
    const auto *smk_1471 = buffer.data(smk + 1471);
    const auto *smk_1472 = buffer.data(smk + 1472);
    const auto *smk_1473 = buffer.data(smk + 1473);
    const auto *smk_1474 = buffer.data(smk + 1474);
    const auto *smk_1475 = buffer.data(smk + 1475);
    const auto *smk_1476 = buffer.data(smk + 1476);
    const auto *smk_1478 = buffer.data(smk + 1478);
    const auto *smk_1479 = buffer.data(smk + 1479);
    const auto *smk_1481 = buffer.data(smk + 1481);
    const auto *smk_1482 = buffer.data(smk + 1482);
    const auto *smk_1485 = buffer.data(smk + 1485);
    const auto *smk_1486 = buffer.data(smk + 1486);
    const auto *smk_1490 = buffer.data(smk + 1490);
    const auto *smk_1491 = buffer.data(smk + 1491);
    const auto *smk_1496 = buffer.data(smk + 1496);
    const auto *smk_1504 = buffer.data(smk + 1504);
    const auto *smk_1505 = buffer.data(smk + 1505);
    const auto *smk_1506 = buffer.data(smk + 1506);
    const auto *smk_1507 = buffer.data(smk + 1507);
    const auto *smk_1508 = buffer.data(smk + 1508);
    const auto *smk_1509 = buffer.data(smk + 1509);
    const auto *smk_1510 = buffer.data(smk + 1510);

#pragma omp simd aligned(t_1763, t_1764, t_1765, pb_x, pc_x, pc_y, sll0_1764, sll0_1765, \
                         slk_1121, slk_1413, slk_1414, sll1_1764, sll1_1765, \
                         smk_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1763[k] = f_19 * slk_1121[k]
                    + f_3 * pc_y[k] * smk_1409[k];

        t_1764[k] = pb_x[k] * sll0_1764[k]
                    + f_19 * slk_1413[k]
                    - f_14 * pc_x[k] * sll1_1764[k];

        t_1765[k] = pb_x[k] * sll0_1765[k]
                    + f_18 * slk_1414[k]
                    - f_14 * pc_x[k] * sll1_1765[k];
    }

#pragma omp simd aligned(t_1766, t_1767, t_1768, pb_x, pc_x, pc_y, pc_z, sll0_1767, slk_1086, \
                         slk_1125, slk_1416, sll1_1767, smk_1410, \
                         smk_1413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1766[k] = f_17 * slk_1086[k]
                    + f_3 * pc_z[k] * smk_1410[k];

        t_1767[k] = pb_x[k] * sll0_1767[k]
                    + f_18 * slk_1416[k]
                    - f_14 * pc_x[k] * sll1_1767[k];

        t_1768[k] = f_19 * slk_1125[k]
                    + f_3 * pc_y[k] * smk_1413[k];
    }

#pragma omp simd aligned(t_1769, t_1770, t_1771, pb_x, pc_x, pc_z, sll0_1769, sll0_1770, \
                         slk_1090, slk_1418, slk_1419, sll1_1769, sll1_1770, \
                         smk_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1769[k] = pb_x[k] * sll0_1769[k]
                    + f_18 * slk_1418[k]
                    - f_14 * pc_x[k] * sll1_1769[k];

        t_1770[k] = pb_x[k] * sll0_1770[k]
                    + f_17 * slk_1419[k]
                    - f_14 * pc_x[k] * sll1_1770[k];

        t_1771[k] = f_17 * slk_1090[k]
                    + f_3 * pc_z[k] * smk_1414[k];
    }

#pragma omp simd aligned(t_1772, t_1773, t_1774, pb_x, pc_x, pc_y, sll0_1772, sll0_1773, \
                         slk_1130, slk_1421, slk_1422, sll1_1772, sll1_1773, \
                         smk_1418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1772[k] = pb_x[k] * sll0_1772[k]
                    + f_17 * slk_1421[k]
                    - f_14 * pc_x[k] * sll1_1772[k];

        t_1773[k] = pb_x[k] * sll0_1773[k]
                    + f_17 * slk_1422[k]
                    - f_14 * pc_x[k] * sll1_1773[k];

        t_1774[k] = f_19 * slk_1130[k]
                    + f_3 * pc_y[k] * smk_1418[k];
    }

#pragma omp simd aligned(t_1775, t_1776, t_1777, pb_x, pc_x, pc_z, sll0_1775, sll0_1776, \
                         slk_1095, slk_1424, slk_1425, sll1_1775, sll1_1776, \
                         smk_1419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1775[k] = pb_x[k] * sll0_1775[k]
                    + f_17 * slk_1424[k]
                    - f_14 * pc_x[k] * sll1_1775[k];

        t_1776[k] = pb_x[k] * sll0_1776[k]
                    + f_16 * slk_1425[k]
                    - f_14 * pc_x[k] * sll1_1776[k];

        t_1777[k] = f_17 * slk_1095[k]
                    + f_3 * pc_z[k] * smk_1419[k];
    }

#pragma omp simd aligned(t_1778, t_1779, t_1780, pb_x, pc_x, sll0_1778, sll0_1779, sll0_1780, \
                         slk_1427, slk_1428, slk_1429, sll1_1778, sll1_1779, \
                         sll1_1780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1778[k] = pb_x[k] * sll0_1778[k]
                    + f_16 * slk_1427[k]
                    - f_14 * pc_x[k] * sll1_1778[k];

        t_1779[k] = pb_x[k] * sll0_1779[k]
                    + f_16 * slk_1428[k]
                    - f_14 * pc_x[k] * sll1_1779[k];

        t_1780[k] = pb_x[k] * sll0_1780[k]
                    + f_16 * slk_1429[k]
                    - f_14 * pc_x[k] * sll1_1780[k];
    }

#pragma omp simd aligned(t_1781, t_1782, t_1783, t_1784, pb_x, pc_x, pc_y, sll0_1782, \
                         slk_1136, slk_1431, slk_1432, slk_1433, sll1_1782, smk_1424, \
                         smk_1432, smk_1433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1781[k] = f_19 * slk_1136[k]
                    + f_3 * pc_y[k] * smk_1424[k];

        t_1782[k] = pb_x[k] * sll0_1782[k]
                    + f_16 * slk_1431[k]
                    - f_14 * pc_x[k] * sll1_1782[k];

        t_1783[k] = f_15 * slk_1432[k]
                    + f_3 * pc_x[k] * smk_1432[k];

        t_1784[k] = f_15 * slk_1433[k]
                    + f_3 * pc_x[k] * smk_1433[k];
    }

#pragma omp simd aligned(t_1785, t_1786, t_1787, t_1788, t_1789, pc_x, slk_1434, slk_1435, \
                         slk_1436, slk_1437, slk_1438, smk_1434, smk_1435, smk_1436, smk_1437, \
                         smk_1438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1785[k] = f_15 * slk_1434[k]
                    + f_3 * pc_x[k] * smk_1434[k];

        t_1786[k] = f_15 * slk_1435[k]
                    + f_3 * pc_x[k] * smk_1435[k];

        t_1787[k] = f_15 * slk_1436[k]
                    + f_3 * pc_x[k] * smk_1436[k];

        t_1788[k] = f_15 * slk_1437[k]
                    + f_3 * pc_x[k] * smk_1437[k];

        t_1789[k] = f_15 * slk_1438[k]
                    + f_3 * pc_x[k] * smk_1438[k];
    }

#pragma omp simd aligned(t_1790, t_1791, t_1792, t_1793, pb_x, pc_x, pc_z, sll0_1791, \
                         sll0_1793, slk_1108, slk_1439, sll1_1791, sll1_1793, smk_1432, \
                         smk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1790[k] = f_15 * slk_1439[k]
                    + f_3 * pc_x[k] * smk_1439[k];

        t_1791[k] = pb_x[k] * sll0_1791[k]
                    - f_14 * pc_x[k] * sll1_1791[k];

        t_1792[k] = f_17 * slk_1108[k]
                    + f_3 * pc_z[k] * smk_1432[k];

        t_1793[k] = pb_x[k] * sll0_1793[k]
                    - f_14 * pc_x[k] * sll1_1793[k];
    }

#pragma omp simd aligned(t_1794, t_1795, t_1796, t_1797, pb_x, pc_x, sll0_1794, sll0_1795, \
                         sll0_1796, sll0_1797, sll1_1794, sll1_1795, sll1_1796, \
                         sll1_1797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1794[k] = pb_x[k] * sll0_1794[k]
                    - f_14 * pc_x[k] * sll1_1794[k];

        t_1795[k] = pb_x[k] * sll0_1795[k]
                    - f_14 * pc_x[k] * sll1_1795[k];

        t_1796[k] = pb_x[k] * sll0_1796[k]
                    - f_14 * pc_x[k] * sll1_1796[k];

        t_1797[k] = pb_x[k] * sll0_1797[k]
                    - f_14 * pc_x[k] * sll1_1797[k];
    }

#pragma omp simd aligned(t_1798, t_1799, t_1800, t_1801, pb_x, pc_x, pc_y, sll0_1799, \
                         sll0_1800, slk_1151, slk_1152, slk_1440, sll1_1799, sll1_1800, \
                         smk_1439, smk_1440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1798[k] = f_19 * slk_1151[k]
                    + f_3 * pc_y[k] * smk_1439[k];

        t_1799[k] = pb_x[k] * sll0_1799[k]
                    - f_14 * pc_x[k] * sll1_1799[k];

        t_1800[k] = pb_x[k] * sll0_1800[k]
                    + f_21 * slk_1440[k]
                    - f_14 * pc_x[k] * sll1_1800[k];

        t_1801[k] = f_18 * slk_1152[k]
                    + f_3 * pc_y[k] * smk_1440[k];
    }

#pragma omp simd aligned(t_1802, t_1803, t_1804, pb_x, pc_x, pc_y, pc_z, sll0_1803, slk_1116, \
                         slk_1154, slk_1443, sll1_1803, smk_1440, \
                         smk_1442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1802[k] = f_18 * slk_1116[k]
                    + f_3 * pc_z[k] * smk_1440[k];

        t_1803[k] = pb_x[k] * sll0_1803[k]
                    + f_20 * slk_1443[k]
                    - f_14 * pc_x[k] * sll1_1803[k];

        t_1804[k] = f_18 * slk_1154[k]
                    + f_3 * pc_y[k] * smk_1442[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, pb_x, pc_x, pc_z, sll0_1805, sll0_1806, \
                         slk_1119, slk_1445, slk_1446, sll1_1805, sll1_1806, \
                         smk_1443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = pb_x[k] * sll0_1805[k]
                    + f_20 * slk_1445[k]
                    - f_14 * pc_x[k] * sll1_1805[k];

        t_1806[k] = pb_x[k] * sll0_1806[k]
                    + f_19 * slk_1446[k]
                    - f_14 * pc_x[k] * sll1_1806[k];

        t_1807[k] = f_18 * slk_1119[k]
                    + f_3 * pc_z[k] * smk_1443[k];
    }

#pragma omp simd aligned(t_1808, t_1809, t_1810, pb_x, pc_x, pc_y, sll0_1809, sll0_1810, \
                         slk_1157, slk_1449, slk_1450, sll1_1809, sll1_1810, \
                         smk_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1808[k] = f_18 * slk_1157[k]
                    + f_3 * pc_y[k] * smk_1445[k];

        t_1809[k] = pb_x[k] * sll0_1809[k]
                    + f_19 * slk_1449[k]
                    - f_14 * pc_x[k] * sll1_1809[k];

        t_1810[k] = pb_x[k] * sll0_1810[k]
                    + f_18 * slk_1450[k]
                    - f_14 * pc_x[k] * sll1_1810[k];
    }

#pragma omp simd aligned(t_1811, t_1812, t_1813, pb_x, pc_x, pc_y, pc_z, sll0_1812, slk_1122, \
                         slk_1161, slk_1452, sll1_1812, smk_1446, \
                         smk_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1811[k] = f_18 * slk_1122[k]
                    + f_3 * pc_z[k] * smk_1446[k];

        t_1812[k] = pb_x[k] * sll0_1812[k]
                    + f_18 * slk_1452[k]
                    - f_14 * pc_x[k] * sll1_1812[k];

        t_1813[k] = f_18 * slk_1161[k]
                    + f_3 * pc_y[k] * smk_1449[k];
    }

#pragma omp simd aligned(t_1814, t_1815, t_1816, pb_x, pc_x, pc_z, sll0_1814, sll0_1815, \
                         slk_1126, slk_1454, slk_1455, sll1_1814, sll1_1815, \
                         smk_1450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1814[k] = pb_x[k] * sll0_1814[k]
                    + f_18 * slk_1454[k]
                    - f_14 * pc_x[k] * sll1_1814[k];

        t_1815[k] = pb_x[k] * sll0_1815[k]
                    + f_17 * slk_1455[k]
                    - f_14 * pc_x[k] * sll1_1815[k];

        t_1816[k] = f_18 * slk_1126[k]
                    + f_3 * pc_z[k] * smk_1450[k];
    }

#pragma omp simd aligned(t_1817, t_1818, t_1819, pb_x, pc_x, pc_y, sll0_1817, sll0_1818, \
                         slk_1166, slk_1457, slk_1458, sll1_1817, sll1_1818, \
                         smk_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1817[k] = pb_x[k] * sll0_1817[k]
                    + f_17 * slk_1457[k]
                    - f_14 * pc_x[k] * sll1_1817[k];

        t_1818[k] = pb_x[k] * sll0_1818[k]
                    + f_17 * slk_1458[k]
                    - f_14 * pc_x[k] * sll1_1818[k];

        t_1819[k] = f_18 * slk_1166[k]
                    + f_3 * pc_y[k] * smk_1454[k];
    }

#pragma omp simd aligned(t_1820, t_1821, t_1822, pb_x, pc_x, pc_z, sll0_1820, sll0_1821, \
                         slk_1131, slk_1460, slk_1461, sll1_1820, sll1_1821, \
                         smk_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1820[k] = pb_x[k] * sll0_1820[k]
                    + f_17 * slk_1460[k]
                    - f_14 * pc_x[k] * sll1_1820[k];

        t_1821[k] = pb_x[k] * sll0_1821[k]
                    + f_16 * slk_1461[k]
                    - f_14 * pc_x[k] * sll1_1821[k];

        t_1822[k] = f_18 * slk_1131[k]
                    + f_3 * pc_z[k] * smk_1455[k];
    }

#pragma omp simd aligned(t_1823, t_1824, t_1825, pb_x, pc_x, sll0_1823, sll0_1824, sll0_1825, \
                         slk_1463, slk_1464, slk_1465, sll1_1823, sll1_1824, \
                         sll1_1825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1823[k] = pb_x[k] * sll0_1823[k]
                    + f_16 * slk_1463[k]
                    - f_14 * pc_x[k] * sll1_1823[k];

        t_1824[k] = pb_x[k] * sll0_1824[k]
                    + f_16 * slk_1464[k]
                    - f_14 * pc_x[k] * sll1_1824[k];

        t_1825[k] = pb_x[k] * sll0_1825[k]
                    + f_16 * slk_1465[k]
                    - f_14 * pc_x[k] * sll1_1825[k];
    }

#pragma omp simd aligned(t_1826, t_1827, t_1828, t_1829, pb_x, pc_x, pc_y, sll0_1827, \
                         slk_1172, slk_1467, slk_1468, slk_1469, sll1_1827, smk_1460, \
                         smk_1468, smk_1469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1826[k] = f_18 * slk_1172[k]
                    + f_3 * pc_y[k] * smk_1460[k];

        t_1827[k] = pb_x[k] * sll0_1827[k]
                    + f_16 * slk_1467[k]
                    - f_14 * pc_x[k] * sll1_1827[k];

        t_1828[k] = f_15 * slk_1468[k]
                    + f_3 * pc_x[k] * smk_1468[k];

        t_1829[k] = f_15 * slk_1469[k]
                    + f_3 * pc_x[k] * smk_1469[k];
    }

#pragma omp simd aligned(t_1830, t_1831, t_1832, t_1833, t_1834, pc_x, slk_1470, slk_1471, \
                         slk_1472, slk_1473, slk_1474, smk_1470, smk_1471, smk_1472, smk_1473, \
                         smk_1474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1830[k] = f_15 * slk_1470[k]
                    + f_3 * pc_x[k] * smk_1470[k];

        t_1831[k] = f_15 * slk_1471[k]
                    + f_3 * pc_x[k] * smk_1471[k];

        t_1832[k] = f_15 * slk_1472[k]
                    + f_3 * pc_x[k] * smk_1472[k];

        t_1833[k] = f_15 * slk_1473[k]
                    + f_3 * pc_x[k] * smk_1473[k];

        t_1834[k] = f_15 * slk_1474[k]
                    + f_3 * pc_x[k] * smk_1474[k];
    }

#pragma omp simd aligned(t_1835, t_1836, t_1837, t_1838, pb_x, pc_x, pc_z, sll0_1836, \
                         sll0_1838, slk_1144, slk_1475, sll1_1836, sll1_1838, smk_1468, \
                         smk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1835[k] = f_15 * slk_1475[k]
                    + f_3 * pc_x[k] * smk_1475[k];

        t_1836[k] = pb_x[k] * sll0_1836[k]
                    - f_14 * pc_x[k] * sll1_1836[k];

        t_1837[k] = f_18 * slk_1144[k]
                    + f_3 * pc_z[k] * smk_1468[k];

        t_1838[k] = pb_x[k] * sll0_1838[k]
                    - f_14 * pc_x[k] * sll1_1838[k];
    }

#pragma omp simd aligned(t_1839, t_1840, t_1841, t_1842, pb_x, pc_x, sll0_1839, sll0_1840, \
                         sll0_1841, sll0_1842, sll1_1839, sll1_1840, sll1_1841, \
                         sll1_1842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1839[k] = pb_x[k] * sll0_1839[k]
                    - f_14 * pc_x[k] * sll1_1839[k];

        t_1840[k] = pb_x[k] * sll0_1840[k]
                    - f_14 * pc_x[k] * sll1_1840[k];

        t_1841[k] = pb_x[k] * sll0_1841[k]
                    - f_14 * pc_x[k] * sll1_1841[k];

        t_1842[k] = pb_x[k] * sll0_1842[k]
                    - f_14 * pc_x[k] * sll1_1842[k];
    }

#pragma omp simd aligned(t_1843, t_1844, t_1845, t_1846, pb_x, pc_x, pc_y, sll0_1844, \
                         sll0_1845, slk_1187, slk_1188, slk_1476, sll1_1844, sll1_1845, \
                         smk_1475, smk_1476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1843[k] = f_18 * slk_1187[k]
                    + f_3 * pc_y[k] * smk_1475[k];

        t_1844[k] = pb_x[k] * sll0_1844[k]
                    - f_14 * pc_x[k] * sll1_1844[k];

        t_1845[k] = pb_x[k] * sll0_1845[k]
                    + f_21 * slk_1476[k]
                    - f_14 * pc_x[k] * sll1_1845[k];

        t_1846[k] = f_17 * slk_1188[k]
                    + f_3 * pc_y[k] * smk_1476[k];
    }

#pragma omp simd aligned(t_1847, t_1848, t_1849, pb_x, pc_x, pc_y, pc_z, sll0_1848, slk_1152, \
                         slk_1190, slk_1479, sll1_1848, smk_1476, \
                         smk_1478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1847[k] = f_19 * slk_1152[k]
                    + f_3 * pc_z[k] * smk_1476[k];

        t_1848[k] = pb_x[k] * sll0_1848[k]
                    + f_20 * slk_1479[k]
                    - f_14 * pc_x[k] * sll1_1848[k];

        t_1849[k] = f_17 * slk_1190[k]
                    + f_3 * pc_y[k] * smk_1478[k];
    }

#pragma omp simd aligned(t_1850, t_1851, t_1852, pb_x, pc_x, pc_z, sll0_1850, sll0_1851, \
                         slk_1155, slk_1481, slk_1482, sll1_1850, sll1_1851, \
                         smk_1479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1850[k] = pb_x[k] * sll0_1850[k]
                    + f_20 * slk_1481[k]
                    - f_14 * pc_x[k] * sll1_1850[k];

        t_1851[k] = pb_x[k] * sll0_1851[k]
                    + f_19 * slk_1482[k]
                    - f_14 * pc_x[k] * sll1_1851[k];

        t_1852[k] = f_19 * slk_1155[k]
                    + f_3 * pc_z[k] * smk_1479[k];
    }

#pragma omp simd aligned(t_1853, t_1854, t_1855, pb_x, pc_x, pc_y, sll0_1854, sll0_1855, \
                         slk_1193, slk_1485, slk_1486, sll1_1854, sll1_1855, \
                         smk_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1853[k] = f_17 * slk_1193[k]
                    + f_3 * pc_y[k] * smk_1481[k];

        t_1854[k] = pb_x[k] * sll0_1854[k]
                    + f_19 * slk_1485[k]
                    - f_14 * pc_x[k] * sll1_1854[k];

        t_1855[k] = pb_x[k] * sll0_1855[k]
                    + f_18 * slk_1486[k]
                    - f_14 * pc_x[k] * sll1_1855[k];
    }

#pragma omp simd aligned(t_1856, t_1857, t_1858, pb_x, pc_x, pc_y, pc_z, sll0_1857, slk_1158, \
                         slk_1197, slk_1488, sll1_1857, smk_1482, \
                         smk_1485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1856[k] = f_19 * slk_1158[k]
                    + f_3 * pc_z[k] * smk_1482[k];

        t_1857[k] = pb_x[k] * sll0_1857[k]
                    + f_18 * slk_1488[k]
                    - f_14 * pc_x[k] * sll1_1857[k];

        t_1858[k] = f_17 * slk_1197[k]
                    + f_3 * pc_y[k] * smk_1485[k];
    }

#pragma omp simd aligned(t_1859, t_1860, t_1861, pb_x, pc_x, pc_z, sll0_1859, sll0_1860, \
                         slk_1162, slk_1490, slk_1491, sll1_1859, sll1_1860, \
                         smk_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1859[k] = pb_x[k] * sll0_1859[k]
                    + f_18 * slk_1490[k]
                    - f_14 * pc_x[k] * sll1_1859[k];

        t_1860[k] = pb_x[k] * sll0_1860[k]
                    + f_17 * slk_1491[k]
                    - f_14 * pc_x[k] * sll1_1860[k];

        t_1861[k] = f_19 * slk_1162[k]
                    + f_3 * pc_z[k] * smk_1486[k];
    }

#pragma omp simd aligned(t_1862, t_1863, t_1864, pb_x, pc_x, pc_y, sll0_1862, sll0_1863, \
                         slk_1202, slk_1493, slk_1494, sll1_1862, sll1_1863, \
                         smk_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1862[k] = pb_x[k] * sll0_1862[k]
                    + f_17 * slk_1493[k]
                    - f_14 * pc_x[k] * sll1_1862[k];

        t_1863[k] = pb_x[k] * sll0_1863[k]
                    + f_17 * slk_1494[k]
                    - f_14 * pc_x[k] * sll1_1863[k];

        t_1864[k] = f_17 * slk_1202[k]
                    + f_3 * pc_y[k] * smk_1490[k];
    }

#pragma omp simd aligned(t_1865, t_1866, t_1867, pb_x, pc_x, pc_z, sll0_1865, sll0_1866, \
                         slk_1167, slk_1496, slk_1497, sll1_1865, sll1_1866, \
                         smk_1491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1865[k] = pb_x[k] * sll0_1865[k]
                    + f_17 * slk_1496[k]
                    - f_14 * pc_x[k] * sll1_1865[k];

        t_1866[k] = pb_x[k] * sll0_1866[k]
                    + f_16 * slk_1497[k]
                    - f_14 * pc_x[k] * sll1_1866[k];

        t_1867[k] = f_19 * slk_1167[k]
                    + f_3 * pc_z[k] * smk_1491[k];
    }

#pragma omp simd aligned(t_1868, t_1869, t_1870, pb_x, pc_x, sll0_1868, sll0_1869, sll0_1870, \
                         slk_1499, slk_1500, slk_1501, sll1_1868, sll1_1869, \
                         sll1_1870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1868[k] = pb_x[k] * sll0_1868[k]
                    + f_16 * slk_1499[k]
                    - f_14 * pc_x[k] * sll1_1868[k];

        t_1869[k] = pb_x[k] * sll0_1869[k]
                    + f_16 * slk_1500[k]
                    - f_14 * pc_x[k] * sll1_1869[k];

        t_1870[k] = pb_x[k] * sll0_1870[k]
                    + f_16 * slk_1501[k]
                    - f_14 * pc_x[k] * sll1_1870[k];
    }

#pragma omp simd aligned(t_1871, t_1872, t_1873, t_1874, pb_x, pc_x, pc_y, sll0_1872, \
                         slk_1208, slk_1503, slk_1504, slk_1505, sll1_1872, smk_1496, \
                         smk_1504, smk_1505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1871[k] = f_17 * slk_1208[k]
                    + f_3 * pc_y[k] * smk_1496[k];

        t_1872[k] = pb_x[k] * sll0_1872[k]
                    + f_16 * slk_1503[k]
                    - f_14 * pc_x[k] * sll1_1872[k];

        t_1873[k] = f_15 * slk_1504[k]
                    + f_3 * pc_x[k] * smk_1504[k];

        t_1874[k] = f_15 * slk_1505[k]
                    + f_3 * pc_x[k] * smk_1505[k];
    }

#pragma omp simd aligned(t_1875, t_1876, t_1877, t_1878, t_1879, pc_x, slk_1506, slk_1507, \
                         slk_1508, slk_1509, slk_1510, smk_1506, smk_1507, smk_1508, smk_1509, \
                         smk_1510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1875[k] = f_15 * slk_1506[k]
                    + f_3 * pc_x[k] * smk_1506[k];

        t_1876[k] = f_15 * slk_1507[k]
                    + f_3 * pc_x[k] * smk_1507[k];

        t_1877[k] = f_15 * slk_1508[k]
                    + f_3 * pc_x[k] * smk_1508[k];

        t_1878[k] = f_15 * slk_1509[k]
                    + f_3 * pc_x[k] * smk_1509[k];

        t_1879[k] = f_15 * slk_1510[k]
                    + f_3 * pc_x[k] * smk_1510[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_1575 = buffer.data(sll0 + 1575);
    const auto *sll0_1580 = buffer.data(sll0 + 1580);
    const auto *sll0_1584 = buffer.data(sll0 + 1584);
    const auto *sll0_1589 = buffer.data(sll0 + 1589);
    const auto *sll0_1595 = buffer.data(sll0 + 1595);
    const auto *sll0_1602 = buffer.data(sll0 + 1602);
    const auto *sll0_1881 = buffer.data(sll0 + 1881);
    const auto *sll0_1883 = buffer.data(sll0 + 1883);
    const auto *sll0_1884 = buffer.data(sll0 + 1884);
    const auto *sll0_1885 = buffer.data(sll0 + 1885);
    const auto *sll0_1886 = buffer.data(sll0 + 1886);
    const auto *sll0_1887 = buffer.data(sll0 + 1887);
    const auto *sll0_1889 = buffer.data(sll0 + 1889);
    const auto *sll0_1890 = buffer.data(sll0 + 1890);
    const auto *sll0_1893 = buffer.data(sll0 + 1893);
    const auto *sll0_1895 = buffer.data(sll0 + 1895);
    const auto *sll0_1896 = buffer.data(sll0 + 1896);
    const auto *sll0_1899 = buffer.data(sll0 + 1899);
    const auto *sll0_1900 = buffer.data(sll0 + 1900);
    const auto *sll0_1902 = buffer.data(sll0 + 1902);
    const auto *sll0_1904 = buffer.data(sll0 + 1904);
    const auto *sll0_1905 = buffer.data(sll0 + 1905);
    const auto *sll0_1907 = buffer.data(sll0 + 1907);
    const auto *sll0_1908 = buffer.data(sll0 + 1908);
    const auto *sll0_1910 = buffer.data(sll0 + 1910);
    const auto *sll0_1911 = buffer.data(sll0 + 1911);
    const auto *sll0_1913 = buffer.data(sll0 + 1913);
    const auto *sll0_1914 = buffer.data(sll0 + 1914);
    const auto *sll0_1915 = buffer.data(sll0 + 1915);
    const auto *sll0_1917 = buffer.data(sll0 + 1917);
    const auto *sll0_1926 = buffer.data(sll0 + 1926);
    const auto *sll0_1928 = buffer.data(sll0 + 1928);
    const auto *sll0_1929 = buffer.data(sll0 + 1929);
    const auto *sll0_1930 = buffer.data(sll0 + 1930);
    const auto *sll0_1931 = buffer.data(sll0 + 1931);
    const auto *sll0_1932 = buffer.data(sll0 + 1932);
    const auto *sll0_1934 = buffer.data(sll0 + 1934);
    const auto *sll0_1938 = buffer.data(sll0 + 1938);
    const auto *sll0_1941 = buffer.data(sll0 + 1941);
    const auto *sll0_1945 = buffer.data(sll0 + 1945);
    const auto *sll0_1947 = buffer.data(sll0 + 1947);
    const auto *sll0_1950 = buffer.data(sll0 + 1950);
    const auto *sll0_1952 = buffer.data(sll0 + 1952);
    const auto *sll0_1953 = buffer.data(sll0 + 1953);
    const auto *sll0_1956 = buffer.data(sll0 + 1956);
    const auto *sll0_1958 = buffer.data(sll0 + 1958);
    const auto *sll0_1959 = buffer.data(sll0 + 1959);
    const auto *sll0_1960 = buffer.data(sll0 + 1960);
    const auto *sll0_1971 = buffer.data(sll0 + 1971);
    const auto *sll0_1973 = buffer.data(sll0 + 1973);
    const auto *sll0_1974 = buffer.data(sll0 + 1974);
    const auto *sll0_1975 = buffer.data(sll0 + 1975);
    const auto *sll0_1976 = buffer.data(sll0 + 1976);
    const auto *sll0_1977 = buffer.data(sll0 + 1977);
    const auto *sll0_1979 = buffer.data(sll0 + 1979);
    const auto *sll0_1980 = buffer.data(sll0 + 1980);
    const auto *sll0_1983 = buffer.data(sll0 + 1983);
    const auto *sll0_1985 = buffer.data(sll0 + 1985);
    const auto *sll0_1986 = buffer.data(sll0 + 1986);
    const auto *sll0_1989 = buffer.data(sll0 + 1989);
    const auto *sll0_1990 = buffer.data(sll0 + 1990);
    const auto *sll0_1992 = buffer.data(sll0 + 1992);
    const auto *sll0_1994 = buffer.data(sll0 + 1994);
    const auto *sll0_1995 = buffer.data(sll0 + 1995);

    const auto *slk_1180 = buffer.data(slk + 1180);
    const auto *slk_1188 = buffer.data(slk + 1188);
    const auto *slk_1191 = buffer.data(slk + 1191);
    const auto *slk_1194 = buffer.data(slk + 1194);
    const auto *slk_1198 = buffer.data(slk + 1198);
    const auto *slk_1203 = buffer.data(slk + 1203);
    const auto *slk_1216 = buffer.data(slk + 1216);
    const auto *slk_1223 = buffer.data(slk + 1223);
    const auto *slk_1224 = buffer.data(slk + 1224);
    const auto *slk_1226 = buffer.data(slk + 1226);
    const auto *slk_1227 = buffer.data(slk + 1227);
    const auto *slk_1229 = buffer.data(slk + 1229);
    const auto *slk_1230 = buffer.data(slk + 1230);
    const auto *slk_1233 = buffer.data(slk + 1233);
    const auto *slk_1234 = buffer.data(slk + 1234);
    const auto *slk_1238 = buffer.data(slk + 1238);
    const auto *slk_1239 = buffer.data(slk + 1239);
    const auto *slk_1244 = buffer.data(slk + 1244);
    const auto *slk_1252 = buffer.data(slk + 1252);
    const auto *slk_1259 = buffer.data(slk + 1259);
    const auto *slk_1260 = buffer.data(slk + 1260);
    const auto *slk_1262 = buffer.data(slk + 1262);
    const auto *slk_1263 = buffer.data(slk + 1263);
    const auto *slk_1265 = buffer.data(slk + 1265);
    const auto *slk_1266 = buffer.data(slk + 1266);
    const auto *slk_1269 = buffer.data(slk + 1269);
    const auto *slk_1270 = buffer.data(slk + 1270);
    const auto *slk_1274 = buffer.data(slk + 1274);
    const auto *slk_1280 = buffer.data(slk + 1280);
    const auto *slk_1295 = buffer.data(slk + 1295);
    const auto *slk_1511 = buffer.data(slk + 1511);
    const auto *slk_1512 = buffer.data(slk + 1512);
    const auto *slk_1515 = buffer.data(slk + 1515);
    const auto *slk_1517 = buffer.data(slk + 1517);
    const auto *slk_1518 = buffer.data(slk + 1518);
    const auto *slk_1521 = buffer.data(slk + 1521);
    const auto *slk_1522 = buffer.data(slk + 1522);
    const auto *slk_1524 = buffer.data(slk + 1524);
    const auto *slk_1526 = buffer.data(slk + 1526);
    const auto *slk_1527 = buffer.data(slk + 1527);
    const auto *slk_1529 = buffer.data(slk + 1529);
    const auto *slk_1530 = buffer.data(slk + 1530);
    const auto *slk_1532 = buffer.data(slk + 1532);
    const auto *slk_1533 = buffer.data(slk + 1533);
    const auto *slk_1535 = buffer.data(slk + 1535);
    const auto *slk_1536 = buffer.data(slk + 1536);
    const auto *slk_1537 = buffer.data(slk + 1537);
    const auto *slk_1539 = buffer.data(slk + 1539);
    const auto *slk_1540 = buffer.data(slk + 1540);
    const auto *slk_1541 = buffer.data(slk + 1541);
    const auto *slk_1542 = buffer.data(slk + 1542);
    const auto *slk_1543 = buffer.data(slk + 1543);
    const auto *slk_1544 = buffer.data(slk + 1544);
    const auto *slk_1545 = buffer.data(slk + 1545);
    const auto *slk_1546 = buffer.data(slk + 1546);
    const auto *slk_1547 = buffer.data(slk + 1547);
    const auto *slk_1551 = buffer.data(slk + 1551);
    const auto *slk_1554 = buffer.data(slk + 1554);
    const auto *slk_1558 = buffer.data(slk + 1558);
    const auto *slk_1560 = buffer.data(slk + 1560);
    const auto *slk_1563 = buffer.data(slk + 1563);
    const auto *slk_1565 = buffer.data(slk + 1565);
    const auto *slk_1566 = buffer.data(slk + 1566);
    const auto *slk_1569 = buffer.data(slk + 1569);
    const auto *slk_1571 = buffer.data(slk + 1571);
    const auto *slk_1572 = buffer.data(slk + 1572);
    const auto *slk_1573 = buffer.data(slk + 1573);
    const auto *slk_1576 = buffer.data(slk + 1576);
    const auto *slk_1577 = buffer.data(slk + 1577);
    const auto *slk_1578 = buffer.data(slk + 1578);
    const auto *slk_1579 = buffer.data(slk + 1579);
    const auto *slk_1580 = buffer.data(slk + 1580);
    const auto *slk_1581 = buffer.data(slk + 1581);
    const auto *slk_1582 = buffer.data(slk + 1582);
    const auto *slk_1583 = buffer.data(slk + 1583);
    const auto *slk_1584 = buffer.data(slk + 1584);
    const auto *slk_1587 = buffer.data(slk + 1587);
    const auto *slk_1589 = buffer.data(slk + 1589);
    const auto *slk_1590 = buffer.data(slk + 1590);
    const auto *slk_1593 = buffer.data(slk + 1593);
    const auto *slk_1594 = buffer.data(slk + 1594);
    const auto *slk_1596 = buffer.data(slk + 1596);
    const auto *slk_1598 = buffer.data(slk + 1598);
    const auto *slk_1599 = buffer.data(slk + 1599);

    const auto *sll1_1575 = buffer.data(sll1 + 1575);
    const auto *sll1_1580 = buffer.data(sll1 + 1580);
    const auto *sll1_1584 = buffer.data(sll1 + 1584);
    const auto *sll1_1589 = buffer.data(sll1 + 1589);
    const auto *sll1_1595 = buffer.data(sll1 + 1595);
    const auto *sll1_1602 = buffer.data(sll1 + 1602);
    const auto *sll1_1881 = buffer.data(sll1 + 1881);
    const auto *sll1_1883 = buffer.data(sll1 + 1883);
    const auto *sll1_1884 = buffer.data(sll1 + 1884);
    const auto *sll1_1885 = buffer.data(sll1 + 1885);
    const auto *sll1_1886 = buffer.data(sll1 + 1886);
    const auto *sll1_1887 = buffer.data(sll1 + 1887);
    const auto *sll1_1889 = buffer.data(sll1 + 1889);
    const auto *sll1_1890 = buffer.data(sll1 + 1890);
    const auto *sll1_1893 = buffer.data(sll1 + 1893);
    const auto *sll1_1895 = buffer.data(sll1 + 1895);
    const auto *sll1_1896 = buffer.data(sll1 + 1896);
    const auto *sll1_1899 = buffer.data(sll1 + 1899);
    const auto *sll1_1900 = buffer.data(sll1 + 1900);
    const auto *sll1_1902 = buffer.data(sll1 + 1902);
    const auto *sll1_1904 = buffer.data(sll1 + 1904);
    const auto *sll1_1905 = buffer.data(sll1 + 1905);
    const auto *sll1_1907 = buffer.data(sll1 + 1907);
    const auto *sll1_1908 = buffer.data(sll1 + 1908);
    const auto *sll1_1910 = buffer.data(sll1 + 1910);
    const auto *sll1_1911 = buffer.data(sll1 + 1911);
    const auto *sll1_1913 = buffer.data(sll1 + 1913);
    const auto *sll1_1914 = buffer.data(sll1 + 1914);
    const auto *sll1_1915 = buffer.data(sll1 + 1915);
    const auto *sll1_1917 = buffer.data(sll1 + 1917);
    const auto *sll1_1926 = buffer.data(sll1 + 1926);
    const auto *sll1_1928 = buffer.data(sll1 + 1928);
    const auto *sll1_1929 = buffer.data(sll1 + 1929);
    const auto *sll1_1930 = buffer.data(sll1 + 1930);
    const auto *sll1_1931 = buffer.data(sll1 + 1931);
    const auto *sll1_1932 = buffer.data(sll1 + 1932);
    const auto *sll1_1934 = buffer.data(sll1 + 1934);
    const auto *sll1_1938 = buffer.data(sll1 + 1938);
    const auto *sll1_1941 = buffer.data(sll1 + 1941);
    const auto *sll1_1945 = buffer.data(sll1 + 1945);
    const auto *sll1_1947 = buffer.data(sll1 + 1947);
    const auto *sll1_1950 = buffer.data(sll1 + 1950);
    const auto *sll1_1952 = buffer.data(sll1 + 1952);
    const auto *sll1_1953 = buffer.data(sll1 + 1953);
    const auto *sll1_1956 = buffer.data(sll1 + 1956);
    const auto *sll1_1958 = buffer.data(sll1 + 1958);
    const auto *sll1_1959 = buffer.data(sll1 + 1959);
    const auto *sll1_1960 = buffer.data(sll1 + 1960);
    const auto *sll1_1971 = buffer.data(sll1 + 1971);
    const auto *sll1_1973 = buffer.data(sll1 + 1973);
    const auto *sll1_1974 = buffer.data(sll1 + 1974);
    const auto *sll1_1975 = buffer.data(sll1 + 1975);
    const auto *sll1_1976 = buffer.data(sll1 + 1976);
    const auto *sll1_1977 = buffer.data(sll1 + 1977);
    const auto *sll1_1979 = buffer.data(sll1 + 1979);
    const auto *sll1_1980 = buffer.data(sll1 + 1980);
    const auto *sll1_1983 = buffer.data(sll1 + 1983);
    const auto *sll1_1985 = buffer.data(sll1 + 1985);
    const auto *sll1_1986 = buffer.data(sll1 + 1986);
    const auto *sll1_1989 = buffer.data(sll1 + 1989);
    const auto *sll1_1990 = buffer.data(sll1 + 1990);
    const auto *sll1_1992 = buffer.data(sll1 + 1992);
    const auto *sll1_1994 = buffer.data(sll1 + 1994);
    const auto *sll1_1995 = buffer.data(sll1 + 1995);

    const auto *smk_1504 = buffer.data(smk + 1504);
    const auto *smk_1511 = buffer.data(smk + 1511);
    const auto *smk_1512 = buffer.data(smk + 1512);
    const auto *smk_1514 = buffer.data(smk + 1514);
    const auto *smk_1515 = buffer.data(smk + 1515);
    const auto *smk_1517 = buffer.data(smk + 1517);
    const auto *smk_1518 = buffer.data(smk + 1518);
    const auto *smk_1521 = buffer.data(smk + 1521);
    const auto *smk_1522 = buffer.data(smk + 1522);
    const auto *smk_1526 = buffer.data(smk + 1526);
    const auto *smk_1527 = buffer.data(smk + 1527);
    const auto *smk_1532 = buffer.data(smk + 1532);
    const auto *smk_1540 = buffer.data(smk + 1540);
    const auto *smk_1541 = buffer.data(smk + 1541);
    const auto *smk_1542 = buffer.data(smk + 1542);
    const auto *smk_1543 = buffer.data(smk + 1543);
    const auto *smk_1544 = buffer.data(smk + 1544);
    const auto *smk_1545 = buffer.data(smk + 1545);
    const auto *smk_1546 = buffer.data(smk + 1546);
    const auto *smk_1547 = buffer.data(smk + 1547);
    const auto *smk_1548 = buffer.data(smk + 1548);
    const auto *smk_1550 = buffer.data(smk + 1550);
    const auto *smk_1551 = buffer.data(smk + 1551);
    const auto *smk_1553 = buffer.data(smk + 1553);
    const auto *smk_1554 = buffer.data(smk + 1554);
    const auto *smk_1557 = buffer.data(smk + 1557);
    const auto *smk_1558 = buffer.data(smk + 1558);
    const auto *smk_1562 = buffer.data(smk + 1562);
    const auto *smk_1563 = buffer.data(smk + 1563);
    const auto *smk_1568 = buffer.data(smk + 1568);
    const auto *smk_1576 = buffer.data(smk + 1576);
    const auto *smk_1577 = buffer.data(smk + 1577);
    const auto *smk_1578 = buffer.data(smk + 1578);
    const auto *smk_1579 = buffer.data(smk + 1579);
    const auto *smk_1580 = buffer.data(smk + 1580);
    const auto *smk_1581 = buffer.data(smk + 1581);
    const auto *smk_1582 = buffer.data(smk + 1582);
    const auto *smk_1583 = buffer.data(smk + 1583);
    const auto *smk_1584 = buffer.data(smk + 1584);
    const auto *smk_1586 = buffer.data(smk + 1586);
    const auto *smk_1587 = buffer.data(smk + 1587);
    const auto *smk_1589 = buffer.data(smk + 1589);
    const auto *smk_1590 = buffer.data(smk + 1590);
    const auto *smk_1593 = buffer.data(smk + 1593);
    const auto *smk_1594 = buffer.data(smk + 1594);

#pragma omp simd aligned(t_1880, t_1881, t_1882, t_1883, pb_x, pc_x, pc_z, sll0_1881, \
                         sll0_1883, slk_1180, slk_1511, sll1_1881, sll1_1883, smk_1504, \
                         smk_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1880[k] = f_15 * slk_1511[k]
                    + f_3 * pc_x[k] * smk_1511[k];

        t_1881[k] = pb_x[k] * sll0_1881[k]
                    - f_14 * pc_x[k] * sll1_1881[k];

        t_1882[k] = f_19 * slk_1180[k]
                    + f_3 * pc_z[k] * smk_1504[k];

        t_1883[k] = pb_x[k] * sll0_1883[k]
                    - f_14 * pc_x[k] * sll1_1883[k];
    }

#pragma omp simd aligned(t_1884, t_1885, t_1886, t_1887, pb_x, pc_x, sll0_1884, sll0_1885, \
                         sll0_1886, sll0_1887, sll1_1884, sll1_1885, sll1_1886, \
                         sll1_1887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1884[k] = pb_x[k] * sll0_1884[k]
                    - f_14 * pc_x[k] * sll1_1884[k];

        t_1885[k] = pb_x[k] * sll0_1885[k]
                    - f_14 * pc_x[k] * sll1_1885[k];

        t_1886[k] = pb_x[k] * sll0_1886[k]
                    - f_14 * pc_x[k] * sll1_1886[k];

        t_1887[k] = pb_x[k] * sll0_1887[k]
                    - f_14 * pc_x[k] * sll1_1887[k];
    }

#pragma omp simd aligned(t_1888, t_1889, t_1890, t_1891, pb_x, pc_x, pc_y, sll0_1889, \
                         sll0_1890, slk_1223, slk_1224, slk_1512, sll1_1889, sll1_1890, \
                         smk_1511, smk_1512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1888[k] = f_17 * slk_1223[k]
                    + f_3 * pc_y[k] * smk_1511[k];

        t_1889[k] = pb_x[k] * sll0_1889[k]
                    - f_14 * pc_x[k] * sll1_1889[k];

        t_1890[k] = pb_x[k] * sll0_1890[k]
                    + f_21 * slk_1512[k]
                    - f_14 * pc_x[k] * sll1_1890[k];

        t_1891[k] = f_16 * slk_1224[k]
                    + f_3 * pc_y[k] * smk_1512[k];
    }

#pragma omp simd aligned(t_1892, t_1893, t_1894, pb_x, pc_x, pc_y, pc_z, sll0_1893, slk_1188, \
                         slk_1226, slk_1515, sll1_1893, smk_1512, \
                         smk_1514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1892[k] = f_20 * slk_1188[k]
                    + f_3 * pc_z[k] * smk_1512[k];

        t_1893[k] = pb_x[k] * sll0_1893[k]
                    + f_20 * slk_1515[k]
                    - f_14 * pc_x[k] * sll1_1893[k];

        t_1894[k] = f_16 * slk_1226[k]
                    + f_3 * pc_y[k] * smk_1514[k];
    }

#pragma omp simd aligned(t_1895, t_1896, t_1897, pb_x, pc_x, pc_z, sll0_1895, sll0_1896, \
                         slk_1191, slk_1517, slk_1518, sll1_1895, sll1_1896, \
                         smk_1515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1895[k] = pb_x[k] * sll0_1895[k]
                    + f_20 * slk_1517[k]
                    - f_14 * pc_x[k] * sll1_1895[k];

        t_1896[k] = pb_x[k] * sll0_1896[k]
                    + f_19 * slk_1518[k]
                    - f_14 * pc_x[k] * sll1_1896[k];

        t_1897[k] = f_20 * slk_1191[k]
                    + f_3 * pc_z[k] * smk_1515[k];
    }

#pragma omp simd aligned(t_1898, t_1899, t_1900, pb_x, pc_x, pc_y, sll0_1899, sll0_1900, \
                         slk_1229, slk_1521, slk_1522, sll1_1899, sll1_1900, \
                         smk_1517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1898[k] = f_16 * slk_1229[k]
                    + f_3 * pc_y[k] * smk_1517[k];

        t_1899[k] = pb_x[k] * sll0_1899[k]
                    + f_19 * slk_1521[k]
                    - f_14 * pc_x[k] * sll1_1899[k];

        t_1900[k] = pb_x[k] * sll0_1900[k]
                    + f_18 * slk_1522[k]
                    - f_14 * pc_x[k] * sll1_1900[k];
    }

#pragma omp simd aligned(t_1901, t_1902, t_1903, pb_x, pc_x, pc_y, pc_z, sll0_1902, slk_1194, \
                         slk_1233, slk_1524, sll1_1902, smk_1518, \
                         smk_1521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1901[k] = f_20 * slk_1194[k]
                    + f_3 * pc_z[k] * smk_1518[k];

        t_1902[k] = pb_x[k] * sll0_1902[k]
                    + f_18 * slk_1524[k]
                    - f_14 * pc_x[k] * sll1_1902[k];

        t_1903[k] = f_16 * slk_1233[k]
                    + f_3 * pc_y[k] * smk_1521[k];
    }

#pragma omp simd aligned(t_1904, t_1905, t_1906, pb_x, pc_x, pc_z, sll0_1904, sll0_1905, \
                         slk_1198, slk_1526, slk_1527, sll1_1904, sll1_1905, \
                         smk_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1904[k] = pb_x[k] * sll0_1904[k]
                    + f_18 * slk_1526[k]
                    - f_14 * pc_x[k] * sll1_1904[k];

        t_1905[k] = pb_x[k] * sll0_1905[k]
                    + f_17 * slk_1527[k]
                    - f_14 * pc_x[k] * sll1_1905[k];

        t_1906[k] = f_20 * slk_1198[k]
                    + f_3 * pc_z[k] * smk_1522[k];
    }

#pragma omp simd aligned(t_1907, t_1908, t_1909, pb_x, pc_x, pc_y, sll0_1907, sll0_1908, \
                         slk_1238, slk_1529, slk_1530, sll1_1907, sll1_1908, \
                         smk_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1907[k] = pb_x[k] * sll0_1907[k]
                    + f_17 * slk_1529[k]
                    - f_14 * pc_x[k] * sll1_1907[k];

        t_1908[k] = pb_x[k] * sll0_1908[k]
                    + f_17 * slk_1530[k]
                    - f_14 * pc_x[k] * sll1_1908[k];

        t_1909[k] = f_16 * slk_1238[k]
                    + f_3 * pc_y[k] * smk_1526[k];
    }

#pragma omp simd aligned(t_1910, t_1911, t_1912, pb_x, pc_x, pc_z, sll0_1910, sll0_1911, \
                         slk_1203, slk_1532, slk_1533, sll1_1910, sll1_1911, \
                         smk_1527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1910[k] = pb_x[k] * sll0_1910[k]
                    + f_17 * slk_1532[k]
                    - f_14 * pc_x[k] * sll1_1910[k];

        t_1911[k] = pb_x[k] * sll0_1911[k]
                    + f_16 * slk_1533[k]
                    - f_14 * pc_x[k] * sll1_1911[k];

        t_1912[k] = f_20 * slk_1203[k]
                    + f_3 * pc_z[k] * smk_1527[k];
    }

#pragma omp simd aligned(t_1913, t_1914, t_1915, pb_x, pc_x, sll0_1913, sll0_1914, sll0_1915, \
                         slk_1535, slk_1536, slk_1537, sll1_1913, sll1_1914, \
                         sll1_1915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1913[k] = pb_x[k] * sll0_1913[k]
                    + f_16 * slk_1535[k]
                    - f_14 * pc_x[k] * sll1_1913[k];

        t_1914[k] = pb_x[k] * sll0_1914[k]
                    + f_16 * slk_1536[k]
                    - f_14 * pc_x[k] * sll1_1914[k];

        t_1915[k] = pb_x[k] * sll0_1915[k]
                    + f_16 * slk_1537[k]
                    - f_14 * pc_x[k] * sll1_1915[k];
    }

#pragma omp simd aligned(t_1916, t_1917, t_1918, t_1919, pb_x, pc_x, pc_y, sll0_1917, \
                         slk_1244, slk_1539, slk_1540, slk_1541, sll1_1917, smk_1532, \
                         smk_1540, smk_1541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1916[k] = f_16 * slk_1244[k]
                    + f_3 * pc_y[k] * smk_1532[k];

        t_1917[k] = pb_x[k] * sll0_1917[k]
                    + f_16 * slk_1539[k]
                    - f_14 * pc_x[k] * sll1_1917[k];

        t_1918[k] = f_15 * slk_1540[k]
                    + f_3 * pc_x[k] * smk_1540[k];

        t_1919[k] = f_15 * slk_1541[k]
                    + f_3 * pc_x[k] * smk_1541[k];
    }

#pragma omp simd aligned(t_1920, t_1921, t_1922, t_1923, t_1924, pc_x, slk_1542, slk_1543, \
                         slk_1544, slk_1545, slk_1546, smk_1542, smk_1543, smk_1544, smk_1545, \
                         smk_1546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1920[k] = f_15 * slk_1542[k]
                    + f_3 * pc_x[k] * smk_1542[k];

        t_1921[k] = f_15 * slk_1543[k]
                    + f_3 * pc_x[k] * smk_1543[k];

        t_1922[k] = f_15 * slk_1544[k]
                    + f_3 * pc_x[k] * smk_1544[k];

        t_1923[k] = f_15 * slk_1545[k]
                    + f_3 * pc_x[k] * smk_1545[k];

        t_1924[k] = f_15 * slk_1546[k]
                    + f_3 * pc_x[k] * smk_1546[k];
    }

#pragma omp simd aligned(t_1925, t_1926, t_1927, t_1928, pb_x, pc_x, pc_z, sll0_1926, \
                         sll0_1928, slk_1216, slk_1547, sll1_1926, sll1_1928, smk_1540, \
                         smk_1547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1925[k] = f_15 * slk_1547[k]
                    + f_3 * pc_x[k] * smk_1547[k];

        t_1926[k] = pb_x[k] * sll0_1926[k]
                    - f_14 * pc_x[k] * sll1_1926[k];

        t_1927[k] = f_20 * slk_1216[k]
                    + f_3 * pc_z[k] * smk_1540[k];

        t_1928[k] = pb_x[k] * sll0_1928[k]
                    - f_14 * pc_x[k] * sll1_1928[k];
    }

#pragma omp simd aligned(t_1929, t_1930, t_1931, t_1932, pb_x, pc_x, sll0_1929, sll0_1930, \
                         sll0_1931, sll0_1932, sll1_1929, sll1_1930, sll1_1931, \
                         sll1_1932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1929[k] = pb_x[k] * sll0_1929[k]
                    - f_14 * pc_x[k] * sll1_1929[k];

        t_1930[k] = pb_x[k] * sll0_1930[k]
                    - f_14 * pc_x[k] * sll1_1930[k];

        t_1931[k] = pb_x[k] * sll0_1931[k]
                    - f_14 * pc_x[k] * sll1_1931[k];

        t_1932[k] = pb_x[k] * sll0_1932[k]
                    - f_14 * pc_x[k] * sll1_1932[k];
    }

#pragma omp simd aligned(t_1933, t_1934, t_1935, t_1936, pb_x, pb_y, pc_x, pc_y, sll0_1575, \
                         sll0_1934, slk_1259, slk_1260, sll1_1575, sll1_1934, smk_1547, \
                         smk_1548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1933[k] = f_16 * slk_1259[k]
                    + f_3 * pc_y[k] * smk_1547[k];

        t_1934[k] = pb_x[k] * sll0_1934[k]
                    - f_14 * pc_x[k] * sll1_1934[k];

        t_1935[k] = pb_y[k] * sll0_1575[k]
                    - f_14 * pc_y[k] * sll1_1575[k];

        t_1936[k] = f_15 * slk_1260[k]
                    + f_3 * pc_y[k] * smk_1548[k];
    }

#pragma omp simd aligned(t_1937, t_1938, t_1939, pb_x, pc_x, pc_y, pc_z, sll0_1938, slk_1224, \
                         slk_1262, slk_1551, sll1_1938, smk_1548, \
                         smk_1550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1937[k] = f_22 * slk_1224[k]
                    + f_3 * pc_z[k] * smk_1548[k];

        t_1938[k] = pb_x[k] * sll0_1938[k]
                    + f_20 * slk_1551[k]
                    - f_14 * pc_x[k] * sll1_1938[k];

        t_1939[k] = f_15 * slk_1262[k]
                    + f_3 * pc_y[k] * smk_1550[k];
    }

#pragma omp simd aligned(t_1940, t_1941, t_1942, pb_x, pb_y, pc_x, pc_y, pc_z, sll0_1580, \
                         sll0_1941, slk_1227, slk_1554, sll1_1580, sll1_1941, \
                         smk_1551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1940[k] = pb_y[k] * sll0_1580[k]
                    - f_14 * pc_y[k] * sll1_1580[k];

        t_1941[k] = pb_x[k] * sll0_1941[k]
                    + f_19 * slk_1554[k]
                    - f_14 * pc_x[k] * sll1_1941[k];

        t_1942[k] = f_22 * slk_1227[k]
                    + f_3 * pc_z[k] * smk_1551[k];
    }

#pragma omp simd aligned(t_1943, t_1944, t_1945, pb_x, pb_y, pc_x, pc_y, sll0_1584, sll0_1945, \
                         slk_1265, slk_1558, sll1_1584, sll1_1945, \
                         smk_1553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1943[k] = f_15 * slk_1265[k]
                    + f_3 * pc_y[k] * smk_1553[k];

        t_1944[k] = pb_y[k] * sll0_1584[k]
                    - f_14 * pc_y[k] * sll1_1584[k];

        t_1945[k] = pb_x[k] * sll0_1945[k]
                    + f_18 * slk_1558[k]
                    - f_14 * pc_x[k] * sll1_1945[k];
    }

#pragma omp simd aligned(t_1946, t_1947, t_1948, pb_x, pc_x, pc_y, pc_z, sll0_1947, slk_1230, \
                         slk_1269, slk_1560, sll1_1947, smk_1554, \
                         smk_1557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1946[k] = f_22 * slk_1230[k]
                    + f_3 * pc_z[k] * smk_1554[k];

        t_1947[k] = pb_x[k] * sll0_1947[k]
                    + f_18 * slk_1560[k]
                    - f_14 * pc_x[k] * sll1_1947[k];

        t_1948[k] = f_15 * slk_1269[k]
                    + f_3 * pc_y[k] * smk_1557[k];
    }

#pragma omp simd aligned(t_1949, t_1950, t_1951, pb_x, pb_y, pc_x, pc_y, pc_z, sll0_1589, \
                         sll0_1950, slk_1234, slk_1563, sll1_1589, sll1_1950, \
                         smk_1558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1949[k] = pb_y[k] * sll0_1589[k]
                    - f_14 * pc_y[k] * sll1_1589[k];

        t_1950[k] = pb_x[k] * sll0_1950[k]
                    + f_17 * slk_1563[k]
                    - f_14 * pc_x[k] * sll1_1950[k];

        t_1951[k] = f_22 * slk_1234[k]
                    + f_3 * pc_z[k] * smk_1558[k];
    }

#pragma omp simd aligned(t_1952, t_1953, t_1954, pb_x, pc_x, pc_y, sll0_1952, sll0_1953, \
                         slk_1274, slk_1565, slk_1566, sll1_1952, sll1_1953, \
                         smk_1562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1952[k] = pb_x[k] * sll0_1952[k]
                    + f_17 * slk_1565[k]
                    - f_14 * pc_x[k] * sll1_1952[k];

        t_1953[k] = pb_x[k] * sll0_1953[k]
                    + f_17 * slk_1566[k]
                    - f_14 * pc_x[k] * sll1_1953[k];

        t_1954[k] = f_15 * slk_1274[k]
                    + f_3 * pc_y[k] * smk_1562[k];
    }

#pragma omp simd aligned(t_1955, t_1956, t_1957, pb_x, pb_y, pc_x, pc_y, pc_z, sll0_1595, \
                         sll0_1956, slk_1239, slk_1569, sll1_1595, sll1_1956, \
                         smk_1563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1955[k] = pb_y[k] * sll0_1595[k]
                    - f_14 * pc_y[k] * sll1_1595[k];

        t_1956[k] = pb_x[k] * sll0_1956[k]
                    + f_16 * slk_1569[k]
                    - f_14 * pc_x[k] * sll1_1956[k];

        t_1957[k] = f_22 * slk_1239[k]
                    + f_3 * pc_z[k] * smk_1563[k];
    }

#pragma omp simd aligned(t_1958, t_1959, t_1960, pb_x, pc_x, sll0_1958, sll0_1959, sll0_1960, \
                         slk_1571, slk_1572, slk_1573, sll1_1958, sll1_1959, \
                         sll1_1960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1958[k] = pb_x[k] * sll0_1958[k]
                    + f_16 * slk_1571[k]
                    - f_14 * pc_x[k] * sll1_1958[k];

        t_1959[k] = pb_x[k] * sll0_1959[k]
                    + f_16 * slk_1572[k]
                    - f_14 * pc_x[k] * sll1_1959[k];

        t_1960[k] = pb_x[k] * sll0_1960[k]
                    + f_16 * slk_1573[k]
                    - f_14 * pc_x[k] * sll1_1960[k];
    }

#pragma omp simd aligned(t_1961, t_1962, t_1963, t_1964, pb_y, pc_x, pc_y, sll0_1602, \
                         slk_1280, slk_1576, slk_1577, sll1_1602, smk_1568, smk_1576, \
                         smk_1577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1961[k] = f_15 * slk_1280[k]
                    + f_3 * pc_y[k] * smk_1568[k];

        t_1962[k] = pb_y[k] * sll0_1602[k]
                    - f_14 * pc_y[k] * sll1_1602[k];

        t_1963[k] = f_15 * slk_1576[k]
                    + f_3 * pc_x[k] * smk_1576[k];

        t_1964[k] = f_15 * slk_1577[k]
                    + f_3 * pc_x[k] * smk_1577[k];
    }

#pragma omp simd aligned(t_1965, t_1966, t_1967, t_1968, t_1969, pc_x, slk_1578, slk_1579, \
                         slk_1580, slk_1581, slk_1582, smk_1578, smk_1579, smk_1580, smk_1581, \
                         smk_1582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1965[k] = f_15 * slk_1578[k]
                    + f_3 * pc_x[k] * smk_1578[k];

        t_1966[k] = f_15 * slk_1579[k]
                    + f_3 * pc_x[k] * smk_1579[k];

        t_1967[k] = f_15 * slk_1580[k]
                    + f_3 * pc_x[k] * smk_1580[k];

        t_1968[k] = f_15 * slk_1581[k]
                    + f_3 * pc_x[k] * smk_1581[k];

        t_1969[k] = f_15 * slk_1582[k]
                    + f_3 * pc_x[k] * smk_1582[k];
    }

#pragma omp simd aligned(t_1970, t_1971, t_1972, t_1973, pb_x, pc_x, pc_z, sll0_1971, \
                         sll0_1973, slk_1252, slk_1583, sll1_1971, sll1_1973, smk_1576, \
                         smk_1583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1970[k] = f_15 * slk_1583[k]
                    + f_3 * pc_x[k] * smk_1583[k];

        t_1971[k] = pb_x[k] * sll0_1971[k]
                    - f_14 * pc_x[k] * sll1_1971[k];

        t_1972[k] = f_22 * slk_1252[k]
                    + f_3 * pc_z[k] * smk_1576[k];

        t_1973[k] = pb_x[k] * sll0_1973[k]
                    - f_14 * pc_x[k] * sll1_1973[k];
    }

#pragma omp simd aligned(t_1974, t_1975, t_1976, t_1977, pb_x, pc_x, sll0_1974, sll0_1975, \
                         sll0_1976, sll0_1977, sll1_1974, sll1_1975, sll1_1976, \
                         sll1_1977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1974[k] = pb_x[k] * sll0_1974[k]
                    - f_14 * pc_x[k] * sll1_1974[k];

        t_1975[k] = pb_x[k] * sll0_1975[k]
                    - f_14 * pc_x[k] * sll1_1975[k];

        t_1976[k] = pb_x[k] * sll0_1976[k]
                    - f_14 * pc_x[k] * sll1_1976[k];

        t_1977[k] = pb_x[k] * sll0_1977[k]
                    - f_14 * pc_x[k] * sll1_1977[k];
    }

#pragma omp simd aligned(t_1978, t_1979, t_1980, t_1981, pb_x, pc_x, pc_y, sll0_1979, \
                         sll0_1980, slk_1295, slk_1584, sll1_1979, sll1_1980, smk_1583, \
                         smk_1584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1978[k] = f_15 * slk_1295[k]
                    + f_3 * pc_y[k] * smk_1583[k];

        t_1979[k] = pb_x[k] * sll0_1979[k]
                    - f_14 * pc_x[k] * sll1_1979[k];

        t_1980[k] = pb_x[k] * sll0_1980[k]
                    + f_21 * slk_1584[k]
                    - f_14 * pc_x[k] * sll1_1980[k];

        t_1981[k] = f_3 * pc_y[k] * smk_1584[k];
    }

#pragma omp simd aligned(t_1982, t_1983, t_1984, pb_x, pc_x, pc_y, pc_z, sll0_1983, slk_1260, \
                         slk_1587, sll1_1983, smk_1584, smk_1586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1982[k] = f_21 * slk_1260[k]
                    + f_3 * pc_z[k] * smk_1584[k];

        t_1983[k] = pb_x[k] * sll0_1983[k]
                    + f_20 * slk_1587[k]
                    - f_14 * pc_x[k] * sll1_1983[k];

        t_1984[k] = f_3 * pc_y[k] * smk_1586[k];
    }

#pragma omp simd aligned(t_1985, t_1986, t_1987, pb_x, pc_x, pc_z, sll0_1985, sll0_1986, \
                         slk_1263, slk_1589, slk_1590, sll1_1985, sll1_1986, \
                         smk_1587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1985[k] = pb_x[k] * sll0_1985[k]
                    + f_20 * slk_1589[k]
                    - f_14 * pc_x[k] * sll1_1985[k];

        t_1986[k] = pb_x[k] * sll0_1986[k]
                    + f_19 * slk_1590[k]
                    - f_14 * pc_x[k] * sll1_1986[k];

        t_1987[k] = f_21 * slk_1263[k]
                    + f_3 * pc_z[k] * smk_1587[k];
    }

#pragma omp simd aligned(t_1988, t_1989, t_1990, pb_x, pc_x, pc_y, sll0_1989, sll0_1990, \
                         slk_1593, slk_1594, sll1_1989, sll1_1990, \
                         smk_1589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1988[k] = f_3 * pc_y[k] * smk_1589[k];

        t_1989[k] = pb_x[k] * sll0_1989[k]
                    + f_19 * slk_1593[k]
                    - f_14 * pc_x[k] * sll1_1989[k];

        t_1990[k] = pb_x[k] * sll0_1990[k]
                    + f_18 * slk_1594[k]
                    - f_14 * pc_x[k] * sll1_1990[k];
    }

#pragma omp simd aligned(t_1991, t_1992, t_1993, pb_x, pc_x, pc_y, pc_z, sll0_1992, slk_1266, \
                         slk_1596, sll1_1992, smk_1590, smk_1593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1991[k] = f_21 * slk_1266[k]
                    + f_3 * pc_z[k] * smk_1590[k];

        t_1992[k] = pb_x[k] * sll0_1992[k]
                    + f_18 * slk_1596[k]
                    - f_14 * pc_x[k] * sll1_1992[k];

        t_1993[k] = f_3 * pc_y[k] * smk_1593[k];
    }

#pragma omp simd aligned(t_1994, t_1995, t_1996, pb_x, pc_x, pc_z, sll0_1994, sll0_1995, \
                         slk_1270, slk_1598, slk_1599, sll1_1994, sll1_1995, \
                         smk_1594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1994[k] = pb_x[k] * sll0_1994[k]
                    + f_18 * slk_1598[k]
                    - f_14 * pc_x[k] * sll1_1994[k];

        t_1995[k] = pb_x[k] * sll0_1995[k]
                    + f_17 * slk_1599[k]
                    - f_14 * pc_x[k] * sll1_1995[k];

        t_1996[k] = f_21 * slk_1270[k]
                    + f_3 * pc_z[k] * smk_1594[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece18(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smi0, const size_t smi1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_1620 = buffer.data(sll0 + 1620);
    const auto *sll0_1623 = buffer.data(sll0 + 1623);
    const auto *sll0_1626 = buffer.data(sll0 + 1626);
    const auto *sll0_1630 = buffer.data(sll0 + 1630);
    const auto *sll0_1635 = buffer.data(sll0 + 1635);
    const auto *sll0_1641 = buffer.data(sll0 + 1641);
    const auto *sll0_1656 = buffer.data(sll0 + 1656);
    const auto *sll0_1658 = buffer.data(sll0 + 1658);
    const auto *sll0_1659 = buffer.data(sll0 + 1659);
    const auto *sll0_1660 = buffer.data(sll0 + 1660);
    const auto *sll0_1661 = buffer.data(sll0 + 1661);
    const auto *sll0_1662 = buffer.data(sll0 + 1662);
    const auto *sll0_1997 = buffer.data(sll0 + 1997);
    const auto *sll0_1998 = buffer.data(sll0 + 1998);
    const auto *sll0_2000 = buffer.data(sll0 + 2000);
    const auto *sll0_2001 = buffer.data(sll0 + 2001);
    const auto *sll0_2003 = buffer.data(sll0 + 2003);
    const auto *sll0_2004 = buffer.data(sll0 + 2004);
    const auto *sll0_2005 = buffer.data(sll0 + 2005);
    const auto *sll0_2007 = buffer.data(sll0 + 2007);
    const auto *sll0_2016 = buffer.data(sll0 + 2016);
    const auto *sll0_2018 = buffer.data(sll0 + 2018);
    const auto *sll0_2019 = buffer.data(sll0 + 2019);
    const auto *sll0_2020 = buffer.data(sll0 + 2020);
    const auto *sll0_2021 = buffer.data(sll0 + 2021);
    const auto *sll0_2022 = buffer.data(sll0 + 2022);
    const auto *sll0_2024 = buffer.data(sll0 + 2024);

    const auto *slk_1275 = buffer.data(slk + 1275);
    const auto *slk_1288 = buffer.data(slk + 1288);
    const auto *slk_1296 = buffer.data(slk + 1296);
    const auto *slk_1298 = buffer.data(slk + 1298);
    const auto *slk_1299 = buffer.data(slk + 1299);
    const auto *slk_1301 = buffer.data(slk + 1301);
    const auto *slk_1302 = buffer.data(slk + 1302);
    const auto *slk_1305 = buffer.data(slk + 1305);
    const auto *slk_1306 = buffer.data(slk + 1306);
    const auto *slk_1310 = buffer.data(slk + 1310);
    const auto *slk_1311 = buffer.data(slk + 1311);
    const auto *slk_1316 = buffer.data(slk + 1316);
    const auto *slk_1324 = buffer.data(slk + 1324);
    const auto *slk_1325 = buffer.data(slk + 1325);
    const auto *slk_1326 = buffer.data(slk + 1326);
    const auto *slk_1327 = buffer.data(slk + 1327);
    const auto *slk_1328 = buffer.data(slk + 1328);
    const auto *slk_1329 = buffer.data(slk + 1329);
    const auto *slk_1330 = buffer.data(slk + 1330);
    const auto *slk_1331 = buffer.data(slk + 1331);
    const auto *slk_1332 = buffer.data(slk + 1332);
    const auto *slk_1334 = buffer.data(slk + 1334);
    const auto *slk_1337 = buffer.data(slk + 1337);
    const auto *slk_1341 = buffer.data(slk + 1341);
    const auto *slk_1346 = buffer.data(slk + 1346);
    const auto *slk_1352 = buffer.data(slk + 1352);
    const auto *slk_1367 = buffer.data(slk + 1367);
    const auto *slk_1368 = buffer.data(slk + 1368);
    const auto *slk_1370 = buffer.data(slk + 1370);
    const auto *slk_1601 = buffer.data(slk + 1601);
    const auto *slk_1602 = buffer.data(slk + 1602);
    const auto *slk_1604 = buffer.data(slk + 1604);
    const auto *slk_1605 = buffer.data(slk + 1605);
    const auto *slk_1607 = buffer.data(slk + 1607);
    const auto *slk_1608 = buffer.data(slk + 1608);
    const auto *slk_1609 = buffer.data(slk + 1609);
    const auto *slk_1611 = buffer.data(slk + 1611);
    const auto *slk_1612 = buffer.data(slk + 1612);
    const auto *slk_1613 = buffer.data(slk + 1613);
    const auto *slk_1614 = buffer.data(slk + 1614);
    const auto *slk_1615 = buffer.data(slk + 1615);
    const auto *slk_1616 = buffer.data(slk + 1616);
    const auto *slk_1617 = buffer.data(slk + 1617);
    const auto *slk_1618 = buffer.data(slk + 1618);
    const auto *slk_1619 = buffer.data(slk + 1619);

    const auto *sll1_1620 = buffer.data(sll1 + 1620);
    const auto *sll1_1623 = buffer.data(sll1 + 1623);
    const auto *sll1_1626 = buffer.data(sll1 + 1626);
    const auto *sll1_1630 = buffer.data(sll1 + 1630);
    const auto *sll1_1635 = buffer.data(sll1 + 1635);
    const auto *sll1_1641 = buffer.data(sll1 + 1641);
    const auto *sll1_1656 = buffer.data(sll1 + 1656);
    const auto *sll1_1658 = buffer.data(sll1 + 1658);
    const auto *sll1_1659 = buffer.data(sll1 + 1659);
    const auto *sll1_1660 = buffer.data(sll1 + 1660);
    const auto *sll1_1661 = buffer.data(sll1 + 1661);
    const auto *sll1_1662 = buffer.data(sll1 + 1662);
    const auto *sll1_1997 = buffer.data(sll1 + 1997);
    const auto *sll1_1998 = buffer.data(sll1 + 1998);
    const auto *sll1_2000 = buffer.data(sll1 + 2000);
    const auto *sll1_2001 = buffer.data(sll1 + 2001);
    const auto *sll1_2003 = buffer.data(sll1 + 2003);
    const auto *sll1_2004 = buffer.data(sll1 + 2004);
    const auto *sll1_2005 = buffer.data(sll1 + 2005);
    const auto *sll1_2007 = buffer.data(sll1 + 2007);
    const auto *sll1_2016 = buffer.data(sll1 + 2016);
    const auto *sll1_2018 = buffer.data(sll1 + 2018);
    const auto *sll1_2019 = buffer.data(sll1 + 2019);
    const auto *sll1_2020 = buffer.data(sll1 + 2020);
    const auto *sll1_2021 = buffer.data(sll1 + 2021);
    const auto *sll1_2022 = buffer.data(sll1 + 2022);
    const auto *sll1_2024 = buffer.data(sll1 + 2024);

    const auto *smi0_1260 = buffer.data(smi0 + 1260);
    const auto *smi0_1263 = buffer.data(smi0 + 1263);
    const auto *smi0_1265 = buffer.data(smi0 + 1265);
    const auto *smi0_1266 = buffer.data(smi0 + 1266);
    const auto *smi0_1269 = buffer.data(smi0 + 1269);
    const auto *smi0_1270 = buffer.data(smi0 + 1270);
    const auto *smi0_1272 = buffer.data(smi0 + 1272);
    const auto *smi0_1274 = buffer.data(smi0 + 1274);
    const auto *smi0_1275 = buffer.data(smi0 + 1275);
    const auto *smi0_1277 = buffer.data(smi0 + 1277);
    const auto *smi0_1278 = buffer.data(smi0 + 1278);
    const auto *smi0_1280 = buffer.data(smi0 + 1280);
    const auto *smi0_1281 = buffer.data(smi0 + 1281);
    const auto *smi0_1283 = buffer.data(smi0 + 1283);
    const auto *smi0_1284 = buffer.data(smi0 + 1284);
    const auto *smi0_1285 = buffer.data(smi0 + 1285);
    const auto *smi0_1286 = buffer.data(smi0 + 1286);
    const auto *smi0_1287 = buffer.data(smi0 + 1287);
    const auto *smi0_1293 = buffer.data(smi0 + 1293);
    const auto *smi0_1297 = buffer.data(smi0 + 1297);
    const auto *smi0_1300 = buffer.data(smi0 + 1300);
    const auto *smi0_1302 = buffer.data(smi0 + 1302);
    const auto *smi0_1305 = buffer.data(smi0 + 1305);
    const auto *smi0_1306 = buffer.data(smi0 + 1306);
    const auto *smi0_1308 = buffer.data(smi0 + 1308);
    const auto *smi0_1311 = buffer.data(smi0 + 1311);
    const auto *smi0_1312 = buffer.data(smi0 + 1312);
    const auto *smi0_1313 = buffer.data(smi0 + 1313);
    const auto *smi0_1315 = buffer.data(smi0 + 1315);
    const auto *smi0_1316 = buffer.data(smi0 + 1316);
    const auto *smi0_1319 = buffer.data(smi0 + 1319);
    const auto *smi0_1321 = buffer.data(smi0 + 1321);

    const auto *smi1_1260 = buffer.data(smi1 + 1260);
    const auto *smi1_1263 = buffer.data(smi1 + 1263);
    const auto *smi1_1265 = buffer.data(smi1 + 1265);
    const auto *smi1_1266 = buffer.data(smi1 + 1266);
    const auto *smi1_1269 = buffer.data(smi1 + 1269);
    const auto *smi1_1270 = buffer.data(smi1 + 1270);
    const auto *smi1_1272 = buffer.data(smi1 + 1272);
    const auto *smi1_1274 = buffer.data(smi1 + 1274);
    const auto *smi1_1275 = buffer.data(smi1 + 1275);
    const auto *smi1_1277 = buffer.data(smi1 + 1277);
    const auto *smi1_1278 = buffer.data(smi1 + 1278);
    const auto *smi1_1280 = buffer.data(smi1 + 1280);
    const auto *smi1_1281 = buffer.data(smi1 + 1281);
    const auto *smi1_1283 = buffer.data(smi1 + 1283);
    const auto *smi1_1284 = buffer.data(smi1 + 1284);
    const auto *smi1_1285 = buffer.data(smi1 + 1285);
    const auto *smi1_1286 = buffer.data(smi1 + 1286);
    const auto *smi1_1287 = buffer.data(smi1 + 1287);
    const auto *smi1_1293 = buffer.data(smi1 + 1293);
    const auto *smi1_1297 = buffer.data(smi1 + 1297);
    const auto *smi1_1300 = buffer.data(smi1 + 1300);
    const auto *smi1_1302 = buffer.data(smi1 + 1302);
    const auto *smi1_1305 = buffer.data(smi1 + 1305);
    const auto *smi1_1306 = buffer.data(smi1 + 1306);
    const auto *smi1_1308 = buffer.data(smi1 + 1308);
    const auto *smi1_1311 = buffer.data(smi1 + 1311);
    const auto *smi1_1312 = buffer.data(smi1 + 1312);
    const auto *smi1_1313 = buffer.data(smi1 + 1313);
    const auto *smi1_1315 = buffer.data(smi1 + 1315);
    const auto *smi1_1316 = buffer.data(smi1 + 1316);
    const auto *smi1_1319 = buffer.data(smi1 + 1319);
    const auto *smi1_1321 = buffer.data(smi1 + 1321);

    const auto *smk_1598 = buffer.data(smk + 1598);
    const auto *smk_1599 = buffer.data(smk + 1599);
    const auto *smk_1604 = buffer.data(smk + 1604);
    const auto *smk_1612 = buffer.data(smk + 1612);
    const auto *smk_1613 = buffer.data(smk + 1613);
    const auto *smk_1614 = buffer.data(smk + 1614);
    const auto *smk_1615 = buffer.data(smk + 1615);
    const auto *smk_1616 = buffer.data(smk + 1616);
    const auto *smk_1617 = buffer.data(smk + 1617);
    const auto *smk_1618 = buffer.data(smk + 1618);
    const auto *smk_1619 = buffer.data(smk + 1619);
    const auto *smk_1620 = buffer.data(smk + 1620);
    const auto *smk_1622 = buffer.data(smk + 1622);
    const auto *smk_1623 = buffer.data(smk + 1623);
    const auto *smk_1625 = buffer.data(smk + 1625);
    const auto *smk_1626 = buffer.data(smk + 1626);
    const auto *smk_1629 = buffer.data(smk + 1629);
    const auto *smk_1630 = buffer.data(smk + 1630);
    const auto *smk_1632 = buffer.data(smk + 1632);
    const auto *smk_1634 = buffer.data(smk + 1634);
    const auto *smk_1635 = buffer.data(smk + 1635);
    const auto *smk_1637 = buffer.data(smk + 1637);
    const auto *smk_1638 = buffer.data(smk + 1638);
    const auto *smk_1640 = buffer.data(smk + 1640);
    const auto *smk_1641 = buffer.data(smk + 1641);
    const auto *smk_1643 = buffer.data(smk + 1643);
    const auto *smk_1644 = buffer.data(smk + 1644);
    const auto *smk_1645 = buffer.data(smk + 1645);
    const auto *smk_1647 = buffer.data(smk + 1647);
    const auto *smk_1648 = buffer.data(smk + 1648);
    const auto *smk_1649 = buffer.data(smk + 1649);
    const auto *smk_1650 = buffer.data(smk + 1650);
    const auto *smk_1651 = buffer.data(smk + 1651);
    const auto *smk_1652 = buffer.data(smk + 1652);
    const auto *smk_1653 = buffer.data(smk + 1653);
    const auto *smk_1654 = buffer.data(smk + 1654);
    const auto *smk_1655 = buffer.data(smk + 1655);
    const auto *smk_1656 = buffer.data(smk + 1656);
    const auto *smk_1658 = buffer.data(smk + 1658);
    const auto *smk_1659 = buffer.data(smk + 1659);
    const auto *smk_1661 = buffer.data(smk + 1661);
    const auto *smk_1662 = buffer.data(smk + 1662);
    const auto *smk_1665 = buffer.data(smk + 1665);
    const auto *smk_1666 = buffer.data(smk + 1666);
    const auto *smk_1668 = buffer.data(smk + 1668);
    const auto *smk_1670 = buffer.data(smk + 1670);
    const auto *smk_1671 = buffer.data(smk + 1671);
    const auto *smk_1673 = buffer.data(smk + 1673);
    const auto *smk_1674 = buffer.data(smk + 1674);
    const auto *smk_1676 = buffer.data(smk + 1676);
    const auto *smk_1679 = buffer.data(smk + 1679);
    const auto *smk_1680 = buffer.data(smk + 1680);
    const auto *smk_1681 = buffer.data(smk + 1681);
    const auto *smk_1683 = buffer.data(smk + 1683);
    const auto *smk_1684 = buffer.data(smk + 1684);
    const auto *smk_1685 = buffer.data(smk + 1685);
    const auto *smk_1686 = buffer.data(smk + 1686);
    const auto *smk_1687 = buffer.data(smk + 1687);
    const auto *smk_1688 = buffer.data(smk + 1688);
    const auto *smk_1689 = buffer.data(smk + 1689);
    const auto *smk_1690 = buffer.data(smk + 1690);
    const auto *smk_1691 = buffer.data(smk + 1691);
    const auto *smk_1692 = buffer.data(smk + 1692);
    const auto *smk_1694 = buffer.data(smk + 1694);
    const auto *smk_1695 = buffer.data(smk + 1695);
    const auto *smk_1697 = buffer.data(smk + 1697);

#pragma omp simd aligned(t_1997, t_1998, t_1999, pb_x, pc_x, pc_y, sll0_1997, sll0_1998, \
                         slk_1601, slk_1602, sll1_1997, sll1_1998, \
                         smk_1598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1997[k] = pb_x[k] * sll0_1997[k]
                    + f_17 * slk_1601[k]
                    - f_14 * pc_x[k] * sll1_1997[k];

        t_1998[k] = pb_x[k] * sll0_1998[k]
                    + f_17 * slk_1602[k]
                    - f_14 * pc_x[k] * sll1_1998[k];

        t_1999[k] = f_3 * pc_y[k] * smk_1598[k];
    }

#pragma omp simd aligned(t_2000, t_2001, t_2002, pb_x, pc_x, pc_z, sll0_2000, sll0_2001, \
                         slk_1275, slk_1604, slk_1605, sll1_2000, sll1_2001, \
                         smk_1599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2000[k] = pb_x[k] * sll0_2000[k]
                    + f_17 * slk_1604[k]
                    - f_14 * pc_x[k] * sll1_2000[k];

        t_2001[k] = pb_x[k] * sll0_2001[k]
                    + f_16 * slk_1605[k]
                    - f_14 * pc_x[k] * sll1_2001[k];

        t_2002[k] = f_21 * slk_1275[k]
                    + f_3 * pc_z[k] * smk_1599[k];
    }

#pragma omp simd aligned(t_2003, t_2004, t_2005, pb_x, pc_x, sll0_2003, sll0_2004, sll0_2005, \
                         slk_1607, slk_1608, slk_1609, sll1_2003, sll1_2004, \
                         sll1_2005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2003[k] = pb_x[k] * sll0_2003[k]
                    + f_16 * slk_1607[k]
                    - f_14 * pc_x[k] * sll1_2003[k];

        t_2004[k] = pb_x[k] * sll0_2004[k]
                    + f_16 * slk_1608[k]
                    - f_14 * pc_x[k] * sll1_2004[k];

        t_2005[k] = pb_x[k] * sll0_2005[k]
                    + f_16 * slk_1609[k]
                    - f_14 * pc_x[k] * sll1_2005[k];
    }

#pragma omp simd aligned(t_2006, t_2007, t_2008, t_2009, pb_x, pc_x, pc_y, sll0_2007, \
                         slk_1611, slk_1612, slk_1613, sll1_2007, smk_1604, smk_1612, \
                         smk_1613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2006[k] = f_3 * pc_y[k] * smk_1604[k];

        t_2007[k] = pb_x[k] * sll0_2007[k]
                    + f_16 * slk_1611[k]
                    - f_14 * pc_x[k] * sll1_2007[k];

        t_2008[k] = f_15 * slk_1612[k]
                    + f_3 * pc_x[k] * smk_1612[k];

        t_2009[k] = f_15 * slk_1613[k]
                    + f_3 * pc_x[k] * smk_1613[k];
    }

#pragma omp simd aligned(t_2010, t_2011, t_2012, t_2013, t_2014, pc_x, slk_1614, slk_1615, \
                         slk_1616, slk_1617, slk_1618, smk_1614, smk_1615, smk_1616, smk_1617, \
                         smk_1618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2010[k] = f_15 * slk_1614[k]
                    + f_3 * pc_x[k] * smk_1614[k];

        t_2011[k] = f_15 * slk_1615[k]
                    + f_3 * pc_x[k] * smk_1615[k];

        t_2012[k] = f_15 * slk_1616[k]
                    + f_3 * pc_x[k] * smk_1616[k];

        t_2013[k] = f_15 * slk_1617[k]
                    + f_3 * pc_x[k] * smk_1617[k];

        t_2014[k] = f_15 * slk_1618[k]
                    + f_3 * pc_x[k] * smk_1618[k];
    }

#pragma omp simd aligned(t_2015, t_2016, t_2017, t_2018, pb_x, pc_x, pc_z, sll0_2016, \
                         sll0_2018, slk_1288, slk_1619, sll1_2016, sll1_2018, smk_1612, \
                         smk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2015[k] = f_15 * slk_1619[k]
                    + f_3 * pc_x[k] * smk_1619[k];

        t_2016[k] = pb_x[k] * sll0_2016[k]
                    - f_14 * pc_x[k] * sll1_2016[k];

        t_2017[k] = f_21 * slk_1288[k]
                    + f_3 * pc_z[k] * smk_1612[k];

        t_2018[k] = pb_x[k] * sll0_2018[k]
                    - f_14 * pc_x[k] * sll1_2018[k];
    }

#pragma omp simd aligned(t_2019, t_2020, t_2021, t_2022, pb_x, pc_x, sll0_2019, sll0_2020, \
                         sll0_2021, sll0_2022, sll1_2019, sll1_2020, sll1_2021, \
                         sll1_2022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2019[k] = pb_x[k] * sll0_2019[k]
                    - f_14 * pc_x[k] * sll1_2019[k];

        t_2020[k] = pb_x[k] * sll0_2020[k]
                    - f_14 * pc_x[k] * sll1_2020[k];

        t_2021[k] = pb_x[k] * sll0_2021[k]
                    - f_14 * pc_x[k] * sll1_2021[k];

        t_2022[k] = pb_x[k] * sll0_2022[k]
                    - f_14 * pc_x[k] * sll1_2022[k];
    }

#pragma omp simd aligned(t_2023, t_2024, t_2025, t_2026, t_2027, pb_x, pc_x, pc_y, pc_z, \
                         sll0_2024, slk_1296, sll1_2024, smi0_1260, smi1_1260, smk_1619, \
                         smk_1620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2023[k] = f_3 * pc_y[k] * smk_1619[k];

        t_2024[k] = pb_x[k] * sll0_2024[k]
                    - f_14 * pc_x[k] * sll1_2024[k];

        t_2025[k] = f_1 * smi0_1260[k]
                    - f_2 * smi1_1260[k]
                    + f_3 * pc_x[k] * smk_1620[k];

        t_2026[k] = f_0 * slk_1296[k]
                    + f_3 * pc_y[k] * smk_1620[k];

        t_2027[k] = f_3 * pc_z[k] * smk_1620[k];
    }

#pragma omp simd aligned(t_2028, t_2029, t_2030, pc_x, pc_y, slk_1298, smi0_1263, smi0_1265, \
                         smi1_1263, smi1_1265, smk_1622, smk_1623, \
                         smk_1625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2028[k] = f_4 * smi0_1263[k]
                    - f_5 * smi1_1263[k]
                    + f_3 * pc_x[k] * smk_1623[k];

        t_2029[k] = f_0 * slk_1298[k]
                    + f_3 * pc_y[k] * smk_1622[k];

        t_2030[k] = f_4 * smi0_1265[k]
                    - f_5 * smi1_1265[k]
                    + f_3 * pc_x[k] * smk_1625[k];
    }

#pragma omp simd aligned(t_2031, t_2032, t_2033, t_2034, pc_x, pc_y, pc_z, slk_1301, \
                         smi0_1266, smi0_1269, smi1_1266, smi1_1269, smk_1623, smk_1625, \
                         smk_1626, smk_1629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2031[k] = f_6 * smi0_1266[k]
                    - f_7 * smi1_1266[k]
                    + f_3 * pc_x[k] * smk_1626[k];

        t_2032[k] = f_3 * pc_z[k] * smk_1623[k];

        t_2033[k] = f_0 * slk_1301[k]
                    + f_3 * pc_y[k] * smk_1625[k];

        t_2034[k] = f_6 * smi0_1269[k]
                    - f_7 * smi1_1269[k]
                    + f_3 * pc_x[k] * smk_1629[k];
    }

#pragma omp simd aligned(t_2035, t_2036, t_2037, t_2038, pc_x, pc_y, pc_z, slk_1305, \
                         smi0_1270, smi0_1272, smi1_1270, smi1_1272, smk_1626, smk_1629, \
                         smk_1630, smk_1632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2035[k] = f_8 * smi0_1270[k]
                    - f_9 * smi1_1270[k]
                    + f_3 * pc_x[k] * smk_1630[k];

        t_2036[k] = f_3 * pc_z[k] * smk_1626[k];

        t_2037[k] = f_8 * smi0_1272[k]
                    - f_9 * smi1_1272[k]
                    + f_3 * pc_x[k] * smk_1632[k];

        t_2038[k] = f_0 * slk_1305[k]
                    + f_3 * pc_y[k] * smk_1629[k];
    }

#pragma omp simd aligned(t_2039, t_2040, t_2041, t_2042, pc_x, pc_z, smi0_1274, smi0_1275, \
                         smi0_1277, smi1_1274, smi1_1275, smi1_1277, smk_1630, smk_1634, \
                         smk_1635, smk_1637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2039[k] = f_8 * smi0_1274[k]
                    - f_9 * smi1_1274[k]
                    + f_3 * pc_x[k] * smk_1634[k];

        t_2040[k] = f_10 * smi0_1275[k]
                    - f_11 * smi1_1275[k]
                    + f_3 * pc_x[k] * smk_1635[k];

        t_2041[k] = f_3 * pc_z[k] * smk_1630[k];

        t_2042[k] = f_10 * smi0_1277[k]
                    - f_11 * smi1_1277[k]
                    + f_3 * pc_x[k] * smk_1637[k];
    }

#pragma omp simd aligned(t_2043, t_2044, t_2045, pc_x, pc_y, slk_1310, smi0_1278, smi0_1280, \
                         smi1_1278, smi1_1280, smk_1634, smk_1638, \
                         smk_1640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2043[k] = f_10 * smi0_1278[k]
                    - f_11 * smi1_1278[k]
                    + f_3 * pc_x[k] * smk_1638[k];

        t_2044[k] = f_0 * slk_1310[k]
                    + f_3 * pc_y[k] * smk_1634[k];

        t_2045[k] = f_10 * smi0_1280[k]
                    - f_11 * smi1_1280[k]
                    + f_3 * pc_x[k] * smk_1640[k];
    }

#pragma omp simd aligned(t_2046, t_2047, t_2048, t_2049, pc_x, pc_z, smi0_1281, smi0_1283, \
                         smi0_1284, smi1_1281, smi1_1283, smi1_1284, smk_1635, smk_1641, \
                         smk_1643, smk_1644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2046[k] = f_12 * smi0_1281[k]
                    - f_13 * smi1_1281[k]
                    + f_3 * pc_x[k] * smk_1641[k];

        t_2047[k] = f_3 * pc_z[k] * smk_1635[k];

        t_2048[k] = f_12 * smi0_1283[k]
                    - f_13 * smi1_1283[k]
                    + f_3 * pc_x[k] * smk_1643[k];

        t_2049[k] = f_12 * smi0_1284[k]
                    - f_13 * smi1_1284[k]
                    + f_3 * pc_x[k] * smk_1644[k];
    }

#pragma omp simd aligned(t_2050, t_2051, t_2052, t_2053, pc_x, pc_y, slk_1316, smi0_1285, \
                         smi0_1287, smi1_1285, smi1_1287, smk_1640, smk_1645, smk_1647, \
                         smk_1648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2050[k] = f_12 * smi0_1285[k]
                    - f_13 * smi1_1285[k]
                    + f_3 * pc_x[k] * smk_1645[k];

        t_2051[k] = f_0 * slk_1316[k]
                    + f_3 * pc_y[k] * smk_1640[k];

        t_2052[k] = f_12 * smi0_1287[k]
                    - f_13 * smi1_1287[k]
                    + f_3 * pc_x[k] * smk_1647[k];

        t_2053[k] = f_3 * pc_x[k] * smk_1648[k];
    }

#pragma omp simd aligned(t_2054, t_2055, t_2056, t_2057, t_2058, t_2059, t_2060, pc_x, \
                         smk_1649, smk_1650, smk_1651, smk_1652, smk_1653, smk_1654, \
                         smk_1655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2054[k] = f_3 * pc_x[k] * smk_1649[k];

        t_2055[k] = f_3 * pc_x[k] * smk_1650[k];

        t_2056[k] = f_3 * pc_x[k] * smk_1651[k];

        t_2057[k] = f_3 * pc_x[k] * smk_1652[k];

        t_2058[k] = f_3 * pc_x[k] * smk_1653[k];

        t_2059[k] = f_3 * pc_x[k] * smk_1654[k];

        t_2060[k] = f_3 * pc_x[k] * smk_1655[k];
    }

#pragma omp simd aligned(t_2061, t_2062, t_2063, pc_y, pc_z, slk_1324, slk_1326, smi0_1281, \
                         smi0_1283, smi1_1281, smi1_1283, smk_1648, \
                         smk_1650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2061[k] = f_0 * slk_1324[k]
                    + f_1 * smi0_1281[k]
                    - f_2 * smi1_1281[k]
                    + f_3 * pc_y[k] * smk_1648[k];

        t_2062[k] = f_3 * pc_z[k] * smk_1648[k];

        t_2063[k] = f_0 * slk_1326[k]
                    + f_4 * smi0_1283[k]
                    - f_5 * smi1_1283[k]
                    + f_3 * pc_y[k] * smk_1650[k];
    }

#pragma omp simd aligned(t_2064, t_2065, t_2066, pc_y, slk_1327, slk_1328, slk_1329, \
                         smi0_1284, smi0_1285, smi0_1286, smi1_1284, smi1_1285, smi1_1286, \
                         smk_1651, smk_1652, smk_1653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2064[k] = f_0 * slk_1327[k]
                    + f_6 * smi0_1284[k]
                    - f_7 * smi1_1284[k]
                    + f_3 * pc_y[k] * smk_1651[k];

        t_2065[k] = f_0 * slk_1328[k]
                    + f_8 * smi0_1285[k]
                    - f_9 * smi1_1285[k]
                    + f_3 * pc_y[k] * smk_1652[k];

        t_2066[k] = f_0 * slk_1329[k]
                    + f_10 * smi0_1286[k]
                    - f_11 * smi1_1286[k]
                    + f_3 * pc_y[k] * smk_1653[k];
    }

#pragma omp simd aligned(t_2067, t_2068, t_2069, t_2070, pb_z, pc_y, pc_z, sll0_1620, \
                         slk_1330, slk_1331, sll1_1620, smi0_1287, smi1_1287, smk_1654, \
                         smk_1655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2067[k] = f_0 * slk_1330[k]
                    + f_12 * smi0_1287[k]
                    - f_13 * smi1_1287[k]
                    + f_3 * pc_y[k] * smk_1654[k];

        t_2068[k] = f_0 * slk_1331[k]
                    + f_3 * pc_y[k] * smk_1655[k];

        t_2069[k] = f_1 * smi0_1287[k]
                    - f_2 * smi1_1287[k]
                    + f_3 * pc_z[k] * smk_1655[k];

        t_2070[k] = pb_z[k] * sll0_1620[k]
                    - f_14 * pc_z[k] * sll1_1620[k];
    }

#pragma omp simd aligned(t_2071, t_2072, t_2073, t_2074, pb_z, pc_y, pc_z, sll0_1623, \
                         slk_1296, slk_1332, slk_1334, sll1_1623, smk_1656, \
                         smk_1658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2071[k] = f_21 * slk_1332[k]
                    + f_3 * pc_y[k] * smk_1656[k];

        t_2072[k] = f_15 * slk_1296[k]
                    + f_3 * pc_z[k] * smk_1656[k];

        t_2073[k] = pb_z[k] * sll0_1623[k]
                    - f_14 * pc_z[k] * sll1_1623[k];

        t_2074[k] = f_21 * slk_1334[k]
                    + f_3 * pc_y[k] * smk_1658[k];
    }

#pragma omp simd aligned(t_2075, t_2076, t_2077, t_2078, pb_z, pc_x, pc_y, pc_z, sll0_1626, \
                         slk_1299, slk_1337, sll1_1626, smi0_1293, smi1_1293, smk_1659, \
                         smk_1661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2075[k] = f_4 * smi0_1293[k]
                    - f_5 * smi1_1293[k]
                    + f_3 * pc_x[k] * smk_1661[k];

        t_2076[k] = pb_z[k] * sll0_1626[k]
                    - f_14 * pc_z[k] * sll1_1626[k];

        t_2077[k] = f_15 * slk_1299[k]
                    + f_3 * pc_z[k] * smk_1659[k];

        t_2078[k] = f_21 * slk_1337[k]
                    + f_3 * pc_y[k] * smk_1661[k];
    }

#pragma omp simd aligned(t_2079, t_2080, t_2081, pb_z, pc_x, pc_z, sll0_1630, slk_1302, \
                         sll1_1630, smi0_1297, smi1_1297, smk_1662, \
                         smk_1665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2079[k] = f_6 * smi0_1297[k]
                    - f_7 * smi1_1297[k]
                    + f_3 * pc_x[k] * smk_1665[k];

        t_2080[k] = pb_z[k] * sll0_1630[k]
                    - f_14 * pc_z[k] * sll1_1630[k];

        t_2081[k] = f_15 * slk_1302[k]
                    + f_3 * pc_z[k] * smk_1662[k];
    }

#pragma omp simd aligned(t_2082, t_2083, t_2084, pc_x, pc_y, slk_1341, smi0_1300, smi0_1302, \
                         smi1_1300, smi1_1302, smk_1665, smk_1668, \
                         smk_1670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2082[k] = f_8 * smi0_1300[k]
                    - f_9 * smi1_1300[k]
                    + f_3 * pc_x[k] * smk_1668[k];

        t_2083[k] = f_21 * slk_1341[k]
                    + f_3 * pc_y[k] * smk_1665[k];

        t_2084[k] = f_8 * smi0_1302[k]
                    - f_9 * smi1_1302[k]
                    + f_3 * pc_x[k] * smk_1670[k];
    }

#pragma omp simd aligned(t_2085, t_2086, t_2087, pb_z, pc_x, pc_z, sll0_1635, slk_1306, \
                         sll1_1635, smi0_1305, smi1_1305, smk_1666, \
                         smk_1673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2085[k] = pb_z[k] * sll0_1635[k]
                    - f_14 * pc_z[k] * sll1_1635[k];

        t_2086[k] = f_15 * slk_1306[k]
                    + f_3 * pc_z[k] * smk_1666[k];

        t_2087[k] = f_10 * smi0_1305[k]
                    - f_11 * smi1_1305[k]
                    + f_3 * pc_x[k] * smk_1673[k];
    }

#pragma omp simd aligned(t_2088, t_2089, t_2090, pc_x, pc_y, slk_1346, smi0_1306, smi0_1308, \
                         smi1_1306, smi1_1308, smk_1670, smk_1674, \
                         smk_1676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2088[k] = f_10 * smi0_1306[k]
                    - f_11 * smi1_1306[k]
                    + f_3 * pc_x[k] * smk_1674[k];

        t_2089[k] = f_21 * slk_1346[k]
                    + f_3 * pc_y[k] * smk_1670[k];

        t_2090[k] = f_10 * smi0_1308[k]
                    - f_11 * smi1_1308[k]
                    + f_3 * pc_x[k] * smk_1676[k];
    }

#pragma omp simd aligned(t_2091, t_2092, t_2093, pb_z, pc_x, pc_z, sll0_1641, slk_1311, \
                         sll1_1641, smi0_1311, smi1_1311, smk_1671, \
                         smk_1679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2091[k] = pb_z[k] * sll0_1641[k]
                    - f_14 * pc_z[k] * sll1_1641[k];

        t_2092[k] = f_15 * slk_1311[k]
                    + f_3 * pc_z[k] * smk_1671[k];

        t_2093[k] = f_12 * smi0_1311[k]
                    - f_13 * smi1_1311[k]
                    + f_3 * pc_x[k] * smk_1679[k];
    }

#pragma omp simd aligned(t_2094, t_2095, t_2096, pc_x, pc_y, slk_1352, smi0_1312, smi0_1313, \
                         smi1_1312, smi1_1313, smk_1676, smk_1680, \
                         smk_1681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2094[k] = f_12 * smi0_1312[k]
                    - f_13 * smi1_1312[k]
                    + f_3 * pc_x[k] * smk_1680[k];

        t_2095[k] = f_12 * smi0_1313[k]
                    - f_13 * smi1_1313[k]
                    + f_3 * pc_x[k] * smk_1681[k];

        t_2096[k] = f_21 * slk_1352[k]
                    + f_3 * pc_y[k] * smk_1676[k];
    }

#pragma omp simd aligned(t_2097, t_2098, t_2099, t_2100, t_2101, t_2102, pc_x, smi0_1315, \
                         smi1_1315, smk_1683, smk_1684, smk_1685, smk_1686, smk_1687, \
                         smk_1688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2097[k] = f_12 * smi0_1315[k]
                    - f_13 * smi1_1315[k]
                    + f_3 * pc_x[k] * smk_1683[k];

        t_2098[k] = f_3 * pc_x[k] * smk_1684[k];

        t_2099[k] = f_3 * pc_x[k] * smk_1685[k];

        t_2100[k] = f_3 * pc_x[k] * smk_1686[k];

        t_2101[k] = f_3 * pc_x[k] * smk_1687[k];

        t_2102[k] = f_3 * pc_x[k] * smk_1688[k];
    }

#pragma omp simd aligned(t_2103, t_2104, t_2105, t_2106, t_2107, pb_z, pc_x, pc_z, sll0_1656, \
                         slk_1324, sll1_1656, smk_1684, smk_1689, smk_1690, \
                         smk_1691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2103[k] = f_3 * pc_x[k] * smk_1689[k];

        t_2104[k] = f_3 * pc_x[k] * smk_1690[k];

        t_2105[k] = f_3 * pc_x[k] * smk_1691[k];

        t_2106[k] = pb_z[k] * sll0_1656[k]
                    - f_14 * pc_z[k] * sll1_1656[k];

        t_2107[k] = f_15 * slk_1324[k]
                    + f_3 * pc_z[k] * smk_1684[k];
    }

#pragma omp simd aligned(t_2108, t_2109, t_2110, pb_z, pc_z, sll0_1658, sll0_1659, sll0_1660, \
                         slk_1325, slk_1326, slk_1327, sll1_1658, sll1_1659, \
                         sll1_1660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2108[k] = pb_z[k] * sll0_1658[k]
                    + f_16 * slk_1325[k]
                    - f_14 * pc_z[k] * sll1_1658[k];

        t_2109[k] = pb_z[k] * sll0_1659[k]
                    + f_17 * slk_1326[k]
                    - f_14 * pc_z[k] * sll1_1659[k];

        t_2110[k] = pb_z[k] * sll0_1660[k]
                    + f_18 * slk_1327[k]
                    - f_14 * pc_z[k] * sll1_1660[k];
    }

#pragma omp simd aligned(t_2111, t_2112, t_2113, pb_z, pc_y, pc_z, sll0_1661, sll0_1662, \
                         slk_1328, slk_1329, slk_1367, sll1_1661, sll1_1662, \
                         smk_1691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2111[k] = pb_z[k] * sll0_1661[k]
                    + f_19 * slk_1328[k]
                    - f_14 * pc_z[k] * sll1_1661[k];

        t_2112[k] = pb_z[k] * sll0_1662[k]
                    + f_20 * slk_1329[k]
                    - f_14 * pc_z[k] * sll1_1662[k];

        t_2113[k] = f_21 * slk_1367[k]
                    + f_3 * pc_y[k] * smk_1691[k];
    }

#pragma omp simd aligned(t_2114, t_2115, t_2116, t_2117, pc_x, pc_y, pc_z, slk_1331, slk_1332, \
                         slk_1368, smi0_1315, smi0_1316, smi1_1315, smi1_1316, smk_1691, \
                         smk_1692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2114[k] = f_15 * slk_1331[k]
                    + f_1 * smi0_1315[k]
                    - f_2 * smi1_1315[k]
                    + f_3 * pc_z[k] * smk_1691[k];

        t_2115[k] = f_1 * smi0_1316[k]
                    - f_2 * smi1_1316[k]
                    + f_3 * pc_x[k] * smk_1692[k];

        t_2116[k] = f_22 * slk_1368[k]
                    + f_3 * pc_y[k] * smk_1692[k];

        t_2117[k] = f_16 * slk_1332[k]
                    + f_3 * pc_z[k] * smk_1692[k];
    }

#pragma omp simd aligned(t_2118, t_2119, t_2120, pc_x, pc_y, slk_1370, smi0_1319, smi0_1321, \
                         smi1_1319, smi1_1321, smk_1694, smk_1695, \
                         smk_1697 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2118[k] = f_4 * smi0_1319[k]
                    - f_5 * smi1_1319[k]
                    + f_3 * pc_x[k] * smk_1695[k];

        t_2119[k] = f_22 * slk_1370[k]
                    + f_3 * pc_y[k] * smk_1694[k];

        t_2120[k] = f_4 * smi0_1321[k]
                    - f_5 * smi1_1321[k]
                    + f_3 * pc_x[k] * smk_1697[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece19(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t slk, const size_t smi0,
                                                           const size_t smi1, const size_t smk,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk_1335 = buffer.data(slk + 1335);
    const auto *slk_1338 = buffer.data(slk + 1338);
    const auto *slk_1342 = buffer.data(slk + 1342);
    const auto *slk_1347 = buffer.data(slk + 1347);
    const auto *slk_1360 = buffer.data(slk + 1360);
    const auto *slk_1367 = buffer.data(slk + 1367);
    const auto *slk_1368 = buffer.data(slk + 1368);
    const auto *slk_1371 = buffer.data(slk + 1371);
    const auto *slk_1373 = buffer.data(slk + 1373);
    const auto *slk_1374 = buffer.data(slk + 1374);
    const auto *slk_1377 = buffer.data(slk + 1377);
    const auto *slk_1378 = buffer.data(slk + 1378);
    const auto *slk_1382 = buffer.data(slk + 1382);
    const auto *slk_1383 = buffer.data(slk + 1383);
    const auto *slk_1388 = buffer.data(slk + 1388);
    const auto *slk_1396 = buffer.data(slk + 1396);
    const auto *slk_1398 = buffer.data(slk + 1398);
    const auto *slk_1399 = buffer.data(slk + 1399);
    const auto *slk_1400 = buffer.data(slk + 1400);
    const auto *slk_1401 = buffer.data(slk + 1401);
    const auto *slk_1402 = buffer.data(slk + 1402);
    const auto *slk_1403 = buffer.data(slk + 1403);
    const auto *slk_1404 = buffer.data(slk + 1404);
    const auto *slk_1406 = buffer.data(slk + 1406);
    const auto *slk_1407 = buffer.data(slk + 1407);
    const auto *slk_1409 = buffer.data(slk + 1409);
    const auto *slk_1410 = buffer.data(slk + 1410);
    const auto *slk_1413 = buffer.data(slk + 1413);
    const auto *slk_1414 = buffer.data(slk + 1414);
    const auto *slk_1418 = buffer.data(slk + 1418);
    const auto *slk_1419 = buffer.data(slk + 1419);
    const auto *slk_1424 = buffer.data(slk + 1424);
    const auto *slk_1432 = buffer.data(slk + 1432);
    const auto *slk_1434 = buffer.data(slk + 1434);
    const auto *slk_1435 = buffer.data(slk + 1435);
    const auto *slk_1436 = buffer.data(slk + 1436);
    const auto *slk_1437 = buffer.data(slk + 1437);
    const auto *slk_1438 = buffer.data(slk + 1438);
    const auto *slk_1439 = buffer.data(slk + 1439);
    const auto *slk_1440 = buffer.data(slk + 1440);
    const auto *slk_1442 = buffer.data(slk + 1442);
    const auto *slk_1445 = buffer.data(slk + 1445);
    const auto *slk_1449 = buffer.data(slk + 1449);
    const auto *slk_1454 = buffer.data(slk + 1454);
    const auto *slk_1460 = buffer.data(slk + 1460);

    const auto *smi0_1322 = buffer.data(smi0 + 1322);
    const auto *smi0_1325 = buffer.data(smi0 + 1325);
    const auto *smi0_1326 = buffer.data(smi0 + 1326);
    const auto *smi0_1328 = buffer.data(smi0 + 1328);
    const auto *smi0_1330 = buffer.data(smi0 + 1330);
    const auto *smi0_1331 = buffer.data(smi0 + 1331);
    const auto *smi0_1333 = buffer.data(smi0 + 1333);
    const auto *smi0_1334 = buffer.data(smi0 + 1334);
    const auto *smi0_1336 = buffer.data(smi0 + 1336);
    const auto *smi0_1337 = buffer.data(smi0 + 1337);
    const auto *smi0_1339 = buffer.data(smi0 + 1339);
    const auto *smi0_1340 = buffer.data(smi0 + 1340);
    const auto *smi0_1341 = buffer.data(smi0 + 1341);
    const auto *smi0_1342 = buffer.data(smi0 + 1342);
    const auto *smi0_1343 = buffer.data(smi0 + 1343);
    const auto *smi0_1344 = buffer.data(smi0 + 1344);
    const auto *smi0_1347 = buffer.data(smi0 + 1347);
    const auto *smi0_1349 = buffer.data(smi0 + 1349);
    const auto *smi0_1350 = buffer.data(smi0 + 1350);
    const auto *smi0_1353 = buffer.data(smi0 + 1353);
    const auto *smi0_1354 = buffer.data(smi0 + 1354);
    const auto *smi0_1356 = buffer.data(smi0 + 1356);
    const auto *smi0_1358 = buffer.data(smi0 + 1358);
    const auto *smi0_1359 = buffer.data(smi0 + 1359);
    const auto *smi0_1361 = buffer.data(smi0 + 1361);
    const auto *smi0_1362 = buffer.data(smi0 + 1362);
    const auto *smi0_1364 = buffer.data(smi0 + 1364);
    const auto *smi0_1365 = buffer.data(smi0 + 1365);
    const auto *smi0_1367 = buffer.data(smi0 + 1367);
    const auto *smi0_1368 = buffer.data(smi0 + 1368);
    const auto *smi0_1369 = buffer.data(smi0 + 1369);
    const auto *smi0_1370 = buffer.data(smi0 + 1370);
    const auto *smi0_1371 = buffer.data(smi0 + 1371);
    const auto *smi0_1372 = buffer.data(smi0 + 1372);
    const auto *smi0_1375 = buffer.data(smi0 + 1375);
    const auto *smi0_1377 = buffer.data(smi0 + 1377);
    const auto *smi0_1378 = buffer.data(smi0 + 1378);
    const auto *smi0_1381 = buffer.data(smi0 + 1381);
    const auto *smi0_1382 = buffer.data(smi0 + 1382);
    const auto *smi0_1384 = buffer.data(smi0 + 1384);
    const auto *smi0_1386 = buffer.data(smi0 + 1386);
    const auto *smi0_1387 = buffer.data(smi0 + 1387);
    const auto *smi0_1389 = buffer.data(smi0 + 1389);
    const auto *smi0_1390 = buffer.data(smi0 + 1390);
    const auto *smi0_1392 = buffer.data(smi0 + 1392);
    const auto *smi0_1393 = buffer.data(smi0 + 1393);
    const auto *smi0_1395 = buffer.data(smi0 + 1395);
    const auto *smi0_1396 = buffer.data(smi0 + 1396);
    const auto *smi0_1397 = buffer.data(smi0 + 1397);
    const auto *smi0_1399 = buffer.data(smi0 + 1399);

    const auto *smi1_1322 = buffer.data(smi1 + 1322);
    const auto *smi1_1325 = buffer.data(smi1 + 1325);
    const auto *smi1_1326 = buffer.data(smi1 + 1326);
    const auto *smi1_1328 = buffer.data(smi1 + 1328);
    const auto *smi1_1330 = buffer.data(smi1 + 1330);
    const auto *smi1_1331 = buffer.data(smi1 + 1331);
    const auto *smi1_1333 = buffer.data(smi1 + 1333);
    const auto *smi1_1334 = buffer.data(smi1 + 1334);
    const auto *smi1_1336 = buffer.data(smi1 + 1336);
    const auto *smi1_1337 = buffer.data(smi1 + 1337);
    const auto *smi1_1339 = buffer.data(smi1 + 1339);
    const auto *smi1_1340 = buffer.data(smi1 + 1340);
    const auto *smi1_1341 = buffer.data(smi1 + 1341);
    const auto *smi1_1342 = buffer.data(smi1 + 1342);
    const auto *smi1_1343 = buffer.data(smi1 + 1343);
    const auto *smi1_1344 = buffer.data(smi1 + 1344);
    const auto *smi1_1347 = buffer.data(smi1 + 1347);
    const auto *smi1_1349 = buffer.data(smi1 + 1349);
    const auto *smi1_1350 = buffer.data(smi1 + 1350);
    const auto *smi1_1353 = buffer.data(smi1 + 1353);
    const auto *smi1_1354 = buffer.data(smi1 + 1354);
    const auto *smi1_1356 = buffer.data(smi1 + 1356);
    const auto *smi1_1358 = buffer.data(smi1 + 1358);
    const auto *smi1_1359 = buffer.data(smi1 + 1359);
    const auto *smi1_1361 = buffer.data(smi1 + 1361);
    const auto *smi1_1362 = buffer.data(smi1 + 1362);
    const auto *smi1_1364 = buffer.data(smi1 + 1364);
    const auto *smi1_1365 = buffer.data(smi1 + 1365);
    const auto *smi1_1367 = buffer.data(smi1 + 1367);
    const auto *smi1_1368 = buffer.data(smi1 + 1368);
    const auto *smi1_1369 = buffer.data(smi1 + 1369);
    const auto *smi1_1370 = buffer.data(smi1 + 1370);
    const auto *smi1_1371 = buffer.data(smi1 + 1371);
    const auto *smi1_1372 = buffer.data(smi1 + 1372);
    const auto *smi1_1375 = buffer.data(smi1 + 1375);
    const auto *smi1_1377 = buffer.data(smi1 + 1377);
    const auto *smi1_1378 = buffer.data(smi1 + 1378);
    const auto *smi1_1381 = buffer.data(smi1 + 1381);
    const auto *smi1_1382 = buffer.data(smi1 + 1382);
    const auto *smi1_1384 = buffer.data(smi1 + 1384);
    const auto *smi1_1386 = buffer.data(smi1 + 1386);
    const auto *smi1_1387 = buffer.data(smi1 + 1387);
    const auto *smi1_1389 = buffer.data(smi1 + 1389);
    const auto *smi1_1390 = buffer.data(smi1 + 1390);
    const auto *smi1_1392 = buffer.data(smi1 + 1392);
    const auto *smi1_1393 = buffer.data(smi1 + 1393);
    const auto *smi1_1395 = buffer.data(smi1 + 1395);
    const auto *smi1_1396 = buffer.data(smi1 + 1396);
    const auto *smi1_1397 = buffer.data(smi1 + 1397);
    const auto *smi1_1399 = buffer.data(smi1 + 1399);

    const auto *smk_1695 = buffer.data(smk + 1695);
    const auto *smk_1697 = buffer.data(smk + 1697);
    const auto *smk_1698 = buffer.data(smk + 1698);
    const auto *smk_1701 = buffer.data(smk + 1701);
    const auto *smk_1702 = buffer.data(smk + 1702);
    const auto *smk_1704 = buffer.data(smk + 1704);
    const auto *smk_1706 = buffer.data(smk + 1706);
    const auto *smk_1707 = buffer.data(smk + 1707);
    const auto *smk_1709 = buffer.data(smk + 1709);
    const auto *smk_1710 = buffer.data(smk + 1710);
    const auto *smk_1712 = buffer.data(smk + 1712);
    const auto *smk_1713 = buffer.data(smk + 1713);
    const auto *smk_1715 = buffer.data(smk + 1715);
    const auto *smk_1716 = buffer.data(smk + 1716);
    const auto *smk_1717 = buffer.data(smk + 1717);
    const auto *smk_1719 = buffer.data(smk + 1719);
    const auto *smk_1720 = buffer.data(smk + 1720);
    const auto *smk_1721 = buffer.data(smk + 1721);
    const auto *smk_1722 = buffer.data(smk + 1722);
    const auto *smk_1723 = buffer.data(smk + 1723);
    const auto *smk_1724 = buffer.data(smk + 1724);
    const auto *smk_1725 = buffer.data(smk + 1725);
    const auto *smk_1726 = buffer.data(smk + 1726);
    const auto *smk_1727 = buffer.data(smk + 1727);
    const auto *smk_1728 = buffer.data(smk + 1728);
    const auto *smk_1730 = buffer.data(smk + 1730);
    const auto *smk_1731 = buffer.data(smk + 1731);
    const auto *smk_1733 = buffer.data(smk + 1733);
    const auto *smk_1734 = buffer.data(smk + 1734);
    const auto *smk_1737 = buffer.data(smk + 1737);
    const auto *smk_1738 = buffer.data(smk + 1738);
    const auto *smk_1740 = buffer.data(smk + 1740);
    const auto *smk_1742 = buffer.data(smk + 1742);
    const auto *smk_1743 = buffer.data(smk + 1743);
    const auto *smk_1745 = buffer.data(smk + 1745);
    const auto *smk_1746 = buffer.data(smk + 1746);
    const auto *smk_1748 = buffer.data(smk + 1748);
    const auto *smk_1749 = buffer.data(smk + 1749);
    const auto *smk_1751 = buffer.data(smk + 1751);
    const auto *smk_1752 = buffer.data(smk + 1752);
    const auto *smk_1753 = buffer.data(smk + 1753);
    const auto *smk_1755 = buffer.data(smk + 1755);
    const auto *smk_1756 = buffer.data(smk + 1756);
    const auto *smk_1757 = buffer.data(smk + 1757);
    const auto *smk_1758 = buffer.data(smk + 1758);
    const auto *smk_1759 = buffer.data(smk + 1759);
    const auto *smk_1760 = buffer.data(smk + 1760);
    const auto *smk_1761 = buffer.data(smk + 1761);
    const auto *smk_1762 = buffer.data(smk + 1762);
    const auto *smk_1763 = buffer.data(smk + 1763);
    const auto *smk_1764 = buffer.data(smk + 1764);
    const auto *smk_1766 = buffer.data(smk + 1766);
    const auto *smk_1767 = buffer.data(smk + 1767);
    const auto *smk_1769 = buffer.data(smk + 1769);
    const auto *smk_1770 = buffer.data(smk + 1770);
    const auto *smk_1773 = buffer.data(smk + 1773);
    const auto *smk_1774 = buffer.data(smk + 1774);
    const auto *smk_1776 = buffer.data(smk + 1776);
    const auto *smk_1778 = buffer.data(smk + 1778);
    const auto *smk_1779 = buffer.data(smk + 1779);
    const auto *smk_1781 = buffer.data(smk + 1781);
    const auto *smk_1782 = buffer.data(smk + 1782);
    const auto *smk_1784 = buffer.data(smk + 1784);
    const auto *smk_1785 = buffer.data(smk + 1785);
    const auto *smk_1787 = buffer.data(smk + 1787);
    const auto *smk_1788 = buffer.data(smk + 1788);
    const auto *smk_1789 = buffer.data(smk + 1789);
    const auto *smk_1791 = buffer.data(smk + 1791);
    const auto *smk_1792 = buffer.data(smk + 1792);
    const auto *smk_1793 = buffer.data(smk + 1793);
    const auto *smk_1794 = buffer.data(smk + 1794);
    const auto *smk_1795 = buffer.data(smk + 1795);
    const auto *smk_1796 = buffer.data(smk + 1796);

#pragma omp simd aligned(t_2121, t_2122, t_2123, pc_x, pc_y, pc_z, slk_1335, slk_1373, \
                         smi0_1322, smi1_1322, smk_1695, smk_1697, \
                         smk_1698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2121[k] = f_6 * smi0_1322[k]
                    - f_7 * smi1_1322[k]
                    + f_3 * pc_x[k] * smk_1698[k];

        t_2122[k] = f_16 * slk_1335[k]
                    + f_3 * pc_z[k] * smk_1695[k];

        t_2123[k] = f_22 * slk_1373[k]
                    + f_3 * pc_y[k] * smk_1697[k];
    }

#pragma omp simd aligned(t_2124, t_2125, t_2126, pc_x, pc_z, slk_1338, smi0_1325, smi0_1326, \
                         smi1_1325, smi1_1326, smk_1698, smk_1701, \
                         smk_1702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2124[k] = f_6 * smi0_1325[k]
                    - f_7 * smi1_1325[k]
                    + f_3 * pc_x[k] * smk_1701[k];

        t_2125[k] = f_8 * smi0_1326[k]
                    - f_9 * smi1_1326[k]
                    + f_3 * pc_x[k] * smk_1702[k];

        t_2126[k] = f_16 * slk_1338[k]
                    + f_3 * pc_z[k] * smk_1698[k];
    }

#pragma omp simd aligned(t_2127, t_2128, t_2129, pc_x, pc_y, slk_1377, smi0_1328, smi0_1330, \
                         smi1_1328, smi1_1330, smk_1701, smk_1704, \
                         smk_1706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2127[k] = f_8 * smi0_1328[k]
                    - f_9 * smi1_1328[k]
                    + f_3 * pc_x[k] * smk_1704[k];

        t_2128[k] = f_22 * slk_1377[k]
                    + f_3 * pc_y[k] * smk_1701[k];

        t_2129[k] = f_8 * smi0_1330[k]
                    - f_9 * smi1_1330[k]
                    + f_3 * pc_x[k] * smk_1706[k];
    }

#pragma omp simd aligned(t_2130, t_2131, t_2132, pc_x, pc_z, slk_1342, smi0_1331, smi0_1333, \
                         smi1_1331, smi1_1333, smk_1702, smk_1707, \
                         smk_1709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2130[k] = f_10 * smi0_1331[k]
                    - f_11 * smi1_1331[k]
                    + f_3 * pc_x[k] * smk_1707[k];

        t_2131[k] = f_16 * slk_1342[k]
                    + f_3 * pc_z[k] * smk_1702[k];

        t_2132[k] = f_10 * smi0_1333[k]
                    - f_11 * smi1_1333[k]
                    + f_3 * pc_x[k] * smk_1709[k];
    }

#pragma omp simd aligned(t_2133, t_2134, t_2135, pc_x, pc_y, slk_1382, smi0_1334, smi0_1336, \
                         smi1_1334, smi1_1336, smk_1706, smk_1710, \
                         smk_1712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2133[k] = f_10 * smi0_1334[k]
                    - f_11 * smi1_1334[k]
                    + f_3 * pc_x[k] * smk_1710[k];

        t_2134[k] = f_22 * slk_1382[k]
                    + f_3 * pc_y[k] * smk_1706[k];

        t_2135[k] = f_10 * smi0_1336[k]
                    - f_11 * smi1_1336[k]
                    + f_3 * pc_x[k] * smk_1712[k];
    }

#pragma omp simd aligned(t_2136, t_2137, t_2138, pc_x, pc_z, slk_1347, smi0_1337, smi0_1339, \
                         smi1_1337, smi1_1339, smk_1707, smk_1713, \
                         smk_1715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2136[k] = f_12 * smi0_1337[k]
                    - f_13 * smi1_1337[k]
                    + f_3 * pc_x[k] * smk_1713[k];

        t_2137[k] = f_16 * slk_1347[k]
                    + f_3 * pc_z[k] * smk_1707[k];

        t_2138[k] = f_12 * smi0_1339[k]
                    - f_13 * smi1_1339[k]
                    + f_3 * pc_x[k] * smk_1715[k];
    }

#pragma omp simd aligned(t_2139, t_2140, t_2141, pc_x, pc_y, slk_1388, smi0_1340, smi0_1341, \
                         smi1_1340, smi1_1341, smk_1712, smk_1716, \
                         smk_1717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2139[k] = f_12 * smi0_1340[k]
                    - f_13 * smi1_1340[k]
                    + f_3 * pc_x[k] * smk_1716[k];

        t_2140[k] = f_12 * smi0_1341[k]
                    - f_13 * smi1_1341[k]
                    + f_3 * pc_x[k] * smk_1717[k];

        t_2141[k] = f_22 * slk_1388[k]
                    + f_3 * pc_y[k] * smk_1712[k];
    }

#pragma omp simd aligned(t_2142, t_2143, t_2144, t_2145, t_2146, t_2147, pc_x, smi0_1343, \
                         smi1_1343, smk_1719, smk_1720, smk_1721, smk_1722, smk_1723, \
                         smk_1724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2142[k] = f_12 * smi0_1343[k]
                    - f_13 * smi1_1343[k]
                    + f_3 * pc_x[k] * smk_1719[k];

        t_2143[k] = f_3 * pc_x[k] * smk_1720[k];

        t_2144[k] = f_3 * pc_x[k] * smk_1721[k];

        t_2145[k] = f_3 * pc_x[k] * smk_1722[k];

        t_2146[k] = f_3 * pc_x[k] * smk_1723[k];

        t_2147[k] = f_3 * pc_x[k] * smk_1724[k];
    }

#pragma omp simd aligned(t_2148, t_2149, t_2150, t_2151, t_2152, pc_x, pc_y, pc_z, slk_1360, \
                         slk_1396, smi0_1337, smi1_1337, smk_1720, smk_1725, smk_1726, \
                         smk_1727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2148[k] = f_3 * pc_x[k] * smk_1725[k];

        t_2149[k] = f_3 * pc_x[k] * smk_1726[k];

        t_2150[k] = f_3 * pc_x[k] * smk_1727[k];

        t_2151[k] = f_22 * slk_1396[k]
                    + f_1 * smi0_1337[k]
                    - f_2 * smi1_1337[k]
                    + f_3 * pc_y[k] * smk_1720[k];

        t_2152[k] = f_16 * slk_1360[k]
                    + f_3 * pc_z[k] * smk_1720[k];
    }

#pragma omp simd aligned(t_2153, t_2154, t_2155, pc_y, slk_1398, slk_1399, slk_1400, \
                         smi0_1339, smi0_1340, smi0_1341, smi1_1339, smi1_1340, smi1_1341, \
                         smk_1722, smk_1723, smk_1724 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2153[k] = f_22 * slk_1398[k]
                    + f_4 * smi0_1339[k]
                    - f_5 * smi1_1339[k]
                    + f_3 * pc_y[k] * smk_1722[k];

        t_2154[k] = f_22 * slk_1399[k]
                    + f_6 * smi0_1340[k]
                    - f_7 * smi1_1340[k]
                    + f_3 * pc_y[k] * smk_1723[k];

        t_2155[k] = f_22 * slk_1400[k]
                    + f_8 * smi0_1341[k]
                    - f_9 * smi1_1341[k]
                    + f_3 * pc_y[k] * smk_1724[k];
    }

#pragma omp simd aligned(t_2156, t_2157, t_2158, pc_y, slk_1401, slk_1402, slk_1403, \
                         smi0_1342, smi0_1343, smi1_1342, smi1_1343, smk_1725, smk_1726, \
                         smk_1727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2156[k] = f_22 * slk_1401[k]
                    + f_10 * smi0_1342[k]
                    - f_11 * smi1_1342[k]
                    + f_3 * pc_y[k] * smk_1725[k];

        t_2157[k] = f_22 * slk_1402[k]
                    + f_12 * smi0_1343[k]
                    - f_13 * smi1_1343[k]
                    + f_3 * pc_y[k] * smk_1726[k];

        t_2158[k] = f_22 * slk_1403[k]
                    + f_3 * pc_y[k] * smk_1727[k];
    }

#pragma omp simd aligned(t_2159, t_2160, t_2161, t_2162, pc_x, pc_y, pc_z, slk_1367, slk_1368, \
                         slk_1404, smi0_1343, smi0_1344, smi1_1343, smi1_1344, smk_1727, \
                         smk_1728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2159[k] = f_16 * slk_1367[k]
                    + f_1 * smi0_1343[k]
                    - f_2 * smi1_1343[k]
                    + f_3 * pc_z[k] * smk_1727[k];

        t_2160[k] = f_1 * smi0_1344[k]
                    - f_2 * smi1_1344[k]
                    + f_3 * pc_x[k] * smk_1728[k];

        t_2161[k] = f_20 * slk_1404[k]
                    + f_3 * pc_y[k] * smk_1728[k];

        t_2162[k] = f_17 * slk_1368[k]
                    + f_3 * pc_z[k] * smk_1728[k];
    }

#pragma omp simd aligned(t_2163, t_2164, t_2165, pc_x, pc_y, slk_1406, smi0_1347, smi0_1349, \
                         smi1_1347, smi1_1349, smk_1730, smk_1731, \
                         smk_1733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2163[k] = f_4 * smi0_1347[k]
                    - f_5 * smi1_1347[k]
                    + f_3 * pc_x[k] * smk_1731[k];

        t_2164[k] = f_20 * slk_1406[k]
                    + f_3 * pc_y[k] * smk_1730[k];

        t_2165[k] = f_4 * smi0_1349[k]
                    - f_5 * smi1_1349[k]
                    + f_3 * pc_x[k] * smk_1733[k];
    }

#pragma omp simd aligned(t_2166, t_2167, t_2168, pc_x, pc_y, pc_z, slk_1371, slk_1409, \
                         smi0_1350, smi1_1350, smk_1731, smk_1733, \
                         smk_1734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2166[k] = f_6 * smi0_1350[k]
                    - f_7 * smi1_1350[k]
                    + f_3 * pc_x[k] * smk_1734[k];

        t_2167[k] = f_17 * slk_1371[k]
                    + f_3 * pc_z[k] * smk_1731[k];

        t_2168[k] = f_20 * slk_1409[k]
                    + f_3 * pc_y[k] * smk_1733[k];
    }

#pragma omp simd aligned(t_2169, t_2170, t_2171, pc_x, pc_z, slk_1374, smi0_1353, smi0_1354, \
                         smi1_1353, smi1_1354, smk_1734, smk_1737, \
                         smk_1738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2169[k] = f_6 * smi0_1353[k]
                    - f_7 * smi1_1353[k]
                    + f_3 * pc_x[k] * smk_1737[k];

        t_2170[k] = f_8 * smi0_1354[k]
                    - f_9 * smi1_1354[k]
                    + f_3 * pc_x[k] * smk_1738[k];

        t_2171[k] = f_17 * slk_1374[k]
                    + f_3 * pc_z[k] * smk_1734[k];
    }

#pragma omp simd aligned(t_2172, t_2173, t_2174, pc_x, pc_y, slk_1413, smi0_1356, smi0_1358, \
                         smi1_1356, smi1_1358, smk_1737, smk_1740, \
                         smk_1742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2172[k] = f_8 * smi0_1356[k]
                    - f_9 * smi1_1356[k]
                    + f_3 * pc_x[k] * smk_1740[k];

        t_2173[k] = f_20 * slk_1413[k]
                    + f_3 * pc_y[k] * smk_1737[k];

        t_2174[k] = f_8 * smi0_1358[k]
                    - f_9 * smi1_1358[k]
                    + f_3 * pc_x[k] * smk_1742[k];
    }

#pragma omp simd aligned(t_2175, t_2176, t_2177, pc_x, pc_z, slk_1378, smi0_1359, smi0_1361, \
                         smi1_1359, smi1_1361, smk_1738, smk_1743, \
                         smk_1745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2175[k] = f_10 * smi0_1359[k]
                    - f_11 * smi1_1359[k]
                    + f_3 * pc_x[k] * smk_1743[k];

        t_2176[k] = f_17 * slk_1378[k]
                    + f_3 * pc_z[k] * smk_1738[k];

        t_2177[k] = f_10 * smi0_1361[k]
                    - f_11 * smi1_1361[k]
                    + f_3 * pc_x[k] * smk_1745[k];
    }

#pragma omp simd aligned(t_2178, t_2179, t_2180, pc_x, pc_y, slk_1418, smi0_1362, smi0_1364, \
                         smi1_1362, smi1_1364, smk_1742, smk_1746, \
                         smk_1748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2178[k] = f_10 * smi0_1362[k]
                    - f_11 * smi1_1362[k]
                    + f_3 * pc_x[k] * smk_1746[k];

        t_2179[k] = f_20 * slk_1418[k]
                    + f_3 * pc_y[k] * smk_1742[k];

        t_2180[k] = f_10 * smi0_1364[k]
                    - f_11 * smi1_1364[k]
                    + f_3 * pc_x[k] * smk_1748[k];
    }

#pragma omp simd aligned(t_2181, t_2182, t_2183, pc_x, pc_z, slk_1383, smi0_1365, smi0_1367, \
                         smi1_1365, smi1_1367, smk_1743, smk_1749, \
                         smk_1751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2181[k] = f_12 * smi0_1365[k]
                    - f_13 * smi1_1365[k]
                    + f_3 * pc_x[k] * smk_1749[k];

        t_2182[k] = f_17 * slk_1383[k]
                    + f_3 * pc_z[k] * smk_1743[k];

        t_2183[k] = f_12 * smi0_1367[k]
                    - f_13 * smi1_1367[k]
                    + f_3 * pc_x[k] * smk_1751[k];
    }

#pragma omp simd aligned(t_2184, t_2185, t_2186, pc_x, pc_y, slk_1424, smi0_1368, smi0_1369, \
                         smi1_1368, smi1_1369, smk_1748, smk_1752, \
                         smk_1753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2184[k] = f_12 * smi0_1368[k]
                    - f_13 * smi1_1368[k]
                    + f_3 * pc_x[k] * smk_1752[k];

        t_2185[k] = f_12 * smi0_1369[k]
                    - f_13 * smi1_1369[k]
                    + f_3 * pc_x[k] * smk_1753[k];

        t_2186[k] = f_20 * slk_1424[k]
                    + f_3 * pc_y[k] * smk_1748[k];
    }

#pragma omp simd aligned(t_2187, t_2188, t_2189, t_2190, t_2191, t_2192, pc_x, smi0_1371, \
                         smi1_1371, smk_1755, smk_1756, smk_1757, smk_1758, smk_1759, \
                         smk_1760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2187[k] = f_12 * smi0_1371[k]
                    - f_13 * smi1_1371[k]
                    + f_3 * pc_x[k] * smk_1755[k];

        t_2188[k] = f_3 * pc_x[k] * smk_1756[k];

        t_2189[k] = f_3 * pc_x[k] * smk_1757[k];

        t_2190[k] = f_3 * pc_x[k] * smk_1758[k];

        t_2191[k] = f_3 * pc_x[k] * smk_1759[k];

        t_2192[k] = f_3 * pc_x[k] * smk_1760[k];
    }

#pragma omp simd aligned(t_2193, t_2194, t_2195, t_2196, t_2197, pc_x, pc_y, pc_z, slk_1396, \
                         slk_1432, smi0_1365, smi1_1365, smk_1756, smk_1761, smk_1762, \
                         smk_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2193[k] = f_3 * pc_x[k] * smk_1761[k];

        t_2194[k] = f_3 * pc_x[k] * smk_1762[k];

        t_2195[k] = f_3 * pc_x[k] * smk_1763[k];

        t_2196[k] = f_20 * slk_1432[k]
                    + f_1 * smi0_1365[k]
                    - f_2 * smi1_1365[k]
                    + f_3 * pc_y[k] * smk_1756[k];

        t_2197[k] = f_17 * slk_1396[k]
                    + f_3 * pc_z[k] * smk_1756[k];
    }

#pragma omp simd aligned(t_2198, t_2199, t_2200, pc_y, slk_1434, slk_1435, slk_1436, \
                         smi0_1367, smi0_1368, smi0_1369, smi1_1367, smi1_1368, smi1_1369, \
                         smk_1758, smk_1759, smk_1760 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2198[k] = f_20 * slk_1434[k]
                    + f_4 * smi0_1367[k]
                    - f_5 * smi1_1367[k]
                    + f_3 * pc_y[k] * smk_1758[k];

        t_2199[k] = f_20 * slk_1435[k]
                    + f_6 * smi0_1368[k]
                    - f_7 * smi1_1368[k]
                    + f_3 * pc_y[k] * smk_1759[k];

        t_2200[k] = f_20 * slk_1436[k]
                    + f_8 * smi0_1369[k]
                    - f_9 * smi1_1369[k]
                    + f_3 * pc_y[k] * smk_1760[k];
    }

#pragma omp simd aligned(t_2201, t_2202, t_2203, pc_y, slk_1437, slk_1438, slk_1439, \
                         smi0_1370, smi0_1371, smi1_1370, smi1_1371, smk_1761, smk_1762, \
                         smk_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2201[k] = f_20 * slk_1437[k]
                    + f_10 * smi0_1370[k]
                    - f_11 * smi1_1370[k]
                    + f_3 * pc_y[k] * smk_1761[k];

        t_2202[k] = f_20 * slk_1438[k]
                    + f_12 * smi0_1371[k]
                    - f_13 * smi1_1371[k]
                    + f_3 * pc_y[k] * smk_1762[k];

        t_2203[k] = f_20 * slk_1439[k]
                    + f_3 * pc_y[k] * smk_1763[k];
    }

#pragma omp simd aligned(t_2204, t_2205, t_2206, t_2207, pc_x, pc_y, pc_z, slk_1403, slk_1404, \
                         slk_1440, smi0_1371, smi0_1372, smi1_1371, smi1_1372, smk_1763, \
                         smk_1764 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2204[k] = f_17 * slk_1403[k]
                    + f_1 * smi0_1371[k]
                    - f_2 * smi1_1371[k]
                    + f_3 * pc_z[k] * smk_1763[k];

        t_2205[k] = f_1 * smi0_1372[k]
                    - f_2 * smi1_1372[k]
                    + f_3 * pc_x[k] * smk_1764[k];

        t_2206[k] = f_19 * slk_1440[k]
                    + f_3 * pc_y[k] * smk_1764[k];

        t_2207[k] = f_18 * slk_1404[k]
                    + f_3 * pc_z[k] * smk_1764[k];
    }

#pragma omp simd aligned(t_2208, t_2209, t_2210, pc_x, pc_y, slk_1442, smi0_1375, smi0_1377, \
                         smi1_1375, smi1_1377, smk_1766, smk_1767, \
                         smk_1769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2208[k] = f_4 * smi0_1375[k]
                    - f_5 * smi1_1375[k]
                    + f_3 * pc_x[k] * smk_1767[k];

        t_2209[k] = f_19 * slk_1442[k]
                    + f_3 * pc_y[k] * smk_1766[k];

        t_2210[k] = f_4 * smi0_1377[k]
                    - f_5 * smi1_1377[k]
                    + f_3 * pc_x[k] * smk_1769[k];
    }

#pragma omp simd aligned(t_2211, t_2212, t_2213, pc_x, pc_y, pc_z, slk_1407, slk_1445, \
                         smi0_1378, smi1_1378, smk_1767, smk_1769, \
                         smk_1770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2211[k] = f_6 * smi0_1378[k]
                    - f_7 * smi1_1378[k]
                    + f_3 * pc_x[k] * smk_1770[k];

        t_2212[k] = f_18 * slk_1407[k]
                    + f_3 * pc_z[k] * smk_1767[k];

        t_2213[k] = f_19 * slk_1445[k]
                    + f_3 * pc_y[k] * smk_1769[k];
    }

#pragma omp simd aligned(t_2214, t_2215, t_2216, pc_x, pc_z, slk_1410, smi0_1381, smi0_1382, \
                         smi1_1381, smi1_1382, smk_1770, smk_1773, \
                         smk_1774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2214[k] = f_6 * smi0_1381[k]
                    - f_7 * smi1_1381[k]
                    + f_3 * pc_x[k] * smk_1773[k];

        t_2215[k] = f_8 * smi0_1382[k]
                    - f_9 * smi1_1382[k]
                    + f_3 * pc_x[k] * smk_1774[k];

        t_2216[k] = f_18 * slk_1410[k]
                    + f_3 * pc_z[k] * smk_1770[k];
    }

#pragma omp simd aligned(t_2217, t_2218, t_2219, pc_x, pc_y, slk_1449, smi0_1384, smi0_1386, \
                         smi1_1384, smi1_1386, smk_1773, smk_1776, \
                         smk_1778 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2217[k] = f_8 * smi0_1384[k]
                    - f_9 * smi1_1384[k]
                    + f_3 * pc_x[k] * smk_1776[k];

        t_2218[k] = f_19 * slk_1449[k]
                    + f_3 * pc_y[k] * smk_1773[k];

        t_2219[k] = f_8 * smi0_1386[k]
                    - f_9 * smi1_1386[k]
                    + f_3 * pc_x[k] * smk_1778[k];
    }

#pragma omp simd aligned(t_2220, t_2221, t_2222, pc_x, pc_z, slk_1414, smi0_1387, smi0_1389, \
                         smi1_1387, smi1_1389, smk_1774, smk_1779, \
                         smk_1781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2220[k] = f_10 * smi0_1387[k]
                    - f_11 * smi1_1387[k]
                    + f_3 * pc_x[k] * smk_1779[k];

        t_2221[k] = f_18 * slk_1414[k]
                    + f_3 * pc_z[k] * smk_1774[k];

        t_2222[k] = f_10 * smi0_1389[k]
                    - f_11 * smi1_1389[k]
                    + f_3 * pc_x[k] * smk_1781[k];
    }

#pragma omp simd aligned(t_2223, t_2224, t_2225, pc_x, pc_y, slk_1454, smi0_1390, smi0_1392, \
                         smi1_1390, smi1_1392, smk_1778, smk_1782, \
                         smk_1784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2223[k] = f_10 * smi0_1390[k]
                    - f_11 * smi1_1390[k]
                    + f_3 * pc_x[k] * smk_1782[k];

        t_2224[k] = f_19 * slk_1454[k]
                    + f_3 * pc_y[k] * smk_1778[k];

        t_2225[k] = f_10 * smi0_1392[k]
                    - f_11 * smi1_1392[k]
                    + f_3 * pc_x[k] * smk_1784[k];
    }

#pragma omp simd aligned(t_2226, t_2227, t_2228, pc_x, pc_z, slk_1419, smi0_1393, smi0_1395, \
                         smi1_1393, smi1_1395, smk_1779, smk_1785, \
                         smk_1787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2226[k] = f_12 * smi0_1393[k]
                    - f_13 * smi1_1393[k]
                    + f_3 * pc_x[k] * smk_1785[k];

        t_2227[k] = f_18 * slk_1419[k]
                    + f_3 * pc_z[k] * smk_1779[k];

        t_2228[k] = f_12 * smi0_1395[k]
                    - f_13 * smi1_1395[k]
                    + f_3 * pc_x[k] * smk_1787[k];
    }

#pragma omp simd aligned(t_2229, t_2230, t_2231, pc_x, pc_y, slk_1460, smi0_1396, smi0_1397, \
                         smi1_1396, smi1_1397, smk_1784, smk_1788, \
                         smk_1789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2229[k] = f_12 * smi0_1396[k]
                    - f_13 * smi1_1396[k]
                    + f_3 * pc_x[k] * smk_1788[k];

        t_2230[k] = f_12 * smi0_1397[k]
                    - f_13 * smi1_1397[k]
                    + f_3 * pc_x[k] * smk_1789[k];

        t_2231[k] = f_19 * slk_1460[k]
                    + f_3 * pc_y[k] * smk_1784[k];
    }

#pragma omp simd aligned(t_2232, t_2233, t_2234, t_2235, t_2236, t_2237, pc_x, smi0_1399, \
                         smi1_1399, smk_1791, smk_1792, smk_1793, smk_1794, smk_1795, \
                         smk_1796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2232[k] = f_12 * smi0_1399[k]
                    - f_13 * smi1_1399[k]
                    + f_3 * pc_x[k] * smk_1791[k];

        t_2233[k] = f_3 * pc_x[k] * smk_1792[k];

        t_2234[k] = f_3 * pc_x[k] * smk_1793[k];

        t_2235[k] = f_3 * pc_x[k] * smk_1794[k];

        t_2236[k] = f_3 * pc_x[k] * smk_1795[k];

        t_2237[k] = f_3 * pc_x[k] * smk_1796[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece20(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t slk, const size_t smi0,
                                                           const size_t smi1, const size_t smk,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_22 = 3.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *slk_1432 = buffer.data(slk + 1432);
    const auto *slk_1439 = buffer.data(slk + 1439);
    const auto *slk_1440 = buffer.data(slk + 1440);
    const auto *slk_1443 = buffer.data(slk + 1443);
    const auto *slk_1446 = buffer.data(slk + 1446);
    const auto *slk_1450 = buffer.data(slk + 1450);
    const auto *slk_1455 = buffer.data(slk + 1455);
    const auto *slk_1468 = buffer.data(slk + 1468);
    const auto *slk_1470 = buffer.data(slk + 1470);
    const auto *slk_1471 = buffer.data(slk + 1471);
    const auto *slk_1472 = buffer.data(slk + 1472);
    const auto *slk_1473 = buffer.data(slk + 1473);
    const auto *slk_1474 = buffer.data(slk + 1474);
    const auto *slk_1475 = buffer.data(slk + 1475);
    const auto *slk_1476 = buffer.data(slk + 1476);
    const auto *slk_1478 = buffer.data(slk + 1478);
    const auto *slk_1479 = buffer.data(slk + 1479);
    const auto *slk_1481 = buffer.data(slk + 1481);
    const auto *slk_1482 = buffer.data(slk + 1482);
    const auto *slk_1485 = buffer.data(slk + 1485);
    const auto *slk_1486 = buffer.data(slk + 1486);
    const auto *slk_1490 = buffer.data(slk + 1490);
    const auto *slk_1491 = buffer.data(slk + 1491);
    const auto *slk_1496 = buffer.data(slk + 1496);
    const auto *slk_1504 = buffer.data(slk + 1504);
    const auto *slk_1506 = buffer.data(slk + 1506);
    const auto *slk_1507 = buffer.data(slk + 1507);
    const auto *slk_1508 = buffer.data(slk + 1508);
    const auto *slk_1509 = buffer.data(slk + 1509);
    const auto *slk_1510 = buffer.data(slk + 1510);
    const auto *slk_1511 = buffer.data(slk + 1511);
    const auto *slk_1512 = buffer.data(slk + 1512);
    const auto *slk_1514 = buffer.data(slk + 1514);
    const auto *slk_1515 = buffer.data(slk + 1515);
    const auto *slk_1517 = buffer.data(slk + 1517);
    const auto *slk_1518 = buffer.data(slk + 1518);
    const auto *slk_1521 = buffer.data(slk + 1521);
    const auto *slk_1526 = buffer.data(slk + 1526);
    const auto *slk_1532 = buffer.data(slk + 1532);
    const auto *slk_1540 = buffer.data(slk + 1540);
    const auto *slk_1542 = buffer.data(slk + 1542);
    const auto *slk_1543 = buffer.data(slk + 1543);
    const auto *slk_1544 = buffer.data(slk + 1544);
    const auto *slk_1545 = buffer.data(slk + 1545);
    const auto *slk_1546 = buffer.data(slk + 1546);
    const auto *slk_1547 = buffer.data(slk + 1547);
    const auto *slk_1548 = buffer.data(slk + 1548);
    const auto *slk_1550 = buffer.data(slk + 1550);
    const auto *slk_1553 = buffer.data(slk + 1553);
    const auto *slk_1557 = buffer.data(slk + 1557);

    const auto *smi0_1393 = buffer.data(smi0 + 1393);
    const auto *smi0_1395 = buffer.data(smi0 + 1395);
    const auto *smi0_1396 = buffer.data(smi0 + 1396);
    const auto *smi0_1397 = buffer.data(smi0 + 1397);
    const auto *smi0_1398 = buffer.data(smi0 + 1398);
    const auto *smi0_1399 = buffer.data(smi0 + 1399);
    const auto *smi0_1400 = buffer.data(smi0 + 1400);
    const auto *smi0_1403 = buffer.data(smi0 + 1403);
    const auto *smi0_1405 = buffer.data(smi0 + 1405);
    const auto *smi0_1406 = buffer.data(smi0 + 1406);
    const auto *smi0_1409 = buffer.data(smi0 + 1409);
    const auto *smi0_1410 = buffer.data(smi0 + 1410);
    const auto *smi0_1412 = buffer.data(smi0 + 1412);
    const auto *smi0_1414 = buffer.data(smi0 + 1414);
    const auto *smi0_1415 = buffer.data(smi0 + 1415);
    const auto *smi0_1417 = buffer.data(smi0 + 1417);
    const auto *smi0_1418 = buffer.data(smi0 + 1418);
    const auto *smi0_1420 = buffer.data(smi0 + 1420);
    const auto *smi0_1421 = buffer.data(smi0 + 1421);
    const auto *smi0_1423 = buffer.data(smi0 + 1423);
    const auto *smi0_1424 = buffer.data(smi0 + 1424);
    const auto *smi0_1425 = buffer.data(smi0 + 1425);
    const auto *smi0_1426 = buffer.data(smi0 + 1426);
    const auto *smi0_1427 = buffer.data(smi0 + 1427);
    const auto *smi0_1428 = buffer.data(smi0 + 1428);
    const auto *smi0_1431 = buffer.data(smi0 + 1431);
    const auto *smi0_1433 = buffer.data(smi0 + 1433);
    const auto *smi0_1434 = buffer.data(smi0 + 1434);
    const auto *smi0_1437 = buffer.data(smi0 + 1437);
    const auto *smi0_1438 = buffer.data(smi0 + 1438);
    const auto *smi0_1440 = buffer.data(smi0 + 1440);
    const auto *smi0_1442 = buffer.data(smi0 + 1442);
    const auto *smi0_1443 = buffer.data(smi0 + 1443);
    const auto *smi0_1445 = buffer.data(smi0 + 1445);
    const auto *smi0_1446 = buffer.data(smi0 + 1446);
    const auto *smi0_1448 = buffer.data(smi0 + 1448);
    const auto *smi0_1449 = buffer.data(smi0 + 1449);
    const auto *smi0_1451 = buffer.data(smi0 + 1451);
    const auto *smi0_1452 = buffer.data(smi0 + 1452);
    const auto *smi0_1453 = buffer.data(smi0 + 1453);
    const auto *smi0_1454 = buffer.data(smi0 + 1454);
    const auto *smi0_1455 = buffer.data(smi0 + 1455);
    const auto *smi0_1456 = buffer.data(smi0 + 1456);
    const auto *smi0_1459 = buffer.data(smi0 + 1459);
    const auto *smi0_1461 = buffer.data(smi0 + 1461);
    const auto *smi0_1462 = buffer.data(smi0 + 1462);
    const auto *smi0_1465 = buffer.data(smi0 + 1465);
    const auto *smi0_1466 = buffer.data(smi0 + 1466);
    const auto *smi0_1468 = buffer.data(smi0 + 1468);
    const auto *smi0_1470 = buffer.data(smi0 + 1470);

    const auto *smi1_1393 = buffer.data(smi1 + 1393);
    const auto *smi1_1395 = buffer.data(smi1 + 1395);
    const auto *smi1_1396 = buffer.data(smi1 + 1396);
    const auto *smi1_1397 = buffer.data(smi1 + 1397);
    const auto *smi1_1398 = buffer.data(smi1 + 1398);
    const auto *smi1_1399 = buffer.data(smi1 + 1399);
    const auto *smi1_1400 = buffer.data(smi1 + 1400);
    const auto *smi1_1403 = buffer.data(smi1 + 1403);
    const auto *smi1_1405 = buffer.data(smi1 + 1405);
    const auto *smi1_1406 = buffer.data(smi1 + 1406);
    const auto *smi1_1409 = buffer.data(smi1 + 1409);
    const auto *smi1_1410 = buffer.data(smi1 + 1410);
    const auto *smi1_1412 = buffer.data(smi1 + 1412);
    const auto *smi1_1414 = buffer.data(smi1 + 1414);
    const auto *smi1_1415 = buffer.data(smi1 + 1415);
    const auto *smi1_1417 = buffer.data(smi1 + 1417);
    const auto *smi1_1418 = buffer.data(smi1 + 1418);
    const auto *smi1_1420 = buffer.data(smi1 + 1420);
    const auto *smi1_1421 = buffer.data(smi1 + 1421);
    const auto *smi1_1423 = buffer.data(smi1 + 1423);
    const auto *smi1_1424 = buffer.data(smi1 + 1424);
    const auto *smi1_1425 = buffer.data(smi1 + 1425);
    const auto *smi1_1426 = buffer.data(smi1 + 1426);
    const auto *smi1_1427 = buffer.data(smi1 + 1427);
    const auto *smi1_1428 = buffer.data(smi1 + 1428);
    const auto *smi1_1431 = buffer.data(smi1 + 1431);
    const auto *smi1_1433 = buffer.data(smi1 + 1433);
    const auto *smi1_1434 = buffer.data(smi1 + 1434);
    const auto *smi1_1437 = buffer.data(smi1 + 1437);
    const auto *smi1_1438 = buffer.data(smi1 + 1438);
    const auto *smi1_1440 = buffer.data(smi1 + 1440);
    const auto *smi1_1442 = buffer.data(smi1 + 1442);
    const auto *smi1_1443 = buffer.data(smi1 + 1443);
    const auto *smi1_1445 = buffer.data(smi1 + 1445);
    const auto *smi1_1446 = buffer.data(smi1 + 1446);
    const auto *smi1_1448 = buffer.data(smi1 + 1448);
    const auto *smi1_1449 = buffer.data(smi1 + 1449);
    const auto *smi1_1451 = buffer.data(smi1 + 1451);
    const auto *smi1_1452 = buffer.data(smi1 + 1452);
    const auto *smi1_1453 = buffer.data(smi1 + 1453);
    const auto *smi1_1454 = buffer.data(smi1 + 1454);
    const auto *smi1_1455 = buffer.data(smi1 + 1455);
    const auto *smi1_1456 = buffer.data(smi1 + 1456);
    const auto *smi1_1459 = buffer.data(smi1 + 1459);
    const auto *smi1_1461 = buffer.data(smi1 + 1461);
    const auto *smi1_1462 = buffer.data(smi1 + 1462);
    const auto *smi1_1465 = buffer.data(smi1 + 1465);
    const auto *smi1_1466 = buffer.data(smi1 + 1466);
    const auto *smi1_1468 = buffer.data(smi1 + 1468);
    const auto *smi1_1470 = buffer.data(smi1 + 1470);

    const auto *smk_1792 = buffer.data(smk + 1792);
    const auto *smk_1794 = buffer.data(smk + 1794);
    const auto *smk_1795 = buffer.data(smk + 1795);
    const auto *smk_1796 = buffer.data(smk + 1796);
    const auto *smk_1797 = buffer.data(smk + 1797);
    const auto *smk_1798 = buffer.data(smk + 1798);
    const auto *smk_1799 = buffer.data(smk + 1799);
    const auto *smk_1800 = buffer.data(smk + 1800);
    const auto *smk_1802 = buffer.data(smk + 1802);
    const auto *smk_1803 = buffer.data(smk + 1803);
    const auto *smk_1805 = buffer.data(smk + 1805);
    const auto *smk_1806 = buffer.data(smk + 1806);
    const auto *smk_1809 = buffer.data(smk + 1809);
    const auto *smk_1810 = buffer.data(smk + 1810);
    const auto *smk_1812 = buffer.data(smk + 1812);
    const auto *smk_1814 = buffer.data(smk + 1814);
    const auto *smk_1815 = buffer.data(smk + 1815);
    const auto *smk_1817 = buffer.data(smk + 1817);
    const auto *smk_1818 = buffer.data(smk + 1818);
    const auto *smk_1820 = buffer.data(smk + 1820);
    const auto *smk_1821 = buffer.data(smk + 1821);
    const auto *smk_1823 = buffer.data(smk + 1823);
    const auto *smk_1824 = buffer.data(smk + 1824);
    const auto *smk_1825 = buffer.data(smk + 1825);
    const auto *smk_1827 = buffer.data(smk + 1827);
    const auto *smk_1828 = buffer.data(smk + 1828);
    const auto *smk_1829 = buffer.data(smk + 1829);
    const auto *smk_1830 = buffer.data(smk + 1830);
    const auto *smk_1831 = buffer.data(smk + 1831);
    const auto *smk_1832 = buffer.data(smk + 1832);
    const auto *smk_1833 = buffer.data(smk + 1833);
    const auto *smk_1834 = buffer.data(smk + 1834);
    const auto *smk_1835 = buffer.data(smk + 1835);
    const auto *smk_1836 = buffer.data(smk + 1836);
    const auto *smk_1838 = buffer.data(smk + 1838);
    const auto *smk_1839 = buffer.data(smk + 1839);
    const auto *smk_1841 = buffer.data(smk + 1841);
    const auto *smk_1842 = buffer.data(smk + 1842);
    const auto *smk_1845 = buffer.data(smk + 1845);
    const auto *smk_1846 = buffer.data(smk + 1846);
    const auto *smk_1848 = buffer.data(smk + 1848);
    const auto *smk_1850 = buffer.data(smk + 1850);
    const auto *smk_1851 = buffer.data(smk + 1851);
    const auto *smk_1853 = buffer.data(smk + 1853);
    const auto *smk_1854 = buffer.data(smk + 1854);
    const auto *smk_1856 = buffer.data(smk + 1856);
    const auto *smk_1857 = buffer.data(smk + 1857);
    const auto *smk_1859 = buffer.data(smk + 1859);
    const auto *smk_1860 = buffer.data(smk + 1860);
    const auto *smk_1861 = buffer.data(smk + 1861);
    const auto *smk_1863 = buffer.data(smk + 1863);
    const auto *smk_1864 = buffer.data(smk + 1864);
    const auto *smk_1865 = buffer.data(smk + 1865);
    const auto *smk_1866 = buffer.data(smk + 1866);
    const auto *smk_1867 = buffer.data(smk + 1867);
    const auto *smk_1868 = buffer.data(smk + 1868);
    const auto *smk_1869 = buffer.data(smk + 1869);
    const auto *smk_1870 = buffer.data(smk + 1870);
    const auto *smk_1871 = buffer.data(smk + 1871);
    const auto *smk_1872 = buffer.data(smk + 1872);
    const auto *smk_1874 = buffer.data(smk + 1874);
    const auto *smk_1875 = buffer.data(smk + 1875);
    const auto *smk_1877 = buffer.data(smk + 1877);
    const auto *smk_1878 = buffer.data(smk + 1878);
    const auto *smk_1881 = buffer.data(smk + 1881);
    const auto *smk_1882 = buffer.data(smk + 1882);
    const auto *smk_1884 = buffer.data(smk + 1884);
    const auto *smk_1886 = buffer.data(smk + 1886);

#pragma omp simd aligned(t_2238, t_2239, t_2240, t_2241, t_2242, pc_x, pc_y, pc_z, slk_1432, \
                         slk_1468, smi0_1393, smi1_1393, smk_1792, smk_1797, smk_1798, \
                         smk_1799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2238[k] = f_3 * pc_x[k] * smk_1797[k];

        t_2239[k] = f_3 * pc_x[k] * smk_1798[k];

        t_2240[k] = f_3 * pc_x[k] * smk_1799[k];

        t_2241[k] = f_19 * slk_1468[k]
                    + f_1 * smi0_1393[k]
                    - f_2 * smi1_1393[k]
                    + f_3 * pc_y[k] * smk_1792[k];

        t_2242[k] = f_18 * slk_1432[k]
                    + f_3 * pc_z[k] * smk_1792[k];
    }

#pragma omp simd aligned(t_2243, t_2244, t_2245, pc_y, slk_1470, slk_1471, slk_1472, \
                         smi0_1395, smi0_1396, smi0_1397, smi1_1395, smi1_1396, smi1_1397, \
                         smk_1794, smk_1795, smk_1796 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2243[k] = f_19 * slk_1470[k]
                    + f_4 * smi0_1395[k]
                    - f_5 * smi1_1395[k]
                    + f_3 * pc_y[k] * smk_1794[k];

        t_2244[k] = f_19 * slk_1471[k]
                    + f_6 * smi0_1396[k]
                    - f_7 * smi1_1396[k]
                    + f_3 * pc_y[k] * smk_1795[k];

        t_2245[k] = f_19 * slk_1472[k]
                    + f_8 * smi0_1397[k]
                    - f_9 * smi1_1397[k]
                    + f_3 * pc_y[k] * smk_1796[k];
    }

#pragma omp simd aligned(t_2246, t_2247, t_2248, pc_y, slk_1473, slk_1474, slk_1475, \
                         smi0_1398, smi0_1399, smi1_1398, smi1_1399, smk_1797, smk_1798, \
                         smk_1799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2246[k] = f_19 * slk_1473[k]
                    + f_10 * smi0_1398[k]
                    - f_11 * smi1_1398[k]
                    + f_3 * pc_y[k] * smk_1797[k];

        t_2247[k] = f_19 * slk_1474[k]
                    + f_12 * smi0_1399[k]
                    - f_13 * smi1_1399[k]
                    + f_3 * pc_y[k] * smk_1798[k];

        t_2248[k] = f_19 * slk_1475[k]
                    + f_3 * pc_y[k] * smk_1799[k];
    }

#pragma omp simd aligned(t_2249, t_2250, t_2251, t_2252, pc_x, pc_y, pc_z, slk_1439, slk_1440, \
                         slk_1476, smi0_1399, smi0_1400, smi1_1399, smi1_1400, smk_1799, \
                         smk_1800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2249[k] = f_18 * slk_1439[k]
                    + f_1 * smi0_1399[k]
                    - f_2 * smi1_1399[k]
                    + f_3 * pc_z[k] * smk_1799[k];

        t_2250[k] = f_1 * smi0_1400[k]
                    - f_2 * smi1_1400[k]
                    + f_3 * pc_x[k] * smk_1800[k];

        t_2251[k] = f_18 * slk_1476[k]
                    + f_3 * pc_y[k] * smk_1800[k];

        t_2252[k] = f_19 * slk_1440[k]
                    + f_3 * pc_z[k] * smk_1800[k];
    }

#pragma omp simd aligned(t_2253, t_2254, t_2255, pc_x, pc_y, slk_1478, smi0_1403, smi0_1405, \
                         smi1_1403, smi1_1405, smk_1802, smk_1803, \
                         smk_1805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2253[k] = f_4 * smi0_1403[k]
                    - f_5 * smi1_1403[k]
                    + f_3 * pc_x[k] * smk_1803[k];

        t_2254[k] = f_18 * slk_1478[k]
                    + f_3 * pc_y[k] * smk_1802[k];

        t_2255[k] = f_4 * smi0_1405[k]
                    - f_5 * smi1_1405[k]
                    + f_3 * pc_x[k] * smk_1805[k];
    }

#pragma omp simd aligned(t_2256, t_2257, t_2258, pc_x, pc_y, pc_z, slk_1443, slk_1481, \
                         smi0_1406, smi1_1406, smk_1803, smk_1805, \
                         smk_1806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2256[k] = f_6 * smi0_1406[k]
                    - f_7 * smi1_1406[k]
                    + f_3 * pc_x[k] * smk_1806[k];

        t_2257[k] = f_19 * slk_1443[k]
                    + f_3 * pc_z[k] * smk_1803[k];

        t_2258[k] = f_18 * slk_1481[k]
                    + f_3 * pc_y[k] * smk_1805[k];
    }

#pragma omp simd aligned(t_2259, t_2260, t_2261, pc_x, pc_z, slk_1446, smi0_1409, smi0_1410, \
                         smi1_1409, smi1_1410, smk_1806, smk_1809, \
                         smk_1810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2259[k] = f_6 * smi0_1409[k]
                    - f_7 * smi1_1409[k]
                    + f_3 * pc_x[k] * smk_1809[k];

        t_2260[k] = f_8 * smi0_1410[k]
                    - f_9 * smi1_1410[k]
                    + f_3 * pc_x[k] * smk_1810[k];

        t_2261[k] = f_19 * slk_1446[k]
                    + f_3 * pc_z[k] * smk_1806[k];
    }

#pragma omp simd aligned(t_2262, t_2263, t_2264, pc_x, pc_y, slk_1485, smi0_1412, smi0_1414, \
                         smi1_1412, smi1_1414, smk_1809, smk_1812, \
                         smk_1814 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2262[k] = f_8 * smi0_1412[k]
                    - f_9 * smi1_1412[k]
                    + f_3 * pc_x[k] * smk_1812[k];

        t_2263[k] = f_18 * slk_1485[k]
                    + f_3 * pc_y[k] * smk_1809[k];

        t_2264[k] = f_8 * smi0_1414[k]
                    - f_9 * smi1_1414[k]
                    + f_3 * pc_x[k] * smk_1814[k];
    }

#pragma omp simd aligned(t_2265, t_2266, t_2267, pc_x, pc_z, slk_1450, smi0_1415, smi0_1417, \
                         smi1_1415, smi1_1417, smk_1810, smk_1815, \
                         smk_1817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2265[k] = f_10 * smi0_1415[k]
                    - f_11 * smi1_1415[k]
                    + f_3 * pc_x[k] * smk_1815[k];

        t_2266[k] = f_19 * slk_1450[k]
                    + f_3 * pc_z[k] * smk_1810[k];

        t_2267[k] = f_10 * smi0_1417[k]
                    - f_11 * smi1_1417[k]
                    + f_3 * pc_x[k] * smk_1817[k];
    }

#pragma omp simd aligned(t_2268, t_2269, t_2270, pc_x, pc_y, slk_1490, smi0_1418, smi0_1420, \
                         smi1_1418, smi1_1420, smk_1814, smk_1818, \
                         smk_1820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2268[k] = f_10 * smi0_1418[k]
                    - f_11 * smi1_1418[k]
                    + f_3 * pc_x[k] * smk_1818[k];

        t_2269[k] = f_18 * slk_1490[k]
                    + f_3 * pc_y[k] * smk_1814[k];

        t_2270[k] = f_10 * smi0_1420[k]
                    - f_11 * smi1_1420[k]
                    + f_3 * pc_x[k] * smk_1820[k];
    }

#pragma omp simd aligned(t_2271, t_2272, t_2273, pc_x, pc_z, slk_1455, smi0_1421, smi0_1423, \
                         smi1_1421, smi1_1423, smk_1815, smk_1821, \
                         smk_1823 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2271[k] = f_12 * smi0_1421[k]
                    - f_13 * smi1_1421[k]
                    + f_3 * pc_x[k] * smk_1821[k];

        t_2272[k] = f_19 * slk_1455[k]
                    + f_3 * pc_z[k] * smk_1815[k];

        t_2273[k] = f_12 * smi0_1423[k]
                    - f_13 * smi1_1423[k]
                    + f_3 * pc_x[k] * smk_1823[k];
    }

#pragma omp simd aligned(t_2274, t_2275, t_2276, pc_x, pc_y, slk_1496, smi0_1424, smi0_1425, \
                         smi1_1424, smi1_1425, smk_1820, smk_1824, \
                         smk_1825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2274[k] = f_12 * smi0_1424[k]
                    - f_13 * smi1_1424[k]
                    + f_3 * pc_x[k] * smk_1824[k];

        t_2275[k] = f_12 * smi0_1425[k]
                    - f_13 * smi1_1425[k]
                    + f_3 * pc_x[k] * smk_1825[k];

        t_2276[k] = f_18 * slk_1496[k]
                    + f_3 * pc_y[k] * smk_1820[k];
    }

#pragma omp simd aligned(t_2277, t_2278, t_2279, t_2280, t_2281, t_2282, pc_x, smi0_1427, \
                         smi1_1427, smk_1827, smk_1828, smk_1829, smk_1830, smk_1831, \
                         smk_1832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2277[k] = f_12 * smi0_1427[k]
                    - f_13 * smi1_1427[k]
                    + f_3 * pc_x[k] * smk_1827[k];

        t_2278[k] = f_3 * pc_x[k] * smk_1828[k];

        t_2279[k] = f_3 * pc_x[k] * smk_1829[k];

        t_2280[k] = f_3 * pc_x[k] * smk_1830[k];

        t_2281[k] = f_3 * pc_x[k] * smk_1831[k];

        t_2282[k] = f_3 * pc_x[k] * smk_1832[k];
    }

#pragma omp simd aligned(t_2283, t_2284, t_2285, t_2286, t_2287, pc_x, pc_y, pc_z, slk_1468, \
                         slk_1504, smi0_1421, smi1_1421, smk_1828, smk_1833, smk_1834, \
                         smk_1835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2283[k] = f_3 * pc_x[k] * smk_1833[k];

        t_2284[k] = f_3 * pc_x[k] * smk_1834[k];

        t_2285[k] = f_3 * pc_x[k] * smk_1835[k];

        t_2286[k] = f_18 * slk_1504[k]
                    + f_1 * smi0_1421[k]
                    - f_2 * smi1_1421[k]
                    + f_3 * pc_y[k] * smk_1828[k];

        t_2287[k] = f_19 * slk_1468[k]
                    + f_3 * pc_z[k] * smk_1828[k];
    }

#pragma omp simd aligned(t_2288, t_2289, t_2290, pc_y, slk_1506, slk_1507, slk_1508, \
                         smi0_1423, smi0_1424, smi0_1425, smi1_1423, smi1_1424, smi1_1425, \
                         smk_1830, smk_1831, smk_1832 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2288[k] = f_18 * slk_1506[k]
                    + f_4 * smi0_1423[k]
                    - f_5 * smi1_1423[k]
                    + f_3 * pc_y[k] * smk_1830[k];

        t_2289[k] = f_18 * slk_1507[k]
                    + f_6 * smi0_1424[k]
                    - f_7 * smi1_1424[k]
                    + f_3 * pc_y[k] * smk_1831[k];

        t_2290[k] = f_18 * slk_1508[k]
                    + f_8 * smi0_1425[k]
                    - f_9 * smi1_1425[k]
                    + f_3 * pc_y[k] * smk_1832[k];
    }

#pragma omp simd aligned(t_2291, t_2292, t_2293, pc_y, slk_1509, slk_1510, slk_1511, \
                         smi0_1426, smi0_1427, smi1_1426, smi1_1427, smk_1833, smk_1834, \
                         smk_1835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2291[k] = f_18 * slk_1509[k]
                    + f_10 * smi0_1426[k]
                    - f_11 * smi1_1426[k]
                    + f_3 * pc_y[k] * smk_1833[k];

        t_2292[k] = f_18 * slk_1510[k]
                    + f_12 * smi0_1427[k]
                    - f_13 * smi1_1427[k]
                    + f_3 * pc_y[k] * smk_1834[k];

        t_2293[k] = f_18 * slk_1511[k]
                    + f_3 * pc_y[k] * smk_1835[k];
    }

#pragma omp simd aligned(t_2294, t_2295, t_2296, t_2297, pc_x, pc_y, pc_z, slk_1475, slk_1476, \
                         slk_1512, smi0_1427, smi0_1428, smi1_1427, smi1_1428, smk_1835, \
                         smk_1836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2294[k] = f_19 * slk_1475[k]
                    + f_1 * smi0_1427[k]
                    - f_2 * smi1_1427[k]
                    + f_3 * pc_z[k] * smk_1835[k];

        t_2295[k] = f_1 * smi0_1428[k]
                    - f_2 * smi1_1428[k]
                    + f_3 * pc_x[k] * smk_1836[k];

        t_2296[k] = f_17 * slk_1512[k]
                    + f_3 * pc_y[k] * smk_1836[k];

        t_2297[k] = f_20 * slk_1476[k]
                    + f_3 * pc_z[k] * smk_1836[k];
    }

#pragma omp simd aligned(t_2298, t_2299, t_2300, pc_x, pc_y, slk_1514, smi0_1431, smi0_1433, \
                         smi1_1431, smi1_1433, smk_1838, smk_1839, \
                         smk_1841 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2298[k] = f_4 * smi0_1431[k]
                    - f_5 * smi1_1431[k]
                    + f_3 * pc_x[k] * smk_1839[k];

        t_2299[k] = f_17 * slk_1514[k]
                    + f_3 * pc_y[k] * smk_1838[k];

        t_2300[k] = f_4 * smi0_1433[k]
                    - f_5 * smi1_1433[k]
                    + f_3 * pc_x[k] * smk_1841[k];
    }

#pragma omp simd aligned(t_2301, t_2302, t_2303, pc_x, pc_y, pc_z, slk_1479, slk_1517, \
                         smi0_1434, smi1_1434, smk_1839, smk_1841, \
                         smk_1842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2301[k] = f_6 * smi0_1434[k]
                    - f_7 * smi1_1434[k]
                    + f_3 * pc_x[k] * smk_1842[k];

        t_2302[k] = f_20 * slk_1479[k]
                    + f_3 * pc_z[k] * smk_1839[k];

        t_2303[k] = f_17 * slk_1517[k]
                    + f_3 * pc_y[k] * smk_1841[k];
    }

#pragma omp simd aligned(t_2304, t_2305, t_2306, pc_x, pc_z, slk_1482, smi0_1437, smi0_1438, \
                         smi1_1437, smi1_1438, smk_1842, smk_1845, \
                         smk_1846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2304[k] = f_6 * smi0_1437[k]
                    - f_7 * smi1_1437[k]
                    + f_3 * pc_x[k] * smk_1845[k];

        t_2305[k] = f_8 * smi0_1438[k]
                    - f_9 * smi1_1438[k]
                    + f_3 * pc_x[k] * smk_1846[k];

        t_2306[k] = f_20 * slk_1482[k]
                    + f_3 * pc_z[k] * smk_1842[k];
    }

#pragma omp simd aligned(t_2307, t_2308, t_2309, pc_x, pc_y, slk_1521, smi0_1440, smi0_1442, \
                         smi1_1440, smi1_1442, smk_1845, smk_1848, \
                         smk_1850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2307[k] = f_8 * smi0_1440[k]
                    - f_9 * smi1_1440[k]
                    + f_3 * pc_x[k] * smk_1848[k];

        t_2308[k] = f_17 * slk_1521[k]
                    + f_3 * pc_y[k] * smk_1845[k];

        t_2309[k] = f_8 * smi0_1442[k]
                    - f_9 * smi1_1442[k]
                    + f_3 * pc_x[k] * smk_1850[k];
    }

#pragma omp simd aligned(t_2310, t_2311, t_2312, pc_x, pc_z, slk_1486, smi0_1443, smi0_1445, \
                         smi1_1443, smi1_1445, smk_1846, smk_1851, \
                         smk_1853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2310[k] = f_10 * smi0_1443[k]
                    - f_11 * smi1_1443[k]
                    + f_3 * pc_x[k] * smk_1851[k];

        t_2311[k] = f_20 * slk_1486[k]
                    + f_3 * pc_z[k] * smk_1846[k];

        t_2312[k] = f_10 * smi0_1445[k]
                    - f_11 * smi1_1445[k]
                    + f_3 * pc_x[k] * smk_1853[k];
    }

#pragma omp simd aligned(t_2313, t_2314, t_2315, pc_x, pc_y, slk_1526, smi0_1446, smi0_1448, \
                         smi1_1446, smi1_1448, smk_1850, smk_1854, \
                         smk_1856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2313[k] = f_10 * smi0_1446[k]
                    - f_11 * smi1_1446[k]
                    + f_3 * pc_x[k] * smk_1854[k];

        t_2314[k] = f_17 * slk_1526[k]
                    + f_3 * pc_y[k] * smk_1850[k];

        t_2315[k] = f_10 * smi0_1448[k]
                    - f_11 * smi1_1448[k]
                    + f_3 * pc_x[k] * smk_1856[k];
    }

#pragma omp simd aligned(t_2316, t_2317, t_2318, pc_x, pc_z, slk_1491, smi0_1449, smi0_1451, \
                         smi1_1449, smi1_1451, smk_1851, smk_1857, \
                         smk_1859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2316[k] = f_12 * smi0_1449[k]
                    - f_13 * smi1_1449[k]
                    + f_3 * pc_x[k] * smk_1857[k];

        t_2317[k] = f_20 * slk_1491[k]
                    + f_3 * pc_z[k] * smk_1851[k];

        t_2318[k] = f_12 * smi0_1451[k]
                    - f_13 * smi1_1451[k]
                    + f_3 * pc_x[k] * smk_1859[k];
    }

#pragma omp simd aligned(t_2319, t_2320, t_2321, pc_x, pc_y, slk_1532, smi0_1452, smi0_1453, \
                         smi1_1452, smi1_1453, smk_1856, smk_1860, \
                         smk_1861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2319[k] = f_12 * smi0_1452[k]
                    - f_13 * smi1_1452[k]
                    + f_3 * pc_x[k] * smk_1860[k];

        t_2320[k] = f_12 * smi0_1453[k]
                    - f_13 * smi1_1453[k]
                    + f_3 * pc_x[k] * smk_1861[k];

        t_2321[k] = f_17 * slk_1532[k]
                    + f_3 * pc_y[k] * smk_1856[k];
    }

#pragma omp simd aligned(t_2322, t_2323, t_2324, t_2325, t_2326, t_2327, pc_x, smi0_1455, \
                         smi1_1455, smk_1863, smk_1864, smk_1865, smk_1866, smk_1867, \
                         smk_1868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2322[k] = f_12 * smi0_1455[k]
                    - f_13 * smi1_1455[k]
                    + f_3 * pc_x[k] * smk_1863[k];

        t_2323[k] = f_3 * pc_x[k] * smk_1864[k];

        t_2324[k] = f_3 * pc_x[k] * smk_1865[k];

        t_2325[k] = f_3 * pc_x[k] * smk_1866[k];

        t_2326[k] = f_3 * pc_x[k] * smk_1867[k];

        t_2327[k] = f_3 * pc_x[k] * smk_1868[k];
    }

#pragma omp simd aligned(t_2328, t_2329, t_2330, t_2331, t_2332, pc_x, pc_y, pc_z, slk_1504, \
                         slk_1540, smi0_1449, smi1_1449, smk_1864, smk_1869, smk_1870, \
                         smk_1871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2328[k] = f_3 * pc_x[k] * smk_1869[k];

        t_2329[k] = f_3 * pc_x[k] * smk_1870[k];

        t_2330[k] = f_3 * pc_x[k] * smk_1871[k];

        t_2331[k] = f_17 * slk_1540[k]
                    + f_1 * smi0_1449[k]
                    - f_2 * smi1_1449[k]
                    + f_3 * pc_y[k] * smk_1864[k];

        t_2332[k] = f_20 * slk_1504[k]
                    + f_3 * pc_z[k] * smk_1864[k];
    }

#pragma omp simd aligned(t_2333, t_2334, t_2335, pc_y, slk_1542, slk_1543, slk_1544, \
                         smi0_1451, smi0_1452, smi0_1453, smi1_1451, smi1_1452, smi1_1453, \
                         smk_1866, smk_1867, smk_1868 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2333[k] = f_17 * slk_1542[k]
                    + f_4 * smi0_1451[k]
                    - f_5 * smi1_1451[k]
                    + f_3 * pc_y[k] * smk_1866[k];

        t_2334[k] = f_17 * slk_1543[k]
                    + f_6 * smi0_1452[k]
                    - f_7 * smi1_1452[k]
                    + f_3 * pc_y[k] * smk_1867[k];

        t_2335[k] = f_17 * slk_1544[k]
                    + f_8 * smi0_1453[k]
                    - f_9 * smi1_1453[k]
                    + f_3 * pc_y[k] * smk_1868[k];
    }

#pragma omp simd aligned(t_2336, t_2337, t_2338, pc_y, slk_1545, slk_1546, slk_1547, \
                         smi0_1454, smi0_1455, smi1_1454, smi1_1455, smk_1869, smk_1870, \
                         smk_1871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2336[k] = f_17 * slk_1545[k]
                    + f_10 * smi0_1454[k]
                    - f_11 * smi1_1454[k]
                    + f_3 * pc_y[k] * smk_1869[k];

        t_2337[k] = f_17 * slk_1546[k]
                    + f_12 * smi0_1455[k]
                    - f_13 * smi1_1455[k]
                    + f_3 * pc_y[k] * smk_1870[k];

        t_2338[k] = f_17 * slk_1547[k]
                    + f_3 * pc_y[k] * smk_1871[k];
    }

#pragma omp simd aligned(t_2339, t_2340, t_2341, t_2342, pc_x, pc_y, pc_z, slk_1511, slk_1512, \
                         slk_1548, smi0_1455, smi0_1456, smi1_1455, smi1_1456, smk_1871, \
                         smk_1872 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2339[k] = f_20 * slk_1511[k]
                    + f_1 * smi0_1455[k]
                    - f_2 * smi1_1455[k]
                    + f_3 * pc_z[k] * smk_1871[k];

        t_2340[k] = f_1 * smi0_1456[k]
                    - f_2 * smi1_1456[k]
                    + f_3 * pc_x[k] * smk_1872[k];

        t_2341[k] = f_16 * slk_1548[k]
                    + f_3 * pc_y[k] * smk_1872[k];

        t_2342[k] = f_22 * slk_1512[k]
                    + f_3 * pc_z[k] * smk_1872[k];
    }

#pragma omp simd aligned(t_2343, t_2344, t_2345, pc_x, pc_y, slk_1550, smi0_1459, smi0_1461, \
                         smi1_1459, smi1_1461, smk_1874, smk_1875, \
                         smk_1877 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2343[k] = f_4 * smi0_1459[k]
                    - f_5 * smi1_1459[k]
                    + f_3 * pc_x[k] * smk_1875[k];

        t_2344[k] = f_16 * slk_1550[k]
                    + f_3 * pc_y[k] * smk_1874[k];

        t_2345[k] = f_4 * smi0_1461[k]
                    - f_5 * smi1_1461[k]
                    + f_3 * pc_x[k] * smk_1877[k];
    }

#pragma omp simd aligned(t_2346, t_2347, t_2348, pc_x, pc_y, pc_z, slk_1515, slk_1553, \
                         smi0_1462, smi1_1462, smk_1875, smk_1877, \
                         smk_1878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2346[k] = f_6 * smi0_1462[k]
                    - f_7 * smi1_1462[k]
                    + f_3 * pc_x[k] * smk_1878[k];

        t_2347[k] = f_22 * slk_1515[k]
                    + f_3 * pc_z[k] * smk_1875[k];

        t_2348[k] = f_16 * slk_1553[k]
                    + f_3 * pc_y[k] * smk_1877[k];
    }

#pragma omp simd aligned(t_2349, t_2350, t_2351, pc_x, pc_z, slk_1518, smi0_1465, smi0_1466, \
                         smi1_1465, smi1_1466, smk_1878, smk_1881, \
                         smk_1882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2349[k] = f_6 * smi0_1465[k]
                    - f_7 * smi1_1465[k]
                    + f_3 * pc_x[k] * smk_1881[k];

        t_2350[k] = f_8 * smi0_1466[k]
                    - f_9 * smi1_1466[k]
                    + f_3 * pc_x[k] * smk_1882[k];

        t_2351[k] = f_22 * slk_1518[k]
                    + f_3 * pc_z[k] * smk_1878[k];
    }

#pragma omp simd aligned(t_2352, t_2353, t_2354, pc_x, pc_y, slk_1557, smi0_1468, smi0_1470, \
                         smi1_1468, smi1_1470, smk_1881, smk_1884, \
                         smk_1886 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2352[k] = f_8 * smi0_1468[k]
                    - f_9 * smi1_1468[k]
                    + f_3 * pc_x[k] * smk_1884[k];

        t_2353[k] = f_16 * slk_1557[k]
                    + f_3 * pc_y[k] * smk_1881[k];

        t_2354[k] = f_8 * smi0_1470[k]
                    - f_9 * smi1_1470[k]
                    + f_3 * pc_x[k] * smk_1886[k];
    }
}

static auto
compute_prim_sml_three_center_electron_repulsion_0_piece21(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sll0,
                                                           const size_t slk, const size_t sll1,
                                                           const size_t smi0, const size_t smi1,
                                                           const size_t smk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;
    const auto f_22 = 3.5 / q;

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
    auto *t_2376 = buffer.data(target + 2376);
    auto *t_2377 = buffer.data(target + 2377);
    auto *t_2378 = buffer.data(target + 2378);
    auto *t_2379 = buffer.data(target + 2379);
    auto *t_2380 = buffer.data(target + 2380);
    auto *t_2381 = buffer.data(target + 2381);
    auto *t_2382 = buffer.data(target + 2382);
    auto *t_2383 = buffer.data(target + 2383);
    auto *t_2384 = buffer.data(target + 2384);
    auto *t_2385 = buffer.data(target + 2385);
    auto *t_2386 = buffer.data(target + 2386);
    auto *t_2387 = buffer.data(target + 2387);
    auto *t_2388 = buffer.data(target + 2388);
    auto *t_2389 = buffer.data(target + 2389);
    auto *t_2390 = buffer.data(target + 2390);
    auto *t_2391 = buffer.data(target + 2391);
    auto *t_2392 = buffer.data(target + 2392);
    auto *t_2393 = buffer.data(target + 2393);
    auto *t_2394 = buffer.data(target + 2394);
    auto *t_2395 = buffer.data(target + 2395);
    auto *t_2396 = buffer.data(target + 2396);
    auto *t_2397 = buffer.data(target + 2397);
    auto *t_2398 = buffer.data(target + 2398);
    auto *t_2399 = buffer.data(target + 2399);
    auto *t_2400 = buffer.data(target + 2400);
    auto *t_2401 = buffer.data(target + 2401);
    auto *t_2402 = buffer.data(target + 2402);
    auto *t_2403 = buffer.data(target + 2403);
    auto *t_2404 = buffer.data(target + 2404);
    auto *t_2405 = buffer.data(target + 2405);
    auto *t_2406 = buffer.data(target + 2406);
    auto *t_2407 = buffer.data(target + 2407);
    auto *t_2408 = buffer.data(target + 2408);
    auto *t_2409 = buffer.data(target + 2409);
    auto *t_2410 = buffer.data(target + 2410);
    auto *t_2411 = buffer.data(target + 2411);
    auto *t_2412 = buffer.data(target + 2412);
    auto *t_2413 = buffer.data(target + 2413);
    auto *t_2414 = buffer.data(target + 2414);
    auto *t_2415 = buffer.data(target + 2415);
    auto *t_2416 = buffer.data(target + 2416);
    auto *t_2417 = buffer.data(target + 2417);
    auto *t_2418 = buffer.data(target + 2418);
    auto *t_2419 = buffer.data(target + 2419);
    auto *t_2420 = buffer.data(target + 2420);
    auto *t_2421 = buffer.data(target + 2421);
    auto *t_2422 = buffer.data(target + 2422);
    auto *t_2423 = buffer.data(target + 2423);
    auto *t_2424 = buffer.data(target + 2424);
    auto *t_2425 = buffer.data(target + 2425);
    auto *t_2426 = buffer.data(target + 2426);
    auto *t_2427 = buffer.data(target + 2427);
    auto *t_2428 = buffer.data(target + 2428);
    auto *t_2429 = buffer.data(target + 2429);
    auto *t_2430 = buffer.data(target + 2430);
    auto *t_2431 = buffer.data(target + 2431);
    auto *t_2432 = buffer.data(target + 2432);
    auto *t_2433 = buffer.data(target + 2433);
    auto *t_2434 = buffer.data(target + 2434);
    auto *t_2435 = buffer.data(target + 2435);
    auto *t_2436 = buffer.data(target + 2436);
    auto *t_2437 = buffer.data(target + 2437);
    auto *t_2438 = buffer.data(target + 2438);
    auto *t_2439 = buffer.data(target + 2439);
    auto *t_2440 = buffer.data(target + 2440);
    auto *t_2441 = buffer.data(target + 2441);
    auto *t_2442 = buffer.data(target + 2442);
    auto *t_2443 = buffer.data(target + 2443);
    auto *t_2444 = buffer.data(target + 2444);
    auto *t_2445 = buffer.data(target + 2445);
    auto *t_2446 = buffer.data(target + 2446);
    auto *t_2447 = buffer.data(target + 2447);
    auto *t_2448 = buffer.data(target + 2448);
    auto *t_2449 = buffer.data(target + 2449);
    auto *t_2450 = buffer.data(target + 2450);
    auto *t_2451 = buffer.data(target + 2451);
    auto *t_2452 = buffer.data(target + 2452);
    auto *t_2453 = buffer.data(target + 2453);
    auto *t_2454 = buffer.data(target + 2454);
    auto *t_2455 = buffer.data(target + 2455);
    auto *t_2456 = buffer.data(target + 2456);
    auto *t_2457 = buffer.data(target + 2457);
    auto *t_2458 = buffer.data(target + 2458);
    auto *t_2459 = buffer.data(target + 2459);
    auto *t_2460 = buffer.data(target + 2460);
    auto *t_2461 = buffer.data(target + 2461);
    auto *t_2462 = buffer.data(target + 2462);
    auto *t_2463 = buffer.data(target + 2463);
    auto *t_2464 = buffer.data(target + 2464);
    auto *t_2465 = buffer.data(target + 2465);
    auto *t_2466 = buffer.data(target + 2466);
    auto *t_2467 = buffer.data(target + 2467);
    auto *t_2468 = buffer.data(target + 2468);
    auto *t_2469 = buffer.data(target + 2469);
    auto *t_2470 = buffer.data(target + 2470);
    auto *t_2471 = buffer.data(target + 2471);
    auto *t_2472 = buffer.data(target + 2472);
    auto *t_2473 = buffer.data(target + 2473);
    auto *t_2474 = buffer.data(target + 2474);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sll0_1980 = buffer.data(sll0 + 1980);
    const auto *sll0_1985 = buffer.data(sll0 + 1985);
    const auto *sll0_1989 = buffer.data(sll0 + 1989);
    const auto *sll0_1994 = buffer.data(sll0 + 1994);
    const auto *sll0_2000 = buffer.data(sll0 + 2000);
    const auto *sll0_2007 = buffer.data(sll0 + 2007);
    const auto *sll0_2016 = buffer.data(sll0 + 2016);
    const auto *sll0_2018 = buffer.data(sll0 + 2018);
    const auto *sll0_2019 = buffer.data(sll0 + 2019);
    const auto *sll0_2020 = buffer.data(sll0 + 2020);
    const auto *sll0_2021 = buffer.data(sll0 + 2021);
    const auto *sll0_2022 = buffer.data(sll0 + 2022);
    const auto *sll0_2024 = buffer.data(sll0 + 2024);

    const auto *slk_1522 = buffer.data(slk + 1522);
    const auto *slk_1527 = buffer.data(slk + 1527);
    const auto *slk_1540 = buffer.data(slk + 1540);
    const auto *slk_1547 = buffer.data(slk + 1547);
    const auto *slk_1548 = buffer.data(slk + 1548);
    const auto *slk_1551 = buffer.data(slk + 1551);
    const auto *slk_1554 = buffer.data(slk + 1554);
    const auto *slk_1558 = buffer.data(slk + 1558);
    const auto *slk_1562 = buffer.data(slk + 1562);
    const auto *slk_1563 = buffer.data(slk + 1563);
    const auto *slk_1568 = buffer.data(slk + 1568);
    const auto *slk_1576 = buffer.data(slk + 1576);
    const auto *slk_1578 = buffer.data(slk + 1578);
    const auto *slk_1579 = buffer.data(slk + 1579);
    const auto *slk_1580 = buffer.data(slk + 1580);
    const auto *slk_1581 = buffer.data(slk + 1581);
    const auto *slk_1582 = buffer.data(slk + 1582);
    const auto *slk_1583 = buffer.data(slk + 1583);
    const auto *slk_1584 = buffer.data(slk + 1584);
    const auto *slk_1586 = buffer.data(slk + 1586);
    const auto *slk_1587 = buffer.data(slk + 1587);
    const auto *slk_1589 = buffer.data(slk + 1589);
    const auto *slk_1590 = buffer.data(slk + 1590);
    const auto *slk_1593 = buffer.data(slk + 1593);
    const auto *slk_1594 = buffer.data(slk + 1594);
    const auto *slk_1598 = buffer.data(slk + 1598);
    const auto *slk_1599 = buffer.data(slk + 1599);
    const auto *slk_1604 = buffer.data(slk + 1604);
    const auto *slk_1612 = buffer.data(slk + 1612);
    const auto *slk_1614 = buffer.data(slk + 1614);
    const auto *slk_1615 = buffer.data(slk + 1615);
    const auto *slk_1616 = buffer.data(slk + 1616);
    const auto *slk_1617 = buffer.data(slk + 1617);
    const auto *slk_1618 = buffer.data(slk + 1618);
    const auto *slk_1619 = buffer.data(slk + 1619);

    const auto *sll1_1980 = buffer.data(sll1 + 1980);
    const auto *sll1_1985 = buffer.data(sll1 + 1985);
    const auto *sll1_1989 = buffer.data(sll1 + 1989);
    const auto *sll1_1994 = buffer.data(sll1 + 1994);
    const auto *sll1_2000 = buffer.data(sll1 + 2000);
    const auto *sll1_2007 = buffer.data(sll1 + 2007);
    const auto *sll1_2016 = buffer.data(sll1 + 2016);
    const auto *sll1_2018 = buffer.data(sll1 + 2018);
    const auto *sll1_2019 = buffer.data(sll1 + 2019);
    const auto *sll1_2020 = buffer.data(sll1 + 2020);
    const auto *sll1_2021 = buffer.data(sll1 + 2021);
    const auto *sll1_2022 = buffer.data(sll1 + 2022);
    const auto *sll1_2024 = buffer.data(sll1 + 2024);

    const auto *smi0_1471 = buffer.data(smi0 + 1471);
    const auto *smi0_1473 = buffer.data(smi0 + 1473);
    const auto *smi0_1474 = buffer.data(smi0 + 1474);
    const auto *smi0_1476 = buffer.data(smi0 + 1476);
    const auto *smi0_1477 = buffer.data(smi0 + 1477);
    const auto *smi0_1479 = buffer.data(smi0 + 1479);
    const auto *smi0_1480 = buffer.data(smi0 + 1480);
    const auto *smi0_1481 = buffer.data(smi0 + 1481);
    const auto *smi0_1482 = buffer.data(smi0 + 1482);
    const auto *smi0_1483 = buffer.data(smi0 + 1483);
    const auto *smi0_1487 = buffer.data(smi0 + 1487);
    const auto *smi0_1490 = buffer.data(smi0 + 1490);
    const auto *smi0_1494 = buffer.data(smi0 + 1494);
    const auto *smi0_1496 = buffer.data(smi0 + 1496);
    const auto *smi0_1499 = buffer.data(smi0 + 1499);
    const auto *smi0_1501 = buffer.data(smi0 + 1501);
    const auto *smi0_1502 = buffer.data(smi0 + 1502);
    const auto *smi0_1505 = buffer.data(smi0 + 1505);
    const auto *smi0_1507 = buffer.data(smi0 + 1507);
    const auto *smi0_1508 = buffer.data(smi0 + 1508);
    const auto *smi0_1509 = buffer.data(smi0 + 1509);
    const auto *smi0_1512 = buffer.data(smi0 + 1512);
    const auto *smi0_1515 = buffer.data(smi0 + 1515);
    const auto *smi0_1517 = buffer.data(smi0 + 1517);
    const auto *smi0_1518 = buffer.data(smi0 + 1518);
    const auto *smi0_1521 = buffer.data(smi0 + 1521);
    const auto *smi0_1522 = buffer.data(smi0 + 1522);
    const auto *smi0_1524 = buffer.data(smi0 + 1524);
    const auto *smi0_1526 = buffer.data(smi0 + 1526);
    const auto *smi0_1527 = buffer.data(smi0 + 1527);
    const auto *smi0_1529 = buffer.data(smi0 + 1529);
    const auto *smi0_1530 = buffer.data(smi0 + 1530);
    const auto *smi0_1532 = buffer.data(smi0 + 1532);
    const auto *smi0_1533 = buffer.data(smi0 + 1533);
    const auto *smi0_1535 = buffer.data(smi0 + 1535);
    const auto *smi0_1536 = buffer.data(smi0 + 1536);
    const auto *smi0_1537 = buffer.data(smi0 + 1537);
    const auto *smi0_1538 = buffer.data(smi0 + 1538);
    const auto *smi0_1539 = buffer.data(smi0 + 1539);

    const auto *smi1_1471 = buffer.data(smi1 + 1471);
    const auto *smi1_1473 = buffer.data(smi1 + 1473);
    const auto *smi1_1474 = buffer.data(smi1 + 1474);
    const auto *smi1_1476 = buffer.data(smi1 + 1476);
    const auto *smi1_1477 = buffer.data(smi1 + 1477);
    const auto *smi1_1479 = buffer.data(smi1 + 1479);
    const auto *smi1_1480 = buffer.data(smi1 + 1480);
    const auto *smi1_1481 = buffer.data(smi1 + 1481);
    const auto *smi1_1482 = buffer.data(smi1 + 1482);
    const auto *smi1_1483 = buffer.data(smi1 + 1483);
    const auto *smi1_1487 = buffer.data(smi1 + 1487);
    const auto *smi1_1490 = buffer.data(smi1 + 1490);
    const auto *smi1_1494 = buffer.data(smi1 + 1494);
    const auto *smi1_1496 = buffer.data(smi1 + 1496);
    const auto *smi1_1499 = buffer.data(smi1 + 1499);
    const auto *smi1_1501 = buffer.data(smi1 + 1501);
    const auto *smi1_1502 = buffer.data(smi1 + 1502);
    const auto *smi1_1505 = buffer.data(smi1 + 1505);
    const auto *smi1_1507 = buffer.data(smi1 + 1507);
    const auto *smi1_1508 = buffer.data(smi1 + 1508);
    const auto *smi1_1509 = buffer.data(smi1 + 1509);
    const auto *smi1_1512 = buffer.data(smi1 + 1512);
    const auto *smi1_1515 = buffer.data(smi1 + 1515);
    const auto *smi1_1517 = buffer.data(smi1 + 1517);
    const auto *smi1_1518 = buffer.data(smi1 + 1518);
    const auto *smi1_1521 = buffer.data(smi1 + 1521);
    const auto *smi1_1522 = buffer.data(smi1 + 1522);
    const auto *smi1_1524 = buffer.data(smi1 + 1524);
    const auto *smi1_1526 = buffer.data(smi1 + 1526);
    const auto *smi1_1527 = buffer.data(smi1 + 1527);
    const auto *smi1_1529 = buffer.data(smi1 + 1529);
    const auto *smi1_1530 = buffer.data(smi1 + 1530);
    const auto *smi1_1532 = buffer.data(smi1 + 1532);
    const auto *smi1_1533 = buffer.data(smi1 + 1533);
    const auto *smi1_1535 = buffer.data(smi1 + 1535);
    const auto *smi1_1536 = buffer.data(smi1 + 1536);
    const auto *smi1_1537 = buffer.data(smi1 + 1537);
    const auto *smi1_1538 = buffer.data(smi1 + 1538);
    const auto *smi1_1539 = buffer.data(smi1 + 1539);

    const auto *smk_1882 = buffer.data(smk + 1882);
    const auto *smk_1886 = buffer.data(smk + 1886);
    const auto *smk_1887 = buffer.data(smk + 1887);
    const auto *smk_1889 = buffer.data(smk + 1889);
    const auto *smk_1890 = buffer.data(smk + 1890);
    const auto *smk_1892 = buffer.data(smk + 1892);
    const auto *smk_1893 = buffer.data(smk + 1893);
    const auto *smk_1895 = buffer.data(smk + 1895);
    const auto *smk_1896 = buffer.data(smk + 1896);
    const auto *smk_1897 = buffer.data(smk + 1897);
    const auto *smk_1899 = buffer.data(smk + 1899);
    const auto *smk_1900 = buffer.data(smk + 1900);
    const auto *smk_1901 = buffer.data(smk + 1901);
    const auto *smk_1902 = buffer.data(smk + 1902);
    const auto *smk_1903 = buffer.data(smk + 1903);
    const auto *smk_1904 = buffer.data(smk + 1904);
    const auto *smk_1905 = buffer.data(smk + 1905);
    const auto *smk_1906 = buffer.data(smk + 1906);
    const auto *smk_1907 = buffer.data(smk + 1907);
    const auto *smk_1908 = buffer.data(smk + 1908);
    const auto *smk_1910 = buffer.data(smk + 1910);
    const auto *smk_1911 = buffer.data(smk + 1911);
    const auto *smk_1913 = buffer.data(smk + 1913);
    const auto *smk_1914 = buffer.data(smk + 1914);
    const auto *smk_1917 = buffer.data(smk + 1917);
    const auto *smk_1918 = buffer.data(smk + 1918);
    const auto *smk_1920 = buffer.data(smk + 1920);
    const auto *smk_1922 = buffer.data(smk + 1922);
    const auto *smk_1923 = buffer.data(smk + 1923);
    const auto *smk_1925 = buffer.data(smk + 1925);
    const auto *smk_1926 = buffer.data(smk + 1926);
    const auto *smk_1928 = buffer.data(smk + 1928);
    const auto *smk_1929 = buffer.data(smk + 1929);
    const auto *smk_1931 = buffer.data(smk + 1931);
    const auto *smk_1932 = buffer.data(smk + 1932);
    const auto *smk_1933 = buffer.data(smk + 1933);
    const auto *smk_1936 = buffer.data(smk + 1936);
    const auto *smk_1937 = buffer.data(smk + 1937);
    const auto *smk_1938 = buffer.data(smk + 1938);
    const auto *smk_1939 = buffer.data(smk + 1939);
    const auto *smk_1940 = buffer.data(smk + 1940);
    const auto *smk_1941 = buffer.data(smk + 1941);
    const auto *smk_1942 = buffer.data(smk + 1942);
    const auto *smk_1943 = buffer.data(smk + 1943);
    const auto *smk_1944 = buffer.data(smk + 1944);
    const auto *smk_1946 = buffer.data(smk + 1946);
    const auto *smk_1947 = buffer.data(smk + 1947);
    const auto *smk_1949 = buffer.data(smk + 1949);
    const auto *smk_1950 = buffer.data(smk + 1950);
    const auto *smk_1953 = buffer.data(smk + 1953);
    const auto *smk_1954 = buffer.data(smk + 1954);
    const auto *smk_1956 = buffer.data(smk + 1956);
    const auto *smk_1958 = buffer.data(smk + 1958);
    const auto *smk_1959 = buffer.data(smk + 1959);
    const auto *smk_1961 = buffer.data(smk + 1961);
    const auto *smk_1962 = buffer.data(smk + 1962);
    const auto *smk_1964 = buffer.data(smk + 1964);
    const auto *smk_1965 = buffer.data(smk + 1965);
    const auto *smk_1967 = buffer.data(smk + 1967);
    const auto *smk_1968 = buffer.data(smk + 1968);
    const auto *smk_1969 = buffer.data(smk + 1969);
    const auto *smk_1971 = buffer.data(smk + 1971);
    const auto *smk_1972 = buffer.data(smk + 1972);
    const auto *smk_1973 = buffer.data(smk + 1973);
    const auto *smk_1974 = buffer.data(smk + 1974);
    const auto *smk_1975 = buffer.data(smk + 1975);
    const auto *smk_1976 = buffer.data(smk + 1976);
    const auto *smk_1977 = buffer.data(smk + 1977);
    const auto *smk_1978 = buffer.data(smk + 1978);
    const auto *smk_1979 = buffer.data(smk + 1979);

#pragma omp simd aligned(t_2355, t_2356, t_2357, pc_x, pc_z, slk_1522, smi0_1471, smi0_1473, \
                         smi1_1471, smi1_1473, smk_1882, smk_1887, \
                         smk_1889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2355[k] = f_10 * smi0_1471[k]
                    - f_11 * smi1_1471[k]
                    + f_3 * pc_x[k] * smk_1887[k];

        t_2356[k] = f_22 * slk_1522[k]
                    + f_3 * pc_z[k] * smk_1882[k];

        t_2357[k] = f_10 * smi0_1473[k]
                    - f_11 * smi1_1473[k]
                    + f_3 * pc_x[k] * smk_1889[k];
    }

#pragma omp simd aligned(t_2358, t_2359, t_2360, pc_x, pc_y, slk_1562, smi0_1474, smi0_1476, \
                         smi1_1474, smi1_1476, smk_1886, smk_1890, \
                         smk_1892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2358[k] = f_10 * smi0_1474[k]
                    - f_11 * smi1_1474[k]
                    + f_3 * pc_x[k] * smk_1890[k];

        t_2359[k] = f_16 * slk_1562[k]
                    + f_3 * pc_y[k] * smk_1886[k];

        t_2360[k] = f_10 * smi0_1476[k]
                    - f_11 * smi1_1476[k]
                    + f_3 * pc_x[k] * smk_1892[k];
    }

#pragma omp simd aligned(t_2361, t_2362, t_2363, pc_x, pc_z, slk_1527, smi0_1477, smi0_1479, \
                         smi1_1477, smi1_1479, smk_1887, smk_1893, \
                         smk_1895 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2361[k] = f_12 * smi0_1477[k]
                    - f_13 * smi1_1477[k]
                    + f_3 * pc_x[k] * smk_1893[k];

        t_2362[k] = f_22 * slk_1527[k]
                    + f_3 * pc_z[k] * smk_1887[k];

        t_2363[k] = f_12 * smi0_1479[k]
                    - f_13 * smi1_1479[k]
                    + f_3 * pc_x[k] * smk_1895[k];
    }

#pragma omp simd aligned(t_2364, t_2365, t_2366, pc_x, pc_y, slk_1568, smi0_1480, smi0_1481, \
                         smi1_1480, smi1_1481, smk_1892, smk_1896, \
                         smk_1897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2364[k] = f_12 * smi0_1480[k]
                    - f_13 * smi1_1480[k]
                    + f_3 * pc_x[k] * smk_1896[k];

        t_2365[k] = f_12 * smi0_1481[k]
                    - f_13 * smi1_1481[k]
                    + f_3 * pc_x[k] * smk_1897[k];

        t_2366[k] = f_16 * slk_1568[k]
                    + f_3 * pc_y[k] * smk_1892[k];
    }

#pragma omp simd aligned(t_2367, t_2368, t_2369, t_2370, t_2371, t_2372, pc_x, smi0_1483, \
                         smi1_1483, smk_1899, smk_1900, smk_1901, smk_1902, smk_1903, \
                         smk_1904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2367[k] = f_12 * smi0_1483[k]
                    - f_13 * smi1_1483[k]
                    + f_3 * pc_x[k] * smk_1899[k];

        t_2368[k] = f_3 * pc_x[k] * smk_1900[k];

        t_2369[k] = f_3 * pc_x[k] * smk_1901[k];

        t_2370[k] = f_3 * pc_x[k] * smk_1902[k];

        t_2371[k] = f_3 * pc_x[k] * smk_1903[k];

        t_2372[k] = f_3 * pc_x[k] * smk_1904[k];
    }

#pragma omp simd aligned(t_2373, t_2374, t_2375, t_2376, t_2377, pc_x, pc_y, pc_z, slk_1540, \
                         slk_1576, smi0_1477, smi1_1477, smk_1900, smk_1905, smk_1906, \
                         smk_1907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2373[k] = f_3 * pc_x[k] * smk_1905[k];

        t_2374[k] = f_3 * pc_x[k] * smk_1906[k];

        t_2375[k] = f_3 * pc_x[k] * smk_1907[k];

        t_2376[k] = f_16 * slk_1576[k]
                    + f_1 * smi0_1477[k]
                    - f_2 * smi1_1477[k]
                    + f_3 * pc_y[k] * smk_1900[k];

        t_2377[k] = f_22 * slk_1540[k]
                    + f_3 * pc_z[k] * smk_1900[k];
    }

#pragma omp simd aligned(t_2378, t_2379, t_2380, pc_y, slk_1578, slk_1579, slk_1580, \
                         smi0_1479, smi0_1480, smi0_1481, smi1_1479, smi1_1480, smi1_1481, \
                         smk_1902, smk_1903, smk_1904 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2378[k] = f_16 * slk_1578[k]
                    + f_4 * smi0_1479[k]
                    - f_5 * smi1_1479[k]
                    + f_3 * pc_y[k] * smk_1902[k];

        t_2379[k] = f_16 * slk_1579[k]
                    + f_6 * smi0_1480[k]
                    - f_7 * smi1_1480[k]
                    + f_3 * pc_y[k] * smk_1903[k];

        t_2380[k] = f_16 * slk_1580[k]
                    + f_8 * smi0_1481[k]
                    - f_9 * smi1_1481[k]
                    + f_3 * pc_y[k] * smk_1904[k];
    }

#pragma omp simd aligned(t_2381, t_2382, t_2383, pc_y, slk_1581, slk_1582, slk_1583, \
                         smi0_1482, smi0_1483, smi1_1482, smi1_1483, smk_1905, smk_1906, \
                         smk_1907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2381[k] = f_16 * slk_1581[k]
                    + f_10 * smi0_1482[k]
                    - f_11 * smi1_1482[k]
                    + f_3 * pc_y[k] * smk_1905[k];

        t_2382[k] = f_16 * slk_1582[k]
                    + f_12 * smi0_1483[k]
                    - f_13 * smi1_1483[k]
                    + f_3 * pc_y[k] * smk_1906[k];

        t_2383[k] = f_16 * slk_1583[k]
                    + f_3 * pc_y[k] * smk_1907[k];
    }

#pragma omp simd aligned(t_2384, t_2385, t_2386, t_2387, pb_y, pc_y, pc_z, sll0_1980, \
                         slk_1547, slk_1548, slk_1584, sll1_1980, smi0_1483, smi1_1483, \
                         smk_1907, smk_1908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2384[k] = f_22 * slk_1547[k]
                    + f_1 * smi0_1483[k]
                    - f_2 * smi1_1483[k]
                    + f_3 * pc_z[k] * smk_1907[k];

        t_2385[k] = pb_y[k] * sll0_1980[k]
                    - f_14 * pc_y[k] * sll1_1980[k];

        t_2386[k] = f_15 * slk_1584[k]
                    + f_3 * pc_y[k] * smk_1908[k];

        t_2387[k] = f_21 * slk_1548[k]
                    + f_3 * pc_z[k] * smk_1908[k];
    }

#pragma omp simd aligned(t_2388, t_2389, t_2390, pb_y, pc_x, pc_y, sll0_1985, slk_1586, \
                         sll1_1985, smi0_1487, smi1_1487, smk_1910, \
                         smk_1911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2388[k] = f_4 * smi0_1487[k]
                    - f_5 * smi1_1487[k]
                    + f_3 * pc_x[k] * smk_1911[k];

        t_2389[k] = f_15 * slk_1586[k]
                    + f_3 * pc_y[k] * smk_1910[k];

        t_2390[k] = pb_y[k] * sll0_1985[k]
                    - f_14 * pc_y[k] * sll1_1985[k];
    }

#pragma omp simd aligned(t_2391, t_2392, t_2393, pc_x, pc_y, pc_z, slk_1551, slk_1589, \
                         smi0_1490, smi1_1490, smk_1911, smk_1913, \
                         smk_1914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2391[k] = f_6 * smi0_1490[k]
                    - f_7 * smi1_1490[k]
                    + f_3 * pc_x[k] * smk_1914[k];

        t_2392[k] = f_21 * slk_1551[k]
                    + f_3 * pc_z[k] * smk_1911[k];

        t_2393[k] = f_15 * slk_1589[k]
                    + f_3 * pc_y[k] * smk_1913[k];
    }

#pragma omp simd aligned(t_2394, t_2395, t_2396, pb_y, pc_x, pc_y, pc_z, sll0_1989, slk_1554, \
                         sll1_1989, smi0_1494, smi1_1494, smk_1914, \
                         smk_1918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2394[k] = pb_y[k] * sll0_1989[k]
                    - f_14 * pc_y[k] * sll1_1989[k];

        t_2395[k] = f_8 * smi0_1494[k]
                    - f_9 * smi1_1494[k]
                    + f_3 * pc_x[k] * smk_1918[k];

        t_2396[k] = f_21 * slk_1554[k]
                    + f_3 * pc_z[k] * smk_1914[k];
    }

#pragma omp simd aligned(t_2397, t_2398, t_2399, pb_y, pc_x, pc_y, sll0_1994, slk_1593, \
                         sll1_1994, smi0_1496, smi1_1496, smk_1917, \
                         smk_1920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2397[k] = f_8 * smi0_1496[k]
                    - f_9 * smi1_1496[k]
                    + f_3 * pc_x[k] * smk_1920[k];

        t_2398[k] = f_15 * slk_1593[k]
                    + f_3 * pc_y[k] * smk_1917[k];

        t_2399[k] = pb_y[k] * sll0_1994[k]
                    - f_14 * pc_y[k] * sll1_1994[k];
    }

#pragma omp simd aligned(t_2400, t_2401, t_2402, pc_x, pc_z, slk_1558, smi0_1499, smi0_1501, \
                         smi1_1499, smi1_1501, smk_1918, smk_1923, \
                         smk_1925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2400[k] = f_10 * smi0_1499[k]
                    - f_11 * smi1_1499[k]
                    + f_3 * pc_x[k] * smk_1923[k];

        t_2401[k] = f_21 * slk_1558[k]
                    + f_3 * pc_z[k] * smk_1918[k];

        t_2402[k] = f_10 * smi0_1501[k]
                    - f_11 * smi1_1501[k]
                    + f_3 * pc_x[k] * smk_1925[k];
    }

#pragma omp simd aligned(t_2403, t_2404, t_2405, pb_y, pc_x, pc_y, sll0_2000, slk_1598, \
                         sll1_2000, smi0_1502, smi1_1502, smk_1922, \
                         smk_1926 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2403[k] = f_10 * smi0_1502[k]
                    - f_11 * smi1_1502[k]
                    + f_3 * pc_x[k] * smk_1926[k];

        t_2404[k] = f_15 * slk_1598[k]
                    + f_3 * pc_y[k] * smk_1922[k];

        t_2405[k] = pb_y[k] * sll0_2000[k]
                    - f_14 * pc_y[k] * sll1_2000[k];
    }

#pragma omp simd aligned(t_2406, t_2407, t_2408, pc_x, pc_z, slk_1563, smi0_1505, smi0_1507, \
                         smi1_1505, smi1_1507, smk_1923, smk_1929, \
                         smk_1931 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2406[k] = f_12 * smi0_1505[k]
                    - f_13 * smi1_1505[k]
                    + f_3 * pc_x[k] * smk_1929[k];

        t_2407[k] = f_21 * slk_1563[k]
                    + f_3 * pc_z[k] * smk_1923[k];

        t_2408[k] = f_12 * smi0_1507[k]
                    - f_13 * smi1_1507[k]
                    + f_3 * pc_x[k] * smk_1931[k];
    }

#pragma omp simd aligned(t_2409, t_2410, t_2411, pc_x, pc_y, slk_1604, smi0_1508, smi0_1509, \
                         smi1_1508, smi1_1509, smk_1928, smk_1932, \
                         smk_1933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2409[k] = f_12 * smi0_1508[k]
                    - f_13 * smi1_1508[k]
                    + f_3 * pc_x[k] * smk_1932[k];

        t_2410[k] = f_12 * smi0_1509[k]
                    - f_13 * smi1_1509[k]
                    + f_3 * pc_x[k] * smk_1933[k];

        t_2411[k] = f_15 * slk_1604[k]
                    + f_3 * pc_y[k] * smk_1928[k];
    }

#pragma omp simd aligned(t_2412, t_2413, t_2414, t_2415, t_2416, t_2417, pb_y, pc_x, pc_y, \
                         sll0_2007, sll1_2007, smk_1936, smk_1937, smk_1938, smk_1939, \
                         smk_1940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2412[k] = pb_y[k] * sll0_2007[k]
                    - f_14 * pc_y[k] * sll1_2007[k];

        t_2413[k] = f_3 * pc_x[k] * smk_1936[k];

        t_2414[k] = f_3 * pc_x[k] * smk_1937[k];

        t_2415[k] = f_3 * pc_x[k] * smk_1938[k];

        t_2416[k] = f_3 * pc_x[k] * smk_1939[k];

        t_2417[k] = f_3 * pc_x[k] * smk_1940[k];
    }

#pragma omp simd aligned(t_2418, t_2419, t_2420, t_2421, pb_y, pc_x, pc_y, sll0_2016, \
                         slk_1612, sll1_2016, smk_1941, smk_1942, \
                         smk_1943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2418[k] = f_3 * pc_x[k] * smk_1941[k];

        t_2419[k] = f_3 * pc_x[k] * smk_1942[k];

        t_2420[k] = f_3 * pc_x[k] * smk_1943[k];

        t_2421[k] = pb_y[k] * sll0_2016[k]
                    + f_21 * slk_1612[k]
                    - f_14 * pc_y[k] * sll1_2016[k];
    }

#pragma omp simd aligned(t_2422, t_2423, t_2424, pb_y, pc_y, pc_z, sll0_2018, sll0_2019, \
                         slk_1576, slk_1614, slk_1615, sll1_2018, sll1_2019, \
                         smk_1936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2422[k] = f_21 * slk_1576[k]
                    + f_3 * pc_z[k] * smk_1936[k];

        t_2423[k] = pb_y[k] * sll0_2018[k]
                    + f_20 * slk_1614[k]
                    - f_14 * pc_y[k] * sll1_2018[k];

        t_2424[k] = pb_y[k] * sll0_2019[k]
                    + f_19 * slk_1615[k]
                    - f_14 * pc_y[k] * sll1_2019[k];
    }

#pragma omp simd aligned(t_2425, t_2426, t_2427, pb_y, pc_y, sll0_2020, sll0_2021, sll0_2022, \
                         slk_1616, slk_1617, slk_1618, sll1_2020, sll1_2021, \
                         sll1_2022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2425[k] = pb_y[k] * sll0_2020[k]
                    + f_18 * slk_1616[k]
                    - f_14 * pc_y[k] * sll1_2020[k];

        t_2426[k] = pb_y[k] * sll0_2021[k]
                    + f_17 * slk_1617[k]
                    - f_14 * pc_y[k] * sll1_2021[k];

        t_2427[k] = pb_y[k] * sll0_2022[k]
                    + f_16 * slk_1618[k]
                    - f_14 * pc_y[k] * sll1_2022[k];
    }

#pragma omp simd aligned(t_2428, t_2429, t_2430, t_2431, pb_y, pc_x, pc_y, sll0_2024, \
                         slk_1619, sll1_2024, smi0_1512, smi1_1512, smk_1943, \
                         smk_1944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2428[k] = f_15 * slk_1619[k]
                    + f_3 * pc_y[k] * smk_1943[k];

        t_2429[k] = pb_y[k] * sll0_2024[k]
                    - f_14 * pc_y[k] * sll1_2024[k];

        t_2430[k] = f_1 * smi0_1512[k]
                    - f_2 * smi1_1512[k]
                    + f_3 * pc_x[k] * smk_1944[k];

        t_2431[k] = f_3 * pc_y[k] * smk_1944[k];
    }

#pragma omp simd aligned(t_2432, t_2433, t_2434, t_2435, pc_x, pc_y, pc_z, slk_1584, \
                         smi0_1515, smi0_1517, smi1_1515, smi1_1517, smk_1944, smk_1946, \
                         smk_1947, smk_1949 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2432[k] = f_0 * slk_1584[k]
                    + f_3 * pc_z[k] * smk_1944[k];

        t_2433[k] = f_4 * smi0_1515[k]
                    - f_5 * smi1_1515[k]
                    + f_3 * pc_x[k] * smk_1947[k];

        t_2434[k] = f_3 * pc_y[k] * smk_1946[k];

        t_2435[k] = f_4 * smi0_1517[k]
                    - f_5 * smi1_1517[k]
                    + f_3 * pc_x[k] * smk_1949[k];
    }

#pragma omp simd aligned(t_2436, t_2437, t_2438, t_2439, pc_x, pc_y, pc_z, slk_1587, \
                         smi0_1518, smi0_1521, smi1_1518, smi1_1521, smk_1947, smk_1949, \
                         smk_1950, smk_1953 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2436[k] = f_6 * smi0_1518[k]
                    - f_7 * smi1_1518[k]
                    + f_3 * pc_x[k] * smk_1950[k];

        t_2437[k] = f_0 * slk_1587[k]
                    + f_3 * pc_z[k] * smk_1947[k];

        t_2438[k] = f_3 * pc_y[k] * smk_1949[k];

        t_2439[k] = f_6 * smi0_1521[k]
                    - f_7 * smi1_1521[k]
                    + f_3 * pc_x[k] * smk_1953[k];
    }

#pragma omp simd aligned(t_2440, t_2441, t_2442, t_2443, pc_x, pc_y, pc_z, slk_1590, \
                         smi0_1522, smi0_1524, smi1_1522, smi1_1524, smk_1950, smk_1953, \
                         smk_1954, smk_1956 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2440[k] = f_8 * smi0_1522[k]
                    - f_9 * smi1_1522[k]
                    + f_3 * pc_x[k] * smk_1954[k];

        t_2441[k] = f_0 * slk_1590[k]
                    + f_3 * pc_z[k] * smk_1950[k];

        t_2442[k] = f_8 * smi0_1524[k]
                    - f_9 * smi1_1524[k]
                    + f_3 * pc_x[k] * smk_1956[k];

        t_2443[k] = f_3 * pc_y[k] * smk_1953[k];
    }

#pragma omp simd aligned(t_2444, t_2445, t_2446, pc_x, pc_z, slk_1594, smi0_1526, smi0_1527, \
                         smi1_1526, smi1_1527, smk_1954, smk_1958, \
                         smk_1959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2444[k] = f_8 * smi0_1526[k]
                    - f_9 * smi1_1526[k]
                    + f_3 * pc_x[k] * smk_1958[k];

        t_2445[k] = f_10 * smi0_1527[k]
                    - f_11 * smi1_1527[k]
                    + f_3 * pc_x[k] * smk_1959[k];

        t_2446[k] = f_0 * slk_1594[k]
                    + f_3 * pc_z[k] * smk_1954[k];
    }

#pragma omp simd aligned(t_2447, t_2448, t_2449, t_2450, pc_x, pc_y, smi0_1529, smi0_1530, \
                         smi0_1532, smi1_1529, smi1_1530, smi1_1532, smk_1958, smk_1961, \
                         smk_1962, smk_1964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2447[k] = f_10 * smi0_1529[k]
                    - f_11 * smi1_1529[k]
                    + f_3 * pc_x[k] * smk_1961[k];

        t_2448[k] = f_10 * smi0_1530[k]
                    - f_11 * smi1_1530[k]
                    + f_3 * pc_x[k] * smk_1962[k];

        t_2449[k] = f_3 * pc_y[k] * smk_1958[k];

        t_2450[k] = f_10 * smi0_1532[k]
                    - f_11 * smi1_1532[k]
                    + f_3 * pc_x[k] * smk_1964[k];
    }

#pragma omp simd aligned(t_2451, t_2452, t_2453, pc_x, pc_z, slk_1599, smi0_1533, smi0_1535, \
                         smi1_1533, smi1_1535, smk_1959, smk_1965, \
                         smk_1967 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2451[k] = f_12 * smi0_1533[k]
                    - f_13 * smi1_1533[k]
                    + f_3 * pc_x[k] * smk_1965[k];

        t_2452[k] = f_0 * slk_1599[k]
                    + f_3 * pc_z[k] * smk_1959[k];

        t_2453[k] = f_12 * smi0_1535[k]
                    - f_13 * smi1_1535[k]
                    + f_3 * pc_x[k] * smk_1967[k];
    }

#pragma omp simd aligned(t_2454, t_2455, t_2456, t_2457, pc_x, pc_y, smi0_1536, smi0_1537, \
                         smi0_1539, smi1_1536, smi1_1537, smi1_1539, smk_1964, smk_1968, \
                         smk_1969, smk_1971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2454[k] = f_12 * smi0_1536[k]
                    - f_13 * smi1_1536[k]
                    + f_3 * pc_x[k] * smk_1968[k];

        t_2455[k] = f_12 * smi0_1537[k]
                    - f_13 * smi1_1537[k]
                    + f_3 * pc_x[k] * smk_1969[k];

        t_2456[k] = f_3 * pc_y[k] * smk_1964[k];

        t_2457[k] = f_12 * smi0_1539[k]
                    - f_13 * smi1_1539[k]
                    + f_3 * pc_x[k] * smk_1971[k];
    }

#pragma omp simd aligned(t_2458, t_2459, t_2460, t_2461, t_2462, t_2463, t_2464, pc_x, \
                         smk_1972, smk_1973, smk_1974, smk_1975, smk_1976, smk_1977, \
                         smk_1978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2458[k] = f_3 * pc_x[k] * smk_1972[k];

        t_2459[k] = f_3 * pc_x[k] * smk_1973[k];

        t_2460[k] = f_3 * pc_x[k] * smk_1974[k];

        t_2461[k] = f_3 * pc_x[k] * smk_1975[k];

        t_2462[k] = f_3 * pc_x[k] * smk_1976[k];

        t_2463[k] = f_3 * pc_x[k] * smk_1977[k];

        t_2464[k] = f_3 * pc_x[k] * smk_1978[k];
    }

#pragma omp simd aligned(t_2465, t_2466, t_2467, t_2468, pc_x, pc_y, pc_z, slk_1612, \
                         smi0_1533, smi0_1535, smi1_1533, smi1_1535, smk_1972, smk_1974, \
                         smk_1979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2465[k] = f_3 * pc_x[k] * smk_1979[k];

        t_2466[k] = f_1 * smi0_1533[k]
                    - f_2 * smi1_1533[k]
                    + f_3 * pc_y[k] * smk_1972[k];

        t_2467[k] = f_0 * slk_1612[k]
                    + f_3 * pc_z[k] * smk_1972[k];

        t_2468[k] = f_4 * smi0_1535[k]
                    - f_5 * smi1_1535[k]
                    + f_3 * pc_y[k] * smk_1974[k];
    }

#pragma omp simd aligned(t_2469, t_2470, t_2471, pc_y, smi0_1536, smi0_1537, smi0_1538, \
                         smi1_1536, smi1_1537, smi1_1538, smk_1975, smk_1976, \
                         smk_1977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2469[k] = f_6 * smi0_1536[k]
                    - f_7 * smi1_1536[k]
                    + f_3 * pc_y[k] * smk_1975[k];

        t_2470[k] = f_8 * smi0_1537[k]
                    - f_9 * smi1_1537[k]
                    + f_3 * pc_y[k] * smk_1976[k];

        t_2471[k] = f_10 * smi0_1538[k]
                    - f_11 * smi1_1538[k]
                    + f_3 * pc_y[k] * smk_1977[k];
    }

#pragma omp simd aligned(t_2472, t_2473, t_2474, pc_y, pc_z, slk_1619, smi0_1539, smi1_1539, \
                         smk_1978, smk_1979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2472[k] = f_12 * smi0_1539[k]
                    - f_13 * smi1_1539[k]
                    + f_3 * pc_y[k] * smk_1978[k];

        t_2473[k] = f_3 * pc_y[k] * smk_1979[k];

        t_2474[k] = f_0 * slk_1619[k]
                    + f_1 * smi0_1539[k]
                    - f_2 * smi1_1539[k]
                    + f_3 * pc_z[k] * smk_1979[k];
    }
}

auto
compute_prim_sml_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sll0, const size_t slk,
                                                   const size_t sll1, const size_t smi0,
                                                   const size_t smi1, const size_t smk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sml_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, sll0, slk,
                                                              sll1, smi0, smi1, smk, ncols,
                                                              gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smi0, smi1, smk,
                                                               ncols, gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smi0, smi1, smk,
                                                               ncols, gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece12(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smi0, smi1, smk,
                                                               ncols, gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece13(buffer, target, pc, slk, smi0,
                                                               smi1, smk, ncols, gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece14(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smi0, smi1, smk,
                                                               ncols, gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece15(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smk, ncols, gamma, p,
                                                               q);

    compute_prim_sml_three_center_electron_repulsion_0_piece16(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smk, ncols, gamma, p,
                                                               q);

    compute_prim_sml_three_center_electron_repulsion_0_piece17(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smk, ncols, gamma, p,
                                                               q);

    compute_prim_sml_three_center_electron_repulsion_0_piece18(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smi0, smi1, smk,
                                                               ncols, gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece19(buffer, target, pc, slk, smi0,
                                                               smi1, smk, ncols, gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece20(buffer, target, pc, slk, smi0,
                                                               smi1, smk, ncols, gamma, p, q);

    compute_prim_sml_three_center_electron_repulsion_0_piece21(buffer, target, pb, pc, sll0,
                                                               slk, sll1, smi0, smi1, smk,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
