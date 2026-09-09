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


#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_0 = buffer.data(msi0 + 0);
    const auto *msi0_3 = buffer.data(msi0 + 3);
    const auto *msi0_5 = buffer.data(msi0 + 5);
    const auto *msi0_6 = buffer.data(msi0 + 6);
    const auto *msi0_9 = buffer.data(msi0 + 9);
    const auto *msi0_10 = buffer.data(msi0 + 10);
    const auto *msi0_14 = buffer.data(msi0 + 14);
    const auto *msi0_21 = buffer.data(msi0 + 21);
    const auto *msi0_27 = buffer.data(msi0 + 27);
    const auto *msi0_31 = buffer.data(msi0 + 31);
    const auto *msi0_34 = buffer.data(msi0 + 34);
    const auto *msi0_38 = buffer.data(msi0 + 38);
    const auto *msi0_56 = buffer.data(msi0 + 56);
    const auto *msi0_61 = buffer.data(msi0 + 61);
    const auto *msi0_65 = buffer.data(msi0 + 65);
    const auto *msi0_68 = buffer.data(msi0 + 68);
    const auto *msi0_70 = buffer.data(msi0 + 70);

    const auto *msh_0 = buffer.data(msh + 0);
    const auto *msh_1 = buffer.data(msh + 1);
    const auto *msh_2 = buffer.data(msh + 2);
    const auto *msh_3 = buffer.data(msh + 3);
    const auto *msh_5 = buffer.data(msh + 5);
    const auto *msh_6 = buffer.data(msh + 6);
    const auto *msh_9 = buffer.data(msh + 9);
    const auto *msh_15 = buffer.data(msh + 15);
    const auto *msh_17 = buffer.data(msh + 17);
    const auto *msh_18 = buffer.data(msh + 18);
    const auto *msh_20 = buffer.data(msh + 20);
    const auto *msh_21 = buffer.data(msh + 21);
    const auto *msh_24 = buffer.data(msh + 24);
    const auto *msh_26 = buffer.data(msh + 26);
    const auto *msh_27 = buffer.data(msh + 27);
    const auto *msh_30 = buffer.data(msh + 30);
    const auto *msh_36 = buffer.data(msh + 36);
    const auto *msh_38 = buffer.data(msh + 38);
    const auto *msh_39 = buffer.data(msh + 39);
    const auto *msh_40 = buffer.data(msh + 40);
    const auto *msh_41 = buffer.data(msh + 41);
    const auto *msh_42 = buffer.data(msh + 42);
    const auto *msh_44 = buffer.data(msh + 44);
    const auto *msh_47 = buffer.data(msh + 47);
    const auto *msh_50 = buffer.data(msh + 50);
    const auto *msh_51 = buffer.data(msh + 51);
    const auto *msh_57 = buffer.data(msh + 57);
    const auto *msh_58 = buffer.data(msh + 58);
    const auto *msh_59 = buffer.data(msh + 59);
    const auto *msh_60 = buffer.data(msh + 60);
    const auto *msh_62 = buffer.data(msh + 62);
    const auto *msh_63 = buffer.data(msh + 63);
    const auto *msh_66 = buffer.data(msh + 66);
    const auto *msh_69 = buffer.data(msh + 69);
    const auto *msh_73 = buffer.data(msh + 73);
    const auto *msh_78 = buffer.data(msh + 78);
    const auto *msh_80 = buffer.data(msh + 80);
    const auto *msh_81 = buffer.data(msh + 81);
    const auto *msh_82 = buffer.data(msh + 82);
    const auto *msh_83 = buffer.data(msh + 83);
    const auto *msh_99 = buffer.data(msh + 99);

    const auto *msi1_0 = buffer.data(msi1 + 0);
    const auto *msi1_3 = buffer.data(msi1 + 3);
    const auto *msi1_5 = buffer.data(msi1 + 5);
    const auto *msi1_6 = buffer.data(msi1 + 6);
    const auto *msi1_9 = buffer.data(msi1 + 9);
    const auto *msi1_10 = buffer.data(msi1 + 10);
    const auto *msi1_14 = buffer.data(msi1 + 14);
    const auto *msi1_21 = buffer.data(msi1 + 21);
    const auto *msi1_27 = buffer.data(msi1 + 27);
    const auto *msi1_31 = buffer.data(msi1 + 31);
    const auto *msi1_34 = buffer.data(msi1 + 34);
    const auto *msi1_38 = buffer.data(msi1 + 38);
    const auto *msi1_56 = buffer.data(msi1 + 56);
    const auto *msi1_61 = buffer.data(msi1 + 61);
    const auto *msi1_65 = buffer.data(msi1 + 65);
    const auto *msi1_68 = buffer.data(msi1 + 68);
    const auto *msi1_70 = buffer.data(msi1 + 70);

    const auto *nsg0_0 = buffer.data(nsg0 + 0);
    const auto *nsg0_1 = buffer.data(nsg0 + 1);
    const auto *nsg0_2 = buffer.data(nsg0 + 2);
    const auto *nsg0_3 = buffer.data(nsg0 + 3);
    const auto *nsg0_5 = buffer.data(nsg0 + 5);
    const auto *nsg0_10 = buffer.data(nsg0 + 10);
    const auto *nsg0_12 = buffer.data(nsg0 + 12);
    const auto *nsg0_13 = buffer.data(nsg0 + 13);
    const auto *nsg0_14 = buffer.data(nsg0 + 14);
    const auto *nsg0_18 = buffer.data(nsg0 + 18);
    const auto *nsg0_25 = buffer.data(nsg0 + 25);
    const auto *nsg0_26 = buffer.data(nsg0 + 26);
    const auto *nsg0_27 = buffer.data(nsg0 + 27);
    const auto *nsg0_32 = buffer.data(nsg0 + 32);
    const auto *nsg0_34 = buffer.data(nsg0 + 34);
    const auto *nsg0_35 = buffer.data(nsg0 + 35);
    const auto *nsg0_41 = buffer.data(nsg0 + 41);
    const auto *nsg0_42 = buffer.data(nsg0 + 42);
    const auto *nsg0_43 = buffer.data(nsg0 + 43);
    const auto *nsg0_44 = buffer.data(nsg0 + 44);
    const auto *nsg0_45 = buffer.data(nsg0 + 45);
    const auto *nsg0_47 = buffer.data(nsg0 + 47);
    const auto *nsg0_48 = buffer.data(nsg0 + 48);
    const auto *nsg0_50 = buffer.data(nsg0 + 50);
    const auto *nsg0_51 = buffer.data(nsg0 + 51);
    const auto *nsg0_55 = buffer.data(nsg0 + 55);
    const auto *nsg0_56 = buffer.data(nsg0 + 56);
    const auto *nsg0_57 = buffer.data(nsg0 + 57);
    const auto *nsg0_59 = buffer.data(nsg0 + 59);

    const auto *nsg1_0 = buffer.data(nsg1 + 0);
    const auto *nsg1_1 = buffer.data(nsg1 + 1);
    const auto *nsg1_2 = buffer.data(nsg1 + 2);
    const auto *nsg1_3 = buffer.data(nsg1 + 3);
    const auto *nsg1_5 = buffer.data(nsg1 + 5);
    const auto *nsg1_10 = buffer.data(nsg1 + 10);
    const auto *nsg1_12 = buffer.data(nsg1 + 12);
    const auto *nsg1_13 = buffer.data(nsg1 + 13);
    const auto *nsg1_14 = buffer.data(nsg1 + 14);
    const auto *nsg1_18 = buffer.data(nsg1 + 18);
    const auto *nsg1_25 = buffer.data(nsg1 + 25);
    const auto *nsg1_26 = buffer.data(nsg1 + 26);
    const auto *nsg1_27 = buffer.data(nsg1 + 27);
    const auto *nsg1_32 = buffer.data(nsg1 + 32);
    const auto *nsg1_34 = buffer.data(nsg1 + 34);
    const auto *nsg1_35 = buffer.data(nsg1 + 35);
    const auto *nsg1_41 = buffer.data(nsg1 + 41);
    const auto *nsg1_42 = buffer.data(nsg1 + 42);
    const auto *nsg1_43 = buffer.data(nsg1 + 43);
    const auto *nsg1_44 = buffer.data(nsg1 + 44);
    const auto *nsg1_45 = buffer.data(nsg1 + 45);
    const auto *nsg1_47 = buffer.data(nsg1 + 47);
    const auto *nsg1_48 = buffer.data(nsg1 + 48);
    const auto *nsg1_50 = buffer.data(nsg1 + 50);
    const auto *nsg1_51 = buffer.data(nsg1 + 51);
    const auto *nsg1_55 = buffer.data(nsg1 + 55);
    const auto *nsg1_56 = buffer.data(nsg1 + 56);
    const auto *nsg1_57 = buffer.data(nsg1 + 57);
    const auto *nsg1_59 = buffer.data(nsg1 + 59);

    const auto *nsh_0 = buffer.data(nsh + 0);
    const auto *nsh_1 = buffer.data(nsh + 1);
    const auto *nsh_2 = buffer.data(nsh + 2);
    const auto *nsh_3 = buffer.data(nsh + 3);
    const auto *nsh_5 = buffer.data(nsh + 5);
    const auto *nsh_6 = buffer.data(nsh + 6);
    const auto *nsh_8 = buffer.data(nsh + 8);
    const auto *nsh_9 = buffer.data(nsh + 9);
    const auto *nsh_10 = buffer.data(nsh + 10);
    const auto *nsh_14 = buffer.data(nsh + 14);
    const auto *nsh_15 = buffer.data(nsh + 15);
    const auto *nsh_17 = buffer.data(nsh + 17);
    const auto *nsh_18 = buffer.data(nsh + 18);
    const auto *nsh_19 = buffer.data(nsh + 19);
    const auto *nsh_20 = buffer.data(nsh + 20);
    const auto *nsh_21 = buffer.data(nsh + 21);
    const auto *nsh_22 = buffer.data(nsh + 22);
    const auto *nsh_24 = buffer.data(nsh + 24);
    const auto *nsh_26 = buffer.data(nsh + 26);
    const auto *nsh_27 = buffer.data(nsh + 27);
    const auto *nsh_28 = buffer.data(nsh + 28);
    const auto *nsh_30 = buffer.data(nsh + 30);
    const auto *nsh_31 = buffer.data(nsh + 31);
    const auto *nsh_36 = buffer.data(nsh + 36);
    const auto *nsh_37 = buffer.data(nsh + 37);
    const auto *nsh_38 = buffer.data(nsh + 38);
    const auto *nsh_39 = buffer.data(nsh + 39);
    const auto *nsh_40 = buffer.data(nsh + 40);
    const auto *nsh_41 = buffer.data(nsh + 41);
    const auto *nsh_42 = buffer.data(nsh + 42);
    const auto *nsh_44 = buffer.data(nsh + 44);
    const auto *nsh_46 = buffer.data(nsh + 46);
    const auto *nsh_47 = buffer.data(nsh + 47);
    const auto *nsh_49 = buffer.data(nsh + 49);
    const auto *nsh_50 = buffer.data(nsh + 50);
    const auto *nsh_51 = buffer.data(nsh + 51);
    const auto *nsh_56 = buffer.data(nsh + 56);
    const auto *nsh_57 = buffer.data(nsh + 57);
    const auto *nsh_58 = buffer.data(nsh + 58);
    const auto *nsh_59 = buffer.data(nsh + 59);
    const auto *nsh_60 = buffer.data(nsh + 60);
    const auto *nsh_61 = buffer.data(nsh + 61);
    const auto *nsh_62 = buffer.data(nsh + 62);
    const auto *nsh_63 = buffer.data(nsh + 63);
    const auto *nsh_64 = buffer.data(nsh + 64);
    const auto *nsh_65 = buffer.data(nsh + 65);
    const auto *nsh_66 = buffer.data(nsh + 66);
    const auto *nsh_68 = buffer.data(nsh + 68);
    const auto *nsh_69 = buffer.data(nsh + 69);
    const auto *nsh_70 = buffer.data(nsh + 70);
    const auto *nsh_72 = buffer.data(nsh + 72);
    const auto *nsh_73 = buffer.data(nsh + 73);
    const auto *nsh_78 = buffer.data(nsh + 78);
    const auto *nsh_79 = buffer.data(nsh + 79);
    const auto *nsh_80 = buffer.data(nsh + 80);
    const auto *nsh_81 = buffer.data(nsh + 81);
    const auto *nsh_82 = buffer.data(nsh + 82);
    const auto *nsh_83 = buffer.data(nsh + 83);
    const auto *nsh_84 = buffer.data(nsh + 84);
    const auto *nsh_86 = buffer.data(nsh + 86);
    const auto *nsh_87 = buffer.data(nsh + 87);
    const auto *nsh_89 = buffer.data(nsh + 89);
    const auto *nsh_90 = buffer.data(nsh + 90);
    const auto *nsh_93 = buffer.data(nsh + 93);
    const auto *nsh_99 = buffer.data(nsh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, msh_0, nsg0_0, \
                         nsg1_0, nsh_0, nsh_1, nsh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * msh_0[k]
                 + f_1 * nsg0_0[k]
                 - f_2 * nsg1_0[k]
                 + f_3 * pc_x[k] * nsh_0[k];

        t_1[k] = f_3 * pc_y[k] * nsh_0[k];

        t_2[k] = f_3 * pc_z[k] * nsh_0[k];

        t_3[k] = f_4 * nsg0_0[k]
                 - f_5 * nsg1_0[k]
                 + f_3 * pc_y[k] * nsh_1[k];

        t_4[k] = f_3 * pc_y[k] * nsh_2[k];

        t_5[k] = f_4 * nsg0_0[k]
                 - f_5 * nsg1_0[k]
                 + f_3 * pc_z[k] * nsh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, nsg0_1, nsg0_2, nsg0_3, nsg1_1, \
                         nsg1_2, nsg1_3, nsh_3, nsh_5, nsh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * nsg0_1[k]
                 - f_7 * nsg1_1[k]
                 + f_3 * pc_y[k] * nsh_3[k];

        t_7[k] = f_3 * pc_z[k] * nsh_3[k];

        t_8[k] = f_3 * pc_y[k] * nsh_5[k];

        t_9[k] = f_6 * nsg0_2[k]
                 - f_7 * nsg1_2[k]
                 + f_3 * pc_z[k] * nsh_5[k];

        t_10[k] = f_8 * nsg0_3[k]
                  - f_9 * nsg1_3[k]
                  + f_3 * pc_y[k] * nsh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, msh_15, nsg0_5, \
                         nsg1_5, nsh_6, nsh_8, nsh_9, nsh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * nsh_6[k];

        t_12[k] = f_4 * nsg0_5[k]
                  - f_5 * nsg1_5[k]
                  + f_3 * pc_y[k] * nsh_8[k];

        t_13[k] = f_3 * pc_y[k] * nsh_9[k];

        t_14[k] = f_8 * nsg0_5[k]
                  - f_9 * nsg1_5[k]
                  + f_3 * pc_z[k] * nsh_9[k];

        t_15[k] = f_0 * msh_15[k]
                  + f_3 * pc_x[k] * nsh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, msh_17, msh_18, \
                         msh_20, nsh_10, nsh_14, nsh_17, nsh_18, \
                         nsh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * nsh_10[k];

        t_17[k] = f_0 * msh_17[k]
                  + f_3 * pc_x[k] * nsh_17[k];

        t_18[k] = f_0 * msh_18[k]
                  + f_3 * pc_x[k] * nsh_18[k];

        t_19[k] = f_3 * pc_y[k] * nsh_14[k];

        t_20[k] = f_0 * msh_20[k]
                  + f_3 * pc_x[k] * nsh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, nsg0_10, nsg0_12, nsg0_13, \
                         nsg1_10, nsg1_12, nsg1_13, nsh_15, nsh_17, \
                         nsh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * nsg0_10[k]
                  - f_2 * nsg1_10[k]
                  + f_3 * pc_y[k] * nsh_15[k];

        t_22[k] = f_3 * pc_z[k] * nsh_15[k];

        t_23[k] = f_8 * nsg0_12[k]
                  - f_9 * nsg1_12[k]
                  + f_3 * pc_y[k] * nsh_17[k];

        t_24[k] = f_6 * nsg0_13[k]
                  - f_7 * nsg1_13[k]
                  + f_3 * pc_y[k] * nsh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, msi0_0, msh_0, \
                         msi1_0, nsg0_14, nsg1_14, nsh_19, nsh_20, \
                         nsh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * nsg0_14[k]
                  - f_5 * nsg1_14[k]
                  + f_3 * pc_y[k] * nsh_19[k];

        t_26[k] = f_3 * pc_y[k] * nsh_20[k];

        t_27[k] = f_1 * nsg0_14[k]
                  - f_2 * nsg1_14[k]
                  + f_3 * pc_z[k] * nsh_20[k];

        t_28[k] = pa_y[k] * msi0_0[k]
                  - f_10 * pc_y[k] * msi1_0[k];

        t_29[k] = f_11 * msh_0[k]
                  + f_3 * pc_y[k] * nsh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, msi0_3, msi0_5, msh_1, \
                         msi1_3, msi1_5, nsh_21, nsh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * nsh_21[k];

        t_31[k] = pa_y[k] * msi0_3[k]
                  + f_12 * msh_1[k]
                  - f_10 * pc_y[k] * msi1_3[k];

        t_32[k] = f_3 * pc_z[k] * nsh_22[k];

        t_33[k] = pa_y[k] * msi0_5[k]
                  - f_10 * pc_y[k] * msi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, msi0_6, msi0_9, msh_3, \
                         msh_5, msi1_6, msi1_9, nsh_24, nsh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * msi0_6[k]
                  + f_13 * msh_3[k]
                  - f_10 * pc_y[k] * msi1_6[k];

        t_35[k] = f_3 * pc_z[k] * nsh_24[k];

        t_36[k] = f_11 * msh_5[k]
                  + f_3 * pc_y[k] * nsh_26[k];

        t_37[k] = pa_y[k] * msi0_9[k]
                  - f_10 * pc_y[k] * msi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, msi0_10, msh_6, msh_9, \
                         msi1_10, nsg0_18, nsg1_18, nsh_27, nsh_28, \
                         nsh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * msi0_10[k]
                  + f_14 * msh_6[k]
                  - f_10 * pc_y[k] * msi1_10[k];

        t_39[k] = f_3 * pc_z[k] * nsh_27[k];

        t_40[k] = f_4 * nsg0_18[k]
                  - f_5 * nsg1_18[k]
                  + f_3 * pc_z[k] * nsh_28[k];

        t_41[k] = f_11 * msh_9[k]
                  + f_3 * pc_y[k] * nsh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, msi0_14, msh_36, \
                         msh_38, msi1_14, nsh_31, nsh_36, nsh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * msi0_14[k]
                  - f_10 * pc_y[k] * msi1_14[k];

        t_43[k] = f_15 * msh_36[k]
                  + f_3 * pc_x[k] * nsh_36[k];

        t_44[k] = f_3 * pc_z[k] * nsh_31[k];

        t_45[k] = f_15 * msh_38[k]
                  + f_3 * pc_x[k] * nsh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, msh_15, msh_39, msh_40, msh_41, \
                         nsg0_25, nsg1_25, nsh_36, nsh_39, nsh_40, \
                         nsh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * msh_39[k]
                  + f_3 * pc_x[k] * nsh_39[k];

        t_47[k] = f_15 * msh_40[k]
                  + f_3 * pc_x[k] * nsh_40[k];

        t_48[k] = f_15 * msh_41[k]
                  + f_3 * pc_x[k] * nsh_41[k];

        t_49[k] = f_11 * msh_15[k]
                  + f_1 * nsg0_25[k]
                  - f_2 * nsg1_25[k]
                  + f_3 * pc_y[k] * nsh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, nsg0_25, nsg0_26, nsg0_27, nsg1_25, \
                         nsg1_26, nsg1_27, nsh_36, nsh_37, nsh_38, \
                         nsh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * nsh_36[k];

        t_51[k] = f_4 * nsg0_25[k]
                  - f_5 * nsg1_25[k]
                  + f_3 * pc_z[k] * nsh_37[k];

        t_52[k] = f_6 * nsg0_26[k]
                  - f_7 * nsg1_26[k]
                  + f_3 * pc_z[k] * nsh_38[k];

        t_53[k] = f_8 * nsg0_27[k]
                  - f_9 * nsg1_27[k]
                  + f_3 * pc_z[k] * nsh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, msi0_0, msi0_27, \
                         msh_20, msi1_0, msi1_27, nsh_41, nsh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * msh_20[k]
                  + f_3 * pc_y[k] * nsh_41[k];

        t_55[k] = pa_y[k] * msi0_27[k]
                  - f_10 * pc_y[k] * msi1_27[k];

        t_56[k] = pa_z[k] * msi0_0[k]
                  - f_10 * pc_z[k] * msi1_0[k];

        t_57[k] = f_3 * pc_y[k] * nsh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, msi0_3, msi0_5, msh_0, \
                         msh_2, msi1_3, msi1_5, nsh_42, nsh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * msh_0[k]
                  + f_3 * pc_z[k] * nsh_42[k];

        t_59[k] = pa_z[k] * msi0_3[k]
                  - f_10 * pc_z[k] * msi1_3[k];

        t_60[k] = f_3 * pc_y[k] * nsh_44[k];

        t_61[k] = pa_z[k] * msi0_5[k]
                  + f_12 * msh_2[k]
                  - f_10 * pc_z[k] * msi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, msi0_6, msi0_9, msh_5, \
                         msi1_6, msi1_9, nsg0_32, nsg1_32, nsh_46, \
                         nsh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * msi0_6[k]
                  - f_10 * pc_z[k] * msi1_6[k];

        t_63[k] = f_4 * nsg0_32[k]
                  - f_5 * nsg1_32[k]
                  + f_3 * pc_y[k] * nsh_46[k];

        t_64[k] = f_3 * pc_y[k] * nsh_47[k];

        t_65[k] = pa_z[k] * msi0_9[k]
                  + f_13 * msh_5[k]
                  - f_10 * pc_z[k] * msi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, msi0_10, msi1_10, nsg0_34, \
                         nsg0_35, nsg1_34, nsg1_35, nsh_49, nsh_50, \
                         nsh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * msi0_10[k]
                  - f_10 * pc_z[k] * msi1_10[k];

        t_67[k] = f_6 * nsg0_34[k]
                  - f_7 * nsg1_34[k]
                  + f_3 * pc_y[k] * nsh_49[k];

        t_68[k] = f_4 * nsg0_35[k]
                  - f_5 * nsg1_35[k]
                  + f_3 * pc_y[k] * nsh_50[k];

        t_69[k] = f_3 * pc_y[k] * nsh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, msi0_14, msh_9, msh_57, \
                         msh_58, msh_59, msi1_14, nsh_57, nsh_58, \
                         nsh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * msi0_14[k]
                  + f_14 * msh_9[k]
                  - f_10 * pc_z[k] * msi1_14[k];

        t_71[k] = f_15 * msh_57[k]
                  + f_3 * pc_x[k] * nsh_57[k];

        t_72[k] = f_15 * msh_58[k]
                  + f_3 * pc_x[k] * nsh_58[k];

        t_73[k] = f_15 * msh_59[k]
                  + f_3 * pc_x[k] * nsh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, msi0_21, msh_60, \
                         msh_62, msi1_21, nsh_56, nsh_60, nsh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_15 * msh_60[k]
                  + f_3 * pc_x[k] * nsh_60[k];

        t_75[k] = f_3 * pc_y[k] * nsh_56[k];

        t_76[k] = f_15 * msh_62[k]
                  + f_3 * pc_x[k] * nsh_62[k];

        t_77[k] = pa_z[k] * msi0_21[k]
                  - f_10 * pc_z[k] * msi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, nsg0_41, nsg0_42, nsg0_43, nsg1_41, nsg1_42, \
                         nsg1_43, nsh_58, nsh_59, nsh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * nsg0_41[k]
                  - f_17 * nsg1_41[k]
                  + f_3 * pc_y[k] * nsh_58[k];

        t_79[k] = f_8 * nsg0_42[k]
                  - f_9 * nsg1_42[k]
                  + f_3 * pc_y[k] * nsh_59[k];

        t_80[k] = f_6 * nsg0_43[k]
                  - f_7 * nsg1_43[k]
                  + f_3 * pc_y[k] * nsh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, msh_20, msh_63, nsg0_44, \
                         nsg0_45, nsg1_44, nsg1_45, nsh_61, nsh_62, \
                         nsh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * nsg0_44[k]
                  - f_5 * nsg1_44[k]
                  + f_3 * pc_y[k] * nsh_61[k];

        t_82[k] = f_3 * pc_y[k] * nsh_62[k];

        t_83[k] = f_11 * msh_20[k]
                  + f_1 * nsg0_44[k]
                  - f_2 * nsg1_44[k]
                  + f_3 * pc_z[k] * nsh_62[k];

        t_84[k] = f_18 * msh_63[k]
                  + f_1 * nsg0_45[k]
                  - f_2 * nsg1_45[k]
                  + f_3 * pc_x[k] * nsh_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, msh_21, msh_66, nsg0_48, \
                         nsg1_48, nsh_63, nsh_64, nsh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * msh_21[k]
                  + f_3 * pc_y[k] * nsh_63[k];

        t_86[k] = f_3 * pc_z[k] * nsh_63[k];

        t_87[k] = f_18 * msh_66[k]
                  + f_8 * nsg0_48[k]
                  - f_9 * nsg1_48[k]
                  + f_3 * pc_x[k] * nsh_66[k];

        t_88[k] = f_3 * pc_z[k] * nsh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, msh_69, nsg0_45, nsg0_51, nsg1_45, \
                         nsg1_51, nsh_65, nsh_66, nsh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * nsg0_45[k]
                  - f_5 * nsg1_45[k]
                  + f_3 * pc_z[k] * nsh_65[k];

        t_90[k] = f_18 * msh_69[k]
                  + f_6 * nsg0_51[k]
                  - f_7 * nsg1_51[k]
                  + f_3 * pc_x[k] * nsh_69[k];

        t_91[k] = f_3 * pc_z[k] * nsh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, msh_26, msh_73, nsg0_47, \
                         nsg0_55, nsg1_47, nsg1_55, nsh_68, nsh_69, \
                         nsh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * msh_26[k]
                  + f_3 * pc_y[k] * nsh_68[k];

        t_93[k] = f_6 * nsg0_47[k]
                  - f_7 * nsg1_47[k]
                  + f_3 * pc_z[k] * nsh_68[k];

        t_94[k] = f_18 * msh_73[k]
                  + f_4 * nsg0_55[k]
                  - f_5 * nsg1_55[k]
                  + f_3 * pc_x[k] * nsh_73[k];

        t_95[k] = f_3 * pc_z[k] * nsh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, msh_30, msh_78, nsg0_48, \
                         nsg0_50, nsg1_48, nsg1_50, nsh_70, nsh_72, \
                         nsh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * nsg0_48[k]
                  - f_5 * nsg1_48[k]
                  + f_3 * pc_z[k] * nsh_70[k];

        t_97[k] = f_12 * msh_30[k]
                  + f_3 * pc_y[k] * nsh_72[k];

        t_98[k] = f_8 * nsg0_50[k]
                  - f_9 * nsg1_50[k]
                  + f_3 * pc_z[k] * nsh_72[k];

        t_99[k] = f_18 * msh_78[k]
                  + f_3 * pc_x[k] * nsh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, msh_80, msh_81, \
                         msh_82, msh_83, nsh_73, nsh_80, nsh_81, nsh_82, \
                         nsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * nsh_73[k];

        t_101[k] = f_18 * msh_80[k]
                   + f_3 * pc_x[k] * nsh_80[k];

        t_102[k] = f_18 * msh_81[k]
                   + f_3 * pc_x[k] * nsh_81[k];

        t_103[k] = f_18 * msh_82[k]
                   + f_3 * pc_x[k] * nsh_82[k];

        t_104[k] = f_18 * msh_83[k]
                   + f_3 * pc_x[k] * nsh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, msh_36, nsg0_55, nsg0_56, \
                         nsg1_55, nsg1_56, nsh_78, nsh_79, nsh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * msh_36[k]
                   + f_1 * nsg0_55[k]
                   - f_2 * nsg1_55[k]
                   + f_3 * pc_y[k] * nsh_78[k];

        t_106[k] = f_3 * pc_z[k] * nsh_78[k];

        t_107[k] = f_4 * nsg0_55[k]
                   - f_5 * nsg1_55[k]
                   + f_3 * pc_z[k] * nsh_79[k];

        t_108[k] = f_6 * nsg0_56[k]
                   - f_7 * nsg1_56[k]
                   + f_3 * pc_z[k] * nsh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, msi0_56, msh_41, \
                         msi1_56, nsg0_57, nsg0_59, nsg1_57, nsg1_59, nsh_81, \
                         nsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * nsg0_57[k]
                   - f_9 * nsg1_57[k]
                   + f_3 * pc_z[k] * nsh_81[k];

        t_110[k] = f_12 * msh_41[k]
                   + f_3 * pc_y[k] * nsh_83[k];

        t_111[k] = f_1 * nsg0_59[k]
                   - f_2 * nsg1_59[k]
                   + f_3 * pc_z[k] * nsh_83[k];

        t_112[k] = pa_y[k] * msi0_56[k]
                   - f_10 * pc_y[k] * msi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, msi0_31, msh_21, \
                         msh_42, msh_44, msi1_31, nsh_84, nsh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * msh_42[k]
                   + f_3 * pc_y[k] * nsh_84[k];

        t_114[k] = f_11 * msh_21[k]
                   + f_3 * pc_z[k] * nsh_84[k];

        t_115[k] = pa_z[k] * msi0_31[k]
                   - f_10 * pc_z[k] * msi1_31[k];

        t_116[k] = f_11 * msh_44[k]
                   + f_3 * pc_y[k] * nsh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, msi0_34, msi0_61, \
                         msh_24, msh_47, msi1_34, msi1_61, nsh_87, \
                         nsh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * msi0_61[k]
                   - f_10 * pc_y[k] * msi1_61[k];

        t_118[k] = pa_z[k] * msi0_34[k]
                   - f_10 * pc_z[k] * msi1_34[k];

        t_119[k] = f_11 * msh_24[k]
                   + f_3 * pc_z[k] * nsh_87[k];

        t_120[k] = f_11 * msh_47[k]
                   + f_3 * pc_y[k] * nsh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, msi0_38, msi0_65, \
                         msh_27, msi1_38, msi1_65, nsh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * msi0_65[k]
                   - f_10 * pc_y[k] * msi1_65[k];

        t_122[k] = pa_z[k] * msi0_38[k]
                   - f_10 * pc_z[k] * msi1_38[k];

        t_123[k] = f_11 * msh_27[k]
                   + f_3 * pc_z[k] * nsh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, msi0_68, msi0_70, \
                         msh_50, msh_51, msh_99, msi1_68, msi1_70, nsh_93, \
                         nsh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * msi0_68[k]
                   + f_12 * msh_50[k]
                   - f_10 * pc_y[k] * msi1_68[k];

        t_125[k] = f_11 * msh_51[k]
                   + f_3 * pc_y[k] * nsh_93[k];

        t_126[k] = pa_y[k] * msi0_70[k]
                   - f_10 * pc_y[k] * msi1_70[k];

        t_127[k] = f_18 * msh_99[k]
                   + f_3 * pc_x[k] * nsh_99[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;

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
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_49 = buffer.data(msi0 + 49);
    const auto *msi0_83 = buffer.data(msi0 + 83);
    const auto *msi0_84 = buffer.data(msi0 + 84);
    const auto *msi0_87 = buffer.data(msi0 + 87);
    const auto *msi0_90 = buffer.data(msi0 + 90);
    const auto *msi0_94 = buffer.data(msi0 + 94);
    const auto *msi0_96 = buffer.data(msi0 + 96);
    const auto *msi0_105 = buffer.data(msi0 + 105);
    const auto *msi0_140 = buffer.data(msi0 + 140);
    const auto *msi0_143 = buffer.data(msi0 + 143);
    const auto *msi0_145 = buffer.data(msi0 + 145);
    const auto *msi0_146 = buffer.data(msi0 + 146);
    const auto *msi0_149 = buffer.data(msi0 + 149);
    const auto *msi0_150 = buffer.data(msi0 + 150);
    const auto *msi0_152 = buffer.data(msi0 + 152);
    const auto *msi0_154 = buffer.data(msi0 + 154);

    const auto *msh_36 = buffer.data(msh + 36);
    const auto *msh_42 = buffer.data(msh + 42);
    const auto *msh_59 = buffer.data(msh + 59);
    const auto *msh_60 = buffer.data(msh + 60);
    const auto *msh_61 = buffer.data(msh + 61);
    const auto *msh_62 = buffer.data(msh + 62);
    const auto *msh_63 = buffer.data(msh + 63);
    const auto *msh_66 = buffer.data(msh + 66);
    const auto *msh_68 = buffer.data(msh + 68);
    const auto *msh_69 = buffer.data(msh + 69);
    const auto *msh_70 = buffer.data(msh + 70);
    const auto *msh_72 = buffer.data(msh + 72);
    const auto *msh_78 = buffer.data(msh + 78);
    const auto *msh_83 = buffer.data(msh + 83);
    const auto *msh_84 = buffer.data(msh + 84);
    const auto *msh_86 = buffer.data(msh + 86);
    const auto *msh_87 = buffer.data(msh + 87);
    const auto *msh_89 = buffer.data(msh + 89);
    const auto *msh_90 = buffer.data(msh + 90);
    const auto *msh_93 = buffer.data(msh + 93);
    const auto *msh_99 = buffer.data(msh + 99);
    const auto *msh_100 = buffer.data(msh + 100);
    const auto *msh_101 = buffer.data(msh + 101);
    const auto *msh_102 = buffer.data(msh + 102);
    const auto *msh_103 = buffer.data(msh + 103);
    const auto *msh_104 = buffer.data(msh + 104);
    const auto *msh_105 = buffer.data(msh + 105);
    const auto *msh_106 = buffer.data(msh + 106);
    const auto *msh_107 = buffer.data(msh + 107);
    const auto *msh_108 = buffer.data(msh + 108);
    const auto *msh_110 = buffer.data(msh + 110);
    const auto *msh_111 = buffer.data(msh + 111);
    const auto *msh_113 = buffer.data(msh + 113);
    const auto *msh_114 = buffer.data(msh + 114);
    const auto *msh_119 = buffer.data(msh + 119);
    const auto *msh_120 = buffer.data(msh + 120);
    const auto *msh_121 = buffer.data(msh + 121);
    const auto *msh_122 = buffer.data(msh + 122);
    const auto *msh_123 = buffer.data(msh + 123);
    const auto *msh_125 = buffer.data(msh + 125);
    const auto *msh_126 = buffer.data(msh + 126);
    const auto *msh_129 = buffer.data(msh + 129);
    const auto *msh_132 = buffer.data(msh + 132);
    const auto *msh_136 = buffer.data(msh + 136);
    const auto *msh_141 = buffer.data(msh + 141);
    const auto *msh_143 = buffer.data(msh + 143);
    const auto *msh_144 = buffer.data(msh + 144);
    const auto *msh_145 = buffer.data(msh + 145);
    const auto *msh_146 = buffer.data(msh + 146);
    const auto *msh_152 = buffer.data(msh + 152);
    const auto *msh_156 = buffer.data(msh + 156);
    const auto *msh_161 = buffer.data(msh + 161);
    const auto *msh_162 = buffer.data(msh + 162);
    const auto *msh_163 = buffer.data(msh + 163);
    const auto *msh_164 = buffer.data(msh + 164);
    const auto *msh_165 = buffer.data(msh + 165);
    const auto *msh_166 = buffer.data(msh + 166);
    const auto *msh_167 = buffer.data(msh + 167);
    const auto *msh_183 = buffer.data(msh + 183);
    const auto *msh_184 = buffer.data(msh + 184);
    const auto *msh_185 = buffer.data(msh + 185);
    const auto *msh_186 = buffer.data(msh + 186);
    const auto *msh_187 = buffer.data(msh + 187);
    const auto *msh_188 = buffer.data(msh + 188);

    const auto *msi1_49 = buffer.data(msi1 + 49);
    const auto *msi1_83 = buffer.data(msi1 + 83);
    const auto *msi1_84 = buffer.data(msi1 + 84);
    const auto *msi1_87 = buffer.data(msi1 + 87);
    const auto *msi1_90 = buffer.data(msi1 + 90);
    const auto *msi1_94 = buffer.data(msi1 + 94);
    const auto *msi1_96 = buffer.data(msi1 + 96);
    const auto *msi1_105 = buffer.data(msi1 + 105);
    const auto *msi1_140 = buffer.data(msi1 + 140);
    const auto *msi1_143 = buffer.data(msi1 + 143);
    const auto *msi1_145 = buffer.data(msi1 + 145);
    const auto *msi1_146 = buffer.data(msi1 + 146);
    const auto *msi1_149 = buffer.data(msi1 + 149);
    const auto *msi1_150 = buffer.data(msi1 + 150);
    const auto *msi1_152 = buffer.data(msi1 + 152);
    const auto *msi1_154 = buffer.data(msi1 + 154);

    const auto *nsg0_72 = buffer.data(nsg0 + 72);
    const auto *nsg0_73 = buffer.data(nsg0 + 73);
    const auto *nsg0_74 = buffer.data(nsg0 + 74);
    const auto *nsg0_75 = buffer.data(nsg0 + 75);
    const auto *nsg0_76 = buffer.data(nsg0 + 76);
    const auto *nsg0_77 = buffer.data(nsg0 + 77);
    const auto *nsg0_78 = buffer.data(nsg0 + 78);
    const auto *nsg0_79 = buffer.data(nsg0 + 79);
    const auto *nsg0_80 = buffer.data(nsg0 + 80);
    const auto *nsg0_84 = buffer.data(nsg0 + 84);
    const auto *nsg0_85 = buffer.data(nsg0 + 85);
    const auto *nsg0_86 = buffer.data(nsg0 + 86);
    const auto *nsg0_87 = buffer.data(nsg0 + 87);
    const auto *nsg0_88 = buffer.data(nsg0 + 88);
    const auto *nsg0_89 = buffer.data(nsg0 + 89);
    const auto *nsg0_90 = buffer.data(nsg0 + 90);
    const auto *nsg0_92 = buffer.data(nsg0 + 92);
    const auto *nsg0_93 = buffer.data(nsg0 + 93);
    const auto *nsg0_95 = buffer.data(nsg0 + 95);
    const auto *nsg0_96 = buffer.data(nsg0 + 96);
    const auto *nsg0_100 = buffer.data(nsg0 + 100);
    const auto *nsg0_101 = buffer.data(nsg0 + 101);
    const auto *nsg0_102 = buffer.data(nsg0 + 102);
    const auto *nsg0_104 = buffer.data(nsg0 + 104);
    const auto *nsg0_110 = buffer.data(nsg0 + 110);
    const auto *nsg0_114 = buffer.data(nsg0 + 114);
    const auto *nsg0_117 = buffer.data(nsg0 + 117);
    const auto *nsg0_118 = buffer.data(nsg0 + 118);
    const auto *nsg0_119 = buffer.data(nsg0 + 119);
    const auto *nsg0_130 = buffer.data(nsg0 + 130);
    const auto *nsg0_132 = buffer.data(nsg0 + 132);

    const auto *nsg1_72 = buffer.data(nsg1 + 72);
    const auto *nsg1_73 = buffer.data(nsg1 + 73);
    const auto *nsg1_74 = buffer.data(nsg1 + 74);
    const auto *nsg1_75 = buffer.data(nsg1 + 75);
    const auto *nsg1_76 = buffer.data(nsg1 + 76);
    const auto *nsg1_77 = buffer.data(nsg1 + 77);
    const auto *nsg1_78 = buffer.data(nsg1 + 78);
    const auto *nsg1_79 = buffer.data(nsg1 + 79);
    const auto *nsg1_80 = buffer.data(nsg1 + 80);
    const auto *nsg1_84 = buffer.data(nsg1 + 84);
    const auto *nsg1_85 = buffer.data(nsg1 + 85);
    const auto *nsg1_86 = buffer.data(nsg1 + 86);
    const auto *nsg1_87 = buffer.data(nsg1 + 87);
    const auto *nsg1_88 = buffer.data(nsg1 + 88);
    const auto *nsg1_89 = buffer.data(nsg1 + 89);
    const auto *nsg1_90 = buffer.data(nsg1 + 90);
    const auto *nsg1_92 = buffer.data(nsg1 + 92);
    const auto *nsg1_93 = buffer.data(nsg1 + 93);
    const auto *nsg1_95 = buffer.data(nsg1 + 95);
    const auto *nsg1_96 = buffer.data(nsg1 + 96);
    const auto *nsg1_100 = buffer.data(nsg1 + 100);
    const auto *nsg1_101 = buffer.data(nsg1 + 101);
    const auto *nsg1_102 = buffer.data(nsg1 + 102);
    const auto *nsg1_104 = buffer.data(nsg1 + 104);
    const auto *nsg1_110 = buffer.data(nsg1 + 110);
    const auto *nsg1_114 = buffer.data(nsg1 + 114);
    const auto *nsg1_117 = buffer.data(nsg1 + 117);
    const auto *nsg1_118 = buffer.data(nsg1 + 118);
    const auto *nsg1_119 = buffer.data(nsg1 + 119);
    const auto *nsg1_130 = buffer.data(nsg1 + 130);
    const auto *nsg1_132 = buffer.data(nsg1 + 132);

    const auto *nsh_99 = buffer.data(nsh + 99);
    const auto *nsh_100 = buffer.data(nsh + 100);
    const auto *nsh_101 = buffer.data(nsh + 101);
    const auto *nsh_102 = buffer.data(nsh + 102);
    const auto *nsh_103 = buffer.data(nsh + 103);
    const auto *nsh_104 = buffer.data(nsh + 104);
    const auto *nsh_105 = buffer.data(nsh + 105);
    const auto *nsh_106 = buffer.data(nsh + 106);
    const auto *nsh_107 = buffer.data(nsh + 107);
    const auto *nsh_108 = buffer.data(nsh + 108);
    const auto *nsh_109 = buffer.data(nsh + 109);
    const auto *nsh_110 = buffer.data(nsh + 110);
    const auto *nsh_111 = buffer.data(nsh + 111);
    const auto *nsh_112 = buffer.data(nsh + 112);
    const auto *nsh_113 = buffer.data(nsh + 113);
    const auto *nsh_114 = buffer.data(nsh + 114);
    const auto *nsh_119 = buffer.data(nsh + 119);
    const auto *nsh_120 = buffer.data(nsh + 120);
    const auto *nsh_121 = buffer.data(nsh + 121);
    const auto *nsh_122 = buffer.data(nsh + 122);
    const auto *nsh_123 = buffer.data(nsh + 123);
    const auto *nsh_124 = buffer.data(nsh + 124);
    const auto *nsh_125 = buffer.data(nsh + 125);
    const auto *nsh_126 = buffer.data(nsh + 126);
    const auto *nsh_127 = buffer.data(nsh + 127);
    const auto *nsh_128 = buffer.data(nsh + 128);
    const auto *nsh_129 = buffer.data(nsh + 129);
    const auto *nsh_131 = buffer.data(nsh + 131);
    const auto *nsh_132 = buffer.data(nsh + 132);
    const auto *nsh_133 = buffer.data(nsh + 133);
    const auto *nsh_135 = buffer.data(nsh + 135);
    const auto *nsh_136 = buffer.data(nsh + 136);
    const auto *nsh_141 = buffer.data(nsh + 141);
    const auto *nsh_142 = buffer.data(nsh + 142);
    const auto *nsh_143 = buffer.data(nsh + 143);
    const auto *nsh_144 = buffer.data(nsh + 144);
    const auto *nsh_145 = buffer.data(nsh + 145);
    const auto *nsh_146 = buffer.data(nsh + 146);
    const auto *nsh_147 = buffer.data(nsh + 147);
    const auto *nsh_149 = buffer.data(nsh + 149);
    const auto *nsh_150 = buffer.data(nsh + 150);
    const auto *nsh_152 = buffer.data(nsh + 152);
    const auto *nsh_153 = buffer.data(nsh + 153);
    const auto *nsh_156 = buffer.data(nsh + 156);
    const auto *nsh_161 = buffer.data(nsh + 161);
    const auto *nsh_162 = buffer.data(nsh + 162);
    const auto *nsh_163 = buffer.data(nsh + 163);
    const auto *nsh_164 = buffer.data(nsh + 164);
    const auto *nsh_165 = buffer.data(nsh + 165);
    const auto *nsh_166 = buffer.data(nsh + 166);
    const auto *nsh_167 = buffer.data(nsh + 167);
    const auto *nsh_168 = buffer.data(nsh + 168);
    const auto *nsh_170 = buffer.data(nsh + 170);
    const auto *nsh_171 = buffer.data(nsh + 171);
    const auto *nsh_173 = buffer.data(nsh + 173);
    const auto *nsh_174 = buffer.data(nsh + 174);
    const auto *nsh_177 = buffer.data(nsh + 177);
    const auto *nsh_183 = buffer.data(nsh + 183);
    const auto *nsh_184 = buffer.data(nsh + 184);
    const auto *nsh_185 = buffer.data(nsh + 185);
    const auto *nsh_186 = buffer.data(nsh + 186);
    const auto *nsh_187 = buffer.data(nsh + 187);
    const auto *nsh_188 = buffer.data(nsh + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, msh_100, msh_101, msh_102, \
                         msh_103, msh_104, nsh_100, nsh_101, nsh_102, nsh_103, \
                         nsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_18 * msh_100[k]
                   + f_3 * pc_x[k] * nsh_100[k];

        t_129[k] = f_18 * msh_101[k]
                   + f_3 * pc_x[k] * nsh_101[k];

        t_130[k] = f_18 * msh_102[k]
                   + f_3 * pc_x[k] * nsh_102[k];

        t_131[k] = f_18 * msh_103[k]
                   + f_3 * pc_x[k] * nsh_103[k];

        t_132[k] = f_18 * msh_104[k]
                   + f_3 * pc_x[k] * nsh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, msi0_49, msh_36, msh_59, \
                         msi1_49, nsg0_72, nsg1_72, nsh_99, nsh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * msi0_49[k]
                   - f_10 * pc_z[k] * msi1_49[k];

        t_134[k] = f_11 * msh_36[k]
                   + f_3 * pc_z[k] * nsh_99[k];

        t_135[k] = f_11 * msh_59[k]
                   + f_8 * nsg0_72[k]
                   - f_9 * nsg1_72[k]
                   + f_3 * pc_y[k] * nsh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, msh_60, msh_61, msh_62, nsg0_73, nsg0_74, \
                         nsg1_73, nsg1_74, nsh_102, nsh_103, nsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * msh_60[k]
                   + f_6 * nsg0_73[k]
                   - f_7 * nsg1_73[k]
                   + f_3 * pc_y[k] * nsh_102[k];

        t_137[k] = f_11 * msh_61[k]
                   + f_4 * nsg0_74[k]
                   - f_5 * nsg1_74[k]
                   + f_3 * pc_y[k] * nsh_103[k];

        t_138[k] = f_11 * msh_62[k]
                   + f_3 * pc_y[k] * nsh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, msi0_83, msh_42, \
                         msh_105, msi1_83, nsg0_75, nsg1_75, nsh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * msi0_83[k]
                   - f_10 * pc_y[k] * msi1_83[k];

        t_140[k] = f_18 * msh_105[k]
                   + f_1 * nsg0_75[k]
                   - f_2 * nsg1_75[k]
                   + f_3 * pc_x[k] * nsh_105[k];

        t_141[k] = f_3 * pc_y[k] * nsh_105[k];

        t_142[k] = f_12 * msh_42[k]
                   + f_3 * pc_z[k] * nsh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, msh_110, nsg0_75, nsg0_80, nsg1_75, \
                         nsg1_80, nsh_106, nsh_107, nsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * nsg0_75[k]
                   - f_5 * nsg1_75[k]
                   + f_3 * pc_y[k] * nsh_106[k];

        t_144[k] = f_3 * pc_y[k] * nsh_107[k];

        t_145[k] = f_18 * msh_110[k]
                   + f_8 * nsg0_80[k]
                   - f_9 * nsg1_80[k]
                   + f_3 * pc_x[k] * nsh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, nsg0_76, nsg0_77, nsg1_76, nsg1_77, \
                         nsh_108, nsh_109, nsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * nsg0_76[k]
                   - f_7 * nsg1_76[k]
                   + f_3 * pc_y[k] * nsh_108[k];

        t_147[k] = f_4 * nsg0_77[k]
                   - f_5 * nsg1_77[k]
                   + f_3 * pc_y[k] * nsh_109[k];

        t_148[k] = f_3 * pc_y[k] * nsh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, msh_114, nsg0_78, nsg0_79, nsg0_84, \
                         nsg1_78, nsg1_79, nsg1_84, nsh_111, nsh_112, \
                         nsh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_18 * msh_114[k]
                   + f_6 * nsg0_84[k]
                   - f_7 * nsg1_84[k]
                   + f_3 * pc_x[k] * nsh_114[k];

        t_150[k] = f_8 * nsg0_78[k]
                   - f_9 * nsg1_78[k]
                   + f_3 * pc_y[k] * nsh_111[k];

        t_151[k] = f_6 * nsg0_79[k]
                   - f_7 * nsg1_79[k]
                   + f_3 * pc_y[k] * nsh_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, msh_119, msh_120, nsg0_80, \
                         nsg0_89, nsg1_80, nsg1_89, nsh_113, nsh_114, nsh_119, \
                         nsh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * nsg0_80[k]
                   - f_5 * nsg1_80[k]
                   + f_3 * pc_y[k] * nsh_113[k];

        t_153[k] = f_3 * pc_y[k] * nsh_114[k];

        t_154[k] = f_18 * msh_119[k]
                   + f_4 * nsg0_89[k]
                   - f_5 * nsg1_89[k]
                   + f_3 * pc_x[k] * nsh_119[k];

        t_155[k] = f_18 * msh_120[k]
                   + f_3 * pc_x[k] * nsh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, msh_121, msh_122, \
                         msh_123, msh_125, nsh_119, nsh_121, nsh_122, nsh_123, \
                         nsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_18 * msh_121[k]
                   + f_3 * pc_x[k] * nsh_121[k];

        t_157[k] = f_18 * msh_122[k]
                   + f_3 * pc_x[k] * nsh_122[k];

        t_158[k] = f_18 * msh_123[k]
                   + f_3 * pc_x[k] * nsh_123[k];

        t_159[k] = f_3 * pc_y[k] * nsh_119[k];

        t_160[k] = f_18 * msh_125[k]
                   + f_3 * pc_x[k] * nsh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, nsg0_85, nsg0_86, nsg0_87, nsg1_85, \
                         nsg1_86, nsg1_87, nsh_120, nsh_121, nsh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * nsg0_85[k]
                   - f_2 * nsg1_85[k]
                   + f_3 * pc_y[k] * nsh_120[k];

        t_162[k] = f_16 * nsg0_86[k]
                   - f_17 * nsg1_86[k]
                   + f_3 * pc_y[k] * nsh_121[k];

        t_163[k] = f_8 * nsg0_87[k]
                   - f_9 * nsg1_87[k]
                   + f_3 * pc_y[k] * nsh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, msh_62, nsg0_88, nsg0_89, \
                         nsg1_88, nsg1_89, nsh_123, nsh_124, nsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * nsg0_88[k]
                   - f_7 * nsg1_88[k]
                   + f_3 * pc_y[k] * nsh_123[k];

        t_165[k] = f_4 * nsg0_89[k]
                   - f_5 * nsg1_89[k]
                   + f_3 * pc_y[k] * nsh_124[k];

        t_166[k] = f_3 * pc_y[k] * nsh_125[k];

        t_167[k] = f_12 * msh_62[k]
                   + f_1 * nsg0_89[k]
                   - f_2 * nsg1_89[k]
                   + f_3 * pc_z[k] * nsh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, msh_63, msh_126, \
                         msh_129, nsg0_90, nsg0_93, nsg1_90, nsg1_93, nsh_126, \
                         nsh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * msh_126[k]
                   + f_1 * nsg0_90[k]
                   - f_2 * nsg1_90[k]
                   + f_3 * pc_x[k] * nsh_126[k];

        t_169[k] = f_13 * msh_63[k]
                   + f_3 * pc_y[k] * nsh_126[k];

        t_170[k] = f_3 * pc_z[k] * nsh_126[k];

        t_171[k] = f_19 * msh_129[k]
                   + f_8 * nsg0_93[k]
                   - f_9 * nsg1_93[k]
                   + f_3 * pc_x[k] * nsh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_x, pc_z, msh_132, nsg0_90, nsg0_96, \
                         nsg1_90, nsg1_96, nsh_127, nsh_128, nsh_129, \
                         nsh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * nsh_127[k];

        t_173[k] = f_4 * nsg0_90[k]
                   - f_5 * nsg1_90[k]
                   + f_3 * pc_z[k] * nsh_128[k];

        t_174[k] = f_19 * msh_132[k]
                   + f_6 * nsg0_96[k]
                   - f_7 * nsg1_96[k]
                   + f_3 * pc_x[k] * nsh_132[k];

        t_175[k] = f_3 * pc_z[k] * nsh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, msh_68, msh_136, \
                         nsg0_92, nsg0_100, nsg1_92, nsg1_100, nsh_131, nsh_132, \
                         nsh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * msh_68[k]
                   + f_3 * pc_y[k] * nsh_131[k];

        t_177[k] = f_6 * nsg0_92[k]
                   - f_7 * nsg1_92[k]
                   + f_3 * pc_z[k] * nsh_131[k];

        t_178[k] = f_19 * msh_136[k]
                   + f_4 * nsg0_100[k]
                   - f_5 * nsg1_100[k]
                   + f_3 * pc_x[k] * nsh_136[k];

        t_179[k] = f_3 * pc_z[k] * nsh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, msh_72, msh_141, \
                         nsg0_93, nsg0_95, nsg1_93, nsg1_95, nsh_133, nsh_135, \
                         nsh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * nsg0_93[k]
                   - f_5 * nsg1_93[k]
                   + f_3 * pc_z[k] * nsh_133[k];

        t_181[k] = f_13 * msh_72[k]
                   + f_3 * pc_y[k] * nsh_135[k];

        t_182[k] = f_8 * nsg0_95[k]
                   - f_9 * nsg1_95[k]
                   + f_3 * pc_z[k] * nsh_135[k];

        t_183[k] = f_19 * msh_141[k]
                   + f_3 * pc_x[k] * nsh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, msh_143, msh_144, \
                         msh_145, msh_146, nsh_136, nsh_143, nsh_144, nsh_145, \
                         nsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * nsh_136[k];

        t_185[k] = f_19 * msh_143[k]
                   + f_3 * pc_x[k] * nsh_143[k];

        t_186[k] = f_19 * msh_144[k]
                   + f_3 * pc_x[k] * nsh_144[k];

        t_187[k] = f_19 * msh_145[k]
                   + f_3 * pc_x[k] * nsh_145[k];

        t_188[k] = f_19 * msh_146[k]
                   + f_3 * pc_x[k] * nsh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_y, pc_z, msh_78, nsg0_100, nsg0_101, \
                         nsg1_100, nsg1_101, nsh_141, nsh_142, \
                         nsh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * msh_78[k]
                   + f_1 * nsg0_100[k]
                   - f_2 * nsg1_100[k]
                   + f_3 * pc_y[k] * nsh_141[k];

        t_190[k] = f_3 * pc_z[k] * nsh_141[k];

        t_191[k] = f_4 * nsg0_100[k]
                   - f_5 * nsg1_100[k]
                   + f_3 * pc_z[k] * nsh_142[k];

        t_192[k] = f_6 * nsg0_101[k]
                   - f_7 * nsg1_101[k]
                   + f_3 * pc_z[k] * nsh_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pc_y, pc_z, msi0_84, msh_83, \
                         msi1_84, nsg0_102, nsg0_104, nsg1_102, nsg1_104, nsh_144, \
                         nsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * nsg0_102[k]
                   - f_9 * nsg1_102[k]
                   + f_3 * pc_z[k] * nsh_144[k];

        t_194[k] = f_13 * msh_83[k]
                   + f_3 * pc_y[k] * nsh_146[k];

        t_195[k] = f_1 * nsg0_104[k]
                   - f_2 * nsg1_104[k]
                   + f_3 * pc_z[k] * nsh_146[k];

        t_196[k] = pa_z[k] * msi0_84[k]
                   - f_10 * pc_z[k] * msi1_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pc_y, pc_z, msi0_87, msh_63, \
                         msh_84, msh_86, msi1_87, nsh_147, nsh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * msh_84[k]
                   + f_3 * pc_y[k] * nsh_147[k];

        t_198[k] = f_11 * msh_63[k]
                   + f_3 * pc_z[k] * nsh_147[k];

        t_199[k] = pa_z[k] * msi0_87[k]
                   - f_10 * pc_z[k] * msi1_87[k];

        t_200[k] = f_12 * msh_86[k]
                   + f_3 * pc_y[k] * nsh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pc_x, pc_z, msi0_90, msh_66, msh_152, \
                         msi1_90, nsg0_110, nsg1_110, nsh_150, \
                         nsh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_19 * msh_152[k]
                   + f_8 * nsg0_110[k]
                   - f_9 * nsg1_110[k]
                   + f_3 * pc_x[k] * nsh_152[k];

        t_202[k] = pa_z[k] * msi0_90[k]
                   - f_10 * pc_z[k] * msi1_90[k];

        t_203[k] = f_11 * msh_66[k]
                   + f_3 * pc_z[k] * nsh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_y, pc_z, msi0_94, msh_89, \
                         msh_156, msi1_94, nsg0_114, nsg1_114, nsh_152, \
                         nsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * msh_89[k]
                   + f_3 * pc_y[k] * nsh_152[k];

        t_205[k] = f_19 * msh_156[k]
                   + f_6 * nsg0_114[k]
                   - f_7 * nsg1_114[k]
                   + f_3 * pc_x[k] * nsh_156[k];

        t_206[k] = pa_z[k] * msi0_94[k]
                   - f_10 * pc_z[k] * msi1_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pc_y, pc_z, msi0_96, msh_69, msh_70, \
                         msh_93, msi1_96, nsh_153, nsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_11 * msh_69[k]
                   + f_3 * pc_z[k] * nsh_153[k];

        t_208[k] = pa_z[k] * msi0_96[k]
                   + f_12 * msh_70[k]
                   - f_10 * pc_z[k] * msi1_96[k];

        t_209[k] = f_12 * msh_93[k]
                   + f_3 * pc_y[k] * nsh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, msh_161, msh_162, msh_163, msh_164, \
                         nsg0_119, nsg1_119, nsh_161, nsh_162, nsh_163, \
                         nsh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_19 * msh_161[k]
                   + f_4 * nsg0_119[k]
                   - f_5 * nsg1_119[k]
                   + f_3 * pc_x[k] * nsh_161[k];

        t_211[k] = f_19 * msh_162[k]
                   + f_3 * pc_x[k] * nsh_162[k];

        t_212[k] = f_19 * msh_163[k]
                   + f_3 * pc_x[k] * nsh_163[k];

        t_213[k] = f_19 * msh_164[k]
                   + f_3 * pc_x[k] * nsh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, msi0_105, msh_165, \
                         msh_166, msh_167, msi1_105, nsh_165, nsh_166, \
                         nsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_19 * msh_165[k]
                   + f_3 * pc_x[k] * nsh_165[k];

        t_215[k] = f_19 * msh_166[k]
                   + f_3 * pc_x[k] * nsh_166[k];

        t_216[k] = f_19 * msh_167[k]
                   + f_3 * pc_x[k] * nsh_167[k];

        t_217[k] = pa_z[k] * msi0_105[k]
                   - f_10 * pc_z[k] * msi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, msh_78, msh_101, msh_102, nsg0_117, \
                         nsg0_118, nsg1_117, nsg1_118, nsh_162, nsh_164, \
                         nsh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * msh_78[k]
                   + f_3 * pc_z[k] * nsh_162[k];

        t_219[k] = f_12 * msh_101[k]
                   + f_8 * nsg0_117[k]
                   - f_9 * nsg1_117[k]
                   + f_3 * pc_y[k] * nsh_164[k];

        t_220[k] = f_12 * msh_102[k]
                   + f_6 * nsg0_118[k]
                   - f_7 * nsg1_118[k]
                   + f_3 * pc_y[k] * nsh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pc_y, pc_z, msi0_140, msh_83, \
                         msh_103, msh_104, msi1_140, nsg0_119, nsg1_119, nsh_166, \
                         nsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * msh_103[k]
                   + f_4 * nsg0_119[k]
                   - f_5 * nsg1_119[k]
                   + f_3 * pc_y[k] * nsh_166[k];

        t_222[k] = f_12 * msh_104[k]
                   + f_3 * pc_y[k] * nsh_167[k];

        t_223[k] = f_11 * msh_83[k]
                   + f_1 * nsg0_119[k]
                   - f_2 * nsg1_119[k]
                   + f_3 * pc_z[k] * nsh_167[k];

        t_224[k] = pa_y[k] * msi0_140[k]
                   - f_10 * pc_y[k] * msi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, msi0_143, msh_84, \
                         msh_105, msh_106, msh_107, msi1_143, nsh_168, \
                         nsh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * msh_105[k]
                   + f_3 * pc_y[k] * nsh_168[k];

        t_226[k] = f_12 * msh_84[k]
                   + f_3 * pc_z[k] * nsh_168[k];

        t_227[k] = pa_y[k] * msi0_143[k]
                   + f_12 * msh_106[k]
                   - f_10 * pc_y[k] * msi1_143[k];

        t_228[k] = f_11 * msh_107[k]
                   + f_3 * pc_y[k] * nsh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_y, pc_z, msi0_145, msi0_146, \
                         msh_87, msh_108, msh_110, msi1_145, msi1_146, nsh_171, \
                         nsh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * msi0_145[k]
                   - f_10 * pc_y[k] * msi1_145[k];

        t_230[k] = pa_y[k] * msi0_146[k]
                   + f_13 * msh_108[k]
                   - f_10 * pc_y[k] * msi1_146[k];

        t_231[k] = f_12 * msh_87[k]
                   + f_3 * pc_z[k] * nsh_171[k];

        t_232[k] = f_11 * msh_110[k]
                   + f_3 * pc_y[k] * nsh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_y, pc_z, msi0_149, msi0_150, msh_90, \
                         msh_111, msi1_149, msi1_150, nsh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * msi0_149[k]
                   - f_10 * pc_y[k] * msi1_149[k];

        t_234[k] = pa_y[k] * msi0_150[k]
                   + f_14 * msh_111[k]
                   - f_10 * pc_y[k] * msi1_150[k];

        t_235[k] = f_12 * msh_90[k]
                   + f_3 * pc_z[k] * nsh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pc_x, pc_y, msi0_152, msi0_154, \
                         msh_113, msh_114, msh_183, msi1_152, msi1_154, nsh_177, \
                         nsh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * msi0_152[k]
                   + f_12 * msh_113[k]
                   - f_10 * pc_y[k] * msi1_152[k];

        t_237[k] = f_11 * msh_114[k]
                   + f_3 * pc_y[k] * nsh_177[k];

        t_238[k] = pa_y[k] * msi0_154[k]
                   - f_10 * pc_y[k] * msi1_154[k];

        t_239[k] = f_19 * msh_183[k]
                   + f_3 * pc_x[k] * nsh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, msh_184, msh_185, msh_186, \
                         msh_187, msh_188, nsh_184, nsh_185, nsh_186, nsh_187, \
                         nsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_19 * msh_184[k]
                   + f_3 * pc_x[k] * nsh_184[k];

        t_241[k] = f_19 * msh_185[k]
                   + f_3 * pc_x[k] * nsh_185[k];

        t_242[k] = f_19 * msh_186[k]
                   + f_3 * pc_x[k] * nsh_186[k];

        t_243[k] = f_19 * msh_187[k]
                   + f_3 * pc_x[k] * nsh_187[k];

        t_244[k] = f_19 * msh_188[k]
                   + f_3 * pc_x[k] * nsh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, msh_99, msh_120, msh_122, nsg0_130, \
                         nsg0_132, nsg1_130, nsg1_132, nsh_183, \
                         nsh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * msh_120[k]
                   + f_1 * nsg0_130[k]
                   - f_2 * nsg1_130[k]
                   + f_3 * pc_y[k] * nsh_183[k];

        t_246[k] = f_12 * msh_99[k]
                   + f_3 * pc_z[k] * nsh_183[k];

        t_247[k] = f_11 * msh_122[k]
                   + f_8 * nsg0_132[k]
                   - f_9 * nsg1_132[k]
                   + f_3 * pc_y[k] * nsh_185[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_167 = buffer.data(msi0 + 167);
    const auto *msi0_168 = buffer.data(msi0 + 168);
    const auto *msi0_171 = buffer.data(msi0 + 171);
    const auto *msi0_174 = buffer.data(msi0 + 174);
    const auto *msi0_178 = buffer.data(msi0 + 178);
    const auto *msi0_180 = buffer.data(msi0 + 180);
    const auto *msi0_189 = buffer.data(msi0 + 189);

    const auto *msh_105 = buffer.data(msh + 105);
    const auto *msh_123 = buffer.data(msh + 123);
    const auto *msh_124 = buffer.data(msh + 124);
    const auto *msh_125 = buffer.data(msh + 125);
    const auto *msh_126 = buffer.data(msh + 126);
    const auto *msh_129 = buffer.data(msh + 129);
    const auto *msh_131 = buffer.data(msh + 131);
    const auto *msh_132 = buffer.data(msh + 132);
    const auto *msh_133 = buffer.data(msh + 133);
    const auto *msh_135 = buffer.data(msh + 135);
    const auto *msh_141 = buffer.data(msh + 141);
    const auto *msh_146 = buffer.data(msh + 146);
    const auto *msh_147 = buffer.data(msh + 147);
    const auto *msh_149 = buffer.data(msh + 149);
    const auto *msh_150 = buffer.data(msh + 150);
    const auto *msh_152 = buffer.data(msh + 152);
    const auto *msh_153 = buffer.data(msh + 153);
    const auto *msh_156 = buffer.data(msh + 156);
    const auto *msh_162 = buffer.data(msh + 162);
    const auto *msh_164 = buffer.data(msh + 164);
    const auto *msh_165 = buffer.data(msh + 165);
    const auto *msh_166 = buffer.data(msh + 166);
    const auto *msh_167 = buffer.data(msh + 167);
    const auto *msh_168 = buffer.data(msh + 168);
    const auto *msh_170 = buffer.data(msh + 170);
    const auto *msh_173 = buffer.data(msh + 173);
    const auto *msh_177 = buffer.data(msh + 177);
    const auto *msh_183 = buffer.data(msh + 183);
    const auto *msh_185 = buffer.data(msh + 185);
    const auto *msh_186 = buffer.data(msh + 186);
    const auto *msh_187 = buffer.data(msh + 187);
    const auto *msh_189 = buffer.data(msh + 189);
    const auto *msh_194 = buffer.data(msh + 194);
    const auto *msh_198 = buffer.data(msh + 198);
    const auto *msh_203 = buffer.data(msh + 203);
    const auto *msh_204 = buffer.data(msh + 204);
    const auto *msh_205 = buffer.data(msh + 205);
    const auto *msh_206 = buffer.data(msh + 206);
    const auto *msh_207 = buffer.data(msh + 207);
    const auto *msh_209 = buffer.data(msh + 209);
    const auto *msh_210 = buffer.data(msh + 210);
    const auto *msh_213 = buffer.data(msh + 213);
    const auto *msh_216 = buffer.data(msh + 216);
    const auto *msh_220 = buffer.data(msh + 220);
    const auto *msh_225 = buffer.data(msh + 225);
    const auto *msh_227 = buffer.data(msh + 227);
    const auto *msh_228 = buffer.data(msh + 228);
    const auto *msh_229 = buffer.data(msh + 229);
    const auto *msh_230 = buffer.data(msh + 230);
    const auto *msh_236 = buffer.data(msh + 236);
    const auto *msh_240 = buffer.data(msh + 240);
    const auto *msh_245 = buffer.data(msh + 245);
    const auto *msh_246 = buffer.data(msh + 246);
    const auto *msh_247 = buffer.data(msh + 247);
    const auto *msh_248 = buffer.data(msh + 248);
    const auto *msh_249 = buffer.data(msh + 249);
    const auto *msh_250 = buffer.data(msh + 250);
    const auto *msh_251 = buffer.data(msh + 251);
    const auto *msh_252 = buffer.data(msh + 252);
    const auto *msh_255 = buffer.data(msh + 255);
    const auto *msh_257 = buffer.data(msh + 257);
    const auto *msh_258 = buffer.data(msh + 258);
    const auto *msh_261 = buffer.data(msh + 261);
    const auto *msh_262 = buffer.data(msh + 262);
    const auto *msh_264 = buffer.data(msh + 264);
    const auto *msh_266 = buffer.data(msh + 266);
    const auto *msh_267 = buffer.data(msh + 267);
    const auto *msh_268 = buffer.data(msh + 268);
    const auto *msh_269 = buffer.data(msh + 269);
    const auto *msh_270 = buffer.data(msh + 270);
    const auto *msh_271 = buffer.data(msh + 271);
    const auto *msh_272 = buffer.data(msh + 272);

    const auto *msi1_167 = buffer.data(msi1 + 167);
    const auto *msi1_168 = buffer.data(msi1 + 168);
    const auto *msi1_171 = buffer.data(msi1 + 171);
    const auto *msi1_174 = buffer.data(msi1 + 174);
    const auto *msi1_178 = buffer.data(msi1 + 178);
    const auto *msi1_180 = buffer.data(msi1 + 180);
    const auto *msi1_189 = buffer.data(msi1 + 189);

    const auto *nsg0_133 = buffer.data(nsg0 + 133);
    const auto *nsg0_134 = buffer.data(nsg0 + 134);
    const auto *nsg0_135 = buffer.data(nsg0 + 135);
    const auto *nsg0_136 = buffer.data(nsg0 + 136);
    const auto *nsg0_137 = buffer.data(nsg0 + 137);
    const auto *nsg0_138 = buffer.data(nsg0 + 138);
    const auto *nsg0_139 = buffer.data(nsg0 + 139);
    const auto *nsg0_140 = buffer.data(nsg0 + 140);
    const auto *nsg0_144 = buffer.data(nsg0 + 144);
    const auto *nsg0_145 = buffer.data(nsg0 + 145);
    const auto *nsg0_146 = buffer.data(nsg0 + 146);
    const auto *nsg0_147 = buffer.data(nsg0 + 147);
    const auto *nsg0_148 = buffer.data(nsg0 + 148);
    const auto *nsg0_149 = buffer.data(nsg0 + 149);
    const auto *nsg0_150 = buffer.data(nsg0 + 150);
    const auto *nsg0_152 = buffer.data(nsg0 + 152);
    const auto *nsg0_153 = buffer.data(nsg0 + 153);
    const auto *nsg0_155 = buffer.data(nsg0 + 155);
    const auto *nsg0_156 = buffer.data(nsg0 + 156);
    const auto *nsg0_160 = buffer.data(nsg0 + 160);
    const auto *nsg0_161 = buffer.data(nsg0 + 161);
    const auto *nsg0_162 = buffer.data(nsg0 + 162);
    const auto *nsg0_164 = buffer.data(nsg0 + 164);
    const auto *nsg0_170 = buffer.data(nsg0 + 170);
    const auto *nsg0_174 = buffer.data(nsg0 + 174);
    const auto *nsg0_177 = buffer.data(nsg0 + 177);
    const auto *nsg0_178 = buffer.data(nsg0 + 178);
    const auto *nsg0_179 = buffer.data(nsg0 + 179);
    const auto *nsg0_180 = buffer.data(nsg0 + 180);
    const auto *nsg0_183 = buffer.data(nsg0 + 183);
    const auto *nsg0_185 = buffer.data(nsg0 + 185);
    const auto *nsg0_186 = buffer.data(nsg0 + 186);
    const auto *nsg0_189 = buffer.data(nsg0 + 189);
    const auto *nsg0_190 = buffer.data(nsg0 + 190);
    const auto *nsg0_192 = buffer.data(nsg0 + 192);
    const auto *nsg0_193 = buffer.data(nsg0 + 193);
    const auto *nsg0_194 = buffer.data(nsg0 + 194);

    const auto *nsg1_133 = buffer.data(nsg1 + 133);
    const auto *nsg1_134 = buffer.data(nsg1 + 134);
    const auto *nsg1_135 = buffer.data(nsg1 + 135);
    const auto *nsg1_136 = buffer.data(nsg1 + 136);
    const auto *nsg1_137 = buffer.data(nsg1 + 137);
    const auto *nsg1_138 = buffer.data(nsg1 + 138);
    const auto *nsg1_139 = buffer.data(nsg1 + 139);
    const auto *nsg1_140 = buffer.data(nsg1 + 140);
    const auto *nsg1_144 = buffer.data(nsg1 + 144);
    const auto *nsg1_145 = buffer.data(nsg1 + 145);
    const auto *nsg1_146 = buffer.data(nsg1 + 146);
    const auto *nsg1_147 = buffer.data(nsg1 + 147);
    const auto *nsg1_148 = buffer.data(nsg1 + 148);
    const auto *nsg1_149 = buffer.data(nsg1 + 149);
    const auto *nsg1_150 = buffer.data(nsg1 + 150);
    const auto *nsg1_152 = buffer.data(nsg1 + 152);
    const auto *nsg1_153 = buffer.data(nsg1 + 153);
    const auto *nsg1_155 = buffer.data(nsg1 + 155);
    const auto *nsg1_156 = buffer.data(nsg1 + 156);
    const auto *nsg1_160 = buffer.data(nsg1 + 160);
    const auto *nsg1_161 = buffer.data(nsg1 + 161);
    const auto *nsg1_162 = buffer.data(nsg1 + 162);
    const auto *nsg1_164 = buffer.data(nsg1 + 164);
    const auto *nsg1_170 = buffer.data(nsg1 + 170);
    const auto *nsg1_174 = buffer.data(nsg1 + 174);
    const auto *nsg1_177 = buffer.data(nsg1 + 177);
    const auto *nsg1_178 = buffer.data(nsg1 + 178);
    const auto *nsg1_179 = buffer.data(nsg1 + 179);
    const auto *nsg1_180 = buffer.data(nsg1 + 180);
    const auto *nsg1_183 = buffer.data(nsg1 + 183);
    const auto *nsg1_185 = buffer.data(nsg1 + 185);
    const auto *nsg1_186 = buffer.data(nsg1 + 186);
    const auto *nsg1_189 = buffer.data(nsg1 + 189);
    const auto *nsg1_190 = buffer.data(nsg1 + 190);
    const auto *nsg1_192 = buffer.data(nsg1 + 192);
    const auto *nsg1_193 = buffer.data(nsg1 + 193);
    const auto *nsg1_194 = buffer.data(nsg1 + 194);

    const auto *nsh_186 = buffer.data(nsh + 186);
    const auto *nsh_187 = buffer.data(nsh + 187);
    const auto *nsh_188 = buffer.data(nsh + 188);
    const auto *nsh_189 = buffer.data(nsh + 189);
    const auto *nsh_190 = buffer.data(nsh + 190);
    const auto *nsh_191 = buffer.data(nsh + 191);
    const auto *nsh_192 = buffer.data(nsh + 192);
    const auto *nsh_193 = buffer.data(nsh + 193);
    const auto *nsh_194 = buffer.data(nsh + 194);
    const auto *nsh_195 = buffer.data(nsh + 195);
    const auto *nsh_196 = buffer.data(nsh + 196);
    const auto *nsh_197 = buffer.data(nsh + 197);
    const auto *nsh_198 = buffer.data(nsh + 198);
    const auto *nsh_203 = buffer.data(nsh + 203);
    const auto *nsh_204 = buffer.data(nsh + 204);
    const auto *nsh_205 = buffer.data(nsh + 205);
    const auto *nsh_206 = buffer.data(nsh + 206);
    const auto *nsh_207 = buffer.data(nsh + 207);
    const auto *nsh_208 = buffer.data(nsh + 208);
    const auto *nsh_209 = buffer.data(nsh + 209);
    const auto *nsh_210 = buffer.data(nsh + 210);
    const auto *nsh_211 = buffer.data(nsh + 211);
    const auto *nsh_212 = buffer.data(nsh + 212);
    const auto *nsh_213 = buffer.data(nsh + 213);
    const auto *nsh_215 = buffer.data(nsh + 215);
    const auto *nsh_216 = buffer.data(nsh + 216);
    const auto *nsh_217 = buffer.data(nsh + 217);
    const auto *nsh_219 = buffer.data(nsh + 219);
    const auto *nsh_220 = buffer.data(nsh + 220);
    const auto *nsh_225 = buffer.data(nsh + 225);
    const auto *nsh_226 = buffer.data(nsh + 226);
    const auto *nsh_227 = buffer.data(nsh + 227);
    const auto *nsh_228 = buffer.data(nsh + 228);
    const auto *nsh_229 = buffer.data(nsh + 229);
    const auto *nsh_230 = buffer.data(nsh + 230);
    const auto *nsh_231 = buffer.data(nsh + 231);
    const auto *nsh_233 = buffer.data(nsh + 233);
    const auto *nsh_234 = buffer.data(nsh + 234);
    const auto *nsh_236 = buffer.data(nsh + 236);
    const auto *nsh_237 = buffer.data(nsh + 237);
    const auto *nsh_240 = buffer.data(nsh + 240);
    const auto *nsh_245 = buffer.data(nsh + 245);
    const auto *nsh_246 = buffer.data(nsh + 246);
    const auto *nsh_247 = buffer.data(nsh + 247);
    const auto *nsh_248 = buffer.data(nsh + 248);
    const auto *nsh_249 = buffer.data(nsh + 249);
    const auto *nsh_250 = buffer.data(nsh + 250);
    const auto *nsh_251 = buffer.data(nsh + 251);
    const auto *nsh_252 = buffer.data(nsh + 252);
    const auto *nsh_254 = buffer.data(nsh + 254);
    const auto *nsh_255 = buffer.data(nsh + 255);
    const auto *nsh_257 = buffer.data(nsh + 257);
    const auto *nsh_258 = buffer.data(nsh + 258);
    const auto *nsh_261 = buffer.data(nsh + 261);
    const auto *nsh_262 = buffer.data(nsh + 262);
    const auto *nsh_264 = buffer.data(nsh + 264);
    const auto *nsh_266 = buffer.data(nsh + 266);
    const auto *nsh_267 = buffer.data(nsh + 267);
    const auto *nsh_268 = buffer.data(nsh + 268);
    const auto *nsh_269 = buffer.data(nsh + 269);
    const auto *nsh_270 = buffer.data(nsh + 270);
    const auto *nsh_271 = buffer.data(nsh + 271);
    const auto *nsh_272 = buffer.data(nsh + 272);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, msh_123, msh_124, msh_125, nsg0_133, \
                         nsg0_134, nsg1_133, nsg1_134, nsh_186, nsh_187, \
                         nsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * msh_123[k]
                   + f_6 * nsg0_133[k]
                   - f_7 * nsg1_133[k]
                   + f_3 * pc_y[k] * nsh_186[k];

        t_249[k] = f_11 * msh_124[k]
                   + f_4 * nsg0_134[k]
                   - f_5 * nsg1_134[k]
                   + f_3 * pc_y[k] * nsh_187[k];

        t_250[k] = f_11 * msh_125[k]
                   + f_3 * pc_y[k] * nsh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pc_x, pc_y, pc_z, msi0_167, \
                         msh_105, msh_189, msi1_167, nsg0_135, nsg1_135, \
                         nsh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * msi0_167[k]
                   - f_10 * pc_y[k] * msi1_167[k];

        t_252[k] = f_19 * msh_189[k]
                   + f_1 * nsg0_135[k]
                   - f_2 * nsg1_135[k]
                   + f_3 * pc_x[k] * nsh_189[k];

        t_253[k] = f_3 * pc_y[k] * nsh_189[k];

        t_254[k] = f_13 * msh_105[k]
                   + f_3 * pc_z[k] * nsh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, msh_194, nsg0_135, nsg0_140, \
                         nsg1_135, nsg1_140, nsh_190, nsh_191, \
                         nsh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * nsg0_135[k]
                   - f_5 * nsg1_135[k]
                   + f_3 * pc_y[k] * nsh_190[k];

        t_256[k] = f_3 * pc_y[k] * nsh_191[k];

        t_257[k] = f_19 * msh_194[k]
                   + f_8 * nsg0_140[k]
                   - f_9 * nsg1_140[k]
                   + f_3 * pc_x[k] * nsh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_y, nsg0_136, nsg0_137, nsg1_136, nsg1_137, \
                         nsh_192, nsh_193, nsh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_6 * nsg0_136[k]
                   - f_7 * nsg1_136[k]
                   + f_3 * pc_y[k] * nsh_192[k];

        t_259[k] = f_4 * nsg0_137[k]
                   - f_5 * nsg1_137[k]
                   + f_3 * pc_y[k] * nsh_193[k];

        t_260[k] = f_3 * pc_y[k] * nsh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_y, msh_198, nsg0_138, nsg0_139, \
                         nsg0_144, nsg1_138, nsg1_139, nsg1_144, nsh_195, nsh_196, \
                         nsh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_19 * msh_198[k]
                   + f_6 * nsg0_144[k]
                   - f_7 * nsg1_144[k]
                   + f_3 * pc_x[k] * nsh_198[k];

        t_262[k] = f_8 * nsg0_138[k]
                   - f_9 * nsg1_138[k]
                   + f_3 * pc_y[k] * nsh_195[k];

        t_263[k] = f_6 * nsg0_139[k]
                   - f_7 * nsg1_139[k]
                   + f_3 * pc_y[k] * nsh_196[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, msh_203, msh_204, nsg0_140, \
                         nsg0_149, nsg1_140, nsg1_149, nsh_197, nsh_198, nsh_203, \
                         nsh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * nsg0_140[k]
                   - f_5 * nsg1_140[k]
                   + f_3 * pc_y[k] * nsh_197[k];

        t_265[k] = f_3 * pc_y[k] * nsh_198[k];

        t_266[k] = f_19 * msh_203[k]
                   + f_4 * nsg0_149[k]
                   - f_5 * nsg1_149[k]
                   + f_3 * pc_x[k] * nsh_203[k];

        t_267[k] = f_19 * msh_204[k]
                   + f_3 * pc_x[k] * nsh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pc_x, pc_y, msh_205, msh_206, \
                         msh_207, msh_209, nsh_203, nsh_205, nsh_206, nsh_207, \
                         nsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_19 * msh_205[k]
                   + f_3 * pc_x[k] * nsh_205[k];

        t_269[k] = f_19 * msh_206[k]
                   + f_3 * pc_x[k] * nsh_206[k];

        t_270[k] = f_19 * msh_207[k]
                   + f_3 * pc_x[k] * nsh_207[k];

        t_271[k] = f_3 * pc_y[k] * nsh_203[k];

        t_272[k] = f_19 * msh_209[k]
                   + f_3 * pc_x[k] * nsh_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_y, nsg0_145, nsg0_146, nsg0_147, nsg1_145, \
                         nsg1_146, nsg1_147, nsh_204, nsh_205, \
                         nsh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * nsg0_145[k]
                   - f_2 * nsg1_145[k]
                   + f_3 * pc_y[k] * nsh_204[k];

        t_274[k] = f_16 * nsg0_146[k]
                   - f_17 * nsg1_146[k]
                   + f_3 * pc_y[k] * nsh_205[k];

        t_275[k] = f_8 * nsg0_147[k]
                   - f_9 * nsg1_147[k]
                   + f_3 * pc_y[k] * nsh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, msh_125, nsg0_148, nsg0_149, \
                         nsg1_148, nsg1_149, nsh_207, nsh_208, \
                         nsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * nsg0_148[k]
                   - f_7 * nsg1_148[k]
                   + f_3 * pc_y[k] * nsh_207[k];

        t_277[k] = f_4 * nsg0_149[k]
                   - f_5 * nsg1_149[k]
                   + f_3 * pc_y[k] * nsh_208[k];

        t_278[k] = f_3 * pc_y[k] * nsh_209[k];

        t_279[k] = f_13 * msh_125[k]
                   + f_1 * nsg0_149[k]
                   - f_2 * nsg1_149[k]
                   + f_3 * pc_z[k] * nsh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, msh_126, msh_210, \
                         msh_213, nsg0_150, nsg0_153, nsg1_150, nsg1_153, nsh_210, \
                         nsh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_20 * msh_210[k]
                   + f_1 * nsg0_150[k]
                   - f_2 * nsg1_150[k]
                   + f_3 * pc_x[k] * nsh_210[k];

        t_281[k] = f_14 * msh_126[k]
                   + f_3 * pc_y[k] * nsh_210[k];

        t_282[k] = f_3 * pc_z[k] * nsh_210[k];

        t_283[k] = f_20 * msh_213[k]
                   + f_8 * nsg0_153[k]
                   - f_9 * nsg1_153[k]
                   + f_3 * pc_x[k] * nsh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_z, msh_216, nsg0_150, nsg0_156, \
                         nsg1_150, nsg1_156, nsh_211, nsh_212, nsh_213, \
                         nsh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * nsh_211[k];

        t_285[k] = f_4 * nsg0_150[k]
                   - f_5 * nsg1_150[k]
                   + f_3 * pc_z[k] * nsh_212[k];

        t_286[k] = f_20 * msh_216[k]
                   + f_6 * nsg0_156[k]
                   - f_7 * nsg1_156[k]
                   + f_3 * pc_x[k] * nsh_216[k];

        t_287[k] = f_3 * pc_z[k] * nsh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, pc_z, msh_131, msh_220, \
                         nsg0_152, nsg0_160, nsg1_152, nsg1_160, nsh_215, nsh_216, \
                         nsh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_14 * msh_131[k]
                   + f_3 * pc_y[k] * nsh_215[k];

        t_289[k] = f_6 * nsg0_152[k]
                   - f_7 * nsg1_152[k]
                   + f_3 * pc_z[k] * nsh_215[k];

        t_290[k] = f_20 * msh_220[k]
                   + f_4 * nsg0_160[k]
                   - f_5 * nsg1_160[k]
                   + f_3 * pc_x[k] * nsh_220[k];

        t_291[k] = f_3 * pc_z[k] * nsh_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, msh_135, msh_225, \
                         nsg0_153, nsg0_155, nsg1_153, nsg1_155, nsh_217, nsh_219, \
                         nsh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * nsg0_153[k]
                   - f_5 * nsg1_153[k]
                   + f_3 * pc_z[k] * nsh_217[k];

        t_293[k] = f_14 * msh_135[k]
                   + f_3 * pc_y[k] * nsh_219[k];

        t_294[k] = f_8 * nsg0_155[k]
                   - f_9 * nsg1_155[k]
                   + f_3 * pc_z[k] * nsh_219[k];

        t_295[k] = f_20 * msh_225[k]
                   + f_3 * pc_x[k] * nsh_225[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, pc_z, msh_227, msh_228, \
                         msh_229, msh_230, nsh_220, nsh_227, nsh_228, nsh_229, \
                         nsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_3 * pc_z[k] * nsh_220[k];

        t_297[k] = f_20 * msh_227[k]
                   + f_3 * pc_x[k] * nsh_227[k];

        t_298[k] = f_20 * msh_228[k]
                   + f_3 * pc_x[k] * nsh_228[k];

        t_299[k] = f_20 * msh_229[k]
                   + f_3 * pc_x[k] * nsh_229[k];

        t_300[k] = f_20 * msh_230[k]
                   + f_3 * pc_x[k] * nsh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_y, pc_z, msh_141, nsg0_160, nsg0_161, \
                         nsg1_160, nsg1_161, nsh_225, nsh_226, \
                         nsh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * msh_141[k]
                   + f_1 * nsg0_160[k]
                   - f_2 * nsg1_160[k]
                   + f_3 * pc_y[k] * nsh_225[k];

        t_302[k] = f_3 * pc_z[k] * nsh_225[k];

        t_303[k] = f_4 * nsg0_160[k]
                   - f_5 * nsg1_160[k]
                   + f_3 * pc_z[k] * nsh_226[k];

        t_304[k] = f_6 * nsg0_161[k]
                   - f_7 * nsg1_161[k]
                   + f_3 * pc_z[k] * nsh_227[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pc_y, pc_z, msi0_168, msh_146, \
                         msi1_168, nsg0_162, nsg0_164, nsg1_162, nsg1_164, nsh_228, \
                         nsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_8 * nsg0_162[k]
                   - f_9 * nsg1_162[k]
                   + f_3 * pc_z[k] * nsh_228[k];

        t_306[k] = f_14 * msh_146[k]
                   + f_3 * pc_y[k] * nsh_230[k];

        t_307[k] = f_1 * nsg0_164[k]
                   - f_2 * nsg1_164[k]
                   + f_3 * pc_z[k] * nsh_230[k];

        t_308[k] = pa_z[k] * msi0_168[k]
                   - f_10 * pc_z[k] * msi1_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_z, pc_y, pc_z, msi0_171, msh_126, \
                         msh_147, msh_149, msi1_171, nsh_231, nsh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_13 * msh_147[k]
                   + f_3 * pc_y[k] * nsh_231[k];

        t_310[k] = f_11 * msh_126[k]
                   + f_3 * pc_z[k] * nsh_231[k];

        t_311[k] = pa_z[k] * msi0_171[k]
                   - f_10 * pc_z[k] * msi1_171[k];

        t_312[k] = f_13 * msh_149[k]
                   + f_3 * pc_y[k] * nsh_233[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_z, pc_x, pc_z, msi0_174, msh_129, msh_236, \
                         msi1_174, nsg0_170, nsg1_170, nsh_234, \
                         nsh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_20 * msh_236[k]
                   + f_8 * nsg0_170[k]
                   - f_9 * nsg1_170[k]
                   + f_3 * pc_x[k] * nsh_236[k];

        t_314[k] = pa_z[k] * msi0_174[k]
                   - f_10 * pc_z[k] * msi1_174[k];

        t_315[k] = f_11 * msh_129[k]
                   + f_3 * pc_z[k] * nsh_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_z, pc_x, pc_y, pc_z, msi0_178, msh_152, \
                         msh_240, msi1_178, nsg0_174, nsg1_174, nsh_236, \
                         nsh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * msh_152[k]
                   + f_3 * pc_y[k] * nsh_236[k];

        t_317[k] = f_20 * msh_240[k]
                   + f_6 * nsg0_174[k]
                   - f_7 * nsg1_174[k]
                   + f_3 * pc_x[k] * nsh_240[k];

        t_318[k] = pa_z[k] * msi0_178[k]
                   - f_10 * pc_z[k] * msi1_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_y, pc_z, msi0_180, msh_132, msh_133, \
                         msh_156, msi1_180, nsh_237, nsh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_11 * msh_132[k]
                   + f_3 * pc_z[k] * nsh_237[k];

        t_320[k] = pa_z[k] * msi0_180[k]
                   + f_12 * msh_133[k]
                   - f_10 * pc_z[k] * msi1_180[k];

        t_321[k] = f_13 * msh_156[k]
                   + f_3 * pc_y[k] * nsh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, msh_245, msh_246, msh_247, msh_248, \
                         nsg0_179, nsg1_179, nsh_245, nsh_246, nsh_247, \
                         nsh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_20 * msh_245[k]
                   + f_4 * nsg0_179[k]
                   - f_5 * nsg1_179[k]
                   + f_3 * pc_x[k] * nsh_245[k];

        t_323[k] = f_20 * msh_246[k]
                   + f_3 * pc_x[k] * nsh_246[k];

        t_324[k] = f_20 * msh_247[k]
                   + f_3 * pc_x[k] * nsh_247[k];

        t_325[k] = f_20 * msh_248[k]
                   + f_3 * pc_x[k] * nsh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pa_z, pc_x, pc_z, msi0_189, msh_249, \
                         msh_250, msh_251, msi1_189, nsh_249, nsh_250, \
                         nsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_20 * msh_249[k]
                   + f_3 * pc_x[k] * nsh_249[k];

        t_327[k] = f_20 * msh_250[k]
                   + f_3 * pc_x[k] * nsh_250[k];

        t_328[k] = f_20 * msh_251[k]
                   + f_3 * pc_x[k] * nsh_251[k];

        t_329[k] = pa_z[k] * msi0_189[k]
                   - f_10 * pc_z[k] * msi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, msh_141, msh_164, msh_165, nsg0_177, \
                         nsg0_178, nsg1_177, nsg1_178, nsh_246, nsh_248, \
                         nsh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * msh_141[k]
                   + f_3 * pc_z[k] * nsh_246[k];

        t_331[k] = f_13 * msh_164[k]
                   + f_8 * nsg0_177[k]
                   - f_9 * nsg1_177[k]
                   + f_3 * pc_y[k] * nsh_248[k];

        t_332[k] = f_13 * msh_165[k]
                   + f_6 * nsg0_178[k]
                   - f_7 * nsg1_178[k]
                   + f_3 * pc_y[k] * nsh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, msh_146, msh_166, msh_167, nsg0_179, \
                         nsg1_179, nsh_250, nsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * msh_166[k]
                   + f_4 * nsg0_179[k]
                   - f_5 * nsg1_179[k]
                   + f_3 * pc_y[k] * nsh_250[k];

        t_334[k] = f_13 * msh_167[k]
                   + f_3 * pc_y[k] * nsh_251[k];

        t_335[k] = f_11 * msh_146[k]
                   + f_1 * nsg0_179[k]
                   - f_2 * nsg1_179[k]
                   + f_3 * pc_z[k] * nsh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, msh_147, msh_168, msh_252, \
                         nsg0_180, nsg1_180, nsh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_20 * msh_252[k]
                   + f_1 * nsg0_180[k]
                   - f_2 * nsg1_180[k]
                   + f_3 * pc_x[k] * nsh_252[k];

        t_337[k] = f_12 * msh_168[k]
                   + f_3 * pc_y[k] * nsh_252[k];

        t_338[k] = f_12 * msh_147[k]
                   + f_3 * pc_z[k] * nsh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, msh_170, msh_255, msh_257, nsg0_183, \
                         nsg0_185, nsg1_183, nsg1_185, nsh_254, nsh_255, \
                         nsh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_20 * msh_255[k]
                   + f_8 * nsg0_183[k]
                   - f_9 * nsg1_183[k]
                   + f_3 * pc_x[k] * nsh_255[k];

        t_340[k] = f_12 * msh_170[k]
                   + f_3 * pc_y[k] * nsh_254[k];

        t_341[k] = f_20 * msh_257[k]
                   + f_8 * nsg0_185[k]
                   - f_9 * nsg1_185[k]
                   + f_3 * pc_x[k] * nsh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, msh_150, msh_173, msh_258, \
                         nsg0_186, nsg1_186, nsh_255, nsh_257, \
                         nsh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_20 * msh_258[k]
                   + f_6 * nsg0_186[k]
                   - f_7 * nsg1_186[k]
                   + f_3 * pc_x[k] * nsh_258[k];

        t_343[k] = f_12 * msh_150[k]
                   + f_3 * pc_z[k] * nsh_255[k];

        t_344[k] = f_12 * msh_173[k]
                   + f_3 * pc_y[k] * nsh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, msh_153, msh_261, msh_262, nsg0_189, \
                         nsg0_190, nsg1_189, nsg1_190, nsh_258, nsh_261, \
                         nsh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_20 * msh_261[k]
                   + f_6 * nsg0_189[k]
                   - f_7 * nsg1_189[k]
                   + f_3 * pc_x[k] * nsh_261[k];

        t_346[k] = f_20 * msh_262[k]
                   + f_4 * nsg0_190[k]
                   - f_5 * nsg1_190[k]
                   + f_3 * pc_x[k] * nsh_262[k];

        t_347[k] = f_12 * msh_153[k]
                   + f_3 * pc_z[k] * nsh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, msh_177, msh_264, msh_266, nsg0_192, \
                         nsg0_194, nsg1_192, nsg1_194, nsh_261, nsh_264, \
                         nsh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_20 * msh_264[k]
                   + f_4 * nsg0_192[k]
                   - f_5 * nsg1_192[k]
                   + f_3 * pc_x[k] * nsh_264[k];

        t_349[k] = f_12 * msh_177[k]
                   + f_3 * pc_y[k] * nsh_261[k];

        t_350[k] = f_20 * msh_266[k]
                   + f_4 * nsg0_194[k]
                   - f_5 * nsg1_194[k]
                   + f_3 * pc_x[k] * nsh_266[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, msh_267, msh_268, msh_269, \
                         msh_270, msh_271, nsh_267, nsh_268, nsh_269, nsh_270, \
                         nsh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_20 * msh_267[k]
                   + f_3 * pc_x[k] * nsh_267[k];

        t_352[k] = f_20 * msh_268[k]
                   + f_3 * pc_x[k] * nsh_268[k];

        t_353[k] = f_20 * msh_269[k]
                   + f_3 * pc_x[k] * nsh_269[k];

        t_354[k] = f_20 * msh_270[k]
                   + f_3 * pc_x[k] * nsh_270[k];

        t_355[k] = f_20 * msh_271[k]
                   + f_3 * pc_x[k] * nsh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, msh_162, msh_183, msh_272, \
                         nsg0_190, nsg1_190, nsh_267, nsh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_20 * msh_272[k]
                   + f_3 * pc_x[k] * nsh_272[k];

        t_357[k] = f_12 * msh_183[k]
                   + f_1 * nsg0_190[k]
                   - f_2 * nsg1_190[k]
                   + f_3 * pc_y[k] * nsh_267[k];

        t_358[k] = f_12 * msh_162[k]
                   + f_3 * pc_z[k] * nsh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, msh_185, msh_186, msh_187, nsg0_192, \
                         nsg0_193, nsg0_194, nsg1_192, nsg1_193, nsg1_194, nsh_269, nsh_270, \
                         nsh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * msh_185[k]
                   + f_8 * nsg0_192[k]
                   - f_9 * nsg1_192[k]
                   + f_3 * pc_y[k] * nsh_269[k];

        t_360[k] = f_12 * msh_186[k]
                   + f_6 * nsg0_193[k]
                   - f_7 * nsg1_193[k]
                   + f_3 * pc_y[k] * nsh_270[k];

        t_361[k] = f_12 * msh_187[k]
                   + f_4 * nsg0_194[k]
                   - f_5 * nsg1_194[k]
                   + f_3 * pc_y[k] * nsh_271[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_252 = buffer.data(msi0 + 252);
    const auto *msi0_255 = buffer.data(msi0 + 255);
    const auto *msi0_257 = buffer.data(msi0 + 257);
    const auto *msi0_258 = buffer.data(msi0 + 258);
    const auto *msi0_261 = buffer.data(msi0 + 261);
    const auto *msi0_262 = buffer.data(msi0 + 262);
    const auto *msi0_264 = buffer.data(msi0 + 264);
    const auto *msi0_266 = buffer.data(msi0 + 266);
    const auto *msi0_279 = buffer.data(msi0 + 279);
    const auto *msi0_280 = buffer.data(msi0 + 280);
    const auto *msi0_283 = buffer.data(msi0 + 283);
    const auto *msi0_286 = buffer.data(msi0 + 286);
    const auto *msi0_290 = buffer.data(msi0 + 290);
    const auto *msi0_292 = buffer.data(msi0 + 292);
    const auto *msi0_301 = buffer.data(msi0 + 301);

    const auto *msh_167 = buffer.data(msh + 167);
    const auto *msh_168 = buffer.data(msh + 168);
    const auto *msh_171 = buffer.data(msh + 171);
    const auto *msh_174 = buffer.data(msh + 174);
    const auto *msh_183 = buffer.data(msh + 183);
    const auto *msh_188 = buffer.data(msh + 188);
    const auto *msh_189 = buffer.data(msh + 189);
    const auto *msh_190 = buffer.data(msh + 190);
    const auto *msh_191 = buffer.data(msh + 191);
    const auto *msh_192 = buffer.data(msh + 192);
    const auto *msh_194 = buffer.data(msh + 194);
    const auto *msh_195 = buffer.data(msh + 195);
    const auto *msh_197 = buffer.data(msh + 197);
    const auto *msh_198 = buffer.data(msh + 198);
    const auto *msh_204 = buffer.data(msh + 204);
    const auto *msh_206 = buffer.data(msh + 206);
    const auto *msh_207 = buffer.data(msh + 207);
    const auto *msh_208 = buffer.data(msh + 208);
    const auto *msh_209 = buffer.data(msh + 209);
    const auto *msh_210 = buffer.data(msh + 210);
    const auto *msh_213 = buffer.data(msh + 213);
    const auto *msh_215 = buffer.data(msh + 215);
    const auto *msh_216 = buffer.data(msh + 216);
    const auto *msh_217 = buffer.data(msh + 217);
    const auto *msh_219 = buffer.data(msh + 219);
    const auto *msh_225 = buffer.data(msh + 225);
    const auto *msh_230 = buffer.data(msh + 230);
    const auto *msh_231 = buffer.data(msh + 231);
    const auto *msh_233 = buffer.data(msh + 233);
    const auto *msh_236 = buffer.data(msh + 236);
    const auto *msh_240 = buffer.data(msh + 240);
    const auto *msh_248 = buffer.data(msh + 248);
    const auto *msh_249 = buffer.data(msh + 249);
    const auto *msh_250 = buffer.data(msh + 250);
    const auto *msh_251 = buffer.data(msh + 251);
    const auto *msh_252 = buffer.data(msh + 252);
    const auto *msh_288 = buffer.data(msh + 288);
    const auto *msh_289 = buffer.data(msh + 289);
    const auto *msh_290 = buffer.data(msh + 290);
    const auto *msh_291 = buffer.data(msh + 291);
    const auto *msh_292 = buffer.data(msh + 292);
    const auto *msh_293 = buffer.data(msh + 293);
    const auto *msh_294 = buffer.data(msh + 294);
    const auto *msh_299 = buffer.data(msh + 299);
    const auto *msh_303 = buffer.data(msh + 303);
    const auto *msh_308 = buffer.data(msh + 308);
    const auto *msh_309 = buffer.data(msh + 309);
    const auto *msh_310 = buffer.data(msh + 310);
    const auto *msh_311 = buffer.data(msh + 311);
    const auto *msh_312 = buffer.data(msh + 312);
    const auto *msh_314 = buffer.data(msh + 314);
    const auto *msh_315 = buffer.data(msh + 315);
    const auto *msh_318 = buffer.data(msh + 318);
    const auto *msh_321 = buffer.data(msh + 321);
    const auto *msh_325 = buffer.data(msh + 325);
    const auto *msh_330 = buffer.data(msh + 330);
    const auto *msh_332 = buffer.data(msh + 332);
    const auto *msh_333 = buffer.data(msh + 333);
    const auto *msh_334 = buffer.data(msh + 334);
    const auto *msh_335 = buffer.data(msh + 335);
    const auto *msh_341 = buffer.data(msh + 341);
    const auto *msh_345 = buffer.data(msh + 345);
    const auto *msh_350 = buffer.data(msh + 350);
    const auto *msh_351 = buffer.data(msh + 351);
    const auto *msh_352 = buffer.data(msh + 352);
    const auto *msh_353 = buffer.data(msh + 353);
    const auto *msh_354 = buffer.data(msh + 354);
    const auto *msh_355 = buffer.data(msh + 355);
    const auto *msh_356 = buffer.data(msh + 356);
    const auto *msh_357 = buffer.data(msh + 357);

    const auto *msi1_252 = buffer.data(msi1 + 252);
    const auto *msi1_255 = buffer.data(msi1 + 255);
    const auto *msi1_257 = buffer.data(msi1 + 257);
    const auto *msi1_258 = buffer.data(msi1 + 258);
    const auto *msi1_261 = buffer.data(msi1 + 261);
    const auto *msi1_262 = buffer.data(msi1 + 262);
    const auto *msi1_264 = buffer.data(msi1 + 264);
    const auto *msi1_266 = buffer.data(msi1 + 266);
    const auto *msi1_279 = buffer.data(msi1 + 279);
    const auto *msi1_280 = buffer.data(msi1 + 280);
    const auto *msi1_283 = buffer.data(msi1 + 283);
    const auto *msi1_286 = buffer.data(msi1 + 286);
    const auto *msi1_290 = buffer.data(msi1 + 290);
    const auto *msi1_292 = buffer.data(msi1 + 292);
    const auto *msi1_301 = buffer.data(msi1 + 301);

    const auto *nsg0_194 = buffer.data(nsg0 + 194);
    const auto *nsg0_205 = buffer.data(nsg0 + 205);
    const auto *nsg0_207 = buffer.data(nsg0 + 207);
    const auto *nsg0_208 = buffer.data(nsg0 + 208);
    const auto *nsg0_209 = buffer.data(nsg0 + 209);
    const auto *nsg0_210 = buffer.data(nsg0 + 210);
    const auto *nsg0_211 = buffer.data(nsg0 + 211);
    const auto *nsg0_212 = buffer.data(nsg0 + 212);
    const auto *nsg0_213 = buffer.data(nsg0 + 213);
    const auto *nsg0_214 = buffer.data(nsg0 + 214);
    const auto *nsg0_215 = buffer.data(nsg0 + 215);
    const auto *nsg0_219 = buffer.data(nsg0 + 219);
    const auto *nsg0_220 = buffer.data(nsg0 + 220);
    const auto *nsg0_221 = buffer.data(nsg0 + 221);
    const auto *nsg0_222 = buffer.data(nsg0 + 222);
    const auto *nsg0_223 = buffer.data(nsg0 + 223);
    const auto *nsg0_224 = buffer.data(nsg0 + 224);
    const auto *nsg0_225 = buffer.data(nsg0 + 225);
    const auto *nsg0_227 = buffer.data(nsg0 + 227);
    const auto *nsg0_228 = buffer.data(nsg0 + 228);
    const auto *nsg0_230 = buffer.data(nsg0 + 230);
    const auto *nsg0_231 = buffer.data(nsg0 + 231);
    const auto *nsg0_235 = buffer.data(nsg0 + 235);
    const auto *nsg0_236 = buffer.data(nsg0 + 236);
    const auto *nsg0_237 = buffer.data(nsg0 + 237);
    const auto *nsg0_239 = buffer.data(nsg0 + 239);
    const auto *nsg0_245 = buffer.data(nsg0 + 245);
    const auto *nsg0_249 = buffer.data(nsg0 + 249);
    const auto *nsg0_252 = buffer.data(nsg0 + 252);
    const auto *nsg0_253 = buffer.data(nsg0 + 253);
    const auto *nsg0_254 = buffer.data(nsg0 + 254);
    const auto *nsg0_255 = buffer.data(nsg0 + 255);

    const auto *nsg1_194 = buffer.data(nsg1 + 194);
    const auto *nsg1_205 = buffer.data(nsg1 + 205);
    const auto *nsg1_207 = buffer.data(nsg1 + 207);
    const auto *nsg1_208 = buffer.data(nsg1 + 208);
    const auto *nsg1_209 = buffer.data(nsg1 + 209);
    const auto *nsg1_210 = buffer.data(nsg1 + 210);
    const auto *nsg1_211 = buffer.data(nsg1 + 211);
    const auto *nsg1_212 = buffer.data(nsg1 + 212);
    const auto *nsg1_213 = buffer.data(nsg1 + 213);
    const auto *nsg1_214 = buffer.data(nsg1 + 214);
    const auto *nsg1_215 = buffer.data(nsg1 + 215);
    const auto *nsg1_219 = buffer.data(nsg1 + 219);
    const auto *nsg1_220 = buffer.data(nsg1 + 220);
    const auto *nsg1_221 = buffer.data(nsg1 + 221);
    const auto *nsg1_222 = buffer.data(nsg1 + 222);
    const auto *nsg1_223 = buffer.data(nsg1 + 223);
    const auto *nsg1_224 = buffer.data(nsg1 + 224);
    const auto *nsg1_225 = buffer.data(nsg1 + 225);
    const auto *nsg1_227 = buffer.data(nsg1 + 227);
    const auto *nsg1_228 = buffer.data(nsg1 + 228);
    const auto *nsg1_230 = buffer.data(nsg1 + 230);
    const auto *nsg1_231 = buffer.data(nsg1 + 231);
    const auto *nsg1_235 = buffer.data(nsg1 + 235);
    const auto *nsg1_236 = buffer.data(nsg1 + 236);
    const auto *nsg1_237 = buffer.data(nsg1 + 237);
    const auto *nsg1_239 = buffer.data(nsg1 + 239);
    const auto *nsg1_245 = buffer.data(nsg1 + 245);
    const auto *nsg1_249 = buffer.data(nsg1 + 249);
    const auto *nsg1_252 = buffer.data(nsg1 + 252);
    const auto *nsg1_253 = buffer.data(nsg1 + 253);
    const auto *nsg1_254 = buffer.data(nsg1 + 254);
    const auto *nsg1_255 = buffer.data(nsg1 + 255);

    const auto *nsh_272 = buffer.data(nsh + 272);
    const auto *nsh_273 = buffer.data(nsh + 273);
    const auto *nsh_275 = buffer.data(nsh + 275);
    const auto *nsh_276 = buffer.data(nsh + 276);
    const auto *nsh_278 = buffer.data(nsh + 278);
    const auto *nsh_279 = buffer.data(nsh + 279);
    const auto *nsh_282 = buffer.data(nsh + 282);
    const auto *nsh_288 = buffer.data(nsh + 288);
    const auto *nsh_289 = buffer.data(nsh + 289);
    const auto *nsh_290 = buffer.data(nsh + 290);
    const auto *nsh_291 = buffer.data(nsh + 291);
    const auto *nsh_292 = buffer.data(nsh + 292);
    const auto *nsh_293 = buffer.data(nsh + 293);
    const auto *nsh_294 = buffer.data(nsh + 294);
    const auto *nsh_295 = buffer.data(nsh + 295);
    const auto *nsh_296 = buffer.data(nsh + 296);
    const auto *nsh_297 = buffer.data(nsh + 297);
    const auto *nsh_298 = buffer.data(nsh + 298);
    const auto *nsh_299 = buffer.data(nsh + 299);
    const auto *nsh_300 = buffer.data(nsh + 300);
    const auto *nsh_301 = buffer.data(nsh + 301);
    const auto *nsh_302 = buffer.data(nsh + 302);
    const auto *nsh_303 = buffer.data(nsh + 303);
    const auto *nsh_308 = buffer.data(nsh + 308);
    const auto *nsh_309 = buffer.data(nsh + 309);
    const auto *nsh_310 = buffer.data(nsh + 310);
    const auto *nsh_311 = buffer.data(nsh + 311);
    const auto *nsh_312 = buffer.data(nsh + 312);
    const auto *nsh_313 = buffer.data(nsh + 313);
    const auto *nsh_314 = buffer.data(nsh + 314);
    const auto *nsh_315 = buffer.data(nsh + 315);
    const auto *nsh_316 = buffer.data(nsh + 316);
    const auto *nsh_317 = buffer.data(nsh + 317);
    const auto *nsh_318 = buffer.data(nsh + 318);
    const auto *nsh_320 = buffer.data(nsh + 320);
    const auto *nsh_321 = buffer.data(nsh + 321);
    const auto *nsh_322 = buffer.data(nsh + 322);
    const auto *nsh_324 = buffer.data(nsh + 324);
    const auto *nsh_325 = buffer.data(nsh + 325);
    const auto *nsh_330 = buffer.data(nsh + 330);
    const auto *nsh_331 = buffer.data(nsh + 331);
    const auto *nsh_332 = buffer.data(nsh + 332);
    const auto *nsh_333 = buffer.data(nsh + 333);
    const auto *nsh_334 = buffer.data(nsh + 334);
    const auto *nsh_335 = buffer.data(nsh + 335);
    const auto *nsh_336 = buffer.data(nsh + 336);
    const auto *nsh_338 = buffer.data(nsh + 338);
    const auto *nsh_339 = buffer.data(nsh + 339);
    const auto *nsh_341 = buffer.data(nsh + 341);
    const auto *nsh_342 = buffer.data(nsh + 342);
    const auto *nsh_345 = buffer.data(nsh + 345);
    const auto *nsh_350 = buffer.data(nsh + 350);
    const auto *nsh_351 = buffer.data(nsh + 351);
    const auto *nsh_352 = buffer.data(nsh + 352);
    const auto *nsh_353 = buffer.data(nsh + 353);
    const auto *nsh_354 = buffer.data(nsh + 354);
    const auto *nsh_355 = buffer.data(nsh + 355);
    const auto *nsh_356 = buffer.data(nsh + 356);
    const auto *nsh_357 = buffer.data(nsh + 357);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_y, pc_y, pc_z, msi0_252, msh_167, \
                         msh_188, msh_189, msi1_252, nsg0_194, nsg1_194, nsh_272, \
                         nsh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * msh_188[k]
                   + f_3 * pc_y[k] * nsh_272[k];

        t_363[k] = f_12 * msh_167[k]
                   + f_1 * nsg0_194[k]
                   - f_2 * nsg1_194[k]
                   + f_3 * pc_z[k] * nsh_272[k];

        t_364[k] = pa_y[k] * msi0_252[k]
                   - f_10 * pc_y[k] * msi1_252[k];

        t_365[k] = f_11 * msh_189[k]
                   + f_3 * pc_y[k] * nsh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pc_y, pc_z, msi0_255, msi0_257, \
                         msh_168, msh_190, msh_191, msi1_255, msi1_257, nsh_273, \
                         nsh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * msh_168[k]
                   + f_3 * pc_z[k] * nsh_273[k];

        t_367[k] = pa_y[k] * msi0_255[k]
                   + f_12 * msh_190[k]
                   - f_10 * pc_y[k] * msi1_255[k];

        t_368[k] = f_11 * msh_191[k]
                   + f_3 * pc_y[k] * nsh_275[k];

        t_369[k] = pa_y[k] * msi0_257[k]
                   - f_10 * pc_y[k] * msi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_y, pc_y, pc_z, msi0_258, msi0_261, \
                         msh_171, msh_192, msh_194, msi1_258, msi1_261, nsh_276, \
                         nsh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_y[k] * msi0_258[k]
                   + f_13 * msh_192[k]
                   - f_10 * pc_y[k] * msi1_258[k];

        t_371[k] = f_13 * msh_171[k]
                   + f_3 * pc_z[k] * nsh_276[k];

        t_372[k] = f_11 * msh_194[k]
                   + f_3 * pc_y[k] * nsh_278[k];

        t_373[k] = pa_y[k] * msi0_261[k]
                   - f_10 * pc_y[k] * msi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_y, pc_y, pc_z, msi0_262, msi0_264, msh_174, \
                         msh_195, msh_197, msi1_262, msi1_264, \
                         nsh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * msi0_262[k]
                   + f_14 * msh_195[k]
                   - f_10 * pc_y[k] * msi1_262[k];

        t_375[k] = f_13 * msh_174[k]
                   + f_3 * pc_z[k] * nsh_279[k];

        t_376[k] = pa_y[k] * msi0_264[k]
                   + f_12 * msh_197[k]
                   - f_10 * pc_y[k] * msi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pc_x, pc_y, msi0_266, msh_198, \
                         msh_288, msh_289, msi1_266, nsh_282, nsh_288, \
                         nsh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * msh_198[k]
                   + f_3 * pc_y[k] * nsh_282[k];

        t_378[k] = pa_y[k] * msi0_266[k]
                   - f_10 * pc_y[k] * msi1_266[k];

        t_379[k] = f_20 * msh_288[k]
                   + f_3 * pc_x[k] * nsh_288[k];

        t_380[k] = f_20 * msh_289[k]
                   + f_3 * pc_x[k] * nsh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, msh_290, msh_291, msh_292, msh_293, \
                         nsh_290, nsh_291, nsh_292, nsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_20 * msh_290[k]
                   + f_3 * pc_x[k] * nsh_290[k];

        t_382[k] = f_20 * msh_291[k]
                   + f_3 * pc_x[k] * nsh_291[k];

        t_383[k] = f_20 * msh_292[k]
                   + f_3 * pc_x[k] * nsh_292[k];

        t_384[k] = f_20 * msh_293[k]
                   + f_3 * pc_x[k] * nsh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, msh_183, msh_204, msh_206, nsg0_205, \
                         nsg0_207, nsg1_205, nsg1_207, nsh_288, \
                         nsh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * msh_204[k]
                   + f_1 * nsg0_205[k]
                   - f_2 * nsg1_205[k]
                   + f_3 * pc_y[k] * nsh_288[k];

        t_386[k] = f_13 * msh_183[k]
                   + f_3 * pc_z[k] * nsh_288[k];

        t_387[k] = f_11 * msh_206[k]
                   + f_8 * nsg0_207[k]
                   - f_9 * nsg1_207[k]
                   + f_3 * pc_y[k] * nsh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, msh_207, msh_208, msh_209, nsg0_208, \
                         nsg0_209, nsg1_208, nsg1_209, nsh_291, nsh_292, \
                         nsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * msh_207[k]
                   + f_6 * nsg0_208[k]
                   - f_7 * nsg1_208[k]
                   + f_3 * pc_y[k] * nsh_291[k];

        t_389[k] = f_11 * msh_208[k]
                   + f_4 * nsg0_209[k]
                   - f_5 * nsg1_209[k]
                   + f_3 * pc_y[k] * nsh_292[k];

        t_390[k] = f_11 * msh_209[k]
                   + f_3 * pc_y[k] * nsh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_x, pc_y, pc_z, msi0_279, \
                         msh_189, msh_294, msi1_279, nsg0_210, nsg1_210, \
                         nsh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_y[k] * msi0_279[k]
                   - f_10 * pc_y[k] * msi1_279[k];

        t_392[k] = f_20 * msh_294[k]
                   + f_1 * nsg0_210[k]
                   - f_2 * nsg1_210[k]
                   + f_3 * pc_x[k] * nsh_294[k];

        t_393[k] = f_3 * pc_y[k] * nsh_294[k];

        t_394[k] = f_14 * msh_189[k]
                   + f_3 * pc_z[k] * nsh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, msh_299, nsg0_210, nsg0_215, \
                         nsg1_210, nsg1_215, nsh_295, nsh_296, \
                         nsh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_4 * nsg0_210[k]
                   - f_5 * nsg1_210[k]
                   + f_3 * pc_y[k] * nsh_295[k];

        t_396[k] = f_3 * pc_y[k] * nsh_296[k];

        t_397[k] = f_20 * msh_299[k]
                   + f_8 * nsg0_215[k]
                   - f_9 * nsg1_215[k]
                   + f_3 * pc_x[k] * nsh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, nsg0_211, nsg0_212, nsg1_211, nsg1_212, \
                         nsh_297, nsh_298, nsh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_6 * nsg0_211[k]
                   - f_7 * nsg1_211[k]
                   + f_3 * pc_y[k] * nsh_297[k];

        t_399[k] = f_4 * nsg0_212[k]
                   - f_5 * nsg1_212[k]
                   + f_3 * pc_y[k] * nsh_298[k];

        t_400[k] = f_3 * pc_y[k] * nsh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_y, msh_303, nsg0_213, nsg0_214, \
                         nsg0_219, nsg1_213, nsg1_214, nsg1_219, nsh_300, nsh_301, \
                         nsh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_20 * msh_303[k]
                   + f_6 * nsg0_219[k]
                   - f_7 * nsg1_219[k]
                   + f_3 * pc_x[k] * nsh_303[k];

        t_402[k] = f_8 * nsg0_213[k]
                   - f_9 * nsg1_213[k]
                   + f_3 * pc_y[k] * nsh_300[k];

        t_403[k] = f_6 * nsg0_214[k]
                   - f_7 * nsg1_214[k]
                   + f_3 * pc_y[k] * nsh_301[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pc_x, pc_y, msh_308, msh_309, nsg0_215, \
                         nsg0_224, nsg1_215, nsg1_224, nsh_302, nsh_303, nsh_308, \
                         nsh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * nsg0_215[k]
                   - f_5 * nsg1_215[k]
                   + f_3 * pc_y[k] * nsh_302[k];

        t_405[k] = f_3 * pc_y[k] * nsh_303[k];

        t_406[k] = f_20 * msh_308[k]
                   + f_4 * nsg0_224[k]
                   - f_5 * nsg1_224[k]
                   + f_3 * pc_x[k] * nsh_308[k];

        t_407[k] = f_20 * msh_309[k]
                   + f_3 * pc_x[k] * nsh_309[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pc_x, pc_y, msh_310, msh_311, \
                         msh_312, msh_314, nsh_308, nsh_310, nsh_311, nsh_312, \
                         nsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_20 * msh_310[k]
                   + f_3 * pc_x[k] * nsh_310[k];

        t_409[k] = f_20 * msh_311[k]
                   + f_3 * pc_x[k] * nsh_311[k];

        t_410[k] = f_20 * msh_312[k]
                   + f_3 * pc_x[k] * nsh_312[k];

        t_411[k] = f_3 * pc_y[k] * nsh_308[k];

        t_412[k] = f_20 * msh_314[k]
                   + f_3 * pc_x[k] * nsh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pc_y, nsg0_220, nsg0_221, nsg0_222, nsg1_220, \
                         nsg1_221, nsg1_222, nsh_309, nsh_310, \
                         nsh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_1 * nsg0_220[k]
                   - f_2 * nsg1_220[k]
                   + f_3 * pc_y[k] * nsh_309[k];

        t_414[k] = f_16 * nsg0_221[k]
                   - f_17 * nsg1_221[k]
                   + f_3 * pc_y[k] * nsh_310[k];

        t_415[k] = f_8 * nsg0_222[k]
                   - f_9 * nsg1_222[k]
                   + f_3 * pc_y[k] * nsh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, msh_209, nsg0_223, nsg0_224, \
                         nsg1_223, nsg1_224, nsh_312, nsh_313, \
                         nsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * nsg0_223[k]
                   - f_7 * nsg1_223[k]
                   + f_3 * pc_y[k] * nsh_312[k];

        t_417[k] = f_4 * nsg0_224[k]
                   - f_5 * nsg1_224[k]
                   + f_3 * pc_y[k] * nsh_313[k];

        t_418[k] = f_3 * pc_y[k] * nsh_314[k];

        t_419[k] = f_14 * msh_209[k]
                   + f_1 * nsg0_224[k]
                   - f_2 * nsg1_224[k]
                   + f_3 * pc_z[k] * nsh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, msh_210, msh_315, \
                         msh_318, nsg0_225, nsg0_228, nsg1_225, nsg1_228, nsh_315, \
                         nsh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_21 * msh_315[k]
                   + f_1 * nsg0_225[k]
                   - f_2 * nsg1_225[k]
                   + f_3 * pc_x[k] * nsh_315[k];

        t_421[k] = f_21 * msh_210[k]
                   + f_3 * pc_y[k] * nsh_315[k];

        t_422[k] = f_3 * pc_z[k] * nsh_315[k];

        t_423[k] = f_21 * msh_318[k]
                   + f_8 * nsg0_228[k]
                   - f_9 * nsg1_228[k]
                   + f_3 * pc_x[k] * nsh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pc_x, pc_z, msh_321, nsg0_225, nsg0_231, \
                         nsg1_225, nsg1_231, nsh_316, nsh_317, nsh_318, \
                         nsh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_z[k] * nsh_316[k];

        t_425[k] = f_4 * nsg0_225[k]
                   - f_5 * nsg1_225[k]
                   + f_3 * pc_z[k] * nsh_317[k];

        t_426[k] = f_21 * msh_321[k]
                   + f_6 * nsg0_231[k]
                   - f_7 * nsg1_231[k]
                   + f_3 * pc_x[k] * nsh_321[k];

        t_427[k] = f_3 * pc_z[k] * nsh_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pc_x, pc_y, pc_z, msh_215, msh_325, \
                         nsg0_227, nsg0_235, nsg1_227, nsg1_235, nsh_320, nsh_321, \
                         nsh_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_21 * msh_215[k]
                   + f_3 * pc_y[k] * nsh_320[k];

        t_429[k] = f_6 * nsg0_227[k]
                   - f_7 * nsg1_227[k]
                   + f_3 * pc_z[k] * nsh_320[k];

        t_430[k] = f_21 * msh_325[k]
                   + f_4 * nsg0_235[k]
                   - f_5 * nsg1_235[k]
                   + f_3 * pc_x[k] * nsh_325[k];

        t_431[k] = f_3 * pc_z[k] * nsh_321[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, pc_y, pc_z, msh_219, msh_330, \
                         nsg0_228, nsg0_230, nsg1_228, nsg1_230, nsh_322, nsh_324, \
                         nsh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_4 * nsg0_228[k]
                   - f_5 * nsg1_228[k]
                   + f_3 * pc_z[k] * nsh_322[k];

        t_433[k] = f_21 * msh_219[k]
                   + f_3 * pc_y[k] * nsh_324[k];

        t_434[k] = f_8 * nsg0_230[k]
                   - f_9 * nsg1_230[k]
                   + f_3 * pc_z[k] * nsh_324[k];

        t_435[k] = f_21 * msh_330[k]
                   + f_3 * pc_x[k] * nsh_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_z, msh_332, msh_333, \
                         msh_334, msh_335, nsh_325, nsh_332, nsh_333, nsh_334, \
                         nsh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_3 * pc_z[k] * nsh_325[k];

        t_437[k] = f_21 * msh_332[k]
                   + f_3 * pc_x[k] * nsh_332[k];

        t_438[k] = f_21 * msh_333[k]
                   + f_3 * pc_x[k] * nsh_333[k];

        t_439[k] = f_21 * msh_334[k]
                   + f_3 * pc_x[k] * nsh_334[k];

        t_440[k] = f_21 * msh_335[k]
                   + f_3 * pc_x[k] * nsh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, msh_225, nsg0_235, nsg0_236, \
                         nsg1_235, nsg1_236, nsh_330, nsh_331, \
                         nsh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_21 * msh_225[k]
                   + f_1 * nsg0_235[k]
                   - f_2 * nsg1_235[k]
                   + f_3 * pc_y[k] * nsh_330[k];

        t_442[k] = f_3 * pc_z[k] * nsh_330[k];

        t_443[k] = f_4 * nsg0_235[k]
                   - f_5 * nsg1_235[k]
                   + f_3 * pc_z[k] * nsh_331[k];

        t_444[k] = f_6 * nsg0_236[k]
                   - f_7 * nsg1_236[k]
                   + f_3 * pc_z[k] * nsh_332[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pc_y, pc_z, msi0_280, msh_230, \
                         msi1_280, nsg0_237, nsg0_239, nsg1_237, nsg1_239, nsh_333, \
                         nsh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * nsg0_237[k]
                   - f_9 * nsg1_237[k]
                   + f_3 * pc_z[k] * nsh_333[k];

        t_446[k] = f_21 * msh_230[k]
                   + f_3 * pc_y[k] * nsh_335[k];

        t_447[k] = f_1 * nsg0_239[k]
                   - f_2 * nsg1_239[k]
                   + f_3 * pc_z[k] * nsh_335[k];

        t_448[k] = pa_z[k] * msi0_280[k]
                   - f_10 * pc_z[k] * msi1_280[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_z, pc_y, pc_z, msi0_283, msh_210, \
                         msh_231, msh_233, msi1_283, nsh_336, nsh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * msh_231[k]
                   + f_3 * pc_y[k] * nsh_336[k];

        t_450[k] = f_11 * msh_210[k]
                   + f_3 * pc_z[k] * nsh_336[k];

        t_451[k] = pa_z[k] * msi0_283[k]
                   - f_10 * pc_z[k] * msi1_283[k];

        t_452[k] = f_14 * msh_233[k]
                   + f_3 * pc_y[k] * nsh_338[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pa_z, pc_x, pc_z, msi0_286, msh_213, msh_341, \
                         msi1_286, nsg0_245, nsg1_245, nsh_339, \
                         nsh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_21 * msh_341[k]
                   + f_8 * nsg0_245[k]
                   - f_9 * nsg1_245[k]
                   + f_3 * pc_x[k] * nsh_341[k];

        t_454[k] = pa_z[k] * msi0_286[k]
                   - f_10 * pc_z[k] * msi1_286[k];

        t_455[k] = f_11 * msh_213[k]
                   + f_3 * pc_z[k] * nsh_339[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pa_z, pc_x, pc_y, pc_z, msi0_290, msh_236, \
                         msh_345, msi1_290, nsg0_249, nsg1_249, nsh_341, \
                         nsh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * msh_236[k]
                   + f_3 * pc_y[k] * nsh_341[k];

        t_457[k] = f_21 * msh_345[k]
                   + f_6 * nsg0_249[k]
                   - f_7 * nsg1_249[k]
                   + f_3 * pc_x[k] * nsh_345[k];

        t_458[k] = pa_z[k] * msi0_290[k]
                   - f_10 * pc_z[k] * msi1_290[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, pc_y, pc_z, msi0_292, msh_216, msh_217, \
                         msh_240, msi1_292, nsh_342, nsh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_11 * msh_216[k]
                   + f_3 * pc_z[k] * nsh_342[k];

        t_460[k] = pa_z[k] * msi0_292[k]
                   + f_12 * msh_217[k]
                   - f_10 * pc_z[k] * msi1_292[k];

        t_461[k] = f_14 * msh_240[k]
                   + f_3 * pc_y[k] * nsh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, msh_350, msh_351, msh_352, msh_353, \
                         nsg0_254, nsg1_254, nsh_350, nsh_351, nsh_352, \
                         nsh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_21 * msh_350[k]
                   + f_4 * nsg0_254[k]
                   - f_5 * nsg1_254[k]
                   + f_3 * pc_x[k] * nsh_350[k];

        t_463[k] = f_21 * msh_351[k]
                   + f_3 * pc_x[k] * nsh_351[k];

        t_464[k] = f_21 * msh_352[k]
                   + f_3 * pc_x[k] * nsh_352[k];

        t_465[k] = f_21 * msh_353[k]
                   + f_3 * pc_x[k] * nsh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pc_x, pc_z, msi0_301, msh_354, \
                         msh_355, msh_356, msi1_301, nsh_354, nsh_355, \
                         nsh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_21 * msh_354[k]
                   + f_3 * pc_x[k] * nsh_354[k];

        t_467[k] = f_21 * msh_355[k]
                   + f_3 * pc_x[k] * nsh_355[k];

        t_468[k] = f_21 * msh_356[k]
                   + f_3 * pc_x[k] * nsh_356[k];

        t_469[k] = pa_z[k] * msi0_301[k]
                   - f_10 * pc_z[k] * msi1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, msh_225, msh_248, msh_249, nsg0_252, \
                         nsg0_253, nsg1_252, nsg1_253, nsh_351, nsh_353, \
                         nsh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * msh_225[k]
                   + f_3 * pc_z[k] * nsh_351[k];

        t_471[k] = f_14 * msh_248[k]
                   + f_8 * nsg0_252[k]
                   - f_9 * nsg1_252[k]
                   + f_3 * pc_y[k] * nsh_353[k];

        t_472[k] = f_14 * msh_249[k]
                   + f_6 * nsg0_253[k]
                   - f_7 * nsg1_253[k]
                   + f_3 * pc_y[k] * nsh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, msh_230, msh_250, msh_251, nsg0_254, \
                         nsg1_254, nsh_355, nsh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * msh_250[k]
                   + f_4 * nsg0_254[k]
                   - f_5 * nsg1_254[k]
                   + f_3 * pc_y[k] * nsh_355[k];

        t_474[k] = f_14 * msh_251[k]
                   + f_3 * pc_y[k] * nsh_356[k];

        t_475[k] = f_11 * msh_230[k]
                   + f_1 * nsg0_254[k]
                   - f_2 * nsg1_254[k]
                   + f_3 * pc_z[k] * nsh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, msh_231, msh_252, msh_357, \
                         nsg0_255, nsg1_255, nsh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_21 * msh_357[k]
                   + f_1 * nsg0_255[k]
                   - f_2 * nsg1_255[k]
                   + f_3 * pc_x[k] * nsh_357[k];

        t_477[k] = f_13 * msh_252[k]
                   + f_3 * pc_y[k] * nsh_357[k];

        t_478[k] = f_12 * msh_231[k]
                   + f_3 * pc_z[k] * nsh_357[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_392 = buffer.data(msi0 + 392);
    const auto *msi0_395 = buffer.data(msi0 + 395);
    const auto *msi0_397 = buffer.data(msi0 + 397);
    const auto *msi0_398 = buffer.data(msi0 + 398);
    const auto *msi0_401 = buffer.data(msi0 + 401);
    const auto *msi0_402 = buffer.data(msi0 + 402);
    const auto *msi0_404 = buffer.data(msi0 + 404);
    const auto *msi0_406 = buffer.data(msi0 + 406);
    const auto *msi0_419 = buffer.data(msi0 + 419);

    const auto *msh_234 = buffer.data(msh + 234);
    const auto *msh_237 = buffer.data(msh + 237);
    const auto *msh_246 = buffer.data(msh + 246);
    const auto *msh_251 = buffer.data(msh + 251);
    const auto *msh_252 = buffer.data(msh + 252);
    const auto *msh_254 = buffer.data(msh + 254);
    const auto *msh_255 = buffer.data(msh + 255);
    const auto *msh_257 = buffer.data(msh + 257);
    const auto *msh_258 = buffer.data(msh + 258);
    const auto *msh_261 = buffer.data(msh + 261);
    const auto *msh_267 = buffer.data(msh + 267);
    const auto *msh_269 = buffer.data(msh + 269);
    const auto *msh_270 = buffer.data(msh + 270);
    const auto *msh_271 = buffer.data(msh + 271);
    const auto *msh_272 = buffer.data(msh + 272);
    const auto *msh_273 = buffer.data(msh + 273);
    const auto *msh_275 = buffer.data(msh + 275);
    const auto *msh_276 = buffer.data(msh + 276);
    const auto *msh_278 = buffer.data(msh + 278);
    const auto *msh_279 = buffer.data(msh + 279);
    const auto *msh_282 = buffer.data(msh + 282);
    const auto *msh_288 = buffer.data(msh + 288);
    const auto *msh_290 = buffer.data(msh + 290);
    const auto *msh_291 = buffer.data(msh + 291);
    const auto *msh_292 = buffer.data(msh + 292);
    const auto *msh_293 = buffer.data(msh + 293);
    const auto *msh_294 = buffer.data(msh + 294);
    const auto *msh_295 = buffer.data(msh + 295);
    const auto *msh_296 = buffer.data(msh + 296);
    const auto *msh_297 = buffer.data(msh + 297);
    const auto *msh_299 = buffer.data(msh + 299);
    const auto *msh_300 = buffer.data(msh + 300);
    const auto *msh_302 = buffer.data(msh + 302);
    const auto *msh_303 = buffer.data(msh + 303);
    const auto *msh_309 = buffer.data(msh + 309);
    const auto *msh_311 = buffer.data(msh + 311);
    const auto *msh_312 = buffer.data(msh + 312);
    const auto *msh_313 = buffer.data(msh + 313);
    const auto *msh_314 = buffer.data(msh + 314);
    const auto *msh_315 = buffer.data(msh + 315);
    const auto *msh_360 = buffer.data(msh + 360);
    const auto *msh_362 = buffer.data(msh + 362);
    const auto *msh_363 = buffer.data(msh + 363);
    const auto *msh_366 = buffer.data(msh + 366);
    const auto *msh_367 = buffer.data(msh + 367);
    const auto *msh_369 = buffer.data(msh + 369);
    const auto *msh_371 = buffer.data(msh + 371);
    const auto *msh_372 = buffer.data(msh + 372);
    const auto *msh_373 = buffer.data(msh + 373);
    const auto *msh_374 = buffer.data(msh + 374);
    const auto *msh_375 = buffer.data(msh + 375);
    const auto *msh_376 = buffer.data(msh + 376);
    const auto *msh_377 = buffer.data(msh + 377);
    const auto *msh_378 = buffer.data(msh + 378);
    const auto *msh_381 = buffer.data(msh + 381);
    const auto *msh_383 = buffer.data(msh + 383);
    const auto *msh_384 = buffer.data(msh + 384);
    const auto *msh_387 = buffer.data(msh + 387);
    const auto *msh_388 = buffer.data(msh + 388);
    const auto *msh_390 = buffer.data(msh + 390);
    const auto *msh_392 = buffer.data(msh + 392);
    const auto *msh_393 = buffer.data(msh + 393);
    const auto *msh_394 = buffer.data(msh + 394);
    const auto *msh_395 = buffer.data(msh + 395);
    const auto *msh_396 = buffer.data(msh + 396);
    const auto *msh_397 = buffer.data(msh + 397);
    const auto *msh_398 = buffer.data(msh + 398);
    const auto *msh_414 = buffer.data(msh + 414);
    const auto *msh_415 = buffer.data(msh + 415);
    const auto *msh_416 = buffer.data(msh + 416);
    const auto *msh_417 = buffer.data(msh + 417);
    const auto *msh_418 = buffer.data(msh + 418);
    const auto *msh_419 = buffer.data(msh + 419);
    const auto *msh_420 = buffer.data(msh + 420);
    const auto *msh_425 = buffer.data(msh + 425);
    const auto *msh_429 = buffer.data(msh + 429);
    const auto *msh_434 = buffer.data(msh + 434);
    const auto *msh_435 = buffer.data(msh + 435);
    const auto *msh_436 = buffer.data(msh + 436);
    const auto *msh_437 = buffer.data(msh + 437);
    const auto *msh_438 = buffer.data(msh + 438);
    const auto *msh_440 = buffer.data(msh + 440);
    const auto *msh_441 = buffer.data(msh + 441);
    const auto *msh_444 = buffer.data(msh + 444);

    const auto *msi1_392 = buffer.data(msi1 + 392);
    const auto *msi1_395 = buffer.data(msi1 + 395);
    const auto *msi1_397 = buffer.data(msi1 + 397);
    const auto *msi1_398 = buffer.data(msi1 + 398);
    const auto *msi1_401 = buffer.data(msi1 + 401);
    const auto *msi1_402 = buffer.data(msi1 + 402);
    const auto *msi1_404 = buffer.data(msi1 + 404);
    const auto *msi1_406 = buffer.data(msi1 + 406);
    const auto *msi1_419 = buffer.data(msi1 + 419);

    const auto *nsg0_258 = buffer.data(nsg0 + 258);
    const auto *nsg0_260 = buffer.data(nsg0 + 260);
    const auto *nsg0_261 = buffer.data(nsg0 + 261);
    const auto *nsg0_264 = buffer.data(nsg0 + 264);
    const auto *nsg0_265 = buffer.data(nsg0 + 265);
    const auto *nsg0_267 = buffer.data(nsg0 + 267);
    const auto *nsg0_268 = buffer.data(nsg0 + 268);
    const auto *nsg0_269 = buffer.data(nsg0 + 269);
    const auto *nsg0_270 = buffer.data(nsg0 + 270);
    const auto *nsg0_273 = buffer.data(nsg0 + 273);
    const auto *nsg0_275 = buffer.data(nsg0 + 275);
    const auto *nsg0_276 = buffer.data(nsg0 + 276);
    const auto *nsg0_279 = buffer.data(nsg0 + 279);
    const auto *nsg0_280 = buffer.data(nsg0 + 280);
    const auto *nsg0_282 = buffer.data(nsg0 + 282);
    const auto *nsg0_283 = buffer.data(nsg0 + 283);
    const auto *nsg0_284 = buffer.data(nsg0 + 284);
    const auto *nsg0_295 = buffer.data(nsg0 + 295);
    const auto *nsg0_297 = buffer.data(nsg0 + 297);
    const auto *nsg0_298 = buffer.data(nsg0 + 298);
    const auto *nsg0_299 = buffer.data(nsg0 + 299);
    const auto *nsg0_300 = buffer.data(nsg0 + 300);
    const auto *nsg0_301 = buffer.data(nsg0 + 301);
    const auto *nsg0_302 = buffer.data(nsg0 + 302);
    const auto *nsg0_303 = buffer.data(nsg0 + 303);
    const auto *nsg0_304 = buffer.data(nsg0 + 304);
    const auto *nsg0_305 = buffer.data(nsg0 + 305);
    const auto *nsg0_309 = buffer.data(nsg0 + 309);
    const auto *nsg0_310 = buffer.data(nsg0 + 310);
    const auto *nsg0_311 = buffer.data(nsg0 + 311);
    const auto *nsg0_312 = buffer.data(nsg0 + 312);
    const auto *nsg0_313 = buffer.data(nsg0 + 313);
    const auto *nsg0_314 = buffer.data(nsg0 + 314);
    const auto *nsg0_315 = buffer.data(nsg0 + 315);
    const auto *nsg0_318 = buffer.data(nsg0 + 318);

    const auto *nsg1_258 = buffer.data(nsg1 + 258);
    const auto *nsg1_260 = buffer.data(nsg1 + 260);
    const auto *nsg1_261 = buffer.data(nsg1 + 261);
    const auto *nsg1_264 = buffer.data(nsg1 + 264);
    const auto *nsg1_265 = buffer.data(nsg1 + 265);
    const auto *nsg1_267 = buffer.data(nsg1 + 267);
    const auto *nsg1_268 = buffer.data(nsg1 + 268);
    const auto *nsg1_269 = buffer.data(nsg1 + 269);
    const auto *nsg1_270 = buffer.data(nsg1 + 270);
    const auto *nsg1_273 = buffer.data(nsg1 + 273);
    const auto *nsg1_275 = buffer.data(nsg1 + 275);
    const auto *nsg1_276 = buffer.data(nsg1 + 276);
    const auto *nsg1_279 = buffer.data(nsg1 + 279);
    const auto *nsg1_280 = buffer.data(nsg1 + 280);
    const auto *nsg1_282 = buffer.data(nsg1 + 282);
    const auto *nsg1_283 = buffer.data(nsg1 + 283);
    const auto *nsg1_284 = buffer.data(nsg1 + 284);
    const auto *nsg1_295 = buffer.data(nsg1 + 295);
    const auto *nsg1_297 = buffer.data(nsg1 + 297);
    const auto *nsg1_298 = buffer.data(nsg1 + 298);
    const auto *nsg1_299 = buffer.data(nsg1 + 299);
    const auto *nsg1_300 = buffer.data(nsg1 + 300);
    const auto *nsg1_301 = buffer.data(nsg1 + 301);
    const auto *nsg1_302 = buffer.data(nsg1 + 302);
    const auto *nsg1_303 = buffer.data(nsg1 + 303);
    const auto *nsg1_304 = buffer.data(nsg1 + 304);
    const auto *nsg1_305 = buffer.data(nsg1 + 305);
    const auto *nsg1_309 = buffer.data(nsg1 + 309);
    const auto *nsg1_310 = buffer.data(nsg1 + 310);
    const auto *nsg1_311 = buffer.data(nsg1 + 311);
    const auto *nsg1_312 = buffer.data(nsg1 + 312);
    const auto *nsg1_313 = buffer.data(nsg1 + 313);
    const auto *nsg1_314 = buffer.data(nsg1 + 314);
    const auto *nsg1_315 = buffer.data(nsg1 + 315);
    const auto *nsg1_318 = buffer.data(nsg1 + 318);

    const auto *nsh_359 = buffer.data(nsh + 359);
    const auto *nsh_360 = buffer.data(nsh + 360);
    const auto *nsh_362 = buffer.data(nsh + 362);
    const auto *nsh_363 = buffer.data(nsh + 363);
    const auto *nsh_366 = buffer.data(nsh + 366);
    const auto *nsh_367 = buffer.data(nsh + 367);
    const auto *nsh_369 = buffer.data(nsh + 369);
    const auto *nsh_371 = buffer.data(nsh + 371);
    const auto *nsh_372 = buffer.data(nsh + 372);
    const auto *nsh_373 = buffer.data(nsh + 373);
    const auto *nsh_374 = buffer.data(nsh + 374);
    const auto *nsh_375 = buffer.data(nsh + 375);
    const auto *nsh_376 = buffer.data(nsh + 376);
    const auto *nsh_377 = buffer.data(nsh + 377);
    const auto *nsh_378 = buffer.data(nsh + 378);
    const auto *nsh_380 = buffer.data(nsh + 380);
    const auto *nsh_381 = buffer.data(nsh + 381);
    const auto *nsh_383 = buffer.data(nsh + 383);
    const auto *nsh_384 = buffer.data(nsh + 384);
    const auto *nsh_387 = buffer.data(nsh + 387);
    const auto *nsh_388 = buffer.data(nsh + 388);
    const auto *nsh_390 = buffer.data(nsh + 390);
    const auto *nsh_392 = buffer.data(nsh + 392);
    const auto *nsh_393 = buffer.data(nsh + 393);
    const auto *nsh_394 = buffer.data(nsh + 394);
    const auto *nsh_395 = buffer.data(nsh + 395);
    const auto *nsh_396 = buffer.data(nsh + 396);
    const auto *nsh_397 = buffer.data(nsh + 397);
    const auto *nsh_398 = buffer.data(nsh + 398);
    const auto *nsh_399 = buffer.data(nsh + 399);
    const auto *nsh_401 = buffer.data(nsh + 401);
    const auto *nsh_402 = buffer.data(nsh + 402);
    const auto *nsh_404 = buffer.data(nsh + 404);
    const auto *nsh_405 = buffer.data(nsh + 405);
    const auto *nsh_408 = buffer.data(nsh + 408);
    const auto *nsh_414 = buffer.data(nsh + 414);
    const auto *nsh_415 = buffer.data(nsh + 415);
    const auto *nsh_416 = buffer.data(nsh + 416);
    const auto *nsh_417 = buffer.data(nsh + 417);
    const auto *nsh_418 = buffer.data(nsh + 418);
    const auto *nsh_419 = buffer.data(nsh + 419);
    const auto *nsh_420 = buffer.data(nsh + 420);
    const auto *nsh_421 = buffer.data(nsh + 421);
    const auto *nsh_422 = buffer.data(nsh + 422);
    const auto *nsh_423 = buffer.data(nsh + 423);
    const auto *nsh_424 = buffer.data(nsh + 424);
    const auto *nsh_425 = buffer.data(nsh + 425);
    const auto *nsh_426 = buffer.data(nsh + 426);
    const auto *nsh_427 = buffer.data(nsh + 427);
    const auto *nsh_428 = buffer.data(nsh + 428);
    const auto *nsh_429 = buffer.data(nsh + 429);
    const auto *nsh_434 = buffer.data(nsh + 434);
    const auto *nsh_435 = buffer.data(nsh + 435);
    const auto *nsh_436 = buffer.data(nsh + 436);
    const auto *nsh_437 = buffer.data(nsh + 437);
    const auto *nsh_438 = buffer.data(nsh + 438);
    const auto *nsh_439 = buffer.data(nsh + 439);
    const auto *nsh_440 = buffer.data(nsh + 440);
    const auto *nsh_441 = buffer.data(nsh + 441);
    const auto *nsh_444 = buffer.data(nsh + 444);

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, msh_254, msh_360, msh_362, nsg0_258, \
                         nsg0_260, nsg1_258, nsg1_260, nsh_359, nsh_360, \
                         nsh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_21 * msh_360[k]
                   + f_8 * nsg0_258[k]
                   - f_9 * nsg1_258[k]
                   + f_3 * pc_x[k] * nsh_360[k];

        t_480[k] = f_13 * msh_254[k]
                   + f_3 * pc_y[k] * nsh_359[k];

        t_481[k] = f_21 * msh_362[k]
                   + f_8 * nsg0_260[k]
                   - f_9 * nsg1_260[k]
                   + f_3 * pc_x[k] * nsh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, msh_234, msh_257, msh_363, \
                         nsg0_261, nsg1_261, nsh_360, nsh_362, \
                         nsh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_21 * msh_363[k]
                   + f_6 * nsg0_261[k]
                   - f_7 * nsg1_261[k]
                   + f_3 * pc_x[k] * nsh_363[k];

        t_483[k] = f_12 * msh_234[k]
                   + f_3 * pc_z[k] * nsh_360[k];

        t_484[k] = f_13 * msh_257[k]
                   + f_3 * pc_y[k] * nsh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, msh_237, msh_366, msh_367, nsg0_264, \
                         nsg0_265, nsg1_264, nsg1_265, nsh_363, nsh_366, \
                         nsh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_21 * msh_366[k]
                   + f_6 * nsg0_264[k]
                   - f_7 * nsg1_264[k]
                   + f_3 * pc_x[k] * nsh_366[k];

        t_486[k] = f_21 * msh_367[k]
                   + f_4 * nsg0_265[k]
                   - f_5 * nsg1_265[k]
                   + f_3 * pc_x[k] * nsh_367[k];

        t_487[k] = f_12 * msh_237[k]
                   + f_3 * pc_z[k] * nsh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, msh_261, msh_369, msh_371, nsg0_267, \
                         nsg0_269, nsg1_267, nsg1_269, nsh_366, nsh_369, \
                         nsh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_21 * msh_369[k]
                   + f_4 * nsg0_267[k]
                   - f_5 * nsg1_267[k]
                   + f_3 * pc_x[k] * nsh_369[k];

        t_489[k] = f_13 * msh_261[k]
                   + f_3 * pc_y[k] * nsh_366[k];

        t_490[k] = f_21 * msh_371[k]
                   + f_4 * nsg0_269[k]
                   - f_5 * nsg1_269[k]
                   + f_3 * pc_x[k] * nsh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, msh_372, msh_373, msh_374, \
                         msh_375, msh_376, nsh_372, nsh_373, nsh_374, nsh_375, \
                         nsh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_21 * msh_372[k]
                   + f_3 * pc_x[k] * nsh_372[k];

        t_492[k] = f_21 * msh_373[k]
                   + f_3 * pc_x[k] * nsh_373[k];

        t_493[k] = f_21 * msh_374[k]
                   + f_3 * pc_x[k] * nsh_374[k];

        t_494[k] = f_21 * msh_375[k]
                   + f_3 * pc_x[k] * nsh_375[k];

        t_495[k] = f_21 * msh_376[k]
                   + f_3 * pc_x[k] * nsh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, msh_246, msh_267, msh_377, \
                         nsg0_265, nsg1_265, nsh_372, nsh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_21 * msh_377[k]
                   + f_3 * pc_x[k] * nsh_377[k];

        t_497[k] = f_13 * msh_267[k]
                   + f_1 * nsg0_265[k]
                   - f_2 * nsg1_265[k]
                   + f_3 * pc_y[k] * nsh_372[k];

        t_498[k] = f_12 * msh_246[k]
                   + f_3 * pc_z[k] * nsh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, msh_269, msh_270, msh_271, nsg0_267, \
                         nsg0_268, nsg0_269, nsg1_267, nsg1_268, nsg1_269, nsh_374, nsh_375, \
                         nsh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * msh_269[k]
                   + f_8 * nsg0_267[k]
                   - f_9 * nsg1_267[k]
                   + f_3 * pc_y[k] * nsh_374[k];

        t_500[k] = f_13 * msh_270[k]
                   + f_6 * nsg0_268[k]
                   - f_7 * nsg1_268[k]
                   + f_3 * pc_y[k] * nsh_375[k];

        t_501[k] = f_13 * msh_271[k]
                   + f_4 * nsg0_269[k]
                   - f_5 * nsg1_269[k]
                   + f_3 * pc_y[k] * nsh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, msh_251, msh_272, msh_378, \
                         nsg0_269, nsg0_270, nsg1_269, nsg1_270, nsh_377, \
                         nsh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * msh_272[k]
                   + f_3 * pc_y[k] * nsh_377[k];

        t_503[k] = f_12 * msh_251[k]
                   + f_1 * nsg0_269[k]
                   - f_2 * nsg1_269[k]
                   + f_3 * pc_z[k] * nsh_377[k];

        t_504[k] = f_21 * msh_378[k]
                   + f_1 * nsg0_270[k]
                   - f_2 * nsg1_270[k]
                   + f_3 * pc_x[k] * nsh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, msh_252, msh_273, \
                         msh_275, msh_381, nsg0_273, nsg1_273, nsh_378, nsh_380, \
                         nsh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * msh_273[k]
                   + f_3 * pc_y[k] * nsh_378[k];

        t_506[k] = f_13 * msh_252[k]
                   + f_3 * pc_z[k] * nsh_378[k];

        t_507[k] = f_21 * msh_381[k]
                   + f_8 * nsg0_273[k]
                   - f_9 * nsg1_273[k]
                   + f_3 * pc_x[k] * nsh_381[k];

        t_508[k] = f_12 * msh_275[k]
                   + f_3 * pc_y[k] * nsh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, msh_255, msh_383, msh_384, nsg0_275, \
                         nsg0_276, nsg1_275, nsg1_276, nsh_381, nsh_383, \
                         nsh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_21 * msh_383[k]
                   + f_8 * nsg0_275[k]
                   - f_9 * nsg1_275[k]
                   + f_3 * pc_x[k] * nsh_383[k];

        t_510[k] = f_21 * msh_384[k]
                   + f_6 * nsg0_276[k]
                   - f_7 * nsg1_276[k]
                   + f_3 * pc_x[k] * nsh_384[k];

        t_511[k] = f_13 * msh_255[k]
                   + f_3 * pc_z[k] * nsh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, msh_278, msh_387, msh_388, nsg0_279, \
                         nsg0_280, nsg1_279, nsg1_280, nsh_383, nsh_387, \
                         nsh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * msh_278[k]
                   + f_3 * pc_y[k] * nsh_383[k];

        t_513[k] = f_21 * msh_387[k]
                   + f_6 * nsg0_279[k]
                   - f_7 * nsg1_279[k]
                   + f_3 * pc_x[k] * nsh_387[k];

        t_514[k] = f_21 * msh_388[k]
                   + f_4 * nsg0_280[k]
                   - f_5 * nsg1_280[k]
                   + f_3 * pc_x[k] * nsh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, msh_258, msh_282, msh_390, \
                         nsg0_282, nsg1_282, nsh_384, nsh_387, \
                         nsh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * msh_258[k]
                   + f_3 * pc_z[k] * nsh_384[k];

        t_516[k] = f_21 * msh_390[k]
                   + f_4 * nsg0_282[k]
                   - f_5 * nsg1_282[k]
                   + f_3 * pc_x[k] * nsh_390[k];

        t_517[k] = f_12 * msh_282[k]
                   + f_3 * pc_y[k] * nsh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, msh_392, msh_393, msh_394, msh_395, \
                         nsg0_284, nsg1_284, nsh_392, nsh_393, nsh_394, \
                         nsh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_21 * msh_392[k]
                   + f_4 * nsg0_284[k]
                   - f_5 * nsg1_284[k]
                   + f_3 * pc_x[k] * nsh_392[k];

        t_519[k] = f_21 * msh_393[k]
                   + f_3 * pc_x[k] * nsh_393[k];

        t_520[k] = f_21 * msh_394[k]
                   + f_3 * pc_x[k] * nsh_394[k];

        t_521[k] = f_21 * msh_395[k]
                   + f_3 * pc_x[k] * nsh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, msh_288, msh_396, msh_397, \
                         msh_398, nsg0_280, nsg1_280, nsh_393, nsh_396, nsh_397, \
                         nsh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_21 * msh_396[k]
                   + f_3 * pc_x[k] * nsh_396[k];

        t_523[k] = f_21 * msh_397[k]
                   + f_3 * pc_x[k] * nsh_397[k];

        t_524[k] = f_21 * msh_398[k]
                   + f_3 * pc_x[k] * nsh_398[k];

        t_525[k] = f_12 * msh_288[k]
                   + f_1 * nsg0_280[k]
                   - f_2 * nsg1_280[k]
                   + f_3 * pc_y[k] * nsh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, msh_267, msh_290, msh_291, nsg0_282, \
                         nsg0_283, nsg1_282, nsg1_283, nsh_393, nsh_395, \
                         nsh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * msh_267[k]
                   + f_3 * pc_z[k] * nsh_393[k];

        t_527[k] = f_12 * msh_290[k]
                   + f_8 * nsg0_282[k]
                   - f_9 * nsg1_282[k]
                   + f_3 * pc_y[k] * nsh_395[k];

        t_528[k] = f_12 * msh_291[k]
                   + f_6 * nsg0_283[k]
                   - f_7 * nsg1_283[k]
                   + f_3 * pc_y[k] * nsh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pa_y, pc_y, pc_z, msi0_392, msh_272, \
                         msh_292, msh_293, msi1_392, nsg0_284, nsg1_284, nsh_397, \
                         nsh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * msh_292[k]
                   + f_4 * nsg0_284[k]
                   - f_5 * nsg1_284[k]
                   + f_3 * pc_y[k] * nsh_397[k];

        t_530[k] = f_12 * msh_293[k]
                   + f_3 * pc_y[k] * nsh_398[k];

        t_531[k] = f_13 * msh_272[k]
                   + f_1 * nsg0_284[k]
                   - f_2 * nsg1_284[k]
                   + f_3 * pc_z[k] * nsh_398[k];

        t_532[k] = pa_y[k] * msi0_392[k]
                   - f_10 * pc_y[k] * msi1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_y, pc_y, pc_z, msi0_395, msh_273, \
                         msh_294, msh_295, msh_296, msi1_395, nsh_399, \
                         nsh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * msh_294[k]
                   + f_3 * pc_y[k] * nsh_399[k];

        t_534[k] = f_14 * msh_273[k]
                   + f_3 * pc_z[k] * nsh_399[k];

        t_535[k] = pa_y[k] * msi0_395[k]
                   + f_12 * msh_295[k]
                   - f_10 * pc_y[k] * msi1_395[k];

        t_536[k] = f_11 * msh_296[k]
                   + f_3 * pc_y[k] * nsh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_y, pc_y, pc_z, msi0_397, msi0_398, \
                         msh_276, msh_297, msh_299, msi1_397, msi1_398, nsh_402, \
                         nsh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pa_y[k] * msi0_397[k]
                   - f_10 * pc_y[k] * msi1_397[k];

        t_538[k] = pa_y[k] * msi0_398[k]
                   + f_13 * msh_297[k]
                   - f_10 * pc_y[k] * msi1_398[k];

        t_539[k] = f_14 * msh_276[k]
                   + f_3 * pc_z[k] * nsh_402[k];

        t_540[k] = f_11 * msh_299[k]
                   + f_3 * pc_y[k] * nsh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pa_y, pc_y, pc_z, msi0_401, msi0_402, msh_279, \
                         msh_300, msi1_401, msi1_402, nsh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pa_y[k] * msi0_401[k]
                   - f_10 * pc_y[k] * msi1_401[k];

        t_542[k] = pa_y[k] * msi0_402[k]
                   + f_14 * msh_300[k]
                   - f_10 * pc_y[k] * msi1_402[k];

        t_543[k] = f_14 * msh_279[k]
                   + f_3 * pc_z[k] * nsh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_x, pc_y, msi0_404, msi0_406, \
                         msh_302, msh_303, msh_414, msi1_404, msi1_406, nsh_408, \
                         nsh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pa_y[k] * msi0_404[k]
                   + f_12 * msh_302[k]
                   - f_10 * pc_y[k] * msi1_404[k];

        t_545[k] = f_11 * msh_303[k]
                   + f_3 * pc_y[k] * nsh_408[k];

        t_546[k] = pa_y[k] * msi0_406[k]
                   - f_10 * pc_y[k] * msi1_406[k];

        t_547[k] = f_21 * msh_414[k]
                   + f_3 * pc_x[k] * nsh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, msh_415, msh_416, msh_417, \
                         msh_418, msh_419, nsh_415, nsh_416, nsh_417, nsh_418, \
                         nsh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_21 * msh_415[k]
                   + f_3 * pc_x[k] * nsh_415[k];

        t_549[k] = f_21 * msh_416[k]
                   + f_3 * pc_x[k] * nsh_416[k];

        t_550[k] = f_21 * msh_417[k]
                   + f_3 * pc_x[k] * nsh_417[k];

        t_551[k] = f_21 * msh_418[k]
                   + f_3 * pc_x[k] * nsh_418[k];

        t_552[k] = f_21 * msh_419[k]
                   + f_3 * pc_x[k] * nsh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, msh_288, msh_309, msh_311, nsg0_295, \
                         nsg0_297, nsg1_295, nsg1_297, nsh_414, \
                         nsh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * msh_309[k]
                   + f_1 * nsg0_295[k]
                   - f_2 * nsg1_295[k]
                   + f_3 * pc_y[k] * nsh_414[k];

        t_554[k] = f_14 * msh_288[k]
                   + f_3 * pc_z[k] * nsh_414[k];

        t_555[k] = f_11 * msh_311[k]
                   + f_8 * nsg0_297[k]
                   - f_9 * nsg1_297[k]
                   + f_3 * pc_y[k] * nsh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, msh_312, msh_313, msh_314, nsg0_298, \
                         nsg0_299, nsg1_298, nsg1_299, nsh_417, nsh_418, \
                         nsh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * msh_312[k]
                   + f_6 * nsg0_298[k]
                   - f_7 * nsg1_298[k]
                   + f_3 * pc_y[k] * nsh_417[k];

        t_557[k] = f_11 * msh_313[k]
                   + f_4 * nsg0_299[k]
                   - f_5 * nsg1_299[k]
                   + f_3 * pc_y[k] * nsh_418[k];

        t_558[k] = f_11 * msh_314[k]
                   + f_3 * pc_y[k] * nsh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_y, pc_x, pc_y, pc_z, msi0_419, \
                         msh_294, msh_420, msi1_419, nsg0_300, nsg1_300, \
                         nsh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pa_y[k] * msi0_419[k]
                   - f_10 * pc_y[k] * msi1_419[k];

        t_560[k] = f_21 * msh_420[k]
                   + f_1 * nsg0_300[k]
                   - f_2 * nsg1_300[k]
                   + f_3 * pc_x[k] * nsh_420[k];

        t_561[k] = f_3 * pc_y[k] * nsh_420[k];

        t_562[k] = f_21 * msh_294[k]
                   + f_3 * pc_z[k] * nsh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, msh_425, nsg0_300, nsg0_305, \
                         nsg1_300, nsg1_305, nsh_421, nsh_422, \
                         nsh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_4 * nsg0_300[k]
                   - f_5 * nsg1_300[k]
                   + f_3 * pc_y[k] * nsh_421[k];

        t_564[k] = f_3 * pc_y[k] * nsh_422[k];

        t_565[k] = f_21 * msh_425[k]
                   + f_8 * nsg0_305[k]
                   - f_9 * nsg1_305[k]
                   + f_3 * pc_x[k] * nsh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_y, nsg0_301, nsg0_302, nsg1_301, nsg1_302, \
                         nsh_423, nsh_424, nsh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_6 * nsg0_301[k]
                   - f_7 * nsg1_301[k]
                   + f_3 * pc_y[k] * nsh_423[k];

        t_567[k] = f_4 * nsg0_302[k]
                   - f_5 * nsg1_302[k]
                   + f_3 * pc_y[k] * nsh_424[k];

        t_568[k] = f_3 * pc_y[k] * nsh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_y, msh_429, nsg0_303, nsg0_304, \
                         nsg0_309, nsg1_303, nsg1_304, nsg1_309, nsh_426, nsh_427, \
                         nsh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_21 * msh_429[k]
                   + f_6 * nsg0_309[k]
                   - f_7 * nsg1_309[k]
                   + f_3 * pc_x[k] * nsh_429[k];

        t_570[k] = f_8 * nsg0_303[k]
                   - f_9 * nsg1_303[k]
                   + f_3 * pc_y[k] * nsh_426[k];

        t_571[k] = f_6 * nsg0_304[k]
                   - f_7 * nsg1_304[k]
                   + f_3 * pc_y[k] * nsh_427[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, msh_434, msh_435, nsg0_305, \
                         nsg0_314, nsg1_305, nsg1_314, nsh_428, nsh_429, nsh_434, \
                         nsh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_4 * nsg0_305[k]
                   - f_5 * nsg1_305[k]
                   + f_3 * pc_y[k] * nsh_428[k];

        t_573[k] = f_3 * pc_y[k] * nsh_429[k];

        t_574[k] = f_21 * msh_434[k]
                   + f_4 * nsg0_314[k]
                   - f_5 * nsg1_314[k]
                   + f_3 * pc_x[k] * nsh_434[k];

        t_575[k] = f_21 * msh_435[k]
                   + f_3 * pc_x[k] * nsh_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, pc_x, pc_y, msh_436, msh_437, \
                         msh_438, msh_440, nsh_434, nsh_436, nsh_437, nsh_438, \
                         nsh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_21 * msh_436[k]
                   + f_3 * pc_x[k] * nsh_436[k];

        t_577[k] = f_21 * msh_437[k]
                   + f_3 * pc_x[k] * nsh_437[k];

        t_578[k] = f_21 * msh_438[k]
                   + f_3 * pc_x[k] * nsh_438[k];

        t_579[k] = f_3 * pc_y[k] * nsh_434[k];

        t_580[k] = f_21 * msh_440[k]
                   + f_3 * pc_x[k] * nsh_440[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, nsg0_310, nsg0_311, nsg0_312, nsg1_310, \
                         nsg1_311, nsg1_312, nsh_435, nsh_436, \
                         nsh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_1 * nsg0_310[k]
                   - f_2 * nsg1_310[k]
                   + f_3 * pc_y[k] * nsh_435[k];

        t_582[k] = f_16 * nsg0_311[k]
                   - f_17 * nsg1_311[k]
                   + f_3 * pc_y[k] * nsh_436[k];

        t_583[k] = f_8 * nsg0_312[k]
                   - f_9 * nsg1_312[k]
                   + f_3 * pc_y[k] * nsh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, msh_314, nsg0_313, nsg0_314, \
                         nsg1_313, nsg1_314, nsh_438, nsh_439, \
                         nsh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * nsg0_313[k]
                   - f_7 * nsg1_313[k]
                   + f_3 * pc_y[k] * nsh_438[k];

        t_585[k] = f_4 * nsg0_314[k]
                   - f_5 * nsg1_314[k]
                   + f_3 * pc_y[k] * nsh_439[k];

        t_586[k] = f_3 * pc_y[k] * nsh_440[k];

        t_587[k] = f_21 * msh_314[k]
                   + f_1 * nsg0_314[k]
                   - f_2 * nsg1_314[k]
                   + f_3 * pc_z[k] * nsh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, msh_315, msh_441, \
                         msh_444, nsg0_315, nsg0_318, nsg1_315, nsg1_318, nsh_441, \
                         nsh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_14 * msh_441[k]
                   + f_1 * nsg0_315[k]
                   - f_2 * nsg1_315[k]
                   + f_3 * pc_x[k] * nsh_441[k];

        t_589[k] = f_20 * msh_315[k]
                   + f_3 * pc_y[k] * nsh_441[k];

        t_590[k] = f_3 * pc_z[k] * nsh_441[k];

        t_591[k] = f_14 * msh_444[k]
                   + f_8 * nsg0_318[k]
                   - f_9 * nsg1_318[k]
                   + f_3 * pc_x[k] * nsh_444[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_420 = buffer.data(msi0 + 420);
    const auto *msi0_423 = buffer.data(msi0 + 423);
    const auto *msi0_426 = buffer.data(msi0 + 426);
    const auto *msi0_430 = buffer.data(msi0 + 430);
    const auto *msi0_432 = buffer.data(msi0 + 432);
    const auto *msi0_441 = buffer.data(msi0 + 441);

    const auto *msh_315 = buffer.data(msh + 315);
    const auto *msh_318 = buffer.data(msh + 318);
    const auto *msh_320 = buffer.data(msh + 320);
    const auto *msh_321 = buffer.data(msh + 321);
    const auto *msh_322 = buffer.data(msh + 322);
    const auto *msh_324 = buffer.data(msh + 324);
    const auto *msh_330 = buffer.data(msh + 330);
    const auto *msh_335 = buffer.data(msh + 335);
    const auto *msh_336 = buffer.data(msh + 336);
    const auto *msh_338 = buffer.data(msh + 338);
    const auto *msh_339 = buffer.data(msh + 339);
    const auto *msh_341 = buffer.data(msh + 341);
    const auto *msh_342 = buffer.data(msh + 342);
    const auto *msh_345 = buffer.data(msh + 345);
    const auto *msh_351 = buffer.data(msh + 351);
    const auto *msh_353 = buffer.data(msh + 353);
    const auto *msh_354 = buffer.data(msh + 354);
    const auto *msh_355 = buffer.data(msh + 355);
    const auto *msh_356 = buffer.data(msh + 356);
    const auto *msh_357 = buffer.data(msh + 357);
    const auto *msh_359 = buffer.data(msh + 359);
    const auto *msh_360 = buffer.data(msh + 360);
    const auto *msh_362 = buffer.data(msh + 362);
    const auto *msh_363 = buffer.data(msh + 363);
    const auto *msh_366 = buffer.data(msh + 366);
    const auto *msh_372 = buffer.data(msh + 372);
    const auto *msh_374 = buffer.data(msh + 374);
    const auto *msh_375 = buffer.data(msh + 375);
    const auto *msh_376 = buffer.data(msh + 376);
    const auto *msh_377 = buffer.data(msh + 377);
    const auto *msh_378 = buffer.data(msh + 378);
    const auto *msh_380 = buffer.data(msh + 380);
    const auto *msh_383 = buffer.data(msh + 383);
    const auto *msh_387 = buffer.data(msh + 387);
    const auto *msh_393 = buffer.data(msh + 393);
    const auto *msh_395 = buffer.data(msh + 395);
    const auto *msh_396 = buffer.data(msh + 396);
    const auto *msh_397 = buffer.data(msh + 397);
    const auto *msh_398 = buffer.data(msh + 398);
    const auto *msh_399 = buffer.data(msh + 399);
    const auto *msh_447 = buffer.data(msh + 447);
    const auto *msh_451 = buffer.data(msh + 451);
    const auto *msh_456 = buffer.data(msh + 456);
    const auto *msh_458 = buffer.data(msh + 458);
    const auto *msh_459 = buffer.data(msh + 459);
    const auto *msh_460 = buffer.data(msh + 460);
    const auto *msh_461 = buffer.data(msh + 461);
    const auto *msh_467 = buffer.data(msh + 467);
    const auto *msh_471 = buffer.data(msh + 471);
    const auto *msh_476 = buffer.data(msh + 476);
    const auto *msh_477 = buffer.data(msh + 477);
    const auto *msh_478 = buffer.data(msh + 478);
    const auto *msh_479 = buffer.data(msh + 479);
    const auto *msh_480 = buffer.data(msh + 480);
    const auto *msh_481 = buffer.data(msh + 481);
    const auto *msh_482 = buffer.data(msh + 482);
    const auto *msh_483 = buffer.data(msh + 483);
    const auto *msh_486 = buffer.data(msh + 486);
    const auto *msh_488 = buffer.data(msh + 488);
    const auto *msh_489 = buffer.data(msh + 489);
    const auto *msh_492 = buffer.data(msh + 492);
    const auto *msh_493 = buffer.data(msh + 493);
    const auto *msh_495 = buffer.data(msh + 495);
    const auto *msh_497 = buffer.data(msh + 497);
    const auto *msh_498 = buffer.data(msh + 498);
    const auto *msh_499 = buffer.data(msh + 499);
    const auto *msh_500 = buffer.data(msh + 500);
    const auto *msh_501 = buffer.data(msh + 501);
    const auto *msh_502 = buffer.data(msh + 502);
    const auto *msh_503 = buffer.data(msh + 503);
    const auto *msh_504 = buffer.data(msh + 504);
    const auto *msh_507 = buffer.data(msh + 507);
    const auto *msh_509 = buffer.data(msh + 509);
    const auto *msh_510 = buffer.data(msh + 510);
    const auto *msh_513 = buffer.data(msh + 513);
    const auto *msh_514 = buffer.data(msh + 514);
    const auto *msh_516 = buffer.data(msh + 516);
    const auto *msh_518 = buffer.data(msh + 518);
    const auto *msh_519 = buffer.data(msh + 519);
    const auto *msh_520 = buffer.data(msh + 520);
    const auto *msh_521 = buffer.data(msh + 521);
    const auto *msh_522 = buffer.data(msh + 522);
    const auto *msh_523 = buffer.data(msh + 523);
    const auto *msh_524 = buffer.data(msh + 524);
    const auto *msh_525 = buffer.data(msh + 525);

    const auto *msi1_420 = buffer.data(msi1 + 420);
    const auto *msi1_423 = buffer.data(msi1 + 423);
    const auto *msi1_426 = buffer.data(msi1 + 426);
    const auto *msi1_430 = buffer.data(msi1 + 430);
    const auto *msi1_432 = buffer.data(msi1 + 432);
    const auto *msi1_441 = buffer.data(msi1 + 441);

    const auto *nsg0_315 = buffer.data(nsg0 + 315);
    const auto *nsg0_317 = buffer.data(nsg0 + 317);
    const auto *nsg0_318 = buffer.data(nsg0 + 318);
    const auto *nsg0_320 = buffer.data(nsg0 + 320);
    const auto *nsg0_321 = buffer.data(nsg0 + 321);
    const auto *nsg0_325 = buffer.data(nsg0 + 325);
    const auto *nsg0_326 = buffer.data(nsg0 + 326);
    const auto *nsg0_327 = buffer.data(nsg0 + 327);
    const auto *nsg0_329 = buffer.data(nsg0 + 329);
    const auto *nsg0_335 = buffer.data(nsg0 + 335);
    const auto *nsg0_339 = buffer.data(nsg0 + 339);
    const auto *nsg0_342 = buffer.data(nsg0 + 342);
    const auto *nsg0_343 = buffer.data(nsg0 + 343);
    const auto *nsg0_344 = buffer.data(nsg0 + 344);
    const auto *nsg0_345 = buffer.data(nsg0 + 345);
    const auto *nsg0_348 = buffer.data(nsg0 + 348);
    const auto *nsg0_350 = buffer.data(nsg0 + 350);
    const auto *nsg0_351 = buffer.data(nsg0 + 351);
    const auto *nsg0_354 = buffer.data(nsg0 + 354);
    const auto *nsg0_355 = buffer.data(nsg0 + 355);
    const auto *nsg0_357 = buffer.data(nsg0 + 357);
    const auto *nsg0_358 = buffer.data(nsg0 + 358);
    const auto *nsg0_359 = buffer.data(nsg0 + 359);
    const auto *nsg0_360 = buffer.data(nsg0 + 360);
    const auto *nsg0_363 = buffer.data(nsg0 + 363);
    const auto *nsg0_365 = buffer.data(nsg0 + 365);
    const auto *nsg0_366 = buffer.data(nsg0 + 366);
    const auto *nsg0_369 = buffer.data(nsg0 + 369);
    const auto *nsg0_370 = buffer.data(nsg0 + 370);
    const auto *nsg0_372 = buffer.data(nsg0 + 372);
    const auto *nsg0_373 = buffer.data(nsg0 + 373);
    const auto *nsg0_374 = buffer.data(nsg0 + 374);
    const auto *nsg0_375 = buffer.data(nsg0 + 375);

    const auto *nsg1_315 = buffer.data(nsg1 + 315);
    const auto *nsg1_317 = buffer.data(nsg1 + 317);
    const auto *nsg1_318 = buffer.data(nsg1 + 318);
    const auto *nsg1_320 = buffer.data(nsg1 + 320);
    const auto *nsg1_321 = buffer.data(nsg1 + 321);
    const auto *nsg1_325 = buffer.data(nsg1 + 325);
    const auto *nsg1_326 = buffer.data(nsg1 + 326);
    const auto *nsg1_327 = buffer.data(nsg1 + 327);
    const auto *nsg1_329 = buffer.data(nsg1 + 329);
    const auto *nsg1_335 = buffer.data(nsg1 + 335);
    const auto *nsg1_339 = buffer.data(nsg1 + 339);
    const auto *nsg1_342 = buffer.data(nsg1 + 342);
    const auto *nsg1_343 = buffer.data(nsg1 + 343);
    const auto *nsg1_344 = buffer.data(nsg1 + 344);
    const auto *nsg1_345 = buffer.data(nsg1 + 345);
    const auto *nsg1_348 = buffer.data(nsg1 + 348);
    const auto *nsg1_350 = buffer.data(nsg1 + 350);
    const auto *nsg1_351 = buffer.data(nsg1 + 351);
    const auto *nsg1_354 = buffer.data(nsg1 + 354);
    const auto *nsg1_355 = buffer.data(nsg1 + 355);
    const auto *nsg1_357 = buffer.data(nsg1 + 357);
    const auto *nsg1_358 = buffer.data(nsg1 + 358);
    const auto *nsg1_359 = buffer.data(nsg1 + 359);
    const auto *nsg1_360 = buffer.data(nsg1 + 360);
    const auto *nsg1_363 = buffer.data(nsg1 + 363);
    const auto *nsg1_365 = buffer.data(nsg1 + 365);
    const auto *nsg1_366 = buffer.data(nsg1 + 366);
    const auto *nsg1_369 = buffer.data(nsg1 + 369);
    const auto *nsg1_370 = buffer.data(nsg1 + 370);
    const auto *nsg1_372 = buffer.data(nsg1 + 372);
    const auto *nsg1_373 = buffer.data(nsg1 + 373);
    const auto *nsg1_374 = buffer.data(nsg1 + 374);
    const auto *nsg1_375 = buffer.data(nsg1 + 375);

    const auto *nsh_442 = buffer.data(nsh + 442);
    const auto *nsh_443 = buffer.data(nsh + 443);
    const auto *nsh_444 = buffer.data(nsh + 444);
    const auto *nsh_446 = buffer.data(nsh + 446);
    const auto *nsh_447 = buffer.data(nsh + 447);
    const auto *nsh_448 = buffer.data(nsh + 448);
    const auto *nsh_450 = buffer.data(nsh + 450);
    const auto *nsh_451 = buffer.data(nsh + 451);
    const auto *nsh_456 = buffer.data(nsh + 456);
    const auto *nsh_457 = buffer.data(nsh + 457);
    const auto *nsh_458 = buffer.data(nsh + 458);
    const auto *nsh_459 = buffer.data(nsh + 459);
    const auto *nsh_460 = buffer.data(nsh + 460);
    const auto *nsh_461 = buffer.data(nsh + 461);
    const auto *nsh_462 = buffer.data(nsh + 462);
    const auto *nsh_464 = buffer.data(nsh + 464);
    const auto *nsh_465 = buffer.data(nsh + 465);
    const auto *nsh_467 = buffer.data(nsh + 467);
    const auto *nsh_468 = buffer.data(nsh + 468);
    const auto *nsh_471 = buffer.data(nsh + 471);
    const auto *nsh_476 = buffer.data(nsh + 476);
    const auto *nsh_477 = buffer.data(nsh + 477);
    const auto *nsh_478 = buffer.data(nsh + 478);
    const auto *nsh_479 = buffer.data(nsh + 479);
    const auto *nsh_480 = buffer.data(nsh + 480);
    const auto *nsh_481 = buffer.data(nsh + 481);
    const auto *nsh_482 = buffer.data(nsh + 482);
    const auto *nsh_483 = buffer.data(nsh + 483);
    const auto *nsh_485 = buffer.data(nsh + 485);
    const auto *nsh_486 = buffer.data(nsh + 486);
    const auto *nsh_488 = buffer.data(nsh + 488);
    const auto *nsh_489 = buffer.data(nsh + 489);
    const auto *nsh_492 = buffer.data(nsh + 492);
    const auto *nsh_493 = buffer.data(nsh + 493);
    const auto *nsh_495 = buffer.data(nsh + 495);
    const auto *nsh_497 = buffer.data(nsh + 497);
    const auto *nsh_498 = buffer.data(nsh + 498);
    const auto *nsh_499 = buffer.data(nsh + 499);
    const auto *nsh_500 = buffer.data(nsh + 500);
    const auto *nsh_501 = buffer.data(nsh + 501);
    const auto *nsh_502 = buffer.data(nsh + 502);
    const auto *nsh_503 = buffer.data(nsh + 503);
    const auto *nsh_504 = buffer.data(nsh + 504);
    const auto *nsh_506 = buffer.data(nsh + 506);
    const auto *nsh_507 = buffer.data(nsh + 507);
    const auto *nsh_509 = buffer.data(nsh + 509);
    const auto *nsh_510 = buffer.data(nsh + 510);
    const auto *nsh_513 = buffer.data(nsh + 513);
    const auto *nsh_514 = buffer.data(nsh + 514);
    const auto *nsh_516 = buffer.data(nsh + 516);
    const auto *nsh_518 = buffer.data(nsh + 518);
    const auto *nsh_519 = buffer.data(nsh + 519);
    const auto *nsh_520 = buffer.data(nsh + 520);
    const auto *nsh_521 = buffer.data(nsh + 521);
    const auto *nsh_522 = buffer.data(nsh + 522);
    const auto *nsh_523 = buffer.data(nsh + 523);
    const auto *nsh_524 = buffer.data(nsh + 524);
    const auto *nsh_525 = buffer.data(nsh + 525);

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, msh_447, nsg0_315, nsg0_321, \
                         nsg1_315, nsg1_321, nsh_442, nsh_443, nsh_444, \
                         nsh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * nsh_442[k];

        t_593[k] = f_4 * nsg0_315[k]
                   - f_5 * nsg1_315[k]
                   + f_3 * pc_z[k] * nsh_443[k];

        t_594[k] = f_14 * msh_447[k]
                   + f_6 * nsg0_321[k]
                   - f_7 * nsg1_321[k]
                   + f_3 * pc_x[k] * nsh_447[k];

        t_595[k] = f_3 * pc_z[k] * nsh_444[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, msh_320, msh_451, \
                         nsg0_317, nsg0_325, nsg1_317, nsg1_325, nsh_446, nsh_447, \
                         nsh_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_20 * msh_320[k]
                   + f_3 * pc_y[k] * nsh_446[k];

        t_597[k] = f_6 * nsg0_317[k]
                   - f_7 * nsg1_317[k]
                   + f_3 * pc_z[k] * nsh_446[k];

        t_598[k] = f_14 * msh_451[k]
                   + f_4 * nsg0_325[k]
                   - f_5 * nsg1_325[k]
                   + f_3 * pc_x[k] * nsh_451[k];

        t_599[k] = f_3 * pc_z[k] * nsh_447[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, pc_z, msh_324, msh_456, \
                         nsg0_318, nsg0_320, nsg1_318, nsg1_320, nsh_448, nsh_450, \
                         nsh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_4 * nsg0_318[k]
                   - f_5 * nsg1_318[k]
                   + f_3 * pc_z[k] * nsh_448[k];

        t_601[k] = f_20 * msh_324[k]
                   + f_3 * pc_y[k] * nsh_450[k];

        t_602[k] = f_8 * nsg0_320[k]
                   - f_9 * nsg1_320[k]
                   + f_3 * pc_z[k] * nsh_450[k];

        t_603[k] = f_14 * msh_456[k]
                   + f_3 * pc_x[k] * nsh_456[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pc_x, pc_z, msh_458, msh_459, \
                         msh_460, msh_461, nsh_451, nsh_458, nsh_459, nsh_460, \
                         nsh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * nsh_451[k];

        t_605[k] = f_14 * msh_458[k]
                   + f_3 * pc_x[k] * nsh_458[k];

        t_606[k] = f_14 * msh_459[k]
                   + f_3 * pc_x[k] * nsh_459[k];

        t_607[k] = f_14 * msh_460[k]
                   + f_3 * pc_x[k] * nsh_460[k];

        t_608[k] = f_14 * msh_461[k]
                   + f_3 * pc_x[k] * nsh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pc_y, pc_z, msh_330, nsg0_325, nsg0_326, \
                         nsg1_325, nsg1_326, nsh_456, nsh_457, \
                         nsh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_20 * msh_330[k]
                   + f_1 * nsg0_325[k]
                   - f_2 * nsg1_325[k]
                   + f_3 * pc_y[k] * nsh_456[k];

        t_610[k] = f_3 * pc_z[k] * nsh_456[k];

        t_611[k] = f_4 * nsg0_325[k]
                   - f_5 * nsg1_325[k]
                   + f_3 * pc_z[k] * nsh_457[k];

        t_612[k] = f_6 * nsg0_326[k]
                   - f_7 * nsg1_326[k]
                   + f_3 * pc_z[k] * nsh_458[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pc_y, pc_z, msi0_420, msh_335, \
                         msi1_420, nsg0_327, nsg0_329, nsg1_327, nsg1_329, nsh_459, \
                         nsh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_8 * nsg0_327[k]
                   - f_9 * nsg1_327[k]
                   + f_3 * pc_z[k] * nsh_459[k];

        t_614[k] = f_20 * msh_335[k]
                   + f_3 * pc_y[k] * nsh_461[k];

        t_615[k] = f_1 * nsg0_329[k]
                   - f_2 * nsg1_329[k]
                   + f_3 * pc_z[k] * nsh_461[k];

        t_616[k] = pa_z[k] * msi0_420[k]
                   - f_10 * pc_z[k] * msi1_420[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pa_z, pc_y, pc_z, msi0_423, msh_315, \
                         msh_336, msh_338, msi1_423, nsh_462, nsh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_21 * msh_336[k]
                   + f_3 * pc_y[k] * nsh_462[k];

        t_618[k] = f_11 * msh_315[k]
                   + f_3 * pc_z[k] * nsh_462[k];

        t_619[k] = pa_z[k] * msi0_423[k]
                   - f_10 * pc_z[k] * msi1_423[k];

        t_620[k] = f_21 * msh_338[k]
                   + f_3 * pc_y[k] * nsh_464[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pa_z, pc_x, pc_z, msi0_426, msh_318, msh_467, \
                         msi1_426, nsg0_335, nsg1_335, nsh_465, \
                         nsh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_14 * msh_467[k]
                   + f_8 * nsg0_335[k]
                   - f_9 * nsg1_335[k]
                   + f_3 * pc_x[k] * nsh_467[k];

        t_622[k] = pa_z[k] * msi0_426[k]
                   - f_10 * pc_z[k] * msi1_426[k];

        t_623[k] = f_11 * msh_318[k]
                   + f_3 * pc_z[k] * nsh_465[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_z, pc_x, pc_y, pc_z, msi0_430, msh_341, \
                         msh_471, msi1_430, nsg0_339, nsg1_339, nsh_467, \
                         nsh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_21 * msh_341[k]
                   + f_3 * pc_y[k] * nsh_467[k];

        t_625[k] = f_14 * msh_471[k]
                   + f_6 * nsg0_339[k]
                   - f_7 * nsg1_339[k]
                   + f_3 * pc_x[k] * nsh_471[k];

        t_626[k] = pa_z[k] * msi0_430[k]
                   - f_10 * pc_z[k] * msi1_430[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_z, pc_y, pc_z, msi0_432, msh_321, msh_322, \
                         msh_345, msi1_432, nsh_468, nsh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_11 * msh_321[k]
                   + f_3 * pc_z[k] * nsh_468[k];

        t_628[k] = pa_z[k] * msi0_432[k]
                   + f_12 * msh_322[k]
                   - f_10 * pc_z[k] * msi1_432[k];

        t_629[k] = f_21 * msh_345[k]
                   + f_3 * pc_y[k] * nsh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, msh_476, msh_477, msh_478, msh_479, \
                         nsg0_344, nsg1_344, nsh_476, nsh_477, nsh_478, \
                         nsh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_14 * msh_476[k]
                   + f_4 * nsg0_344[k]
                   - f_5 * nsg1_344[k]
                   + f_3 * pc_x[k] * nsh_476[k];

        t_631[k] = f_14 * msh_477[k]
                   + f_3 * pc_x[k] * nsh_477[k];

        t_632[k] = f_14 * msh_478[k]
                   + f_3 * pc_x[k] * nsh_478[k];

        t_633[k] = f_14 * msh_479[k]
                   + f_3 * pc_x[k] * nsh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_z, pc_x, pc_z, msi0_441, msh_480, \
                         msh_481, msh_482, msi1_441, nsh_480, nsh_481, \
                         nsh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_14 * msh_480[k]
                   + f_3 * pc_x[k] * nsh_480[k];

        t_635[k] = f_14 * msh_481[k]
                   + f_3 * pc_x[k] * nsh_481[k];

        t_636[k] = f_14 * msh_482[k]
                   + f_3 * pc_x[k] * nsh_482[k];

        t_637[k] = pa_z[k] * msi0_441[k]
                   - f_10 * pc_z[k] * msi1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, msh_330, msh_353, msh_354, nsg0_342, \
                         nsg0_343, nsg1_342, nsg1_343, nsh_477, nsh_479, \
                         nsh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * msh_330[k]
                   + f_3 * pc_z[k] * nsh_477[k];

        t_639[k] = f_21 * msh_353[k]
                   + f_8 * nsg0_342[k]
                   - f_9 * nsg1_342[k]
                   + f_3 * pc_y[k] * nsh_479[k];

        t_640[k] = f_21 * msh_354[k]
                   + f_6 * nsg0_343[k]
                   - f_7 * nsg1_343[k]
                   + f_3 * pc_y[k] * nsh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, msh_335, msh_355, msh_356, nsg0_344, \
                         nsg1_344, nsh_481, nsh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_21 * msh_355[k]
                   + f_4 * nsg0_344[k]
                   - f_5 * nsg1_344[k]
                   + f_3 * pc_y[k] * nsh_481[k];

        t_642[k] = f_21 * msh_356[k]
                   + f_3 * pc_y[k] * nsh_482[k];

        t_643[k] = f_11 * msh_335[k]
                   + f_1 * nsg0_344[k]
                   - f_2 * nsg1_344[k]
                   + f_3 * pc_z[k] * nsh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, msh_336, msh_357, msh_483, \
                         nsg0_345, nsg1_345, nsh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_14 * msh_483[k]
                   + f_1 * nsg0_345[k]
                   - f_2 * nsg1_345[k]
                   + f_3 * pc_x[k] * nsh_483[k];

        t_645[k] = f_14 * msh_357[k]
                   + f_3 * pc_y[k] * nsh_483[k];

        t_646[k] = f_12 * msh_336[k]
                   + f_3 * pc_z[k] * nsh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, msh_359, msh_486, msh_488, nsg0_348, \
                         nsg0_350, nsg1_348, nsg1_350, nsh_485, nsh_486, \
                         nsh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_14 * msh_486[k]
                   + f_8 * nsg0_348[k]
                   - f_9 * nsg1_348[k]
                   + f_3 * pc_x[k] * nsh_486[k];

        t_648[k] = f_14 * msh_359[k]
                   + f_3 * pc_y[k] * nsh_485[k];

        t_649[k] = f_14 * msh_488[k]
                   + f_8 * nsg0_350[k]
                   - f_9 * nsg1_350[k]
                   + f_3 * pc_x[k] * nsh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, msh_339, msh_362, msh_489, \
                         nsg0_351, nsg1_351, nsh_486, nsh_488, \
                         nsh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_14 * msh_489[k]
                   + f_6 * nsg0_351[k]
                   - f_7 * nsg1_351[k]
                   + f_3 * pc_x[k] * nsh_489[k];

        t_651[k] = f_12 * msh_339[k]
                   + f_3 * pc_z[k] * nsh_486[k];

        t_652[k] = f_14 * msh_362[k]
                   + f_3 * pc_y[k] * nsh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, msh_342, msh_492, msh_493, nsg0_354, \
                         nsg0_355, nsg1_354, nsg1_355, nsh_489, nsh_492, \
                         nsh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_14 * msh_492[k]
                   + f_6 * nsg0_354[k]
                   - f_7 * nsg1_354[k]
                   + f_3 * pc_x[k] * nsh_492[k];

        t_654[k] = f_14 * msh_493[k]
                   + f_4 * nsg0_355[k]
                   - f_5 * nsg1_355[k]
                   + f_3 * pc_x[k] * nsh_493[k];

        t_655[k] = f_12 * msh_342[k]
                   + f_3 * pc_z[k] * nsh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, msh_366, msh_495, msh_497, nsg0_357, \
                         nsg0_359, nsg1_357, nsg1_359, nsh_492, nsh_495, \
                         nsh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_14 * msh_495[k]
                   + f_4 * nsg0_357[k]
                   - f_5 * nsg1_357[k]
                   + f_3 * pc_x[k] * nsh_495[k];

        t_657[k] = f_14 * msh_366[k]
                   + f_3 * pc_y[k] * nsh_492[k];

        t_658[k] = f_14 * msh_497[k]
                   + f_4 * nsg0_359[k]
                   - f_5 * nsg1_359[k]
                   + f_3 * pc_x[k] * nsh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, msh_498, msh_499, msh_500, \
                         msh_501, msh_502, nsh_498, nsh_499, nsh_500, nsh_501, \
                         nsh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_14 * msh_498[k]
                   + f_3 * pc_x[k] * nsh_498[k];

        t_660[k] = f_14 * msh_499[k]
                   + f_3 * pc_x[k] * nsh_499[k];

        t_661[k] = f_14 * msh_500[k]
                   + f_3 * pc_x[k] * nsh_500[k];

        t_662[k] = f_14 * msh_501[k]
                   + f_3 * pc_x[k] * nsh_501[k];

        t_663[k] = f_14 * msh_502[k]
                   + f_3 * pc_x[k] * nsh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, msh_351, msh_372, msh_503, \
                         nsg0_355, nsg1_355, nsh_498, nsh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_14 * msh_503[k]
                   + f_3 * pc_x[k] * nsh_503[k];

        t_665[k] = f_14 * msh_372[k]
                   + f_1 * nsg0_355[k]
                   - f_2 * nsg1_355[k]
                   + f_3 * pc_y[k] * nsh_498[k];

        t_666[k] = f_12 * msh_351[k]
                   + f_3 * pc_z[k] * nsh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, msh_374, msh_375, msh_376, nsg0_357, \
                         nsg0_358, nsg0_359, nsg1_357, nsg1_358, nsg1_359, nsh_500, nsh_501, \
                         nsh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * msh_374[k]
                   + f_8 * nsg0_357[k]
                   - f_9 * nsg1_357[k]
                   + f_3 * pc_y[k] * nsh_500[k];

        t_668[k] = f_14 * msh_375[k]
                   + f_6 * nsg0_358[k]
                   - f_7 * nsg1_358[k]
                   + f_3 * pc_y[k] * nsh_501[k];

        t_669[k] = f_14 * msh_376[k]
                   + f_4 * nsg0_359[k]
                   - f_5 * nsg1_359[k]
                   + f_3 * pc_y[k] * nsh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, msh_356, msh_377, msh_504, \
                         nsg0_359, nsg0_360, nsg1_359, nsg1_360, nsh_503, \
                         nsh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * msh_377[k]
                   + f_3 * pc_y[k] * nsh_503[k];

        t_671[k] = f_12 * msh_356[k]
                   + f_1 * nsg0_359[k]
                   - f_2 * nsg1_359[k]
                   + f_3 * pc_z[k] * nsh_503[k];

        t_672[k] = f_14 * msh_504[k]
                   + f_1 * nsg0_360[k]
                   - f_2 * nsg1_360[k]
                   + f_3 * pc_x[k] * nsh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, msh_357, msh_378, \
                         msh_380, msh_507, nsg0_363, nsg1_363, nsh_504, nsh_506, \
                         nsh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * msh_378[k]
                   + f_3 * pc_y[k] * nsh_504[k];

        t_674[k] = f_13 * msh_357[k]
                   + f_3 * pc_z[k] * nsh_504[k];

        t_675[k] = f_14 * msh_507[k]
                   + f_8 * nsg0_363[k]
                   - f_9 * nsg1_363[k]
                   + f_3 * pc_x[k] * nsh_507[k];

        t_676[k] = f_13 * msh_380[k]
                   + f_3 * pc_y[k] * nsh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, msh_360, msh_509, msh_510, nsg0_365, \
                         nsg0_366, nsg1_365, nsg1_366, nsh_507, nsh_509, \
                         nsh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_14 * msh_509[k]
                   + f_8 * nsg0_365[k]
                   - f_9 * nsg1_365[k]
                   + f_3 * pc_x[k] * nsh_509[k];

        t_678[k] = f_14 * msh_510[k]
                   + f_6 * nsg0_366[k]
                   - f_7 * nsg1_366[k]
                   + f_3 * pc_x[k] * nsh_510[k];

        t_679[k] = f_13 * msh_360[k]
                   + f_3 * pc_z[k] * nsh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, msh_383, msh_513, msh_514, nsg0_369, \
                         nsg0_370, nsg1_369, nsg1_370, nsh_509, nsh_513, \
                         nsh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * msh_383[k]
                   + f_3 * pc_y[k] * nsh_509[k];

        t_681[k] = f_14 * msh_513[k]
                   + f_6 * nsg0_369[k]
                   - f_7 * nsg1_369[k]
                   + f_3 * pc_x[k] * nsh_513[k];

        t_682[k] = f_14 * msh_514[k]
                   + f_4 * nsg0_370[k]
                   - f_5 * nsg1_370[k]
                   + f_3 * pc_x[k] * nsh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, msh_363, msh_387, msh_516, \
                         nsg0_372, nsg1_372, nsh_510, nsh_513, \
                         nsh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * msh_363[k]
                   + f_3 * pc_z[k] * nsh_510[k];

        t_684[k] = f_14 * msh_516[k]
                   + f_4 * nsg0_372[k]
                   - f_5 * nsg1_372[k]
                   + f_3 * pc_x[k] * nsh_516[k];

        t_685[k] = f_13 * msh_387[k]
                   + f_3 * pc_y[k] * nsh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, msh_518, msh_519, msh_520, msh_521, \
                         nsg0_374, nsg1_374, nsh_518, nsh_519, nsh_520, \
                         nsh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_14 * msh_518[k]
                   + f_4 * nsg0_374[k]
                   - f_5 * nsg1_374[k]
                   + f_3 * pc_x[k] * nsh_518[k];

        t_687[k] = f_14 * msh_519[k]
                   + f_3 * pc_x[k] * nsh_519[k];

        t_688[k] = f_14 * msh_520[k]
                   + f_3 * pc_x[k] * nsh_520[k];

        t_689[k] = f_14 * msh_521[k]
                   + f_3 * pc_x[k] * nsh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, msh_393, msh_522, msh_523, \
                         msh_524, nsg0_370, nsg1_370, nsh_519, nsh_522, nsh_523, \
                         nsh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_14 * msh_522[k]
                   + f_3 * pc_x[k] * nsh_522[k];

        t_691[k] = f_14 * msh_523[k]
                   + f_3 * pc_x[k] * nsh_523[k];

        t_692[k] = f_14 * msh_524[k]
                   + f_3 * pc_x[k] * nsh_524[k];

        t_693[k] = f_13 * msh_393[k]
                   + f_1 * nsg0_370[k]
                   - f_2 * nsg1_370[k]
                   + f_3 * pc_y[k] * nsh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, msh_372, msh_395, msh_396, nsg0_372, \
                         nsg0_373, nsg1_372, nsg1_373, nsh_519, nsh_521, \
                         nsh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * msh_372[k]
                   + f_3 * pc_z[k] * nsh_519[k];

        t_695[k] = f_13 * msh_395[k]
                   + f_8 * nsg0_372[k]
                   - f_9 * nsg1_372[k]
                   + f_3 * pc_y[k] * nsh_521[k];

        t_696[k] = f_13 * msh_396[k]
                   + f_6 * nsg0_373[k]
                   - f_7 * nsg1_373[k]
                   + f_3 * pc_y[k] * nsh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, msh_377, msh_397, msh_398, nsg0_374, \
                         nsg1_374, nsh_523, nsh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * msh_397[k]
                   + f_4 * nsg0_374[k]
                   - f_5 * nsg1_374[k]
                   + f_3 * pc_y[k] * nsh_523[k];

        t_698[k] = f_13 * msh_398[k]
                   + f_3 * pc_y[k] * nsh_524[k];

        t_699[k] = f_13 * msh_377[k]
                   + f_1 * nsg0_374[k]
                   - f_2 * nsg1_374[k]
                   + f_3 * pc_z[k] * nsh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, msh_378, msh_399, msh_525, \
                         nsg0_375, nsg1_375, nsh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_14 * msh_525[k]
                   + f_1 * nsg0_375[k]
                   - f_2 * nsg1_375[k]
                   + f_3 * pc_x[k] * nsh_525[k];

        t_701[k] = f_12 * msh_399[k]
                   + f_3 * pc_y[k] * nsh_525[k];

        t_702[k] = f_14 * msh_378[k]
                   + f_3 * pc_z[k] * nsh_525[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_560 = buffer.data(msi0 + 560);
    const auto *msi0_563 = buffer.data(msi0 + 563);
    const auto *msi0_565 = buffer.data(msi0 + 565);
    const auto *msi0_566 = buffer.data(msi0 + 566);
    const auto *msi0_569 = buffer.data(msi0 + 569);
    const auto *msi0_570 = buffer.data(msi0 + 570);
    const auto *msi0_572 = buffer.data(msi0 + 572);
    const auto *msi0_574 = buffer.data(msi0 + 574);
    const auto *msi0_587 = buffer.data(msi0 + 587);
    const auto *msi0_588 = buffer.data(msi0 + 588);
    const auto *msi0_591 = buffer.data(msi0 + 591);
    const auto *msi0_594 = buffer.data(msi0 + 594);

    const auto *msh_381 = buffer.data(msh + 381);
    const auto *msh_384 = buffer.data(msh + 384);
    const auto *msh_393 = buffer.data(msh + 393);
    const auto *msh_398 = buffer.data(msh + 398);
    const auto *msh_399 = buffer.data(msh + 399);
    const auto *msh_401 = buffer.data(msh + 401);
    const auto *msh_402 = buffer.data(msh + 402);
    const auto *msh_404 = buffer.data(msh + 404);
    const auto *msh_405 = buffer.data(msh + 405);
    const auto *msh_408 = buffer.data(msh + 408);
    const auto *msh_414 = buffer.data(msh + 414);
    const auto *msh_416 = buffer.data(msh + 416);
    const auto *msh_417 = buffer.data(msh + 417);
    const auto *msh_418 = buffer.data(msh + 418);
    const auto *msh_419 = buffer.data(msh + 419);
    const auto *msh_420 = buffer.data(msh + 420);
    const auto *msh_421 = buffer.data(msh + 421);
    const auto *msh_422 = buffer.data(msh + 422);
    const auto *msh_423 = buffer.data(msh + 423);
    const auto *msh_425 = buffer.data(msh + 425);
    const auto *msh_426 = buffer.data(msh + 426);
    const auto *msh_428 = buffer.data(msh + 428);
    const auto *msh_429 = buffer.data(msh + 429);
    const auto *msh_435 = buffer.data(msh + 435);
    const auto *msh_437 = buffer.data(msh + 437);
    const auto *msh_438 = buffer.data(msh + 438);
    const auto *msh_439 = buffer.data(msh + 439);
    const auto *msh_440 = buffer.data(msh + 440);
    const auto *msh_441 = buffer.data(msh + 441);
    const auto *msh_444 = buffer.data(msh + 444);
    const auto *msh_446 = buffer.data(msh + 446);
    const auto *msh_450 = buffer.data(msh + 450);
    const auto *msh_456 = buffer.data(msh + 456);
    const auto *msh_461 = buffer.data(msh + 461);
    const auto *msh_462 = buffer.data(msh + 462);
    const auto *msh_464 = buffer.data(msh + 464);
    const auto *msh_528 = buffer.data(msh + 528);
    const auto *msh_530 = buffer.data(msh + 530);
    const auto *msh_531 = buffer.data(msh + 531);
    const auto *msh_534 = buffer.data(msh + 534);
    const auto *msh_535 = buffer.data(msh + 535);
    const auto *msh_537 = buffer.data(msh + 537);
    const auto *msh_539 = buffer.data(msh + 539);
    const auto *msh_540 = buffer.data(msh + 540);
    const auto *msh_541 = buffer.data(msh + 541);
    const auto *msh_542 = buffer.data(msh + 542);
    const auto *msh_543 = buffer.data(msh + 543);
    const auto *msh_544 = buffer.data(msh + 544);
    const auto *msh_545 = buffer.data(msh + 545);
    const auto *msh_561 = buffer.data(msh + 561);
    const auto *msh_562 = buffer.data(msh + 562);
    const auto *msh_563 = buffer.data(msh + 563);
    const auto *msh_564 = buffer.data(msh + 564);
    const auto *msh_565 = buffer.data(msh + 565);
    const auto *msh_566 = buffer.data(msh + 566);
    const auto *msh_567 = buffer.data(msh + 567);
    const auto *msh_572 = buffer.data(msh + 572);
    const auto *msh_576 = buffer.data(msh + 576);
    const auto *msh_581 = buffer.data(msh + 581);
    const auto *msh_582 = buffer.data(msh + 582);
    const auto *msh_583 = buffer.data(msh + 583);
    const auto *msh_584 = buffer.data(msh + 584);
    const auto *msh_585 = buffer.data(msh + 585);
    const auto *msh_587 = buffer.data(msh + 587);
    const auto *msh_588 = buffer.data(msh + 588);
    const auto *msh_591 = buffer.data(msh + 591);
    const auto *msh_594 = buffer.data(msh + 594);
    const auto *msh_598 = buffer.data(msh + 598);
    const auto *msh_603 = buffer.data(msh + 603);
    const auto *msh_605 = buffer.data(msh + 605);
    const auto *msh_606 = buffer.data(msh + 606);
    const auto *msh_607 = buffer.data(msh + 607);
    const auto *msh_608 = buffer.data(msh + 608);
    const auto *msh_614 = buffer.data(msh + 614);

    const auto *msi1_560 = buffer.data(msi1 + 560);
    const auto *msi1_563 = buffer.data(msi1 + 563);
    const auto *msi1_565 = buffer.data(msi1 + 565);
    const auto *msi1_566 = buffer.data(msi1 + 566);
    const auto *msi1_569 = buffer.data(msi1 + 569);
    const auto *msi1_570 = buffer.data(msi1 + 570);
    const auto *msi1_572 = buffer.data(msi1 + 572);
    const auto *msi1_574 = buffer.data(msi1 + 574);
    const auto *msi1_587 = buffer.data(msi1 + 587);
    const auto *msi1_588 = buffer.data(msi1 + 588);
    const auto *msi1_591 = buffer.data(msi1 + 591);
    const auto *msi1_594 = buffer.data(msi1 + 594);

    const auto *nsg0_378 = buffer.data(nsg0 + 378);
    const auto *nsg0_380 = buffer.data(nsg0 + 380);
    const auto *nsg0_381 = buffer.data(nsg0 + 381);
    const auto *nsg0_384 = buffer.data(nsg0 + 384);
    const auto *nsg0_385 = buffer.data(nsg0 + 385);
    const auto *nsg0_387 = buffer.data(nsg0 + 387);
    const auto *nsg0_388 = buffer.data(nsg0 + 388);
    const auto *nsg0_389 = buffer.data(nsg0 + 389);
    const auto *nsg0_400 = buffer.data(nsg0 + 400);
    const auto *nsg0_402 = buffer.data(nsg0 + 402);
    const auto *nsg0_403 = buffer.data(nsg0 + 403);
    const auto *nsg0_404 = buffer.data(nsg0 + 404);
    const auto *nsg0_405 = buffer.data(nsg0 + 405);
    const auto *nsg0_406 = buffer.data(nsg0 + 406);
    const auto *nsg0_407 = buffer.data(nsg0 + 407);
    const auto *nsg0_408 = buffer.data(nsg0 + 408);
    const auto *nsg0_409 = buffer.data(nsg0 + 409);
    const auto *nsg0_410 = buffer.data(nsg0 + 410);
    const auto *nsg0_414 = buffer.data(nsg0 + 414);
    const auto *nsg0_415 = buffer.data(nsg0 + 415);
    const auto *nsg0_416 = buffer.data(nsg0 + 416);
    const auto *nsg0_417 = buffer.data(nsg0 + 417);
    const auto *nsg0_418 = buffer.data(nsg0 + 418);
    const auto *nsg0_419 = buffer.data(nsg0 + 419);
    const auto *nsg0_420 = buffer.data(nsg0 + 420);
    const auto *nsg0_422 = buffer.data(nsg0 + 422);
    const auto *nsg0_423 = buffer.data(nsg0 + 423);
    const auto *nsg0_425 = buffer.data(nsg0 + 425);
    const auto *nsg0_426 = buffer.data(nsg0 + 426);
    const auto *nsg0_430 = buffer.data(nsg0 + 430);
    const auto *nsg0_431 = buffer.data(nsg0 + 431);
    const auto *nsg0_432 = buffer.data(nsg0 + 432);
    const auto *nsg0_434 = buffer.data(nsg0 + 434);
    const auto *nsg0_440 = buffer.data(nsg0 + 440);

    const auto *nsg1_378 = buffer.data(nsg1 + 378);
    const auto *nsg1_380 = buffer.data(nsg1 + 380);
    const auto *nsg1_381 = buffer.data(nsg1 + 381);
    const auto *nsg1_384 = buffer.data(nsg1 + 384);
    const auto *nsg1_385 = buffer.data(nsg1 + 385);
    const auto *nsg1_387 = buffer.data(nsg1 + 387);
    const auto *nsg1_388 = buffer.data(nsg1 + 388);
    const auto *nsg1_389 = buffer.data(nsg1 + 389);
    const auto *nsg1_400 = buffer.data(nsg1 + 400);
    const auto *nsg1_402 = buffer.data(nsg1 + 402);
    const auto *nsg1_403 = buffer.data(nsg1 + 403);
    const auto *nsg1_404 = buffer.data(nsg1 + 404);
    const auto *nsg1_405 = buffer.data(nsg1 + 405);
    const auto *nsg1_406 = buffer.data(nsg1 + 406);
    const auto *nsg1_407 = buffer.data(nsg1 + 407);
    const auto *nsg1_408 = buffer.data(nsg1 + 408);
    const auto *nsg1_409 = buffer.data(nsg1 + 409);
    const auto *nsg1_410 = buffer.data(nsg1 + 410);
    const auto *nsg1_414 = buffer.data(nsg1 + 414);
    const auto *nsg1_415 = buffer.data(nsg1 + 415);
    const auto *nsg1_416 = buffer.data(nsg1 + 416);
    const auto *nsg1_417 = buffer.data(nsg1 + 417);
    const auto *nsg1_418 = buffer.data(nsg1 + 418);
    const auto *nsg1_419 = buffer.data(nsg1 + 419);
    const auto *nsg1_420 = buffer.data(nsg1 + 420);
    const auto *nsg1_422 = buffer.data(nsg1 + 422);
    const auto *nsg1_423 = buffer.data(nsg1 + 423);
    const auto *nsg1_425 = buffer.data(nsg1 + 425);
    const auto *nsg1_426 = buffer.data(nsg1 + 426);
    const auto *nsg1_430 = buffer.data(nsg1 + 430);
    const auto *nsg1_431 = buffer.data(nsg1 + 431);
    const auto *nsg1_432 = buffer.data(nsg1 + 432);
    const auto *nsg1_434 = buffer.data(nsg1 + 434);
    const auto *nsg1_440 = buffer.data(nsg1 + 440);

    const auto *nsh_527 = buffer.data(nsh + 527);
    const auto *nsh_528 = buffer.data(nsh + 528);
    const auto *nsh_530 = buffer.data(nsh + 530);
    const auto *nsh_531 = buffer.data(nsh + 531);
    const auto *nsh_534 = buffer.data(nsh + 534);
    const auto *nsh_535 = buffer.data(nsh + 535);
    const auto *nsh_537 = buffer.data(nsh + 537);
    const auto *nsh_539 = buffer.data(nsh + 539);
    const auto *nsh_540 = buffer.data(nsh + 540);
    const auto *nsh_541 = buffer.data(nsh + 541);
    const auto *nsh_542 = buffer.data(nsh + 542);
    const auto *nsh_543 = buffer.data(nsh + 543);
    const auto *nsh_544 = buffer.data(nsh + 544);
    const auto *nsh_545 = buffer.data(nsh + 545);
    const auto *nsh_546 = buffer.data(nsh + 546);
    const auto *nsh_548 = buffer.data(nsh + 548);
    const auto *nsh_549 = buffer.data(nsh + 549);
    const auto *nsh_551 = buffer.data(nsh + 551);
    const auto *nsh_552 = buffer.data(nsh + 552);
    const auto *nsh_555 = buffer.data(nsh + 555);
    const auto *nsh_561 = buffer.data(nsh + 561);
    const auto *nsh_562 = buffer.data(nsh + 562);
    const auto *nsh_563 = buffer.data(nsh + 563);
    const auto *nsh_564 = buffer.data(nsh + 564);
    const auto *nsh_565 = buffer.data(nsh + 565);
    const auto *nsh_566 = buffer.data(nsh + 566);
    const auto *nsh_567 = buffer.data(nsh + 567);
    const auto *nsh_568 = buffer.data(nsh + 568);
    const auto *nsh_569 = buffer.data(nsh + 569);
    const auto *nsh_570 = buffer.data(nsh + 570);
    const auto *nsh_571 = buffer.data(nsh + 571);
    const auto *nsh_572 = buffer.data(nsh + 572);
    const auto *nsh_573 = buffer.data(nsh + 573);
    const auto *nsh_574 = buffer.data(nsh + 574);
    const auto *nsh_575 = buffer.data(nsh + 575);
    const auto *nsh_576 = buffer.data(nsh + 576);
    const auto *nsh_581 = buffer.data(nsh + 581);
    const auto *nsh_582 = buffer.data(nsh + 582);
    const auto *nsh_583 = buffer.data(nsh + 583);
    const auto *nsh_584 = buffer.data(nsh + 584);
    const auto *nsh_585 = buffer.data(nsh + 585);
    const auto *nsh_586 = buffer.data(nsh + 586);
    const auto *nsh_587 = buffer.data(nsh + 587);
    const auto *nsh_588 = buffer.data(nsh + 588);
    const auto *nsh_589 = buffer.data(nsh + 589);
    const auto *nsh_590 = buffer.data(nsh + 590);
    const auto *nsh_591 = buffer.data(nsh + 591);
    const auto *nsh_593 = buffer.data(nsh + 593);
    const auto *nsh_594 = buffer.data(nsh + 594);
    const auto *nsh_595 = buffer.data(nsh + 595);
    const auto *nsh_597 = buffer.data(nsh + 597);
    const auto *nsh_598 = buffer.data(nsh + 598);
    const auto *nsh_603 = buffer.data(nsh + 603);
    const auto *nsh_604 = buffer.data(nsh + 604);
    const auto *nsh_605 = buffer.data(nsh + 605);
    const auto *nsh_606 = buffer.data(nsh + 606);
    const auto *nsh_607 = buffer.data(nsh + 607);
    const auto *nsh_608 = buffer.data(nsh + 608);
    const auto *nsh_609 = buffer.data(nsh + 609);
    const auto *nsh_611 = buffer.data(nsh + 611);
    const auto *nsh_612 = buffer.data(nsh + 612);
    const auto *nsh_614 = buffer.data(nsh + 614);

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, msh_401, msh_528, msh_530, nsg0_378, \
                         nsg0_380, nsg1_378, nsg1_380, nsh_527, nsh_528, \
                         nsh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_14 * msh_528[k]
                   + f_8 * nsg0_378[k]
                   - f_9 * nsg1_378[k]
                   + f_3 * pc_x[k] * nsh_528[k];

        t_704[k] = f_12 * msh_401[k]
                   + f_3 * pc_y[k] * nsh_527[k];

        t_705[k] = f_14 * msh_530[k]
                   + f_8 * nsg0_380[k]
                   - f_9 * nsg1_380[k]
                   + f_3 * pc_x[k] * nsh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, msh_381, msh_404, msh_531, \
                         nsg0_381, nsg1_381, nsh_528, nsh_530, \
                         nsh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_14 * msh_531[k]
                   + f_6 * nsg0_381[k]
                   - f_7 * nsg1_381[k]
                   + f_3 * pc_x[k] * nsh_531[k];

        t_707[k] = f_14 * msh_381[k]
                   + f_3 * pc_z[k] * nsh_528[k];

        t_708[k] = f_12 * msh_404[k]
                   + f_3 * pc_y[k] * nsh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, msh_384, msh_534, msh_535, nsg0_384, \
                         nsg0_385, nsg1_384, nsg1_385, nsh_531, nsh_534, \
                         nsh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_14 * msh_534[k]
                   + f_6 * nsg0_384[k]
                   - f_7 * nsg1_384[k]
                   + f_3 * pc_x[k] * nsh_534[k];

        t_710[k] = f_14 * msh_535[k]
                   + f_4 * nsg0_385[k]
                   - f_5 * nsg1_385[k]
                   + f_3 * pc_x[k] * nsh_535[k];

        t_711[k] = f_14 * msh_384[k]
                   + f_3 * pc_z[k] * nsh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, msh_408, msh_537, msh_539, nsg0_387, \
                         nsg0_389, nsg1_387, nsg1_389, nsh_534, nsh_537, \
                         nsh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_14 * msh_537[k]
                   + f_4 * nsg0_387[k]
                   - f_5 * nsg1_387[k]
                   + f_3 * pc_x[k] * nsh_537[k];

        t_713[k] = f_12 * msh_408[k]
                   + f_3 * pc_y[k] * nsh_534[k];

        t_714[k] = f_14 * msh_539[k]
                   + f_4 * nsg0_389[k]
                   - f_5 * nsg1_389[k]
                   + f_3 * pc_x[k] * nsh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, msh_540, msh_541, msh_542, \
                         msh_543, msh_544, nsh_540, nsh_541, nsh_542, nsh_543, \
                         nsh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_14 * msh_540[k]
                   + f_3 * pc_x[k] * nsh_540[k];

        t_716[k] = f_14 * msh_541[k]
                   + f_3 * pc_x[k] * nsh_541[k];

        t_717[k] = f_14 * msh_542[k]
                   + f_3 * pc_x[k] * nsh_542[k];

        t_718[k] = f_14 * msh_543[k]
                   + f_3 * pc_x[k] * nsh_543[k];

        t_719[k] = f_14 * msh_544[k]
                   + f_3 * pc_x[k] * nsh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, msh_393, msh_414, msh_545, \
                         nsg0_385, nsg1_385, nsh_540, nsh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_14 * msh_545[k]
                   + f_3 * pc_x[k] * nsh_545[k];

        t_721[k] = f_12 * msh_414[k]
                   + f_1 * nsg0_385[k]
                   - f_2 * nsg1_385[k]
                   + f_3 * pc_y[k] * nsh_540[k];

        t_722[k] = f_14 * msh_393[k]
                   + f_3 * pc_z[k] * nsh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, msh_416, msh_417, msh_418, nsg0_387, \
                         nsg0_388, nsg0_389, nsg1_387, nsg1_388, nsg1_389, nsh_542, nsh_543, \
                         nsh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * msh_416[k]
                   + f_8 * nsg0_387[k]
                   - f_9 * nsg1_387[k]
                   + f_3 * pc_y[k] * nsh_542[k];

        t_724[k] = f_12 * msh_417[k]
                   + f_6 * nsg0_388[k]
                   - f_7 * nsg1_388[k]
                   + f_3 * pc_y[k] * nsh_543[k];

        t_725[k] = f_12 * msh_418[k]
                   + f_4 * nsg0_389[k]
                   - f_5 * nsg1_389[k]
                   + f_3 * pc_y[k] * nsh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_y, pc_y, pc_z, msi0_560, msh_398, \
                         msh_419, msh_420, msi1_560, nsg0_389, nsg1_389, nsh_545, \
                         nsh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * msh_419[k]
                   + f_3 * pc_y[k] * nsh_545[k];

        t_727[k] = f_14 * msh_398[k]
                   + f_1 * nsg0_389[k]
                   - f_2 * nsg1_389[k]
                   + f_3 * pc_z[k] * nsh_545[k];

        t_728[k] = pa_y[k] * msi0_560[k]
                   - f_10 * pc_y[k] * msi1_560[k];

        t_729[k] = f_11 * msh_420[k]
                   + f_3 * pc_y[k] * nsh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_y, pc_y, pc_z, msi0_563, msi0_565, \
                         msh_399, msh_421, msh_422, msi1_563, msi1_565, nsh_546, \
                         nsh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_21 * msh_399[k]
                   + f_3 * pc_z[k] * nsh_546[k];

        t_731[k] = pa_y[k] * msi0_563[k]
                   + f_12 * msh_421[k]
                   - f_10 * pc_y[k] * msi1_563[k];

        t_732[k] = f_11 * msh_422[k]
                   + f_3 * pc_y[k] * nsh_548[k];

        t_733[k] = pa_y[k] * msi0_565[k]
                   - f_10 * pc_y[k] * msi1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pa_y, pc_y, pc_z, msi0_566, msi0_569, \
                         msh_402, msh_423, msh_425, msi1_566, msi1_569, nsh_549, \
                         nsh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pa_y[k] * msi0_566[k]
                   + f_13 * msh_423[k]
                   - f_10 * pc_y[k] * msi1_566[k];

        t_735[k] = f_21 * msh_402[k]
                   + f_3 * pc_z[k] * nsh_549[k];

        t_736[k] = f_11 * msh_425[k]
                   + f_3 * pc_y[k] * nsh_551[k];

        t_737[k] = pa_y[k] * msi0_569[k]
                   - f_10 * pc_y[k] * msi1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pa_y, pc_y, pc_z, msi0_570, msi0_572, msh_405, \
                         msh_426, msh_428, msi1_570, msi1_572, \
                         nsh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pa_y[k] * msi0_570[k]
                   + f_14 * msh_426[k]
                   - f_10 * pc_y[k] * msi1_570[k];

        t_739[k] = f_21 * msh_405[k]
                   + f_3 * pc_z[k] * nsh_552[k];

        t_740[k] = pa_y[k] * msi0_572[k]
                   + f_12 * msh_428[k]
                   - f_10 * pc_y[k] * msi1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pa_y, pc_x, pc_y, msi0_574, msh_429, \
                         msh_561, msh_562, msi1_574, nsh_555, nsh_561, \
                         nsh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * msh_429[k]
                   + f_3 * pc_y[k] * nsh_555[k];

        t_742[k] = pa_y[k] * msi0_574[k]
                   - f_10 * pc_y[k] * msi1_574[k];

        t_743[k] = f_14 * msh_561[k]
                   + f_3 * pc_x[k] * nsh_561[k];

        t_744[k] = f_14 * msh_562[k]
                   + f_3 * pc_x[k] * nsh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, msh_563, msh_564, msh_565, msh_566, \
                         nsh_563, nsh_564, nsh_565, nsh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_14 * msh_563[k]
                   + f_3 * pc_x[k] * nsh_563[k];

        t_746[k] = f_14 * msh_564[k]
                   + f_3 * pc_x[k] * nsh_564[k];

        t_747[k] = f_14 * msh_565[k]
                   + f_3 * pc_x[k] * nsh_565[k];

        t_748[k] = f_14 * msh_566[k]
                   + f_3 * pc_x[k] * nsh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, msh_414, msh_435, msh_437, nsg0_400, \
                         nsg0_402, nsg1_400, nsg1_402, nsh_561, \
                         nsh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * msh_435[k]
                   + f_1 * nsg0_400[k]
                   - f_2 * nsg1_400[k]
                   + f_3 * pc_y[k] * nsh_561[k];

        t_750[k] = f_21 * msh_414[k]
                   + f_3 * pc_z[k] * nsh_561[k];

        t_751[k] = f_11 * msh_437[k]
                   + f_8 * nsg0_402[k]
                   - f_9 * nsg1_402[k]
                   + f_3 * pc_y[k] * nsh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, msh_438, msh_439, msh_440, nsg0_403, \
                         nsg0_404, nsg1_403, nsg1_404, nsh_564, nsh_565, \
                         nsh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * msh_438[k]
                   + f_6 * nsg0_403[k]
                   - f_7 * nsg1_403[k]
                   + f_3 * pc_y[k] * nsh_564[k];

        t_753[k] = f_11 * msh_439[k]
                   + f_4 * nsg0_404[k]
                   - f_5 * nsg1_404[k]
                   + f_3 * pc_y[k] * nsh_565[k];

        t_754[k] = f_11 * msh_440[k]
                   + f_3 * pc_y[k] * nsh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pa_y, pc_x, pc_y, pc_z, msi0_587, \
                         msh_420, msh_567, msi1_587, nsg0_405, nsg1_405, \
                         nsh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pa_y[k] * msi0_587[k]
                   - f_10 * pc_y[k] * msi1_587[k];

        t_756[k] = f_14 * msh_567[k]
                   + f_1 * nsg0_405[k]
                   - f_2 * nsg1_405[k]
                   + f_3 * pc_x[k] * nsh_567[k];

        t_757[k] = f_3 * pc_y[k] * nsh_567[k];

        t_758[k] = f_20 * msh_420[k]
                   + f_3 * pc_z[k] * nsh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, msh_572, nsg0_405, nsg0_410, \
                         nsg1_405, nsg1_410, nsh_568, nsh_569, \
                         nsh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_4 * nsg0_405[k]
                   - f_5 * nsg1_405[k]
                   + f_3 * pc_y[k] * nsh_568[k];

        t_760[k] = f_3 * pc_y[k] * nsh_569[k];

        t_761[k] = f_14 * msh_572[k]
                   + f_8 * nsg0_410[k]
                   - f_9 * nsg1_410[k]
                   + f_3 * pc_x[k] * nsh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_y, nsg0_406, nsg0_407, nsg1_406, nsg1_407, \
                         nsh_570, nsh_571, nsh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_6 * nsg0_406[k]
                   - f_7 * nsg1_406[k]
                   + f_3 * pc_y[k] * nsh_570[k];

        t_763[k] = f_4 * nsg0_407[k]
                   - f_5 * nsg1_407[k]
                   + f_3 * pc_y[k] * nsh_571[k];

        t_764[k] = f_3 * pc_y[k] * nsh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_y, msh_576, nsg0_408, nsg0_409, \
                         nsg0_414, nsg1_408, nsg1_409, nsg1_414, nsh_573, nsh_574, \
                         nsh_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_14 * msh_576[k]
                   + f_6 * nsg0_414[k]
                   - f_7 * nsg1_414[k]
                   + f_3 * pc_x[k] * nsh_576[k];

        t_766[k] = f_8 * nsg0_408[k]
                   - f_9 * nsg1_408[k]
                   + f_3 * pc_y[k] * nsh_573[k];

        t_767[k] = f_6 * nsg0_409[k]
                   - f_7 * nsg1_409[k]
                   + f_3 * pc_y[k] * nsh_574[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, msh_581, msh_582, nsg0_410, \
                         nsg0_419, nsg1_410, nsg1_419, nsh_575, nsh_576, nsh_581, \
                         nsh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * nsg0_410[k]
                   - f_5 * nsg1_410[k]
                   + f_3 * pc_y[k] * nsh_575[k];

        t_769[k] = f_3 * pc_y[k] * nsh_576[k];

        t_770[k] = f_14 * msh_581[k]
                   + f_4 * nsg0_419[k]
                   - f_5 * nsg1_419[k]
                   + f_3 * pc_x[k] * nsh_581[k];

        t_771[k] = f_14 * msh_582[k]
                   + f_3 * pc_x[k] * nsh_582[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pc_x, pc_y, msh_583, msh_584, \
                         msh_585, msh_587, nsh_581, nsh_583, nsh_584, nsh_585, \
                         nsh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_14 * msh_583[k]
                   + f_3 * pc_x[k] * nsh_583[k];

        t_773[k] = f_14 * msh_584[k]
                   + f_3 * pc_x[k] * nsh_584[k];

        t_774[k] = f_14 * msh_585[k]
                   + f_3 * pc_x[k] * nsh_585[k];

        t_775[k] = f_3 * pc_y[k] * nsh_581[k];

        t_776[k] = f_14 * msh_587[k]
                   + f_3 * pc_x[k] * nsh_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, nsg0_415, nsg0_416, nsg0_417, nsg1_415, \
                         nsg1_416, nsg1_417, nsh_582, nsh_583, \
                         nsh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_1 * nsg0_415[k]
                   - f_2 * nsg1_415[k]
                   + f_3 * pc_y[k] * nsh_582[k];

        t_778[k] = f_16 * nsg0_416[k]
                   - f_17 * nsg1_416[k]
                   + f_3 * pc_y[k] * nsh_583[k];

        t_779[k] = f_8 * nsg0_417[k]
                   - f_9 * nsg1_417[k]
                   + f_3 * pc_y[k] * nsh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, msh_440, nsg0_418, nsg0_419, \
                         nsg1_418, nsg1_419, nsh_585, nsh_586, \
                         nsh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * nsg0_418[k]
                   - f_7 * nsg1_418[k]
                   + f_3 * pc_y[k] * nsh_585[k];

        t_781[k] = f_4 * nsg0_419[k]
                   - f_5 * nsg1_419[k]
                   + f_3 * pc_y[k] * nsh_586[k];

        t_782[k] = f_3 * pc_y[k] * nsh_587[k];

        t_783[k] = f_20 * msh_440[k]
                   + f_1 * nsg0_419[k]
                   - f_2 * nsg1_419[k]
                   + f_3 * pc_z[k] * nsh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pc_x, pc_y, pc_z, msh_441, msh_588, \
                         msh_591, nsg0_420, nsg0_423, nsg1_420, nsg1_423, nsh_588, \
                         nsh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_13 * msh_588[k]
                   + f_1 * nsg0_420[k]
                   - f_2 * nsg1_420[k]
                   + f_3 * pc_x[k] * nsh_588[k];

        t_785[k] = f_19 * msh_441[k]
                   + f_3 * pc_y[k] * nsh_588[k];

        t_786[k] = f_3 * pc_z[k] * nsh_588[k];

        t_787[k] = f_13 * msh_591[k]
                   + f_8 * nsg0_423[k]
                   - f_9 * nsg1_423[k]
                   + f_3 * pc_x[k] * nsh_591[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, t_791, pc_x, pc_z, msh_594, nsg0_420, nsg0_426, \
                         nsg1_420, nsg1_426, nsh_589, nsh_590, nsh_591, \
                         nsh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_3 * pc_z[k] * nsh_589[k];

        t_789[k] = f_4 * nsg0_420[k]
                   - f_5 * nsg1_420[k]
                   + f_3 * pc_z[k] * nsh_590[k];

        t_790[k] = f_13 * msh_594[k]
                   + f_6 * nsg0_426[k]
                   - f_7 * nsg1_426[k]
                   + f_3 * pc_x[k] * nsh_594[k];

        t_791[k] = f_3 * pc_z[k] * nsh_591[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pc_x, pc_y, pc_z, msh_446, msh_598, \
                         nsg0_422, nsg0_430, nsg1_422, nsg1_430, nsh_593, nsh_594, \
                         nsh_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_19 * msh_446[k]
                   + f_3 * pc_y[k] * nsh_593[k];

        t_793[k] = f_6 * nsg0_422[k]
                   - f_7 * nsg1_422[k]
                   + f_3 * pc_z[k] * nsh_593[k];

        t_794[k] = f_13 * msh_598[k]
                   + f_4 * nsg0_430[k]
                   - f_5 * nsg1_430[k]
                   + f_3 * pc_x[k] * nsh_598[k];

        t_795[k] = f_3 * pc_z[k] * nsh_594[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pc_x, pc_y, pc_z, msh_450, msh_603, \
                         nsg0_423, nsg0_425, nsg1_423, nsg1_425, nsh_595, nsh_597, \
                         nsh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_4 * nsg0_423[k]
                   - f_5 * nsg1_423[k]
                   + f_3 * pc_z[k] * nsh_595[k];

        t_797[k] = f_19 * msh_450[k]
                   + f_3 * pc_y[k] * nsh_597[k];

        t_798[k] = f_8 * nsg0_425[k]
                   - f_9 * nsg1_425[k]
                   + f_3 * pc_z[k] * nsh_597[k];

        t_799[k] = f_13 * msh_603[k]
                   + f_3 * pc_x[k] * nsh_603[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, pc_x, pc_z, msh_605, msh_606, \
                         msh_607, msh_608, nsh_598, nsh_605, nsh_606, nsh_607, \
                         nsh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_3 * pc_z[k] * nsh_598[k];

        t_801[k] = f_13 * msh_605[k]
                   + f_3 * pc_x[k] * nsh_605[k];

        t_802[k] = f_13 * msh_606[k]
                   + f_3 * pc_x[k] * nsh_606[k];

        t_803[k] = f_13 * msh_607[k]
                   + f_3 * pc_x[k] * nsh_607[k];

        t_804[k] = f_13 * msh_608[k]
                   + f_3 * pc_x[k] * nsh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, pc_y, pc_z, msh_456, nsg0_430, nsg0_431, \
                         nsg1_430, nsg1_431, nsh_603, nsh_604, \
                         nsh_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_19 * msh_456[k]
                   + f_1 * nsg0_430[k]
                   - f_2 * nsg1_430[k]
                   + f_3 * pc_y[k] * nsh_603[k];

        t_806[k] = f_3 * pc_z[k] * nsh_603[k];

        t_807[k] = f_4 * nsg0_430[k]
                   - f_5 * nsg1_430[k]
                   + f_3 * pc_z[k] * nsh_604[k];

        t_808[k] = f_6 * nsg0_431[k]
                   - f_7 * nsg1_431[k]
                   + f_3 * pc_z[k] * nsh_605[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, t_812, pa_z, pc_y, pc_z, msi0_588, msh_461, \
                         msi1_588, nsg0_432, nsg0_434, nsg1_432, nsg1_434, nsh_606, \
                         nsh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_8 * nsg0_432[k]
                   - f_9 * nsg1_432[k]
                   + f_3 * pc_z[k] * nsh_606[k];

        t_810[k] = f_19 * msh_461[k]
                   + f_3 * pc_y[k] * nsh_608[k];

        t_811[k] = f_1 * nsg0_434[k]
                   - f_2 * nsg1_434[k]
                   + f_3 * pc_z[k] * nsh_608[k];

        t_812[k] = pa_z[k] * msi0_588[k]
                   - f_10 * pc_z[k] * msi1_588[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, t_816, pa_z, pc_y, pc_z, msi0_591, msh_441, \
                         msh_462, msh_464, msi1_591, nsh_609, nsh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_20 * msh_462[k]
                   + f_3 * pc_y[k] * nsh_609[k];

        t_814[k] = f_11 * msh_441[k]
                   + f_3 * pc_z[k] * nsh_609[k];

        t_815[k] = pa_z[k] * msi0_591[k]
                   - f_10 * pc_z[k] * msi1_591[k];

        t_816[k] = f_20 * msh_464[k]
                   + f_3 * pc_y[k] * nsh_611[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, pa_z, pc_x, pc_z, msi0_594, msh_444, msh_614, \
                         msi1_594, nsg0_440, nsg1_440, nsh_612, \
                         nsh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_13 * msh_614[k]
                   + f_8 * nsg0_440[k]
                   - f_9 * nsg1_440[k]
                   + f_3 * pc_x[k] * nsh_614[k];

        t_818[k] = pa_z[k] * msi0_594[k]
                   - f_10 * pc_z[k] * msi1_594[k];

        t_819[k] = f_11 * msh_444[k]
                   + f_3 * pc_z[k] * nsh_612[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_598 = buffer.data(msi0 + 598);
    const auto *msi0_600 = buffer.data(msi0 + 600);
    const auto *msi0_609 = buffer.data(msi0 + 609);

    const auto *msh_447 = buffer.data(msh + 447);
    const auto *msh_448 = buffer.data(msh + 448);
    const auto *msh_456 = buffer.data(msh + 456);
    const auto *msh_461 = buffer.data(msh + 461);
    const auto *msh_462 = buffer.data(msh + 462);
    const auto *msh_465 = buffer.data(msh + 465);
    const auto *msh_467 = buffer.data(msh + 467);
    const auto *msh_468 = buffer.data(msh + 468);
    const auto *msh_471 = buffer.data(msh + 471);
    const auto *msh_477 = buffer.data(msh + 477);
    const auto *msh_479 = buffer.data(msh + 479);
    const auto *msh_480 = buffer.data(msh + 480);
    const auto *msh_481 = buffer.data(msh + 481);
    const auto *msh_482 = buffer.data(msh + 482);
    const auto *msh_483 = buffer.data(msh + 483);
    const auto *msh_485 = buffer.data(msh + 485);
    const auto *msh_486 = buffer.data(msh + 486);
    const auto *msh_488 = buffer.data(msh + 488);
    const auto *msh_489 = buffer.data(msh + 489);
    const auto *msh_492 = buffer.data(msh + 492);
    const auto *msh_498 = buffer.data(msh + 498);
    const auto *msh_500 = buffer.data(msh + 500);
    const auto *msh_501 = buffer.data(msh + 501);
    const auto *msh_502 = buffer.data(msh + 502);
    const auto *msh_503 = buffer.data(msh + 503);
    const auto *msh_504 = buffer.data(msh + 504);
    const auto *msh_506 = buffer.data(msh + 506);
    const auto *msh_507 = buffer.data(msh + 507);
    const auto *msh_509 = buffer.data(msh + 509);
    const auto *msh_510 = buffer.data(msh + 510);
    const auto *msh_513 = buffer.data(msh + 513);
    const auto *msh_519 = buffer.data(msh + 519);
    const auto *msh_521 = buffer.data(msh + 521);
    const auto *msh_522 = buffer.data(msh + 522);
    const auto *msh_523 = buffer.data(msh + 523);
    const auto *msh_524 = buffer.data(msh + 524);
    const auto *msh_525 = buffer.data(msh + 525);
    const auto *msh_527 = buffer.data(msh + 527);
    const auto *msh_530 = buffer.data(msh + 530);
    const auto *msh_534 = buffer.data(msh + 534);
    const auto *msh_540 = buffer.data(msh + 540);
    const auto *msh_542 = buffer.data(msh + 542);
    const auto *msh_543 = buffer.data(msh + 543);
    const auto *msh_544 = buffer.data(msh + 544);
    const auto *msh_545 = buffer.data(msh + 545);
    const auto *msh_618 = buffer.data(msh + 618);
    const auto *msh_623 = buffer.data(msh + 623);
    const auto *msh_624 = buffer.data(msh + 624);
    const auto *msh_625 = buffer.data(msh + 625);
    const auto *msh_626 = buffer.data(msh + 626);
    const auto *msh_627 = buffer.data(msh + 627);
    const auto *msh_628 = buffer.data(msh + 628);
    const auto *msh_629 = buffer.data(msh + 629);
    const auto *msh_630 = buffer.data(msh + 630);
    const auto *msh_633 = buffer.data(msh + 633);
    const auto *msh_635 = buffer.data(msh + 635);
    const auto *msh_636 = buffer.data(msh + 636);
    const auto *msh_639 = buffer.data(msh + 639);
    const auto *msh_640 = buffer.data(msh + 640);
    const auto *msh_642 = buffer.data(msh + 642);
    const auto *msh_644 = buffer.data(msh + 644);
    const auto *msh_645 = buffer.data(msh + 645);
    const auto *msh_646 = buffer.data(msh + 646);
    const auto *msh_647 = buffer.data(msh + 647);
    const auto *msh_648 = buffer.data(msh + 648);
    const auto *msh_649 = buffer.data(msh + 649);
    const auto *msh_650 = buffer.data(msh + 650);
    const auto *msh_651 = buffer.data(msh + 651);
    const auto *msh_654 = buffer.data(msh + 654);
    const auto *msh_656 = buffer.data(msh + 656);
    const auto *msh_657 = buffer.data(msh + 657);
    const auto *msh_660 = buffer.data(msh + 660);
    const auto *msh_661 = buffer.data(msh + 661);
    const auto *msh_663 = buffer.data(msh + 663);
    const auto *msh_665 = buffer.data(msh + 665);
    const auto *msh_666 = buffer.data(msh + 666);
    const auto *msh_667 = buffer.data(msh + 667);
    const auto *msh_668 = buffer.data(msh + 668);
    const auto *msh_669 = buffer.data(msh + 669);
    const auto *msh_670 = buffer.data(msh + 670);
    const auto *msh_671 = buffer.data(msh + 671);
    const auto *msh_672 = buffer.data(msh + 672);
    const auto *msh_675 = buffer.data(msh + 675);
    const auto *msh_677 = buffer.data(msh + 677);
    const auto *msh_678 = buffer.data(msh + 678);
    const auto *msh_681 = buffer.data(msh + 681);
    const auto *msh_682 = buffer.data(msh + 682);
    const auto *msh_684 = buffer.data(msh + 684);
    const auto *msh_686 = buffer.data(msh + 686);
    const auto *msh_687 = buffer.data(msh + 687);
    const auto *msh_688 = buffer.data(msh + 688);
    const auto *msh_689 = buffer.data(msh + 689);
    const auto *msh_690 = buffer.data(msh + 690);
    const auto *msh_691 = buffer.data(msh + 691);
    const auto *msh_692 = buffer.data(msh + 692);
    const auto *msh_693 = buffer.data(msh + 693);

    const auto *msi1_598 = buffer.data(msi1 + 598);
    const auto *msi1_600 = buffer.data(msi1 + 600);
    const auto *msi1_609 = buffer.data(msi1 + 609);

    const auto *nsg0_444 = buffer.data(nsg0 + 444);
    const auto *nsg0_447 = buffer.data(nsg0 + 447);
    const auto *nsg0_448 = buffer.data(nsg0 + 448);
    const auto *nsg0_449 = buffer.data(nsg0 + 449);
    const auto *nsg0_450 = buffer.data(nsg0 + 450);
    const auto *nsg0_453 = buffer.data(nsg0 + 453);
    const auto *nsg0_455 = buffer.data(nsg0 + 455);
    const auto *nsg0_456 = buffer.data(nsg0 + 456);
    const auto *nsg0_459 = buffer.data(nsg0 + 459);
    const auto *nsg0_460 = buffer.data(nsg0 + 460);
    const auto *nsg0_462 = buffer.data(nsg0 + 462);
    const auto *nsg0_463 = buffer.data(nsg0 + 463);
    const auto *nsg0_464 = buffer.data(nsg0 + 464);
    const auto *nsg0_465 = buffer.data(nsg0 + 465);
    const auto *nsg0_468 = buffer.data(nsg0 + 468);
    const auto *nsg0_470 = buffer.data(nsg0 + 470);
    const auto *nsg0_471 = buffer.data(nsg0 + 471);
    const auto *nsg0_474 = buffer.data(nsg0 + 474);
    const auto *nsg0_475 = buffer.data(nsg0 + 475);
    const auto *nsg0_477 = buffer.data(nsg0 + 477);
    const auto *nsg0_478 = buffer.data(nsg0 + 478);
    const auto *nsg0_479 = buffer.data(nsg0 + 479);
    const auto *nsg0_480 = buffer.data(nsg0 + 480);
    const auto *nsg0_483 = buffer.data(nsg0 + 483);
    const auto *nsg0_485 = buffer.data(nsg0 + 485);
    const auto *nsg0_486 = buffer.data(nsg0 + 486);
    const auto *nsg0_489 = buffer.data(nsg0 + 489);
    const auto *nsg0_490 = buffer.data(nsg0 + 490);
    const auto *nsg0_492 = buffer.data(nsg0 + 492);
    const auto *nsg0_493 = buffer.data(nsg0 + 493);
    const auto *nsg0_494 = buffer.data(nsg0 + 494);
    const auto *nsg0_495 = buffer.data(nsg0 + 495);

    const auto *nsg1_444 = buffer.data(nsg1 + 444);
    const auto *nsg1_447 = buffer.data(nsg1 + 447);
    const auto *nsg1_448 = buffer.data(nsg1 + 448);
    const auto *nsg1_449 = buffer.data(nsg1 + 449);
    const auto *nsg1_450 = buffer.data(nsg1 + 450);
    const auto *nsg1_453 = buffer.data(nsg1 + 453);
    const auto *nsg1_455 = buffer.data(nsg1 + 455);
    const auto *nsg1_456 = buffer.data(nsg1 + 456);
    const auto *nsg1_459 = buffer.data(nsg1 + 459);
    const auto *nsg1_460 = buffer.data(nsg1 + 460);
    const auto *nsg1_462 = buffer.data(nsg1 + 462);
    const auto *nsg1_463 = buffer.data(nsg1 + 463);
    const auto *nsg1_464 = buffer.data(nsg1 + 464);
    const auto *nsg1_465 = buffer.data(nsg1 + 465);
    const auto *nsg1_468 = buffer.data(nsg1 + 468);
    const auto *nsg1_470 = buffer.data(nsg1 + 470);
    const auto *nsg1_471 = buffer.data(nsg1 + 471);
    const auto *nsg1_474 = buffer.data(nsg1 + 474);
    const auto *nsg1_475 = buffer.data(nsg1 + 475);
    const auto *nsg1_477 = buffer.data(nsg1 + 477);
    const auto *nsg1_478 = buffer.data(nsg1 + 478);
    const auto *nsg1_479 = buffer.data(nsg1 + 479);
    const auto *nsg1_480 = buffer.data(nsg1 + 480);
    const auto *nsg1_483 = buffer.data(nsg1 + 483);
    const auto *nsg1_485 = buffer.data(nsg1 + 485);
    const auto *nsg1_486 = buffer.data(nsg1 + 486);
    const auto *nsg1_489 = buffer.data(nsg1 + 489);
    const auto *nsg1_490 = buffer.data(nsg1 + 490);
    const auto *nsg1_492 = buffer.data(nsg1 + 492);
    const auto *nsg1_493 = buffer.data(nsg1 + 493);
    const auto *nsg1_494 = buffer.data(nsg1 + 494);
    const auto *nsg1_495 = buffer.data(nsg1 + 495);

    const auto *nsh_614 = buffer.data(nsh + 614);
    const auto *nsh_615 = buffer.data(nsh + 615);
    const auto *nsh_618 = buffer.data(nsh + 618);
    const auto *nsh_623 = buffer.data(nsh + 623);
    const auto *nsh_624 = buffer.data(nsh + 624);
    const auto *nsh_625 = buffer.data(nsh + 625);
    const auto *nsh_626 = buffer.data(nsh + 626);
    const auto *nsh_627 = buffer.data(nsh + 627);
    const auto *nsh_628 = buffer.data(nsh + 628);
    const auto *nsh_629 = buffer.data(nsh + 629);
    const auto *nsh_630 = buffer.data(nsh + 630);
    const auto *nsh_632 = buffer.data(nsh + 632);
    const auto *nsh_633 = buffer.data(nsh + 633);
    const auto *nsh_635 = buffer.data(nsh + 635);
    const auto *nsh_636 = buffer.data(nsh + 636);
    const auto *nsh_639 = buffer.data(nsh + 639);
    const auto *nsh_640 = buffer.data(nsh + 640);
    const auto *nsh_642 = buffer.data(nsh + 642);
    const auto *nsh_644 = buffer.data(nsh + 644);
    const auto *nsh_645 = buffer.data(nsh + 645);
    const auto *nsh_646 = buffer.data(nsh + 646);
    const auto *nsh_647 = buffer.data(nsh + 647);
    const auto *nsh_648 = buffer.data(nsh + 648);
    const auto *nsh_649 = buffer.data(nsh + 649);
    const auto *nsh_650 = buffer.data(nsh + 650);
    const auto *nsh_651 = buffer.data(nsh + 651);
    const auto *nsh_653 = buffer.data(nsh + 653);
    const auto *nsh_654 = buffer.data(nsh + 654);
    const auto *nsh_656 = buffer.data(nsh + 656);
    const auto *nsh_657 = buffer.data(nsh + 657);
    const auto *nsh_660 = buffer.data(nsh + 660);
    const auto *nsh_661 = buffer.data(nsh + 661);
    const auto *nsh_663 = buffer.data(nsh + 663);
    const auto *nsh_665 = buffer.data(nsh + 665);
    const auto *nsh_666 = buffer.data(nsh + 666);
    const auto *nsh_667 = buffer.data(nsh + 667);
    const auto *nsh_668 = buffer.data(nsh + 668);
    const auto *nsh_669 = buffer.data(nsh + 669);
    const auto *nsh_670 = buffer.data(nsh + 670);
    const auto *nsh_671 = buffer.data(nsh + 671);
    const auto *nsh_672 = buffer.data(nsh + 672);
    const auto *nsh_674 = buffer.data(nsh + 674);
    const auto *nsh_675 = buffer.data(nsh + 675);
    const auto *nsh_677 = buffer.data(nsh + 677);
    const auto *nsh_678 = buffer.data(nsh + 678);
    const auto *nsh_681 = buffer.data(nsh + 681);
    const auto *nsh_682 = buffer.data(nsh + 682);
    const auto *nsh_684 = buffer.data(nsh + 684);
    const auto *nsh_686 = buffer.data(nsh + 686);
    const auto *nsh_687 = buffer.data(nsh + 687);
    const auto *nsh_688 = buffer.data(nsh + 688);
    const auto *nsh_689 = buffer.data(nsh + 689);
    const auto *nsh_690 = buffer.data(nsh + 690);
    const auto *nsh_691 = buffer.data(nsh + 691);
    const auto *nsh_692 = buffer.data(nsh + 692);
    const auto *nsh_693 = buffer.data(nsh + 693);

#pragma omp simd aligned(t_820, t_821, t_822, pa_z, pc_x, pc_y, pc_z, msi0_598, msh_467, \
                         msh_618, msi1_598, nsg0_444, nsg1_444, nsh_614, \
                         nsh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_820[k] = f_20 * msh_467[k]
                   + f_3 * pc_y[k] * nsh_614[k];

        t_821[k] = f_13 * msh_618[k]
                   + f_6 * nsg0_444[k]
                   - f_7 * nsg1_444[k]
                   + f_3 * pc_x[k] * nsh_618[k];

        t_822[k] = pa_z[k] * msi0_598[k]
                   - f_10 * pc_z[k] * msi1_598[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pa_z, pc_y, pc_z, msi0_600, msh_447, msh_448, \
                         msh_471, msi1_600, nsh_615, nsh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_11 * msh_447[k]
                   + f_3 * pc_z[k] * nsh_615[k];

        t_824[k] = pa_z[k] * msi0_600[k]
                   + f_12 * msh_448[k]
                   - f_10 * pc_z[k] * msi1_600[k];

        t_825[k] = f_20 * msh_471[k]
                   + f_3 * pc_y[k] * nsh_618[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, msh_623, msh_624, msh_625, msh_626, \
                         nsg0_449, nsg1_449, nsh_623, nsh_624, nsh_625, \
                         nsh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_13 * msh_623[k]
                   + f_4 * nsg0_449[k]
                   - f_5 * nsg1_449[k]
                   + f_3 * pc_x[k] * nsh_623[k];

        t_827[k] = f_13 * msh_624[k]
                   + f_3 * pc_x[k] * nsh_624[k];

        t_828[k] = f_13 * msh_625[k]
                   + f_3 * pc_x[k] * nsh_625[k];

        t_829[k] = f_13 * msh_626[k]
                   + f_3 * pc_x[k] * nsh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pa_z, pc_x, pc_z, msi0_609, msh_627, \
                         msh_628, msh_629, msi1_609, nsh_627, nsh_628, \
                         nsh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_13 * msh_627[k]
                   + f_3 * pc_x[k] * nsh_627[k];

        t_831[k] = f_13 * msh_628[k]
                   + f_3 * pc_x[k] * nsh_628[k];

        t_832[k] = f_13 * msh_629[k]
                   + f_3 * pc_x[k] * nsh_629[k];

        t_833[k] = pa_z[k] * msi0_609[k]
                   - f_10 * pc_z[k] * msi1_609[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_y, pc_z, msh_456, msh_479, msh_480, nsg0_447, \
                         nsg0_448, nsg1_447, nsg1_448, nsh_624, nsh_626, \
                         nsh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_11 * msh_456[k]
                   + f_3 * pc_z[k] * nsh_624[k];

        t_835[k] = f_20 * msh_479[k]
                   + f_8 * nsg0_447[k]
                   - f_9 * nsg1_447[k]
                   + f_3 * pc_y[k] * nsh_626[k];

        t_836[k] = f_20 * msh_480[k]
                   + f_6 * nsg0_448[k]
                   - f_7 * nsg1_448[k]
                   + f_3 * pc_y[k] * nsh_627[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, msh_461, msh_481, msh_482, nsg0_449, \
                         nsg1_449, nsh_628, nsh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_20 * msh_481[k]
                   + f_4 * nsg0_449[k]
                   - f_5 * nsg1_449[k]
                   + f_3 * pc_y[k] * nsh_628[k];

        t_838[k] = f_20 * msh_482[k]
                   + f_3 * pc_y[k] * nsh_629[k];

        t_839[k] = f_11 * msh_461[k]
                   + f_1 * nsg0_449[k]
                   - f_2 * nsg1_449[k]
                   + f_3 * pc_z[k] * nsh_629[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, pc_y, pc_z, msh_462, msh_483, msh_630, \
                         nsg0_450, nsg1_450, nsh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_13 * msh_630[k]
                   + f_1 * nsg0_450[k]
                   - f_2 * nsg1_450[k]
                   + f_3 * pc_x[k] * nsh_630[k];

        t_841[k] = f_21 * msh_483[k]
                   + f_3 * pc_y[k] * nsh_630[k];

        t_842[k] = f_12 * msh_462[k]
                   + f_3 * pc_z[k] * nsh_630[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, pc_y, msh_485, msh_633, msh_635, nsg0_453, \
                         nsg0_455, nsg1_453, nsg1_455, nsh_632, nsh_633, \
                         nsh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_13 * msh_633[k]
                   + f_8 * nsg0_453[k]
                   - f_9 * nsg1_453[k]
                   + f_3 * pc_x[k] * nsh_633[k];

        t_844[k] = f_21 * msh_485[k]
                   + f_3 * pc_y[k] * nsh_632[k];

        t_845[k] = f_13 * msh_635[k]
                   + f_8 * nsg0_455[k]
                   - f_9 * nsg1_455[k]
                   + f_3 * pc_x[k] * nsh_635[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pc_x, pc_y, pc_z, msh_465, msh_488, msh_636, \
                         nsg0_456, nsg1_456, nsh_633, nsh_635, \
                         nsh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_13 * msh_636[k]
                   + f_6 * nsg0_456[k]
                   - f_7 * nsg1_456[k]
                   + f_3 * pc_x[k] * nsh_636[k];

        t_847[k] = f_12 * msh_465[k]
                   + f_3 * pc_z[k] * nsh_633[k];

        t_848[k] = f_21 * msh_488[k]
                   + f_3 * pc_y[k] * nsh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, msh_468, msh_639, msh_640, nsg0_459, \
                         nsg0_460, nsg1_459, nsg1_460, nsh_636, nsh_639, \
                         nsh_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_13 * msh_639[k]
                   + f_6 * nsg0_459[k]
                   - f_7 * nsg1_459[k]
                   + f_3 * pc_x[k] * nsh_639[k];

        t_850[k] = f_13 * msh_640[k]
                   + f_4 * nsg0_460[k]
                   - f_5 * nsg1_460[k]
                   + f_3 * pc_x[k] * nsh_640[k];

        t_851[k] = f_12 * msh_468[k]
                   + f_3 * pc_z[k] * nsh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, pc_x, pc_y, msh_492, msh_642, msh_644, nsg0_462, \
                         nsg0_464, nsg1_462, nsg1_464, nsh_639, nsh_642, \
                         nsh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_13 * msh_642[k]
                   + f_4 * nsg0_462[k]
                   - f_5 * nsg1_462[k]
                   + f_3 * pc_x[k] * nsh_642[k];

        t_853[k] = f_21 * msh_492[k]
                   + f_3 * pc_y[k] * nsh_639[k];

        t_854[k] = f_13 * msh_644[k]
                   + f_4 * nsg0_464[k]
                   - f_5 * nsg1_464[k]
                   + f_3 * pc_x[k] * nsh_644[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, pc_x, msh_645, msh_646, msh_647, \
                         msh_648, msh_649, nsh_645, nsh_646, nsh_647, nsh_648, \
                         nsh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_13 * msh_645[k]
                   + f_3 * pc_x[k] * nsh_645[k];

        t_856[k] = f_13 * msh_646[k]
                   + f_3 * pc_x[k] * nsh_646[k];

        t_857[k] = f_13 * msh_647[k]
                   + f_3 * pc_x[k] * nsh_647[k];

        t_858[k] = f_13 * msh_648[k]
                   + f_3 * pc_x[k] * nsh_648[k];

        t_859[k] = f_13 * msh_649[k]
                   + f_3 * pc_x[k] * nsh_649[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pc_x, pc_y, pc_z, msh_477, msh_498, msh_650, \
                         nsg0_460, nsg1_460, nsh_645, nsh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_13 * msh_650[k]
                   + f_3 * pc_x[k] * nsh_650[k];

        t_861[k] = f_21 * msh_498[k]
                   + f_1 * nsg0_460[k]
                   - f_2 * nsg1_460[k]
                   + f_3 * pc_y[k] * nsh_645[k];

        t_862[k] = f_12 * msh_477[k]
                   + f_3 * pc_z[k] * nsh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_y, msh_500, msh_501, msh_502, nsg0_462, \
                         nsg0_463, nsg0_464, nsg1_462, nsg1_463, nsg1_464, nsh_647, nsh_648, \
                         nsh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_21 * msh_500[k]
                   + f_8 * nsg0_462[k]
                   - f_9 * nsg1_462[k]
                   + f_3 * pc_y[k] * nsh_647[k];

        t_864[k] = f_21 * msh_501[k]
                   + f_6 * nsg0_463[k]
                   - f_7 * nsg1_463[k]
                   + f_3 * pc_y[k] * nsh_648[k];

        t_865[k] = f_21 * msh_502[k]
                   + f_4 * nsg0_464[k]
                   - f_5 * nsg1_464[k]
                   + f_3 * pc_y[k] * nsh_649[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, msh_482, msh_503, msh_651, \
                         nsg0_464, nsg0_465, nsg1_464, nsg1_465, nsh_650, \
                         nsh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_21 * msh_503[k]
                   + f_3 * pc_y[k] * nsh_650[k];

        t_867[k] = f_12 * msh_482[k]
                   + f_1 * nsg0_464[k]
                   - f_2 * nsg1_464[k]
                   + f_3 * pc_z[k] * nsh_650[k];

        t_868[k] = f_13 * msh_651[k]
                   + f_1 * nsg0_465[k]
                   - f_2 * nsg1_465[k]
                   + f_3 * pc_x[k] * nsh_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, pc_z, msh_483, msh_504, \
                         msh_506, msh_654, nsg0_468, nsg1_468, nsh_651, nsh_653, \
                         nsh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * msh_504[k]
                   + f_3 * pc_y[k] * nsh_651[k];

        t_870[k] = f_13 * msh_483[k]
                   + f_3 * pc_z[k] * nsh_651[k];

        t_871[k] = f_13 * msh_654[k]
                   + f_8 * nsg0_468[k]
                   - f_9 * nsg1_468[k]
                   + f_3 * pc_x[k] * nsh_654[k];

        t_872[k] = f_14 * msh_506[k]
                   + f_3 * pc_y[k] * nsh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_z, msh_486, msh_656, msh_657, nsg0_470, \
                         nsg0_471, nsg1_470, nsg1_471, nsh_654, nsh_656, \
                         nsh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_13 * msh_656[k]
                   + f_8 * nsg0_470[k]
                   - f_9 * nsg1_470[k]
                   + f_3 * pc_x[k] * nsh_656[k];

        t_874[k] = f_13 * msh_657[k]
                   + f_6 * nsg0_471[k]
                   - f_7 * nsg1_471[k]
                   + f_3 * pc_x[k] * nsh_657[k];

        t_875[k] = f_13 * msh_486[k]
                   + f_3 * pc_z[k] * nsh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_y, msh_509, msh_660, msh_661, nsg0_474, \
                         nsg0_475, nsg1_474, nsg1_475, nsh_656, nsh_660, \
                         nsh_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * msh_509[k]
                   + f_3 * pc_y[k] * nsh_656[k];

        t_877[k] = f_13 * msh_660[k]
                   + f_6 * nsg0_474[k]
                   - f_7 * nsg1_474[k]
                   + f_3 * pc_x[k] * nsh_660[k];

        t_878[k] = f_13 * msh_661[k]
                   + f_4 * nsg0_475[k]
                   - f_5 * nsg1_475[k]
                   + f_3 * pc_x[k] * nsh_661[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, msh_489, msh_513, msh_663, \
                         nsg0_477, nsg1_477, nsh_657, nsh_660, \
                         nsh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * msh_489[k]
                   + f_3 * pc_z[k] * nsh_657[k];

        t_880[k] = f_13 * msh_663[k]
                   + f_4 * nsg0_477[k]
                   - f_5 * nsg1_477[k]
                   + f_3 * pc_x[k] * nsh_663[k];

        t_881[k] = f_14 * msh_513[k]
                   + f_3 * pc_y[k] * nsh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pc_x, msh_665, msh_666, msh_667, msh_668, \
                         nsg0_479, nsg1_479, nsh_665, nsh_666, nsh_667, \
                         nsh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_13 * msh_665[k]
                   + f_4 * nsg0_479[k]
                   - f_5 * nsg1_479[k]
                   + f_3 * pc_x[k] * nsh_665[k];

        t_883[k] = f_13 * msh_666[k]
                   + f_3 * pc_x[k] * nsh_666[k];

        t_884[k] = f_13 * msh_667[k]
                   + f_3 * pc_x[k] * nsh_667[k];

        t_885[k] = f_13 * msh_668[k]
                   + f_3 * pc_x[k] * nsh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pc_x, pc_y, msh_519, msh_669, msh_670, \
                         msh_671, nsg0_475, nsg1_475, nsh_666, nsh_669, nsh_670, \
                         nsh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_13 * msh_669[k]
                   + f_3 * pc_x[k] * nsh_669[k];

        t_887[k] = f_13 * msh_670[k]
                   + f_3 * pc_x[k] * nsh_670[k];

        t_888[k] = f_13 * msh_671[k]
                   + f_3 * pc_x[k] * nsh_671[k];

        t_889[k] = f_14 * msh_519[k]
                   + f_1 * nsg0_475[k]
                   - f_2 * nsg1_475[k]
                   + f_3 * pc_y[k] * nsh_666[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pc_y, pc_z, msh_498, msh_521, msh_522, nsg0_477, \
                         nsg0_478, nsg1_477, nsg1_478, nsh_666, nsh_668, \
                         nsh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * msh_498[k]
                   + f_3 * pc_z[k] * nsh_666[k];

        t_891[k] = f_14 * msh_521[k]
                   + f_8 * nsg0_477[k]
                   - f_9 * nsg1_477[k]
                   + f_3 * pc_y[k] * nsh_668[k];

        t_892[k] = f_14 * msh_522[k]
                   + f_6 * nsg0_478[k]
                   - f_7 * nsg1_478[k]
                   + f_3 * pc_y[k] * nsh_669[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, pc_z, msh_503, msh_523, msh_524, nsg0_479, \
                         nsg1_479, nsh_670, nsh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_14 * msh_523[k]
                   + f_4 * nsg0_479[k]
                   - f_5 * nsg1_479[k]
                   + f_3 * pc_y[k] * nsh_670[k];

        t_894[k] = f_14 * msh_524[k]
                   + f_3 * pc_y[k] * nsh_671[k];

        t_895[k] = f_13 * msh_503[k]
                   + f_1 * nsg0_479[k]
                   - f_2 * nsg1_479[k]
                   + f_3 * pc_z[k] * nsh_671[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_x, pc_y, pc_z, msh_504, msh_525, msh_672, \
                         nsg0_480, nsg1_480, nsh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_13 * msh_672[k]
                   + f_1 * nsg0_480[k]
                   - f_2 * nsg1_480[k]
                   + f_3 * pc_x[k] * nsh_672[k];

        t_897[k] = f_13 * msh_525[k]
                   + f_3 * pc_y[k] * nsh_672[k];

        t_898[k] = f_14 * msh_504[k]
                   + f_3 * pc_z[k] * nsh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, msh_527, msh_675, msh_677, nsg0_483, \
                         nsg0_485, nsg1_483, nsg1_485, nsh_674, nsh_675, \
                         nsh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_13 * msh_675[k]
                   + f_8 * nsg0_483[k]
                   - f_9 * nsg1_483[k]
                   + f_3 * pc_x[k] * nsh_675[k];

        t_900[k] = f_13 * msh_527[k]
                   + f_3 * pc_y[k] * nsh_674[k];

        t_901[k] = f_13 * msh_677[k]
                   + f_8 * nsg0_485[k]
                   - f_9 * nsg1_485[k]
                   + f_3 * pc_x[k] * nsh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, msh_507, msh_530, msh_678, \
                         nsg0_486, nsg1_486, nsh_675, nsh_677, \
                         nsh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_13 * msh_678[k]
                   + f_6 * nsg0_486[k]
                   - f_7 * nsg1_486[k]
                   + f_3 * pc_x[k] * nsh_678[k];

        t_903[k] = f_14 * msh_507[k]
                   + f_3 * pc_z[k] * nsh_675[k];

        t_904[k] = f_13 * msh_530[k]
                   + f_3 * pc_y[k] * nsh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, msh_510, msh_681, msh_682, nsg0_489, \
                         nsg0_490, nsg1_489, nsg1_490, nsh_678, nsh_681, \
                         nsh_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_13 * msh_681[k]
                   + f_6 * nsg0_489[k]
                   - f_7 * nsg1_489[k]
                   + f_3 * pc_x[k] * nsh_681[k];

        t_906[k] = f_13 * msh_682[k]
                   + f_4 * nsg0_490[k]
                   - f_5 * nsg1_490[k]
                   + f_3 * pc_x[k] * nsh_682[k];

        t_907[k] = f_14 * msh_510[k]
                   + f_3 * pc_z[k] * nsh_678[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, msh_534, msh_684, msh_686, nsg0_492, \
                         nsg0_494, nsg1_492, nsg1_494, nsh_681, nsh_684, \
                         nsh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_13 * msh_684[k]
                   + f_4 * nsg0_492[k]
                   - f_5 * nsg1_492[k]
                   + f_3 * pc_x[k] * nsh_684[k];

        t_909[k] = f_13 * msh_534[k]
                   + f_3 * pc_y[k] * nsh_681[k];

        t_910[k] = f_13 * msh_686[k]
                   + f_4 * nsg0_494[k]
                   - f_5 * nsg1_494[k]
                   + f_3 * pc_x[k] * nsh_686[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pc_x, msh_687, msh_688, msh_689, \
                         msh_690, msh_691, nsh_687, nsh_688, nsh_689, nsh_690, \
                         nsh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_13 * msh_687[k]
                   + f_3 * pc_x[k] * nsh_687[k];

        t_912[k] = f_13 * msh_688[k]
                   + f_3 * pc_x[k] * nsh_688[k];

        t_913[k] = f_13 * msh_689[k]
                   + f_3 * pc_x[k] * nsh_689[k];

        t_914[k] = f_13 * msh_690[k]
                   + f_3 * pc_x[k] * nsh_690[k];

        t_915[k] = f_13 * msh_691[k]
                   + f_3 * pc_x[k] * nsh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pc_x, pc_y, pc_z, msh_519, msh_540, msh_692, \
                         nsg0_490, nsg1_490, nsh_687, nsh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_13 * msh_692[k]
                   + f_3 * pc_x[k] * nsh_692[k];

        t_917[k] = f_13 * msh_540[k]
                   + f_1 * nsg0_490[k]
                   - f_2 * nsg1_490[k]
                   + f_3 * pc_y[k] * nsh_687[k];

        t_918[k] = f_14 * msh_519[k]
                   + f_3 * pc_z[k] * nsh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pc_y, msh_542, msh_543, msh_544, nsg0_492, \
                         nsg0_493, nsg0_494, nsg1_492, nsg1_493, nsg1_494, nsh_689, nsh_690, \
                         nsh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_13 * msh_542[k]
                   + f_8 * nsg0_492[k]
                   - f_9 * nsg1_492[k]
                   + f_3 * pc_y[k] * nsh_689[k];

        t_920[k] = f_13 * msh_543[k]
                   + f_6 * nsg0_493[k]
                   - f_7 * nsg1_493[k]
                   + f_3 * pc_y[k] * nsh_690[k];

        t_921[k] = f_13 * msh_544[k]
                   + f_4 * nsg0_494[k]
                   - f_5 * nsg1_494[k]
                   + f_3 * pc_y[k] * nsh_691[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pc_x, pc_y, pc_z, msh_524, msh_545, msh_693, \
                         nsg0_494, nsg0_495, nsg1_494, nsg1_495, nsh_692, \
                         nsh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * msh_545[k]
                   + f_3 * pc_y[k] * nsh_692[k];

        t_923[k] = f_14 * msh_524[k]
                   + f_1 * nsg0_494[k]
                   - f_2 * nsg1_494[k]
                   + f_3 * pc_z[k] * nsh_692[k];

        t_924[k] = f_13 * msh_693[k]
                   + f_1 * nsg0_495[k]
                   - f_2 * nsg1_495[k]
                   + f_3 * pc_x[k] * nsh_693[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_756 = buffer.data(msi0 + 756);
    const auto *msi0_759 = buffer.data(msi0 + 759);
    const auto *msi0_761 = buffer.data(msi0 + 761);
    const auto *msi0_762 = buffer.data(msi0 + 762);
    const auto *msi0_765 = buffer.data(msi0 + 765);
    const auto *msi0_766 = buffer.data(msi0 + 766);
    const auto *msi0_768 = buffer.data(msi0 + 768);
    const auto *msi0_770 = buffer.data(msi0 + 770);
    const auto *msi0_783 = buffer.data(msi0 + 783);
    const auto *msi0_784 = buffer.data(msi0 + 784);
    const auto *msi0_787 = buffer.data(msi0 + 787);
    const auto *msi0_790 = buffer.data(msi0 + 790);

    const auto *msh_525 = buffer.data(msh + 525);
    const auto *msh_528 = buffer.data(msh + 528);
    const auto *msh_531 = buffer.data(msh + 531);
    const auto *msh_540 = buffer.data(msh + 540);
    const auto *msh_545 = buffer.data(msh + 545);
    const auto *msh_546 = buffer.data(msh + 546);
    const auto *msh_548 = buffer.data(msh + 548);
    const auto *msh_549 = buffer.data(msh + 549);
    const auto *msh_551 = buffer.data(msh + 551);
    const auto *msh_552 = buffer.data(msh + 552);
    const auto *msh_555 = buffer.data(msh + 555);
    const auto *msh_561 = buffer.data(msh + 561);
    const auto *msh_563 = buffer.data(msh + 563);
    const auto *msh_564 = buffer.data(msh + 564);
    const auto *msh_565 = buffer.data(msh + 565);
    const auto *msh_566 = buffer.data(msh + 566);
    const auto *msh_567 = buffer.data(msh + 567);
    const auto *msh_568 = buffer.data(msh + 568);
    const auto *msh_569 = buffer.data(msh + 569);
    const auto *msh_570 = buffer.data(msh + 570);
    const auto *msh_572 = buffer.data(msh + 572);
    const auto *msh_573 = buffer.data(msh + 573);
    const auto *msh_575 = buffer.data(msh + 575);
    const auto *msh_576 = buffer.data(msh + 576);
    const auto *msh_582 = buffer.data(msh + 582);
    const auto *msh_584 = buffer.data(msh + 584);
    const auto *msh_585 = buffer.data(msh + 585);
    const auto *msh_586 = buffer.data(msh + 586);
    const auto *msh_587 = buffer.data(msh + 587);
    const auto *msh_588 = buffer.data(msh + 588);
    const auto *msh_591 = buffer.data(msh + 591);
    const auto *msh_593 = buffer.data(msh + 593);
    const auto *msh_597 = buffer.data(msh + 597);
    const auto *msh_603 = buffer.data(msh + 603);
    const auto *msh_608 = buffer.data(msh + 608);
    const auto *msh_609 = buffer.data(msh + 609);
    const auto *msh_611 = buffer.data(msh + 611);
    const auto *msh_696 = buffer.data(msh + 696);
    const auto *msh_698 = buffer.data(msh + 698);
    const auto *msh_699 = buffer.data(msh + 699);
    const auto *msh_702 = buffer.data(msh + 702);
    const auto *msh_703 = buffer.data(msh + 703);
    const auto *msh_705 = buffer.data(msh + 705);
    const auto *msh_707 = buffer.data(msh + 707);
    const auto *msh_708 = buffer.data(msh + 708);
    const auto *msh_709 = buffer.data(msh + 709);
    const auto *msh_710 = buffer.data(msh + 710);
    const auto *msh_711 = buffer.data(msh + 711);
    const auto *msh_712 = buffer.data(msh + 712);
    const auto *msh_713 = buffer.data(msh + 713);
    const auto *msh_729 = buffer.data(msh + 729);
    const auto *msh_730 = buffer.data(msh + 730);
    const auto *msh_731 = buffer.data(msh + 731);
    const auto *msh_732 = buffer.data(msh + 732);
    const auto *msh_733 = buffer.data(msh + 733);
    const auto *msh_734 = buffer.data(msh + 734);
    const auto *msh_735 = buffer.data(msh + 735);
    const auto *msh_740 = buffer.data(msh + 740);
    const auto *msh_744 = buffer.data(msh + 744);
    const auto *msh_749 = buffer.data(msh + 749);
    const auto *msh_750 = buffer.data(msh + 750);
    const auto *msh_751 = buffer.data(msh + 751);
    const auto *msh_752 = buffer.data(msh + 752);
    const auto *msh_753 = buffer.data(msh + 753);
    const auto *msh_755 = buffer.data(msh + 755);
    const auto *msh_756 = buffer.data(msh + 756);
    const auto *msh_759 = buffer.data(msh + 759);
    const auto *msh_762 = buffer.data(msh + 762);
    const auto *msh_766 = buffer.data(msh + 766);
    const auto *msh_771 = buffer.data(msh + 771);
    const auto *msh_773 = buffer.data(msh + 773);
    const auto *msh_774 = buffer.data(msh + 774);
    const auto *msh_775 = buffer.data(msh + 775);
    const auto *msh_776 = buffer.data(msh + 776);
    const auto *msh_782 = buffer.data(msh + 782);

    const auto *msi1_756 = buffer.data(msi1 + 756);
    const auto *msi1_759 = buffer.data(msi1 + 759);
    const auto *msi1_761 = buffer.data(msi1 + 761);
    const auto *msi1_762 = buffer.data(msi1 + 762);
    const auto *msi1_765 = buffer.data(msi1 + 765);
    const auto *msi1_766 = buffer.data(msi1 + 766);
    const auto *msi1_768 = buffer.data(msi1 + 768);
    const auto *msi1_770 = buffer.data(msi1 + 770);
    const auto *msi1_783 = buffer.data(msi1 + 783);
    const auto *msi1_784 = buffer.data(msi1 + 784);
    const auto *msi1_787 = buffer.data(msi1 + 787);
    const auto *msi1_790 = buffer.data(msi1 + 790);

    const auto *nsg0_498 = buffer.data(nsg0 + 498);
    const auto *nsg0_500 = buffer.data(nsg0 + 500);
    const auto *nsg0_501 = buffer.data(nsg0 + 501);
    const auto *nsg0_504 = buffer.data(nsg0 + 504);
    const auto *nsg0_505 = buffer.data(nsg0 + 505);
    const auto *nsg0_507 = buffer.data(nsg0 + 507);
    const auto *nsg0_508 = buffer.data(nsg0 + 508);
    const auto *nsg0_509 = buffer.data(nsg0 + 509);
    const auto *nsg0_520 = buffer.data(nsg0 + 520);
    const auto *nsg0_522 = buffer.data(nsg0 + 522);
    const auto *nsg0_523 = buffer.data(nsg0 + 523);
    const auto *nsg0_524 = buffer.data(nsg0 + 524);
    const auto *nsg0_525 = buffer.data(nsg0 + 525);
    const auto *nsg0_526 = buffer.data(nsg0 + 526);
    const auto *nsg0_527 = buffer.data(nsg0 + 527);
    const auto *nsg0_528 = buffer.data(nsg0 + 528);
    const auto *nsg0_529 = buffer.data(nsg0 + 529);
    const auto *nsg0_530 = buffer.data(nsg0 + 530);
    const auto *nsg0_534 = buffer.data(nsg0 + 534);
    const auto *nsg0_535 = buffer.data(nsg0 + 535);
    const auto *nsg0_536 = buffer.data(nsg0 + 536);
    const auto *nsg0_537 = buffer.data(nsg0 + 537);
    const auto *nsg0_538 = buffer.data(nsg0 + 538);
    const auto *nsg0_539 = buffer.data(nsg0 + 539);
    const auto *nsg0_540 = buffer.data(nsg0 + 540);
    const auto *nsg0_542 = buffer.data(nsg0 + 542);
    const auto *nsg0_543 = buffer.data(nsg0 + 543);
    const auto *nsg0_545 = buffer.data(nsg0 + 545);
    const auto *nsg0_546 = buffer.data(nsg0 + 546);
    const auto *nsg0_550 = buffer.data(nsg0 + 550);
    const auto *nsg0_551 = buffer.data(nsg0 + 551);
    const auto *nsg0_552 = buffer.data(nsg0 + 552);
    const auto *nsg0_554 = buffer.data(nsg0 + 554);
    const auto *nsg0_560 = buffer.data(nsg0 + 560);

    const auto *nsg1_498 = buffer.data(nsg1 + 498);
    const auto *nsg1_500 = buffer.data(nsg1 + 500);
    const auto *nsg1_501 = buffer.data(nsg1 + 501);
    const auto *nsg1_504 = buffer.data(nsg1 + 504);
    const auto *nsg1_505 = buffer.data(nsg1 + 505);
    const auto *nsg1_507 = buffer.data(nsg1 + 507);
    const auto *nsg1_508 = buffer.data(nsg1 + 508);
    const auto *nsg1_509 = buffer.data(nsg1 + 509);
    const auto *nsg1_520 = buffer.data(nsg1 + 520);
    const auto *nsg1_522 = buffer.data(nsg1 + 522);
    const auto *nsg1_523 = buffer.data(nsg1 + 523);
    const auto *nsg1_524 = buffer.data(nsg1 + 524);
    const auto *nsg1_525 = buffer.data(nsg1 + 525);
    const auto *nsg1_526 = buffer.data(nsg1 + 526);
    const auto *nsg1_527 = buffer.data(nsg1 + 527);
    const auto *nsg1_528 = buffer.data(nsg1 + 528);
    const auto *nsg1_529 = buffer.data(nsg1 + 529);
    const auto *nsg1_530 = buffer.data(nsg1 + 530);
    const auto *nsg1_534 = buffer.data(nsg1 + 534);
    const auto *nsg1_535 = buffer.data(nsg1 + 535);
    const auto *nsg1_536 = buffer.data(nsg1 + 536);
    const auto *nsg1_537 = buffer.data(nsg1 + 537);
    const auto *nsg1_538 = buffer.data(nsg1 + 538);
    const auto *nsg1_539 = buffer.data(nsg1 + 539);
    const auto *nsg1_540 = buffer.data(nsg1 + 540);
    const auto *nsg1_542 = buffer.data(nsg1 + 542);
    const auto *nsg1_543 = buffer.data(nsg1 + 543);
    const auto *nsg1_545 = buffer.data(nsg1 + 545);
    const auto *nsg1_546 = buffer.data(nsg1 + 546);
    const auto *nsg1_550 = buffer.data(nsg1 + 550);
    const auto *nsg1_551 = buffer.data(nsg1 + 551);
    const auto *nsg1_552 = buffer.data(nsg1 + 552);
    const auto *nsg1_554 = buffer.data(nsg1 + 554);
    const auto *nsg1_560 = buffer.data(nsg1 + 560);

    const auto *nsh_693 = buffer.data(nsh + 693);
    const auto *nsh_695 = buffer.data(nsh + 695);
    const auto *nsh_696 = buffer.data(nsh + 696);
    const auto *nsh_698 = buffer.data(nsh + 698);
    const auto *nsh_699 = buffer.data(nsh + 699);
    const auto *nsh_702 = buffer.data(nsh + 702);
    const auto *nsh_703 = buffer.data(nsh + 703);
    const auto *nsh_705 = buffer.data(nsh + 705);
    const auto *nsh_707 = buffer.data(nsh + 707);
    const auto *nsh_708 = buffer.data(nsh + 708);
    const auto *nsh_709 = buffer.data(nsh + 709);
    const auto *nsh_710 = buffer.data(nsh + 710);
    const auto *nsh_711 = buffer.data(nsh + 711);
    const auto *nsh_712 = buffer.data(nsh + 712);
    const auto *nsh_713 = buffer.data(nsh + 713);
    const auto *nsh_714 = buffer.data(nsh + 714);
    const auto *nsh_716 = buffer.data(nsh + 716);
    const auto *nsh_717 = buffer.data(nsh + 717);
    const auto *nsh_719 = buffer.data(nsh + 719);
    const auto *nsh_720 = buffer.data(nsh + 720);
    const auto *nsh_723 = buffer.data(nsh + 723);
    const auto *nsh_729 = buffer.data(nsh + 729);
    const auto *nsh_730 = buffer.data(nsh + 730);
    const auto *nsh_731 = buffer.data(nsh + 731);
    const auto *nsh_732 = buffer.data(nsh + 732);
    const auto *nsh_733 = buffer.data(nsh + 733);
    const auto *nsh_734 = buffer.data(nsh + 734);
    const auto *nsh_735 = buffer.data(nsh + 735);
    const auto *nsh_736 = buffer.data(nsh + 736);
    const auto *nsh_737 = buffer.data(nsh + 737);
    const auto *nsh_738 = buffer.data(nsh + 738);
    const auto *nsh_739 = buffer.data(nsh + 739);
    const auto *nsh_740 = buffer.data(nsh + 740);
    const auto *nsh_741 = buffer.data(nsh + 741);
    const auto *nsh_742 = buffer.data(nsh + 742);
    const auto *nsh_743 = buffer.data(nsh + 743);
    const auto *nsh_744 = buffer.data(nsh + 744);
    const auto *nsh_749 = buffer.data(nsh + 749);
    const auto *nsh_750 = buffer.data(nsh + 750);
    const auto *nsh_751 = buffer.data(nsh + 751);
    const auto *nsh_752 = buffer.data(nsh + 752);
    const auto *nsh_753 = buffer.data(nsh + 753);
    const auto *nsh_754 = buffer.data(nsh + 754);
    const auto *nsh_755 = buffer.data(nsh + 755);
    const auto *nsh_756 = buffer.data(nsh + 756);
    const auto *nsh_757 = buffer.data(nsh + 757);
    const auto *nsh_758 = buffer.data(nsh + 758);
    const auto *nsh_759 = buffer.data(nsh + 759);
    const auto *nsh_761 = buffer.data(nsh + 761);
    const auto *nsh_762 = buffer.data(nsh + 762);
    const auto *nsh_763 = buffer.data(nsh + 763);
    const auto *nsh_765 = buffer.data(nsh + 765);
    const auto *nsh_766 = buffer.data(nsh + 766);
    const auto *nsh_771 = buffer.data(nsh + 771);
    const auto *nsh_772 = buffer.data(nsh + 772);
    const auto *nsh_773 = buffer.data(nsh + 773);
    const auto *nsh_774 = buffer.data(nsh + 774);
    const auto *nsh_775 = buffer.data(nsh + 775);
    const auto *nsh_776 = buffer.data(nsh + 776);
    const auto *nsh_777 = buffer.data(nsh + 777);
    const auto *nsh_779 = buffer.data(nsh + 779);
    const auto *nsh_780 = buffer.data(nsh + 780);
    const auto *nsh_782 = buffer.data(nsh + 782);

#pragma omp simd aligned(t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, msh_525, msh_546, \
                         msh_548, msh_696, nsg0_498, nsg1_498, nsh_693, nsh_695, \
                         nsh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_12 * msh_546[k]
                   + f_3 * pc_y[k] * nsh_693[k];

        t_926[k] = f_21 * msh_525[k]
                   + f_3 * pc_z[k] * nsh_693[k];

        t_927[k] = f_13 * msh_696[k]
                   + f_8 * nsg0_498[k]
                   - f_9 * nsg1_498[k]
                   + f_3 * pc_x[k] * nsh_696[k];

        t_928[k] = f_12 * msh_548[k]
                   + f_3 * pc_y[k] * nsh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pc_x, pc_z, msh_528, msh_698, msh_699, nsg0_500, \
                         nsg0_501, nsg1_500, nsg1_501, nsh_696, nsh_698, \
                         nsh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_13 * msh_698[k]
                   + f_8 * nsg0_500[k]
                   - f_9 * nsg1_500[k]
                   + f_3 * pc_x[k] * nsh_698[k];

        t_930[k] = f_13 * msh_699[k]
                   + f_6 * nsg0_501[k]
                   - f_7 * nsg1_501[k]
                   + f_3 * pc_x[k] * nsh_699[k];

        t_931[k] = f_21 * msh_528[k]
                   + f_3 * pc_z[k] * nsh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pc_x, pc_y, msh_551, msh_702, msh_703, nsg0_504, \
                         nsg0_505, nsg1_504, nsg1_505, nsh_698, nsh_702, \
                         nsh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * msh_551[k]
                   + f_3 * pc_y[k] * nsh_698[k];

        t_933[k] = f_13 * msh_702[k]
                   + f_6 * nsg0_504[k]
                   - f_7 * nsg1_504[k]
                   + f_3 * pc_x[k] * nsh_702[k];

        t_934[k] = f_13 * msh_703[k]
                   + f_4 * nsg0_505[k]
                   - f_5 * nsg1_505[k]
                   + f_3 * pc_x[k] * nsh_703[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pc_x, pc_y, pc_z, msh_531, msh_555, msh_705, \
                         nsg0_507, nsg1_507, nsh_699, nsh_702, \
                         nsh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_21 * msh_531[k]
                   + f_3 * pc_z[k] * nsh_699[k];

        t_936[k] = f_13 * msh_705[k]
                   + f_4 * nsg0_507[k]
                   - f_5 * nsg1_507[k]
                   + f_3 * pc_x[k] * nsh_705[k];

        t_937[k] = f_12 * msh_555[k]
                   + f_3 * pc_y[k] * nsh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pc_x, msh_707, msh_708, msh_709, msh_710, \
                         nsg0_509, nsg1_509, nsh_707, nsh_708, nsh_709, \
                         nsh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_13 * msh_707[k]
                   + f_4 * nsg0_509[k]
                   - f_5 * nsg1_509[k]
                   + f_3 * pc_x[k] * nsh_707[k];

        t_939[k] = f_13 * msh_708[k]
                   + f_3 * pc_x[k] * nsh_708[k];

        t_940[k] = f_13 * msh_709[k]
                   + f_3 * pc_x[k] * nsh_709[k];

        t_941[k] = f_13 * msh_710[k]
                   + f_3 * pc_x[k] * nsh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, msh_561, msh_711, msh_712, \
                         msh_713, nsg0_505, nsg1_505, nsh_708, nsh_711, nsh_712, \
                         nsh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_13 * msh_711[k]
                   + f_3 * pc_x[k] * nsh_711[k];

        t_943[k] = f_13 * msh_712[k]
                   + f_3 * pc_x[k] * nsh_712[k];

        t_944[k] = f_13 * msh_713[k]
                   + f_3 * pc_x[k] * nsh_713[k];

        t_945[k] = f_12 * msh_561[k]
                   + f_1 * nsg0_505[k]
                   - f_2 * nsg1_505[k]
                   + f_3 * pc_y[k] * nsh_708[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pc_y, pc_z, msh_540, msh_563, msh_564, nsg0_507, \
                         nsg0_508, nsg1_507, nsg1_508, nsh_708, nsh_710, \
                         nsh_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_21 * msh_540[k]
                   + f_3 * pc_z[k] * nsh_708[k];

        t_947[k] = f_12 * msh_563[k]
                   + f_8 * nsg0_507[k]
                   - f_9 * nsg1_507[k]
                   + f_3 * pc_y[k] * nsh_710[k];

        t_948[k] = f_12 * msh_564[k]
                   + f_6 * nsg0_508[k]
                   - f_7 * nsg1_508[k]
                   + f_3 * pc_y[k] * nsh_711[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pa_y, pc_y, pc_z, msi0_756, msh_545, \
                         msh_565, msh_566, msi1_756, nsg0_509, nsg1_509, nsh_712, \
                         nsh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_12 * msh_565[k]
                   + f_4 * nsg0_509[k]
                   - f_5 * nsg1_509[k]
                   + f_3 * pc_y[k] * nsh_712[k];

        t_950[k] = f_12 * msh_566[k]
                   + f_3 * pc_y[k] * nsh_713[k];

        t_951[k] = f_21 * msh_545[k]
                   + f_1 * nsg0_509[k]
                   - f_2 * nsg1_509[k]
                   + f_3 * pc_z[k] * nsh_713[k];

        t_952[k] = pa_y[k] * msi0_756[k]
                   - f_10 * pc_y[k] * msi1_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pa_y, pc_y, pc_z, msi0_759, msh_546, \
                         msh_567, msh_568, msh_569, msi1_759, nsh_714, \
                         nsh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_11 * msh_567[k]
                   + f_3 * pc_y[k] * nsh_714[k];

        t_954[k] = f_20 * msh_546[k]
                   + f_3 * pc_z[k] * nsh_714[k];

        t_955[k] = pa_y[k] * msi0_759[k]
                   + f_12 * msh_568[k]
                   - f_10 * pc_y[k] * msi1_759[k];

        t_956[k] = f_11 * msh_569[k]
                   + f_3 * pc_y[k] * nsh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pa_y, pc_y, pc_z, msi0_761, msi0_762, \
                         msh_549, msh_570, msh_572, msi1_761, msi1_762, nsh_717, \
                         nsh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pa_y[k] * msi0_761[k]
                   - f_10 * pc_y[k] * msi1_761[k];

        t_958[k] = pa_y[k] * msi0_762[k]
                   + f_13 * msh_570[k]
                   - f_10 * pc_y[k] * msi1_762[k];

        t_959[k] = f_20 * msh_549[k]
                   + f_3 * pc_z[k] * nsh_717[k];

        t_960[k] = f_11 * msh_572[k]
                   + f_3 * pc_y[k] * nsh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pa_y, pc_y, pc_z, msi0_765, msi0_766, msh_552, \
                         msh_573, msi1_765, msi1_766, nsh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pa_y[k] * msi0_765[k]
                   - f_10 * pc_y[k] * msi1_765[k];

        t_962[k] = pa_y[k] * msi0_766[k]
                   + f_14 * msh_573[k]
                   - f_10 * pc_y[k] * msi1_766[k];

        t_963[k] = f_20 * msh_552[k]
                   + f_3 * pc_z[k] * nsh_720[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pa_y, pc_x, pc_y, msi0_768, msi0_770, \
                         msh_575, msh_576, msh_729, msi1_768, msi1_770, nsh_723, \
                         nsh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pa_y[k] * msi0_768[k]
                   + f_12 * msh_575[k]
                   - f_10 * pc_y[k] * msi1_768[k];

        t_965[k] = f_11 * msh_576[k]
                   + f_3 * pc_y[k] * nsh_723[k];

        t_966[k] = pa_y[k] * msi0_770[k]
                   - f_10 * pc_y[k] * msi1_770[k];

        t_967[k] = f_13 * msh_729[k]
                   + f_3 * pc_x[k] * nsh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, msh_730, msh_731, msh_732, \
                         msh_733, msh_734, nsh_730, nsh_731, nsh_732, nsh_733, \
                         nsh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_13 * msh_730[k]
                   + f_3 * pc_x[k] * nsh_730[k];

        t_969[k] = f_13 * msh_731[k]
                   + f_3 * pc_x[k] * nsh_731[k];

        t_970[k] = f_13 * msh_732[k]
                   + f_3 * pc_x[k] * nsh_732[k];

        t_971[k] = f_13 * msh_733[k]
                   + f_3 * pc_x[k] * nsh_733[k];

        t_972[k] = f_13 * msh_734[k]
                   + f_3 * pc_x[k] * nsh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pc_y, pc_z, msh_561, msh_582, msh_584, nsg0_520, \
                         nsg0_522, nsg1_520, nsg1_522, nsh_729, \
                         nsh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_11 * msh_582[k]
                   + f_1 * nsg0_520[k]
                   - f_2 * nsg1_520[k]
                   + f_3 * pc_y[k] * nsh_729[k];

        t_974[k] = f_20 * msh_561[k]
                   + f_3 * pc_z[k] * nsh_729[k];

        t_975[k] = f_11 * msh_584[k]
                   + f_8 * nsg0_522[k]
                   - f_9 * nsg1_522[k]
                   + f_3 * pc_y[k] * nsh_731[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_y, msh_585, msh_586, msh_587, nsg0_523, \
                         nsg0_524, nsg1_523, nsg1_524, nsh_732, nsh_733, \
                         nsh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_11 * msh_585[k]
                   + f_6 * nsg0_523[k]
                   - f_7 * nsg1_523[k]
                   + f_3 * pc_y[k] * nsh_732[k];

        t_977[k] = f_11 * msh_586[k]
                   + f_4 * nsg0_524[k]
                   - f_5 * nsg1_524[k]
                   + f_3 * pc_y[k] * nsh_733[k];

        t_978[k] = f_11 * msh_587[k]
                   + f_3 * pc_y[k] * nsh_734[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pa_y, pc_x, pc_y, pc_z, msi0_783, \
                         msh_567, msh_735, msi1_783, nsg0_525, nsg1_525, \
                         nsh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pa_y[k] * msi0_783[k]
                   - f_10 * pc_y[k] * msi1_783[k];

        t_980[k] = f_13 * msh_735[k]
                   + f_1 * nsg0_525[k]
                   - f_2 * nsg1_525[k]
                   + f_3 * pc_x[k] * nsh_735[k];

        t_981[k] = f_3 * pc_y[k] * nsh_735[k];

        t_982[k] = f_19 * msh_567[k]
                   + f_3 * pc_z[k] * nsh_735[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_x, pc_y, msh_740, nsg0_525, nsg0_530, \
                         nsg1_525, nsg1_530, nsh_736, nsh_737, \
                         nsh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_4 * nsg0_525[k]
                   - f_5 * nsg1_525[k]
                   + f_3 * pc_y[k] * nsh_736[k];

        t_984[k] = f_3 * pc_y[k] * nsh_737[k];

        t_985[k] = f_13 * msh_740[k]
                   + f_8 * nsg0_530[k]
                   - f_9 * nsg1_530[k]
                   + f_3 * pc_x[k] * nsh_740[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_y, nsg0_526, nsg0_527, nsg1_526, nsg1_527, \
                         nsh_738, nsh_739, nsh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_6 * nsg0_526[k]
                   - f_7 * nsg1_526[k]
                   + f_3 * pc_y[k] * nsh_738[k];

        t_987[k] = f_4 * nsg0_527[k]
                   - f_5 * nsg1_527[k]
                   + f_3 * pc_y[k] * nsh_739[k];

        t_988[k] = f_3 * pc_y[k] * nsh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_y, msh_744, nsg0_528, nsg0_529, \
                         nsg0_534, nsg1_528, nsg1_529, nsg1_534, nsh_741, nsh_742, \
                         nsh_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_13 * msh_744[k]
                   + f_6 * nsg0_534[k]
                   - f_7 * nsg1_534[k]
                   + f_3 * pc_x[k] * nsh_744[k];

        t_990[k] = f_8 * nsg0_528[k]
                   - f_9 * nsg1_528[k]
                   + f_3 * pc_y[k] * nsh_741[k];

        t_991[k] = f_6 * nsg0_529[k]
                   - f_7 * nsg1_529[k]
                   + f_3 * pc_y[k] * nsh_742[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, pc_x, pc_y, msh_749, msh_750, nsg0_530, \
                         nsg0_539, nsg1_530, nsg1_539, nsh_743, nsh_744, nsh_749, \
                         nsh_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_4 * nsg0_530[k]
                   - f_5 * nsg1_530[k]
                   + f_3 * pc_y[k] * nsh_743[k];

        t_993[k] = f_3 * pc_y[k] * nsh_744[k];

        t_994[k] = f_13 * msh_749[k]
                   + f_4 * nsg0_539[k]
                   - f_5 * nsg1_539[k]
                   + f_3 * pc_x[k] * nsh_749[k];

        t_995[k] = f_13 * msh_750[k]
                   + f_3 * pc_x[k] * nsh_750[k];
    }

#pragma omp simd aligned(t_996, t_997, t_998, t_999, t_1000, pc_x, pc_y, msh_751, msh_752, \
                         msh_753, msh_755, nsh_749, nsh_751, nsh_752, nsh_753, \
                         nsh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_996[k] = f_13 * msh_751[k]
                   + f_3 * pc_x[k] * nsh_751[k];

        t_997[k] = f_13 * msh_752[k]
                   + f_3 * pc_x[k] * nsh_752[k];

        t_998[k] = f_13 * msh_753[k]
                   + f_3 * pc_x[k] * nsh_753[k];

        t_999[k] = f_3 * pc_y[k] * nsh_749[k];

        t_1000[k] = f_13 * msh_755[k]
                    + f_3 * pc_x[k] * nsh_755[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, pc_y, nsg0_535, nsg0_536, nsg0_537, nsg1_535, \
                         nsg1_536, nsg1_537, nsh_750, nsh_751, \
                         nsh_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_1 * nsg0_535[k]
                    - f_2 * nsg1_535[k]
                    + f_3 * pc_y[k] * nsh_750[k];

        t_1002[k] = f_16 * nsg0_536[k]
                    - f_17 * nsg1_536[k]
                    + f_3 * pc_y[k] * nsh_751[k];

        t_1003[k] = f_8 * nsg0_537[k]
                    - f_9 * nsg1_537[k]
                    + f_3 * pc_y[k] * nsh_752[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, msh_587, nsg0_538, \
                         nsg0_539, nsg1_538, nsg1_539, nsh_753, nsh_754, \
                         nsh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_6 * nsg0_538[k]
                    - f_7 * nsg1_538[k]
                    + f_3 * pc_y[k] * nsh_753[k];

        t_1005[k] = f_4 * nsg0_539[k]
                    - f_5 * nsg1_539[k]
                    + f_3 * pc_y[k] * nsh_754[k];

        t_1006[k] = f_3 * pc_y[k] * nsh_755[k];

        t_1007[k] = f_19 * msh_587[k]
                    + f_1 * nsg0_539[k]
                    - f_2 * nsg1_539[k]
                    + f_3 * pc_z[k] * nsh_755[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, msh_588, msh_756, \
                         msh_759, nsg0_540, nsg0_543, nsg1_540, nsg1_543, nsh_756, \
                         nsh_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_12 * msh_756[k]
                    + f_1 * nsg0_540[k]
                    - f_2 * nsg1_540[k]
                    + f_3 * pc_x[k] * nsh_756[k];

        t_1009[k] = f_18 * msh_588[k]
                    + f_3 * pc_y[k] * nsh_756[k];

        t_1010[k] = f_3 * pc_z[k] * nsh_756[k];

        t_1011[k] = f_12 * msh_759[k]
                    + f_8 * nsg0_543[k]
                    - f_9 * nsg1_543[k]
                    + f_3 * pc_x[k] * nsh_759[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, t_1015, pc_x, pc_z, msh_762, nsg0_540, \
                         nsg0_546, nsg1_540, nsg1_546, nsh_757, nsh_758, nsh_759, \
                         nsh_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_3 * pc_z[k] * nsh_757[k];

        t_1013[k] = f_4 * nsg0_540[k]
                    - f_5 * nsg1_540[k]
                    + f_3 * pc_z[k] * nsh_758[k];

        t_1014[k] = f_12 * msh_762[k]
                    + f_6 * nsg0_546[k]
                    - f_7 * nsg1_546[k]
                    + f_3 * pc_x[k] * nsh_762[k];

        t_1015[k] = f_3 * pc_z[k] * nsh_759[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, pc_z, msh_593, msh_766, \
                         nsg0_542, nsg0_550, nsg1_542, nsg1_550, nsh_761, nsh_762, \
                         nsh_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_18 * msh_593[k]
                    + f_3 * pc_y[k] * nsh_761[k];

        t_1017[k] = f_6 * nsg0_542[k]
                    - f_7 * nsg1_542[k]
                    + f_3 * pc_z[k] * nsh_761[k];

        t_1018[k] = f_12 * msh_766[k]
                    + f_4 * nsg0_550[k]
                    - f_5 * nsg1_550[k]
                    + f_3 * pc_x[k] * nsh_766[k];

        t_1019[k] = f_3 * pc_z[k] * nsh_762[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, pc_x, pc_y, pc_z, msh_597, msh_771, \
                         nsg0_543, nsg0_545, nsg1_543, nsg1_545, nsh_763, nsh_765, \
                         nsh_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_4 * nsg0_543[k]
                    - f_5 * nsg1_543[k]
                    + f_3 * pc_z[k] * nsh_763[k];

        t_1021[k] = f_18 * msh_597[k]
                    + f_3 * pc_y[k] * nsh_765[k];

        t_1022[k] = f_8 * nsg0_545[k]
                    - f_9 * nsg1_545[k]
                    + f_3 * pc_z[k] * nsh_765[k];

        t_1023[k] = f_12 * msh_771[k]
                    + f_3 * pc_x[k] * nsh_771[k];
    }

#pragma omp simd aligned(t_1024, t_1025, t_1026, t_1027, t_1028, pc_x, pc_z, msh_773, msh_774, \
                         msh_775, msh_776, nsh_766, nsh_773, nsh_774, nsh_775, \
                         nsh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1024[k] = f_3 * pc_z[k] * nsh_766[k];

        t_1025[k] = f_12 * msh_773[k]
                    + f_3 * pc_x[k] * nsh_773[k];

        t_1026[k] = f_12 * msh_774[k]
                    + f_3 * pc_x[k] * nsh_774[k];

        t_1027[k] = f_12 * msh_775[k]
                    + f_3 * pc_x[k] * nsh_775[k];

        t_1028[k] = f_12 * msh_776[k]
                    + f_3 * pc_x[k] * nsh_776[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pc_y, pc_z, msh_603, nsg0_550, \
                         nsg0_551, nsg1_550, nsg1_551, nsh_771, nsh_772, \
                         nsh_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = f_18 * msh_603[k]
                    + f_1 * nsg0_550[k]
                    - f_2 * nsg1_550[k]
                    + f_3 * pc_y[k] * nsh_771[k];

        t_1030[k] = f_3 * pc_z[k] * nsh_771[k];

        t_1031[k] = f_4 * nsg0_550[k]
                    - f_5 * nsg1_550[k]
                    + f_3 * pc_z[k] * nsh_772[k];

        t_1032[k] = f_6 * nsg0_551[k]
                    - f_7 * nsg1_551[k]
                    + f_3 * pc_z[k] * nsh_773[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, t_1036, pa_z, pc_y, pc_z, msi0_784, msh_608, \
                         msi1_784, nsg0_552, nsg0_554, nsg1_552, nsg1_554, nsh_774, \
                         nsh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_8 * nsg0_552[k]
                    - f_9 * nsg1_552[k]
                    + f_3 * pc_z[k] * nsh_774[k];

        t_1034[k] = f_18 * msh_608[k]
                    + f_3 * pc_y[k] * nsh_776[k];

        t_1035[k] = f_1 * nsg0_554[k]
                    - f_2 * nsg1_554[k]
                    + f_3 * pc_z[k] * nsh_776[k];

        t_1036[k] = pa_z[k] * msi0_784[k]
                    - f_10 * pc_z[k] * msi1_784[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, t_1040, pa_z, pc_y, pc_z, msi0_787, msh_588, \
                         msh_609, msh_611, msi1_787, nsh_777, nsh_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_19 * msh_609[k]
                    + f_3 * pc_y[k] * nsh_777[k];

        t_1038[k] = f_11 * msh_588[k]
                    + f_3 * pc_z[k] * nsh_777[k];

        t_1039[k] = pa_z[k] * msi0_787[k]
                    - f_10 * pc_z[k] * msi1_787[k];

        t_1040[k] = f_19 * msh_611[k]
                    + f_3 * pc_y[k] * nsh_779[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, pa_z, pc_x, pc_z, msi0_790, msh_591, msh_782, \
                         msi1_790, nsg0_560, nsg1_560, nsh_780, \
                         nsh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_12 * msh_782[k]
                    + f_8 * nsg0_560[k]
                    - f_9 * nsg1_560[k]
                    + f_3 * pc_x[k] * nsh_782[k];

        t_1042[k] = pa_z[k] * msi0_790[k]
                    - f_10 * pc_z[k] * msi1_790[k];

        t_1043[k] = f_11 * msh_591[k]
                    + f_3 * pc_z[k] * nsh_780[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t msi0,
                                                          const size_t msh, const size_t msi1,
                                                          const size_t nsg0, const size_t nsg1,
                                                          const size_t nsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_794 = buffer.data(msi0 + 794);
    const auto *msi0_796 = buffer.data(msi0 + 796);
    const auto *msi0_805 = buffer.data(msi0 + 805);

    const auto *msh_594 = buffer.data(msh + 594);
    const auto *msh_595 = buffer.data(msh + 595);
    const auto *msh_603 = buffer.data(msh + 603);
    const auto *msh_608 = buffer.data(msh + 608);
    const auto *msh_609 = buffer.data(msh + 609);
    const auto *msh_612 = buffer.data(msh + 612);
    const auto *msh_614 = buffer.data(msh + 614);
    const auto *msh_615 = buffer.data(msh + 615);
    const auto *msh_618 = buffer.data(msh + 618);
    const auto *msh_624 = buffer.data(msh + 624);
    const auto *msh_626 = buffer.data(msh + 626);
    const auto *msh_627 = buffer.data(msh + 627);
    const auto *msh_628 = buffer.data(msh + 628);
    const auto *msh_629 = buffer.data(msh + 629);
    const auto *msh_630 = buffer.data(msh + 630);
    const auto *msh_632 = buffer.data(msh + 632);
    const auto *msh_633 = buffer.data(msh + 633);
    const auto *msh_635 = buffer.data(msh + 635);
    const auto *msh_636 = buffer.data(msh + 636);
    const auto *msh_639 = buffer.data(msh + 639);
    const auto *msh_645 = buffer.data(msh + 645);
    const auto *msh_647 = buffer.data(msh + 647);
    const auto *msh_648 = buffer.data(msh + 648);
    const auto *msh_649 = buffer.data(msh + 649);
    const auto *msh_650 = buffer.data(msh + 650);
    const auto *msh_651 = buffer.data(msh + 651);
    const auto *msh_653 = buffer.data(msh + 653);
    const auto *msh_654 = buffer.data(msh + 654);
    const auto *msh_656 = buffer.data(msh + 656);
    const auto *msh_657 = buffer.data(msh + 657);
    const auto *msh_660 = buffer.data(msh + 660);
    const auto *msh_666 = buffer.data(msh + 666);
    const auto *msh_668 = buffer.data(msh + 668);
    const auto *msh_669 = buffer.data(msh + 669);
    const auto *msh_670 = buffer.data(msh + 670);
    const auto *msh_671 = buffer.data(msh + 671);
    const auto *msh_672 = buffer.data(msh + 672);
    const auto *msh_674 = buffer.data(msh + 674);
    const auto *msh_677 = buffer.data(msh + 677);
    const auto *msh_681 = buffer.data(msh + 681);
    const auto *msh_687 = buffer.data(msh + 687);
    const auto *msh_689 = buffer.data(msh + 689);
    const auto *msh_690 = buffer.data(msh + 690);
    const auto *msh_691 = buffer.data(msh + 691);
    const auto *msh_692 = buffer.data(msh + 692);
    const auto *msh_786 = buffer.data(msh + 786);
    const auto *msh_791 = buffer.data(msh + 791);
    const auto *msh_792 = buffer.data(msh + 792);
    const auto *msh_793 = buffer.data(msh + 793);
    const auto *msh_794 = buffer.data(msh + 794);
    const auto *msh_795 = buffer.data(msh + 795);
    const auto *msh_796 = buffer.data(msh + 796);
    const auto *msh_797 = buffer.data(msh + 797);
    const auto *msh_798 = buffer.data(msh + 798);
    const auto *msh_801 = buffer.data(msh + 801);
    const auto *msh_803 = buffer.data(msh + 803);
    const auto *msh_804 = buffer.data(msh + 804);
    const auto *msh_807 = buffer.data(msh + 807);
    const auto *msh_808 = buffer.data(msh + 808);
    const auto *msh_810 = buffer.data(msh + 810);
    const auto *msh_812 = buffer.data(msh + 812);
    const auto *msh_813 = buffer.data(msh + 813);
    const auto *msh_814 = buffer.data(msh + 814);
    const auto *msh_815 = buffer.data(msh + 815);
    const auto *msh_816 = buffer.data(msh + 816);
    const auto *msh_817 = buffer.data(msh + 817);
    const auto *msh_818 = buffer.data(msh + 818);
    const auto *msh_819 = buffer.data(msh + 819);
    const auto *msh_822 = buffer.data(msh + 822);
    const auto *msh_824 = buffer.data(msh + 824);
    const auto *msh_825 = buffer.data(msh + 825);
    const auto *msh_828 = buffer.data(msh + 828);
    const auto *msh_829 = buffer.data(msh + 829);
    const auto *msh_831 = buffer.data(msh + 831);
    const auto *msh_833 = buffer.data(msh + 833);
    const auto *msh_834 = buffer.data(msh + 834);
    const auto *msh_835 = buffer.data(msh + 835);
    const auto *msh_836 = buffer.data(msh + 836);
    const auto *msh_837 = buffer.data(msh + 837);
    const auto *msh_838 = buffer.data(msh + 838);
    const auto *msh_839 = buffer.data(msh + 839);
    const auto *msh_840 = buffer.data(msh + 840);
    const auto *msh_843 = buffer.data(msh + 843);
    const auto *msh_845 = buffer.data(msh + 845);
    const auto *msh_846 = buffer.data(msh + 846);
    const auto *msh_849 = buffer.data(msh + 849);
    const auto *msh_850 = buffer.data(msh + 850);
    const auto *msh_852 = buffer.data(msh + 852);
    const auto *msh_854 = buffer.data(msh + 854);
    const auto *msh_855 = buffer.data(msh + 855);
    const auto *msh_856 = buffer.data(msh + 856);
    const auto *msh_857 = buffer.data(msh + 857);
    const auto *msh_858 = buffer.data(msh + 858);
    const auto *msh_859 = buffer.data(msh + 859);
    const auto *msh_860 = buffer.data(msh + 860);
    const auto *msh_861 = buffer.data(msh + 861);

    const auto *msi1_794 = buffer.data(msi1 + 794);
    const auto *msi1_796 = buffer.data(msi1 + 796);
    const auto *msi1_805 = buffer.data(msi1 + 805);

    const auto *nsg0_564 = buffer.data(nsg0 + 564);
    const auto *nsg0_567 = buffer.data(nsg0 + 567);
    const auto *nsg0_568 = buffer.data(nsg0 + 568);
    const auto *nsg0_569 = buffer.data(nsg0 + 569);
    const auto *nsg0_570 = buffer.data(nsg0 + 570);
    const auto *nsg0_573 = buffer.data(nsg0 + 573);
    const auto *nsg0_575 = buffer.data(nsg0 + 575);
    const auto *nsg0_576 = buffer.data(nsg0 + 576);
    const auto *nsg0_579 = buffer.data(nsg0 + 579);
    const auto *nsg0_580 = buffer.data(nsg0 + 580);
    const auto *nsg0_582 = buffer.data(nsg0 + 582);
    const auto *nsg0_583 = buffer.data(nsg0 + 583);
    const auto *nsg0_584 = buffer.data(nsg0 + 584);
    const auto *nsg0_585 = buffer.data(nsg0 + 585);
    const auto *nsg0_588 = buffer.data(nsg0 + 588);
    const auto *nsg0_590 = buffer.data(nsg0 + 590);
    const auto *nsg0_591 = buffer.data(nsg0 + 591);
    const auto *nsg0_594 = buffer.data(nsg0 + 594);
    const auto *nsg0_595 = buffer.data(nsg0 + 595);
    const auto *nsg0_597 = buffer.data(nsg0 + 597);
    const auto *nsg0_598 = buffer.data(nsg0 + 598);
    const auto *nsg0_599 = buffer.data(nsg0 + 599);
    const auto *nsg0_600 = buffer.data(nsg0 + 600);
    const auto *nsg0_603 = buffer.data(nsg0 + 603);
    const auto *nsg0_605 = buffer.data(nsg0 + 605);
    const auto *nsg0_606 = buffer.data(nsg0 + 606);
    const auto *nsg0_609 = buffer.data(nsg0 + 609);
    const auto *nsg0_610 = buffer.data(nsg0 + 610);
    const auto *nsg0_612 = buffer.data(nsg0 + 612);
    const auto *nsg0_613 = buffer.data(nsg0 + 613);
    const auto *nsg0_614 = buffer.data(nsg0 + 614);
    const auto *nsg0_615 = buffer.data(nsg0 + 615);

    const auto *nsg1_564 = buffer.data(nsg1 + 564);
    const auto *nsg1_567 = buffer.data(nsg1 + 567);
    const auto *nsg1_568 = buffer.data(nsg1 + 568);
    const auto *nsg1_569 = buffer.data(nsg1 + 569);
    const auto *nsg1_570 = buffer.data(nsg1 + 570);
    const auto *nsg1_573 = buffer.data(nsg1 + 573);
    const auto *nsg1_575 = buffer.data(nsg1 + 575);
    const auto *nsg1_576 = buffer.data(nsg1 + 576);
    const auto *nsg1_579 = buffer.data(nsg1 + 579);
    const auto *nsg1_580 = buffer.data(nsg1 + 580);
    const auto *nsg1_582 = buffer.data(nsg1 + 582);
    const auto *nsg1_583 = buffer.data(nsg1 + 583);
    const auto *nsg1_584 = buffer.data(nsg1 + 584);
    const auto *nsg1_585 = buffer.data(nsg1 + 585);
    const auto *nsg1_588 = buffer.data(nsg1 + 588);
    const auto *nsg1_590 = buffer.data(nsg1 + 590);
    const auto *nsg1_591 = buffer.data(nsg1 + 591);
    const auto *nsg1_594 = buffer.data(nsg1 + 594);
    const auto *nsg1_595 = buffer.data(nsg1 + 595);
    const auto *nsg1_597 = buffer.data(nsg1 + 597);
    const auto *nsg1_598 = buffer.data(nsg1 + 598);
    const auto *nsg1_599 = buffer.data(nsg1 + 599);
    const auto *nsg1_600 = buffer.data(nsg1 + 600);
    const auto *nsg1_603 = buffer.data(nsg1 + 603);
    const auto *nsg1_605 = buffer.data(nsg1 + 605);
    const auto *nsg1_606 = buffer.data(nsg1 + 606);
    const auto *nsg1_609 = buffer.data(nsg1 + 609);
    const auto *nsg1_610 = buffer.data(nsg1 + 610);
    const auto *nsg1_612 = buffer.data(nsg1 + 612);
    const auto *nsg1_613 = buffer.data(nsg1 + 613);
    const auto *nsg1_614 = buffer.data(nsg1 + 614);
    const auto *nsg1_615 = buffer.data(nsg1 + 615);

    const auto *nsh_782 = buffer.data(nsh + 782);
    const auto *nsh_783 = buffer.data(nsh + 783);
    const auto *nsh_786 = buffer.data(nsh + 786);
    const auto *nsh_791 = buffer.data(nsh + 791);
    const auto *nsh_792 = buffer.data(nsh + 792);
    const auto *nsh_793 = buffer.data(nsh + 793);
    const auto *nsh_794 = buffer.data(nsh + 794);
    const auto *nsh_795 = buffer.data(nsh + 795);
    const auto *nsh_796 = buffer.data(nsh + 796);
    const auto *nsh_797 = buffer.data(nsh + 797);
    const auto *nsh_798 = buffer.data(nsh + 798);
    const auto *nsh_800 = buffer.data(nsh + 800);
    const auto *nsh_801 = buffer.data(nsh + 801);
    const auto *nsh_803 = buffer.data(nsh + 803);
    const auto *nsh_804 = buffer.data(nsh + 804);
    const auto *nsh_807 = buffer.data(nsh + 807);
    const auto *nsh_808 = buffer.data(nsh + 808);
    const auto *nsh_810 = buffer.data(nsh + 810);
    const auto *nsh_812 = buffer.data(nsh + 812);
    const auto *nsh_813 = buffer.data(nsh + 813);
    const auto *nsh_814 = buffer.data(nsh + 814);
    const auto *nsh_815 = buffer.data(nsh + 815);
    const auto *nsh_816 = buffer.data(nsh + 816);
    const auto *nsh_817 = buffer.data(nsh + 817);
    const auto *nsh_818 = buffer.data(nsh + 818);
    const auto *nsh_819 = buffer.data(nsh + 819);
    const auto *nsh_821 = buffer.data(nsh + 821);
    const auto *nsh_822 = buffer.data(nsh + 822);
    const auto *nsh_824 = buffer.data(nsh + 824);
    const auto *nsh_825 = buffer.data(nsh + 825);
    const auto *nsh_828 = buffer.data(nsh + 828);
    const auto *nsh_829 = buffer.data(nsh + 829);
    const auto *nsh_831 = buffer.data(nsh + 831);
    const auto *nsh_833 = buffer.data(nsh + 833);
    const auto *nsh_834 = buffer.data(nsh + 834);
    const auto *nsh_835 = buffer.data(nsh + 835);
    const auto *nsh_836 = buffer.data(nsh + 836);
    const auto *nsh_837 = buffer.data(nsh + 837);
    const auto *nsh_838 = buffer.data(nsh + 838);
    const auto *nsh_839 = buffer.data(nsh + 839);
    const auto *nsh_840 = buffer.data(nsh + 840);
    const auto *nsh_842 = buffer.data(nsh + 842);
    const auto *nsh_843 = buffer.data(nsh + 843);
    const auto *nsh_845 = buffer.data(nsh + 845);
    const auto *nsh_846 = buffer.data(nsh + 846);
    const auto *nsh_849 = buffer.data(nsh + 849);
    const auto *nsh_850 = buffer.data(nsh + 850);
    const auto *nsh_852 = buffer.data(nsh + 852);
    const auto *nsh_854 = buffer.data(nsh + 854);
    const auto *nsh_855 = buffer.data(nsh + 855);
    const auto *nsh_856 = buffer.data(nsh + 856);
    const auto *nsh_857 = buffer.data(nsh + 857);
    const auto *nsh_858 = buffer.data(nsh + 858);
    const auto *nsh_859 = buffer.data(nsh + 859);
    const auto *nsh_860 = buffer.data(nsh + 860);
    const auto *nsh_861 = buffer.data(nsh + 861);

#pragma omp simd aligned(t_1044, t_1045, t_1046, pa_z, pc_x, pc_y, pc_z, msi0_794, msh_614, \
                         msh_786, msi1_794, nsg0_564, nsg1_564, nsh_782, \
                         nsh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_19 * msh_614[k]
                    + f_3 * pc_y[k] * nsh_782[k];

        t_1045[k] = f_12 * msh_786[k]
                    + f_6 * nsg0_564[k]
                    - f_7 * nsg1_564[k]
                    + f_3 * pc_x[k] * nsh_786[k];

        t_1046[k] = pa_z[k] * msi0_794[k]
                    - f_10 * pc_z[k] * msi1_794[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pa_z, pc_y, pc_z, msi0_796, msh_594, msh_595, \
                         msh_618, msi1_796, nsh_783, nsh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = f_11 * msh_594[k]
                    + f_3 * pc_z[k] * nsh_783[k];

        t_1048[k] = pa_z[k] * msi0_796[k]
                    + f_12 * msh_595[k]
                    - f_10 * pc_z[k] * msi1_796[k];

        t_1049[k] = f_19 * msh_618[k]
                    + f_3 * pc_y[k] * nsh_786[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pc_x, msh_791, msh_792, msh_793, \
                         msh_794, nsg0_569, nsg1_569, nsh_791, nsh_792, nsh_793, \
                         nsh_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_12 * msh_791[k]
                    + f_4 * nsg0_569[k]
                    - f_5 * nsg1_569[k]
                    + f_3 * pc_x[k] * nsh_791[k];

        t_1051[k] = f_12 * msh_792[k]
                    + f_3 * pc_x[k] * nsh_792[k];

        t_1052[k] = f_12 * msh_793[k]
                    + f_3 * pc_x[k] * nsh_793[k];

        t_1053[k] = f_12 * msh_794[k]
                    + f_3 * pc_x[k] * nsh_794[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, pa_z, pc_x, pc_z, msi0_805, msh_795, \
                         msh_796, msh_797, msi1_805, nsh_795, nsh_796, \
                         nsh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_12 * msh_795[k]
                    + f_3 * pc_x[k] * nsh_795[k];

        t_1055[k] = f_12 * msh_796[k]
                    + f_3 * pc_x[k] * nsh_796[k];

        t_1056[k] = f_12 * msh_797[k]
                    + f_3 * pc_x[k] * nsh_797[k];

        t_1057[k] = pa_z[k] * msi0_805[k]
                    - f_10 * pc_z[k] * msi1_805[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_y, pc_z, msh_603, msh_626, msh_627, \
                         nsg0_567, nsg0_568, nsg1_567, nsg1_568, nsh_792, nsh_794, \
                         nsh_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_11 * msh_603[k]
                    + f_3 * pc_z[k] * nsh_792[k];

        t_1059[k] = f_19 * msh_626[k]
                    + f_8 * nsg0_567[k]
                    - f_9 * nsg1_567[k]
                    + f_3 * pc_y[k] * nsh_794[k];

        t_1060[k] = f_19 * msh_627[k]
                    + f_6 * nsg0_568[k]
                    - f_7 * nsg1_568[k]
                    + f_3 * pc_y[k] * nsh_795[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pc_y, pc_z, msh_608, msh_628, msh_629, \
                         nsg0_569, nsg1_569, nsh_796, nsh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_19 * msh_628[k]
                    + f_4 * nsg0_569[k]
                    - f_5 * nsg1_569[k]
                    + f_3 * pc_y[k] * nsh_796[k];

        t_1062[k] = f_19 * msh_629[k]
                    + f_3 * pc_y[k] * nsh_797[k];

        t_1063[k] = f_11 * msh_608[k]
                    + f_1 * nsg0_569[k]
                    - f_2 * nsg1_569[k]
                    + f_3 * pc_z[k] * nsh_797[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, pc_x, pc_y, pc_z, msh_609, msh_630, msh_798, \
                         nsg0_570, nsg1_570, nsh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_12 * msh_798[k]
                    + f_1 * nsg0_570[k]
                    - f_2 * nsg1_570[k]
                    + f_3 * pc_x[k] * nsh_798[k];

        t_1065[k] = f_20 * msh_630[k]
                    + f_3 * pc_y[k] * nsh_798[k];

        t_1066[k] = f_12 * msh_609[k]
                    + f_3 * pc_z[k] * nsh_798[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pc_x, pc_y, msh_632, msh_801, msh_803, \
                         nsg0_573, nsg0_575, nsg1_573, nsg1_575, nsh_800, nsh_801, \
                         nsh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_12 * msh_801[k]
                    + f_8 * nsg0_573[k]
                    - f_9 * nsg1_573[k]
                    + f_3 * pc_x[k] * nsh_801[k];

        t_1068[k] = f_20 * msh_632[k]
                    + f_3 * pc_y[k] * nsh_800[k];

        t_1069[k] = f_12 * msh_803[k]
                    + f_8 * nsg0_575[k]
                    - f_9 * nsg1_575[k]
                    + f_3 * pc_x[k] * nsh_803[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, msh_612, msh_635, msh_804, \
                         nsg0_576, nsg1_576, nsh_801, nsh_803, \
                         nsh_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_12 * msh_804[k]
                    + f_6 * nsg0_576[k]
                    - f_7 * nsg1_576[k]
                    + f_3 * pc_x[k] * nsh_804[k];

        t_1071[k] = f_12 * msh_612[k]
                    + f_3 * pc_z[k] * nsh_801[k];

        t_1072[k] = f_20 * msh_635[k]
                    + f_3 * pc_y[k] * nsh_803[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, pc_z, msh_615, msh_807, msh_808, \
                         nsg0_579, nsg0_580, nsg1_579, nsg1_580, nsh_804, nsh_807, \
                         nsh_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_12 * msh_807[k]
                    + f_6 * nsg0_579[k]
                    - f_7 * nsg1_579[k]
                    + f_3 * pc_x[k] * nsh_807[k];

        t_1074[k] = f_12 * msh_808[k]
                    + f_4 * nsg0_580[k]
                    - f_5 * nsg1_580[k]
                    + f_3 * pc_x[k] * nsh_808[k];

        t_1075[k] = f_12 * msh_615[k]
                    + f_3 * pc_z[k] * nsh_804[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, pc_y, msh_639, msh_810, msh_812, \
                         nsg0_582, nsg0_584, nsg1_582, nsg1_584, nsh_807, nsh_810, \
                         nsh_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_12 * msh_810[k]
                    + f_4 * nsg0_582[k]
                    - f_5 * nsg1_582[k]
                    + f_3 * pc_x[k] * nsh_810[k];

        t_1077[k] = f_20 * msh_639[k]
                    + f_3 * pc_y[k] * nsh_807[k];

        t_1078[k] = f_12 * msh_812[k]
                    + f_4 * nsg0_584[k]
                    - f_5 * nsg1_584[k]
                    + f_3 * pc_x[k] * nsh_812[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, t_1083, pc_x, msh_813, msh_814, \
                         msh_815, msh_816, msh_817, nsh_813, nsh_814, nsh_815, nsh_816, \
                         nsh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_12 * msh_813[k]
                    + f_3 * pc_x[k] * nsh_813[k];

        t_1080[k] = f_12 * msh_814[k]
                    + f_3 * pc_x[k] * nsh_814[k];

        t_1081[k] = f_12 * msh_815[k]
                    + f_3 * pc_x[k] * nsh_815[k];

        t_1082[k] = f_12 * msh_816[k]
                    + f_3 * pc_x[k] * nsh_816[k];

        t_1083[k] = f_12 * msh_817[k]
                    + f_3 * pc_x[k] * nsh_817[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, pc_x, pc_y, pc_z, msh_624, msh_645, msh_818, \
                         nsg0_580, nsg1_580, nsh_813, nsh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_12 * msh_818[k]
                    + f_3 * pc_x[k] * nsh_818[k];

        t_1085[k] = f_20 * msh_645[k]
                    + f_1 * nsg0_580[k]
                    - f_2 * nsg1_580[k]
                    + f_3 * pc_y[k] * nsh_813[k];

        t_1086[k] = f_12 * msh_624[k]
                    + f_3 * pc_z[k] * nsh_813[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, msh_647, msh_648, msh_649, nsg0_582, \
                         nsg0_583, nsg0_584, nsg1_582, nsg1_583, nsg1_584, nsh_815, nsh_816, \
                         nsh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_20 * msh_647[k]
                    + f_8 * nsg0_582[k]
                    - f_9 * nsg1_582[k]
                    + f_3 * pc_y[k] * nsh_815[k];

        t_1088[k] = f_20 * msh_648[k]
                    + f_6 * nsg0_583[k]
                    - f_7 * nsg1_583[k]
                    + f_3 * pc_y[k] * nsh_816[k];

        t_1089[k] = f_20 * msh_649[k]
                    + f_4 * nsg0_584[k]
                    - f_5 * nsg1_584[k]
                    + f_3 * pc_y[k] * nsh_817[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, msh_629, msh_650, msh_819, \
                         nsg0_584, nsg0_585, nsg1_584, nsg1_585, nsh_818, \
                         nsh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_20 * msh_650[k]
                    + f_3 * pc_y[k] * nsh_818[k];

        t_1091[k] = f_12 * msh_629[k]
                    + f_1 * nsg0_584[k]
                    - f_2 * nsg1_584[k]
                    + f_3 * pc_z[k] * nsh_818[k];

        t_1092[k] = f_12 * msh_819[k]
                    + f_1 * nsg0_585[k]
                    - f_2 * nsg1_585[k]
                    + f_3 * pc_x[k] * nsh_819[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pc_x, pc_y, pc_z, msh_630, msh_651, \
                         msh_653, msh_822, nsg0_588, nsg1_588, nsh_819, nsh_821, \
                         nsh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_21 * msh_651[k]
                    + f_3 * pc_y[k] * nsh_819[k];

        t_1094[k] = f_13 * msh_630[k]
                    + f_3 * pc_z[k] * nsh_819[k];

        t_1095[k] = f_12 * msh_822[k]
                    + f_8 * nsg0_588[k]
                    - f_9 * nsg1_588[k]
                    + f_3 * pc_x[k] * nsh_822[k];

        t_1096[k] = f_21 * msh_653[k]
                    + f_3 * pc_y[k] * nsh_821[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_z, msh_633, msh_824, msh_825, \
                         nsg0_590, nsg0_591, nsg1_590, nsg1_591, nsh_822, nsh_824, \
                         nsh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_12 * msh_824[k]
                    + f_8 * nsg0_590[k]
                    - f_9 * nsg1_590[k]
                    + f_3 * pc_x[k] * nsh_824[k];

        t_1098[k] = f_12 * msh_825[k]
                    + f_6 * nsg0_591[k]
                    - f_7 * nsg1_591[k]
                    + f_3 * pc_x[k] * nsh_825[k];

        t_1099[k] = f_13 * msh_633[k]
                    + f_3 * pc_z[k] * nsh_822[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_y, msh_656, msh_828, msh_829, \
                         nsg0_594, nsg0_595, nsg1_594, nsg1_595, nsh_824, nsh_828, \
                         nsh_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_21 * msh_656[k]
                    + f_3 * pc_y[k] * nsh_824[k];

        t_1101[k] = f_12 * msh_828[k]
                    + f_6 * nsg0_594[k]
                    - f_7 * nsg1_594[k]
                    + f_3 * pc_x[k] * nsh_828[k];

        t_1102[k] = f_12 * msh_829[k]
                    + f_4 * nsg0_595[k]
                    - f_5 * nsg1_595[k]
                    + f_3 * pc_x[k] * nsh_829[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, pc_y, pc_z, msh_636, msh_660, msh_831, \
                         nsg0_597, nsg1_597, nsh_825, nsh_828, \
                         nsh_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_13 * msh_636[k]
                    + f_3 * pc_z[k] * nsh_825[k];

        t_1104[k] = f_12 * msh_831[k]
                    + f_4 * nsg0_597[k]
                    - f_5 * nsg1_597[k]
                    + f_3 * pc_x[k] * nsh_831[k];

        t_1105[k] = f_21 * msh_660[k]
                    + f_3 * pc_y[k] * nsh_828[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, msh_833, msh_834, msh_835, \
                         msh_836, nsg0_599, nsg1_599, nsh_833, nsh_834, nsh_835, \
                         nsh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_12 * msh_833[k]
                    + f_4 * nsg0_599[k]
                    - f_5 * nsg1_599[k]
                    + f_3 * pc_x[k] * nsh_833[k];

        t_1107[k] = f_12 * msh_834[k]
                    + f_3 * pc_x[k] * nsh_834[k];

        t_1108[k] = f_12 * msh_835[k]
                    + f_3 * pc_x[k] * nsh_835[k];

        t_1109[k] = f_12 * msh_836[k]
                    + f_3 * pc_x[k] * nsh_836[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pc_x, pc_y, msh_666, msh_837, \
                         msh_838, msh_839, nsg0_595, nsg1_595, nsh_834, nsh_837, nsh_838, \
                         nsh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_12 * msh_837[k]
                    + f_3 * pc_x[k] * nsh_837[k];

        t_1111[k] = f_12 * msh_838[k]
                    + f_3 * pc_x[k] * nsh_838[k];

        t_1112[k] = f_12 * msh_839[k]
                    + f_3 * pc_x[k] * nsh_839[k];

        t_1113[k] = f_21 * msh_666[k]
                    + f_1 * nsg0_595[k]
                    - f_2 * nsg1_595[k]
                    + f_3 * pc_y[k] * nsh_834[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_y, pc_z, msh_645, msh_668, msh_669, \
                         nsg0_597, nsg0_598, nsg1_597, nsg1_598, nsh_834, nsh_836, \
                         nsh_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_13 * msh_645[k]
                    + f_3 * pc_z[k] * nsh_834[k];

        t_1115[k] = f_21 * msh_668[k]
                    + f_8 * nsg0_597[k]
                    - f_9 * nsg1_597[k]
                    + f_3 * pc_y[k] * nsh_836[k];

        t_1116[k] = f_21 * msh_669[k]
                    + f_6 * nsg0_598[k]
                    - f_7 * nsg1_598[k]
                    + f_3 * pc_y[k] * nsh_837[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, pc_y, pc_z, msh_650, msh_670, msh_671, \
                         nsg0_599, nsg1_599, nsh_838, nsh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_21 * msh_670[k]
                    + f_4 * nsg0_599[k]
                    - f_5 * nsg1_599[k]
                    + f_3 * pc_y[k] * nsh_838[k];

        t_1118[k] = f_21 * msh_671[k]
                    + f_3 * pc_y[k] * nsh_839[k];

        t_1119[k] = f_13 * msh_650[k]
                    + f_1 * nsg0_599[k]
                    - f_2 * nsg1_599[k]
                    + f_3 * pc_z[k] * nsh_839[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, pc_x, pc_y, pc_z, msh_651, msh_672, msh_840, \
                         nsg0_600, nsg1_600, nsh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = f_12 * msh_840[k]
                    + f_1 * nsg0_600[k]
                    - f_2 * nsg1_600[k]
                    + f_3 * pc_x[k] * nsh_840[k];

        t_1121[k] = f_14 * msh_672[k]
                    + f_3 * pc_y[k] * nsh_840[k];

        t_1122[k] = f_14 * msh_651[k]
                    + f_3 * pc_z[k] * nsh_840[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, pc_x, pc_y, msh_674, msh_843, msh_845, \
                         nsg0_603, nsg0_605, nsg1_603, nsg1_605, nsh_842, nsh_843, \
                         nsh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_12 * msh_843[k]
                    + f_8 * nsg0_603[k]
                    - f_9 * nsg1_603[k]
                    + f_3 * pc_x[k] * nsh_843[k];

        t_1124[k] = f_14 * msh_674[k]
                    + f_3 * pc_y[k] * nsh_842[k];

        t_1125[k] = f_12 * msh_845[k]
                    + f_8 * nsg0_605[k]
                    - f_9 * nsg1_605[k]
                    + f_3 * pc_x[k] * nsh_845[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, pc_x, pc_y, pc_z, msh_654, msh_677, msh_846, \
                         nsg0_606, nsg1_606, nsh_843, nsh_845, \
                         nsh_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_12 * msh_846[k]
                    + f_6 * nsg0_606[k]
                    - f_7 * nsg1_606[k]
                    + f_3 * pc_x[k] * nsh_846[k];

        t_1127[k] = f_14 * msh_654[k]
                    + f_3 * pc_z[k] * nsh_843[k];

        t_1128[k] = f_14 * msh_677[k]
                    + f_3 * pc_y[k] * nsh_845[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pc_x, pc_z, msh_657, msh_849, msh_850, \
                         nsg0_609, nsg0_610, nsg1_609, nsg1_610, nsh_846, nsh_849, \
                         nsh_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_12 * msh_849[k]
                    + f_6 * nsg0_609[k]
                    - f_7 * nsg1_609[k]
                    + f_3 * pc_x[k] * nsh_849[k];

        t_1130[k] = f_12 * msh_850[k]
                    + f_4 * nsg0_610[k]
                    - f_5 * nsg1_610[k]
                    + f_3 * pc_x[k] * nsh_850[k];

        t_1131[k] = f_14 * msh_657[k]
                    + f_3 * pc_z[k] * nsh_846[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, pc_x, pc_y, msh_681, msh_852, msh_854, \
                         nsg0_612, nsg0_614, nsg1_612, nsg1_614, nsh_849, nsh_852, \
                         nsh_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_12 * msh_852[k]
                    + f_4 * nsg0_612[k]
                    - f_5 * nsg1_612[k]
                    + f_3 * pc_x[k] * nsh_852[k];

        t_1133[k] = f_14 * msh_681[k]
                    + f_3 * pc_y[k] * nsh_849[k];

        t_1134[k] = f_12 * msh_854[k]
                    + f_4 * nsg0_614[k]
                    - f_5 * nsg1_614[k]
                    + f_3 * pc_x[k] * nsh_854[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, t_1138, t_1139, pc_x, msh_855, msh_856, \
                         msh_857, msh_858, msh_859, nsh_855, nsh_856, nsh_857, nsh_858, \
                         nsh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = f_12 * msh_855[k]
                    + f_3 * pc_x[k] * nsh_855[k];

        t_1136[k] = f_12 * msh_856[k]
                    + f_3 * pc_x[k] * nsh_856[k];

        t_1137[k] = f_12 * msh_857[k]
                    + f_3 * pc_x[k] * nsh_857[k];

        t_1138[k] = f_12 * msh_858[k]
                    + f_3 * pc_x[k] * nsh_858[k];

        t_1139[k] = f_12 * msh_859[k]
                    + f_3 * pc_x[k] * nsh_859[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, pc_x, pc_y, pc_z, msh_666, msh_687, msh_860, \
                         nsg0_610, nsg1_610, nsh_855, nsh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_12 * msh_860[k]
                    + f_3 * pc_x[k] * nsh_860[k];

        t_1141[k] = f_14 * msh_687[k]
                    + f_1 * nsg0_610[k]
                    - f_2 * nsg1_610[k]
                    + f_3 * pc_y[k] * nsh_855[k];

        t_1142[k] = f_14 * msh_666[k]
                    + f_3 * pc_z[k] * nsh_855[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, pc_y, msh_689, msh_690, msh_691, nsg0_612, \
                         nsg0_613, nsg0_614, nsg1_612, nsg1_613, nsg1_614, nsh_857, nsh_858, \
                         nsh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_14 * msh_689[k]
                    + f_8 * nsg0_612[k]
                    - f_9 * nsg1_612[k]
                    + f_3 * pc_y[k] * nsh_857[k];

        t_1144[k] = f_14 * msh_690[k]
                    + f_6 * nsg0_613[k]
                    - f_7 * nsg1_613[k]
                    + f_3 * pc_y[k] * nsh_858[k];

        t_1145[k] = f_14 * msh_691[k]
                    + f_4 * nsg0_614[k]
                    - f_5 * nsg1_614[k]
                    + f_3 * pc_y[k] * nsh_859[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pc_x, pc_y, pc_z, msh_671, msh_692, msh_861, \
                         nsg0_614, nsg0_615, nsg1_614, nsg1_615, nsh_860, \
                         nsh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_14 * msh_692[k]
                    + f_3 * pc_y[k] * nsh_860[k];

        t_1147[k] = f_14 * msh_671[k]
                    + f_1 * nsg0_614[k]
                    - f_2 * nsg1_614[k]
                    + f_3 * pc_z[k] * nsh_860[k];

        t_1148[k] = f_12 * msh_861[k]
                    + f_1 * nsg0_615[k]
                    - f_2 * nsg1_615[k]
                    + f_3 * pc_x[k] * nsh_861[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msi0,
                                                           const size_t msh, const size_t msi1,
                                                           const size_t nsg0, const size_t nsg1,
                                                           const size_t nsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_980 = buffer.data(msi0 + 980);
    const auto *msi0_983 = buffer.data(msi0 + 983);
    const auto *msi0_985 = buffer.data(msi0 + 985);
    const auto *msi0_986 = buffer.data(msi0 + 986);
    const auto *msi0_989 = buffer.data(msi0 + 989);
    const auto *msi0_990 = buffer.data(msi0 + 990);
    const auto *msi0_992 = buffer.data(msi0 + 992);
    const auto *msi0_994 = buffer.data(msi0 + 994);
    const auto *msi0_1007 = buffer.data(msi0 + 1007);

    const auto *msh_672 = buffer.data(msh + 672);
    const auto *msh_675 = buffer.data(msh + 675);
    const auto *msh_678 = buffer.data(msh + 678);
    const auto *msh_687 = buffer.data(msh + 687);
    const auto *msh_692 = buffer.data(msh + 692);
    const auto *msh_693 = buffer.data(msh + 693);
    const auto *msh_695 = buffer.data(msh + 695);
    const auto *msh_696 = buffer.data(msh + 696);
    const auto *msh_698 = buffer.data(msh + 698);
    const auto *msh_699 = buffer.data(msh + 699);
    const auto *msh_702 = buffer.data(msh + 702);
    const auto *msh_708 = buffer.data(msh + 708);
    const auto *msh_710 = buffer.data(msh + 710);
    const auto *msh_711 = buffer.data(msh + 711);
    const auto *msh_712 = buffer.data(msh + 712);
    const auto *msh_713 = buffer.data(msh + 713);
    const auto *msh_714 = buffer.data(msh + 714);
    const auto *msh_716 = buffer.data(msh + 716);
    const auto *msh_717 = buffer.data(msh + 717);
    const auto *msh_719 = buffer.data(msh + 719);
    const auto *msh_720 = buffer.data(msh + 720);
    const auto *msh_723 = buffer.data(msh + 723);
    const auto *msh_729 = buffer.data(msh + 729);
    const auto *msh_731 = buffer.data(msh + 731);
    const auto *msh_732 = buffer.data(msh + 732);
    const auto *msh_733 = buffer.data(msh + 733);
    const auto *msh_734 = buffer.data(msh + 734);
    const auto *msh_735 = buffer.data(msh + 735);
    const auto *msh_736 = buffer.data(msh + 736);
    const auto *msh_737 = buffer.data(msh + 737);
    const auto *msh_738 = buffer.data(msh + 738);
    const auto *msh_740 = buffer.data(msh + 740);
    const auto *msh_741 = buffer.data(msh + 741);
    const auto *msh_743 = buffer.data(msh + 743);
    const auto *msh_744 = buffer.data(msh + 744);
    const auto *msh_750 = buffer.data(msh + 750);
    const auto *msh_752 = buffer.data(msh + 752);
    const auto *msh_753 = buffer.data(msh + 753);
    const auto *msh_754 = buffer.data(msh + 754);
    const auto *msh_755 = buffer.data(msh + 755);
    const auto *msh_864 = buffer.data(msh + 864);
    const auto *msh_866 = buffer.data(msh + 866);
    const auto *msh_867 = buffer.data(msh + 867);
    const auto *msh_870 = buffer.data(msh + 870);
    const auto *msh_871 = buffer.data(msh + 871);
    const auto *msh_873 = buffer.data(msh + 873);
    const auto *msh_875 = buffer.data(msh + 875);
    const auto *msh_876 = buffer.data(msh + 876);
    const auto *msh_877 = buffer.data(msh + 877);
    const auto *msh_878 = buffer.data(msh + 878);
    const auto *msh_879 = buffer.data(msh + 879);
    const auto *msh_880 = buffer.data(msh + 880);
    const auto *msh_881 = buffer.data(msh + 881);
    const auto *msh_882 = buffer.data(msh + 882);
    const auto *msh_885 = buffer.data(msh + 885);
    const auto *msh_887 = buffer.data(msh + 887);
    const auto *msh_888 = buffer.data(msh + 888);
    const auto *msh_891 = buffer.data(msh + 891);
    const auto *msh_892 = buffer.data(msh + 892);
    const auto *msh_894 = buffer.data(msh + 894);
    const auto *msh_896 = buffer.data(msh + 896);
    const auto *msh_897 = buffer.data(msh + 897);
    const auto *msh_898 = buffer.data(msh + 898);
    const auto *msh_899 = buffer.data(msh + 899);
    const auto *msh_900 = buffer.data(msh + 900);
    const auto *msh_901 = buffer.data(msh + 901);
    const auto *msh_902 = buffer.data(msh + 902);
    const auto *msh_918 = buffer.data(msh + 918);
    const auto *msh_919 = buffer.data(msh + 919);
    const auto *msh_920 = buffer.data(msh + 920);
    const auto *msh_921 = buffer.data(msh + 921);
    const auto *msh_922 = buffer.data(msh + 922);
    const auto *msh_923 = buffer.data(msh + 923);
    const auto *msh_924 = buffer.data(msh + 924);
    const auto *msh_929 = buffer.data(msh + 929);
    const auto *msh_933 = buffer.data(msh + 933);
    const auto *msh_938 = buffer.data(msh + 938);
    const auto *msh_939 = buffer.data(msh + 939);
    const auto *msh_940 = buffer.data(msh + 940);
    const auto *msh_941 = buffer.data(msh + 941);
    const auto *msh_942 = buffer.data(msh + 942);
    const auto *msh_944 = buffer.data(msh + 944);

    const auto *msi1_980 = buffer.data(msi1 + 980);
    const auto *msi1_983 = buffer.data(msi1 + 983);
    const auto *msi1_985 = buffer.data(msi1 + 985);
    const auto *msi1_986 = buffer.data(msi1 + 986);
    const auto *msi1_989 = buffer.data(msi1 + 989);
    const auto *msi1_990 = buffer.data(msi1 + 990);
    const auto *msi1_992 = buffer.data(msi1 + 992);
    const auto *msi1_994 = buffer.data(msi1 + 994);
    const auto *msi1_1007 = buffer.data(msi1 + 1007);

    const auto *nsg0_618 = buffer.data(nsg0 + 618);
    const auto *nsg0_620 = buffer.data(nsg0 + 620);
    const auto *nsg0_621 = buffer.data(nsg0 + 621);
    const auto *nsg0_624 = buffer.data(nsg0 + 624);
    const auto *nsg0_625 = buffer.data(nsg0 + 625);
    const auto *nsg0_627 = buffer.data(nsg0 + 627);
    const auto *nsg0_628 = buffer.data(nsg0 + 628);
    const auto *nsg0_629 = buffer.data(nsg0 + 629);
    const auto *nsg0_630 = buffer.data(nsg0 + 630);
    const auto *nsg0_633 = buffer.data(nsg0 + 633);
    const auto *nsg0_635 = buffer.data(nsg0 + 635);
    const auto *nsg0_636 = buffer.data(nsg0 + 636);
    const auto *nsg0_639 = buffer.data(nsg0 + 639);
    const auto *nsg0_640 = buffer.data(nsg0 + 640);
    const auto *nsg0_642 = buffer.data(nsg0 + 642);
    const auto *nsg0_643 = buffer.data(nsg0 + 643);
    const auto *nsg0_644 = buffer.data(nsg0 + 644);
    const auto *nsg0_655 = buffer.data(nsg0 + 655);
    const auto *nsg0_657 = buffer.data(nsg0 + 657);
    const auto *nsg0_658 = buffer.data(nsg0 + 658);
    const auto *nsg0_659 = buffer.data(nsg0 + 659);
    const auto *nsg0_660 = buffer.data(nsg0 + 660);
    const auto *nsg0_661 = buffer.data(nsg0 + 661);
    const auto *nsg0_662 = buffer.data(nsg0 + 662);
    const auto *nsg0_663 = buffer.data(nsg0 + 663);
    const auto *nsg0_664 = buffer.data(nsg0 + 664);
    const auto *nsg0_665 = buffer.data(nsg0 + 665);
    const auto *nsg0_669 = buffer.data(nsg0 + 669);
    const auto *nsg0_670 = buffer.data(nsg0 + 670);
    const auto *nsg0_671 = buffer.data(nsg0 + 671);
    const auto *nsg0_672 = buffer.data(nsg0 + 672);
    const auto *nsg0_673 = buffer.data(nsg0 + 673);
    const auto *nsg0_674 = buffer.data(nsg0 + 674);

    const auto *nsg1_618 = buffer.data(nsg1 + 618);
    const auto *nsg1_620 = buffer.data(nsg1 + 620);
    const auto *nsg1_621 = buffer.data(nsg1 + 621);
    const auto *nsg1_624 = buffer.data(nsg1 + 624);
    const auto *nsg1_625 = buffer.data(nsg1 + 625);
    const auto *nsg1_627 = buffer.data(nsg1 + 627);
    const auto *nsg1_628 = buffer.data(nsg1 + 628);
    const auto *nsg1_629 = buffer.data(nsg1 + 629);
    const auto *nsg1_630 = buffer.data(nsg1 + 630);
    const auto *nsg1_633 = buffer.data(nsg1 + 633);
    const auto *nsg1_635 = buffer.data(nsg1 + 635);
    const auto *nsg1_636 = buffer.data(nsg1 + 636);
    const auto *nsg1_639 = buffer.data(nsg1 + 639);
    const auto *nsg1_640 = buffer.data(nsg1 + 640);
    const auto *nsg1_642 = buffer.data(nsg1 + 642);
    const auto *nsg1_643 = buffer.data(nsg1 + 643);
    const auto *nsg1_644 = buffer.data(nsg1 + 644);
    const auto *nsg1_655 = buffer.data(nsg1 + 655);
    const auto *nsg1_657 = buffer.data(nsg1 + 657);
    const auto *nsg1_658 = buffer.data(nsg1 + 658);
    const auto *nsg1_659 = buffer.data(nsg1 + 659);
    const auto *nsg1_660 = buffer.data(nsg1 + 660);
    const auto *nsg1_661 = buffer.data(nsg1 + 661);
    const auto *nsg1_662 = buffer.data(nsg1 + 662);
    const auto *nsg1_663 = buffer.data(nsg1 + 663);
    const auto *nsg1_664 = buffer.data(nsg1 + 664);
    const auto *nsg1_665 = buffer.data(nsg1 + 665);
    const auto *nsg1_669 = buffer.data(nsg1 + 669);
    const auto *nsg1_670 = buffer.data(nsg1 + 670);
    const auto *nsg1_671 = buffer.data(nsg1 + 671);
    const auto *nsg1_672 = buffer.data(nsg1 + 672);
    const auto *nsg1_673 = buffer.data(nsg1 + 673);
    const auto *nsg1_674 = buffer.data(nsg1 + 674);

    const auto *nsh_861 = buffer.data(nsh + 861);
    const auto *nsh_863 = buffer.data(nsh + 863);
    const auto *nsh_864 = buffer.data(nsh + 864);
    const auto *nsh_866 = buffer.data(nsh + 866);
    const auto *nsh_867 = buffer.data(nsh + 867);
    const auto *nsh_870 = buffer.data(nsh + 870);
    const auto *nsh_871 = buffer.data(nsh + 871);
    const auto *nsh_873 = buffer.data(nsh + 873);
    const auto *nsh_875 = buffer.data(nsh + 875);
    const auto *nsh_876 = buffer.data(nsh + 876);
    const auto *nsh_877 = buffer.data(nsh + 877);
    const auto *nsh_878 = buffer.data(nsh + 878);
    const auto *nsh_879 = buffer.data(nsh + 879);
    const auto *nsh_880 = buffer.data(nsh + 880);
    const auto *nsh_881 = buffer.data(nsh + 881);
    const auto *nsh_882 = buffer.data(nsh + 882);
    const auto *nsh_884 = buffer.data(nsh + 884);
    const auto *nsh_885 = buffer.data(nsh + 885);
    const auto *nsh_887 = buffer.data(nsh + 887);
    const auto *nsh_888 = buffer.data(nsh + 888);
    const auto *nsh_891 = buffer.data(nsh + 891);
    const auto *nsh_892 = buffer.data(nsh + 892);
    const auto *nsh_894 = buffer.data(nsh + 894);
    const auto *nsh_896 = buffer.data(nsh + 896);
    const auto *nsh_897 = buffer.data(nsh + 897);
    const auto *nsh_898 = buffer.data(nsh + 898);
    const auto *nsh_899 = buffer.data(nsh + 899);
    const auto *nsh_900 = buffer.data(nsh + 900);
    const auto *nsh_901 = buffer.data(nsh + 901);
    const auto *nsh_902 = buffer.data(nsh + 902);
    const auto *nsh_903 = buffer.data(nsh + 903);
    const auto *nsh_905 = buffer.data(nsh + 905);
    const auto *nsh_906 = buffer.data(nsh + 906);
    const auto *nsh_908 = buffer.data(nsh + 908);
    const auto *nsh_909 = buffer.data(nsh + 909);
    const auto *nsh_912 = buffer.data(nsh + 912);
    const auto *nsh_918 = buffer.data(nsh + 918);
    const auto *nsh_919 = buffer.data(nsh + 919);
    const auto *nsh_920 = buffer.data(nsh + 920);
    const auto *nsh_921 = buffer.data(nsh + 921);
    const auto *nsh_922 = buffer.data(nsh + 922);
    const auto *nsh_923 = buffer.data(nsh + 923);
    const auto *nsh_924 = buffer.data(nsh + 924);
    const auto *nsh_925 = buffer.data(nsh + 925);
    const auto *nsh_926 = buffer.data(nsh + 926);
    const auto *nsh_927 = buffer.data(nsh + 927);
    const auto *nsh_928 = buffer.data(nsh + 928);
    const auto *nsh_929 = buffer.data(nsh + 929);
    const auto *nsh_930 = buffer.data(nsh + 930);
    const auto *nsh_931 = buffer.data(nsh + 931);
    const auto *nsh_932 = buffer.data(nsh + 932);
    const auto *nsh_933 = buffer.data(nsh + 933);
    const auto *nsh_938 = buffer.data(nsh + 938);
    const auto *nsh_939 = buffer.data(nsh + 939);
    const auto *nsh_940 = buffer.data(nsh + 940);
    const auto *nsh_941 = buffer.data(nsh + 941);
    const auto *nsh_942 = buffer.data(nsh + 942);
    const auto *nsh_943 = buffer.data(nsh + 943);
    const auto *nsh_944 = buffer.data(nsh + 944);

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, msh_672, msh_693, \
                         msh_695, msh_864, nsg0_618, nsg1_618, nsh_861, nsh_863, \
                         nsh_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_13 * msh_693[k]
                    + f_3 * pc_y[k] * nsh_861[k];

        t_1150[k] = f_21 * msh_672[k]
                    + f_3 * pc_z[k] * nsh_861[k];

        t_1151[k] = f_12 * msh_864[k]
                    + f_8 * nsg0_618[k]
                    - f_9 * nsg1_618[k]
                    + f_3 * pc_x[k] * nsh_864[k];

        t_1152[k] = f_13 * msh_695[k]
                    + f_3 * pc_y[k] * nsh_863[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, pc_x, pc_z, msh_675, msh_866, msh_867, \
                         nsg0_620, nsg0_621, nsg1_620, nsg1_621, nsh_864, nsh_866, \
                         nsh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_12 * msh_866[k]
                    + f_8 * nsg0_620[k]
                    - f_9 * nsg1_620[k]
                    + f_3 * pc_x[k] * nsh_866[k];

        t_1154[k] = f_12 * msh_867[k]
                    + f_6 * nsg0_621[k]
                    - f_7 * nsg1_621[k]
                    + f_3 * pc_x[k] * nsh_867[k];

        t_1155[k] = f_21 * msh_675[k]
                    + f_3 * pc_z[k] * nsh_864[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, pc_x, pc_y, msh_698, msh_870, msh_871, \
                         nsg0_624, nsg0_625, nsg1_624, nsg1_625, nsh_866, nsh_870, \
                         nsh_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_13 * msh_698[k]
                    + f_3 * pc_y[k] * nsh_866[k];

        t_1157[k] = f_12 * msh_870[k]
                    + f_6 * nsg0_624[k]
                    - f_7 * nsg1_624[k]
                    + f_3 * pc_x[k] * nsh_870[k];

        t_1158[k] = f_12 * msh_871[k]
                    + f_4 * nsg0_625[k]
                    - f_5 * nsg1_625[k]
                    + f_3 * pc_x[k] * nsh_871[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, pc_x, pc_y, pc_z, msh_678, msh_702, msh_873, \
                         nsg0_627, nsg1_627, nsh_867, nsh_870, \
                         nsh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_21 * msh_678[k]
                    + f_3 * pc_z[k] * nsh_867[k];

        t_1160[k] = f_12 * msh_873[k]
                    + f_4 * nsg0_627[k]
                    - f_5 * nsg1_627[k]
                    + f_3 * pc_x[k] * nsh_873[k];

        t_1161[k] = f_13 * msh_702[k]
                    + f_3 * pc_y[k] * nsh_870[k];
    }

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, pc_x, msh_875, msh_876, msh_877, \
                         msh_878, nsg0_629, nsg1_629, nsh_875, nsh_876, nsh_877, \
                         nsh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = f_12 * msh_875[k]
                    + f_4 * nsg0_629[k]
                    - f_5 * nsg1_629[k]
                    + f_3 * pc_x[k] * nsh_875[k];

        t_1163[k] = f_12 * msh_876[k]
                    + f_3 * pc_x[k] * nsh_876[k];

        t_1164[k] = f_12 * msh_877[k]
                    + f_3 * pc_x[k] * nsh_877[k];

        t_1165[k] = f_12 * msh_878[k]
                    + f_3 * pc_x[k] * nsh_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, pc_x, pc_y, msh_708, msh_879, \
                         msh_880, msh_881, nsg0_625, nsg1_625, nsh_876, nsh_879, nsh_880, \
                         nsh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_12 * msh_879[k]
                    + f_3 * pc_x[k] * nsh_879[k];

        t_1167[k] = f_12 * msh_880[k]
                    + f_3 * pc_x[k] * nsh_880[k];

        t_1168[k] = f_12 * msh_881[k]
                    + f_3 * pc_x[k] * nsh_881[k];

        t_1169[k] = f_13 * msh_708[k]
                    + f_1 * nsg0_625[k]
                    - f_2 * nsg1_625[k]
                    + f_3 * pc_y[k] * nsh_876[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, pc_y, pc_z, msh_687, msh_710, msh_711, \
                         nsg0_627, nsg0_628, nsg1_627, nsg1_628, nsh_876, nsh_878, \
                         nsh_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_21 * msh_687[k]
                    + f_3 * pc_z[k] * nsh_876[k];

        t_1171[k] = f_13 * msh_710[k]
                    + f_8 * nsg0_627[k]
                    - f_9 * nsg1_627[k]
                    + f_3 * pc_y[k] * nsh_878[k];

        t_1172[k] = f_13 * msh_711[k]
                    + f_6 * nsg0_628[k]
                    - f_7 * nsg1_628[k]
                    + f_3 * pc_y[k] * nsh_879[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, pc_y, pc_z, msh_692, msh_712, msh_713, \
                         nsg0_629, nsg1_629, nsh_880, nsh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_13 * msh_712[k]
                    + f_4 * nsg0_629[k]
                    - f_5 * nsg1_629[k]
                    + f_3 * pc_y[k] * nsh_880[k];

        t_1174[k] = f_13 * msh_713[k]
                    + f_3 * pc_y[k] * nsh_881[k];

        t_1175[k] = f_21 * msh_692[k]
                    + f_1 * nsg0_629[k]
                    - f_2 * nsg1_629[k]
                    + f_3 * pc_z[k] * nsh_881[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, pc_x, pc_y, pc_z, msh_693, msh_714, msh_882, \
                         nsg0_630, nsg1_630, nsh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_12 * msh_882[k]
                    + f_1 * nsg0_630[k]
                    - f_2 * nsg1_630[k]
                    + f_3 * pc_x[k] * nsh_882[k];

        t_1177[k] = f_12 * msh_714[k]
                    + f_3 * pc_y[k] * nsh_882[k];

        t_1178[k] = f_20 * msh_693[k]
                    + f_3 * pc_z[k] * nsh_882[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pc_x, pc_y, msh_716, msh_885, msh_887, \
                         nsg0_633, nsg0_635, nsg1_633, nsg1_635, nsh_884, nsh_885, \
                         nsh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = f_12 * msh_885[k]
                    + f_8 * nsg0_633[k]
                    - f_9 * nsg1_633[k]
                    + f_3 * pc_x[k] * nsh_885[k];

        t_1180[k] = f_12 * msh_716[k]
                    + f_3 * pc_y[k] * nsh_884[k];

        t_1181[k] = f_12 * msh_887[k]
                    + f_8 * nsg0_635[k]
                    - f_9 * nsg1_635[k]
                    + f_3 * pc_x[k] * nsh_887[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pc_x, pc_y, pc_z, msh_696, msh_719, msh_888, \
                         nsg0_636, nsg1_636, nsh_885, nsh_887, \
                         nsh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_12 * msh_888[k]
                    + f_6 * nsg0_636[k]
                    - f_7 * nsg1_636[k]
                    + f_3 * pc_x[k] * nsh_888[k];

        t_1183[k] = f_20 * msh_696[k]
                    + f_3 * pc_z[k] * nsh_885[k];

        t_1184[k] = f_12 * msh_719[k]
                    + f_3 * pc_y[k] * nsh_887[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_x, pc_z, msh_699, msh_891, msh_892, \
                         nsg0_639, nsg0_640, nsg1_639, nsg1_640, nsh_888, nsh_891, \
                         nsh_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_12 * msh_891[k]
                    + f_6 * nsg0_639[k]
                    - f_7 * nsg1_639[k]
                    + f_3 * pc_x[k] * nsh_891[k];

        t_1186[k] = f_12 * msh_892[k]
                    + f_4 * nsg0_640[k]
                    - f_5 * nsg1_640[k]
                    + f_3 * pc_x[k] * nsh_892[k];

        t_1187[k] = f_20 * msh_699[k]
                    + f_3 * pc_z[k] * nsh_888[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, pc_x, pc_y, msh_723, msh_894, msh_896, \
                         nsg0_642, nsg0_644, nsg1_642, nsg1_644, nsh_891, nsh_894, \
                         nsh_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_12 * msh_894[k]
                    + f_4 * nsg0_642[k]
                    - f_5 * nsg1_642[k]
                    + f_3 * pc_x[k] * nsh_894[k];

        t_1189[k] = f_12 * msh_723[k]
                    + f_3 * pc_y[k] * nsh_891[k];

        t_1190[k] = f_12 * msh_896[k]
                    + f_4 * nsg0_644[k]
                    - f_5 * nsg1_644[k]
                    + f_3 * pc_x[k] * nsh_896[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, t_1194, t_1195, pc_x, msh_897, msh_898, \
                         msh_899, msh_900, msh_901, nsh_897, nsh_898, nsh_899, nsh_900, \
                         nsh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = f_12 * msh_897[k]
                    + f_3 * pc_x[k] * nsh_897[k];

        t_1192[k] = f_12 * msh_898[k]
                    + f_3 * pc_x[k] * nsh_898[k];

        t_1193[k] = f_12 * msh_899[k]
                    + f_3 * pc_x[k] * nsh_899[k];

        t_1194[k] = f_12 * msh_900[k]
                    + f_3 * pc_x[k] * nsh_900[k];

        t_1195[k] = f_12 * msh_901[k]
                    + f_3 * pc_x[k] * nsh_901[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, pc_z, msh_708, msh_729, msh_902, \
                         nsg0_640, nsg1_640, nsh_897, nsh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_12 * msh_902[k]
                    + f_3 * pc_x[k] * nsh_902[k];

        t_1197[k] = f_12 * msh_729[k]
                    + f_1 * nsg0_640[k]
                    - f_2 * nsg1_640[k]
                    + f_3 * pc_y[k] * nsh_897[k];

        t_1198[k] = f_20 * msh_708[k]
                    + f_3 * pc_z[k] * nsh_897[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_y, msh_731, msh_732, msh_733, nsg0_642, \
                         nsg0_643, nsg0_644, nsg1_642, nsg1_643, nsg1_644, nsh_899, nsh_900, \
                         nsh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_12 * msh_731[k]
                    + f_8 * nsg0_642[k]
                    - f_9 * nsg1_642[k]
                    + f_3 * pc_y[k] * nsh_899[k];

        t_1200[k] = f_12 * msh_732[k]
                    + f_6 * nsg0_643[k]
                    - f_7 * nsg1_643[k]
                    + f_3 * pc_y[k] * nsh_900[k];

        t_1201[k] = f_12 * msh_733[k]
                    + f_4 * nsg0_644[k]
                    - f_5 * nsg1_644[k]
                    + f_3 * pc_y[k] * nsh_901[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, pa_y, pc_y, pc_z, msi0_980, msh_713, \
                         msh_734, msh_735, msi1_980, nsg0_644, nsg1_644, nsh_902, \
                         nsh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_12 * msh_734[k]
                    + f_3 * pc_y[k] * nsh_902[k];

        t_1203[k] = f_20 * msh_713[k]
                    + f_1 * nsg0_644[k]
                    - f_2 * nsg1_644[k]
                    + f_3 * pc_z[k] * nsh_902[k];

        t_1204[k] = pa_y[k] * msi0_980[k]
                    - f_10 * pc_y[k] * msi1_980[k];

        t_1205[k] = f_11 * msh_735[k]
                    + f_3 * pc_y[k] * nsh_903[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pa_y, pc_y, pc_z, msi0_983, msi0_985, \
                         msh_714, msh_736, msh_737, msi1_983, msi1_985, nsh_903, \
                         nsh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_19 * msh_714[k]
                    + f_3 * pc_z[k] * nsh_903[k];

        t_1207[k] = pa_y[k] * msi0_983[k]
                    + f_12 * msh_736[k]
                    - f_10 * pc_y[k] * msi1_983[k];

        t_1208[k] = f_11 * msh_737[k]
                    + f_3 * pc_y[k] * nsh_905[k];

        t_1209[k] = pa_y[k] * msi0_985[k]
                    - f_10 * pc_y[k] * msi1_985[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pa_y, pc_y, pc_z, msi0_986, msi0_989, \
                         msh_717, msh_738, msh_740, msi1_986, msi1_989, nsh_906, \
                         nsh_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = pa_y[k] * msi0_986[k]
                    + f_13 * msh_738[k]
                    - f_10 * pc_y[k] * msi1_986[k];

        t_1211[k] = f_19 * msh_717[k]
                    + f_3 * pc_z[k] * nsh_906[k];

        t_1212[k] = f_11 * msh_740[k]
                    + f_3 * pc_y[k] * nsh_908[k];

        t_1213[k] = pa_y[k] * msi0_989[k]
                    - f_10 * pc_y[k] * msi1_989[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, pa_y, pc_y, pc_z, msi0_990, msi0_992, \
                         msh_720, msh_741, msh_743, msi1_990, msi1_992, \
                         nsh_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pa_y[k] * msi0_990[k]
                    + f_14 * msh_741[k]
                    - f_10 * pc_y[k] * msi1_990[k];

        t_1215[k] = f_19 * msh_720[k]
                    + f_3 * pc_z[k] * nsh_909[k];

        t_1216[k] = pa_y[k] * msi0_992[k]
                    + f_12 * msh_743[k]
                    - f_10 * pc_y[k] * msi1_992[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, pa_y, pc_x, pc_y, msi0_994, msh_744, \
                         msh_918, msh_919, msi1_994, nsh_912, nsh_918, \
                         nsh_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_11 * msh_744[k]
                    + f_3 * pc_y[k] * nsh_912[k];

        t_1218[k] = pa_y[k] * msi0_994[k]
                    - f_10 * pc_y[k] * msi1_994[k];

        t_1219[k] = f_12 * msh_918[k]
                    + f_3 * pc_x[k] * nsh_918[k];

        t_1220[k] = f_12 * msh_919[k]
                    + f_3 * pc_x[k] * nsh_919[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, t_1224, pc_x, msh_920, msh_921, msh_922, \
                         msh_923, nsh_920, nsh_921, nsh_922, nsh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_12 * msh_920[k]
                    + f_3 * pc_x[k] * nsh_920[k];

        t_1222[k] = f_12 * msh_921[k]
                    + f_3 * pc_x[k] * nsh_921[k];

        t_1223[k] = f_12 * msh_922[k]
                    + f_3 * pc_x[k] * nsh_922[k];

        t_1224[k] = f_12 * msh_923[k]
                    + f_3 * pc_x[k] * nsh_923[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, pc_y, pc_z, msh_729, msh_750, msh_752, \
                         nsg0_655, nsg0_657, nsg1_655, nsg1_657, nsh_918, \
                         nsh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_11 * msh_750[k]
                    + f_1 * nsg0_655[k]
                    - f_2 * nsg1_655[k]
                    + f_3 * pc_y[k] * nsh_918[k];

        t_1226[k] = f_19 * msh_729[k]
                    + f_3 * pc_z[k] * nsh_918[k];

        t_1227[k] = f_11 * msh_752[k]
                    + f_8 * nsg0_657[k]
                    - f_9 * nsg1_657[k]
                    + f_3 * pc_y[k] * nsh_920[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, pc_y, msh_753, msh_754, msh_755, nsg0_658, \
                         nsg0_659, nsg1_658, nsg1_659, nsh_921, nsh_922, \
                         nsh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = f_11 * msh_753[k]
                    + f_6 * nsg0_658[k]
                    - f_7 * nsg1_658[k]
                    + f_3 * pc_y[k] * nsh_921[k];

        t_1229[k] = f_11 * msh_754[k]
                    + f_4 * nsg0_659[k]
                    - f_5 * nsg1_659[k]
                    + f_3 * pc_y[k] * nsh_922[k];

        t_1230[k] = f_11 * msh_755[k]
                    + f_3 * pc_y[k] * nsh_923[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pa_y, pc_x, pc_y, pc_z, msi0_1007, \
                         msh_735, msh_924, msi1_1007, nsg0_660, nsg1_660, \
                         nsh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = pa_y[k] * msi0_1007[k]
                    - f_10 * pc_y[k] * msi1_1007[k];

        t_1232[k] = f_12 * msh_924[k]
                    + f_1 * nsg0_660[k]
                    - f_2 * nsg1_660[k]
                    + f_3 * pc_x[k] * nsh_924[k];

        t_1233[k] = f_3 * pc_y[k] * nsh_924[k];

        t_1234[k] = f_18 * msh_735[k]
                    + f_3 * pc_z[k] * nsh_924[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pc_x, pc_y, msh_929, nsg0_660, nsg0_665, \
                         nsg1_660, nsg1_665, nsh_925, nsh_926, \
                         nsh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_4 * nsg0_660[k]
                    - f_5 * nsg1_660[k]
                    + f_3 * pc_y[k] * nsh_925[k];

        t_1236[k] = f_3 * pc_y[k] * nsh_926[k];

        t_1237[k] = f_12 * msh_929[k]
                    + f_8 * nsg0_665[k]
                    - f_9 * nsg1_665[k]
                    + f_3 * pc_x[k] * nsh_929[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, pc_y, nsg0_661, nsg0_662, nsg1_661, nsg1_662, \
                         nsh_927, nsh_928, nsh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_6 * nsg0_661[k]
                    - f_7 * nsg1_661[k]
                    + f_3 * pc_y[k] * nsh_927[k];

        t_1239[k] = f_4 * nsg0_662[k]
                    - f_5 * nsg1_662[k]
                    + f_3 * pc_y[k] * nsh_928[k];

        t_1240[k] = f_3 * pc_y[k] * nsh_929[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pc_y, msh_933, nsg0_663, nsg0_664, \
                         nsg0_669, nsg1_663, nsg1_664, nsg1_669, nsh_930, nsh_931, \
                         nsh_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_12 * msh_933[k]
                    + f_6 * nsg0_669[k]
                    - f_7 * nsg1_669[k]
                    + f_3 * pc_x[k] * nsh_933[k];

        t_1242[k] = f_8 * nsg0_663[k]
                    - f_9 * nsg1_663[k]
                    + f_3 * pc_y[k] * nsh_930[k];

        t_1243[k] = f_6 * nsg0_664[k]
                    - f_7 * nsg1_664[k]
                    + f_3 * pc_y[k] * nsh_931[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, pc_x, pc_y, msh_938, msh_939, \
                         nsg0_665, nsg0_674, nsg1_665, nsg1_674, nsh_932, nsh_933, nsh_938, \
                         nsh_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_4 * nsg0_665[k]
                    - f_5 * nsg1_665[k]
                    + f_3 * pc_y[k] * nsh_932[k];

        t_1245[k] = f_3 * pc_y[k] * nsh_933[k];

        t_1246[k] = f_12 * msh_938[k]
                    + f_4 * nsg0_674[k]
                    - f_5 * nsg1_674[k]
                    + f_3 * pc_x[k] * nsh_938[k];

        t_1247[k] = f_12 * msh_939[k]
                    + f_3 * pc_x[k] * nsh_939[k];
    }

#pragma omp simd aligned(t_1248, t_1249, t_1250, t_1251, t_1252, pc_x, pc_y, msh_940, msh_941, \
                         msh_942, msh_944, nsh_938, nsh_940, nsh_941, nsh_942, \
                         nsh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1248[k] = f_12 * msh_940[k]
                    + f_3 * pc_x[k] * nsh_940[k];

        t_1249[k] = f_12 * msh_941[k]
                    + f_3 * pc_x[k] * nsh_941[k];

        t_1250[k] = f_12 * msh_942[k]
                    + f_3 * pc_x[k] * nsh_942[k];

        t_1251[k] = f_3 * pc_y[k] * nsh_938[k];

        t_1252[k] = f_12 * msh_944[k]
                    + f_3 * pc_x[k] * nsh_944[k];
    }

#pragma omp simd aligned(t_1253, t_1254, t_1255, pc_y, nsg0_670, nsg0_671, nsg0_672, nsg1_670, \
                         nsg1_671, nsg1_672, nsh_939, nsh_940, \
                         nsh_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1253[k] = f_1 * nsg0_670[k]
                    - f_2 * nsg1_670[k]
                    + f_3 * pc_y[k] * nsh_939[k];

        t_1254[k] = f_16 * nsg0_671[k]
                    - f_17 * nsg1_671[k]
                    + f_3 * pc_y[k] * nsh_940[k];

        t_1255[k] = f_8 * nsg0_672[k]
                    - f_9 * nsg1_672[k]
                    + f_3 * pc_y[k] * nsh_941[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pc_y, pc_z, msh_755, nsg0_673, \
                         nsg0_674, nsg1_673, nsg1_674, nsh_942, nsh_943, \
                         nsh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_6 * nsg0_673[k]
                    - f_7 * nsg1_673[k]
                    + f_3 * pc_y[k] * nsh_942[k];

        t_1257[k] = f_4 * nsg0_674[k]
                    - f_5 * nsg1_674[k]
                    + f_3 * pc_y[k] * nsh_943[k];

        t_1258[k] = f_3 * pc_y[k] * nsh_944[k];

        t_1259[k] = f_18 * msh_755[k]
                    + f_1 * nsg0_674[k]
                    - f_2 * nsg1_674[k]
                    + f_3 * pc_z[k] * nsh_944[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msi0,
                                                           const size_t msh, const size_t msi1,
                                                           const size_t nsg0, const size_t nsg1,
                                                           const size_t nsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.5 / q;
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_1008 = buffer.data(msi0 + 1008);
    const auto *msi0_1011 = buffer.data(msi0 + 1011);
    const auto *msi0_1014 = buffer.data(msi0 + 1014);
    const auto *msi0_1018 = buffer.data(msi0 + 1018);
    const auto *msi0_1260 = buffer.data(msi0 + 1260);
    const auto *msi0_1263 = buffer.data(msi0 + 1263);
    const auto *msi0_1266 = buffer.data(msi0 + 1266);
    const auto *msi0_1270 = buffer.data(msi0 + 1270);
    const auto *msi0_1281 = buffer.data(msi0 + 1281);
    const auto *msi0_1283 = buffer.data(msi0 + 1283);
    const auto *msi0_1284 = buffer.data(msi0 + 1284);
    const auto *msi0_1285 = buffer.data(msi0 + 1285);
    const auto *msi0_1287 = buffer.data(msi0 + 1287);
    const auto *msi0_1293 = buffer.data(msi0 + 1293);
    const auto *msi0_1297 = buffer.data(msi0 + 1297);
    const auto *msi0_1300 = buffer.data(msi0 + 1300);
    const auto *msi0_1302 = buffer.data(msi0 + 1302);
    const auto *msi0_1309 = buffer.data(msi0 + 1309);
    const auto *msi0_1311 = buffer.data(msi0 + 1311);
    const auto *msi0_1312 = buffer.data(msi0 + 1312);
    const auto *msi0_1313 = buffer.data(msi0 + 1313);
    const auto *msi0_1315 = buffer.data(msi0 + 1315);
    const auto *msi0_1316 = buffer.data(msi0 + 1316);
    const auto *msi0_1319 = buffer.data(msi0 + 1319);
    const auto *msi0_1321 = buffer.data(msi0 + 1321);
    const auto *msi0_1322 = buffer.data(msi0 + 1322);
    const auto *msi0_1325 = buffer.data(msi0 + 1325);
    const auto *msi0_1326 = buffer.data(msi0 + 1326);
    const auto *msi0_1328 = buffer.data(msi0 + 1328);
    const auto *msi0_1330 = buffer.data(msi0 + 1330);
    const auto *msi0_1337 = buffer.data(msi0 + 1337);
    const auto *msi0_1339 = buffer.data(msi0 + 1339);
    const auto *msi0_1340 = buffer.data(msi0 + 1340);
    const auto *msi0_1341 = buffer.data(msi0 + 1341);
    const auto *msi0_1343 = buffer.data(msi0 + 1343);
    const auto *msi0_1344 = buffer.data(msi0 + 1344);
    const auto *msi0_1347 = buffer.data(msi0 + 1347);
    const auto *msi0_1349 = buffer.data(msi0 + 1349);
    const auto *msi0_1350 = buffer.data(msi0 + 1350);
    const auto *msi0_1353 = buffer.data(msi0 + 1353);
    const auto *msi0_1354 = buffer.data(msi0 + 1354);
    const auto *msi0_1356 = buffer.data(msi0 + 1356);
    const auto *msi0_1358 = buffer.data(msi0 + 1358);
    const auto *msi0_1365 = buffer.data(msi0 + 1365);
    const auto *msi0_1367 = buffer.data(msi0 + 1367);
    const auto *msi0_1368 = buffer.data(msi0 + 1368);
    const auto *msi0_1369 = buffer.data(msi0 + 1369);
    const auto *msi0_1371 = buffer.data(msi0 + 1371);
    const auto *msi0_1372 = buffer.data(msi0 + 1372);
    const auto *msi0_1375 = buffer.data(msi0 + 1375);
    const auto *msi0_1377 = buffer.data(msi0 + 1377);
    const auto *msi0_1378 = buffer.data(msi0 + 1378);
    const auto *msi0_1381 = buffer.data(msi0 + 1381);
    const auto *msi0_1382 = buffer.data(msi0 + 1382);

    const auto *msh_756 = buffer.data(msh + 756);
    const auto *msh_759 = buffer.data(msh + 759);
    const auto *msh_761 = buffer.data(msh + 761);
    const auto *msh_762 = buffer.data(msh + 762);
    const auto *msh_765 = buffer.data(msh + 765);
    const auto *msh_771 = buffer.data(msh + 771);
    const auto *msh_776 = buffer.data(msh + 776);
    const auto *msh_777 = buffer.data(msh + 777);
    const auto *msh_779 = buffer.data(msh + 779);
    const auto *msh_780 = buffer.data(msh + 780);
    const auto *msh_782 = buffer.data(msh + 782);
    const auto *msh_783 = buffer.data(msh + 783);
    const auto *msh_786 = buffer.data(msh + 786);
    const auto *msh_792 = buffer.data(msh + 792);
    const auto *msh_797 = buffer.data(msh + 797);
    const auto *msh_798 = buffer.data(msh + 798);
    const auto *msh_800 = buffer.data(msh + 800);
    const auto *msh_801 = buffer.data(msh + 801);
    const auto *msh_803 = buffer.data(msh + 803);
    const auto *msh_804 = buffer.data(msh + 804);
    const auto *msh_807 = buffer.data(msh + 807);
    const auto *msh_813 = buffer.data(msh + 813);
    const auto *msh_818 = buffer.data(msh + 818);
    const auto *msh_819 = buffer.data(msh + 819);
    const auto *msh_821 = buffer.data(msh + 821);
    const auto *msh_822 = buffer.data(msh + 822);
    const auto *msh_824 = buffer.data(msh + 824);
    const auto *msh_828 = buffer.data(msh + 828);
    const auto *msh_839 = buffer.data(msh + 839);
    const auto *msh_840 = buffer.data(msh + 840);
    const auto *msh_842 = buffer.data(msh + 842);
    const auto *msh_845 = buffer.data(msh + 845);
    const auto *msh_945 = buffer.data(msh + 945);
    const auto *msh_948 = buffer.data(msh + 948);
    const auto *msh_951 = buffer.data(msh + 951);
    const auto *msh_955 = buffer.data(msh + 955);
    const auto *msh_960 = buffer.data(msh + 960);
    const auto *msh_962 = buffer.data(msh + 962);
    const auto *msh_963 = buffer.data(msh + 963);
    const auto *msh_964 = buffer.data(msh + 964);
    const auto *msh_965 = buffer.data(msh + 965);
    const auto *msh_971 = buffer.data(msh + 971);
    const auto *msh_975 = buffer.data(msh + 975);
    const auto *msh_978 = buffer.data(msh + 978);
    const auto *msh_980 = buffer.data(msh + 980);
    const auto *msh_981 = buffer.data(msh + 981);
    const auto *msh_982 = buffer.data(msh + 982);
    const auto *msh_983 = buffer.data(msh + 983);
    const auto *msh_984 = buffer.data(msh + 984);
    const auto *msh_985 = buffer.data(msh + 985);
    const auto *msh_986 = buffer.data(msh + 986);
    const auto *msh_987 = buffer.data(msh + 987);
    const auto *msh_990 = buffer.data(msh + 990);
    const auto *msh_992 = buffer.data(msh + 992);
    const auto *msh_993 = buffer.data(msh + 993);
    const auto *msh_996 = buffer.data(msh + 996);
    const auto *msh_997 = buffer.data(msh + 997);
    const auto *msh_999 = buffer.data(msh + 999);
    const auto *msh_1001 = buffer.data(msh + 1001);
    const auto *msh_1002 = buffer.data(msh + 1002);
    const auto *msh_1003 = buffer.data(msh + 1003);
    const auto *msh_1004 = buffer.data(msh + 1004);
    const auto *msh_1005 = buffer.data(msh + 1005);
    const auto *msh_1006 = buffer.data(msh + 1006);
    const auto *msh_1007 = buffer.data(msh + 1007);
    const auto *msh_1008 = buffer.data(msh + 1008);
    const auto *msh_1011 = buffer.data(msh + 1011);
    const auto *msh_1013 = buffer.data(msh + 1013);
    const auto *msh_1014 = buffer.data(msh + 1014);
    const auto *msh_1017 = buffer.data(msh + 1017);
    const auto *msh_1018 = buffer.data(msh + 1018);
    const auto *msh_1020 = buffer.data(msh + 1020);
    const auto *msh_1022 = buffer.data(msh + 1022);
    const auto *msh_1023 = buffer.data(msh + 1023);
    const auto *msh_1024 = buffer.data(msh + 1024);
    const auto *msh_1025 = buffer.data(msh + 1025);
    const auto *msh_1026 = buffer.data(msh + 1026);
    const auto *msh_1027 = buffer.data(msh + 1027);
    const auto *msh_1028 = buffer.data(msh + 1028);
    const auto *msh_1029 = buffer.data(msh + 1029);
    const auto *msh_1032 = buffer.data(msh + 1032);
    const auto *msh_1034 = buffer.data(msh + 1034);
    const auto *msh_1035 = buffer.data(msh + 1035);
    const auto *msh_1038 = buffer.data(msh + 1038);
    const auto *msh_1039 = buffer.data(msh + 1039);

    const auto *msi1_1008 = buffer.data(msi1 + 1008);
    const auto *msi1_1011 = buffer.data(msi1 + 1011);
    const auto *msi1_1014 = buffer.data(msi1 + 1014);
    const auto *msi1_1018 = buffer.data(msi1 + 1018);
    const auto *msi1_1260 = buffer.data(msi1 + 1260);
    const auto *msi1_1263 = buffer.data(msi1 + 1263);
    const auto *msi1_1266 = buffer.data(msi1 + 1266);
    const auto *msi1_1270 = buffer.data(msi1 + 1270);
    const auto *msi1_1281 = buffer.data(msi1 + 1281);
    const auto *msi1_1283 = buffer.data(msi1 + 1283);
    const auto *msi1_1284 = buffer.data(msi1 + 1284);
    const auto *msi1_1285 = buffer.data(msi1 + 1285);
    const auto *msi1_1287 = buffer.data(msi1 + 1287);
    const auto *msi1_1293 = buffer.data(msi1 + 1293);
    const auto *msi1_1297 = buffer.data(msi1 + 1297);
    const auto *msi1_1300 = buffer.data(msi1 + 1300);
    const auto *msi1_1302 = buffer.data(msi1 + 1302);
    const auto *msi1_1309 = buffer.data(msi1 + 1309);
    const auto *msi1_1311 = buffer.data(msi1 + 1311);
    const auto *msi1_1312 = buffer.data(msi1 + 1312);
    const auto *msi1_1313 = buffer.data(msi1 + 1313);
    const auto *msi1_1315 = buffer.data(msi1 + 1315);
    const auto *msi1_1316 = buffer.data(msi1 + 1316);
    const auto *msi1_1319 = buffer.data(msi1 + 1319);
    const auto *msi1_1321 = buffer.data(msi1 + 1321);
    const auto *msi1_1322 = buffer.data(msi1 + 1322);
    const auto *msi1_1325 = buffer.data(msi1 + 1325);
    const auto *msi1_1326 = buffer.data(msi1 + 1326);
    const auto *msi1_1328 = buffer.data(msi1 + 1328);
    const auto *msi1_1330 = buffer.data(msi1 + 1330);
    const auto *msi1_1337 = buffer.data(msi1 + 1337);
    const auto *msi1_1339 = buffer.data(msi1 + 1339);
    const auto *msi1_1340 = buffer.data(msi1 + 1340);
    const auto *msi1_1341 = buffer.data(msi1 + 1341);
    const auto *msi1_1343 = buffer.data(msi1 + 1343);
    const auto *msi1_1344 = buffer.data(msi1 + 1344);
    const auto *msi1_1347 = buffer.data(msi1 + 1347);
    const auto *msi1_1349 = buffer.data(msi1 + 1349);
    const auto *msi1_1350 = buffer.data(msi1 + 1350);
    const auto *msi1_1353 = buffer.data(msi1 + 1353);
    const auto *msi1_1354 = buffer.data(msi1 + 1354);
    const auto *msi1_1356 = buffer.data(msi1 + 1356);
    const auto *msi1_1358 = buffer.data(msi1 + 1358);
    const auto *msi1_1365 = buffer.data(msi1 + 1365);
    const auto *msi1_1367 = buffer.data(msi1 + 1367);
    const auto *msi1_1368 = buffer.data(msi1 + 1368);
    const auto *msi1_1369 = buffer.data(msi1 + 1369);
    const auto *msi1_1371 = buffer.data(msi1 + 1371);
    const auto *msi1_1372 = buffer.data(msi1 + 1372);
    const auto *msi1_1375 = buffer.data(msi1 + 1375);
    const auto *msi1_1377 = buffer.data(msi1 + 1377);
    const auto *msi1_1378 = buffer.data(msi1 + 1378);
    const auto *msi1_1381 = buffer.data(msi1 + 1381);
    const auto *msi1_1382 = buffer.data(msi1 + 1382);

    const auto *nsg0_675 = buffer.data(nsg0 + 675);
    const auto *nsg0_677 = buffer.data(nsg0 + 677);
    const auto *nsg0_678 = buffer.data(nsg0 + 678);
    const auto *nsg0_680 = buffer.data(nsg0 + 680);

    const auto *nsg1_675 = buffer.data(nsg1 + 675);
    const auto *nsg1_677 = buffer.data(nsg1 + 677);
    const auto *nsg1_678 = buffer.data(nsg1 + 678);
    const auto *nsg1_680 = buffer.data(nsg1 + 680);

    const auto *nsh_945 = buffer.data(nsh + 945);
    const auto *nsh_946 = buffer.data(nsh + 946);
    const auto *nsh_947 = buffer.data(nsh + 947);
    const auto *nsh_948 = buffer.data(nsh + 948);
    const auto *nsh_950 = buffer.data(nsh + 950);
    const auto *nsh_951 = buffer.data(nsh + 951);
    const auto *nsh_952 = buffer.data(nsh + 952);
    const auto *nsh_954 = buffer.data(nsh + 954);
    const auto *nsh_955 = buffer.data(nsh + 955);
    const auto *nsh_960 = buffer.data(nsh + 960);
    const auto *nsh_962 = buffer.data(nsh + 962);
    const auto *nsh_963 = buffer.data(nsh + 963);
    const auto *nsh_964 = buffer.data(nsh + 964);
    const auto *nsh_965 = buffer.data(nsh + 965);
    const auto *nsh_966 = buffer.data(nsh + 966);
    const auto *nsh_968 = buffer.data(nsh + 968);
    const auto *nsh_969 = buffer.data(nsh + 969);
    const auto *nsh_971 = buffer.data(nsh + 971);
    const auto *nsh_972 = buffer.data(nsh + 972);
    const auto *nsh_975 = buffer.data(nsh + 975);
    const auto *nsh_981 = buffer.data(nsh + 981);
    const auto *nsh_982 = buffer.data(nsh + 982);
    const auto *nsh_983 = buffer.data(nsh + 983);
    const auto *nsh_984 = buffer.data(nsh + 984);
    const auto *nsh_985 = buffer.data(nsh + 985);
    const auto *nsh_986 = buffer.data(nsh + 986);
    const auto *nsh_987 = buffer.data(nsh + 987);
    const auto *nsh_989 = buffer.data(nsh + 989);
    const auto *nsh_990 = buffer.data(nsh + 990);
    const auto *nsh_992 = buffer.data(nsh + 992);
    const auto *nsh_993 = buffer.data(nsh + 993);
    const auto *nsh_996 = buffer.data(nsh + 996);
    const auto *nsh_1002 = buffer.data(nsh + 1002);
    const auto *nsh_1003 = buffer.data(nsh + 1003);
    const auto *nsh_1004 = buffer.data(nsh + 1004);
    const auto *nsh_1005 = buffer.data(nsh + 1005);
    const auto *nsh_1006 = buffer.data(nsh + 1006);
    const auto *nsh_1007 = buffer.data(nsh + 1007);
    const auto *nsh_1008 = buffer.data(nsh + 1008);
    const auto *nsh_1010 = buffer.data(nsh + 1010);
    const auto *nsh_1011 = buffer.data(nsh + 1011);
    const auto *nsh_1013 = buffer.data(nsh + 1013);
    const auto *nsh_1014 = buffer.data(nsh + 1014);
    const auto *nsh_1017 = buffer.data(nsh + 1017);
    const auto *nsh_1023 = buffer.data(nsh + 1023);
    const auto *nsh_1024 = buffer.data(nsh + 1024);
    const auto *nsh_1025 = buffer.data(nsh + 1025);
    const auto *nsh_1026 = buffer.data(nsh + 1026);
    const auto *nsh_1027 = buffer.data(nsh + 1027);
    const auto *nsh_1028 = buffer.data(nsh + 1028);
    const auto *nsh_1029 = buffer.data(nsh + 1029);
    const auto *nsh_1031 = buffer.data(nsh + 1031);
    const auto *nsh_1032 = buffer.data(nsh + 1032);
    const auto *nsh_1034 = buffer.data(nsh + 1034);

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pa_x, pc_x, pc_y, pc_z, msi0_1260, \
                         msi0_1263, msh_756, msh_945, msh_948, msi1_1260, msi1_1263, \
                         nsh_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = pa_x[k] * msi0_1260[k]
                    + f_20 * msh_945[k]
                    - f_10 * pc_x[k] * msi1_1260[k];

        t_1261[k] = f_15 * msh_756[k]
                    + f_3 * pc_y[k] * nsh_945[k];

        t_1262[k] = f_3 * pc_z[k] * nsh_945[k];

        t_1263[k] = pa_x[k] * msi0_1263[k]
                    + f_14 * msh_948[k]
                    - f_10 * pc_x[k] * msi1_1263[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, t_1267, pa_x, pc_x, pc_z, msi0_1266, msh_951, \
                         msi1_1266, nsg0_675, nsg1_675, nsh_946, nsh_947, \
                         nsh_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_3 * pc_z[k] * nsh_946[k];

        t_1265[k] = f_4 * nsg0_675[k]
                    - f_5 * nsg1_675[k]
                    + f_3 * pc_z[k] * nsh_947[k];

        t_1266[k] = pa_x[k] * msi0_1266[k]
                    + f_13 * msh_951[k]
                    - f_10 * pc_x[k] * msi1_1266[k];

        t_1267[k] = f_3 * pc_z[k] * nsh_948[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, t_1271, pa_x, pc_x, pc_y, pc_z, msi0_1270, \
                         msh_761, msh_955, msi1_1270, nsg0_677, nsg1_677, nsh_950, \
                         nsh_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_15 * msh_761[k]
                    + f_3 * pc_y[k] * nsh_950[k];

        t_1269[k] = f_6 * nsg0_677[k]
                    - f_7 * nsg1_677[k]
                    + f_3 * pc_z[k] * nsh_950[k];

        t_1270[k] = pa_x[k] * msi0_1270[k]
                    + f_12 * msh_955[k]
                    - f_10 * pc_x[k] * msi1_1270[k];

        t_1271[k] = f_3 * pc_z[k] * nsh_951[k];
    }

#pragma omp simd aligned(t_1272, t_1273, t_1274, t_1275, pc_x, pc_y, pc_z, msh_765, msh_960, \
                         nsg0_678, nsg0_680, nsg1_678, nsg1_680, nsh_952, nsh_954, \
                         nsh_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1272[k] = f_4 * nsg0_678[k]
                    - f_5 * nsg1_678[k]
                    + f_3 * pc_z[k] * nsh_952[k];

        t_1273[k] = f_15 * msh_765[k]
                    + f_3 * pc_y[k] * nsh_954[k];

        t_1274[k] = f_8 * nsg0_680[k]
                    - f_9 * nsg1_680[k]
                    + f_3 * pc_z[k] * nsh_954[k];

        t_1275[k] = f_11 * msh_960[k]
                    + f_3 * pc_x[k] * nsh_960[k];
    }

#pragma omp simd aligned(t_1276, t_1277, t_1278, t_1279, t_1280, pc_x, pc_z, msh_962, msh_963, \
                         msh_964, msh_965, nsh_955, nsh_962, nsh_963, nsh_964, \
                         nsh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1276[k] = f_3 * pc_z[k] * nsh_955[k];

        t_1277[k] = f_11 * msh_962[k]
                    + f_3 * pc_x[k] * nsh_962[k];

        t_1278[k] = f_11 * msh_963[k]
                    + f_3 * pc_x[k] * nsh_963[k];

        t_1279[k] = f_11 * msh_964[k]
                    + f_3 * pc_x[k] * nsh_964[k];

        t_1280[k] = f_11 * msh_965[k]
                    + f_3 * pc_x[k] * nsh_965[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, t_1284, pa_x, pc_x, pc_z, msi0_1281, \
                         msi0_1283, msi0_1284, msi1_1281, msi1_1283, msi1_1284, \
                         nsh_960 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = pa_x[k] * msi0_1281[k]
                    - f_10 * pc_x[k] * msi1_1281[k];

        t_1282[k] = f_3 * pc_z[k] * nsh_960[k];

        t_1283[k] = pa_x[k] * msi0_1283[k]
                    - f_10 * pc_x[k] * msi1_1283[k];

        t_1284[k] = pa_x[k] * msi0_1284[k]
                    - f_10 * pc_x[k] * msi1_1284[k];
    }

#pragma omp simd aligned(t_1285, t_1286, t_1287, pa_x, pc_x, pc_y, msi0_1285, msi0_1287, \
                         msh_776, msi1_1285, msi1_1287, nsh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1285[k] = pa_x[k] * msi0_1285[k]
                    - f_10 * pc_x[k] * msi1_1285[k];

        t_1286[k] = f_15 * msh_776[k]
                    + f_3 * pc_y[k] * nsh_965[k];

        t_1287[k] = pa_x[k] * msi0_1287[k]
                    - f_10 * pc_x[k] * msi1_1287[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pa_z, pc_y, pc_z, msi0_1008, \
                         msi0_1011, msh_756, msh_777, msi1_1008, msi1_1011, \
                         nsh_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = pa_z[k] * msi0_1008[k]
                    - f_10 * pc_z[k] * msi1_1008[k];

        t_1289[k] = f_18 * msh_777[k]
                    + f_3 * pc_y[k] * nsh_966[k];

        t_1290[k] = f_11 * msh_756[k]
                    + f_3 * pc_z[k] * nsh_966[k];

        t_1291[k] = pa_z[k] * msi0_1011[k]
                    - f_10 * pc_z[k] * msi1_1011[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, pa_x, pa_z, pc_x, pc_y, pc_z, msi0_1014, \
                         msi0_1293, msh_779, msh_971, msi1_1014, msi1_1293, \
                         nsh_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_18 * msh_779[k]
                    + f_3 * pc_y[k] * nsh_968[k];

        t_1293[k] = pa_x[k] * msi0_1293[k]
                    + f_14 * msh_971[k]
                    - f_10 * pc_x[k] * msi1_1293[k];

        t_1294[k] = pa_z[k] * msi0_1014[k]
                    - f_10 * pc_z[k] * msi1_1014[k];
    }

#pragma omp simd aligned(t_1295, t_1296, t_1297, pa_x, pc_x, pc_y, pc_z, msi0_1297, msh_759, \
                         msh_782, msh_975, msi1_1297, nsh_969, \
                         nsh_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1295[k] = f_11 * msh_759[k]
                    + f_3 * pc_z[k] * nsh_969[k];

        t_1296[k] = f_18 * msh_782[k]
                    + f_3 * pc_y[k] * nsh_971[k];

        t_1297[k] = pa_x[k] * msi0_1297[k]
                    + f_13 * msh_975[k]
                    - f_10 * pc_x[k] * msi1_1297[k];
    }

#pragma omp simd aligned(t_1298, t_1299, t_1300, pa_x, pa_z, pc_x, pc_z, msi0_1018, msi0_1300, \
                         msh_762, msh_978, msi1_1018, msi1_1300, \
                         nsh_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1298[k] = pa_z[k] * msi0_1018[k]
                    - f_10 * pc_z[k] * msi1_1018[k];

        t_1299[k] = f_11 * msh_762[k]
                    + f_3 * pc_z[k] * nsh_972[k];

        t_1300[k] = pa_x[k] * msi0_1300[k]
                    + f_12 * msh_978[k]
                    - f_10 * pc_x[k] * msi1_1300[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, pa_x, pc_x, pc_y, msi0_1302, msh_786, \
                         msh_980, msh_981, msh_982, msi1_1302, nsh_975, nsh_981, \
                         nsh_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_18 * msh_786[k]
                    + f_3 * pc_y[k] * nsh_975[k];

        t_1302[k] = pa_x[k] * msi0_1302[k]
                    + f_12 * msh_980[k]
                    - f_10 * pc_x[k] * msi1_1302[k];

        t_1303[k] = f_11 * msh_981[k]
                    + f_3 * pc_x[k] * nsh_981[k];

        t_1304[k] = f_11 * msh_982[k]
                    + f_3 * pc_x[k] * nsh_982[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, t_1308, pc_x, msh_983, msh_984, msh_985, \
                         msh_986, nsh_983, nsh_984, nsh_985, nsh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = f_11 * msh_983[k]
                    + f_3 * pc_x[k] * nsh_983[k];

        t_1306[k] = f_11 * msh_984[k]
                    + f_3 * pc_x[k] * nsh_984[k];

        t_1307[k] = f_11 * msh_985[k]
                    + f_3 * pc_x[k] * nsh_985[k];

        t_1308[k] = f_11 * msh_986[k]
                    + f_3 * pc_x[k] * nsh_986[k];
    }

#pragma omp simd aligned(t_1309, t_1310, t_1311, t_1312, pa_x, pc_x, pc_z, msi0_1309, \
                         msi0_1311, msi0_1312, msh_771, msi1_1309, msi1_1311, msi1_1312, \
                         nsh_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1309[k] = pa_x[k] * msi0_1309[k]
                    - f_10 * pc_x[k] * msi1_1309[k];

        t_1310[k] = f_11 * msh_771[k]
                    + f_3 * pc_z[k] * nsh_981[k];

        t_1311[k] = pa_x[k] * msi0_1311[k]
                    - f_10 * pc_x[k] * msi1_1311[k];

        t_1312[k] = pa_x[k] * msi0_1312[k]
                    - f_10 * pc_x[k] * msi1_1312[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, t_1316, pa_x, pc_x, pc_y, msi0_1313, \
                         msi0_1315, msi0_1316, msh_797, msh_987, msi1_1313, msi1_1315, \
                         msi1_1316, nsh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = pa_x[k] * msi0_1313[k]
                    - f_10 * pc_x[k] * msi1_1313[k];

        t_1314[k] = f_18 * msh_797[k]
                    + f_3 * pc_y[k] * nsh_986[k];

        t_1315[k] = pa_x[k] * msi0_1315[k]
                    - f_10 * pc_x[k] * msi1_1315[k];

        t_1316[k] = pa_x[k] * msi0_1316[k]
                    + f_20 * msh_987[k]
                    - f_10 * pc_x[k] * msi1_1316[k];
    }

#pragma omp simd aligned(t_1317, t_1318, t_1319, t_1320, pa_x, pc_x, pc_y, pc_z, msi0_1319, \
                         msh_777, msh_798, msh_800, msh_990, msi1_1319, nsh_987, \
                         nsh_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1317[k] = f_19 * msh_798[k]
                    + f_3 * pc_y[k] * nsh_987[k];

        t_1318[k] = f_12 * msh_777[k]
                    + f_3 * pc_z[k] * nsh_987[k];

        t_1319[k] = pa_x[k] * msi0_1319[k]
                    + f_14 * msh_990[k]
                    - f_10 * pc_x[k] * msi1_1319[k];

        t_1320[k] = f_19 * msh_800[k]
                    + f_3 * pc_y[k] * nsh_989[k];
    }

#pragma omp simd aligned(t_1321, t_1322, t_1323, pa_x, pc_x, pc_z, msi0_1321, msi0_1322, \
                         msh_780, msh_992, msh_993, msi1_1321, msi1_1322, \
                         nsh_990 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1321[k] = pa_x[k] * msi0_1321[k]
                    + f_14 * msh_992[k]
                    - f_10 * pc_x[k] * msi1_1321[k];

        t_1322[k] = pa_x[k] * msi0_1322[k]
                    + f_13 * msh_993[k]
                    - f_10 * pc_x[k] * msi1_1322[k];

        t_1323[k] = f_12 * msh_780[k]
                    + f_3 * pc_z[k] * nsh_990[k];
    }

#pragma omp simd aligned(t_1324, t_1325, t_1326, pa_x, pc_x, pc_y, msi0_1325, msi0_1326, \
                         msh_803, msh_996, msh_997, msi1_1325, msi1_1326, \
                         nsh_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1324[k] = f_19 * msh_803[k]
                    + f_3 * pc_y[k] * nsh_992[k];

        t_1325[k] = pa_x[k] * msi0_1325[k]
                    + f_13 * msh_996[k]
                    - f_10 * pc_x[k] * msi1_1325[k];

        t_1326[k] = pa_x[k] * msi0_1326[k]
                    + f_12 * msh_997[k]
                    - f_10 * pc_x[k] * msi1_1326[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, pa_x, pc_x, pc_y, pc_z, msi0_1328, msh_783, \
                         msh_807, msh_999, msi1_1328, nsh_993, \
                         nsh_996 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = f_12 * msh_783[k]
                    + f_3 * pc_z[k] * nsh_993[k];

        t_1328[k] = pa_x[k] * msi0_1328[k]
                    + f_12 * msh_999[k]
                    - f_10 * pc_x[k] * msi1_1328[k];

        t_1329[k] = f_19 * msh_807[k]
                    + f_3 * pc_y[k] * nsh_996[k];
    }

#pragma omp simd aligned(t_1330, t_1331, t_1332, t_1333, pa_x, pc_x, msi0_1330, msh_1001, \
                         msh_1002, msh_1003, msh_1004, msi1_1330, nsh_1002, nsh_1003, \
                         nsh_1004 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1330[k] = pa_x[k] * msi0_1330[k]
                    + f_12 * msh_1001[k]
                    - f_10 * pc_x[k] * msi1_1330[k];

        t_1331[k] = f_11 * msh_1002[k]
                    + f_3 * pc_x[k] * nsh_1002[k];

        t_1332[k] = f_11 * msh_1003[k]
                    + f_3 * pc_x[k] * nsh_1003[k];

        t_1333[k] = f_11 * msh_1004[k]
                    + f_3 * pc_x[k] * nsh_1004[k];
    }

#pragma omp simd aligned(t_1334, t_1335, t_1336, t_1337, pa_x, pc_x, msi0_1337, msh_1005, \
                         msh_1006, msh_1007, msi1_1337, nsh_1005, nsh_1006, \
                         nsh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1334[k] = f_11 * msh_1005[k]
                    + f_3 * pc_x[k] * nsh_1005[k];

        t_1335[k] = f_11 * msh_1006[k]
                    + f_3 * pc_x[k] * nsh_1006[k];

        t_1336[k] = f_11 * msh_1007[k]
                    + f_3 * pc_x[k] * nsh_1007[k];

        t_1337[k] = pa_x[k] * msi0_1337[k]
                    - f_10 * pc_x[k] * msi1_1337[k];
    }

#pragma omp simd aligned(t_1338, t_1339, t_1340, t_1341, pa_x, pc_x, pc_z, msi0_1339, \
                         msi0_1340, msi0_1341, msh_792, msi1_1339, msi1_1340, msi1_1341, \
                         nsh_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1338[k] = f_12 * msh_792[k]
                    + f_3 * pc_z[k] * nsh_1002[k];

        t_1339[k] = pa_x[k] * msi0_1339[k]
                    - f_10 * pc_x[k] * msi1_1339[k];

        t_1340[k] = pa_x[k] * msi0_1340[k]
                    - f_10 * pc_x[k] * msi1_1340[k];

        t_1341[k] = pa_x[k] * msi0_1341[k]
                    - f_10 * pc_x[k] * msi1_1341[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, t_1345, pa_x, pc_x, pc_y, msi0_1343, \
                         msi0_1344, msh_818, msh_819, msh_1008, msi1_1343, msi1_1344, \
                         nsh_1007, nsh_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_19 * msh_818[k]
                    + f_3 * pc_y[k] * nsh_1007[k];

        t_1343[k] = pa_x[k] * msi0_1343[k]
                    - f_10 * pc_x[k] * msi1_1343[k];

        t_1344[k] = pa_x[k] * msi0_1344[k]
                    + f_20 * msh_1008[k]
                    - f_10 * pc_x[k] * msi1_1344[k];

        t_1345[k] = f_20 * msh_819[k]
                    + f_3 * pc_y[k] * nsh_1008[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pa_x, pc_x, pc_y, pc_z, msi0_1347, msh_798, \
                         msh_821, msh_1011, msi1_1347, nsh_1008, \
                         nsh_1010 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = f_13 * msh_798[k]
                    + f_3 * pc_z[k] * nsh_1008[k];

        t_1347[k] = pa_x[k] * msi0_1347[k]
                    + f_14 * msh_1011[k]
                    - f_10 * pc_x[k] * msi1_1347[k];

        t_1348[k] = f_20 * msh_821[k]
                    + f_3 * pc_y[k] * nsh_1010[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pa_x, pc_x, pc_z, msi0_1349, msi0_1350, \
                         msh_801, msh_1013, msh_1014, msi1_1349, msi1_1350, \
                         nsh_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = pa_x[k] * msi0_1349[k]
                    + f_14 * msh_1013[k]
                    - f_10 * pc_x[k] * msi1_1349[k];

        t_1350[k] = pa_x[k] * msi0_1350[k]
                    + f_13 * msh_1014[k]
                    - f_10 * pc_x[k] * msi1_1350[k];

        t_1351[k] = f_13 * msh_801[k]
                    + f_3 * pc_z[k] * nsh_1011[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, pa_x, pc_x, pc_y, msi0_1353, msi0_1354, \
                         msh_824, msh_1017, msh_1018, msi1_1353, msi1_1354, \
                         nsh_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_20 * msh_824[k]
                    + f_3 * pc_y[k] * nsh_1013[k];

        t_1353[k] = pa_x[k] * msi0_1353[k]
                    + f_13 * msh_1017[k]
                    - f_10 * pc_x[k] * msi1_1353[k];

        t_1354[k] = pa_x[k] * msi0_1354[k]
                    + f_12 * msh_1018[k]
                    - f_10 * pc_x[k] * msi1_1354[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pa_x, pc_x, pc_y, pc_z, msi0_1356, msh_804, \
                         msh_828, msh_1020, msi1_1356, nsh_1014, \
                         nsh_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_13 * msh_804[k]
                    + f_3 * pc_z[k] * nsh_1014[k];

        t_1356[k] = pa_x[k] * msi0_1356[k]
                    + f_12 * msh_1020[k]
                    - f_10 * pc_x[k] * msi1_1356[k];

        t_1357[k] = f_20 * msh_828[k]
                    + f_3 * pc_y[k] * nsh_1017[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, t_1361, pa_x, pc_x, msi0_1358, msh_1022, \
                         msh_1023, msh_1024, msh_1025, msi1_1358, nsh_1023, nsh_1024, \
                         nsh_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = pa_x[k] * msi0_1358[k]
                    + f_12 * msh_1022[k]
                    - f_10 * pc_x[k] * msi1_1358[k];

        t_1359[k] = f_11 * msh_1023[k]
                    + f_3 * pc_x[k] * nsh_1023[k];

        t_1360[k] = f_11 * msh_1024[k]
                    + f_3 * pc_x[k] * nsh_1024[k];

        t_1361[k] = f_11 * msh_1025[k]
                    + f_3 * pc_x[k] * nsh_1025[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, t_1365, pa_x, pc_x, msi0_1365, msh_1026, \
                         msh_1027, msh_1028, msi1_1365, nsh_1026, nsh_1027, \
                         nsh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = f_11 * msh_1026[k]
                    + f_3 * pc_x[k] * nsh_1026[k];

        t_1363[k] = f_11 * msh_1027[k]
                    + f_3 * pc_x[k] * nsh_1027[k];

        t_1364[k] = f_11 * msh_1028[k]
                    + f_3 * pc_x[k] * nsh_1028[k];

        t_1365[k] = pa_x[k] * msi0_1365[k]
                    - f_10 * pc_x[k] * msi1_1365[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, t_1369, pa_x, pc_x, pc_z, msi0_1367, \
                         msi0_1368, msi0_1369, msh_813, msi1_1367, msi1_1368, msi1_1369, \
                         nsh_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_13 * msh_813[k]
                    + f_3 * pc_z[k] * nsh_1023[k];

        t_1367[k] = pa_x[k] * msi0_1367[k]
                    - f_10 * pc_x[k] * msi1_1367[k];

        t_1368[k] = pa_x[k] * msi0_1368[k]
                    - f_10 * pc_x[k] * msi1_1368[k];

        t_1369[k] = pa_x[k] * msi0_1369[k]
                    - f_10 * pc_x[k] * msi1_1369[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, t_1373, pa_x, pc_x, pc_y, msi0_1371, \
                         msi0_1372, msh_839, msh_840, msh_1029, msi1_1371, msi1_1372, \
                         nsh_1028, nsh_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_20 * msh_839[k]
                    + f_3 * pc_y[k] * nsh_1028[k];

        t_1371[k] = pa_x[k] * msi0_1371[k]
                    - f_10 * pc_x[k] * msi1_1371[k];

        t_1372[k] = pa_x[k] * msi0_1372[k]
                    + f_20 * msh_1029[k]
                    - f_10 * pc_x[k] * msi1_1372[k];

        t_1373[k] = f_21 * msh_840[k]
                    + f_3 * pc_y[k] * nsh_1029[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, pa_x, pc_x, pc_y, pc_z, msi0_1375, msh_819, \
                         msh_842, msh_1032, msi1_1375, nsh_1029, \
                         nsh_1031 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_14 * msh_819[k]
                    + f_3 * pc_z[k] * nsh_1029[k];

        t_1375[k] = pa_x[k] * msi0_1375[k]
                    + f_14 * msh_1032[k]
                    - f_10 * pc_x[k] * msi1_1375[k];

        t_1376[k] = f_21 * msh_842[k]
                    + f_3 * pc_y[k] * nsh_1031[k];
    }

#pragma omp simd aligned(t_1377, t_1378, t_1379, pa_x, pc_x, pc_z, msi0_1377, msi0_1378, \
                         msh_822, msh_1034, msh_1035, msi1_1377, msi1_1378, \
                         nsh_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1377[k] = pa_x[k] * msi0_1377[k]
                    + f_14 * msh_1034[k]
                    - f_10 * pc_x[k] * msi1_1377[k];

        t_1378[k] = pa_x[k] * msi0_1378[k]
                    + f_13 * msh_1035[k]
                    - f_10 * pc_x[k] * msi1_1378[k];

        t_1379[k] = f_14 * msh_822[k]
                    + f_3 * pc_z[k] * nsh_1032[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, pa_x, pc_x, pc_y, msi0_1381, msi0_1382, \
                         msh_845, msh_1038, msh_1039, msi1_1381, msi1_1382, \
                         nsh_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_21 * msh_845[k]
                    + f_3 * pc_y[k] * nsh_1034[k];

        t_1381[k] = pa_x[k] * msi0_1381[k]
                    + f_13 * msh_1038[k]
                    - f_10 * pc_x[k] * msi1_1381[k];

        t_1382[k] = pa_x[k] * msi0_1382[k]
                    + f_12 * msh_1039[k]
                    - f_10 * pc_x[k] * msi1_1382[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msi0,
                                                           const size_t msh, const size_t msi1,
                                                           const size_t nsh, const size_t ncols,
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
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_1232 = buffer.data(msi0 + 1232);
    const auto *msi0_1237 = buffer.data(msi0 + 1237);
    const auto *msi0_1241 = buffer.data(msi0 + 1241);
    const auto *msi0_1246 = buffer.data(msi0 + 1246);
    const auto *msi0_1384 = buffer.data(msi0 + 1384);
    const auto *msi0_1386 = buffer.data(msi0 + 1386);
    const auto *msi0_1393 = buffer.data(msi0 + 1393);
    const auto *msi0_1395 = buffer.data(msi0 + 1395);
    const auto *msi0_1396 = buffer.data(msi0 + 1396);
    const auto *msi0_1397 = buffer.data(msi0 + 1397);
    const auto *msi0_1399 = buffer.data(msi0 + 1399);
    const auto *msi0_1400 = buffer.data(msi0 + 1400);
    const auto *msi0_1403 = buffer.data(msi0 + 1403);
    const auto *msi0_1405 = buffer.data(msi0 + 1405);
    const auto *msi0_1406 = buffer.data(msi0 + 1406);
    const auto *msi0_1409 = buffer.data(msi0 + 1409);
    const auto *msi0_1410 = buffer.data(msi0 + 1410);
    const auto *msi0_1412 = buffer.data(msi0 + 1412);
    const auto *msi0_1414 = buffer.data(msi0 + 1414);
    const auto *msi0_1421 = buffer.data(msi0 + 1421);
    const auto *msi0_1423 = buffer.data(msi0 + 1423);
    const auto *msi0_1424 = buffer.data(msi0 + 1424);
    const auto *msi0_1425 = buffer.data(msi0 + 1425);
    const auto *msi0_1427 = buffer.data(msi0 + 1427);
    const auto *msi0_1428 = buffer.data(msi0 + 1428);
    const auto *msi0_1431 = buffer.data(msi0 + 1431);
    const auto *msi0_1433 = buffer.data(msi0 + 1433);
    const auto *msi0_1434 = buffer.data(msi0 + 1434);
    const auto *msi0_1437 = buffer.data(msi0 + 1437);
    const auto *msi0_1438 = buffer.data(msi0 + 1438);
    const auto *msi0_1440 = buffer.data(msi0 + 1440);
    const auto *msi0_1442 = buffer.data(msi0 + 1442);
    const auto *msi0_1449 = buffer.data(msi0 + 1449);
    const auto *msi0_1451 = buffer.data(msi0 + 1451);
    const auto *msi0_1452 = buffer.data(msi0 + 1452);
    const auto *msi0_1453 = buffer.data(msi0 + 1453);
    const auto *msi0_1455 = buffer.data(msi0 + 1455);
    const auto *msi0_1456 = buffer.data(msi0 + 1456);
    const auto *msi0_1459 = buffer.data(msi0 + 1459);
    const auto *msi0_1461 = buffer.data(msi0 + 1461);
    const auto *msi0_1462 = buffer.data(msi0 + 1462);
    const auto *msi0_1465 = buffer.data(msi0 + 1465);
    const auto *msi0_1466 = buffer.data(msi0 + 1466);
    const auto *msi0_1468 = buffer.data(msi0 + 1468);
    const auto *msi0_1470 = buffer.data(msi0 + 1470);
    const auto *msi0_1477 = buffer.data(msi0 + 1477);
    const auto *msi0_1479 = buffer.data(msi0 + 1479);
    const auto *msi0_1480 = buffer.data(msi0 + 1480);
    const auto *msi0_1481 = buffer.data(msi0 + 1481);
    const auto *msi0_1483 = buffer.data(msi0 + 1483);
    const auto *msi0_1487 = buffer.data(msi0 + 1487);
    const auto *msi0_1490 = buffer.data(msi0 + 1490);
    const auto *msi0_1494 = buffer.data(msi0 + 1494);
    const auto *msi0_1496 = buffer.data(msi0 + 1496);

    const auto *msh_825 = buffer.data(msh + 825);
    const auto *msh_834 = buffer.data(msh + 834);
    const auto *msh_840 = buffer.data(msh + 840);
    const auto *msh_843 = buffer.data(msh + 843);
    const auto *msh_846 = buffer.data(msh + 846);
    const auto *msh_849 = buffer.data(msh + 849);
    const auto *msh_855 = buffer.data(msh + 855);
    const auto *msh_860 = buffer.data(msh + 860);
    const auto *msh_861 = buffer.data(msh + 861);
    const auto *msh_863 = buffer.data(msh + 863);
    const auto *msh_864 = buffer.data(msh + 864);
    const auto *msh_866 = buffer.data(msh + 866);
    const auto *msh_867 = buffer.data(msh + 867);
    const auto *msh_870 = buffer.data(msh + 870);
    const auto *msh_876 = buffer.data(msh + 876);
    const auto *msh_881 = buffer.data(msh + 881);
    const auto *msh_882 = buffer.data(msh + 882);
    const auto *msh_884 = buffer.data(msh + 884);
    const auto *msh_885 = buffer.data(msh + 885);
    const auto *msh_887 = buffer.data(msh + 887);
    const auto *msh_888 = buffer.data(msh + 888);
    const auto *msh_891 = buffer.data(msh + 891);
    const auto *msh_897 = buffer.data(msh + 897);
    const auto *msh_902 = buffer.data(msh + 902);
    const auto *msh_903 = buffer.data(msh + 903);
    const auto *msh_905 = buffer.data(msh + 905);
    const auto *msh_906 = buffer.data(msh + 906);
    const auto *msh_908 = buffer.data(msh + 908);
    const auto *msh_909 = buffer.data(msh + 909);
    const auto *msh_912 = buffer.data(msh + 912);
    const auto *msh_923 = buffer.data(msh + 923);
    const auto *msh_924 = buffer.data(msh + 924);
    const auto *msh_926 = buffer.data(msh + 926);
    const auto *msh_929 = buffer.data(msh + 929);
    const auto *msh_933 = buffer.data(msh + 933);
    const auto *msh_1041 = buffer.data(msh + 1041);
    const auto *msh_1043 = buffer.data(msh + 1043);
    const auto *msh_1044 = buffer.data(msh + 1044);
    const auto *msh_1045 = buffer.data(msh + 1045);
    const auto *msh_1046 = buffer.data(msh + 1046);
    const auto *msh_1047 = buffer.data(msh + 1047);
    const auto *msh_1048 = buffer.data(msh + 1048);
    const auto *msh_1049 = buffer.data(msh + 1049);
    const auto *msh_1050 = buffer.data(msh + 1050);
    const auto *msh_1053 = buffer.data(msh + 1053);
    const auto *msh_1055 = buffer.data(msh + 1055);
    const auto *msh_1056 = buffer.data(msh + 1056);
    const auto *msh_1059 = buffer.data(msh + 1059);
    const auto *msh_1060 = buffer.data(msh + 1060);
    const auto *msh_1062 = buffer.data(msh + 1062);
    const auto *msh_1064 = buffer.data(msh + 1064);
    const auto *msh_1065 = buffer.data(msh + 1065);
    const auto *msh_1066 = buffer.data(msh + 1066);
    const auto *msh_1067 = buffer.data(msh + 1067);
    const auto *msh_1068 = buffer.data(msh + 1068);
    const auto *msh_1069 = buffer.data(msh + 1069);
    const auto *msh_1070 = buffer.data(msh + 1070);
    const auto *msh_1071 = buffer.data(msh + 1071);
    const auto *msh_1074 = buffer.data(msh + 1074);
    const auto *msh_1076 = buffer.data(msh + 1076);
    const auto *msh_1077 = buffer.data(msh + 1077);
    const auto *msh_1080 = buffer.data(msh + 1080);
    const auto *msh_1081 = buffer.data(msh + 1081);
    const auto *msh_1083 = buffer.data(msh + 1083);
    const auto *msh_1085 = buffer.data(msh + 1085);
    const auto *msh_1086 = buffer.data(msh + 1086);
    const auto *msh_1087 = buffer.data(msh + 1087);
    const auto *msh_1088 = buffer.data(msh + 1088);
    const auto *msh_1089 = buffer.data(msh + 1089);
    const auto *msh_1090 = buffer.data(msh + 1090);
    const auto *msh_1091 = buffer.data(msh + 1091);
    const auto *msh_1092 = buffer.data(msh + 1092);
    const auto *msh_1095 = buffer.data(msh + 1095);
    const auto *msh_1097 = buffer.data(msh + 1097);
    const auto *msh_1098 = buffer.data(msh + 1098);
    const auto *msh_1101 = buffer.data(msh + 1101);
    const auto *msh_1102 = buffer.data(msh + 1102);
    const auto *msh_1104 = buffer.data(msh + 1104);
    const auto *msh_1106 = buffer.data(msh + 1106);
    const auto *msh_1107 = buffer.data(msh + 1107);
    const auto *msh_1108 = buffer.data(msh + 1108);
    const auto *msh_1109 = buffer.data(msh + 1109);
    const auto *msh_1110 = buffer.data(msh + 1110);
    const auto *msh_1111 = buffer.data(msh + 1111);
    const auto *msh_1112 = buffer.data(msh + 1112);
    const auto *msh_1116 = buffer.data(msh + 1116);
    const auto *msh_1119 = buffer.data(msh + 1119);
    const auto *msh_1123 = buffer.data(msh + 1123);
    const auto *msh_1125 = buffer.data(msh + 1125);
    const auto *msh_1128 = buffer.data(msh + 1128);
    const auto *msh_1129 = buffer.data(msh + 1129);
    const auto *msh_1130 = buffer.data(msh + 1130);

    const auto *msi1_1232 = buffer.data(msi1 + 1232);
    const auto *msi1_1237 = buffer.data(msi1 + 1237);
    const auto *msi1_1241 = buffer.data(msi1 + 1241);
    const auto *msi1_1246 = buffer.data(msi1 + 1246);
    const auto *msi1_1384 = buffer.data(msi1 + 1384);
    const auto *msi1_1386 = buffer.data(msi1 + 1386);
    const auto *msi1_1393 = buffer.data(msi1 + 1393);
    const auto *msi1_1395 = buffer.data(msi1 + 1395);
    const auto *msi1_1396 = buffer.data(msi1 + 1396);
    const auto *msi1_1397 = buffer.data(msi1 + 1397);
    const auto *msi1_1399 = buffer.data(msi1 + 1399);
    const auto *msi1_1400 = buffer.data(msi1 + 1400);
    const auto *msi1_1403 = buffer.data(msi1 + 1403);
    const auto *msi1_1405 = buffer.data(msi1 + 1405);
    const auto *msi1_1406 = buffer.data(msi1 + 1406);
    const auto *msi1_1409 = buffer.data(msi1 + 1409);
    const auto *msi1_1410 = buffer.data(msi1 + 1410);
    const auto *msi1_1412 = buffer.data(msi1 + 1412);
    const auto *msi1_1414 = buffer.data(msi1 + 1414);
    const auto *msi1_1421 = buffer.data(msi1 + 1421);
    const auto *msi1_1423 = buffer.data(msi1 + 1423);
    const auto *msi1_1424 = buffer.data(msi1 + 1424);
    const auto *msi1_1425 = buffer.data(msi1 + 1425);
    const auto *msi1_1427 = buffer.data(msi1 + 1427);
    const auto *msi1_1428 = buffer.data(msi1 + 1428);
    const auto *msi1_1431 = buffer.data(msi1 + 1431);
    const auto *msi1_1433 = buffer.data(msi1 + 1433);
    const auto *msi1_1434 = buffer.data(msi1 + 1434);
    const auto *msi1_1437 = buffer.data(msi1 + 1437);
    const auto *msi1_1438 = buffer.data(msi1 + 1438);
    const auto *msi1_1440 = buffer.data(msi1 + 1440);
    const auto *msi1_1442 = buffer.data(msi1 + 1442);
    const auto *msi1_1449 = buffer.data(msi1 + 1449);
    const auto *msi1_1451 = buffer.data(msi1 + 1451);
    const auto *msi1_1452 = buffer.data(msi1 + 1452);
    const auto *msi1_1453 = buffer.data(msi1 + 1453);
    const auto *msi1_1455 = buffer.data(msi1 + 1455);
    const auto *msi1_1456 = buffer.data(msi1 + 1456);
    const auto *msi1_1459 = buffer.data(msi1 + 1459);
    const auto *msi1_1461 = buffer.data(msi1 + 1461);
    const auto *msi1_1462 = buffer.data(msi1 + 1462);
    const auto *msi1_1465 = buffer.data(msi1 + 1465);
    const auto *msi1_1466 = buffer.data(msi1 + 1466);
    const auto *msi1_1468 = buffer.data(msi1 + 1468);
    const auto *msi1_1470 = buffer.data(msi1 + 1470);
    const auto *msi1_1477 = buffer.data(msi1 + 1477);
    const auto *msi1_1479 = buffer.data(msi1 + 1479);
    const auto *msi1_1480 = buffer.data(msi1 + 1480);
    const auto *msi1_1481 = buffer.data(msi1 + 1481);
    const auto *msi1_1483 = buffer.data(msi1 + 1483);
    const auto *msi1_1487 = buffer.data(msi1 + 1487);
    const auto *msi1_1490 = buffer.data(msi1 + 1490);
    const auto *msi1_1494 = buffer.data(msi1 + 1494);
    const auto *msi1_1496 = buffer.data(msi1 + 1496);

    const auto *nsh_1035 = buffer.data(nsh + 1035);
    const auto *nsh_1038 = buffer.data(nsh + 1038);
    const auto *nsh_1044 = buffer.data(nsh + 1044);
    const auto *nsh_1045 = buffer.data(nsh + 1045);
    const auto *nsh_1046 = buffer.data(nsh + 1046);
    const auto *nsh_1047 = buffer.data(nsh + 1047);
    const auto *nsh_1048 = buffer.data(nsh + 1048);
    const auto *nsh_1049 = buffer.data(nsh + 1049);
    const auto *nsh_1050 = buffer.data(nsh + 1050);
    const auto *nsh_1052 = buffer.data(nsh + 1052);
    const auto *nsh_1053 = buffer.data(nsh + 1053);
    const auto *nsh_1055 = buffer.data(nsh + 1055);
    const auto *nsh_1056 = buffer.data(nsh + 1056);
    const auto *nsh_1059 = buffer.data(nsh + 1059);
    const auto *nsh_1065 = buffer.data(nsh + 1065);
    const auto *nsh_1066 = buffer.data(nsh + 1066);
    const auto *nsh_1067 = buffer.data(nsh + 1067);
    const auto *nsh_1068 = buffer.data(nsh + 1068);
    const auto *nsh_1069 = buffer.data(nsh + 1069);
    const auto *nsh_1070 = buffer.data(nsh + 1070);
    const auto *nsh_1071 = buffer.data(nsh + 1071);
    const auto *nsh_1073 = buffer.data(nsh + 1073);
    const auto *nsh_1074 = buffer.data(nsh + 1074);
    const auto *nsh_1076 = buffer.data(nsh + 1076);
    const auto *nsh_1077 = buffer.data(nsh + 1077);
    const auto *nsh_1080 = buffer.data(nsh + 1080);
    const auto *nsh_1086 = buffer.data(nsh + 1086);
    const auto *nsh_1087 = buffer.data(nsh + 1087);
    const auto *nsh_1088 = buffer.data(nsh + 1088);
    const auto *nsh_1089 = buffer.data(nsh + 1089);
    const auto *nsh_1090 = buffer.data(nsh + 1090);
    const auto *nsh_1091 = buffer.data(nsh + 1091);
    const auto *nsh_1092 = buffer.data(nsh + 1092);
    const auto *nsh_1094 = buffer.data(nsh + 1094);
    const auto *nsh_1095 = buffer.data(nsh + 1095);
    const auto *nsh_1097 = buffer.data(nsh + 1097);
    const auto *nsh_1098 = buffer.data(nsh + 1098);
    const auto *nsh_1101 = buffer.data(nsh + 1101);
    const auto *nsh_1107 = buffer.data(nsh + 1107);
    const auto *nsh_1108 = buffer.data(nsh + 1108);
    const auto *nsh_1109 = buffer.data(nsh + 1109);
    const auto *nsh_1110 = buffer.data(nsh + 1110);
    const auto *nsh_1111 = buffer.data(nsh + 1111);
    const auto *nsh_1112 = buffer.data(nsh + 1112);
    const auto *nsh_1113 = buffer.data(nsh + 1113);
    const auto *nsh_1115 = buffer.data(nsh + 1115);
    const auto *nsh_1116 = buffer.data(nsh + 1116);
    const auto *nsh_1118 = buffer.data(nsh + 1118);
    const auto *nsh_1119 = buffer.data(nsh + 1119);
    const auto *nsh_1122 = buffer.data(nsh + 1122);
    const auto *nsh_1128 = buffer.data(nsh + 1128);
    const auto *nsh_1129 = buffer.data(nsh + 1129);
    const auto *nsh_1130 = buffer.data(nsh + 1130);

#pragma omp simd aligned(t_1383, t_1384, t_1385, pa_x, pc_x, pc_y, pc_z, msi0_1384, msh_825, \
                         msh_849, msh_1041, msi1_1384, nsh_1035, \
                         nsh_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1383[k] = f_14 * msh_825[k]
                    + f_3 * pc_z[k] * nsh_1035[k];

        t_1384[k] = pa_x[k] * msi0_1384[k]
                    + f_12 * msh_1041[k]
                    - f_10 * pc_x[k] * msi1_1384[k];

        t_1385[k] = f_21 * msh_849[k]
                    + f_3 * pc_y[k] * nsh_1038[k];
    }

#pragma omp simd aligned(t_1386, t_1387, t_1388, t_1389, pa_x, pc_x, msi0_1386, msh_1043, \
                         msh_1044, msh_1045, msh_1046, msi1_1386, nsh_1044, nsh_1045, \
                         nsh_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1386[k] = pa_x[k] * msi0_1386[k]
                    + f_12 * msh_1043[k]
                    - f_10 * pc_x[k] * msi1_1386[k];

        t_1387[k] = f_11 * msh_1044[k]
                    + f_3 * pc_x[k] * nsh_1044[k];

        t_1388[k] = f_11 * msh_1045[k]
                    + f_3 * pc_x[k] * nsh_1045[k];

        t_1389[k] = f_11 * msh_1046[k]
                    + f_3 * pc_x[k] * nsh_1046[k];
    }

#pragma omp simd aligned(t_1390, t_1391, t_1392, t_1393, pa_x, pc_x, msi0_1393, msh_1047, \
                         msh_1048, msh_1049, msi1_1393, nsh_1047, nsh_1048, \
                         nsh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1390[k] = f_11 * msh_1047[k]
                    + f_3 * pc_x[k] * nsh_1047[k];

        t_1391[k] = f_11 * msh_1048[k]
                    + f_3 * pc_x[k] * nsh_1048[k];

        t_1392[k] = f_11 * msh_1049[k]
                    + f_3 * pc_x[k] * nsh_1049[k];

        t_1393[k] = pa_x[k] * msi0_1393[k]
                    - f_10 * pc_x[k] * msi1_1393[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, t_1397, pa_x, pc_x, pc_z, msi0_1395, \
                         msi0_1396, msi0_1397, msh_834, msi1_1395, msi1_1396, msi1_1397, \
                         nsh_1044 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_14 * msh_834[k]
                    + f_3 * pc_z[k] * nsh_1044[k];

        t_1395[k] = pa_x[k] * msi0_1395[k]
                    - f_10 * pc_x[k] * msi1_1395[k];

        t_1396[k] = pa_x[k] * msi0_1396[k]
                    - f_10 * pc_x[k] * msi1_1396[k];

        t_1397[k] = pa_x[k] * msi0_1397[k]
                    - f_10 * pc_x[k] * msi1_1397[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, t_1401, pa_x, pc_x, pc_y, msi0_1399, \
                         msi0_1400, msh_860, msh_861, msh_1050, msi1_1399, msi1_1400, \
                         nsh_1049, nsh_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_21 * msh_860[k]
                    + f_3 * pc_y[k] * nsh_1049[k];

        t_1399[k] = pa_x[k] * msi0_1399[k]
                    - f_10 * pc_x[k] * msi1_1399[k];

        t_1400[k] = pa_x[k] * msi0_1400[k]
                    + f_20 * msh_1050[k]
                    - f_10 * pc_x[k] * msi1_1400[k];

        t_1401[k] = f_14 * msh_861[k]
                    + f_3 * pc_y[k] * nsh_1050[k];
    }

#pragma omp simd aligned(t_1402, t_1403, t_1404, pa_x, pc_x, pc_y, pc_z, msi0_1403, msh_840, \
                         msh_863, msh_1053, msi1_1403, nsh_1050, \
                         nsh_1052 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1402[k] = f_21 * msh_840[k]
                    + f_3 * pc_z[k] * nsh_1050[k];

        t_1403[k] = pa_x[k] * msi0_1403[k]
                    + f_14 * msh_1053[k]
                    - f_10 * pc_x[k] * msi1_1403[k];

        t_1404[k] = f_14 * msh_863[k]
                    + f_3 * pc_y[k] * nsh_1052[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, pa_x, pc_x, pc_z, msi0_1405, msi0_1406, \
                         msh_843, msh_1055, msh_1056, msi1_1405, msi1_1406, \
                         nsh_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = pa_x[k] * msi0_1405[k]
                    + f_14 * msh_1055[k]
                    - f_10 * pc_x[k] * msi1_1405[k];

        t_1406[k] = pa_x[k] * msi0_1406[k]
                    + f_13 * msh_1056[k]
                    - f_10 * pc_x[k] * msi1_1406[k];

        t_1407[k] = f_21 * msh_843[k]
                    + f_3 * pc_z[k] * nsh_1053[k];
    }

#pragma omp simd aligned(t_1408, t_1409, t_1410, pa_x, pc_x, pc_y, msi0_1409, msi0_1410, \
                         msh_866, msh_1059, msh_1060, msi1_1409, msi1_1410, \
                         nsh_1055 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1408[k] = f_14 * msh_866[k]
                    + f_3 * pc_y[k] * nsh_1055[k];

        t_1409[k] = pa_x[k] * msi0_1409[k]
                    + f_13 * msh_1059[k]
                    - f_10 * pc_x[k] * msi1_1409[k];

        t_1410[k] = pa_x[k] * msi0_1410[k]
                    + f_12 * msh_1060[k]
                    - f_10 * pc_x[k] * msi1_1410[k];
    }

#pragma omp simd aligned(t_1411, t_1412, t_1413, pa_x, pc_x, pc_y, pc_z, msi0_1412, msh_846, \
                         msh_870, msh_1062, msi1_1412, nsh_1056, \
                         nsh_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1411[k] = f_21 * msh_846[k]
                    + f_3 * pc_z[k] * nsh_1056[k];

        t_1412[k] = pa_x[k] * msi0_1412[k]
                    + f_12 * msh_1062[k]
                    - f_10 * pc_x[k] * msi1_1412[k];

        t_1413[k] = f_14 * msh_870[k]
                    + f_3 * pc_y[k] * nsh_1059[k];
    }

#pragma omp simd aligned(t_1414, t_1415, t_1416, t_1417, pa_x, pc_x, msi0_1414, msh_1064, \
                         msh_1065, msh_1066, msh_1067, msi1_1414, nsh_1065, nsh_1066, \
                         nsh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1414[k] = pa_x[k] * msi0_1414[k]
                    + f_12 * msh_1064[k]
                    - f_10 * pc_x[k] * msi1_1414[k];

        t_1415[k] = f_11 * msh_1065[k]
                    + f_3 * pc_x[k] * nsh_1065[k];

        t_1416[k] = f_11 * msh_1066[k]
                    + f_3 * pc_x[k] * nsh_1066[k];

        t_1417[k] = f_11 * msh_1067[k]
                    + f_3 * pc_x[k] * nsh_1067[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, t_1421, pa_x, pc_x, msi0_1421, msh_1068, \
                         msh_1069, msh_1070, msi1_1421, nsh_1068, nsh_1069, \
                         nsh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_11 * msh_1068[k]
                    + f_3 * pc_x[k] * nsh_1068[k];

        t_1419[k] = f_11 * msh_1069[k]
                    + f_3 * pc_x[k] * nsh_1069[k];

        t_1420[k] = f_11 * msh_1070[k]
                    + f_3 * pc_x[k] * nsh_1070[k];

        t_1421[k] = pa_x[k] * msi0_1421[k]
                    - f_10 * pc_x[k] * msi1_1421[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, t_1425, pa_x, pc_x, pc_z, msi0_1423, \
                         msi0_1424, msi0_1425, msh_855, msi1_1423, msi1_1424, msi1_1425, \
                         nsh_1065 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_21 * msh_855[k]
                    + f_3 * pc_z[k] * nsh_1065[k];

        t_1423[k] = pa_x[k] * msi0_1423[k]
                    - f_10 * pc_x[k] * msi1_1423[k];

        t_1424[k] = pa_x[k] * msi0_1424[k]
                    - f_10 * pc_x[k] * msi1_1424[k];

        t_1425[k] = pa_x[k] * msi0_1425[k]
                    - f_10 * pc_x[k] * msi1_1425[k];
    }

#pragma omp simd aligned(t_1426, t_1427, t_1428, t_1429, pa_x, pc_x, pc_y, msi0_1427, \
                         msi0_1428, msh_881, msh_882, msh_1071, msi1_1427, msi1_1428, \
                         nsh_1070, nsh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1426[k] = f_14 * msh_881[k]
                    + f_3 * pc_y[k] * nsh_1070[k];

        t_1427[k] = pa_x[k] * msi0_1427[k]
                    - f_10 * pc_x[k] * msi1_1427[k];

        t_1428[k] = pa_x[k] * msi0_1428[k]
                    + f_20 * msh_1071[k]
                    - f_10 * pc_x[k] * msi1_1428[k];

        t_1429[k] = f_13 * msh_882[k]
                    + f_3 * pc_y[k] * nsh_1071[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, pa_x, pc_x, pc_y, pc_z, msi0_1431, msh_861, \
                         msh_884, msh_1074, msi1_1431, nsh_1071, \
                         nsh_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_20 * msh_861[k]
                    + f_3 * pc_z[k] * nsh_1071[k];

        t_1431[k] = pa_x[k] * msi0_1431[k]
                    + f_14 * msh_1074[k]
                    - f_10 * pc_x[k] * msi1_1431[k];

        t_1432[k] = f_13 * msh_884[k]
                    + f_3 * pc_y[k] * nsh_1073[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pa_x, pc_x, pc_z, msi0_1433, msi0_1434, \
                         msh_864, msh_1076, msh_1077, msi1_1433, msi1_1434, \
                         nsh_1074 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = pa_x[k] * msi0_1433[k]
                    + f_14 * msh_1076[k]
                    - f_10 * pc_x[k] * msi1_1433[k];

        t_1434[k] = pa_x[k] * msi0_1434[k]
                    + f_13 * msh_1077[k]
                    - f_10 * pc_x[k] * msi1_1434[k];

        t_1435[k] = f_20 * msh_864[k]
                    + f_3 * pc_z[k] * nsh_1074[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, pa_x, pc_x, pc_y, msi0_1437, msi0_1438, \
                         msh_887, msh_1080, msh_1081, msi1_1437, msi1_1438, \
                         nsh_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_13 * msh_887[k]
                    + f_3 * pc_y[k] * nsh_1076[k];

        t_1437[k] = pa_x[k] * msi0_1437[k]
                    + f_13 * msh_1080[k]
                    - f_10 * pc_x[k] * msi1_1437[k];

        t_1438[k] = pa_x[k] * msi0_1438[k]
                    + f_12 * msh_1081[k]
                    - f_10 * pc_x[k] * msi1_1438[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, pa_x, pc_x, pc_y, pc_z, msi0_1440, msh_867, \
                         msh_891, msh_1083, msi1_1440, nsh_1077, \
                         nsh_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = f_20 * msh_867[k]
                    + f_3 * pc_z[k] * nsh_1077[k];

        t_1440[k] = pa_x[k] * msi0_1440[k]
                    + f_12 * msh_1083[k]
                    - f_10 * pc_x[k] * msi1_1440[k];

        t_1441[k] = f_13 * msh_891[k]
                    + f_3 * pc_y[k] * nsh_1080[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, t_1445, pa_x, pc_x, msi0_1442, msh_1085, \
                         msh_1086, msh_1087, msh_1088, msi1_1442, nsh_1086, nsh_1087, \
                         nsh_1088 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = pa_x[k] * msi0_1442[k]
                    + f_12 * msh_1085[k]
                    - f_10 * pc_x[k] * msi1_1442[k];

        t_1443[k] = f_11 * msh_1086[k]
                    + f_3 * pc_x[k] * nsh_1086[k];

        t_1444[k] = f_11 * msh_1087[k]
                    + f_3 * pc_x[k] * nsh_1087[k];

        t_1445[k] = f_11 * msh_1088[k]
                    + f_3 * pc_x[k] * nsh_1088[k];
    }

#pragma omp simd aligned(t_1446, t_1447, t_1448, t_1449, pa_x, pc_x, msi0_1449, msh_1089, \
                         msh_1090, msh_1091, msi1_1449, nsh_1089, nsh_1090, \
                         nsh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1446[k] = f_11 * msh_1089[k]
                    + f_3 * pc_x[k] * nsh_1089[k];

        t_1447[k] = f_11 * msh_1090[k]
                    + f_3 * pc_x[k] * nsh_1090[k];

        t_1448[k] = f_11 * msh_1091[k]
                    + f_3 * pc_x[k] * nsh_1091[k];

        t_1449[k] = pa_x[k] * msi0_1449[k]
                    - f_10 * pc_x[k] * msi1_1449[k];
    }

#pragma omp simd aligned(t_1450, t_1451, t_1452, t_1453, pa_x, pc_x, pc_z, msi0_1451, \
                         msi0_1452, msi0_1453, msh_876, msi1_1451, msi1_1452, msi1_1453, \
                         nsh_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1450[k] = f_20 * msh_876[k]
                    + f_3 * pc_z[k] * nsh_1086[k];

        t_1451[k] = pa_x[k] * msi0_1451[k]
                    - f_10 * pc_x[k] * msi1_1451[k];

        t_1452[k] = pa_x[k] * msi0_1452[k]
                    - f_10 * pc_x[k] * msi1_1452[k];

        t_1453[k] = pa_x[k] * msi0_1453[k]
                    - f_10 * pc_x[k] * msi1_1453[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, t_1457, pa_x, pc_x, pc_y, msi0_1455, \
                         msi0_1456, msh_902, msh_903, msh_1092, msi1_1455, msi1_1456, \
                         nsh_1091, nsh_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_13 * msh_902[k]
                    + f_3 * pc_y[k] * nsh_1091[k];

        t_1455[k] = pa_x[k] * msi0_1455[k]
                    - f_10 * pc_x[k] * msi1_1455[k];

        t_1456[k] = pa_x[k] * msi0_1456[k]
                    + f_20 * msh_1092[k]
                    - f_10 * pc_x[k] * msi1_1456[k];

        t_1457[k] = f_12 * msh_903[k]
                    + f_3 * pc_y[k] * nsh_1092[k];
    }

#pragma omp simd aligned(t_1458, t_1459, t_1460, pa_x, pc_x, pc_y, pc_z, msi0_1459, msh_882, \
                         msh_905, msh_1095, msi1_1459, nsh_1092, \
                         nsh_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1458[k] = f_19 * msh_882[k]
                    + f_3 * pc_z[k] * nsh_1092[k];

        t_1459[k] = pa_x[k] * msi0_1459[k]
                    + f_14 * msh_1095[k]
                    - f_10 * pc_x[k] * msi1_1459[k];

        t_1460[k] = f_12 * msh_905[k]
                    + f_3 * pc_y[k] * nsh_1094[k];
    }

#pragma omp simd aligned(t_1461, t_1462, t_1463, pa_x, pc_x, pc_z, msi0_1461, msi0_1462, \
                         msh_885, msh_1097, msh_1098, msi1_1461, msi1_1462, \
                         nsh_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1461[k] = pa_x[k] * msi0_1461[k]
                    + f_14 * msh_1097[k]
                    - f_10 * pc_x[k] * msi1_1461[k];

        t_1462[k] = pa_x[k] * msi0_1462[k]
                    + f_13 * msh_1098[k]
                    - f_10 * pc_x[k] * msi1_1462[k];

        t_1463[k] = f_19 * msh_885[k]
                    + f_3 * pc_z[k] * nsh_1095[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, pa_x, pc_x, pc_y, msi0_1465, msi0_1466, \
                         msh_908, msh_1101, msh_1102, msi1_1465, msi1_1466, \
                         nsh_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_12 * msh_908[k]
                    + f_3 * pc_y[k] * nsh_1097[k];

        t_1465[k] = pa_x[k] * msi0_1465[k]
                    + f_13 * msh_1101[k]
                    - f_10 * pc_x[k] * msi1_1465[k];

        t_1466[k] = pa_x[k] * msi0_1466[k]
                    + f_12 * msh_1102[k]
                    - f_10 * pc_x[k] * msi1_1466[k];
    }

#pragma omp simd aligned(t_1467, t_1468, t_1469, pa_x, pc_x, pc_y, pc_z, msi0_1468, msh_888, \
                         msh_912, msh_1104, msi1_1468, nsh_1098, \
                         nsh_1101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1467[k] = f_19 * msh_888[k]
                    + f_3 * pc_z[k] * nsh_1098[k];

        t_1468[k] = pa_x[k] * msi0_1468[k]
                    + f_12 * msh_1104[k]
                    - f_10 * pc_x[k] * msi1_1468[k];

        t_1469[k] = f_12 * msh_912[k]
                    + f_3 * pc_y[k] * nsh_1101[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, pa_x, pc_x, msi0_1470, msh_1106, \
                         msh_1107, msh_1108, msh_1109, msi1_1470, nsh_1107, nsh_1108, \
                         nsh_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = pa_x[k] * msi0_1470[k]
                    + f_12 * msh_1106[k]
                    - f_10 * pc_x[k] * msi1_1470[k];

        t_1471[k] = f_11 * msh_1107[k]
                    + f_3 * pc_x[k] * nsh_1107[k];

        t_1472[k] = f_11 * msh_1108[k]
                    + f_3 * pc_x[k] * nsh_1108[k];

        t_1473[k] = f_11 * msh_1109[k]
                    + f_3 * pc_x[k] * nsh_1109[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, pa_x, pc_x, msi0_1477, msh_1110, \
                         msh_1111, msh_1112, msi1_1477, nsh_1110, nsh_1111, \
                         nsh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_11 * msh_1110[k]
                    + f_3 * pc_x[k] * nsh_1110[k];

        t_1475[k] = f_11 * msh_1111[k]
                    + f_3 * pc_x[k] * nsh_1111[k];

        t_1476[k] = f_11 * msh_1112[k]
                    + f_3 * pc_x[k] * nsh_1112[k];

        t_1477[k] = pa_x[k] * msi0_1477[k]
                    - f_10 * pc_x[k] * msi1_1477[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, t_1481, pa_x, pc_x, pc_z, msi0_1479, \
                         msi0_1480, msi0_1481, msh_897, msi1_1479, msi1_1480, msi1_1481, \
                         nsh_1107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_19 * msh_897[k]
                    + f_3 * pc_z[k] * nsh_1107[k];

        t_1479[k] = pa_x[k] * msi0_1479[k]
                    - f_10 * pc_x[k] * msi1_1479[k];

        t_1480[k] = pa_x[k] * msi0_1480[k]
                    - f_10 * pc_x[k] * msi1_1480[k];

        t_1481[k] = pa_x[k] * msi0_1481[k]
                    - f_10 * pc_x[k] * msi1_1481[k];
    }

#pragma omp simd aligned(t_1482, t_1483, t_1484, t_1485, pa_x, pa_y, pc_x, pc_y, msi0_1232, \
                         msi0_1483, msh_923, msh_924, msi1_1232, msi1_1483, nsh_1112, \
                         nsh_1113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1482[k] = f_12 * msh_923[k]
                    + f_3 * pc_y[k] * nsh_1112[k];

        t_1483[k] = pa_x[k] * msi0_1483[k]
                    - f_10 * pc_x[k] * msi1_1483[k];

        t_1484[k] = pa_y[k] * msi0_1232[k]
                    - f_10 * pc_y[k] * msi1_1232[k];

        t_1485[k] = f_11 * msh_924[k]
                    + f_3 * pc_y[k] * nsh_1113[k];
    }

#pragma omp simd aligned(t_1486, t_1487, t_1488, pa_x, pc_x, pc_y, pc_z, msi0_1487, msh_903, \
                         msh_926, msh_1116, msi1_1487, nsh_1113, \
                         nsh_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1486[k] = f_18 * msh_903[k]
                    + f_3 * pc_z[k] * nsh_1113[k];

        t_1487[k] = pa_x[k] * msi0_1487[k]
                    + f_14 * msh_1116[k]
                    - f_10 * pc_x[k] * msi1_1487[k];

        t_1488[k] = f_11 * msh_926[k]
                    + f_3 * pc_y[k] * nsh_1115[k];
    }

#pragma omp simd aligned(t_1489, t_1490, t_1491, pa_x, pa_y, pc_x, pc_y, pc_z, msi0_1237, \
                         msi0_1490, msh_906, msh_1119, msi1_1237, msi1_1490, \
                         nsh_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1489[k] = pa_y[k] * msi0_1237[k]
                    - f_10 * pc_y[k] * msi1_1237[k];

        t_1490[k] = pa_x[k] * msi0_1490[k]
                    + f_13 * msh_1119[k]
                    - f_10 * pc_x[k] * msi1_1490[k];

        t_1491[k] = f_18 * msh_906[k]
                    + f_3 * pc_z[k] * nsh_1116[k];
    }

#pragma omp simd aligned(t_1492, t_1493, t_1494, pa_x, pa_y, pc_x, pc_y, msi0_1241, msi0_1494, \
                         msh_929, msh_1123, msi1_1241, msi1_1494, \
                         nsh_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1492[k] = f_11 * msh_929[k]
                    + f_3 * pc_y[k] * nsh_1118[k];

        t_1493[k] = pa_y[k] * msi0_1241[k]
                    - f_10 * pc_y[k] * msi1_1241[k];

        t_1494[k] = pa_x[k] * msi0_1494[k]
                    + f_12 * msh_1123[k]
                    - f_10 * pc_x[k] * msi1_1494[k];
    }

#pragma omp simd aligned(t_1495, t_1496, t_1497, pa_x, pc_x, pc_y, pc_z, msi0_1496, msh_909, \
                         msh_933, msh_1125, msi1_1496, nsh_1119, \
                         nsh_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1495[k] = f_18 * msh_909[k]
                    + f_3 * pc_z[k] * nsh_1119[k];

        t_1496[k] = pa_x[k] * msi0_1496[k]
                    + f_12 * msh_1125[k]
                    - f_10 * pc_x[k] * msi1_1496[k];

        t_1497[k] = f_11 * msh_933[k]
                    + f_3 * pc_y[k] * nsh_1122[k];
    }

#pragma omp simd aligned(t_1498, t_1499, t_1500, t_1501, pa_y, pc_x, pc_y, msi0_1246, \
                         msh_1128, msh_1129, msh_1130, msi1_1246, nsh_1128, nsh_1129, \
                         nsh_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1498[k] = pa_y[k] * msi0_1246[k]
                    - f_10 * pc_y[k] * msi1_1246[k];

        t_1499[k] = f_11 * msh_1128[k]
                    + f_3 * pc_x[k] * nsh_1128[k];

        t_1500[k] = f_11 * msh_1129[k]
                    + f_3 * pc_x[k] * nsh_1129[k];

        t_1501[k] = f_11 * msh_1130[k]
                    + f_3 * pc_x[k] * nsh_1130[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msi0,
                                                           const size_t msh, const size_t msi1,
                                                           const size_t nsg0, const size_t nsg1,
                                                           const size_t nsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.0 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_1260 = buffer.data(msi0 + 1260);
    const auto *msi0_1261 = buffer.data(msi0 + 1261);
    const auto *msi0_1263 = buffer.data(msi0 + 1263);
    const auto *msi0_1266 = buffer.data(msi0 + 1266);
    const auto *msi0_1270 = buffer.data(msi0 + 1270);
    const auto *msi0_1281 = buffer.data(msi0 + 1281);
    const auto *msi0_1283 = buffer.data(msi0 + 1283);
    const auto *msi0_1284 = buffer.data(msi0 + 1284);
    const auto *msi0_1285 = buffer.data(msi0 + 1285);
    const auto *msi0_1505 = buffer.data(msi0 + 1505);
    const auto *msi0_1507 = buffer.data(msi0 + 1507);
    const auto *msi0_1508 = buffer.data(msi0 + 1508);
    const auto *msi0_1509 = buffer.data(msi0 + 1509);
    const auto *msi0_1511 = buffer.data(msi0 + 1511);
    const auto *msi0_1512 = buffer.data(msi0 + 1512);
    const auto *msi0_1517 = buffer.data(msi0 + 1517);
    const auto *msi0_1521 = buffer.data(msi0 + 1521);
    const auto *msi0_1526 = buffer.data(msi0 + 1526);
    const auto *msi0_1533 = buffer.data(msi0 + 1533);
    const auto *msi0_1534 = buffer.data(msi0 + 1534);
    const auto *msi0_1535 = buffer.data(msi0 + 1535);
    const auto *msi0_1536 = buffer.data(msi0 + 1536);
    const auto *msi0_1537 = buffer.data(msi0 + 1537);
    const auto *msi0_1539 = buffer.data(msi0 + 1539);

    const auto *msh_918 = buffer.data(msh + 918);
    const auto *msh_924 = buffer.data(msh + 924);
    const auto *msh_944 = buffer.data(msh + 944);
    const auto *msh_960 = buffer.data(msh + 960);
    const auto *msh_961 = buffer.data(msh + 961);
    const auto *msh_962 = buffer.data(msh + 962);
    const auto *msh_963 = buffer.data(msh + 963);
    const auto *msh_965 = buffer.data(msh + 965);
    const auto *msh_981 = buffer.data(msh + 981);
    const auto *msh_986 = buffer.data(msh + 986);
    const auto *msh_1002 = buffer.data(msh + 1002);
    const auto *msh_1004 = buffer.data(msh + 1004);
    const auto *msh_1005 = buffer.data(msh + 1005);
    const auto *msh_1006 = buffer.data(msh + 1006);
    const auto *msh_1007 = buffer.data(msh + 1007);
    const auto *msh_1131 = buffer.data(msh + 1131);
    const auto *msh_1132 = buffer.data(msh + 1132);
    const auto *msh_1133 = buffer.data(msh + 1133);
    const auto *msh_1134 = buffer.data(msh + 1134);
    const auto *msh_1139 = buffer.data(msh + 1139);
    const auto *msh_1143 = buffer.data(msh + 1143);
    const auto *msh_1148 = buffer.data(msh + 1148);
    const auto *msh_1149 = buffer.data(msh + 1149);
    const auto *msh_1150 = buffer.data(msh + 1150);
    const auto *msh_1151 = buffer.data(msh + 1151);
    const auto *msh_1152 = buffer.data(msh + 1152);
    const auto *msh_1154 = buffer.data(msh + 1154);

    const auto *msi1_1260 = buffer.data(msi1 + 1260);
    const auto *msi1_1261 = buffer.data(msi1 + 1261);
    const auto *msi1_1263 = buffer.data(msi1 + 1263);
    const auto *msi1_1266 = buffer.data(msi1 + 1266);
    const auto *msi1_1270 = buffer.data(msi1 + 1270);
    const auto *msi1_1281 = buffer.data(msi1 + 1281);
    const auto *msi1_1283 = buffer.data(msi1 + 1283);
    const auto *msi1_1284 = buffer.data(msi1 + 1284);
    const auto *msi1_1285 = buffer.data(msi1 + 1285);
    const auto *msi1_1505 = buffer.data(msi1 + 1505);
    const auto *msi1_1507 = buffer.data(msi1 + 1507);
    const auto *msi1_1508 = buffer.data(msi1 + 1508);
    const auto *msi1_1509 = buffer.data(msi1 + 1509);
    const auto *msi1_1511 = buffer.data(msi1 + 1511);
    const auto *msi1_1512 = buffer.data(msi1 + 1512);
    const auto *msi1_1517 = buffer.data(msi1 + 1517);
    const auto *msi1_1521 = buffer.data(msi1 + 1521);
    const auto *msi1_1526 = buffer.data(msi1 + 1526);
    const auto *msi1_1533 = buffer.data(msi1 + 1533);
    const auto *msi1_1534 = buffer.data(msi1 + 1534);
    const auto *msi1_1535 = buffer.data(msi1 + 1535);
    const auto *msi1_1536 = buffer.data(msi1 + 1536);
    const auto *msi1_1537 = buffer.data(msi1 + 1537);
    const auto *msi1_1539 = buffer.data(msi1 + 1539);

    const auto *nsg0_810 = buffer.data(nsg0 + 810);
    const auto *nsg0_811 = buffer.data(nsg0 + 811);
    const auto *nsg0_812 = buffer.data(nsg0 + 812);
    const auto *nsg0_813 = buffer.data(nsg0 + 813);
    const auto *nsg0_814 = buffer.data(nsg0 + 814);
    const auto *nsg0_815 = buffer.data(nsg0 + 815);
    const auto *nsg0_825 = buffer.data(nsg0 + 825);
    const auto *nsg0_826 = buffer.data(nsg0 + 826);
    const auto *nsg0_828 = buffer.data(nsg0 + 828);
    const auto *nsg0_830 = buffer.data(nsg0 + 830);
    const auto *nsg0_831 = buffer.data(nsg0 + 831);
    const auto *nsg0_833 = buffer.data(nsg0 + 833);
    const auto *nsg0_834 = buffer.data(nsg0 + 834);
    const auto *nsg0_835 = buffer.data(nsg0 + 835);
    const auto *nsg0_836 = buffer.data(nsg0 + 836);
    const auto *nsg0_837 = buffer.data(nsg0 + 837);
    const auto *nsg0_838 = buffer.data(nsg0 + 838);
    const auto *nsg0_839 = buffer.data(nsg0 + 839);
    const auto *nsg0_842 = buffer.data(nsg0 + 842);
    const auto *nsg0_844 = buffer.data(nsg0 + 844);
    const auto *nsg0_845 = buffer.data(nsg0 + 845);
    const auto *nsg0_847 = buffer.data(nsg0 + 847);
    const auto *nsg0_848 = buffer.data(nsg0 + 848);
    const auto *nsg0_849 = buffer.data(nsg0 + 849);
    const auto *nsg0_851 = buffer.data(nsg0 + 851);
    const auto *nsg0_852 = buffer.data(nsg0 + 852);
    const auto *nsg0_853 = buffer.data(nsg0 + 853);
    const auto *nsg0_854 = buffer.data(nsg0 + 854);
    const auto *nsg0_855 = buffer.data(nsg0 + 855);
    const auto *nsg0_856 = buffer.data(nsg0 + 856);
    const auto *nsg0_857 = buffer.data(nsg0 + 857);
    const auto *nsg0_858 = buffer.data(nsg0 + 858);
    const auto *nsg0_859 = buffer.data(nsg0 + 859);
    const auto *nsg0_860 = buffer.data(nsg0 + 860);
    const auto *nsg0_861 = buffer.data(nsg0 + 861);
    const auto *nsg0_862 = buffer.data(nsg0 + 862);
    const auto *nsg0_863 = buffer.data(nsg0 + 863);
    const auto *nsg0_864 = buffer.data(nsg0 + 864);
    const auto *nsg0_865 = buffer.data(nsg0 + 865);
    const auto *nsg0_866 = buffer.data(nsg0 + 866);
    const auto *nsg0_867 = buffer.data(nsg0 + 867);
    const auto *nsg0_868 = buffer.data(nsg0 + 868);
    const auto *nsg0_869 = buffer.data(nsg0 + 869);
    const auto *nsg0_870 = buffer.data(nsg0 + 870);

    const auto *nsg1_810 = buffer.data(nsg1 + 810);
    const auto *nsg1_811 = buffer.data(nsg1 + 811);
    const auto *nsg1_812 = buffer.data(nsg1 + 812);
    const auto *nsg1_813 = buffer.data(nsg1 + 813);
    const auto *nsg1_814 = buffer.data(nsg1 + 814);
    const auto *nsg1_815 = buffer.data(nsg1 + 815);
    const auto *nsg1_825 = buffer.data(nsg1 + 825);
    const auto *nsg1_826 = buffer.data(nsg1 + 826);
    const auto *nsg1_828 = buffer.data(nsg1 + 828);
    const auto *nsg1_830 = buffer.data(nsg1 + 830);
    const auto *nsg1_831 = buffer.data(nsg1 + 831);
    const auto *nsg1_833 = buffer.data(nsg1 + 833);
    const auto *nsg1_834 = buffer.data(nsg1 + 834);
    const auto *nsg1_835 = buffer.data(nsg1 + 835);
    const auto *nsg1_836 = buffer.data(nsg1 + 836);
    const auto *nsg1_837 = buffer.data(nsg1 + 837);
    const auto *nsg1_838 = buffer.data(nsg1 + 838);
    const auto *nsg1_839 = buffer.data(nsg1 + 839);
    const auto *nsg1_842 = buffer.data(nsg1 + 842);
    const auto *nsg1_844 = buffer.data(nsg1 + 844);
    const auto *nsg1_845 = buffer.data(nsg1 + 845);
    const auto *nsg1_847 = buffer.data(nsg1 + 847);
    const auto *nsg1_848 = buffer.data(nsg1 + 848);
    const auto *nsg1_849 = buffer.data(nsg1 + 849);
    const auto *nsg1_851 = buffer.data(nsg1 + 851);
    const auto *nsg1_852 = buffer.data(nsg1 + 852);
    const auto *nsg1_853 = buffer.data(nsg1 + 853);
    const auto *nsg1_854 = buffer.data(nsg1 + 854);
    const auto *nsg1_855 = buffer.data(nsg1 + 855);
    const auto *nsg1_856 = buffer.data(nsg1 + 856);
    const auto *nsg1_857 = buffer.data(nsg1 + 857);
    const auto *nsg1_858 = buffer.data(nsg1 + 858);
    const auto *nsg1_859 = buffer.data(nsg1 + 859);
    const auto *nsg1_860 = buffer.data(nsg1 + 860);
    const auto *nsg1_861 = buffer.data(nsg1 + 861);
    const auto *nsg1_862 = buffer.data(nsg1 + 862);
    const auto *nsg1_863 = buffer.data(nsg1 + 863);
    const auto *nsg1_864 = buffer.data(nsg1 + 864);
    const auto *nsg1_865 = buffer.data(nsg1 + 865);
    const auto *nsg1_866 = buffer.data(nsg1 + 866);
    const auto *nsg1_867 = buffer.data(nsg1 + 867);
    const auto *nsg1_868 = buffer.data(nsg1 + 868);
    const auto *nsg1_869 = buffer.data(nsg1 + 869);
    const auto *nsg1_870 = buffer.data(nsg1 + 870);

    const auto *nsh_1128 = buffer.data(nsh + 1128);
    const auto *nsh_1131 = buffer.data(nsh + 1131);
    const auto *nsh_1132 = buffer.data(nsh + 1132);
    const auto *nsh_1133 = buffer.data(nsh + 1133);
    const auto *nsh_1134 = buffer.data(nsh + 1134);
    const auto *nsh_1135 = buffer.data(nsh + 1135);
    const auto *nsh_1136 = buffer.data(nsh + 1136);
    const auto *nsh_1137 = buffer.data(nsh + 1137);
    const auto *nsh_1138 = buffer.data(nsh + 1138);
    const auto *nsh_1139 = buffer.data(nsh + 1139);
    const auto *nsh_1140 = buffer.data(nsh + 1140);
    const auto *nsh_1141 = buffer.data(nsh + 1141);
    const auto *nsh_1142 = buffer.data(nsh + 1142);
    const auto *nsh_1143 = buffer.data(nsh + 1143);
    const auto *nsh_1148 = buffer.data(nsh + 1148);
    const auto *nsh_1149 = buffer.data(nsh + 1149);
    const auto *nsh_1150 = buffer.data(nsh + 1150);
    const auto *nsh_1151 = buffer.data(nsh + 1151);
    const auto *nsh_1152 = buffer.data(nsh + 1152);
    const auto *nsh_1154 = buffer.data(nsh + 1154);
    const auto *nsh_1155 = buffer.data(nsh + 1155);
    const auto *nsh_1156 = buffer.data(nsh + 1156);
    const auto *nsh_1158 = buffer.data(nsh + 1158);
    const auto *nsh_1160 = buffer.data(nsh + 1160);
    const auto *nsh_1161 = buffer.data(nsh + 1161);
    const auto *nsh_1163 = buffer.data(nsh + 1163);
    const auto *nsh_1164 = buffer.data(nsh + 1164);
    const auto *nsh_1165 = buffer.data(nsh + 1165);
    const auto *nsh_1167 = buffer.data(nsh + 1167);
    const auto *nsh_1168 = buffer.data(nsh + 1168);
    const auto *nsh_1169 = buffer.data(nsh + 1169);
    const auto *nsh_1170 = buffer.data(nsh + 1170);
    const auto *nsh_1171 = buffer.data(nsh + 1171);
    const auto *nsh_1172 = buffer.data(nsh + 1172);
    const auto *nsh_1173 = buffer.data(nsh + 1173);
    const auto *nsh_1174 = buffer.data(nsh + 1174);
    const auto *nsh_1175 = buffer.data(nsh + 1175);
    const auto *nsh_1178 = buffer.data(nsh + 1178);
    const auto *nsh_1180 = buffer.data(nsh + 1180);
    const auto *nsh_1181 = buffer.data(nsh + 1181);
    const auto *nsh_1183 = buffer.data(nsh + 1183);
    const auto *nsh_1184 = buffer.data(nsh + 1184);
    const auto *nsh_1185 = buffer.data(nsh + 1185);
    const auto *nsh_1187 = buffer.data(nsh + 1187);
    const auto *nsh_1188 = buffer.data(nsh + 1188);
    const auto *nsh_1189 = buffer.data(nsh + 1189);
    const auto *nsh_1190 = buffer.data(nsh + 1190);
    const auto *nsh_1191 = buffer.data(nsh + 1191);
    const auto *nsh_1192 = buffer.data(nsh + 1192);
    const auto *nsh_1193 = buffer.data(nsh + 1193);
    const auto *nsh_1194 = buffer.data(nsh + 1194);
    const auto *nsh_1195 = buffer.data(nsh + 1195);
    const auto *nsh_1196 = buffer.data(nsh + 1196);
    const auto *nsh_1197 = buffer.data(nsh + 1197);
    const auto *nsh_1198 = buffer.data(nsh + 1198);
    const auto *nsh_1199 = buffer.data(nsh + 1199);
    const auto *nsh_1200 = buffer.data(nsh + 1200);
    const auto *nsh_1201 = buffer.data(nsh + 1201);
    const auto *nsh_1202 = buffer.data(nsh + 1202);
    const auto *nsh_1203 = buffer.data(nsh + 1203);
    const auto *nsh_1204 = buffer.data(nsh + 1204);
    const auto *nsh_1205 = buffer.data(nsh + 1205);
    const auto *nsh_1206 = buffer.data(nsh + 1206);
    const auto *nsh_1207 = buffer.data(nsh + 1207);
    const auto *nsh_1208 = buffer.data(nsh + 1208);
    const auto *nsh_1209 = buffer.data(nsh + 1209);
    const auto *nsh_1210 = buffer.data(nsh + 1210);
    const auto *nsh_1211 = buffer.data(nsh + 1211);
    const auto *nsh_1212 = buffer.data(nsh + 1212);
    const auto *nsh_1213 = buffer.data(nsh + 1213);
    const auto *nsh_1214 = buffer.data(nsh + 1214);
    const auto *nsh_1215 = buffer.data(nsh + 1215);
    const auto *nsh_1216 = buffer.data(nsh + 1216);
    const auto *nsh_1217 = buffer.data(nsh + 1217);
    const auto *nsh_1218 = buffer.data(nsh + 1218);

#pragma omp simd aligned(t_1502, t_1503, t_1504, t_1505, pa_x, pc_x, msi0_1505, msh_1131, \
                         msh_1132, msh_1133, msi1_1505, nsh_1131, nsh_1132, \
                         nsh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = f_11 * msh_1131[k]
                    + f_3 * pc_x[k] * nsh_1131[k];

        t_1503[k] = f_11 * msh_1132[k]
                    + f_3 * pc_x[k] * nsh_1132[k];

        t_1504[k] = f_11 * msh_1133[k]
                    + f_3 * pc_x[k] * nsh_1133[k];

        t_1505[k] = pa_x[k] * msi0_1505[k]
                    - f_10 * pc_x[k] * msi1_1505[k];
    }

#pragma omp simd aligned(t_1506, t_1507, t_1508, t_1509, pa_x, pc_x, pc_z, msi0_1507, \
                         msi0_1508, msi0_1509, msh_918, msi1_1507, msi1_1508, msi1_1509, \
                         nsh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1506[k] = f_18 * msh_918[k]
                    + f_3 * pc_z[k] * nsh_1128[k];

        t_1507[k] = pa_x[k] * msi0_1507[k]
                    - f_10 * pc_x[k] * msi1_1507[k];

        t_1508[k] = pa_x[k] * msi0_1508[k]
                    - f_10 * pc_x[k] * msi1_1508[k];

        t_1509[k] = pa_x[k] * msi0_1509[k]
                    - f_10 * pc_x[k] * msi1_1509[k];
    }

#pragma omp simd aligned(t_1510, t_1511, t_1512, t_1513, pa_x, pc_x, pc_y, msi0_1511, \
                         msi0_1512, msh_944, msh_1134, msi1_1511, msi1_1512, nsh_1133, \
                         nsh_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1510[k] = f_11 * msh_944[k]
                    + f_3 * pc_y[k] * nsh_1133[k];

        t_1511[k] = pa_x[k] * msi0_1511[k]
                    - f_10 * pc_x[k] * msi1_1511[k];

        t_1512[k] = pa_x[k] * msi0_1512[k]
                    + f_20 * msh_1134[k]
                    - f_10 * pc_x[k] * msi1_1512[k];

        t_1513[k] = f_3 * pc_y[k] * nsh_1134[k];
    }

#pragma omp simd aligned(t_1514, t_1515, t_1516, pc_y, pc_z, msh_924, nsg0_810, nsg1_810, \
                         nsh_1134, nsh_1135, nsh_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1514[k] = f_15 * msh_924[k]
                    + f_3 * pc_z[k] * nsh_1134[k];

        t_1515[k] = f_4 * nsg0_810[k]
                    - f_5 * nsg1_810[k]
                    + f_3 * pc_y[k] * nsh_1135[k];

        t_1516[k] = f_3 * pc_y[k] * nsh_1136[k];
    }

#pragma omp simd aligned(t_1517, t_1518, t_1519, pa_x, pc_x, pc_y, msi0_1517, msh_1139, \
                         msi1_1517, nsg0_811, nsg0_812, nsg1_811, nsg1_812, nsh_1137, \
                         nsh_1138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1517[k] = pa_x[k] * msi0_1517[k]
                    + f_14 * msh_1139[k]
                    - f_10 * pc_x[k] * msi1_1517[k];

        t_1518[k] = f_6 * nsg0_811[k]
                    - f_7 * nsg1_811[k]
                    + f_3 * pc_y[k] * nsh_1137[k];

        t_1519[k] = f_4 * nsg0_812[k]
                    - f_5 * nsg1_812[k]
                    + f_3 * pc_y[k] * nsh_1138[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, pa_x, pc_x, pc_y, msi0_1521, msh_1143, \
                         msi1_1521, nsg0_813, nsg1_813, nsh_1139, \
                         nsh_1140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_3 * pc_y[k] * nsh_1139[k];

        t_1521[k] = pa_x[k] * msi0_1521[k]
                    + f_13 * msh_1143[k]
                    - f_10 * pc_x[k] * msi1_1521[k];

        t_1522[k] = f_8 * nsg0_813[k]
                    - f_9 * nsg1_813[k]
                    + f_3 * pc_y[k] * nsh_1140[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pc_y, nsg0_814, nsg0_815, nsg1_814, nsg1_815, \
                         nsh_1141, nsh_1142, nsh_1143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_6 * nsg0_814[k]
                    - f_7 * nsg1_814[k]
                    + f_3 * pc_y[k] * nsh_1141[k];

        t_1524[k] = f_4 * nsg0_815[k]
                    - f_5 * nsg1_815[k]
                    + f_3 * pc_y[k] * nsh_1142[k];

        t_1525[k] = f_3 * pc_y[k] * nsh_1143[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, t_1529, pa_x, pc_x, msi0_1526, msh_1148, \
                         msh_1149, msh_1150, msh_1151, msi1_1526, nsh_1149, nsh_1150, \
                         nsh_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = pa_x[k] * msi0_1526[k]
                    + f_12 * msh_1148[k]
                    - f_10 * pc_x[k] * msi1_1526[k];

        t_1527[k] = f_11 * msh_1149[k]
                    + f_3 * pc_x[k] * nsh_1149[k];

        t_1528[k] = f_11 * msh_1150[k]
                    + f_3 * pc_x[k] * nsh_1150[k];

        t_1529[k] = f_11 * msh_1151[k]
                    + f_3 * pc_x[k] * nsh_1151[k];
    }

#pragma omp simd aligned(t_1530, t_1531, t_1532, t_1533, pa_x, pc_x, pc_y, msi0_1533, \
                         msh_1152, msh_1154, msi1_1533, nsh_1148, nsh_1152, \
                         nsh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1530[k] = f_11 * msh_1152[k]
                    + f_3 * pc_x[k] * nsh_1152[k];

        t_1531[k] = f_3 * pc_y[k] * nsh_1148[k];

        t_1532[k] = f_11 * msh_1154[k]
                    + f_3 * pc_x[k] * nsh_1154[k];

        t_1533[k] = pa_x[k] * msi0_1533[k]
                    - f_10 * pc_x[k] * msi1_1533[k];
    }

#pragma omp simd aligned(t_1534, t_1535, t_1536, t_1537, pa_x, pc_x, msi0_1534, msi0_1535, \
                         msi0_1536, msi0_1537, msi1_1534, msi1_1535, msi1_1536, \
                         msi1_1537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1534[k] = pa_x[k] * msi0_1534[k]
                    - f_10 * pc_x[k] * msi1_1534[k];

        t_1535[k] = pa_x[k] * msi0_1535[k]
                    - f_10 * pc_x[k] * msi1_1535[k];

        t_1536[k] = pa_x[k] * msi0_1536[k]
                    - f_10 * pc_x[k] * msi1_1536[k];

        t_1537[k] = pa_x[k] * msi0_1537[k]
                    - f_10 * pc_x[k] * msi1_1537[k];
    }

#pragma omp simd aligned(t_1538, t_1539, t_1540, t_1541, pa_x, pc_x, pc_y, msi0_1539, \
                         msi1_1539, nsg0_825, nsg0_826, nsg1_825, nsg1_826, nsh_1154, \
                         nsh_1155, nsh_1156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1538[k] = f_3 * pc_y[k] * nsh_1154[k];

        t_1539[k] = pa_x[k] * msi0_1539[k]
                    - f_10 * pc_x[k] * msi1_1539[k];

        t_1540[k] = f_1 * nsg0_825[k]
                    - f_2 * nsg1_825[k]
                    + f_3 * pc_x[k] * nsh_1155[k];

        t_1541[k] = f_16 * nsg0_826[k]
                    - f_17 * nsg1_826[k]
                    + f_3 * pc_x[k] * nsh_1156[k];
    }

#pragma omp simd aligned(t_1542, t_1543, t_1544, t_1545, pc_x, pc_z, nsg0_828, nsg0_830, \
                         nsg1_828, nsg1_830, nsh_1155, nsh_1156, nsh_1158, \
                         nsh_1160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1542[k] = f_3 * pc_z[k] * nsh_1155[k];

        t_1543[k] = f_8 * nsg0_828[k]
                    - f_9 * nsg1_828[k]
                    + f_3 * pc_x[k] * nsh_1158[k];

        t_1544[k] = f_3 * pc_z[k] * nsh_1156[k];

        t_1545[k] = f_8 * nsg0_830[k]
                    - f_9 * nsg1_830[k]
                    + f_3 * pc_x[k] * nsh_1160[k];
    }

#pragma omp simd aligned(t_1546, t_1547, t_1548, t_1549, pc_x, pc_z, nsg0_831, nsg0_833, \
                         nsg0_834, nsg1_831, nsg1_833, nsg1_834, nsh_1158, nsh_1161, nsh_1163, \
                         nsh_1164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1546[k] = f_6 * nsg0_831[k]
                    - f_7 * nsg1_831[k]
                    + f_3 * pc_x[k] * nsh_1161[k];

        t_1547[k] = f_3 * pc_z[k] * nsh_1158[k];

        t_1548[k] = f_6 * nsg0_833[k]
                    - f_7 * nsg1_833[k]
                    + f_3 * pc_x[k] * nsh_1163[k];

        t_1549[k] = f_6 * nsg0_834[k]
                    - f_7 * nsg1_834[k]
                    + f_3 * pc_x[k] * nsh_1164[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, t_1553, pc_x, pc_z, nsg0_835, nsg0_837, \
                         nsg0_838, nsg1_835, nsg1_837, nsg1_838, nsh_1161, nsh_1165, nsh_1167, \
                         nsh_1168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = f_4 * nsg0_835[k]
                    - f_5 * nsg1_835[k]
                    + f_3 * pc_x[k] * nsh_1165[k];

        t_1551[k] = f_3 * pc_z[k] * nsh_1161[k];

        t_1552[k] = f_4 * nsg0_837[k]
                    - f_5 * nsg1_837[k]
                    + f_3 * pc_x[k] * nsh_1167[k];

        t_1553[k] = f_4 * nsg0_838[k]
                    - f_5 * nsg1_838[k]
                    + f_3 * pc_x[k] * nsh_1168[k];
    }

#pragma omp simd aligned(t_1554, t_1555, t_1556, t_1557, t_1558, t_1559, pc_x, nsg0_839, \
                         nsg1_839, nsh_1169, nsh_1170, nsh_1171, nsh_1172, nsh_1173, \
                         nsh_1174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1554[k] = f_4 * nsg0_839[k]
                    - f_5 * nsg1_839[k]
                    + f_3 * pc_x[k] * nsh_1169[k];

        t_1555[k] = f_3 * pc_x[k] * nsh_1170[k];

        t_1556[k] = f_3 * pc_x[k] * nsh_1171[k];

        t_1557[k] = f_3 * pc_x[k] * nsh_1172[k];

        t_1558[k] = f_3 * pc_x[k] * nsh_1173[k];

        t_1559[k] = f_3 * pc_x[k] * nsh_1174[k];
    }

#pragma omp simd aligned(t_1560, t_1561, t_1562, t_1563, pc_x, pc_y, pc_z, msh_960, nsg0_835, \
                         nsg1_835, nsh_1170, nsh_1171, nsh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1560[k] = f_3 * pc_x[k] * nsh_1175[k];

        t_1561[k] = f_0 * msh_960[k]
                    + f_1 * nsg0_835[k]
                    - f_2 * nsg1_835[k]
                    + f_3 * pc_y[k] * nsh_1170[k];

        t_1562[k] = f_3 * pc_z[k] * nsh_1170[k];

        t_1563[k] = f_4 * nsg0_835[k]
                    - f_5 * nsg1_835[k]
                    + f_3 * pc_z[k] * nsh_1171[k];
    }

#pragma omp simd aligned(t_1564, t_1565, t_1566, t_1567, pc_y, pc_z, msh_965, nsg0_836, \
                         nsg0_837, nsg0_839, nsg1_836, nsg1_837, nsg1_839, nsh_1172, nsh_1173, \
                         nsh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1564[k] = f_6 * nsg0_836[k]
                    - f_7 * nsg1_836[k]
                    + f_3 * pc_z[k] * nsh_1172[k];

        t_1565[k] = f_8 * nsg0_837[k]
                    - f_9 * nsg1_837[k]
                    + f_3 * pc_z[k] * nsh_1173[k];

        t_1566[k] = f_0 * msh_965[k]
                    + f_3 * pc_y[k] * nsh_1175[k];

        t_1567[k] = f_1 * nsg0_839[k]
                    - f_2 * nsg1_839[k]
                    + f_3 * pc_z[k] * nsh_1175[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, t_1571, pa_z, pc_x, pc_z, msi0_1260, \
                         msi0_1261, msi0_1263, msi1_1260, msi1_1261, msi1_1263, nsg0_842, \
                         nsg1_842, nsh_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = pa_z[k] * msi0_1260[k]
                    - f_10 * pc_z[k] * msi1_1260[k];

        t_1569[k] = pa_z[k] * msi0_1261[k]
                    - f_10 * pc_z[k] * msi1_1261[k];

        t_1570[k] = f_16 * nsg0_842[k]
                    - f_17 * nsg1_842[k]
                    + f_3 * pc_x[k] * nsh_1178[k];

        t_1571[k] = pa_z[k] * msi0_1263[k]
                    - f_10 * pc_z[k] * msi1_1263[k];
    }

#pragma omp simd aligned(t_1572, t_1573, t_1574, pa_z, pc_x, pc_z, msi0_1266, msi1_1266, \
                         nsg0_844, nsg0_845, nsg1_844, nsg1_845, nsh_1180, \
                         nsh_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1572[k] = f_8 * nsg0_844[k]
                    - f_9 * nsg1_844[k]
                    + f_3 * pc_x[k] * nsh_1180[k];

        t_1573[k] = f_8 * nsg0_845[k]
                    - f_9 * nsg1_845[k]
                    + f_3 * pc_x[k] * nsh_1181[k];

        t_1574[k] = pa_z[k] * msi0_1266[k]
                    - f_10 * pc_z[k] * msi1_1266[k];
    }

#pragma omp simd aligned(t_1575, t_1576, t_1577, pc_x, nsg0_847, nsg0_848, nsg0_849, nsg1_847, \
                         nsg1_848, nsg1_849, nsh_1183, nsh_1184, \
                         nsh_1185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1575[k] = f_6 * nsg0_847[k]
                    - f_7 * nsg1_847[k]
                    + f_3 * pc_x[k] * nsh_1183[k];

        t_1576[k] = f_6 * nsg0_848[k]
                    - f_7 * nsg1_848[k]
                    + f_3 * pc_x[k] * nsh_1184[k];

        t_1577[k] = f_6 * nsg0_849[k]
                    - f_7 * nsg1_849[k]
                    + f_3 * pc_x[k] * nsh_1185[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pa_z, pc_x, pc_z, msi0_1270, msi1_1270, \
                         nsg0_851, nsg0_852, nsg1_851, nsg1_852, nsh_1187, \
                         nsh_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = pa_z[k] * msi0_1270[k]
                    - f_10 * pc_z[k] * msi1_1270[k];

        t_1579[k] = f_4 * nsg0_851[k]
                    - f_5 * nsg1_851[k]
                    + f_3 * pc_x[k] * nsh_1187[k];

        t_1580[k] = f_4 * nsg0_852[k]
                    - f_5 * nsg1_852[k]
                    + f_3 * pc_x[k] * nsh_1188[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, t_1584, t_1585, pc_x, nsg0_853, nsg0_854, \
                         nsg1_853, nsg1_854, nsh_1189, nsh_1190, nsh_1191, nsh_1192, \
                         nsh_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_4 * nsg0_853[k]
                    - f_5 * nsg1_853[k]
                    + f_3 * pc_x[k] * nsh_1189[k];

        t_1582[k] = f_4 * nsg0_854[k]
                    - f_5 * nsg1_854[k]
                    + f_3 * pc_x[k] * nsh_1190[k];

        t_1583[k] = f_3 * pc_x[k] * nsh_1191[k];

        t_1584[k] = f_3 * pc_x[k] * nsh_1192[k];

        t_1585[k] = f_3 * pc_x[k] * nsh_1193[k];
    }

#pragma omp simd aligned(t_1586, t_1587, t_1588, t_1589, t_1590, pa_z, pc_x, pc_z, msi0_1281, \
                         msh_960, msi1_1281, nsh_1191, nsh_1194, nsh_1195, \
                         nsh_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1586[k] = f_3 * pc_x[k] * nsh_1194[k];

        t_1587[k] = f_3 * pc_x[k] * nsh_1195[k];

        t_1588[k] = f_3 * pc_x[k] * nsh_1196[k];

        t_1589[k] = pa_z[k] * msi0_1281[k]
                    - f_10 * pc_z[k] * msi1_1281[k];

        t_1590[k] = f_11 * msh_960[k]
                    + f_3 * pc_z[k] * nsh_1191[k];
    }

#pragma omp simd aligned(t_1591, t_1592, t_1593, pa_z, pc_z, msi0_1283, msi0_1284, msi0_1285, \
                         msh_961, msh_962, msh_963, msi1_1283, msi1_1284, \
                         msi1_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1591[k] = pa_z[k] * msi0_1283[k]
                    + f_12 * msh_961[k]
                    - f_10 * pc_z[k] * msi1_1283[k];

        t_1592[k] = pa_z[k] * msi0_1284[k]
                    + f_13 * msh_962[k]
                    - f_10 * pc_z[k] * msi1_1284[k];

        t_1593[k] = pa_z[k] * msi0_1285[k]
                    + f_14 * msh_963[k]
                    - f_10 * pc_z[k] * msi1_1285[k];
    }

#pragma omp simd aligned(t_1594, t_1595, t_1596, pc_x, pc_y, pc_z, msh_965, msh_986, nsg0_854, \
                         nsg0_855, nsg1_854, nsg1_855, nsh_1196, \
                         nsh_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1594[k] = f_15 * msh_986[k]
                    + f_3 * pc_y[k] * nsh_1196[k];

        t_1595[k] = f_11 * msh_965[k]
                    + f_1 * nsg0_854[k]
                    - f_2 * nsg1_854[k]
                    + f_3 * pc_z[k] * nsh_1196[k];

        t_1596[k] = f_1 * nsg0_855[k]
                    - f_2 * nsg1_855[k]
                    + f_3 * pc_x[k] * nsh_1197[k];
    }

#pragma omp simd aligned(t_1597, t_1598, t_1599, pc_x, nsg0_856, nsg0_857, nsg0_858, nsg1_856, \
                         nsg1_857, nsg1_858, nsh_1198, nsh_1199, \
                         nsh_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1597[k] = f_16 * nsg0_856[k]
                    - f_17 * nsg1_856[k]
                    + f_3 * pc_x[k] * nsh_1198[k];

        t_1598[k] = f_16 * nsg0_857[k]
                    - f_17 * nsg1_857[k]
                    + f_3 * pc_x[k] * nsh_1199[k];

        t_1599[k] = f_8 * nsg0_858[k]
                    - f_9 * nsg1_858[k]
                    + f_3 * pc_x[k] * nsh_1200[k];
    }

#pragma omp simd aligned(t_1600, t_1601, t_1602, pc_x, nsg0_859, nsg0_860, nsg0_861, nsg1_859, \
                         nsg1_860, nsg1_861, nsh_1201, nsh_1202, \
                         nsh_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1600[k] = f_8 * nsg0_859[k]
                    - f_9 * nsg1_859[k]
                    + f_3 * pc_x[k] * nsh_1201[k];

        t_1601[k] = f_8 * nsg0_860[k]
                    - f_9 * nsg1_860[k]
                    + f_3 * pc_x[k] * nsh_1202[k];

        t_1602[k] = f_6 * nsg0_861[k]
                    - f_7 * nsg1_861[k]
                    + f_3 * pc_x[k] * nsh_1203[k];
    }

#pragma omp simd aligned(t_1603, t_1604, t_1605, pc_x, nsg0_862, nsg0_863, nsg0_864, nsg1_862, \
                         nsg1_863, nsg1_864, nsh_1204, nsh_1205, \
                         nsh_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1603[k] = f_6 * nsg0_862[k]
                    - f_7 * nsg1_862[k]
                    + f_3 * pc_x[k] * nsh_1204[k];

        t_1604[k] = f_6 * nsg0_863[k]
                    - f_7 * nsg1_863[k]
                    + f_3 * pc_x[k] * nsh_1205[k];

        t_1605[k] = f_6 * nsg0_864[k]
                    - f_7 * nsg1_864[k]
                    + f_3 * pc_x[k] * nsh_1206[k];
    }

#pragma omp simd aligned(t_1606, t_1607, t_1608, pc_x, nsg0_865, nsg0_866, nsg0_867, nsg1_865, \
                         nsg1_866, nsg1_867, nsh_1207, nsh_1208, \
                         nsh_1209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1606[k] = f_4 * nsg0_865[k]
                    - f_5 * nsg1_865[k]
                    + f_3 * pc_x[k] * nsh_1207[k];

        t_1607[k] = f_4 * nsg0_866[k]
                    - f_5 * nsg1_866[k]
                    + f_3 * pc_x[k] * nsh_1208[k];

        t_1608[k] = f_4 * nsg0_867[k]
                    - f_5 * nsg1_867[k]
                    + f_3 * pc_x[k] * nsh_1209[k];
    }

#pragma omp simd aligned(t_1609, t_1610, t_1611, t_1612, t_1613, pc_x, nsg0_868, nsg0_869, \
                         nsg1_868, nsg1_869, nsh_1210, nsh_1211, nsh_1212, nsh_1213, \
                         nsh_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1609[k] = f_4 * nsg0_868[k]
                    - f_5 * nsg1_868[k]
                    + f_3 * pc_x[k] * nsh_1210[k];

        t_1610[k] = f_4 * nsg0_869[k]
                    - f_5 * nsg1_869[k]
                    + f_3 * pc_x[k] * nsh_1211[k];

        t_1611[k] = f_3 * pc_x[k] * nsh_1212[k];

        t_1612[k] = f_3 * pc_x[k] * nsh_1213[k];

        t_1613[k] = f_3 * pc_x[k] * nsh_1214[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, t_1617, t_1618, pc_x, pc_y, pc_z, msh_981, \
                         msh_1002, nsg0_865, nsg1_865, nsh_1212, nsh_1215, nsh_1216, \
                         nsh_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = f_3 * pc_x[k] * nsh_1215[k];

        t_1615[k] = f_3 * pc_x[k] * nsh_1216[k];

        t_1616[k] = f_3 * pc_x[k] * nsh_1217[k];

        t_1617[k] = f_18 * msh_1002[k]
                    + f_1 * nsg0_865[k]
                    - f_2 * nsg1_865[k]
                    + f_3 * pc_y[k] * nsh_1212[k];

        t_1618[k] = f_12 * msh_981[k]
                    + f_3 * pc_z[k] * nsh_1212[k];
    }

#pragma omp simd aligned(t_1619, t_1620, t_1621, pc_y, msh_1004, msh_1005, msh_1006, nsg0_867, \
                         nsg0_868, nsg0_869, nsg1_867, nsg1_868, nsg1_869, nsh_1214, nsh_1215, \
                         nsh_1216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1619[k] = f_18 * msh_1004[k]
                    + f_8 * nsg0_867[k]
                    - f_9 * nsg1_867[k]
                    + f_3 * pc_y[k] * nsh_1214[k];

        t_1620[k] = f_18 * msh_1005[k]
                    + f_6 * nsg0_868[k]
                    - f_7 * nsg1_868[k]
                    + f_3 * pc_y[k] * nsh_1215[k];

        t_1621[k] = f_18 * msh_1006[k]
                    + f_4 * nsg0_869[k]
                    - f_5 * nsg1_869[k]
                    + f_3 * pc_y[k] * nsh_1216[k];
    }

#pragma omp simd aligned(t_1622, t_1623, t_1624, pc_x, pc_y, pc_z, msh_986, msh_1007, \
                         nsg0_869, nsg0_870, nsg1_869, nsg1_870, nsh_1217, \
                         nsh_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1622[k] = f_18 * msh_1007[k]
                    + f_3 * pc_y[k] * nsh_1217[k];

        t_1623[k] = f_12 * msh_986[k]
                    + f_1 * nsg0_869[k]
                    - f_2 * nsg1_869[k]
                    + f_3 * pc_z[k] * nsh_1217[k];

        t_1624[k] = f_1 * nsg0_870[k]
                    - f_2 * nsg1_870[k]
                    + f_3 * pc_x[k] * nsh_1218[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t msh, const size_t nsg0,
                                                           const size_t nsg1, const size_t nsh,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 2.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msh_1002 = buffer.data(msh + 1002);
    const auto *msh_1007 = buffer.data(msh + 1007);
    const auto *msh_1023 = buffer.data(msh + 1023);
    const auto *msh_1025 = buffer.data(msh + 1025);
    const auto *msh_1026 = buffer.data(msh + 1026);
    const auto *msh_1027 = buffer.data(msh + 1027);
    const auto *msh_1028 = buffer.data(msh + 1028);
    const auto *msh_1044 = buffer.data(msh + 1044);
    const auto *msh_1046 = buffer.data(msh + 1046);
    const auto *msh_1047 = buffer.data(msh + 1047);
    const auto *msh_1048 = buffer.data(msh + 1048);
    const auto *msh_1049 = buffer.data(msh + 1049);
    const auto *msh_1065 = buffer.data(msh + 1065);
    const auto *msh_1067 = buffer.data(msh + 1067);
    const auto *msh_1068 = buffer.data(msh + 1068);
    const auto *msh_1069 = buffer.data(msh + 1069);
    const auto *msh_1070 = buffer.data(msh + 1070);
    const auto *msh_1086 = buffer.data(msh + 1086);
    const auto *msh_1088 = buffer.data(msh + 1088);
    const auto *msh_1089 = buffer.data(msh + 1089);
    const auto *msh_1090 = buffer.data(msh + 1090);
    const auto *msh_1091 = buffer.data(msh + 1091);

    const auto *nsg0_871 = buffer.data(nsg0 + 871);
    const auto *nsg0_872 = buffer.data(nsg0 + 872);
    const auto *nsg0_873 = buffer.data(nsg0 + 873);
    const auto *nsg0_874 = buffer.data(nsg0 + 874);
    const auto *nsg0_875 = buffer.data(nsg0 + 875);
    const auto *nsg0_876 = buffer.data(nsg0 + 876);
    const auto *nsg0_877 = buffer.data(nsg0 + 877);
    const auto *nsg0_878 = buffer.data(nsg0 + 878);
    const auto *nsg0_879 = buffer.data(nsg0 + 879);
    const auto *nsg0_880 = buffer.data(nsg0 + 880);
    const auto *nsg0_881 = buffer.data(nsg0 + 881);
    const auto *nsg0_882 = buffer.data(nsg0 + 882);
    const auto *nsg0_883 = buffer.data(nsg0 + 883);
    const auto *nsg0_884 = buffer.data(nsg0 + 884);
    const auto *nsg0_885 = buffer.data(nsg0 + 885);
    const auto *nsg0_886 = buffer.data(nsg0 + 886);
    const auto *nsg0_887 = buffer.data(nsg0 + 887);
    const auto *nsg0_888 = buffer.data(nsg0 + 888);
    const auto *nsg0_889 = buffer.data(nsg0 + 889);
    const auto *nsg0_890 = buffer.data(nsg0 + 890);
    const auto *nsg0_891 = buffer.data(nsg0 + 891);
    const auto *nsg0_892 = buffer.data(nsg0 + 892);
    const auto *nsg0_893 = buffer.data(nsg0 + 893);
    const auto *nsg0_894 = buffer.data(nsg0 + 894);
    const auto *nsg0_895 = buffer.data(nsg0 + 895);
    const auto *nsg0_896 = buffer.data(nsg0 + 896);
    const auto *nsg0_897 = buffer.data(nsg0 + 897);
    const auto *nsg0_898 = buffer.data(nsg0 + 898);
    const auto *nsg0_899 = buffer.data(nsg0 + 899);
    const auto *nsg0_900 = buffer.data(nsg0 + 900);
    const auto *nsg0_901 = buffer.data(nsg0 + 901);
    const auto *nsg0_902 = buffer.data(nsg0 + 902);
    const auto *nsg0_903 = buffer.data(nsg0 + 903);
    const auto *nsg0_904 = buffer.data(nsg0 + 904);
    const auto *nsg0_905 = buffer.data(nsg0 + 905);
    const auto *nsg0_906 = buffer.data(nsg0 + 906);
    const auto *nsg0_907 = buffer.data(nsg0 + 907);
    const auto *nsg0_908 = buffer.data(nsg0 + 908);
    const auto *nsg0_909 = buffer.data(nsg0 + 909);
    const auto *nsg0_910 = buffer.data(nsg0 + 910);
    const auto *nsg0_911 = buffer.data(nsg0 + 911);
    const auto *nsg0_912 = buffer.data(nsg0 + 912);
    const auto *nsg0_913 = buffer.data(nsg0 + 913);
    const auto *nsg0_914 = buffer.data(nsg0 + 914);
    const auto *nsg0_915 = buffer.data(nsg0 + 915);
    const auto *nsg0_916 = buffer.data(nsg0 + 916);
    const auto *nsg0_917 = buffer.data(nsg0 + 917);
    const auto *nsg0_918 = buffer.data(nsg0 + 918);
    const auto *nsg0_919 = buffer.data(nsg0 + 919);
    const auto *nsg0_920 = buffer.data(nsg0 + 920);
    const auto *nsg0_921 = buffer.data(nsg0 + 921);
    const auto *nsg0_922 = buffer.data(nsg0 + 922);
    const auto *nsg0_923 = buffer.data(nsg0 + 923);
    const auto *nsg0_924 = buffer.data(nsg0 + 924);
    const auto *nsg0_925 = buffer.data(nsg0 + 925);
    const auto *nsg0_926 = buffer.data(nsg0 + 926);
    const auto *nsg0_927 = buffer.data(nsg0 + 927);
    const auto *nsg0_928 = buffer.data(nsg0 + 928);
    const auto *nsg0_929 = buffer.data(nsg0 + 929);
    const auto *nsg0_930 = buffer.data(nsg0 + 930);
    const auto *nsg0_931 = buffer.data(nsg0 + 931);
    const auto *nsg0_932 = buffer.data(nsg0 + 932);
    const auto *nsg0_933 = buffer.data(nsg0 + 933);

    const auto *nsg1_871 = buffer.data(nsg1 + 871);
    const auto *nsg1_872 = buffer.data(nsg1 + 872);
    const auto *nsg1_873 = buffer.data(nsg1 + 873);
    const auto *nsg1_874 = buffer.data(nsg1 + 874);
    const auto *nsg1_875 = buffer.data(nsg1 + 875);
    const auto *nsg1_876 = buffer.data(nsg1 + 876);
    const auto *nsg1_877 = buffer.data(nsg1 + 877);
    const auto *nsg1_878 = buffer.data(nsg1 + 878);
    const auto *nsg1_879 = buffer.data(nsg1 + 879);
    const auto *nsg1_880 = buffer.data(nsg1 + 880);
    const auto *nsg1_881 = buffer.data(nsg1 + 881);
    const auto *nsg1_882 = buffer.data(nsg1 + 882);
    const auto *nsg1_883 = buffer.data(nsg1 + 883);
    const auto *nsg1_884 = buffer.data(nsg1 + 884);
    const auto *nsg1_885 = buffer.data(nsg1 + 885);
    const auto *nsg1_886 = buffer.data(nsg1 + 886);
    const auto *nsg1_887 = buffer.data(nsg1 + 887);
    const auto *nsg1_888 = buffer.data(nsg1 + 888);
    const auto *nsg1_889 = buffer.data(nsg1 + 889);
    const auto *nsg1_890 = buffer.data(nsg1 + 890);
    const auto *nsg1_891 = buffer.data(nsg1 + 891);
    const auto *nsg1_892 = buffer.data(nsg1 + 892);
    const auto *nsg1_893 = buffer.data(nsg1 + 893);
    const auto *nsg1_894 = buffer.data(nsg1 + 894);
    const auto *nsg1_895 = buffer.data(nsg1 + 895);
    const auto *nsg1_896 = buffer.data(nsg1 + 896);
    const auto *nsg1_897 = buffer.data(nsg1 + 897);
    const auto *nsg1_898 = buffer.data(nsg1 + 898);
    const auto *nsg1_899 = buffer.data(nsg1 + 899);
    const auto *nsg1_900 = buffer.data(nsg1 + 900);
    const auto *nsg1_901 = buffer.data(nsg1 + 901);
    const auto *nsg1_902 = buffer.data(nsg1 + 902);
    const auto *nsg1_903 = buffer.data(nsg1 + 903);
    const auto *nsg1_904 = buffer.data(nsg1 + 904);
    const auto *nsg1_905 = buffer.data(nsg1 + 905);
    const auto *nsg1_906 = buffer.data(nsg1 + 906);
    const auto *nsg1_907 = buffer.data(nsg1 + 907);
    const auto *nsg1_908 = buffer.data(nsg1 + 908);
    const auto *nsg1_909 = buffer.data(nsg1 + 909);
    const auto *nsg1_910 = buffer.data(nsg1 + 910);
    const auto *nsg1_911 = buffer.data(nsg1 + 911);
    const auto *nsg1_912 = buffer.data(nsg1 + 912);
    const auto *nsg1_913 = buffer.data(nsg1 + 913);
    const auto *nsg1_914 = buffer.data(nsg1 + 914);
    const auto *nsg1_915 = buffer.data(nsg1 + 915);
    const auto *nsg1_916 = buffer.data(nsg1 + 916);
    const auto *nsg1_917 = buffer.data(nsg1 + 917);
    const auto *nsg1_918 = buffer.data(nsg1 + 918);
    const auto *nsg1_919 = buffer.data(nsg1 + 919);
    const auto *nsg1_920 = buffer.data(nsg1 + 920);
    const auto *nsg1_921 = buffer.data(nsg1 + 921);
    const auto *nsg1_922 = buffer.data(nsg1 + 922);
    const auto *nsg1_923 = buffer.data(nsg1 + 923);
    const auto *nsg1_924 = buffer.data(nsg1 + 924);
    const auto *nsg1_925 = buffer.data(nsg1 + 925);
    const auto *nsg1_926 = buffer.data(nsg1 + 926);
    const auto *nsg1_927 = buffer.data(nsg1 + 927);
    const auto *nsg1_928 = buffer.data(nsg1 + 928);
    const auto *nsg1_929 = buffer.data(nsg1 + 929);
    const auto *nsg1_930 = buffer.data(nsg1 + 930);
    const auto *nsg1_931 = buffer.data(nsg1 + 931);
    const auto *nsg1_932 = buffer.data(nsg1 + 932);
    const auto *nsg1_933 = buffer.data(nsg1 + 933);

    const auto *nsh_1219 = buffer.data(nsh + 1219);
    const auto *nsh_1220 = buffer.data(nsh + 1220);
    const auto *nsh_1221 = buffer.data(nsh + 1221);
    const auto *nsh_1222 = buffer.data(nsh + 1222);
    const auto *nsh_1223 = buffer.data(nsh + 1223);
    const auto *nsh_1224 = buffer.data(nsh + 1224);
    const auto *nsh_1225 = buffer.data(nsh + 1225);
    const auto *nsh_1226 = buffer.data(nsh + 1226);
    const auto *nsh_1227 = buffer.data(nsh + 1227);
    const auto *nsh_1228 = buffer.data(nsh + 1228);
    const auto *nsh_1229 = buffer.data(nsh + 1229);
    const auto *nsh_1230 = buffer.data(nsh + 1230);
    const auto *nsh_1231 = buffer.data(nsh + 1231);
    const auto *nsh_1232 = buffer.data(nsh + 1232);
    const auto *nsh_1233 = buffer.data(nsh + 1233);
    const auto *nsh_1234 = buffer.data(nsh + 1234);
    const auto *nsh_1235 = buffer.data(nsh + 1235);
    const auto *nsh_1236 = buffer.data(nsh + 1236);
    const auto *nsh_1237 = buffer.data(nsh + 1237);
    const auto *nsh_1238 = buffer.data(nsh + 1238);
    const auto *nsh_1239 = buffer.data(nsh + 1239);
    const auto *nsh_1240 = buffer.data(nsh + 1240);
    const auto *nsh_1241 = buffer.data(nsh + 1241);
    const auto *nsh_1242 = buffer.data(nsh + 1242);
    const auto *nsh_1243 = buffer.data(nsh + 1243);
    const auto *nsh_1244 = buffer.data(nsh + 1244);
    const auto *nsh_1245 = buffer.data(nsh + 1245);
    const auto *nsh_1246 = buffer.data(nsh + 1246);
    const auto *nsh_1247 = buffer.data(nsh + 1247);
    const auto *nsh_1248 = buffer.data(nsh + 1248);
    const auto *nsh_1249 = buffer.data(nsh + 1249);
    const auto *nsh_1250 = buffer.data(nsh + 1250);
    const auto *nsh_1251 = buffer.data(nsh + 1251);
    const auto *nsh_1252 = buffer.data(nsh + 1252);
    const auto *nsh_1253 = buffer.data(nsh + 1253);
    const auto *nsh_1254 = buffer.data(nsh + 1254);
    const auto *nsh_1255 = buffer.data(nsh + 1255);
    const auto *nsh_1256 = buffer.data(nsh + 1256);
    const auto *nsh_1257 = buffer.data(nsh + 1257);
    const auto *nsh_1258 = buffer.data(nsh + 1258);
    const auto *nsh_1259 = buffer.data(nsh + 1259);
    const auto *nsh_1260 = buffer.data(nsh + 1260);
    const auto *nsh_1261 = buffer.data(nsh + 1261);
    const auto *nsh_1262 = buffer.data(nsh + 1262);
    const auto *nsh_1263 = buffer.data(nsh + 1263);
    const auto *nsh_1264 = buffer.data(nsh + 1264);
    const auto *nsh_1265 = buffer.data(nsh + 1265);
    const auto *nsh_1266 = buffer.data(nsh + 1266);
    const auto *nsh_1267 = buffer.data(nsh + 1267);
    const auto *nsh_1268 = buffer.data(nsh + 1268);
    const auto *nsh_1269 = buffer.data(nsh + 1269);
    const auto *nsh_1270 = buffer.data(nsh + 1270);
    const auto *nsh_1271 = buffer.data(nsh + 1271);
    const auto *nsh_1272 = buffer.data(nsh + 1272);
    const auto *nsh_1273 = buffer.data(nsh + 1273);
    const auto *nsh_1274 = buffer.data(nsh + 1274);
    const auto *nsh_1275 = buffer.data(nsh + 1275);
    const auto *nsh_1276 = buffer.data(nsh + 1276);
    const auto *nsh_1277 = buffer.data(nsh + 1277);
    const auto *nsh_1278 = buffer.data(nsh + 1278);
    const auto *nsh_1279 = buffer.data(nsh + 1279);
    const auto *nsh_1280 = buffer.data(nsh + 1280);
    const auto *nsh_1281 = buffer.data(nsh + 1281);
    const auto *nsh_1282 = buffer.data(nsh + 1282);
    const auto *nsh_1283 = buffer.data(nsh + 1283);
    const auto *nsh_1284 = buffer.data(nsh + 1284);
    const auto *nsh_1285 = buffer.data(nsh + 1285);
    const auto *nsh_1286 = buffer.data(nsh + 1286);
    const auto *nsh_1287 = buffer.data(nsh + 1287);
    const auto *nsh_1288 = buffer.data(nsh + 1288);
    const auto *nsh_1289 = buffer.data(nsh + 1289);
    const auto *nsh_1290 = buffer.data(nsh + 1290);
    const auto *nsh_1291 = buffer.data(nsh + 1291);
    const auto *nsh_1292 = buffer.data(nsh + 1292);
    const auto *nsh_1293 = buffer.data(nsh + 1293);
    const auto *nsh_1294 = buffer.data(nsh + 1294);
    const auto *nsh_1295 = buffer.data(nsh + 1295);
    const auto *nsh_1296 = buffer.data(nsh + 1296);
    const auto *nsh_1297 = buffer.data(nsh + 1297);
    const auto *nsh_1298 = buffer.data(nsh + 1298);
    const auto *nsh_1299 = buffer.data(nsh + 1299);
    const auto *nsh_1300 = buffer.data(nsh + 1300);
    const auto *nsh_1301 = buffer.data(nsh + 1301);
    const auto *nsh_1302 = buffer.data(nsh + 1302);
    const auto *nsh_1303 = buffer.data(nsh + 1303);
    const auto *nsh_1304 = buffer.data(nsh + 1304);
    const auto *nsh_1305 = buffer.data(nsh + 1305);

#pragma omp simd aligned(t_1625, t_1626, t_1627, pc_x, nsg0_871, nsg0_872, nsg0_873, nsg1_871, \
                         nsg1_872, nsg1_873, nsh_1219, nsh_1220, \
                         nsh_1221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1625[k] = f_16 * nsg0_871[k]
                    - f_17 * nsg1_871[k]
                    + f_3 * pc_x[k] * nsh_1219[k];

        t_1626[k] = f_16 * nsg0_872[k]
                    - f_17 * nsg1_872[k]
                    + f_3 * pc_x[k] * nsh_1220[k];

        t_1627[k] = f_8 * nsg0_873[k]
                    - f_9 * nsg1_873[k]
                    + f_3 * pc_x[k] * nsh_1221[k];
    }

#pragma omp simd aligned(t_1628, t_1629, t_1630, pc_x, nsg0_874, nsg0_875, nsg0_876, nsg1_874, \
                         nsg1_875, nsg1_876, nsh_1222, nsh_1223, \
                         nsh_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1628[k] = f_8 * nsg0_874[k]
                    - f_9 * nsg1_874[k]
                    + f_3 * pc_x[k] * nsh_1222[k];

        t_1629[k] = f_8 * nsg0_875[k]
                    - f_9 * nsg1_875[k]
                    + f_3 * pc_x[k] * nsh_1223[k];

        t_1630[k] = f_6 * nsg0_876[k]
                    - f_7 * nsg1_876[k]
                    + f_3 * pc_x[k] * nsh_1224[k];
    }

#pragma omp simd aligned(t_1631, t_1632, t_1633, pc_x, nsg0_877, nsg0_878, nsg0_879, nsg1_877, \
                         nsg1_878, nsg1_879, nsh_1225, nsh_1226, \
                         nsh_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1631[k] = f_6 * nsg0_877[k]
                    - f_7 * nsg1_877[k]
                    + f_3 * pc_x[k] * nsh_1225[k];

        t_1632[k] = f_6 * nsg0_878[k]
                    - f_7 * nsg1_878[k]
                    + f_3 * pc_x[k] * nsh_1226[k];

        t_1633[k] = f_6 * nsg0_879[k]
                    - f_7 * nsg1_879[k]
                    + f_3 * pc_x[k] * nsh_1227[k];
    }

#pragma omp simd aligned(t_1634, t_1635, t_1636, pc_x, nsg0_880, nsg0_881, nsg0_882, nsg1_880, \
                         nsg1_881, nsg1_882, nsh_1228, nsh_1229, \
                         nsh_1230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1634[k] = f_4 * nsg0_880[k]
                    - f_5 * nsg1_880[k]
                    + f_3 * pc_x[k] * nsh_1228[k];

        t_1635[k] = f_4 * nsg0_881[k]
                    - f_5 * nsg1_881[k]
                    + f_3 * pc_x[k] * nsh_1229[k];

        t_1636[k] = f_4 * nsg0_882[k]
                    - f_5 * nsg1_882[k]
                    + f_3 * pc_x[k] * nsh_1230[k];
    }

#pragma omp simd aligned(t_1637, t_1638, t_1639, t_1640, t_1641, pc_x, nsg0_883, nsg0_884, \
                         nsg1_883, nsg1_884, nsh_1231, nsh_1232, nsh_1233, nsh_1234, \
                         nsh_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1637[k] = f_4 * nsg0_883[k]
                    - f_5 * nsg1_883[k]
                    + f_3 * pc_x[k] * nsh_1231[k];

        t_1638[k] = f_4 * nsg0_884[k]
                    - f_5 * nsg1_884[k]
                    + f_3 * pc_x[k] * nsh_1232[k];

        t_1639[k] = f_3 * pc_x[k] * nsh_1233[k];

        t_1640[k] = f_3 * pc_x[k] * nsh_1234[k];

        t_1641[k] = f_3 * pc_x[k] * nsh_1235[k];
    }

#pragma omp simd aligned(t_1642, t_1643, t_1644, t_1645, t_1646, pc_x, pc_y, pc_z, msh_1002, \
                         msh_1023, nsg0_880, nsg1_880, nsh_1233, nsh_1236, nsh_1237, \
                         nsh_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1642[k] = f_3 * pc_x[k] * nsh_1236[k];

        t_1643[k] = f_3 * pc_x[k] * nsh_1237[k];

        t_1644[k] = f_3 * pc_x[k] * nsh_1238[k];

        t_1645[k] = f_19 * msh_1023[k]
                    + f_1 * nsg0_880[k]
                    - f_2 * nsg1_880[k]
                    + f_3 * pc_y[k] * nsh_1233[k];

        t_1646[k] = f_13 * msh_1002[k]
                    + f_3 * pc_z[k] * nsh_1233[k];
    }

#pragma omp simd aligned(t_1647, t_1648, t_1649, pc_y, msh_1025, msh_1026, msh_1027, nsg0_882, \
                         nsg0_883, nsg0_884, nsg1_882, nsg1_883, nsg1_884, nsh_1235, nsh_1236, \
                         nsh_1237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1647[k] = f_19 * msh_1025[k]
                    + f_8 * nsg0_882[k]
                    - f_9 * nsg1_882[k]
                    + f_3 * pc_y[k] * nsh_1235[k];

        t_1648[k] = f_19 * msh_1026[k]
                    + f_6 * nsg0_883[k]
                    - f_7 * nsg1_883[k]
                    + f_3 * pc_y[k] * nsh_1236[k];

        t_1649[k] = f_19 * msh_1027[k]
                    + f_4 * nsg0_884[k]
                    - f_5 * nsg1_884[k]
                    + f_3 * pc_y[k] * nsh_1237[k];
    }

#pragma omp simd aligned(t_1650, t_1651, t_1652, pc_x, pc_y, pc_z, msh_1007, msh_1028, \
                         nsg0_884, nsg0_885, nsg1_884, nsg1_885, nsh_1238, \
                         nsh_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1650[k] = f_19 * msh_1028[k]
                    + f_3 * pc_y[k] * nsh_1238[k];

        t_1651[k] = f_13 * msh_1007[k]
                    + f_1 * nsg0_884[k]
                    - f_2 * nsg1_884[k]
                    + f_3 * pc_z[k] * nsh_1238[k];

        t_1652[k] = f_1 * nsg0_885[k]
                    - f_2 * nsg1_885[k]
                    + f_3 * pc_x[k] * nsh_1239[k];
    }

#pragma omp simd aligned(t_1653, t_1654, t_1655, pc_x, nsg0_886, nsg0_887, nsg0_888, nsg1_886, \
                         nsg1_887, nsg1_888, nsh_1240, nsh_1241, \
                         nsh_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1653[k] = f_16 * nsg0_886[k]
                    - f_17 * nsg1_886[k]
                    + f_3 * pc_x[k] * nsh_1240[k];

        t_1654[k] = f_16 * nsg0_887[k]
                    - f_17 * nsg1_887[k]
                    + f_3 * pc_x[k] * nsh_1241[k];

        t_1655[k] = f_8 * nsg0_888[k]
                    - f_9 * nsg1_888[k]
                    + f_3 * pc_x[k] * nsh_1242[k];
    }

#pragma omp simd aligned(t_1656, t_1657, t_1658, pc_x, nsg0_889, nsg0_890, nsg0_891, nsg1_889, \
                         nsg1_890, nsg1_891, nsh_1243, nsh_1244, \
                         nsh_1245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1656[k] = f_8 * nsg0_889[k]
                    - f_9 * nsg1_889[k]
                    + f_3 * pc_x[k] * nsh_1243[k];

        t_1657[k] = f_8 * nsg0_890[k]
                    - f_9 * nsg1_890[k]
                    + f_3 * pc_x[k] * nsh_1244[k];

        t_1658[k] = f_6 * nsg0_891[k]
                    - f_7 * nsg1_891[k]
                    + f_3 * pc_x[k] * nsh_1245[k];
    }

#pragma omp simd aligned(t_1659, t_1660, t_1661, pc_x, nsg0_892, nsg0_893, nsg0_894, nsg1_892, \
                         nsg1_893, nsg1_894, nsh_1246, nsh_1247, \
                         nsh_1248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1659[k] = f_6 * nsg0_892[k]
                    - f_7 * nsg1_892[k]
                    + f_3 * pc_x[k] * nsh_1246[k];

        t_1660[k] = f_6 * nsg0_893[k]
                    - f_7 * nsg1_893[k]
                    + f_3 * pc_x[k] * nsh_1247[k];

        t_1661[k] = f_6 * nsg0_894[k]
                    - f_7 * nsg1_894[k]
                    + f_3 * pc_x[k] * nsh_1248[k];
    }

#pragma omp simd aligned(t_1662, t_1663, t_1664, pc_x, nsg0_895, nsg0_896, nsg0_897, nsg1_895, \
                         nsg1_896, nsg1_897, nsh_1249, nsh_1250, \
                         nsh_1251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1662[k] = f_4 * nsg0_895[k]
                    - f_5 * nsg1_895[k]
                    + f_3 * pc_x[k] * nsh_1249[k];

        t_1663[k] = f_4 * nsg0_896[k]
                    - f_5 * nsg1_896[k]
                    + f_3 * pc_x[k] * nsh_1250[k];

        t_1664[k] = f_4 * nsg0_897[k]
                    - f_5 * nsg1_897[k]
                    + f_3 * pc_x[k] * nsh_1251[k];
    }

#pragma omp simd aligned(t_1665, t_1666, t_1667, t_1668, t_1669, pc_x, nsg0_898, nsg0_899, \
                         nsg1_898, nsg1_899, nsh_1252, nsh_1253, nsh_1254, nsh_1255, \
                         nsh_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1665[k] = f_4 * nsg0_898[k]
                    - f_5 * nsg1_898[k]
                    + f_3 * pc_x[k] * nsh_1252[k];

        t_1666[k] = f_4 * nsg0_899[k]
                    - f_5 * nsg1_899[k]
                    + f_3 * pc_x[k] * nsh_1253[k];

        t_1667[k] = f_3 * pc_x[k] * nsh_1254[k];

        t_1668[k] = f_3 * pc_x[k] * nsh_1255[k];

        t_1669[k] = f_3 * pc_x[k] * nsh_1256[k];
    }

#pragma omp simd aligned(t_1670, t_1671, t_1672, t_1673, t_1674, pc_x, pc_y, pc_z, msh_1023, \
                         msh_1044, nsg0_895, nsg1_895, nsh_1254, nsh_1257, nsh_1258, \
                         nsh_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1670[k] = f_3 * pc_x[k] * nsh_1257[k];

        t_1671[k] = f_3 * pc_x[k] * nsh_1258[k];

        t_1672[k] = f_3 * pc_x[k] * nsh_1259[k];

        t_1673[k] = f_20 * msh_1044[k]
                    + f_1 * nsg0_895[k]
                    - f_2 * nsg1_895[k]
                    + f_3 * pc_y[k] * nsh_1254[k];

        t_1674[k] = f_14 * msh_1023[k]
                    + f_3 * pc_z[k] * nsh_1254[k];
    }

#pragma omp simd aligned(t_1675, t_1676, t_1677, pc_y, msh_1046, msh_1047, msh_1048, nsg0_897, \
                         nsg0_898, nsg0_899, nsg1_897, nsg1_898, nsg1_899, nsh_1256, nsh_1257, \
                         nsh_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1675[k] = f_20 * msh_1046[k]
                    + f_8 * nsg0_897[k]
                    - f_9 * nsg1_897[k]
                    + f_3 * pc_y[k] * nsh_1256[k];

        t_1676[k] = f_20 * msh_1047[k]
                    + f_6 * nsg0_898[k]
                    - f_7 * nsg1_898[k]
                    + f_3 * pc_y[k] * nsh_1257[k];

        t_1677[k] = f_20 * msh_1048[k]
                    + f_4 * nsg0_899[k]
                    - f_5 * nsg1_899[k]
                    + f_3 * pc_y[k] * nsh_1258[k];
    }

#pragma omp simd aligned(t_1678, t_1679, t_1680, pc_x, pc_y, pc_z, msh_1028, msh_1049, \
                         nsg0_899, nsg0_900, nsg1_899, nsg1_900, nsh_1259, \
                         nsh_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1678[k] = f_20 * msh_1049[k]
                    + f_3 * pc_y[k] * nsh_1259[k];

        t_1679[k] = f_14 * msh_1028[k]
                    + f_1 * nsg0_899[k]
                    - f_2 * nsg1_899[k]
                    + f_3 * pc_z[k] * nsh_1259[k];

        t_1680[k] = f_1 * nsg0_900[k]
                    - f_2 * nsg1_900[k]
                    + f_3 * pc_x[k] * nsh_1260[k];
    }

#pragma omp simd aligned(t_1681, t_1682, t_1683, pc_x, nsg0_901, nsg0_902, nsg0_903, nsg1_901, \
                         nsg1_902, nsg1_903, nsh_1261, nsh_1262, \
                         nsh_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1681[k] = f_16 * nsg0_901[k]
                    - f_17 * nsg1_901[k]
                    + f_3 * pc_x[k] * nsh_1261[k];

        t_1682[k] = f_16 * nsg0_902[k]
                    - f_17 * nsg1_902[k]
                    + f_3 * pc_x[k] * nsh_1262[k];

        t_1683[k] = f_8 * nsg0_903[k]
                    - f_9 * nsg1_903[k]
                    + f_3 * pc_x[k] * nsh_1263[k];
    }

#pragma omp simd aligned(t_1684, t_1685, t_1686, pc_x, nsg0_904, nsg0_905, nsg0_906, nsg1_904, \
                         nsg1_905, nsg1_906, nsh_1264, nsh_1265, \
                         nsh_1266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1684[k] = f_8 * nsg0_904[k]
                    - f_9 * nsg1_904[k]
                    + f_3 * pc_x[k] * nsh_1264[k];

        t_1685[k] = f_8 * nsg0_905[k]
                    - f_9 * nsg1_905[k]
                    + f_3 * pc_x[k] * nsh_1265[k];

        t_1686[k] = f_6 * nsg0_906[k]
                    - f_7 * nsg1_906[k]
                    + f_3 * pc_x[k] * nsh_1266[k];
    }

#pragma omp simd aligned(t_1687, t_1688, t_1689, pc_x, nsg0_907, nsg0_908, nsg0_909, nsg1_907, \
                         nsg1_908, nsg1_909, nsh_1267, nsh_1268, \
                         nsh_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1687[k] = f_6 * nsg0_907[k]
                    - f_7 * nsg1_907[k]
                    + f_3 * pc_x[k] * nsh_1267[k];

        t_1688[k] = f_6 * nsg0_908[k]
                    - f_7 * nsg1_908[k]
                    + f_3 * pc_x[k] * nsh_1268[k];

        t_1689[k] = f_6 * nsg0_909[k]
                    - f_7 * nsg1_909[k]
                    + f_3 * pc_x[k] * nsh_1269[k];
    }

#pragma omp simd aligned(t_1690, t_1691, t_1692, pc_x, nsg0_910, nsg0_911, nsg0_912, nsg1_910, \
                         nsg1_911, nsg1_912, nsh_1270, nsh_1271, \
                         nsh_1272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1690[k] = f_4 * nsg0_910[k]
                    - f_5 * nsg1_910[k]
                    + f_3 * pc_x[k] * nsh_1270[k];

        t_1691[k] = f_4 * nsg0_911[k]
                    - f_5 * nsg1_911[k]
                    + f_3 * pc_x[k] * nsh_1271[k];

        t_1692[k] = f_4 * nsg0_912[k]
                    - f_5 * nsg1_912[k]
                    + f_3 * pc_x[k] * nsh_1272[k];
    }

#pragma omp simd aligned(t_1693, t_1694, t_1695, t_1696, t_1697, pc_x, nsg0_913, nsg0_914, \
                         nsg1_913, nsg1_914, nsh_1273, nsh_1274, nsh_1275, nsh_1276, \
                         nsh_1277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1693[k] = f_4 * nsg0_913[k]
                    - f_5 * nsg1_913[k]
                    + f_3 * pc_x[k] * nsh_1273[k];

        t_1694[k] = f_4 * nsg0_914[k]
                    - f_5 * nsg1_914[k]
                    + f_3 * pc_x[k] * nsh_1274[k];

        t_1695[k] = f_3 * pc_x[k] * nsh_1275[k];

        t_1696[k] = f_3 * pc_x[k] * nsh_1276[k];

        t_1697[k] = f_3 * pc_x[k] * nsh_1277[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, t_1701, t_1702, pc_x, pc_y, pc_z, msh_1044, \
                         msh_1065, nsg0_910, nsg1_910, nsh_1275, nsh_1278, nsh_1279, \
                         nsh_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_3 * pc_x[k] * nsh_1278[k];

        t_1699[k] = f_3 * pc_x[k] * nsh_1279[k];

        t_1700[k] = f_3 * pc_x[k] * nsh_1280[k];

        t_1701[k] = f_21 * msh_1065[k]
                    + f_1 * nsg0_910[k]
                    - f_2 * nsg1_910[k]
                    + f_3 * pc_y[k] * nsh_1275[k];

        t_1702[k] = f_21 * msh_1044[k]
                    + f_3 * pc_z[k] * nsh_1275[k];
    }

#pragma omp simd aligned(t_1703, t_1704, t_1705, pc_y, msh_1067, msh_1068, msh_1069, nsg0_912, \
                         nsg0_913, nsg0_914, nsg1_912, nsg1_913, nsg1_914, nsh_1277, nsh_1278, \
                         nsh_1279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1703[k] = f_21 * msh_1067[k]
                    + f_8 * nsg0_912[k]
                    - f_9 * nsg1_912[k]
                    + f_3 * pc_y[k] * nsh_1277[k];

        t_1704[k] = f_21 * msh_1068[k]
                    + f_6 * nsg0_913[k]
                    - f_7 * nsg1_913[k]
                    + f_3 * pc_y[k] * nsh_1278[k];

        t_1705[k] = f_21 * msh_1069[k]
                    + f_4 * nsg0_914[k]
                    - f_5 * nsg1_914[k]
                    + f_3 * pc_y[k] * nsh_1279[k];
    }

#pragma omp simd aligned(t_1706, t_1707, t_1708, pc_x, pc_y, pc_z, msh_1049, msh_1070, \
                         nsg0_914, nsg0_915, nsg1_914, nsg1_915, nsh_1280, \
                         nsh_1281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1706[k] = f_21 * msh_1070[k]
                    + f_3 * pc_y[k] * nsh_1280[k];

        t_1707[k] = f_21 * msh_1049[k]
                    + f_1 * nsg0_914[k]
                    - f_2 * nsg1_914[k]
                    + f_3 * pc_z[k] * nsh_1280[k];

        t_1708[k] = f_1 * nsg0_915[k]
                    - f_2 * nsg1_915[k]
                    + f_3 * pc_x[k] * nsh_1281[k];
    }

#pragma omp simd aligned(t_1709, t_1710, t_1711, pc_x, nsg0_916, nsg0_917, nsg0_918, nsg1_916, \
                         nsg1_917, nsg1_918, nsh_1282, nsh_1283, \
                         nsh_1284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1709[k] = f_16 * nsg0_916[k]
                    - f_17 * nsg1_916[k]
                    + f_3 * pc_x[k] * nsh_1282[k];

        t_1710[k] = f_16 * nsg0_917[k]
                    - f_17 * nsg1_917[k]
                    + f_3 * pc_x[k] * nsh_1283[k];

        t_1711[k] = f_8 * nsg0_918[k]
                    - f_9 * nsg1_918[k]
                    + f_3 * pc_x[k] * nsh_1284[k];
    }

#pragma omp simd aligned(t_1712, t_1713, t_1714, pc_x, nsg0_919, nsg0_920, nsg0_921, nsg1_919, \
                         nsg1_920, nsg1_921, nsh_1285, nsh_1286, \
                         nsh_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1712[k] = f_8 * nsg0_919[k]
                    - f_9 * nsg1_919[k]
                    + f_3 * pc_x[k] * nsh_1285[k];

        t_1713[k] = f_8 * nsg0_920[k]
                    - f_9 * nsg1_920[k]
                    + f_3 * pc_x[k] * nsh_1286[k];

        t_1714[k] = f_6 * nsg0_921[k]
                    - f_7 * nsg1_921[k]
                    + f_3 * pc_x[k] * nsh_1287[k];
    }

#pragma omp simd aligned(t_1715, t_1716, t_1717, pc_x, nsg0_922, nsg0_923, nsg0_924, nsg1_922, \
                         nsg1_923, nsg1_924, nsh_1288, nsh_1289, \
                         nsh_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1715[k] = f_6 * nsg0_922[k]
                    - f_7 * nsg1_922[k]
                    + f_3 * pc_x[k] * nsh_1288[k];

        t_1716[k] = f_6 * nsg0_923[k]
                    - f_7 * nsg1_923[k]
                    + f_3 * pc_x[k] * nsh_1289[k];

        t_1717[k] = f_6 * nsg0_924[k]
                    - f_7 * nsg1_924[k]
                    + f_3 * pc_x[k] * nsh_1290[k];
    }

#pragma omp simd aligned(t_1718, t_1719, t_1720, pc_x, nsg0_925, nsg0_926, nsg0_927, nsg1_925, \
                         nsg1_926, nsg1_927, nsh_1291, nsh_1292, \
                         nsh_1293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1718[k] = f_4 * nsg0_925[k]
                    - f_5 * nsg1_925[k]
                    + f_3 * pc_x[k] * nsh_1291[k];

        t_1719[k] = f_4 * nsg0_926[k]
                    - f_5 * nsg1_926[k]
                    + f_3 * pc_x[k] * nsh_1292[k];

        t_1720[k] = f_4 * nsg0_927[k]
                    - f_5 * nsg1_927[k]
                    + f_3 * pc_x[k] * nsh_1293[k];
    }

#pragma omp simd aligned(t_1721, t_1722, t_1723, t_1724, t_1725, pc_x, nsg0_928, nsg0_929, \
                         nsg1_928, nsg1_929, nsh_1294, nsh_1295, nsh_1296, nsh_1297, \
                         nsh_1298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1721[k] = f_4 * nsg0_928[k]
                    - f_5 * nsg1_928[k]
                    + f_3 * pc_x[k] * nsh_1294[k];

        t_1722[k] = f_4 * nsg0_929[k]
                    - f_5 * nsg1_929[k]
                    + f_3 * pc_x[k] * nsh_1295[k];

        t_1723[k] = f_3 * pc_x[k] * nsh_1296[k];

        t_1724[k] = f_3 * pc_x[k] * nsh_1297[k];

        t_1725[k] = f_3 * pc_x[k] * nsh_1298[k];
    }

#pragma omp simd aligned(t_1726, t_1727, t_1728, t_1729, t_1730, pc_x, pc_y, pc_z, msh_1065, \
                         msh_1086, nsg0_925, nsg1_925, nsh_1296, nsh_1299, nsh_1300, \
                         nsh_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1726[k] = f_3 * pc_x[k] * nsh_1299[k];

        t_1727[k] = f_3 * pc_x[k] * nsh_1300[k];

        t_1728[k] = f_3 * pc_x[k] * nsh_1301[k];

        t_1729[k] = f_14 * msh_1086[k]
                    + f_1 * nsg0_925[k]
                    - f_2 * nsg1_925[k]
                    + f_3 * pc_y[k] * nsh_1296[k];

        t_1730[k] = f_20 * msh_1065[k]
                    + f_3 * pc_z[k] * nsh_1296[k];
    }

#pragma omp simd aligned(t_1731, t_1732, t_1733, pc_y, msh_1088, msh_1089, msh_1090, nsg0_927, \
                         nsg0_928, nsg0_929, nsg1_927, nsg1_928, nsg1_929, nsh_1298, nsh_1299, \
                         nsh_1300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1731[k] = f_14 * msh_1088[k]
                    + f_8 * nsg0_927[k]
                    - f_9 * nsg1_927[k]
                    + f_3 * pc_y[k] * nsh_1298[k];

        t_1732[k] = f_14 * msh_1089[k]
                    + f_6 * nsg0_928[k]
                    - f_7 * nsg1_928[k]
                    + f_3 * pc_y[k] * nsh_1299[k];

        t_1733[k] = f_14 * msh_1090[k]
                    + f_4 * nsg0_929[k]
                    - f_5 * nsg1_929[k]
                    + f_3 * pc_y[k] * nsh_1300[k];
    }

#pragma omp simd aligned(t_1734, t_1735, t_1736, pc_x, pc_y, pc_z, msh_1070, msh_1091, \
                         nsg0_929, nsg0_930, nsg1_929, nsg1_930, nsh_1301, \
                         nsh_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1734[k] = f_14 * msh_1091[k]
                    + f_3 * pc_y[k] * nsh_1301[k];

        t_1735[k] = f_20 * msh_1070[k]
                    + f_1 * nsg0_929[k]
                    - f_2 * nsg1_929[k]
                    + f_3 * pc_z[k] * nsh_1301[k];

        t_1736[k] = f_1 * nsg0_930[k]
                    - f_2 * nsg1_930[k]
                    + f_3 * pc_x[k] * nsh_1302[k];
    }

#pragma omp simd aligned(t_1737, t_1738, t_1739, pc_x, nsg0_931, nsg0_932, nsg0_933, nsg1_931, \
                         nsg1_932, nsg1_933, nsh_1303, nsh_1304, \
                         nsh_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1737[k] = f_16 * nsg0_931[k]
                    - f_17 * nsg1_931[k]
                    + f_3 * pc_x[k] * nsh_1303[k];

        t_1738[k] = f_16 * nsg0_932[k]
                    - f_17 * nsg1_932[k]
                    + f_3 * pc_x[k] * nsh_1304[k];

        t_1739[k] = f_8 * nsg0_933[k]
                    - f_9 * nsg1_933[k]
                    + f_3 * pc_x[k] * nsh_1305[k];
    }
}

static auto
compute_prim_nsi_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t msi0,
                                                           const size_t msh, const size_t msi1,
                                                           const size_t nsg0, const size_t nsg1,
                                                           const size_t nsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 4.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 4.0 / q;
    const auto f_19 = 3.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *msi0_1512 = buffer.data(msi0 + 1512);
    const auto *msi0_1514 = buffer.data(msi0 + 1514);
    const auto *msi0_1517 = buffer.data(msi0 + 1517);
    const auto *msi0_1521 = buffer.data(msi0 + 1521);
    const auto *msi0_1526 = buffer.data(msi0 + 1526);
    const auto *msi0_1533 = buffer.data(msi0 + 1533);
    const auto *msi0_1535 = buffer.data(msi0 + 1535);
    const auto *msi0_1536 = buffer.data(msi0 + 1536);
    const auto *msi0_1537 = buffer.data(msi0 + 1537);
    const auto *msi0_1539 = buffer.data(msi0 + 1539);

    const auto *msh_1086 = buffer.data(msh + 1086);
    const auto *msh_1091 = buffer.data(msh + 1091);
    const auto *msh_1107 = buffer.data(msh + 1107);
    const auto *msh_1109 = buffer.data(msh + 1109);
    const auto *msh_1110 = buffer.data(msh + 1110);
    const auto *msh_1111 = buffer.data(msh + 1111);
    const auto *msh_1112 = buffer.data(msh + 1112);
    const auto *msh_1128 = buffer.data(msh + 1128);
    const auto *msh_1130 = buffer.data(msh + 1130);
    const auto *msh_1131 = buffer.data(msh + 1131);
    const auto *msh_1132 = buffer.data(msh + 1132);
    const auto *msh_1133 = buffer.data(msh + 1133);
    const auto *msh_1149 = buffer.data(msh + 1149);
    const auto *msh_1151 = buffer.data(msh + 1151);
    const auto *msh_1152 = buffer.data(msh + 1152);
    const auto *msh_1153 = buffer.data(msh + 1153);
    const auto *msh_1154 = buffer.data(msh + 1154);

    const auto *msi1_1512 = buffer.data(msi1 + 1512);
    const auto *msi1_1514 = buffer.data(msi1 + 1514);
    const auto *msi1_1517 = buffer.data(msi1 + 1517);
    const auto *msi1_1521 = buffer.data(msi1 + 1521);
    const auto *msi1_1526 = buffer.data(msi1 + 1526);
    const auto *msi1_1533 = buffer.data(msi1 + 1533);
    const auto *msi1_1535 = buffer.data(msi1 + 1535);
    const auto *msi1_1536 = buffer.data(msi1 + 1536);
    const auto *msi1_1537 = buffer.data(msi1 + 1537);
    const auto *msi1_1539 = buffer.data(msi1 + 1539);

    const auto *nsg0_934 = buffer.data(nsg0 + 934);
    const auto *nsg0_935 = buffer.data(nsg0 + 935);
    const auto *nsg0_936 = buffer.data(nsg0 + 936);
    const auto *nsg0_937 = buffer.data(nsg0 + 937);
    const auto *nsg0_938 = buffer.data(nsg0 + 938);
    const auto *nsg0_939 = buffer.data(nsg0 + 939);
    const auto *nsg0_940 = buffer.data(nsg0 + 940);
    const auto *nsg0_941 = buffer.data(nsg0 + 941);
    const auto *nsg0_942 = buffer.data(nsg0 + 942);
    const auto *nsg0_943 = buffer.data(nsg0 + 943);
    const auto *nsg0_944 = buffer.data(nsg0 + 944);
    const auto *nsg0_945 = buffer.data(nsg0 + 945);
    const auto *nsg0_946 = buffer.data(nsg0 + 946);
    const auto *nsg0_947 = buffer.data(nsg0 + 947);
    const auto *nsg0_948 = buffer.data(nsg0 + 948);
    const auto *nsg0_949 = buffer.data(nsg0 + 949);
    const auto *nsg0_950 = buffer.data(nsg0 + 950);
    const auto *nsg0_951 = buffer.data(nsg0 + 951);
    const auto *nsg0_952 = buffer.data(nsg0 + 952);
    const auto *nsg0_953 = buffer.data(nsg0 + 953);
    const auto *nsg0_954 = buffer.data(nsg0 + 954);
    const auto *nsg0_955 = buffer.data(nsg0 + 955);
    const auto *nsg0_956 = buffer.data(nsg0 + 956);
    const auto *nsg0_957 = buffer.data(nsg0 + 957);
    const auto *nsg0_958 = buffer.data(nsg0 + 958);
    const auto *nsg0_959 = buffer.data(nsg0 + 959);
    const auto *nsg0_961 = buffer.data(nsg0 + 961);
    const auto *nsg0_963 = buffer.data(nsg0 + 963);
    const auto *nsg0_964 = buffer.data(nsg0 + 964);
    const auto *nsg0_966 = buffer.data(nsg0 + 966);
    const auto *nsg0_967 = buffer.data(nsg0 + 967);
    const auto *nsg0_968 = buffer.data(nsg0 + 968);
    const auto *nsg0_970 = buffer.data(nsg0 + 970);
    const auto *nsg0_971 = buffer.data(nsg0 + 971);
    const auto *nsg0_972 = buffer.data(nsg0 + 972);
    const auto *nsg0_973 = buffer.data(nsg0 + 973);
    const auto *nsg0_975 = buffer.data(nsg0 + 975);
    const auto *nsg0_977 = buffer.data(nsg0 + 977);
    const auto *nsg0_978 = buffer.data(nsg0 + 978);
    const auto *nsg0_980 = buffer.data(nsg0 + 980);
    const auto *nsg0_981 = buffer.data(nsg0 + 981);
    const auto *nsg0_982 = buffer.data(nsg0 + 982);
    const auto *nsg0_984 = buffer.data(nsg0 + 984);
    const auto *nsg0_985 = buffer.data(nsg0 + 985);
    const auto *nsg0_986 = buffer.data(nsg0 + 986);
    const auto *nsg0_987 = buffer.data(nsg0 + 987);
    const auto *nsg0_988 = buffer.data(nsg0 + 988);
    const auto *nsg0_989 = buffer.data(nsg0 + 989);

    const auto *nsg1_934 = buffer.data(nsg1 + 934);
    const auto *nsg1_935 = buffer.data(nsg1 + 935);
    const auto *nsg1_936 = buffer.data(nsg1 + 936);
    const auto *nsg1_937 = buffer.data(nsg1 + 937);
    const auto *nsg1_938 = buffer.data(nsg1 + 938);
    const auto *nsg1_939 = buffer.data(nsg1 + 939);
    const auto *nsg1_940 = buffer.data(nsg1 + 940);
    const auto *nsg1_941 = buffer.data(nsg1 + 941);
    const auto *nsg1_942 = buffer.data(nsg1 + 942);
    const auto *nsg1_943 = buffer.data(nsg1 + 943);
    const auto *nsg1_944 = buffer.data(nsg1 + 944);
    const auto *nsg1_945 = buffer.data(nsg1 + 945);
    const auto *nsg1_946 = buffer.data(nsg1 + 946);
    const auto *nsg1_947 = buffer.data(nsg1 + 947);
    const auto *nsg1_948 = buffer.data(nsg1 + 948);
    const auto *nsg1_949 = buffer.data(nsg1 + 949);
    const auto *nsg1_950 = buffer.data(nsg1 + 950);
    const auto *nsg1_951 = buffer.data(nsg1 + 951);
    const auto *nsg1_952 = buffer.data(nsg1 + 952);
    const auto *nsg1_953 = buffer.data(nsg1 + 953);
    const auto *nsg1_954 = buffer.data(nsg1 + 954);
    const auto *nsg1_955 = buffer.data(nsg1 + 955);
    const auto *nsg1_956 = buffer.data(nsg1 + 956);
    const auto *nsg1_957 = buffer.data(nsg1 + 957);
    const auto *nsg1_958 = buffer.data(nsg1 + 958);
    const auto *nsg1_959 = buffer.data(nsg1 + 959);
    const auto *nsg1_961 = buffer.data(nsg1 + 961);
    const auto *nsg1_963 = buffer.data(nsg1 + 963);
    const auto *nsg1_964 = buffer.data(nsg1 + 964);
    const auto *nsg1_966 = buffer.data(nsg1 + 966);
    const auto *nsg1_967 = buffer.data(nsg1 + 967);
    const auto *nsg1_968 = buffer.data(nsg1 + 968);
    const auto *nsg1_970 = buffer.data(nsg1 + 970);
    const auto *nsg1_971 = buffer.data(nsg1 + 971);
    const auto *nsg1_972 = buffer.data(nsg1 + 972);
    const auto *nsg1_973 = buffer.data(nsg1 + 973);
    const auto *nsg1_975 = buffer.data(nsg1 + 975);
    const auto *nsg1_977 = buffer.data(nsg1 + 977);
    const auto *nsg1_978 = buffer.data(nsg1 + 978);
    const auto *nsg1_980 = buffer.data(nsg1 + 980);
    const auto *nsg1_981 = buffer.data(nsg1 + 981);
    const auto *nsg1_982 = buffer.data(nsg1 + 982);
    const auto *nsg1_984 = buffer.data(nsg1 + 984);
    const auto *nsg1_985 = buffer.data(nsg1 + 985);
    const auto *nsg1_986 = buffer.data(nsg1 + 986);
    const auto *nsg1_987 = buffer.data(nsg1 + 987);
    const auto *nsg1_988 = buffer.data(nsg1 + 988);
    const auto *nsg1_989 = buffer.data(nsg1 + 989);

    const auto *nsh_1306 = buffer.data(nsh + 1306);
    const auto *nsh_1307 = buffer.data(nsh + 1307);
    const auto *nsh_1308 = buffer.data(nsh + 1308);
    const auto *nsh_1309 = buffer.data(nsh + 1309);
    const auto *nsh_1310 = buffer.data(nsh + 1310);
    const auto *nsh_1311 = buffer.data(nsh + 1311);
    const auto *nsh_1312 = buffer.data(nsh + 1312);
    const auto *nsh_1313 = buffer.data(nsh + 1313);
    const auto *nsh_1314 = buffer.data(nsh + 1314);
    const auto *nsh_1315 = buffer.data(nsh + 1315);
    const auto *nsh_1316 = buffer.data(nsh + 1316);
    const auto *nsh_1317 = buffer.data(nsh + 1317);
    const auto *nsh_1318 = buffer.data(nsh + 1318);
    const auto *nsh_1319 = buffer.data(nsh + 1319);
    const auto *nsh_1320 = buffer.data(nsh + 1320);
    const auto *nsh_1321 = buffer.data(nsh + 1321);
    const auto *nsh_1322 = buffer.data(nsh + 1322);
    const auto *nsh_1323 = buffer.data(nsh + 1323);
    const auto *nsh_1324 = buffer.data(nsh + 1324);
    const auto *nsh_1325 = buffer.data(nsh + 1325);
    const auto *nsh_1326 = buffer.data(nsh + 1326);
    const auto *nsh_1327 = buffer.data(nsh + 1327);
    const auto *nsh_1328 = buffer.data(nsh + 1328);
    const auto *nsh_1329 = buffer.data(nsh + 1329);
    const auto *nsh_1330 = buffer.data(nsh + 1330);
    const auto *nsh_1331 = buffer.data(nsh + 1331);
    const auto *nsh_1332 = buffer.data(nsh + 1332);
    const auto *nsh_1333 = buffer.data(nsh + 1333);
    const auto *nsh_1334 = buffer.data(nsh + 1334);
    const auto *nsh_1335 = buffer.data(nsh + 1335);
    const auto *nsh_1336 = buffer.data(nsh + 1336);
    const auto *nsh_1337 = buffer.data(nsh + 1337);
    const auto *nsh_1338 = buffer.data(nsh + 1338);
    const auto *nsh_1339 = buffer.data(nsh + 1339);
    const auto *nsh_1340 = buffer.data(nsh + 1340);
    const auto *nsh_1341 = buffer.data(nsh + 1341);
    const auto *nsh_1342 = buffer.data(nsh + 1342);
    const auto *nsh_1343 = buffer.data(nsh + 1343);
    const auto *nsh_1345 = buffer.data(nsh + 1345);
    const auto *nsh_1347 = buffer.data(nsh + 1347);
    const auto *nsh_1348 = buffer.data(nsh + 1348);
    const auto *nsh_1350 = buffer.data(nsh + 1350);
    const auto *nsh_1351 = buffer.data(nsh + 1351);
    const auto *nsh_1352 = buffer.data(nsh + 1352);
    const auto *nsh_1354 = buffer.data(nsh + 1354);
    const auto *nsh_1355 = buffer.data(nsh + 1355);
    const auto *nsh_1356 = buffer.data(nsh + 1356);
    const auto *nsh_1357 = buffer.data(nsh + 1357);
    const auto *nsh_1359 = buffer.data(nsh + 1359);
    const auto *nsh_1360 = buffer.data(nsh + 1360);
    const auto *nsh_1361 = buffer.data(nsh + 1361);
    const auto *nsh_1362 = buffer.data(nsh + 1362);
    const auto *nsh_1363 = buffer.data(nsh + 1363);
    const auto *nsh_1364 = buffer.data(nsh + 1364);
    const auto *nsh_1365 = buffer.data(nsh + 1365);
    const auto *nsh_1367 = buffer.data(nsh + 1367);
    const auto *nsh_1368 = buffer.data(nsh + 1368);
    const auto *nsh_1370 = buffer.data(nsh + 1370);
    const auto *nsh_1371 = buffer.data(nsh + 1371);
    const auto *nsh_1372 = buffer.data(nsh + 1372);
    const auto *nsh_1374 = buffer.data(nsh + 1374);
    const auto *nsh_1375 = buffer.data(nsh + 1375);
    const auto *nsh_1376 = buffer.data(nsh + 1376);
    const auto *nsh_1377 = buffer.data(nsh + 1377);
    const auto *nsh_1379 = buffer.data(nsh + 1379);
    const auto *nsh_1380 = buffer.data(nsh + 1380);
    const auto *nsh_1381 = buffer.data(nsh + 1381);
    const auto *nsh_1382 = buffer.data(nsh + 1382);
    const auto *nsh_1383 = buffer.data(nsh + 1383);
    const auto *nsh_1384 = buffer.data(nsh + 1384);
    const auto *nsh_1385 = buffer.data(nsh + 1385);

#pragma omp simd aligned(t_1740, t_1741, t_1742, pc_x, nsg0_934, nsg0_935, nsg0_936, nsg1_934, \
                         nsg1_935, nsg1_936, nsh_1306, nsh_1307, \
                         nsh_1308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1740[k] = f_8 * nsg0_934[k]
                    - f_9 * nsg1_934[k]
                    + f_3 * pc_x[k] * nsh_1306[k];

        t_1741[k] = f_8 * nsg0_935[k]
                    - f_9 * nsg1_935[k]
                    + f_3 * pc_x[k] * nsh_1307[k];

        t_1742[k] = f_6 * nsg0_936[k]
                    - f_7 * nsg1_936[k]
                    + f_3 * pc_x[k] * nsh_1308[k];
    }

#pragma omp simd aligned(t_1743, t_1744, t_1745, pc_x, nsg0_937, nsg0_938, nsg0_939, nsg1_937, \
                         nsg1_938, nsg1_939, nsh_1309, nsh_1310, \
                         nsh_1311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1743[k] = f_6 * nsg0_937[k]
                    - f_7 * nsg1_937[k]
                    + f_3 * pc_x[k] * nsh_1309[k];

        t_1744[k] = f_6 * nsg0_938[k]
                    - f_7 * nsg1_938[k]
                    + f_3 * pc_x[k] * nsh_1310[k];

        t_1745[k] = f_6 * nsg0_939[k]
                    - f_7 * nsg1_939[k]
                    + f_3 * pc_x[k] * nsh_1311[k];
    }

#pragma omp simd aligned(t_1746, t_1747, t_1748, pc_x, nsg0_940, nsg0_941, nsg0_942, nsg1_940, \
                         nsg1_941, nsg1_942, nsh_1312, nsh_1313, \
                         nsh_1314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1746[k] = f_4 * nsg0_940[k]
                    - f_5 * nsg1_940[k]
                    + f_3 * pc_x[k] * nsh_1312[k];

        t_1747[k] = f_4 * nsg0_941[k]
                    - f_5 * nsg1_941[k]
                    + f_3 * pc_x[k] * nsh_1313[k];

        t_1748[k] = f_4 * nsg0_942[k]
                    - f_5 * nsg1_942[k]
                    + f_3 * pc_x[k] * nsh_1314[k];
    }

#pragma omp simd aligned(t_1749, t_1750, t_1751, t_1752, t_1753, pc_x, nsg0_943, nsg0_944, \
                         nsg1_943, nsg1_944, nsh_1315, nsh_1316, nsh_1317, nsh_1318, \
                         nsh_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1749[k] = f_4 * nsg0_943[k]
                    - f_5 * nsg1_943[k]
                    + f_3 * pc_x[k] * nsh_1315[k];

        t_1750[k] = f_4 * nsg0_944[k]
                    - f_5 * nsg1_944[k]
                    + f_3 * pc_x[k] * nsh_1316[k];

        t_1751[k] = f_3 * pc_x[k] * nsh_1317[k];

        t_1752[k] = f_3 * pc_x[k] * nsh_1318[k];

        t_1753[k] = f_3 * pc_x[k] * nsh_1319[k];
    }

#pragma omp simd aligned(t_1754, t_1755, t_1756, t_1757, t_1758, pc_x, pc_y, pc_z, msh_1086, \
                         msh_1107, nsg0_940, nsg1_940, nsh_1317, nsh_1320, nsh_1321, \
                         nsh_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1754[k] = f_3 * pc_x[k] * nsh_1320[k];

        t_1755[k] = f_3 * pc_x[k] * nsh_1321[k];

        t_1756[k] = f_3 * pc_x[k] * nsh_1322[k];

        t_1757[k] = f_13 * msh_1107[k]
                    + f_1 * nsg0_940[k]
                    - f_2 * nsg1_940[k]
                    + f_3 * pc_y[k] * nsh_1317[k];

        t_1758[k] = f_19 * msh_1086[k]
                    + f_3 * pc_z[k] * nsh_1317[k];
    }

#pragma omp simd aligned(t_1759, t_1760, t_1761, pc_y, msh_1109, msh_1110, msh_1111, nsg0_942, \
                         nsg0_943, nsg0_944, nsg1_942, nsg1_943, nsg1_944, nsh_1319, nsh_1320, \
                         nsh_1321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1759[k] = f_13 * msh_1109[k]
                    + f_8 * nsg0_942[k]
                    - f_9 * nsg1_942[k]
                    + f_3 * pc_y[k] * nsh_1319[k];

        t_1760[k] = f_13 * msh_1110[k]
                    + f_6 * nsg0_943[k]
                    - f_7 * nsg1_943[k]
                    + f_3 * pc_y[k] * nsh_1320[k];

        t_1761[k] = f_13 * msh_1111[k]
                    + f_4 * nsg0_944[k]
                    - f_5 * nsg1_944[k]
                    + f_3 * pc_y[k] * nsh_1321[k];
    }

#pragma omp simd aligned(t_1762, t_1763, t_1764, pc_x, pc_y, pc_z, msh_1091, msh_1112, \
                         nsg0_944, nsg0_945, nsg1_944, nsg1_945, nsh_1322, \
                         nsh_1323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1762[k] = f_13 * msh_1112[k]
                    + f_3 * pc_y[k] * nsh_1322[k];

        t_1763[k] = f_19 * msh_1091[k]
                    + f_1 * nsg0_944[k]
                    - f_2 * nsg1_944[k]
                    + f_3 * pc_z[k] * nsh_1322[k];

        t_1764[k] = f_1 * nsg0_945[k]
                    - f_2 * nsg1_945[k]
                    + f_3 * pc_x[k] * nsh_1323[k];
    }

#pragma omp simd aligned(t_1765, t_1766, t_1767, pc_x, nsg0_946, nsg0_947, nsg0_948, nsg1_946, \
                         nsg1_947, nsg1_948, nsh_1324, nsh_1325, \
                         nsh_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1765[k] = f_16 * nsg0_946[k]
                    - f_17 * nsg1_946[k]
                    + f_3 * pc_x[k] * nsh_1324[k];

        t_1766[k] = f_16 * nsg0_947[k]
                    - f_17 * nsg1_947[k]
                    + f_3 * pc_x[k] * nsh_1325[k];

        t_1767[k] = f_8 * nsg0_948[k]
                    - f_9 * nsg1_948[k]
                    + f_3 * pc_x[k] * nsh_1326[k];
    }

#pragma omp simd aligned(t_1768, t_1769, t_1770, pc_x, nsg0_949, nsg0_950, nsg0_951, nsg1_949, \
                         nsg1_950, nsg1_951, nsh_1327, nsh_1328, \
                         nsh_1329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1768[k] = f_8 * nsg0_949[k]
                    - f_9 * nsg1_949[k]
                    + f_3 * pc_x[k] * nsh_1327[k];

        t_1769[k] = f_8 * nsg0_950[k]
                    - f_9 * nsg1_950[k]
                    + f_3 * pc_x[k] * nsh_1328[k];

        t_1770[k] = f_6 * nsg0_951[k]
                    - f_7 * nsg1_951[k]
                    + f_3 * pc_x[k] * nsh_1329[k];
    }

#pragma omp simd aligned(t_1771, t_1772, t_1773, pc_x, nsg0_952, nsg0_953, nsg0_954, nsg1_952, \
                         nsg1_953, nsg1_954, nsh_1330, nsh_1331, \
                         nsh_1332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1771[k] = f_6 * nsg0_952[k]
                    - f_7 * nsg1_952[k]
                    + f_3 * pc_x[k] * nsh_1330[k];

        t_1772[k] = f_6 * nsg0_953[k]
                    - f_7 * nsg1_953[k]
                    + f_3 * pc_x[k] * nsh_1331[k];

        t_1773[k] = f_6 * nsg0_954[k]
                    - f_7 * nsg1_954[k]
                    + f_3 * pc_x[k] * nsh_1332[k];
    }

#pragma omp simd aligned(t_1774, t_1775, t_1776, pc_x, nsg0_955, nsg0_956, nsg0_957, nsg1_955, \
                         nsg1_956, nsg1_957, nsh_1333, nsh_1334, \
                         nsh_1335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1774[k] = f_4 * nsg0_955[k]
                    - f_5 * nsg1_955[k]
                    + f_3 * pc_x[k] * nsh_1333[k];

        t_1775[k] = f_4 * nsg0_956[k]
                    - f_5 * nsg1_956[k]
                    + f_3 * pc_x[k] * nsh_1334[k];

        t_1776[k] = f_4 * nsg0_957[k]
                    - f_5 * nsg1_957[k]
                    + f_3 * pc_x[k] * nsh_1335[k];
    }

#pragma omp simd aligned(t_1777, t_1778, t_1779, t_1780, t_1781, pc_x, nsg0_958, nsg0_959, \
                         nsg1_958, nsg1_959, nsh_1336, nsh_1337, nsh_1338, nsh_1339, \
                         nsh_1340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1777[k] = f_4 * nsg0_958[k]
                    - f_5 * nsg1_958[k]
                    + f_3 * pc_x[k] * nsh_1336[k];

        t_1778[k] = f_4 * nsg0_959[k]
                    - f_5 * nsg1_959[k]
                    + f_3 * pc_x[k] * nsh_1337[k];

        t_1779[k] = f_3 * pc_x[k] * nsh_1338[k];

        t_1780[k] = f_3 * pc_x[k] * nsh_1339[k];

        t_1781[k] = f_3 * pc_x[k] * nsh_1340[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, t_1785, t_1786, pc_x, pc_y, pc_z, msh_1107, \
                         msh_1128, nsg0_955, nsg1_955, nsh_1338, nsh_1341, nsh_1342, \
                         nsh_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = f_3 * pc_x[k] * nsh_1341[k];

        t_1783[k] = f_3 * pc_x[k] * nsh_1342[k];

        t_1784[k] = f_3 * pc_x[k] * nsh_1343[k];

        t_1785[k] = f_12 * msh_1128[k]
                    + f_1 * nsg0_955[k]
                    - f_2 * nsg1_955[k]
                    + f_3 * pc_y[k] * nsh_1338[k];

        t_1786[k] = f_18 * msh_1107[k]
                    + f_3 * pc_z[k] * nsh_1338[k];
    }

#pragma omp simd aligned(t_1787, t_1788, t_1789, pc_y, msh_1130, msh_1131, msh_1132, nsg0_957, \
                         nsg0_958, nsg0_959, nsg1_957, nsg1_958, nsg1_959, nsh_1340, nsh_1341, \
                         nsh_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1787[k] = f_12 * msh_1130[k]
                    + f_8 * nsg0_957[k]
                    - f_9 * nsg1_957[k]
                    + f_3 * pc_y[k] * nsh_1340[k];

        t_1788[k] = f_12 * msh_1131[k]
                    + f_6 * nsg0_958[k]
                    - f_7 * nsg1_958[k]
                    + f_3 * pc_y[k] * nsh_1341[k];

        t_1789[k] = f_12 * msh_1132[k]
                    + f_4 * nsg0_959[k]
                    - f_5 * nsg1_959[k]
                    + f_3 * pc_y[k] * nsh_1342[k];
    }

#pragma omp simd aligned(t_1790, t_1791, t_1792, pa_y, pc_y, pc_z, msi0_1512, msh_1112, \
                         msh_1133, msi1_1512, nsg0_959, nsg1_959, \
                         nsh_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1790[k] = f_12 * msh_1133[k]
                    + f_3 * pc_y[k] * nsh_1343[k];

        t_1791[k] = f_18 * msh_1112[k]
                    + f_1 * nsg0_959[k]
                    - f_2 * nsg1_959[k]
                    + f_3 * pc_z[k] * nsh_1343[k];

        t_1792[k] = pa_y[k] * msi0_1512[k]
                    - f_10 * pc_y[k] * msi1_1512[k];
    }

#pragma omp simd aligned(t_1793, t_1794, t_1795, pa_y, pc_x, pc_y, msi0_1514, msi1_1514, \
                         nsg0_961, nsg0_963, nsg1_961, nsg1_963, nsh_1345, \
                         nsh_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1793[k] = f_16 * nsg0_961[k]
                    - f_17 * nsg1_961[k]
                    + f_3 * pc_x[k] * nsh_1345[k];

        t_1794[k] = pa_y[k] * msi0_1514[k]
                    - f_10 * pc_y[k] * msi1_1514[k];

        t_1795[k] = f_8 * nsg0_963[k]
                    - f_9 * nsg1_963[k]
                    + f_3 * pc_x[k] * nsh_1347[k];
    }

#pragma omp simd aligned(t_1796, t_1797, t_1798, pa_y, pc_x, pc_y, msi0_1517, msi1_1517, \
                         nsg0_964, nsg0_966, nsg1_964, nsg1_966, nsh_1348, \
                         nsh_1350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1796[k] = f_8 * nsg0_964[k]
                    - f_9 * nsg1_964[k]
                    + f_3 * pc_x[k] * nsh_1348[k];

        t_1797[k] = pa_y[k] * msi0_1517[k]
                    - f_10 * pc_y[k] * msi1_1517[k];

        t_1798[k] = f_6 * nsg0_966[k]
                    - f_7 * nsg1_966[k]
                    + f_3 * pc_x[k] * nsh_1350[k];
    }

#pragma omp simd aligned(t_1799, t_1800, t_1801, pa_y, pc_x, pc_y, msi0_1521, msi1_1521, \
                         nsg0_967, nsg0_968, nsg1_967, nsg1_968, nsh_1351, \
                         nsh_1352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1799[k] = f_6 * nsg0_967[k]
                    - f_7 * nsg1_967[k]
                    + f_3 * pc_x[k] * nsh_1351[k];

        t_1800[k] = f_6 * nsg0_968[k]
                    - f_7 * nsg1_968[k]
                    + f_3 * pc_x[k] * nsh_1352[k];

        t_1801[k] = pa_y[k] * msi0_1521[k]
                    - f_10 * pc_y[k] * msi1_1521[k];
    }

#pragma omp simd aligned(t_1802, t_1803, t_1804, pc_x, nsg0_970, nsg0_971, nsg0_972, nsg1_970, \
                         nsg1_971, nsg1_972, nsh_1354, nsh_1355, \
                         nsh_1356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1802[k] = f_4 * nsg0_970[k]
                    - f_5 * nsg1_970[k]
                    + f_3 * pc_x[k] * nsh_1354[k];

        t_1803[k] = f_4 * nsg0_971[k]
                    - f_5 * nsg1_971[k]
                    + f_3 * pc_x[k] * nsh_1355[k];

        t_1804[k] = f_4 * nsg0_972[k]
                    - f_5 * nsg1_972[k]
                    + f_3 * pc_x[k] * nsh_1356[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, t_1808, t_1809, pa_y, pc_x, pc_y, msi0_1526, \
                         msi1_1526, nsg0_973, nsg1_973, nsh_1357, nsh_1359, nsh_1360, \
                         nsh_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = f_4 * nsg0_973[k]
                    - f_5 * nsg1_973[k]
                    + f_3 * pc_x[k] * nsh_1357[k];

        t_1806[k] = pa_y[k] * msi0_1526[k]
                    - f_10 * pc_y[k] * msi1_1526[k];

        t_1807[k] = f_3 * pc_x[k] * nsh_1359[k];

        t_1808[k] = f_3 * pc_x[k] * nsh_1360[k];

        t_1809[k] = f_3 * pc_x[k] * nsh_1361[k];
    }

#pragma omp simd aligned(t_1810, t_1811, t_1812, t_1813, pa_y, pc_x, pc_y, msi0_1533, \
                         msh_1149, msi1_1533, nsh_1362, nsh_1363, \
                         nsh_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1810[k] = f_3 * pc_x[k] * nsh_1362[k];

        t_1811[k] = f_3 * pc_x[k] * nsh_1363[k];

        t_1812[k] = f_3 * pc_x[k] * nsh_1364[k];

        t_1813[k] = pa_y[k] * msi0_1533[k]
                    + f_20 * msh_1149[k]
                    - f_10 * pc_y[k] * msi1_1533[k];
    }

#pragma omp simd aligned(t_1814, t_1815, t_1816, pa_y, pc_y, pc_z, msi0_1535, msi0_1536, \
                         msh_1128, msh_1151, msh_1152, msi1_1535, msi1_1536, \
                         nsh_1359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1814[k] = f_15 * msh_1128[k]
                    + f_3 * pc_z[k] * nsh_1359[k];

        t_1815[k] = pa_y[k] * msi0_1535[k]
                    + f_14 * msh_1151[k]
                    - f_10 * pc_y[k] * msi1_1535[k];

        t_1816[k] = pa_y[k] * msi0_1536[k]
                    + f_13 * msh_1152[k]
                    - f_10 * pc_y[k] * msi1_1536[k];
    }

#pragma omp simd aligned(t_1817, t_1818, t_1819, pa_y, pc_y, msi0_1537, msi0_1539, msh_1153, \
                         msh_1154, msi1_1537, msi1_1539, nsh_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1817[k] = pa_y[k] * msi0_1537[k]
                    + f_12 * msh_1153[k]
                    - f_10 * pc_y[k] * msi1_1537[k];

        t_1818[k] = f_11 * msh_1154[k]
                    + f_3 * pc_y[k] * nsh_1364[k];

        t_1819[k] = pa_y[k] * msi0_1539[k]
                    - f_10 * pc_y[k] * msi1_1539[k];
    }

#pragma omp simd aligned(t_1820, t_1821, t_1822, t_1823, t_1824, pc_x, pc_y, nsg0_975, \
                         nsg0_977, nsg0_978, nsg1_975, nsg1_977, nsg1_978, nsh_1365, nsh_1367, \
                         nsh_1368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1820[k] = f_1 * nsg0_975[k]
                    - f_2 * nsg1_975[k]
                    + f_3 * pc_x[k] * nsh_1365[k];

        t_1821[k] = f_3 * pc_y[k] * nsh_1365[k];

        t_1822[k] = f_16 * nsg0_977[k]
                    - f_17 * nsg1_977[k]
                    + f_3 * pc_x[k] * nsh_1367[k];

        t_1823[k] = f_8 * nsg0_978[k]
                    - f_9 * nsg1_978[k]
                    + f_3 * pc_x[k] * nsh_1368[k];

        t_1824[k] = f_3 * pc_y[k] * nsh_1367[k];
    }

#pragma omp simd aligned(t_1825, t_1826, t_1827, t_1828, pc_x, pc_y, nsg0_980, nsg0_981, \
                         nsg0_982, nsg1_980, nsg1_981, nsg1_982, nsh_1370, nsh_1371, \
                         nsh_1372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1825[k] = f_8 * nsg0_980[k]
                    - f_9 * nsg1_980[k]
                    + f_3 * pc_x[k] * nsh_1370[k];

        t_1826[k] = f_6 * nsg0_981[k]
                    - f_7 * nsg1_981[k]
                    + f_3 * pc_x[k] * nsh_1371[k];

        t_1827[k] = f_6 * nsg0_982[k]
                    - f_7 * nsg1_982[k]
                    + f_3 * pc_x[k] * nsh_1372[k];

        t_1828[k] = f_3 * pc_y[k] * nsh_1370[k];
    }

#pragma omp simd aligned(t_1829, t_1830, t_1831, pc_x, nsg0_984, nsg0_985, nsg0_986, nsg1_984, \
                         nsg1_985, nsg1_986, nsh_1374, nsh_1375, \
                         nsh_1376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1829[k] = f_6 * nsg0_984[k]
                    - f_7 * nsg1_984[k]
                    + f_3 * pc_x[k] * nsh_1374[k];

        t_1830[k] = f_4 * nsg0_985[k]
                    - f_5 * nsg1_985[k]
                    + f_3 * pc_x[k] * nsh_1375[k];

        t_1831[k] = f_4 * nsg0_986[k]
                    - f_5 * nsg1_986[k]
                    + f_3 * pc_x[k] * nsh_1376[k];
    }

#pragma omp simd aligned(t_1832, t_1833, t_1834, t_1835, t_1836, pc_x, pc_y, nsg0_987, \
                         nsg0_989, nsg1_987, nsg1_989, nsh_1374, nsh_1377, nsh_1379, nsh_1380, \
                         nsh_1381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1832[k] = f_4 * nsg0_987[k]
                    - f_5 * nsg1_987[k]
                    + f_3 * pc_x[k] * nsh_1377[k];

        t_1833[k] = f_3 * pc_y[k] * nsh_1374[k];

        t_1834[k] = f_4 * nsg0_989[k]
                    - f_5 * nsg1_989[k]
                    + f_3 * pc_x[k] * nsh_1379[k];

        t_1835[k] = f_3 * pc_x[k] * nsh_1380[k];

        t_1836[k] = f_3 * pc_x[k] * nsh_1381[k];
    }

#pragma omp simd aligned(t_1837, t_1838, t_1839, t_1840, t_1841, pc_x, pc_y, nsg0_985, \
                         nsg1_985, nsh_1380, nsh_1382, nsh_1383, nsh_1384, \
                         nsh_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1837[k] = f_3 * pc_x[k] * nsh_1382[k];

        t_1838[k] = f_3 * pc_x[k] * nsh_1383[k];

        t_1839[k] = f_3 * pc_x[k] * nsh_1384[k];

        t_1840[k] = f_3 * pc_x[k] * nsh_1385[k];

        t_1841[k] = f_1 * nsg0_985[k]
                    - f_2 * nsg1_985[k]
                    + f_3 * pc_y[k] * nsh_1380[k];
    }

#pragma omp simd aligned(t_1842, t_1843, t_1844, pc_y, nsg0_986, nsg0_987, nsg0_988, nsg1_986, \
                         nsg1_987, nsg1_988, nsh_1381, nsh_1382, \
                         nsh_1383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1842[k] = f_16 * nsg0_986[k]
                    - f_17 * nsg1_986[k]
                    + f_3 * pc_y[k] * nsh_1381[k];

        t_1843[k] = f_8 * nsg0_987[k]
                    - f_9 * nsg1_987[k]
                    + f_3 * pc_y[k] * nsh_1382[k];

        t_1844[k] = f_6 * nsg0_988[k]
                    - f_7 * nsg1_988[k]
                    + f_3 * pc_y[k] * nsh_1383[k];
    }

#pragma omp simd aligned(t_1845, t_1846, t_1847, pc_y, pc_z, msh_1154, nsg0_989, nsg1_989, \
                         nsh_1384, nsh_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1845[k] = f_4 * nsg0_989[k]
                    - f_5 * nsg1_989[k]
                    + f_3 * pc_y[k] * nsh_1384[k];

        t_1846[k] = f_3 * pc_y[k] * nsh_1385[k];

        t_1847[k] = f_0 * msh_1154[k]
                    + f_1 * nsg0_989[k]
                    - f_2 * nsg1_989[k]
                    + f_3 * pc_z[k] * nsh_1385[k];
    }
}

auto
compute_prim_nsi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t msi0, const size_t msh,
                                                   const size_t msi1, const size_t nsg0,
                                                   const size_t nsg1, const size_t nsh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_nsi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, msi0, msh,
                                                              msi1, nsg0, nsg1, nsh, ncols,
                                                              gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, msi0,
                                                               msh, msi1, nsg0, nsg1, nsh,
                                                               ncols, gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece11(buffer, target, pa, pc, msi0,
                                                               msh, msi1, nsg0, nsg1, nsh,
                                                               ncols, gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece12(buffer, target, pa, pc, msi0,
                                                               msh, msi1, nsh, ncols, gamma, p,
                                                               q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece13(buffer, target, pa, pc, msi0,
                                                               msh, msi1, nsg0, nsg1, nsh,
                                                               ncols, gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece14(buffer, target, pc, msh, nsg0,
                                                               nsg1, nsh, ncols, gamma, p, q);

    compute_prim_nsi_three_center_electron_repulsion_0_piece15(buffer, target, pa, pc, msi0,
                                                               msh, msi1, nsg0, nsg1, nsh,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
