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


#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.0 / q;
    const auto f_10 = 4.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *sng0_0 = buffer.data(sng0 + 0);
    const auto *sng0_3 = buffer.data(sng0 + 3);
    const auto *sng0_5 = buffer.data(sng0 + 5);
    const auto *sng0_10 = buffer.data(sng0 + 10);
    const auto *sng0_14 = buffer.data(sng0 + 14);
    const auto *sng0_18 = buffer.data(sng0 + 18);
    const auto *sng0_25 = buffer.data(sng0 + 25);
    const auto *sng0_30 = buffer.data(sng0 + 30);
    const auto *sng0_35 = buffer.data(sng0 + 35);
    const auto *sng0_44 = buffer.data(sng0 + 44);
    const auto *sng0_45 = buffer.data(sng0 + 45);
    const auto *sng0_48 = buffer.data(sng0 + 48);
    const auto *sng0_55 = buffer.data(sng0 + 55);
    const auto *sng0_75 = buffer.data(sng0 + 75);
    const auto *sng0_78 = buffer.data(sng0 + 78);

    const auto *snf_0 = buffer.data(snf + 0);
    const auto *snf_1 = buffer.data(snf + 1);
    const auto *snf_2 = buffer.data(snf + 2);
    const auto *snf_3 = buffer.data(snf + 3);
    const auto *snf_5 = buffer.data(snf + 5);
    const auto *snf_6 = buffer.data(snf + 6);
    const auto *snf_7 = buffer.data(snf + 7);
    const auto *snf_8 = buffer.data(snf + 8);
    const auto *snf_9 = buffer.data(snf + 9);
    const auto *snf_10 = buffer.data(snf + 10);
    const auto *snf_12 = buffer.data(snf + 12);
    const auto *snf_16 = buffer.data(snf + 16);
    const auto *snf_17 = buffer.data(snf + 17);
    const auto *snf_18 = buffer.data(snf + 18);
    const auto *snf_19 = buffer.data(snf + 19);
    const auto *snf_20 = buffer.data(snf + 20);
    const auto *snf_22 = buffer.data(snf + 22);
    const auto *snf_26 = buffer.data(snf + 26);
    const auto *snf_27 = buffer.data(snf + 27);
    const auto *snf_28 = buffer.data(snf + 28);
    const auto *snf_29 = buffer.data(snf + 29);
    const auto *snf_30 = buffer.data(snf + 30);
    const auto *snf_32 = buffer.data(snf + 32);
    const auto *snf_33 = buffer.data(snf + 33);
    const auto *snf_35 = buffer.data(snf + 35);
    const auto *snf_36 = buffer.data(snf + 36);
    const auto *snf_37 = buffer.data(snf + 37);
    const auto *snf_38 = buffer.data(snf + 38);
    const auto *snf_39 = buffer.data(snf + 39);
    const auto *snf_40 = buffer.data(snf + 40);
    const auto *snf_42 = buffer.data(snf + 42);
    const auto *snf_46 = buffer.data(snf + 46);
    const auto *snf_47 = buffer.data(snf + 47);
    const auto *snf_48 = buffer.data(snf + 48);
    const auto *snf_49 = buffer.data(snf + 49);
    const auto *snf_50 = buffer.data(snf + 50);
    const auto *snf_51 = buffer.data(snf + 51);
    const auto *snf_52 = buffer.data(snf + 52);
    const auto *snf_53 = buffer.data(snf + 53);
    const auto *snf_55 = buffer.data(snf + 55);
    const auto *snf_56 = buffer.data(snf + 56);
    const auto *snf_57 = buffer.data(snf + 57);
    const auto *snf_58 = buffer.data(snf + 58);
    const auto *snf_59 = buffer.data(snf + 59);
    const auto *snf_60 = buffer.data(snf + 60);
    const auto *snf_63 = buffer.data(snf + 63);
    const auto *snf_65 = buffer.data(snf + 65);
    const auto *snf_66 = buffer.data(snf + 66);
    const auto *snf_67 = buffer.data(snf + 67);
    const auto *snf_68 = buffer.data(snf + 68);
    const auto *snf_69 = buffer.data(snf + 69);
    const auto *snf_75 = buffer.data(snf + 75);
    const auto *snf_76 = buffer.data(snf + 76);
    const auto *snf_77 = buffer.data(snf + 77);
    const auto *snf_78 = buffer.data(snf + 78);
    const auto *snf_79 = buffer.data(snf + 79);

    const auto *sng1_0 = buffer.data(sng1 + 0);
    const auto *sng1_3 = buffer.data(sng1 + 3);
    const auto *sng1_5 = buffer.data(sng1 + 5);
    const auto *sng1_10 = buffer.data(sng1 + 10);
    const auto *sng1_14 = buffer.data(sng1 + 14);
    const auto *sng1_18 = buffer.data(sng1 + 18);
    const auto *sng1_25 = buffer.data(sng1 + 25);
    const auto *sng1_30 = buffer.data(sng1 + 30);
    const auto *sng1_35 = buffer.data(sng1 + 35);
    const auto *sng1_44 = buffer.data(sng1 + 44);
    const auto *sng1_45 = buffer.data(sng1 + 45);
    const auto *sng1_48 = buffer.data(sng1 + 48);
    const auto *sng1_55 = buffer.data(sng1 + 55);
    const auto *sng1_75 = buffer.data(sng1 + 75);
    const auto *sng1_78 = buffer.data(sng1 + 78);

    const auto *sod0_0 = buffer.data(sod0 + 0);
    const auto *sod0_3 = buffer.data(sod0 + 3);
    const auto *sod0_5 = buffer.data(sod0 + 5);
    const auto *sod0_9 = buffer.data(sod0 + 9);
    const auto *sod0_11 = buffer.data(sod0 + 11);
    const auto *sod0_17 = buffer.data(sod0 + 17);
    const auto *sod0_18 = buffer.data(sod0 + 18);
    const auto *sod0_21 = buffer.data(sod0 + 21);
    const auto *sod0_23 = buffer.data(sod0 + 23);
    const auto *sod0_29 = buffer.data(sod0 + 29);
    const auto *sod0_30 = buffer.data(sod0 + 30);
    const auto *sod0_33 = buffer.data(sod0 + 33);
    const auto *sod0_35 = buffer.data(sod0 + 35);
    const auto *sod0_36 = buffer.data(sod0 + 36);
    const auto *sod0_39 = buffer.data(sod0 + 39);
    const auto *sod0_41 = buffer.data(sod0 + 41);
    const auto *sod0_47 = buffer.data(sod0 + 47);

    const auto *sod1_0 = buffer.data(sod1 + 0);
    const auto *sod1_3 = buffer.data(sod1 + 3);
    const auto *sod1_5 = buffer.data(sod1 + 5);
    const auto *sod1_9 = buffer.data(sod1 + 9);
    const auto *sod1_11 = buffer.data(sod1 + 11);
    const auto *sod1_17 = buffer.data(sod1 + 17);
    const auto *sod1_18 = buffer.data(sod1 + 18);
    const auto *sod1_21 = buffer.data(sod1 + 21);
    const auto *sod1_23 = buffer.data(sod1 + 23);
    const auto *sod1_29 = buffer.data(sod1 + 29);
    const auto *sod1_30 = buffer.data(sod1 + 30);
    const auto *sod1_33 = buffer.data(sod1 + 33);
    const auto *sod1_35 = buffer.data(sod1 + 35);
    const auto *sod1_36 = buffer.data(sod1 + 36);
    const auto *sod1_39 = buffer.data(sod1 + 39);
    const auto *sod1_41 = buffer.data(sod1 + 41);
    const auto *sod1_47 = buffer.data(sod1 + 47);

    const auto *sof_0 = buffer.data(sof + 0);
    const auto *sof_2 = buffer.data(sof + 2);
    const auto *sof_3 = buffer.data(sof + 3);
    const auto *sof_5 = buffer.data(sof + 5);
    const auto *sof_6 = buffer.data(sof + 6);
    const auto *sof_7 = buffer.data(sof + 7);
    const auto *sof_8 = buffer.data(sof + 8);
    const auto *sof_9 = buffer.data(sof + 9);
    const auto *sof_10 = buffer.data(sof + 10);
    const auto *sof_12 = buffer.data(sof + 12);
    const auto *sof_16 = buffer.data(sof + 16);
    const auto *sof_17 = buffer.data(sof + 17);
    const auto *sof_18 = buffer.data(sof + 18);
    const auto *sof_19 = buffer.data(sof + 19);
    const auto *sof_20 = buffer.data(sof + 20);
    const auto *sof_22 = buffer.data(sof + 22);
    const auto *sof_26 = buffer.data(sof + 26);
    const auto *sof_27 = buffer.data(sof + 27);
    const auto *sof_28 = buffer.data(sof + 28);
    const auto *sof_29 = buffer.data(sof + 29);
    const auto *sof_30 = buffer.data(sof + 30);
    const auto *sof_32 = buffer.data(sof + 32);
    const auto *sof_33 = buffer.data(sof + 33);
    const auto *sof_35 = buffer.data(sof + 35);
    const auto *sof_36 = buffer.data(sof + 36);
    const auto *sof_37 = buffer.data(sof + 37);
    const auto *sof_38 = buffer.data(sof + 38);
    const auto *sof_39 = buffer.data(sof + 39);
    const auto *sof_40 = buffer.data(sof + 40);
    const auto *sof_42 = buffer.data(sof + 42);
    const auto *sof_46 = buffer.data(sof + 46);
    const auto *sof_47 = buffer.data(sof + 47);
    const auto *sof_48 = buffer.data(sof + 48);
    const auto *sof_49 = buffer.data(sof + 49);
    const auto *sof_50 = buffer.data(sof + 50);
    const auto *sof_52 = buffer.data(sof + 52);
    const auto *sof_53 = buffer.data(sof + 53);
    const auto *sof_55 = buffer.data(sof + 55);
    const auto *sof_56 = buffer.data(sof + 56);
    const auto *sof_57 = buffer.data(sof + 57);
    const auto *sof_58 = buffer.data(sof + 58);
    const auto *sof_59 = buffer.data(sof + 59);
    const auto *sof_60 = buffer.data(sof + 60);
    const auto *sof_62 = buffer.data(sof + 62);
    const auto *sof_63 = buffer.data(sof + 63);
    const auto *sof_65 = buffer.data(sof + 65);
    const auto *sof_66 = buffer.data(sof + 66);
    const auto *sof_67 = buffer.data(sof + 67);
    const auto *sof_68 = buffer.data(sof + 68);
    const auto *sof_69 = buffer.data(sof + 69);
    const auto *sof_70 = buffer.data(sof + 70);
    const auto *sof_72 = buffer.data(sof + 72);
    const auto *sof_75 = buffer.data(sof + 75);
    const auto *sof_76 = buffer.data(sof + 76);
    const auto *sof_77 = buffer.data(sof + 77);
    const auto *sof_78 = buffer.data(sof + 78);
    const auto *sof_79 = buffer.data(sof + 79);
    const auto *sof_80 = buffer.data(sof + 80);
    const auto *sof_82 = buffer.data(sof + 82);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, snf_0, snf_3, sod0_0, sod0_3, \
                         sod1_0, sod1_3, sof_0, sof_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * snf_0[k]
                 + f_1 * sod0_0[k]
                 - f_2 * sod1_0[k]
                 + f_3 * pc_x[k] * sof_0[k];

        t_1[k] = f_3 * pc_y[k] * sof_0[k];

        t_2[k] = f_3 * pc_z[k] * sof_0[k];

        t_3[k] = f_0 * snf_3[k]
                 + f_4 * sod0_3[k]
                 - f_5 * sod1_3[k]
                 + f_3 * pc_x[k] * sof_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, snf_5, snf_6, snf_7, sod0_5, sod1_5, \
                         sof_2, sof_5, sof_6, sof_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sof_2[k];

        t_5[k] = f_0 * snf_5[k]
                 + f_4 * sod0_5[k]
                 - f_5 * sod1_5[k]
                 + f_3 * pc_x[k] * sof_5[k];

        t_6[k] = f_0 * snf_6[k]
                 + f_3 * pc_x[k] * sof_6[k];

        t_7[k] = f_0 * snf_7[k]
                 + f_3 * pc_x[k] * sof_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, snf_8, snf_9, sod0_3, sod1_3, \
                         sof_6, sof_8, sof_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * snf_8[k]
                 + f_3 * pc_x[k] * sof_8[k];

        t_9[k] = f_0 * snf_9[k]
                 + f_3 * pc_x[k] * sof_9[k];

        t_10[k] = f_1 * sod0_3[k]
                  - f_2 * sod1_3[k]
                  + f_3 * pc_y[k] * sof_6[k];

        t_11[k] = f_3 * pc_z[k] * sof_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, sng0_0, snf_0, \
                         sng1_0, sod0_5, sod1_5, sof_8, sof_9, sof_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * sod0_5[k]
                  - f_5 * sod1_5[k]
                  + f_3 * pc_y[k] * sof_8[k];

        t_13[k] = f_3 * pc_y[k] * sof_9[k];

        t_14[k] = f_1 * sod0_5[k]
                  - f_2 * sod1_5[k]
                  + f_3 * pc_z[k] * sof_9[k];

        t_15[k] = pb_y[k] * sng0_0[k]
                  - f_6 * pc_y[k] * sng1_0[k];

        t_16[k] = f_7 * snf_0[k]
                  + f_3 * pc_y[k] * sof_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, sng0_3, sng0_5, snf_1, \
                         snf_2, sng1_3, sng1_5, sof_10, sof_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * sof_10[k];

        t_18[k] = pb_y[k] * sng0_3[k]
                  + f_8 * snf_1[k]
                  - f_6 * pc_y[k] * sng1_3[k];

        t_19[k] = f_7 * snf_2[k]
                  + f_3 * pc_y[k] * sof_12[k];

        t_20[k] = pb_y[k] * sng0_5[k]
                  - f_6 * pc_y[k] * sng1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, snf_16, snf_17, snf_18, snf_19, sof_16, \
                         sof_17, sof_18, sof_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * snf_16[k]
                  + f_3 * pc_x[k] * sof_16[k];

        t_22[k] = f_9 * snf_17[k]
                  + f_3 * pc_x[k] * sof_17[k];

        t_23[k] = f_9 * snf_18[k]
                  + f_3 * pc_x[k] * sof_18[k];

        t_24[k] = f_9 * snf_19[k]
                  + f_3 * pc_x[k] * sof_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, snf_6, snf_8, snf_9, sod0_9, \
                         sod0_11, sod1_9, sod1_11, sof_16, sof_18, \
                         sof_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * snf_6[k]
                  + f_1 * sod0_9[k]
                  - f_2 * sod1_9[k]
                  + f_3 * pc_y[k] * sof_16[k];

        t_26[k] = f_3 * pc_z[k] * sof_16[k];

        t_27[k] = f_7 * snf_8[k]
                  + f_4 * sod0_11[k]
                  - f_5 * sod1_11[k]
                  + f_3 * pc_y[k] * sof_18[k];

        t_28[k] = f_7 * snf_9[k]
                  + f_3 * pc_y[k] * sof_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, sng0_0, sng0_14, \
                         snf_0, sng1_0, sng1_14, sof_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * sng0_14[k]
                  - f_6 * pc_y[k] * sng1_14[k];

        t_30[k] = pb_z[k] * sng0_0[k]
                  - f_6 * pc_z[k] * sng1_0[k];

        t_31[k] = f_3 * pc_y[k] * sof_20[k];

        t_32[k] = f_7 * snf_0[k]
                  + f_3 * pc_z[k] * sof_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, sng0_3, sng0_5, \
                         snf_2, snf_26, sng1_3, sng1_5, sof_22, \
                         sof_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * sng0_3[k]
                  - f_6 * pc_z[k] * sng1_3[k];

        t_34[k] = f_3 * pc_y[k] * sof_22[k];

        t_35[k] = pb_z[k] * sng0_5[k]
                  + f_8 * snf_2[k]
                  - f_6 * pc_z[k] * sng1_5[k];

        t_36[k] = f_9 * snf_26[k]
                  + f_3 * pc_x[k] * sof_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, sng0_10, snf_27, snf_28, \
                         snf_29, sng1_10, sof_27, sof_28, sof_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * snf_27[k]
                  + f_3 * pc_x[k] * sof_27[k];

        t_38[k] = f_9 * snf_28[k]
                  + f_3 * pc_x[k] * sof_28[k];

        t_39[k] = f_9 * snf_29[k]
                  + f_3 * pc_x[k] * sof_29[k];

        t_40[k] = pb_z[k] * sng0_10[k]
                  - f_6 * pc_z[k] * sng1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, snf_6, snf_9, sod0_17, sod1_17, \
                         sof_26, sof_28, sof_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * snf_6[k]
                  + f_3 * pc_z[k] * sof_26[k];

        t_42[k] = f_4 * sod0_17[k]
                  - f_5 * sod1_17[k]
                  + f_3 * pc_y[k] * sof_28[k];

        t_43[k] = f_3 * pc_y[k] * sof_29[k];

        t_44[k] = f_7 * snf_9[k]
                  + f_1 * sod0_17[k]
                  - f_2 * sod1_17[k]
                  + f_3 * pc_z[k] * sof_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, snf_10, snf_30, snf_33, \
                         sod0_18, sod0_21, sod1_18, sod1_21, sof_30, \
                         sof_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_10 * snf_30[k]
                  + f_1 * sod0_18[k]
                  - f_2 * sod1_18[k]
                  + f_3 * pc_x[k] * sof_30[k];

        t_46[k] = f_8 * snf_10[k]
                  + f_3 * pc_y[k] * sof_30[k];

        t_47[k] = f_3 * pc_z[k] * sof_30[k];

        t_48[k] = f_10 * snf_33[k]
                  + f_4 * sod0_21[k]
                  - f_5 * sod1_21[k]
                  + f_3 * pc_x[k] * sof_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, snf_12, snf_35, snf_36, snf_37, \
                         sod0_23, sod1_23, sof_32, sof_35, sof_36, \
                         sof_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * snf_12[k]
                  + f_3 * pc_y[k] * sof_32[k];

        t_50[k] = f_10 * snf_35[k]
                  + f_4 * sod0_23[k]
                  - f_5 * sod1_23[k]
                  + f_3 * pc_x[k] * sof_35[k];

        t_51[k] = f_10 * snf_36[k]
                  + f_3 * pc_x[k] * sof_36[k];

        t_52[k] = f_10 * snf_37[k]
                  + f_3 * pc_x[k] * sof_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, snf_16, snf_38, snf_39, \
                         sod0_21, sod1_21, sof_36, sof_38, sof_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * snf_38[k]
                  + f_3 * pc_x[k] * sof_38[k];

        t_54[k] = f_10 * snf_39[k]
                  + f_3 * pc_x[k] * sof_39[k];

        t_55[k] = f_8 * snf_16[k]
                  + f_1 * sod0_21[k]
                  - f_2 * sod1_21[k]
                  + f_3 * pc_y[k] * sof_36[k];

        t_56[k] = f_3 * pc_z[k] * sof_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, sng0_30, snf_18, snf_19, \
                         sng1_30, sod0_23, sod1_23, sof_38, sof_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * snf_18[k]
                  + f_4 * sod0_23[k]
                  - f_5 * sod1_23[k]
                  + f_3 * pc_y[k] * sof_38[k];

        t_58[k] = f_8 * snf_19[k]
                  + f_3 * pc_y[k] * sof_39[k];

        t_59[k] = f_1 * sod0_23[k]
                  - f_2 * sod1_23[k]
                  + f_3 * pc_z[k] * sof_39[k];

        t_60[k] = pb_y[k] * sng0_30[k]
                  - f_6 * pc_y[k] * sng1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sng0_18, snf_10, snf_20, \
                         snf_22, sng1_18, sof_40, sof_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * snf_20[k]
                  + f_3 * pc_y[k] * sof_40[k];

        t_62[k] = f_7 * snf_10[k]
                  + f_3 * pc_z[k] * sof_40[k];

        t_63[k] = pb_z[k] * sng0_18[k]
                  - f_6 * pc_z[k] * sng1_18[k];

        t_64[k] = f_7 * snf_22[k]
                  + f_3 * pc_y[k] * sof_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, sng0_35, snf_46, snf_47, \
                         snf_48, sng1_35, sof_46, sof_47, sof_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sng0_35[k]
                  - f_6 * pc_y[k] * sng1_35[k];

        t_66[k] = f_10 * snf_46[k]
                  + f_3 * pc_x[k] * sof_46[k];

        t_67[k] = f_10 * snf_47[k]
                  + f_3 * pc_x[k] * sof_47[k];

        t_68[k] = f_10 * snf_48[k]
                  + f_3 * pc_x[k] * sof_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, pc_x, pc_z, sng0_25, snf_16, snf_49, sng1_25, \
                         sof_46, sof_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * snf_49[k]
                  + f_3 * pc_x[k] * sof_49[k];

        t_70[k] = pb_z[k] * sng0_25[k]
                  - f_6 * pc_z[k] * sng1_25[k];

        t_71[k] = f_7 * snf_16[k]
                  + f_3 * pc_z[k] * sof_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, pc_y, sng0_44, snf_28, snf_29, sng1_44, \
                         sod0_29, sod1_29, sof_48, sof_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * snf_28[k]
                  + f_4 * sod0_29[k]
                  - f_5 * sod1_29[k]
                  + f_3 * pc_y[k] * sof_48[k];

        t_73[k] = f_7 * snf_29[k]
                  + f_3 * pc_y[k] * sof_49[k];

        t_74[k] = pb_y[k] * sng0_44[k]
                  - f_6 * pc_y[k] * sng1_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, snf_20, snf_50, snf_53, \
                         sod0_30, sod0_33, sod1_30, sod1_33, sof_50, \
                         sof_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * snf_50[k]
                  + f_1 * sod0_30[k]
                  - f_2 * sod1_30[k]
                  + f_3 * pc_x[k] * sof_50[k];

        t_76[k] = f_3 * pc_y[k] * sof_50[k];

        t_77[k] = f_8 * snf_20[k]
                  + f_3 * pc_z[k] * sof_50[k];

        t_78[k] = f_10 * snf_53[k]
                  + f_4 * sod0_33[k]
                  - f_5 * sod1_33[k]
                  + f_3 * pc_x[k] * sof_53[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, snf_55, snf_56, snf_57, sod0_35, \
                         sod1_35, sof_52, sof_55, sof_56, sof_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_y[k] * sof_52[k];

        t_80[k] = f_10 * snf_55[k]
                  + f_4 * sod0_35[k]
                  - f_5 * sod1_35[k]
                  + f_3 * pc_x[k] * sof_55[k];

        t_81[k] = f_10 * snf_56[k]
                  + f_3 * pc_x[k] * sof_56[k];

        t_82[k] = f_10 * snf_57[k]
                  + f_3 * pc_x[k] * sof_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, snf_26, snf_58, snf_59, \
                         sod0_33, sod1_33, sof_56, sof_58, sof_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_10 * snf_58[k]
                  + f_3 * pc_x[k] * sof_58[k];

        t_84[k] = f_10 * snf_59[k]
                  + f_3 * pc_x[k] * sof_59[k];

        t_85[k] = f_1 * sod0_33[k]
                  - f_2 * sod1_33[k]
                  + f_3 * pc_y[k] * sof_56[k];

        t_86[k] = f_8 * snf_26[k]
                  + f_3 * pc_z[k] * sof_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, snf_29, snf_60, sod0_35, \
                         sod0_36, sod1_35, sod1_36, sof_58, sof_59, \
                         sof_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * sod0_35[k]
                  - f_5 * sod1_35[k]
                  + f_3 * pc_y[k] * sof_58[k];

        t_88[k] = f_3 * pc_y[k] * sof_59[k];

        t_89[k] = f_8 * snf_29[k]
                  + f_1 * sod0_35[k]
                  - f_2 * sod1_35[k]
                  + f_3 * pc_z[k] * sof_59[k];

        t_90[k] = f_11 * snf_60[k]
                  + f_1 * sod0_36[k]
                  - f_2 * sod1_36[k]
                  + f_3 * pc_x[k] * sof_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, snf_30, snf_32, snf_63, \
                         sod0_39, sod1_39, sof_60, sof_62, sof_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_12 * snf_30[k]
                  + f_3 * pc_y[k] * sof_60[k];

        t_92[k] = f_3 * pc_z[k] * sof_60[k];

        t_93[k] = f_11 * snf_63[k]
                  + f_4 * sod0_39[k]
                  - f_5 * sod1_39[k]
                  + f_3 * pc_x[k] * sof_63[k];

        t_94[k] = f_12 * snf_32[k]
                  + f_3 * pc_y[k] * sof_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, snf_65, snf_66, snf_67, snf_68, \
                         sod0_41, sod1_41, sof_65, sof_66, sof_67, \
                         sof_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * snf_65[k]
                  + f_4 * sod0_41[k]
                  - f_5 * sod1_41[k]
                  + f_3 * pc_x[k] * sof_65[k];

        t_96[k] = f_11 * snf_66[k]
                  + f_3 * pc_x[k] * sof_66[k];

        t_97[k] = f_11 * snf_67[k]
                  + f_3 * pc_x[k] * sof_67[k];

        t_98[k] = f_11 * snf_68[k]
                  + f_3 * pc_x[k] * sof_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_x, pc_y, pc_z, snf_36, snf_69, sod0_39, \
                         sod1_39, sof_66, sof_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_11 * snf_69[k]
                  + f_3 * pc_x[k] * sof_69[k];

        t_100[k] = f_12 * snf_36[k]
                   + f_1 * sod0_39[k]
                   - f_2 * sod1_39[k]
                   + f_3 * pc_y[k] * sof_66[k];

        t_101[k] = f_3 * pc_z[k] * sof_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, sng0_45, snf_38, \
                         snf_39, sng1_45, sod0_41, sod1_41, sof_68, \
                         sof_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * snf_38[k]
                   + f_4 * sod0_41[k]
                   - f_5 * sod1_41[k]
                   + f_3 * pc_y[k] * sof_68[k];

        t_103[k] = f_12 * snf_39[k]
                   + f_3 * pc_y[k] * sof_69[k];

        t_104[k] = f_1 * sod0_41[k]
                   - f_2 * sod1_41[k]
                   + f_3 * pc_z[k] * sof_69[k];

        t_105[k] = pb_z[k] * sng0_45[k]
                   - f_6 * pc_z[k] * sng1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, sng0_48, snf_30, \
                         snf_40, snf_42, sng1_48, sof_70, sof_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * snf_40[k]
                   + f_3 * pc_y[k] * sof_70[k];

        t_107[k] = f_7 * snf_30[k]
                   + f_3 * pc_z[k] * sof_70[k];

        t_108[k] = pb_z[k] * sng0_48[k]
                   - f_6 * pc_z[k] * sng1_48[k];

        t_109[k] = f_8 * snf_42[k]
                   + f_3 * pc_y[k] * sof_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, snf_75, snf_76, snf_77, snf_78, \
                         sod0_47, sod1_47, sof_75, sof_76, sof_77, \
                         sof_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_11 * snf_75[k]
                   + f_4 * sod0_47[k]
                   - f_5 * sod1_47[k]
                   + f_3 * pc_x[k] * sof_75[k];

        t_111[k] = f_11 * snf_76[k]
                   + f_3 * pc_x[k] * sof_76[k];

        t_112[k] = f_11 * snf_77[k]
                   + f_3 * pc_x[k] * sof_77[k];

        t_113[k] = f_11 * snf_78[k]
                   + f_3 * pc_x[k] * sof_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_z, pc_x, pc_z, sng0_55, snf_36, snf_79, \
                         sng1_55, sof_76, sof_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_11 * snf_79[k]
                   + f_3 * pc_x[k] * sof_79[k];

        t_115[k] = pb_z[k] * sng0_55[k]
                   - f_6 * pc_z[k] * sng1_55[k];

        t_116[k] = f_7 * snf_36[k]
                   + f_3 * pc_z[k] * sof_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pc_y, pc_z, sng0_75, snf_39, \
                         snf_48, snf_49, sng1_75, sod0_47, sod1_47, sof_78, \
                         sof_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * snf_48[k]
                   + f_4 * sod0_47[k]
                   - f_5 * sod1_47[k]
                   + f_3 * pc_y[k] * sof_78[k];

        t_118[k] = f_8 * snf_49[k]
                   + f_3 * pc_y[k] * sof_79[k];

        t_119[k] = f_7 * snf_39[k]
                   + f_1 * sod0_47[k]
                   - f_2 * sod1_47[k]
                   + f_3 * pc_z[k] * sof_79[k];

        t_120[k] = pb_y[k] * sng0_75[k]
                   - f_6 * pc_y[k] * sng1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pc_y, pc_z, sng0_78, snf_40, \
                         snf_50, snf_51, snf_52, sng1_78, sof_80, \
                         sof_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * snf_50[k]
                   + f_3 * pc_y[k] * sof_80[k];

        t_122[k] = f_8 * snf_40[k]
                   + f_3 * pc_z[k] * sof_80[k];

        t_123[k] = pb_y[k] * sng0_78[k]
                   + f_8 * snf_51[k]
                   - f_6 * pc_y[k] * sng1_78[k];

        t_124[k] = f_7 * snf_52[k]
                   + f_3 * pc_y[k] * sof_82[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
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
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_80 = buffer.data(sng0 + 80);
    const auto *sng0_89 = buffer.data(sng0 + 89);
    const auto *sng0_90 = buffer.data(sng0 + 90);
    const auto *sng0_93 = buffer.data(sng0 + 93);
    const auto *sng0_100 = buffer.data(sng0 + 100);
    const auto *sng0_135 = buffer.data(sng0 + 135);
    const auto *sng0_138 = buffer.data(sng0 + 138);
    const auto *sng0_140 = buffer.data(sng0 + 140);
    const auto *sng0_149 = buffer.data(sng0 + 149);
    const auto *sng0_150 = buffer.data(sng0 + 150);
    const auto *sng0_153 = buffer.data(sng0 + 153);

    const auto *snf_46 = buffer.data(snf + 46);
    const auto *snf_50 = buffer.data(snf + 50);
    const auto *snf_56 = buffer.data(snf + 56);
    const auto *snf_58 = buffer.data(snf + 58);
    const auto *snf_59 = buffer.data(snf + 59);
    const auto *snf_60 = buffer.data(snf + 60);
    const auto *snf_62 = buffer.data(snf + 62);
    const auto *snf_66 = buffer.data(snf + 66);
    const auto *snf_68 = buffer.data(snf + 68);
    const auto *snf_69 = buffer.data(snf + 69);
    const auto *snf_70 = buffer.data(snf + 70);
    const auto *snf_72 = buffer.data(snf + 72);
    const auto *snf_76 = buffer.data(snf + 76);
    const auto *snf_78 = buffer.data(snf + 78);
    const auto *snf_79 = buffer.data(snf + 79);
    const auto *snf_80 = buffer.data(snf + 80);
    const auto *snf_82 = buffer.data(snf + 82);
    const auto *snf_86 = buffer.data(snf + 86);
    const auto *snf_87 = buffer.data(snf + 87);
    const auto *snf_88 = buffer.data(snf + 88);
    const auto *snf_89 = buffer.data(snf + 89);
    const auto *snf_90 = buffer.data(snf + 90);
    const auto *snf_91 = buffer.data(snf + 91);
    const auto *snf_92 = buffer.data(snf + 92);
    const auto *snf_93 = buffer.data(snf + 93);
    const auto *snf_95 = buffer.data(snf + 95);
    const auto *snf_96 = buffer.data(snf + 96);
    const auto *snf_97 = buffer.data(snf + 97);
    const auto *snf_98 = buffer.data(snf + 98);
    const auto *snf_99 = buffer.data(snf + 99);
    const auto *snf_100 = buffer.data(snf + 100);
    const auto *snf_102 = buffer.data(snf + 102);
    const auto *snf_103 = buffer.data(snf + 103);
    const auto *snf_105 = buffer.data(snf + 105);
    const auto *snf_106 = buffer.data(snf + 106);
    const auto *snf_107 = buffer.data(snf + 107);
    const auto *snf_108 = buffer.data(snf + 108);
    const auto *snf_109 = buffer.data(snf + 109);
    const auto *snf_110 = buffer.data(snf + 110);
    const auto *snf_112 = buffer.data(snf + 112);
    const auto *snf_115 = buffer.data(snf + 115);
    const auto *snf_116 = buffer.data(snf + 116);
    const auto *snf_117 = buffer.data(snf + 117);
    const auto *snf_118 = buffer.data(snf + 118);
    const auto *snf_119 = buffer.data(snf + 119);
    const auto *snf_120 = buffer.data(snf + 120);
    const auto *snf_123 = buffer.data(snf + 123);
    const auto *snf_125 = buffer.data(snf + 125);
    const auto *snf_126 = buffer.data(snf + 126);
    const auto *snf_127 = buffer.data(snf + 127);
    const auto *snf_128 = buffer.data(snf + 128);
    const auto *snf_129 = buffer.data(snf + 129);
    const auto *snf_136 = buffer.data(snf + 136);
    const auto *snf_137 = buffer.data(snf + 137);
    const auto *snf_138 = buffer.data(snf + 138);
    const auto *snf_139 = buffer.data(snf + 139);
    const auto *snf_140 = buffer.data(snf + 140);
    const auto *snf_143 = buffer.data(snf + 143);
    const auto *snf_145 = buffer.data(snf + 145);
    const auto *snf_146 = buffer.data(snf + 146);
    const auto *snf_147 = buffer.data(snf + 147);
    const auto *snf_148 = buffer.data(snf + 148);
    const auto *snf_149 = buffer.data(snf + 149);
    const auto *snf_150 = buffer.data(snf + 150);
    const auto *snf_153 = buffer.data(snf + 153);
    const auto *snf_155 = buffer.data(snf + 155);
    const auto *snf_156 = buffer.data(snf + 156);
    const auto *snf_157 = buffer.data(snf + 157);
    const auto *snf_158 = buffer.data(snf + 158);
    const auto *snf_159 = buffer.data(snf + 159);

    const auto *sng1_80 = buffer.data(sng1 + 80);
    const auto *sng1_89 = buffer.data(sng1 + 89);
    const auto *sng1_90 = buffer.data(sng1 + 90);
    const auto *sng1_93 = buffer.data(sng1 + 93);
    const auto *sng1_100 = buffer.data(sng1 + 100);
    const auto *sng1_135 = buffer.data(sng1 + 135);
    const auto *sng1_138 = buffer.data(sng1 + 138);
    const auto *sng1_140 = buffer.data(sng1 + 140);
    const auto *sng1_149 = buffer.data(sng1 + 149);
    const auto *sng1_150 = buffer.data(sng1 + 150);
    const auto *sng1_153 = buffer.data(sng1 + 153);

    const auto *sod0_51 = buffer.data(sod0 + 51);
    const auto *sod0_53 = buffer.data(sod0 + 53);
    const auto *sod0_54 = buffer.data(sod0 + 54);
    const auto *sod0_57 = buffer.data(sod0 + 57);
    const auto *sod0_59 = buffer.data(sod0 + 59);
    const auto *sod0_60 = buffer.data(sod0 + 60);
    const auto *sod0_63 = buffer.data(sod0 + 63);
    const auto *sod0_65 = buffer.data(sod0 + 65);
    const auto *sod0_71 = buffer.data(sod0 + 71);
    const auto *sod0_72 = buffer.data(sod0 + 72);
    const auto *sod0_75 = buffer.data(sod0 + 75);
    const auto *sod0_77 = buffer.data(sod0 + 77);
    const auto *sod0_81 = buffer.data(sod0 + 81);
    const auto *sod0_83 = buffer.data(sod0 + 83);
    const auto *sod0_84 = buffer.data(sod0 + 84);
    const auto *sod0_87 = buffer.data(sod0 + 87);
    const auto *sod0_89 = buffer.data(sod0 + 89);
    const auto *sod0_90 = buffer.data(sod0 + 90);
    const auto *sod0_93 = buffer.data(sod0 + 93);
    const auto *sod0_95 = buffer.data(sod0 + 95);

    const auto *sod1_51 = buffer.data(sod1 + 51);
    const auto *sod1_53 = buffer.data(sod1 + 53);
    const auto *sod1_54 = buffer.data(sod1 + 54);
    const auto *sod1_57 = buffer.data(sod1 + 57);
    const auto *sod1_59 = buffer.data(sod1 + 59);
    const auto *sod1_60 = buffer.data(sod1 + 60);
    const auto *sod1_63 = buffer.data(sod1 + 63);
    const auto *sod1_65 = buffer.data(sod1 + 65);
    const auto *sod1_71 = buffer.data(sod1 + 71);
    const auto *sod1_72 = buffer.data(sod1 + 72);
    const auto *sod1_75 = buffer.data(sod1 + 75);
    const auto *sod1_77 = buffer.data(sod1 + 77);
    const auto *sod1_81 = buffer.data(sod1 + 81);
    const auto *sod1_83 = buffer.data(sod1 + 83);
    const auto *sod1_84 = buffer.data(sod1 + 84);
    const auto *sod1_87 = buffer.data(sod1 + 87);
    const auto *sod1_89 = buffer.data(sod1 + 89);
    const auto *sod1_90 = buffer.data(sod1 + 90);
    const auto *sod1_93 = buffer.data(sod1 + 93);
    const auto *sod1_95 = buffer.data(sod1 + 95);

    const auto *sof_86 = buffer.data(sof + 86);
    const auto *sof_87 = buffer.data(sof + 87);
    const auto *sof_88 = buffer.data(sof + 88);
    const auto *sof_89 = buffer.data(sof + 89);
    const auto *sof_90 = buffer.data(sof + 90);
    const auto *sof_92 = buffer.data(sof + 92);
    const auto *sof_93 = buffer.data(sof + 93);
    const auto *sof_95 = buffer.data(sof + 95);
    const auto *sof_96 = buffer.data(sof + 96);
    const auto *sof_97 = buffer.data(sof + 97);
    const auto *sof_98 = buffer.data(sof + 98);
    const auto *sof_99 = buffer.data(sof + 99);
    const auto *sof_100 = buffer.data(sof + 100);
    const auto *sof_102 = buffer.data(sof + 102);
    const auto *sof_103 = buffer.data(sof + 103);
    const auto *sof_105 = buffer.data(sof + 105);
    const auto *sof_106 = buffer.data(sof + 106);
    const auto *sof_107 = buffer.data(sof + 107);
    const auto *sof_108 = buffer.data(sof + 108);
    const auto *sof_109 = buffer.data(sof + 109);
    const auto *sof_110 = buffer.data(sof + 110);
    const auto *sof_112 = buffer.data(sof + 112);
    const auto *sof_115 = buffer.data(sof + 115);
    const auto *sof_116 = buffer.data(sof + 116);
    const auto *sof_117 = buffer.data(sof + 117);
    const auto *sof_118 = buffer.data(sof + 118);
    const auto *sof_119 = buffer.data(sof + 119);
    const auto *sof_120 = buffer.data(sof + 120);
    const auto *sof_122 = buffer.data(sof + 122);
    const auto *sof_123 = buffer.data(sof + 123);
    const auto *sof_125 = buffer.data(sof + 125);
    const auto *sof_126 = buffer.data(sof + 126);
    const auto *sof_127 = buffer.data(sof + 127);
    const auto *sof_128 = buffer.data(sof + 128);
    const auto *sof_129 = buffer.data(sof + 129);
    const auto *sof_130 = buffer.data(sof + 130);
    const auto *sof_132 = buffer.data(sof + 132);
    const auto *sof_136 = buffer.data(sof + 136);
    const auto *sof_137 = buffer.data(sof + 137);
    const auto *sof_138 = buffer.data(sof + 138);
    const auto *sof_139 = buffer.data(sof + 139);
    const auto *sof_140 = buffer.data(sof + 140);
    const auto *sof_142 = buffer.data(sof + 142);
    const auto *sof_143 = buffer.data(sof + 143);
    const auto *sof_145 = buffer.data(sof + 145);
    const auto *sof_146 = buffer.data(sof + 146);
    const auto *sof_147 = buffer.data(sof + 147);
    const auto *sof_148 = buffer.data(sof + 148);
    const auto *sof_149 = buffer.data(sof + 149);
    const auto *sof_150 = buffer.data(sof + 150);
    const auto *sof_152 = buffer.data(sof + 152);
    const auto *sof_153 = buffer.data(sof + 153);
    const auto *sof_155 = buffer.data(sof + 155);
    const auto *sof_156 = buffer.data(sof + 156);
    const auto *sof_157 = buffer.data(sof + 157);
    const auto *sof_158 = buffer.data(sof + 158);
    const auto *sof_159 = buffer.data(sof + 159);
    const auto *sof_160 = buffer.data(sof + 160);
    const auto *sof_162 = buffer.data(sof + 162);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, sng0_80, snf_86, \
                         snf_87, snf_88, sng1_80, sof_86, sof_87, \
                         sof_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * sng0_80[k]
                   - f_6 * pc_y[k] * sng1_80[k];

        t_126[k] = f_11 * snf_86[k]
                   + f_3 * pc_x[k] * sof_86[k];

        t_127[k] = f_11 * snf_87[k]
                   + f_3 * pc_x[k] * sof_87[k];

        t_128[k] = f_11 * snf_88[k]
                   + f_3 * pc_x[k] * sof_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, pc_z, snf_46, snf_56, snf_89, \
                         sod0_51, sod1_51, sof_86, sof_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * snf_89[k]
                   + f_3 * pc_x[k] * sof_89[k];

        t_130[k] = f_7 * snf_56[k]
                   + f_1 * sod0_51[k]
                   - f_2 * sod1_51[k]
                   + f_3 * pc_y[k] * sof_86[k];

        t_131[k] = f_8 * snf_46[k]
                   + f_3 * pc_z[k] * sof_86[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pc_y, sng0_89, snf_58, snf_59, sng1_89, \
                         sod0_53, sod1_53, sof_88, sof_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_7 * snf_58[k]
                   + f_4 * sod0_53[k]
                   - f_5 * sod1_53[k]
                   + f_3 * pc_y[k] * sof_88[k];

        t_133[k] = f_7 * snf_59[k]
                   + f_3 * pc_y[k] * sof_89[k];

        t_134[k] = pb_y[k] * sng0_89[k]
                   - f_6 * pc_y[k] * sng1_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, snf_50, snf_90, snf_93, \
                         sod0_54, sod0_57, sod1_54, sod1_57, sof_90, \
                         sof_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_11 * snf_90[k]
                   + f_1 * sod0_54[k]
                   - f_2 * sod1_54[k]
                   + f_3 * pc_x[k] * sof_90[k];

        t_136[k] = f_3 * pc_y[k] * sof_90[k];

        t_137[k] = f_12 * snf_50[k]
                   + f_3 * pc_z[k] * sof_90[k];

        t_138[k] = f_11 * snf_93[k]
                   + f_4 * sod0_57[k]
                   - f_5 * sod1_57[k]
                   + f_3 * pc_x[k] * sof_93[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, snf_95, snf_96, snf_97, \
                         sod0_59, sod1_59, sof_92, sof_95, sof_96, \
                         sof_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * pc_y[k] * sof_92[k];

        t_140[k] = f_11 * snf_95[k]
                   + f_4 * sod0_59[k]
                   - f_5 * sod1_59[k]
                   + f_3 * pc_x[k] * sof_95[k];

        t_141[k] = f_11 * snf_96[k]
                   + f_3 * pc_x[k] * sof_96[k];

        t_142[k] = f_11 * snf_97[k]
                   + f_3 * pc_x[k] * sof_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, snf_56, snf_98, snf_99, \
                         sod0_57, sod1_57, sof_96, sof_98, sof_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * snf_98[k]
                   + f_3 * pc_x[k] * sof_98[k];

        t_144[k] = f_11 * snf_99[k]
                   + f_3 * pc_x[k] * sof_99[k];

        t_145[k] = f_1 * sod0_57[k]
                   - f_2 * sod1_57[k]
                   + f_3 * pc_y[k] * sof_96[k];

        t_146[k] = f_12 * snf_56[k]
                   + f_3 * pc_z[k] * sof_96[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, snf_59, snf_100, \
                         sod0_59, sod0_60, sod1_59, sod1_60, sof_98, sof_99, \
                         sof_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * sod0_59[k]
                   - f_5 * sod1_59[k]
                   + f_3 * pc_y[k] * sof_98[k];

        t_148[k] = f_3 * pc_y[k] * sof_99[k];

        t_149[k] = f_12 * snf_59[k]
                   + f_1 * sod0_59[k]
                   - f_2 * sod1_59[k]
                   + f_3 * pc_z[k] * sof_99[k];

        t_150[k] = f_13 * snf_100[k]
                   + f_1 * sod0_60[k]
                   - f_2 * sod1_60[k]
                   + f_3 * pc_x[k] * sof_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, snf_60, snf_62, \
                         snf_103, sod0_63, sod1_63, sof_100, sof_102, \
                         sof_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_14 * snf_60[k]
                   + f_3 * pc_y[k] * sof_100[k];

        t_152[k] = f_3 * pc_z[k] * sof_100[k];

        t_153[k] = f_13 * snf_103[k]
                   + f_4 * sod0_63[k]
                   - f_5 * sod1_63[k]
                   + f_3 * pc_x[k] * sof_103[k];

        t_154[k] = f_14 * snf_62[k]
                   + f_3 * pc_y[k] * sof_102[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, snf_105, snf_106, snf_107, snf_108, \
                         sod0_65, sod1_65, sof_105, sof_106, sof_107, \
                         sof_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_13 * snf_105[k]
                   + f_4 * sod0_65[k]
                   - f_5 * sod1_65[k]
                   + f_3 * pc_x[k] * sof_105[k];

        t_156[k] = f_13 * snf_106[k]
                   + f_3 * pc_x[k] * sof_106[k];

        t_157[k] = f_13 * snf_107[k]
                   + f_3 * pc_x[k] * sof_107[k];

        t_158[k] = f_13 * snf_108[k]
                   + f_3 * pc_x[k] * sof_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, pc_z, snf_66, snf_109, sod0_63, \
                         sod1_63, sof_106, sof_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_13 * snf_109[k]
                   + f_3 * pc_x[k] * sof_109[k];

        t_160[k] = f_14 * snf_66[k]
                   + f_1 * sod0_63[k]
                   - f_2 * sod1_63[k]
                   + f_3 * pc_y[k] * sof_106[k];

        t_161[k] = f_3 * pc_z[k] * sof_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_z, pc_y, pc_z, sng0_90, snf_68, \
                         snf_69, sng1_90, sod0_65, sod1_65, sof_108, \
                         sof_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_14 * snf_68[k]
                   + f_4 * sod0_65[k]
                   - f_5 * sod1_65[k]
                   + f_3 * pc_y[k] * sof_108[k];

        t_163[k] = f_14 * snf_69[k]
                   + f_3 * pc_y[k] * sof_109[k];

        t_164[k] = f_1 * sod0_65[k]
                   - f_2 * sod1_65[k]
                   + f_3 * pc_z[k] * sof_109[k];

        t_165[k] = pb_z[k] * sng0_90[k]
                   - f_6 * pc_z[k] * sng1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, sng0_93, snf_60, \
                         snf_70, snf_72, sng1_93, sof_110, sof_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_12 * snf_70[k]
                   + f_3 * pc_y[k] * sof_110[k];

        t_167[k] = f_7 * snf_60[k]
                   + f_3 * pc_z[k] * sof_110[k];

        t_168[k] = pb_z[k] * sng0_93[k]
                   - f_6 * pc_z[k] * sng1_93[k];

        t_169[k] = f_12 * snf_72[k]
                   + f_3 * pc_y[k] * sof_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, snf_115, snf_116, snf_117, snf_118, \
                         sod0_71, sod1_71, sof_115, sof_116, sof_117, \
                         sof_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_13 * snf_115[k]
                   + f_4 * sod0_71[k]
                   - f_5 * sod1_71[k]
                   + f_3 * pc_x[k] * sof_115[k];

        t_171[k] = f_13 * snf_116[k]
                   + f_3 * pc_x[k] * sof_116[k];

        t_172[k] = f_13 * snf_117[k]
                   + f_3 * pc_x[k] * sof_117[k];

        t_173[k] = f_13 * snf_118[k]
                   + f_3 * pc_x[k] * sof_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_z, pc_x, pc_z, sng0_100, snf_66, snf_119, \
                         sng1_100, sof_116, sof_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * snf_119[k]
                   + f_3 * pc_x[k] * sof_119[k];

        t_175[k] = pb_z[k] * sng0_100[k]
                   - f_6 * pc_z[k] * sng1_100[k];

        t_176[k] = f_7 * snf_66[k]
                   + f_3 * pc_z[k] * sof_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pc_y, pc_z, snf_69, snf_78, snf_79, sod0_71, \
                         sod1_71, sof_118, sof_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_12 * snf_78[k]
                   + f_4 * sod0_71[k]
                   - f_5 * sod1_71[k]
                   + f_3 * pc_y[k] * sof_118[k];

        t_178[k] = f_12 * snf_79[k]
                   + f_3 * pc_y[k] * sof_119[k];

        t_179[k] = f_7 * snf_69[k]
                   + f_1 * sod0_71[k]
                   - f_2 * sod1_71[k]
                   + f_3 * pc_z[k] * sof_119[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pc_x, pc_y, pc_z, snf_70, snf_80, snf_120, \
                         sod0_72, sod1_72, sof_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_13 * snf_120[k]
                   + f_1 * sod0_72[k]
                   - f_2 * sod1_72[k]
                   + f_3 * pc_x[k] * sof_120[k];

        t_181[k] = f_8 * snf_80[k]
                   + f_3 * pc_y[k] * sof_120[k];

        t_182[k] = f_8 * snf_70[k]
                   + f_3 * pc_z[k] * sof_120[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_x, pc_y, snf_82, snf_123, snf_125, sod0_75, \
                         sod0_77, sod1_75, sod1_77, sof_122, sof_123, \
                         sof_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_13 * snf_123[k]
                   + f_4 * sod0_75[k]
                   - f_5 * sod1_75[k]
                   + f_3 * pc_x[k] * sof_123[k];

        t_184[k] = f_8 * snf_82[k]
                   + f_3 * pc_y[k] * sof_122[k];

        t_185[k] = f_13 * snf_125[k]
                   + f_4 * sod0_77[k]
                   - f_5 * sod1_77[k]
                   + f_3 * pc_x[k] * sof_125[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, snf_126, snf_127, snf_128, snf_129, \
                         sof_126, sof_127, sof_128, sof_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_13 * snf_126[k]
                   + f_3 * pc_x[k] * sof_126[k];

        t_187[k] = f_13 * snf_127[k]
                   + f_3 * pc_x[k] * sof_127[k];

        t_188[k] = f_13 * snf_128[k]
                   + f_3 * pc_x[k] * sof_128[k];

        t_189[k] = f_13 * snf_129[k]
                   + f_3 * pc_x[k] * sof_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_y, pc_z, snf_76, snf_86, snf_88, sod0_75, \
                         sod0_77, sod1_75, sod1_77, sof_126, sof_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * snf_86[k]
                   + f_1 * sod0_75[k]
                   - f_2 * sod1_75[k]
                   + f_3 * pc_y[k] * sof_126[k];

        t_191[k] = f_8 * snf_76[k]
                   + f_3 * pc_z[k] * sof_126[k];

        t_192[k] = f_8 * snf_88[k]
                   + f_4 * sod0_77[k]
                   - f_5 * sod1_77[k]
                   + f_3 * pc_y[k] * sof_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_y, pc_y, pc_z, sng0_135, snf_79, \
                         snf_89, snf_90, sng1_135, sod0_77, sod1_77, sof_129, \
                         sof_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * snf_89[k]
                   + f_3 * pc_y[k] * sof_129[k];

        t_194[k] = f_8 * snf_79[k]
                   + f_1 * sod0_77[k]
                   - f_2 * sod1_77[k]
                   + f_3 * pc_z[k] * sof_129[k];

        t_195[k] = pb_y[k] * sng0_135[k]
                   - f_6 * pc_y[k] * sng1_135[k];

        t_196[k] = f_7 * snf_90[k]
                   + f_3 * pc_y[k] * sof_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_y, pc_y, pc_z, sng0_138, sng0_140, \
                         snf_80, snf_91, snf_92, sng1_138, sng1_140, sof_130, \
                         sof_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * snf_80[k]
                   + f_3 * pc_z[k] * sof_130[k];

        t_198[k] = pb_y[k] * sng0_138[k]
                   + f_8 * snf_91[k]
                   - f_6 * pc_y[k] * sng1_138[k];

        t_199[k] = f_7 * snf_92[k]
                   + f_3 * pc_y[k] * sof_132[k];

        t_200[k] = pb_y[k] * sng0_140[k]
                   - f_6 * pc_y[k] * sng1_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, snf_136, snf_137, snf_138, snf_139, \
                         sof_136, sof_137, sof_138, sof_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_13 * snf_136[k]
                   + f_3 * pc_x[k] * sof_136[k];

        t_202[k] = f_13 * snf_137[k]
                   + f_3 * pc_x[k] * sof_137[k];

        t_203[k] = f_13 * snf_138[k]
                   + f_3 * pc_x[k] * sof_138[k];

        t_204[k] = f_13 * snf_139[k]
                   + f_3 * pc_x[k] * sof_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pc_y, pc_z, snf_86, snf_96, snf_98, sod0_81, \
                         sod0_83, sod1_81, sod1_83, sof_136, sof_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_7 * snf_96[k]
                   + f_1 * sod0_81[k]
                   - f_2 * sod1_81[k]
                   + f_3 * pc_y[k] * sof_136[k];

        t_206[k] = f_12 * snf_86[k]
                   + f_3 * pc_z[k] * sof_136[k];

        t_207[k] = f_7 * snf_98[k]
                   + f_4 * sod0_83[k]
                   - f_5 * sod1_83[k]
                   + f_3 * pc_y[k] * sof_138[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pc_x, pc_y, sng0_149, snf_99, \
                         snf_140, sng1_149, sod0_84, sod1_84, sof_139, \
                         sof_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * snf_99[k]
                   + f_3 * pc_y[k] * sof_139[k];

        t_209[k] = pb_y[k] * sng0_149[k]
                   - f_6 * pc_y[k] * sng1_149[k];

        t_210[k] = f_13 * snf_140[k]
                   + f_1 * sod0_84[k]
                   - f_2 * sod1_84[k]
                   + f_3 * pc_x[k] * sof_140[k];

        t_211[k] = f_3 * pc_y[k] * sof_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, snf_90, snf_143, sod0_87, \
                         sod1_87, sof_140, sof_142, sof_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_14 * snf_90[k]
                   + f_3 * pc_z[k] * sof_140[k];

        t_213[k] = f_13 * snf_143[k]
                   + f_4 * sod0_87[k]
                   - f_5 * sod1_87[k]
                   + f_3 * pc_x[k] * sof_143[k];

        t_214[k] = f_3 * pc_y[k] * sof_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, snf_145, snf_146, snf_147, snf_148, \
                         sod0_89, sod1_89, sof_145, sof_146, sof_147, \
                         sof_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_13 * snf_145[k]
                   + f_4 * sod0_89[k]
                   - f_5 * sod1_89[k]
                   + f_3 * pc_x[k] * sof_145[k];

        t_216[k] = f_13 * snf_146[k]
                   + f_3 * pc_x[k] * sof_146[k];

        t_217[k] = f_13 * snf_147[k]
                   + f_3 * pc_x[k] * sof_147[k];

        t_218[k] = f_13 * snf_148[k]
                   + f_3 * pc_x[k] * sof_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pc_x, pc_y, pc_z, snf_96, snf_149, \
                         sod0_87, sod0_89, sod1_87, sod1_89, sof_146, sof_148, \
                         sof_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_13 * snf_149[k]
                   + f_3 * pc_x[k] * sof_149[k];

        t_220[k] = f_1 * sod0_87[k]
                   - f_2 * sod1_87[k]
                   + f_3 * pc_y[k] * sof_146[k];

        t_221[k] = f_14 * snf_96[k]
                   + f_3 * pc_z[k] * sof_146[k];

        t_222[k] = f_4 * sod0_89[k]
                   - f_5 * sod1_89[k]
                   + f_3 * pc_y[k] * sof_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pc_x, pc_y, pc_z, snf_99, snf_100, \
                         snf_150, sod0_89, sod0_90, sod1_89, sod1_90, sof_149, \
                         sof_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_y[k] * sof_149[k];

        t_224[k] = f_14 * snf_99[k]
                   + f_1 * sod0_89[k]
                   - f_2 * sod1_89[k]
                   + f_3 * pc_z[k] * sof_149[k];

        t_225[k] = f_15 * snf_150[k]
                   + f_1 * sod0_90[k]
                   - f_2 * sod1_90[k]
                   + f_3 * pc_x[k] * sof_150[k];

        t_226[k] = f_16 * snf_100[k]
                   + f_3 * pc_y[k] * sof_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, snf_102, snf_153, sod0_93, \
                         sod1_93, sof_150, sof_152, sof_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_3 * pc_z[k] * sof_150[k];

        t_228[k] = f_15 * snf_153[k]
                   + f_4 * sod0_93[k]
                   - f_5 * sod1_93[k]
                   + f_3 * pc_x[k] * sof_153[k];

        t_229[k] = f_16 * snf_102[k]
                   + f_3 * pc_y[k] * sof_152[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, snf_155, snf_156, snf_157, snf_158, \
                         sod0_95, sod1_95, sof_155, sof_156, sof_157, \
                         sof_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_15 * snf_155[k]
                   + f_4 * sod0_95[k]
                   - f_5 * sod1_95[k]
                   + f_3 * pc_x[k] * sof_155[k];

        t_231[k] = f_15 * snf_156[k]
                   + f_3 * pc_x[k] * sof_156[k];

        t_232[k] = f_15 * snf_157[k]
                   + f_3 * pc_x[k] * sof_157[k];

        t_233[k] = f_15 * snf_158[k]
                   + f_3 * pc_x[k] * sof_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, pc_z, snf_106, snf_159, sod0_93, \
                         sod1_93, sof_156, sof_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_15 * snf_159[k]
                   + f_3 * pc_x[k] * sof_159[k];

        t_235[k] = f_16 * snf_106[k]
                   + f_1 * sod0_93[k]
                   - f_2 * sod1_93[k]
                   + f_3 * pc_y[k] * sof_156[k];

        t_236[k] = f_3 * pc_z[k] * sof_156[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pb_z, pc_y, pc_z, sng0_150, snf_108, \
                         snf_109, sng1_150, sod0_95, sod1_95, sof_158, \
                         sof_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_16 * snf_108[k]
                   + f_4 * sod0_95[k]
                   - f_5 * sod1_95[k]
                   + f_3 * pc_y[k] * sof_158[k];

        t_238[k] = f_16 * snf_109[k]
                   + f_3 * pc_y[k] * sof_159[k];

        t_239[k] = f_1 * sod0_95[k]
                   - f_2 * sod1_95[k]
                   + f_3 * pc_z[k] * sof_159[k];

        t_240[k] = pb_z[k] * sng0_150[k]
                   - f_6 * pc_z[k] * sng1_150[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_z, pc_y, pc_z, sng0_153, snf_100, \
                         snf_110, snf_112, sng1_153, sof_160, sof_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_14 * snf_110[k]
                   + f_3 * pc_y[k] * sof_160[k];

        t_242[k] = f_7 * snf_100[k]
                   + f_3 * pc_z[k] * sof_160[k];

        t_243[k] = pb_z[k] * sng0_153[k]
                   - f_6 * pc_z[k] * sng1_153[k];

        t_244[k] = f_14 * snf_112[k]
                   + f_3 * pc_y[k] * sof_162[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
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
    const auto f_12 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_160 = buffer.data(sng0 + 160);
    const auto *sng0_210 = buffer.data(sng0 + 210);
    const auto *sng0_213 = buffer.data(sng0 + 213);
    const auto *sng0_215 = buffer.data(sng0 + 215);
    const auto *sng0_224 = buffer.data(sng0 + 224);
    const auto *sng0_225 = buffer.data(sng0 + 225);
    const auto *sng0_228 = buffer.data(sng0 + 228);
    const auto *sng0_235 = buffer.data(sng0 + 235);

    const auto *snf_106 = buffer.data(snf + 106);
    const auto *snf_109 = buffer.data(snf + 109);
    const auto *snf_110 = buffer.data(snf + 110);
    const auto *snf_116 = buffer.data(snf + 116);
    const auto *snf_118 = buffer.data(snf + 118);
    const auto *snf_119 = buffer.data(snf + 119);
    const auto *snf_120 = buffer.data(snf + 120);
    const auto *snf_122 = buffer.data(snf + 122);
    const auto *snf_126 = buffer.data(snf + 126);
    const auto *snf_128 = buffer.data(snf + 128);
    const auto *snf_129 = buffer.data(snf + 129);
    const auto *snf_130 = buffer.data(snf + 130);
    const auto *snf_132 = buffer.data(snf + 132);
    const auto *snf_136 = buffer.data(snf + 136);
    const auto *snf_138 = buffer.data(snf + 138);
    const auto *snf_139 = buffer.data(snf + 139);
    const auto *snf_140 = buffer.data(snf + 140);
    const auto *snf_141 = buffer.data(snf + 141);
    const auto *snf_142 = buffer.data(snf + 142);
    const auto *snf_146 = buffer.data(snf + 146);
    const auto *snf_148 = buffer.data(snf + 148);
    const auto *snf_149 = buffer.data(snf + 149);
    const auto *snf_150 = buffer.data(snf + 150);
    const auto *snf_152 = buffer.data(snf + 152);
    const auto *snf_156 = buffer.data(snf + 156);
    const auto *snf_158 = buffer.data(snf + 158);
    const auto *snf_159 = buffer.data(snf + 159);
    const auto *snf_160 = buffer.data(snf + 160);
    const auto *snf_162 = buffer.data(snf + 162);
    const auto *snf_165 = buffer.data(snf + 165);
    const auto *snf_166 = buffer.data(snf + 166);
    const auto *snf_167 = buffer.data(snf + 167);
    const auto *snf_168 = buffer.data(snf + 168);
    const auto *snf_169 = buffer.data(snf + 169);
    const auto *snf_170 = buffer.data(snf + 170);
    const auto *snf_172 = buffer.data(snf + 172);
    const auto *snf_173 = buffer.data(snf + 173);
    const auto *snf_175 = buffer.data(snf + 175);
    const auto *snf_176 = buffer.data(snf + 176);
    const auto *snf_177 = buffer.data(snf + 177);
    const auto *snf_178 = buffer.data(snf + 178);
    const auto *snf_179 = buffer.data(snf + 179);
    const auto *snf_180 = buffer.data(snf + 180);
    const auto *snf_183 = buffer.data(snf + 183);
    const auto *snf_185 = buffer.data(snf + 185);
    const auto *snf_186 = buffer.data(snf + 186);
    const auto *snf_187 = buffer.data(snf + 187);
    const auto *snf_188 = buffer.data(snf + 188);
    const auto *snf_189 = buffer.data(snf + 189);
    const auto *snf_196 = buffer.data(snf + 196);
    const auto *snf_197 = buffer.data(snf + 197);
    const auto *snf_198 = buffer.data(snf + 198);
    const auto *snf_199 = buffer.data(snf + 199);
    const auto *snf_200 = buffer.data(snf + 200);
    const auto *snf_203 = buffer.data(snf + 203);
    const auto *snf_205 = buffer.data(snf + 205);
    const auto *snf_206 = buffer.data(snf + 206);
    const auto *snf_207 = buffer.data(snf + 207);
    const auto *snf_208 = buffer.data(snf + 208);
    const auto *snf_209 = buffer.data(snf + 209);
    const auto *snf_210 = buffer.data(snf + 210);
    const auto *snf_213 = buffer.data(snf + 213);
    const auto *snf_215 = buffer.data(snf + 215);
    const auto *snf_216 = buffer.data(snf + 216);
    const auto *snf_217 = buffer.data(snf + 217);
    const auto *snf_218 = buffer.data(snf + 218);
    const auto *snf_219 = buffer.data(snf + 219);
    const auto *snf_225 = buffer.data(snf + 225);
    const auto *snf_226 = buffer.data(snf + 226);
    const auto *snf_227 = buffer.data(snf + 227);
    const auto *snf_228 = buffer.data(snf + 228);
    const auto *snf_229 = buffer.data(snf + 229);
    const auto *snf_230 = buffer.data(snf + 230);
    const auto *snf_233 = buffer.data(snf + 233);
    const auto *snf_235 = buffer.data(snf + 235);
    const auto *snf_236 = buffer.data(snf + 236);
    const auto *snf_237 = buffer.data(snf + 237);
    const auto *snf_238 = buffer.data(snf + 238);
    const auto *snf_239 = buffer.data(snf + 239);
    const auto *snf_240 = buffer.data(snf + 240);

    const auto *sng1_160 = buffer.data(sng1 + 160);
    const auto *sng1_210 = buffer.data(sng1 + 210);
    const auto *sng1_213 = buffer.data(sng1 + 213);
    const auto *sng1_215 = buffer.data(sng1 + 215);
    const auto *sng1_224 = buffer.data(sng1 + 224);
    const auto *sng1_225 = buffer.data(sng1 + 225);
    const auto *sng1_228 = buffer.data(sng1 + 228);
    const auto *sng1_235 = buffer.data(sng1 + 235);

    const auto *sod0_101 = buffer.data(sod0 + 101);
    const auto *sod0_102 = buffer.data(sod0 + 102);
    const auto *sod0_105 = buffer.data(sod0 + 105);
    const auto *sod0_107 = buffer.data(sod0 + 107);
    const auto *sod0_108 = buffer.data(sod0 + 108);
    const auto *sod0_111 = buffer.data(sod0 + 111);
    const auto *sod0_113 = buffer.data(sod0 + 113);
    const auto *sod0_117 = buffer.data(sod0 + 117);
    const auto *sod0_119 = buffer.data(sod0 + 119);
    const auto *sod0_120 = buffer.data(sod0 + 120);
    const auto *sod0_123 = buffer.data(sod0 + 123);
    const auto *sod0_125 = buffer.data(sod0 + 125);
    const auto *sod0_126 = buffer.data(sod0 + 126);
    const auto *sod0_129 = buffer.data(sod0 + 129);
    const auto *sod0_131 = buffer.data(sod0 + 131);
    const auto *sod0_137 = buffer.data(sod0 + 137);
    const auto *sod0_138 = buffer.data(sod0 + 138);
    const auto *sod0_141 = buffer.data(sod0 + 141);
    const auto *sod0_143 = buffer.data(sod0 + 143);
    const auto *sod0_144 = buffer.data(sod0 + 144);

    const auto *sod1_101 = buffer.data(sod1 + 101);
    const auto *sod1_102 = buffer.data(sod1 + 102);
    const auto *sod1_105 = buffer.data(sod1 + 105);
    const auto *sod1_107 = buffer.data(sod1 + 107);
    const auto *sod1_108 = buffer.data(sod1 + 108);
    const auto *sod1_111 = buffer.data(sod1 + 111);
    const auto *sod1_113 = buffer.data(sod1 + 113);
    const auto *sod1_117 = buffer.data(sod1 + 117);
    const auto *sod1_119 = buffer.data(sod1 + 119);
    const auto *sod1_120 = buffer.data(sod1 + 120);
    const auto *sod1_123 = buffer.data(sod1 + 123);
    const auto *sod1_125 = buffer.data(sod1 + 125);
    const auto *sod1_126 = buffer.data(sod1 + 126);
    const auto *sod1_129 = buffer.data(sod1 + 129);
    const auto *sod1_131 = buffer.data(sod1 + 131);
    const auto *sod1_137 = buffer.data(sod1 + 137);
    const auto *sod1_138 = buffer.data(sod1 + 138);
    const auto *sod1_141 = buffer.data(sod1 + 141);
    const auto *sod1_143 = buffer.data(sod1 + 143);
    const auto *sod1_144 = buffer.data(sod1 + 144);

    const auto *sof_165 = buffer.data(sof + 165);
    const auto *sof_166 = buffer.data(sof + 166);
    const auto *sof_167 = buffer.data(sof + 167);
    const auto *sof_168 = buffer.data(sof + 168);
    const auto *sof_169 = buffer.data(sof + 169);
    const auto *sof_170 = buffer.data(sof + 170);
    const auto *sof_172 = buffer.data(sof + 172);
    const auto *sof_173 = buffer.data(sof + 173);
    const auto *sof_175 = buffer.data(sof + 175);
    const auto *sof_176 = buffer.data(sof + 176);
    const auto *sof_177 = buffer.data(sof + 177);
    const auto *sof_178 = buffer.data(sof + 178);
    const auto *sof_179 = buffer.data(sof + 179);
    const auto *sof_180 = buffer.data(sof + 180);
    const auto *sof_182 = buffer.data(sof + 182);
    const auto *sof_183 = buffer.data(sof + 183);
    const auto *sof_185 = buffer.data(sof + 185);
    const auto *sof_186 = buffer.data(sof + 186);
    const auto *sof_187 = buffer.data(sof + 187);
    const auto *sof_188 = buffer.data(sof + 188);
    const auto *sof_189 = buffer.data(sof + 189);
    const auto *sof_190 = buffer.data(sof + 190);
    const auto *sof_192 = buffer.data(sof + 192);
    const auto *sof_196 = buffer.data(sof + 196);
    const auto *sof_197 = buffer.data(sof + 197);
    const auto *sof_198 = buffer.data(sof + 198);
    const auto *sof_199 = buffer.data(sof + 199);
    const auto *sof_200 = buffer.data(sof + 200);
    const auto *sof_202 = buffer.data(sof + 202);
    const auto *sof_203 = buffer.data(sof + 203);
    const auto *sof_205 = buffer.data(sof + 205);
    const auto *sof_206 = buffer.data(sof + 206);
    const auto *sof_207 = buffer.data(sof + 207);
    const auto *sof_208 = buffer.data(sof + 208);
    const auto *sof_209 = buffer.data(sof + 209);
    const auto *sof_210 = buffer.data(sof + 210);
    const auto *sof_212 = buffer.data(sof + 212);
    const auto *sof_213 = buffer.data(sof + 213);
    const auto *sof_215 = buffer.data(sof + 215);
    const auto *sof_216 = buffer.data(sof + 216);
    const auto *sof_217 = buffer.data(sof + 217);
    const auto *sof_218 = buffer.data(sof + 218);
    const auto *sof_219 = buffer.data(sof + 219);
    const auto *sof_220 = buffer.data(sof + 220);
    const auto *sof_222 = buffer.data(sof + 222);
    const auto *sof_225 = buffer.data(sof + 225);
    const auto *sof_226 = buffer.data(sof + 226);
    const auto *sof_227 = buffer.data(sof + 227);
    const auto *sof_228 = buffer.data(sof + 228);
    const auto *sof_229 = buffer.data(sof + 229);
    const auto *sof_230 = buffer.data(sof + 230);
    const auto *sof_232 = buffer.data(sof + 232);
    const auto *sof_233 = buffer.data(sof + 233);
    const auto *sof_235 = buffer.data(sof + 235);
    const auto *sof_236 = buffer.data(sof + 236);
    const auto *sof_237 = buffer.data(sof + 237);
    const auto *sof_238 = buffer.data(sof + 238);
    const auto *sof_239 = buffer.data(sof + 239);
    const auto *sof_240 = buffer.data(sof + 240);

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, snf_165, snf_166, snf_167, snf_168, \
                         sod0_101, sod1_101, sof_165, sof_166, sof_167, \
                         sof_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_15 * snf_165[k]
                   + f_4 * sod0_101[k]
                   - f_5 * sod1_101[k]
                   + f_3 * pc_x[k] * sof_165[k];

        t_246[k] = f_15 * snf_166[k]
                   + f_3 * pc_x[k] * sof_166[k];

        t_247[k] = f_15 * snf_167[k]
                   + f_3 * pc_x[k] * sof_167[k];

        t_248[k] = f_15 * snf_168[k]
                   + f_3 * pc_x[k] * sof_168[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pb_z, pc_x, pc_z, sng0_160, snf_106, snf_169, \
                         sng1_160, sof_166, sof_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_15 * snf_169[k]
                   + f_3 * pc_x[k] * sof_169[k];

        t_250[k] = pb_z[k] * sng0_160[k]
                   - f_6 * pc_z[k] * sng1_160[k];

        t_251[k] = f_7 * snf_106[k]
                   + f_3 * pc_z[k] * sof_166[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pc_y, pc_z, snf_109, snf_118, snf_119, sod0_101, \
                         sod1_101, sof_168, sof_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * snf_118[k]
                   + f_4 * sod0_101[k]
                   - f_5 * sod1_101[k]
                   + f_3 * pc_y[k] * sof_168[k];

        t_253[k] = f_14 * snf_119[k]
                   + f_3 * pc_y[k] * sof_169[k];

        t_254[k] = f_7 * snf_109[k]
                   + f_1 * sod0_101[k]
                   - f_2 * sod1_101[k]
                   + f_3 * pc_z[k] * sof_169[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, pc_z, snf_110, snf_120, snf_170, \
                         sod0_102, sod1_102, sof_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_15 * snf_170[k]
                   + f_1 * sod0_102[k]
                   - f_2 * sod1_102[k]
                   + f_3 * pc_x[k] * sof_170[k];

        t_256[k] = f_12 * snf_120[k]
                   + f_3 * pc_y[k] * sof_170[k];

        t_257[k] = f_8 * snf_110[k]
                   + f_3 * pc_z[k] * sof_170[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, snf_122, snf_173, snf_175, sod0_105, \
                         sod0_107, sod1_105, sod1_107, sof_172, sof_173, \
                         sof_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_15 * snf_173[k]
                   + f_4 * sod0_105[k]
                   - f_5 * sod1_105[k]
                   + f_3 * pc_x[k] * sof_173[k];

        t_259[k] = f_12 * snf_122[k]
                   + f_3 * pc_y[k] * sof_172[k];

        t_260[k] = f_15 * snf_175[k]
                   + f_4 * sod0_107[k]
                   - f_5 * sod1_107[k]
                   + f_3 * pc_x[k] * sof_175[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pc_x, snf_176, snf_177, snf_178, snf_179, \
                         sof_176, sof_177, sof_178, sof_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_15 * snf_176[k]
                   + f_3 * pc_x[k] * sof_176[k];

        t_262[k] = f_15 * snf_177[k]
                   + f_3 * pc_x[k] * sof_177[k];

        t_263[k] = f_15 * snf_178[k]
                   + f_3 * pc_x[k] * sof_178[k];

        t_264[k] = f_15 * snf_179[k]
                   + f_3 * pc_x[k] * sof_179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pc_y, pc_z, snf_116, snf_126, snf_128, sod0_105, \
                         sod0_107, sod1_105, sod1_107, sof_176, \
                         sof_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_12 * snf_126[k]
                   + f_1 * sod0_105[k]
                   - f_2 * sod1_105[k]
                   + f_3 * pc_y[k] * sof_176[k];

        t_266[k] = f_8 * snf_116[k]
                   + f_3 * pc_z[k] * sof_176[k];

        t_267[k] = f_12 * snf_128[k]
                   + f_4 * sod0_107[k]
                   - f_5 * sod1_107[k]
                   + f_3 * pc_y[k] * sof_178[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_x, pc_y, pc_z, snf_119, snf_129, snf_180, \
                         sod0_107, sod0_108, sod1_107, sod1_108, sof_179, \
                         sof_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_12 * snf_129[k]
                   + f_3 * pc_y[k] * sof_179[k];

        t_269[k] = f_8 * snf_119[k]
                   + f_1 * sod0_107[k]
                   - f_2 * sod1_107[k]
                   + f_3 * pc_z[k] * sof_179[k];

        t_270[k] = f_15 * snf_180[k]
                   + f_1 * sod0_108[k]
                   - f_2 * sod1_108[k]
                   + f_3 * pc_x[k] * sof_180[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, snf_120, snf_130, \
                         snf_132, snf_183, sod0_111, sod1_111, sof_180, sof_182, \
                         sof_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_8 * snf_130[k]
                   + f_3 * pc_y[k] * sof_180[k];

        t_272[k] = f_12 * snf_120[k]
                   + f_3 * pc_z[k] * sof_180[k];

        t_273[k] = f_15 * snf_183[k]
                   + f_4 * sod0_111[k]
                   - f_5 * sod1_111[k]
                   + f_3 * pc_x[k] * sof_183[k];

        t_274[k] = f_8 * snf_132[k]
                   + f_3 * pc_y[k] * sof_182[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pc_x, snf_185, snf_186, snf_187, snf_188, \
                         sod0_113, sod1_113, sof_185, sof_186, sof_187, \
                         sof_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_15 * snf_185[k]
                   + f_4 * sod0_113[k]
                   - f_5 * sod1_113[k]
                   + f_3 * pc_x[k] * sof_185[k];

        t_276[k] = f_15 * snf_186[k]
                   + f_3 * pc_x[k] * sof_186[k];

        t_277[k] = f_15 * snf_187[k]
                   + f_3 * pc_x[k] * sof_187[k];

        t_278[k] = f_15 * snf_188[k]
                   + f_3 * pc_x[k] * sof_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pc_x, pc_y, pc_z, snf_126, snf_136, snf_189, \
                         sod0_111, sod1_111, sof_186, sof_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_15 * snf_189[k]
                   + f_3 * pc_x[k] * sof_189[k];

        t_280[k] = f_8 * snf_136[k]
                   + f_1 * sod0_111[k]
                   - f_2 * sod1_111[k]
                   + f_3 * pc_y[k] * sof_186[k];

        t_281[k] = f_12 * snf_126[k]
                   + f_3 * pc_z[k] * sof_186[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_y, pc_y, pc_z, sng0_210, snf_129, \
                         snf_138, snf_139, sng1_210, sod0_113, sod1_113, sof_188, \
                         sof_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_8 * snf_138[k]
                   + f_4 * sod0_113[k]
                   - f_5 * sod1_113[k]
                   + f_3 * pc_y[k] * sof_188[k];

        t_283[k] = f_8 * snf_139[k]
                   + f_3 * pc_y[k] * sof_189[k];

        t_284[k] = f_12 * snf_129[k]
                   + f_1 * sod0_113[k]
                   - f_2 * sod1_113[k]
                   + f_3 * pc_z[k] * sof_189[k];

        t_285[k] = pb_y[k] * sng0_210[k]
                   - f_6 * pc_y[k] * sng1_210[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, sng0_213, snf_130, \
                         snf_140, snf_141, snf_142, sng1_213, sof_190, \
                         sof_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_7 * snf_140[k]
                   + f_3 * pc_y[k] * sof_190[k];

        t_287[k] = f_14 * snf_130[k]
                   + f_3 * pc_z[k] * sof_190[k];

        t_288[k] = pb_y[k] * sng0_213[k]
                   + f_8 * snf_141[k]
                   - f_6 * pc_y[k] * sng1_213[k];

        t_289[k] = f_7 * snf_142[k]
                   + f_3 * pc_y[k] * sof_192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_x, pc_y, sng0_215, snf_196, \
                         snf_197, snf_198, sng1_215, sof_196, sof_197, \
                         sof_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_y[k] * sng0_215[k]
                   - f_6 * pc_y[k] * sng1_215[k];

        t_291[k] = f_15 * snf_196[k]
                   + f_3 * pc_x[k] * sof_196[k];

        t_292[k] = f_15 * snf_197[k]
                   + f_3 * pc_x[k] * sof_197[k];

        t_293[k] = f_15 * snf_198[k]
                   + f_3 * pc_x[k] * sof_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pc_x, pc_y, pc_z, snf_136, snf_146, snf_199, \
                         sod0_117, sod1_117, sof_196, sof_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_15 * snf_199[k]
                   + f_3 * pc_x[k] * sof_199[k];

        t_295[k] = f_7 * snf_146[k]
                   + f_1 * sod0_117[k]
                   - f_2 * sod1_117[k]
                   + f_3 * pc_y[k] * sof_196[k];

        t_296[k] = f_14 * snf_136[k]
                   + f_3 * pc_z[k] * sof_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_y, pc_y, sng0_224, snf_148, snf_149, \
                         sng1_224, sod0_119, sod1_119, sof_198, \
                         sof_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_7 * snf_148[k]
                   + f_4 * sod0_119[k]
                   - f_5 * sod1_119[k]
                   + f_3 * pc_y[k] * sof_198[k];

        t_298[k] = f_7 * snf_149[k]
                   + f_3 * pc_y[k] * sof_199[k];

        t_299[k] = pb_y[k] * sng0_224[k]
                   - f_6 * pc_y[k] * sng1_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, snf_140, snf_200, \
                         snf_203, sod0_120, sod0_123, sod1_120, sod1_123, sof_200, \
                         sof_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_15 * snf_200[k]
                   + f_1 * sod0_120[k]
                   - f_2 * sod1_120[k]
                   + f_3 * pc_x[k] * sof_200[k];

        t_301[k] = f_3 * pc_y[k] * sof_200[k];

        t_302[k] = f_16 * snf_140[k]
                   + f_3 * pc_z[k] * sof_200[k];

        t_303[k] = f_15 * snf_203[k]
                   + f_4 * sod0_123[k]
                   - f_5 * sod1_123[k]
                   + f_3 * pc_x[k] * sof_203[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, snf_205, snf_206, snf_207, \
                         sod0_125, sod1_125, sof_202, sof_205, sof_206, \
                         sof_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_3 * pc_y[k] * sof_202[k];

        t_305[k] = f_15 * snf_205[k]
                   + f_4 * sod0_125[k]
                   - f_5 * sod1_125[k]
                   + f_3 * pc_x[k] * sof_205[k];

        t_306[k] = f_15 * snf_206[k]
                   + f_3 * pc_x[k] * sof_206[k];

        t_307[k] = f_15 * snf_207[k]
                   + f_3 * pc_x[k] * sof_207[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pc_x, pc_y, pc_z, snf_146, snf_208, \
                         snf_209, sod0_123, sod1_123, sof_206, sof_208, \
                         sof_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_15 * snf_208[k]
                   + f_3 * pc_x[k] * sof_208[k];

        t_309[k] = f_15 * snf_209[k]
                   + f_3 * pc_x[k] * sof_209[k];

        t_310[k] = f_1 * sod0_123[k]
                   - f_2 * sod1_123[k]
                   + f_3 * pc_y[k] * sof_206[k];

        t_311[k] = f_16 * snf_146[k]
                   + f_3 * pc_z[k] * sof_206[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, snf_149, snf_210, \
                         sod0_125, sod0_126, sod1_125, sod1_126, sof_208, sof_209, \
                         sof_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * sod0_125[k]
                   - f_5 * sod1_125[k]
                   + f_3 * pc_y[k] * sof_208[k];

        t_313[k] = f_3 * pc_y[k] * sof_209[k];

        t_314[k] = f_16 * snf_149[k]
                   + f_1 * sod0_125[k]
                   - f_2 * sod1_125[k]
                   + f_3 * pc_z[k] * sof_209[k];

        t_315[k] = f_16 * snf_210[k]
                   + f_1 * sod0_126[k]
                   - f_2 * sod1_126[k]
                   + f_3 * pc_x[k] * sof_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, snf_150, snf_152, \
                         snf_213, sod0_129, sod1_129, sof_210, sof_212, \
                         sof_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_15 * snf_150[k]
                   + f_3 * pc_y[k] * sof_210[k];

        t_317[k] = f_3 * pc_z[k] * sof_210[k];

        t_318[k] = f_16 * snf_213[k]
                   + f_4 * sod0_129[k]
                   - f_5 * sod1_129[k]
                   + f_3 * pc_x[k] * sof_213[k];

        t_319[k] = f_15 * snf_152[k]
                   + f_3 * pc_y[k] * sof_212[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, snf_215, snf_216, snf_217, snf_218, \
                         sod0_131, sod1_131, sof_215, sof_216, sof_217, \
                         sof_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_16 * snf_215[k]
                   + f_4 * sod0_131[k]
                   - f_5 * sod1_131[k]
                   + f_3 * pc_x[k] * sof_215[k];

        t_321[k] = f_16 * snf_216[k]
                   + f_3 * pc_x[k] * sof_216[k];

        t_322[k] = f_16 * snf_217[k]
                   + f_3 * pc_x[k] * sof_217[k];

        t_323[k] = f_16 * snf_218[k]
                   + f_3 * pc_x[k] * sof_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pc_x, pc_y, pc_z, snf_156, snf_219, sod0_129, \
                         sod1_129, sof_216, sof_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_16 * snf_219[k]
                   + f_3 * pc_x[k] * sof_219[k];

        t_325[k] = f_15 * snf_156[k]
                   + f_1 * sod0_129[k]
                   - f_2 * sod1_129[k]
                   + f_3 * pc_y[k] * sof_216[k];

        t_326[k] = f_3 * pc_z[k] * sof_216[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pb_z, pc_y, pc_z, sng0_225, snf_158, \
                         snf_159, sng1_225, sod0_131, sod1_131, sof_218, \
                         sof_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_15 * snf_158[k]
                   + f_4 * sod0_131[k]
                   - f_5 * sod1_131[k]
                   + f_3 * pc_y[k] * sof_218[k];

        t_328[k] = f_15 * snf_159[k]
                   + f_3 * pc_y[k] * sof_219[k];

        t_329[k] = f_1 * sod0_131[k]
                   - f_2 * sod1_131[k]
                   + f_3 * pc_z[k] * sof_219[k];

        t_330[k] = pb_z[k] * sng0_225[k]
                   - f_6 * pc_z[k] * sng1_225[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_z, pc_y, pc_z, sng0_228, snf_150, \
                         snf_160, snf_162, sng1_228, sof_220, sof_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_16 * snf_160[k]
                   + f_3 * pc_y[k] * sof_220[k];

        t_332[k] = f_7 * snf_150[k]
                   + f_3 * pc_z[k] * sof_220[k];

        t_333[k] = pb_z[k] * sng0_228[k]
                   - f_6 * pc_z[k] * sng1_228[k];

        t_334[k] = f_16 * snf_162[k]
                   + f_3 * pc_y[k] * sof_222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pc_x, snf_225, snf_226, snf_227, snf_228, \
                         sod0_137, sod1_137, sof_225, sof_226, sof_227, \
                         sof_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_16 * snf_225[k]
                   + f_4 * sod0_137[k]
                   - f_5 * sod1_137[k]
                   + f_3 * pc_x[k] * sof_225[k];

        t_336[k] = f_16 * snf_226[k]
                   + f_3 * pc_x[k] * sof_226[k];

        t_337[k] = f_16 * snf_227[k]
                   + f_3 * pc_x[k] * sof_227[k];

        t_338[k] = f_16 * snf_228[k]
                   + f_3 * pc_x[k] * sof_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pb_z, pc_x, pc_z, sng0_235, snf_156, snf_229, \
                         sng1_235, sof_226, sof_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_16 * snf_229[k]
                   + f_3 * pc_x[k] * sof_229[k];

        t_340[k] = pb_z[k] * sng0_235[k]
                   - f_6 * pc_z[k] * sng1_235[k];

        t_341[k] = f_7 * snf_156[k]
                   + f_3 * pc_z[k] * sof_226[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_y, pc_z, snf_159, snf_168, snf_169, sod0_137, \
                         sod1_137, sof_228, sof_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_16 * snf_168[k]
                   + f_4 * sod0_137[k]
                   - f_5 * sod1_137[k]
                   + f_3 * pc_y[k] * sof_228[k];

        t_343[k] = f_16 * snf_169[k]
                   + f_3 * pc_y[k] * sof_229[k];

        t_344[k] = f_7 * snf_159[k]
                   + f_1 * sod0_137[k]
                   - f_2 * sod1_137[k]
                   + f_3 * pc_z[k] * sof_229[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_y, pc_z, snf_160, snf_170, snf_230, \
                         sod0_138, sod1_138, sof_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_16 * snf_230[k]
                   + f_1 * sod0_138[k]
                   - f_2 * sod1_138[k]
                   + f_3 * pc_x[k] * sof_230[k];

        t_346[k] = f_14 * snf_170[k]
                   + f_3 * pc_y[k] * sof_230[k];

        t_347[k] = f_8 * snf_160[k]
                   + f_3 * pc_z[k] * sof_230[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, snf_172, snf_233, snf_235, sod0_141, \
                         sod0_143, sod1_141, sod1_143, sof_232, sof_233, \
                         sof_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_16 * snf_233[k]
                   + f_4 * sod0_141[k]
                   - f_5 * sod1_141[k]
                   + f_3 * pc_x[k] * sof_233[k];

        t_349[k] = f_14 * snf_172[k]
                   + f_3 * pc_y[k] * sof_232[k];

        t_350[k] = f_16 * snf_235[k]
                   + f_4 * sod0_143[k]
                   - f_5 * sod1_143[k]
                   + f_3 * pc_x[k] * sof_235[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, pc_x, snf_236, snf_237, snf_238, snf_239, \
                         sof_236, sof_237, sof_238, sof_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_16 * snf_236[k]
                   + f_3 * pc_x[k] * sof_236[k];

        t_352[k] = f_16 * snf_237[k]
                   + f_3 * pc_x[k] * sof_237[k];

        t_353[k] = f_16 * snf_238[k]
                   + f_3 * pc_x[k] * sof_238[k];

        t_354[k] = f_16 * snf_239[k]
                   + f_3 * pc_x[k] * sof_239[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, pc_y, pc_z, snf_166, snf_176, snf_178, sod0_141, \
                         sod0_143, sod1_141, sod1_143, sof_236, \
                         sof_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_14 * snf_176[k]
                   + f_1 * sod0_141[k]
                   - f_2 * sod1_141[k]
                   + f_3 * pc_y[k] * sof_236[k];

        t_356[k] = f_8 * snf_166[k]
                   + f_3 * pc_z[k] * sof_236[k];

        t_357[k] = f_14 * snf_178[k]
                   + f_4 * sod0_143[k]
                   - f_5 * sod1_143[k]
                   + f_3 * pc_y[k] * sof_238[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pc_x, pc_y, pc_z, snf_169, snf_179, snf_240, \
                         sod0_143, sod0_144, sod1_143, sod1_144, sof_239, \
                         sof_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_14 * snf_179[k]
                   + f_3 * pc_y[k] * sof_239[k];

        t_359[k] = f_8 * snf_169[k]
                   + f_1 * sod0_143[k]
                   - f_2 * sod1_143[k]
                   + f_3 * pc_z[k] * sof_239[k];

        t_360[k] = f_16 * snf_240[k]
                   + f_1 * sod0_144[k]
                   - f_2 * sod1_144[k]
                   + f_3 * pc_x[k] * sof_240[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
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
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_300 = buffer.data(sng0 + 300);
    const auto *sng0_303 = buffer.data(sng0 + 303);
    const auto *sng0_305 = buffer.data(sng0 + 305);
    const auto *sng0_314 = buffer.data(sng0 + 314);
    const auto *sng0_315 = buffer.data(sng0 + 315);
    const auto *sng0_318 = buffer.data(sng0 + 318);
    const auto *sng0_325 = buffer.data(sng0 + 325);

    const auto *snf_170 = buffer.data(snf + 170);
    const auto *snf_176 = buffer.data(snf + 176);
    const auto *snf_179 = buffer.data(snf + 179);
    const auto *snf_180 = buffer.data(snf + 180);
    const auto *snf_182 = buffer.data(snf + 182);
    const auto *snf_186 = buffer.data(snf + 186);
    const auto *snf_188 = buffer.data(snf + 188);
    const auto *snf_189 = buffer.data(snf + 189);
    const auto *snf_190 = buffer.data(snf + 190);
    const auto *snf_192 = buffer.data(snf + 192);
    const auto *snf_196 = buffer.data(snf + 196);
    const auto *snf_198 = buffer.data(snf + 198);
    const auto *snf_199 = buffer.data(snf + 199);
    const auto *snf_200 = buffer.data(snf + 200);
    const auto *snf_201 = buffer.data(snf + 201);
    const auto *snf_202 = buffer.data(snf + 202);
    const auto *snf_206 = buffer.data(snf + 206);
    const auto *snf_208 = buffer.data(snf + 208);
    const auto *snf_209 = buffer.data(snf + 209);
    const auto *snf_210 = buffer.data(snf + 210);
    const auto *snf_212 = buffer.data(snf + 212);
    const auto *snf_216 = buffer.data(snf + 216);
    const auto *snf_218 = buffer.data(snf + 218);
    const auto *snf_219 = buffer.data(snf + 219);
    const auto *snf_220 = buffer.data(snf + 220);
    const auto *snf_222 = buffer.data(snf + 222);
    const auto *snf_226 = buffer.data(snf + 226);
    const auto *snf_228 = buffer.data(snf + 228);
    const auto *snf_229 = buffer.data(snf + 229);
    const auto *snf_230 = buffer.data(snf + 230);
    const auto *snf_232 = buffer.data(snf + 232);
    const auto *snf_236 = buffer.data(snf + 236);
    const auto *snf_238 = buffer.data(snf + 238);
    const auto *snf_239 = buffer.data(snf + 239);
    const auto *snf_240 = buffer.data(snf + 240);
    const auto *snf_242 = buffer.data(snf + 242);
    const auto *snf_243 = buffer.data(snf + 243);
    const auto *snf_245 = buffer.data(snf + 245);
    const auto *snf_246 = buffer.data(snf + 246);
    const auto *snf_247 = buffer.data(snf + 247);
    const auto *snf_248 = buffer.data(snf + 248);
    const auto *snf_249 = buffer.data(snf + 249);
    const auto *snf_250 = buffer.data(snf + 250);
    const auto *snf_253 = buffer.data(snf + 253);
    const auto *snf_255 = buffer.data(snf + 255);
    const auto *snf_256 = buffer.data(snf + 256);
    const auto *snf_257 = buffer.data(snf + 257);
    const auto *snf_258 = buffer.data(snf + 258);
    const auto *snf_259 = buffer.data(snf + 259);
    const auto *snf_266 = buffer.data(snf + 266);
    const auto *snf_267 = buffer.data(snf + 267);
    const auto *snf_268 = buffer.data(snf + 268);
    const auto *snf_269 = buffer.data(snf + 269);
    const auto *snf_270 = buffer.data(snf + 270);
    const auto *snf_273 = buffer.data(snf + 273);
    const auto *snf_275 = buffer.data(snf + 275);
    const auto *snf_276 = buffer.data(snf + 276);
    const auto *snf_277 = buffer.data(snf + 277);
    const auto *snf_278 = buffer.data(snf + 278);
    const auto *snf_279 = buffer.data(snf + 279);
    const auto *snf_280 = buffer.data(snf + 280);
    const auto *snf_283 = buffer.data(snf + 283);
    const auto *snf_285 = buffer.data(snf + 285);
    const auto *snf_286 = buffer.data(snf + 286);
    const auto *snf_287 = buffer.data(snf + 287);
    const auto *snf_288 = buffer.data(snf + 288);
    const auto *snf_289 = buffer.data(snf + 289);
    const auto *snf_295 = buffer.data(snf + 295);
    const auto *snf_296 = buffer.data(snf + 296);
    const auto *snf_297 = buffer.data(snf + 297);
    const auto *snf_298 = buffer.data(snf + 298);
    const auto *snf_299 = buffer.data(snf + 299);
    const auto *snf_300 = buffer.data(snf + 300);
    const auto *snf_303 = buffer.data(snf + 303);
    const auto *snf_305 = buffer.data(snf + 305);
    const auto *snf_306 = buffer.data(snf + 306);
    const auto *snf_307 = buffer.data(snf + 307);
    const auto *snf_308 = buffer.data(snf + 308);
    const auto *snf_309 = buffer.data(snf + 309);
    const auto *snf_310 = buffer.data(snf + 310);
    const auto *snf_313 = buffer.data(snf + 313);
    const auto *snf_315 = buffer.data(snf + 315);
    const auto *snf_316 = buffer.data(snf + 316);
    const auto *snf_317 = buffer.data(snf + 317);
    const auto *snf_318 = buffer.data(snf + 318);
    const auto *snf_319 = buffer.data(snf + 319);

    const auto *sng1_300 = buffer.data(sng1 + 300);
    const auto *sng1_303 = buffer.data(sng1 + 303);
    const auto *sng1_305 = buffer.data(sng1 + 305);
    const auto *sng1_314 = buffer.data(sng1 + 314);
    const auto *sng1_315 = buffer.data(sng1 + 315);
    const auto *sng1_318 = buffer.data(sng1 + 318);
    const auto *sng1_325 = buffer.data(sng1 + 325);

    const auto *sod0_147 = buffer.data(sod0 + 147);
    const auto *sod0_149 = buffer.data(sod0 + 149);
    const auto *sod0_150 = buffer.data(sod0 + 150);
    const auto *sod0_153 = buffer.data(sod0 + 153);
    const auto *sod0_155 = buffer.data(sod0 + 155);
    const auto *sod0_159 = buffer.data(sod0 + 159);
    const auto *sod0_161 = buffer.data(sod0 + 161);
    const auto *sod0_162 = buffer.data(sod0 + 162);
    const auto *sod0_165 = buffer.data(sod0 + 165);
    const auto *sod0_167 = buffer.data(sod0 + 167);
    const auto *sod0_168 = buffer.data(sod0 + 168);
    const auto *sod0_171 = buffer.data(sod0 + 171);
    const auto *sod0_173 = buffer.data(sod0 + 173);
    const auto *sod0_179 = buffer.data(sod0 + 179);
    const auto *sod0_180 = buffer.data(sod0 + 180);
    const auto *sod0_183 = buffer.data(sod0 + 183);
    const auto *sod0_185 = buffer.data(sod0 + 185);
    const auto *sod0_186 = buffer.data(sod0 + 186);
    const auto *sod0_189 = buffer.data(sod0 + 189);
    const auto *sod0_191 = buffer.data(sod0 + 191);

    const auto *sod1_147 = buffer.data(sod1 + 147);
    const auto *sod1_149 = buffer.data(sod1 + 149);
    const auto *sod1_150 = buffer.data(sod1 + 150);
    const auto *sod1_153 = buffer.data(sod1 + 153);
    const auto *sod1_155 = buffer.data(sod1 + 155);
    const auto *sod1_159 = buffer.data(sod1 + 159);
    const auto *sod1_161 = buffer.data(sod1 + 161);
    const auto *sod1_162 = buffer.data(sod1 + 162);
    const auto *sod1_165 = buffer.data(sod1 + 165);
    const auto *sod1_167 = buffer.data(sod1 + 167);
    const auto *sod1_168 = buffer.data(sod1 + 168);
    const auto *sod1_171 = buffer.data(sod1 + 171);
    const auto *sod1_173 = buffer.data(sod1 + 173);
    const auto *sod1_179 = buffer.data(sod1 + 179);
    const auto *sod1_180 = buffer.data(sod1 + 180);
    const auto *sod1_183 = buffer.data(sod1 + 183);
    const auto *sod1_185 = buffer.data(sod1 + 185);
    const auto *sod1_186 = buffer.data(sod1 + 186);
    const auto *sod1_189 = buffer.data(sod1 + 189);
    const auto *sod1_191 = buffer.data(sod1 + 191);

    const auto *sof_240 = buffer.data(sof + 240);
    const auto *sof_242 = buffer.data(sof + 242);
    const auto *sof_243 = buffer.data(sof + 243);
    const auto *sof_245 = buffer.data(sof + 245);
    const auto *sof_246 = buffer.data(sof + 246);
    const auto *sof_247 = buffer.data(sof + 247);
    const auto *sof_248 = buffer.data(sof + 248);
    const auto *sof_249 = buffer.data(sof + 249);
    const auto *sof_250 = buffer.data(sof + 250);
    const auto *sof_252 = buffer.data(sof + 252);
    const auto *sof_253 = buffer.data(sof + 253);
    const auto *sof_255 = buffer.data(sof + 255);
    const auto *sof_256 = buffer.data(sof + 256);
    const auto *sof_257 = buffer.data(sof + 257);
    const auto *sof_258 = buffer.data(sof + 258);
    const auto *sof_259 = buffer.data(sof + 259);
    const auto *sof_260 = buffer.data(sof + 260);
    const auto *sof_262 = buffer.data(sof + 262);
    const auto *sof_266 = buffer.data(sof + 266);
    const auto *sof_267 = buffer.data(sof + 267);
    const auto *sof_268 = buffer.data(sof + 268);
    const auto *sof_269 = buffer.data(sof + 269);
    const auto *sof_270 = buffer.data(sof + 270);
    const auto *sof_272 = buffer.data(sof + 272);
    const auto *sof_273 = buffer.data(sof + 273);
    const auto *sof_275 = buffer.data(sof + 275);
    const auto *sof_276 = buffer.data(sof + 276);
    const auto *sof_277 = buffer.data(sof + 277);
    const auto *sof_278 = buffer.data(sof + 278);
    const auto *sof_279 = buffer.data(sof + 279);
    const auto *sof_280 = buffer.data(sof + 280);
    const auto *sof_282 = buffer.data(sof + 282);
    const auto *sof_283 = buffer.data(sof + 283);
    const auto *sof_285 = buffer.data(sof + 285);
    const auto *sof_286 = buffer.data(sof + 286);
    const auto *sof_287 = buffer.data(sof + 287);
    const auto *sof_288 = buffer.data(sof + 288);
    const auto *sof_289 = buffer.data(sof + 289);
    const auto *sof_290 = buffer.data(sof + 290);
    const auto *sof_292 = buffer.data(sof + 292);
    const auto *sof_295 = buffer.data(sof + 295);
    const auto *sof_296 = buffer.data(sof + 296);
    const auto *sof_297 = buffer.data(sof + 297);
    const auto *sof_298 = buffer.data(sof + 298);
    const auto *sof_299 = buffer.data(sof + 299);
    const auto *sof_300 = buffer.data(sof + 300);
    const auto *sof_302 = buffer.data(sof + 302);
    const auto *sof_303 = buffer.data(sof + 303);
    const auto *sof_305 = buffer.data(sof + 305);
    const auto *sof_306 = buffer.data(sof + 306);
    const auto *sof_307 = buffer.data(sof + 307);
    const auto *sof_308 = buffer.data(sof + 308);
    const auto *sof_309 = buffer.data(sof + 309);
    const auto *sof_310 = buffer.data(sof + 310);
    const auto *sof_312 = buffer.data(sof + 312);
    const auto *sof_313 = buffer.data(sof + 313);
    const auto *sof_315 = buffer.data(sof + 315);
    const auto *sof_316 = buffer.data(sof + 316);
    const auto *sof_317 = buffer.data(sof + 317);
    const auto *sof_318 = buffer.data(sof + 318);
    const auto *sof_319 = buffer.data(sof + 319);

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pc_x, pc_y, pc_z, snf_170, snf_180, \
                         snf_182, snf_243, sod0_147, sod1_147, sof_240, sof_242, \
                         sof_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_12 * snf_180[k]
                   + f_3 * pc_y[k] * sof_240[k];

        t_362[k] = f_12 * snf_170[k]
                   + f_3 * pc_z[k] * sof_240[k];

        t_363[k] = f_16 * snf_243[k]
                   + f_4 * sod0_147[k]
                   - f_5 * sod1_147[k]
                   + f_3 * pc_x[k] * sof_243[k];

        t_364[k] = f_12 * snf_182[k]
                   + f_3 * pc_y[k] * sof_242[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, snf_245, snf_246, snf_247, snf_248, \
                         sod0_149, sod1_149, sof_245, sof_246, sof_247, \
                         sof_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_16 * snf_245[k]
                   + f_4 * sod0_149[k]
                   - f_5 * sod1_149[k]
                   + f_3 * pc_x[k] * sof_245[k];

        t_366[k] = f_16 * snf_246[k]
                   + f_3 * pc_x[k] * sof_246[k];

        t_367[k] = f_16 * snf_247[k]
                   + f_3 * pc_x[k] * sof_247[k];

        t_368[k] = f_16 * snf_248[k]
                   + f_3 * pc_x[k] * sof_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pc_x, pc_y, pc_z, snf_176, snf_186, snf_249, \
                         sod0_147, sod1_147, sof_246, sof_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_16 * snf_249[k]
                   + f_3 * pc_x[k] * sof_249[k];

        t_370[k] = f_12 * snf_186[k]
                   + f_1 * sod0_147[k]
                   - f_2 * sod1_147[k]
                   + f_3 * pc_y[k] * sof_246[k];

        t_371[k] = f_12 * snf_176[k]
                   + f_3 * pc_z[k] * sof_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, snf_179, snf_188, snf_189, sod0_149, \
                         sod1_149, sof_248, sof_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_12 * snf_188[k]
                   + f_4 * sod0_149[k]
                   - f_5 * sod1_149[k]
                   + f_3 * pc_y[k] * sof_248[k];

        t_373[k] = f_12 * snf_189[k]
                   + f_3 * pc_y[k] * sof_249[k];

        t_374[k] = f_12 * snf_179[k]
                   + f_1 * sod0_149[k]
                   - f_2 * sod1_149[k]
                   + f_3 * pc_z[k] * sof_249[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_y, pc_z, snf_180, snf_190, snf_250, \
                         sod0_150, sod1_150, sof_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_16 * snf_250[k]
                   + f_1 * sod0_150[k]
                   - f_2 * sod1_150[k]
                   + f_3 * pc_x[k] * sof_250[k];

        t_376[k] = f_8 * snf_190[k]
                   + f_3 * pc_y[k] * sof_250[k];

        t_377[k] = f_14 * snf_180[k]
                   + f_3 * pc_z[k] * sof_250[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pc_x, pc_y, snf_192, snf_253, snf_255, sod0_153, \
                         sod0_155, sod1_153, sod1_155, sof_252, sof_253, \
                         sof_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_16 * snf_253[k]
                   + f_4 * sod0_153[k]
                   - f_5 * sod1_153[k]
                   + f_3 * pc_x[k] * sof_253[k];

        t_379[k] = f_8 * snf_192[k]
                   + f_3 * pc_y[k] * sof_252[k];

        t_380[k] = f_16 * snf_255[k]
                   + f_4 * sod0_155[k]
                   - f_5 * sod1_155[k]
                   + f_3 * pc_x[k] * sof_255[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, snf_256, snf_257, snf_258, snf_259, \
                         sof_256, sof_257, sof_258, sof_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_16 * snf_256[k]
                   + f_3 * pc_x[k] * sof_256[k];

        t_382[k] = f_16 * snf_257[k]
                   + f_3 * pc_x[k] * sof_257[k];

        t_383[k] = f_16 * snf_258[k]
                   + f_3 * pc_x[k] * sof_258[k];

        t_384[k] = f_16 * snf_259[k]
                   + f_3 * pc_x[k] * sof_259[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, snf_186, snf_196, snf_198, sod0_153, \
                         sod0_155, sod1_153, sod1_155, sof_256, \
                         sof_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_8 * snf_196[k]
                   + f_1 * sod0_153[k]
                   - f_2 * sod1_153[k]
                   + f_3 * pc_y[k] * sof_256[k];

        t_386[k] = f_14 * snf_186[k]
                   + f_3 * pc_z[k] * sof_256[k];

        t_387[k] = f_8 * snf_198[k]
                   + f_4 * sod0_155[k]
                   - f_5 * sod1_155[k]
                   + f_3 * pc_y[k] * sof_258[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_y, pc_y, pc_z, sng0_300, snf_189, \
                         snf_199, snf_200, sng1_300, sod0_155, sod1_155, sof_259, \
                         sof_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_8 * snf_199[k]
                   + f_3 * pc_y[k] * sof_259[k];

        t_389[k] = f_14 * snf_189[k]
                   + f_1 * sod0_155[k]
                   - f_2 * sod1_155[k]
                   + f_3 * pc_z[k] * sof_259[k];

        t_390[k] = pb_y[k] * sng0_300[k]
                   - f_6 * pc_y[k] * sng1_300[k];

        t_391[k] = f_7 * snf_200[k]
                   + f_3 * pc_y[k] * sof_260[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pb_y, pc_y, pc_z, sng0_303, sng0_305, \
                         snf_190, snf_201, snf_202, sng1_303, sng1_305, sof_260, \
                         sof_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_16 * snf_190[k]
                   + f_3 * pc_z[k] * sof_260[k];

        t_393[k] = pb_y[k] * sng0_303[k]
                   + f_8 * snf_201[k]
                   - f_6 * pc_y[k] * sng1_303[k];

        t_394[k] = f_7 * snf_202[k]
                   + f_3 * pc_y[k] * sof_262[k];

        t_395[k] = pb_y[k] * sng0_305[k]
                   - f_6 * pc_y[k] * sng1_305[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pc_x, snf_266, snf_267, snf_268, snf_269, \
                         sof_266, sof_267, sof_268, sof_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_16 * snf_266[k]
                   + f_3 * pc_x[k] * sof_266[k];

        t_397[k] = f_16 * snf_267[k]
                   + f_3 * pc_x[k] * sof_267[k];

        t_398[k] = f_16 * snf_268[k]
                   + f_3 * pc_x[k] * sof_268[k];

        t_399[k] = f_16 * snf_269[k]
                   + f_3 * pc_x[k] * sof_269[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pc_y, pc_z, snf_196, snf_206, snf_208, sod0_159, \
                         sod0_161, sod1_159, sod1_161, sof_266, \
                         sof_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_7 * snf_206[k]
                   + f_1 * sod0_159[k]
                   - f_2 * sod1_159[k]
                   + f_3 * pc_y[k] * sof_266[k];

        t_401[k] = f_16 * snf_196[k]
                   + f_3 * pc_z[k] * sof_266[k];

        t_402[k] = f_7 * snf_208[k]
                   + f_4 * sod0_161[k]
                   - f_5 * sod1_161[k]
                   + f_3 * pc_y[k] * sof_268[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pb_y, pc_x, pc_y, sng0_314, snf_209, \
                         snf_270, sng1_314, sod0_162, sod1_162, sof_269, \
                         sof_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_7 * snf_209[k]
                   + f_3 * pc_y[k] * sof_269[k];

        t_404[k] = pb_y[k] * sng0_314[k]
                   - f_6 * pc_y[k] * sng1_314[k];

        t_405[k] = f_16 * snf_270[k]
                   + f_1 * sod0_162[k]
                   - f_2 * sod1_162[k]
                   + f_3 * pc_x[k] * sof_270[k];

        t_406[k] = f_3 * pc_y[k] * sof_270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pc_x, pc_y, pc_z, snf_200, snf_273, sod0_165, \
                         sod1_165, sof_270, sof_272, sof_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_15 * snf_200[k]
                   + f_3 * pc_z[k] * sof_270[k];

        t_408[k] = f_16 * snf_273[k]
                   + f_4 * sod0_165[k]
                   - f_5 * sod1_165[k]
                   + f_3 * pc_x[k] * sof_273[k];

        t_409[k] = f_3 * pc_y[k] * sof_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, snf_275, snf_276, snf_277, snf_278, \
                         sod0_167, sod1_167, sof_275, sof_276, sof_277, \
                         sof_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_16 * snf_275[k]
                   + f_4 * sod0_167[k]
                   - f_5 * sod1_167[k]
                   + f_3 * pc_x[k] * sof_275[k];

        t_411[k] = f_16 * snf_276[k]
                   + f_3 * pc_x[k] * sof_276[k];

        t_412[k] = f_16 * snf_277[k]
                   + f_3 * pc_x[k] * sof_277[k];

        t_413[k] = f_16 * snf_278[k]
                   + f_3 * pc_x[k] * sof_278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, snf_206, snf_279, \
                         sod0_165, sod0_167, sod1_165, sod1_167, sof_276, sof_278, \
                         sof_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_16 * snf_279[k]
                   + f_3 * pc_x[k] * sof_279[k];

        t_415[k] = f_1 * sod0_165[k]
                   - f_2 * sod1_165[k]
                   + f_3 * pc_y[k] * sof_276[k];

        t_416[k] = f_15 * snf_206[k]
                   + f_3 * pc_z[k] * sof_276[k];

        t_417[k] = f_4 * sod0_167[k]
                   - f_5 * sod1_167[k]
                   + f_3 * pc_y[k] * sof_278[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pc_x, pc_y, pc_z, snf_209, snf_210, \
                         snf_280, sod0_167, sod0_168, sod1_167, sod1_168, sof_279, \
                         sof_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_3 * pc_y[k] * sof_279[k];

        t_419[k] = f_15 * snf_209[k]
                   + f_1 * sod0_167[k]
                   - f_2 * sod1_167[k]
                   + f_3 * pc_z[k] * sof_279[k];

        t_420[k] = f_14 * snf_280[k]
                   + f_1 * sod0_168[k]
                   - f_2 * sod1_168[k]
                   + f_3 * pc_x[k] * sof_280[k];

        t_421[k] = f_13 * snf_210[k]
                   + f_3 * pc_y[k] * sof_280[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, pc_x, pc_y, pc_z, snf_212, snf_283, sod0_171, \
                         sod1_171, sof_280, sof_282, sof_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_3 * pc_z[k] * sof_280[k];

        t_423[k] = f_14 * snf_283[k]
                   + f_4 * sod0_171[k]
                   - f_5 * sod1_171[k]
                   + f_3 * pc_x[k] * sof_283[k];

        t_424[k] = f_13 * snf_212[k]
                   + f_3 * pc_y[k] * sof_282[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, snf_285, snf_286, snf_287, snf_288, \
                         sod0_173, sod1_173, sof_285, sof_286, sof_287, \
                         sof_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_14 * snf_285[k]
                   + f_4 * sod0_173[k]
                   - f_5 * sod1_173[k]
                   + f_3 * pc_x[k] * sof_285[k];

        t_426[k] = f_14 * snf_286[k]
                   + f_3 * pc_x[k] * sof_286[k];

        t_427[k] = f_14 * snf_287[k]
                   + f_3 * pc_x[k] * sof_287[k];

        t_428[k] = f_14 * snf_288[k]
                   + f_3 * pc_x[k] * sof_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_x, pc_y, pc_z, snf_216, snf_289, sod0_171, \
                         sod1_171, sof_286, sof_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_14 * snf_289[k]
                   + f_3 * pc_x[k] * sof_289[k];

        t_430[k] = f_13 * snf_216[k]
                   + f_1 * sod0_171[k]
                   - f_2 * sod1_171[k]
                   + f_3 * pc_y[k] * sof_286[k];

        t_431[k] = f_3 * pc_z[k] * sof_286[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pb_z, pc_y, pc_z, sng0_315, snf_218, \
                         snf_219, sng1_315, sod0_173, sod1_173, sof_288, \
                         sof_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_13 * snf_218[k]
                   + f_4 * sod0_173[k]
                   - f_5 * sod1_173[k]
                   + f_3 * pc_y[k] * sof_288[k];

        t_433[k] = f_13 * snf_219[k]
                   + f_3 * pc_y[k] * sof_289[k];

        t_434[k] = f_1 * sod0_173[k]
                   - f_2 * sod1_173[k]
                   + f_3 * pc_z[k] * sof_289[k];

        t_435[k] = pb_z[k] * sng0_315[k]
                   - f_6 * pc_z[k] * sng1_315[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pb_z, pc_y, pc_z, sng0_318, snf_210, \
                         snf_220, snf_222, sng1_318, sof_290, sof_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_15 * snf_220[k]
                   + f_3 * pc_y[k] * sof_290[k];

        t_437[k] = f_7 * snf_210[k]
                   + f_3 * pc_z[k] * sof_290[k];

        t_438[k] = pb_z[k] * sng0_318[k]
                   - f_6 * pc_z[k] * sng1_318[k];

        t_439[k] = f_15 * snf_222[k]
                   + f_3 * pc_y[k] * sof_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pc_x, snf_295, snf_296, snf_297, snf_298, \
                         sod0_179, sod1_179, sof_295, sof_296, sof_297, \
                         sof_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_14 * snf_295[k]
                   + f_4 * sod0_179[k]
                   - f_5 * sod1_179[k]
                   + f_3 * pc_x[k] * sof_295[k];

        t_441[k] = f_14 * snf_296[k]
                   + f_3 * pc_x[k] * sof_296[k];

        t_442[k] = f_14 * snf_297[k]
                   + f_3 * pc_x[k] * sof_297[k];

        t_443[k] = f_14 * snf_298[k]
                   + f_3 * pc_x[k] * sof_298[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pb_z, pc_x, pc_z, sng0_325, snf_216, snf_299, \
                         sng1_325, sof_296, sof_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_14 * snf_299[k]
                   + f_3 * pc_x[k] * sof_299[k];

        t_445[k] = pb_z[k] * sng0_325[k]
                   - f_6 * pc_z[k] * sng1_325[k];

        t_446[k] = f_7 * snf_216[k]
                   + f_3 * pc_z[k] * sof_296[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pc_y, pc_z, snf_219, snf_228, snf_229, sod0_179, \
                         sod1_179, sof_298, sof_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_15 * snf_228[k]
                   + f_4 * sod0_179[k]
                   - f_5 * sod1_179[k]
                   + f_3 * pc_y[k] * sof_298[k];

        t_448[k] = f_15 * snf_229[k]
                   + f_3 * pc_y[k] * sof_299[k];

        t_449[k] = f_7 * snf_219[k]
                   + f_1 * sod0_179[k]
                   - f_2 * sod1_179[k]
                   + f_3 * pc_z[k] * sof_299[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, pc_x, pc_y, pc_z, snf_220, snf_230, snf_300, \
                         sod0_180, sod1_180, sof_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_14 * snf_300[k]
                   + f_1 * sod0_180[k]
                   - f_2 * sod1_180[k]
                   + f_3 * pc_x[k] * sof_300[k];

        t_451[k] = f_16 * snf_230[k]
                   + f_3 * pc_y[k] * sof_300[k];

        t_452[k] = f_8 * snf_220[k]
                   + f_3 * pc_z[k] * sof_300[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pc_x, pc_y, snf_232, snf_303, snf_305, sod0_183, \
                         sod0_185, sod1_183, sod1_185, sof_302, sof_303, \
                         sof_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_14 * snf_303[k]
                   + f_4 * sod0_183[k]
                   - f_5 * sod1_183[k]
                   + f_3 * pc_x[k] * sof_303[k];

        t_454[k] = f_16 * snf_232[k]
                   + f_3 * pc_y[k] * sof_302[k];

        t_455[k] = f_14 * snf_305[k]
                   + f_4 * sod0_185[k]
                   - f_5 * sod1_185[k]
                   + f_3 * pc_x[k] * sof_305[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pc_x, snf_306, snf_307, snf_308, snf_309, \
                         sof_306, sof_307, sof_308, sof_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * snf_306[k]
                   + f_3 * pc_x[k] * sof_306[k];

        t_457[k] = f_14 * snf_307[k]
                   + f_3 * pc_x[k] * sof_307[k];

        t_458[k] = f_14 * snf_308[k]
                   + f_3 * pc_x[k] * sof_308[k];

        t_459[k] = f_14 * snf_309[k]
                   + f_3 * pc_x[k] * sof_309[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, snf_226, snf_236, snf_238, sod0_183, \
                         sod0_185, sod1_183, sod1_185, sof_306, \
                         sof_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_16 * snf_236[k]
                   + f_1 * sod0_183[k]
                   - f_2 * sod1_183[k]
                   + f_3 * pc_y[k] * sof_306[k];

        t_461[k] = f_8 * snf_226[k]
                   + f_3 * pc_z[k] * sof_306[k];

        t_462[k] = f_16 * snf_238[k]
                   + f_4 * sod0_185[k]
                   - f_5 * sod1_185[k]
                   + f_3 * pc_y[k] * sof_308[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_x, pc_y, pc_z, snf_229, snf_239, snf_310, \
                         sod0_185, sod0_186, sod1_185, sod1_186, sof_309, \
                         sof_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_16 * snf_239[k]
                   + f_3 * pc_y[k] * sof_309[k];

        t_464[k] = f_8 * snf_229[k]
                   + f_1 * sod0_185[k]
                   - f_2 * sod1_185[k]
                   + f_3 * pc_z[k] * sof_309[k];

        t_465[k] = f_14 * snf_310[k]
                   + f_1 * sod0_186[k]
                   - f_2 * sod1_186[k]
                   + f_3 * pc_x[k] * sof_310[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pc_x, pc_y, pc_z, snf_230, snf_240, \
                         snf_242, snf_313, sod0_189, sod1_189, sof_310, sof_312, \
                         sof_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_14 * snf_240[k]
                   + f_3 * pc_y[k] * sof_310[k];

        t_467[k] = f_12 * snf_230[k]
                   + f_3 * pc_z[k] * sof_310[k];

        t_468[k] = f_14 * snf_313[k]
                   + f_4 * sod0_189[k]
                   - f_5 * sod1_189[k]
                   + f_3 * pc_x[k] * sof_313[k];

        t_469[k] = f_14 * snf_242[k]
                   + f_3 * pc_y[k] * sof_312[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, snf_315, snf_316, snf_317, snf_318, \
                         sod0_191, sod1_191, sof_315, sof_316, sof_317, \
                         sof_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_14 * snf_315[k]
                   + f_4 * sod0_191[k]
                   - f_5 * sod1_191[k]
                   + f_3 * pc_x[k] * sof_315[k];

        t_471[k] = f_14 * snf_316[k]
                   + f_3 * pc_x[k] * sof_316[k];

        t_472[k] = f_14 * snf_317[k]
                   + f_3 * pc_x[k] * sof_317[k];

        t_473[k] = f_14 * snf_318[k]
                   + f_3 * pc_x[k] * sof_318[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, pc_z, snf_236, snf_246, snf_319, \
                         sod0_189, sod1_189, sof_316, sof_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_14 * snf_319[k]
                   + f_3 * pc_x[k] * sof_319[k];

        t_475[k] = f_14 * snf_246[k]
                   + f_1 * sod0_189[k]
                   - f_2 * sod1_189[k]
                   + f_3 * pc_y[k] * sof_316[k];

        t_476[k] = f_12 * snf_236[k]
                   + f_3 * pc_z[k] * sof_316[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
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
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_405 = buffer.data(sng0 + 405);
    const auto *sng0_408 = buffer.data(sng0 + 408);
    const auto *sng0_410 = buffer.data(sng0 + 410);
    const auto *sng0_419 = buffer.data(sng0 + 419);
    const auto *sng0_420 = buffer.data(sng0 + 420);
    const auto *sng0_423 = buffer.data(sng0 + 423);
    const auto *sng0_430 = buffer.data(sng0 + 430);

    const auto *snf_239 = buffer.data(snf + 239);
    const auto *snf_240 = buffer.data(snf + 240);
    const auto *snf_246 = buffer.data(snf + 246);
    const auto *snf_248 = buffer.data(snf + 248);
    const auto *snf_249 = buffer.data(snf + 249);
    const auto *snf_250 = buffer.data(snf + 250);
    const auto *snf_252 = buffer.data(snf + 252);
    const auto *snf_256 = buffer.data(snf + 256);
    const auto *snf_258 = buffer.data(snf + 258);
    const auto *snf_259 = buffer.data(snf + 259);
    const auto *snf_260 = buffer.data(snf + 260);
    const auto *snf_262 = buffer.data(snf + 262);
    const auto *snf_266 = buffer.data(snf + 266);
    const auto *snf_268 = buffer.data(snf + 268);
    const auto *snf_269 = buffer.data(snf + 269);
    const auto *snf_270 = buffer.data(snf + 270);
    const auto *snf_271 = buffer.data(snf + 271);
    const auto *snf_272 = buffer.data(snf + 272);
    const auto *snf_276 = buffer.data(snf + 276);
    const auto *snf_278 = buffer.data(snf + 278);
    const auto *snf_279 = buffer.data(snf + 279);
    const auto *snf_280 = buffer.data(snf + 280);
    const auto *snf_282 = buffer.data(snf + 282);
    const auto *snf_286 = buffer.data(snf + 286);
    const auto *snf_288 = buffer.data(snf + 288);
    const auto *snf_289 = buffer.data(snf + 289);
    const auto *snf_290 = buffer.data(snf + 290);
    const auto *snf_292 = buffer.data(snf + 292);
    const auto *snf_296 = buffer.data(snf + 296);
    const auto *snf_298 = buffer.data(snf + 298);
    const auto *snf_299 = buffer.data(snf + 299);
    const auto *snf_300 = buffer.data(snf + 300);
    const auto *snf_302 = buffer.data(snf + 302);
    const auto *snf_306 = buffer.data(snf + 306);
    const auto *snf_308 = buffer.data(snf + 308);
    const auto *snf_309 = buffer.data(snf + 309);
    const auto *snf_310 = buffer.data(snf + 310);
    const auto *snf_312 = buffer.data(snf + 312);
    const auto *snf_320 = buffer.data(snf + 320);
    const auto *snf_323 = buffer.data(snf + 323);
    const auto *snf_325 = buffer.data(snf + 325);
    const auto *snf_326 = buffer.data(snf + 326);
    const auto *snf_327 = buffer.data(snf + 327);
    const auto *snf_328 = buffer.data(snf + 328);
    const auto *snf_329 = buffer.data(snf + 329);
    const auto *snf_330 = buffer.data(snf + 330);
    const auto *snf_333 = buffer.data(snf + 333);
    const auto *snf_335 = buffer.data(snf + 335);
    const auto *snf_336 = buffer.data(snf + 336);
    const auto *snf_337 = buffer.data(snf + 337);
    const auto *snf_338 = buffer.data(snf + 338);
    const auto *snf_339 = buffer.data(snf + 339);
    const auto *snf_346 = buffer.data(snf + 346);
    const auto *snf_347 = buffer.data(snf + 347);
    const auto *snf_348 = buffer.data(snf + 348);
    const auto *snf_349 = buffer.data(snf + 349);
    const auto *snf_350 = buffer.data(snf + 350);
    const auto *snf_353 = buffer.data(snf + 353);
    const auto *snf_355 = buffer.data(snf + 355);
    const auto *snf_356 = buffer.data(snf + 356);
    const auto *snf_357 = buffer.data(snf + 357);
    const auto *snf_358 = buffer.data(snf + 358);
    const auto *snf_359 = buffer.data(snf + 359);
    const auto *snf_360 = buffer.data(snf + 360);
    const auto *snf_363 = buffer.data(snf + 363);
    const auto *snf_365 = buffer.data(snf + 365);
    const auto *snf_366 = buffer.data(snf + 366);
    const auto *snf_367 = buffer.data(snf + 367);
    const auto *snf_368 = buffer.data(snf + 368);
    const auto *snf_369 = buffer.data(snf + 369);
    const auto *snf_375 = buffer.data(snf + 375);
    const auto *snf_376 = buffer.data(snf + 376);
    const auto *snf_377 = buffer.data(snf + 377);
    const auto *snf_378 = buffer.data(snf + 378);
    const auto *snf_379 = buffer.data(snf + 379);
    const auto *snf_380 = buffer.data(snf + 380);
    const auto *snf_383 = buffer.data(snf + 383);
    const auto *snf_385 = buffer.data(snf + 385);
    const auto *snf_386 = buffer.data(snf + 386);
    const auto *snf_387 = buffer.data(snf + 387);
    const auto *snf_388 = buffer.data(snf + 388);
    const auto *snf_389 = buffer.data(snf + 389);
    const auto *snf_390 = buffer.data(snf + 390);
    const auto *snf_393 = buffer.data(snf + 393);
    const auto *snf_395 = buffer.data(snf + 395);
    const auto *snf_396 = buffer.data(snf + 396);
    const auto *snf_397 = buffer.data(snf + 397);
    const auto *snf_398 = buffer.data(snf + 398);

    const auto *sng1_405 = buffer.data(sng1 + 405);
    const auto *sng1_408 = buffer.data(sng1 + 408);
    const auto *sng1_410 = buffer.data(sng1 + 410);
    const auto *sng1_419 = buffer.data(sng1 + 419);
    const auto *sng1_420 = buffer.data(sng1 + 420);
    const auto *sng1_423 = buffer.data(sng1 + 423);
    const auto *sng1_430 = buffer.data(sng1 + 430);

    const auto *sod0_191 = buffer.data(sod0 + 191);
    const auto *sod0_192 = buffer.data(sod0 + 192);
    const auto *sod0_195 = buffer.data(sod0 + 195);
    const auto *sod0_197 = buffer.data(sod0 + 197);
    const auto *sod0_198 = buffer.data(sod0 + 198);
    const auto *sod0_201 = buffer.data(sod0 + 201);
    const auto *sod0_203 = buffer.data(sod0 + 203);
    const auto *sod0_207 = buffer.data(sod0 + 207);
    const auto *sod0_209 = buffer.data(sod0 + 209);
    const auto *sod0_210 = buffer.data(sod0 + 210);
    const auto *sod0_213 = buffer.data(sod0 + 213);
    const auto *sod0_215 = buffer.data(sod0 + 215);
    const auto *sod0_216 = buffer.data(sod0 + 216);
    const auto *sod0_219 = buffer.data(sod0 + 219);
    const auto *sod0_221 = buffer.data(sod0 + 221);
    const auto *sod0_227 = buffer.data(sod0 + 227);
    const auto *sod0_228 = buffer.data(sod0 + 228);
    const auto *sod0_231 = buffer.data(sod0 + 231);
    const auto *sod0_233 = buffer.data(sod0 + 233);
    const auto *sod0_234 = buffer.data(sod0 + 234);
    const auto *sod0_237 = buffer.data(sod0 + 237);
    const auto *sod0_239 = buffer.data(sod0 + 239);

    const auto *sod1_191 = buffer.data(sod1 + 191);
    const auto *sod1_192 = buffer.data(sod1 + 192);
    const auto *sod1_195 = buffer.data(sod1 + 195);
    const auto *sod1_197 = buffer.data(sod1 + 197);
    const auto *sod1_198 = buffer.data(sod1 + 198);
    const auto *sod1_201 = buffer.data(sod1 + 201);
    const auto *sod1_203 = buffer.data(sod1 + 203);
    const auto *sod1_207 = buffer.data(sod1 + 207);
    const auto *sod1_209 = buffer.data(sod1 + 209);
    const auto *sod1_210 = buffer.data(sod1 + 210);
    const auto *sod1_213 = buffer.data(sod1 + 213);
    const auto *sod1_215 = buffer.data(sod1 + 215);
    const auto *sod1_216 = buffer.data(sod1 + 216);
    const auto *sod1_219 = buffer.data(sod1 + 219);
    const auto *sod1_221 = buffer.data(sod1 + 221);
    const auto *sod1_227 = buffer.data(sod1 + 227);
    const auto *sod1_228 = buffer.data(sod1 + 228);
    const auto *sod1_231 = buffer.data(sod1 + 231);
    const auto *sod1_233 = buffer.data(sod1 + 233);
    const auto *sod1_234 = buffer.data(sod1 + 234);
    const auto *sod1_237 = buffer.data(sod1 + 237);
    const auto *sod1_239 = buffer.data(sod1 + 239);

    const auto *sof_318 = buffer.data(sof + 318);
    const auto *sof_319 = buffer.data(sof + 319);
    const auto *sof_320 = buffer.data(sof + 320);
    const auto *sof_322 = buffer.data(sof + 322);
    const auto *sof_323 = buffer.data(sof + 323);
    const auto *sof_325 = buffer.data(sof + 325);
    const auto *sof_326 = buffer.data(sof + 326);
    const auto *sof_327 = buffer.data(sof + 327);
    const auto *sof_328 = buffer.data(sof + 328);
    const auto *sof_329 = buffer.data(sof + 329);
    const auto *sof_330 = buffer.data(sof + 330);
    const auto *sof_332 = buffer.data(sof + 332);
    const auto *sof_333 = buffer.data(sof + 333);
    const auto *sof_335 = buffer.data(sof + 335);
    const auto *sof_336 = buffer.data(sof + 336);
    const auto *sof_337 = buffer.data(sof + 337);
    const auto *sof_338 = buffer.data(sof + 338);
    const auto *sof_339 = buffer.data(sof + 339);
    const auto *sof_340 = buffer.data(sof + 340);
    const auto *sof_342 = buffer.data(sof + 342);
    const auto *sof_346 = buffer.data(sof + 346);
    const auto *sof_347 = buffer.data(sof + 347);
    const auto *sof_348 = buffer.data(sof + 348);
    const auto *sof_349 = buffer.data(sof + 349);
    const auto *sof_350 = buffer.data(sof + 350);
    const auto *sof_352 = buffer.data(sof + 352);
    const auto *sof_353 = buffer.data(sof + 353);
    const auto *sof_355 = buffer.data(sof + 355);
    const auto *sof_356 = buffer.data(sof + 356);
    const auto *sof_357 = buffer.data(sof + 357);
    const auto *sof_358 = buffer.data(sof + 358);
    const auto *sof_359 = buffer.data(sof + 359);
    const auto *sof_360 = buffer.data(sof + 360);
    const auto *sof_362 = buffer.data(sof + 362);
    const auto *sof_363 = buffer.data(sof + 363);
    const auto *sof_365 = buffer.data(sof + 365);
    const auto *sof_366 = buffer.data(sof + 366);
    const auto *sof_367 = buffer.data(sof + 367);
    const auto *sof_368 = buffer.data(sof + 368);
    const auto *sof_369 = buffer.data(sof + 369);
    const auto *sof_370 = buffer.data(sof + 370);
    const auto *sof_372 = buffer.data(sof + 372);
    const auto *sof_375 = buffer.data(sof + 375);
    const auto *sof_376 = buffer.data(sof + 376);
    const auto *sof_377 = buffer.data(sof + 377);
    const auto *sof_378 = buffer.data(sof + 378);
    const auto *sof_379 = buffer.data(sof + 379);
    const auto *sof_380 = buffer.data(sof + 380);
    const auto *sof_382 = buffer.data(sof + 382);
    const auto *sof_383 = buffer.data(sof + 383);
    const auto *sof_385 = buffer.data(sof + 385);
    const auto *sof_386 = buffer.data(sof + 386);
    const auto *sof_387 = buffer.data(sof + 387);
    const auto *sof_388 = buffer.data(sof + 388);
    const auto *sof_389 = buffer.data(sof + 389);
    const auto *sof_390 = buffer.data(sof + 390);
    const auto *sof_392 = buffer.data(sof + 392);
    const auto *sof_393 = buffer.data(sof + 393);
    const auto *sof_395 = buffer.data(sof + 395);
    const auto *sof_396 = buffer.data(sof + 396);
    const auto *sof_397 = buffer.data(sof + 397);
    const auto *sof_398 = buffer.data(sof + 398);

#pragma omp simd aligned(t_477, t_478, t_479, pc_y, pc_z, snf_239, snf_248, snf_249, sod0_191, \
                         sod1_191, sof_318, sof_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_14 * snf_248[k]
                   + f_4 * sod0_191[k]
                   - f_5 * sod1_191[k]
                   + f_3 * pc_y[k] * sof_318[k];

        t_478[k] = f_14 * snf_249[k]
                   + f_3 * pc_y[k] * sof_319[k];

        t_479[k] = f_12 * snf_239[k]
                   + f_1 * sod0_191[k]
                   - f_2 * sod1_191[k]
                   + f_3 * pc_z[k] * sof_319[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, pc_x, pc_y, pc_z, snf_240, snf_250, snf_320, \
                         sod0_192, sod1_192, sof_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_14 * snf_320[k]
                   + f_1 * sod0_192[k]
                   - f_2 * sod1_192[k]
                   + f_3 * pc_x[k] * sof_320[k];

        t_481[k] = f_12 * snf_250[k]
                   + f_3 * pc_y[k] * sof_320[k];

        t_482[k] = f_14 * snf_240[k]
                   + f_3 * pc_z[k] * sof_320[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pc_x, pc_y, snf_252, snf_323, snf_325, sod0_195, \
                         sod0_197, sod1_195, sod1_197, sof_322, sof_323, \
                         sof_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_14 * snf_323[k]
                   + f_4 * sod0_195[k]
                   - f_5 * sod1_195[k]
                   + f_3 * pc_x[k] * sof_323[k];

        t_484[k] = f_12 * snf_252[k]
                   + f_3 * pc_y[k] * sof_322[k];

        t_485[k] = f_14 * snf_325[k]
                   + f_4 * sod0_197[k]
                   - f_5 * sod1_197[k]
                   + f_3 * pc_x[k] * sof_325[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, pc_x, snf_326, snf_327, snf_328, snf_329, \
                         sof_326, sof_327, sof_328, sof_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_14 * snf_326[k]
                   + f_3 * pc_x[k] * sof_326[k];

        t_487[k] = f_14 * snf_327[k]
                   + f_3 * pc_x[k] * sof_327[k];

        t_488[k] = f_14 * snf_328[k]
                   + f_3 * pc_x[k] * sof_328[k];

        t_489[k] = f_14 * snf_329[k]
                   + f_3 * pc_x[k] * sof_329[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, pc_y, pc_z, snf_246, snf_256, snf_258, sod0_195, \
                         sod0_197, sod1_195, sod1_197, sof_326, \
                         sof_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = f_12 * snf_256[k]
                   + f_1 * sod0_195[k]
                   - f_2 * sod1_195[k]
                   + f_3 * pc_y[k] * sof_326[k];

        t_491[k] = f_14 * snf_246[k]
                   + f_3 * pc_z[k] * sof_326[k];

        t_492[k] = f_12 * snf_258[k]
                   + f_4 * sod0_197[k]
                   - f_5 * sod1_197[k]
                   + f_3 * pc_y[k] * sof_328[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, pc_x, pc_y, pc_z, snf_249, snf_259, snf_330, \
                         sod0_197, sod0_198, sod1_197, sod1_198, sof_329, \
                         sof_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_12 * snf_259[k]
                   + f_3 * pc_y[k] * sof_329[k];

        t_494[k] = f_14 * snf_249[k]
                   + f_1 * sod0_197[k]
                   - f_2 * sod1_197[k]
                   + f_3 * pc_z[k] * sof_329[k];

        t_495[k] = f_14 * snf_330[k]
                   + f_1 * sod0_198[k]
                   - f_2 * sod1_198[k]
                   + f_3 * pc_x[k] * sof_330[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pc_x, pc_y, pc_z, snf_250, snf_260, \
                         snf_262, snf_333, sod0_201, sod1_201, sof_330, sof_332, \
                         sof_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_8 * snf_260[k]
                   + f_3 * pc_y[k] * sof_330[k];

        t_497[k] = f_16 * snf_250[k]
                   + f_3 * pc_z[k] * sof_330[k];

        t_498[k] = f_14 * snf_333[k]
                   + f_4 * sod0_201[k]
                   - f_5 * sod1_201[k]
                   + f_3 * pc_x[k] * sof_333[k];

        t_499[k] = f_8 * snf_262[k]
                   + f_3 * pc_y[k] * sof_332[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pc_x, snf_335, snf_336, snf_337, snf_338, \
                         sod0_203, sod1_203, sof_335, sof_336, sof_337, \
                         sof_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_14 * snf_335[k]
                   + f_4 * sod0_203[k]
                   - f_5 * sod1_203[k]
                   + f_3 * pc_x[k] * sof_335[k];

        t_501[k] = f_14 * snf_336[k]
                   + f_3 * pc_x[k] * sof_336[k];

        t_502[k] = f_14 * snf_337[k]
                   + f_3 * pc_x[k] * sof_337[k];

        t_503[k] = f_14 * snf_338[k]
                   + f_3 * pc_x[k] * sof_338[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pc_x, pc_y, pc_z, snf_256, snf_266, snf_339, \
                         sod0_201, sod1_201, sof_336, sof_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_14 * snf_339[k]
                   + f_3 * pc_x[k] * sof_339[k];

        t_505[k] = f_8 * snf_266[k]
                   + f_1 * sod0_201[k]
                   - f_2 * sod1_201[k]
                   + f_3 * pc_y[k] * sof_336[k];

        t_506[k] = f_16 * snf_256[k]
                   + f_3 * pc_z[k] * sof_336[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, pb_y, pc_y, pc_z, sng0_405, snf_259, \
                         snf_268, snf_269, sng1_405, sod0_203, sod1_203, sof_338, \
                         sof_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_8 * snf_268[k]
                   + f_4 * sod0_203[k]
                   - f_5 * sod1_203[k]
                   + f_3 * pc_y[k] * sof_338[k];

        t_508[k] = f_8 * snf_269[k]
                   + f_3 * pc_y[k] * sof_339[k];

        t_509[k] = f_16 * snf_259[k]
                   + f_1 * sod0_203[k]
                   - f_2 * sod1_203[k]
                   + f_3 * pc_z[k] * sof_339[k];

        t_510[k] = pb_y[k] * sng0_405[k]
                   - f_6 * pc_y[k] * sng1_405[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pb_y, pc_y, pc_z, sng0_408, snf_260, \
                         snf_270, snf_271, snf_272, sng1_408, sof_340, \
                         sof_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_7 * snf_270[k]
                   + f_3 * pc_y[k] * sof_340[k];

        t_512[k] = f_15 * snf_260[k]
                   + f_3 * pc_z[k] * sof_340[k];

        t_513[k] = pb_y[k] * sng0_408[k]
                   + f_8 * snf_271[k]
                   - f_6 * pc_y[k] * sng1_408[k];

        t_514[k] = f_7 * snf_272[k]
                   + f_3 * pc_y[k] * sof_342[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pb_y, pc_x, pc_y, sng0_410, snf_346, \
                         snf_347, snf_348, sng1_410, sof_346, sof_347, \
                         sof_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pb_y[k] * sng0_410[k]
                   - f_6 * pc_y[k] * sng1_410[k];

        t_516[k] = f_14 * snf_346[k]
                   + f_3 * pc_x[k] * sof_346[k];

        t_517[k] = f_14 * snf_347[k]
                   + f_3 * pc_x[k] * sof_347[k];

        t_518[k] = f_14 * snf_348[k]
                   + f_3 * pc_x[k] * sof_348[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_x, pc_y, pc_z, snf_266, snf_276, snf_349, \
                         sod0_207, sod1_207, sof_346, sof_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_14 * snf_349[k]
                   + f_3 * pc_x[k] * sof_349[k];

        t_520[k] = f_7 * snf_276[k]
                   + f_1 * sod0_207[k]
                   - f_2 * sod1_207[k]
                   + f_3 * pc_y[k] * sof_346[k];

        t_521[k] = f_15 * snf_266[k]
                   + f_3 * pc_z[k] * sof_346[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, pb_y, pc_y, sng0_419, snf_278, snf_279, \
                         sng1_419, sod0_209, sod1_209, sof_348, \
                         sof_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_7 * snf_278[k]
                   + f_4 * sod0_209[k]
                   - f_5 * sod1_209[k]
                   + f_3 * pc_y[k] * sof_348[k];

        t_523[k] = f_7 * snf_279[k]
                   + f_3 * pc_y[k] * sof_349[k];

        t_524[k] = pb_y[k] * sng0_419[k]
                   - f_6 * pc_y[k] * sng1_419[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, pc_x, pc_y, pc_z, snf_270, snf_350, \
                         snf_353, sod0_210, sod0_213, sod1_210, sod1_213, sof_350, \
                         sof_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_14 * snf_350[k]
                   + f_1 * sod0_210[k]
                   - f_2 * sod1_210[k]
                   + f_3 * pc_x[k] * sof_350[k];

        t_526[k] = f_3 * pc_y[k] * sof_350[k];

        t_527[k] = f_13 * snf_270[k]
                   + f_3 * pc_z[k] * sof_350[k];

        t_528[k] = f_14 * snf_353[k]
                   + f_4 * sod0_213[k]
                   - f_5 * sod1_213[k]
                   + f_3 * pc_x[k] * sof_353[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pc_x, pc_y, snf_355, snf_356, snf_357, \
                         sod0_215, sod1_215, sof_352, sof_355, sof_356, \
                         sof_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_3 * pc_y[k] * sof_352[k];

        t_530[k] = f_14 * snf_355[k]
                   + f_4 * sod0_215[k]
                   - f_5 * sod1_215[k]
                   + f_3 * pc_x[k] * sof_355[k];

        t_531[k] = f_14 * snf_356[k]
                   + f_3 * pc_x[k] * sof_356[k];

        t_532[k] = f_14 * snf_357[k]
                   + f_3 * pc_x[k] * sof_357[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, pc_z, snf_276, snf_358, \
                         snf_359, sod0_213, sod1_213, sof_356, sof_358, \
                         sof_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_14 * snf_358[k]
                   + f_3 * pc_x[k] * sof_358[k];

        t_534[k] = f_14 * snf_359[k]
                   + f_3 * pc_x[k] * sof_359[k];

        t_535[k] = f_1 * sod0_213[k]
                   - f_2 * sod1_213[k]
                   + f_3 * pc_y[k] * sof_356[k];

        t_536[k] = f_13 * snf_276[k]
                   + f_3 * pc_z[k] * sof_356[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, snf_279, snf_360, \
                         sod0_215, sod0_216, sod1_215, sod1_216, sof_358, sof_359, \
                         sof_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * sod0_215[k]
                   - f_5 * sod1_215[k]
                   + f_3 * pc_y[k] * sof_358[k];

        t_538[k] = f_3 * pc_y[k] * sof_359[k];

        t_539[k] = f_13 * snf_279[k]
                   + f_1 * sod0_215[k]
                   - f_2 * sod1_215[k]
                   + f_3 * pc_z[k] * sof_359[k];

        t_540[k] = f_12 * snf_360[k]
                   + f_1 * sod0_216[k]
                   - f_2 * sod1_216[k]
                   + f_3 * pc_x[k] * sof_360[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, snf_280, snf_282, \
                         snf_363, sod0_219, sod1_219, sof_360, sof_362, \
                         sof_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_11 * snf_280[k]
                   + f_3 * pc_y[k] * sof_360[k];

        t_542[k] = f_3 * pc_z[k] * sof_360[k];

        t_543[k] = f_12 * snf_363[k]
                   + f_4 * sod0_219[k]
                   - f_5 * sod1_219[k]
                   + f_3 * pc_x[k] * sof_363[k];

        t_544[k] = f_11 * snf_282[k]
                   + f_3 * pc_y[k] * sof_362[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pc_x, snf_365, snf_366, snf_367, snf_368, \
                         sod0_221, sod1_221, sof_365, sof_366, sof_367, \
                         sof_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_12 * snf_365[k]
                   + f_4 * sod0_221[k]
                   - f_5 * sod1_221[k]
                   + f_3 * pc_x[k] * sof_365[k];

        t_546[k] = f_12 * snf_366[k]
                   + f_3 * pc_x[k] * sof_366[k];

        t_547[k] = f_12 * snf_367[k]
                   + f_3 * pc_x[k] * sof_367[k];

        t_548[k] = f_12 * snf_368[k]
                   + f_3 * pc_x[k] * sof_368[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pc_x, pc_y, pc_z, snf_286, snf_369, sod0_219, \
                         sod1_219, sof_366, sof_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_12 * snf_369[k]
                   + f_3 * pc_x[k] * sof_369[k];

        t_550[k] = f_11 * snf_286[k]
                   + f_1 * sod0_219[k]
                   - f_2 * sod1_219[k]
                   + f_3 * pc_y[k] * sof_366[k];

        t_551[k] = f_3 * pc_z[k] * sof_366[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pb_z, pc_y, pc_z, sng0_420, snf_288, \
                         snf_289, sng1_420, sod0_221, sod1_221, sof_368, \
                         sof_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_11 * snf_288[k]
                   + f_4 * sod0_221[k]
                   - f_5 * sod1_221[k]
                   + f_3 * pc_y[k] * sof_368[k];

        t_553[k] = f_11 * snf_289[k]
                   + f_3 * pc_y[k] * sof_369[k];

        t_554[k] = f_1 * sod0_221[k]
                   - f_2 * sod1_221[k]
                   + f_3 * pc_z[k] * sof_369[k];

        t_555[k] = pb_z[k] * sng0_420[k]
                   - f_6 * pc_z[k] * sng1_420[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, pb_z, pc_y, pc_z, sng0_423, snf_280, \
                         snf_290, snf_292, sng1_423, sof_370, sof_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_13 * snf_290[k]
                   + f_3 * pc_y[k] * sof_370[k];

        t_557[k] = f_7 * snf_280[k]
                   + f_3 * pc_z[k] * sof_370[k];

        t_558[k] = pb_z[k] * sng0_423[k]
                   - f_6 * pc_z[k] * sng1_423[k];

        t_559[k] = f_13 * snf_292[k]
                   + f_3 * pc_y[k] * sof_372[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, pc_x, snf_375, snf_376, snf_377, snf_378, \
                         sod0_227, sod1_227, sof_375, sof_376, sof_377, \
                         sof_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_12 * snf_375[k]
                   + f_4 * sod0_227[k]
                   - f_5 * sod1_227[k]
                   + f_3 * pc_x[k] * sof_375[k];

        t_561[k] = f_12 * snf_376[k]
                   + f_3 * pc_x[k] * sof_376[k];

        t_562[k] = f_12 * snf_377[k]
                   + f_3 * pc_x[k] * sof_377[k];

        t_563[k] = f_12 * snf_378[k]
                   + f_3 * pc_x[k] * sof_378[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pb_z, pc_x, pc_z, sng0_430, snf_286, snf_379, \
                         sng1_430, sof_376, sof_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_12 * snf_379[k]
                   + f_3 * pc_x[k] * sof_379[k];

        t_565[k] = pb_z[k] * sng0_430[k]
                   - f_6 * pc_z[k] * sng1_430[k];

        t_566[k] = f_7 * snf_286[k]
                   + f_3 * pc_z[k] * sof_376[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pc_y, pc_z, snf_289, snf_298, snf_299, sod0_227, \
                         sod1_227, sof_378, sof_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_13 * snf_298[k]
                   + f_4 * sod0_227[k]
                   - f_5 * sod1_227[k]
                   + f_3 * pc_y[k] * sof_378[k];

        t_568[k] = f_13 * snf_299[k]
                   + f_3 * pc_y[k] * sof_379[k];

        t_569[k] = f_7 * snf_289[k]
                   + f_1 * sod0_227[k]
                   - f_2 * sod1_227[k]
                   + f_3 * pc_z[k] * sof_379[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pc_x, pc_y, pc_z, snf_290, snf_300, snf_380, \
                         sod0_228, sod1_228, sof_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_12 * snf_380[k]
                   + f_1 * sod0_228[k]
                   - f_2 * sod1_228[k]
                   + f_3 * pc_x[k] * sof_380[k];

        t_571[k] = f_15 * snf_300[k]
                   + f_3 * pc_y[k] * sof_380[k];

        t_572[k] = f_8 * snf_290[k]
                   + f_3 * pc_z[k] * sof_380[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, pc_x, pc_y, snf_302, snf_383, snf_385, sod0_231, \
                         sod0_233, sod1_231, sod1_233, sof_382, sof_383, \
                         sof_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_12 * snf_383[k]
                   + f_4 * sod0_231[k]
                   - f_5 * sod1_231[k]
                   + f_3 * pc_x[k] * sof_383[k];

        t_574[k] = f_15 * snf_302[k]
                   + f_3 * pc_y[k] * sof_382[k];

        t_575[k] = f_12 * snf_385[k]
                   + f_4 * sod0_233[k]
                   - f_5 * sod1_233[k]
                   + f_3 * pc_x[k] * sof_385[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pc_x, snf_386, snf_387, snf_388, snf_389, \
                         sof_386, sof_387, sof_388, sof_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_12 * snf_386[k]
                   + f_3 * pc_x[k] * sof_386[k];

        t_577[k] = f_12 * snf_387[k]
                   + f_3 * pc_x[k] * sof_387[k];

        t_578[k] = f_12 * snf_388[k]
                   + f_3 * pc_x[k] * sof_388[k];

        t_579[k] = f_12 * snf_389[k]
                   + f_3 * pc_x[k] * sof_389[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, pc_y, pc_z, snf_296, snf_306, snf_308, sod0_231, \
                         sod0_233, sod1_231, sod1_233, sof_386, \
                         sof_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_15 * snf_306[k]
                   + f_1 * sod0_231[k]
                   - f_2 * sod1_231[k]
                   + f_3 * pc_y[k] * sof_386[k];

        t_581[k] = f_8 * snf_296[k]
                   + f_3 * pc_z[k] * sof_386[k];

        t_582[k] = f_15 * snf_308[k]
                   + f_4 * sod0_233[k]
                   - f_5 * sod1_233[k]
                   + f_3 * pc_y[k] * sof_388[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, pc_x, pc_y, pc_z, snf_299, snf_309, snf_390, \
                         sod0_233, sod0_234, sod1_233, sod1_234, sof_389, \
                         sof_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_15 * snf_309[k]
                   + f_3 * pc_y[k] * sof_389[k];

        t_584[k] = f_8 * snf_299[k]
                   + f_1 * sod0_233[k]
                   - f_2 * sod1_233[k]
                   + f_3 * pc_z[k] * sof_389[k];

        t_585[k] = f_12 * snf_390[k]
                   + f_1 * sod0_234[k]
                   - f_2 * sod1_234[k]
                   + f_3 * pc_x[k] * sof_390[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pc_x, pc_y, pc_z, snf_300, snf_310, \
                         snf_312, snf_393, sod0_237, sod1_237, sof_390, sof_392, \
                         sof_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_16 * snf_310[k]
                   + f_3 * pc_y[k] * sof_390[k];

        t_587[k] = f_12 * snf_300[k]
                   + f_3 * pc_z[k] * sof_390[k];

        t_588[k] = f_12 * snf_393[k]
                   + f_4 * sod0_237[k]
                   - f_5 * sod1_237[k]
                   + f_3 * pc_x[k] * sof_393[k];

        t_589[k] = f_16 * snf_312[k]
                   + f_3 * pc_y[k] * sof_392[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pc_x, snf_395, snf_396, snf_397, snf_398, \
                         sod0_239, sod1_239, sof_395, sof_396, sof_397, \
                         sof_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_12 * snf_395[k]
                   + f_4 * sod0_239[k]
                   - f_5 * sod1_239[k]
                   + f_3 * pc_x[k] * sof_395[k];

        t_591[k] = f_12 * snf_396[k]
                   + f_3 * pc_x[k] * sof_396[k];

        t_592[k] = f_12 * snf_397[k]
                   + f_3 * pc_x[k] * sof_397[k];

        t_593[k] = f_12 * snf_398[k]
                   + f_3 * pc_x[k] * sof_398[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
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
    const auto f_10 = 4.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_525 = buffer.data(sng0 + 525);
    const auto *sng0_528 = buffer.data(sng0 + 528);
    const auto *sng0_530 = buffer.data(sng0 + 530);
    const auto *sng0_539 = buffer.data(sng0 + 539);
    const auto *sng0_540 = buffer.data(sng0 + 540);
    const auto *sng0_543 = buffer.data(sng0 + 543);
    const auto *sng0_550 = buffer.data(sng0 + 550);

    const auto *snf_306 = buffer.data(snf + 306);
    const auto *snf_309 = buffer.data(snf + 309);
    const auto *snf_310 = buffer.data(snf + 310);
    const auto *snf_316 = buffer.data(snf + 316);
    const auto *snf_318 = buffer.data(snf + 318);
    const auto *snf_319 = buffer.data(snf + 319);
    const auto *snf_320 = buffer.data(snf + 320);
    const auto *snf_322 = buffer.data(snf + 322);
    const auto *snf_326 = buffer.data(snf + 326);
    const auto *snf_328 = buffer.data(snf + 328);
    const auto *snf_329 = buffer.data(snf + 329);
    const auto *snf_330 = buffer.data(snf + 330);
    const auto *snf_332 = buffer.data(snf + 332);
    const auto *snf_336 = buffer.data(snf + 336);
    const auto *snf_338 = buffer.data(snf + 338);
    const auto *snf_339 = buffer.data(snf + 339);
    const auto *snf_340 = buffer.data(snf + 340);
    const auto *snf_342 = buffer.data(snf + 342);
    const auto *snf_346 = buffer.data(snf + 346);
    const auto *snf_348 = buffer.data(snf + 348);
    const auto *snf_349 = buffer.data(snf + 349);
    const auto *snf_350 = buffer.data(snf + 350);
    const auto *snf_351 = buffer.data(snf + 351);
    const auto *snf_352 = buffer.data(snf + 352);
    const auto *snf_356 = buffer.data(snf + 356);
    const auto *snf_358 = buffer.data(snf + 358);
    const auto *snf_359 = buffer.data(snf + 359);
    const auto *snf_360 = buffer.data(snf + 360);
    const auto *snf_362 = buffer.data(snf + 362);
    const auto *snf_366 = buffer.data(snf + 366);
    const auto *snf_368 = buffer.data(snf + 368);
    const auto *snf_369 = buffer.data(snf + 369);
    const auto *snf_370 = buffer.data(snf + 370);
    const auto *snf_372 = buffer.data(snf + 372);
    const auto *snf_378 = buffer.data(snf + 378);
    const auto *snf_379 = buffer.data(snf + 379);
    const auto *snf_380 = buffer.data(snf + 380);
    const auto *snf_399 = buffer.data(snf + 399);
    const auto *snf_400 = buffer.data(snf + 400);
    const auto *snf_403 = buffer.data(snf + 403);
    const auto *snf_405 = buffer.data(snf + 405);
    const auto *snf_406 = buffer.data(snf + 406);
    const auto *snf_407 = buffer.data(snf + 407);
    const auto *snf_408 = buffer.data(snf + 408);
    const auto *snf_409 = buffer.data(snf + 409);
    const auto *snf_410 = buffer.data(snf + 410);
    const auto *snf_413 = buffer.data(snf + 413);
    const auto *snf_415 = buffer.data(snf + 415);
    const auto *snf_416 = buffer.data(snf + 416);
    const auto *snf_417 = buffer.data(snf + 417);
    const auto *snf_418 = buffer.data(snf + 418);
    const auto *snf_419 = buffer.data(snf + 419);
    const auto *snf_420 = buffer.data(snf + 420);
    const auto *snf_423 = buffer.data(snf + 423);
    const auto *snf_425 = buffer.data(snf + 425);
    const auto *snf_426 = buffer.data(snf + 426);
    const auto *snf_427 = buffer.data(snf + 427);
    const auto *snf_428 = buffer.data(snf + 428);
    const auto *snf_429 = buffer.data(snf + 429);
    const auto *snf_436 = buffer.data(snf + 436);
    const auto *snf_437 = buffer.data(snf + 437);
    const auto *snf_438 = buffer.data(snf + 438);
    const auto *snf_439 = buffer.data(snf + 439);
    const auto *snf_440 = buffer.data(snf + 440);
    const auto *snf_443 = buffer.data(snf + 443);
    const auto *snf_445 = buffer.data(snf + 445);
    const auto *snf_446 = buffer.data(snf + 446);
    const auto *snf_447 = buffer.data(snf + 447);
    const auto *snf_448 = buffer.data(snf + 448);
    const auto *snf_449 = buffer.data(snf + 449);
    const auto *snf_450 = buffer.data(snf + 450);
    const auto *snf_453 = buffer.data(snf + 453);
    const auto *snf_455 = buffer.data(snf + 455);
    const auto *snf_456 = buffer.data(snf + 456);
    const auto *snf_457 = buffer.data(snf + 457);
    const auto *snf_458 = buffer.data(snf + 458);
    const auto *snf_459 = buffer.data(snf + 459);
    const auto *snf_465 = buffer.data(snf + 465);
    const auto *snf_466 = buffer.data(snf + 466);
    const auto *snf_467 = buffer.data(snf + 467);
    const auto *snf_468 = buffer.data(snf + 468);
    const auto *snf_469 = buffer.data(snf + 469);
    const auto *snf_470 = buffer.data(snf + 470);

    const auto *sng1_525 = buffer.data(sng1 + 525);
    const auto *sng1_528 = buffer.data(sng1 + 528);
    const auto *sng1_530 = buffer.data(sng1 + 530);
    const auto *sng1_539 = buffer.data(sng1 + 539);
    const auto *sng1_540 = buffer.data(sng1 + 540);
    const auto *sng1_543 = buffer.data(sng1 + 543);
    const auto *sng1_550 = buffer.data(sng1 + 550);

    const auto *sod0_237 = buffer.data(sod0 + 237);
    const auto *sod0_239 = buffer.data(sod0 + 239);
    const auto *sod0_240 = buffer.data(sod0 + 240);
    const auto *sod0_243 = buffer.data(sod0 + 243);
    const auto *sod0_245 = buffer.data(sod0 + 245);
    const auto *sod0_246 = buffer.data(sod0 + 246);
    const auto *sod0_249 = buffer.data(sod0 + 249);
    const auto *sod0_251 = buffer.data(sod0 + 251);
    const auto *sod0_252 = buffer.data(sod0 + 252);
    const auto *sod0_255 = buffer.data(sod0 + 255);
    const auto *sod0_257 = buffer.data(sod0 + 257);
    const auto *sod0_261 = buffer.data(sod0 + 261);
    const auto *sod0_263 = buffer.data(sod0 + 263);
    const auto *sod0_264 = buffer.data(sod0 + 264);
    const auto *sod0_267 = buffer.data(sod0 + 267);
    const auto *sod0_269 = buffer.data(sod0 + 269);
    const auto *sod0_270 = buffer.data(sod0 + 270);
    const auto *sod0_273 = buffer.data(sod0 + 273);
    const auto *sod0_275 = buffer.data(sod0 + 275);
    const auto *sod0_281 = buffer.data(sod0 + 281);
    const auto *sod0_282 = buffer.data(sod0 + 282);

    const auto *sod1_237 = buffer.data(sod1 + 237);
    const auto *sod1_239 = buffer.data(sod1 + 239);
    const auto *sod1_240 = buffer.data(sod1 + 240);
    const auto *sod1_243 = buffer.data(sod1 + 243);
    const auto *sod1_245 = buffer.data(sod1 + 245);
    const auto *sod1_246 = buffer.data(sod1 + 246);
    const auto *sod1_249 = buffer.data(sod1 + 249);
    const auto *sod1_251 = buffer.data(sod1 + 251);
    const auto *sod1_252 = buffer.data(sod1 + 252);
    const auto *sod1_255 = buffer.data(sod1 + 255);
    const auto *sod1_257 = buffer.data(sod1 + 257);
    const auto *sod1_261 = buffer.data(sod1 + 261);
    const auto *sod1_263 = buffer.data(sod1 + 263);
    const auto *sod1_264 = buffer.data(sod1 + 264);
    const auto *sod1_267 = buffer.data(sod1 + 267);
    const auto *sod1_269 = buffer.data(sod1 + 269);
    const auto *sod1_270 = buffer.data(sod1 + 270);
    const auto *sod1_273 = buffer.data(sod1 + 273);
    const auto *sod1_275 = buffer.data(sod1 + 275);
    const auto *sod1_281 = buffer.data(sod1 + 281);
    const auto *sod1_282 = buffer.data(sod1 + 282);

    const auto *sof_396 = buffer.data(sof + 396);
    const auto *sof_398 = buffer.data(sof + 398);
    const auto *sof_399 = buffer.data(sof + 399);
    const auto *sof_400 = buffer.data(sof + 400);
    const auto *sof_402 = buffer.data(sof + 402);
    const auto *sof_403 = buffer.data(sof + 403);
    const auto *sof_405 = buffer.data(sof + 405);
    const auto *sof_406 = buffer.data(sof + 406);
    const auto *sof_407 = buffer.data(sof + 407);
    const auto *sof_408 = buffer.data(sof + 408);
    const auto *sof_409 = buffer.data(sof + 409);
    const auto *sof_410 = buffer.data(sof + 410);
    const auto *sof_412 = buffer.data(sof + 412);
    const auto *sof_413 = buffer.data(sof + 413);
    const auto *sof_415 = buffer.data(sof + 415);
    const auto *sof_416 = buffer.data(sof + 416);
    const auto *sof_417 = buffer.data(sof + 417);
    const auto *sof_418 = buffer.data(sof + 418);
    const auto *sof_419 = buffer.data(sof + 419);
    const auto *sof_420 = buffer.data(sof + 420);
    const auto *sof_422 = buffer.data(sof + 422);
    const auto *sof_423 = buffer.data(sof + 423);
    const auto *sof_425 = buffer.data(sof + 425);
    const auto *sof_426 = buffer.data(sof + 426);
    const auto *sof_427 = buffer.data(sof + 427);
    const auto *sof_428 = buffer.data(sof + 428);
    const auto *sof_429 = buffer.data(sof + 429);
    const auto *sof_430 = buffer.data(sof + 430);
    const auto *sof_432 = buffer.data(sof + 432);
    const auto *sof_436 = buffer.data(sof + 436);
    const auto *sof_437 = buffer.data(sof + 437);
    const auto *sof_438 = buffer.data(sof + 438);
    const auto *sof_439 = buffer.data(sof + 439);
    const auto *sof_440 = buffer.data(sof + 440);
    const auto *sof_442 = buffer.data(sof + 442);
    const auto *sof_443 = buffer.data(sof + 443);
    const auto *sof_445 = buffer.data(sof + 445);
    const auto *sof_446 = buffer.data(sof + 446);
    const auto *sof_447 = buffer.data(sof + 447);
    const auto *sof_448 = buffer.data(sof + 448);
    const auto *sof_449 = buffer.data(sof + 449);
    const auto *sof_450 = buffer.data(sof + 450);
    const auto *sof_452 = buffer.data(sof + 452);
    const auto *sof_453 = buffer.data(sof + 453);
    const auto *sof_455 = buffer.data(sof + 455);
    const auto *sof_456 = buffer.data(sof + 456);
    const auto *sof_457 = buffer.data(sof + 457);
    const auto *sof_458 = buffer.data(sof + 458);
    const auto *sof_459 = buffer.data(sof + 459);
    const auto *sof_460 = buffer.data(sof + 460);
    const auto *sof_462 = buffer.data(sof + 462);
    const auto *sof_465 = buffer.data(sof + 465);
    const auto *sof_466 = buffer.data(sof + 466);
    const auto *sof_467 = buffer.data(sof + 467);
    const auto *sof_468 = buffer.data(sof + 468);
    const auto *sof_469 = buffer.data(sof + 469);
    const auto *sof_470 = buffer.data(sof + 470);

#pragma omp simd aligned(t_594, t_595, t_596, pc_x, pc_y, pc_z, snf_306, snf_316, snf_399, \
                         sod0_237, sod1_237, sof_396, sof_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_12 * snf_399[k]
                   + f_3 * pc_x[k] * sof_399[k];

        t_595[k] = f_16 * snf_316[k]
                   + f_1 * sod0_237[k]
                   - f_2 * sod1_237[k]
                   + f_3 * pc_y[k] * sof_396[k];

        t_596[k] = f_12 * snf_306[k]
                   + f_3 * pc_z[k] * sof_396[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, pc_y, pc_z, snf_309, snf_318, snf_319, sod0_239, \
                         sod1_239, sof_398, sof_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_16 * snf_318[k]
                   + f_4 * sod0_239[k]
                   - f_5 * sod1_239[k]
                   + f_3 * pc_y[k] * sof_398[k];

        t_598[k] = f_16 * snf_319[k]
                   + f_3 * pc_y[k] * sof_399[k];

        t_599[k] = f_12 * snf_309[k]
                   + f_1 * sod0_239[k]
                   - f_2 * sod1_239[k]
                   + f_3 * pc_z[k] * sof_399[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pc_x, pc_y, pc_z, snf_310, snf_320, snf_400, \
                         sod0_240, sod1_240, sof_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_12 * snf_400[k]
                   + f_1 * sod0_240[k]
                   - f_2 * sod1_240[k]
                   + f_3 * pc_x[k] * sof_400[k];

        t_601[k] = f_14 * snf_320[k]
                   + f_3 * pc_y[k] * sof_400[k];

        t_602[k] = f_14 * snf_310[k]
                   + f_3 * pc_z[k] * sof_400[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, pc_x, pc_y, snf_322, snf_403, snf_405, sod0_243, \
                         sod0_245, sod1_243, sod1_245, sof_402, sof_403, \
                         sof_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = f_12 * snf_403[k]
                   + f_4 * sod0_243[k]
                   - f_5 * sod1_243[k]
                   + f_3 * pc_x[k] * sof_403[k];

        t_604[k] = f_14 * snf_322[k]
                   + f_3 * pc_y[k] * sof_402[k];

        t_605[k] = f_12 * snf_405[k]
                   + f_4 * sod0_245[k]
                   - f_5 * sod1_245[k]
                   + f_3 * pc_x[k] * sof_405[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, t_609, pc_x, snf_406, snf_407, snf_408, snf_409, \
                         sof_406, sof_407, sof_408, sof_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = f_12 * snf_406[k]
                   + f_3 * pc_x[k] * sof_406[k];

        t_607[k] = f_12 * snf_407[k]
                   + f_3 * pc_x[k] * sof_407[k];

        t_608[k] = f_12 * snf_408[k]
                   + f_3 * pc_x[k] * sof_408[k];

        t_609[k] = f_12 * snf_409[k]
                   + f_3 * pc_x[k] * sof_409[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, pc_y, pc_z, snf_316, snf_326, snf_328, sod0_243, \
                         sod0_245, sod1_243, sod1_245, sof_406, \
                         sof_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_14 * snf_326[k]
                   + f_1 * sod0_243[k]
                   - f_2 * sod1_243[k]
                   + f_3 * pc_y[k] * sof_406[k];

        t_611[k] = f_14 * snf_316[k]
                   + f_3 * pc_z[k] * sof_406[k];

        t_612[k] = f_14 * snf_328[k]
                   + f_4 * sod0_245[k]
                   - f_5 * sod1_245[k]
                   + f_3 * pc_y[k] * sof_408[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, pc_x, pc_y, pc_z, snf_319, snf_329, snf_410, \
                         sod0_245, sod0_246, sod1_245, sod1_246, sof_409, \
                         sof_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_14 * snf_329[k]
                   + f_3 * pc_y[k] * sof_409[k];

        t_614[k] = f_14 * snf_319[k]
                   + f_1 * sod0_245[k]
                   - f_2 * sod1_245[k]
                   + f_3 * pc_z[k] * sof_409[k];

        t_615[k] = f_12 * snf_410[k]
                   + f_1 * sod0_246[k]
                   - f_2 * sod1_246[k]
                   + f_3 * pc_x[k] * sof_410[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pc_x, pc_y, pc_z, snf_320, snf_330, \
                         snf_332, snf_413, sod0_249, sod1_249, sof_410, sof_412, \
                         sof_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = f_12 * snf_330[k]
                   + f_3 * pc_y[k] * sof_410[k];

        t_617[k] = f_16 * snf_320[k]
                   + f_3 * pc_z[k] * sof_410[k];

        t_618[k] = f_12 * snf_413[k]
                   + f_4 * sod0_249[k]
                   - f_5 * sod1_249[k]
                   + f_3 * pc_x[k] * sof_413[k];

        t_619[k] = f_12 * snf_332[k]
                   + f_3 * pc_y[k] * sof_412[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, pc_x, snf_415, snf_416, snf_417, snf_418, \
                         sod0_251, sod1_251, sof_415, sof_416, sof_417, \
                         sof_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_12 * snf_415[k]
                   + f_4 * sod0_251[k]
                   - f_5 * sod1_251[k]
                   + f_3 * pc_x[k] * sof_415[k];

        t_621[k] = f_12 * snf_416[k]
                   + f_3 * pc_x[k] * sof_416[k];

        t_622[k] = f_12 * snf_417[k]
                   + f_3 * pc_x[k] * sof_417[k];

        t_623[k] = f_12 * snf_418[k]
                   + f_3 * pc_x[k] * sof_418[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pc_x, pc_y, pc_z, snf_326, snf_336, snf_419, \
                         sod0_249, sod1_249, sof_416, sof_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_12 * snf_419[k]
                   + f_3 * pc_x[k] * sof_419[k];

        t_625[k] = f_12 * snf_336[k]
                   + f_1 * sod0_249[k]
                   - f_2 * sod1_249[k]
                   + f_3 * pc_y[k] * sof_416[k];

        t_626[k] = f_16 * snf_326[k]
                   + f_3 * pc_z[k] * sof_416[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pc_y, pc_z, snf_329, snf_338, snf_339, sod0_251, \
                         sod1_251, sof_418, sof_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_12 * snf_338[k]
                   + f_4 * sod0_251[k]
                   - f_5 * sod1_251[k]
                   + f_3 * pc_y[k] * sof_418[k];

        t_628[k] = f_12 * snf_339[k]
                   + f_3 * pc_y[k] * sof_419[k];

        t_629[k] = f_16 * snf_329[k]
                   + f_1 * sod0_251[k]
                   - f_2 * sod1_251[k]
                   + f_3 * pc_z[k] * sof_419[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, pc_x, pc_y, pc_z, snf_330, snf_340, snf_420, \
                         sod0_252, sod1_252, sof_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_12 * snf_420[k]
                   + f_1 * sod0_252[k]
                   - f_2 * sod1_252[k]
                   + f_3 * pc_x[k] * sof_420[k];

        t_631[k] = f_8 * snf_340[k]
                   + f_3 * pc_y[k] * sof_420[k];

        t_632[k] = f_15 * snf_330[k]
                   + f_3 * pc_z[k] * sof_420[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, snf_342, snf_423, snf_425, sod0_255, \
                         sod0_257, sod1_255, sod1_257, sof_422, sof_423, \
                         sof_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_12 * snf_423[k]
                   + f_4 * sod0_255[k]
                   - f_5 * sod1_255[k]
                   + f_3 * pc_x[k] * sof_423[k];

        t_634[k] = f_8 * snf_342[k]
                   + f_3 * pc_y[k] * sof_422[k];

        t_635[k] = f_12 * snf_425[k]
                   + f_4 * sod0_257[k]
                   - f_5 * sod1_257[k]
                   + f_3 * pc_x[k] * sof_425[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pc_x, snf_426, snf_427, snf_428, snf_429, \
                         sof_426, sof_427, sof_428, sof_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_12 * snf_426[k]
                   + f_3 * pc_x[k] * sof_426[k];

        t_637[k] = f_12 * snf_427[k]
                   + f_3 * pc_x[k] * sof_427[k];

        t_638[k] = f_12 * snf_428[k]
                   + f_3 * pc_x[k] * sof_428[k];

        t_639[k] = f_12 * snf_429[k]
                   + f_3 * pc_x[k] * sof_429[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pc_y, pc_z, snf_336, snf_346, snf_348, sod0_255, \
                         sod0_257, sod1_255, sod1_257, sof_426, \
                         sof_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_8 * snf_346[k]
                   + f_1 * sod0_255[k]
                   - f_2 * sod1_255[k]
                   + f_3 * pc_y[k] * sof_426[k];

        t_641[k] = f_15 * snf_336[k]
                   + f_3 * pc_z[k] * sof_426[k];

        t_642[k] = f_8 * snf_348[k]
                   + f_4 * sod0_257[k]
                   - f_5 * sod1_257[k]
                   + f_3 * pc_y[k] * sof_428[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pb_y, pc_y, pc_z, sng0_525, snf_339, \
                         snf_349, snf_350, sng1_525, sod0_257, sod1_257, sof_429, \
                         sof_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_8 * snf_349[k]
                   + f_3 * pc_y[k] * sof_429[k];

        t_644[k] = f_15 * snf_339[k]
                   + f_1 * sod0_257[k]
                   - f_2 * sod1_257[k]
                   + f_3 * pc_z[k] * sof_429[k];

        t_645[k] = pb_y[k] * sng0_525[k]
                   - f_6 * pc_y[k] * sng1_525[k];

        t_646[k] = f_7 * snf_350[k]
                   + f_3 * pc_y[k] * sof_430[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, pb_y, pc_y, pc_z, sng0_528, sng0_530, \
                         snf_340, snf_351, snf_352, sng1_528, sng1_530, sof_430, \
                         sof_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_13 * snf_340[k]
                   + f_3 * pc_z[k] * sof_430[k];

        t_648[k] = pb_y[k] * sng0_528[k]
                   + f_8 * snf_351[k]
                   - f_6 * pc_y[k] * sng1_528[k];

        t_649[k] = f_7 * snf_352[k]
                   + f_3 * pc_y[k] * sof_432[k];

        t_650[k] = pb_y[k] * sng0_530[k]
                   - f_6 * pc_y[k] * sng1_530[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, pc_x, snf_436, snf_437, snf_438, snf_439, \
                         sof_436, sof_437, sof_438, sof_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_12 * snf_436[k]
                   + f_3 * pc_x[k] * sof_436[k];

        t_652[k] = f_12 * snf_437[k]
                   + f_3 * pc_x[k] * sof_437[k];

        t_653[k] = f_12 * snf_438[k]
                   + f_3 * pc_x[k] * sof_438[k];

        t_654[k] = f_12 * snf_439[k]
                   + f_3 * pc_x[k] * sof_439[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pc_y, pc_z, snf_346, snf_356, snf_358, sod0_261, \
                         sod0_263, sod1_261, sod1_263, sof_436, \
                         sof_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_7 * snf_356[k]
                   + f_1 * sod0_261[k]
                   - f_2 * sod1_261[k]
                   + f_3 * pc_y[k] * sof_436[k];

        t_656[k] = f_13 * snf_346[k]
                   + f_3 * pc_z[k] * sof_436[k];

        t_657[k] = f_7 * snf_358[k]
                   + f_4 * sod0_263[k]
                   - f_5 * sod1_263[k]
                   + f_3 * pc_y[k] * sof_438[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pb_y, pc_x, pc_y, sng0_539, snf_359, \
                         snf_440, sng1_539, sod0_264, sod1_264, sof_439, \
                         sof_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_7 * snf_359[k]
                   + f_3 * pc_y[k] * sof_439[k];

        t_659[k] = pb_y[k] * sng0_539[k]
                   - f_6 * pc_y[k] * sng1_539[k];

        t_660[k] = f_12 * snf_440[k]
                   + f_1 * sod0_264[k]
                   - f_2 * sod1_264[k]
                   + f_3 * pc_x[k] * sof_440[k];

        t_661[k] = f_3 * pc_y[k] * sof_440[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, pc_x, pc_y, pc_z, snf_350, snf_443, sod0_267, \
                         sod1_267, sof_440, sof_442, sof_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_11 * snf_350[k]
                   + f_3 * pc_z[k] * sof_440[k];

        t_663[k] = f_12 * snf_443[k]
                   + f_4 * sod0_267[k]
                   - f_5 * sod1_267[k]
                   + f_3 * pc_x[k] * sof_443[k];

        t_664[k] = f_3 * pc_y[k] * sof_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pc_x, snf_445, snf_446, snf_447, snf_448, \
                         sod0_269, sod1_269, sof_445, sof_446, sof_447, \
                         sof_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_12 * snf_445[k]
                   + f_4 * sod0_269[k]
                   - f_5 * sod1_269[k]
                   + f_3 * pc_x[k] * sof_445[k];

        t_666[k] = f_12 * snf_446[k]
                   + f_3 * pc_x[k] * sof_446[k];

        t_667[k] = f_12 * snf_447[k]
                   + f_3 * pc_x[k] * sof_447[k];

        t_668[k] = f_12 * snf_448[k]
                   + f_3 * pc_x[k] * sof_448[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, t_672, pc_x, pc_y, pc_z, snf_356, snf_449, \
                         sod0_267, sod0_269, sod1_267, sod1_269, sof_446, sof_448, \
                         sof_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_12 * snf_449[k]
                   + f_3 * pc_x[k] * sof_449[k];

        t_670[k] = f_1 * sod0_267[k]
                   - f_2 * sod1_267[k]
                   + f_3 * pc_y[k] * sof_446[k];

        t_671[k] = f_11 * snf_356[k]
                   + f_3 * pc_z[k] * sof_446[k];

        t_672[k] = f_4 * sod0_269[k]
                   - f_5 * sod1_269[k]
                   + f_3 * pc_y[k] * sof_448[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, snf_359, snf_360, \
                         snf_450, sod0_269, sod0_270, sod1_269, sod1_270, sof_449, \
                         sof_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_3 * pc_y[k] * sof_449[k];

        t_674[k] = f_11 * snf_359[k]
                   + f_1 * sod0_269[k]
                   - f_2 * sod1_269[k]
                   + f_3 * pc_z[k] * sof_449[k];

        t_675[k] = f_8 * snf_450[k]
                   + f_1 * sod0_270[k]
                   - f_2 * sod1_270[k]
                   + f_3 * pc_x[k] * sof_450[k];

        t_676[k] = f_10 * snf_360[k]
                   + f_3 * pc_y[k] * sof_450[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_y, pc_z, snf_362, snf_453, sod0_273, \
                         sod1_273, sof_450, sof_452, sof_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_3 * pc_z[k] * sof_450[k];

        t_678[k] = f_8 * snf_453[k]
                   + f_4 * sod0_273[k]
                   - f_5 * sod1_273[k]
                   + f_3 * pc_x[k] * sof_453[k];

        t_679[k] = f_10 * snf_362[k]
                   + f_3 * pc_y[k] * sof_452[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, snf_455, snf_456, snf_457, snf_458, \
                         sod0_275, sod1_275, sof_455, sof_456, sof_457, \
                         sof_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_8 * snf_455[k]
                   + f_4 * sod0_275[k]
                   - f_5 * sod1_275[k]
                   + f_3 * pc_x[k] * sof_455[k];

        t_681[k] = f_8 * snf_456[k]
                   + f_3 * pc_x[k] * sof_456[k];

        t_682[k] = f_8 * snf_457[k]
                   + f_3 * pc_x[k] * sof_457[k];

        t_683[k] = f_8 * snf_458[k]
                   + f_3 * pc_x[k] * sof_458[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, pc_x, pc_y, pc_z, snf_366, snf_459, sod0_273, \
                         sod1_273, sof_456, sof_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_8 * snf_459[k]
                   + f_3 * pc_x[k] * sof_459[k];

        t_685[k] = f_10 * snf_366[k]
                   + f_1 * sod0_273[k]
                   - f_2 * sod1_273[k]
                   + f_3 * pc_y[k] * sof_456[k];

        t_686[k] = f_3 * pc_z[k] * sof_456[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, pb_z, pc_y, pc_z, sng0_540, snf_368, \
                         snf_369, sng1_540, sod0_275, sod1_275, sof_458, \
                         sof_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_10 * snf_368[k]
                   + f_4 * sod0_275[k]
                   - f_5 * sod1_275[k]
                   + f_3 * pc_y[k] * sof_458[k];

        t_688[k] = f_10 * snf_369[k]
                   + f_3 * pc_y[k] * sof_459[k];

        t_689[k] = f_1 * sod0_275[k]
                   - f_2 * sod1_275[k]
                   + f_3 * pc_z[k] * sof_459[k];

        t_690[k] = pb_z[k] * sng0_540[k]
                   - f_6 * pc_z[k] * sng1_540[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, pb_z, pc_y, pc_z, sng0_543, snf_360, \
                         snf_370, snf_372, sng1_543, sof_460, sof_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_11 * snf_370[k]
                   + f_3 * pc_y[k] * sof_460[k];

        t_692[k] = f_7 * snf_360[k]
                   + f_3 * pc_z[k] * sof_460[k];

        t_693[k] = pb_z[k] * sng0_543[k]
                   - f_6 * pc_z[k] * sng1_543[k];

        t_694[k] = f_11 * snf_372[k]
                   + f_3 * pc_y[k] * sof_462[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, pc_x, snf_465, snf_466, snf_467, snf_468, \
                         sod0_281, sod1_281, sof_465, sof_466, sof_467, \
                         sof_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_8 * snf_465[k]
                   + f_4 * sod0_281[k]
                   - f_5 * sod1_281[k]
                   + f_3 * pc_x[k] * sof_465[k];

        t_696[k] = f_8 * snf_466[k]
                   + f_3 * pc_x[k] * sof_466[k];

        t_697[k] = f_8 * snf_467[k]
                   + f_3 * pc_x[k] * sof_467[k];

        t_698[k] = f_8 * snf_468[k]
                   + f_3 * pc_x[k] * sof_468[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, pb_z, pc_x, pc_z, sng0_550, snf_366, snf_469, \
                         sng1_550, sof_466, sof_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_8 * snf_469[k]
                   + f_3 * pc_x[k] * sof_469[k];

        t_700[k] = pb_z[k] * sng0_550[k]
                   - f_6 * pc_z[k] * sng1_550[k];

        t_701[k] = f_7 * snf_366[k]
                   + f_3 * pc_z[k] * sof_466[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, pc_y, pc_z, snf_369, snf_378, snf_379, sod0_281, \
                         sod1_281, sof_468, sof_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = f_11 * snf_378[k]
                   + f_4 * sod0_281[k]
                   - f_5 * sod1_281[k]
                   + f_3 * pc_y[k] * sof_468[k];

        t_703[k] = f_11 * snf_379[k]
                   + f_3 * pc_y[k] * sof_469[k];

        t_704[k] = f_7 * snf_369[k]
                   + f_1 * sod0_281[k]
                   - f_2 * sod1_281[k]
                   + f_3 * pc_z[k] * sof_469[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, pc_x, pc_y, pc_z, snf_370, snf_380, snf_470, \
                         sod0_282, sod1_282, sof_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_8 * snf_470[k]
                   + f_1 * sod0_282[k]
                   - f_2 * sod1_282[k]
                   + f_3 * pc_x[k] * sof_470[k];

        t_706[k] = f_13 * snf_380[k]
                   + f_3 * pc_y[k] * sof_470[k];

        t_707[k] = f_8 * snf_370[k]
                   + f_3 * pc_z[k] * sof_470[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
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
    const auto f_10 = 4.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_660 = buffer.data(sng0 + 660);
    const auto *sng0_663 = buffer.data(sng0 + 663);
    const auto *sng0_665 = buffer.data(sng0 + 665);
    const auto *sng0_674 = buffer.data(sng0 + 674);

    const auto *snf_376 = buffer.data(snf + 376);
    const auto *snf_379 = buffer.data(snf + 379);
    const auto *snf_380 = buffer.data(snf + 380);
    const auto *snf_382 = buffer.data(snf + 382);
    const auto *snf_386 = buffer.data(snf + 386);
    const auto *snf_388 = buffer.data(snf + 388);
    const auto *snf_389 = buffer.data(snf + 389);
    const auto *snf_390 = buffer.data(snf + 390);
    const auto *snf_392 = buffer.data(snf + 392);
    const auto *snf_396 = buffer.data(snf + 396);
    const auto *snf_398 = buffer.data(snf + 398);
    const auto *snf_399 = buffer.data(snf + 399);
    const auto *snf_400 = buffer.data(snf + 400);
    const auto *snf_402 = buffer.data(snf + 402);
    const auto *snf_406 = buffer.data(snf + 406);
    const auto *snf_408 = buffer.data(snf + 408);
    const auto *snf_409 = buffer.data(snf + 409);
    const auto *snf_410 = buffer.data(snf + 410);
    const auto *snf_412 = buffer.data(snf + 412);
    const auto *snf_416 = buffer.data(snf + 416);
    const auto *snf_418 = buffer.data(snf + 418);
    const auto *snf_419 = buffer.data(snf + 419);
    const auto *snf_420 = buffer.data(snf + 420);
    const auto *snf_422 = buffer.data(snf + 422);
    const auto *snf_426 = buffer.data(snf + 426);
    const auto *snf_428 = buffer.data(snf + 428);
    const auto *snf_429 = buffer.data(snf + 429);
    const auto *snf_430 = buffer.data(snf + 430);
    const auto *snf_432 = buffer.data(snf + 432);
    const auto *snf_436 = buffer.data(snf + 436);
    const auto *snf_438 = buffer.data(snf + 438);
    const auto *snf_439 = buffer.data(snf + 439);
    const auto *snf_440 = buffer.data(snf + 440);
    const auto *snf_441 = buffer.data(snf + 441);
    const auto *snf_442 = buffer.data(snf + 442);
    const auto *snf_446 = buffer.data(snf + 446);
    const auto *snf_448 = buffer.data(snf + 448);
    const auto *snf_449 = buffer.data(snf + 449);
    const auto *snf_473 = buffer.data(snf + 473);
    const auto *snf_475 = buffer.data(snf + 475);
    const auto *snf_476 = buffer.data(snf + 476);
    const auto *snf_477 = buffer.data(snf + 477);
    const auto *snf_478 = buffer.data(snf + 478);
    const auto *snf_479 = buffer.data(snf + 479);
    const auto *snf_480 = buffer.data(snf + 480);
    const auto *snf_483 = buffer.data(snf + 483);
    const auto *snf_485 = buffer.data(snf + 485);
    const auto *snf_486 = buffer.data(snf + 486);
    const auto *snf_487 = buffer.data(snf + 487);
    const auto *snf_488 = buffer.data(snf + 488);
    const auto *snf_489 = buffer.data(snf + 489);
    const auto *snf_490 = buffer.data(snf + 490);
    const auto *snf_493 = buffer.data(snf + 493);
    const auto *snf_495 = buffer.data(snf + 495);
    const auto *snf_496 = buffer.data(snf + 496);
    const auto *snf_497 = buffer.data(snf + 497);
    const auto *snf_498 = buffer.data(snf + 498);
    const auto *snf_499 = buffer.data(snf + 499);
    const auto *snf_500 = buffer.data(snf + 500);
    const auto *snf_503 = buffer.data(snf + 503);
    const auto *snf_505 = buffer.data(snf + 505);
    const auto *snf_506 = buffer.data(snf + 506);
    const auto *snf_507 = buffer.data(snf + 507);
    const auto *snf_508 = buffer.data(snf + 508);
    const auto *snf_509 = buffer.data(snf + 509);
    const auto *snf_510 = buffer.data(snf + 510);
    const auto *snf_513 = buffer.data(snf + 513);
    const auto *snf_515 = buffer.data(snf + 515);
    const auto *snf_516 = buffer.data(snf + 516);
    const auto *snf_517 = buffer.data(snf + 517);
    const auto *snf_518 = buffer.data(snf + 518);
    const auto *snf_519 = buffer.data(snf + 519);
    const auto *snf_520 = buffer.data(snf + 520);
    const auto *snf_523 = buffer.data(snf + 523);
    const auto *snf_525 = buffer.data(snf + 525);
    const auto *snf_526 = buffer.data(snf + 526);
    const auto *snf_527 = buffer.data(snf + 527);
    const auto *snf_528 = buffer.data(snf + 528);
    const auto *snf_529 = buffer.data(snf + 529);
    const auto *snf_536 = buffer.data(snf + 536);
    const auto *snf_537 = buffer.data(snf + 537);
    const auto *snf_538 = buffer.data(snf + 538);
    const auto *snf_539 = buffer.data(snf + 539);
    const auto *snf_540 = buffer.data(snf + 540);
    const auto *snf_543 = buffer.data(snf + 543);
    const auto *snf_545 = buffer.data(snf + 545);
    const auto *snf_546 = buffer.data(snf + 546);
    const auto *snf_547 = buffer.data(snf + 547);
    const auto *snf_548 = buffer.data(snf + 548);
    const auto *snf_549 = buffer.data(snf + 549);

    const auto *sng1_660 = buffer.data(sng1 + 660);
    const auto *sng1_663 = buffer.data(sng1 + 663);
    const auto *sng1_665 = buffer.data(sng1 + 665);
    const auto *sng1_674 = buffer.data(sng1 + 674);

    const auto *sod0_285 = buffer.data(sod0 + 285);
    const auto *sod0_287 = buffer.data(sod0 + 287);
    const auto *sod0_288 = buffer.data(sod0 + 288);
    const auto *sod0_291 = buffer.data(sod0 + 291);
    const auto *sod0_293 = buffer.data(sod0 + 293);
    const auto *sod0_294 = buffer.data(sod0 + 294);
    const auto *sod0_297 = buffer.data(sod0 + 297);
    const auto *sod0_299 = buffer.data(sod0 + 299);
    const auto *sod0_300 = buffer.data(sod0 + 300);
    const auto *sod0_303 = buffer.data(sod0 + 303);
    const auto *sod0_305 = buffer.data(sod0 + 305);
    const auto *sod0_306 = buffer.data(sod0 + 306);
    const auto *sod0_309 = buffer.data(sod0 + 309);
    const auto *sod0_311 = buffer.data(sod0 + 311);
    const auto *sod0_312 = buffer.data(sod0 + 312);
    const auto *sod0_315 = buffer.data(sod0 + 315);
    const auto *sod0_317 = buffer.data(sod0 + 317);
    const auto *sod0_321 = buffer.data(sod0 + 321);
    const auto *sod0_323 = buffer.data(sod0 + 323);
    const auto *sod0_324 = buffer.data(sod0 + 324);
    const auto *sod0_327 = buffer.data(sod0 + 327);
    const auto *sod0_329 = buffer.data(sod0 + 329);

    const auto *sod1_285 = buffer.data(sod1 + 285);
    const auto *sod1_287 = buffer.data(sod1 + 287);
    const auto *sod1_288 = buffer.data(sod1 + 288);
    const auto *sod1_291 = buffer.data(sod1 + 291);
    const auto *sod1_293 = buffer.data(sod1 + 293);
    const auto *sod1_294 = buffer.data(sod1 + 294);
    const auto *sod1_297 = buffer.data(sod1 + 297);
    const auto *sod1_299 = buffer.data(sod1 + 299);
    const auto *sod1_300 = buffer.data(sod1 + 300);
    const auto *sod1_303 = buffer.data(sod1 + 303);
    const auto *sod1_305 = buffer.data(sod1 + 305);
    const auto *sod1_306 = buffer.data(sod1 + 306);
    const auto *sod1_309 = buffer.data(sod1 + 309);
    const auto *sod1_311 = buffer.data(sod1 + 311);
    const auto *sod1_312 = buffer.data(sod1 + 312);
    const auto *sod1_315 = buffer.data(sod1 + 315);
    const auto *sod1_317 = buffer.data(sod1 + 317);
    const auto *sod1_321 = buffer.data(sod1 + 321);
    const auto *sod1_323 = buffer.data(sod1 + 323);
    const auto *sod1_324 = buffer.data(sod1 + 324);
    const auto *sod1_327 = buffer.data(sod1 + 327);
    const auto *sod1_329 = buffer.data(sod1 + 329);

    const auto *sof_472 = buffer.data(sof + 472);
    const auto *sof_473 = buffer.data(sof + 473);
    const auto *sof_475 = buffer.data(sof + 475);
    const auto *sof_476 = buffer.data(sof + 476);
    const auto *sof_477 = buffer.data(sof + 477);
    const auto *sof_478 = buffer.data(sof + 478);
    const auto *sof_479 = buffer.data(sof + 479);
    const auto *sof_480 = buffer.data(sof + 480);
    const auto *sof_482 = buffer.data(sof + 482);
    const auto *sof_483 = buffer.data(sof + 483);
    const auto *sof_485 = buffer.data(sof + 485);
    const auto *sof_486 = buffer.data(sof + 486);
    const auto *sof_487 = buffer.data(sof + 487);
    const auto *sof_488 = buffer.data(sof + 488);
    const auto *sof_489 = buffer.data(sof + 489);
    const auto *sof_490 = buffer.data(sof + 490);
    const auto *sof_492 = buffer.data(sof + 492);
    const auto *sof_493 = buffer.data(sof + 493);
    const auto *sof_495 = buffer.data(sof + 495);
    const auto *sof_496 = buffer.data(sof + 496);
    const auto *sof_497 = buffer.data(sof + 497);
    const auto *sof_498 = buffer.data(sof + 498);
    const auto *sof_499 = buffer.data(sof + 499);
    const auto *sof_500 = buffer.data(sof + 500);
    const auto *sof_502 = buffer.data(sof + 502);
    const auto *sof_503 = buffer.data(sof + 503);
    const auto *sof_505 = buffer.data(sof + 505);
    const auto *sof_506 = buffer.data(sof + 506);
    const auto *sof_507 = buffer.data(sof + 507);
    const auto *sof_508 = buffer.data(sof + 508);
    const auto *sof_509 = buffer.data(sof + 509);
    const auto *sof_510 = buffer.data(sof + 510);
    const auto *sof_512 = buffer.data(sof + 512);
    const auto *sof_513 = buffer.data(sof + 513);
    const auto *sof_515 = buffer.data(sof + 515);
    const auto *sof_516 = buffer.data(sof + 516);
    const auto *sof_517 = buffer.data(sof + 517);
    const auto *sof_518 = buffer.data(sof + 518);
    const auto *sof_519 = buffer.data(sof + 519);
    const auto *sof_520 = buffer.data(sof + 520);
    const auto *sof_522 = buffer.data(sof + 522);
    const auto *sof_523 = buffer.data(sof + 523);
    const auto *sof_525 = buffer.data(sof + 525);
    const auto *sof_526 = buffer.data(sof + 526);
    const auto *sof_527 = buffer.data(sof + 527);
    const auto *sof_528 = buffer.data(sof + 528);
    const auto *sof_529 = buffer.data(sof + 529);
    const auto *sof_530 = buffer.data(sof + 530);
    const auto *sof_532 = buffer.data(sof + 532);
    const auto *sof_536 = buffer.data(sof + 536);
    const auto *sof_537 = buffer.data(sof + 537);
    const auto *sof_538 = buffer.data(sof + 538);
    const auto *sof_539 = buffer.data(sof + 539);
    const auto *sof_540 = buffer.data(sof + 540);
    const auto *sof_542 = buffer.data(sof + 542);
    const auto *sof_543 = buffer.data(sof + 543);
    const auto *sof_545 = buffer.data(sof + 545);
    const auto *sof_546 = buffer.data(sof + 546);
    const auto *sof_547 = buffer.data(sof + 547);
    const auto *sof_548 = buffer.data(sof + 548);
    const auto *sof_549 = buffer.data(sof + 549);

#pragma omp simd aligned(t_708, t_709, t_710, pc_x, pc_y, snf_382, snf_473, snf_475, sod0_285, \
                         sod0_287, sod1_285, sod1_287, sof_472, sof_473, \
                         sof_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_8 * snf_473[k]
                   + f_4 * sod0_285[k]
                   - f_5 * sod1_285[k]
                   + f_3 * pc_x[k] * sof_473[k];

        t_709[k] = f_13 * snf_382[k]
                   + f_3 * pc_y[k] * sof_472[k];

        t_710[k] = f_8 * snf_475[k]
                   + f_4 * sod0_287[k]
                   - f_5 * sod1_287[k]
                   + f_3 * pc_x[k] * sof_475[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, t_714, pc_x, snf_476, snf_477, snf_478, snf_479, \
                         sof_476, sof_477, sof_478, sof_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_8 * snf_476[k]
                   + f_3 * pc_x[k] * sof_476[k];

        t_712[k] = f_8 * snf_477[k]
                   + f_3 * pc_x[k] * sof_477[k];

        t_713[k] = f_8 * snf_478[k]
                   + f_3 * pc_x[k] * sof_478[k];

        t_714[k] = f_8 * snf_479[k]
                   + f_3 * pc_x[k] * sof_479[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, pc_y, pc_z, snf_376, snf_386, snf_388, sod0_285, \
                         sod0_287, sod1_285, sod1_287, sof_476, \
                         sof_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_13 * snf_386[k]
                   + f_1 * sod0_285[k]
                   - f_2 * sod1_285[k]
                   + f_3 * pc_y[k] * sof_476[k];

        t_716[k] = f_8 * snf_376[k]
                   + f_3 * pc_z[k] * sof_476[k];

        t_717[k] = f_13 * snf_388[k]
                   + f_4 * sod0_287[k]
                   - f_5 * sod1_287[k]
                   + f_3 * pc_y[k] * sof_478[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, pc_x, pc_y, pc_z, snf_379, snf_389, snf_480, \
                         sod0_287, sod0_288, sod1_287, sod1_288, sof_479, \
                         sof_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_13 * snf_389[k]
                   + f_3 * pc_y[k] * sof_479[k];

        t_719[k] = f_8 * snf_379[k]
                   + f_1 * sod0_287[k]
                   - f_2 * sod1_287[k]
                   + f_3 * pc_z[k] * sof_479[k];

        t_720[k] = f_8 * snf_480[k]
                   + f_1 * sod0_288[k]
                   - f_2 * sod1_288[k]
                   + f_3 * pc_x[k] * sof_480[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pc_x, pc_y, pc_z, snf_380, snf_390, \
                         snf_392, snf_483, sod0_291, sod1_291, sof_480, sof_482, \
                         sof_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_15 * snf_390[k]
                   + f_3 * pc_y[k] * sof_480[k];

        t_722[k] = f_12 * snf_380[k]
                   + f_3 * pc_z[k] * sof_480[k];

        t_723[k] = f_8 * snf_483[k]
                   + f_4 * sod0_291[k]
                   - f_5 * sod1_291[k]
                   + f_3 * pc_x[k] * sof_483[k];

        t_724[k] = f_15 * snf_392[k]
                   + f_3 * pc_y[k] * sof_482[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pc_x, snf_485, snf_486, snf_487, snf_488, \
                         sod0_293, sod1_293, sof_485, sof_486, sof_487, \
                         sof_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_8 * snf_485[k]
                   + f_4 * sod0_293[k]
                   - f_5 * sod1_293[k]
                   + f_3 * pc_x[k] * sof_485[k];

        t_726[k] = f_8 * snf_486[k]
                   + f_3 * pc_x[k] * sof_486[k];

        t_727[k] = f_8 * snf_487[k]
                   + f_3 * pc_x[k] * sof_487[k];

        t_728[k] = f_8 * snf_488[k]
                   + f_3 * pc_x[k] * sof_488[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_y, pc_z, snf_386, snf_396, snf_489, \
                         sod0_291, sod1_291, sof_486, sof_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_8 * snf_489[k]
                   + f_3 * pc_x[k] * sof_489[k];

        t_730[k] = f_15 * snf_396[k]
                   + f_1 * sod0_291[k]
                   - f_2 * sod1_291[k]
                   + f_3 * pc_y[k] * sof_486[k];

        t_731[k] = f_12 * snf_386[k]
                   + f_3 * pc_z[k] * sof_486[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pc_y, pc_z, snf_389, snf_398, snf_399, sod0_293, \
                         sod1_293, sof_488, sof_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_15 * snf_398[k]
                   + f_4 * sod0_293[k]
                   - f_5 * sod1_293[k]
                   + f_3 * pc_y[k] * sof_488[k];

        t_733[k] = f_15 * snf_399[k]
                   + f_3 * pc_y[k] * sof_489[k];

        t_734[k] = f_12 * snf_389[k]
                   + f_1 * sod0_293[k]
                   - f_2 * sod1_293[k]
                   + f_3 * pc_z[k] * sof_489[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pc_x, pc_y, pc_z, snf_390, snf_400, snf_490, \
                         sod0_294, sod1_294, sof_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_8 * snf_490[k]
                   + f_1 * sod0_294[k]
                   - f_2 * sod1_294[k]
                   + f_3 * pc_x[k] * sof_490[k];

        t_736[k] = f_16 * snf_400[k]
                   + f_3 * pc_y[k] * sof_490[k];

        t_737[k] = f_14 * snf_390[k]
                   + f_3 * pc_z[k] * sof_490[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pc_x, pc_y, snf_402, snf_493, snf_495, sod0_297, \
                         sod0_299, sod1_297, sod1_299, sof_492, sof_493, \
                         sof_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_8 * snf_493[k]
                   + f_4 * sod0_297[k]
                   - f_5 * sod1_297[k]
                   + f_3 * pc_x[k] * sof_493[k];

        t_739[k] = f_16 * snf_402[k]
                   + f_3 * pc_y[k] * sof_492[k];

        t_740[k] = f_8 * snf_495[k]
                   + f_4 * sod0_299[k]
                   - f_5 * sod1_299[k]
                   + f_3 * pc_x[k] * sof_495[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pc_x, snf_496, snf_497, snf_498, snf_499, \
                         sof_496, sof_497, sof_498, sof_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_8 * snf_496[k]
                   + f_3 * pc_x[k] * sof_496[k];

        t_742[k] = f_8 * snf_497[k]
                   + f_3 * pc_x[k] * sof_497[k];

        t_743[k] = f_8 * snf_498[k]
                   + f_3 * pc_x[k] * sof_498[k];

        t_744[k] = f_8 * snf_499[k]
                   + f_3 * pc_x[k] * sof_499[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, pc_y, pc_z, snf_396, snf_406, snf_408, sod0_297, \
                         sod0_299, sod1_297, sod1_299, sof_496, \
                         sof_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_16 * snf_406[k]
                   + f_1 * sod0_297[k]
                   - f_2 * sod1_297[k]
                   + f_3 * pc_y[k] * sof_496[k];

        t_746[k] = f_14 * snf_396[k]
                   + f_3 * pc_z[k] * sof_496[k];

        t_747[k] = f_16 * snf_408[k]
                   + f_4 * sod0_299[k]
                   - f_5 * sod1_299[k]
                   + f_3 * pc_y[k] * sof_498[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, pc_x, pc_y, pc_z, snf_399, snf_409, snf_500, \
                         sod0_299, sod0_300, sod1_299, sod1_300, sof_499, \
                         sof_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_16 * snf_409[k]
                   + f_3 * pc_y[k] * sof_499[k];

        t_749[k] = f_14 * snf_399[k]
                   + f_1 * sod0_299[k]
                   - f_2 * sod1_299[k]
                   + f_3 * pc_z[k] * sof_499[k];

        t_750[k] = f_8 * snf_500[k]
                   + f_1 * sod0_300[k]
                   - f_2 * sod1_300[k]
                   + f_3 * pc_x[k] * sof_500[k];
    }

#pragma omp simd aligned(t_751, t_752, t_753, t_754, pc_x, pc_y, pc_z, snf_400, snf_410, \
                         snf_412, snf_503, sod0_303, sod1_303, sof_500, sof_502, \
                         sof_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_751[k] = f_14 * snf_410[k]
                   + f_3 * pc_y[k] * sof_500[k];

        t_752[k] = f_16 * snf_400[k]
                   + f_3 * pc_z[k] * sof_500[k];

        t_753[k] = f_8 * snf_503[k]
                   + f_4 * sod0_303[k]
                   - f_5 * sod1_303[k]
                   + f_3 * pc_x[k] * sof_503[k];

        t_754[k] = f_14 * snf_412[k]
                   + f_3 * pc_y[k] * sof_502[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pc_x, snf_505, snf_506, snf_507, snf_508, \
                         sod0_305, sod1_305, sof_505, sof_506, sof_507, \
                         sof_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_8 * snf_505[k]
                   + f_4 * sod0_305[k]
                   - f_5 * sod1_305[k]
                   + f_3 * pc_x[k] * sof_505[k];

        t_756[k] = f_8 * snf_506[k]
                   + f_3 * pc_x[k] * sof_506[k];

        t_757[k] = f_8 * snf_507[k]
                   + f_3 * pc_x[k] * sof_507[k];

        t_758[k] = f_8 * snf_508[k]
                   + f_3 * pc_x[k] * sof_508[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, pc_z, snf_406, snf_416, snf_509, \
                         sod0_303, sod1_303, sof_506, sof_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_8 * snf_509[k]
                   + f_3 * pc_x[k] * sof_509[k];

        t_760[k] = f_14 * snf_416[k]
                   + f_1 * sod0_303[k]
                   - f_2 * sod1_303[k]
                   + f_3 * pc_y[k] * sof_506[k];

        t_761[k] = f_16 * snf_406[k]
                   + f_3 * pc_z[k] * sof_506[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_y, pc_z, snf_409, snf_418, snf_419, sod0_305, \
                         sod1_305, sof_508, sof_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_14 * snf_418[k]
                   + f_4 * sod0_305[k]
                   - f_5 * sod1_305[k]
                   + f_3 * pc_y[k] * sof_508[k];

        t_763[k] = f_14 * snf_419[k]
                   + f_3 * pc_y[k] * sof_509[k];

        t_764[k] = f_16 * snf_409[k]
                   + f_1 * sod0_305[k]
                   - f_2 * sod1_305[k]
                   + f_3 * pc_z[k] * sof_509[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_y, pc_z, snf_410, snf_420, snf_510, \
                         sod0_306, sod1_306, sof_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_8 * snf_510[k]
                   + f_1 * sod0_306[k]
                   - f_2 * sod1_306[k]
                   + f_3 * pc_x[k] * sof_510[k];

        t_766[k] = f_12 * snf_420[k]
                   + f_3 * pc_y[k] * sof_510[k];

        t_767[k] = f_15 * snf_410[k]
                   + f_3 * pc_z[k] * sof_510[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_x, pc_y, snf_422, snf_513, snf_515, sod0_309, \
                         sod0_311, sod1_309, sod1_311, sof_512, sof_513, \
                         sof_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_8 * snf_513[k]
                   + f_4 * sod0_309[k]
                   - f_5 * sod1_309[k]
                   + f_3 * pc_x[k] * sof_513[k];

        t_769[k] = f_12 * snf_422[k]
                   + f_3 * pc_y[k] * sof_512[k];

        t_770[k] = f_8 * snf_515[k]
                   + f_4 * sod0_311[k]
                   - f_5 * sod1_311[k]
                   + f_3 * pc_x[k] * sof_515[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, pc_x, snf_516, snf_517, snf_518, snf_519, \
                         sof_516, sof_517, sof_518, sof_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_8 * snf_516[k]
                   + f_3 * pc_x[k] * sof_516[k];

        t_772[k] = f_8 * snf_517[k]
                   + f_3 * pc_x[k] * sof_517[k];

        t_773[k] = f_8 * snf_518[k]
                   + f_3 * pc_x[k] * sof_518[k];

        t_774[k] = f_8 * snf_519[k]
                   + f_3 * pc_x[k] * sof_519[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, pc_y, pc_z, snf_416, snf_426, snf_428, sod0_309, \
                         sod0_311, sod1_309, sod1_311, sof_516, \
                         sof_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_12 * snf_426[k]
                   + f_1 * sod0_309[k]
                   - f_2 * sod1_309[k]
                   + f_3 * pc_y[k] * sof_516[k];

        t_776[k] = f_15 * snf_416[k]
                   + f_3 * pc_z[k] * sof_516[k];

        t_777[k] = f_12 * snf_428[k]
                   + f_4 * sod0_311[k]
                   - f_5 * sod1_311[k]
                   + f_3 * pc_y[k] * sof_518[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, pc_x, pc_y, pc_z, snf_419, snf_429, snf_520, \
                         sod0_311, sod0_312, sod1_311, sod1_312, sof_519, \
                         sof_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_12 * snf_429[k]
                   + f_3 * pc_y[k] * sof_519[k];

        t_779[k] = f_15 * snf_419[k]
                   + f_1 * sod0_311[k]
                   - f_2 * sod1_311[k]
                   + f_3 * pc_z[k] * sof_519[k];

        t_780[k] = f_8 * snf_520[k]
                   + f_1 * sod0_312[k]
                   - f_2 * sod1_312[k]
                   + f_3 * pc_x[k] * sof_520[k];
    }

#pragma omp simd aligned(t_781, t_782, t_783, t_784, pc_x, pc_y, pc_z, snf_420, snf_430, \
                         snf_432, snf_523, sod0_315, sod1_315, sof_520, sof_522, \
                         sof_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_781[k] = f_8 * snf_430[k]
                   + f_3 * pc_y[k] * sof_520[k];

        t_782[k] = f_13 * snf_420[k]
                   + f_3 * pc_z[k] * sof_520[k];

        t_783[k] = f_8 * snf_523[k]
                   + f_4 * sod0_315[k]
                   - f_5 * sod1_315[k]
                   + f_3 * pc_x[k] * sof_523[k];

        t_784[k] = f_8 * snf_432[k]
                   + f_3 * pc_y[k] * sof_522[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, pc_x, snf_525, snf_526, snf_527, snf_528, \
                         sod0_317, sod1_317, sof_525, sof_526, sof_527, \
                         sof_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_8 * snf_525[k]
                   + f_4 * sod0_317[k]
                   - f_5 * sod1_317[k]
                   + f_3 * pc_x[k] * sof_525[k];

        t_786[k] = f_8 * snf_526[k]
                   + f_3 * pc_x[k] * sof_526[k];

        t_787[k] = f_8 * snf_527[k]
                   + f_3 * pc_x[k] * sof_527[k];

        t_788[k] = f_8 * snf_528[k]
                   + f_3 * pc_x[k] * sof_528[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, pc_x, pc_y, pc_z, snf_426, snf_436, snf_529, \
                         sod0_315, sod1_315, sof_526, sof_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = f_8 * snf_529[k]
                   + f_3 * pc_x[k] * sof_529[k];

        t_790[k] = f_8 * snf_436[k]
                   + f_1 * sod0_315[k]
                   - f_2 * sod1_315[k]
                   + f_3 * pc_y[k] * sof_526[k];

        t_791[k] = f_13 * snf_426[k]
                   + f_3 * pc_z[k] * sof_526[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pb_y, pc_y, pc_z, sng0_660, snf_429, \
                         snf_438, snf_439, sng1_660, sod0_317, sod1_317, sof_528, \
                         sof_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_8 * snf_438[k]
                   + f_4 * sod0_317[k]
                   - f_5 * sod1_317[k]
                   + f_3 * pc_y[k] * sof_528[k];

        t_793[k] = f_8 * snf_439[k]
                   + f_3 * pc_y[k] * sof_529[k];

        t_794[k] = f_13 * snf_429[k]
                   + f_1 * sod0_317[k]
                   - f_2 * sod1_317[k]
                   + f_3 * pc_z[k] * sof_529[k];

        t_795[k] = pb_y[k] * sng0_660[k]
                   - f_6 * pc_y[k] * sng1_660[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pb_y, pc_y, pc_z, sng0_663, snf_430, \
                         snf_440, snf_441, snf_442, sng1_663, sof_530, \
                         sof_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_7 * snf_440[k]
                   + f_3 * pc_y[k] * sof_530[k];

        t_797[k] = f_11 * snf_430[k]
                   + f_3 * pc_z[k] * sof_530[k];

        t_798[k] = pb_y[k] * sng0_663[k]
                   + f_8 * snf_441[k]
                   - f_6 * pc_y[k] * sng1_663[k];

        t_799[k] = f_7 * snf_442[k]
                   + f_3 * pc_y[k] * sof_532[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, pb_y, pc_x, pc_y, sng0_665, snf_536, \
                         snf_537, snf_538, sng1_665, sof_536, sof_537, \
                         sof_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = pb_y[k] * sng0_665[k]
                   - f_6 * pc_y[k] * sng1_665[k];

        t_801[k] = f_8 * snf_536[k]
                   + f_3 * pc_x[k] * sof_536[k];

        t_802[k] = f_8 * snf_537[k]
                   + f_3 * pc_x[k] * sof_537[k];

        t_803[k] = f_8 * snf_538[k]
                   + f_3 * pc_x[k] * sof_538[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pc_x, pc_y, pc_z, snf_436, snf_446, snf_539, \
                         sod0_321, sod1_321, sof_536, sof_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_8 * snf_539[k]
                   + f_3 * pc_x[k] * sof_539[k];

        t_805[k] = f_7 * snf_446[k]
                   + f_1 * sod0_321[k]
                   - f_2 * sod1_321[k]
                   + f_3 * pc_y[k] * sof_536[k];

        t_806[k] = f_11 * snf_436[k]
                   + f_3 * pc_z[k] * sof_536[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, pb_y, pc_y, sng0_674, snf_448, snf_449, \
                         sng1_674, sod0_323, sod1_323, sof_538, \
                         sof_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = f_7 * snf_448[k]
                   + f_4 * sod0_323[k]
                   - f_5 * sod1_323[k]
                   + f_3 * pc_y[k] * sof_538[k];

        t_808[k] = f_7 * snf_449[k]
                   + f_3 * pc_y[k] * sof_539[k];

        t_809[k] = pb_y[k] * sng0_674[k]
                   - f_6 * pc_y[k] * sng1_674[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pc_x, pc_y, pc_z, snf_440, snf_540, \
                         snf_543, sod0_324, sod0_327, sod1_324, sod1_327, sof_540, \
                         sof_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_8 * snf_540[k]
                   + f_1 * sod0_324[k]
                   - f_2 * sod1_324[k]
                   + f_3 * pc_x[k] * sof_540[k];

        t_811[k] = f_3 * pc_y[k] * sof_540[k];

        t_812[k] = f_10 * snf_440[k]
                   + f_3 * pc_z[k] * sof_540[k];

        t_813[k] = f_8 * snf_543[k]
                   + f_4 * sod0_327[k]
                   - f_5 * sod1_327[k]
                   + f_3 * pc_x[k] * sof_543[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, t_817, pc_x, pc_y, snf_545, snf_546, snf_547, \
                         sod0_329, sod1_329, sof_542, sof_545, sof_546, \
                         sof_547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_3 * pc_y[k] * sof_542[k];

        t_815[k] = f_8 * snf_545[k]
                   + f_4 * sod0_329[k]
                   - f_5 * sod1_329[k]
                   + f_3 * pc_x[k] * sof_545[k];

        t_816[k] = f_8 * snf_546[k]
                   + f_3 * pc_x[k] * sof_546[k];

        t_817[k] = f_8 * snf_547[k]
                   + f_3 * pc_x[k] * sof_547[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, pc_x, pc_y, pc_z, snf_446, snf_548, \
                         snf_549, sod0_327, sod1_327, sof_546, sof_548, \
                         sof_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_8 * snf_548[k]
                   + f_3 * pc_x[k] * sof_548[k];

        t_819[k] = f_8 * snf_549[k]
                   + f_3 * pc_x[k] * sof_549[k];

        t_820[k] = f_1 * sod0_327[k]
                   - f_2 * sod1_327[k]
                   + f_3 * pc_y[k] * sof_546[k];

        t_821[k] = f_10 * snf_446[k]
                   + f_3 * pc_z[k] * sof_546[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
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
    const auto f_9 = 5.0 / q;
    const auto f_10 = 4.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_675 = buffer.data(sng0 + 675);
    const auto *sng0_678 = buffer.data(sng0 + 678);
    const auto *sng0_825 = buffer.data(sng0 + 825);
    const auto *sng0_828 = buffer.data(sng0 + 828);
    const auto *sng0_830 = buffer.data(sng0 + 830);
    const auto *sng0_835 = buffer.data(sng0 + 835);
    const auto *sng0_837 = buffer.data(sng0 + 837);
    const auto *sng0_839 = buffer.data(sng0 + 839);
    const auto *sng0_845 = buffer.data(sng0 + 845);
    const auto *sng0_850 = buffer.data(sng0 + 850);
    const auto *sng0_852 = buffer.data(sng0 + 852);
    const auto *sng0_854 = buffer.data(sng0 + 854);
    const auto *sng0_855 = buffer.data(sng0 + 855);
    const auto *sng0_858 = buffer.data(sng0 + 858);
    const auto *sng0_860 = buffer.data(sng0 + 860);
    const auto *sng0_865 = buffer.data(sng0 + 865);
    const auto *sng0_867 = buffer.data(sng0 + 867);
    const auto *sng0_869 = buffer.data(sng0 + 869);
    const auto *sng0_870 = buffer.data(sng0 + 870);
    const auto *sng0_873 = buffer.data(sng0 + 873);
    const auto *sng0_875 = buffer.data(sng0 + 875);
    const auto *sng0_880 = buffer.data(sng0 + 880);
    const auto *sng0_882 = buffer.data(sng0 + 882);
    const auto *sng0_884 = buffer.data(sng0 + 884);
    const auto *sng0_885 = buffer.data(sng0 + 885);
    const auto *sng0_888 = buffer.data(sng0 + 888);
    const auto *sng0_890 = buffer.data(sng0 + 890);
    const auto *sng0_895 = buffer.data(sng0 + 895);
    const auto *sng0_897 = buffer.data(sng0 + 897);
    const auto *sng0_899 = buffer.data(sng0 + 899);
    const auto *sng0_900 = buffer.data(sng0 + 900);
    const auto *sng0_903 = buffer.data(sng0 + 903);
    const auto *sng0_905 = buffer.data(sng0 + 905);
    const auto *sng0_910 = buffer.data(sng0 + 910);
    const auto *sng0_912 = buffer.data(sng0 + 912);
    const auto *sng0_914 = buffer.data(sng0 + 914);
    const auto *sng0_915 = buffer.data(sng0 + 915);
    const auto *sng0_918 = buffer.data(sng0 + 918);
    const auto *sng0_920 = buffer.data(sng0 + 920);
    const auto *sng0_925 = buffer.data(sng0 + 925);
    const auto *sng0_927 = buffer.data(sng0 + 927);
    const auto *sng0_929 = buffer.data(sng0 + 929);
    const auto *sng0_930 = buffer.data(sng0 + 930);
    const auto *sng0_933 = buffer.data(sng0 + 933);
    const auto *sng0_935 = buffer.data(sng0 + 935);
    const auto *sng0_940 = buffer.data(sng0 + 940);
    const auto *sng0_942 = buffer.data(sng0 + 942);
    const auto *sng0_944 = buffer.data(sng0 + 944);
    const auto *sng0_945 = buffer.data(sng0 + 945);

    const auto *snf_449 = buffer.data(snf + 449);
    const auto *snf_450 = buffer.data(snf + 450);
    const auto *snf_452 = buffer.data(snf + 452);
    const auto *snf_456 = buffer.data(snf + 456);
    const auto *snf_459 = buffer.data(snf + 459);
    const auto *snf_460 = buffer.data(snf + 460);
    const auto *snf_462 = buffer.data(snf + 462);
    const auto *snf_466 = buffer.data(snf + 466);
    const auto *snf_469 = buffer.data(snf + 469);
    const auto *snf_470 = buffer.data(snf + 470);
    const auto *snf_472 = buffer.data(snf + 472);
    const auto *snf_476 = buffer.data(snf + 476);
    const auto *snf_479 = buffer.data(snf + 479);
    const auto *snf_480 = buffer.data(snf + 480);
    const auto *snf_482 = buffer.data(snf + 482);
    const auto *snf_486 = buffer.data(snf + 486);
    const auto *snf_489 = buffer.data(snf + 489);
    const auto *snf_490 = buffer.data(snf + 490);
    const auto *snf_492 = buffer.data(snf + 492);
    const auto *snf_496 = buffer.data(snf + 496);
    const auto *snf_499 = buffer.data(snf + 499);
    const auto *snf_500 = buffer.data(snf + 500);
    const auto *snf_502 = buffer.data(snf + 502);
    const auto *snf_506 = buffer.data(snf + 506);
    const auto *snf_509 = buffer.data(snf + 509);
    const auto *snf_510 = buffer.data(snf + 510);
    const auto *snf_512 = buffer.data(snf + 512);
    const auto *snf_516 = buffer.data(snf + 516);
    const auto *snf_519 = buffer.data(snf + 519);
    const auto *snf_520 = buffer.data(snf + 520);
    const auto *snf_522 = buffer.data(snf + 522);
    const auto *snf_529 = buffer.data(snf + 529);
    const auto *snf_530 = buffer.data(snf + 530);
    const auto *snf_550 = buffer.data(snf + 550);
    const auto *snf_553 = buffer.data(snf + 553);
    const auto *snf_555 = buffer.data(snf + 555);
    const auto *snf_556 = buffer.data(snf + 556);
    const auto *snf_557 = buffer.data(snf + 557);
    const auto *snf_558 = buffer.data(snf + 558);
    const auto *snf_559 = buffer.data(snf + 559);
    const auto *snf_565 = buffer.data(snf + 565);
    const auto *snf_566 = buffer.data(snf + 566);
    const auto *snf_567 = buffer.data(snf + 567);
    const auto *snf_568 = buffer.data(snf + 568);
    const auto *snf_569 = buffer.data(snf + 569);
    const auto *snf_570 = buffer.data(snf + 570);
    const auto *snf_573 = buffer.data(snf + 573);
    const auto *snf_575 = buffer.data(snf + 575);
    const auto *snf_576 = buffer.data(snf + 576);
    const auto *snf_577 = buffer.data(snf + 577);
    const auto *snf_578 = buffer.data(snf + 578);
    const auto *snf_579 = buffer.data(snf + 579);
    const auto *snf_580 = buffer.data(snf + 580);
    const auto *snf_583 = buffer.data(snf + 583);
    const auto *snf_585 = buffer.data(snf + 585);
    const auto *snf_586 = buffer.data(snf + 586);
    const auto *snf_587 = buffer.data(snf + 587);
    const auto *snf_588 = buffer.data(snf + 588);
    const auto *snf_589 = buffer.data(snf + 589);
    const auto *snf_590 = buffer.data(snf + 590);
    const auto *snf_593 = buffer.data(snf + 593);
    const auto *snf_595 = buffer.data(snf + 595);
    const auto *snf_596 = buffer.data(snf + 596);
    const auto *snf_597 = buffer.data(snf + 597);
    const auto *snf_598 = buffer.data(snf + 598);
    const auto *snf_599 = buffer.data(snf + 599);
    const auto *snf_600 = buffer.data(snf + 600);
    const auto *snf_603 = buffer.data(snf + 603);
    const auto *snf_605 = buffer.data(snf + 605);
    const auto *snf_606 = buffer.data(snf + 606);
    const auto *snf_607 = buffer.data(snf + 607);
    const auto *snf_608 = buffer.data(snf + 608);
    const auto *snf_609 = buffer.data(snf + 609);
    const auto *snf_610 = buffer.data(snf + 610);
    const auto *snf_613 = buffer.data(snf + 613);
    const auto *snf_615 = buffer.data(snf + 615);
    const auto *snf_616 = buffer.data(snf + 616);
    const auto *snf_617 = buffer.data(snf + 617);
    const auto *snf_618 = buffer.data(snf + 618);
    const auto *snf_619 = buffer.data(snf + 619);
    const auto *snf_620 = buffer.data(snf + 620);
    const auto *snf_623 = buffer.data(snf + 623);
    const auto *snf_625 = buffer.data(snf + 625);
    const auto *snf_626 = buffer.data(snf + 626);
    const auto *snf_627 = buffer.data(snf + 627);
    const auto *snf_628 = buffer.data(snf + 628);
    const auto *snf_629 = buffer.data(snf + 629);
    const auto *snf_630 = buffer.data(snf + 630);

    const auto *sng1_675 = buffer.data(sng1 + 675);
    const auto *sng1_678 = buffer.data(sng1 + 678);
    const auto *sng1_825 = buffer.data(sng1 + 825);
    const auto *sng1_828 = buffer.data(sng1 + 828);
    const auto *sng1_830 = buffer.data(sng1 + 830);
    const auto *sng1_835 = buffer.data(sng1 + 835);
    const auto *sng1_837 = buffer.data(sng1 + 837);
    const auto *sng1_839 = buffer.data(sng1 + 839);
    const auto *sng1_845 = buffer.data(sng1 + 845);
    const auto *sng1_850 = buffer.data(sng1 + 850);
    const auto *sng1_852 = buffer.data(sng1 + 852);
    const auto *sng1_854 = buffer.data(sng1 + 854);
    const auto *sng1_855 = buffer.data(sng1 + 855);
    const auto *sng1_858 = buffer.data(sng1 + 858);
    const auto *sng1_860 = buffer.data(sng1 + 860);
    const auto *sng1_865 = buffer.data(sng1 + 865);
    const auto *sng1_867 = buffer.data(sng1 + 867);
    const auto *sng1_869 = buffer.data(sng1 + 869);
    const auto *sng1_870 = buffer.data(sng1 + 870);
    const auto *sng1_873 = buffer.data(sng1 + 873);
    const auto *sng1_875 = buffer.data(sng1 + 875);
    const auto *sng1_880 = buffer.data(sng1 + 880);
    const auto *sng1_882 = buffer.data(sng1 + 882);
    const auto *sng1_884 = buffer.data(sng1 + 884);
    const auto *sng1_885 = buffer.data(sng1 + 885);
    const auto *sng1_888 = buffer.data(sng1 + 888);
    const auto *sng1_890 = buffer.data(sng1 + 890);
    const auto *sng1_895 = buffer.data(sng1 + 895);
    const auto *sng1_897 = buffer.data(sng1 + 897);
    const auto *sng1_899 = buffer.data(sng1 + 899);
    const auto *sng1_900 = buffer.data(sng1 + 900);
    const auto *sng1_903 = buffer.data(sng1 + 903);
    const auto *sng1_905 = buffer.data(sng1 + 905);
    const auto *sng1_910 = buffer.data(sng1 + 910);
    const auto *sng1_912 = buffer.data(sng1 + 912);
    const auto *sng1_914 = buffer.data(sng1 + 914);
    const auto *sng1_915 = buffer.data(sng1 + 915);
    const auto *sng1_918 = buffer.data(sng1 + 918);
    const auto *sng1_920 = buffer.data(sng1 + 920);
    const auto *sng1_925 = buffer.data(sng1 + 925);
    const auto *sng1_927 = buffer.data(sng1 + 927);
    const auto *sng1_929 = buffer.data(sng1 + 929);
    const auto *sng1_930 = buffer.data(sng1 + 930);
    const auto *sng1_933 = buffer.data(sng1 + 933);
    const auto *sng1_935 = buffer.data(sng1 + 935);
    const auto *sng1_940 = buffer.data(sng1 + 940);
    const auto *sng1_942 = buffer.data(sng1 + 942);
    const auto *sng1_944 = buffer.data(sng1 + 944);
    const auto *sng1_945 = buffer.data(sng1 + 945);

    const auto *sod0_329 = buffer.data(sod0 + 329);

    const auto *sod1_329 = buffer.data(sod1 + 329);

    const auto *sof_548 = buffer.data(sof + 548);
    const auto *sof_549 = buffer.data(sof + 549);
    const auto *sof_550 = buffer.data(sof + 550);
    const auto *sof_552 = buffer.data(sof + 552);
    const auto *sof_556 = buffer.data(sof + 556);
    const auto *sof_557 = buffer.data(sof + 557);
    const auto *sof_558 = buffer.data(sof + 558);
    const auto *sof_559 = buffer.data(sof + 559);
    const auto *sof_560 = buffer.data(sof + 560);
    const auto *sof_562 = buffer.data(sof + 562);
    const auto *sof_566 = buffer.data(sof + 566);
    const auto *sof_567 = buffer.data(sof + 567);
    const auto *sof_568 = buffer.data(sof + 568);
    const auto *sof_569 = buffer.data(sof + 569);
    const auto *sof_570 = buffer.data(sof + 570);
    const auto *sof_572 = buffer.data(sof + 572);
    const auto *sof_576 = buffer.data(sof + 576);
    const auto *sof_577 = buffer.data(sof + 577);
    const auto *sof_578 = buffer.data(sof + 578);
    const auto *sof_579 = buffer.data(sof + 579);
    const auto *sof_580 = buffer.data(sof + 580);
    const auto *sof_582 = buffer.data(sof + 582);
    const auto *sof_586 = buffer.data(sof + 586);
    const auto *sof_587 = buffer.data(sof + 587);
    const auto *sof_588 = buffer.data(sof + 588);
    const auto *sof_589 = buffer.data(sof + 589);
    const auto *sof_590 = buffer.data(sof + 590);
    const auto *sof_592 = buffer.data(sof + 592);
    const auto *sof_596 = buffer.data(sof + 596);
    const auto *sof_597 = buffer.data(sof + 597);
    const auto *sof_598 = buffer.data(sof + 598);
    const auto *sof_599 = buffer.data(sof + 599);
    const auto *sof_600 = buffer.data(sof + 600);
    const auto *sof_602 = buffer.data(sof + 602);
    const auto *sof_606 = buffer.data(sof + 606);
    const auto *sof_607 = buffer.data(sof + 607);
    const auto *sof_608 = buffer.data(sof + 608);
    const auto *sof_609 = buffer.data(sof + 609);
    const auto *sof_610 = buffer.data(sof + 610);
    const auto *sof_612 = buffer.data(sof + 612);
    const auto *sof_616 = buffer.data(sof + 616);
    const auto *sof_617 = buffer.data(sof + 617);
    const auto *sof_618 = buffer.data(sof + 618);
    const auto *sof_619 = buffer.data(sof + 619);
    const auto *sof_620 = buffer.data(sof + 620);
    const auto *sof_622 = buffer.data(sof + 622);
    const auto *sof_626 = buffer.data(sof + 626);
    const auto *sof_627 = buffer.data(sof + 627);
    const auto *sof_628 = buffer.data(sof + 628);
    const auto *sof_629 = buffer.data(sof + 629);
    const auto *sof_630 = buffer.data(sof + 630);

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pb_x, pc_x, pc_y, pc_z, sng0_825, \
                         snf_449, snf_550, sng1_825, sod0_329, sod1_329, sof_548, \
                         sof_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_4 * sod0_329[k]
                   - f_5 * sod1_329[k]
                   + f_3 * pc_y[k] * sof_548[k];

        t_823[k] = f_3 * pc_y[k] * sof_549[k];

        t_824[k] = f_10 * snf_449[k]
                   + f_1 * sod0_329[k]
                   - f_2 * sod1_329[k]
                   + f_3 * pc_z[k] * sof_549[k];

        t_825[k] = pb_x[k] * sng0_825[k]
                   + f_14 * snf_550[k]
                   - f_6 * pc_x[k] * sng1_825[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pb_x, pc_x, pc_y, pc_z, sng0_828, \
                         snf_450, snf_452, snf_553, sng1_828, sof_550, \
                         sof_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_9 * snf_450[k]
                   + f_3 * pc_y[k] * sof_550[k];

        t_827[k] = f_3 * pc_z[k] * sof_550[k];

        t_828[k] = pb_x[k] * sng0_828[k]
                   + f_8 * snf_553[k]
                   - f_6 * pc_x[k] * sng1_828[k];

        t_829[k] = f_9 * snf_452[k]
                   + f_3 * pc_y[k] * sof_552[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pb_x, pc_x, sng0_830, snf_555, snf_556, \
                         snf_557, snf_558, sng1_830, sof_556, sof_557, \
                         sof_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = pb_x[k] * sng0_830[k]
                   + f_8 * snf_555[k]
                   - f_6 * pc_x[k] * sng1_830[k];

        t_831[k] = f_7 * snf_556[k]
                   + f_3 * pc_x[k] * sof_556[k];

        t_832[k] = f_7 * snf_557[k]
                   + f_3 * pc_x[k] * sof_557[k];

        t_833[k] = f_7 * snf_558[k]
                   + f_3 * pc_x[k] * sof_558[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, pb_x, pc_x, pc_z, sng0_835, sng0_837, \
                         snf_559, sng1_835, sng1_837, sof_556, \
                         sof_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_7 * snf_559[k]
                   + f_3 * pc_x[k] * sof_559[k];

        t_835[k] = pb_x[k] * sng0_835[k]
                   - f_6 * pc_x[k] * sng1_835[k];

        t_836[k] = f_3 * pc_z[k] * sof_556[k];

        t_837[k] = pb_x[k] * sng0_837[k]
                   - f_6 * pc_x[k] * sng1_837[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, pb_x, pb_z, pc_x, pc_y, pc_z, sng0_675, \
                         sng0_839, snf_459, sng1_675, sng1_839, \
                         sof_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_9 * snf_459[k]
                   + f_3 * pc_y[k] * sof_559[k];

        t_839[k] = pb_x[k] * sng0_839[k]
                   - f_6 * pc_x[k] * sng1_839[k];

        t_840[k] = pb_z[k] * sng0_675[k]
                   - f_6 * pc_z[k] * sng1_675[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pb_z, pc_y, pc_z, sng0_678, snf_450, \
                         snf_460, snf_462, sng1_678, sof_560, sof_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_10 * snf_460[k]
                   + f_3 * pc_y[k] * sof_560[k];

        t_842[k] = f_7 * snf_450[k]
                   + f_3 * pc_z[k] * sof_560[k];

        t_843[k] = pb_z[k] * sng0_678[k]
                   - f_6 * pc_z[k] * sng1_678[k];

        t_844[k] = f_10 * snf_462[k]
                   + f_3 * pc_y[k] * sof_562[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, pb_x, pc_x, sng0_845, snf_565, snf_566, \
                         snf_567, snf_568, sng1_845, sof_566, sof_567, \
                         sof_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = pb_x[k] * sng0_845[k]
                   + f_8 * snf_565[k]
                   - f_6 * pc_x[k] * sng1_845[k];

        t_846[k] = f_7 * snf_566[k]
                   + f_3 * pc_x[k] * sof_566[k];

        t_847[k] = f_7 * snf_567[k]
                   + f_3 * pc_x[k] * sof_567[k];

        t_848[k] = f_7 * snf_568[k]
                   + f_3 * pc_x[k] * sof_568[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, pb_x, pc_x, pc_z, sng0_850, sng0_852, \
                         snf_456, snf_569, sng1_850, sng1_852, sof_566, \
                         sof_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_7 * snf_569[k]
                   + f_3 * pc_x[k] * sof_569[k];

        t_850[k] = pb_x[k] * sng0_850[k]
                   - f_6 * pc_x[k] * sng1_850[k];

        t_851[k] = f_7 * snf_456[k]
                   + f_3 * pc_z[k] * sof_566[k];

        t_852[k] = pb_x[k] * sng0_852[k]
                   - f_6 * pc_x[k] * sng1_852[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, pb_x, pc_x, pc_y, sng0_854, sng0_855, \
                         snf_469, snf_470, snf_570, sng1_854, sng1_855, sof_569, \
                         sof_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_10 * snf_469[k]
                   + f_3 * pc_y[k] * sof_569[k];

        t_854[k] = pb_x[k] * sng0_854[k]
                   - f_6 * pc_x[k] * sng1_854[k];

        t_855[k] = pb_x[k] * sng0_855[k]
                   + f_14 * snf_570[k]
                   - f_6 * pc_x[k] * sng1_855[k];

        t_856[k] = f_11 * snf_470[k]
                   + f_3 * pc_y[k] * sof_570[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, pb_x, pc_x, pc_y, pc_z, sng0_858, snf_460, \
                         snf_472, snf_573, sng1_858, sof_570, sof_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = f_8 * snf_460[k]
                   + f_3 * pc_z[k] * sof_570[k];

        t_858[k] = pb_x[k] * sng0_858[k]
                   + f_8 * snf_573[k]
                   - f_6 * pc_x[k] * sng1_858[k];

        t_859[k] = f_11 * snf_472[k]
                   + f_3 * pc_y[k] * sof_572[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, pb_x, pc_x, sng0_860, snf_575, snf_576, \
                         snf_577, snf_578, sng1_860, sof_576, sof_577, \
                         sof_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = pb_x[k] * sng0_860[k]
                   + f_8 * snf_575[k]
                   - f_6 * pc_x[k] * sng1_860[k];

        t_861[k] = f_7 * snf_576[k]
                   + f_3 * pc_x[k] * sof_576[k];

        t_862[k] = f_7 * snf_577[k]
                   + f_3 * pc_x[k] * sof_577[k];

        t_863[k] = f_7 * snf_578[k]
                   + f_3 * pc_x[k] * sof_578[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, t_867, pb_x, pc_x, pc_z, sng0_865, sng0_867, \
                         snf_466, snf_579, sng1_865, sng1_867, sof_576, \
                         sof_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = f_7 * snf_579[k]
                   + f_3 * pc_x[k] * sof_579[k];

        t_865[k] = pb_x[k] * sng0_865[k]
                   - f_6 * pc_x[k] * sng1_865[k];

        t_866[k] = f_8 * snf_466[k]
                   + f_3 * pc_z[k] * sof_576[k];

        t_867[k] = pb_x[k] * sng0_867[k]
                   - f_6 * pc_x[k] * sng1_867[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, t_871, pb_x, pc_x, pc_y, sng0_869, sng0_870, \
                         snf_479, snf_480, snf_580, sng1_869, sng1_870, sof_579, \
                         sof_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = f_11 * snf_479[k]
                   + f_3 * pc_y[k] * sof_579[k];

        t_869[k] = pb_x[k] * sng0_869[k]
                   - f_6 * pc_x[k] * sng1_869[k];

        t_870[k] = pb_x[k] * sng0_870[k]
                   + f_14 * snf_580[k]
                   - f_6 * pc_x[k] * sng1_870[k];

        t_871[k] = f_13 * snf_480[k]
                   + f_3 * pc_y[k] * sof_580[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pb_x, pc_x, pc_y, pc_z, sng0_873, snf_470, \
                         snf_482, snf_583, sng1_873, sof_580, sof_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_12 * snf_470[k]
                   + f_3 * pc_z[k] * sof_580[k];

        t_873[k] = pb_x[k] * sng0_873[k]
                   + f_8 * snf_583[k]
                   - f_6 * pc_x[k] * sng1_873[k];

        t_874[k] = f_13 * snf_482[k]
                   + f_3 * pc_y[k] * sof_582[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, pb_x, pc_x, sng0_875, snf_585, snf_586, \
                         snf_587, snf_588, sng1_875, sof_586, sof_587, \
                         sof_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = pb_x[k] * sng0_875[k]
                   + f_8 * snf_585[k]
                   - f_6 * pc_x[k] * sng1_875[k];

        t_876[k] = f_7 * snf_586[k]
                   + f_3 * pc_x[k] * sof_586[k];

        t_877[k] = f_7 * snf_587[k]
                   + f_3 * pc_x[k] * sof_587[k];

        t_878[k] = f_7 * snf_588[k]
                   + f_3 * pc_x[k] * sof_588[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, t_882, pb_x, pc_x, pc_z, sng0_880, sng0_882, \
                         snf_476, snf_589, sng1_880, sng1_882, sof_586, \
                         sof_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_7 * snf_589[k]
                   + f_3 * pc_x[k] * sof_589[k];

        t_880[k] = pb_x[k] * sng0_880[k]
                   - f_6 * pc_x[k] * sng1_880[k];

        t_881[k] = f_12 * snf_476[k]
                   + f_3 * pc_z[k] * sof_586[k];

        t_882[k] = pb_x[k] * sng0_882[k]
                   - f_6 * pc_x[k] * sng1_882[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pb_x, pc_x, pc_y, sng0_884, sng0_885, \
                         snf_489, snf_490, snf_590, sng1_884, sng1_885, sof_589, \
                         sof_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_13 * snf_489[k]
                   + f_3 * pc_y[k] * sof_589[k];

        t_884[k] = pb_x[k] * sng0_884[k]
                   - f_6 * pc_x[k] * sng1_884[k];

        t_885[k] = pb_x[k] * sng0_885[k]
                   + f_14 * snf_590[k]
                   - f_6 * pc_x[k] * sng1_885[k];

        t_886[k] = f_15 * snf_490[k]
                   + f_3 * pc_y[k] * sof_590[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, pb_x, pc_x, pc_y, pc_z, sng0_888, snf_480, \
                         snf_492, snf_593, sng1_888, sof_590, sof_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_14 * snf_480[k]
                   + f_3 * pc_z[k] * sof_590[k];

        t_888[k] = pb_x[k] * sng0_888[k]
                   + f_8 * snf_593[k]
                   - f_6 * pc_x[k] * sng1_888[k];

        t_889[k] = f_15 * snf_492[k]
                   + f_3 * pc_y[k] * sof_592[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pb_x, pc_x, sng0_890, snf_595, snf_596, \
                         snf_597, snf_598, sng1_890, sof_596, sof_597, \
                         sof_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = pb_x[k] * sng0_890[k]
                   + f_8 * snf_595[k]
                   - f_6 * pc_x[k] * sng1_890[k];

        t_891[k] = f_7 * snf_596[k]
                   + f_3 * pc_x[k] * sof_596[k];

        t_892[k] = f_7 * snf_597[k]
                   + f_3 * pc_x[k] * sof_597[k];

        t_893[k] = f_7 * snf_598[k]
                   + f_3 * pc_x[k] * sof_598[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pb_x, pc_x, pc_z, sng0_895, sng0_897, \
                         snf_486, snf_599, sng1_895, sng1_897, sof_596, \
                         sof_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_7 * snf_599[k]
                   + f_3 * pc_x[k] * sof_599[k];

        t_895[k] = pb_x[k] * sng0_895[k]
                   - f_6 * pc_x[k] * sng1_895[k];

        t_896[k] = f_14 * snf_486[k]
                   + f_3 * pc_z[k] * sof_596[k];

        t_897[k] = pb_x[k] * sng0_897[k]
                   - f_6 * pc_x[k] * sng1_897[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, pb_x, pc_x, pc_y, sng0_899, sng0_900, \
                         snf_499, snf_500, snf_600, sng1_899, sng1_900, sof_599, \
                         sof_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * snf_499[k]
                   + f_3 * pc_y[k] * sof_599[k];

        t_899[k] = pb_x[k] * sng0_899[k]
                   - f_6 * pc_x[k] * sng1_899[k];

        t_900[k] = pb_x[k] * sng0_900[k]
                   + f_14 * snf_600[k]
                   - f_6 * pc_x[k] * sng1_900[k];

        t_901[k] = f_16 * snf_500[k]
                   + f_3 * pc_y[k] * sof_600[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pb_x, pc_x, pc_y, pc_z, sng0_903, snf_490, \
                         snf_502, snf_603, sng1_903, sof_600, sof_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_16 * snf_490[k]
                   + f_3 * pc_z[k] * sof_600[k];

        t_903[k] = pb_x[k] * sng0_903[k]
                   + f_8 * snf_603[k]
                   - f_6 * pc_x[k] * sng1_903[k];

        t_904[k] = f_16 * snf_502[k]
                   + f_3 * pc_y[k] * sof_602[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pb_x, pc_x, sng0_905, snf_605, snf_606, \
                         snf_607, snf_608, sng1_905, sof_606, sof_607, \
                         sof_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = pb_x[k] * sng0_905[k]
                   + f_8 * snf_605[k]
                   - f_6 * pc_x[k] * sng1_905[k];

        t_906[k] = f_7 * snf_606[k]
                   + f_3 * pc_x[k] * sof_606[k];

        t_907[k] = f_7 * snf_607[k]
                   + f_3 * pc_x[k] * sof_607[k];

        t_908[k] = f_7 * snf_608[k]
                   + f_3 * pc_x[k] * sof_608[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, pb_x, pc_x, pc_z, sng0_910, sng0_912, \
                         snf_496, snf_609, sng1_910, sng1_912, sof_606, \
                         sof_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_7 * snf_609[k]
                   + f_3 * pc_x[k] * sof_609[k];

        t_910[k] = pb_x[k] * sng0_910[k]
                   - f_6 * pc_x[k] * sng1_910[k];

        t_911[k] = f_16 * snf_496[k]
                   + f_3 * pc_z[k] * sof_606[k];

        t_912[k] = pb_x[k] * sng0_912[k]
                   - f_6 * pc_x[k] * sng1_912[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, pb_x, pc_x, pc_y, sng0_914, sng0_915, \
                         snf_509, snf_510, snf_610, sng1_914, sng1_915, sof_609, \
                         sof_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = f_16 * snf_509[k]
                   + f_3 * pc_y[k] * sof_609[k];

        t_914[k] = pb_x[k] * sng0_914[k]
                   - f_6 * pc_x[k] * sng1_914[k];

        t_915[k] = pb_x[k] * sng0_915[k]
                   + f_14 * snf_610[k]
                   - f_6 * pc_x[k] * sng1_915[k];

        t_916[k] = f_14 * snf_510[k]
                   + f_3 * pc_y[k] * sof_610[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pb_x, pc_x, pc_y, pc_z, sng0_918, snf_500, \
                         snf_512, snf_613, sng1_918, sof_610, sof_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_15 * snf_500[k]
                   + f_3 * pc_z[k] * sof_610[k];

        t_918[k] = pb_x[k] * sng0_918[k]
                   + f_8 * snf_613[k]
                   - f_6 * pc_x[k] * sng1_918[k];

        t_919[k] = f_14 * snf_512[k]
                   + f_3 * pc_y[k] * sof_612[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pb_x, pc_x, sng0_920, snf_615, snf_616, \
                         snf_617, snf_618, sng1_920, sof_616, sof_617, \
                         sof_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = pb_x[k] * sng0_920[k]
                   + f_8 * snf_615[k]
                   - f_6 * pc_x[k] * sng1_920[k];

        t_921[k] = f_7 * snf_616[k]
                   + f_3 * pc_x[k] * sof_616[k];

        t_922[k] = f_7 * snf_617[k]
                   + f_3 * pc_x[k] * sof_617[k];

        t_923[k] = f_7 * snf_618[k]
                   + f_3 * pc_x[k] * sof_618[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pb_x, pc_x, pc_z, sng0_925, sng0_927, \
                         snf_506, snf_619, sng1_925, sng1_927, sof_616, \
                         sof_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_7 * snf_619[k]
                   + f_3 * pc_x[k] * sof_619[k];

        t_925[k] = pb_x[k] * sng0_925[k]
                   - f_6 * pc_x[k] * sng1_925[k];

        t_926[k] = f_15 * snf_506[k]
                   + f_3 * pc_z[k] * sof_616[k];

        t_927[k] = pb_x[k] * sng0_927[k]
                   - f_6 * pc_x[k] * sng1_927[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, t_931, pb_x, pc_x, pc_y, sng0_929, sng0_930, \
                         snf_519, snf_520, snf_620, sng1_929, sng1_930, sof_619, \
                         sof_620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_14 * snf_519[k]
                   + f_3 * pc_y[k] * sof_619[k];

        t_929[k] = pb_x[k] * sng0_929[k]
                   - f_6 * pc_x[k] * sng1_929[k];

        t_930[k] = pb_x[k] * sng0_930[k]
                   + f_14 * snf_620[k]
                   - f_6 * pc_x[k] * sng1_930[k];

        t_931[k] = f_12 * snf_520[k]
                   + f_3 * pc_y[k] * sof_620[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pb_x, pc_x, pc_y, pc_z, sng0_933, snf_510, \
                         snf_522, snf_623, sng1_933, sof_620, sof_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_13 * snf_510[k]
                   + f_3 * pc_z[k] * sof_620[k];

        t_933[k] = pb_x[k] * sng0_933[k]
                   + f_8 * snf_623[k]
                   - f_6 * pc_x[k] * sng1_933[k];

        t_934[k] = f_12 * snf_522[k]
                   + f_3 * pc_y[k] * sof_622[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, pb_x, pc_x, sng0_935, snf_625, snf_626, \
                         snf_627, snf_628, sng1_935, sof_626, sof_627, \
                         sof_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = pb_x[k] * sng0_935[k]
                   + f_8 * snf_625[k]
                   - f_6 * pc_x[k] * sng1_935[k];

        t_936[k] = f_7 * snf_626[k]
                   + f_3 * pc_x[k] * sof_626[k];

        t_937[k] = f_7 * snf_627[k]
                   + f_3 * pc_x[k] * sof_627[k];

        t_938[k] = f_7 * snf_628[k]
                   + f_3 * pc_x[k] * sof_628[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pb_x, pc_x, pc_z, sng0_940, sng0_942, \
                         snf_516, snf_629, sng1_940, sng1_942, sof_626, \
                         sof_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_7 * snf_629[k]
                   + f_3 * pc_x[k] * sof_629[k];

        t_940[k] = pb_x[k] * sng0_940[k]
                   - f_6 * pc_x[k] * sng1_940[k];

        t_941[k] = f_13 * snf_516[k]
                   + f_3 * pc_z[k] * sof_626[k];

        t_942[k] = pb_x[k] * sng0_942[k]
                   - f_6 * pc_x[k] * sng1_942[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, pb_x, pc_x, pc_y, sng0_944, sng0_945, \
                         snf_529, snf_530, snf_630, sng1_944, sng1_945, sof_629, \
                         sof_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_12 * snf_529[k]
                   + f_3 * pc_y[k] * sof_629[k];

        t_944[k] = pb_x[k] * sng0_944[k]
                   - f_6 * pc_x[k] * sng1_944[k];

        t_945[k] = pb_x[k] * sng0_945[k]
                   + f_14 * snf_630[k]
                   - f_6 * pc_x[k] * sng1_945[k];

        t_946[k] = f_8 * snf_530[k]
                   + f_3 * pc_y[k] * sof_630[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.0 / q;
    const auto f_10 = 4.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_810 = buffer.data(sng0 + 810);
    const auto *sng0_815 = buffer.data(sng0 + 815);
    const auto *sng0_825 = buffer.data(sng0 + 825);
    const auto *sng0_828 = buffer.data(sng0 + 828);
    const auto *sng0_835 = buffer.data(sng0 + 835);
    const auto *sng0_837 = buffer.data(sng0 + 837);
    const auto *sng0_948 = buffer.data(sng0 + 948);
    const auto *sng0_950 = buffer.data(sng0 + 950);
    const auto *sng0_955 = buffer.data(sng0 + 955);
    const auto *sng0_957 = buffer.data(sng0 + 957);
    const auto *sng0_959 = buffer.data(sng0 + 959);
    const auto *sng0_963 = buffer.data(sng0 + 963);
    const auto *sng0_970 = buffer.data(sng0 + 970);
    const auto *sng0_972 = buffer.data(sng0 + 972);
    const auto *sng0_974 = buffer.data(sng0 + 974);
    const auto *sng0_975 = buffer.data(sng0 + 975);
    const auto *sng0_978 = buffer.data(sng0 + 978);
    const auto *sng0_980 = buffer.data(sng0 + 980);
    const auto *sng0_985 = buffer.data(sng0 + 985);
    const auto *sng0_987 = buffer.data(sng0 + 987);
    const auto *sng0_989 = buffer.data(sng0 + 989);

    const auto *snf_520 = buffer.data(snf + 520);
    const auto *snf_526 = buffer.data(snf + 526);
    const auto *snf_530 = buffer.data(snf + 530);
    const auto *snf_532 = buffer.data(snf + 532);
    const auto *snf_536 = buffer.data(snf + 536);
    const auto *snf_539 = buffer.data(snf + 539);
    const auto *snf_540 = buffer.data(snf + 540);
    const auto *snf_542 = buffer.data(snf + 542);
    const auto *snf_546 = buffer.data(snf + 546);
    const auto *snf_549 = buffer.data(snf + 549);
    const auto *snf_550 = buffer.data(snf + 550);
    const auto *snf_552 = buffer.data(snf + 552);
    const auto *snf_556 = buffer.data(snf + 556);
    const auto *snf_557 = buffer.data(snf + 557);
    const auto *snf_558 = buffer.data(snf + 558);
    const auto *snf_559 = buffer.data(snf + 559);
    const auto *snf_560 = buffer.data(snf + 560);
    const auto *snf_562 = buffer.data(snf + 562);
    const auto *snf_566 = buffer.data(snf + 566);
    const auto *snf_569 = buffer.data(snf + 569);
    const auto *snf_570 = buffer.data(snf + 570);
    const auto *snf_572 = buffer.data(snf + 572);
    const auto *snf_576 = buffer.data(snf + 576);
    const auto *snf_578 = buffer.data(snf + 578);
    const auto *snf_579 = buffer.data(snf + 579);
    const auto *snf_580 = buffer.data(snf + 580);
    const auto *snf_582 = buffer.data(snf + 582);
    const auto *snf_586 = buffer.data(snf + 586);
    const auto *snf_588 = buffer.data(snf + 588);
    const auto *snf_589 = buffer.data(snf + 589);
    const auto *snf_590 = buffer.data(snf + 590);
    const auto *snf_592 = buffer.data(snf + 592);
    const auto *snf_596 = buffer.data(snf + 596);
    const auto *snf_598 = buffer.data(snf + 598);
    const auto *snf_599 = buffer.data(snf + 599);
    const auto *snf_600 = buffer.data(snf + 600);
    const auto *snf_602 = buffer.data(snf + 602);
    const auto *snf_606 = buffer.data(snf + 606);
    const auto *snf_633 = buffer.data(snf + 633);
    const auto *snf_635 = buffer.data(snf + 635);
    const auto *snf_636 = buffer.data(snf + 636);
    const auto *snf_637 = buffer.data(snf + 637);
    const auto *snf_638 = buffer.data(snf + 638);
    const auto *snf_639 = buffer.data(snf + 639);
    const auto *snf_643 = buffer.data(snf + 643);
    const auto *snf_646 = buffer.data(snf + 646);
    const auto *snf_647 = buffer.data(snf + 647);
    const auto *snf_648 = buffer.data(snf + 648);
    const auto *snf_649 = buffer.data(snf + 649);
    const auto *snf_650 = buffer.data(snf + 650);
    const auto *snf_653 = buffer.data(snf + 653);
    const auto *snf_655 = buffer.data(snf + 655);
    const auto *snf_656 = buffer.data(snf + 656);
    const auto *snf_657 = buffer.data(snf + 657);
    const auto *snf_658 = buffer.data(snf + 658);
    const auto *snf_659 = buffer.data(snf + 659);

    const auto *sng1_810 = buffer.data(sng1 + 810);
    const auto *sng1_815 = buffer.data(sng1 + 815);
    const auto *sng1_825 = buffer.data(sng1 + 825);
    const auto *sng1_828 = buffer.data(sng1 + 828);
    const auto *sng1_835 = buffer.data(sng1 + 835);
    const auto *sng1_837 = buffer.data(sng1 + 837);
    const auto *sng1_948 = buffer.data(sng1 + 948);
    const auto *sng1_950 = buffer.data(sng1 + 950);
    const auto *sng1_955 = buffer.data(sng1 + 955);
    const auto *sng1_957 = buffer.data(sng1 + 957);
    const auto *sng1_959 = buffer.data(sng1 + 959);
    const auto *sng1_963 = buffer.data(sng1 + 963);
    const auto *sng1_970 = buffer.data(sng1 + 970);
    const auto *sng1_972 = buffer.data(sng1 + 972);
    const auto *sng1_974 = buffer.data(sng1 + 974);
    const auto *sng1_975 = buffer.data(sng1 + 975);
    const auto *sng1_978 = buffer.data(sng1 + 978);
    const auto *sng1_980 = buffer.data(sng1 + 980);
    const auto *sng1_985 = buffer.data(sng1 + 985);
    const auto *sng1_987 = buffer.data(sng1 + 987);
    const auto *sng1_989 = buffer.data(sng1 + 989);

    const auto *sod0_396 = buffer.data(sod0 + 396);
    const auto *sod0_399 = buffer.data(sod0 + 399);
    const auto *sod0_401 = buffer.data(sod0 + 401);
    const auto *sod0_407 = buffer.data(sod0 + 407);
    const auto *sod0_408 = buffer.data(sod0 + 408);
    const auto *sod0_411 = buffer.data(sod0 + 411);
    const auto *sod0_413 = buffer.data(sod0 + 413);
    const auto *sod0_414 = buffer.data(sod0 + 414);
    const auto *sod0_417 = buffer.data(sod0 + 417);
    const auto *sod0_419 = buffer.data(sod0 + 419);
    const auto *sod0_420 = buffer.data(sod0 + 420);
    const auto *sod0_423 = buffer.data(sod0 + 423);
    const auto *sod0_425 = buffer.data(sod0 + 425);
    const auto *sod0_426 = buffer.data(sod0 + 426);
    const auto *sod0_429 = buffer.data(sod0 + 429);
    const auto *sod0_431 = buffer.data(sod0 + 431);

    const auto *sod1_396 = buffer.data(sod1 + 396);
    const auto *sod1_399 = buffer.data(sod1 + 399);
    const auto *sod1_401 = buffer.data(sod1 + 401);
    const auto *sod1_407 = buffer.data(sod1 + 407);
    const auto *sod1_408 = buffer.data(sod1 + 408);
    const auto *sod1_411 = buffer.data(sod1 + 411);
    const auto *sod1_413 = buffer.data(sod1 + 413);
    const auto *sod1_414 = buffer.data(sod1 + 414);
    const auto *sod1_417 = buffer.data(sod1 + 417);
    const auto *sod1_419 = buffer.data(sod1 + 419);
    const auto *sod1_420 = buffer.data(sod1 + 420);
    const auto *sod1_423 = buffer.data(sod1 + 423);
    const auto *sod1_425 = buffer.data(sod1 + 425);
    const auto *sod1_426 = buffer.data(sod1 + 426);
    const auto *sod1_429 = buffer.data(sod1 + 429);
    const auto *sod1_431 = buffer.data(sod1 + 431);

    const auto *sof_630 = buffer.data(sof + 630);
    const auto *sof_632 = buffer.data(sof + 632);
    const auto *sof_636 = buffer.data(sof + 636);
    const auto *sof_637 = buffer.data(sof + 637);
    const auto *sof_638 = buffer.data(sof + 638);
    const auto *sof_639 = buffer.data(sof + 639);
    const auto *sof_640 = buffer.data(sof + 640);
    const auto *sof_642 = buffer.data(sof + 642);
    const auto *sof_646 = buffer.data(sof + 646);
    const auto *sof_647 = buffer.data(sof + 647);
    const auto *sof_648 = buffer.data(sof + 648);
    const auto *sof_649 = buffer.data(sof + 649);
    const auto *sof_650 = buffer.data(sof + 650);
    const auto *sof_652 = buffer.data(sof + 652);
    const auto *sof_656 = buffer.data(sof + 656);
    const auto *sof_657 = buffer.data(sof + 657);
    const auto *sof_658 = buffer.data(sof + 658);
    const auto *sof_659 = buffer.data(sof + 659);
    const auto *sof_660 = buffer.data(sof + 660);
    const auto *sof_662 = buffer.data(sof + 662);
    const auto *sof_663 = buffer.data(sof + 663);
    const auto *sof_665 = buffer.data(sof + 665);
    const auto *sof_666 = buffer.data(sof + 666);
    const auto *sof_667 = buffer.data(sof + 667);
    const auto *sof_668 = buffer.data(sof + 668);
    const auto *sof_669 = buffer.data(sof + 669);
    const auto *sof_670 = buffer.data(sof + 670);
    const auto *sof_672 = buffer.data(sof + 672);
    const auto *sof_675 = buffer.data(sof + 675);
    const auto *sof_676 = buffer.data(sof + 676);
    const auto *sof_677 = buffer.data(sof + 677);
    const auto *sof_678 = buffer.data(sof + 678);
    const auto *sof_679 = buffer.data(sof + 679);
    const auto *sof_680 = buffer.data(sof + 680);
    const auto *sof_682 = buffer.data(sof + 682);
    const auto *sof_683 = buffer.data(sof + 683);
    const auto *sof_685 = buffer.data(sof + 685);
    const auto *sof_686 = buffer.data(sof + 686);
    const auto *sof_687 = buffer.data(sof + 687);
    const auto *sof_688 = buffer.data(sof + 688);
    const auto *sof_689 = buffer.data(sof + 689);
    const auto *sof_690 = buffer.data(sof + 690);
    const auto *sof_692 = buffer.data(sof + 692);
    const auto *sof_693 = buffer.data(sof + 693);
    const auto *sof_695 = buffer.data(sof + 695);
    const auto *sof_696 = buffer.data(sof + 696);
    const auto *sof_697 = buffer.data(sof + 697);
    const auto *sof_698 = buffer.data(sof + 698);
    const auto *sof_699 = buffer.data(sof + 699);
    const auto *sof_700 = buffer.data(sof + 700);
    const auto *sof_702 = buffer.data(sof + 702);
    const auto *sof_703 = buffer.data(sof + 703);
    const auto *sof_705 = buffer.data(sof + 705);
    const auto *sof_706 = buffer.data(sof + 706);
    const auto *sof_707 = buffer.data(sof + 707);
    const auto *sof_708 = buffer.data(sof + 708);
    const auto *sof_709 = buffer.data(sof + 709);
    const auto *sof_710 = buffer.data(sof + 710);
    const auto *sof_712 = buffer.data(sof + 712);
    const auto *sof_713 = buffer.data(sof + 713);
    const auto *sof_715 = buffer.data(sof + 715);
    const auto *sof_716 = buffer.data(sof + 716);
    const auto *sof_717 = buffer.data(sof + 717);
    const auto *sof_718 = buffer.data(sof + 718);
    const auto *sof_719 = buffer.data(sof + 719);

#pragma omp simd aligned(t_947, t_948, t_949, pb_x, pc_x, pc_y, pc_z, sng0_948, snf_520, \
                         snf_532, snf_633, sng1_948, sof_630, sof_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_11 * snf_520[k]
                   + f_3 * pc_z[k] * sof_630[k];

        t_948[k] = pb_x[k] * sng0_948[k]
                   + f_8 * snf_633[k]
                   - f_6 * pc_x[k] * sng1_948[k];

        t_949[k] = f_8 * snf_532[k]
                   + f_3 * pc_y[k] * sof_632[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, pb_x, pc_x, sng0_950, snf_635, snf_636, \
                         snf_637, snf_638, sng1_950, sof_636, sof_637, \
                         sof_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = pb_x[k] * sng0_950[k]
                   + f_8 * snf_635[k]
                   - f_6 * pc_x[k] * sng1_950[k];

        t_951[k] = f_7 * snf_636[k]
                   + f_3 * pc_x[k] * sof_636[k];

        t_952[k] = f_7 * snf_637[k]
                   + f_3 * pc_x[k] * sof_637[k];

        t_953[k] = f_7 * snf_638[k]
                   + f_3 * pc_x[k] * sof_638[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, t_957, pb_x, pc_x, pc_z, sng0_955, sng0_957, \
                         snf_526, snf_639, sng1_955, sng1_957, sof_636, \
                         sof_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_7 * snf_639[k]
                   + f_3 * pc_x[k] * sof_639[k];

        t_955[k] = pb_x[k] * sng0_955[k]
                   - f_6 * pc_x[k] * sng1_955[k];

        t_956[k] = f_11 * snf_526[k]
                   + f_3 * pc_z[k] * sof_636[k];

        t_957[k] = pb_x[k] * sng0_957[k]
                   - f_6 * pc_x[k] * sng1_957[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, t_961, pb_x, pb_y, pc_x, pc_y, sng0_810, \
                         sng0_959, snf_539, snf_540, sng1_810, sng1_959, sof_639, \
                         sof_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = f_8 * snf_539[k]
                   + f_3 * pc_y[k] * sof_639[k];

        t_959[k] = pb_x[k] * sng0_959[k]
                   - f_6 * pc_x[k] * sng1_959[k];

        t_960[k] = pb_y[k] * sng0_810[k]
                   - f_6 * pc_y[k] * sng1_810[k];

        t_961[k] = f_7 * snf_540[k]
                   + f_3 * pc_y[k] * sof_640[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, pb_x, pc_x, pc_y, pc_z, sng0_963, snf_530, \
                         snf_542, snf_643, sng1_963, sof_640, sof_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_10 * snf_530[k]
                   + f_3 * pc_z[k] * sof_640[k];

        t_963[k] = pb_x[k] * sng0_963[k]
                   + f_8 * snf_643[k]
                   - f_6 * pc_x[k] * sng1_963[k];

        t_964[k] = f_7 * snf_542[k]
                   + f_3 * pc_y[k] * sof_642[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, pb_y, pc_x, pc_y, sng0_815, snf_646, \
                         snf_647, snf_648, sng1_815, sof_646, sof_647, \
                         sof_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = pb_y[k] * sng0_815[k]
                   - f_6 * pc_y[k] * sng1_815[k];

        t_966[k] = f_7 * snf_646[k]
                   + f_3 * pc_x[k] * sof_646[k];

        t_967[k] = f_7 * snf_647[k]
                   + f_3 * pc_x[k] * sof_647[k];

        t_968[k] = f_7 * snf_648[k]
                   + f_3 * pc_x[k] * sof_648[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, pb_x, pc_x, pc_z, sng0_970, sng0_972, \
                         snf_536, snf_649, sng1_970, sng1_972, sof_646, \
                         sof_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_7 * snf_649[k]
                   + f_3 * pc_x[k] * sof_649[k];

        t_970[k] = pb_x[k] * sng0_970[k]
                   - f_6 * pc_x[k] * sng1_970[k];

        t_971[k] = f_10 * snf_536[k]
                   + f_3 * pc_z[k] * sof_646[k];

        t_972[k] = pb_x[k] * sng0_972[k]
                   - f_6 * pc_x[k] * sng1_972[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, pb_x, pc_x, pc_y, sng0_974, sng0_975, \
                         snf_549, snf_650, sng1_974, sng1_975, sof_649, \
                         sof_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_7 * snf_549[k]
                   + f_3 * pc_y[k] * sof_649[k];

        t_974[k] = pb_x[k] * sng0_974[k]
                   - f_6 * pc_x[k] * sng1_974[k];

        t_975[k] = pb_x[k] * sng0_975[k]
                   + f_14 * snf_650[k]
                   - f_6 * pc_x[k] * sng1_975[k];

        t_976[k] = f_3 * pc_y[k] * sof_650[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, pb_x, pc_x, pc_y, pc_z, sng0_978, snf_540, \
                         snf_653, sng1_978, sof_650, sof_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_9 * snf_540[k]
                   + f_3 * pc_z[k] * sof_650[k];

        t_978[k] = pb_x[k] * sng0_978[k]
                   + f_8 * snf_653[k]
                   - f_6 * pc_x[k] * sng1_978[k];

        t_979[k] = f_3 * pc_y[k] * sof_652[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, pb_x, pc_x, sng0_980, snf_655, snf_656, \
                         snf_657, snf_658, sng1_980, sof_656, sof_657, \
                         sof_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = pb_x[k] * sng0_980[k]
                   + f_8 * snf_655[k]
                   - f_6 * pc_x[k] * sng1_980[k];

        t_981[k] = f_7 * snf_656[k]
                   + f_3 * pc_x[k] * sof_656[k];

        t_982[k] = f_7 * snf_657[k]
                   + f_3 * pc_x[k] * sof_657[k];

        t_983[k] = f_7 * snf_658[k]
                   + f_3 * pc_x[k] * sof_658[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, t_987, pb_x, pc_x, pc_z, sng0_985, sng0_987, \
                         snf_546, snf_659, sng1_985, sng1_987, sof_656, \
                         sof_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_7 * snf_659[k]
                   + f_3 * pc_x[k] * sof_659[k];

        t_985[k] = pb_x[k] * sng0_985[k]
                   - f_6 * pc_x[k] * sng1_985[k];

        t_986[k] = f_9 * snf_546[k]
                   + f_3 * pc_z[k] * sof_656[k];

        t_987[k] = pb_x[k] * sng0_987[k]
                   - f_6 * pc_x[k] * sng1_987[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, t_991, t_992, pb_x, pc_x, pc_y, pc_z, sng0_989, \
                         snf_550, sng1_989, sod0_396, sod1_396, sof_659, \
                         sof_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_3 * pc_y[k] * sof_659[k];

        t_989[k] = pb_x[k] * sng0_989[k]
                   - f_6 * pc_x[k] * sng1_989[k];

        t_990[k] = f_1 * sod0_396[k]
                   - f_2 * sod1_396[k]
                   + f_3 * pc_x[k] * sof_660[k];

        t_991[k] = f_0 * snf_550[k]
                   + f_3 * pc_y[k] * sof_660[k];

        t_992[k] = f_3 * pc_z[k] * sof_660[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, pc_x, pc_y, snf_552, sod0_399, sod0_401, \
                         sod1_399, sod1_401, sof_662, sof_663, sof_665, \
                         sof_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = f_4 * sod0_399[k]
                   - f_5 * sod1_399[k]
                   + f_3 * pc_x[k] * sof_663[k];

        t_994[k] = f_0 * snf_552[k]
                   + f_3 * pc_y[k] * sof_662[k];

        t_995[k] = f_4 * sod0_401[k]
                   - f_5 * sod1_401[k]
                   + f_3 * pc_x[k] * sof_665[k];

        t_996[k] = f_3 * pc_x[k] * sof_666[k];
    }

#pragma omp simd aligned(t_997, t_998, t_999, t_1000, t_1001, pc_x, pc_y, pc_z, snf_556, \
                         sod0_399, sod1_399, sof_666, sof_667, sof_668, \
                         sof_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = f_3 * pc_x[k] * sof_667[k];

        t_998[k] = f_3 * pc_x[k] * sof_668[k];

        t_999[k] = f_3 * pc_x[k] * sof_669[k];

        t_1000[k] = f_0 * snf_556[k]
                    + f_1 * sod0_399[k]
                    - f_2 * sod1_399[k]
                    + f_3 * pc_y[k] * sof_666[k];

        t_1001[k] = f_3 * pc_z[k] * sof_666[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, t_1005, pb_z, pc_y, pc_z, sng0_825, snf_558, \
                         snf_559, sng1_825, sod0_401, sod1_401, sof_668, \
                         sof_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_0 * snf_558[k]
                    + f_4 * sod0_401[k]
                    - f_5 * sod1_401[k]
                    + f_3 * pc_y[k] * sof_668[k];

        t_1003[k] = f_0 * snf_559[k]
                    + f_3 * pc_y[k] * sof_669[k];

        t_1004[k] = f_1 * sod0_401[k]
                    - f_2 * sod1_401[k]
                    + f_3 * pc_z[k] * sof_669[k];

        t_1005[k] = pb_z[k] * sng0_825[k]
                    - f_6 * pc_z[k] * sng1_825[k];
    }

#pragma omp simd aligned(t_1006, t_1007, t_1008, t_1009, pb_z, pc_y, pc_z, sng0_828, snf_550, \
                         snf_560, snf_562, sng1_828, sof_670, sof_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1006[k] = f_9 * snf_560[k]
                    + f_3 * pc_y[k] * sof_670[k];

        t_1007[k] = f_7 * snf_550[k]
                    + f_3 * pc_z[k] * sof_670[k];

        t_1008[k] = pb_z[k] * sng0_828[k]
                    - f_6 * pc_z[k] * sng1_828[k];

        t_1009[k] = f_9 * snf_562[k]
                    + f_3 * pc_y[k] * sof_672[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, t_1013, t_1014, pc_x, sod0_407, sod1_407, \
                         sof_675, sof_676, sof_677, sof_678, sof_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_4 * sod0_407[k]
                    - f_5 * sod1_407[k]
                    + f_3 * pc_x[k] * sof_675[k];

        t_1011[k] = f_3 * pc_x[k] * sof_676[k];

        t_1012[k] = f_3 * pc_x[k] * sof_677[k];

        t_1013[k] = f_3 * pc_x[k] * sof_678[k];

        t_1014[k] = f_3 * pc_x[k] * sof_679[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, t_1018, pb_z, pc_y, pc_z, sng0_835, sng0_837, \
                         snf_556, snf_557, snf_569, sng1_835, sng1_837, sof_676, \
                         sof_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = pb_z[k] * sng0_835[k]
                    - f_6 * pc_z[k] * sng1_835[k];

        t_1016[k] = f_7 * snf_556[k]
                    + f_3 * pc_z[k] * sof_676[k];

        t_1017[k] = pb_z[k] * sng0_837[k]
                    + f_8 * snf_557[k]
                    - f_6 * pc_z[k] * sng1_837[k];

        t_1018[k] = f_9 * snf_569[k]
                    + f_3 * pc_y[k] * sof_679[k];
    }

#pragma omp simd aligned(t_1019, t_1020, t_1021, t_1022, pc_x, pc_y, pc_z, snf_559, snf_560, \
                         snf_570, sod0_407, sod0_408, sod1_407, sod1_408, sof_679, \
                         sof_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1019[k] = f_7 * snf_559[k]
                    + f_1 * sod0_407[k]
                    - f_2 * sod1_407[k]
                    + f_3 * pc_z[k] * sof_679[k];

        t_1020[k] = f_1 * sod0_408[k]
                    - f_2 * sod1_408[k]
                    + f_3 * pc_x[k] * sof_680[k];

        t_1021[k] = f_10 * snf_570[k]
                    + f_3 * pc_y[k] * sof_680[k];

        t_1022[k] = f_8 * snf_560[k]
                    + f_3 * pc_z[k] * sof_680[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, t_1026, pc_x, pc_y, snf_572, sod0_411, \
                         sod0_413, sod1_411, sod1_413, sof_682, sof_683, sof_685, \
                         sof_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = f_4 * sod0_411[k]
                    - f_5 * sod1_411[k]
                    + f_3 * pc_x[k] * sof_683[k];

        t_1024[k] = f_10 * snf_572[k]
                    + f_3 * pc_y[k] * sof_682[k];

        t_1025[k] = f_4 * sod0_413[k]
                    - f_5 * sod1_413[k]
                    + f_3 * pc_x[k] * sof_685[k];

        t_1026[k] = f_3 * pc_x[k] * sof_686[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, t_1031, pc_x, pc_y, pc_z, snf_566, \
                         snf_576, sod0_411, sod1_411, sof_686, sof_687, sof_688, \
                         sof_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_3 * pc_x[k] * sof_687[k];

        t_1028[k] = f_3 * pc_x[k] * sof_688[k];

        t_1029[k] = f_3 * pc_x[k] * sof_689[k];

        t_1030[k] = f_10 * snf_576[k]
                    + f_1 * sod0_411[k]
                    - f_2 * sod1_411[k]
                    + f_3 * pc_y[k] * sof_686[k];

        t_1031[k] = f_8 * snf_566[k]
                    + f_3 * pc_z[k] * sof_686[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, pc_y, pc_z, snf_569, snf_578, snf_579, \
                         sod0_413, sod1_413, sof_688, sof_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = f_10 * snf_578[k]
                    + f_4 * sod0_413[k]
                    - f_5 * sod1_413[k]
                    + f_3 * pc_y[k] * sof_688[k];

        t_1033[k] = f_10 * snf_579[k]
                    + f_3 * pc_y[k] * sof_689[k];

        t_1034[k] = f_8 * snf_569[k]
                    + f_1 * sod0_413[k]
                    - f_2 * sod1_413[k]
                    + f_3 * pc_z[k] * sof_689[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, pc_x, pc_y, pc_z, snf_570, snf_580, \
                         sod0_414, sod0_417, sod1_414, sod1_417, sof_690, \
                         sof_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = f_1 * sod0_414[k]
                    - f_2 * sod1_414[k]
                    + f_3 * pc_x[k] * sof_690[k];

        t_1036[k] = f_11 * snf_580[k]
                    + f_3 * pc_y[k] * sof_690[k];

        t_1037[k] = f_12 * snf_570[k]
                    + f_3 * pc_z[k] * sof_690[k];

        t_1038[k] = f_4 * sod0_417[k]
                    - f_5 * sod1_417[k]
                    + f_3 * pc_x[k] * sof_693[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, t_1042, t_1043, pc_x, pc_y, snf_582, \
                         sod0_419, sod1_419, sof_692, sof_695, sof_696, sof_697, \
                         sof_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_11 * snf_582[k]
                    + f_3 * pc_y[k] * sof_692[k];

        t_1040[k] = f_4 * sod0_419[k]
                    - f_5 * sod1_419[k]
                    + f_3 * pc_x[k] * sof_695[k];

        t_1041[k] = f_3 * pc_x[k] * sof_696[k];

        t_1042[k] = f_3 * pc_x[k] * sof_697[k];

        t_1043[k] = f_3 * pc_x[k] * sof_698[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, pc_x, pc_y, pc_z, snf_576, snf_586, sod0_417, \
                         sod1_417, sof_696, sof_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_3 * pc_x[k] * sof_699[k];

        t_1045[k] = f_11 * snf_586[k]
                    + f_1 * sod0_417[k]
                    - f_2 * sod1_417[k]
                    + f_3 * pc_y[k] * sof_696[k];

        t_1046[k] = f_12 * snf_576[k]
                    + f_3 * pc_z[k] * sof_696[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pc_y, pc_z, snf_579, snf_588, snf_589, \
                         sod0_419, sod1_419, sof_698, sof_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = f_11 * snf_588[k]
                    + f_4 * sod0_419[k]
                    - f_5 * sod1_419[k]
                    + f_3 * pc_y[k] * sof_698[k];

        t_1048[k] = f_11 * snf_589[k]
                    + f_3 * pc_y[k] * sof_699[k];

        t_1049[k] = f_12 * snf_579[k]
                    + f_1 * sod0_419[k]
                    - f_2 * sod1_419[k]
                    + f_3 * pc_z[k] * sof_699[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pc_x, pc_y, pc_z, snf_580, snf_590, \
                         sod0_420, sod0_423, sod1_420, sod1_423, sof_700, \
                         sof_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_1 * sod0_420[k]
                    - f_2 * sod1_420[k]
                    + f_3 * pc_x[k] * sof_700[k];

        t_1051[k] = f_13 * snf_590[k]
                    + f_3 * pc_y[k] * sof_700[k];

        t_1052[k] = f_14 * snf_580[k]
                    + f_3 * pc_z[k] * sof_700[k];

        t_1053[k] = f_4 * sod0_423[k]
                    - f_5 * sod1_423[k]
                    + f_3 * pc_x[k] * sof_703[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, t_1058, pc_x, pc_y, snf_592, \
                         sod0_425, sod1_425, sof_702, sof_705, sof_706, sof_707, \
                         sof_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_13 * snf_592[k]
                    + f_3 * pc_y[k] * sof_702[k];

        t_1055[k] = f_4 * sod0_425[k]
                    - f_5 * sod1_425[k]
                    + f_3 * pc_x[k] * sof_705[k];

        t_1056[k] = f_3 * pc_x[k] * sof_706[k];

        t_1057[k] = f_3 * pc_x[k] * sof_707[k];

        t_1058[k] = f_3 * pc_x[k] * sof_708[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, pc_x, pc_y, pc_z, snf_586, snf_596, sod0_423, \
                         sod1_423, sof_706, sof_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = f_3 * pc_x[k] * sof_709[k];

        t_1060[k] = f_13 * snf_596[k]
                    + f_1 * sod0_423[k]
                    - f_2 * sod1_423[k]
                    + f_3 * pc_y[k] * sof_706[k];

        t_1061[k] = f_14 * snf_586[k]
                    + f_3 * pc_z[k] * sof_706[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, pc_y, pc_z, snf_589, snf_598, snf_599, \
                         sod0_425, sod1_425, sof_708, sof_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = f_13 * snf_598[k]
                    + f_4 * sod0_425[k]
                    - f_5 * sod1_425[k]
                    + f_3 * pc_y[k] * sof_708[k];

        t_1063[k] = f_13 * snf_599[k]
                    + f_3 * pc_y[k] * sof_709[k];

        t_1064[k] = f_14 * snf_589[k]
                    + f_1 * sod0_425[k]
                    - f_2 * sod1_425[k]
                    + f_3 * pc_z[k] * sof_709[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, pc_x, pc_y, pc_z, snf_590, snf_600, \
                         sod0_426, sod0_429, sod1_426, sod1_429, sof_710, \
                         sof_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_1 * sod0_426[k]
                    - f_2 * sod1_426[k]
                    + f_3 * pc_x[k] * sof_710[k];

        t_1066[k] = f_15 * snf_600[k]
                    + f_3 * pc_y[k] * sof_710[k];

        t_1067[k] = f_16 * snf_590[k]
                    + f_3 * pc_z[k] * sof_710[k];

        t_1068[k] = f_4 * sod0_429[k]
                    - f_5 * sod1_429[k]
                    + f_3 * pc_x[k] * sof_713[k];
    }

#pragma omp simd aligned(t_1069, t_1070, t_1071, t_1072, t_1073, pc_x, pc_y, snf_602, \
                         sod0_431, sod1_431, sof_712, sof_715, sof_716, sof_717, \
                         sof_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1069[k] = f_15 * snf_602[k]
                    + f_3 * pc_y[k] * sof_712[k];

        t_1070[k] = f_4 * sod0_431[k]
                    - f_5 * sod1_431[k]
                    + f_3 * pc_x[k] * sof_715[k];

        t_1071[k] = f_3 * pc_x[k] * sof_716[k];

        t_1072[k] = f_3 * pc_x[k] * sof_717[k];

        t_1073[k] = f_3 * pc_x[k] * sof_718[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, pc_x, pc_y, pc_z, snf_596, snf_606, sod0_429, \
                         sod1_429, sof_716, sof_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = f_3 * pc_x[k] * sof_719[k];

        t_1075[k] = f_15 * snf_606[k]
                    + f_1 * sod0_429[k]
                    - f_2 * sod1_429[k]
                    + f_3 * pc_y[k] * sof_716[k];

        t_1076[k] = f_16 * snf_596[k]
                    + f_3 * pc_z[k] * sof_716[k];
    }
}

static auto
compute_prim_sog_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sng0,
                                                          const size_t snf, const size_t sng1,
                                                          const size_t sod0, const size_t sod1,
                                                          const size_t sof, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.0 / q;
    const auto f_10 = 4.5 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sng0_975 = buffer.data(sng0 + 975);
    const auto *sng0_980 = buffer.data(sng0 + 980);
    const auto *sng0_985 = buffer.data(sng0 + 985);
    const auto *sng0_987 = buffer.data(sng0 + 987);
    const auto *sng0_989 = buffer.data(sng0 + 989);

    const auto *snf_599 = buffer.data(snf + 599);
    const auto *snf_600 = buffer.data(snf + 600);
    const auto *snf_606 = buffer.data(snf + 606);
    const auto *snf_608 = buffer.data(snf + 608);
    const auto *snf_609 = buffer.data(snf + 609);
    const auto *snf_610 = buffer.data(snf + 610);
    const auto *snf_612 = buffer.data(snf + 612);
    const auto *snf_616 = buffer.data(snf + 616);
    const auto *snf_618 = buffer.data(snf + 618);
    const auto *snf_619 = buffer.data(snf + 619);
    const auto *snf_620 = buffer.data(snf + 620);
    const auto *snf_622 = buffer.data(snf + 622);
    const auto *snf_626 = buffer.data(snf + 626);
    const auto *snf_628 = buffer.data(snf + 628);
    const auto *snf_629 = buffer.data(snf + 629);
    const auto *snf_630 = buffer.data(snf + 630);
    const auto *snf_632 = buffer.data(snf + 632);
    const auto *snf_636 = buffer.data(snf + 636);
    const auto *snf_638 = buffer.data(snf + 638);
    const auto *snf_639 = buffer.data(snf + 639);
    const auto *snf_640 = buffer.data(snf + 640);
    const auto *snf_642 = buffer.data(snf + 642);
    const auto *snf_646 = buffer.data(snf + 646);
    const auto *snf_648 = buffer.data(snf + 648);
    const auto *snf_649 = buffer.data(snf + 649);
    const auto *snf_650 = buffer.data(snf + 650);
    const auto *snf_652 = buffer.data(snf + 652);
    const auto *snf_656 = buffer.data(snf + 656);
    const auto *snf_658 = buffer.data(snf + 658);
    const auto *snf_659 = buffer.data(snf + 659);

    const auto *sng1_975 = buffer.data(sng1 + 975);
    const auto *sng1_980 = buffer.data(sng1 + 980);
    const auto *sng1_985 = buffer.data(sng1 + 985);
    const auto *sng1_987 = buffer.data(sng1 + 987);
    const auto *sng1_989 = buffer.data(sng1 + 989);

    const auto *sod0_431 = buffer.data(sod0 + 431);
    const auto *sod0_432 = buffer.data(sod0 + 432);
    const auto *sod0_435 = buffer.data(sod0 + 435);
    const auto *sod0_437 = buffer.data(sod0 + 437);
    const auto *sod0_438 = buffer.data(sod0 + 438);
    const auto *sod0_441 = buffer.data(sod0 + 441);
    const auto *sod0_443 = buffer.data(sod0 + 443);
    const auto *sod0_444 = buffer.data(sod0 + 444);
    const auto *sod0_447 = buffer.data(sod0 + 447);
    const auto *sod0_449 = buffer.data(sod0 + 449);
    const auto *sod0_450 = buffer.data(sod0 + 450);
    const auto *sod0_453 = buffer.data(sod0 + 453);
    const auto *sod0_455 = buffer.data(sod0 + 455);
    const auto *sod0_459 = buffer.data(sod0 + 459);
    const auto *sod0_462 = buffer.data(sod0 + 462);
    const auto *sod0_465 = buffer.data(sod0 + 465);
    const auto *sod0_467 = buffer.data(sod0 + 467);

    const auto *sod1_431 = buffer.data(sod1 + 431);
    const auto *sod1_432 = buffer.data(sod1 + 432);
    const auto *sod1_435 = buffer.data(sod1 + 435);
    const auto *sod1_437 = buffer.data(sod1 + 437);
    const auto *sod1_438 = buffer.data(sod1 + 438);
    const auto *sod1_441 = buffer.data(sod1 + 441);
    const auto *sod1_443 = buffer.data(sod1 + 443);
    const auto *sod1_444 = buffer.data(sod1 + 444);
    const auto *sod1_447 = buffer.data(sod1 + 447);
    const auto *sod1_449 = buffer.data(sod1 + 449);
    const auto *sod1_450 = buffer.data(sod1 + 450);
    const auto *sod1_453 = buffer.data(sod1 + 453);
    const auto *sod1_455 = buffer.data(sod1 + 455);
    const auto *sod1_459 = buffer.data(sod1 + 459);
    const auto *sod1_462 = buffer.data(sod1 + 462);
    const auto *sod1_465 = buffer.data(sod1 + 465);
    const auto *sod1_467 = buffer.data(sod1 + 467);

    const auto *sof_718 = buffer.data(sof + 718);
    const auto *sof_719 = buffer.data(sof + 719);
    const auto *sof_720 = buffer.data(sof + 720);
    const auto *sof_722 = buffer.data(sof + 722);
    const auto *sof_723 = buffer.data(sof + 723);
    const auto *sof_725 = buffer.data(sof + 725);
    const auto *sof_726 = buffer.data(sof + 726);
    const auto *sof_727 = buffer.data(sof + 727);
    const auto *sof_728 = buffer.data(sof + 728);
    const auto *sof_729 = buffer.data(sof + 729);
    const auto *sof_730 = buffer.data(sof + 730);
    const auto *sof_732 = buffer.data(sof + 732);
    const auto *sof_733 = buffer.data(sof + 733);
    const auto *sof_735 = buffer.data(sof + 735);
    const auto *sof_736 = buffer.data(sof + 736);
    const auto *sof_737 = buffer.data(sof + 737);
    const auto *sof_738 = buffer.data(sof + 738);
    const auto *sof_739 = buffer.data(sof + 739);
    const auto *sof_740 = buffer.data(sof + 740);
    const auto *sof_742 = buffer.data(sof + 742);
    const auto *sof_743 = buffer.data(sof + 743);
    const auto *sof_745 = buffer.data(sof + 745);
    const auto *sof_746 = buffer.data(sof + 746);
    const auto *sof_747 = buffer.data(sof + 747);
    const auto *sof_748 = buffer.data(sof + 748);
    const auto *sof_749 = buffer.data(sof + 749);
    const auto *sof_750 = buffer.data(sof + 750);
    const auto *sof_752 = buffer.data(sof + 752);
    const auto *sof_753 = buffer.data(sof + 753);
    const auto *sof_755 = buffer.data(sof + 755);
    const auto *sof_756 = buffer.data(sof + 756);
    const auto *sof_757 = buffer.data(sof + 757);
    const auto *sof_758 = buffer.data(sof + 758);
    const auto *sof_759 = buffer.data(sof + 759);
    const auto *sof_760 = buffer.data(sof + 760);
    const auto *sof_762 = buffer.data(sof + 762);
    const auto *sof_763 = buffer.data(sof + 763);
    const auto *sof_766 = buffer.data(sof + 766);
    const auto *sof_767 = buffer.data(sof + 767);
    const auto *sof_768 = buffer.data(sof + 768);
    const auto *sof_769 = buffer.data(sof + 769);
    const auto *sof_770 = buffer.data(sof + 770);
    const auto *sof_772 = buffer.data(sof + 772);
    const auto *sof_773 = buffer.data(sof + 773);
    const auto *sof_775 = buffer.data(sof + 775);
    const auto *sof_776 = buffer.data(sof + 776);
    const auto *sof_777 = buffer.data(sof + 777);
    const auto *sof_778 = buffer.data(sof + 778);
    const auto *sof_779 = buffer.data(sof + 779);

#pragma omp simd aligned(t_1077, t_1078, t_1079, pc_y, pc_z, snf_599, snf_608, snf_609, \
                         sod0_431, sod1_431, sof_718, sof_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = f_15 * snf_608[k]
                    + f_4 * sod0_431[k]
                    - f_5 * sod1_431[k]
                    + f_3 * pc_y[k] * sof_718[k];

        t_1078[k] = f_15 * snf_609[k]
                    + f_3 * pc_y[k] * sof_719[k];

        t_1079[k] = f_16 * snf_599[k]
                    + f_1 * sod0_431[k]
                    - f_2 * sod1_431[k]
                    + f_3 * pc_z[k] * sof_719[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, pc_x, pc_y, pc_z, snf_600, snf_610, \
                         sod0_432, sod0_435, sod1_432, sod1_435, sof_720, \
                         sof_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_1 * sod0_432[k]
                    - f_2 * sod1_432[k]
                    + f_3 * pc_x[k] * sof_720[k];

        t_1081[k] = f_16 * snf_610[k]
                    + f_3 * pc_y[k] * sof_720[k];

        t_1082[k] = f_15 * snf_600[k]
                    + f_3 * pc_z[k] * sof_720[k];

        t_1083[k] = f_4 * sod0_435[k]
                    - f_5 * sod1_435[k]
                    + f_3 * pc_x[k] * sof_723[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, t_1087, t_1088, pc_x, pc_y, snf_612, \
                         sod0_437, sod1_437, sof_722, sof_725, sof_726, sof_727, \
                         sof_728 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_16 * snf_612[k]
                    + f_3 * pc_y[k] * sof_722[k];

        t_1085[k] = f_4 * sod0_437[k]
                    - f_5 * sod1_437[k]
                    + f_3 * pc_x[k] * sof_725[k];

        t_1086[k] = f_3 * pc_x[k] * sof_726[k];

        t_1087[k] = f_3 * pc_x[k] * sof_727[k];

        t_1088[k] = f_3 * pc_x[k] * sof_728[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pc_x, pc_y, pc_z, snf_606, snf_616, sod0_435, \
                         sod1_435, sof_726, sof_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_3 * pc_x[k] * sof_729[k];

        t_1090[k] = f_16 * snf_616[k]
                    + f_1 * sod0_435[k]
                    - f_2 * sod1_435[k]
                    + f_3 * pc_y[k] * sof_726[k];

        t_1091[k] = f_15 * snf_606[k]
                    + f_3 * pc_z[k] * sof_726[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pc_y, pc_z, snf_609, snf_618, snf_619, \
                         sod0_437, sod1_437, sof_728, sof_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_16 * snf_618[k]
                    + f_4 * sod0_437[k]
                    - f_5 * sod1_437[k]
                    + f_3 * pc_y[k] * sof_728[k];

        t_1093[k] = f_16 * snf_619[k]
                    + f_3 * pc_y[k] * sof_729[k];

        t_1094[k] = f_15 * snf_609[k]
                    + f_1 * sod0_437[k]
                    - f_2 * sod1_437[k]
                    + f_3 * pc_z[k] * sof_729[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, t_1098, pc_x, pc_y, pc_z, snf_610, snf_620, \
                         sod0_438, sod0_441, sod1_438, sod1_441, sof_730, \
                         sof_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = f_1 * sod0_438[k]
                    - f_2 * sod1_438[k]
                    + f_3 * pc_x[k] * sof_730[k];

        t_1096[k] = f_14 * snf_620[k]
                    + f_3 * pc_y[k] * sof_730[k];

        t_1097[k] = f_13 * snf_610[k]
                    + f_3 * pc_z[k] * sof_730[k];

        t_1098[k] = f_4 * sod0_441[k]
                    - f_5 * sod1_441[k]
                    + f_3 * pc_x[k] * sof_733[k];
    }

#pragma omp simd aligned(t_1099, t_1100, t_1101, t_1102, t_1103, pc_x, pc_y, snf_622, \
                         sod0_443, sod1_443, sof_732, sof_735, sof_736, sof_737, \
                         sof_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1099[k] = f_14 * snf_622[k]
                    + f_3 * pc_y[k] * sof_732[k];

        t_1100[k] = f_4 * sod0_443[k]
                    - f_5 * sod1_443[k]
                    + f_3 * pc_x[k] * sof_735[k];

        t_1101[k] = f_3 * pc_x[k] * sof_736[k];

        t_1102[k] = f_3 * pc_x[k] * sof_737[k];

        t_1103[k] = f_3 * pc_x[k] * sof_738[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, pc_x, pc_y, pc_z, snf_616, snf_626, sod0_441, \
                         sod1_441, sof_736, sof_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_3 * pc_x[k] * sof_739[k];

        t_1105[k] = f_14 * snf_626[k]
                    + f_1 * sod0_441[k]
                    - f_2 * sod1_441[k]
                    + f_3 * pc_y[k] * sof_736[k];

        t_1106[k] = f_13 * snf_616[k]
                    + f_3 * pc_z[k] * sof_736[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, pc_y, pc_z, snf_619, snf_628, snf_629, \
                         sod0_443, sod1_443, sof_738, sof_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = f_14 * snf_628[k]
                    + f_4 * sod0_443[k]
                    - f_5 * sod1_443[k]
                    + f_3 * pc_y[k] * sof_738[k];

        t_1108[k] = f_14 * snf_629[k]
                    + f_3 * pc_y[k] * sof_739[k];

        t_1109[k] = f_13 * snf_619[k]
                    + f_1 * sod0_443[k]
                    - f_2 * sod1_443[k]
                    + f_3 * pc_z[k] * sof_739[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pc_x, pc_y, pc_z, snf_620, snf_630, \
                         sod0_444, sod0_447, sod1_444, sod1_447, sof_740, \
                         sof_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_1 * sod0_444[k]
                    - f_2 * sod1_444[k]
                    + f_3 * pc_x[k] * sof_740[k];

        t_1111[k] = f_12 * snf_630[k]
                    + f_3 * pc_y[k] * sof_740[k];

        t_1112[k] = f_11 * snf_620[k]
                    + f_3 * pc_z[k] * sof_740[k];

        t_1113[k] = f_4 * sod0_447[k]
                    - f_5 * sod1_447[k]
                    + f_3 * pc_x[k] * sof_743[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, t_1117, t_1118, pc_x, pc_y, snf_632, \
                         sod0_449, sod1_449, sof_742, sof_745, sof_746, sof_747, \
                         sof_748 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_12 * snf_632[k]
                    + f_3 * pc_y[k] * sof_742[k];

        t_1115[k] = f_4 * sod0_449[k]
                    - f_5 * sod1_449[k]
                    + f_3 * pc_x[k] * sof_745[k];

        t_1116[k] = f_3 * pc_x[k] * sof_746[k];

        t_1117[k] = f_3 * pc_x[k] * sof_747[k];

        t_1118[k] = f_3 * pc_x[k] * sof_748[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, pc_x, pc_y, pc_z, snf_626, snf_636, sod0_447, \
                         sod1_447, sof_746, sof_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = f_3 * pc_x[k] * sof_749[k];

        t_1120[k] = f_12 * snf_636[k]
                    + f_1 * sod0_447[k]
                    - f_2 * sod1_447[k]
                    + f_3 * pc_y[k] * sof_746[k];

        t_1121[k] = f_11 * snf_626[k]
                    + f_3 * pc_z[k] * sof_746[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, pc_y, pc_z, snf_629, snf_638, snf_639, \
                         sod0_449, sod1_449, sof_748, sof_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = f_12 * snf_638[k]
                    + f_4 * sod0_449[k]
                    - f_5 * sod1_449[k]
                    + f_3 * pc_y[k] * sof_748[k];

        t_1123[k] = f_12 * snf_639[k]
                    + f_3 * pc_y[k] * sof_749[k];

        t_1124[k] = f_11 * snf_629[k]
                    + f_1 * sod0_449[k]
                    - f_2 * sod1_449[k]
                    + f_3 * pc_z[k] * sof_749[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, t_1128, pc_x, pc_y, pc_z, snf_630, snf_640, \
                         sod0_450, sod0_453, sod1_450, sod1_453, sof_750, \
                         sof_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = f_1 * sod0_450[k]
                    - f_2 * sod1_450[k]
                    + f_3 * pc_x[k] * sof_750[k];

        t_1126[k] = f_8 * snf_640[k]
                    + f_3 * pc_y[k] * sof_750[k];

        t_1127[k] = f_10 * snf_630[k]
                    + f_3 * pc_z[k] * sof_750[k];

        t_1128[k] = f_4 * sod0_453[k]
                    - f_5 * sod1_453[k]
                    + f_3 * pc_x[k] * sof_753[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, t_1132, t_1133, pc_x, pc_y, snf_642, \
                         sod0_455, sod1_455, sof_752, sof_755, sof_756, sof_757, \
                         sof_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_8 * snf_642[k]
                    + f_3 * pc_y[k] * sof_752[k];

        t_1130[k] = f_4 * sod0_455[k]
                    - f_5 * sod1_455[k]
                    + f_3 * pc_x[k] * sof_755[k];

        t_1131[k] = f_3 * pc_x[k] * sof_756[k];

        t_1132[k] = f_3 * pc_x[k] * sof_757[k];

        t_1133[k] = f_3 * pc_x[k] * sof_758[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, pc_x, pc_y, pc_z, snf_636, snf_646, sod0_453, \
                         sod1_453, sof_756, sof_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_3 * pc_x[k] * sof_759[k];

        t_1135[k] = f_8 * snf_646[k]
                    + f_1 * sod0_453[k]
                    - f_2 * sod1_453[k]
                    + f_3 * pc_y[k] * sof_756[k];

        t_1136[k] = f_10 * snf_636[k]
                    + f_3 * pc_z[k] * sof_756[k];
    }

#pragma omp simd aligned(t_1137, t_1138, t_1139, t_1140, pb_y, pc_y, pc_z, sng0_975, snf_639, \
                         snf_648, snf_649, sng1_975, sod0_455, sod1_455, sof_758, \
                         sof_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = f_8 * snf_648[k]
                    + f_4 * sod0_455[k]
                    - f_5 * sod1_455[k]
                    + f_3 * pc_y[k] * sof_758[k];

        t_1138[k] = f_8 * snf_649[k]
                    + f_3 * pc_y[k] * sof_759[k];

        t_1139[k] = f_10 * snf_639[k]
                    + f_1 * sod0_455[k]
                    - f_2 * sod1_455[k]
                    + f_3 * pc_z[k] * sof_759[k];

        t_1140[k] = pb_y[k] * sng0_975[k]
                    - f_6 * pc_y[k] * sng1_975[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, t_1144, pc_x, pc_y, pc_z, snf_640, snf_650, \
                         snf_652, sod0_459, sod1_459, sof_760, sof_762, \
                         sof_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = f_7 * snf_650[k]
                    + f_3 * pc_y[k] * sof_760[k];

        t_1142[k] = f_9 * snf_640[k]
                    + f_3 * pc_z[k] * sof_760[k];

        t_1143[k] = f_4 * sod0_459[k]
                    - f_5 * sod1_459[k]
                    + f_3 * pc_x[k] * sof_763[k];

        t_1144[k] = f_7 * snf_652[k]
                    + f_3 * pc_y[k] * sof_762[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, t_1148, t_1149, pb_y, pc_x, pc_y, sng0_980, \
                         sng1_980, sof_766, sof_767, sof_768, sof_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = pb_y[k] * sng0_980[k]
                    - f_6 * pc_y[k] * sng1_980[k];

        t_1146[k] = f_3 * pc_x[k] * sof_766[k];

        t_1147[k] = f_3 * pc_x[k] * sof_767[k];

        t_1148[k] = f_3 * pc_x[k] * sof_768[k];

        t_1149[k] = f_3 * pc_x[k] * sof_769[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, pb_y, pc_y, pc_z, sng0_985, sng0_987, \
                         snf_646, snf_656, snf_658, sng1_985, sng1_987, \
                         sof_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = pb_y[k] * sng0_985[k]
                    + f_14 * snf_656[k]
                    - f_6 * pc_y[k] * sng1_985[k];

        t_1151[k] = f_9 * snf_646[k]
                    + f_3 * pc_z[k] * sof_766[k];

        t_1152[k] = pb_y[k] * sng0_987[k]
                    + f_8 * snf_658[k]
                    - f_6 * pc_y[k] * sng1_987[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, t_1156, pb_y, pc_x, pc_y, sng0_989, snf_659, \
                         sng1_989, sod0_462, sod1_462, sof_769, \
                         sof_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_7 * snf_659[k]
                    + f_3 * pc_y[k] * sof_769[k];

        t_1154[k] = pb_y[k] * sng0_989[k]
                    - f_6 * pc_y[k] * sng1_989[k];

        t_1155[k] = f_1 * sod0_462[k]
                    - f_2 * sod1_462[k]
                    + f_3 * pc_x[k] * sof_770[k];

        t_1156[k] = f_3 * pc_y[k] * sof_770[k];
    }

#pragma omp simd aligned(t_1157, t_1158, t_1159, t_1160, pc_x, pc_y, pc_z, snf_650, sod0_465, \
                         sod0_467, sod1_465, sod1_467, sof_770, sof_772, sof_773, \
                         sof_775 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1157[k] = f_0 * snf_650[k]
                    + f_3 * pc_z[k] * sof_770[k];

        t_1158[k] = f_4 * sod0_465[k]
                    - f_5 * sod1_465[k]
                    + f_3 * pc_x[k] * sof_773[k];

        t_1159[k] = f_3 * pc_y[k] * sof_772[k];

        t_1160[k] = f_4 * sod0_467[k]
                    - f_5 * sod1_467[k]
                    + f_3 * pc_x[k] * sof_775[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, t_1165, t_1166, pc_x, pc_y, pc_z, \
                         snf_656, sod0_465, sod1_465, sof_776, sof_777, sof_778, \
                         sof_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_3 * pc_x[k] * sof_776[k];

        t_1162[k] = f_3 * pc_x[k] * sof_777[k];

        t_1163[k] = f_3 * pc_x[k] * sof_778[k];

        t_1164[k] = f_3 * pc_x[k] * sof_779[k];

        t_1165[k] = f_1 * sod0_465[k]
                    - f_2 * sod1_465[k]
                    + f_3 * pc_y[k] * sof_776[k];

        t_1166[k] = f_0 * snf_656[k]
                    + f_3 * pc_z[k] * sof_776[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, pc_y, pc_z, snf_659, sod0_467, sod1_467, \
                         sof_778, sof_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = f_4 * sod0_467[k]
                    - f_5 * sod1_467[k]
                    + f_3 * pc_y[k] * sof_778[k];

        t_1168[k] = f_3 * pc_y[k] * sof_779[k];

        t_1169[k] = f_0 * snf_659[k]
                    + f_1 * sod0_467[k]
                    - f_2 * sod1_467[k]
                    + f_3 * pc_z[k] * sof_779[k];
    }
}

auto
compute_prim_sog_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sng0, const size_t snf,
                                                   const size_t sng1, const size_t sod0,
                                                   const size_t sod1, const size_t sof,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sog_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);

    compute_prim_sog_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, sng0, snf,
                                                              sng1, sod0, sod1, sof, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
