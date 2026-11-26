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


#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sfi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdi0,
                                                          const size_t sdh, const size_t sdi1,
                                                          const size_t sfg0, const size_t sfg1,
                                                          const size_t sfh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdi0_0 = buffer.data(sdi0 + 0);
    const auto *sdi0_3 = buffer.data(sdi0 + 3);
    const auto *sdi0_5 = buffer.data(sdi0 + 5);
    const auto *sdi0_6 = buffer.data(sdi0 + 6);
    const auto *sdi0_9 = buffer.data(sdi0 + 9);
    const auto *sdi0_10 = buffer.data(sdi0 + 10);
    const auto *sdi0_12 = buffer.data(sdi0 + 12);
    const auto *sdi0_14 = buffer.data(sdi0 + 14);
    const auto *sdi0_21 = buffer.data(sdi0 + 21);
    const auto *sdi0_27 = buffer.data(sdi0 + 27);
    const auto *sdi0_31 = buffer.data(sdi0 + 31);
    const auto *sdi0_34 = buffer.data(sdi0 + 34);
    const auto *sdi0_38 = buffer.data(sdi0 + 38);
    const auto *sdi0_56 = buffer.data(sdi0 + 56);
    const auto *sdi0_61 = buffer.data(sdi0 + 61);
    const auto *sdi0_65 = buffer.data(sdi0 + 65);
    const auto *sdi0_84 = buffer.data(sdi0 + 84);
    const auto *sdi0_87 = buffer.data(sdi0 + 87);
    const auto *sdi0_89 = buffer.data(sdi0 + 89);
    const auto *sdi0_90 = buffer.data(sdi0 + 90);
    const auto *sdi0_93 = buffer.data(sdi0 + 93);
    const auto *sdi0_94 = buffer.data(sdi0 + 94);
    const auto *sdi0_96 = buffer.data(sdi0 + 96);
    const auto *sdi0_98 = buffer.data(sdi0 + 98);
    const auto *sdi0_105 = buffer.data(sdi0 + 105);
    const auto *sdi0_107 = buffer.data(sdi0 + 107);
    const auto *sdi0_108 = buffer.data(sdi0 + 108);
    const auto *sdi0_109 = buffer.data(sdi0 + 109);
    const auto *sdi0_111 = buffer.data(sdi0 + 111);

    const auto *sdh_0 = buffer.data(sdh + 0);
    const auto *sdh_1 = buffer.data(sdh + 1);
    const auto *sdh_2 = buffer.data(sdh + 2);
    const auto *sdh_3 = buffer.data(sdh + 3);
    const auto *sdh_5 = buffer.data(sdh + 5);
    const auto *sdh_6 = buffer.data(sdh + 6);
    const auto *sdh_7 = buffer.data(sdh + 7);
    const auto *sdh_8 = buffer.data(sdh + 8);
    const auto *sdh_9 = buffer.data(sdh + 9);
    const auto *sdh_10 = buffer.data(sdh + 10);
    const auto *sdh_12 = buffer.data(sdh + 12);
    const auto *sdh_14 = buffer.data(sdh + 14);
    const auto *sdh_15 = buffer.data(sdh + 15);
    const auto *sdh_16 = buffer.data(sdh + 16);
    const auto *sdh_17 = buffer.data(sdh + 17);
    const auto *sdh_18 = buffer.data(sdh + 18);
    const auto *sdh_19 = buffer.data(sdh + 19);
    const auto *sdh_20 = buffer.data(sdh + 20);
    const auto *sdh_21 = buffer.data(sdh + 21);
    const auto *sdh_23 = buffer.data(sdh + 23);
    const auto *sdh_24 = buffer.data(sdh + 24);
    const auto *sdh_26 = buffer.data(sdh + 26);
    const auto *sdh_27 = buffer.data(sdh + 27);
    const auto *sdh_30 = buffer.data(sdh + 30);
    const auto *sdh_36 = buffer.data(sdh + 36);
    const auto *sdh_37 = buffer.data(sdh + 37);
    const auto *sdh_38 = buffer.data(sdh + 38);
    const auto *sdh_39 = buffer.data(sdh + 39);
    const auto *sdh_40 = buffer.data(sdh + 40);
    const auto *sdh_41 = buffer.data(sdh + 41);
    const auto *sdh_42 = buffer.data(sdh + 42);
    const auto *sdh_44 = buffer.data(sdh + 44);
    const auto *sdh_47 = buffer.data(sdh + 47);
    const auto *sdh_57 = buffer.data(sdh + 57);
    const auto *sdh_58 = buffer.data(sdh + 58);
    const auto *sdh_59 = buffer.data(sdh + 59);
    const auto *sdh_60 = buffer.data(sdh + 60);
    const auto *sdh_61 = buffer.data(sdh + 61);
    const auto *sdh_62 = buffer.data(sdh + 62);
    const auto *sdh_63 = buffer.data(sdh + 63);
    const auto *sdh_66 = buffer.data(sdh + 66);
    const auto *sdh_68 = buffer.data(sdh + 68);
    const auto *sdh_69 = buffer.data(sdh + 69);
    const auto *sdh_72 = buffer.data(sdh + 72);
    const auto *sdh_73 = buffer.data(sdh + 73);
    const auto *sdh_75 = buffer.data(sdh + 75);
    const auto *sdh_77 = buffer.data(sdh + 77);
    const auto *sdh_78 = buffer.data(sdh + 78);
    const auto *sdh_79 = buffer.data(sdh + 79);
    const auto *sdh_80 = buffer.data(sdh + 80);
    const auto *sdh_81 = buffer.data(sdh + 81);
    const auto *sdh_82 = buffer.data(sdh + 82);
    const auto *sdh_83 = buffer.data(sdh + 83);

    const auto *sdi1_0 = buffer.data(sdi1 + 0);
    const auto *sdi1_3 = buffer.data(sdi1 + 3);
    const auto *sdi1_5 = buffer.data(sdi1 + 5);
    const auto *sdi1_6 = buffer.data(sdi1 + 6);
    const auto *sdi1_9 = buffer.data(sdi1 + 9);
    const auto *sdi1_10 = buffer.data(sdi1 + 10);
    const auto *sdi1_12 = buffer.data(sdi1 + 12);
    const auto *sdi1_14 = buffer.data(sdi1 + 14);
    const auto *sdi1_21 = buffer.data(sdi1 + 21);
    const auto *sdi1_27 = buffer.data(sdi1 + 27);
    const auto *sdi1_31 = buffer.data(sdi1 + 31);
    const auto *sdi1_34 = buffer.data(sdi1 + 34);
    const auto *sdi1_38 = buffer.data(sdi1 + 38);
    const auto *sdi1_56 = buffer.data(sdi1 + 56);
    const auto *sdi1_61 = buffer.data(sdi1 + 61);
    const auto *sdi1_65 = buffer.data(sdi1 + 65);
    const auto *sdi1_84 = buffer.data(sdi1 + 84);
    const auto *sdi1_87 = buffer.data(sdi1 + 87);
    const auto *sdi1_89 = buffer.data(sdi1 + 89);
    const auto *sdi1_90 = buffer.data(sdi1 + 90);
    const auto *sdi1_93 = buffer.data(sdi1 + 93);
    const auto *sdi1_94 = buffer.data(sdi1 + 94);
    const auto *sdi1_96 = buffer.data(sdi1 + 96);
    const auto *sdi1_98 = buffer.data(sdi1 + 98);
    const auto *sdi1_105 = buffer.data(sdi1 + 105);
    const auto *sdi1_107 = buffer.data(sdi1 + 107);
    const auto *sdi1_108 = buffer.data(sdi1 + 108);
    const auto *sdi1_109 = buffer.data(sdi1 + 109);
    const auto *sdi1_111 = buffer.data(sdi1 + 111);

    const auto *sfg0_0 = buffer.data(sfg0 + 0);
    const auto *sfg0_3 = buffer.data(sfg0 + 3);
    const auto *sfg0_5 = buffer.data(sfg0 + 5);
    const auto *sfg0_6 = buffer.data(sfg0 + 6);
    const auto *sfg0_9 = buffer.data(sfg0 + 9);
    const auto *sfg0_10 = buffer.data(sfg0 + 10);
    const auto *sfg0_12 = buffer.data(sfg0 + 12);
    const auto *sfg0_13 = buffer.data(sfg0 + 13);
    const auto *sfg0_14 = buffer.data(sfg0 + 14);
    const auto *sfg0_25 = buffer.data(sfg0 + 25);
    const auto *sfg0_27 = buffer.data(sfg0 + 27);
    const auto *sfg0_28 = buffer.data(sfg0 + 28);
    const auto *sfg0_29 = buffer.data(sfg0 + 29);
    const auto *sfg0_42 = buffer.data(sfg0 + 42);
    const auto *sfg0_43 = buffer.data(sfg0 + 43);
    const auto *sfg0_44 = buffer.data(sfg0 + 44);

    const auto *sfg1_0 = buffer.data(sfg1 + 0);
    const auto *sfg1_3 = buffer.data(sfg1 + 3);
    const auto *sfg1_5 = buffer.data(sfg1 + 5);
    const auto *sfg1_6 = buffer.data(sfg1 + 6);
    const auto *sfg1_9 = buffer.data(sfg1 + 9);
    const auto *sfg1_10 = buffer.data(sfg1 + 10);
    const auto *sfg1_12 = buffer.data(sfg1 + 12);
    const auto *sfg1_13 = buffer.data(sfg1 + 13);
    const auto *sfg1_14 = buffer.data(sfg1 + 14);
    const auto *sfg1_25 = buffer.data(sfg1 + 25);
    const auto *sfg1_27 = buffer.data(sfg1 + 27);
    const auto *sfg1_28 = buffer.data(sfg1 + 28);
    const auto *sfg1_29 = buffer.data(sfg1 + 29);
    const auto *sfg1_42 = buffer.data(sfg1 + 42);
    const auto *sfg1_43 = buffer.data(sfg1 + 43);
    const auto *sfg1_44 = buffer.data(sfg1 + 44);

    const auto *sfh_0 = buffer.data(sfh + 0);
    const auto *sfh_2 = buffer.data(sfh + 2);
    const auto *sfh_3 = buffer.data(sfh + 3);
    const auto *sfh_5 = buffer.data(sfh + 5);
    const auto *sfh_6 = buffer.data(sfh + 6);
    const auto *sfh_9 = buffer.data(sfh + 9);
    const auto *sfh_10 = buffer.data(sfh + 10);
    const auto *sfh_12 = buffer.data(sfh + 12);
    const auto *sfh_14 = buffer.data(sfh + 14);
    const auto *sfh_15 = buffer.data(sfh + 15);
    const auto *sfh_16 = buffer.data(sfh + 16);
    const auto *sfh_17 = buffer.data(sfh + 17);
    const auto *sfh_18 = buffer.data(sfh + 18);
    const auto *sfh_19 = buffer.data(sfh + 19);
    const auto *sfh_20 = buffer.data(sfh + 20);
    const auto *sfh_21 = buffer.data(sfh + 21);
    const auto *sfh_23 = buffer.data(sfh + 23);
    const auto *sfh_24 = buffer.data(sfh + 24);
    const auto *sfh_26 = buffer.data(sfh + 26);
    const auto *sfh_27 = buffer.data(sfh + 27);
    const auto *sfh_30 = buffer.data(sfh + 30);
    const auto *sfh_36 = buffer.data(sfh + 36);
    const auto *sfh_37 = buffer.data(sfh + 37);
    const auto *sfh_38 = buffer.data(sfh + 38);
    const auto *sfh_39 = buffer.data(sfh + 39);
    const auto *sfh_40 = buffer.data(sfh + 40);
    const auto *sfh_41 = buffer.data(sfh + 41);
    const auto *sfh_42 = buffer.data(sfh + 42);
    const auto *sfh_44 = buffer.data(sfh + 44);
    const auto *sfh_45 = buffer.data(sfh + 45);
    const auto *sfh_47 = buffer.data(sfh + 47);
    const auto *sfh_48 = buffer.data(sfh + 48);
    const auto *sfh_51 = buffer.data(sfh + 51);
    const auto *sfh_57 = buffer.data(sfh + 57);
    const auto *sfh_58 = buffer.data(sfh + 58);
    const auto *sfh_59 = buffer.data(sfh + 59);
    const auto *sfh_60 = buffer.data(sfh + 60);
    const auto *sfh_61 = buffer.data(sfh + 61);
    const auto *sfh_62 = buffer.data(sfh + 62);
    const auto *sfh_63 = buffer.data(sfh + 63);
    const auto *sfh_65 = buffer.data(sfh + 65);
    const auto *sfh_66 = buffer.data(sfh + 66);
    const auto *sfh_68 = buffer.data(sfh + 68);
    const auto *sfh_69 = buffer.data(sfh + 69);
    const auto *sfh_72 = buffer.data(sfh + 72);
    const auto *sfh_78 = buffer.data(sfh + 78);
    const auto *sfh_79 = buffer.data(sfh + 79);
    const auto *sfh_80 = buffer.data(sfh + 80);
    const auto *sfh_81 = buffer.data(sfh + 81);
    const auto *sfh_82 = buffer.data(sfh + 82);
    const auto *sfh_83 = buffer.data(sfh + 83);
    const auto *sfh_84 = buffer.data(sfh + 84);
    const auto *sfh_86 = buffer.data(sfh + 86);
    const auto *sfh_87 = buffer.data(sfh + 87);
    const auto *sfh_89 = buffer.data(sfh + 89);
    const auto *sfh_90 = buffer.data(sfh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sdh_0, sdh_3, sfg0_0, sfg0_3, \
                         sfg1_0, sfg1_3, sfh_0, sfh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdh_0[k]
                 + f_1 * sfg0_0[k]
                 - f_2 * sfg1_0[k]
                 + f_3 * pc_x[k] * sfh_0[k];

        t_1[k] = f_3 * pc_y[k] * sfh_0[k];

        t_2[k] = f_3 * pc_z[k] * sfh_0[k];

        t_3[k] = f_0 * sdh_3[k]
                 + f_4 * sfg0_3[k]
                 - f_5 * sfg1_3[k]
                 + f_3 * pc_x[k] * sfh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sdh_5, sdh_6, sfg0_5, sfg0_6, sfg1_5, \
                         sfg1_6, sfh_2, sfh_5, sfh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sfh_2[k];

        t_5[k] = f_0 * sdh_5[k]
                 + f_4 * sfg0_5[k]
                 - f_5 * sfg1_5[k]
                 + f_3 * pc_x[k] * sfh_5[k];

        t_6[k] = f_0 * sdh_6[k]
                 + f_6 * sfg0_6[k]
                 - f_7 * sfg1_6[k]
                 + f_3 * pc_x[k] * sfh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sdh_9, sfg0_9, sfg1_9, sfh_3, sfh_5, \
                         sfh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sfh_3[k];

        t_8[k] = f_3 * pc_y[k] * sfh_5[k];

        t_9[k] = f_0 * sdh_9[k]
                 + f_6 * sfg0_9[k]
                 - f_7 * sfg1_9[k]
                 + f_3 * pc_x[k] * sfh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sdh_10, sdh_12, sfg0_10, sfg0_12, \
                         sfg1_10, sfg1_12, sfh_6, sfh_10, sfh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sdh_10[k]
                  + f_8 * sfg0_10[k]
                  - f_9 * sfg1_10[k]
                  + f_3 * pc_x[k] * sfh_10[k];

        t_11[k] = f_3 * pc_z[k] * sfh_6[k];

        t_12[k] = f_0 * sdh_12[k]
                  + f_8 * sfg0_12[k]
                  - f_9 * sfg1_12[k]
                  + f_3 * pc_x[k] * sfh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, sdh_14, sdh_15, sdh_16, sfg0_14, \
                         sfg1_14, sfh_9, sfh_14, sfh_15, sfh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sfh_9[k];

        t_14[k] = f_0 * sdh_14[k]
                  + f_8 * sfg0_14[k]
                  - f_9 * sfg1_14[k]
                  + f_3 * pc_x[k] * sfh_14[k];

        t_15[k] = f_0 * sdh_15[k]
                  + f_3 * pc_x[k] * sfh_15[k];

        t_16[k] = f_0 * sdh_16[k]
                  + f_3 * pc_x[k] * sfh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, sdh_17, sdh_18, sdh_19, sdh_20, sfh_17, \
                         sfh_18, sfh_19, sfh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * sdh_17[k]
                  + f_3 * pc_x[k] * sfh_17[k];

        t_18[k] = f_0 * sdh_18[k]
                  + f_3 * pc_x[k] * sfh_18[k];

        t_19[k] = f_0 * sdh_19[k]
                  + f_3 * pc_x[k] * sfh_19[k];

        t_20[k] = f_0 * sdh_20[k]
                  + f_3 * pc_x[k] * sfh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, sfg0_10, sfg0_12, sfg0_13, \
                         sfg1_10, sfg1_12, sfg1_13, sfh_15, sfh_17, \
                         sfh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * sfg0_10[k]
                  - f_2 * sfg1_10[k]
                  + f_3 * pc_y[k] * sfh_15[k];

        t_22[k] = f_3 * pc_z[k] * sfh_15[k];

        t_23[k] = f_4 * sfg0_12[k]
                  - f_5 * sfg1_12[k]
                  + f_3 * pc_y[k] * sfh_17[k];

        t_24[k] = f_6 * sfg0_13[k]
                  - f_7 * sfg1_13[k]
                  + f_3 * pc_y[k] * sfh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, sdi0_0, sdh_0, \
                         sdi1_0, sfg0_14, sfg1_14, sfh_19, sfh_20, \
                         sfh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * sfg0_14[k]
                  - f_9 * sfg1_14[k]
                  + f_3 * pc_y[k] * sfh_19[k];

        t_26[k] = f_3 * pc_y[k] * sfh_20[k];

        t_27[k] = f_1 * sfg0_14[k]
                  - f_2 * sfg1_14[k]
                  + f_3 * pc_z[k] * sfh_20[k];

        t_28[k] = pb_y[k] * sdi0_0[k]
                  - f_10 * pc_y[k] * sdi1_0[k];

        t_29[k] = f_11 * sdh_0[k]
                  + f_3 * pc_y[k] * sfh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, sdi0_3, sdi0_5, sdh_1, \
                         sdh_2, sdi1_3, sdi1_5, sfh_21, sfh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * sfh_21[k];

        t_31[k] = pb_y[k] * sdi0_3[k]
                  + f_12 * sdh_1[k]
                  - f_10 * pc_y[k] * sdi1_3[k];

        t_32[k] = f_11 * sdh_2[k]
                  + f_3 * pc_y[k] * sfh_23[k];

        t_33[k] = pb_y[k] * sdi0_5[k]
                  - f_10 * pc_y[k] * sdi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, sdi0_6, sdi0_9, sdh_3, \
                         sdh_5, sdi1_6, sdi1_9, sfh_24, sfh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * sdi0_6[k]
                  + f_0 * sdh_3[k]
                  - f_10 * pc_y[k] * sdi1_6[k];

        t_35[k] = f_3 * pc_z[k] * sfh_24[k];

        t_36[k] = f_11 * sdh_5[k]
                  + f_3 * pc_y[k] * sfh_26[k];

        t_37[k] = pb_y[k] * sdi0_9[k]
                  - f_10 * pc_y[k] * sdi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sdi0_10, sdi0_12, sdh_6, \
                         sdh_8, sdh_9, sdi1_10, sdi1_12, sfh_27, \
                         sfh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * sdi0_10[k]
                  + f_13 * sdh_6[k]
                  - f_10 * pc_y[k] * sdi1_10[k];

        t_39[k] = f_3 * pc_z[k] * sfh_27[k];

        t_40[k] = pb_y[k] * sdi0_12[k]
                  + f_12 * sdh_8[k]
                  - f_10 * pc_y[k] * sdi1_12[k];

        t_41[k] = f_11 * sdh_9[k]
                  + f_3 * pc_y[k] * sfh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, sdi0_14, sdh_36, sdh_37, \
                         sdh_38, sdi1_14, sfh_36, sfh_37, sfh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * sdi0_14[k]
                  - f_10 * pc_y[k] * sdi1_14[k];

        t_43[k] = f_12 * sdh_36[k]
                  + f_3 * pc_x[k] * sfh_36[k];

        t_44[k] = f_12 * sdh_37[k]
                  + f_3 * pc_x[k] * sfh_37[k];

        t_45[k] = f_12 * sdh_38[k]
                  + f_3 * pc_x[k] * sfh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, sdh_15, sdh_39, sdh_40, sdh_41, \
                         sfg0_25, sfg1_25, sfh_36, sfh_39, sfh_40, \
                         sfh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * sdh_39[k]
                  + f_3 * pc_x[k] * sfh_39[k];

        t_47[k] = f_12 * sdh_40[k]
                  + f_3 * pc_x[k] * sfh_40[k];

        t_48[k] = f_12 * sdh_41[k]
                  + f_3 * pc_x[k] * sfh_41[k];

        t_49[k] = f_11 * sdh_15[k]
                  + f_1 * sfg0_25[k]
                  - f_2 * sfg1_25[k]
                  + f_3 * pc_y[k] * sfh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, sdh_17, sdh_18, sfg0_27, sfg0_28, \
                         sfg1_27, sfg1_28, sfh_36, sfh_38, sfh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * sfh_36[k];

        t_51[k] = f_11 * sdh_17[k]
                  + f_4 * sfg0_27[k]
                  - f_5 * sfg1_27[k]
                  + f_3 * pc_y[k] * sfh_38[k];

        t_52[k] = f_11 * sdh_18[k]
                  + f_6 * sfg0_28[k]
                  - f_7 * sfg1_28[k]
                  + f_3 * pc_y[k] * sfh_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, sdi0_27, sdh_19, sdh_20, sdi1_27, \
                         sfg0_29, sfg1_29, sfh_40, sfh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * sdh_19[k]
                  + f_8 * sfg0_29[k]
                  - f_9 * sfg1_29[k]
                  + f_3 * pc_y[k] * sfh_40[k];

        t_54[k] = f_11 * sdh_20[k]
                  + f_3 * pc_y[k] * sfh_41[k];

        t_55[k] = pb_y[k] * sdi0_27[k]
                  - f_10 * pc_y[k] * sdi1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, sdi0_0, sdi0_3, \
                         sdh_0, sdi1_0, sdi1_3, sfh_42, sfh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * sdi0_0[k]
                  - f_10 * pc_z[k] * sdi1_0[k];

        t_57[k] = f_3 * pc_y[k] * sfh_42[k];

        t_58[k] = f_11 * sdh_0[k]
                  + f_3 * pc_z[k] * sfh_42[k];

        t_59[k] = pb_z[k] * sdi0_3[k]
                  - f_10 * pc_z[k] * sdi1_3[k];

        t_60[k] = f_3 * pc_y[k] * sfh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sdi0_5, sdi0_6, sdh_2, \
                         sdh_3, sdi1_5, sdi1_6, sfh_45, sfh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * sdi0_5[k]
                  + f_12 * sdh_2[k]
                  - f_10 * pc_z[k] * sdi1_5[k];

        t_62[k] = pb_z[k] * sdi0_6[k]
                  - f_10 * pc_z[k] * sdi1_6[k];

        t_63[k] = f_11 * sdh_3[k]
                  + f_3 * pc_z[k] * sfh_45[k];

        t_64[k] = f_3 * pc_y[k] * sfh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, sdi0_9, sdi0_10, sdi0_12, sdh_5, \
                         sdh_6, sdh_7, sdi1_9, sdi1_10, sdi1_12, \
                         sfh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * sdi0_9[k]
                  + f_0 * sdh_5[k]
                  - f_10 * pc_z[k] * sdi1_9[k];

        t_66[k] = pb_z[k] * sdi0_10[k]
                  - f_10 * pc_z[k] * sdi1_10[k];

        t_67[k] = f_11 * sdh_6[k]
                  + f_3 * pc_z[k] * sfh_48[k];

        t_68[k] = pb_z[k] * sdi0_12[k]
                  + f_12 * sdh_7[k]
                  - f_10 * pc_z[k] * sdi1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, sdi0_14, sdh_9, \
                         sdh_57, sdh_58, sdi1_14, sfh_51, sfh_57, \
                         sfh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * sfh_51[k];

        t_70[k] = pb_z[k] * sdi0_14[k]
                  + f_13 * sdh_9[k]
                  - f_10 * pc_z[k] * sdi1_14[k];

        t_71[k] = f_12 * sdh_57[k]
                  + f_3 * pc_x[k] * sfh_57[k];

        t_72[k] = f_12 * sdh_58[k]
                  + f_3 * pc_x[k] * sfh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, sdh_59, sdh_60, sdh_61, sdh_62, sfh_59, \
                         sfh_60, sfh_61, sfh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_12 * sdh_59[k]
                  + f_3 * pc_x[k] * sfh_59[k];

        t_74[k] = f_12 * sdh_60[k]
                  + f_3 * pc_x[k] * sfh_60[k];

        t_75[k] = f_12 * sdh_61[k]
                  + f_3 * pc_x[k] * sfh_61[k];

        t_76[k] = f_12 * sdh_62[k]
                  + f_3 * pc_x[k] * sfh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, sdi0_21, sdh_15, sdi1_21, \
                         sfg0_42, sfg1_42, sfh_57, sfh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * sdi0_21[k]
                  - f_10 * pc_z[k] * sdi1_21[k];

        t_78[k] = f_11 * sdh_15[k]
                  + f_3 * pc_z[k] * sfh_57[k];

        t_79[k] = f_4 * sfg0_42[k]
                  - f_5 * sfg1_42[k]
                  + f_3 * pc_y[k] * sfh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, sdh_20, sfg0_43, sfg0_44, \
                         sfg1_43, sfg1_44, sfh_60, sfh_61, sfh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * sfg0_43[k]
                  - f_7 * sfg1_43[k]
                  + f_3 * pc_y[k] * sfh_60[k];

        t_81[k] = f_8 * sfg0_44[k]
                  - f_9 * sfg1_44[k]
                  + f_3 * pc_y[k] * sfh_61[k];

        t_82[k] = f_3 * pc_y[k] * sfh_62[k];

        t_83[k] = f_11 * sdh_20[k]
                  + f_1 * sfg0_44[k]
                  - f_2 * sfg1_44[k]
                  + f_3 * pc_z[k] * sfh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pc_x, pc_y, pc_z, sdi0_84, sdi0_87, \
                         sdh_21, sdh_63, sdh_66, sdi1_84, sdi1_87, \
                         sfh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pb_x[k] * sdi0_84[k]
                  + f_14 * sdh_63[k]
                  - f_10 * pc_x[k] * sdi1_84[k];

        t_85[k] = f_12 * sdh_21[k]
                  + f_3 * pc_y[k] * sfh_63[k];

        t_86[k] = f_3 * pc_z[k] * sfh_63[k];

        t_87[k] = pb_x[k] * sdi0_87[k]
                  + f_13 * sdh_66[k]
                  - f_10 * pc_x[k] * sdi1_87[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pb_x, pc_x, pc_y, sdi0_89, sdi0_90, sdh_23, sdh_68, \
                         sdh_69, sdi1_89, sdi1_90, sfh_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * sdh_23[k]
                  + f_3 * pc_y[k] * sfh_65[k];

        t_89[k] = pb_x[k] * sdi0_89[k]
                  + f_13 * sdh_68[k]
                  - f_10 * pc_x[k] * sdi1_89[k];

        t_90[k] = pb_x[k] * sdi0_90[k]
                  + f_0 * sdh_69[k]
                  - f_10 * pc_x[k] * sdi1_90[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_x, pc_x, pc_y, pc_z, sdi0_93, sdh_26, sdh_72, \
                         sdi1_93, sfh_66, sfh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * sfh_66[k];

        t_92[k] = f_12 * sdh_26[k]
                  + f_3 * pc_y[k] * sfh_68[k];

        t_93[k] = pb_x[k] * sdi0_93[k]
                  + f_0 * sdh_72[k]
                  - f_10 * pc_x[k] * sdi1_93[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_x, pc_x, pc_z, sdi0_94, sdi0_96, sdh_73, sdh_75, \
                         sdi1_94, sdi1_96, sfh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = pb_x[k] * sdi0_94[k]
                  + f_12 * sdh_73[k]
                  - f_10 * pc_x[k] * sdi1_94[k];

        t_95[k] = f_3 * pc_z[k] * sfh_69[k];

        t_96[k] = pb_x[k] * sdi0_96[k]
                  + f_12 * sdh_75[k]
                  - f_10 * pc_x[k] * sdi1_96[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pc_x, pc_y, sdi0_98, sdh_30, sdh_77, \
                         sdh_78, sdh_79, sdi1_98, sfh_72, sfh_78, \
                         sfh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * sdh_30[k]
                  + f_3 * pc_y[k] * sfh_72[k];

        t_98[k] = pb_x[k] * sdi0_98[k]
                  + f_12 * sdh_77[k]
                  - f_10 * pc_x[k] * sdi1_98[k];

        t_99[k] = f_11 * sdh_78[k]
                  + f_3 * pc_x[k] * sfh_78[k];

        t_100[k] = f_11 * sdh_79[k]
                   + f_3 * pc_x[k] * sfh_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, sdh_80, sdh_81, sdh_82, sdh_83, \
                         sfh_80, sfh_81, sfh_82, sfh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_11 * sdh_80[k]
                   + f_3 * pc_x[k] * sfh_80[k];

        t_102[k] = f_11 * sdh_81[k]
                   + f_3 * pc_x[k] * sfh_81[k];

        t_103[k] = f_11 * sdh_82[k]
                   + f_3 * pc_x[k] * sfh_82[k];

        t_104[k] = f_11 * sdh_83[k]
                   + f_3 * pc_x[k] * sfh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, pc_x, pc_z, sdi0_105, sdi0_107, \
                         sdi0_108, sdi1_105, sdi1_107, sdi1_108, \
                         sfh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pb_x[k] * sdi0_105[k]
                   - f_10 * pc_x[k] * sdi1_105[k];

        t_106[k] = f_3 * pc_z[k] * sfh_78[k];

        t_107[k] = pb_x[k] * sdi0_107[k]
                   - f_10 * pc_x[k] * sdi1_107[k];

        t_108[k] = pb_x[k] * sdi0_108[k]
                   - f_10 * pc_x[k] * sdi1_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_x, pb_y, pc_x, pc_y, sdi0_56, \
                         sdi0_109, sdi0_111, sdh_41, sdi1_56, sdi1_109, sdi1_111, \
                         sfh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * sdi0_109[k]
                   - f_10 * pc_x[k] * sdi1_109[k];

        t_110[k] = f_12 * sdh_41[k]
                   + f_3 * pc_y[k] * sfh_83[k];

        t_111[k] = pb_x[k] * sdi0_111[k]
                   - f_10 * pc_x[k] * sdi1_111[k];

        t_112[k] = pb_y[k] * sdi0_56[k]
                   - f_10 * pc_y[k] * sdi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_z, pc_y, pc_z, sdi0_31, sdh_21, \
                         sdh_42, sdh_44, sdi1_31, sfh_84, sfh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * sdh_42[k]
                   + f_3 * pc_y[k] * sfh_84[k];

        t_114[k] = f_11 * sdh_21[k]
                   + f_3 * pc_z[k] * sfh_84[k];

        t_115[k] = pb_z[k] * sdi0_31[k]
                   - f_10 * pc_z[k] * sdi1_31[k];

        t_116[k] = f_11 * sdh_44[k]
                   + f_3 * pc_y[k] * sfh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pb_z, pc_y, pc_z, sdi0_34, sdi0_61, \
                         sdh_24, sdh_47, sdi1_34, sdi1_61, sfh_87, \
                         sfh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pb_y[k] * sdi0_61[k]
                   - f_10 * pc_y[k] * sdi1_61[k];

        t_118[k] = pb_z[k] * sdi0_34[k]
                   - f_10 * pc_z[k] * sdi1_34[k];

        t_119[k] = f_11 * sdh_24[k]
                   + f_3 * pc_z[k] * sfh_87[k];

        t_120[k] = f_11 * sdh_47[k]
                   + f_3 * pc_y[k] * sfh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, sdi0_38, sdi0_65, \
                         sdh_27, sdi1_38, sdi1_65, sfh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pb_y[k] * sdi0_65[k]
                   - f_10 * pc_y[k] * sdi1_65[k];

        t_122[k] = pb_z[k] * sdi0_38[k]
                   - f_10 * pc_z[k] * sdi1_38[k];

        t_123[k] = f_11 * sdh_27[k]
                   + f_3 * pc_z[k] * sfh_90[k];
    }
}

static auto
compute_prim_sfi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdi0,
                                                          const size_t sdh, const size_t sdi1,
                                                          const size_t sfg0, const size_t sfg1,
                                                          const size_t sfh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdi0_70 = buffer.data(sdi0 + 70);
    const auto *sdi0_84 = buffer.data(sdi0 + 84);
    const auto *sdi0_87 = buffer.data(sdi0 + 87);
    const auto *sdi0_90 = buffer.data(sdi0 + 90);
    const auto *sdi0_94 = buffer.data(sdi0 + 94);
    const auto *sdi0_105 = buffer.data(sdi0 + 105);
    const auto *sdi0_107 = buffer.data(sdi0 + 107);
    const auto *sdi0_108 = buffer.data(sdi0 + 108);
    const auto *sdi0_109 = buffer.data(sdi0 + 109);
    const auto *sdi0_124 = buffer.data(sdi0 + 124);
    const auto *sdi0_133 = buffer.data(sdi0 + 133);
    const auto *sdi0_135 = buffer.data(sdi0 + 135);
    const auto *sdi0_136 = buffer.data(sdi0 + 136);
    const auto *sdi0_137 = buffer.data(sdi0 + 137);
    const auto *sdi0_139 = buffer.data(sdi0 + 139);
    const auto *sdi0_140 = buffer.data(sdi0 + 140);
    const auto *sdi0_143 = buffer.data(sdi0 + 143);
    const auto *sdi0_145 = buffer.data(sdi0 + 145);
    const auto *sdi0_146 = buffer.data(sdi0 + 146);
    const auto *sdi0_149 = buffer.data(sdi0 + 149);
    const auto *sdi0_150 = buffer.data(sdi0 + 150);
    const auto *sdi0_152 = buffer.data(sdi0 + 152);
    const auto *sdi0_154 = buffer.data(sdi0 + 154);
    const auto *sdi0_161 = buffer.data(sdi0 + 161);
    const auto *sdi0_163 = buffer.data(sdi0 + 163);
    const auto *sdi0_164 = buffer.data(sdi0 + 164);
    const auto *sdi0_165 = buffer.data(sdi0 + 165);
    const auto *sdi0_167 = buffer.data(sdi0 + 167);

    const auto *sdh_36 = buffer.data(sdh + 36);
    const auto *sdh_42 = buffer.data(sdh + 42);
    const auto *sdh_45 = buffer.data(sdh + 45);
    const auto *sdh_48 = buffer.data(sdh + 48);
    const auto *sdh_51 = buffer.data(sdh + 51);
    const auto *sdh_57 = buffer.data(sdh + 57);
    const auto *sdh_62 = buffer.data(sdh + 62);
    const auto *sdh_63 = buffer.data(sdh + 63);
    const auto *sdh_65 = buffer.data(sdh + 65);
    const auto *sdh_66 = buffer.data(sdh + 66);
    const auto *sdh_68 = buffer.data(sdh + 68);
    const auto *sdh_69 = buffer.data(sdh + 69);
    const auto *sdh_72 = buffer.data(sdh + 72);
    const auto *sdh_78 = buffer.data(sdh + 78);
    const auto *sdh_79 = buffer.data(sdh + 79);
    const auto *sdh_80 = buffer.data(sdh + 80);
    const auto *sdh_81 = buffer.data(sdh + 81);
    const auto *sdh_82 = buffer.data(sdh + 82);
    const auto *sdh_83 = buffer.data(sdh + 83);
    const auto *sdh_84 = buffer.data(sdh + 84);
    const auto *sdh_86 = buffer.data(sdh + 86);
    const auto *sdh_87 = buffer.data(sdh + 87);
    const auto *sdh_89 = buffer.data(sdh + 89);
    const auto *sdh_90 = buffer.data(sdh + 90);
    const auto *sdh_93 = buffer.data(sdh + 93);
    const auto *sdh_96 = buffer.data(sdh + 96);
    const auto *sdh_99 = buffer.data(sdh + 99);
    const auto *sdh_100 = buffer.data(sdh + 100);
    const auto *sdh_101 = buffer.data(sdh + 101);
    const auto *sdh_102 = buffer.data(sdh + 102);
    const auto *sdh_103 = buffer.data(sdh + 103);
    const auto *sdh_104 = buffer.data(sdh + 104);
    const auto *sdh_105 = buffer.data(sdh + 105);
    const auto *sdh_107 = buffer.data(sdh + 107);
    const auto *sdh_108 = buffer.data(sdh + 108);
    const auto *sdh_110 = buffer.data(sdh + 110);
    const auto *sdh_111 = buffer.data(sdh + 111);
    const auto *sdh_114 = buffer.data(sdh + 114);
    const auto *sdh_115 = buffer.data(sdh + 115);
    const auto *sdh_117 = buffer.data(sdh + 117);
    const auto *sdh_119 = buffer.data(sdh + 119);
    const auto *sdh_120 = buffer.data(sdh + 120);
    const auto *sdh_121 = buffer.data(sdh + 121);
    const auto *sdh_122 = buffer.data(sdh + 122);
    const auto *sdh_123 = buffer.data(sdh + 123);
    const auto *sdh_124 = buffer.data(sdh + 124);
    const auto *sdh_125 = buffer.data(sdh + 125);

    const auto *sdi1_70 = buffer.data(sdi1 + 70);
    const auto *sdi1_84 = buffer.data(sdi1 + 84);
    const auto *sdi1_87 = buffer.data(sdi1 + 87);
    const auto *sdi1_90 = buffer.data(sdi1 + 90);
    const auto *sdi1_94 = buffer.data(sdi1 + 94);
    const auto *sdi1_105 = buffer.data(sdi1 + 105);
    const auto *sdi1_107 = buffer.data(sdi1 + 107);
    const auto *sdi1_108 = buffer.data(sdi1 + 108);
    const auto *sdi1_109 = buffer.data(sdi1 + 109);
    const auto *sdi1_124 = buffer.data(sdi1 + 124);
    const auto *sdi1_133 = buffer.data(sdi1 + 133);
    const auto *sdi1_135 = buffer.data(sdi1 + 135);
    const auto *sdi1_136 = buffer.data(sdi1 + 136);
    const auto *sdi1_137 = buffer.data(sdi1 + 137);
    const auto *sdi1_139 = buffer.data(sdi1 + 139);
    const auto *sdi1_140 = buffer.data(sdi1 + 140);
    const auto *sdi1_143 = buffer.data(sdi1 + 143);
    const auto *sdi1_145 = buffer.data(sdi1 + 145);
    const auto *sdi1_146 = buffer.data(sdi1 + 146);
    const auto *sdi1_149 = buffer.data(sdi1 + 149);
    const auto *sdi1_150 = buffer.data(sdi1 + 150);
    const auto *sdi1_152 = buffer.data(sdi1 + 152);
    const auto *sdi1_154 = buffer.data(sdi1 + 154);
    const auto *sdi1_161 = buffer.data(sdi1 + 161);
    const auto *sdi1_163 = buffer.data(sdi1 + 163);
    const auto *sdi1_164 = buffer.data(sdi1 + 164);
    const auto *sdi1_165 = buffer.data(sdi1 + 165);
    const auto *sdi1_167 = buffer.data(sdi1 + 167);

    const auto *sfg0_90 = buffer.data(sfg0 + 90);
    const auto *sfg0_93 = buffer.data(sfg0 + 93);
    const auto *sfg0_95 = buffer.data(sfg0 + 95);
    const auto *sfg0_96 = buffer.data(sfg0 + 96);
    const auto *sfg0_99 = buffer.data(sfg0 + 99);
    const auto *sfg0_100 = buffer.data(sfg0 + 100);
    const auto *sfg0_102 = buffer.data(sfg0 + 102);
    const auto *sfg0_103 = buffer.data(sfg0 + 103);
    const auto *sfg0_104 = buffer.data(sfg0 + 104);
    const auto *sfg0_110 = buffer.data(sfg0 + 110);
    const auto *sfg0_114 = buffer.data(sfg0 + 114);
    const auto *sfg0_117 = buffer.data(sfg0 + 117);
    const auto *sfg0_119 = buffer.data(sfg0 + 119);
    const auto *sfg0_123 = buffer.data(sfg0 + 123);
    const auto *sfg0_126 = buffer.data(sfg0 + 126);
    const auto *sfg0_130 = buffer.data(sfg0 + 130);
    const auto *sfg0_132 = buffer.data(sfg0 + 132);

    const auto *sfg1_90 = buffer.data(sfg1 + 90);
    const auto *sfg1_93 = buffer.data(sfg1 + 93);
    const auto *sfg1_95 = buffer.data(sfg1 + 95);
    const auto *sfg1_96 = buffer.data(sfg1 + 96);
    const auto *sfg1_99 = buffer.data(sfg1 + 99);
    const auto *sfg1_100 = buffer.data(sfg1 + 100);
    const auto *sfg1_102 = buffer.data(sfg1 + 102);
    const auto *sfg1_103 = buffer.data(sfg1 + 103);
    const auto *sfg1_104 = buffer.data(sfg1 + 104);
    const auto *sfg1_110 = buffer.data(sfg1 + 110);
    const auto *sfg1_114 = buffer.data(sfg1 + 114);
    const auto *sfg1_117 = buffer.data(sfg1 + 117);
    const auto *sfg1_119 = buffer.data(sfg1 + 119);
    const auto *sfg1_123 = buffer.data(sfg1 + 123);
    const auto *sfg1_126 = buffer.data(sfg1 + 126);
    const auto *sfg1_130 = buffer.data(sfg1 + 130);
    const auto *sfg1_132 = buffer.data(sfg1 + 132);

    const auto *sfh_93 = buffer.data(sfh + 93);
    const auto *sfh_99 = buffer.data(sfh + 99);
    const auto *sfh_100 = buffer.data(sfh + 100);
    const auto *sfh_101 = buffer.data(sfh + 101);
    const auto *sfh_102 = buffer.data(sfh + 102);
    const auto *sfh_103 = buffer.data(sfh + 103);
    const auto *sfh_104 = buffer.data(sfh + 104);
    const auto *sfh_105 = buffer.data(sfh + 105);
    const auto *sfh_107 = buffer.data(sfh + 107);
    const auto *sfh_108 = buffer.data(sfh + 108);
    const auto *sfh_110 = buffer.data(sfh + 110);
    const auto *sfh_111 = buffer.data(sfh + 111);
    const auto *sfh_114 = buffer.data(sfh + 114);
    const auto *sfh_120 = buffer.data(sfh + 120);
    const auto *sfh_121 = buffer.data(sfh + 121);
    const auto *sfh_122 = buffer.data(sfh + 122);
    const auto *sfh_123 = buffer.data(sfh + 123);
    const auto *sfh_124 = buffer.data(sfh + 124);
    const auto *sfh_125 = buffer.data(sfh + 125);
    const auto *sfh_126 = buffer.data(sfh + 126);
    const auto *sfh_128 = buffer.data(sfh + 128);
    const auto *sfh_129 = buffer.data(sfh + 129);
    const auto *sfh_131 = buffer.data(sfh + 131);
    const auto *sfh_132 = buffer.data(sfh + 132);
    const auto *sfh_135 = buffer.data(sfh + 135);
    const auto *sfh_136 = buffer.data(sfh + 136);
    const auto *sfh_138 = buffer.data(sfh + 138);
    const auto *sfh_140 = buffer.data(sfh + 140);
    const auto *sfh_141 = buffer.data(sfh + 141);
    const auto *sfh_142 = buffer.data(sfh + 142);
    const auto *sfh_143 = buffer.data(sfh + 143);
    const auto *sfh_144 = buffer.data(sfh + 144);
    const auto *sfh_145 = buffer.data(sfh + 145);
    const auto *sfh_146 = buffer.data(sfh + 146);
    const auto *sfh_147 = buffer.data(sfh + 147);
    const auto *sfh_149 = buffer.data(sfh + 149);
    const auto *sfh_150 = buffer.data(sfh + 150);
    const auto *sfh_152 = buffer.data(sfh + 152);
    const auto *sfh_153 = buffer.data(sfh + 153);
    const auto *sfh_156 = buffer.data(sfh + 156);
    const auto *sfh_159 = buffer.data(sfh + 159);
    const auto *sfh_161 = buffer.data(sfh + 161);
    const auto *sfh_162 = buffer.data(sfh + 162);
    const auto *sfh_163 = buffer.data(sfh + 163);
    const auto *sfh_164 = buffer.data(sfh + 164);
    const auto *sfh_165 = buffer.data(sfh + 165);
    const auto *sfh_166 = buffer.data(sfh + 166);
    const auto *sfh_167 = buffer.data(sfh + 167);
    const auto *sfh_168 = buffer.data(sfh + 168);
    const auto *sfh_170 = buffer.data(sfh + 170);
    const auto *sfh_171 = buffer.data(sfh + 171);
    const auto *sfh_173 = buffer.data(sfh + 173);
    const auto *sfh_174 = buffer.data(sfh + 174);
    const auto *sfh_177 = buffer.data(sfh + 177);
    const auto *sfh_178 = buffer.data(sfh + 178);
    const auto *sfh_180 = buffer.data(sfh + 180);
    const auto *sfh_183 = buffer.data(sfh + 183);
    const auto *sfh_184 = buffer.data(sfh + 184);
    const auto *sfh_185 = buffer.data(sfh + 185);
    const auto *sfh_186 = buffer.data(sfh + 186);
    const auto *sfh_187 = buffer.data(sfh + 187);
    const auto *sfh_188 = buffer.data(sfh + 188);

#pragma omp simd aligned(t_124, t_125, t_126, pb_x, pb_y, pc_x, pc_y, sdi0_70, sdi0_124, \
                         sdh_51, sdh_96, sdi1_70, sdi1_124, sfh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_x[k] * sdi0_124[k]
                   + f_12 * sdh_96[k]
                   - f_10 * pc_x[k] * sdi1_124[k];

        t_125[k] = f_11 * sdh_51[k]
                   + f_3 * pc_y[k] * sfh_93[k];

        t_126[k] = pb_y[k] * sdi0_70[k]
                   - f_10 * pc_y[k] * sdi1_70[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pc_x, sdh_99, sdh_100, sdh_101, \
                         sdh_102, sdh_103, sfh_99, sfh_100, sfh_101, sfh_102, \
                         sfh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_11 * sdh_99[k]
                   + f_3 * pc_x[k] * sfh_99[k];

        t_128[k] = f_11 * sdh_100[k]
                   + f_3 * pc_x[k] * sfh_100[k];

        t_129[k] = f_11 * sdh_101[k]
                   + f_3 * pc_x[k] * sfh_101[k];

        t_130[k] = f_11 * sdh_102[k]
                   + f_3 * pc_x[k] * sfh_102[k];

        t_131[k] = f_11 * sdh_103[k]
                   + f_3 * pc_x[k] * sfh_103[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, pc_x, pc_z, sdi0_133, sdi0_135, \
                         sdh_36, sdh_104, sdi1_133, sdi1_135, sfh_99, \
                         sfh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_11 * sdh_104[k]
                   + f_3 * pc_x[k] * sfh_104[k];

        t_133[k] = pb_x[k] * sdi0_133[k]
                   - f_10 * pc_x[k] * sdi1_133[k];

        t_134[k] = f_11 * sdh_36[k]
                   + f_3 * pc_z[k] * sfh_99[k];

        t_135[k] = pb_x[k] * sdi0_135[k]
                   - f_10 * pc_x[k] * sdi1_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pb_x, pc_x, pc_y, sdi0_136, sdi0_137, \
                         sdi0_139, sdh_62, sdi1_136, sdi1_137, sdi1_139, \
                         sfh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * sdi0_136[k]
                   - f_10 * pc_x[k] * sdi1_136[k];

        t_137[k] = pb_x[k] * sdi0_137[k]
                   - f_10 * pc_x[k] * sdi1_137[k];

        t_138[k] = f_11 * sdh_62[k]
                   + f_3 * pc_y[k] * sfh_104[k];

        t_139[k] = pb_x[k] * sdi0_139[k]
                   - f_10 * pc_x[k] * sdi1_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_x, pc_x, pc_y, pc_z, sdi0_140, \
                         sdi0_143, sdh_42, sdh_105, sdh_108, sdi1_140, sdi1_143, \
                         sfh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pb_x[k] * sdi0_140[k]
                   + f_14 * sdh_105[k]
                   - f_10 * pc_x[k] * sdi1_140[k];

        t_141[k] = f_3 * pc_y[k] * sfh_105[k];

        t_142[k] = f_12 * sdh_42[k]
                   + f_3 * pc_z[k] * sfh_105[k];

        t_143[k] = pb_x[k] * sdi0_143[k]
                   + f_13 * sdh_108[k]
                   - f_10 * pc_x[k] * sdi1_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pb_x, pc_x, pc_y, sdi0_145, sdi0_146, sdh_110, \
                         sdh_111, sdi1_145, sdi1_146, sfh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_3 * pc_y[k] * sfh_107[k];

        t_145[k] = pb_x[k] * sdi0_145[k]
                   + f_13 * sdh_110[k]
                   - f_10 * pc_x[k] * sdi1_145[k];

        t_146[k] = pb_x[k] * sdi0_146[k]
                   + f_0 * sdh_111[k]
                   - f_10 * pc_x[k] * sdi1_146[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_x, pc_x, pc_y, pc_z, sdi0_149, sdh_45, \
                         sdh_114, sdi1_149, sfh_108, sfh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_12 * sdh_45[k]
                   + f_3 * pc_z[k] * sfh_108[k];

        t_148[k] = f_3 * pc_y[k] * sfh_110[k];

        t_149[k] = pb_x[k] * sdi0_149[k]
                   + f_0 * sdh_114[k]
                   - f_10 * pc_x[k] * sdi1_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pb_x, pc_x, pc_z, sdi0_150, sdi0_152, sdh_48, \
                         sdh_115, sdh_117, sdi1_150, sdi1_152, \
                         sfh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_x[k] * sdi0_150[k]
                   + f_12 * sdh_115[k]
                   - f_10 * pc_x[k] * sdi1_150[k];

        t_151[k] = f_12 * sdh_48[k]
                   + f_3 * pc_z[k] * sfh_111[k];

        t_152[k] = pb_x[k] * sdi0_152[k]
                   + f_12 * sdh_117[k]
                   - f_10 * pc_x[k] * sdi1_152[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pb_x, pc_x, pc_y, sdi0_154, sdh_119, \
                         sdh_120, sdh_121, sdi1_154, sfh_114, sfh_120, \
                         sfh_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_3 * pc_y[k] * sfh_114[k];

        t_154[k] = pb_x[k] * sdi0_154[k]
                   + f_12 * sdh_119[k]
                   - f_10 * pc_x[k] * sdi1_154[k];

        t_155[k] = f_11 * sdh_120[k]
                   + f_3 * pc_x[k] * sfh_120[k];

        t_156[k] = f_11 * sdh_121[k]
                   + f_3 * pc_x[k] * sfh_121[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, sdh_122, sdh_123, sdh_124, sdh_125, \
                         sfh_122, sfh_123, sfh_124, sfh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_11 * sdh_122[k]
                   + f_3 * pc_x[k] * sfh_122[k];

        t_158[k] = f_11 * sdh_123[k]
                   + f_3 * pc_x[k] * sfh_123[k];

        t_159[k] = f_11 * sdh_124[k]
                   + f_3 * pc_x[k] * sfh_124[k];

        t_160[k] = f_11 * sdh_125[k]
                   + f_3 * pc_x[k] * sfh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_x, pc_x, pc_z, sdi0_161, sdi0_163, \
                         sdi0_164, sdh_57, sdi1_161, sdi1_163, sdi1_164, \
                         sfh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_x[k] * sdi0_161[k]
                   - f_10 * pc_x[k] * sdi1_161[k];

        t_162[k] = f_12 * sdh_57[k]
                   + f_3 * pc_z[k] * sfh_120[k];

        t_163[k] = pb_x[k] * sdi0_163[k]
                   - f_10 * pc_x[k] * sdi1_163[k];

        t_164[k] = pb_x[k] * sdi0_164[k]
                   - f_10 * pc_x[k] * sdi1_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pb_x, pc_x, pc_y, sdi0_165, sdi0_167, \
                         sdi1_165, sdi1_167, sfg0_90, sfg1_90, sfh_125, \
                         sfh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pb_x[k] * sdi0_165[k]
                   - f_10 * pc_x[k] * sdi1_165[k];

        t_166[k] = f_3 * pc_y[k] * sfh_125[k];

        t_167[k] = pb_x[k] * sdi0_167[k]
                   - f_10 * pc_x[k] * sdi1_167[k];

        t_168[k] = f_1 * sfg0_90[k]
                   - f_2 * sfg1_90[k]
                   + f_3 * pc_x[k] * sfh_126[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pc_x, pc_y, pc_z, sdh_63, sdh_65, \
                         sfg0_93, sfg1_93, sfh_126, sfh_128, sfh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_0 * sdh_63[k]
                   + f_3 * pc_y[k] * sfh_126[k];

        t_170[k] = f_3 * pc_z[k] * sfh_126[k];

        t_171[k] = f_4 * sfg0_93[k]
                   - f_5 * sfg1_93[k]
                   + f_3 * pc_x[k] * sfh_129[k];

        t_172[k] = f_0 * sdh_65[k]
                   + f_3 * pc_y[k] * sfh_128[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pc_x, pc_y, pc_z, sdh_68, sfg0_95, \
                         sfg0_96, sfg1_95, sfg1_96, sfh_129, sfh_131, \
                         sfh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_4 * sfg0_95[k]
                   - f_5 * sfg1_95[k]
                   + f_3 * pc_x[k] * sfh_131[k];

        t_174[k] = f_6 * sfg0_96[k]
                   - f_7 * sfg1_96[k]
                   + f_3 * pc_x[k] * sfh_132[k];

        t_175[k] = f_3 * pc_z[k] * sfh_129[k];

        t_176[k] = f_0 * sdh_68[k]
                   + f_3 * pc_y[k] * sfh_131[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pc_x, pc_z, sfg0_99, sfg0_100, sfg0_102, \
                         sfg1_99, sfg1_100, sfg1_102, sfh_132, sfh_135, sfh_136, \
                         sfh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_6 * sfg0_99[k]
                   - f_7 * sfg1_99[k]
                   + f_3 * pc_x[k] * sfh_135[k];

        t_178[k] = f_8 * sfg0_100[k]
                   - f_9 * sfg1_100[k]
                   + f_3 * pc_x[k] * sfh_136[k];

        t_179[k] = f_3 * pc_z[k] * sfh_132[k];

        t_180[k] = f_8 * sfg0_102[k]
                   - f_9 * sfg1_102[k]
                   + f_3 * pc_x[k] * sfh_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, pc_x, pc_y, sdh_72, sfg0_104, \
                         sfg1_104, sfh_135, sfh_140, sfh_141, sfh_142, \
                         sfh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_0 * sdh_72[k]
                   + f_3 * pc_y[k] * sfh_135[k];

        t_182[k] = f_8 * sfg0_104[k]
                   - f_9 * sfg1_104[k]
                   + f_3 * pc_x[k] * sfh_140[k];

        t_183[k] = f_3 * pc_x[k] * sfh_141[k];

        t_184[k] = f_3 * pc_x[k] * sfh_142[k];

        t_185[k] = f_3 * pc_x[k] * sfh_143[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pc_x, pc_y, pc_z, sdh_78, \
                         sfg0_100, sfg1_100, sfh_141, sfh_144, sfh_145, \
                         sfh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_3 * pc_x[k] * sfh_144[k];

        t_187[k] = f_3 * pc_x[k] * sfh_145[k];

        t_188[k] = f_3 * pc_x[k] * sfh_146[k];

        t_189[k] = f_0 * sdh_78[k]
                   + f_1 * sfg0_100[k]
                   - f_2 * sfg1_100[k]
                   + f_3 * pc_y[k] * sfh_141[k];

        t_190[k] = f_3 * pc_z[k] * sfh_141[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pc_y, sdh_80, sdh_81, sdh_82, sfg0_102, \
                         sfg0_103, sfg0_104, sfg1_102, sfg1_103, sfg1_104, sfh_143, sfh_144, \
                         sfh_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_0 * sdh_80[k]
                   + f_4 * sfg0_102[k]
                   - f_5 * sfg1_102[k]
                   + f_3 * pc_y[k] * sfh_143[k];

        t_192[k] = f_0 * sdh_81[k]
                   + f_6 * sfg0_103[k]
                   - f_7 * sfg1_103[k]
                   + f_3 * pc_y[k] * sfh_144[k];

        t_193[k] = f_0 * sdh_82[k]
                   + f_8 * sfg0_104[k]
                   - f_9 * sfg1_104[k]
                   + f_3 * pc_y[k] * sfh_145[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pb_z, pc_y, pc_z, sdi0_84, sdh_83, \
                         sdh_84, sdi1_84, sfg0_104, sfg1_104, sfh_146, \
                         sfh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_0 * sdh_83[k]
                   + f_3 * pc_y[k] * sfh_146[k];

        t_195[k] = f_1 * sfg0_104[k]
                   - f_2 * sfg1_104[k]
                   + f_3 * pc_z[k] * sfh_146[k];

        t_196[k] = pb_z[k] * sdi0_84[k]
                   - f_10 * pc_z[k] * sdi1_84[k];

        t_197[k] = f_12 * sdh_84[k]
                   + f_3 * pc_y[k] * sfh_147[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_z, pc_y, pc_z, sdi0_87, sdh_63, sdh_86, \
                         sdi1_87, sfh_147, sfh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_11 * sdh_63[k]
                   + f_3 * pc_z[k] * sfh_147[k];

        t_199[k] = pb_z[k] * sdi0_87[k]
                   - f_10 * pc_z[k] * sdi1_87[k];

        t_200[k] = f_12 * sdh_86[k]
                   + f_3 * pc_y[k] * sfh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pb_z, pc_x, pc_y, pc_z, sdi0_90, sdh_66, \
                         sdh_89, sdi1_90, sfg0_110, sfg1_110, sfh_150, \
                         sfh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_4 * sfg0_110[k]
                   - f_5 * sfg1_110[k]
                   + f_3 * pc_x[k] * sfh_152[k];

        t_202[k] = pb_z[k] * sdi0_90[k]
                   - f_10 * pc_z[k] * sdi1_90[k];

        t_203[k] = f_11 * sdh_66[k]
                   + f_3 * pc_z[k] * sfh_150[k];

        t_204[k] = f_12 * sdh_89[k]
                   + f_3 * pc_y[k] * sfh_152[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pb_z, pc_x, pc_z, sdi0_94, sdh_69, sdi1_94, \
                         sfg0_114, sfg1_114, sfh_153, sfh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_6 * sfg0_114[k]
                   - f_7 * sfg1_114[k]
                   + f_3 * pc_x[k] * sfh_156[k];

        t_206[k] = pb_z[k] * sdi0_94[k]
                   - f_10 * pc_z[k] * sdi1_94[k];

        t_207[k] = f_11 * sdh_69[k]
                   + f_3 * pc_z[k] * sfh_153[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_x, pc_y, sdh_93, sfg0_117, sfg0_119, \
                         sfg1_117, sfg1_119, sfh_156, sfh_159, sfh_161, \
                         sfh_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_8 * sfg0_117[k]
                   - f_9 * sfg1_117[k]
                   + f_3 * pc_x[k] * sfh_159[k];

        t_209[k] = f_12 * sdh_93[k]
                   + f_3 * pc_y[k] * sfh_156[k];

        t_210[k] = f_8 * sfg0_119[k]
                   - f_9 * sfg1_119[k]
                   + f_3 * pc_x[k] * sfh_161[k];

        t_211[k] = f_3 * pc_x[k] * sfh_162[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, sdi0_105, \
                         sdi1_105, sfh_163, sfh_164, sfh_165, sfh_166, \
                         sfh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_x[k] * sfh_163[k];

        t_213[k] = f_3 * pc_x[k] * sfh_164[k];

        t_214[k] = f_3 * pc_x[k] * sfh_165[k];

        t_215[k] = f_3 * pc_x[k] * sfh_166[k];

        t_216[k] = f_3 * pc_x[k] * sfh_167[k];

        t_217[k] = pb_z[k] * sdi0_105[k]
                   - f_10 * pc_z[k] * sdi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pb_z, pc_z, sdi0_107, sdi0_108, sdh_78, sdh_79, \
                         sdh_80, sdi1_107, sdi1_108, sfh_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * sdh_78[k]
                   + f_3 * pc_z[k] * sfh_162[k];

        t_219[k] = pb_z[k] * sdi0_107[k]
                   + f_12 * sdh_79[k]
                   - f_10 * pc_z[k] * sdi1_107[k];

        t_220[k] = pb_z[k] * sdi0_108[k]
                   + f_0 * sdh_80[k]
                   - f_10 * pc_z[k] * sdi1_108[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pb_z, pc_y, pc_z, sdi0_109, sdh_81, sdh_83, \
                         sdh_104, sdi1_109, sfg0_119, sfg1_119, \
                         sfh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pb_z[k] * sdi0_109[k]
                   + f_13 * sdh_81[k]
                   - f_10 * pc_z[k] * sdi1_109[k];

        t_222[k] = f_12 * sdh_104[k]
                   + f_3 * pc_y[k] * sfh_167[k];

        t_223[k] = f_11 * sdh_83[k]
                   + f_1 * sfg0_119[k]
                   - f_2 * sfg1_119[k]
                   + f_3 * pc_z[k] * sfh_167[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pb_y, pc_x, pc_y, pc_z, sdi0_140, sdh_84, \
                         sdh_105, sdi1_140, sfg0_123, sfg1_123, sfh_168, \
                         sfh_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pb_y[k] * sdi0_140[k]
                   - f_10 * pc_y[k] * sdi1_140[k];

        t_225[k] = f_11 * sdh_105[k]
                   + f_3 * pc_y[k] * sfh_168[k];

        t_226[k] = f_12 * sdh_84[k]
                   + f_3 * pc_z[k] * sfh_168[k];

        t_227[k] = f_4 * sfg0_123[k]
                   - f_5 * sfg1_123[k]
                   + f_3 * pc_x[k] * sfh_171[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pb_y, pc_x, pc_y, sdi0_145, sdh_107, sdi1_145, \
                         sfg0_126, sfg1_126, sfh_170, sfh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_11 * sdh_107[k]
                   + f_3 * pc_y[k] * sfh_170[k];

        t_229[k] = pb_y[k] * sdi0_145[k]
                   - f_10 * pc_y[k] * sdi1_145[k];

        t_230[k] = f_6 * sfg0_126[k]
                   - f_7 * sfg1_126[k]
                   + f_3 * pc_x[k] * sfh_174[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pb_y, pc_y, pc_z, sdi0_149, sdh_87, sdh_110, \
                         sdi1_149, sfh_171, sfh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_12 * sdh_87[k]
                   + f_3 * pc_z[k] * sfh_171[k];

        t_232[k] = f_11 * sdh_110[k]
                   + f_3 * pc_y[k] * sfh_173[k];

        t_233[k] = pb_y[k] * sdi0_149[k]
                   - f_10 * pc_y[k] * sdi1_149[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_z, sdh_90, sfg0_130, sfg0_132, \
                         sfg1_130, sfg1_132, sfh_174, sfh_178, \
                         sfh_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_8 * sfg0_130[k]
                   - f_9 * sfg1_130[k]
                   + f_3 * pc_x[k] * sfh_178[k];

        t_235[k] = f_12 * sdh_90[k]
                   + f_3 * pc_z[k] * sfh_174[k];

        t_236[k] = f_8 * sfg0_132[k]
                   - f_9 * sfg1_132[k]
                   + f_3 * pc_x[k] * sfh_180[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pb_y, pc_x, pc_y, sdi0_154, \
                         sdh_114, sdi1_154, sfh_177, sfh_183, sfh_184, \
                         sfh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_11 * sdh_114[k]
                   + f_3 * pc_y[k] * sfh_177[k];

        t_238[k] = pb_y[k] * sdi0_154[k]
                   - f_10 * pc_y[k] * sdi1_154[k];

        t_239[k] = f_3 * pc_x[k] * sfh_183[k];

        t_240[k] = f_3 * pc_x[k] * sfh_184[k];

        t_241[k] = f_3 * pc_x[k] * sfh_185[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pb_y, pc_x, pc_y, sdi0_161, sdh_120, \
                         sdi1_161, sfh_186, sfh_187, sfh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_3 * pc_x[k] * sfh_186[k];

        t_243[k] = f_3 * pc_x[k] * sfh_187[k];

        t_244[k] = f_3 * pc_x[k] * sfh_188[k];

        t_245[k] = pb_y[k] * sdi0_161[k]
                   + f_14 * sdh_120[k]
                   - f_10 * pc_y[k] * sdi1_161[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pb_y, pc_y, pc_z, sdi0_163, sdi0_164, sdh_99, \
                         sdh_122, sdh_123, sdi1_163, sdi1_164, \
                         sfh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * sdh_99[k]
                   + f_3 * pc_z[k] * sfh_183[k];

        t_247[k] = pb_y[k] * sdi0_163[k]
                   + f_13 * sdh_122[k]
                   - f_10 * pc_y[k] * sdi1_163[k];

        t_248[k] = pb_y[k] * sdi0_164[k]
                   + f_0 * sdh_123[k]
                   - f_10 * pc_y[k] * sdi1_164[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pb_y, pc_y, sdi0_165, sdi0_167, sdh_124, \
                         sdh_125, sdi1_165, sdi1_167, sfh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pb_y[k] * sdi0_165[k]
                   + f_12 * sdh_124[k]
                   - f_10 * pc_y[k] * sdi1_165[k];

        t_250[k] = f_11 * sdh_125[k]
                   + f_3 * pc_y[k] * sfh_188[k];

        t_251[k] = pb_y[k] * sdi0_167[k]
                   - f_10 * pc_y[k] * sdi1_167[k];
    }
}

static auto
compute_prim_sfi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t sdh, const size_t sfg0,
                                                          const size_t sfg1, const size_t sfh,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.5 / gamma;
    const auto f_5 = 1.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdh_105 = buffer.data(sdh + 105);
    const auto *sdh_108 = buffer.data(sdh + 108);
    const auto *sdh_111 = buffer.data(sdh + 111);
    const auto *sdh_120 = buffer.data(sdh + 120);
    const auto *sdh_125 = buffer.data(sdh + 125);

    const auto *sfg0_135 = buffer.data(sfg0 + 135);
    const auto *sfg0_138 = buffer.data(sfg0 + 138);
    const auto *sfg0_140 = buffer.data(sfg0 + 140);
    const auto *sfg0_141 = buffer.data(sfg0 + 141);
    const auto *sfg0_144 = buffer.data(sfg0 + 144);
    const auto *sfg0_145 = buffer.data(sfg0 + 145);
    const auto *sfg0_147 = buffer.data(sfg0 + 147);
    const auto *sfg0_148 = buffer.data(sfg0 + 148);
    const auto *sfg0_149 = buffer.data(sfg0 + 149);

    const auto *sfg1_135 = buffer.data(sfg1 + 135);
    const auto *sfg1_138 = buffer.data(sfg1 + 138);
    const auto *sfg1_140 = buffer.data(sfg1 + 140);
    const auto *sfg1_141 = buffer.data(sfg1 + 141);
    const auto *sfg1_144 = buffer.data(sfg1 + 144);
    const auto *sfg1_145 = buffer.data(sfg1 + 145);
    const auto *sfg1_147 = buffer.data(sfg1 + 147);
    const auto *sfg1_148 = buffer.data(sfg1 + 148);
    const auto *sfg1_149 = buffer.data(sfg1 + 149);

    const auto *sfh_189 = buffer.data(sfh + 189);
    const auto *sfh_191 = buffer.data(sfh + 191);
    const auto *sfh_192 = buffer.data(sfh + 192);
    const auto *sfh_194 = buffer.data(sfh + 194);
    const auto *sfh_195 = buffer.data(sfh + 195);
    const auto *sfh_198 = buffer.data(sfh + 198);
    const auto *sfh_199 = buffer.data(sfh + 199);
    const auto *sfh_201 = buffer.data(sfh + 201);
    const auto *sfh_203 = buffer.data(sfh + 203);
    const auto *sfh_204 = buffer.data(sfh + 204);
    const auto *sfh_205 = buffer.data(sfh + 205);
    const auto *sfh_206 = buffer.data(sfh + 206);
    const auto *sfh_207 = buffer.data(sfh + 207);
    const auto *sfh_208 = buffer.data(sfh + 208);
    const auto *sfh_209 = buffer.data(sfh + 209);

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, pc_x, pc_y, pc_z, sdh_105, \
                         sfg0_135, sfg0_138, sfg1_135, sfg1_138, sfh_189, sfh_191, \
                         sfh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_1 * sfg0_135[k]
                   - f_2 * sfg1_135[k]
                   + f_3 * pc_x[k] * sfh_189[k];

        t_253[k] = f_3 * pc_y[k] * sfh_189[k];

        t_254[k] = f_0 * sdh_105[k]
                   + f_3 * pc_z[k] * sfh_189[k];

        t_255[k] = f_4 * sfg0_138[k]
                   - f_5 * sfg1_138[k]
                   + f_3 * pc_x[k] * sfh_192[k];

        t_256[k] = f_3 * pc_y[k] * sfh_191[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pc_x, pc_y, pc_z, sdh_108, sfg0_140, \
                         sfg0_141, sfg1_140, sfg1_141, sfh_192, sfh_194, \
                         sfh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_4 * sfg0_140[k]
                   - f_5 * sfg1_140[k]
                   + f_3 * pc_x[k] * sfh_194[k];

        t_258[k] = f_6 * sfg0_141[k]
                   - f_7 * sfg1_141[k]
                   + f_3 * pc_x[k] * sfh_195[k];

        t_259[k] = f_0 * sdh_108[k]
                   + f_3 * pc_z[k] * sfh_192[k];

        t_260[k] = f_3 * pc_y[k] * sfh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, sdh_111, sfg0_144, sfg0_145, \
                         sfg1_144, sfg1_145, sfh_195, sfh_198, \
                         sfh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_6 * sfg0_144[k]
                   - f_7 * sfg1_144[k]
                   + f_3 * pc_x[k] * sfh_198[k];

        t_262[k] = f_8 * sfg0_145[k]
                   - f_9 * sfg1_145[k]
                   + f_3 * pc_x[k] * sfh_199[k];

        t_263[k] = f_0 * sdh_111[k]
                   + f_3 * pc_z[k] * sfh_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, pc_x, pc_y, sfg0_147, sfg0_149, \
                         sfg1_147, sfg1_149, sfh_198, sfh_201, sfh_203, sfh_204, \
                         sfh_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_8 * sfg0_147[k]
                   - f_9 * sfg1_147[k]
                   + f_3 * pc_x[k] * sfh_201[k];

        t_265[k] = f_3 * pc_y[k] * sfh_198[k];

        t_266[k] = f_8 * sfg0_149[k]
                   - f_9 * sfg1_149[k]
                   + f_3 * pc_x[k] * sfh_203[k];

        t_267[k] = f_3 * pc_x[k] * sfh_204[k];

        t_268[k] = f_3 * pc_x[k] * sfh_205[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, pc_x, pc_y, sfg0_145, sfg1_145, \
                         sfh_204, sfh_206, sfh_207, sfh_208, sfh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_3 * pc_x[k] * sfh_206[k];

        t_270[k] = f_3 * pc_x[k] * sfh_207[k];

        t_271[k] = f_3 * pc_x[k] * sfh_208[k];

        t_272[k] = f_3 * pc_x[k] * sfh_209[k];

        t_273[k] = f_1 * sfg0_145[k]
                   - f_2 * sfg1_145[k]
                   + f_3 * pc_y[k] * sfh_204[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, pc_y, pc_z, sdh_120, sfg0_147, sfg0_148, \
                         sfg1_147, sfg1_148, sfh_204, sfh_206, \
                         sfh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_0 * sdh_120[k]
                   + f_3 * pc_z[k] * sfh_204[k];

        t_275[k] = f_4 * sfg0_147[k]
                   - f_5 * sfg1_147[k]
                   + f_3 * pc_y[k] * sfh_206[k];

        t_276[k] = f_6 * sfg0_148[k]
                   - f_7 * sfg1_148[k]
                   + f_3 * pc_y[k] * sfh_207[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, pc_y, pc_z, sdh_125, sfg0_149, sfg1_149, \
                         sfh_208, sfh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_8 * sfg0_149[k]
                   - f_9 * sfg1_149[k]
                   + f_3 * pc_y[k] * sfh_208[k];

        t_278[k] = f_3 * pc_y[k] * sfh_209[k];

        t_279[k] = f_0 * sdh_125[k]
                   + f_1 * sfg0_149[k]
                   - f_2 * sfg1_149[k]
                   + f_3 * pc_z[k] * sfh_209[k];
    }
}

auto
compute_prim_sfi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sdi0, const size_t sdh,
                                                   const size_t sdi1, const size_t sfg0,
                                                   const size_t sfg1, const size_t sfh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sfi_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sdi0, sdh,
                                                              sdi1, sfg0, sfg1, sfh, ncols,
                                                              gamma, p, q);

    compute_prim_sfi_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sdi0, sdh,
                                                              sdi1, sfg0, sfg1, sfh, ncols,
                                                              gamma, p, q);

    compute_prim_sfi_three_center_electron_repulsion_0_piece2(buffer, target, pc, sdh, sfg0,
                                                              sfg1, sfh, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
